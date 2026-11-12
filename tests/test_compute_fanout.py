"""Offline tests for compute-first execution, failure modes, and provenance."""
import json
from pathlib import Path
import shlex
import sys
import tempfile
import unittest

from scripts import compute_fanout as compute

ROOT = Path(__file__).resolve().parents[1]
DEMO = f"{shlex.quote(sys.executable)} {shlex.quote(str(ROOT / 'scripts/demo_compute_worker.py'))}"


class ComputeFanoutTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        (self.root / 'tasks').mkdir()
        (self.root / 'catalog/sources').mkdir(parents=True)
        (self.root / 'tasks/cancer-task.json').write_text(json.dumps({
            'id': 'cancer-task', 'input_source_ids': ['civic', 'gdc']}))
        for source in ('civic', 'gdc'):
            (self.root / 'catalog/sources' / (source + '.json')).write_text('{}')
        self.spec = {
            'schema_version': 1, 'task_id': 'cancer-task',
            'question': 'Compare public source metadata without inferring any clinical results.',
            'source_ids': ['gdc', 'civic'],
            'axes': {'gene': ['EGFR', 'KRAS'], 'source-check': ['denominator', 'terminology']},
        }
        self.jobs = compute.plan(self.spec, 4, self.root)

    def test_deterministic_plan_and_bound(self):
        self.assertEqual(4, len(self.jobs))
        reverse = dict(self.spec)
        reverse['axes'] = dict(reversed(list(self.spec['axes'].items())))
        reverse['source_ids'] = list(reversed(self.spec['source_ids']))
        self.assertEqual(self.jobs, compute.plan(reverse, 4, self.root))
        with self.assertRaisesRegex(ValueError, 'exceeds --max-jobs'):
            compute.plan(self.spec, 3, self.root)
        invalid = dict(self.spec, source_ids=['private-source'])
        with self.assertRaisesRegex(ValueError, 'checked-in task'):
            compute.plan(invalid, 100, self.root)

    def test_worker_output_and_resume_with_hash_binding(self):
        out = self.root / 'out'
        def run(resume=False):
            return compute.run_batch(self.jobs, out, DEMO, DEMO, 3, 5, 0, 4, 0, 1, resume)
        summary = run()
        self.assertEqual(4, summary['completed'])
        self.assertEqual(0, summary['failed'])
        self.assertEqual({'synthetic': 4}, summary['triage'])
        for job in self.jobs:
            record = compute.read_envelope(out / 'records' / (job['id'] + '.json'))
            self.assertEqual('unreviewed', record['review_status'])
            self.assertEqual('not-published', record['publication_status'])
            self.assertEqual('not-assessed', record['clinical_validation_status'])
            self.assertEqual('insufficient-evidence', record['critique']['verdict'])
            self.assertTrue(record['candidate']['synthetic'])
        first = {p.name: p.read_bytes() for p in (out / 'records').iterdir()}
        self.assertEqual(summary, run(resume=True))
        self.assertEqual(first, {p.name: p.read_bytes() for p in (out / 'records').iterdir()})
        with self.assertRaisesRegex(ValueError, 'pass --resume'):
            run(resume=False)
        with self.assertRaisesRegex(ValueError, 'configuration differs'):
            compute.run_batch(self.jobs, out, DEMO, DEMO, 3, 6, 0, 4, 0, 1, True)
        record_path = next((out / 'records').iterdir())
        content = json.loads(record_path.read_text())
        content['body']['triage'] = 'tampered'
        record_path.write_text(json.dumps(content))
        with self.assertRaisesRegex(ValueError, 'modified checkpoint'):
            run(resume=True)

    def test_shards_do_not_overlap_and_cover_all(self):
        shard_ids = []
        for i in range(3):
            out = self.root / f'shard{i}'
            summary = compute.run_batch(self.jobs, out, DEMO, DEMO, 2, 5, 0, 4, i, 3, False)
            self.assertEqual(0, summary['failed'])
            shard_ids += [r['id'] for r in summary['records']]
        self.assertEqual(sorted(j['id'] for j in self.jobs), sorted(shard_ids))
        self.assertEqual(len(shard_ids), len(set(shard_ids)))

    def test_missing_or_bad_worker_fails_closed(self):
        invalid_worker = f"{shlex.quote(sys.executable)} -c {shlex.quote('print(\"not json\")')}"
        summary = compute.run_batch(self.jobs[:1], self.root / 'bad', invalid_worker,
                                    DEMO, 1, 5, 1, 1, 0, 1, False)
        self.assertEqual(1, summary['failed'])
        record = compute.read_envelope(next((self.root / 'bad/records').iterdir()))
        self.assertIn('invalid-worker-response', record['error'])
        self.assertEqual('not-published', record['publication_status'])
        self.assertNotIn('critique', record)

    def test_validation_rejects_untrusted_protocol_and_duplicates(self):
        job = dict(self.jobs[0]); job['id'] = 'cancer-task-fake'
        with self.assertRaisesRegex(ValueError, 'digest'):
            compute.validate_job(job)
        jobs_path = self.root / 'jobs.jsonl'
        jobs_path.write_text(json.dumps(self.jobs[0]) + '\n' + json.dumps(self.jobs[0]) + '\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate job'):
            compute.load_jobs(jobs_path)
        evidence = {'summary': 'A sufficiently long but wholly invented result.',
                    'evidence': [{'source_id': 'private-source', 'url': 'https://bad.invalid',
                                  'relation': 'supports'}],
                    'limitations': ['This is a placeholder limitation only.'], 'synthetic': False}
        with self.assertRaisesRegex(ValueError, 'Evidence needs'):
            compute.check_response('investigate', evidence, self.jobs[0])
        with self.assertRaisesRegex(ValueError, 'explain issues'):
            compute.check_response('critique', {'verdict': 'inconsistent', 'issues': []}, self.jobs[0])


if __name__ == '__main__':
    unittest.main()
