"""Test real inference adapter using fake HTTP, never a live model or source service."""
from http.client import BadStatusLine
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock

from scripts import compute_fanout
from scripts import civic_evidence, grounded_claims
from scripts import openai_compatible_worker as worker


class Response:
    def __init__(self, body):
        self.body = json.dumps({"choices": [{"message": {"content": json.dumps(body)}}]}).encode()
    def __enter__(self): return self
    def __exit__(self, *args): return False
    def read(self, n): return self.body[:n]


class ModelAdapterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / 'catalog/sources').mkdir(parents=True)
        (self.root / 'catalog/sources/civic.json').write_text(json.dumps({
            'id': 'civic', 'name': 'CIViC', 'summary': 'Public evidence directory',
            'homepage': 'https://civicdb.org/',
            'documentation_url': 'https://docs.civicdb.org/en/latest/api.html'}))
        self.job = {'protocol': compute_fanout.PROTOCOL, 'task_id': 'test-task',
                    'question': 'Compare public source terminology without drawing clinical conclusions.',
                    'source_ids': ['civic'], 'dimensions': {'variant': 'EGFR L858R'}}
        self.job['id'] = compute_fanout.stable_id(self.job)
        self.env = {'OPEN_CANCER_LLM_URL': 'http://127.0.0.1:8000/v1/chat/completions',
                    'OPEN_CANCER_LLM_MODEL': 'local-model', 'OPEN_CANCER_LLM_KEY': 'secret'}

    def test_investigator_validates_context_and_never_rewrites_sources(self):
        model = {'summary': 'An unverified computational question about cancer variant terminology.',
                 'evidence': [{'source_id': 'civic', 'url': 'https://civicdb.org/',
                               'relation': 'context'}],
                 'limitations': ['No actual source evidence items were fetched or examined.'],
                 'synthetic': False}
        seen = []
        def opener(request, timeout):
            seen.append(request)
            return Response(model)
        result = worker.work({'protocol': compute_fanout.PROTOCOL, 'stage': 'investigate',
                              'job': self.job}, self.env, self.root, opener)
        self.assertEqual(model, result)
        payload = json.loads(seen[0].data)
        self.assertEqual('local-model', payload['model'])
        self.assertIn('DIRECTORY', payload['messages'][0]['content'])
        self.assertEqual('Bearer secret', seen[0].headers['Authorization'])
        self.assertNotIn('secret', json.dumps(result))

    def test_model_can_not_claim_verified_record_from_directory_only(self):
        record = {'summary': 'Model invented an evidence-supported biological finding.',
                  'evidence': [{'source_id': 'civic', 'url': 'https://civicdb.org/links/evidence_items/999',
                                'relation': 'supports'}],
                  'limitations': ['The retrieved record was not independently verified.'],
                  'synthetic': False}
        with self.assertRaisesRegex(ValueError, 'directory only'):
            worker.work({'protocol': compute_fanout.PROTOCOL, 'stage': 'investigate',
                         'job': self.job}, self.env, self.root, lambda *a, **k: Response(record))

    def test_critic_is_always_downgraded_when_records_absent(self):
        verdict = {'verdict': 'no-internal-issues-found', 'issues': []}
        result = worker.work({'protocol': compute_fanout.PROTOCOL, 'stage': 'critique',
                              'job': self.job, 'candidate': {'summary': 'unverified'}},
                             self.env, self.root, lambda *a, **k: Response(verdict))
        self.assertEqual('insufficient-evidence', result['verdict'])
        self.assertTrue(result['issues'])

    def test_endpoint_rejects_plain_http_off_machine_and_credentials_in_url(self):
        for url in ['http://example.org/v1/chat/completions',
                    'https://user:pw@example.org/v1/chat/completions',
                    'file:///tmp/test', 'http://localhost/v1/not-chat']:
            with self.subTest(url=url), self.assertRaises(ValueError):
                worker.validate_endpoint(url)

    def test_bad_model_response_fails_not_publishes(self):
        with self.assertRaisesRegex(ValueError, 'JSON message content'):
            worker.extract_response({'choices': [{'message': {'content': 'not-json'}}]})

    def test_inference_status_line_and_redirect_are_not_forwarded(self):
        failing = Mock(side_effect=BadStatusLine('untrusted status line'))
        with self.assertRaisesRegex(RuntimeError, '^Inference endpoint request failed$'):
            worker.work({'protocol': compute_fanout.PROTOCOL, 'stage': 'investigate',
                         'job': self.job}, self.env, self.root, failing)
        redirect = worker.NoRedirect()
        self.assertIsNone(redirect.redirect_request(None, None, 302, '', {}, 'https://elsewhere/'))


class GroundedModelAdapterTests(unittest.TestCase):
    def setUp(self):
        self.catalog = json.loads((worker.ROOT / 'catalog/sources/civic.json').read_text())
        self.bundle = civic_evidence.retrieve_fixture(self.catalog)
        self.job = {'protocol': compute_fanout.PROTOCOL, 'task_id': 'egfr-luad-evidence-coverage',
                    'question': 'Which public CIViC annotation directions occur in this sample?',
                    'source_ids': ['civic'], 'dimensions': {'alteration': 'EGFR L858R'}}
        self.job['id'] = compute_fanout.stable_id(self.job)
        self.request = {'protocol': grounded_claims.GROUNDED_PROTOCOL, 'stage': 'investigate',
                        'job': self.job, 'evidence_bundle': self.bundle}
        self.env = {'OPEN_CANCER_LLM_URL': 'http://127.0.0.1:8000/v1/chat/completions',
                    'OPEN_CANCER_INVESTIGATOR_MODEL': 'investigator',
                    'OPEN_CANCER_CRITIC_MODEL': 'critic'}

    def test_models_receive_same_scoped_records_and_critic_can_cite_counterexample(self):
        records = civic_evidence.records_for_job(self.bundle, self.job)
        support = next(r for r in records if r['evidence_direction'] == 'SUPPORTS')
        conflict = next(r for r in records if r['evidence_direction'] == 'DOES_NOT_SUPPORT')
        def citation(record, relation):
            return {'record_id': record['record_id'], 'sha256': record['sha256'],
                    'relation': relation}
        proposal = {
            'synthetic': True, 'outcome': 'record-observation',
            'observation': {'field': 'evidence_direction', 'value': 'SUPPORTS',
                            'evidence_type': support['evidence_type']},
            'evidence': [citation(support, 'supports')], 'limitations': ['record-annotations-only']
        }
        challenge = {'verdict': 'inconsistent', 'issues': ['annotation-counterexample'],
                     'evidence': [citation(conflict, 'contradicts')]}
        opener = Mock(side_effect=[Response(proposal), Response(challenge)])
        candidate = worker.work(self.request, self.env, opener=opener)
        critique = worker.work(dict(self.request, stage='critique', candidate=candidate),
                               self.env, opener=opener)
        self.assertEqual(candidate, proposal)
        self.assertEqual(critique, challenge)
        calls = [json.loads(call.args[0].data) for call in opener.call_args_list]
        self.assertEqual(['investigator', 'critic'], [call['model'] for call in calls])
        contexts = [json.loads(call['messages'][1]['content']) for call in calls]
        for context in contexts:
            self.assertEqual(records, context['records'])
            self.assertTrue(context['synthetic'])
            self.assertNotIn('source_directory', context)
        self.assertNotIn('candidate', contexts[0])
        self.assertEqual(candidate, contexts[1]['candidate'])

        proposal['evidence'][0]['record_id'] = 'EID999999999'
        with self.assertRaises(ValueError):
            worker.work(self.request, self.env, opener=Mock(return_value=Response(proposal)))

    def test_empty_scope_short_circuits_both_stages_without_credentials(self):
        self.job['dimensions']['alteration'] = 'EGFR exon 20 insertion'
        self.job['id'] = compute_fanout.stable_id(self.job)
        opener = Mock(side_effect=AssertionError('inference should be skipped'))
        candidate = worker.work(self.request, {}, opener=opener)
        self.assertEqual('insufficient-evidence', candidate['outcome'])
        critique = worker.work(dict(self.request, stage='critique', candidate=candidate),
                               {}, opener=opener)
        self.assertEqual('insufficient-evidence', critique['verdict'])
        opener.assert_not_called()

    def test_invalid_bundle_and_invalid_candidate_never_reach_the_model(self):
        opener = Mock(side_effect=AssertionError('inference should be skipped'))
        with self.assertRaises(ValueError):
            worker.work(dict(self.request, stage='critique', candidate={'summary': 'untrusted'}),
                        self.env, opener=opener)
        self.bundle['records'][0]['description'] = 'unexpected upstream free text'
        with self.assertRaises(ValueError):
            worker.work(self.request, self.env, opener=opener)
        opener.assert_not_called()

    def test_version_mismatch_cannot_silently_fall_back_to_directory(self):
        with self.assertRaisesRegex(ValueError, 'Unsupported worker request'):
            worker.work(dict(self.request, protocol=compute_fanout.PROTOCOL), {})
        del self.request['evidence_bundle']
        with self.assertRaisesRegex(ValueError, 'Unsupported grounded worker request'):
            worker.work(self.request, {})


if __name__ == '__main__':
    unittest.main()
