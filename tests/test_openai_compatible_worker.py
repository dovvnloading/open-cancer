"""Test real inference adapter using fake HTTP, never a live model or source service."""
import json
from pathlib import Path
import tempfile
import unittest

from scripts import compute_fanout
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


if __name__ == '__main__':
    unittest.main()
