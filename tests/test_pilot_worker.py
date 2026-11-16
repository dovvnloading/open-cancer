"""Exercise the production message and transport protocol with software test doubles."""
from http.server import BaseHTTPRequestHandler, HTTPServer
import io
import json
from threading import Thread
import unittest
from unittest.mock import Mock, patch
from urllib.error import HTTPError

from scripts import civic_evidence, compute_fanout, openai_compatible_worker as adapter
from scripts import pilot_evidence, pilot_worker


def reply(value, **meta):
    message = {"model": "served-build-72", "choices": [{"message": {"content": json.dumps(value)}}],
               **meta}
    return io.BytesIO(json.dumps(message).encode())


class PilotWorkerTests(unittest.TestCase):
    def setUp(self):
        self.env = {
            "OPEN_CANCER_LLM_URL": "https://investigator.invalid/v1/chat/completions",
            "OPEN_CANCER_INVESTIGATOR_MODEL": "research-v1", "OPEN_CANCER_LLM_KEY": "first-secret",
            "OPEN_CANCER_CRITIC_URL": "https://critic.invalid/v1/chat/completions",
            "OPEN_CANCER_CRITIC_MODEL": "challenge-v3", "OPEN_CANCER_CRITIC_KEY": "second-secret"}
        jobs = compute_fanout.plan(compute_fanout.read_json(
            compute_fanout.ROOT / "workloads/egfr-luad-compute.json"), 9)
        self.job = next(j for j in jobs if j["dimensions"]["alteration"] == "EGFR L858R")
        with patch.object(civic_evidence, "_fetch_page", side_effect=AssertionError("Unexpected GraphQL")), \
                patch.object(pilot_evidence, "fetch", side_effect=AssertionError("Unexpected TSV")):
            self.fixture = pilot_evidence.evidence_for_job(
                pilot_evidence.prepare_snapshot("fixture"), self.job)

        # Stub at the prepared packet boundary. No source or inference service is contacted.
        self.packet = dict(self.fixture, synthetic=False)
        self.model_answer = pilot_worker._fixture_response("investigate", self.packet)

    def call(self, stage, opener, candidate=None):
        return pilot_worker.invoke(stage, self.job, self.packet, candidate, "real", self.env, opener)

    def test_separate_models_keys_and_identical_packet_and_reproducible_request(self):
        out = reply(self.model_answer, system_fingerprint="weights-30",
                    usage={"prompt_tokens": 110, "completion_tokens": 70, "total_tokens": 180})
        investigator = Mock(return_value=out)
        first = self.call("investigate", investigator)
        self.assertIsNone(first["error"])
        self.assertEqual({"prompt_tokens": 110, "completion_tokens": 70, "total_tokens": 180}, first["usage"])
        self.assertTrue(first["endpoint_attempted"])
        self.assertEqual("served-build-72", first["returned_model"])
        self.assertEqual("weights-30", first["system_fingerprint"])
        self.assertFalse(first["response"]["synthetic"])
        again = self.call("investigate", Mock(return_value=reply(self.model_answer)))
        self.assertEqual(first["request_sha256"], again["request_sha256"])

        critic_answer = pilot_worker._fixture_response("critique", self.packet)
        critic = Mock(return_value=reply(critic_answer))
        second = self.call("critique", critic, first["response"])
        self.assertIsNone(second["error"])
        self.assertNotEqual(first["request_sha256"], second["request_sha256"])
        for stage, op, role, model, bearer in [
            ("investigate", investigator, "LLM", "research-v1", "first-secret"),
            ("critique", critic, "CRITIC", "challenge-v3", "second-secret")
        ]:
            with self.subTest(stage=stage):
                req = op.call_args.args[0]
                self.assertEqual(self.env[f"OPEN_CANCER_{role}_URL"], req.full_url)
                self.assertEqual("Bearer " + bearer, req.get_header("Authorization"))
                body = json.loads(req.data)
                self.assertEqual(model, body["model"])
                self.assertEqual(0, body["temperature"])
                self.assertEqual(2400, body["max_tokens"])
                user = json.loads(body["messages"][1]["content"])
                self.assertEqual(self.packet, user["packet"])
                self.assertEqual(self.job, user["job"])
                if stage == "critique":
                    self.assertEqual(first["response"], user["candidate"])
                self.assertNotIn(bearer.encode(), req.data)
        exposed = json.dumps([first, second, pilot_worker.configuration("real", self.env)])
        for private in ("first-secret", "second-secret", "https://critic.invalid"):
            self.assertNotIn(private, exposed)

    def test_loopback_http_completion_and_refused_redirect(self):
        """Actual HTTP POST to a local stub; no inference hardware or external calls."""
        paths = []
        response = reply(self.model_answer, usage={"prompt_tokens": 8, "completion_tokens": 4}).getvalue()

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_POST(self):
                paths.append(self.path)
                self.rfile.read(int(self.headers["Content-Length"]))
                if self.path.startswith("/redirect/"):
                    self.send_response(302)
                    self.send_header("Location", "/unexpected")
                    self.end_headers()
                else:
                    self.send_response(200)
                    self.end_headers()
                    self.wfile.write(response)

            def do_GET(self):
                paths.append(self.path)
                self.send_error(400)

        with HTTPServer(("127.0.0.1", 0), Handler) as server:
            thread = Thread(target=server.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True)
            thread.start()
            try:
                local = f"http://127.0.0.1:{server.server_port}"
                self.env["OPEN_CANCER_LLM_URL"] = local + "/v1/chat/completions"
                first = self.call("investigate", None)
                self.assertEqual(4, first["usage"]["completion_tokens"])
                self.assertIsNone(first["error"])
                self.env["OPEN_CANCER_LLM_URL"] = local + "/redirect/chat/completions"
                second = self.call("investigate", None)
                self.assertEqual("Inference endpoint HTTP 302", second["error"])
                self.assertEqual(["/v1/chat/completions", "/redirect/chat/completions"], paths)
            finally:
                server.shutdown()
                thread.join()

    def test_separate_endpoint_without_own_key_does_not_receive_shared_key(self):
        del self.env["OPEN_CANCER_CRITIC_KEY"]
        candidate = dict(self.model_answer, synthetic=False)
        opener = Mock(return_value=reply(pilot_worker._fixture_response("critique", self.packet)))
        self.assertIsNone(self.call("critique", opener, candidate)["error"])
        self.assertIsNone(opener.call_args.args[0].get_header("Authorization"))
        self.assertIsNone(adapter.NoRedirect().redirect_request(
            opener.call_args.args[0], None, 302, None, {}, "https://elsewhere.invalid"))

    def test_failures_preserve_metering_without_exposing_provider_data(self):
        for choice in ({"message": {"content": "Private debug data"}},
                       {"message": {"content": json.dumps(self.model_answer)}, "finish_reason": "length"},
                       {"message": {"content": '{"issues": "Private debug data"}'}}):
            with self.subTest(choice=choice.get("finish_reason")):
                opener = Mock(return_value=reply(None, choices=[choice],
                                                 usage={"prompt_tokens": 81, "completion_tokens": False}))
                failed = self.call("investigate", opener)
                self.assertIsNone(failed["response"])
                self.assertIsNotNone(failed["error"])
                self.assertEqual(81, failed["usage"]["prompt_tokens"])
                self.assertIsNone(failed["usage"]["completion_tokens"])
                self.assertEqual("served-build-72", failed["returned_model"])
                self.assertEqual(64, len(failed["provider_response_sha256"]))
                self.assertNotIn("Private debug data", json.dumps(failed))
        http_error = HTTPError("https://private.invalid", 429, "secret", {}, io.BytesIO(b"secret"))
        failed = self.call("investigate", Mock(side_effect=http_error))
        self.assertTrue(failed["endpoint_attempted"])
        self.assertEqual("Inference endpoint HTTP 429", failed["error"])
        self.assertIsNone(failed["usage"])

    def test_fixture_never_calls_endpoint_and_rejects_mode_mixing(self):
        opener = Mock(side_effect=AssertionError("Fixture attempted an HTTP request"))
        attempt = pilot_worker.invoke("investigate", self.job, self.fixture, None, "fixture",
                                      {}, opener)
        self.assertTrue(attempt["response"]["synthetic"])
        self.assertFalse(attempt["endpoint_attempted"])
        self.assertIsNone(attempt["usage"])
        self.assertEqual(1, len(attempt["response"]["contradictory_evidence"]))
        for stage, packet, candidate in [("investigate", self.fixture, None),
                                          ("critique", self.packet, attempt["response"])]:
            failed = pilot_worker.invoke(stage, self.job, packet, candidate, "real", self.env, opener)
            self.assertEqual("invalid-input", failed["error"])
            self.assertFalse(failed["endpoint_attempted"])
        opener.assert_not_called()

    def test_output_schema_enforces_citations_and_critic_job_binding(self):
        self.model_answer["observations"][0]["citations"] = []
        self.assertEqual("invalid-model-schema",
                         self.call("investigate", Mock(return_value=reply(self.model_answer)))["error"])
        self.model_answer.update(observations=[], contradictory_evidence=[])
        empty = self.call("investigate", Mock(return_value=reply(self.model_answer,
                                                                  usage={"completion_tokens": 8})))
        self.assertEqual("invalid-model-schema", empty["error"])
        self.assertIsNone(empty["response"])
        self.assertEqual(8, empty["usage"]["completion_tokens"])
        candidate = dict(pilot_worker._fixture_response("investigate", self.packet), synthetic=False)
        critique = pilot_worker._fixture_response("critique", self.packet)
        critique["issues"][0]["observation_id"] = "from-another-job"
        self.assertEqual("invalid-model-schema",
                         self.call("critique", Mock(return_value=reply(critique)), candidate)["error"])

    def test_empty_selection_is_a_gap(self):
        self.packet["records"] = []
        answer = self.call("investigate", Mock(return_value=reply(
            pilot_worker._fixture_response("investigate", self.packet))))["response"]
        self.assertEqual([], answer["observations"])
        incorrect = {"verdict": "no-internal-issues-found", "issues": []}
        second = self.call("critique", Mock(return_value=reply(incorrect)), answer)
        self.assertEqual("invalid-model-schema", second["error"])

    def test_passing_critic_verdict_cannot_contradict_its_listed_issues(self):
        candidate = dict(self.model_answer, synthetic=False)
        contradictory = pilot_worker._fixture_response("critique", self.packet)
        contradictory["verdict"] = "no-internal-issues-found"
        answer = self.call("critique", Mock(return_value=reply(contradictory)), candidate)
        self.assertIsNone(answer["response"])
        self.assertEqual("invalid-model-schema", answer["error"])

    def test_truncated_retrieval_prevents_passing_verdict_with_retained_records(self):
        self.packet["sources"]["civic"]["retrieval_status"] = "truncated"
        candidate = dict(self.model_answer, synthetic=False)
        passing = {"verdict": "no-internal-issues-found", "issues": []}
        rejected = self.call("critique", Mock(return_value=reply(passing)), candidate)
        self.assertEqual("invalid-model-schema", rejected["error"])
        incomplete = pilot_worker._fixture_response("critique", self.packet)
        accepted = self.call("critique", Mock(return_value=reply(incomplete)), candidate)
        self.assertEqual("insufficient-evidence", accepted["response"]["verdict"])

    def test_config_binds_role_model_and_endpoint_but_not_key(self):
        config = pilot_worker.configuration("real", self.env)
        self.env["OPEN_CANCER_CRITIC_KEY"] = "rotate"
        self.assertEqual(config, pilot_worker.configuration("real", self.env))
        self.env["OPEN_CANCER_CRITIC_URL"] = "https://another.invalid/v1/chat/completions"
        updated = pilot_worker.configuration("real", self.env)
        self.assertNotEqual(config["critique"]["endpoint_sha256"], updated["critique"]["endpoint_sha256"])
        self.assertNotEqual(config["critique"]["prompt_sha256"], config["investigate"]["prompt_sha256"])
        for setting, value in (("OPEN_CANCER_LLM_TIMEOUT", "0"), ("OPEN_CANCER_MAX_TOKENS", "infinite"),
                               ("OPEN_CANCER_CRITIC_URL", "https://other.invalid:88889/v1/chat/completions")):
            with self.subTest(setting=setting), self.assertRaises(ValueError):
                pilot_worker.configuration("real", dict(self.env, **{setting: value}))


if __name__ == "__main__":
    unittest.main()
