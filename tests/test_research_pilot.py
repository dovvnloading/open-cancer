"""Integration tests for the nine-job driver. All sources and inference are synthetic or stubbed."""
from collections import Counter
from http.server import BaseHTTPRequestHandler, HTTPServer
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Thread
import unittest
from unittest.mock import patch

from scripts import civic_evidence, compute_fanout as compute, pilot_evidence as evidence, pilot_worker
from scripts import research_pilot as pilot


class ResearchPilotTests(unittest.TestCase):
    def setUp(self):
        base = pilot.ROOT / ".compute"
        base.mkdir(exist_ok=True)
        tmp = TemporaryDirectory(prefix="driver-test-", dir=base)
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.out = self.tmp / "test-run"
        # Every test opts in to exact stub responses; neither source may use an external socket.
        self.primary = self.enterContext(patch.object(civic_evidence, "_fetch_page",
            side_effect=AssertionError("Unexpected GraphQL retrieval")))
        self.tsv = self.enterContext(patch.object(evidence, "fetch",
            side_effect=AssertionError("Unexpected TSV retrieval")))

    def test_nine_jobs_checkpointed_two_independent_executions_and_no_rebilling(self):
        first = pilot.run("fixture", self.out, env={})
        self.assertEqual("finished", first["status"])
        self.assertTrue(first["synthetic"])
        self.assertEqual({"total": 9, "completed": 9, "failed": 0,
                          "with_candidate": 9, "with_critique": 9}, first["metrics"]["jobs"])
        resources = first["metrics"]["execution"]
        self.assertEqual(18, resources["attempts"])
        self.assertEqual(0, resources["endpoint_attempts"])
        self.assertIsNone(resources["usage"]["total_tokens"])
        self.assertIsNone(resources["estimated_cost_usd"])
        self.assertIn("synthetic", (self.out / "report.md").read_text().lower())
        config = first["configuration"]
        self.assertEqual(compute.digest(compute.read_json(self.out / "plan.json")), config["plan_sha256"])
        self.assertIn("scripts/pilot_worker.py", config["file_sha256"])
        self.assertEqual(2, len(config["models"]))
        self.assertEqual(first, compute.read_json(self.out / "benchmark.json"))
        files = sorted((self.out / "records").glob("*.json"))
        self.assertEqual(9, len(files))
        checkpoint_bytes = [path.read_bytes() for path in files]
        for p in files:
            record = compute.read_envelope(p)
            self.assertEqual(["investigate", "critique"], [a["stage"] for a in record["attempts"]])
            self.assertEqual("unreviewed", record["review_status"])
            self.assertEqual("not-published", record["publication_status"])
            self.assertTrue(record["nonclinical"])
            self.assertTrue(record["candidate"]["synthetic"] and record["critique"]["synthetic"])

        with self.assertRaisesRegex(ValueError, "Run directory exists"):
            pilot.run("fixture", self.out, env={})
        with patch.object(pilot_worker, "invoke", side_effect=AssertionError("Repeated finished stage")):
            resumed = pilot.run("fixture", self.out, resume=True, env={})
        self.assertEqual(first["metrics"], resumed["metrics"])
        self.assertEqual(checkpoint_bytes, [p.read_bytes() for p in files])

        other = self.tmp / "repeat"
        repeat = pilot.run("fixture", other, supplied=self.out / "evidence.json", env={})
        self.assertNotEqual(first["run_id"], repeat["run_id"])
        self.assertEqual(first["config_sha256"], repeat["config_sha256"])
        compare = pilot.compare(self.out, other)["consistency"]
        self.assertEqual(9, compare["paired_completed"])
        self.assertEqual(9, compare["both_outputs_exact"])
        self.assertEqual(9, compare["verdict_matches"])
        self.assertEqual(1.0, compare["mean_cited_record_jaccard"])
        with self.assertRaisesRegex(ValueError, "independently launched"):
            pilot.compare(self.out, self.out)
        # A lost comparison input has no output to compare. It stays in the denominator.
        next((other / "records").glob("*.json")).unlink()
        partial = pilot.compare(self.out, other)["consistency"]
        self.assertEqual(8, partial["paired_completed"])
        self.assertEqual(1, partial["missing_or_failed_pairs"])

    def test_checkpoint_corruption_rejected_before_dispatch(self):
        pilot.run("fixture", self.out, env={})
        checkpoint = next((self.out / "records").glob("*.json"))
        envelope = compute.read_json(checkpoint)
        envelope["body"]["candidate"]["assumptions"] = ["tampered"]
        compute.atomic_json(checkpoint, envelope)
        with patch.object(pilot_worker, "invoke", side_effect=AssertionError("unexpected call")):
            with self.assertRaisesRegex(ValueError, "modified checkpoint"):
                pilot.run("fixture", self.out, resume=True, env={})
            with self.assertRaisesRegex(ValueError, "modified checkpoint"):
                pilot.read_run(self.out)

    def test_missing_checkpoint_cannot_rebill_a_job(self):
        pilot.run("fixture", self.out, env={})
        next((self.out / "records").glob("*.json")).unlink()
        with patch.object(pilot_worker, "invoke", side_effect=AssertionError("Repeated missing checkpoint")):
            with self.assertRaises((ValueError, FileNotFoundError)):
                pilot.run("fixture", self.out, resume=True, env={})

    def test_mode_configuration_and_snapshot_mismatch_rejected(self):
        first = pilot.run("fixture", self.out, env={})
        with self.assertRaisesRegex(ValueError, "mode"):
            pilot.run("real", self.out, resume=True, env={})
        snap = compute.read_json(self.out / "evidence.json")
        snap["records"][0]["molecular_profile"] = "unrelated"
        altered = self.tmp / "altered.json"
        compute.atomic_json(altered, snap)
        with self.assertRaisesRegex(ValueError, "evidence differs"):
            pilot.run("fixture", self.out, supplied=altered, resume=True, env={})
        changed = dict(first["configuration"]["models"])
        changed["investigate"] = dict(changed["investigate"], model="different")
        with patch.object(pilot_worker, "configuration", return_value=changed):
            with self.assertRaisesRegex(ValueError, "configuration differs"):
                pilot.run("fixture", self.out, resume=True, env={})
            other = self.tmp / "different-config"
            pilot.run("fixture", other, env={})
        with self.assertRaisesRegex(ValueError, "configurations or evidence differ"):
            pilot.compare(self.out, other)

    def test_critic_failure_retains_candidate_and_resume_retries_only_that_critic(self):
        invoke = pilot_worker.invoke
        stages = []

        def fail_first_critic(stage, *args):
            stages.append(stage)
            attempt = invoke(stage, *args)
            if stage == "critique" and stages.count("critique") == 1:
                attempt.update(error="invalid-model-schema", response=None)
            return attempt

        with patch.object(pilot_worker, "invoke", side_effect=fail_first_critic):
            report = pilot.run("fixture", self.out, env={})
            self.assertEqual("finished-with-failures", report["status"])
            self.assertEqual(1, report["metrics"]["jobs"]["failed"])
            self.assertEqual(9, report["metrics"]["jobs"]["with_candidate"])
            self.assertEqual(8, report["metrics"]["jobs"]["with_critique"])
            failed = next(b for b in pilot.read_run(self.out)[1] if b["status"] == "failed")
            self.assertIsInstance(failed["candidate"]["observations"], list)
            self.assertIsNone(failed["critique"])
            self.assertEqual("invalid-model-schema", failed["attempts"][-1]["error"])
            resumed = pilot.run("fixture", self.out, resume=True, env={})
            self.assertEqual("finished", resumed["status"])
            self.assertEqual(Counter(investigate=9, critique=10), Counter(stages))
            self.assertEqual(1, resumed["metrics"]["execution"]["failed_attempts"])
            repaired = compute.read_envelope(pilot.record_path(self.out, failed["job"]))
            self.assertEqual(failed["candidate"], repaired["candidate"])
            self.assertEqual([1, 1, 2], [a["number"] for a in repaired["attempts"]])

    def test_interrupted_request_reservations_survive_and_bound_restarts(self):
        with patch.object(pilot_worker, "invoke", side_effect=KeyboardInterrupt) as stopped:
            for resume in (False, True):
                with self.assertRaises(KeyboardInterrupt):
                    pilot.run("fixture", self.out, resume=resume, env={})
            self.assertEqual(2, stopped.call_count)
        records = [compute.read_envelope(path) for path in (self.out / "records").glob("*.json")]
        saved = next(record for record in records if record["attempts"])
        self.assertIsNone(saved["candidate"])
        self.assertEqual([1, 2], [a["number"] for a in saved["attempts"]])
        for attempt in saved["attempts"]:
            self.assertEqual("interrupted-attempt-outcome-unknown", attempt["error"])
            self.assertIsNone(attempt["usage"])
            self.assertIsNone(attempt["elapsed_seconds"])
        with patch.object(pilot_worker, "invoke", wraps=pilot_worker.invoke) as calls:
            resumed = pilot.run("fixture", self.out, resume=True, env={})
            self.assertEqual(16, calls.call_count)  # eight other complete jobs, one exhausted
        self.assertEqual(8, resumed["metrics"]["jobs"]["completed"])
        self.assertEqual(1, resumed["metrics"]["jobs"]["failed"])
        self.assertEqual(18, resumed["metrics"]["execution"]["attempts"])
        with patch.object(pilot_worker, "invoke", side_effect=AssertionError("retry over bound")):
            pilot.run("fixture", self.out, resume=True, env={})

    def test_preflight_block_is_inspectable_and_safe_to_resume(self):
        with patch.object(pilot, "prepare_snapshot", side_effect=ValueError("private upstream error")), \
                patch.object(pilot_worker, "invoke", side_effect=AssertionError("unexpected call")):
            result = pilot.run("real", self.out, env={})
        self.assertEqual("blocked-before-inference", result["status"])
        self.assertEqual((9, 0, 0), (result["jobs_blocked"], result["jobs_attempted"], result["jobs_failed"]))
        self.assertIsNone(result["sources"])
        self.assertNotIn("private upstream error", (self.out / "report.md").read_text())
        self.assertFalse((self.out / "records").exists())
        self.assertFalse((self.out / "evidence.json").exists())
        # Source access recovers; without any configured endpoint zero model calls are made.
        self.primary.side_effect = None
        self.primary.return_value = civic_evidence.FIXTURE.read_bytes()
        with patch.object(pilot_worker, "invoke", side_effect=AssertionError("unexpected call")):
            result = pilot.run("real", self.out, resume=True, env={})
        self.assertEqual("blocked-before-inference", result["status"])
        self.assertEqual(9, result["jobs_blocked"])
        self.assertEqual("not-retrieved", result["sources"]["gdc"]["status"])
        self.assertTrue((self.out / "evidence.json").exists())
        self.assertFalse((self.out / "records").exists())

    def test_import_unavailable_native_bundle_preserves_failed_query_without_dispatch(self):
        native = evidence.civic_evidence.retrieve_fixture(None)
        native.update(synthetic=False, status="unavailable", failure="rate-limit",
                      failed_variables={"first": native["page_size"], "after": None},
                      pages=[], records=[])
        native["sha256"] = compute.digest({k: v for k, v in native.items() if k != "sha256"})
        path = self.tmp / "unavailable.json"
        compute.atomic_json(path, native)
        with patch.object(evidence, "fetch", side_effect=AssertionError("unrequested fallback")), \
                patch.object(pilot_worker, "invoke", side_effect=AssertionError("unexpected dispatch")):
            report = pilot.run("real", self.out, supplied=path, env={})
        self.assertEqual("blocked-before-inference", report["status"])
        self.assertEqual(0, report["jobs_attempted"])
        self.assertEqual("rate-limit", report["sources"]["civic"]["agent1"]["failure"])
        self.assertEqual(native["sha256"], report["sources"]["civic"]["agent1"]["sha256"])
        self.assertFalse(report["sources"]["civic"]["nightly_fallback_attempted"])

    def test_prohibits_public_output_and_overlapping_writers(self):
        targets = ["results/egfr-luad-evidence-coverage/test-output", "catalog/test-output",
                   "explorer/test-output", "_site/test-output", ".compute"]
        for target in targets:
            with self.subTest(target=target), self.assertRaisesRegex(ValueError, "run directory"):
                pilot.run("fixture", pilot.ROOT / target, env={})
        with pilot.exclusive(self.out):
            with self.assertRaisesRegex(ValueError, "another pilot process|Another pilot process"):
                pilot.run("fixture", self.out, resume=True, env={})

    def test_real_mode_wiring_against_mock_source_and_mock_inference_server(self):
        """The real-mode code path is executed in isolation; these are NOT real data or model results."""
        seen = []
        unavailable = []

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_POST(self):
                req = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
                prompt = json.loads(req["messages"][1]["content"])
                stage = "critique" if "candidate" in prompt else "investigate"
                seen.append((stage, self.path, req["model"], self.headers.get("Authorization")))
                unavailable.extend(r for r in prompt["packet"]["records"] if r["source_id"] != "civic")
                answer = pilot_worker._fixture_response(stage, prompt["packet"])
                if stage == "investigate":
                    answer["assumptions"] = ["This is a synthetic transport-stub test."]
                response = {"model": "stub-build-" + stage, "usage": {"prompt_tokens": 12,
                            "completion_tokens": 8, "total_tokens": 20}, "choices": [
                                {"message": {"content": json.dumps(answer)}, "finish_reason": "stop"}]}
                self.send_response(200)
                self.end_headers()
                self.wfile.write(json.dumps(response).encode())

        with HTTPServer(("127.0.0.1", 0), Handler) as local:
            thread = Thread(target=local.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True)
            thread.start()
            try:
                url = f"http://127.0.0.1:{local.server_port}"
                env = {}
                for role in ("INVESTIGATOR", "CRITIC"):
                    prefix = "OPEN_CANCER_" + role
                    env.update({prefix + "_MODEL": "stub-" + role, prefix + "_KEY": "test-key-" + role,
                                prefix + "_URL": url + "/" + role + "/chat/completions",
                                prefix + "_INPUT_USD_PER_MILLION": "1", prefix + "_OUTPUT_USD_PER_MILLION": "2"})
                self.primary.side_effect = None
                self.primary.return_value = civic_evidence.FIXTURE.read_bytes()
                report = pilot.run("real", self.out, env=env)
                self.primary.assert_called_once()
                self.tsv.assert_not_called()
            finally:
                local.shutdown()
                thread.join()
        self.assertEqual(18, len(seen))
        self.assertFalse(unavailable)
        for stage, path, model, key in seen:
            role = "INVESTIGATOR" if stage == "investigate" else "CRITIC"
            self.assertEqual((f"/{role}/chat/completions", "stub-" + role, "Bearer test-key-" + role),
                             (path, model, key))
        self.assertEqual(9, report["metrics"]["jobs"]["completed"])
        costs = report["metrics"]["execution"]
        self.assertEqual(18, costs["endpoint_attempts"])
        self.assertEqual(360, costs["usage"]["total_tokens"])
        self.assertAlmostEqual(0.000504, costs["estimated_cost_usd"])
        self.assertFalse(report["synthetic"])  # tests the mode marker, not proof of live evidence
        for path in self.out.rglob("*"):
            if path.is_file():
                self.assertNotIn(b"test-key-", path.read_bytes())


if __name__ == "__main__":
    unittest.main()
