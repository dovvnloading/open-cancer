"""Offline end-to-end checks through the coordinator and real subprocess protocol."""
from contextlib import redirect_stdout, redirect_stderr
import io
import json
from pathlib import Path
import shlex
import sys
import tempfile
import unittest
from unittest.mock import patch

from scripts import civic_evidence as civic
from scripts import compute_fanout as compute

DEMO = shlex.join([sys.executable, str(compute.ROOT / "scripts/demo_compute_worker.py")])


class GroundedFanoutTest(unittest.TestCase):
    def setUp(self):
        staging = compute.ROOT / ".compute"
        staging.mkdir(exist_ok=True)
        tmp = tempfile.TemporaryDirectory(dir=staging)
        self.addCleanup(tmp.cleanup)
        self.directory = Path(tmp.name)
        self.snapshot = self.directory / "civic.json"
        self.plan = self.directory / "jobs.jsonl"
        self.out = self.directory / "run"
        self.assertEqual(0, self.cli("retrieve", "--source", "civic", "--fixture", "--output", self.snapshot))
        spec = compute.ROOT / "workloads/civic-egfr-luad-annotations.json"
        self.assertEqual(0, self.cli("plan", "--spec", spec, "--output", self.plan))

    @staticmethod
    def cli(*args):
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            return compute.main([str(x) for x in args])

    def batch(self, *extra, investigator=DEMO, critic=DEMO):
        return self.cli("run", "--jobs", self.plan, "--evidence", self.snapshot, "--output",
                        self.out, "--investigator", investigator, "--critic", critic, "--retries", "0",
                        "--timeout", "10", *extra)

    def test_fixture_citation_lineage_through_both_stages_and_resume(self):
        with patch.object(civic, "_fetch_page", side_effect=AssertionError("Live API called")):
            self.assertEqual(0, self.batch())
            config = compute.read_json(self.out / "run-config.json")
            bundle = config["evidence_bundle"]
            self.assertEqual(compute.load_evidence(self.snapshot), bundle)
            before = {p: p.read_bytes() for p in (self.out / "records").glob("*.json")}
            by_profile = {}
            for path in before:
                body = compute.read_envelope(path)
                self.assertEqual(bundle["sha256"], body["evidence_sha256"])
                self.assertEqual("not-published", body["publication_status"])
                self.assertTrue(body["candidate"]["synthetic"])
                by_profile[body["job"]["dimensions"]["alteration"]] = body
                index = {r["record_id"]: r for r in civic.records_for_job(bundle, body["job"])}
                for answer in (body["candidate"], body["critique"]):
                    for citation in answer["evidence"]:
                        record = index[citation["record_id"]]
                        self.assertEqual(record["sha256"], citation["sha256"])
                        self.assertIn("response_sha256", bundle["pages"][record["page"]])
            contested = by_profile["EGFR L858R"]
            self.assertEqual("inconsistent", contested["critique"]["verdict"])
            self.assertEqual("contradicts", contested["critique"]["evidence"][0]["relation"])
            missing = by_profile["EGFR exon 20 insertion"]
            self.assertEqual({"investigate": 0, "critique": 0}, missing["attempts"])
            self.assertEqual("insufficient-evidence", missing["triage"])
            self.assertEqual(0, self.batch("--resume"))
            self.assertEqual(before, {p: p.read_bytes() for p in before})

            with patch.dict("os.environ", OPEN_CANCER_CRITIC_MODEL="changed-model"):
                self.assertEqual(2, self.batch("--resume"))
            self.assertEqual(2, self.cli("retrieve", "--source", "civic", "--fixture",
                                         "--output", self.snapshot))
            self.assertEqual(bundle, compute.load_evidence(self.snapshot))

    def test_retrieval_failure_cached_explicitly_no_subprocess_or_fabricated_zero(self):
        self.snapshot.unlink()
        with patch.object(civic, "_fetch_page", return_value=b'{"errors":[{"message":"private"}]}'):
            self.assertEqual(1, self.cli("retrieve", "--source", "civic", "--output", self.snapshot))
        bundle = compute.load_evidence(self.snapshot)
        self.assertEqual("unavailable", bundle["status"])
        self.assertNotIn("private", self.snapshot.read_text())
        with patch.object(compute, "invoke", side_effect=AssertionError("Worker was called")):
            self.assertEqual(0, self.batch())
        summary = compute.read_json(self.out / "summary.json")
        self.assertEqual({"insufficient-evidence": 3}, summary["triage"])
        for p in (self.out / "records").glob("*.json"):
            body = compute.read_envelope(p)
            self.assertEqual("insufficient-evidence", body["candidate"]["outcome"])
            self.assertEqual([], body["candidate"]["evidence"])

    def test_changed_snapshot_tampered_citation_and_new_permission_policy_refuse_resume(self):
        self.assertEqual(0, self.batch())
        bundle = compute.load_evidence(self.snapshot)
        other = dict(bundle, retrieved_at="2026-02-01T00:00:00Z")
        other["sha256"] = compute.digest({k: v for k, v in other.items() if k != "sha256"})
        compute.atomic_json(self.snapshot, other)
        self.assertEqual(2, self.batch("--resume"))
        compute.atomic_json(self.snapshot, bundle)
        catalog = dict(compute.civic_catalog(), access="controlled")
        with patch.object(compute, "civic_catalog", return_value=catalog):
            self.assertEqual(2, self.batch("--resume"))
        for p in (self.out / "records").glob("*.json"):
            body = compute.read_envelope(p)
            if body["candidate"]["evidence"]:
                body["candidate"]["evidence"][0]["sha256"] = "0" * 64
                compute.atomic_json(p, compute.envelope(body))
                break
        self.assertEqual(2, self.batch("--resume"))

    def test_arbitrary_subprocess_cannot_forge_investigator_or_critic_output(self):
        jobs = [j for j in compute.load_jobs(self.plan) if j["dimensions"]["alteration"] == "EGFR L858R"]
        self.plan.write_text(json.dumps(jobs[0]) + "\n")
        adapter = self.directory / "worker.py"
        adapter.write_text(
            "import json, sys\n"
            f"sys.path.insert(0, {str(compute.ROOT)!r})\n"
            "from scripts.grounded_claims import fixture_response\n"
            "r = json.load(sys.stdin)\n"
            "a = fixture_response(r['stage'], r['evidence_bundle'], r['job'], r.get('candidate'))\n"
            "a['evidence'][0]['record_id'] = 'EID1999999999'\n"
            "print(json.dumps(a))\n")
        bad = shlex.join([sys.executable, str(adapter)])
        for stage in ("investigate", "critique"):
            with self.subTest(stage=stage):
                self.out = self.directory / stage
                self.assertEqual(1, self.batch(investigator=bad if stage == "investigate" else DEMO,
                                               critic=bad if stage == "critique" else DEMO))
                body = compute.read_envelope(next((self.out / "records").glob("*.json")))
                self.assertEqual("failed", body["status"])
                self.assertIn(stage + ": invalid-worker-response", body["error"])
                self.assertNotIn("1999999999", json.dumps(body))

    def test_reject_public_output_and_unrelated_plan_before_workers(self):
        for directory in ("results", "explorer", "docs", "_site", "catalog"):
            with self.subTest(directory=directory), self.assertRaisesRegex(ValueError, "under .compute"):
                compute.private_output(compute.ROOT / directory / "snapshot.json")
        jobs = compute.plan(compute.read_json(compute.ROOT / "workloads/egfr-luad-compute.json"), 9)
        self.plan.write_text(json.dumps(jobs[0]) + "\n")
        self.assertEqual(2, self.batch())
        self.assertFalse(self.out.exists())

    def test_legacy_custom_worker_cannot_assert_or_invent_record_links(self):
        job = compute.load_jobs(self.plan)[0]
        payload = {"summary": "Source metadata supports a research question only.",
                   "synthetic": False, "limitations": ["No records were supplied."],
                   "evidence": [{"source_id": "civic", "url": "https://civicdb.org/",
                                 "relation": "supports"}]}
        for url, relation in (("https://civicdb.org/", "supports"), ("https://example.org", "context")):
            payload["evidence"][0].update(url=url, relation=relation)
            with self.assertRaisesRegex(ValueError, "directory only"):
                compute.check_directory_response("investigate", payload, job)
        for verdict in ("inconsistent", "no-internal-issues-found"):
            answer = compute.check_directory_response("critique", {
                "verdict": verdict, "issues": ["The interpretation lacks actual records."]}, job)
            self.assertEqual("inconsistent" if verdict == "inconsistent" else "insufficient-evidence",
                             answer["verdict"])
            self.assertIn("The interpretation lacks actual records.", answer["issues"])
