"""Scientific-review schema and version-lock regression tests."""
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from jsonschema import Draft202012Validator, FormatChecker
from scripts.validate import ROOT, validate_repository


class ReviewContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema = json.loads((ROOT / "schemas/review.schema.json").read_text())

    def test_synthetic_reviews_validate_but_are_not_live_evidence(self):
        for path in (ROOT / "examples/reviews").glob("*.json"):
            item = json.loads(path.read_text())
            self.assertTrue(item["synthetic"])
            self.assertEqual([], list(Draft202012Validator(
                self.schema, format_checker=FormatChecker()).iter_errors(item)))

    def test_agent_cannot_attest_scientific_interpretation(self):
        item = json.loads((ROOT / "examples/reviews/synthetic-methods-example.json").read_text())
        item["reviewer"]["type"] = "agent"
        item["scope"] = "scientific-interpretation"
        item["decision"] = "scientific-assessment-completed"
        self.assertTrue(list(Draft202012Validator(self.schema).iter_errors(item)))

    def test_method_reproduction_must_have_command(self):
        item = json.loads((ROOT / "examples/reviews/synthetic-methods-example.json").read_text())
        item["reproduction"]["command"] = ""
        self.assertTrue(list(Draft202012Validator(self.schema).iter_errors(item)))

    def test_review_does_not_apply_to_changed_result_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for dirname in ("schemas", "catalog", "tasks", "examples"):
                shutil.copytree(ROOT / dirname, root / dirname)
            name = "synthetic-check"
            result_dir = root / "results/egfr-luad-evidence-coverage" / name
            result_dir.mkdir(parents=True)
            (result_dir / "analyze.py").write_text("print('synthetic')\n")
            (result_dir / "summary.csv").write_text("name,count\nexample,1\n")
            art_sha = hashlib.sha256((result_dir / "summary.csv").read_bytes()).hexdigest()
            manifest = {
                "schema_version": 1, "task_id": "egfr-luad-evidence-coverage",
                "run_id": name, "submitted_by": {"type": "human", "name": "fictional-test"},
                "created_at": "2026-11-10T00:00:00Z",
                "inputs": [{"source_id": "civic", "accessed_on": "2026-11-10",
                            "query_or_filter": "synthetic unit test fixture only",
                            "url": "https://civicdb.org/"}],
                "execution": {"code_path": "analyze.py", "command": "python analyze.py",
                              "environment": "Python 3.12"},
                "artifacts": [{"path": "summary.csv", "sha256": art_sha,
                               "description": "Fictional unit test results"}],
                "findings": [], "limitations": "All synthetic data and not real medical evidence."
            }
            manifest_path = result_dir / "manifest.json"
            manifest_path.write_text(json.dumps(manifest))
            review_dir = root / "reviews/egfr-luad-evidence-coverage" / name
            review_dir.mkdir(parents=True)
            review = json.loads((ROOT / "examples/reviews/synthetic-methods-example.json").read_text())
            review.update({
                "id": "mock-review",
                "task_id": "egfr-luad-evidence-coverage",
                "run_id": name,
                "reviewed_manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                "synthetic": False
            })
            (review_dir / "mock-review.json").write_text(json.dumps(review))
            self.assertEqual([], validate_repository(root))
            # Change the result manifest without updating the review attestation:
            manifest["limitations"] += " Changed after review."
            manifest_path.write_text(json.dumps(manifest))
            self.assertTrue(any("stale or mismatched reviewed manifest" in e
                                for e in validate_repository(root)))


if __name__ == "__main__":
    unittest.main()
