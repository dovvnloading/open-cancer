"""Offline regression tests for the contribution contract."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from jsonschema import Draft202012Validator

from scripts.validate import ROOT, validate_repository


class ValidationTests(unittest.TestCase):
    def test_current_repository_passes(self):
        self.assertEqual([], validate_repository(ROOT))

    def test_source_requires_reuse_provenance(self):
        schema = json.loads((ROOT / "schemas/source.schema.json").read_text())
        source = json.loads((ROOT / "catalog/sources/civic.json").read_text())
        del source["reuse"]
        self.assertTrue(list(Draft202012Validator(schema).iter_errors(source)))

    def copy_catalog_to(self, root):
        for name in ("schemas", "catalog", "tasks"):
            shutil.copytree(ROOT / name, root / name)

    def test_unknown_source_reference_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.copy_catalog_to(root)
            task_path = root / "tasks/egfr-luad-evidence-coverage.json"
            task = json.loads(task_path.read_text())
            task["input_source_ids"] = ["not-a-real-source"]
            task_path.write_text(json.dumps(task))
            self.assertTrue(any("unknown source id" in err for err in validate_repository(root)))

    def test_artifact_checksum_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.copy_catalog_to(root)
            run = root / "results/egfr-luad-evidence-coverage/test-run"
            run.mkdir(parents=True)
            (run / "analyze.py").write_text("print('example')\n")
            (run / "summary.csv").write_text("source,count\ncivic,0\n")
            manifest = {
                "schema_version": 1,
                "task_id": "egfr-luad-evidence-coverage",
                "run_id": "test-run",
                "submitted_by": {"type": "agent", "name": "test-agent"},
                "created_at": "2026-11-10T00:00:00Z",
                "inputs": [{"source_id": "civic", "accessed_on": "2026-11-10",
                            "query_or_filter": "synthetic unit test only",
                            "url": "https://civicdb.org/"}],
                "execution": {"code_path": "analyze.py", "command": "python analyze.py",
                              "environment": "Python 3.12"},
                "artifacts": [{"path": "summary.csv", "sha256": "0" * 64,
                               "description": "Synthetic unit-test artifact"}],
                "findings": [], "limitations": "Synthetic data only; not a biomedical finding."
            }
            (run / "manifest.json").write_text(json.dumps(manifest))
            self.assertTrue(any("SHA-256 mismatch" in err for err in validate_repository(root)))


if __name__ == "__main__":
    unittest.main()
