"""Offline dataset catalog discovery and integrity tests."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from scripts.search_catalog import search_datasets
from scripts.validate import ROOT, validate_repository


class DatasetCatalogTests(unittest.TestCase):
    def setUp(self):
        self.data = [json.loads(path.read_text())
                     for path in sorted((ROOT / "catalog/datasets").glob("*.json"))]

    def test_exact_disease_filter(self):
        found = search_datasets(self.data, disease="DOID:3910")
        self.assertEqual([r["id"] for r in found], [
            "cbioportal-luad-tcga-pancan-atlas-2018", "gdc-tcga-luad-project"
        ])

    def test_source_and_search(self):
        self.assertEqual(len(search_datasets(self.data, source="civic")), 1)
        self.assertEqual(len(search_datasets(self.data, text="lung")), 2)

    def test_unknown_source_fails_validation(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for name in ("schemas", "catalog", "tasks"):
                shutil.copytree(ROOT / name, root / name)
            file = root / "catalog/datasets/civic-nightly-clinical-evidence.json"
            data = json.loads(file.read_text())
            data["source_id"] = "invented-source"
            file.write_text(json.dumps(data))
            errors = validate_repository(root)
            self.assertTrue(any("dataset references unknown source" in e for e in errors))

    def test_mutable_release_is_explicit(self):
        self.assertTrue(all(r["release"]["mutable"] is True for r in self.data))


if __name__ == "__main__":
    unittest.main()
