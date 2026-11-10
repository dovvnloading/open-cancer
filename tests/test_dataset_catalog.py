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
        ids = {row["id"] for row in found}
        self.assertTrue({
            "cbioportal-luad-tcga-pancan-atlas-2018",
            "gdc-tcga-luad-project",
            "gdc-tcga-luad-egfr-ssm-occurrences",
        }.issubset(ids))
        self.assertTrue(all(row["disease"]["ontology_id"] == "DOID:3910"
                            for row in found))

    def test_source_and_search(self):
        civic = search_datasets(self.data, source="civic")
        self.assertTrue(civic)
        self.assertTrue(all(row["source_id"] == "civic" for row in civic))
        lung = search_datasets(self.data, text="lung")
        self.assertTrue(lung)
        self.assertTrue(all("lung" in " ".join(
            [row.get("id", ""), row.get("title", ""), row.get("description", ""),
             row.get("cohort_definition", ""), row.get("disease", {}).get("name", "")]
        ).casefold() for row in lung))

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
        self.assertTrue(all(isinstance(r["release"]["mutable"], bool)
                            for r in self.data))


if __name__ == "__main__":
    unittest.main()
