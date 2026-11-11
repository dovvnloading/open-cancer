"""Regression checks for public metadata discovery without network calls."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts import refresh_catalog as refresh


class RefreshCatalogTests(unittest.TestCase):
    def test_gdc_project_modalities_follow_upstream_categories(self):
        project = {
            "project_id": "TCGA-BRCA", "name": "Breast Invasive Carcinoma",
            "program": {"name": "TCGA"},
            "summary": {"data_categories": [{"data_category": "Transcriptome Profiling"}]},
        }
        with patch.object(refresh, "get_json", side_effect=[
            {"data": project}, {"status": "OK", "data_release": "Data Release 46.0"},
        ]) as request:
            record, = refresh.gdc_records("TCGA-BRCA", "2026-01-01")
        self.assertEqual(record["modalities"], ["transcriptomics"])
        self.assertEqual(record["disease"], {"name": "Breast Invasive Carcinoma"})
        self.assertEqual(record["access"]["level"], "mixed")
        self.assertEqual(record["retrieval_url"], request.call_args_list[0].args[0])

    def test_truncated_gdc_list_and_wrong_accession_are_rejected(self):
        with patch.object(refresh, "get_json", return_value={
            "data": {"pagination": {"total": 33}, "hits": []}
        }), self.assertRaisesRegex(ValueError, "incomplete"):
            refresh.gdc_records(None, "2026-01-01")
        with patch.object(refresh, "get_json", return_value={"data": {"project_id": "TCGA-LUAD"}}), \
                self.assertRaisesRegex(ValueError, "different project"):
            refresh.gdc_records("TCGA-BRCA", "2026-01-01")

    def test_cbio_scope_public_status_and_modalities(self):
        study = {
            "studyId": "brca_tcga_pan_can_atlas_2018", "publicStudy": True,
            "name": "Breast Invasive Carcinoma (TCGA, PanCancer Atlas)", "pmid": "123,456,123",
            "cancerType": {"name": "Breast Cancer"}, "importDate": "2018-12-01",
            "cnaSampleCount": 1, "sequencedSampleCount": 0, "rppaSampleCount": 0,
            "massSpectrometrySampleCount": 0, "mrnaRnaSeqV2SampleCount": 0,
            "mrnaRnaSeqSampleCount": 0, "mrnaMicroarraySampleCount": 0,
        }
        with patch.object(refresh, "get_json", return_value=[study, {"studyId": "other"}]), \
                patch.object(refresh, "cbio_total", return_value=2):
            record, = refresh.cbio_records(None, "2026-01-01")
        self.assertEqual(record["modalities"], ["genomics"])
        self.assertEqual(record["reuse"]["status"], "verify-upstream")
        self.assertEqual(record["citation_url"], "https://pubmed.ncbi.nlm.nih.gov/?term=123%5Buid%5D+OR+456%5Buid%5D")
        self.assertIn(study["importDate"], record["release"]["label"])
        with self.assertRaisesRegex(ValueError, "public"):
            refresh.cbio_record(study | {"publicStudy": False}, "2026-01-01")

    def test_incomplete_cbio_collection_and_gdc_warnings_fail_closed(self):
        with patch.object(refresh, "get_json", return_value=[{"studyId": "brca"}]), \
                patch.object(refresh, "cbio_total", return_value=600), \
                self.assertRaisesRegex(ValueError, "incomplete"):
            refresh.cbio_records(None, "2026-01-01")
        payload = {"warnings": {"fields": "unrecognized field"}, "data": {"project_id": "TCGA-LUAD"}}
        with patch.object(refresh, "get_json", return_value=payload), \
                self.assertRaisesRegex(ValueError, "warnings"):
            refresh.gdc_records("TCGA-LUAD", "2026-01-01")

    def test_civic_checks_named_export_header(self):
        filename = "nightly-MolecularProfileSummaries.tsv"
        with patch.object(refresh, "civic_header",
                          return_value=["molecular_profile_id", "last_review_date", "is_flagged"]) as header:
            record, = refresh.civic_records(filename, "2026-01-01")
        self.assertEqual(record["accession"]["value"], filename)
        self.assertEqual(record["retrieval_url"], header.call_args.args[0])
        with patch.object(refresh, "civic_header", return_value=["<html>"]), \
                self.assertRaisesRegex(ValueError, "unexpected columns"):
            refresh.civic_records(filename, "2026-01-01")

    def test_rejects_duplicate_results_and_bad_lookup(self):
        with patch.dict(refresh.CONNECTORS, {"gdc": lambda *args: [{"id": "same"}, {"id": "same"}]}), \
                self.assertRaisesRegex(ValueError, "duplicate"):
            refresh.discover("gdc")
        with self.assertRaisesRegex(ValueError, "no matching"):
            refresh.discover("civic", "unknown.tsv")
        with self.assertRaisesRegex(ValueError, "select --source"):
            refresh.discover(accession="TCGA-BRCA")

    def test_staging_preserves_existing_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "output"
            refresh.save_new_directory(path, {"one": {"id": "one"}})
            with self.assertRaisesRegex(ValueError, "new output"):
                refresh.save_new_directory(path, {"one": {"id": "different"}})
            self.assertEqual(json.loads((path / "one.json").read_text()), {"id": "one"})
