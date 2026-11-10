"""Offline fixtures guard against false matching and leaking patient identifiers."""
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

from scripts.egfr_luad_pilot import (
    CIVIC, CBIO_STUDY_ID, SOURCE_URLS, as_csv, civic_analysis,
    cbio_analysis, gdc_analysis, run,
)

HEADERS = ["molecular_profile", "evidence_id", "evidence_status", "disease",
           "doid", "evidence_type", "evidence_civic_url"]


def civic_fixture():
    rows = [
        ["EGFR L858R", "101", "accepted", "Lung Adenocarcinoma", "3910",
         "Predictive", "https://civicdb.org/evidence/101/summary"],
        ["EGFR T790M", "102", "ACCEPTED", "Non-Small Cell Lung Carcinoma",
         "3908", "Predictive", "https://civicdb.org/evidence/102/summary"],
        ["EGFR L858R", "103", "SUBMITTED", "Lung Adenocarcinoma",
         "3910", "Predictive", "https://civicdb.org/evidence/103/summary"],
        ["NOTEGFR L858R", "104", "accepted", "Lung Adenocarcinoma",
         "3910", "Predictive", "https://civicdb.org/evidence/104/summary"],
    ]
    f = io.StringIO(newline="")
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(HEADERS)
    writer.writerows(rows)
    return f.getvalue().encode()


def cbio_fixtures():
    study = {"studyId": CBIO_STUDY_ID, "name": "TCGA LUAD"}
    records = [
        {"entrezGeneId": 1956, "studyId": CBIO_STUDY_ID, "proteinChange": "L858R",
         "mutationType": "Missense_Mutation", "sampleId": "PRIVATE-NOT-PUBLISHED"},
        {"entrezGeneId": 1956, "studyId": CBIO_STUDY_ID, "proteinChange": "L858R",
         "mutationType": "Missense_Mutation", "patientId": "PRIVATE-NOT-PUBLISHED"},
    ]
    return json.dumps(study).encode(), json.dumps(records).encode()


class ResearchPilotTests(unittest.TestCase):
    def test_civic_narrow_doid_filter_and_status(self):
        counts, rows = civic_analysis(civic_fixture())
        self.assertEqual(counts["all_export_rows"], 4)
        self.assertEqual(counts["egfr_accepted_evidence_exact_luad"], 1)
        self.assertEqual(counts["egfr_accepted_other_lung_labels_excluded"], 1)
        self.assertEqual(rows[0]["record_id"], "101")

    def test_cbio_aggregates_without_record_identifiers(self):
        counts, rows = cbio_analysis(*cbio_fixtures())
        self.assertEqual(counts["egfr_mutation_call_rows"], 2)
        self.assertEqual(rows[0]["mutation_call_rows"], 2)
        self.assertNotIn("PRIVATE-NOT-PUBLISHED", json.dumps([counts, rows]))

    def test_fails_on_unexpected_gene(self):
        study, records = cbio_fixtures()
        changed = json.loads(records)
        changed[0]["entrezGeneId"] = 123
        with self.assertRaises(ValueError):
            cbio_analysis(study, json.dumps(changed).encode())

    def test_gdc_requires_project_id(self):
        with self.assertRaises(ValueError):
            gdc_analysis(json.dumps({"data": {"project_id": "TCGA-BRCA"}}).encode())

    def test_csv_formula_sanitization(self):
        self.assertIn("'+injection", as_csv([{"a": "+injection"}], ["a"]))

    def test_end_to_end_offline_no_patient_data(self):
        study, muts = cbio_fixtures()
        raw = {
            "civic": civic_fixture(), "cbio_study": study,
            "cbio_mutations": muts,
            "gdc": json.dumps({"data": {"project_id": "TCGA-LUAD", "summary": {"case_count": 585}}}).encode(),
        }
        self.assertEqual(set(raw), set(SOURCE_URLS))
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "pilot-20261110"
            run(raw, run_dir, "pilot-20261110", "2026-11-10", "2026-11-10T22:00:00Z")
            snapshot = json.loads((run_dir / "coverage.json").read_text())
            manifest = json.loads((run_dir / "manifest.json").read_text())
            self.assertEqual(snapshot["status"], "unreviewed")
            self.assertEqual(len(manifest["artifacts"]), 3)
            self.assertNotIn("PRIVATE-NOT-PUBLISHED", "".join(
                p.read_text() for p in run_dir.glob("*.csv")))
            self.assertEqual(manifest["findings"], [])


if __name__ == "__main__":
    unittest.main()
