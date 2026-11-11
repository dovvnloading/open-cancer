"""Synthetic offline fixtures for live public EGFR source-drift extraction.

These do NOT represent newly verified biological evidence or source acceptance.
"""
from __future__ import annotations

import csv
import io
import json
import hashlib
import unittest

from scripts.egfr_live_source_audit import (
    ROOT, RESULT, parse_civic_snapshot, compare_civic, make_live_report,
    is_strict_luad, canonical_doid,
)

HEADERS = ("evidence_id", "molecular_profile", "evidence_status",
           "disease", "doid", "evidence_type", "evidence_civic_url")


def fake_tsv():
    lines = [
        HEADERS,
        ("879", "EGFR L858R", "accepted", "Lung Adenocarcinoma",
         "DOID:3910", "Predictive", "https://civicdb.org/links/evidence_items/879"),
        ("880", "EGFR Exon 19 Deletion", "submitted", "Lung Adenocarcinoma",
         "DOID:3910", "Predictive", "https://civicdb.org/links/evidence_items/880"),
        ("999", "EGFR Mutation", "accepted", "Non-Small Cell Lung Cancer",
         "DOID:3908", "Predictive", "https://civicdb.org/links/evidence_items/999"),
        ("1000", "KRAS G12C", "accepted", "Lung Adenocarcinoma",
         "DOID:3910", "Predictive", "https://civicdb.org/links/evidence_items/1000"),
    ]
    out = io.StringIO(newline="")
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerows(lines)
    return out.getvalue().encode("utf-8")


class LiveCIViCAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = (RESULT / "manifest.json").read_bytes()
        cls.original_csv = (RESULT / "civic_evidence_ids.csv").read_bytes()
        cls.coverage = json.loads((RESULT / "coverage.json").read_bytes())

    def test_canonical_doid_namespace_does_not_create_fabricated_drift(self):
        self.assertEqual(canonical_doid("DOID:3910"), "3910")
        self.assertEqual(canonical_doid("doid: 3910"), "3910")
        self.assertEqual(canonical_doid("3910"), "3910")
        self.assertNotEqual(canonical_doid("3908"), canonical_doid("DOID:3910"))
        current = parse_civic_snapshot(fake_tsv())
        original = list(csv.DictReader(io.StringIO(self.original_csv.decode())))
        selected, _ = compare_civic(original, current)
        original_879 = next(x for x in selected if x["record_id"] == "879")
        self.assertEqual(original_879["original_doid"], "DOID:3910")
        self.assertEqual(original_879["current_doid"], "DOID:3910")
        self.assertEqual(original_879["changed_metadata_fields"], "none")
        # A genuinely different ontology term must still trigger review.
        amended = parse_civic_snapshot(fake_tsv())
        amended["all_by_id"]["879"]["doid"] = "DOID:3908"
        revised, _ = compare_civic(original, amended)
        updated = next(x for x in revised if x["record_id"] == "879")
        self.assertIn("doid", updated["changed_metadata_fields"])
        self.assertIn("current_inclusion", updated["changed_metadata_fields"])

    def test_disease_selectors_distinguish_exact_from_broad(self):
        self.assertTrue(is_strict_luad({"doid": "DOID:3910", "disease": "Lung Adenocarcinoma"}))
        self.assertTrue(is_strict_luad({"doid": "", "disease": "adenocarcinoma of the lung"}))
        self.assertFalse(is_strict_luad({"doid": "DOID:3908", "disease": "Non-Small Cell Lung Cancer"}))
        self.assertFalse(is_strict_luad({"doid": "DOID:3908", "disease": "Lung Adenocarcinoma"}))

    def test_classify_and_compare_preserves_historical_evidence_ids(self):
        current = parse_civic_snapshot(fake_tsv())
        self.assertEqual(current["total"], 4)
        self.assertEqual(set(current["strict"]), {"879"})
        self.assertEqual(set(current["broader"]), {"999"})
        original = list(csv.DictReader(io.StringIO(self.original_csv.decode())))
        selected, broader = compare_civic(original, current)
        self.assertEqual(len(selected), 56)
        self.assertEqual(len(broader), 1)
        by_id = {x["record_id"]: x for x in selected}
        self.assertEqual(by_id["879"]["changed_metadata_fields"], "none")
        self.assertIn("current_inclusion", by_id["880"]["changed_metadata_fields"])
        self.assertEqual(by_id["880"]["current_status"], "submitted")
        self.assertEqual(by_id["1396"]["current_status"], "not-present")
        self.assertEqual(broader[0]["source_record_independently_assessed"], "no")
        self.assertTrue(all(x["source_record_independently_assessed"] == "no"
                            for x in selected))
        self.assertTrue(all(x["record_url"].startswith(
            "https://civicdb.org/links/evidence_items/") for x in selected))

    def test_derived_report_tracks_version_drift_without_patient_data(self):
        raw = fake_tsv()
        cbio = b'{"studyId":"luad_tcga_pan_can_atlas_2018","name":"LUAD","sequencedSampleCount":400}'
        gdc = b'{"data":{"project_id":"TCGA-LUAD","summary":{"case_count":600}}}'
        packet = make_live_report(
            civic_raw=raw, cbio_study_raw=cbio, gdc_project_raw=gdc,
            original_csv=self.original_csv, manifest_bytes=self.manifest,
            coverage=self.coverage, checked_at="2026-11-11T13:00:00Z",
        )
        self.assertEqual(set(packet), {
            "civic-selected-record-drift.csv", "civic-broader-lung-candidates.csv",
            "public-source-metadata-summary.json", "SHA256SUMS",
        })
        data = json.loads(packet["public-source-metadata-summary.json"])
        self.assertEqual(data["civic"]["original_selected"], 56)
        self.assertEqual(data["civic"]["live_strict_selected"], 1)
        self.assertEqual(data["civic"]["live_broader_lung_candidates"], 1)
        self.assertTrue(data["civic"]["civic_snapshot_bytes_changed"])
        self.assertIsNone(data["cbioportal"]["egfr_mutation_profile_denominator"])
        self.assertIsNone(data["gdc"]["egfr_profiled_cohort_denominator"])
        self.assertEqual(data["cbioportal"]["study_level_sequenced_sample_count_if_provided"], 400)
        self.assertEqual(data["gdc"]["project_level_case_count_if_provided"], 600)
        self.assertFalse(data["cross_source_prevalence_calculable"])
        self.assertFalse(data["scientific_interpretation_assessed"])
        for line in packet["SHA256SUMS"].decode().splitlines():
            digest, name = line.split("  ", 1)
            self.assertEqual(digest, hashlib.sha256(packet[name]).hexdigest())

    def test_invalid_and_unsafe_sources_fail_without_false_zero(self):
        with self.assertRaisesRegex(ValueError, "columns changed"):
            parse_civic_snapshot(b"bad-header\n")
        rows = fake_tsv() + b"879\tEGFR Mutation\taccepted\tLung Adenocarcinoma\tDOID:3910\tPredictive\turl\n"
        with self.assertRaisesRegex(ValueError, "duplicate"):
            parse_civic_snapshot(rows)
        cbio_wrong = b'{"studyId":"someone_else"}'
        valid_gdc = b'{"data":{"project_id":"TCGA-LUAD"}}'
        args = dict(civic_raw=fake_tsv(), cbio_study_raw=cbio_wrong,
                    gdc_project_raw=valid_gdc, original_csv=self.original_csv,
                    manifest_bytes=self.manifest, coverage=self.coverage,
                    checked_at="2026-11-11T13:00:00Z")
        with self.assertRaisesRegex(ValueError, "Unexpected cBioPortal"):
            make_live_report(**args)
        args["cbio_study_raw"] = b'{"studyId":"luad_tcga_pan_can_atlas_2018"}'
        args["gdc_project_raw"] = b'{"data":{"project_id":"TCGA-BRCA"}}'
        with self.assertRaisesRegex(ValueError, "Unexpected GDC"):
            make_live_report(**args)
        args["gdc_project_raw"] = valid_gdc
        args["manifest_bytes"] = self.manifest + b" "
        with self.assertRaisesRegex(ValueError, "manifest has changed"):
            make_live_report(**args)


if __name__ == "__main__":
    unittest.main()
