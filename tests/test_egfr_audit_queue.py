"""Offline EGFR/LUAD QA queue tests: metadata arithmetic, not human review."""
from __future__ import annotations

from collections import Counter
import copy
import csv
import hashlib
import io
import json
import unittest

from scripts.build_review_packet import ROOT, RESULT
from scripts.egfr_audit_queue import make_civic_triage, make_cross_source_units, text_only_flags


class EGFRSourceAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = ROOT / RESULT
        cls.manifest = (directory / "manifest.json").read_bytes()
        cls.coverage = json.loads((directory / "coverage.json").read_bytes())
        cls.civic = list(csv.DictReader(io.StringIO(
            (directory / "civic_evidence_ids.csv").read_text(encoding="utf-8"))))
        cls.cbiorows = list(csv.DictReader(io.StringIO(
            (directory / "cbio_protein_changes.csv").read_text(encoding="utf-8"))))

    def test_full_selected_civic_inventory_not_an_upstream_validation(self):
        content, totals = make_civic_triage(self.civic)
        records = list(csv.DictReader(io.StringIO(content.decode("utf-8"))))
        self.assertEqual(len(records), 56)
        self.assertEqual(totals, {"Oncogenic": 2, "Predictive": 52, "Prognostic": 2})
        self.assertEqual(len({r["record_id"] for r in records}), 56)
        self.assertEqual([int(r["record_id"]) for r in records],
                         sorted(int(r["record_id"]) for r in records))
        self.assertTrue(all(r["upstream_record_checked"] == "no" and
                            r["scientific_interpretation_checked"] == "no"
                            for r in records))
        by_id = {r["record_id"]: r for r in records}
        self.assertIn("multiple-terms-AND-OR", by_id["1396"]["text_only_flags"])
        self.assertIn("fusion-or-partner-profile", by_id["5909"]["text_only_flags"])
        self.assertIn("categorical-unspecified-mutation", by_id["1701"]["text_only_flags"])
        self.assertEqual(by_id["879"]["text_only_flags"], "none-detected")
        self.assertTrue(all(r["record_url"].endswith("/" + r["record_id"])
                            for r in records))

    def test_source_profile_flags_are_heuristics_not_classifications(self):
        self.assertEqual(text_only_flags("EGFR L858R"), [])
        self.assertEqual(text_only_flags("EGFR Mutation"), ["categorical-unspecified-mutation"])
        self.assertIn("multiple-terms-AND-OR",
                      text_only_flags("EGFR T790M AND EGFR C797S"))
        self.assertIn("fusion-or-partner-profile", text_only_flags("KIF5B::EGFR Fusion"))
        with self.assertRaisesRegex(ValueError, "EGFR-containing"):
            text_only_flags("KRAS G12C")

    def test_original_snapshot_integrity_and_distinct_units(self):
        data = make_cross_source_units(
            self.coverage, self.manifest, len(self.civic),
            sum(int(r["mutation_call_rows"]) for r in self.cbiorows))
        audit = json.loads(data)
        self.assertEqual(audit["original_manifest_sha256"],
                         "9fce9beaaf873b5f0b0a14e338c3b372d8aa5d19c03be11b911ea1eb2b92931e")
        self.assertEqual(audit["sources"]["civic"]["count"], 56)
        self.assertEqual(audit["sources"]["cbioportal"]["count"], 86)
        self.assertEqual(audit["sources"]["gdc"]["occurrence_rows"], 93)
        self.assertEqual(audit["sources"]["gdc"]["distinct_case_ids_aggregated_in_memory"], 78)
        self.assertIsNone(audit["sources"]["gdc"]["project_profiled_denominator"])
        self.assertIsNone(audit["sources"]["cbioportal"]["molecular_profiled_denominator"])
        self.assertFalse(audit["cross_source_prevalence_calculable"])
        self.assertFalse(audit["independent_domain_review_completed"])
        self.assertFalse(audit["sources"]["civic"]["excluded_record_ids_available_in_original_artifact"])
        self.assertEqual(audit["sources"]["civic"]["upstream_item_verification"], "not-performed")
        self.assertNotIn("case.case_id", data.decode("utf-8"))

    def test_corruption_and_false_zero_rejected(self):
        bad = [dict(item) for item in self.civic]
        bad[1]["record_id"] = bad[0]["record_id"]
        with self.assertRaisesRegex(ValueError, "duplicated"):
            make_civic_triage(bad)
        bad = [dict(item) for item in self.civic]
        bad[0]["record_url"] = "https://other.org/not-civic"
        with self.assertRaisesRegex(ValueError, "record link"):
            make_civic_triage(bad)
        bad = [dict(item) for item in self.civic]
        bad[0]["doid"] = "DOID:123"
        with self.assertRaisesRegex(ValueError, "strict LUAD"):
            make_civic_triage(bad)
        with self.assertRaisesRegex(ValueError, "CIViC selected count"):
            make_cross_source_units(self.coverage, self.manifest, 55, 86)
        with self.assertRaisesRegex(ValueError, "cBioPortal aggregate"):
            make_cross_source_units(self.coverage, self.manifest, 56, 85)
        coverage = copy.deepcopy(self.coverage)
        coverage["gdc"]["egfr_occurrences"]["distinct_case_ids_with_recorded_occurrence"] = 94
        with self.assertRaisesRegex(ValueError, "Inconsistent GDC"):
            make_cross_source_units(coverage, self.manifest, 56, 86)
        coverage = copy.deepcopy(self.coverage)
        coverage["gdc"]["egfr_occurrences"]["status"] = "unavailable"
        coverage["gdc"]["egfr_occurrences"]["somatic_mutation_occurrence_rows"] = 0
        coverage["gdc"]["egfr_occurrences"]["distinct_case_ids_with_recorded_occurrence"] = None
        with self.assertRaisesRegex(ValueError, "remain null"):
            make_cross_source_units(coverage, self.manifest, 56, 86)


if __name__ == "__main__":
    unittest.main()
