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
    is_strict_luad, canonical_doid, audit_historical_source_queries,
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


    def test_historical_cancer_source_query_provenance_is_exact_and_nonclinical(self):
        record = audit_historical_source_queries(self.manifest, self.coverage)
        self.assertEqual(record["original_access_dates"], ["2026-11-10"])
        self.assertFalse(record["source_fingerprints_are_immutable_releases"])
        cbio = record["cbioportal"]
        self.assertEqual(cbio["study_id"], "luad_tcga_pan_can_atlas_2018")
        self.assertEqual(cbio["molecular_profile_id"], "luad_tcga_pan_can_atlas_2018_mutations")
        self.assertEqual(cbio["sample_list_id"], "luad_tcga_pan_can_atlas_2018_all")
        self.assertEqual(cbio["entrez_gene_id"], 1956)
        self.assertEqual(cbio["historical_mutation_call_rows"], 86)
        self.assertIsNone(cbio["egfr_profiled_sample_denominator"])
        self.assertFalse(cbio["profile_and_assay_coverage_verified"])
        self.assertIn("not-mutation-assayed-denominator", cbio["sample_list_scope"])
        gdc = record["gdc"]
        self.assertEqual(gdc["ensembl_gene_id"], "ENSG00000146648")
        self.assertEqual(gdc["historical_occurrence_rows"], 93)
        self.assertEqual(gdc["historical_distinct_case_ids_with_occurrence"], 78)
        self.assertEqual(gdc["requested_fields"], ["case.case_id", "ssm.ssm_id"])
        self.assertEqual(gdc["requested_from"], 0)
        self.assertEqual(gdc["requested_page_size"], 1000)
        self.assertFalse(gdc["consequence_type_filter_applied"])
        self.assertFalse(gdc["completeness_reverified_live"])
        self.assertIsNone(gdc["egfr_profiled_case_denominator"])
        self.assertEqual(gdc["filters"], {"op": "and", "content": [
            {"op": "in", "content": {"field": "ssm.consequence.transcript.gene.gene_id",
                                      "value": ["ENSG00000146648"]}},
            {"op": "in", "content": {"field": "case.project.project_id",
                                      "value": ["TCGA-LUAD"]}},
        ]})

    def test_historical_query_pin_rejects_tampered_filter_or_source_digest(self):
        altered = json.loads(json.dumps(self.coverage))
        altered["source_urls"]["gdc_occurrences"] = altered["source_urls"]["gdc_occurrences"].replace(
            "size=1000", "size=1"
        )
        with self.assertRaisesRegex(ValueError, "Historical source URLs"):
            audit_historical_source_queries(self.manifest, altered)
        altered = json.loads(json.dumps(self.coverage))
        altered["source_content_sha256"]["cbio_mutations"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "fingerprint mismatch"):
            audit_historical_source_queries(self.manifest, altered)
        with self.assertRaisesRegex(ValueError, "manifest has changed"):
            audit_historical_source_queries(self.manifest + b" ", self.coverage)

    def test_metadata_fingerprint_comparison_does_not_expose_sample_identifiers(self):
        cbio = b'{"studyId":"luad_tcga_pan_can_atlas_2018","sampleIds":["NEVER_OUTPUT_ID"]}'
        gdc = b'{"data":{"project_id":"TCGA-LUAD","case_ids":["NEVER_OUTPUT_GDC_ID"]}}'
        fixture = json.loads(json.dumps(self.coverage))
        fixture["source_content_sha256"]["cbio_study"] = hashlib.sha256(cbio).hexdigest()
        results = make_live_report(
            civic_raw=fake_tsv(), cbio_study_raw=cbio, gdc_project_raw=gdc,
            original_csv=self.original_csv, manifest_bytes=self.manifest,
            coverage=fixture, checked_at="2026-11-11T13:00:00Z",
        )
        summary = results["public-source-metadata-summary.json"]
        data = json.loads(summary)
        matches = data["metadata_response_byte_identity_with_original"]
        self.assertTrue(matches["cbio_study_metadata"])
        self.assertFalse(matches["gdc_project_metadata"])
        self.assertFalse(matches["civic"])
        self.assertEqual(data["historical_query_provenance"]["audit_type"],
                         "historical-query-provenance-only-not-new-row-retrieval")
        self.assertNotIn(b"NEVER_OUTPUT_ID", summary)
        self.assertNotIn(b"NEVER_OUTPUT_GDC_ID", summary)



if __name__ == "__main__":
    unittest.main()
