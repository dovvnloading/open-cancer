"""Offline tests for native Agent 1 bundles, job-specific citations and the bounded TSV fallback."""
from copy import deepcopy
import csv
from hashlib import sha256
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase, mock

from scripts import civic_evidence as agent1
from scripts.compute_fanout import digest
from scripts import pilot_evidence as evidence


def job(alteration):
    return {"task_id": evidence.TASK, "source_ids": ["civic", "cbioportal", "gdc"],
            "dimensions": {"alteration": alteration}}


def fixture_page():
    return agent1.FIXTURE.read_bytes()


def unavailable(reason):
    with mock.patch.object(agent1, "_fetch_page", side_effect=agent1._Unavailable(reason)):
        return agent1.retrieve(None)


def nightly_rows():
    path = evidence.ROOT / "tests/fixtures/egfr_workload_civic_synthetic.tsv"
    rows = list(csv.DictReader(io.StringIO(path.read_text()), delimiter="\t"))
    for row in rows:
        row["evidence_civic_url"] = agent1.LINK + row["evidence_id"]
    return rows


def pack_rows(rows):
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buf.getvalue().encode()


def resign(s):
    s["snapshot_sha256"] = digest({k: v for k, v in s.items() if k != "snapshot_sha256"})


class EvidenceTests(TestCase):
    def test_native_fixture_through_agent1_without_network(self):
        with mock.patch.object(evidence, "fetch", side_effect=AssertionError), \
             mock.patch.object(agent1, "_fetch_page", side_effect=AssertionError):
            s = evidence.prepare_snapshot("fixture")
            self.assertEqual(s, evidence.from_agent1_bundle(agent1.retrieve_fixture(None), "fixture"))
            evidence.validate_snapshot(s, "fixture")
        c = s["sources"]["civic"]
        self.assertEqual(c["connector"], "agent1-graphql")
        self.assertEqual(c["agent1"]["pages"][0]["response_sha256"], sha256(fixture_page()).hexdigest())
        self.assertEqual(c["sha256"], s["agent1_bundle"]["sha256"])
        self.assertEqual(c["matching_by_alteration"], dict(zip(evidence.ALTERATIONS, [2, 1, 0])))
        original = {r["record_id"]: r for r in s["agent1_bundle"]["records"]}
        for alt in evidence.ALTERATIONS:
            pkt = evidence.evidence_for_job(s, job(alt))
            self.assertEqual(pkt["snapshot_sha256"], s["snapshot_sha256"])
            self.assertNotIn("agent1_bundle", pkt)
            self.assertNotIn("records", pkt["sources"]["civic"]["agent1"])
            for r in pkt["records"]:
                self.assertTrue(r["url"].startswith("https://fixture.invalid/"))
                self.assertEqual(r["record_sha256"], original[r["record_id"]]["sha256"])
                self.assertEqual(r["molecular_profile"].casefold(), alt.casefold())
            for missing in ("gdc", "cbioportal"):
                self.assertEqual(pkt["sources"][missing]["status"], "not-retrieved")
                self.assertNotIn("matching_records", pkt["sources"][missing])

    def test_primary_real_transport_validates_and_preserves_page_bytes(self):
        raw = fixture_page()
        with mock.patch.object(agent1, "_fetch_page", return_value=raw) as graphql, \
             mock.patch.object(evidence, "fetch", side_effect=AssertionError):
            s = evidence.prepare_snapshot("real")
        graphql.assert_called_once_with({"first": 20, "after": None})
        c = s["sources"]["civic"]
        self.assertEqual(c["retrieval_status"], "complete")
        self.assertEqual(c["url"], agent1.ENDPOINT)
        self.assertEqual(c["agent1"]["pages"][0]["response_sha256"], sha256(raw).hexdigest())
        row = s["records"][0]
        self.assertEqual(row["url"], agent1.LINK + row["record_id"][3:])
        with self.assertRaisesRegex(ValueError, "mode"):
            evidence.validate_snapshot(s, "fixture")
        with self.assertRaisesRegex(ValueError, "mode"):
            evidence.from_agent1_bundle(s["agent1_bundle"], "fixture")
        with self.assertRaisesRegex(ValueError, "mode"):
            evidence.from_agent1_bundle(agent1.retrieve_fixture(None), "real")

    def test_truncated_native_input_remains_partial_without_a_fallback(self):
        page = json.loads(fixture_page())
        data = page["data"]["evidenceItems"]
        data.update(totalCount=4)
        data["pageInfo"].update(hasNextPage=True)
        with mock.patch.object(agent1, "_fetch_page", return_value=json.dumps(page).encode()):
            bundle = agent1.retrieve(None, max_pages=1, page_size=3)
        with mock.patch.object(agent1, "retrieve", return_value=bundle), \
             mock.patch.object(evidence, "fetch", side_effect=AssertionError):
            s = evidence.prepare_snapshot("real")
        p = evidence.evidence_for_job(s, job("EGFR exon 20 insertion"))
        self.assertEqual(p["sources"]["civic"]["retrieval_status"], "truncated")
        self.assertIn("Additional matches", p["sources"]["civic"]["reason"])
        self.assertEqual(p["records"], [])

    def test_unavailable_native_import_rate_limit_and_schema_drift_fail_closed(self):
        for reason in ("rate-limit", "network-error", "malformed-response", "source-changed"):
            bundle = unavailable(reason)
            with self.assertRaisesRegex(evidence.EvidenceUnavailable, reason) as err:
                evidence.from_agent1_bundle(bundle, "real")
            self.assertEqual(err.exception.sources["civic"]["agent1"]["sha256"], bundle["sha256"])
            if reason not in evidence.FALLBACK_ERRORS:
                with mock.patch.object(agent1, "retrieve", return_value=bundle), \
                     mock.patch.object(evidence, "fetch", side_effect=AssertionError):
                    with self.assertRaisesRegex(evidence.EvidenceUnavailable, reason):
                        evidence.prepare_snapshot("real")

    def test_nightly_fallback_preserves_graphql_failure_and_both_fingerprints(self):
        bundle = unavailable("network-error")
        rows = nightly_rows()
        for row in rows:
            row["patientId"] = "PRIVATE-SENTINEL"
            row["evidence_statement"] = "UNPROJECTED-NARRATIVE"
        rows[1]["is_flagged"] = "TRUE"
        rows[3]["doid"] = ""  # Do not infer strict ontology from its label.
        raw = pack_rows(rows)
        with mock.patch.object(agent1, "retrieve", return_value=bundle), \
             mock.patch.object(evidence, "fetch", return_value=raw) as fetch:
            s = evidence.prepare_snapshot("real")
        fetch.assert_called_once_with(evidence.CIVIC, retries=0)
        c = s["sources"]["civic"]
        self.assertEqual(c["connector"], "nightly-tsv-fallback")
        self.assertEqual(c["agent1"]["failed_variables"], {"first": 20, "after": None})
        self.assertEqual(c["agent1"]["sha256"], bundle["sha256"])
        self.assertEqual(c["sha256"], sha256(raw).hexdigest())
        self.assertEqual({r["record_id"] for r in s["records"]}, {"EID91", "EID93"})
        for sentinel in ("PRIVATE-SENTINEL", "UNPROJECTED-NARRATIVE"):
            self.assertNotIn(sentinel, str(s))
        evidence.validate_snapshot(s, "real")
        with self.assertRaisesRegex(ValueError, "mode"):
            evidence.validate_snapshot(s, "fixture")
        bad = deepcopy(s)
        bad["records"][0]["url"] += "0"
        resign(bad)
        with self.assertRaisesRegex(ValueError, "scope"):
            evidence.validate_snapshot(bad, "real")

    def test_no_raw_exports_or_fixtures_substituted_on_fallback_failure(self):
        bundle = unavailable("http-error")
        for invalid in (b"new header\n", pack_rows(nightly_rows()).replace(b"evidence_direction", b"x"),
                        pack_rows(nightly_rows()).replace(agent1.LINK.encode(), b"https://fixture.invalid/")):
            with mock.patch.object(agent1, "retrieve", return_value=bundle), \
                 mock.patch.object(evidence, "fetch", return_value=invalid):
                with self.assertRaisesRegex(evidence.EvidenceUnavailable, "Nightly fallback") as err:
                    evidence.prepare_snapshot("real")
                self.assertTrue(err.exception.sources["civic"]["nightly_fallback_attempted"])

    def test_catalog_policy_checked_before_either_network_request(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = root / "catalog/sources/civic.json"
            src.parent.mkdir(parents=True)
            src.write_bytes(agent1.CATALOG.read_bytes())
            rel = "catalog/datasets/civic-nightly-clinical-evidence.json"
            dataset = json.loads((evidence.ROOT / rel).read_text())
            dataset["reuse"]["status"] = "verify-upstream"
            target = root / rel
            target.parent.mkdir(parents=True)
            target.write_text(json.dumps(dataset))
            bundle = unavailable("http-error")
            with mock.patch.object(evidence, "ROOT", root), mock.patch.object(evidence, "fetch") as fetch:
                with mock.patch.object(agent1, "retrieve", return_value=bundle):
                    with self.assertRaises(evidence.EvidenceUnavailable):
                        evidence.prepare_snapshot("real")
                source = json.loads(src.read_bytes())
                source["reuse"]["status"] = "verify-upstream"
                src.write_text(json.dumps(source))
                with mock.patch.object(agent1, "_fetch_page", side_effect=AssertionError):
                    with self.assertRaisesRegex(ValueError, "permission"):
                        evidence.prepare_snapshot("real")
            fetch.assert_not_called()

    def test_cap_does_not_hide_matching_count_or_rare_direction(self):
        page = json.loads(fixture_page())
        data = page["data"]["evidenceItems"]
        data["nodes"] = [dict(data["nodes"][0], id=i) for i in range(1, 12)]
        data["nodes"][-1].update(evidenceDirection="DOES_NOT_SUPPORT")
        data.update(totalCount=11, pageInfo={"hasNextPage": False, "endCursor": "MTE"})
        with mock.patch.object(agent1, "_fetch_page", return_value=json.dumps(page).encode()):
            s = evidence.prepare_snapshot("real")
        pkt = evidence.evidence_for_job(s, job("EGFR L858R"))
        self.assertEqual(pkt["sources"]["civic"]["matching_records"], 11)
        self.assertEqual(pkt["sources"]["civic"]["included_records"], 8)
        self.assertIn("EID11", [r["record_id"] for r in pkt["records"]])

    def test_pinned_snapshot_binds_native_provenance_and_projected_values(self):
        base = evidence.prepare_snapshot("fixture")
        for mutate in (
            lambda s: s["records"][0].update(url="https://fixture.invalid/made-up"),
            lambda s: s["records"][0].update(patient_id="not permitted"),
            lambda s: s["sources"]["civic"].update(sha256="0" * 64),
            lambda s: s["agent1_bundle"]["pages"][0].update(response_sha256="0" * 64),
            lambda s: s["agent1_bundle"]["records"][0].update(page=1),
            lambda s: s["sources"]["gdc"].update(status="retrieved"),
        ):
            s = deepcopy(base)
            mutate(s)
            resign(s)
            with self.assertRaises(ValueError):
                evidence.validate_snapshot(s, "fixture")
        with self.assertRaisesRegex(ValueError, "scope"):
            evidence.evidence_for_job(base, job("EGFR"))
