"""Offline investigator/critic contract tests using CIViC-shaped synthetic response data."""
import json
import unittest

from scripts import civic_evidence as civic
from scripts import grounded_claims as claims


class GroundedClaimsTests(unittest.TestCase):
    def setUp(self):
        self.bundle = civic.retrieve_fixture(None)
        self.job = {"task_id": civic.TASK, "source_ids": ["civic"],
                    "dimensions": {"alteration": "EGFR L858R"}}
        self.a, self.b, self.other = self.bundle["records"]
        self.candidate = claims.fixture_response("investigate", self.bundle, self.job)

    def alternate_snapshot(self, edit):
        """Exercise variant cases through the same parsing, linkage, and fingerprint checks."""
        response = json.loads(civic.FIXTURE.read_bytes())
        edit(response["data"]["evidenceItems"])
        return civic._snapshot(None, 1, 20, lambda _: json.dumps(response).encode(),
                               self.bundle["retrieved_at"], True)

    def check(self, candidate):
        return claims.check_grounded_response("investigate", candidate, self.bundle, self.job)

    def critique(self, answer):
        return claims.check_grounded_response("critique", answer, self.bundle, self.job, self.candidate)

    @staticmethod
    def cite(row, relation):
        return {"record_id": row["record_id"], "sha256": row["sha256"], "relation": relation}

    def test_critic_finds_counterexample_omitted_by_investigator(self):
        self.assertEqual([self.cite(self.a, "supports")], self.candidate["evidence"])
        answer = claims.fixture_response("critique", self.bundle, self.job, self.candidate)
        self.assertEqual({"verdict": "inconsistent", "issues": ["annotation-counterexample"],
                          "evidence": [self.cite(self.b, "contradicts")]}, answer)
        for bad in (
            dict(answer, verdict="no-internal-issues-found", issues=[]),
            dict(answer, verdict="insufficient-evidence", issues=["insufficient-record-evidence"]),
            dict(answer, evidence=[self.cite(self.a, "supports")]),
            dict(answer, evidence=[self.cite(self.other, "contradicts")]),
            dict(answer, issues=["The model's unrestricted commentary."]),
        ):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                self.critique(bad)

    def test_ids_hashes_relations_and_provenance_must_match_scope(self):
        valid = self.candidate["evidence"][0]
        for evidence in (
            [],
            [dict(valid, record_id="EID99")],
            [dict(valid, record_id=2000000001)],
            [dict(valid, sha256=self.b["sha256"])],
            [dict(valid, url=self.a["url"])],
            [self.cite(self.other, "supports")],
            [self.cite(self.a, "context")],
            [self.cite(self.b, "supports")],
            [self.cite(self.b, "contradicts")],
            [valid, valid],
        ):
            with self.subTest(evidence=evidence), self.assertRaises(ValueError):
                self.check(dict(self.candidate, evidence=evidence))
        both = dict(self.candidate, evidence=[valid, self.cite(self.b, "contradicts")])
        self.assertEqual(both, self.check(both))

    def test_abstain_on_missing_out_of_scope_and_unavailable_records(self):
        unavailable = self.alternate_snapshot(lambda page: page.update(nodes=None))
        cases = ((self.bundle, dict(self.job, dimensions={"alteration": "EGFR exon 20 insertion"})),
                 (self.bundle, dict(self.job, source_ids=["gdc"])),
                 (unavailable, self.job))
        for bundle, job in cases:
            with self.subTest(job=job, status=bundle["status"]):
                candidate = claims.fixture_response("investigate", bundle, job)
                self.assertEqual(claims.insufficient_response("investigate", bundle), candidate)
                response = claims.fixture_response("critique", bundle, job, candidate)
                self.assertEqual("insufficient-evidence", response["verdict"])
                for try_stage, answer in (("investigate", self.candidate), ("critique", dict(
                        response, verdict="no-internal-issues-found"))):
                    with self.assertRaises(ValueError):
                        claims.check_grounded_response(try_stage, answer, bundle, job, candidate)

    def test_other_evidence_types_provide_context_only(self):
        # Select an alternative categorical type for the second annotation.
        self.bundle = self.alternate_snapshot(lambda page: page["nodes"][1].update(
            evidenceType="DIAGNOSTIC", significance="POSITIVE"))
        other_type = self.bundle["records"][1]
        candidate = dict(self.candidate, evidence=[self.cite(self.a, "supports"),
                                                   self.cite(other_type, "context")])
        self.assertEqual(candidate, self.check(candidate))
        with self.assertRaises(ValueError):
            self.check(dict(candidate, evidence=[self.cite(self.a, "supports"),
                                                 self.cite(other_type, "contradicts")]))

    def test_only_complete_uncontested_annotations_can_pass_critique(self):
        self.candidate["observation"].update(field="significance", value="RESISTANCE")
        answer = {"verdict": "no-internal-issues-found", "issues": [],
                  "evidence": [self.cite(self.b, "supports")]}
        self.assertEqual(answer, self.critique(answer))
        with self.assertRaises(ValueError):
            self.critique(dict(answer, evidence=[]))
        def simulate_next_page(page):
            page["totalCount"] += 1
            page["pageInfo"]["hasNextPage"] = True

        self.bundle = self.alternate_snapshot(simulate_next_page)
        self.assertEqual(self.bundle["status"], "truncated")
        with self.assertRaises(ValueError):
            self.check(self.candidate)
        self.candidate["limitations"].append("bounded-retrieval")
        self.assertEqual(self.candidate, self.check(self.candidate))
        with self.assertRaises(ValueError):
            self.critique(answer)
        insufficient = {"verdict": "insufficient-evidence", "issues": ["bounded-retrieval"],
                        "evidence": []}
        self.assertEqual(insufficient, self.critique(insufficient))

    def test_closed_vocabulary_and_synthetic_marker(self):
        for changes in (
            {"outcome": "verified"}, {"synthetic": False}, {"synthetic": 1},
            {"summary": "Without an evidence-backed schema a model could put anything here."},
            {"limitations": ["https://example.invalid/private"]},
            {"limitations": ["bounded-retrieval", "record-annotations-only"]},
            {"observation": dict(self.candidate["observation"], field="summary")},
            {"observation": dict(self.candidate["observation"], value="verified")},
            {"observation": dict(self.candidate["observation"], evidence_type={})},
            {"observation": dict(self.candidate["observation"], record_url=self.a["url"])},
        ):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.check(dict(self.candidate, **changes))
        with self.assertRaisesRegex(ValueError, "synthetic fixture"):
            claims.fixture_response("investigate", dict(self.bundle, synthetic=False), self.job)

    def test_explicit_abstention_with_available_data(self):
        # A model may find the retrieved fields insufficient for its research question.
        self.candidate = claims.insufficient_response("investigate", self.bundle)
        self.assertEqual(self.candidate, self.check(self.candidate))
        with self.assertRaises(ValueError):
            self.check(dict(self.candidate, evidence=[self.cite(self.a, "context")]))
        response = claims.insufficient_response("critique", self.bundle)
        self.assertEqual(response, self.critique(response))
        with self.assertRaises(ValueError):
            self.critique(dict(response, evidence=[self.cite(self.a, "context")]))


if __name__ == "__main__":
    unittest.main()
