"""Synthetic packet and response pairs for measuring the pilot scorer itself."""
import copy
import unittest

from scripts import pilot_evaluation as score


class PilotEvaluationTests(unittest.TestCase):
    def setUp(self):
        self.job = {"id": "j1", "source_ids": ["civic", "gdc", "cbioportal"]}
        self.ref = {"source_id": "civic", "record_id": "1", "url": "https://fixture.invalid/1"}
        self.packet = {
            "records": [{**self.ref, "disease": "Lung Adenocarcinoma",
                         "molecular_profile": "EGFR L858R", "significance": ""}],
            "sources": {"civic": {"status": "fixture", "retrieval_status": "complete", "matching_records": 12},
                        "gdc": {"status": "not-retrieved"},
                        "cbioportal": {"status": "not-retrieved"}}}
        self.obs = {"id": "o1", "text": "The supplied packet contains 1 CIViC evidence item.",
                    "citations": [self.ref], "checks": [self.count(1)]}
        self.candidate = {"observations": [self.obs], "contradictory_evidence": [],
                          "assumptions": ["Synthetic examples."], "limitations": ["No cases."]}
        self.critique = {"verdict": "insufficient-evidence", "issues": [
            {"category": "missing-denominator", "observation_id": None,
             "detail": "No cohort denominator.", "citations": []}]}

    @staticmethod
    def count(n, source="civic", unit="evidence-items"):
        return {"kind": "packet-count", "source_id": source, "value": n, "unit": unit}

    @staticmethod
    def field(field, value, rid="1"):
        return {"kind": "record-field", "record_id": rid, "field": field, "value": value}

    def body(self, complete=True, attempts=()):
        critique = self.critique if complete else None
        return {"job": self.job, "candidate": self.candidate, "critique": critique,
                "status": "completed" if complete else "failed", "attempts": list(attempts),
                "evaluation": score.evaluate_job(self.job, self.packet, self.candidate, critique)}

    def test_citations_match_exact_record_id_url_and_source_in_each_role(self):
        self.obs["citations"] += [{**self.ref, "source_id": "gdc"},
                                  {**self.ref, "url": "https://fixture.invalid/"},
                                  {**self.ref, "record_id": "2"}]
        self.critique["issues"][0]["citations"] = [self.ref]
        result = self.body()["evaluation"]
        self.assertEqual([True, False, False, False, True],
                         [c["valid"] for c in result["citations"]])
        totals = score.aggregate([self.body(), self.body(complete=False)])
        self.assertEqual({"valid": 3, "unsupported": 6}, totals["citations"])
        self.assertEqual(2, totals["observations"]["flagged_unsupported"])
        self.assertEqual(1, totals["jobs"]["with_critique"])
        self.assertEqual(1, totals["jobs"]["failed"])
        self.assertEqual(1, totals["source_coverage"]["civic"]["unique_records_supplied"])
        self.assertEqual(1, totals["source_coverage"]["civic"]["unique_records_cited_by_investigator"])
        self.assertEqual(1, totals["source_coverage"]["civic"]["fraction_supplied_cited"])
        self.assertEqual(2, totals["source_coverage"]["gdc"]["jobs_without_retrieval"])
        self.assertIsNone(totals["source_coverage"]["cbioportal"]["fraction_supplied_cited"])

    def test_recomputed_counts_fields_units_and_unavailable_denominators(self):
        self.obs["checks"] = [self.count(1), self.count(12), self.count(1, unit="patients"),
                              self.count(0, source="gdc"), self.field("disease", "Lung Adenocarcinoma"),
                              self.field("molecular_profile", "KRAS"), self.field("significance", ""),
                              self.field("disease", "Lung Adenocarcinoma", rid="99"),
                              self.field("cohort_size", "100")]
        body = self.body()
        self.assertEqual(["passed", "failed", "failed", "unassessed", "passed", "failed",
                          "unassessed", "unassessed", "unassessed"],
                         [c["status"] for c in body["evaluation"]["checks"]])
        self.assertEqual({"passed": 2, "failed": 3, "unassessed": 4},
                         score.aggregate([body])["deterministic_checks"])

    def test_unstructured_count_conflicting_with_structured_check_and_cohort_review(self):
        self.obs["text"] = "The supplied packet contains 56 evidence items, from 31 patients."
        body = self.body()
        self.assertEqual([("passed", "model"), ("failed", "prose")],
                         [(c["status"], c["origin"]) for c in body["evaluation"]["checks"]])
        summary = score.aggregate([body])
        self.assertEqual(1, summary["observations"]["with_cohort_quantity_to_review"])
        self.assertEqual(1, summary["observations"]["flagged_unsupported"])
        self.obs.update(text="Examining EGFR L858R and exon 19 deletion labels.", checks=[])
        summary = score.aggregate([self.body()])
        self.assertEqual(1, summary["observations"]["without_checkable_assertions"])
        self.assertEqual(0, summary["observations"]["with_cohort_quantity_to_review"])

    def test_direct_claims_of_examination_require_source_records(self):
        self.obs["text"] = "We examined the public GDC occurrence rows and we inspected cBioPortal records."
        result = self.body()["evaluation"]
        self.assertEqual(["gdc", "cbioportal"],
                         [check["assertion"]["source_id"] for check in result["checks"]
                          if check["status"] == "failed"])
        self.obs["text"] = "We did not inspect GDC records. Historical cBioPortal records are out of scope."
        self.assertFalse(self.body()["evaluation"]["observations"][0]["unsupported_reasons"])

    def test_retrieved_zero_match_separate_from_missing_source_and_fixture(self):
        self.packet["records"] = []
        self.packet["sources"]["civic"].update(status="retrieved", matching_records=0)
        self.candidate["observations"] = []
        result = score.aggregate([self.body()], mode="real")
        self.assertEqual(1, result["source_coverage"]["civic"]["jobs_without_matches"])
        self.assertEqual(0, result["source_coverage"]["civic"]["jobs_without_retrieval"])
        self.assertEqual(1, result["source_coverage"]["civic"]["jobs_live_retrieved"])
        self.assertEqual(0, result["source_coverage"]["gdc"]["jobs_without_matches"])
        self.assertIsNone(result["execution"]["estimated_cost_usd"])

    def test_truncated_search_with_no_observed_match_preserves_unknown_scope(self):
        self.packet["records"] = []
        self.packet["sources"]["civic"].update(status="retrieved", retrieval_status="truncated",
                                              matching_records=0)
        self.candidate["observations"] = []
        body = self.body()
        coverage = score.aggregate([body], mode="real")["source_coverage"]["civic"]
        self.assertEqual(0, coverage["jobs_without_matches"])
        self.assertEqual(1, coverage["jobs_with_truncated_retrieval"])
        self.assertEqual(1, coverage["jobs_without_matches_before_truncation"])
        self.assertEqual(0, coverage["jobs_without_retrieval"])
        self.assertIsNone(coverage["fraction_supplied_cited"])

    def test_critic_flags_and_counterevidence_counted_separately(self):
        self.candidate["contradictory_evidence"] = [{**self.obs, "id": "c1"}]
        self.critique["issues"] += [{"category": "contradiction", "observation_id": "c1",
                                     "detail": "The candidate's interpretations conflict.",
                                     "citations": [self.ref]},
                                    {"category": "unsupported-claim", "observation_id": "o1",
                                     "detail": "The interpretation needs independent review.",
                                     "citations": []}]
        summary = score.aggregate([self.body()])
        self.assertEqual(1, summary["observations"]["contradictory_evidence_entries"])
        self.assertEqual(1, summary["critic_issues_by_category"]["contradiction"])
        self.assertEqual(1, summary["critic_issues_by_category"]["unsupported-claim"])
        self.assertEqual(1, summary["critic_issues_by_category"]["missing-denominator"])
        self.assertEqual(0, summary["observations"]["flagged_unsupported"])

    @staticmethod
    def attempt(stage, prompt=None, completion=None, error=None):
        return {"stage": stage, "usage": {"prompt_tokens": prompt, "completion_tokens": completion},
                "elapsed_seconds": 1.5, "error": error, "synthetic": False, "endpoint_attempted": True}

    def test_cost_counts_failed_billed_attempts_and_preserves_partial_usage(self):
        prices = {"investigate": {"input_per_million_usd": 2.50, "output_per_million_usd": 9},
                  "critique": {"input_per_million_usd": 1, "output_per_million_usd": 5}}
        attempts = [self.attempt("investigate", 10, 7, "bad-schema"),
                    self.attempt("investigate", 20, 5),
                    self.attempt("critique", 40, 15)]
        result = score.aggregate([self.body(attempts=attempts)], prices, mode="real")["execution"]
        self.assertEqual(97, result["usage"]["total_tokens"])
        self.assertEqual(1, result["failed_attempts"])
        self.assertEqual(4.5, result["elapsed_attempt_seconds"])
        self.assertAlmostEqual(.000298, result["estimated_cost_usd"])
        attempts[-1]["usage"]["completion_tokens"] = None
        partial = score.aggregate([self.body(attempts=attempts)], prices, mode="real")["execution"]
        self.assertIsNone(partial["estimated_cost_usd"])
        self.assertAlmostEqual(.000223, partial["known_estimated_cost_subtotal_usd"])
        self.assertEqual(12, partial["known_token_subtotals"]["completion_tokens"])
        self.assertIsNone(partial["usage"]["completion_tokens"])
        self.assertIsNone(partial["usage"]["total_tokens"])
        self.assertEqual(70, partial["usage"]["prompt_tokens"])
        self.assertEqual(2, partial["attempts_with_token_counts"]["total_tokens"])
        self.assertEqual(5, partial["priced_input_or_output_components"])

    def test_interruption_preserves_unknown_billing_and_elapsed_time(self):
        prices = {"investigate": {"input_per_million_usd": 1, "output_per_million_usd": 2}}
        complete = self.attempt("investigate", 100, 50)
        interrupted = {"stage": "critique", "elapsed_seconds": None, "endpoint_attempted": None,
                       "usage": None, "error": "interrupted-attempt-outcome-unknown", "synthetic": False}
        body = self.body(complete=False, attempts=[complete, interrupted])
        execution = score.aggregate([body], prices, mode="real")["execution"]
        self.assertEqual(2, execution["attempts"])
        self.assertEqual(1, execution["endpoint_attempts"])
        self.assertEqual(1, execution["attempts_of_unknown_dispatch"])
        self.assertEqual(1, execution["attempts_with_elapsed_seconds"])
        self.assertIsNone(execution["elapsed_attempt_seconds"])
        self.assertEqual(1.5, execution["known_elapsed_attempt_seconds"])
        self.assertIsNone(execution["estimated_cost_usd"])
        self.assertAlmostEqual(0.0002, execution["known_estimated_cost_subtotal_usd"])
        self.assertIsNone(execution["usage"]["total_tokens"])
        self.assertEqual(150, execution["known_token_subtotals"]["total_tokens"])
        self.assertIsNone(execution["by_stage"]["critique"]["known_elapsed_attempt_seconds"])

    def test_metering_preserves_zero_and_excludes_non_endpoint_or_synthetic_costs(self):
        prices = {"investigate": {"input_per_million_usd": 0, "output_per_million_usd": 0}}
        a = self.attempt("investigate", 0, 0)
        skipped = {**self.attempt("critique", error="invalid-input"), "endpoint_attempted": False}
        result = score.aggregate([self.body(attempts=[a, skipped])], prices, mode="real")["execution"]
        self.assertEqual(1, result["endpoint_attempts"])
        self.assertEqual(0, result["usage"]["total_tokens"])
        self.assertEqual(0.0, result["estimated_cost_usd"])
        with self.assertRaisesRegex(ValueError, "synthetic"):
            score.aggregate([self.body(attempts=[a])], prices, mode="fixture")
        a.update(synthetic=True, endpoint_attempted=False)
        fixture = score.aggregate([self.body(attempts=[a])], prices)["execution"]
        self.assertIsNone(fixture["estimated_cost_usd"])
        self.assertIsNone(fixture["usage"]["total_tokens"])
        self.assertEqual(0, fixture["endpoint_attempts"])
        for value in [-1, "1", float("inf"), float("nan"), True]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                score.validate_prices({"investigate": {"input_per_million_usd": value}})

    def test_output_consistency_excludes_failed_and_missing_pairs(self):
        first = self.body()
        second = copy.deepcopy(first)
        failed = self.body(complete=False)
        failed["job"] = {**self.job, "id": "j2"}
        other = copy.deepcopy(failed)
        other["job"]["id"] = "j3"
        other["status"] = "completed"
        comparison = score.compare_outputs([first, failed], [second, other])
        self.assertEqual(2, comparison["missing_or_failed_pairs"])
        self.assertEqual(1, comparison["paired_completed"])
        self.assertEqual(1.0, comparison["mean_cited_record_jaccard"])
        self.assertEqual(1.0, comparison["exact_fraction"])
        second["candidate"]["limitations"].append("A revised limitation.")
        self.assertEqual(0.0, score.compare_outputs([first], [second])["exact_fraction"])
        no_pairs = score.compare_outputs([failed], [])
        self.assertIsNone(no_pairs["exact_fraction"])
        self.assertIsNone(no_pairs["mean_cited_record_jaccard"])
        self.assertEqual("missing", no_pairs["pairs"][0]["second_status"])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            score.compare_outputs([failed, failed], [])

    def test_comparison_detects_served_build_drift_with_stable_request_alias(self):
        stage = self.attempt("investigate")
        served = {**stage, "returned_model": "served-a", "system_fingerprint": None}
        # A failed attempt's provider identity is unrelated to the accepted response.
        failed = {**served, "returned_model": "failed-build", "error": "invalid-model-schema"}
        critic = {**self.attempt("critique"), "returned_model": "critic",
                  "system_fingerprint": "critic-build-a"}
        first = self.body(attempts=[failed, served, critic])
        second = copy.deepcopy(first)
        second["attempts"][1].update(returned_model="served-b", system_fingerprint="investigator-b")
        second["attempts"][2]["system_fingerprint"] = "critic-build-b"
        # The compare CLI has already checked that both runs requested the same model configuration.
        row = score.compare_outputs([first], [second])["pairs"][0]
        identifiers = row["provider_identifiers"]
        self.assertTrue(row["candidate_exact"])
        self.assertEqual({"first": "served-a", "second": "served-b", "differ": True},
                         identifiers["investigate"]["returned_model"])
        self.assertIsNone(identifiers["investigate"]["system_fingerprint"]["differ"])
        self.assertFalse(identifiers["critique"]["returned_model"]["differ"])
        self.assertTrue(identifiers["critique"]["system_fingerprint"]["differ"])
        for body in (first, second):
            for attempt in body["attempts"]:
                attempt.update(synthetic=True, endpoint_attempted=False)
        fixture = score.compare_outputs([first], [second])["pairs"][0]["provider_identifiers"]
        for stage in fixture.values():
            for identity in stage.values():
                self.assertEqual({"first": None, "second": None, "differ": None}, identity)

    def test_source_failure_has_no_fake_zero_counts_or_role_outputs(self):
        evaluation = score.evaluate_job(self.job, None)
        self.assertIsNone(evaluation["source_coverage"]["civic"]["matching_records"])
        self.assertEqual({"investigate": None, "critique": None}, evaluation["output_sha256"])
        totals = score.aggregate([{"job": self.job, "evaluation": evaluation,
                                   "attempts": [], "status": "failed"}])
        self.assertEqual(1, totals["jobs"]["failed"])
        self.assertEqual(0, totals["jobs"]["with_candidate"])
        self.assertIsNone(totals["execution"]["usage"]["total_tokens"])
