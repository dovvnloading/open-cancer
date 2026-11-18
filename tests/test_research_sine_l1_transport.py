"""Offline, cross-platform stdlib integrity tests for frozen ranked research job 5."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import socket
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "sine-l1-transport"
SCRIPT = RESEARCH / "analysis.py"
SPEC = importlib.util.spec_from_file_location("research_sine_l1_transport_analysis", SCRIPT)
ANALYSIS = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(ANALYSIS)

# Numerical snapshots: absolute+relative 3e-7. Optimizer stopping traces are
# deliberately tested against protocol bounds, NEVER snapshot exact floats.
TOLERANCE = 3e-7
EXPECTED = {
    "baseline": 0.45381839464321067,
    "sine": 0.4433046516276496,
    "gain": 0.010513743015561072,
    "jackknife_center": 0.02513552622558116,
    "jackknife_se": 0.010670012463819434,
    "lower": 0.0013612569039949396,
    "upper": 0.04890979554716738,
    "candidate_gain_min": 0.009535222478026628,
    "candidate_gain_max": 0.011529692807449643,
    "candidate_lower_min": 0.00005099412229885719,
    "candidate_upper_max": 0.05256644352134192,
}


def assert_close(test: unittest.TestCase, actual: float, expected: float):
    test.assertTrue(math.isclose(actual, expected, rel_tol=TOLERANCE, abs_tol=TOLERANCE),
                    (actual, expected))


class SourceAndMath(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = ANALYSIS.load_panel()

    def test_pinned_minimal_derived_panel_and_count_gate(self):
        rows = self.rows
        self.assertEqual((len(rows), len(set(x["order"] for x in rows))), (55, 11))
        self.assertEqual(sum(r["n"] for r in rows), 3533)
        self.assertEqual(sum(r["neo_k"] for r in rows), 691)
        self.assertEqual(sum(r["mal_k"] for r in rows), 405)
        self.assertEqual(sum(r["x"][0] == 0 for r in rows), 9)
        self.assertEqual(sum(not r["ratio_ok"] for r in rows), 31)
        self.assertEqual(math.prod(len(r["neo_candidates"]) for r in rows), 45)
        self.assertEqual([r["species"] for r in rows if len(r["neo_candidates"]) > 1],
                         ["Cercopithecus_mitis", "Giraffa_camelopardalis", "Panthera_pardus"])
        self.assertEqual([r["species"] for r in rows], sorted(r["species"] for r in rows))
        self.assertFalse(any(r["neo_k"] == 0 or r["mal_k"] == 0 for r in rows),
                         "upstream panel excludes zero-recorded-case species")
        self.assertEqual(len(list(ANALYSIS.all_assignments(rows))), 45)
        self.assertEqual(len({tuple(sorted(a.items())) for a in ANALYSIS.all_assignments(rows)}), 45)
        self.assertTrue(all(r["neo_k"] in r["neo_candidates"] for r in rows))
        self.assertTrue(all(r["mal_k"] in r["mal_candidates"] for r in rows))

    def test_source_licence_audit_and_secondary_abstention(self):
        audit = json.loads((RESEARCH / "data" / "source_audit.json").read_text(encoding="utf-8"))
        self.assertIn("CC BY 4.0", audit["data_origin"]["article_license"])
        self.assertIn("null", audit["data_origin"]["repository_license"])
        self.assertEqual(audit["data_origin"]["pinned_commit"],
                         "792af204b4113238276b3ac43e511b62872dd368")
        self.assertEqual(audit["audit"]["total_neoplasia_nominal"], 691)
        self.assertEqual(audit["audit"]["total_malignancy_nominal"], 405)
        self.assertEqual(audit["audit"]["secondary_ratio_precision_mismatch_count"], 31)
        self.assertFalse(audit["audit"]["secondary_gate_pass"])
        self.assertEqual(audit["audit"]["primary_compatible_assignment_count"], 45)

    def test_endpoints_and_numerical_primitives(self):
        self.assertEqual(ANALYSIS.cross_entropy(0, 10, 0.1)[0], -math.log1p(-0.1))
        self.assertEqual(ANALYSIS.cross_entropy(10, 10, 0.1)[0], -math.log(0.1))
        self.assertTrue(ANALYSIS.cross_entropy(0, 10, 0.0)[1])
        self.assertEqual(ANALYSIS.softplus(1000), 1000)
        self.assertEqual(ANALYSIS.softplus(-1000), 0)
        self.assertEqual(ANALYSIS.sigmoid(1000), 1.0)
        self.assertEqual(ANALYSIS.sigmoid(-1000), 0.0)
        sol = ANALYSIS.cholesky_solve([[4.0, 2.0], [2.0, 3.0]], [6.0, 5.0])
        assert_close(self, 4 * sol[0] + 2 * sol[1], 6)
        assert_close(self, 2 * sol[0] + 3 * sol[1], 5)
        with self.assertRaises(ANALYSIS.AnalysisError):
            ANALYSIS.cholesky_solve([[1, 2], [2, 1]], [1, 2])

    def test_training_only_transforms_and_heldout_label_isolation(self):
        rows = self.rows
        held = [r for r in rows if r["order"] == "Cingulata"]
        train = [r for r in rows if r["order"] != "Cingulata"]
        a, b, scale1 = ANALYSIS.transform(train, held, True, True)
        self.assertEqual(len(a[0]), 6)
        self.assertEqual(len(b[0]), 6)
        extreme = dict(held[0], x=(0.0, 100000.0, 0.0, 1e12, 0.0), neo_k=0)
        _, _, scale2 = ANALYSIS.transform(train, [extreme], True, True)
        self.assertEqual(scale1, scale2, "held-out features must not change training scaling")
        original = ANALYSIS.score_orders(rows, detailed=True)
        replaced_rows = [dict(r, neo_k=0) if r["order"] == "Cingulata" else r for r in rows]
        changed = ANALYSIS.score_orders(replaced_rows, detailed=True)
        orig_pred = [p for p in original["predictions"] if p["held_out_order"] == "Cingulata"][0]
        other_pred = [p for p in changed["predictions"] if p["held_out_order"] == "Cingulata"][0]
        assert_close(self, orig_pred["baseline_p"], other_pred["baseline_p"])
        assert_close(self, orig_pred["sine_p"], other_pred["sine_p"])
        self.assertNotEqual(orig_pred["baseline_loss"], other_pred["baseline_loss"])

    def test_newline_normalized_hashes_and_posix_keys(self):
        with tempfile.TemporaryDirectory() as td:
            temp = Path(td)
            a, b = temp / "lf.txt", temp / "crlf.txt"
            a.write_bytes(b"A\nB\n")
            b.write_bytes(b"A\r\nB\r\n")
            self.assertEqual(ANALYSIS.normalized_hash(a), ANALYSIS.normalized_hash(b))
            self.assertEqual(ANALYSIS.normalized_hash(a), hashlib.sha256(b"A\nB\n").hexdigest())
        self.assertEqual((RESEARCH / "data" / "species_aggregates.csv").relative_to(RESEARCH).as_posix(),
                         "data/species_aggregates.csv")
        self.assertEqual(ANALYSIS.normalized_hash(ANALYSIS.PANEL),
                         "79736cfcf782900caceec5d7eeb2f7b5db89a1497b158f5f7088faa8ebf0a0e2")


class EntireFrozenRefit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # This explicit network denial makes test failure immediate if a later
        # revision accidentally introduces a download into the offline tests.
        with patch.object(socket, "socket", side_effect=AssertionError("offline-only test")):
            cls.outcome = ANALYSIS.run_all(exhaustive=True)

    def test_nominal_gain_jackknife_independently_recovered(self):
        o = self.outcome["nominal"]
        assert_close(self, o["baseline_loss"], EXPECTED["baseline"])
        assert_close(self, o["sine_loss"], EXPECTED["sine"])
        assert_close(self, o["gain"], EXPECTED["gain"])
        self.assertEqual(o["orders_improving"], 8)
        self.assertEqual(len(o["per_order"]), 11)
        self.assertEqual(len(o["predictions"]), 55)
        self.assertEqual(len(o["diagnostics"]), 22)
        self.assertEqual(o["probability_clips"], 0)
        j = o["jackknife"]
        self.assertEqual(len(j["deletions"]), 11)
        for key, expected in (("center", "jackknife_center"), ("se", "jackknife_se"),
                              ("lower", "lower"), ("upper", "upper")):
            assert_close(self, j[key], EXPECTED[expected])
        q = [11 * o["gain"] - 10 * d["gain"] for d in j["deletions"]]
        center = sum(q) / 11
        se = math.sqrt(sum((v - center)**2 for v in q) / (11*10))
        assert_close(self, j["center"], center)
        assert_close(self, j["se"], se)
        assert_close(self, j["lower"], center - ANALYSIS.T95_10*se)
        assert_close(self, j["upper"], center + ANALYSIS.T95_10*se)
        self.assertTrue(all(d["gain"] > 0 for d in j["deletions"]))

    def test_independent_prediction_loss_and_order_balancing(self):
        o = self.outcome["nominal"]
        groups = {}
        for r in o["predictions"]:
            p = r["baseline_p"]
            q = r["sine_p"]
            k, n = r["k"], r["n"]
            self.assertTrue(0 < p < 1 and 0 < q < 1)
            calc = lambda value: -k*math.log(value)/n - (n-k)*math.log1p(-value)/n
            assert_close(self, r["baseline_loss"], calc(p))
            assert_close(self, r["sine_loss"], calc(q))
            groups.setdefault(r["held_out_order"], []).append(r)
        self.assertEqual(set(groups), {o["order"] for o in o["per_order"]})
        for order in o["per_order"]:
            members = groups[order["order"]]
            assert_close(self, order["baseline_loss"], sum(x["baseline_loss"] for x in members)/len(members))
            assert_close(self, order["sine_loss"], sum(x["sine_loss"] for x in members)/len(members))
        assert_close(self, o["baseline_loss"], sum(x["baseline_loss"] for x in o["per_order"])/11)
        assert_close(self, o["sine_loss"], sum(x["sine_loss"] for x in o["per_order"])/11)

    def test_gradient_and_optimisation_within_declared_bounds_not_exact_floats(self):
        for f in self.outcome["nominal"]["diagnostics"]:
            self.assertTrue(f["converged"])
            self.assertTrue(0 <= f["iterations"] <= 100)
            self.assertTrue(0 <= f["line_search_reductions"] <= 44*f["iterations"])
            self.assertGreater(f["effective_n"], 0)
            self.assertGreaterEqual(f["objective"], 0)
            self.assertLessEqual(f["max_gradient"], f["tolerance"]*(1+1e-8))
            self.assertLessEqual(f["gradient_tolerance_ratio"], 1+1e-8)
            self.assertTrue(all(math.isfinite(b) for b in f["beta"]))
            self.assertTrue(all(s > 0 for s in f["train_sd"]))
            self.assertEqual(f["train_species"] + f["test_species"], 55)

    def test_exhaustive_grid_sensitivities_and_frozen_verdict(self):
        data = self.outcome
        self.assertEqual(data["verdict"], "consistent")
        self.assertTrue(data["nonclinical_unreviewed"])
        self.assertEqual(len(data["candidate_assignments"]), 45)
        cases = data["candidate_assignments"]
        self.assertEqual(set(x["verdict"] for x in cases), {"consistent"})
        assert_close(self, min(v["gain"] for v in cases), EXPECTED["candidate_gain_min"])
        assert_close(self, max(v["gain"] for v in cases), EXPECTED["candidate_gain_max"])
        assert_close(self, min(v["interval"][0] for v in cases), EXPECTED["candidate_lower_min"])
        assert_close(self, max(v["interval"][1] for v in cases), EXPECTED["candidate_upper_max"])
        for case in cases:
            self.assertGreater(case["gain"], 0)
            self.assertGreater(case["interval"][0], 0)
            self.assertEqual(len(case["deletion_gains"]), 11)
            self.assertTrue(all(v > 0 for v in case["deletion_gains"].values()))
            self.assertGreater(case["orders_improving"], 5)
            self.assertGreater(case["assembly_quality_gain"], 0)
            self.assertTrue(all(v > 0 for v in case["sensitivity_gains"].values()))
        s = data["sensitivities"]
        self.assertEqual((s["assembly_quality"]["species"], s["assembly_quality"]["orders"]), (38, 9))
        self.assertGreater(s["assembly_quality"]["score"]["gain"], 0)
        self.assertTrue(all(s[k]["gain"] > 0 for k in
                            ("equal_species_training", "ridge_0_1", "ridge_10", "training_winsor_n50")))
        self.assertEqual(data["audit"]["secondary_gate"],
                         "failed: 31 source Ratio precision mismatches; not scored")
        self.assertTrue(all("\\" not in k for k in data["input_hashes"]))
        self.assertEqual(data["input_hashes"]["data/species_aggregates.csv"],
                         ANALYSIS.normalized_hash(ANALYSIS.PANEL))
        self.assertEqual(data["input_hashes"]["analysis.py"],
                         ANALYSIS.normalized_hash(SCRIPT))

    def test_negative_verdict_boundaries_are_exhaustive(self):
        fake = {"gain": -0.01, "jackknife": {"lower": -0.03, "upper": -0.001},
                "orders_improving": 0}
        self.assertEqual(ANALYSIS.verdict(fake, self.outcome["sensitivities"]), "falsified")
        fake["jackknife"]["upper"] = 0.0
        self.assertEqual(ANALYSIS.verdict(fake, self.outcome["sensitivities"]), "inconclusive")
        fake["gain"] = 0.01
        fake["jackknife"] = {"lower": -0.001, "upper": 0.02}
        self.assertEqual(ANALYSIS.verdict(fake, self.outcome["sensitivities"]), "inconclusive")


if __name__ == "__main__":
    unittest.main()
