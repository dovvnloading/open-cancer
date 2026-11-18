"""Nonclinical aggregate CH study: independent calculations and integrity checks.

Stdlib only; no network, source downloads, external executables, or patient
data. This test checks reproducibility and mathematical implementations, not
source truth, biological novelty, or individual clinical prediction.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import math
from pathlib import Path
import random
import statistics
import sys
import unittest


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "research" / "clonal-hematopoiesis-progression"


def _load_modules():
    # Script modules have short import names. Isolate these names so discovery
    # neither consumes nor changes modules imported by unrelated repo tests.
    names = ("model", "simulate", "analysis", "verify")
    old = {name: sys.modules.get(name) for name in names}
    loaded = {}
    try:
        for name in names:
            spec = importlib.util.spec_from_file_location(name, ROOT / (name + ".py"))
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            spec.loader.exec_module(module)
            loaded[name] = module
    finally:
        for name, previous in old.items():
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous
    return (loaded[name] for name in names)


MODEL, SIMULATE, ANALYSIS, VERIFY = _load_modules()


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def simpson(function, horizon, intervals=600):
    if horizon == 0:
        return 0.0
    step = horizon / intervals
    total = function(0.0) + function(horizon)
    total += sum((4 if i % 2 else 2) * function(i * step)
                 for i in range(1, intervals))
    return total * step / 3


def numerical_survival(birth, death, hit, horizon, steps=2048):
    """Independent RK4 for the backward equation, not the closed-form roots."""
    def derivative(u):
        return birth * u * u - (birth + death + hit) * u + death
    step, u = horizon / steps, 1.0
    for _ in range(steps):
        k1 = derivative(u)
        k2 = derivative(u + step * k1 / 2)
        k3 = derivative(u + step * k2 / 2)
        k4 = derivative(u + step * k3)
        u += step * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return u


class ClonalHematopoiesisProgressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = ANALYSIS.run(ROOT)
        cls.fitness = read_csv(ROOT / "data" / "fitness.csv")
        cls.outcomes = read_csv(ROOT / "data" / "outcomes.csv")

    def test_all_reported_results_recompute_and_all_files_have_normalized_hashes(self):
        checked = VERIFY.verify_manifest(ROOT)
        self.assertIn("tests/test_research_clonal_hematopoiesis_progression.py", checked)
        self.assertTrue(all(p.startswith("research/clonal-hematopoiesis-progression/") or
                            p == "tests/test_research_clonal_hematopoiesis_progression.py"
                            for p in checked))
        self.assertEqual(VERIFY.digest("a\nb\nc\n"), VERIFY.digest(b"a\r\nb\rc\n"))
        text = (ROOT / "README.md").read_text(encoding="utf-8")
        actual = ANALYSIS.BEGIN + text.split(ANALYSIS.BEGIN, 1)[1].split(ANALYSIS.END, 1)[0] + ANALYSIS.END
        self.assertEqual(actual, ANALYSIS.render_report(self.report))
        # README contains the complete proof text, not only a link to it.
        proofs = (ROOT / "THEOREMS.md").read_text(encoding="utf-8").strip()
        self.assertIn(proofs, text)
        commit = VERIFY.normalized_bytes((ROOT / "PRESPEC_COMMIT.txt").read_bytes())
        git_object = b"commit " + str(len(commit)).encode("ascii") + b"\0" + commit
        self.assertEqual(hashlib.sha1(git_object).hexdigest(),
                         "7409aa2e15be9215d337ae8e27cef14b628e6c45")

    def test_provenance_units_and_complete_source_gene_inclusion(self):
        VERIFY.verify_provenance(ROOT)
        self.assertEqual(VERIFY.verify_derived_outcomes(ROOT), 22)
        source_genes = {
            "ASXL1", "BCOR", "CBL", "DNMT3A", "IDH1", "IDH2", "JAK2", "KDM6A",
            "KMT2C", "KMT2D", "KRAS", "NF1", "NRAS", "PHF6", "PTPN11", "RAD21",
            "RUNX1", "SF3B1", "SRSF2", "TET2", "TP53", "U2AF1"}
        self.assertEqual({r["gene"] for r in self.outcomes if r["raw_log_estimate"]}, source_genes)
        self.assertEqual(len(self.fitness), 23)
        self.assertEqual(len({r["gene"] for r in self.fitness}), 17)
        for row in self.fitness:
            lo, center, hi = (float(row[k]) for k in ("lower", "estimate", "upper"))
            self.assertLessEqual(lo, center)
            self.assertLessEqual(center, hi)
            self.assertEqual(row["units"], "year^-1")
            self.assertEqual(row["interval_type"], "95% highest posterior density interval")
            for key in ("estimate_cell", "lower_cell", "upper_cell", "sheet", "interval_audit_url"):
                self.assertTrue(row[key])
        included = {r["gene"] for r in self.report["primary"]["rows"]}
        excluded = {r["gene"] for r in self.report["excluded"]}
        self.assertFalse(included & excluded)
        self.assertEqual(included | excluded, source_genes)
        self.assertEqual(included, {"ASXL1", "IDH1", "KRAS", "PTPN11", "SF3B1", "SRSF2", "U2AF1"})

    def test_independent_loo_and_every_sensitivity_from_csvs(self):
        # Select by source measurement scope independently of ANALYSIS helper.
        groups = {}
        for row in self.fitness:
            groups.setdefault(row["gene"], []).append(row)
        evidence = []
        for outcome in self.outcomes:
            if outcome["study_id"] != "fabre2022_abelson_refit":
                continue
            fits = groups.get(outcome["gene"], [])
            if len(fits) == 1 and fits[0]["estimator_scope"] == "gene_class_baseline_effect":
                evidence.append((outcome["gene"], fits[0], outcome))
        scenarios = read_csv(ROOT / "analysis_scenarios.csv")
        for scenario, computed in zip(scenarios, self.report["scenarios"]):
            horizon = float(scenario["horizon_years"])
            def log_opportunity(s):
                fraction = 2 * float(scenario["vaf"]) if scenario["vaf"] else 0.0
                def path(t):
                    growth = math.exp(s * t)
                    return growth / (1 - fraction + fraction * growth)
                return math.log(simpson(path, horizon, 1000))
            independent = []
            for gene, fit, outcome in evidence:
                beta = float(outcome["raw_log_estimate"])
                se = float(outcome["log_standard_error"])
                training = [(f, o) for g, f, o in evidence if g != gene]
                intercept = statistics.fmean(float(o["raw_log_estimate"]) -
                                             log_opportunity(float(f["estimate"])) for f, o in training)
                lower_intercept = statistics.fmean(float(o["raw_log_estimate"]) - 1.96 * float(o["log_standard_error"]) -
                                                   log_opportunity(float(f["upper"])) for f, o in training)
                upper_intercept = statistics.fmean(float(o["raw_log_estimate"]) + 1.96 * float(o["log_standard_error"]) -
                                                   log_opportunity(float(f["lower"])) for f, o in training)
                log_prediction = intercept + log_opportunity(float(fit["estimate"]))
                low = lower_intercept + log_opportunity(float(fit["lower"]))
                high = upper_intercept + log_opportunity(float(fit["upper"]))
                row = next(r for r in computed["rows"] if r["gene"] == gene)
                self.assertAlmostEqual(row["prediction"], math.exp(log_prediction), delta=1e-7)
                self.assertAlmostEqual(row["prediction_lower"], math.exp(low), delta=1e-7)
                self.assertAlmostEqual(row["prediction_upper"], math.exp(high), delta=1e-6)
                self.assertEqual(row["envelopes_overlap"], max(beta - 1.96 * se, low) <= min(beta + 1.96 * se, high))
                self.assertAlmostEqual(row["observed_over_predicted"], math.exp(beta - log_prediction), delta=1e-7)
                null = statistics.fmean(float(o["raw_log_estimate"]) for _, o in training)
                independent.append((beta - log_prediction, beta - null))
            rmse = math.sqrt(statistics.fmean(e[0]**2 for e in independent))
            null_rmse = math.sqrt(statistics.fmean(e[1]**2 for e in independent))
            self.assertAlmostEqual(computed["log_rmse"], rmse, places=9)
            self.assertAlmostEqual(computed["null_log_rmse"], null_rmse, places=12)
            self.assertAlmostEqual(computed["median_absolute_log_error"],
                                   statistics.median(abs(e[0]) for e in independent), places=9)
        self.assertEqual(self.report["primary"]["failed_genes"], ["U2AF1"])
        self.assertEqual(next(r for r in self.report["scenarios"] if r["scenario_id"] == "long_20y")["failed_genes"], [])
        self.assertEqual(self.report["mechanistic_verdict"], "untestable")

    def test_heldout_outcome_never_changes_its_own_prediction(self):
        data = ANALYSIS.load_inputs(ROOT)
        pairs, _ = ANALYSIS.select_proxy_pairs(data)
        scenario = data["scenarios"][0]
        before = ANALYSIS.loo(pairs, scenario)
        held = pairs[-1]
        held["outcome"] = dict(held["outcome"])
        held["outcome"]["estimate"] = str(float(held["outcome"]["estimate"]) * 100)
        after = ANALYSIS.loo(pairs, scenario)
        old = next(r for r in before["rows"] if r["gene"] == held["gene"])
        new = next(r for r in after["rows"] if r["gene"] == held["gene"])
        for key in ("prediction", "prediction_lower", "prediction_upper"):
            self.assertEqual(old[key], new[key])

    def test_post_hoc_intercept_feasibility_and_class_sensitivity(self):
        data = ANALYSIS.load_inputs(ROOT)
        pairs, _ = ANALYSIS.select_proxy_pairs(data)
        for scenario, row in zip(data["scenarios"], self.report["joint_box_audit"]):
            intervals = []
            for pair in pairs:
                o, f = pair["outcome"], pair["fitness"]
                lo = math.log(float(o["lower"])) - math.log(ANALYSIS.opportunity(float(f["upper"]), scenario))
                hi = math.log(float(o["upper"])) - math.log(ANALYSIS.opportunity(float(f["lower"]), scenario))
                intervals.append((lo, hi))
            lower, upper = max(x[0] for x in intervals), min(x[1] for x in intervals)
            self.assertEqual(row["feasible"], lower <= upper)
            self.assertAlmostEqual(row["gap"], lower - upper)
            self.assertEqual(row["status"], "post_hoc_sensitivity_box_check")
        primary = data["scenarios"][0]
        alpha, alo, ahi = ANALYSIS.calibration(pairs, primary)
        for row in self.report["excluded_gene_sensitivity"]:
            fits = [f for f in self.fitness if f["gene"] == row["gene"]]
            points = [math.exp(alpha) * MODEL.exposure(float(f["estimate"]), 10) for f in fits]
            self.assertAlmostEqual(row["point_prediction_min"], min(points))
            self.assertAlmostEqual(row["point_prediction_max"], max(points))
            self.assertAlmostEqual(row["envelope_min"], min(math.exp(alo) * MODEL.exposure(float(f["lower"]), 10) for f in fits))
            self.assertAlmostEqual(row["envelope_max"], max(math.exp(ahi) * MODEL.exposure(float(f["upper"]), 10) for f in fits))

    def test_reference_groups_reconstruct_published_and_derived_odds(self):
        c = read_csv(ROOT / "data" / "cohorts_outcomes.csv")[0]
        total, cases, carriers, carrier_cases, sweeps, sweep_cases = (
            int(c[k]) for k in ("total_n", "total_cases", "carrier_n", "carrier_cases", "sweep_n", "sweep_cases"))
        all_comparison = (sweep_cases / (sweeps - sweep_cases)) / (
            (cases - sweep_cases) / ((total - sweeps) - (cases - sweep_cases)))
        carrier_comparison = (sweep_cases / (sweeps - sweep_cases)) / (
            (carrier_cases - sweep_cases) / ((carriers - sweeps) - (carrier_cases - sweep_cases)))
        report = self.report["reference_group_audit"]
        self.assertAlmostEqual(report["all_participant_or"], all_comparison)
        self.assertAlmostEqual(report["within_carriers_or"], carrier_comparison)
        self.assertAlmostEqual(report["comparator_or_factor"], all_comparison / carrier_comparison)
        source = next(r for r in self.outcomes if r["study_id"] == "watson2026_ukb")
        self.assertEqual(round(all_comparison), int(source["estimate"]))
        self.assertEqual(int(source["n_controls"]), total - cases)
        self.assertNotEqual(source["endpoint"], "AML")

    def test_analytic_branching_against_independent_backward_ode(self):
        rng = random.Random(260117)
        cases = [(0, 0, .03, 5), (0, .2, .03, 5), (.2, 0, .03, 5),
                 (.2, .2, .03, 5), (.1, .2, .03, 5), (.2, .1, 0, 5), (.2, .1, .03, 0)]
        cases += [(rng.random(), rng.random(), rng.random() / 5, rng.random() * 6) for _ in range(24)]
        for birth, death, hit, horizon in cases:
            exact = math.exp(MODEL.log_survival_one(birth, death, hit, horizon))
            numerical = numerical_survival(birth, death, hit, horizon)
            self.assertAlmostEqual(exact, numerical, delta=1e-10)
            for n in (0, 1, 7):
                self.assertAlmostEqual(MODEL.progression_probability(n, birth, death, hit, horizon), 1 - exact**n, delta=1e-12)
            if hit and horizon:
                eps = 1e-5
                derivative = (MODEL.cumulative_hazard(3, birth, death, hit, horizon + eps) -
                              MODEL.cumulative_hazard(3, birth, death, hit, horizon - eps)) / (2 * eps)
                self.assertAlmostEqual(MODEL.instantaneous_hazard(3, birth, death, hit, horizon), derivative, delta=1e-8)

    def test_second_moment_and_proved_error_bound(self):
        for birth, death, horizon in ((.1, .1, 3), (.1, .100000001, 3), (.2, .1, 4), (.1, .3, 4), (0, .2, 4)):
            s = birth - death
            first = simpson(lambda t: math.exp(s * t), horizon)
            def variance_integrand(r):
                variance_z = (birth + death) * (r if s == 0 else math.exp(s*r) * math.expm1(s*r) / s)
                tail = simpson(lambda u: math.exp(s * u), horizon - r, 80)
                return 2 * variance_z * tail
            second = first**2 + simpson(variance_integrand, horizon, 400)
            self.assertAlmostEqual(MODEL.integrated_second_moment(birth, death, horizon), second, delta=1e-7)
            for n in (1, 5, 100):
                for hit in (0, 1e-10, .01, .5):
                    exact = MODEL.cumulative_hazard(n, birth, death, hit, horizon)
                    upper = n * hit * first
                    lower = max(0, upper - n * hit**2 * second / 2)
                    self.assertLessEqual(exact, upper + 1e-9)
                    self.assertGreaterEqual(exact, lower - 1e-9)
        # Finite-risk limit: n*q fixed, and error shrinks as the proved 1/n bound.
        for n in (10, 100, 1000):
            target = .2 * MODEL.exposure(.05, 5)
            exact = MODEL.cumulative_hazard(n, .12, .07, .2 / n, 5)
            bound = .2**2 * MODEL.integrated_second_moment(.12, .07, 5) / (2 * n)
            self.assertLessEqual(abs(exact - target), bound + 1e-12)

    def test_logistic_exposure_kernel_bounds_and_sufficiency_counterexamples(self):
        for s in (-.2, -1e-10, 0, 1e-10, .2):
            a = simpson(lambda t: math.exp(s * t), 10)
            self.assertAlmostEqual(MODEL.exposure(s, 10), a, delta=1e-9)
            for x in (0, .04, .2, .8, 1):
                value = MODEL.logistic_exposure(s, x, 10)
                check = simpson(lambda t: math.exp(s*t) / (1 - x + x*math.exp(s*t)), 10)
                self.assertAlmostEqual(value, check, delta=1e-9)
                self.assertGreaterEqual(value, min(a, 10) - 1e-9)
                self.assertLessEqual(value, max(a, 10) + 1e-9)
        # Endpoint masses attain both bounds; arbitrary shared discrete kernels obey them.
        for points, weights in (([0], [1]), ([10], [1]), ([0, 2, 7, 10], [2, 3, 1, 4])):
            k1 = sum(w * math.exp(.2 * t) for t, w in zip(points, weights))
            k0 = sum(w * math.exp(.05 * t) for t, w in zip(points, weights))
            self.assertGreaterEqual(k1 / k0, 1)
            self.assertLessEqual(k1 / k0, math.exp(1.5) + 1e-12)
        self.assertLess(MODEL.cumulative_hazard(4, .4, .4, .1, 5), MODEL.cumulative_hazard(4, 0, 0, .1, 5))
        self.assertLess(MODEL.progression_probability(4, .12, .07, .03, 5), MODEL.progression_probability(4, .2, .07, .03, 5))

    def test_multistage_delay_hidden_states_and_marginal_mixture(self):
        rate, horizon = .13, 5
        for stages in (1, 2, 3, 5):
            cdf = 1 - math.exp(-rate*horizon) * sum((rate*horizon)**j / math.factorial(j) for j in range(stages))
            self.assertAlmostEqual(MODEL.erlang_completion(stages, rate, horizon), cdf, delta=1e-14)
            counts = [3] + [0] * (stages - 1)
            survival = MODEL.multistage_survival(counts, [0]*stages, [0]*stages, [rate]*stages, horizon, steps=1024)
            self.assertAlmostEqual(survival, (1-cdf)**3, delta=1e-11)
        lower = MODEL.multistage_survival([3, 0], [0, 0], [0, 0], [rate, rate], horizon, steps=1024)
        advanced = MODEL.multistage_survival([0, 3], [0, 0], [0, 0], [rate, rate], horizon, steps=1024)
        self.assertGreater(lower, advanced)
        n, weights = (1, 10, 100), (.2, .4, .4)
        for kl, kh in ((.01, .02), (.02, .2), (.2, 2)):
            def marginal(k):
                return -math.log(sum(w * math.exp(-count*k) for w, count in zip(weights, n)))
            self.assertGreaterEqual(marginal(kh)/marginal(kl), 1)
            self.assertLess(marginal(kh)/marginal(kl), kh/kl)

    def test_seeded_simulation_diagnostics_are_recomputed(self):
        with self.assertRaises(ValueError):
            SIMULATE.one_run(.5, .1, .1, .1, 1, random.Random(1))
        with self.assertRaises(ValueError):
            SIMULATE.one_run(1, math.inf, .1, .1, 1, random.Random(1))
        repeated = SIMULATE.run_scenarios(ROOT)
        self.assertEqual(repeated, self.report["simulations"])
        for row in repeated:
            self.assertLess(abs(row["standardized_error"]), 6)
            self.assertGreaterEqual(row["exact_cumulative_hazard"], row["rare_bound_lower"] - 1e-12)
            self.assertLessEqual(row["exact_cumulative_hazard"], row["rare_bound_upper"] + 1e-12)


if __name__ == "__main__":
    unittest.main()
