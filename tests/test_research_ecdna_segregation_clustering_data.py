"""Nonclinical aggregate ecDNA calibration; independent stdlib-only checks."""

import ast
import csv
from fractions import Fraction
import importlib.util
import math
from pathlib import Path
import sys
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "research" / "ecdna-segregation-clustering-data"
F = Fraction


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, DIRECTORY / filename)
    instance = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(instance)
    return instance


def independent_bb(n, p, rho):
    """Separate rising-factorial calculation of a small rational kernel."""
    if rho == 0:
        return [F(math.comb(n, j))*p**j*(1-p)**(n-j) for j in range(n+1)]
    a, b = p*(1/rho-1), (1-p)*(1/rho-1)
    denominator = math.prod(a+b+i for i in range(n))
    return [F(math.comb(n, j))*math.prod(a+i for i in range(j))
            * math.prod(b+i for i in range(n-j))/denominator for j in range(n+1)]


def lgamma_observation_distribution(R, B, log_size, log_mean, maximum=100):
    """Independent synthetic observation law; never source observations.

    Direct lgamma probabilities, all ordered allocations and explicit sorting
    exercise zero floors and ties independently of the production score code.
    """
    p = (1-math.sqrt(B))/2
    rho = (R-B)/(1-B)
    a, b = p*(1/rho-1), (1-p)*(1/rho-1)
    size, mean = math.exp(log_size), math.exp(log_mean)
    q0 = (size/(size+mean))**size
    out = {}
    logbeta = math.lgamma(a)+math.lgamma(b)-math.lgamma(a+b)
    for n in range(1, maximum+1):
        ln_nb = (math.lgamma(size+n)-math.lgamma(size)-math.lgamma(n+1)
                 + size*math.log(size/(size+mean))+n*math.log(mean/(size+mean))
                 - math.log1p(-q0))
        for j in range(n+1):
            ln_bb = (math.lgamma(n+1)-math.lgamma(j+1)-math.lgamma(n-j+1)
                     + math.lgamma(a+j)+math.lgamma(b+n-j)-math.lgamma(a+b+n)
                     - logbeta)
            key = tuple(sorted((max(j, 1), max(n-j, 1))))
            out[key] = out.get(key, 0.0)+math.exp(ln_nb+ln_bb)
    return out


class TestEcDNASegregationClusteringData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.start = time.perf_counter()
        cls.analysis = module("segregation_data_analysis_test", "analysis.py")
        cls.verifier = cls.analysis.verifier_module()
        cls.result = cls.analysis.results()
        cls.result_seconds = time.perf_counter()-cls.start

    def test_every_reported_number_recomputed_from_csv(self):
        text = (DIRECTORY / "README.md").read_text(encoding="utf-8")
        actual = text.split("<!-- BEGIN RECOMPUTED RESULTS -->\n", 1)[1].split(
            "<!-- END RECOMPUTED RESULTS -->", 1)[0]
        self.assertEqual(actual, self.analysis.render_report()+"\n")
        self.assertLess(self.result_seconds, 60)

    def test_source_provenance_all_rows_and_no_raw_payload(self):
        expected = {"shibata_fits.csv", "shibata_models.csv", "lundberg.csv",
                    "other_aggregates.csv", "literature_audit.csv", "analysis_design.csv"}
        self.assertEqual({p.name for p in (DIRECTORY / "data").iterdir()}, expected)
        for path in (DIRECTORY / "data").glob("*.csv"):
            self.assertLess(path.stat().st_size, 100_000)
            rows = self.analysis.load_csv(path.name)
            self.assertTrue(rows)
            for row in rows:
                for key in ("source_url", "locator", "access_date", "reuse_basis"):
                    self.assertTrue(row[key].strip(), (path.name, key, row))
                self.assertEqual(row["access_date"], "2026-11-17")
                if path.name != "analysis_design.csv":
                    self.assertTrue(row["doi"].startswith("10."))
                    self.assertTrue(row["source_url"].startswith("https://"))
                if "value" in row and path.name != "analysis_design.csv":
                    F(row["value"])
            id_field = "key" if path.name == "analysis_design.csv" else "id"
            self.assertEqual(len({row[id_field] for row in rows}), len(rows))
        allowed = {"patient_id", "patient_name", "cell_id", "genome", "sequence"}
        for path in (DIRECTORY / "data").glob("*.csv"):
            self.assertFalse(allowed.intersection(self.analysis.load_csv(path.name)[0]))
        self.assertFalse(any(DIRECTORY.rglob("*.xlsx")))
        self.assertFalse(any(DIRECTORY.rglob("*.pdf")))

    def test_stdlib_only_and_no_runtime_retrieval(self):
        for filename in ("analysis.py", "verifier.py"):
            tree = ast.parse((DIRECTORY / filename).read_text())
            for node in ast.walk(tree):
                names = []
                if isinstance(node, ast.Import):
                    names = [alias.name.split(".")[0] for alias in node.names]
                elif isinstance(node, ast.ImportFrom):
                    names = [(node.module or "").split(".")[0]]
                for name in names:
                    self.assertIn(name, sys.stdlib_module_names)
                    self.assertNotIn(name, {"urllib", "http", "socket", "subprocess"})

    def test_all_clones_controls_and_conditions_preserved(self):
        rows = self.analysis.load_csv("shibata_fits.csv")
        clones = {row["id"] for row in rows if row["id"].startswith("shibata_clone")}
        self.assertEqual(clones, {f"shibata_clone{i}" for i in range(1, 13)})
        controls = {row["id"] for row in rows if row["scope"] == "HSR_control"}
        self.assertEqual(controls, {"shibata_clone4", "shibata_clone6", "shibata_clone7"})
        self.assertEqual(len(rows)-len(clones), 4)
        self.assertTrue(controls.isdisjoint(row["id"] for row in self.result["fits"]))
        self.assertEqual(self.result["primary_summary"]["n_fits"], 9)
        audit = {row["id"]: row for row in self.analysis.load_csv("literature_audit.csv")}
        for seed in ("lundberg_2008", "lange_2022", "shibata_2026", "scanu_2025",
                     "sankar_2026", "nichols_2025", "bruckner_2026"):
            self.assertIn(seed, audit)
        self.assertIn("preprint", audit["scanu_2025"]["publication_status"])

    def test_exact_kernels_and_conditional_moments(self):
        for row in self.analysis.load_csv("shibata_fits.csv"):
            p, rho = F(row["p"]), F(row["rho"])
            R = self.verifier.effective_r(p, rho)
            for k in (1, 2, 7):
                nums, den = self.verifier.integer_kernel(k, p, rho)
                self.assertEqual(sum(nums), den)
                symmetric = [F(nums[j]+nums[2*k-j], 2*den) for j in range(2*k+1)]
                self.assertEqual(sum(symmetric), 1)
                self.assertEqual(sum(j*q for j, q in enumerate(symmetric)), k)
                variance = sum((j-k)**2*q for j, q in enumerate(symmetric))
                self.assertEqual(variance, R*k*k+(1-R)*k/2)
                self.assertEqual(variance/k, (1-R)/2+R*k)
        for p, rho in ((F(2, 5), F(1, 7)), (F(1, 2), F(0))):
            nums, den = self.verifier.integer_kernel(4, p, rho)
            self.assertEqual([F(num, den) for num in nums], independent_bb(8, p, rho))

    def test_fold_alias_and_R_does_not_determine_kernel(self):
        R = F(1, 5)
        for k in (1, 3, 8):
            one, den1 = self.verifier.integer_kernel(k, F(1, 2), R)
            two, den2 = self.verifier.integer_kernel(k, (1-R)/2, R/(1+R))
            self.assertEqual([F(one[j]+one[2*k-j], den1) for j in range(2*k+1)],
                             [F(two[j]+two[2*k-j], den2) for j in range(2*k+1)])
        self.assertEqual(self.verifier.effective_r(F(1, 2), F(1, 4)),
                         self.verifier.effective_r(F(3, 4), F(0)))
        a = self.verifier.threshold(F(1, 2), F(1, 4))
        b = self.verifier.threshold(F(3, 4), F(0))
        self.assertLess(a["upper"], b["lower"])

    def test_all_point_certificates_rechecked_and_rounding_encloses(self):
        originals = {row["id"]: row for row in self.analysis.load_csv("shibata_fits.csv")}
        for row in self.result["fits"]:
            original = originals[row["id"]]
            p, rho = F(original["p"]), F(original["rho"])
            point = row["point_threshold"]
            cap = point["cap"]
            kernels = tuple(self.verifier.integer_kernel(k, p, rho) for k in range(1, cap+1))
            self.assertIsNotNone(self.verifier._lower_witness(kernels, F(point["lower_exact"]), cap))
            self.assertIsNotNone(self.verifier._upper_witness(kernels, F(point["upper_exact"]), cap))
            rounding = row["rounding_threshold"]
            self.assertLessEqual(F(rounding["lower_exact"]), F(point["lower_exact"]))
            self.assertGreaterEqual(F(rounding["upper_exact"]), F(point["upper_exact"]))
            self.assertLessEqual(point["width"], float(self.analysis.settings()["threshold_width"]))
            self.assertLess(rounding["upper"], self.result["baseline"]["lower"])

    def test_tiny_parameter_box_contains_all_kernel_entries(self):
        row = self.analysis.load_csv("shibata_fits.csv")[0]
        p, rho = F(row["p"]), F(row["rho"])
        eps = F(1, 200_000_000)
        lower = self.verifier._small_box_rows(p-eps, p+eps, rho-eps, rho+eps, 4, False)
        upper = self.verifier._small_box_rows(p-eps, p+eps, rho-eps, rho+eps, 4, True)
        for pp in (p-eps, p, p+eps):
            for rr in (rho-eps, rho, rho+eps):
                for k in range(1, 5):
                    nums, den = self.verifier.integer_kernel(k, pp, rr)
                    low, ld = lower[k-1]
                    high, hd = upper[k-1]
                    for j, value in enumerate(nums):
                        self.assertLessEqual(low[j]*den, value*ld)
                        self.assertLessEqual(value*hd, high[j]*den)
        # The alternative Polya-factor bounds also cover a broad interior box.
        lo, ld = self.verifier.interval_kernel(3, F(2, 5), F(1, 2), F(0), F(1, 5), False)
        hi, hd = self.verifier.interval_kernel(3, F(2, 5), F(1, 2), F(0), F(1, 5), True)
        for pp, rr in ((F(2, 5), F(0)), (F(9, 20), F(1, 10)), (F(1, 2), F(1, 5))):
            nums, den = self.verifier.integer_kernel(3, pp, rr)
            for j, value in enumerate(nums):
                self.assertLessEqual(lo[j]*den, value*ld)
                self.assertLessEqual(value*hd, hi[j]*den)

    def test_false_truncation_certificate_is_rejected(self):
        rows = (self.verifier.integer_kernel(1, F(1, 2), F(0)),)
        # Type-one principal mean is <1 even where the infinite process persists.
        self.assertFalse(self.verifier._certificate(rows, F(1, 10), (1,), tail=1))
        self.assertFalse(self.verifier._certificate(rows, F(1, 10), (0,), tail=1))
        self.assertGreater(self.result["baseline"]["lower"], 0.1)
        for p, rho in ((0, F(1, 3)), (1, F(0)), (F(1, 2), 1)):
            self.assertEqual(self.verifier.threshold(p, rho)["upper"], 0)

    def test_uniform_moment_bound_algebra_and_coverage_containment(self):
        for R in (F(0), F(1, 10), F(1, 5), F(1, 3), F(1, 2), F(1)):
            bound = self.verifier.moment_threshold_bounds(R, R)
            self.assertEqual(F(bound["upper_exact"]), (1-R)/2)
            c = F(bound["lower_exact"])/2
            if c == 0:
                continue
            for T in (R*R, (R*R+R)/2, R):
                aa = (1-R)/(1+c)
                bb = (1+R)/(2*(1+c))
                cc = (1-T)/(2*(1+2*c))
                dd = 3*(1-2*R+T)/(4*(1+2*c))
                polynomial = (c+R)*(1+8*c+6*R-3*T)-(1+R)*(1-T)
                self.assertEqual(((1-aa)*(1-dd)-bb*cc)*4*(1+c)*(1+2*c), polynomial)
                self.assertLess(polynomial, 0)
        for row in self.result["fits"]:
            band, bound = row["outward_sensitivity_R"], row["propagated_sensitivity"]
            self.assertLessEqual(band[0], row["R"])
            self.assertLessEqual(row["R"], band[1])
            self.assertLessEqual(bound["lower"], row["point_threshold"]["lower"])
            self.assertGreaterEqual(bound["upper"], row["point_threshold"]["upper"])
            self.assertLessEqual(bound["lower"], self.result["baseline"]["lower"])
            self.assertGreaterEqual(bound["upper"], self.result["baseline"]["upper"])
        self.assertEqual(self.result["unidentified_R_threshold"]["lower"], 0)
        self.assertEqual(self.result["unidentified_R_threshold"]["upper"], 0.5)

    def test_rigorous_enclosures_display_outward(self):
        for row in self.result["fits"]:
            for key in ("rounding_threshold", "propagated_sensitivity"):
                interval = row[key]
                display = self.analysis.display_enclosure(interval)
                low, high = map(F, display.strip("[]").split(", "))
                self.assertLessEqual(low, F(interval["lower_exact"]))
                self.assertGreaterEqual(high, F(interval["upper_exact"]))
                base = self.result["baseline"]
                shift = self.analysis.display_relative_enclosure(interval, base, 4)
                low, high = map(F, shift.strip("[]").split(", "))
                self.assertLessEqual(low, 100*(F(interval["lower_exact"])/F(base["upper_exact"])-1))
                self.assertGreaterEqual(high, 100*(F(interval["upper_exact"])/F(base["lower_exact"])-1))
        self.assertEqual(self.analysis.outward(F(-1, 3), 2), "-0.34")
        self.assertEqual(self.analysis.outward(F(-1, 3), 2, True), "-0.33")

    def test_normalized_information_against_independent_finite_differences(self):
        p, rho, size, mean = 0.4, 0.2, 2.0, 3.0
        R, B = self.analysis.effective_r(p, rho), (2*p-1)**2
        theta = [R, B, math.log(size), math.log(mean)]
        base = lgamma_observation_distribution(*theta)
        self.assertAlmostEqual(sum(base.values()), 1, places=11)
        epsilon = 1e-5
        derivatives = []
        for dimension in range(4):
            plus, minus = theta[:], theta[:]
            plus[dimension] += epsilon
            minus[dimension] -= epsilon
            hi = lgamma_observation_distribution(*plus)
            lo = lgamma_observation_distribution(*minus)
            derivatives.append({key: (hi[key]-lo[key])/(2*epsilon) for key in base})
        expected = [[math.fsum(derivatives[i][key]*derivatives[j][key]/prob
                              for key, prob in base.items()) for j in range(4)] for i in range(4)]
        actual = self.analysis.design_information(p, rho, size, mean)
        for i in range(4):
            for j in range(4):
                self.assertAlmostEqual(actual["information"][i][j], expected[i][j], delta=2e-6)
        self.assertAlmostEqual(actual["included_mass"], 1, delta=2e-11)
        self.assertLess(max(map(abs, actual["score_sum"])), 1e-7)
        self.assertGreater(base[(1, 1)], 0)

    def test_sensitivity_not_mislabeled_as_empirical_confidence(self):
        for row in self.result["fits"]:
            info = row["local_information"]
            self.assertIn("no validated sampling coverage", info["label"])
            self.assertNotIn("empirical_p", info)
            self.assertLess(info["omitted_nb_mass"], 2e-11)
        clone10 = next(row for row in self.result["fits"] if row["id"] == "shibata_clone10")
        self.assertGreater(clone10["local_information"]["bias2_se"], 1)
        readme = (DIRECTORY / "README.md").read_text()
        self.assertIn("stronger requested empirical discovery is not established", readme)
        self.assertIn("post-freeze local Fisher sensitivity", readme)

    def test_lundberg_counts_exact_tails_and_dependence_limit(self):
        for row in self.result["lundberg"]["lines"]:
            n, count = row["n"], min(row["compatible_counts"])
            exact = sum(F(math.comb(n, j))*F(1, 100)**j*F(99, 100)**(n-j)
                        for j in range(count, n+1))
            self.assertAlmostEqual(row["conditional_p"], float(exact), places=13)
            self.assertGreaterEqual(row["markov_p_bound"], 0.05)
        self.assertEqual(self.analysis.compatible_counts(100, "2"), [2])
        self.assertEqual(self.analysis.compatible_counts(50, "2"), [1])
        self.assertEqual(self.analysis.compatible_counts(1000, "2"), list(range(15, 26)))
        self.assertEqual(self.result["lundberg"]["combined_count"], 13)

    def test_source_AIC_and_negative_results_preserved(self):
        rows = self.analysis.load_csv("shibata_models.csv")
        for row in rows:
            residual = F(row["aic"])-(-2*F(row["reported_negative_log_likelihood"])
                                      + 2*int(row["n_parameters"]))
            self.assertLessEqual(abs(residual), F(15, 10_000))
        best = min(self.result["source_models"], key=lambda x: x["aic"])
        self.assertEqual(best["model"], "M4")
        lange = [row for row in self.analysis.load_csv("other_aggregates.csv")
                 if row["study"] == "Lange et al. 2022" and row["quantity"] == "source_ks_p"]
        self.assertEqual(len(lange), 8)
        self.assertTrue(all(float(row["value"]) > 0.05 for row in lange))

    def test_future_prediction_from_exact_kernel_and_exact_null(self):
        design = self.analysis.settings()
        future = self.result["prediction"]
        n, cut = future["total"], future["minority_cutoff"]
        q0 = F(sum(math.comb(n, j) for j in range(n+1) if min(j, n-j) <= cut), 2**n)
        self.assertEqual(float(q0), future["binomial_extreme_probability"])
        source = next(row for row in self.analysis.load_csv("shibata_fits.csv")
                      if row["id"] == future["source_id"])
        mass = independent_bb(n, F(source["p"]), F(source["rho"]))
        q1 = sum(mass[j] for j in range(n+1) if min(j, n-j) <= cut)
        self.assertEqual(float(q1), future["fitted_extreme_probability"])
        m, critical = future["independent_pairs"], future["critical_count"]
        exact_tail = sum(F(math.comb(m, j))*q0**j*(1-q0)**(m-j)
                         for j in range(critical, m+1))
        previous = exact_tail+F(math.comb(m, critical-1))*q0**(critical-1)*(1-q0)**(m-critical+1)
        self.assertLessEqual(exact_tail, F(design["prediction_test_alpha"]))
        self.assertGreater(previous, F(design["prediction_test_alpha"]))
        self.assertAlmostEqual(float(exact_tail), future["null_rejection_probability"], places=14)
        lower, upper = future["conditional_prediction_interval"]
        probability_inside = sum(F(math.comb(m, j))*q1**j*(1-q1)**(m-j)
                                 for j in range(lower, upper+1))
        self.assertGreaterEqual(probability_inside, 1-F(design["nominal_alpha"]))
        self.assertGreater(future["conditional_power"], 0.99)


if __name__ == "__main__":
    unittest.main()
