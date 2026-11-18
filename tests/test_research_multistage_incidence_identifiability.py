"""Offline, stdlib-only checks of the incidence identifiability research.

Nonclinical computational research; no medical advice. Finite checks and
independent simulations corroborate, but cannot replace, the README proofs.
"""

from __future__ import annotations

import csv
from decimal import Decimal, localcontext
from fractions import Fraction as F
import importlib.util
import json
from math import factorial
from pathlib import Path
import sys
import tempfile
import unittest


RESEARCH = Path(__file__).resolve().parents[1] / "research" / "multistage-incidence-identifiability"


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, RESEARCH / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


V = _load("research_msce_verifier", "verifier.py")
# Resolve the analysis CLI's local dependency without leaking a generic module
# name into the rest of the repository's unittest discovery process.
_prior_verifier = sys.modules.get("verifier")
sys.modules["verifier"] = V
try:
    A = _load("research_msce_analysis", "analysis.py")
finally:
    if _prior_verifier is None:
        del sys.modules["verifier"]
    else:
        sys.modules["verifier"] = _prior_verifier
S = _load("research_msce_simulation", "simulation.py")


class MultistageIncidenceIdentifiabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = A.analyze()
        cls.simulation = S.run_checks(quick=True)

    def test_01_exact_branching_onset_and_boundaries(self):
        checked = V.exact_checks()
        self.assertEqual(checked["onset_parameter_cases"], 128)
        for k in range(1, 13):
            rates = tuple(F(i + 1, 13) for i in range(k - 1))
            p = V.Branching(2, (0,) * (k - 1), (0,) * (k - 1),
                            rates, (True,) * (k - 1))
            jet = V.branching_jet(p, k + 2)
            product = F(2)
            for mu in rates:
                product *= mu
            self.assertEqual(jet[:k - 1], (F(0),) * (k - 1))
            self.assertEqual(jet[k - 1] * factorial(k - 1), product)
        with self.assertRaises(ValueError):
            V.Branching(1, (0,), (0,), (0,), (False,))
        with self.assertRaises(ValueError):
            V.Branching(1, (0,), (0,), (1,), (2,))

    def test_02_complete_tsce_invariants_and_feasible_class(self):
        # Fix (a,b,c)=(2,1,3/2), and vary the unidentified initiation rate.
        reference = None
        for lam in (F(2), F(3), F(7), F(13)):
            mu, birth = F(2) / lam, F(3, 2) * lam
            death = birth - mu - 1
            self.assertGreaterEqual(death, 0)
            p = V.Branching(lam, (birth,), (death,), (mu,), (False,))
            jet = V.branching_jet(p, 20)
            self.assertEqual(V.recover_tsce(jet), (F(2), F(1), F(3, 2)))
            if reference is None:
                reference = jet
            self.assertEqual(jet, reference)
        no_birth = V.Branching(2, (0,), (1,), (2,), (False,))
        self.assertEqual(V.recover_tsce(V.branching_jet(no_birth, 3)), (F(4), F(-3), F(0)))
        self.assertLess(self.simulation["tsce_symmetry"]["max_hazard_residual"], 1e-10)

    def test_03_rare_marking_scalar_inequality_and_spectral_symmetry(self):
        # Independent rational enclosures of 0<=x-(1-e^-x)<=x^2/2.
        for x in (F(0), F(1, 1000), F(1, 3), F(1), F(3)):
            lo, hi = V.exp_interval(-x, 100)
            self.assertGreaterEqual(x - 1 + lo, 0)
            self.assertLessEqual(x - 1 + hi, x * x / 2)
        growth = (F(1, 5), F(0), F(-2, 7), F(1, 5))
        self.assertEqual(V.linear_jet(F(3), growth, 22),
                         V.linear_jet(F(3), growth[::-1], 22))
        numerator, denominator = V.spectral_signature(F(3), growth)
        self.assertEqual(numerator, 3)
        self.assertEqual(len(denominator) - 1, 5)
        self.assertEqual(denominator[-1], 1)
        self.assertEqual(denominator[:2], (F(0), F(0)))
        for row in self.simulation["rare_marking_limit"]:
            for finite in row["finite_epsilon"]:
                self.assertLessEqual(finite["gap_to_mc_limit"], finite["moment_bound_mc"] + 1e-12)

    def test_04_exact_observation_pole_assignments(self):
        observed = V.spectral_signature(F(3), (F(1, 4), F(0), F(-2), F(-3)))
        counts = set()
        for mask in range(4):
            biological = [F(1, 4), F(0)]
            delays = []
            amplitude = F(3)
            for i, rate in enumerate((F(2), F(3))):
                if mask & (1 << i):
                    delays.append(rate)
                    amplitude /= rate
                else:
                    biological.append(-rate)
            counts.add(len(biological) + 1)
            self.assertEqual(V.spectral_signature(amplitude, tuple(biological), tuple(delays)), observed)
        self.assertEqual(counts, {3, 4, 5})
        # No biological intermediate stages: k=1, both poles are observation.
        self.assertEqual(V.spectral_signature(F(1, 2), (), (F(2), F(3))),
                         V.spectral_signature(F(3), (F(-2), F(-3))))
        values = [r["observed_hazard"] for r in self.simulation["delay_pole_assignments"]]
        self.assertLess(max(values) - min(values), 1e-9)

    def test_05_exact_converter_filter_and_deficit(self):
        p = V.Branching(3, (F(2, 5), F(1, 3)), (F(1, 7), F(1, 2)),
                        (F(2, 7), F(3, 11)), (False, True))
        jet = V.branching_jet(p, 24)
        speed = F(9, 2)
        extra = V.Branching(p.initiation, (F(0),) + p.birth,
                            (F(0),) + p.death, (speed,) + p.mutation,
                            (True,) + p.conversion)
        self.assertEqual(V.branching_jet(extra, 24), V.exponential_delay_jet(jet, speed))
        # k=1 identity: integrated deficit A*(1-exp(-LT))/L, strictly < A/L.
        lo, hi = V.exp_interval(F(-6), 100)
        deficit_lower, deficit_upper = F(3, 2) * (1 - hi), F(3, 2) * (1 - lo)
        self.assertGreater(deficit_lower, 0)
        self.assertLess(deficit_upper, F(3, 2))
        for row in self.simulation["fast_precursor"]:
            self.assertLessEqual(row["integrated_hazard_deficit"], row["bound_h_T_over_L"] + 1e-7)

    def test_06_poisson_affinity_from_discrete_probabilities_and_exact_risk(self):
        # sqrt(Pois(1)[z]*Pois(4)[z]) = exp(-5/2)*2^z/z!.
        # Summing the actual discrete coefficients independently checks the
        # closed-form affinity exp(-1/2), including its omitted-series tail.
        n = 45
        series = sum((F(2) ** z / factorial(z) for z in range(n + 1)), F(0))
        tail = F(2) ** (n + 1) / factorial(n + 1) / (1 - F(2, n + 2))
        elo, ehi = V.exp_interval(F(-5, 2), 120)
        rho_lower, rho_upper = elo * series, ehi * (series + tail)
        closed_lo, closed_hi = V.exp_interval(F(-1, 2), 120)
        self.assertLessEqual(rho_lower, closed_hi)
        self.assertGreaterEqual(rho_upper, closed_lo)
        self.assertLess(rho_upper - rho_lower, F(1, 10**40))
        # For means 1 vs 4, the likelihood ratio changes sign between 2 and 3.
        e3lo, e3hi = V.exp_interval(F(3), 120)
        self.assertGreater(e3lo, 16)
        self.assertLess(e3hi, 64)
        # Exact optimal equal-prior error = (1-2.5e^-1+13e^-4)/2.
        a, b = V.exp_interval(F(-1), 120)
        c, d = V.exp_interval(F(-4), 120)
        risk_lower = (1 - F(5, 2) * b + 13 * c) / 2
        risk_upper = (1 - F(5, 2) * a + 13 * d) / 2
        bounds = V.poisson_bounds_exact((F(1),), (F(4),))
        self.assertLess(bounds["affinity_bound_upper"], risk_lower)
        self.assertLess(bounds["coupling_bound_upper"], risk_lower)
        self.assertLess(risk_upper, F(1, 2))
        equal = V.poisson_bounds_exact((F(0), F(2)), (F(0), F(2)))
        self.assertEqual(equal["coupling_bound_lower"], F(1, 2))
        self.assertGreater(equal["affinity_bound_lower"], F(499999, 1000000))
        zero = V.poisson_bounds_exact((F(0),), (F(1),))
        self.assertGreater(zero["coupling_bound_lower"], 0)
        experiment = self.simulation["poisson_decision"]
        self.assertGreaterEqual(experiment["equal_prior_error"] + experiment["six_standard_error_upper_bound"],
                                experiment["hellinger_lower_bound"])

    def test_07_calibration_is_sharp_and_positive_tolerance_hides_fast_steps(self):
        y, shape = (F(2), F(8), F(4)), (F(1), F(1), F(1))
        self.assertEqual(V.calibration_ratio(y, y, shape), 4)
        self.assertEqual(V.amplitude_interval(y, y, shape, F(2)), (F(4), F(4)))
        for rho in (F(3, 2), F(19999, 10000)):
            left, right = V.amplitude_interval(y, y, shape, rho)
            self.assertGreater(left, right)
        self.assertEqual(V.calibration_ratio((F(1), F(2)), (F(3), F(5)), (F(1), F(2))), 1)
        self.assertEqual(V.flexible_shape_calibration(y, y, ((F(1), F(4)),) * 3), 1)
        # Exact general deficit logic of (13), avoiding a fitted toy curve:
        for rho in (F(1001, 1000), F(11, 10), F(2)):
            means, wh = (F(1, 4), F(2), F(9)), F(5)
            speed = wh / (min(means) * (1 - 1 / rho ** 2))
            delta = wh / speed
            worst_shapes = tuple(m - delta for m in means)
            self.assertGreater(min(worst_shapes), 0)
            self.assertLessEqual(V.calibration_ratio(means, means, worst_shapes), rho ** 2)
        self.assertTrue(all(not r["factor_below_threshold_is_feasible"]
                            for r in self.simulation["bounded_observation_multiplier"]))

    def test_08_monotone_falsification_uses_actual_source_values(self):
        groups = A.load_data()
        breast = groups["female_breast"]
        self.assertEqual([r["rate_per_100000"] for r in breast[-3:]], ["470.2", "458.9", "413.5"])
        mono = self.results["cancers"]["female_breast"]["monotonicity"]
        self.assertFalse(mono["homogeneous_monotone_mean_compatible"])
        self.assertEqual(mono["minimum_rho_for_any_nondecreasing_mean"], "1.066238")
        for site in ("pancreas", "colon_and_rectum", "lung_and_bronchus"):
            self.assertTrue(self.results["cancers"][site]["monotonicity"]["homogeneous_monotone_mean_compatible"])
        # Construct a nondecreasing envelope at a rational factor above the
        # exact threshold; its intervals must intersect all source intervals.
        lower, upper = A._rounding_intervals(breast)
        rho = F(1067, 1000)
        previous = F(0)
        for l, u in zip(lower, upper):
            previous = max(previous, l / rho)
            self.assertLessEqual(previous, rho * u)

    def test_09_independent_closed_form_for_signed_bin_series(self):
        self.assertIsInstance(V.power_bin(3, 30, 35), F)
        self.assertEqual(V.power_bin(3, 30, 35), V.power_bin(3, F(30), F(35)))
        with localcontext() as ctx:
            ctx.prec = 180
            def dec(x):
                return Decimal(x.numerator) / Decimal(x.denominator)
            # Integrate the alternative exp-minus-polynomial formula rather
            # than reimplementing the same Taylor series as the verifier.
            for k in (2, 3, 6, 10):
                for g in (F(-3, 50), F(-1, 200), F(1, 200), F(3, 50)):
                    for a, b in ((F(0), F(5)), (F(80), F(85))):
                        lo, hi = V.one_growth_bin_interval(k, g, a, b, 40)
                        gd, ad, bd = dec(g), dec(a), dec(b)
                        exp_integral = ((gd * bd).exp() - (gd * ad).exp()) / gd
                        polynomial = sum((gd ** j * (bd ** (j + 1) - ad ** (j + 1)) /
                                          Decimal(factorial(j + 1)) for j in range(k - 1)), Decimal(0))
                        independent = (exp_integral - polynomial) / (gd ** (k - 1) * (bd - ad))
                        self.assertLessEqual(dec(lo), independent)
                        self.assertGreaterEqual(dec(hi), independent)

    def test_10_data_provenance_and_deterministic_full_reproduction(self):
        self.assertEqual(self.results["source_sha256"], "1a0b4c66b6267c9841b05dfb47cf2831e66e9d02a03bc3d3017b30cfc4fa0198")
        # Git core.autocrlf=true may produce CRLF, or an older checkout CR.
        # A provenance digest must not depend on platform-specific newlines.
        canonical = A.DATA.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
        with tempfile.TemporaryDirectory() as temporary:
            for idx, ending in enumerate((b"\n", b"\r\n", b"\r")):
                alternate = Path(temporary) / f"alternate-newlines-{idx}.csv"
                alternate.write_bytes(canonical.replace(b"\n", ending))
                self.assertEqual(A.source_sha256(alternate), self.results["source_sha256"])
                self.assertEqual(sum(map(len, A.load_data(alternate).values())), 44)
        stored = json.loads((RESEARCH / "analysis_results.json").read_text(encoding="utf-8"))
        self.assertEqual(self.results, stored)
        groups = A.load_data()
        self.assertEqual(sum(map(len, groups.values())), 44)
        sums = {site: str(sum((F(r["rate_per_100000"]) for r in rows), F(0)))
                for site, rows in groups.items()}
        self.assertEqual(sums, {"pancreas": "733/2", "colon_and_rectum": "4873/5",
                               "lung_and_bronchus": "16197/10", "female_breast": "15192/5"})
        for rows in groups.values():
            for r in rows:
                self.assertEqual(r["access_date"], "2026-11-17")
                self.assertEqual(r["reuse_terms_url"], "https://seer.cancer.gov/archive/csr/1975_2017/")
                self.assertEqual((r["case_count"], r["person_years"]), ("", ""))
                self.assertIn("seer.cancer.gov/archive/csr/1975_2017/", r["source_pdf_url"])
        for record in self.results["cancers"].values():
            self.assertIsNone(record["identified_biological_k"])
            self.assertTrue(record["joint_grid_k_minimum_certified_unique"])
            self.assertTrue(all(row["grid_minimum_certified_unique"] for row in record["by_k"]))

    def test_11_reject_overlaps_and_inconsistent_source_strata(self):
        with A.DATA.open(newline="", encoding="utf-8") as stream:
            original = list(csv.DictReader(stream))
        with tempfile.TemporaryDirectory(dir=RESEARCH) as temporary:
            path = Path(temporary) / "bad.csv"
            for bad in ("overlap", "metadata", "label", "invented_denominator"):
                rows = [dict(r) for r in original]
                if bad == "overlap":
                    rows[1].update(age_start="32", age_end_exclusive="37", age_group="32-36")
                elif bad == "metadata":
                    rows[1]["sex"] = "female"
                elif bad == "label":
                    rows[1]["age_group"] = "40-44"
                else:
                    rows[1]["person_years"] = "100000"
                with path.open("w", newline="", encoding="utf-8") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(original[0]))
                    writer.writeheader()
                    writer.writerows(rows)
                with self.assertRaises(ValueError, msg=bad):
                    A.load_data(path)

    def test_12_synthetic_error_bounds_are_downward_rounded(self):
        for record in self.results["synthetic_poisson_exposure_certificates"]:
            speed = F(record["extra_conversion_rate_per_year"])
            lo, hi = V.exp_interval(-F(100) / speed, 100)
            printed = F(record["coupling_error_lower_bound"])
            self.assertLessEqual(printed, lo / 2)
            self.assertLess(hi / 2 - printed, F(1, 10**8))

    def test_13_independent_simulation_covers_all_theorem_families(self):
        self.assertTrue(self.simulation["passed"])
        self.assertEqual(self.simulation["seed"], 19270617)
        self.assertEqual(len(self.simulation["event_level_branching"]), 6)
        self.assertEqual(len(self.simulation["onset_order"]), 6)
        self.assertEqual(len(self.simulation["rare_marking_limit"]), 2)
        self.assertEqual(len(self.simulation["spectral_limit"]), 4)
        self.assertEqual(len(self.simulation["fast_precursor"]), 3)
        self.assertEqual(self.simulation["poisson_decision"]["trials_per_hypothesis"], 5000)


if __name__ == "__main__":
    unittest.main()
