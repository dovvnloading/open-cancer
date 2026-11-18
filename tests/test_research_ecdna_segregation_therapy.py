"""Nonclinical ecDNA mathematics: exact independent checks and seeded simulation.

No patient data or treatment advice. Finite tests do not prove an infinite-state
theorem; README.md supplies the proofs and identifies each observable.
"""

from collections import Counter, defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
import importlib.util
import math
from pathlib import Path
import random
import unittest


RESEARCH = Path(__file__).resolve().parents[1] / "research" / "ecdna-segregation-therapy"


def load_module(name):
    spec = importlib.util.spec_from_file_location("ecdna_" + name, RESEARCH / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


exact = load_module("verifier")
simulation = load_module("simulation")
F = Fraction


@lru_cache(maxsize=None)
def enumerate_offspring(k, u, s):
    """Independent bit-mask enumeration, without the verifier's binomial kernel."""
    probability = u * s**k
    result = defaultdict(Fraction)
    result[()] += 1 - probability
    # Only small generations use this enumerator. No model is truncated.
    for bits in range(1 << (2 * k)):
        left = bits.bit_count()
        children = tuple(sorted((left, 2 * k - left)))
        result[children] += probability / (1 << (2 * k))
    return dict(result)


def population_law(initial, generations, u, s):
    """Enumerate whole populations, retaining sisters and extinct populations."""
    law = {tuple(sorted(initial)): F(1)}
    for _ in range(generations):
        next_law = defaultdict(Fraction)
        for state, state_probability in law.items():
            children_law = {(): F(1)}
            for k in state:
                expanded = defaultdict(Fraction)
                for accumulated, accumulated_probability in children_law.items():
                    for children, probability in enumerate_offspring(k, u, s).items():
                        expanded[tuple(sorted(accumulated + children))] += (
                            accumulated_probability * probability
                        )
                children_law = expanded
            for children, probability in children_law.items():
                next_law[children] += state_probability * probability
        law = dict(next_law)
    return law


def decimal_positive_expectation(initial, n, u, s):
    """Independent high-precision finite-horizon evaluation; not a certificate."""
    with localcontext() as context:
        context.prec = 80
        sd = Decimal(s.numerator) / Decimal(s.denominator)
        ud = Decimal(u.numerator) / Decimal(u.denominator)
        a, b = Decimal(1), Decimal(0)
        for _ in range(n):
            a, b = sd * (1 + a)**2 / 4, sd * (1 + b)**2 / 4
        return float((2 * ud)**n * sum((a**k - b**k for k in initial), Decimal(0)))


class EcDNASegregationTherapyTests(unittest.TestCase):
    """Tests are deterministic, offline, stdlib-only, and bounded in runtime."""

    def assert_mc(self, observed, expected, standard_error, label):
        self.assertIsNotNone(observed, label)
        self.assertIsNotNone(standard_error, label)
        # A fixed six-standard-error check is a falsification screen, not a
        # confidence claim or proof. Numeric allowance covers float conversion.
        self.assertLessEqual(
            abs(observed - float(expected)),
            6 * standard_error + 1e-10 * max(1.0, abs(float(expected))),
            label,
        )

    def test_standalone_exact_verifier(self):
        report = exact.verify()
        self.assertGreaterEqual(report["total_checks"], 700)

    def test_exact_binomial_and_complementary_sisters_at_large_types(self):
        for k in (0, 1, 2, 5, 17, 64, 137):
            with self.subTest(k=k):
                probabilities = exact.binomial_kernel(k)
                self.assertEqual(sum(probabilities), 1)
                mean = sum(j * p for j, p in enumerate(probabilities))
                variance = sum((j-k)**2 * p for j, p in enumerate(probabilities))
                covariance = sum((j-k)*(k-j)*p for j, p in enumerate(probabilities))
                self.assertEqual(mean, k)
                self.assertEqual(variance, F(k, 2))
                self.assertEqual(covariance, -F(k, 2))
                self.assertEqual(
                    sum((j + (2*k-j))**2*p for j, p in enumerate(probabilities)),
                    4*k*k,
                )

    def test_complete_joint_population_laws_extinction_and_all_raw_moments(self):
        u, s = F(3, 4), F(4, 5)
        statistics = {
            "cells": len,
            "zero": lambda state: state.count(0),
            "positive": lambda state: sum(k > 0 for k in state),
            "copies": sum,
            "copy_squares": lambda state: sum(k*k for k in state),
        }
        for initial, n in (((), 2), ((0,), 2), ((1,), 3), ((2,), 2), ((0, 1), 2), ((1, 2), 2)):
            with self.subTest(initial=initial, n=n):
                law = population_law(initial, n, u, s)
                self.assertEqual(sum(law.values()), 1)
                extinction = sum(p for state, p in law.items() if not any(state))
                self.assertEqual(extinction, exact.population_extinction(initial, n, u, s))
                counts = defaultdict(Fraction)
                for state, probability in law.items():
                    for k, count in Counter(state).items():
                        counts[k] += probability * count
                self.assertEqual(
                    exact.summary_from_counts(counts),
                    exact.pgf_summary(initial, n, u, s),
                )
                for name, statistic in statistics.items():
                    mean = sum(p*statistic(state) for state, p in law.items())
                    second = sum(p*statistic(state)**2 for state, p in law.items())
                    self.assertEqual(
                        exact.population_raw_moments(initial, n, u, s, name),
                        {"mean": mean, "second": second, "variance": second-mean**2},
                    )

    def test_independent_sister_shortcut_is_detectably_wrong(self):
        # Same daughter marginals can preserve means while changing extinction.
        # Actual sisters from k=1 cannot both have zero copies.
        u, s = F(4, 5), F(3, 4)
        true = exact.extinction_by_generation(1, 1, u, s)
        wrong_independent_sisters = 1 - u*s + u*s*F(1, 4)**2
        self.assertEqual(true, 1-u*s)
        self.assertNotEqual(true, wrong_independent_sisters)
        self.assertEqual(exact.extinction_by_generation(1, 2, 1, F(3, 4)), F(7, 16))

    def test_seeded_exact_copy_mass_conservation_without_killing(self):
        for k in (1, 3, 10):
            observed = simulation.run_ensemble((k,), 2, 1.0, 1.0,
                                               trials=1000, seed=501+k)
            self.assertEqual(observed["means"]["copy_sum"], 4*k)
            self.assertEqual(observed["sample_variances"]["copy_sum"], 0)

    def test_endpoints_and_neutral_moment_law(self):
        for u, s in ((0, 0), (0, 1), (1, 0), (1, 1)):
            for initial in ((), (0,), (1,), (0, 2, 3)):
                for n in range(4):
                    self.assertEqual(
                        exact.pgf_summary(initial, n, u, s),
                        exact.summary_from_counts(exact.forward_expectations(initial, n, u, s)),
                    )
                    if n and (u == 0 or s == 0):
                        self.assertEqual(exact.population_extinction(initial, n, u, s), 1)
                    if u == s == 1 and any(initial):
                        self.assertEqual(exact.population_extinction(initial, n, u, s), 0)
        initial = (0, 1, 3, 8)
        mean = F(sum(initial), len(initial))
        variance = F(sum(k*k for k in initial), len(initial)) - mean**2
        for u in (F(1, 3), F(1, 2), F(1)):
            for n in range(7):
                result = exact.pgf_summary(initial, n, u, 1)
                self.assertEqual(result["pooled_mean"], mean)
                self.assertEqual(result["pooled_variance"], variance + n*mean/2)
        self.assertIsNone(exact.pgf_summary((), 3, 1, F(3, 4))["pooled_mean"])
        self.assertIsNone(exact.pgf_summary((0,), 3, 1, F(3, 4))["positive_mean"])

    def test_critical_boundary_exact_sign_and_general_copy_certificates(self):
        for u in (F(1), F(3, 4), F(2, 3), F(51, 100)):
            critical_s = (4*u-1)/(4*u*u)
            epsilon = (1-critical_s)/10**10
            self.assertEqual(exact.threshold_regime(u, critical_s), "critical")
            self.assertEqual(exact.threshold_regime(u, critical_s-epsilon), "eradication")
            self.assertEqual(exact.threshold_regime(u, critical_s+epsilon), "persistence")
            a = 1/(4*u-1)
            self.assertEqual(exact.threshold_certificate(u, a)["lambda"], 1)
        for u in (F(0), F(1, 4), F(1, 2)):
            self.assertEqual(exact.threshold_regime(u, F(999999, 1000000)), "eradication")
        for u, a in ((F(1), F(1, 4)), (F(1), F(1, 3)), (F(1), F(1, 2)), (F(4, 5), F(3, 4))):
            for k in (1, 2, 11, 53, 128):
                self.assertTrue(exact.check_eigen_identity(k, u, a))

    def test_positive_mean_bound_and_survival_supersolution(self):
        for u, a in ((F(1), F(1, 4)), (F(1), F(1, 3)), (F(1), F(1, 2)), (F(4, 5), F(3, 4))):
            cert = exact.threshold_certificate(u, a)
            s, lam = cert["s"], cert["lambda"]
            for n in range(8):
                result = exact.pgf_summary((0, 1, 9), n, u, s)
                self.assertLessEqual(result["positive"], 3*10*lam**n)
            for initial in ((1,), (3,), (0, 1, 4), (7, 7)):
                bound = exact.survival_lower_bound(initial, u, a)
                self.assertLessEqual(bound, 1-exact.population_extinction(initial, 3, u, s))
                self.assertEqual(bound > 0, lam > 1)

    def test_monomorphic_fano_collapse_and_two_generation_inverse(self):
        for s in (F(1, 3), F(3, 4), F(8, 9), F(1)):
            for n in range(6):
                expected = exact.pgf_summary((1,), n, 1, s)["pooled_fano"]
                for k in (2, 8, 31):
                    actual = exact.pgf_summary((k, k, k), n, F(3, 5), s)["pooled_fano"]
                    self.assertEqual(actual, expected)
            fano2 = exact.pgf_summary((7,), 2, F(4, 5), s)["pooled_fano"]
            self.assertEqual(fano2, (3+s)/(2*(1+s)))
            self.assertEqual((3-2*fano2)/(2*fano2-1), s)
        self.assertEqual(exact.pgf_summary((5,), 2, 1, F(3, 4))["pooled_fano"], F(15, 14))

    def test_certified_fano_series_rounding_and_long_run_mixture_limit(self):
        for s in (F(1, 2), F(3, 4), F(8, 9)):
            partial = F(1)
            for n in range(6):
                a, d, _ = exact.iterate_jet(s, 1, n)
                partial += d/(1+a)
            lower, upper = exact.fano_limit_interval(s, steps=6, digits=18)
            tail = s**6/(1-s)
            self.assertLessEqual(lower, partial)
            self.assertGreaterEqual(upper-tail, partial)
            lower, upper = exact.fano_limit_interval(s, steps=250, digits=36)
            self.assertLess(upper-lower, F(1, 10**11))
            for initial in ((1,), (0, 1, 9), (2, 4, 7, 11)):
                # An independent precision algorithm checks the certified enclosure;
                # the interval's proof is outward rounding plus the analytic tail.
                decimal_result = exact.decimal_pooled_moments(initial, 300, F(4, 5), s)
                limit_approximation = F(decimal_result["pooled_fano"])
                self.assertLessEqual(lower, limit_approximation)
                self.assertLessEqual(limit_approximation, upper)
                self.assertLess(decimal_result["pooled_mean"], Decimal("1e-40"))
        lower, upper = exact.fano_limit_interval(F(3, 4), steps=130, digits=35)
        self.assertGreater(lower, F("2.6793354938974588"))
        self.assertLess(upper, F("2.6793354938974591"))

    def test_one_step_selection_variance_and_response_identity(self):
        initial = (0, 1, 1, 2, 4, 8)
        for s in (F(1, 4), F(3, 4), F(9, 10)):
            weights = [s**k for k in initial]
            denominator = sum(weights)
            mean = sum(k*w for k, w in zip(initial, weights))/denominator
            variance = sum(k*k*w for k, w in zip(initial, weights))/denominator-mean**2
            result = exact.pgf_summary(initial, 1, F(4, 5), s)
            self.assertEqual(result["pooled_mean"], mean)
            self.assertEqual(result["pooled_variance"], variance+mean/2)
            self.assertEqual(mean/2-result["pooled_variance"], -variance)

    def test_generic_comparison_and_linear_killing_certificates(self):
        u, s = F(9, 10), F(4, 5)
        p = lambda k: u*s**k
        lower = lambda k: p(k)*F(k+2, k+3)
        upper = lambda k: p(k)+(1-p(k))/F(k+3)
        for initial in ((1,), (0, 1, 3), (5,)):
            for n in range(4):
                q = exact.population_extinction(initial, n, u, s)
                self.assertEqual(exact.generic_extinction(initial, n, p), q)
                self.assertGreaterEqual(exact.generic_extinction(initial, n, lower), q)
                self.assertLessEqual(exact.generic_extinction(initial, n, upper), q)
        for k in (1, 2, 3, 5, 17, 64, 151):
            for c in (F(1, 4), F(1, 3), F(1, 2), F(2)):
                self.assertGreaterEqual(1/(1+c*k), (1/(1+c))**k)
                if c >= F(1, 2):
                    self.assertLessEqual(2*(1-F(1, 4)**k)/(1+c*k), 1)

    def test_seeded_finite_population_extinction_and_copy_moments(self):
        initial, n, u, s, trials = (0, 1, 3), 4, F(9, 10), F(4, 5), 12000
        observed = simulation.run_ensemble(initial, n, float(u), float(s), trials=trials)
        theory = exact.pgf_summary(initial, n, u, s)
        mapping = {
            "population": "cells", "zero": "zero", "positive": "positive",
            "copy_sum": "copies", "copy_square_sum": "copy_squares",
        }
        for simulated_name, exact_name in mapping.items():
            raw = exact.population_raw_moments(initial, n, u, s, exact_name)
            se = math.sqrt(float(raw["variance"])/trials)
            self.assert_mc(observed["means"][simulated_name], theory[exact_name], se, exact_name)
        q = exact.population_extinction(initial, n, u, s)
        self.assert_mc(observed["means"]["extinct_positive"], q,
                       math.sqrt(float(q*(1-q))/trials), "extinction")
        for name in ("mean", "variance", "fano"):
            self.assert_mc(observed["pooled"][name], theory["pooled_"+name],
                           observed["pooled"][name+"_standard_error"], "pooled "+name)
        self.assert_mc(observed["pooled"]["zero_fraction"], theory["zero"]/theory["cells"],
                       observed["pooled"]["zero_fraction_standard_error"], "zero fraction")

    def test_seeded_eigenfunction_and_threshold_finite_horizon(self):
        for index, a in enumerate((F(1, 4), F(1, 3), F(1, 2))):
            cert = exact.threshold_certificate(1, a)
            s, lam = cert["s"], cert["lambda"]
            n = 10
            observed = simulation.run_ensemble(
                (1,), n, 1.0, float(s), trials=10000,
                seed=4000+index, keep_zeros=False, weight_base=float(a),
            )
            self.assert_mc(observed["means"]["h_sum"], lam**n,
                           observed["standard_errors"]["h_sum"], "weighted growth")
            self.assert_mc(observed["means"]["positive"],
                           decimal_positive_expectation((1,), n, F(1), s),
                           observed["standard_errors"]["positive"], "positive growth")
            if lam > 1:
                finite_survival = 1-observed["means"]["extinct_positive"]
                bound = float(exact.survival_lower_bound((1,), 1, a))
                self.assertGreaterEqual(
                    finite_survival+6*observed["standard_errors"]["extinct_positive"], bound,
                )

    def test_seeded_fano_collapse_limit_and_zero_dominance(self):
        for k in (1, 8):
            observed = simulation.run_ensemble((k,), 4, 1.0, 0.9,
                                               trials=12000, seed=5000+k)
            theory = exact.pgf_summary((k,), 4, 1, F(9, 10))["pooled_fano"]
            self.assert_mc(observed["pooled"]["fano"], theory,
                           observed["pooled"]["fano_standard_error"], "Fano collapse")
        observed = simulation.run_ensemble((0, 1, 3), 12, 1.0, 0.75,
                                           trials=15000, seed=6011)
        theory = exact.decimal_pooled_moments((0, 1, 3), 12, 1, F(3, 4))
        self.assert_mc(observed["pooled"]["fano"], theory["pooled_fano"],
                       observed["pooled"]["fano_standard_error"], "finite Fano approaching limit")
        self.assertGreater(observed["pooled"]["zero_fraction"], 0.99)
        # A surviving zero founder guarantees a nonzero total denominator here.
        # This is a finite illustration of Theorem 4, not a proof of its limit.
        for seed in range(20):
            state = simulation.simulate_one((0, 1), 18, 1.0, 8/9, random.Random(seed))
            self.assertGreaterEqual(state["zero_count"], 2**18)

    def test_seeded_general_linear_hazard_finite_predictions(self):
        for index, c in enumerate((F(1, 4), F(1, 2))):
            p = lambda k: 1/(1+c*k)
            initial, n = (1, 3), 4
            counts = exact.generic_forward_expectations(initial, n, p)
            theory = exact.summary_from_counts(counts)
            q = exact.generic_extinction(initial, n, p)
            observed = simulation.run_ensemble(
                initial, n, 1.0, 1.0, trials=10000, seed=7000+index,
                division_probability=lambda k: float(p(k)),
            )
            for metric, key in (("population", "cells"), ("positive", "positive"), ("copy_sum", "copies")):
                self.assert_mc(observed["means"][metric], theory[key],
                               observed["standard_errors"][metric], "linear "+metric)
            self.assert_mc(observed["means"]["extinct_positive"], q,
                           math.sqrt(float(q*(1-q))/10000), "linear extinction")

    def test_input_validation_and_uncensored_simulation(self):
        for initial in ((-1,), (True,), (1.5,)):
            with self.assertRaises((ValueError, TypeError)):
                exact.pgf_summary(initial, 1, 1, F(3, 4))
        with self.assertRaises(ValueError):
            exact.extinction_by_generation(1, -1, 1, F(3, 4))
        with self.assertRaises(TypeError):
            exact.pgf_summary((1,), 1, 1, 0.75)
        for s in (0, 1):
            with self.assertRaises(ValueError):
                exact.fano_limit_interval(s)
        with self.assertRaises(ValueError):
            exact.generic_extinction((1,), 1, lambda k: F(3, 2))
        with self.assertRaises(RuntimeError):
            simulation.simulate_one((1, 1), 1, 1.0, 1.0, random.Random(1), max_cells=1)
        for u, s in ((0.0, 0.0), (0.0, 1.0), (1.0, 0.0)):
            state = simulation.simulate_one((1, 2), 1, u, s, random.Random(1))
            self.assertFalse(state["positive_copies"])


if __name__ == "__main__":
    unittest.main()
