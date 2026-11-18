"""Nonclinical ecDNA mathematics: independent exact and stochastic checks.

No patient data or treatment advice. These offline stdlib tests check finite
certificates and falsifiable predictions; the README proves infinite-state
and limiting assertions. Only independent replicate populations set MC errors.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import importlib.util
import math
from pathlib import Path
import random
import unittest


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "ecdna-killing-laws-and-clustering"
F = Fraction


def load_module(name):
    spec = importlib.util.spec_from_file_location("killing_clustering_"+name,
                                                RESEARCH / (name+".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


exact = load_module("verifier")
simulation = load_module("simulation")


@lru_cache(maxsize=None)
def independent_kernel(k, mode="binomial", parameter=F(1, 2)):
    """Sequential Bernoulli/Pólya-urn DP, not the verifier's PMF formula.

    For beta mode parameter is rho, and the urn starts with a balls of each
    color, a=(1-rho)/(2*rho). Rational initial weights are valid urn laws.
    """
    law = {0: F(1)}
    if mode == "beta" and parameter == 1:
        if not k:
            return (F(1),)
        return (F(1, 2),)+(F(0),)*(2*k-1)+(F(1, 2),)
    if mode == "beta" and parameter:
        a = (1-parameter)/(2*parameter)
    else:
        a = None
    for used in range(2*k):
        following = defaultdict(Fraction)
        for left, chance in law.items():
            if a is not None:
                probability = (a+left)/(2*a+used)
            else:
                probability = parameter if mode == "bias" else F(1, 2)
            following[left+1] += chance*probability
            following[left] += chance*(1-probability)
        law = dict(following)
    return tuple(law.get(j, F(0)) for j in range(2*k+1))


def whole_population_law(initial, generations, p, mode="binomial", parameter=F(1, 2)):
    """Enumerate complete populations; retain sisters and empty outcomes."""
    law = {tuple(sorted(initial)): F(1)}
    for _ in range(generations):
        following = defaultdict(Fraction)
        for state, state_chance in law.items():
            descendants = {(): F(1)}
            for k in state:
                offspring = defaultdict(Fraction)
                offspring[()] += 1-p(k)
                for j, chance in enumerate(independent_kernel(k, mode, parameter)):
                    offspring[tuple(sorted((j, 2*k-j)))] += p(k)*chance
                expanded = defaultdict(Fraction)
                for current, current_chance in descendants.items():
                    for children, chance in offspring.items():
                        expanded[tuple(sorted(current+children))] += current_chance*chance
                descendants = dict(expanded)
            for state2, chance in descendants.items():
                following[state2] += state_chance*chance
        law = dict(following)
    return law


def falling(value, order):
    result = 1
    for offset in range(order):
        result *= value-offset
    return result


class EcDNAKillingAndClusteringTests(unittest.TestCase):
    """No network, nonstdlib packages, fitted seeds, or silent simulation caps."""

    def assert_mc(self, observed, expected, se, label):
        self.assertIsNotNone(observed, label)
        self.assertIsNotNone(se, label)
        self.assertLessEqual(abs(observed-float(expected)),
                             6*se+1e-10*max(1, abs(float(expected))), label)

    def test_standalone_exact_verifier_and_rational_thresholds(self):
        report = exact.verify()
        self.assertGreaterEqual(report["exact_checks"], 1500)
        bracket = report["linear_threshold"]
        self.assertEqual(bracket["width"], F(1, 10**12))
        self.assertGreater(bracket["lower"], F(3, 8))
        self.assertLess(bracket["upper"], F(2, 5))
        pair = report["segregation_thresholds"]["same_second_moment"]
        self.assertEqual(pair["c"], F(13, 50))
        self.assertEqual(pair["beta_regime"], "eradication")
        self.assertEqual(pair["bias_regime"], "persistence")

    def test_independent_urn_and_bernoulli_kernels(self):
        for k in (0, 1, 2, 5, 11):
            for rho in (F(0), F(1, 41), F(1, 4), F(1, 3), F(3, 4), F(1)):
                self.assertEqual(exact.beta_binomial_kernel(k, rho),
                                 independent_kernel(k, "beta", rho))
            for bias in (F(0), F(1, 5), F(1, 2), F(3, 4), F(1)):
                self.assertEqual(exact.binomial_kernel(k, bias),
                                 independent_kernel(k, "bias", bias))
        # Direct bit-sequence enumeration is a third construction at small k.
        for k in range(5):
            counts = Counter(bits.bit_count() for bits in range(1 << (2*k)))
            self.assertEqual(exact.binomial_kernel(k),
                             tuple(F(counts[j], 1 << (2*k)) for j in range(2*k+1)))

    def test_complete_population_extinction_and_raw_second_moments(self):
        statistics = {
            "cells": len, "positive": lambda state: sum(k > 0 for k in state),
            "zero": lambda state: state.count(0), "copies": sum,
            "copy_squares": lambda state: sum(k*k for k in state),
        }
        p = lambda k: 1/(1+F(2, 5)*k)
        for mode, parameter, kernel in (
            ("binomial", F(1, 2), exact.binomial_kernel),
            ("beta", F(1, 3), lambda k: exact.beta_binomial_kernel(k, F(1, 3))),
            ("bias", F(3, 4), lambda k: exact.binomial_kernel(k, F(3, 4))),
        ):
            for initial, n in (((), 2), ((0,), 2), ((1,), 2), ((2,), 2), ((0, 1), 2)):
                with self.subTest(mode=mode, initial=initial):
                    law = whole_population_law(initial, n, p, mode, parameter)
                    self.assertEqual(sum(law.values()), 1)
                    q = sum(chance for state, chance in law.items() if not any(state))
                    self.assertEqual(q, exact.population_extinction(initial, n, p, kernel))
                    counts = defaultdict(Fraction)
                    for state, chance in law.items():
                        for k, number in Counter(state).items():
                            counts[k] += chance*number
                    forward = exact.forward_expectations(initial, n, p, kernel)
                    self.assertEqual(dict(counts), forward)
                    for name, statistic in statistics.items():
                        mean = sum(chance*statistic(state) for state, chance in law.items())
                        second = sum(chance*statistic(state)**2 for state, chance in law.items())
                        self.assertEqual(exact.population_raw_moments(initial, n, p, kernel, name),
                                         {"mean": mean, "second": second,
                                          "variance": second-mean*mean})

    def test_certificates_have_independently_checked_head_and_tail_inequalities(self):
        # Independent Pascal rows and direct sister sums check the exact witnesses.
        for c, cap, wanted in ((F("0.378474994554"), 25, "persistence"),
                               (F("0.378474994555"), 48, "eradication")):
            cert = exact.linear_certificate(c, cap)
            self.assertEqual(cert["regime"], wanted)
            weights = cert["weights"]
            strict = False
            for k in range(1, cap+1):
                coefficients = [1]
                for _ in range(2*k):
                    coefficients = [1]+[coefficients[j-1]+coefficients[j]
                                        for j in range(1, len(coefficients))]+[1]
                def weight(j):
                    if not j:
                        return F(0)
                    if j <= cap:
                        return weights[j-1]
                    return F(wanted == "eradication")
                mw = sum(F(number, 1 << (2*k))*(weight(j)+weight(2*k-j))
                         for j, number in enumerate(coefficients))/(1+c*k)
                if wanted == "persistence":
                    self.assertGreaterEqual(mw, weights[k-1])
                    strict |= mw > weights[k-1]
                else:
                    self.assertEqual(mw, weights[k-1])
            if wanted == "persistence":
                self.assertTrue(strict)
            else:
                self.assertGreater(min(weights), 0)
                lhs = 2*max(F(1), max(weights))/(1+c*(cap+1))
                self.assertEqual(lhs, cert["tail_lhs"])
                self.assertLess(lhs, F(1, 4))
                for k in (cap+1, 2*cap, 10**6):
                    self.assertLessEqual(2*max(F(1), max(weights))/(1+c*k), lhs)
        # Independently rebuild the centerpiece's own beta/bias witnesses,
        # using the sequential urn/Bernoulli law instead of either exact PMF.
        c = F(13, 50)
        for mode, parameter, cap, wanted in (
            ("beta", F(1, 4), 14, "eradication"),
            ("bias", F(3, 4), 8, "persistence"),
        ):
            kernel = (lambda k: exact.beta_binomial_kernel(k, parameter)) if mode == "beta" \
                else (lambda k: exact.binomial_kernel(k, parameter))
            cert = exact.linear_certificate(c, cap, kernel)
            self.assertEqual(cert["regime"], wanted)
            weights = cert["weights"]
            strict = False
            def weight(j):
                if not j:
                    return F(0)
                return weights[j-1] if j <= cap else F(wanted == "eradication")
            for k in range(1, cap+1):
                law = independent_kernel(k, mode, parameter)
                mw = sum(chance*(weight(j)+weight(2*k-j))
                         for j, chance in enumerate(law))/(1+c*k)
                if wanted == "eradication":
                    self.assertEqual(mw, weights[k-1])
                else:
                    self.assertGreaterEqual(mw, weights[k-1])
                    strict |= mw > weights[k-1]
            if wanted == "eradication":
                self.assertLess(2*max(F(1), max(weights))/(1+c*(cap+1)), 1)
            else:
                self.assertTrue(strict)

    def test_small_truncation_is_not_a_false_full_model_extinction_certificate(self):
        # The full binomial model persists at its certified lower endpoint,
        # while the one-type truncation is subcritical.
        c = F("0.378474994554")
        small = exact.linear_certificate(c, 1)
        self.assertEqual(small["regime"], "undetermined")
        self.assertGreater(small["tail_lhs"], 1)
        self.assertEqual(exact.linear_certificate(c, 25)["regime"], "persistence")
        # Singular finite eigenproblem, not an arbitrary floating sign choice.
        singular = exact.certificate(lambda k: F(1), 1, F(1))
        self.assertEqual(singular["regime"], "undetermined")

    def test_power_laws_refuge_and_proved_large_exponent_scaling(self):
        for c, cap, wanted in ((F(1, 5), 2, "persistence"),
                               (F(1, 4), 4, "eradication")):
            cert = exact.certificate(lambda k: 1/(1+c*k*k), cap,
                                     1/(1+c*(cap+1)**2))
            self.assertEqual(cert["regime"], wanted)
        for alpha in (16, 24, 40):
            for constant, wanted in ((F(49, 100), "persistence"),
                                      (F(51, 100), "eradication")):
                c = constant/F(2**(alpha//2))
                analytic = exact.two_type_power_witness(c, 2**alpha, 3**alpha)
                general = exact.certificate(lambda k: 1/(1+c*k**alpha), 2,
                                            1/(1+c*3**alpha))
                self.assertEqual(analytic["regime"], wanted)
                self.assertEqual(general["solution"], analytic["solution"])
                self.assertEqual(general["regime"], wanted)
        for c, band in ((F(10), 64), (F(100), 512)):
            p_min = 1/(1+c/(band+1))  # decreasing hazard c/k
            self.assertGreater(2*p_min*exact.refuge_bound(band), 1)
        # Independent evaluation of the finite sine-band inequality. The
        # certificate itself is rational; floats here are an extra cross-check.
        band = 24
        for k in range(band+1, 2*band):
            h = math.sin(math.pi*(k-band)/band)
            mean = sum(float(chance)*math.sin(math.pi*(j-band)/band)
                       for j, chance in enumerate(exact.binomial_kernel(k))
                       if band < j < 2*band)
            self.assertGreater(mean/h, float(exact.refuge_bound(band)))

    def test_saturation_bounds_and_a_whole_family_counterexample(self):
        for h in (F(1, 2), F(1), F(3), F(7)):
            for c in (F(1, 2), F(1), F(3, 2), h+1):
                p1, p2 = (h+1)/(h+1+c), (h+2)/(h+2+2*c)
                determinant = (1-p1)*(1-3*p2/4)-p1*p2/4
                polynomial = 8*c*c+(h+2)*c-(h+1)*(h+2)
                self.assertEqual(determinant,
                                 polynomial/(4*(h+1+c)*(h+2+2*c)))
                if h > 2 and c == 1:
                    self.assertLess(determinant, 0)
                if c >= h+1:
                    for k in (1, 2, 13, 1000):
                        self.assertLessEqual(1/(1+c*k/(h+k)), F(1, 2))
        # Dose-capped hazards: every finite D has plateau below b.
        for dose in (F(0), F(1), F(100), F(10**9)):
            ceiling = dose/(2+dose)
            self.assertLess(ceiling, 1)

    def test_cluster_count_bounds_and_allocation_symmetry(self):
        for rho in (F(0), F(1, 41), F(1, 4), F(1, 2), F(99, 100), F(1)):
            a1 = (1-rho)/2
            for k in (1, 2, 9, 37, 80):
                law = exact.beta_binomial_kernel(k, rho)
                a_k = 1-law[0]-law[-1]
                self.assertLessEqual(a_k, min(k*a1, 1))
                if a1:
                    self.assertLessEqual((1+a_k)/(1+a1*k), 1)
        p = lambda k: 1/(1+F(3, 10)*k)
        for bias in (F(1, 10), F(1, 4), F(1, 2)):
            self.assertEqual(exact.mean_matrix(8, p, lambda k: exact.binomial_kernel(k, bias)),
                             exact.mean_matrix(8, p, lambda k: exact.binomial_kernel(k, 1-bias)))
        # The strict binomial lower bound certifies a whole R interval.
        c_lower = F("0.378474994554")
        for rho in (F(1, 4), F(1, 3), F(1, 2), F(9, 10)):
            self.assertLess((1-rho)/2, c_lower)

    def test_complete_clumping_and_zero_death_boundaries(self):
        kernel = lambda k: exact.beta_binomial_kernel(k, F(1))
        c = F(2, 3)
        for alpha in (-1, 0, 1, 2):
            p = lambda k: F(1) if k == 0 else 1/(1+c*F(k)**alpha)
            for n in range(6):
                survival = F(1)
                for generation in range(n):
                    survival *= p(3*2**generation)
                self.assertEqual(exact.population_extinction((3,), n, p, kernel),
                                 1-survival)
        for kernel in (exact.binomial_kernel,
                       lambda k: exact.beta_binomial_kernel(k, F(1, 2)),
                       lambda k: exact.binomial_kernel(k, F(9, 10))):
            for initial in ((), (0,), (1,), (0, 2)):
                q = exact.population_extinction(initial, 3, lambda k: F(1), kernel)
                self.assertEqual(q, 0 if any(initial) else 1)
        self.assertIsNone(exact.summary_from_counts({})["pooled_fano"])
        self.assertIsNone(exact.summary_from_counts({0: F(3)})["pooled_fano"])

    def test_conditional_scaling_fourth_moment_and_second_moment_obstruction(self):
        for bias in (F(3, 5), F(3, 4), F(9, 10)):
            r = (2*bias-1)**2
            beta = lambda k: exact.beta_binomial_kernel(k, r)
            asymmetric = lambda k: exact.binomial_kernel(k, bias)
            for k in (2, 3, 9, 21):
                for kernel, expected_fourth in ((beta, 3*r*r/(1+2*r)),
                                                 (asymmetric, r*r)):
                    law = kernel(k)
                    pooled_second = sum((j-k)**2*chance for j, chance in enumerate(law))
                    self.assertEqual(pooled_second, (1-r)*k/2+r*k*k)
                    fourth = sum(chance*(falling(j, 4)+falling(2*k-j, 4))/2
                                 for j, chance in enumerate(law))
                    statistic = 16*fourth/F(falling(2*k, 4))-1-6*r
                    self.assertEqual(statistic, expected_fourth)
            for c in (F(0), F(1, 10), F(13, 50), F(2), F(100)):
                p = lambda k: 1/(1+c*k)
                first = exact.summary_from_counts(exact.forward_expectations((1,), 2, p, beta))
                second = exact.summary_from_counts(exact.forward_expectations((1,), 2, p, asymmetric))
                theory = exact.two_generation_prediction(c, r)
                for metric, name in (("pooled_mean", "mean"), ("pooled_fano", "fano")):
                    self.assertEqual(first[metric], second[metric])
                    self.assertEqual(first[metric], theory[name])
                self.assertEqual(first["pooled_zero_fraction"]-second["pooled_zero_fraction"],
                                 theory["beta_minus_bias_zero_fraction"])
                self.assertGreater(theory["beta_minus_bias_zero_fraction"], 0)

    def test_dose_derivative_algebra_and_fixed_input_zero_fraction(self):
        for r in (F(0), F(1, 4), F(1, 3), F(99, 100), F(1)):
            for c in (F(0), F(1, 100), F(13, 50), F(1), F(100)):
                a, d = 2+(3-r)*c, 2+4*c+(1+r)*c*c
                polynomial = (2+2*r*r)+(8-4*r+4*r*r)*c+(8-7*r+2*r*r+r**3)*c*c
                self.assertEqual(a*a-(1-r*r)*d, polynomial)
                prediction = exact.two_generation_prediction(c, r)
                self.assertLess(prediction["mean_derivative"], 0)
                self.assertGreater(prediction["fano_derivative"], 0)
        initial = (0, 1, 1, 2, 5)
        for kernel in (exact.binomial_kernel,
                       lambda k: exact.beta_binomial_kernel(k, F(1, 3)),
                       lambda k: exact.binomial_kernel(k, F(3, 4))):
            previous = F(-1)
            for c in (F(0), F(1, 5), F(2, 5), F(1), F(5)):
                weights = [1/(1+c*k*k) for k in initial]
                z = [(kernel(k)[0]+kernel(k)[-1])/2 for k in initial]
                score = [F(k*k)/(1+c*k*k) for k in initial]
                total = sum(weights)
                average_z = sum(w*x for w, x in zip(weights, z))/total
                average_score = sum(w*x for w, x in zip(weights, score))/total
                covariance = sum(w*x*y for w, x, y in zip(weights, z, score))/total
                covariance -= average_z*average_score
                self.assertLess(covariance, 0)
                result = exact.summary_from_counts(exact.forward_expectations(
                    initial, 1, lambda k: 1/(1+c*k*k), kernel))
                self.assertEqual(average_z, result["pooled_zero_fraction"])
                self.assertGreater(average_z, previous)
                previous = average_z

    def test_exact_falsification_of_finite_time_clustering_order(self):
        p = lambda k: 1/(1+F(3, 8)*k)
        clustered = lambda k: exact.beta_binomial_kernel(k, F(1, 41))
        differences = []
        for k in (1, 5):
            differences.append(exact.population_extinction((k,), 2, p, clustered)
                               -exact.population_extinction((k,), 2, p))
        self.assertEqual(differences[0], F(1200, 381997))
        self.assertEqual(differences[1],
                         -F(28567481394946251, 16939722098833101568))

    def test_seeded_simulation_matches_all_law_families_at_finite_horizons(self):
        cases = (
            ("linear below", (1,), 4, lambda k: 1/(1+F(3, 8)*k), "binomial", F(0), F(1, 2)),
            ("linear above", (1,), 4, lambda k: 1/(1+F(2, 5)*k), "binomial", F(0), F(1, 2)),
            ("beta", (1,), 4, lambda k: 1/(1+F(13, 50)*k), "beta", F(1, 4), F(1, 2)),
            ("bias", (1,), 4, lambda k: 1/(1+F(13, 50)*k), "biased", F(0), F(3, 4)),
            ("quadratic", (1,), 3, lambda k: 1/(1+F(1, 4)*k*k), "binomial", F(0), F(1, 2)),
            ("saturating", (1,), 3, lambda k: 1/(1+F(k, 3+k)), "binomial", F(0), F(1, 2)),
            ("constant critical", (1,), 3, lambda k: F(1, 2) if k else F(1), "binomial", F(0), F(1, 2)),
            ("declining refuge", (65,), 2, lambda k: F(k, k+10) if k else F(1), "binomial", F(0), F(1, 2)),
            ("clumped declining", (2,), 4, lambda k: F(k, k+1) if k else F(1), "beta", F(1), F(1, 2)),
            ("beta declining", (3,), 3, lambda k: F(k, k+1) if k else F(1), "beta", F(1, 3), F(1, 2)),
            ("biased declining", (3,), 3, lambda k: F(k, k+1) if k else F(1), "biased", F(0), F(3, 4)),
            ("steep power lower", (1,), 3, lambda k: 1/(1+F(49, 1600)*k**8), "binomial", F(0), F(1, 2)),
            ("steep power upper", (1,), 3, lambda k: 1/(1+F(51, 1600)*k**8), "binomial", F(0), F(1, 2)),
            ("clumped linear", (1,), 5, lambda k: 1/(1+F(1, 5)*k), "beta", F(1), F(1, 2)),
        )
        trials = 7000
        for index, (name, initial, n, p, mode, rho, bias) in enumerate(cases):
            with self.subTest(case=name):
                if mode == "beta":
                    kernel = lambda k, value=rho: exact.beta_binomial_kernel(k, value)
                elif mode == "biased":
                    kernel = lambda k, value=bias: exact.binomial_kernel(k, value)
                else:
                    kernel = exact.binomial_kernel
                observed = simulation.run_ensemble(
                    initial, n, lambda k: float(p(k)), trials=trials,
                    seed=2026111800+index, mode=mode, rho=float(rho), bias=float(bias))
                counts = exact.forward_expectations(initial, n, p, kernel)
                theory = exact.summary_from_counts(counts)
                for metric in ("cells", "positive", "zero", "copies", "copy_squares"):
                    self.assert_mc(observed["means"][metric], theory[metric],
                                   observed["standard_errors"][metric], name+" "+metric)
                q = exact.population_extinction(initial, n, p, kernel)
                se = math.sqrt(float(q*(1-q))/trials)
                self.assert_mc(observed["means"]["extinct_positive"], q, se, name+" extinction")
                for metric in ("mean", "variance", "fano", "zero_fraction"):
                    self.assert_mc(observed["pooled"][metric], theory["pooled_"+metric],
                                   observed["pooled"][metric+"_standard_error"], name+" "+metric)

    def test_seeded_conditional_variance_and_fourth_moment_predictions(self):
        trials = 10000
        for index, (mode, rho, bias, r) in enumerate((
            ("binomial", 0.0, 0.5, F(0)),
            ("beta", 0.25, 0.5, F(1, 4)),
            ("biased", 0.0, 0.75, F(1, 4)),
        )):
            for k in (2, 9):
                rng = random.Random(2026111900+10*index+k)
                samples = []
                fourths = []
                for _ in range(trials):
                    j, sibling = simulation.sample_split(k, rng, mode, rho, bias)
                    self.assertEqual(j+sibling, 2*k)
                    samples.append((j-k)**2)
                    fourths.append((falling(j, 4)+falling(sibling, 4))/2)
                mean = sum(samples)/trials
                se = math.sqrt(sum((x-mean)**2 for x in samples)/(trials*(trials-1)))
                self.assert_mc(mean, (1-r)*k/2+r*k*k, se, "conditional variance")
                moment4 = sum(fourths)/trials
                se4 = math.sqrt(sum((x-moment4)**2 for x in fourths)/(trials*(trials-1)))
                x4 = 3*r*r/(1+2*r) if mode == "beta" else r*r
                expected = F(falling(2*k, 4))*(1+6*r+x4)/16
                self.assert_mc(moment4, expected, se4, "factorial fourth moment")

    def test_seeded_equal_fano_curves_and_copy_conservation(self):
        for index, (mode, rho, bias) in enumerate((
            ("beta", 0.25, 0.5), ("biased", 0.0, 0.75),
        )):
            c = F(13, 50)
            observed = simulation.run_ensemble(
                (1,), 2, lambda k: 1/(1+float(c)*k), trials=14000,
                seed=2026112000+index, mode=mode, rho=rho, bias=bias)
            theory = exact.two_generation_prediction(c, F(1, 4))
            for metric in ("mean", "fano"):
                self.assert_mc(observed["pooled"][metric], theory[metric],
                               observed["pooled"][metric+"_standard_error"], "same-R "+metric)
            neutral = simulation.run_ensemble(
                (0, 1, 3), 3, lambda k: 1, trials=100, seed=301+index,
                mode=mode, rho=rho, bias=bias)
            self.assertEqual(neutral["means"]["copies"], 32)
            self.assertEqual(neutral["sample_variances"]["copies"], 0)
            self.assertEqual(neutral["means"]["cells"], 24)

    def test_validation_caching_and_uncensored_resource_guards(self):
        exact.binomial_kernel(2, F(1, 2))
        exact.beta_binomial_kernel(2, F(1, 4))
        for call in (lambda: exact.binomial_kernel(2, 0.5),
                     lambda: exact.beta_binomial_kernel(2, 0.25),
                     lambda: exact.binomial_kernel(True),
                     lambda: exact.forward_expectations((-1,), 1, lambda k: F(1)),
                     lambda: exact.forward_expectations((1,), -1, lambda k: F(1)),
                     lambda: exact.population_extinction((1,), 1, lambda k: F(3, 2)),
                     lambda: exact.two_generation_prediction(-1, F(1, 4)),
                     lambda: exact.refuge_bound(1)):
            with self.assertRaises((TypeError, ValueError)):
                call()
        with self.assertRaises(RuntimeError):
            simulation.simulate_one((1, 1), 1, lambda k: 1, random.Random(1), max_cells=1)
        with self.assertRaises(RuntimeError):
            simulation.simulate_one((5,), 1, lambda k: 1, random.Random(1), max_copies=8)
        dead = simulation.run_ensemble((1,), 2, lambda k: 0, trials=20, seed=17)
        self.assertEqual(dead["means"]["extinct_positive"], 1)
        self.assertIsNone(dead["pooled"]["fano"])
        with self.assertRaises(ValueError):
            simulation.run_ensemble((1,), 1, lambda k: float("nan"), trials=2)


if __name__ == "__main__":
    unittest.main()
