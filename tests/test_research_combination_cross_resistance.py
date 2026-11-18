"""Independent checks of synthetic, nonclinical branching-process mathematics.

The algebra oracle expands the population-state CTMC, not the verifier's
four-dimensional backward ODE. Monte Carlo is a separate implementation and
is a consistency check, not a proof or biological validation.
"""

from collections import defaultdict
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
import importlib.util
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "combination-cross-resistance"
Q = Fraction
ZERO = (0, 0, 0, 0)


def _load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, RESEARCH / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


verifier = _load_module("combination_cross_resistance_verifier", "verifier.py")
simulation = _load_module("combination_cross_resistance_simulation", "simulation.py")
algebra = _load_module("combination_cross_resistance_algebra", "algebra_certificate.py")


@lru_cache(maxsize=None)
def _population_generator_row(model, fraction_a, state):
    """Construct a CTMC row directly from cell-level events, with exact rates.

    No verifier rate builder or ODE function is called. A mutant offspring
    event adds its daughter and retains the parent; a death removes its parent.
    """
    b, d, k, c = model.b, model.d, model.k, model.c
    mu, p, w, nu = model.mu, model.p, model.w, model.nu
    death_rates = (d + k, d + c * k * (1 - fraction_a),
                   d + k * fraction_a, d)
    mutant_channels = (
        ((1, mu * (1 - p) * w), (2, mu * (1 - p) * (1 - w)), (3, mu * p)),
        ((3, nu),),
        ((3, nu),),
        (),
    )
    row = defaultdict(Q)
    total = Q(0)
    for parent, count in enumerate(state):
        if not count:
            continue
        death = count * death_rates[parent]
        if death:
            target = list(state)
            target[parent] -= 1
            row[tuple(target)] += death
            total += death
        for daughter, rate in ((parent, b),) + mutant_channels[parent]:
            if rate:
                target = list(state)
                target[daughter] += 1
                rate *= count
                row[tuple(target)] += rate
                total += rate
    if total:
        row[state] -= total
    return tuple((state, rate) for state, rate in row.items() if rate)


def _apply_generator(distribution, model, fraction_a):
    result = defaultdict(Q)
    for state, weight in distribution.items():
        for target, rate in _population_generator_row(model, fraction_a, state):
            result[target] += weight * rate
    return {state: weight for state, weight in result.items() if weight}


def _ctmc_extinction_coefficients(model, phases, order=4):
    """Coefficients of P(N(T)=0), using chronological row-vector exponentials.

    Phase lengths are fractions of T. Applying a generator at most ``order``
    times reaches every state relevant to that derivative. There is no
    population cutoff approximation. Diagonal terms are retained throughout.
    """
    series = [{(1, 0, 0, 0): Q(1)}] + [{} for _ in range(order)]
    for duration_fraction, fraction_a in phases:
        updated = [defaultdict(Q) for _ in range(order + 1)]
        for degree, distribution in enumerate(series):
            term = distribution
            factor = Q(1)
            for extra in range(order - degree + 1):
                for state, weight in term.items():
                    updated[degree + extra][state] += factor * weight
                if degree + extra < order:
                    term = _apply_generator(term, model, fraction_a)
                    factor *= duration_fraction / (extra + 1)
        series = [dict((state, weight) for state, weight in row.items() if weight)
                  for row in updated]
    return tuple(row.get(ZERO, Q(0)) for row in series)


def _independent_phases(cycles, first="A"):
    first_dose = Q(1) if first == "A" else Q(0)
    return ((Q(1, 2 * cycles), first_dose),
            (Q(1, 2 * cycles), 1 - first_dose)) * cycles


@lru_cache(maxsize=None)
def _linear_sensitivity_generator_row(model, fraction_a, state):
    """Q(c=1+delta) = Q(c=1) + delta * (Q(c=2)-Q(c=1))."""
    base = dict(_population_generator_row(replace(model, c=Q(1)), fraction_a, state))
    shifted = dict(_population_generator_row(replace(model, c=Q(2)), fraction_a, state))
    return tuple((target, base.get(target, Q(0)),
                  shifted.get(target, Q(0)) - base.get(target, Q(0)))
                 for target in base.keys() | shifted.keys())


def _ctmc_sensitivity_coefficients(model, phases, order=6):
    """Exact population CTMC series in T, with coefficients cubic in c-1.

    Only delta powers through three can contribute to the critical curve
    through T**3, because its gap starts at T**3 and delta itself is O(T).
    This uses the independently constructed population generator above, never
    a backward-PGF recurrence or the production symbolic algebra certificate.
    """
    zero_polynomial = (Q(0),) * 4

    def apply(distribution, fraction_a):
        result = {}
        for state, polynomial in distribution.items():
            for target, base, slope in _linear_sensitivity_generator_row(model, fraction_a, state):
                row = result.setdefault(target, [Q(0)] * 4)
                for power in range(4):
                    row[power] += base * polynomial[power]
                    if power:
                        row[power] += slope * polynomial[power - 1]
        return {state: tuple(row) for state, row in result.items() if any(row)}

    series = [{(1, 0, 0, 0): (Q(1), Q(0), Q(0), Q(0))}] + [{} for _ in range(order)]
    for duration_fraction, fraction_a in phases:
        updated = [{} for _ in range(order + 1)]
        for degree, distribution in enumerate(series):
            term = distribution
            factor = Q(1)
            for extra in range(order - degree + 1):
                for state, polynomial in term.items():
                    row = updated[degree + extra].setdefault(state, [Q(0)] * 4)
                    for power in range(4):
                        row[power] += factor * polynomial[power]
                if degree + extra < order:
                    term = apply(term, fraction_a)
                    factor *= duration_fraction / (extra + 1)
        series = [{state: tuple(row) for state, row in distribution.items() if any(row)}
                  for distribution in updated]
    return tuple(distribution.get(ZERO, zero_polynomial) for distribution in series)


def _ctmc_critical_curve_coefficients(model):
    """Solve the CTMC gap coefficients for c=1+alpha*T+beta*T²+gamma*T³."""
    sequence = _ctmc_sensitivity_coefficients(model, _independent_phases(1))
    combination = _ctmc_sensitivity_coefficients(model, ((Q(1), Q(1, 2)),))
    gap = tuple(tuple(x - y for x, y in zip(left, right))
                for left, right in zip(sequence, combination))
    if any(any(polynomial) for polynomial in gap[:3]) or gap[3][0]:
        raise AssertionError("The independently computed gap does not start as expected")
    c3, c4, c5, c6 = gap[3:]
    alpha = -c4[0] / c3[1]
    beta = -(c3[2] * alpha**2 + c4[1] * alpha + c5[0]) / c3[1]
    gamma = -(2 * c3[2] * alpha * beta + c3[3] * alpha**3 + c4[1] * beta
              + c4[2] * alpha**2 + c5[1] * alpha + c6[0]) / c3[1]
    return alpha, beta, gamma


@lru_cache(maxsize=None)
def _bounds(model, schedule):
    return verifier.extinction_bounds(model, schedule, order=18, bits=72)


def _decimal(value):
    return Decimal(value.numerator) / Decimal(value.denominator)


class CombinationCrossResistanceTests(unittest.TestCase):
    def assert_tight_probability(self, interval):
        self.assertGreaterEqual(interval.lower, 0)
        self.assertLessEqual(interval.upper, 1)
        self.assertLessEqual(interval.width, Q(1, 10**8))

    def test_independent_population_generator_matches_local_coefficients(self):
        models = (
            verifier.Model(),
            verifier.Model(b=Q(2, 7), d=Q(1, 3), k=Q(4, 5), c=Q(3, 2),
                           mu=Q(2, 9), p=Q(2, 5), w=Q(2, 5), nu=Q(1, 11)),
            verifier.Model(b=Q(7, 9), d=Q(0), k=Q(3, 4), c=Q(7, 3),
                           mu=Q(3, 8), p=Q(0), w=Q(1, 4), nu=Q(2, 13)),
            verifier.Model(b=Q(1, 5), d=Q(2, 7), k=Q(5, 6), c=Q(5, 2),
                           mu=Q(4, 9), p=Q(7, 10), w=Q(3, 4), nu=Q(3, 11)),
        )
        for model in models:
            reference = _ctmc_extinction_coefficients(model, ((Q(1), Q(1, 2)),))
            for cycles in (1, 2, 3):
                for first in ("A", "B"):
                    with self.subTest(model=model, cycles=cycles, first=first):
                        actual = _ctmc_extinction_coefficients(
                            model, _independent_phases(cycles, first))
                        gap = tuple(a - b for a, b in zip(actual, reference))
                        self.assertEqual(gap[:3], (0, 0, 0))
                        self.assertEqual(gap[3:], verifier.local_coefficients(model, cycles, first))

    def test_boundary_has_strict_negative_fourth_order_term(self):
        for w in (Q(1, 2), Q(2, 5), Q(1, 3), Q(1, 4), Q(1, 5)):
            c = (1 - w) / w
            model = verifier.Model(c=c, w=w, p=Q(2, 5), d=Q(0))
            a = model.mu * (1 - model.p) * w
            reference = _ctmc_extinction_coefficients(model, ((Q(1), Q(1, 2)),))
            for cycles in (1, 2, 3):
                with self.subTest(w=w, cycles=cycles):
                    actual = _ctmc_extinction_coefficients(model, _independent_phases(cycles))
                    gap = tuple(x - y for x, y in zip(actual, reference))
                    expected = -a * c * model.k**2 * (model.d + model.k) * (
                        (c - 1) * (8 * cycles**2 - 1) + 2 * cycles * (c + 1)
                    ) / (384 * cycles**3)
                    self.assertEqual(gap[:4], (0, 0, 0, 0))
                    self.assertEqual(gap[4], expected)
                    self.assertLess(expected, 0)
                    reverse = _ctmc_extinction_coefficients(model, _independent_phases(cycles, "B"))
                    reverse_gap = tuple(x - y for x, y in zip(reverse, reference))
                    reverse_expected = a * c * model.k**2 * (model.d + model.k) * (
                        (c - 1) * (8 * cycles**2 - 1) - 2 * cycles * (c + 1)
                    ) / (384 * cycles**3)
                    self.assertEqual(reverse_gap[:4], (0, 0, 0, 0))
                    self.assertEqual(reverse_gap[4], reverse_expected)
                    self.assertEqual(verifier.local_coefficients(model, cycles, "B"),
                                     (0, reverse_expected))

    def test_boundary_existence_can_require_more_than_one_cycle(self):
        # The boundary can defeat BOTH orders at a specified m. Existence
        # over some m is different from existence for every fixed m.
        for c, w, cycles, sign in ((Q(3, 2), Q(2, 5), 1, -1),
                                   (Q(3, 2), Q(2, 5), 2, 1),
                                   (Q(4), Q(1, 5), 1, 1),
                                   (Q(1), Q(1, 2), 3, -1)):
            model = verifier.Model(c=c, w=w)
            reference = _ctmc_extinction_coefficients(model, ((Q(1), Q(1, 2)),))
            reverse = _ctmc_extinction_coefficients(model, _independent_phases(cycles, "B"))
            gap = reverse[4] - reference[4]
            self.assertEqual(reverse[:4], reference[:4])
            self.assertEqual((gap > 0) - (gap < 0), sign)
            self.assertLess(verifier.local_coefficients(model, cycles, "A")[1], 0)

    def test_threshold_sign_and_cross_resistance_scaling(self):
        for cycles in (1, 2, 3):
            for c, expected_sign in ((Q(3, 2), -1), (Q(2), 0), (Q(5, 2), 1)):
                model = verifier.Model(c=c, w=Q(1, 3), p=Q(0))
                cubic, _ = verifier.local_coefficients(model, cycles)
                self.assertEqual((cubic > 0) - (cubic < 0), expected_sign)
                attenuated, _ = verifier.local_coefficients(replace(model, p=Q(4, 5)), cycles)
                self.assertEqual(attenuated, cubic / 5)
                reverse, _ = verifier.local_coefficients(model, cycles, "B")
                self.assertEqual(reverse, -cubic)

    def test_schedules_preserve_total_time_and_each_drug_exposure(self):
        duration = Q(7, 5)
        for cycles in (1, 2, 3):
            for first in ("A", "B"):
                schedule = verifier.alternating(duration, cycles, first)
                self.assertEqual(sum(t for t, _ in schedule), duration)
                self.assertEqual(sum(t * u for t, u in schedule), duration / 2)
                self.assertEqual(sum(t * (1 - u) for t, u in schedule), duration / 2)
                self.assertEqual(verifier.exposure(schedule), verifier.exposure(verifier.combination(duration)))
                self.assertEqual(schedule, tuple((duration * t, u) for t, u in _independent_phases(cycles, first)))
        with self.assertRaises(ValueError):
            verifier.comparison_bounds(verifier.Model(), ((Q(1), Q(1)),), verifier.combination(Q(1)))

    def test_rational_certificates_show_improvement_and_longer_course_reversal(self):
        model = verifier.Model()
        for duration, sign in ((Q(1, 2), 1), (Q(2), -1)):
            reference = _bounds(model, verifier.combination(duration))[0]
            self.assert_tight_probability(reference)
            for cycles in (1, 2):
                with self.subTest(duration=duration, cycles=cycles):
                    candidate = _bounds(model, verifier.alternating(duration, cycles))[0]
                    self.assert_tight_probability(candidate)
                    gap = candidate - reference
                    self.assertGreater(gap.lower if sign > 0 else -gap.upper, 0)
        # With the same parameters and course, one pre-existing A-resistant
        # founder reverses the pure-S result. The initial-state assumption matters.
        resistant = _bounds(model, verifier.alternating(Q(1, 2)))[1]
        resistant_reference = _bounds(model, verifier.combination(Q(1, 2)))[1]
        self.assertLess((resistant - resistant_reference).upper, 0)

    def test_balanced_boundary_and_population_product(self):
        model = verifier.Model(c=Q(1))
        duration = Q(1, 2)
        candidate = _bounds(model, verifier.alternating(duration))
        reference = _bounds(model, verifier.combination(duration))
        self.assertLess((candidate[0] - reference[0]).upper, 0)
        for founders in (1, 2, 7):
            left = verifier.population_bounds(candidate, (founders, 0, 0, 0))
            right = verifier.population_bounds(reference, (founders, 0, 0, 0))
            self.assertEqual(left.lower, candidate[0].lower**founders)
            self.assertEqual(left.upper, candidate[0].upper**founders)
            self.assertLess((left - right).upper, 0)
        self.assertEqual(verifier.population_bounds(candidate, ZERO), verifier.Interval(Q(1), Q(1)))

    def test_certified_finite_course_cross_resistance_reversal(self):
        # A constructed near-boundary witness, with all mutation channels
        # positive. The tiny effect requires narrower exact enclosures.
        # These signs are not claimed to be detectable by the Monte Carlo test.
        duration = Q(1, 10)
        collateral = Q(5084894333, 5000000000)
        for p, sign in ((Q(1, 10), -1), (Q(2, 5), 1)):
            model = verifier.Model(c=collateral, p=p)
            gap = verifier.comparison_bounds(
                model, verifier.alternating(duration), verifier.combination(duration),
                order=28, bits=112)
            self.assertLess(gap.width, Q(1, 10**18))
            self.assertGreater(gap.lower if sign > 0 else -gap.upper, 0)

    def test_degenerate_exact_ties_have_overlapping_enclosures(self):
        for model in (verifier.Model(p=Q(1)), verifier.Model(mu=Q(0)), verifier.Model(k=Q(0))):
            with self.subTest(model=model):
                self.assertEqual(verifier.local_coefficients(model), (0, 0))
                reference = _bounds(model, verifier.combination(Q(1, 3)))[0]
                candidate = _bounds(model, verifier.alternating(Q(1, 3)))[0]
                self.assert_tight_probability(candidate)
                self.assert_tight_probability(reference)
                # Overlap is only a numerical consistency check; the model's
                # reachable rates prove equality, not an interval containing 0.
                gap = candidate - reference
                self.assertLessEqual(gap.lower, 0)
                self.assertGreaterEqual(gap.upper, 0)
        model = verifier.Model(p=Q(1))
        for fraction_a in (Q(0), Q(1, 2), Q(1)):
            distribution = {(1, 0, 0, 0): Q(1)}
            for _ in range(4):
                distribution = _apply_generator(distribution, model, fraction_a)
                self.assertTrue(all(state[1] == state[2] == 0 for state in distribution))

    def test_doubly_resistant_component_matches_scalar_birth_death_solution(self):
        duration = Q(2, 3)
        for birth, death in ((Q(2, 5), Q(1, 10)), (Q(1, 10), Q(2, 5)),
                             (Q(1, 3), Q(1, 3)), (Q(0), Q(1, 3))):
            model = verifier.Model(b=birth, d=death, k=Q(0), mu=Q(0), nu=Q(0))
            bounds = _bounds(model, verifier.combination(duration))
            with localcontext() as context:
                context.prec = 60
                b, d, t = map(_decimal, (birth, death, duration))
                if birth == death:
                    exact = d * t / (1 + d * t)
                else:
                    decay = (-(b - d) * t).exp()
                    exact = d * (1 - decay) / (b - d * decay)
                for founders in (1, 3):
                    interval = verifier.population_bounds(bounds, (0, 0, 0, founders))
                    self.assert_tight_probability(interval)
                    self.assertLessEqual(_decimal(interval.lower), exact**founders)
                    self.assertGreaterEqual(_decimal(interval.upper), exact**founders)

    def test_enclosure_refinement_preserves_exact_critical_birth_death_probability(self):
        # With birth=death=1 and no mutation or drug killing, each founder's
        # exact extinction probability at T=1/2 is T/(1+T)=1/3.
        model = verifier.Model(b=Q(1), d=Q(1), k=Q(0), mu=Q(0), nu=Q(0))
        schedule = verifier.combination(Q(1, 2))
        coarse = verifier.extinction_bounds(model, schedule, order=4, bits=16)
        refined = verifier.extinction_bounds(model, schedule, order=18, bits=72)
        truth = Q(1, 3)
        for earlier, later in zip(coarse, refined):
            self.assertLessEqual(earlier.lower, truth)
            self.assertGreaterEqual(earlier.upper, truth)
            self.assertLessEqual(later.lower, truth)
            self.assertGreaterEqual(later.upper, truth)
            self.assertLess(later.width, earlier.width / 1000)
        # The theorem guarantees convergence; it does not require nested bounds.

    def test_independent_sixth_order_cross_resistance_effect_on_critical_curve(self):
        models = (
            verifier.Model(),
            verifier.Model(b=Q(1, 20), d=Q(2), k=Q(1), c=Q(1),
                           mu=Q(6), nu=Q(1, 50)),
            verifier.Model(b=Q(2), d=Q(3), k=Q(5), c=Q(1), mu=Q(7), nu=Q(11)),
        )
        for model in models:
            cubic_coefficients = []
            for p in (Q(0), Q(1, 2)):
                with self.subTest(model=model, p=p):
                    alpha, beta, gamma = _ctmc_critical_curve_coefficients(replace(model, p=p))
                    expected_beta = model.k * (-115 * model.b + 2 * model.d + 55 * model.k
                                                - 2 * model.mu + 61 * model.nu) / 720
                    self.assertEqual(alpha, model.k / 6)
                    self.assertEqual(beta, expected_beta)
                    cubic_coefficients.append(gamma)
            self.assertEqual(cubic_coefficients[1] - cubic_coefficients[0],
                             -model.mu * model.k**2 * Q(1, 2) / 1440)

    def test_generic_rational_polynomial_proof_certificate(self):
        # verify() actually computes Lie coefficients and requires every
        # parameter-wide polynomial difference to have zero coefficients.
        # The independent forward-CTMC test above uses a different algorithm.
        certificate = algebra.verify()
        self.assertEqual(certificate["identities"], 11)
        self.assertEqual(certificate["maximum_time_degree"], 6)
        self.assertEqual(certificate["threshold_p_slope"], "-mu*k**2/1440")

    def test_zero_horizon_and_invalid_parameters(self):
        zero_bounds = verifier.extinction_bounds(verifier.Model(), verifier.combination(Q(0)))
        self.assertEqual(zero_bounds, (verifier.Interval(Q(0), Q(0)),) * 4)
        for kwargs in ({"mu": Q(-1)}, {"p": Q(2)}, {"c": Q(1, 2)}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                verifier.Model(**kwargs)
        with self.assertRaises(TypeError):
            verifier.Model(p=0.1)
        with self.assertRaises(ValueError):
            verifier.alternating(Q(1), cycles=0)

    def test_independent_fixed_seed_gillespie_agrees_with_certificates(self):
        model = verifier.Model(b=Q(1, 20), d=Q(2), k=Q(1), c=Q(8),
                               mu=Q(6), p=Q(1, 10), w=Q(1, 2), nu=Q(1, 50))
        parameters = {name: float(getattr(model, name))
                      for name in ("b", "d", "k", "c", "mu", "p", "w", "nu")}
        duration = Q(3, 10)
        schedules = (verifier.combination(duration), verifier.alternating(duration),
                     verifier.alternating(duration, cycles=2))
        trials = 150_000
        estimates = []
        certificates = []
        for offset, schedule in enumerate(schedules):
            certificate = _bounds(model, schedule)[0]
            self.assert_tight_probability(certificate)
            result = simulation.simulate(parameters, [(float(t), float(u)) for t, u in schedule],
                                         trials=trials, seed=20261117 + offset)
            self.assertEqual(result.trials, trials)
            truth = float(certificate.midpoint)
            standard_error = math.sqrt(truth * (1 - truth) / trials)
            self.assertLessEqual(abs(result.probability - truth),
                                 6 * standard_error + float(certificate.width))
            estimates.append(result.probability)
            certificates.append(certificate)
        for index in (1, 2):
            self.assertGreater((certificates[index] - certificates[0]).lower, 0)
            self.assertGreater(estimates[index], estimates[0])
        # Neither the fixed-seed sign nor a six-SE agreement is a theorem.
        # Strict inequalities are supplied independently by rational enclosures.


if __name__ == "__main__":
    unittest.main()
