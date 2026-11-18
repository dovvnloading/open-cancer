"""Stdlib checks of nonclinical, synthetic immune-escape mathematics.

Exact certificates carry the mathematical assertions. Independent Gillespie
and birth/death simulations are stochastic diagnostics, not biological data.
No network, numpy, patient records, or treatment recommendations are used.
"""

from collections import Counter
from fractions import Fraction as F
import importlib.util
from itertools import product
import math
from pathlib import Path
import statistics
import sys
import unittest


DIRECTORY = (Path(__file__).resolve().parents[1] / "research"
             / "immune-escape-neoantigen-clonality")


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, DIRECTORY / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


V = load("immune_escape_exact_verifier", "verifier.py")
S = load("immune_escape_independent_simulation", "simulation.py")


def simulated_rates(p):
    return S.Parameters(b=float(p.b), d=float(p.d), kappa=float(p.kappa),
                        gain=float(p.u), loss=float(p.ell), hla=float(p.h),
                        b_hla=float(p.bH), d_hla=float(p.dH))


def interval_sum(terms):
    """Independent signed interval arithmetic for linear identities."""
    lower, upper = F(0), F(0)
    for coefficient, bounds in terms:
        a, b = coefficient * bounds[0], coefficient * bounds[1]
        lower += min(a, b)
        upper += max(a, b)
    return lower, upper


def solve_linear_two_rhs(matrix, rhs):
    """Exact Gaussian elimination of a forward population absorption system."""
    n = len(matrix)
    rows = [list(row) + list(value) for row, value in zip(matrix, rhs)]
    for column in range(n):
        pivot = next(i for i in range(column, n) if rows[i][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        divisor = rows[column][column]
        rows[column] = [value / divisor for value in rows[column]]
        for i in range(column + 1, n):
            multiplier = rows[i][column]
            if multiplier:
                rows[i] = [x - multiplier * y
                           for x, y in zip(rows[i], rows[column])]
    answer = [[F(0), F(0)] for _ in range(n)]
    for i in reversed(range(n)):
        for end in (0, 1):
            answer[i][end] = rows[i][n + end] - sum(
                rows[i][j] * answer[j][end] for j in range(i + 1, n))
    return answer


def forward_population_brackets(p, population_cap=3, antigen_cap=3):
    """Independent linear forward-chain enclosure, not a branching PGF.

    A state records ALL simultaneous source cell counts. A surviving H family
    is an absorbing failure, and an extinct H family removes its source cell.
    Population/type overflow is unresolved and assigned both possible fates.
    The two finite Dirichlet problems give lower/upper extinction bounds for
    the original unbounded process; the cap is not silently called extinction.
    """
    states = [s for s in product(range(population_cap + 1),
                                repeat=antigen_cap + 1)
              if 0 < sum(s) <= population_cap]
    index = {s: i for i, s in enumerate(states)}
    matrix = [[F(0)] * len(states) for _ in states]
    rhs = [[F(0), F(0)] for _ in states]

    def transfer(row, state, rate):
        if not rate:
            return
        if state is None or sum(state) > population_cap:
            rhs[row][1] += rate
        elif not any(state):
            rhs[row][0] += rate
            rhs[row][1] += rate
        else:
            matrix[row][index[tuple(state)]] -= rate

    for state, row in index.items():
        for k, count in enumerate(state):
            if not count:
                continue
            death = p.d + p.kappa * k
            matrix[row][row] += count * (p.b + death + p.u + p.ell * k + p.h)
            born = list(state)
            born[k] += 1
            transfer(row, born, count * p.b)
            removed = list(state)
            removed[k] -= 1
            transfer(row, removed, count * (death + p.h * (1 - p.sigma)))
            if k == antigen_cap:
                transfer(row, None, count * p.u)
            else:
                gained = list(removed)
                gained[k + 1] += 1
                transfer(row, gained, count * p.u)
            if k:
                lost = list(removed)
                lost[k - 1] += 1
                transfer(row, lost, count * p.ell * k)
            # The remaining rate count*h*sigma enters successful-H absorption,
            # whose extinction value is zero in both Dirichlet problems.
    answer = solve_linear_two_rhs(matrix, rhs)
    return {state: tuple(answer[i]) for state, i in index.items()}


class ExactModelTests(unittest.TestCase):
    def test_theorem_1_generator_and_critical_surface(self):
        for b, loss, gain in product((F(1, 2), F(1), F(7, 4)),
                                     (F(1, 10), F(1)), (F(1, 2), F(2))):
            p = V.Parameters(b=b, ell=loss, u=gain)
            for k in (0, 1, 2, 11):
                self.assertEqual(*V.eigenfunction_identity(p, k))
                for x in (F(1, 7), F(1, 2), F(1)):
                    g = p.c + (1 - p.c) * x
                    derivative = (p.r + p.u * p.kappa * x / p.a
                                  - k * p.kappa * x / g)
                    generator = (p.b - p.d - p.kappa * k + p.u * (g - 1)
                                 + p.ell * k * (1 / g - 1))
                    self.assertEqual(derivative, generator)
        p = V.Parameters(b=1, d="1/2", u=1, ell="1/4", h=0)
        critical_kappa = (p.b - p.d) * p.ell / (p.u - p.b + p.d)
        self.assertEqual(critical_kappa, F(1, 4))
        self.assertEqual(p.replace(kappa=critical_kappa).r, 0)
        self.assertGreater(p.replace(kappa=critical_kappa / 2).r, 0)
        self.assertLess(p.replace(kappa=critical_kappa * 2).r, 0)
        self.assertEqual(V.solve_extinction(p.replace(kappa=critical_kappa), 3),
                         (F(1), F(1), 0))

    def test_theorem_2_against_independent_forward_population_chain(self):
        for p in (V.Parameters(), V.Parameters(b=2, u="1/4", h="3/5")):
            forward = forward_population_brackets(p)
            backward = V.event_tree_bounds(p, 3, depth=128)
            for k in range(3):
                state = tuple(int(i == k) for i in range(4))
                low, high = forward[state]
                self.assertLessEqual(low, backward[k][0])
                self.assertLessEqual(backward[k][1], high)
                self.assertLess(high - low, F(4, 5))
            state = (0, 1, 1, 0)
            low, high = V.founder_product_bounds(backward, {1: 1, 2: 1})
            self.assertLessEqual(forward[state][0], low)
            self.assertLessEqual(high, forward[state][1])

    def test_theorem_2_refinement_rounding_and_known_quadratic(self):
        p = V.Parameters(b=1, d=1, u=0, ell=0, h=1, bH=1, dH=0)
        # At k=0: q^2-3q+1=0; isolate its smaller root without sqrt/float.
        low, high = V.extinction_bounds(p, 0, depth=64)
        polynomial = lambda x: x * x - 3 * x + 1
        self.assertGreaterEqual(polynomial(low), 0)
        self.assertLessEqual(polynomial(high), 0)
        self.assertLess(high - low, F(1, 10 ** 15))
        p = V.Parameters()
        older = V.event_tree_bounds(p, 4, depth=24)
        finer = V.event_tree_bounds(p, 4, depth=72)
        for k, (low, high) in enumerate(finer):
            self.assertLessEqual(older[k][0], low)
            self.assertLessEqual(high, older[k][1])
            self.assertLessEqual(high - low, V.event_depth_error_bound(p, k, 72))
        low, high, depth = V.solve_extinction(p, 1, F(1, 10 ** 8))
        self.assertLessEqual(high - low, F(1, 10 ** 8))
        self.assertGreater(depth, 0)

    def test_theorem_3_positive_resolvent_and_linear_equations(self):
        p = V.Parameters()
        for z in (F(0), F(1, 10), F(21, 10)):
            bounds = [V.resolvent_bounds(p, k, z, terms=24) for k in range(10)]
            for k in range(1, 8):
                diagonal = p.d + p.kappa * k + p.u + p.ell * k + z - p.b
                low, high = interval_sum(((diagonal, bounds[k]),
                                         (-p.u, bounds[k + 1]),
                                         (-p.ell * k, bounds[k - 1])))
                self.assertLessEqual(low, 1)
                self.assertGreaterEqual(high, 1)
                self.assertLess(bounds[k + 1][1], bounds[k][0])
                self.assertGreater(bounds[k - 1][0] + bounds[k + 1][0],
                                   2 * bounds[k][1])
                self.assertGreater(bounds[k - 1][0] * bounds[k + 1][0],
                                   bounds[k][1] ** 2)
            for order in range(1, 5):
                terms = [(F((-1) ** j * math.comb(order, j)), bounds[j])
                         for j in range(order + 1)]
                self.assertGreater(interval_sum(terms)[0], 0)
            finer = V.resolvent_bounds(p, 3, z, terms=32)
            self.assertLessEqual(bounds[3][0], finer[0])
            self.assertLessEqual(finer[1], bounds[3][1])

    def test_theorem_4_uniform_rare_h_error_and_poisson_limit(self):
        base = V.Parameters()
        t = V.resolvent_bounds(base, 1)
        for n in (20, 40, 80, 160):
            p = base.replace(h=F(1, n))
            low, high = V.extinction_bounds(p, 1, depth=96)
            risk = (1 - high, 1 - low)
            linear = (p.sigma * p.h * t[0], p.sigma * p.h * t[1])
            self.assertLessEqual(risk[1], linear[0])
            self.assertLessEqual(linear[1] - risk[0], V.rare_h_error_bound(p, 1))
            exact = (low ** n, high ** n)
            # Theorem 4's log-product enclosure, checked after exponentiation.
            t0 = V.resolvent_bounds(p, 0)[1]
            lam_low, lam_high = p.sigma * t[0], p.sigma * t[1]
            err_up = p.h * (t0 + p.b * t0 ** 2) * lam_high
            x = p.sigma * p.h * t0
            err_down = n * (p.sigma * p.h * t[1]) ** 2 / (2 * (1 - x))
            predicted_low = V.exp_bounds(-lam_high - err_down)[0]
            predicted_high = V.exp_bounds(-lam_low + err_up)[1]
            self.assertLessEqual(predicted_low, exact[0])
            self.assertLessEqual(exact[1], predicted_high)
        limiting = V.exp_bounds(-base.sigma * t[0])
        self.assertLess(abs(float(exact[0]) - float(limiting[0])), 0.08)

    def test_theorem_5_clonal_advantage_has_a_finite_h_certificate(self):
        p = V.Parameters(h="1/100")
        n = 20
        t = [V.resolvent_bounds(p, k) for k in range(3)]
        a = n * t[1][1]
        b_lower = n * (t[0][0] + t[2][0]) / 2
        b_upper = n * (t[0][1] + t[2][1]) / 2
        constant = t[0][1] + p.b * t[0][1] ** 2
        x = p.sigma * p.h * t[0][1]
        lower_log_a = (-p.sigma * p.h * a
                       - n * (p.sigma * p.h * t[1][1]) ** 2 / (2 * (1 - x)))
        upper_log_b = -p.sigma * p.h * b_lower + p.sigma * p.h ** 2 * constant * b_upper
        self.assertGreater(lower_log_a, upper_log_b)
        intervals = V.event_tree_bounds(p, 2, depth=112)
        uniform = V.founder_product_bounds(intervals, {1: n})
        patchy = V.founder_product_bounds(intervals, {0: n // 2, 2: n // 2})
        self.assertGreater(uniform[0], patchy[1])
        self.assertGreater(F(1, 2) * t[0][0], t[1][1])

    def test_theorem_6_fixed_h_inverse_and_harmonic_limits(self):
        p = V.Parameters()
        intervals = V.event_tree_bounds(p, 192, depth=96)
        coefficient = p.sigma * p.h / p.kappa
        previous_error = F(1)
        for k in (8, 16, 32, 64, 128):
            low, high = intervals[k]
            sandwich = V.finite_h_escape_bounds(p, k)
            self.assertLessEqual(sandwich[0], 1 - high)
            self.assertLessEqual(1 - low, sandwich[1])
            scaled = (k * (1 - high), k * (1 - low))
            error = max(abs(scaled[0] - coefficient), abs(scaled[1] - coefficient))
            self.assertLess(error, previous_error)  # This synthetic sequence only.
            previous_error = error
        self.assertLess(previous_error, F(1, 2000))
        n = 128
        homogeneous = V.founder_product_bounds(intervals, {n: n})
        heterogeneous = V.founder_product_bounds(
            intervals, {n // 2: n // 2, 3 * n // 2: n // 2})
        first_limit = V.exp_bounds(-coefficient)
        second_limit = V.exp_bounds(-coefficient * F(4, 3))
        self.assertLess(abs(float(homogeneous[0] - first_limit[0])), 0.001)
        self.assertLess(abs(float(heterogeneous[0] - second_limit[0])), 0.002)
        self.assertLess(heterogeneous[1], homogeneous[0])

    def test_boundary_cases_and_rejected_overclaims(self):
        p = V.Parameters(h=0)
        self.assertEqual(V.solve_extinction(p, 0), (F(1), F(1), 0))
        for hidden_death in (1, 2):
            self.assertEqual(V.solve_extinction(p.replace(h=1, dH=hidden_death), 1),
                             (F(1), F(1), 0))
        with self.assertRaises(ValueError):
            V.solve_extinction(p.replace(b=3), 1)
        with self.assertRaises(ValueError):
            V.resolvent_bounds(p.replace(b=3), 0)
        with self.assertRaises(TypeError):
            V.Parameters(h=0.1)
        with self.assertRaises(ValueError):
            V.event_tree_bounds(p, -1)
        with self.assertRaises(ValueError):
            S.sample_escape(simulated_rates(p), [-1], 10, 1)
        # Without individual losses the initial count matters to the threshold.
        no_losses = p.replace(b=2, d="1/2", u="1/4", ell=0)
        self.assertGreater(no_losses.r, 0)
        self.assertEqual(V.solve_extinction(no_losses, 2), (F(1), F(1), 0))
        self.assertGreater(V.resolvent_bounds(no_losses, 2)[0], 0)
        # A hostile counterexample: susceptibility convexity cannot be promoted
        # to exact finite-h risk convexity over the whole parameter space.
        hostile = V.Parameters(b=10, d="1/10", u=10, ell="1/10",
                               kappa=1, h="1/10", bH=1, dH="1/4")
        q = V.event_tree_bounds(hostile, 2, depth=176)
        second_difference_lower = -q[0][1] + 2 * q[1][0] - q[2][1]
        second_difference_upper = -q[0][0] + 2 * q[1][1] - q[2][0]
        self.assertGreater(second_difference_lower, F("-0.004580404"))
        self.assertLess(second_difference_upper, F("-0.004580403"))
        # Saturated killing would give this positive k-independent lower bound,
        # so the inverse-load law cannot be asserted for bounded killing.
        saturation = F(5)
        lower_for_all_k = p.sigma * F(1, 10) / (F(1, 10) + p.b + p.d + saturation)
        self.assertGreater(lower_for_all_k, 0)

    def test_prediction_algebra_and_published_aggregate_arithmetic(self):
        # Public published 2x2 summary only; these are not patient-level records.
        self.assertEqual(F(12 * 16, 1 * 2), 96)
        # Under the EXPLICIT ideal observation bridge Q(k)=exp(-A/k), doubling
        # k raises odds by exp(A/(2k))+1 > 2. This is not a clinical theorem.
        for exponent in (F(1, 20), F(1, 2), F(3)):
            elo, ehi = V.exp_bounds(exponent)
            self.assertGreater(elo + 1, 2)
            # Check the polynomial identity exactly; a rounded Decimal value
            # must not be compared with a much narrower rational certificate.
            e = (elo + ehi) / 2
            odds_ratio = (e * e - 1) / (e - 1)
            self.assertEqual(odds_ratio, e + 1)
            self.assertLessEqual(elo + 1, odds_ratio)
            self.assertLessEqual(odds_ratio, ehi + 1)


class IndependentStochasticTests(unittest.TestCase):
    def assert_probability(self, result, bounds):
        # At most 32 bounded Bernoulli comparisons; union error <= 1e-6.
        tolerance = math.sqrt(math.log(64 / 1e-6) / (2 * result["trials"]))
        error = result["hla_false_positive_bound"]
        unresolved = result["unresolved"] / result["trials"]
        observed = result["eradication_frequency"]
        self.assertGreaterEqual(observed + unresolved + error + tolerance,
                                float(bounds[0]))
        self.assertLessEqual(observed - tolerance, float(bounds[1]))

    def assert_empirical_mean(self, values, bounds):
        self.assertTrue(all(value is not None for value in values))
        # A fixed-seed six-standard-error diagnostic, not a distribution-free
        # confidence claim for these unbounded random variables.
        error = 6 * math.sqrt(statistics.variance(values) / len(values)) + 0.002
        observed = statistics.mean(values)
        self.assertGreaterEqual(observed + error, float(bounds[0]))
        self.assertLessEqual(observed - error, float(bounds[1]))

    def test_theorem_1_means_below_at_and_above_threshold(self):
        for i, birth in enumerate((F(1), F(13, 10), F(8, 5))):
            p = V.Parameters(b=birth, h=0)
            result = S.sample_source_at_time(simulated_rates(p), 0, 3.0,
                                             10000, 20261117 + i)
            self.assertEqual(result["censored"], 0)
            self.assert_empirical_mean(result["counts"], V.mean_bounds(p, 0, 3))

    def test_theorems_2_and_4_full_escape_and_rare_h(self):
        for i, (p, loads) in enumerate((
            (V.Parameters(), [0]), (V.Parameters(), [1, 2, 3]),
            (V.Parameters(h="1/100"), [1] * 20),
        )):
            result = S.sample_escape(simulated_rates(p), loads, 8000,
                                     20261217 + i, hla_threshold=24)
            self.assertEqual(result["unresolved"], 0)
            bounds = V.population_extinction_bounds(p, Counter(loads), depth=112)
            self.assert_probability(result, bounds)
            if i == 2:
                risk_low, risk_high = V.rare_h_escape_bounds(p, 1)
                self.assert_probability(result, ((1 - risk_high) ** 20,
                                                 (1 - risk_low) ** 20))

    def test_theorem_3_independent_cell_time_and_resolvent(self):
        p = V.Parameters(h=0)
        for i, k in enumerate((0, 1, 4)):
            result = S.source_lifetime_samples(simulated_rates(p), [k],
                                               10000, 20261317 + i)
            self.assertEqual(result["censored"], 0)
            self.assert_empirical_mean(result["areas"], V.resolvent_bounds(p, k))
            self.assertGreaterEqual(result["largest_type"], k + 1)

    def test_theorem_5_coverage_and_labelled_inheritance(self):
        p = V.Parameters(h="1/100")
        uniform = S.sample_escape(simulated_rates(p), [1] * 20, 8000, 20261417)
        patchy = S.sample_escape(simulated_rates(p), [0] * 10 + [2] * 10,
                                 8000, 20261418)
        self.assert_probability(uniform, V.population_extinction_bounds(p, {1: 20}, 112))
        self.assert_probability(patchy, V.population_extinction_bounds(p, {0: 10, 2: 10}, 112))
        self.assertGreater(uniform["eradication_frequency"], patchy["eradication_frequency"])
        shared = [S.LabelledCell(frozenset({7})) for _ in range(20)]
        private = [S.LabelledCell(frozenset({j})) for j in range(20)]
        first = S.sample_escape(simulated_rates(p), [len(c.antigens) for c in shared],
                                200, 7)
        second = S.sample_escape(simulated_rates(p), [len(c.antigens) for c in private],
                                 200, 7)
        self.assertEqual(first, second)
        parent = S.LabelledCell(frozenset({1, 2}))
        left, right = S.divide_labelled(parent)
        left = S.lose_labelled(S.gain_labelled(left, 3), 1)
        self.assertEqual(left.antigens, frozenset({2, 3}))
        self.assertEqual(right.antigens, frozenset({1, 2}))
        self.assertEqual(S.hide_labelled(right).presented_antigens, frozenset())
        growing = S.grow_labelled(S.Parameters(b=2, d=0.1, kappa=0.01,
                                               gain=1, loss=0.01),
                                  64, 123, initial_antigens=(0,))
        self.assertEqual(growing["status"], "target")
        self.assertEqual(len(growing["cells"]), 64)
        self.assertGreater(growing["next_antigen"], 1)

    def test_theorem_6_large_populations_and_harmonic_penalty(self):
        p = V.Parameters()
        for i, counts in enumerate(({64: 64}, {32: 32, 96: 32})):
            loads = [k for k, number in counts.items() for _ in range(number)]
            result = S.sample_escape(simulated_rates(p), loads, 6000,
                                     20261517 + i, hla_threshold=24)
            self.assertEqual(result["unresolved"], 0)
            bounds = V.population_extinction_bounds(p, counts, depth=96)
            self.assert_probability(result, bounds)
            exponent = p.sigma * p.h / p.kappa * sum(F(n, k) for k, n in counts.items())
            limiting = V.exp_bounds(-exponent)
            self.assertLess(abs(float(bounds[0] - limiting[0])), 0.004)

    def test_hidden_family_and_censoring_are_not_silently_extinction(self):
        p = V.Parameters()
        hidden = S.sample_escape(simulated_rates(p), [], 8000, 20261617,
                                 initial_hla=1, hla_threshold=24)
        self.assert_probability(hidden, (1 - p.sigma, 1 - p.sigma))
        censored = S.sample_escape(simulated_rates(p), [0] * 20, 100, 8,
                                   max_events=1)
        self.assertEqual(censored["unresolved"], 100)
        self.assertEqual(censored["eradicated"], 0)
        area = S.source_lifetime_samples(simulated_rates(p), [0] * 20,
                                         100, 9, max_events=1)
        self.assertEqual(area["censored"], 100)
        self.assertIsNone(area["mean"])


if __name__ == "__main__":
    unittest.main()
