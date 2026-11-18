"""Independent checks of nonclinical metastatic-seeding mathematics.

Only synthetic laws and published aggregate counts are used. Fraction outcome
enumeration, a separately implemented moment linear program, high-precision
Decimal probabilities and fixed-seed event simulations check the implementations.
These checks do not establish biological validity, clinical utility, or priority.
Standard library only; no existing repository file is modified by this test.
"""

from decimal import Decimal, localcontext
from fractions import Fraction
import importlib.util
from itertools import combinations
import math
from pathlib import Path
import random
import sys
import unittest


RESEARCH = (Path(__file__).resolve().parents[1] / "research"
            / "metastatic-seeding-timing")
F = Fraction
D = Decimal


def load_module(name):
    spec = importlib.util.spec_from_file_location(
        "metastatic_seeding_timing_" + name, RESEARCH / (name + ".py")
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


exact = load_module("verifier")
simulation = load_module("simulation")


def dec(value):
    value = F(value)
    return D(value.numerator) / D(value.denominator)


def poisson_cdf_decimal(mean, k):
    """Independent direct finite polynomial and Decimal exponential."""
    if k < 0:
        return D(0)
    mean = dec(mean)
    if mean == 0:
        return D(1)
    return (-mean).exp() * sum(
        (mean**j / D(math.factorial(j)) for j in range(k + 1)), D(0)
    )


def polynomial_value(coefficients, x):
    value = F(0)
    for coefficient in reversed(coefficients):
        value = value * x + coefficient
    return value


def multiply_polynomials(left, right):
    result = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def trim_polynomial(coefficients):
    result = list(coefficients)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def categorical_coverage(probabilities, draws):
    """State recursion on observed label masks; no inclusion-exclusion."""
    probabilities = tuple(map(F, probabilities))
    missing = 1 - sum(probabilities)
    states = {0: F(1)}
    for _ in range(draws):
        updated = {}
        for mask, mass in states.items():
            updated[mask] = updated.get(mask, F(0)) + mass * missing
            for index, probability in enumerate(probabilities):
                target = mask | (1 << index)
                updated[target] = updated.get(target, F(0)) + mass * probability
        states = updated
    return states.get((1 << len(probabilities)) - 1, F(0))


def finite_pool_by_subsets(size, invisible, regions, sensitivity):
    """Enumerate distinct region subsets and condition on their visible count."""
    total = F(0)
    denominator = math.comb(size, regions)
    for selected in combinations(range(size), regions):
        visible = sum(index >= invisible for index in selected)
        total += (1 - sensitivity)**visible
    return total / denominator


def solve_square(matrix, rhs):
    """Independent Fraction Gaussian elimination with row pivoting."""
    rows = [list(map(F, row)) + [F(value)] for row, value in zip(matrix, rhs)]
    width = len(rows)
    for column in range(width):
        pivot = next((i for i in range(column, width)
                      if rows[i][column] != 0), None)
        if pivot is None:
            return None
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [x / scale for x in rows[column]]
        for i in range(width):
            if i != column:
                scale = rows[i][column]
                rows[i] = [x - scale * y for x, y in zip(rows[i], rows[column])]
    return tuple(row[-1] for row in rows)


def finite_pool_lp_oracle(size, first, second, regions, sensitivity):
    """Enumerate basic feasible measures, using independently solved systems."""
    moments = lambda u: (F(1), F(u, size), F(u * (u - 1), size * (size - 1)))
    target = (F(1), F(first), F(second))
    feasible = []
    for support_size in (1, 2, 3):
        for support in combinations(range(size + 1), support_size):
            columns = tuple(moments(u) for u in support)
            # Use every square selection of moment rows, then check all rows.
            for selected_rows in combinations(range(3), support_size):
                matrix = [[column[row] for column in columns] for row in selected_rows]
                weights = solve_square(matrix, [target[row] for row in selected_rows])
                if weights is None or min(weights) < 0:
                    continue
                if any(sum(weights[j] * columns[j][i] for j in range(support_size))
                       != target[i] for i in range(3)):
                    continue
                objective = sum(
                    (weight * finite_pool_by_subsets(size, u, regions, sensitivity)
                     for u, weight in zip(support, weights)), F(0)
                )
                feasible.append(objective)
                break
    if not feasible:
        raise ValueError("infeasible moment vector")
    return min(feasible), max(feasible)


class MetastaticSeedingTimingTests(unittest.TestCase):
    """Boundary cases and independent numerical checks of every theorem."""

    def assert_decimal_enclosed(self, interval, value):
        self.assertIsInstance(interval.lower, F)
        self.assertIsInstance(interval.upper, F)
        self.assertLessEqual(interval.lower, interval.upper)
        # Decimal evaluates the independent oracle at 80-digit precision;
        # this allowance is far smaller than the verifier's rational widths.
        tolerance = D("1e-65")
        self.assertLessEqual(dec(interval.lower) - tolerance, value)
        self.assertLessEqual(value, dec(interval.upper) + tolerance)

    def test_exact_exponential_and_poisson_against_decimal_oracle(self):
        with localcontext() as context:
            context.prec = 80
            for mean in (F(0), F(1, 17), F(1), F(7, 2), F(25)):
                interval = exact.exp_neg_bounds(mean)
                self.assert_decimal_enclosed(interval, (-dec(mean)).exp())
                self.assertLessEqual(interval.upper - interval.lower, F(1, 10**30))
                for k in (-1, 0, 1, 5, 20, 40):
                    self.assert_decimal_enclosed(
                        exact.poisson_cdf_bounds(mean, k), poisson_cdf_decimal(mean, k)
                    )
                    self.assert_decimal_enclosed(
                        exact.poisson_tail_bounds(mean, k),
                        1 - poisson_cdf_decimal(mean, k - 1)
                    )

    def test_t1_survival_cox_law_and_shared_edge_counts(self):
        with localcontext() as context:
            context.prec = 80
            for birth, death in ((2, 1), (1, 0), (1, 1), (1, 2), (0, 3), (0, 0)):
                for age in (F(0), F(1, 7), F(2), F(30)):
                    if birth == death:
                        expected = 1 / (1 + D(birth) * dec(age))
                    else:
                        rate = D(birth - death)
                        expected = rate / (birth - death * (-rate * dec(age)).exp())
                    self.assert_decimal_enclosed(
                        exact.bd_survival_bounds(birth, death, age), expected
                    )
            integrated = exact.piecewise_seed_mean(
                F(2, 7), ((F(1, 3), 2), (F(1, 2), 5), (F(1, 6), 0))
            )
            self.assertEqual(integrated, F(19, 21))
            for z in (F(0), F(1, 3), F(1)):
                self.assert_decimal_enclosed(
                    exact.cox_count_pgf_bounds(integrated, z),
                    (-dec(integrated * (1 - z))).exp()
                )
            edges = ((F(1), ("P", "M")), (F(2), ("P",)),
                     (F(3), ("M",)), (F(1, 2), ("M",)))
            means = exact.mutation_pattern_means(edges, F(2), exposure=3)
            self.assertEqual(means, {frozenset(("P", "M")): F(6),
                                     frozenset(("P",)): F(12),
                                     frozenset(("M",)): F(21)})
            z_p, z_m = F(2, 3), F(1, 4)
            # The shared edge has z_P*z_M, not two independent leaf factors.
            exponent = 6 * (z_p * z_m - 1) + 12 * (z_p - 1) + 21 * (z_m - 1)
            self.assert_decimal_enclosed(
                exact.mutation_pattern_pgf_bounds(means, {"P": z_p, "M": z_m}),
                dec(exponent).exp()
            )

    def test_t2_exact_confounding_and_direction_of_least_favorable_rates(self):
        for c, h in ((F(1), F(1)), (F(5), F(1)), (F(4), F(4))):
            endpoint = c + h
            for lag_bound in (F(0), h, 2 * h, 10 * h):
                result = exact.anchor_clock_rates(F(3, 2), endpoint, c, h, lag_bound)
                if lag_bound >= 2 * h:
                    self.assertEqual(result["status"], "confounded")
                    early_s, early_lag = result["early_pair"]
                    late_s, late_lag = result["late_pair"]
                    self.assertLessEqual(early_s, c - h)
                    self.assertGreaterEqual(late_s, c + h)
                    self.assertTrue(0 <= early_lag <= min(lag_bound, early_s))
                    self.assertTrue(0 <= late_lag <= min(lag_bound, late_s))
                    self.assertEqual(early_s - early_lag, late_s - late_lag)
                else:
                    self.assertEqual(result["status"], "separated")
                    self.assertGreater(result["lambda_early"], result["lambda_late"])
                    for seed in (F(0), (c - h) / 2, c - h):
                        for lag in (F(0), min(seed, lag_bound)):
                            self.assertGreaterEqual(F(3, 2) * (endpoint - seed + lag),
                                                    result["lambda_early"])
                    for seed in (c + h, endpoint):
                        for lag in (F(0), min(seed, lag_bound)):
                            self.assertLessEqual(F(3, 2) * (endpoint - seed + lag),
                                                 result["lambda_late"])

    def test_t2_minimum_deterministic_exposure_against_exhaustive_decimal_search(self):
        with localcontext() as context:
            context.prec = 80
            for early, late, alpha in ((F(4), F(1), F(1, 20)),
                                       (F(4), F(2), F(1, 20)),
                                       (F(1), F(0), F(1, 20))):
                result = exact.minimal_threshold_exposure(early, late, alpha, 20)
                self.assertEqual(result["status"], "attained")
                found = None
                for exposure in range(1, 21):
                    for threshold in range(1, 151):
                        early_error = poisson_cdf_decimal(exposure * early, threshold - 1)
                        late_error = 1 - poisson_cdf_decimal(exposure * late, threshold - 1)
                        if max(early_error, late_error) <= dec(alpha):
                            found = (exposure, threshold)
                            break
                    if found is not None:
                        break
                self.assertEqual((result["exposure"], result["threshold"]), found)
                self.assertEqual(len(result["rejected_exposures"]), found[0] - 1)
                for certificate in result["rejected_exposures"]:
                    self.assertFalse(certificate["feasible"])
                    self.assertGreater(certificate["early"].lower, alpha)
                self.assertLessEqual(result["early"].upper, alpha)
                self.assertLessEqual(result["late"].upper, alpha)
                self.assertGreater(result["previous_threshold_late"].lower, alpha)
                delta = (dec(early).sqrt() - dec(late).sqrt())**2
                sufficient = 2 * (1 / dec(alpha)).ln() / delta
                necessary = (1 / (4 * dec(alpha) * (1 - dec(alpha)))).ln() / delta
                self.assertLessEqual(necessary, D(found[0]))
                self.assertLessEqual(found[0], math.ceil(sufficient))
            self.assertEqual(exact.minimal_threshold_exposure(1, 1, F(1, 10))["status"],
                             "confounded")
            self.assertEqual(exact.minimal_threshold_exposure(4, 2, F(1, 20), 2)["status"],
                             "search_exhausted")

    def test_t3_sharp_moments_and_all_degenerate_boundaries(self):
        laws = (
            ((F(0), F(1)),),
            ((F(1), F(1)),),
            ((F(2, 5), F(1)),),
            ((F(0), F(1, 3)), (F(1), F(2, 3))),
            ((F(0), F(1, 5)), (F(1, 3), F(2, 5)), (F(1), F(2, 5))),
            ((F(1, 7), F(2, 7)), (F(3, 5), F(3, 7)), (F(4, 5), F(2, 7))),
        )
        for law in laws:
            first = sum(f * weight for f, weight in law)
            second = sum(f * f * weight for f, weight in law)
            for regions in (0, 1, 2, 3, 7, 25):
                for sensitivity in (F(0), F(1, 3), F(1)):
                    result = exact.sharp_miss_bounds(first, second, regions, sensitivity)
                    observed = sum(weight * (1 - sensitivity * f)**regions
                                   for f, weight in law)
                    self.assertLessEqual(result["lower"], observed)
                    self.assertLessEqual(observed, result["upper"])
                    for side in ("lower", "upper"):
                        atoms = result[side + "_atoms"]
                        self.assertEqual(sum(weight for _, weight in atoms), 1)
                        self.assertTrue(all(0 <= f <= 1 and weight >= 0 for f, weight in atoms))
                        self.assertEqual(sum(f * weight for f, weight in atoms), first)
                        self.assertEqual(sum(f * f * weight for f, weight in atoms), second)
                        self.assertEqual(sum(weight * (1 - sensitivity * f)**regions
                                             for f, weight in atoms), result[side])
                    if regions <= 2 or sensitivity == 0:
                        self.assertEqual(result["lower"], result["upper"])

    def test_t3_maximum_zero_visibility_fraction(self):
        # The upper-extremizer's mass at F=0 is the sharp Cauchy bound;
        # an arbitrary compatible population may have much less zero mass.
        laws = (
            ((F(0), F(1, 6)), (F(3, 5), F(5, 6))),
            ((F(2, 5), F(1)),),
            ((F(0), F(1, 5)), (F(1, 3), F(2, 5)), (F(1), F(2, 5))),
            ((F(1), F(1)),),
        )
        for law in laws:
            m = sum(value * weight for value, weight in law)
            s = sum(value * value * weight for value, weight in law)
            maximum_invisible = 1 - m*m/s
            observed_invisible = sum(weight for value, weight in law if value == 0)
            self.assertLessEqual(observed_invisible, maximum_invisible)
            upper = exact.sharp_miss_bounds(m, s, 7, 1)["upper_atoms"]
            self.assertEqual(sum(weight for value, weight in upper if value == 0),
                             maximum_invisible)
        self.assertEqual(exact.sharp_miss_bounds(0, 0, 7, 1)["upper"], 1)

    def test_t3_dual_polynomial_coefficients_and_global_sign_certificates(self):
        for first, second in ((F(1, 2), F(3, 10)), (F(1, 3), F(1, 9)),
                              (F(2, 3), F(2, 3))):
            for regions in (3, 4, 7, 12):
                for sensitivity in (F(1, 3), F(1)):
                    g = tuple(F(math.comb(regions, j)) * (-sensitivity)**j
                              for j in range(regions + 1))
                    for side in ("upper", "lower"):
                        certificate = exact.dual_certificate(first, second, regions,
                                                             sensitivity, side)
                        node = second / first if side == "upper" else (first - second) / (1 - first)
                        endpoint = F(0) if side == "upper" else F(1)
                        coefficients = solve_square(
                            ((1, node, node * node), (0, 1, 2 * node),
                             (1, endpoint, endpoint * endpoint)),
                            ((1 - sensitivity * node)**regions,
                             -regions * sensitivity * (1 - sensitivity * node)**(regions - 1),
                             (1 - sensitivity * endpoint)**regions)
                        )
                        self.assertEqual(certificate["dual"], trim_polynomial(coefficients))
                        difference = []
                        for j in range(regions + 1):
                            qj = coefficients[j] if j < 3 else F(0)
                            difference.append(qj - g[j] if side == "upper" else g[j] - qj)
                        self.assertTrue(certificate["identity"])
                        self.assertEqual(trim_polynomial(difference), multiply_polynomials(
                            certificate["factor"], certificate["quotient"]
                        ))
                        self.assertTrue(all(weight >= 0 for weight in certificate["positive_weights"]))
                        # Coefficient identity plus nonnegative basis on [0,1]
                        # is a global certificate; these points also check coding.
                        for x in (F(0), F(1, 7), node, F(4, 5), F(1)):
                            self.assertEqual(polynomial_value(certificate["quotient"], x),
                                             sum(weight * (1 - sensitivity * x)**j
                                                 for j, weight in enumerate(certificate["positive_weights"])))
                            self.assertGreaterEqual(polynomial_value(difference, x), 0)

    def test_t4_depth_region_minima_and_irreducible_error_floor(self):
        cases = ((F(1, 2), F(3, 10), F(1, 2), 4, F(1, 5)),
                 (F(1, 2), F(1, 4), F(1, 3), 3, F(1, 20)),
                 (F(2, 3), F(2, 3), F(1), 1, F(1, 3)),
                 (F(1), F(1), F(1), 1, F(0)))
        for first, second, vaf, depth, alpha in cases:
            result = exact.minimum_regions(first, second, vaf, depth, alpha, 100)
            sensitivity = 1 - (1 - vaf)**depth
            brute = next(r for r in range(101)
                         if exact.sharp_miss_bounds(first, second, r, sensitivity)["upper"] <= alpha)
            self.assertEqual(result["status"], "attained")
            self.assertEqual(result["regions"], brute)
            self.assertLessEqual(result["upper"], alpha)
            if brute:
                self.assertGreater(result["previous_upper"], alpha)
        depth_result = exact.minimum_depth(F(1, 2), F(1, 4), F(1, 5), 5, F(1, 10), 100)
        brute_depth = next(d for d in range(101) if
                           (1 - F(1, 2) * (1 - F(4, 5)**d))**5 <= F(1, 10))
        self.assertEqual(depth_result["depth"], brute_depth)
        for alpha in (F(1, 7), F(1, 6)):
            self.assertEqual(exact.minimum_regions(F(1, 2), F(3, 10), F(1, 2), 4,
                                                    alpha)["status"], "unreachable")
        self.assertEqual(exact.minimum_depth(F(1, 2), F(1, 4), F(1, 2), 2,
                                              F(1, 4))["status"], "unreachable")
        self.assertEqual(exact.minimum_depth(F(1, 2), F(1, 4), F(1), 2,
                                              F(1, 4))["depth"], 1)
        self.assertEqual(exact.minimum_regions(F(1), F(1), F(0), 100, F(1, 10))["status"],
                         "unreachable")
        self.assertEqual(exact.minimum_regions(F(0), F(0), F(1), 1, F(1))["regions"], 0)

    def test_t4_pilot_inversion_robustness_and_union_bound(self):
        first, second, c0 = F(1, 2), F(3, 10), F(3, 4)
        miss1 = 1 - c0 * first
        miss2 = 1 - 2 * c0 * first + c0 * c0 * second
        self.assertEqual(exact.pilot_moments(miss1, miss2, c0), (first, second))
        h1, h2 = c0 * first, c0 * c0 * second
        epsilon = F(1, 25)
        result = exact.pilot_upper_bound(h1, h2, c0, epsilon, epsilon, 7, F(4, 5))
        self.assertEqual(result["status"], "bounded")
        # Exhaust many feasible latent moment pairs inside the one-sided set.
        for numerator in range(21):
            m = F(numerator, 20)
            for numerator2 in range(41):
                s = F(numerator2, 40)
                if m*m <= s <= m and c0*m >= h1-epsilon and c0*c0*s <= h2+epsilon:
                    self.assertLessEqual(exact.sharp_miss_bounds(m, s, 7, F(4, 5))["upper"],
                                         result["upper"])
        self.assertEqual(exact.pilot_upper_bound(F(9, 10), F(1, 10), 1,
                                                F(1, 100), F(1, 100), 3, 1)["status"],
                         "inconsistent")
        self.assertEqual(exact.joint_error_bound(F(1, 5), F(1, 10), F(1, 20)), F(7, 20))
        self.assertEqual(exact.joint_error_bound(F(3, 4), F(1, 2)), 1)
        with localcontext() as context:
            context.prec = 80
            expected = (-2 * D(100) * dec(epsilon)**2).exp()
            expected += (-2 * D(200) * dec(epsilon)**2).exp()
            self.assert_decimal_enclosed(exact.pilot_failure_bounds(100, 200, epsilon, epsilon),
                                         min(D(1), expected))

    def test_t5_finite_pool_law_by_literal_subset_enumeration(self):
        for size in (2, 3, 5, 8):
            for invisible in range(size + 1):
                for regions in range(size + 1):
                    for sensitivity in (F(0), F(1, 3), F(1)):
                        self.assertEqual(
                            exact.finite_pool_miss(size, invisible, regions, sensitivity),
                            finite_pool_by_subsets(size, invisible, regions, sensitivity)
                        )

    def test_t5_sharp_lp_against_independent_fraction_elimination(self):
        for size in (3, 5, 7):
            for law in (((0, F(1, 3)), (size, F(2, 3))),
                        ((1, F(1, 4)), (size // 2, F(1, 2)), (size - 1, F(1, 4))),
                        ((size // 2, F(1)),)):
                first = sum(weight * F(u, size) for u, weight in law)
                second = sum(weight * F(u*(u-1), size*(size-1)) for u, weight in law)
                for regions, sensitivity in ((0, F(1)), (2, F(1, 2)), (size, F(1))):
                    expected = finite_pool_lp_oracle(size, first, second, regions, sensitivity)
                    result = exact.finite_pool_extrema(size, first, second, regions, sensitivity)
                    self.assertEqual(result["status"], "attained")
                    self.assertEqual((result["lower"], result["upper"]), expected)
                    for side in ("lower", "upper"):
                        atoms = result[side + "_atoms"]
                        self.assertLessEqual(len(atoms), 3)
                        self.assertEqual(sum(weight for _, weight in atoms), 1)
                        self.assertEqual(sum(weight * F(u, size) for u, weight in atoms), first)
                        self.assertEqual(sum(weight * F(u*(u-1), size*(size-1))
                                             for u, weight in atoms), second)
        self.assertEqual(exact.finite_pool_extrema(4, F(1, 2), F(0), 2, 1)["status"],
                         "infeasible")

    def test_t6_exact_clone_coverage_by_state_recursion(self):
        for probabilities in ((), (F(0),), (F(1),), (F(1, 3), F(2, 3)),
                              (F(1, 10), F(1, 5), F(3, 10)),
                              (F(0), F(1, 4), F(1, 2))):
            for draws in range(10):
                self.assertEqual(exact.occupancy_all_seen(probabilities, draws),
                                 categorical_coverage(probabilities, draws))
        rng = random.Random(2026111706)
        draws = 4000
        observed = sum(simulation.clone_trial(rng, [0.5], 1, 1, 1) for _ in range(draws)) / draws
        allowance = math.sqrt(math.log(2e10) / (2 * draws))
        self.assertLessEqual(abs(observed - .5), allowance)
        self.assertEqual(simulation.clone_trial(rng, [], 10, 1, 1), 0)
        self.assertEqual(simulation.clone_trial(rng, [0], 10, 1, 1), 0)
        with self.assertRaises(ValueError):
            simulation.clone_trial(rng, [.6, .6], 1, 1, 1)

    def test_t6_clone_support_bound_and_smallest_guaranteed_budget(self):
        probabilities = (F(1, 10), F(1, 5), F(3, 10))
        for draws in (0, 1, 2, 5, 10, 20):
            exact_miss = 1 - categorical_coverage(probabilities, draws)
            bound = exact.occupancy_miss_bound(3, F(1, 10), F(1), draws)
            self.assertLessEqual(exact_miss, bound)
        result = exact.occupancy_budget(3, F(1, 10), F(1, 2), 4, F(1, 20), 100)
        expected = next(n for n in range(101)
                        if 3 * (1 - F(1, 10) * F(15, 16))**n <= F(1, 20))
        self.assertEqual(result["regions"], expected)
        self.assertEqual(exact.occupancy_budget(0, 0, 0, 0, 0)["regions"], 0)
        self.assertEqual(exact.occupancy_budget(1, 0, 1, 1, F(1, 10))["status"], "unreachable")

    def test_published_aggregate_tail_and_necessary_not_sufficient_region_floor(self):
        aggregate = exact.published_aggregate_check()
        numerator = sum(math.comb(90, count) for count in range(75, 91))
        expected_tail = F(numerator, 2**90)
        self.assertEqual(aggregate["n"], 90)
        self.assertEqual(aggregate["observed"], 75)
        self.assertEqual(aggregate["null_p"], F(1, 2))
        self.assertEqual(aggregate["model_check_upper_tail"], expected_tail)
        self.assertEqual(expected_tail, F(1774990439854189, 38685626227668133590597632))
        # Enumerate binary outcome vectors as a separate small-n tail oracle.
        for size in (0, 1, 4, 7):
            for probability in (F(0), F(1, 3), F(1, 2), F(1)):
                for threshold in range(-1, size + 3):
                    expected = F(0)
                    for mask in range(1 << size):
                        successes = mask.bit_count()
                        if successes >= threshold:
                            expected += probability**successes * (1 - probability)**(size - successes)
                    self.assertEqual(exact.binomial_upper_tail(size, threshold, probability), expected)
        self.assertEqual(aggregate["observed_single_region_miss"], F(5, 6))
        self.assertEqual(aggregate["homogeneous_floor_R5percent"], 17)
        self.assertEqual(aggregate["homogeneous_floor"], F(5, 6)**17)
        self.assertEqual(aggregate["previous_homogeneous_floor"], F(5, 6)**16)
        self.assertLessEqual(F(5, 6)**17, aggregate["target"])
        self.assertGreater(F(5, 6)**16, aggregate["target"])
        # A common first moment has radically different repetition behavior:
        # constant miss probability 5/6 versus all-or-none patient visibility.
        # This verifies why the computed 17 is a necessary floor only.
        for regions in (1, 2, 17, 100):
            homogeneous_miss = F(5, 6)**regions
            heterogeneous_miss = F(1, 6) * F(0)**regions + F(5, 6) * F(1)**regions
            self.assertLessEqual(homogeneous_miss, heterogeneous_miss)
            self.assertEqual(heterogeneous_miss, F(5, 6))
            self.assertGreater(heterogeneous_miss, aggregate["target"])

    def test_rejects_infeasible_or_unstated_inputs(self):
        calls = (
            lambda: exact.exp_neg_bounds(F(-1)),
            lambda: exact.exp_neg_bounds(F(1), F(0)),
            lambda: exact.poisson_cdf_bounds(1, True),
            lambda: exact.sharp_miss_bounds(F(1, 2), F(1, 5), 3, 1),
            lambda: exact.sharp_miss_bounds(F(1, 2), F(3, 5), 3, 1),
            lambda: exact.detection_probability(F(1, 2), -1),
            lambda: exact.detection_probability(F(2), 3),
            lambda: exact.finite_pool_miss(4, 5, 2, 1),
            lambda: exact.finite_pool_miss(4, 2, 5, 1),
            lambda: exact.occupancy_all_seen((F(2, 3), F(2, 3)), 3),
            lambda: exact.anchor_clock_rates(0, 10, 5, 1, 0),
            lambda: exact.anchor_clock_rates(1, 10, 5, 0, 0),
            lambda: exact.anchor_clock_rates(1, 10, 1, 2, 0),
            lambda: exact.pilot_moments(F(1, 2), F(1, 4), 0),
            lambda: exact.minimal_threshold_exposure(1, 2, F(1, 10)),
            lambda: exact.binomial_upper_tail(-1, 0, F(1, 2)),
            lambda: exact.binomial_upper_tail(3, True, F(1, 2)),
            lambda: exact.binomial_upper_tail(3, 1, F(2)),
        )
        for call in calls:
            with self.assertRaises(ValueError):
                call()
        with self.assertRaises(TypeError):
            exact.exp_neg_bounds(0.1)

    def test_independent_fixed_seed_stochastic_agreement_for_all_theorems(self):
        report = simulation.run_simulations(seed=20261117, replicates=4000)
        checks = report["checks"]
        self.assertEqual({record["theorem"] for record in checks},
                         {"T1", "T2", "T3", "T4", "T5", "T6"})
        self.assertGreaterEqual(len(checks), 20)
        for record in checks:
            # Predeclared 7 analytic standard errors plus a .001 numerical
            # allowance. This is a simulation diagnostic, not a confidence
            # interval, formal proof, or evidence of biological validity.
            allowance = 7 * math.sqrt(record["variance"] / record["trials"]) + .001
            self.assertLessEqual(abs(record["observed"] - record["expected"]), allowance,
                                 record["theorem"] + ": " + record["name"])
        path = report["primary_path"]
        self.assertGreater(path["population_events"], 1000)
        self.assertGreater(path["final_cells"], 1000)
        self.assertEqual(path["extinct_runs_discarded"], 0)
        self.assertIsNone(path["population_cap"])
        # Cross-check the independent simulator's analytical targets against
        # exact formulas, rather than comparing the simulator only to itself.
        clock_errors = exact.threshold_errors(5, 7, 4, 27)
        indexed = {record["name"]: record for record in checks}
        self.assertAlmostEqual(indexed["early_boundary_error"]["expected"],
                               float(clock_errors["early"].lower), places=13)
        self.assertAlmostEqual(indexed["late_boundary_error"]["expected"],
                               float(clock_errors["late"].lower), places=13)
        sampling = exact.sharp_miss_bounds(F(1, 2), F(5, 16), 8, F(2101, 3125))
        for side in ("lower", "upper"):
            self.assertAlmostEqual(report["sampling_fixture"][side],
                                   float(sampling[side]), places=13)
        finite = exact.finite_pool_extrema(6, F(3, 8), F(1, 5), 3, F(3, 4))
        for side in ("lower", "upper"):
            self.assertEqual(F(report["finite_pool_fixture"][side]), finite[side])
        clone_recovery = categorical_coverage((F(3, 8), F(9, 40), F(3, 20)), 12)
        self.assertAlmostEqual(report["clone_fixture"]["exact_recovery"],
                               float(clone_recovery), places=13)
        self.assertTrue(report["timing_fixture"]["confounding"]["coupled_observations_equal"])
        self.assertGreater(report["sampling_fixture"]["linked_region_expected"],
                           3 * report["sampling_fixture"]["iid_expected_at_fixed_F_half"])
        pilot = report["pilot_fixture"]
        pilot_bound = exact.pilot_failure_bounds(50, 50, F(3, 20), F(3, 20))
        self.assertAlmostEqual(pilot["hoeffding_union_upper"], float(pilot_bound.upper), places=13)
        self.assertLessEqual(pilot["exact_failure_via_binomial_oracle"], float(pilot_bound.upper))


if __name__ == "__main__":
    unittest.main()
