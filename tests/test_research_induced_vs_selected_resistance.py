"""Independent computational checks of a nonclinical count-identifiability model.

All inputs and random observations are synthetic. Passing tests establishes
neither biological validity nor treatment guidance. Uses the Python stdlib only.
"""

from decimal import Decimal, localcontext
from fractions import Fraction as F
import importlib.util
from itertools import product
import math
from pathlib import Path
import random
import unittest


RESEARCH = (Path(__file__).resolve().parents[1] / "research"
            / "induced-vs-selected-resistance")


def load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, RESEARCH / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


V = load_module("induced_selected_verifier", "verifier.py")
S = load_module("induced_selected_independent_simulation", "simulation.py")


def decimal_fraction(value):
    return Decimal(value.numerator) / Decimal(value.denominator)


def independent_thinning(kernel, retention):
    """Multiply individual Bernoulli polynomials; no binomial coefficients."""
    result = [F(0)] * (max(kernel) + 1)
    for size, probability in kernel.items():
        polynomial = [F(1)]
        for _ in range(size):
            next_polynomial = [F(0)] * (len(polynomial) + 1)
            for index, coefficient in enumerate(polynomial):
                next_polynomial[index] += (1 - retention) * coefficient
                next_polynomial[index + 1] += retention * coefficient
            polynomial = next_polynomial
        for index, coefficient in enumerate(polynomial):
            result[index] += probability * coefficient
    return result


def independent_binomial(n, probability):
    """A Bernoulli-convolution calculation independent of the LR implementation."""
    coefficients = [F(1)]
    for _ in range(n):
        updated = [F(0)] * (len(coefficients) + 1)
        for index, coefficient in enumerate(coefficients):
            updated[index] += (1 - probability) * coefficient
            updated[index + 1] += probability * coefficient
        coefficients = updated
    return coefficients


def independent_best_tests(n, p0, p1, alpha):
    """Enumerate threshold tests and separately solve the finite NP knapsack."""
    null = independent_binomial(n, p0)
    alternative = independent_binomial(n, p1)
    deterministic_power = max(
        sum(alternative[k:], F(0)) for k in range(n + 2)
        if sum(null[k:], F(0)) <= alpha
    )
    budget, randomized_power = alpha, F(0)
    for k in sorted(range(n + 1), key=lambda i: alternative[i] / null[i],
                    reverse=True):
        allowed = min(budget, null[k])
        randomized_power += allowed * alternative[k] / null[k]
        budget -= allowed
        if budget == 0:
            break
    return deterministic_power, randomized_power


class ExactDistributionTests(unittest.TestCase):
    def test_finite_opportunity_convolution_obeys_poisson_scaling_bound(self):
        # Independently expand N Bernoulli opportunities with marked D/I/P
        # events. Truncation leaves an exact finite-law tail mass; the CP tail
        # is enclosed by the verifier. Their sum bounds unenumerated TV mass.
        kernel = {1: F(1, 2), 3: F(1, 2)}
        a, b, c, retention = F(1, 4), F(1, 4), F(1, 2), F(2, 3)
        total_intensity = a + b + c
        thinned = independent_thinning(kernel, retention)
        levy = V.levy_coefficients(kernel, a, b, c, retention)
        limit = 32
        cp_intervals, cp_tail = V.pmf_bounds(levy, limit, terms=24)
        for opportunities in (8, 32, 128):
            with self.subTest(opportunities=opportunities):
                one = [F(0)] * (max(kernel) + 1)
                one[0] = 1 - total_intensity / opportunities
                one[1] += b / opportunities
                for size, probability in kernel.items():
                    one[size] += a * probability / opportunities
                for size, probability in enumerate(thinned):
                    one[size] += c * probability / opportunities
                self.assertEqual(sum(one), 1)
                finite = [F(1)] + [F(0)] * limit
                for _ in range(opportunities):
                    updated = [F(0)] * (limit + 1)
                    for size, probability in enumerate(one):
                        for count in range(limit - size + 1):
                            updated[count + size] += finite[count] * probability
                    finite = updated
                finite_tail = 1 - sum(finite)
                prefix_error_upper = sum(
                    max(abs(probability - lo), abs(probability - hi))
                    for probability, (lo, hi) in zip(finite, cp_intervals)
                )
                tv_upper = (prefix_error_upper + finite_tail + cp_tail[1]) / 2
                self.assertLessEqual(tv_upper, total_intensity**2 / opportunities)
                for z in (F(0), F(1, 3), F(2, 3)):
                    finite_pgf = sum(probability * z**size
                                     for size, probability in enumerate(one))**opportunities
                    lo, hi = V.exp_minus_bounds(-V.log_pgf(levy, z), terms=24)
                    self.assertLessEqual(max(abs(finite_pgf - lo), abs(finite_pgf - hi)),
                                         total_intensity**2 / opportunities)

    def test_thinning_by_individual_cells_and_coefficient_distribution(self):
        kernels = ({1: F(1, 4), 2: F(1, 4), 5: F(1, 2)},
                   {3: F(1, 3), 8: F(2, 3)})
        for kernel, q in product(kernels, (F(0), F(1, 3), F(1))):
            with self.subTest(kernel=kernel, q=q):
                thinned = independent_thinning(kernel, q)
                self.assertEqual(V.thin_kernel(kernel, q), dict(enumerate(thinned)))
                self.assertEqual(sum(thinned), 1)
                a, b, c = F(3, 7), F(2, 5), F(4, 9)
                independent_levy = {
                    k: a * kernel.get(k, 0) + c * thinned[k]
                    + (b if k == 1 else 0)
                    for k in range(1, len(thinned))
                }
                independent_levy = {k: value for k, value in independent_levy.items()
                                    if value}
                levy = V.levy_coefficients(kernel, a, b, c, q)
                self.assertEqual(levy, independent_levy)
                weights = V.reduced_pmf(levy, 22)
                self.assertEqual(weights, V.partition_weights(levy, 22))
                self.assertEqual(V.recover_levy_from_weights(weights), levy)
                intervals, tail = V.pmf_bounds(levy, 22, terms=24)
                self.assertTrue(all(0 <= lo <= hi for lo, hi in intervals))
                self.assertTrue(0 <= tail[0] <= tail[1] <= 1)
                self.assertLessEqual(sum(lo for lo, _ in intervals) + tail[0], 1)
                self.assertGreaterEqual(sum(hi for _, hi in intervals) + tail[1], 1)

    def test_rational_exponential_enclosures_against_independent_decimal(self):
        with localcontext() as context:
            context.prec = 110
            for x, terms in product((F(0), F(1, 7), F(1), F(11, 3), F(27)),
                                    (8, 17, 30)):
                with self.subTest(x=x, terms=terms):
                    lo, hi = V.exp_minus_bounds(x, terms)
                    independent_value = (-decimal_fraction(x)).exp()
                    self.assertLessEqual(decimal_fraction(lo), independent_value)
                    self.assertLessEqual(independent_value, decimal_fraction(hi))
                    self.assertTrue(0 <= lo <= hi <= 1)
            levy = {1: F(3, 5), 3: F(2, 7), 6: F(1, 9)}
            weights = V.partition_weights(levy, 24)
            intervals, tail = V.pmf_bounds(levy, 24, terms=30)
            scale = (-decimal_fraction(sum(levy.values(), F(0)))).exp()
            actual_prefix = Decimal(0)
            for weight, (lo, hi) in zip(weights, intervals):
                value = scale * decimal_fraction(weight)
                self.assertLessEqual(decimal_fraction(lo), value)
                self.assertLessEqual(value, decimal_fraction(hi))
                actual_prefix += value
            self.assertLessEqual(decimal_fraction(tail[0]), 1 - actual_prefix)
            self.assertLessEqual(1 - actual_prefix, decimal_fraction(tail[1]))


class ExactIdentifiabilityTests(unittest.TestCase):
    def test_single_positive_arm_triangular_alias_at_arbitrary_clone_support(self):
        # Back substitution is implemented here independently of any inverse
        # in verifier.py. The equality compares all positive log coefficients,
        # so it proves equality of the full distributions on each tested grid.
        for maximum in (3, 7, 17):
            kernel = {k: F(1, maximum) for k in range(1, maximum + 1)}
            x, intensity, b, p = F(1, 3), F(2), F(1, 4), F(1, 2)
            a, c = x * intensity, (1 - x) * intensity
            observed = V.levy_coefficients(kernel, a, b, c, p)
            changed_p = p + F(1, 100 * maximum * maximum)
            changed_marks = {}
            for j in range(maximum, 0, -1):
                rhs = observed[j] - (b if j == 1 else 0)
                rhs -= (1 - x) * sum(
                    changed_marks[k] * math.comb(k, j) * changed_p**j
                    * (1 - changed_p)**(k - j)
                    for k in range(j + 1, maximum + 1)
                )
                changed_marks[j] = rhs / (x + (1 - x) * changed_p**j)
            self.assertTrue(all(value > 0 for value in changed_marks.values()))
            changed_intensity = sum(changed_marks.values())
            changed_kernel = {k: value / changed_intensity
                              for k, value in changed_marks.items()}
            changed_a, changed_c = x * changed_intensity, (1 - x) * changed_intensity
            self.assertEqual(sum(changed_kernel.values()), 1)
            self.assertEqual(V.levy_coefficients(changed_kernel, changed_a, b,
                                                changed_c, changed_p), observed)
            self.assertNotEqual(p, changed_p)
            self.assertNotEqual(V.levy_coefficients(changed_kernel, changed_a, b,
                                                   changed_c),
                                V.levy_coefficients(kernel, a, b, c))

    def test_two_arm_reconstruction_general_grid_and_boundary_intensities(self):
        kernels = ({1: F(1, 4), 2: F(1, 4), 4: F(1, 2)},
                   {1: F(1, 7), 3: F(6, 7)}, {64: F(1)})
        for kernel, a, b, c, q in product(
            kernels, (F(0), F(1, 2), F(3)), (F(0), F(2, 7)),
            (F(1, 3), F(2)), (F(0), F(1, 5), F(2, 3), F(9, 10))
        ):
            with self.subTest(kernel=kernel, a=a, b=b, c=c, q=q):
                baseline = V.levy_coefficients(kernel, a, b, c)
                holiday = V.levy_coefficients(kernel, a, b, c, q)
                inferred = V.recover_two_arm(V.cumulants(baseline),
                                            V.cumulants(holiday), sum(baseline.values()))
                self.assertEqual(inferred, {"a": a, "b": b, "c": c, "p": q})
                self.assertEqual(V.recover_kernel(baseline, a, b, c), kernel)

    def test_three_arm_two_way_inverse_needs_only_second_cumulants(self):
        kernels = ({1: F(1, 2), 2: F(1, 2)},
                   {1: F(1, 7), 3: F(2, 7), 5: F(4, 7)}, {64: F(1)})
        parameters = ((F(0), F(0), F(1)), (F(1, 3), F(2, 5), F(2, 3)))
        for kernel, (a, b, c), pi, r in product(
            kernels, parameters, (F(0), F(1, 4), F(9, 10)),
            (F(1, 7), F(1, 2), F(4, 5))
        ):
            with self.subTest(kernel=kernel, a=a, b=b, c=c, pi=pi, r=r):
                arms = [V.levy_coefficients(kernel, a, b, c,
                                          pi + (1 - pi) * r**j) for j in range(3)]
                cumulants = [V.cumulants(arm, order=2) for arm in arms]
                inferred = V.recover_three_arm(*cumulants, sum(arms[0].values()))
                self.assertEqual(inferred, {"a": a, "b": b, "c": c,
                                           "p": pi + (1 - pi) * r,
                                           "pi": pi, "r": r})
                self.assertEqual(V.recover_kernel(arms[0], a, b, c), kernel)

    def test_baseline_preexisting_and_reversible_clones_have_same_entire_law(self):
        for maximum in (3, 9, 32, 64):
            kernel = {k: F(1, maximum) for k in range(1, maximum + 1)}
            selected = V.levy_coefficients(kernel, F(7, 5), 0, 0)
            reversible = V.levy_coefficients(kernel, 0, 0, F(7, 5))
            self.assertEqual(selected, reversible)
            self.assertEqual(V.reduced_pmf(selected, 12),
                             V.reduced_pmf(reversible, 12))
            after_holiday = V.levy_coefficients(kernel, 0, 0, F(7, 5), F(1, 2))
            self.assertNotEqual(selected, after_holiday)

    def test_doubleton_continuum_alias_and_third_arm_separation(self):
        kernel = {1: F(1, 2), 2: F(1, 2)}
        baseline = V.levy_coefficients(kernel, 1, 1, 1)
        holiday = V.levy_coefficients(kernel, 1, 1, 1, F(1, 2))
        self.assertEqual(baseline, {1: F(2), 2: F(1)})
        self.assertEqual(holiday, {1: F(2), 2: F(5, 8)})
        for q in (F(1, 3), F(2, 5), F(3, 5)):
            alias = V.doubleton_alias(baseline, holiday, q)
            for retention, expected in ((F(1), baseline), (q, holiday)):
                actual = V.levy_coefficients(alias["kernel"], alias["a"],
                                             alias["b"], alias["c"], retention)
                self.assertEqual(actual, expected)
            self.assertNotEqual(
                V.levy_coefficients(kernel, 1, 1, 1, F(1, 4)),
                V.levy_coefficients(alias["kernel"], alias["a"], alias["b"],
                                    alias["c"], q * q),
            )
        alias = V.doubleton_alias(baseline, holiday, F(1, 3))
        self.assertEqual(alias, {"a": F(185, 192), "b": F(4, 3),
                                 "c": F(45, 64), "q": F(1, 3),
                                 "kernel": {1: F(2, 5), 2: F(3, 5)}})

    def test_singleton_four_arm_recovery_and_two_way_all_time_ambiguity(self):
        for a, b, c, q in product((F(0), F(2)), (F(0), F(3, 7)),
                                  (F(1, 2), F(4)), (F(0), F(1, 3), F(4, 5))):
            means = [a + b + c * q**j for j in range(3)] + [a + 2 * b + c]
            self.assertEqual(V.recover_singleton_four(means),
                             {"a": a, "b": b, "c": c, "q": q})
        # Both parameter sets have a+c*pi=5/2 and c*(1-pi)=3/2.
        # Consequently every real-time occupancy law is the same exponential
        # curve; the rational points below are exact checks of that identity.
        for exposure, decay in product((F(0), F(1), F(3)),
                                       (F(0), F(1, 9), F(1, 2), F(1))):
            first = 2 + F(1, 2) * exposure + 2 * (F(1, 4) + F(3, 4) * decay)
            second = 1 + F(1, 2) * exposure + 3 * (F(1, 2) + F(1, 2) * decay)
            self.assertEqual(first, second)

    def test_two_switching_laws_agree_at_one_holiday_but_not_two(self):
        kernel = {1: F(1, 3), 3: F(2, 3)}
        first = (F(1, 4), F(1, 3))
        second = (F(1, 3), F(1, 4))
        laws = []
        for pi, r in (first, second):
            laws.append([V.levy_coefficients(kernel, F(2, 3), F(1, 4), F(3, 2),
                                            pi + (1 - pi) * r**j)
                         for j in range(3)])
        self.assertEqual(laws[0][0], laws[1][0])
        self.assertEqual(laws[0][1], laws[1][1])
        self.assertNotEqual(laws[0][2], laws[1][2])

    def test_two_positive_holiday_arms_do_not_replace_the_immediate_arm(self):
        # An unbounded geometric-kernel alias: equality holds as a rational
        # function of z at both H and 2H, despite different switching rates.
        def geometric_pgf(z, probability):
            return probability * z / (1 - (1 - probability) * z)

        for holiday, z in product((1, 2), (F(0), F(1, 4), F(1, 2), F(3, 4))):
            p = F(1, 4) + F(3, 4) * F(1, 2)**holiday
            p_other = F(11, 136) + F(125, 136) * F(33, 50)**holiday
            first = geometric_pgf(1 - p + p * z, F(1, 2)) - 1
            second = F(20, 21) * (
                geometric_pgf(1 - p_other + p_other * z, F(11, 21)) - 1
            )
            self.assertEqual(first, second)
        self.assertNotEqual(-F(1), -F(20, 21))  # Immediate zero probabilities differ.

    def test_nonidentifiable_boundaries_are_rejected(self):
        # Admissible-looking higher-cumulant ratios do not excuse a mean
        # increase: positive clone sizes and loss require a strict mean drop.
        with self.assertRaises(ValueError):
            V.recover_two_arm((1, 1, 1), (2, F(1, 4), F(1, 8)), 4)
        kernel = {1: F(1, 4), 2: F(1, 4), 4: F(1, 2)}
        for c, p in ((F(0), F(1, 2)), (F(1), F(1))):
            baseline = V.levy_coefficients(kernel, 1, 1, c)
            holiday = V.levy_coefficients(kernel, 1, 1, c, p)
            with self.assertRaises(ValueError):
                V.recover_two_arm(V.cumulants(baseline), V.cumulants(holiday),
                                  sum(baseline.values()))
        for kernel in ({1: F(1)}, {1: F(1, 2), 2: F(1, 2)}):
            baseline = V.levy_coefficients(kernel, 1, 1, 1)
            holiday = V.levy_coefficients(kernel, 1, 1, 1, F(1, 2))
            with self.assertRaises(ValueError):
                V.recover_two_arm(V.cumulants(baseline), V.cumulants(holiday),
                                  sum(baseline.values()))


class ExactPowerTests(unittest.TestCase):
    def test_certified_tv_enclosures_obey_both_coupling_obstructions(self):
        kernel = {1: F(1, 3), 3: F(2, 3)}

        def pgf(z):
            return sum(mass * z**size for size, mass in kernel.items())

        def tv_upper(first, second):
            first_prefix, first_tail = V.pmf_bounds(first, 30, terms=24)
            second_prefix, second_tail = V.pmf_bounds(second, 30, terms=24)
            prefix = sum(max(abs(lo0 - hi1), abs(hi0 - lo1))
                         for (lo0, hi0), (lo1, hi1) in zip(first_prefix, second_prefix))
            return (prefix + first_tail[1] + second_tail[1]) / 2

        common_a, common_b = F(1, 2), F(1, 3)
        no_persisters = V.levy_coefficients(kernel, common_a, common_b, 0)
        previous_lower = F(1)
        for c in (F(1, 2), F(1, 8), F(1, 32), F(1, 128)):
            p = F(3, 5)
            added = V.levy_coefficients(kernel, common_a, common_b, c, p)
            exponent = c * (1 - pgf(1 - p))
            exp_lo, exp_hi = V.exp_minus_bounds(exponent, terms=24)
            coupling_lower, coupling_upper = 1 - exp_hi, 1 - exp_lo
            self.assertLessEqual(tv_upper(no_persisters, added), coupling_lower)
            self.assertLessEqual(coupling_upper, exponent)
            self.assertLess(coupling_upper, previous_lower)
            previous_lower = coupling_lower

        c = F(1, 2)
        permanent = V.levy_coefficients(kernel, common_a + c, common_b, 0)
        previous_lower = F(1)
        for p in (F(1, 2), F(7, 8), F(31, 32), F(127, 128)):
            reversible = V.levy_coefficients(kernel, common_a, common_b, c, p)
            # This exponent uses G(p), whereas addition used G(1-p).
            exponent = c * (1 - pgf(p))
            exp_lo, exp_hi = V.exp_minus_bounds(exponent, terms=24)
            coupling_lower, coupling_upper = 1 - exp_hi, 1 - exp_lo
            self.assertLessEqual(tv_upper(permanent, reversible), coupling_lower)
            self.assertLessEqual(coupling_upper, exponent)
            self.assertLess(coupling_upper, previous_lower)
            previous_lower = coupling_lower

    def test_bounded_pgf_sample_bound_has_an_exact_exponential_certificate(self):
        kernel = {1: F(1, 3), 3: F(2, 3)}
        null = V.levy_coefficients(kernel, F(1, 2), F(1, 3), 0)
        alternative = V.levy_coefficients(kernel, 0, F(1, 3), F(1, 2), F(1, 2))
        z = F(1, 2)
        f0_lo, f0_hi = V.exp_minus_bounds(-V.log_pgf(null, z), terms=16)
        f1_lo, f1_hi = V.exp_minus_bounds(-V.log_pgf(alternative, z), terms=16)
        delta_lower = f1_lo - f0_hi
        delta_upper = f1_hi - f0_lo
        self.assertTrue(0 < delta_lower <= delta_upper <= 1)
        alpha, beta = F(1, 20), F(1, 10)
        # Floating arithmetic proposes an integer; only the subsequent exact
        # rational enclosure certifies the claimed sufficient sample count.
        n = math.ceil(2 * max(math.log(1 / float(alpha)), math.log(1 / float(beta)))
                      / float(delta_lower)**2)
        while V.exp_minus_bounds(n * delta_lower**2 / 2, terms=16)[1] > min(alpha, beta):
            n += 1
        certified_error_upper = V.exp_minus_bounds(n * delta_lower**2 / 2, terms=16)[1]
        self.assertLessEqual(certified_error_upper, alpha)
        self.assertLessEqual(certified_error_upper, beta)
        self.assertLess(n, 10000)

    def test_minimum_replicate_certificates_and_all_smaller_designs(self):
        p0, p1, alpha, target = F(1, 4), F(1, 2), F(1, 20), F(4, 5)
        deterministic = V.minimum_zero_test(p0, p1)
        randomized = V.minimum_zero_test(p0, p1, randomized=True)
        self.assertEqual((deterministic["n"], deterministic["boundary"],
                          deterministic["gamma"]), (26, 10, F(0)))
        self.assertEqual(deterministic["size"], F(22565567197793, 562949953421312))
        self.assertEqual(deterministic["power"], F(1754331, 2097152))
        self.assertEqual(randomized["n"], 23)
        self.assertEqual(randomized["boundary"], 9)
        self.assertEqual(randomized["gamma"], F(1622209696738, 9771486092775))
        self.assertEqual(randomized["size"], alpha)
        self.assertEqual(randomized["power"], F(40811361573509, 50153065021440))
        for n in range(1, 27):
            deterministic_power, randomized_power = independent_best_tests(n, p0, p1, alpha)
            self.assertEqual(deterministic_power, V.best_zero_test(n, p0, p1)["power"])
            self.assertEqual(randomized_power,
                             V.best_zero_test(n, p0, p1, randomized=True)["power"])
            if n < 26:
                self.assertLess(deterministic_power, target)
            if n < 23:
                self.assertLess(randomized_power, target)


class ExactAssumptionFailureTests(unittest.TestCase):
    def test_distinct_unknown_clone_laws_allow_stable_singleton_reallocation(self):
        stable = {1: F(1, 4), 4: F(3, 4)}
        reversible = {1: F(1, 3), 2: F(1, 3), 7: F(1, 3)}
        a, b, c, moved = F(2), F(1, 3), F(5, 4), F(1, 4)
        changed_a, changed_b = a - moved, b + moved
        changed_stable = {k: (a * mass - (moved if k == 1 else 0)) / changed_a
                          for k, mass in stable.items()}
        self.assertEqual(changed_stable, {1: F(1, 7), 4: F(6, 7)})
        for p in (F(0), F(1, 7), F(1, 2), F(9, 10), F(1)):
            thinned = independent_thinning(reversible, p)
            original_levy = {
                k: a * stable.get(k, 0) + c * thinned[k] + (b if k == 1 else 0)
                for k in range(1, len(thinned))
            }
            changed_levy = {
                k: changed_a * changed_stable.get(k, 0) + c * thinned[k]
                + (changed_b if k == 1 else 0)
                for k in range(1, len(thinned))
            }
            self.assertEqual(original_levy, changed_levy)

    def test_coherent_whole_clone_switching_has_all_time_floor_alias(self):
        kernel = {1: F(1, 7), 5: F(2, 7), 32: F(4, 7)}
        # For whole-clone switching, each reversible event contributes its
        # entire clone or zero. Its intensity is c*p, not cell-thinned marks.
        for decay in (F(0), F(1, 13), F(1, 3), F(4, 5), F(1)):
            p = F(1, 4) + F(3, 4) * decay
            changed_p = F(1, 2) + F(1, 2) * decay
            original = {k: (2 + 2 * p) * mass + (F(1, 2) if k == 1 else 0)
                        for k, mass in kernel.items()}
            changed = {k: (1 + 3 * changed_p) * mass + (F(1, 2) if k == 1 else 0)
                       for k, mass in kernel.items()}
            self.assertEqual(original, changed)
        independently_switched = V.levy_coefficients(kernel, 2, F(1, 2), 2, F(1, 2))
        coherently_switched = {k: 3 * mass + (F(1, 2) if k == 1 else 0)
                              for k, mass in kernel.items()}
        self.assertNotEqual(independently_switched, coherently_switched)

    def test_common_unknown_plating_recovers_only_effective_intensities(self):
        kernel = {1: F(1, 5), 3: F(2, 5), 7: F(2, 5)}
        a, b, c = F(3, 4), F(1, 2), F(1, 3)
        for efficiency in (F(1, 5), F(2, 3), F(1)):
            observed_mark = independent_thinning(kernel, efficiency)
            survival = 1 - observed_mark[0]
            effective_kernel = {k: mass / survival
                                for k, mass in enumerate(observed_mark) if k and mass}
            effective_a, effective_b, effective_c = (
                a * survival, b * efficiency, c * survival
            )
            baseline = V.levy_coefficients(effective_kernel, effective_a,
                                           effective_b, effective_c)
            for retention in (F(1, 4), F(7, 10), F(1)):
                observed_reversible = independent_thinning(kernel, efficiency * retention)
                direct = {
                    k: a * observed_mark[k] + c * observed_reversible[k]
                    + (b * efficiency if k == 1 else 0)
                    for k in range(1, len(observed_mark))
                }
                direct = {k: value for k, value in direct.items() if value}
                effective = V.levy_coefficients(effective_kernel, effective_a,
                                                effective_b, effective_c, retention)
                self.assertEqual(direct, effective)
                if retention < 1:
                    inferred = V.recover_two_arm(V.cumulants(baseline),
                                                V.cumulants(direct), sum(baseline.values()))
                    self.assertEqual(inferred, {"a": effective_a, "b": effective_b,
                                               "c": effective_c, "p": retention})
            if efficiency < 1:
                self.assertNotEqual((effective_a, effective_b, effective_c), (a, b, c))


class IndependentStochasticTests(unittest.TestCase):
    def assert_bounded_pgf(self, counts, log_probability):
        # Union bound over at most 24 bounded statistics, family error <=1e-8.
        tolerance = math.sqrt(math.log(48 / 1e-8) / (2 * len(counts)))
        for z in (0.0, 0.5, 2 / 3):
            empirical = math.fsum(z**count for count in counts) / len(counts)
            expected = math.exp(log_probability(z))
            self.assertLessEqual(abs(empirical - expected), tolerance,
                                 (z, empirical, expected, tolerance))

    def test_event_samples_match_rational_mixture_pgfs(self):
        kernel = {1: F(1, 4), 2: F(1, 4), 4: F(1, 2)}
        a, b, c = F(3, 5), F(1, 3), F(4, 5)
        for index, q in enumerate((F(0), F(1, 2), F(1))):
            counts = S.sample_counts(30000, 20261117 + index, a, b, c,
                                     clone_pmf=kernel, retention=q)
            levy = V.levy_coefficients(kernel, a, b, c, q)
            self.assert_bounded_pgf(counts, lambda z: sum(
                float(value) * (z**k - 1) for k, value in levy.items()))

    def test_direct_switching_clocks_match_three_endpoint_pgfs(self):
        kernel = {1: F(1, 4), 2: F(1, 4), 4: F(1, 2)}
        a, b, c = 0.6, 1 / 3, 0.8

        def mark_pgf(z):
            return math.fsum(float(value) * z**k for k, value in kernel.items())

        for index, duration in enumerate((0.0, 1.0, 2.0)):
            counts = S.sample_counts(30000, 20261217 + index, a, b, c,
                                     clone_pmf=kernel, holiday_time=duration,
                                     forward_rate=0.25, reverse_rate=0.75)
            # This formula is deliberately used only in the test. The simulator
            # independently follows random exponential switching clocks.
            p = 0.25 + 0.75 * math.exp(-duration)
            self.assert_bounded_pgf(counts, lambda z: (
                a * (mark_pgf(z) - 1) + b * (z - 1)
                + c * (mark_pgf(1 - p + p * z) - 1)
            ))

    def test_unbounded_clone_simulation_and_actual_likelihood_test_power(self):
        experiments, per_experiment = 4000, 26
        n = experiments * per_experiment
        intensity = 2 * math.log(2)
        selected = S.sample_counts(n, 10120261117, intensity, 0, 0, geometric_p=0.5)
        reversible = S.sample_counts(n, 10220261117, 0, 0, intensity,
                                     geometric_p=0.5, retention=1 / 3)
        certificate = V.minimum_zero_test(F(1, 4), F(1, 2))
        random_certificate = V.minimum_zero_test(F(1, 4), F(1, 2), randomized=True)
        rng = random.Random(10320261117)
        for counts, zero_probability, expected, randomized_expected in (
            (selected, 0.25, float(certificate["size"]),
             float(random_certificate["size"])),
            (reversible, 0.5, float(certificate["power"]),
             float(random_certificate["power"])),
        ):
            self.assertLess(abs(sum(k == 0 for k in counts) / n - zero_probability), 0.01)
            rejections, randomized_rejections = 0, 0
            for start in range(0, n, per_experiment):
                group = counts[start:start + per_experiment]
                zero_count = sum(k == 0 for k in group)
                rejections += zero_count >= 11
                zero_count = sum(k == 0 for k in group[:23])
                randomized_rejections += (
                    zero_count > 9 or (zero_count == 9
                                      and rng.random() < float(random_certificate["gamma"]))
                )
            self.assertLess(abs(rejections / experiments - expected), 0.04)
            self.assertLess(abs(randomized_rejections / experiments - randomized_expected),
                            0.04)


if __name__ == "__main__":
    unittest.main()
