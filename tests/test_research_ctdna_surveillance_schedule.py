"""Independent checks of NONCLINICAL ctDNA surveillance mathematics.

Synthetic distributions only; no patient data or monitoring advice. Exact
outcome enumeration, independently integrated growth paths and seeded event
simulations can falsify implementations, but do not replace the README proofs
or establish clinical validity or research priority. Standard library only.
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
            / "ctdna-surveillance-schedule")
F = Fraction
D = Decimal


def load_module(name):
    spec = importlib.util.spec_from_file_location(
        "ctdna_schedule_" + name, RESEARCH / (name + ".py")
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


exact = load_module("verifier")
simulation = load_module("simulation")


def all_subsets(width):
    return [indices for size in range(width + 1)
            for indices in combinations(range(width), size)]


def enumerate_outcomes(weights, rows):
    """Whole Bernoulli outcome vectors, not the verifier's product-of-failures."""
    width = len(rows[0])
    distribution = {}
    for mask in range(1 << width):
        probability = F(0)
        for weight, row in zip(weights, rows):
            outcome = weight
            for j, p in enumerate(row):
                outcome *= p if mask & (1 << j) else 1 - p
            probability += outcome
        distribution[mask] = probability
    return distribution


def event_probability(distribution, selected):
    selected_mask = sum(1 << j for j in selected)
    return sum((p for outcome, p in distribution.items()
                if outcome & selected_mask), F(0))


def adaptive_integral(function, left, right, tolerance=2e-12):
    """Independent float adaptive Simpson integration of a smooth integrand."""
    if left == right:
        return 0.0
    midpoint = (left + right) / 2
    fa, fm, fb = function(left), function(midpoint), function(right)
    whole = (right - left) * (fa + 4 * fm + fb) / 6

    def refine(a, b, f_a, f_m, f_b, previous, error, depth):
        c = (a + b) / 2
        fl, fr = function((a + c) / 2), function((c + b) / 2)
        first = (c - a) * (f_a + 4 * fl + f_m) / 6
        second = (b - c) * (f_m + 4 * fr + f_b) / 6
        correction = first + second - previous
        if depth == 0 or abs(correction) <= 15 * error:
            return first + second + correction / 15
        return (refine(a, c, f_a, fl, f_m, first, error / 2, depth - 1)
                + refine(c, b, f_m, fr, f_b, second, error / 2, depth - 1))

    return refine(left, right, fa, fm, fb, whole, tolerance, 25)


def gamma_cdf_independent(shape, value):
    """Direct density quadrature, with y=x**shape only for singular shapes."""
    if value <= 0:
        return 0.0
    if shape >= 1:
        return adaptive_integral(
            lambda x: x**(shape - 1) * math.exp(-x), 0.0, value
        ) / math.gamma(shape)
    return adaptive_integral(
        lambda y: math.exp(-y ** (1 / shape)), 0.0, value**shape
    ) / (shape * math.gamma(shape))


def lognormal_cdf(mu, sigma):
    def cdf(rate):
        if rate <= 0:
            return 0.0
        return (1 + math.erf((math.log(rate) - mu)
                             / (sigma * math.sqrt(2)))) / 2
    return cdf


def gamma_integer_cdf(shape, scale):
    def cdf(rate):
        x = max(0.0, rate / scale)
        return 1 - math.exp(-x) * sum(
            x**j / math.factorial(j) for j in range(shape)
        )
    return cdf


def independent_success(schedule, a, b, cdf, endpoint=None):
    """Partition all event boundaries and evaluate direct burden inequalities.

    Unlike the production union merger, this keeps overlaps and determines the
    success event at each cell's midpoint before integrating that cell's mass.
    """
    times = tuple(map(float, schedule))
    boundaries = {0.0}
    for time in times:
        boundaries.update((a / time, b / time))
    if endpoint is not None:
        boundaries.add(b / endpoint)
    ordered = sorted(boundaries)
    value = 0.0
    for left, right in zip(ordered, ordered[1:]):
        rate = (left + right) / 2
        scored = endpoint is None or b / rate <= endpoint
        if scored and any(a <= rate * time < b for time in times):
            value += cdf(right) - cdf(left)
    return value


def independent_lognormal_lead(schedule, a, b, mu, sigma):
    """Integrate first-detection lead directly over log growth rate."""
    times = sorted(map(float, schedule))
    boundaries = sorted({math.log(v / t) for t in times for v in (a, b)})
    total = 0.0
    for left, right in zip(boundaries, boundaries[1:]):
        midpoint_rate = math.exp((left + right) / 2)
        first = next((t for t in times if a <= midpoint_rate * t < b), None)
        if first is not None:
            total += adaptive_integral(
                lambda x: (b * math.exp(-x) - first)
                * math.exp(-((x - mu) / sigma)**2 / 2)
                / (sigma * math.sqrt(2 * math.pi)), left, right
            )
    return total


def uniform_loss_by_cells(schedule, horizon, window):
    """Exact midpoint integration of piecewise affine event-path loss."""
    times = sorted(map(F, schedule))
    boundaries = sorted({F(0), F(horizon), *times})
    result = F(0)
    for left, right in zip(boundaries, boundaries[1:]):
        onset = (left + right) / 2
        first = next((t for t in times if onset <= t < onset + window), None)
        loss = F(window) if first is None else first - onset
        result += (right - left) * loss / horizon
    return result


class CtDNASurveillanceScheduleTests(unittest.TestCase):
    """Every numbered theorem has independent numerical and seeded checks."""

    @staticmethod
    def bounded_allowance(draws, width):
        # Hoeffding: P(|sample mean-E mean|>epsilon)<=1e-10 for each
        # comparison of independent observations in an interval of this width.
        # Fewer than 100 bounded comparisons give family failure<=1e-8.
        return width * math.sqrt(math.log(2e10) / (2 * draws)) + 2e-11

    def assert_bounded_mc(self, observed, expected, draws, width, label):
        allowance = self.bounded_allowance(draws, width)
        self.assertLessEqual(abs(observed - float(expected)), allowance, label)

    def assert_second_moment_mc(self, observed, expected, second_moment_bound,
                                draws, label):
        # Chebyshev applies to the unbounded lognormal loss because its second
        # moment is bounded analytically by (b-a)^2 E[R^-2]. Each of the three
        # comparisons has failure probability<=.001, hence combined<=.003.
        allowance = math.sqrt(second_moment_bound / (draws * .001)) + 2e-11
        self.assertLessEqual(abs(observed - float(expected)), allowance, label)

    def assert_probability_mc(self, result, expected, label):
        self.assert_bounded_mc(result["probability"], expected, result["draws"], 1, label)

    def test_assay_exact_individual_molecule_enumeration_and_lod_gate(self):
        for molecules in range(1, 8):
            for fraction in (F(0), F(1, 5), F(1, 2), F(4, 5), F(1)):
                for lod in (F(1, 4), F(1, 2), F(4, 5)):
                    probability = F(0)
                    if fraction >= lod:
                        for mask in range(1 << molecules):
                            count = mask.bit_count()
                            if F(count, molecules) >= lod:
                                probability += (fraction**count
                                                * (1 - fraction)**(molecules - count))
                    self.assertEqual(
                        exact.binomial_detection(fraction, lod, molecules), probability
                    )
        self.assertEqual(exact.binomial_detection(F(1, 2), F(1, 2), 2), F(3, 4))
        self.assertEqual(exact.binomial_detection(F(49, 100), F(1, 2), 2), 0)

    def test_t1_whole_outcomes_marginals_submodularity_and_greedy_bound(self):
        weights = (F(1, 10), F(2, 10), F(3, 10), F(4, 10))
        rows = ((F(0), F(1, 2), F(3, 4), F(0), F(1)),
                (F(1, 2), F(0), F(1, 3), F(1), F(0)),
                (F(0), F(2, 3), F(1, 4), F(1, 2), F(1, 5)),
                (F(1, 3), F(1, 2), F(1, 3), F(1, 4), F(0)))
        law = enumerate_outcomes(weights, rows)
        self.assertEqual(sum(law.values()), 1)
        subsets = all_subsets(5)
        values = {s: event_probability(law, s) for s in subsets}
        for selected in subsets:
            self.assertEqual(exact.grid_value(weights, rows, selected), values[selected])
            # Formula (7), checked against whole outcome vectors. Here the
            # synthetic clinical time is 6 and candidate times are j+1.
            enumerated_lead = sum(
                probability * max((5 - j for j in selected if mask & (1 << j)), default=0)
                for mask, probability in law.items()
            )
            first_detection_formula = F(0)
            for weight, row in zip(weights, rows):
                survived = F(1)
                for j in selected:
                    first_detection_formula += weight * (5 - j) * row[j] * survived
                    survived *= 1 - row[j]
            self.assertEqual(first_detection_formula, enumerated_lead)
            for candidate in range(5):
                enlarged = tuple(sorted(set(selected) | {candidate}))
                self.assertEqual(exact.grid_marginal(weights, rows, selected, candidate),
                                 values[enlarged] - values[selected])
        for small in subsets:
            for large in subsets:
                if set(small) <= set(large):
                    for j in set(range(5)) - set(large):
                        gain_small = values[tuple(sorted((*small, j)))] - values[small]
                        gain_large = values[tuple(sorted((*large, j)))] - values[large]
                        self.assertGreaterEqual(gain_small, gain_large)
        for budget in range(6):
            chosen, optimum = exact.grid_optimum(weights, rows, budget)
            self.assertEqual(optimum, max(value for s, value in values.items()
                                          if len(s) == budget))
            self.assertEqual(values[tuple(sorted(chosen))], optimum)
            greedy, achieved = exact.grid_greedy(weights, rows, budget)
            self.assertEqual(values[tuple(sorted(greedy))], achieved)
            if budget:
                self.assertGreaterEqual(achieved,
                                        (1 - F(budget - 1, budget)**budget) * optimum)

    def test_t1_seeded_scenario_and_individual_test_simulation(self):
        weights = (F(1, 4), F(3, 4))
        rows = ((F(0), F(1, 2), F(3, 4)), (F(1, 3), F(1, 5), F(1, 2)))
        rng, draws, successes, gains, lead_sum = random.Random(40117), 18000, 0, 0, 0
        for _ in range(draws):
            row = rows[0 if rng.random() < .25 else 1]
            outcomes = [rng.random() < float(p) for p in row]
            prior = outcomes[0] or outcomes[1]
            successes += prior
            gains += (not prior) and outcomes[2]
            lead_sum += max((4 - j for j, positive in enumerate(outcomes) if positive), default=0)
        probability = float(exact.grid_value(weights, rows, (0, 1)))
        marginal = float(exact.grid_marginal(weights, rows, (0, 1), 2))
        self.assert_bounded_mc(successes / draws, probability, draws, 1,
                               "T1 unconditional detection")
        self.assert_bounded_mc(gains / draws, marginal, draws, 1, "T1 one additional test")
        lead_reference = sum(
            p * max((4 - j for j in range(3) if mask & (1 << j)), default=0)
            for mask, p in enumerate_outcomes(weights, rows).items()
        )
        self.assert_bounded_mc(lead_sum / draws, lead_reference, draws, 4,
                               "T1 noisy first-detection lead")

    def test_independent_normal_and_general_gamma_numeric_references(self):
        with localcontext() as context:
            context.prec = 70
            self.assertEqual(exact.normal_cdf(0), D("0.5"))
            known = D("0.84134474606854294858523254563203792247791296672660")
            self.assertLess(abs(exact.normal_cdf(1) - known), D("1e-49"))
            for value in ("0.1", "0.7", "2.3", "7.1"):
                self.assertLess(abs(exact.normal_cdf(value)
                                    + exact.normal_cdf(-D(value)) - 1), D("1e-65"))
        for shape in (.5, .8, 1.7, 3.5, 6.4):
            for value in (.03, .5, 2.0, 9.0):
                with self.subTest(shape=shape, value=value):
                    reference = gamma_cdf_independent(shape, value)
                    self.assertAlmostEqual(float(exact.gamma_cdf(shape, value)),
                                           reference, delta=2e-10)
        for shape in (1, 2, 5):
            for value in (.01, 1, 7, 40):
                self.assertAlmostEqual(float(exact.gamma_cdf(shape, value)),
                                       gamma_integer_cdf(shape, 1)(value), delta=3e-15)
        for value in (.000001, .1, 2, 20):
            reference = math.erf(math.sqrt(value))
            self.assertAlmostEqual(float(exact.gamma_cdf(.5, value)), reference, delta=3e-15)
            reference -= 2 * math.sqrt(value / math.pi) * math.exp(-value)
            self.assertAlmostEqual(float(exact.gamma_cdf(1.5, value)), reference, delta=3e-15)
            reference -= 4 * value**1.5 * math.exp(-value) / (3 * math.sqrt(math.pi))
            self.assertAlmostEqual(float(exact.gamma_cdf(2.5, value)), reference, delta=3e-15)

    def test_t2_direct_growth_event_integration_and_overlapping_intervals(self):
        cdf = lognormal_cdf(.2, .8)
        decimal_cdf = lambda r: exact.normal_cdf((r.ln() - D(".2")) / D(".8"))
        for times in ((), (.1,), (.1, .2, .4), (.7, .8, 1.0),
                      (.1, 1.2, 1.3, 7.0), (.2, .2, .4)):
            for endpoint in (None, .4, 1.5, 9.0):
                expected = independent_success(times, 1.0, 2.0, cdf, endpoint)
                actual = exact.ideal_value(times, 1, 2, decimal_cdf, endpoint)
                self.assertAlmostEqual(float(actual), expected, delta=3e-15)

    def test_t2_global_optima_against_independent_competitor_search(self):
        cases = []
        for n, horizon, endpoint in ((1, 4.0, None), (3, 8.0, None),
                                     (3, .6, None), (3, 2.0, 2.0)):
            cases.append((exact.optimal_lognormal(n, 1, 2, horizon, .1, .8, endpoint),
                          n, horizon, endpoint, lognormal_cdf(.1, .8)))
            cases.append((exact.optimal_gamma(n, 1, 2, horizon, 2, .6, endpoint),
                          n, horizon, endpoint, gamma_integer_cdf(2, .6)))
        for solution, n, horizon, endpoint, cdf in cases:
            with self.subTest(n=n, horizon=horizon, endpoint=endpoint,
                              lower=str(solution["lower_rate"])):
                times = tuple(map(float, solution["schedule"]))
                self.assertEqual(len(times), n)
                self.assertTrue(all(0 < t <= horizon * (1 + 1e-14) for t in times))
                for left, right in zip(times, times[1:]):
                    self.assertAlmostEqual(right / left, 2.0, places=13)
                expected = independent_success(times, 1, 2, cdf, endpoint)
                self.assertAlmostEqual(float(solution["value"]), expected, delta=3e-14)
                denominator = 1.0 if endpoint is None else 1 - cdf(2 / endpoint)
                self.assertAlmostEqual(float(solution["denominator"]), denominator,
                                       delta=3e-14)
                self.assertAlmostEqual(float(solution["conditional_value"]),
                                       expected / denominator, delta=3e-13)
                candidates = sorted({horizon * math.exp(-4 + j / 3)
                                     for j in range(13)} | set(times))
                for competitor in combinations(candidates, n):
                    self.assertLessEqual(independent_success(competitor, 1, 2, cdf, endpoint),
                                         expected + 5e-14)

    def test_t2_nonbinary_ratios_near_one_and_initial_wait_boundary(self):
        # Cover q<2, q>2, and a narrow detectable-to-clinical window. A ratio
        # accidentally hardcoded to two would pass the main reference fixtures.
        for a, b in ((2, 3), (1, 3), (1, 1.001)):
            for solution, cdf in (
                (exact.optimal_lognormal(3, a, b, 20, .1, .8), lognormal_cdf(.1, .8)),
                (exact.optimal_gamma(3, a, b, 20, 2, .6), gamma_integer_cdf(2, .6)),
            ):
                times = tuple(map(float, solution["schedule"]))
                q = b / a
                self.assertAlmostEqual(times[1] / times[0], q, delta=3e-14)
                self.assertAlmostEqual(times[2] / times[1], q, delta=3e-14)
                expected = independent_success(times, a, b, cdf)
                self.assertAlmostEqual(float(solution["value"]), expected, delta=3e-14)
                arithmetic = (times[0], (times[0] + times[-1]) / 2, times[-1])
                self.assertGreater(expected, independent_success(arithmetic, a, b, cdf) + 1e-14)
                if q < 2:
                    self.assertLess(times[1] - times[0], times[0])
                candidates = sorted({20 * math.exp(-5 + j / 2) for j in range(11)} | set(times))
                for competitor in combinations(candidates, 3):
                    self.assertLessEqual(independent_success(competitor, a, b, cdf),
                                         expected + 5e-14)

    def test_t2_horizon_estimands_have_different_boundary_schedules(self):
        eventual = exact.optimal_lognormal(3, 1, 2, 1, 0, 1)
        within = exact.optimal_lognormal(3, 1, 2, 1, 0, 1, endpoint=1)
        self.assertEqual(eventual["schedule"][-1], 1)
        self.assertEqual(within["schedule"][-1], D(".5"))
        self.assertEqual(eventual["denominator"], 1)
        self.assertLess(within["denominator"], 1)
        self.assertLess(within["value"], within["conditional_value"])
        self.assertEqual(exact.ideal_value((1,), 1, 2,
                                          lambda r: exact.normal_cdf(r.ln()), endpoint=1), 0)

    def test_t2_t3_seeded_lognormal_paths_lead_and_additional_test_value(self):
        mu, sigma, a, b, horizon = 0.0, .8, 1.0, 2.0, 10.0
        solutions = {"three": exact.optimal_lognormal(3, a, b, horizon, mu, sigma),
                     "four": exact.optimal_lognormal(4, a, b, horizon, mu, sigma)}
        schedules = {name: tuple(map(float, item["schedule"]))
                     for name, item in solutions.items()}
        schedules["equally_spaced"] = (horizon / 3, 2 * horizon / 3, horizon)
        observations = simulation.simulate_growth(
            schedules, 22000, 48121, lambda rng: rng.lognormvariate(mu, sigma), a=a, b=b
        )
        cdf = lognormal_cdf(mu, sigma)
        for name, times in schedules.items():
            probability = independent_success(times, a, b, cdf)
            self.assert_probability_mc(observations[name], probability, "T2 " + name)
            lead = independent_lognormal_lead(times, a, b, mu, sigma)
            # A successful draw at t has R>=a/t, hence P-t<=((b/a)-1)H.
            self.assert_bounded_mc(observations[name]["mean_lead"], lead, 22000,
                                   (b / a - 1) * horizon, "first-detection lead")
            loss = (b - a) * math.exp(-mu + sigma**2 / 2) - lead
            second_moment_bound = (b - a)**2 * math.exp(-2 * mu + 2 * sigma**2)
            self.assert_second_moment_mc(observations[name]["mean_loss"], loss,
                                        second_moment_bound, 22000, "lead-time loss")
        paired = observations["four"]["paired"]["three"]
        delta = float(solutions["four"]["value"] - solutions["three"]["value"])
        self.assert_bounded_mc(paired["probability_difference"], delta, 22000, 2,
                               "T3 reoptimized delta")
        self.assertGreater(float(solutions["three"]["value"]),
                           independent_success(schedules["equally_spaced"], a, b, cdf))

    def test_t2_t3_seeded_gamma_and_finite_endpoint(self):
        solutions = {str(n): exact.optimal_gamma(n, 1, 2, 12, "2.5", ".4")
                     for n in (2, 3)}
        schedules = {name: tuple(map(float, item["schedule"]))
                     for name, item in solutions.items()}
        results = simulation.simulate_growth(schedules, 18000, 55219,
                                             lambda rng: rng.gammavariate(2.5, .4), a=1, b=2)
        independent_cdf = lambda r: gamma_cdf_independent(2.5, r / .4)
        for name, times in schedules.items():
            self.assert_probability_mc(results[name],
                                       independent_success(times, 1, 2, independent_cdf),
                                       "T2 gamma " + name)
        paired = results["3"]["paired"]["2"]
        delta = float(solutions["3"]["value"] - solutions["2"]["value"])
        self.assert_bounded_mc(paired["probability_difference"], delta, 18000, 2,
                               "T3 gamma delta")
        bounded = exact.optimal_lognormal(3, 1, 2, 2, 0, .8, endpoint=2)
        observed = simulation.simulate_growth(
            {"bounded": tuple(map(float, bounded["schedule"]))}, 18000, 65009,
            lambda rng: rng.lognormvariate(0, .8), a=1, b=2, endpoint=2
        )["bounded"]
        self.assert_probability_mc(observed, bounded["value"], "T2 endpoint unconditional")
        denominator = float(bounded["denominator"])
        self.assert_bounded_mc(observed["scored"] / observed["draws"], denominator,
                               observed["draws"], 1, "endpoint denominator")

    def test_t3_heterogeneity_parameters_do_not_determine_front_loading(self):
        normal = [exact.optimal_lognormal(3, 1, 2, 100, ".2", sigma)
                  for sigma in (".3", ".8", "1.4")]
        self.assertEqual(normal[0]["schedule"], normal[1]["schedule"])
        self.assertEqual(normal[1]["schedule"], normal[2]["schedule"])
        self.assertGreater(normal[0]["value"], normal[1]["value"])
        self.assertGreater(normal[1]["value"], normal[2]["value"])
        gamma = [exact.optimal_gamma(3, 1, 2, 100, shape, F(1, 1) / F(shape))
                 for shape in (F(1, 2), F(1), F(2), F(5))]
        with localcontext() as context:
            context.prec = 70
            for left, right in zip(gamma, gamma[1:]):
                for x, y in zip(left["schedule"], right["schedule"]):
                    self.assertLess(abs(x - y), D("1e-60"))
                self.assertLess(left["value"], right["value"])
        # Fixed arithmetic mean: more lognormal dispersion shifts every optimum
        # later by exp((sigma_2**2 - sigma_1**2)/2) when the horizon is inactive.
        first = exact.optimal_lognormal(3, 1, 2, 100, "-.08", ".4")
        second = exact.optimal_lognormal(3, 1, 2, 100, "-.5", 1)
        for x, y in zip(first["schedule"], second["schedule"]):
            self.assertAlmostEqual(float(y / x), math.exp(.42), places=13)

    def test_t3_exact_extra_test_formula_and_diminishing_optimal_values(self):
        with localcontext() as context:
            context.prec = 70
            for horizon in (".4", "2", "100"):
                for factory in (
                    lambda n: exact.optimal_lognormal(n, 1, 2, horizon, 0, ".9"),
                    lambda n: exact.optimal_gamma(n, 1, 2, horizon, "1.7", ".7"),
                ):
                    values = [D(0)] + [factory(n)["value"] for n in range(1, 8)]
                    differences = [right - left for left, right in zip(values, values[1:])]
                    self.assertTrue(all(delta >= 0 for delta in differences))
                    for first, second in zip(differences, differences[1:]):
                        self.assertLessEqual(second, first + D("1e-60"))
            width = math.log(2)
            previous = 0.0
            for n in range(1, 7):
                value = float(exact.optimal_lognormal(n, 1, 2, 1000, 0, 1)["value"])
                reference_delta = (math.erf(n * width / (2 * math.sqrt(2)))
                                   - math.erf((n - 1) * width / (2 * math.sqrt(2))))
                self.assertAlmostEqual(value - previous, reference_delta, delta=4e-15)
                previous = value

    def test_t3_additional_test_value_has_an_interior_heterogeneity_peak(self):
        for budget in (1, 3, 6):
            width = math.log(2)
            peak = width * math.sqrt((2 * budget + 1)
                                     / (8 * math.log1p(1 / budget)))
            increments = []
            for sigma in (peak / 2, peak, peak * 2):
                before = exact.optimal_lognormal(budget, 1, 2, 1000, 0, sigma)
                after = exact.optimal_lognormal(budget + 1, 1, 2, 1000, 0, sigma)
                increments.append(float(after["value"] - before["value"]))
            self.assertGreater(increments[1], increments[0])
            self.assertGreater(increments[1], increments[2])
            # Independent differentiation of the normal CDF expression checks
            # the critical-point equation, not only neighboring objective values.
            lower, upper = budget * width / 2, (budget + 1) * width / 2
            derivative_numerator = (lower * math.exp(-lower**2 / (2 * peak**2))
                                    - upper * math.exp(-upper**2 / (2 * peak**2)))
            self.assertAlmostEqual(derivative_numerator, 0, delta=4e-15)

    def test_t4_bound_terms_limits_and_finite_endpoint(self):
        for finite in (False, True):
            result = exact.robustness_bound(3, 10, 20, ".01", ".5", 800, ".08", ".05", finite)
            shift = math.log(.13 * .92 / (.08 * .87))
            burden = .5 * (3 * (math.log(10 / 9.99) + math.log(20 / 19.99))
                             + (math.log(20 / 19.99) if finite else 0))
            noise = 3 * .5 * math.log1p(shift / 9.99) + math.exp(-2 * 800 * .05**2)
            self.assertAlmostEqual(float(result["burden_error"]), burden, delta=3e-16)
            self.assertAlmostEqual(float(result["noise_error"]), noise, delta=3e-16)
            self.assertAlmostEqual(float(result["regret_bound"]), 2 * burden + noise, delta=6e-16)
            self.assertLess(result["regret_bound"], D(".11"))
        errors = []
        for m in (10000, 1000000, 100000000):
            delta = m ** (-1 / 3)
            result = exact.robustness_bound(3, 1, 2, 0, ".5", m, ".1", delta)
            self.assertEqual(result["burden_error"], 0)
            errors.append(result["noise_error"])
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        with self.assertRaises(ValueError):
            exact.robustness_bound(3, 1, 2, 1, ".5", 100, ".1", ".01")

    def test_t4_seeded_random_burden_and_individual_molecule_noise(self):
        optimum = exact.optimal_lognormal(3, 10, 20, 50, 0, ".8", endpoint=30)
        schedules = {"geometric": tuple(map(float, optimum["schedule"])),
                     "competitor": (4.0, 8.0, 14.0)}
        parameters = dict(a=10, b=20, endpoint=30,
                          burden_sampler=lambda rng: rng.uniform(-.01, .01))
        perfect = simulation.simulate_growth(schedules, 14000, 77521,
                                              lambda rng: rng.lognormvariate(0, .8), **parameters)
        noisy = simulation.simulate_growth(schedules, 14000, 77521,
                                            lambda rng: rng.lognormvariate(0, .8),
                                            molecule_count=800, lod=.08, **parameters)
        bound = exact.robustness_bound(3, 10, 20, ".01", ".5", 800, ".08", ".05", True)
        p0 = float(optimum["value"])
        tolerance = self.bounded_allowance(14000, 1)
        self.assertLessEqual(abs(perfect["geometric"]["probability"] - p0),
                             float(bound["burden_error"]) + tolerance)
        for name in schedules:
            # Same seed and separated patient/assay streams yield the same
            # patients. Noise can only remove successes under the hard gate.
            self.assertEqual(noisy[name]["scored"], perfect[name]["scored"])
            self.assertLessEqual(noisy[name]["successes"], perfect[name]["successes"])
            self.assertLessEqual(perfect[name]["probability"] - noisy[name]["probability"],
                                 float(bound["noise_error"]) + tolerance)
        paired = noisy["competitor"]["paired"]["geometric"]
        self.assertLessEqual(paired["probability_difference"], float(bound["regret_bound"])
                             + self.bounded_allowance(14000, 2))

    def test_falsification_geometric_exactness_fails_with_noisy_unimodal_rates(self):
        ideal = exact.optimal_lognormal(2, 1, 2, 2, 0, ".001")
        schedules = {"geometric": tuple(map(float, ideal["schedule"])),
                     "overlap": (1.5, 1.6)}
        observed = simulation.simulate_growth(schedules, 12000, 89231,
                                              lambda rng: rng.lognormvariate(0, .001),
                                              a=1, b=2, molecule_count=1, lod=.5)
        # Each geometric window hits at most once. All eligible burden excesses
        # are <1, so success<=logistic(1)<11/15. For overlap, |log R|<=.01 has
        # probability>=.99 by Chebyshev and both tests succeed with p>=1/2.
        self.assertLess(math.e / (1 + math.e), 11 / 15)
        self.assertGreater(.99 * .75, 11 / 15)
        self.assertLess(observed["geometric"]["probability"],
                        11 / 15 + self.bounded_allowance(12000, 1))
        self.assertGreater(observed["overlap"]["probability"],
                           .99 * .75 - self.bounded_allowance(12000, 1))
        paired = observed["overlap"]["paired"]["geometric"]
        self.assertGreater(paired["probability_difference"], self.bounded_allowance(12000, 2))

    def test_falsification_unimodality_is_needed_for_a_single_geometric_block(self):
        # Synthetic rate atoms 1 and 100. Each ratio-two schedule spans a total
        # rate ratio of four and cannot cover both atoms. Two separated tests do.
        observed = simulation.simulate_growth(
            {"separated": (.015, 1.5), "one_block": (.75, 1.5)}, 4000, 90017,
            lambda rng: 1.0 if rng.getrandbits(1) else 100.0, a=1, b=2
        )
        self.assertEqual(observed["separated"]["successes"], 4000)
        self.assert_probability_mc(observed["one_block"], .5, "bimodal counterexample")

    def test_t5_equal_spacing_exact_loss_and_additional_test_value(self):
        horizon, window = F(12), F(15)
        for n in range(1, 5):
            optimum = tuple(i * horizon / n for i in range(1, n + 1))
            expected = horizon / (2 * n)
            self.assertEqual(uniform_loss_by_cells(optimum, horizon, window), expected)
            self.assertEqual(exact.uniform_onset_loss(optimum, horizon, window), expected)
            next_times = tuple(i * horizon / (n + 1) for i in range(1, n + 2))
            improvement = expected - uniform_loss_by_cells(next_times, horizon, window)
            self.assertEqual(improvement, horizon / (2 * n * (n + 1)))
            grid = sorted(set(F(i) for i in range(1, 13)) | set(optimum))
            for times in combinations(grid, n):
                independent = uniform_loss_by_cells(times, horizon, window)
                self.assertEqual(exact.uniform_onset_loss(times, horizon, window), independent)
                self.assertGreaterEqual(independent, expected)
        self.assertEqual(exact.uniform_onset_loss((), horizon, window), window)
        self.assertEqual(exact.uniform_onset_loss((horizon,), horizon, horizon), horizon / 2)

    def test_t5_seeded_uniform_onset_paths(self):
        schedules = {"equal": (2.5, 5.0, 7.5, 10.0), "irregular": (1.0, 4.0, 8.0, 9.0),
                     "none": ()}
        observed = simulation.sample_uniform_onsets(schedules, 18000, 99173, H=10, W=12)
        for name, times in schedules.items():
            expected = uniform_loss_by_cells(times, F(10), F(12))
            self.assert_bounded_mc(observed[name]["mean_loss"], expected, 18000, 12, "T5 " + name)
        self.assertEqual(observed["none"]["mean_loss"], 12)
        self.assertLess(observed["equal"]["mean_loss"], observed["irregular"]["mean_loss"])

    def test_t6_sharp_mean_support_bounds_exact_finite_distributions(self):
        for leads in ((F(0), F(1), F(3), F(9)), (F(0), F(0), F(4), F(8)),
                      (F(2), F(3), F(7)), (F(5),), (F(0),)):
            for step in (F(1), F(3), F(12)):
                mean = sum(leads) / len(leads)
                retained = sum((lead + max(F(0), lead - step)) / 2
                               for lead in leads) / len(leads)
                lower, upper = exact.thinning_bounds(mean, min(leads), max(leads), step)
                self.assertLessEqual(lower, retained)
                self.assertLessEqual(retained, upper)
        mean, low, high, step = F(87, 10), F(4, 5), F(33, 2), F(3)
        lower, upper = exact.thinning_bounds(mean, low, high, step)
        self.assertEqual(lower, F(36, 5))
        self.assertEqual(upper, F(6081, 785))
        # Exact extremizers of the relaxed moment/support problem. These are
        # synthetic distributions, not reconstructed cohort observations.
        self.assertEqual((mean + max(0, mean - step)) / 2, lower)
        high_weight = (mean - low) / (high - low)
        endpoint_value = ((1 - high_weight) * (low + max(0, low - step)) / 2
                          + high_weight * (high + max(0, high - step)) / 2)
        self.assertEqual(endpoint_value, upper)
        rounding_lower = exact.thinning_bounds(F(173, 20), F(3, 4), F(331, 20), 3)[0]
        rounding_upper = exact.thinning_bounds(F(35, 4), F(3, 4), F(331, 20), 3)[1]
        self.assertEqual(rounding_lower, F(143, 20))
        self.assertEqual(rounding_upper, F(4933, 632))

    def test_t6_seeded_thinning_and_zero_lead_misses_keep_denominator(self):
        for leads in ((8.7,), (0.0, 1.0, 3.0, 9.0), tuple([.8] * 78 + [16.5] * 79)):
            observed = simulation.thinning_simulation(leads, 3, 20000, 102017)
            expected = sum((lead + max(0, lead - 3)) / 2 for lead in leads) / len(leads)
            reduction = sum(min(lead, 3) / 2 for lead in leads) / len(leads)
            self.assert_bounded_mc(observed["thinned_mean"], expected, 20000, max(leads),
                                   "T6 retained-cohort mean")
            self.assert_bounded_mc(observed["mean_reduction"], reduction, 20000, 3, "T6 reduction")
        zeros = simulation.thinning_simulation((0,), 3, 100, 102019)
        self.assertEqual(zeros["successes"], 0)
        self.assertEqual(zeros["thinned_mean"], 0)

    def test_invalid_inputs_are_not_silently_given_theorem_interpretations(self):
        for arguments in ((0, 1, 2, 10, 0, 1), (2, 1, 1, 10, 0, 1),
                          (2, 1, 2, 10, 0, 0), (2, 1, 2, 0, 0, 1)):
            with self.assertRaises(ValueError):
                exact.optimal_lognormal(*arguments)
        with self.assertRaises(ValueError):
            exact.optimal_gamma(2, 1, 2, 10, 0, 1)
        with self.assertRaises(ValueError):
            exact.uniform_onset_loss((1,), 3, 2)
        with self.assertRaises(ValueError):
            exact.grid_value((F(1),), ((F(1, 2),),), (0, 0))
        with self.assertRaises(ValueError):
            exact.thinning_bounds(2, 3, 4, 1)


if __name__ == "__main__":
    unittest.main()
