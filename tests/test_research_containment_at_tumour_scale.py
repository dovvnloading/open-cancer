"""Independent checks of synthetic, nonclinical mathematics; no treatment advice.

Only the standard library is used. Exact certificates are compared with an
independently integrated branching backward equation, and the population
simulator generates its own birth/death/mutation events. Fixed-seed checks
test finite-N witnesses and scaling slopes, not equality of finite-N means
with their asymptotic coefficients. The workload is fixed and bounded.
"""

from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction as F
import importlib.util
import math
from pathlib import Path
import random
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / "research" / "containment-at-tumour-scale"


def load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, NOTE / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


V = load_module("tumour_scale_verifier", "verifier.py")
S = load_module("tumour_scale_simulation", "simulation.py")


def decimal_fraction(x):
    return Decimal(x.numerator) / Decimal(x.denominator)


def backward_branching(s0=1.0, step=0.05, finish=800.0):
    """Independent numerical route: backward Riccati + integrated arrivals.

    Hard-coded default rates intentionally do not call the verifier's series,
    quadrature, parameter object, or survival function. At the remote terminal
    time competition is <5e-12; the terminal survival probability tends to
    (birth-death)/birth. A step-halving check controls numerical evidence here;
    this floating calculation is not described as an exact certificate.
    """
    birth, released, decline, competition = 0.8, 0.2, 0.03, 0.12
    t, survival, integral = finish, released / birth, 0.0
    count = math.ceil(finish / step)
    dt = -finish / count

    def derivative(time, p):
        sensitive = s0 * math.exp(-decline * time)
        return (-(released - competition * sensitive) * p + birth * p * p,
                -sensitive * p)

    for _ in range(count):
        a, ai = derivative(t, survival)
        b, bi = derivative(t + dt / 2, survival + dt * a / 2)
        c, ci = derivative(t + dt / 2, survival + dt * b / 2)
        d, di = derivative(t + dt, survival + dt * c)
        survival += dt * (a + 2 * b + 2 * c + d) / 6
        integral += dt * (ai + 2 * bi + 2 * ci + di) / 6
        t += dt
    # The omitted mutation integral is between zero and this analytic bound.
    tail = s0 * math.exp(-decline * finish) * released / (birth * decline)
    return survival, integral, tail


class ExactMathematics(unittest.TestCase):
    def test_assumptions_and_strict_boundaries(self):
        p = V.DEFAULT
        V.validate(p)
        self.assertEqual((p.delta + p.a_sr * p.B) * p.L, F(108, 125))
        for bad in (
            replace(p, alpha=-1),
            replace(p, k_max=F(1, 3)),
            replace(p, L=F(5)),  # L=1/g0 is excluded.
            replace(p, alpha=p.g0 - p.beta * p.epsilon),
            replace(p, alpha=p.g0 - p.beta * p.epsilon / 2),
            replace(p, a_sr=F(1, 4)),
        ):
            with self.assertRaises(ValueError):
                V.validate(bad)
        with self.assertRaises(TypeError):
            V.RationalParameters(alpha=0.12)
        with self.assertRaises(ValueError):
            V.coupling_bound(100, 100)
        with self.assertRaises(ValueError):
            V.supply_interval(panels=15)

    def test_elementary_intervals_against_decimal(self):
        with localcontext() as context:
            context.prec = 80
            for value in (F(1, 10**40), F(1, 10), F(1), F(42, 25),
                          F(2), F(10**11), F(10**40)):
                lo, hi = V.log_interval(value)
                actual = decimal_fraction(value).ln()
                self.assertLessEqual(decimal_fraction(lo), actual)
                self.assertGreaterEqual(decimal_fraction(hi), actual)
                self.assertLess(hi - lo, F(1, 10**23))
            for value in (F(0), F(1, 1000), F(1), F(627, 100), F(100)):
                lo, hi = V.exp_neg_interval(value)
                actual = (-decimal_fraction(value)).exp()
                self.assertLessEqual(decimal_fraction(lo), actual)
                self.assertGreaterEqual(decimal_fraction(hi), actual)
                self.assertLess(hi - lo, F(1, 10**23))

    def test_survival_and_supply_against_backward_equation(self):
        coarse_p, coarse_a, _ = backward_branching(step=0.1)
        fine_p, fine_a, tail = backward_branching(step=0.05)
        self.assertLess(abs(coarse_p - fine_p), 2e-9)
        self.assertLess(abs(coarse_a - fine_a), 2e-8)
        plo, phi = V.survival_interval(1)
        self.assertAlmostEqual(float((plo + phi) / 2), fine_p, delta=2e-9)
        alo, ahi = V.supply_interval()
        self.assertLessEqual(float(alo), fine_a + tail)
        self.assertGreaterEqual(float(ahi), fine_a)
        self.assertLess(ahi - alo, F(1, 10**6))
        # Independent backward solve at other founder-time densities.
        for density in (F(0), F(1, 4), F(2)):
            expected, _, _ = backward_branching(float(density), step=0.1)
            lo, hi = V.survival_interval(density)
            self.assertAlmostEqual(float((lo + hi) / 2), expected, delta=3e-8)

    def test_exact_zero_competition_special_case(self):
        p = replace(V.DEFAULT, alpha=F(0))
        self.assertEqual(V.survival_interval(17, p), (F(1, 4), F(1, 4)))
        self.assertEqual(V.supply_interval(p), (F(25, 3), F(25, 3)))
        self.assertIsNone(V.phase_interval(p))
        for supply in (F(1, 1000), F(1, 10), F(1), F(4)):
            values = V.asymptotic_mean_intervals(supply, p)
            self.assertGreater(values["maximum"][0], values["containment"][1])

    def test_phase_enclosure_sign_and_uniqueness(self):
        low, high = V.phase_interval()
        self.assertGreater(low, F("0.08273970"))
        self.assertLess(high, F("0.08273972"))
        self.assertLess(high - low, F(1, 10**7))
        below = V.asymptotic_mean_intervals(low - F(1, 10**5))
        above = V.asymptotic_mean_intervals(high + F(1, 10**5))
        self.assertGreater(below["maximum"][0], below["containment"][1])
        self.assertGreater(above["containment"][0], above["maximum"][1])
        # Independent algebraic signs: delta(lambda) has positive derivative
        # (J-d)*A*exp(-lambda*A) and the limits c-J<0<c-d.
        d, c, J = F(5), F(25, 2), F(88, 5)
        self.assertLess(c - J, 0)
        self.assertGreater(c - d, 0)
        last = None
        for supply in (F(1, 100), F(1, 20), F(1, 10), F(1, 2), F(1)):
            values = V.asymptotic_mean_intervals(supply)
            gap = (c - values["maximum"][1], c - values["maximum"][0])
            if last is not None:
                self.assertGreater(gap[0], last[1])
            last = gap

    def test_zero_supply_critical_cost_and_strong_competition(self):
        zero = V.asymptotic_mean_intervals(0)
        self.assertEqual(zero, {"containment": (F(16), F(16)),
                                "maximum": (F(88, 5), F(88, 5))})
        positive = V.asymptotic_mean_intervals(F(1, 10**8))
        self.assertEqual(positive["containment"], (F(25, 2), F(25, 2)))
        for death in (F(4, 5), F(9, 10)):
            p = replace(V.DEFAULT, d_r=death)
            self.assertIsNone(V.phase_interval(p))
            for supply in (F(0), F(1, 100), F(2)):
                self.assertEqual(V.asymptotic_mean_intervals(supply, p), zero)
        for competition in (F(1, 5), F(1, 4)):
            p = replace(V.DEFAULT, alpha=competition)
            self.assertEqual(V.asymptotic_mean_intervals(1, p)["containment"],
                             (F(16), F(16)))
            self.assertIsNotNone(V.phase_interval(p))

    def test_shared_hill_law_and_steepness_invariance(self):
        for exponent in (1, 2, 4, 8):
            p = replace(V.DEFAULT, hill_h=exponent)
            V.validate(p)
            self.assertEqual(V.hill_kill(0, p), 0)
            self.assertEqual(V.hill_kill(1, p), F(33, 100))
            # Same realized maximum kill => same demographic limit.
            self.assertEqual(V.phase_interval(p), V.phase_interval())
            floating = replace(S.Parameters(), hill=float(exponent))
            for kill in (0.0, 0.12, 0.3, 0.33):
                dose = S.inverse_hill(kill, floating)
                self.assertGreaterEqual(dose, 0.0)
                self.assertLessEqual(dose, 1.0)
                self.assertAlmostEqual(S.hill_kill(dose, floating), kill, places=13)
        # Away from u_max=EC50, steepness changes realized maximum kill and
        # therefore genuinely enters the computed comparison, when feasible.
        boundaries = []
        for exponent in (1, 4):
            p = replace(V.DEFAULT, hill_h=exponent, hill_ec50=F(101, 100))
            p = replace(p, k_max=V.hill_kill(p.u_max, p))
            V.validate(p)
            boundaries.append(V.phase_interval(p))
        self.assertTrue(boundaries[0][1] < boundaries[1][0]
                        or boundaries[1][1] < boundaries[0][0])

    def test_quantitative_branching_coupling(self):
        previous = None
        for n in (10**8, 10**9, 10**10, 10**11):
            bound = V.coupling_bound(n, F(1, 10))
            m = bound["m"]
            self.assertGreaterEqual(m**4, n)
            self.assertLess((m - 1)**4, n)
            self.assertLess(bound["maximum"], F(1, 1000))
            self.assertGreater(bound["maximum"], bound["containment"])
            if previous is not None:
                self.assertLess(bound["maximum"], previous)
            previous = bound["maximum"]


class IndependentPopulationSimulation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Same CLI seed convention as the README's four-scale experiment.
        # Two endpoint sizes keep this test comfortably below 60 seconds.
        cls.results = {}
        for population, offset in ((10**8, 0), (10**11, 12)):
            for supply, supply_offset in ((0.005, 0), (1.0, 2)):
                for policy, policy_offset in (("C", 0), ("M", 1)):
                    cls.results[population, supply, policy] = S.simulate(
                        policy, population, supply, 400,
                        20261117 + offset + supply_offset + policy_offset)

    def test_tumour_scale_ranking_witnesses(self):
        for population in (10**8, 10**11):
            for supply, sign in ((0.005, -1), (1.0, 1)):
                c = self.results[population, supply, "C"]
                m = self.results[population, supply, "M"]
                difference = c["mean_time"] - m["mean_time"]
                se = math.hypot(c["se_time"], m["se_time"])
                self.assertGreater(sign * difference, 6 * se)
                limits = V.asymptotic_mean_intervals(F(str(supply)))
                if sign > 0:
                    self.assertGreater(limits["containment"][0], limits["maximum"][1])
                else:
                    self.assertGreater(limits["maximum"][0], limits["containment"][1])

    def test_admissibility_and_pathwise_bound(self):
        p = S.Parameters()
        for (population, supply, _), result in self.results.items():
            horizon = p.course_scale * math.log(population)
            r_sensitive = (1 - supply / population) * p.b_s - p.d_s
            bound = horizon + (math.log(p.threshold / p.s0)
                     + (p.k_max - r_sensitive + p.a_sr * p.threshold) * horizon
                     ) / (r_sensitive - p.a_sr * p.threshold)
            self.assertLess(result["mean_time"], bound)
            self.assertLessEqual(result["max_dose"], p.u_max + 1e-12)
            self.assertLessEqual(result["max_exposure"], p.u_max * horizon + 1e-8)
            self.assertGreater(result["min_sensitive_count"], 1000)

    def test_establishment_weight_from_actual_events(self):
        alo, ahi = V.supply_interval()
        a = float((alo + ahi) / 2)
        for population in (10**8, 10**11):
            for supply in (0.005, 1.0):
                result = self.results[population, supply, "M"]
                expected = -math.expm1(-supply * a)
                # A theoretical Bernoulli variance avoids zero-SE assertions
                # when all sampled trajectories switch. The 0.003 tolerance
                # also allows finite switching-time/tail effects; no exact
                # finite-N identity with establishment probability is claimed.
                se = math.sqrt(expected * (1 - expected) / result["samples"])
                coupling = float(V.coupling_bound(population, F(str(supply)))["maximum"])
                self.assertLess(abs(result["switch_probability"] - expected),
                                6 * se + coupling + 0.003)

    def test_growth_slopes_not_false_finite_size_equality(self):
        log_ratio = math.log(10**11 / 10**8)
        for policy, slope in (("C", 12.5), ("M", 5.0)):
            small = self.results[10**8, 1.0, policy]
            large = self.results[10**11, 1.0, policy]
            estimate = (large["mean_time"] - small["mean_time"]) / log_ratio
            se = math.hypot(small["se_time"], large["se_time"]) / log_ratio
            # M's actual limiting mean slope is ~5.024, not exactly 5,
            # because a small noescape tail remains at lambda=1.
            self.assertLess(abs(estimate - slope), 6 * se + 0.25)
        c_small = self.results[10**8, 1.0, "C"]["mean_time_over_log_n"]
        c_large = self.results[10**11, 1.0, "C"]["mean_time_over_log_n"]
        self.assertLess(c_small, c_large)
        self.assertLess(c_large, 12.5)
        # Startup at small positive supply remains material, as documented.
        self.assertGreater(self.results[10**8, 0.005, "C"]["mean_time_over_log_n"], 15)

    def test_zero_mutation_exact_rebound(self):
        p = S.Parameters()
        for population in (10**8, 10**9, 10**10, 10**11):
            horizon = p.course_scale * math.log(population)
            for policy in ("C", "M"):
                result = S.trial(policy, population, 0, random.Random(7))
                expected = horizon + math.log(p.threshold / p.s0) / p.r_s
                if policy == "M":
                    expected += (p.k_max - p.r_s) * horizon / p.r_s
                self.assertAlmostEqual(result.time, expected, places=9)
                self.assertFalse(result.switched)
                self.assertEqual(result.events, 0)

    def test_critical_resistance_and_strong_competition(self):
        population = 10**8
        base = S.Parameters()
        horizon = base.course_scale * math.log(population)
        for p in (replace(base, d_r=0.8), replace(base, d_r=0.9)):
            for policy in ("C", "M"):
                result = S.simulate(policy, population, 1.0, 12, 319, p)
                mutation_free = horizon + math.log(p.threshold / p.s0) / p.r_s
                if policy == "M":
                    mutation_free += (p.k_max - p.r_s) * horizon / p.r_s
                self.assertAlmostEqual(result["mean_time"], mutation_free, delta=0.01)
        strong = replace(base, alpha=0.25)
        result = S.simulate("C", population, 1.0, 12, 719, strong)
        self.assertAlmostEqual(result["mean_time"],
                               horizon + math.log(base.threshold / base.s0) / base.r_s,
                               delta=0.01)

    def test_event_hazard_and_ode_refinement(self):
        for growth in (-0.3, -0.01, 0.0, 0.05, 0.3):
            for constant, varying in ((0.0, 2.0), (0.7, 0.0), (2.0, 7.0)):
                for exponential in (0.01, 0.3, 2.0, 10.0):
                    cap = 12.0
                    wait = S._event_wait(exponential, constant, varying, growth, cap)
                    def integrated(t):
                        if growth == 0:
                            return (constant + varying) * t
                        return constant * t + varying * math.expm1(growth * t) / growth
                    if math.isfinite(wait):
                        self.assertLess(wait, cap + 1e-10)
                        self.assertAlmostEqual(integrated(wait), exponential, delta=1e-10)
                    else:
                        self.assertGreaterEqual(exponential + 1e-10, integrated(cap))
        for policy in ("C", "M"):
            for seed in range(5):
                coarse = S.trial(policy, 10**11, 1.0, random.Random(seed), ode_step=0.5)
                fine = S.trial(policy, 10**11, 1.0, random.Random(seed), ode_step=0.25)
                self.assertEqual(coarse.events, fine.events)
                self.assertEqual(coarse.switch_time, fine.switch_time)
                self.assertLess(abs(coarse.time - fine.time), 0.001)

    def test_simulation_independence(self):
        source = (NOTE / "simulation.py").read_text()
        self.assertNotIn("import verifier", source)
        self.assertNotIn("supply_interval(", source)
        self.assertNotIn("phase_interval(", source)
        self.assertNotIn("asymptotic_mean_intervals(", source)


if __name__ == "__main__":
    unittest.main()
