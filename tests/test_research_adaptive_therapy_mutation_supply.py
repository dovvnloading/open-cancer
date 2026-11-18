"""Independent checks of nonclinical finite-state mathematics; stdlib only."""

from fractions import Fraction as F
import importlib.util
from pathlib import Path
import sys
import unittest


DIRECTORY = (Path(__file__).resolve().parents[1] / "research"
             / "adaptive-therapy-mutation-supply")


def load_file(name):
    spec = importlib.util.spec_from_file_location(
        "adaptive_mutation_research_" + name, DIRECTORY / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


verifier = load_file("verifier")
simulation = load_file("simulation")


def independent_jump_bounds(mu, policy, tolerance=F(1, 10**25)):
    """Exact absorbing jump-chain sum, NOT the verifier's linear solve.

    Hand-enumerated states and rates from the model include an independent
    censor clock of rate 1.  The unabsorbed probability bounds the entire
    remaining series, so the result is an exact enclosure, not a tolerance
    comparison between two floating-point solvers.
    """
    # R1/R2 are the permanently sensitive-free states.  -1 is progression.
    if policy == "C":
        # State order C0,C1,R2,R1; cure/censor destinations are omitted.
        rows = [
            (2 * mu, {1: 2 * mu}),
            (mu + 4, {0: F(3), 2: mu + 1}),
            (F(8), {3: F(6), -1: F(2)}),
            (F(103), {2: F(100)}),
        ]
    else:
        # State order M20,M11,M10,R2,R1.
        rows = [
            (2 * mu + 4, {1: 2 * mu, 2: F(4)}),
            (mu + 6, {-1: F(1), 2: F(3), 3: mu, 4: F(2)}),
            (mu + 2, {4: mu}),
            (F(8), {-1: F(2), 4: F(6)}),
            (F(103), {3: F(100)}),
        ]
    mass = [F(1)] + [F(0)] * (len(rows) - 1)
    lower = F(0)
    for _ in range(2000):
        next_mass = [F(0)] * len(rows)
        for i, (outgoing, destinations) in enumerate(rows):
            for target, rate in destinations.items():
                weight = mass[i] * rate / (1 + outgoing)
                if target == -1:
                    lower += weight
                else:
                    next_mass[target] += weight
        remaining = sum(next_mass, F(0))
        if remaining < tolerance:
            return lower, lower + remaining
        mass = next_mass
    raise AssertionError("jump-chain enclosure did not converge")


class MutationSupplyMathematicsTests(unittest.TestCase):
    def test_general_formula_against_exact_occupancy_generator(self):
        report = verifier.verify()
        self.assertEqual(report["exact_generator_comparisons"], 96)

    def test_independent_exact_jump_chain_encloses_closed_forms(self):
        for mu in (F(0), F(1, 4), F(5, 2), F(25)):
            for policy in ("C", "M"):
                with self.subTest(mu=mu, policy=policy):
                    lower, upper = independent_jump_bounds(mu, policy)
                    closed = verifier.formula_laplace(mu, policy)
                    self.assertLessEqual(lower, closed)
                    self.assertLessEqual(closed, upper)
                    self.assertLess(upper - lower, F(1, 10**25))

    def test_polynomial_identity_and_complete_positive_root_count(self):
        polynomial = verifier.MAIN_POLYNOMIAL
        # Coefficientwise identity, valid for all mu, not a fitted polynomial.
        self.assertEqual(verifier.comparison_polynomial(),
                         tuple(-coefficient / 26 for coefficient in polynomial))
        self.assertEqual(verifier.sign_variations(polynomial), 2)
        self.assertGreater(verifier.poly_value(polynomial, F(0)), 0)
        self.assertLess(verifier.poly_value(polynomial, F(1)), 0)
        self.assertGreater(verifier.poly_value(polynomial, F(13)), 0)
        for left, right in ((F(0), F(1)), (F(1), F(13))):
            lo, hi = verifier.isolate_root(polynomial, left, right)
            self.assertLessEqual(hi - lo, F(1, 10**45))
            self.assertLess(verifier.poly_value(polynomial, lo)
                            * verifier.poly_value(polynomial, hi), 0)
        # Exact rational, decimal-endpoint certificates quoted in the note.
        for lo, hi in (
            ("0.637350672391423878287347192554", "0.637350672391423878287347192555"),
            ("12.442615642162591815264324420886", "12.442615642162591815264324420887"),
        ):
            self.assertLess(verifier.poly_value(polynomial, F(lo))
                            * verifier.poly_value(polynomial, F(hi)), 0)

    def test_exact_ranking_witnesses_and_sign_direction(self):
        for mu, expected_c_minus_m in (
            (F(1, 4), F(16606, 4441437)),
            (F(5, 2), -F(101, 7524)),
            (F(25), F(38275, 14281344)),
        ):
            difference = (verifier.formula_control_time(mu, "C")
                          - verifier.formula_control_time(mu, "M"))
            self.assertEqual(difference, expected_c_minus_m)

    def test_no_sensitive_resurrection_after_extinction(self):
        for state in ((0, 0), (0, 1), (0, 2)):
            c = verifier.transitions("C", state, F(7))
            self.assertEqual(c, verifier.transitions("M", state, F(7)))
            self.assertTrue(all(dest is None or dest[0] == 0 for dest in c))
        self.assertEqual(verifier.transitions("C", (0, 2), F(7)),
                         {None: F(2), (0, 1): F(6)})

    def test_progression_is_total_above_two_not_first_two_resistant(self):
        c = verifier.transitions("C", (1, 1), F(2))
        m = verifier.transitions("M", (1, 1), F(2))
        self.assertNotIn(None, c)
        self.assertEqual(c[(0, 2)], 3)  # Birth + conversion, still total 2.
        self.assertEqual(m[None], 1)  # Uncompensated R birth, total 3.
        self.assertEqual(m[(0, 2)], 2)  # Conversion conserves total burden.

    def test_zero_mutation_and_eventual_probability_discontinuity(self):
        for policy in ("C", "M"):
            self.assertEqual(verifier.formula_laplace(F(0), policy), 0)
            self.assertEqual(verifier.formula_control_time(F(0), policy), 1)
            self.assertEqual(verifier.markov_laplace(F(0), policy, eventual=True), 0)
        self.assertEqual(verifier.markov_laplace(F(1, 10**8), "C", eventual=True),
                         F(103, 112))

    def test_eventual_progression_always_favors_mtd_in_main_example(self):
        for mu in (F(1, 10**6), F(1, 4), F(5, 2), F(25), F(10**6)):
            c = verifier.markov_laplace(mu, "C", eventual=True)
            m = verifier.markov_laplace(mu, "M", eventual=True)
            self.assertLess(m, c)
            self.assertEqual(1 - c, F(9, 112))
            self.assertGreaterEqual(1 - m, (2 / (2 + mu))**2)

    def test_falsification_regimes_have_no_two_crossings(self):
        # Every coefficient has the stated sign, proving it for all mu>0.
        for rates, desired_sign in (
            (verifier.Rates(a=1), 1),
            (verifier.Rates(z=F(1, 10)), 1),
            (verifier.Rates(z=10), -1),
            (verifier.Rates(k=20), 1),
        ):
            polynomial = verifier.comparison_polynomial(rates)
            self.assertTrue(all(c * desired_sign > 0 for c in polynomial))
            for mu in (F(1, 4), F(5, 2), F(25)):
                difference = (verifier.markov_laplace(mu, "C", rates)
                              - verifier.markov_laplace(mu, "M", rates))
                self.assertGreater(difference * desired_sign, 0)

    def test_cancelled_high_mutation_cubic_is_not_a_tie(self):
        # Real leading-coefficient cancellation: both signs of the quadratic
        # are possible. This guards the boundary excluded from the cubic
        # high-rate winner rule and the direction U_M-U_C.
        witnesses = (
            (verifier.Rates(a=2, b=1, d=1, k=F(1, 2), z=1),
             (-F(9, 4), -F(15, 4), -F(3, 2)), -1),
            (verifier.Rates(a=100, b=1, d=3, k=F(13, 2), z=1),
             (F(28765, 52), F(45197, 52), F(2633, 26)), 1),
        )
        for rates, expected, winner in witnesses:
            with self.subTest(rates=rates):
                q, h, *_ = verifier.reduced_constants(rates)
                self.assertEqual(4 * rates.k * (1 - q) - rates.z, 0)
                polynomial = verifier.comparison_polynomial(rates)
                self.assertEqual(polynomial, expected)
                self.assertEqual(len(polynomial), 3)
                for mu in (F(1, 1000), F(1), F(10**8)):
                    gap = (verifier.formula_control_time(mu, "M", rates)
                           - verifier.formula_control_time(mu, "C", rates))
                    independent = (verifier.markov_laplace(mu, "C", rates)
                                   - verifier.markov_laplace(mu, "M", rates)) / rates.z
                    self.assertEqual(gap, independent)
                    self.assertGreater(winner * gap, 0)
                # With cubic term zero the decay is mu^-2, not mu^-1.
                mu = F(10**8)
                observed = mu**2 * (verifier.formula_control_time(mu, "M", rates)
                                    - verifier.formula_control_time(mu, "C", rates))
                self.assertLess(abs(observed - h * polynomial[2] / (2 * rates.z)),
                                F(1, 100))

    def test_fixed_horizon_changes_ranking_with_exact_error_bounds(self):
        c = verifier.fixed_horizon_interval(F(5, 2), "C", F(1))
        m = verifier.fixed_horizon_interval(F(5, 2), "M", F(1))
        self.assertGreater(c[0] - m[1], F(58, 1000))
        self.assertLess((c[1] - c[0]) + (m[1] - m[0]), F(1, 10**75))
        self.assertLess(verifier.formula_control_time(F(5, 2), "C")
                        - verifier.formula_control_time(F(5, 2), "M"), 0)
        self.assertEqual(verifier.fixed_horizon_interval(F(5, 2), "C", F(0)),
                         (F(0), F(0)))
        with self.assertRaises(ValueError):
            verifier.fixed_horizon_interval(F(5, 2), "C", F(1), terms=1)

    def test_invalid_inputs_fail_and_conversion_rates_are_not_probabilities(self):
        for kwargs in ({"d": 0}, {"k": -1}, {"z": 0}, {"a": 0}, {"a": 1, "b": 2}):
            with self.assertRaises(ValueError):
                verifier.Rates(**kwargs)
        with self.assertRaises(TypeError):
            verifier.Rates(a=1.0)
        with self.assertRaises(ValueError):
            verifier.formula_laplace(F(-1), "C")
        with self.assertRaises(TypeError):
            verifier.formula_laplace(0.1, "C")
        self.assertGreater(verifier.formula_laplace(F(1000), "C"), 0)

    def test_independent_fixed_seed_event_simulations(self):
        # Shared clocks couple the two arms only to reduce the variance of
        # their difference. Neither simulation arm imports deterministic code.
        for offset, mu in enumerate((F(1, 4), F(5, 2), F(25))):
            with self.subTest(mu=mu):
                result = simulation.simulate_pair(
                    float(mu), samples=300000, seed=20261117 + offset)
                exact_c = float(verifier.formula_control_time(mu, "C"))
                exact_m = float(verifier.formula_control_time(mu, "M"))
                exact_gap = exact_c - exact_m
                self.assertLess(abs(result["gap"] - exact_gap),
                                6 * result["gap_standard_error"])
                self.assertGreater(result["gap"] * exact_gap, 0)
                self.assertGreater(abs(result["gap"]),
                                   6 * result["gap_standard_error"])
                for policy, exact in (("c", exact_c), ("m", exact_m)):
                    self.assertLess(abs(result["mean_" + policy] - exact),
                                    6 * result["standard_error_" + policy])


if __name__ == "__main__":
    unittest.main()
