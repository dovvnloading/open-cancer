"""Independent finite-vertex checks for fixed-predictor target deletion."""
from fractions import Fraction
import importlib.util
import itertools
import math
from pathlib import Path
import random
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "barcode_target_robustness", ROOT / "research/barcode-lineage-memory/robustness.py")
BOUND = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BOUND)


def vertices(probabilities, contrasts, *, budget=None, weighted_removal=None):
    """Enumerate all box/intersection vertices, without sorting the contrasts."""
    n = len(probabilities)
    for free in range(n):
        rest = [i for i in range(n) if i != free]
        for bits in itertools.product((0, 1), repeat=n-1):
            removed = [Fraction(0)] * n
            for index, bit in zip(rest, bits):
                removed[index] = probabilities[index] * bit
            if budget is not None:
                removed[free] = budget - sum(removed)
            elif contrasts[free] != 0:
                removed[free] = (weighted_removal-sum(e*d for e, d in zip(removed, contrasts))) / contrasts[free]
            else:
                continue
            if 0 <= removed[free] <= probabilities[free]:
                yield removed


def exact_bound(p, d, epsilon):
    return min(sum((pi-ei)*di for pi, ei, di in zip(p, e, d)) / (1-epsilon)
               for e in vertices(p, d, budget=epsilon))


def exact_tip(p, d):
    gain = sum(pi*di for pi, di in zip(p, d))
    if gain <= 0:
        return Fraction(0)
    candidates = [sum(e) for e in vertices(p, d, weighted_removal=gain) if sum(e) < 1]
    return min(candidates) if candidates else None


class TargetDeletionTests(unittest.TestCase):
    def test_exact_vertex_optimization_on_seeded_rational_tables(self):
        rng = random.Random(11811)
        for n in range(2, 6):
            for _ in range(12):
                weights = [rng.randrange(1, 7) for _ in range(n)]
                p = [Fraction(value, sum(weights)) for value in weights]
                d = [Fraction(rng.randrange(-5, 6)) for _ in range(n)]
                for k in range(8):
                    epsilon = Fraction(k, 8)
                    actual = BOUND.adverse_deletion_bound(p, d, epsilon)["lower_gain"]
                    self.assertAlmostEqual(actual, float(exact_bound(p, d, epsilon)), places=11)

    def test_returned_witness_attains_bound_and_respects_each_capacity(self):
        p, d, epsilon = [0.1, 0.2, 0.3, 0.4], [-2, 4, 4, -1], 0.37
        result = BOUND.adverse_deletion_bound(p, d, epsilon, return_witness=True)
        removed, retained = result["deleted_probabilities"], result["retained_probabilities"]
        self.assertAlmostEqual(math.fsum(removed), epsilon)
        self.assertAlmostEqual(math.fsum(retained), 1-epsilon)
        for before, erase, keep in zip(p, removed, retained):
            self.assertGreaterEqual(erase, 0)
            self.assertLessEqual(erase, before)
            self.assertAlmostEqual(erase+keep, before)
        attained = math.fsum(mass*value for mass, value in zip(retained, d)) / math.fsum(retained)
        self.assertAlmostEqual(result["lower_gain"], attained)

    def test_fractional_boundary_is_required_and_attained(self):
        result = BOUND.adverse_deletion_bound([0.2, 0.8], [3, -1], 0.1, return_witness=True)
        self.assertEqual(result["deleted_probabilities"], [0.1, 0.0])
        self.assertAlmostEqual(result["lower_gain"], -0.5/0.9)

    def test_tipping_matches_independent_rational_vertex_optimization(self):
        rng = random.Random(619)
        for n in range(2, 6):
            for _ in range(40):
                weights = [rng.randrange(1, 8) for _ in range(n)]
                p = [Fraction(value, sum(weights)) for value in weights]
                d = [Fraction(rng.randrange(-4, 6)) for _ in range(n)]
                expected, actual = exact_tip(p, d), BOUND.tipping_mass(p, d)
                if expected is None:
                    self.assertEqual(actual["status"], "unreachable_below_one")
                    self.assertIsNone(actual["tipping_mass"])
                else:
                    self.assertAlmostEqual(actual["tipping_mass"], float(expected), places=12)

    def test_tipping_cases_keep_zero_budget_and_zero_tail_distinct(self):
        cases = [([0.2, 0.3, 0.5], [4, 1, -1], 0.15, "attained"),
                 ([0.2, 0.8], [1, 0], 0.2, "attained"),
                 ([0.1, 0.9], [1, -1], 0.0, "already_nonpositive"),
                 ([1], [3], None, "unreachable_below_one")]
        for p, d, expected, status in cases:
            found = BOUND.tipping_mass(p, d)
            self.assertEqual(found["status"], status)
            if expected is None:
                self.assertIsNone(found["tipping_mass"])
            else:
                self.assertAlmostEqual(found["tipping_mass"], expected)

    def test_unoccupied_categories_and_permutation_do_not_change_answer(self):
        p, d = [0.3, 0.7, 0], [3, -1, -1e9]
        baseline = BOUND.adverse_deletion_bound(p, d, 0.1)["lower_gain"]
        for order in itertools.permutations(range(3)):
            got = BOUND.adverse_deletion_bound([p[i] for i in order], [d[i] for i in order], .1)
            self.assertAlmostEqual(got["lower_gain"], baseline)
        self.assertEqual(BOUND.tipping_mass([1, 0], [2, -1e9])["status"], "unreachable_below_one")

    def test_count_scaling_and_contrast_translation(self):
        original = BOUND.adverse_deletion_bound([0.2, 0.8], [2, -1], 0.17)["lower_gain"]
        changed = BOUND.adverse_deletion_bound([20, 80], [9, 6], 0.17)["lower_gain"]
        self.assertAlmostEqual(changed, original+7)

    def test_larger_adverse_budget_never_improves_the_retained_mean(self):
        values = [BOUND.adverse_deletion_bound([.1, .2, .3, .4], [4, -2, 3, 1], e)["lower_gain"]
                  for e in (0, .01, .1, .4, .6, .8, .99)]
        self.assertTrue(all(a >= b-1e-12 for a, b in zip(values, values[1:])))

    def test_constant_contrasts_survive_every_valid_budget(self):
        for value in (-3, 0, 2):
            for epsilon in (0, .001, .5, .999999):
                self.assertEqual(BOUND.adverse_deletion_bound([2, 3, 5], [value]*3, epsilon)["lower_gain"], value)

    def test_invalid_probability_inputs_are_rejected(self):
        for p in ([], [0, 0], [-1, 2], [math.inf, 1], [math.nan, 1]):
            with self.assertRaises(ValueError):
                BOUND.adverse_deletion_bound(p, [1]*len(p), .1)
        with self.assertRaises(ValueError):
            BOUND.tipping_mass([1], [1, 2])

    def test_invalid_budgets_and_nonfinite_contrasts_are_rejected(self):
        for epsilon in (-.1, 1, 2, math.nan, math.inf):
            with self.assertRaises(ValueError):
                BOUND.adverse_deletion_bound([1], [1], epsilon)
        for contrast in (math.nan, math.inf, -math.inf):
            with self.assertRaises(ValueError):
                BOUND.tipping_mass([1], [contrast])

    def test_near_one_budget_preserves_a_real_positive_tail(self):
        epsilon = 1-1e-8
        lower = BOUND.adverse_deletion_bound([.4, .6], [8, 2], epsilon)["lower_gain"]
        self.assertAlmostEqual(lower, 2, places=7)
        self.assertEqual(BOUND.tipping_mass([.4, .6], [8, 2])["status"], "unreachable_below_one")


if __name__ == "__main__":
    unittest.main()
