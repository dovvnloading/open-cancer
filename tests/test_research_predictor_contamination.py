"""Independent high-precision checks of explicitly synthetic probability claims."""
from decimal import Decimal, localcontext
import importlib.util
import json
import math
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "research/barcode-lineage-memory" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SYNTHETIC = load("barcode_synthetic_counterexample", "predictor_counterexample.py")
BOUND = load("barcode_synthetic_target_check", "robustness.py")


class PredictorContaminationTests(unittest.TestCase):
    def test_decimal_fourth_roots_reproduce_numeric_examples(self):
        with localcontext() as ctx:
            ctx.prec = 80
            w, epsilon = Decimal("0.2"), Decimal("0.0001")
            for text in ("0.001", "0.00001", "0.00000001", "1e-12", "1e-20"):
                m = Decimal(text)
                x = (w, m, 1-w-m)
                z = (m, 1-2*m, m)
                dirty = ((1-epsilon)*m+epsilon, (1-epsilon)*(1-2*m), (1-epsilon)*m)
                actual = SYNTHETIC.counterexample(m=float(m))
                for donor, key in ((z, "clean"), (dirty, "contaminated")):
                    weights = [a*b.sqrt().sqrt() for a, b in zip(x, donor)]
                    probability = weights[0]/sum(weights)
                    gain = (probability/w).ln()
                    self.assertAlmostEqual(actual[key+"_target_prediction"], float(probability), places=13)
                    self.assertAlmostEqual(actual[key+"_gain"], float(gain), places=12)

    def test_fixed_budget_limit_has_clean_zero_and_contaminated_finite_gain(self):
        rows = SYNTHETIC.default_examples()
        self.assertLess(abs(rows[-1]["clean_gain"]), 1e-12)
        self.assertAlmostEqual(rows[-1]["clean_target_prediction"], .2, places=12)
        self.assertLess(1-rows[-1]["contaminated_target_prediction"], .0005)
        self.assertLess(abs(rows[-1]["contaminated_gain"]+math.log(.2)), .0005)
        self.assertGreater(rows[-1]["contaminated_gain"], rows[0]["contaminated_gain"])

    def test_arbitrarily_smaller_donor_budgets_admit_positive_synthetic_witnesses(self):
        for epsilon in (1e-3, 1e-6, 1e-9, 1e-12):
            row = SYNTHETIC.counterexample(epsilon=epsilon, m=epsilon**8)
            self.assertLess(abs(row["clean_gain"]), 1e-8)
            self.assertGreater(row["contaminated_gain"], 1.6)
            self.assertLessEqual(row["donor_total_variation"], epsilon*(1+1e-12))

    def test_smaller_baseline_probability_allows_larger_ordered_limit(self):
        gains = []
        for w in (.1, .01, .001, .000001):
            row = SYNTHETIC.counterexample(w=w, m=w**16*1e-32)
            gains.append(row["contaminated_gain"])
            self.assertAlmostEqual(row["contaminated_gain"], -math.log(w), places=8)
        self.assertTrue(all(a < b for a, b in zip(gains, gains[1:])))

    def test_every_target_only_budget_below_one_keeps_the_same_pointmass_gain(self):
        row = SYNTHETIC.counterexample()
        base = [math.log(value) for value in row["own"]]
        dirty = SYNTHETIC.product_log_prediction(row["own"], row["contaminated_donor"])
        contrast = [x-y for x, y in zip(dirty, base)]
        for epsilon in (0, .001, .01, .5, .99, 1-1e-12):
            got = BOUND.adverse_deletion_bound([1, 0, 0], contrast, epsilon)
            self.assertAlmostEqual(got["lower_gain"], row["contaminated_gain"], places=13)
        self.assertEqual(BOUND.tipping_mass([1, 0, 0], contrast)["status"], "unreachable_below_one")

    def test_fixed_floor_bound_covers_each_category_and_tends_to_zero(self):
        own, donor = [.05, .25, .7], [.2, .3, .5]
        mixtures = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [.1, .8, .1]]
        for b in (0, .25, 1, 2):
            clean = SYNTHETIC.product_log_prediction(own, donor, 1.5, b)
            for epsilon in (0, 1e-8, .001, .1, .9):
                bound = SYNTHETIC.fixed_donor_mixture_score_bound(epsilon, min(donor), b)
                for mixture in mixtures:
                    contaminated = [(1-epsilon)*z+epsilon*v for z, v in zip(donor, mixture)]
                    got = SYNTHETIC.product_log_prediction(own, contaminated, 1.5, b)
                    self.assertLessEqual(max(abs(a-c) for a, c in zip(got, clean)), bound+5e-15)
        self.assertLess(SYNTHETIC.fixed_donor_mixture_score_bound(1e-9, .2), 2e-9)
        self.assertTrue(math.isfinite(SYNTHETIC.fixed_donor_mixture_score_bound(.9, 1e-320)))

    def test_raw_contamination_before_smoothing_has_the_same_limits(self):
        row = SYNTHETIC.before_smoothing_example(delta=1e-24)
        self.assertAlmostEqual(row["clean_target_prediction"], .2, places=12)
        self.assertLess(1-row["contaminated_target_prediction"], 1e-4)
        self.assertGreater(row["contaminated_gain"], 1.6)
        self.assertAlmostEqual(row["smoothed_donor_total_variation"], row["epsilon"]*(1-row["delta"]))

    def test_invalid_domains_and_nonfinite_output_are_prevented(self):
        for kwargs in ({"w": 0}, {"w": 1}, {"m": 0}, {"m": .3}, {"b": 0}, {"b": 1},
                       {"epsilon": 0}, {"epsilon": 1}, {"m": math.nan}):
            with self.assertRaises(ValueError):
                SYNTHETIC.counterexample(**kwargs)
        for args in ((-1, .1), (1, .1), (.1, 0), (.1, 2), (.1, .1, -1)):
            with self.assertRaises(ValueError):
                SYNTHETIC.fixed_donor_mixture_score_bound(*args)
        with self.assertRaises(ValueError):
            SYNTHETIC.product_log_prediction([1, 0], [.5, .5])
        json.dumps(SYNTHETIC.default_examples(), allow_nan=False)


if __name__ == "__main__":
    unittest.main()
