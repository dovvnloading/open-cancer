"""Offline, stdlib-only checks of the NONCLINICAL Merkel aggregate audit.

The reference calculations below use separate least-squares and probability
implementations. Reproducibility checks do not establish source truth, assay
calibration, biological causation, or clinical usefulness.
"""

import copy
import csv
from decimal import Decimal
import hashlib
import importlib.util
import io
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "merkel-viral-antigen-test"
ORIGINAL = ROOT / "research" / "immune-escape-pancancer-test"
PROTOCOL_SHA256 = "22923506eff606d81284a476a56813ec0773d73d238a5dcf633067344e711ea1"
CALIBRATION_SHA256 = "ef597c8dd3ff770848fb06092853e620136a014c1fdd76c9e3873618cbf3af6e"


def csv_reference(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256_reference(path):
    data = Path(path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


def clipped(value):
    return min(1.0, max(0.0, value))


def reciprocal_reference(rows):
    """Independent dense search with a projected analytic ceiling.

    Unlike production's logarithmic grid plus golden-section search, refine
    every locally minimal grid bracket by ternary minimization. This checks
    the original deterministic fit without fitting to any JAVELIN outcome.
    """
    pairs = [(float(row["tmb"]), float(row["orr"])) for row in rows]

    def objective(b):
        shape = [math.exp(-b / tmb) for tmb, _ in pairs]
        squares = sum(value * value for value in shape)
        a = clipped(sum(value * pair[1] for value, pair in zip(shape, pairs))
                    / squares) if squares else 0.0
        loss = sum((pair[1] - a * value) ** 2
                   for value, pair in zip(shape, pairs))
        return loss, a, b

    minimum = min(tmb for tmb, _ in pairs) * 1e-9
    maximum = max(tmb for tmb, _ in pairs) * 100
    grid = [0.0] + [minimum * (maximum / minimum) ** (index / 240)
                    for index in range(241)]
    evaluations = [objective(b) for b in grid]
    candidates = list(evaluations)
    for index in range(1, len(grid) - 1):
        if (evaluations[index][0] <= evaluations[index - 1][0]
                and evaluations[index][0] <= evaluations[index + 1][0]):
            left, right = grid[index - 1], grid[index + 1]
            for _ in range(110):
                first = left + (right - left) / 3
                second = right - (right - left) / 3
                if objective(first)[0] < objective(second)[0]:
                    right = second
                else:
                    left = first
            candidates.extend((objective(left), objective(right),
                               objective((left + right) / 2)))
    loss, a, b = min(candidates)
    return {"A": a, "B": b, "sse": loss}


def log_reference(rows):
    """OLS from pairwise differences, independent of production centering."""
    points = [(math.log(float(row["tmb"])), float(row["orr"])) for row in rows]
    numerator = sum((x1 - x2) * (y1 - y2)
                    for x1, y1 in points for x2, y2 in points)
    denominator = sum((x1 - x2) ** 2
                      for x1, _ in points for x2, _ in points)
    slope = numerator / denominator if denominator else 0.0
    intercept = sum(y - slope * x for x, y in points) / len(points)
    return {"a": intercept, "b": slope}


def probability_reference(curve, value):
    if "A" in curve:
        if value <= 0:
            return 0.0 if curve["B"] > 0 else curve["A"]
        return curve["A"] * math.exp(-curve["B"] / value)
    if value <= 0:
        return 0.0
    return clipped(curve["a"] + curve["b"] * math.log(value))


def inverse_reference(curve, median, lower, upper):
    """Preimage of a closed probability interval, in finite positive scales."""
    if lower > upper:
        return None
    if "A" in curve:
        a, b = curve["A"], curve["B"]
        if b == 0 or a == 0:
            return (0.0, math.inf) if lower <= a <= upper else None
        if upper <= 0 or lower >= a:
            return None
        low = 0.0 if lower <= 0 else b / (median * math.log(a / lower))
        high = math.inf if upper >= a else b / (median * math.log(a / upper))
    else:
        a, b = curve["a"], curve["b"]
        if b == 0:
            return (0.0, math.inf) if lower <= clipped(a) <= upper else None
        if b < 0:
            raise ValueError("Reference restricts the frozen log fit to positive slope")
        low = 0.0 if lower <= 0 else math.exp((lower - a) / b) / median
        high = math.inf if upper >= 1 else math.exp((upper - a) / b) / median
    return None if low > high or math.isinf(low) or high <= 0 else (low, high)


def intersection_reference(*intervals):
    if not intervals or any(interval is None for interval in intervals):
        return None
    low = max(interval[0] for interval in intervals)
    high = min(interval[1] for interval in intervals)
    return None if low > high or math.isinf(low) or high <= 0 else (low, high)


def rounding_reference(value, percentage=False):
    value = Decimal(str(value))
    half_unit = Decimal(5).scaleb(value.as_tuple().exponent - 1)
    return float(half_unit / 100 if percentage else half_unit)


def interval_distance_reference(prediction, lower, upper):
    return max(lower - prediction, prediction - upper, 0.0)


def clopper_pearson_reference(successes, trials):
    """Invert finite binomial tails only where source publishes both counts."""
    def tail(p, start, stop):
        return sum(math.comb(trials, k) * p ** k * (1 - p) ** (trials - k)
                   for k in range(start, stop))

    def bisect(function, increasing):
        left, right = 0.0, 1.0
        for _ in range(70):
            midpoint = (left + right) / 2
            if (function(midpoint) < 0.025) == increasing:
                left = midpoint
            else:
                right = midpoint
        return (left + right) / 2

    lower = bisect(lambda p: tail(p, successes, trials + 1), True) if successes else 0.0
    upper = bisect(lambda p: tail(p, 0, successes + 1), False) if successes < trials else 1.0
    return lower, upper


def panel_reference(curve, inputs, mode):
    result = {}
    for cohort in ("A", "B"):
        groups = {}
        for status in ("positive", "negative"):
            prefix = cohort.lower()
            marker = "_subset_orr_" if mode == "subset_a" and cohort == "A" else "_orr_"
            outcome, tmb = inputs[prefix + marker + status], inputs[prefix + "_tmb_" + status]
            median = float(tmb["value"])
            lower, upper = float(outcome["ci_low"]) / 100, float(outcome["ci_high"]) / 100
            if mode == "rounded":
                error = rounding_reference(tmb["value"])
                lower = max(0, lower - rounding_reference(outcome["ci_low"], True))
                upper = min(1, upper + rounding_reference(outcome["ci_high"], True))
                at_low = inverse_reference(curve, median - error, lower, upper)
                at_high = inverse_reference(curve, median + error, lower, upper)
                region = None if at_high is None else (at_high[0], at_low[1])
            else:
                if mode == "median_envelope":
                    lower, upper = max(0, 2 * lower - curve["A"]), min(curve["A"], 2 * upper)
                region = inverse_reference(curve, median, lower, upper)
            groups[status] = region
        groups["joint"] = intersection_reference(groups["positive"], groups["negative"])
        result[cohort] = groups
    result["joint"] = intersection_reference(result["A"]["joint"], result["B"]["joint"])
    result["bands"] = {
        "counting_assumption": intersection_reference(result["joint"], (1.0, 2.0)),
        "wide_assumption": intersection_reference(result["joint"], (0.5, 5.0)),
    }
    return result


def score_reference(predictions, inputs):
    outcomes = [inputs["b_orr_" + status] for status in ("positive", "negative")]
    residuals = [float(row["value"]) / 100 - prediction for row, prediction in zip(outcomes, predictions)]
    distances = [interval_distance_reference(prediction, float(row["ci_low"]) / 100,
                                             float(row["ci_high"]) / 100)
                 for row, prediction in zip(outcomes, predictions)]
    mse = sum(value * value for value in residuals) / 2
    return {"mse": mse, "rmse": math.sqrt(mse),
            "mae": sum(abs(value) for value in residuals) / 2,
            "misses": sum(value > 0 for value in distances),
            "interval_mse": sum(value * value for value in distances) / 2,
            "residuals": residuals, "distances": distances}


def viral_reference(curve, median, target, scale):
    baseline = probability_reference(curve, median * scale)
    if target <= baseline:
        return 0.0
    if target >= curve["A"]:
        return None
    return curve["B"] / math.log(curve["A"] / target) - median * scale


def load_analysis():
    spec = importlib.util.spec_from_file_location(
        "merkel_viral_antigen_analysis", RESEARCH / "analysis.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MerkelSyntheticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.analysis = load_analysis()

    def assertRegion(self, actual, expected, places=10):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        for observed, reference in zip(actual, expected):
            if math.isinf(reference):
                self.assertEqual(observed, reference)
            else:
                self.assertAlmostEqual(observed, reference, places=places)

    def test_frozen_protocol_fingerprint(self):
        self.assertEqual(sha256_reference(RESEARCH / "PROTOCOL.md"), PROTOCOL_SHA256)

    def test_line_endings_have_identical_fingerprints(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory) / name for name in ("lf", "crlf", "cr")]
            for path, separator in zip(paths, (b"\n", b"\r\n", b"\r")):
                path.write_bytes(separator.join((b"alpha", b"beta", b"")))
            expected = hashlib.sha256(b"alpha\nbeta\n").hexdigest()
            for path in paths:
                self.assertEqual(self.analysis.normalized_bytes(path), b"alpha\nbeta\n")
                self.assertEqual(self.analysis.digest(path), expected)

    def test_reciprocal_inverse_and_unattained_ceilings(self):
        curve = {"kind": "reciprocal", "A": 0.6, "B": 2.7}
        median = 0.8
        cases = ((0.1, 0.4), (0.0, 0.6), (0.6, 0.6),
                 (0.7, 1.0), (0.0, 0.0), (0.0, 1.0), (0.3, 1.0))
        for lower, upper in cases:
            with self.subTest(lower=lower, upper=upper):
                actual = self.analysis.compatibility(curve, median, lower, upper)
                expected = inverse_reference(curve, median, lower, upper)
                self.assertRegion(actual, expected)
                if actual is None:
                    continue
                if actual[0] > 0:
                    self.assertAlmostEqual(self.analysis.predict(curve, median * actual[0]),
                                           lower, places=13)
                if math.isfinite(actual[1]):
                    self.assertAlmostEqual(self.analysis.predict(curve, median * actual[1]),
                                           upper, places=13)
        for scale in (1e-4, 0.5, 1.0, 2.0, 1000.0):
            self.assertAlmostEqual(self.analysis.predict(curve, scale * median),
                                   probability_reference(curve, scale * median), places=14)

    def test_log_inverse_includes_attained_clipping_boundaries(self):
        curve = {"kind": "log", "a": 0.15, "b": 0.12}
        for median in (0.2, 1.0, 4.0):
            for lower, upper in ((0, 0), (1, 1), (0, 1), (0.2, 0.8), (0.4, 1)):
                with self.subTest(median=median, lower=lower, upper=upper):
                    self.assertRegion(
                        self.analysis.compatibility(curve, median, lower, upper),
                        inverse_reference(curve, median, lower, upper))
        self.assertEqual(self.analysis.predict(curve, math.exp(-20)), 0)
        self.assertEqual(self.analysis.predict(curve, math.exp(20)), 1)

    def test_constant_curve_boundaries(self):
        for curve in ({"kind": "reciprocal", "A": 0.6, "B": 0.0},
                      {"kind": "reciprocal", "A": 0.0, "B": 2.0},
                      {"kind": "log", "a": 0.6, "b": 0.0}):
            value = probability_reference(curve, 1.0)
            self.assertRegion(self.analysis.compatibility(curve, 1.0, value, value),
                              (0.0, math.inf))
            self.assertIsNone(self.analysis.compatibility(curve, 1.0, value + 0.1, 1.0))

    def test_intersection_requires_one_shared_finite_positive_scale(self):
        cases = (((1, 2), (2, 5)), ((0, math.inf), (4, 7)),
                 ((1, 2), (3, math.inf)), ((1, 2), None),
                 ((math.inf, math.inf),), ((0, 0),))
        for regions in cases:
            with self.subTest(regions=regions):
                self.assertRegion(self.analysis.intersect(*regions),
                                  intersection_reference(*regions))

    def test_median_envelope_has_correct_half_mass_limits(self):
        curve = {"kind": "reciprocal", "A": 0.6, "B": 2.7}
        for median in (0.2, 1.0, 4.0):
            for lower, upper in ((0.1, 0.4), (0.4, 0.5), (0.6, 0.8), (0.7, 0.8), (0, 0)):
                expected = inverse_reference(curve, median,
                                             max(0, 2 * lower - curve["A"]),
                                             min(curve["A"], 2 * upper))
                self.assertRegion(self.analysis.median_compatibility(
                    curve, median, lower, upper), expected)
            # Every finite odd sample below has the specified sample median.
            for values in ((0, 0, median, median, median),
                           (median, median, median, 1000, 1000),
                           (0, median / 2, median, 2 * median, 1000)):
                midpoint = probability_reference(curve, median)
                average = sum(probability_reference(curve, value) for value in values) / len(values)
                self.assertGreaterEqual(average, midpoint / 2)
                self.assertLessEqual(average, (curve["A"] + midpoint) / 2)

    def test_viral_parameter_boundaries_do_not_invent_finite_ceiling_fit(self):
        curve = {"kind": "reciprocal", "A": 0.6, "B": 2.7}
        median, scale = 0.8, 2.0
        baseline = probability_reference(curve, median * scale)
        for target in (0.0, baseline / 2, baseline):
            actual = self.analysis.viral_calibration(curve, median, target, scale)
            self.assertEqual(actual["v"], 0.0)
            self.assertEqual(actual["status"], "zero_boundary")
            self.assertAlmostEqual(actual["fitted_orr"], baseline, places=14)
        for target in (curve["A"], 0.8, 1.0):
            actual = self.analysis.viral_calibration(curve, median, target, scale)
            self.assertIsNone(actual["v"])
            self.assertIsNone(actual["fitted_orr"])
            self.assertEqual(actual["status"], "unattained_ceiling")
        target = (baseline + curve["A"]) / 2
        fit = self.analysis.viral_calibration(curve, median, target, scale)
        expected = curve["B"] / math.log(curve["A"] / target) - median * scale
        self.assertAlmostEqual(fit["v"], expected, places=13)
        self.assertAlmostEqual(fit["fitted_orr"], target, places=13)


class MerkelAggregateTests(unittest.TestCase):
    assertRegion = MerkelSyntheticTests.assertRegion

    @classmethod
    def setUpClass(cls):
        cls.analysis = load_analysis()
        cls.inputs = cls.analysis.load_inputs()
        cls.curves = cls.analysis.frozen_curves()
        cls.report = cls.analysis.build_report()
        cls.original_rows = [row for row in csv_reference(ORIGINAL / "data" / "tmb_orr.csv")
                             if row["tmb"]]

    def test_original_fits_reconstructed_by_independent_optimizers(self):
        self.assertEqual(set(self.curves), {"reciprocal_holdout", "reciprocal_all",
                                          "log_holdout", "log_all"})
        for suffix in ("holdout", "all"):
            rows = [row for row in self.original_rows
                    if suffix == "all" or row["parent_type"] != "merkel"]
            references = {"reciprocal": reciprocal_reference(rows), "log": log_reference(rows)}
            for kind, reference in references.items():
                actual = self.curves[kind + "_" + suffix]
                self.assertEqual(actual["n"], len(rows))
                self.assertEqual(actual["parents"], len({r["parent_type"] for r in rows}))
                for key in (("A", "B") if kind == "reciprocal" else ("a", "b")):
                    self.assertAlmostEqual(actual[key], reference[key], delta=2e-6)
                for row in rows:
                    self.assertAlmostEqual(self.analysis.predict(actual, float(row["tmb"])),
                                           probability_reference(actual, float(row["tmb"])), places=13)

    def test_frozen_curves_never_read_javelin_inputs(self):
        original_reader = self.analysis.read_csv

        def guarded_reader(path):
            self.assertIn(Path(path), (ORIGINAL / "data" / "tmb_orr.csv",
                                      RESEARCH / "data" / "frozen_curves.csv"))
            return original_reader(path)

        self.analysis.frozen_curves.cache_clear()
        with patch.object(self.analysis, "read_csv", side_effect=guarded_reader):
            self.assertEqual(self.analysis.frozen_curves(), self.curves)

    def test_source_snapshots_are_stable_under_small_replay_drift_and_reject_mismatch(self):
        snapshot_rows = csv_reference(RESEARCH / "data" / "frozen_curves.csv")
        self.assertEqual(len(snapshot_rows), len(self.curves))
        for row in snapshot_rows:
            for key in ("source_title", "source_url", "locator", "access_date", "reuse_basis", "license_url"):
                self.assertTrue(row[key].strip())
            self.assertEqual(row["pr67_code_sha256"], sha256_reference(ORIGINAL / "analysis.py"))
            self.assertEqual(row["pr67_csv_sha256"], sha256_reference(ORIGINAL / "data" / "tmb_orr.csv"))
        small_drift = copy.deepcopy(self.curves)
        for curve in small_drift.values():
            for key in (("A", "B") if curve["kind"] == "reciprocal" else ("a", "b")):
                curve[key] += 1e-7
        try:
            self.analysis.frozen_curves.cache_clear()
            with patch.object(self.analysis, "replay_curves", return_value=small_drift):
                self.assertEqual(self.analysis.frozen_curves(), self.curves)
                self.assertEqual(self.analysis.build_report(), self.report)
            large_drift = copy.deepcopy(self.curves)
            large_drift["reciprocal_holdout"]["B"] += 1e-3
            self.analysis.frozen_curves.cache_clear()
            with patch.object(self.analysis, "replay_curves", return_value=large_drift):
                with self.assertRaisesRegex(ValueError, "replay exceeds numerical tolerance"):
                    self.analysis.frozen_curves()
            original_reader = self.analysis.read_csv
            for field in ("pr67_code_sha256", "pr67_csv_sha256"):
                def bad_snapshot_reader(path):
                    rows = original_reader(path)
                    if Path(path) == RESEARCH / "data" / "frozen_curves.csv":
                        rows[0][field] = "0" * 64
                    return rows
                self.analysis.frozen_curves.cache_clear()
                with patch.object(self.analysis, "read_csv", side_effect=bad_snapshot_reader):
                    with self.assertRaisesRegex(ValueError, "source fingerprint mismatch"):
                        self.analysis.frozen_curves()
        finally:
            self.analysis.frozen_curves.cache_clear()
            self.analysis.frozen_curves()

    def test_published_rows_have_provenance_and_separate_marginal_denominators(self):
        required = ("source_title", "doi", "source_url", "data_url", "locator",
                    "access_date", "reuse_basis", "license_url", "inspection", "notes")
        for row in self.inputs.values():
            with self.subTest(row=row["id"]):
                for key in required:
                    self.assertTrue(row[key].strip(), key)
                for key in ("source_url", "data_url", "license_url"):
                    self.assertTrue(row[key].startswith("https://"))
                self.assertEqual(row["access_date"], "2026-11-18")
                self.assertIn(row["cohort"], ("A", "B"))
                self.assertTrue(row["id"].startswith(row["cohort"].lower() + "_"))
                if row["metric"] == "orr_percent" and row["value"]:
                    self.assertLessEqual(float(row["ci_low"]), float(row["value"]))
                    self.assertLessEqual(float(row["value"]), float(row["ci_high"]))
                    self.assertTrue(row["n"])
                if row["metric"] == "tmb_median":
                    self.assertEqual(row["unit"], "NSSV/Mb")
                    self.assertEqual(row["analysis_set"], "tmb_evaluable")
        for cohort in ("a", "b"):
            for status in ("positive", "negative"):
                orr, tmb = (self.inputs[cohort + marker + status]
                            for marker in ("_orr_", "_tmb_"))
                self.assertGreater(int(orr["n"]), int(tmb["n"]))
                self.assertEqual(orr["analysis_set"], "all_treated")
        for status in ("positive", "negative"):
            absent = self.inputs["b_subset_orr_" + status]
            self.assertEqual([absent[k] for k in ("value", "ci_low", "ci_high", "n", "responders")],
                             [""] * 5)
            self.assertEqual(self.inputs["a_subset_orr_" + status]["n"],
                             self.inputs["a_tmb_" + status]["n"])
        self.assertNotEqual(self.inputs["a_orr_all"]["line"], self.inputs["b_orr_all"]["line"])
        self.assertNotEqual(self.inputs["a_orr_positive"]["data_cutoff"],
                            self.inputs["a_tmb_positive"]["data_cutoff"])
        self.assertEqual(self.inputs["a_orr_positive"]["responders"], "")
        self.assertEqual(self.inputs["a_orr_negative"]["responders"], "")

    def test_count_identities_and_published_count_supported_intervals(self):
        data = self.inputs
        self.assertEqual(int(data["a_orr_positive"]["n"]) + int(data["a_orr_negative"]["n"])
                         + int(data["a_count_viral_unevaluable"]["value"]),
                         int(data["a_orr_all"]["n"]))
        self.assertEqual(sum(int(data["b_orr_" + status]["n"])
                             for status in ("positive", "negative", "unknown")),
                         int(data["b_orr_all"]["n"]))
        self.assertEqual(int(data["a_tmb_positive"]["n"]) + int(data["a_tmb_negative"]["n"]),
                         int(data["a_tmb_all"]["n"]))
        self.assertLess(int(data["b_tmb_positive"]["n"]) + int(data["b_tmb_negative"]["n"]),
                        int(data["b_tmb_all"]["n"]))
        for row in data.values():
            if row["metric"] != "orr_percent" or not row["responders"]:
                continue
            with self.subTest(row=row["id"]):
                successes, trials = int(row["responders"]), int(row["n"])
                self.assertLessEqual(abs(100 * successes / trials - float(row["value"])),
                                     rounding_reference(row["value"]) + 1e-9)
                lower, upper = clopper_pearson_reference(successes, trials)
                for expected, field in ((lower, "ci_low"), (upper, "ci_high")):
                    self.assertLessEqual(abs(100 * expected - float(row[field])),
                                         rounding_reference(row[field]) + 1e-9)

    def test_primary_and_subset_compatibility_endpoints_independently(self):
        for curve_name, curve in self.curves.items():
            for row in self.inputs.values():
                if (row["metric"] != "orr_percent" or not row["value"]
                        or row["viral_status"] not in ("positive", "negative")):
                    continue
                tmb = self.inputs[row["cohort"].lower() + "_tmb_" + row["viral_status"]]
                median = float(tmb["value"])
                lower, upper = float(row["ci_low"]) / 100, float(row["ci_high"]) / 100
                with self.subTest(curve=curve_name, row=row["id"]):
                    expected = inverse_reference(curve, median, lower, upper)
                    actual = self.analysis.compatibility(curve, median, lower, upper)
                    self.assertRegion(actual, expected)
                    for scale in (0.5, 1.0, 4 / 3, 2.0, 5.0, 10.0, 100.0, 10000.0):
                        membership = actual is not None and actual[0] <= scale <= actual[1]
                        prediction = probability_reference(curve, scale * median)
                        self.assertEqual(membership, lower <= prediction <= upper)

    def test_rounding_union_and_median_relaxation_for_every_eligible_row(self):
        for value in ("0.49", "10.52", "16.0", "48.4", "1"):
            self.assertEqual(self.analysis.half_unit(value), rounding_reference(value))
        for row in self.inputs.values():
            if (row["metric"] != "orr_percent" or not row["value"]
                    or row["viral_status"] not in ("positive", "negative")):
                continue
            tmb = self.inputs[row["cohort"].lower() + "_tmb_" + row["viral_status"]]
            median, error = float(tmb["value"]), rounding_reference(tmb["value"])
            lower = max(0, float(row["ci_low"]) / 100 - rounding_reference(row["ci_low"], True))
            upper = min(1, float(row["ci_high"]) / 100 + rounding_reference(row["ci_high"], True))
            for name, curve in self.curves.items():
                low_m = inverse_reference(curve, median - error, lower, upper)
                high_m = inverse_reference(curve, median + error, lower, upper)
                expected = None if high_m is None else (high_m[0], low_m[1])
                with self.subTest(curve=name, row=row["id"]):
                    self.assertRegion(self.analysis.rounded_compatibility(curve, tmb, row), expected)
            curve = self.curves["reciprocal_holdout"]
            lower, upper = float(row["ci_low"]) / 100, float(row["ci_high"]) / 100
            expected = inverse_reference(curve, median, max(0, 2 * lower - curve["A"]),
                                         min(curve["A"], 2 * upper))
            self.assertRegion(self.analysis.median_compatibility(curve, median, lower, upper), expected)

    def test_every_committed_part_a_calibration_is_recomputed(self):
        self.assertEqual(sha256_reference(RESEARCH / "data" / "part_a_calibration.csv"),
                         CALIBRATION_SHA256)
        self.assertTrue(self.analysis.check_calibration())
        actual = self.analysis.calibrate_part_a(self.inputs, self.curves)
        self.assertEqual(len(actual), 2 * len(self.analysis.SCALES) * 3)
        curve = self.curves["reciprocal_holdout"]
        for row in actual:
            orr = self.inputs[row["source_orr_id"]]
            tmb = self.inputs[row["source_tmb_id"]]
            field = "value" if row["target"] == "point" else row["target"]
            y, median = float(orr[field]) / 100, float(tmb["value"])
            baseline = probability_reference(curve, row["scale"] * median)
            self.assertEqual(row["target_orr"], y)
            self.assertEqual(row["median"], median)
            if y <= baseline:
                self.assertEqual(row["v"], 0)
                self.assertEqual(row["status"], "zero_boundary")
                self.assertAlmostEqual(row["fitted_orr"], baseline, places=13)
            elif y >= curve["A"]:
                self.assertIsNone(row["v"])
                self.assertEqual(row["status"], "unattained_ceiling")
            else:
                expected = curve["B"] / math.log(curve["A"] / y) - row["scale"] * median
                self.assertAlmostEqual(row["v"], expected, places=12)
                self.assertAlmostEqual(row["fitted_orr"], y, places=13)
            self.assertEqual((row["A"], row["B"]), (curve["A"], curve["B"]))
        saved = csv_reference(RESEARCH / "data" / "part_a_calibration.csv")
        numerical = {"scale", "A", "B", "target_orr", "median", "v", "fitted_orr"}
        recorded = self.analysis.recorded_calibration()
        for row, recomputed, fixed in zip(saved, actual, recorded):
            self.assertEqual(row["part_a_sha256"], sha256_reference(RESEARCH / "data" / "part_a.csv"))
            self.assertEqual(row["pr67_code_sha256"], sha256_reference(ORIGINAL / "analysis.py"))
            self.assertEqual(row["pr67_csv_sha256"], sha256_reference(ORIGINAL / "data" / "tmb_orr.csv"))
            self.assertEqual(row["protocol_sha256"], PROTOCOL_SHA256)
            for key in self.analysis.CALIBRATION_FIELDS[:13]:
                if key in numerical and row[key]:
                    self.assertTrue(math.isclose(float(row[key]), recomputed[key],
                                                 rel_tol=5e-12, abs_tol=5e-12))
                    self.assertEqual(fixed[key], float(row[key]))
                elif key in numerical:
                    self.assertIsNone(recomputed[key])
                    self.assertIsNone(fixed[key])
                else:
                    self.assertEqual(row[key], recomputed[key])
                    self.assertEqual(fixed[key], row[key])
        original_text = self.analysis.calibration_text()
        drifted = list(csv.DictReader(io.StringIO(original_text)))
        drifted[0]["v"] = str(float(drifted[0]["v"]) + 1e-13)
        small_text = self.analysis.csv_text(drifted, self.analysis.CALIBRATION_FIELDS)
        with patch.object(self.analysis, "calibration_text", return_value=small_text):
            self.assertTrue(self.analysis.check_calibration())
        drifted[0]["v"] = str(float(drifted[0]["v"]) + 1e-6)
        large_text = self.analysis.csv_text(drifted, self.analysis.CALIBRATION_FIELDS)
        with patch.object(self.analysis, "calibration_text", return_value=large_text):
            with self.assertRaisesRegex(ValueError, "calibration replay mismatch"):
                self.analysis.check_calibration()

    def test_part_b_cannot_change_any_part_a_parameter(self):
        expected = self.analysis.calibrate_part_a(self.inputs, self.curves)
        perturbed = copy.deepcopy(self.inputs)
        for key, row in perturbed.items():
            if key.startswith("b_"):
                for field in ("value", "ci_low", "ci_high", "n", "responders"):
                    row[field] = "unusable_part_b_sentinel"
        self.assertEqual(self.analysis.calibrate_part_a(perturbed, self.curves), expected)
        only_a = {key: row for key, row in self.inputs.items() if key.startswith("a_")}
        self.assertEqual(self.analysis.calibrate_part_a(only_a, self.curves), expected)
        original_reader = self.analysis.read_csv

        def guarded_reader(path):
            self.assertEqual(Path(path), RESEARCH / "data" / "part_a.csv")
            return original_reader(path)

        with patch.object(self.analysis, "read_csv", side_effect=guarded_reader):
            self.assertEqual(self.analysis.calibrate_part_a(curves=self.curves), expected)

    def test_manifest_uses_explicit_posix_relative_paths(self):
        expected_paths = [ORIGINAL / "analysis.py", ORIGINAL / "data" / "tmb_orr.csv",
                          RESEARCH / "PROTOCOL.md", RESEARCH / "analysis.py", Path(__file__).resolve()]
        expected_paths += [RESEARCH / "data" / name
                           for name in ("part_a.csv", "part_b.csv", "part_a_calibration.csv", "frozen_curves.csv")]
        expected = {path.relative_to(ROOT).as_posix(): sha256_reference(path)
                    for path in sorted(expected_paths, key=lambda path: path.relative_to(ROOT).as_posix())}
        self.assertEqual(self.analysis.manifest(), expected)
        self.assertEqual(list(self.analysis.manifest()), sorted(expected))
        self.assertTrue(all("\\" not in key and not Path(key).is_absolute() for key in expected))

    def test_bad_source_rows_fail_without_implicit_zero_imputation(self):
        all_rows = [copy.deepcopy(row) for row in self.inputs.values()]
        for field, value in (("source_url", ""), ("value", "nan"),
                             ("ci_low", "99.9"), ("n", "")):
            rows = copy.deepcopy(all_rows)
            target = next(row for row in rows if row["id"] == "b_orr_positive")
            target[field] = value

            def reader(path):
                cohort = "A" if Path(path).name == "part_a.csv" else "B"
                return [row for row in rows if row["cohort"] == cohort]

            with self.subTest(field=field), patch.object(self.analysis, "read_csv", side_effect=reader):
                with self.assertRaises(ValueError):
                    self.analysis.load_inputs()

    def test_all_reported_compatibility_regions_and_shared_scale_intersections(self):
        for name, curve in self.curves.items():
            for mode in ("printed", "rounded", "subset_a"):
                actual = self.report["regions"][name][mode]
                expected = panel_reference(curve, self.inputs, mode)
                with self.subTest(curve=name, mode=mode):
                    self.assertRegion(actual["joint"], expected["joint"])
                    for cohort in ("A", "B"):
                        for group in ("positive", "negative", "joint"):
                            self.assertRegion(actual[cohort][group], expected[cohort][group])
                    for band in ("counting_assumption", "wide_assumption"):
                        self.assertRegion(actual["bands"][band], expected["bands"][band])
        expected = panel_reference(self.curves["reciprocal_holdout"], self.inputs, "median_envelope")
        actual = self.report["median_envelope"]
        for cohort in ("A", "B"):
            for group in ("positive", "negative", "joint"):
                self.assertRegion(actual[cohort][group], expected[cohort][group])
        self.assertRegion(actual["joint"], expected["joint"])
        for band in ("counting_assumption", "wide_assumption"):
            self.assertRegion(actual["bands"][band], expected["bands"][band])

    def assertScore(self, actual, expected):
        self.assertEqual(set(actual), set(expected))
        for field in ("mse", "rmse", "mae", "interval_mse"):
            self.assertAlmostEqual(actual[field], expected[field], places=14)
        self.assertEqual(actual["misses"], expected["misses"])
        for field in ("residuals", "distances"):
            self.assertEqual(len(actual[field]), len(expected[field]))
            for first, second in zip(actual[field], expected[field]):
                self.assertAlmostEqual(first, second, places=14)

    def test_every_baseline_prediction_score_and_historical_residual(self):
        grid = (("0.5", 0.5), ("1", 1.0), ("4/3", 4 / 3), ("2", 2.0), ("5", 5.0))
        self.assertEqual(self.analysis.SCALES, grid)
        rows = self.report["baseline"]
        self.assertEqual(len(rows), len(self.curves) * len(grid))
        self.assertEqual({(row["curve"], row["scale_label"]) for row in rows},
                         {(curve, label) for curve in self.curves for label, _ in grid})
        for row in rows:
            curve = self.curves[row["curve"]]
            for cohort in ("A", "B"):
                predictions = [probability_reference(curve, row["scale"] *
                               float(self.inputs[cohort.lower() + "_tmb_" + status]["value"]))
                               for status in ("positive", "negative")]
                for actual, expected in zip(row["predictions"][cohort], predictions):
                    self.assertAlmostEqual(actual, expected, places=14)
            self.assertScore(row["score"], score_reference(row["predictions"]["B"], self.inputs))
        historic = self.report["historic"]
        raw = next(row for row in self.original_rows if row["parent_type"] == "merkel")
        self.assertEqual((historic["tmb"], historic["orr"]), (float(raw["tmb"]), float(raw["orr"])))
        expected = probability_reference(self.curves["reciprocal_holdout"], float(raw["tmb"]))
        self.assertAlmostEqual(historic["prediction"], expected, places=14)
        self.assertAlmostEqual(historic["residual"], float(raw["orr"]) - expected, places=14)

    def test_every_transferred_extension_endpoint_score_and_repair_verdict(self):
        rows = self.report["extensions"]
        self.assertEqual(len(rows), 2 * len(self.analysis.SCALES))
        self.assertEqual({(row["calibration_set"], row["scale_label"]) for row in rows},
                         {(group, label) for group in ("broad", "tmb_subset")
                          for label, _ in self.analysis.SCALES})
        primary = self.curves["reciprocal_holdout"]
        median_a = float(self.inputs["a_tmb_positive"]["value"])
        median_b = float(self.inputs["b_tmb_positive"]["value"])
        negative_b = float(self.inputs["b_tmb_negative"]["value"])
        for row in rows:
            scale = row["scale"]
            orr_id = "a_orr_positive" if row["calibration_set"] == "broad" else "a_subset_orr_positive"
            outcome = self.inputs[orr_id]
            parameter = viral_reference(primary, median_a, float(outcome["value"]) / 100, scale)
            if parameter is None:
                self.assertIsNone(row["predictions"])
                self.assertIsNone(row["score"])
                self.assertFalse(row["repair"])
                continue
            self.assertAlmostEqual(row["v"], parameter, places=12)
            predictions = [probability_reference(primary, scale * median_b + parameter),
                           probability_reference(primary, scale * negative_b)]
            for actual, expected in zip(row["predictions"], predictions):
                self.assertAlmostEqual(actual, expected, places=14)
            score = score_reference(predictions, self.inputs)
            baseline_predictions = [probability_reference(primary, scale * median_b),
                                    probability_reference(primary, scale * negative_b)]
            baseline_score = score_reference(baseline_predictions, self.inputs)
            self.assertEqual(row["predictions"][1], baseline_predictions[1])
            self.assertScore(row["score"], score)
            self.assertAlmostEqual(row["delta_mse"], score["mse"] - baseline_score["mse"], places=14)
            self.assertEqual(row["repair"], score["mse"] < baseline_score["mse"] and score["misses"] == 0)
            for endpoint in ("ci_low", "ci_high"):
                value = viral_reference(primary, median_a, float(outcome[endpoint]) / 100, scale)
                expected = None if value is None else probability_reference(primary, scale * median_b + value)
                if expected is None:
                    self.assertIsNone(row["endpoint_predictions"][endpoint])
                else:
                    self.assertAlmostEqual(row["endpoint_predictions"][endpoint], expected, places=13)

    def test_generated_readme_and_saved_fingerprints_match_all_recomputed_numbers(self):
        self.assertTrue(self.analysis.check_manifest())
        self.assertTrue(self.analysis.check_readme(self.report))
        text = (RESEARCH / "README.md").read_text(encoding="utf-8")
        actual = text.split(self.analysis.BEGIN + "\n", 1)[1].split(self.analysis.END, 1)[0]
        self.assertEqual(actual, "\n" + self.analysis.render_report(self.report) + "\n")


if __name__ == "__main__":
    unittest.main()
