"""Offline standard-library checks of the NONCLINICAL Merkel TMB-bin audit.

Independent bisection, least-squares, finite-distribution and provenance checks
test the implementation. They do not establish source truth, assay calibration,
biological attribution, clinical usefulness or independent scientific review.
"""

import copy
import csv
from decimal import Decimal
import hashlib
import importlib.util
import math
from pathlib import Path
import socket
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "merkel-tmb-bin-test"
PR72 = ROOT / "research" / "merkel-viral-antigen-test"
PR67 = ROOT / "research" / "immune-escape-pancancer-test"
PROTOCOL_SHA256 = "4789e3562692bf563d0eed4754e519b5f5ae45424092516e074117fd7a1ae1d6"
COEFFICIENT_TOLERANCE = 2e-6
GRADIENT_BOUND = 1e-7


def csv_rows(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def digest_reference(path):
    canonical = Path(path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(canonical).hexdigest()


def load_analysis():
    spec = importlib.util.spec_from_file_location(
        "merkel_tmb_bin_analysis_test", RESEARCH / "analysis.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def response_reference(curve, load):
    if not load:
        return 0.0
    if curve["kind"] == "reciprocal":
        return curve["A"] * math.exp(-curve["B"] / load)
    return max(0.0, min(1.0, curve["a"] + curve["b"] * math.log(load)))


def target_reference(row, mode):
    if mode == "point":
        value = float(Decimal(row["orr_percent"]) / 100)
        return value, value
    low, high = (Decimal(row[key]) for key in ("ci_low_percent", "ci_high_percent"))
    if mode == "rounded":
        low -= Decimal(5).scaleb(low.as_tuple().exponent - 1)
        high += Decimal(5).scaleb(high.as_tuple().exponent - 1)
    return float(max(Decimal(0), low) / 100), float(min(Decimal(100), high) / 100)


def inverse_reference(curve, response, cutoff):
    """Invert numerically in log scale, independently of analytic production."""
    low, high = -40.0, 40.0
    for _ in range(150):
        midpoint = (low + high) / 2
        if response_reference(curve, cutoff * math.exp(midpoint)) < response:
            low = midpoint
        else:
            high = midpoint
    return math.exp((low + high) / 2)


def region_reference(curve, bin_name, lower, upper, cutoff=2.0):
    """Return an independent (lo, hi, left_closed, right_closed) tuple."""
    if curve["kind"] == "reciprocal":
        ceiling = curve["A"]
        if lower >= ceiling or (bin_name == "high" and upper <= 0):
            return None
        if bin_name == "low":
            return (0.0, math.inf, False, False) if lower == 0 else (
                inverse_reference(curve, lower, cutoff), math.inf, False, False)
        return (0.0, math.inf, False, False) if upper >= ceiling else (
            0.0, inverse_reference(curve, upper, cutoff), False, True)
    if bin_name == "low":
        return (0.0, math.inf, False, False) if lower == 0 else (
            inverse_reference(curve, lower, cutoff), math.inf, False, False)
    if upper == 1:
        return (0.0, math.inf, False, False)
    # The zero plateau has a finite right edge, unlike reciprocal TMB.
    if upper == 0:
        return (0.0, math.exp(-curve["a"] / curve["b"]) / cutoff, False, True)
    return (0.0, inverse_reference(curve, upper, cutoff), False, True)


def intersection_reference(regions):
    if any(region is None for region in regions):
        return None
    lower = max(region[0] for region in regions)
    upper = min(region[1] for region in regions)
    left_closed = all(region[2] for region in regions if region[0] == lower)
    right_closed = all(region[3] for region in regions if region[1] == upper)
    if lower > upper or (lower == upper and not (left_closed and right_closed)):
        return None
    return lower, upper, left_closed, right_closed


def contains_reference(region, value):
    return region is not None and (value > region[0] or (
        value == region[0] and region[2])) and (value < region[1] or (
        value == region[1] and region[3]))


def means_reference(curve, bin_name, scale, cutoff=2.0):
    threshold = response_reference(curve, cutoff * scale)
    if curve["kind"] == "reciprocal":
        return (0.0, threshold, True, False) if bin_name == "low" else (
            threshold, curve["A"], True, False)
    if bin_name == "high":
        return threshold, 1.0, True, True
    x0 = math.exp(-curve["a"] / curve["b"])
    x1 = math.exp((1.0 - curve["a"]) / curve["b"])
    if cutoff * scale <= x0:
        return 0.0, 0.0, True, True
    return 0.0, threshold, True, cutoff * scale > x1


def independent_reciprocal_fit(rows):
    """Projected ceiling and derivative-root search, with explicit bounds.

    This is independent of the original golden-section optimizer. Gradients
    and objective comparisons are bounded diagnostics, never saved fixtures.
    """
    pairs = [(float(row["tmb"]), float(row["orr"])) for row in rows]

    def profile(b):
        shapes = [math.exp(-b / tmb) for tmb, _ in pairs]
        denominator = math.fsum(z * z for z in shapes)
        a = min(1.0, max(0.0, math.fsum(z * y for z, (_, y) in zip(shapes, pairs))
                         / denominator)) if denominator else 0.0
        residuals = [a * z - y for z, (_, y) in zip(shapes, pairs)]
        gradient = -2 * math.fsum(r * a * z / tmb
                                for r, z, (tmb, _) in zip(residuals, shapes, pairs))
        return math.fsum(r * r for r in residuals), a, b, gradient

    maximum = max(tmb for tmb, _ in pairs) * 50
    grid = [0.0] + [math.exp(math.log(1e-7) + i / 240 * (
        math.log(maximum) - math.log(1e-7))) for i in range(241)]
    candidates = [profile(b) for b in grid]
    for left, right in zip(grid, grid[1:]):
        if profile(left)[3] < 0 < profile(right)[3]:
            for _ in range(100):
                midpoint = (left + right) / 2
                if profile(midpoint)[3] < 0:
                    left = midpoint
                else:
                    right = midpoint
            candidates.append(profile((left + right) / 2))
    return min(candidates), profile


class MerkelTmbBinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with patch.object(socket, "create_connection", side_effect=AssertionError("network forbidden")):
            cls.analysis = load_analysis()
            cls.rows = cls.analysis.load_inputs()
            cls.curves = cls.analysis.load_curves()
            cls.report = cls.analysis.build_report()

    def assertRegion(self, actual, expected):
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertEqual((actual.left_closed, actual.right_closed), expected[2:])
        for value, reference in zip((actual.lo, actual.hi), expected[:2]):
            if math.isinf(reference):
                self.assertEqual(value, reference)
            else:
                self.assertTrue(math.isclose(value, reference, rel_tol=1e-12, abs_tol=1e-12),
                                (value, reference))

    def test_printed_aggregate_values_and_denominator_provenance(self):
        expected = {
            "A_low": ("A", "low", "28.0", "12.1", "49.4", "36"),
            "A_high": ("A", "high", "45.5", "16.7", "76.6", "36"),
            "B_low": ("B", "low", "41.2", "24.6", "59.3", "52"),
            "B_high": ("B", "high", "50.0", "26.0", "74.0", "52"),
        }
        self.assertEqual({row["id"] for row in self.rows}, set(expected))
        self.assertEqual(len(self.rows), len(expected))
        required = ("population", "endpoint", "source_title", "doi", "source_url", "html_url",
                    "section", "paragraph", "denominator_locator", "access_date", "access_timezone",
                    "reuse_basis", "license_url", "inspection", "limitations", "notes")
        for row in self.rows:
            self.assertEqual(tuple(row[key] for key in ("part", "bin", "orr_percent", "ci_low_percent",
                                                      "ci_high_percent", "tmb_evaluable_cohort_n")),
                             expected[row["id"]])
            self.assertEqual(float(row["cutoff_nssv_mb"]), 2.0)
            self.assertEqual(row["bin_n"], "")
            self.assertEqual(row["responders"], "")
            self.assertTrue(all(row[key].strip() for key in required))
            self.assertIn("HTML", row["inspection"].upper())
            self.assertNotIn(".pdf", row["html_url"].lower())
            self.assertEqual(row["access_timezone"], "America/New_York")
            self.assertIn("pmc.ncbi.nlm.nih.gov/articles/PMC", row["source_url"])
            self.assertIn("by-nc" if row["part"] == "B" else "/by/", row["license_url"])
        self.assertEqual(self.rows, csv_rows(self.analysis.DATA))

    def test_original_curve_snapshots_and_independent_training_checks(self):
        snapshot = csv_rows(self.analysis.CURVES_PATH)
        self.assertEqual(set(self.curves), {row["name"] for row in snapshot})
        training = [row for row in csv_rows(PR67 / "data" / "tmb_orr.csv") if row["tmb"]]
        for raw in snapshot:
            curve = self.curves[raw["name"]]
            for key in ("A", "B") if curve["kind"] == "reciprocal" else ("a", "b"):
                self.assertEqual(curve[key], float(raw[key]))
            rows = [row for row in training if raw["name"].endswith("_all")
                    or row["parent_type"] != "merkel"]
            self.assertEqual(curve["n"], len(rows))
            self.assertEqual(curve["parents"], len({row["parent_type"] for row in rows}))
            if curve["kind"] == "reciprocal":
                (loss, a, b, gradient), profile = independent_reciprocal_fit(rows)
                self.assertAlmostEqual(curve["A"], a, delta=COEFFICIENT_TOLERANCE)
                self.assertAlmostEqual(curve["B"], b, delta=COEFFICIENT_TOLERANCE)
                self.assertLessEqual(abs(gradient), GRADIENT_BOUND)
                self.assertLessEqual(abs(profile(curve["B"])[3]), GRADIENT_BOUND)
                self.assertLessEqual(profile(curve["B"])[0] - loss, 1e-12)
                shapes = [math.exp(-curve["B"] / float(row["tmb"])) for row in rows]
                gradient_a = 2 * math.fsum((curve["A"] * z - float(row["orr"])) * z
                                          for row, z in zip(rows, shapes))
                self.assertLessEqual(abs(gradient_a), GRADIENT_BOUND)
            else:
                points = [(math.log(float(row["tmb"])), float(row["orr"])) for row in rows]
                denominator = math.fsum((x - xx) ** 2 for x, _ in points for xx, _ in points)
                slope = math.fsum((x - xx) * (y - yy) for x, y in points for xx, yy in points) / denominator
                intercept = math.fsum(y - slope * x for x, y in points) / len(points)
                self.assertAlmostEqual(curve["a"], intercept, delta=COEFFICIENT_TOLERANCE)
                self.assertAlmostEqual(curve["b"], slope, delta=COEFFICIENT_TOLERANCE)
            self.assertEqual(raw["pr67_code_sha256"], digest_reference(PR67 / "analysis.py"))
            self.assertEqual(raw["pr67_csv_sha256"], digest_reference(PR67 / "data" / "tmb_orr.csv"))

    def test_endpoint_aware_intersections(self):
        interval = self.analysis.Interval
        cases = [
            ([interval(0, 1, True, False), interval(1, 2, True, True)], None),
            ([interval(0, 1, True, True), interval(1, 2, True, True)], (1, 1, True, True)),
            ([interval(0, 1, True, True), interval(1, 2, False, True)], None),
            ([interval(0, 2, False, False), interval(1, 3, True, False)], (1, 2, True, False)),
            ([interval(0, math.inf, False, False), interval(0, 1, True, True)], (0, 1, False, True)),
        ]
        for regions, expected in cases:
            self.assertRegion(self.analysis.intersect(*regions), expected)
        self.assertIsNone(self.analysis.intersect(interval(0, 1, True, True), None))

    def test_reciprocal_full_range_and_unattained_endpoints(self):
        curve = {"kind": "reciprocal", "A": 0.6, "B": 2.0}
        threshold = response_reference(curve, 2)
        self.assertRegion(self.analysis.attainable(curve, "low", 1), (0, threshold, True, False))
        self.assertRegion(self.analysis.attainable(curve, "high", 1), (threshold, 0.6, True, False))
        for bin_name in ("low", "high"):
            for lower, upper in ((0.6, 0.6), (0.6, 1), (0.7, 0.8)):
                self.assertIsNone(self.analysis.scale_region(curve, bin_name, lower, upper))
        self.assertIsNone(self.analysis.scale_region(curve, "high", 0, 0))
        self.assertRegion(self.analysis.scale_region(curve, "low", 0, 0),
                          (0, math.inf, False, False))
        for bin_name, endpoint in (("low", threshold), ("high", 0.6)):
            region = self.analysis.attainable(curve, bin_name, 1)
            self.assertFalse(region.contains(endpoint))
            self.assertEqual(max(region.lo - endpoint, endpoint - region.hi, 0.0), 0.0)

    def test_clipped_log_plateau_endpoint_rules(self):
        curve = {"kind": "log", "a": 0.0, "b": 1.0}
        x0, x1 = 1.0, math.e
        for cutoff in (1.0, 2.0, 3.0):
            self.assertRegion(self.analysis.attainable(curve, "low", x0 / cutoff, cutoff),
                              (0, 0, True, True))
            self.assertRegion(self.analysis.attainable(curve, "low", x1 / cutoff, cutoff),
                              (0, 1, True, False))
            self.assertRegion(self.analysis.attainable(curve, "low", 2 * x1 / cutoff, cutoff),
                              (0, 1, True, True))
            self.assertRegion(self.analysis.scale_region(curve, "high", 0, 0, cutoff),
                              (0, x0 / cutoff, False, True))
            self.assertRegion(self.analysis.scale_region(curve, "low", 1, 1, cutoff),
                              (x1 / cutoff, math.inf, False, False))
            self.assertRegion(self.analysis.scale_region(curve, "high", 1, 1, cutoff),
                              (0, math.inf, False, False))

    def test_every_target_and_scale_region_recomputed_from_csv(self):
        for row in self.rows:
            for mode in ("interval", "point", "rounded"):
                expected_target = target_reference(row, mode)
                for actual, expected in zip(self.analysis.target(row, mode), expected_target):
                    self.assertAlmostEqual(actual, expected, delta=1e-12)
                for curve in self.curves.values():
                    expected = region_reference(curve, row["bin"], *expected_target,
                                                float(row["cutoff_nssv_mb"]))
                    actual = self.analysis.scale_region(curve, row["bin"], *expected_target,
                                                        float(row["cutoff_nssv_mb"]))
                    self.assertRegion(actual, expected)
                    if actual is not None and actual.lo > 0:
                        self.assertFalse(actual.contains(actual.lo))
                        self.assertTrue(actual.contains(actual.lo * (1 + 1e-8)))
                    if actual is not None and math.isfinite(actual.hi):
                        self.assertTrue(actual.contains(actual.hi))
                        self.assertFalse(actual.contains(actual.hi * (1 + 1e-8)))

    def test_witness_does_not_close_an_excluded_inverse_endpoint(self):
        for curve in self.curves.values():
            for row in self.rows:
                if row["bin"] != "low":
                    continue
                region = self.analysis.scale_region(curve, "low", *target_reference(row, "interval"))
                self.assertIsNotNone(region)
                self.assertFalse(region.left_closed)
                # A forward log/exp round-off can place f(2*c) one ulp above
                # the target. It must not create a forbidden T=2 low-bin witness.
                self.assertIsNone(self.analysis.witness(curve, row, region.lo))

    def test_printed_decimal_precision_and_clipped_rounding(self):
        for row in (
                {"orr_percent": "40.0", "ci_low_percent": "24.60", "ci_high_percent": "59.3"},
                {"orr_percent": "50.0", "ci_low_percent": "0.00", "ci_high_percent": "100.0"}):
            for mode in ("interval", "point", "rounded"):
                self.assertEqual(self.analysis.target(row, mode), target_reference(row, mode))

    def test_source_loader_rejects_unprinted_counts_and_broken_rows(self):
        mutations = []
        for field, value in (("bin_n", "25"), ("responders", "7"), ("source_url", ""),
                             ("orr_percent", "-1"), ("ci_low_percent", "nan"),
                             ("cutoff_nssv_mb", "3"), ("tmb_evaluable_cohort_n", "999")):
            changed = copy.deepcopy(self.rows)
            changed[0][field] = value
            mutations.append(changed)
        mutations.extend((self.rows[:-1], self.rows + [self.rows[0]]))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.csv"
            for rows in mutations:
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(self.rows[0]))
                    writer.writeheader()
                    writer.writerows(rows)
                with self.assertRaises(ValueError):
                    self.analysis.load_inputs(path)

    def test_every_reported_region_and_band_intersection(self):
        self.assertEqual(self.report["inputs"], self.rows)
        self.assertEqual(self.report["curves"], self.curves)
        band_references = {"counting_assumption": (1, 2, True, True),
                           "wide_assumption": (0.5, 5, True, True)}
        for mode in ("interval", "point", "rounded"):
            for name, curve in self.curves.items():
                actual = self.report["modes"][mode][name]
                references = {row["id"]: region_reference(curve, row["bin"],
                    *target_reference(row, mode), float(row["cutoff_nssv_mb"])) for row in self.rows}
                self.assertEqual(set(actual["bins"]), set(references))
                for row_id, expected in references.items():
                    self.assertRegion(actual["bins"][row_id], expected)
                for part in ("A", "B"):
                    self.assertRegion(actual["cohorts"][part], intersection_reference(
                        [references[row["id"]] for row in self.rows if row["part"] == part]))
                joint = intersection_reference(list(references.values()))
                self.assertRegion(actual["joint"], joint)
                self.assertEqual(set(actual["bands"]), set(band_references))
                for band, limits in band_references.items():
                    self.assertRegion(actual["bands"][band], intersection_reference([joint, limits]))

    def test_every_grid_range_score_and_witness(self):
        scales = [("0.5", 0.5), ("1", 1.0), ("4/3", 4 / 3), ("2", 2.0), ("5", 5.0)]
        for name, curve in self.curves.items():
            grid = self.report["grid"][name]
            self.assertEqual([(row["scale_label"], row["scale"]) for row in grid], scales)
            for entry in grid:
                interval_passes = point_passes = 0
                gaps = []
                for row in self.rows:
                    row_id, cutoff = row["id"], float(row["cutoff_nssv_mb"])
                    expected = means_reference(curve, row["bin"], entry["scale"], cutoff)
                    self.assertRegion(entry["ranges"][row_id], expected)
                    self.assertRegion(self.analysis.attainable(curve, row["bin"], entry["scale"], cutoff), expected)
                    observed = float(Decimal(row["orr_percent"]) / 100)
                    gap = 100 * max(expected[0] - observed, observed - expected[1], 0)
                    gaps.append(gap)
                    self.assertAlmostEqual(entry["point_gaps_pp"][row_id], gap, delta=1e-10)
                    interval_region = region_reference(curve, row["bin"], *target_reference(row, "interval"), cutoff)
                    compatible = contains_reference(interval_region, entry["scale"])
                    interval_passes += compatible
                    point_passes += contains_reference(expected, observed)
                    witness = entry["witnesses"][row_id]
                    if compatible:
                        self.assertIsNotNone(witness)
                        mean, tmb = witness["mean"], witness["tmb"]
                        lower, upper = target_reference(row, "interval")
                        self.assertLessEqual(lower, mean)
                        self.assertLessEqual(mean, upper)
                        self.assertTrue(math.isfinite(tmb))
                        self.assertGreaterEqual(tmb, 0)
                        self.assertLess(tmb, cutoff) if row["bin"] == "low" else self.assertGreaterEqual(tmb, cutoff)
                        self.assertAlmostEqual(response_reference(curve, entry["scale"] * tmb), mean, delta=1e-12)
                    else:
                        self.assertIsNone(witness)
                self.assertEqual(entry["interval_passes"], interval_passes)
                self.assertEqual(entry["point_passes"], point_passes)
                self.assertAlmostEqual(entry["rms_gap_pp"], math.sqrt(math.fsum(g * g for g in gaps) / len(gaps)),
                                       delta=1e-10)

    def test_arbitrary_finite_distributions_respect_the_sharp_ranges(self):
        for curve in self.curves.values():
            for scale in (0.5, 1.0, 2.0, 5.0, 100.0):
                for bin_name, supports in (("low", [(0.0,), (0.2, 0.8, 1.7), (1.999,)]),
                                           ("high", [(2.0,), (2.0, 5.0, 500.0), (10000.0,)])):
                    actual = self.analysis.attainable(curve, bin_name, scale)
                    for support in supports:
                        mean = math.fsum(response_reference(curve, scale * tmb) for tmb in support) / len(support)
                        self.assertTrue(actual.contains(mean), (curve, scale, bin_name, support, mean, actual))

    def test_part_b_perturbation_does_not_change_frozen_curves_or_part_a(self):
        changed = copy.deepcopy(self.rows)
        for row in changed:
            if row["part"] == "B":
                row.update(orr_percent="1.0", ci_low_percent="0.5", ci_high_percent="2.0")
        altered = self.analysis.build_report(rows=changed)
        self.assertEqual(altered["curves"], self.report["curves"])
        for mode in ("interval", "point", "rounded"):
            for name in self.curves:
                self.assertEqual(altered["modes"][mode][name]["cohorts"]["A"],
                                 self.report["modes"][mode][name]["cohorts"]["A"])
        self.assertNotEqual(altered["modes"], self.report["modes"])

    def test_protocol_fingerprint_newline_normalization_and_posix_keys(self):
        self.assertEqual(digest_reference(RESEARCH / "PROTOCOL.md"), PROTOCOL_SHA256)
        readme = (RESEARCH / "README.md").read_text(encoding="utf-8")
        preserved = readme.split("<!-- BEGIN ORIGINAL PROTOCOL COMMIT -->", 1)[1].split(
            "<!-- END ORIGINAL PROTOCOL COMMIT -->", 1)[0]
        payload = preserved.split("```text\n", 1)[1].split("```", 1)[0].encode("utf-8")
        header = ("commit " + str(len(payload)) + "\0").encode("ascii")
        self.assertEqual(hashlib.sha1(header + payload, usedforsecurity=False).hexdigest(),
                         self.analysis.LOCAL_PROTOCOL_COMMIT)
        hashes = self.analysis.file_hashes()
        self.assertEqual(hashes, self.report["hashes"])
        for key, value in hashes.items():
            self.assertNotIn("\\", key)
            self.assertEqual(Path(key).as_posix(), key)
            self.assertFalse(Path(key).is_absolute())
            self.assertEqual(value, digest_reference(ROOT / Path(key)))
        for path in (RESEARCH / "PROTOCOL.md", RESEARCH / "analysis.py", self.analysis.DATA,
                     self.analysis.CURVES_PATH, Path(__file__).resolve()):
            self.assertIn(path.relative_to(ROOT).as_posix(), hashes)
        self.assertNotIn((RESEARCH / "README.md").relative_to(ROOT).as_posix(), hashes)
        text = "TMB <2\nπ = nonclinical\n"
        with tempfile.TemporaryDirectory() as directory:
            for index, newline in enumerate(("\n", "\r\n", "\r")):
                path = Path(directory) / (str(index) + ".txt")
                path.write_bytes(text.replace("\n", newline).encode("utf-8"))
                self.assertEqual(self.analysis.normalized_bytes(path), text.encode("utf-8"))
                self.assertEqual(self.analysis.digest(path), hashlib.sha256(text.encode("utf-8")).hexdigest())

    def test_readme_numeric_block_is_recomputed_without_network(self):
        with patch.object(socket, "create_connection", side_effect=AssertionError("network forbidden")), \
             patch.object(socket, "socket", side_effect=AssertionError("network forbidden")):
            regenerated = self.analysis.render_report(self.analysis.build_report())
        readme = (RESEARCH / "README.md").read_text(encoding="utf-8")
        self.assertEqual(readme.count(self.analysis.BEGIN), 1)
        self.assertEqual(readme.count(self.analysis.END), 1)
        stored = readme.split(self.analysis.BEGIN, 1)[1].split(self.analysis.END, 1)[0]
        self.assertEqual(stored.strip(), regenerated.strip())


if __name__ == "__main__":
    unittest.main()
