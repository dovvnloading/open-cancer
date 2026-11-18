"""Offline stdlib checks for the stem-cell-divisions aggregate reanalysis.

Nonclinical computational research; no medical advice. The complete empirical
report is recomputed from the cited CSVs. Independent QR/pairwise regressions,
Decimal survivor probabilities, rational identities and fold perturbations
check material risks beyond comparing a cached output with its producer.
Finite checks corroborate the implementations, not the general proofs.
"""

from __future__ import annotations

from collections import defaultdict
import csv
from datetime import date
from decimal import Decimal, localcontext
from functools import lru_cache
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import random
import statistics
import sys
import tempfile
import unittest
from unittest.mock import patch
from urllib.parse import urlsplit


RESEARCH = Path(__file__).resolve().parents[1] / "research" / "stem-cell-divisions-risk"
_spec = importlib.util.spec_from_file_location("research_scd_analysis", RESEARCH / "analysis.py")
A = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = A
_spec.loader.exec_module(A)
V = A.V

NUMERIC = ("risk", "stem_cells", "divisions_per_year", "divisions_per_lifetime", "total_divisions")
PROVENANCE = ("source", "source_url", "source_locator", "source_version", "access_date",
              "reuse_basis", "reuse_terms_url")
HIGH_RISK = {"colorectal_fap", "colorectal_lynch", "duodenum_fap", "head_neck_hpv16",
             "liver_hcv", "lung_smokers"}
VERSION_CHANGES = {"glioblastoma", "medulloblastoma", "ovarian_germ"}


def raw_rows(filename):
    with (RESEARCH / "data" / filename).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def numeric_rows(filename):
    """Separate parser, deliberately independent of analysis.load_rows."""
    rows = raw_rows(filename)
    for row in rows:
        for field in NUMERIC:
            row[field] = float(Decimal(row[field]))
        row["high_risk_subgroup"] = int(row["high_risk_subgroup"])
    return rows


def primary(rows):
    return [row for row in rows if row["tissue_id"] not in HIGH_RISK]


def exposure(row, kind):
    if kind == "effective":
        return row["total_divisions"] / row["stem_cells"]
    return row["divisions_per_lifetime"] + math.log(row["stem_cells"], 2)


def log_h(row):
    return math.log10(-math.log1p(-row["risk"]))


def pairwise_slope(xs, ys):
    """OLS via pair differences, independently of centered normal equations."""
    pairs = tuple(itertools.combinations(range(len(xs)), 2))
    numerator = math.fsum((xs[j] - xs[i]) * (ys[j] - ys[i]) for i, j in pairs)
    denominator = math.fsum((xs[j] - xs[i]) ** 2 for i, j in pairs)
    if not denominator:
        raise ValueError("constant predictor")
    return numerator / denominator


def independent_qr(design, response):
    """Modified Gram--Schmidt/back-substitution, not the analysis's 2x2 solve."""
    columns = [list(column) for column in zip(*design)]
    count = len(columns)
    q = []
    upper = [[0.0] * count for _ in columns]
    for j, column in enumerate(columns):
        work = column[:]
        for i in range(j):
            upper[i][j] = math.fsum(a * b for a, b in zip(q[i], work))
            work = [a - upper[i][j] * b for a, b in zip(work, q[i])]
        upper[j][j] = math.sqrt(math.fsum(a * a for a in work))
        if upper[j][j] < 1e-12:
            raise ValueError("dependent design")
        q.append([a / upper[j][j] for a in work])
    target = [math.fsum(a * b for a, b in zip(column, response)) for column in q]
    coefficients = [0.0] * count
    for j in range(count - 1, -1, -1):
        coefficients[j] = (target[j] - math.fsum(upper[j][i] * coefficients[i]
                                                for i in range(j + 1, count))) / upper[j][j]
    return coefficients


@lru_cache(maxsize=None)
def independent_log_phi(k, z):
    """Direct fixed positive series; data/grid z<=17 makes 160 terms ample."""
    if z == 0:
        return 0.0
    if not 0 < z < 20:
        raise AssertionError("independent series only justified on the tested grid")
    value = math.fsum(z ** n / math.factorial(k + n) for n in range(160)) * math.factorial(k)
    return math.log(value)


def independent_shape(row, k, kind, gamma=0.0):
    v = exposure(row, kind)
    # A different log base/formula avoids sharing analysis.log_shape.
    return (math.log10(row["stem_cells"]) + k * math.log10(v)
            - math.log10(math.factorial(k)) + independent_log_phi(k, gamma * v) / math.log(10))


def independent_fit(rows, family, kind="effective", k=None):
    y = [log_h(row) for row in rows]
    if family == "descriptive":
        x = [math.log10(row["total_divisions"]) for row in rows]
        b = pairwise_slope(x, y)
        a = statistics.fmean(y) - b * statistics.fmean(x)
        return (a, b)
    if family == "separable":
        design = [(1, math.log10(row["stem_cells"]), math.log10(exposure(row, kind))) for row in rows]
        return independent_qr(design, y)
    grid = (0, 0.25, 0.5, 1, 2, 4, 8) if family == "expansion" else (0,)
    candidates = []
    for c in grid:
        gamma = c / max(exposure(row, kind) for row in rows)
        shape = [independent_shape(row, k, kind, gamma) for row in rows]
        a = math.fsum(actual - expected for actual, expected in zip(y, shape)) / len(y)
        sse = math.fsum((actual - a - expected) ** 2 for actual, expected in zip(y, shape))
        candidates.append((sse, c, a, gamma))
    _, c, a, gamma = min(candidates)
    return (a, gamma, c)


def independent_prediction(parameters, row, family, kind="effective", k=None):
    if family == "descriptive":
        return parameters[0] + parameters[1] * math.log10(row["total_divisions"])
    if family == "separable":
        return (parameters[0] + parameters[1] * math.log10(row["stem_cells"])
                + parameters[2] * math.log10(exposure(row, kind)))
    return parameters[0] + independent_shape(row, k, kind, parameters[1])


def decimal_psi(k, x):
    """Independent high-precision *survival* computation, not float tail code."""
    survival = (-x).exp() * sum((x ** j / math.factorial(j) for j in range(k)), Decimal(0))
    return -survival.ln()


class StemCellDivisionsRiskTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = numeric_rows("tv2015.csv")
        cls.revised = numeric_rows("little2016.csv")
        cls.corrected = numeric_rows("tv2015_nunney2020.csv")
        cls.results = A.analyze()
        cls.saved = json.loads((RESEARCH / "analysis_results.json").read_text(encoding="utf-8"))
        cls.panels = {
            "primary": (primary(cls.original), "effective"),
            "all_finite": ([r for r in cls.original if r["risk"] < 1], "effective"),
            "primary_developmental": (primary(cls.original), "developmental"),
            "later_effective": (primary(cls.revised), "effective"),
            "later_developmental": (primary(cls.revised), "developmental"),
            "corrected_effective": (primary(cls.corrected), "effective"),
        }

    def assert_close(self, actual, expected, *, rel=3e-10, absolute=0.0, message=None):
        self.assertTrue(math.isclose(actual, expected, rel_tol=rel, abs_tol=absolute),
                        message or f"{actual!r} != {expected!r}")

    def assert_tree_equal(self, actual, expected, path="results"):
        if isinstance(expected, dict):
            self.assertIsInstance(actual, dict, path)
            self.assertEqual(set(actual), set(expected), path)
            for key in expected:
                self.assert_tree_equal(actual[key], expected[key], f"{path}.{key}")
        elif isinstance(expected, list):
            self.assertIsInstance(actual, list, path)
            self.assertEqual(len(actual), len(expected), path)
            for i, (left, right) in enumerate(zip(actual, expected)):
                self.assert_tree_equal(left, right, f"{path}[{i}]")
        elif isinstance(expected, float):
            self.assertIsInstance(actual, (float, int), path)
            self.assertTrue(math.isfinite(actual), path)
            # Only cancellation diagnostics may have an absolute tolerance.
            # Tiny rates/interval endpoints retain a relative error criterion.
            cancellation = ("max_exact_logH_discrepancy" in path
                            or "max_relative_exact_vs_rare_rate" in path)
            self.assert_close(actual, expected, rel=5e-11,
                              absolute=2e-13 if cancellation else 0, message=path)
        else:
            self.assertEqual(actual, expected, path)

    def test_01_all_reported_results_and_readme_are_recomputed(self):
        # Archive the complete published summary tree. Auxiliary per-fold
        # predictions and per-tissue rates remain in the fresh full results
        # used by the independent checks below, avoiding a duplicated payload.
        self.assert_tree_equal(A.compact(A.artifact_results(self.results)), self.saved)
        text = (RESEARCH / "README.md").read_text(encoding="utf-8")
        begin, end = "<!-- BEGIN RECOMPUTED RESULTS -->", "<!-- END RECOMPUTED RESULTS -->"
        self.assertEqual(text.count(begin), 1)
        self.assertEqual(text.count(end), 1)
        block = text[text.index(begin):text.index(end) + len(end)]
        # Compare both the archived numbers and a newly recomputed rendering.
        self.assertEqual(block, A.render_results(self.saved))
        self.assertEqual(block, A.render_results(A.compact(self.results)))
        self.assertIn("nonclinical", text.lower())
        self.assertIn("no medical advice", text.lower())
        self.assertEqual(self.saved["preregistration_commit"], "f5c1f4e87bfd03c85eb015fa14e72a28d33844ad")
        self.assertEqual(self.saved["prospective_amendment_commit"], "c7d255d7bfafe0ac720e1752a4b55058d136ecd3")

    def test_02_csv_provenance_selection_and_anatomical_groups(self):
        expected_groups = {
            "hematopoietic", "skin", "colorectal", "small_intestine", "esophagus", "gallbladder",
            "brain", "head_neck", "liver", "lung", "bone", "ovary", "pancreas", "testis", "thyroid",
        }
        for filename in ("tv2015.csv", "little2016.csv", "tv2015_nunney2020.csv"):
            rows = raw_rows(filename)
            self.assertEqual(len(rows), 31)
            self.assertEqual(len({r["tissue_id"] for r in rows}), 31)
            self.assertEqual({r["tissue_group"] for r in rows}, expected_groups)
            self.assertEqual({r["tissue_id"] for r in rows if r["high_risk_subgroup"] == "1"}, HIGH_RISK)
            for row in rows:
                for field in PROVENANCE:
                    self.assertTrue(row[field].strip(), (filename, row["tissue_id"], field))
                for field in ("source_url", "reuse_terms_url", "original_source_url"):
                    self.assertEqual(urlsplit(row[field]).scheme, "https")
                date.fromisoformat(row["access_date"])
                self.assertIn("CC0", row["reuse_basis"])
                self.assertRegex(row["source_locator"], r"row \d+")
                self.assertIn(row["high_risk_subgroup"], {"0", "1"})
                self.assertTrue(0 < Decimal(row["risk"]) <= 1)
                self.assertGreater(Decimal(row["stem_cells"]), 0)
                self.assertGreater(Decimal(row["total_divisions"]), 0)
                self.assertGreaterEqual(Decimal(row["divisions_per_year"]), 0)
                self.assertGreaterEqual(Decimal(row["divisions_per_lifetime"]), 0)
        self.assertEqual(len(primary(self.original)), 25)
        self.assertEqual(len({r["tissue_group"] for r in primary(self.original)}), 15)
        by_id = {r["tissue_id"]: r for r in self.original}
        self.assertEqual(by_id["duodenum"]["tissue_group"], by_id["small_intestine"]["tissue_group"])
        self.assertEqual(by_id["aml"]["tissue_group"], by_id["cll"]["tissue_group"])
        self.assertIn("myeloid", by_id["aml"]["tissue"].lower())
        self.assertIn("lymphocytic", by_id["cll"]["tissue"].lower())
        # Visually audited original values where circulating transcription
        # shortcuts are known to disagree with the cited Table S1.
        self.assertEqual(by_id["esophagus_squamous"]["stem_cells"], 864000)
        self.assertEqual(by_id["thyroid_papillary_follicular"]["total_divisions"], 585000000)
        self.assertEqual(by_id["thyroid_medullary"]["total_divisions"], 58500000)

    def test_03_source_versions_and_exact_decimal_correction(self):
        original, revised, corrected = map(raw_rows, ("tv2015.csv", "little2016.csv", "tv2015_nunney2020.csv"))
        changed = []
        for left, right, fixed in zip(original, revised, corrected):
            self.assertEqual(left["tissue_id"], right["tissue_id"])
            self.assertEqual(left["tissue_id"], fixed["tissue_id"])
            for field in NUMERIC:
                a, b = Decimal(left[field]), Decimal(right[field])
                if a != b:
                    changed.append((left["tissue_id"], field))
                if field != "total_divisions":
                    self.assertEqual(Decimal(fixed[field]), a)
            n, lifetime = Decimal(left["stem_cells"]), Decimal(left["divisions_per_lifetime"])
            self.assertEqual(Decimal(fixed["total_divisions"]), n * (1 + lifetime) - 1)
            self.assertIn("2020.2291", fixed["source_url"] + fixed["original_source_url"]
                          + fixed["source_locator"] + fixed["notes"])
        self.assertEqual(len(original) * len(NUMERIC), 155)
        self.assertEqual(set(changed), {(t, f) for t in VERSION_CHANGES
                                      for f in ("divisions_per_year", "divisions_per_lifetime")})
        self.assertEqual(len(changed), 6)
        for row in original:
            if row["tissue_id"] in VERSION_CHANGES:
                self.assertEqual(Decimal(row["divisions_per_lifetime"]), 0)
        for row in revised:
            if row["tissue_id"] in VERSION_CHANGES:
                self.assertEqual(Decimal(row["divisions_per_year"]), Decimal("0.025"))
                self.assertEqual(Decimal(row["divisions_per_lifetime"]), 2)
        self.assertNotIn("expansion", self.results["corrected_effective"])
        self.assertNotIn("stress_frontier", self.results["primary_developmental"])
        self.assertNotIn("stress_frontier", self.results["later_developmental"])

    def test_04_normalized_hashes_preserve_every_other_byte(self):
        expected_paths = {"PROTOCOL.md", "AMENDMENTS.md"}
        expected_paths.update(str(p.relative_to(RESEARCH)) for p in (RESEARCH / "data").glob("*.csv"))
        self.assertEqual(set(self.results["hashes"]), expected_paths)
        for relative, digest in self.results["hashes"].items():
            payload = (RESEARCH / relative).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
            self.assertEqual(digest, hashlib.sha256(payload).hexdigest())
            self.assertEqual(digest, V.normalized_sha256(RESEARCH / relative))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "newlines.txt"
            reference = "a,b\né,2\nlast".encode("utf-8")
            expected = hashlib.sha256(reference).hexdigest()
            for ending in (b"\n", b"\r\n", b"\r"):
                path.write_bytes(reference.replace(b"\n", ending))
                self.assertEqual(A.normalized_hash(path), expected)
                self.assertEqual(V.normalized_sha256(path), expected)
            path.write_bytes(reference + b"\n")
            self.assertNotEqual(A.normalized_hash(path), expected)

    def test_05_independent_regressions_and_all_fixed_model_predictions(self):
        historical = self.results["historical_all"]
        x = [math.log10(r["total_divisions"]) for r in self.original]
        y = [math.log10(r["risk"]) for r in self.original]
        b = pairwise_slope(x, y)
        self.assert_close(historical["slope"], b)
        self.assert_close(historical["r_squared"], pairwise_slope(x, y) * pairwise_slope(y, x))
        for name, (rows, kind) in self.panels.items():
            panel = self.results[name]
            models = [("descriptive", None, panel["descriptive"]),
                      ("separable", None, panel["separable"])]
            models.extend(("sequential", k, value) for k, value in enumerate(panel["sequential"], 1))
            models.extend(("expansion", value["fit"]["k"], value) for value in panel.get("expansion", []))
            for family, k, result in models:
                reference = independent_fit(rows, family, kind, k)
                predictions = [independent_prediction(reference, r, family, kind, k) for r in rows]
                for row, predicted in zip(rows, predictions):
                    self.assert_close(A.predict(result["fit"], row), predicted, absolute=1e-11,
                                      message=(name, family, k, row["tissue_id"]))
                expected_rmse = math.sqrt(statistics.fmean((log_h(row) - predicted) ** 2
                                                          for row, predicted in zip(rows, predictions)))
                self.assert_close(result["training_rmse"], expected_rmse)
                if family in ("sequential", "expansion"):
                    self.assertEqual(result["fit"]["grid_c"], reference[2])
                if family == "sequential":
                    xs = [math.log10(r["total_divisions"]) for r in rows]
                    shape_slope = (pairwise_slope(xs, [math.log10(r["stem_cells"]) for r in rows])
                                   + k * pairwise_slope(xs, [math.log10(exposure(r, kind)) for r in rows]))
                    self.assert_close(result["predicted_logH_slope"], shape_slope)

    def test_06_independent_leave_group_and_row_out_scores(self):
        for name, (rows, kind) in self.panels.items():
            panel = self.results[name]
            models = [("descriptive", None, panel["descriptive"]),
                      ("separable", None, panel["separable"])]
            models.extend(("sequential", k, value) for k, value in enumerate(panel["sequential"], 1))
            models.extend(("expansion", value["fit"]["k"], value) for value in panel.get("expansion", []))
            for family, k, result in models:
                for key, cv_name in (("tissue_group", "leave_group_out"), ("tissue_id", "leave_row_out")):
                    cv = result[cv_name]
                    actual = {p["tissue_id"]: p for p in cv["predictions"]}
                    squared = defaultdict(list)
                    self.assertEqual(set(actual), {r["tissue_id"] for r in rows})
                    for fold in sorted({r[key] for r in rows}):
                        training = [r for r in rows if r[key] != fold]
                        parameters = independent_fit(training, family, kind, k)
                        for row in (r for r in rows if r[key] == fold):
                            prediction = independent_prediction(parameters, row, family, kind, k)
                            residual = log_h(row) - prediction
                            measured = actual[row["tissue_id"]]
                            self.assert_close(measured["prediction"], prediction, absolute=2e-11,
                                              message=(name, family, k, cv_name, row["tissue_id"]))
                            self.assert_close(measured["residual"], residual, absolute=2e-11)
                            self.assert_close(measured["amplitude_ratio"], 10 ** residual)
                            if family in ("sequential", "expansion"):
                                self.assertEqual(measured["training_grid_c"], parameters[2])
                            squared[fold].append(residual ** 2)
                    self.assertEqual(cv["fold_count"], len(squared))
                    group_rmse = math.sqrt(statistics.fmean(statistics.fmean(v) for v in squared.values()))
                    row_rmse = math.sqrt(statistics.fmean(v for values in squared.values() for v in values))
                    self.assert_close(cv["group_balanced_rmse"], group_rmse)
                    self.assert_close(cv["row_rmse"], row_rmse)

    def test_07_training_grid_cannot_use_held_out_outcomes_or_clocks(self):
        rows = primary(self.original)
        target = "colorectal"  # Contains the largest observed clock.
        fit_function = A.fit

        def captured(data):
            calls = []

            def checked(training, family, kind="effective", k=None):
                model = fit_function(training, family, kind, k)
                calls.append((set(r["tissue_id"] for r in training), model))
                self.assert_close(model["gamma"], model["grid_c"] / max(exposure(r, kind) for r in training))
                return model

            with patch.object(A, "fit", side_effect=checked):
                cv = A.cross_validate(data, "expansion", k=3)
            groups = sorted({r["tissue_group"] for r in data})
            for group, (training_ids, _) in zip(groups, calls):
                self.assertEqual(training_ids, {r["tissue_id"] for r in data if r["tissue_group"] != group})
            return cv, dict(zip(groups, (model for _, model in calls)))

        base, base_models = captured(rows)
        altered_risk = [dict(row, risk=0.8) if row["tissue_group"] == target else dict(row) for row in rows]
        changed, changed_models = captured(altered_risk)
        self.assertEqual(base_models[target], changed_models[target])
        self.assertEqual([p["prediction"] for p in base["predictions"] if p["group"] == target],
                         [p["prediction"] for p in changed["predictions"] if p["group"] == target])
        altered_clock = [dict(row, total_divisions=row["total_divisions"] * 1000)
                         if row["tissue_group"] == target else dict(row) for row in rows]
        _, clock_models = captured(altered_clock)
        self.assertEqual(base_models[target], clock_models[target])
        # The test is a real perturbation: the largest held-out clock exceeds
        # every training clock by orders of magnitude, but gamma is unchanged.
        self.assertGreater(max(exposure(r, "effective") for r in altered_clock),
                           100 * max(exposure(r, "effective") for r in rows))

    def test_08_group_bootstrap_is_reconstructed_from_csv_rows(self):
        for rows, recorded in ((self.original, self.results["historical_all"]["bootstrap"]),
                               (primary(self.original), self.results["primary"]["slope_bootstrap"])):
            rng = random.Random(1260825)
            groups = sorted({r["tissue_group"] for r in rows})
            members = {g: [r for r in rows if r["tissue_group"] == g] for g in groups}
            slopes, skipped = [], 0
            for _ in range(2000):
                sampled = [rng.choice(groups) for _ in range(len(groups))]
                selected = [r for g in sampled for r in members[g]]
                x = [math.log10(r["total_divisions"]) for r in selected]
                y = [math.log10(r["risk"]) for r in selected]
                try:
                    slopes.append(pairwise_slope(x, y))
                except ValueError:
                    skipped += 1
            slopes.sort()
            for q, key in ((0.025, "lower"), (0.975, "upper")):
                position = q * (len(slopes) - 1)
                lo = int(position)
                fraction = position - lo
                value = (1 - fraction) * slopes[lo] + fraction * slopes[math.ceil(position)]
                self.assert_close(recorded[key], value)
            self.assertEqual(recorded["draws"], 2000)
            self.assertEqual(recorded["seed"], 1260825)
            self.assertEqual(recorded["skipped"], skipped)

    def test_09_decimal_survivor_checks_all_empirical_point_inverses(self):
        with localcontext() as ctx:
            ctx.prec = 90
            for name, (rows, kind) in self.panels.items():
                by_id = {r["tissue_id"]: r for r in rows}
                for certificate in self.results[name]["exact_frontier"]:
                    k = certificate["k"]
                    for point in certificate["rates"]:
                        row = by_id[point["tissue_id"]]
                        risk, n = Decimal(str(row["risk"])), Decimal(str(row["stem_cells"]))
                        expected_h = -(1 - risk).ln() / n
                        x = Decimal(str(point["rate"])) * Decimal(str(exposure(row, kind)))
                        evaluated_h = decimal_psi(k, x)
                        self.assert_close(float(evaluated_h / expected_h), 1.0, rel=2e-12,
                                          message=(name, k, point["tissue_id"]))
                    intervals = [point["interval"] for point in certificate["rates"]]
                    lower, upper = max(v[0] for v in intervals), min(v[1] for v in intervals)
                    self.assert_close(certificate["rho"] ** 2,
                                      max(1, lower / upper))
                    self.assert_close(certificate["baseline_interval"][0],
                                      certificate["baseline_interval"][1])
                    self.assertEqual(certificate["lower_constraint"],
                                     [p["tissue_id"] for p in certificate["rates"] if p["interval"][0] == lower])
                    self.assertEqual(certificate["upper_constraint"],
                                     [p["tissue_id"] for p in certificate["rates"] if p["interval"][1] == upper])
                    self.assertEqual(certificate["compatible_allowances"],
                                     [rho for rho in (1, 2, 5, 10, 100) if rho >= certificate["rho"]])
                if kind == "effective":
                    for point, stress in zip(self.results[name]["exact_frontier"],
                                             self.results[name]["stress_frontier"]):
                        self.assertLessEqual(stress["rho"], point["rho"] * (1 + 1e-13))
                        intervals = [r["interval"] for r in stress["rates"]]
                        lower, upper = max(v[0] for v in intervals), min(v[1] for v in intervals)
                        self.assert_close(stress["rho"], math.sqrt(max(1, lower / upper)))
                        self.assertLessEqual(lower / stress["rho"], stress["rho"] * upper * (1 + 1e-13))
                        for interval in stress["rates"]:
                            row = by_id[interval["tissue_id"]]
                            r, n, d = row["risk"], row["stem_cells"], row["total_divisions"]
                            corners = ((r / 1.25, n / 2, 2 * d),
                                       (min(1 - 1e-12, 1.25 * r), 2 * n, d / 2))
                            for u, corner in zip(interval["interval"], corners):
                                cr, cn, cd = map(Decimal.from_float, corner)
                                expected_h = -(1 - cr).ln() / cn
                                x = Decimal(str(u)) * cd / cn
                                self.assert_close(float(decimal_psi(stress["k"], x) / expected_h),
                                                  1, rel=2e-12,
                                                  message=(name, stress["k"], interval["tissue_id"], corner))
                        if stress["rho"] > 1:
                            too_small = stress["rho"] * (1 - 1e-7)
                            self.assertGreater(lower / too_small, too_small * upper)

    def test_10_independent_synthetic_math_checks_and_domain_boundaries(self):
        checks = V.exact_checks()
        self.assertEqual(checks["fraction_algebra_checks"], 122)
        self.assertEqual(checks["independent_decimal_inverse_cases"], 40)
        self.assertEqual(checks["elasticity_cases"], 72)
        self.assertEqual(checks["coupled_rectangle_cases"], 12)
        self.assertEqual(checks["same_input_ratio_cases"], 12)
        self.assertLess(checks["max_relative_psi_error"], 2e-13)
        self.assertLess(checks["max_relative_inverse_error"], 5e-14)
        self.assertLess(checks["max_phi_log_error"], 2e-13)
        for bad in (0, 1, -0.1, 1.1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                V.rate(2, bad, 100, 1000)
            with self.assertRaises(ValueError):
                A.hazard({"risk": bad})
        for bad_k in (0, -1, 1.5, True):
            with self.assertRaises(ValueError):
                V.inverse_psi(bad_k, 0.01)
        with self.assertRaises(ValueError):
            A.simple_fit([1, 1, 1], [2, 3, 4])
        with self.assertRaises(ValueError):
            A.separated_fit([self.original[0]] * 4, "effective")
        endpoint = next(r for r in self.original if r["tissue_id"] == "colorectal_fap")
        self.assertEqual(endpoint["risk"], 1)
        self.assertEqual(self.results["historical_all"]["n"], 31)
        self.assertEqual(self.results["all_finite"]["n"], 30)
        for name in self.panels:
            self.assertNotIn("colorectal_fap", self.results[name]["ids"])

    def test_11_every_identical_input_pair_and_hazard_rate_contrast(self):
        expected = {}
        for a, b in itertools.combinations(self.original, 2):
            if (a["stem_cells"], a["total_divisions"]) != (b["stem_cells"], b["total_divisions"]):
                continue
            if a["risk"] == b["risk"]:
                continue
            low, high = sorted((a, b), key=lambda r: (r["risk"], r["tissue_id"]))
            expected[(low["tissue_id"], high["tissue_id"])] = (low, high)
        actual = {(p["lower"], p["higher"]): p for p in self.results["duplicate_contrasts"]}
        point_rates = {c["k"]: {r["tissue_id"]: r["rate"] for r in c["rates"]}
                       for c in self.results["all_finite"]["exact_frontier"]}
        self.assertEqual(set(actual), set(expected))
        for pair, (low, high) in expected.items():
            contrast = actual[pair]
            self.assert_close(contrast["risk_ratio"], high["risk"] / low["risk"])
            self.assertEqual(low["tissue_group"], high["tissue_group"])
            if high["risk"] == 1:
                self.assertIsNone(contrast["hazard_ratio"])
                self.assertIsNone(contrast["rate_ratios"])
                continue
            ratio = math.log1p(-high["risk"]) / math.log1p(-low["risk"])
            self.assert_close(contrast["hazard_ratio"], ratio)
            self.assertEqual(set(contrast["rate_ratios"]), {str(k) for k in range(1, 9)})
            for k in range(1, 9):
                u_ratio = contrast["rate_ratios"][str(k)]
                self.assert_close(u_ratio, point_rates[k][high["tissue_id"]] / point_rates[k][low["tissue_id"]])
                self.assertGreaterEqual(u_ratio * (1 + 1e-12), ratio ** (1 / k))
                self.assertLessEqual(u_ratio, ratio * (1 + 1e-12))
            self.assert_close(contrast["rate_ratios"]["1"], ratio)

    def test_12_outliers_saturation_and_ordering_conditions(self):
        for name, (rows, kind) in self.panels.items():
            panel = self.results[name]
            predictions = panel["descriptive"]["leave_group_out"]["predictions"]
            positives = sorted(predictions, key=lambda p: (-p["residual"], p["tissue_id"]))[:3]
            negatives = sorted(predictions, key=lambda p: (p["residual"], p["tissue_id"]))[:3]
            self.assertEqual(panel["outliers_positive"], positives)
            self.assertEqual(panel["outliers_negative"], negatives)
            maximum = max(r["risk"] for r in rows)
            expected_attenuation = -(1 - maximum) * math.log1p(-maximum) / maximum
            self.assert_close(panel["saturation"]["sigma_at_max_risk"], expected_attenuation)
            count = 0
            for a, b in itertools.combinations(rows, 2):
                delta = a["total_divisions"] - b["total_divisions"]
                if delta * (a["stem_cells"] - b["stem_cells"]) < 0 or delta * (exposure(a, kind) - exposure(b, kind)) < 0:
                    count += 1
            self.assertEqual(panel["saturation"]["ordering_violations"], count)
            self.assertEqual(panel["saturation"]["ordered_cloud_bound_applicable"], count == 0)
            if kind == "effective":
                one = panel["sequential"][0]
                u = one["implied_equal_rate"]
                x = [math.log10(r["total_divisions"]) for r in rows]
                log_r = [math.log10(-math.expm1(-u * r["total_divisions"])) for r in rows]
                self.assert_close(one["predicted_logH_slope"], 1)
                self.assert_close(one["predicted_logR_slope"], pairwise_slope(x, log_r))

    def test_13_later_ranked_risk_audit_cannot_be_a_coherent_panel(self):
        """Optional source audit is inspected separately, never fitted as data."""
        audit_path = RESEARCH / "data" / "liu2022_ranked_risk_audit.csv"
        if not audit_path.exists():
            self.skipTest("The source-coverage audit is documented without a numeric CSV")
        rows = raw_rows(audit_path.name)
        self.assertEqual(len(rows), 17)
        self.assertEqual(len({r["tissue_id"] for r in rows}), 17)
        reported = self.results["international_endpoint_audit"]
        self.assertEqual(reported["n"], len(rows))
        self.assertIn("no coherent risk-plus-divisions fit", reported["status"])
        self.assertEqual(len(reported["rows"]), len(rows))
        for source, result in zip(rows, reported["rows"]):
            self.assertEqual(result["tissue_id"], source["tissue_id"])
            self.assertEqual(result["tissue"], source["tissue"])
            self.assertEqual(result["registry_count"], int(source["registry_count"]))
            self.assertEqual(result["age_endpoint"], source["age_endpoint"])
            for field in ("risk_min", "risk_median", "risk_max"):
                self.assert_close(result[field], float(Decimal(source[field])))
        for row in rows:
            for field in PROVENANCE:
                self.assertTrue(row.get(field, "").strip(), (row, field))
            self.assertNotIn("total_divisions", row)
            self.assertNotIn("stem_cells", row)
            self.assertNotIn("divisions_per_year", row)
            self.assertEqual(int(row["registry_count"]), 423)
            low, median, high = (Decimal(row[k]) for k in ("risk_min", "risk_median", "risk_max"))
            self.assertTrue(0 <= low <= median <= high < 1)
            self.assertIn("rank 1 minimum", row["source_locator"])
            self.assertIn("rank 212 median", row["source_locator"])
            self.assertIn("rank 423 maximum", row["source_locator"])
            self.assertIn("age 74", row["age_endpoint"])
            self.assertIn("not full lifetime", row["age_endpoint"])
            note = " ".join(row.values()).lower()
            self.assertRegex(note, r"independent|unmatched|not matched|not a matched")
        self.assertFalse(any("ranked" in name or "liu2022" in name for name in self.panels))

    def test_14_unordered_target_inverse_and_limiting_amplitude_equivalence(self):
        # Independently reconstruct exact positive-endpoint target inverses in
        # Decimal arithmetic. This is a different process from ordered Erlang
        # waits, despite its equal-amplitude rare-power limit.
        with localcontext() as ctx:
            ctx.prec = 90
            for k, n, v, raw_r in itertools.product(
                    range(1, 9), (1, 3, 100000000), (Decimal("0.2"), Decimal(4)),
                    ("0.0001", "0.03", "0.9")):
                risk = Decimal(raw_r)
                h = -(1 - risk).ln()
                completed = 1 - (-h / n).exp()
                target_hit = (completed.ln() / k).exp()
                intensity = -(1 - target_hit).ln() / v
                forward_q = (1 - (-intensity * v).exp()) ** k
                forward_r = 1 - (1 - forward_q) ** n
                self.assertLess(abs(forward_r / risk - 1), Decimal("1e-70"))
        for k in range(1, 9):
            errors = []
            for x in (0.01, 0.0001, 0.000001):
                scaled = x / math.factorial(k) ** (1 / k)
                unordered_h = -math.log1p(-(-math.expm1(-scaled)) ** k)
                ordered_h = V.psi(k, x)
                errors.append(abs(unordered_h / ordered_h - 1))
            self.assertLess(errors[-1], 2e-6 * k)
            self.assertLessEqual(errors[-1], errors[0] + 1e-14)
        x = 2.0
        target_h = -math.log1p(-(-math.expm1(-x / math.sqrt(2))) ** 2)
        self.assertGreater(abs(target_h - V.psi(2, x)), 0.001)

    def test_15_source_consistent_rectangles_reverse_the_coupled_n_effect(self):
        def source_rate(k, r, n, lifetime, c):
            # All three arguments of the effective-clock inverse change under
            # this source constraint; holding D fixed would test another map.
            d = n * (c + lifetime) - c
            return V.rate(k, r, n, d)

        for k, c in itertools.product(range(1, 9), (0, 1, 2)):
            lower = source_rate(k, 0.0001, 80, 3, c)
            upper = source_rate(k, 0.2, 8, 0.1, c)
            corners = [source_rate(k, r, n, lifetime, c)
                       for r, n, lifetime in itertools.product((0.0001, 0.2), (8, 80), (0.1, 3))]
            self.assert_close(lower, min(corners))
            self.assert_close(upper, max(corners))
            for r, n, lifetime in itertools.product((0.0001, 0.03, 0.2), (8, 25, 80), (0.1, 1, 3)):
                rate = source_rate(k, r, n, lifetime, c)
                self.assertGreaterEqual(rate * (1 + 1e-13), lower)
                self.assertLessEqual(rate, upper * (1 + 1e-13))
            self.assertLess(source_rate(k, 0.03, 80, 1, c), source_rate(k, 0.03, 8, 1, c))
            if k == 1:
                self.assert_close(V.rate(k, 0.03, 80, 1000), V.rate(k, 0.03, 8, 1000))
            else:
                self.assertGreater(V.rate(k, 0.03, 80, 1000), V.rate(k, 0.03, 8, 1000))


if __name__ == "__main__":
    unittest.main()
