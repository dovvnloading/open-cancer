"""Independent offline checks of the published-trend compatibility audit.

Nonclinical, unreviewed computational research. No screening or medical advice.
These standard-library checks establish arithmetic and file consistency, not
source truth, simultaneous statistical coverage, or clinical validity.
"""

from __future__ import annotations

import copy
import csv
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
import math
from pathlib import Path, PurePosixPath, PureWindowsPath
import subprocess
import sys
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
RELATIVE = Path("research") / "breast-age-period-stage"
ROOT = REPO / RELATIVE
CSV_PATH = ROOT / "data" / "lee_argov_table2.csv"
PROTOCOL_COMMIT = "04ad47e180bb6ee1444784f2e0f96e0ef449aafb"
AGES = ("65-74", "75-84", "85+")
STAGES = ("localized", "regional", "distant", "unknown")
MODES = ("primary", "printed_ci", "full_unit", "point_rounding")
NUMERIC_FIELDS = ("printed_aapc", "printed_ci_lower", "printed_ci_upper")
FROZEN_HASHES = {
    (RELATIVE / "PROTOCOL.md").as_posix():
        "d8ad921768f0e387687f17b9a82faf0ba22e0789b9adf833e4584019ddd58d40",
    (RELATIVE / "data" / "lee_argov_table2.csv").as_posix():
        "856f7b1ab0c7bab4a47c005247efb70abc75dacf430380eba0868ce6c58a5f47",
}

# A separate transcription check, not an assertion that hashes verify the source.
SOURCE_TRIPLES = {
    "localized": (("0.95", "0.78", "1.18"),
                  ("0.38", "0.20", "0.74"),
                  ("-0.66", "-0.94", "-0.36")),
    "regional": (("-0.80", "-1.07", "-0.48"),
                 ("-0.92", "-1.23", "-0.60"),
                 ("-0.89", "-1.23", "-0.53")),
    "distant": (("1.13", "0.92", "1.34"),
                ("1.64", "1.33", "1.96"),
                ("1.64", "1.36", "1.95")),
    "unknown": (("-3.99", "-4.61", "-3.34"),
                ("-3.64", "-4.37", "-3.07"),
                ("-3.80", "-4.32", "-3.11")),
}


def load_module():
    spec = importlib.util.spec_from_file_location(
        "research_breast_age_period_stage", ROOT / "analysis.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def raw_cells():
    """Independently read the sole analytical input, retaining printed strings."""
    with CSV_PATH.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    return rows, {(row["stage"], row["age_bin"]): row for row in rows}


def independent_hash(path):
    text = Path(path).read_bytes().decode("utf-8")
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def reference_half_unit(value):
    digits = len(value.split(".")[1])
    return Decimal(5).scaleb(-digits - 1)


def reference_interval(row, mode="primary"):
    """Exact percentage arithmetic, independently implementing the protocol."""
    if mode == "point_rounding":
        left = right = row["printed_aapc"]
        scale = 1
    else:
        left, right = row["printed_ci_lower"], row["printed_ci_upper"]
        scale = {"primary": 1, "printed_ci": 0, "full_unit": 2}[mode]
    return (Decimal(left) - scale * reference_half_unit(left),
            Decimal(right) + scale * reference_half_unit(right))


def reference_log(percent):
    with localcontext() as context:
        context.prec = 70
        return (Decimal(1) + Decimal(percent) / 100).ln()


def reference_summary(cells, stage, ages=AGES, mode="primary"):
    intervals = {age: reference_interval(cells[stage, age], mode) for age in ages}
    lower = max(pair[0] for pair in intervals.values())
    upper = min(pair[1] for pair in intervals.values())
    with localcontext() as context:
        context.prec = 70
        log_lower, log_upper = reference_log(lower), reference_log(upper)
        gap = max(Decimal(0), log_lower - log_upper)
        allowance = gap / 2
        reference = (log_lower + log_upper) / 2
    return {
        "ages": list(ages),
        "verdict": "PASS" if lower <= upper else "FAIL",
        "lower_percent": str(lower),
        "upper_percent": str(upper),
        "lower_log": float(log_lower),
        "upper_log": float(log_upper),
        "lower_ages": [age for age in ages if intervals[age][0] == lower],
        "upper_ages": [age for age in ages if intervals[age][1] == upper],
        "gap": float(gap),
        "allowance": float(allowance),
        "reference": float(reference),
    }


def reference_results(cells):
    cell_results = []
    for stage in STAGES:
        for age in AGES:
            row = cells[stage, age]
            low, high = reference_interval(row)
            point_low, point_high = reference_interval(row, "point_rounding")
            cell_results.append({
                "stage": stage,
                "age_bin": age,
                **{key: row[key] for key in NUMERIC_FIELDS},
                "point_log": float(reference_log(Decimal(row["printed_aapc"]))),
                "point_lower_log": float(reference_log(point_low)),
                "point_upper_log": float(reference_log(point_high)),
                "ci_lower_percent": str(low),
                "ci_upper_percent": str(high),
                "ci_lower_log": float(reference_log(low)),
                "ci_upper_log": float(reference_log(high)),
            })
    modes = {
        mode: {stage: reference_summary(cells, stage, mode=mode) for stage in STAGES}
        for mode in MODES
    }
    primary = modes["primary"]
    failed = [stage for stage in STAGES if primary[stage]["verdict"] == "FAIL"]
    passed = [stage for stage in STAGES if primary[stage]["verdict"] == "PASS"]
    relative_files = (
        Path("PROTOCOL.md"), Path("data") / "lee_argov_table2.csv", Path("analysis.py")
    )
    return {
        "cells": cell_results,
        "modes": modes,
        "leave_one_out": {
            stage: {
                omitted: reference_summary(cells, stage,
                                           tuple(age for age in AGES if age != omitted))
                for omitted in AGES
            } for stage in STAGES
        },
        "overall": {
            "verdict": "FAIL" if failed else "PASS",
            "failed_stages": failed,
            "passing_stages": passed,
            "allowance": max(primary[stage]["allowance"] for stage in STAGES),
        },
        "hashes": {(RELATIVE / file).as_posix(): independent_hash(ROOT / file)
                   for file in relative_files},
        "row_count": len(cells),
    }


class BreastAgePeriodStageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.analysis = load_module()
        cls.rows, cls.cells = raw_cells()
        cls.actual = cls.analysis.analyze()
        cls.expected = reference_results(cls.cells)

    def assert_nested_close(self, actual, expected, location="result"):
        if isinstance(expected, dict):
            self.assertEqual(set(actual), set(expected), location)
            for key, value in expected.items():
                self.assert_nested_close(actual[key], value, f"{location}.{key}")
        elif isinstance(expected, list):
            self.assertEqual(len(actual), len(expected), location)
            for index, value in enumerate(expected):
                self.assert_nested_close(actual[index], value, f"{location}[{index}]")
        elif isinstance(expected, float):
            self.assertTrue(math.isfinite(actual), location)
            self.assertAlmostEqual(actual, expected, delta=2e-15, msg=location)
        else:
            self.assertEqual(actual, expected, location)

    def test_complete_panel_and_all_source_cells(self):
        self.assertEqual(len(self.rows), len(STAGES) * len(AGES))
        self.assertEqual(set(self.cells), {(s, a) for s in STAGES for a in AGES})
        self.assertEqual(tuple(self.analysis.AGES), AGES)
        self.assertEqual(tuple(self.analysis.STAGES), STAGES)
        self.assertEqual(tuple(self.analysis.MODES), MODES)
        for stage in STAGES:
            for age, triple in zip(AGES, SOURCE_TRIPLES[stage]):
                with self.subTest(stage=stage, age=age):
                    row = self.cells[stage, age]
                    self.assertEqual(tuple(row[key] for key in NUMERIC_FIELDS), triple)
                    for value in triple:
                        self.assertRegex(value, r"^-?\d+\.\d{2}$")

    def test_every_row_has_specific_source_and_reuse_provenance(self):
        columns = {"65-74": "65-74 y", "75-84": "75-84 y", "85+": "≥85 y"}
        numeric_or_labels = {"stage", "age_bin", *NUMERIC_FIELDS}
        for row in self.rows:
            with self.subTest(stage=row["stage"], age=row["age_bin"]):
                for key, value in row.items():
                    self.assertTrue(value and value.strip(), key)
                self.assertIn("Lee Argov", row["source"])
                self.assertIn("JAMA Network Open", row["source"])
                self.assertEqual(row["doi"], "10.1001/jamanetworkopen.2025.16947")
                self.assertEqual(row["source_url"],
                                 "https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2835519")
                self.assertEqual(row["table"], "Table 2")
                self.assertEqual(row["source_section"], "Stage at diagnosis")
                expected_row = row["stage"].capitalize()
                if row["stage"] == "unknown":
                    expected_row += " (footnote a)"
                self.assertEqual(row["source_row"], expected_row)
                self.assertEqual(row["source_column"], columns[row["age_bin"]])
                self.assertIn(row["source_row"], row["locator"])
                self.assertIn(row["source_column"], row["locator"])
                self.assertEqual((row["period_start"], row["period_end"]), ("2001", "2019"))
                self.assertEqual(row["measure"], "AAPC")
                self.assertEqual(row["interval_level"], "95%")
                self.assertEqual(row["interval_type"], "published marginal confidence interval")
                self.assertEqual(row["access_date"], "2026-11-18")
                self.assertEqual(row["access_timezone"], "America/New_York")
                self.assertIn("CC-BY", row["reuse_basis"])
                self.assertIn("version not stated", row["reuse_basis"])
                self.assertEqual(row["reuse_url"],
                                 "https://jamanetwork.com/pages/cc-by-license-permissions")
                self.assertNotIn("4.0", row["reuse_basis"])
                self.assertGreater(len(set(row) - numeric_or_labels), 0)
        for absent in ("counts", "annual_counts", "population", "covariance",
                       "midpoint", "age_midpoint", "reclassification"):
            self.assertNotIn(absent, self.rows[0])

    def test_every_reported_numeric_result_against_decimal_reference(self):
        self.assert_nested_close(self.actual, self.expected)

    def test_entire_generated_readme_uses_independently_recomputed_results(self):
        text = (ROOT / "README.md").read_text(encoding="utf-8")
        begin, end = self.analysis.BEGIN, self.analysis.END
        self.assertEqual(text.count(begin), 1)
        self.assertEqual(text.count(end), 1)
        block = text.split(begin, 1)[1].split(end, 1)[0].strip()
        generated = self.analysis.render_results(self.expected).strip()
        # Formatting is shared; every source value and arithmetic field above is not.
        self.assertEqual(block, generated)
        self.assertEqual(generated, self.analysis.render_results(self.actual).strip())
        self.assertIn(PROTOCOL_COMMIT, text)

    def test_primary_verdicts_and_conjunction_retain_passing_stages(self):
        verdicts = {stage: self.actual["modes"]["primary"][stage]["verdict"]
                    for stage in STAGES}
        self.assertEqual(verdicts, {"localized": "FAIL", "regional": "PASS",
                                    "distant": "FAIL", "unknown": "PASS"})
        self.assertEqual(self.actual["overall"]["verdict"], "FAIL")
        self.assertEqual(self.actual["overall"]["passing_stages"], ["regional", "unknown"])
        for stage in self.actual["overall"]["passing_stages"]:
            self.assertEqual(self.actual["modes"]["primary"][stage]["allowance"], 0)

    def test_exact_rounding_boundary_and_all_prespecified_sensitivities(self):
        # Full-unit distant endpoints meet exactly; floating-point equality is not used.
        full = self.actual["modes"]["full_unit"]["distant"]
        self.assertEqual(Decimal(full["lower_percent"]), Decimal(full["upper_percent"]))
        self.assertEqual(full["verdict"], "PASS")
        self.assertEqual(full["allowance"], 0)
        for stage in STAGES:
            for mode in MODES:
                self.assert_nested_close(self.actual["modes"][mode][stage],
                                         reference_summary(self.cells, stage, mode=mode))
            self.assertEqual(set(self.actual["leave_one_out"][stage]), set(AGES))
            for omitted in AGES:
                kept = tuple(age for age in AGES if age != omitted)
                self.assert_nested_close(self.actual["leave_one_out"][stage][omitted],
                                         reference_summary(self.cells, stage, kept))

    def test_sharp_allowance_attainment_and_positive_lower_bound(self):
        # Decimal arithmetic supplies an independent witness and strict lower bound.
        with localcontext() as context:
            context.prec = 70
            for stage in STAGES:
                intervals = [tuple(reference_log(x) for x in reference_interval(self.cells[stage, age]))
                             for age in AGES]
                lower, upper = max(x[0] for x in intervals), min(x[1] for x in intervals)
                allowance = max(Decimal(0), (lower - upper) / 2)
                witness = (lower + upper) / 2
                for low, high in intervals:
                    self.assertLessEqual(low - allowance - Decimal("1e-68"), witness)
                    self.assertLessEqual(witness, high + allowance + Decimal("1e-68"))
                if allowance:
                    smaller = allowance / 2
                    self.assertGreater(lower - smaller, upper + smaller)
                    self.assertAlmostEqual(float(2 * allowance),
                                           self.actual["modes"]["primary"][stage]["gap"],
                                           delta=2e-15)

    def test_rounding_precedes_log_and_point_interval_is_distinct(self):
        for row in self.rows:
            for field in NUMERIC_FIELDS:
                self.assertEqual(self.analysis.half_unit(row[field]), reference_half_unit(row[field]))
            for mode in MODES:
                self.assertEqual(self.analysis.percent_interval(row, mode), reference_interval(row, mode))
            for endpoint in reference_interval(row):
                self.assertAlmostEqual(self.analysis.log_trend(endpoint),
                                       float(reference_log(endpoint)), delta=2e-15)
        self.assertEqual(self.analysis.half_unit("1.2"), Decimal("0.05"))
        self.assertEqual(self.analysis.half_unit("1.200"), Decimal("0.0005"))

    def test_synthetic_closed_intersection_ties_and_permutation(self):
        # Synthetic percentages, not additional observations.
        intervals = {"a": (Decimal("-1"), Decimal("0")),
                     "b": (Decimal("0"), Decimal("1")),
                     "c": (Decimal("0"), Decimal("0"))}
        result = self.analysis.summarize_intervals(intervals)
        self.assertEqual(result["verdict"], "PASS")
        self.assertEqual(result["allowance"], 0)
        self.assertEqual(result["lower_ages"], ["b", "c"])
        self.assertEqual(result["upper_ages"], ["a", "c"])
        self.assertEqual(result["lower_log"], result["upper_log"])
        data = self.analysis.load_data()
        shuffled = dict(reversed(list(data.items())))
        self.assert_nested_close(self.analysis.analyze(shuffled), self.actual)

    def test_synthetic_common_log_translation_cannot_repair_separation(self):
        original = {"a": (Decimal("0"), Decimal("1")),
                    "b": (Decimal("2"), Decimal("3"))}
        # Multiplying each annual factor by 2 is a common translation by log(2).
        shifted = {age: tuple((100 + value) * 2 - 100 for value in interval)
                   for age, interval in original.items()}
        first = self.analysis.summarize_intervals(original)
        second = self.analysis.summarize_intervals(shifted)
        self.assertEqual(first["verdict"], "FAIL")
        self.assertEqual(second["verdict"], first["verdict"])
        self.assertAlmostEqual(second["allowance"], first["allowance"], delta=2e-15)
        self.assertAlmostEqual(second["reference"] - first["reference"],
                               float(Decimal(2).ln()), delta=2e-15)

    def test_synthetic_equal_averages_do_not_establish_full_factorisation(self):
        # Hypothetical exact rational rates: equal endpoint changes, unequal interior ratio.
        one, two = (1, 2, 4), (2, 6, 8)
        self.assertEqual(one[-1] * two[0], two[-1] * one[0])
        self.assertNotEqual(one[1] * two[0], two[1] * one[0])

    def test_log_domain_boundary_is_not_clipped(self):
        for invalid in (Decimal("-100"), Decimal("-100.01"),
                        Decimal("NaN"), Decimal("Infinity")):
            with self.subTest(value=str(invalid)), self.assertRaises((ValueError, ArithmeticError)):
                self.analysis.log_trend(invalid)
        near = Decimal("-99.999")
        self.assertAlmostEqual(self.analysis.log_trend(near), float(reference_log(near)), delta=2e-10)

    def write_fixture(self, directory, rows):
        path = Path(directory) / "fixture.csv"
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(self.rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_malformed_or_missing_inputs_are_errors_without_imputation(self):
        variants = []
        variants.append(("missing cell", copy.deepcopy(self.rows[:-1])))
        variants.append(("duplicate cell", copy.deepcopy(self.rows + [self.rows[0]])))
        mutations = (
            ("stage", "in_situ"), ("age_bin", "all ages"),
            ("printed_aapc", "NA"), ("printed_ci_lower", "NaN"),
            ("printed_ci_upper", "Infinity"), ("printed_ci_lower", "99.00"),
            ("printed_ci_upper", "-99.00"), ("printed_aapc", "5.00"),
            ("source_url", ""), ("reuse_basis", ""),
            ("period_end", "2020"), ("source_column", "75-84 y"),
        )
        for field, value in mutations:
            rows = copy.deepcopy(self.rows)
            rows[0][field] = value
            variants.append((f"{field}={value}", rows))
        invalid_domain = copy.deepcopy(self.rows)
        invalid_domain[0].update(printed_aapc="-99.99", printed_ci_lower="-100.00",
                                 printed_ci_upper="-99.90")
        variants.append(("widened lower endpoint outside log domain", invalid_domain))
        with tempfile.TemporaryDirectory() as temporary:
            for label, rows in variants:
                with self.subTest(case=label):
                    path = self.write_fixture(temporary, rows)
                    with self.assertRaises((ValueError, ArithmeticError)):
                        self.analysis.load_data(path)

    def test_frozen_protocol_and_source_hashes(self):
        for key, expected in FROZEN_HASHES.items():
            with self.subTest(path=key):
                self.assertEqual(independent_hash(REPO / Path(key)), expected)
                self.assertEqual(self.analysis.normalized_sha256(REPO / Path(key)), expected)
        self.assertEqual(self.analysis.file_hashes(), self.expected["hashes"])

    def test_hashes_and_parsing_are_newline_portable(self):
        filenames = (Path("PROTOCOL.md"), Path("data") / "lee_argov_table2.csv", Path("analysis.py"))
        with tempfile.TemporaryDirectory() as temporary:
            for name, newline in (("lf", "\n"), ("crlf", "\r\n"), ("cr", "\r")):
                root = Path(temporary) / name
                for filename in filenames:
                    target = root / filename
                    target.parent.mkdir(parents=True, exist_ok=True)
                    text = (ROOT / filename).read_text(encoding="utf-8")
                    target.write_bytes(text.replace("\n", newline).encode("utf-8"))
                self.assertEqual(self.analysis.file_hashes(root), self.expected["hashes"])
                self.assertEqual(self.analysis.load_data(root / "data" / "lee_argov_table2.csv"),
                                 self.analysis.load_data())

    def test_path_keys_are_posix_on_both_path_flavours(self):
        expected = {key for key in self.expected["hashes"]}
        components = [tuple(Path(key).parts) for key in expected]
        windows = {PureWindowsPath(*parts).as_posix() for parts in components}
        posix = {PurePosixPath(*parts).as_posix() for parts in components}
        native = {Path(*parts).as_posix() for parts in components}
        self.assertEqual(expected, windows)
        self.assertEqual(expected, posix)
        self.assertEqual(expected, native)
        self.assertTrue(all("\\" not in key for key in expected))

    def test_offline_command_line_reproduction(self):
        command = [sys.executable, "-B", str(ROOT / "analysis.py")]
        checked = subprocess.run(command + ["--check"], cwd=REPO, text=True,
                                 capture_output=True, timeout=15, check=False)
        self.assertEqual(checked.returncode, 0, checked.stdout + checked.stderr)
        encoded = subprocess.run(command + ["--format", "json"], cwd=REPO, text=True,
                                 capture_output=True, timeout=15, check=False)
        self.assertEqual(encoded.returncode, 0, encoded.stdout + encoded.stderr)
        self.assert_nested_close(json.loads(encoded.stdout), self.expected)


if __name__ == "__main__":
    unittest.main()
