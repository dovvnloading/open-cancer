"""Independent stdlib/offline guards for an unreviewed nonclinical BMJ audit.

Passing checks cannot confer clinical validity, simultaneous interval coverage
or source truth. Cases/PY MUST NOT be used to reconstruct adjusted IRRs.
"""

from __future__ import annotations

import ast
import copy
import csv
from fractions import Fraction
import hashlib
import importlib.util
import io
import math
from pathlib import Path, PurePosixPath, PureWindowsPath
import socket
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch


REPO = Path(__file__).resolve().parents[1]
REL = Path("research") / "hpv-age-duration-factorisation"
ROOT = REPO / REL
DATA = ROOT / "data" / "wu_figure3.csv"
SPEC = importlib.util.spec_from_file_location("hpv_age_duration_factorisation", ROOT / "analysis.py")
ANALYSIS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ANALYSIS
SPEC.loader.exec_module(ANALYSIS)
F = Fraction
FROZEN_FACTS = {
    ("7-9", "10-16"): ("0.24", "0.10", "0.58", 5, 797583),
    ("7-9", "17+"): ("0.77", "0.52", "1.15", 26, 264588),
    ("10-12", "10-16"): ("0.21", "0.12", "0.40", 11, 668369),
    ("10-12", "17+"): ("0.54", "0.33", "0.86", 18, 195405),
    ("13-15", "10-16"): ("0.23", "0.11", "0.46", 8, 212315),
    ("13-15", "17+"): ("0.23", "0.08", "0.60", 4, 87765),
}


def independently_read_csv():
    with DATA.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def canonical_hash(path):
    text = Path(path).read_bytes().decode("utf-8")
    return hashlib.sha256(text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")).hexdigest()


def reference_boxes(facts, widening):
    """Independent rational implementation directly from the six printed facts."""
    out = {}
    for duration in ("7-9", "10-12", "13-15"):
        e = facts[(duration, "10-16")]
        l = facts[(duration, "17+")]
        if widening == "point":
            e0, e1 = F(e[0]) - F(1, 200), F(e[0]) + F(1, 200)
            l0, l1 = F(l[0]) - F(1, 200), F(l[0]) + F(1, 200)
        else:
            d = {"printed": F(0), "half": F(1, 200),
                 "full": F(1, 100)}[widening]
            e0, e1 = F(e[1])-d, F(e[2])+d
            l0, l1 = F(l[1])-d, F(l[2])+d
        out[duration] = (e0/l1, e1/l0)
    return out


def independent_cdf(k, mean):
    if k < 0:
        return 0.0
    return math.fsum(math.exp(-mean)*mean**j/math.factorial(j) for j in range(k+1))


class TestAgeDurationFactorisation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.perf_counter()
        cls.rows = ANALYSIS.load_data()
        cls.result = ANALYSIS.analyze(cls.rows, include_hashes=False)

    def test_exact_six_transcribed_inputs_and_rights(self):
        rows = independently_read_csv()
        self.assertEqual(len(rows), 6)
        self.assertEqual(set(FROZEN_FACTS), set(self.rows))
        for record in rows:
            key = (record["duration_bin"], record["age_at_first_dose"])
            self.assertEqual((record["adjusted_irr"], record["ci95_lower"],
                              record["ci95_upper"], int(record["vaccinated_cases"]),
                              int(record["vaccinated_person_years"])), FROZEN_FACTS[key])
            self.assertEqual(record["doi"], "10.1136/bmj-2025-087326")
            self.assertEqual(record["license"], "CC BY 4.0")
            self.assertEqual(record["source_url"], ANALYSIS.SOURCE_URL)
            self.assertIn(key[0], record["figure_locator"])
            self.assertEqual(record["interval_kind"], "published marginal 95% CI")

    def test_primary_exact_fractions_independently(self):
        primary = reference_boxes(FROZEN_FACTS, "half")
        self.assertEqual(self.result["boxes"]["primary"], primary)
        expected_lower = F("0.105") / F("0.605")
        expected_upper = F("0.585") / F("0.515")
        summary = self.result["summaries"]["primary"]
        self.assertEqual(summary["lower"], expected_lower)
        self.assertEqual(summary["upper"], expected_upper)
        self.assertEqual(summary["verdict"], "CONSISTENT")
        self.assertEqual(summary["limiting_lower"], ["13-15"])
        self.assertEqual(summary["limiting_upper"], ["7-9"])
        self.assertEqual(summary["allowance"], 0)
        self.assertLess(float(expected_lower), float(expected_upper))
        for lower, upper in primary.values():
            self.assertLessEqual(lower, expected_upper)
            self.assertGreaterEqual(upper, expected_lower)

    def test_rounding_sensitivities_and_exact_point_intervals(self):
        for mode, label in (("no_rounding", "printed"),
                            ("full_unit", "full"), ("point_only", "point")):
            self.assertEqual(self.result["boxes"][mode],
                             reference_boxes(FROZEN_FACTS, label))
        self.assertEqual(self.result["summaries"]["no_rounding"]["verdict"],
                         "CONSISTENT")
        self.assertEqual(self.result["summaries"]["full_unit"]["verdict"],
                         "CONSISTENT")
        point = self.result["summaries"]["point_only"]
        self.assertEqual(point["verdict"], "FALSIFIED")
        self.assertGreater(point["lower"], point["upper"])
        expected = math.log(float(point["lower"] / point["upper"])) / 2
        self.assertAlmostEqual(point["allowance"], expected, delta=1e-12)
        self.assertAlmostEqual(point["allowance"], 0.547566255671, delta=1e-11)
        witness = point["log_witness"]
        # Sharp: witness attains both touching expanded outer endpoints.
        self.assertAlmostEqual(witness,
                               math.log(float(point["lower"])) - point["allowance"],
                               delta=1e-12)
        self.assertAlmostEqual(witness,
                               math.log(float(point["upper"])) + point["allowance"],
                               delta=1e-12)
        self.assertGreater(math.log(float(point["lower"])) - (point["allowance"]-1e-5),
                           math.log(float(point["upper"])) + (point["allowance"]-1e-5))

    def test_all_three_leave_one_bin_out_and_multiplicity_variant(self):
        for omitted in ANALYSIS.BINS:
            rest = {k: v for k, v in self.result["boxes"]["primary"].items()
                    if k != omitted}
            independent = (max(t[0] for t in rest.values()),
                           min(t[1] for t in rest.values()))
            got = self.result["omissions"][omitted]
            self.assertEqual((got["lower"], got["upper"]), independent)
            self.assertEqual(got["verdict"], "CONSISTENT")
            self.assertEqual(got["allowance"], 0)
        approx = self.result["summaries"]["bonferroni_approx"]
        self.assertEqual(approx["verdict"], "CONSISTENT")
        self.assertAlmostEqual(approx["lower"], 0.075008455197, delta=1e-11)
        self.assertAlmostEqual(approx["upper"], 1.996076523107, delta=1e-11)
        self.assertEqual(approx["allowance"], 0)
        self.assertEqual(self.result["inferential_verdict"].split(":")[0], "INCONCLUSIVE")

    def test_counts_are_descriptive_and_poisson_tails_are_exact(self):
        self.assertEqual(len(self.result["counts"]), 6)
        for cell in self.result["counts"]:
            x, years = cell["cases"], cell["years"]
            lo, hi = ANALYSIS.poisson_mean_interval(x)
            interval = cell["crude_rate_interval_per_100k"]
            self.assertAlmostEqual(interval[0], lo*100000/years, delta=1e-10)
            self.assertAlmostEqual(interval[1], hi*100000/years, delta=1e-10)
            self.assertAlmostEqual(1-independent_cdf(x-1, lo), 0.025, delta=2e-12)
            self.assertAlmostEqual(independent_cdf(x, hi), 0.025, delta=2e-12)
        self.assertEqual(ANALYSIS.poisson_mean_interval(0)[0], 0)
        self.assertAlmostEqual(ANALYSIS.poisson_mean_interval(0)[1],
                               -math.log(0.025), delta=1e-13)
        self.assertIn("UNTESTABLE", self.result["zero_bin_adjusted_ratio"])
        self.assertIn("infinity", self.result["zero_bin_adjusted_ratio"])
        # Reweighting descriptive denominators may change only crude diagnostics.
        modified = copy.deepcopy(self.rows)
        modified[("7-9", "10-16")]["vaccinated_person_years"] = "900001"
        modified[("7-9", "10-16")]["vaccinated_cases"] = "4"
        recalculated = ANALYSIS.analyze(modified, include_hashes=False)
        self.assertEqual(self.result["boxes"], recalculated["boxes"])
        self.assertNotEqual(self.result["counts"], recalculated["counts"])

    def test_malformed_data_and_missing_adjusted_irr_fail_closed(self):
        for field, value in (("adjusted_irr", "-0.24"),
                             ("adjusted_irr", "-"),
                             ("ci95_lower", "0.00"),
                             ("vaccinated_cases", "0"),
                             ("license", "CC BY-NC"),
                             ("source_url", "https://example.invalid"),
                             ("figure_locator", "Fig. 3 other bin")):
            rows = copy.deepcopy(self.rows)
            rows[("7-9", "10-16")][field] = value
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "altered.csv"
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=ANALYSIS.COLUMNS)
                    writer.writeheader()
                    writer.writerows(rows.values())
                with self.assertRaisesRegex(ValueError, "UNTESTABLE"):
                    ANALYSIS.load_data(path)
        for value in ("-", "0.0", "0.000", "garbled"):
            with self.assertRaises(ValueError):
                ANALYSIS.decimal_fraction(value)
        with self.assertRaises(ValueError):
            ANALYSIS.poisson_mean_interval(-1)
        with self.assertRaises(ValueError):
            ANALYSIS.poisson_mean_interval(True)
        with self.assertRaises(ValueError):
            ANALYSIS.irr_range(self.rows[("7-9", "10-16")], "unknown")

    def test_normalized_sha_and_path_keys_are_cross_platform(self):
        paths = (Path("PROTOCOL.md"), Path("data") / "wu_figure3.csv", Path("analysis.py"))
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            for item in paths:
                src = ROOT / item
                if not src.exists():
                    # The local container has no GitHub network/clone access;
                    # CI checks the real checked-out frozen protocol.
                    self.assertEqual(item, Path("PROTOCOL.md"))
                    continue
                data = src.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
                dest = folder / item
                dest.parent.mkdir(parents=True, exist_ok=True)
                for ending in (b"\n", b"\r\n", b"\r"):
                    dest.write_bytes(data.replace(b"\n", ending))
                    self.assertEqual(canonical_hash(dest), canonical_hash(src))
                    self.assertEqual(ANALYSIS.normalized_sha256(dest), canonical_hash(src))
            for path in (REL / name for name in paths):
                parts = path.parts
                self.assertEqual(PureWindowsPath(*parts).as_posix(),
                                 PurePosixPath(*parts).as_posix())
                self.assertNotIn("\\", path.as_posix())

    def test_no_network_stdlib_ast_and_count_isolation(self):
        tree = ast.parse((ROOT / "analysis.py").read_text(encoding="utf-8"))
        imports = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name.split(".")[0] for alias in node.names)
            if isinstance(node, ast.ImportFrom) and node.module:
                imports.add(node.module.split(".")[0])
        self.assertTrue(imports <= sys.stdlib_module_names, imports-sys.stdlib_module_names)
        self.assertFalse(imports & {"socket", "urllib", "requests", "http", "subprocess"})
        for func_name in ("irr_range", "ratio_ranges", "bonferroni_irr_range"):
            function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                            and node.name == func_name)
            strings = {n.value for n in ast.walk(function)
                       if isinstance(n, ast.Constant) and isinstance(n.value, str)}
            self.assertNotIn("vaccinated_cases", strings)
            self.assertNotIn("vaccinated_person_years", strings)
        with patch.object(socket, "socket", side_effect=AssertionError("network blocked")), \
             patch.object(socket, "create_connection", side_effect=AssertionError("network blocked")):
            replay = ANALYSIS.analyze(self.rows, include_hashes=False)
        self.assertEqual(replay["boxes"], self.result["boxes"])

    def test_cli_and_report_recomputation(self):
        block = ANALYSIS.format_results(self.result)
        self.assertIn("CONSISTENT", block)
        self.assertIn("FALSIFIED", block)
        self.assertIn("INCONCLUSIVE", block)
        self.assertIn("UNTESTABLE", block)
        self.assertIn("0.547566255671", block)
        self.assertIn("CC BY", (ROOT / "README.md").read_text(encoding="utf-8"))
        source = (ROOT / "REPORT.md").read_text(encoding="utf-8")
        self.assertEqual(source.count(ANALYSIS.BEGIN), 1)
        self.assertEqual(source.count(ANALYSIS.END), 1)
        if (ROOT / "PROTOCOL.md").exists():
            # All actual GitHub/CI checkouts must take the full offline CLI path.
            complete = ANALYSIS.analyze(self.rows)
            embedded = source.split(ANALYSIS.BEGIN)[1].split(ANALYSIS.END)[0].strip()
            self.assertEqual(embedded, ANALYSIS.format_results(complete).strip())
            result = subprocess.run([sys.executable, "-B", str(ROOT/"analysis.py"), "--check"],
                                    cwd=REPO, text=True, capture_output=True, timeout=15)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        else:
            self.assertNotIn("PROTOCOL.md", ANALYSIS.file_hashes.__defaults__ or ())

    def test_zz_under_60_seconds(self):
        self.assertLess(time.perf_counter()-self.started, 60)


if __name__ == "__main__":
    unittest.main()
