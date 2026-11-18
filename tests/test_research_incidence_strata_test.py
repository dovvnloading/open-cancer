"""Offline checks of the frozen published-strata reanalysis; stdlib only.

Nonclinical, unreviewed computational research. No medical advice. These tests
recompute the reported results from the aggregate CSVs, independently evaluate
the exponential shapes, and check content integrity. They cannot establish
source truth, population sampling uncertainty, or biological validity.
"""

from __future__ import annotations

from collections import defaultdict
import copy
import csv
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import product
from math import factorial
from pathlib import Path, PurePosixPath, PureWindowsPath
import sys
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
RESEARCH = REPO / "research" / "incidence-strata-test"
PROTOCOL_COMMIT = "f06822e9f2427cbec8fef16c48589290561b40c3"
PROTOCOL_HASH = "daa9c8b3d838d79043c253cc9c16d679ec3dabfb39938ae16d3ce3cc786ef822"
DATA_NAMES = (
    "female_breast.csv", "colon_and_rectum.csv", "lung_and_bronchus.csv",
    "pancreas.csv",
)
AGES = tuple(f"{start}-{start + 4}" for start in range(30, 85, 5))
RACES = ("All Races", "Whites", "Blacks")
SEXES = ("Both Sexes", "Males", "Females")
KS = tuple(range(2, 11))
GROWTHS = tuple(F(i, 200) for i in range(-12, 13))
HALF = F(1, 20)
BASE_URL = "https://seer.cancer.gov/archive/csr/1975_2017/"
SOURCES = {
    "female_breast": ("4.11", "sect_04_table.11", "4"),
    "colon_and_rectum": ("6.10", "sect_06_table.10", "6"),
    "lung_and_bronchus": ("15.10", "sect_15_table.10", "15"),
    "pancreas": ("22.7", "sect_22_table.07", "22"),
}


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


A = load_module("research_incidence_strata_analysis", RESEARCH / "analysis.py")


def raw_groups(data_dir=RESEARCH / "data"):
    """Separate parser: no production loader or production grouping helper."""
    groups = defaultdict(list)
    for name in DATA_NAMES:
        with (data_dir / name).open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                groups[row["group_id"]].append(row)
    return {
        key: sorted(rows, key=lambda row: int(row["age_start"]))
        for key, rows in groups.items()
    }


def fraction_decimal(value):
    value = F(value)
    return Decimal(value.numerator) / Decimal(value.denominator)


def display_rho(ratio):
    with localcontext() as context:
        context.prec = 160
        return format(fraction_decimal(ratio).sqrt(), ".6f")


def intervals(rows, half=HALF):
    values = [F(row["rate_per_100000"]) for row in rows]
    return [value - half for value in values], [value + half for value in values]


def monotone_reference(rows, half=HALF):
    """Use a prefix maximum for the threshold, then independently list pairs."""
    lower, upper = intervals(rows, half)
    prefix = F(0)
    ratio = F(1)
    for left, right in zip(lower, upper):
        prefix = max(prefix, left)
        ratio = max(ratio, prefix / right)
    violations = {}
    witnesses = set()
    for earlier in range(len(rows)):
        for later in range(earlier + 1, len(rows)):
            key = rows[earlier]["age_group"], rows[later]["age_group"]
            decline = lower[earlier] - upper[later]
            if decline > 0:
                violations[key] = decline
            if ratio > 1 and lower[earlier] / upper[later] == ratio:
                witnesses.add(key)
    return ratio, violations, witnesses


def deletion_reference(rows):
    """Exhaust all retained subsequences using constructive interval feasibility.

    This does not construct the production algorithm's violation graph or its
    minimum vertex covers. Work in exact twentieth-of-a-rate integer units.
    """
    lower, upper = intervals(rows)
    lower = [int(value * 20) for value in lower]
    upper = [int(value * 20) for value in upper]
    ages = [row["age_group"] for row in rows]
    best_size = -1
    answers = set()
    for mask in range(1 << len(rows)):
        size = bin(mask).count("1")
        if size < best_size:
            continue
        prefix = 0
        feasible = True
        for index, (left, right) in enumerate(zip(lower, upper)):
            if mask & (1 << index):
                prefix = max(prefix, left)
                if prefix > right:
                    feasible = False
                    break
        if feasible:
            if size > best_size:
                best_size, answers = size, set()
            answers.add(tuple(age for i, age in enumerate(ages) if not mask & (1 << i)))
    return answers


def power_reference(rows, k, flexible=False):
    """Exact integrals, or endpoint-closure ranges, independent of #64 helpers."""
    lower, upper = intervals(rows)
    shape_lower, shape_upper = [], []
    for row in rows:
        start, stop = F(row["age_start"]), F(row["age_end_exclusive"])
        if flexible:
            shape_lower.append(start ** (k - 1) / factorial(k - 1))
            shape_upper.append(stop ** (k - 1) / factorial(k - 1))
        else:
            value = (stop ** k - start ** k) / ((stop - start) * factorial(k))
            shape_lower.append(value)
            shape_upper.append(value)
    left = [value / shape for value, shape in zip(lower, shape_upper)]
    right = [value / shape for value, shape in zip(upper, shape_lower)]
    L, U = max(left), min(right)
    return (
        max(F(1), L / U),
        {row["age_group"] for row, value in zip(rows, left) if value == L},
        {row["age_group"] for row, value in zip(rows, right) if value == U},
        L, U,
    )


def decimal_bin_shapes():
    """Independent exponential-minus-polynomial integral, not a Taylor loop.

    All signs of growth, every grid index, and every age bin are evaluated at
    160 digits. Exponentials are shared across k; they are not fitted to data.
    """
    shapes = {}
    with localcontext() as context:
        context.prec = 160
        for growth in GROWTHS:
            gd = fraction_decimal(growth)
            exponentials = {
                age: (gd * Decimal(age)).exp() for age in range(30, 90, 5)
            } if growth else {}
            for k in KS:
                values = []
                for start in range(30, 85, 5):
                    stop = start + 5
                    if not growth:
                        value = Decimal(stop ** k - start ** k) / Decimal(5 * factorial(k))
                    else:
                        integrated_exponential = (exponentials[stop] - exponentials[start]) / gd
                        polynomial = sum(
                            (gd ** degree) * Decimal(stop ** (degree + 1) - start ** (degree + 1))
                            / Decimal(factorial(degree + 1))
                            for degree in range(k - 1)
                        )
                        value = (integrated_exponential - polynomial) / (Decimal(5) * gd ** (k - 1))
                    if value <= 0:
                        raise AssertionError("Independent shape must be positive")
                    values.append(value)
                shapes[k, growth] = tuple(values)
    return shapes


def decimal_frontier(rows, shape):
    lower, upper = intervals(rows)
    left = [fraction_decimal(value) / f for value, f in zip(lower, shape)]
    right = [fraction_decimal(value) / f for value, f in zip(upper, shape)]
    L, U = max(left), min(right)
    tolerance = Decimal("1e-110")
    return (
        max(Decimal(1), L / U),
        {row["age_group"] for row, value in zip(rows, left)
         if abs(value - L) <= tolerance * max(abs(L), Decimal(1))},
        {row["age_group"] for row, value in zip(rows, right)
         if abs(value - U) <= tolerance * max(abs(U), Decimal(1))},
    )


class IncidenceStrataTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # The rational shapes and the full 30-curve report are computed once.
        cls.results = A.analyze()
        cls.rows = raw_groups()

    def test_01_complete_frozen_inputs_and_per_row_provenance(self):
        expected_ids = set()
        for site in SOURCES:
            sexes = ("Females",) if site == "female_breast" else SEXES
            for race, sex in product(RACES, sexes):
                expected_ids.add("__".join((site, race.lower().replace(" ", "_"),
                                              sex.lower().replace(" ", "_"))))
        self.assertEqual(set(self.rows), expected_ids)
        self.assertEqual(len(self.rows), 30)
        self.assertEqual(sum(map(len, self.rows.values())), 330)
        self.assertEqual(self.rows, A.load_data())
        self.assertEqual(set(self.results["groups"]), expected_ids)
        for group_id, rows in self.rows.items():
            self.assertEqual([row["age_group"] for row in rows], list(AGES))
            self.assertEqual(len(rows), 11)
            for row in rows:
                table, stem, section = SOURCES[row["cancer_site"]]
                self.assertEqual(row["group_id"], group_id)
                self.assertEqual(row["behavior"], "invasive")
                self.assertEqual(row["registry_set"], "SEER 21 areas")
                self.assertEqual((row["period_start"], row["period_end"]), ("2013", "2017"))
                self.assertEqual(row["rate_units"], "cases per 100000 population per year")
                self.assertEqual(row["reported_decimal_places"], "1")
                self.assertEqual((row["case_count"], row["person_years"]), ("", ""))
                self.assertRegex(row["rate_per_100000"], r"^\d+\.\d$")
                self.assertGreater(F(row["rate_per_100000"]), HALF)
                start, stop = int(row["age_start"]), int(row["age_end_exclusive"])
                self.assertEqual(stop, start + 5)
                self.assertEqual(row["age_group"], f"{start}-{stop - 1}")
                self.assertEqual(row["source_table"], table)
                self.assertEqual(row["source_section"], "Age-Specific SEER Incidence Rates, 2013-2017")
                self.assertEqual(row["source_column"], f"{row['race']} / {row['sex']}")
                self.assertEqual(row["source_page"], "1")
                self.assertEqual(row["source_url"], f"{BASE_URL}download_csr_datafile.php/{stem}.csv")
                self.assertEqual(row["source_pdf_url"], f"{BASE_URL}results_single/{stem}.pdf")
                self.assertEqual(row["source_html_url"],
                                 f"{BASE_URL}browse_csr.php?pageSEL={stem}&sectionSEL={section}")
                self.assertEqual(row["access_date"], "2026-11-18")
                self.assertEqual(row["reuse_terms_url"], BASE_URL)
                self.assertIn("public domain", row["reuse_basis"])
                self.assertIn("National Cancer Institute", row["source_name"])
                self.assertIn("November 2019", row["source_version"])
                self.assertIn("2020-04-15", row["source_version"])
            result = self.results["groups"][group_id]
            self.assertEqual(result["status"], "testable")
            self.assertEqual(result["n"], len(rows))
            self.assertEqual(result["site"], rows[0]["cancer_site"])
            self.assertEqual(result["race"], rows[0]["race"])
            self.assertEqual(result["sex"], rows[0]["sex"])

    def test_02_frozen_protocol_hashes_and_windows_path_keys(self):
        self.assertEqual(self.results["protocol_commit"], PROTOCOL_COMMIT)
        expected_paths = {Path("PROTOCOL.md").as_posix()}
        expected_paths.update((Path("data") / name).as_posix() for name in DATA_NAMES)
        self.assertEqual(set(self.results["hashes"]), expected_paths)
        self.assertEqual(self.results["hashes"][Path("PROTOCOL.md").as_posix()], PROTOCOL_HASH)
        for key, digest in self.results["hashes"].items():
            self.assertNotIn("\\", key)
            path = RESEARCH / key
            normalized = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
            self.assertEqual(digest, hashlib.sha256(normalized).hexdigest())
            self.assertEqual(digest, A.source_sha256(path))
        inherited = self.results["inherited_method"]
        inherited_path = Path("research") / "multistage-incidence-identifiability" / "verifier.py"
        self.assertEqual(inherited["path"], inherited_path.as_posix())
        normalized = (REPO / inherited_path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
        self.assertEqual(inherited["sha256"], hashlib.sha256(normalized).hexdigest())
        for root in (PureWindowsPath(r"C:\work\incidence-strata-test"),
                     PurePosixPath("/work/incidence-strata-test")):
            for name in DATA_NAMES:
                canonical = (Path("data") / name).as_posix()
                self.assertEqual(A.path_key(root / "data" / name, root), canonical)
                self.assertIn(canonical, self.results["hashes"])
            self.assertEqual(A.path_key(root / "PROTOCOL.md", root), Path("PROTOCOL.md").as_posix())
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "newlines.csv"
            original = "column,rate\né,29.4\nlast".encode("utf-8")
            digest = hashlib.sha256(original).hexdigest()
            for ending in (b"\n", b"\r\n", b"\r"):
                path.write_bytes(original.replace(b"\n", ending))
                self.assertEqual(A.source_sha256(path), digest)
            path.write_bytes(original + b"\n")
            self.assertNotEqual(A.source_sha256(path), digest)
            path.write_bytes(original.replace(b"29.4", b"29.5"))
            self.assertNotEqual(A.source_sha256(path), digest)

    def assert_monotonicity(self, rows, actual, half=HALF):
        ratio, declines, witnesses = monotone_reference(rows, half)
        self.assertEqual(actual["rho_squared"], ratio)
        self.assertEqual(actual["rho"], display_rho(ratio))
        self.assertEqual(actual["compatible"], ratio == 1)
        actual_declines = {
            (value["earlier_age"], value["later_age"]): value["minimum_decline"]
            for value in actual["violations"]
        }
        self.assertEqual(len(actual_declines), len(actual["violations"]))
        self.assertEqual(actual_declines, declines)
        by_age = {row["age_group"]: row["rate_per_100000"] for row in rows}
        for pair in actual["violations"]:
            earlier, later = by_age[pair["earlier_age"]], by_age[pair["later_age"]]
            self.assertEqual(pair["earlier_rate"], earlier)
            self.assertEqual(pair["later_rate"], later)
            self.assertEqual(pair["rho_squared"], (F(earlier) - half) / (F(later) + half))
        self.assertEqual({(value["earlier_age"], value["later_age"])
                          for value in actual["limiting_pairs"]}, witnesses)
        self.assertEqual(len(actual["limiting_pairs"]), len(witnesses))
        # Multiply the feasible isotonic envelope by rho to work in exact
        # squared coordinates. No rounded square root decides feasibility.
        lower, upper = intervals(rows, half)
        prefix = F(0)
        for left, right in zip(lower, upper):
            prefix = max(prefix, left)
            self.assertLessEqual(prefix, ratio * right)
        if ratio > 1:
            below = (ratio + 1) / 2
            self.assertTrue(any(max(lower[:i + 1]) > below * right
                                for i, right in enumerate(upper)))

    def test_03_every_rounding_and_point_monotonicity_result(self):
        for group_id, rows in self.rows.items():
            with self.subTest(group=group_id):
                result = self.results["groups"][group_id]
                self.assert_monotonicity(rows, result["monotonicity"])
                self.assert_monotonicity(rows, result["point_monotonicity"], F(0))
                self.assertLessEqual(result["monotonicity"]["rho_squared"],
                                     result["point_monotonicity"]["rho_squared"])

    def test_04_every_age_prefix_deletion_and_minimum_deletion_set(self):
        for group_id, rows in self.rows.items():
            with self.subTest(group=group_id):
                actual = self.results["groups"][group_id]["age_sensitivity"]
                prefixes = actual["prefixes"]
                self.assertEqual([entry["ending_age"] for entry in prefixes], list(AGES))
                first_failure = None
                previous = F(1)
                for count, entry in enumerate(prefixes, start=1):
                    ratio = monotone_reference(rows[:count])[0]
                    self.assertEqual(entry["rho_squared"], ratio)
                    self.assertEqual(entry["rho"], display_rho(ratio))
                    self.assertGreaterEqual(ratio, previous)
                    previous = ratio
                    if first_failure is None and ratio > 1:
                        first_failure = entry["ending_age"]
                self.assertEqual(actual["first_failure_age"], first_failure)
                deletions = actual["leave_one_out"]
                self.assertEqual([entry["removed_age"] for entry in deletions], list(AGES))
                for index, entry in enumerate(deletions):
                    retained = rows[:index] + rows[index + 1:]
                    ratio = monotone_reference(retained)[0]
                    self.assertEqual(entry["rho_squared"], ratio)
                    self.assertEqual(entry["rho"], display_rho(ratio))
                    self.assertLessEqual(ratio, prefixes[-1]["rho_squared"])
                expected = deletion_reference(rows)
                self.assertEqual({tuple(value) for value in actual["minimum_deletion_sets"]}, expected)
                self.assertEqual(len(actual["minimum_deletion_sets"]), len(expected))
                for deleted in expected:
                    kept = [row for row in rows if row["age_group"] not in deleted]
                    self.assertEqual(monotone_reference(kept)[0], 1)

    def test_05_every_exact_power_and_unknown_weight_frontier(self):
        for group_id, rows in self.rows.items():
            result = self.results["groups"][group_id]
            self.assertEqual([entry["k"] for entry in result["by_k"]], list(KS))
            pure_ratios = {}
            for entry in result["by_k"]:
                k = entry["k"]
                for field, flexible in (("power", False), ("flexible_power", True)):
                    with self.subTest(group=group_id, k=k, family=field):
                        ratio, lower_ages, upper_ages, L, U = power_reference(rows, k, flexible)
                        actual = entry[field]
                        self.assertEqual(actual["rho_squared"], ratio)
                        self.assertEqual(actual["rho"], display_rho(ratio))
                        self.assertEqual(set(actual["lower_ages"]), lower_ages)
                        self.assertEqual(set(actual["upper_ages"]), upper_ages)
                        self.assertLessEqual(L, ratio * U)
                        if ratio > 1:
                            self.assertGreater(L, ((ratio + 1) / 2) * U)
                        self.assertGreaterEqual(ratio, result["monotonicity"]["rho_squared"])
                pure_ratios[k] = entry["power"]["rho_squared"]
                self.assertLessEqual(entry["flexible_power"]["rho_squared"],
                                     entry["power"]["rho_squared"])
            best = min(pure_ratios.values())
            expected = {k for k, value in pure_ratios.items() if value == best}
            self.assertEqual(set(result["best_power_contenders"]), expected)
            self.assertIn(result["best_power_k"], expected)

    def test_06_all_one_growth_candidates_against_independent_closed_forms(self):
        shapes = decimal_bin_shapes()
        with localcontext() as context:
            context.prec = 160
            tolerance = Decimal("1e-110")
            for group_id, rows in self.rows.items():
                result = self.results["groups"][group_id]
                global_values = {}
                for entry in result["by_k"]:
                    k, actual = entry["k"], entry["grid"]
                    independent = {
                        g: decimal_frontier(rows, shapes[k, g]) for g in GROWTHS
                    }
                    optimum = min(value[0] for value in independent.values())
                    winners = {g for g, value in independent.items()
                               if abs(value[0] - optimum) <= tolerance * max(optimum, Decimal(1))}
                    with self.subTest(group=group_id, k=k):
                        self.assertEqual(set(actual["contenders"]), winners)
                        self.assertIn(actual["selected_g"], winners)
                        self.assertEqual(actual["unique"], len(winners) == 1)
                        lo, hi = actual["rho_squared_interval"]
                        self.assertLessEqual(lo, hi)
                        self.assertLess(hi - lo, F(1, 10 ** 24))
                        if lo == hi:
                            # At g=0 the enclosure is a single exact rational.
                            # Successive Decimal divisions can differ by a last
                            # digit from dividing that rational once. Check the
                            # exact reference, and bound presentation roundoff;
                            # do not widen the production numerical certificate.
                            if actual["selected_g"] == 0:
                                self.assertEqual(lo, power_reference(rows, k)[0])
                            else:
                                self.assertEqual(lo, 1)
                            self.assertLessEqual(abs(fraction_decimal(lo) - optimum),
                                                 Decimal("1e-140") * max(optimum, Decimal(1)))
                        else:
                            self.assertLessEqual(fraction_decimal(lo), optimum)
                            self.assertGreaterEqual(fraction_decimal(hi), optimum)
                        expected_rho = format(optimum.sqrt(), ".6f")
                        self.assertEqual(actual["rho"], expected_rho)
                        self.assertEqual(display_rho(lo), expected_rho)
                        self.assertEqual(display_rho(hi), expected_rho)
                        _, lower_ages, upper_ages = independent[actual["selected_g"]]
                        self.assertEqual(set(actual["lower_ages"]), lower_ages)
                        self.assertEqual(set(actual["upper_ages"]), upper_ages)
                        self.assertLessEqual(lo, entry["power"]["rho_squared"])
                        self.assertGreaterEqual(hi, result["monotonicity"]["rho_squared"])
                    for growth, reference in independent.items():
                        global_values[k, growth] = reference[0]
                best = min(global_values.values())
                winners = {key for key, value in global_values.items()
                           if abs(value - best) <= tolerance * max(best, Decimal(1))}
                actual_winners = {(entry["k"], entry["g"])
                                  for entry in result["best_grid_contenders"]}
                self.assertEqual(actual_winners, winners)
                self.assertIn(result["best_grid_k"], {k for k, _ in winners})

    def test_07_primary_verdict_and_prior_44_rate_benchmark(self):
        breast = {
            race: monotone_reference(self.rows[f"female_breast__{race}__females"])[0]
            for race in ("all_races", "whites", "blacks")
        }
        if breast["all_races"] <= 1:
            expected = "untestable"
        elif breast["whites"] > 1 or breast["blacks"] > 1:
            expected = "shape failure persists within strata"
        else:
            expected = "shape failure disappears within the prespecified strata"
        self.assertEqual(self.results["verdict"], expected)
        self.assertEqual(display_rho(breast["all_races"]), "1.066238")
        original_path = REPO / "research" / "multistage-incidence-identifiability" / "data" / "seer_age_specific_2013_2017.csv"
        with original_path.open(encoding="utf-8", newline="") as stream:
            original_rows = list(csv.DictReader(stream))
        self.assertEqual(len(original_rows), 44)
        for row in original_rows:
            sex = "females" if row["cancer_site"] == "female_breast" else "both_sexes"
            group_id = f"{row['cancer_site']}__all_races__{sex}"
            current = {entry["age_group"]: entry["rate_per_100000"]
                       for entry in self.rows[group_id]}
            self.assertEqual(current[row["age_group"]], row["rate_per_100000"])

    def test_08_rounding_equality_nonadjacent_violations_and_ties(self):
        template = self.rows["female_breast__all_races__females"]
        for printed in (("1.1", "1.0"), ("1.2", "1.1", "1.0"),
                        ("2.0", "2.0", "1.0", "1.0"), ("1.0",)):
            rows = copy.deepcopy(template[:len(printed)])
            for row, rate in zip(rows, printed):
                row["rate_per_100000"] = rate
            self.assert_monotonicity(rows, A.monotonicity(rows))
        equal = copy.deepcopy(template[:2])
        equal[0]["rate_per_100000"], equal[1]["rate_per_100000"] = "1.1", "1.0"
        self.assertEqual(A.monotonicity(equal)["rho_squared"], 1)
        nonadjacent = copy.deepcopy(template[:3])
        for row, rate in zip(nonadjacent, ("1.2", "1.1", "1.0")):
            row["rate_per_100000"] = rate
        self.assertEqual(len(A.monotonicity(nonadjacent)["violations"]), 1)
        self.assertGreater(A.monotonicity(nonadjacent)["rho_squared"], 1)
        with self.assertRaises(ValueError):
            A.monotonicity([])
        with self.assertRaises(ValueError):
            A.monotonicity(template, half=-HALF)
        zero_lower = copy.deepcopy(template[:1])
        zero_lower[0]["rate_per_100000"] = "0.05"
        with self.assertRaises(ValueError):
            A.monotonicity(zero_lower)

    def test_09_missing_suppressed_and_nonpositive_lower_cells_abstain(self):
        original = self.rows["female_breast__all_races__females"]
        self.assertIsNone(A.eligibility_reason(original))
        for label in ("missing_age", "suppressed", "blank", "nonnumeric", "zero", "lower_zero"):
            rows = copy.deepcopy(original)
            if label == "missing_age":
                rows.pop(4)
            else:
                rows[4]["rate_per_100000"] = {
                    "suppressed": "-", "blank": "", "nonnumeric": "unavailable",
                    "zero": "0.0", "lower_zero": "0.05",
                }[label]
            with self.subTest(case=label):
                reason = A.eligibility_reason(rows)
                self.assertIsInstance(reason, str)
                self.assertTrue(reason)

    def test_10_duplicate_and_inconsistent_metadata_fail_loading(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "data"
            target.mkdir()
            originals = {}
            for name in DATA_NAMES:
                with (RESEARCH / "data" / name).open(encoding="utf-8", newline="") as stream:
                    originals[name] = list(csv.DictReader(stream))
            for label in ("duplicate", "period", "column", "page", "denominator", "age_label"):
                for name, original in originals.items():
                    rows = copy.deepcopy(original)
                    if name == "female_breast.csv":
                        if label == "duplicate":
                            rows.append(dict(rows[0]))
                        else:
                            field, value = {
                                "period": ("period_end", "2018"),
                                "column": ("source_column", "All Races / Males"),
                                "page": ("source_page", ""),
                                "denominator": ("person_years", "100000"),
                                "age_label": ("age_group", "31-35"),
                            }[label]
                            rows[0][field] = value
                    with (target / name).open("w", encoding="utf-8", newline="") as stream:
                        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                        writer.writeheader()
                        writer.writerows(rows)
                with self.subTest(case=label), self.assertRaises(ValueError):
                    A.load_data(target)

    def test_11_entire_published_result_block_is_regenerated(self):
        text = (RESEARCH / "README.md").read_text(encoding="utf-8")
        begin, end = "<!-- BEGIN RECOMPUTED RESULTS -->", "<!-- END RECOMPUTED RESULTS -->"
        self.assertEqual(text.count(begin), 1)
        self.assertEqual(text.count(end), 1)
        actual = text.split(begin, 1)[1].split(end, 1)[0].strip()
        expected = A.render_results(self.results).strip()
        if begin in expected:
            expected = expected.split(begin, 1)[1].split(end, 1)[0].strip()
        self.assertEqual(actual, expected)
        self.assertIn(PROTOCOL_COMMIT, text)
        self.assertIn("nonclinical", text.lower())
        self.assertIn("no medical advice", text.lower())

    def test_12_verdict_branches_preserve_partial_and_unavailable_outcomes(self):
        races = ("all_races", "whites", "blacks")
        for passing in product((False, True), repeat=3):
            for unavailable in product((False, True), repeat=3):
                groups = {
                    f"female_breast__{race}__females": {
                        "status": "untestable" if missing else "testable",
                        "monotonicity": {"compatible": compatible},
                    }
                    for race, compatible, missing in zip(races, passing, unavailable)
                }
                # Prespecified contrast outcomes cannot change the primary rule.
                groups["pancreas__all_races__both_sexes"] = {
                    "status": "testable", "monotonicity": {"compatible": False},
                }
                if any(unavailable) or passing[0]:
                    expected = "untestable"
                elif passing[1] and passing[2]:
                    expected = "shape failure disappears within the prespecified strata"
                else:
                    expected = "shape failure persists within strata"
                with self.subTest(passing=passing, unavailable=unavailable):
                    self.assertEqual(A.primary_verdict(groups), expected)

    def test_13_absent_panel_returns_abstentions_without_fitting(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "data"
            target.mkdir()
            result = A.analyze(target)
        self.assertEqual(result["verdict"], "untestable")
        self.assertEqual(set(result["groups"]), set(self.rows))
        for group in result["groups"].values():
            self.assertEqual(group["status"], "untestable")
            self.assertEqual(group["n"], 0)
            self.assertTrue(group["reason"])
            self.assertNotIn("by_k", group)
            self.assertNotIn("monotonicity", group)


if __name__ == "__main__":
    unittest.main()
