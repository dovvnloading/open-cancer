"""Independent, offline checks of the frozen aggregate breast-incidence audit.

Nonclinical, unreviewed computational research. No screening or medical advice.
These checks recompute the displayed analysis from factual aggregate CSV cells.
They establish neither source truth, a causal effect, nor clinical validity.
Only Python's standard library is used.
"""

from __future__ import annotations

import copy
import csv
from fractions import Fraction as F
import hashlib
import importlib.util
import math
from pathlib import Path, PurePosixPath, PureWindowsPath
import sys
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
RESEARCH = REPO / "research" / "breast-incidence-decline-stage"
CSV_PATH = RESEARCH / "data" / "turbow_table1.csv"
PROTOCOL_COMMIT = "523133aeb2949e82aa5a0411a61ce5cce9f486a8"
RELATIVE = Path("research") / "breast-incidence-decline-stage"
FROZEN_HASHES = {
    (RELATIVE / "PROTOCOL.md").as_posix():
        "ff029700849020859911e9f8a1ace73eaa9acd6939e1cba912489a1f96930b14",
    (RELATIVE / "data" / "turbow_table1.csv").as_posix():
        "49d90dd3ea6eb9c7cffb5bf24bc0ea45491676469d9624676fe250d03797eeb0",
}
AGES = ("65-69", "70-74", "75-79", "80-84", "85+")
STAGES = ("overall", "localized", "regional", "distant", "unknown")
NAMED = ("localized", "regional", "distant")
EDGES = tuple(zip(AGES[:3], AGES[1:4]))
ALPHA = 0.05
FAMILY = 15


def load_module():
    spec = importlib.util.spec_from_file_location(
        "research_breast_incidence_decline_stage", RESEARCH / "analysis.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


A = load_module()


def raw_cells():
    """Independent CSV parser; never calls the production data loader."""
    with CSV_PATH.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    return rows, {(row["age_bin"], row["stage"]): row for row in rows}


def crude(row):
    return F(100000 * int(row["cases"]), int(row["population"]))


def count(row):
    return int(row["cases"])


def exposure(row):
    return int(row["population"])


def residual_reference(cells):
    return {
        age: count(cells[age, "overall"])
        - sum(count(cells[age, stage]) for stage in STAGES[1:])
        for age in AGES
    }


def _lower_tail_below_mean(n, k, p):
    """Direct binomial mass recurrence, independent of incomplete beta.

    The caller chooses the lower tail below the mean, so successive terms
    decrease. The omitted positive terms are bounded by a geometric series.
    math.lgamma is used only for the first probability mass; the remaining
    masses are related by exact distributional probability ratios.
    """
    log_mass = (math.lgamma(n + 1) - math.lgamma(k + 1)
                - math.lgamma(n - k + 1)
                + k * math.log(p) + (n - k) * math.log1p(-p))
    mass = math.exp(log_mass)
    terms = [mass]
    j = k
    while j > 0:
        ratio = j * (1 - p) / ((n - j + 1) * p)
        mass *= ratio
        terms.append(mass)
        j -= 1
        if mass == 0:
            break
        next_ratio = j * (1 - p) / ((n - j + 1) * p)
        if next_ratio < 1:
            remainder = mass * next_ratio / (1 - next_ratio)
            if remainder <= 1e-15 * math.fsum(terms):
                break
    return math.fsum(terms)


def binomial_leq(n, k, p):
    """P[Binomial(n,p) <= k] by a direct, numerically small tail sum."""
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    if p == 0:
        return 1.0
    if p == 1:
        return 0.0
    if k <= n * p:
        return _lower_tail_below_mean(n, k, p)
    return 1 - _lower_tail_below_mean(n, n - k - 1, 1 - p)


def binomial_geq(n, k, p):
    return binomial_leq(n, n - k, 1 - p)


def probability_from_rr(rr, younger_exposure, older_exposure):
    if math.isinf(rr):
        return 1.0
    return rr * older_exposure / (younger_exposure + rr * older_exposure)


def newline_hash(path):
    text = path.read_bytes().decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


class BreastIncidenceStageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.raw = raw_cells()
        cls.data = A.load_data()
        cls.report = A.analyze()
        cls.residual = residual_reference(cls.raw)

    def assert_close(self, actual, expected, rel=3e-8, absolute=3e-11):
        self.assertTrue(
            math.isclose(float(actual), float(expected), rel_tol=rel, abs_tol=absolute),
            f"{actual!r} differs from independent reference {expected!r}",
        )

    def assert_interval_tails(self, result, younger, older, ey, eo):
        """Check mathematical defining equations, not production inversion."""
        n = younger + older
        tail = ALPHA / (2 * FAMILY)
        if older:
            p = probability_from_rr(result["lower"], ey, eo)
            self.assert_close(binomial_geq(n, older, p), tail)
        else:
            self.assertEqual(result["lower"], 0)
        if younger:
            p = probability_from_rr(result["upper"], ey, eo)
            self.assert_close(binomial_leq(n, older, p), tail)
        else:
            self.assertTrue(math.isinf(result["upper"]))
        if "one_sided_upper" in result and younger:
            p = probability_from_rr(result["one_sided_upper"], ey, eo)
            self.assert_close(binomial_leq(n, older, p), ALPHA / FAMILY)
            self.assertLessEqual(result["one_sided_upper"], result["upper"])

    def test_01_panel_provenance_and_frozen_protocol(self):
        self.assertEqual(len(self.rows), 25)
        expected_keys = {(age, stage) for age in AGES for stage in STAGES}
        self.assertEqual(set(self.raw), expected_keys)
        self.assertEqual(set(self.data), expected_keys)
        self.assertEqual(tuple(A.AGES), AGES)
        self.assertEqual(tuple(A.STAGES), STAGES)
        self.assertEqual(tuple(A.CONTRASTS), EDGES)
        self.assertEqual(A.FAMILY_SIZE, FAMILY)
        self.assertEqual(A.ALPHA, ALPHA)
        for key, row in self.raw.items():
            age, stage = key
            with self.subTest(cell=key):
                self.assertTrue(all(row.values()))
                self.assertEqual(row["doi"], "10.1007/s10549-021-06160-4")
                self.assertEqual(row["source_url"],
                                 "https://pmc.ncbi.nlm.nih.gov/articles/PMC10846538/#T1")
                self.assertEqual(row["table"], "Table 1")
                self.assertEqual(row["source_row"],
                                 "Overall" if stage == "overall" else "Stage: " + stage.title())
                self.assertEqual(row["source_column"], ">=85" if age == "85+" else age)
                self.assertEqual(row["population_locator"], "Table 1 age-column header n")
                self.assertEqual(row["access_date"], "2026-11-18")
                self.assertEqual(row["access_timezone"], "America/New_York")
                self.assertIn("19 age groups", row["rate_adjustment"])
                self.assertIn("2016 USCS", row["dataset_version"])
                self.assertIn("factual aggregate statistics", row["reuse_basis"])
                self.assertEqual(row["reuse_url"], "https://www.copyright.gov/help/faq/faq-protect.html")
                self.assertEqual(self.data[key]["cases"], count(row))
                self.assertEqual(self.data[key]["population"], exposure(row))
                for field in ("printed_rate", "printed_ci_lower", "printed_ci_upper"):
                    self.assertEqual(self.data[key][field], row[field])
                self.assertGreater(exposure(row), 0)
                self.assertGreaterEqual(count(row), 0)
                self.assertLess(F(row["printed_ci_lower"]), F(row["printed_rate"]))
                self.assertLess(F(row["printed_rate"]), F(row["printed_ci_upper"]))
        for age in AGES:
            self.assertEqual(len({exposure(self.raw[age, s]) for s in STAGES}), 1)
        for relative, expected in FROZEN_HASHES.items():
            self.assertEqual(newline_hash(REPO / relative), expected)

    def test_02_every_reconciliation_cell(self):
        report = self.report["reconciliation"]
        self.assertEqual(len(report), 25)
        self.assertEqual({(r["age"], r["stage"]) for r in report}, set(self.raw))
        for row in report:
            key = row["age"], row["stage"]
            raw = self.raw[key]
            observed = crude(raw)
            printed = F(raw["printed_rate"])
            with self.subTest(cell=key):
                self.assertEqual(row["cases"], count(raw))
                self.assertEqual(row["population"], exposure(raw))
                self.assertEqual(A.rate(self.data[key]), observed)
                self.assertEqual(row["computed_rate"], float(observed))
                self.assertEqual(row["printed_rate"], float(printed))
                self.assertEqual(row["printed_ci_lower"], float(raw["printed_ci_lower"]))
                self.assertEqual(row["printed_ci_upper"], float(raw["printed_ci_upper"]))
                self.assertEqual(row["rate_difference"], float(observed - printed))
                self.assertEqual(row["rounding_matches"], abs(observed - printed) <= F(1, 20))
                self.assertEqual(row["printed_ci_contains_rate"],
                                 F(raw["printed_ci_lower"]) <= observed <= F(raw["printed_ci_upper"]))

    def test_03_accounting_retains_unassigned_cases(self):
        self.assertEqual(A.residuals(self.data), self.residual)
        self.assertEqual(self.residual, {"65-69": 2, "70-74": 0, "75-79": 0, "80-84": 1, "85+": 1})
        self.assertEqual(len(self.report["accounting"]), 5)
        for row in self.report["accounting"]:
            age = row["age"]
            overall = count(self.raw[age, "overall"])
            stage_sum = sum(count(self.raw[age, stage]) for stage in STAGES[1:])
            self.assertEqual(row["overall"], overall)
            self.assertEqual(row["stage_sum"], stage_sum)
            self.assertEqual(row["residual"], overall - stage_sum)
            self.assertEqual(row["unknown_pool"], count(self.raw[age, "unknown"]) + self.residual[age])

    def test_04_independent_binomial_reference_on_small_exact_examples(self):
        for n in (1, 2, 5, 12, 30):
            for rational in (F(1, 10), F(1, 2), F(4, 5)):
                p = float(rational)
                for k in range(-1, n + 1):
                    exact = sum((F(math.comb(n, j)) * rational ** j * (1 - rational) ** (n - j)
                                 for j in range(max(0, k + 1))), F(0))
                    self.assert_close(binomial_leq(n, k, p), exact, rel=1e-12, absolute=1e-13)

    def test_05_all_primary_exact_intervals_p_values_and_decisions(self):
        rows = self.report["contrasts"]
        self.assertEqual(len(rows), FAMILY)
        self.assertEqual({(r["younger"], r["older"], r["stage"]) for r in rows},
                         {(young, old, stage) for young, old in EDGES for stage in STAGES})
        for row in rows:
            young = self.raw[row["younger"], row["stage"]]
            old = self.raw[row["older"], row["stage"]]
            y, o, ey, eo = count(young), count(old), exposure(young), exposure(old)
            p0 = eo / (ey + eo)
            p_decline = binomial_leq(y + o, o, p0)
            p_increase = binomial_geq(y + o, o, p0)
            p_two = min(1, 2 * min(p_decline, p_increase))
            with self.subTest(age=(row["younger"], row["older"]), stage=row["stage"]):
                self.assert_close(row["rr"], F(o * ey, y * eo), rel=1e-14)
                self.assertEqual(row["rate_younger"], float(crude(young)))
                self.assertEqual(row["rate_older"], float(crude(old)))
                self.assertEqual(row["change"], float(crude(old) - crude(young)))
                self.assert_interval_tails(row, y, o, ey, eo)
                self.assert_close(row["p_one_sided"], p_decline, absolute=1e-300)
                self.assert_close(row["p_two_sided"], p_two, absolute=1e-300)
                self.assert_close(row["adjusted_two_sided_p"], min(1, FAMILY * p_two), absolute=1e-300)
                decision = ("decline" if row["upper"] < 1
                            else "increase" if row["lower"] > 1 else "unresolved")
                self.assertEqual(row["decision"], decision)
                self.assertEqual(row["one_sided_decline"], row["one_sided_upper"] < 1)
                self.assertEqual(row["decision"] != "unresolved", p_two < ALPHA / FAMILY)

    def test_06_every_signed_decomposition_and_localized_share_bound(self):
        self.assertEqual(len(self.report["decomposition"]), len(EDGES))
        overall_decisions = {(r["younger"], r["older"]): r["decision"]
                             for r in self.report["contrasts"] if r["stage"] == "overall"}
        for row in self.report["decomposition"]:
            young, old = row["younger"], row["older"]
            self.assertIn((young, old), EDGES)
            losses = {stage: crude(self.raw[young, stage]) - crude(self.raw[old, stage])
                      for stage in STAGES}
            ey, eo = exposure(self.raw[young, "overall"]), exposure(self.raw[old, "overall"])
            losses["residual"] = F(100000 * self.residual[young], ey) - F(100000 * self.residual[old], eo)
            self.assertEqual(row["losses"], {key: float(value) for key, value in losses.items()})
            self.assertEqual(sum(losses[s] for s in STAGES[1:]) + losses["residual"], losses["overall"])
            if losses["overall"] > 0:
                shares = {s: losses[s] / losses["overall"] for s in (*STAGES[1:], "residual")}
                self.assertEqual(row["shares"], {key: float(value) for key, value in shares.items()})
                self.assertEqual(sum(shares.values()), 1)
                pool_y = count(self.raw[young, "unknown"]) + self.residual[young]
                pool_o = count(self.raw[old, "unknown"]) + self.residual[old]
                minimum = (losses["localized"] - F(100000 * pool_o, eo)) / losses["overall"]
                maximum = (losses["localized"] + F(100000 * pool_y, ey)) / losses["overall"]
                self.assertEqual(row["localized_share_min"], float(minimum))
                self.assertEqual(row["localized_share_max"], float(maximum))
                robust_overall = overall_decisions[young, old] == "decline"
                self.assertEqual(row["recorded_majority"], robust_overall and shares["localized"] > F(1, 2))
                self.assertEqual(row["allocation_majority"], robust_overall and minimum > F(1, 2))
            else:
                self.assertIsNone(row["shares"])
                self.assertIsNone(row["localized_share_min"])
                self.assertIsNone(row["localized_share_max"])
                self.assertFalse(row["recorded_majority"])
                self.assertFalse(row["allocation_majority"])

    def test_07_all_adverse_completion_intervals(self):
        rows = self.report["allocation"]
        self.assertEqual(len(rows), 9)
        self.assertEqual({(r["younger"], r["older"], r["stage"]) for r in rows},
                         {(young, old, stage) for young, old in EDGES for stage in NAMED})
        for row in rows:
            young, old, stage = row["younger"], row["older"], row["stage"]
            y, o = count(self.raw[young, stage]), count(self.raw[old, stage])
            ey, eo = exposure(self.raw[young, stage]), exposure(self.raw[old, stage])
            uy = count(self.raw[young, "unknown"]) + self.residual[young]
            uo = count(self.raw[old, "unknown"]) + self.residual[old]
            with self.subTest(age=(young, old), stage=stage):
                self.assert_close(row["rr_min"], F(o * ey, (y + uy) * eo))
                self.assert_close(row["rr_max"], F((o + uo) * ey, y * eo))
                lower_p = probability_from_rr(row["lower"], ey, eo)
                upper_p = probability_from_rr(row["upper"], ey, eo)
                self.assert_close(binomial_geq(y + uy + o, o, lower_p), ALPHA / (2 * FAMILY))
                self.assert_close(binomial_leq(y + o + uo, o + uo, upper_p), ALPHA / (2 * FAMILY))
                self.assertEqual(row["decline_survives"], row["upper"] < 1)
                nominal = next(r for r in self.report["contrasts"]
                               if (r["younger"], r["older"], r["stage"]) == (young, old, stage))
                self.assertLessEqual(row["lower"], nominal["lower"])
                self.assertGreaterEqual(row["upper"], nominal["upper"])

    def test_08_zero_count_equal_rate_and_reversal_boundaries(self):
        for y, o, ey, eo in ((0, 9, 100, 100), (9, 0, 100, 100),
                             (4, 4, 100, 100), (12, 18, 200, 300), (1, 2, 100, 700)):
            with self.subTest(counts=(y, o), exposures=(ey, eo)):
                result = A.conditional_rr(y, o, ey, eo)
                self.assert_interval_tails(result, y, o, ey, eo)
                reversed_result = A.conditional_rr(o, y, eo, ey)
                if result["lower"] > 0:
                    self.assert_close(reversed_result["upper"], 1 / result["lower"], rel=1e-10)
                if math.isfinite(result["upper"]):
                    self.assert_close(reversed_result["lower"], 1 / result["upper"], rel=1e-10)
                if y * eo == o * ey:
                    self.assertEqual(result["decision"], "unresolved")
                    self.assert_close(result["rr"], 1)
        empty = A.conditional_rr(0, 0, 100, 100)
        self.assertEqual(empty["lower"], 0)
        self.assertTrue(math.isinf(empty["upper"]))
        self.assertEqual(empty["decision"], "unresolved")

    def test_09_newline_hashing_and_portable_path_keys(self):
        actual = A.file_hashes()
        expected_paths = (RESEARCH / "PROTOCOL.md", CSV_PATH, RESEARCH / "analysis.py")
        expected = {path.relative_to(REPO).as_posix(): newline_hash(path) for path in expected_paths}
        self.assertEqual(actual, expected)
        self.assertEqual(self.report["hashes"], expected)
        for path in expected_paths:
            relative = path.relative_to(REPO)
            self.assertEqual(PureWindowsPath(*relative.parts).as_posix(), relative.as_posix())
            self.assertEqual(PurePosixPath(*relative.parts).as_posix(), relative.as_posix())
            self.assertNotIn("\\", relative.as_posix())
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "line-endings.txt"
            original = "alpha\nbeta \nnon-ASCII: \u03bb\n".encode("utf-8")
            expected_hash = hashlib.sha256(original).hexdigest()
            for separator in (b"\n", b"\r\n", b"\r"):
                path.write_bytes(original.replace(b"\n", separator))
                self.assertEqual(A.normalized_sha256(path), expected_hash)
            path.write_bytes(original + b"\n")
            self.assertNotEqual(A.normalized_sha256(path), expected_hash)
            path.write_bytes(original.replace(b"beta ", b"beta"))
            self.assertNotEqual(A.normalized_sha256(path), expected_hash)

    def test_10_entire_displayed_numerical_block_recomputes(self):
        text = (RESEARCH / "README.md").read_text(encoding="utf-8")
        begin, end = "<!-- BEGIN RECOMPUTED RESULTS -->", "<!-- END RECOMPUTED RESULTS -->"
        self.assertEqual(text.count(begin), 1)
        self.assertEqual(text.count(end), 1)
        actual = text.split(begin, 1)[1].split(end, 1)[0].strip()
        self.assertEqual(actual, A.render_results(self.report).strip())
        self.assertIn(PROTOCOL_COMMIT, text)
        # The app publishes an identical tree after an unauthenticated shell
        # push fails; preserve, rather than backdate, the original local object.
        original = text.split("<!-- BEGIN ORIGINAL PROTOCOL COMMIT -->", 1)[1]
        original = original.split("<!-- END ORIGINAL PROTOCOL COMMIT -->", 1)[0]
        commit_bytes = original.split("```text\n", 1)[1].split("```", 1)[0].encode("utf-8")
        header = b"commit " + str(len(commit_bytes)).encode("ascii") + b"\0"
        self.assertEqual(hashlib.sha1(header + commit_bytes).hexdigest(), PROTOCOL_COMMIT)
        self.assertIn("nonclinical", text.lower())
        self.assertIn("no screening or medical advice", text.lower())

    def test_11_stage_proportions_and_separate_one_sided_family(self):
        proportions = self.report["stage_proportions"]
        self.assertEqual(len(proportions), len(AGES))
        self.assertEqual({row["age"] for row in proportions}, set(AGES))
        for row in proportions:
            age = row["age"]
            total = count(self.raw[age, "overall"])
            self.assertEqual(row["distant_percent"], float(F(100 * count(self.raw[age, "distant"]), total)))
            self.assertEqual(row["unknown_percent"], float(F(100 * count(self.raw[age, "unknown"]), total)))
        expected_changes = [row for row in self.report["contrasts"]
                            if (row["decision"] == "decline") != row["one_sided_decline"]]
        self.assertEqual(self.report["one_sided_changes"], expected_changes)
        # These regression targets come from independent exact-tail decisions,
        # not from a rule promoting a sensitivity result to the primary family.
        declines = {(row["younger"], row["older"], row["stage"])
                    for row in self.report["contrasts"] if row["decision"] == "decline"}
        self.assertEqual(declines, {(y, o, stage) for y, o in EDGES[1:]
                                    for stage in ("overall", "localized")})
        changed = {(r["younger"], r["older"], r["stage"]) for r in expected_changes}
        self.assertEqual(changed, {("75-79", "80-84", "regional")})
        self.assertNotIn(("75-79", "80-84", "regional"), declines)

    def test_12_negative_residual_preserves_accounting_and_abstains_allocation(self):
        data = copy.deepcopy(self.data)
        # A failed partition in the middle age affects both touching edges.
        age = "75-79"
        data[age, "overall"]["cases"] = sum(data[age, s]["cases"] for s in STAGES[1:]) - 1
        report = A.analyze(data)
        self.assertEqual(A.residuals(data)[age], -1)
        for row in report["allocation"]:
            affected = age in (row["younger"], row["older"])
            self.assertEqual(row["allocation_valid"], not affected)
            if affected:
                for key in ("lower", "upper", "rr_min", "rr_max", "decline_survives"):
                    self.assertIsNone(row[key])
        for row in report["decomposition"]:
            if age in (row["younger"], row["older"]):
                self.assertIsNone(row["localized_share_min"])
                self.assertIsNone(row["localized_share_max"])
                self.assertFalse(row["allocation_majority"])
        self.assertEqual(len(report["contrasts"]), FAMILY)

    def test_13_changed_printed_rate_cannot_change_count_based_contrasts(self):
        data = copy.deepcopy(self.data)
        data["65-69", "overall"]["printed_rate"] = "417.9"
        report = A.analyze(data)
        cell = next(r for r in report["reconciliation"]
                    if (r["age"], r["stage"]) == ("65-69", "overall"))
        self.assertFalse(cell["rounding_matches"])
        self.assertEqual(report["contrasts"], self.report["contrasts"])
        self.assertEqual(report["allocation"], self.report["allocation"])
        self.assertEqual(report["decomposition"], self.report["decomposition"])

    def test_14_malformed_source_inputs_are_not_imputed_or_repaired(self):
        variants = {}
        for name, field, value in (
            ("missing_count", "cases", ""),
            ("fractional_count", "cases", "37.5"),
            ("negative_count", "cases", "-1"),
            ("nonpositive_exposure", "population", "0"),
            ("different_age_denominator", "population", "8889000"),
            ("reversed_interval", "printed_ci_lower", "999.9"),
            ("missing_source", "source_url", ""),
        ):
            rows = copy.deepcopy(self.rows)
            rows[0][field] = value
            variants[name] = rows
        variants["duplicate_cell"] = copy.deepcopy(self.rows) + [copy.deepcopy(self.rows[0])]
        variants["missing_cell"] = copy.deepcopy(self.rows[1:])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.csv"
            for name, rows in variants.items():
                with path.open("w", encoding="utf-8", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(self.rows[0]))
                    writer.writeheader()
                    writer.writerows(rows)
                with self.subTest(case=name), self.assertRaises(ValueError):
                    A.load_data(path)
        for args in ((-1, 1, 10, 10), (1.5, 1, 10, 10), (True, 1, 10, 10),
                     (1, 1, 0, 10), (1, 1, 10, -1)):
            with self.subTest(arguments=args), self.assertRaises(ValueError):
                A.conditional_rr(*args)


if __name__ == "__main__":
    unittest.main()
