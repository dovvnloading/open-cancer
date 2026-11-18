"""Independent offline checks of the nonclinical paired-call compatibility audit.

Printed counts are read from the CSV.  The hard-coded transcription guard is
not an alternative data source.  A Decimal binomial recurrence independently
inverts confidence tails, and integer weights independently recover exact
tests.  Passing software checks do not establish source truth or biology.
"""

from __future__ import annotations

import ast
import csv
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib.util
import io
from itertools import product
import json
import math
from pathlib import Path, PureWindowsPath
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch


REPO = Path(__file__).resolve().parents[1]
NOTE = REPO / "research" / "ctdna-paired-calls"
DATA = NOTE / "data" / "aggregates.csv"
SPEC = importlib.util.spec_from_file_location(
    "research_ctdna_paired_calls", NOTE / "analysis.py")
ANALYSIS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ANALYSIS
SPEC.loader.exec_module(ANALYSIS)

F = Fraction
ALPHA = F(1, 40)
MARGINAL_ALPHA = F(1, 20)
CAPS = (F(0), F(1, 200), F(1, 100), F(1, 50), F(3, 100), F(1, 20))
REVERSAL_CAPS = (F(0), F(1, 200), F(1, 100), F(1, 50), F(3, 100))
CELLS = ("++", "+-", "-+", "--")
PRINTED_COUNTS = {
    "guided_itt": 294,
    "successful_analysis": 291,
    "both_positive": 31,
    "week4_only": 8,
    "week7_only": 6,
    "any_positive": 45,
    "reported_negative": 246,
    "itt_excluded": 14,
    "itt_missing_week7": 6,
    "itt_ineligible": 3,
    "itt_withdrew": 5,
}
NORMALIZED_SHA256 = {
    Path("PROTOCOL.md").as_posix():
        "6e615dcb5eef9a745d4cbd68e03865455f880d2f56ba382b37e15ec9cbc43e93",
    (Path("data") / "aggregates.csv").as_posix():
        "c96bd1f345900bc5987fe4de3e5be7d957ca96fa75ea53198aaedaec09080310",
}


def normalized_bytes(path):
    return Path(path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def read_source_rows():
    return list(csv.DictReader(io.StringIO(
        normalized_bytes(DATA).decode("utf-8"), newline="")))


def integer_tail(x, n, probability):
    """Complementary lower-tail integer recurrence, independent of the helper."""
    p = F(probability)
    if not (0 <= x <= n+1 and 0 <= p <= 1):
        raise ValueError("Invalid binomial arguments")
    a, d = p.numerator, p.denominator
    coefficient, lower_numerator = 1, 0
    for j in range(x):
        lower_numerator += coefficient * a**j * (d-a)**(n-j)
        coefficient = coefficient*(n-j)//(j+1)
    denominator = d**n
    return F(denominator-lower_numerator, denominator)


def decimal_tail(x, n, p):
    """Recurrence starting at P(X=0), used only as the independent oracle."""
    if x == 0:
        return Decimal(1)
    if p == 0:
        return Decimal(0)
    if p == 1:
        return Decimal(1)
    mass = (1-p)**n
    total = Decimal(0)
    for j in range(n+1):
        if j >= x:
            total += mass
        if j < n:
            mass *= Decimal(n-j) * p / (Decimal(j+1) * (1-p))
    return total


@lru_cache(maxsize=None)
def independent_lower(x, n, alpha=ALPHA):
    if x == 0:
        return 0.0
    with localcontext() as context:
        context.prec = 65
        target = Decimal(alpha.numerator) / Decimal(alpha.denominator)
        low, high = Decimal(0), Decimal(1)
        for _ in range(96):
            mid = (low+high)/2
            if decimal_tail(x, n, mid) < target:
                low = mid
            else:
                high = mid
        return float((low+high)/2)


def independent_ci(x, n):
    return [independent_lower(x, n), 1-independent_lower(n-x, n)]


def independent_discordants(x, y):
    return integer_tail(x, x+y, F(1, 2))


def independent_completions(m):
    return [row for row in product(range(m+1), repeat=4) if sum(row) == m]


def expected_numeric_results(counts):
    """All inferential output is rebuilt from the CSV, never copied decimals."""
    n, total = counts["successful_analysis"], counts["guided_itt"]
    pp, x, y, nn = (counts[key] for key in (
        "both_positive", "week4_only", "week7_only", "reported_negative"))
    m, d = total-n, x+y
    q = F(x, n)
    q_ci = independent_ci(x, n)
    q_low = q_ci[0]
    p = independent_discordants(x, y)
    table = dict(zip(CELLS, (pp, x, y, nn)))

    successful = {
        "n": n,
        "table": table,
        "q": {
            "fraction": str(q), "value": float(q),
            "lower_95": independent_lower(x, n, MARGINAL_ALPHA),
            "lower_975": q_low, "ci_95": q_ci,
        },
        "marginals": {
            "week4_positive_count": pp+x,
            "week7_positive_count": pp+y,
            "week4_positive_fraction": str(F(pp+x, n)),
            "week7_positive_fraction": str(F(pp+y, n)),
            "week4_positive_rate": (pp+x)/n,
            "week7_positive_rate": (pp+y)/n,
            "difference_fraction": str(F(y-x, n)),
            "difference": (y-x)/n,
        },
        "discordant": {
            "n": d, "reverse": x, "forward": y,
            "fraction": str(F(x, d)) if d else None,
            "value": x/d if d else None,
            "ci_95": independent_ci(x, d) if d else [0.0, 1.0],
            "p_fraction": str(p), "p_value": float(p),
            "reject_primary": p < ALPHA,
            "reject_marginal": p < MARGINAL_ALPHA,
        },
        "caps": [],
        "tradeoff": [],
    }
    for cap in CAPS:
        tail = integer_tail(x, n, cap)
        successful["caps"].append({
            "cap_fraction": str(cap), "cap": float(cap),
            "tail_fraction": str(tail), "tail": float(tail),
            "exact_table_compatible": cap >= q,
            "reject": tail < ALPHA,
        })
    for b in REVERSAL_CAPS:
        minimum = max(F(0), q-b)
        lower = max(0.0, q_low-float(b))
        successful["tradeoff"].append({
            "b_fraction": str(b), "b": float(b),
            "combined_error_fraction": str(minimum),
            "combined_error": float(minimum),
            "equal_error_fraction": str(minimum/2),
            "equal_error": float(minimum/2),
            "lower_combined_error_975": lower,
            "lower_equal_error_975": lower/2,
        })

    assignments = independent_completions(m)
    completion_rows = []
    for allocation in assignments:
        completed_x, completed_y = x+allocation[1], y+allocation[2]
        c_q = F(completed_x, total)
        c_p = independent_discordants(completed_x, completed_y)
        completion_rows.append({
            "allocation": dict(zip(CELLS, allocation)),
            "q_fraction": str(c_q), "q": float(c_q),
            "lower_975": independent_lower(completed_x, total),
            "ci_95": independent_ci(completed_x, total),
            "discordant_p_fraction": str(c_p),
            "discordant_p": float(c_p),
            "discordant_reject": c_p < ALPHA,
        })
    low_q, high_q = F(x, total), F(x+m, total)
    tails = [F(row["discordant_p_fraction"]) for row in completion_rows]
    min_p, max_p = min(tails), max(tails)
    missing = {
        "n_total": total, "m": m,
        "q_bounds": {"fractions": [str(low_q), str(high_q)],
                     "values": [float(low_q), float(high_q)]},
        "lower_975_bounds": [independent_lower(x, total),
                             independent_lower(x+m, total)],
        "ci_95_outer": [independent_lower(x, total),
                        1-independent_lower(total-x-m, total)],
        "cap_sensitivity": [], "tradeoff": [],
        "completions": completion_rows,
        "directional_extrema": {
            "min_fraction": str(min_p), "max_fraction": str(max_p),
            "min_value": float(min_p), "max_value": float(max_p),
            "min_attainers": [row["allocation"] for row in completion_rows
                              if F(row["discordant_p_fraction"]) == min_p],
            "max_attainers": [row["allocation"] for row in completion_rows
                              if F(row["discordant_p_fraction"]) == max_p],
            "rejecting_completions": sum(t < ALPHA for t in tails),
        },
    }
    for cap in CAPS:
        candidate_tails = [integer_tail(x+a[1], total, cap) for a in assignments]
        rejected = sum(t < ALPHA for t in candidate_tails)
        scope = "every" if rejected == len(assignments) else (
            "some" if rejected else "none")
        minimum, maximum = min(candidate_tails), max(candidate_tails)
        missing["cap_sensitivity"].append({
            "cap_fraction": str(cap), "cap": float(cap),
            "rejection_scope": scope,
            "rejecting_completions": rejected,
            "total_completions": len(assignments),
            "min_tail_fraction": str(minimum), "max_tail_fraction": str(maximum),
            "min_tail": float(minimum), "max_tail": float(maximum),
        })
    for b in REVERSAL_CAPS:
        errors = [max(F(0), value-b) for value in (low_q, high_q)]
        missing["tradeoff"].append({
            "b_fraction": str(b), "b": float(b),
            "combined_error_bounds_fractions": list(map(str, errors)),
            "combined_error_bounds": list(map(float, errors)),
            "equal_error_bounds_fractions": [str(e/2) for e in errors],
            "equal_error_bounds": [float(e/2) for e in errors],
        })
    return {
        "counts": counts,
        "design": {
            "family_alpha": float(2*ALPHA), "per_test_alpha": float(ALPHA),
            "primary_cap_fraction": str(F(1, 100)),
            "source_access_date": "2026-11-18",
        },
        "paired_classification_supported": True,
        "successful": successful, "missing": missing,
        "exclusion_audit": {
            "printed_total": counts["itt_excluded"],
            "parts_total": sum(counts[key] for key in (
                "itt_missing_week7", "itt_ineligible", "itt_withdrew")),
            "matches": True,
        },
    }


def numeric_leaves(value, location=()):
    """Find every numerical/boolean leaf so added output cannot evade the oracle."""
    if isinstance(value, dict):
        answer = {}
        for key, child in value.items():
            answer.update(numeric_leaves(child, location+(key,)))
        return answer
    if isinstance(value, (list, tuple)):
        answer = {}
        for index, child in enumerate(value):
            answer.update(numeric_leaves(child, location+(index,)))
        return answer
    if isinstance(value, (int, float)):
        return {location: value}
    if isinstance(value, str):
        try:
            F(value)
        except (ValueError, ZeroDivisionError):
            return {}
        return {location: value}
    return {}


class PairedCallsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.perf_counter()
        cls.rows = read_source_rows()
        cls.counts = {row["id"]: int(row["count"]) for row in cls.rows}
        cls.results = ANALYSIS.compute_results()
        cls.expected = expected_numeric_results(cls.counts)

    def assert_tree(self, expected, actual, location="result"):
        """Compare the independent numerical tree, allowing prose annotations."""
        if isinstance(expected, dict):
            self.assertIsInstance(actual, dict, location)
            for key, value in expected.items():
                self.assertIn(key, actual, location)
                self.assert_tree(value, actual[key], f"{location}.{key}")
        elif isinstance(expected, (list, tuple)):
            self.assertEqual(len(expected), len(actual), location)
            for index, (left, right) in enumerate(zip(expected, actual)):
                self.assert_tree(left, right, f"{location}[{index}]")
        elif isinstance(expected, float):
            self.assertIsInstance(actual, (int, float), location)
            self.assertTrue(math.isfinite(actual), location)
            self.assertAlmostEqual(expected, actual, delta=4e-13, msg=location)
        else:
            self.assertEqual(expected, actual, location)

    def test_csv_transcription_guards_and_provenance(self):
        self.assertEqual(len(self.rows), len(self.counts))
        self.assertEqual(self.counts, PRINTED_COUNTS)
        self.assertEqual(ANALYSIS.load_counts(), self.counts)
        for row in self.rows:
            with self.subTest(row=row["id"]):
                self.assertTrue(all(str(value).strip() for value in row.values()))
                self.assertEqual(row["source_url"],
                                 "https://pmc.ncbi.nlm.nih.gov/articles/PMC12974608/")
                self.assertEqual(row["doi"], "10.1038/s41591-025-03579-w")
                self.assertEqual(row["access_date"], "2026-11-18")
                self.assertIn("printed aggregate", row["unit"])
                self.assertIn("exclusive", row["article_rights"])
                self.assertIn("no Creative Commons", row["article_rights"])
                self.assertIn("factual counts", row["reuse_basis"])
                self.assertIn("not a publisher reuse grant", row["reuse_basis"])
                self.assertTrue(row["rights_url"].startswith("https://www.nature.com/"))
                self.assertTrue(row["reuse_basis_url"].startswith(
                    "https://www.copyright.gov/"))
                self.assertIn("Figure 1 caption" if row["id"].startswith("itt_")
                              else ("Extended Data Table 1" if row["id"] ==
                                    "reported_negative" else "first paragraph"),
                              row["locator"])
        negative = next(row for row in self.rows if row["id"] == "reported_negative")
        self.assertIn("not a directly printed negative/negative", negative["notes"])
        self.assertIn("both-negative definition", negative["notes"])
        self.assertEqual(self.counts["any_positive"], sum(
            self.counts[key] for key in ("both_positive", "week4_only", "week7_only")))
        self.assertEqual(self.counts["successful_analysis"],
                         self.counts["any_positive"] + self.counts["reported_negative"])
        self.assertEqual(self.counts["itt_excluded"], sum(self.counts[key] for key in (
            "itt_missing_week7", "itt_ineligible", "itt_withdrew")))

    def test_every_reported_numeric_result_from_csv(self):
        self.assert_tree(self.expected, self.results)
        self.assertEqual(set(numeric_leaves(self.expected)),
                         set(numeric_leaves(self.results)),
                         "Every numerical output must have an independent CSV-derived oracle")

    def test_cap_verdicts_preserve_nonrejection_and_point_mismatch(self):
        successful = self.results["successful"]
        for expected, actual in zip(self.expected["successful"]["caps"],
                                    successful["caps"]):
            verdict = "falsified under stated caps" if expected["reject"] else (
                "compatible with stated latitude")
            self.assertEqual(actual["verdict"], verdict)
        primary = next(row for row in successful["caps"]
                       if F(row["cap_fraction"]) == F(1, 100))
        self.assertEqual(successful["primary_verdict"], primary["verdict"])
        self.assertTrue(any(not row["exact_table_compatible"] and not row["reject"]
                            for row in successful["caps"]),
                        "An empirical point mismatch must not be promoted to rejection")
        self.assertEqual(successful["discordant"]["verdict"],
                         "compatible with stated latitude")
        for key in ("individual_availability_audit",
                    "unconditional_biological_no_intervention_claim",
                    "independently_calibrated_error_or_reversal_rates",
                    "assay_error_shedding_or_biological_clearance_attribution",
                    "recurrence_sensitivity", "schedule_optimum", "treatment_decision"):
            self.assertEqual(self.results["verdicts"][key], "untestable")

    def test_source_gate_withholds_unsupported_full_table(self):
        result = ANALYSIS.compute_results(paired_classification_supported=False)
        self.assertIsNone(result["successful"])
        self.assertIsNone(result["missing"])
        self.assertEqual(result["verdicts"]["paired_reporting_gate"], "untestable")
        discordant = result["source_supported_discordants"]
        self.assertEqual(discordant["verdict"], "untestable")
        self.assert_tree(self.expected["successful"]["discordant"], discordant)
        text = ANALYSIS.format_results(result)
        self.assertIn("Paired-availability gate: untestable", text)
        self.assertNotIn("Exact plug-in total-latitude frontier", text)
        with self.assertRaises(ValueError):
            ANALYSIS.compute_results(paired_classification_supported="yes")

    def test_absent_malformed_and_inconsistent_aggregates_are_rejected(self):
        examples = []
        for key, bad_value in (("guided_itt", 0), ("successful_analysis", 0),
                               ("reported_negative", -1), ("week4_only", 8.0),
                               ("week4_only", True), ("any_positive", 0),
                               ("itt_excluded", 0)):
            changed = dict(self.counts)
            changed[key] = bad_value
            examples.append(changed)
        absent = dict(self.counts)
        del absent["reported_negative"]
        examples.append(absent)
        extra = dict(self.counts, invented_negative=0)
        examples.append(extra)
        wrong_denominator = dict(self.counts)
        wrong_denominator["successful_analysis"] = self.counts["any_positive"]
        examples.append(wrong_denominator)
        for counts in examples:
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                ANALYSIS.compute_results(counts)

        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "bad.csv"
            variants = [self.rows+self.rows[:1], self.rows[:-1]]
            for field, value in (("count", "-1"), ("count", "1.5"),
                                 ("count", "١"), ("locator", ""), ("notes", " ")):
                changed_rows = [dict(row) for row in self.rows]
                changed_rows[0][field] = value
                variants.append(changed_rows)
            for rows in variants:
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(self.rows[0]))
                    writer.writeheader()
                    writer.writerows(rows)
                with self.subTest(rows=rows), self.assertRaises(ValueError):
                    ANALYSIS.load_counts(path)
            path.write_text("id,count\nguided_itt,294\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                ANALYSIS.load_counts(path)

    def test_exact_tail_oracle_and_confidence_endpoints(self):
        n, x = self.counts["successful_analysis"], self.counts["week4_only"]
        for cap in CAPS:
            self.assertEqual(ANALYSIS.binomial_tail(x, n, cap), integer_tail(x, n, cap))
        for nn in (1, 4, 11):
            for xx in range(nn+1):
                for probability in (F(0), F(1, 7), F(1, 2), F(6, 7), F(1)):
                    with self.subTest(x=xx, n=nn, p=probability):
                        self.assertEqual(ANALYSIS.binomial_tail(xx, nn, probability),
                                         integer_tail(xx, nn, probability))
        cases = ((0, n), (x, n), (n, n), (x, x+self.counts["week7_only"]))
        for xx, nn in cases:
            self.assert_tree(independent_ci(xx, nn), ANALYSIS.cp_interval(xx, nn))
            for alpha in (ALPHA, MARGINAL_ALPHA):
                lower = ANALYSIS.lower_bound(xx, nn, alpha)
                self.assertAlmostEqual(lower, independent_lower(xx, nn, alpha),
                                       delta=4e-13)
                if xx:
                    with localcontext() as context:
                        context.prec = 65
                        observed = decimal_tail(xx, nn, Decimal(str(lower)))
                    self.assertAlmostEqual(float(observed), float(alpha), delta=2e-11)
                else:
                    self.assertEqual(lower, 0.0)

    def test_missingness_is_exhaustive_not_missing_at_random(self):
        m = self.counts["guided_itt"]-self.counts["successful_analysis"]
        rows = independent_completions(m)
        self.assertEqual(sorted(ANALYSIS.completions(m)), rows)
        self.assertEqual(len(rows), math.comb(m+3, 3))
        bounds = self.results["missing"]["q_bounds"]["fractions"]
        q_values = [F(row["q_fraction"]) for row in self.results["missing"]["completions"]]
        self.assertEqual([min(q_values), max(q_values)], list(map(F, bounds)))
        for row in self.results["missing"]["completions"]:
            allocation = row["allocation"]
            self.assertEqual(sum(allocation.values()), m)
            self.assertTrue(all(isinstance(value, int) and value >= 0
                                for value in allocation.values()))
            self.assertGreaterEqual(row["ci_95"][0],
                                    self.results["missing"]["ci_95_outer"][0]-1e-14)
            self.assertLessEqual(row["ci_95"][1],
                                 self.results["missing"]["ci_95_outer"][1]+1e-14)

    def test_paired_direction_uses_both_discordants(self):
        x, y = self.counts["week4_only"], self.counts["week7_only"]
        p = independent_discordants(x, y)
        self.assertEqual(ANALYSIS.discordant_pvalue(x, y), p)
        self.assertEqual(ANALYSIS.discordant_pvalue(0, 0), F(1))
        self.assertFalse(self.results["successful"]["discordant"]["reject_primary"])
        self.assertFalse(self.results["successful"]["discordant"]["reject_marginal"])
        self.assertEqual(self.results["missing"]["directional_extrema"]
                         ["rejecting_completions"], 0)
        # The null is on the discordant probability, never x/n or x/(++ + +-).
        self.assertNotEqual(p, integer_tail(x, self.counts["successful_analysis"], F(1, 2)))
        self.assertNotEqual(F(x, x+y), F(x, self.counts["both_positive"]+x))

    def test_sharpness_cap_constructions_preserve_the_full_table(self):
        n = self.counts["successful_analysis"]
        probabilities = {key: F(value, n) for key, value in
                         self.results["successful"]["table"].items()}
        q = probabilities["+-"]
        cases = ((q, F(0), F(0)), (F(0), q, F(0)), (F(0), F(0), q),
                 (q/3, q/3, q/3), (q/2, q/2, F(0)), (F(1), F(1), F(1)))
        for caps in cases:
            with self.subTest(caps=caps):
                a4, a7, ab = ANALYSIS.sharpness_allocation(q, *caps)
                self.assertEqual(a4+a7+ab, q)
                self.assertTrue(all(0 <= value <= cap for value, cap in
                                    zip((a4, a7, ab), caps)))
                joint = [((1, 0), (0, 0), a4), ((1, 0), (1, 1), a7),
                         ((1, 0), (1, 0), ab)]
                for cell, yy in (("++", (1, 1)), ("-+", (0, 1)), ("--", (0, 0))):
                    joint.append((yy, yy, probabilities[cell]))
                observed = {cell: sum(weight for yy, zz, weight in joint
                                     if yy == tuple(1 if sign == "+" else 0
                                                    for sign in cell))
                            for cell in CELLS}
                self.assertEqual(observed, probabilities)
                actual = (sum(weight for yy, zz, weight in joint if yy[0] != zz[0]),
                          sum(weight for yy, zz, weight in joint if yy[1] != zz[1]),
                          sum(weight for yy, zz, weight in joint if zz == (1, 0)))
                self.assertEqual(actual, (a4, a7, ab))
        self.assertIsNone(ANALYSIS.sharpness_allocation(q, q/4, q/4, q/4))
        self.assertEqual(ANALYSIS.sharpness_allocation(F(0), F(0), F(0), F(0)),
                         (F(0), F(0), F(0)))

    def test_synthetic_observation_counterexamples(self):
        # All latent states persist at 11, but always calling 10 makes q>r.
        # Thus persistence alone cannot justify the directional restriction.
        latent, observed = (1, 1), (1, 0)
        self.assertLessEqual(latent[0], latent[1])
        e4, e7 = int(observed[0] != latent[0]), int(observed[1] != latent[1])
        self.assertEqual((e4, e7), (0, 1))
        self.assertGreater(int(observed == (1, 0)), int(observed == (0, 1)))
        # A shared increasing noisy response can have positive reversal mass.
        for p4, p7 in ((F(1, 4), F(1, 2)), (F(1, 2), F(1, 2))):
            reverse = p4*(1-p7)
            forward = (1-p4)*p7
            self.assertGreater(reverse, 0)
            self.assertLessEqual(reverse, forward)
            self.assertEqual(forward-reverse, p7-p4)

    def test_newline_hashes_and_posix_path_keys(self):
        repository_keys = {(NOTE / key).relative_to(REPO).as_posix(): value
                           for key, value in NORMALIZED_SHA256.items()}
        self.assertEqual(self.results["hashes"], repository_keys)
        for key, expected in NORMALIZED_SHA256.items():
            source = NOTE / key
            self.assertNotIn("\\", key)
            self.assertEqual(Path(*PureWindowsPath(key).parts).as_posix(), key)
            self.assertEqual(hashlib.sha256(normalized_bytes(source)).hexdigest(), expected)
            self.assertEqual(ANALYSIS.normalized_sha256(source), expected)
            with tempfile.TemporaryDirectory() as temp:
                path = Path(temp) / "text.txt"
                raw = normalized_bytes(source)
                for ending in (b"\n", b"\r\n", b"\r"):
                    path.write_bytes(raw.replace(b"\n", ending))
                    self.assertEqual(ANALYSIS.normalized_sha256(path), expected)
                    if source == DATA:
                        self.assertEqual(ANALYSIS.load_counts(path), self.counts)

    def test_analysis_has_only_offline_standard_library_dependencies(self):
        tree = ast.parse((NOTE / "analysis.py").read_text(encoding="utf-8"))
        roots = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                roots.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                roots.add(node.module.split(".")[0])
        self.assertTrue(roots <= sys.stdlib_module_names, roots-sys.stdlib_module_names)
        self.assertFalse(roots & {"socket", "urllib", "http", "ftplib", "requests"})
        with patch("socket.socket", side_effect=AssertionError("Network prohibited")), \
                patch("socket.create_connection", side_effect=AssertionError("Network prohibited")):
            self.assert_tree(self.expected, ANALYSIS.compute_results())

    def test_readme_result_block_and_cli_recompute(self):
        readme = (NOTE / "README.md").read_text(encoding="utf-8")
        start, end = "<!-- BEGIN RECOMPUTED RESULTS -->", "<!-- END RECOMPUTED RESULTS -->"
        self.assertEqual(readme.count(start), 1)
        self.assertEqual(readme.count(end), 1)
        retained = readme.split(start, 1)[1].split(end, 1)[0].strip()
        self.assertEqual(retained, ANALYSIS.format_results(self.results).strip())
        for option in ("--check", "--json"):
            run = subprocess.run([sys.executable, "-B", str(NOTE / "analysis.py"), option],
                                 cwd=REPO, capture_output=True, text=True, timeout=20,
                                 check=False)
            self.assertEqual(run.returncode, 0, run.stdout+run.stderr)
            if option == "--json":
                self.assert_tree(self.expected, json.loads(run.stdout))

    def test_zz_execution_budget(self):
        self.assertLess(time.perf_counter()-self.started, 60.0,
                        "The entire offline test file must finish within one minute")


if __name__ == "__main__":
    unittest.main()
