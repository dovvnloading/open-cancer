"""Offline checks of a nonclinical, preclinical aggregate analysis.

The independent oracle uses factorial multinomial weights, not the analysis
module's column-wise combination formula.  Synthetic mathematical fixtures
test the software only; passing checks do not adjudicate biological causes.
"""

from __future__ import annotations

import ast
from collections import Counter
import csv
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib.util
import io
from itertools import product
import math
from pathlib import Path, PureWindowsPath
import sys
import tempfile
import time
import unittest
from unittest.mock import patch


REPO = Path(__file__).resolve().parents[1]
NOTE = REPO / "research" / "schedule-terminal-cause"
SPEC = importlib.util.spec_from_file_location(
    "research_schedule_terminal_cause", NOTE / "analysis.py")
ANALYSIS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ANALYSIS
SPEC.loader.exec_module(ANALYSIS)

F = Fraction
ALPHA = F(1, 20)
CELL_ALPHA = F(1, 240)
CAUSES = ("toxicity", "cancer", "other")
GROUPS = ("vehicle", "ping_pong", "tandem", "combo_mtd")
PAIRS = (("combo_mtd", "ping_pong"), ("combo_mtd", "tandem"),
         ("combo_mtd", "vehicle"), ("tandem", "ping_pong"),
         ("ping_pong", "vehicle"), ("tandem", "vehicle"))

# These are source-transcription guards, not substitutes for the CSV inputs.
PRINTED_COUNTS = {
    "cap_adaptive": (4, 5, 4), "cap_mtd": (2, 4, 0),
    "gem_adaptive": (7, 1, 5), "gem_mtd": (5, 1, 0),
    "vehicle": (1, 3, 0), "ping_pong": (3, 6, 3),
    "tandem": (5, 4, 3), "combo_mtd": (6, 0, 0),
}
PRINTED_DOSES = {
    ("gem_dm", "GEM"): ("256", "5.8"),
    ("gem_intermittent", "GEM"): ("33.3", "5.5"),
    ("gem_mtd", "GEM"): ("783.3", "12.8"),
    ("cap_dm", "CAP"): ("533.8", "11.8"),
    ("cap_intermittent", "CAP"): ("1226.6", "26"),
    ("cap_mtd", "CAP"): ("1373.3", "28.4"),
    ("ping_dm", "GEM"): ("379.1", "4.9"),
    ("ping_dm", "CAP"): ("475.5", "6.1"),
    ("tandem_dm", "GEM"): ("229", "4.3"),
    ("tandem_dm", "CAP"): ("492.1", "8.7"),
    ("ping_intermittent", "GEM"): ("550", "6.4"),
    ("ping_intermittent", "CAP"): ("660", "8"),
    ("tandem_intermittent", "GEM"): ("341.6", "9"),
    ("tandem_intermittent", "CAP"): ("606.6", "12.5"),
    ("combo_mtd", "GEM"): ("433.3", "20.6"),
    ("combo_mtd", "CAP"): ("766.6", "30.6"),
}
NORMALIZED_SHA256 = {
    Path("PROTOCOL.md").as_posix():
        "0f228ecb01a7d187b37a57f43d163fb8b4cedb954ee0258c24ccd94363283b90",
    (Path("data") / "allocations.csv").as_posix():
        "fb9f938da5d3e61672f6ff8fccc2cab8c278530aff7a5449105c880d97eb9ac4",
    (Path("data") / "counts.csv").as_posix():
        "f79fc1a13ebf7b8c36cf480491822c5d523c40b9488f2acf260818e74a97940b",
    (Path("data") / "dose_conflicts.csv").as_posix():
        "c9af39a48d9254fc6c63debf8cf48963024dae84aa37b759f4e72c2e9a5ac3c5",
    (Path("data") / "doses.csv").as_posix():
        "826de1f5b7c456745922a69644339006239e12b8a24b0c8db96d85cf53686093",
}


def normalized_bytes(path):
    return Path(path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def read_csvs():
    return {path.stem: list(csv.DictReader(io.StringIO(
                normalized_bytes(path).decode("utf-8-sig"), newline="")))
            for path in sorted((NOTE / "data").glob("*.csv"))}


def allocations(total, caps):
    """All bounded integer compositions, including empty-column cases."""
    if len(caps) == 1:
        if 0 <= total <= caps[0]:
            yield (total,)
        return
    for first in range(max(0, total - sum(caps[1:])), min(total, caps[0]) + 1):
        for rest in allocations(total - first, caps[1:]):
            yield (first,) + rest


@lru_cache(maxsize=None)
def factorial_exact(first, second):
    """Independent conditional probability-ordering test, using row weights."""
    n1, n2 = sum(first), sum(second)
    columns = tuple(a + b for a, b in zip(first, second))
    denominator = math.factorial(n1 + n2) // math.prod(
        math.factorial(x) for x in columns)
    row_numerator = math.factorial(n1) * math.factorial(n2)

    def weight(row):
        return row_numerator // math.prod(
            math.factorial(x) * math.factorial(c - x)
            for x, c in zip(row, columns))

    observed = weight(first)
    weights = [weight(row) for row in allocations(n1, columns)]
    if sum(weights) != denominator:
        raise AssertionError("Independent exact distribution does not normalize")
    return F(sum(w for w in weights if w <= observed), denominator)


def independent_tests(first, second):
    return (factorial_exact(first, second),
            factorial_exact((first[0], sum(first[1:])),
                            (second[0], sum(second[1:]))))


def independent_holm(values):
    """Closed-form prefix maximum, allowing arbitrary dependent hypotheses."""
    ranked = sorted(range(len(values)), key=values.__getitem__)
    answer = [None] * len(values)
    for rank, index in enumerate(ranked):
        answer[index] = min(F(1), max(
            (len(values) - j) * values[ranked[j]] for j in range(rank + 1)))
    return tuple(answer)


def independent_row_states(counts, include_other):
    t, c, o = counts
    return tuple(((t - dt, c + dt + do, o - do), dt + do)
                 for dt in range(t + 1)
                 for do in range((o if include_other else 0) + 1))


@lru_cache(maxsize=None)
def independent_evaluation(state):
    rows = dict(zip(GROUPS, state))
    raw = tuple(value for a, b in PAIRS
                for value in independent_tests(rows[a], rows[b]))
    adjusted = independent_holm(raw)
    differences = tuple(tuple(F(x, sum(rows[a])) - F(y, sum(rows[b]))
                              for x, y in zip(rows[a], rows[b]))
                        for a, b in PAIRS)
    return differences, raw, adjusted


def independent_states(source, include_other):
    alternatives = [independent_row_states(source[g], include_other) for g in GROUPS]
    for selection in product(*alternatives):
        state = tuple(item[0] for item in selection)
        yield state, sum(item[1] for item in selection), independent_evaluation(state)


def binomial_mass(k, n, p):
    return math.factorial(n) / (math.factorial(k) * math.factorial(n - k)) * p**k * (1 - p)**(n - k)


def binomial_tail(k, n, p, upper):
    support = range(k, n + 1) if upper else range(k + 1)
    return math.fsum(binomial_mass(j, n, p) for j in support)


@lru_cache(maxsize=None)
def independent_interval(k, n):
    """Invert binomial tails using independently computed factorial masses."""
    target = float(CELL_ALPHA / 2)
    lower, upper = 0.0, 1.0
    if k:
        lo, hi = 0.0, 1.0
        for _ in range(64):
            middle = (lo + hi) / 2
            if binomial_tail(k, n, middle, True) < target:
                lo = middle
            else:
                hi = middle
        lower = (lo + hi) / 2
    if k < n:
        lo, hi = 0.0, 1.0
        for _ in range(64):
            middle = (lo + hi) / 2
            if binomial_tail(k, n, middle, False) > target:
                lo = middle
            else:
                hi = middle
        upper = (lo + hi) / 2
    return lower, upper


def independent_metrics(state, evaluation=None):
    differences, raw, adjusted = evaluation or independent_evaluation(state)
    rows = dict(zip(GROUPS, state))
    result = []
    for i, (first, second) in enumerate(PAIRS):
        metrics = {f"{cause}_difference": value
                   for cause, value in zip(CAUSES, differences[i])}
        metrics.update(composition_p=raw[2 * i], toxicity_p=raw[2 * i + 1],
                       composition_holm=adjusted[2 * i], toxicity_holm=adjusted[2 * i + 1])
        for j, cause in enumerate(CAUSES):
            left = independent_interval(rows[first][j], sum(rows[first]))
            right = independent_interval(rows[second][j], sum(rows[second]))
            metrics[f"{cause}_interval_lower"] = left[0] - right[1]
            metrics[f"{cause}_interval_upper"] = left[1] - right[0]
        result.append(metrics)
    return result


def lost_criteria(metrics, baseline):
    direction = baseline["toxicity_difference"]
    lost_direction = not (direction * metrics["toxicity_difference"] > 0)
    return {
        "raw_composition": metrics["composition_p"] > ALPHA,
        "raw_toxicity": metrics["toxicity_p"] > ALPHA,
        "holm_composition": metrics["composition_holm"] > ALPHA,
        "holm_toxicity": metrics["toxicity_holm"] > ALPHA,
        "direction": lost_direction,
        "combined": lost_direction or metrics["composition_holm"] > ALPHA
                    or metrics["toxicity_holm"] > ALPHA,
    }


class ScheduleTerminalCauseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.monotonic()
        cls.tables = read_csvs()
        cls.count_rows = {row["group_id"]: row for row in cls.tables["counts"]}
        cls.counts = {key: tuple(int(row[c]) for c in CAUSES)
                      for key, row in cls.count_rows.items()}
        with patch("socket.socket", side_effect=AssertionError("Network access forbidden")):
            cls.report = ANALYSIS.build_report()
        cls.oracle_states = {
            include_other: list(independent_states(cls.counts, include_other))
            for include_other in (False, True)
        }
        cls.oracle_metrics = {state: independent_metrics(state, evaluation)
                              for state, _, evaluation in cls.oracle_states[True]}
        cls.baseline = tuple(cls.counts[group] for group in GROUPS)
        cls.baseline_metrics = cls.oracle_metrics[cls.baseline]

    def test_aggregate_transcriptions_and_denominators_are_not_repaired(self):
        self.assertEqual(self.counts, PRINTED_COUNTS)
        arms = {row["arm_id"]: row for row in self.tables["allocations"]}
        self.assertEqual(len(arms), len(self.tables["allocations"]))
        self.assertEqual(len(arms), 12)
        totals = Counter()
        for row in arms.values():
            expected = 4 if row["arm_id"] == "vehicle" else 6
            self.assertEqual(int(row["animals"]), expected)
            totals[row["pooled_group_id"]] += int(row["animals"])
        self.assertEqual(sum(totals.values()), 70)
        self.assertEqual(sum(totals[group] for group in GROUPS), 34)
        self.assertEqual({key: sum(value) - totals[key] for key, value in self.counts.items()},
                         {key: int(key in {"cap_adaptive", "gem_adaptive"}) for key in self.counts})
        for key, row in self.count_rows.items():
            self.assertEqual({a for a in row["allocation_arm_ids"].split("|")},
                             {a for a, record in arms.items() if record["pooled_group_id"] == key})
        self.assertEqual(sum(row["group_id"] == "combo_mtd" for row in self.tables["counts"]), 1)
        self.assertEqual(sum(row["group_id"] == "vehicle" for row in self.tables["counts"]), 1)

    def test_printed_dose_means_and_source_conflicts_are_preserved(self):
        observed = {(row["arm_id"], row["drug"]):
                    (Decimal(row["mean_cumulative_mg_kg"]), Decimal(row["mean_daily_mg_kg_day"]))
                    for row in self.tables["doses"]}
        expected = {key: tuple(Decimal(value) for value in pair) for key, pair in PRINTED_DOSES.items()}
        self.assertEqual(observed, expected)
        self.assertEqual(len(observed), len(self.tables["doses"]))
        conflicts = {(row["arm_id"], row["drug"]): Decimal(row["mean_cumulative_mg_kg"])
                     for row in self.tables["dose_conflicts"]}
        self.assertEqual(conflicts, {("gem_dm", "GEM"): Decimal("256.9"),
                                     ("gem_intermittent", "GEM"): Decimal("333.3")})
        for key, value in conflicts.items():
            self.assertNotEqual(observed[key][0], value)

    def test_every_csv_row_has_exact_locator_and_reuse_provenance(self):
        forbidden = {"mouse_id", "patient_id", "subject_id", "trajectory", "individual_id"}
        for table, rows in self.tables.items():
            self.assertTrue(rows)
            for index, row in enumerate(rows, 2):
                with self.subTest(table=table, row=index):
                    self.assertNotIn(None, row)
                    self.assertFalse(forbidden & set(row))
                    for key in ("source_url", "source_table", "source_row", "access_date", "reuse_basis"):
                        self.assertTrue(row[key].strip())
                    self.assertTrue(row["source_url"].startswith("https://"))
                    self.assertEqual(row["access_date"], "2026-11-18")
                    self.assertIn("CC BY 4.0", row["reuse_basis"])
                    self.assertTrue(row.get("pdf_page") or row.get("source_page"))
                    self.assertTrue(any(value for key, value in row.items() if key.endswith("column")))
        for row in self.tables["counts"]:
            self.assertIn(row["source_table"], {"S2", "S3", "S4"})
            self.assertIn(row["pdf_page"], {"32", "33"})
            self.assertEqual(tuple(row[f"{cause}_column"] for cause in CAUSES),
                             ("Died of toxicity", "Died of cancer", "Died of something else"))
        for row in self.tables["doses"]:
            self.assertEqual((row["source_table"], row["pdf_page"]), ("S5", "34"))

    def test_probability_ordering_exact_tests_use_ties_and_zero_columns(self):
        fixtures = (((0, 3, 0), (0, 7, 0)), ((2, 0, 1), (0, 2, 1)),
                    ((1, 0, 0), (0, 1, 0)), ((0, 0, 5), (4, 0, 0)),
                    ((2, 2, 2), (2, 2, 2)), ((3, 6, 3), (6, 0, 0)))
        for first, second in fixtures:
            with self.subTest(first=first, second=second):
                expected = factorial_exact(first, second)
                self.assertEqual(ANALYSIS.exact_conditional(first, second), expected)
                self.assertEqual(ANALYSIS.exact_conditional(second, first), expected)
                self.assertEqual(ANALYSIS.exact_conditional(first[::-1], second[::-1]), expected)
                self.assertEqual(tuple(ANALYSIS.pair_tests(first, second)), independent_tests(first, second))
        self.assertEqual(ANALYSIS.exact_conditional((1, 0), (0, 1)), F(1))
        self.assertEqual(ANALYSIS.exact_conditional((0, 3, 0), (0, 7, 0)), F(1))
        for first in product(range(3), repeat=3):
            if not sum(first):
                continue
            second = (2, 1, 2)
            self.assertEqual(ANALYSIS.exact_conditional(first, second), factorial_exact(first, second))

    def test_holm_is_a_fixed_twelve_test_family_with_arbitrary_dependence(self):
        fixtures = ((F(1, 100), F(1, 100), F(1), F(0)),
                    tuple(F(k, 200) for k in (4, 1, 9, 7, 2, 3, 8, 6, 12, 10, 11, 5)),
                    (F(1, 240),) * 12, (F(0),) * 12, (F(1),) * 12)
        for raw in fixtures:
            self.assertEqual(tuple(ANALYSIS.holm(raw)), independent_holm(raw))
        baseline = tuple(self.counts[group] for group in GROUPS)
        _, raw, adjusted = independent_evaluation(baseline)
        self.assertEqual(len(raw), 12)
        self.assertTrue(all(value >= p for value, p in zip(adjusted, raw)))

    def test_clopper_pearson_bounds_solve_the_correct_tails(self):
        target = float(CELL_ALPHA / 2)
        for n in sorted({sum(self.counts[group]) for group in GROUPS}):
            for k in range(n + 1):
                with self.subTest(k=k, n=n):
                    lower, upper = ANALYSIS.clopper_pearson(k, n, CELL_ALPHA)
                    self.assertLessEqual(lower, k / n)
                    self.assertLessEqual(k / n, upper)
                    if k == 0:
                        self.assertEqual(lower, 0)
                        self.assertAlmostEqual(upper, 1 - target**(1 / n), places=12)
                    else:
                        self.assertAlmostEqual(binomial_tail(k, n, lower, True), target, places=12)
                    if k == n:
                        self.assertEqual(upper, 1)
                        self.assertAlmostEqual(lower, target**(1 / n), places=12)
                    else:
                        self.assertAlmostEqual(binomial_tail(k, n, upper, False), target, places=12)
                    mirrored = ANALYSIS.clopper_pearson(n - k, n, CELL_ALPHA)
                    self.assertAlmostEqual(lower, 1 - mirrored[1], places=12)
                    self.assertAlmostEqual(upper, 1 - mirrored[0], places=12)

    def assert_number(self, actual, expected):
        if isinstance(expected, float):
            self.assertAlmostEqual(actual, expected, places=12)
        else:
            self.assertEqual(actual, expected)

    def assert_metric_map(self, actual, expected):
        self.assertEqual(set(actual), set(expected))
        for name, value in expected.items():
            self.assert_number(actual[name], value)

    def test_every_baseline_number_is_recomputed_from_the_csvs(self):
        report, s4 = self.report, self.report["s4"]
        self.assertEqual(report["protocol_commit"], "dbf3c54746ba853a81a65f177108531081460978")
        self.assertEqual((report["alpha"], report["cell_alpha"]), (ALPHA, CELL_ALPHA))
        self.assertEqual(s4["state"], self.baseline)
        self.assertEqual(s4["total_n"], sum(map(sum, self.baseline)))
        for report_key, table in (("counts_rows", "counts"), ("allocation_rows", "allocations"),
                                  ("doses", "doses"), ("dose_conflicts", "dose_conflicts")):
            self.assertEqual(report[report_key], self.tables[table])
        allocations_by_group = Counter()
        for row in self.tables["allocations"]:
            allocations_by_group[row["pooled_group_id"]] += int(row["animals"])
        self.assertEqual(report["allocated_total"], sum(allocations_by_group.values()))
        self.assertEqual(report["printed_total"], sum(map(sum, self.counts.values())))
        expected_denominators = [
            {"group_id": key, "counts": self.counts[key], "allocated": allocations_by_group[key],
             "printed": sum(self.counts[key]),
             "excess": sum(self.counts[key]) - allocations_by_group[key],
             "table": row["source_table"], "page": row["pdf_page"]}
            for key, row in self.count_rows.items()
        ]
        self.assertEqual(report["denominators"], expected_denominators)
        self.assertEqual(len(s4["arms"]), len(GROUPS))
        for arm, group in zip(s4["arms"], GROUPS):
            counts = self.counts[group]
            self.assertEqual((arm["group_id"], arm["counts"], arm["n"]),
                             (group, counts, sum(counts)))
            self.assertEqual(arm["proportions"], tuple(F(x, sum(counts)) for x in counts))
            for actual, k in zip(arm["intervals"], counts):
                expected = independent_interval(k, sum(counts))
                self.assert_number(actual[0], expected[0])
                self.assert_number(actual[1], expected[1])
        self.assertEqual(len(s4["baseline"]), len(PAIRS))
        for record, pair, expected in zip(s4["baseline"], PAIRS, self.baseline_metrics):
            self.assertEqual((record["first"], record["second"]), pair)
            self.assert_metric_map(record["metrics"], expected)
        primary = self.baseline_metrics[0]
        self.assertEqual(s4["printed_primary_pass"], not lost_criteria(primary, primary)["combined"])

    def assert_s4_witness(self, witness, include_other, budget=None, pair=None):
        state = tuple(tuple(row) for row in witness)
        self.assertEqual(len(state), len(GROUPS))
        changes = []
        for i, (before, after) in enumerate(zip(self.baseline, state)):
            self.assertEqual(sum(before), sum(after))
            dt, do = before[0] - after[0], before[2] - after[2]
            self.assertGreaterEqual(dt, 0)
            self.assertGreaterEqual(do, 0)
            self.assertLessEqual(dt, before[0])
            self.assertLessEqual(do, before[2])
            self.assertEqual(after[1], before[1] + dt + do)
            if not include_other:
                self.assertEqual(do, 0)
            if pair is not None and GROUPS[i] not in pair:
                self.assertEqual(after, before)
            changes.append(dt + do)
        if budget is not None:
            self.assertLessEqual(sum(changes), budget)
        self.assertIn(state, self.oracle_metrics)
        return state, sum(changes)

    def assert_attained_ranges(self, actual, cases, pair_index, include_other, budget=None):
        expected_metrics = self.baseline_metrics[pair_index]
        self.assertEqual(set(actual), set(expected_metrics))
        for metric, bounds in actual.items():
            values = [self.oracle_metrics[state][pair_index][metric] for state, _, _ in cases]
            self.assert_number(bounds["min"], min(values))
            self.assert_number(bounds["max"], max(values))
            for side in ("min", "max"):
                witness, _ = self.assert_s4_witness(bounds[side + "_witness"], include_other, budget)
                self.assert_number(bounds[side], self.oracle_metrics[witness][pair_index][metric])

    def test_every_envelope_extremum_and_witness_uses_a_coherent_shared_control(self):
        for name, include_other in (("toxicity_only", False), ("toxicity_and_other", True)):
            envelope = self.report["s4"]["envelopes"][name]
            cases = self.oracle_states[include_other]
            eligible = sum(row[0] + (row[2] if include_other else 0) for row in self.baseline)
            self.assertEqual(envelope["state_count"], len(cases))
            self.assertEqual(envelope["eligible_labels"], eligible)
            self.assertEqual(envelope["total_n"], sum(map(sum, self.baseline)))
            self.assertEqual(len({state for state, _, _ in cases}), len(cases))
            for i, baseline_record in enumerate(self.report["s4"]["baseline"]):
                pair_id = baseline_record["id"]
                with self.subTest(envelope=name, pair=pair_id):
                    self.assert_attained_ranges(envelope["ranges"][pair_id], cases, i, include_other)
                    survives = all(not lost_criteria(self.oracle_metrics[state][i],
                                                     self.baseline_metrics[i])["combined"]
                                   for state, _, _ in cases)
                    self.assertEqual(envelope["verdicts"][pair_id], "survives" if survives else "fails")
        primary_id = self.report["s4"]["baseline"][0]["id"]
        self.assertEqual(self.report["s4"]["verdict"],
                         self.report["s4"]["envelopes"]["toxicity_and_other"]["verdicts"][primary_id])

    def assert_breakpoint(self, actual, cases, pair_index, criterion, include_other, pair=None):
        eligible_cases = cases
        if pair is not None:
            eligible_cases = [case for case in cases
                              if all(case[0][i] == self.baseline[i]
                                     for i, group in enumerate(GROUPS) if group not in pair)]
        failing = [(state, k) for state, k, _ in eligible_cases
                   if lost_criteria(self.oracle_metrics[state][pair_index],
                                    self.baseline_metrics[pair_index])[criterion]]
        expected = min((k for _, k in failing), default=None)
        if expected is None:
            self.assertIsNone(actual)
            return
        self.assertIsNotNone(actual)
        self.assertEqual(actual["k"], expected)
        scope = pair if pair is not None else GROUPS
        total = sum(sum(self.counts[group]) for group in scope)
        eligible = sum(self.counts[group][0] + (self.counts[group][2] if include_other else 0)
                       for group in scope)
        self.assertEqual(actual["fraction_animals"], F(expected, total))
        self.assertEqual(actual["fraction_eligible"], F(expected, eligible) if eligible else F(0))
        witness, changed = self.assert_s4_witness(actual["witness"], include_other, expected, pair)
        self.assertEqual(changed, expected)
        self.assertTrue(lost_criteria(self.oracle_metrics[witness][pair_index],
                                      self.baseline_metrics[pair_index])[criterion])

    def test_every_pair_local_and_global_breakpoint_is_minimal(self):
        criteria = set(lost_criteria(self.baseline_metrics[0], self.baseline_metrics[0]))
        for name, include_other in (("toxicity_only", False), ("toxicity_and_other", True)):
            envelope = self.report["s4"]["envelopes"][name]
            cases = self.oracle_states[include_other]
            for i, (record, pair) in enumerate(zip(self.report["s4"]["baseline"], PAIRS)):
                points = envelope["pair_local_breakpoints"][record["id"]]
                self.assertEqual(set(points), criteria)
                for criterion, actual in points.items():
                    with self.subTest(envelope=name, pair=pair, criterion=criterion):
                        self.assert_breakpoint(actual, cases, i, criterion, include_other, pair)
            points = envelope["primary_global_breakpoints"]
            self.assertEqual(set(points), criteria)
            for criterion, actual in points.items():
                self.assert_breakpoint(actual, cases, 0, criterion, include_other)

    def test_all_cumulative_budgets_recompute_every_pair_metric(self):
        for name, include_other in (("toxicity_only", False), ("toxicity_and_other", True)):
            envelope = self.report["s4"]["envelopes"][name]
            eligible = sum(row[0] + (row[2] if include_other else 0) for row in self.baseline)
            self.assertEqual([b["k"] for b in envelope["budgets"]], list(range(eligible + 1)))
            for budget in envelope["budgets"]:
                k = budget["k"]
                cases = [case for case in self.oracle_states[include_other] if case[1] <= k]
                self.assertEqual(budget["states"], len(cases))
                self.assertEqual(budget["fraction_animals"], F(k, sum(map(sum, self.baseline))))
                self.assertEqual(budget["fraction_eligible"], F(k, eligible))
                for i, record in enumerate(self.report["s4"]["baseline"]):
                    with self.subTest(envelope=name, k=k, pair=record["id"]):
                        self.assert_attained_ranges(budget["pairs"][record["id"]], cases,
                                                    i, include_other, k)

    def test_single_drug_completions_and_all_conditional_bounds_are_exhaustive(self):
        for drug in ("cap", "gem"):
            adaptive, mtd = self.counts[drug + "_adaptive"], self.counts[drug + "_mtd"]
            actual = self.report["singles"][drug]
            completions = tuple(tuple(value - (j == i) for j, value in enumerate(adaptive))
                                for i in range(len(CAUSES)) if adaptive[i])
            self.assertEqual(actual["printed_adaptive"], adaptive)
            self.assertEqual(actual["printed_mtd"], mtd)
            self.assertEqual(actual["completions"], completions)
            self.assertEqual(actual["unconditional_verdict"], "untestable")
            self.assertTrue(all(sum(row) == 12 for row in completions))
            for model, include_other in (("as_printed_completion", None),
                                         ("toxicity_only", False), ("toxicity_and_other", True)):
                cases = []
                for completion in completions:
                    left = ((completion, 0),) if include_other is None else independent_row_states(completion, include_other)
                    right = ((mtd, 0),) if include_other is None else independent_row_states(mtd, include_other)
                    for a, b in product(left, right):
                        first, second = a[0], b[0]
                        comp, tox = independent_tests(first, second)
                        metrics = {f"{cause}_difference": F(first[i], sum(first)) - F(second[i], sum(second))
                                   for i, cause in enumerate(CAUSES)}
                        metrics.update(composition_p=comp, toxicity_p=tox,
                                       composition_bonferroni=min(F(1), 4 * comp),
                                       toxicity_bonferroni=min(F(1), 4 * tox))
                        cases.append((completion, (first, second), metrics))
                bounds = actual["bounds"][model]
                self.assertEqual(bounds["completion_scenarios"], len(cases))
                self.assertEqual(bounds["unique_states"], len({state for _, state, _ in cases}))
                self.assertEqual(set(bounds["ranges"]), set(cases[0][2]))
                indexed = {(completion, state): metrics for completion, state, metrics in cases}
                for metric, observed in bounds["ranges"].items():
                    with self.subTest(drug=drug, model=model, metric=metric):
                        values = [metrics[metric] for _, _, metrics in cases]
                        self.assertEqual(observed["min"], min(values))
                        self.assertEqual(observed["max"], max(values))
                        for side in ("min", "max"):
                            witness = observed[side + "_witness"]
                            key = (tuple(witness["completion"]), tuple(tuple(row) for row in witness["state"]))
                            self.assertIn(key, indexed)
                            self.assertEqual(observed[side], indexed[key][metric])

    def test_generated_markdown_and_all_hashes_match_the_checked_report(self):
        readme = (NOTE / "README.md").read_text(encoding="utf-8")
        # The authenticated connector recreated publication commits. Retain
        # enough metadata to verify the original local pre-inference commit
        # without requiring Git, credentials, or network access in this test.
        commit_payload = (readme.split("```text\n", 1)[1].split("\n```", 1)[0]
                          + "\n").encode("utf-8")
        commit_object = (b"commit " + str(len(commit_payload)).encode("ascii")
                         + b"\0" + commit_payload)
        self.assertEqual(hashlib.sha1(commit_object).hexdigest(), self.report["protocol_commit"])
        self.assertEqual(readme.count(ANALYSIS.BEGIN), 1)
        self.assertEqual(readme.count(ANALYSIS.END), 1)
        rendered = ANALYSIS.render_report(self.report)
        actual = readme.split(ANALYSIS.BEGIN, 1)[1].split(ANALYSIS.END, 1)[0]
        self.assertEqual(actual.strip(), rendered.strip())
        expected = dict(NORMALIZED_SHA256)
        expected[Path("analysis.py").as_posix()] = hashlib.sha256(normalized_bytes(NOTE / "analysis.py")).hexdigest()
        self.assertEqual(self.report["hashes"], expected)
        for key, value in self.report["hashes"].items():
            self.assertEqual(key, Path(key).as_posix())
            self.assertEqual(key, PureWindowsPath(key).as_posix())
            self.assertEqual(ANALYSIS.normalized_sha256(NOTE / key), value)

    def test_lf_crlf_and_cr_checkouts_give_identical_hashes_and_results(self):
        paths = [NOTE / "PROTOCOL.md", NOTE / "analysis.py", *sorted((NOTE / "data").glob("*.csv"))]
        with tempfile.TemporaryDirectory(prefix="terminal-causes-newlines-") as directory:
            temporary = Path(directory)
            (temporary / "data").mkdir()
            for newline in (b"\n", b"\r\n", b"\r"):
                for source in paths:
                    target = temporary / source.relative_to(NOTE)
                    target.write_bytes(normalized_bytes(source).replace(b"\n", newline))
                    self.assertEqual(ANALYSIS.normalized_sha256(source), ANALYSIS.normalized_sha256(target))
                with patch("socket.socket", side_effect=AssertionError("Network access forbidden")):
                    alternate = ANALYSIS.build_report(temporary)
                self.assertEqual(alternate, self.report)

    def test_normalized_hashes_and_portable_path_keys(self):
        for key, expected in NORMALIZED_SHA256.items():
            self.assertEqual(key, Path(key).as_posix())
            self.assertEqual(key, PureWindowsPath(key).as_posix())
            self.assertNotIn("\\", key)
            self.assertEqual(hashlib.sha256(normalized_bytes(NOTE / key)).hexdigest(), expected)

    def test_stdlib_offline_scope_and_runtime(self):
        for path in (NOTE / "analysis.py", Path(__file__)):
            imports = set()
            for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
                if isinstance(node, ast.Import):
                    imports.update(alias.name.split(".")[0] for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imports.add(node.module.split(".")[0])
            self.assertTrue(imports <= sys.stdlib_module_names | {"__future__"})
            self.assertFalse(imports & {"socket", "urllib", "http", "requests", "subprocess"})
        self.assertLess(time.monotonic() - self.started, 60)


if __name__ == "__main__":
    unittest.main()
