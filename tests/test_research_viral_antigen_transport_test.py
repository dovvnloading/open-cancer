"""Offline NONCLINICAL checks; source truth and reuse are separate manual audits.

All reference probabilities use Decimal exponentials rather than the production
math implementation. Synthetic cases exercise unavailable empirical branches;
they are mathematical fixtures, never new observations.
"""

import csv
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import itertools
from pathlib import Path, PureWindowsPath
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "viral-antigen-transport-test"
SPEC = importlib.util.spec_from_file_location("viral_transport", RESEARCH / "analysis.py")
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def read_rows(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def hash_reference(path):
    # Text universal-newline reading is independent of the byte replacement path.
    text = Path(path).read_text(encoding="utf-8")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def probability_reference(load, constants):
    if load == 0:
        return 0.0
    with localcontext() as context:
        context.prec = 55
        a, b, x = (Decimal(str(value)) for value in (constants["A"], constants["B"], load))
        return float(a * (-b / x).exp())


def half_unit_reference(text):
    if "." in text:
        return 0.5 * 10 ** (-len(text.split(".")[1]))
    return 0.5


def interval_reference(row, rounding=False):
    if not row["ci_low_percent"] or row["ci_level_percent"] != "95":
        return None
    lower, upper = float(row["ci_low_percent"]), float(row["ci_high_percent"])
    if rounding:
        lower = max(0, lower - half_unit_reference(row["ci_low_percent"]))
        upper = min(100, upper + half_unit_reference(row["ci_high_percent"]))
    return lower / 100, upper / 100


def range_reference(floor, ceiling, bounds):
    if bounds is None:
        return "no_published_95ci"
    if bounds[1] < floor:
        return "below_floor"
    if bounds[0] > ceiling:
        return "above_ceiling"
    if bounds[0] == ceiling:
        return "ceiling_boundary_only"
    return "compatible"


class ViralAntigenTransportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = ANALYSIS.analyze()
        cls.constants = cls.report["constants"]
        cls.outcomes = cls.report["outcomes"]
        cls.predictors = {r["outcome_id"]: r for r in cls.report["tmb"] if r["usable"] == "yes"}

    def test_immutable_protocol_and_inherited_files(self):
        self.assertEqual(ANALYSIS.PROTOCOL_SHA256,
                         "75c9e5a87ad2f3f132efd3f52491ccba5a142494c0f34f7c006ae047bd701ce9")
        for key, expected in ANALYSIS.PINNED.items():
            self.assertEqual(hash_reference(ROOT / key), expected, key)
        text = (RESEARCH / "PROTOCOL.md").read_text(encoding="utf-8")
        self.assertIn("e1c48e377bd438ac13ea94dc53d1a9e7db088de3", text)
        self.assertIn("No browser file downloads", text)

    def test_exact_carryover_not_a_target_fit(self):
        prior = ROOT / "research" / "merkel-viral-antigen-test" / "data"
        curves = {r["name"]: r for r in read_rows(prior / "frozen_curves.csv")}
        for name in ("A", "B"):
            self.assertEqual(self.constants[name], float(curves["reciprocal_holdout"][name]))
        for name in ("a", "b"):
            self.assertEqual(self.constants[name], float(curves["log_holdout"][name]))
        calibration = [r for r in read_rows(prior / "part_a_calibration.csv")
                       if r["calibration_set"] == "broad"]
        actual = [e["scenario"] for e in self.report["evaluations"]]
        self.assertEqual(len(actual), 15)
        for frozen, inherited in zip(actual, calibration):
            self.assertEqual(frozen["v"], float(inherited["v"]))
            self.assertEqual(frozen["c"], float(inherited["scale"]))
            self.assertEqual(frozen["target"], inherited["target"])
        # Changing a target observation can change its error, never the inputs.
        constants_before, scenarios_before = ANALYSIS.load_frozen()
        if self.predictors:
            chosen = next(r for r in self.outcomes if r["id"] in self.predictors)
            changed = dict(chosen, orr_percent="0")
            s = self.report["primary"]["scenario"]
            first = ANALYSIS.evaluate_stratum(chosen, self.predictors[chosen["id"]], s, self.constants)
            second = ANALYSIS.evaluate_stratum(changed, self.predictors[chosen["id"]], s, self.constants)
            for key in ("baseline", "offset", "log", "floor"):
                self.assertEqual(first[key], second[key])
        self.assertEqual(ANALYSIS.load_frozen(), (constants_before, scenarios_before))

    def test_every_row_has_auditable_provenance_and_missingness(self):
        for rows in (self.outcomes, self.report["tmb"], self.report["sources"]):
            self.assertEqual(len(rows), len({r["id"] for r in rows}))
            for row in rows:
                for key in ANALYSIS.PROVENANCE:
                    self.assertTrue(row[key].strip(), (row["id"], key))
                self.assertTrue(row["source_url"].startswith("https://"))
                self.assertTrue(row["terms_url"].startswith("https://"))
                self.assertEqual(row["nonclinical"], "yes")
                self.assertEqual(row["access_date"], "2026-11-18")
        for row in self.outcomes:
            self.assertEqual(bool(row["ci_low_percent"]), bool(row["ci_high_percent"]))
            if row["responders"] and row["n_orr"] and row["orr_percent"]:
                actual = 100 * int(row["responders"]) / int(row["n_orr"])
                consistent = abs(actual - float(row["orr_percent"])) <= half_unit_reference(row["orr_percent"]) + 1e-10
                self.assertTrue(consistent or row["count_rate_note"], row["id"])
        for row in self.report["tmb"]:
            if row["usable"] == "yes":
                self.assertEqual(row["unit"], "mutations/Mb")
                self.assertIn(row["population_match"], {"same", "within_cohort_subset"})
                self.assertIn(row["outcome_id"], {r["id"] for r in self.outcomes})

    def test_no_primary_duplicate_followups(self):
        ledger_families = {r["source_doi"]: r["family_id"] for r in self.report["sources"]}
        for row in self.outcomes + self.report["tmb"]:
            self.assertEqual(ledger_families[row["source_doi"]], row["family_id"])
        families = {}
        for row in self.outcomes:
            if row["role"] == "primary":
                families.setdefault(row["family_id"], set()).add(row["cohort_id"])
        for family, cohorts in families.items():
            self.assertEqual(len(cohorts), 1, (family, cohorts))

    def test_all_predictions_and_errors_independently(self):
        for evaluation in self.report["evaluations"]:
            s = evaluation["scenario"]
            for row in self.outcomes:
                value = evaluation["strata"][row["id"]]
                positive = row["viral_status"] == "positive"
                floor = probability_reference(s["v"], self.constants) if positive else 0.0
                self.assertAlmostEqual(value["floor"], floor, places=13)
                if (row["id"] not in self.predictors or row["viral_status"] == "unknown"
                        or row["eligible_monotherapy"] != "yes" or not row["orr_percent"]):
                    self.assertFalse(value["scorable"])
                    self.assertNotIn("baseline", value)
                    continue
                m = float(self.predictors[row["id"]]["median"])
                base = probability_reference(s["c"] * m, self.constants)
                offset = probability_reference(s["c"] * m + (s["v"] if positive else 0), self.constants)
                y = float(row["orr_percent"]) / 100
                self.assertAlmostEqual(value["baseline"], base, places=13)
                self.assertAlmostEqual(value["offset"], offset, places=13)
                with localcontext() as context:
                    context.prec = 50
                    x = Decimal(str(s["c"] * m))
                    log = 0.0 if x == 0 else float(Decimal(str(self.constants["a"])) +
                          Decimal(str(self.constants["b"])) * x.ln())
                self.assertAlmostEqual(value["log"], max(0, min(1, log)), places=13)
                self.assertAlmostEqual(value["baseline_residual"], y - base, places=13)
                self.assertAlmostEqual(value["offset_residual"], y - offset, places=13)
                self.assertAlmostEqual(value["baseline_se"], (y - base) ** 2, places=13)
                self.assertAlmostEqual(value["offset_se"], (y - offset) ** 2, places=13)
                self.assertAlmostEqual(value["delta"], (y - offset) ** 2 - (y - base) ** 2, places=13)
                if not positive:
                    self.assertEqual(value["baseline"], value["offset"])
                    self.assertEqual(value["delta"], 0)
                    self.assertEqual(value["rounding_delta_bounds"], (0, 0))
                base_envelope = (base / 2, (base + self.constants["A"]) / 2)
                offset_envelope = ((floor + offset) / 2, (offset + self.constants["A"]) / 2)
                for model, envelope in (("baseline", base_envelope), ("offset", offset_envelope)):
                    for actual, expected in zip(value[model + "_median_envelope"], envelope):
                        self.assertAlmostEqual(actual, expected, places=13)
                bounds = interval_reference(row)
                if bounds is not None:
                    for model, envelope in (("baseline", base_envelope), ("offset", offset_envelope)):
                        compatible = not (envelope[1] < bounds[0] or envelope[0] > bounds[1])
                        self.assertEqual(value[model + "_median_compatible"], compatible)
                    expanded = interval_reference(row, True)
                    dm = half_unit_reference(self.predictors[row["id"]]["median"])
                    pl = probability_reference(s["c"] * max(0, m - dm) + (s["v"] if positive else 0), self.constants)
                    pu = probability_reference(s["c"] * (m + dm) + (s["v"] if positive else 0), self.constants)
                    self.assertEqual(value["rounding_can_include"], not (pu < expanded[0] or pl > expanded[1]))

    def test_all_range_checks_independently(self):
        for evaluation in self.report["evaluations"]:
            s = evaluation["scenario"]
            for row in self.outcomes:
                v = evaluation["strata"][row["id"]]
                floor = probability_reference(s["v"], self.constants) if row["viral_status"] == "positive" else 0
                self.assertEqual(v["baseline_range"], range_reference(0, self.constants["A"], interval_reference(row)))
                self.assertEqual(v["offset_range"], range_reference(floor, self.constants["A"], interval_reference(row)))
                self.assertEqual(v["rounded_offset_range"], range_reference(floor, self.constants["A"], interval_reference(row, True)))
                expanded = interval_reference(row, True)
                if expanded is None:
                    self.assertIsNone(v["rounded_range_gap"])
                else:
                    gap = max(floor - expanded[1], expanded[0] - self.constants["A"], 0)
                    self.assertAlmostEqual(v["rounded_range_gap"], gap, places=13)

    def test_range_boundaries_and_rounding_are_not_confidence(self):
        self.assertEqual(ANALYSIS.range_status(.2, .5, (.05, .19)), "below_floor")
        self.assertEqual(ANALYSIS.range_status(.2, .5, (.51, .8)), "above_ceiling")
        self.assertEqual(ANALYSIS.range_status(.2, .5, (.5, .8)), "ceiling_boundary_only")
        self.assertEqual(ANALYSIS.range_status(.2, .5, (.1, .2)), "compatible")
        self.assertEqual(ANALYSIS.range_status(.2, .5, None), "no_published_95ci")
        for text in ("0", "100", "27.5", "5.06", "0.000"):
            self.assertEqual(ANALYSIS.half_unit(text), half_unit_reference(text))
        self.assertEqual(ANALYSIS.rounded_bounds("0", 100, 1), (0, .005))
        self.assertEqual(ANALYSIS.rounded_bounds("100", 100, 1), (.995, 1))

    def test_uncertainty_boxes_use_exact_affine_extrema(self):
        # Synthetic nonclinical probability fixtures, including reversed order.
        for p0, p1, bounds in ((.1, .4, (.2, .6)), (.6, .2, (.1, .5)), (.3, .3, (0, 1))):
            lower, upper = ANALYSIS.delta_box(p0, p1, bounds)
            endpoint_values = [(p1 - y) ** 2 - (p0 - y) ** 2 for y in bounds]
            self.assertEqual((lower, upper), (min(endpoint_values), max(endpoint_values)))
            for i in range(21):
                y = bounds[0] + (bounds[1] - bounds[0]) * i / 20
                delta = (p1 - y) ** 2 - (p0 - y) ** 2
                self.assertLessEqual(lower - 1e-14, delta)
                self.assertLessEqual(delta, upper + 1e-14)
        for evaluation in self.report["evaluations"]:
            for row in self.outcomes:
                v = evaluation["strata"][row["id"]]
                if "delta_bounds" in v:
                    endpoint_values = [(v["offset"] - y) ** 2 - (v["baseline"] - y) ** 2
                                       for y in interval_reference(row)]
                    self.assertEqual(tuple(v["delta_bounds"]), (min(endpoint_values), max(endpoint_values)))

    def test_rounding_bounds_enclose_correlated_predictions(self):
        for evaluation in self.report["evaluations"]:
            s = evaluation["scenario"]
            for row in self.outcomes:
                if row["id"] not in self.predictors:
                    continue
                r = self.predictors[row["id"]]
                m = float(r["median"])
                dm = half_unit_reference(r["median"])
                y = float(row["orr_percent"]) / 100
                dy = half_unit_reference(row["orr_percent"]) / 100
                lows, highs = evaluation["strata"][row["id"]]["rounding_delta_bounds"]
                for median, target in itertools.product((max(0, m - dm), m, m + dm),
                                                        (max(0, y - dy), y, min(1, y + dy))):
                    base = probability_reference(s["c"] * median, self.constants)
                    offset = probability_reference(s["c"] * median +
                        (s["v"] if row["viral_status"] == "positive" else 0), self.constants)
                    delta = (offset - target) ** 2 - (base - target) ** 2
                    self.assertLessEqual(lows - 1e-13, delta)
                    self.assertLessEqual(delta, highs + 1e-13)

    def test_median_envelope_respects_changed_lower_endpoint(self):
        c = self.constants
        s = self.report["primary"]["scenario"]
        floor = probability_reference(s["v"], c)
        p = probability_reference(s["v"] + 3, c)
        envelope = ANALYSIS.median_envelope(p, floor, c["A"])
        self.assertEqual(envelope, ((floor + p) / 2, (p + c["A"]) / 2))
        self.assertGreater(envelope[0], p / 2)
        # A symmetric three-atom distribution has median 3 and lies in the bound.
        mean = (probability_reference(s["v"], c) + p + probability_reference(s["v"] + 30, c)) / 3
        self.assertLessEqual(envelope[0], mean)
        self.assertLessEqual(mean, envelope[1])

    def test_pair_scoring_and_empty_panel_behavior(self):
        self.assertIsNone(ANALYSIS.balanced_summary([]))
        for evaluation in self.report["evaluations"]:
            for pair in evaluation["pairs"]:
                rows = [r for r in self.outcomes if r["pair_id"] == pair["pair_id"] and r["role"] == "primary"]
                if not pair["point_scorable"]:
                    self.assertNotIn("baseline", pair)
                    self.assertIsNone(pair["pass"])
                    continue
                self.assertEqual({r["viral_status"] for r in rows}, {"positive", "negative"})
                for model in ("baseline", "offset", "log"):
                    errors = [evaluation["strata"][r["id"]][model] - float(r["orr_percent"]) / 100 for r in rows]
                    mse = sum(e * e for e in errors) / 2
                    self.assertAlmostEqual(pair[model]["mse"], mse, places=14)
                    self.assertAlmostEqual(pair[model]["rmse_pp"], 100 * mse ** .5, places=12)
                    self.assertAlmostEqual(pair[model]["mae_pp"], 50 * sum(abs(e) for e in errors), places=12)
                if pair["full_interval_rule"]:
                    expected = pair["delta"] < 0 and all(
                        interval_reference(r)[0] <= evaluation["strata"][r["id"]]["offset"] <= interval_reference(r)[1]
                        for r in rows)
                    self.assertEqual(pair["pass"], expected)

    def test_verdict_recomputed_without_selecting_sensitivity(self):
        primary = self.report["primary"]
        self.assertEqual(primary["scenario"]["scale_label"], "1")
        self.assertEqual(primary["scenario"]["target"], "point")
        required = [r for r in self.outcomes if r["role"] == "primary" and r["cancer_group"] in ANALYSIS.REQUIRED
                    and r["eligible_monotherapy"] == "yes" and r["viral_status"] in {"positive", "negative"}]
        failures = [r["id"] for r in required if primary["strata"][r["id"]]["rounded_offset_range"]
                    in {"below_floor", "above_ceiling"}]
        paired_failures = [p["pair_id"] for p in primary["pairs"] if p["cancer_group"] in ANALYSIS.REQUIRED
                           and p["full_interval_rule"] and not p["pass"]]
        self.assertEqual(self.report["verdict"]["range_failures"], failures)
        self.assertEqual(self.report["verdict"]["pair_failures"], paired_failures)
        if failures or paired_failures:
            expected = "viral offset fails to transport"
        else:
            passed = {p["cancer_group"] for p in primary["pairs"] if p["full_interval_rule"] and p["pass"] and p["independent"]}
            incomplete = any(p["point_scorable"] and not p["full_interval_rule"]
                             for p in primary["pairs"] if p["cancer_group"] in ANALYSIS.REQUIRED)
            expected = "viral offset transports" if ANALYSIS.REQUIRED.issubset(passed) and not incomplete else "untestable"
        self.assertEqual(self.report["verdict"]["verdict"], expected)

    def test_secondary_context_or_overlap_cannot_manufacture_primary_success(self):
        # Synthetic records isolate the decision rule independently of the audit.
        row = {"id": "synthetic", "role": "secondary", "cancer_group": "hpv_hnscc", "historical_overlap": "possible",
               "eligible_monotherapy": "yes", "viral_status": "positive"}
        strata = {"synthetic": {"rounded_offset_range": "below_floor"}}
        self.assertEqual(ANALYSIS.verdict([row], strata, [])["verdict"], "untestable")
        row["role"] = "primary"
        self.assertEqual(ANALYSIS.verdict([row], strata, [])["verdict"], "viral offset fails to transport")
        self.assertEqual(ANALYSIS.verdict([row], strata, [], True)["verdict"], "untestable")
        row["eligible_monotherapy"] = "no"
        self.assertEqual(ANALYSIS.verdict([row], strata, [])["verdict"], "untestable")
        pairs = [{"pair_id": group, "cancer_group": group, "full_interval_rule": True,
                  "point_scorable": True, "pass": True, "independent": False} for group in ANALYSIS.REQUIRED]
        self.assertEqual(ANALYSIS.verdict([], {}, pairs)["verdict"], "untestable")
        for pair in pairs:
            pair["independent"] = True
        self.assertEqual(ANALYSIS.verdict([], {}, pairs)["verdict"], "viral offset transports")
        pairs.append({"pair_id": "synthetic_missing_interval", "cancer_group": "hpv_hnscc",
                      "point_scorable": True, "full_interval_rule": False, "pass": None, "independent": True})
        self.assertEqual(ANALYSIS.verdict([], {}, pairs)["verdict"], "untestable")

    def test_canonical_hashes_and_posix_keys(self):
        hashes = ANALYSIS.fingerprints()
        self.assertEqual(list(hashes), sorted(hashes))
        for key, value in hashes.items():
            self.assertNotIn("\\", key)
            self.assertEqual(Path(key).as_posix(), key)
            self.assertEqual(PureWindowsPath(key).as_posix(), key)
            self.assertEqual(value, hash_reference(ROOT / key))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "newline_fixture.txt"
            expected = hashlib.sha256("alpha\nbeta\nγ\n".encode()).hexdigest()
            for content in ("alpha\nbeta\nγ\n", "alpha\r\nbeta\r\nγ\r\n", "alpha\rbeta\rγ\r"):
                path.write_bytes(content.encode())
                self.assertEqual(ANALYSIS.sha256(path), expected)
                self.assertEqual(hash_reference(path), expected)

    def test_every_generated_readme_number_and_hash(self):
        actual = (RESEARCH / "README.md").read_text(encoding="utf-8")
        self.assertEqual(actual, ANALYSIS.report_text(self.report))
        self.assertEqual(actual.count(ANALYSIS.BEGIN), 1)
        self.assertEqual(actual.count(ANALYSIS.END), 1)
        self.assertIn("Nonclinical", actual)
        self.assertIn("f47119556e7d0ad6ceb546f11a3860863047b3d6", actual)

    def test_write_report_preserves_prose_and_is_idempotent(self):
        # Regression: opening the destination before reading the template
        # previously truncated the narrative. Exercise the actual CLI write path.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "README.md"
            before, after = "Nonclinical preface\n\n", "\n\nRetained limitations.\n"
            original = before + ANALYSIS.BEGIN + "\nold block\n" + ANALYSIS.END + after
            path.write_text(original, encoding="utf-8")
            expected = before + ANALYSIS.BEGIN + "\n\nSynthetic replacement\n\n" + ANALYSIS.END + after
            with patch.object(ANALYSIS, "HERE", root), \
                    patch.object(ANALYSIS, "analyze", return_value={}), \
                    patch.object(ANALYSIS, "render_report", return_value="Synthetic replacement\n"), \
                    patch("sys.argv", ["analysis.py", "--write-report"]):
                ANALYSIS.main()
                self.assertEqual(path.read_text(encoding="utf-8"), expected)
                ANALYSIS.main()
                self.assertEqual(path.read_text(encoding="utf-8"), expected)


if __name__ == "__main__":
    unittest.main()
