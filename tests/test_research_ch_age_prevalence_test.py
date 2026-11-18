"""Nonclinical CH aggregate audit: reproducibility, mathematics and leakage checks.

Standard library only; no network, external executables, downloads or patient
records. Passing tests is not source verification or biological validation.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import importlib.util
import json
import math
from collections import defaultdict
from datetime import date
from pathlib import Path, PurePosixPath, PureWindowsPath
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
NOTE = REPO / "research" / "ch-age-prevalence-test"
DATA = NOTE / "data"
THIS_TEST = Path(__file__).resolve()


def _load(name: str):
    # Unique module names avoid consuming/changing another research test's
    # short 'model', 'analysis', or 'verify' module names during discovery.
    path = NOTE / (name + ".py")
    spec = importlib.util.spec_from_file_location("test_ch_age_prevalence_" + name, path)
    if spec is None or spec.loader is None:
        raise AssertionError(f"cannot load {path.name}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ANALYSIS = _load("analysis")
VERIFY = _load("verify")


def _read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _normalized(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def _hash(path: Path) -> str:
    return hashlib.sha256(_normalized(path)).hexdigest()


def _data_key(filename: str) -> str:
    return Path(filename).as_posix()


class CHAgePrevalenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Exactly one full report recomputation checks every saved numerical
        # leaf and the complete generated README section. Reuse its result.
        cls.report = ANALYSIS.check_or_write()
        cls.mathematics = VERIFY.verify_math()
        cls.tables = {path.relative_to(DATA).as_posix(): _read_csv(path)
                      for path in sorted(DATA.glob("*.csv"))}
        cls.observations = cls.tables[_data_key("prevalence.csv")]
        cls.context = cls.tables[_data_key("context_aggregates.csv")]
        cls.cohorts = {row["cohort_id"]: row for row in cls.tables[_data_key("cohorts.csv")]}
        cls.parameters = ANALYSIS.model.load_parameters()
        # LOCO must use unrounded intensities: serialized very small numbers
        # are intentionally rounded. Production model caches are already warm.
        cls.raw_scored = ANALYSIS.score_rows(cls.observations, cls.parameters, cls.cohorts)
        cls.external = [row for row in cls.raw_scored if row["role"] == "external_age"]

    def assertNumericClose(self, actual, expected):
        if expected is None:
            self.assertIsNone(actual)
        else:
            self.assertTrue(math.isclose(actual, expected,
                                        rel_tol=ANALYSIS.REPORT_FLOAT_REL_TOL,
                                        abs_tol=ANALYSIS.REPORT_FLOAT_ABS_TOL),
                            (actual, expected))

    def test_complete_report_and_readme_recompute(self):
        # setUpClass has already run the default checking entry point once.
        saved = json.loads((NOTE / "report.json").read_text(encoding="utf-8"))
        ANALYSIS.assert_reports_close(saved, self.report)
        ANALYSIS.assert_reports_close(ANALYSIS.numeric_clean(self.raw_scored),
                                      self.report["scored_rows"])
        text = (NOTE / "README.md").read_text(encoding="utf-8")
        self.assertEqual(text.count(ANALYSIS.START), 1)
        self.assertEqual(text.count(ANALYSIS.END), 1)
        actual = (ANALYSIS.START + text.split(ANALYSIS.START, 1)[1]
                  .split(ANALYSIS.END, 1)[0] + ANALYSIS.END)
        self.assertEqual(actual, ANALYSIS.render_results(self.report))
        self.assertIn("nonclinical", text.lower())
        self.assertIn("no medical advice", text.lower())
        expected_hashes = {path.relative_to(NOTE).as_posix(): _hash(path)
                           for path in sorted(DATA.glob("*.csv"))}
        self.assertEqual(self.report["normalized_input_sha256"], expected_hashes)
        self.assertEqual(self.report["frozen_protocol_sha256"], _hash(NOTE / "PROTOCOL.md"))

    def test_report_comparison_does_not_ignore_structural_changes(self):
        # Exercise the comparator's safeguards without repeating a model fit.
        expected = {"float": 1.0, "count": 2, "label": "retained", "rows": [1, 2]}
        for altered in (
            {"float": 1.01, "count": 2, "label": "retained", "rows": [1, 2]},
            {"float": 1.0, "count": 3, "label": "retained", "rows": [1, 2]},
            {"float": 1.0, "count": 2, "label": "changed", "rows": [1, 2]},
            {"float": 1.0, "count": 2, "label": "retained", "rows": [2, 1]},
            {"float": 1.0, "count": 2, "label": "retained"},
        ):
            with self.assertRaises(AssertionError):
                ANALYSIS.assert_reports_close(altered, expected)

    def test_manifest_covers_every_payload_with_normalized_posix_hashes(self):
        manifest = json.loads((NOTE / "manifest.json").read_text(encoding="utf-8"))
        files = manifest["files"]
        expected_paths = {path.relative_to(REPO).as_posix(): path
                          for path in NOTE.rglob("*")
                          if path.is_file() and path.name != "manifest.json"
                          and "__pycache__" not in path.parts and path.suffix != ".pyc"}
        expected_paths[THIS_TEST.relative_to(REPO).as_posix()] = THIS_TEST
        self.assertEqual(set(files), set(expected_paths))
        for key, path in expected_paths.items():
            with self.subTest(path=key):
                pure = PurePosixPath(key)
                self.assertFalse(pure.is_absolute())
                self.assertNotIn("..", pure.parts)
                self.assertNotIn("\\", key)
                self.assertEqual(key, Path(*pure.parts).as_posix())
                self.assertEqual(key, PureWindowsPath(key).as_posix())
                self.assertEqual(files[key], _hash(path))

    def test_crlf_and_cr_normalization_is_platform_independent(self):
        with tempfile.TemporaryDirectory(dir=NOTE) as temporary:
            directory = Path(temporary)
            lf, mixed = directory / "lf.txt", directory / "mixed.txt"
            lf.write_bytes(b"alpha\nbeta\ngamma\n")
            mixed.write_bytes(b"alpha\r\nbeta\rgamma\r\n")
            self.assertEqual(_normalized(lf), _normalized(mixed))
            self.assertEqual(ANALYSIS.normalized_bytes(lf), ANALYSIS.normalized_bytes(mixed))
            self.assertEqual(_hash(lf), _hash(mixed))
            self.assertEqual(ANALYSIS.digest(lf), ANALYSIS.digest(mixed))
        windows = PureWindowsPath("research", "ch-age-prevalence-test", "data", "prevalence.csv")
        native = Path("research", "ch-age-prevalence-test", "data", "prevalence.csv")
        self.assertEqual(windows.as_posix(), native.as_posix())

    def test_protocol_commit_payload_and_frozen_files(self):
        evidence = json.loads((NOTE / "freeze_evidence.json").read_text(encoding="utf-8"))
        payload = evidence["local_commit_payload"].encode("utf-8")
        git_object = b"commit " + str(len(payload)).encode("ascii") + b"\0" + payload
        self.assertEqual(hashlib.sha1(git_object).hexdigest(), evidence["local_commit_sha"])
        self.assertIn(("parent " + evidence["base_main_sha"] + "\n").encode("ascii"), payload)
        self.assertRegex(evidence["published_protocol_commit_sha"], r"^[0-9a-f]{40}$")
        frozen_paths = [NOTE / "PROTOCOL.md", DATA / "model_parameters.csv", DATA / "mutation_supply.csv"]
        expected = {path.relative_to(REPO).as_posix(): _hash(path) for path in frozen_paths}
        self.assertEqual(evidence["normalized_sha256"], expected)
        # This verifies the recorded object identity and file integrity. It
        # cannot independently attest a self-authored timestamp or blinding.
        self.assertIn("not", evidence["chronology"].lower())
        self.assertIn("third-party", evidence["chronology"].lower())

    def test_all_csv_rows_have_provenance_and_unique_keys(self):
        id_fields = {
            _data_key("cohorts.csv"): ("cohort_id",),
            _data_key("model_parameters.csv"): ("model_id",),
            _data_key("mutation_supply.csv"): ("gene",),
            _data_key("prevalence.csv"): ("record_id",),
            _data_key("context_aggregates.csv"): ("record_id",),
            _data_key("math_scenarios.csv"): ("fixture_id",),
            _data_key("prediction_grid.csv"): ("model_id", "age", "vaf"),
            _data_key("scenarios.csv"): ("scenario_id",),
            _data_key("source_audit.csv"): ("study_id",),
        }
        self.assertEqual(set(self.tables), set(id_fields))
        for name, rows in self.tables.items():
            with self.subTest(csv=name):
                self.assertTrue(rows)
                keys = [tuple(row[field] for field in id_fields[name]) for row in rows]
                self.assertEqual(len(set(keys)), len(keys))
                for row in rows:
                    for field in ("source_url", "locator", "access_date", "reuse_basis"):
                        self.assertTrue(row[field].strip(), (name, field))
                    date.fromisoformat(row["access_date"])
                    self.assertTrue(row["source_url"].startswith("https://"))
        scored_ids = {row["record_id"] for row in self.observations}
        context_ids = {row["record_id"] for row in self.context}
        self.assertFalse(scored_ids & context_ids)
        self.assertEqual(scored_ids, {row["record_id"] for row in self.report["scored_rows"]})
        for row in self.observations + self.context:
            self.assertIn(row["cohort_id"], self.cohorts)
            if row["cases"]:
                self.assertTrue(row["total"])
                self.assertLessEqual(int(row["cases"]), int(row["total"]))
            if row["noncases"]:
                self.assertEqual(int(row["cases"])+int(row["noncases"]), int(row["total"]))
        for row in self.observations:
            self.assertIn(row["model_id"], self.parameters)

    def test_frozen_supply_joins_and_source_denominators(self):
        supplies = {row["gene"]: float(row["mu_per_cell_year"])
                    for row in self.tables[_data_key("mutation_supply.csv")]}
        for name, params in self.parameters.items():
            self.assertEqual(params["mu_per_cell_year"], supplies[name])
        self.assertTrue(math.isclose(supplies["POOLED10"],
                                    math.fsum(value for key, value in supplies.items() if key != "POOLED10"),
                                    rel_tol=1e-12, abs_tol=1e-18))
        # Independent checksum of the complete printed CNSR-III age table,
        # not a vector of expected fitted predictions. Table 1 supplies both
        # column totals; the source URLs and locators are on each CSV row.
        cn = [row for row in self.observations if row["cohort_id"] == "qiu2023_cnsr"]
        self.assertEqual(sum(int(row["cases"]) for row in cn), 220)
        self.assertEqual(sum(int(row["total"]) for row in cn), 6016)
        zhang = [row for row in self.observations if row["cohort_id"] == "zhang2025_ukb"]
        self.assertTrue(zhang)
        for row in zhang:
            self.assertTrue(row["noncases"])
            self.assertEqual(int(row["total"]), int(row["cases"])+int(row["noncases"]))
            self.assertIn("sum", row["denominator_rule"].lower())

    def test_uncertainty_missing_counts_and_source_conflicts_are_retained(self):
        conflicts = []
        by_id = {row["record_id"]: row for row in self.raw_scored}
        for row in self.observations:
            observation = by_id[row["record_id"]]["observation"]
            self.assertLessEqual(observation["lower"], observation["point"])
            self.assertLessEqual(observation["point"], observation["upper"])
            if not row["cases"]:
                self.assertIsNone(observation["cases"])
                self.assertIsNone(by_id[row["record_id"]]["source_comparison"]["binomial_deviance"])
                if not row["ci_lower_percent"]:
                    self.assertIn("rounding", observation["uncertainty"])
                    self.assertIn("sampling uncertainty unavailable", observation["uncertainty"])
            if row["cases"] and row["prevalence_percent"]:
                point = int(row["cases"])/int(row["total"])
                alternative = float(row["prevalence_percent"])/100
                radius = 0.5*10**(-int(row["percentage_decimals"]))/100
                conflict = abs(point-alternative) > radius+1e-15
                self.assertEqual(observation["source_conflict"], conflict)
                if conflict:
                    conflicts.append(row["record_id"])
                    self.assertLessEqual(observation["lower"], max(0, alternative-radius))
                    self.assertGreaterEqual(observation["upper"], min(1, alternative+radius))
                    self.assertIn("hull", observation["uncertainty"])
                    self.assertFalse(by_id[row["record_id"]]["strict_eligible"])
        self.assertTrue(conflicts, "the audited printed-count conflicts must remain visible")
        reported_conflicts = [row["record_id"] for row in self.report["scored_rows"]
                              if row["observation"]["source_conflict"]]
        self.assertEqual(conflicts, reported_conflicts)

    def test_equal_parent_scores_and_complete_sensitivity_coverage(self):
        summary = self.report["external_substitute_summary"]
        groups = defaultdict(list)
        for row in self.external:
            prediction = row["source_prediction"]["uniform_age_mean"]
            if prediction is not None:
                groups[row["parent_cohort"]].append(100*(row["observation"]["point"]-prediction))
        self.assertEqual(summary["rows"], len(self.external))
        self.assertEqual(summary["parents"], len({row["parent_cohort"] for row in self.external}))
        self.assertEqual(summary["point_rows"], sum(map(len, groups.values())))
        self.assertEqual(summary["point_parents"], len(groups))
        if groups:
            mae = sum(sum(map(abs, errors))/len(errors) for errors in groups.values())/len(groups)
            mse = sum(sum(value*value for value in errors)/len(errors) for errors in groups.values())/len(groups)
            self.assertNumericClose(summary["macro_mae_pp"], mae)
            self.assertNumericClose(summary["macro_rmse_pp"], math.sqrt(mse))
        self.assertEqual(summary["nonoverlap_rows"],
                         sum(not row["source_comparison"]["interval_compatible"] for row in self.external))
        self.assertEqual({row["scenario_id"] for row in self.report["sensitivities"]},
                         {row["scenario_id"] for row in self.tables[_data_key("scenarios.csv")]})
        grid_keys = {(row["model_id"], float(row["age"]), float(row["vaf"]))
                     for row in self.tables[_data_key("prediction_grid.csv")]}
        self.assertEqual(grid_keys, {(row["model_id"], row["age"], row["vaf"])
                                     for row in self.report["future_predictions"]})
        expected_detection = {key for key, row in self.cohorts.items()
                              if row["depth_for_illustration"] and row["min_alt_reads"]}
        self.assertEqual(expected_detection, {row["cohort_id"] for row in self.report["detection_illustrations"]})

    def test_whole_parent_holdouts_and_heldout_outcome_perturbations(self):
        row_parents = {row["record_id"]: row["parent_cohort"] for row in self.raw_scored}
        finite_parents = sorted({row["parent_cohort"] for row in self.external
                                 if row["source_prediction"]["midpoint_intensity"] is not None})
        self.assertGreaterEqual(len(finite_parents), 2)
        for parent in finite_parents:
            with self.subTest(parent=parent):
                before = ANALYSIS.loco(self.external, parent)
                self.assertIn("multiplier", before)
                self.assertNotIn(parent, before["training_parents"])
                self.assertTrue(all(row_parents[key] != parent for key in before["training_ids"]))
                altered = copy.deepcopy(self.external)
                for row in altered:
                    if row["parent_cohort"] == parent:
                        point = row["observation"]["point"]
                        row["observation"]["point"] = point/2 if point else 0.25
                after = ANALYSIS.loco(altered, parent)
                for key in ("multiplier", "constant_baseline", "training_ids", "training_parents"):
                    self.assertEqual(before[key], after[key])
                self.assertEqual([(r["record_id"], r["prediction"], r["constant_prediction"])
                                  for r in before["predictions"]],
                                 [(r["record_id"], r["prediction"], r["constant_prediction"])
                                  for r in after["predictions"]])
        # Fixed-model deletion changes aggregate scores only; re-score the
        # retained rows without changing any other source or model input.
        parent = finite_parents[0]
        retained_source = [row for row in self.observations
                           if self.cohorts[row["cohort_id"]]["parent_cohort"] != parent]
        retained = ANALYSIS.score_rows(retained_source, self.parameters, self.cohorts)
        original = {row["record_id"]: row for row in self.raw_scored}
        for row in retained:
            self.assertEqual(row["source_prediction"], original[row["record_id"]]["source_prediction"])

    def test_loco_finite_parent_guards_and_training_exclusions(self):
        finite = [row for row in self.external
                  if row["source_prediction"]["midpoint_intensity"] is not None]
        parent = finite[0]["parent_cohort"]
        one_parent = [row for row in finite if row["parent_cohort"] == parent]
        open_only = copy.deepcopy(finite[0])
        open_only["parent_cohort"] = "synthetic_open_only_parent"
        open_only["record_id"] = "synthetic_open_only_row"
        open_only["source_prediction"]["midpoint_intensity"] = None
        insufficient = ANALYSIS.loco_collection(one_parent + [open_only])
        self.assertEqual(insufficient["finite_parents"], 1)
        self.assertEqual(insufficient["results"], [])
        self.assertIn("untestable", insufficient["status"])
        no_test_age = ANALYSIS.loco(self.external + [open_only], open_only["parent_cohort"])
        self.assertIn("held parent has no finite-band rows", no_test_age["status"])
        training = next(row for row in finite if row["parent_cohort"] != parent)
        for endpoint in (0.0, 1.0):
            altered = copy.deepcopy(self.external)
            next(row for row in altered if row["record_id"] == training["record_id"])["observation"]["point"] = endpoint
            result = ANALYSIS.loco(altered, parent)
            self.assertNotIn(training["record_id"], result["training_ids"])
            self.assertTrue(any(row["record_id"] == training["record_id"]
                                and "no pseudocount" in row["reason"]
                                for row in result["training_exclusions"]))

    def test_strict_eligibility_age_assay_and_parent_independence_gates(self):
        statuses = defaultdict(set)
        for cohort in self.cohorts.values():
            statuses[cohort["parent_cohort"]].add(cohort["independence"])
        source_by_id = {row["record_id"]: row for row in self.observations}
        for row in self.raw_scored:
            source = source_by_id[row["record_id"]]
            cohort = self.cohorts[row["cohort_id"]]
            expected = (source["analysis_role"] in ("external_age", "external_gene_age")
                        and bool(source["age_lower"]) and bool(source["age_upper"])
                        and source["endpoint_match"] == "yes"
                        and cohort["assay_calibrated"] == "yes"
                        and statuses[cohort["parent_cohort"]] == {"independent"}
                        and not row["observation"]["source_conflict"])
            self.assertEqual(row["strict_eligible"], expected)
        strict = [row for row in self.raw_scored if row["strict_eligible"]]
        parents = {row["parent_cohort"] for row in strict}
        self.assertEqual(self.report["strict_eligible_rows"], len(strict))
        self.assertEqual(self.report["strict_eligible_parents"], len(parents))
        gene_age = any(row["role"] == "external_gene_age" and row["model_id"] != "POOLED10"
                       for row in strict)
        if len(parents) < 2 or not gene_age:
            self.assertEqual(self.report["strict_verdict"], "untestable")
        # Positive and negative controls exercise the guards even when all
        # real endpoints are unmatched. Modified copies are never saved.
        finite = next(copy.deepcopy(row) for row in self.observations
                      if row["analysis_role"] == "external_age" and row["age_lower"] and row["age_upper"]
                      and not ANALYSIS.observed(row)["source_conflict"])
        finite["endpoint_match"] = "yes"
        cohorts = copy.deepcopy(self.cohorts)
        linked_parent = cohorts[finite["cohort_id"]]["parent_cohort"]
        for cohort in cohorts.values():
            if cohort["parent_cohort"] == linked_parent:
                cohort["independence"] = "independent"
                cohort["assay_calibrated"] = "yes"
        matched = ANALYSIS.score_rows([finite], self.parameters, cohorts)[0]
        self.assertTrue(matched["strict_eligible"])
        open_band = dict(finite, age_upper="")
        self.assertFalse(ANALYSIS.score_rows([open_band], self.parameters, cohorts)[0]["strict_eligible"])
        linked = copy.deepcopy(cohorts[finite["cohort_id"]])
        linked["cohort_id"] = "synthetic_linked_release"
        linked["independence"] = "overlap_uncertain"
        cohorts[linked["cohort_id"]] = linked
        contaminated = ANALYSIS.score_rows([finite], self.parameters, cohorts)[0]
        self.assertFalse(contaminated["strict_eligible"])
        self.assertEqual(contaminated["parent_independence_status"], "not established independent")

    def test_endpoint_separation_and_common_supply_intersections(self):
        expected_groups = defaultdict(list)
        for row in self.raw_scored:
            key = (row["cohort_id"], row["model_id"], row["endpoint_id"], row["analysis_role"])
            expected_groups[key].append(row)
        stored = self.report["supply_feasibility"]
        self.assertEqual(len(stored), len(expected_groups))
        for result in stored:
            key = tuple(result[field] for field in ("cohort_id", "model_id", "endpoint_id", "analysis_role"))
            rows = expected_groups[key]
            low = max(row["supply_feasible"][0] for row in rows)
            high = min(row["supply_feasible"][1] for row in rows)
            self.assertEqual(result["row_ids"], [row["record_id"] for row in rows])
            self.assertEqual(result["lower"], ANALYSIS.numeric_clean(low))
            self.assertEqual(result["upper"], ANALYSIS.numeric_clean(high))
            self.assertEqual(result["compatible"], low <= high)
            self.assertEqual(result["includes_frozen_supply"], low <= 1 <= high)
        independent = {row["parent_cohort"] for row in self.external
                       if all(cohort["independence"] == "independent"
                              for cohort in self.cohorts.values()
                              if cohort["parent_cohort"] == row["parent_cohort"])}
        self.assertEqual(set(self.report["external_independent_summary"]["parent_ids"]), independent)

    def test_gene_pool_event_inclusion_audit_recomputes_from_frozen_grid(self):
        # The event inclusion inequality does not require independent genes.
        # Only same-age/same-threshold candidate forecasts are comparable, and
        # source mutation-supply membership establishes the intended nesting.
        grid = self.tables[_data_key("prediction_grid.csv")]
        supplies = {row["gene"]: row for row in self.tables[_data_key("mutation_supply.csv")]}
        panel = set(supplies)-{"POOLED10"}
        pooled_grid = {(float(row["age"]), float(row["vaf"])): row
                       for row in grid if row["model_id"] == "POOLED10"}
        comparisons = 0
        violations = []
        for row in grid:
            gene = row["model_id"]
            if gene == "POOLED10":
                continue
            self.assertIn(gene, panel)
            self.assertEqual(self.parameters[gene]["mu_per_cell_year"],
                             float(supplies[gene]["mu_per_cell_year"]))
            age, vaf = float(row["age"]), float(row["vaf"])
            if (age, vaf) not in pooled_grid:
                continue
            comparisons += 1
            gene_probability = ANALYSIS.model.prevalence(self.parameters[gene], age, vaf)
            pooled_probability = ANALYSIS.model.prevalence(self.parameters["POOLED10"], age, vaf)
            if gene_probability > pooled_probability:
                violations.append(dict(gene=gene, age=age, vaf=vaf,
                                       gene_prediction=gene_probability,
                                       pooled_prediction=pooled_probability,
                                       gap_pp=100*(gene_probability-pooled_probability)))
        stored = self.report["subset_coherence"]
        self.assertEqual(stored["comparisons"], comparisons)
        self.assertEqual(stored["violation_count"], len(violations))
        ANALYSIS.assert_reports_close(stored["violations"], ANALYSIS.numeric_clean(violations))
        maximum = max((row["gap_pp"] for row in violations), default=0.0)
        self.assertNumericClose(stored["maximum_gap_pp"], maximum)
        self.assertIn("separate", stored["interpretation"])
        self.assertIn("joint", stored["interpretation"])

    def test_independent_mathematical_verifier(self):
        self.assertEqual(self.mathematics["status"], "passed")
        fixtures = self.tables[_data_key("math_scenarios.csv")]
        self.assertEqual(self.mathematics["fixture_rows"], len(fixtures))
        self.assertEqual(self.mathematics["ctmc_scenarios"], sum(row["kind"] == "bd_ctmc" for row in fixtures))
        self.assertEqual(self.mathematics["series_bound_scenarios"], sum(row["kind"] == "bd_series" for row in fixtures))
        self.assertEqual(self.mathematics["parameter_models"], len(self.parameters))
        path = NOTE / "math_verification.json"
        if path.exists():
            ANALYSIS.assert_reports_close(json.loads(path.read_text(encoding="utf-8")), self.mathematics)


if __name__ == "__main__":
    unittest.main()
