"""Nonclinical aggregate-evidence reproducibility; no network or patient data.

The tests check transcription snapshots and calculations, not the truth of
papers or biological validity. All empirical inputs come from the CSVs.
Small mathematical edge fixtures below are explicitly synthetic unit tests.
"""

from __future__ import annotations

import ast
import copy
import csv
from decimal import Decimal
import hashlib
import importlib.util
import io
import math
from pathlib import Path
import sys
import tempfile
import time
import unittest


REPO = Path(__file__).resolve().parents[1]
EVIDENCE = REPO / "research" / "containment-evidence"
SPEC = importlib.util.spec_from_file_location("containment_evidence", EVIDENCE / "analysis.py")
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)

# Populated after the source audit. Logical newline normalization is part of
# this certificate: git checkout's LF/CRLF policy cannot change these hashes.
NORMALIZED_SHA256 = {
    # HASHES_BEGIN
    'PROTOCOL.md': 'fbe2c9bc37c58a90a4e08a17eaceccd3cb3576814b9b375316b14661fe524890',
    'data/studies.csv': '43aabec151b9995d11a5bf864eb257873ef3633db64487866d43ab36d19b3f80',
    'data/aggregates.csv': '0e3476f71fc696fbb980d708c3cfafc300d472898ffc89c2a717e32ee9d4ff2a',
    'data/comparisons.csv': '65aec53ce60b3221835a986d3b00ab9e0a1480fc1a47f18200cb313a74a7a57a',
    'data/mapping.csv': '31b6b1bc73d9b8eaef087a9694d49bb5bad73860cf0bf9f216fe03fb57e86cea',
    'data/gates.csv': '447499ea9db9d09e3affeacb05a56f16d88c5f8accf3ee3c26618f110dc5ede4',
    'data/searches.csv': '01e3116aeb94b29a79963f463cfb6f167ab907ddee5e2edd6a01f9b8540998f9',
    # HASHES_END
}


def independent_tables():
    tables = {}
    for name in ANALYSIS.TABLES:
        with (EVIDENCE / "data" / f"{name}.csv").open(encoding="utf-8", newline="") as stream:
            tables[name] = list(csv.DictReader(stream))
    return tables


class ContainmentEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tables = independent_tables()
        cls.result = ANALYSIS.compute(cls.tables)
        cls.aggregates = {r["row_id"]: r for r in cls.tables["aggregates"]}
        cls.studies = {s["study_id"]: s for s in cls.tables["studies"]}
        cls.contrasts = {r["comparison_id"]: r for r in cls.result["contrasts"]}
        cls.verdicts = {r["prediction"]: r for r in cls.result["predictions"]}

    def test_frozen_protocol_and_data_hashes_use_normalized_text(self):
        self.assertEqual(set(NORMALIZED_SHA256), {"PROTOCOL.md", *(
            f"data/{name}.csv" for name in ANALYSIS.TABLES)})
        for relative, expected in NORMALIZED_SHA256.items():
            with self.subTest(file=relative):
                data = (EVIDENCE / relative).read_bytes()
                normalized = data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
                self.assertEqual(hashlib.sha256(normalized).hexdigest(), expected)
                self.assertEqual(ANALYSIS.sha256(EVIDENCE / relative), expected)

    def test_windows_and_linux_read_identically(self):
        with tempfile.TemporaryDirectory(prefix="evidence-newlines-") as directory:
            temporary = Path(directory)
            for newline in (b"\r\n", b"\r", b"\n"):
                for name in ANALYSIS.TABLES:
                    source = EVIDENCE / "data" / f"{name}.csv"
                    logical = source.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
                    target = temporary / f"{name}.csv"
                    target.write_bytes(logical.replace(b"\n", newline))
                    self.assertEqual(ANALYSIS.sha256(source), ANALYSIS.sha256(target))
                self.assertEqual(ANALYSIS.read_tables(temporary), self.tables)
                self.assertEqual(ANALYSIS.compute(ANALYSIS.read_tables(temporary)), self.result)

    def test_provenance_rights_and_public_aggregate_schema(self):
        for name, rows in self.tables.items():
            for row in rows:
                with self.subTest(table=name, row=row.get("row_id", row.get("study_id"))):
                    for field in ("source_url", "locator", "access_date", "reuse_basis"):
                        self.assertTrue(row[field].strip())
                    self.assertTrue(row["source_url"].startswith("https://"))
                    self.assertEqual(row["access_date"], "2026-11-17")
                    self.assertFalse(set(row) & {"patient_id", "subject_id", "sample_id", "mouse_id",
                                                 "date_of_birth", "individual_trajectory", "genome"})
        for row in self.aggregates.values():
            study = self.studies[row["study_id"]]
            self.assertIn(study["disposition"], {"included", "included_model", "included_single_arm"})
            self.assertEqual(study["license_url"], "https://creativecommons.org/licenses/by/4.0/")
            self.assertIn("CC BY 4.0", row["reuse_basis"])
        excluded_ids = {s["study_id"] for s in self.studies.values()
                        if s["disposition"] not in {"included", "included_model", "included_single_arm"}}
        self.assertFalse(excluded_ids & {r["study_id"] for r in self.aggregates.values()})

    def test_requested_sources_and_negative_candidates_are_visible(self):
        required = {"zhang2017", "zhang2022", "gatenby_2009", "enriquez_navas_2016", "smalley_2019",
                    "anzadapt", "strobl2021", "strobl2022", "west2020", "brady2020", "brady2022_rbat",
                    "pressley_2026", "algazi2020", "interim2024", "gonzalez_cao2021", "bfr14_2024",
                    "coin2011", "improve2025", "improve2026"}
        self.assertTrue(required <= set(self.studies))
        self.assertEqual(self.studies["anzadapt"]["disposition"], "no_posted_results_observed")
        self.assertIn("not inspectable", self.studies["adaptive_capecitabine"]["limitations"].lower())

    def test_every_contrast_recomputed_independently_with_decimal(self):
        for row in self.tables["comparisons"]:
            actual = self.contrasts[row["comparison_id"]]
            with self.subTest(contrast=row["comparison_id"]):
                if row["experimental_row"]:
                    x, y = (self.aggregates[row[key]] for key in ("experimental_row", "reference_row"))
                    self.assertEqual(x["unit"], y["unit"])
                    self.assertEqual(x["statistic"], y["statistic"])
                    self.assertEqual(actual["experimental_interval"], ANALYSIS.interval(x))
                    self.assertEqual(actual["reference_interval"], ANALYSIS.interval(y))
                    if x["relation"] in {"exact", "at_least"} and y["relation"] == "exact" and x["value"] and y["value"]:
                        left, right = Decimal(x["value"]), Decimal(y["value"])
                        self.assertAlmostEqual(actual["difference"], float(left - right), places=11)
                        self.assertAlmostEqual(actual["ratio"], float(left / right), places=11)
                        self.assertEqual(actual["bound_type"], ">=" if x["relation"] == "at_least" else "")
                    else:
                        self.assertIsNone(actual["difference"])
                        self.assertIsNone(actual["ratio"])
                else:
                    self.assertIsNone(actual["difference"])
                    self.assertIsNone(actual["ratio"])
                if row["hr_row"]:
                    hr = self.aggregates[row["hr_row"]]
                    self.assertEqual(actual["hr"], float(Decimal(hr["value"])))
                    self.assertEqual(actual["hr_low"], float(Decimal(hr["lower"])))
                    self.assertEqual(actual["hr_high"], float(Decimal(hr["upper"])))
                    excludes = Decimal(hr["upper"]) < 1 or Decimal(hr["lower"]) > 1
                    self.assertEqual(actual["interval_excludes_null"], excludes)
                else:
                    self.assertIsNone(actual["hr"])

    def test_all_reported_hr_contrasts_retained(self):
        source_hrs = {r["row_id"] for r in self.aggregates.values() if r["statistic"] == "hazard_ratio"}
        displayed_hrs = {c["hr_row"] for c in self.tables["comparisons"] if c["hr_row"]}
        self.assertEqual(source_hrs, displayed_hrs)
        self.assertEqual(len(source_hrs), 19)

    def test_every_input_and_reported_table_is_reproduced(self):
        block = ANALYSIS.render(self.tables, self.result)
        ANALYSIS.update_readme(block, check=True)
        for row_id in self.aggregates:
            self.assertIn(row_id, block)
        for comparison_id in self.contrasts:
            self.assertIn(comparison_id, block)
        for study_id in self.studies:
            self.assertIn(study_id, block)
        self.assertEqual(self.result["counts"]["aggregate_rows"], len(self.aggregates))
        self.assertEqual(self.result["counts"]["study_records"], len(self.studies))
        self.assertEqual(self.result["counts"]["search_records"], len(self.tables["searches"]))

    def test_censored_medians_and_limits_are_not_point_estimates(self):
        for key in ("int24_os_c", "adt22_rttp", "adt22_psattp", "hockings_sensitive_at", "hockings_mixed_at"):
            self.assertEqual(self.aggregates[key]["relation"], "not_reached")
            self.assertEqual(self.aggregates[key]["value"], "")
        for key in ("interim_os", "hockings_sensitive", "hockings_mixed"):
            self.assertIsNone(self.contrasts[key]["ratio"])
        self.assertEqual(self.aggregates["z17_ttp_a"]["relation"], "at_least")
        self.assertEqual(Decimal(self.aggregates["z17_ttp_c"]["value"]), Decimal("14"))
        self.assertEqual(self.contrasts["zhang2017_ttp"]["bound_type"], ">=")
        bound = self.aggregates["pressley_omnibus_p"]
        self.assertEqual(bound["relation"], "less_than")
        self.assertEqual(bound["lower"], "")
        self.assertEqual(bound["upper"], "")

    def test_prostate_endpoints_exposure_and_overlap(self):
        self.assertEqual(self.aggregates["z22_os_hr"]["measure"], "OS_HR")
        self.assertIsNone(self.contrasts["zhang2022_ttp"]["hr"])
        self.assertEqual(Decimal(self.aggregates["z22_dose_rate"]["value"])
                         + Decimal(self.aggregates["z22_time_off"]["value"]), Decimal(1))
        self.assertIn("Discussion", self.aggregates["z22_dose_rate"]["notes"])
        self.assertEqual(self.studies["zhang2017"]["cohort_id"], self.studies["zhang2022"]["cohort_id"])
        self.assertEqual(self.contrasts["zhang2017_ttp"]["role"], "superseded_same_cohort")
        self.assertEqual(self.studies["interim2024"]["cohort_id"], "ISRCTN18183156")

    def test_confidence_levels_iqr_and_negative_evidence(self):
        self.assertEqual(self.aggregates["coin11_os_hr"]["uncertainty_type"], "80% CI")
        for key in ("coin11_os_a", "coin11_os_c"):
            self.assertIn("IQR", self.aggregates[key]["uncertainty_type"])
        for key in ("gonzalez_cao_pfs", "interim_pfs", "coin_os", "pressley_adaptive_continuous"):
            self.assertLess(self.contrasts[key]["difference"], 0)
        self.assertEqual(self.contrasts["pressley_adaptive_continuous"]["difference"], -94)
        self.assertLess(self.contrasts["interim_pfs"]["hr_low"], 1)
        self.assertGreater(self.contrasts["interim_pfs"]["hr_high"], 1)

    def test_shared_controls_neutral_arms_and_toxic_comparator_remain_visible(self):
        rows = [r for r in self.tables["comparisons"] if r["study_id"] == "seyedi_2024"]
        self.assertEqual({r["independence_group"] for r in rows}, {"seyedi_shared_experiment"})
        for key in ("seyedi_pp_dm_mtd", "seyedi_pp_int_mtd", "seyedi_tandem_dm_mtd", "seyedi_tandem_int_mtd"):
            self.assertIn("toxicity", self.contrasts[key]["notes"].lower())
        neutral = self.aggregates["seyedi_gem_survival_missing"]
        self.assertEqual(neutral["value"], "")
        self.assertIn("no significant", neutral["notes"].lower())
        self.assertEqual(self.contrasts["hockings_resistant_st_vehicle"]["role"], "contextual_or_secondary")
        self.assertEqual(self.contrasts["pressley_adaptive_control"]["role"], "contextual_or_secondary")
        self.assertEqual(self.aggregates["hockings_n_adaptive_range"]["value"], "")
        self.assertEqual(self.aggregates["hockings_n_adaptive_range"]["lower"], "3")
        self.assertEqual(self.aggregates["hockings_n_adaptive_range"]["upper"], "5")

    def test_fitted_counterfactuals_are_not_independent_clinical_outcomes(self):
        for c in self.result["contrasts"]:
            if c["study_id"] == "brady2020":
                self.assertEqual(c["role"], "model_counterfactual")
                self.assertEqual(c["design"], "fitted_model")
                self.assertNotIn(c["comparison_id"], self.verdicts["P6a"]["evidence"])
        for row in self.aggregates.values():
            if row["statistic"].startswith("source_inter_individual"):
                self.assertFalse(row["lower"] or row["upper"])
        self.assertEqual(self.aggregates["bacevic_monolayer_fitness"]["statistic"], "reported_range")

    def test_all_parameter_statuses_and_conditional_transformations(self):
        for study_id in self.studies:
            rows = [r for r in self.tables["mapping"] if r["study_id"] == study_id]
            self.assertEqual(sorted(r["parameter"] for r in rows), sorted(("d", "c", "J", "A", "lambda")))
            for row in rows:
                if row["parameter"] in {"c", "J", "A", "lambda"}:
                    self.assertEqual(row["status"], "unidentifiable")
        transformed = self.result["conditional_parameters"]
        self.assertEqual({r["study_id"] for r in transformed}, {"zhang2022", "wang_2021_reachability"})
        for row in transformed:
            product = Decimal(1)
            for key in row["input_ids"].split(";"):
                product *= Decimal(self.aggregates[key]["value"])
            self.assertAlmostEqual(row["value"], float(Decimal(1) / product), places=12)
            self.assertEqual(row["status"], "conditional")
            self.assertIn("not identifiable", row["interval"])

    def test_mechanistic_eligibility_is_computed_not_favourable_outcome_inferred(self):
        for gate in self.tables["gates"]:
            self.assertNotEqual(gate["complete_parameters"], "yes")
            self.assertNotEqual(gate["measured_supply"], "yes")
        for key in ("P1", "P2", "P2a", "P3", "P4", "P5"):
            self.assertEqual(self.verdicts[key]["verdict"], "untestable")
            self.assertEqual(self.verdicts[key]["eligible"], 0)
        self.assertEqual(self.result["formal_mechanistic_tests"], 0)
        self.assertEqual(self.result["pooled_effects"], 0)
        changed = copy.deepcopy(self.tables)
        for key in ("supply_series", "mean_burden_endpoint", "model_assumptions"):
            changed["gates"][0][key] = "yes"
        with self.assertRaisesRegex(ValueError, "New eligible P1"):
            ANALYSIS.compute(changed)

    def test_original_broad_rule_and_post_extraction_sensitivities_both_reported(self):
        # Independently apply the source-interval rule; no normal-approximation
        # CI is manufactured from COIN's 80% CI or a p-value.
        adverse = []
        for row in self.tables["comparisons"]:
            if row["role"] == "observed_latest" and row["design"] == "randomized" and row["hr_row"]:
                source = self.aggregates[row["hr_row"]]
                if Decimal(source["lower"]) > 1:
                    adverse.append(row["comparison_id"])
        self.assertEqual(adverse, ["coin_os"])
        self.assertEqual(self.verdicts["P6b"]["evidence"], ";".join(adverse))
        self.assertEqual(self.verdicts["P6b"]["verdict"], "falsified")
        self.assertIn("80%", self.verdicts["P6b"]["reason"])
        self.assertIn("NOT a 95%", self.verdicts["P6b"]["reason"])
        for key in ("P6b-control", "P6b-control95"):
            self.assertEqual(self.verdicts[key]["verdict"], "untestable")
            self.assertIn("POST-EXTRACTION", self.verdicts[key]["reason"])
        self.assertEqual(self.verdicts["P6a"]["verdict"], "consistent")

    def test_missing_and_invalid_measurements_fail_closed(self):
        bad = copy.deepcopy(self.tables)
        bad["aggregates"][0]["source_url"] = ""
        with self.assertRaises(ValueError):
            ANALYSIS.compute(bad)
        bad = copy.deepcopy(self.tables)
        bad["aggregates"][0]["upper"] = "1"
        bad["aggregates"][0]["lower"] = ""
        with self.assertRaises(ValueError):
            ANALYSIS.compute(bad)
        for args in ((1, 1, 2, 1), (1, 2, 2, 1), (1, 2, 3, 0)):
            with self.assertRaises(ValueError):
                ANALYSIS.threshold(*args)

    def test_synthetic_formula_edges_do_not_generate_empirical_estimates(self):
        # Synthetic algebra fixtures only, never reported as biological data.
        self.assertAlmostEqual(ANALYSIS.threshold(1, 2, 3, 4), math.log(2) / 4)
        # An independent positive-series quadrature tests the saturated-c
        # counterexample in an admissible SYNTHETIC parameter regime. This
        # is not a tumour calibration or an empirical threshold estimate.
        def establishment_integral(alpha):
            panels = 256
            values = []
            for i in range(panels):
                density = (i + 0.5) / panels
                term = total = 1.0
                for k in range(1, 121):
                    term *= alpha * density / (0.2 + k * 0.03)
                    total += term
                values.append(0.2 / (0.8 * total))
            return math.fsum(values) / (panels * 0.03)
        # rS=.3, delta=.03, L=16, B=1.2, aSR=beta=.02, g0=.2.
        self.assertLess((0.03 + 0.02 * 1.2) * 16, 1)
        self.assertGreater(0.3, 0.02 * 1.2)
        self.assertGreater(16, 1 / 0.2)
        self.assertGreater(0.2, 0.02 * 1.2)
        alphas = (0.24, 0.36)
        self.assertTrue(all(0.2 - a <= 0 for a in alphas))
        low, high = (establishment_integral(a) for a in alphas)
        self.assertLess(high, low)
        d, c, j, supply = 5, 16, 16 * (1 + 0.03 / 0.3), 0.2
        advantage_low = c - d - (j - d) * math.exp(-supply * low)
        advantage_high = c - d - (j - d) * math.exp(-supply * high)
        self.assertLess(advantage_high, advantage_low)
        self.assertGreater(ANALYSIS.threshold(d, c, j, high), ANALYSIS.threshold(d, c, j, low))
        # No aggregate lineage counts qualify, so these are API boundary
        # checks rather than an analysis of clinical nonprogression counts.
        with self.assertRaises(ValueError):
            ANALYSIS.wilson(2, 1)
        low, high = ANALYSIS.wilson(0, 10)
        self.assertAlmostEqual(low, 0)
        self.assertGreater(high, 0)
        self.assertTrue(math.isinf(ANALYSIS.supply_from_zero_fraction(0, 10)["theta"]))

    def test_standard_library_only_and_bounded_execution(self):
        tree = ast.parse((EVIDENCE / "analysis.py").read_text(encoding="utf-8"))
        imports = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.add(node.module.split(".")[0])
        self.assertTrue(imports <= sys.stdlib_module_names | {"__future__"})
        self.assertFalse(imports & {"socket", "urllib", "http", "subprocess"})
        start = time.monotonic()
        tables = ANALYSIS.read_tables()
        result = ANALYSIS.compute(tables)
        ANALYSIS.update_readme(ANALYSIS.render(tables, result), check=True)
        self.assertLess(time.monotonic() - start, 60)


if __name__ == "__main__":
    unittest.main()
