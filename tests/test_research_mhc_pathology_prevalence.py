"""Offline checks of the nonclinical MHC-I pathology-prevalence audit.

No medical or veterinary advice. Numerical consistency does not establish
source truth, population independence, immune causation, or biological validity.
"""

from contextlib import ExitStack
import copy
import csv
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path, PureWindowsPath
import tempfile
import time
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "mhc-pathology-prevalence"
FREEZE = "196f9970a80381bcf2cf0e2603ee33f843045e7c"
PROTOCOL_SHA256 = "e25aca9e7037e303739b2cdb0b2342fe2a6db1fd57f7d7a86156d3a767f7c06d"
N, K, SECONDARY = "records", "neoplasia", "malignant"
MASS, LIFE, MHC, DIET = (
    "log_body_mass_kg", "log_adult_life_expectancy_days", "log1p_mhc_i", "mammal_prey")
BEGIN, END = "<!-- BEGIN GENERATED RESULTS -->", "<!-- END GENERATED RESULTS -->"
DIAGNOSTICS = {
    "maximum_gradient_ratio", "maximum_absolute_gradient", "iterations",
    "maximum_fit_iterations", "line_search_halvings",
}
MAJOR = {"equal_species_training", "penalty_0_1", "penalty_10", "exclude_named_rodents"}
EXPECTED_SOURCES = {
    Path("vincze.csv").as_posix(): (
        "2c8a5c5b5fda6fb34c86b0e5a2b6a58dcc85cc3e",
        "8a937a2c8d4f1749ae33c1b12d2ed52e895852f9",
        "b9d6dedc9447e9f99b275cbfbdb54e3b0e318d08c12f0384572a65101d057744"),
    Path("compton.csv").as_posix(): (
        "4819e0933358c5fc7eab483a0c05161ea8c074b7",
        "d5e18150557a80b2dff1f53779182abd99e7d5ec",
        "70f8e2f6675e4c42dcf87acb0b826bb7c4ea32fed27a46900c1cfdf0bce70390"),
    Path("dictionary.txt").as_posix(): (
        "4819e0933358c5fc7eab483a0c05161ea8c074b7",
        "293705259ec5d0e4fabad000e183e162057e14c5",
        "39eec5bd11bd14b09ddbf6f3248322869711dd4968e853dafe0a238724a36140"),
}


def no_network():
    stack = ExitStack()
    for target in ("urllib.request.urlopen", "socket.create_connection", "socket.socket"):
        stack.enter_context(mock.patch(target, side_effect=AssertionError(
            "Network access is forbidden in the offline pathology test")))
    return stack


def load_module(name, path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


with no_network():
    A = load_module("mhc_pathology_analysis_test", RESEARCH / "analysis.py")
    F = load_module("mhc_pathology_fetch_test", RESEARCH / "data" / "fetch_verify.py")


def canonical_text(path):
    # Universal-newline decoding is independent of the explicit replacements
    # used by the analysis and preserves all other whitespace.
    with Path(path).open(encoding="utf-8", newline=None) as stream:
        return stream.read()


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def independent_rows():
    rows = read_csv(RESEARCH / "data" / "species_aggregates.csv")
    for row in rows:
        for field in (N, K, SECONDARY, DIET):
            row[field] = int(row[field])
        for field in (MASS, LIFE, MHC):
            row[field] = float(row[field])
    return rows


def normalize_name(name):
    return " ".join(name.replace("_", " ").split()).casefold()


def independent_probability(row, fit):
    features = [1.0, (row[MASS] - fit["means"][MASS]) / fit["scales"][MASS],
                (row[LIFE] - fit["means"][LIFE]) / fit["scales"][LIFE], float(row[DIET])]
    if fit["include_mhc"]:
        features.append((row[MHC] - fit["means"][MHC]) / fit["scales"][MHC])
    eta = math.fsum(x * b for x, b in zip(features, fit["coefficients"]))
    if eta < 0:
        exponential = math.exp(eta)
        probability = exponential / (1 + exponential)
    else:
        probability = 1 / (1 + math.exp(-eta))
    return features, probability


def independent_loss(k, n, probability):
    with localcontext() as context:
        context.prec = 50
        p = Decimal.from_float(min(1 - 1e-12, max(1e-12, probability)))
        fraction = Decimal(k) / Decimal(n)
        return float(-fraction * p.ln() - (1 - fraction) * (1 - p).ln())


class MhcPathologyPrevalenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.perf_counter()
        cls.network_guard = no_network()
        cls.addClassCleanup(cls.network_guard.close)
        cls.rows = independent_rows()
        cls.evaluations = []
        evaluate = A.evaluate

        def capture(rows, *args, **kwargs):
            result = evaluate(rows, *args, **kwargs)
            cls.evaluations.append((tuple(copy.deepcopy(rows)), result))
            return result

        # One full rebuild retains every sensitivity/deletion's fits before
        # build_report reduces those evaluations to report summaries.
        with mock.patch.object(A, "evaluate", side_effect=capture):
            cls.report = A.build_report()
        cls.rebuild_seconds = time.perf_counter() - cls.started

    def assert_close(self, actual, expected, label="", rel=2e-10, absolute=3e-12):
        self.assertTrue(math.isfinite(actual), label)
        self.assertTrue(math.isfinite(expected), label)
        self.assertTrue(math.isclose(actual, expected, rel_tol=rel, abs_tol=absolute),
                        label or f"{actual!r} != {expected!r}")

    def check_diagnostic(self, key, record, label):
        value = record[key]
        self.assertNotIsInstance(value, bool, label)
        self.assertIsInstance(value, (int, float), label)
        self.assertTrue(math.isfinite(value), label)
        self.assertGreaterEqual(value, 0, label)
        if key == "maximum_gradient_ratio":
            self.assertLessEqual(value, 1 + 1e-12, label)
        elif key == "maximum_absolute_gradient":
            self.assertGreater(record["gradient_tolerance"], 0, label)
            self.assertLessEqual(value, record["gradient_tolerance"] * (1 + 1e-12), label)
        elif key in ("iterations", "maximum_fit_iterations"):
            self.assertIsInstance(value, int, label)
            self.assertLessEqual(value, 100, label)
        elif key == "line_search_halvings":
            self.assertIsInstance(value, int, label)
            self.assertLessEqual(value, 44 * record["iterations"], label)

    def assert_tree(self, actual, expected, label="report", snapshot=False):
        if isinstance(expected, dict):
            self.assertIsInstance(actual, dict, label)
            self.assertEqual(set(actual), set(expected), label)
            for key in expected:
                if key in DIAGNOSTICS:
                    # A stopping certificate is not a portable exact value.
                    self.check_diagnostic(key, actual, label + "." + key)
                    self.check_diagnostic(key, expected, label + "." + key)
                else:
                    self.assert_tree(actual[key], expected[key], label + "." + key, snapshot)
        elif isinstance(expected, list):
            self.assertIsInstance(actual, list, label)
            self.assertEqual(len(actual), len(expected), label)
            for index, (left, right) in enumerate(zip(actual, expected)):
                self.assert_tree(left, right, f"{label}[{index}]", snapshot)
        elif isinstance(expected, float):
            self.assertNotIsInstance(actual, bool, label)
            self.assertIsInstance(actual, (float, int), label)
            tolerance = 3e-7 if snapshot else 3e-12
            self.assert_close(actual, expected, label, rel=3e-7 if snapshot else 2e-10,
                              absolute=tolerance)
        else:
            self.assertIs(type(actual), type(expected), label)
            self.assertEqual(actual, expected, label)

    def test_01_complete_saved_report_and_archived_markdown(self):
        def reject_nonfinite(token):
            raise AssertionError("Nonfinite JSON token: " + token)

        saved = json.loads(canonical_text(RESEARCH / "data" / "results.json"),
                           parse_constant=reject_nonfinite)
        self.assert_tree(saved, self.report, snapshot=True)
        self.assertEqual(A.compare_reports(self.report, saved), [])
        readme = canonical_text(RESEARCH / "README.md")
        self.assertEqual(readme.count(BEGIN), 1)
        self.assertEqual(readme.count(END), 1)
        block = readme[readme.index(BEGIN):readme.index(END) + len(END)] + "\n"
        self.assertEqual(block, A.render_results(saved))
        self.assertIn(FREEZE, readme)
        text = " ".join(readme.casefold().split())
        self.assertIn("nonclinical", text)
        self.assertIn("no medical or veterinary advice", text)
        self.assertIn("functional-load theorem", text)
        self.assertRegex(text, r"possible (?:individual-)?animal overlap")
        # Stored primary and secondary fits receive an independent stationarity
        # check as well as the tolerant fresh-versus-archived comparison.
        lookup = {row["species"]: row for row in self.rows}
        for key in ("primary", "secondary"):
            for fold in saved[key]["folds"]:
                train = [lookup[name] for name in fold["training_species_ids"]]
                self.check_fit(train, fold["baseline_fit"])
                self.check_fit(train, fold["extended_fit"])

    def test_02_recount_panel_reconciliation_and_endpoint_audit(self):
        rows = self.rows
        names = [r["species"] for r in rows]
        self.assertEqual(names, sorted(set(names)))
        self.assertEqual(A.load_rows(), rows)
        self.assertEqual(len(rows), 42)
        self.assertEqual({r["species"] for r in rows if r[K] == 0},
                         {"Pithecia_pithecia", "Saimiri_boliviensis"})
        self.assertIn("Heterocephalus_glaber", names)
        self.assertNotIn("Dolichotis_patagonum", names)
        for row in rows:
            self.assertGreater(row[N], 0)
            self.assertLessEqual(0, row[SECONDARY])
            self.assertLessEqual(row[SECONDARY], row[K])
            self.assertLessEqual(row[K], row[N])
            self.assertIn(row[DIET], (0, 1))
            self.assertGreaterEqual(row[MHC], 0)
            self.assert_close(math.expm1(row[MHC]), round(math.expm1(row[MHC])))
            self.assertFalse({"CMR", "ICM", "RDI", "PropMalignant"} & row.keys())
        ledger = read_csv(RESEARCH / "data" / "reconciliation.csv")
        audit = read_csv(RESEARCH / "data" / "endpoint_audit.csv")
        v_names, c_names, included, excluded, disagreements = set(), set(), [], [], []
        for item in ledger:
            key = item["normalized_species"]
            for field, keys in (("vincze_species", v_names), ("compton_species", c_names)):
                if item[field]:
                    self.assertEqual(normalize_name(item[field]), key)
                    self.assertNotIn(key, keys)
                    keys.add(key)
            if item["status"] == "included":
                self.assertEqual(item["reason"], "")
                included.append(item)
                if item["vincze_order"] != item["compton_order"]:
                    disagreements.append({"species": item["vincze_species"],
                                          "vincze_order": item["vincze_order"],
                                          "compton_order": item["compton_order"]})
            else:
                self.assertEqual(item["status"], "excluded")
                self.assertTrue(item["reason"])
                if item["vincze_species"] and item["compton_species"]:
                    excluded.append(item["vincze_species"])
        selected = {normalize_name(name) for name in names}
        self.assertEqual(selected, {r["normalized_species"] for r in included})
        self.assertEqual({r["normalized_species"] for r in ledger}, v_names | c_names)
        self.assertEqual({r["normalized_species"] for r in audit}, c_names)
        self.assertEqual(len(audit), len(c_names))
        mismatch = [0, 0, 0, 0]
        for item in audit:
            self.assertEqual(normalize_name(item["species"]), item["normalized_species"])
            values = [Decimal(item[k]) for k in (
                "neoplasia_residual", "malignant_residual", "malignancyknown_residual")]
            tolerances = [Decimal(item["neoplasia_tolerance"]),
                          Decimal(item["malignancy_tolerance"])]
            self.assertTrue(all(value.is_finite() for value in values + tolerances))
            self.assertTrue(all(value >= Decimal("2e-15") for value in tolerances))
            flags = [abs(values[0]) > tolerances[0], abs(values[1]) > tolerances[1],
                     abs(values[2]) > tolerances[1]]
            for index, flag in enumerate(flags):
                mismatch[index] += flag
            mismatch[3] += flags[2] and item["normalized_species"] in selected
        species_per_order = {g: sum(r["order"] == g for r in rows)
                             for g in sorted({r["order"] for r in rows})}
        counts = {
            "vincze_species": len(v_names), "compton_species": len(c_names),
            "normalized_matches": len(v_names & c_names), "union_species": len(ledger),
            "primary_species": len(rows), "orders": len(species_per_order),
            "species_per_order": species_per_order,
            "records": sum(r[N] for r in rows), "neoplasia": sum(r[K] for r in rows),
            "malignant": sum(r[SECONDARY] for r in rows),
            "zero_neoplasia_species": sum(r[K] == 0 for r in rows),
            "zero_malignant_species": sum(r[SECONDARY] == 0 for r in rows),
            "matched_exclusions": excluded,
            "primary_classes": sorted({r["compton_class"] for r in included}),
            "order_label_disagreements": disagreements,
            "endpoint_audit": dict(zip(("audited_compton_species",
                "neoplasia_prevalence_mismatches", "malignant_prevalence_mismatches",
                "malignancyknown_prevalence_mismatches",
                "panel_malignancyknown_prevalence_mismatches"), [len(audit), *mismatch])),
        }
        self.assertEqual(counts, self.report["prefit_summary"])
        self.assertEqual(counts, self.report["source_audit"]["prefit_summary"])
        self.assertEqual(counts, F.summarize(rows, ledger, audit))
        self.assertEqual((counts["records"], counts["neoplasia"], counts["malignant"]),
                         (2215, 312, 209))
        panel = self.report["panel"]
        for field in ("orders", "records", "neoplasia", "malignant", "species_per_order",
                      "zero_neoplasia_species", "zero_malignant_species"):
            self.assertEqual(panel[field], counts[field])
        self.assertEqual(panel["species"], len(rows))
        self.assertEqual(panel["named_rodents_present"], ["Heterocephalus_glaber"])
        self.assertEqual(panel["named_rodents_absent"], ["Dolichotis_patagonum"])
        by_name = {row["species"]: row for row in rows}
        for item in included:
            self.assertEqual(by_name[item["vincze_species"]]["order"], item["vincze_order"])

    def check_fit(self, train, fit):
        fields = [MASS, LIFE] + ([MHC] if fit["include_mhc"] else [])
        means = {field: math.fsum(row[field] for row in train) / len(train) for field in fields}
        scales = {field: math.sqrt(math.fsum((row[field] - means[field]) ** 2
                   for row in train) / len(train)) or 1.0 for field in fields}
        self.assert_tree(fit["means"], means)
        self.assert_tree(fit["scales"], scales)
        self.assertEqual(fit["feature_names"], ["intercept", MASS, LIFE, DIET] +
                         ([MHC] if fit["include_mhc"] else []))
        beta, ridge = fit["coefficients"], fit["penalty"]
        self.assertEqual(len(beta), len(fit["feature_names"]))
        self.assertGreater(ridge, 0)
        self.assertIn(fit["weighting"], ("counts", "species"))
        self.assertIn(fit["outcome"], (K, SECONDARY))
        mean_n = math.fsum(row[N] for row in train) / len(train)
        terms, likelihood, effective = [[] for _ in beta], [], []
        for row in train:
            features, p = independent_probability(row, fit)
            multiplier = mean_n / row[N] if fit["weighting"] == "species" else 1
            n, k = row[N] * multiplier, row[fit["outcome"]] * multiplier
            effective.append(n)
            self.assertGreater(p, 0)
            self.assertLess(p, 1)
            for index, feature in enumerate(features):
                terms[index].append((n * p - k) * feature)
            likelihood.append(-k * math.log(p) - (n - k) * math.log1p(-p))
        gradient = [math.fsum(values) + (ridge * beta[j] if j else 0)
                    for j, values in enumerate(terms)]
        total_n = math.fsum(effective)
        bound = max(1e-8, 1e-10 * total_n)
        self.assert_close(fit["gradient_tolerance"], bound)
        # The independent fsum/log formulation adds only rounding allowance;
        # the optimizer's own stored gradient must satisfy the declared bound.
        self.assertLessEqual(max(abs(g) for g in gradient), bound + 2e-11)
        self.check_diagnostic("maximum_absolute_gradient", fit, "fit gradient")
        for key in ("iterations", "line_search_halvings"):
            self.check_diagnostic(key, fit, key)
        objective = math.fsum(likelihood) + ridge * math.fsum(b * b for b in beta[1:]) / 2
        self.assert_close(fit["objective"], objective, absolute=2e-9)
        self.assert_close(fit["effective_training_denominator"], total_n)
        self.assertEqual(fit["training_records"], sum(r[N] for r in train))
        self.assertEqual(fit["training_species"], len(train))
        self.assertTrue(fit["converged"])

    def test_03_all_complete_refits_have_training_only_scaling_and_stationarity(self):
        total_fits = 0
        for rows, result in self.evaluations:
            lookup = {row["species"]: row for row in rows}
            seen = []
            self.assertEqual(len(result["folds"]), len({r["order"] for r in rows}))
            for fold in result["folds"]:
                train_ids, test_ids = fold["training_species_ids"], fold["test_species_ids"]
                self.assertFalse(set(train_ids) & set(test_ids))
                self.assertEqual(set(train_ids) | set(test_ids), set(lookup))
                train = [lookup[name] for name in train_ids]
                test = [lookup[name] for name in test_ids]
                self.assertEqual({r["order"] for r in test}, {fold["heldout_group"]})
                self.assertNotIn(fold["heldout_group"], {r["order"] for r in train})
                self.assertEqual(fold["training_orders"], sorted({r["order"] for r in train}))
                self.assertEqual(fold["test_orders"], [fold["heldout_group"]])
                seen.extend(test_ids)
                for key, mhc in (("baseline_fit", False), ("extended_fit", True)):
                    fit = fold[key]
                    self.assertEqual(fit["include_mhc"], mhc)
                    for field in ("penalty", "weighting", "outcome"):
                        self.assertEqual(fit[field], result[field])
                    self.check_fit(train, fit)
                    total_fits += 1
            self.assertEqual(sorted(seen), sorted(lookup))
        self.assertEqual(total_fits, self.report["numerical_audit"]["fits"])

    def test_04_every_prediction_and_all_weighted_scores(self):
        for rows, result in self.evaluations:
            lookup = {r["species"]: r for r in rows}
            folds = {f["heldout_group"]: f for f in result["folds"]}
            groups, weighted, clipping = {}, [[], []], 0
            self.assertEqual([p["species"] for p in result["predictions"]], sorted(lookup))
            for prediction in result["predictions"]:
                row = lookup[prediction["species"]]
                self.assertEqual(prediction["order"], row["order"])
                losses = []
                for index, model in enumerate(("baseline", "extended")):
                    _, p = independent_probability(row, folds[row["order"]][model + "_fit"])
                    loss = independent_loss(row[result["outcome"]], row[N], p)
                    self.assert_close(prediction[model + "_probability"], p)
                    self.assert_close(prediction[model + "_loss"], loss)
                    clipping += p < 1e-12 or p > 1 - 1e-12
                    losses.append(loss)
                    weighted[index].append(row[N] * loss)
                self.assertEqual(prediction["events"], row[result["outcome"]])
                self.assertEqual(prediction[N], row[N])
                self.assert_close(prediction["delta"], losses[0] - losses[1])
                groups.setdefault(row["order"], []).append((row, losses))
            group_means, improving = [[], []], 0
            for item in result["group_scores"]:
                members = groups[item["group"]]
                means = [math.fsum(losses[j] for _, losses in members) / len(members)
                         for j in range(2)]
                self.assertEqual(item["species"], len(members))
                self.assertEqual(item[N], sum(r[N] for r, _ in members))
                self.assertEqual(item["events"], sum(r[result["outcome"]] for r, _ in members))
                self.assertEqual(item["zero_case_species"],
                                 sum(r[result["outcome"]] == 0 for r, _ in members))
                for j, model in enumerate(("baseline", "extended")):
                    self.assert_close(item[model + "_loss"], means[j])
                    group_means[j].append(means[j])
                self.assert_close(item["delta"], means[0] - means[1])
                improving += means[0] > means[1]
            means = [math.fsum(values) / len(groups) for values in group_means]
            record_means = [math.fsum(values) / sum(r[N] for r in rows) for values in weighted]
            for index, model in enumerate(("baseline", "extended")):
                self.assert_close(result[model + "_loss"], means[index])
                self.assert_close(result["record_weighted_" + model + "_loss"], record_means[index])
            self.assert_close(result["delta"], means[0] - means[1])
            self.assert_close(result["record_weighted_delta"], record_means[0] - record_means[1])
            self.assertEqual(result["improving_groups"], improving)
            self.assertEqual(result["probabilities_clipped_for_scoring"], clipping)
            self.assertEqual(result["species"], len(rows))
            self.assertEqual(result["groups"], len(groups))
            self.assertEqual(result["orders"], len(groups))
            self.assertEqual(result[N], sum(r[N] for r in rows))
            self.assertEqual(result["events"], sum(r[result["outcome"]] for r in rows))
            self.assertEqual(result["zero_case_species"], sum(r[result["outcome"]] == 0 for r in rows))
            self.assertEqual(result["fit_count"], 2 * len(groups))

    def test_05_two_complete_refit_jackknives_and_t8_quantile(self):
        orders = sorted({row["order"] for row in self.rows})
        for name, key in (("primary", "uncertainty"), ("secondary", "secondary_uncertainty")):
            uncertainty = self.report[key]
            deletions = uncertainty["deletions"]
            self.assertEqual([item["deleted_order"] for item in deletions], orders)
            self.assertEqual(uncertainty["degrees_of_freedom"], 8)
            self.assertEqual(uncertainty["groups"], 9)
            self.assertEqual(uncertainty["outcome"], self.report[name]["outcome"])
            self.assertEqual(uncertainty["full_panel_point_estimate"], self.report[name]["delta"])
            with localcontext() as context:
                context.prec = 60
                full = Decimal.from_float(self.report[name]["delta"])
                pseudo = [9 * full - 8 * Decimal.from_float(item["delta"]) for item in deletions]
                center = sum(pseudo) / 9
                se = (sum((q - center) ** 2 for q in pseudo) / 72).sqrt()
                critical = Decimal("2.306004135204166")
                for item, expected in zip(deletions, pseudo):
                    self.assert_close(item["pseudovalue"], float(expected))
                    retained = tuple(r for r in self.rows if r["order"] != item["deleted_order"])
                    matches = [result for rows, result in self.evaluations
                               if rows == retained and result["outcome"] == item["outcome"]
                               and result["penalty"] == 1 and result["weighting"] == "counts"]
                    self.assertEqual(len(matches), 1)
                    actual_summary = {k: v for k, v in item.items()
                                      if k not in ("deleted_order", "pseudovalue")}
                    expected_summary = {k: v for k, v in matches[0].items()
                                        if k not in ("folds", "predictions")}
                    self.assert_tree(actual_summary, expected_summary)
                    self.assertEqual(item["orders"], 8)
                for field, expected in (("bias_corrected_center", center), ("standard_error", se),
                                        ("ci_lower", center - critical * se),
                                        ("ci_upper", center + critical * se)):
                    self.assert_close(uncertainty[field], float(expected))
            self.assert_close(uncertainty["delete_order_delta_min"], min(d["delta"] for d in deletions))
            self.assert_close(uncertainty["delete_order_delta_max"], max(d["delta"] for d in deletions))
            # Closed-form integral of the t8 density, independent of the stored
            # quantile constant and the jackknife code.
            t = uncertainty["t_critical"]
            s = t / math.sqrt(8 + t * t)
            cdf = 0.5 + 35 / 32 * (s - s ** 3 + 3 * s ** 5 / 5 - s ** 7 / 7)
            self.assert_close(cdf, 0.975)

    def test_06_sensitivities_and_secondary_preserve_declared_design(self):
        sensitivities = self.report["sensitivities"]
        self.assertEqual(set(sensitivities), MAJOR | {"exclude_recorded_zeros"})
        for name, summary in sensitivities.items():
            rows = self.rows
            penalty, weighting = 1, "counts"
            if name == "exclude_named_rodents":
                rows = [r for r in rows if r["species"] != "Heterocephalus_glaber"]
            elif name == "exclude_recorded_zeros":
                rows = [r for r in rows if r[K] > 0]
            elif name == "equal_species_training":
                weighting = "species"
            elif name == "penalty_0_1":
                penalty = 0.1
            elif name == "penalty_10":
                penalty = 10
            matches = [result for retained, result in self.evaluations
                       if retained == tuple(rows) and result["outcome"] == K
                       and result["weighting"] == weighting and result["penalty"] == penalty]
            self.assertEqual(len(matches), 1, name)
            self.assert_tree(summary, {k: v for k, v in matches[0].items()
                                      if k not in ("folds", "predictions")})
        self.assertEqual(sensitivities["exclude_recorded_zeros"]["zero_case_species"], 0)
        self.assertEqual({p["species"] for p in self.report["secondary"]["predictions"]},
                         {r["species"] for r in self.rows})
        self.assertEqual(self.report["secondary"]["outcome"], SECONDARY)
        audit = self.report["numerical_audit"]
        self.assertEqual(audit["complete_evaluations"], len(self.evaluations))
        self.assertEqual(audit["fits"], sum(e["fit_count"] for _, e in self.evaluations))
        self.assertTrue(audit["all_fits_converged"])
        self.assertEqual(audit["probabilities_clipped_for_scoring"],
                         sum(e["probabilities_clipped_for_scoring"] for _, e in self.evaluations))
        for key in ("maximum_fit_iterations", "maximum_gradient_ratio"):
            self.check_diagnostic(key, audit, key)

    def test_07_heldout_labels_and_covariates_cannot_enter_training(self):
        order = "Rodentia"
        before = next(f for f in self.report["primary"]["folds"] if f["heldout_group"] == order)
        labels = [dict(r, **{K: r[N] - r[K]}) if r["order"] == order else dict(r)
                  for r in self.rows]
        changed = A.evaluate(labels)
        after = next(f for f in changed["folds"] if f["heldout_group"] == order)
        self.assertEqual(before, after)
        old_predictions = {p["species"]: p for p in self.report["primary"]["predictions"]
                           if p["order"] == order}
        new_predictions = {p["species"]: p for p in changed["predictions"] if p["order"] == order}
        for species in old_predictions:
            for model in ("baseline", "extended"):
                self.assertEqual(old_predictions[species][model + "_probability"],
                                 new_predictions[species][model + "_probability"])
        self.assertTrue(any(old_predictions[s]["baseline_loss"] != new_predictions[s]["baseline_loss"]
                            for s in old_predictions))
        covariates = [dict(r, **{MASS: r[MASS] + 2, LIFE: r[LIFE] - 1,
                                MHC: r[MHC] + 1.5, DIET: 1 - r[DIET]})
                      if r["order"] == order else dict(r) for r in self.rows]
        changed = A.evaluate(covariates)
        after = next(f for f in changed["folds"] if f["heldout_group"] == order)
        self.assertEqual(before, after)
        self.assertTrue(any(p["baseline_probability"] != old_predictions[p["species"]]["baseline_probability"]
                            for p in changed["predictions"] if p["order"] == order))

    def test_08_numerical_boundaries_degenerate_fits_and_zero_loss(self):
        self.assert_close(A.pathology_loss(0, 100, 0.2), -math.log(0.8))
        self.assertGreater(A.pathology_loss(0, 100, 0.2), 0)
        for k, n, p in ((0, 10, 0), (0, 10, 1), (10, 10, 0), (10, 10, 1), (4, 19, 0.3)):
            self.assert_close(A.pathology_loss(k, n, p), independent_loss(k, n, p))
        for k, n, p in ((1, 0, 0.2), (-1, 10, 0.2), (11, 10, 0.2),
                        (2, 10, math.nan), (2, 10, -0.1), (2, 10, 1.1)):
            with self.assertRaises(ValueError):
                A.pathology_loss(k, n, p)
        for endpoint in (0, 1):
            constant = [dict(r, **{K: endpoint * r[N]}) for r in self.rows[:8]]
            with self.assertRaises(A.FitError):
                A.fit_model(constant, True)
        toy = [{MASS: 1.0, LIFE: 2.0, MHC: 3.0, DIET: i % 2, N: 10,
                K: 10 * (i % 2), SECONDARY: 0} for i in range(8)]
        fit = A.fit_model(toy, True)
        self.assertEqual(set(fit["scales"].values()), {1.0})
        self.check_fit(toy, fit)
        self.assertTrue(all(math.isfinite(value) for value in fit["coefficients"]))
        with mock.patch.object(A, "MAX_ITERATIONS", 0), self.assertRaises(A.FitError):
            A.fit_model(toy, True)

    def test_09_verdict_boundaries_and_secondary_cannot_rescue_primary(self):
        primary = {"orders": 9, "delta": 0.1, "improving_groups": 5,
                   "all_fits_converged": True,
                   "group_scores": [{"group": str(i)} for i in range(9)]}
        uncertainty = {"ci_lower": 0.01, "ci_upper": 0.2,
                       "deletions": [{"deleted_order": str(i), "all_fits_converged": True}
                                     for i in range(9)]}
        sensitivities = {key: {"delta": 0.1, "all_fits_converged": True} for key in MAJOR}
        sensitivities["exclude_recorded_zeros"] = {"delta": -99, "all_fits_converged": False}

        def decide(p=None, u=None, s=None, panel_valid=True):
            return A.verdict(primary if p is None else p, uncertainty if u is None else u,
                             sensitivities if s is None else s, panel_valid=panel_valid)

        self.assertEqual(decide()["code"], "transports")
        self.assertFalse(decide()["secondary_can_rescue_primary"])
        self.assertEqual(decide(panel_valid=False)["code"], "inconclusive")
        for change in ({"orders": 7}, {"all_fits_converged": False}, {"improving_groups": 4},
                       {"delta": 0}, {"delta": -0.01}, {"group_scores": []}):
            self.assertEqual(decide(p=dict(primary, **change))["code"], "inconclusive")
        for change in ({"ci_lower": 0}, {"ci_lower": -0.1}, {"deletions": []}):
            self.assertEqual(decide(u=dict(uncertainty, **change))["code"], "inconclusive")
        for key in MAJOR:
            for change in ({"delta": 0}, {"delta": -0.1}, {"all_fits_converged": False}):
                altered = copy.deepcopy(sensitivities)
                altered[key].update(change)
                self.assertEqual(decide(s=altered)["code"], "inconclusive")
            altered = copy.deepcopy(sensitivities)
            del altered[key]
            self.assertEqual(decide(s=altered)["code"], "inconclusive")
        altered = copy.deepcopy(uncertainty)
        altered["deletions"][0]["all_fits_converged"] = False
        self.assertEqual(decide(u=altered)["code"], "inconclusive")
        altered = copy.deepcopy(uncertainty)
        altered["deletions"][0]["deleted_order"] = "1"
        self.assertEqual(decide(u=altered)["code"], "inconclusive")
        for delta in (0, -0.1):
            negative = dict(primary, delta=delta)
            for upper, code in ((-0.01, "does_not_transport"), (0, "inconclusive"),
                                (0.01, "inconclusive")):
                answer = decide(p=negative, u=dict(uncertainty, ci_lower=-0.2, ci_upper=upper))
                self.assertEqual(answer["code"], code)
                self.assertTrue(answer["point_improvement_failure"])
                self.assertEqual(answer["interval_contains_zero"], upper >= 0)
        self.assertEqual(decide(u=dict(uncertainty, ci_lower=-0.2, ci_upper=-0.01))["code"],
                         "inconclusive")

    def test_10_newline_hashes_source_pins_and_portable_path_keys(self):
        paths = [Path("PROTOCOL.md"), Path("analysis.py")]
        paths += [Path("data") / name for name in (
            "species_aggregates.csv", "reconciliation.csv", "endpoint_audit.csv",
            "source_audit.json", "fetch_verify.py")]
        contents = {p.as_posix(): canonical_text(RESEARCH / p).encode("utf-8") for p in paths}
        expected = {key: hashlib.sha256(value).hexdigest() for key, value in contents.items()}
        self.assertEqual(expected, self.report["file_hashes_lf"])
        self.assertEqual(expected, A.input_hashes())
        self.assertEqual(expected[Path("PROTOCOL.md").as_posix()], PROTOCOL_SHA256)
        for key in expected:
            self.assertNotIn("\\", key)
            self.assertEqual(PureWindowsPath(key).as_posix(), key)
        for newline in (b"\n", b"\r\n", b"\r"):
            with tempfile.TemporaryDirectory() as folder:
                base = Path(folder)
                for key, content in contents.items():
                    target = base / key
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(content.replace(b"\n", newline))
                self.assertEqual(A.input_hashes(base), expected)
                self.assertEqual(A.load_rows(base), self.rows)
                self.assertEqual(A.prefit_summary(self.rows, base), self.report["prefit_summary"])
        self.assertEqual(F.normalized_sha256("a \r\n\r"), hashlib.sha256(b"a \n\n").hexdigest())
        self.assertNotEqual(F.normalized_sha256("a \r\n\r"), F.normalized_sha256("a"))
        audit = self.report["source_audit"]
        self.assertEqual(set(audit["sources"]), set(EXPECTED_SOURCES))
        self.assertEqual(set(F.SOURCE_FILES), set(EXPECTED_SOURCES))
        for name, expected_pins in EXPECTED_SOURCES.items():
            for metadata in (audit["sources"][name], F.SOURCE_FILES[name]):
                self.assertEqual(tuple(metadata[key] for key in
                                 ("commit", "git_blob_sha1", "sha256_lf")), expected_pins)
                self.assertEqual(Path(metadata["path"]).as_posix(), metadata["path"])
            self.assertEqual(PureWindowsPath(name).as_posix(), name)
        for name, digest in audit["derived_file_hashes_lf"].items():
            self.assertEqual(Path(name).as_posix(), name)
            self.assertNotIn("\\", name)
            self.assertEqual(hashlib.sha256(canonical_text(RESEARCH / "data" / name)
                             .encode("utf-8")).hexdigest(), digest)
        self.assertTrue(audit["licence_audit"]["checked_before_source_data_read"])
        self.assertIsNone(audit["licence_audit"]["vincze_repository_license"])
        self.assertIsNone(audit["licence_audit"]["compton_repository_license"])

    def test_11_synthetic_parser_rounding_and_no_mortality_eligibility(self):
        self.assertEqual(F.normalize_species("  ALPHA__beta\t Gamma  "), "alpha beta gamma")
        self.assertEqual(F.normalize_species("Straße_one"), "strasse one")
        self.assertNotEqual(F.normalize_species("Panthera leo"), F.normalize_species("Panthera leo persica"))
        for text in ("Species,x\nAlpha_one,1\n alpha ONE ,2\n",
                     "Species,x\n,1\n", "Species,x\nAlpha_one,1,2\n",
                     "Species,x\nAlpha_one\n"):
            with self.assertRaises(ValueError):
                F.parse_csv(text, "synthetic")
        for token in ("0", "0.0", "1", "1.0000"):
            self.assertEqual(F.rounding_tolerance(token), Decimal("2e-15"))
        self.assertEqual(F.rounding_tolerance("0.35416667"), Decimal("0.000000005000002"))
        self.assertEqual(F.rounding_tolerance("0.1"), Decimal("0.050000000000002"))
        for token in ("-0.1", "1.1", "NaN", "Infinity", "missing"):
            with self.assertRaises(ValueError):
                F.rounding_tolerance(token)
        for token, expected in (("0", 0), ("20.0", 20), ("2e1", 20)):
            self.assertEqual(F.integer(token, "synthetic"), expected)
        for token in ("1.5", "-1", "NaN", "Infinity", ""):
            with self.assertRaises(ValueError):
                F.integer(token, "synthetic")
        with self.assertRaises(ValueError):
            F.integer("0", "denominator", 1)
        toy_vincze = ("Species,Taxonomic order,Body mass (kg),Adult life expectancy,"
                      "Consuming mammal prey,MHC I gene family size\nAlpha_one,Order_A,2,10,0,3\n")
        toy_compton = ("Species,Class,Orders,RecordsWithDenominators,NeoplasiaWithDenominators,"
                       "Malignant,MalignancyKnown,NeoplasiaPrevalence,MalignancyPrevalence\n"
                       "alpha one,Mammalia,Order_A,20,0,0,0,0,0\n")
        texts = {Path("vincze.csv").as_posix(): toy_vincze,
                 Path("compton.csv").as_posix(): toy_compton,
                 Path("dictionary.txt").as_posix(): "Synthetic schema test only.\n"}
        metadata = copy.deepcopy(F.SOURCE_FILES)
        for name, text in texts.items():
            metadata[name]["sha256_lf"] = hashlib.sha256(text.encode("utf-8")).hexdigest()
        with mock.patch.object(F, "SOURCE_FILES", metadata):
            outputs = F.derive_outputs(texts)
            aggregates = list(csv.DictReader(io.StringIO(outputs[Path("species_aggregates.csv").as_posix()])))
            self.assertEqual(len(aggregates), 1)
            self.assertEqual(aggregates[0][K], "0")
            self.assertEqual(aggregates[0]["species"], "Alpha_one")
            self.assertFalse(any(field in toy_vincze for field in ("CMR", "RDI", "ICM")))
            for newline in ("\n", "\r\n", "\r"):
                with tempfile.TemporaryDirectory() as folder:
                    for name, text in texts.items():
                        (Path(folder) / name).write_bytes(text.replace("\n", newline).encode("utf-8"))
                    self.assertEqual(F.load_sources(Path(folder)), texts)
                    self.assertEqual(F.derive_outputs(F.load_sources(Path(folder))), outputs)
                    (Path(folder) / Path("vincze.csv")).write_text("corrupted", encoding="utf-8")
                    with self.assertRaises(ValueError):
                        F.load_sources(Path(folder))
            with tempfile.TemporaryDirectory() as folder:
                for name, content in outputs.items():
                    self.assertEqual(Path(name).as_posix(), name)
                    (Path(folder) / name).write_bytes(content.replace("\n", "\r\n").encode("utf-8"))
                F.compare_outputs(outputs, Path(folder))
            altered = dict(texts)
            altered[Path("vincze.csv").as_posix()] += "corrupted"
            with self.assertRaises(ValueError):
                F.derive_outputs(altered)
        with self.assertRaises(ValueError):
            F.load_sources()
        with self.assertRaises(ValueError):
            F.load_sources(Path("unused"), fetch=True)
        fixture = {"Species": "Alpha one", "RecordsWithDenominators": "20",
                   "NeoplasiaWithDenominators": "2", "Malignant": "1",
                   "MalignancyKnown": "2", "NeoplasiaPrevalence": "0.1",
                   "MalignancyPrevalence": "0.05"}
        audit = F.endpoint_audit_row(fixture)
        self.assertEqual(Decimal(audit["neoplasia_residual"]), 0)
        self.assertEqual(Decimal(audit["malignant_residual"]), 0)
        self.assertEqual(Decimal(audit["malignancyknown_residual"]), Decimal("-0.05"))
        with self.assertRaises(ValueError):
            F.endpoint_audit_row(dict(fixture, Malignant="3"))

    def test_12_diagnostic_changes_need_bounds_not_saved_stopping_values(self):
        changed = copy.deepcopy(self.report)
        changed["numerical_audit"]["maximum_fit_iterations"] = 100
        changed["numerical_audit"]["maximum_gradient_ratio"] = 0.99
        fit = changed["primary"]["folds"][0]["baseline_fit"]
        fit["iterations"] = 99
        fit["line_search_halvings"] = 44 * 99
        fit["maximum_absolute_gradient"] = fit["gradient_tolerance"] * 0.9
        self.assertEqual(A.compare_reports(self.report, changed), [])
        self.assert_tree(changed, self.report, snapshot=True)
        fit["maximum_absolute_gradient"] = fit["gradient_tolerance"] * 2
        self.assertTrue(A.compare_reports(self.report, changed))
        changed = copy.deepcopy(self.report)
        changed["primary"]["delta"] += 1e-3
        self.assertTrue(A.compare_reports(self.report, changed))

    def test_zz_dedicated_workload_under_sixty_seconds(self):
        self.assertLess(time.perf_counter() - self.started, 60)


if __name__ == "__main__":
    unittest.main()
