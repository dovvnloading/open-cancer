"""Offline reproducibility of the nonclinical MHC mortality transport audit.

No medical or veterinary advice. These checks establish numerical consistency,
not source truth, population independence, functional repertoires or causation.
"""

import ast
from contextlib import ExitStack
import copy
import csv
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path, PureWindowsPath
import sys
import tempfile
import time
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "mhc-cancer-mortality-transport"
FREEZE = "3695a1bfcb379a2ed37a654f2440995a9d488c95"
PROTOCOL_SHA256 = "aa0b6e05590d2016129b2207d0842ff89117ce00e8458c7219fb9fd9c1046c50"
SOURCE_COMMIT = "2c8a5c5b5fda6fb34c86b0e5a2b6a58dcc85cc3e"
EXPECTED_SOURCE_HASHES = {
    Path("SourceData_FINAL.csv").as_posix(): {
        "git_blob_sha1": "8a937a2c8d4f1749ae33c1b12d2ed52e895852f9",
        "sha256_lf": "b9d6dedc9447e9f99b275cbfbdb54e3b0e318d08c12f0384572a65101d057744",
    },
    Path("phylogeny.tre").as_posix(): {
        "git_blob_sha1": "64970f63fc167b7fdbadb5ebabbe3833dfb07c82",
        "sha256_lf": "23dab5f42ff4c6e736486fd54d115fff67b4b82d12a57a755f6577607b596c6a",
    },
}
N = "known_cause_deaths"
K = "nominal_cancer_deaths"
MASS = "log_body_mass_kg"
LIFE = "log_adult_life_expectancy_days"
MHC = "log1p_mhc_i"
DIET = "mammal_prey"
BEGIN = "<!-- BEGIN GENERATED RESULTS -->"
END = "<!-- END GENERATED RESULTS -->"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = load_module("research_mhc_transport_analysis", RESEARCH / "analysis.py")
F = load_module("research_mhc_transport_source", RESEARCH / "data" / "fetch_verify.py")


def canonical_text(path):
    # Text universal-newline decoding is independent of the analysis's explicit
    # CRLF/CR replacement and preserves trailing spaces and final-newline state.
    with path.open(encoding="utf-8", newline=None) as handle:
        return handle.read()


def independent_rows():
    with (RESEARCH / "data" / "species_aggregates.csv").open(
            encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for field in (N, K, "count_lower", "count_upper", DIET):
            row[field] = int(row[field])
        for field in (MASS, LIFE, MHC, "phylo_fair_proportion"):
            row[field] = float(row[field])
    return rows


def enumerated_counts(token, denominator):
    """Enumerate integer fractions; do not reuse either ceil/floor algorithm."""
    with localcontext() as context:
        context.prec = 65
        c = Decimal(token)
        n = Decimal(denominator)
        if c == 0 or c == 1:
            compatible = [int(c * n)]
        else:
            unit = Decimal(10) ** c.as_tuple().exponent
            half = unit / 2 + Decimal("0.000000000000002")
            lower, upper = max(Decimal(0), c - half), min(Decimal(1), c + half)
            compatible = [k for k in range(denominator + 1)
                          if lower <= Decimal(k) / n <= upper]
        nominal = int((c * n).to_integral_value(rounding=ROUND_HALF_EVEN))
        return nominal, compatible


def independent_probability(row, fit):
    features = [1.0, (row[MASS] - fit["means"][MASS]) / fit["scales"][MASS],
                (row[LIFE] - fit["means"][LIFE]) / fit["scales"][LIFE], float(row[DIET])]
    if fit["include_mhc"]:
        features.append((row[MHC] - fit["means"][MHC]) / fit["scales"][MHC])
    eta = math.fsum(x * b for x, b in zip(features, fit["coefficients"]))
    # The empirical range is moderate. The alternate branch keeps this helper
    # valid for the separate numerical boundary fixture as well.
    if eta < 0:
        exponential = math.exp(eta)
        p = exponential / (1 + exponential)
    else:
        p = 1 / (1 + math.exp(-eta))
    return features, p


def independent_loss(k, n, probability):
    """Compute the score with high-precision Decimal logarithms."""
    with localcontext() as context:
        context.prec = 50
        p = Decimal.from_float(min(1 - 1e-12, max(1e-12, probability)))
        event = Decimal(k) / Decimal(n)
        return float(-event * p.ln() - (1 - event) * (1 - p).ln())


def no_network():
    stack = ExitStack()
    for target in ("urllib.request.urlopen", "socket.create_connection", "socket.socket"):
        stack.enter_context(mock.patch(target, side_effect=AssertionError(
            "Network access is prohibited in this offline research test")))
    return stack


class MhcCancerMortalityTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.perf_counter()
        cls.rows = independent_rows()
        cls.evaluations = []
        evaluate = A.evaluate

        def retain_evaluation(rows, *args, **kwargs):
            result = evaluate(rows, *args, **kwargs)
            # Capture each complete run before build_report reduces some to
            # summaries. This permits independent checks of all 1188 fits and
            # every prediction without doing a second full fit or storing more
            # source material. Inputs are retained without modification.
            cls.evaluations.append((tuple(rows), result))
            return result

        with no_network(), mock.patch.object(A, "evaluate", side_effect=retain_evaluation):
            cls.report = A.build_report()
        cls.rebuild_seconds = time.perf_counter() - cls.started

    def assert_close(self, actual, expected, label="", rel=2e-10, absolute=3e-12):
        self.assertTrue(math.isfinite(actual), label)
        self.assertTrue(math.isfinite(expected), label)
        self.assertTrue(math.isclose(actual, expected, rel_tol=rel, abs_tol=absolute),
                        label or "{!r} != {!r}".format(actual, expected))

    def assert_tree(self, actual, expected, label="report"):
        if isinstance(expected, dict):
            self.assertIsInstance(actual, dict, label)
            self.assertEqual(set(actual), set(expected), label)
            for key in expected:
                self.assert_tree(actual[key], expected[key], label + "." + key)
        elif isinstance(expected, list):
            self.assertIsInstance(actual, list, label)
            self.assertEqual(len(actual), len(expected), label)
            for index, (left, right) in enumerate(zip(actual, expected)):
                self.assert_tree(left, right, label + "[" + str(index) + "]")
        elif isinstance(expected, float):
            self.assertIsInstance(actual, (float, int), label)
            self.assert_close(actual, expected, label)
        else:
            self.assertEqual(actual, expected, label)

    def test_01_entire_saved_report_and_exact_readme_block(self):
        def reject_nonfinite(value):
            raise AssertionError("Non-JSON numeric constant: " + value)

        saved = json.loads((RESEARCH / "data" / "results.json").read_text(encoding="utf-8"),
                           parse_constant=reject_nonfinite)
        self.assert_tree(saved, self.report)
        readme = canonical_text(RESEARCH / "README.md")
        self.assertEqual(readme.count(BEGIN), 1)
        self.assertEqual(readme.count(END), 1)
        block = readme[readme.index(BEGIN):readme.index(END) + len(END)] + "\n"
        self.assertEqual(block, A.render_results(self.report))
        self.assertIn(FREEZE, readme)
        self.assertIn("nonclinical", readme.lower())
        self.assertIn("no medical or veterinary advice", " ".join(readme.lower().split()))
        self.assertIn("not independent mortality validation", " ".join(readme.lower().split()))
        self.assertEqual(self.report["verdict"]["code"], "untestable")

    def test_02_counts_panel_and_source_audit_are_independently_recovered(self):
        rows = self.rows
        self.assertEqual(rows, A.load_rows())
        names = [r["species"] for r in rows]
        self.assertEqual(names, sorted(set(names)))
        self.assertTrue({"Heterocephalus_glaber", "Dolichotis_patagonum"} <= set(names))
        self.assertEqual(len(rows), 125)
        self.assertEqual(len({r["order"] for r in rows}), 16)
        self.assertEqual(sum(r[K] == 0 for r in rows), 17)
        ambiguous = {}
        residuals = []
        for row in rows:
            self.assertGreater(row[N], 0)
            nominal, candidates = enumerated_counts(row["cmr_decimal"], row[N])
            self.assertTrue(candidates, row["species"])
            self.assertIn(nominal, candidates, row["species"])
            expected = (nominal, min(candidates), max(candidates))
            self.assertEqual(expected, (row[K], row["count_lower"], row["count_upper"]))
            self.assertEqual(A.count_interval(row["cmr_decimal"], row[N]), expected)
            self.assertEqual(F.count_interval(row["cmr_decimal"], row[N]), expected)
            if len(candidates) > 1:
                ambiguous[row["species"]] = candidates
            residuals.append((abs(Decimal(row["cmr_decimal"]) * row[N] - row[K]),
                              row["species"]))
            self.assertIn(row[DIET], (0, 1))
            self.assertGreater(row["phylo_fair_proportion"], 0)
            self.assertGreaterEqual(row[MHC], 0)
            self.assert_close(math.expm1(row[MHC]), round(math.expm1(row[MHC])))
            self.assertNotIn("icm", " ".join(row).lower())
        self.assertEqual(ambiguous, {"Eulemur_macaco": list(range(2, 7)),
                                     "Myrmecophaga_tridactyla": list(range(1, 4))})
        audit = self.report["source_audit"]
        panel = self.report["panel"]
        for field, value in ((N, sum(r[N] for r in rows)), (K, sum(r[K] for r in rows)),
                             ("count_lower_total", sum(r["count_lower"] for r in rows)),
                             ("count_upper_total", sum(r["count_upper"] for r in rows))):
            self.assertEqual(panel[field], value)
            self.assertEqual(audit[field], value)
        self.assertEqual((panel[N], panel[K], panel["count_lower_total"],
                          panel["count_upper_total"]), (12033, 697, 694, 700))
        self.assertEqual(panel["largest_integer_residual"], str(max(residuals)[0]))
        self.assertEqual(panel["largest_integer_residual_species"], max(residuals)[1])
        self.assertEqual(audit["count_audit"]["unique_count_species"], len(rows) - len(ambiguous))
        with (RESEARCH / "data" / "exclusions.csv").open(encoding="utf-8", newline="") as stream:
            excluded = list(csv.DictReader(stream))
        self.assertEqual(len(excluded), 14)
        self.assertFalse(set(names) & {r["species"] for r in excluded})
        self.assertEqual(sum("MHC" in r["reason"] for r in excluded), 12)
        self.assertEqual({r["species"] for r in excluded if "MHC" not in r["reason"]},
                         {"Bos_gaurus", "Ovis_orientalis"})
        self.assertEqual(audit["source_rows"], len(rows) + len(excluded))
        for order in {r["order"] for r in rows}:
            self.assertEqual(len({r["phylo_block"] for r in rows if r["order"] == order}), 1)
        blocks = {b: sum(r["phylo_block"] == b for r in rows)
                  for b in sorted({r["phylo_block"] for r in rows})}
        self.assertEqual(blocks, {"clade_001": 119, "clade_002": 6})
        self.assertEqual(blocks, panel["species_per_phylo_block"])
        self.assertEqual(blocks, audit["phylogeny"]["clade_blocks"]["species_per_block"])
        self.assert_close(sum(r["phylo_fair_proportion"] for r in rows),
                          audit["phylogeny"]["retained_branch_length"])

    def check_fit(self, train, fit):
        """Certify each strict-convex ridge optimum from its gradient, independently."""
        fields = [MASS, LIFE] + ([MHC] if fit["include_mhc"] else [])
        means = {f: math.fsum(r[f] for r in train) / len(train) for f in fields}
        scales = {f: math.sqrt(math.fsum((r[f] - means[f]) ** 2 for r in train) /
                                   len(train)) or 1.0 for f in fields}
        self.assert_tree(fit["means"], means)
        self.assert_tree(fit["scales"], scales)
        self.assertEqual(fit["feature_names"], ["intercept", MASS, LIFE, DIET] +
                         ([MHC] if fit["include_mhc"] else []))
        beta, ridge = fit["coefficients"], fit["penalty"]
        self.assertEqual(len(beta), len(fit["feature_names"]))
        self.assertGreater(ridge, 0)
        mean_n = math.fsum(r[N] for r in train) / len(train)
        mean_fp = math.fsum(r["phylo_fair_proportion"] for r in train) / len(train)
        gradient_terms = [[] for _ in beta]
        likelihood = []
        effective_n = []
        for row in train:
            features, probability = independent_probability(row, fit)
            if fit["weighting"] == "counts":
                multiplier = 1.0
            elif fit["weighting"] == "species":
                multiplier = mean_n / row[N]
            else:
                self.assertEqual(fit["weighting"], "phylogeny")
                multiplier = row["phylo_fair_proportion"] / mean_fp
            n, k = row[N] * multiplier, row[K] * multiplier
            effective_n.append(n)
            for j, x in enumerate(features):
                gradient_terms[j].append((n * probability - k) * x)
            self.assertGreater(probability, 0)
            self.assertLess(probability, 1)
            likelihood.append(-k * math.log(probability) - (n - k) * math.log1p(-probability))
        gradient = [math.fsum(terms) + (ridge * beta[j] if j else 0)
                    for j, terms in enumerate(gradient_terms)]
        norm = max(abs(g) for g in gradient)
        self.assertLessEqual(norm, fit["gradient_tolerance"] * 1.01 + 1e-9)
        self.assert_close(norm, fit["maximum_absolute_gradient"], absolute=2e-9)
        objective = math.fsum(likelihood) + ridge * math.fsum(b * b for b in beta[1:]) / 2
        self.assert_close(objective, fit["objective"], absolute=2e-9)
        self.assertEqual(fit["training_species"], len(train))
        self.assertEqual(fit["training_known_cause_deaths"], sum(r[N] for r in train))
        self.assert_close(fit["effective_training_denominator"], math.fsum(effective_n))
        self.assertTrue(fit["converged"])
        self.assertLessEqual(fit["iterations"], 100)

    def test_03_every_fit_has_training_only_scaling_weights_and_stationarity(self):
        fit_count = 0
        for rows, result in self.evaluations:
            lookup = {r["species"]: r for r in rows}
            seen_test = []
            for fold in result["folds"]:
                train_ids, test_ids = fold["training_species_ids"], fold["test_species_ids"]
                self.assertFalse(set(train_ids) & set(test_ids))
                self.assertEqual(set(train_ids) | set(test_ids), set(lookup))
                train = [lookup[s] for s in train_ids]
                test = [lookup[s] for s in test_ids]
                self.assertFalse({r["order"] for r in train} & {r["order"] for r in test})
                self.assertEqual({r[result["group_key"]] for r in test}, {fold["heldout_group"]})
                self.assertNotIn(fold["heldout_group"], {r[result["group_key"]] for r in train})
                self.assertEqual(fold["training_orders"], sorted({r["order"] for r in train}))
                self.assertEqual(fold["test_orders"], sorted({r["order"] for r in test}))
                seen_test.extend(test_ids)
                for key, mhc in (("baseline_fit", False), ("extended_fit", True)):
                    fit = fold[key]
                    self.assertEqual(fit["include_mhc"], mhc)
                    self.assertEqual(fit["penalty"], result["penalty"])
                    self.assertEqual(fit["weighting"], result["weighting"])
                    self.check_fit(train, fit)
                    fit_count += 1
            self.assertEqual(sorted(seen_test), sorted(lookup))
        self.assertEqual(fit_count, self.report["numerical_audit"]["fits"])

    def test_04_every_prediction_and_all_score_aggregations(self):
        for rows, result in self.evaluations:
            lookup = {r["species"]: r for r in rows}
            folds = {f["heldout_group"]: f for f in result["folds"]}
            groups = {}
            weighted = [[], []]
            for prediction in result["predictions"]:
                row = lookup[prediction["species"]]
                fold = folds[row[result["group_key"]]]
                losses = []
                for index, model in enumerate(("baseline", "extended")):
                    _, p = independent_probability(row, fold[model + "_fit"])
                    loss = independent_loss(row[K], row[N], p)
                    self.assert_close(prediction[model + "_probability"], p)
                    self.assert_close(prediction[model + "_loss"], loss)
                    losses.append(loss)
                    weighted[index].append(row[N] * loss)
                self.assertEqual((prediction[N], prediction[K]), (row[N], row[K]))
                self.assert_close(prediction["delta"], losses[0] - losses[1])
                groups.setdefault(prediction["heldout_group"], []).append((row, losses))
            group_means = [[], []]
            positive = 0
            for record in result["group_scores"]:
                members = groups[record["group"]]
                means = [math.fsum(item[1][i] for item in members) / len(members)
                         for i in range(2)]
                self.assertEqual(record["species"], len(members))
                self.assertEqual(record[N], sum(r[N] for r, _ in members))
                self.assertEqual(record[K], sum(r[K] for r, _ in members))
                self.assertEqual(record["zero_case_species"], sum(r[K] == 0 for r, _ in members))
                self.assert_close(record["baseline_loss"], means[0])
                self.assert_close(record["extended_loss"], means[1])
                self.assert_close(record["delta"], means[0] - means[1])
                positive += int(means[0] > means[1])
                for i in range(2):
                    group_means[i].append(means[i])
            means = [math.fsum(values) / len(groups) for values in group_means]
            death_means = [math.fsum(values) / sum(r[N] for r in rows) for values in weighted]
            for index, model in enumerate(("baseline", "extended")):
                self.assert_close(result[model + "_loss"], means[index])
                self.assert_close(result["death_weighted_" + model + "_loss"], death_means[index])
            self.assert_close(result["delta"], means[0] - means[1])
            self.assert_close(result["death_weighted_delta"], death_means[0] - death_means[1])
            self.assertEqual(result["improving_groups"], positive)
            self.assertEqual(result["species"], len(rows))
            self.assertEqual(result["groups"], len(groups))
            self.assertEqual(result[N], sum(r[N] for r in rows))
            self.assertEqual(result[K], sum(r[K] for r in rows))
            self.assertEqual(result["fit_count"], 2 * len(groups))

    def test_05_complete_refit_jackknife_and_student_quantile(self):
        uncertainty = self.report["uncertainty"]
        deletions = uncertainty["deletions"]
        orders = sorted({r["order"] for r in self.rows})
        self.assertEqual([d["deleted_order"] for d in deletions], orders)
        with localcontext() as context:
            context.prec = 60
            g = Decimal(len(orders))
            delta = Decimal(str(self.report["primary"]["delta"]))
            pseudo = [g * delta - (g - 1) * Decimal(str(d["delta"])) for d in deletions]
            center = sum(pseudo) / g
            se = (sum((q - center) ** 2 for q in pseudo) / (g * (g - 1))).sqrt()
            critical = Decimal("2.131449545559323")
            for row, expected in zip(deletions, pseudo):
                self.assert_close(row["pseudovalue"], float(expected))
                self.assertEqual(row["groups"], 15)
                self.assertNotIn(row["deleted_order"], {s["group"] for s in row["group_scores"]})
                self.assertEqual(row["species"], sum(r["order"] != row["deleted_order"]
                                                     for r in self.rows))
            for key, expected in (("bias_corrected_center", center), ("standard_error", se),
                                  ("ci_lower", center - critical * se),
                                  ("ci_upper", center + critical * se)):
                self.assert_close(uncertainty[key], float(expected))
        self.assert_close(uncertainty["delete_order_delta_min"], min(d["delta"] for d in deletions))
        self.assert_close(uncertainty["delete_order_delta_max"], max(d["delta"] for d in deletions))
        # Numerically integrate the t15 density instead of trusting the stored
        # critical value solely because two constants agree.
        degrees = 15
        t = uncertainty["t_critical"]
        normalizer = math.gamma((degrees + 1) / 2) / (
            math.sqrt(degrees * math.pi) * math.gamma(degrees / 2))
        density = lambda x: normalizer * (1 + x * x / degrees) ** (-(degrees + 1) / 2)
        subdivisions = 2000
        step = t / subdivisions
        area = step / 3 * (density(0) + density(t) + math.fsum(
            (4 if i % 2 else 2) * density(i * step) for i in range(1, subdivisions)))
        self.assert_close(0.5 + area, 0.975, absolute=2e-12)

    def test_06_rounding_and_every_sensitivity_preserve_the_declared_design(self):
        grid = self.report["rounding"]
        observed = {(r["assignment"]["Eulemur_macaco"],
                     r["assignment"]["Myrmecophaga_tridactyla"]) for r in grid}
        self.assertEqual(observed, set(itertools.product(range(2, 7), range(1, 4))))
        nominal = [r for r in grid if r["assignment"] ==
                   {"Eulemur_macaco": 4, "Myrmecophaga_tridactyla": 2}]
        self.assertEqual(len(nominal), 1)
        self.assert_close(nominal[0]["delta"], self.report["primary"]["delta"])
        for item in grid:
            self.assertEqual(item[K], sum(r[K] for r in self.rows) - 6 + sum(item["assignment"].values()))
            self.assertEqual(item["species"], len(self.rows))
            self.assertEqual(item["groups"], 16)
        summary = self.report["rounding_summary"]
        self.assertEqual(summary["scenarios"], len(grid))
        self.assert_close(summary["minimum_delta"], min(r["delta"] for r in grid))
        self.assert_close(summary["maximum_delta"], max(r["delta"] for r in grid))
        self.assertEqual(summary["positive_scenarios"], sum(r["delta"] > 0 for r in grid))
        sensitivities = self.report["sensitivities"]
        self.assertEqual(set(sensitivities), {"exclude_named_rodents", "exclude_recorded_zeros",
                         "equal_species_training", "penalty_0_1", "penalty_10",
                         "phylogeny_weighting", "deeper_phylogeny_holdout"})
        self.assertEqual(sensitivities["exclude_named_rodents"]["species"], len(self.rows) - 2)
        self.assertEqual(sensitivities["exclude_recorded_zeros"]["species"],
                         sum(r[K] > 0 for r in self.rows))
        self.assertEqual(sensitivities["exclude_recorded_zeros"]["zero_case_species"], 0)
        for key, weighting in (("equal_species_training", "species"),
                               ("phylogeny_weighting", "phylogeny")):
            self.assertEqual(sensitivities[key]["weighting"], weighting)
            self.assertEqual(sensitivities[key]["species"], len(self.rows))
            self.assertEqual(sensitivities[key]["penalty"], 1)
        self.assertEqual(sensitivities["penalty_0_1"]["penalty"], 0.1)
        self.assertEqual(sensitivities["penalty_10"]["penalty"], 10)
        deep = sensitivities["deeper_phylogeny_holdout"]
        self.assertEqual(deep["group_key"], "phylo_block")
        self.assertEqual(deep["groups"], 2)
        self.assertFalse(self.report["verdict"]["positive_claim_gates"]["deeper_tree_gate_applies"])
        audit = self.report["numerical_audit"]
        self.assertEqual(audit["complete_evaluations"], len(self.evaluations))
        self.assertEqual(audit["fits"], sum(r["fit_count"] for _, r in self.evaluations))
        self.assertEqual(audit["maximum_fit_iterations"],
                         max(r["maximum_fit_iterations"] for _, r in self.evaluations))
        self.assert_close(audit["maximum_gradient_ratio"],
                          max(r["maximum_gradient_ratio"] for _, r in self.evaluations))
        self.assertEqual(audit["probabilities_clipped_for_scoring"],
                         sum(r["probabilities_clipped_for_scoring"] for _, r in self.evaluations))

    def test_07_heldout_outcome_perturbation_cannot_change_its_predictions(self):
        heldout = "Rodentia"
        altered = [dict(r, **{K: r[N] - r[K]}) if r["order"] == heldout else dict(r)
                   for r in self.rows]
        with no_network():
            result = A.evaluate(altered)
        before_fold = next(f for f in self.report["primary"]["folds"] if f["heldout_group"] == heldout)
        after_fold = next(f for f in result["folds"] if f["heldout_group"] == heldout)
        self.assertEqual(before_fold, after_fold)
        original = {p["species"]: p for p in self.report["primary"]["predictions"]
                    if p["order"] == heldout}
        changed = {p["species"]: p for p in result["predictions"] if p["order"] == heldout}
        self.assertEqual(set(original), set(changed))
        for species in original:
            for field in ("baseline_probability", "extended_probability"):
                self.assertEqual(original[species][field], changed[species][field])
        self.assertTrue(any(original[s]["baseline_loss"] != changed[s]["baseline_loss"]
                            for s in original))

    def test_08_zero_events_numerical_boundaries_and_degenerate_fits(self):
        self.assert_close(A.cancer_death_loss(0, 100, 0.2), -math.log(0.8))
        self.assert_close(A.cancer_death_loss(100, 100, 0.2), -math.log(0.2))
        self.assertGreater(A.cancer_death_loss(0, 100, 0.2), 0)
        for k, n, p in ((0, 10, 0), (0, 10, 1), (10, 10, 0), (10, 10, 1), (4, 19, 0.3)):
            self.assert_close(A.cancer_death_loss(k, n, p), independent_loss(k, n, p))
        for k, n in ((1, 0), (-1, 10), (11, 10)):
            with self.assertRaises(ValueError):
                A.cancer_death_loss(k, n, 0.2)
        for endpoint in (0, 1):
            constant = [dict(r, **{K: endpoint * r[N]}) for r in self.rows[:8]]
            with self.assertRaises(A.FitError):
                A.fit_model(constant, True)
        # All continuous columns constant and binary diet separates the outcomes.
        # The unpenalized intercept is identifiable; ridge prevents infinite slopes.
        toy = [{MASS: 1.0, LIFE: 2.0, MHC: 3.0, DIET: i % 2,
                N: 10, K: 10 * (i % 2), "phylo_fair_proportion": 1.0}
               for i in range(8)]
        fit = A.fit_model(toy, True)
        self.assertEqual(set(fit["scales"].values()), {1.0})
        self.check_fit(toy, fit)
        self.assertTrue(all(math.isfinite(x) for x in fit["coefficients"]))

    def test_09_verdict_rules_pass_fail_and_abstain_without_vacuous_success(self):
        primary = {"orders": 16, "groups": 16, "delta": 1.0,
                   "improving_groups": 9, "all_fits_converged": True}
        uncertainty = {"ci_lower": 0.5, "ci_upper": 1.5,
                       "deletions": [{"all_fits_converged": True} for _ in range(16)]}
        rounding = [{"delta": 0.2, "all_fits_converged": True} for _ in range(15)]
        sensitivities = {key: {"delta": 0.1, "all_fits_converged": True}
                         for key in ("exclude_named_rodents", "equal_species_training",
                                     "penalty_0_1", "penalty_10", "phylogeny_weighting")}
        sensitivities["deeper_phylogeny_holdout"] = {
            "groups": 2, "delta": -10.0, "all_fits_converged": True}
        sensitivities["exclude_recorded_zeros"] = {"delta": -10.0, "all_fits_converged": True}

        def decide(p=None, u=None, r=None, s=None, panel_valid=True):
            return A.verdict(primary if p is None else p, uncertainty if u is None else u,
                             rounding if r is None else r, sensitivities if s is None else s,
                             panel_valid=panel_valid)["code"]

        self.assertEqual(decide(), "adds")
        self.assertEqual(decide(panel_valid=False), "untestable")
        for changed in (dict(primary, orders=7), dict(primary, all_fits_converged=False),
                        dict(primary, improving_groups=8)):
            self.assertEqual(decide(p=changed), "untestable")
        self.assertEqual(decide(u=dict(uncertainty, ci_lower=-0.1)), "untestable")
        self.assertEqual(decide(u=dict(uncertainty, deletions=[])), "untestable")
        failed = copy.deepcopy(uncertainty)
        failed["deletions"][0]["all_fits_converged"] = False
        self.assertEqual(decide(u=failed), "untestable")
        self.assertEqual(decide(r=[]), "untestable")
        self.assertEqual(decide(r=rounding[:-1]), "untestable")
        for change in ({"delta": 0}, {"all_fits_converged": False}):
            altered = copy.deepcopy(rounding)
            altered[0].update(change)
            self.assertEqual(decide(r=altered), "untestable")
        for key in ("exclude_named_rodents", "equal_species_training", "penalty_0_1",
                    "penalty_10", "phylogeny_weighting"):
            for change in ({"delta": 0}, {"all_fits_converged": False}):
                altered = copy.deepcopy(sensitivities)
                altered[key].update(change)
                self.assertEqual(decide(s=altered), "untestable")
        altered = copy.deepcopy(sensitivities)
        altered["deeper_phylogeny_holdout"].update(groups=4)
        self.assertEqual(decide(s=altered), "untestable")
        altered["deeper_phylogeny_holdout"].update(delta=0.1)
        self.assertEqual(decide(s=altered), "adds")
        altered["deeper_phylogeny_holdout"].update(all_fits_converged=False)
        self.assertEqual(decide(s=altered), "untestable")
        negative = dict(primary, delta=-0.2)
        self.assertEqual(decide(p=negative, u=dict(uncertainty, ci_lower=-0.3, ci_upper=-0.1)),
                         "does_not_add")
        self.assertEqual(decide(p=negative, u=dict(uncertainty, ci_lower=-0.3, ci_upper=0.1)),
                         "untestable")

    def test_10_normalized_hashes_windows_paths_and_full_crlf_recomputation(self):
        paths = [Path("PROTOCOL.md"), Path("analysis.py")]
        paths += [Path("data") / name for name in (
            "species_aggregates.csv", "exclusions.csv", "source_audit.json", "fetch_verify.py")]
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
                    target = base / Path(key)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(content.replace(b"\n", newline))
                self.assertEqual(A.input_hashes(base), expected)
                if newline == b"\r\n":
                    with no_network():
                        regenerated = A.build_report(base)
                    self.assert_tree(regenerated, self.report)
        self.assertEqual(F.normalized_sha256("a \r\n\r"), hashlib.sha256(b"a \n\n").hexdigest())
        self.assertNotEqual(F.normalized_sha256("a \r\n\r"), F.normalized_sha256("a"))
        audit = self.report["source_audit"]
        self.assertEqual(audit["source_commit"], SOURCE_COMMIT)
        self.assertEqual(F.SOURCE_COMMIT, SOURCE_COMMIT)
        self.assertEqual(audit["source_files"], EXPECTED_SOURCE_HASHES)
        self.assertEqual(F.SOURCE_FILES, EXPECTED_SOURCE_HASHES)
        for key, expected_hash in audit["derived_file_hashes_lf"].items():
            self.assertEqual(Path(key).as_posix(), key)
            self.assertEqual(hashlib.sha256(canonical_text(RESEARCH / "data" / key).encode("utf-8"))
                             .hexdigest(), expected_hash)

    def test_11_synthetic_nexus_fair_proportion_and_order_respecting_blocks(self):
        fixture = """#NEXUS
BEGIN TAXA;
 DIMENSIONS NTAX=4;
END;
BEGIN TREES;
 TRANSLATE 1 Alpha_one, 2 Beta_two, 3 Gamma_three, 4 Delta_four;
 TREE test = [&R] ((1:1,2:1)Internal:2,(3:2,4:2):1):0;
END;
"""
        tree, names = F.parse_nexus(fixture)
        self.assertEqual(set(names.values()), {"Alpha_one", "Beta_two", "Gamma_three", "Delta_four"})
        fp, length = F.fair_proportion(tree, list(names.values()))
        self.assertEqual(fp, {"Alpha_one": 2.0, "Beta_two": 2.0,
                              "Gamma_three": 2.5, "Delta_four": 2.5})
        self.assertEqual(length, 9.0)
        reduced, length = F.fair_proportion(tree, ["Alpha_one", "Gamma_three", "Delta_four"])
        self.assertEqual(reduced, {"Alpha_one": 3.0, "Gamma_three": 2.5, "Delta_four": 2.5})
        self.assertEqual(length, 8.0)
        orders = {"Alpha_one": "Order_A", "Beta_two": "Order_B",
                  "Gamma_three": "Order_C", "Delta_four": "Order_C"}
        blocks, audit = F.phylogenetic_blocks(tree, orders)
        self.assertEqual(audit["half_height_threshold"], 1.5)
        self.assertEqual(audit["initial_crossing_edge_groups"], 3)
        self.assertEqual(audit["merged_block_count"], 2)
        self.assertEqual(blocks["Alpha_one"], blocks["Beta_two"])
        self.assertEqual(blocks["Gamma_three"], blocks["Delta_four"])
        self.assertNotEqual(blocks["Alpha_one"], blocks["Gamma_three"])
        orders.update(Gamma_three="Order_A", Delta_four="Order_B")
        blocks, audit = F.phylogenetic_blocks(tree, orders)
        self.assertEqual(len(set(blocks.values())), 1)  # Transitive order unions.
        self.assertTrue(audit["orders_are_not_split"])
        for bad in (fixture.replace("1:1", "1:-1"),
                    fixture.replace("Delta_four", "Gamma_three"),
                    fixture.replace("):0;", "):1;")):
            with self.assertRaises(ValueError):
                F.parse_nexus(bad)
        with self.assertRaises(ValueError):
            F.fair_proportion(tree, ["Missing_species"])
        with self.assertRaises(ValueError):
            F.load_sources()  # Retrieval is explicitly opt-in.

    def test_12_standard_library_and_frozen_scientific_boundaries(self):
        paths = [RESEARCH / "analysis.py", RESEARCH / "data" / "fetch_verify.py", Path(__file__)]
        for path in paths:
            tree = ast.parse(canonical_text(path))
            imports = []
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imports.extend(alias.name.split(".")[0] for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imports.append(node.module.split(".")[0])
            self.assertTrue(set(imports) <= sys.stdlib_module_names, path.as_posix())
        protocol = canonical_text(RESEARCH / "PROTOCOL.md")
        self.assertEqual(hashlib.sha256(protocol.encode("utf-8")).hexdigest(), PROTOCOL_SHA256)
        plain_protocol = " ".join(protocol.split())
        for phrase in ("not independent mortality validation", "functional", "immune causation",
                       "Never recover an integer count from ICM", "No medical or veterinary advice"):
            self.assertIn(phrase, plain_protocol)
        self.assertIsNone(self.report["source_audit"]["reuse"]["repository_license"])
        self.assertEqual(self.report["source_audit"]["reuse"]["article_license"], "CC-BY-NC-ND-4.0")
        self.assertFalse((RESEARCH / "data" / "SourceData_FINAL.csv").exists())
        self.assertFalse((RESEARCH / "data" / "phylogeny.tre").exists())

    def test_99_entire_dedicated_suite_remains_below_one_minute(self):
        self.assertLess(self.rebuild_seconds, 60)
        self.assertLess(time.perf_counter() - self.started, 60)


if __name__ == "__main__":
    unittest.main()
