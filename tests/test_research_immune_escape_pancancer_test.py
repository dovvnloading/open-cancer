"""Offline checks of NONCLINICAL cancer-type aggregate calculations.

The references use pairwise-sum OLS, direct probability predictions, parent
weights, and nested convex minimization rather than the production active-face
solver. Source-field checks cannot establish source truth, permissions,
independent biological replication or clinical usefulness. No patient records,
network requests or third-party packages are used.
"""

import csv
from decimal import Decimal
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import tempfile
import time
import unittest
from unittest.mock import patch


RESEARCH = (Path(__file__).resolve().parents[1] / "research"
            / "immune-escape-pancancer-test")
SPEC = importlib.util.spec_from_file_location(
    "immune_escape_pancancer_analysis", RESEARCH / "analysis.py")
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)


def bounded(value, lower=0.0, upper=1.0):
    return min(upper, max(lower, value))


def average(values):
    values = list(values)
    return math.fsum(values) / len(values)


def percentile(values, probability):
    """Explicit interpolation order statistic, without production quantiles."""
    ordered = sorted(values)
    if not ordered:
        return None
    index = probability * (len(ordered) - 1)
    left, right = math.floor(index), math.ceil(index)
    return ordered[left] + (index-left) * (ordered[right]-ordered[left])


def ols_reference(rows, model):
    """Pairwise covariance identity, independently of centered OLS code."""
    x = [0.0 if model == "mean" else
         math.log(row["tmb"]) if model == "log_tmb" else row["tmb"]
         for row in rows]
    y = [row["orr"] for row in rows]
    denominator = math.fsum((x[i]-x[j])**2
                            for i in range(len(x)) for j in range(i))
    numerator = math.fsum((x[i]-x[j])*(y[i]-y[j])
                          for i in range(len(x)) for j in range(i))
    slope = numerator / denominator if denominator else 0.0
    return {"model": model, "intercept": average(y)-slope*average(x),
            "slope": slope}


def prediction_reference(parameters, row):
    model = parameters["model"]
    if model in ("mean", "linear_tmb", "log_tmb"):
        x = (0.0 if model == "mean" else math.log(row["tmb"])
             if model == "log_tmb" else row["tmb"])
        return bounded(parameters["intercept"] + parameters["slope"] * x)
    load = row["tmb"]
    if model == "power_tmb":
        load **= parameters["q"]
    elif model == "power_clonal":
        load *= row["fraction"] ** parameters["theta"]
    elif model in ("reciprocal_clonal", "reciprocal_full", "free_full"):
        load *= row["fraction"]
    uses_hla = model in ("reciprocal_hla", "reciprocal_full", "free_hla",
                         "free_full") or parameters.get("with_hla", False)
    prevalence = row["loh_fraction"] if uses_hla else 0.0
    no_loh = parameters["A0"] * math.exp(-parameters["B"] / load)
    with_loh = parameters["A1"] * math.exp(-parameters["B"] / load)
    return (1-prevalence)*no_loh + prevalence*with_loh


def metric_reference(rows, predicted):
    """Explicit row weights giving every parent cancer total weight 1/G."""
    counts = {}
    for row in rows:
        group = row["parent_type"]
        counts[group] = counts.get(group, 0) + 1
    errors = [row["orr"]-prediction for row, prediction in zip(rows, predicted)]
    weights = [1/(len(counts)*counts[row["parent_type"]]) for row in rows]
    mse = math.fsum(w*e*e for w, e in zip(weights, errors))
    return {"mse": mse, "rmse": math.sqrt(mse),
            "mae": math.fsum(w*abs(e) for w, e in zip(weights, errors)),
            "mean_residual": math.fsum(w*e for w, e in zip(weights, errors)),
            "row_rmse": math.sqrt(average(e*e for e in errors)),
            "group_mse": {group: average(errors[i]**2 for i, row in enumerate(rows)
                                          if row["parent_type"] == group)
                          for group in counts}}


def ceiling_reference(z, outcomes, loh=None, adverse=True):
    """Nested convex profile, independent of enumerating feasible faces.

    At fixed A0, the least-squares A1 is a clipped scalar regression. Its
    profile in A0 is convex on [0,1], so ternary minimization plus endpoints
    covers the entire triangle/square without active-set enumeration.
    """
    if loh is None:
        denominator = math.fsum(value**2 for value in z)
        a = bounded(math.fsum(x*y for x, y in zip(z, outcomes))/denominator) if denominator else 0.0
        return math.fsum((y-a*x)**2 for x, y in zip(z, outcomes)), a, a
    u = [value*(1-prevalence) for value, prevalence in zip(z, loh)]
    v = [value*prevalence for value, prevalence in zip(z, loh)]
    uu, vv = math.fsum(x*x for x in u), math.fsum(x*x for x in v)
    uv = math.fsum(x*y for x, y in zip(u, v))
    uy = math.fsum(x*y for x, y in zip(u, outcomes))
    vy = math.fsum(x*y for x, y in zip(v, outcomes))
    yy = math.fsum(y*y for y in outcomes)

    def value(a0):
        a1 = bounded((vy-uv*a0)/vv, 0.0, a0 if adverse else 1.0) if vv else 0.0
        loss = yy - 2*a0*uy - 2*a1*vy + a0*a0*uu + 2*a0*a1*uv + a1*a1*vv
        return loss, a0, a1

    low, high = 0.0, 1.0
    for _ in range(64):
        first, second = (2*low+high)/3, (low+2*high)/3
        if value(first)[0] < value(second)[0]:
            high = second
        else:
            low = first
    _, a0, a1 = min((value(0.0), value(1.0), value((low+high)/2)),
                    key=lambda item: item[0])
    loss = math.fsum((y-a0*x-a1*w)**2 for y, x, w in zip(outcomes, u, v))
    return loss, a0, a1


@lru_cache(maxsize=None)
def reciprocal_reference(loads, outcomes, loh=None, adverse=True):
    """A wider independent log-B grid and ternary refinement of each basin."""
    def objective(b):
        z = [math.exp(-b/load) for load in loads]
        return ceiling_reference(z, outcomes, loh, adverse)

    first, last = math.log(min(loads)*1e-9), math.log(100*max(loads))
    logs = [first+(last-first)*i/100 for i in range(101)]
    values = [objective(math.exp(value)) for value in logs]
    candidates = [(objective(0.0), 0.0)] + list(zip(values, map(math.exp, logs)))
    for i in range(1, len(logs)-1):
        if (values[i][0] <= values[i-1][0] and values[i][0] <= values[i+1][0]
                and (values[i][0] < values[i-1][0] or values[i][0] < values[i+1][0])):
            low, high = logs[i-1], logs[i+1]
            for _ in range(46):
                left, right = (2*low+high)/3, (low+2*high)/3
                if objective(math.exp(left))[0] < objective(math.exp(right))[0]:
                    high = right
                else:
                    low = left
            b = math.exp((low+high)/2)
            candidates.append((objective(b), b))
    (loss, a0, a1), b = min(candidates, key=lambda item: item[0][0])
    return {"A0": a0, "A1": a1, "B": b, "sse": loss}


def nonlinear_inputs(rows, model, parameters=None):
    parameters = parameters or {}
    loads = []
    for row in rows:
        load = row["tmb"]
        if model == "power_tmb":
            load **= parameters["q"]
        elif model == "power_clonal":
            load *= row["fraction"] ** parameters["theta"]
        elif model in ("reciprocal_clonal", "reciprocal_full", "free_full"):
            load *= row["fraction"]
        loads.append(load)
    uses_hla = model in ("reciprocal_hla", "reciprocal_full", "free_hla",
                         "free_full") or parameters.get("with_hla", False)
    loh = tuple(row["loh_fraction"] for row in rows) if uses_hla else None
    return tuple(loads), tuple(row["orr"] for row in rows), loh, not model.startswith("free_")


def synthetic_rows():
    values = [0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 13.0, 17.0]
    prevalence = [0.05, 0.30, 0.20, 0.40, 0.10, 0.25,
                  0.15, 0.35, 0.45, 0.17, 0.27, 0.37]
    return [{"id": str(i), "parent_type": str(i//2), "tmb": load,
             "fraction": 0.4+0.02*i, "loh_fraction": prevalence[i],
             "orr": (0.7*(1-prevalence[i])+0.2*prevalence[i])*math.exp(-1.2/load)}
            for i, load in enumerate(values)]


class TestPanCancerSyntheticGuards(unittest.TestCase):
    def test_probability_bridge_identity_and_direction(self):
        parameters = {"model": "reciprocal_full", "A0": 0.7, "A1": 0.2, "B": 1.2}
        row = {"tmb": 5.0, "fraction": 0.6, "loh_fraction": 0.3}
        first = analysis.predict(parameters, row)
        doubled = analysis.predict(parameters, dict(row, tmb=10.0))
        ceiling = 0.7*0.7+0.2*0.3
        self.assertAlmostEqual(doubled/ceiling, math.sqrt(first/ceiling), places=14)
        self.assertGreater(analysis.predict(parameters, dict(row, fraction=0.8)), first)
        self.assertLess(analysis.predict(parameters, dict(row, loh_fraction=0.6)), first)

    def test_boundary_projection_independent_reference(self):
        generator = random.Random(917)
        for adverse in (False, True):
            for _ in range(12):
                z = [generator.random() for _ in range(9)]
                y = [generator.random() for _ in range(9)]
                loh = [generator.random() for _ in range(9)]
                got = analysis._bounded_ceiling(z, y, loh, adverse)
                reference = ceiling_reference(z, y, loh, adverse)
                self.assertAlmostEqual(got[0], reference[0], places=10)
                self.assertTrue(0 <= got[1] <= 1 and 0 <= got[2] <= 1)
                if adverse:
                    self.assertLessEqual(got[2], got[1])
        for y in ([0.0]*4, [1.0]*4):
            got = analysis._bounded_ceiling([1.0]*4, y, [0.0, 0.25, 0.75, 1.0])
            self.assertLess(got[0], 1e-20)

    def test_exact_synthetic_curve_and_unidentified_designs(self):
        rows = synthetic_rows()
        fitted = analysis.fit(rows, "reciprocal_hla")
        for key, expected in (("A0", 0.7), ("A1", 0.2), ("B", 1.2)):
            self.assertAlmostEqual(fitted[key], expected, places=7)
        fixed = [dict(row, loh_fraction=0.2) for row in rows]
        self.assertTrue(analysis.fit(fixed, "reciprocal_hla")["singular"])
        self.assertTrue(analysis.reciprocal_fit([1.0, 1.0], [0.2, 0.3])["singular"])
        zero = analysis.reciprocal_fit([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
        self.assertEqual(zero["sse"], 0.0)
        self.assertTrue(zero["singular"])

    def test_group_holdout_and_equal_parent_scoring(self):
        rows = synthetic_rows()
        seen = []
        original = analysis.fit

        def record_fit(training, model, **kwargs):
            seen.append({row["parent_type"] for row in training})
            return original(training, model, **kwargs)

        with patch.object(analysis, "fit", side_effect=record_fit):
            result = analysis.cross_validate(rows, "log_tmb")
        groups = list(dict.fromkeys(row["parent_type"] for row in rows))
        self.assertEqual(len(seen), len(groups))
        for omitted, train_groups in zip(groups, seen):
            self.assertNotIn(omitted, train_groups)
            self.assertEqual(train_groups, set(groups)-{omitted})
        predicted = [result["predictions"][row["id"]] for row in rows]
        expected = metric_reference(rows, predicted)
        for key in ("mse", "rmse", "mae", "mean_residual", "row_rmse"):
            self.assertAlmostEqual(result[key], expected[key], places=13)
        unequal = rows[:3] + rows[4:]
        zero_predictions = [0.0]*len(unequal)
        result = analysis.scores(unequal, zero_predictions)
        self.assertAlmostEqual(result["mse"], metric_reference(unequal, zero_predictions)["mse"], places=14)

    def test_missing_values_invalid_probabilities_and_design_guard(self):
        rows = synthetic_rows()
        with self.assertRaises(ValueError):
            analysis.fit([dict(rows[0], fraction=None)], "reciprocal_clonal")
        with self.assertRaises(ValueError):
            analysis.fit([dict(rows[0], loh_fraction=None)], "reciprocal_hla")
        for loads, outcomes, loh in (([0.0], [0.2], None), ([math.nan], [0.2], None),
                                     ([1.0], [1.01], None), ([1.0], [math.nan], None),
                                     ([1.0], [0.2], [-0.1])):
            with self.assertRaises(ValueError):
                analysis.reciprocal_fit(loads, outcomes, loh)
        with self.assertRaises(ValueError):
            analysis.scores(rows, [0.0])
        with self.assertRaises(ValueError):
            analysis.cross_validate(rows[:2], "linear_tmb")
        first = {"group_mse": {"x": 0.1}, "predictions": {"a": 0.1}}
        second = {"group_mse": {"x": 0.1}, "predictions": {"b": 0.1}}
        with self.assertRaises(ValueError):
            analysis.paired_error_interval(first, second, replicates=10)
        missing = [dict(row, loh_fraction=None) if row["parent_type"] == "0" else row for row in rows]
        with self.assertRaises(ValueError):
            analysis.fit(missing, "power_clonal", clonal_with_hla=True)
        no_missing = [row for row in missing if row["parent_type"] != "0"]
        result = analysis.fit(no_missing, "power_clonal", clonal_with_hla=False)
        self.assertFalse(result["with_hla"])
        self.assertTrue(-2 <= result["theta"] <= 3)

    def test_all_hashes_are_newline_normalized(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory)/name for name in ("lf.csv", "crlf.csv", "cr.csv")]
            for path, data in zip(paths, (b"a,b\n1,2\n", b"a,b\r\n1,2\r\n", b"a,b\r1,2\r")):
                path.write_bytes(data)
            expected = hashlib.sha256(b"a,b\n1,2\n").hexdigest()
            self.assertEqual({analysis.digest(path) for path in paths}, {expected})


class TestPanCancerRecomputedEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.monotonic()
        cls.report = analysis.build_report()
        cls.rows = cls.report["rows"]
        cls.by_id = {row["id"]: row for row in cls.rows}
        cls.raw = {}
        for path in sorted((RESEARCH/"data").glob("*.csv")):
            with path.open(encoding="utf-8", newline="") as stream:
                cls.raw[path.name] = list(csv.DictReader(stream))

    def panel_rows(self, panel):
        return panel.get("rows") or [self.by_id[row_id] for row_id in panel["row_ids"]]

    def fitted_panels(self):
        return [panel for panel in self.report["panels"]+self.report["exploratory_panels"]
                if "models" in panel]

    def assert_metrics(self, expected, actual):
        for key in ("mse", "rmse", "mae", "mean_residual", "row_rmse"):
            self.assertAlmostEqual(expected[key], actual[key], places=11, msg=key)
        self.assertEqual(set(expected["group_mse"]), set(actual["group_mse"]))
        for group in expected["group_mse"]:
            self.assertAlmostEqual(expected["group_mse"][group], actual["group_mse"][group], places=11)

    def comparison_reference(self, rows, model):
        """Refit each training fold; independently form predictions and scores."""
        predictions = {model: {}, "log_tmb": {}}
        for omitted in {row["parent_type"] for row in rows}:
            training = [row for row in rows if row["parent_type"] != omitted]
            parameters = {model: analysis.fit(training, model),
                          "log_tmb": ols_reference(training, "log_tmb")}
            for row in rows:
                if row["parent_type"] == omitted:
                    for name in parameters:
                        predictions[name][row["id"]] = prediction_reference(parameters[name], row)
        metrics = {name: metric_reference(rows, [values[row["id"]] for row in rows])
                   for name, values in predictions.items()}
        parameters = analysis.fit(rows, model)
        return {"n_rows": len(rows), "n_groups": len({row["parent_type"] for row in rows}),
                "theory_rmse": metrics[model]["rmse"], "log_rmse": metrics["log_tmb"]["rmse"],
                "delta_mse": metrics[model]["mse"]-metrics["log_tmb"]["mse"],
                **{key: parameters[key] for key in ("A0", "A1", "B", "eta")}}

    def assert_comparison(self, rows, model, reported):
        expected = self.comparison_reference(rows, model)
        for field, value in expected.items():
            if value is None:
                self.assertIsNone(reported[field])
            else:
                self.assertAlmostEqual(value, reported[field], places=10, msg=field)

    def test_source_provenance_original_units_and_missingness(self):
        self.assertEqual(len(self.by_id), len(self.rows))
        for filename, rows in self.raw.items():
            self.assertTrue(rows, filename)
            self.assertEqual(len({row["id"] for row in rows}), len(rows), filename)
            for row in rows:
                self.assertTrue(row["id"])
                if filename == "tmb_orr.csv":
                    for prefix in ("orr_", "tmb_"):
                        if prefix == "tmb_" and not row.get("tmb"):
                            continue
                        for field in ("source_url", "source_title", "locator", "access_date", "reuse_basis", "terms_url"):
                            self.assertTrue(row.get(prefix+field, "").strip(), (row["id"], prefix+field))
                    self.assertEqual(Decimal(row["orr"]), Decimal(row["orr_original_value"])/100)
                    if row.get("tmb"):
                        self.assertEqual(Decimal(row["tmb"]), Decimal(row["tmb_original_value"]))
                    self.assertGreaterEqual(Decimal(row["orr_half_width"]), 0)
                else:
                    for field in ("source_url", "source_title", "access_date", "reuse_basis", "terms_url"):
                        self.assertTrue(row.get(field, "").strip(), (row["id"], field))
                    self.assertTrue(row.get("locator", row.get("source_locator", "")).strip())
                if filename == "clonal_fractions.csv":
                    self.assertEqual(row["measure_kind"], "neoantigen_clonal_fraction")
                    if row["original_value"] == "NA":
                        self.assertEqual(row["fraction"], "")
                    else:
                        self.assertEqual(Decimal(row["fraction"]), Decimal(row["original_value"])/100)
                        self.assertEqual(Decimal(row["original_value"])+Decimal(row["original_subclonal_percent"]), 100)
                        expected_half = Decimal(5).scaleb(Decimal(row["original_value"]).as_tuple().exponent-3)
                        self.assertEqual(Decimal(row["fraction_half_width"]), expected_half)
                    self.assertIn("not_patient_median", row["aggregation"])
                    self.assertIn("indexed", row["extraction_method"])
                if filename == "hla_loh.csv":
                    self.assertTrue(0 <= float(row["loh_fraction"]) <= 1)
                    self.assertGreaterEqual(float(row["loh_half_width"]), 0)
                    if row.get("pixel_bar_right"):
                        # Published aggregate bars, not individual patient points.
                        digitized = (float(row["pixel_bar_right"])-497.7)/728.5
                        self.assertAlmostEqual(float(row["loh_fraction"]), round(digitized, 2), places=12)
                    for pixel_field, fraction_field in (("pixel_ci_left", "source_ci95_low"),
                                                        ("pixel_ci_right", "source_ci95_high")):
                        if row.get(pixel_field):
                            digitized = (float(row[pixel_field])-497.7)/728.5
                            self.assertAlmostEqual(float(row[fraction_field]), round(bounded(digitized), 2), places=12)
                    if row.get("loh_numerator") and row.get("loh_denominator"):
                        exact = int(row["loh_numerator"])/int(row["loh_denominator"])
                        self.assertLessEqual(abs(float(row["loh_fraction"])-exact),
                                             float(row["loh_half_width"])+1e-10)
        paired = {row["id"]: row for row in self.raw["tmb_orr.csv"]}
        clonal = {row["id"]: row for row in self.raw["clonal_fractions.csv"]}
        hla = {row["id"]: row for row in self.raw["hla_loh.csv"]}
        self.assertEqual(set(self.by_id), {key for key, row in paired.items() if row["tmb"]})
        self.assertEqual({row["id"] for row in self.report["excluded_rows"]},
                         {key for key, row in paired.items() if not row["tmb"]})
        for row in self.rows:
            self.assertEqual(row["tmb"], float(paired[row["id"]]["tmb"]))
            self.assertEqual(row["orr"], float(paired[row["id"]]["orr"]))
            self.assertGreater(row["tmb"], 0)
            self.assertTrue(0 <= row["orr"] <= 1)
            self.assertTrue(row["parent_type"])
            for field in ("fraction", "loh_fraction"):
                if row.get(field) is not None:
                    self.assertTrue(0 <= row[field] <= 1, (row["id"], field))
            if row.get("fraction") is not None:
                self.assertEqual(row["fraction"], float(clonal[row["clonal_source_id"]]["fraction"]))
            if row.get("loh_fraction") is not None:
                self.assertEqual(row["loh_fraction"], float(hla[row["hla_source_id"]]["loh_fraction"]))

    def test_gates_use_distinct_parents_and_matched_subsets(self):
        expected_ids = {
            "all": {row["id"] for row in self.rows},
            "clonal": {row["id"] for row in self.rows if row.get("fraction") is not None},
            "hla": {row["id"] for row in self.rows if row.get("loh_fraction") is not None},
            "full": {row["id"] for row in self.rows if row.get("fraction") is not None and row.get("loh_fraction") is not None},
        }
        gates = {"all": 3, "clonal": 5, "hla": 10, "full": 10}
        for panel in self.report["panels"]:
            rows = self.panel_rows(panel)
            self.assertEqual(set(panel["row_ids"]), expected_ids[panel["name"]])
            groups = {row["parent_type"] for row in rows}
            self.assertEqual(panel["n_rows"], len(rows))
            self.assertEqual(panel["n_groups"], len(groups))
            self.assertEqual(panel["gate"], gates[panel["name"]])
            self.assertEqual(panel["eligible"], len(groups) >= panel["gate"])
            if panel["eligible"]:
                for result in panel["models"].values():
                    self.assertEqual(set(result["cv"]["predictions"]), expected_ids[panel["name"]])
                    self.assertEqual(set(result["cv"]["folds"]), groups)
            else:
                self.assertNotIn("models", panel)
        for panel in self.report["exploratory_panels"]:
            self.assertEqual(panel["status"], "exploratory")
            self.assertNotEqual(self.report["primary_panel"], panel["name"])
            rows = self.panel_rows(panel)
            self.assertEqual(panel["n_rows"], len(rows))
            self.assertEqual(panel["n_groups"], len({row["parent_type"] for row in rows}))
            for result in panel["models"].values():
                self.assertEqual(set(result["cv"]["predictions"]), set(panel["row_ids"]))

    def test_every_baseline_fit_heldout_prediction_and_score(self):
        for panel in self.fitted_panels():
            rows = self.panel_rows(panel)
            for model in ("mean", "linear_tmb", "log_tmb"):
                result = panel["models"][model]
                parameters = ols_reference(rows, model)
                for key in ("intercept", "slope"):
                    self.assertAlmostEqual(parameters[key], result["parameters"][key], places=10)
                fitted = [prediction_reference(parameters, row) for row in rows]
                loss = math.fsum((row["orr"]-value)**2 for row, value in zip(rows, fitted))
                tss = math.fsum((row["orr"]-average(row["orr"] for row in rows))**2 for row in rows)
                self.assertAlmostEqual(loss, result["training_sse"], places=10)
                if tss:
                    self.assertAlmostEqual(1-loss/tss, result["training_r2"], places=11)
                predictions = []
                for row in rows:
                    train = [other for other in rows if other["parent_type"] != row["parent_type"]]
                    parameters = ols_reference(train, model)
                    for key in ("intercept", "slope"):
                        self.assertAlmostEqual(parameters[key], result["cv"]["folds"][row["parent_type"]][key], places=10)
                    predicted = prediction_reference(parameters, row)
                    self.assertAlmostEqual(predicted, result["cv"]["predictions"][row["id"]], places=10)
                    predictions.append(predicted)
                self.assert_metrics(metric_reference(rows, predictions), result["cv"])

    def test_every_nonlinear_fit_and_fold_against_independent_wide_search(self):
        for panel in self.fitted_panels():
            rows = self.panel_rows(panel)
            results = list(panel["models"].items())
            if "sign_diagnostic" in panel:
                sign = panel["sign_diagnostic"]
                results.append((sign["parameters"]["model"], sign))
            for model, result in results:
                if model in ("mean", "linear_tmb", "log_tmb"):
                    continue
                all_fits = [(rows, result["parameters"])]
                all_fits.extend(([row for row in rows if row["parent_type"] != omitted], parameters)
                                for omitted, parameters in result["cv"]["folds"].items())
                for training, parameters in all_fits:
                    reference = reciprocal_reference(*nonlinear_inputs(training, model, parameters))
                    predictions = [prediction_reference(parameters, row) for row in training]
                    loss = math.fsum((row["orr"]-predicted)**2 for row, predicted in zip(training, predictions))
                    self.assertAlmostEqual(loss, parameters["sse"], places=11)
                    self.assertLessEqual(loss, reference["sse"]+2e-8, model)
                    self.assertLessEqual(reference["sse"], loss+2e-8, model)
                    self.assertTrue(0 <= parameters["A0"] <= 1 and 0 <= parameters["A1"] <= 1)
                    self.assertGreaterEqual(parameters["B"], 0)
                    if not model.startswith("free_"):
                        self.assertLessEqual(parameters["A1"], parameters["A0"]+1e-12)
                    if parameters["eta"] is not None:
                        self.assertAlmostEqual(parameters["eta"], 1-parameters["A1"]/parameters["A0"], places=12)
                    if "ceiling_difference" in parameters:
                        self.assertAlmostEqual(parameters["ceiling_difference"],
                                               parameters["A1"]-parameters["A0"], places=12)
                predicted = []
                for row in rows:
                    parameters = result["cv"]["folds"][row["parent_type"]]
                    value = prediction_reference(parameters, row)
                    self.assertAlmostEqual(value, result["cv"]["predictions"][row["id"]], places=12)
                    predicted.append(value)
                self.assert_metrics(metric_reference(rows, predicted), result["cv"])
                fitted = [prediction_reference(result["parameters"], row) for row in rows]
                sse = math.fsum((row["orr"]-value)**2 for row, value in zip(rows, fitted))
                tss = math.fsum((row["orr"]-average(row["orr"] for row in rows))**2 for row in rows)
                self.assertAlmostEqual(sse, result["training_sse"], places=12)
                if tss:
                    self.assertAlmostEqual(1-sse/tss, result["training_r2"], places=11)

    def test_dense_grid_check_summary_recomputes(self):
        differences = []
        for panel in self.fitted_panels():
            rows = self.panel_rows(panel)
            for model, output in panel["models"].items():
                if model.startswith("reciprocal"):
                    refined = analysis.fit(rows, model, dense=True)
                    differences.append(output["parameters"]["sse"]-refined["sse"])
        self.assertAlmostEqual(self.report["grid_check_max_difference"], max(differences), places=12)
        self.assertLessEqual(max(differences), 1e-9)

    def test_paired_loss_bootstrap_recomputed_from_parent_errors(self):
        for panel in self.fitted_panels():
            theory = panel["models"][panel["model"]]["cv"]
            groups = sorted(theory["group_mse"])
            for comparator, comparison in panel["comparisons"].items():
                baseline = panel["models"][comparator]["cv"]
                differences = [theory["group_mse"][group]-baseline["group_mse"][group]
                               for group in groups]
                generator = random.Random(analysis.SEED)
                resampled = [average(generator.choices(differences, k=len(groups)))
                             for _ in range(comparison["replicates"])]
                self.assertEqual(comparison["replicates"], 10000)
                self.assertEqual(comparison["groups"], len(groups))
                self.assertAlmostEqual(comparison["delta_mse"], average(differences), places=13)
                for probability, actual in zip((0.025, 0.975), comparison["interval"]):
                    self.assertAlmostEqual(percentile(resampled, probability), actual, places=12)

    def test_parameter_bootstrap_recomputed_and_degenerate_draws_reported(self):
        for panel in self.fitted_panels():
            rows = self.panel_rows(panel)
            specifications = [("parameter_uncertainty", panel["model"])]
            if "sign_uncertainty" in panel:
                specifications.append(("sign_uncertainty", panel["sign_diagnostic"]["parameters"]["model"]))
            for field, model in specifications:
                reported = panel[field]
                groups = sorted({row["parent_type"] for row in rows})
                generator = random.Random(analysis.SEED)
                values, skipped, boundaries = {}, 0, 0
                for _ in range(reported["replicates"]):
                    chosen = generator.choices(groups, k=len(groups))
                    if len(set(chosen)) < 3:
                        skipped += 1
                        continue
                    sample = [row for group in chosen for row in rows if row["parent_type"] == group]
                    parameters = analysis.fit(sample, model)
                    if parameters.get("singular", False):
                        skipped += 1
                        continue
                    for key in ("A0", "A1", "B", "eta", "ceiling_difference", "intercept", "slope"):
                        if parameters.get(key) is not None:
                            values.setdefault(key, []).append(parameters[key])
                    # Distance to each feasible boundary; count a fit once.
                    distances = [parameters.get("B", 1), parameters.get("A0", .5),
                                 1-parameters.get("A0", .5)]
                    if model in ("reciprocal_hla", "reciprocal_full", "free_hla", "free_full"):
                        distances.extend([parameters["A1"], 1-parameters["A1"]])
                        if not model.startswith("free_"):
                            distances.append(abs(parameters["A0"]-parameters["A1"]))
                    boundaries += min(distances) < 1e-8
                self.assertEqual(reported["replicates"], 1000)
                self.assertEqual(reported["skipped"], skipped)
                self.assertEqual(reported["boundary_fits"], boundaries)
                self.assertEqual(reported["usable"], {key: len(value) for key, value in values.items()})
                for key, value in values.items():
                    expected = [percentile(value, probability) for probability in (0.025, 0.975)]
                    for first, second in zip(expected, reported["intervals"][key]):
                        self.assertAlmostEqual(first, second, places=10)

    def test_exponents_are_secondary_constrained_profiles(self):
        def check_profile(rows, parameters, field, lower, upper):
            # Refit the secondary coefficient, then compare its conditional
            # minimum with an independent projected search, boundary profiles
            # and both sides of the reported coefficient.
            model = parameters["model"]
            options = {"clonal_with_hla": parameters["with_hla"]} if model == "power_clonal" else {}
            repeated = analysis.fit(rows, model, **options)
            for key in ("A0", "A1", "B", "sse", field):
                self.assertAlmostEqual(repeated[key], parameters[key], places=10, msg=key)
            reference = reciprocal_reference(*nonlinear_inputs(rows, model, parameters))
            self.assertAlmostEqual(reference["sse"], parameters["sse"], places=8)
            for candidate in {lower, upper, 1.0, bounded(parameters[field]-.01, lower, upper),
                              bounded(parameters[field]+.01, lower, upper)}:
                alternative = dict(parameters, **{field: candidate})
                reference = reciprocal_reference(*nonlinear_inputs(rows, model, alternative))
                self.assertLessEqual(parameters["sse"], reference["sse"]+2e-8)

        diagnostic = self.report["power_diagnostic"]
        parameters = diagnostic["parameters"]
        self.assertTrue(0.05 <= parameters["q"] <= 4)
        check_profile(self.rows, parameters, "q", .05, 4.0)
        self.assertEqual(set(diagnostic["jackknife_q"]), {row["parent_type"] for row in self.rows})
        self.assertTrue(all(0.05 <= value <= 4 for value in diagnostic["jackknife_q"].values()))
        for omitted, value in diagnostic["jackknife_q"].items():
            remaining = [row for row in self.rows if row["parent_type"] != omitted]
            self.assertAlmostEqual(value, analysis.fit(remaining, "power_tmb")["q"], places=10)
        fixed = analysis.fit(self.rows, "reciprocal_tmb")
        self.assertLessEqual(parameters["sse"], fixed["sse"]+2e-8)
        if "clonal_power_diagnostic" in self.report:
            clonal = [row for row in self.rows if row.get("fraction") is not None]
            diagnostic = self.report["clonal_power_diagnostic"]
            parameters = diagnostic["parameters"]
            self.assertTrue(-2 <= parameters["theta"] <= 3)
            check_profile(clonal, parameters, "theta", -2.0, 3.0)
            self.assertEqual(set(diagnostic["jackknife_theta"]), {row["parent_type"] for row in clonal})
            self.assertTrue(all(-2 <= value <= 3 for value in diagnostic["jackknife_theta"].values()))
            for omitted, value in diagnostic["jackknife_theta"].items():
                remaining = [row for row in clonal if row["parent_type"] != omitted]
                repeated = analysis.fit(remaining, "power_clonal", clonal_with_hla=parameters["with_hla"])
                self.assertAlmostEqual(value, repeated["theta"], places=10)
            fixed_model = "reciprocal_full" if parameters["with_hla"] else "reciprocal_clonal"
            self.assertLessEqual(parameters["sse"], analysis.fit(clonal, fixed_model)["sse"]+2e-8)

        for panel in self.report["exploratory_panels"]:
            if "clonal_power" in panel:
                parameters = panel["clonal_power"]
                self.assertTrue(-2 <= parameters["theta"] <= 3)
                self.assertTrue(parameters["with_hla"])
                check_profile(self.panel_rows(panel), parameters, "theta", -2.0, 3.0)
                self.assertLessEqual(parameters["sse"], panel["models"]["reciprocal_full"]["parameters"]["sse"]+2e-8)
                self.assertEqual(set(panel["clonal_power_jackknife"]),
                                 {row["parent_type"] for row in self.panel_rows(panel)})
                for omitted, value in panel["clonal_power_jackknife"].items():
                    remaining = [row for row in self.panel_rows(panel) if row["parent_type"] != omitted]
                    self.assertAlmostEqual(value, analysis.fit(remaining, "power_clonal", clonal_with_hla=True)["theta"], places=10)

    def test_sensitivity_values_recomputed_with_fixed_mapping_and_perturbations(self):
        primary = next(panel for panel in self.report["panels"]
                       if panel["name"] == self.report["primary_panel"])
        full = [row for row in self.rows
                if row.get("fraction") is not None and row.get("loh_fraction") is not None]
        specifications = [("primary", self.panel_rows(primary), primary["model"])]
        if "coherent_full" in self.report["sensitivities"]:
            specifications.append(("coherent_full", full, "reciprocal_full"))
        for name, rows, model in specifications:
            reported = self.report["sensitivities"][name]
            exact = [row for row in rows if row["mapping_quality"] == "exact"
                     and (model not in ("reciprocal_full", "reciprocal_clonal")
                          or row["clonal_mapping_quality"] == "exact")
                     and (model not in ("reciprocal_full", "reciprocal_hla")
                          or row["hla_mapping_quality"] == "exact")]
            if len({row["parent_type"] for row in exact}) >= 3:
                self.assert_comparison(exact, model, reported["exact_mapping"])
            else:
                self.assertTrue(reported["exact_mapping"]["unavailable"])
                self.assertEqual(reported["exact_mapping"]["n_rows"], len(exact))
            groups = {row["parent_type"] for row in rows}
            self.assertEqual(set(reported["single_parent_deletions"]), groups)
            for omitted, comparison in reported["single_parent_deletions"].items():
                remaining = [row for row in rows if row["parent_type"] != omitted]
                self.assert_comparison(remaining, model, comparison)
            fields = [("tmb", "tmb_half_width"), ("orr", "orr_half_width")]
            if model in ("reciprocal_full", "reciprocal_clonal"):
                fields.append(("fraction", "fraction_half_width"))
            if model in ("reciprocal_full", "reciprocal_hla"):
                fields.append(("loh_fraction", "loh_half_width"))
            generator = random.Random(analysis.SEED+1)
            differences = []
            for _ in range(reported["rounding_envelope"]["draws"]):
                changed = []
                for row in rows:
                    copy = dict(row)
                    for field, width in fields:
                        copy[field] = row[field] + (-1+2*generator.random())*row[width]
                        if field != "tmb":
                            copy[field] = bounded(copy[field], 1e-12 if field == "fraction" else 0.0)
                    changed.append(copy)
                differences.append(self.comparison_reference(changed, model)["delta_mse"])
            for expected, actual in zip((min(differences), max(differences)),
                                        reported["rounding_envelope"]["delta_mse_range"]):
                self.assertAlmostEqual(expected, actual, places=10)
            selection = reported["biomarker_selection"]
            self.assertEqual(selection["known_selected"], sum(row["biomarker_selected"] == "yes" for row in rows))
            self.assertEqual(selection["unknown"], sum(row["biomarker_selected"] == "unknown" for row in rows))
            if "hla_source_interval_scenarios" in reported:
                for label, field in (("all_lower", "source_ci95_low"), ("all_upper", "source_ci95_high")):
                    changed = [dict(row, loh_fraction=row[field] if row.get(field) is not None else row["loh_fraction"])
                               for row in rows]
                    self.assert_comparison(changed, model, reported["hla_source_interval_scenarios"][label])
        if "figure_instead_of_printed" in self.report["sensitivities"]:
            rows, _, _ = analysis.load_rows(figure_values=True)
            rows = [row for row in rows if row["fraction"] is not None and row["loh_fraction"] is not None]
            self.assert_comparison(rows, "reciprocal_full", self.report["sensitivities"]["figure_instead_of_printed"])
        eligible_without_indexed = [row for row in self.raw["clonal_fractions.csv"]
                                   if row["fraction"] and "indexed" not in row["extraction_method"]]
        self.assertEqual(self.report["source_quality"]["clonal_without_indexed_only"], len(eligible_without_indexed))

    def test_readme_every_generated_number_hash_and_json_recomputes(self):
        text = (RESEARCH/"README.md").read_text(encoding="utf-8")
        self.assertEqual(text.count(analysis.BEGIN), 1)
        self.assertEqual(text.count(analysis.END), 1)
        actual = text.split(analysis.BEGIN, 1)[1].split(analysis.END, 1)[0]
        self.assertEqual(actual, "\n\n"+analysis.render_report(self.report)+"\n\n")
        json.dumps(self.report, allow_nan=False)
        for filename, checksum in self.report["hashes"].items():
            path = RESEARCH/filename
            if not path.exists():
                path = RESEARCH/"data"/filename
            contents = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
            self.assertEqual(checksum, hashlib.sha256(contents).hexdigest(), filename)
        protocol = (RESEARCH/"PROTOCOL.md").read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
        self.assertEqual(self.report["protocol_sha256"], hashlib.sha256(protocol).hexdigest())
        primary = next(panel for panel in self.report["panels"]
                       if panel["name"] == self.report["primary_panel"])
        for model in (primary["model"], "log_tmb"):
            predicted = primary["models"][model]["cv"]["predictions"]
            parents = {row["parent_type"] for row in self.rows}
            expected = sorted(parents, key=lambda parent: (
                -average(abs(row["orr"]-predicted[row["id"]])
                         for row in self.rows if row["parent_type"] == parent), parent))[:5]
            self.assertEqual([parent for parent, _ in analysis.ranked_outliers(self.rows, predicted)], expected)
        self.assertIn("nonclinical", text.lower())
        self.assertIn("ecological", text.lower())

    def test_zz_dedicated_recomputation_finishes_within_one_minute(self):
        self.assertLess(time.monotonic()-self.started, 60, "Dedicated recomputation must finish within one minute")


if __name__ == "__main__":
    unittest.main()
