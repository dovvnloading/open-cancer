"""Reproduce the Peto aggregate note offline; stdlib only and nonclinical.

These checks establish reproducibility and consistency, not biological validity.
Independent calculations use the transcribed CSVs, pair-difference regression,
the binomial score equation, and Decimal likelihoods. No medical advice.
"""

import csv
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path, PureWindowsPath
import sys
import tempfile
import unittest
from urllib.parse import urlsplit


RESEARCH = Path(__file__).resolve().parents[1] / "research" / "petos-paradox-mutation-rates"
_spec = importlib.util.spec_from_file_location("research_petos_mutation_analysis", RESEARCH / "analysis.py")
A = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = A
_spec.loader.exec_module(A)

KS = tuple(range(1, 9))
PANELS = ("primary", "maximum_longevity", "slope_estimator")
SCENARIOS = ("point", "cmr95", "sex_mass", "slope_estimator", "slope_ci", "factor_two")
PROVENANCE = ("source", "source_url", "source_table", "locator", "access_date",
              "reuse_basis", "reuse_terms_url")
FREEZE = "b98d28e5258c09132d197e26529ff4f58d18b8a9"


def source_rows(filename):
    with (RESEARCH / "data" / filename).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def score_interval(events, total):
    """Solve n*(observed-p)**2 = z**2*p*(1-p), independently in Decimal."""
    with localcontext() as context:
        context.prec = 65
        n, e = Decimal(total), Decimal(events)
        z = Decimal("1.959963984540054")
        aa, bb, cc = n + z * z, -(2 * e + z * z), e * e / n
        root = (bb * bb - 4 * aa * cc).sqrt()
        return float((-bb - root) / (2 * aa)), float((-bb + root) / (2 * aa))


def pairwise_line(xs, ys):
    """OLS from all pair differences, not centered normal equations."""
    pairs = list(itertools.combinations(range(len(xs)), 2))
    denominator = math.fsum((xs[j] - xs[i]) ** 2 for i, j in pairs)
    if denominator == 0:
        raise ValueError("rank-deficient independent design")
    slope = math.fsum((xs[j] - xs[i]) * (ys[j] - ys[i]) for i, j in pairs) / denominator
    intercept = math.fsum(y - slope * x for x, y in zip(xs, ys)) / len(xs)
    return intercept, slope


def model_specs():
    yield "intercept", "intercept", 0
    yield "mass", "mass", 0
    for family in ("burden", "cell", "free_mass"):
        for k in KS:
            yield f"{family}_k{k}", family, k


def independent_parameters(rows, family, k):
    exponent = 0 if family in ("intercept", "mass") else k
    y = [math.log(r["H"]) - exponent * math.log(r["B"]) for r in rows]
    if family == "free_mass":
        intercept, mass = pairwise_line([math.log(r["M"]) for r in rows], y)
    else:
        mass = 1.0 if family in ("cell", "mass") else 0.0
        intercept = math.fsum(v - mass * math.log(r["M"]) for v, r in zip(y, rows)) / len(rows)
    return {"intercept": intercept, "mass_exponent": mass, "burden_exponent": exponent}


def decimal_prediction(log_h, risk):
    """Use direct high-precision survival; no analysis probability/NLL helpers."""
    with localcontext() as context:
        context.prec = 85
        h = Decimal(str(log_h)).exp()
        observed = Decimal(str(risk))
        # For h>1000 the omitted log-risk is below 1e-434, far below
        # double precision at every score in this calculation.
        p = Decimal(1) if h > 1000 else 1 - (-h).exp()
        loss = -observed * p.ln() + (1 - observed) * h
        return float(p), float(loss)


def independent_cv(rows, family, k, grouping="species"):
    predictions = []
    for group in sorted({r[grouping] for r in rows}):
        train = [r for r in rows if r[grouping] != group]
        if family == "free_mass":
            parameters = independent_parameters(train, family, k)
        else:
            mass = 1 if family in ("mass", "cell") else 0
            exponent = 0 if family in ("intercept", "mass") else k
            # Closed-form deleted-group intercept: total adjusted response
            # minus held-out adjusted response, divided by training size.
            adjusted = [(r, math.log(r["H"]) - mass * math.log(r["M"])
                         - exponent * math.log(r["B"])) for r in rows]
            total = math.fsum(v for _, v in adjusted)
            removed = math.fsum(v for r, v in adjusted if r[grouping] == group)
            parameters = {"intercept": (total - removed) / len(train),
                          "mass_exponent": mass, "burden_exponent": exponent}
        for row in rows:
            if row[grouping] != group:
                continue
            log_h = (parameters["intercept"] + parameters["mass_exponent"] * math.log(row["M"])
                     + parameters["burden_exponent"] * math.log(row["B"]))
            residual = math.log10(row["H"]) - log_h / math.log(10)
            predicted_risk, nll = decimal_prediction(log_h, row["R"])
            predictions.append({"species": row["species"], "held_out_group": group,
                                "predicted_log10_h": log_h / math.log(10),
                                "predicted_risk": predicted_risk,
                                "residual_log10_h": residual, "squared_error": residual ** 2,
                                "nll": nll, "training_mass_exponent": parameters["mass_exponent"]})
    predictions.sort(key=lambda r: r["species"])
    return {"rmse_log10_h": math.sqrt(math.fsum(p["squared_error"] for p in predictions) / len(predictions)),
            "mean_nll": math.fsum(p["nll"] for p in predictions) / len(predictions),
            "predictions": predictions}


def independent_box(row, scenario):
    if scenario == "point":
        hazards = (row["H"], row["H"])
    else:
        hazards = tuple(row["H_box"])
    masses = tuple(row["M_box"]) if scenario in ("sex_mass", "factor_two") else (row["M"], row["M"])
    if scenario == "factor_two":
        burdens = (row["m"] * row["L"] / 4, row["m"] * row["L"] * 4)
    elif scenario in ("slope_estimator", "slope_ci"):
        rates = (row["slope"], row["slope"])
        if scenario == "slope_ci" and row["slope_ci"] is not None:
            rates = (max(0, row["slope_ci"][0]), row["slope_ci"][1])
        burdens = (rates[0] * row["L"], rates[1] * row["L"])
    else:
        burdens = (row["m"] * row["L"], row["m"] * row["L"])
    return {"H": hazards, "M": masses, "B": burdens}


def independent_pair(first, second, k, scenario):
    a, b = independent_box(first, scenario), independent_box(second, scenario)
    yl = (a["H"][0] / b["H"][1]) / (a["M"][1] / b["M"][0])
    yu = (a["H"][1] / b["H"][0]) / (a["M"][0] / b["M"][1])
    bl = a["B"][0] / b["B"][1]
    bu = a["B"][1] / b["B"][0] if b["B"][0] else math.inf
    ql, qu = min(bl, bl ** k), max(bu, bu ** k)
    low, high = yl / qu, yu / ql if ql else math.inf
    fold = max(1.0, low, 1 / high)
    return {"y_interval": [yl, yu], "burden_ratio_interval": [bl, bu],
            "multiplier_low": low, "multiplier_high": high,
            "fold_infimum": fold, "rho_infimum": math.sqrt(fold)}


class PetosParadoxMutationRatesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mutation = {r["species"]: r for r in source_rows("cagan_mutation_rates.csv")}
        cls.traits = {r["species"]: r for r in source_rows("compton_life_history.csv")}
        cls.mortality = source_rows("vincze_mortality.csv")
        cls.rows = []
        for source in cls.mortality:
            species = source["species"]
            if species not in cls.mutation or species not in cls.traits:
                continue
            mutation, traits = cls.mutation[species], cls.traits[species]
            e, n = int(source["cancer_deaths"]), int(source["known_cause_deaths"])
            r = e / n
            lower, upper = score_interval(e, n)
            female, male = float(source["female_mass_kg"]), float(source["male_mass_kg"])
            row = {"species": species, "order": source["order"],
                   "m": float(mutation["mean_sbs_per_genome_year"]),
                   "slope": float(mutation["zero_intercept_slope_sbs_per_genome_year"]),
                   "slope_ci": ([float(mutation["slope_ci95_lower"]), float(mutation["slope_ci95_upper"])]
                                if mutation["slope_ci95_lower"] else None),
                   "M": (female + male) / 2, "M_box": sorted([female, male]),
                   "L": float(source["adult_life_expectancy_days"]) / 365.25,
                   "L_max": float(traits["max_longevity_months"]) / 12,
                   "events": e, "n": n, "total_animals": int(source["total_animals"]),
                   "total_deaths": int(source["total_deaths"]), "R": r, "R_box": [lower, upper],
                   "H": -math.log1p(-r), "H_box": [-math.log1p(-lower), -math.log1p(-upper)],
                   "ICM": float(source["icm"]) if source["icm"] else None}
            cls.rows.append(row)
        cls.rows.sort(key=lambda r: r["species"])
        cls.panels = {}
        for panel in PANELS:
            rows = [dict(r) for r in cls.rows]
            for row in rows:
                if panel == "maximum_longevity":
                    row["L"] = row["L_max"]
                if panel == "slope_estimator":
                    row["m"] = row["slope"]
                row["B"] = row["m"] * row["L"]
            cls.panels[panel] = rows
        cls.results = A.analyze()

    def assert_close(self, actual, expected, message="", rel=3e-10, absolute=0.0):
        if math.isinf(expected):
            self.assertEqual(actual, expected, message)
            return
        self.assertTrue(math.isfinite(actual), message)
        self.assertTrue(math.isclose(actual, expected, rel_tol=rel,
                                     abs_tol=absolute if expected else max(absolute, 2e-13)),
                        message or f"{actual!r} != {expected!r}")

    def assert_tree(self, actual, expected, path="results"):
        if isinstance(expected, dict):
            self.assertIsInstance(actual, dict, path)
            self.assertEqual(set(actual), set(expected), path)
            for key in expected:
                self.assert_tree(actual[key], expected[key], f"{path}.{key}")
        elif isinstance(expected, (list, tuple)):
            self.assertIsInstance(actual, (list, tuple), path)
            self.assertEqual(len(actual), len(expected), path)
            for i, (left, right) in enumerate(zip(actual, expected)):
                self.assert_tree(left, right, f"{path}[{i}]")
        elif isinstance(expected, float):
            self.assertIsInstance(actual, (int, float), path)
            self.assert_close(actual, expected, path)
        else:
            self.assertEqual(actual, expected, path)

    def test_01_complete_saved_report_and_generated_readme(self):
        def reject_constant(value):
            raise AssertionError(f"non-JSON numeric constant: {value}")

        saved = json.loads((RESEARCH / "analysis_results.json").read_text(encoding="utf-8"),
                           parse_constant=reject_constant)
        self.assert_tree(saved, A.serializable(self.results))
        self.assertEqual(saved["protocol_freeze_commit"], FREEZE)
        self.assertEqual(saved["biological_verdict"], "not_identifiable")
        self.assertEqual(saved["status"], "nonclinical_unreviewed")
        text = (RESEARCH / "README.md").read_text(encoding="utf-8")
        begin, end = "<!-- BEGIN COMPUTED RESULTS -->", "<!-- END COMPUTED RESULTS -->"
        self.assertEqual(text.count(begin), 1)
        self.assertEqual(text.count(end), 1)
        block = text[text.index(begin):text.index(end) + len(end)]
        self.assertEqual(block, A.render_results(self.results))
        prose = " ".join(text.lower().split())
        self.assertIn("nonclinical", prose)
        self.assertIn("no medical advice", prose)
        self.assertIn(FREEZE, text)

    def test_02_hashes_are_independent_of_checkout_newlines_and_platform_paths(self):
        paths = sorted((RESEARCH / "data").glob("*.csv")) + [RESEARCH / "PROTOCOL.md"]
        expected, contents = {}, {}
        for path in paths:
            key = path.relative_to(RESEARCH).as_posix()
            # Universal-newline text reading is independent of the analysis's
            # byte-replacement implementation and preserves a final newline.
            with path.open(encoding="utf-8", newline=None) as stream:
                canonical = stream.read().encode("utf-8")
            contents[key] = canonical
            expected[key] = hashlib.sha256(canonical).hexdigest()
            self.assertNotIn("\\", key)
            self.assertEqual(PureWindowsPath(key).as_posix(), key)
        self.assertEqual(expected, self.results["input_sha256"])
        self.assertEqual(expected, A.input_hashes())
        for newline in (b"\n", b"\r\n", b"\r"):
            with tempfile.TemporaryDirectory() as folder:
                base = Path(folder)
                for key, content in contents.items():
                    target = base / key
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(content.replace(b"\n", newline))
                self.assertEqual(A.input_hashes(base), expected)
                self.assert_tree(A.serializable(A.analyze(base)), A.serializable(self.results))

    def test_03_source_metadata_missingness_and_join_contract(self):
        for path in sorted((RESEARCH / "data").glob("*.csv")):
            rows = source_rows(path.name)
            self.assertTrue(rows, path.name)
            self.assertEqual(len(rows), len({r["species"] for r in rows}), path.name)
            for row in rows:
                self.assertNotIn(None, row, path.name)
                for field in PROVENANCE:
                    self.assertTrue(row[field].strip(), (path.name, row["species"], field))
                for field in ("source_url", "reuse_terms_url"):
                    parsed = urlsplit(row[field])
                    self.assertEqual(parsed.scheme, "https")
                    self.assertTrue(parsed.netloc)
                self.assertRegex(row["access_date"], r"^\d{4}-\d{2}-\d{2}$")
                self.assertIn("by/4.0", row["reuse_terms_url"])
                if row.get("source_commit"):
                    self.assertRegex(row["source_commit"], r"^[0-9a-f]{40}$")
        mutation_rows = list(self.mutation.values())
        for row in mutation_rows:
            self.assertEqual(row["mutation_individuals_n"], "")
            self.assertEqual(row["mutation_crypts_n"], "")
            self.assertIn("Not reported", row["counts_status"])
            self.assertGreater(float(row["mean_sbs_per_genome_year"]), 0)
            self.assertGreater(float(row["zero_intercept_slope_sbs_per_genome_year"]), 0)
            self.assertEqual(bool(row["slope_ci95_lower"]), bool(row["slope_ci95_upper"]))
            self.assertNotIn("mean_ci95_lower", row)
            self.assertNotIn("mean_ci95_upper", row)
        self.assertLess(float(self.mutation["Panthera_tigris"]["slope_ci95_lower"]), 0)
        self.assertEqual(self.mutation["Giraffa_camelopardalis"]["slope_ci95_lower"], "")
        self.assertEqual(self.mutation["Heterocephalus_glaber"]["source_common_name"], "Naked mole rat")
        expected = set(self.mutation) & set(self.traits) & {r["species"] for r in self.mortality}
        self.assertEqual({r["species"] for r in self.results["observations"]}, expected)
        self.assert_tree(self.results["observations"], self.panels["primary"])
        loaded, _, _ = A.load_data()
        for name in PANELS:
            self.assert_tree(A.panel_version(loaded, name), self.panels[name])
        with self.assertRaises(ValueError):
            A.panel_version(loaded, "unregistered_panel")

    def test_04_death_denominators_score_intervals_and_missing_icm(self):
        for source in self.mortality:
            events, known, deaths, animals = (int(source[k]) for k in
                                            ("cancer_deaths", "known_cause_deaths", "total_deaths", "total_animals"))
            self.assertTrue(0 < events <= known <= deaths <= animals)
            self.assert_close(float(source["cmr"]), events / known, absolute=2e-14)
            self.assertNotEqual(known, deaths)
            self.assertNotEqual(known, animals)
            self.assertEqual(source["published_cmr_ci_lower"], "")
            self.assertEqual(source["published_cmr_ci_upper"], "")
            self.assertEqual(source["ci_status"], "absent_in_source_table")
            low, high = A.wilson(events, known)
            self.assert_tree((low, high), score_interval(events, known))
            for p in (low, high):
                score = known * (events / known - p) ** 2 / (p * (1 - p))
                self.assert_close(score, 1.959963984540054 ** 2)
        self.assertIsNone(next(r for r in self.rows if r["species"] == "Heterocephalus_glaber")["ICM"])
        for events, total in ((0, 1), (1, 1), (0, 97), (97, 97), (3, 100)):
            self.assert_tree(A.wilson(events, total), score_interval(events, total))
        for events, total in ((-1, 5), (6, 5), (1, 0), (0.5, 10)):
            with self.assertRaises(ValueError):
                A.wilson(events, total)
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            (base / "data").mkdir()
            for path in (RESEARCH / "data").glob("*.csv"):
                (base / "data" / path.name).write_bytes(path.read_bytes())
            rows = [dict(r) for r in self.mortality]
            rows[0]["cmr"] = str(int(rows[0]["cancer_deaths"]) / int(rows[0]["total_animals"]))
            with (base / "data" / "vincze_mortality.csv").open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "denominator"):
                A.load_data(base)

    def test_05_every_family_fixed_k_panel_and_fold_is_independently_recomputed(self):
        self.assertEqual(tuple(A.K_VALUES), KS)
        self.assertEqual(set(self.results["panels"]), set(PANELS))
        for panel_name, rows in self.panels.items():
            panel = self.results["panels"][panel_name]
            self.assertEqual(set(panel), {name for name, _, _ in model_specs()})
            for name, family, k in model_specs():
                with self.subTest(panel=panel_name, model=name):
                    self.assert_tree(panel[name]["parameters"], independent_parameters(rows, family, k))
                    self.assert_tree(panel[name]["loso"], independent_cv(rows, family, k))
                    if panel_name == "primary":
                        self.assert_tree(panel[name]["leave_order_out"], independent_cv(rows, family, k, "order"))
                    else:
                        self.assertNotIn("leave_order_out", panel[name])

    def test_06_score_difference_ranges_and_outliers_use_every_held_out_fold(self):
        for panel_name, panel in self.results["panels"].items():
            reference = {p["species"]: p for p in independent_cv(self.panels[panel_name], "intercept", 0)["predictions"]}
            for name, family, k in model_specs():
                predictions = independent_cv(self.panels[panel_name], family, k)["predictions"]
                for score in ("nll", "squared_error"):
                    differences = [p[score] - reference[p["species"]][score] for p in predictions]
                    total, n = math.fsum(differences), len(differences)
                    deleted = [(total - value) / (n - 1) for value in differences]
                    expected = {"mean": total / n, "delete_fold_min": min(deleted), "delete_fold_max": max(deleted)}
                    self.assert_tree(panel[name]["vs_intercept"][score], expected)
        cv = independent_cv(self.panels["primary"], "cell", 3)["predictions"]
        outliers = [{"species": p["species"], "residual_log10_h": p["residual_log10_h"],
                     "observed_over_predicted_h": math.exp(math.log(10) * p["residual_log10_h"])} for p in cv]
        outliers.sort(key=lambda p: (-p["residual_log10_h"], p["species"]))
        self.assert_tree(self.results["outliers_k3"], outliers)

    def test_07_held_out_outcomes_cannot_leak_into_predictions_or_order_folds(self):
        rows = self.panels["primary"]
        groups = {order: {r["species"] for r in rows if r["order"] == order} for order in {r["order"] for r in rows}}
        self.assertEqual(groups["Primates"], {"Colobus_guereza", "Lemur_catta"})
        self.assertEqual(groups["Carnivora"], {"Panthera_leo", "Panthera_tigris"})
        for grouping in ("species", "order"):
            for name, family, k in model_specs():
                original = {p["species"]: p for p in A.cross_validate(rows, family, k, grouping)["predictions"]}
                for group in {r[grouping] for r in rows}:
                    changed = [dict(r) for r in rows]
                    held = set()
                    for i, row in enumerate(changed):
                        if row[grouping] == group:
                            held.add(row["species"])
                            row["R"] = 0.45 + 0.01 * i
                            row["H"] = -math.log1p(-row["R"])
                    predicted = {p["species"]: p for p in A.cross_validate(changed, family, k, grouping)["predictions"]}
                    for species in held:
                        for field in ("predicted_log10_h", "predicted_risk", "training_mass_exponent"):
                            self.assert_close(predicted[species][field], original[species][field], f"{grouping}/{name}/{species}/{field}")
                        self.assertEqual(predicted[species]["held_out_group"], group)
                        self.assertNotEqual(predicted[species]["residual_log10_h"], original[species]["residual_log10_h"])

    def test_08_all_pair_envelopes_multiplier_directions_and_worst_pairs(self):
        rows = self.panels["primary"]
        expected_pairs = {}
        for first, second in itertools.combinations(rows, 2):
            first, second = sorted((first, second), key=lambda r: (-r["M"], r["species"]))
            expected_pairs[(first["species"], second["species"])] = (first, second)
        actual = {(p["larger"], p["smaller"]): p for p in self.results["exact_pairs"]}
        self.assertEqual(set(actual), set(expected_pairs))
        for key, (first, second) in expected_pairs.items():
            self.assertEqual(set(actual[key]["scenarios"]), set(SCENARIOS))
            for scenario in SCENARIOS:
                self.assertEqual(set(actual[key]["scenarios"][scenario]), {str(k) for k in KS})
                for k in KS:
                    expected = independent_pair(first, second, k, scenario)
                    self.assert_tree(actual[key]["scenarios"][scenario][str(k)], expected)
                    reverse = A.pair_certificate(second, first, k, scenario)
                    self.assert_close(reverse["fold_infimum"], expected["fold_infimum"])
                    self.assert_close(reverse["multiplier_low"], 1 / expected["multiplier_high"])
                    self.assert_close(reverse["multiplier_high"], 1 / expected["multiplier_low"] if expected["multiplier_low"] else math.inf)
        main = actual[("Giraffa_camelopardalis", "Heterocephalus_glaber")]
        self.assert_tree(self.results["main_pair"], main)
        for k in KS:
            choices = [{"larger": key[0], "smaller": key[1],
                        "fold_infimum": independent_pair(*pair, k, "sex_mass")["fold_infimum"]}
                       for key, pair in expected_pairs.items()]
            expected = max(choices, key=lambda p: (p["fold_infimum"], p["larger"], p["smaller"]))
            self.assert_tree(self.results["exact_worst_pair_by_k"][str(k)], expected)

    def test_09_sensitivity_boxes_preserve_estimators_mass_and_factor_two_scope(self):
        for row in self.panels["primary"]:
            for scenario in SCENARIOS:
                self.assert_tree(A.box(row, scenario), independent_box(row, scenario))
            self.assert_tree(A.box(row, "factor_two")["B"], (row["B"] / 4, row["B"] * 4))
            for scenario in ("slope_estimator", "slope_ci"):
                self.assertEqual(A.box(row, scenario)["M"], (row["M"], row["M"]))
            self.assertEqual(A.box(row, "sex_mass")["M"], tuple(row["M_box"]))
            if row["slope_ci"] is None:
                self.assertEqual(A.box(row, "slope_ci")["B"], (row["slope"] * row["L"],) * 2)
            else:
                source = self.mutation[row["species"]]
                self.assertEqual(A.box(row, "slope_ci")["B"][1], float(source["slope_ci95_upper"]) * row["L"])
                if row["slope_ci"][0] < 0:
                    self.assertEqual(A.box(row, "slope_ci")["B"][0], 0)
                    self.assertGreater(row["m"], 0)
                    self.assertGreater(A.box(row, "cmr95")["B"][0], 0)
        for pair in self.results["exact_pairs"]:
            point = pair["scenarios"]["point"]["1"]["burden_ratio_interval"][0]
            stress = pair["scenarios"]["factor_two"]["1"]["burden_ratio_interval"]
            self.assert_close(stress[0], point / 16)
            self.assert_close(stress[1], point * 16)
        for k in KS:
            for scenario in ("point", "cmr95", "factor_two"):
                intervals = []
                for row in self.panels["primary"]:
                    box = independent_box(row, scenario)
                    lower = box["H"][0] / (box["M"][1] * box["B"][1] ** k)
                    upper = box["H"][1] / (box["M"][0] * box["B"][0] ** k)
                    intervals.append((row["species"], lower, upper))
                low = min(intervals, key=lambda v: (-v[1], v[0]))
                high = min(intervals, key=lambda v: (v[2], v[0]))
                rho = max(1.0, math.sqrt(low[1] / high[2]))
                expected = {"rho": rho, "log_gap": 2 * math.log(rho),
                            "lower_constraint_species": low[0], "upper_constraint_species": high[0],
                            "common_amplitude_in_box": low[1] <= high[2]}
                self.assert_tree(self.results["rare_certificates"][str(k)][scenario], expected)

    def test_10_printed_pair_is_separate_and_recomputed_with_correct_horizons(self):
        printed = source_rows("cagan_printed_pair.csv")
        expected_values = {"Giraffa_camelopardalis": (99, 24), "Heterocephalus_glaber": (93, 25)}
        self.assertEqual({r["species"] for r in printed}, set(expected_values))
        lookup = {r["species"]: dict(r) for r in self.rows}
        for source in printed:
            species = source["species"]
            rate, years = float(source["mean_sbs_per_genome_year"]), float(source["lifespan80_years"])
            self.assertEqual((rate, years), expected_values[species])
            self.assertIn("Rounded", source["precision_note"])
            self.assertIn("80th percentile", source["precision_note"])
            self.assertIn("Page 522", source["locator"])
            self.assertNotIn("adult_mass_g", source)
            lookup[species].update(m=rate, L=years)
        first, second = lookup["Giraffa_camelopardalis"], lookup["Heterocephalus_glaber"]
        expected = {"burden_giraffe": first["m"] * first["L"],
                    "burden_naked_mole_rat": second["m"] * second["L"],
                    "exact_pair_cmr95": {str(k): independent_pair(first, second, k, "sex_mass") for k in KS}}
        self.assert_tree(self.results["printed_pair_sensitivity"], expected)

    def test_11_maximum_longevity_scaling_and_all_prediction_comparators(self):
        species = sorted(set(self.mutation) & set(self.traits))
        xs = [math.log(float(self.traits[s]["max_longevity_months"]) / 12) for s in species]
        ys = [math.log(float(self.mutation[s]["mean_sbs_per_genome_year"])) for s in species]
        intercept, slope = pairwise_line(xs, ys)
        expected = {"n_species": len(species), "intercept": intercept, "slope": slope, "models": {}}
        fold_slopes = []
        for mode in ("free", "inverse", "intercept"):
            predictions = []
            for j, species_name in enumerate(species):
                tx = [v for i, v in enumerate(xs) if i != j]
                ty = [v for i, v in enumerate(ys) if i != j]
                if mode == "free":
                    a, b = pairwise_line(tx, ty)
                    fold_slopes.append(b)
                else:
                    b = -1 if mode == "inverse" else 0
                    a = math.fsum(y - b * x for x, y in zip(tx, ty)) / len(tx)
                predicted = (a + b * xs[j]) / math.log(10)
                predictions.append({"species": species_name, "predicted_log10_rate": predicted,
                                    "error_log10_rate": ys[j] / math.log(10) - predicted})
            expected["models"][mode] = {"predictions": predictions,
                                        "rmse_log10_rate": math.sqrt(math.fsum(p["error_log10_rate"] ** 2
                                                                               for p in predictions) / len(species))}
        expected["leave_species_out_slope_range"] = [min(fold_slopes), max(fold_slopes)]
        self.assert_tree(self.results["scaling"], expected)

    def test_12_coverage_integrity_and_numerical_boundaries(self):
        matched = {r["species"] for r in self.rows}
        expected = {"mutation_rows": len(self.mutation), "life_history_rows": len(self.traits),
                    "matched_species": len(self.rows), "orders": len({r["order"] for r in self.rows}),
                    "cancer_deaths": sum(r["events"] for r in self.rows),
                    "known_cause_deaths": sum(r["n"] for r in self.rows),
                    "mutation_rows_without_mortality": sorted(set(self.mutation) - matched),
                    "matched_slope_intervals": sum(r["slope_ci"] is not None for r in self.rows),
                    "matched_slope_intervals_reaching_zero": sum(r["slope_ci"] is not None and r["slope_ci"][0] <= 0
                                                                 for r in self.rows)}
        self.assert_tree(self.results["coverage"], expected)
        checks = A.V.run_checks()
        self.assertEqual(checks["status"], "passed")
        self.assertGreater(checks["checks"], 0)
        self.assertIn("no biological validation", checks["scope"])
        for log_h, observed in itertools.product((-40, -10, -1, 0, 1, 4, 10), (0, 0.01, 0.5, 1)):
            expected_risk, expected_nll = decimal_prediction(log_h, observed)
            self.assert_close(A.risk_from_log_h(log_h), expected_risk)
            self.assert_close(A.bernoulli_nll(log_h, observed), expected_nll)
        for observed in (0, 0.3, 1):
            self.assert_close(A.bernoulli_nll(-1000, observed), 1000 * observed)
            self.assertEqual(A.bernoulli_nll(1000, observed), 0 if observed == 1 else math.inf)
        self.assertEqual(A.risk_from_log_h(-1000), 0)
        self.assertEqual(A.risk_from_log_h(1000), 1)
        self.assertEqual(A.bernoulli_nll(-math.inf, 0), 0)
        self.assertEqual(A.bernoulli_nll(-math.inf, 0.3), math.inf)
        for invalid in (-0.1, 1.1, math.nan):
            with self.assertRaises(ValueError):
                A.hazard(invalid)
            with self.assertRaises(ValueError):
                A.bernoulli_nll(0, invalid)
        with self.assertRaises(ValueError):
            A.serializable(float("nan"))
        with self.assertRaises(ValueError):
            A.ols([1, 1, 1], [1, 2, 3])


if __name__ == "__main__":
    unittest.main()
