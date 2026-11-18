"""Offline, independent checks of a nonclinical colony-count reanalysis.

The finite probability-law fixtures below are explicitly synthetic. Source
observations are loaded only from the committed aggregates. Reproducibility,
hash checks and an algebraic alias establish neither source truth nor clinical
validity. No treatment advice. Python standard library only; no network use.
"""

from collections import defaultdict
from contextlib import ExitStack, contextmanager
import csv
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import product
import json
import math
from pathlib import Path
import shutil
import socket
import sys
import tempfile
import unittest
from unittest import mock
import urllib.request
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "persister-colony-reanalysis"
DATA = RESEARCH / "data"
FROZEN_PROTOCOL_SHA256 = (
    "a2a2b44d24977626c5cb395e8832e45378d0db2f9e898c9db23a2cda5b6b35e2"
)


@contextmanager
def offline():
    """Fail immediately if loading or recomputing attempts a network request."""
    error = AssertionError("The persister-colony regression test is offline")
    with ExitStack() as stack:
        for target in (
            "socket.create_connection", "socket.getaddrinfo",
            "socket.socket.connect", "socket.socket.connect_ex",
            "urllib.request.urlopen", "urllib.request.OpenerDirector.open",
        ):
            stack.enter_context(mock.patch(target, side_effect=error))
        yield


def load_analysis(filename="analysis.py", name="persister_colony_reanalysis_under_test"):
    spec = importlib.util.spec_from_file_location(name, RESEARCH / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


with offline():
    A = load_analysis()
    E = load_analysis("extract_supplement.py", "persister_colony_extraction_under_test")


def normalized_hash(path):
    text = Path(path).read_bytes().decode("utf-8")
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def rising(value, count):
    result = F(1)
    for offset in range(count):
        result *= value + offset
    return result


def exact_multinomial_predictive(counts, alpha):
    """Finite rising products, independent of the analysis log-gamma formula."""
    n = sum(counts)
    coefficient = F(math.factorial(n), math.prod(math.factorial(k) for k in counts))
    return coefficient * math.prod(
        rising(a, k) for a, k in zip(alpha, counts)
    ) / rising(sum(alpha), n)


def exact_beta_binomial(n, a, b):
    denominator = rising(a + b, n)
    return [F(math.comb(n, k)) * rising(a, k) * rising(b, n - k)
            / denominator for k in range(n + 1)]


def independent_quantile(probabilities, threshold):
    mass = 0
    for count, probability in enumerate(probabilities):
        mass += probability
        if mass >= threshold:
            return count
    return len(probabilities) - 1


def pgf(kernel, arguments):
    return sum(
        (mass * math.prod(z ** count for z, count in zip(arguments, mark))
         for mark, mass in kernel.items()), F(0)
    )


class FileNormalizationTests(unittest.TestCase):
    def test_text_digests_ignore_only_newline_encoding_and_paths_use_posix(self):
        content = "source,measurement\nDiFi,1.25\nWiDr,2.5\n"
        expected = hashlib.sha256(content.encode("utf-8")).hexdigest()
        with tempfile.TemporaryDirectory(prefix="persister-newline-test-") as temporary:
            for index, newline in enumerate(("\n", "\r\n", "\r")):
                path = Path(temporary) / f"example_{index}.csv"
                path.write_bytes(content.replace("\n", newline).encode("utf-8"))
                self.assertEqual(A.normalized_text_hash(path), expected)
                self.assertEqual(normalized_hash(path), expected)
        self.assertEqual((Path("data") / "colony_counts.csv").as_posix(),
                         "data/colony_counts.csv")

    def test_source_exponents_survive_rich_text_extraction(self):
        node = ET.fromstring(
            '<si xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            '<r><t>10</t></r><r><rPr><vertAlign val="superscript"/></rPr>'
            '<t>5</t></r></si>'
        )
        plain, rich = E.string_value(node)
        self.assertEqual((plain, rich), ("105", "10^5"))
        self.assertEqual(E.number_from_text(rich), 100000)
        self.assertEqual(E.number_from_text(plain), 105)
        for text, expected in (("1.2 ∗ 10^8", 120000000),
                               ("2.5 × 10^−6", 0.0000025),
                               ("7.2999999999999995E-2", 0.073)):
            self.assertAlmostEqual(E.number_from_text(text), expected, places=14)
        for malformed in ("NaN", "Infinity", "_1010", "10^−^1", "10^x"):
            with self.assertRaises(ValueError):
                E.number_from_text(malformed)


class ProbabilityLawTests(unittest.TestCase):
    def test_full_vector_predictive_matches_exact_finite_rational_law(self):
        # Exhaust the support, including zero counts and the multinomial factor.
        alpha = (F(1, 2), F(3, 2), F(5, 2))
        for n in (0, 1, 4, 8):
            total = F(0)
            for early in range(n + 1):
                for late in range(n - early + 1):
                    counts = (early, late, n - early - late)
                    exact = exact_multinomial_predictive(counts, alpha)
                    observed = math.exp(A.dirichlet_multinomial_logpmf(counts, alpha))
                    self.assertAlmostEqual(observed, float(exact), places=13)
                    total += exact
            self.assertEqual(total, 1)

    def test_beta_binomial_probabilities_and_quantiles_match_exact_products(self):
        for n, a, b in ((0, F(1, 2), F(1, 2)),
                        (12, F(3, 2), F(5, 2)),
                        (19, F(1, 4), F(9, 4))):
            exact = exact_beta_binomial(n, a, b)
            self.assertEqual(sum(exact), 1)
            observed = A.beta_binomial_pmf(n, a, b)
            for x, y in zip(observed, exact):
                self.assertAlmostEqual(x, float(y), places=13)
            expected_interval = [independent_quantile(exact, F(1, 40)),
                                 independent_quantile(exact, F(39, 40))]
            self.assertEqual(A.beta_binomial_interval(n, a, b), expected_interval)
            lo, hi = expected_interval
            self.assertLessEqual(sum(exact[:lo]), F(1, 40))
            self.assertLessEqual(sum(exact[hi + 1:]), F(1, 40))
            self.assertGreaterEqual(sum(exact[lo:hi + 1]), F(19, 20))

    def test_beta_binomial_on_assay_sized_support_by_independent_recurrence(self):
        # This initial mass is safely representable; recurrence uses no gamma.
        n, a, b = 1920, 0.5, 1.5
        initial = math.prod((b + j) / (a + b + j) for j in range(n))
        recurrence = [initial]
        for k in range(n):
            recurrence.append(recurrence[-1] * (n - k) * (a + k)
                              / ((k + 1) * (b + n - k - 1)))
        self.assertAlmostEqual(math.fsum(recurrence), 1.0, places=12)
        actual = A.beta_binomial_pmf(n, a, b)
        self.assertAlmostEqual(math.fsum(actual), 1.0, places=13)
        self.assertLess(max(abs(x - y) for x, y in zip(actual, recurrence)), 1e-11)
        self.assertAlmostEqual(math.fsum(k * p for k, p in enumerate(actual)),
                               n * a / (a + b), places=7)
        self.assertEqual(A.beta_binomial_interval(n, a, b),
                         [independent_quantile(recurrence, 0.025),
                          independent_quantile(recurrence, 0.975)])

    def test_concentration_caps_preserve_training_mean_without_adding_wells(self):
        rows = [{"total_wells": 20, "early": 2, "late": 3},
                {"total_wells": 20, "early": 1, "late": 6}]
        expected = (F(7, 2), F(19, 2), F(57, 2))
        self.assertEqual(A.posterior_alpha(rows), tuple(map(float, expected)))
        self.assertEqual(A.posterior_alpha(rows, cap=100), tuple(map(float, expected)))
        capped = A.posterior_alpha(rows, cap=10)
        for actual, value in zip(capped, expected):
            self.assertAlmostEqual(actual, float(10 * value / sum(expected)), places=14)
        self.assertAlmostEqual(sum(capped), 10, places=13)


class ObservationalAliasTests(unittest.TestCase):
    def test_vector_alias_with_unequal_kernels_and_invisible_events(self):
        components = [
            (F(2, 3), {(0, 0): F(1, 4), (0, 1): F(1, 4), (1, 2): F(1, 2)}),
            (F(5, 7), {(0, 0): F(2, 5), (0, 2): F(1, 5), (2, 1): F(2, 5)}),
        ]
        effective, delayed_kernel = A.positive_mark_alias(components)
        expected_effective = sum(rate * (1 - kernel[(0, 0)])
                                 for rate, kernel in components)
        self.assertEqual(effective, expected_effective)
        self.assertEqual(sum(delayed_kernel.values()), 1)
        self.assertNotIn((0, 0), delayed_kernel)
        self.assertTrue(all(value >= 0 for value in delayed_kernel.values()))
        # A common unrelated pre-existing component stays fixed in both models.
        background = {(1, 0): F(1, 3), (2, 1): F(2, 3)}
        for arguments in product((F(0), F(1, 3), F(2, 3), F(1)), repeat=2):
            acquired_log_pgf = sum(rate * (pgf(kernel, arguments) - 1)
                                   for rate, kernel in components)
            delayed_log_pgf = effective * (pgf(delayed_kernel, arguments) - 1)
            self.assertEqual(acquired_log_pgf, delayed_log_pgf)
            fixed_log_pgf = F(2, 5) * (pgf(background, arguments) - 1)
            for theta in (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)):
                allocated = (fixed_log_pgf + (1 - theta) * delayed_log_pgf
                             + theta * acquired_log_pgf)
                self.assertEqual(allocated, fixed_log_pgf + acquired_log_pgf)

    def test_all_invisible_events_have_no_identifiable_positive_kernel(self):
        self.assertEqual(A.positive_mark_alias([(F(3, 4), {(0, 0): F(1)})]),
                         (F(0), {}))
        for invalid in (
            [(F(-1), {(0, 1): F(1)})],
            [(F(1), {(0, 1): F(3, 4)})],
            [(F(1), {(1,): F(1, 2), (1, 2): F(1, 2)})],
        ):
            with self.assertRaises(ValueError):
                A.positive_mark_alias(invalid)


class PopulationCalibrationTests(unittest.TestCase):
    def test_net_slope_does_not_determine_gross_birth_and_death(self):
        # A synthetic declining proxy curve with unequal measurement times.
        points = [(day, 123 * math.exp(-0.06 * day)) for day in (14, 19, 25, 32)]
        slope, intercept = A.log_slope(points)
        self.assertAlmostEqual(slope, -0.06, places=14)
        self.assertAlmostEqual(intercept, math.log(123), places=13)
        rows = [{"cell_line": "synthetic", "clone": "test", "assay": "proxy",
                 "measurement": "arbitrary_units", "biological_replicate": "1",
                 "day": day, "value": value} for day, value in points]
        result = A.calibrate_dynamics(rows)
        self.assertEqual(len(result["fits"]), 2)
        for fit in result["fits"]:
            rates = fit["turnover_sensitivity"]
            self.assertEqual(len(rates), 3)
            for pair in rates:
                self.assertGreaterEqual(pair["birth_per_day"], 0)
                self.assertGreaterEqual(pair["death_per_day"], 0)
                self.assertAlmostEqual(pair["birth_per_day"] - pair["death_per_day"],
                                       -0.06, places=14)
        self.assertEqual(A.calibrate_dynamics([])["fits"], [])
        self.assertTrue(A.calibrate_dynamics([])["gaps"])


class CommittedReanalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with offline():
            cls.result = A.analyze(DATA)
            cls.colonies = A.load_colonies(DATA / "colony_counts.csv")
        cls.groups = defaultdict(list)
        for row in cls.colonies:
            cls.groups[(row["cell_line"], row["clone"])].append(row)

    def test_original_local_freeze_objects_bind_protocol_and_prefit_deviation(self):
        evidence = json.loads((DATA / "protocol_freeze.json").read_text(encoding="utf-8"))

        def object_hash(kind, payload):
            header = f"{kind} {len(payload)}\0".encode("ascii")
            return hashlib.sha1(header + payload).hexdigest()

        def tree_entries(payload):
            entries, position = {}, 0
            while position < len(payload):
                space = payload.index(b" ", position)
                nul = payload.index(b"\0", space)
                mode = payload[position:space].decode("ascii")
                name = payload[space + 1:nul].decode("utf-8")
                digest = payload[nul + 1:nul + 21].hex()
                self.assertEqual(len(digest), 40)
                entries[name] = (mode, digest)
                position = nul + 21
            self.assertEqual(position, len(payload))
            return entries

        protocol = (RESEARCH / "PROTOCOL.md").read_text(encoding="utf-8")
        protocol = protocol.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")
        for label in ("protocol", "source_deviation"):
            record = evidence[label]
            commit_payload = record["payload"].encode("utf-8")
            self.assertEqual(object_hash("commit", commit_payload), record["sha1"])
            trees = record["trees"]
            self.assertTrue(record["payload"].startswith("tree " + trees["."]["sha1"] + "\n"))
            decoded = {}
            for name, tree in trees.items():
                raw = bytes.fromhex(tree["payload_hex"])
                self.assertEqual(object_hash("tree", raw), tree["sha1"])
                decoded[name] = tree_entries(raw)
            self.assertEqual(decoded["."]["research"][1], trees["research"]["sha1"])
            self.assertEqual(decoded["research"]["persister-colony-reanalysis"][1],
                             trees["research/persister-colony-reanalysis"]["sha1"])
            leaf = decoded["research/persister-colony-reanalysis"]
            self.assertEqual(leaf["PROTOCOL.md"], ("100644", object_hash("blob", protocol)))
            if label == "protocol":
                self.assertIn("parent " + evidence["base_commit"] + "\n", record["payload"])
                self.assertEqual(set(leaf), {"PROTOCOL.md"})
            else:
                self.assertIn("parent " + evidence["protocol"]["sha1"] + "\n", record["payload"])
                self.assertEqual(set(leaf), {"PROTOCOL.md", "README.md"})
                self.assertEqual(leaf["README.md"][1], object_hash(
                    "blob", record["prefit_readme"].encode("utf-8")))
                self.assertIn("Deviation D1, recorded before fitting", record["prefit_readme"])
        self.assertTrue(evidence["protocol"]["sha1"].startswith("cb76f03"))
        self.assertTrue(evidence["source_deviation"]["sha1"].startswith("431881a"))

    def test_frozen_protocol_and_normalized_manifest_hashes(self):
        manifest_path = DATA / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        expected = manifest["normalized_sha256"]
        self.assertEqual(manifest["protocol_sha256"], FROZEN_PROTOCOL_SHA256)
        self.assertEqual(normalized_hash(RESEARCH / "PROTOCOL.md"),
                         FROZEN_PROTOCOL_SHA256)
        actual_paths = {
            path.relative_to(ROOT).as_posix()
            for path in RESEARCH.rglob("*")
            if path.is_file() and "__pycache__" not in path.parts
            and path not in (manifest_path, RESEARCH / "README.md")
        }
        self_key = Path(__file__).resolve().relative_to(ROOT).as_posix()
        if self_key in expected:
            actual_paths.add(self_key)
        self.assertEqual(set(expected), actual_paths)
        for relative, digest in expected.items():
            with self.subTest(path=relative):
                path = Path(relative)
                self.assertEqual(path.as_posix(), relative)
                self.assertNotIn("\\", relative)
                self.assertFalse(path.is_absolute())
                self.assertNotIn("..", path.parts)
                self.assertEqual(normalized_hash(ROOT / path), digest)
        for relative, digest in self.result["input_hashes"].items():
            self.assertEqual(Path(relative).as_posix(), relative)
            self.assertEqual(normalized_hash(RESEARCH / relative), digest)

    def test_readme_numerical_report_recomputes_offline_exactly(self):
        readme = (RESEARCH / "README.md").read_text(encoding="utf-8")
        begin, end = "<!-- BEGIN COMPUTED REPORT -->", "<!-- END COMPUTED REPORT -->"
        self.assertEqual(readme.count(begin), 1)
        self.assertEqual(readme.count(end), 1)
        observed = begin + readme.split(begin, 1)[1].split(end, 1)[0] + end
        self.assertEqual(observed, A.render_markdown(self.result).strip())
        self.assertIn("nonclinical", readme.lower())
        self.assertIn("no treatment advice", readme.lower())

    def test_every_replicate_is_held_out_whole_and_training_sums_are_independent(self):
        expected_settings = {
            "primary": (0.5, None), "cap_10": (0.5, 10),
            "cap_100": (0.5, 100), "cap_1000": (0.5, 1000),
            "prior_0.25": (0.25, None), "prior_1": (1.0, None),
        }
        scenarios = {row["label"]: row for row in self.result["predictive_scenarios"]}
        self.assertEqual(set(scenarios), set(expected_settings))
        for label, scenario in scenarios.items():
            prior, cap = expected_settings[label]
            expected_keys = {(line, clone, row["biological_replicate"])
                             for (line, clone), rows in self.groups.items()
                             if len(rows) >= 2 for row in rows}
            observed_keys = {(row["cell_line"], row["clone"], row["held_out_replicate"])
                             for row in scenario["folds"]}
            self.assertEqual(expected_keys, observed_keys)
            self.assertEqual(len(expected_keys), scenario["n_scored_folds"])
            for fold in scenario["folds"]:
                key = fold["cell_line"], fold["clone"]
                held_id = fold["held_out_replicate"]
                training = [row for row in self.groups[key]
                            if row["biological_replicate"] != held_id]
                held = next(row for row in self.groups[key]
                            if row["biological_replicate"] == held_id)
                self.assertEqual(set(fold["training_replicates"]),
                                 {row["biological_replicate"] for row in training})
                self.assertNotIn(held_id, fold["training_replicates"])
                alpha = [prior + sum(row["early"] for row in training),
                         prior + sum(row["late"] for row in training),
                         prior + sum(row["total_wells"] - row["early"] - row["late"]
                                     for row in training)]
                if cap is not None and sum(alpha) > cap:
                    alpha = [value * cap / sum(alpha) for value in alpha]
                for actual, expected in zip(fold["alpha"], alpha):
                    self.assertAlmostEqual(actual, expected, places=11)
                n, early, late = held["total_wells"], held["early"], held["late"]
                self.assertEqual(fold["counts"], [early, late, n - early - late])
                self.assertEqual(fold["total_wells"], n)
                for actual, value in zip(fold["predictive_means"], alpha):
                    self.assertAlmostEqual(actual, n * value / sum(alpha), places=10)
                self.assertEqual(fold["log_score_m0"], fold["log_score_m1"])
                self.assertEqual(fold["early_covered"],
                                 fold["early_interval"][0] <= early <= fold["early_interval"][1])
                self.assertEqual(fold["late_covered"],
                                 fold["late_interval"][0] <= late <= fold["late_interval"][1])
            self.assertAlmostEqual(scenario["sum_log_score_m0"],
                                   math.fsum(row["log_score_m0"] for row in scenario["folds"]),
                                   places=12)
            self.assertEqual(scenario["sum_log_score_m0"], scenario["sum_log_score_m1"])
            self.assertEqual(scenario["score_difference"], 0)
            self.assertEqual(scenario["early_covered"],
                             sum(row["early_covered"] for row in scenario["folds"]))
            self.assertEqual(scenario["late_covered"],
                             sum(row["late_covered"] for row in scenario["folds"]))

    def test_mutating_a_held_out_outcome_does_not_change_its_predictive_distribution(self):
        original = self.colonies[0]
        identity = (original["cell_line"], original["clone"],
                    original["biological_replicate"])
        self.assertGreater(original["total_wells"] - original["early"] - original["late"], 3)
        with tempfile.TemporaryDirectory(prefix="persister-colony-test-") as temporary:
            copied = Path(temporary) / "data"
            shutil.copytree(DATA, copied)
            csv_path = copied / "colony_counts.csv"
            with csv_path.open(encoding="utf-8", newline="") as stream:
                reader = csv.DictReader(stream)
                fields, rows = reader.fieldnames, list(reader)
            for row in rows:
                if (row["cell_line"], row["clone"], row["biological_replicate"]) == identity:
                    row["early"] = str(int(row["early"]) + 1)
                    row["late"] = str(int(row["late"]) + 2)
            with csv_path.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            with offline():
                changed = A.analyze(copied)
        for before, after in zip(self.result["predictive_scenarios"],
                                 changed["predictive_scenarios"]):
            self.assertEqual(before["label"], after["label"])
            def selected(scenario):
                return next(fold for fold in scenario["folds"]
                            if (fold["cell_line"], fold["clone"],
                                fold["held_out_replicate"]) == identity)
            left, right = selected(before), selected(after)
            for name in ("alpha", "predictive_means", "early_interval", "late_interval",
                         "training_replicates"):
                self.assertEqual(left[name], right[name])
            self.assertNotEqual(left["counts"], right["counts"])
        self.assertEqual(self.result["dynamics"], changed["dynamics"])

    def test_published_cells_preserve_numeric_values_units_and_exact_locators(self):
        with (DATA / "published_parameters.csv").open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 32)
        self.assertEqual({int(row["source_table"]) for row in rows}, {1, 2, 4, 5})
        source_cells = {(row["sheet"], row["cell"]) for row in rows}
        self.assertEqual(len(source_cells), len(rows))
        for row in rows:
            with self.subTest(locator=row["locator"]):
                self.assertEqual(row["locator"], row["sheet"] + "!" + row["cell"])
                self.assertEqual(row["source_file"], "41588_2022_1105_MOESM4_ESM.xlsx")
                self.assertEqual(row["licence"], "CC-BY-4.0")
                # These selected cells are numeric, not ambiguous rich strings.
                # This checks internal transcription consistency, not source truth;
                # the complete workbook is deliberately not bundled in offline CI.
                self.assertEqual(float(row["raw_cell_value"]), float(row["value"]))
                self.assertEqual(E.number_from_text(row["rich_text_value"]), float(row["value"]))
                if int(row["source_table"]) == 5:
                    expected_scale = "1e-7" if row["parameter"].startswith("sensitive_") else "1e-6"
                    self.assertEqual(row["units"], expected_scale + "_per_cell_per_day")
                if row["statistic"].startswith("hdi_"):
                    self.assertIn("3%/97%", row["note"])
                    self.assertIn("2.5%/97.5%", row["note"])

    def test_published_rate_scales_ratios_and_all_turnover_pairs_recompute(self):
        with (DATA / "published_parameters.csv").open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        lookup = {(int(row["source_table"]), row["cell_line"], row["clone"],
                   row["parameter"], row["statistic"]): row for row in rows}
        context = self.result["published_context"]
        self.assertFalse(context["used_to_fit_colony_counts"])
        self.assertFalse(context["raw_population_recalibration"])
        self.assertEqual(len(context["rate_comparisons"]), 2)
        for comparison in context["rate_comparisons"]:
            values = {}
            for population, statistic in product(("sensitive", "persister"),
                                                  ("estimate", "uncertainty_sd")):
                source = lookup[(5, comparison["cell_line"], comparison["clone"],
                                 f"{population}_chronological_resistance_rate", statistic)]
                value = F(source["value"]) * F(source["units"].split("_", 1)[0])
                key = f"{population}_{statistic}_per_cell_per_day"
                values[key] = value
                self.assertAlmostEqual(comparison[key], float(value), places=18)
            ratio = values["persister_estimate_per_cell_per_day"] / values["sensitive_estimate_per_cell_per_day"]
            self.assertAlmostEqual(comparison["ratio_of_published_point_estimates"],
                                   float(ratio), places=12)
        anchors = context["fallback_net_loss_anchors"]
        self.assertEqual(len(anchors), 6)
        self.assertEqual(sum(len(row["turnover_sensitivity"]) for row in anchors), 18)
        for anchor in anchors:
            source = lookup[(4, anchor["cell_line"], anchor["clone"],
                             "D_p", anchor["source_statistic"])]
            loss = F(source["value"])
            self.assertEqual(anchor["source_locator"], source["locator"])
            self.assertEqual(anchor["published_D_p_per_day"], float(loss))
            self.assertEqual(anchor["assumed_net_growth_per_day"], -float(loss))
            self.assertEqual([row["birth_per_day"] for row in anchor["turnover_sensitivity"]],
                             [0, 0.5, 1])
            for rates in anchor["turnover_sensitivity"]:
                birth = F(str(rates["birth_per_day"]))
                self.assertAlmostEqual(rates["death_per_day"], float(birth + loss), places=14)
                self.assertAlmostEqual(rates["birth_per_day"] - rates["death_per_day"],
                                       -float(loss), places=14)
                self.assertEqual(rates["net_per_day"], -float(loss))

    def test_changing_published_context_does_not_change_colony_fits_or_scores(self):
        with tempfile.TemporaryDirectory(prefix="persister-context-test-") as temporary:
            copied = Path(temporary) / "data"
            shutil.copytree(DATA, copied)
            csv_path = copied / "published_parameters.csv"
            with csv_path.open(encoding="utf-8", newline="") as stream:
                reader = csv.DictReader(stream)
                fields, rows = reader.fieldnames, list(reader)
            for row in rows:
                # Deliberate synthetic perturbation in a disposable copy only.
                row["value"] = str(2 * float(row["value"]))
            with csv_path.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            with offline():
                changed = A.analyze(copied)
        self.assertNotEqual(self.result["published_context"], changed["published_context"])
        for key in ("clone_summaries", "predictive_scenarios", "dynamics"):
            self.assertEqual(self.result[key], changed[key])

    def test_source_classification_binds_report_to_incomplete_raw_audit(self):
        metadata = json.loads((DATA / "sources.json").read_text(encoding="utf-8"))
        status = self.result["source_status"]
        self.assertIs(metadata["raw_archive_verified"], False)
        self.assertIs(metadata["raw_population_calibration_completed"], False)
        self.assertEqual(status["analysis_status"], metadata["analysis_status"])
        self.assertTrue(status["analysis_status"].startswith("partial_"))
        self.assertIs(status["raw_archive_verified"], False)
        self.assertIs(status["raw_population_calibration_completed"], False)
        source_types = {row["source_type"] for row in self.colonies}
        self.assertEqual(source_types, {"published_figure_annotation_not_raw_deposit"})
        self.assertEqual(status["colony_source_types"], sorted(source_types))
        self.assertEqual(self.result["population_observation_rows"], 0)
        self.assertEqual(self.result["dynamics"]["fits"], [])
        self.assertTrue(self.result["dynamics"]["gaps"])
        self.assertIn("requested raw audit incomplete", self.result["verdict"])
        report = A.render_markdown(self.result)
        self.assertIn(self.result["verdict"], report)
        self.assertIn("Raw archive checksum verified: **false**", report)
        self.assertIn("Raw population calibration completed: **false**", report)
        self.assertIn("published_figure_annotation_not_raw_deposit", report)

    def test_empirical_totals_and_all_origin_allocations_preserve_observables(self):
        for summary in self.result["clone_summaries"]:
            rows = self.groups[summary["cell_line"], summary["clone"]]
            n = sum(row["total_wells"] for row in rows)
            early, late = (sum(row[field] for row in rows) for field in ("early", "late"))
            self.assertEqual(summary["total_wells"], n)
            self.assertEqual(summary["early"], early)
            self.assertEqual(summary["late"], late)
            self.assertEqual(summary["negative"], n - early - late)
            self.assertEqual(summary["late_risk_wells"], n - early)
            self.assertAlmostEqual(summary["empirical_late_fraction_of_risk"], late / (n - early))
            expected_p = [(early + 0.5) / (n + 1.5), (late + 0.5) / (n + 1.5),
                          (n - early - late + 0.5) / (n + 1.5)]
            self.assertEqual([row["theta"] for row in summary["aliases"]],
                             [0, 0.25, 0.5, 0.75, 1])
            a, b = summary["effective_early_intensity"], summary["effective_late_intensity"]
            for alias in summary["aliases"]:
                theta, delayed, acquired = alias["theta"], alias["b_delayed"], alias["b_acquired"]
                self.assertGreaterEqual(delayed, 0)
                self.assertGreaterEqual(acquired, 0)
                self.assertAlmostEqual(delayed, (1 - theta) * b, places=14)
                self.assertAlmostEqual(acquired, theta * b, places=14)
                self.assertAlmostEqual(delayed + acquired, b, places=14)
                probabilities = A.alias_probabilities(a, delayed, acquired)
                self.assertAlmostEqual(sum(probabilities), 1, places=14)
                for actual, expected in zip(probabilities, expected_p):
                    self.assertAlmostEqual(actual, expected, places=13)
                for actual, expected in zip(alias["probabilities"], expected_p):
                    self.assertAlmostEqual(actual, expected, places=13)
            conversions = summary["latent_intensity_scalings"]
            self.assertEqual({(row["detection_efficiency"], row["kernel_scale"])
                              for row in conversions},
                             set(product((0.1, 0.5, 1.0), (0.25, 1.0, 4.0))))
            for conversion in conversions:
                # Binary occupancy depends on Pr(mark != 0), not positive size.
                q = conversion["detection_efficiency"]
                latent = conversion["required_latent_intensity_per_well"]
                self.assertTrue(0 < q <= 1)
                self.assertAlmostEqual(latent * q, b, places=13)


if __name__ == "__main__":
    unittest.main()
