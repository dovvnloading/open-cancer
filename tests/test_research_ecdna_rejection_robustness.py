"""Nonclinical rejection-count sensitivity: independent offline stdlib checks.

These calculations test mathematical/software claims, not the truth of the
published observations or an allocation mechanism. No treatment advice.
"""

import ast
import csv
from fractions import Fraction
import hashlib
import importlib.util
from itertools import combinations, product
import math
from pathlib import Path
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "research" / "ecdna-rejection-robustness"
SOURCE = ROOT / "research" / "ecdna-segregation-clustering-data" / "data" / "lundberg.csv"
F = Fraction
LEVEL = F(1, 20)
SOURCE_SHA256 = "d6d89d7beb6757db0a93174674c3e0baa38e363f98c5835f66007ad4b917ebea"
PROTOCOL_SHA256 = "bed446347c21936e5561d712b03cc361b026fc8fb9e96b02d2df27684b414185"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    instance = importlib.util.module_from_spec(spec)
    sys.modules[name] = instance
    spec.loader.exec_module(instance)
    return instance


def digest(path):
    content = path.read_bytes().decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def bernoulli_convolution(n, p):
    """Independent probability-generating-function multiplication."""
    masses = [F(1)]
    for _ in range(n):
        following = [F(0)] * (len(masses) + 1)
        for s, mass in enumerate(masses):
            following[s] += mass * (1 - p)
            following[s + 1] += mass * p
        masses = following
    return masses


def polya_urn(n, p, rho):
    """Independent sequential predictive construction of the count law."""
    if rho == 0:
        return bernoulli_convolution(n, p)
    if rho == 1:
        return [1 - p] + [F(0)] * (n - 1) + [p]
    if p in (0, 1):
        return [F(s == n * p) for s in range(n + 1)]
    concentration = (1 - rho) / rho
    alpha = p * concentration
    masses = [F(1)]
    for t in range(n):
        following = [F(0)] * (t + 2)
        for s, mass in enumerate(masses):
            chance = (alpha + s) / (concentration + t)
            following[s] += mass * (1 - chance)
            following[s + 1] += mass * chance
        masses = following
    return masses


def feasible_vertices(n, p, rho, support_limit=None):
    """Enumerate the primal LP independently of the production dual formula.

    With normalization and fixed mean, a vertex with slack second moment has
    at most two atoms. An active second moment permits at most three atoms.
    A three-atom solution is solved by Lagrange interpolation of moments.
    """
    mean = n * p
    second = n * (n - 1) * (p * p + rho * p * (1 - p))
    out = []
    maximum = n if support_limit is None else support_limit
    if not 0 <= maximum <= n:
        raise ValueError("Support limit must lie between zero and n")

    def retain(law):
        if any(value < 0 for value in law.values()):
            return
        if sum(law.values()) != 1 or sum(s * q for s, q in law.items()) != mean:
            return
        if sum(s * (s - 1) * q for s, q in law.items()) > second:
            return
        out.append(law)

    for s in range(maximum + 1):
        retain({s: F(1)})
    for a, b in combinations(range(maximum + 1), 2):
        retain({a: (b - mean) / (b - a), b: (mean - a) / (b - a)})
    for support in combinations(range(maximum + 1), 3):
        law = {}
        for s in support:
            a, b = [t for t in support if t != s]
            law[s] = (second + (1 - a - b) * mean + a * b) / ((s - a) * (s - b))
        retain(law)
    if not out:
        raise AssertionError("The declared nonnegative-correlation cap must be feasible")
    return out


def holm_closed(values):
    """Independent closed Bonferroni expression; no sorted-step algorithm."""
    keys = tuple(values)
    answer = {}
    for key in keys:
        candidates = []
        for size in range(1, len(keys) + 1):
            for subset in combinations(keys, size):
                if key in subset:
                    candidates.append(size * min(values[item] for item in subset))
        answer[key] = min(F(1), max(candidates))
    return answer


class TestEcDNARejectionRobustness(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.start = time.perf_counter()
        cls.analysis = load_module("ecdna_rejection_robustness_analysis_test", HERE / "analysis.py")
        cls.verifier = load_module("ecdna_rejection_robustness_verifier_test", HERE / "verifier.py")
        cls.result = cls.analysis.results()
        cls.result_seconds = time.perf_counter() - cls.start

    def test_every_reported_number_is_recomputed(self):
        text = (HERE / "README.md").read_text(encoding="utf-8")
        marker = "<!-- BEGIN RECOMPUTED RESULTS -->\n"
        end = "<!-- END RECOMPUTED RESULTS -->"
        self.assertEqual(text.count(marker), 1)
        self.assertEqual(text.count(end), 1)
        actual = text.split(marker, 1)[1].split(end, 1)[0]
        self.assertEqual(actual, self.analysis.render_report(self.result) + "\n")
        self.assertLess(self.result_seconds, 60)

    def test_frozen_protocol_and_source_hashes(self):
        self.assertEqual(digest(SOURCE), SOURCE_SHA256)
        self.assertEqual(digest(HERE / "PROTOCOL.md"), PROTOCOL_SHA256)
        self.assertEqual(self.analysis.PROTOCOL_COMMIT, "d58107d841841dfb42577e10f1f6065c11c909bc")
        paths = self.analysis.hashes()
        required = (SOURCE, HERE / "PROTOCOL.md", HERE / "analysis.py", HERE / "verifier.py")
        for path in required:
            key = path.relative_to(ROOT).as_posix()
            self.assertIn(key, paths)
            self.assertEqual(paths[key], digest(path))
        for key, value in paths.items():
            self.assertNotIn("\\", key)
            self.assertEqual(key, Path(key).as_posix())
            self.assertFalse(Path(key).is_absolute())
            self.assertEqual(len(value), 64)

    def test_hashes_normalize_crlf_cr_and_lf(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "line-endings.txt"
            expected = hashlib.sha256("alpha\nbeta\ngamma\n".encode("utf-8")).hexdigest()
            for content in (b"alpha\nbeta\ngamma\n", b"alpha\r\nbeta\r\ngamma\r\n",
                            b"alpha\rbeta\rgamma\r", b"alpha\r\nbeta\rgamma\n"):
                path.write_bytes(content)
                self.assertEqual(self.analysis.normalized_hash(path), expected)
                self.assertEqual(self.analysis.normalized_text(path), "alpha\nbeta\ngamma\n")

    def test_source_selection_and_rounding_are_independent(self):
        with SOURCE.open(encoding="utf-8", newline="") as handle:
            source_rows = list(csv.DictReader(handle))
        source = {row["id"]: row for row in source_rows}
        self.assertEqual(len(source_rows), len(source))
        chosen = self.analysis.load_data()
        self.assertEqual(tuple(row["line"] for row in chosen), ("CHP-212", "LA-N-5", "MC-IX"))
        identifiers = {"CHP-212": "chp212", "LA-N-5": "lan5", "MC-IX": "mcix"}
        expected = {"CHP-212": (100, 2), "LA-N-5": (100, 10), "MC-IX": (50, 1)}
        selected_ids = {"lundberg_alpha"}
        for row in chosen:
            prefix = "lundberg_" + identifiers[row["line"]]
            nr, pr = source[prefix + "_n"], source[prefix + "_pct"]
            selected_ids.update((nr["id"], pr["id"]))
            n, percent = int(nr["value"]), F(pr["value"])
            compatible = tuple(k for k in range(n + 1)
                               if percent - F(1, 2) <= F(100 * k, n) <= percent + F(1, 2))
            self.assertEqual(tuple(row["counts"]), compatible)
            self.assertEqual((row["n"], row["k"]), expected[row["line"]])
            self.assertEqual(row["k"], min(compatible))
            self.assertEqual(row["nominal"], F(1, 100))
            self.assertEqual(nr["quantity"], "scored_divisions_n")
            self.assertEqual(pr["quantity"], "rejected_divisions_percent")
        self.assertEqual(len(selected_ids), 7)
        for identifier in selected_ids:
            row = source[identifier]
            self.assertEqual(row["condition"], "untreated")
            self.assertEqual(row["doi"], "10.1371/journal.pone.0003099")
            self.assertEqual(row["access_date"], "2026-11-17")
            for key in ("source_url", "locator", "reuse_basis", "notes"):
                self.assertTrue(row[key].strip())
            self.assertIn("Creative Commons Attribution", row["reuse_basis"])
            self.assertIn("Statistical distribution", row["locator"])

    def test_independent_binomial_tails_and_published_endpoints(self):
        for n, p in product((1, 3, 7, 12), (F(0), F(1, 100), F(1, 3), F(1))):
            law = bernoulli_convolution(n, p)
            self.assertEqual(sum(law), 1)
            for k in range(1, n + 1):
                self.assertEqual(self.analysis.binomial_tail(n, k, p), sum(law[k:]))
        for row in self.analysis.load_data():
            n, k, p = row["n"], row["k"], row["nominal"]
            direct = sum(F(math.comb(n, j)) * p ** j * (1 - p) ** (n - j)
                         for j in range(k, n + 1))
            self.assertEqual(self.analysis.binomial_tail(n, k, p), direct)
        nominal = self.analysis.scenario("independent", F(1, 100), F(0))
        self.assertEqual(tuple(nominal["rejected"]), ("LA-N-5",))
        self.assertEqual(nominal["adjusted"], holm_closed(nominal["raw"]))

    def test_beta_tails_against_polya_urn_and_moments(self):
        for n, p, rho in product((1, 4, 8), (F(0), F(1, 100), F(2, 5), F(1)),
                                 (F(0), F(1, 1000), F(1, 7), F(3, 4), F(1))):
            law = polya_urn(n, p, rho)
            self.assertEqual(sum(law), 1)
            self.assertEqual(sum(s * q for s, q in enumerate(law)), n * p)
            self.assertEqual(sum(s * (s - 1) * q for s, q in enumerate(law)),
                             n * (n - 1) * (p * p + rho * p * (1 - p)))
            for k in range(1, n + 1):
                self.assertEqual(self.analysis.beta_tail(n, k, p, rho), sum(law[k:]))

    def test_sharp_bound_against_full_small_primal_vertex_enumeration(self):
        for n, p, rho in product((2, 4, 6), (F(0), F(1, 100), F(1, 4), F(1, 2), F(1)),
                                 (F(0), F(1, 10), F(1))):
            vertices = feasible_vertices(n, p, rho)
            for k in range(1, n + 1):
                optimum = max(sum(q for s, q in law.items() if s >= k) for law in vertices)
                self.assertEqual(self.analysis.sharp_tail(n, k, p, rho), optimum,
                                 (n, k, p, rho))

    def test_primal_dual_certificates_and_integer_ties(self):
        # Includes mean at/beyond threshold, degenerate marginal endpoints,
        # and exact active-branch ties. A displayed bound alone is insufficient.
        cases = [(n, k, p, rho) for n in (2, 5, 8) for k in range(1, n + 1)
                 for p in (F(0), F(1, 100), F(1, 4), F(1, 2), F(1))
                 for rho in (F(0), F(1, 5), F(1))]
        cases += [(row["n"], row["k"], F(1, 100), rho)
                  for row in self.analysis.load_data()
                  for rho in (F(0), F(1, 1000), F(1, 100), F(1))]
        for n, k, p, rho in cases:
            certificate = self.analysis.sharp_certificate(n, k, p, rho)
            law, (c0, c1, c2) = certificate["law"], certificate["dual"]
            self.assertTrue(all(isinstance(s, int) and 0 <= s <= n for s in law))
            self.assertTrue(all(q >= 0 for q in law.values()))
            self.assertEqual(sum(law.values()), 1)
            self.assertEqual(sum(s * q for s, q in law.items()), n * p)
            second = n * (n - 1) * (p * p + rho * p * (1 - p))
            self.assertLessEqual(sum(s * (s - 1) * q for s, q in law.items()), second)
            self.assertEqual(certificate["mean"], n * p)
            self.assertEqual(certificate["second_cap"], second)
            self.assertGreaterEqual(c2, 0)
            for s in range(n + 1):
                self.assertGreaterEqual(c0 + c1 * s + c2 * s * (s - 1), F(s >= k))
            primal = sum(q for s, q in law.items() if s >= k)
            dual = c0 + c1 * n * p + c2 * second
            self.assertEqual(primal, dual)
            self.assertEqual(primal, certificate["bound"])
            self.assertEqual(primal, self.analysis.sharp_tail(n, k, p, rho))

    def test_dependence_free_bound_is_attained_and_contains_other_models(self):
        for n, k, p in product((3, 8), (1, 2, 3), (F(0), F(1, 100), F(1, 5), F(1))):
            free = min(F(1), n * p / k)
            self.assertEqual(self.analysis.markov_tail(n, k, p), free)
            self.assertEqual(self.analysis.sharp_tail(n, k, p, F(1)), free)
            for rho in (F(0), F(1, 100), F(1, 3), F(1)):
                point = self.analysis.beta_tail(n, k, p, rho)
                envelope = self.analysis.sharp_tail(n, k, p, rho)
                self.assertLessEqual(point, envelope)
                self.assertLessEqual(envelope, free)

    def test_marginal_cap_augmentation_is_feasible(self):
        n, low, high, cap = 5, F(1, 5), F(2, 5), F(1, 10)
        flip = (high - low) / (1 - low)
        for law in feasible_vertices(n, low, cap):
            augmented = [F(0)] * (n + 1)
            for s, probability in law.items():
                extras = bernoulli_convolution(n - s, flip)
                for extra, probability_extra in enumerate(extras):
                    augmented[s + extra] += probability * probability_extra
            self.assertEqual(sum(augmented), 1)
            self.assertEqual(sum(s * q for s, q in enumerate(augmented)), n * high)
            factorial = sum(s * (s - 1) * q for s, q in enumerate(augmented))
            self.assertLessEqual(factorial, n * (n - 1) * (high * high + cap * high * (1 - high)))
            for k in range(1, n + 1):
                self.assertGreaterEqual(sum(augmented[k:]), sum(q for s, q in law.items() if s >= k))

    def test_holm_matches_closed_testing_with_ties_and_reordering(self):
        for values in product((F(0), F(1, 100), F(1, 60), F(1, 40), F(1, 20), F(1, 2), F(1)), repeat=3):
            raw = dict(zip(("A", "B", "C"), values))
            adjusted = self.analysis.holm(raw)
            self.assertEqual(adjusted, holm_closed(raw))
            self.assertEqual(any(p <= LEVEL for p in adjusted.values()), 3 * min(values) <= LEVEL)
        self.assertEqual(self.analysis.holm({"A": F(1, 60), "B": F(1), "C": F(1)})["A"], LEVEL)

    def test_worst_case_tail_pvalues_are_superuniform(self):
        n, p, cap = 5, F(1, 5), F(1, 10)
        pvalues = [F(1)] + [self.analysis.sharp_tail(n, k, p, cap) for k in range(1, n + 1)]
        for law in feasible_vertices(n, p, cap):
            for threshold in set(pvalues + [F(1, 20), F(1, 60)]):
                probability = sum(q for s, q in law.items() if pvalues[s] <= threshold)
                self.assertLessEqual(probability, threshold)

    def test_zero_pairwise_correlation_does_not_imply_independence(self):
        n, k, p = 100, 10, F(1, 100)
        certificate = self.analysis.sharp_certificate(n, k, p, F(0))
        law = certificate["law"]
        factorial = sum(s * (s - 1) * q for s, q in law.items())
        self.assertEqual(factorial, n * (n - 1) * p * p)
        self.assertGreater(certificate["bound"], self.analysis.binomial_tail(n, k, p))
        self.assertLess(len([q for q in law.values() if q]), n + 1)

    def test_beta_dependence_is_not_a_monotone_robustness_parameter(self):
        p, n, k = F(1, 100), 100, 10
        at_zero = self.analysis.beta_tail(n, k, p, F(0))
        interior = self.analysis.beta_tail(n, k, p, F(1, 10))
        at_one = self.analysis.beta_tail(n, k, p, F(1))
        self.assertGreater(interior, at_zero)
        self.assertGreater(interior, at_one)
        self.assertEqual(at_one, p)
        self.assertTrue(self.analysis.scenario("beta", p, F(1))["all"])
        self.assertFalse(self.analysis.scenario("exchangeable", p, F(1))["any"])

    def test_line_specific_caps_and_holm_family_are_not_pooled(self):
        rows = self.analysis.load_data()
        caps = {"CHP-212": (F(1, 100), F(0)), "LA-N-5": (F(1, 50), F(1, 100)),
                "MC-IX": (F(1, 200), F(1, 10))}
        raw = {row["line"]: self.analysis.sharp_tail(row["n"], row["k"], *caps[row["line"]])
               for row in rows}
        self.assertEqual(self.analysis.holm(raw), holm_closed(raw))
        for model in ("independent", "beta", "exchangeable", "markov"):
            result = self.analysis.scenario(model, F(1, 100), F(1, 100))
            self.assertEqual(set(result["raw"]), set(caps))
            self.assertEqual(result["global"], min(F(1), 3 * min(result["raw"].values())))
            self.assertEqual(result["adjusted"], holm_closed(result["raw"]))
            self.assertEqual(result["any"], any(p <= LEVEL for p in result["adjusted"].values()))
            self.assertEqual(result["all"], all(p <= LEVEL for p in result["adjusted"].values()))

    def test_monotone_size_and_robust_correlation_caps(self):
        for model in ("independent", "beta", "exchangeable", "markov"):
            for rho in (F(0), F(1, 1000), F(1, 10), F(1)):
                prior = {line: F(0) for line in ("CHP-212", "LA-N-5", "MC-IX")}
                for p in (F(0), F(1, 1000), F(1, 100), F(1, 20), F(1)):
                    current = self.analysis.scenario(model, p, rho)["adjusted"]
                    for line in prior:
                        self.assertLessEqual(prior[line], current[line])
                    prior = current
        for p in (F(1, 1000), F(1, 100), F(1, 20)):
            previous = F(0)
            for rho in (F(0), F(1, 1000), F(1, 100), F(1, 10), F(1)):
                current = self.analysis.scenario("exchangeable", p, rho)["global"]
                self.assertLessEqual(previous, current)
                previous = current

    def test_every_reported_size_frontier_has_exact_decision_brackets(self):
        targets = ("CHP-212", "LA-N-5", "MC-IX", "any")
        table = [(model, row["rho"], row["boundary"])
                 for model, rows in self.result["frontiers"].items() for row in rows]
        table += [("independent", F(0), self.result["independent_boundary"]),
                  ("markov", F(0), self.result["markov_boundary"])]
        for model, rho, frontiers in table:
            self.assertEqual(set(frontiers), set(targets))
            for target, bracket in frontiers.items():
                lo, hi = bracket
                self.assertLessEqual(F(0), lo)
                self.assertLess(lo, hi)
                self.assertLessEqual(hi, F(1))
                self.assertLessEqual(hi - lo, F(1, 10**12))
                low = self.analysis.scenario(model, lo, rho)
                high = self.analysis.scenario(model, hi, rho)
                low_value = low["global"] if target == "any" else low["adjusted"][target]
                high_value = high["global"] if target == "any" else high["adjusted"][target]
                self.assertLessEqual(low_value, LEVEL, (model, rho, target))
                self.assertGreater(high_value, LEVEL, (model, rho, target))
                self.assertEqual(tuple(bracket), self.analysis.boundary(model, rho, target))
                # Printed decimal endpoints must enclose, rather than round
                # inward from, the rational computation.
                rendered = self.analysis.bracket_text(bracket)
                shown = tuple(map(F, rendered.strip("[]").split(", ")))
                self.assertLessEqual(shown[0], lo)
                self.assertGreaterEqual(shown[1], hi)
        self.assertEqual(tuple(row["rho"] for row in self.result["frontiers"]["beta"]),
                         self.analysis.SLICES)
        self.assertEqual(tuple(row["rho"] for row in self.result["frontiers"]["exchangeable"]),
                         self.analysis.SLICES)

    def test_exact_nominal_robust_frontier_and_equality_convention(self):
        transitions = self.result["robust_crossings"]
        self.assertEqual(len(transitions), 1)
        transition = transitions[0]
        self.assertEqual(transition["rho"], F(17, 3267))
        self.assertEqual(tuple(transition["left"]), ("LA-N-5",))
        self.assertEqual(tuple(transition["at"]), ("LA-N-5",))
        self.assertEqual(tuple(transition["right"]), ())
        exact = self.analysis.scenario("exchangeable", F(1, 100), transition["rho"])
        self.assertEqual(exact["adjusted"]["LA-N-5"], LEVEL)
        self.assertTrue(exact["any"])
        above = self.analysis.scenario("exchangeable", F(1, 100), transition["rho"] + F(1, 10**12))
        self.assertFalse(above["any"])
        for name, cap in (("uncorrelated", F(0)), ("nominal_boundary", transition["rho"])):
            witness = self.result["witnesses"][name]
            self.verifier.verify_moment_law(100, 10, F(1, 100), cap,
                                            witness["bound"], witness["law"])
        self.assertEqual(self.result["witnesses"]["uncorrelated"]["bound"], F(11, 1000))

    def test_joint_closed_form_against_independent_primal_laws_and_holm(self):
        examples = self.result["joint_examples"]
        self.assertEqual(tuple(row["A"] for row in examples),
                         (F(1, 100), F(23, 2000), F(3, 250), F(7, 500)))
        sizes = {row["A"] for row in examples}
        sizes.update((F(23, 2000) - F(1, 10**8), F(23, 2000) + F(1, 10**8)))
        for size in sorted(sizes):
            correlation = self.analysis.nominal_or_worse_boundary(size)
            self.assertIsNotNone(correlation)
            for cap, reject in ((correlation / 2, True), (correlation, True),
                                (correlation + F(1, 10**10), False)):
                raw = {}
                for row in self.analysis.load_data():
                    n, k = row["n"], row["k"]
                    # Enumerate all low-support primal vertices independently.
                    # This supplies an attaining lower bound. The separately
                    # checked polynomial majorant holds on ALL 0,...,n and
                    # supplies the matching upper bound, so no count cap is
                    # mistaken for a proof about the original model.
                    laws = feasible_vertices(n, size, cap, support_limit=k)
                    best = max(sum(q for s, q in law.items() if s >= k) for law in laws)
                    certificate = self.analysis.sharp_certificate(n, k, size, cap)
                    c0, c1, c2 = certificate["dual"]
                    self.assertGreaterEqual(c2, 0)
                    for s in range(n + 1):
                        self.assertGreaterEqual(c0 + c1*s + c2*s*(s-1), F(s >= k))
                    ceiling = c0 + c1*n*size + c2*n*(n-1)*(size*size + cap*size*(1-size))
                    self.assertEqual(best, ceiling)
                    raw[row["line"]] = best
                adjusted = holm_closed(raw)
                self.assertEqual(any(value <= LEVEL for value in adjusted.values()), reject)
                if cap == correlation:
                    self.assertEqual(adjusted["LA-N-5"], LEVEL)
            for row in examples:
                if row["A"] == size:
                    self.assertEqual(row["R_boundary"], correlation)
        tie_size = F(23, 2000)
        tie_cap = self.analysis.nominal_or_worse_boundary(tie_size)
        laws = feasible_vertices(100, tie_size, tie_cap, support_limit=10)
        attaining = [{s: q for s, q in law.items() if q}
                     for law in laws if law.get(10, F(0)) == F(1, 60)]
        self.assertTrue(attaining)
        self.assertTrue(all(law == {1: F(59, 60), 10: F(1, 60)} for law in attaining))

    def test_algebraic_size_endpoint_and_specialization_domain(self):
        lo, hi = self.result["frontiers"]["exchangeable"][0]["boundary"]["any"]
        # Equivalent algebraic definitions of the upper root; no float sqrt.
        polynomial = lambda a: 9900*a*a - 200*a + F(4, 5)
        self.assertLessEqual(polynomial(lo), 0)
        self.assertGreater(polynomial(hi), 0)
        self.assertGreater(2475*lo - 25, 0)
        self.assertLessEqual((2475*lo - 25)**2, 130)
        self.assertGreater((2475*hi - 25)**2, 130)
        self.assertIsNotNone(self.analysis.nominal_or_worse_boundary(lo))
        self.assertIsNone(self.analysis.nominal_or_worse_boundary(hi))
        for size in (F(1, 50), F(1, 20), F(1)):
            self.assertIsNone(self.analysis.nominal_or_worse_boundary(size))
            self.assertFalse(self.analysis.scenario("exchangeable", size, F(0))["any"])
        with self.assertRaises(ValueError):
            self.analysis.nominal_or_worse_boundary(F(9, 1000))

    def test_calibration_budget_and_independent_minimal_design_count(self):
        budget = self.result["calibration_budget"]
        self.assertEqual(budget["epsilon"], F(1, 100))
        self.assertEqual(budget["test_level"], F(1, 25))
        self.assertEqual(budget["epsilon"] + budget["test_level"], LEVEL)
        design = budget["design"]
        self.assertEqual((design["n"], design["A"], design["R"], design["k"]),
                         (100, F(1, 100), F(1, 200), 12))
        self.assertEqual(design["threshold"], budget["test_level"] / 3)
        for k in range(1, design["k"] + 1):
            laws = feasible_vertices(design["n"], design["A"], design["R"], support_limit=k)
            optimum = max(sum(q for s, q in law.items() if s >= k) for law in laws)
            certificate = self.analysis.sharp_certificate(design["n"], k, design["A"], design["R"])
            c0, c1, c2 = certificate["dual"]
            self.assertGreaterEqual(c2, 0)
            for s in range(design["n"] + 1):
                self.assertGreaterEqual(c0+c1*s+c2*s*(s-1), F(s >= k))
            self.assertEqual(optimum, c0+c1*certificate["mean"]+c2*certificate["second_cap"])
            if k == design["k"]:
                self.assertEqual(optimum, design["tail"])
                self.assertLessEqual(optimum, design["threshold"])
            else:
                self.assertGreater(optimum, design["threshold"])
                if k == design["k"] - 1:
                    self.assertEqual(optimum, design["previous_tail"])

    def test_beta_polynomials_equal_probability_differences(self):
        # A polynomial sign certificate is useful only if its algebra agrees
        # with the original probability law, especially after canceling 1-rho.
        for n, k in ((1, 1), (7, 1), (7, 3), (7, 7), (50, 1), (100, 2), (100, 10)):
            for threshold in (F(1, 20), F(1, 40), F(1, 60)):
                polynomial = self.verifier.beta_tail_polynomial(n, k, F(1, 100), threshold)
                for rho in (F(0), F(1, 1000), F(1, 10), F(3, 4), F(1)):
                    numerator = self.verifier.evaluate_polynomial(polynomial.coefficients, rho)
                    denominator = self.verifier.evaluate_polynomial(polynomial.denominator, rho)
                    self.assertGreater(denominator, 0)
                    self.assertEqual(polynomial.scale * numerator / denominator,
                                     self.analysis.beta_tail(n, k, F(1, 100), rho) - threshold)

    def test_complete_root_isolator_on_known_polynomials_and_failure(self):
        polynomials = [((0, -1, 1), (F(0), F(1))),
                       ((1, -4, 4), (F(1, 2),)),
                       ((24, -250, 875, -1250, 625),
                        (F(1, 5), F(2, 5), F(3, 5), F(4, 5))),
                       ((1, 0, 1), ())]
        for polynomial, expected in polynomials:
            certificate = self.verifier.isolate_unit_roots_certificate(polynomial)
            self.assertEqual(len(certificate.roots), len(expected))
            self.assertEqual(certificate.zero_variation_leaves + certificate.one_variation_leaves,
                             certificate.splits + 1)
            for (lo, hi), root in zip(certificate.roots, expected):
                self.assertLessEqual(lo, root)
                self.assertLessEqual(root, hi)
                self.assertLessEqual(hi - lo, F(1, 10**13))
        with self.assertRaises(ValueError):
            self.verifier.isolate_unit_roots((0,))
        # A non-dyadic repeated root cannot acquire a one-variation leaf;
        # the bounded algorithm must stop honestly instead of missing it.
        with self.assertRaises(ArithmeticError):
            self.verifier.isolate_unit_roots((1, -6, 9))

    def test_bernstein_conversion_and_subdivision_preserve_polynomial(self):
        def bernstein(coefficients, x):
            degree = len(coefficients) - 1
            return sum(F(math.comb(degree, j)) * coefficient * x**j * (1-x)**(degree-j)
                       for j, coefficient in enumerate(coefficients))

        for polynomial in ((2, -5, 0, 3), (1, -7, 2, 9), (24, -250, 875, -1250, 625)):
            coefficients = self.verifier._power_to_bernstein(polynomial)
            coefficients += (0,) * (len(polynomial) - len(coefficients))
            anchor = F(1, 3)
            scale = bernstein(coefficients, anchor) / self.verifier.evaluate_polynomial(polynomial, anchor)
            self.assertGreater(scale, 0)
            for x in (F(0), F(1, 10), F(1, 2), F(4, 5), F(1)):
                self.assertEqual(bernstein(coefficients, x),
                                 scale * self.verifier.evaluate_polynomial(polynomial, x))
            left, right = self.verifier._bisect_bernstein(coefficients)
            for child, offset in ((left, F(0)), (right, F(1, 2))):
                child_scale = bernstein(child, anchor) / bernstein(coefficients, offset + anchor / 2)
                self.assertGreater(child_scale, 0)
                for x in (F(0), F(1, 10), F(1, 2), F(4, 5), F(1)):
                    self.assertEqual(bernstein(child, x),
                                     child_scale * bernstein(coefficients, offset + x / 2))

    def test_false_moment_certificates_are_rejected(self):
        witness = self.analysis.sharp_certificate(100, 10, F(1, 100), F(0))
        with self.assertRaises(AssertionError):
            self.verifier.verify_moment_law(100, 10, F(1, 100), F(0), F(1), witness["law"])
        with self.assertRaises(AssertionError):
            self.verifier.verify_moment_law(100, 10, F(1, 100), F(0), F(0), {0: F(1)})
        with self.assertRaises(AssertionError):
            self.verifier.verify_moment_law(100, 10, F(1, 100), F(0), F(1), {101: F(1)})

    def test_all_nominal_beta_crossings_and_regimes(self):
        result = self.result["beta_crossings"]
        observations = {row["line"]: (row["n"], row["k"]) for row in self.analysis.load_data()}
        certificates = result["polynomial_certificates"]
        self.assertEqual(len(certificates), 9)
        self.assertEqual(sum(row["root_count"] for row in certificates), len(result["raw_roots"]))
        self.assertEqual(len(result["raw_roots"]), 10)
        self.assertEqual(len(result["transitions"]), 2)
        for row in certificates:
            self.assertEqual(row["zero_variation_leaves"] + row["one_variation_leaves"], row["splits"] + 1)
            self.assertEqual(row["root_count"], row["one_variation_leaves"] + row["exact_roots"])
        for root in result["raw_roots"]:
            n, k = observations[root["line"]]
            lo, hi = root["lo"], root["hi"]
            self.assertLessEqual(hi - lo, F(1, 10**12))
            left = self.analysis.beta_tail(n, k, F(1, 100), lo) - LEVEL / root["divisor"]
            right = self.analysis.beta_tail(n, k, F(1, 100), hi) - LEVEL / root["divisor"]
            self.assertLess(left * right, 0)
        for regime in result["regimes"]:
            middle = (regime["lo"] + regime["hi"]) / 2
            observed = self.analysis.scenario("beta", F(1, 100), middle)["rejected"]
            self.assertEqual(tuple(sorted(observed)), tuple(regime["rejected"]))
        self.assertEqual(tuple(result["at_zero"]), ("LA-N-5",))
        self.assertEqual(set(result["at_one"]), set(observations))
        for root in result["transitions"]:
            before = self.analysis.scenario("beta", F(1, 100), root["lo"])["rejected"]
            after = self.analysis.scenario("beta", F(1, 100), root["hi"])["rejected"]
            self.assertEqual(tuple(sorted(before)), tuple(root["before"]))
            self.assertEqual(tuple(sorted(after)), tuple(root["after"]))
            self.assertEqual(set(root["at_boundary"]), set(before) | set(after))

    def test_experiment_counts_are_minimal_and_not_observed_power(self):
        for row in self.result["experiments"]:
            self.assertEqual(row["n"], 100)
            self.assertEqual(row["threshold"], LEVEL / 3)
            self.assertIsNotNone(row["k"])
            self.assertLessEqual(row["tail"], row["threshold"])
            self.assertGreater(row["previous_tail"], row["threshold"])
            if row["model"] == "independent":
                n, k, p = row["n"], row["k"], row["A"]
                independent = sum(F(math.comb(n, s)) * p**s * (1 - p)**(n - s)
                                  for s in range(k, n + 1))
                self.assertEqual(row["tail"], independent)
            elif row["model"] == "markov":
                self.assertEqual(row["tail"], min(F(1), row["n"] * row["A"] / row["k"]))
            else:
                certificate = self.analysis.sharp_certificate(row["n"], row["k"], row["A"], row["R"])
                self.verifier.verify_moment_law(row["n"], row["k"], row["A"], row["R"],
                                                row["tail"], certificate["law"])

    def test_exact_input_validation_and_tail_boundaries(self):
        for function, extra in ((self.analysis.binomial_tail, ()),
                                (self.analysis.beta_tail, (F(1, 10),)),
                                (self.analysis.sharp_tail, (F(1, 10),)),
                                (self.analysis.markov_tail, ())):
            self.assertEqual(function(5, 0, F(1, 5), *extra), 1)
            self.assertEqual(function(5, 6, F(1, 5), *extra), 0)
            for n, k in ((0, 1), (True, 1), (5, True), (F(5), 1), (5, F(1))):
                with self.assertRaises(ValueError):
                    function(n, k, F(1, 5), *extra)
            for invalid in (0.1, True, False):
                with self.assertRaises(TypeError):
                    function(5, 1, invalid, *extra)
            for invalid in (F(-1, 100), F(101, 100)):
                with self.assertRaises(ValueError):
                    function(5, 1, invalid, *extra)
        for function in (self.analysis.beta_tail, self.analysis.sharp_tail):
            for invalid in (0.1, True):
                with self.assertRaises(TypeError):
                    function(5, 1, F(1, 5), invalid)
            for invalid in (F(-1, 100), F(101, 100)):
                with self.assertRaises(ValueError):
                    function(5, 1, F(1, 5), invalid)
        with self.assertRaises(ValueError):
            self.analysis.holm({})
        with self.assertRaises(ValueError):
            self.analysis.scenario("unrecognized", F(1, 100), F(0))

    def test_stdlib_only_no_runtime_retrieval_and_nonclinical_labels(self):
        forbidden = {"urllib", "http", "socket", "requests", "subprocess"}
        for filename in ("analysis.py", "verifier.py"):
            tree = ast.parse((HERE / filename).read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    names = [alias.name.split(".")[0] for alias in node.names]
                elif isinstance(node, ast.ImportFrom):
                    names = [(node.module or "").split(".")[0]]
                else:
                    continue
                for name in names:
                    self.assertIn(name, sys.stdlib_module_names)
                    self.assertNotIn(name, forbidden)
        for filename in ("README.md", "PROTOCOL.md"):
            text = (HERE / filename).read_text(encoding="utf-8").lower()
            self.assertIn("nonclinical", text)
            self.assertIn("no treatment advice", text)


if __name__ == "__main__":
    unittest.main()
