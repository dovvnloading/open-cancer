"""Offline, nonclinical verification of the frozen separability audit.

Decimal-based reference calculations below read the cited CSVs independently.
Synthetic controls are arithmetic fixtures, not new biological observations.
Passing these tests establishes neither source truth nor clinical validity.
"""

import csv
from datetime import date
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import itertools
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research" / "mutation-clock-separability"
SPEC = importlib.util.spec_from_file_location(
    "mutation_clock_separability", RESEARCH / "analysis.py"
)
ANALYSIS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ANALYSIS
SPEC.loader.exec_module(ANALYSIS)

CSV_NAMES = ("blood_rates.csv", "colon_rates.csv", "neuron_rates.csv", "assay_rates.csv")
INPUT_PATHS = (Path("PROTOCOL.md"),) + tuple(Path("data") / name for name in CSV_NAMES)
PROTOCOL_SHA256 = "d3dc484ae6ee879c09dcdfc209311fc6ef66ebf2fa49350f1b1cdb8148b18191"
CONSISTENT = "Separable clock consistent"
REJECTED = "Rejected under the stated measurement assumptions"
UNTESTABLE = "Untestable"
PROTOCOL_COMMIT = "3dea146499a321f6fed5e175b2610c9884ac89ca"


def read_rows(base=RESEARCH):
    """Independent source reader; preserve the actual printed decimal strings."""
    rows = {}
    for name in CSV_NAMES:
        with (Path(base) / "data" / name).open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if row["id"] in rows:
                    raise AssertionError("duplicate source ID: " + row["id"])
                rows[row["id"]] = row
    return rows


def number(row, key="estimate"):
    return Decimal(row[key])


def reference_hash(path):
    # Universal-newline decoding provides an independent normalization route.
    return hashlib.sha256(Path(path).read_text(encoding="utf-8").encode("utf-8")).hexdigest()


def reference_departures(q):
    with localcontext() as context:
        context.prec = 60
        q = Decimal(str(q))
        fold = max(q, Decimal(1) / q)
        return {
            "Q": float(q),
            "log_Q": float(q.ln()),
            "D": float(fold),
            "delta": float(fold.sqrt()),
            "rho": float(fold.sqrt().sqrt()),
        }


def reference_point(bm, bh, cm, ch):
    with localcontext() as context:
        context.prec = 60
        bm, bh, cm, ch = (Decimal(str(x)) for x in (bm, bh, cm, ch))
        result = reference_departures((bm / bh) / (cm / ch))
        result.update(blood_or_neuron_ratio=float(bm / bh), colon_ratio=float(cm / ch))
        return result


def reference_compatibility(lower, upper):
    with localcontext() as context:
        context.prec = 60
        lower, upper = Decimal(str(lower)), Decimal(str(upper))
        inside = lower <= 1 <= upper
        nearest_fold = Decimal(1) if inside else (lower if lower > 1 else 1 / upper)
        return {
            "q_low": float(lower), "q_high": float(upper),
            "log_distance_min": float(nearest_fold.ln()),
            "D_min": float(nearest_fold),
            "delta_min": float(nearest_fold.sqrt()),
            "rho_min": float(nearest_fold.sqrt().sqrt()),
            "contains_one": inside,
            "verdict": CONSISTENT if inside else REJECTED,
        }


def reference_box(boxes):
    """Enumerate all corners, independent of production monotonicity code."""
    with localcontext() as context:
        context.prec = 60
        corners = []
        for bm, bh, cm, ch in itertools.product(*boxes):
            corners.append((bm * ch) / (bh * cm))
        return reference_compatibility(min(corners), max(corners))


def printed_box(row, rounding=False):
    if row["ci_lower"] and row["ci_upper"]:
        lower, upper = number(row, "ci_lower"), number(row, "ci_upper")
    else:
        lower = upper = number(row)
    if rounding:
        half_unit = number(row, "printed_resolution") / 2
        lower, upper = lower - half_unit, upper + half_unit
    return lower, upper


def sem_box(row):
    estimate = number(row)
    sem = number(row, "sem") if row["sem"] else Decimal(0)
    return estimate - sem, estimate + sem


def reference_calibration(nh, nm, ch, cm, cf):
    with localcontext() as context:
        context.prec = 60
        nh, nm, ch, cm, cf = (Decimal(str(x)) for x in (nh, nm, ch, cm, cf))
        hr, mr = nh / ch, nm / cm
        multiplier = (hr * mr).sqrt()
        alpha = multiplier.ln()
        rh, rm = hr.ln() - alpha, mr.ln() - alpha
        return {
            "human_ratio": float(hr), "mouse_ratio": float(mr),
            "alpha": float(alpha), "multiplier": float(multiplier),
            "training_log_residuals": {"human": float(rh), "mouse": float(rm)},
            "training_log_rmse": float(((rh * rh + rm * rm) / 2).sqrt()),
            "predicted_ferret": float(cf * multiplier),
        }


def reference_score(observed, predicted):
    with localcontext() as context:
        context.prec = 60
        observed, predicted = Decimal(str(observed)), Decimal(str(predicted))
        ratio = observed / predicted
        return {
            "observed": float(observed), "predicted": float(predicted),
            "observed_over_predicted": float(ratio),
            "signed_log_error": float(ratio.ln()),
            "absolute_fold_error": float(max(ratio, 1 / ratio)),
        }


class MutationClockSeparabilityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = read_rows()
        # A network call cannot silently turn the offline fixture into a live test.
        with patch("socket.create_connection", side_effect=AssertionError("network forbidden")):
            cls.report = ANALYSIS.analyze()

    def assert_recursive_close(self, actual, expected, path="result"):
        """Check keys as well as values, so new reported numbers need an oracle."""
        if isinstance(expected, dict):
            self.assertEqual(set(actual), set(expected), path)
            for key, value in expected.items():
                self.assert_recursive_close(actual[key], value, path + "." + key)
        elif isinstance(expected, (list, tuple)):
            self.assertEqual(len(actual), len(expected), path)
            for i, (left, right) in enumerate(zip(actual, expected)):
                self.assert_recursive_close(left, right, path + "." + str(i))
        elif isinstance(expected, bool) or expected is None or isinstance(expected, str):
            self.assertEqual(actual, expected, path)
        else:
            self.assertTrue(math.isfinite(actual), path)
            self.assertTrue(math.isclose(actual, expected, rel_tol=3e-12, abs_tol=2e-14),
                            (path, actual, expected))

    def test_protocol_is_immutable_and_file_hashes_are_portable(self):
        self.assertEqual(reference_hash(RESEARCH / "PROTOCOL.md"), PROTOCOL_SHA256)
        # Publishing through the GitHub connector creates a new commit identity.
        # This is the original local pre-analysis Git object, also retained in
        # the public protocol-publication commit message. No Git/network needed.
        freeze_payload = (
            "tree f5ddac4c38de1e9672cddf568d5b0160c2d10e96\n"
            "parent 09d916cdea7053f25643c6d0bd9c45827436e315\n"
            "author Codex <codex@openai.com> 1794989366 -0500\n"
            "committer Codex <codex@openai.com> 1794989366 -0500\n"
            "\nFreeze mutation-clock separability protocol before calculations\n"
        ).encode("utf-8")
        git_header = ("commit " + str(len(freeze_payload)) + "\0").encode("ascii")
        self.assertEqual(hashlib.sha1(git_header + freeze_payload).hexdigest(), PROTOCOL_COMMIT)
        expected = {path.as_posix(): reference_hash(RESEARCH / path) for path in INPUT_PATHS}
        self.assertEqual(ANALYSIS.input_hashes(), expected)
        self.assertEqual(self.report["hashes"], expected)
        for key in self.report["hashes"]:
            self.assertEqual(key, Path(key).as_posix())
            self.assertNotIn("\\", key)
        protocol = (RESEARCH / "PROTOCOL.md").read_text(encoding="utf-8")
        self.assertIn("before executing any empirical arithmetic", protocol)
        self.assertIn("No browser file downloads", protocol)

    def test_hashes_and_analysis_survive_lf_crlf_and_lone_cr(self):
        original = self.report
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            (base / "data").mkdir()
            for ending in (b"\n", b"\r\n", b"\r"):
                for relative in INPUT_PATHS:
                    text = (RESEARCH / relative).read_text(encoding="utf-8")
                    (base / relative).write_bytes(text.replace("\n", ending.decode("ascii")).encode("utf-8"))
                self.assertEqual(ANALYSIS.input_hashes(base), original["hashes"])
                self.assert_recursive_close(ANALYSIS.analyze(base), original)
        fixture = "nonclinical α\nsecond line\n".encode("utf-8")
        for ending in (b"\n", b"\r\n", b"\r"):
            self.assertEqual(ANALYSIS.normalized_bytes(fixture.replace(b"\n", ending)), fixture)
        with self.assertRaises(UnicodeDecodeError):
            ANALYSIS.normalized_bytes(b"\xff")

    def test_every_input_has_provenance_and_distinct_units_and_uncertainty(self):
        expected_ids = {
            "mouse_hsc", "human_hsc", "human_colon", "mouse_colon", "ferret_colon_mean",
            "ferret_colon_slope", "human_neuron", "mouse_neuron", "ferret_neuron",
            "mouse_hsc_per_base", "mouse_phylogeny_per_base", "mouse_duplex_per_base",
        }
        self.assertEqual(set(self.rows), expected_ids)
        required = (
            "id", "species", "source_taxon", "tissue", "estimate", "unit", "original_unit",
            "estimator", "assay", "genotype", "printed_resolution", "uncertainty_status",
            "approximate", "role", "source", "source_url", "version", "locator",
            "access_date", "reuse_basis", "reuse_terms_url", "notes",
        )
        for row in self.rows.values():
            for field in required:
                self.assertTrue(row[field].strip(), (row["id"], field))
            for field in ("source_url", "reuse_terms_url"):
                self.assertTrue(row[field].startswith("https://"), (row["id"], field))
            self.assertEqual(date.fromisoformat(row["access_date"]), date(2026, 11, 18))
            self.assertGreater(number(row), 0)
            self.assertGreater(number(row, "printed_resolution"), 0)
            self.assertEqual(bool(row["ci_lower"]), bool(row["ci_upper"]))
            self.assertFalse(row["sem"] and row["ci_lower"], row["id"])
            if row["ci_lower"]:
                self.assertLessEqual(number(row, "ci_lower"), number(row))
                self.assertGreaterEqual(number(row, "ci_upper"), number(row))
            if row["sem"]:
                self.assertGreater(number(row, "sem"), 0)
        for key in ("mouse_hsc", "human_hsc", "human_colon", "mouse_colon",
                    "ferret_colon_mean", "ferret_colon_slope"):
            self.assertEqual(self.rows[key]["unit"], "SBS_per_cell_or_genome_year")
        for key in ("human_neuron", "mouse_neuron", "ferret_neuron"):
            self.assertEqual(self.rows[key]["unit"], "SNV_per_genome_year")
            self.assertFalse(self.rows[key]["ci_lower"])
            self.assertFalse(self.rows[key]["ci_upper"])
        for key in ("mouse_hsc", "human_hsc"):
            self.assertEqual(self.rows[key]["ci_level"], "95")
        for key in ("human_colon", "mouse_colon", "ferret_colon_mean", "ferret_colon_slope"):
            self.assertFalse(self.rows[key]["ci_lower"])
            self.assertFalse(self.rows[key]["sem"])
        self.assertFalse(self.rows["human_neuron"]["sem"])
        for key in ("mouse_neuron", "ferret_neuron"):
            self.assertTrue(self.rows[key]["sem"])
        self.assertEqual(self.rows["mouse_hsc_per_base"]["unit"], "SBS_per_bp_year")
        for key in ("mouse_phylogeny_per_base", "mouse_duplex_per_base"):
            self.assertEqual(self.rows[key]["unit"], "SBS_per_bp_year_interpreted")
            self.assertIn("bp per year", self.rows[key]["original_unit"])
        self.assertNotEqual(self.rows["ferret_colon_mean"]["estimate"],
                            self.rows["ferret_colon_slope"]["estimate"])
        self.assertIn("audit", self.rows["ferret_colon_slope"]["role"].lower())
        # The frozen IDs survive the audit correction; the actual source tissue
        # must never silently turn the ferret predictor into measured colon.
        for key in ("ferret_colon_mean", "ferret_colon_slope"):
            self.assertEqual(self.rows[key]["tissue"], "small_intestinal_crypt")
            self.assertIn("3b4db10c7c38c63aeebfc9a48ce4bfa376dfd2bc", self.rows[key]["source_url"])
            self.assertIn("3b4db10c7c38c63aeebfc9a48ce4bfa376dfd2bc", self.rows[key]["version"])
        for key in ("human_neuron", "mouse_neuron", "ferret_neuron"):
            self.assertIn("Version 1", self.rows[key]["version"])
            self.assertIn("not peer reviewed", self.rows[key]["version"])
        self.assertNotEqual(self.rows["mouse_neuron"]["genotype"], self.rows["mouse_hsc"]["genotype"])

    def test_reported_inputs_are_exact_csv_transcriptions(self):
        self.assertEqual(set(self.report["inputs"]), set(self.rows))
        for key, row in self.rows.items():
            actual = self.report["inputs"][key]
            expected = {field: row[field] for field in (
                "species", "source_taxon", "tissue", "unit", "original_unit", "estimator",
                "assay", "genotype", "uncertainty_status", "approximate", "role", "source",
                "source_url", "version", "locator", "access_date", "reuse_basis",
                "reuse_terms_url", "notes",
            )}
            for field in ("estimate", "ci_lower", "ci_upper", "ci_level", "sem", "printed_resolution"):
                expected[field] = float(row[field]) if row[field] else None
                expected[field + "_printed"] = row[field]
            self.assertEqual(actual, expected, key)

    def test_every_primary_and_assay_number_independently(self):
        rows = self.rows
        ordered = [rows[key] for key in ("mouse_hsc", "human_hsc", "mouse_colon", "human_colon")]
        point = reference_point(*(number(row) for row in ordered))
        box = reference_box([printed_box(row) for row in ordered])
        rounding = reference_box([printed_box(row, rounding=True) for row in ordered])
        self.assert_recursive_close(self.report["primary"], {
            "point": point, "printed_interval_box": box, "rounding_box": rounding,
            "fixed_missing_intervals": ["mouse_colon", "human_colon"],
            "unconditional_verdict": UNTESTABLE,
        })
        duplex, phylogeny = rows["mouse_duplex_per_base"], rows["mouse_phylogeny_per_base"]
        with localcontext() as context:
            context.prec = 60
            assay_ratio = number(duplex) / number(phylogeny)
            ratio_low = number(duplex, "ci_lower") / number(phylogeny, "ci_upper")
            ratio_high = number(duplex, "ci_upper") / number(phylogeny, "ci_lower")
            observed_q = number(ordered[0]) * number(ordered[3]) / (number(ordered[1]) * number(ordered[2]))
            q_low = (number(ordered[0], "ci_lower") * number(ordered[3]) /
                     (number(ordered[1], "ci_upper") * number(ordered[2])))
            q_high = (number(ordered[0], "ci_upper") * number(ordered[3]) /
                      (number(ordered[1], "ci_lower") * number(ordered[2])))
            bands = [dict(K=k, **reference_compatibility(q_low / k, q_high * k))
                     for k in (1, 2, 4, 8)]
            expected = {
                "break_even_point": float(observed_q),
                "break_even_box": {"low": float(q_low), "high": float(q_high)},
                "hypothetical_bands": bands,
                "duplex_phylogenetic_ratio": float(assay_ratio),
                "duplex_phylogenetic_box": {"low": float(ratio_low), "high": float(ratio_high)},
                "bridged_point": reference_departures(observed_q * assay_ratio),
                "bridged_box": reference_compatibility(q_low * ratio_low, q_high * ratio_high),
            }
        self.assert_recursive_close(self.report["assay"], expected)

    def test_every_neuron_calibration_score_and_sem_number_independently(self):
        rows = self.rows
        nh, nm, nf = (number(rows[key]) for key in ("human_neuron", "mouse_neuron", "ferret_neuron"))
        ch, cm, cf = (number(rows[key]) for key in ("human_colon", "mouse_colon", "ferret_colon_mean"))
        calibration = reference_calibration(nh, nm, ch, cm, cf)
        hlo, hhi = sem_box(rows["human_neuron"])
        mlo, mhi = sem_box(rows["mouse_neuron"])
        flo, fhi = sem_box(rows["ferret_neuron"])
        with localcontext() as context:
            context.prec = 60
            scale = ((nh / ch) * (nm / cm)).sqrt()
            lower_scale = ((hlo / ch) * (mlo / cm)).sqrt()
            upper_scale = ((hhi / ch) * (mhi / cm)).sqrt()
            lower_prediction, upper_prediction = cf * lower_scale, cf * upper_scale
            transfer = {
                "available": True,
                "alpha_low": float(lower_scale.ln()), "alpha_high": float(upper_scale.ln()),
                "multiplier_low": float(lower_scale), "multiplier_high": float(upper_scale),
                "predicted_low": float(lower_prediction), "predicted_high": float(upper_prediction),
                "observed_low": float(flo), "observed_high": float(fhi),
                "ratio_low": float(flo / upper_prediction), "ratio_high": float(fhi / lower_prediction),
                "intervals_intersect": max(lower_prediction, flo) <= min(upper_prediction, fhi),
            }
            expected = {
                "transfer_scope": {
                    "predictor_tissue": "small_intestinal_crypt",
                    "matched_colon_transfer_verdict": UNTESTABLE,
                    "assumption": "For numerical proxy scores only: ferret small-intestine rate substitutes for unmeasured ferret colon rate; protocol ID names are not tissue claims.",
                },
                "point": reference_point(nm, nh, cm, ch),
                "one_sem_box": dict(available=True, **reference_box([(mlo, mhi), (hlo, hhi), (cm, cm), (ch, ch)])),
                "unpropagated_neuronal_uncertainty": ["human_neuron"],
                "calibration": calibration,
                "ferret_score": reference_score(nf, cf * scale),
                "anchor_predictions": {"human": float(cf * nh / ch), "mouse": float(cf * nm / cm)},
                "one_sem_transfer": transfer,
            }
        self.assert_recursive_close(self.report["neuron"], expected)

    def test_schema_covers_all_reported_numbers(self):
        self.assertEqual(set(self.report), {
            "scope", "protocol_freeze_commit", "inputs", "hashes", "primary", "assay", "neuron",
        })
        self.assertEqual(self.report["protocol_freeze_commit"], PROTOCOL_COMMIT)
        self.assertIn("nonclinical", self.report["scope"].lower())
        self.assertIn("no medical advice", self.report["scope"].lower())

    def test_sharp_per_cell_bound_has_constructive_attainment(self):
        rows = self.rows
        examples = [
            (3, 1, 7, 2), (7, 2, 3, 1), (6, 2, 9, 3),
            tuple(float(rows[key]["estimate"]) for key in ("mouse_hsc", "human_hsc", "mouse_colon", "human_colon")),
        ]
        for bm, bh, cm, ch in examples:
            actual = ANALYSIS.ratio_of_ratios(bm, bh, cm, ch)
            self.assert_recursive_close(actual, reference_point(bm, bh, cm, ch))
            z = math.log(bm) - math.log(bh) - math.log(cm) + math.log(ch)
            log_residuals = (z / 4, -z / 4, -z / 4, z / 4)
            residuals = tuple(math.exp(x) for x in log_residuals)
            fitted = tuple(x / u for x, u in zip((bm, bh, cm, ch), residuals))
            self.assertAlmostEqual((fitted[0] * fitted[3]) / (fitted[1] * fitted[2]), 1, places=13)
            self.assertAlmostEqual(max(max(u, 1 / u) for u in residuals), actual["rho"], places=12)
            self.assertAlmostEqual(4 * math.log(actual["rho"]), abs(z), places=12)
            common_ratio = math.sqrt((bm / bh) * (cm / ch))
            allowances = [max(r / common_ratio, common_ratio / r) for r in (bm / bh, cm / ch)]
            self.assertAlmostEqual(max(allowances), actual["delta"], places=12)
            if abs(z) > 1e-12:
                self.assertLess(4 * math.log(actual["rho"] * (1 - 1e-6)), abs(z))

    def test_positive_boxes_and_synthetic_assay_rescue(self):
        cases = [
            ((Decimal(6), Decimal(6)), (Decimal(2), Decimal(2)), (Decimal(9), Decimal(9)), (Decimal(3), Decimal(3))),
            ((Decimal(2), Decimal(4)), (Decimal(1), Decimal(1)), (Decimal(3), Decimal(3)), (Decimal(1), Decimal(1))),
            ((Decimal(1), Decimal(2)), (Decimal(2), Decimal(3)), (Decimal(8), Decimal(9)), (Decimal(1), Decimal(2))),
            ((Decimal(8), Decimal(9)), (Decimal(1), Decimal(2)), (Decimal(1), Decimal(2)), (Decimal(2), Decimal(3))),
        ]
        for boxes in cases:
            floats = [tuple(float(x) for x in box) for box in boxes]
            self.assert_recursive_close(ANALYSIS.ratio_box(*floats), reference_box(boxes))
        for lower, upper in ((1, 1), (.5, 1), (1, 2), (.2, .3), (2, 3)):
            self.assert_recursive_close(ANALYSIS.compatibility(lower, upper), reference_compatibility(lower, upper))
        before = ANALYSIS.compatibility(.2, .3)
        after = ANALYSIS.compatibility(.2 / 4, .3 * 4)
        self.assertEqual(before["verdict"], REJECTED)
        self.assertEqual(after["verdict"], CONSISTENT)
        self.assertEqual(after["rho_min"], 1)
        # A measurement interaction equal to observed Q exactly rescues latent Q.
        self.assertEqual(.25 / .25, 1)

    def test_invalid_rates_intervals_and_sem_are_not_silently_clipped(self):
        for value in (0, -1, float("inf"), float("nan")):
            with self.assertRaises(ValueError):
                ANALYSIS.departures(value)
        for lower, upper in ((2, 1), (0, 1), (-1, 2)):
            with self.assertRaises(ValueError):
                ANALYSIS.compatibility(lower, upper)
        row = dict(self.rows["mouse_neuron"])
        row["sem"] = row["estimate"]
        self.assertIsNone(ANALYSIS.sem_bounds(row))
        row["sem"] = "999"
        self.assertIsNone(ANALYSIS.sem_bounds(row))
        row["sem"] = "-1"
        with self.assertRaises(ValueError):
            ANALYSIS.sem_bounds(row)
        row["sem"] = ""
        self.assertEqual(ANALYSIS.sem_bounds(row), (float(row["estimate"]), float(row["estimate"])))
        broken_interval = dict(self.rows["mouse_hsc"], ci_upper="")
        with self.assertRaises(ValueError):
            ANALYSIS.printed_bounds(broken_interval)
        # Synthetic endpoint failures keep their point score and refuse only
        # the invalid log sensitivity. No clipping to a small positive value.
        for row_id in ("mouse_neuron", "ferret_neuron"):
            changed = {key: dict(value) for key, value in self.rows.items()}
            changed[row_id]["sem"] = changed[row_id]["estimate"]
            with patch.object(ANALYSIS, "load_rates", return_value=changed):
                result = ANALYSIS.analyze()["neuron"]
            self.assertFalse(result["one_sem_transfer"]["available"])
            self.assertIn("Nonpositive", result["one_sem_transfer"]["reason"])
            self.assertEqual(result["calibration"], self.report["neuron"]["calibration"])
            self.assertEqual(result["ferret_score"], self.report["neuron"]["ferret_score"])
            self.assertEqual(result["one_sem_box"]["available"], row_id == "ferret_neuron")
        changed = {key: dict(value) for key, value in self.rows.items()}
        changed["mouse_hsc"]["unit"] = "SBS_per_bp_year"
        with patch.object(ANALYSIS, "load_rates", return_value=changed):
            with self.assertRaises(ValueError):
                ANALYSIS.analyze()

    def test_species_and_tissue_scale_factors_cancel(self):
        values = [float(self.rows[key]["estimate"]) for key in
                  ("mouse_hsc", "human_hsc", "mouse_colon", "human_colon")]
        baseline = ANALYSIS.ratio_of_ratios(*values)
        for mouse, human, blood, colon in ((2.7, 3.1, 1, 1), (1, 1, 7, .25), (2.7, 3.1, 7, .25)):
            bm, bh, cm, ch = values
            changed = ANALYSIS.ratio_of_ratios(bm * mouse * blood, bh * human * blood,
                                             cm * mouse * colon, ch * human * colon)
            self.assert_recursive_close({key: changed[key] for key in ("Q", "log_Q", "D", "delta", "rho")},
                                        {key: baseline[key] for key in ("Q", "log_Q", "D", "delta", "rho")})

    def test_ferret_outcome_and_audit_slope_cannot_influence_calibration(self):
        original = self.report["neuron"]
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            (base / "data").mkdir()
            for relative in INPUT_PATHS:
                (base / relative).write_bytes((RESEARCH / relative).read_bytes())
            changes = {"ferret_neuron": "12345", "ferret_colon_slope": "54321"}
            for filename in ("neuron_rates.csv", "colon_rates.csv"):
                path = base / "data" / filename
                with path.open(encoding="utf-8", newline="") as handle:
                    records = list(csv.DictReader(handle))
                fields = list(records[0])
                for row in records:
                    if row["id"] in changes:
                        row["estimate"] = changes[row["id"]]
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows(records)
            changed = ANALYSIS.analyze(base)["neuron"]
        self.assertEqual(changed["calibration"], original["calibration"])
        self.assertEqual(changed["anchor_predictions"], original["anchor_predictions"])
        self.assertEqual(changed["point"], original["point"])
        self.assertNotEqual(changed["ferret_score"], original["ferret_score"])
        self.assertEqual(changed["one_sem_transfer"]["predicted_low"], original["one_sem_transfer"]["predicted_low"])
        self.assertEqual(changed["one_sem_transfer"]["predicted_high"], original["one_sem_transfer"]["predicted_high"])

    def test_generated_readme_and_nonwriting_cli_check(self):
        readme_path = RESEARCH / "README.md"
        before = readme_path.read_bytes()
        text = before.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
        begin, end = "<!-- BEGIN COMPUTED RESULTS -->", "<!-- END COMPUTED RESULTS -->"
        self.assertEqual(text.count(begin), 1)
        self.assertEqual(text.count(end), 1)
        actual = text[text.index(begin):text.index(end) + len(end)]
        self.assertEqual(actual, ANALYSIS.render_results(self.report).strip())
        self.assertIn("nonclinical", text.lower())
        self.assertIn("no medical advice", text.lower())
        self.assertIn("3dea146", text)
        result = subprocess.run([sys.executable, "-B", str(RESEARCH / "analysis.py"), "--check"],
                                cwd=ROOT, capture_output=True, text=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(readme_path.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
