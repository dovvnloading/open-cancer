"""Offline tests for the nonclinical, unreviewed priority triage and DNA-loss analysis."""
import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "research/expression-loss-sl-priority-triage"
OLD = ROOT / "research/expression-loss-synthetic-lethality"
SPEC = importlib.util.spec_from_file_location("priority_analysis", HERE / "analysis.py")
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)


class PriorityTriageTests(unittest.TestCase):
    def test_scope_and_search_accounting(self):
        with (OLD / "derived/strict_candidate_queue.csv").open(encoding="utf-8", newline="") as stream:
            pairs = {f"{row['A']}<-{row['B']}" for row in csv.DictReader(stream)
                     if not row["proxy_or_known_explanation"].strip()}
        with (HERE / "triage.csv").open(encoding="utf-8", newline="") as stream:
            items = list(csv.DictReader(stream))
        self.assertEqual(len(pairs), 18)
        self.assertEqual(len(items), 20)
        self.assertEqual({row["item"] for row in items}, pairs | {"H1", "H2"})
        for row in items:
            with self.subTest(item=row["item"]):
                self.assertIn(row["verdict"], {"PUBLISHED", "PARTIALLY KNOWN", "NOT FOUND"})
                self.assertEqual(row["search_date"], "2026-11-22")
                self.assertTrue(row["evidence"].strip())
                self.assertTrue(row["queries"].strip())
                self.assertTrue(row["proxy"].strip())
                if row["verdict"] != "NOT FOUND":
                    self.assertTrue(row["citations"].strip())
        indexed = {row["item"]: row["verdict"] for row in items}
        self.assertEqual(indexed["H1"], "PUBLISHED")
        self.assertEqual(indexed["H2"], "PUBLISHED")

    def test_frozen_offline_result_and_diagnostics(self):
        saved = json.loads((HERE / "derived/decisive_test.json").read_text(encoding="utf-8"))
        output = subprocess.check_output([sys.executable, str(HERE / "analysis.py")], timeout=15)
        replay = json.loads(output)
        for report in (saved, replay):
            self.assertEqual(report["n_screen"], 1103)
            self.assertEqual(sum(report["eligibility"].values()), report["n_screen"])
            self.assertEqual(sum(report["exposure"].values()), report["eligibility"]["representative"])
            self.assertEqual((report["matched_loss"], report["matched_retained"]), (6, 202))
            self.assertEqual(sum(s["loss"] for s in report["strata"]), 6)
            self.assertEqual(sum(s["retained"] for s in report["strata"]), 202)
            self.assertEqual(len(report["strata"]), 6)
            self.assertGreater(report["retained_minus_loss"], 0.2051)
            self.assertLess(report["retained_minus_loss"], 0.2053)
            self.assertEqual(report["diagnostics"]["loss"]["STS"], {"positive": 4, "denominator": 6})
            self.assertEqual(report["diagnostics"]["loss"]["PNPLA4"]["positive"], 3)
            self.assertEqual(report["diagnostics"]["loss"]["PUDP_RNA_off"]["positive"], 5)
            self.assertEqual(report["verdict"], "INCONCLUSIVE")
            self.assertEqual(report["reason"], "insufficient exposure")
            self.assertEqual(report["draws"], 0)
            self.assertIsNone(report["p_zero_upper"])
            self.assertIsNone(report["p_margin_lower"])
        paths = (analysis.SCREEN, analysis.EXTRACT, HERE / "PROTOCOL.md")
        hashes = {path.relative_to(ROOT).as_posix(): analysis.normalized_sha256(path) for path in paths}
        self.assertEqual(saved["sha256_lf_utf8"], hashes)
        self.assertEqual(replay["sha256_lf_utf8"], hashes)

    def test_portable_hash_and_fail_closed_join(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "copy.csv"
            raw = analysis.EXTRACT.read_bytes().replace(b"\r\n", b"\n")
            path.write_bytes(raw.replace(b"\n", b"\r\n"))
            self.assertEqual(analysis.normalized_sha256(path), analysis.normalized_sha256(analysis.EXTRACT))
            path.write_bytes(b"ModelID\r\nACH-1\r\n")
            with self.assertRaisesRegex(ValueError, "ModelIDs"):
                analysis.analyze(path)
            path.write_bytes(b"ModelID\r\nACH-1\r\nACH-1\r\n")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                analysis.read_models(path)

    def test_exposure_boundaries_and_tails(self):
        for invalid in (-0.1, float("nan"), float("inf")):
            with self.assertRaisesRegex(ValueError, "dosage"):
                analysis.exposure(invalid)
        for ratio, expected in ((0.20, "loss"), (0.201, "intermediate_or_gain"), (0.75, "retained"),
                                (1.25, "retained"), (1.251, "intermediate_or_gain")):
            self.assertEqual(analysis.exposure(ratio), expected)
        # Two artificial strata exercise both tails, which the real exposure-count gate skips.
        for decrement, expected in ((1.0, "SUPPORTED"), (0.0, "FALSIFIED")):
            stratum = [("loss", -decrement)] * 4 + [("retained", 0.0)] * 4
            upper, lower = analysis.permutation([stratum, stratum], decrement)
            self.assertEqual(analysis.decision(decrement, upper, lower), expected)
        self.assertEqual(analysis.decision(0.19, 0.1, 0.2), "INCONCLUSIVE")

    def test_corrupted_extract_does_not_become_a_negative_result(self):
        rows = list(analysis.read_models(analysis.EXTRACT).values())
        case = next(i for i, row in enumerate(rows) if row["eligibility"] == "representative")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "invalid.csv"
            for field, value in (("eligibility", "typo"), ("PUDP", "NaN"), ("sex", "")):
                with path.open("w", encoding="utf-8", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=rows[0])
                    writer.writeheader()
                    for i, row in enumerate(rows):
                        writer.writerow({**row, field: value} if i == case else row)
                with self.subTest(field=field), self.assertRaises(ValueError):
                    analysis.analyze(path)


if __name__ == "__main__":
    unittest.main()
