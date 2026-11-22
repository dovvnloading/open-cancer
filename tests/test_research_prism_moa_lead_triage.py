"""Offline checks for nonclinical, unreviewed PRISM lead triage."""
import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from datetime import date
from pathlib import Path

STUDY = Path(__file__).resolve().parents[1] / "research" / "prism-moa-lead-triage"
DERIVED = STUDY / "derived"
SPEC = importlib.util.spec_from_file_location("prism_lead_triage", STUDY / "analysis.py")
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)


def rows(path):
    with path.open(encoding="utf-8", newline="") as data:
        return list(csv.DictReader(data))


class LeadTriageTests(unittest.TestCase):
    def test_six_distinct_claims_with_search_trail(self):
        triage = rows(STUDY / "triage.csv")
        self.assertEqual(len(triage), 6)
        self.assertEqual({r["id"]: r["verdict"] for r in triage}, {
            "selpercatinib": "PARTIALLY KNOWN", "zanubrutinib": "PUBLISHED",
            "pbt1033": "PARTIALLY KNOWN", "ruboxistaurin": "PARTIALLY KNOWN",
            "nim811": "PARTIALLY KNOWN", "dars1": "PARTIALLY KNOWN",
        })
        log = rows(DERIVED / "search_log.csv")
        self.assertEqual({r["id"] for r in log}, {r["id"] for r in triage})
        for claim in triage:
            self.assertTrue(all(claim.values()))
            self.assertIn("doi.org/10.", claim["citations"])
            evidence = [r for r in log if r["id"] == claim["id"]]
            coverage = " ".join(" ".join(r.values()) for r in evidence).casefold()
            for resource in ("pubmed", "europe pmc", "biorxiv", "prism"):
                self.assertIn(resource, coverage)
            for record in evidence:
                self.assertTrue(all(record.values()))
                self.assertEqual(date.fromisoformat(record["date"]), date(2026, 11, 22))

    def test_censoring_ties_and_tied_ranks(self):
        for text in ("", "<50", ">100", "0", "-3", "nan", "inf", "abc"):
            self.assertIsNone(analysis.capped_ic50(text))
        for text in (">2500", "> 31600", "2500", "8000.7"):
            self.assertLess(abs(analysis.capped_ic50(text) - 2500), 1e-12)
        self.assertLess(abs(analysis.capped_ic50("100.5") - 100.5), 1e-12)
        ranks = analysis.centered_ranks([9, 3, 3, 1])
        for actual, expected in zip(ranks, [0.375, 0, 0, -0.375], strict=True):
            self.assertLess(abs(actual - expected), 1e-12)

    def test_identifiers_ret_exclusions_and_insufficient_coverage(self):
        panel = rows(DERIVED / "panel.csv")
        names = {r["cell_line"] for r in panel}
        self.assertEqual(len(names), len(panel))
        self.assertIn("5637", names)
        groups, reasons = analysis.eligibility(panel)
        self.assertEqual(reasons["TT"], reasons["LC-2/ad"])
        self.assertIn("RET", reasons["TT"])
        self.assertIn("ambiguous", reasons["Rh30"])
        ids = [member[0] for group in groups.values() for member in group]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertTrue(all(len(group) >= 2 for group in groups.values()))
        small = [dict(panel[0], cell_line="A", model_id="dup"),
                 dict(panel[0], cell_line="B", model_id="dup")]
        result = analysis.analyze(small)
        self.assertEqual(result["verdict"], "INCONCLUSIVE")
        self.assertEqual(result["n"], 0)
        self.assertNotIn("p", result)
        self.assertEqual(set(result["excluded"].values()), {"missing or reused model ID"})

    def test_cli_replay_and_normalized_hashes(self):
        saved = json.loads((DERIVED / "summary.json").read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            subprocess.run([sys.executable, str(STUDY / "analysis.py"), "--out", temp],
                           check=True, timeout=45)
            replay = json.loads((out / "summary.json").read_text(encoding="utf-8"))
            for key in ("excluded", "exclusion_counts", "lineage_counts", "sha256_lf"):
                self.assertEqual(replay[key], saved[key])
            for result in (saved, replay):
                self.assertEqual(result["source_rows"], 134)
                self.assertEqual(result["n"], 101)
                self.assertEqual(result["below_cap"], 19)
                self.assertEqual(result["informative_lineages"], 9)
                self.assertEqual(result["n"] + len(result["excluded"]), result["source_rows"])
                self.assertEqual(Counter(result["excluded"].values()), result["exclusion_counts"])
                self.assertEqual(sum(result["lineage_counts"].values()), result["n"])
                self.assertEqual(result["permutations"], 9999)
                self.assertEqual(result["at_least_observed"], 491)
                self.assertTrue(0.17 < result["r"] < 0.18)
                self.assertTrue(0.045 < result["p"] <= 0.05)
                self.assertEqual(result["verdict"], "FAIL")
            keys = {Path("derived/panel.csv").as_posix(), Path("PROTOCOL.md").as_posix()}
            self.assertEqual(set(saved["sha256_lf"]), keys)
            for key, expected in saved["sha256_lf"].items():
                source = STUDY / key
                self.assertEqual(analysis.digest(source), expected)
                crlf = out / source.name
                crlf.write_bytes(source.read_text(encoding="utf-8").replace("\n", "\r\n").encode())
                self.assertEqual(analysis.digest(crlf), expected)


if __name__ == "__main__":
    unittest.main()
