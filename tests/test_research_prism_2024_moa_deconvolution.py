"""Offline checks for research/prism-2024-moa-deconvolution (stdlib only)."""

import csv
import json
import math
import unittest
from pathlib import Path

DER = Path(__file__).resolve().parents[1] / "research" / "prism-2024-moa-deconvolution" / "derived"


def rows(name):
    with open(DER / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def f(x):
    return float(x) if x not in ("", None) else math.nan


class TestMoa(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = json.loads((DER / "summary.json").read_text(encoding="utf-8"))
        cls.sim = rows("drug_gene_similarity.csv")
        cls.disc = rows("new_compound_discoveries.csv")

    def test_calibration_recomputed_and_gate_failed(self):
        ann = [r for r in self.sim if int(r["n_targets_screened"] or 0) > 0]
        self.assertEqual(len(ann), self.s["calibration"]["n_annotated_evaluable"])
        top10 = sum(f(r["best_target_rank"]) <= 10 for r in ann) / len(ann)
        self.assertAlmostEqual(top10, self.s["calibration"]["top10_rate"], places=9)
        self.assertLess(top10, 0.10)                       # frozen gate not met
        self.assertGreater(top10 / self.s["calibration"]["chance_top10_rate"], 20)
        self.assertFalse(self.s["calibration"]["pass"])

    def test_discovery_rules(self):
        self.assertEqual(len(self.disc), self.s["n_discovery_pairs"])
        for r in self.disc:
            self.assertIn(r["screen"], ("REP.1M", "REP.300"))
            self.assertGreaterEqual(f(r["top_r"]), f(r["r_threshold"]) - 1e-9)
            self.assertNotIn(r["top_gene"], r["targets"].split(","))
            rep = (f(r["q_sanger"]) <= 0.05) or (f(r["q_rnai"]) <= 0.05)
            self.assertEqual(rep, r["replicated"] == "True")
        self.assertEqual(sum(r["replicated"] == "True" for r in self.disc), self.s["n_replicated"])

    def test_known_target_recoveries(self):
        top = {r["drug"]: r for r in self.disc}
        for drug, gene in [("MILADEMETAN", "MDM2"), ("AZD5991", "MCL1"), ("MRTX1133", "KRAS"),
                           ("RMC-4550", "PTPN11"), ("GDC-0077", "PIK3CA")]:
            self.assertEqual(top[drug]["top_gene"], gene)
            self.assertEqual(top[drug]["replicated"], "True")
            self.assertGreater(f(top[drug]["r_sanger"]), 0.15)

    def test_no_sanger_values_or_ids(self):
        for p in DER.iterdir():
            self.assertNotIn("SIDM0", p.read_text(encoding="utf-8"), p.as_posix())


if __name__ == "__main__":
    unittest.main()
