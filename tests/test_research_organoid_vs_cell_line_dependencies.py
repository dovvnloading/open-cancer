"""Offline checks for research/organoid-vs-cell-line-dependencies (stdlib only)."""

import csv
import json
import math
import unittest
from pathlib import Path

DER = Path(__file__).resolve().parents[1] / "research" / "organoid-vs-cell-line-dependencies" / "derived"


def load(name):
    return json.loads((DER / name).read_text(encoding="utf-8"))


def rows(name):
    with open(DER / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


class TestOrganoid(unittest.TestCase):
    def test_primary_counts_and_platform_signature(self):
        s = load("summary.json")
        self.assertEqual(s["crc"]["n_organoids"], 84)
        self.assertEqual(s["crc"]["n_cell_lines"], 29)
        hits = rows("crc_hits_with_esca_replication.csv")
        self.assertEqual(len(hits), s["n_crc_hits"])
        self.assertEqual(sum(h["replicated"] == "True" for h in hits), s["n_replicated"])
        for h in hits:
            self.assertLessEqual(float(h["q"]), 0.01)
            self.assertGreaterEqual(abs(float(h["diff"])), 0.3)
        genes = {h["gene"] for h in hits}
        self.assertIn("KRTAP4-11", genes)          # multi-copy artefact flags a platform effect
        crc = {r["gene"]: float(r["diff"]) for r in rows("crc_all_genes.csv") if r["diff"]}
        esc = {r["gene"]: float(r["diff"]) for r in rows("esca_all_genes.csv") if r["diff"]}
        common = [g for g in crc if g in esc and math.isfinite(crc[g]) and math.isfinite(esc[g])]
        x = [crc[g] for g in common]; y = [esc[g] for g in common]
        mx, my = sum(x) / len(x), sum(y) / len(y)
        r = sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(
            sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
        self.assertGreater(r, 0.7)                   # genome-wide shared axis between cancer types

    def test_within_platform_positive_control(self):
        w = load("within_platform_msi.json")
        self.assertTrue(w["WRN_positive_control_pass"])
        self.assertLess(w["WRN"]["diff_org"], 0)
        self.assertLess(w["WRN"]["diff_line"], w["WRN"]["diff_org"])

    def test_wrn_attenuation_not_specific(self):
        a = load("wrn_attenuation.json")
        self.assertEqual(a["verdict"], "general compression / not specific")
        self.assertTrue(0.3 < a["WRN_MSI"]["ratio"] < 0.6)
        self.assertTrue(0.3 < a["KRAS_mut"]["ratio"] < 0.7)
        lo_k = a["KRAS_mut"]["ratio_ci95"][0]
        self.assertGreater(a["WRN_MSI"]["ratio_ci95"][1], lo_k)   # CIs overlap
        self.assertLess(a["BRAF_V600E"]["ratio"], 0.3)

    def test_no_per_model_sanger_values(self):
        for p in DER.iterdir():
            self.assertNotIn("SIDM0", p.read_text(encoding="utf-8"), p.as_posix())


if __name__ == "__main__":
    unittest.main()
