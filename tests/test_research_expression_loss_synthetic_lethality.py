"""Offline checks for research/expression-loss-synthetic-lethality (stdlib only)."""

import csv
import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "research" / "expression-loss-synthetic-lethality"
DER = ROOT / "derived"


def rows(name):
    with open(DER / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def load_json(name):
    return json.loads((DER / name).read_text(encoding="utf-8"))


def fnum(x):
    return float(x) if x not in ("", None) else math.nan


class TestExtractRecomputation(unittest.TestCase):
    """Recompute dichotomised penetrance from the committed CC BY 4.0 DepMap extract."""

    @classmethod
    def setUpClass(cls):
        cls.x = rows("depmap24q4_extract_selected_pairs.csv")
        cls.hits = {(r["A"], r["B"]): r for r in rows("discovery_hits.csv")}

    def frac(self, a, b, off):
        sel = [r for r in self.x if (fnum(r[b + "_log2tpm1"]) < 1.0 if off else fnum(r[b + "_log2tpm1"]) > 3.0)]
        dep = [fnum(r[a + "_gene_effect"]) < -0.5 for r in sel if not math.isnan(fnum(r[a + "_gene_effect"]))]
        return len(sel), sum(dep) / len(dep)

    def test_line_count(self):
        self.assertEqual(len(self.x), 1103)

    def test_penetrance_matches_pipeline(self):
        for a, b in [("FAM50A", "FAM50B"), ("ITPA", "PUDP"), ("SRD5A3", "HTATIP2"), ("PFAS", "GJA1")]:
            n_off, f_off = self.frac(a, b, True)
            _, f_on = self.frac(a, b, False)
            h = self.hits[(a, b)]
            self.assertEqual(n_off, int(h["n_b_off"]))
            self.assertAlmostEqual(f_off, fnum(h["frac_dependent_b_off"]), places=3)
            self.assertAlmostEqual(f_on, fnum(h["frac_dependent_b_on"]), places=3)
            self.assertGreaterEqual(f_off - f_on, 0.20)  # annotate.MIN_PENETRANCE_GAP

    def test_within_lineage_association_positive(self):
        # Lineage-centred correlation between "B expressed" and A effect (coarser than the
        # pipeline's full covariate model) must agree in sign and be substantial for FAM50A/B.
        groups = {}
        for r in self.x:
            e = fnum(r["FAM50A_gene_effect"])
            if math.isnan(e):
                continue
            groups.setdefault(r["OncotreeLineage"], []).append((1.0 if fnum(r["FAM50B_log2tpm1"]) >= 1 else 0.0, e))
        xs, ys = [], []
        for vals in groups.values():
            mx = sum(v[0] for v in vals) / len(vals); my = sum(v[1] for v in vals) / len(vals)
            xs += [v[0] - mx for v in vals]; ys += [v[1] - my for v in vals]
        sxy = sum(a * b for a, b in zip(xs, ys)); sxx = sum(a * a for a in xs); syy = sum(b * b for b in ys)
        self.assertGreater(sxy / math.sqrt(sxx * syy), 0.6)


class TestScreenSummary(unittest.TestCase):
    def test_v2_counts_and_threshold(self):
        s = load_json("run_summary.json")
        self.assertEqual(s["n_lines_crispr"], 1103)
        self.assertEqual(s["n_tests"], s["n_A"] * s["n_B"])
        self.assertEqual(s["n_discovery_hits"], len(rows("discovery_hits.csv")))
        self.assertTrue(0.17 < s["r_threshold"] < 0.19)
        self.assertEqual(s["params"]["predictor"], "binary_off")
        self.assertEqual(s["params"]["expr_pcs"], 20)
        for r in rows("discovery_hits.csv"):
            self.assertGreaterEqual(fnum(r["r_crispr"]), s["r_threshold"] - 1e-9)
            self.assertNotEqual(r["A"], r["B"])

    def test_v1_retained_as_failure(self):
        v1 = load_json("v1_run_summary.json")
        self.assertGreater(v1["n_discovery_hits"], 50 * load_json("run_summary.json")["n_discovery_hits"])

    def test_benchmark_recall(self):
        b = rows("benchmark_recall.csv")
        self.assertEqual(len(b), 18)
        screenable = [r for r in b if r["A_screened"] == "True" and r["B_screened"] == "True"]
        self.assertEqual(len(screenable), 4)
        self.assertTrue(all(r["discovery_hit"] == "True" for r in screenable))


class TestHypotheses(unittest.TestCase):
    def test_h1_h2_h3_verdicts(self):
        h = load_json("hypothesis_tests.json")
        self.assertTrue(h["H1"]["pass"])
        self.assertLess(h["H1"]["crispr_set"]["empirical_p_one_sided"], 0.01)
        self.assertFalse(h["H2"]["pass"])  # RNAi did not replicate; recorded, not hidden
        self.assertLess(h["H2"]["crispr_set"]["empirical_p_one_sided"], 0.01)
        self.assertFalse(h["H3"]["pass"])
        self.assertGreater(h["H3"]["crispr_set"]["empirical_p_one_sided"], 0.5)
        sp = h["H1_growth_pattern_split"]
        self.assertGreater(sp["adherent"]["set_mean_r_including_seeds"], 3 * sp["suspension"]["set_mean_r_including_seeds"])
        self.assertLess(abs(h["H1_prism"]["mean_r"]), 0.02)

    def test_followup_controls_and_sanger(self):
        f = load_json("followup_tests.json")
        for k in ["ribosome", "proteasome", "spliceosome_core", "rna_pol_II"]:
            self.assertTrue(f["H1_specificity"][k]["pass"], k)
        s = f["sanger_set_replication"]
        self.assertLess(s["H1"]["empirical_p_one_sided"], 0.01)
        self.assertLess(s["H2"]["empirical_p_one_sided"], 0.01)
        self.assertLess(s["H2"]["negative_control"]["mean_r"], 0.02)

    def test_paired_bootstrap_and_enrichment(self):
        b = load_json("h1_paired_bootstrap.json")
        self.assertGreater(b["difference_ci95"][0], 0)
        self.assertLess(b["drug_mean_ci95"][0], 0)
        e = load_json("h1_enzyme_enrichment.json")
        self.assertLess(e["mannwhitney_p_one_sided"], 0.01)
        self.assertGreater(e["top50_gja1_protected_common_depleted_enzyme_fraction"],
                           3 * e["enzyme_fraction_all_common_depleted"])


class TestQueue(unittest.TestCase):
    def test_strict_queue_rule(self):
        q = rows("strict_candidate_queue.csv")
        self.assertEqual(len(q), 51)
        self.assertEqual(sum(1 for r in q if not r["proxy_or_known_explanation"]), 18)
        for r in q:
            rnai_ok = fnum(r["q_rnai"]) <= 0.05
            sanger_ok = fnum(r["q_sanger"]) <= 0.05 and fnum(r["n_b_off_sanger"]) >= 5
            self.assertTrue(rnai_ok or sanger_ok, (r["A"], r["B"]))
            self.assertGreaterEqual(fnum(r["lolo_min_r"]), 0.5 * fnum(r["r_crispr"]) - 1e-6)

    def test_no_sanger_per_line_values(self):
        header = (DER / "candidates_with_sanger.csv").read_text(encoding="utf-8").splitlines()[0]
        self.assertNotIn("SIDM", header)
        for p in DER.iterdir():
            self.assertNotIn("SIDM0", p.read_text(encoding="utf-8", errors="ignore")[:20000], p.as_posix())


if __name__ == "__main__":
    unittest.main()
