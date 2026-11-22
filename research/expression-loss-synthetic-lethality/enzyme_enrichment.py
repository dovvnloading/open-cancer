"""Addendum C: are GJA1-protected dependencies enriched for metabolic enzymes?

    python enzyme_enrichment.py --data <dir> --out derived
"""
import argparse, json
from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu
import pipeline as P

ap = argparse.ArgumentParser(); ap.add_argument("--data", type=Path, required=True)
ap.add_argument("--out", type=Path, default=Path(__file__).with_name("derived")); a = ap.parse_args()
C = P.read_matrix(a.data / "CRISPRGeneEffect.csv"); C.columns = [P.symbol(c) for c in C.columns]
r = pd.read_csv(a.out / "gja1_genomewide_partial_r.csv", index_col=0).iloc[:, 0]
mean = C.mean().reindex(r.index)
hg = pd.read_csv(a.data / "hgnc_complete_set.txt", sep="\t", dtype=str, usecols=["symbol", "enzyme_id"])
ec = set(hg.symbol[hg.enzyme_id.notna()])
ess = r[mean < -0.5]
enz, non = ess[ess.index.isin(ec)], ess[~ess.index.isin(ec)]
u = mannwhitneyu(enz, non, alternative="greater")
res = {"n_common_depleted": int(len(ess)), "n_enzymes": int(len(enz)), "n_non_enzymes": int(len(non)),
       "median_r_enzymes": float(enz.median()), "median_r_non_enzymes": float(non.median()),
       "mannwhitney_p_one_sided": float(u.pvalue),
       "top50_gja1_protected_common_depleted_enzyme_fraction": float(ess.sort_values(ascending=False).head(50).index.isin(ec).mean()),
       "enzyme_fraction_all_common_depleted": float(len(enz) / len(ess))}
(a.out / "h1_enzyme_enrichment.json").write_text(json.dumps(res, indent=2)); print(res)
