"""Addendum B: WRN-in-MSI attenuation in organoids vs comparator oncogene addictions.

    python wrn_attenuation.py --sanger <dir> --out derived
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
from compare import load_fc, normalise, per_model

ap = argparse.ArgumentParser(); ap.add_argument("--sanger", type=Path, required=True)
ap.add_argument("--out", type=Path, default=Path(__file__).with_name("derived")); a = ap.parse_args(); s = a.sanger
ml = pd.read_csv(s / "model_list_20261013.csv", low_memory=False).set_index("model_id")
sup = pd.read_csv(s / "CRISPR_supplemental_information_20260723.csv")
core = set(sup.symbol[sup.core_fitness_pancan.astype(str).str.upper() == "TRUE"])
mut = pd.read_csv(s / "mutations_summary_latest.csv.gz", low_memory=False)
drv = mut[mut.cancer_driver.astype(str) == "t"]
def mutant(gene, prot=None):
    x = drv[drv.gene_symbol == gene]
    if prot: x = x[x.protein_mutation.astype(str).str.startswith(prot)]
    return set(x.model_id)
Xo, mo = load_fc(s / "Organoids_fitness_scores_fold_change_values_20260723.tsv")
Xl, mlm = load_fc(next(s.glob("Project_score_combined_*fold_change_values_*.tsv")))
P = {"organoid": per_model(normalise(Xo, core), mo, (mo.qc == "TRUE").to_numpy()),
     "cell_line": per_model(normalise(Xl, core), mlm, ((mlm.source == "Sanger") & (mlm.qc == "TRUE")).to_numpy())}
mtype = {"organoid": "Organoid", "cell_line": "Cell Line"}
screened = {k: set(v.columns) for k, v in P.items()}
contrasts = {"WRN_MSI": ("WRN", lambda m: ml.loc[m, "msi_status"] == "MSI", lambda m: ml.loc[m, "msi_status"] == "MSS"),
             "KRAS_mut": ("KRAS", lambda m: m in mutant("KRAS"), lambda m: m not in mutant("KRAS")),
             "BRAF_V600E": ("BRAF", lambda m: m in mutant("BRAF", "p.V600E"), lambda m: m not in mutant("BRAF")),
             "CTNNB1_mut": ("CTNNB1", lambda m: m in mutant("CTNNB1"), lambda m: m not in mutant("CTNNB1"))}
rng = np.random.default_rng(20261122); res = {}
for name, (gene, is_pos, is_neg) in contrasts.items():
    vals = {}
    for plat, M in P.items():
        ids = [m for m in M.columns if m in ml.index and ml.loc[m, "model_type"] == mtype[plat]
               and ml.loc[m, "cancer_type"] == "Colorectal Carcinoma" and pd.notna(M.at[gene, m]) if gene in M.index]
        if name == "WRN_MSI":
            ids = [m for m in ids if ml.loc[m, "msi_status"] in ("MSI", "MSS")]
        pos = np.array([M.at[gene, m] for m in ids if is_pos(m)]); neg = np.array([M.at[gene, m] for m in ids if is_neg(m)])
        vals[plat] = (pos, neg)
    (po, no), (pl, nl) = vals["organoid"], vals["cell_line"]
    entry = {"gene": gene, "n_pos_organoid": len(po), "n_neg_organoid": len(no), "n_pos_line": len(pl), "n_neg_line": len(nl)}
    if min(len(po), len(no), len(pl), len(nl)) >= 3:
        d_o, d_l = po.mean() - no.mean(), pl.mean() - nl.mean()
        boots = []
        for _ in range(2000):
            bo = rng.choice(po, len(po)).mean() - rng.choice(no, len(no)).mean()
            bl = rng.choice(pl, len(pl)).mean() - rng.choice(nl, len(nl)).mean()
            if bl < 0: boots.append(bo / bl)
        entry.update({"diff_organoid": float(d_o), "diff_line": float(d_l), "ratio": float(d_o / d_l),
                      "ratio_ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
                      "frac_pos_effect_lt_-0.5_organoid": float((po < -0.5).mean()),
                      "frac_pos_effect_lt_-0.5_line": float((pl < -0.5).mean())})
    res[name] = entry
eligible = [k for k in res if k != "WRN_MSI" and "ratio_ci95" in res[k]]
w = res["WRN_MSI"]["ratio_ci95"]
res["verdict"] = ("WRN-specific attenuation" if eligible and all(w[1] < res[k]["ratio_ci95"][0] for k in eligible)
                  else "general compression / not specific")
res["eligible_comparators"] = eligible
(a.out / "wrn_attenuation.json").write_text(json.dumps(res, indent=2)); print(json.dumps(res, indent=1))
