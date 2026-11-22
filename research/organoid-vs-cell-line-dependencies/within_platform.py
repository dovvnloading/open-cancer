"""Addendum A: within-platform MSI vs MSS dependency contrast, organoids vs cell lines.

    python within_platform.py --sanger <dir> --out derived
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import ttest_ind, spearmanr
from compare import load_fc, normalise, per_model, bh

ap = argparse.ArgumentParser(); ap.add_argument("--sanger", type=Path, required=True)
ap.add_argument("--out", type=Path, default=Path(__file__).with_name("derived")); a = ap.parse_args(); s = a.sanger
ml = pd.read_csv(s / "model_list_20261013.csv", low_memory=False).set_index("model_id")
sup = pd.read_csv(s / "CRISPR_supplemental_information_20260723.csv")
core = set(sup.symbol[sup.core_fitness_pancan.astype(str).str.upper() == "TRUE"])
Xo, mo = load_fc(s / "Organoids_fitness_scores_fold_change_values_20260723.tsv")
Xl, mlm = load_fc(next(s.glob("Project_score_combined_*fold_change_values_*.tsv")))
O = per_model(normalise(Xo, core), mo, (mo.qc == "TRUE").to_numpy())
L = per_model(normalise(Xl, core), mlm, ((mlm.source == "Sanger") & (mlm.qc == "TRUE")).to_numpy())

def contrast(M, mtype):
    ids = [m for m in M.columns if m in ml.index and ml.loc[m, "model_type"] == mtype
           and ml.loc[m, "cancer_type"] == "Colorectal Carcinoma" and ml.loc[m, "msi_status"] in ("MSI", "MSS")]
    msi = [m for m in ids if ml.loc[m, "msi_status"] == "MSI"]; mss = [m for m in ids if ml.loc[m, "msi_status"] == "MSS"]
    X = M.loc[:, ids]; X = X[X.notna().mean(axis=1) >= 0.8]
    t = ttest_ind(X[msi].to_numpy(), X[mss].to_numpy(), axis=1, equal_var=False, nan_policy="omit")
    r = pd.DataFrame({"diff": X[msi].mean(axis=1) - X[mss].mean(axis=1), "t": np.asarray(t.statistic),
                      "p_two": np.asarray(t.pvalue)}, index=X.index)
    r["p_msi_more_essential"] = np.where(r["diff"] < 0, r.p_two / 2, 1 - r.p_two / 2)
    r["q"] = bh(r.p_two.fillna(1).to_numpy())
    return r, len(msi), len(mss)

ro, no_msi, no_mss = contrast(O, "Organoid")
rl, nl_msi, nl_mss = contrast(L, "Cell Line")
j = ro.join(rl, lsuffix="_org", rsuffix="_line", how="inner").dropna(subset=["diff_org", "diff_line"])
rho = spearmanr(j.diff_org, j.diff_line).statistic
line_sp = j[(j.q_line <= 0.05) & (j.diff_line <= -0.3)]
transfer = line_sp[(line_sp.diff_org <= 0) & (line_sp.p_msi_more_essential_org <= 0.05)]
org_only = j[(j.q_org <= 0.05) & (j.diff_org <= -0.3) & (j.p_two_line > 0.2)].sort_values("diff_org")
res = {"n_organoid_msi": no_msi, "n_organoid_mss": no_mss, "n_line_msi": nl_msi, "n_line_mss": nl_mss,
       "n_genes": int(len(j)), "spearman_diff_org_vs_line": float(rho),
       "WRN": {k: float(j.loc["WRN", k]) for k in ["diff_org", "p_msi_more_essential_org", "diff_line", "p_msi_more_essential_line"]},
       "WRN_positive_control_pass": bool(j.loc["WRN", "p_msi_more_essential_org"] <= 0.01 and j.loc["WRN", "p_msi_more_essential_line"] <= 0.01),
       "n_line_msi_specific": int(len(line_sp)), "n_transfer_to_organoids": int(len(transfer)),
       "transfer_fraction": float(len(transfer) / len(line_sp)) if len(line_sp) else None,
       "line_msi_specific_genes": sorted(line_sp.index), "transferred_genes": sorted(transfer.index),
       "n_organoid_only_msi": int(len(org_only)), "organoid_only_top": org_only.head(25).index.tolist()}
j.index.name = "gene"
j.to_csv(a.out / "crc_msi_contrast_by_platform.csv", float_format="%.5g")
(a.out / "within_platform_msi.json").write_text(json.dumps(res, indent=2)); print(json.dumps(res, indent=1))
print(org_only.head(25)[["diff_org", "q_org", "diff_line", "p_two_line"]].round(3).to_string())
