"""Paired cell-line bootstrap: CRISPR H1 set mean r vs active-PRISM-drug mean r on the
same adherent lines (Addendum B in HYPOTHESES.md), plus genome-wide GJA1 partial r.

    python bootstrap_h1.py --data <dir> --out derived
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
import pipeline as P
from pathways import HYP, covars, partial_r
from followup import load

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path(__file__).with_name("derived")); a = ap.parse_args()
    model, lineage, C, E, pcs, lines, ess, strength = load(a.data)
    X = covars(lines, lineage, strength, pcs)
    rg = partial_r(C, E.loc[lines, "GJA1"], X).sort_values(ascending=False)
    rg.rename("partial_r_GJA1_expressed").to_csv(a.out / "gja1_genomewide_partial_r.csv", float_format="%.5g")
    act = pd.read_csv(a.out / "h1_prism_active_same_lines.csv")
    cl = pd.read_csv(a.data / "PRISM_Extended_Primary_Compound_List.csv")
    M = pd.read_csv(a.data / "PRISM_Extended_Primary_Data_Matrix.csv", index_col=0)
    gp = model["GrowthPattern"].fillna("")
    adherent = set(gp.index[gp.str.contains("Adherent", case=False) & ~gp.str.contains("Suspension", case=False)])
    same = sorted(set(lines) & set(M.columns) & adherent)
    members = [g for g in HYP["H1"]["set"] if g not in HYP["H1"]["seeds"] and g in C.columns]
    rows = [i for i, r in cl.iterrows() if (r["Drug.Name"], r["dose"]) in set(zip(act.drug, act.dose))]
    D = M.iloc[rows].T.reindex(same)
    rng = np.random.default_rng(20261124)
    diffs, cr, dr = [], [], []
    for b in range(1000):
        s = list(rng.choice(same, len(same), replace=True))
        Xs = covars(s, lineage, strength, pcs)
        c = partial_r(C.loc[s, members].reset_index(drop=True), E.loc[s, "GJA1"].reset_index(drop=True), Xs).mean()
        ds = []
        for j in range(D.shape[1]):
            v = D.iloc[:, j].loc[s].reset_index(drop=True)
            ok = v.notna().to_numpy()
            Xd = covars([x for x, k in zip(s, ok) if k], lineage, None, pcs)
            ds.append(partial_r(pd.DataFrame({"v": v[ok].to_numpy()}), E.loc[s, "GJA1"].reset_index(drop=True)[ok].reset_index(drop=True), Xd).iloc[0])
        cr.append(c); dr.append(float(np.mean(ds))); diffs.append(c - float(np.mean(ds)))
    q = lambda x: [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]
    out = {"n_lines": len(same), "n_drug_conditions": D.shape[1], "B": 1000,
           "crispr_set_mean_ci95": q(cr), "drug_mean_ci95": q(dr), "difference_ci95": q(diffs),
           "ratio_drug_over_crispr_median": float(np.median(np.array(dr) / np.array(cr))),
           "frac_boot_difference_le_0": float(np.mean(np.array(diffs) <= 0))}
    (a.out / "h1_paired_bootstrap.json").write_text(json.dumps(out, indent=2)); print(out)

if __name__ == "__main__":
    main()
