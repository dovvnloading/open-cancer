"""Addendum-A checks from HYPOTHESES.md (H1 specificity, same-line contrast,
PRISM activity sensitivity, Sanger set-level replication of H1/H2).

    python followup.py --data <dir> --sanger <dir> --out derived
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import pipeline as P
from pathways import HYP, H1_DRUG_TARGETS, N_NULL, covars, partial_r, set_test

SEED = 20261123
CONTROLS = {
    "ribosome": None,  # filled from symbols RPL*/RPS* (cytosolic, excluding MRP*)
    "proteasome": [f"PSMA{i}" for i in range(1, 8)] + [f"PSMB{i}" for i in range(1, 8)],
    "spliceosome_core": ["SF3B1", "SF3B2", "SF3B3", "SF3B4", "SF3B5", "SNRPD1", "SNRPD2", "SNRPD3", "PRPF8"],
    "rna_pol_II": [f"POLR2{c}" for c in "ABCDEFGHIJKL"],
}


def load(d: Path):
    model = pd.read_csv(d / "Model.csv", index_col="ModelID", low_memory=False)
    lineage = model["OncotreeLineage"]
    C = P.read_matrix(d / "CRISPRGeneEffect.csv"); C.columns = [P.symbol(c) for c in C.columns]
    E = P.read_matrix(d / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"); E.columns = [P.symbol(c) for c in E.columns]
    C = C.loc[:, ~C.columns.duplicated()]; E = E.loc[:, ~E.columns.duplicated()]
    pcs = P.expression_pcs(E, lineage)
    lines = sorted(set(C.index) & set(E.index) & set(lineage.dropna().index))
    C = C.loc[lines]; C = C.loc[:, C.isna().mean() <= P.PARAMS["max_nan_fraction"]]
    ess = C.mean().nsmallest(P.PARAMS["essential_strength_genes"]).index
    strength = C[ess].mean(axis=1)
    return model, lineage, C, E, pcs, lines, ess, strength


def load_sanger(s: Path, E, lineage, ess):
    f = next(s.glob("*fold_change_values_*.tsv"))
    head = pd.read_csv(f, sep="\t", header=None, nrows=4)
    src = head.iloc[2, 3:].to_numpy(); sid = head.iloc[1, 3:].to_numpy(); qc = head.iloc[3, 3:].to_numpy()
    keep = np.where((src == "Sanger") & (qc.astype(str) == "TRUE"))[0]
    S = pd.read_csv(f, sep="\t", header=None, skiprows=5, low_memory=False)
    genes = S.iloc[:, 1].to_numpy()
    S = S.iloc[:, 3 + keep].astype(float); S.index = genes; S.columns = sid[keep]
    S = S.T.groupby(level=0).mean()
    S = S.loc[:, ~S.columns.duplicated()]
    ml = pd.read_csv(s / "model_list_20261013.csv", usecols=["model_id", "BROAD_ID"]).dropna()
    ml = ml[ml.BROAD_ID.str.startswith("ACH-")].drop_duplicates("model_id").set_index("model_id").BROAD_ID
    S = S[S.index.isin(ml.index)]; S.index = ml.reindex(S.index).values; S = S[~S.index.duplicated()]
    sl = sorted(set(S.index) & set(E.index) & set(lineage.dropna().index))
    S = S.loc[sl]
    S = S.loc[:, S.isna().mean() <= 0.05]
    return S, S[[g for g in ess if g in S.columns]].mean(axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--sanger", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path(__file__).with_name("derived"))
    a = ap.parse_args()
    rng = np.random.default_rng(SEED)
    model, lineage, C, E, pcs, lines, ess, strength = load(a.data)
    X = covars(lines, lineage, strength, pcs)
    res = {}

    # 1. specificity
    r_g = partial_r(C, E.loc[lines, "GJA1"], X)
    h1 = set_test(r_g, [g for g in HYP["H1"]["set"] if g not in HYP["H1"]["seeds"]], rng)
    CONTROLS["ribosome"] = sorted(g for g in C.columns if (g.startswith("RPL") or g.startswith("RPS"))
                                  and g[3:4].isdigit() and not g.endswith("L1") and "Y" not in g)
    spec = {"H1_set_mean_r": h1["mean_r"]}
    for k, genes in CONTROLS.items():
        t = set_test(r_g, genes, rng)
        t["pass"] = bool(t["mean_r"] < 0.5 * h1["mean_r"] and t["empirical_p_one_sided"] > 0.01)
        t.pop("per_gene_r", None)
        spec[k] = t
    res["H1_specificity"] = spec
    print("specificity", {k: (round(v["mean_r"], 4), v["empirical_p_one_sided"], v["pass"]) for k, v in spec.items() if isinstance(v, dict)}, flush=True)

    # 2/3. PRISM: same-line contrast and active-drug sensitivity
    cl = pd.read_csv(a.data / "PRISM_Extended_Primary_Compound_List.csv")
    M = pd.read_csv(a.data / "PRISM_Extended_Primary_Data_Matrix.csv", index_col=0)
    # The compound list is NOT in matrix row order: align metadata to matrix rows by Broad ID
    # (bug found 2026-11-22 by the MOA calibration in research/prism-2024-moa-deconvolution).
    cl = cl.drop_duplicates("IDs").set_index("IDs").reindex(M.index)
    cl.index.name = "IDs"; cl = cl.reset_index()
    gp = model["GrowthPattern"].fillna("")
    adherent = set(gp.index[gp.str.contains("Adherent", case=False) & ~gp.str.contains("Suspension", case=False)])
    prism_lines = set(M.columns)
    same = sorted(set(lines) & prism_lines & adherent)
    Xs = covars(same, lineage, strength, pcs)
    members = [g for g in HYP["H1"]["set"] if g not in HYP["H1"]["seeds"] and g in C.columns]
    r_same = partial_r(C.loc[same, members], E.loc[same, "GJA1"], Xs)
    res["H1_crispr_same_adherent_prism_lines"] = {"n_lines": len(same),
        "n_gja1_off": int((E.loc[same, "GJA1"] < P.PARAMS["b_off_logtpm"]).sum()),
        "set_mean_r": float(r_same.mean()), "n_positive": int((r_same > 0).sum()), "n_genes": len(members)}
    drows = []
    for i, row in cl.iterrows():
        tg = set(str(row["repurposing_target"]).replace(" ", "").split(",")) & H1_DRUG_TARGETS
        if not tg:
            continue
        v = M.iloc[i].reindex(same).dropna()
        if len(v) < 150 or (v < -1).mean() < 0.10:
            continue
        sub = list(v.index)
        Xd = covars(sub, lineage, None, pcs)
        r = float(partial_r(pd.DataFrame({"lfc": v}), E.loc[sub, "GJA1"], Xd).iloc[0])
        drows.append({"drug": row["Drug.Name"], "dose": row["dose"], "targets": ",".join(sorted(tg)),
                      "n_lines": len(sub), "frac_lfc_lt_-1": round(float((v < -1).mean()), 3), "r": round(r, 4)})
    dr = pd.DataFrame(drows)
    dr.to_csv(a.out / "h1_prism_active_same_lines.csv", index=False)
    boot = [dr.r.sample(len(dr), replace=True, random_state=int(s)).mean() for s in rng.integers(0, 2**31, 5000)] if len(dr) else []
    res["H1_prism_active_same_lines"] = {"n_conditions": len(dr), "mean_r": float(dr.r.mean()) if len(dr) else None,
                                         "bootstrap95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))] if len(dr) else None}
    print("same-line CRISPR", res["H1_crispr_same_adherent_prism_lines"], flush=True)
    print("PRISM active", res["H1_prism_active_same_lines"], flush=True)

    # 4. Sanger set-level replication
    S, s_strength = load_sanger(a.sanger, E, lineage, ess)
    sl = list(S.index)
    Xsg = covars(sl, lineage, s_strength, pcs)
    sang = {"n_models": len(sl)}
    for hid in ["H1", "H2"]:
        B = HYP[hid]["B"]
        r = partial_r(S, E.loc[sl, B], Xsg)
        t = set_test(r, [g for g in HYP[hid]["set"] if g not in HYP[hid]["seeds"]], rng)
        t["seed_r"] = {g: round(float(r[g]), 4) for g in HYP[hid]["seeds"] if g in r.index}
        t["n_b_off"] = int((E.loc[sl, B] < P.PARAMS["b_off_logtpm"]).sum())
        if hid == "H2":
            nc = set_test(r, HYP["H2"]["negative_control"], rng); nc.pop("per_gene_r", None)
            t["negative_control"] = nc
        sang[hid] = t
    res["sanger_set_replication"] = sang
    print("Sanger sets", json.dumps(sang)[:1500], flush=True)
    (a.out / "followup_tests.json").write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
