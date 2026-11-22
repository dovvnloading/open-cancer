"""Tests of the frozen hypotheses in HYPOTHESES.md plus Sanger replication.

    python pathways.py --data <DepMap/RNAi/PRISM dir> --sanger <Sanger dir> --out derived

Needs numpy/pandas/scipy. Sanger Project Score values are used only to compute
summary statistics; none are written to disk.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import pipeline as P

SEED = 20261122
N_NULL = 10_000

HYP = {
    "H1": {"B": "GJA1", "seeds": ["PFAS", "ATIC", "CTPS1", "YRDC"],
           "set": ["PPAT", "GART", "PAICS", "ADSL", "IMPDH2", "GMPS", "ADSS2", "CAD", "DHODH",
                   "UMPS", "CTPS2", "PRPS1", "TYMS", "DHFR", "RRM1", "RRM2"]},
    "H2": {"B": "ALDH2", "seeds": ["BRIP1"],
           "set": ["FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCL",
                   "FANCM", "UBE2T", "FAAP24", "FAAP100", "SLX4", "ERCC4", "ERCC1"],
           "negative_control": ["MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7"]},
    "H3": {"B": "ACOT9", "seeds": ["SUCLG1", "SUCLA2"], "set": ["SUCLG2", "OGDH", "DLST", "DLD"]},
}
H1_DRUG_TARGETS = {"PPAT", "GART", "PAICS", "ADSL", "ATIC", "PFAS", "IMPDH1", "IMPDH2", "GMPS",
                   "CAD", "DHODH", "UMPS", "CTPS1", "TYMS", "DHFR", "RRM1", "RRM2"}


def covars(lines, lineage, strength, pcs):
    lin = lineage.reindex(lines).fillna("Unknown")
    counts = lin.value_counts()
    lin = lin.where(~lin.isin(counts[counts < P.PARAMS["min_lineage_size"]].index), "Other")
    cols = [np.ones(len(lines)), pd.get_dummies(lin, drop_first=True).astype(float).to_numpy()]
    if strength is not None:
        s = strength.reindex(lines).to_numpy(float)
        cols.append((s - s.mean()) / s.std())
    cols.append(pcs.reindex(lines).to_numpy(float))
    return np.column_stack(cols)


def partial_r(y: pd.DataFrame, b: pd.Series, X: np.ndarray) -> pd.Series:
    """Partial r between each column of y and predictor(b) given X (rows aligned)."""
    zy = P.residual_z(y.to_numpy(), X)
    zb = P.residual_z(P.predictor(b.to_numpy()[:, None]), X)[:, 0]
    return pd.Series((zy * zb[:, None]).mean(axis=0), index=y.columns)


def set_test(r: pd.Series, members: list[str], rng) -> dict:
    present = [g for g in members if g in r.index and np.isfinite(r[g])]
    if not present:
        return {"n": 0}
    obs = float(r[present].mean())
    pool = r.dropna().to_numpy()
    null = np.array([rng.choice(pool, len(present), replace=False).mean() for _ in range(N_NULL)])
    return {"n": len(present), "members": present, "mean_r": obs,
            "per_gene_r": {g: round(float(r[g]), 4) for g in present},
            "empirical_p_one_sided": float((1 + (null >= obs).sum()) / (N_NULL + 1))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--sanger", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path(__file__).with_name("derived"))
    a = ap.parse_args()
    d = a.data
    rng = np.random.default_rng(SEED)

    model = pd.read_csv(d / "Model.csv", index_col="ModelID", low_memory=False)
    lineage = model["OncotreeLineage"]
    C = P.read_matrix(d / "CRISPRGeneEffect.csv"); C.columns = [P.symbol(c) for c in C.columns]
    E = P.read_matrix(d / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"); E.columns = [P.symbol(c) for c in E.columns]
    C = C.loc[:, ~C.columns.duplicated()]; E = E.loc[:, ~E.columns.duplicated()]
    pcs = P.expression_pcs(E, lineage)
    lines = sorted(set(C.index) & set(E.index) & set(lineage.dropna().index))
    C = C.loc[lines]
    C = C.loc[:, C.isna().mean() <= P.PARAMS["max_nan_fraction"]]
    ess = C.mean().nsmallest(P.PARAMS["essential_strength_genes"]).index
    strength = C[ess].mean(axis=1)
    X = covars(lines, lineage, strength, pcs)

    D = P.read_matrix(d / "D2_combined_gene_dep_scores.csv", index_is_rows=False)
    D.columns = [P.symbol(c) for c in D.columns]; D = D.loc[:, ~D.columns.duplicated()]
    cc = model["CCLEName"].dropna(); m = pd.Series(cc.index, index=cc.values)
    m = m[~m.index.duplicated()]
    D = D[D.index.isin(m.index)]; D.index = m.reindex(D.index).values; D = D[~D.index.duplicated()]
    rl = sorted(set(D.index) & set(E.index) & set(lineage.dropna().index))
    D = D.loc[rl]
    d_strength = D[[g for g in ess if g in D.columns]].mean(axis=1)

    results = {}
    for hid, h in HYP.items():
        B = h["B"]
        r = partial_r(C, E.loc[lines, B], X)
        res = {"B": B, "n_b_off_crispr_lines": int((E.loc[lines, B] < P.PARAMS["b_off_logtpm"]).sum()),
               "seed_r": {g: round(float(r[g]), 4) for g in h["seeds"] if g in r.index},
               "crispr_set": set_test(r, [g for g in h["set"] if g not in h["seeds"]], rng)}
        # RNAi
        genes = [g for g in h["set"] + h["seeds"] if g in D.columns]
        rr = {}
        for g in genes:
            ok = D[g].notna()
            sub = list(D.index[ok])
            if len(sub) < P.PARAMS["replication_min_lines"]:
                continue
            Xs = covars(sub, lineage, d_strength, pcs)
            rr[g] = round(float(partial_r(D.loc[sub, [g]], E.loc[sub, B], Xs).iloc[0]), 4)
        set_rnai = [rr[g] for g in h["set"] if g in rr and g not in h["seeds"]]
        res["rnai_per_gene_r"] = rr
        res["rnai_set_mean_r"] = float(np.mean(set_rnai)) if set_rnai else None
        if "negative_control" in h:
            res["negative_control_set"] = set_test(r, h["negative_control"], rng)
        cs = res["crispr_set"]
        res["pass"] = bool(cs.get("n") and cs["empirical_p_one_sided"] < 0.01
                           and (res["rnai_set_mean_r"] or 0) > 0)
        results[hid] = res
        print(hid, json.dumps({k: v for k, v in res.items() if k != "rnai_per_gene_r"})[:600], flush=True)

    # H1 growth-pattern split
    gp = model["GrowthPattern"].reindex(lines).fillna("Unknown")
    split = {}
    members = HYP["H1"]["set"] + HYP["H1"]["seeds"]
    for label, mask in {"adherent": gp.str.contains("Adherent", case=False) & ~gp.str.contains("Suspension", case=False),
                        "suspension": gp.str.contains("Suspension", case=False) & ~gp.str.contains("Adherent", case=False)}.items():
        sub = list(gp.index[mask])
        Xs = covars(sub, lineage, strength, pcs)
        b = E.loc[sub, "GJA1"]
        r = partial_r(C.loc[sub, [g for g in members if g in C.columns]], b, Xs)
        split[label] = {"n_lines": len(sub), "n_gja1_off": int((b < P.PARAMS["b_off_logtpm"]).sum()),
                        "set_mean_r_including_seeds": float(r.mean()),
                        "per_gene_r": {g: round(float(v), 4) for g, v in r.items()}}
    results["H1_growth_pattern_split"] = split
    print("H1 split", {k: (v["n_lines"], v["n_gja1_off"], round(v["set_mean_r_including_seeds"], 4)) for k, v in split.items()}, flush=True)

    # H1 PRISM drug-level test
    cl = pd.read_csv(d / "PRISM_Extended_Primary_Compound_List.csv")
    M = pd.read_csv(d / "PRISM_Extended_Primary_Data_Matrix.csv", index_col=0)
    M.index = cl["IDs"].values if (M.index == cl["IDs"]).all() else M.index
    drug_rows = []
    for i, row in cl.iterrows():
        tg = set(str(row["repurposing_target"]).replace(" ", "").split(","))
        if tg & H1_DRUG_TARGETS:
            v = M.iloc[i].dropna()
            sub = sorted(set(v.index) & set(E.index) & set(lineage.dropna().index))
            if len(sub) < 200:
                continue
            Xs = covars(sub, lineage, None, pcs)
            b = E.loc[sub, "GJA1"]
            r = float(partial_r(pd.DataFrame({"lfc": v.reindex(sub)}), b, Xs).iloc[0])
            dfree = len(sub) - Xs.shape[1] - 1
            drug_rows.append({"drug": row["Drug.Name"], "screen": row["screen"], "dose": row["dose"],
                              "targets": ",".join(sorted(tg & H1_DRUG_TARGETS)), "n_lines": len(sub),
                              "n_gja1_off": int((b < P.PARAMS["b_off_logtpm"]).sum()),
                              "r": round(r, 4), "p_one_sided": float(P.r_to_p_one_sided(np.array([r]), dfree)[0])})
    dr = pd.DataFrame(drug_rows)
    dr.to_csv(a.out / "h1_prism_antimetabolites.csv", index=False)
    results["H1_prism"] = {"n_drug_conditions": len(dr), "mean_r": float(dr.r.mean()) if len(dr) else None,
                           "n_nominal_p_lt_0.05": int((dr.p_one_sided < 0.05).sum()) if len(dr) else 0}
    print("H1 PRISM", results["H1_prism"], flush=True)

    # Sanger Project Score replication (Sanger-screened columns only)
    s = a.sanger
    f = next(s.glob("*fold_change_values_*.tsv"))
    head = pd.read_csv(f, sep="\t", header=None, nrows=4, index_col=None)
    src = head.iloc[2, 3:].to_numpy(); sid = head.iloc[1, 3:].to_numpy(); qc = head.iloc[3, 3:].to_numpy()
    keep = np.where((src == "Sanger") & (qc.astype(str) == "TRUE"))[0]
    S = pd.read_csv(f, sep="\t", header=None, skiprows=5, low_memory=False)
    genes = S.iloc[:, 1].to_numpy()
    S = S.iloc[:, 3 + keep].astype(float); S.index = genes; S.columns = sid[keep]
    S = S.T.groupby(level=0).mean()  # models x genes
    ml = pd.read_csv(s / "model_list_20261013.csv", usecols=["model_id", "BROAD_ID"]).dropna()
    ml = ml[ml.BROAD_ID.str.startswith("ACH-")].drop_duplicates("model_id").set_index("model_id").BROAD_ID
    S = S[S.index.isin(ml.index)]; S.index = ml.reindex(S.index).values
    S = S[~S.index.duplicated()]
    sl = sorted(set(S.index) & set(E.index) & set(lineage.dropna().index))
    S = S.loc[sl]
    s_strength = S[[g for g in ess if g in S.columns]].mean(axis=1)
    Xs = covars(sl, lineage, s_strength, pcs)
    ann = pd.read_csv(a.out / "discovery_hits_annotated.csv")
    cand = ann[ann["class"].isin(["non_paralog_candidate", "paralog", "chrY_loss_proxy"])].copy()
    rows = []
    for A, B in zip(cand.A, cand.B):
        if A not in S.columns or B not in E.columns:
            rows.append((A, B, np.nan, 0, 0)); continue
        y = S[[A]]
        ok = y[A].notna().to_numpy()
        sub = [l for l, k in zip(sl, ok) if k]
        r = float(partial_r(S.loc[sub, [A]], E.loc[sub, B], Xs[ok]).iloc[0])
        rows.append((A, B, r, len(sub), int((E.loc[sub, B] < P.PARAMS["b_off_logtpm"]).sum())))
    sr = pd.DataFrame(rows, columns=["A", "B", "r_sanger", "n_sanger", "n_b_off_sanger"])
    dfree = sr.n_sanger - Xs.shape[1] - 1
    sr["p_sanger"] = [float(P.r_to_p_one_sided(np.array([r]), int(k))[0]) if np.isfinite(r) and k > 0 else np.nan
                      for r, k in zip(sr.r_sanger, dfree)]
    testable = sr.p_sanger.notna()
    sr["q_sanger"] = np.nan
    sr.loc[testable, "q_sanger"] = P.bh(sr.loc[testable, "p_sanger"].to_numpy())
    out = cand.merge(sr, on=["A", "B"], how="left")
    out.to_csv(a.out / "candidates_with_sanger.csv", index=False, float_format="%.6g")
    results["sanger"] = {"n_models": len(sl), "n_pairs_testable": int(testable.sum()),
                         "n_q_le_0.05": int((sr.q_sanger <= 0.05).sum())}
    for c in ["paralog", "non_paralog_candidate", "chrY_loss_proxy"]:
        sub = out[(out["class"] == c) & out.q_sanger.notna()]
        results["sanger"][c] = {"testable": int(len(sub)), "q_le_0.05": int((sub.q_sanger <= 0.05).sum())}
    print("Sanger", results["sanger"], flush=True)
    (a.out / "hypothesis_tests.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
