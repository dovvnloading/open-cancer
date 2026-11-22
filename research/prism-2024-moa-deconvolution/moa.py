"""PRISM 24Q2 compound x CRISPR knockout similarity screen (see PROTOCOL.md).

    python moa.py --data <DepMap/PRISM/RNAi dir> --sanger <Sanger dir> --out derived

Needs numpy, pandas, scipy. Input files are not in the repository (CC BY 4.0 figshare
releases listed in README.md; Sanger data used for summary statistics only).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as tdist

MIN_LINEAGE = 10
N_PCS = 20
PC_GENES = 5000
ESS_GENES = 300
MIN_ACTIVE_LINES = 20
ALPHA = 0.05


def sym(c: str) -> str:
    return c.split(" (")[0]


def read_matrix(path: Path) -> pd.DataFrame:
    return pd.concat(c.astype(np.float32) for c in pd.read_csv(path, index_col=0, chunksize=100, low_memory=False))


def expression_pcs(E: pd.DataFrame, lineage: pd.Series) -> pd.DataFrame:
    E = E[E.index.isin(lineage.dropna().index)]
    M = E[E.std().nlargest(PC_GENES).index].to_numpy(np.float64)
    lin = lineage.reindex(E.index).to_numpy()
    for L in np.unique(lin):
        M[lin == L] -= M[lin == L].mean(axis=0)
    M /= M.std(axis=0) + 1e-9
    U, S, _ = np.linalg.svd(M, full_matrices=False)
    return pd.DataFrame(U[:, :N_PCS] * S[:N_PCS], index=E.index)


def covars(lines, lineage, pcs, strength=None):
    lin = lineage.reindex(lines).fillna("Unknown")
    vc = lin.value_counts()
    lin = lin.where(~lin.isin(vc[vc < MIN_LINEAGE].index), "Other")
    cols = [np.ones(len(lines)), pd.get_dummies(lin, drop_first=True).astype(float).to_numpy(),
            pcs.reindex(lines).to_numpy(float)]
    if strength is not None:
        s = strength.reindex(lines).to_numpy(float)
        cols.append((s - s.mean()) / s.std())
    return np.column_stack(cols)


def resid_z(M: np.ndarray, X: np.ndarray) -> np.ndarray:
    M = np.array(M, dtype=np.float64)
    if np.isnan(M).any():
        mu = np.nanmean(M, axis=0)
        i = np.where(np.isnan(M))
        M[i] = np.take(mu, i[1])
    Q, _ = np.linalg.qr(X)
    R = M - Q @ (Q.T @ M)
    sd = R.std(axis=0)
    sd[sd == 0] = np.nan
    return R / sd


def r_threshold(df: int, n_tests: int) -> float:
    lo, hi = 0.0, 0.99
    for _ in range(60):
        mid = (lo + hi) / 2
        tt = mid * np.sqrt(df / (1 - mid * mid))
        lo, hi = (mid, hi) if tdist.sf(tt, df) > ALPHA / n_tests else (lo, mid)
    return hi


def p_one(r: float, df: int) -> float:
    return float(tdist.sf(r * np.sqrt(df / max(1e-12, 1 - r * r)), df))


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def targets(s) -> set[str]:
    return set() if pd.isna(s) else {x.strip() for x in str(s).split(",") if x.strip()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--sanger", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path(__file__).with_name("derived"))
    a = ap.parse_args(); d = a.data; a.out.mkdir(exist_ok=True)

    model = pd.read_csv(d / "Model.csv", index_col="ModelID", low_memory=False)
    lineage = model["OncotreeLineage"]
    C = read_matrix(d / "CRISPRGeneEffect.csv"); C.columns = [sym(c) for c in C.columns]
    C = C.loc[:, ~C.columns.duplicated()]
    E = read_matrix(d / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"); E.columns = [sym(c) for c in E.columns]
    pcs = expression_pcs(E, lineage); del E
    C = C[C.index.isin(pcs.index) & C.index.isin(lineage.dropna().index)]
    C = C.loc[:, C.isna().mean() <= 0.05]
    strength = C[C.mean().nsmallest(ESS_GENES).index].mean(axis=1)
    genes = np.array(C.columns)

    cl = pd.read_csv(d / "PRISM_Extended_Primary_Compound_List.csv")
    M = pd.read_csv(d / "PRISM_Extended_Primary_Data_Matrix.csv", index_col=0)
    # The compound list is NOT in matrix row order: align metadata to matrix rows by Broad ID.
    cl = cl.drop_duplicates("IDs").set_index("IDs").reindex(M.index)
    cl.index.name = "IDs"; cl = cl.reset_index()
    assert (cl.IDs.values == M.index.values).all()
    shared_cols = [c for c in M.columns if c in C.index]
    M = M[shared_cols]

    # group active conditions by their set of non-missing shared lines
    groups: dict[tuple, list[int]] = {}
    for i in range(len(M)):
        v = M.iloc[i]
        ok = v.notna()
        if (v[ok] < -1).sum() < MIN_ACTIVE_LINES or ok.sum() < 100:
            continue
        groups.setdefault(tuple(np.flatnonzero(ok.to_numpy())), []).append(i)
    print(f"active conditions={sum(map(len, groups.values()))} groups={len(groups)} genes={len(genes)}", flush=True)

    rows, top_tables = [], []
    for key, idx in groups.items():
        lines = [shared_cols[j] for j in key]
        X = covars(lines, lineage, pcs, strength)
        Xd = covars(lines, lineage, pcs)
        Zc = resid_z(C.loc[lines].to_numpy(), X)
        Zd = resid_z(M.iloc[idx][lines].to_numpy().T, Xd)
        R = (Zd.T @ np.nan_to_num(Zc)) / len(lines)          # conditions x genes
        df = len(lines) - X.shape[1] - 1
        thr = r_threshold(df, len(genes))
        for k, i in enumerate(idx):
            r = R[k]
            order = np.argsort(-r)
            tg = targets(cl.loc[i, "repurposing_target"])
            ranks = {g: int(np.where(genes[order] == g)[0][0]) + 1 for g in tg if g in set(genes)}
            rows.append({
                "row": i, "drug": cl.loc[i, "Drug.Name"], "screen": cl.loc[i, "screen"],
                "broad_id": cl.loc[i, "IDs"], "targets": ",".join(sorted(tg)), "moa": cl.loc[i, "MOA"],
                "n_lines": len(lines), "n_lfc_lt_-1": int((M.iloc[i][lines] < -1).sum()), "df": df, "r_threshold": thr,
                "top_gene": genes[order[0]], "top_r": float(r[order[0]]),
                "top5": ";".join(f"{genes[j]}:{r[j]:.3f}" for j in order[:5]),
                "best_target_rank": min(ranks.values()) if ranks else None,
                "best_target": min(ranks, key=ranks.get) if ranks else "",
                "best_target_r": float(r[np.where(genes == min(ranks, key=ranks.get))[0][0]]) if ranks else None,
                "n_targets_screened": len(ranks),
            })
        print(f"group n_lines={len(lines)} conditions={len(idx)} thr={thr:.3f}", flush=True)
    res = pd.DataFrame(rows)
    res["top_is_hit"] = res.top_r >= res.r_threshold
    res.to_csv(a.out / "drug_gene_similarity.csv", index=False, float_format="%.5g")

    # calibration
    ann = res[res.n_targets_screened > 0]
    n_genes = len(genes)
    chance_top10 = float((1 - (1 - 10 / n_genes) ** ann.n_targets_screened).mean())
    cal = {"n_annotated_evaluable": int(len(ann)), "n_genes": n_genes,
           "top1_rate": float((ann.best_target_rank == 1).mean()),
           "top10_rate": float((ann.best_target_rank <= 10).mean()),
           "chance_top10_rate": chance_top10}
    cal["fold_over_chance"] = cal["top10_rate"] / chance_top10
    cal["pass"] = bool(cal["top10_rate"] >= 0.10 and cal["fold_over_chance"] >= 20)
    print("calibration", cal, flush=True)

    # discovery set
    new = res[res.screen.isin(["REP.1M", "REP.300"]) & res.top_is_hit].copy()
    new["top_is_annotated_target"] = [g in set(t.split(",")) for g, t in zip(new.top_gene, new.targets)]
    disc = new[~new.top_is_annotated_target].copy()

    # replication: Sanger + RNAi for (drug, top_gene)
    f = next(a.sanger.glob("*fold_change_values_*.tsv"))
    head = pd.read_csv(f, sep="\t", header=None, nrows=4)
    src, sid, qc = head.iloc[2, 3:].to_numpy(), head.iloc[1, 3:].to_numpy(), head.iloc[3, 3:].astype(str).to_numpy()
    keep = np.where((src == "Sanger") & (qc == "TRUE"))[0]
    S = pd.read_csv(f, sep="\t", header=None, skiprows=5, low_memory=False)
    sg = S.iloc[:, 1].to_numpy(); S = S.iloc[:, 3 + keep].astype(float); S.index = sg; S.columns = sid[keep]
    S = S.T.groupby(level=0).mean(); S = S.loc[:, ~S.columns.duplicated()]
    ml = pd.read_csv(a.sanger / "model_list_20261013.csv", usecols=["model_id", "BROAD_ID"]).dropna()
    ml = ml[ml.BROAD_ID.str.startswith("ACH-")].drop_duplicates("model_id").set_index("model_id").BROAD_ID
    S = S[S.index.isin(ml.index)]; S.index = ml.reindex(S.index).values; S = S[~S.index.duplicated()]
    D = read_matrix(d / "D2_combined_gene_dep_scores.csv").T
    D.columns = [sym(c) for c in D.columns]; D = D.loc[:, ~D.columns.duplicated()]
    cc = model["CCLEName"].dropna(); mp = pd.Series(cc.index, index=cc.values); mp = mp[~mp.index.duplicated()]
    D = D[D.index.isin(mp.index)]; D.index = mp.reindex(D.index).values; D = D[~D.index.duplicated()]

    def rep(src_mat, i, g):
        if g not in src_mat.columns:
            return np.nan, 0
        v = M.iloc[i]
        lines = [l for l in src_mat.index if l in v.index and pd.notna(v[l]) and pd.notna(src_mat.at[l, g])
                 and l in pcs.index and pd.notna(lineage.get(l))]
        if len(lines) < 60:
            return np.nan, 0
        X = covars(lines, lineage, pcs)
        zd = resid_z(v[lines].to_numpy()[:, None], X)[:, 0]
        zg = resid_z(src_mat.loc[lines, g].to_numpy()[:, None], X)[:, 0]
        return float(np.mean(zd * zg)), len(lines) - X.shape[1] - 1

    for name, mat in [("sanger", S), ("rnai", D)]:
        out = [rep(mat, i, g) for i, g in zip(disc.row, disc.top_gene)]
        disc[f"r_{name}"] = [o[0] for o in out]; disc[f"df_{name}"] = [o[1] for o in out]
        p = np.array([p_one(r, k) if np.isfinite(r) and k > 0 else np.nan for r, k in out])
        q = np.full(len(p), np.nan); ok = np.isfinite(p)
        if ok.any():
            q[ok] = bh(p[ok])
        disc[f"q_{name}"] = q
    disc["replicated"] = (disc.q_sanger <= 0.05) | (disc.q_rnai <= 0.05)
    disc = disc.sort_values("top_r", ascending=False)
    disc.to_csv(a.out / "new_compound_discoveries.csv", index=False, float_format="%.5g")
    summ = {"calibration": cal,
            "n_active_conditions": int(len(res)),
            "n_new_active": int(res.screen.isin(["REP.1M", "REP.300"]).sum()),
            "n_new_with_hit": int(len(new)),
            "n_new_hit_top_is_annotated_target": int(new.top_is_annotated_target.sum()),
            "n_discovery_pairs": int(len(disc)),
            "n_replicated": int(disc.replicated.sum())}
    (a.out / "summary.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
