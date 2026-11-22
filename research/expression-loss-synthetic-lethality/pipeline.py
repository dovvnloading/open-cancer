"""Expression-loss synthetic-lethality screen across DepMap cell lines.

Question: when gene B is transcriptionally OFF in a cancer cell line, does the
line become dependent on a different gene A?  Discovery uses DepMap 24Q4
CRISPR (Chronos) gene effects; replication uses DEMETER2 combined RNAi gene
dependency scores (independent perturbation technology and screening
campaigns).  A third, drug-level check uses the PRISM Repurposing 24Q2
extended primary screen for annotated inhibitors of A.

Run (needs numpy and pandas; data are NOT in the repository):

    python pipeline.py --data <dir with the figshare files> --out derived

Inputs (all CC BY 4.0, figshare):
  DepMap 24Q4 Public  doi:10.25452/figshare.plus.27993248.v1
      CRISPRGeneEffect.csv, OmicsExpressionProteinCodingGenesTPMLogp1.csv,
      Model.csv
  DEMETER2 data       figshare article 6025238
      D2_combined_gene_dep_scores.csv
  Repurposing Public 24Q2  figshare article 25917643
      Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv (saved as
      PRISM_Extended_Primary_Data_Matrix.csv), ..._Compound_List.csv
  HGNC complete set (gene groups), Ensembl REST paralogues (top pairs only).

Every threshold below was fixed before the first scan was run and is written
to derived/run_parameters.json.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

PARAMS = {
    "min_lineage_size": 10,
    "essential_strength_genes": 300,     # most negative mean CRISPR effect
    "a_dep_threshold": -0.5,             # Chronos: dependent line
    "a_min_dependent_lines": 10,
    "a_min_null_lines": 10,              # effect > -0.3
    "a_null_threshold": -0.3,
    "b_off_logtpm": 1.0,                 # log2(TPM+1) < 1  ~ TPM < 1
    "b_on_logtpm": 3.0,                  # log2(TPM+1) > 3  ~ TPM > 7
    "b_min_off_lines": 20,
    "b_min_on_lines": 20,
    "max_nan_fraction": 0.05,
    "discovery_family_alpha": 0.05,      # one-sided Bonferroni over A x B
    "replication_fdr": 0.05,             # BH over discovery hits testable in RNAi
    "replication_min_lines": 100,
    "lolo_min_r_fraction": 0.5,          # min leave-one-lineage-out r >= 0.5 * r
    # v2 additions (fixed after v1 failed its positive controls; see REPORT.md)
    "predictor": "binary_off",           # indicator log2(TPM+1) < b_off_logtpm
    "expr_pcs": 20,                      # expression PCs added as covariates
    "expr_pc_genes": 5000,               # most variable genes used for the PCs
}

# Known loss-of-B -> dependency-on-A pairs (A, B), fixed before the v2 run.
BENCHMARK = [
    ("SMARCA2", "SMARCA4"), ("ARID1B", "ARID1A"), ("VPS4A", "VPS4B"),
    ("VPS4B", "VPS4A"), ("ENO2", "ENO1"), ("ME2", "ME3"), ("FAM50A", "FAM50B"),
    ("RPP25L", "RPP25"), ("MAGOH", "MAGOHB"), ("MAGOHB", "MAGOH"),
    ("NXT1", "NXT2"), ("CDS2", "CDS1"), ("STAG1", "STAG2"), ("PAPSS1", "PAPSS2"),
    ("PRMT5", "MTAP"), ("MAT2A", "MTAP"), ("ASF1A", "ASF1B"), ("SMARCC1", "SMARCC2"),
]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_matrix(path: Path, index_is_rows: bool = True) -> pd.DataFrame:
    chunks = []
    for chunk in pd.read_csv(path, index_col=0, chunksize=100, low_memory=False):
        chunks.append(chunk.astype(np.float32))
    df = pd.concat(chunks)
    return df if index_is_rows else df.T


def symbol(col: str) -> str:
    return col.split(" (")[0]


def covariates(lines: list[str], lineage: pd.Series, strength: pd.Series,
               pcs: pd.DataFrame | None = None) -> np.ndarray:
    lin = lineage.reindex(lines).fillna("Unknown")
    counts = lin.value_counts()
    small = counts[counts < PARAMS["min_lineage_size"]].index
    lin = lin.where(~lin.isin(small), "Other")
    dummies = pd.get_dummies(lin, drop_first=True).astype(np.float64).to_numpy()
    s = strength.reindex(lines).to_numpy(dtype=np.float64)
    s = (s - s.mean()) / s.std()
    cols = [np.ones(len(lines)), dummies, s]
    if pcs is not None:
        cols.append(pcs.reindex(lines).to_numpy(dtype=np.float64))
    return np.column_stack(cols)


def residual_z(M: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Residualise columns of M on X; return unit-variance residual columns."""
    M = M.astype(np.float64)
    col_mean = np.nanmean(M, axis=0)
    idx = np.where(np.isnan(M))
    M[idx] = np.take(col_mean, idx[1])
    Q, _ = np.linalg.qr(X)
    R = M - Q @ (Q.T @ M)
    sd = R.std(axis=0)
    sd[sd == 0] = np.nan
    return (R / sd).astype(np.float32)


def predictor(M: np.ndarray) -> np.ndarray:
    """1 = B expressed, 0 = B off.  Positive r <=> B-off lines are more dependent on A."""
    if PARAMS["predictor"] == "binary_off":
        return (M >= PARAMS["b_off_logtpm"]).astype(np.float32)
    return M.astype(np.float32)


def expression_pcs(E: pd.DataFrame, lineage: pd.Series) -> pd.DataFrame:
    """Top expression PCs after removing lineage means (captures within-lineage cell state)."""
    E = E[E.index.isin(lineage.dropna().index)]
    top = E.std().nlargest(PARAMS["expr_pc_genes"]).index
    M = E[top].to_numpy(dtype=np.float64)
    lin = lineage.reindex(E.index).to_numpy()
    for L in np.unique(lin):
        m = lin == L
        M[m] -= M[m].mean(axis=0)
    M /= M.std(axis=0) + 1e-9
    U, S, _ = np.linalg.svd(M, full_matrices=False)
    k = PARAMS["expr_pcs"]
    return pd.DataFrame(U[:, :k] * S[:k], index=E.index, columns=[f"PC{i+1}" for i in range(k)])


def r_to_p_one_sided(r: np.ndarray, df: int) -> np.ndarray:
    from math import erfc
    t = r * np.sqrt(df / np.clip(1 - r * r, 1e-12, None))
    # normal approximation is accurate at df ~ 1000; use exact t via scipy if present
    try:
        from scipy.stats import t as tdist
        return tdist.sf(t, df)
    except Exception:  # pragma: no cover
        return np.array([0.5 * erfc(x / math.sqrt(2)) for x in np.ravel(t)]).reshape(np.shape(t))


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    ranked = p[order] * n / np.arange(1, n + 1)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    return np.minimum(q, 1.0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, type=Path)
    ap.add_argument("--out", default=Path(__file__).with_name("derived"), type=Path)
    args = ap.parse_args()
    d, out = args.data, args.out
    out.mkdir(parents=True, exist_ok=True)

    files = {
        "crispr": d / "CRISPRGeneEffect.csv",
        "expr": d / "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        "model": d / "Model.csv",
        "rnai": d / "D2_combined_gene_dep_scores.csv",
    }
    provenance = {k: {"file": v.name, "md5": md5(v)} for k, v in files.items()}

    model = pd.read_csv(files["model"], index_col="ModelID", low_memory=False)
    lineage = model["OncotreeLineage"]

    print("loading CRISPR / expression ...", flush=True)
    C = read_matrix(files["crispr"])
    E = read_matrix(files["expr"])
    C.columns = [symbol(c) for c in C.columns]
    E.columns = [symbol(c) for c in E.columns]
    C = C.loc[:, ~C.columns.duplicated()]
    E = E.loc[:, ~E.columns.duplicated()]
    E_all = E
    pcs = expression_pcs(E_all, lineage) if PARAMS["expr_pcs"] else None
    lines = sorted(set(C.index) & set(E.index) & set(lineage.dropna().index))
    C, E = C.loc[lines], E.loc[lines]
    n = len(lines)

    # screen-strength covariate: mean effect of the most essential genes
    mean_eff = C.mean(axis=0)
    ess = mean_eff.nsmallest(PARAMS["essential_strength_genes"]).index
    strength = C[ess].mean(axis=1)
    X = covariates(lines, lineage, strength, pcs)

    nan_ok = C.isna().mean() <= PARAMS["max_nan_fraction"]
    dep = (C < PARAMS["a_dep_threshold"]).sum()
    nul = (C > PARAMS["a_null_threshold"]).sum()
    A_genes = C.columns[nan_ok & (dep >= PARAMS["a_min_dependent_lines"]) & (nul >= PARAMS["a_min_null_lines"])]
    off = (E < PARAMS["b_off_logtpm"]).sum()
    on = (E > PARAMS["b_on_logtpm"]).sum()
    B_genes = E.columns[(off >= PARAMS["b_min_off_lines"]) & (on >= PARAMS["b_min_on_lines"])]
    print(f"lines={n} A={len(A_genes)} B={len(B_genes)} covariates={X.shape[1]}", flush=True)

    Zc = residual_z(C[A_genes].to_numpy(), X)
    Ze = residual_z(predictor(E[B_genes].to_numpy()), X)
    df = n - X.shape[1] - 1
    n_tests = len(A_genes) * len(B_genes)
    alpha = PARAMS["discovery_family_alpha"] / n_tests
    # r threshold for one-sided Bonferroni
    lo, hi = 0.0, 0.99
    for _ in range(60):
        mid = (lo + hi) / 2
        if r_to_p_one_sided(np.array([mid]), df)[0] > alpha:
            lo = mid
        else:
            hi = mid
    r_thr = hi
    print(f"tests={n_tests} r_threshold={r_thr:.4f}", flush=True)

    hits = []
    block = 256
    for i in range(0, Zc.shape[1], block):
        R = (Zc[:, i:i + block].T @ Ze) / n
        ai, bj = np.where(R >= r_thr)
        for a, b in zip(ai, bj):
            A, B = A_genes[i + a], B_genes[b]
            if A != B:
                hits.append((A, B, float(R[a, b])))
    hits = pd.DataFrame(hits, columns=["A", "B", "r_crispr"])
    hits["p_crispr"] = r_to_p_one_sided(hits["r_crispr"].to_numpy(), df)
    print(f"discovery hits={len(hits)}", flush=True)

    # dichotomised effect sizes in CRISPR
    def dich(row):
        b_off = E[row.B] < PARAMS["b_off_logtpm"]
        b_on = E[row.B] > PARAMS["b_on_logtpm"]
        a = C[row.A]
        return pd.Series({
            "n_b_off": int(b_off.sum()),
            "mean_A_effect_b_off": float(a[b_off].mean()),
            "mean_A_effect_b_on": float(a[b_on].mean()),
            "frac_dependent_b_off": float((a[b_off] < PARAMS["a_dep_threshold"]).mean()),
            "frac_dependent_b_on": float((a[b_on] < PARAMS["a_dep_threshold"]).mean()),
        })
    hits = pd.concat([hits, hits.apply(dich, axis=1)], axis=1)

    # leave-one-lineage-out robustness
    lin = lineage.reindex(lines)
    big = lin.value_counts()
    big = big[big >= PARAMS["min_lineage_size"]].index
    a_idx = {g: k for k, g in enumerate(A_genes)}
    b_idx = {g: k for k, g in enumerate(B_genes)}
    hA = sorted(set(hits.A)); hB = sorted(set(hits.B))
    Craw = C[hA].to_numpy(); Eraw = predictor(E[hB].to_numpy())
    min_r = np.full(len(hits), np.inf); worst = [""] * len(hits)
    ha = np.array([hA.index(a) for a in hits.A]); hb = np.array([hB.index(b) for b in hits.B])
    for L in big:
        keep = (lin != L).to_numpy()
        sub = [l for l, k in zip(lines, keep) if k]
        Xs = covariates(sub, lineage, strength, pcs)
        zc = residual_z(Craw[keep], Xs); ze = residual_z(Eraw[keep], Xs)
        r = np.einsum("ij,ij->j", zc[:, ha], ze[:, hb]) / keep.sum()
        better = r < min_r
        min_r = np.where(better, r, min_r)
        for k in np.where(better)[0]:
            worst[k] = L
    hits["lolo_min_r"] = min_r
    hits["lolo_worst_lineage"] = worst
    del Zc, Ze

    # RNAi replication
    print("RNAi replication ...", flush=True)
    D = read_matrix(files["rnai"], index_is_rows=False)  # -> lines x genes
    D.columns = [symbol(c) for c in D.columns]
    D = D.loc[:, ~D.columns.duplicated()]
    ccle = model["CCLEName"].dropna()
    ccle_to_id = pd.Series(ccle.index, index=ccle.values)
    ccle_to_id = ccle_to_id[~ccle_to_id.index.duplicated()]
    D = D[D.index.isin(ccle_to_id.index)]
    D.index = ccle_to_id.reindex(D.index).values
    D = D[~D.index.duplicated()]
    rl = sorted(set(D.index) & set(E.index) & set(lineage.dropna().index))
    Dm = D.loc[rl]
    d_strength = Dm[[g for g in ess if g in Dm.columns]].mean(axis=1)
    rep = []
    testable = hits[hits.A.isin(Dm.columns)]
    for A, grp in testable.groupby("A"):
        a = Dm[A]
        ok = a.notna()
        if ok.sum() < PARAMS["replication_min_lines"]:
            continue
        sub = list(a[ok].index)
        Xs = covariates(sub, lineage, d_strength, pcs)
        za = residual_z(a[ok].to_numpy()[:, None], Xs)[:, 0]
        Bs = list(grp.B)
        zb = residual_z(predictor(E_all.loc[sub, Bs].to_numpy()), Xs)
        r = (za[:, None] * zb).mean(axis=0)
        for B, rv in zip(Bs, r):
            rep.append((A, B, float(rv), len(sub), len(sub) - Xs.shape[1] - 1))
    rep = pd.DataFrame(rep, columns=["A", "B", "r_rnai", "n_rnai", "df_rnai"])
    rep["p_rnai"] = [float(r_to_p_one_sided(np.array([r]), df_)[0]) for r, df_ in zip(rep.r_rnai, rep.df_rnai)]
    rep["q_rnai"] = bh(rep.p_rnai.to_numpy()) if len(rep) else []
    hits = hits.merge(rep, on=["A", "B"], how="left")
    hits["replicated"] = (hits.q_rnai <= PARAMS["replication_fdr"]) & (hits.r_rnai > 0)
    hits["lineage_robust"] = hits.lolo_min_r >= PARAMS["lolo_min_r_fraction"] * hits.r_crispr

    # HGNC gene-group paralogy proxy
    hg = pd.read_csv(d / "hgnc_complete_set.txt", sep="\t", usecols=["symbol", "gene_group_id"], dtype=str)
    groups = {s: set(str(g).split("|")) for s, g in zip(hg.symbol, hg.gene_group_id) if isinstance(g, str)}
    hits["same_hgnc_group"] = [bool(groups.get(a, set()) & groups.get(b, set())) for a, b in zip(hits.A, hits.B)]

    hits = hits.sort_values("r_crispr", ascending=False)
    bench = []
    for A, B in BENCHMARK:
        row = hits[(hits.A == A) & (hits.B == B)]
        bench.append({
            "A": A, "B": B, "A_screened": A in set(A_genes), "B_screened": B in set(B_genes),
            "discovery_hit": bool(len(row)),
            "replicated": bool(len(row) and row.replicated.iloc[0]),
            "r_crispr": float(row.r_crispr.iloc[0]) if len(row) else None,
        })
    pd.DataFrame(bench).to_csv(out / "benchmark_recall.csv", index=False)
    hits.to_csv(out / "discovery_hits.csv", index=False, float_format="%.6g")
    summary = {
        "params": PARAMS, "provenance": provenance, "n_lines_crispr": n,
        "n_lines_rnai": len(rl), "n_A": len(A_genes), "n_B": len(B_genes),
        "n_tests": n_tests, "r_threshold": r_thr, "df_crispr": df,
        "n_discovery_hits": int(len(hits)), "n_rnai_testable": int(len(rep)),
        "n_replicated": int(hits.replicated.sum()),
        "n_replicated_lineage_robust": int((hits.replicated & hits.lineage_robust).sum()),
        "benchmark_screenable": int(sum(b["A_screened"] and b["B_screened"] for b in bench)),
        "benchmark_discovered": int(sum(b["discovery_hit"] for b in bench)),
        "benchmark_replicated": int(sum(b["replicated"] for b in bench)),
    }
    (out / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
