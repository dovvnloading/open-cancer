"""Organoid vs 2D cell-line CRISPR dependency comparison (see PROTOCOL.md).

    python compare.py --sanger <dir with Sanger Cell Model Passports files> --out derived

Needs numpy, pandas, scipy. Sanger data are used under their research-use licence; only
gene-level summary statistics are written.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

CLASSES = {
    "niche_signalling": ["LGR4", "LGR5", "RNF43", "ZNRF3", "CTNNB1", "TCF7L2", "EGFR", "ERBB2", "ERBB3", "KRAS",
                         "MAPK1", "PTPN11", "SOS1", "GRB2"],
    "adhesion_integrin": ["ITGAV", "ITGB1", "ITGB5", "PTK2", "TLN1", "FERMT2", "ILK", "PXN", "CRK", "CRKL", "RAC1"],
    "oxygen_metabolism": ["HIF1A", "VHL", "EGLN1", "ARNT", "PDK1", "SLC2A1"],
    "interferon_innate": ["STAT1", "IRF1", "JAK1", "IFNGR1", "MAVS"],
}


def load_fc(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    head = pd.read_csv(path, sep="\t", header=None, nrows=4)
    meta = pd.DataFrame({"model_id": head.iloc[1, 3:].to_numpy(), "source": head.iloc[2, 3:].to_numpy(),
                         "qc": head.iloc[3, 3:].astype(str).to_numpy()})
    X = pd.read_csv(path, sep="\t", header=None, skiprows=5, low_memory=False)
    genes = X.iloc[:, 1].to_numpy()
    X = X.iloc[:, 3:].astype(float)
    X.index = genes
    X.columns = range(X.shape[1])
    return X, meta


def normalise(X: pd.DataFrame, core: set[str]) -> pd.DataFrame:
    med_all = X.median()
    med_core = X[X.index.isin(core)].median()
    return (X - med_all) / (med_all - med_core)


def per_model(X: pd.DataFrame, meta: pd.DataFrame, keep: np.ndarray) -> pd.DataFrame:
    X = X.loc[:, keep]
    X.columns = meta.model_id.to_numpy()[keep]
    X = X.T.groupby(level=0).mean().T
    return X.loc[~X.index.duplicated()]


def compare(org: pd.DataFrame, lin: pd.DataFrame) -> pd.DataFrame:
    genes = org.index.intersection(lin.index)
    o, l = org.loc[genes], lin.loc[genes]
    ok = (o.notna().mean(axis=1) >= 0.8) & (l.notna().mean(axis=1) >= 0.8)
    o, l = o[ok], l[ok]
    t = ttest_ind(o.to_numpy(), l.to_numpy(), axis=1, equal_var=False, nan_policy="omit")
    res = pd.DataFrame({"mean_organoid": o.mean(axis=1), "mean_cell_line": l.mean(axis=1),
                        "diff": o.mean(axis=1) - l.mean(axis=1), "t": np.asarray(t.statistic),
                        "p": np.asarray(t.pvalue)}, index=o.index)
    return res


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def lomo_stable(org, lin, gene, sign) -> bool:
    """Leave-one-model-out in each group: the sign of the difference never flips."""
    o, l = org.loc[gene].dropna(), lin.loc[gene].dropna()
    for grp, other, s in [(o, l, 1), (l, o, -1)]:
        for m in grp.index:
            d = (grp.drop(m).mean() - other.mean()) * s
            if np.sign(d) != sign:
                return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sanger", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path(__file__).with_name("derived"))
    a = ap.parse_args(); s = a.sanger; a.out.mkdir(exist_ok=True)
    ml = pd.read_csv(s / "model_list_20261013.csv", low_memory=False).set_index("model_id")
    sup = pd.read_csv(s / "CRISPR_supplemental_information_20260723.csv")
    core = set(sup.symbol[sup.core_fitness_pancan.astype(str).str.upper() == "TRUE"])

    Xo, mo = load_fc(s / "Organoids_fitness_scores_fold_change_values_20260723.tsv")
    Xl, mlm = load_fc(next(s.glob("Project_score_combined_*fold_change_values_*.tsv")))
    Xo, Xl = normalise(Xo, core), normalise(Xl, core)
    O = per_model(Xo, mo, (mo.qc == "TRUE").to_numpy())
    L = per_model(Xl, mlm, ((mlm.source == "Sanger") & (mlm.qc == "TRUE")).to_numpy())
    O = O.loc[:, [m for m in O.columns if ml.loc[m, "model_type"] == "Organoid"]] if len(O.columns) else O
    L = L.loc[:, [m for m in L.columns if m in ml.index and ml.loc[m, "model_type"] == "Cell Line"]]

    def group(M, ctype):
        return M.loc[:, [m for m in M.columns if ml.loc[m, "cancer_type"] == ctype]]

    out, summ = {}, {"n_core_genes": len(core)}
    for label, ctype in [("crc", "Colorectal Carcinoma"), ("esca", "Esophageal Carcinoma")]:
        o, l = group(O, ctype), group(L, ctype)
        r = compare(o, l)
        r["q"] = bh(r.p.to_numpy())
        out[label] = (r, o, l)
        msi = lambda ids: ml.loc[ids, "msi_status"].fillna("NA").value_counts().to_dict()
        summ[label] = {"n_organoids": o.shape[1], "n_cell_lines": l.shape[1], "n_genes_tested": len(r),
                       "msi_organoids": msi(list(o.columns)), "msi_cell_lines": msi(list(l.columns))}
        print(label, summ[label], flush=True)

    crc, o_crc, l_crc = out["crc"]
    esc, _, _ = out["esca"]
    hits = crc[(crc.q <= 0.01) & (crc["diff"].abs() >= 0.3)].copy()
    hits = hits.join(esc[["diff", "t", "p"]].rename(columns={"diff": "diff_esca", "t": "t_esca", "p": "p_esca_two_sided"}))
    same = np.sign(hits.diff_esca) == np.sign(hits["diff"])
    p1 = np.where(same, hits.p_esca_two_sided / 2, 1 - hits.p_esca_two_sided / 2)
    hits["p_esca_one_sided"] = p1
    ok = np.isfinite(p1)
    hits["q_esca"] = np.nan
    hits.loc[ok, "q_esca"] = bh(p1[ok])
    hits["replicated"] = (hits.q_esca <= 0.05) & same
    hits["direction"] = np.where(hits["diff"] < 0, "more_essential_in_organoids", "more_essential_in_cell_lines")
    hits["lomo_stable_crc"] = [lomo_stable(o_crc, l_crc, g, np.sign(d)) for g, d in zip(hits.index, hits["diff"])]
    hits = hits.sort_values("t")
    hits.index.name = "gene"
    hits.to_csv(a.out / "crc_hits_with_esca_replication.csv", float_format="%.5g")
    crc.index.name = esc.index.name = "gene"
    crc.to_csv(a.out / "crc_all_genes.csv", float_format="%.5g")
    esc.to_csv(a.out / "esca_all_genes.csv", float_format="%.5g")

    rep = hits[hits.replicated]
    summ["n_crc_hits"] = int(len(hits))
    summ["n_replicated"] = int(len(rep))
    summ["replicated_by_direction"] = rep.direction.value_counts().to_dict()
    summ["replicated_lomo_unstable"] = int((~rep.lomo_stable_crc).sum())
    cls = {}
    for k, genes in CLASSES.items():
        present = [g for g in genes if g in crc.index]
        cls[k] = {"n_present": len(present), "crc_hits": sorted(set(present) & set(hits.index)),
                  "replicated": sorted(set(present) & set(rep.index)),
                  "mean_diff_crc": float(crc.loc[present, "diff"].mean()) if present else None,
                  "mean_diff_esca": float(esc.loc[[g for g in present if g in esc.index], "diff"].mean()) if present else None}
    summ["classes"] = cls
    (a.out / "summary.json").write_text(json.dumps(summ, indent=2))
    print(json.dumps({k: v for k, v in summ.items() if k != "classes"}, indent=2))
    print(json.dumps(cls, indent=1))


if __name__ == "__main__":
    main()
