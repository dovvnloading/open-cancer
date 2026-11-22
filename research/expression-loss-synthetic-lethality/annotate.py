"""Annotate v2 hits: Ensembl paralogy, chromosome of B, penetrance, candidate class.

Run after pipeline.py:  python annotate.py --data <dir> --out derived
Ensembl REST (https://rest.ensembl.org) paralogue calls are cached in
derived/ensembl_paralogues.json so reruns are offline and deterministic.
"""
from __future__ import annotations
import argparse, json, time, urllib.request
from pathlib import Path
import pandas as pd

MIN_PENETRANCE_GAP = 0.20      # frac dependent (B off) - frac dependent (B on)
MIN_FRAC_DEP_OFF = 0.25

def paralogues(sym, cache):
    if sym in cache:
        return cache[sym]
    url = f"https://rest.ensembl.org/homology/symbol/human/{sym}?type=paralogues;format=condensed;content-type=application/json"
    ids = []
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                d = json.load(r)
            ids = sorted({h["id"] for x in d.get("data", []) for h in x.get("homologies", [])})
            break
        except Exception as e:  # 400 = unknown symbol
            if "400" in str(e):
                break
            time.sleep(2)
    cache[sym] = ids
    time.sleep(0.08)
    return ids

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path(__file__).with_name("derived")); a = ap.parse_args()
    h = pd.read_csv(a.out / "discovery_hits.csv")
    hg = pd.read_csv(a.data / "hgnc_complete_set.txt", sep="\t", dtype=str,
                     usecols=["symbol", "location", "ensembl_gene_id"])
    loc = dict(zip(hg.symbol, hg.location)); ens = dict(zip(hg.symbol, hg.ensembl_gene_id))
    cpath = a.out / "ensembl_paralogues.json"
    cache = json.loads(cpath.read_text()) if cpath.exists() else {}
    from concurrent.futures import ThreadPoolExecutor
    todo = [g for g in sorted(set(h.A) | set(h.B)) if g not in cache]
    with ThreadPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(lambda g: (g, paralogues(g, {})), todo))
    cache.update(dict(res))
    cpath.write_text(json.dumps(cache, indent=0, sort_keys=True))
    def is_par(A, B):
        return (ens.get(B) in set(cache.get(A, []))) or (ens.get(A) in set(cache.get(B, [])))
    h["ensembl_paralog"] = [is_par(A, B) for A, B in zip(h.A, h.B)]
    h["B_location"] = [loc.get(b, "") for b in h.B]
    h["B_on_chrY"] = h.B_location.fillna("").str.startswith("Y")
    h["penetrance_gap"] = h.frac_dependent_b_off - h.frac_dependent_b_on
    pen = (h.penetrance_gap >= MIN_PENETRANCE_GAP) & (h.frac_dependent_b_off >= MIN_FRAC_DEP_OFF)
    cls = []
    for _, r in h.iterrows():
        if not pen[_]:
            cls.append("weak_or_non_penetrant")
        elif r.ensembl_paralog or r.same_hgnc_group:
            cls.append("paralog")
        elif r.B_on_chrY:
            cls.append("chrY_loss_proxy")
        else:
            cls.append("non_paralog_candidate")
    h["class"] = cls
    h["evidence_tier"] = [
        "replicated_rnai+lineage_robust" if (r.replicated and r.lineage_robust) else
        "lineage_robust_rnai_not_replicated" if (r.lineage_robust and r.n_rnai == r.n_rnai) else
        "lineage_robust_rnai_untestable" if r.lineage_robust else "lineage_fragile"
        for r in h.itertuples()]
    h.to_csv(a.out / "discovery_hits_annotated.csv", index=False, float_format="%.6g")
    print(h.groupby(["class", "evidence_tier"]).size().to_string())

if __name__ == "__main__":
    main()
