"""Build the strict candidate queue from derived/candidates_with_sanger.csv.

Strict tier: class == non_paralog_candidate, lineage_robust, and replication either in
RNAi (BH q<=0.05) or in Sanger Project Score (BH q<=0.05 with >=5 B-off Sanger models).
Rule-based proxy labels mark hits explained by already-known biology.
"""
from pathlib import Path
import pandas as pd

OUT = Path(__file__).with_name("derived")
CTA_PREFIX = ("SSX", "MAGE", "XAGE", "GAGE", "PAGE", "CTAG", "VCX", "SAGE", "DCAF4L2", "COX7B2",
              "CSMD1", "NPFFR2", "FAM133A", "C5orf58")

def proxy(a, b):
    if a == "NAA10" and b.startswith(CTA_PREFIX):
        return "known: co-activation with NAA11 (NAA10 paralog) on DNA hypomethylation"
    if a == "WRN":
        return "known: MSI / MLH1-promoter silencing marker (WRN dependency in MSI)"
    if b == "KLHL9":
        return "known: 9p21.3 co-deletion with FOCAD (FOCAD loss -> HBS1L/PELO)"
    if b == "GJA1":
        return "H1: gap-junction cross-feeding in pooled screens (see REPORT)"
    return ""

def main():
    c = pd.read_csv(OUT / "candidates_with_sanger.csv")
    strict = c[(c["class"] == "non_paralog_candidate") & c.lineage_robust.astype(bool) &
               ((c.replicated.astype(bool)) | ((c.q_sanger <= 0.05) & (c.n_b_off_sanger >= 5)))].copy()
    strict["proxy_or_known_explanation"] = [proxy(a, b) for a, b in zip(strict.A, strict.B)]
    cols = ["A", "B", "B_location", "r_crispr", "n_b_off", "frac_dependent_b_off", "frac_dependent_b_on",
            "lolo_min_r", "r_rnai", "n_rnai", "q_rnai", "r_sanger", "n_sanger", "n_b_off_sanger", "q_sanger",
            "proxy_or_known_explanation"]
    strict = strict[cols].sort_values(["proxy_or_known_explanation", "r_crispr"], ascending=[True, False])
    strict.to_csv(OUT / "strict_candidate_queue.csv", index=False, float_format="%.4g")
    print(len(strict), "strict;", (strict.proxy_or_known_explanation == "").sum(), "without a known/proxy explanation")
    print(strict[strict.proxy_or_known_explanation == ""][["A", "B", "r_crispr", "n_b_off", "r_rnai", "r_sanger", "n_b_off_sanger", "q_sanger"]].round(3).to_string())

if __name__ == "__main__":
    main()
