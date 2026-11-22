# Report: methods, numbers, limitations

**Nonclinical, unreviewed computational research.** Run 2026-11-22, Python 3.14, numpy 2.4,
pandas 3.0, scipy 1.17.

## Data (all downloaded 2026-11-22; MD5 verified against the repository API)

| Source | Files | Licence | Used for |
| --- | --- | --- | --- |
| DepMap 24Q4 Public, doi:10.25452/figshare.plus.27993248.v1 | CRISPRGeneEffect.csv (md5 6edf7ade…), OmicsExpressionProteinCodingGenesTPMLogp1.csv (71794802…), Model.csv (675210d1…) | CC BY 4.0 | Discovery (Chronos gene effect, 1,103 lines with expression and lineage) |
| DEMETER2, figshare 6025238 (McFarland et al. 2018) | D2_combined_gene_dep_scores.csv (e82566a4…) | CC BY 4.0 | RNAi replication (541 lines with expression) |
| PRISM Repurposing Public 24Q2, figshare 25917643 | Extended primary matrix (042b15d3…) and compound list | CC BY 4.0 | Drug-level test of H1 |
| Sanger Cell Model Passports, Project Score combined Sanger v2/Broad 21Q2 (2025-06-24), model list 2026-10-13 | fold-change matrix, model list | Sanger DepMap data usage policy (research use, no redistribution) | Independent CRISPR replication; **only Sanger-screened, QC-passed columns**; 257 models mapped to DepMap IDs; **no per-line values committed** |
| HGNC complete set; Ensembl REST paralogues | | public | Gene groups, EC numbers, paralogy |

The DepMap portal download API is behind an interactive verification page and was **not** used;
the identical public release on figshare was used instead.

## Screen design

* Lines: CRISPR ∩ expression ∩ lineage annotation (n=1,103).
* A (dependency) genes: ≥10 lines with effect < −0.5 and ≥10 lines with effect > −0.3, ≤5% missing
  (3,954 genes). B (biomarker) genes: ≥20 lines with log2(TPM+1) < 1 and ≥20 with > 3 (7,194 genes).
* Covariates: intercept, lineage indicators (lineages with <10 lines pooled as "Other"),
  screen strength (mean effect of the 300 most depleted genes), and 20 expression principal
  components computed from the 5,000 most variable genes after removing lineage means.
* Predictor (v2): indicator that B is expressed (log2(TPM+1) ≥ 1). Partial r between
  residualised A effect and residualised indicator; positive = B-off lines more dependent.
* Discovery threshold: one-sided Bonferroni α=0.05 over 28,445,076 tests (r ≥ 0.1803, df 1,056).
* Leave-one-lineage-out: min r over dropping each lineage with ≥10 lines must be ≥ 0.5 × r.
* RNAi replication: same model on DEMETER2, one-sided, BH 5% over 279 testable hits.
* Penetrance: (fraction of B-off lines with effect < −0.5) − (same for B-on) ≥ 0.20 and
  B-off fraction ≥ 0.25.
* Classes: Ensembl paralog or shared HGNC gene group → *paralog*; B on chrY → *chrY loss proxy*;
  else *non-paralog candidate*.

## v1 (retained negative result)

Continuous expression, no expression PCs: 29,356 hits; 12,104 replicated in RNAi; replicated
hits were dominated by state-marker B genes (TGFBI 98 partners, ITGA3 82, HRH1 76, GBP3 73,
IFI16 68 …; `derived/v1_hub_B.csv`) and adhesion/OXPHOS A genes; benchmark pairs SMARCA2←SMARCA4,
ARID1B←ARID1A, VPS4A/B, ENO2←ENO1, MAT2A←MTAP were not recovered. Lineage adjustment alone does
not remove within-lineage cell-state confounding. v2 was designed **after** this failure; the
18-pair benchmark (`pipeline.py: BENCHMARK`) was fixed before v2 ran.

## v2 screen numbers (`derived/run_summary.json`)

296 discovery hits; 279 testable in RNAi; 86 replicated (BH 5%); 64 replicated and
lineage-robust. Benchmark: 4/18 pairs screenable under the expression filters (most benchmark
B genes are lost by mutation/deletion with residual mRNA, not transcriptional silencing);
4/4 screenable pairs discovered (FAM50A←FAM50B r=0.736, RPP25L←RPP25 0.658, CDS2←CDS1 0.328,
PRMT5←MTAP 0.249); 1/4 replicated in RNAi. Sanger replication (BH 5%): paralog 38/42,
chrY proxy 18/18, non-paralog candidate 95/132 testable.

Strict queue (`queue.py`): 51 pairs; 33 carry a rule-based known/proxy explanation; 18 do not.

## Hypothesis tests (`HYPOTHESES.md` written before each step)

| Test | Statistic | Value |
| --- | --- | --- |
| H1 CRISPR set (15 held-out enzymes) | mean partial r; empirical p (10,000 random sets) | 0.106; p ≤ 1e-4 |
| H1 RNAi | mean r over set members present | 0.064 |
| H1 adherent / suspension | set incl. seeds, mean r (lines; GJA1-off) | 0.144 (890; 313) / 0.021 (156; 122) |
| H1 controls (GJA1) | ribosome / proteasome / spliceosome / Pol II mean r (p) | 0.004 (0.27) / −0.023 (1.00) / 0.002 (0.46) / 0.005 (0.37); all pass |
| H1 PRISM primary (frozen) | 60 drug conditions annotated to the enzymes, mean r | 0.006 (4 nominal p<0.05) |
| H1 same adherent lines | CRISPR set mean r on 625 lines with PRISM | 0.148 (14/15 positive) |
| H1 PRISM active-drug sensitivity (post hoc) | 28 active conditions, mean r, condition bootstrap | 0.032 [0.018, 0.044] |
| H1 paired line bootstrap (1,000; 28 conditions) | CRISPR − drug; drug alone | [0.044, 0.178]; [−0.022, 0.090]; median drug/CRISPR ratio 0.24 |
| H1 Sanger | set mean r; p (106 GJA1-off models) | 0.135; ≤ 1e-4 |
| H1 enzyme enrichment | EC enzymes vs others among 1,572 commonly depleted genes | median 0.0072 vs −0.0040; Mann-Whitney p=1.1e-4; top-50 enzyme fraction 0.52 vs 0.11 base |
| H2 CRISPR set (16 FA/ICL genes) | mean r; p | 0.085; ≤ 1e-4 (15/16 positive); MCM2-7 control −0.006 (p 0.67) |
| H2 RNAi | mean r | −0.001 (**not replicated**) |
| H2 Sanger | mean r; p (30 ALDH2-off models) | 0.111; ≤ 1e-4; MCM2-7 control −0.045 (p 0.96) |
| H3 CRISPR set | mean r; p | −0.014; 0.81 (**falsified**) |

## Correction (2026-11-22, same day)

The PRISM compound list is **not** stored in the same row order as the PRISM data matrix. The first
versions of `pathways.py`, `followup.py` and `bootstrap_h1.py` paired matrix rows with list rows by
position, so the H1 drug-level statistics were computed on the wrong compounds. The bug was found
when the calibration of the follow-up compound-mechanism screen (research/prism-2024-moa-deconvolution)
came out near chance. All three scripts now align metadata to matrix rows by Broad ID. Every H1
drug-level number above is from the corrected rerun. First-version values were: primary 0.001;
active 0.040 [0.005, 0.075] over 11 conditions; paired difference [0.064, 0.174], drug alone
[−0.006, 0.060], ratio 0.18. The non-PRISM results reproduced exactly. The qualitative conclusion
is unchanged: drug-level attenuation holds and a drug-level effect is not established. The numbers
quoted inside HYPOTHESES.md Addenda B–C are the first-version values, kept as written at the time.

## Limitations

* "B off" (TPM < 1) mixes deletion, promoter silencing and lineage non-expression; copy number
  was not used. Most classic paralog SL pairs (SMARCA4, ARID1A, VPS4B, ENO1) are lost by
  mutation or deletion with mRNA still present, so they are invisible to this screen.
* RNAi and Sanger replications reuse many of the same cell lines. They replicate perturbation
  technology, library and laboratory, not tumour sampling.
* RNAi (DEMETER2) has known weaker signal for many genes; RNAi non-replication is reported, not
  hidden. H2 rests on two CRISPR datasets.
* H1 alternative explanations considered: generic growth/quality confounding (argued against by
  null protein-machine controls and the screen-strength covariate); a GJA1-linked cell-intrinsic
  metabolic state (argued against by the paired drug attenuation, but a small drug effect is not
  excluded). PRISM is itself a pooled assay, but every cell is exposed to the drug. No
  experiment was performed. The cross-feeding mechanism is an inference.
* Post hoc elements are declared: v2 design after v1; H1–H3 were derived from v2 hits and tested
  leave-seed-out; Addenda A–C were sequential, and the active-drug restriction is post hoc.
* Literature novelty has not yet been checked. Some of these may be published.
