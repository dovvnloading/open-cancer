# Protocol (frozen before any organoid-vs-cell-line statistic was computed)

**Nonclinical, unreviewed computational research.** Written 2026-11-22.

## Question
Which gene knockouts have systematically different fitness effects in patient-derived tumour
**organoids** than in established **2D cell lines** of the same cancer type, screened by the same
laboratory (Wellcome Sanger Institute)? Genes that are essential only in organoids may mark tumour
biology that cell lines lost in culture. Genes essential only in cell lines may mark culture
adaptation.

## Data (Sanger Cell Model Passports; research-use licence; only gene-level summary statistics are written)
* Organoids: `Organoids_fitness_scores_fold_change_values_20260723.tsv` (QC-passed models).
* Cell lines: `Project_score_combined_Sanger_v2_Broad_21Q2_fitness_scores_fold_change_values_20250624.tsv`,
  **Sanger-screened, QC-passed columns only** (the same laboratory and pipeline).
* Annotation: `model_list_20261013.csv` (tissue, cancer_type, model_type);
  `CRISPR_supplemental_information_20260723.csv` (core_fitness_pancan).

## Normalisation (per model)
Gene-level fold change scaled so that the median of all genes = 0 and the median of pan-cancer core-fitness
genes = −1. This removes differences in screen strength between models. Duplicate screens of one model are averaged.

## Discovery: colorectal (CRC)
Organoids with cancer_type "Colorectal Carcinoma" vs Sanger cell lines with cancer_type "Colorectal Carcinoma".
Per gene: difference in means (organoid − line) and Welch t-test, two-sided; BH FDR over genes with
data in ≥80% of models of both groups. Discovery hit: q ≤ 0.01 and |difference| ≥ 0.3.

## Replication: oesophageal (independent patients and models)
The same statistic for oesophageal carcinoma organoids vs oesophageal carcinoma cell lines. A CRC hit
**replicates** if the oesophageal difference has the same sign with one-sided p ≤ 0.05 after BH over the
discovery hits.

## Interpretation guards (decided now)
* Report how replicated hits split by direction. Then test enrichment of pre-declared classes:
  (a) WNT/R-spondin/EGF niche signalling (organoid medium supplies niche factors: LGR4, LGR5, RNF43, ZNRF3,
  CTNNB1, TCF7L2, EGFR, ERBB2, ERBB3, KRAS, MAPK1, PTPN11, SOS1, GRB2); (b) cell adhesion and integrin genes
  (2D plastic vs 3D matrix: ITGAV, ITGB1, ITGB5, PTK2, TLN1, FERMT2, ILK, PXN, CRK, CRKL, RAC1); (c) oxygen and
  metabolism (HIF1A, VHL, EGLN1, ARNT, PDK1, SLC2A1); (d) interferon/innate immunity (STAT1, IRF1, JAK1, IFNGR1, MAVS).
* Any hit driven by a single model (leave-one-model-out sign change in either group) is flagged.
* Organoid-vs-line differences can reflect medium, matrix, passage, tumour stage, MSI/CRIS composition or
  library batch. These are named limitations, and no causal claim is made. MSI composition is reported for both groups.

## Addendum A (written after the primary run, before `within_platform.py`)

Primary result: 623 CRC hits, 481 replicated in oesophageal models. But both directions are dominated
by core-essential genes (spliceosome more essential in organoids; ribosome more essential in lines),
the organoid−line differences correlate genome-wide between the two cancer types (r = 0.75), and a
multi-copy artefact gene (KRTAP4-11) is a hit. This is the signature of a screening-campaign or
library effect. The primary design **cannot separate biology from platform**, so its hit list is
reported as platform-confounded and not interpreted gene by gene.

Declared rescue design, which is immune to between-platform library effects: **within-platform MSI vs MSS contrast.**
Per gene, MSI − MSS mean difference among CRC organoids, and separately among CRC cell lines (same
normalisation). Positive control fixed now: WRN must be more essential in MSI in **both** platforms
(Welch one-sided p ≤ 0.01 in each). Transfer statistic: Spearman correlation of the two genome-wide
difference vectors, and, among cell-line MSI-specific dependencies (q ≤ 0.05 and difference ≤ −0.3),
the fraction that also shows difference ≤ 0 with one-sided p ≤ 0.05 in organoids. Genes significant
(q ≤ 0.05, difference ≤ −0.3) in organoids but not in lines (p > 0.2) are listed as organoid-only MSI
dependencies. They are exploratory; priority is not checked.

## Addendum B (written after Addendum A, before `wrn_attenuation.py`)

Addendum A: WRN passes the positive control in both platforms, but the MSI−MSS WRN difference is
−0.67 in organoids vs −1.44 in cell lines (ratio 0.46). Genome-wide transfer is underpowered (10 MSI lines).
Question fixed now: **is WRN-in-MSI attenuation in organoids specific, or a general compression of selective
dependencies in organoid screens?** Comparators: oncogene addiction contrasts in CRC with mutation status
from Sanger `mutations_summary_latest` (cancer_driver = t): KRAS (KRAS-mutant vs wild-type), BRAF
(BRAF V600E vs wild-type), CTNNB1 (CTNNB1 driver-mutant vs wild-type), each computed per platform with the
same normalisation. Statistic: organoid/line ratio of the mutant−WT (or MSI−MSS) difference, with a
bootstrap 95% CI over models (2,000 resamples within each group). **WRN-specific attenuation** is supported
only if the WRN ratio CI lies entirely below the CI of every comparator with ≥3 mutant models in each platform.
Otherwise the result is "general compression / not specific". Also reported: the fraction of MSI models with a
WRN effect < −0.5 in each platform.
