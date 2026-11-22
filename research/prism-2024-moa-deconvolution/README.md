# PRISM 24Q2 compounds vs CRISPR knockouts: mechanism-of-action screen

**Nonclinical, unreviewed computational research. No treatment advice.** Record-level public
cell-line data: PRISM Repurposing Public 24Q2 (figshare 25917643), DepMap 24Q4 (doi:10.25452/figshare.plus.27993248.v1),
DEMETER2 (figshare 6025238), all CC BY 4.0. Sanger Project Score is used for summary statistics only,
under its research-use licence. No literature priority check has been made, so nothing here is claimed as new.

## Verdict: the frozen calibration gate FAILED, so no discovery is claimed

`PROTOCOL.md` was frozen before any statistic. It required the annotated target to be in the top 10 of
17,205 knockout profiles for ≥10% of annotated active compounds, and to be ≥20× chance.

| | value |
| --- | --- |
| active compound conditions (≥20 lines with LFC < −1) | 2,325 |
| annotated and evaluable | 1,457 |
| annotated target ranked 1st | 3.7% |
| annotated target in top 10 | **8.0%** (chance 0.16%, **50×**) |
| gate (≥10% and ≥20×) | **fail** |

The specificity is far above chance, but the frozen bar was not met, so the discovery list below is
**exploratory**. An earlier run gave 0.4% because the PRISM compound list is not stored in matrix
row order. The calibration exposed that bug, and it also affected research/expression-loss-synthetic-lethality
(corrected there in PR #94).

## Why the gate misses: annotation, not signal (post hoc observation)

Most new-screen compounds (REP.1M/REP.300) carry **no** `repurposing_target`, so they cannot count
towards calibration. Yet the screen recovers their known pharmacology, with replication in independent
Sanger CRISPR **and** RNAi:

| compound (target field empty) | top knockout match | DepMap r | Sanger r | RNAi r |
| --- | --- | --- | --- | --- |
| milademetan | MDM2 | 0.62 | 0.44 | 0.55 |
| MIK665 / AMG-397 / AZD5991 / AMG-176 | MCL1 | 0.47 / 0.47 / 0.45 / 0.35 | 0.50 / 0.47 / 0.53 / 0.28 | 0.29 / 0.35 / 0.33 / 0.26 |
| MRTX1133 | KRAS | 0.33 | 0.37 | 0.26 |
| mobocertinib / pyrotinib | ERBB2 | 0.31 / 0.30 | 0.43 / 0.39 | 0.32 / 0.30 |
| GDC-0077 (inavolisib) | PIK3CA | 0.28 | 0.26 | 0.27 |
| RMC-4550 (MOA field says "SHIP2 phosphatase inhibitor") | **PTPN11 (SHP2)** | 0.35 | 0.20 | 0.35 |

RMC-4550 is a SHP2 (PTPN11) inhibitor, and the knockout match corrects the dataset's MOA label.

## Exploratory, unclaimed leads (priority not checked)

236 new-compound pairs pass per-compound Bonferroni with a non-annotated top gene. 124 replicate
(Sanger 102, RNAi 67, both 45; BH 5%). Most are recoveries of known targets like those above. Leads
that are not obvious recoveries:

* **Selpercatinib (RET inhibitor) ↔ heparan-sulfate synthesis and FGFR2** (EXT1 0.28, SLC35B2 0.26, EXT2 0.23,
  HS2ST1 0.22, FGFR2 0.20). Heparan sulfate is an FGF co-receptor, which suggests FGFR-pathway
  activity at 2.5 µM. Replicated in Sanger (0.26), not in RNAi.
* **Zanubrutinib (BTK inhibitor) ↔ ERBB2/EGFR** (0.23; Sanger 0.33, RNAi 0.27). This fits an
  off-target HER2/EGFR effect at 2.5 µM, of the kind known for other covalent BTK inhibitors.
* **PBT-1033 (metal chelator) ↔ iron-sulfur cluster assembly** (HSCB 0.22, IREB2, LYRM4; Sanger 0.30,
  RNAi 0.16).
* **Ruboxistaurin (PKC inhibitor) ↔ p53-wild-type signature** (PPM1D 0.24, MDM2 0.16; RNAi 0.13).
* **NIM811 (cyclophilin inhibitor) ↔ SLC33A1**, the ER acetyl-CoA transporter (0.25; Sanger 0.33, RNAi 0.24).
* **Microtubule agents ↔ DARS1.** 7 of 9 compounds whose top match is DARS1 are tubulin binders
  (eribulin, DM4, soblidotin, mivobulin, batabulin, crolibulin, ON-01500). This may be a
  mitotic-rate proxy that screen strength does not absorb; it is unexplained.

## Files
`PROTOCOL.md` (frozen), `moa.py` (numpy/pandas/scipy; data not in repo), `derived/summary.json`,
`derived/drug_gene_similarity.csv` (every active condition: top-5 genes and annotated-target rank),
`derived/new_compound_discoveries.csv` (exploratory pairs with replication). Test:
`tests/test_research_prism_2024_moa_deconvolution.py` (stdlib, offline).
