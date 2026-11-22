# Frozen pathway-level hypotheses (written before running `pathways.py`)

Written after the v2 pair screen and its annotation (`derived/discovery_hits_annotated.csv`),
**before** any of the tests below were computed. Each hypothesis generalises one or more
screen hits to a pre-declared gene set, so the test is no longer the pair that produced it.
The seed pairs are excluded from every set-level statistic ("leave-seed-out").

Convention: predictor = 1 if B is expressed (log2(TPM+1) >= 1), 0 if off. Positive partial r
(after lineage, screen-strength and 20 expression-PC adjustment) = B-off lines are more
dependent on A.  Set statistic = mean partial r over set members screened; null = 10,000
random gene sets of equal size from all screened CRISPR genes (same B).  Pass = one-sided
empirical p < 0.01 **and** the mean is positive in RNAi (DEMETER2) over set members present.

| ID | B (loss) | Seed hits | Pre-declared A gene set (seeds excluded from test) | Extra falsification test |
|----|----------|-----------|-----------------------------------------------------|--------------------------|
| H1 | GJA1 (connexin-43) | PFAS, ATIC, CTPS1, YRDC | De novo nucleotide synthesis: PPAT, GART, PAICS, ADSL, IMPDH2, GMPS, ADSS2, CAD, DHODH, UMPS, CTPS2, PRPS1, TYMS, DHFR, RRM1, RRM2 | **Pooled-screen cross-feeding test.** If the signal is gap-junction transfer of nucleotides from neighbouring cells in pooled CRISPR screens, it (a) must be present in adherent lines and weaker/absent in suspension lines, and (b) must be **absent** for drugs that inhibit the same enzymes in PRISM (every cell is inhibited, no healthy donors). A GJA1 effect of the same sign for those drugs falsifies the artefact explanation and supports a cell-intrinsic mechanism instead. |
| H2 | ALDH2 | BRIP1 | Fanconi anaemia / ICL repair: FANCA, FANCB, FANCC, FANCD2, FANCE, FANCF, FANCG, FANCI, FANCL, FANCM, UBE2T, FAAP24, FAAP100, SLX4, ERCC4, ERCC1 | Biology: ALDH2 loss raises endogenous acetaldehyde, which forms DNA interstrand crosslinks (Aldh2/Fancd2 double-mutant mice; Langevin et al. 2011). Prediction: FA pathway dependency in ALDH2-off lines; negative control set = homologous-recombination-independent essential DNA replication genes (MCM2-7) shows no shift. |
| H3 | ACOT9 | SUCLG1, SUCLA2 | Succinyl-CoA synthetase: SUCLG2 (only remaining subunit) + succinyl-CoA producing/consuming neighbours OGDH, DLST, DLD | Weak by construction (small set). Report the pair-level replication of SUCLG1/SUCLA2 in Sanger Project Score as the main independent test. |
| H4 | PUDP | ITPA | Pair only (no set): replicate in Sanger Project Score | Both are noncanonical-nucleotide "house-cleaning" enzymes; PUDP lies in the Xp22.31 (STS) region. |

Additional, independent replication for every non-paralog candidate pair: Sanger Project
Score (Sanger v2 screens only; different library and laboratory), one-sided partial r with
the same covariate scheme, BH FDR 5% over candidates. Sanger data are used under their
research-use licence; **no Sanger per-line values are committed**, only summary statistics.

## Addendum A (written after the first `pathways.py` run, before `followup.py`)

First-run outcome: H1 passes (CRISPR set p<1e-4, RNAi mean positive), with the predicted
adherent/suspension split and no PRISM antimetabolite signal; H2 is coherent in CRISPR but not
in RNAi; H3 fails. The following checks were fixed before running `followup.py`:

1. **H1 specificity (negative controls).** Gap junctions pass small metabolites (< ~1 kDa), not
   proteins. GJA1-off must **not** shift dependency on essential protein machines whose
   loss cannot be rescued by a neighbour: cytosolic ribosome (RPL/RPS), proteasome
   (PSMA1-7, PSMB1-7), spliceosome core (SF3B1-5, SNRPD1-3, PRPF8), RNA polymerase II (POLR2A-L).
   Pass = each control-set mean partial r is below half of the H1 set mean and its empirical
   p > 0.01.
2. **H1 same-line contrast.** Recompute the CRISPR H1 set statistic on exactly the adherent
   lines that also have PRISM data, so that the CRISPR-vs-drug difference is not a
   line-composition artefact.
3. **H1 PRISM activity sensitivity (declared post hoc).** Restrict PRISM to drug conditions
   whose target annotation is a human nucleotide-synthesis enzyme **and** which are active
   (>= 10% of lines with log2 fold-change < -1). Report the mean r with a bootstrap 95% CI over
   drug conditions. This is a sensitivity analysis, not a replacement for the frozen primary.
4. **Sanger set-level replication of H1 and H2** using the same set statistic (members present
   in Sanger, seeds excluded, 10,000 random equal-size sets).

## Addendum B (written after `followup.py`, before `bootstrap_h1.py`)

Active-drug sensitivity gave a small positive GJA1 association (mean r 0.040, drug-condition
bootstrap 95% CI 0.005-0.075) versus 0.148 for CRISPR on the same 625 adherent lines. Because
the drug-condition bootstrap ignores the shared cell lines, the decisive comparison is a
**paired bootstrap over cell lines** (1,000 resamples, both statistics recomputed on each
resample). Decision rule fixed now: "attenuated at drug level" if the 95% CI of
(CRISPR set mean r - active-drug mean r) excludes 0; "cell-intrinsic component present" if
the 95% CI of the drug mean r excludes 0; both can hold.

## Addendum C (written after `bootstrap_h1.py`, before `enzyme_enrichment.py`)

Paired line bootstrap: CRISPR minus active-drug mean r 95% CI [0.064, 0.174] (attenuated at drug
level); drug mean r CI [-0.006, 0.060] (cell-intrinsic component not established). The
genome-wide GJA1 ranking (`derived/gja1_genomewide_partial_r.csv`) was then inspected; its top is
dominated by small-metabolite enzymes. To test this as a pattern rather than by eye, fixed now:
among **commonly depleted genes** (mean CRISPR effect < -0.5 over the screened lines), genes with
an HGNC EC number (`enzyme_id`) have higher GJA1 partial r than genes without one
(one-sided Mann-Whitney, p < 0.01), and the median difference is reported. Mitochondrially
encoded/translation machinery is not excluded, so this test is conservative for the enzyme class.

> **Correction note:** the PRISM numbers quoted in Addenda B and C came from a row-misaligned compound
> list. Corrected values are in REPORT.md, section "Correction". The decision rules above were not changed.
