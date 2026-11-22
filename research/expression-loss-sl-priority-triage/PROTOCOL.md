# Can deep DNA loss at PUDP predict ITPA dependency?

Nonclinical, unreviewed computational research. Frozen on 2026-11-22 before parsing copy-number records, assigning exposures, joining outcomes or computing the contrast. Feasibility checks covered access, release metadata and column headers. No browser download was requested.

## Choice and estimand

H1 and H2 have direct published precedents. Among the unlabeled pairs, ITPA ← PUDP has the strongest discovery partial correlation (0.4005), a second CRISPR signal (0.4243), and weak RNAi evidence. DepMap has already publicly linked PUDP expression to ITPA predictability. The question worth testing is whether an independently measured DNA lesion predicts a materially larger fitness cost of ITPA knockout. This is one test with a new exposure in previously screened cell lines. It cannot establish genetic causality, exclude regional co-deletion, or establish generalization to new models.

## Data and eligibility

Use Broad DepMap Public 24Q4 v1, CC BY 4.0, DOI [10.25452/figshare.plus.27993248.v1](https://doi.org/10.25452/figshare.plus.27993248.v1). Use only model-indexed relative DNA copy ratios from OmicsCNGene.csv (Figshare file 51065324, original MD5 4696565fa74dc805983bc278b4a41802) and sex and donor grouping from Model.csv (51065297, 675210d17675f3517b0ce39a3c274f16). Relative copy ratios in this release are linear. Absolute PureCN calls exclude sex chromosomes. Use ITPA Chronos effects, lineage and, for descriptive concordance only, PUDP expression from PR #94's committed DepMap extract at 77d53bedf6bf93528f14077a3789f5c8533acb60. Join strictly on ModelID and reject duplicate keys. No Sanger records enter this analysis.

Start with all records in the screen extract. Require finite ITPA effect, nonblank lineage, annotated Female/Male sex, known source donor, finite nonnegative PUDP CN and all five reference X-gene CN values, and positive reference median. Within the eligible records select the lexicographically first ModelID per donor, before classifying exposure. Retain only de-identified cell model accessions, sex and representative status from this grouping. No patient identifiers or clinical records are exported.

## Exposure, outcome and single contrast

For each model divide the PUDP (8226) relative DNA copy ratio by the median for MID1 (4281), KDM5C (8242), HUWE1 (10075), RPS4X (6191) and G6PD (2539) in that same model. These loci span Xp and Xq. This denominator and exact stratification by source sex address baseline X dosage. Deep relative loss is a ratio ≤0.20. A retained locus is 0.75–1.25 inclusive. Exclude intermediate and gain calls. This definition measures regional DNA loss against the rest of X, without assuming absolute zero copies.

Stratify by the exact OncotreeLineage and sex, requiring at least one loss and two retained donors per stratum. Calculate D = mean across loss donors of (mean retained ITPA Chronos effect in that stratum minus that donor's ITPA Chronos effect). D > 0 means greater dependency with loss. Every loss donor has equal weight. Report stratum sizes and contrasts. As descriptive checks only, report in each exposure group of the matched cohort the proportion with PUDP RNA off (log2[TPM+1] <1), and the fractions with deep loss of each of the adjacent genes STS and PNPLA4 by the same X reference. Exclude missing values separately from each descriptive denominator. These checks cannot change the verdict.

## Frozen decision rule

The sharp null is no effect of loss conditional on the strata. The minimum material decrement is 0.20 Chronos units, chosen as one fifth of the approximate nonessential-to-common-essential reference interval. Require at least eight distinct matched loss donors across at least two strata, otherwise report INCONCLUSIVE (insufficient exposure). Missing or corrupt release inputs, ambiguous keys or schema failures are data failures, never negative evidence.

Use 19,999 random within-stratum permutations (Python Random seed 9401, sorted strata and models), preserving group sizes. For the sharp null, obtain the upper-tail p0 by comparing permuted D to observed D, inclusive. For the 0.20 constant-decrement boundary, add 0.20 to each observed loss outcome, then obtain the lower-tail p0.20 using the same label permutations and the shifted observed D. Use the add-one correction and a 1e-12 tolerance at ties. Both calculations invert the same contrast under alternative sharp nulls. They assume conditional exchangeability and the boundary calculation assumes a constant additive effect.

- SUPPORTED if D ≥0.20 and p0 ≤0.05.
- FALSIFIED for the prespecified material DNA-loss prediction if D <0.20 and p0.20 ≤0.05.
- Otherwise INCONCLUSIVE. Failure to reject zero alone does not falsify the material prediction.

Do not revise the thresholds, switch genes, add outcomes or rerun a more favorable subgroup after unblinding. Screen selection, model reuse, measurement error, chromosome-arm changes, regional passengers, residual subtype or screen quality can affect interpretation. A result for this lesion does not resolve all mechanisms of transcriptional silencing.
