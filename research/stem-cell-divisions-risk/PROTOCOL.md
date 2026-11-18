# Frozen analysis protocol

**Nonclinical computational research; no medical advice.**

Protocol date: **2026-11-17, America/New_York**. The first commit containing
this file is to be pushed before running any regressions, cross-validation,
inverse-rate calculations, or outcome-based model selection for this project.
The final README will identify that commit and any amendments. This is a
prospective analysis specification for a reanalysis of already-public data,
not a blinded study or an independent clinical-trial registration. Published
correlations, familiar outliers, and the literature were known beforehand.
Source acquisition, rights checks, transcription, and arithmetic auditing may
precede this freeze. None of those checks is a fit.

## Question and literature-driven pivot

Noble et al. (2015), Little et al. (2016), the Tomasetti--Vogelstein reply
(2017, PLOS ONE e0175535), Nowak--Waclaw (2017), and Nunney--Thai (2020)
already address subunit slopes, heterogeneous tissue protection/step counts,
and the failure of cross-tissue correlation to identify attributable causes.
Rediscovering those points is not the proposed novelty. The candidate new
contribution is an exact, computable **compatibility frontier** for a specified
sequential model with bounded tissue-specific per-division rates, including
the coupling of uncertain stem-cell number N and the clock D/N. Priority is
unestablished. The project will retain negative novelty findings.

## Inputs and versions

1. The 31 tissue/cancer-category rows of Tomasetti--Vogelstein (TV) 2015
   Table S1: printed lifetime risk R, stem-cell number N, annual and lifetime
   divisions per stem cell, and total lifetime divisions D. Use a source
   permitting small-aggregate redistribution, with upstream table checking.
2. Little--Hendry--Puskin 2016 S2 public-domain aggregate data as a distinct
   later reanalysis version. Do not call reanalysis changes formal TV errata.
3. A separately identified TV 2017 or permissively licensed published
   reanalysis aggregate panel, if source access and reuse terms permit.
   Do not pool an age-truncated cumulative incidence with lifetime risk, or
   fill unreported N/rates by guessing. Run comparable analyses only for
   complete fields and label endpoint differences. Record inaccessible or
   nonredistributable coverage explicitly.
4. Every numeric-data row must carry source URL, source/version, precise
   locator, access date, reuse basis and reuse URL. Never include person-level
   data. Preserve printed values and explain corrections instead of silently
   repairing them. Hash UTF-8 files after CRLF/CR-to-LF normalization.

## Fixed panels and grouping

The historical replication uses all 31 original rows for log R versus log D.
The **primary mechanistic panel** excludes the six explicitly elevated-risk
categories: colorectal FAP, colorectal Lynch, duodenum FAP, head/neck HPV16,
liver HCV, and lung smokers. This exclusion follows the published subgroup
definitions, before fitting; it does not imply the retained rows have no
environmental or inherited causes. No outlier is removed on a residual rule.

Sensitivity panels are all rows with 0<R<1, and the later source version with
the same exclusion rule. A printed risk of exactly one remains usable in the
historical log R replication but is excluded from log cumulative-hazard and
inverse-rate calculations; do not manufacture a finite hazard for it.

Freeze anatomical groups from labels before fitting: all osteosarcoma rows
form one bone group; pancreatic subtypes one pancreas group; medulloblastoma
and glioblastoma one brain group; AML and CLL one hematopoietic group;
basal-cell skin and melanoma one skin group; thyroid subtypes one thyroid
group. Group each exposure/hereditary subgroup with its corresponding organ.
Duodenum and small intestine share one group because one is a subsite of the
other; holding them out separately could also retain an overlapping input.
Other sites are their own groups. Grouping is conservative and is not a claim
of statistical independence between the remaining groups. Publish every row's
group and inclusion flag.

## Observations, models and fits

Use base-10 logs and H=-log(1-R). The primary operational division clock is
v=D/N. This is an **effective allocation surrogate**, not an assertion that
TV's total developmental divisions measure a terminal lineage's depth. The
prespecified alternative clock is v=L+log2(N), with L the source lifetime
divisions per stem cell; this is a balanced-development path proxy, also not
an observed lineage history. Models assume a common observation horizon and
direct onset risk; real lifetime risks can violate these assumptions.

Fit these fixed model families, with equal row weights in each training set:

| Model | Training rule |
|---|---|
| Historical descriptive regression | OLS log R = a+b log D |
| Hazard descriptive regression | OLS log H = a+b log D |
| Separable descriptive regression | OLS log H = a+bN log N+bv log v |
| Homogeneous sequential rare-event family, k=1,...,8 | log H = a+log N+k log v-log(k!); fit only a by its residual mean |
| One expanding intermediate, k=2,...,8 | Add log Phi_k(gamma v), Phi_k(z)=sum(n>=0) k! z^n/(n+k)!; fit amplitude and select gamma by training SSE |

For the expansion grid, set gamma=c/max(v_training), with
c in {0, 0.25, 0.5, 1, 2, 4, 8}. Use the same fitted gamma for held-out rows;
never use their clocks to set the training grid. Break exact SSE ties toward
smaller c. This is a finite-grid comparison in a proved rare-marking limit,
not optimization over all clonal-expansion biology. At a fixed k, expansion
in calendar time with a common horizon is absorbed into amplitude and has
the same limiting predictions as the no-expansion shape.

Do not claim a best k is an identified driver count. The code will check the
finite equal-rate chain at the rates implied by each no-expansion limiting
fit and report the maximum log-hazard approximation discrepancy. The fits
remain explicitly fits of the limiting family, whether that check is small
or large. The exact-chain compatibility calculation below does not use the
rare-event approximation.

For the alternative developmental clock, repeat the two descriptive hazard
models, all fixed-k no-expansion fits and their prediction scores, and exact
point-input rate frontiers. Do not repeat the expansion grid or rectangular
N/D uncertainty calculation for this alternative clock. The historical log-R
regression depends only on D and is unchanged. The saturation-only test is
the k=1 model H=uD fitted as above, plus its fitted log-R versus log-D slope
and the attenuation bound proved for an ordered N/v cloud; explicitly test
whether that ordering condition holds before using the bound.

## Prediction and uncertainty

Re-fit everything in each training fold. Report both leave-one-row-out and
**leave-one-anatomical-group-out** predictions. The primary score is the
square root of the mean, across held-out groups, of each group's mean squared
log H error. Also report ordinary row-weighted RMSE. No held-out outcome is
used for k/gamma selection. Report each fixed k; an optional selected-k
score must choose k by training SSE inside the fold and be labeled as such.
The implementation will report fixed-k scores only, avoiding an extra
selected-k estimand. A rank-deficient OLS design is reported unavailable;
it is never repaired by deleting predictors. For bootstrap intervals of the
univariate slope, omit a draw only if log D is constant and report the count
of such draws. Do not perform a multivariable-coefficient bootstrap.

For descriptive slope sensitivity use 2,000 fixed-seed (seed 1260825)
anatomical-group bootstrap draws, retaining all rows of each sampled group.
Report the 2.5th and 97.5th percentiles with linearly interpolated empirical
quantiles. These are conditional resampling intervals, not registry sampling
confidence intervals or a probability of a biological hypothesis. Do not
invent denominators, standard errors for N or D, or a population sampling
design. All source-version comparisons and deterministic uncertainty
scenarios must be shown, including disagreements.

## Exact compatibility certificate

For each k=1,...,8, the fixed independent-lineage model is
R=1-Q_k(u v)^N, Q_k(x)=exp(-x) sum(j=0,...,k-1) x^j/j!, with a continuous
per-division transition **intensity** u>0. All k waiting stages have this
intensity. This is not a Bernoulli probability or a universal model of tissue
development. Derive and numerically invert psi_k=-log Q_k, giving
u_i=psi_k^{-1}(H_i/N_i)/v_i.

Report the exact smallest symmetric rate-heterogeneity factor rho allowing
all u_i in [u0/rho,u0*rho] for some common u0>0, its limiting-form comparison,
and the tissue(s) attaining the upper and lower constraints. Also report
compatibility for a prespecified allowance rho in {1,2,5,10,100}.

Use two uncertainty scenarios for the effective clock certificate:

* Exact printed point inputs.
* An explicitly hypothetical rectangular stress test:
  R in [R/1.25,min(1-1e-12,1.25R)], N in [N/2,2N], D in [D/2,2D].

The stress test is not a confidence region. D and N are independent interval
coordinates, but v=D/N must retain their algebraic coupling in the inverse.
This box can over-cover source-coherent biological histories. Prove its
extremizers rather than sampling corners and assuming they are exhaustive.
For integer-target models, source estimates are operational real-valued
scales; distinguish the continuous extension from literal integer counts.

## Falsification, outliers and interpretation

1. State which homogeneous models fail exact compatibility and which merely
   predict poorly. Do not turn lack of sampling uncertainty into a p-value.
2. Find exact duplicate (N,D) source inputs across subgroup rows. A homogeneous
   same-k/same-rate model predicts equal risks. Report their observed risk
   contrasts, finite inverse-rate contrasts where defined, and why these
   falsify the identical-parameter specification without identifying causes.
3. Name the three largest positive and three largest negative held-out
   log-H residuals from the primary descriptive hazard regression; use the
   group-held-out predictions, deterministic tissue-ID tie breaking, and no
   residual-based deletion. The implied hazard-amplitude factor is exactly
   H_observed/H_predicted (use its reciprocal to describe a deficit). This
   ratio is not automatically a probability or an environmental fraction.
4. State at least one mechanistic hypothesis with an independently measurable
   prediction and a condition that would refute it. A fitted free tissue
   multiplier alone is explanatory bookkeeping, not a tested mechanism.
5. Test saturation-only, alternative clock, source revisions, high-risk-row
   exclusion, grouping and expansion as attempts to defeat the preferred
   interpretation. Preserve all negative findings. Distinguish mutation
   attribution, case attribution and preventability.

## Verification and publication

Numbered theorems require all assumptions and complete proofs, including
general-parameter or proved-limit statements. Finite computational examples
do not prove a theorem. Check stable inverses, Taylor bounds, exact algebra,
inequality directions, edge cases and independently implemented numerical
identities. One new stdlib-only test file must reproduce every reported
empirical number and normalized data hash in under 60 seconds offline.

All changes are confined to research/stem-cell-divisions-risk/ and
tests/test_research_stem_cell_divisions_risk.py. Run existing repository gates
without modifying them. Hostile self-review must examine proofs, source
transcriptions, model implementation, scope and novelty before the final
push. Open a PR to main and report its actual CI status. No clinical use,
medical advice, biological validation or independent peer review is claimed.
