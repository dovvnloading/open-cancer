# Frozen protocol: transport of published CH fitness to age prevalence

**Nonclinical, unreviewed research. No medical advice.**

Protocol drafted on 2026-11-18 (America/New_York), from current `main`
`5993e6cc474958bd7645721ce1771b6bc7848f08`. This file and its parameter
CSVs must be committed **before target prevalence extraction or scoring** and
remain unchanged thereafter. The later README will identify the commit and
any deviations. This is an in-session frozen analysis, not a registered study.

## Prior art, exposure, and the question

Watson et al., Science 2020, DOI
[10.1126/science.aay9333](https://doi.org/10.1126/science.aay9333), already
predicted age- and detection-threshold-dependent prevalence (Figure 3) and
checked the combined DNMT3A R882H/C age trend (Figure S18). Jaiswal 2014 and
Genovese 2014 contributed to its development collection and cannot be
independent holdouts. Watson estimated a distribution of fitness effects
(DFE) within genes, not one interchangeable whole-gene scalar. Fabre 2022,
Mitchell 2022, and newer competition and CH-prediction work further preclude
claiming the general idea as new. Bibliographic priority is unestablished.

The pivot is a reproducible **unchanged-parameter external transport audit**:
do published fitness distributions and mutation supplies predict independently
published, matched age-prevalence aggregates? If the measurement bridge cannot
be established, test the closest declared descriptive transport and prove the
precise identification obstruction. Do not rename a proxy a biological test.

Before this commit we read the research index, DATA_POLICY, PR #68's model,
proofs and coefficients, Watson's model-development article/notebooks, and
literature abstracts/methods. Broad search snippets incidentally exposed a UK
Biobank overall prevalence/sample-size/age-range summary, Fabre overall and
oldest-age prevalence summaries, and Bick's abstract totals. These sources are
not fully blinded. No independent target age table has been deliberately
extracted or scored. Watson's original development comparisons were read as
prior art. The final audit must retain these disclosures.

## Frozen model and inputs

Use Watson's published truncated stretched-exponential DFE,

\[
w_g(s)=Z_g^{-1}\exp[-(s/d_g)^{b_g}],\quad
Z_g=\int_0^{s_{\max,g}}\exp[-(s/d_g)^{b_g}]\,ds.
\]

Every numerical fitness value is copied unchanged in
`data/model_parameters.csv` from author code pinned to commit
`3a250fd3459660812e0ed4392f4e3bf98ea1ce7d`. The separately recovered genes
are DNMT3A, TET2, ASXL1 and TP53. The source also supplies a pooled ten-gene
DFE. Its genes are DNMT3A, TET2, ASXL1, JAK2, TP53, CBL, SF3B1, SRSF2,
IDH2 and KRAS. The pooled DFE is its own published fit, not a sum of the four
gene fits. Do not invent distributions for the other six genes or average
PR #68's distinct mutation classes.

Set the effective population-times-turnover parameter to the published code
value **A=N*tau=94016.661 years**. Fix the per-cell annual nonsynonymous-SNV
mutation supplies in `data/mutation_supply.csv`; the pooled value is
**0.000038011236665 per year**. No target prevalence calibrates these inputs.
Their genomic panels and mutation classes remain attached to the rows.

For a diploid, heterozygous, copy-neutral variant, use Watson's finite-fraction
coordinate `x=v/(1-2v)` and

\[
\phi(s,t)=\frac{\exp(st)-1}{2As},\qquad
\Lambda_g(t,v)=2A\mu_g\int_{0.0001}^{s_{\max,g}}
  w_g(s)E_1\{x(v)/\phi(s,t)\}\,ds,\quad
p_g(t,v)=1-\exp[-\Lambda_g(t,v)].
\]

Here `E1(z)=integral_z^infinity exp(-u)/u du`, and `phi(0,t)=t/(2A)`.
Independent mutation arrivals and independent clones are assumptions. Fixed
fitness, constant mutation supply from birth, no pre-existing prenatal clone,
no competition, no mortality selection, and a common unbiased VAF mapping
are transport assumptions, not findings. The finite-fraction density is a
Watson model convention, not an exact finite-population branching theorem.

The source integrates fitness from 0.0001 while normalizing over 0..s_max;
retain that convention. Use accurate integration of the stated normalization
instead of copying the notebook's numerical rectangle sum. Integrate VAF up
to its mathematical limit 0.5; separately check the notebook cutoff 0.4999.
These numerical choices are fixed here before scoring. Derive the E1
reduction and, separately, the exact integer birth-death clone threshold law
as numbered theorems with proofs; do not conflate cross-sectional detection
with having ever crossed a threshold.

## Detection and observation rules

The reference prediction is **true VAF >= 0.02** under perfect detection.
Use a cohort's explicitly printed VAF cutoff for its source-threshold
comparison; if absent, retain 0.02 as a labelled reference scenario only.
Predeclare sensitivity cutoffs **0.01, 0.03, 0.05 and 0.10** and population
scenarios **A=70000 and A=130000**. All are reported regardless of agreement;
none is selected as a fitted effective threshold or confidence limit.

Record each cohort's printed sequencing depth, depth statistic (mean/median),
variant classes/panel, read support, caller rules, and accessible sensitivity
information. A nominal 2% definition or mean depth alone does not establish
perfect detection at 2%. If depth D and minimum alternate-read count r are
printed, compute the illustrative binomial probability P[Bin(D,v)>=r] at the
reference thresholds, rounding a noninteger reported mean depth to the nearest
integer only for this explicit illustration. It is not a calibrated caller.
Do not infer a depth distribution, false-negative rate, hidden denominator or
unprinted age-by-gene count. Missing assay calibration is a primary-test gate.

The target must match the model's nonsynonymous-SNV panel to receive a strict
mechanistic score. A canonical-driver CHIP endpoint including indels or other
genes is a **descriptive endpoint substitution**, even if the same gene names
appear. Neither endpoint is assumed to contain the other. Preserve genotype
and caller mismatches even if the apparent prediction is good.

## Extraction and cohort independence

After committing, inspect Jaiswal 2014 and Genovese 2014 for context; Bick 2020
TOPMed, later UK Biobank, All of Us, Lifelines and other located publications
for external targets; Fabre 2022 and Mitchell 2022 for model/measurement
counterevidence. Search newer primary publications through 2026-11-18.
Include every eligible printed age row in an encountered eligible table,
without filtering on fit, direction, significance or gene ranking.

Transcribe only small **published population aggregates** in article text or
openly readable tables: age bounds, gene, prevalence/count, printed uncertainty,
and denominator where printed. No plot digitization, individual points,
participant records, raw genomes, controlled data or browser file downloads.
If a source offers only a download, stop there, retain the link and limitation,
and proceed with readable material. Do not reconstruct data from a figure's
axis or estimate bar heights. A missing value remains missing.

Every CSV row must carry source URL, precise text/table/figure or notebook
locator, access date, and a defensible reuse basis. Attribute isolated factual
numeric extracts without reproducing source prose, figure artwork or entire
source tables; retain explicit source licenses where supplied. Public access
alone is not a redistribution license. When terms do not support the intended
extract, retain bibliographic metadata and exclusion only. Repository MIT does
not relicense third-party material.

Record development membership, shared parent cohorts and known/uncertain
subject overlap. Jaiswal/Genovese and Watson development studies are never
independent. Publication after Watson does not by itself prove independence.
Collapse repeated UKB releases and other shared cohorts into a single parent
for holdout/influence analysis. No overlapping subgroup is an extra cohort.

## Scoring, age uncertainty and verdict

For a finite band [a,b], report both `p(a)..p(b)` (distribution-free age
bounds under monotonicity) and an explicitly assumed uniform-age mean for
descriptive point scoring. Preserve half-open boundaries in the source data.
An open-ended band has no point score without a published finite upper bound;
its conservative upper prediction is 1. A reported age mean is not substituted
for an unreported age distribution. Overall prevalence without age support is
context only. Birth is age zero; never choose an onset offset to improve fit.

Observed uncertainty priority: printed numerical confidence limits; otherwise
Wilson 95% intervals when cases and N are both printed; otherwise the rounding
interval at the printed precision. Label rounding and sampling uncertainty
separately. A rounded percentage does not reveal its binomial denominator.
All model parameter sensitivities are scenarios, not 95% prediction intervals.
Unknown covariance/selection and the fitted DFE's missing joint uncertainty
preclude an unconditional significance claim.

Report signed errors in percentage points, MAE/RMSE on finite-band uniform-age
predictions, and the minimum gap between observed and predicted age intervals.
For printed counts also report binomial deviance against the fixed prediction.
Use equal weight within each parent cohort and equal weight across parents;
report source/age/gene-specific rows and the largest deviations. Development,
external-but-overlap-uncertain, strictly eligible, and substituted endpoints
must have separate summaries. Descriptive scores never upgrade eligibility.

The strict verdict is **fitness predicts age-prevalence** only when at least
two independent parent cohorts have matched panel/caller aggregates, every
row's observed uncertainty interval overlaps the frozen age prediction range,
and there are usable gene-specific rows; this is compatibility at stated
resolution, not proof of biological truth. A nonoverlap is **fails** for the
specified compound frozen prediction, not evidence that fitness is irrelevant.
If these measurement/independence/data gates fail, the strict verdict is
**untestable**. Still calculate all permissible descriptive scores and state
whether their frozen interval-compatibility criterion fails. Do not invent a
gene deviation ranking from pooled prevalence when gene-age rows are absent.

## Leave-one-cohort-out and identifiability checks

For the unchanged model, deleting any parent does not change a prediction.
Report recomputed aggregate errors after each parent deletion: this measures
influence and is not fitted cross-validation.

As a separately labelled secondary transport exercise, permit only one common
multiplier c on mutation supply. On training finite-band rows with 0<p_obs<1,
estimate log(c) by the mean of
`log(-log(1-p_obs)) - log(Lambda_at_band_midpoint)`, giving each parent equal
weight. Exclude a held-out parent and every linked release completely. Predict
its row at the midpoint with `1-exp(-c*Lambda)` and compare against a
training-only constant-prevalence baseline. Exclude endpoint probabilities
0/1 from this transform, with a reason; no pseudocounts. Do not refit fitness,
thresholds, ages, N*tau, gene-specific supply, or a slope. Run only if at least
two distinct external parent cohorts provide finite-band rows. Any unresolved
overlap retains a non-independent label. Keep all-substitute and strictly
matched versions distinct. This is not the unchanged-parameter primary test.

Prove the sharp common-c feasibility interval implied by each age bin and its
observed bounds, then intersect across bins. Prove the mutation-supply versus
constant assay-thinning alias and explain what additional published aggregate
would break it. Do not claim that unrestricted assay confounding makes every
possible age-shape test impossible. State a quantitative falsifiable prediction
for matched future aggregate reporting, including genes, ages and threshold.

## Verification, review and delivery

Use only standard-library Python in analysis/verifier and the single new test.
Recompute generated numerical results from the small CSVs. Check independent
E1/quadrature calculations, birth-death generating equations and boundary
cases, age-bin sharpness, nuisance aliases, no held-out leakage and all
exclusions. Normalize CRLF and CR to LF before SHA-256; every manifest path
key is constructed with `Path.as_posix()`. The dedicated test must run in
under 60 seconds on ordinary CI hardware and work on Windows and Linux.

Keep changes only in `research/ch-age-prevalence-test/` and
`tests/test_research_ch_age_prevalence_test.py`. Do not modify the research
index, other notes or workflows. Hostilely review math, sources, scope and
claims; downgrade anything unproved. Push the branch, open a PR, and inspect
the existing CI. Neither passed tests nor AI review is external scientific or
clinical validation.
