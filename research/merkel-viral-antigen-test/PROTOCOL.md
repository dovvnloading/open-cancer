# Frozen protocol: Merkel viral-antigen and assay-scale test

**Nonclinical, unreviewed computational research. No treatment, diagnosis,
surveillance or patient-selection advice.**

Protocol date: 2026-11-18, America/New_York. Branch base is current main
`16c7596637a60836056593c1fb22b1f5c9d07501`. This file is the sole addition in
the first branch commit and must remain unchanged thereafter. No model fit,
curve evaluation, assay compatibility calculation or outcome score for this
task precedes that commit. Reading prior results, source methods and published
aggregate values to establish eligibility preceded the freeze. Consequently
this is a pre-computation plan for a retrospective analysis, **not blinded
preregistration, an independent registry or an unseen-outcome experiment**.

Read first: [ranked job and negative results](../README.md),
[#67](../immune-escape-pancancer-test/),
[#59](../immune-escape-neoantigen-clonality/), [AGENTS.md](../../AGENTS.md),
and [DATA_POLICY.md](../../DATA_POLICY.md). All later decisions, extraction
corrections and departures belong in the new README, with their timing.

## Question and estimand

Does the frozen TMB-only response surrogate that produced #67's largest
Merkel-cell residual remain descriptively compatible with published viral
strata after a **common** assay conversion? Can one nonnegative viral-load
parameter fixed on later-line Part A improve prediction of distinct first-line
Part B without altering the virus-negative control?

The observational target is each report's published aggregate objective
response rate (ORR, RECIST independent review), conditional on its reported
MCPyV classification and treatment population. Predictor is the published
median NSSV/Mb in its TMB-evaluable viral subgroup. These are separate marginal
summaries, often with different denominators and cutoffs. Substitution of a
median into a nonlinear response curve is an explicit observation assumption:
`E[f(T)] != f(median(T))` in general. There is no individual-data fit,
causal viral effect, functional-antigen measurement or clinical validation.

Part A is the previously chemotherapy-treated cohort in
[D'Angelo et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7239697/),
DOI `10.1136/jitc-2020-000674`; Part B is the previously untreated metastatic
cohort in [D'Angelo et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311489/),
DOI `10.1136/jitc-2021-002646`. They are distinct patient cohorts in the same
trial, with shared drug, infrastructure and selection conventions. They are
not independent programmes, randomized viral strata or interchangeable
treatment lines. Earlier subsets and later follow-ups of either cohort are
not extra replications. Independence from all historical response inputs in
#67 is not presumed; primary removal of the Merkel parent avoids direct
reuse of that cancer's response in the frozen training curve.

## Frozen curves: identity before evaluation

Use the exact `analysis.py` and `data/tmb_orr.csv` at the base commit above.
Read nonmissing positive TMB rows from that CSV and the original ORRs; do not
import clonal or HLA features. Recreate the original deterministic `fit`
algorithm, without any JAVELIN data in its training inputs:

* **Primary:** `reciprocal_tmb` after excluding every `parent_type=merkel`
  row, the original fold used to calculate the largest Merkel residual.
  Its form is `f(x)=A*exp(-B/x)` with the original constrained least squares.
* **Required sensitivities:** the all-panel reciprocal fit and the
  `log_tmb` fits on these same two frozen training sets. The latter use
  `clip(a+b*ln(x),0,1)` as in #67. Do not change the fitting methods or select
  a curve after examining JAVELIN compatibility.

The base code and input CSV, not rounded printed parameters, fix the curves.
Report their full-precision regenerated coefficients and newline-normalized
SHA-256 fingerprints. No refitting of A, B, a or b on either JAVELIN cohort.
Parameter-estimation uncertainty is not added or optimized: this question is
about the fixed curves, not whether a different pan-cancer curve could fit.

## Admissible assay bridges, fixed without outcome calibration

Let `m` be JAVELIN NSSV/Mb and `x=c*m` be the input on the Chalmers scale.
The admissible algebraic family is **one finite positive multiplicative c,
shared across both viral strata and both cohorts**, with no intercept,
status-specific multiplier, cohort-specific multiplier or nonlinear warp.

[Chalmers methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC5395719/)
count coding substitutions and indels, including synonymous changes, with
driver and germline filtering; JAVELIN reports nonsynonymous somatic variants
from matched tumour/blood whole-exome profiles. Differences in territory,
filtering and variant ascertainment prevent a universal synonymous-count
conversion. In particular, `c>=1` is not guaranteed across these assays.
No published paired calibration for the exact assays is established here.

Report exact compatibility over **all c>0** and its intersections with these
predeclared, explicitly assumed sensitivity bands:

* `[1,2]`: a counting-convention scenario, not a calibrated interval;
* `[0.5,5]`: a deliberately wider assay scenario, also not calibrated.

Use the fixed grid `c in {0.5,1,4/3,2,5}` for predictions and scores. The
headline reference is `c=1`, merely numerical identity of the reported units;
`4/3` is illustrative, not an estimated synonymous fraction. Report every
grid point. Never choose a scale from Part B to improve the score, interpret
an unbounded compatible scale as an assay calibration, or treat the bounded
bands as empirical confidence intervals. All-positive feasibility is a
falsification check, not evidence of a plausible assay mechanism.

## Extraction, sources and denominators

Transcribe only published cohort-level factual summaries into small CSVs.
Each row must contain source title/DOI/URL, exact figure/table/text locator,
access date, reuse basis and terms URL, population, endpoint, units and notes.
Use published numeric labels/text before reading plot geometry. No patient
points, reconstructed records, individual sequencing data or inferred
unreported responder counts/denominators. No source PDFs are committed.
Missing values remain blank, with reasons; missing does not mean zero.

Primary rows are MCPyV+ and MCPyV- broad ORRs and their published 95% intervals,
paired only ecologically with the corresponding TMB medians. Keep total and
viral-stratum ORR denominators separate from total and viral-stratum TMB
denominators. Preserve source definitions of evaluability: an all-treated ORR
denominator may include people without an evaluable best response. Record
unclassified viral groups where published, and retain cutoffs for both endpoints.

Part A response/biomarker Results and Figs. 1/3, and Part B Fig. 1D,
biomarker Results and Supplementary Fig. 3B are required anchors. Part A
Fig. 4 provides a prespecified **secondary** ORR sensitivity restricted to its
TMB subset; it cannot silently replace the broad primary ORR. If corresponding
Part B viral-specific TMB-subset ORRs are not published, the denominator-matched
held-out test is unavailable. Do not manufacture it from response-group and
virus-group margins. Preserve Part A one-versus-multiple-prior-treatment ORRs
as context; these are not virus-by-treatment-line cells. No pooling of Parts A/B.

Verify reuse notices individually. Repository MIT does not relicense source
observations. Record any source accessible only as indexed text and distinguish
it from direct visual verification. Unknown rights mean metadata only.

## Compatibility rule and sensitivity to median substitution

For each frozen curve and cohort/status pair, calculate the entire set
`C={c>0: L <= f(c*m) <= U}`, using published ORR interval endpoints `[L,U]`
on the probability scale. Invert the monotone curve analytically, including
zero, clipping and ceiling boundaries. A reciprocal ceiling A is approached
but never reached at finite positive input when B>0. Infinite-only solutions
are empty, not compatible. Check inverse endpoints by forward evaluation.

Intersect the two statuses within each cohort, then intersect Parts A and B.
Report individual sets, cohort sets, joint set and both assumed-band
intersections. Do this separately for each frozen curve, without optimization
or model selection. These are deterministic **compatibility regions**, not
joint confidence regions, equivalence tests, calibrated hypothesis tests or
p-values. Broad published marginal intervals do not supply a joint likelihood.

Primary values are printed medians and endpoints. Prespecified rounding
sensitivity expands each published number by half its last printed unit:
one-decimal percentages by 0.05 percentage points, two-decimal medians by
0.005 NSSV/Mb (other precision follows its printed digits). Clip probability
limits to [0,1]. Report the union of scales compatible with some rounded
median/interval realization; do not call rounding error sampling uncertainty.

Also relax median substitution for the **primary reciprocal curve**. For any
nonnegative TMB distribution with a population median m and bounded monotone
f ranging from 0 to A, its mean predicted response lies in the conservative
closure `[f(c*m)/2, (A+f(c*m))/2]`. Compute compatibility if this envelope
intersects `[L,U]`. This is a distributional stress test, not a recovered
distribution or empirical confidence interval. It still assumes the measured
subset's median represents the ORR population; missingness and sample-median
uncertainty are not identified. Do not apply this formula to the added viral
model without accounting for its changed lower endpoint.

## Exactly one added viral parameter and held-out score

For the primary reciprocal curve only, define
`p(m,z;c,v)=A*exp[-B/(c*m+v*z)]`, `v>=0`, where z=1 for MCPyV+ and zero for
MCPyV-. The added v is a response-equivalent load offset in Chalmers-scale
units, **not** a measured number of viral antigens. It does not follow as a
theorem from #59; equal presentation, quality and common response bridge are
additional assumptions. Do not add a cohort intercept, new ceiling, treatment
effect or viral-specific scale.

At each fixed grid c, fix v using **only Part A broad MCPyV+** point ORR yA
and median mA. If `0<yA<A`, set
`v=max(0, B/ln(A/yA)-c*mA)`. This minimizes squared calibration error for
that stratum; report any zero boundary and the remaining residual. If
`yA>=A`, no finite exact match exists: report the unattained boundary and
the extension as unavailable instead of choosing an arbitrary large v.
For `yA<=f(c*mA)`, use v=0. The virus-negative prediction is unchanged.
The five fits are separate prespecified scale sensitivities with one added
parameter each, **not** a fit of both c and v or a search for a winning c.

After protocol freeze and Part A transcription, compute and commit a small
Part-A calibration CSV with v and frozen coefficients for every grid point
**before running any Part B prediction or score**. Include the Part A
TMB-subset calibration as separately labelled sensitivity at the same grid;
it will be scored against broad Part B only as an expressly mismatched
stress test if matched Part B ORRs are absent. Preserve this calibration
commit. Automated recomputation must establish that altering Part B cannot
change any fitted v. Reading published Part B aggregates before this stage
does not make subsequent scoring blinded; disclose it.

On Part B report, for both statuses and every c, baseline/extension predicted
ORRs, residuals (observed minus predicted), interval inclusion and distance
to the nearest interval endpoint (zero inside). Primary descriptive score is
`MSE=((p_plus-y_plus)^2+(p_minus-y_minus)^2)/2`; also report RMSE, MAE,
number of interval misses, and mean squared interval distance. Equal status
weight avoids conflating different ORR and TMB denominator precision. No
binomial likelihood, reconstructed response count, patient-level Brier score,
bootstrap significance or p-value is justified by the ecological pairing.

A numerical repair must improve status-balanced Part B MSE **and** put both
predictions inside their published intervals at the stated fixed c; report
failure of either condition. An improvement confined to the positive stratum
cannot repair a failing negative control. Report the full fixed-grid result
and the reference result even when adverse. Fit the Part A interval endpoints
as a separate sensitivity for v when finite; these envelopes are conditional
training uncertainty scenarios, not calibrated predictive intervals.

## Verdict, falsification and novelty rules

Give an explicit verdict using the requested categories: explained by viral
antigen supply / explained by assay scale / unexplained / untestable. Attach
its estimand and assumptions. A compatible arbitrary c or improved v does
not establish a causal explanation. An assay-scale explanation requires
independent evidence that the compatible scale is admissible for these
assays; a viral explanation requires more than fitting a low-TMB residual.
If the available marginals/uncalibrated bridge cannot distinguish explanations,
the attribution verdict is **untestable**, with numerical surrogate successes
and failures still reported. Use **unexplained** for a failed declared
surrogate repair, without upgrading it to a biological rejection.

List attempted falsifications: common scale within and across cohorts;
bounded versus all-positive scales; leave-Merkel-out versus full-panel curves;
log-TMB comparator; rounding; changed ORR denominator/cutoff; median-envelope
relaxation; the unaltered virus-negative control; transport of a Part-A-only
viral offset; source/assay classification uncertainty; and historical overlap.
Keep failed, null and unavailable checks visible. No silent dropping of a
cohort, endpoint or adverse scale. Corrections of transcription mistakes are
allowed only with an explicit source-audit record, never outcome-driven editing.

For #59 distinguish its functional per-cell antigen load from mutation-derived
TMB. A conceptual `k_total=k_mutational+k_nonmutational` term could include
presented viral antigens, but the equal-quality, common-presentation and
large-load assumptions still need justification. This exercise neither
estimates harmonic functional exposure nor proves/refutes #59's theorem.
Credit existing MCPyV biology, trial biomarker analyses and ecological TMB
response regressions. Novelty is at most a repository-specific frozen-curve
transport/provenance audit; no new mechanism, biomarker or priority claim.

## Reproducibility and scope

Create only this directory (README, this protocol, small CSVs in data/ and
analysis code) and `tests/test_research_merkel_viral_antigen_test.py`.
No existing file changes, including the research index or other research
directories. Use Python standard library only, offline, with a bounded test
workload below 60 seconds. Recompute every reported analysis number from the
published CSVs, including the unchanged #67 input CSV; verify generated README
numeric blocks, source schema, missingness, denominator separation, independent
inverse/fit checks, Part B leakage and rounding/median diagnostics.

Hash UTF-8 text after CRLF/CR-to-LF normalization. Every manifest/path key must
be constructed with `Path.as_posix()` for Windows/Linux equivalence. Preserve
protocol and calibration fingerprints and commits. Before pushing, independently
recheck each transcribed number against its source, hostilely review scientific
claims and scope, and downgrade unproven assertions. Run the repository gates,
push the branch, open a PR, observe its CI and report actual PR/status. CI and
AI review test code/bookkeeping; they do not confer scientific validation.
