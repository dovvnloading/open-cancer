# Frozen protocol: viral-antigen offset transport

**Nonclinical, unreviewed computational research. No treatment, diagnosis,
surveillance or patient-selection advice.**

Date: 2026-11-18, America/New_York. Branch from current main
`e1c48e377bd438ac13ea94dc53d1a9e7db088de3`.
This file must be the sole addition in the first branch commit and remain
unchanged. Commit it before searching for, opening or transcribing any new
HPV, EBV or other target-cohort response outcomes, and before scoring them.
Only the requested existing repository notes and their already recorded
inputs were read before this freeze. This timing does not erase prior
literature knowledge, establish an independent registry or make published
data a blinded prospective experiment. Record the commit and any subsequent
decisions, corrections or deviations in README.md; never rewrite this plan.

Read first: [research programme](../README.md),
[PR72](../merkel-viral-antigen-test/),
[PR67](../immune-escape-pancancer-test/),
[PR59](../immune-escape-neoantigen-clonality/),
[DATA_POLICY.md](../../DATA_POLICY.md) and [AGENTS.md](../../AGENTS.md).

## Question, estimand and scope of transport

Does the exact response-equivalent viral-load offset calibrated on JAVELIN
Merkel 200 Part A improve the frozen reciprocal-TMB curve's predictions in
other virus-associated cancers? The observable estimand is a report's
published viral-status-specific objective response rate (ORR) to anti-PD-1
or anti-PD-L1 monotherapy, in its stated treatment population and at its
stated cutoff. RECIST response is not eradication. Viral status is not
randomized and is not a measurement of presented functional antigen load.

For a complete cohort pair the predictor is its published median TMB in each
viral stratum, preferably in exactly the ORR population. A sequencing subset
in the same trial can supply a separately flagged ecological pairing; a
different untreated cohort cannot supply the primary predictor. Different
ORR and TMB denominators must remain separate. Missing values stay missing.
The plug-in assumption `E[f(T)] = f(median(T))` is generally false and is an
explicit limitation, not an identity established by PR59.

"Out of sample" means no new target outcome can alter the inherited curve,
offset, scale grid or decision rules. It does not automatically mean new
patients relative to every historical PR67 source. Audit the original
Mo/Yarchoan response lineage. Mark known or possible overlapping historical
trials, repeated follow-ups, pooled datasets and biomarker subsets. Report a
strict independence sensitivity excluding known/possible training overlap;
if independence cannot be established, do not claim independent validation.
The esophagogastric PR67 row had missing TMB and was excluded from its fit;
the head/neck parent was included. Do not refit to remove either parent.

## Exact inherited model and offsets: no fitting

Use the full-precision snapshot `reciprocal_holdout` in PR72
`data/frozen_curves.csv`, representing PR67's original leave-Merkel-parent-out
curve. The constants are

```
A = 0.48718323459062135
B = 5.5353933555249
p0(m; c)    = A * exp[-B / (c*m)]
p1(m,z; c)  = A * exp[-B / (c*m + v_c*z)]
```

The continuous value at zero effective load is zero; TMB must be nonnegative.
`z=1` for the source-defined viral-positive category and `z=0` for negative.
Unknown, equivocal or unassessed status is not negative. Preserve p16, direct
HPV testing, EBER, PCR, serology and viral-load definitions distinctly.
No new coefficient, intercept, ceiling, status correction or conversion is
estimated from target outcomes. Do not select a winning curve or scale.

Carry the following exact broad-population Part A rows of PR72's
`data/part_a_calibration.csv` forward unchanged:

| c | Point v_c | Part-A CI-low v_c | Part-A CI-high v_c |
| --- | --- | --- | --- |
| 0.5 | 9.945464925200476 | 4.726315626607673 | 48.61354519512048 |
| 1 | 9.700464925200475 | 4.481315626607673 | 48.368545195120475 |
| 4/3 | 9.53713159186714 | 4.31798229327434 | 48.20521186178714 |
| 2 | 9.210464925200474 | 3.991315626607673 | 47.87854519512048 |
| 5 | 7.740464925200475 | 2.521315626607673 | 46.408545195120475 |

Primary scenario: `c=1`, point `v=9.700464925200475`. All five point scenarios
are required sensitivities. The two inherited Part-A CI-endpoint offsets at
each scale are training-uncertainty scenarios, not confidence or predictive
intervals for a transported parameter. Do not fit the new outcomes, calibrate
a new viral species, re-use the Part B outcome, or promote a favorable
endpoint sensitivity to the primary model. The Part-A TMB-subset calibration
is not used: its different estimand is not another selectable offset.

For complete pairs also show the inherited PR72 `log_holdout` TMB-only
comparator, `clip(a+b*ln(c*m),0,1)`, with
`a=-0.03576193945209477`, `b=0.12639822421480662` and zero-input value zero.
It is context for PR67's negative baseline comparison, not a newly fitted
viral extension and not part of the primary pass/fail comparison.

## Assay-scale assumptions

The offset is in response-equivalent Chalmers-scale mutations/Mb, not a count
of viral antigens. PR72 mapped NSSV/Mb to that scale with `x=c*m`. Treat
published per-Mb target TMB medians analogously only as an explicit numerical
bridge: count definition, panel territory, germline/driver filters and
synonymous inclusion differ. The same fixed c applies to both statuses and
every cohort within a reported scenario. None of the five c values is an
empirically calibrated conversion for these new assays. A validated external
assay mapping, if found, is described separately and cannot replace the grid.

Do not convert mutations/exome, mutations/tumour, mutation counts, means,
neoantigen counts or TMB thresholds to median mutations/Mb using an assumed
exome size. Retain such published aggregates as context with their actual
units; they are ineligible for plug-in scoring. An exact per-Mb equivalent
printed by the source may be used and must have its own locator. No medians
are inferred from boxplots, patient points, thresholded categories or pooled
medians. No status-specific assay multiplier is allowed.

## Target search, eligibility, duplication and stopping

Required anchors, including unfavorable or inaccessible results:

1. HPV-positive versus HPV-negative recurrent/metastatic head and neck
   squamous cell carcinoma under anti-PD-1: KEYNOTE-012, KEYNOTE-040,
   KEYNOTE-048, CheckMate 141, and later pooled/biomarker analyses including
   KEYNOTE-012/055 when found. Keep p16-defined oropharyngeal positivity
   distinct from direct HPV detection and non-oropharyngeal classification.
2. EBV-positive versus EBV-negative advanced gastric/GEJ cancer: Kim et al.
   2018 Nature Medicine and later anti-PD-(L)1 monotherapy cohorts. A
   non-EBV molecular subgroup such as MSS/EBV-negative is not automatically
   the entire EBV-negative cohort; MSI and other selections remain explicit.
3. Targeted searches for other virus-associated tumours with printed
   viral-status-stratified ORR: hepatocellular carcinoma (HBV/HCV categories),
   nasopharyngeal, cervical and anal carcinomas. Infection etiology or plasma
   viral DNA is not silently equated with tumour viral-antigen presentation.
   These are secondary transport challenges with their assay/biology limits.

Search primary articles and openly readable main/supplemental tables and
article text. Reviews may locate primary reports but do not supply the
primary numeric rows. Stop after all named anchors, their relevant accessible
updates and targeted missing-aggregate/contrary-result searches have an
entry in the source ledger. This is a bounded targeted audit, not an
exhaustive systematic review. No arbitrary minimum response selects a report.

Include published aggregate ORR points even when intervals, denominators or
TMB are absent; label each unavailable computation. Exclude survival-only,
disease-control-only, neoadjuvant pathological-response and combination arms
from the ORR transport score. Keep treatment lines, PD-L1 selection and
response assessment method. A trial's monotherapy arm remains eligible when
other arms use chemotherapy, but their ORRs are not mixed.

For overlapping reports rank availability first: both viral-stratum ORRs,
published 95% intervals, denominators and usable within-cohort TMB medians;
then both ORRs with intervals/denominators; then a single viral stratum with
interval/denominator; then less complete ORR reports. Within a rank choose
the latest publication date, then lexical DOI order to break an exact tie.
This hierarchy uses availability, never outcome magnitude. One selected
report per overlap family is primary; retain other versions as context.
A positive-only primary report may support the range check but not paired
scoring. If a pooled report supplies otherwise unavailable TMB and ORR,
score it as one overlap family and do not also count its components. Never
pool status medians, reconstruct an unreported stratum by subtracting
rounded ORRs or treat repeated trial follow-ups as independent replications.

## Transcription, provenance, access and reuse

Only printed factual aggregate numbers in article text or openly readable
tables are transcribed. Do not digitize plots or transcribe individual
patient rows. A numerical aggregate label accessible as text may be used
with an exact figure/table locator and inspection limitation. No inferred
responder counts or denominators from rounded percentages. Preserve printed
precision, interval level/method, source population and missingness reasons.

Every row, including missingness/source-ledger rows, must include source
title, DOI or identifier, URL, precise text/table/figure locator, access date
in America/New_York, reuse basis and terms URL. Record source license terms
per paper; repository MIT does not relicense third-party facts. Prefer
explicitly reusable open-access sources. If reuse rights cannot be
established, retain metadata and the limitation without copying the numbers.
No original articles, figures, source workbooks or patient records are
committed. Distinguish direct article/table reading from search-indexed text.

**No browser file downloads.** If a download prompt appears, proceed without
it and record the limitation. Use readable HTML/text or a readable table
instead; do not work around a download prompt or access wall. No source-file
downloads are needed to build the aggregate CSVs. No restricted/controlled
access data, authentication circumvention, patient-level reconstruction,
private credentials or clinical advice.

## Scoring and honest uncertainty

For each complete pair, compute both status predictions and residuals
`observed minus predicted`. The primary score is the equally status-weighted
MSE on the probability scale. Report RMSE and MAE in percentage points,
interval inclusion, number of interval misses and mean squared distance to
the reported interval (zero inside). Keep denominator-matched and
within-cohort ecological/subset pairings separately identified. Do not use a
binomial likelihood, patient-level Brier score or bootstrap of pseudo-patients.

For summaries, average within an overlap family first, then equally over
independent cohort families. If complete families occur in both required
cancer groups, also report the equally cancer-group-weighted result. Never
call an empty panel's score zero. Point-only pairs can have a point score but
cannot pass the full interval rule. No unreported CI is invented. Published
95% ORR intervals are required for a complete interval-rule verdict.

For each available ORR interval box calculate the exact range of
`MSE_offset-MSE_baseline` as the true ORRs range over the box. At fixed
predictions each squared-error difference is affine in the unknown ORR, so
its extrema occur at interval endpoints. This is a descriptive uncertainty
box, not a joint confidence interval or multiplicity-controlled significance
test. Positive/negative direction at a point alone is not statistically
demonstrated improvement. Published counts and denominators are used only
to check consistency; unknown counts stay unknown.

Required stress tests: every fixed scale, both inherited CI-endpoint
offsets, exclusion of possible historical overlap, retention of the
unchanged negative control, and source rounding. Expand printed ORR points
and interval endpoints by half their last printed unit and per-Mb medians
likewise (truncate at zero and probabilities at one); use monotonic bounds
for predictions and disclose any changed inclusion or score direction.
Unknown TMB sampling uncertainty and parameter-estimation uncertainty cannot
be recovered from these summaries.

Also show the conservative population-median envelope. For a monotone curve
with lower limit l and upper limit A, a nonnegative population distribution
with median m has mean prediction in the closure
`[(l+f(m))/2, (f(m)+A)/2]`. For the baseline/negative stratum l=0;
for the offset-positive stratum l=`A*exp(-B/v_c)`. Report whether that envelope
intersects the ORR interval. This relaxes plug-in substitution but still
assumes the sample median describes the response population; it supplies no
fitted distribution, calibrated confidence band or biological explanation.

### Checks possible without a usable TMB median

For finite TMB, both models have ceiling A, approached but not reached.
The baseline ranges from zero to A. With positive v the viral-positive
extension has lower endpoint `A*exp(-B/v)` at zero TMB and the same upper
limit A. Consequently an ORR interval wholly below this floor or with lower
endpoint at/above A is incompatible with the specified model for every
nonnegative TMB and every positive multiplicative assay scale. Check these
range restrictions at every frozen v, including inherited endpoint scenarios.
Use closures for conservative rounding/uncertainty compatibility; exact
boundary-only cases must be disclosed, not called a robust contradiction.
An above-ceiling contradiction challenges the common ceiling/response
bridge of both models; it is not evidence specifically against viral antigens.

This necessary-condition check may falsify universal numerical transport
without identifying a paired MSE. It must never be presented as a prediction
score at an invented TMB, or as proof that the baseline beats the offset.

## Frozen verdict rule

Per complete cohort pair, a **descriptive transport pass** at the primary
scenario requires lower offset MSE and both offset predictions inside their
published 95% ORR intervals. A missing component is untestable, not a pass.
Report all sensitivities; favorable scales do not replace c=1.

Give one explicit top-level verdict with this priority order:

* **Viral offset fails to transport** under the declared numerical rule if
  any eligible required-cancer primary cohort fails the complete-pair rule,
  or an eligible required-cancer viral stratum has a rounding-robust
  TMB-independent range contradiction. Only selected primary reports at
  c=1 with the broad-point offset can trigger this top-level failure;
  other scales, CI-endpoint offsets and overlapping context rows are
  sensitivities. Attribute an above-ceiling failure
  to the shared response bridge, not uniquely to the offset. Label possible
  historical overlap and show the strict independence sensitivity.
* **Viral offset transports**, conditionally and descriptively, only if at
  least one nonoverlapping complete pair in each of HPV/HNSCC and EBV/gastric
  passes, every eligible complete primary pair passes, and no required-cancer
  range contradiction remains. Documented target-outcome independence from
  the inherited fit is required for this overall positive verdict; numerical
  passes with unresolved overlap remain conditional descriptive results.
* **Untestable** otherwise: insufficient complete pairs/intervals or
  unresolved independence, without a decisive declared-rule contradiction.

Always report separately whether comparative improvement over TMB alone is
testable and supported. An absolute range failure can coexist with an
untestable comparative score. Secondary tumour categories cannot manufacture
the required HPV/EBV coverage. These are frozen descriptive falsification
rules, not calibrated hypothesis tests, causal findings or clinical rules.

## Falsification, novelty, reproducibility and delivery

List every attempted falsification and unavailable test: fixed-size transfer
across viral species; baseline/offset errors; unchanged negative control;
ceiling and offset floor; alternative inherited scales and uncertainty;
nonlinear median substitution; assay and status definitions; different lines
and selected populations; missing medians/intervals; historical overlap and
pooled-report duplication; and counterexamples in later cohorts.

Explain the distinction between a conceptual
`k_total=k_mutational+k_nonmutational` contribution to PR59's functional
per-cell load and a universal equal-sized response-equivalent offset. The
theorem assumes equal quality, presentation, source suppression and a
large-minimum-load limit; TMB and viral presence do not measure those inputs.
Audit primary prior work on viral status, TMB and immunotherapy response.
Novelty is at most this repository's frozen-parameter cross-virus
transport/identifiability audit, not discovery of viral immunogenicity,
a validated biomarker, or a bibliographic priority claim.

Create only `research/viral-antigen-transport-test/` (this protocol, README,
small aggregate/provenance CSVs in data/ and analysis code) and
`tests/test_research_viral_antigen_transport_test.py`. Do not modify any
existing file, including research/README.md or other research directories.
Use Python stdlib only, offline, and keep the dedicated test under 60 seconds.
Recompute every reported analytical number from CSVs; generate and verify
README numeric tables, missingness, overlap grouping and all verdicts.
Read inherited CSVs without editing them. Pin their newline-normalized
SHA-256 hashes and the protocol hash. Normalize CRLF/CR to LF before hashing;
build every file/path key with `Path.as_posix()` and sort canonical keys.
Test newline/Windows-path portability and absence of any target-driven fit.

Before publishing results, independently recheck every transcribed numeric
field against its source, hostilely review the scientific claims and downgrade
anything unproven. Run the dedicated test and applicable repository gates.
Push the branch, open a PR, observe its actual CI and report PR number and
status. CI, hashes and AI critique establish neither source truth nor
external scientific or clinical validation.
