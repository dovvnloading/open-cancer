# Frozen protocol: outcome-matched Merkel TMB bins

**Nonclinical, unreviewed computational research. No treatment, diagnosis,
surveillance or patient-selection advice.**

Protocol date: 2026-11-18, America/New_York. Branch from current main
`37c37e2b4cfc711ce938201123ef3085b8868f2f`. This file alone must be committed
before any new curve evaluation, fit replay, scale bound or outcome score.
Preserve that commit and this file unchanged. Reading existing results,
mathematical derivation and article text to establish eligibility precedes
the freeze; source outcomes are already visible. This is a retrospective
pre-computation specification, **not blinding, independent preregistration
or an unseen-outcome experiment**. Record later corrections or departures,
with timing, in the new README.

Read first: [ranked job 3 and its cautions](../README.md),
[PR #72](../merkel-viral-antigen-test/),
[PR #67](../immune-escape-pancancer-test/),
[AGENTS.md](../../AGENTS.md), and [DATA_POLICY.md](../../DATA_POLICY.md).
The user-authorized scope is this directory and the single new test named
below; no task catalog, result bundle, index or existing file is modified.

## Question, population and additional observation assumption

Does each frozen monotone TMB-only response curve survive the published
outcome-matched low/high TMB bins in JAVELIN Merkel 200 Parts A and B,
without replacing a distribution by its median?

Let T be a finite, nonnegative NSSV/Mb value in the **TMB-evaluable** cohort.
The low bin is `0 <= T < 2`; the high bin is `2 <= T < infinity`.
For each cohort j and bin b, explicitly assume the mean-response bridge

`ORR_(j,b) = E[f(c*T) | cohort j, bin b]`.

Here ORR is the population probability corresponding to the published
aggregate objective-response endpoint. Sampling uncertainty is represented
only by the printed marginal intervals as specified below. **Applying an
ecologically fitted cancer-type curve inside a within-bin expectation is
an additional assumption, not a consequence of its original regression or
PR #59's theorem.** RECIST response is not eradication.

Use one finite positive multiplicative assay scale `c`, common to both bins
and both cohorts. Permit arbitrary different within-bin distributions in
each cohort, finite T almost surely, with no upper support bound, median,
mean, variance, parametric family or inferred patient records. The bin
threshold is the source-defined exact classification rule, not a rounded
measurement to perturb. Monotonicity implies low mean at most `f(2*c)` and
high mean at least `f(2*c)`; exact attainability, not these weak inequalities
alone, determines feasibility.

Part A is the previously chemotherapy-treated metastatic cohort in
[D'Angelo et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7239697/),
DOI `10.1136/jitc-2020-000674`, Results "Exploratory biomarker analyses",
the TMB-evaluability and high/low-TMB ORR paragraphs. Part B is the first-line
metastatic cohort in
[D'Angelo et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311489/),
DOI `10.1136/jitc-2021-002646`, Results "Tumor mutational burden".
Use the reported ORRs for these TMB bins, never broad viral-status ORRs.
The cohorts differ in treatment history and share trial infrastructure;
both were used in PR #72. This is reused evidence, not a new independent
trial replication. Removing the Merkel parent from the primary training
curve does not establish independence of all historical sources.

## Frozen curves and scale domains

Read the four existing full-precision rows in
`research/merkel-viral-antigen-test/data/frozen_curves.csv` at the base commit.
Use those recorded coefficients unchanged, with their original definitions:

* Primary `reciprocal_holdout`: `f(x)=A*exp(-B/x)` for x>0, `f(0)=0`,
  trained with all Merkel-parent rows excluded.
* Required sensitivity `reciprocal_all`: the frozen full-panel reciprocal.
* Required sensitivities `log_holdout` and `log_all`:
  `f(x)=clip(a+b*ln(x),0,1)` for x>0 and `f(0)=0`.

These snapshots have positive A, B and b. Pin the original CSV and relevant
PR #72 / #67 code and protocol by newline-normalized SHA-256. Replay the
historical training calculation only as a numerical identity check; it must
not replace the snapshots or use any new JAVELIN outcomes. Retain the
inherited absolute coefficient comparison tolerance `2e-6`. Any independent
optimizer or gradient diagnostic is checked against a declared error bound,
never exact saved optimizer output. No new curve, viral offset, ceiling,
cohort intercept, Part B outcome parameter or preferred-curve selection.

Compute the full feasible set over **all finite c>0**, then intersect it with
the already-declared assumed bands `[1,2]` and `[0.5,5]`. These are sensitivity
scenarios, not measured assay-calibration intervals. Also retain the inherited
fixed grid `c in {0.5,1,4/3,2,5}`, with identity `c=1` only a numerical
reference. Do not optimize or select a scale using Part B. A feasible set is
a mathematical identified set under the bridge, not an estimated calibration.

## HTML-only aggregate extraction and reuse

Read only article HTML text, including an explicitly identified indexed
HTML rendering if direct article text is blocked. Do not read PDFs, download
files through a browser, retrieve supplementary archives, digitize figures,
reconstruct participant records or circumvent access controls. If a browser
file-download prompt appears, deny it, proceed without that file and record
the limitation. Record access failures and the actual successful text route.

Create one small aggregate CSV with four outcome rows, one per cohort/bin.
Each row records printed ORR percent, printed lower/upper 95% limits, bin
rule/units, the printed total TMB-evaluable cohort denominator, and any
explicitly printed bin denominator. Missing bin counts remain blank with
their reason. Do not infer a denominator or responder count from rounded
percentages, confidence limits, differences or neighbouring viral groups.
The repeated total denominator is cohort context, never the bin denominator
or a weighting factor. Do not sum repeated totals or pool cohorts.

Every row must carry article title, DOI, canonical source URL, inspected HTML
URL/route, Results section, paragraph locator (including short opening anchor
or paragraph description), denominator locator, access date/timezone,
inspection limitation, reuse basis and licence URL. Preserve printed decimal
precision as strings. Only factual numerical summaries and minimal locators
are transcribed; do not mirror source prose or figures. Verify article reuse
notices independently. Part A's CC BY and Part B's noncommercial conditions
remain attached; repository MIT does not relicense third-party observations.
Unknown rights or irreconcilable source/denominator ambiguity means metadata
only and an untestable affected check, not imputation.

## Sharp attainable ranges and scale inversion

Write `s=f(2*c)`. For positive-B reciprocal curves the attainable bin means
are exactly `low=[0,s)` and `high=[s,A)`. The asymptotic ceiling A is never
attained by any distribution on finite T, even with unbounded support. Zero
is attainable in the low bin. For clipped log curves define
`x0=exp(-a/b)` and `x1=exp((1-a)/b)`. The high mean range is `[s,1]`.
The low range is `{0}` if `2*c<=x0`, `[0,s)` if `x0<2*c<=x1`, and `[0,1]`
if `2*c>x1`. In particular, `2*c=x1` does not attain low-bin mean one.

Prove necessity from the response ranges and strict endpoint expectations.
Prove sufficiency by finite-support distributions (a single permitted T
already achieves every interior response for these continuous curves).
These are mathematical feasibility witnesses, never reconstructed patients.
For every curve, each bin and cohort, invert the intersection of its exact
mean range with a target response interval analytically. Track all open and
closed finite endpoints. Zero and infinity in a scale set are excluded domain
limits. Infinite-only, unattained singleton and out-of-response-range
solutions are empty. Intersect bins within cohort, then across cohorts.
Report the per-bin, per-cohort, common and assumed-band sets for every curve.

## Frozen scoring rule, uncertainty and secondary diagnostics

**Primary score:** at a fixed c, count how many of the four exact attainable
mean ranges have nonempty intersection with their printed marginal 95% ORR
intervals. A shared-c compatibility pass requires all four. The primary
domain score is whether the common exact scale set intersects that domain;
compute it for all positive c and both inherited bands, without a grid
approximation. No optimization over patient distributions is necessary.

Use printed interval endpoints unchanged for the headline result. Their
intersection is a deterministic compatibility rule, **not a simultaneous
95% confidence test, joint confidence region, p-value, equivalence test or
calibrated falsification probability**. No independence assumption, binomial
likelihood, invented counts, meta-analysis, bootstrap or multiplicity claim.

Required secondary checks, all reported regardless of direction:

1. Treat the printed point ORRs as exact target means and compute the same
   sets separately. An empty point set does not establish population failure
   when the printed interval rule passes.
2. Widen each printed interval endpoint by half its last printed percentage
   unit, clipped to [0,100], then recompute the sets. This is only a source
   rounding sensitivity; it neither replaces the primary intervals nor adds
   sampling precision. The bin boundary and frozen coefficients stay fixed.
3. At every inherited grid c, report attainable ranges, interval-compatible
   bin count, point-compatible bin count and the infimum distance from each
   point ORR to its attainable range. Report the root mean of the four squared
   distances in percentage points. A zero distance at an excluded endpoint
   is not attainment; report exact compatibility separately. This lower-bound
   descriptive loss grants arbitrary distributions, is not an observed
   prediction error, and cannot be called a fitted or validated response.
4. Check forward/inverse agreement, finite feasible witnesses, strict lower
   endpoints, reciprocal ceilings, log clipping plateaus, closed singletons
   and empty intersections with independent tests. At a compatible grid
   scale a witness may select an interior mean in each observed interval;
   that is an existence certificate, not an outcome fit or recovered TMB.

No extra empirical sensitivity is added after scores to rescue a result.
Numerical comparisons use documented tolerances; structural endpoint flags
must remain exact and cannot be closed by a floating-point tolerance.

## Verdict, falsification attempts and novelty

Give a primary verdict for `reciprocal_holdout` and explicit verdicts for all
required sensitivities. If an eligible common all-positive region is empty,
say **observation bridge falsified under the declared interval-compatibility
rule**, with the domain and non-significance qualification. If nonempty,
say **survives with stated latitude**, listing admitted scales, any rejected
bands/grid reference, arbitrary-distribution freedom and marginal uncertainty.
If rights, eligibility or essential observations prevent scoring, say
**untestable**, identifying what is missing. Band-specific failure is not
all-positive failure; survival is not confirmation. Missing bin-specific
counts alone do not block this count-free compatibility calculation.

Challenge the common-scale restriction, both inherited bands, the primary
identity reference, all four fixed curves, exact points versus intervals,
source rounding, curve ceilings/clipping, open-bin boundaries, source
denominators, ecological-to-within-bin transfer, selection and evidence reuse.
Keep every adverse, null and unavailable check visible. Record source-audit
corrections and hostile mathematical/software review in README, distinguishing
review before computation from review after results. AI review and CI do not
provide external scientific validation.

This test removes the particular median substitution and broad-ORR/TMB-subset
denominator mismatch in PR #72 by using the printed TMB-bin ORRs. It does
not remove the ecological response bridge, sequencing-subset selection,
assay uncertainty or different treatment contexts. It **does not isolate
viral antigens and does not test PR #59's functional-load theorem**. PR #67's
pan-cancer comparison is not rerun or reversed. Credit the trial authors'
existing bin analyses and prior ecological regressions. Novelty is at most a
repository-specific frozen-curve/observation-map audit with elementary sharp
bounds; no new biology, biomarker, general theorem or bibliographic priority.

## Reproducibility and delivery

Create only `research/merkel-tmb-bin-test/` (README.md, this PROTOCOL.md,
the small CSV under data/, analysis code) and
`tests/test_research_merkel_tmb_bin_test.py`. Read other files without editing.
Use the Python standard library only, offline, with the new test comfortably
below 60 seconds on Python 3.12 and 3.14, Windows and Linux. Recompute every
reported analysis number from the new aggregate CSV and the unchanged frozen
curve/training CSVs; verify a generated README numeric block. Do not require
network access or git history to run tests in a shallow CI checkout.

Hash UTF-8 text after CRLF and CR normalization to LF. Construct **all path
keys with `Path.as_posix()`**, including fingerprints and diagnostic mappings.
Include protocol, aggregate CSV, frozen source dependencies, analysis and the
new test in reported fingerprints; exclude README from its own hash to avoid
circularity. Independent replay coefficient tolerance is `2e-6`; inverse and
score comparisons use absolute probability tolerance `1e-12` or documented
equivalent scaled tolerances. Any optimizer diagnostics have bounded
acceptance checks rather than exact saved-value comparisons.

Commit the protocol before computing; preserve it unchanged thereafter.
Hostilely self-review source transcription, proof, numeric reproducibility,
endpoint handling, claims and exact changed-file scope before pushing.
Run relevant repository validation, Python and Node checks. Push the branch
to the requested repository, open a PR targeting main, observe its CI, and
report the actual PR number and CI state. Do not merge the PR.
