# Prospective analysis plan: cancer-type immune-escape test

**Nonclinical, unreviewed computational research. No treatment advice.**

Frozen before new source retrieval or numerical extraction for this task.
Protocol date: 2026-11-17, America/New_York.
Base main: `2b03737d05c5647587a70671ab10bb4bc969cff3`.
The first branch commit contains only this plan. Read first: PR59
(`research/immune-escape-neoantigen-clonality/`), PR63
(`research/prediction-confrontation/`), AGENTS.md and DATA_POLICY.md.
Those notes already disclose published associations and contrary HLA findings;
this is a prospective plan for a retrospective ecological analysis, not a
blinded experiment or registration in an independent study registry. Publish
this protocol-only commit before the new search. Preserve it unchanged; record
subsequent decisions and their timing in AMENDMENTS.md.

## 1. Prediction fixed before collection

PR59 Theorem 6 gives, under source suppression, fixed rates, no initial
presentation-null cells and increasing minimum functional per-cell load,

`Q ~= exp[-(sigma*h/kappa) * sum_i(1/k_i)]`.

For homogeneous load k and effective founder count N this is `exp(-D*N/k)`.
To obtain a testable cancer-type surrogate, explicitly ASSUME:

* Comparable median TMB T is proportional to potential antigen burden.
* A reported clonal/truncal neoantigen fraction f supplies a coverage proxy
  `k = c*T*f`; N/c and sigma*h/kappa are constant across types.
* ORR in a non-LOH population is `A0*Q`, with a common ceiling A0 and zero
  false-positive response bridge. This is a restrictive working surrogate,
  not a claim that RECIST response equals eradication.
* A type's HLA-LOH prevalence L describes a population mixture. The LOH
  component has the same burden shape but ceiling A1, with
  `0 <= A1 <= A0 <= 1`. L is NOT a transition hazard or an initial null-cell
  count. Ordinary HLA LOH is NOT complete antigen-presentation loss.

The full, three-parameter working prediction is

`p(T,f,L) = [A0*(1-L) + A1*L] * exp[-B/(T*f)]`, `B >= 0`.

Equivalently `p=A0*(1-eta*L)*exp[-B/(T*f)]`,
`eta=1-A1/A0 in [0,1]` when A0>0.
The mathematical theorem fixes neither A0 nor eta and does not quantify an
ordinary-LOH penalty. B is a cross-type fitted composite, not an estimate of
sigma, h, kappa or N.

**Sharpest falsifiable shape consequence:** after division by the fitted
mixture ceiling C(L), `log[p/C(L)] = -B/(T*f)`; doubling T*f at fixed L gives
`p(2*T*f,L)/C(L) = sqrt[p(T*f,L)/C(L)]`. At fixed T and L, increasing f
cannot decrease p; at fixed T and f, increasing L cannot increase p.
The reciprocal-power exponent is one. These are consequences of the
theorem PLUS the measurement, common-parameter and endpoint assumptions.
PR59 Theorem 5(i) already proves that clonal-ID fraction alone does not
identify the harmonic per-cell load; f is not a sufficient statistic.
A contradiction rejects this ecological surrogate, not the theorem.

## 2. Data, source ordering and stopping

Use only small, already-published cancer-type aggregates. No patient-level,
controlled-access, re-identifiable or reconstructed individual data, no
patient-level downloads to compute aggregates, no raw source workbook of
unknown/mixed granularity. Read primary papers and aggregate-only supplements
or explicitly licensed published aggregate republications. No inference of
unreported denominators from rounded percentages.

Search the named anchors: Yarchoan 2017 and its cited response sources,
Chalmers 2017, McGranahan 2016/2017, Litchfield 2021, Montesion 2021 and
TRACERx; follow relevant primary updates and contrary findings. Source
priority for the main paired TMB-ORR panel: the original Yarchoan 2017
published numeric panel if accessible with an affirmative reuse basis;
otherwise a published, reusable aggregate update/republication with explicit
provenance back to the original trial summaries. Prefer a coherent panel
over picking a favorable study per type. Choose by source availability,
endpoint comparability and aggregate coverage, never by fitted performance.
Keep source lineage and overlap explicit. A later independent panel, if
available under the same rules, is a sensitivity, not an independent patient
replication if trial cohorts overlap.

Eligible outcomes are objective response to anti-PD-1/PD-L1 monotherapy,
in the original source's stated populations. Exclude combination-treatment,
anti-CTLA-4-only and survival-only endpoints from the ORR model.
Keep PD-L1/MSI/viral selection and line of therapy visible. Retain zero ORRs.
Do not manufacture an unselected population from selected cohorts.

For TMB, retain source assay, units and median definition; mean and median
are not interchangeable. Prefer paired source-panel medians, and audit
Chalmers as a source/cross-check. For f prefer reported median clonal/truncal
neoantigen fractions per type, then reported aggregate fraction with the
different estimand labelled. Clonal mutation fractions are a separately
labelled exploratory proxy, never silently equated to neoantigen fractions.
For L use reported type-specific any-HLA-I-LOH prevalence in the broad
Montesion pan-cancer series when reusable/available, then a coherent
alternative pan-cancer series. Lung-only TRACERx is contextual/sensitivity.
Do not mix any-HLA LOH with HLA-A-only frequencies in the same primary panel.

Include every eligible numeric type from the selected paired panel;
missing f/L does not remove it from the TMB-only analysis.
Use exact cancer categories where possible. Record and flag any coarsening,
histology substitution, shared parent-type predictor or biomarker-stratum
mismatch. A pan-cancer pooled f or L cannot substitute for type-specific
values. No imputation of missing f or L. Mapping decisions must precede
fitting. A sensitivity excludes all approximate mappings.
Related strata/histologies of the same parent cancer are grouped for
leave-one-cancer-type-out evaluation, preventing leakage across MSI strata
or lung histologies. Also report row count and independent parent-group count.

Printed numeric tables/text take priority over aggregate figure digitization.
Only aggregate summary markers may be digitized; record method/resolution,
rounding/reading uncertainty and exact panel. Never digitize patient points.
Keep an inclusion/exclusion/search ledger covering all named anchors,
rights/access failures, incompatible endpoints and missing aggregates.
Stop after these anchors, their relevant aggregate references, and targeted
searches for missing cross-type clonal and HLA summaries have documented
outcomes. This is a targeted audit, not an exhaustive systematic review.

Every source-data CSV row must have source title/DOI or identifier, URL,
exact locator, access date, reuse basis/terms URL, type/population, definition,
units and notes. Keep source-specific license attribution; repository MIT
does not relicense third-party data. Without established rights, retain
bibliographic metadata only. Source checks and rights checks are manual
evidence audits; unit tests and hashes cannot certify their truth.

## 3. Models and prespecified sample gates

Equal weight per published cancer-type row for fitting. Equal weight per
independent parent cancer in prediction scoring (average within parent first).
Report descriptive ungrouped row scores only as a sensitivity.

* Intercept-only reference.
* Linear TMB: `p=a+b*T`.
* Log-linear TMB: `p=a+b*log(T)` (natural logarithm).
* Reduced reciprocal shape: `p=A*exp(-B/T)`, `0<=A<=1, B>=0`.
  Here f is assumed constant and absorbed into B; this does not test clonality.
* Clonal reciprocal: `p=A*exp[-B/(T*f)]`, same constraints, if at least
  five distinct parent types have eligible f.
* HLA reciprocal: `p=[A0*(1-L)+A1*L]*exp(-B/T)`, with
  `0<=A1<=A0<=1`, if at least ten distinct parent types have eligible L.
* Full prediction in section 1 if at least ten distinct parent types have
  both eligible f and L. Baselines must be refitted on the identical subset.

Do not describe an unavailable full model as tested or substitute f=1.
If a gate fails, report exact coverage and missing requirements while still
performing the real reduced shape test and all eligible reduced extensions.
A small available table of f remains part of the evidence, not a completed
pan-cancer coefficient test.

Estimate by ordinary least squares on the ORR probability scale, including
zeros. For linear baselines use ordinary least squares then clip predictions
to [0,1] for comparable probability scoring; disclose clipping. Use bounded
least squares for reciprocal models, variable projection of the ceiling
coefficients for fixed B, a deterministic broad one-dimensional B search
and local refinement. Include boundary solutions and verify the minimum
against a wider grid. No outcome-driven model or hyperparameter selection.

Prespecified diagnostics, clearly secondary: relax A1<=A0 but retain
0<=A0,A1<=1 to see whether the HLA sign is opposite; fit the reciprocal
shape exponent q in `exp(-B/T^q)`, q in [0.05,4], to assess the fixed q=1
restriction; when at least ten eligible f types exist, allow a clonal
exponent theta in `T*f^theta`, theta in [-2,3], to assess positive sign and
the fixed theta=1 magnitude. A boundary estimate or broad/unstable interval
means weak identification. Never select a preferred exponent after inspecting
the outcomes. These diagnostics cannot rescue a failed prespecified model.

## 4. Prediction, uncertainty and decision

Refit all parameters from scratch holding out one complete parent cancer
group at a time. No held-out data in training or source mapping. Primary
score: parent-balanced mean squared error and its square root (RMSE).
Secondary: parent-balanced MAE and calibration mean residual.
Primary comparison is full theory versus log-linear TMB on identical
complete cases if its gate is met; otherwise the reduced reciprocal shape
versus log-linear TMB is the explicitly limited primary available comparison.
Linear TMB remains a required simpler comparator. Show all models, training
fit and held-out predictions; training R-squared is not predictive evidence.

Uncertainty: cancer-group bootstrap with a fixed seed, resampling whole
parent types, for fitted parameters where estimable; target at least 1000
replicates, disclose skipped singular samples. For paired held-out errors,
use 10000 deterministic bootstrap resamples of parent error differences
and report percentile intervals. These intervals describe variation across
this selected set of cancer types conditional on the extraction; overlapping
training folds, few types and source selection limit their inferential
interpretation. They are not patient-level confidence intervals or formal
external validation. Do not infer binomial sampling precision from unknown
pooled-trial denominators. Use published ORR intervals/denominators only
where the actual aggregate construction supports them.

Report source rounding/digitization uncertainty separately from sampling
uncertainty. Prespecified sensitivity: deterministic perturbations within
reported rounding/reading envelopes, exact mappings only, removal of
biomarker-selected strata, and all single-parent deletions.
If a sensitivity is impossible from available metadata, say so.

A predictive improvement requires lower held-out MSE than BOTH simple TMB
baselines and a paired-bootstrap interval favoring the theory against the
log-linear primary comparator. A negative/unstable difference means no
demonstrated added predictive value. Do not turn failure to reject into
equivalence or validation. This is a descriptive ecological criterion, not
a multiplicity-controlled causal significance test. HLA/clonality diagnostics
are secondary; report all, avoid families of opportunistic p-values, and make
no significance claims from nominal intervals across multiple diagnostics.

## 5. Outliers, contradictions, novelty and reproducibility

Rank parent cancers by absolute held-out residual under the primary available
theory and log-linear baseline; report the five largest (or all if fewer).
Never drop an outlier because a mechanism can be proposed after fitting.
For each, propose only a minimal sourced mechanistic hypothesis, flag that it
is post hoc, and state an aggregate measurement that could distinguish it.
Separate model contradictions from unavailable tests and bridge failures.
Explicitly discuss ecological fallacy, cross-cohort TMB/ORR/f/L mismatch,
antigen quality/expression, viral antigens, immune infiltration/exclusion,
HLA partial loss, tumour burden N, selection and shared-trial dependence.

Novelty audit must credit earlier ecological TMB-ORR regressions and
clonality/HLA-response work. A different nonlinear curve and reproducible
negative result are not a novel biomarker discovery. No AI agreement,
reproducibility or successful CI is biological validation.

Only create `research/immune-escape-pancancer-test/` and
`tests/test_research_immune_escape_pancancer_test.py`; no other repository
file may change. Analysis and new tests use the Python standard library,
offline, and the dedicated test runs in under 60 seconds. Recompute all
reported numeric results from the CSVs, compare generated README numeric
blocks, check source/range/missingness/mapping guardrails, and normalize
CRLF/CR to LF before every file hash. Run the repository's required checks.
Hostilely verify every transcribed number against its primary source,
preserve exclusions and contrary evidence, then publish the reviewed results
on this branch, open a PR, and observe CI.
