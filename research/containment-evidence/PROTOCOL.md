# Analysis plan fixed before evidence extraction

**Nonclinical, unreviewed evidence assessment; no treatment advice.**

Protocol date: 2026-11-17 (America/New_York). Base: current `main`,
`4f2098f543eb36d2810809af600c4372673794e7`. This file is committed before
extracting study results or running this analysis. It is an analysis-plan
commit, not prospective registration of already published experiments.
The named papers and PR #61's conclusions were known when planning.
Any later departure must be documented explicitly; this file is not rewritten.

## Search and inclusion

Search primary papers, publisher full texts, supplements and trial registries
through the search date. Start with Zhang 2017/2022, Enriquez-Navas 2016,
Gatenby 2009, Smalley 2019, ANZadapt, and Strobl/West/Brady fitted cohorts;
search their references/citing studies and negative intermittent-treatment
comparisons. Include clinical, animal and in-vitro comparisons of adaptive,
evolutionary or containment treatment, and dose-modulation/intermittent
comparisons relevant to the claim. Include fitted group-level quantities but
keep fitted counterfactual outcomes distinct from observed comparisons.
Generic dose finding, toxicity-only reductions, reviews and pure simulations
do not supply independent experimental outcomes. Log exclusions, unavailable
full texts, absent results, licensing failures, duplicate cohorts and queries.
Search retrieval cannot prove that every worldwide publication was found;
report actual coverage and unresolved candidates instead of claiming that.

Only published aggregates are eligible. No individual human/animal trajectories,
patient-level supplements, reconstructed survival records, raw sequencing,
controlled data or graph digitisation. No new aggregation of person-level
data. Use small factual tables only where reuse terms permit; otherwise retain
metadata and a link with the exclusion reason. Record source, exact locator,
access date, reuse basis and uncertainty for every numerical row. Do not
reinterpret censored/not-reached medians as exact values. Deduplicate by cohort
and use the latest eligible report for an endpoint while retaining earlier
reports visibly. No outcome-driven selection of studies or endpoints.

## Parameter mapping and hard eligibility gates

PR #61 uses d=1/g0, c=min(L,1/gC) (or L for gC<=0),
J=L*(1+delta/rS), A=integral bS*sM(t)*pM(t) dt and lambda=N*mu.
For each study classify each parameter as measured, estimated (with its
reported uncertainty), assumed/conditional, or unidentifiable. Net tumour
growth does not identify birth/death turnover, mutation probability,
absolute viable cell number, establishment or competition. TMB is not mu;
PSA/volume is not N without a validated calibration. A fitted Lotka-Volterra
coefficient is not automatically PR #61's dimensional alpha.

Assess zero initial resistance, constant resistant drug insensitivity,
common bounded kill law, sensitive-reservoir control, common finite course,
log-N regime and ordinary uncensored mean total-burden progression. Do not
convert medians, PFS hazards, PSA endpoints or survival into that mean.
Do not infer high lambda from a favourable outcome and then use it to test
the same prediction. A biological prediction is untestable if its required
measurements or comparable endpoint are missing, even when the software runs.

## Fixed prediction family and estimands

1. **P1, mutation-supply slope:** at otherwise comparable model parameters,
   the mean containment advantage increases with independently measured
   lambda. Require at least three independent supply levels with common
   design, mean outcomes and uncertainty; estimate a slope with an interval.
2. **P2, competition:** containment advantage increases with competition,
   at fixed supply and other model parameters. Treat this as an additional
   biological hypothesis: alpha changes both c and A, so PR #61 does not
   assert universal monotonicity in alpha. Require independently measured
   competition contrasts and comparable outcomes. Separately assess P2a:
   no containment-winning region at measured zero competition under A1-A6.
3. **P3, quantitative threshold:** compare independently estimated lambda
   and lambda_star=log((J-d)/(c-d))/A only when d<c<J, A>0 and model gates
   hold. Propagate reported joint uncertainty or conservative bounds;
   never borrow synthetic PR #61 values. Report unknown if any necessary
   parameter is unidentified. Finite-N discrepancies alone cannot falsify
   the limit theorem without an independently checked finite-size model.
4. **P4, establishment:** under maximum dose, log(no established new
   resistant lineage)=-lambda*A. Require aggregate lineage counts, initial
   resistance exclusion and independent lambda/A; use binomial intervals.
5. **P5, one crossing:** require a controlled supply series with at least
   three informative levels and independent rate measurements; assess
   whether two reliably opposite crossings contradict the fixed-parameter
   prediction. Do not treat different cancers as a mutation-dose series.
6. **P6, empirical transport check:** describe observed control/progression
   outcomes and dose exposure for all eligible comparisons. A favourable
   direction is only consistent with possible containment benefit, never
   threshold validation. Explicitly test the broader universal claim that
   all adaptive/intermittent schedules improve control; any eligible
   randomized comparison with an interval wholly favouring continuous
   treatment contradicts that extension, not necessarily PR #61.

Each verdict is exactly `consistent`, `falsified` or `untestable`, with
scope and the cited supporting rows. `Consistent` means compatibility,
not proof, equivalence, clinical efficacy or independent validation.

## Calculations, uncertainty and multiplicity

Report within-comparison median differences and ratios descriptively; retain
published HRs and confidence intervals with their stated orientation. Do
not infer an HR from medians or construct a CI from a p-value. If separate
median confidence limits exist, endpoint-bound ratios are sensitivity bounds,
not confidence intervals for the ratio. Missing uncertainty stays missing.
For fixed-time aggregate counts, use Wilson intervals and a conservative
difference interval; do not treat censoring as observed failures. Dose-time
fractions and cumulative dose are separate estimands. Do not pool unlike
endpoints, species, schedules or overlapping cohorts. No meta-regression on
cancer labels, TMB or inferred mutation supply; no overall efficacy p-value.

The mechanistic family P1-P5 is fixed. If adequate independent data support
formal tests, use Holm correction across tested mechanistic hypotheses at
family alpha 0.05 and show both unadjusted and adjusted results. Otherwise
report no hypothesis-test p-values and explain why. P6 effect intervals are
source-reported descriptive intervals, not a new multiple-testing exercise.
Do not count multiple endpoints/cohort updates as replications.

## Execution and audit

Use standard-library-only deterministic analysis reading only the small CSVs.
Recompute every result table and verdict, and test data provenance, identities,
cohort exclusions, missingness gates, numerical calculations, README freshness
and newline-normalised hashes. The new dedicated test must finish in under
60 seconds and make no network calls. Recheck each extracted value against
its primary source before pushing; record corrections, negative evidence,
source-access limits and inability to identify parameters. Passing CI or
agreement between AI reviewers is not scientific validation.
