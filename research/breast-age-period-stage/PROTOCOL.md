# Frozen protocol: age-by-period compatibility within breast-cancer stage

**Nonclinical, unreviewed computational research. No screening or medical
advice.** This is a secondary audit of published aggregate trends, not a
clinical prediction or causal analysis.

Protocol date: **2026-11-18, America/New_York**. Branch from current `main`
`2a6d357224a286a05dbae03d6d3ad88ed622d552`. PR #84 is open at design time;
the user's brief supplies the specification for its highest-ranked job.
Read [the research programme](../README.md),
[PR #80's stage analysis](../breast-incidence-decline-stage/),
[PR #73's strata analysis](../incidence-strata-test/),
[DATA_POLICY.md](../../DATA_POLICY.md), and [AGENTS.md](../../AGENTS.md).

Commit this file alone before any transformed-trend, compatibility, allowance,
or sensitivity calculation or execution of result-generating code. Preserve
it unchanged; identify the freeze commit in the README. Source reading,
licence checking, factual transcription, and mathematical design may precede
the commit. The published trends and predecessor results are already visible.
This is a committed computational specification, not blinding, independent
registration, independent replication, or a claim of independent timestamps.

## Question, restriction, and observation bridge

For each fixed recorded stage `s`, ask whether the incidence surface could
satisfy

`I_s(a,t) = A_s(a) * P_s(t)`, with positive age and period factors.

Each stage has its own unrestricted period factor; no equality across stages
is imposed. Under this restriction, all ages share the same log-period
change over the same calendar interval. With fixed nonnegative age weights
within a reported bin, the bin's standardized rate also factorises:
`J_sg(t) = B_sg * P_s(t)`. The source's fixed standard-population convention
motivates this bridge but does not verify all underlying composition and
observation assumptions.

The necessary full-period restriction is a common average log-period change
`beta_sg = beta_s` across the three reported age groups, using

`beta_sg = log(1 + AAPC_sg / 100)`.

The transformed published AAPCs are estimates of those target summaries,
not exact population slopes. The factorisation implication concerns the
targets; the interval rule below audits their published estimates.

Use natural logarithms, in inverse calendar years. The time span is 2001-2019.
This is a transform of the published joinpoint AAPC summary, not a new fit,
an assumption that every annual change is constant, or a reconstruction of
annual observations. Equality of these averages is necessary, not sufficient,
for common annual trajectories or the full pointwise factorisation. Different
within-period paths can have the same average. Estimation and model selection
uncertainty are only represented to the extent captured by the printed
marginal intervals; the observation bridge remains conditional.

## Fixed source, licence, panel, and transcription

Only numerical input: Lee Argov EJ, Lui ML, Karr AG, Tehranifar P, Kehm RD
(2025), *Breast Cancer Incidence Trends in Older US Women by Race, Ethnicity,
Geography, and Stage*, JAMA Network Open,
[doi:10.1001/jamanetworkopen.2025.16947](https://doi.org/10.1001/jamanetworkopen.2025.16947),
[main article HTML](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2835519).
Use main **Table 2**, **Stage at diagnosis**, exactly the **Localized,
Regional, Distant, and Unknown** rows, and **65-74 y, 75-84 y, and >=85 y**
columns. Retain every selected stage and age regardless of the result.
Exclude the overlapping All ages column, In situ row, all demographic and
geographic rows, and supplementary tables. In situ is outside the requested
four-stage panel, not an outcome-dependent exclusion.

Read **Table 1** for rate-summary context and the main Methods for the
definitions. Table 1's period-pooled age-adjusted rates and rate ratios are
not annual series and will not be used to reconstruct slopes. The article
uses female incidence, a fixed US standard population, and a combined
summary-stage variable. All regional forms are combined. Unknown includes
unknown, unstaged, unspecified, and death-certificate-only cases. The HTML
text reader does not expose the table footnote prose; do not claim to have
read it. The Methods supply the stated definitions. Changes in staging
systems are described incompletely for this purpose; no stage-reclassification
history will be supplied by assumption.

The article's **Article Information / Open Access** notice was checked before
transcription: **CC-BY**. Its linked
[JAMA permissions page](https://jamanetwork.com/pages/cc-by-license-permissions)
permits reuse with author and journal credit. Neither readable notice states
a licence version; record **CC-BY, version not stated**, not an invented
version. Preserve attribution and note that the CSV reformats selected facts,
normalizes minus signs and age labels, and adds provenance. The repository's
MIT licence does not relicense third-party material.

Read the source only as article HTML text. Never initiate a browser file
download. If a prompt appears, deny it, continue without that file, and record
the limitation. No PDFs, supplementary files, patient-level records, raw USCS
data, or file downloads are needed. Record access failures or discrepancies.

The small CSV will contain exactly the selected printed AAPCs and both
marginal 95% interval endpoints as decimal strings, with stage and age labels,
source citation, DOI, URL, table, row, column, human-readable locator, period,
units, access date/timezone, and reuse basis/URL on every row. Do not add
invented annual counts, denominators, covariance, age midpoints, reclassification
histories, or derived numerical fields. No Table 1 numbers enter the CSV.

Missing, suppressed, duplicated, nonnumeric, inconsistent, or malformed cells
are errors that stop the run with an explicit untestable/error message; do
not impute or silently select a smaller primary panel. Require ordered
interval endpoints containing the printed estimate and a valid log domain
after rounding widening. All selected cells must be audited against source
HTML before publishing the analysis.

## Frozen rounding rule and primary compatibility decision

For any printed decimal `x` with `d` digits after the decimal point, let
`q(x) = 0.5 * 10^(-d)` in **percentage points of AAPC**. Treat its unknown
unrounded value as lying in the closed interval `[x-q(x), x+q(x)]`.
The selected table uses two decimal places, hence a half-unit of `0.005`
percentage points. This is an assumed last-digit rounding convention, not
a verified publisher algorithm or a probability distribution.

For printed estimate `x`, lower endpoint `l`, and upper endpoint `u`, retain
the point-rounding interval, and form the primary marginal interval hull by
widening the **printed endpoints before transforming**:

`C_sg = [log(1 + (l-q(l))/100), log(1 + (u+q(u))/100)]`.

The point-rounding interval is
`[log(1+(x-q(x))/100), log(1+(x+q(x))/100)]`; it does not replace the printed
sampling-uncertainty interval. Rounding uncertainty and the source's
marginal 95% intervals have distinct meanings. Do not take symmetric
standard errors from the interval widths or round a derived result before
deciding compatibility.

Within each stage, define `L_s = max_g lower(C_sg)` and
`U_s = min_g upper(C_sg)`. **PASS** means `L_s <= U_s`; **FAIL** means
`L_s > U_s`, an empty intersection. Boundary equality is compatible.
For decisions, exploit the monotonicity of log to compare the exact decimal
percent endpoints; use sufficient numerical precision for the reported log
values. Report both limiting endpoints and all age groups attaining each.
For a pass, `[L_s,U_s]` is the common-average compatibility set. For a fail,
the reversed limits are not an interval: report the separating bounds and
the positive gap explicitly.

This **interval intersection is a descriptive compatibility rule, not a
simultaneous significance test**. There is no calibrated overall error rate,
joint confidence region, p-value, or declaration of statistical significance
for a cross-age contrast. Marginal source coverage is not simultaneous
coverage over ages or stages, and joint covariance is unavailable.

## Sharp additive allowance and overall rule

Define the effect size as the smallest common nonnegative allowance `epsilon`
added symmetrically to **each** log-trend interval, permitting each age's
trend to lie within `epsilon` of a common reference `b_s`:

`epsilon_s = inf_b max_g distance(b, C_sg)
           = max(0, (L_s-U_s)/2)`.

Report `epsilon_s` in natural-log units per calendar year, the full positive
gap `max(0,L_s-U_s)`, and an attaining reference `b_s=(L_s+U_s)/2`.
The allowance is half the gap because both opposing interval ends may move;
it is not a one-sided shift or an AAPC percentage-point difference. A pass
has zero allowance even if printed point estimates differ. A single common
translation of all intervals cannot repair separation.

Proof obligation: after widening by `epsilon`, compatibility is exactly
`L_s-epsilon <= U_s+epsilon`. This proves both the lower bound and attainment.
For a failing stage the midpoint attains the bound; any smaller symmetric
allowance fails. This is a set-distance summary conditional on the chosen
intervals, not a measured biological interaction, confidence bound on true
heterogeneity, or fitted screening effect.

**Overall verdict:** FAIL the joint descriptive restriction if any selected
stage fails; PASS only if all selected stages pass. This conjunction allows
a different reference trend for every stage. Report each individual verdict,
including passing stages, and the overall verdict. The smallest allowance
shared across stages is `max_s epsilon_s`, with separate `b_s`. An overall
pass would still not validate the full factorised model. An overall fail
identifies an incompatible published-summary restriction under this rule,
not an independently calibrated falsification of biological separability.

## Prespecified attempts to falsify the conclusion

1. Independently re-read every selected numeric cell, signs, age headers,
   period, and stage label against the main article HTML. Compare duplicated
   HTML representations when exposed; they remain one source, not replication.
2. Repeat all four-stage calculations using the printed CI endpoints with
   no widening, and using **one full last-digit unit** on each endpoint
   (`2*q`) as a rounding stress. The latter exceeds the primary assumption;
   it does not assert publisher error. Retain every changed result.
3. Calculate compatibility of the point-rounding intervals alone for every
   stage. This intentionally omits the source uncertainty and cannot overrule
   the primary interval verdict or support significance.
4. Leave out each age group in turn for every stage, recomputing the primary
   interval rule on the remaining pair. Report all results. These are
   diagnostics of which ages constrain the intersection, not exclusions from
   the primary panel or a causal localization of an age effect.
5. Check allowance attainment and failure below it using independent
   high-precision calculations; check equality, ties, a common shift, malformed
   data, and logarithm-domain boundaries on clearly labelled synthetic inputs.
6. Try to invalidate causal/replication/novelty language by comparison with
   the article's own findings and the predecessor notes. Report limitations
   rather than manufacture missing observations or broader claims.

No result-dependent extra source, stage, age, weighting, outcome, significance
test, or rounding rule will replace the frozen primary specification.

## Interpretation, novelty, and reproduction commitments

The published incidence trends are **prior findings**. Any increment is an
auditable restriction/rounding/allowance calculation, not discovery of those
trends, a new mathematical theorem, a systematic literature review, or a
claim of bibliographic priority. Record adverse and compatible findings.

Period slopes here differ from PR #80's **cross-sectional age slopes**:
that note compares age groups within one calendar year using finite
five-year bins and a conditional count model. This audit compares calendar
trends within broader age groups. It imports neither that note's Poisson
model nor its simultaneous calibration. The open-ended 85+ group can supply
its reported period trend without any midpoint or assumed bin width.
Shared USCS populations, overlapping calendar coverage, and different age
bins prevent an independent replication claim or attribution to age, cohort,
screening, detection, competing risks, stage migration, or biological causes.
Unknown stage is retained as an observation category, not resolved into true
stages. No annual counts, joint covariance, 85+ midpoints, or stage-
reclassification histories may be invented.

Create only this directory (README, this frozen PROTOCOL, small aggregate CSV
under data/, analysis code) and
`tests/test_research_breast_age_period_stage.py`. Do not modify existing files,
the central research index, workflows, or other research directories. All
reported analytical numbers must be recomputed from the CSV. Use standard-
library Python, no network in analysis/tests, and a dedicated test suite
that completes in under 60 seconds. Include independent arithmetic, source
and panel checks, all generated README numbers, frozen protocol/data hashes,
and newline/platform portability. Hash UTF-8 content after CRLF and CR are
normalized to LF; construct **all path keys with `Path.as_posix()`**, with
portable Windows/POSIX checks. Support Python 3.12 and 3.14 without optional
packages. Distinguish actually executed environments from compatibility
reasoning. Hostilely review the complete work before pushing it, open a PR
against main, and report the actual CI status. Software checks and AI review
do not establish source truth, independent science, or clinical validity.
