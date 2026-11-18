# Frozen protocol: older-age breast incidence and stage decomposition

**Nonclinical, unreviewed computational research. No screening or medical
advice.** This protocol specifies a secondary analysis of published aggregates,
not a clinical prediction tool or a causal screening study.

Protocol date: **2026-11-18, America/New_York**. Base: current `main`
`09d916cdea7053f25643c6d0bd9c45827436e315`, which includes merged PR #77.
Commit this file alone before executing any reconciliation arithmetic,
statistical contrast, stage decomposition, or result-generating code. Source
reading, transcription, and method review precede this commit. Preserve this
file unchanged; the README will identify the commit and disclose deviations.
The published age pattern and earlier repository results are already known:
this is a prospective computational specification, not blinding, independent
registration, or independent data replication.

## Question and permitted inputs

Run ranked job 2 in [the research programme](../README.md), following
[PR #64](../multistage-incidence-identifiability/),
[PR #73](../incidence-strata-test/), [DATA_POLICY.md](../../DATA_POLICY.md),
and [AGENTS.md](../../AGENTS.md). Ask which adjacent-age incidence declines
survive sampling uncertainty and whether recorded stage explains their
arithmetic contribution. Biological transition counts remain unidentified.

Use only Turbow SD, White MC, Breslau ES, Sabatino SA (2021),
*Mammography use and breast cancer incidence among older U.S. women*,
Breast Cancer Research and Treatment 188:307-316,
[doi:10.1007/s10549-021-06160-4](https://doi.org/10.1007/s10549-021-06160-4),
[PMC HTML Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10846538/#T1).
The table uses 2016 invasive-breast USCS incidence, November 2018 submission,
released June 2019. Read HTML text; never use a browser file download. Deny
any download prompt and record that limitation. Search/PDF snippets are not
authoritative transcription evidence. Do not acquire case-level USCS data.

Transcribe the printed overall population denominators and the Overall,
Localized, Regional, Distant, and Unknown-stage case counts, rates, and
95% interval endpoints for ages 65-69, 70-74, 75-79, 80-84, and 85+.
Store only those aggregate facts in a small independently structured CSV,
with source, DOI/URL, table, row, age column, denominator locator, access
date/timezone, units, adjustment label, dataset version, and reuse basis on
every row. Do not include derived counts or silently repair source cells.
Reuse factual statistics with attribution, without reproducing article prose,
figures, or table layout, claiming a CC/public-domain article license, or
relicensing third-party material under the repository MIT license.

Exclude race/ethnicity and receptor strata: their matching denominators and
race-by-stage cross-classifications are not printed. Do not substitute overall
denominators or infer denominators by dividing counts by rounded rates. The
separate NHIS screening survey is not an analytical input. No new data source
or outcome is selected after seeing the computed contrasts.

## Reconciliation and stopping rules

For every selected cell calculate `r = 100000 * cases / population` using
the printed population as the fixed one-calendar-year exposure proxy. Compare
it with the published rate, its closed half-last-digit rounding interval
`[printed - 0.05, printed + 0.05]`, and the printed confidence interval.
Record the numerical difference and whether it fits rounding. Do not use
published interval overlap as an adjacent-age test or treat rounding as
sampling uncertainty. The original intervals are retained for audit; their
calculation method is not assumed or reverse-engineered.

Retain Table 1's age-adjusted label and its footnote: 2000 U.S. standard,
19 age groups, Census P25-1130. Check the official standard and document the
single-standard-cell interpretation where supported. Our estimand remains
the count/population bin rate. Agreement after rounding supports compatible
arithmetic, not independent verification of SEER*Stat's unpublished session.
If a cell differs beyond rounding, report the discrepancy and restrict that
cell's inference to the explicitly substituted crude rate; do not call it a
test of the published adjusted rate. If exposure is absent/nonpositive or a
count is missing/noninteger, mark affected contrasts untestable without
imputation. Invalid confidence endpoints or duplicated cells are errors.

Compute `residual = overall - localized - regional - distant - unknown`
for every age. Retain the signed residual as an unassigned accounting term.
A positive residual is not a verified new stage; a negative residual means
the partition fails and the completed-stage allocation analysis abstains.
Report source prose/count contradictions explicitly and prefer the checked
counts for arithmetic without correcting the article.

## Frozen contrasts, sampling law, and simultaneous uncertainty

Finite bins are `[65,70), [70,75), [75,80), [80,85)`. The primary contrasts,
all oriented **older / younger**, are exactly:

1. 65-69 -> 70-74;
2. 70-74 -> 75-79;
3. 75-79 -> 80-84.

Apply each to Overall, Localized, Regional, Distant, and Unknown: **15
contrasts**, with no selection of only observed decreases. Keep **85+ purely
descriptive and outside this inferential family**. Do not assign it a midpoint,
pretend it has a five-year width, or use it to claim a finite-bin replication.

Assume independent Poisson counts across disjoint age groups, conditional on
fixed published exposures: `Y[a,s] ~ Poisson(E[a] * lambda[a,s])`. Recorded
stage categories may be represented by a marked Poisson process; overall
overlaps stage counts. Do not assume independence across overlapping endpoints
or adjacent contrasts. Bonferroni below does not require that independence.
This is process-model uncertainty for registry aggregates, not survey sampling
error: the registries aim to enumerate cases. It ignores population-estimate
error, coding/registration error, clustering, extra-Poisson variation, and
unobserved heterogeneity; no inference is unconditional on those assumptions.
Population is not a measured cancer-free person-time risk set.

For younger count `y`, older count `o`, and exposures `E_y,E_o`, condition on
`n=y+o`. Then `o | n ~ Binomial(n,p)` where
`p = E_o*lambda_o/(E_y*lambda_y+E_o*lambda_o)`. Transform the exact
Clopper-Pearson bounds for p to the rate ratio by
`RR = (E_y/E_o) * p/(1-p)`.

Use **two-sided Bonferroni simultaneous coverage at least 95%** over the
15 comparisons: each interval has error `0.05/15`, split equally between
tails. The primary decline criterion is strictly `RR_upper < 1`; increase
is `RR_lower > 1`; otherwise direction is unresolved. Equality is unresolved.
Report RR, both simultaneous endpoints, signed older-minus-younger rate
difference, and direction for all comparisons. Optionally report the matching
central exact two-sided p-value `min(1, 2*min(P[X<=o],P[X>=o]))` and its
Bonferroni adjustment, using `p0=E_o/(E_y+E_o)`.

As a prespecified **separate sensitivity**, repeat decline-only upper bounds
with all `0.05/15` in the upper tail (one-sided Bonferroni). Report every
changed decision and preserve the two-sided primary verdict. This sensitivity
has its own familywise guarantee, not joint coverage with the primary family.
Do not change alpha, multiplicity, or tail choice in response to results.
Use exact-tail inversion in standard-library Python, with an independent
numerical tail check; no asymptotic normal approximation is primary.

## Stage accounting and unknown-stage sensitivity

Define loss `D_s = r_younger,s - r_older,s`; positive values contribute to an
overall decline. By retaining the residual rate,
`D_overall = D_localized + D_regional + D_distant + D_unknown + D_residual`.
For positive total loss, report each signed share `D_s / D_overall`.
Shares may be negative or greater than 100%; do not renormalize them or
discard offsets. If total loss is nonpositive, do not call a share a fraction
of decline. Define **mainly localized** as localized contributing strictly
more than half the observed net loss, conditional on a robust overall decline.
This is an arithmetic, descriptive majority rule, not a confidence bound on
the share. Localized stage is anatomic extent, not observed screen detection.

For each localized/regional/distant stage, treat
`U[a] = unknown[a] + nonnegative residual[a]` as a possible unassigned pool.
Under the additional assumption that the overall and known-stage counts are
truthful and these unassigned cases belong among the named stages, a completed
stage count lies in `[known[a,s], known[a,s]+U[a]]`. No missing-at-random
assumption is made. Each stage is stress-tested separately: do not allocate
the entire pool to multiple stages in the same alleged completed table.

For every finite adjacent pair, form the union of exact RR intervals over all
permitted completions. Monotonicity gives its lower endpoint by assigning all
younger unassigned cases to the stage and none older, and its upper endpoint
by assigning none younger and all older. Retain the same two-sided tail
allocation `0.05/(2*15)`. A decline survives allocation plus sampling only
when this **worst-case upper endpoint is below one**. This is a separate
conditional completed-stage sensitivity family; it does not verify that
unreconciled residuals are real omitted cases or correct source mistakes.

For a positive total loss also report the sharp descriptive localized-share
range obtained from those two assignments, keeping the observed total fixed.
Call the majority allocation-robust only if its lower share is above one-half.
These share ranges handle allocation, **not sampling uncertainty**. Report
nominal stage patterns even if allocation can erase them; a failure to reject
is neither equality nor proof that a decline is absent.

## Falsification, novelty, reproducibility, and scope

Try to falsify each proposed conclusion with cell-by-cell source checking,
rate/interval reconciliation, exact stage summation, all frozen edges, the
two tail conventions, adverse missing-stage allocation, and separation of the
open-ended bin. Record negative and inconclusive outcomes. Source conflicts
cannot be fixed by a favorable model or passing CI. Distinguish finite checks
of software from biological evidence.

Compare explicitly with Turbow et al.'s own age and stage findings. A new
discovery of the decline is not claimed; any increment is a reproducible
uncertainty/decomposition audit. USCS includes SEER and is not an independent
replication of the older SEER-based notes. Different period/geography and
missing race denominators prevent literal replication of their full panel.

No result identifies screening effects causally, birth-cohort effects,
competing-mortality mechanisms, latent onset hazards, or biological transition
counts. The separate incidence and screening datasets cannot be linked
ecologically to supply those missing observations. Give no screening, medical,
diagnostic, surveillance, or treatment advice.

Only create this research directory (README, frozen protocol, aggregate CSV
under data/, analysis code) and
`tests/test_research_breast_incidence_decline_stage.py`. Existing files,
including the central index and all other research notes, stay unchanged.
Generate every displayed analytical number from the CSV. Include SHA-256
file hashes of UTF-8 content with CRLF and CR normalized to LF; construct every
path key with `Path.as_posix()`. Tests must be offline, standard-library only,
and complete in under 60 seconds. Check the full generated README block,
frozen-protocol hash, numerical inversion, decomposition, allocation extrema,
source provenance, scope, and newline/platform portability. Hostilely review
before pushing; open a PR against main and report actual CI status.
