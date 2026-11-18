# Swedish HPV age-at-first-dose / time-since-dose factorisation

**NONCLINICAL, UNREVIEWED SECONDARY RESEARCH. No vaccination, screening,
diagnosis, treatment, or clinical advice.** This note is not a new clinical
finding, external scientific review, or a causal study.

## Question and frozen approach

Could Wu et al.'s **fully adjusted** relative-to-unvaccinated invasive
cervical-cancer IRRs from initiation ages **10-16** versus **17+** have a
single constant ratio over the matched **7-9, 10-12, and 13-15 years since
first dose** bins? This is a necessary *printed-estimate compatibility*
condition, not a formal regression interaction test. Every printed 95%
interval is marginal and the adjusted covariance is unpublished here.

The [frozen protocol](PROTOCOL.md) was the only file added in first commit
[`b1401db41cbca4b0fd015656c216fb1a38456b9f`](https://github.com/dovvnloading/open-cancer/commit/b1401db41cbca4b0fd015656c216fb1a38456b9f),
prior to any result calculation. Branch base was main at
`307aa3c336548046c344772f3bedd4f325a1b1a0` on 2026-11-18.
[Attribution and raw transcribed facts](data/wu_figure3.csv),
[offline script](analysis.py), [full report](REPORT.md), and
[independent tests](../../tests/test_research_hpv_age_duration_factorisation.py)
follow that freeze. We did not modify existing files.

**Primary descriptive verdict: CONSISTENT.** The ratio interval boxes,
widened for every printed last digit, intersect in **[0.173553719008,
1.135922330097]**. The sharp symmetric log-ratio allowance is **0**.
This is not a finding of one true biological constant. Printed point ratios
differ; their rounding-only bands do not overlap (point-only **FALSIFIED**,
sharp log allowance **0.547566255671**). With marginal confidence limits
and unknown joint covariance, the *statistical interaction inference* is
**INCONCLUSIVE**. No-rounding and full-unit widened 95% CI projections,
all three leave-one-bin-out checks, and an explicitly approximate six-IRR
Bonferroni/log-Wald sensitivity all remain **CONSISTENT**. For full exact
fractions, bounds, sparse-cell limits and limitations, see [REPORT.md](REPORT.md).

## What the literature already establishes and what is new here

Source: Shiqiang Wu, Yunyang Deng, Tiia Lepp et al.,
"Extended follow-up of invasive cervical cancer risk after quadrivalent HPV
vaccination: nationwide, register based study," *BMJ* (2026);392:e087326,
DOI [10.1136/bmj-2025-087326](https://doi.org/10.1136/bmj-2025-087326),
25 February 2026. Only the [article HTML](https://www.bmj.com/content/392/bmj-2025-087326.full),
**Results > Invasive cervical cancer by time since vaccination**, **Fig. 3**,
and its Methods, supplies the adjusted IRRs, printed marginal CIs and
contextual vaccinated count/person-year pairs. Article's footer explicitly
states [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); numerical
facts are sparse, attributed, and not raw register data. Our transformations
are new and not supplied by the source; repo MIT does not relicense data.

Wu et al. already publish the duration-specific estimates, lower invasive
cervical-cancer incidence among vaccinees, age-at-vaccination differences,
and no observed waning. The earlier Swedish study (Jiayao Lei et al.,
*NEJM* 2020, [DOI 10.1056/NEJMoa1917338](https://doi.org/10.1056/NEJMoa1917338))
already reported early-versus-later vaccination differences. The checked
main Methods/Results do **not** explicitly report a three-bin common-ratio
homogeneity/interaction test, so this small necessary-factorisation
compatibility question is the incremental descriptive exercise. This
literature survey does not establish absolute bibliographic novelty; if
a prior formal constancy test is identified, the only defensible remaining
contribution is the narrower printed-CI/rounding compatibility audit.

**THIS DOES NOT TEST WHETHER VACCINATION WORKS OR WHETHER PROTECTION WANES.**
Those subjects have already been evaluated in the vaccination literature;
Wu et al. report sustained association without indication of waning. The
present ratio is a ratio of adjusted observational comparisons, not a
head-to-head vaccine-age trial. Programme, birth cohort, prior HPV exposure,
attained age/calendrical support, selection and screening are confounded
with any biological interpretation. Nothing here justifies a vaccination,
screening, diagnosis or treatment recommendation.

## Source and browser/data boundary

Read the searchable article **HTML text** for the fully adjusted IRRs, model
adjustments, duration bins and printed CI endpoints; the HTML article footer
also independently establishes the exact CC BY 4.0 licence. The direct
HTML open attempt returned HTTP 403 in the available reader; indexed HTML
text remained readable. The figure's native rendered numerical labels were
not independently machine-readable here: the six case/person-year pairs
specified in the research queue and user brief are retained *only* as
attributed descriptive context, not substituted adjusted estimates.
No browser files were downloaded; no prompts were accepted, and no PDF,
supplement, raw registers or individual participants were ingested. This
access limitation prevents claiming an independent audit of the figure
image. Only source-text-selected aggregates with locators were stored.

## Offline reproduction

From the repository root (Python 3.12 and 3.14, stdlib only, no network):

    python -B research/hpv-age-duration-factorisation/analysis.py --check
    python -B research/hpv-age-duration-factorisation/analysis.py --json
    python -B -m unittest discover -s tests -p 'test_research_hpv_age_duration_factorisation.py' -v

The script recomputes the frozen result block of REPORT.md (without writing)
and normalised SHA-256 file hashes. Use `--write-report` only when deliberately
regenerating this **newly added** report. The exact fraction decisions are
independent of floating point; log allowances and approximate Bonferroni and
Poisson-tail roots are bounded numerical diagnostics. Passing CI is not
source verification or scientific review. See [DATA_POLICY.md](../../DATA_POLICY.md).
