# Source audit and factorisation compatibility report

**NONCLINICAL, UNREVIEWED RESEARCH. Not clinical advice.**

## Literature screen and attribution

Wu et al. (2026), *The BMJ* 392:e087326, DOI
[10.1136/bmj-2025-087326](https://doi.org/10.1136/bmj-2025-087326),
[article HTML](https://www.bmj.com/content/392/bmj-2025-087326.full),
Results "Invasive cervical cancer by time since vaccination," Fig. 3,
Methods "Statistical analyses"; checked 2026-11-18 (America/New_York).
The article's HTML licence notice says [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Only six aggregate adjusted IRRs with their printed marginal intervals and
six descriptive vaccinated case/person-year pairs were transcribed in the
[attributed sparse source file](data/wu_figure3.csv). The adjusted IRRs come
from the published fully adjusted Poisson model with age, calendar time,
county, maternal origin and history, education and household income. No
unvaccinated person-years, joint covariance, individual records, source
images, figure pixels or supplementary material were used.

The paper's stated findings on vaccination associations, age at initiation
and no detected waning, and Lei et al. (2020), DOI
[10.1056/NEJMoa1917338](https://doi.org/10.1056/NEJMoa1917338),
already cover those general topics. In the surveyed main Methods/Results,
no explicit three-bin common age-at-initiation-IRR ratio test was identified.
If one is subsequently found, this note is confined to the narrower
rounded-published-box compatibility exercise; no new biological discovery
or absolute novelty claim follows.

## Frozen assumptions, partial identification and limits

The [pre-computation protocol](PROTOCOL.md) was committed alone as
[b1401db](https://github.com/dovvnloading/open-cancer/commit/b1401db41cbca4b0fd015656c216fb1a38456b9f).
The exact fractional projection of each **marginal** CI pair supplies a
necessary compatible set, not a fitted joint likelihood, an inferential CI
for a ratio, or a biological age-duration effect. An empty box intersection
would falsify only that descriptive compatibility rule, not vaccination.
A nonempty intersection does not validate constant effects. The two-sided
Poisson count-model intervals in the generated block below characterize
**crude vaccinated rates only**, never adjusted IRRs. Missing or zero-event
adjusted age-bin comparisons remain untestable, not zero ratio estimates.

## Recomputed results (do not edit generated block by hand)

<!-- BEGIN RECOMPUTED RESULTS -->

### Primary exact-rational rounded-interval projection

Each ratio is the source's adjusted IRR for dose age 10-16 divided by its adjusted IRR for dose age 17+ **in the same duration bin**. Ratios below are projections of **marginal** intervals, not joint 95% confidence intervals or direct vaccinated-only comparisons.

| Years since first dose | Printed point ratio | Exact point fraction | Primary lower | Primary upper | Exact rounded lower / upper |
| --- | ---: | --- | ---: | ---: | --- |
| 7-9 | 0.311688311688 | 24/77 | 0.082251082251 | 1.135922330097 | 19/231 / 117/103 |
| 10-12 | 0.388888888889 | 7/18 | 0.132947976879 | 1.246153846154 | 23/173 / 81/65 |
| 13-15 | 1.000000000000 | 1/1 | 0.173553719008 | 6.200000000000 | 21/121 / 31/5 |

### Frozen intersection and sensitivity results

The two middle bounds are max(lower) and min(upper), not a valid ordered interval when max(lower) exceeds min(upper).

| Prespecified variant | Descriptive verdict | max lower | min upper | Sharp symmetric log allowance | Limiting lower / upper bins |
| --- | --- | ---: | ---: | ---: | --- |
| Primary; 0.005-rounded endpoints | **CONSISTENT** | 0.173553719008 | 1.135922330097 | 0.000000000000 | 13-15 / 7-9 |
| Printed CI; no widening | **CONSISTENT** | 0.183333333333 | 1.115384615385 | 0.000000000000 | 13-15 / 7-9 |
| Full 0.01-unit widening | **CONSISTENT** | 0.163934426230 | 1.156862745098 | 0.000000000000 | 13-15 / 7-9 |
| Rounded points only (not CI) | **FALSIFIED** | 0.957446808511 | 0.320261437908 | 0.547566255671 | 13-15 / 7-9 |
| Approximate six-IRR Bonferroni/Wald | **CONSISTENT** | 0.075008455197 | 1.996076523107 | 0.000000000000 | 10-12 / 7-9 |
| Primary minus 7-9 years | **CONSISTENT** | 0.173553719008 | 1.246153846154 | 0.000000000000 | 13-15 / 10-12 |
| Primary minus 10-12 years | **CONSISTENT** | 0.173553719008 | 1.135922330097 | 0.000000000000 | 13-15 / 7-9 |
| Primary minus 13-15 years | **CONSISTENT** | 0.132947976879 | 1.135922330097 | 0.000000000000 | 10-12 / 7-9 |

**Primary descriptive verdict: CONSISTENT.** Common q in [0.173553719008, 1.135922330097]. Exact intersection extrema: 21/121 through 117/103. Minimum symmetric log allowance epsilon*=0.000000000000; multiplicative allowance exp(epsilon*)=1.000000000000. Attaining log midpoint=-0.811911580653 (q=0.444008496420).

**Point-only sensitivity: FALSIFIED**, epsilon*=0.547566255671, exp(epsilon*)=1.729039852128, attaining log reference=-0.591051367611, q=0.553744789244; this is a point-rounding diagnostic and does not override the marginal-CI primary rule.

The approximate Bonferroni/Wald variant expands six source 95% marginal log-CI widths using a normal quantile at alpha/6. Because the actual IRR covariance and underlying interval method are unavailable, this is **not** an established 95% family-wise test. Statistical inference: **INCONCLUSIVE**.

### Descriptive exact-count Poisson identification intervals

These invert exact Poisson tails conditional on each **vaccinated** person-year denominator. They describe crude rates only; they never reconstruct an adjusted IRR or its standard error.

| Years | Dose age | Vaccinated cases / person-years | 95% exact-count crude rate per 100,000 PY |
| --- | --- | ---: | --- |
| 7-9 | 10-16 | 5 / 797583 | [0.203550776548, 1.462961482294] |
| 7-9 | 17+ | 26 / 264588 | [6.419060280737, 14.398243337992] |
| 10-12 | 10-16 | 11 / 668369 | [0.821576160360, 2.944786265267] |
| 10-12 | 17+ | 18 / 195405 | [5.459400107674, 14.558358418427] |
| 13-15 | 10-16 | 8 / 212315 | [1.626749017615, 7.424435023523] |
| 13-15 | 17+ | 4 / 87765 | [1.241799548369, 11.669331368317] |

For excluded early zero-case cells (1-3, 4-6, 16-18 years), a two-sided exact-count Poisson bound has mean lower 0 and upper -ln(0.025), conditional on known exposure. Their **adjusted age-ratio identification set is [0,+infinity] / UNTESTABLE**, because the early adjusted IRR and common unvaccinated adjustment/reference are missing. No continuity corrections or invented point ratios are used.

### Normalized reproducibility hashes

UTF-8 CRLF and CR are converted to LF before SHA-256; keys use Path.as_posix(). Hashes document text identity, not source truth.

| Input path | SHA-256 |
| --- | --- |
| research/hpv-age-duration-factorisation/PROTOCOL.md | bda8a6cf82cf54bf1c0bc7da700130c8a29d1b782b1943995817bf694ea3253a |
| research/hpv-age-duration-factorisation/data/wu_figure3.csv | 459d3b99f33240bc1090a2a73f64d66963f7ddcee7be354772e03d405ef798f5 |
| research/hpv-age-duration-factorisation/analysis.py | 6756bb424f1d88394e5b340b967af9147733cccb15d8c70a3671103bf6ad36c4 |

<!-- END RECOMPUTED RESULTS -->

## Hostile self-review: attempted falsification of the interpretation

1. **False statistical certainty.** Six marginal 95% intervals are not a
   simultaneous 95% confidence region for three ratio contrasts. Their
   dependence, cross-row regression covariance, CI construction, and
   populations' shared reference effects are unknown. The primary result
   is a deterministic interval-box audit only. Even the prespecified
   normal-approximate Bonferroni variant cannot supply a calibrated p-value
   or a definitive interaction rejection/nonrejection.
2. **Point estimates can mislead.** The point-only rounded ratio bands
   fail to intersect, with a positive sharp allowance, whereas the
   uncertainty-aware primary bands intersect. Reporting the point-only
   outcome as an inferential rejection would contradict the protocol.
3. **Sparse denominator trap.** The late 13-15-year group contains only four
   reported cases; early 7-9 contains five. The six vaccinated person-year
   denominators are not the unvaccinated reference and cannot reconstruct
   fully adjusted IRRs. Exact-count *crude* intervals are separate, not a
   second confidence sample used to update adjusted source CIs.
4. **Zero-event trap.** In 1-3, 4-6 and 16-18 years, the early age group
   lacks a fitted adjusted IRR. Treating the figure dash or zero cases as
   an IRR of zero, or adding 0.5 cases to manufacture one, would silently
   expand the panel and introduce unvalidated adjustment assumptions.
   Their adjusted age ratio is **UNTESTABLE**, even with a Poisson crude
   zero-event upper rate bound.
5. **Noncomparable biology.** The two vaccination-age groups have different
   attainable ages, historical vaccination programmes, birth cohorts,
   screening opportunities, prior HPV exposure and potential risk behaviour.
   Covariate adjustment does not make them exchangeable. A constant ratio
   here is not a biological absence of age-duration interaction; a changing
   ratio would not establish a cellular, immunity, waning or causal effect.
6. **Source and legal boundary.** Live direct opening of BMJ's HTML returned
   403; its indexed HTML text exposed the adjusted results and the article
   licence statement. The figure image's numeric labels could not be
   independently rechecked in this environment; the six queued/context
   pairs are disclosed, attributed but not independently figure-audited.
   No browser download prompt was accepted; no PDFs or source archives were
   downloaded. This is a reproducibility limitation, not licence permission
   to mirror registers. Article CC BY 4.0 does not waive register privacy.
7. **Software and novelty boundary.** A normal-tail approximation is not an
   exact simultaneous interaction test; exact fractions and count-tail
   inversion have their declared input assumptions. All sensitivity modes
   were frozen before computation, and no bin was dropped post hoc. The
   earlier Swedish report and Wu et al. already cover efficacy associations,
   age differences and no observed waning. There is no clinical conclusion,
   no independent replication and no claim of bibliographic priority.

**Final interpretation:** A shared early-versus-later adjusted-IRR multiplier
is **CONSISTENT** with the selected published marginal interval boxes;
statistical evidence for or against genuine age-by-duration effect
modification is **INCONCLUSIVE**. This exercise **does not test whether HPV
vaccination works or wanes**; those issues have already been studied, with
no indication of waning reported in this source. Programme, birth cohort,
prior HPV exposure and screening are confounded with any biological reading.
No vaccination, screening or clinical advice is provided.
