# Frozen protocol: somatic mutation rates and Peto's paradox

**Nonclinical, unreviewed computational research. No medical advice.**

Protocol date: **2026-11-18, America/New_York**. This file will be committed
and pushed on a branch of current `main` before any empirical regression,
cross-validation, model scoring, inverse-rate calculation, or residual-bound
calculation. The README will link the public freeze commit. This is a
prospective analysis specification for already published, partly familiar
data, not blinded research or independent trial registration. Source reading,
rights checks, literal transcription and checking column definitions precede
the freeze. The protocol will not be silently rewritten after fitting.

## 1. Literature-first question and pivot

Read first: `research/README.md`, the full model/observation discussions in
`research/multistage-incidence-identifiability/` (PR #64),
`research/stem-cell-divisions-risk/` (PR #70), `AGENTS.md`, and `DATA_POLICY.md`.

The classical multistage scaling is not new. Peto (1975/1977), Caulin and
Maley (2011), Nunney (2020), and Cagan et al. (2022), especially Supplementary
Note 2, already discuss rate reduction, cell numbers, and additional genetic
controls. Compton et al. (2023 preprint/2025 Cancer Discovery) already joins
mutation rates with cancer prevalence. Zheng (2025 JITC) explicitly joins
Cagan and Vincze and claims quantitative resolution. Butler et al. (2025),
the ensuing PNAS exchange, and Dujon et al. (2025) demonstrate that endpoint,
effect magnitude, phylogeny, and sampling uncertainty matter. Firsanov et al.
(2025) supplies counterevidence to a universal extra-oncogenic-hits mechanism.
Recent 2026 work will be described without converting mechanistic hypotheses
or germline conservation into measured somatic transition rates.

**Pivot:** audit the claimed resolution against a specified cell-aware model,
derive a sharp two-species finite-risk compatibility envelope without fixing
the absolute target-cell count or passenger-to-driver conversion, and confront
that envelope and held-out prediction with eligible published aggregates.
The candidate contribution is this composed certificate and its transparent
application, not first discovery of inverse-lifespan scaling, rate/step
tradeoffs, or nonidentifiability. Bibliographic priority is unestablished.

## 2. Permitted sources, fields, and missingness

Use only small published per-species aggregates in article text or openly
readable tables. An article-linked GitHub CSV rendered as a table qualifies.
Do not use browser file downloads, binary spreadsheet/archive extraction,
digitized point positions, raw genomes, or individual human/animal records.
If a download prompt or inaccessible supplement blocks a source, proceed
without it and document the limitation. No access-control bypasses.

Eligible planned sources:

1. Cagan 2022 annual SBS/genome mutation rates, as printed in its text or
   reproduced in the openly readable `Mutation_Data.csv` linked by Compton's
   CC BY 4.0 2023 preprint. Pin the table's 2022 creation commit and verify
   its Cagan columns against the later table version. Retain mean rates and
   zero-intercept slopes separately. A confidence interval for a regression
   slope is **not** an interval for the mean-rate estimator. Preserve negative
   printed lower limits; do not silently clip them to a small positive value.
2. Vincze 2022 article-linked `data.csv`: CMR, neoplasia-attributed deaths,
   known-cause deaths, total deaths, total animals, female/male adult mass,
   adult remaining life expectancy, and ICM where reported. Pin its public
   publication commit. CMR's denominator is `knownDeaths`, not `Sample` or
   all `Deaths`. ICM is a different, censoring-adjusted endpoint.
3. The same Compton preprint-linked table's published maximum longevity and
   adult mass for the nine mutation-overlap species, provided attribution and
   reuse permit. These life traits are not Cagan's 80%-mortality lifespan.
4. Cagan's printed giraffe/naked-mole-rat lifespan/rate comparison as a
   separately labelled, rounded-input sensitivity, if its numbers are readable.

Every numeric-data row must have source citation, source URL, table name,
exact row/column locator, access date, reuse basis and terms URL. Retain source
versions and original units. Record absent intervals and sample sizes as
missing, never zero, estimated from a figure, or invented. Species-level
mutation sample sizes may remain missing if only inaccessible supplements
provide them. Do not reproduce Compton's cancer outcome data under an assumed
license from its more restrictive final article. Original Cagan/Vincze
CC BY data and the separately documented preprint route retain attribution;
the repository's MIT license does not relicense third-party data.

No species is removed because it is an outlier. The planned matched panel is
the intersection of the nine accessible mutation rows with Vincze's table;
six overlaps are known before fitting. Unavailable original Cagan rows and
incompatible taxa are recorded. No substitute species or breed is inserted.
Species-level joins collapse source subspecies/common labels only to the
explicit binomial printed by the source table; disclose that loss of detail.

## 3. Models, target, and endpoints

Let `m` denote measured genome-wide passenger SBS/year and `mu=c*m` a
postulated functional transition intensity with unknown `c>0`. They are not
the same measured quantity. For common integer `k>=1`, independent initially
clean targets, no growth/loss/competition, equal exponential stage rates and
immediate observation, use

`R=1-exp(-H)`, `H=N*psi_k(mu*L)`,
`psi_k(x)=x-log(sum(j=0,...,k-1) x^j/j!)`.

Prove the rare-lineage limit, the per-division clock version (including the
distinction from literal discrete Bernoulli trials), the required mass
exponent, endpoint equivalences, and the observation obstruction from
competing death/detection. Distinguish an equal-rate modeled transition from
an independently measured biological bottleneck. General expanding models
from PR #64 are not silently replaced by an expected malignant-cell count.

Primary **proxy** panel: `R=CMR`, `M=(female_mass+male_mass)/2` kg,
`L=adult_remaining_days/365.25`, and `m=published mean annual SBS`.
Setting a random remaining lifetime to its mean and equating a death fraction
with fixed-horizon onset are explicitly unverified bridges. The clean-start
adult clock is a mathematical surrogate, not a claim that adults have no
prior mutations. Secondary horizon: published maximum longevity/12 years.
Secondary mutation estimator: zero-intercept regression slope. These are
source-defined sensitivity versions, not interchangeable estimates.

The biological target is whether rates/lifespans alone explain cancer under
known target scaling. The proxy observations cannot by themselves identify
that target. A conditional failure of the shared model must not be described
as proof that large animals need additional biological steps.

## 4. Frozen fitting and prediction rules

Use natural logs internally; report log10-hazard RMSE. For each panel with
at least four positive finite CMRs, fit the following rare-lineage families,
with equal species weight and `B=m*L`:

| Family | log hazard specification; fitted quantities |
|---|---|
| Intercept comparator | `a`; fit residual mean |
| Mass comparator | `a+log(M)`; fit residual mean |
| Burden-only comparator, each k=1,...,8 | `a+k*log(B)`; fit residual mean |
| Cell-aware model, each k=1,...,8 | `a+log(M)+k*log(B)`; fit residual mean |
| Free mass exponent, each k=1,...,8 | `a+b*log(M)+k*log(B)`; OLS a,b |

The common `1/k!` factor is absorbed in each independently fitted amplitude.
The burden-only comparator omits the target-cell factor and is descriptive;
better prediction by it does not prove mutation scaling alone suffices.
The free mass exponent represents conditional aggregate residual scaling,
not a measured target-cell exponent or an evolutionary causal effect.

Report **every fixed k**, not a favorable post-hoc index. Refit within each
leave-one-species-out (LOSO) training set. Primary score: equal-species RMSE
of log10 H. Also report mean per-known-death Bernoulli negative log score,
averaging species equally; this is conditional on independent binomial death
labels, not an established sampling design. Compute extreme probabilities
stably, without arbitrary clipping. Report held-out risk for every species.
No held-out outcome enters fitting or model selection. No extra selected-k
score is planned. A rank-deficient fit is unavailable, not repaired ad hoc.

Also repeat the primary-panel comparison leaving out each taxonomic order
(four planned orders; both primates and both Panthera species held together).
This is a relatedness sensitivity, not a phylogenetic model or independent
evolutionary replication. State the small, selected sample explicitly.

For the nine-species life-trait table, descriptively fit log m on log maximum
longevity with intercept and compare LOSO predictions with a fixed inverse
slope (-1) and an intercept-only comparator. Do not call this replication of
Cagan's original 16-species, 80%-mortality-lifespan analysis. Report the free
slope's full leave-one-species-out range; no population p-value.

## 5. Uncertainty and exact compatibility

Use published intervals when they match their estimator. Where per-species
CMR intervals are absent, compute Wilson 95% intervals from the published
neoplasia/known-death counts with z=1.959963984540054. Label these as **derived**
binomial intervals, conditional on exchangeable independent known-cause
deaths. Missing causes, zoo sampling, species dependence and mortality/onset
mismatch are not repaired by these intervals. Do not call marginal interval
boxes a simultaneous confidence region or a calibrated hypothesis test.

For each fixed k report the minimum symmetric multiplicative hazard allowance
needed to fit the cell-aware rare family: for log residual intervals [l_i,u_i],
`rho=exp(max(0,max(l_i)-min(u_i))/2)`. Give both point and CMR-interval versions,
the constraining species, and whether rho=1 is feasible. This is a deterministic
compatibility certificate under stated boxes, not a significance test.

Derive and apply the **exact** pair envelope
`H_i/H_j in n*[min(b,b^k),max(b,b^k)]`, with `n=M_i/M_j`, `b=B_i/B_j`.
For k>1,b!=1 the finite-parameter feasible set has open endpoints; closure
endpoints are infima. Prove two-species sufficiency in the positive effective-
target extension; do not assert pairwise feasibility is sufficient for an
entire multi-species panel. Enumerate all pairs but predesignate
giraffe/naked mole rat (already highlighted by Cagan) as the main diagnostic.
Report the residual hazard multiplier interval implied by CMR intervals and
the exact envelope, so the amount of conditional additional suppression is
quantified without guessing an absolute driver rate or target count.

Sensitivity boxes: (a) both sexes' printed mass endpoints instead of their
mean; (b) mean annual rates replaced by published slope estimates; (c) all
available positive slope confidence limits applied only to those slopes,
with other inputs fixed and missing uncertainty explicitly unresolved;
(d) hypothetical independent factor-two boxes for each m and L, and both-sex
mass endpoints. Scenario (d) is a stress test, not a confidence interval.
If a slope interval reaches zero, a log lower bound is unbounded; report it.
Use monotonic extrema in the pair envelope and rare residuals, not random
corners presented as a proof. No biological extra-step count will be inferred
without an absolute effective driver clock and independent target definition.

For uncertainty in score comparisons, report all fold errors and the full
range of paired mean-score differences after deleting each species' scored
fold. This is a finite-sample sensitivity range, not a confidence interval.
No bootstrap that treats the six selected species as an evolutionary random
sample and no fabricated mutation-rate standard errors.

## 6. Verdict, falsification and review

Use an explicit two-level verdict:

* Mutation-rate scaling alone **suffices** only conditional on a specified
  target/driver/observation map that is compatible; predictive association
  alone cannot establish this biological verdict.
* **Additional suppression required** means the specified common-k/common-c,
  target-proportional model fails the deterministic compatibility criterion;
  quantify the necessary multiplier and identify the assumptions. It does
  not uniquely identify extra steps, apoptosis, repair, immunity or target
  counts. Source-data structural uncertainty can weaken that conditional claim.
* **Not identifiable** is the biological verdict when these mechanisms or
  observation maps remain observationally equivalent or the required matching
  information is absent. An untestable mechanism is not confirmed.

Name the largest positive and negative primary cell-aware k=3 LOSO residuals
(a fixed reporting index, not selected after fitting), and their ratios.
Retain all species. Give a prospective, falsifiable pair prediction conditional
on independently matched targets, driver clocks, tissue and age window; do not
offer clinical or experimental treatment advice.

Record attempts to defeat the conclusion: prior-art duplication, finite-risk
saturation, mortality/incidence conflation, different horizons, rate estimators,
reported versus missing uncertainty, related-species leakage, unfavorable
held-out scores, unbounded calibration/step aliases, and bowhead counterevidence.
Downgrade every unproved novelty or biological interpretation before pushing
the final analysis. AI critique and CI are not independent scientific review.

## 7. Reproducibility and delivery

Only add files within `research/petos-paradox-mutation-rates/` and one new
test file `tests/test_research_petos_paradox_mutation_rates.py`. Do not change
any existing file, workflow, research index or other research directory.
Use standard-library Python only. Recompute the complete reported numerical
result tree and README numerical block from CSVs. Hash UTF-8 content after
CRLF and lone CR become LF; construct every relative path key using
`Path.as_posix()` for Windows and Linux. Independent numerical checks must
exercise stable exact-chain formulas, scale invariance, pair bounds, interval
directions, fitting, held-out leakage and source denominators. The new test
file must run offline in under 60 seconds. Run repository validation and
existing tests as required. Push the completed branch, open a PR to `main`,
observe its actual CI, and report the PR number and CI status.
