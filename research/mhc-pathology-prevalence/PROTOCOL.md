# Frozen protocol: MHC-I proxy and pathology prevalence

**Nonclinical, unreviewed computational research. No medical or veterinary advice.**

## Status, question, and scope

Freeze this protocol and the derived common-panel files in a separate commit
before fitting any model or computing predictive scores. Inspection of source
licences, schema, names, missingness, endpoint arithmetic, and descriptive
panel counts is allowed before that commit. The published MHC association and
PR #79's inconclusive mortality result are already known. This is a retrospective
evaluation, not blinded testing, independent preregistration, or a prospective
validation study. Preserve this text after fitting; any correction must be a
dated, explicit amendment in a separate file and must not overwrite the freeze.

Base repository: `dovvnloading/open-cancer`, current main at task start
`620c1df128ce3a30dbaaa190179cb5404048174a`. This implements ranked job 2 in
[the pinned research index](https://github.com/dovvnloading/open-cancer/blob/620c1df128ce3a30dbaaa190179cb5404048174a/research/README.md).
The governing data policy is [DATA_POLICY.md at that revision](https://github.com/dovvnloading/open-cancer/blob/620c1df128ce3a30dbaaa190179cb5404048174a/DATA_POLICY.md).

Question: does one `log(1 + MHC-I gene-family size)` term improve prediction
of denominator-linked neoplasia prevalence in a wholly held-out order beyond
body mass, adult life expectancy, and diet? We transport PR #79's **learning
rules**, refitting all coefficients on the new endpoint; we do not reuse its
mortality-fitted coefficients. The MHC source covers mammals, so the realized
test concerns the matched mammalian subset of a vertebrate pathology dataset,
not transport across all vertebrate classes.

## Sources and licence audit before reading data

Access/reuse audit: **2026-11-18 (America/New_York)**. Both repository licences
and both final-article licences were checked before reading either source CSV.

| Source | Immutable revision / publication | Licence finding |
|---|---|---|
| Vincze et al., *Immunological surveillance against cancer across mammals* | [DOI 10.1038/s41467-025-65286-x](https://doi.org/10.1038/s41467-025-65286-x), [PMC12644480](https://pmc.ncbi.nlm.nih.gov/articles/PMC12644480/) | Final article: CC BY-NC-ND 4.0. |
| Vincze repository | `OrsolyaVincze/VinczeEtalNatComm_MHC-cancer`, commit `2c8a5c5b5fda6fb34c86b0e5a2b6a58dcc85cc3e` | Complete pinned tree has no licence file or README; current GitHub repository metadata reports no licence. |
| Compton et al., *Cancer Prevalence across Vertebrates* | [DOI 10.1158/2159-8290.CD-24-0573](https://doi.org/10.1158/2159-8290.CD-24-0573), [PMC11726020](https://pmc.ncbi.nlm.nih.gov/articles/PMC11726020/) | Final article: CC BY-NC-ND 4.0. Do not substitute the preprint's different licence. |
| Compton repository | `zacharycompton/cancerAcrossVertebrates`, commit `4819e0933358c5fc7eab483a0c05161ea8c074b7` | Complete pinned tree has no licence file; its README gives no redistribution grant; current GitHub metadata reports no licence. |

Audit locators: [Vincze pinned tree](https://api.github.com/repos/OrsolyaVincze/VinczeEtalNatComm_MHC-cancer/git/trees/2c8a5c5b5fda6fb34c86b0e5a2b6a58dcc85cc3e?recursive=1),
[Compton pinned tree](https://api.github.com/repos/zacharycompton/cancerAcrossVertebrates/git/trees/4819e0933358c5fc7eab483a0c05161ea8c074b7?recursive=1),
[Compton pinned README](https://github.com/zacharycompton/cancerAcrossVertebrates/blob/4819e0933358c5fc7eab483a0c05161ea8c074b7/README.md),
and [CC BY-NC-ND 4.0 terms](https://creativecommons.org/licenses/by-nc-nd/4.0/).
Article access included official indexed article text when the direct PMC
reader returned a bot-check page; that limitation did not affect the GitHub
plain-text source reads. No access control was bypassed.

Public availability and article access do not establish a right to mirror the
repository assets. Apply the user's explicit aggregate-only fallback, as in
PR #79: commit original code/documentation, narrowly derived species-level
numerical facts needed to reproduce the report, reconciliation metadata, and
pinned fetch-and-verify instructions with attribution. Do not commit either
upstream CSV, the dictionary, upstream R code, article prose, or individual
records. Do not represent the articles' ND licence as permission to distribute
adapted licensed material. The project MIT licence does not relicense any
third-party material. This is a conservative factual-aggregate route, not a
claim that either upstream repository supplies an unrestricted reuse grant.

Read only these data inputs, as UTF-8 text through GitHub API/raw endpoints:

- [Vincze `SourceData_FINAL.csv`](https://github.com/OrsolyaVincze/VinczeEtalNatComm_MHC-cancer/blob/2c8a5c5b5fda6fb34c86b0e5a2b6a58dcc85cc3e/SourceData_FINAL.csv).
- [Compton `species-cancer-prevalence-data.csv`](https://github.com/zacharycompton/cancerAcrossVertebrates/blob/4819e0933358c5fc7eab483a0c05161ea8c074b7/species-cancer-prevalence-data.csv).
- [Compton data dictionary](https://github.com/zacharycompton/cancerAcrossVertebrates/blob/4819e0933358c5fc7eab483a0c05161ea8c074b7/data-description-species-cancer-prevalence-data.txt).

No browser file download or download prompt was used. Any future browser
download prompt must be denied and the limitation recorded; it is not a
reason to retrieve a file by browser. The optional fetch script reads pinned
response text into memory and writes only derived outputs. Restricted
Species360 records and individual animal records are not requested.

## Species reconciliation and the frozen common panel

For each source species name, replace underscores with spaces, strip leading
and trailing whitespace, collapse internal whitespace to one space, and apply
Unicode `casefold()`. Join only exact equality of these normalized strings.
Do not fuzzy-match, consult additional taxonomies, map synonyms, truncate
subspecies, or combine different names. Duplicate normalized names in either
source fail derivation rather than trigger aggregation or a hand-picked match.
Store both source names and every exclusion in `data/reconciliation.csv`.
An unmatched name is a limitation of this strict panel, not proof of distinct
biological species.

Use Vincze's `Taxonomic order` consistently for every fold. In particular,
Compton's missing order for `Antilocapra americana` does not exclude it or
create an artificial group. A disagreement in supplied order metadata is
recorded; it never changes this source rule after fitting.

Eligibility requires a matched name; nonmissing finite positive Vincze
`Body mass (kg)` and `Adult life expectancy`; binary `Consuming mammal prey`;
nonnegative integer `MHC I gene family size`; a nonempty Vincze order; and
Compton integer `RecordsWithDenominators > 0` with integer
`0 <= NeoplasiaWithDenominators <= RecordsWithDenominators`. Invalid observed
values stop derivation; missing required values produce explicit exclusions.
Both models always use the same rows. Retain all eligible zero-neoplasm rows.
Do not require Vincze CMR, RDI, ICM, any mortality outcome, or secondary-endpoint
availability. Do not substitute Compton maximum longevity for adult life
expectancy. Do not exclude a species for an inconvenient outcome or fit.

The prefit audit gives 139 Vincze and 327 Compton species, 47 exact normalized
matches, and **42 eligible species in 9 orders**, including **2 zero-neoplasm
species**, with **312 recorded neoplasms among 2215 denominator-linked records**.
The sorted species identifiers and values in the freeze commit's
`data/species_aggregates.csv` define the common panel. These are descriptive
data audits, not model results. All included Compton classes are Mammalia.
The strict join retains `Heterocephalus_glaber`; `Dolichotis_patagonum` has no
strict match and is not manufactured or substituted.

## Endpoints and the malignancy dictionary conflict

Primary: `k = NeoplasiaWithDenominators`, `n = RecordsWithDenominators`.
Use the observed integer counts, not a reconstruction from a rounded
prevalence column, and not counts from `RecordsFromAllSources` or
`NeoplasiaFromAllSources`. Neoplasia includes benign and malignant findings.
The ratio is a recorded pathology prevalence among available denominator-linked
records, not onset incidence, lifetime risk, or cause-specific death.

Before fitting, audit the printed neoplasia and malignancy prevalence columns
against candidate count ratios. For a nonendpoint decimal token with exponent
e, accept absolute discrepancy at most `0.5 * 10**e + 2e-15`. Literal 0 and 1
are exact endpoint proportions with only the fixed `2e-15` export allowance;
integer display does not create a half-probability tolerance. This respects
different displayed precisions. Store derived residuals and tolerances per
species in `data/endpoint_audit.csv`, including rows outside the common panel.

The dictionary describes `MalignancyKnown` as the numerator of
`MalignancyPrevalence`, but the pinned numerical column is reproduced by
**`Malignant / RecordsWithDenominators`**. The alternative dictionary ratio
fails for many species, including `Carollia perspicillata` in the common panel.
The source [bootstrap script](https://github.com/zacharycompton/cancerAcrossVertebrates/blob/4819e0933358c5fc7eab483a0c05161ea8c074b7/bootstrap.R)
also sums `Malignant` indicators when constructing its malignant endpoint.
This resolves which count reproduces the published endpoint; it does not
reconstruct undocumented pooling of records or classify every unknown tumour.
`MalignancyKnown` can exceed denominator-qualified neoplasia and must not be
silently relabelled as classified neoplasms in the same cohort.

Secondary, only after this arithmetic audit passes: recorded malignant
pathology, `k = Malignant`, with the same `n` and the same primary panel.
Require integer `0 <= Malignant <= NeoplasiaWithDenominators` and agreement
with the printed malignancy prevalence under the frozen tolerance. If any
included row fails, abstain from the secondary analysis and report the failure;
do not alter the primary panel. The observed panel passes these prefit gates.
Unknown malignant status is not known benign status. This secondary endpoint
does not estimate latent true malignancy, and cannot rescue the primary verdict.
Do not use `PropMalignant` anywhere in fitting or scoring, and **never interpret
it as a progression probability**. Its denominator/provenance conflicts remain
unresolved.

## Paired learning rules and numerical bounds

Both models are binomial ridge logistic regressions with an unpenalized
intercept. Baseline slopes: natural log of body mass (kg), natural log of
adult life expectancy (days), and mammal-prey diet (0/1). The extension adds
exactly one natural-log `log1p(MHC-I family size)` slope. Adult life expectancy
is mean remaining life after sexual maturity, not maximum lifespan or life
from birth. Gene-family size is an inferred genomic proxy, not expressed,
functional, or peptide-binding repertoire size.

In each training fold, standardize every continuous predictor using the
unweighted training-species mean and population standard deviation. Give a
zero-variance column scale 1. Leave binary diet unstandardized. Apply those
training transformations to held-out species. No global scaling, interactions,
polynomials, taxonomic indicators, variable selection, or held-out tuning.

Minimize

    sum_i [n_i * softplus(eta_i) - k_i * eta_i]
        + (lambda / 2) * sum_{j > 0} beta_j**2

with fixed `lambda = 1`. Diet is penalized like all other slopes; the intercept
is not. The primary data likelihood is count-weighted, although the primary
test score below weights orders equally. A half-event intercept initialization
does not change counts, likelihood, or fitted outcomes.

Use the stable damped Newton/Cholesky algorithm from PR #79, with at most
100 accepted steps, up to 45 line-search candidates per step, Armijo factor
`1e-4`, and objective roundoff allowance `1e-12 * max(1, abs(objective))`.
Convergence requires maximum absolute penalized gradient no larger than
`max(1e-8, 1e-10 * effective_training_denominator)`. These are declared numerical
bounds, not post-hoc values selected to accommodate an observed fit. A fit
with no finite intercept optimum, invalid Hessian, failed line search, or
unmet gradient bound fails explicitly. Never omit that fold. Store diagnostics;
test their bounds rather than exact saved iteration or gradient values.

## Complete held-out-order scoring

For each of the 9 Vincze orders, fit both models on all other orders and predict
every species in the held-out order. The held-out order contributes no outcomes
to its training, scaling, or learning-rule choices. Its covariates are used only
for its predictions. Preserve membership, training transformations, fitted
coefficients, diagnostics, and predictions in the structured results.

Species log loss is

    L_i = -[k_i * log(p_i) + (n_i-k_i) * log(1-p_i)] / n_i.

It includes zero-event species: when k=0 the loss is `-log(1-p)`, not zero.
For numerical scoring only, clip probabilities to `[1e-12, 1-1e-12]` and report
the number clipped. Binomial combinatorial constants cancel between models.
Average losses equally across species within each order, then average the
order means equally. The primary gain is

    delta = equal_order_baseline_loss - equal_order_extended_loss.

Positive delta favors adding MHC-I. Units are natural-log cross-entropy per
denominator-linked record, with the stated species/order averaging. Also
report record-weighted losses and gain separately, each order's counts and
losses, the number of improving orders, and every held-out prediction.
Record-weighted gain cannot replace the primary estimand.

## Order-level uncertainty

Delete each order in turn and rerun the **entire** held-out-order evaluation
on the remaining orders, including both model fits and training transforms
for every remaining holdout. Deleting fixed prediction scores is insufficient.
For G=9, full-panel delta, and complete-refit `delta_(-g)`, calculate

    q_g = G * delta - (G-1) * delta_(-g)
    center = mean(q_g)
    SE = sqrt(sum_g (q_g-center)**2 / [G*(G-1)])
    interval = center +/- 2.306004135204166 * SE.

The critical value is the two-sided 95% Student t quantile with 8 degrees of
freedom, not PR #79's t15. Report the unadjusted full-panel gain separately
from the bias-corrected jackknife center; the interval is centered on the
latter. Report every deletion result, its pseudovalue, and the deletion range.
Repeat the same procedure for the same-panel secondary endpoint, explicitly
exploratory with no multiplicity-adjusted confirmatory claim.

These intervals are approximate for a small, selected, uneven and
phylogenetically related set of orders. They are not independent-animal
uncertainty, a source-ascertainment correction, a phylogenetic covariance
adjustment, or a guarantee of population coverage.

## Prespecified falsification attempts

Report every attempt. Never select a favorable sensitivity after seeing gains.
Keep the primary held-out-order scoring rule and refit both models for:

1. **Order influence:** all complete delete-order refits above, including any
   sign reversal or dependence on a singleton order.
2. **Equal-species training:** give every training species effective denominator
   equal to that training fold's mean n, preserving its k/n and total training
   likelihood scale; retain lambda=1 and original test counts.
3. **Penalty:** separately use lambda=0.1 and lambda=10, with no selection.
4. **Named rodents:** remove `Heterocephalus_glaber` and
   `Dolichotis_patagonum` if present. The latter is already absent under the
   frozen join; record which removal actually occurs.
5. **Recorded zeros:** remove k=0 species as an explicitly outcome-selected
   stress test. This changes the target and can neither rescue nor veto the
   primary claim.
6. **Secondary endpoint:** the dictionary-resolved malignant count, with all
   primary-panel species retained, complete refits, and its own order jackknife.
   It does not determine the primary verdict.

The four major robustness checks for a positive primary verdict are
equal-species training, the two penalties, and named-rodent removal. Their
failures or failed fits block a positive claim. The secondary endpoint and
outcome-selected zero removal do not. No new tree source is added: full
phylogenetic covariance adjustment, alternative trees, alternative synonym
maps, and independently ascertained populations remain untested. Mortality
count-rounding sensitivities from PR #79 do not apply to these integer counts.

## Frozen verdict and failure observable

First require the frozen panel/count/hash gates, at least 8 complete orders,
and successful primary and all complete-refit uncertainty fits. Otherwise
report **inconclusive (not evaluable under the protocol)** with the reason.
Never call a failed computation evidence of no association.

For an evaluable primary analysis:

- **MHC-I transports to pathology prevalence within this evaluation** only
  if delta>0, the jackknife lower bound>0, strictly more than half the held-out
  orders improve, and all four major sensitivity gains are positive with
  successful fits.
- **MHC-I does not transport within this evaluation** if delta<=0 and the
  jackknife upper bound is strictly below zero. This contradicts improvement
  for these learning rules, this panel, and this scoring procedure only.
- **Inconclusive** otherwise, including any interval containing zero, a zero
  interval boundary, conflicting point/interval direction, or failure of a
  required positive-claim sensitivity. An interval spanning zero is not
  evidence of equivalence or of no MHC effect.

Independently report the **point-gain failure observable, delta<=0**. A
nonpositive point gain fails improvement but does not by itself establish the
directional non-improvement verdict. Use strict upper<0 for contradiction;
upper=0 is not below zero. A lower>0 interval cannot override a nonpositive
full-panel gain. No secondary finding can change the primary verdict.

## Interpretation, novelty, and limits

This tests **endpoint portability of a measured gene-family-size proxy after
PR #79**, not PR #59's functional-load theorem, expressed antigen repertoire,
immune surveillance as a mechanism, or complete presentation loss. Functional
antigen load, peptide quality, antigen presentation, and immune killing are
not measured. The prior MHC association is credited to Vincze et al.; pathology
prevalence comparisons are credited to Compton et al.; the held-order audit
design is inherited from PR #79. The bounded contribution is the pinned join
and frozen endpoint-prediction audit. Do not claim historical priority, a new
immune mechanism, or independent biological validation. A limited literature
check cannot establish novelty by absence.

Pathology prevalence is neither onset risk nor cause-specific death.
Ascertainment, differing necropsy practices, ages, sampling, tissue coverage,
annotation quality, and record denominators can explain associations.
Shared zoos and possible animal overlap with PR #79 remain limitations;
species overlap is not evidence for or against individual-animal overlap.
The new endpoint does not make observations independent. No result establishes
immune causation, applicability to wild populations, or medical/veterinary
utility. Negative and inconclusive results do not refute immune surveillance.

## Reproduction, hostile review, and immutable scope

Use only the Python standard library. The dedicated test must run offline
in under 60 seconds, recompute every reported numerical analysis from committed
data, verify zero retention and reconciliation/audit arithmetic, check training
membership and label/covariate isolation, independently recalculate gradients,
probabilities, scores and jackknife formulas, and exercise verdict boundaries.
Patch network entry points during the analysis rebuild. Compare scientific
floating results with declared absolute/relative tolerance 3e-7; require exact
integers, labels, structural membership and hashes. Check optimizer diagnostics
against their declared bounds in both fresh and archived results rather than
exact saved stopping values. Check the README numerical block against rendering
the archived results, while comparing fresh scientific results numerically;
do not demand identical platform-specific convergence traces.

Hash UTF-8 file content after replacing CRLF and CR with LF, without trimming
whitespace or changing final-newline presence. Build **every file-path key**
with `Path.as_posix()`, including source and derived-file provenance maps.
Exercise LF/CRLF/CR and Windows-style path behavior. Python 3.12 and 3.14,
Windows and Linux must use the same bounded numerical contract. Record which
execution environments were actually available; portability tests do not
claim an unrun CI matrix.

Only add `research/mhc-pathology-prevalence/` and
`tests/test_research_mhc_pathology_prevalence.py`. Modify no existing file,
other research directory, index, workflow or data policy. Before pushing,
perform hostile review of licensing, joins, endpoint semantics, numerical
claims, leakage, uncertainty, reproducibility and novelty; record corrections,
failed attempts and remaining limitations. Preserve the prefit commit as an
ancestor of results. Push the branch, open the requested PR, and inspect its
actual CI. Do not merge. Passing tests and AI review establish computational
consistency only, not source truth, biological validity or clinical utility.
