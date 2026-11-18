# Does the MHC-I proxy transport to pathology prevalence?

**Nonclinical, unreviewed computational research. No medical or veterinary advice.**

**Verdict: inconclusive.** Adding MHC-I gene-family size worsens the primary
held-out-order point score, so the frozen point-improvement observable fails.
The approximate order-level interval contains zero, so the stricter directional
non-improvement verdict is not established. The major prespecified sensitivity
checks also fail to produce positive gains. This analysis therefore supplies
no positive evidence of portable incremental prediction on the evaluated panel,
while leaving the directional claim unresolved under its frozen uncertainty
procedure.

This implements ranked job 2 in [the research index at the starting revision](https://github.com/dovvnloading/open-cancer/blob/620c1df128ce3a30dbaaa190179cb5404048174a/research/README.md).
The endpoint is denominator-linked pathology prevalence, distinct from the
attributed cancer death tested inconclusively in [PR #79](https://github.com/dovvnloading/open-cancer/pull/79).
The actual overlap is mammalian: using a vertebrate source does not extend the
result to birds, reptiles, amphibians, or fish.

## Prefit freeze and what was tested

The [protocol](PROTOCOL.md), common-panel aggregates, name reconciliation,
endpoint arithmetic audit, licence record, and fetch script were committed and
published **before any model fitting** in
[freeze commit `196f997`](https://github.com/dovvnloading/open-cancer/commit/196f9970a80381bcf2cf0e2603ee33f843045e7c).
Its parent is current main at task start,
[`620c1df`](https://github.com/dovvnloading/open-cancer/commit/620c1df128ce3a30dbaaa190179cb5404048174a).
The published association, the earlier mortality result, and the source data
definitions were already visible. This was a retrospective protocol freeze,
not blinded testing or independent preregistration.

Both learning rules from the prior mortality audit are refit on pathology:
an intercept plus log body mass, log adult remaining life expectancy, and
mammal-prey diet; the extended model adds one `log1p(MHC-I family size)` term.
Ridge regularization, training-only centering/scaling, whole-order holdouts,
equal-order scoring, and complete-refit order jackknife follow the frozen
specification. No mortality coefficient is transported unchanged, and no
held-out label selects a transform or penalty. The primary unit is a species
aggregate, with species averaged within order and orders weighted equally.

Species matching replaces underscores with spaces, collapses whitespace, and
case-folds before exact equality. There are no synonyms, fuzzy matches, or
subspecies truncations. This conservative join can miss true taxonomic synonyms;
it is a fixed reproducible panel, not a claim of complete biological overlap.
Vincze supplies every order label and every predictor, including adult remaining
life after maturity rather than Compton maximum longevity. Missing mortality
outcomes never determine pathology eligibility. Zero-neoplasm species remain
in the common primary panel.

## Numerical results

The following block is generated from [the complete structured report](data/results.json).
It includes all endpoint scores, order contributions, uncertainty refits,
prespecified sensitivities, decision gates, and numerical checks. Every
reported analysis quantity is recomputed from the committed derived inputs
by the offline test.

<!-- BEGIN GENERATED RESULTS -->
### Primary result and scope

**Verdict: Inconclusive on the available pathology-prevalence panel.**

The pinned sources contain **139 Vincze species** and
**327 Compton species**; the normalized union has
**419 identifiers**, with **47 matches**.
The predictor-completeness rule excludes **5 matched species**.

The common panel contains **42 mammalian species in 9 orders**,
with **312 denominator-linked neoplasms among 2215 records**.
It retains **2 zero-neoplasm species**. The same panel has
**209 recorded malignant cases** and **4 species**
with zero recorded malignancies. It does not cover the other vertebrate classes.

### Endpoint arithmetic audit

Across **327 source species**, printed neoplasia prevalence
has **0 discrepancies** against the denominator-linked
neoplasia counts; printed malignancy prevalence has **0 discrepancies**
against `Malignant / RecordsWithDenominators`, under token-specific rounding allowances.
The dictionary's `MalignancyKnown` numerator disagrees for **182 source species**
and **34 common-panel species**.
This resolves which recorded count reproduces the published prevalence; it does not
reconstruct individual-record provenance. `PropMalignant` is unused and is never a progression probability.

### Paired held-out-order scores

Losses are natural-log cross-entropy per denominator-linked record.
The primary target averages species within each order, then orders equally.
Delta is baseline minus extended loss; positive values favor adding MHC-I.

| Endpoint / scoring target | Baseline | Baseline + MHC-I | Delta |
|---|---:|---:|---:|
| Primary neoplasia: equal orders | 0.432671423 | 0.434589384 | -0.001917961 |
| Primary neoplasia: record weighted | 0.404380234 | 0.406757816 | -0.002377582 |
| Secondary recorded malignancy: equal orders | 0.334369264 | 0.329874570 | +0.004494694 |
| Secondary recorded malignancy: record weighted | 0.320088149 | 0.322045791 | -0.001957642 |

The primary extension improves **2 of 9 orders**.
The complete-refit order jackknife reports the following uncertainty.
Its bias-corrected center is distinct from the full-panel point estimate.

| Endpoint | Full-panel delta | Jackknife center | Standard error | Approximate 95% interval |
|---|---:|---:|---:|---|
| Primary neoplasia | -0.001917961 | +0.011057504 | 0.007236988 | [-0.005631021, +0.027746029] |
| Secondary recorded malignancy | +0.004494694 | +0.011381839 | 0.008439095 | [-0.008078748, +0.030842427] |

The t critical value is **2.306004135204** with
**8 degrees of freedom** for each endpoint.
These are approximate intervals over a small selected set of related orders.
The secondary endpoint is exploratory, has no separate multiplicity-adjusted claim,
and cannot rescue the primary verdict. Unknown malignancy status is not known benign status.

### Primary neoplasia: every held-out order

| Order | Species | Records | Events | Zero species | Baseline | + MHC-I | Delta |
|---|---:|---:|---:|---:|---:|---:|---:|
| Artiodactyla | 8 | 287 | 28 | 0 | 0.381053156 | 0.381633154 | -0.000579998 |
| Carnivora | 14 | 674 | 129 | 0 | 0.526703389 | 0.526089913 | +0.000613476 |
| Chiroptera | 4 | 387 | 27 | 0 | 0.259689230 | 0.260011042 | -0.000321812 |
| Diprotodontia | 1 | 84 | 22 | 0 | 0.664166650 | 0.664436874 | -0.000270224 |
| Hyracoidea | 1 | 148 | 12 | 0 | 0.283588411 | 0.283447509 | +0.000140901 |
| Pilosa | 1 | 20 | 1 | 0 | 0.261986330 | 0.262260492 | -0.000274162 |
| Primates | 8 | 336 | 63 | 2 | 0.445308381 | 0.447221341 | -0.001912960 |
| Proboscidea | 1 | 22 | 9 | 0 | 0.688721865 | 0.697295171 | -0.008573306 |
| Rodentia | 4 | 257 | 21 | 0 | 0.382825395 | 0.388908956 | -0.006083561 |

### Secondary recorded malignancy: every held-out order

| Order | Species | Records | Events | Zero species | Baseline | + MHC-I | Delta |
|---|---:|---:|---:|---:|---:|---:|---:|
| Artiodactyla | 8 | 287 | 14 | 2 | 0.271258485 | 0.273920298 | -0.002661813 |
| Carnivora | 14 | 674 | 88 | 0 | 0.434520453 | 0.429282453 | +0.005238000 |
| Chiroptera | 4 | 387 | 19 | 0 | 0.213137030 | 0.212583793 | +0.000553236 |
| Diprotodontia | 1 | 84 | 21 | 0 | 0.710510209 | 0.713384724 | -0.002874515 |
| Hyracoidea | 1 | 148 | 8 | 0 | 0.216423345 | 0.211895488 | +0.004527857 |
| Pilosa | 1 | 20 | 1 | 0 | 0.214626995 | 0.214663942 | -0.000036947 |
| Primates | 8 | 336 | 40 | 2 | 0.304247742 | 0.302232988 | +0.002014754 |
| Proboscidea | 1 | 22 | 1 | 0 | 0.293376140 | 0.248202686 | +0.045173454 |
| Rodentia | 4 | 257 | 17 | 0 | 0.351222974 | 0.362704753 | -0.011481779 |

### Prespecified primary sensitivities and falsification attempts

| Sensitivity | Species | Orders | Baseline | + MHC-I | Delta |
|---|---:|---:|---:|---:|---:|
| Equal-species training weights | 42 | 9 | 0.434379708 | 0.436969208 | -0.002589500 |
| Ridge penalty 0.1 | 42 | 9 | 0.432803417 | 0.434701153 | -0.001897737 |
| Ridge penalty 10 | 42 | 9 | 0.431659963 | 0.433756312 | -0.002096349 |
| Remove eligible named rodents | 41 | 9 | 0.442977973 | 0.452322071 | -0.009344099 |
| Remove recorded zeros (outcome-selected) | 40 | 9 | 0.445585579 | 0.445454085 | +0.000131495 |

Named rodents actually removed: `Heterocephalus_glaber`.
Named rodents already absent from the eligible panel: `Dolichotis_patagonum`.
Only the first four sensitivities are positive-verdict gates. Zero removal is
deliberately outcome-selected and cannot rescue or veto the primary claim.

### Primary neoplasia: complete delete-order refits

The delete-order delta range is **[-0.010264705, -0.001203180]**.
Every row reruns all remaining holdouts, their model fits, and their training transformations.

| Deleted order | Remaining species | Remaining orders | Refit delta | Pseudovalue |
|---|---:|---:|---:|---:|
| Artiodactyla | 34 | 8 | -0.002611060 | +0.003626836 |
| Carnivora | 28 | 8 | -0.002814951 | +0.005257959 |
| Chiroptera | 38 | 8 | -0.003114914 | +0.007657664 |
| Diprotodontia | 41 | 8 | -0.002103086 | -0.000436955 |
| Hyracoidea | 41 | 8 | -0.002836360 | +0.005429233 |
| Pilosa | 41 | 8 | -0.002024504 | -0.001065616 |
| Primates | 34 | 8 | -0.004886283 | +0.021828622 |
| Proboscidea | 41 | 8 | -0.001203180 | -0.007636206 |
| Rodentia | 38 | 8 | -0.010264705 | +0.064855995 |

### Secondary recorded malignancy: complete delete-order refits

The delete-order delta range is **[-0.000500394, +0.008042737]**.
Every row reruns all remaining holdouts, their model fits, and their training transformations.

| Deleted order | Remaining species | Remaining orders | Refit delta | Pseudovalue |
|---|---:|---:|---:|---:|
| Artiodactyla | 34 | 8 | +0.007533936 | -0.019819241 |
| Carnivora | 28 | 8 | +0.003194366 | +0.014897315 |
| Chiroptera | 38 | 8 | -0.000153129 | +0.041677278 |
| Diprotodontia | 41 | 8 | +0.005675968 | -0.004955494 |
| Hyracoidea | 41 | 8 | +0.002626935 | +0.019436765 |
| Pilosa | 41 | 8 | +0.005163920 | -0.000859110 |
| Primates | 34 | 8 | +0.001119870 | +0.031493289 |
| Proboscidea | 41 | 8 | -0.000500394 | +0.044455400 |
| Rodentia | 38 | 8 | +0.008042737 | -0.023889647 |

### Frozen decision gates

| Gate | Outcome |
|---|---|
| all four major sensitivities positive | false |
| all four major sensitivities valid | true |
| positive jackknife lower bound | false |
| positive primary delta | false |
| strict majority orders improve | false |
| valid panel and primary uncertainty fits | true |
| point improvement failure | true |
| interval contains zero | true |
| nonimprovement gate | false |

A nonpositive point gain fails the point-improvement observable. An interval
containing zero leaves the directional claim inconclusive; an upper bound strictly
below zero contradicts improvement under this frozen procedure.

### Numerical reproducibility

Across **25 complete evaluations and 414 fits**:
all fits converged **true**;
**0 probabilities** were clipped for scoring.
Each fit must satisfy the declared maximum-gradient bound
`max(1e-8, 1e-10 * effective_training_denominator)` within 100 Newton steps.
Saved stopping diagnostics are checked against these bounds, not their exact
platform-specific values. Scientific snapshots use absolute/relative tolerance `3e-7`.
All primary and secondary predictions, training-fold parameters, coefficients,
fit certificates, refit summaries, and newline-normalized file hashes are in `data/results.json`.

**Nonclinical, unreviewed computational research. No medical or veterinary advice.**
<!-- END GENERATED RESULTS -->

## How to interpret the failed and unresolved tests

The primary point estimate, the record-weighted contrast, every primary
complete delete-order refit, and all major sensitivity gains are negative.
That pattern fails the proposed improvement test on this panel. It is not a
test of equivalence and does not establish an absence of biological MHC effects.

The primary jackknife center is positive even though the full-panel gain is
negative. These are different estimators: the frozen interval uses
`G * delta - (G - 1) * mean(delete_order_delta)` as its bias-corrected center.
Deleting an order changes both the training populations and the scoring target;
the substantial correction reflects instability of this small-order procedure.
We report the center separately and retain the full-panel gain as the point
estimate. We do not move the interval to the negative point estimate, switch
uncertainty methods after seeing the answer, or promote its positive center
into an improvement claim. Approximate jackknife coverage is not established
for these selected and related orders.

Removing the named rodent that is actually present makes the gain more
negative. The other prespecified rodent has no exact match in the source join.
Removing recorded zeros alone produces a small positive gain. That operation
conditions on the outcome and changes the target; the protocol explicitly
prevents it from rescuing the result. This is a concrete falsification attempt
against a positive claim obtained through outcome-selected exclusions.

The recorded-malignancy secondary point score is positive under equal-order
weighting, but its interval contains zero and its record-weighted contrast is
negative. Its deletion refits also include sign reversals. This remains an
exploratory secondary result, without a separate multiplicity-adjusted claim,
and cannot change the primary verdict.

## Sources, endpoint definitions, and reuse

| Input | Source and exact version | Use in this analysis |
|---|---|---|
| Genomic and life-history predictors | Vincze et al., [*Immunological surveillance against cancer across mammals*](https://doi.org/10.1038/s41467-025-65286-x); [pinned source CSV](https://github.com/OrsolyaVincze/VinczeEtalNatComm_MHC-cancer/blob/2c8a5c5b5fda6fb34c86b0e5a2b6a58dcc85cc3e/SourceData_FINAL.csv) | MHC-I family size, body mass, adult remaining lifespan, diet, and order. Mortality fields are unused. |
| Pathology aggregates | Compton et al., [*Cancer Prevalence across Vertebrates*](https://doi.org/10.1158/2159-8290.CD-24-0573); [pinned species CSV](https://github.com/zacharycompton/cancerAcrossVertebrates/blob/4819e0933358c5fc7eab483a0c05161ea8c074b7/species-cancer-prevalence-data.csv) | Primary `NeoplasiaWithDenominators / RecordsWithDenominators`; secondary recorded `Malignant / RecordsWithDenominators`. |
| Endpoint documentation | [Pinned Compton dictionary](https://github.com/zacharycompton/cancerAcrossVertebrates/blob/4819e0933358c5fc7eab483a0c05161ea8c074b7/data-description-species-cancer-prevalence-data.txt) and [source bootstrap code](https://github.com/zacharycompton/cancerAcrossVertebrates/blob/4819e0933358c5fc7eab483a0c05161ea8c074b7/bootstrap.R) | Audit conflicting numerator definitions before permitting the secondary analysis. |

### Malignancy: an operational resolution with remaining provenance limits

The dictionary proposes `MalignancyKnown` as the malignancy-prevalence numerator.
The pinned numerical prevalence instead agrees with `Malignant`, as the
per-species discrepancy audit shows. The source bootstrap code also constructs
malignant outcomes from `Malignant` indicators. Thus the secondary analysis
uses the recorded malignant count that reproduces the published prevalence.
We have resolved which column defines that observable, not the full history
of pooled records or the classification of every tumour. `MalignancyKnown`
must not be silently reinterpreted as known malignant cases or a consistently
denominator-qualified classification count. Unknown malignancy status is not
known benign status.

`PropMalignant` is unused. It is a cross-sectional source aggregate with
unresolved denominator/provenance issues and **is never read as a progression
probability**. Neither a ratio of pathological diagnoses nor this predictive
analysis supplies observations of benign-to-malignant transitions.

### Licence and access decisions

Both repository licences and both final-article licences were checked before
reading either source CSV. Complete Git-tree
inspection at both pinned revisions found no licence file; Compton's pinned
README has no reuse grant. GitHub repository metadata reported no licence for
either source. Both final articles specify
[CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/), as recorded
in the [source audit](data/source_audit.json) and prefit protocol. The separate
preprint's licence was not substituted for the final article or repository.

The original CSVs, dictionary, source R code, article text, and individual
records are not mirrored. The committed files contain only the narrow derived
species facts and audit metadata necessary for the reported results, plus
original code/documentation, following the user-authorized factual-aggregate
fallback and [DATA_POLICY.md](../../DATA_POLICY.md). This does not claim an
upstream unrestricted redistribution grant or permission under ND to share
adapted licensed material. The repository's MIT licence does not relicense
third-party assets.

All source data were read as plain UTF-8 text through GitHub's API. No browser
file download or prompt occurred. Direct PMC reads sometimes returned a
bot-check page; official indexed article text and publisher text supplied the
article licence evidence without bypassing that check. Restricted individual
animal records were neither requested nor accessed.

## What this cannot establish

This tests **endpoint portability of a measured genomic proxy after PR #79**.
It does **not** test [PR #59's functional-load theorem](https://github.com/dovvnloading/open-cancer/pull/59),
which concerns functional presented antigen loads and specified immune-escape
dynamics. A gene-family count does not measure antigen expression, peptide
binding, presentation quality, functional load, immune killing, or complete
presentation loss. The measurement bridge required to test that theorem is
absent.

**Pathology prevalence is neither onset risk nor cause-specific death.**
The source observations depend on which animals and tissues receive pathology,
their ages and survival, necropsy practices, detection and classification,
available denominators, institutional sources, and genome annotation. The
analysis does not harmonize those observation processes. Zero recorded cases
are not proof of cancer resistance.

**Ascertainment, shared zoos, and possible individual-animal overlap with
PR #79 remain limitations.** A distinct endpoint and different published
dataset do not establish disjoint animals, institutions, or observations.
Species overlap cannot resolve animal overlap. No individual identifiers were
available or sought. Order holdout prevents that order from entering its own
training fold, but does not remove all deep phylogenetic dependence or shared
ascertainment. No phylogenetic covariance model or alternative-tree uncertainty
was fit here.

No result establishes immune causation, application to wild populations, or
medical/veterinary utility. An inconclusive or negative predictive contrast
does not falsify immune surveillance. Successful computation, AI critique,
and passing CI establish neither source truth nor biological validity.

## Prior work and honest novelty assessment

The cross-species MHC-I/cancer association is prior work by Vincze and colleagues.
Compton and colleagues supplied the pathology endpoint and published comparative
analyses. PR #79 supplied the held-out-order prediction and uncertainty framework.
This contribution is the pinned endpoint join, explicit dictionary audit,
prefit comparison, and fully replayable negative/inconclusive evidence.

A focused search checked the two final articles and the relation between MHC-I
family size and comparative cancer endpoints. It did not establish a prior
identical frozen join-and-prediction experiment, but that limited search cannot
establish bibliographic priority. No first-association, new-mechanism, theorem
validation, or independent-replication claim is made. The apparent attractiveness
of a new endpoint does not turn an unresolved prediction into a discovery.

## Hostile self-review and retained limitations

| Challenge | Resolution or retained limitation |
|---|---|
| Did the panel inherit mortality eligibility or discard recorded zeros? | No mortality outcomes enter eligibility. The frozen primary panel retains zeros and uses identical rows for both models. |
| Was name reconciliation chosen to improve the score? | Exact normalization and all exclusions were published before fitting; no synonym search or post-fit expansion was used. Strict matching remains incomplete as taxonomy. |
| Did the secondary numerator follow the conflicting dictionary uncritically? | Every printed prevalence was audited at its own decimal precision. Only the arithmetic-compatible recorded malignant count is used; pooled-source classification remains unresolved. |
| Was an apparent arithmetic discrepancy manufactured by display precision? | A preliminary fixed-decimal comparison was rejected before the freeze. Token-specific precision and exact endpoint handling remove that false discrepancy. |
| Could held-out outcomes leak into transforms, penalties, or refits? | Tests inspect every training/test partition and independently recalculate transforms, gradients, probabilities, scores and jackknife quantities, including perturbed held-out labels and covariates. |
| Could negative gains be rewritten as a positive transport result? | The point-failure flag, all order scores, every complete deletion, negative major sensitivities, and the positive outcome-selected zero-removal result are all retained. The frozen verdict stays inconclusive. |
| Is the jackknife center being mistaken for the primary gain? | Both estimators are shown separately. Their differing signs and the approximate small-order uncertainty are discussed explicitly. |
| Does source arithmetic establish complete provenance? | No. The immutable input panel has valid secondary counts; its regeneration is fail-closed on altered or invalid pinned inputs. No new incomplete-secondary dataset or individual classification history was evaluated. |
| Are platform-specific optimizer traces treated as scientific constants? | Fresh and archived diagnostics must meet declared bounds. Scientific snapshots use numerical tolerances; text hashes normalize line endings and file keys use POSIX form. |
| Do licences or tests warrant a biological claim? | No raw assets are mirrored, no blanket reuse permission is asserted, and all outputs remain nonclinical and unreviewed. |

The protocol, panel, model specifications, scoring, uncertainty and verdict rule
were not changed after fitting. Documentation and executable checks were
reviewed against the frozen decisions; no favorable sensitivity replaced the
primary endpoint. The prefit commit remains in branch history.

## Reproduce and inspect

From the repository root, all analysis/test commands below run offline and use
only the Python standard library:

```bash
python -B research/mhc-pathology-prevalence/analysis.py --check
python -B -m unittest discover -s tests -p 'test_research_mhc_pathology_prevalence.py' -v
```

The analysis recomputes the complete report. Its check compares scientific
values numerically, validates optimization diagnostics against bounds, and
checks the README block against the archived report's rendering. The dedicated
test independently checks the fit equations and scoring, rather than relying
only on matching a saved snapshot. It also checks normalized hashes, source
audit summaries, zero retention, malformed input handling, endpoint arithmetic,
constant columns, finite ridge fits, held-out isolation, and decision boundaries.

Optional online source regeneration reads the pinned GitHub raw text in memory,
verifies newline-normalized hashes, and compares the derived files without
overwriting them:

```bash
python -B research/mhc-pathology-prevalence/data/fetch_verify.py --fetch --check
```

If source text is already available outside the repository as `vincze.csv`,
`compton.csv`, and `dictionary.txt`, use:

```bash
python -B research/mhc-pathology-prevalence/data/fetch_verify.py --source-dir PATH_TO_EXTERNAL_SOURCE_TEXT --check
```

The source script never opens a browser or writes upstream assets into the
repository. Network access is opt-in, separate from the offline test. Hashes
normalize CRLF/CR to LF without trimming or changing final-newline presence.
Path keys are constructed with `Path.as_posix()`. Tiny differences in regenerated
logarithms across C math libraries are checked numerically while archived file
content and its normalized hashes remain exact.

Validation environments and actual PR CI are recorded in the pull request.
The existing repository workflow discovers the new test without workflow
changes. A Windows CI run must not be inferred from line-ending/path fixtures
or a Linux run. No existing file, other research directory, or research index
was modified.

**Nonclinical, unreviewed computational research. No medical or veterinary advice.**
