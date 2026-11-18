# Can published CH fitness predict age prevalence?

**Nonclinical, unreviewed research. No medical advice.** This directory concerns
population aggregate observations and a mathematical model. It cannot estimate
an individual's disease risk, diagnose CH, or guide treatment. Passed code checks
and AI critique are not external scientific or clinical validation.

## Question, scope, and current interpretation

**Verdict: UNTESTABLE for fitness alone with these published aggregates.** The
specified compound prediction transported to the available substituted
endpoints **fails** the frozen interval-compatibility criterion. This is a
quantitative negative transport result, with a reproducible test, rather than
evidence that fitness is irrelevant. Missing panel and assay equivalence and
missing matched gene-by-age rows prevent the stronger biological verdict.

An additional model-only check fails: the unchanged separate gene and pooled
fits sometimes violate the necessary inequality that a gene's carrier
prevalence cannot exceed prevalence for any gene in its containing panel.
The tabulated curves are therefore separate candidate forecasts, not one
coherent joint model. This was not repaired by refitting or changing supply.

The intended test is unchanged-parameter prediction of the age prevalence of
clones carrying a specified nonsynonymous SNV, above a specified true VAF
threshold, using Watson's previously published fitness distributions and
mutation supplies. The strict result and the complete descriptive scores are
generated below from the CSVs. Curated-driver CHIP measured by different assays
is an explicitly substituted endpoint. Its disagreement tests that compound
transport prediction; it does not isolate an error in fitness.

This is the closest feasible open variant of the question under the permitted
aggregate-only access: a frozen transport audit with a sharp statement of what
age bands, incomplete assays, and mutation-supply confounding can identify.
There is no claim of a new branching model or the first external comparison.
The gene aggregates recovered here are broad cohort marginals, not a quantitative
validation of gene-specific age curves. The missing matched curves remain an
explicit empirical gap.

## What was already known

| Primary work | What it already established or supplied | Consequence for this test |
|---|---|---|
| [Watson et al., Science 2020](https://doi.org/10.1126/science.aay9333), Figure 3 and Figure S18 | Predicted prevalence as a function of age and detection threshold; compared the DNMT3A R882H/C age trend; inferred within-gene fitness distributions from VAF spectra. | The age-prevalence idea and these comparisons are prior art. Jaiswal and Genovese were development inputs. A whole-gene average is not a substitute for the published distribution. |
| [Jaiswal et al., NEJM 2014](https://doi.org/10.1056/NEJMoa1408617) and [Genovese et al., NEJM 2014](https://doi.org/10.1056/NEJMoa1409405) | Published strong age associations and some exact age-band prevalence values. | Useful development diagnostics, never independent holdouts. Their panels, variant classes and detection rules differ from the frozen target. |
| [Bick et al., Nature 2020](https://doi.org/10.1038/s41586-020-2819-2) | TOPMed age/gene prevalence and empirical low-VAF detection limitations. | Relevant external evidence, but exact age rows were not recovered in readable text. Plots were not digitized. Some parent-cohort overlap remains unresolved. |
| [Schenz et al., 2022](https://doi.org/10.3389/fimmu.2022.968778), Results and Figure 2 | Already compared observed hospitalized-COVID age groups with Watson expectations. | This defeats a claim to the first independent external comparison. Recompute the frozen model, preserve disease selection, and do not borrow their median-age test as ours. |
| [Fabre et al., Nature 2022](https://doi.org/10.1038/s41586-022-04785-z) | Longitudinal clonal growth depends on gene and age, with early and late growth patterns. | PR #68's logistic/class coefficients are not lifelong intrinsic per-gene fitness estimates. Fixed fitness is a transport assumption that may fail. |
| [Mitchell et al., Nature 2022](https://doi.org/10.1038/s41586-022-04786-y) | Haematopoietic phylogenies, age-related oligoclonality and a selection model. | A colony/phylogeny sample is not a population prevalence table at a common VAF threshold. No participant data are extracted. |
| [CH-Predict, Blood Advances 2026](https://doi.org/10.1182/bloodadvances.2025017940) | External prediction using demographic and blood-count features across UKB, All of Us and MSK. | External CH prediction already exists. It is not a test of unchanged Watson fitness and supply alone. |
| [Mon Père et al., Cancer Discovery 2026](https://doi.org/10.1158/2159-8290.CD-25-0990) and [Wang et al., 2026 preprint](https://doi.org/10.64898/2026.01.27.702030) | Newer competition and branching-inference work addresses selection, mutation and environment. | Intrinsic fitness and observed net expansion need not coincide. These works prevent an unsupported novelty claim. |

The source audit also includes newer timing/structure work and the recent
[Watson pre-AML paper](https://doi.org/10.1038/s41586-026-11269-x). Searches ran
through **2026-11-18**. This is a targeted primary-source survey, not an exhaustive
systematic review or proof of priority. [Source audit](data/source_audit.csv)
records accessible material and exclusions. The original
[CH progression note](../clonal-hematopoiesis-progression/README.md) is read-only
background; it and the research index are unchanged.

## Frozen protocol and audit trail

[PROTOCOL.md](PROTOCOL.md), [model parameters](data/model_parameters.csv) and
[mutation supplies](data/mutation_supply.csv) were locally committed as
`6d00600fc57919ab89e9408717a65526fe069a66` on **2026-11-18 at 02:18:27 EST**,
before deliberate target age-table extraction and before any target scoring.
The source base was main `5993e6cc474958bd7645721ce1771b6bc7848f08`.
The published branch begins with the normalized equivalent protocol commit
[`bdbf5d30cd25f3c8e6739148f922bde1d22b2a80`](https://github.com/dovvnloading/open-cancer/commit/bdbf5d30cd25f3c8e6739148f922bde1d22b2a80).
The local CSV writer used CRLF; publication used LF, so tree/commit IDs differ
while the newline-normalized contents agree. [Freeze evidence](freeze_evidence.json)
retains the original Git commit payload and content hashes. A self-authored
commit timestamp is an audit trail, not third-party preregistration evidence.

The aggregate extraction and conservative source-handling decisions were then
locally committed as `16ac1c5` and published as
[`066b4489bdf190a61691df1c1715df38404ea6c7`](https://github.com/dovvnloading/open-cancer/commit/066b4489bdf190a61691df1c1715df38404ea6c7)
before the first target score. Later implementation fixes preserve these
inputs. The finite-band-only supply-intersection check was added after the
first score as a labelled robustness check; it changes no primary result.

The freeze was **not fully blinded**: broad search snippets had already shown
some overall/oldest-age summaries, and development comparisons had been read as
prior art. The protocol identifies that exposure. No target fit, fitness refit,
effective threshold selection, onset-age shift or mutation-supply calibration
is hidden in the primary result. The explicitly secondary supply calibration
uses training parents only.

[DECISIONS.md](DECISIONS.md) records source-handling decisions made after
extraction but before scoring, including the conservative source-conflict hull,
completed-year age bands, and the newly found Schenz prior comparison.
These are disclosed amendments/clarifications; the frozen files are immutable.

## Data and measurement bridge

[Prevalence CSV](data/prevalence.csv) contains only small printed population
aggregates. [Cohort CSV](data/cohorts.csv) records panel, VAF rule, depth statistic,
read support, ascertainment, development membership and shared parents.
[Context aggregates](data/context_aggregates.csv) retain nested definitions and
genes without a frozen model separately. Every row, including design fixtures,
has a source, locator, access date and reuse basis.

- **Jaiswal:** use explicit main-text prevalence/intervals, with cases only
  when aligned in that text. The figure-only numeric table was not visually
  verified, so its additional rows and counts are not admitted. Different SNV
  and indel detection limits prevent a single source-threshold score; the
  reference comparison is labelled accordingly.
- **Genovese:** candidate-driver age rows are the development diagnostic.
  Candidate-or-unknown and mutation-count proxy definitions share participants
  and are not pooled into that endpoint's score or feasibility intersection.
- **Bick/TOPMed:** age/gene plots and low-VAF sensitivity information are
  meaningful, but an exact readable age table was not recovered. The nominal
  threshold is not its effective sensitivity. No plotted height or unprinted
  age denominator was inferred.
- **Zhang UKB:** derive each denominator by summing printed CHIP-positive and
  CHIP-negative age cells. The abstract total disagrees with Methods/Table 1.
  Intended recruitment ages are insufficient to assert this analytic subset's
  exact support; the extreme bands remain open. The variant filter's use of
  age association is an additional ascertainment caveat.
- **Wu and Liu UKB:** Wu supplies broad gene marginals over its printed age
  support; TET2 and ASXL1 count/percentage contradictions remain visible. Liu
  supplies gene-associated counts including TP53, but unique-carrier
  deduplication is not explicit. Its ambiguous counts remain context with no
  prevalence score and do not repair Wu. All UKB releases share one parent.
- **Qiu/CNSR-III and Schenz/COVID:** printed age aggregates are available in
  distinct independent cohorts. Both are disease-selected. Their panels and
  assays differ from the Watson target; the conditional descriptive scores
  do not become population calibration or a test of fitness alone.
- **Lifelines:** the oldest-age cohort supplies an exact pooled count and some
  gene aggregates, but no verified maximum age. Preserve an open upper band,
  no uniform-age point score, and unresolved development-overlap status.

Independence here is an audit of the published named cohort provenance against
Watson's development collections, not a participant-identity deduplication
proof. No participant-level data were sought to establish it.

Source URLs and exact locators live with every row. Direct PMC/publisher access
sometimes returned reCAPTCHA or access errors; readable indexed primary text
or public author-hosted article text was used where available. No challenge was
bypassed. No browser file download, supplemental spreadsheet, raw genome,
participant record, individual figure point, or plotted coordinate was used.
Unavailable supplemental/readout details remain missing. No download prompt
was accepted.

The original CSV schema contains isolated numeric facts, not mirrored article
text, figure artwork or source-table formatting. Copyright and source licenses
are retained in the provenance fields. NC/ND licenses are not asserted to permit
adaptation of protected expression, public access is not a blanket reuse
license, and this repository's MIT license does not relicense third-party
material. The work is a noncommercial factual research audit.

## Model and numbered results

The frozen model integrates each published truncated stretched-exponential
fitness distribution against Watson's threshold kernel and Poisson mutation
supply. It uses the finite-fraction coordinate `x=v/(1-2v)` and the published
population-times-turnover parameter, as specified in the protocol. Numerical
normalization is evaluated accurately; the source fitness cutoff is retained
without renormalizing it away. The pooled model is a separately published fit,
not a sum of the separately recovered gene fits.

[THEOREMS.md](THEOREMS.md) supplies numbered statements and proofs:

1. **DFE normalization and Watson E1 reduction.** Reduce the VAF integral to
   a one-dimensional integral and prove monotonicity of the resulting age curve.
   Exactness is relative to the stated approximate Watson density.
2. **Exact integer birth-death threshold law.** Obtain the cross-sectional
   probability that a clone beginning with one cell has at least a specified
   number of cells at a later time, including critical/pure-birth boundaries.
3. **Integrated threshold kernel.** Integrate that probability over arrival
   times, with a positive-series remainder bound and a late-time asymptotic.
4. **Poisson arrivals and observation.** Derive the probability of at least
   one detected clone under independent arrivals/clones and a specified caller.
5. **Sharp age-bin and common-supply bounds.** Unknown age distributions in the stated reported-age support give
   sharp endpoint bounds; each observed interval restricts a common mutation
   multiplier, and their intersection tests whether one multiplier can work.
6. **Identifiability limits.** One point can be matched by free supply;
   constant independent assay thinning aliases supply. Net fitness alone does
   not fix exact branching tails without turnover. These statements do not
   rule out every restricted age-shape test.
7. **Necessary subset consistency.** A gene-specific carrier event is contained
   in the event of a clone in the pooled panel. Its probability cannot exceed
   the pooled probability, regardless of dependence among genes. The frozen
   separate fits violate this condition when interpreted jointly.

The exact birth-death law does not prove the finite-population VAF mapping,
constant fitness, independent clones, or transportability of mutation supply.
A survival/first-crossing probability is not substituted for cross-sectional
clone presence. Prenatal initiation, competition, survivorship, inflammation,
and gene/caller selection can violate the compound transport model.

## Computed results

<!-- BEGIN RECOMPUTED RESULTS -->

**Strict biological verdict: UNTESTABLE.**
Eligible matched rows: **0**; independent eligible parents: **0**. The descriptive substituted-endpoint result **fails frozen interval-compatibility criterion**.

External substituted endpoints: **12 rows**, **7 finite-band point scores**, **4 represented parent cohorts**, of which **2 contribute point scores**. With each source's printed threshold (or the explicitly labelled reference when missing), equal-parent MAE is **6.7480 percentage points**, RMSE **9.7242 pp**; **6 rows** have nonoverlapping observation and age-prediction intervals.

These errors use an assumed uniform age distribution within finite bands. The intervals below retain unknown age composition. Rounding intervals do not supply sampling uncertainty. Model sensitivities do not supply fitted-parameter confidence intervals.

### External independence strata

A parent is in the independent stratum only if every audited linked release is labelled independent. The other stratum retains unresolved or known overlap; the combined summary is descriptive, not independent validation.

| Independence stratum | Rows | Represented parents | Point-score parents | MAE pp | RMSE pp | Nonoverlap rows |
|---|---:|---:|---:|---:|---:|---:|
| Independent as audited | 9 | 2 | 2 | 6.7480 | 9.7242 | 5 |
| Overlap unresolved or linked | 3 | 2 | 0 | unavailable | unavailable | 1 |

### Complete age-band comparison

Predictions use the pooled ten-gene fit for any-CH endpoints. These endpoints are substitutions for Watson's variant panel. Error = observed minus uniform-age prediction.

Rows marked reference use 2% because a uniform source cutoff is unavailable. Published/Wilson labels denote marginal 95% intervals; rounding-only rows have no computed sampling interval. Every interval is conditional on reported-age support.

| Row | VAF % / basis | Observed % | Observation interval % / kind | Prediction %: age range; uniform mean | Error pp | Gap pp |
|---|---|---:|---|---:|---:|---:|
| jaiswal2014_60_69 (60-69) | 2 / reference only | 5.6000 | 5.0000–6.3000 / published 95% | 9.0248–21.5668; 14.9330 | -9.3330 | 2.7248 |
| jaiswal2014_70_79 (70-79) | 2 / reference only | 9.5217 | 8.4000–10.8000 / published 95% | 21.5668–36.9361; 29.1429 | -19.6211 | 10.7668 |
| jaiswal2014_80_89 (80-89) | 2 / reference only | 11.6719 | 8.6000–15.7000 / published 95% | 36.9361–52.3226; 44.7236 | -33.0517 | 21.2361 |
| jaiswal2014_90_108_text (90-108) | 2 / reference only | 18.4466 | 12.1000–27.0000 / published 95% | 52.3226–74.7117; 64.2377 | -45.7911 | 25.3226 |
| genovese2014_candidate_drivers_lt50 (<50) | 2 / reference only | 0.7000 | 0.6500–0.7500 / rounding only | 0.0000–1.9193; 0.2255 | 0.4745 | 0.0000 |
| genovese2014_candidate_drivers_gt65 (>65) | 2 / reference only | 5.7000 | 5.6500–5.7500 / rounding only | 14.7510–56.6306; 35.4486 | -29.7486 | 9.0010 |
| zhang2025_ukb_lt60 (<60) | 2 / source | 1.8037 | 1.7522–1.8566 / Wilson 95% | 0.0000–9.0248; unavailable | unavailable | 0.0000 |
| zhang2025_ukb_ge60 (>=60) | 2 / source | 4.6765 | 4.5858–4.7689 / Wilson 95% | 9.0248–100.0000; unavailable | unavailable | 4.2559 |
| schenz2022_covid_lt45 (<45) | 1 / source | 0.0000 | 0.0000–32.4408 / Wilson 95% | 0.0000–2.3618; 0.4177 | -0.4177 | 0.0000 |
| schenz2022_covid_45_65 (45-65) | 1 / source | 32.0000 | 31.5000–32.5000 / rounding only | 2.3618–24.6795; 11.6208 | 20.3792 | 6.8205 |
| schenz2022_covid_gt65 (>65) | 1 / source | 56.3000 | 56.2500–56.3500 / rounding only | 24.6795–70.3790; 48.7088 | 7.5912 | 0.0000 |
| qiu2023_cnsr_lt40 (<40) | 5 / source | 0.5848 | 0.1033–3.2377 / Wilson 95% | 0.0000–0.0003; unavailable | unavailable | 0.1030 |
| qiu2023_cnsr_40_49 (40-49) | 5 / source | 1.8571 | 1.0885–3.1513 / Wilson 95% | 0.0003–0.2055; 0.0450 | 1.8121 | 0.8829 |
| qiu2023_cnsr_50_59 (50-59) | 5 / source | 2.1699 | 1.5643–3.0027 / Wilson 95% | 0.2055–2.7594; 1.1369 | 1.0330 | 0.0000 |
| qiu2023_cnsr_60_69 (60-69) | 5 / source | 3.5201 | 2.8001–4.4168 / Wilson 95% | 2.7594–10.3335; 6.1009 | -2.5808 | 0.0000 |
| qiu2023_cnsr_70_79 (70-79) | 5 / source | 5.3963 | 4.2484–6.8322 / Wilson 95% | 10.3335–22.4648; 16.1040 | -10.7077 | 3.5013 |
| qiu2023_cnsr_ge80 (>=80) | 5 / source | 10.9422 | 8.0092–14.7768 / Wilson 95% | 22.4648–100.0000; unavailable | unavailable | 7.6879 |
| lifelines2021_all_source_genes (>=80 whole cohort) | 1 / source | 61.5137 | 57.6272–65.2586 / Wilson 95% | 49.3574–100.0000; unavailable | unavailable | 0.0000 |

### Cohort scores

Alternative endpoints, model panels and whole-cohort gene marginals remain separate.

| Cohort / model / endpoint / role | Point rows | MAE pp | RMSE pp | Nonoverlap rows |
|---|---:|---:|---:|---:|
| genovese2014 / POOLED10 / candidate_drivers / development | 2 | 15.1116 | 21.0381 | 1 |
| jaiswal2014 / POOLED10 / curated_160_gene_SNV_indel / development | 4 | 26.9492 | 30.2545 | 4 |
| lifelines2021 / ASXL1 / targeted_27_gene_CH / external_gene_marginal | 0 | unavailable | unavailable | 1 |
| lifelines2021 / DNMT3A / targeted_27_gene_CH / external_gene_marginal | 0 | unavailable | unavailable | 1 |
| lifelines2021 / POOLED10 / targeted_27_gene_CH / external_age | 0 | unavailable | unavailable | 0 |
| lifelines2021 / TET2 / targeted_27_gene_CH / external_gene_marginal | 0 | unavailable | unavailable | 0 |
| lifelines2021 / TP53 / targeted_27_gene_CH / external_gene_marginal | 0 | unavailable | unavailable | 1 |
| qiu2023_cnsr / POOLED10 / curated_74_gene_CHIP / external_age | 4 | 4.0334 | 5.6051 | 4 |
| schenz2022_covid / DNMT3A / targeted_driver_CHIP / external_gene_marginal | 1 | 7.6797 | 7.6797 | 0 |
| schenz2022_covid / POOLED10 / targeted_driver_CHIP / external_age | 3 | 9.4627 | 12.5580 | 1 |
| schenz2022_covid / TET2 / targeted_driver_CHIP / external_gene_marginal | 1 | 7.0421 | 7.0421 | 0 |
| wu2025_ukb / ASXL1 / curated_58_gene_CHIP / external_gene_marginal | 1 | 2.3706 | 2.3706 | 0 |
| wu2025_ukb / DNMT3A / curated_58_gene_CHIP / external_gene_marginal | 1 | 7.7177 | 7.7177 | 0 |
| wu2025_ukb / TET2 / curated_58_gene_CHIP / external_gene_marginal | 1 | 5.2451 | 5.2451 | 0 |
| zhang2025_ukb / POOLED10 / curated_55_gene_CHIP / external_age | 0 | unavailable | unavailable | 1 |

### Leave-one-parent-out checks

With fixed external inputs, deleting a cohort changes the aggregate error, never the remaining predictions. This is an influence check.

| Parent removed | Remaining point rows | MAE pp | RMSE pp |
|---|---:|---:|---:|
| CNSR_III | 3 | 9.4627 | 12.5580 |
| German_COVID_hospital_2020_2021 | 4 | 4.0334 | 5.6051 |
| Lifelines | 7 | 6.7480 | 9.7242 |
| UK_Biobank | 7 | 6.7480 | 9.7242 |

The secondary calculation changes only the common mutation-supply multiplier using other parents. All rows sharing the held-out parent are excluded from calibration. This is a descriptive transport check across unmatched endpoints, not primary validation. The combined calculation retains a non-independent interpretation whenever any parent has unresolved or linked overlap.

Combined-calculation status: **descriptive secondary calibration**.

| Held-out parent | Independence interpretation | Training-only supply multiplier | Held-out MAE pp | Held-out RMSE pp | Constant-baseline MAE pp | Constant-baseline RMSE pp |
|---|---|---:|---:|---:|---:|---:|
| CNSR_III | non-independent interpretation: unresolved or linked cohort overlap | 2.041687 | 8.6563 | 12.9269 | 40.9142 | 40.9379 |
| German_COVID_hospital_2020_2021 | non-independent interpretation: unresolved or linked cohort overlap | 2.562237 | 11.1162 | 15.6524 | 28.3547 | 34.8982 |

The same calibration is also computed separately within each independence stratum; a stratum with fewer than two finite-band parents has no LOCO fit.

| Independence stratum | Finite-band parents | Status |
|---|---:|---|
| independent | 2 | descriptive secondary calibration |
| overlap_uncertain | 0 | untestable: fewer than two parents with finite-band rows |

| Independence stratum / held-out parent | Supply multiplier | MAE pp | RMSE pp | Constant MAE pp | Constant RMSE pp |
|---|---:|---:|---:|---:|---:|
| independent / CNSR_III | 2.041687 | 8.6563 | 12.9269 | 40.9142 | 40.9379 |
| independent / German_COVID_hospital_2020_2021 | 2.562237 | 11.1162 | 15.6524 | 28.3547 | 34.8982 |

### Common-supply compatibility

The intersection uses every admitted age interval within the same endpoint, model and role. Finite-band bounds are sharp under their stated assumptions. Open or missing age support instead gives a conservative closure; a zero lower multiplier can be an unattained infimum. An empty intersection rules out a common multiplier under the stated kernel/measurement assumptions. These interval boxes are not a simultaneous confidence test.

| Cohort / model / endpoint / role | Lower multiplier | Upper multiplier | Nonempty? | Contains frozen multiplier 1? |
|---|---:|---:|---|---|
| genovese2014 / POOLED10 / candidate_drivers / development | 0.336507 | 0.371062 | True | False |
| jaiswal2014 / POOLED10 / curated_160_gene_SNV_indel / development | 0.211150 | 0.370456 | True | False |
| lifelines2021 / ASXL1 / targeted_27_gene_CH / external_gene_marginal | 0.000000 | 0.553891 | True | False |
| lifelines2021 / DNMT3A / targeted_27_gene_CH / external_gene_marginal | 0.000000 | 0.702930 | True | False |
| lifelines2021 / POOLED10 / targeted_27_gene_CH / external_age | 0.000000 | 1.553904 | True | True |
| lifelines2021 / TET2 / targeted_27_gene_CH / external_gene_marginal | 0.000000 | 1.171964 | True | True |
| lifelines2021 / TP53 / targeted_27_gene_CH / external_gene_marginal | 0.000000 | 0.761782 | True | False |
| qiu2023_cnsr / POOLED10 / curated_74_gene_CHIP / external_age | 297.545536 | 0.628432 | False | False |
| schenz2022_covid / DNMT3A / targeted_driver_CHIP / external_gene_marginal | 0.097804 | 31444.874122 | True | True |
| schenz2022_covid / POOLED10 / targeted_driver_CHIP / external_age | 1.334904 | 2.924887 | True | False |
| schenz2022_covid / TET2 / targeted_driver_CHIP / external_gene_marginal | 0.061693 | 9829.415910 | True | True |
| wu2025_ukb / ASXL1 / curated_58_gene_CHIP / external_gene_marginal | 0.034093 | 2.787067 | True | True |
| wu2025_ukb / DNMT3A / curated_58_gene_CHIP / external_gene_marginal | 0.069365 | 10.649562 | True | True |
| wu2025_ukb / TET2 / curated_58_gene_CHIP / external_gene_marginal | 0.027772 | 2.522105 | True | True |
| zhang2025_ukb / POOLED10 / curated_55_gene_CHIP / external_age | 0.186901 | 0.516619 | True | False |

A post-score robustness check repeats the same intersection using only finite external age bands. It changes no model input or primary score and cannot be used to select which bands to retain.

| Finite-band cohort | Lower multiplier | Upper multiplier | Nonempty? |
|---|---:|---:|---|
| qiu2023_cnsr | 5.319347 | 0.648814 | False |
| schenz2022_covid | 1.334904 | 2.924887 | True |

### Broad gene marginals

These are whole-cohort gene counts or printed percentages over their full age support, not gene-by-decade curves. Their uniform-age residuals cannot identify the most discordant gene-specific age trend. Conflicting source percentages are retained in the observation hull.

| Gene / source / ages | Observed % | Printed alternative % | Prediction age range %; uniform mean % | Uniform-assumption error pp | Interval gap pp | Source conflict? |
|---|---:|---:|---:|---:|---:|---|
| DNMT3A / wu2025_ukb / 37-73 (whole cohort marginal) | 2.3966 | 2.4000 | 0.2318–29.0518; 10.1143 | -7.7177 | 0.0000 | False |
| TET2 / wu2025_ukb / 37-73 (whole cohort marginal) | 0.5619 | 0.4700 | 0.2319–15.4497; 5.8070 | -5.2451 | 0.0000 | True |
| ASXL1 / wu2025_ukb / 37-73 (whole cohort marginal) | 0.3512 | 0.2500 | 0.1325–6.9423; 2.7218 | -2.3706 | 0.0000 | True |
| DNMT3A / schenz2022_covid / 24-95 whole cohort | 17.7778 | 17.8000 | 0.0010–70.4838; 25.4574 | -7.6797 | 0.0000 | False |
| TET2 / schenz2022_covid / 24-95 whole cohort | 6.6667 | 6.7000 | 0.0015–39.8864; 13.7087 | -7.0421 | 0.0000 | False |
| DNMT3A / lifelines2021 / >=80 whole cohort | 35.2657 | unavailable | 50.6195–100.0000; unavailable | unavailable | 11.5158 | False |
| TET2 / lifelines2021 / >=80 whole cohort | 26.7311 | unavailable | 26.5511–100.0000; unavailable | unavailable | 0.0000 | False |
| ASXL1 / lifelines2021 / >=80 whole cohort | 6.0000 | unavailable | 11.4267–100.0000; unavailable | unavailable | 4.9267 | False |
| TP53 / lifelines2021 / >=80 whole cohort | 3.0000 | unavailable | 4.5691–100.0000; unavailable | unavailable | 1.0691 | False |

### Largest conditional deviations

Finite-band age errors ranked by absolute error: schenz2022_covid_45_65: 20.3792 pp; qiu2023_cnsr_70_79: -10.7077 pp; schenz2022_covid_gt65: 7.5912 pp.

Largest age-interval gaps, including open bands: qiu2023_cnsr_ge80: 7.6879 pp; schenz2022_covid_45_65: 6.8205 pp; zhang2025_ukb_ge60: 4.2559 pp.

Largest gene-marginal interval gaps: DNMT3A (lifelines2021, >=80 whole cohort): 11.5158 pp; ASXL1 (lifelines2021, >=80 whole cohort): 4.9267 pp; TP53 (lifelines2021, >=80 whole cohort): 1.0691 pp. These identify disagreements of substituted endpoints under their stated bounds. They do not establish which gene's fitness or age trend is wrong.

### All frozen parameter sensitivities

Each scenario applies its common threshold to all external age rows. It is a hypothetical sensitivity, not a choice of the best cohort-specific cutoff.

| Scenario | MAE pp | RMSE pp | Nonoverlap rows |
|---|---:|---:|---:|
| reference | 12.5077 | 16.2949 | 6 |
| vaf_01 | 13.0475 | 17.4608 | 6 |
| vaf_03 | 12.3783 | 17.0423 | 8 |
| vaf_05 | 12.5730 | 18.7235 | 6 |
| vaf_10 | 13.2388 | 21.3884 | 6 |
| A_70000 | 12.6794 | 16.5159 | 6 |
| A_130000 | 12.2854 | 16.2713 | 6 |

### Read-support illustrations

These are binomial read-sampling probabilities, not caller sensitivity estimates. Other filters, variable depths and variant classes remain unmodelled.

| Cohort | Illustrative depth / alternate reads | Detection at true VAF 2% | 5% | 10% |
|---|---:|---:|---:|---:|
| lifelines2021 | 9726 / 10 | 100.0000% | 100.0000% | 100.0000% |

### Falsifiable matched-assay predictions

Fixed-age prevalence (%) for the source variant panels under perfect detection at true VAF >= 2%. These are predictions, not observations or clinical probabilities. A future matched aggregate should supply age support, carrier counts, denominators and caller calibration. For an age band, use its two endpoint predictions as the age-only bounds. These are separately fitted candidate forecasts, not a coherent joint carrier model; the subset-order check immediately below exposes the inconsistency.

| Model | Age 40 | Age 50 | Age 60 | Age 70 | Age 80 | Age 90 |
|---|---:|---:|---:|---:|---:|---:|
| DNMT3A | 0.5898 | 4.1495 | 12.4010 | 24.8061 | 39.3469 | 53.7837 |
| TET2 | 0.5255 | 2.8157 | 7.2000 | 13.3401 | 20.7281 | 28.8421 |
| ASXL1 | 0.2893 | 1.4136 | 3.3976 | 6.0461 | 9.1790 | 12.6527 |
| TP53 | 0.1076 | 0.5401 | 1.3170 | 2.3712 | 3.6417 | 5.0812 |
| POOLED10 | 0.0749 | 1.9193 | 9.0248 | 21.5668 | 36.9361 | 52.3226 |

### Internal consistency of the separate published fits

If a gene's variants are a subset of the pooled endpoint, its carrier prevalence cannot exceed the prevalence of any pooled-gene clone. This necessary event-inclusion inequality does not assume independent genes. The unchanged separate DFE fits violate it in **10 of 24** grid comparisons; the largest excess is **3.3762 pp**. Therefore these fits cannot jointly describe one common population/panel as exact probabilities. No fit is repaired or selected after seeing targets. This structural incompatibility concerns the assembled candidate predictions; it is not evidence against every individual fitness estimate or all branching models.

| Gene | Age | Gene prevalence % | Pooled prevalence % | Excess pp |
|---|---:|---:|---:|---:|
| DNMT3A | 40 | 0.5898 | 0.0749 | 0.5150 |
| DNMT3A | 50 | 4.1495 | 1.9193 | 2.2302 |
| DNMT3A | 60 | 12.4010 | 9.0248 | 3.3762 |
| DNMT3A | 70 | 24.8061 | 21.5668 | 3.2392 |
| DNMT3A | 80 | 39.3469 | 36.9361 | 2.4108 |
| DNMT3A | 90 | 53.7837 | 52.3226 | 1.4611 |
| TET2 | 40 | 0.5255 | 0.0749 | 0.4506 |
| TET2 | 50 | 2.8157 | 1.9193 | 0.8964 |
| ASXL1 | 40 | 0.2893 | 0.0749 | 0.2145 |
| TP53 | 40 | 0.1076 | 0.0749 | 0.0328 |

The maximum difference from the source notebook's VAF upper cutoff over this prediction grid is **0.00158410 pp**. This checks a numerical convention; it does not measure biological model error.

Report verification allows floating-point differences of at most relative **5e-09** or absolute **5e-10** for platform math-library and numerical-integration variation. Integer values, metadata, keys and input hashes are exact; this displayed results block must reproduce exactly. These are computational tolerances, not statistical uncertainty.

<!-- END RECOMPUTED RESULTS -->

## What would make the test decisive?

Publish independent, age-stratified **carrier counts and denominators** for
DNMT3A, TET2, ASXL1 and TP53, plus the exact Watson nonsynonymous-SNV coding
panels and the pooled panel, with declared age boundaries and true-VAF cutoffs.
Provide aggregate validation of detection probability by VAF, read depth,
variant class and relevant gene coverage; mean depth alone is insufficient.
A finer age distribution within bins would replace the descriptive uniform-age
assumption. Development-membership and parent-cohort overlap must be declared.

The generated future-prediction table is a falsifiable unchanged-parameter
prediction at the stated ages and reference threshold. For a finite reporting
band, use the generated curve's endpoint range, not the centre-age value as an
exact bin average. A reproducible mismatch with an appropriately calibrated,
matched endpoint rejects that compound frozen prediction. Shared failure
across independent genes/ages would be stronger evidence against the model
than a selected cohort or a nominal-threshold mismatch.

A second threshold measured in the same age/gene bin, with known relative
sensitivity, constrains shape beyond a constant supply multiplier. Independent
caller calibration or an externally measured mutation supply can break the
constant-thinning/supply alias. Arbitrary unmeasured age-dependent sensitivity
would need stronger observation assumptions; it does not license fitting it
post hoc to every discrepancy.

## Reproduce and review

From the repository root, with standard-library Python only:

```sh
python research/ch-age-prevalence-test/analysis.py
python research/ch-age-prevalence-test/verify.py
python -m unittest discover -s tests -p test_research_ch_age_prevalence_test.py -v
```

`analysis.py --write` explicitly regenerates the numerical report, only the
marked results section of this README, and the payload manifest. It does not
alter frozen inputs or the separate freeze-evidence hashes.
The committed report exposes per-row predictions, intervals, deviances,
scenario results, held-out predictions, exclusions and CSV hashes. The verifier
uses independent mathematical calculations, including a stopped continuous-time
Markov chain with controlled tail omission and separate quadrature formulas.
The mathematical fixtures are labelled original/synthetic and are not empirical
CH observations. Algorithmic agreement is not an estimate of biological or
parameter uncertainty.

`manifest.json` hashes all owned code, data, documents and the single test file,
excluding itself. All keys are built with `Path.as_posix()` and all hashes
normalize CRLF/CR to LF. Float comparisons use documented tolerances so platform
math-library differences do not cause spurious byte-equality failures. Rendered
rounded results must reproduce exactly. No network access or external packages
are needed to run the checks.

Hostile review checked prior-art claims, source count conflicts, nested endpoint
mixing, uncertainty labels, parent-level holdout leakage, finite-age eligibility,
mathematical boundary cases, and Windows newline/path handling. It corrected
proof-domain assumptions and conservative source handling. Any failed primary
data gate stays a failed gate; descriptive agreement cannot upgrade it.
