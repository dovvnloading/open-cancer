# Outcome-source audit and extraction record

Nonclinical, unreviewed research. No medical advice. Access date for this
record: 2026-11-17, America/New_York. This is a source and estimand audit,
not independent verification of the underlying cohorts or a clinical validation.

## Included: Fabre's published AML coefficient file

**Primary source:** Fabre et al., *Nature* (2022),
[article and methods](https://doi.org/10.1038/s41586-022-04785-z), especially
“Association between clonal haematopoiesis dynamics and AML risk” and Figure 5a.
The article's Data availability statement links
[de Almeida, Fabre, Gerstung and Vassiliou, Data, version 1](https://doi.org/10.6084/m9.figshare.15029118.v1).
The [Figshare metadata API](https://api.figshare.com/v2/articles/15029118)
identifies the dataset license as
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The article separately carries CC BY 4.0. Source attribution and license remain
attached to every extracted row; the repository's MIT license does not replace
those source terms.

The file inventory allowed retrieval of **only** the small aggregate file
[`aml-risk-coef.csv`](https://ndownloader.figshare.com/files/28917537).
The other archive entries include individual sequencing, phenotyping and
phylogeny records and were not downloaded. No bulk archive was requested.
All source rows whose `group` equals `Genes` are retained in `data/outcomes.csv`;
demographic coefficients are not gene prediction targets and are omitted.
No sign, significance, precision or model agreement filter was applied.
The CSV retains the original `coef` and `sd` strings as `raw_log_estimate`
and `log_standard_error`. A second fetch of this aggregate source was compared
with every retained coefficient and standard error after initial extraction.

### Scale and uncertainty

The source methods describe Cox proportional-hazards coefficients, adjusted for
age, sex and cohort, **agnostic of clone size**. Coefficients are interpreted
on the Cox log-hazard scale; negative source coefficients also rule out treating
the raw column as a hazard ratio. The exported ratio is `exp(coef)`.
The source also contains a standardized-statistic column consistent, to source
rounding, with coefficient divided by `sd`.

The displayed lower and upper limits in the extracted CSV are our Wald
transformation `exp(coef - 1.96*sd)` and `exp(coef + 1.96*sd)`. They are labeled
**derived from published coefficient and standard error**, not transcribed
published confidence limits. The public analysis
[notebook](https://github.com/josegcpa/clonal_dynamics/blob/master/Notebook_GrowthCoefficients_AgeAtOnset_PossibleAssociations.Rmd)
reads this exact file and plots `coef +/- 1.96*sd`. Figure 5's caption instead
describes AML error bars as standard errors. That caption/code discrepancy is
preserved rather than silently resolved by calling the plotted bars a
particular published confidence level. The original coefficient and `sd`
allow a reader to use a different explicitly labeled multiplier.

### What this estimand does not provide

- The file is a refit agnostic of clone size, not the original Abelson
  per-VAF coefficient table. A per-VAF coefficient cannot replace a
  gene-presence coefficient.
- The exact analytic denominator for this particular refit is not supplied
  in the coefficient file or its Fabre methods paragraph. The corresponding
  cohort-size fields are blank with a reason. Related Abelson sample and
  participant totals are not silently substituted for a verified refit N.
- The complete Cox formula, design matrix, any coefficient covariance and
  any regularization details are not in this aggregate file. The methods
  establish adjustment for age, sex and cohort but do not unambiguously say
  whether all gene indicators were included jointly in this particular
  refit. We therefore do not assert either a separate univariate model per
  gene or a fully specified common-reference multivariable model. The CSV
  reference text means the indicated gene absent under the fitted model;
  it does not establish that all gene-negative reference populations have
  identical biological risk.
- No gene-specific VAF distribution, current stem-cell burden, progression
  latency, age-specific growth trajectory, hidden-hit distribution or
  baseline hazard is supplied. A Cox coefficient is not a fixed-horizon
  cumulative-hazard ratio under the branching model.
- Fabre's original Figure 5 uses pooled gene, site and unknown-cause growth
  effects. Supplementary Table 6's baseline mutation-class coefficients are
  not that pooled quantity. The present aggregate stress test cannot claim
  to reproduce Fabre's original correlation or their exact predictor.

These limits remain even if a fixed-slope stress test fits well. A failure
does not identify which transport assumption, covariate imbalance or
biological mechanism caused it.

## Included: Watson's published high-VAF aggregate

**Source:** Watson, Poon, MacGregor et al., *Nature* (2026),
[Evolutionary dynamics in the years preceding acute myeloid leukaemia](https://doi.org/10.1038/s41586-026-11269-x),
Figure 5f and Methods “High VAF variants in the UK Biobank.” The article's
Rights and permissions section explicitly grants CC BY 4.0. Only its main-text
published aggregates were extracted; its participant-level supplementary
tables, sequencing data, Zenodo records and controlled-access resources were
not accessed.

`data/outcomes.csv` preserves the reported sweep odds ratio and confidence
limits. `data/cohorts_outcomes.csv` preserves the published nested counts
needed to audit the reference group. A sweep means at least one qualifying
hotspot variant **strictly above** the stated VAF threshold. The endpoint is
prospective **myeloid leukaemia, ICD-10 C92**, not AML alone. Eligibility
excluded previous cancers as specified by the source and other subsequent
haematological malignancies. The source tests sweeps against **all other
eligible participants**, including participants with no detected hotspot CH.

The scatterplot depicts hotspot carriers only. Its carrier controls must not
be mistaken for the population-wide controls used for the published odds
ratio. The analysis reconstructs the population-wide cross-product odds
ratio from the published counts and compares it with the rounded reported
estimate. It also computes a separate within-carrier sweep-versus-nonsweep
odds ratio. The latter is **derived here**, with no newly asserted confidence
interval and no claim that it was a published source estimate. No individual
records or independently inferred subgroups are needed for these count
identities.

Neither comparison identifies the clone-size proportionality theorem:
within-bin VAF moments, driver/age composition, conversion rates and a common
prediction horizon are unavailable. The VAF aggregate is a direct check of
reference-group transport, not an extra gene-level observation in the fitness
regression.

## Prespecified or searched outcome sources not transcribed into the fit

| Source | What was checked | Reason for exclusion from numerical regression |
|---|---|---|
| [Abelson et al. 2018](https://doi.org/10.1038/s41586-018-0317-6) | Article, methods and aggregate Supplementary Table 3 model definitions; main paper and supplement contain different per-VAF and model coefficients | Source-specific permission for redistribution of an extracted table was not established from the research-use manuscript notice. The explicitly CC-licensed Fabre/author aggregate refit is used with its own estimand instead. Original per-VAF coefficients are not silently substituted. |
| [Desai et al. 2018](https://doi.org/10.1038/s41591-018-0081-z) | Main article gene-group odds-ratio and VAF findings; supplement descriptions | Unrestricted or otherwise applicable table-redistribution terms were not verified. Gene-group estimates also do not provide separate fitness-matched estimates for each member. Bibliographic/qualitative comparison only; no source workbook downloaded. |
| [Weeks et al. 2023](https://doi.org/10.1056/EVIDoa2200310) | Article/abstract and CHRS estimand, mutation/VAF/age and blood-index variables | Reuse permission for extracted numerical tables was not verified; availability of an author manuscript in PMC is not a redistribution license. The MN/CHRS endpoint and reference groups also differ from AML gene-presence coefficients. No numerical series is represented as independent validation. |
| [Gu et al. 2023](https://doi.org/10.1038/s41588-023-01472-1) | CC BY 4.0 article and Supplementary Table 7, identified in the supplement description as ranked final-model HRs | This aggregate table supplies point estimates but not their numerical uncertainty limits; the confidence intervals are graphical. The frozen specification prohibits inventing limits or digitizing plots. The models also distinguish AML, MDS and MPN and split some mutation classes; no endpoint pooling is performed. |

These exclusions concern the requirements of this repository and this frozen
analysis. They do not claim that reuse of isolated statistical facts is
universally prohibited, nor that the omitted studies are scientifically
invalid. No source was excluded because its reported direction was
inconvenient. Missing permissions or estimands are explicit abstentions.

### Mixed-workbook screening incident

The Gu supplemental workbook is described as containing both aggregate model
tables and external validation results. An in-memory schema/preview operation
unexpectedly encountered participant-level validation worksheets after the
aggregate tables. This was a screening mistake. Those records were not used
for any computation, and no workbook or participant rows were written to the
repository. Further processing of those worksheets stopped. No automated
source refresh or download of this mixed workbook is part of the rerun path.
The included CSVs contain only the separate rights-cleared aggregates named
above. This incident is recorded so the audit does not imply that every source
screening step was free of error.

## Remaining source and transport risks

The observation units in the two included outcome sources differ. Their
effects must not be pooled. Shared subjects, overlapping genes and jointly
fitted effects may produce dependent uncertainties, so the deterministic
envelopes are not joint statistical prediction intervals. Published aggregate
values can be transcribed accurately while their underlying estimates remain
confounded, biased or nontransportable. Passing the computational checks does
not resolve any of those biological or statistical limitations.
