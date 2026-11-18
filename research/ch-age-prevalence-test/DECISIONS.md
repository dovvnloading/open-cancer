# Extraction and implementation decisions after the freeze

**Nonclinical, unreviewed research. No medical advice.**

Recorded after inspecting source aggregates but **before the first target model
score**. These decisions do not modify the frozen protocol, fitness, mutation
supply, thresholds, or scoring weights. They are disclosed to avoid presenting
unanticipated source handling as preregistered.

1. **Source conflicts.** Wu's TET2 and ASXL1 carrier counts disagree with its
   printed percentages. Retain both in CSV. The primary diagnostic point uses
   the directly printed count/denominator; the reported observation interval is
   the hull of its Wilson interval and the alternative percentage's rounding
   interval. This conservative hull is an amendment, not a sampling confidence
   interval. Show both points and prohibit strict adjudication of these rows.
   Their residuals are conditional on the count interpretation; they are not an
   empirical gene-by-age ranking.
2. **Completed-year decade labels.** Interpret a label such as 40-49 as the
   half-open completed-year band [40,50), using its closure for conservative
   bounds. Preserve explicit inequalities and whole-cohort printed minima and
   maxima. Zhang's intended recruitment range is not established analytic
   support: leave its missing extreme-age bounds blank. Do the same for CNSR
   extreme bands and Lifelines' unreported maximum. Do not insert an assumed
   maximum or a cohort mean to obtain a point score.
   Whole-cohort integer ranges are evaluated as reported-age support using the
   printed endpoints. Sharp bounds are conditional on that reported-age model;
   unknown rounding of exact biological age is additional uncertainty.
3. **Nested and unsupported endpoints.** Only Genovese's candidate-driver
   definition enters its development diagnostic; candidate-or-unknown and
   mutation-count proxies are separately retained as context. Gene marginals
   with no frozen DFE remain context. Whole-cohort pooled context is not an
   extra age band when age rows already describe that cohort. Group scores
   and multiplier intersections by endpoint, model and evidence role.
   Liu's gene-associated n values do not explicitly establish unique-carrier
   deduplication. They are retained as ambiguous context, with no prevalence
   score, rather than assuming variant counts equal carrier counts.
4. **Denominators without cases.** A printed percentage and subgroup N do not
   authorize reverse engineering an unprinted integer case count. Retain N,
   use rounding uncertainty only, and omit binomial deviance. Zhang denominators
   are exact sums of its printed positive and negative age-table cells; neither
   within-CHIP column percentage is age prevalence. Retain the abstract-total
   discrepancy without substituting it for the table cells.
5. **Assay metadata.** Minimum allowed read depth is not mean depth; a target
   coverage is not an observed depth. CNSR's average greater-than bound and
   minimum coverage do not define a single observed D. Liu's settings referenced
   to Kar are not silently imported as independently verified settings. Only a
   printed exact mean/median plus read-count rule receives the single-depth
   illustration; it is not a calibrated caller or correction to prevalence.
6. **Independence.** German post-Watson COVID recruitment and the distinct
   Chinese stroke registry supply the two documented independent parents here.
   UKB releases share one parent and retain unresolved overlap against every
   development collection. Lifelines is a distinct Dutch parent but retains
   unresolved development-overlap status. No publication date alone proves
   participant independence. Parent-level linkage is applied before every
   holdout, including releases with different panels.
7. **Prior-art correction.** Schenz already compared its COVID age bands with
   Watson expectations. It was found during target extraction. This audit makes
   no first external-test or new-biological-theory claim. Its proposed value is
   explicit freezing, reproducible aggregate scoring, sharp age/supply bounds,
   and disclosure of the observation/identifiability gap.
8. **Numerical portability.** CSV/code/document hashes use CRLF/CR-to-LF
   normalization and POSIX path keys. Floating-point report leaves are checked
   with declared numerical tolerances; keys, strings, input hashes and displayed
   rounded report values must agree exactly. Numerical implementation agreement
   does not quantify DFE uncertainty or establish biological calibration.

No source was excluded because of model agreement. No graphical coordinates,
individual participant points, downloadable supplements, or controlled datasets
were used. Later bug fixes and hostile-review changes are recorded in the README.
