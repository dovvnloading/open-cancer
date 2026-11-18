# Source-availability decisions frozen before fitting

**Nonclinical, unreviewed research. No treatment advice.**

This document was frozen after source extraction and before the first
empirical model fit. It supplements, and does not replace, [PROTOCOL.md](PROTOCOL.md).
The original protocol was published before new retrieval as
[5ce3b2261a05a9763f79a7942b3f7fe30e108033](https://github.com/dovvnloading/open-cancer/commit/5ce3b2261a05a9763f79a7942b3f7fe30e108033).
Its GitHub tree was fetched into the local branch before implementation.
The earlier local protocol-only commit had an extra final blank line;
the public protocol is authoritative and is unchanged by this work.
The original local commit was `fee90f0f68b2956e78a4d541c9bbb8a99c3ad0eb`.
All checksum calculations normalize CRLF and CR to LF, without stripping
other whitespace. No backdated publication or independent registry is claimed.

## Source choices and categorical matching

The source audit followed the protocol's availability/reuse hierarchy:

* The original Yarchoan numeric panel did not have established affirmative
  redistribution terms. Mo et al. (2023) provides a coherent licensed
  published response table with Yarchoan lineage. Chalmers supplies exact
  published median TMB values. This is a cross-source pairing, not a
  reproduction of Mo's MSK-IMPACT analysis or of the original Yarchoan fit.
* The Chalmers category rule prefers the named entity or NOS category when
  one exists. Explicit morphology/site approximations remain flagged.
  Cervical squamous carcinoma, bladder urothelial carcinoma, pleural
  mesothelioma, prostate acinar adenocarcinoma, soft-tissue sarcoma NOS and
  testicular nonseminoma are approximations to broader response labels.
  Anal, head-and-neck, ovarian and endometrial morphology mismatches are
  conservatively flagged too. No synthetic medians are constructed.
  MSI/MSS colorectal, MSI noncolorectal and pooled esophagogastric outcomes
  remain retained with missing TMB.
* Charoentong's Table S3 supplies a cancer-type fraction of predicted
  neoantigens assigned clonal origin, not a patient median and not a
  calibrated harmonic functional load. Its indexed primary table is
  inspectable, but the live PDF currently fails. This single-source
  limitation is retained, not presented as visual or byte-level verification.
* DASH supplies the coherent primary HLA panel. Its printed percentages
  take priority over aggregate figure measurements. Montesion's licensed
  republication and TRACERx remain separately identified panels.

The exact per-row TMB mappings and original labels appear in `tmb_orr.csv`.
Clonal joining uses the stated TCGA code and keeps OV/STAD missing. Combined
NSCLC receives no artificial average of lung-specific clonal fractions.
The KIRC-to-general-renal, PAAD-to-general-pancreatic and analogous broader
mapping assumptions are separately flagged in the analysis.

Primary HLA categorical joins, fixed before fitting:

| ORR row ID | DASH source row ID |
| --- | --- |
| urothelial | dash2022_bladder |
| breast | dash2022_breast |
| cervical | dash2022_cervical |
| head_neck | dash2022_hnsc |
| renal_cell | dash2022_kidney |
| hepatocellular | dash2022_liver |
| ovarian | dash2022_ovarian |
| pancreatic | dash2022_pancreatic |
| skin_melanoma | dash2022_melanoma |
| endometrial | dash2022_uterine |

Kidney, liver and uterine are organ labels, not proof of histologic identity.
All these broadenings remain explicit and are removed by the relevant
exact-mapping sensitivity. Bladder versus general urothelial and skin versus
unspecified melanoma require the same caution. NSCLC and SCLC share a
parent lung holdout group. No colorectal biomarker strata receive a generic
colorectal predictor. Unreported selection is unknown, never unselected.

## The full-model gate is preserved

The coherent primary sources yield nine complete parent cancer types. This
falls below the protocol's ten-parent full-model gate. **The preregistered
full test is unavailable; the reduced reciprocal-TMB comparison remains the
primary available test.** The clonal-only and HLA-only gates remain unchanged.

To make the available full-model information reviewable, specify these
additional analyses now, before fitting, and report all of them:

1. Fit the full model on the coherent nine-parent subset as **exploratory,
   below the preregistered gate**, with matched baselines, grouped prediction,
   coefficient uncertainty and the relaxed HLA sign diagnostic. This does
   not waive the gate or supply a confirmatory verdict.
2. Construct one **exploratory mixed-source** HLA sensitivity using the
   Montesion republication wherever it has a matching category and DASH
   otherwise. This fixed source-priority rule is applied to every match,
   not just to a cancer that changes a coefficient. Montesion matches are
   breast, head-neck, pancreatic, prostate, ovarian, combined NSCLC and
   mesothelioma; all remaining eligible HLA rows use DASH. Its complete
   full-model subset contains ten parents because prostate is added.
   Different stages, assays and study selection remain confounded with
   cancer type. This panel is not promoted to the primary analysis.
3. In that mixed complete-case panel, fit the free clonal exponent jointly
   with HLA and report its full fit and all parent-deletion estimates.
   This is a secondary, potentially weakly identified diagnostic.
   The coherent clonal-only exponent diagnostic is also retained.

Matched reciprocal-TMB ablations and matched single-term reciprocal models
are reported alongside the required linear/log-linear comparisons. They
distinguish shape changes from covariate additions. All are descriptive;
no opportunistic p-values, selected winning panel, or multiple-testing
significance claim is introduced.

## Numerical implementation clarifications before fitting

* The response scale is probability, and natural log is used for log-linear
  TMB. Zero responses remain valid observations. Linear baseline predictions
  are clipped only when scored/predicted; the fitted coefficients are OLS.
* Nonlinear fitting profiles out the bounded ceiling parameters and searches
  the entire prespecified nonnegative coefficient range, including boundaries.
  Every grid-local minimum is refined, and a denser grid checks final fits.
  Synthetic examples were used to test the optimizer before empirical fits;
  those examples are not evidence and are not included as source data.
* Parameter uncertainty uses 1,000 group-bootstrap draws, seed 20261117.
  Singular resamples and per-parameter usable counts are disclosed.
  Paired held-out error differences use 10,000 group resamples with that seed.
  Their percentile intervals are conditional descriptive intervals, not
  calibrated external-prediction or patient-sampling confidence intervals.
* Each free-exponent diagnostic holds its feature set fixed across parent
  deletions. It never adds HLA merely because deletion removes a missing row.
  The deletion range is an influence diagnostic, not a confidence interval.
* Rounding/reading sensitivity uses 100 deterministic uniform perturbations
  within the row-specific extraction envelopes, seed 20261118. The reported
  min/max is the range over those scenarios, not an exhaustive bound or CI.
  Percentages printed as zero without decimals receive the more conservative
  half-percentage-point envelope, clipped at zero.
* DASH's published marginal prevalence intervals are preserved separately
  from digitization envelopes. Two additional scenarios set all available
  HLA prevalences to their published lower or upper endpoints. These extremes
  are stress scenarios, not a joint confidence region or error propagation.
  Their construction is not inferred to be binomial/Wilson without a source.
* Selection-exclusion sensitivity cannot establish an unselected panel when
  the original pooled trial rows have unknown selection. Its unavailable
  status will be shown explicitly. Removing the indexed-only clonal source
  leaves no clonal/full panel; it cannot be reported as a negative test.

No fit, coefficient, held-out error or outlier ranking was used to choose
these sources, mappings, gates or additions. This chronology does not erase
the investigators' access to published source values, the prior notes or
known ecological findings.
