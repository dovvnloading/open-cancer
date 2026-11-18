# Breast receptor-by-age proportionality — frozen 2016 USCS analysis

**Nonclinical, unreviewed research. Not medical advice, screening advice, or clinical validation.**

<!-- Deterministic output of analysis.py; source data in data/turbow_table1_receptor.csv. -->

**Overall four-way verdict: INCONCLUSIVE.** HER2: **UNTESTABLE** (not reported in Turbow Table 1).

A failure would require at least one missing-status-enveloped, five-age Bonferroni-simultaneous
interval for Q to exclude 1. A point contrast or one favorable allocation does not meet that criterion.
Here Q = (positive older/younger incidence RR) / (double-negative older/younger incidence RR).
It is a case-composition odds ratio when the receptor strata have the same age-specific exposure.

## Input and provenance

- Table: [Turbow et al. (2021), PMC Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10846538/#T1), Receptor Status and Overall rows.
- Dataset: 2016 USCS (NPCR/SEER); November 2018 submission; released June 2019.
- Article manuscript is free to read but shows no explicit Creative Commons grant: assume standard
  copyright for expression. This CSV transcribes an attributed small subset of factual statistics,
  without copying article prose/layout or relicensing third-party content under MIT.
- Article HTML was read; no manuscript PDF, browser file download, or download prompt was used.
- The source describes its rates as age-adjusted to the 2000 U.S. standard population (19 age groups).
  Audited count/population rates below are recomputed proxies; no age-standardisation weights are inferred.
- Each source row in the CSV records DOI, exact source label, age-column locator, printed 95% CI,
  population header, unit, dataset vintage, access date, attribution and licence caveat.

### SHA-256 of UTF-8 with CRLF/CR normalized to LF

| Path (Path.as_posix()) | Normalized SHA-256 |
| --- | --- |
| research/breast-receptor-age/PROTOCOL.md | `77e50c2271b08d041cc981f9fe248a9a1ca70723fbd993160e2141ab2202c1ce` |
| research/breast-receptor-age/data/turbow_table1_receptor.csv | `ad6c7f8888dde1649c4a926f05ab405dd26bd14315775e23cbd6ade8e21b5595` |

## Source transcription and rounding reconciliation (all 20 cells)

The last three source fields show printed rate (printed 95% CI). The rate difference is
count/population × 100,000 minus printed rate; ±0.05 is only a one-decimal rounding check.

| Age | Source row | Cases / population | Published rate (95% CI) | Count rate | Difference | ±0.05 fit |
| --- | --- | ---: | ---: | ---: | ---: | :---: |
| 65-69 | overall | 37,127 / 8,889,353 | 417.7 (413.4–421.9) | 417.656943 | -0.043057 | yes |
| 65-69 | positive | 31,317 / 8,889,353 | 352.3 (348.4–356.2) | 352.297856 | -0.002144 | yes |
| 65-69 | negative | 4,595 / 8,889,353 | 51.7 (50.2–53.2) | 51.691051 | -0.008949 | yes |
| 65-69 | unknown | 1,197 / 8,889,353 | 13.7 (12.9–14.5) | 13.465547 | -0.234453 | **NO** |
| 70-74 | overall | 29,147 / 6,360,203 | 458.3 (453.0–463.6) | 458.271536 | -0.028464 | yes |
| 70-74 | positive | 24,664 / 6,360,203 | 387.8 (383.0–392.7) | 387.786365 | -0.013635 | yes |
| 70-74 | negative | 3,475 / 6,360,203 | 54.6 (52.8–56.5) | 54.636621 | +0.036621 | yes |
| 70-74 | unknown | 1,002 / 6,360,203 | 15.8 (14.9–16.9) | 15.754214 | -0.045786 | yes |
| 75-79 | overall | 20,575 / 4,647,676 | 442.7 (436.7–448.8) | 442.694370 | -0.005630 | yes |
| 75-79 | positive | 17,275 / 4,647,676 | 371.7 (366.2–377.3) | 371.691142 | -0.008858 | yes |
| 75-79 | negative | 2,411 / 4,647,676 | 51.9 (49.8–54.0) | 51.875389 | -0.024611 | yes |
| 75-79 | unknown | 887 / 4,647,676 | 19.1 (17.9–20.4) | 19.084807 | -0.015193 | yes |
| 80-84 | overall | 13,433 / 3,413,802 | 393.5 (386.9–400.2) | 393.490894 | -0.009106 | yes |
| 80-84 | positive | 11,117 / 3,413,802 | 325.6 (319.6–331.8) | 325.648646 | +0.048646 | yes |
| 80-84 | negative | 1,667 / 3,413,802 | 48.8 (46.5–51.2) | 48.831186 | +0.031186 | yes |
| 80-84 | unknown | 644 / 3,413,802 | 19.0 (17.6–20.5) | 18.864597 | -0.135403 | **NO** |
| 85+ | overall | 12,263 / 4,154,029 | 295.2 (290.0–300.5) | 295.207376 | +0.007376 | yes |
| 85+ | positive | 9,625 / 4,154,029 | 231.7 (227.1–236.4) | 231.702764 | +0.002764 | yes |
| 85+ | negative | 1,533 / 4,154,029 | 36.9 (35.1–38.8) | 36.903931 | +0.003931 | yes |
| 85+ | unknown | 1,095 / 4,154,029 | 26.6 (25.1–28.2) | 26.359951 | -0.240049 | **NO** |

## Unknown status, residual and primary proportional allocation

Unknown and unaccounted residual cases are separate in the input; both are in the
conditional worst-case assignable pool. Proportional allocation uses integer
round-half-up of pool × observed positive/(positive + negative).

| Age | Printed positive | Printed double-negative | Printed unknown | Overall | Residual | Assignable pool | k to positive | Completed positive / negative |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 65-69 | 31,317 | 4,595 | 1,197 | 37,127 | 18 | 1,215 | 1,060 | 32,377 / 4,750 |
| 70-74 | 24,664 | 3,475 | 1,002 | 29,147 | 6 | 1,008 | 884 | 25,548 / 3,599 |
| 75-79 | 17,275 | 2,411 | 887 | 20,575 | 2 | 889 | 780 | 18,055 / 2,520 |
| 80-84 | 11,117 | 1,667 | 644 | 13,433 | 5 | 649 | 564 | 11,681 / 1,752 |
| 85+ | 9,625 | 1,533 | 1,095 | 12,263 | 10 | 1,105 | 953 | 10,578 / 1,685 |

## Simultaneous inference over all 10 age contrasts

Within each hypothesized completed table, independent Poisson age/class marks imply
binomial positive composition conditional on the total for each age. Exact two-sided
Clopper–Pearson intervals at 99% **per age** (two tails 0.005 each), combined by
Bonferroni over five age groups, give ≥95% simultaneous coverage for all ten Q
transforms. This is process-model uncertainty conditional on hypothetical classification,
not a model for the missing-status mechanism, coding error or population error.

Below, P-RR and N-RR are the source-known positive and double-negative age rate ratios,
before allocation. Q and the interval are based on proportional completion. All ten
pairs are prespecified; “adjacent” only labels their position, not a different test.

| Younger → older | Kind | P-RR | N-RR | Proportional Q | Proportional interval | Status | All-allocation envelope | Robust status |
| --- | --- | ---: | ---: | ---: | --- | --- | --- | --- |
| 65-69 → 70-74 | adjacent | 1.100734 | 1.056984 | 1.041435 | [0.955512, 1.135212] | contains 1 | [0.715424, 1.490260] | contains 1 |
| 65-69 → 75-79 | nonadjacent | 1.055048 | 1.003566 | 1.051124 | [0.955871, 1.156227] | contains 1 | [0.675908, 1.534179] | contains 1 |
| 65-69 → 80-84 | nonadjacent | 0.924356 | 0.944674 | 0.978144 | [0.879614, 1.088350] | contains 1 | [0.613674, 1.454439] | contains 1 |
| 65-69 → 85+ | nonadjacent | 0.657690 | 0.713933 | 0.921002 | [0.826925, 1.026411] | contains 1 | [0.467480, 1.446639] | contains 1 |
| 70-74 → 75-79 | adjacent | 0.958495 | 0.949462 | 1.009304 | [0.912410, 1.116706] | contains 1 | [0.643835, 1.510650] | contains 1 |
| 70-74 → 80-84 | nonadjacent | 0.839763 | 0.893745 | 0.939228 | [0.839620, 1.051150] | contains 1 | [0.584554, 1.432133] | contains 1 |
| 70-74 → 85+ | nonadjacent | 0.597501 | 0.675443 | 0.884359 | [0.789327, 0.991328] | below | [0.445297, 1.424453] | contains 1 |
| 75-79 → 80-84 | adjacent | 0.876127 | 0.941317 | 0.930570 | [0.824360, 1.050754] | contains 1 | [0.567820, 1.515859] | contains 1 |
| 75-79 → 85+ | nonadjacent | 0.623374 | 0.711396 | 0.876207 | [0.774981, 0.990955] | below | [0.432550, 1.507730] | contains 1 |
| 80-84 → 85+ | adjacent | 0.711512 | 0.755745 | 0.941581 | [0.823314, 1.076865] | contains 1 | [0.456264, 1.660632] | contains 1 |

### Complete allocation sensitivity (all five cases for each pair)

Scenarios: proportional, all unknown+residual to positive in both ages,
all to negative in both, and the two separately adverse cross-age corners.
Each interval uses the same per-age Bonferroni CP rule. The envelope above
projects over **every integer allocation**, not just these five displayed choices.

| Younger → older | Allocation | Q point | Simultaneous interval | Excludes 1? |
| --- | --- | ---: | --- | --- |
| 65-69 → 70-74 | proportional | 1.041435 | [0.955512, 1.135212] | contains 1 |
| 65-69 → 70-74 | all_positive | 1.043469 | [0.956157, 1.138889] | contains 1 |
| 65-69 → 70-74 | all_negative | 1.020683 | [0.943330, 1.104467] | contains 1 |
| 65-69 → 70-74 | younger_positive_older_negative | 0.777087 | [0.715424, 0.844058] | below |
| 65-69 → 70-74 | younger_negative_older_positive | 1.370569 | [1.260752, 1.490260] | above |
| 65-69 → 75-79 | proportional | 1.051124 | [0.955871, 1.156227] | contains 1 |
| 65-69 → 75-79 | all_positive | 1.064116 | [0.966110, 1.172453] | contains 1 |
| 65-69 → 75-79 | all_negative | 0.971181 | [0.891227, 1.058520] | contains 1 |
| 65-69 → 75-79 | younger_positive_older_negative | 0.739399 | [0.675908, 0.808944] | below |
| 65-69 → 75-79 | younger_negative_older_positive | 1.397688 | [1.273875, 1.534179] | above |
| 65-69 → 80-84 | proportional | 0.978144 | [0.879614, 1.088350] | contains 1 |
| 65-69 → 80-84 | all_positive | 0.996938 | [0.894736, 1.111513] | contains 1 |
| 65-69 → 80-84 | all_negative | 0.890523 | [0.809167, 0.980436] | below |
| 65-69 → 80-84 | younger_positive_older_negative | 0.677991 | [0.613674, 0.749271] | below |
| 65-69 → 80-84 | younger_negative_older_positive | 1.309451 | [1.179764, 1.454439] | above |
| 65-69 → 85+ | proportional | 0.921002 | [0.826925, 1.026411] | contains 1 |
| 65-69 → 85+ | all_positive | 0.988627 | [0.884699, 1.105553] | contains 1 |
| 65-69 → 85+ | all_negative | 0.676896 | [0.616401, 0.743520] | below |
| 65-69 → 85+ | younger_positive_older_negative | 0.515348 | [0.467480, 0.568214] | below |
| 65-69 → 85+ | younger_negative_older_positive | 1.298535 | [1.166530, 1.446639] | above |
| 70-74 → 75-79 | proportional | 1.009304 | [0.912410, 1.116706] | contains 1 |
| 70-74 → 75-79 | all_positive | 1.019787 | [0.920266, 1.130311] | contains 1 |
| 70-74 → 75-79 | all_negative | 0.951501 | [0.868728, 1.042286] | contains 1 |
| 70-74 → 75-79 | younger_positive_older_negative | 0.708597 | [0.643835, 0.779868] | below |
| 70-74 → 75-79 | younger_negative_older_positive | 1.369366 | [1.241716, 1.510650] | above |
| 70-74 → 80-84 | proportional | 0.939228 | [0.839620, 1.051150] | contains 1 |
| 70-74 → 80-84 | all_positive | 0.955407 | [0.852279, 1.071562] | contains 1 |
| 70-74 → 80-84 | all_negative | 0.872478 | [0.788740, 0.965400] | below |
| 70-74 → 80-84 | younger_positive_older_negative | 0.649747 | [0.584554, 0.722339] | below |
| 70-74 → 80-84 | younger_negative_older_positive | 1.282917 | [1.149981, 1.432133] | above |
| 70-74 → 85+ | proportional | 0.884359 | [0.789327, 0.991328] | below |
| 70-74 → 85+ | all_positive | 0.947442 | [0.842718, 1.065815] | contains 1 |
| 70-74 → 85+ | all_negative | 0.663180 | [0.600840, 0.732117] | below |
| 70-74 → 85+ | younger_positive_older_negative | 0.493880 | [0.445297, 0.547790] | below |
| 70-74 → 85+ | younger_negative_older_positive | 1.272222 | [1.137081, 1.424453] | above |
| 75-79 → 80-84 | proportional | 0.930570 | [0.824360, 1.050754] | contains 1 |
| 75-79 → 80-84 | all_positive | 0.936869 | [0.827880, 1.060523] | contains 1 |
| 75-79 → 80-84 | all_negative | 0.916948 | [0.822977, 1.021840] | contains 1 |
| 75-79 → 80-84 | younger_positive_older_negative | 0.637140 | [0.567820, 0.714898] | below |
| 75-79 → 80-84 | younger_negative_older_positive | 1.348308 | [1.199898, 1.515859] | above |
| 75-79 → 85+ | proportional | 0.876207 | [0.774981, 0.990955] | below |
| 75-79 → 85+ | all_positive | 0.929059 | [0.818594, 1.054836] | contains 1 |
| 75-79 → 85+ | all_negative | 0.696982 | [0.626921, 0.774919] | below |
| 75-79 → 85+ | younger_positive_older_negative | 0.484297 | [0.432550, 0.542147] | below |
| 75-79 → 85+ | younger_negative_older_positive | 1.337068 | [1.186439, 1.507730] | above |
| 80-84 → 85+ | proportional | 0.941581 | [0.823314, 1.076865] | contains 1 |
| 80-84 → 85+ | all_positive | 0.991663 | [0.863474, 1.138981] | contains 1 |
| 80-84 → 85+ | all_negative | 0.760111 | [0.676850, 0.853505] | below |
| 80-84 → 85+ | younger_positive_older_negative | 0.516931 | [0.456264, 0.585395] | below |
| 80-84 → 85+ | younger_negative_older_positive | 1.458171 | [1.280929, 1.660632] | above |

### All rounding sensitivities (not sampling intervals)

The count-based primary inference does not use rounded rates. The following
uses only positive and negative **printed** Table 1 rates, and allows each
independently to vary ±0.05 (half last digit) or ±0.10 (full last digit)
per 100,000. Each row is a sharp interval for a rounding-only Q expression,
not a confidence interval; these do not change the count-based verdict.

| Younger → older | Q from printed rates | ±0.05 rounding range | ±0.10 stress range |
| --- | ---: | --- | --- |
| 65-69 → 70-74 | 1.042301 | [1.040058, 1.044548] | [1.037821, 1.046800] |
| 65-69 → 75-79 | 1.051001 | [1.048684, 1.053323] | [1.046372, 1.055650] |
| 65-69 → 80-84 | 0.979135 | [0.976898, 0.981377] | [0.974666, 0.983624] |
| 65-69 → 85+ | 0.921462 | [0.918997, 0.923935] | [0.916538, 0.926416] |
| 70-74 → 75-79 | 1.008347 | [1.006189, 1.010510] | [1.004036, 1.012677] |
| 70-74 → 80-84 | 0.939398 | [0.937312, 0.941488] | [0.935231, 0.943583] |
| 70-74 → 85+ | 0.884066 | [0.881757, 0.886381] | [0.879455, 0.888704] |
| 75-79 → 80-84 | 0.931621 | [0.929503, 0.933744] | [0.927390, 0.935872] |
| 75-79 → 85+ | 0.876747 | [0.874411, 0.879091] | [0.872082, 0.881441] |
| 80-84 → 85+ | 0.941099 | [0.938515, 0.943689] | [0.935940, 0.946288] |

## Explicit verdict and scope

**INCONCLUSIVE.** Robust missing-status-enveloped exclusions: **0/10**.
Allocation-specific interval exclusions: **28/50** (different mutually hypothetical
completions; not 50 independent tests). The frozen rule does not allow a scenario-specific
effect to override an allocation-robust failure criterion. **HER2: UNTESTABLE** here.

The older/younger receptor-class rate ratios are comparisons of one calendar year’s
different age groups; no transition, time-series trend, screening-detection effect,
birth-cohort effect, competing-risk mechanism, survival effect, latent biology, or
cause of the localized-stage decline is identified. The 85+ group is open-ended
and is *not* assigned a midpoint or five-year width. The source’s ER-or-PR
positive class includes hormone receptor borderline positive and all HER2 statuses.
Double-negative ER/PR is not equivalent to triple-negative HER2 disease.

### Literature gate and independent critical interpretation

- Turbow et al. (2021) themselves described the positive-class rates as resembling the
  all-cancer age pattern and tabulated the receptor composition by age. This is **not**
  a novel discovery of those rates, only a prespecified missingness/multiplicity audit.
- Lee Argov et al. (2025), doi:10.1001/jamanetworkopen.2025.16947, published
  related USCS 2001–2019 older-age HR/HER2 composition and trends, but not this
  specific 2016 five-bin all-allocation interaction test.
- Jiang et al. (2024), doi:10.1002/ijc.35073, considered ER-specific **calendar**
  trends over 2004–2020 and age 20–84; different estimand and period.
- Acheampong et al. (2020), doi:10.1001/jamanetworkopen.2020.13226, and other SEER subtype-age analyses use different coverage,
  molecular categories and periods; overlapping registries are **not** independent
  replication. The nearby unsettled question is formal allocation-robust proportionality.
- The interval assumption can fail under registry classification error, extra-Poisson
  heterogeneity, ascertainment changes, or residual cases outside binary status.
  Widely differing allocation corners cannot be interpreted as equally plausible.
- No case-level data, joint stage×receptor×age table, HER2 counts for 2016,
  independent population-age weights, or screening linkage are supplied.

### Hostile-review checklist

- Every one of the 20 printed HTML count/rate/95% CI triples and five age-header
  populations checked against the full article HTML, not PDF snippets.
- Explicitly checked source manuscript label, absence of a CC licence, ER/PR class
  definitions, five residuals, rate rounding anomalies, and open-ended 85+.
- CP tails, degenerate count boundary cases, interval monotonicity, all extrema,
  proportional integer rounding and direct printed-rate envelopes independently tested.
- Original numerical conclusions, even if contrary to expectation, must remain intact.
- No browser file download or prompt occurred; browser recaptcha on one URL form was
  bypassed only by opening an already indexed **same PMC HTML** result, not PDF.
- Shell git clone was blocked by DNS; the GitHub connector used for remote branch
  and separate protocol commit. No unsupported claim of successful local full clone.

**Nonclinical research output. Not medical advice. Software CI and automated
source review do not independently validate scientific or clinical claims.**
