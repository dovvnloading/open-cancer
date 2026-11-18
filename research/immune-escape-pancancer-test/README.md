# Cancer-type test of the immune-escape reciprocal-burden prediction

**Nonclinical, unreviewed computational research. No treatment, diagnosis,
surveillance or patient-selection advice.** These are published cancer-type
aggregates, not individual observations. Passing tests, source fingerprints
and agreement among AI agents do not validate biology or clinical utility.

## Outcome

**The proposed form does not demonstrate predictive value beyond TMB alone
in this published aggregate panel.** The primary reciprocal-TMB surrogate
predicts held-out parent cancers worse than the log-linear TMB baseline.
The fixed clonal-fraction extension performs worse on its matched subset;
its freely estimated clonal exponent has the opposite direction to the
proposed positive-coverage bridge. The HLA penalty has the expected point
direction, but its uncertainty permits no penalty and the opposite direction.

The coherent complete-case panel does not meet the preregistered full-model
sample gate. Its below-gate fit and an explicitly exploratory mixed-source
fit are reported anyway, with matched baselines; neither rescues a positive
predictive claim. The clonal finding is **descriptive contrary evidence to
this ecological observation bridge**, not evidence that clonal neoantigens
cause worse response. The conditional stochastic theorem in PR59 is neither
empirically validated nor mathematically refuted by these data.

All empirical numbers below are generated from the included CSVs by
[analysis.py](analysis.py). The result is a real, reproducible test of a
specified cancer-type surrogate, with explicit limits on what was testable.

## Prediction derived before new retrieval

Read first:
[PR59's immune-escape note](../immune-escape-neoantigen-clonality/README.md),
[PR63's prediction confrontation](../prediction-confrontation/README.md),
and [DATA_POLICY.md](../../DATA_POLICY.md).

PR59 Theorem 6, conditional on source suppression, fixed positive escape
rates, no pre-existing presentation-null cells and increasing minimum
functional antigen load, gives

\[
Q \simeq \exp\!\left[-\frac{\sigma h}{\kappa}
                    \sum_{i=1}^{N}\frac{1}{k_i}\right].
\]

For homogeneous load this becomes \(Q\simeq\exp[-\sigma hN/(\kappa k)]\).
The sharp variable is the harmonic-load exposure, or \(N/k\) in that special
case. **Neither median TMB nor a clonal-ID fraction identifies that exposure.**
Theorem 5 explicitly permits identical per-cell loads and identical outcomes
despite different clonal fractions. This analysis therefore proposes and
tests an additional measurement/observation model; it does not rename TMB
as the theorem's measured input.

For cancer type \(c\), let \(T_c\) be median mutations/Mb, \(f_c\) the
published clonal neoantigen fraction and \(L_c\) any-HLA-I-LOH prevalence.
Make the restrictive working assumptions that functional load is proportional
to \(T_cf_c\), effective founder burden and the rate/conversion factors combine
into a common \(B\), and objective response is a common multiplicative
surrogate for eradication. Distinguish the population with and without LOH:

\[
\boxed{p_c=[A_0(1-L_c)+A_1L_c]\exp[-B/(T_cf_c)]},
\quad 0\le A_1\le A_0\le1,\quad B\ge0.
\]

Equivalently, for positive \(A_0\),
\(p_c=A_0(1-\eta L_c)e^{-B/(T_cf_c)}\), where
\(\eta=1-A_1/A_0\). LOH prevalence is a population-mixture weight,
**not the transition hazard \(h\)**. Ordinary HLA LOH leaves residual
presentation pathways; it is not PR59's complete loss of all modeled targets.
The theorem supplies no numerical prediction for an ordinary-LOH penalty.

### Sharpest predeclared consequence

With \(C(L)=A_0(1-L)+A_1L\), the adjusted log-response must be linear in
reciprocal clonal burden:

\[
\log\{p/C(L)\}=-B/(Tf),\qquad
\frac{p(2Tf,L)}{C(L)}=\sqrt{\frac{p(Tf,L)}{C(L)}}.
\]

The reciprocal exponent is fixed at unity. At matched TMB and LOH, increasing
the clonal fraction cannot lower response; at matched TMB and clonality,
increasing LOH prevalence cannot raise it. A free clonal exponent tests that
direction and the fixed exponent as a secondary diagnostic. A relaxed HLA
model allows the ceiling ordering to reverse rather than enforcing the
desired sign and then calling it evidence.

These are falsifiable consequences of the theorem **plus** strong endpoint,
measurement and common-parameter assumptions. ORR includes partial shrinkage;
it is not eventual eradication. The multiplicative bridge assumes no
background response unrelated to modeled eradication. The large-load error
is not calibrated here. Unknown founder counts, presentation quality and
immune strength cannot be estimated separately from the fitted composite.
Also \(E[e^{-B/K}]\ne e^{-B/\operatorname{median}(K)}\): substituting medians
from unrelated cohorts into a response probability is an ecological
approximation, not a mathematical averaging identity.

## Protocol chronology and what changed

The branch starts from main
`2b03737d05c5647587a70671ab10bb4bc969cff3`.
[PROTOCOL.md](PROTOCOL.md) was published in the
[protocol-only commit](https://github.com/dovvnloading/open-cancer/commit/5ce3b2261a05a9763f79a7942b3f7fe30e108033)
before any new source retrieval. It specifies the equations, signs, source
hierarchy, missing-data rules, model gates, baselines, prediction units,
uncertainty and improvement criterion. This is prospective planning for
retrospective public data, **not a blinded study or independent registry**:
the requested prior notes already exposed motivating published findings.

After source extraction, the
[pre-fit amendment commit](https://github.com/dovvnloading/open-cancer/commit/4155af90522219d8a3b2d038e6c329dd75336b09)
froze the exact source/category choices and exploratory additions before the
first empirical fit. [AMENDMENTS.md](AMENDMENTS.md) preserves those decisions.
It does not lower the full-model gate. The below-gate coherent fit and the
Montesion-priority/DASH-fallback fit remain exploratory throughout. No
coefficient, held-out error or outlier ranking selected the source panel,
category mappings or these additions. Source access dates use
America/New_York and are recorded on each row.

## Published aggregates and provenance

| CSV | Numeric source and exact location | Meaning and reuse |
| --- | --- | --- |
| [tmb_orr.csv](data/tmb_orr.csv) | [Mo et al. 2023](https://doi.org/10.3389/fcell.2023.1209243), Supplement Table S1, printed PDF page 5; [Chalmers et al. 2017](https://doi.org/10.1186/s13073-017-0424-2), Additional file 3, Sheet1, exact cells recorded per row | Anti-PD-1/PD-L1 aggregate ORRs with historical Yarchoan lineage, paired to published median TMB. Mo CC BY; Chalmers CC BY with a CC0 data waiver unless otherwise stated. |
| [clonal_fractions.csv](data/clonal_fractions.csv) | [Charoentong et al. 2017](https://doi.org/10.1016/j.celrep.2016.12.019), Supplementary Table S3, printed Table 3, Clonal/Subclonal columns | Published cancer-type fractions of predicted neoantigens assigned clonal origin; not patient medians. CC BY-NC-ND; unchanged isolated factual values with attribution and unit conversion. |
| [hla_loh.csv](data/hla_loh.csv) | [Marty Pyke et al. 2022](https://doi.org/10.1038/s41467-022-29203-w), Results and Figure 5a; [Lozac'hmeur et al. 2024](https://doi.org/10.1038/s41698-024-00665-z), Table 1 Montesion columns; [McGranahan et al. 2017](https://doi.org/10.1016/j.cell.2017.10.001), prevalence/timing Results | Any-HLA-I-LOH prevalence; coherent DASH panel, separately labeled Montesion republication and overlapping TRACERx context. DASH/McGranahan CC BY; Lozac'hmeur CC BY-NC-ND factual excerpts. HLA-A-only estimates are excluded. |

Each CSV row carries its source, locator, access date, definitions, units,
population and reuse basis. TMB and ORR have separate provenance on the same
row. Original values and extraction precision are preserved; missing is
blank, never zero. The source papers, figures and original workbooks are not
redistributed. Repository MIT does not relicense third-party observations.
The compilation is noncommercial research under the attached source terms.

The detailed audits are [TMB/ORR](TMB_ORR_SOURCES.md),
[clonality](CLONAL_SOURCES.md) and [HLA](HLA_SOURCES.md). They cover all
user-named anchors, usable follow-ups, rejected source/endpoint substitutions,
access failures, rights, search stopping and row-level verification.

### Material source limitations

* **Response lineage:** Mo's table is the numerical source; it is not claimed
  to reconstruct all original trials. Its stated lineage overlaps Yarchoan,
  so this is not an independent replication of the historical response
  compilation. The source's adjacent miscitation is retained in the audit.
  Reported trial denominators do not justify reconstructing responder counts
  or treating heterogeneous pooled outcomes as a homogeneous binomial sample.
* **TMB matching:** Chalmers is a different cohort and assay family from the
  response trials. Median TMB is not averaged across histologies. Biomarker
  strata and pooled esophagogastric outcomes with no defensible median remain
  visible but are not fitted. Named/NOS matches do not establish identical
  histology distributions, treatment lines or eligibility.
* **Clonal access:** every numerical clonal cell was rechecked against the
  complete indexed text of the primary author-hosted supplement. Direct
  retrieval currently fails; the live TCIA site reports maintenance. This
  is **not visual verification or a byte-pinned source version**. The
  estimand's per-patient weighting is not established by the table. Removing
  this source makes clonality/full-model analyses unavailable. No patient
  points or mixed patient-level workbooks were used to fill the gap.
* **HLA extraction:** only published aggregate bars and confidence-interval
  endpoints were digitized, with coordinates, calibration and conservative
  reading envelopes retained. Printed values take priority. The source's
  lung-squamous text/figure discrepancy is preserved as an alternative row,
  not counted as another cancer. Published marginal prevalence intervals
  remain distinct from reading precision. Their statistical construction
  was not assumed. No denominators were inferred from count-bar lengths.
* **Definition and cohort mismatch:** ordinary HLA LOH, HLA-A-only loss,
  complete presentation loss and an evolving loss hazard are different
  quantities. DASH, Montesion and TRACERx differ in stage, assay sensitivity,
  region sampling and eligibility. The mixed-source sensitivity cannot
  separate source identity from cancer type. Combined NSCLC is never given
  an arbitrary average of separate lung-histology predictors.

## Analysis and decision rules

| Model | Fitted form | Role |
| --- | --- | --- |
| Mean only | Constant ORR | Reference |
| Linear TMB | \(a+bT\) | Required simple baseline |
| Log-linear TMB | \(a+b\log T\) | Required primary comparator; natural logarithm |
| Reciprocal TMB | \(A\exp(-B/T)\) | Reduced shape test, with unmeasured constant clonality absorbed into B |
| Reciprocal clonal | \(A\exp[-B/(Tf)]\) | Fixed clonal-proxy extension on matched rows |
| Reciprocal TMB + HLA | \([A_0(1-L)+A_1L]\exp(-B/T)\) | HLA extension on matched rows |
| Full reciprocal | \([A_0(1-L)+A_1L]\exp[-B/(Tf)]\) | Full specification, gated; additional fits explicitly exploratory |

Fit on the ORR probability scale with equal weight per source cancer row,
including zeros. Linear baselines use OLS and clipped probability predictions.
The reciprocal models use bounded least squares: for each B, solve the
ceiling problem exactly, search a broad nonnegative B grid, refine every
grid-local minimum, and check the final fit against a denser grid. Boundary
solutions are retained and counted. The exponential is not fitted by taking
logs of observed zeros or adding pseudocounts.

**Prediction:** refit every parameter after leaving out a whole parent cancer.
Related lung categories share a fold. Score mean squared error and absolute
error within each parent, then give each parent equal weight. Matched subsets
get freshly fitted baselines; changing the subset is not an improvement.
Training R² is descriptive. The primary decision requires lower held-out MSE
than both simple TMB baselines and a favorable descriptive paired-error
interval against log-linear TMB. No p-values or multiple-comparison
significance claims are manufactured.

**Uncertainty:** group bootstrap intervals describe variation across this
selected cancer panel conditional on the source values. They do not capture
patient-level sampling, trial heterogeneity, assay error, unmeasured covariates
or a random sample of all cancers. Held-out losses share training data, so
their resampling intervals are not calibrated independent-test CIs. Source
rounding/reading scenarios and published HLA interval endpoints are separate
stress tests. A box of marginal source intervals is not a joint confidence
region. Parent-deletion exponent ranges are influence diagnostics, not CIs.

Every model, alternative source reading and specified sensitivity is retained.
The free-exponent models are diagnostics; they do not replace the failed fixed
form. Their numerical search bounds are fixed in the protocol. The full
coherent model's gate and the unavailable selection/exact-mapping sensitivities
remain visible rather than being relaxed after seeing results.

## Recomputed results

<!-- BEGIN RECOMPUTED RESULTS -->

### Coverage and preregistered gates

| Panel | Rows | Parent cancers | Required parents | Status |
| --- | ---: | ---: | ---: | --- |
| all | 21 | 20 | 3 | eligible |
| clonal | 11 | 11 | 5 | eligible |
| hla | 10 | 10 | 10 | eligible |
| full | 9 | 9 | 10 | **unavailable under frozen gate** |
| full_below_gate | 9 | 9 | unchanged | **exploratory** |
| full_mixed_sources | 10 | 10 | unchanged | **exploratory** |

The CSVs retain 73 source rows: clonal_fractions.csv: 19, hla_loh.csv: 29, tmb_orr.csv: 25.
TMB is missing for 4 original ORR rows. The source ORR table reports 6348 trial participants across its rows; this is a transcribed table total, not a verified unique-person count.

### Source inputs and coherent-panel joins

ORR and predictor fractions below are percentages; TMB is mutations/Mb. A dash means missing. HLA is the coherent DASH panel. Source-specific reading/rounding envelopes and provenance are in the CSVs.

| Cancer category | TMB | ORR % | Clonal % | HLA-LOH % | TMB category match |
| --- | ---: | ---: | ---: | ---: | --- |
| Anal carcinoma | 5.40 | 24.3 | — | — | approximate |
| Urothelial carcinoma | 7.20 | 18.1 | 66.54 | 23.00 | approximate |
| Breast carcinoma | 3.80 | 5.7 | 62.75 | 9.00 | exact |
| Cervical carcinoma | 5.40 | 14.9 | 49.96 | 39.00 | approximate |
| Microsatellite-instable colorectal cancer | — | 30.4 | — | — | missing |
| Microsatellite-stable colorectal cancer | — | 0.0 | — | — | missing |
| Cutaneous squamous-cell carcinoma | 45.20 | 52.2 | — | — | exact |
| Esophago-gastric carcinoma | — | 10.9 | — | — | missing |
| Glioblastoma | 2.70 | 8.7 | 44.88 | — | exact |
| Head and neck carcinoma | 5.00 | 14.6 | 58.57 | 40.00 | approximate |
| Renal-cell carcinoma | 2.70 | 25.1 | 48.46 | 16.00 | exact |
| Hepatocellular carcinoma | 3.60 | 18.1 | 55.70 | 4.00 | exact |
| Merkel-cell carcinoma | 4.30 | 56.0 | — | — | exact |
| Mesothelioma | 1.80 | 13.8 | — | — | approximate |
| Microsatellite-instable non-colorectal cancer | — | 37.7 | — | — | missing |
| Non-small-cell lung cancer | 8.10 | 19.8 | — | — | exact |
| Ovarian cancer | 3.60 | 9.9 | — | 14.00 | approximate |
| Pancreatic cancer | 2.50 | 0.0 | 52.08 | 26.00 | exact |
| Prostate carcinoma | 2.70 | 6.8 | 56.50 | — | approximate |
| Sarcoma | 2.50 | 9.0 | — | — | approximate |
| Small-cell lung cancer | 9.90 | 11.0 | — | — | exact |
| Skin melanoma | 14.40 | 37.2 | 60.86 | 7.00 | exact |
| Germ-cell cancer | 2.70 | 0.0 | — | — | approximate |
| Endometrial carcinoma | 4.50 | 13.3 | 60.59 | 7.00 | approximate |
| Uveal melanoma | 1.30 | 3.6 | — | — | exact |

### Held-out prediction

Each model is refitted after holding out an entire parent cancer. RMSE and MAE give equal weight to parents. Training R² uses row-level squared residuals and is descriptive. Models within a panel use identical rows; comparisons across different panels are not performance improvements.

| Panel | Model | Training R² | Held-out RMSE (pp) | Held-out MAE (pp) | Mean residual (pp) |
| --- | --- | ---: | ---: | ---: | ---: |
| all | Mean only | 0.0000 | 15.775 | 11.608 | 0.092 |
| all | Linear TMB | 0.3908 | 14.884 | 9.861 | -1.459 |
| all | Log-linear TMB | 0.4213 | 12.084 | 8.579 | 0.789 |
| all | Reciprocal TMB | 0.3801 | 13.173 | 9.183 | 1.445 |
| clonal | Mean only | 0.0000 | 10.670 | 7.907 | -0.000 |
| clonal | Linear TMB | 0.5942 | 7.321 | 5.193 | 0.449 |
| clonal | Log-linear TMB | 0.5269 | 8.792 | 6.715 | 0.749 |
| clonal | Reciprocal TMB | 0.5036 | 9.223 | 7.041 | 1.289 |
| clonal | Reciprocal clonal | 0.3986 | 10.349 | 8.059 | 1.374 |
| hla | Mean only | 0.0000 | 10.829 | 7.942 | -0.000 |
| hla | Linear TMB | 0.5648 | 8.153 | 6.062 | 0.814 |
| hla | Log-linear TMB | 0.4809 | 9.943 | 7.819 | 1.065 |
| hla | Reciprocal TMB | 0.4584 | 10.045 | 7.826 | 1.601 |
| hla | Reciprocal TMB + HLA | 0.5139 | 9.708 | 7.651 | 1.672 |
| full_below_gate | Mean only | 0.0000 | 11.329 | 8.292 | -0.000 |
| full_below_gate | Linear TMB | 0.5525 | 8.837 | 6.803 | 0.985 |
| full_below_gate | Log-linear TMB | 0.4668 | 10.703 | 8.700 | 1.200 |
| full_below_gate | Reciprocal TMB | 0.4392 | 10.787 | 8.814 | 1.711 |
| full_below_gate | Reciprocal clonal | 0.3307 | 12.040 | 9.802 | 1.653 |
| full_below_gate | Reciprocal TMB + HLA | 0.5030 | 10.361 | 8.340 | 1.777 |
| full_below_gate | Full reciprocal | 0.3832 | 11.950 | 9.783 | 1.941 |
| full_mixed_sources | Mean only | 0.0000 | 11.080 | 8.218 | -0.000 |
| full_mixed_sources | Linear TMB | 0.5786 | 7.917 | 5.865 | 0.588 |
| full_mixed_sources | Log-linear TMB | 0.5078 | 9.539 | 7.542 | 0.896 |
| full_mixed_sources | Reciprocal TMB | 0.4850 | 9.834 | 7.670 | 1.462 |
| full_mixed_sources | Reciprocal clonal | 0.3841 | 10.945 | 8.470 | 1.570 |
| full_mixed_sources | Reciprocal TMB + HLA | 0.5491 | 9.083 | 6.964 | 1.365 |
| full_mixed_sources | Full reciprocal | 0.4380 | 10.703 | 8.359 | 1.746 |

### Paired prediction-error differences

ΔMSE = designated theory model minus comparator, in squared percentage points; negative favors theory. Intervals are conditional descriptive percentile intervals from 10000 parent resamples, seed 20261117. They are not a formal significance test: folds share training data, cancers are not a random sample, and these are multiple descriptive comparisons.

| Panel | Comparator | ΔMSE (pp²) | Descriptive 95% interval (pp²) |
| --- | --- | ---: | --- |
| all | Log-linear TMB | 27.491 | [-0.144, 66.098] |
| all | Linear TMB | -48.003 | [-161.917, 29.539] |
| clonal | Log-linear TMB | 29.813 | [5.092, 66.631] |
| clonal | Linear TMB | 53.513 | [3.412, 125.515] |
| clonal | Reciprocal TMB | 22.048 | [4.984, 43.748] |
| hla | Log-linear TMB | -4.631 | [-24.913, 13.670] |
| hla | Linear TMB | 27.768 | [-10.263, 74.631] |
| hla | Reciprocal TMB | -6.662 | [-22.316, 9.802] |
| full_below_gate | Linear TMB | 64.701 | [2.963, 140.827] |
| full_below_gate | Log-linear TMB | 28.251 | [-5.496, 69.661] |
| full_below_gate | Reciprocal TMB | 26.446 | [-0.618, 55.378] |
| full_below_gate | Reciprocal clonal | -2.174 | [-20.310, 20.550] |
| full_below_gate | Reciprocal TMB + HLA | 35.438 | [11.356, 65.025] |
| full_mixed_sources | Linear TMB | 51.874 | [-6.098, 129.740] |
| full_mixed_sources | Log-linear TMB | 23.559 | [-5.730, 64.348] |
| full_mixed_sources | Reciprocal TMB | 17.847 | [-3.375, 41.907] |
| full_mixed_sources | Reciprocal clonal | -5.248 | [-17.458, 6.904] |
| full_mixed_sources | Reciprocal TMB + HLA | 32.054 | [6.675, 68.006] |

**Primary descriptive improvement criterion: not met.** The available primary panel is `all`; the preregistered full model remains below its complete-case gate.

### Fitted composite parameters and uncertainty

A₀ and A₁ are fitted endpoint ceilings, B is an uncalibrated composite. η is a population-mixture penalty, not the presentation-loss hazard. Intervals use 1000 cancer-group bootstrap draws, seed 20261117; they condition on the source values and do not include between-study measurement error or cure/response calibration.

| Panel | Parameter | Estimate | Descriptive 95% interval | Usable draws |
| --- | --- | ---: | --- | ---: |
| all | A0 | 0.4491 | [0.2170, 0.5998] | 1000 |
| all | B | 4.1992 | [1.5544, 7.3884] | 1000 |
| clonal | A0 | 0.3787 | [0.1347, 0.5996] | 1000 |
| clonal | B | 2.2681 | [0.0000, 4.4435] | 1000 |
| hla | A0 | 0.5026 | [0.1448, 0.6791] | 1000 |
| hla | B | 4.5241 | [0.0000, 8.2609] | 1000 |
| hla | A1 | 0.0507 | [0.0000, 0.6298] | 1000 |
| hla | eta | 0.8992 | [0.0000, 1.0000] | 1000 |
| full_below_gate | A0 | 0.4254 | [0.1429, 0.6579] | 1000 |
| full_below_gate | B | 2.0808 | [0.0000, 5.0684] | 1000 |
| full_below_gate | A1 | 0.0517 | [0.0000, 0.6336] | 1000 |
| full_below_gate | eta | 0.8786 | [0.0000, 1.0000] | 1000 |
| full_mixed_sources | A0 | 0.4520 | [0.1359, 0.6107] | 1000 |
| full_mixed_sources | B | 2.2863 | [0.0000, 4.3385] | 1000 |
| full_mixed_sources | A1 | 0.0000 | [0.0000, 0.5904] | 1000 |
| full_mixed_sources | eta | 1.0000 | [0.0000, 1.0000] | 1000 |

| Panel | Skipped singular draws | Fits on any active constraint boundary |
| --- | ---: | ---: |
| all | 0 | 1 |
| clonal | 0 | 67 |
| hla | 0 | 504 |
| full_below_gate | 0 | 573 |
| full_mixed_sources | 0 | 646 |

### HLA sign and exponent diagnostics

The HLA diagnostic removes A₁ ≤ A₀ while keeping both ceilings in [0,1]. Positive A₁ − A₀ (negative η) would reverse the adverse bridge. The relaxed intervals, rather than an enforced nonnegative sign, show uncertainty.

| Panel | A₁ − A₀ | Descriptive 95% interval | η | Descriptive 95% interval | Held-out RMSE (pp) |
| --- | ---: | --- | ---: | --- | ---: |
| hla | -0.4519 | [-0.6050, 0.2149] | 0.8992 | [-1.0453, 1.0000] | 9.708 |
| full_below_gate | -0.3737 | [-0.5675, 0.3918] | 0.8786 | [-1.7533, 1.0000] | 11.950 |
| full_mixed_sources | -0.4520 | [-0.5759, 0.3686] | 1.0000 | [-1.9844, 1.0000] | 10.703 |

| Shape diagnostic | Fitted exponent | Range after each parent deletion | Boundary deletions |
| --- | ---: | --- | ---: |
| TMB exponent q; specified q = 1 | 0.3989 | [0.3795, 0.9059] | 0 |
| Clonal exponent θ without HLA; specified θ = 1 | -1.6857 | [-2.0000, -0.8042] | 3 |
| Joint θ with mixed-source HLA, exploratory | -2.0000 | [-2.0000, -1.3565] | 9 |

Deletion ranges are influence diagnostics, **not confidence intervals**. A fitted boundary is not a precise magnitude estimate. The negative clonal exponent is descriptive contrary evidence to the positive ecological coverage proxy, not a causal claim that clonality is harmful.

### Every primary held-out prediction

| Cancer category | Observed ORR % | Reciprocal % | Log-linear % | Linear % | Reciprocal residual (pp) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Anal carcinoma | 24.300 | 20.376 | 19.430 | 15.579 | +3.924 |
| Urothelial carcinoma | 18.100 | 25.662 | 23.650 | 17.799 | -7.562 |
| Breast carcinoma | 5.700 | 15.626 | 15.775 | 14.870 | -9.926 |
| Cervical carcinoma | 14.900 | 21.041 | 19.932 | 16.058 | -6.141 |
| Cutaneous squamous-cell carcinoma | 52.200 | 30.917 | 40.639 | 89.729 | +21.283 |
| Glioblastoma | 8.700 | 9.552 | 11.161 | 13.558 | -0.852 |
| Head and neck carcinoma | 14.600 | 19.734 | 18.923 | 15.651 | -5.134 |
| Renal-cell carcinoma | 25.100 | 8.178 | 9.964 | 12.577 | +16.922 |
| Hepatocellular carcinoma | 18.100 | 13.632 | 14.408 | 13.972 | +4.468 |
| Merkel-cell carcinoma | 56.000 | 13.447 | 14.860 | 12.698 | +42.553 |
| Mesothelioma | 13.800 | 3.904 | 4.893 | 12.286 | +9.896 |
| Non-small-cell lung cancer | 19.800 | 30.017 | 26.940 | 19.195 | -10.217 |
| Ovarian cancer | 9.900 | 14.339 | 14.851 | 14.429 | -4.439 |
| Pancreatic cancer | 0.000 | 9.122 | 10.841 | 13.879 | -9.122 |
| Prostate carcinoma | 6.800 | 9.723 | 11.300 | 13.672 | -2.923 |
| Sarcoma | 9.000 | 8.320 | 10.113 | 13.331 | +0.680 |
| Small-cell lung cancer | 11.000 | 33.258 | 29.772 | 21.040 | -22.258 |
| Skin melanoma | 37.200 | 32.500 | 30.959 | 23.975 | +4.700 |
| Germ-cell cancer | 0.000 | 10.354 | 11.796 | 14.079 | -10.354 |
| Endometrial carcinoma | 13.300 | 17.990 | 17.597 | 15.191 | -4.690 |
| Uveal melanoma | 3.600 | 1.742 | 1.473 | 12.439 | +1.858 |

### Largest parent-cancer outliers

Rank by parent mean absolute held-out residual; preserve each row sign within a parent. Positive means the observed response exceeds prediction.

| Model | Rank | Parent / constituent category | Observed % | Predicted % | Residual (pp) |
| --- | ---: | --- | ---: | ---: | ---: |
| Reciprocal TMB | 1 | merkel: Merkel-cell carcinoma | 56.000 | 13.447 | +42.553 |
| Reciprocal TMB | 2 | cutaneous_scc: Cutaneous squamous-cell carcinoma | 52.200 | 30.917 | +21.283 |
| Reciprocal TMB | 3 | renal_cell: Renal-cell carcinoma | 25.100 | 8.178 | +16.922 |
| Reciprocal TMB | 4 | lung: Non-small-cell lung cancer | 19.800 | 30.017 | -10.217 |
| Reciprocal TMB | 4 | lung: Small-cell lung cancer | 11.000 | 33.258 | -22.258 |
| Reciprocal TMB | 5 | germ_cell: Germ-cell cancer | 0.000 | 10.354 | -10.354 |
| Log-linear TMB | 1 | merkel: Merkel-cell carcinoma | 56.000 | 14.860 | +41.140 |
| Log-linear TMB | 2 | renal_cell: Renal-cell carcinoma | 25.100 | 9.964 | +15.136 |
| Log-linear TMB | 3 | lung: Non-small-cell lung cancer | 19.800 | 26.940 | -7.140 |
| Log-linear TMB | 3 | lung: Small-cell lung cancer | 11.000 | 29.772 | -18.772 |
| Log-linear TMB | 4 | germ_cell: Germ-cell cancer | 0.000 | 11.796 | -11.796 |
| Log-linear TMB | 5 | cutaneous_scc: Cutaneous squamous-cell carcinoma | 52.200 | 40.639 | +11.561 |

### Sensitivities and failed falsification attempts

All differences below compare the designated theory with log-linear TMB on identical rows. Scenario ranges are neither CIs nor exhaustive bounds.

| Panel / sensitivity | Parent cancers | Theory RMSE (pp) | Log-linear RMSE (pp) | ΔMSE (pp²) |
| --- | ---: | ---: | ---: | ---: |
| primary: exact mappings only | 10 | 18.494 | 16.387 | 73.505 |
| coherent_full: exact mappings only | 1 | unavailable | unavailable | unavailable |
| coherent_full: HLA source intervals all_lower | 9 | 11.862 | 10.703 | 26.152 |
| coherent_full: HLA source intervals all_upper | 9 | 11.985 | 10.703 | 29.100 |
| coherent_full: figure instead of printed HLA | 9 | 11.923 | 10.703 | 27.621 |

| Panel | Extraction scenarios | ΔMSE scenario range (pp²) | Parent-deletion ΔMSE range (pp²) | Deletions favoring theory |
| --- | ---: | --- | --- | ---: |
| primary | 100 | [26.749, 28.244] | [5.692, 34.617] | 0/20 |
| coherent_full | 100 | [25.462, 31.456] | [-14.894, 36.022] | 2/9 |

Excluding indexed-only clonal evidence leaves 0 usable source clonal categories: clonal/full tests are unavailable. Trial selection is unknown for 21 included ORR rows, so no unselected-only sensitivity is claimed.

### Normalized provenance fingerprints

| File | SHA-256 of LF-normalized content |
| --- | --- |
| data/clonal_fractions.csv | `242b0f239462b52b3e70b0e726dca996ea9e322735e1c2b4de6c21b458e669f8` |
| data/hla_loh.csv | `c3e44d0440ee5430d769dc340b1e037fb9fc8f2ed847100588f106d6279d166e` |
| data/tmb_orr.csv | `b34c1dd8d8d2c403f7251cf02611283d5326cbdc2caf2c72a6afadcf5ddc540b` |
| PROTOCOL.md | `65c9cd84cf10b26f1f6236adc0f3d9b9d622d210dd5376ba98cdc9ed41f36ab3` |

<!-- END RECOMPUTED RESULTS -->

## Verdicts and falsification record

| Question / attempted falsification | Verdict and interpretation |
| --- | --- |
| Does the predeclared reciprocal shape outperform TMB alone? | **No demonstrated predictive gain.** It does not outperform the log-linear primary comparator and fails the declared improvement criterion. The paired interval crosses zero; do not call the difference a calibrated significant rejection. |
| Does the fixed positive clonal correction help? | **Contrary descriptive evidence.** The clonal correction increases held-out error against both simple baselines and the matched reciprocal-TMB ablation. The free clonal exponent reverses the proposed direction, including the parent-deletion diagnostics. This contradicts the cancer-type proxy assumption in this panel, not the biological benefit of clonal targets within a tumour. |
| Do HLA sign and magnitude match? | **Point direction compatible; magnitude unestablished.** The adverse direction remains when the sign is relaxed, but the relaxed uncertainty includes the opposite direction. A fitted maximum penalty is a boundary/extrapolation, not evidence that ordinary HLA LOH eliminates response. The theorem predicts no ordinary-LOH coefficient. |
| Is the complete preregistered model tested? | **Unavailable under its frozen sample gate.** The coherent below-gate and mixed-source full fits are explicitly exploratory and do not improve on the required baselines. Additional heterogeneous source coverage does not restore a preregistered test. |
| Can dropping an inconvenient cancer rescue the primary claim? | Every parent deletion is reported; the primary reciprocal-vs-log-linear disadvantage retains its direction. The exploratory coherent full comparison is more fragile, changing direction under some deletions. No cancer is dropped from the main fit. |
| Is the result a digitization/rounding artifact? | The specified extraction perturbations and printed-versus-figure HLA comparison retain the main direction. These finite scenarios do not bound shared-source systematic error. Removing the indexed-only clonal source prevents the clonal/full test altogether. |
| Does using exact category matches establish the result? | The reduced exact-match sensitivity still disfavors the reciprocal model. Exact matching leaves too few full-model categories to analyze. Histology proxies therefore remain a material limitation. |
| Does this test the sharp coefficient or N/k scaling itself? | **No.** Functional load, effective founder number, calibrated eradication/response mapping and rates are absent. Neither B nor HLA prevalence measures the theorem's mechanistic coefficient. |

The fitted HLA ceilings concern hypothetical all-LOH versus no-LOH populations,
outside much of the observed prevalence range. A large eta can therefore be
weakly identified even when predictions are stable. Constrained intervals
touching the entire permitted penalty range are not evidence for a precise
penalty. The relaxed ceiling difference is shown to avoid relying only on
a ratio that becomes unstable near a zero ceiling.

The free TMB-exponent fit also reaches the upper endpoint-ceiling constraint.
Its exponent and ceiling are coupled, so the exponent profile is a shape
diagnostic, not a separately calibrated mechanistic scaling estimate.

The negative clonal exponent is especially useful as a falsification attempt
because it was not prohibited by the diagnostic. It could reflect assay and
stage differences, restricted histology coverage, antigen quality, immune
infiltration, mutation-class differences, or between-cancer confounding.
These explanations are alternatives to test, not licenses to relabel a
contradiction as support. A positive within-cancer effect can coexist with
a negative between-cancer coefficient: **ecological associations cannot
determine individual treatment response or causal benefit**.

## Largest outliers: minimal post hoc explanations

The ranking is generated above from held-out residuals and was not used to
exclude data or tune the model. The same parent cancers dominate both the
theory and log-linear baseline, in different orders. The following are
mechanistic hypotheses or sampling explanations, not established causes.

| Outlier | Minimal explanation consistent with existing evidence | Smallest useful distinguishing aggregate |
| --- | --- | --- |
| Merkel-cell carcinoma | Mutation-derived neoantigens omit viral targets. [Nghiem et al. 2016](https://doi.org/10.1056/NEJMoa1603702) explicitly discusses immune recognition of polyomavirus-associated tumours despite low mutation burden. | Published ORR and median TMB by viral status, plus an aggregate viral-antigen expression/presentation measure. Do not infer individual viral status from TMB. |
| Cutaneous squamous-cell carcinoma | This holdout extrapolates to the highest TMB region; the shared fitted ceiling is restrictive. Immune competence and trial selection are possible ceiling modifiers. [Migden et al. 2018](https://doi.org/10.1056/NEJMoa1805131) describes UV hypermutation and the immune-surveillance/eligibility context. TMB itself is already in the model and cannot explain its own residual. | An independent published high-TMB response panel, with comparable eligibility and treatment line and aggregates stratified by immune suppression. Fix any new ceiling rule before inspecting those responses. |
| Renal-cell carcinoma | Mutation classes can differ in the number and quality of targets they generate. [Turajlic et al. 2017](https://doi.org/10.1016/S1470-2045(17)30516-8) reports indel enrichment and predicted frameshift-derived immunogenicity in renal cancers. | Type/stratum-level frameshift-indel or expressed-neoantigen burden alongside TMB and ORR, with comparable histology. This tests effective antigen yield rather than post hoc adjustment of f. |
| Lung, especially small-cell | High mutation burden can coexist with impaired MHC-I presentation. [Mahadevan et al. 2021](https://doi.org/10.1158/2159-8290.CD-20-0913) distinguishes MHC-I-high and low neuroendocrine states and ICB benefit in small-cell lung cancer. This explanation is specific to SCLC, not automatically all NSCLC. | Subtype-specific aggregate MHC-I-high prevalence, TMB and ORR, measured in comparable treatment cohorts. Keep NSCLC and SCLC endpoints visible inside the parent holdout. |
| Germ-cell cancer | A small heavily pretreated nonseminomatous cohort can observe no responses without requiring a mechanistic change to a positive population prediction. [Adra et al. 2018](https://doi.org/10.1093/annonc/mdx680) documents that trial context. Sampling and cohort mismatch are sufficient competing explanations. | Independent published responder numerators and denominators separated by seminoma/nonseminoma and prior therapy, with matched antigen/TMB summaries. No causal immune-failure mechanism is inferred from a sample zero. |

The minimal general repair would let functional target yield, presentation
and immune effectiveness vary between cancers. That is a hypothesis for a
new independently specified study. Allowing every cancer its own free
coefficient now would erase falsifiability and cannot count as a successful
test of the shared-parameter form.

## Novelty assessment

This is a repository-specific confrontation of PR59's proposed surrogate,
**not the first ecological TMB–response analysis, multivariable predictor,
clonality biomarker, HLA-loss analysis or use of cross-validation**.

| Prior work | What was already done | What this note can claim |
| --- | --- | --- |
| [Yarchoan et al. 2017](https://doi.org/10.1056/NEJMc1713444) | Ecological association and log-linear prediction of checkpoint response from median TMB, including known cancer-specific exceptions. | A different, predeclared reciprocal functional form is compared on a later published aggregate table sharing that response lineage. This is not independent outcome replication. |
| [Lee and Ruppin 2019](https://doi.org/10.1001/jamaoncol.2019.2311) | Cross-cancer immune/neoantigen/checkpoint predictors and leave-one-out prediction beyond TMB. | Grouped prediction, source/gate transparency and a negative result for this specific constrained curve; no priority claim for multivariable ecological prediction or cross-validation. |
| [Charoentong et al. 2017](https://doi.org/10.1016/j.celrep.2016.12.019), [McGranahan et al. 2016](https://doi.org/10.1126/science.aaf1490), [Litchfield et al. 2021](https://doi.org/10.1016/j.cell.2021.01.002) | Pan-cancer immunogenomics, clonal neoantigen biology and response-associated clonal burden. | Reuse of existing type-level fractions to challenge an explicit observation bridge, not discovery of those data or the biological clonality association. |
| [Montesion et al. 2021](https://doi.org/10.1158/2159-8290.CD-20-0672) | HLA loss and refinement of TMB-based response interpretation. | An explicitly restricted prevalence mixture and an uncertainty/sign audit; HLA prevalence is never relabeled a transition hazard. |
| [Zhang et al. 2025](https://doi.org/10.17305/bb.2024.11069) | Ecological response prediction combining TMB, HLA expression and immune infiltration. | No claim to be the first antigen-presentation-aware ecological model. Expression and loss prevalence are different predictors; its source table was not used to repair this fit. |

The useful contribution is the auditable specification, aggregate extraction,
held-out comparison, preserved contradiction and explicit missing-test record.
The probability law belongs to PR59 under its stated mathematical assumptions;
no new theorem is claimed here. Publication-level novelty remains unestablished.
This targeted source audit is not an exhaustive systematic review.

## Reproduce and audit

From the repository root, using only the Python standard library:

```sh
python -B research/immune-escape-pancancer-test/analysis.py --check
python -B research/immune-escape-pancancer-test/analysis.py --json
python -B -m unittest discover -s tests -p 'test_research_immune_escape_pancancer_test.py' -v
```

The JSON command prints every fit, fold, parameter, bootstrap interval,
exclusion and sensitivity for inspection; no external fetch is required.
`analysis.py --write` regenerates only this README's marked numeric block.
The dedicated test recomputes the report from the CSVs, independently checks
linear fits, bounded nonlinear solutions, held-out predictions, error metrics,
bootstrap arithmetic, source conversions and joins, gate behavior and the
generated block. File hashes use newline-normalized content. No package
installation, network access, patient records or clinical credentials are
needed to rerun this analysis.

Repository integration checks, with the repository's existing validation
dependencies already available, are:

```sh
python scripts/validate.py
python -B -m unittest discover -s tests -v
node --test tests/*.test.*
```

### Hostile review and boundaries of verification

All ORR entries were checked against the rendered original aggregate PDF
and its extracted text. Every retained TMB median, source count and label
was compared with its recorded original workbook cell. Every clonal value
was re-read from the indexed primary table, with its complementary percentage
and access limitation preserved. The HLA printed text, republished table
cells and aggregate figure readings were checked separately; inconsistent
printed/figure values were not silently repaired.

Code review challenged constrained sign claims, singular coefficient designs,
feature changes during deletion, unmatched comparisons and boundary counts.
The implementation preserves zero outcomes, separates missing data, fixes
features across deletion, refits held-out folds, counts constraint boundaries,
and reports coefficient uncertainty without requiring the desired sign.
Source review and software tests have different meanings: neither proves
that a source's biological measurements or cross-cohort assumptions are true.

Only this new directory and
`tests/test_research_immune_escape_pancancer_test.py` are added. No existing
research note, index, catalog, workflow, requirement, result or script is
modified. The existing PR workflow discovers the new test. This research
remains nonclinical and provides no treatment advice.
