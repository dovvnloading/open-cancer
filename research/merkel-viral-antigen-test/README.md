# Merkel-cell response: frozen TMB curve, assay scale and viral antigens

**Nonclinical, unreviewed computational research. No treatment, diagnosis,
surveillance or patient-selection advice.** Published aggregate summaries
only. Code checks, cryptographic hashes and AI critique do not establish
biological truth, clinical validity or external scientific review.

## Verdict

**Untestable as an attribution claim:** these published marginals cannot
establish that viral antigen supply or assay scale explains the original
Merkel-cell residual. The restricted surrogate itself is testable, and its
failures and conditional repairs are informative:

* The frozen reciprocal curve with **median substitution** fails both
  predeclared assay bands. There is nevertheless a common finite-positive
  scale region outside those bands that satisfies both cohorts. Its large
  magnitude has no independent calibration here, so this is mathematical
  feasibility rather than an assay-scale explanation.
* A single nonnegative viral-load offset fixed on Part A improves the
  status-balanced Part B point score at every prespecified scale. At the
  identity reference it still fails the complete repair rule because the
  unchanged virus-negative control falls below its published interval.
  Other fixed scales pass the descriptive rule. This conditional transport
  result supports studying non-mutational antigen supply; it does not
  identify the cause of the improvement.
* **Median substitution is a decisive limit on attribution.**
  Allowing arbitrary TMB distributions with the reported medians removes
  incompatibility throughout both assumed bands for the primary reciprocal
  curve. This broad bound does not show that the real distributions fit;
  it shows that the median summaries alone do not force a viral term.
* The denominator-matched Part A-to-Part B viral test is **unavailable**.
  The inspected Part B report gives TMB by viral status and TMB by response
  separately, without the needed viral-status ORR within its TMB subset.
  No individual points or margins were combined to manufacture that result.

At the reference scale, the declared viral repair therefore remains
**unexplained under its complete interval rule**, despite a better point
score. The top-level attribution verdict remains **untestable**, rather than
“explained by viral antigen supply” or “explained by assay scale.”

## What was frozen, and when

This executes ranked job 2 in the unchanged [research programme](../README.md),
following its negative-results guidance. The branch starts at
[`16c7596`](https://github.com/dovvnloading/open-cancer/commit/16c7596637a60836056593c1fb22b1f5c9d07501).

| Stage | Immutable record | Timing and meaning |
| --- | --- | --- |
| Question, curve identities, conversion family, scores, sensitivity rules and verdict rules | [PROTOCOL.md](PROTOCOL.md); original local commit `09ecbb18e6b159d8cf74d30eac82059dc54e472e`; [published identical-tree replica `b5d8eb4`](https://github.com/dovvnloading/open-cancer/commit/b5d8eb4e19e91f1decdcf258dda8a7a022f06d6f) | Locally committed at 2026-11-18 00:37:27 America/New_York, before any new fit, curve evaluation, conversion region or score. Source methods and published outcome summaries had been read to establish feasibility. |
| Part A coefficients and viral offsets | [part_a_calibration.csv](data/part_a_calibration.csv); original local commit `1d4f3682aa1d42b601a3ecac669f40b34ded65a4`; [published identical-tree replica `4c6d018`](https://github.com/dovvnloading/open-cancer/commit/4c6d0182d728671c753a1a935ead928d3bb21699) | Locally committed at 2026-11-18 00:40:02 America/New_York, before any Part B prediction or score. This commit includes Part A inputs and the analysis implementation that reads only Part A for calibration; it contains no Part B input file. |
| Part B confrontation and source cross-audit | This note, code and dedicated test | Performed after the calibration commit. Neither frozen artifact was edited. All predeclared scales and sensitivities are retained. |

This is a prospective computation plan applied to already published, partly
seen outcomes. It is not blinded validation or registration in an independent
registry. There was no outcome-driven protocol amendment. The exact base code
and source CSV define the curves; printed rounded coefficients were never used
as exact calibration inputs.

Publication used the authenticated GitHub connector after ordinary HTTPS Git
push failed because the local CLI had no credentials. The connector recreated
the ordered commits with new identifiers and publication timestamps. Both
freeze trees were verified identical to their original local Git trees:
`e54ad44bbfedecd114febe3a3220f44e4e187197` for the protocol and
`a69feb57696d45efe7474eb77e02303028c8e10d` for calibration. The original local
commit objects remain in the working repository; the published commit messages
also record their identifiers and dates. These are local computational freezes
published later, not independent public timestamps preceding computation.
Neither frozen file was changed during publication.

The post-scoring portability review added [frozen_curves.csv](data/frozen_curves.csv)
as a derived snapshot of the already specified original fits. This is an
implementation correction, with no JAVELIN refitting, changed scale or changed
decision rule. Replays must match the recorded coefficients within an absolute
numerical tolerance of `2e-6`, after exact canonical hashes verify the original
code and CSV. Prediction uses the recorded full-precision coefficients and
pre-Part-B viral offsets. Calibration replay uses numerical tolerance only
for floating-point fields (`5e-12` relative/absolute); its immutable file hash
and all provenance fields remain exact. This avoids assuming identical
platform math-library rounding at a flat optimizer minimum. Manifest order
is sorted by canonical POSIX strings, not filesystem-specific path ordering.

The primary curve is #67's **original leave-Merkel-parent-out reciprocal fit**,
which generated its largest held-out residual. Its full-panel fit and the
corresponding clipped log-TMB curves are required sensitivities. They are
replayed from the unchanged [original code](../immune-escape-pancancer-test/analysis.py)
and [TMB/ORR CSV](../immune-escape-pancancer-test/data/tmb_orr.csv), without
JAVELIN fitting of their coefficients or selection among their results.

## Sources, reuse and extraction audit

Access date for the new published aggregates: **2026-11-18,
America/New_York**. Every data row retains its own source, locator, population,
cutoff, units, inspection method and reuse basis. Transformations are factual
aggregate transcription, categorical labels and percent-to-probability
conversion in code. Derived calibration rows are explicitly labelled derived.

| Input | Primary evidence and locator | Reuse and inspection |
| --- | --- | --- |
| [part_a.csv](data/part_a.csv) | [D'Angelo et al., JITC 2020, e000674](https://pmc.ncbi.nlm.nih.gov/articles/PMC7239697/), Results “Antitumor activity”; Table 1; Results “Exploratory biomarker analyses”; Figs. 3B and 4; Supplementary Table S1 for unavailable viral status | Article CC BY 4.0. Main paper text and aggregate figures inspected, with a separate cross-audit of every populated numeric cell. Supplementary Table S1 was inspected as indexed publisher text; no direct visual audit of that supplement is claimed. |
| [part_b.csv](data/part_b.csv) | [D'Angelo et al., JITC 2021, e002646](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311489/), Table 2; Results “Efficacy” and “Tumor mutational burden”; Fig. 1D; Supplementary Fig. 3B | Article CC BY-NC 4.0. Main paper and Fig. 1D inspected; every numerical field cross-audited. Supplement labels accessed as indexed publisher text, with medians and TMB denominators independently printed in the accessible main Results. |
| Frozen mutation-scale context | [Chalmers et al., Genome Medicine 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5395719/), Methods “Tumor mutational burden”; existing #67 CSV and [source audit](../immune-escape-pancancer-test/TMB_ORR_SOURCES.md) | Existing attributed aggregate input used read-only. Chalmers counts synonymous changes with stated filters; no exact assay-to-assay calibration is supplied by that definition. |

The original papers and figures are not redistributed. Source licences remain
attached to the observations: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
for Part A and [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)
for Part B. This is noncommercial research use; repository MIT does not
relicense third-party data.

### What hostile source checking found

The independent cross-audit found **no numerical transcription correction**.
The primary Part A figure labels were visually rechecked. Part B's author-hosted
PDF has a low-resolution figure, so its smaller subgroup labels were additionally
corroborated in the indexed [primary PMC PDF](https://rcastoragev2.blob.core.windows.net/c00ff4d1f8a3cfecc6f0026ec61eccb6/PMC8311489.pdf).
The publisher's [Part B supplement](https://jitc.bmj.com/content/jitc/suppl/2021/07/23/jitc-2021-002646.DC1/jitc-2021-002646supp002_data_supplement.pdf)
returned HTTP 403 for direct retrieval and screenshots. Its indexed Fig. 3B
labels and printed page locator were checked, and all TMB numbers used here
also appear in the main Results. No digitized patient scatter points enter
the data. The indexed supplementary patient-level HLA table was not extracted
or used to infer the missing viral-by-response table.

The frozen Part A row `a_count_viral_unevaluable` uses `not_stratified` in its
categorical field; its identifier, population and metric explicitly mean the
published count with unavailable viral status. It is excluded from model
scoring. Its inspection phrase means indexed-text inspection. These metadata
qualifications are recorded here without rewriting the frozen calibration
input. They change no numerical field or analysis.

### Population and endpoint restrictions

Part A and Part B are distinct participant cohorts, previously treated and
first-line for metastatic disease, respectively. They share a trial programme,
agent, assay context and analysis infrastructure. They supply a cohort transfer
check, with no claim of independent research-programme replication. Earlier
subsets and later follow-ups are overlapping reports, not further cohorts.
The historical #67 panel has Mo/Yarchoan response lineage and cross-cohort
Chalmers medians; it is not reconstructed into JAVELIN data. Removing the Merkel
parent from the primary fit avoids directly training on that cancer's original
response, while historical overlap elsewhere is not presumed absent.

Part B explicitly identifies IHC for the viral groups used in these Results.
Part A describes both PCR and IHC, without a clear row-specific combination
rule in the report. We preserve the published groups and this uncertainty;
we do not equate the methods by assumption. Part A's prior-treatment margins
refer to systemic treatments at any disease stage; they cannot be crossed
with the viral margins. First-line Part B means first-line metastatic therapy.

ORR denominators include the source-defined all-treated populations, including
participants without an evaluable best response. TMB denominators describe
the smaller sequencing subsets. TMB and response dates differ within Part A;
baseline TMB is a measurement associated with the report cutoff, not a
measurement taken on that date. No source denominators were inferred from
rounded rates, and no median of separate medians was pooled.

## Methods and interpretation

For each frozen reciprocal curve,

\[
f(x)=A\exp(-B/x),\qquad x=c\,m,
\]

where `m` is a JAVELIN NSSV/Mb median and a single finite positive `c` is
shared across cohorts and viral statuses. The admissible family is a
multiplicative sensitivity bridge. Assumed bands and fixed score scenarios
are stated in [PROTOCOL.md](PROTOCOL.md); none is an estimated calibration.
Synonymous counting alone cannot guarantee the direction or magnitude of a
panel-versus-exome conversion after territory and filtering differences.

For a published response interval `[L,U]`, invert the increasing curve to
obtain every compatible `c`. When an endpoint lies strictly between zero and
the reciprocal ceiling, its inverse is `B/(m*ln(A/p))`. Boundary handling
retains unattainable ceiling and zero limits. Intersect both status sets
within each cohort, then both cohorts. Log-TMB fits use their corresponding
analytic inverse with clipping. These interval intersections are deterministic
descriptions of published marginal uncertainty, not joint statistical tests.

The one-parameter extension is

\[
p(m,z)=A\exp\!\left[-\frac{B}{c\,m+v z}\right],\qquad v\ge0,
\]

with `z` indicating published MCPyV-positive status. At each prespecified
fixed `c`, Part A's positive-stratum ORR fixes `v`; neither Part B nor the
negative stratum enters that fit. Zero and unattainable boundaries are handled
explicitly. A status-balanced mean squared error scores Part B's two published
ORR point estimates. Smaller loss and inclusion of both predictions in their
published intervals are both required for a numerical repair. The reported
CI-endpoint calibration sensitivities do not provide predictive intervals.

### Why the median-envelope sensitivity matters

For any nonnegative TMB distribution with median `m`, monotonicity and the
reciprocal response ceiling give the conservative closure

\[
\frac{f(cm)}{2}\ \le\ E[f(cT)]\ \le\ \frac{A+f(cm)}{2}.
\]

At least half the mass is at or above the median, giving the lower bound;
at least half is at or below it, giving the upper bound. Values at unattained
limits are treated as closures rather than realized finite-TMB solutions.
The analysis tests whether these envelopes overlap the published ORR intervals.
This is deliberately weak: observed sample medians are being treated as
population constraints, subset representativeness is assumed, and no sampling
uncertainty or distribution shape is identified. Its compatibility shows why
a nonlinear plug-in failure alone cannot select a biological mechanism.
The envelope permits arbitrarily high upper tails and does not impose the
observed sample ranges as population support. It is a permissive robustness
bound, not evidence that any actual cohort has the required distribution.

## Recomputed results

All empirical numeric results in the block below are regenerated from the
included CSVs and unchanged #67 CSV. Positive and negative groups receive
equal score weight; the source denominators are preserved for interpretation,
not substituted into a joint binomial likelihood.

<!-- BEGIN RECOMPUTED RESULTS -->

### Published marginal inputs

ORR is percent [published 95% CI]; TMB is median NSSV/Mb. Adjacent cells are an ecological pairing, not individual records. Unknown responder counts remain unreported.

| Cohort / viral stratum | ORR % [95% CI] | ORR denominator | Published responders | TMB median | TMB denominator |
| --- | --- | ---: | ---: | ---: | ---: |
| A / all | 33.0 [23.3, 43.8] | 88 | 29 | 0.58 | 36 |
| A / positive | 28.3 [16.0, 43.5] | 46 | unreported | 0.49 | 23 |
| A / negative | 35.5 [19.2, 54.6] | 31 | unreported | 2.72 | 13 |
| B / all | 39.7 [30.7, 49.2] | 116 | 46 | 0.34 | 52 |
| B / positive | 34.3 [23.3, 46.6] | 70 | 24 | 0.22 | 31 |
| B / negative | 48.6 [31.9, 65.6] | 37 | 18 | 10.52 | 19 |
| B / unknown | 44.4 [13.7, 78.8] | 9 | 4 | missing | missing |

Additional published margins remain separate:

| Part A context | ORR % [95% CI] | Denominator | Published responders |
| --- | --- | ---: | ---: |
| One prior systemic anticancer treatment in any disease stage | 40.4 [27.0, 54.9] | 52 | unreported |
| At least two prior systemic anticancer treatments in any disease stage | 22.2 [10.1, 39.2] | 36 | unreported |
| All TMB-evaluable participants | 33.3 [18.6, 51.0] | 36 | 12 |
| MCPyV-positive TMB-evaluable participants | 26.1 [10.2, 48.4] | 23 | 6 |
| MCPyV-negative TMB-evaluable participants | 46.2 [19.2, 74.9] | 13 | 6 |

Part A explicitly reports 11 participants with unevaluable viral status and 18 with unevaluable best response; these categories must not be equated.

| Cohort / endpoint | Source cutoff |
| --- | --- |
| A broad ORR | 2018-09-14 |
| A TMB and TMB-subset ORR | 2017-09-24 |
| B ORR and TMB | 2019-05-02 |

Transcribed source rows: Part A 13; Part B 9. The Part B count includes the explicitly missing viral-specific TMB-subset ORR rows.

### Exact frozen-curve replay

| Curve | Training rows / parents | A or a | B or b |
| --- | --- | ---: | ---: |
| Reciprocal, Merkel held out (primary) | 20 / 19 | 0.48718323459062135 | 5.5353933555249002 |
| Reciprocal, full panel | 21 / 20 | 0.44912947779826107 | 4.199233915543866 |
| Log TMB, Merkel held out | 20 / 19 | -0.035761939452094771 | 0.12639822421480662 |
| Log TMB, full panel | 21 / 20 | -0.014484116854756118 | 0.12526723878661711 |

Replayed historical Merkel input: median TMB 4.30, ORR 56.000%; held-out reciprocal prediction 13.447%, residual +42.553 pp. This is the original mixed-source motivating residual, not a JAVELIN viral contrast.

### Common-scale compatibility

All intervals below are sets of finite positive conversion factors c. They use printed median substitution and published ORR intervals. They are not calibrated confidence regions. Infinity denotes an open limit.

| Curve | Cohort | MCPyV+ scales | MCPyV- scales | Both statuses |
| --- | --- | --- | --- | --- |
| Reciprocal, Merkel held out (primary) | A | [10.146, 99.711] | [2.186, infinity) | [10.146, 99.711] |
| Reciprocal, Merkel held out (primary) | B | [34.112, 565.990] | [1.243, infinity) | [34.112, 565.990] |
| Reciprocal, full panel | A | [8.303, 268.100] | [1.817, infinity) | [8.303, 268.100] |
| Reciprocal, full panel | B | [29.085, infinity) | [1.167, infinity) | [29.085, infinity) |
| Log TMB, Merkel held out | A | [9.603, 84.587] | [2.228, 36.671] | [9.603, 36.671] |
| Log TMB, Merkel held out | B | [38.109, 240.764] | [1.574, 22.637] | empty |
| Log TMB, full panel | A | [8.217, 73.814] | [1.911, 32.255] | [8.217, 32.255] |
| Log TMB, full panel | B | [32.779, 210.566] | [1.362, 20.069] | empty |

| Curve | Both cohorts / both statuses | Within assumed [1, 2] | Within assumed [0.5, 5] |
| --- | --- | --- | --- |
| Reciprocal, Merkel held out (primary) | [34.112, 99.711] | empty | empty |
| Reciprocal, full panel | [29.085, 268.100] | empty | empty |
| Log TMB, Merkel held out | empty | empty | empty |
| Log TMB, full panel | empty | empty | empty |

Rounding expands medians and interval endpoints by half the last printed unit; it does not account for sampling error in the medians:

| Curve | Rounded A intersection | Rounded B intersection | Rounded common intersection |
| --- | --- | --- | --- |
| Reciprocal, Merkel held out (primary) | [10.015, 101.771] | [33.257, 593.468] | [33.257, 101.771] |
| Reciprocal, full panel | [8.194, 280.961] | [28.345, infinity) | [28.345, 280.961] |
| Log TMB, Merkel held out | [9.469, 36.884] | empty | empty |
| Log TMB, full panel | [8.102, 32.444] | empty | empty |

Substituting the published Part A TMB-subset ORRs (secondary) gives:

| Curve | A subset MCPyV+ | A subset MCPyV- | A subset intersection | With broad B intersection |
| --- | --- | --- | --- | --- |
| Reciprocal, Merkel held out (primary) | [7.225, 1723.270] | [2.186, infinity) | [7.225, 1723.270] | [34.112, 565.990] |
| Reciprocal, full panel | [5.781, infinity) | [1.817, infinity) | [5.781, infinity) | [29.085, infinity) |
| Log TMB, Merkel held out | [6.069, 124.642] | [2.228, 182.731] | [6.069, 124.642] | empty |
| Log TMB, full panel | [5.172, 109.149] | [1.911, 163.075] | [5.172, 109.149] | empty |

**The matched Part A-to-Part B test is untestable:** viral-specific TMB-subset ORRs in Part B are unpublished in the inspected sources. The last column above still uses broad B and therefore retains the denominator mismatch.

### Relaxing median substitution

For the primary reciprocal curve, replace f(c*m) by the conservative mean-response closure [f(c*m)/2, (A+f(c*m))/2]. This assumes the TMB median represents the ORR population; it does not identify its distribution or remove selection bias.

| Cohort | MCPyV+ scales | MCPyV- scales | Both statuses |
| --- | --- | --- | --- |
| A | (0, infinity) | (0, infinity) | (0, infinity) |
| B | (0, infinity) | [0.449, infinity) | [0.449, infinity) |

Common envelope-compatible scales: **[0.449, infinity)**. Assumed-band intersections: [1, 2] -> [1.000, 2.000]; [0.5, 5] -> [0.500, 5.000].

### Predictions of every frozen curve at every fixed scale

Percent predictions; Part B RMSE is status-balanced, in percentage points. Misses count predictions outside the published marginal intervals. Very small predictions use scientific notation to avoid implying an exact zero.

| Curve | c | A+ % | A- % | B+ % | B- % | B RMSE pp | B misses |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Reciprocal, Merkel held out (primary) | 0.5 | 7.507e-09 | 0.832 | 6.811e-21 | 17.008 | 32.974 | 2 |
| Reciprocal, Merkel held out (primary) | 1 | 0.001 | 6.366 | 5.761e-10 | 28.786 | 28.010 | 2 |
| Reciprocal, Merkel held out (primary) | 4/3 | 0.010 | 10.588 | 3.106e-07 | 32.833 | 26.694 | 1 |
| Reciprocal, Merkel held out (primary) | 2 | 0.172 | 17.611 | 1.675e-04 | 37.448 | 25.503 | 1 |
| Reciprocal, Merkel held out (primary) | 5 | 5.087 | 32.429 | 0.318 | 43.852 | 24.262 | 1 |
| Reciprocal, full panel | 0.5 | 1.617e-06 | 2.048 | 1.184e-15 | 20.214 | 31.482 | 2 |
| Reciprocal, full panel | 1 | 0.009 | 9.592 | 2.306e-07 | 30.131 | 27.546 | 2 |
| Reciprocal, full panel | 4/3 | 0.073 | 14.110 | 2.724e-05 | 33.293 | 26.559 | 1 |
| Reciprocal, full panel | 2 | 0.619 | 20.755 | 0.003 | 36.787 | 25.650 | 1 |
| Reciprocal, full panel | 5 | 8.091 | 32.982 | 0.987 | 41.467 | 24.090 | 1 |
| Log TMB, Merkel held out | 0.5 | 0.000e+00 | 0.310 | 0.000e+00 | 17.408 | 32.783 | 2 |
| Log TMB, Merkel held out | 1 | 0.000e+00 | 9.072 | 0.000e+00 | 26.169 | 28.980 | 2 |
| Log TMB, Merkel held out | 4/3 | 0.000e+00 | 12.708 | 0.000e+00 | 29.805 | 27.656 | 2 |
| Log TMB, Merkel held out | 2 | 0.000e+00 | 17.833 | 0.000e+00 | 34.930 | 26.109 | 1 |
| Log TMB, Merkel held out | 5 | 7.750 | 29.415 | 0.000e+00 | 46.512 | 24.299 | 1 |
| Log TMB, full panel | 0.5 | 0.000e+00 | 2.403 | 0.000e+00 | 19.348 | 31.876 | 2 |
| Log TMB, full panel | 1 | 0.000e+00 | 11.086 | 0.000e+00 | 28.030 | 28.281 | 2 |
| Log TMB, full panel | 4/3 | 0.000e+00 | 14.690 | 0.000e+00 | 31.634 | 27.059 | 2 |
| Log TMB, full panel | 2 | 0.000e+00 | 19.769 | 0.000e+00 | 36.713 | 25.669 | 1 |
| Log TMB, full panel | 5 | 9.777 | 31.247 | 0.000e+00 | 48.191 | 24.255 | 1 |

### Part-A-fixed viral extension scored on Part B

Primary calibration uses broad Part A MCPyV+ ORR. Only v is fitted; the scale is fixed by the scenario. Virus-negative predictions remain unchanged. Residuals are observed minus predicted, in pp.

| c | Frozen v | B+ predicted % | B- predicted % | B+ residual pp | B- residual pp | B+ interval distance pp | B- interval distance pp |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5 | 9.945465 | 28.094 | 17.008 | +6.206 | +31.592 | 0.000 | 14.892 |
| 1 | 9.700465 | 27.885 | 28.786 | +6.415 | +19.814 | 0.000 | 3.114 |
| 4/3 | 9.537132 | 27.743 | 32.833 | +6.557 | +15.767 | 0.000 | 0.000 |
| 2 | 9.210465 | 27.453 | 37.448 | +6.847 | +11.152 | 0.000 | 0.000 |
| 5 | 7.740465 | 26.047 | 43.852 | +8.253 | +4.748 | 0.000 | 0.000 |

MSE and interval MSE are in squared percentage points; RMSE/MAE are pp. Repair requires smaller MSE and no interval misses. These are descriptive aggregate scores, not patient-level prediction scores or significance tests.

| c | Baseline MSE | Extension MSE | Change MSE | Baseline RMSE | Extension RMSE | Extension MAE | Extension interval MSE | Misses | Repair rule |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0.5 | 1087.264 | 518.274 | -568.990 | 32.974 | 22.766 | 18.899 | 110.882 | 1 | fails |
| 1 | 784.549 | 216.882 | -567.667 | 28.010 | 14.727 | 13.115 | 4.850 | 1 | fails |
| 4/3 | 712.551 | 145.805 | -566.745 | 26.694 | 12.075 | 11.162 | 0.000 | 0 | passes |
| 2 | 650.417 | 85.620 | -564.797 | 25.503 | 9.253 | 8.999 | 0.000 | 0 | passes |
| 5 | 588.665 | 45.326 | -543.339 | 24.262 | 6.732 | 6.500 | 0.000 | 0 | passes |

Required calibration sensitivities (all scored against broad B, not a missing matched B outcome). The CI-endpoint columns are conditional Part A calibration scenarios for predicted B+; they are not predictive confidence intervals.

| Calibration | c | v at point ORR | B+ % | B RMSE pp | B misses | B+ % from A lower CI | B+ % from A upper CI |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| broad | 0.5 | 9.945465 | 28.094 | 22.766 | 1 | 15.510 | 43.486 |
| broad | 1 | 9.700465 | 27.885 | 14.727 | 1 | 15.009 | 43.473 |
| broad | 4/3 | 9.537132 | 27.743 | 12.075 | 0 | 14.668 | 43.463 |
| broad | 2 | 9.210465 | 27.453 | 9.253 | 0 | 13.970 | 43.445 |
| broad | 5 | 7.740465 | 26.047 | 6.732 | 0 | 10.564 | 43.360 |
| tmb_subset | 0.5 | 8.624119 | 25.849 | 23.124 | 1 | 9.587 | 48.400 |
| tmb_subset | 1 | 8.379119 | 25.594 | 15.304 | 1 | 8.965 | 48.400 |
| tmb_subset | 4/3 | 8.215785 | 25.420 | 12.796 | 0 | 8.545 | 48.400 |
| tmb_subset | 2 | 7.889119 | 25.065 | 10.238 | 0 | 7.698 | 48.400 |
| tmb_subset | 5 | 6.419119 | 23.333 | 8.450 | 0 | 3.890 | 48.399 |

### Canonical input and implementation fingerprints

SHA-256 uses UTF-8 text after CRLF/CR normalization to LF. Path keys are repository-relative `Path.as_posix()` values. README and the manifest itself are excluded to avoid a circular fingerprint; the test regenerates this entire block.

| Path | SHA-256 |
| --- | --- |
| `research/immune-escape-pancancer-test/analysis.py` | `06a7e15a24fa23b52f4c74b3659c678dd2778a88fdb429b0712e170e4b564b1c` |
| `research/immune-escape-pancancer-test/data/tmb_orr.csv` | `b34c1dd8d8d2c403f7251cf02611283d5326cbdc2caf2c72a6afadcf5ddc540b` |
| `research/merkel-viral-antigen-test/PROTOCOL.md` | `22923506eff606d81284a476a56813ec0773d73d238a5dcf633067344e711ea1` |
| `research/merkel-viral-antigen-test/analysis.py` | `39c6d7be8b81a02c8b6ee5c9ffacff37df30e9d705a6745039424af4cebcde31` |
| `research/merkel-viral-antigen-test/data/frozen_curves.csv` | `1acc4eff8810b2afd59c6e5cac250705f2e4e4f50b1c5f079a0c18308bb893f9` |
| `research/merkel-viral-antigen-test/data/part_a.csv` | `3594998062d9e59450066241db1613c59f00af8ec6f3773f365a0acde92b93fb` |
| `research/merkel-viral-antigen-test/data/part_a_calibration.csv` | `ef597c8dd3ff770848fb06092853e620136a014c1fdd76c9e3873618cbf3af6e` |
| `research/merkel-viral-antigen-test/data/part_b.csv` | `c3b197377b000c9a79816d75f1cdc09cfe5f8ac2e86ba007913eda6e44167eb3` |
| `tests/test_research_merkel_viral_antigen_test.py` | `b6d00ee08d32563ac3d5e1b97ea930f091aa28f906bf7681c18613709930fa9a` |

<!-- END RECOMPUTED RESULTS -->

## What this teaches about the functional-load model

The TMB-only observation bridge needs to distinguish mutation-derived antigen
supply from other presented targets. A useful conceptual decomposition is
`k_i = k_mutational,i + k_nonmutational,i`, with viral antigens among the latter.
[#59](../immune-escape-neoantigen-clonality/README.md) already defines `k_i` as
functional, equally weighted presented load. Its theorem's harmonic exposure
depends on all those loads; a mutation assay is an incomplete proxy whenever
non-mutational targets contribute. The viral offset tested here acts only as
a proxy for the missing contribution in that empirical mapping.

The fitted offset has mutation-scale units and cannot be read as a count of
viral peptides, their antigenic quality, presentation, coverage or immune
killing strength. Viral status also correlates with tumour biology and immune
state; treatment history differs across cohorts. The Part B report itself
describes viral-group immune-expression differences, and neither cohort
provides the calibrated functional load or response-to-eradication bridge
required to test #59 directly. Common antigen quality and loss dynamics would
also need justification before adding viral and mutational counts.

The transferred offset's point-score improvement is a conditional empirical
reason to measure that contribution. Its scale-dependent repair, the unchanged
negative-control failure, broad Part A interval sensitivity, missing matched
Part B response table and successful median-envelope relaxation prevent
claiming that this analysis identified it. #59's conditional theorem remains
untested by these marginals; #67's failed pan-cancer improvement criterion is
not reversed by this focused follow-up.

## Falsification attempts and retained counterevidence

| Attempt | Outcome and limit |
| --- | --- |
| Reproduce the actual motivating held-out curve | Exact replay from pinned input/code; the full-panel constants were not substituted for the held-out fold. |
| Fit both viral groups with a shared scale in each cohort and across cohorts | Reciprocal curves admit common scales only beyond the assumed bands under median substitution. Therefore neither “no scale exists” nor “ordinary conversion explains it” is supported. |
| Use the primary simpler log-TMB comparator | Its Part B status-scale sets are disjoint even over the full positive domain. This distinguishes the fixed shapes under the plug-in bridge; it is not evidence that the reciprocal curve is biologically correct. |
| Blame printed rounding | Expanded regions retain the same qualitative bounded-band failure and log-comparator inconsistency. Rounding tolerance is not a confidence interval for the median. |
| Change the original training-panel choice | Both the primary held-out curve and full-panel reciprocal sensitivity retain the large-scale qualification. No curve is selected by its new score. |
| Replace broad Part A ORRs with its published TMB-subset ORRs | Reported as a secondary sensitivity; it changes the compatible region and calibration but cannot supply the missing matched Part B outcome. |
| Remove median substitution while retaining median constraints | The conservative envelope becomes compatible throughout both assumed bands. This defeats an unconditional claim that these summaries require a viral offset. It does not fit or observe an actual distribution. |
| Transfer only one Part-A-fitted viral parameter | Point loss improves at every fixed scale; the complete repair rule fails at the identity reference because the virus-negative control is unchanged and incompatible. Passing scenarios are all displayed. |
| Stress Part A calibration with its published interval endpoints | Wide transferred positive-stratum predictions remain possible. These are descriptive calibration scenarios, without a predictive coverage guarantee. |
| Infer missing viral-by-response cells from separate biomarker margins | Rejected. Missing Part B subset-ORR rows remain blank and the matched test remains unavailable. |
| Treat overlapping reports or a sequencing subset as independent replications | Rejected. Cohort lineage, cutoffs and denominator differences remain explicit. |
| Make causal or clinical claims from a better score | Rejected. No randomized viral comparison, exact assay calibration, antigen measurement or clinical decision rule is supplied. |

## Honest novelty assessment

The broad biology is established prior work. [Goh et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4823115/)
described contrasting virus-associated and mutation-rich Merkel genomes and
their antigen implications. [Iyer et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3207011/)
demonstrated virus-specific T-cell recognition and oncoprotein presentation.
The two JAVELIN papers already report the viral/TMB/response patterns used
here. The [original #67 source audit](../immune-escape-pancancer-test/TMB_ORR_SOURCES.md)
credits prior ecological TMB-response regressions and alternative antigen
quality explanations.

The contribution here is limited to a repository-specific replay of a frozen
curve, complete common-scale compatibility regions, a temporally frozen
single-parameter cohort transfer, and an explicit counterexample to treating
median substitution as an identified mechanism. The median bounds use
elementary monotonicity; no new probability theorem is claimed. This is not a
new viral mechanism, novel biomarker, first ecological association, efficacy
estimate or demonstrated bibliographic priority. The targeted primary-source
check covered the named trial reports, their assay definitions and the cited
biological antecedents; it was not a comprehensive systematic review.

## Reproduce, verify and review

From the repository root:

```bash
python -B research/merkel-viral-antigen-test/analysis.py --check
python -B research/merkel-viral-antigen-test/analysis.py --json
python -B -m unittest discover -s tests -p 'test_research_merkel_viral_antigen_test.py' -v
```

The JSON report also supplies every baseline residual, interval distance, MAE
and interval MSE; the baseline Markdown table displays predictions, RMSE and
miss counts for compact comparison.

The test uses only Python's standard library and a bounded workload well
below the requested minute. It independently checks frozen fits, analytic
inversion and intersections, all scores, rounding and median sensitivities,
published-count consistency where counts are actually reported, denominator
separation, absence of Part B calibration leakage, canonical text hashing,
POSIX path keys, the frozen calibration and every generated README result.
Canonical fingerprints are also in [fingerprints.csv](data/fingerprints.csv).
An interval or hash check cannot verify publication truth or source permissions.

`--write-report` updates only generated results/fingerprints after a reviewed
code change. It never changes the frozen Part A artifact. `--freeze-part-a`
was used only before the calibration commit and refuses to overwrite an
existing artifact. Its inputs never include Part B.

Required repository gates:

```bash
python -B scripts/validate.py
python -B -m unittest discover -s tests -v
node --test tests/*.test.*
```

The existing PR workflow discovers the dedicated test; no workflow, dependency,
index, catalog, result bundle or pre-existing research file is changed. The
only additions are this directory and the single named test file. CI checks
software and bookkeeping, without scientific or clinical endorsement.
