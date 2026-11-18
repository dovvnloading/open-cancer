# Prediction confrontation with published aggregate evidence

**Nonclinical, unreviewed research. No diagnosis, treatment or surveillance advice.**

This directory confronts the mathematical predictions in merged PRs
[#55](https://github.com/dovvnloading/open-cancer/pull/55),
[#57](https://github.com/dovvnloading/open-cancer/pull/57),
[#58](https://github.com/dovvnloading/open-cancer/pull/58),
[#59](https://github.com/dovvnloading/open-cancer/pull/59) and
[#60](https://github.com/dovvnloading/open-cancer/pull/60) with published
aggregate observations. It does not change any of those notes.

**Outcome:** no core mathematical prediction is empirically falsified by an
eligible aggregate test in this audit. None is empirically validated either.
The precise escape, schedule, generational moment, threshold and holiday-law
predictions lack at least one necessary public aggregate or design control.
A weaker clonal-burden response association is consistent with the published
estimates. A blanket adverse HLA-LOH association has nominal counterevidence
but no rejection after the specified multiplicity correction. An exploratory
ecDNA calculation identifies a strong conditional challenge to universally
fair independent allocation; its sampling and test-calibration assumptions
are not established by the paper.

These are substantive limits on what the programme can currently claim.
The next useful work is to obtain the missing **aggregate summaries and
calibrations**, not to relabel a correlated clinical endpoint as mechanistic
escape or to fit a threshold from an uncalibrated drug concentration.

## Scope, timing and provenance

Read first: [research index](../README.md), [data policy](../../DATA_POLICY.md),
and all five merged notes. Base main commit:
`864bd6aa09cbd2bd100a141750a99deb2688921c`.
The [protocol](PROTOCOL.md) was committed locally as `29bedbd` before new literature
retrieval and extraction. It is a prospective analysis plan for a
retrospective targeted audit, **not a registered or blinded study**: the
merged notes already reveal some of the motivating observations.
[Amendments](AMENDMENTS.md) preserve post-extraction clarifications and the
explicitly exploratory segregation test. They do not overwrite the protocol.
The [publication-transport record](AMENDMENTS.md#publication-transport-record)
links the matching remote protocol tree and distinguishes local freeze time
from the later GitHub publication time.
Access/search cut-off: **2026-11-17, America/New_York**.

The five small [data CSVs](data/) contain only scalar published aggregates.
Each row has a source URL, DOI, exact text/table/figure locator, population,
unit, access date, reuse basis and interpretation notes. Fitted parameters
are labelled as fits. No patient, cell or culture records, individual figure
points, reconstructed survival records, source workbooks, source PDFs or
source images are included. Mixed source-data workbooks containing patient
records were not downloaded for this analysis. Small aggregate panels are retained
even when their direction is inconvenient or their eligibility fails.

This is a targeted falsification search, not an exhaustive systematic review.
The [search and source audit](SOURCES.md) includes the user-named studies,
null/adverse/contrary results, follow-up references, rights exclusions,
access limitations, overlap, and source inconsistencies. Absence from the
located aggregates is not a claim that no suitable aggregate exists anywhere.

## Sharp predictions and tests specified before extraction

Symbols below follow their source notes. A population size, hazard or dose
parameter in a theorem is not automatically the similarly named clinical
measurement. All assumptions are part of the testable conjunction.

### Immune escape and neoantigen clonality — PR59

Source: [immune-escape note](../immune-escape-neoantigen-clonality/README.md),
Theorem 6 and its observation-bridge section.

**I1: reciprocal functional-load escape.** For fixed rates, a suppressed
presenting source, no pre-existing presentation-null cells, and increasing
minimum functional load,

\[
 k p_{\rm escape}(k)\longrightarrow\sigma h/\kappa,
 \qquad Q\longrightarrow\exp[-\sigma hN/(\kappa k)].
\]

For heterogeneous founders the exponent depends on the reciprocal-load sum;
equal groups at half and one-and-a-half of the common load give the stated
factor of four-thirds. A test needs **functional presented load**, effective
founder count, independently measured rates and actual escape/extinction
counts. A response-based alternative needs a separately calibrated
false-positive/sensitivity bridge, justified finite-load error, and at least
three disjoint finite positive load bins with response numerators and
denominators. The planned test is the note's conservative binomial
acceptance-interval intersection for a common coefficient, including dose
and bridge uncertainty. An empty intersection rejects the calibrated
model-plus-bridge; a fitted coefficient with no independent bridge does not.

The retained Litchfield/Ravi estimates are response associations on
standardized or log-feature scales. They supply none of the requisite dose
bins, founder counts or bridge calibration. Total TMB, an open-ended
high-TMB bin, HLA-LOH prevalence and a survival hazard ratio cannot fill those
fields. **Verdict: untestable with current public aggregates.**

**I2: weaker clonal-burden association.** Compare each eligible published
clonal-response OR with one, retaining its source interval and all selected
contrasts. Litchfield and Ravi are consistent with a positive association.
This does not identify a reciprocal curve, its coefficient, a causal benefit,
or a matched-background effect. **I2-H**, the lower-heterogeneity component
at fixed burden, remains **untestable**: comparing separate clonal/subclonal
ORs without their joint design and covariance cannot estimate it. McGranahan
was screened, but affirmative redistribution rights for the inspected copy
were not established; its numbers are not copied here.

**I3: blanket adverse HLA-LOH observation bridge.** This intentionally tests
a tempting stronger extrapolation that PR59 itself cautions against, not
the theorem about complete presentation loss. Response OR below one and
death HR above one are adverse; effects are kept on their own endpoints.
Test each eligible reported contrast with the frozen joint Holm rule.
The nominal contrary Yang and Lozac'hmeur intervals are retained. None gives
a corrected rejection. Under the frozen rules a non-rejection is **not** evidence of consistency
with a universal adverse effect: the estimates have opposing directions,
are largely unadjusted or otherwise differently adjusted, and do not
calibrate LOH as complete loss of presentation. **Verdict: untestable with
current public aggregates (conflicting directions; no corrected rejection).**
Nominal counterevidence is retained; non-significance is not equivalence.

### ctDNA lead time and sampling — PR58

Source: [surveillance note](../ctdna-surveillance-schedule/README.md),
Theorem 6, its sharp support bounds and geometric-schedule result.

**C1: alternate-grid thinning.** On an exact interval-\(d\) grid, persistent
detectability and unchanged recurrence dates, average both retained parities
and give missed detections zero lead:

\[
 E[L_{\rm thin}\mid L]=\tfrac12\{L+(L-d)_+\}.
\]

For observed original mean \(m\) and support \([a,b]\), define

\[
 T_{\min}=\frac{b-m}{b-a}\min(a,d)+\frac{m-a}{b-a}\min(b,d),
 \qquad T_{\max}=\min(m,d).
\]

Then the thinned mean lies in
\([m-T_{\max}/2,m-T_{\min}/2]\). The prespecified test compares an
**observed same-cohort both-parity, zero-imputed thinned mean** with the
outward rounding envelope. Reinert permits a conditional forecast only:
planned quarterly sampling does not document an exact realized grid, and
the reported mean is conditional on detection. Henriksen's restriction to
blood draws at CT dates is the closest actual retrospective comparison, but
uses a different selection and reports a median, not the target mean.
**Verdict: untestable.** The computed envelope is not an observed success of
the theorem and is not a population confidence interval.

**C2: geometric schedule.** With fixed residual burden, a perfect detection
threshold, unimodal log-growth distribution and the note's horizon
conditions, the optimal ratio is
\(q=\log(C/B_*)/\log(D/B_*)\). A test needs independently calibrated
threshold/residual quantities and aggregate capture counts for geometric
and equal-calendar schedules with the same test budget and horizon; compare
the capture risk difference with its interval. Data-defined slow/fast growth
groups cannot independently test log-growth unimodality. Gains obtained by
extrapolating a growth model cannot independently validate that model's
schedule prediction. **Verdict: untestable.** No adaptive policy is tested.

**C3: restricted observation-model stress test.** A common detection window,
exact periodic perfect testing and continuous recurrence observation imply
lead-range width at most \(d\). A known imaging-lag bound \(I\) changes the
bound to \(d+I\). The retained range is wider than the blood-test interval,
but the actual imaging endpoint and common-window conditions do not satisfy
the restricted conjunction. **Verdict: untestable**, not a falsification of
PR58's general heterogeneous model. A planned imaging interval is not a
verified bound on every realized diagnostic delay.

### ecDNA segregation, killing and clustering — PR55 and PR60

Sources: [segregation note](../ecdna-segregation-therapy/README.md) and
[killing/clustering note](../ecdna-killing-laws-and-clustering/README.md).

**E1: generation-two Fano under exponential killing.** Monomorphic founders,
all generation-two cells including zeros, and pooling by expected cell
number give
\(F_2=(3+s)/(2(1+s))\), strictly between one and three-halves for interior
\(s\), independent of founder copy count and \(u\). Test a measured Fano
and its replicate interval against the prediction with independently
calibrated \(s\). Cross-sectional calendar-time cells, positive-only counts
and sister-fraction fits are different estimands. **Verdict: untestable.**

**E2: killing thresholds.** PR55 predicts
\(2u(1-\sqrt{1-s})=1\), giving \(s_*=3/4\) when \(u=1\).
PR60's linear-hazard fair-binomial model brackets
\(c_*=\gamma/b\) between \(0.378474994554\) and \(0.378474994555\).
A test needs independent copy-specific death/division calibration, founder
state, segregation kernel and aggregate surviving-well counts at known
finite genealogical generations. Compare those counts to the model's
finite-generation probabilities with binomial intervals. Finite persistence
cannot by itself falsify eventual extinction. Drug concentration is not
\(\gamma/b\). **Verdict: untestable.** These numbers are theoretical constants
from the merged notes, not newly measured estimates.

**E3: conditional segregation moments.** Given maternal pre-replication
copy count \(k\), fair independent segregation predicts daughter mean
\(k\), variance \(k/2\), and Fano \(1/2\). The symmetric clustered/bias
family predicts

\[
 F(k)=(1-R)/2+Rk,\qquad
 R=(2p-1)^2+4p(1-p)\rho.
\]

The required observable is conditional daughter moments by maternal count
with uncertainty, or a valid published test/parameter interval. In particular,
\(R\ne\rho\) when allocation is biased and \(\rho<1\). Shibata's fitted parameters have no
intervals and the linear relationship is imposed by their fitted family;
transforming them does not independently test PR60. Its plotted KS tests
have sister dependence and calibration concerns; two HSR controls and an
internally disputed clone must not be treated as unambiguous ecDNA. Lange's
main/supplement discrepancy is retained. **Verdict: untestable for both the
precise moment constraint and an unconditional rejection of its baseline.**

**E4: generation-two dose curves.** For **one-copy founders**, all cells
including zeros, fixed \(R\), and calibrated linear dose \(c\), set
\(A=2+(3-R)c\), \(D=2+4c+(1+R)c^2\), \(B=3+R+4c\).
The mean is \(A/D\) and Fano is
\((1+R)B/A+(1-R)/2-A/D\): mean falls and Fano rises with dose.
Compare measured moments and their uncertainty with these curves. Equal
second-moment curves cannot identify the killing threshold; the note's
matched beta/bias example requires a fourth-moment distinction. No located
aggregate has the joint founder, generation, dose and uncertainty inputs.
**Verdict: untestable.**

**E3-X: exploratory baseline challenge, specified after seeing the source
percentages.** Lundberg reports the fraction of divisions whose allocation
departs from a fair-binomial test, with known denominators and test level.
Retain every untreated cell line, enumerate all compatible integer rejected
counts under outward whole-percent rounding, compute a conservative exact
binomial upper tail and Wilson interval, and adjust all panels and the
combined family. These are rejection-frequency effects, not measured
allocation bias or Fano inflation. The independent-binomial calculation
rejects the fair-allocation-plus-observation conjunction for one cell line.
However, the paper does not document the individual test's tail construction
or sampled divisions' allocation across independent cultures. A
Markov bound removes independence but still assumes calibrated marginal
tests; it cannot reject at the chosen level and does not prove compatibility.
**Empirical verdict: untestable; conditional inconsistency only.** This
calculation cannot falsify or validate the clustered extension or dose law.

### Induced versus selected resistance — PR57

Source: [resistance note](../induced-vs-selected-resistance/README.md),
Theorem 3, assumptions A3/A4 and stated failure strata.

**R1: holiday moments and a held-out prediction.** With independent cellular
retention, a shared positive clone-size law with finite second moment and
nonzero probability of nonsingleton clones, unchanged exposure/readout and
no holiday growth, death or acquisition,

\[
 \mu_j=B+Cr^j,\quad C>0,\;0<r<1,\qquad
 W_j=\operatorname{Var}(Y_j)-\mu_j=A_0+A_1r^j+A_2r^{2j},
\]

where \(A_1\ge0,A_2>0\). Independent cultures at holidays
\(0,H,2H,3H\) should permit

\[
 r=\frac{\mu_1-\mu_2}{\mu_0-\mu_1},\qquad
 \widehat\mu_3=\mu_2-\frac{(\mu_1-\mu_2)^2}{\mu_0-\mu_1}.
\]

The planned test propagates a simultaneous uncertainty box through this
held-out prediction and the three-arm inverse admissibility constraints;
ordinary marginal intervals are not a simultaneous box. Reject only if no
compatible inputs satisfy the constraints. The located fluctuation,
resensitization and lineage studies do not publish this matched quartet.
**Verdict: untestable.** A singleton-clone failure stratum or different
clone-size laws is already outside the stated theorem.

**R2: independent-retention assumption.** Conditioned on initial state and
founder count, residual within-lineage retention covariance should be zero.
A test needs conditioned paired aggregate counts or covariance and uncertainty.
Iyer's approximate lineage dispersion is a useful contextual challenge but
is not the required conditional statistic. Even independent binomial
retention with variable founder number \(N\) gives

\[
 \operatorname{Var}(Y)/E[Y]=1-p+p\operatorname{Var}(N)/E[N].
\]

The excess over a fixed-founder comparator therefore does not identify
correlated switching. **Verdict: untestable.** No p-value or confidence
interval is invented from an approximate VMR and a lineage count.

## Numerical methods and error control

[analysis.py](analysis.py) reads the CSVs offline and generates the entire
numerical block below. The only dependencies are Python's standard library.
The formal family contains all eligible clonal/HLA contrasts specified in
[AMENDMENTS.md](AMENDMENTS.md). Overlapping studies and HLA loci are not
independent replications; Holm remains valid under arbitrary dependence
when its component p-values are valid. No meta-analysis is constructed.
Total and subclonal TMB are context, not extra selectively interpreted tests.

Fisher's two-sided p uses exact integer hypergeometric probability ordering.
The central exact conditional OR interval inverts noncentral hypergeometric
tails at fixed margins; it need not invert that two-sided p definition.
The point estimate is the crude cross-product OR, not a conditional MLE.
Group risks use Wilson intervals, risk difference uses Newcombe's
Wilson-based construction, and risk ratio uses a log-Wald interval without
pseudocounts. All are marginal intervals; formal rejection uses Holm, not
an isolated displayed interval excluding the null.

Source regression/meta-analysis p-values are transcribed, not regenerated
from unavailable records. Ravi's raw p is approximated from the rounded
log-OR interval. A separate sensitivity uses the upper rounding endpoint of
the source BH-adjusted p, conditional on ordinary BH adjustment and rounding
to the nearest displayed last digit. This changes no verdict. The
CI-derived family is not advertised as a numerically exact error guarantee.
No subgroup effect is converted to a causal or calibrated mechanistic effect.

The ctDNA support calculation is algebraic. The outward envelope allows
nearest-tenth rounding of the mean and endpoints, holding the theorem's
exact cadence assumption fixed. On this rectangle the upper bound increases
with mean and upper support and decreases with lower support, so corner
evaluation gives the outer envelope. It covers decimal rounding only,
not sampling variation, unknown imaging delay, missed visits or assay error.

The exploratory segregation tail uses rational arithmetic. Reported
percentages are converted only using explicitly reported denominators;
unknown denominators are never inferred. The combined immune/segregation
Holm calculation is an additional conditional sensitivity, not a repair
for the segregation study's undocumented calibration and dependence.

<!-- BEGIN RECOMPUTED RESULTS -->
### Recomputed results

Sampling intervals below are marginal 95% intervals; model/forecast bounds, ranges, IQRs and rounding envelopes are not confidence intervals. Holm uses all 14 specified immune contrasts; they are not 14 independent studies. Source regression/KS test inputs are absent from the retained aggregates, so their p-values are transcribed; all subsequent transformations are checked.

| CSV | Published scalar rows |
| --- | ---: |
| `immune.csv` | 30 |
| `hla_response.csv` | 36 |
| `ctdna.csv` | 34 |
| `ecdna.csv` | 87 |
| `resistance.csv` | 3 |

#### Immune associations

OR > 1 favors response; HR > 1 means greater death hazard. These are different endpoints. Litchfield ORs are per within-study SD, including the binary HLA feature; Ravi uses a log-transformed feature with no calibrated doubling interpretation. The q-bound column uses an upper bound on the source BH-adjusted value conditional on nearest-last-digit rounding; it is not a claimed raw p.

| Contrast | Effect | Source CI | p used | p provenance | Holm 14 | Holm, Ravi q bound |
| --- | ---: | --- | ---: | --- | ---: | ---: |
| Clonal TMB; Litchfield | OR 1.74 | [1.41, 2.15] | 2.9e-07 | source reported | 4.06e-06 | 4.06e-06 |
| Total TMB; Litchfield | OR 1.7 | [1.33, 2.17] | 1.9e-05 | source reported | excluded | excluded |
| Subclonal TMB; Litchfield | OR 1.18 | [0.99, 1.41] | 0.07 | source reported | excluded | excluded |
| HLA LOH; Litchfield | OR 1.14 | [0.95, 1.36] | 0.16 | source reported | 1 | 1 |
| Clonal neoantigens; Ravi | OR 5.4 | [2.7, 11] | 2.52351e-06 | approximate from rounded CI | 3.28057e-05 | 0.00325 |
| HLA-A LOH; Lozac’hmeur LUAD | HR 0.6 | [0.37, 0.96] | 0.032 | source reported | 0.352 | 0.352 |
| HLA-A LOH; Lozac’hmeur LUSC | HR 1.64 | [0.79, 3.41] | 0.183 | source reported | 1 | 1 |

Ravi source BH q = 0.0002; rounding-conditional upper bound = 0.00025; reconstructed log-OR SE = 0.358334, z = 4.70622. Total/subclonal TMB are contextual comparators, excluded from formal inference; no difference between their correlated ORs is tested.

| Yang cohort / HLA locus | LOH R, NR | Intact R, NR | Crude response OR | Exact conditional OR CI | Fisher two-sided p | Holm 14 |
| --- | --- | --- | ---: | --- | ---: | ---: |
| Urothelial / anti-PDL1 / A | 17, 26 | 29, 107 | 2.41247 | [1.06981, 5.33246] | 0.0265574 | 0.318689 |
| Urothelial / anti-PDL1 / B | 17, 39 | 30, 113 | 1.64188 | [0.75998, 3.46304] | 0.194092 | 1 |
| Urothelial / anti-PDL1 / C | 12, 30 | 35, 109 | 1.24571 | [0.522269, 2.8334] | 0.552851 | 1 |
| Melanoma / anti-PD1 / A | 8, 3 | 5, 5 | 2.66667 | [0.319033, 24.4908] | 0.386997 | 1 |
| Melanoma / anti-PD1 / B | 4, 3 | 7, 7 | 1.33333 | [0.153345, 12.5985] | 1 | 1 |
| Melanoma / anti-PD1 / C | 6, 5 | 6, 4 | 0.8 | [0.101087, 6.14873] | 1 | 1 |
| Melanoma / anti-CTLA4 / A | 3, 17 | 4, 13 | 0.573529 | [0.0721764, 4.11672] | 0.679995 | 1 |
| Melanoma / anti-CTLA4 / B | 3, 10 | 4, 22 | 1.65 | [0.20036, 11.7492] | 0.666416 | 1 |
| Melanoma / anti-CTLA4 / C | 4, 12 | 3, 16 | 1.77778 | [0.244215, 14.278] | 0.677217 | 1 |

| Yang cohort / locus | LOH response risk (Wilson CI) | Intact response risk (Wilson CI) | Risk difference (Newcombe CI) | Risk ratio (log-Wald CI) |
| --- | --- | --- | --- | --- |
| Urothelial / anti-PDL1 / A | 0.395349 [0.263653, 0.544209] | 0.213235 [0.152771, 0.289455] | 0.182114 [0.029952, 0.342785] | 1.85405 [1.13497, 3.02871] |
| Urothelial / anti-PDL1 / B | 0.303571 [0.199003, 0.433359] | 0.20979 [0.151091, 0.283673] | 0.0937812 [-0.0342548, 0.236225] | 1.44702 [0.87025, 2.40607] |
| Urothelial / anti-PDL1 / C | 0.285714 [0.17167, 0.435672] | 0.243056 [0.180269, 0.319194] | 0.0426587 [-0.0944659, 0.20523] | 1.17551 [0.672581, 2.05451] |
| Melanoma / anti-PD1 / A | 0.727273 [0.434355, 0.902539] | 0.5 [0.236593, 0.763407] | 0.227273 [-0.166661, 0.543661] | 1.45455 [0.709628, 2.98143] |
| Melanoma / anti-PD1 / B | 0.571429 [0.250458, 0.84178] | 0.5 [0.267992, 0.732008] | 0.0714286 [-0.324614, 0.427684] | 1.14286 [0.499222, 2.61632] |
| Melanoma / anti-PD1 / C | 0.545455 [0.280092, 0.787287] | 0.6 [0.312674, 0.83182] | -0.0545455 [-0.406906, 0.321007] | 0.909091 [0.433882, 1.90477] |
| Melanoma / anti-CTLA4 / A | 0.15 [0.0523687, 0.360419] | 0.235294 [0.0955504, 0.472618] | -0.0852941 [-0.341915, 0.167301] | 0.6375 [0.165247, 2.45939] |
| Melanoma / anti-CTLA4 / B | 0.230769 [0.0817953, 0.502564] | 0.153846 [0.0615003, 0.335312] | 0.0769231 [-0.15786, 0.363977] | 1.5 [0.392475, 5.73285] |
| Melanoma / anti-CTLA4 / C | 0.25 [0.101821, 0.494983] | 0.157895 [0.0552047, 0.375655] | 0.0921053 [-0.171289, 0.35774] | 1.58333 [0.414129, 6.05354] |

#### ctDNA: forecasts and endpoint audit

| Quantity | Recomputed value | Interpretation |
| --- | --- | --- |
| Reinert original detected-only mean / range / planned cadence | 8.7 / [0.8, 16.5] / 3 months | Range is not a CI; exact grid not documented |
| Both-parity thinning mean forecast | [7.2, 7.7465] months | Conditional forecast, no observed target aggregate |
| Forecast loss of lead | [0.953503, 1.5] months | Sharp plug-in bounds, not a CI |
| Outward decimal-rounding forecast envelope | [7.15, 7.80538] months | Assumes nearest tenth, not sampling uncertainty |
| Original lead range width | 15.7 months; rounding [15.6, 15.8] | Cannot reject C3 without its observation assumptions |
| Reinert detected / recurrent | 14 / 16 | Descriptive; no cross-cohort rate comparison |
| Henriksen 2022 CT-date restriction | 6 / 18 earlier; among earlier detections median 10, IQR [8, 14] months | Not alternate-grid, both-parity, zero-imputed mean |
| Log of Henriksen 2022 group mean multiplier | slow 0.223144, fast 0.887891 per month | Not mean log-growth; data-defined groups, no bimodality test |
| Henriksen 2024 frequency extrapolation | 9 / 30 delayed; delayed subgroup median gain 2, IQR [1, 2] months | Model-derived, not an independently observed monthly schedule |
| TRACERx median / range / recurrent n | 119 / [0, 1137] days / 63 | Zeros included for absent/post-relapse detection; median not mean |

#### ecDNA: fit transformations and excluded sensitivity tests

Shibata R values below are transformations of fitted parameters with no uncertainty; n is printed in the source, with fractions plotted per daughter, not asserted independent mothers. KS adjustments are sensitivity arithmetic only. Including the disputed clone cannot resolve its classification or test calibration.

| Clone | Source p allocation | Source rho | Derived R | Source KS p | n | Scope | Holm 9 | Holm 10, ambiguous 7 included |
| --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| 1 | 0.468947 | 0.0884258 | 0.0919418 | 0.32 | 486 | ecDNA | 1 | 1 |
| 2 | 0.458063 | 0.112367 | 0.118612 | 0.96 | 108 | ecDNA | 1 | 1 |
| 3 | 0.463154 | 0.098527 | 0.103423 | 0.024 | 554 | ecDNA | 0.216 | 0.216 |
| 4 | 0.464332 | 0.097302 | 0.101896 | 0.049 | 484 | HSR control | excluded | excluded |
| 5 | 0.446876 | 0.138734 | 0.148457 | 0.46 | 494 | ecDNA | 1 | 1 |
| 6 | 0.452968 | 0.124863 | 0.132606 | 8e-06 | 452 | HSR control | excluded | excluded |
| 7 | 0.471172 | 0.0865586 | 0.089595 | 3.9e-10 | 558 | classification conflict | excluded | 3.9e-09 |
| 8 | 0.466738 | 0.0910859 | 0.0951082 | 0.81 | 508 | ecDNA | 1 | 1 |
| 9 | 0.473924 | 0.102812 | 0.105252 | 0.98 | 480 | ecDNA | 1 | 1 |
| 10 | 0.468222 | 0.0860257 | 0.0897175 | 0.045 | 608 | ecDNA | 0.36 | 0.36 |
| 11 | 0.446612 | 0.137344 | 0.14718 | 0.99 | 268 | ecDNA | 1 | 1 |
| 12 | 0.463886 | 0.0999508 | 0.104646 | 0.44 | 516 | ecDNA | 1 | 1 |

The eligible-by-classification fitted R range is [0.0897175, 0.148457]; this is not a confidence interval.

| Shibata arm ID | p allocation | rho | Derived R |
| --- | ---: | ---: | ---: |
| ctrl1_dmso | 0.428609 | 0.17924 | 0.195972 |
| ctrl1_jq1 | 0.313859 | 0.0877144 | 0.214151 |
| ctrlen2_dmso | 0.449921 | 0.130785 | 0.139505 |
| ctrlen2_jq1 | 0.438837 | 0.164261 | 0.176767 |

Lange supplementary Table SI1 is retained only for the source-version audit; conflicting main/SI representations and missing conditional moments prevent formal inference.

| SI1 panel | Source KS D | Source KS p | Source n |
| --- | ---: | ---: | ---: |
| pc3_myc | 0.065 | 0.375 | 200 |
| snu16_myc | 0.039 | 0.918 | 194 |
| snu16_fgfr2 | 0.063 | 0.415 | 196 |
| gbm39_egfr | 0.072 | 0.221 | 210 |
| colo_myc | 0.033 | 0.973 | 206 |
| tr14_mycn | 0.075 | 0.315 | 162 |
| tr14_cdk4 | 0.086 | 0.239 | 142 |
| colo320_dm | 0.049 | 0.996 | 66 |

#### E3-X: exploratory after extraction, conditional on valid tests and independent mitoses

The percentage-to-count conversion uses the published denominator and every compatible integer numerator. The Wilson intervals and binomial tail require independent sampling. Markov removes independence but still requires marginal test size at most the source cutoff. Excess/fold refer to the rejection-rate ceiling, not a biological allocation effect.

| Cell line | Published % / n | Compatible k | Rate, Wilson CI | Excess / fold vs nominal | Conditional p | Holm 3 | Conditional Holm combined 17 | Markov upper bound |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| chp212 | 2% / 100 | 2 | 0.02 [0.00550197, 0.0700118] | 0.01 / 2 | 0.264238 | 0.528476 | 1 | 0.5 |
| la_n5 | 10% / 100 | 10 | 0.1 [0.0552291, 0.174366] | 0.09 / 10 | 7.63159e-08 | 2.28948e-07 | 1.29737e-06 | 0.1 |
| mcix | 2% / 50 | 1 | 0.02 [0.00353926, 0.104954] | 0.01 / 2 | 0.394994 | 0.528476 | 1 | 0.5 |

#### Resistance: contextual dispersion, no holiday-law test

Iyer reports approximately p = 0.36, VMR = 2.4, with 65 retained lineages. The fixed-founder independent-binomial comparator is 1-p = 0.64; excess = 1.76, ratio = 3.75. No CI or p-value is calculable for the required conditioned contrast from these aggregates. Variable founder counts alone change this comparator.

#### Verdicts under the frozen rules

| ID | Verdict |
| --- | --- |
| I1 | **untestable with current public aggregates** |
| I2 | **consistent (clonal-burden association only)** |
| I2-H | **untestable with current public aggregates** |
| I3 | **untestable with current public aggregates (conflicting directions; no corrected rejection)** |
| C1 | **untestable with current public aggregates** |
| C2 | **untestable with current public aggregates** |
| C3 | **untestable with current public aggregates** |
| E1 | **untestable with current public aggregates** |
| E2 | **untestable with current public aggregates** |
| E3 | **untestable with current public aggregates** |
| E4 | **untestable with current public aggregates** |
| R1 | **untestable with current public aggregates** |
| R2 | **untestable with current public aggregates** |
| E3-X | **untestable with current public aggregates** |

I3 remains unresolved: nominal opposite HLA-LOH contrasts do not survive Holm correction, and the public aggregates cannot establish the universal adverse bridge. E3-X gives a conditional challenge, not an unconditional empirical falsification.
<!-- END RECOMPUTED RESULTS -->

## What the contrary observations would change

No unconditional corrected rejection of a core theorem was obtained.
Consequently the following are **minimal candidate repairs for challenged
applications**, not models fitted or validated by these aggregates.

| Challenged application | Minimal change | New aggregate needed to distinguish it |
| --- | --- | --- |
| Every HLA LOH is complete presentation loss with the same adverse response effect | Keep allele-specific antigen coverage; distinguish partial loss from the theorem's complete-loss state; permit a histology/regimen-specific response bridge | Matched aggregate response or survival contrasts by functional lost coverage, burden and regimen, with adjustment and uncertainty |
| Fair independent allocation is universal across ecDNA copy numbers and states | Permit correlated allocation or bias through a nonzero, possibly copy/state-dependent R; include counting/truncation error | Conditional second and fourth moments by maternal count, independent-mother counts and replicate uncertainty; measured counting-error summary |
| A blood-test interval alone bounds observed clinical lead time | Add the imaging observation lag to the threshold-crossing process; allow heterogeneous growth/windows | Same-cohort both-parity zero-imputed mean and actual imaging-lag summary; unchanged recurrence definition |
| Excess lineage dispersion uniquely proves dependent switching | Condition on founder number and initial state; if residual covariance remains, add one shared clone-level latent state | State/founder-stratified paired counts or residual covariance with culture-cluster uncertainty |
| Response increases with clonal burden, therefore escape scales reciprocally | Retain a separately calibrated observation bridge and functional burden; test finite-load error before transporting the limit | Disjoint finite dose-bin counts, founder/rate calibration and bridge bounds |

These changes are small mechanistically but may destroy an exact inverse or
threshold. They require a new derivation and a new aggregate test. For example,
a saturating immune killing law is already a stated failure mode in PR59;
this audit provides no evidence for choosing it over the original law.
Likewise, positive fitted R does not establish a fixed-R dose curve.

## Data provenance and reuse

| File | Primary sources and locators | Reuse basis retained in every row |
| --- | --- | --- |
| [immune.csv](data/immune.csv) | [Litchfield 2021](https://doi.org/10.1016/j.cell.2021.01.002), Results/Fig. 2A and STAR Methods; [Ravi 2023](https://doi.org/10.1038/s41588-023-01355-5), Fig. 1f-associated Results/Methods; [Lozac'hmeur 2024](https://doi.org/10.1038/s41698-024-00665-z), Fig. 3 caption and survival Results | Litchfield/Ravi CC BY; Lozac'hmeur CC BY-NC-ND; unchanged scalar facts, source rights preserved |
| [hla_response.csv](data/hla_response.csv) | [Yang, Kim & Kim 2022](https://doi.org/10.1007/s13258-021-01207-8), printed response/nonresponse labels in every main Fig. 1c–e panel | CC BY; counts transcribed from aggregate labels; no underlying records |
| [ctdna.csv](data/ctdna.csv) | [Reinert 2019](https://doi.org/10.1001/jamaoncol.2019.0528), corrected Results/Methods; [Henriksen 2022](https://doi.org/10.1158/1078-0432.CCR-21-2404), lead-time/growth Results; [Henriksen 2024](https://doi.org/10.1016/j.annonc.2023.11.009), sampling-frequency Results; [Abbosh/TRACERx 2023](https://doi.org/10.1038/s41586-023-05776-4), lead-time Results and Methods | CC BY except Henriksen 2022 CC BY-NC-ND; unchanged scalar facts with original attribution |
| [ecdna.csv](data/ecdna.csv) | [Lange 2022](https://doi.org/10.1038/s41588-022-01177-x), Table SI1; [Shibata 2026](https://doi.org/10.1038/s41467-026-76566-5), Fig. 3d and aggregate Supplementary Data 5/8; [Lundberg 2008](https://doi.org/10.1371/journal.pone.0003099), Statistical distribution of DMs at mitosis | Lange CC BY; Shibata CC BY-NC-ND; Lundberg Creative Commons Attribution as stated by the article; no invented license version |
| [resistance.csv](data/resistance.csv) | [Iyer 2025](https://doi.org/10.1371/journal.pcbi.1013446), Results around Fig. 6e and Table 1 caption | CC BY; approximate experimental summaries, excluding fitted simulation outputs |

The repository's MIT license does not relicense these source facts or remove
source restrictions. This noncommercial research audit transcribes unchanged
small numerical excerpts with attribution; it does not redistribute or adapt
source articles, expressive figures/tables or underlying databases. The
[Creative Commons FAQ](https://creativecommons.org/faq/) distinguishes
unchanged excerpts/format changes from adaptation. Sources without an
established affirmative reuse basis remain bibliographic links only.

## Hostile review, limitations and rerun

Transcriptions were checked at primary text, figure-label or aggregate
supplement-cell locators. Review preserved all Yang panels and both LOH
histologies, the Lange source discrepancy, all Shibata clone panels with
scope flags, and every untreated Lundberg panel. It caught the HSR/classification
trap, per-SD effect-scale trap, response-versus-extinction bridge, dependence
and marginal-versus-simultaneous interval problems. The detailed audit is in
[SOURCES.md](SOURCES.md). It records what was inspected and what was not;
source truth is not proved by a test suite or agreement among AI agents.

Remaining limitations include targeted-search selection and publication
bias, unresolved rights exclusions, overlapping source cohorts, confounding,
rounded source estimates, source-dependent assays/denominators, absent joint
uncertainty, and empirical tests whose calibration cannot be reconstructed
without disallowed records. A missing-data eligibility verdict is an explicit
abstention, not a negative experimental result. The code checks that its
manual eligibility annotations have evidence anchors; it cannot prove the
literature search complete or verify a primary-source transcription by itself.

Run from the repository root:

```bash
python -B research/prediction-confrontation/analysis.py --check
python -B research/prediction-confrontation/analysis.py --json
python -B -m unittest discover -s tests -p test_research_prediction_confrontation.py -v
```

Without flags, the script prints Markdown result tables. `--check` rejects a
stale README block. It does not fetch sources or modify files. The sole new
[test file](../../tests/test_research_prediction_confrontation.py) independently
recomputes all reported derived statistics from the CSVs, verifies interval
tail equations, source-field completeness, effect orientation, multiplicity,
source exclusions, eligibility guards and exact README synchronization. It
uses only the standard library and completes well below the requested
minute. The existing repository Python/Node/validation and PR CI gates are
run without changing workflows, requirements, indexes or any existing file.

To make the untestable predictions testable, publish the listed **aggregate
sufficient summaries**: calibrated disjoint bins for I1; a same-cohort
zero-imputed parity-thinning mean for C1; comparable schedule capture counts
for C2; independent generational moments/rate measurements for E1–E4; and
matched holiday count moments with simultaneous uncertainty and growth/readout
controls for R1/R2. For count-variance uncertainty, independent-culture counts
and aggregate power sums through the fourth moment can be sufficient;
individual records are not required. No recommendation here concerns care.
