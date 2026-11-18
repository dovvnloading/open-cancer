# ecDNA segregation clustering: aggregate calibration and its limits

**Nonclinical computational research. Unreviewed. No treatment advice.**
This directory contains published cell-line aggregates, original analysis,
mathematical proofs and reproducible checks. It contains no patient-level,
controlled-access, individual-cell, or re-identifiable records. A mathematical
certificate establishes a model statement; it does not establish biological
validity, therapeutic efficacy, or scientific priority.

## Finding and scientific verdict

Published beta-binomial segregation fits imply appreciably lower **model**
thresholds for eradication of an ecDNA-positive lineage than the independent
binomial model. The numerical shift, source-parameter rounding sensitivity,
and conservative propagation of segregation uncertainty are computed below.
The strongest falsifiable output is a fixed-copy sister-pair experiment with
an exact test against unbiased binomial segregation.

**The stronger requested empirical discovery is not established by the
permitted aggregates.** The usable source fits lack a sampling covariance,
profile likelihood or bootstrap distribution. The published observation
likelihood and sampling units also have consequential ambiguities. There is
therefore no defensible new empirical beta-binomial refit or calibrated
confidence distribution for true segregation dispersion in this submission.
We reuse the published maximum-likelihood fits, compute their implications,
and expose this limitation instead of manufacturing observations, confidence
intervals, or a successful discovery claim. The conditional sensitivity bands
below must not be quoted as empirical confidence intervals.

The starting theory is [PR55's segregation note](../ecdna-segregation-therapy/README.md),
[PR60's killing-law and clustering note](../ecdna-killing-laws-and-clustering/README.md),
and [PR63's evidence confrontation](../prediction-confrontation/README.md).
The present work neither modifies those notes nor treats their prior
transcriptions as authoritative primary sources. Every retained numerical
source extract was checked against the current primary text, rendered figure,
or supplement cell.

## Scope, provenance and search

### What was searched

The search date is **2026-11-17, America/New_York**. This is a targeted primary
literature and citation audit, not a systematic-review registration or a
guarantee that every published measurement has been located. Searches combined
`ecDNA`, `extrachromosomal DNA`, `double minute`, `segregation`, `partition`,
`mitosis`, `daughter`, `sister`, `binomial`, `beta-binomial`, `overdispersion`,
`clustering`, `retention`, `hitchhiking`, `live imaging`, and named investigators.
Backward references from Lundberg, Lange, Kanda, Hung, Nichols and Sankar and
forward/current searches were followed. Historical and contemporary sources
were searched even when their findings opposed clustering. Exact searches
also included `Scanu 2025 ecDNA segregation quantitative`, `double minute
segregation variance`, and `ecDNA segregation beta-binomial threshold linear
killing clustering`. Novelty searches included earlier amplification and
infinite-type branching models.

[data/literature_audit.csv](data/literature_audit.csv) records each inspected
publication, version/status, source URL, inspection depth, access limits,
reuse basis and inclusion/exclusion reason. Superseded preprints are labelled
duplicates, not independent studies. Some historical papers were accessible
only as abstracts: their undiscovered tables are **missing evidence**, not
proof that no quantitative data exist. Search indexing and access limitations
prevent an exhaustive completeness claim.

### Main source decisions

| Source | Measurement and role | Eligibility for numerical dispersion inference |
| --- | --- | --- |
| [Lundberg et al., PLoS ONE 2008](https://doi.org/10.1371/journal.pone.0003099) | Opposing anaphase/telophase pole counts; fractions of divisions rejecting fair binomial allocation; sister correlations; copy-dependent clustering and irradiation summaries | Rejection-frequency test conditional on calibration/independence. Rejection fractions and marginal sister correlations cannot identify rho without conditional totals and moments. |
| [Lange et al., Nature Genetics 2022](https://doi.org/10.1038/s41588-022-01177-x) | Sister-fraction distributions and all rows of Supplementary Table SI 1; live imaging | Retain every published KS summary as opposing evidence. Fractions, unspecified sample-unit dependence and simulation design do not supply a new conditional likelihood or dispersion interval. |
| [Shibata et al., Nature Communications 2026](https://doi.org/10.1038/s41467-026-76566-5) | Published folded beta-binomial and negative-binomial MLEs; aggregate model comparisons; figure sample sizes | Source-fit calibration only. No source sampling intervals. Keep baseline clones, chromosomal controls and perturbation arms distinct. |
| [Hung et al., Nature 2024](https://doi.org/10.1038/s41586-024-07861-8) | Correlation of different ecDNA species in daughter cells | Cross-species correlation is not within-species rho or PR60's R. The article describes correlated species alongside binomial marginal inheritance. |
| [Scanu, Werner and Huang, 2025 preprint version](https://doi.org/10.1101/2024.10.22.619675) | Computational ecDNA variant maintenance and de-mixing | No new sister-cell segregation measurement; assumed binomial allocation cannot validate itself. The preprint is not relabelled as a peer-reviewed experiment. |
| [Sankar et al., Nature 2026; online 2025](https://doi.org/10.1038/s41586-025-09764-8) | Retention elements, endogenous chromosome colocalization, introduced-plasmid hitchhiking failure | Retention/loss and fractional DNA signal are different estimands. Endogenous colocalization reuses Lange images; do not count it as independent replication or convert a signal fraction into a failed-cell count. |
| [Nichols et al., Molecular Cell 2025](https://doi.org/10.1016/j.molcel.2025.06.013), [Xie et al., Nature Cell Biology 2026](https://doi.org/10.1038/s41556-026-01982-0) | Mitotic clustering, chromosome attachment and cytosolic loss; differing mechanistic interpretations | Preserve live-imaging evidence and mechanistic disagreements. Nuclear loss violates the conserved-copy threshold model. No plot coordinates or cell-level arrays were extracted. |
| [Brueckner et al., Nature Genetics 2026](https://doi.org/10.1038/s41588-026-02782-w) | Real paired-daughter FISH and ABC inference of biased-versus-random segregation, Supplementary Fig. 12 | A relevant additional segregation study. The supplement plots empirical posterior densities/boxes but prints no usable empirical parameter estimates, interval endpoints or bin counts. Retain printed cell-line sample sizes; do not digitize posterior curves or mistake synthetic-validation parameters for measured fits. |
| [Wang et al., Research Square 2026 preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC13142634/) | Approximate unequal inheritance fitted jointly with population recovery and growth/death; sister imaging | Descriptive model fit, confounded with other dynamics and lacking usable parameter uncertainty. Imaging measures total DNA. Not an independently calibrated beta-binomial fit. |
| Historical work by Levan; Barker; Hamkalo; Jack; Takayama; Kanda; Tanaka; Shimizu; Pauletti; Kimmel and colleagues | Double-minute pairing, chromosome hitchhiking, retention, unequal allocation and branching models | Full details and access qualifications are in the audit. Neither clustering nor unequal segregation, nor an amplification branching model, is newly discovered here. |

For Brueckner's copy-conserving simulation mixture, if its numerical
parameters were available, R would be q(2p1-1)^2, where q is the probability
of biased allocation and p1 its allocation probability. This is not the
paper's micronuclear loss rate. Lack of printable fitted values, rather
than an assumption that the paper lacks segregation measurements, is the
reason no numerical R is imported from that source.

### Data files and reuse

| File | Contents and rights |
| --- | --- |
| [shibata_fits.csv](data/shibata_fits.csv) | Tiny unchanged fitted numerical aggregates and figure labels; full source names and exact spreadsheet cell/figure locators. Shibata article is **CC BY-NC-ND 4.0**; attributed unchanged numerical facts are retained for noncommercial research, with no article, figure, table layout or adapted source narrative redistributed. |
| [shibata_models.csv](data/shibata_models.csv) | Complete small published model-comparison table as numerical facts, with its signed-likelihood discrepancy preserved; same source rights. |
| [lundberg.csv](data/lundberg.csv) | Tiny unchanged quantitative facts; original article grants Creative Commons Attribution redistribution with attribution without specifying a version in its copyright paragraph. |
| [other_aggregates.csv](data/other_aggregates.csv) | Eligible Lange, Hung, Sankar, Wang and Brueckner aggregate facts, including negative/contextual measurements. Each row retains its own source-specific license/reuse basis. |
| [literature_audit.csv](data/literature_audit.csv) | Primary metadata and original inspection/eligibility prose. Where quantitative redistribution rights were not established, only metadata and links are included. |
| [analysis_design.csv](data/analysis_design.csv) | Clearly labelled original mathematical inputs, numerical tolerances and future-experiment design choices. **These are not published observations.** |

Every source row identifies its population/condition, source, DOI, locator,
access date and reuse basis; numerical rows also identify units. Shibata's
wide rows give separate locators/URLs for the fit sheet and sample-size
figure. Missing entries remain blank. The repository MIT license applies
to original code/documentation as appropriate and **does not relicense
third-party source facts or materials**. A reader seeking broader reuse
must consult each original source and its terms. No source workbook/PDF,
source figure, patient-related sheet, raw array, individual-cell row,
reconstructed cell list, or plot digitization is committed.

## Pre-specified plan and amendments

[PROTOCOL.md](PROTOCOL.md) was committed locally at
`2af64b4866b78dd8242ffa875418b90099dfaa5c` before this project's new numerical
fitting/threshold calculations, based on current main
`4f2098f543eb36d2810809af600c4372673794e7`. The first source inspection and
the prior notes already exposed some results; this was never blinded.
The protocol prefers a new aggregate likelihood or sufficient conditional
moments, permits explicitly conditional model reconstruction when justified,
and requires abstention where inputs are insufficient. No compatible
published aggregate likelihood was recovered for a new empirical refit.

The shell lacked Git push credentials, so publication used the connected
GitHub Git-data API. The [public protocol-preservation commit](https://github.com/dovvnloading/open-cancer/commit/fd427d4a91baf9a16fe98298ccb16c2ddfda2bd2)
preserves the same frozen content before the results commit, on main after
PR64. Its publication time is later than the local freeze; it is **not an
independent preregistration timestamp**. The local freeze hash above records
the actual computational sequence, not a claim of external registration.

**Computational amendments, 2026-11-17:**

1. After discovering the source observation-likelihood and sample-unit
   ambiguities, the expected-information calculation was implemented as a
   normalized hypothetical reconstruction. Hostile review of its first
   calculations found very weak squared-bias information. It was therefore
   **downgraded to post-freeze local Fisher sensitivity**, with no claimed
   sampling coverage. It is not the protocol's primary empirical analysis
   and is not presented as a recovered source confidence interval.
2. The exact-certificate width, rounding grid, outward sensitivity rounding
   and fixed-copy prediction design were selected after source inspection.
   They are recorded in the design CSV. Clone12 was selected for the future
   experiment because it is an ecDNA baseline clone with a larger inferred
   total-copy count and relatively precise local information. This selection
   is exploratory; the experiment and rejection rule remain unperformed.
3. Current primary text and rendered Fig. 3b consistently classify clones
   4, 6 and 7 as HSR. All three are excluded from primary ecDNA calibration;
   the earlier repo note's disputed-clone treatment is not carried forward.
   The prior directory remains unchanged.
4. An initial broad screening decision treated Brueckner principally as a
   nuclear-loss study. A hostile follow-up checked its actual supplement,
   recognized the genuine sister-pair ABC study, and corrected its audit
   entry. The missing printed posterior numbers are now the explicit limit.

## Model and interpretation

Let a type-k cell have division rate b and whole-cell death rate gamma*k.
It divides before dying with probability

$$p_k=\frac{1}{1+ck},\qquad c=\gamma/b.$$

At division k copies replicate exactly to N=2k. Draw an allocation
probability P independently for each mother, then J conditional on P is
Binomial(2k,P). The daughters are J and 2k-J. They are complementary,
not independent observations. Zero-copy cells cannot reacquire ecDNA.
The process has no nuclear loss, reintegration, immigration, variable
replication, density regulation, heritable environmental state or
copy-dependent division rate. These are assumptions, not measured facts.

For the source beta-binomial family,

$$P\sim\operatorname{Beta}(p(1/\rho-1),(1-p)(1/\rho-1)).$$

Here p is the unfurled Beta mean, constrained to at most one-half by
sister-label symmetry; it is **not** the mean minority fraction. For a
uniformly selected sister the effective allocation second moment is

$$R=E[(2P-1)^2]=(2p-1)^2+4p(1-p)\rho,$$

and the conditional daughter variance and Fano factor are

$$\operatorname{Var}(J\mid k)=Rk^2+(1-R)k/2,
\qquad F(k)=Rk+(1-R)/2.$$

This R is neither a marginal sister-count correlation across mothers of
different sizes nor a correlation between different ecDNA species.
The NB source fits describe the total counts of sampled divisions; they
are used only in the optional information reconstruction. The branching
process does **not** redraw its maternal count from that distribution.
Transporting the source segregation kernel to every maternal k, particularly
the low-copy states that control persistence, is an additional assumption.

The threshold c* separates certain eventual loss of the **entire
ecDNA-positive genealogy** from positive probability of its persistence.
It is not tumor extinction, an ecDNA-positive fraction, a finite-time cell
count, or a pharmacological dose. No experiment here measured gamma/b or
established linear whole-cell killing. The exponential-killing law in PR55
is a different model; its threshold is not numerically substituted for this
linear-killing threshold.

## Recomputed results

This block is generated by `analysis.py --markdown` and matched exactly by
the new test. All calculations use the checked-in aggregate/design CSVs.
The report deliberately separates published measurements, source-fit
transport, numerical sensitivity, uncalibrated information bands and
conditional future-experiment sampling variation.

<!-- BEGIN RECOMPUTED RESULTS -->
### Recomputed evidence inventory and calibration

The audit contains **31 literature records** (including explicitly labelled duplicates), **33 Lundberg scalar rows**, **41 other-study scalar rows**, **16 Shibata fitted-parameter rows**, and **6 Shibata model-comparison rows**. Analysis-design rows are original choices, not published observations.

Point statistics and probabilities are rounded decimal displays. Rigorous enclosure endpoints are displayed with outward rounding; JSON output also supplies exact rational endpoints.

The independently recomputed binomial threshold is enclosed by **[0.37847450, 0.37847525]**. These are exact-certificate endpoints, not a statistical interval.

Across the **9 ecDNA baseline clone fits**, R ranges from **0.089717 to 0.148457**, with descriptive median **0.104646**. Their conditional model thresholds are **11.19%–18.76% lower** than binomial. This is heterogeneity of published point fits in one selected cell-line system, not a population distribution of true R or a demonstrated treatment effect.

The threshold intervals below include numerical error and a declared half-last-decimal rounding box for both source parameters. **They exclude unknown sampling and model error.** The full asymmetric beta-binomial kernel is used in every row.

| Source fit | p | rho | Effective R | c* rounding/numerical enclosure | Relative shift (%) |
| --- | ---: | ---: | ---: | --- | --- |
| clone1 | 0.46894739 | 0.08842581 | 0.091942 | [0.335031, 0.335034] | [-11.4788, -11.4777] |
| clone2 | 0.45806315 | 0.11236728 | 0.118612 | [0.322076, 0.322080] | [-14.9018, -14.9004] |
| clone3 | 0.46315360 | 0.09852702 | 0.103423 | [0.329467, 0.329472] | [-12.9489, -12.9473] |
| clone5 | 0.44687637 | 0.13873426 | 0.148457 | [0.307467, 0.307471] | [-18.7617, -18.7604] |
| clone8 | 0.46673800 | 0.09108586 | 0.095108 | [0.333496, 0.333505] | [-11.8844, -11.8817] |
| clone9 | 0.47392428 | 0.10281156 | 0.105252 | [0.328581, 0.328585] | [-13.1830, -13.1817] |
| clone10 | 0.46822233 | 0.08602567 | 0.089717 | [0.336106, 0.336110] | [-11.1948, -11.1934] |
| clone11 | 0.44661171 | 0.13734416 | 0.147180 | [0.308094, 0.308098] | [-18.5960, -18.5947] |
| clone12 | 0.46388591 | 0.09995085 | 0.104646 | [0.328873, 0.328877] | [-13.1059, -13.1045] |
| condition1 (siCtrl_1, DMSO) | 0.42860890 | 0.17923980 | 0.195972 | [0.284045, 0.284052] | [-24.9502, -24.9481] |
| condition2 (siCtr1_1, JQ1) | 0.31385949 | 0.08771439 | 0.214151 | [0.277741, 0.277745] | [-26.6159, -26.6146] |
| condition3 (siCtrl_EN2, DMSO) | 0.44992150 | 0.13078506 | 0.139505 | [0.311859, 0.311863] | [-17.6013, -17.5999] |
| condition4 (siCtr1_EN2, JQ1) | 0.43883707 | 0.16426146 | 0.176767 | [0.293531, 0.293534] | [-22.4439, -22.4428] |

The perturbation arms are separate conditions from the same study; their drug names identify source experiments and do not calibrate a drug dose to c. No averaging across these kernels is interpreted as a mixed-population threshold.

### Uncertainty that can and cannot be propagated

The next table is a **post-freeze local-information sensitivity analysis, not empirical confidence intervals**. The R bands use the nominal normal multiplier from the design CSV and the normalized hypothetical observation model described below. Source sampling units, filtering, likelihood normalization, and finite-sample regularity are unresolved. The mathematical c envelope is uniform over every fixed allocation-mixture law with R inside the displayed band, using Theorem 1. It does not grant the input band statistical coverage.

| Source fit | Assumed independent divisions | Local R sensitivity band | Uniform c envelope | Relative shift envelope (%) |
| --- | ---: | --- | --- | --- |
| clone1 | 243 | [0.067977, 0.115907] | [0.209966, 0.466012] | [-44.53, 23.13] |
| clone2 | 54 | [0.072359, 0.164864] | [0.173044, 0.463821] | [-54.28, 22.56] |
| clone3 | 277 | [0.077453, 0.129392] | [0.199817, 0.461274] | [-47.21, 21.88] |
| clone5 | 247 | [0.119469, 0.177444] | [0.163521, 0.440266] | [-56.80, 16.33] |
| clone8 | 254 | [0.045494, 0.144722] | [0.188261, 0.477253] | [-50.26, 26.10] |
| clone9 | 240 | [0.053290, 0.157214] | [0.178828, 0.473355] | [-52.76, 25.07] |
| clone10 | 304 | [0.043465, 0.135970] | [0.194861, 0.478268] | [-48.52, 26.37] |
| clone11 | 134 | [0.073250, 0.221109] | [0.130353, 0.463375] | [-65.56, 22.44] |
| clone12 | 258 | [0.084728, 0.124565] | [0.203452, 0.457636] | [-46.25, 20.92] |
| condition1 | 50 | [0.105312, 0.286633] | [0.080238, 0.447344] | [-78.80, 18.20] |
| condition2 | 123 | [0.154588, 0.273714] | [0.090152, 0.422706] | [-76.19, 11.69] |
| condition3 | 74 | [0.084801, 0.194208] | [0.150808, 0.457600] | [-60.16, 20.91] |
| condition4 | 55 | [0.105459, 0.248075] | [0.109779, 0.447271] | [-71.00, 18.18] |

A diagnostic of failed local regularity: clone10 has squared bias B=0.004039, but its local B standard error is 159.63, larger than the entire admissible domain. This is why even a positive R band is **not used to reject R=0 empirically**. Every propagated envelope above contains the binomial reference; the direction of the shift is unresolved by these conservative bounds.

Without additional validated sampling/measurement assumptions, the aggregates supply no nontrivial calibrated confidence set for true R. Taking its full admissible range [0, 1] gives c in [0, 0.5] and a relative-shift enclosure **[-100.00, 32.11]%** under the branching assumptions. This is an outer bound, not a sharp identified set or a confidence interval; its positive upper endpoint does not prove that an increase is attainable.

### Published model comparisons and heterogeneity

All values below are from the source aggregate comparison across the full clone set, including HSR controls. AIC is a source-model comparison, not a p value or a probability that a biological mechanism is true. The source's negative-valued column headed 'Negative log-likelihood' behaves as a log-likelihood in its AIC calculation; the CSV preserves the signed values and the label discrepancy.

| Model | Parameters | Allocation model | Source AIC | Recomputed delta AIC |
| --- | ---: | --- | ---: | ---: |
| M1 | 25 | binomial; p common; rho zero | 155052.775 | 6892.480 |
| M2 | 36 | binomial; p clone_specific; rho zero | 154882.441 | 6722.146 |
| M3 | 26 | beta_binomial; p common; rho common | 148238.831 | 78.536 |
| M4 | 37 | beta_binomial; p clone_specific; rho common | 148160.295 | 0.000 |
| M5 | 37 | beta_binomial; p common; rho clone_specific | 148181.077 | 20.782 |
| M6 | 48 | beta_binomial; p clone_specific; rho clone_specific | 148182.080 | 21.785 |

The favored source model has clone-specific allocation means and common rho. The separate clone fits therefore do **not** establish heterogeneous rho. No compatible second study supplies a new likelihood for the same estimand, so a cross-study hierarchical dispersion estimate is not fit.

HSR controls are retained in the audit but excluded from ecDNA calibration:

| HSR control | Source rho | Transformed R (descriptive only) |
| --- | ---: | ---: |
| clone4 | 0.09730201 | 0.101896 |
| clone6 | 0.12486279 | 0.132606 |
| clone7 | 0.08655862 | 0.089595 |

Nonzero fitted dispersion in these chromosomal controls warns that the fitted parameter can absorb counting effects or model mismatch; it is not by itself a measurement of physical ecDNA clustering.

### Binomial-baseline challenges and opposing summaries

Lundberg's integer counts are derived by enumerating all counts compatible with the printed whole-percent rounding and published denominators. The exact tails and Holm adjustment assume calibrated per-division tests and independent sampled mothers; Wilson intervals have the same sampling qualification. The Markov bound dispenses with independence but still needs calibrated marginal tests.

| Line | Divisions | Compatible rejected counts | Conditional upper-tail p | Holm p | Wilson rejection-fraction interval | Dependence-robust Markov p bound |
| --- | ---: | --- | ---: | ---: | --- | ---: |
| CHP-212 | 100 | [2] | 0.26423802 | 0.52847604 | [0.005501, 0.070012] | 0.500000 |
| LA-N-5 | 100 | [10] | 7.6315875e-08 | 2.2894763e-07 | [0.055229, 0.174366] | 0.100000 |
| MC-IX | 50 | [1] | 0.39499393 | 0.52847604 | [0.003539, 0.104955] | 0.500000 |

The secondary pooled count is 13/250, with conditional p=1.9358638e-06 and Markov bound 0.192308. It is not an additional independent experiment or a fourth Holm-family member. This preserves PR63's verdict: conditional inconsistency, not an unconditional biological rejection and not an estimate of rho.

Lange's complete retained Table SI 1 is opposing evidence:

| Cell line / amplicon | Source KS distance | Source KS p | Source #samples |
| --- | ---: | ---: | ---: |
| PC3_Myc | 0.065 | 0.375 | 200 |
| SNU16_Myc | 0.039 | 0.918 | 194 |
| SNU16_fgfr2 | 0.063 | 0.415 | 196 |
| GBM39_EGFR | 0.072 | 0.221 | 210 |
| COLO_Myc | 0.033 | 0.973 | 206 |
| TR14_MYCN | 0.075 | 0.315 | 162 |
| TR14_CDK4 | 0.086 | 0.239 | 142 |
| COLO320_DM | 0.049 | 0.996 | 66 |

These source KS tests do not reject their simulated distribution. Complementary sister dependence, variable maternal totals and reported main-figure/supplement differences prevent recovering a calibrated rho confidence interval from them. Failure to reject does not establish zero overdispersion. No study or source panel is dropped for having an inconvenient sign.

### Single sharp future prediction

Use untreated **Clone12**, verified maternal pre-replication count **k=10**, and conserved post-replication total **N=20**. In **100 independent mother divisions**, count H, the number with a minority daughter receiving **at most 4 copies**. Count molecular copies and true zeros, not unresolved hubs. Eligibility is set from the pre-division mother; nuclear loss is recorded separately and invalidates the conservation model rather than being silently excluded.

| Quantity | Unbiased binomial | Source-fit transported prediction |
| --- | ---: | ---: |
| Probability of an extreme split | 0.011817932 | 0.164996542 |
| Expected H | 1.181793 | 16.499654 |
| Conditional central 95% predictive interval for H | Not used for the rejection rule | [10, 24] |

Reject the unbiased-binomial experiment model when **H >= 5**. Its exact null size is **0.006869524**, and its power under the fixed published-fit kernel is **0.999890254**. Conversely, H < 10 or H > 24 falsifies this transported point-kernel prediction at the stated two-tail predictive rule (tail probabilities 0.023536384 and 0.019350592). Both models can fail.

The same fixed-kernel prediction has conditional daughter Fano **1.494140091** and variance **14.941400913**. These are alternative summaries of this one experiment, not additional independent tests. The prediction interval and power condition on the published kernel; they do not include unknown parameter, assay, clone-transfer or culture-dependence uncertainty.

<!-- END RECOMPUTED RESULTS -->

## Theorems and verification

### Theorem 1: a uniform threshold envelope from an interval for R

Assume a fixed allocation-mixture law for P at every maternal type, the
copy-conserving branching assumptions above, and linear killing. Randomly
orient the sisters and put X=2P-1, R=E[X^2], T=E[X^4]. If a justified
interval contains R, say R_L <= R <= R_U, define z1(r), z2(r) as the
positive roots of

$$8c^2+(1+11r)c+r+4r^2-1=0,$$

$$8c^2+(1+14r-3r^2)c+7r^2-2r^3-1=0.$$

Set a root to zero when its constant term is nonnegative, and let
L(r)=min(z1(r),z2(r)). Then

$$\boxed{L(R_U)\le c_*\le (1-R_L)/2.}$$

This is a uniform identification envelope, not a formula identifying c*
from R. Its explicit lower bound is the proposed new mathematical result.
The upper bound is PR60's Theorem 5 and is credited to that note.

**Proof.** Jensen's inequality and |X|<=1 give R^2<=T<=R. The positive
two-type principal mean matrix, for retained types one and two, is

$$M_2=\begin{pmatrix}
\frac{1-R}{1+c}&\frac{1+R}{2(1+c)}\\
\frac{1-T}{2(1+2c)}&\frac{3(1-2R+T)}{4(1+2c)}
\end{pmatrix}.$$

Multiplication by the positive denominator shows that det(I-M2) has the
sign of

$$D(c,R,T)=(c+R)(1+8c+6R-3T)-(1+R)(1-T).$$

This is affine in T; its maximum over [R^2,R] is at an endpoint. At
T=R and T=R^2 it equals the two polynomials defining z1 and z2.
For c>=0 and R in [0,1], their derivatives with respect to R are
11c+1+8R and (14-6R)c+14R-6R^2, both nonnegative. Thus c<L(R_U)
makes both polynomials negative for every R in the interval, and makes
det(I-M2)<0. A nonnegative two-type matrix then has Perron root above
one, certifying persistence in a finite subprocess and hence the full
process. For R<1 the allocation law has positive mass in the interior,
so positive types communicate; for R_U=1 the lower bound is zero and
the persistence statement below it is vacuous. Taking the supremum of
c values with certified persistence gives the claimed lower threshold bound.

For the upper bound, the probability A_k of two positive daughters obeys
A_k<=k(1-R)/2 (PR60). Hence the mean positive offspring count is at
most one when c>=(1-R)/2. The uniform positive death floor and the
finite-population visitation argument in PR60 prove extinction, including
the boundary. Replacing R by its lower endpoint proves the uniform bound.
Complete clumping has c*=0 for c>0, while death-free c=0 preserves a
positive lineage. These cases are consistent with the envelope. QED.

**Coverage corollary.** If an R interval has valid coverage at least
1-alpha under an explicit statistical model, this uniform transformation
has at least the same coverage for c* under the branching assumptions.
This follows because interval containment for R implies threshold
containment pointwise. A deterministic proof cannot confer coverage on
an unvalidated local Fisher band. The broad bounds in this report are
conservative; their failure to prove a shift does not demonstrate that
all values inside the envelope are attainable.

### Theorem 2: certificates control every copy number

Let K_k(j) be the ordered beta-binomial mass and

$$M_c(k,j)=\frac{K_k(j)+K_k(2k-j)}{1+ck},\qquad j\ge1.$$

PR60's compact-tail criterion applies for interior beta parameters and
c>0, including critical extinction. At a finite head size m, choose
positive integer weights w1,...,wm. If every head row of the principal
matrix satisfies M_c w>w, its Perron root exceeds one, proving persistence
in the retained finite subprocess. This supplies a lower threshold endpoint.

For an upper endpoint, use w_k=t>0 for every k>m. Check every head row,
including all its exits above m, satisfies M_c w<=w. Also check

$$\frac{2\max(t,\max_{j\le m}w_j)}{1+c(m+1)}<t.$$

Every omitted row satisfies this inequality because its mean positive
offspring count is at most 2/(1+ck). The bounded weights have a positive
minimum, so expected population size remains bounded. There is a uniform
death probability at least c/(1+c). For any finite population bound B,
the probability of immediate extinction when the population is between
one and B is bounded below by that death probability to the power B.
Repeated visits therefore cannot occur forever on survival; population
size tends to infinity on survival. Fatou's lemma contradicts bounded
expected size unless survival probability is zero. Strict tail slack and
PR60's positive left eigenvector give a strict upper threshold endpoint.

For rational source inputs, every K row has a common integer denominator:
if Beta shapes are A/D and B/D and their sum is S/D, then

$$K(j)=\binom{n}{j}
\frac{\prod_{i<j}(A+iD)\prod_{i<n-j}(B+iD)}
{\prod_{i<n}(S+iD)}.$$

The verifier proposes weights using floating arithmetic, rounds weights to
integers, and accepts them only after exact integer residual comparisons
and the all-type tail inequality. Floating error can cause a failed
proposal; it cannot make an invalid inequality pass. The numerical head
size is not a biological copy cap.

### Uniform rounding and parameter-box propagation

The point kernel uses the published decimals as exact inputs. Separately,
the rounding sensitivity allows each parameter to vary by half its last
nonzero printed decimal unit. This is a declared rounding scenario, not
author-supplied estimation uncertainty. Positive rising-factorial bounds
on Beta shape numerators and the concentration denominator enclose every
kernel throughout each small box; exact lower/upper certificate checks
then enclose every threshold in that box.

For broader boxes, the optional `rectangle_threshold` routine uses the
Polya-urn factor

$$f(q,r;i,t)=\frac{q(1-r)+ir}{1-r+tr}.$$

It increases in q; for fixed q its derivative in r has constant sign
i-qt. Its extrema on a rectangle are therefore at corners. Products of
factor minima/maxima bound each kernel entry, which in turn bounds its
mean operator. These envelopes need not sum to one: they are matrix
bounds. Subdivision followed by the same exact head/tail checks supplies
uniform enclosures, with the universal bound as fallback. **No monotonicity
of the true threshold in p or rho is assumed.**

### Lemma: folded observations do not globally identify bias and rho

For a>0, symmetrizing Beta(a,a+1) gives Beta(a,a): factor the densities
and use x+(1-x)=1 and B(a,a+1)=B(a,a)/2. Equivalently, for any interior
R, parameter pairs

$$ (p,\rho)=(1/2,R),\qquad
(p,\rho)=((1-R)/2,R/(1+R)) $$

generate identical unordered sister-count laws for every total count.
They have the same R and the same threshold. Thus even perfect unordered
counts do not globally identify the bias/rho decomposition. This does
not imply R is unidentifiable from sufficiently rich conditional moments.
Conversely, equal R need not imply equal full kernels or thresholds:
PR60's beta-versus-fixed-bias example already demonstrates that distinction.

## Why the optional information calculation is only sensitivity

For a transparent hypothetical experiment, let the total N follow the
positive-truncated NB distribution with the published size and mean.
Draw J conditional on N from the published beta-binomial fit. Observe
the sorted pair (max(J,1),max(N-J,1)). The probability of observing
(1,m) combines latent outcomes (N,J)=(m,0) and (m+1,1), together with
reflections; ties are counted only once. This construction normalizes.

It differs from a literal reconstruction of Supplementary Note 3 in
material ways: its beta-binomial zero-floor branch prints inconsistent
totals, and its likelihood equations omit the positive-NB normalization
used in its simulations. We resolve these for the declared hypothetical
model, **not** by asserting what the source code actually computed.
The source's image-analysis method also removes values above each group's
upper percentile cutoff, whose effects cannot be recovered from these
aggregates. No source algorithm is silently repaired and then relabelled
an empirical reanalysis.

Analytic probability derivatives are accumulated in coordinates
(R, B=(2p-1)^2, log NB-size, log NB-mean). Summing score outer products
gives the full expected Fisher matrix, including nuisance parameters and
the zero-floor aliases. The inverse matrix's R entry, divided by the
assumed independent pair count, supplies a local scale. The nominal
normal multiplier makes a sensitivity band, then its endpoints are
rounded outward before applying Theorem 1. The calculation reports its
NB mass truncation and tests normalization and score identities.

Baseline plotted cell counts are divided by two only as a reconstruction
assumption. Perturbation minority-fraction counts are treated as pair
counts only as an assumption. Neither assumption establishes independent
cultures or effective sample size. Very weak information in B and the
exact aliases above defeat any claim of demonstrated finite-sample Wald
coverage. No biological null test is based on a positive endpoint of this
band; no synthetic sample is used as if it were an observed experiment.

## Falsification attempts and hostile review

| Attempted claim or shortcut | Result and consequence |
| --- | --- |
| Use the published R range as a population confidence interval | Rejected. These are selected clone point fits, not independently sampled cancer types or a sampling distribution. |
| Refit a beta-binomial likelihood from fitted parameters alone | Rejected. The original likelihood surface and sufficient observations are absent. The code labels source fits as source fits. |
| Interpret the local Fisher bands as empirical confidence intervals | Rejected after information diagnostics, missing sampling units and observation-model discrepancies. Bands are secondary sensitivity only. |
| Infer established rho heterogeneity from separate-clone point estimates | Rejected. The complete source model comparison favors common rho with varying allocation means; the comparison also includes HSR controls. |
| Include the disputed clone as ecDNA | Rejected by the current primary text and rendered metaphase classification; all HSR controls are retained separately. |
| Treat every newer clustering paper as a direct rho estimate | Rejected. Cross-species correlation, interphase hubs, nuclear loss and conditional same-species allocation are different estimands. Brueckner's actual paired-daughter supplement was specifically rechecked. |
| Digitize an ABC posterior plot to obtain missing uncertainty | Rejected by the aggregate-only scope. No printed numeric posterior summary was recovered. |
| Treat Lundberg's rejection fractions as literal clustering rates or rho | Rejected. Their conversion requires count totals, calibrated tests and an observation model. Dependence-robust sensitivity does not reject. |
| Ignore Lange's nonsignificant comparisons | Rejected. Every table row is retained. Differences in copy ranges, assays and simulations prevent treating this as a pooled estimate of the same R. |
| Assume a constant allocation law across copy numbers | Challenged by Lundberg's association of clustering with high copy counts and selection of high-count mitoses for cluster inspection. Constant R is a falsifiable transport assumption. |
| Infer eradication from a subcritical finite matrix | Rejected. All upper certificates include exact head exits and an inequality covering every omitted copy number. |
| Use c*=(1-R)/2 or compare only endpoints of an R interval | Rejected. This is an upper bound; R does not identify the kernel, and general threshold monotonicity is unproved. Uniform envelopes replace interpolation. |
| Treat the four perturbation arms as a dose-response test of linear killing | Rejected. They measure segregation changes, not independent gamma/b calibration. |
| Claim a first model of ecDNA clustering or unequal inheritance | Rejected by historical double-minute biology, prior amplification branching models, Hung, Nichols, Shibata and related sources. |
| Infer eventual extinction from a short surviving-well assay | Rejected. The chosen sharp prediction is a finite, directly observable sister-count experiment with explicit model-falsification outcomes. |

Independent source checks re-read the spreadsheet scalars and rendered figure
labels, including the rights statement. Mathematical checks independently
rebuilt beta-binomial masses, differentiated a log-gamma implementation of
the observation likelihood, checked zero-floor tie/alias accounting, and
tested false-certificate rejection. These checks provide computational and
transcription scrutiny; agreement among AI agents is not independent
biological validation or external scientific review.

## Novelty assessment

The biological observation of unequal or clustered double-minute inheritance
is old. Beta-binomial segregation inference is already in Shibata; ecDNA
species co-inheritance and tethering mechanisms are already in the cited
literature. Amplification branching models, variable segregation and
deamplification predate modern ecDNA terminology, including
[Kimmel and Axelrod 1990](https://doi.org/10.1093/genetics/125.3.633) and
[Kimmel, Axelrod and Wahl 1992](https://doi.org/10.1016/0165-1110(92)90010-7).

The contribution here is narrower: traceable transport of published fitted
kernels into PR60's infinite-type killing model; efficient exact
certificates for those source-fit thresholds and their rounding boxes;
the explicit uniform two-type lower envelope paired with the existing
upper bound; an observation-model/identifiability audit; and a sharp
prospective count experiment. No located source supplied this exact
combined calibration and envelope. That is a provisional targeted novelty
assessment, not a proof of priority. The Beta identity is elementary and
not claimed as a newly discovered general statistical identity.

## Limitations and what would resolve them

The most useful additional publication would supply, for each independently
sampled culture and maternal total, **already-aggregated** counts of unordered
daughter pairs or sufficient conditional second/fourth moments, with the
number of independent divisions, zero handling, imaging error calibration,
filtering and lost-copy outcomes. Published profile likelihoods or bootstrap
replicate summaries with joint parameter uncertainty would also help.
These can preserve the aggregate-only and privacy boundaries; no patient
records are required.

More cohorts are required before between-study heterogeneity can be
estimated. Fluorescence foci may contain multiple molecular copies.
Compensating replication, micronuclear sequestration, chromosome retention,
lineage memory and copy-dependent division can all invalidate the
transported kernel or killing law. Accurate estimates in a culture would
still not establish the therapy model. Neither a simulation, source-fit
threshold nor passed CI supports clinical advice.

## Reproduce and review

From the repository root, using only the Python standard library:

```sh
python research/ecdna-segregation-clustering-data/analysis.py --json
python research/ecdna-segregation-clustering-data/analysis.py --check-readme
python -m unittest discover -s tests -p 'test_research_ecdna_segregation_clustering_data.py' -v
```

`analysis.py --markdown` regenerates the exact numerical report block;
`--design` prints only the optional local-information sensitivity details.
The test recomputes every displayed source/derived result from CSVs,
checks source-row provenance and permitted scope, verifies independent
kernel identities and exact head/tail certificates, and exercises the
uncertainty and prediction calculations. It requires no network or external
packages and is designed to finish within the requested runtime budget.
The source policy and existing files/workflows are unchanged. The existing
pull-request validation workflow discovers the new test automatically.

Only this new directory and the new test file are part of the change.
This work remains **nonclinical and unreviewed**, regardless of CI status.
