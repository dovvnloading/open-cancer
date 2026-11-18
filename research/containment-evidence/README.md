# Containment theory confronted with published aggregate evidence

**Nonclinical literature computation. Research status: unreviewed. No treatment
advice, patient-level data, new experiment, or claim of clinical validation.**
Search/extraction date: **17 November 2026 (America/New_York)**. This work starts
from `main` at `4f2098f543eb36d2810809af600c4372673794e7` and adds only this
directory and its dedicated test. The source studies include clinical research;
this repository analysis is not a clinical decision tool.

## Finding

The published aggregates support **context-dependent outcomes under adaptive
and intermittent treatment**, and expose a decisive identification gap:
**they do not locate a real tumour above or below PR #61's mutation-supply
threshold**. Favourable prostate observations, negative randomized melanoma
point estimates, beneficial and neutral breast xenograft comparisons, and an
adverse lung xenograft comparison all remain in the evidence map. Numerical
effects, their published uncertainty, source links and exclusions are
recomputed in the tables below. A favourable treatment comparison cannot
identify the mechanism that produced it.

The [tumour-scale theorem](../containment-at-tumour-scale/README.md) predicts
an **ordinary uncensored mean**, under specific population, mutation and
treatment assumptions. The available median progression times, survival
hazards, PSA fits and animal humane-endpoint summaries generally target
different quantities. There is no justified numerical estimate or confidence
interval for a real tumour's `lambda/lambda_star` in this extraction. Missing
parameters are **unidentifiable**, not zero. This is an evidence-grounded
negative result about identifiability, not a biological confirmation or
refutation of the theorem.

A second correction matters before any experiment is interpreted: a globally
increasing containment advantage with **competition strength** is not a
prediction of PR #61. The theorem proves monotonicity in **mutation supply
with the other parameters fixed**, and a separate no-competition result.
Competition affects both suppression and establishment; the derivation below
shows why its overall effect need not be monotone.

## Provenance, search coverage and reuse

The immutable [protocol](PROTOCOL.md) was committed locally as `2d46668` before this
evidence extraction and the analysis. That local audit SHA is not a remote
GitHub ref: authenticated publication recreated the Git commits after the
command-line push lacked credentials, as recorded in AMENDMENTS.md.
This is prospective specification of
this computation on already published papers, **not** prospective registration
of the experiments or an analysis blinded to familiar literature. Later
clarifications and source corrections are recorded in [AMENDMENTS.md](AMENDMENTS.md).

Searches used two web search indexes, primary publisher/author manuscripts,
the official Europe PMC article service, trial registries and sponsor records.
`data/searches.csv` records the executed queries, search date, engine and
screened candidates; it is an audit log, not a reproducible snapshot of a
changing ranked search index. Seeds were all studies requested in the issue,
then references and related/citing comparisons, newer experiments, fitted
cohorts, trial results and negative intermittent-treatment comparisons.
No login barrier, CAPTCHA or controlled-data access was circumvented.

The operational scope is adaptive/evolutionary/containment treatment and
intermittent or dose-modulated comparisons relevant to competitive release,
resistance or treatment holidays. Generic dose finding and tolerability-only
dose reductions are outside this scope. We include broader intermittent
comparisons as **context**, not as equivalent implementations of containment.
This is a search-audited evidence map, **not a complete worldwide systematic
review**: subscription sources, unresolved candidates, conference updates,
uninspectable registry results tabs and reuse exclusions remain. The user's
requested universal census cannot honestly be certified from this search.
Every located candidate is retained with its disposition; inability to
extract is not evidence of no benefit or no publication.

Every CSV row has a source URL, locator, access date and reuse basis. Numeric
rows are small published aggregates from sources with an explicit
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) grant. Attribution
is through the study title, author-year identifier, source and locator.
Changes consist of selecting aggregate facts, transcribing them to CSV,
standardizing units/labels and computing the explicitly described contrasts.
These source materials retain their original terms; the repository's MIT
license does not relicense them. Full articles, figures and underlying
datasets are not redistributed. No patient/animal-level trajectories, subject
identifiers, individual fitted parameters, raw genomes or reconstructed
survival records are included or used in the analysis.

Some seminal sources have an author manuscript without a clear reuse grant;
others use noncommercial/no-derivatives terms whose applicability to a derived
CSV was not resolved. Under [DATA_POLICY.md](../../DATA_POLICY.md), those
records remain metadata and links only. This conservative choice is a **reuse
limitation**, not a claim that every numerical fact is copyrightable or that
those studies are scientifically ineligible. In particular, Gatenby,
Enriquez-Navas, Smalley, S1320 and West are visibly screened rather than silently
omitted. The study ledger also distinguishes **explicitly no posted results
observed** from **results not located with incomplete access** for registered
trials. A registered protocol is not an efficacy result.

### Data files

| File | Meaning |
| --- | --- |
| `data/studies.csv` | Source/cohort identity, design, setting, inclusion or exclusion, rights and limitations. Reports sharing a cohort are not replications. |
| `data/aggregates.csv` | Published arm summaries, HRs and intervals, reported exposure, and labelled group-level fitted/model quantities. Empty values are missing or not reached, never zero. |
| `data/comparisons.csv` | Explicit links between compatible arm summaries or a published HR, endpoint class, design, overlapping-control group, observed versus model role. |
| `data/mapping.csv` | Five parameter statuses per source, with reasons and any conditional transformation. A fit is not relabelled a direct measurement. |
| `data/gates.csv` | Study-specific eligibility for each mechanistic test; source-linked reasons for failure or uncertainty. |
| `data/searches.csv` | Executed search queries and candidate trail, including negatives, rights and access gaps. |

`exact` in a numeric row means an exact transcription of a **reported point
summary**, not exact biological truth. Source rounding remains. `at_least`,
`less_than`, `approximate` and `not_reached` preserve the source's qualification.
The source-labelled confidence level is retained. IQRs, ranges, standard
deviations, confidence intervals and missing uncertainty are distinct.

## Mapping the biological studies to PR #61

With `g0=bR-dR`, `gC=g0-alpha*s0`, `delta=kmax-rS` and
`H_N=L log(N)`, the theorem uses:

| Parameter | Definition and units | What must be measured independently | What the published aggregate usually supplies |
| --- | --- | --- | --- |
| `d` | `1/g0`, time | Resistant net expansion after competitive release, in the same system and units | A PSA/net-growth fit or no resistant growth estimate. A fitted number is conditional on a different observation/dynamics model. |
| `c` | `min(L,1/gC)` in the admissible positive-growth regime; `L` if `gC<=0`, time | Resistant invasion growth at the maintained sensitive density and the common course scale | Total burden control or a fitted competition coefficient; not a fixed sensitive reservoir and finite common `L`. |
| `J` | `L*(1+delta/rS)`, time | Sensitive net growth, realized maximum-dose killing and post-course rebound under a common finite course | PSA decline or dose fraction; no transferable kill/rebound decomposition. |
| `A` | `integral bS*sM(t)*pM(t) dt`, dimensionless | Sensitive divisions, time-varying density, resistant birth/death turnover and establishment | Net growth alone, which does not separate births/deaths or identify establishment. |
| `lambda` | `N*mu`, dimensionless | Absolute viable tumour-cell scale and acquisition of the relevant heritable resistance per completed sensitive division | Relative PSA/volume, existing resistance or mutation burden; none is `N*mu`. |

The cell-number normalization must be consistent. A carrying capacity in
arbitrary PSA units is not a measured count of viable cells. In particular,
the fitted eLife model uses a growth multiplier as well as intrinsic PSA-slope
estimates. The conditional inverse-growth calculation in the tables
includes that multiplier; treating the reciprocal of the reported resistant
slope alone as `d` would be wrong even within that source model. The conditional
calculation has no reported transferable CI and cannot identify the threshold.

Turnover causes a structural identification problem, not merely wide error
bars: adding the same positive rate to resistant birth and death preserves
`g0` while changing establishment probability. Increasing sensitive birth
and death together preserves `rS` while changing mutation supply per unit
time. Thus the same net-growth summaries can imply different `A`. Neither
PSA, mutational burden, initial resistant fraction, nor a source model's
parameter also named `lambda` repairs this missing information.

The eligibility gates additionally require the theorem's zero initial
resistance, irreversible mutation per division, absence of direct drug effect
on R, common bounded kill law, sensitive-reservoir feedback, common finite
course and exposure ceiling, and total-burden progression followed after
withdrawal. Pre-existing resistant mixtures, repeated treatment cycles,
drug dependence, phenotype switching, immune effects, PSA response, overall
survival and humane euthanasia can each change the comparison. Failure of a
gate does not make a study uninformative; it limits which claim it can test.

### Why the competition-strength extension fails mathematically

For fixed `lambda`, write

\[
\Delta=(c-d)-(J-d)e^{-\lambda A}.
\]

Holding the other demographic and treatment parameters fixed, PR #61's
establishment integral decreases with `alpha`: `A'(alpha)<0`. Therefore

\[
\frac{d\Delta}{d\alpha}
=c'(\alpha)+(J-d)\lambda e^{-\lambda A} A'(\alpha).
\]

The two terms have opposing signs. In the admissible saturated regime
`c=L`, the first is zero and the derivative is negative. Stronger competition
then reduces this asymptotic advantage and raises the supply threshold.
This is an algebraic correction to an **extension** of the theorem, not an
empirical falsification or a new fitted biological result. It leaves intact
the separate result that zero competition has no containment-winning region
under the theorem's assumptions. The dedicated test checks the sign logic;
no synthetic numerical example is presented as published data.

## Pre-specified tests and interpretation

The protocol fixes the mutation-supply trend (P1), competition trend as an
auxiliary hypothesis (P2), zero-competition result (P2a), quantitative boundary
(P3), establishment law (P4), one-crossing prediction (P5), and descriptive
transport checks (P6). Each mechanistic test first checks whether its required
independent exposure/rate measurements and endpoint exist. Cross-cancer labels
are not substituted for a controlled mutation or competition series.

`consistent` means the observed aggregate direction is compatible with the
stated claim in its limited scope, not validation. `falsified` requires the
specified contradictory evidence in that scope. `untestable` means the claim
cannot be adjudicated from these allowed aggregates: required variables,
comparable endpoints or sufficiently informative contrast uncertainty are
missing. It is not evidence that the theorem is true. No p-value is assigned
to an untestable prediction.

P6b's universal phrase is ambiguous about overall survival versus tumour
control and does not prescribe one confidence level. We therefore report
**both readings**: the original broad outcome/any-published-interval rule is
met by COIN's adverse OS interval, while explicitly **post-extraction**
sensitivities restricted to tumour control, and then to published 95% tumour
control intervals, are inconclusive. The broad verdict is a descriptive
contradiction at the source's 80% level, **not** a conventional 95% rejection,
a tumour-control result or a falsification of PR #61. The ambiguity prevents
a unique fully pre-specified P6b verdict. Failure of noninferiority alone is
not proof of harm. The original negative result is retained rather than
replaced by the stricter sensitivities; see AMENDMENTS.md.

No new pooled efficacy estimate, meta-regression, inferred hazard ratio or
multiple endpoint significance search is performed. Source intervals are
shown descriptively, including null-spanning intervals and shared controls.
The planned Holm adjustment applies only if mechanistic tests pass the data
gates. Aggregate median ratios are descriptive ratios of medians, not ratios
of individual benefit or estimates of the theorem's mean coefficient. A
source's arm-median CIs are not a CI for their ratio. Animal group counts and
published fit summaries are not made into patient observations.

## Contradictions and hostile source audit

The review checked the numerical claims against primary sources, independently
recomputed all transformations, and attempted to overturn favourable
interpretations. Corrections are recorded with source locators in
[AMENDMENTS.md](AMENDMENTS.md); the final test guards the important mistakes.

- The initial prostate report is a lower-bound interim result, and its
  reported comparator median differs between updates. The reports share
  a trial; fitted reanalyses of that trial do not add clinical replications.
- In the eLife update, the reported survival HR is for **overall survival**,
  not time to progression. Its Results define an average dosing rate across follow-up and
  complementary time off; a Discussion sentence reverses the time-off value.
  We retain the explicit Results definitions and flag the internal discrepancy.
  An average dosing rate is not total lifetime cumulative exposure.
- Randomized intermittent melanoma comparisons have adverse point estimates;
  the available HR interval spans the null and the other study reports
  arm-median intervals. Neither is discarded or converted into a fabricated
  precise adverse effect. Fixed schedules are not clonal-feedback containment.
- COIN's source-reported HR interval has a lower confidence level, while its
  median-arm brackets are IQRs. They are never relabelled as 95% CIs. Its
  endpoint is overall survival. Negative S1320 and BFR14 sources also remain
  linked despite the numeric reuse restriction.
- The breast experiment's neutral gemcitabine result is retained even where
  no numerical HR was printed. Its multiarm comparisons share controls;
  apparent repeated benefit is not independent replication. Exposure versus
  survival associations may also reflect different survival/accumulation windows.
  The authors explicitly warn that the constant combined GEM/CAP comparator
  may have exceeded tolerable dosing and that toxicity could explain much of
  the multidrug advantage. This is not an isolated competition effect.
- The newer lung xenograft's adverse mean survival contrast is retained. Its
  reported omnibus test is not a pairwise CI. Failure to detect density
  dependence does not establish `alpha=0`; the fitted growth model is not
  an independent competition-removal experiment.
- Ovarian adaptive-arm medians that were not reached/undefined remain missing.
  Neither the experimental stopping horizon nor follow-up is substituted for
  a median, and no Kaplan-Meier curve is digitized. In that study, standard
  treatment ended after a fixed course while adaptive dosing continued; the
  authors explicitly could not separate evolutionary benefit from increased
  cumulative exposure. This cannot establish a competition mechanism.
- Model-derived improvements in fitted cohorts are separate from observed
  comparisons. Initial resistance, stem renewal, switching and arbitrary
  burden scales are not silently recoded as PR #61's de novo mutation model.

These are reproducibility and interpretation checks. AI assistance, repeated
calculation, a checksum and passing CI do not establish source truth,
independent scientific review or biological validity.

## The single most decision-relevant open measurement

**Measure the integrated supply of newly established resistant lineages during
maximum-dose depletion**, independently of the comparative treatment outcome:
`theta=lambda*A`. Matched nonclinical replicate systems with common parameters and independently
checked sensitive starting lineages could report the aggregate fraction `p0`
without an established new resistant lineage; the model gives
`theta=-log(p0)`. An aggregate binomial interval can be transformed with reversed
endpoints. With heterogeneous supplies, pooled `p0=E[exp(-theta_i)]`, so
`-log(p0)` is not the mean supply. This avoids treating tumour mutational burden as a mutation rate
or inferring high supply from a favourable adaptive outcome.

The decision boundary can then be written

\[
\theta > \log\!\left(\frac{J-d}{c-d}\right),
\qquad d<c<J.
\]

The measurement needs a pre-defined lineage-establishment criterion, an
explicit observation horizon and bounds on undetected initial resistance.
A finite-course zero-lineage fraction is not automatically the theorem's
infinite-depletion-limit probability. This is the most discriminatory missing
measurement **conditional on** subsequently identifying `d,c,J` and testing
the growth/competition and finite-size assumptions; it cannot by itself make
the threshold usable. No such measurement was collected or commissioned here.

## How to rerun

From the repository root, using Python's standard library only:

```sh
python research/containment-evidence/analysis.py --check
python research/containment-evidence/analysis.py --json
python -m unittest discover -s tests -p test_research_containment_evidence.py -v
```

To update the generated block after an explicitly documented source correction:

```sh
python research/containment-evidence/analysis.py --write-readme
```

`analysis.py` makes no network calls. The dedicated test checks every input
join, provenance field, numerical contrast, missingness rule, cohort role,
parameter mapping, verdict and generated README table. The test includes
independent arithmetic and normalized-content checks; all file hashes use
UTF-8 content with CRLF and CR converted to LF, so Windows and Linux agree.
The test is bounded well below the requested one-minute ceiling. Re-running
does not recheck a publisher's current source or permissions: source auditing
is a separate requirement, documented in the row locators and audit ledger.

## Limits that remain

The search has no guarantee of worldwide completeness, and restrictions on
reuse/access limit the quantitative subset. Small animal experiments,
nonrandomized clinical controls, responder selection, censoring, unequal
exposure windows, model-dependent fits and publication bias limit causal
inference. A study's reported uncertainty need not cover protocol selection,
measurement error, structural model error or cross-system transport.

Most importantly, the theorem's large-population limit is not an exact
finite-tumour decision boundary. PR #61 itself retains near-boundary
finite-size discrepancies. A real experimental disagreement would need an
independently parameterized finite-size prediction and matched assumptions
before it could refute the model. No synthetic values from PR #61 are used
to fill a biological parameter gap here.

<!-- BEGIN RECOMPUTED EVIDENCE -->

## Recomputed evidence tables

The ledger contains **69 source/study records**, **15 reports with licensed aggregates**, **179 aggregate rows**, **44 contrasts**, and **165 search records. These are extraction counts, not independent studies or participants.**

Disposition counts: context_only=4; excluded_scope=6; included=10; included_model=4; included_single_arm=1; metadata_rights_or_access=26; model_only=6; no_posted_results_observed=1; overlap_or_review=3; results_not_located=5; unresolved_candidate=3.

### Observed outcome contrasts

Differences and ratios compare published arm summaries; they are not hazard ratios or uncensored means. No CI for a difference/ratio is invented. `>=` is a lower bound; `not reached/undefined` is not a numerical median. Published HRs retain experimental/reference orientation and confidence level. Earlier reports from the same cohort are explicitly marked. Shared-control rows are dependent.

| Contrast | Endpoint / units | Experimental | Reference | Difference | Ratio | HR | HR interval | Role | Source / input rows |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| zhang2017_ttp | Radiographic TTP (months) | >=27 | 14 | >=13 | >=1.92857 | not reported | not reported | superseded_same_cohort | [zhang2017](https://jessica-cunningham.com/wp-content/uploads/2019/11/cunningham_et_al-2017-nature_communications.pdf) / `z17_ttp_a;z17_ttp_c` |
| zhang2022_ttp | Radiographic TTP (months) | 33.5 | 14.3 | 19.2 | 2.34266 | not reported | not reported | observed_latest | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) / `z22_ttp_a;z22_ttp_c` |
| zhang2022_os | Overall survival (months) | 58.5 | 31.3 | 27.2 | 1.86901 | 0.41 | [0.2, 0.83] (95% CI) | observed_latest | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) / `z22_os_a;z22_os_c;z22_os_hr` |
| gonzalez_cao_pfs | PFS (months) | 6.9 | 16.2 | -9.3 | 0.425926 | not reported | not reported | observed_latest | [gonzalez_cao2021](https://www.nature.com/articles/s41467-021-26572-6) / `gc21_pfs_a;gc21_pfs_c` |
| gonzalez_cao_orr | ORR (fraction; rounded) | 0.6 | 0.714 | -0.114 | 0.840336 | not reported | not reported | observed_latest | [gonzalez_cao2021](https://www.nature.com/articles/s41467-021-26572-6) / `gc21_orr_a;gc21_orr_c` |
| interim_pfs | PFS (months) | 8.5 | 10.7 | -2.2 | 0.794393 | 1.39 | [0.79, 2.45] (95% CI) | observed_latest | [interim2024](https://eprints.whiterose.ac.uk/id/eprint/206919/7/EJC-D-23-02078_R1.pdf) / `int24_pfs_a;int24_pfs_c;int24_pfs_hr` |
| interim_os | Overall survival (months) | 18.1 | not reached/undefined | not estimable | not estimable | 1.69 | [0.87, 3.28] (95% CI) | observed_latest | [interim2024](https://eprints.whiterose.ac.uk/id/eprint/206919/7/EJC-D-23-02078_R1.pdf) / `int24_os_a;int24_os_c;int24_os_hr` |
| interim_orr | ORR (fraction; rounded) | 0.57 | 0.77 | -0.2 | 0.74026 | not reported | not reported | observed_latest | [interim2024](https://eprints.whiterose.ac.uk/id/eprint/206919/7/EJC-D-23-02078_R1.pdf) / `int24_orr_a;int24_orr_c` |
| coin_os | Overall survival ITT (months) | 14.4 | 15.8 | -1.4 | 0.911392 | 1.084 | [1.008, 1.165] (80% CI) | observed_latest | [coin2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3159416/) / `coin11_os_a;coin11_os_c;coin11_os_hr` |
| seyedi_cap_dm_mtd | OS: CAP dose modulation / CAP MTD | not reported | not reported | not estimable | not estimable | 0.22 | [0.043, 1.1] (95% CI) | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_cap_dm_mtd` |
| seyedi_cap_int_mtd | OS: CAP intermittent / CAP MTD | not reported | not reported | not estimable | not estimable | 0.78 | [0.236, 2.6] (95% CI) | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_cap_int_mtd` |
| seyedi_cap_dm_none | OS: CAP dose modulation / no treatment | not reported | not reported | not estimable | not estimable | 0.24 | [0.042, 1.4] (95% CI) | contextual_or_secondary | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_cap_dm_none` |
| seyedi_cap_int_none | OS: CAP intermittent / no treatment | not reported | not reported | not estimable | not estimable | 0.96 | [0.247, 3.7] (95% CI) | contextual_or_secondary | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_cap_int_none` |
| seyedi_cap_dm_int | OS: CAP dose modulation / CAP intermittent | not reported | not reported | not estimable | not estimable | 0.26 | [0.049, 1.4] (95% CI) | contextual_or_secondary | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_cap_dm_int` |
| seyedi_cap_mtd_none | OS: CAP MTD / no treatment | not reported | not reported | not estimable | not estimable | 1.2 | [0.323, 4.5] (95% CI) | contextual_or_secondary | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_cap_mtd_none` |
| seyedi_pp_dm_mtd | OS: GEM+CAP ping-pong dose modulation / GEM+CAP MTD | not reported | not reported | not estimable | not estimable | 0.11 | [0.024, 0.55] (95% CI) | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_pp_dm_mtd` |
| seyedi_pp_int_mtd | OS: GEM+CAP ping-pong intermittent / GEM+CAP MTD | not reported | not reported | not estimable | not estimable | 0.07 | [0.013, 0.42] (95% CI) | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_pp_int_mtd` |
| seyedi_pp_int_none | OS: GEM+CAP ping-pong intermittent / no treatment | not reported | not reported | not estimable | not estimable | 0.13 | [0.02, 0.84] (95% CI) | contextual_or_secondary | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_pp_int_none` |
| seyedi_pp_dm_none | OS: GEM+CAP ping-pong dose modulation / no treatment | not reported | not reported | not estimable | not estimable | 0.2 | [0.038, 1.09] (95% CI) | contextual_or_secondary | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_pp_dm_none` |
| seyedi_tandem_dm_mtd | OS: GEM+CAP tandem dose modulation / GEM+CAP MTD | not reported | not reported | not estimable | not estimable | 0.22 | [0.05, 0.94] (95% CI) | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_tandem_dm_mtd` |
| seyedi_tandem_int_mtd | OS: GEM+CAP tandem intermittent / GEM+CAP MTD | not reported | not reported | not estimable | not estimable | 0.34 | [0.07, 1.6] (95% CI) | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_tandem_int_mtd` |
| seyedi_tandem_dm_none | OS: GEM+CAP tandem dose modulation / no treatment | not reported | not reported | not estimable | not estimable | 0.34 | [0.07, 1.6] (95% CI) | contextual_or_secondary | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_tandem_dm_none` |
| seyedi_tandem_int_none | OS: GEM+CAP tandem intermittent / no treatment | not reported | not reported | not estimable | not estimable | 0.73 | [0.17, 3] (95% CI) | contextual_or_secondary | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_tandem_int_none` |
| seyedi_combo_mtd_none | OS: GEM+CAP MTD / no treatment | not reported | not reported | not estimable | not estimable | 1.9 | [0.48, 7.4] (95% CI) | contextual_or_secondary | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_combo_mtd_none` |
| hockings_sensitive | Sensitive survival from inoculation (weeks) | not reached/undefined | 18.75 | not estimable | not estimable | not reported | not reported | observed_latest | [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) / `hockings_sensitive_at;hockings_sensitive_st` |
| hockings_mixed | Mixed survival from inoculation (weeks) | not reached/undefined | 11.25 | not estimable | not estimable | not reported | not reported | observed_latest | [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) / `hockings_mixed_at;hockings_mixed_st` |
| hockings_resistant_st_vehicle | Resistant ST / vehicle survival (weeks) | 5.75 | 15.25 | -9.5 | 0.377049 | not reported | not reported | contextual_or_secondary | [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) / `hockings_resistant_st;hockings_resistant_vehicle` |
| pressley_adaptive_continuous | Mean overall survival (days) | 127 | 221 | -94 | 0.574661 | not reported | not reported | observed_latest | [pressley_2026](https://www.nature.com/articles/s41598-026-63987-x) / `pressley_adaptive_os;pressley_continuous_os` |
| pressley_continuous_control | Mean OS continuous / untreated (days) | 221 | 52 | 169 | 4.25 | not reported | not reported | contextual_or_secondary | [pressley_2026](https://www.nature.com/articles/s41598-026-63987-x) / `pressley_continuous_os;pressley_control_os` |
| pressley_adaptive_control | Mean OS adaptive / untreated (days) | 127 | 52 | 75 | 2.44231 | not reported | not reported | contextual_or_secondary | [pressley_2026](https://www.nature.com/articles/s41598-026-63987-x) / `pressley_adaptive_os;pressley_control_os` |

Published arm-summary intervals (not intervals for the derived ratio or difference):

| Contrast | Experimental uncertainty | Reference uncertainty |
| --- | --- | --- |
| gonzalez_cao_pfs | [5.2, 9.3] (95% CI) | [9.5, 24.1] (95% CI) |
| coin_os | [8, 24.7] (IQR) | [9.4, 26.1] (IQR) |

### Exposure contrasts

Cumulative amount, mean daily intensity and on-treatment time are different quantities. Longer survival changes the accumulation window. Doses below have no published sampling interval unless listed in the inputs.

| Contrast | Endpoint / units | Experimental | Reference | Difference | Ratio | HR | HR interval | Role | Source / input rows |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| seyedi_cap_dose_modulation_cumulative | CAP cumulative dose (mg/kg) | 533.8 | 1373.3 | -839.5 | 0.388699 | not reported | not reported | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_cap_dose_modulation_dose;seyedi_cap_mtd_dose` |
| seyedi_cap_dose_modulation_daily | CAP alternating mean daily dose (mg/kg/day) | 6.1 | 30.6 | -24.5 | 0.199346 | not reported | not reported | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_daily_cap_ping-pong_dose_modulation;seyedi_daily_cap_MTD` |
| seyedi_cap_intermittent_cumulative | CAP cumulative dose (mg/kg) | 1226.6 | 1373.3 | -146.7 | 0.893177 | not reported | not reported | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_cap_intermittent_dose;seyedi_cap_mtd_dose` |
| seyedi_cap_intermittent_daily | CAP alternating mean daily dose (mg/kg/day) | 8 | 30.6 | -22.6 | 0.261438 | not reported | not reported | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_daily_cap_ping-pong_intermittent;seyedi_daily_cap_MTD` |
| seyedi_gem_dose_modulation_cumulative | GEM cumulative dose (mg/kg) | 256.9 | 783.3 | -526.4 | 0.327971 | not reported | not reported | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_gem_dose_modulation_dose;seyedi_gem_mtd_dose` |
| seyedi_gem_dose_modulation_daily | GEM alternating mean daily dose (mg/kg/day) | 4.9 | 20.6 | -15.7 | 0.237864 | not reported | not reported | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_daily_gem_ping-pong_dose_modulation;seyedi_daily_gem_MTD` |
| seyedi_gem_intermittent_cumulative | GEM cumulative dose (mg/kg) | 333.3 | 783.3 | -450 | 0.425507 | not reported | not reported | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_gem_intermittent_dose;seyedi_gem_mtd_dose` |
| seyedi_gem_intermittent_daily | GEM alternating mean daily dose (mg/kg/day) | 6.4 | 20.6 | -14.2 | 0.31068 | not reported | not reported | observed_latest | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) / `seyedi_daily_gem_ping-pong_intermittent;seyedi_daily_gem_MTD` |

### Published fitted-model counterfactuals

These summaries came from published simulations fitted to clinical cohorts. They are not observed treatment comparisons or new independent clinical replications.

| Contrast | Endpoint / units | Experimental | Reference | Difference | Ratio | HR | HR interval | Role | Source / input rows |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| brady2020_iadt_continuous | Fitted-model mean TTP (months) | 67.11 | 57.25 | 9.86 | 1.17223 | not reported | not reported | model_counterfactual | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) / `brady20_ttp_iadt;brady20_ttp_continuous` |
| brady2020_no_induction_continuous | Fitted-model mean TTP (months) | 71.87 | 57.25 | 14.62 | 1.25537 | not reported | not reported | model_counterfactual | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) / `brady20_ttp_no_induction;brady20_ttp_continuous` |
| brady2020_no_induction_iadt | Fitted-model mean TTP (months) | 71.87 | 67.11 | 4.76 | 1.07093 | not reported | not reported | model_counterfactual | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) / `brady20_ttp_no_induction;brady20_ttp_iadt` |
| brady2020_trigger90 | Fitted-model mean TTP (months) | 66.79 | 67.11 | -0.32 | 0.995232 | not reported | not reported | model_counterfactual | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) / `brady20_ttp_trigger90;brady20_ttp_iadt` |
| brady2020_trigger50 | Fitted-model mean TTP (months) | 72.52 | 67.11 | 5.41 | 1.08061 | not reported | not reported | model_counterfactual | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) / `brady20_ttp_trigger50;brady20_ttp_iadt` |
| brady2020_trigger30 | Fitted-model mean TTP (months) | 76.6 | 67.11 | 9.49 | 1.14141 | not reported | not reported | model_counterfactual | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) / `brady20_ttp_trigger30;brady20_ttp_iadt` |

### Other published aggregates

| Input row | Measure | Value | Unit | Reported uncertainty | Statistic | Source |
| --- | --- | --- | --- | --- | --- | --- |
| z17_n_a | n | 11 | people | not reported | count | [zhang2017](https://jessica-cunningham.com/wp-content/uploads/2019/11/cunningham_et_al-2017-nature_communications.pdf) |
| z17_n_c | n | 16 | people | not reported | count | [zhang2017](https://jessica-cunningham.com/wp-content/uploads/2019/11/cunningham_et_al-2017-nature_communications.pdf) |
| z17_dose_a | cumulative_dose_fraction_vs_SOC | 0.47 | fraction | not reported | mean | [zhang2017](https://jessica-cunningham.com/wp-content/uploads/2019/11/cunningham_et_al-2017-nature_communications.pdf) |
| z17_progressed_a | radiographic_progressions_at_report | 1 | people | not reported | count | [zhang2017](https://jessica-cunningham.com/wp-content/uploads/2019/11/cunningham_et_al-2017-nature_communications.pdf) |
| z17_progressed_c | radiographic_progressions_at_report | 14 | people | not reported | count | [zhang2017](https://jessica-cunningham.com/wp-content/uploads/2019/11/cunningham_et_al-2017-nature_communications.pdf) |
| z22_n_a | n | 17 | people | not reported | count | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) |
| z22_n_c | n | 16 | people | not reported | count | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) |
| z22_dose_rate | dosing_rate_fraction_vs_SOC | 0.54 | fraction | not reported | mean | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) |
| z22_time_off | time_off_treatment_fraction | 0.46 | fraction | not reported | mean | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) |
| gc21_n_total | n | 70 | people | not reported | count | [gonzalez_cao2021](https://www.nature.com/articles/s41467-021-26572-6) |
| gc21_response_a | objective_responders | 21 | people | not reported | count | [gonzalez_cao2021](https://www.nature.com/articles/s41467-021-26572-6) |
| gc21_response_c | objective_responders | 25 | people | not reported | count | [gonzalez_cao2021](https://www.nature.com/articles/s41467-021-26572-6) |
| int24_n_a | n | 39 | people | not reported | count | [interim2024](https://eprints.whiterose.ac.uk/id/eprint/206919/7/EJC-D-23-02078_R1.pdf) |
| int24_n_c | n | 40 | people | not reported | count | [interim2024](https://eprints.whiterose.ac.uk/id/eprint/206919/7/EJC-D-23-02078_R1.pdf) |
| int24_response_n_a | response_evaluable_n | 37 | people | not reported | count | [interim2024](https://eprints.whiterose.ac.uk/id/eprint/206919/7/EJC-D-23-02078_R1.pdf) |
| int24_response_n_c | response_evaluable_n | 35 | people | not reported | count | [interim2024](https://eprints.whiterose.ac.uk/id/eprint/206919/7/EJC-D-23-02078_R1.pdf) |
| coin11_n_a | n | 815 | people | not reported | count | [coin2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3159416/) |
| coin11_n_c | n | 815 | people | not reported | count | [coin2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3159416/) |
| adt22_n | n | 16 | people | not reported | count | [zhang_adt2022](https://www.mdpi.com/2072-6694/14/21/5225) |
| adt22_nonprogression12 | no_imaging_progression_at_12_months | 16 | people | not reported | count | [zhang_adt2022](https://www.mdpi.com/2072-6694/14/21/5225) |
| adt22_followup | followup | 26 | months | not reported | median | [zhang_adt2022](https://www.mdpi.com/2072-6694/14/21/5225) |
| adt22_rttp | radiographic_TTP | not reached/undefined | months | not reported | median | [zhang_adt2022](https://www.mdpi.com/2072-6694/14/21/5225) |
| adt22_psattp | PSA_TTP | not reached/undefined | months | not reported | median | [zhang_adt2022](https://www.mdpi.com/2072-6694/14/21/5225) |
| seyedi_cap_dm_mtd_p | overall_survival_hazard_ratio_p_value | 0.065 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_cap_int_mtd_p | overall_survival_hazard_ratio_p_value | 0.68 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_cap_dm_none_p | overall_survival_hazard_ratio_p_value | 0.1 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_cap_int_none_p | overall_survival_hazard_ratio_p_value | 0.95 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_cap_dm_int_p | overall_survival_hazard_ratio_p_value | 0.09 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_cap_mtd_none_p | overall_survival_hazard_ratio_p_value | 0.78 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_pp_dm_mtd_p | overall_survival_hazard_ratio_p_value | 0.007 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_pp_int_mtd_p | overall_survival_hazard_ratio_p_value | 0.003 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_pp_int_none_p | overall_survival_hazard_ratio_p_value | 0.032 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_pp_dm_none_p | overall_survival_hazard_ratio_p_value | 0.06 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_tandem_dm_mtd_p | overall_survival_hazard_ratio_p_value | 0.04 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_tandem_int_mtd_p | overall_survival_hazard_ratio_p_value | 0.17 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_tandem_dm_none_p | overall_survival_hazard_ratio_p_value | 0.17 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_tandem_int_none_p | overall_survival_hazard_ratio_p_value | 0.66 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_combo_mtd_none_p | overall_survival_hazard_ratio_p_value | 0.36 | probability | not reported | reported_p | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_gem_survival_missing | overall_survival_comparison | not estimable |  | not reported | not_quantified | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| hockings_sensitive_p | survival_logrank_p | 0.01 | probability | not reported | reported_p | [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) |
| hockings_mixed_p | survival_logrank_p | 0.06 | probability | not reported | reported_p | [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) |
| hockings_resistant_st_p | survival_logrank_p | 0.09 | probability | not reported | reported_p | [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) |
| hockings_resistant_at_p | survival_logrank_p | 0.04 | probability | not reported | reported_p | [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) |
| pressley_omnibus_p | survival_kruskal_wallis_p | <0.002 | probability | not reported | reported_upper_bound | [pressley_2026](https://www.nature.com/articles/s41598-026-63987-x) |
| wang_reach_logrank_1 | time_to_volume_endpoint_logrank_p | 0.1671 | probability | not reported | reported_p | [wang_2021_reachability](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) |
| wang_reach_logrank_2 | time_to_volume_endpoint_logrank_p | 0.0035 | probability | not reported | reported_p | [wang_2021_reachability](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) |
| wang_reach_logrank_3 | time_to_volume_endpoint_logrank_p | 0.0035 | probability | not reported | reported_p | [wang_2021_reachability](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) |
| wang_reach_logrank_4 | time_to_volume_endpoint_logrank_p | 0.0035 | probability | not reported | reported_p | [wang_2021_reachability](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) |
| wang_reach_logrank_5 | time_to_volume_endpoint_logrank_p | 0.0035 | probability | not reported | reported_p | [wang_2021_reachability](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) |
| wang_reach_logrank_6 | time_to_volume_endpoint_logrank_p | 0.0035 | probability | not reported | reported_p | [wang_2021_reachability](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) |
| wang_reach_logrank_7 | time_to_volume_endpoint_logrank_p | 0.01 | probability | not reported | reported_p | [wang_2021_reachability](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) |
| wang_reach_logrank_8 | time_to_volume_endpoint_logrank_p | 0.0035 | probability | not reported | reported_p | [wang_2021_reachability](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) |
| wang_reach_logrank_9 | time_to_volume_endpoint_logrank_p | 0.0035 | probability | not reported | reported_p | [wang_2021_reachability](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) |
| wang_reach_logrank_10 | time_to_volume_endpoint_logrank_p | 0.0065 | probability | not reported | reported_p | [wang_2021_reachability](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) |
| bacevic_monolayer_fitness | relative_fitness | not estimable | ratio | [0.9, 0.95] (descriptive_range_not_confidence_interval) | reported_range | [bacevic_2017](https://www.nature.com/articles/s41467-017-01516-1) |
| hockings_low_resource_growth_difference | sensitive_minus_resistant_growth | 0.08 | doublings/day | not reported | fitted_group_growth_difference | [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) |
| pressley_growth_fixed | net_tumor_growth_rate | 0.076 | per_day | not reported | common_fixed_fitted_estimate | [pressley_2026](https://www.nature.com/articles/s41598-026-63987-x) |
| strobl22_n_fitted | n | 65 | people | not reported | count | [strobl2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9053239/) |
| strobl22_n_excluded_means | n_excluded_from_mean | 8 | people | not reported | count | [strobl2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9053239/) |
| strobl22_cost | resistance_cost_c_R | 0.78 | fraction | not reported | fitted_group_mean | [strobl2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9053239/) |
| strobl22_turnover | normalized_turnover_d_T | 0.14 | fraction | not reported | fitted_group_mean | [strobl2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9053239/) |
| strobl22_density | initial_density_n0 | 0.59 | fraction | not reported | fitted_group_mean | [strobl2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9053239/) |
| strobl22_resistance | initial_resistance_fraction_f_R | 0.04 | fraction | not reported | fitted_group_mean | [strobl2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9053239/) |
| strobl22_r_s_fixed | sensitive_proliferation_r_S | 0.027 | day^-1 | not reported | assumed_model_constant | [strobl2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9053239/) |
| brady20_n | n | 70 | people | not reported | count | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) |
| brady20_n_training | n | 35 | people | not reported | count | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) |
| brady20_n_testing | n | 35 | people | not reported | count | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) |
| brady20_ps | stem_self_renewal_p_s | 0.0278 | dimensionless | [0.0265, 0.0611] (95% CI) | fitted_group_median | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) |
| brady20_alpha | ADT_cytotoxicity_alpha | 0.036 | day^-1 | [0.0258, 0.1379] (95% CI) | fitted_group_median | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) |
| brady20_rho | PSA_production_rho | 0.000187 | microgram_per_litre_per_day | not reported | fitted_common_parameter | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) |
| brady20_phi | PSA_decay_phi | 0.0856 | day^-1 | not reported | fitted_common_parameter | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) |
| brady20_dose_trigger90 | time_adjusted_dose_fraction_vs_IADT | ~0.6 | fraction | not reported | source_approximate_model_projection | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) |
| brady20_dose_trigger50 | time_adjusted_dose_fraction_vs_IADT | ~0.25 | fraction | not reported | source_approximate_model_projection | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) |
| brady20_dose_trigger30 | time_adjusted_dose_fraction_vs_IADT | ~0.22 | fraction | not reported | source_approximate_model_projection | [brady2020](https://www.nature.com/articles/s41467-020-15424-4) |
| brady22_n | n | 16 | people | not reported | count | [brady2022_rbat](https://pmc.ncbi.nlm.nih.gov/articles/PMC9657943/) |
| brady22_n_ttp_improved | n_simulations_with_longer_TTP | 11 | people_simulated | not reported | count | [brady2022_rbat](https://pmc.ncbi.nlm.nih.gov/articles/PMC9657943/) |
| brady22_n_dose_reduced | n_simulations_with_reduced_cumulative_dose | 15 | people_simulated | not reported | count | [brady2022_rbat](https://pmc.ncbi.nlm.nih.gov/articles/PMC9657943/) |
| brady22_n_both_controlled | n_controlled_through_horizon | 4 | people_simulated | not reported | count | [brady2022_rbat](https://pmc.ncbi.nlm.nih.gov/articles/PMC9657943/) |
| brady22_horizon | simulation_horizon | 10 | years | not reported | assumed_design | [brady2022_rbat](https://pmc.ncbi.nlm.nih.gov/articles/PMC9657943/) |
| brady22_n_ttp_worse | n_simulations_with_shorter_TTP | 1 | people_simulated | not reported | count | [brady2022_rbat](https://pmc.ncbi.nlm.nih.gov/articles/PMC9657943/) |
| brady22_max_dose_reduction | maximum_cumulative_dose_reduction | 0.4 | fraction | not reported | simulated_maximum | [brady2022_rbat](https://pmc.ncbi.nlm.nih.gov/articles/PMC9657943/) |
| kim21_n | n | 8 | people | not reported | count | [kim2021](https://www.mdpi.com/2072-6694/13/4/823) |
| kim21_fits_per_subject | fits_per_subject | 50 | fits | not reported | simulation_design | [kim2021](https://www.mdpi.com/2072-6694/13/4/823) |
| kim21_min_ttp_gain | TTP_gain_range_min | 6 | months | not reported | model_sensitivity_range_endpoint | [kim2021](https://www.mdpi.com/2072-6694/13/4/823) |
| kim21_max_ttp_gain | TTP_gain_range_max | 25 | months | not reported | model_sensitivity_range_endpoint | [kim2021](https://www.mdpi.com/2072-6694/13/4/823) |
| kim21_min_time_on | time_on_treatment_fraction_min | 0.06 | fraction | not reported | model_sensitivity_range_endpoint | [kim2021](https://www.mdpi.com/2072-6694/13/4/823) |
| kim21_max_time_on | time_on_treatment_fraction_max | 0.74 | fraction | not reported | model_sensitivity_range_endpoint | [kim2021](https://www.mdpi.com/2072-6694/13/4/823) |
| kim21_corr_r | correlation_TTP_fitted_rR | -0.77 | correlation | not reported | fitted_parameter_association | [kim2021](https://www.mdpi.com/2072-6694/13/4/823) |
| kim21_corr_c | correlation_TTP_fitted_C | -0.64 | correlation | not reported | fitted_parameter_association | [kim2021](https://www.mdpi.com/2072-6694/13/4/823) |
| kim21_corr_r_switch | correlation_TTP_fitted_rR_switch_model | -0.67 | correlation | not reported | fitted_parameter_association | [kim2021](https://www.mdpi.com/2072-6694/13/4/823) |
| z22_r_s | sensitive_net_growth_r_S | 0.0156 | day^-1 | not reported | PSA_slope_group_mean | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) |
| z22_r_r | resistant_net_growth_r_R | 0.0091 | day^-1 | not reported | PSA_slope_group_mean | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) |
| z22_alpha_rs | competition_alpha_RS | 6 | dimensionless | not reported | optimized_shared_parameter | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) |
| z22_beta_sc | growth_scale_beta_SC | 8 | dimensionless | not reported | optimized_shared_parameter | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) |
| z22_alpha_sr | competition_alpha_SR | 1 | dimensionless | not reported | assumed_model_constant | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) |
| z22_K | carrying_capacity_scale_K | 10000 | arbitrary_scalable_burden_units | not reported | assumed_model_constant | [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) |
| wang21_fit_r_a | resistant_A549_intrinsic_growth_r_a | 0.184 | day^-1 | not reported | fitted_population_typical_value | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| wang21_fit_r_a_iiv | resistant_A549_intrinsic_growth_r_a_inter_individual_variation | 0.183 | source_log_normal_variation_scale | not reported | source_inter_individual_variation | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| wang21_fit_r_h | sensitive_H1975_intrinsic_growth_r_h | 0.187 | day^-1 | not reported | fitted_population_typical_value | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| wang21_fit_r_h_iiv | sensitive_H1975_intrinsic_growth_r_h_inter_individual_variation | 0.254 | source_log_normal_variation_scale | not reported | source_inter_individual_variation | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| wang21_fit_K_a | resistant_A549_carrying_capacity_K_a | 1350 | mm3 | not reported | fitted_population_typical_value | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| wang21_fit_K_a_iiv | resistant_A549_carrying_capacity_K_a_inter_individual_variation | 0.316 | source_log_normal_variation_scale | not reported | source_inter_individual_variation | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| wang21_fit_K_h | sensitive_H1975_carrying_capacity_K_h | 1410 | mm3 | not reported | fitted_population_typical_value | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| wang21_fit_K_h_iiv | sensitive_H1975_carrying_capacity_K_h_inter_individual_variation | 0.352 | source_log_normal_variation_scale | not reported | source_inter_individual_variation | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| wang21_fit_a | competition_A549_on_H1975_a | 0.109 | dimensionless | not reported | fitted_population_typical_value | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| wang21_fit_a_iiv | competition_A549_on_H1975_a_inter_individual_variation | 2.05 | source_log_normal_variation_scale | not reported | source_inter_individual_variation | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| wang21_fit_h | competition_H1975_on_A549_h | 0.554 | dimensionless | not reported | fitted_population_typical_value | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| wang21_fit_h_iiv | competition_H1975_on_A549_h_inter_individual_variation | 1.56 | source_log_normal_variation_scale | not reported | source_inter_individual_variation | [wang_2021_reachability](https://www.mdpi.com/article/10.3390/cancers13215262/s1) |
| seyedi_n_treatment | initial_group_size | 6 | mice | not reported | count | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_n_vehicle | initial_group_size | 4 | mice | not reported | count | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| seyedi_n_total | total_randomized | 70 | mice | not reported | count | [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) |
| hockings_n_st | initial_group_size | 2 | mice | not reported | count | [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) |
| hockings_n_vehicle | initial_group_size | 2 | mice | not reported | count | [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) |
| hockings_n_adaptive_range | initial_group_size | not estimable | mice | [3, 5] (descriptive range; not a CI) | reported_range | [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) |
| pressley_n_arm | initial_group_size | 5 | mice | not reported | count | [pressley_2026](https://www.nature.com/articles/s41598-026-63987-x) |

### Parameter identifiability by report

`conditional` denotes a point transformation inside the source's fitted model, with no transferable uncertainty interval. It does not satisfy PR #61's eligibility gates. Full reasons and source locators are in `data/mapping.csv`.

| Report | d | c | J | A | lambda |
| --- | --- | --- | --- | --- | --- |
| [zhang2017](https://jessica-cunningham.com/wp-content/uploads/2019/11/cunningham_et_al-2017-nature_communications.pdf) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [zhang2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) | conditional | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [gonzalez_cao2021](https://www.nature.com/articles/s41467-021-26572-6) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [interim2024](https://eprints.whiterose.ac.uk/id/eprint/206919/7/EJC-D-23-02078_R1.pdf) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [coin2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3159416/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [zhang_adt2022](https://www.mdpi.com/2072-6694/14/21/5225) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [algazi2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC8063889/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [improve2025](https://ascopubs.org/doi/10.1200/JCO.24.00979) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [improve2026](https://www.annalsofoncology.org/article/S0923-7534%2826%2900257-7/fulltext) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [bfr14_2024](https://pubmed.ncbi.nlm.nih.gov/39127063/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [crook2012](https://doi.org/10.1056/NEJMoa1201546) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [hussain2013](https://www.nejm.org/doi/full/10.1056/NEJMoa1212299) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [iceland2016](https://pubmed.ncbi.nlm.nih.gov/26520703/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [optimox2_2009](https://doi.org/10.1200/JCO.2009.23.4344) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [anzadapt](https://www.anzctr.org.au/Trial/Registration/TrialReview.aspx?id=21010&isClinicalTrial=True) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [actov](https://pmc.ncbi.nlm.nih.gov/articles/PMC11667365/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [bcc_vismodegib](https://www.cancer.gov/research/participate/clinical-trials-search/v?id=NCI-2022-10569&r=1) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [adaptive_capecitabine](https://clinicaltrials.gov/study/NCT06525766) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [adaptive_rcc](https://www.cancer.gov/research/participate/clinical-trials-search/v?id=NCI-2025-01856) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [adaptive_adt_docetaxel](https://www.cancer.gov/research/participate/clinical-trials-search/v?id=NCI-2024-10406) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [mistry2021](https://www.nature.com/articles/s41467-020-20174-4) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [seyedi_2024](https://www.mdpi.com/2072-6694/16/2/257) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [hockings_2025](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [pressley_2026](https://www.nature.com/articles/s41598-026-63987-x) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [wang_2021_reachability](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) | conditional | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [bacevic_2017](https://www.nature.com/articles/s41467-017-01516-1) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [gatenby_2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC3728826/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [enriquez_navas_2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4962860/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [smalley_2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6838387/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [strobl_2024_parpi](https://pmc.ncbi.nlm.nih.gov/articles/PMC11190943/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [silva_2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3525750/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [gallaher_2018](https://aacrjournals.org/cancerres/article/78/8/2127/633375/Spatial-Heterogeneity-and-Evolutionary-Dynamics) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [das_thakur_2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3930354/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [maltas_2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10918730/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [wang_2021_fixed_periodic](https://pmc.ncbi.nlm.nih.gov/articles/PMC8643150/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [kam_2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4262548/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [kam_2014_review](https://pmc.ncbi.nlm.nih.gov/articles/PMC4258899/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [antifragile_2026](https://pmc.ncbi.nlm.nih.gov/articles/PMC13283037/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [plasticity_2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12515544/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [microbial_2026](https://www.nature.com/articles/s41559-026-03178-z) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [gevertz_2026](https://www.nature.com/articles/s41540-025-00635-6) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [enriquez_2014_ismrm](https://cds.ismrm.org/protected/14MProceedings/PDFfiles/1130.pdf) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [enriquez_2015_aacr](https://aacrjournals.org/mct/article/14/12_Supplement_1/B11/285815/Abstract-B11-Adaptive-therapy-modifies-tumor) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [colorectal_2026_abstract](https://aacrjournals.org/cancerres/article/86/7_Supplement/702/776389/Abstract-702-Harnessing-clonal-dynamics-to) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [qi_2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7874629/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [seyedi_2023_preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC10541126/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [hockings_2023_preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC10401956/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [antifragile_2025_preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458151/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [chattopadhyay_2024_abstract](https://doi.org/10.1158/1538-7445.CANEvol23-PR003) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [seyedi_2022_model](https://pmc.ncbi.nlm.nih.gov/articles/PMC9179496/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [acar_2020_steering](https://pmc.ncbi.nlm.nih.gov/articles/PMC7174377/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [strobl2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9053239/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [brady2020](https://www.nature.com/articles/s41467-020-15424-4) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [brady2022_rbat](https://pmc.ncbi.nlm.nih.gov/articles/PMC9657943/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [kim2021](https://www.mdpi.com/2072-6694/13/4/823) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [strobl2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8455086/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [west2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7307613/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [brady2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8322456/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [gallagher2026](https://jamanetwork.com/journals/jamaoncology/fullarticle/2852635) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [pasetto2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8604892/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [west2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6665681/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [hansen2020_modifying](https://pmc.ncbi.nlm.nih.gov/articles/PMC7761372/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [salvioli2024](https://doi.org/10.1371/journal.pone.0308173) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [model_deep_learning2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11148552/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [deep_rl2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11174533/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [uncertainty2026](https://arxiv.org/abs/2608.18387) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [gallaher2025_psa_testosterone](https://pmc.ncbi.nlm.nih.gov/articles/PMC12248141/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [gallaher2022_sum](https://jeffreywest.github.io/pdf/gallaher2022sum.pdf) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |
| [gallagher2025_timing](https://pmc.ncbi.nlm.nih.gov/articles/PMC12417256/) | unidentifiable | unidentifiable | unidentifiable | unidentifiable | unidentifiable |

Conditional transformations (not tumour threshold estimates):

| Report | Parameter | Value (days) | Input rows | Uncertainty |
| --- | --- | --- | --- | --- |
| zhang2022 | d | 13.7363 | z22_beta_sc;z22_r_r | not identifiable from reported uncertainty |
| wang_2021_reachability | d | 5.43478 | wang21_fit_r_a | not identifiable from reported uncertainty |

### Prediction verdicts and declared sensitivities

Eligibility counts are report counts for mechanistic rows and distinct comparison groups for P6. No mechanistic slope, threshold, establishment or crossing test passed its data gate; **0 formal mechanistic hypothesis tests and 0 pooled effects** were calculated. Holm correction therefore has no p-values to adjust. P6 uses descriptive published intervals only.

| Prediction | Verdict | Eligible | Evidence rows | Interpretation |
| --- | --- | --- | --- | --- |
| P1 | untestable | 0 | data/gates.csv and data/mapping.csv (all study-specific source locators) | No independently measured controlled mutation-supply series with the theorem's mean endpoint. |
| P2 | untestable | 0 | data/gates.csv and data/mapping.csv (all study-specific source locators) | No fixed-supply measured competition contrast; increasing competition is an auxiliary hypothesis, not a universal consequence of PR #61. |
| P2a | untestable | 0 | data/gates.csv and data/mapping.csv (all study-specific source locators) | No measured zero-competition comparison satisfying the theorem's assumptions; absence of detected density dependence is not alpha=0. |
| P3 | untestable | 0 | data/gates.csv and data/mapping.csv (all study-specific source locators) | No report identifies d,c,J,A,lambda jointly; clinical medians and HRs do not identify the uncensored mean coefficients. |
| P4 | untestable | 0 | data/gates.csv and data/mapping.csv (all study-specific source locators) | No aggregate count of new established lineages with independent lambda and A; treatment failures are not founder counts. |
| P5 | untestable | 0 | data/gates.csv and data/mapping.csv (all study-specific source locators) | No otherwise comparable independently measured mutation-supply series; different cancers are not supply levels. |
| P6a | consistent | 8 | zhang2022_ttp;zhang2022_os;seyedi_cap_dm_mtd;seyedi_cap_int_mtd;seyedi_pp_dm_mtd;seyedi_pp_int_mtd;seyedi_tandem_dm_mtd;seyedi_tandem_int_mtd | Some observed comparisons favour reduced/adaptive exposure. This is compatibility with possible benefit, not a test of the mutation threshold; no pooled efficacy claim. |
| P6b | falsified | 6 | coin_os | BROAD outcome reading of the ambiguous universal extension only: COIN OS has an adverse published 80% CI. Thus the original any-reported-interval rule is met. This is a descriptive contradiction at that level, NOT a 95% rejection, proof of harm from noninferiority failure, tumour-control falsification, or PR #61 falsification. |
| P6b-control | untestable | 2 | gonzalez_cao_pfs;interim_pfs | Explicit POST-EXTRACTION sensitivity restricted to tumour-control endpoints: INTERIM's PFS CI spans 1; Gonzalez-Cao supplies arm-median CIs, not a contrast CI. Original endpoint ambiguity prevents a unique fully pre-specified P6b verdict. Adverse point estimates remain visible. |
| P6b-control95 | untestable | 2 | gonzalez_cao_pfs;interim_pfs | Explicit POST-EXTRACTION sensitivity restricted to tumour-control endpoints and published 95% intervals: INTERIM's PFS CI spans 1; Gonzalez-Cao supplies arm-median CIs, not a contrast CI. Original endpoint ambiguity prevents a unique fully pre-specified P6b verdict. Adverse point estimates remain visible. |

### Search coverage and exclusions

Every screened source below is retained regardless of outcome or extraction permission. Rights/access exclusions and missing posted results are not negative efficacy outcomes. Original dispositions, full reuse basis and exact source locators remain in the CSVs.

| Source | Design / setting | Disposition | Limitation / reason |
| --- | --- | --- | --- |
| [zhang2017: Integrating evolutionary dynamics into treatment of metastatic castrate-resistant prostate cancer](https://jessica-cunningham.com/wp-content/uploads/2019/11/cunningham_et_al-2017-nature_communications.pdf) | single_arm_with_nonrandomized_contemporaneous_control; human_mCRPC | included | Preliminary cohort update; use2022 for same-cohort primary comparisons; lower-bound median; control median14 months not2022 value14.3 |
| [zhang2022: Evolution-based mathematical models significantly prolong response to abiraterone in metastatic castrate-resistant prostate cancer and identify strategies to further improve outcomes](https://pmc.ncbi.nlm.nih.gov/articles/PMC9239688/) | single_arm_with_nonrandomized_contemporaneous_control; human_mCRPC | included | Latest full comparative report located; selected initial PSA responders; self-selected SOC chart-review comparator; postprogression and on-treatment endpoint rules differ from PR61 |
| [gonzalez_cao2021: Intermittent BRAF inhibition in advanced BRAF mutated melanoma results of a phase II randomized trial](https://www.nature.com/articles/s41467-021-26572-6) | randomized_phase2_fixed_intermittent_vs_continuous; human_melanoma | included | Experimental intermittent vemurafenib/cobimetinib after standard induction; comparator standard cobimetinib itself three weeks on one week off; medians not uncensored means |
| [interim2024: A randomised phase 2 study of intermittent versus continuous dosing of dabrafenib plus trametinib in patients with BRAFV600 mutant advanced melanoma (INTERIM)](https://eprints.whiterose.ac.uk/id/eprint/206919/7/EJC-D-23-02078_R1.pdf) | randomized_phase2_fixed_intermittent_vs_continuous; human_melanoma | included | Randomized79; no induction selection; intermittent dabrafenibd1-21/trametinibd1-14 in28-day cycles; published HR CIs include null Registry identifiers in author manuscript: ISRCTN18183156; EudraCT2016-005228-27. Earlier scratch NCT03352947 was unverified and corrected. |
| [coin2011: Intermittent versus continuous oxaliplatin and fluoropyrimidine combination chemotherapy for first-line treatment of advanced colorectal cancer: results of the randomised phase3 MRC COIN trial](https://pmc.ncbi.nlm.nih.gov/articles/PMC3159416/) | randomized_phase3_noninferiority_fixed_intermittent_vs_continuous; human_colorectal | included | Use ITT armsA/C only; third cetuximab arm irrelevant. Survival interval after medians is IQR not confidence interval; HR uses80%CI. Treatment-holiday design endpoint differs from AT |
| [zhang_adt2022: A Phase1b Adaptive Androgen Deprivation Therapy Trial in Metastatic Castration Sensitive Prostate Cancer](https://www.mdpi.com/2072-6694/14/21/5225) | single_arm_feasibility; human_mCSPC | included_single_arm | No continuous comparator. Methods correctly identifyNCT03511196; Introduction erroneously repeatsNCT02415621. Both median progression endpoints not reached; do not set equal to26-month median followup |
| [algazi2020: Continuous versus intermittent BRAF and MEK inhibition in patients with BRAF-mutated melanoma: a randomized phase2 trial (S1320)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8063889/) | randomized_phase2_fixed_intermittent_vs_continuous; human_melanoma | metadata_rights_or_access | Known negative comparator retained visibly; published HR interval is80% and prespecified two-sidedalpha.20; do not label95% or conventionalP<.05. Registry results access incomplete |
| [improve2025: Intermittent or Continuous Panitumumab Plus Fluorouracil Leucovorin and Irinotecan for First-Line Treatment of RAS and BRAF Wild-Type Metastatic Colorectal Cancer: The IMPROVE Trial](https://ascopubs.org/doi/10.1200/JCO.24.00979) | randomized_phase2_noncomparative_intermittent_vs_continuous; human_colorectal | metadata_rights_or_access | Quantitative comparison exists but formal design noncomparative; PFS-on-treatment excludes treatment holidays and must not equate to ordinary PFS or PR61 progression |
| [improve2026: IMPROVE survival update and ctDNA analysis abstract72P](https://www.annalsofoncology.org/article/S0923-7534%2826%2900257-7/fulltext) | conference_update_same_cohort; human_colorectal | metadata_rights_or_access | 2026 update located; same cohort asIMPROVE2025 not independent study; ctDNA emergent resistance mutations are not independently measured per-divisionmu |
| [bfr14_2024: Discontinuation versus continuation of imatinib in patients with advanced gastrointestinal stromal tumours (BFR14): exploratory long-term follow-up](https://pubmed.ncbi.nlm.nih.gov/39127063/) | randomized_phase3_interruption_vs_continuation_followup; human_GIST | metadata_rights_or_access | Strong negative interruption context; prespecified interruption after sustained disease control differs from sensitive-reservoir AT; multiple timed randomization cohorts must not assume independence |
| [crook2012: Intermittent androgen suppression for rising PSA level after radiotherapy](https://doi.org/10.1056/NEJMoa1201546) | randomized_phase3_noninferiority_intermittent_vs_continuous; human_prostate_biochemical_recurrence | metadata_rights_or_access | Landmark broader intermittent comparison; no original PR61 mechanistic rate/lineage measurements |
| [hussain2013: Intermittent versus Continuous Androgen Deprivation in Prostate Cancer](https://www.nejm.org/doi/full/10.1056/NEJMoa1212299) | randomized_phase3_noninferiority_intermittent_vs_continuous; human_metastatic_hormone_sensitive_prostate | metadata_rights_or_access | SWOG9346; noninferiority not established; seven-month induction selection differs from PR61 sensitive maintenance. Identifier independently verified by primary secondary report PMC2684851. |
| [iceland2016: Intermittent Versus Continuous Androgen Deprivation Therapy in Patients with Relapsing or Locally Advanced Prostate Cancer (ICELAND)](https://pubmed.ncbi.nlm.nih.gov/26520703/) | randomized_phase3b_intermittent_vs_continuous; human_nonmetastatic_prostate | metadata_rights_or_access | Additional null intermittent comparison; source reports similar efficacy; group injection aggregate present but no rights to derivedCSV established |
| [optimox2_2009: Can chemotherapy be discontinued in unresectable metastatic colorectal cancer? The GERCOR OPTIMOX2 Study](https://doi.org/10.1200/JCO.2009.23.4344) | randomized_phase2_maintenance_vs_chemotherapy_free_interval; human_colorectal | metadata_rights_or_access | Broader chemotherapy-free interval comparator screened; maintenance arm is not common maximum-dose arm; do not copy secondary-review or conference values in place of final primary report |
| [anzadapt: ANZadapt phaseII randomized adaptive versus continuous abiraterone or enzalutamide](https://www.anzctr.org.au/Trial/Registration/TrialReview.aspx?id=21010&isClinicalTrial=True) | registered_randomized_phase2; human_mCRPC | no_posted_results_observed | ANZCTR Results section explicitly says results not provided and no documents uploaded. ANZUP sponsor pagehttps://anzup.org.au/clinical-trial/anzadapt/ says In Follow-up whereas older registry text still Recruiting. Do not invent results |
| [actov: Adaptive ChemoTherapy for Ovarian cancer (ACTOv)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11667365/) | registered_randomized_phase2; human_ovarian | results_not_located | 2024 BMJOpen protocol and sponsor2026 open recruitment listing found; no primary aggregate results located. Registry browser shell/API unavailable, so absence is not a confirmed global registry negative |
| [bcc_vismodegib: Adaptive Therapy of Vismodegib in Advanced Basal Cell Carcinoma](https://www.cancer.gov/research/participate/clinical-trials-search/v?id=NCI-2022-10569&r=1) | registered_randomized_early_phase1; human_BCC | results_not_located | AACR2025 abstract6242 describes interim modeling without quantitative group outcome comparison. NCI active; direct registry unreadable in tool, secondary registry says No Results but not treated as primary proof |
| [adaptive_capecitabine: Adaptive Therapy With Capecitabine for Treatment of Metastatic ER Positive HER2 Negative Breast Cancer](https://clinicaltrials.gov/study/NCT06525766) | registered_phase2; human_breast | results_not_located | No primary aggregate results located in query; secondary registry snapshotRecruiting but primary results tab not inspectable. Do not assert definitive no-results |
| [adaptive_rcc: Adaptive Pembrolizumab and Axitinib for Treatment of Advanced Clear Cell Renal Cell Carcinoma](https://www.cancer.gov/research/participate/clinical-trials-search/v?id=NCI-2025-01856) | registered_phase2; human_RCC | results_not_located | Primary NCI lists closed to accrual and intervention; secondary registry copies Suspended. No primary quantitative results located; registry API inaccessible |
| [adaptive_adt_docetaxel: Adaptive Androgen Deprivation Therapy and Docetaxel for Treatment of Metastatic Castration Sensitive Prostate Cancer](https://www.cancer.gov/research/participate/clinical-trials-search/v?id=NCI-2024-10406) | registered_phase2a; human_mCSPC | results_not_located | Current NCI trial description located; no published aggregate result found in targeted search |
| [mistry2021: On the reporting and analysis of a cancer evolutionary adaptive dosing trial](https://www.nature.com/articles/s41467-020-20174-4) | methodological_criticism; human_mCRPC_commentary | excluded_scope | Retain published statistical critique of2017 trial and responsehttps://www.nature.com/articles/s41467-020-20175-3; not independent treatment comparison Primary indexed text now retrieved: critique concerns unreported prognostic covariates and responder selection versus historical PhaseIII; author response admits imperfect nonrandomized controls but says contemporaneous controls also met responder eligibility. Do not apply historical responder mismatch to contemporaneous cohort. |
| [seyedi_2024: Testing Adaptive Therapy Protocols Using Gemcitabine and Capecitabine in a Preclinical Model of Endocrine-Resistant Breast Cancer](https://www.mdpi.com/2072-6694/16/2/257) | Random allocation reported; MCF7 endocrine-resistant breast xenograft; capecitabine and gemcitabine single and combination protocols; shared MTD and untreated controls; mouse_xenograft | included | Include all directly printed active comparator HRs and exposure means, retaining null single-GEM result. No pooling of shared arms. All HR, confidence limits, p-values and dose values independently crosschecked against official EuropePMC fullTextXML PMC10813385 on access date. Authors Discussion4.1 warn combined GEM/CAP nominal MTD may exceed tolerable dosing; toxicity may explain much of the multidrug advantage. |
| [hockings_2025: Adaptive Therapy Exploits Fitness Deficits in Chemotherapy-Resistant Ovarian Cancer to Achieve Long-Term Tumor Control](https://aacrjournals.org/cancerres/article/85/18/3503/764448/Adaptive-Therapy-Exploits-Fitness-Deficits-in) | Random allocation reported; OVCAR4 sensitive; Ov4Carbo-Luc resistant; and mixed 80:20 sensitive:resistant cohorts; carboplatin adaptive vs standard treatment and vehicle; mouse_xenograft | included | Censored adaptive medians remain undefined; no ratio to an invented endpoint. Distinct resistant-only control comparison kept separate. Authors explicitly cannot distinguish evolutionary benefit from higher cumulative AT dose in this experiment. |
| [pressley_2026: Density dependence and evolvability limit adaptive therapy in non-small cell lung cancer mouse model](https://www.nature.com/articles/s41598-026-63987-x) | H3122 NSCLC xenograft in NSG mice; alectinib continuous versus response-triggered adaptive versus untreated; mouse_xenograft | included | Include adverse adaptive direction; no mean variance or confidence interval printed. Omnibus p is not an arm-specific interval. |
| [wang_2021_reachability: Optimizing Adaptive Therapy Based on the Reachability to Tumor Resistant Subpopulation](https://pmc.ncbi.nlm.nih.gov/articles/PMC8582524/) | Random allocation reported; H1975/A549 mixed NSCLC, osimertinib; classic AT1, reachability AT2, MTD, late-start positive control and untreated; mouse_xenograft | included | Table1 gives pairwise corrected logrank p values and narrative direction; no printed numeric median or HR. Curves not digitized. Companion Wang2021 CMAR same authors/cells but cohort overlap unresolved. |
| [bacevic_2017: Spatial competition constrains resistance to targeted cancer therapy](https://www.nature.com/articles/s41467-017-01516-1) | HCT116 WT and NU6102-resistant R50; adaptive and high-dose monolayer, low versus high-dose spheroids; in_vitro | included | Experimental efficacy curves lack exact endpoint aggregates in text; numerical median increases in Figure5 are simulations, excluded. Preserve monolayer no-advantage vs spheroid advantage narrative. |
| [gatenby_2009: Adaptive Therapy](https://pmc.ncbi.nlm.nih.gov/articles/PMC3728826/) | OVCAR3 ovarian xenograft carboplatin adaptive vs fixed standard; mouse_xenograft | metadata_rights_or_access | Do not mistake simulation durations for observed xenograft outcomes; main experimental comparisons graph-based. |
| [enriquez_navas_2016: Exploiting evolutionary principles to prolong tumor control in preclinical models of breast cancer](https://pmc.ncbi.nlm.nih.gov/articles/PMC4962860/) | MDA-MB-231/luc and MCF7 orthotopic breast; paclitaxel AT1 modulation AT2 skipping vs standard; mouse_xenograft | metadata_rights_or_access | Both favorable dose modulation and weaker skipping retained in metadata; continuing AT versus finite standard course mismatches PR61. |
| [smalley_2019: Leveraging transcriptional dynamics to improve BRAF inhibitor responses in melanoma](https://pmc.ncbi.nlm.nih.gov/articles/PMC6838387/) | WM164 melanoma; model-informed adaptive versus continuous and fixed intermittent; mouse_xenograft | metadata_rights_or_access | Phenotypic/transcriptional dynamics provide alternative mechanism; no data copied due outcome-independent rights gate. |
| [strobl_2024_parpi: To modulate or to skip: De-escalating PARP inhibitor maintenance therapy in ovarian cancer using adaptive therapy](https://pmc.ncbi.nlm.nih.gov/articles/PMC11190943/) | Ovarian PARP inhibitor time-lapse assay, fitted simulations and in vivo pilot modulation/skipping; in_vitro_and_mouse | metadata_rights_or_access | Dose modulation favorable, skipping less successful. Public PMC access alone is not permission; no numeric rows. |
| [silva_2012: Evolutionary approaches to prolong progression-free survival in breast cancer](https://pmc.ncbi.nlm.nih.gov/articles/PMC3525750/) | Breast resistance fitness assays parameterize simulated therapy comparisons; in_vitro_and_model | metadata_rights_or_access | Headline progression gain derives from computational models, not independent observed AT animal/clinical survival. Fitness experiments relevant but no data copied. |
| [gallaher_2018: Spatial Heterogeneity and Evolutionary Dynamics Modulate Time to Recurrence in Continuous and Adaptive Cancer Therapies](https://aacrjournals.org/cancerres/article/78/8/2127/633375/Spatial-Heterogeneity-and-Evolutionary-Dynamics) | MCF7/MCF7Dox competition assays plus spatial ABM comparisons; in_vitro_and_model | metadata_rights_or_access | Continuous/adaptive outcome advantages mainly simulated; assay empirical but no numeric transcription due rights. |
| [das_thakur_2013: Modelling vemurafenib resistance in melanoma reveals a strategy to forestall drug resistance](https://pmc.ncbi.nlm.nih.gov/articles/PMC3930354/) | Two melanoma PDX models; intermittent vs continuous BRAF inhibition; mouse_xenograft | metadata_rights_or_access | Positive intermittent result linked to drug addiction, not independently identified competition; number transcription excluded by rights rule. |
| [maltas_2024: Drug dependence in cancer is exploitable by optimally constructed treatment holidays](https://pmc.ncbi.nlm.nih.gov/articles/PMC10918730/) | Melanoma drug-dependent/sensitive mixtures; intermittent laboratory schedules and model-derived adaptive schedules; in_vitro_and_model | metadata_rights_or_access | Author correction DOI10.1038/s41559-023-02300-9 exists; no numeric extraction. Distinguish experimentally tested holidays from modeled adaptive optimization. |
| [wang_2021_fixed_periodic: Is the Fixed Periodic Treatment Effective for the Tumor System without Complete Information?](https://pmc.ncbi.nlm.nih.gov/articles/PMC8643150/) | H1975/A549 osimertinib/FK866 periodic cycles; nonlinear mixed-effects model; mouse_and_model | metadata_rights_or_access | Same cells and authors as reachability paper, cohort overlap unresolved; some survival comparisons are simulations, not animals. |
| [kam_2015: Sweat but no gain: Inhibiting proliferation of multidrug resistant cancer cells with ersatzdroges](https://pmc.ncbi.nlm.nih.gov/articles/PMC4262548/) | P-glycoprotein substrate energetic-cost experiments; in_vitro_and_mouse | context_only | Manipulates resistance fitness cost rather than a direct adaptive-versus-high-dose scheduling comparison. |
| [kam_2014_review: Evolutionary strategy for systemic therapy of metastatic breast cancer: balancing response with suppression of resistance](https://pmc.ncbi.nlm.nih.gov/articles/PMC4258899/) | Review with early breast xenograft programme illustrations; review | overlap_or_review | Author manuscript2015; final publication2014; potential overlap with Enriquez2016, not counted as independent. |
| [antifragile_2026: Evolutionary antifragile therapy is a treatment strategy to suppress drug resistance by exploiting dose response convexity](https://pmc.ncbi.nlm.nih.gov/articles/PMC13283037/) | Cancer Research86(17):4434; DOI10.1158/0008-5472.CAN-25-4262; variable-dose antifragility; experimental_candidate | metadata_rights_or_access | Retained current experimental candidate; no invented values or presumed comparison direction. Official EuropePMC fullTextXML returned500; no numeric values or reuse license verified. |
| [plasticity_2025: Impacts of competition and phenotypic plasticity on the viability of adaptive therapy](https://pmc.ncbi.nlm.nih.gov/articles/PMC12515544/) | Mathematical Biosciences DOI10.1016/j.mbs.2025.109522; preprintPMC11974694; model | model_only | No independent experimental scheduling outcome identified; do not count fitted/model conclusions as additional experiments. |
| [microbial_2026: A microbial therapy-mimicry assay shows how spatial resource dynamics control resistance escape](https://www.nature.com/articles/s41559-026-03178-z) | Microbial spatial resources and therapy schedules; yeast | context_only | Mechanistically informative but yeast is not a real tumour comparison. |
| [gevertz_2026: Delaying cancer progression by integrating toxicity constraints in a model of adaptive therapy](https://www.nature.com/articles/s41540-025-00635-6) | Toxicity-aware mathematical adaptive treatment comparisons; model | model_only | No independent empirical cohort or experiment. |
| [enriquez_2014_ismrm: Early adaptive breast therapy experiment proceedings](https://cds.ismrm.org/protected/14MProceedings/PDFfiles/1130.pdf) | ISMRM2014 breast adaptive experiments; conference_abstract | overlap_or_review | Early breast programme; potential overlap with2016 paper, not independent. |
| [enriquez_2015_aacr: Abstract B11: Adaptive therapy modifies tumor](https://aacrjournals.org/mct/article/14/12_Supplement_1/B11/285815/Abstract-B11-Adaptive-therapy-modifies-tumor) | AACR2015 early breast adaptive experiments; conference_abstract | overlap_or_review | Early breast programme; potential overlap with2016 paper, not independent. |
| [colorectal_2026_abstract: Harnessing clonal dynamics to overcome resistance: An adaptive therapy framework in colorectal cancer](https://aacrjournals.org/cancerres/article/86/7_Supplement/702/776389/Abstract-702-Harnessing-clonal-dynamics-to) | AACR2026 abstract0702; immunocompetent multi-drug adaptive model; conference_abstract | unresolved_candidate | Current candidate found in broad search; numerical comparison not yet verified. |
| [qi_2021: Adriamycin-resistant cells are significantly less fit than adriamycin-sensitive cells in cervical cancer](https://pmc.ncbi.nlm.nih.gov/articles/PMC7874629/) | HeLa/ADR versus parental growth and competition assays; in_vitro_and_mouse | context_only | No adaptive versus high-dose schedule experiment; growth/competition results graph-based and not digitized. |
| [seyedi_2023_preprint: Testing Adaptive Therapy Protocols using Gemcitabine and Capecitabine on a Mouse Model of Endocrine-Resistant Breast Cancer](https://pmc.ncbi.nlm.nih.gov/articles/PMC10541126/) | Earlier report of Seyedi2024 experiment; preprint | excluded_scope | 2024 final article supersedes preprint; not independent. |
| [hockings_2023_preprint: Adaptive therapy achieves long-term control of chemotherapy resistance in high grade ovarian cancer](https://pmc.ncbi.nlm.nih.gov/articles/PMC10401956/) | Earlier report of Hockings2025 programme; preprint | excluded_scope | Same programme and experimental design; do not count as independent. |
| [antifragile_2025_preprint: Evolutionary antifragile therapy](https://pmc.ncbi.nlm.nih.gov/articles/PMC12458151/) | High/low intermittent versus even dosing; precedes West2026 CancerResearch; preprint | excluded_scope | Final2026 report candidate retained separately; no independent count. |
| [chattopadhyay_2024_abstract: Abstract PR003: Multi-drug adaptive therapy protocol in vitro: An exercise in optimization for neuroblastoma](https://doi.org/10.1158/1538-7445.CANEvol23-PR003) | Neuroblastoma sknbe2c/lan1 cisplatin/vincristine assay plus stochastic model; conference_abstract | excluded_scope | Abstract says comparative protocol efficacy still being tested; no quantitative efficacy comparison reported. |
| [seyedi_2022_model: In Silico Investigations of Multi-Drug Adaptive Therapy Protocols](https://pmc.ncbi.nlm.nih.gov/articles/PMC9179496/) | Hybrid agent-based multidrug adaptive protocols; model | model_only | Simulation precursor; not an independent experimental comparison. |
| [acar_2020_steering: Exploiting evolutionary steering to induce collateral drug sensitivity in cancer](https://pmc.ncbi.nlm.nih.gov/articles/PMC7174377/) | Cancer drug-sequence selection and collateral sensitivity; in_vitro | context_only | Evolutionary steering with different drugs; not directly matched containment versus maximum-dose comparison, retained contextual screen. |
| [strobl2022: Spatial structure impacts adaptive therapy by shaping intra-tumoral competition](https://pmc.ncbi.nlm.nih.gov/articles/PMC9053239/) | spatial_ABM_fitted_to_IADT_cohort; human_prostate_fit | included_model | Observed PSA used for fitting; adaptive-versus-continuous comparisons synthetic; overlap Brady2020/Strobl2021 |
| [brady2020: Prostate-specific antigen dynamics predict individual responses to intermittent androgen deprivation](https://www.nature.com/articles/s41467-020-15424-4) | stem_cell_PSA_model_fitted_to_IADT_cohort; human_prostate_fit | included_model | Published fit aggregates and counterfactual policy means; no independent observed adaptive/continuous arms |
| [brady2022_rbat: Range-Bounded Adaptive Therapy in Metastatic Prostate Cancer](https://pmc.ncbi.nlm.nih.gov/articles/PMC9657943/) | fitted_counterfactual_schedule_simulations; human_mCRPC_fit | included_model | Same Zhang cohort; RBAT alternatives simulated only |
| [kim2021: Adaptive Therapy for Metastatic Melanoma: Predictions from Patient Calibrated Mathematical Models](https://www.mdpi.com/2072-6694/13/4/823) | fitted_counterfactual_LV_and_switching_models; human_melanoma_fit | included_model | Eight actual continuous-treatment subjects; no observed adaptive arm; 50 fits per subject are not independent participants |
| [strobl2021: Turnover modulates the need for a cost of resistance in adaptive therapy](https://pmc.ncbi.nlm.nih.gov/articles/PMC8455086/) | fitted_counterfactual_model; human_prostate_fit | metadata_rights_or_access | NIH author manuscript no explicit redistributable license located; publisher denied; overlap Strobl2022/Brady2020; cost and turnover not directly measured. DOI10.1158/0008-5472.CAN-20-0806. |
| [west2020: Towards multi-drug adaptive therapy](https://pmc.ncbi.nlm.nih.gov/articles/PMC7307613/) | selected_case_counterfactual_model; human_mCRPC_fit | metadata_rights_or_access | NIH manuscript without explicit permissive license located; four selected case studies; model constants fixed except initial conditions; DOI10.1158/0008-5472.CAN-19-2669. |
| [brady2021: Predicting patient-specific response to adaptive therapy in metastatic castration-resistant prostate cancer using prostate-specific antigen dynamics](https://pmc.ncbi.nlm.nih.gov/articles/PMC8322456/) | prognostic_model; human_mCRPC_fit | metadata_rights_or_access | CC BY-NC-ND4.0; same16 cohort; conservative policy excludes derived CSV; DOI10.1016/j.neo.2021.06.013. |
| [gallagher2026: Mathematical Biomarkers of Adaptive Therapy Outcomes in Prostate Cancer](https://jamanetwork.com/journals/jamaoncology/fullarticle/2852635) | prognostic_model; human_prostate_fit | metadata_rights_or_access | CC BY-NC-ND and explicit restriction on alterations/commercial TDM/AI; no numerical rows reused; DOI10.1001/jamaoncol.2026.2781. |
| [pasetto2022: Intermittent Hormone Therapy Models Analysis and Bayesian Model Comparison for Prostate Cancer](https://pmc.ncbi.nlm.nih.gov/articles/PMC8604892/) | model_comparison; human_prostate_fit | excluded_scope | CC BY4.0; model comparison existing data, no new observed treatment comparison or identifiable PR61 quantities; DOI10.1007/s11538-021-00953-w. |
| [west2019: Multidrug Cancer Therapy in Metastatic Castrate-Resistant Prostate Cancer: An Evolution-Based Strategy](https://pmc.ncbi.nlm.nih.gov/articles/PMC6665681/) | selected_case_counterfactual_model; human_mCRPC_fit | metadata_rights_or_access | Selected individual model cases; no group comparative outcome; rights not further audited; DOI10.1158/1078-0432.CCR-19-0006. |
| [hansen2020_modifying: Modifying Adaptive Therapy to Enhance Competitive Suppression](https://pmc.ncbi.nlm.nih.gov/articles/PMC7761372/) | mathematical_model; synthetic | model_only | Borrowed model rates; no new empirical cohort fit. DOI10.3390/cancers12123556. |
| [salvioli2024: The effect of tumor composition on the success of adaptive therapy: The case of metastatic Castrate-Resistant Prostate Cancer](https://doi.org/10.1371/journal.pone.0308173) | mathematical_model; synthetic | model_only | Prior assumed parameter matrices/PSA contributions; no independent empirical fit. |
| [model_deep_learning2024: Mathematical Model-Driven Deep Learning Enables Personalized Adaptive Therapy](https://pmc.ncbi.nlm.nih.gov/articles/PMC11148552/) | fitted_counterfactual_algorithm; human_prostate_fit | unresolved_candidate | Main full text blocked during audit; unresolved algorithm/model candidate, no numerical reuse. |
| [deep_rl2024: Deep reinforcement learning identifies personalized intermittent androgen deprivation therapy for prostate cancer](https://pmc.ncbi.nlm.nih.gov/articles/PMC11174533/) | fitted_counterfactual_algorithm; human_prostate_fit | metadata_rights_or_access | CC BY-NC4.0; algorithmic counterfactuals no observed treatment comparison; not fully extracted. |
| [uncertainty2026: Adaptive therapy under parametric, structural, and measurement uncertainty](https://arxiv.org/abs/2608.18387) | preprint_model_uncertainty; human_prostate_fit | unresolved_candidate | Preprint; Bayesian model uncertainty, no new observed comparison; not fully extracted. |
| [gallaher2025_psa_testosterone: PSA/testosterone ratio as a potential biomarker to identify early progressors of adaptive therapy in metastatic castration sensitive prostate cancer](https://pmc.ncbi.nlm.nih.gov/articles/PMC12248141/) | single_arm_prognostic_model_preprint; human_mCSPC_fit | metadata_rights_or_access | CC BY-NC4.0; same single-arm feasibility cohort no independent continuous comparator; DOI10.1101/2025.05.08.25327179. |
| [gallaher2022_sum: The sum and the parts: dynamics of multiple and individual metastases during adaptive therapy](https://jeffreywest.github.io/pdf/gallaher2022sum.pdf) | fitted_or_analytical_model; model | metadata_rights_or_access | Preprint; all rights reserved/no reuse; model of multiple metastases using existing cohorts, no new independent observed comparison. |
| [gallagher2025_timing: Deriving Optimal Treatment Timing for Adaptive Therapy: Matching the Model to the Tumor Dynamics](https://pmc.ncbi.nlm.nih.gov/articles/PMC12417256/) | fitted_or_analytical_model; model | model_only | Analytical modelling and selected existing-cohort cases; no new observed treatment comparison. |

<!-- END RECOMPUTED EVIDENCE -->
