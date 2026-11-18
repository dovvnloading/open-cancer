# Which older-age breast-incidence declines survive uncertainty and stage accounting?

**Nonclinical, unreviewed computational research. No screening or medical advice.
No patient-level data, clinical prediction, or biological validation.**

**Verdict:** the overall incidence declines after the age group with the
highest rate survive the frozen simultaneous count-based test. The recorded
localized-stage losses account for most of those net declines. A decline is
not established in regional or distant disease by the primary two-sided rule;
the regional result on the later edge is sensitive to the prespecified tail
choice. Distant-stage incidence increases on the earlier declining-total edge.
The results do not support treating all stages as sharing the same decline.

The strongest missing-stage-robust finding is the localized decline from
75-79 to 80-84. It survives both Poisson uncertainty and the adverse allocation
of every older unknown-stage case and unassigned residual to localized
disease. The earlier localized decline does not survive that combined stress.
Its apparent majority contribution can also be erased by allocation alone.
The later localized majority survives allocation **as a descriptive share**;
we do not establish a sampling confidence statement for the majority fraction.
Localized stage is anatomic extent and does not establish screen detection.

This executes ranked job 2 from the [research programme](../README.md), after
the integration in [PR #77](https://github.com/dovvnloading/open-cancer/pull/77).
It extends the count-free observations in
[PR #64](../multistage-incidence-identifiability/) and
[PR #73](../incidence-strata-test/), with different, partly overlapping
aggregates. Race-specific sampling inference remains unavailable because the
matching race denominators are not printed here.

## Frozen design and audit trail

The branch starts at current-main commit
[`09d916c`](https://github.com/dovvnloading/open-cancer/commit/09d916cdea7053f25643c6d0bd9c45827436e315).
[PROTOCOL.md](PROTOCOL.md) was committed alone locally as
`523133aeb2949e82aa5a0411a61ce5cce9f486a8`
**before any numerical reconciliation, inference, decomposition, or execution
of result-generating code**. Source reading, factual transcription, and method
review occurred before the freeze. No protocol amendment or analytical
choice changed. There is one operational deviation for hypothetical malformed
inputs: the loader stops the whole run with an explicit error rather than
printing unaffected contrasts alongside per-contrast abstentions. All required
cells are present in the audited CSV, so this affects no reported result.
This is a committed computational specification for an already published
pattern, not blinded work or independent registration.

**Publication transport:** the shell's Git push had no configured GitHub
credential. The connected GitHub app therefore uploaded the same protocol
tree in [commit `1955c13`](https://github.com/dovvnloading/open-cancer/commit/1955c13859a28ef718c9dbc80ebf447788adb4e2),
followed by the analysis. The imported tree is identical to the original
pre-analysis tree. Its later server timestamp records upload, not the original
freeze, and the original local commit is not claimed to be a remotely
registered commit. The original Git commit object is preserved here; the test
recomputes its SHA-1 as well as the separate normalized protocol-content hash.
These metadata document the workflow and do not supply independent timestamp
attestation or scientific validation.

<!-- BEGIN ORIGINAL PROTOCOL COMMIT -->
```text
tree eacf57b706c3ab22b411cc785b635f20e97e6678
parent 09d916cdea7053f25643c6d0bd9c45827436e315
author Open Cancer Research <research@users.noreply.github.com> 1794989241 -0500
committer Open Cancer Research <research@users.noreply.github.com> 1794989241 -0500

research: freeze breast incidence stage analysis protocol
```
<!-- END ORIGINAL PROTOCOL COMMIT -->

The analysis fixes all adjacent comparisons among the finite age bins and
all printed stage endpoints, including increases and inconclusive outcomes.
The open-ended age group is descriptive. The exact two-sided primary family,
decline-only sensitivity, and hypothetical stage-completion sensitivity each
have the coverage interpretation specified in the protocol; there is no
claim of one common confidence guarantee spanning all sensitivity families.

## Data, source access, and reuse

The only numerical input is
[data/turbow_table1.csv](data/turbow_table1.csv), a small transcription of
selected factual cells from **Turbow SD, White MC, Breslau ES, Sabatino SA
(2021), Mammography use and breast cancer incidence among older U.S. women,
Breast Cancer Research and Treatment 188:307-316,
[doi:10.1007/s10549-021-06160-4](https://doi.org/10.1007/s10549-021-06160-4),
[PMC Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10846538/#T1)** [1].
The table describes invasive breast-cancer incidence in USCS for 2016,
using the November 2018 submission released in June 2019. The population
values come from the age-column headers, not the separate survey table.
Every CSV row includes the source, table, row, age column, denominator locator,
access date/timezone, units, adjustment label, database vintage, and reuse basis.

**Access on 2026-11-18, America/New_York:** the full PMC article HTML exposed
the table. An initial article request and a standalone-table request showed
a browser-check page; the full article text was subsequently readable. No
browser was used, no browser file download was attempted, and no download
prompt appeared. Search-engine PDF text snippets disagreed with several HTML
cells and were discarded as transcription evidence. We did not download or
parse the manuscript PDF or acquire the underlying USCS case records.

The CSV is an independently structured record of factual statistics with
attribution. It reproduces no article prose, figures, or table layout.
The [U.S. Copyright Office's facts/expression distinction][2] is the reuse
basis; [PMC's copyright notice][3] makes clear that free access does not
establish a broad redistribution license. We found no article-specific
Creative Commons grant and make no blanket public-domain claim. The
repository's MIT license does not relicense third-party material. This
limited factual transcription follows [DATA_POLICY.md](../../DATA_POLICY.md).

### Reconciling the age-adjustment label

Table 1 and its methods explicitly label the rates age adjusted to the
2000 U.S. standard population, with the 19-group Census P25-1130 standard [1].
The [NCI standard-population table][4] has a separate cell for each age group
used here. [Historical CDC guidance][5] explains that restricting the age
range renormalizes the standard weights. **Mathematical inference:** when
only one such standardization cell is selected, its normalized weight is
one, so its adjusted rate can equal its count/population rate. The complete
rounding reconciliation below is consistent with that interpretation.

We retain the source's label. We do not claim access to its original
SEER*Stat session or reproduce its interval algorithm. The actual inferential
estimand is the printed-count/population bin rate under a fixed exposure
model. Matching arithmetic does not authenticate the source or permit
substitution of raw counts for arbitrary age-adjusted rates in other tables.

### Source inconsistencies that remain visible

The overall counts and the sum of printed stage counts leave small nonzero
residuals at some ages. They are not silently added to unknown stage in the
data or base decomposition. Only the expressly conditional completion
sensitivity allows them to be assigned along with unknown stage.

The Results prose describes a sequence of percentages as distant-stage
proportions, although the checked counts reproduce that sequence for unknown
stage. The generated proportion table below shows the disagreement. The
source's later Discussion itself associates the oldest-group percentage with
unknown stage. This appears to be a category-label error; we retain it as an
unresolved source inconsistency and use the printed counts for arithmetic.

Race-specific population denominators and race-by-stage counts are absent.
The total female population cannot be used as a denominator for White or
Black cases. Excluded race, ethnicity, and receptor cells are not silently
repaired or used to complete the selected data. This note therefore supplies
no race-specific sampling verdict.

<!-- BEGIN RECOMPUTED RESULTS -->

## Recomputed results

The CSV contains **25 cells** across **5 age groups**. **25/25** count/population rates agree with the printed rates within half a last digit; **25/25** also lie within their printed intervals. These checks establish arithmetic compatibility, not source truth. The largest absolute rate discrepancy is **0.045420** per 100,000 (80-84, unknown).

### Every transcribed count/rate reconciliation

Rates are cases per 100,000 women. The source interval is its printed marginal 95% interval, not our simultaneous contrast interval. Gap = count/population rate minus printed rate. The CSV retains the original decimal strings and provenance.

| Age | Endpoint | Cases | Population | Computed rate | Printed rate | Printed interval | Gap | Rounding match |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 65-69 | overall | 37,127 | 8,889,353 | 417.656943 | 417.7 | [413.4, 421.9] | -0.043057 | yes |
| 65-69 | localized | 26,034 | 8,889,353 | 292.867209 | 292.9 | [289.3, 296.4] | -0.032791 | yes |
| 65-69 | regional | 8,303 | 8,889,353 | 93.403873 | 93.4 | [91.4, 95.4] | +0.003873 | yes |
| 65-69 | distant | 1,955 | 8,889,353 | 21.992602 | 22.0 | [21.0, 23.0] | -0.007398 | yes |
| 65-69 | unknown | 833 | 8,889,353 | 9.370761 | 9.4 | [8.7, 10.0] | -0.029239 | yes |
| 70-74 | overall | 29,147 | 6,360,203 | 458.271536 | 458.3 | [453.0, 463.6] | -0.028464 | yes |
| 70-74 | localized | 20,970 | 6,360,203 | 329.706457 | 329.7 | [325.3, 334.2] | +0.006457 | yes |
| 70-74 | regional | 6,012 | 6,360,203 | 94.525285 | 94.5 | [92.2, 96.9] | +0.025285 | yes |
| 70-74 | distant | 1,490 | 6,360,203 | 23.426925 | 23.4 | [22.3, 24.6] | +0.026925 | yes |
| 70-74 | unknown | 675 | 6,360,203 | 10.612869 | 10.6 | [9.8, 11.4] | +0.012869 | yes |
| 75-79 | overall | 20,575 | 4,647,676 | 442.694370 | 442.7 | [436.7, 448.8] | -0.005630 | yes |
| 75-79 | localized | 14,496 | 4,647,676 | 311.897817 | 311.9 | [306.8, 317.0] | -0.002183 | yes |
| 75-79 | regional | 4,250 | 4,647,676 | 91.443552 | 91.4 | [88.7, 94.2] | +0.043552 | yes |
| 75-79 | distant | 1,277 | 4,647,676 | 27.476098 | 27.5 | [26.0, 29.0] | -0.023902 | yes |
| 75-79 | unknown | 552 | 4,647,676 | 11.876904 | 11.9 | [10.9, 12.9] | -0.023096 | yes |
| 80-84 | overall | 13,433 | 3,413,802 | 393.490894 | 393.5 | [386.9, 400.2] | -0.009106 | yes |
| 80-84 | localized | 9,100 | 3,413,802 | 266.564962 | 266.6 | [261.1, 272.1] | -0.035038 | yes |
| 80-84 | regional | 2,921 | 3,413,802 | 85.564423 | 85.6 | [82.5, 88.7] | -0.035577 | yes |
| 80-84 | distant | 952 | 3,413,802 | 27.886796 | 27.9 | [26.1, 29.7] | -0.013204 | yes |
| 80-84 | unknown | 459 | 3,413,802 | 13.445420 | 13.4 | [12.2, 14.7] | +0.045420 | yes |
| 85+ | overall | 12,263 | 4,154,029 | 295.207376 | 295.2 | [290.0, 300.5] | +0.007376 | yes |
| 85+ | localized | 7,598 | 4,154,029 | 182.906764 | 182.9 | [178.8, 187.1] | +0.006764 | yes |
| 85+ | regional | 2,817 | 4,154,029 | 67.813682 | 67.8 | [65.3, 70.4] | +0.013682 | yes |
| 85+ | distant | 1,036 | 4,154,029 | 24.939643 | 24.9 | [23.4, 26.5] | +0.039643 | yes |
| 85+ | unknown | 811 | 4,154,029 | 19.523215 | 19.5 | [18.2, 20.9] | +0.023215 | yes |

### Unassigned residuals

Residual = overall minus the four printed stage counts. The residual is kept separate from unknown stage. The possible-allocation pool is unknown plus residual, conditional on these being assignable cases.

| Age | Overall | Stage sum | Residual | Unknown + residual |
| --- | --- | --- | --- | --- |
| 65-69 | 37127 | 37125 | 2 | 835 |
| 70-74 | 29147 | 29147 | 0 | 675 |
| 75-79 | 20575 | 20575 | 0 | 552 |
| 80-84 | 13433 | 13432 | 1 | 460 |
| 85+ | 12263 | 12262 | 1 | 812 |

### Primary simultaneous adjacent-age inference

There are **15** frozen comparisons. The exact conditional Poisson RR intervals use error **0.0033333333** per comparison (**0.0016666667** in each tail), giving at least **95%** simultaneous coverage within this family. RR is older/younger; change is older minus younger per 100,000. Decisions use unrounded endpoints. 'Unresolved' does not mean equal rates.

| Age contrast | Endpoint | Rate change | RR | Simultaneous RR interval | Adjusted two-sided p | Primary direction |
| --- | --- | --- | --- | --- | --- | --- |
| 65-69 -> 70-74 | overall | +40.614593 | 1.097244 | [1.072284, 1.122764] | 4.07115e-31 | increase |
| 65-69 -> 70-74 | localized | +36.839249 | 1.125788 | [1.095480, 1.156907] | 5.74762e-36 | increase |
| 65-69 -> 70-74 | regional | +1.121412 | 1.012006 | [0.962737, 1.063671] | 1 | unresolved |
| 65-69 -> 70-74 | distant | +1.434323 | 1.065218 | [0.962161, 1.178827] | 1 | unresolved |
| 65-69 -> 70-74 | unknown | +1.242108 | 1.132551 | [0.971132, 1.319831] | 0.264733 | unresolved |
| 70-74 -> 75-79 | overall | -15.577166 | 0.966009 | [0.940476, 0.992198] | 0.0021985 | decline |
| 70-74 -> 75-79 | localized | -17.808640 | 0.945986 | [0.916386, 0.976490] | 4.04727e-06 | decline |
| 70-74 -> 75-79 | regional | -3.081733 | 0.967398 | [0.911867, 1.026128] | 1 | unresolved |
| 70-74 -> 75-79 | distant | +4.049173 | 1.172843 | [1.047710, 1.312540] | 0.000494808 | increase |
| 70-74 -> 75-79 | unknown | +1.264035 | 1.119104 | [0.943596, 1.326107] | 0.805317 | unresolved |
| 75-79 -> 80-84 | overall | -49.203477 | 0.888855 | [0.860299, 0.918298] | 2.43114e-25 | decline |
| 75-79 -> 80-84 | localized | -45.332855 | 0.854655 | [0.821637, 0.888907] | 5.7857e-31 | decline |
| 75-79 -> 80-84 | regional | -5.879128 | 0.935708 | [0.871596, 1.004258] | 0.0873122 | unresolved |
| 75-79 -> 80-84 | distant | +0.410698 | 1.014947 | [0.893920, 1.151560] | 1 | unresolved |
| 75-79 -> 80-84 | unknown | +1.568516 | 1.132064 | [0.938062, 1.364869] | 0.805692 | unresolved |

**Primary verdict:** 2 overall declines survive: 70-74 -> 75-79; 75-79 -> 80-84. Localized disease has 2 robust declines on 70-74 -> 75-79; 75-79 -> 80-84. Regional/distant disease has 0 robust declines under this two-sided simultaneous rule. This last finding is limited evidence, not proof of no regional or distant decline.

### Signed decomposition of the observed change

Loss = younger minus older rate. Positive values contribute to decline; negative values offset it. Residuals stay explicit. Shares divide by the net overall loss only when it is positive. These are arithmetic contributions, not causal fractions or confidence bounds.

| Age contrast | Overall loss | Localized | Regional | Distant | Unknown | Residual |
| --- | --- | --- | --- | --- | --- | --- |
| 65-69 -> 70-74 | -40.614593 | -36.839249 | -1.121412 | -1.434323 | -1.242108 | +0.022499 |
| 70-74 -> 75-79 | +15.577166 | +17.808640 | +3.081733 | -4.049173 | -1.264035 | +0.000000 |
| 75-79 -> 80-84 | +49.203477 | +45.332855 | +5.879128 | -0.410698 | -1.568516 | -0.029293 |

| Declining age contrast | Localized share | Regional share | Distant share | Unknown share | Residual share | Mainly localized (recorded) |
| --- | --- | --- | --- | --- | --- | --- |
| 70-74 -> 75-79 | +114.325% | +19.784% | -25.994% | -8.115% | +0.000% | yes |
| 75-79 -> 80-84 | +92.133% | +11.949% | -0.835% | -3.188% | -0.060% | yes |

### Adverse allocation of unknown stage and residuals

For each named stage separately, assign all older unknown/residual cases to it and none younger to maximize the RR; reverse for the minimum. The interval envelope also includes the frozen sampling uncertainty. These are conditional completion bounds, using the same conservative multiplicity factor in a separate sensitivity family. They do not presume all adverse stage allocations occur together.

| Age contrast | Stage | Possible point RR range | Simultaneous RR envelope | Decline survives |
| --- | --- | --- | --- | --- |
| 65-69 -> 70-74 | localized | [1.090802, 1.162026] | [1.061636, 1.193866] | no |
| 65-69 -> 70-74 | regional | [0.919532, 1.125629] | [0.875601, 1.181365] | no |
| 65-69 -> 70-74 | distant | [0.746417, 1.547784] | [0.678695, 1.697170] | no |
| 70-74 -> 75-79 | localized | [0.916486, 0.982009] | [0.887987, 1.013326] | no |
| 70-74 -> 75-79 | regional | [0.869747, 1.093046] | [0.820830, 1.157095] | no |
| 70-74 -> 75-79 | distant | [0.807176, 1.679819] | [0.726996, 1.862472] | no |
| 75-79 -> 80-84 | localized | [0.823304, 0.897857] | [0.791716, 0.933298] | yes |
| 75-79 -> 80-84 | regional | [0.828146, 1.083063] | [0.772679, 1.159080] | no |
| 75-79 -> 80-84 | distant | [0.708632, 1.505363] | [0.629268, 1.687498] | no |

| Declining age contrast | Localized share across allocations | Majority survives allocation |
| --- | --- | --- |
| 70-74 -> 75-79 | [38.080%, 182.456%] | no |
| 75-79 -> 80-84 | [64.748%, 116.272%] | yes |

The share ranges condition on the observed totals and quantify allocation uncertainty only. They are **not sampling confidence intervals** for the majority claim. The RR envelopes, separately, incorporate the assumed Poisson variability.

### Prespecified one-sided sensitivity

These decline-only upper limits spend the full per-comparison error in one tail. They do not replace the primary two-sided results. Every primary comparison is shown so the sensitivity does not select favorable endpoints.

| Age contrast | Endpoint | One-sided simultaneous RR upper | Decline-only decision |
| --- | --- | --- | --- |
| 65-69 -> 70-74 | overall | 1.120815 | unresolved |
| 65-69 -> 70-74 | localized | 1.154527 | unresolved |
| 65-69 -> 70-74 | regional | 1.059687 | unresolved |
| 65-69 -> 70-74 | distant | 1.169890 | unresolved |
| 65-69 -> 70-74 | unknown | 1.304794 | unresolved |
| 70-74 -> 75-79 | overall | 0.990196 | decline |
| 70-74 -> 75-79 | localized | 0.974154 | decline |
| 70-74 -> 75-79 | regional | 1.021584 | unresolved |
| 70-74 -> 75-79 | distant | 1.301491 | unresolved |
| 70-74 -> 75-79 | unknown | 1.309376 | unresolved |
| 75-79 -> 80-84 | overall | 0.916042 | decline |
| 75-79 -> 80-84 | localized | 0.886277 | decline |
| 75-79 -> 80-84 | regional | 0.998931 | decline |
| 75-79 -> 80-84 | distant | 1.140714 | unresolved |
| 75-79 -> 80-84 | unknown | 1.345921 | unresolved |

Changed decline decisions: **1**. 75-79 -> 80-84, regional

### Open-ended age group and source-prose check

The 85+ cells are reported descriptively in the reconciliation, residual, and proportion tables; no hypothesis test or finite-bin midpoint is assigned to them. Recomputing stage proportions exposes the article's prose/category conflict: its sequence described as distant-stage percentages matches the unknown-stage count shares after rounding. Here both are calculated from the transcribed counts.

| Age | Distant / overall | Unknown / overall |
| --- | --- | --- |
| 65-69 | 5.266% | 2.244% |
| 70-74 | 5.112% | 2.316% |
| 75-79 | 6.207% | 2.683% |
| 80-84 | 7.087% | 3.417% |
| 85+ | 8.448% | 6.613% |

### Newline-normalized source and code hashes

SHA-256 uses UTF-8 after CRLF/CR -> LF normalization; path keys use `Path.as_posix()`. The README is excluded to avoid self-reference. A hash pins content, not source accuracy or scientific validity.

| Repository path | SHA-256 |
| --- | --- |
| research/breast-incidence-decline-stage/PROTOCOL.md | `ff029700849020859911e9f8a1ace73eaa9acd6939e1cba912489a1f96930b14` |
| research/breast-incidence-decline-stage/data/turbow_table1.csv | `49d90dd3ea6eb9c7cffb5bf24bc0ea45491676469d9624676fe250d03797eeb0` |
| research/breast-incidence-decline-stage/analysis.py | `08d31bb17cfc02f0f9e993c7200004f1d52fb1546b2b4b963e76b5f032ace6ce` |

<!-- END RECOMPUTED RESULTS -->

## What the uncertainty calculations mean

For each fixed age pair, independent Poisson counts conditional on fixed
population exposures imply a binomial distribution for the older count
conditional on the sum. We invert the binomial tails and transform their
Clopper-Pearson endpoints to older/younger incidence-rate ratios. Bonferroni
controls the primary family's simultaneous coverage without requiring
independence between adjacent comparisons or between stage and overall
endpoints. The calculation uses exact distributional tails evaluated in
floating point, not an asymptotic normal approximation or interval-certified
arithmetic. The test checks the endpoints using an independent probability
recurrence.

The registries aim to enumerate cases. Thus “sampling uncertainty” here means
the uncertainty in a hypothetical repeated count-generating process, not
uncertainty from a sampled survey. Fixed population estimates are a proxy for
annual exposure; they are not measured cancer-free person-years. Population
error, registry omissions, staging error, dependence, and extra-Poisson
variation can invalidate the calibration. They cannot be estimated from this
table, and the reported coverage is conditional on the declared model.

The source's printed marginal intervals and the half-digit rounding checks
serve different purposes from our simultaneous intervals. Nonoverlap of
source intervals is not the decision rule. A narrower decline-only family
changes the regional decision, which is why the primary claim stays tied to
the frozen two-sided rule. A nonsignificant regional result does not erase
its observed negative change or demonstrate no underlying decline.

The completion sensitivity is a form of partial identification. It assumes
that the total and known-stage counts are correct, and that all currently
unassigned cases could be distributed among the named stages. For any actual
completion within those restrictions, its interval is contained in the
computed envelope. Taking the union is conservative; it does not require
missing-at-random staging. If residuals instead represent source errors, the
completion interpretation is only hypothetical. The printed counts alone
cannot choose among those explanations.

## Falsification attempts and outcomes

| Challenge to the proposed result | What was checked and what survived |
|---|---|
| A transcription or exposure error creates the decline | Every selected numeric CSV field was independently re-read against the full PMC HTML, including the population headers and interval endpoints. The complete exact-rational rounding audit passes. This is a transcription check, not validation of the underlying registry. |
| The published age-adjustment label prevents count inference | The explicit standard's age cells and the raw-rate reconciliation support the single-cell interpretation. Inference remains explicitly conditional on the printed-count/population estimand; unpublished settings are not claimed as reproduced. |
| Multiplicity removes the apparent age decline | All frozen edges and endpoints remain in the primary family. Overall and localized declines survive on both declining-total edges; regional and distant declines do not satisfy that primary rule. |
| The decline is a common decrease across stages | The regional result is unresolved under the primary rule and sensitive to tail choice. Distant disease increases on the earlier edge and is unresolved on the later edge. The common-decrease description is unsupported. |
| Unknown staging and residual cases overturn localized attribution | Adverse assignment removes the earlier localized decline's sampling robustness and its descriptive majority guarantee. The later localized decline and descriptive majority survive the specified allocation stress. The majority fraction has no sampling confidence bound. |
| The increase in unknown stage can be ignored | Unknown stage and signed residuals remain explicit in both the identity and the sensitivity pool. Neither is renormalized away. |
| The open-ended oldest group strengthens a finite-bin result | It is excluded from hypothesis testing and supplied only as descriptive source/reconciliation information. |
| Stage labels or the paper's prose prove a screening mechanism | Localized is not measured detection mode. The prose/count discrepancy is exposed rather than copied. Neither the stage partition nor the separate survey supplies a causal contrast. |

The hostile review checks both adverse and favorable outcomes. The source
audit and mathematical/software review were conducted by separate AI
investigations with access to the same evidence. Their agreement is not
independent biological replication, human peer review, or clinical validation.

## What this can and cannot establish

**It can establish a conditional observable restriction failure.** Under the
fixed-exposure Poisson model, the primary family rejects nondecreasing
finite-bin incidence means at its declared error level. This gives a statistical
analogue of the earlier rounded-mean failure in different aggregates. It
tests a restriction on observed diagnosis rates, not a biological mechanism
alone. The recorded-stage decomposition locates which counts supply the net
rate loss, and the completion analysis shows how much that statement depends
on missing stage.

**It cannot identify a screening effect.** Screening history, mode of
detection, counterfactual diagnoses without screening, lead time, and
overdiagnosis are not observed. Localized cancers can be symptom-detected;
screen-detected cancers need not all be localized. The paper's screening
survey uses 2018 NHIS, whereas incidence is from 2016 USCS, with no individual
linkage [1]. Their coexisting age patterns do not identify the consequences
of screening cessation or provide an intervention comparison.

**It cannot separate aging, cohort, and competing mortality.** A single
calendar year's age bins contain different birth cohorts. Population
denominators do not reveal earlier exposures, individual risk trajectories,
survival selection, diagnosis delays, or cancer-free risk time. The analysis
cannot recover an unobserved onset hazard or attribute the decline to death,
selection, reduced detection, or changing disease occurrence. Those
explanations need additional, appropriately linked observations.

**It cannot count biological transitions.** Localized/regional/distant stage
at diagnosis is not the number of mutations or model transitions. PR #64's
monotone-mean theorem has homogeneous initiation, clone, and observation
assumptions. A violation of its compound observable restriction does not
choose the failed assumption or identify an alternative transition count.
PR #73's race-stratum pattern remains background evidence; its missing
count/exposure uncertainty is not retroactively supplied by pooled USCS data.

The cohorts also overlap: USCS incorporates SEER, and its incidence year is
within the period of the earlier SEER summaries. Broader geographic coverage
does not make this a fully independent replication [1; PR #64; PR #73].

## Novelty assessment

**This is an incremental reproducibility and uncertainty audit.** Turbow et
al. already reported the overall peak and subsequent decline, the dominance
of localized disease, the higher distant-stage incidence in the older finite
groups, and increasing unknown-stage incidence [1, Results and Discussion].
Rediscovering those patterns is not a new epidemiological finding.

The addition here is the committed contrast family; a complete count/exposure
and residual audit; exact simultaneous adjacent-age inference; signed
stage-loss accounting; and explicit worst-case stage-completion bounds. The
paper's described significance tests concern screening-use differences in
NHIS; it does not report this simultaneous USCS adjacent-age family or this
allocation analysis [1, Analysis]. That comparison is against the identified
paper, not a systematic literature search or a claim of bibliographic
priority. The exact conditional method and Bonferroni principle are
established statistics, not newly developed methods.

## Reproduce and verify

From the repository root, with Python's standard library only:

```bash
python -B research/breast-incidence-decline-stage/analysis.py --check
python -B research/breast-incidence-decline-stage/analysis.py --format json
python -B -m unittest discover -s tests -p 'test_research_breast_incidence_decline_stage.py' -v
```

[analysis.py](analysis.py) derives the numerical section directly from the
CSV; `--write-readme` regenerates that section deliberately. No network,
untracked result file, downloaded source document, or optional Python package
is required. The [dedicated test](../../tests/test_research_breast_incidence_decline_stage.py)
independently checks source invariants, rate and decomposition arithmetic,
exact interval tails, allocation extrema, frozen hashes, the generated
section, and newline/platform portability. It is a bounded workload intended
to finish within the requested runtime limit.

Only this new directory and the dedicated test are part of the change.
The existing PR workflow discovers the test without a workflow edit.
Repository validation and CI are software checks, not scientific peer review.
No publication-registry entry, clinical recommendation, or claim of reviewed
scientific status is created.

## Sources

1. Turbow et al. (2021), [article and Table 1][1]. Source for all aggregate
   values, stage definitions, database vintage, source labels, and comparison
   with the original analysis.
2. U.S. Copyright Office, [What Does Copyright Protect?][2]. Factual-statistic
   reuse basis; no claim of a license to reproduce article expression.
3. NLM, [PMC Copyright Notice][3]. Distinguishes free reading and applicable
   reuse rights from open licensing.
4. NCI, [Standard Populations: 19 Age Groups][4]. Historical standard's cell
   definitions; no switch to a later database's standard.
5. CDC, [historical USCS incidence documentation][5], age-adjusted-rate
   calculation and age-restriction notes. Used for the normalization principle,
   not as an additional data source or a substitute for the study's release.

All linked source pages were read as HTML on 2026-11-18. The original paper
is the evidence source; the other links clarify provenance and arithmetic.

[1]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10846538/#T1
[2]: https://www.copyright.gov/help/faq/faq-protect.html
[3]: https://pmc.ncbi.nlm.nih.gov/about/copyright/
[4]: https://seer.cancer.gov/stdpopulations/stdpop.19ages.html
[5]: https://wonder.cdc.gov/wonder/help/cancer-v2009.html
