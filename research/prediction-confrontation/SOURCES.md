# Search, eligibility and transcription audit

**Nonclinical, unreviewed. No treatment advice.** Access date and targeted
search cut-off: 2026-11-17 (America/New_York). The frozen protocol is in
local commit `29bedbd`; amendments and the remote publication mapping are
explicit in [AMENDMENTS.md](AMENDMENTS.md).
This ledger records a targeted falsification search, not systematic-review
completeness. Search results were used for discovery; retained numerical
rows were checked against primary text, printed aggregate labels or
identified aggregate supplement cells.

## Search coverage and stopping rule

Read the requested notes before searching. Search families covered clonal
neoantigen/heterogeneity response, HLA LOH and checkpoint outcomes including
null/contrary directions, TRACERx/Reinert/Henriksen lead time and thinning,
ecDNA conditional segregation/Fano/killing/dose, and resistance fluctuation
assays, drug holidays, count moments and lineage retention. Named author,
exact-title and DOI searches resolved primary copies and source terms.
References led to additional mechanistic and contrary studies below.

Examples of query families actually used include `cancer fluctuation assay
drug holiday variance colonies`, `resistance holiday factorial cumulants`,
`drug holiday mean variance colony cancer fluctuation`, `drug holiday Luria
zero`, `ecDNA segregation Fano factor conditional variance daughter mother
Lange Scanu Shibata`, and `ecDNA segregation generation two Fano dose
threshold measured cell killing aggregate`. HLA searches explicitly sought
favorable, null and adverse associations. No outcome-dependent cutoff,
histology or endpoint was selected to force agreement or rejection.

Both general search systems and primary-source follow-up were used. Some
PMC pages returned CAPTCHA/identity pages, some publisher endpoints returned
access errors, and one supplementary redirect failed. Ordinary open primary
publisher or institutional copies were used where available; no access
controls were bypassed. An inaccessible supplement is not claimed inspected.
No underlying patient datasets, mixed patient-level source-data workbooks,
individual figure points or reconstructed survival curves were used.

Stop after the user-named studies and relevant primary follow-up had either
supplied eligible aggregates or had an explicit mismatch/access/rights
reason, with a final search for the missing generation/rate and holiday
moment designs. The retained evidence has publication and search selection
bias. This audit cannot exclude the existence of a suitable unlocated table.

## Included numerical sources

### Immune response associations

**Litchfield et al. (2021), Cell**, [DOI](https://doi.org/10.1016/j.cell.2021.01.002),
[primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7933824/).
Locators: Results, “Landscape of published biomarkers”, Figure 2A; STAR
Methods, “Biomarker analysis”. Article-specific CC BY license checked.
All extracted clonal, total, subclonal and HLA point/interval/p fields were
re-read from the primary Results. Every feature, including a binary HLA-LOH
indicator, was standardized within study before the random-effects
meta-analysis. Thus the HLA magnitude is not the raw LOH-versus-intact OR.
Clonal mutations are not equivalent to functional presented neoantigens.
CR/PR versus SD/PD/NE is the response endpoint. The source includes cohorts
also used in earlier clonal/HLA papers; no independent replication count or
pooled effect is claimed. Total/subclonal comparators are retained without a
significance claim or an unsupported comparison of correlated coefficients.

**Ravi et al. (2023), Nature Genetics**, [DOI](https://doi.org/10.1038/s41588-023-01355-5),
[publisher](https://www.nature.com/articles/s41588-023-01355-5),
[primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10181943/).
Locators: clonal-neoantigen response Results near Figure 1f; Methods,
“Response association testing”. CC BY article license checked. Point/CI and
BH-adjusted q were re-read in the publisher text. Primary Methods confirm
univariate logistic CR/PR versus SD/PD, log-transformed mutation-burden
features and Benjamini–Hochberg correction. No base/offset is assumed for
a burden-doubling interpretation. A separate metric's natural-log statement
is not used to calibrate this feature. Pretreatment HLA LOH and related
presentation variables have no reported response association in the text;
no effect was inferred from a graph or an uninspected mixed workbook.
The source q is not copied into a raw-p field. Rounding and CI approximation
are explicit in the analysis.

**Lozac'hmeur et al. (2024), npj Precision Oncology**,
[DOI](https://doi.org/10.1038/s41698-024-00665-z),
[primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11300791/).
Locators: survival Results, Figure 3 caption, Cox Methods. CC BY-NC-ND
article license checked; only unchanged scalar facts are transcribed.
Both LUAD and LUSC are retained. HR is death hazard under HLA-A LOH versus
stable, with source univariable Cox models stratified by regimen and line;
this is not a TMB-adjusted causal estimate. Regimens include checkpoint
therapy with/without chemotherapy. The LUSC lower CI differs between abstract
and figure caption; the CSV consistently uses Figure 3c and flags that
choice. LOH-associated TMB differences and treatment composition prevent
interpreting the opposite LUAD estimate as causal protection or as a direct
refutation of Montesion's differently defined cohort.

**Yang, Kim & Kim (2022), Genes & Genomics**,
[DOI](https://doi.org/10.1007/s13258-021-01207-8),
[publisher](https://link.springer.com/article/10.1007/s13258-021-01207-8),
[main Figure 1](https://link.springer.com/article/10.1007/s13258-021-01207-8/figures/1).
The publisher's Rights and Permissions explicitly states CC BY. All printed
response/nonresponse segments in Figure 1c–e were visually read and rechecked;
all loci/cohorts are retained. These are unadjusted contingency counts.
The underlying cohorts/loci are dependent and overlap the broader CPI
literature; no pooled replication analysis is performed. Supplement S1,
which contains patient data, was not accessed.

Three source conflicts remain unresolved and are recorded on every row:
the abstract's total differs from Methods/Figure 1a; HLA-C response-panel
cell sums disagree with the text/Figure 1b total; and one Results sentence
swaps the melanoma treatment/sample-count mapping. Use the printed response
cells without repairing denominators. Methods and the figure caption agree
that panel d is anti-PD1 and panel e is anti-CTLA4, which is the mapping used.
A significant result in an internally inconsistent source is not silently
promoted to stronger evidence.

### ctDNA lead time and schedules

**Reinert et al. (2019), JAMA Oncology**,
[DOI](https://doi.org/10.1001/jamaoncol.2019.0528),
[corrected primary article](https://jamanetwork.com/journals/jamaoncology/fullarticle/2733132).
Locators: “Association of Longitudinal ctDNA Analysis With Patient Outcome”,
paragraph citing Figure 3A, and sampling Methods. Article Information states
CC BY; this is an article-specific license, not an assumption about all JAMA
content. The corrected publication was used. Mean, range, planned draw
interval and detected/total recurrence counts were checked. Lead mean is
conditional on detection; no mean CI is supplied. The planned cadence does
not establish exact realized spacing. No source patient trajectories were
read off the figure or supplement.

**Henriksen et al. (2022), Clinical Cancer Research**,
[DOI](https://doi.org/10.1158/1078-0432.CCR-21-2404),
[permitted institutional manuscript](https://findresearcher.sdu.dk/ws/portalfiles/portal/200565953/507.pdf).
Locators: printed page 513 lead-time and growth Results, and cohort imaging/
sampling Methods. CC BY-NC-ND license checked; unchanged numerical excerpts
only. Fast/slow multiplicative growth summaries are means and SDs from
**data-defined groups of fitted slopes**. Transforming their means with a
logarithm is not estimating a mean log-rate or independently testing a
mixture. For the CT-date restriction, the median/IQR refer to earlier-detected
cases among the eligible restricted cohort. This is the nearest observed
schedule-restriction check, but it is not parity thinning. Planned imaging
schedules are retained as design context, never asserted to bound actual
individual delays. The Danish cohort overlaps earlier work; no pooling.

**Henriksen et al. (2024), Annals of Oncology**,
[DOI](https://doi.org/10.1016/j.annonc.2023.11.009),
[primary manuscript](https://findresearcher.sdu.dk/ws/portalfiles/portal/255744575/1-s2.0-S0923753423050731-main.pdf).
Locator: printed page 233, sampling-frequency section. CC BY license checked.
The frequency-gain estimate is a growth-model extrapolation; it is not an
independently observed monthly schedule. The median/IQR refer to the subgroup
whose detection was estimated to be delayed under the observed cadence.
Retaining it tests the temptation to call model-derived outcomes independent
validation, not the schedule optimum itself.

**Abbosh et al. (2023), TRACERx, Nature**,
[DOI](https://doi.org/10.1038/s41586-023-05776-4),
[permitted accepted manuscript](https://eprints.whiterose.ac.uk/id/eprint/203072/1/Tracking%20early%20lung%20cancer%20metastatic%20dissemination%20in%20TRACERx%20using%20ctDNA.pdf).
Locators: lead-time Results, extended figure context and Methods “Lead time”
plus sampling schedule. CC BY license checked. Median/range/recurrence
count and early/late planned cadences were checked. Zeros are assigned when
no pre-relapse detection occurs, unlike Reinert's detected-only mean. It is
not legitimate to subtract the two studies' summaries as a schedule effect.

### ecDNA

**Lange et al. (2022), Nature Genetics**,
[DOI](https://doi.org/10.1038/s41588-022-01177-x),
[publisher](https://www.nature.com/articles/s41588-022-01177-x),
[primary supplement](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-022-01177-x/MediaObjects/41588_2022_1177_MOESM1_ESM.pdf).
CC BY license checked. Table SI1, printed page 5/PDF page 6, contains the
random-segregation KS statistic, p and plotted sample size. All panels were
transcribed and re-read. The final main Figure 1c has different p-values for
some apparently matching panels, including a panel with the same n. The
primary-hosted and author-hosted supplements agree; the discrepancy with the
main figure remains unresolved. No preferred representation is selected to
produce a verdict. The retained supplement is descriptive/source-audit
context only. Maternal-total mixtures and pooled sister fractions are not
conditional moments or a generation-two Fano. Patient-specific tumour
summaries and underlying individual cell observations were not extracted.

**Shibata et al. (2026), Nature Communications**,
[DOI](https://doi.org/10.1038/s41467-026-76566-5),
[publisher](https://www.nature.com/articles/s41467-026-76566-5),
[primary PDF](https://www.nature.com/articles/s41467-026-76566-5.pdf).
CC BY-NC-ND license is explicit on the article's final page. Only unchanged
scalar fitted parameters and printed aggregate test labels are retained.

Aggregate source cells and locators:

- [Supplementary Data 5](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-76566-5/MediaObjects/41467_2026_76566_MOESM7_ESM.xlsx): fitted p in B5:B16, rho in C5:C16, sheet of the same name. The source decimal precision is preserved, removing only spreadsheet binary-floating representation artifacts.
- [Supplementary Data 8](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-76566-5/MediaObjects/41467_2026_76566_MOESM10_ESM.xlsx): p in B5:B8, rho in C5:C8. All source arms retained. Source spelling distinctions such as siCtr1/siCtrl are preserved in CSV notes.
- Main Figure 3d, PDF page 5: every printed n/p label was visually rechecked. The caption calls these two-sided KS comparisons to simulated distributions. Both sisters' fractions are plotted; no independent-mother denominator or pair-aware calibration is asserted from n.
- [Supplementary Notes](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-76566-5/MediaObjects/41467_2026_76566_MOESM1_ESM.pdf), Note 3: a separate per-mother Monte Carlo procedure is described, but its numerical p-results were not found. The plotted KS values are not relabelled as that Monte Carlo test.

All fitted scalar cells were compared with the primary aggregate workbook
cells. No CI/SE/covariance for the fitted p/rho appears in the inspected
aggregate tables. The conditional linear Fano form is imposed by the fitted
beta-binomial family, not independently observed by this transformation.
The image-analysis Methods also trim high observations, further preventing
an unqualified exact-moment interpretation.

**Critical classification conflict:** main text on PDF page 3 calls clones
4, 6 and 7 predominantly HSR-like. Figure 3b clearly depicts 4/6 as HSR but
7 predominantly ecDNA-colored. Treat 4/6 as controls, and 7 as unresolved;
report the sensitivity with and without 7, never silently choose whichever
produces the desired p-value. Source classification and sister dependence
are not resolved by multiplicity correction. Other fitted model-comparison
tables were inspected for context but unnecessary values were not copied.

**Lundberg et al. (2008), PLOS ONE**,
[DOI](https://doi.org/10.1371/journal.pone.0003099),
[primary article](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003099).
The Copyright statement permits attributed redistribution under Creative
Commons Attribution; no license version is invented. Locator: Results,
“Statistical distribution of DMs at mitosis”, first paragraphs and Figure
3E–G context. Every untreated cell-line denominator, rejection percentage
and significance cutoff was re-read from the primary HTML. The CSV stores
those published quantities; integer rejected counts are derived only at
analysis time from the known denominators and an outward rounding envelope.
No individual mitoses or plotted points were digitized.

The paper does not document whether each per-mitosis p is a calibrated
cumulative-tail test or its one/two-sided construction. Divisions' assignment
to independent culture replicates is absent. Detached double minutes are
excluded, and clustering may complicate counting. These omissions block an
unconditional binomial-rejection-frequency interpretation. The irradiated
percentage summaries lack explicit corresponding denominators and were not
used. Colony tests based on presumed sisters/inferred division counts were
also excluded. The exploratory Markov sensitivity removes dependence only;
it retains the calibration assumption.

### Resistance

**Iyer et al. (2025), PLOS Computational Biology**,
[DOI](https://doi.org/10.1371/journal.pcbi.1013446),
[primary article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013446).
CC BY license checked. Locators: Results around Figure 6e and Table 1 caption.
The approximate experimental persister fraction, VMR and retained lineage
count were independently re-read. These are not four independent holiday
culture moments. Repeated random daughter choice reuses the same lineages
and is not biological replication. Fixed initial state and founder size are
not conditioned out. Model-generated VMR and fitted fate parameters are not
included. The word “approximately” does not justify a narrow decimal-rounding
confidence interval. A supplementary redirect was inaccessible, so no unseen
supplementary covariance or uncertainty is claimed.

## Metadata-only exclusions and negative eligibility checks

No numerical facts from rights-uncertain sources are copied to CSVs. Rights
exclusion is a conservative project-policy decision, not a claim that isolated
facts are universally copyrightable or that an inaccessible study has no
scientific evidence. Qualitative adverse/null evidence is retained to avoid
an apparent one-sided literature selection.

| Source | Primary locator/design screened | Why it cannot supply the planned aggregate test here |
| --- | --- | --- |
| [McGranahan 2016](https://doi.org/10.1126/science.aaf1490); [inspected manuscript](https://backend.orbit.dtu.dk/ws/portalfiles/portal/140198020/befo_nihms806693.pdf) | Clonal/heterogeneity response Results and figure; manuscript cover terms | No affirmative redistribution permission established for inspected copy; metadata only. Its joint cutoff groups are also not finite calibrated functional-dose bins for I1. Earlier literature is already known to the merged note. |
| [Rosenthal 2019](https://doi.org/10.1038/s41586-019-1032-7); [primary manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6954100/) | Early untreated lung cancer; neoantigen expression/loss/restriction | Mechanistic context, not calibrated checkpoint response versus escape hazard. No necessary numeric extraction or established affirmative redistribution basis for inspected copy. |
| [Montesion 2021](https://doi.org/10.1158/2159-8290.CD-20-0672) | HLA-LOH and checkpoint survival/response; publisher indexed Results | Adverse evidence retained qualitatively. Direct full access/affirmative reuse unresolved; no numerical extraction. Some source HRs are intact-versus-LOH, so blindly importing their sign would invert the comparison. Second-line monotherapy differs from the Lozac'hmeur cohort. |
| [Samstein 2019](https://doi.org/10.1038/s41588-018-0312-8) | Total TMB, histology-dependent percentile cutoffs and survival | Not functional clonal load, calibrated founder count or escape/response bridge; no numerical extraction. |
| [Valero 2021](https://doi.org/10.1001/jamaoncol.2020.7684); [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7893543/) | Total-TMB response thresholds, nested cutoffs | Article copyright is All Rights Reserved; metadata only. Nested/open-ended thresholds must not become disjoint dose bins or arbitrary midpoints. |
| [Morris 2023 online / 2024](https://doi.org/10.1002/lary.30761); [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10687312/) | HLA LOH, clinical benefit; unadjusted and adjusted models | Opposite unadjusted association becomes imprecise after adjustment. No affirmative reuse license found; no numbers copied. Different benefit endpoint, not I1 escape. |
| [Mismatch-repair deficiency / immunogenicity study, 2023](https://doi.org/10.1038/s41588-023-01499-4) | Primary Results on mixed mouse clones, antigen-bearing cell mass and T-cell responses; encountered during final source recheck | Relevant mechanistic warning about relative clonality versus absolute antigen-bearing cell mass. Mixtures include nonimmunogenic clones and do not provide I1's positive minimum functional-load/rate/escape calibration or a matched I2-H checkpoint-response contrast. No additional numeric extraction or after-the-fact formal test. |
| [Scanu, Werner & Huang, preprint](https://doi.org/10.1101/2024.10.22.619675); [version inspected](https://pmc.ncbi.nlm.nih.gov/articles/PMC12330626/) | Model/simulation of ecDNA maintenance and de-mixing | Explicitly a preprint; no new independent calibrated generation-two moments or killing test. Simulations under a binomial assumption cannot independently verify it. |
| [Hung 2024](https://doi.org/10.1038/s41586-024-07861-8) | Coordinated inheritance of different ecDNA species | A joint cross-species correlation does not by itself reject the single-species marginal binomial law. No required conditional moment/rate aggregates. |
| [Yi et al., PMID34819316](https://pubmed.ncbi.nlm.nih.gov/34819316/) | Live-cell ecTag inheritance abstract | Uneven inheritance is compatible with stochastic binomial segregation; no matched founder/generation/rate moment test located. |
| [Sugimoto 2026](https://doi.org/10.1093/nar/gkag005); [correction](https://doi.org/10.1093/nar/gkag498) | ecDNA engineering, fitted damage survival and copy changes; primary indexed content | Fitted survival law is not independent rate calibration for E2. Publisher read failed; no numerical values extracted or correction-dependent claim made. |
| [MMEJ ecDNA study, 2026](https://doi.org/10.1038/s41586-026-11048-8) | Primary indexed growth/proliferation/ecDNA-loss results | No required conditional moments or independent threshold calibration located in screening; no numbers extracted. |
| [Shaffer 2017](https://doi.org/10.1038/nature22794) | Fluctuation Methods and holiday example | Baseline fluctuation plus exposure/recovery example, not matched growth-controlled holiday quartet with count-moment uncertainty. |
| [Russo 2022](https://doi.org/10.1038/s41588-022-01105-z) | Two-step fluctuation and early/late resistance | Continuous exposure, growth/mutation fitting; not holiday retention moments. |
| [Whiting 2025](https://doi.org/10.1038/s41467-025-59479-7) | Lineage/barcode dynamics and recovery-period confounding | Growth-dependent resensitization, no matched R1 moment quartet. Fitted parameter distributions are not new independent data. |
| [Pisco 2013](https://doi.org/10.1038/ncomms3467) | MDR phenotypic transitions | Does not supply the specified holiday count laws or conditioned covariance. |
| [Roux 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6848270/) | TRAIL repeat-exposure/recovery experiment | Recovery and signaling measurements, no matched four-arm count-moment law. |
| [Treatment scheduling 2021](https://doi.org/10.1038/s41523-021-00270-4) | Heterogeneous populations, holiday schedules | Ongoing growth, changed exposure and unequal holidays/viability endpoints; not R1's design. |
| [H3K4me3/OGT resistance study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10719180/) | Holiday after differing pretreatments | Changed pretreatment, not equal holiday steps after the same exposure. |

## Verification limits

Primary-source transcriptions and reuse checks were performed for every
retained row. Independent arithmetic checks used alternate hypergeometric
weighting/tail inversion, rational binomial sums and direct algebra before
push. These are distinct tasks: the unit tests recompute arithmetic and
report synchronization; they do not authenticate source truth, correct the
source papers or certify a complete literature search. Agreement among AI
agents is not biological evidence. No new treatment, diagnostic or
surveillance recommendation follows from this audit.

## Independent PR #63 pre-merge random-row check (2026-11-17)

A reviewer parsed all five CSVs and checked **all 190 scalar records** for a
unique ID, numeric aggregate field, study/population, locator, source DOI/URL,
access date, source-specific reuse statement and interpretation notes. No
patient-level fields or uncontrolled-access records were present. All 12
named primary sources have an affirmative article-specific Creative Commons
reuse basis recorded, including the NC-ND restrictions on unchanged facts
where applicable; third-party source materials are not relicensed by the
repository. This is a provenance/schema and rights *review*, not a mechanical
proof of a publisher's license or independent replication of every paper.

The random draw was performed *without replacement* across the 190 CSV data
rows using Fisher–Yates and mulberry32 with seed 63, before selecting any
source for verification. The first **20** rows were checked. Ten values were
rechecked directly against accessible primary-source text, figure labels or
PDF pages:

| Drawn CSV row | Source locator | Published scalar checked |
| --- | --- | ---: |
| `henriksen22_draw_interval` | Henriksen 2022, study-design section, printed p508 | 3 months |
| `yang_2022_melanoma_pd1_hla_b_response_loh` | Yang 2022 Fig. 1d, B Loss, R | 4 |
| `yang_2022_urothelial_hla_c_nonresponse_loh` | Yang 2022 Fig. 1c, C Loss, NR | 30 |
| `henriksen24_gain_q1` | Henriksen 2024, printed p233 sampling-frequency Results | 1 month |
| `lozac_luad_p` | Lozac'hmeur 2024, Fig. 3c caption | 0.032 |
| `yang_2022_melanoma_ctla4_hla_c_nonresponse_loh` | Yang 2022 Fig. 1e, C Loss, NR | 12 |
| `shibata_2026_clone9_ks_p` | Shibata 2026, Fig. 3d, Clone 9 | 0.98 |
| `yang_2022_melanoma_ctla4_hla_c_response_loh` | Yang 2022 Fig. 1e, C Loss, R | 4 |
| `shibata_2026_clone2_ks_p` | Shibata 2026, Fig. 3d, Clone 2 | 0.96 |
| `henriksen24_gain_q3` | Henriksen 2024, printed p233 sampling-frequency Results | 2 months |

The other **ten** selected rows are
`lange_2022_tr14_cdk4_n`,
`shibata_2026_ctrlen2_jq1_p`,
`lange_2022_gbm39_egfr_ks_p`,
`lange_2022_snu16_fgfr2_ks_p`,
`lange_2022_colo_myc_n`,
`shibata_2026_clone1_rho`,
`shibata_2026_ctrl1_dmso_rho`,
`shibata_2026_clone11_p`,
`lange_2022_tr14_mycn_ks_p`, and
`shibata_2026_clone10_p`.
Their exact cell/table locators and values remain in the CSVs, and the earlier
author's original source-cell audit is documented above; their Springer
supplementary XLSX/PDF content could not be independently fetched during this
pre-merge recheck. **Do not interpret the 10/20 directly checked as full
independent source verification of all 190 rows.** No scraping workaround,
unlicensed reproduction or imputation was used. This limitation is documented
rather than silently treating inaccessible supplements as inspected.

The review also separately compared the protocol's three verdict definitions
to the output and corrected I3's overstated non-rejection label as recorded in
AMENDMENTS.md. The exploratory E3-X result remains explicitly conditional on
mitosis independence and calibrated tests, and cannot certify any clustered
segregation extension or dose threshold.
