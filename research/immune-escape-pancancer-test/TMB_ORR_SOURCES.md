# TMB and objective-response source audit

**Nonclinical, unreviewed computational research. No treatment advice.**

Access date throughout this audit and its source CSV is **2026-11-17,
America/New_York**. Sources were retrieved after the published protocol-only
commit `5ce3b2261a05a9763f79a7942b3f7fe30e108033`. Source choice and category
mapping preceded fitting. This is a targeted source audit, not a systematic
review of every checkpoint-inhibitor trial.

## Selected sources and affirmative reuse basis

The objective-response panel is the complete aggregate **Table S1** of
Mo, Cai, Kuai, Li and Chen, *Universal cutoff for tumor mutational burden in
predicting the efficacy of anti-PD-(L)1 therapy for advanced cancers*,
Frontiers in Cell and Developmental Biology (2023),
[DOI: 10.3389/fcell.2023.1209243](https://doi.org/10.3389/fcell.2023.1209243).
The [publisher article](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2023.1209243/full)
has an explicit **CC BY 4.0** copyright notice attributing Mo and colleagues.
No separate exclusion from that notice was found on the aggregate supplement.
The PDF `DataSheet1.PDF`, printed page **5/6**, Table S1, was retrieved from the
official [Europe PMC supplementary-files endpoint](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10248461/supplementaryFiles).
The corresponding deposit's [full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10248461/fullTextXML)
identifies the supplement and reproduces the article's license notice.
The table contains only cancer-type summary counts and response rates.

All published rows of that table are retained in `data/tmb_orr.csv`, including
zero responses and biomarker-defined strata. The source's methods identify
anti-PD-1/PD-L1 monotherapy and attribute the response summaries to Yarchoan.
Its data-availability statement explicitly identifies
**10.1056/NEJMc1713444**. The methods contain a bibliographic inconsistency:
the text names Yarchoan but the adjacent parenthetical citation points to
Wang (2023). The explicit data-availability DOI establishes the stated lineage;
it does not prove that the later authors exactly reproduced every original
trial-pooling decision. We therefore identify the numeric source as **Mo's
published table**, not as our reconstruction of the original Yarchoan panel.
This is a reuse of earlier outcome evidence, not an independent response cohort.

The TMB values are the numerical medians in **Additional file 3, Table S1** of
Chalmers et al., *Analysis of 100,000 human cancer genomes reveals the landscape
of tumor mutational burden*, Genome Medicine (2017),
[DOI: 10.1186/s13073-017-0424-2](https://doi.org/10.1186/s13073-017-0424-2).
The [publisher article](https://link.springer.com/article/10.1186/s13073-017-0424-2)
expressly licenses the article **CC BY 4.0** and applies the **CC0 data waiver**
to published data unless otherwise stated. No contrary notice was found in
the selected aggregate file. The exact workbook is
[`13073_2017_424_MOESM3_ESM.xlsx`](https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-017-0424-2/MediaObjects/13073_2017_424_MOESM3_ESM.xlsx),
sheet `Sheet1`; the CSV records the original disease label, specimen count,
median, and exact cell locator for each selected disease.
The article's description of this particular file explicitly says it contains
summary TMB properties by disease. Its contents were inspected before use;
it has no patient records. We did not access the underlying clinical specimens
or individual genomic records. The article's statement that its underlying
patient data are unavailable for privacy reasons remains applicable to those
underlying records; it does not turn the separately published, licensed
disease-level table into patient data.

These licenses remain attached to the third-party numbers. Repository MIT
licensing does not replace them. The relevant terms are
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) and
[CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).
Changes consist of selecting the relevant published aggregate rows, normalizing
category labels, and converting ORR percentages to probabilities. TMB medians
and counts are transcribed without recalculation or synthetic pooling. Each
CSV row carries separate source title, DOI, article URL, data URL, locator,
access date, reuse basis, terms URL, definition, population, and transformation
fields for TMB and ORR. Missing TMB rows retain metadata and the reason for
missingness, with an empty numerical field.

## Source priority and excluded alternatives

| Source inspected | Decision before fitting | Reason and locator |
| --- | --- | --- |
| Yarchoan, Hopkins and Jaffee (2017), *Tumor Mutational Burden and Response Rate to PD-1 Inhibition*, DOI `10.1056/NEJMc1713444` | Bibliographic lineage and methods context; no direct numerical panel import | [NEJM primary record](https://www.nejm.org/doi/abs/10.1056/NEJMc1713444) states copyright and directs permission requests to NEJM; the [PMC author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6549688/) supplies access but no affirmative CC reuse license was established. The licensed Mo republication supplies the response facts used here. |
| Yarchoan et al. (2019), *PD-L1 expression and tumor mutational burden are independent biomarkers in most cancers*, DOI `10.1172/jci.insight.126908` | Bibliographic/methods context only | The [article](https://insight.jci.org/articles/view/126908) and its aggregate supplemental-table descriptions were inspected. Article copyright belongs to ASCI. [JCI terms](https://token.jci.org/kiosks/terms/) distinguish articles before August 20, 2020 from later CC BY articles; current site-wide open-access language does not retrospectively establish reuse rights for this article. |
| Lee and Ruppin (2019), *Multiomics Prediction of Response Rates to Therapies to Inhibit Programmed Cell Death 1 and Programmed Cell Death 1 Ligand 1*, DOI `10.1001/jamaoncol.2019.2311` | Prior ecological multi-predictor and cross-validation work; no numerical import | The [publisher record](https://jamanetwork.com/journals/jamaoncology/fullarticle/2748890) and [primary PDF](https://jamanetwork.com/journals/jamaoncology/articlepdf/2748890/jamaoncology_lee_2019_br_190007.pdf) identify prior mutation/immune/checkpoint prediction. The PDF says copyright American Medical Association, all rights reserved. NIH author affiliations were not treated as permission to redistribute the table. |
| Mo et al. (2023), DOI `10.3389/fcell.2023.1209243` | Selected complete ORR panel | Explicit license, cancer-type-only numerical Table S1, monotherapy endpoint, and stated Yarchoan lineage. The paper's TMB analysis uses a different cohort and its own TMB medians are displayed as boxplot markers in Figure 2A rather than printed as numerical medians in Table S1. The protocol's numerical-table priority favors Chalmers's printed medians; no patient data were obtained to calculate Mo-cohort medians. |
| Zhang et al. (2025), *The MHC-I-dependent neoantigen presentation pathway predicts response rate to PD-1/PD-L1 blockade*, DOI `10.17305/bb.2024.11069` | Relevant prior work; not selected for the response panel | [Primary article](https://www.bjbms.org/ojs/index.php/bjbms/article/view/11069), [official deposited XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12042685/fullTextXML), and [supplement PDF](https://www.bjbms.org/ojs/index.php/bjbms/article/download/11069/3629/91031) have CC BY material. Table S1 has response summaries but uses a rectal category whose relation to the original colorectal MSI/MMR strata is unresolved, and its sample-size column is not established as response-trial denominators. Table S2 reports correlations, not the missing per-cancer feature panel. These endpoint/category ambiguities disfavor it relative to Mo before fitting. Publisher access intermittently failed; official deposited XML and the aggregate supplement were readable. |
| Lahman et al. (2021), *Quality Is King: Fundamental Insights into Tumor Antigenicity from Virus-Associated Merkel Cell Carcinoma*, DOI `10.1016/j.jid.2020.12.037` | Mechanistic/novelty context only | [PubMed](https://pubmed.ncbi.nlm.nih.gov/33863500/) states all rights reserved; the [publisher record](https://www.sciencedirect.com/science/article/pii/S0022202X21009994) identifies an Elsevier user license/open archive. Free access was not treated as an affirmative blanket redistribution license. This is a review, not an independent aggregate monotherapy response panel. |

No alternative was chosen or rejected by its fitted correlation, residuals,
or apparent agreement with the immune-escape prediction. All retrieval problems
were addressed through publisher or official archival copies, without bypassing
authentication or access controls. No source PDF, workbook, individual record,
or extracted individual datum is committed to the repository.

## Measurement and population definitions

**Response.** The response endpoint is the published cancer-type aggregate
objective response to anti-PD-1/PD-L1 monotherapy. The response table preserves
its own reported trial sample sizes; these are not the same as Mo's genomic
cohort sizes. The source's trial denominators are not enough to reconstruct
all trial-specific response definitions, weights, treatment-line distributions,
or selection rules. No responder counts are reverse-engineered from rounded
percentages. In particular, no binomial likelihood or patient-level confidence
interval is justified by treating a pooled response estimate as one homogeneous
binomial sample.

MSI/MMR-defined rows are explicitly selected by that biomarker. For other rows,
row-specific selection status is recorded as `unknown`, not as a claim of
an unselected population. The original Yarchoan methods describe selection
rules, but Mo's table does not establish complete row-specific PD-L1, viral,
receptor, or treatment-line metadata. An analysis excluding only known selected
rows cannot eliminate unknown selection in the retained types.

**TMB.** Chalmers defines TMB using coding somatic base substitutions and indels
per sequenced coding megabase. Synonymous substitutions are initially counted;
known driver alterations and likely germline variants are filtered as described
in the original methods. This assay-derived count is not a measured clonal
neoantigen burden. It does not identify antigen expression, binding, or
presentation.

The clinical profiling cohort contains mostly substantially pretreated,
advanced/metastatic cancers. The source methods state that **all FoundationOne
assay versions** were included. The frequently cited approximately 1.1 Mb
coding territory is the validation panel description, not an assurance of one
identical assay version for every specimen contributing to every median.
The CSV retains that distinction. The response and TMB values come from
different patient populations. No matching at patient level is claimed.

## Category mapping fixed before fitting

The mapping follows named disease/NOS categories where available. A median of
medians is never used to fabricate a pooled population median. An exact label
mapping says only that the published disease labels match; it does not establish
that a NOS subgroup represents the same histological mixture, stage, treatment
history, or population as the response trials. These cross-cohort limits also
apply to rows classified `exact`.

The CSV's `mapping_notes` enumerate every source category and substitution.
In particular:

* Urothelial response is paired with bladder urothelial TMB; cervical and anal
  labels use explicitly flagged squamous histology substitutes; head/neck uses
  the explicitly flagged HNSCC median.
* Ovarian response uses the epithelial NOS median with an approximate flag,
  not a silently assumed high-grade-serous or all-histology median. Endometrial
  response uses an explicitly flagged adenocarcinoma NOS category.
* Renal-cell and pancreatic categories use the published renal-cell NOS and
  pancreatic-carcinoma NOS medians. In particular, pancreatic TMB is not changed
  to the separate ductal-adenocarcinoma value to better match a hypothesized
  response pattern.
* Mesothelioma, prostate, sarcoma, and germ-cell categories use explicitly
  flagged pleural, acinar-adenocarcinoma, soft-tissue NOS, and testicular
  non-seminoma substitutes, respectively.
* The NSCLC response row uses the NSCLC NOS median. It is not split into
  synthetic adenocarcinoma and squamous response outcomes. NSCLC and SCLC have
  the same held-out parent group, `lung`.
* MSI and MSS colorectal response rows retain missing TMB because Chalmers's
  broad colon/rectum rows do not identify those strata. The heterogeneous
  non-colorectal MSI group also has no matching TMB. The pooled esophagogastric
  response has no corresponding pooled median; separate stomach, esophageal,
  and gastroesophageal-junction medians are not averaged. All these rows remain
  in the source CSV.

`tcga_code` is a later biomarker-joining anchor, not a claim of identical
histology. For example, `KIRC` beside renal-cell NOS or `OV` beside broad ovarian
response does not make either mapping exact for a clear-cell or high-grade-serous
clonality estimate. Those subsequent predictor mappings require their own
explicit approximation flags. A blank TCGA code for NSCLC prevents assigning
one lung histology's biomarker to the mixed outcome by accident.

## Numerical extraction audit and uncertainty

The ORR supplement was read as text and its original table page was rendered
and visually inspected. Every response percentage and every source trial count
in the CSV was checked against that page. A second extraction from the PDF's
text layer was then compared programmatically with the saved CSV. Each TMB
median, specimen count, and original category label was read directly from
the licensed XLSX's XML cells, and then re-read and compared against the CSV
using its recorded cell locator. This confirms faithful transcription from
these sources; it does not independently certify the original publications'
upstream calculations or biological claims.

The conversion is `orr = original percentage / 100`. The CSV preserves the
original displayed percentage, including `0%`. A percentage shown with one
decimal is assigned a rounding half-width of `0.0005` in probability units.
The conservatively wider half-width for a percentage printed simply as `0%`
is `0.005`; perturbations must be clipped to the probability interval. TMB
medians have a `0.05 mutations/Mb` rounding half-width. These are reading/rounding
envelopes, not sampling confidence intervals. There is no figure digitization
in the final ORR/TMB input.

The audit deliberately leaves all missing categories and zero outcomes visible.
Automated tests, hashes, source agreement, and AI critique can check bookkeeping
and computational reproducibility; they do not validate immune-escape biology,
remove ecological confounding, or justify clinical use.
