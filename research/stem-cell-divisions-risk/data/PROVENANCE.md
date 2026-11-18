# Aggregate-data provenance and source audit

**Nonclinical research; no medical advice.** These files contain published
tissue-level aggregates and elementary calculations from those aggregates.
They contain no patient records, identifiers, individual genomes, or
controlled-access data. The calculations are unreviewed research outputs.

All sources below were accessed on **2026-11-17, America/New_York**. The
corresponding workspace clock was already 2026-11-18 UTC. Every CSV row
records its source, locator, access date, and reuse basis. The repository's
MIT license does not relicense third-party data.

## 1. Inventory and separation of versions

| CSV | Meaning | Permitted numeric source or derivation |
| --- | --- | --- |
| `tv2015.csv` | Original Tomasetti–Vogelstein 2015 Table S1 values for the five requested fields and 31 cancer entries | Little et al. 2016 CC0 republication, with the three original developmental-only zero-rate conventions restored as explained below |
| `little2016.csv` | Literal later Little et al. 2016 values for those same entries and fields | Unmodified factual columns of the CC0 S2 dataset |
| `tv2015_nunney2020.csv` | Our derived application of the developmental-division correction to the original 2015 inputs | Own arithmetic on the permitted original-input reconstruction; mathematical correction attributed to Nunney and Thai 2020 |
| `liu2022_ranked_risk_audit.csv` | A later endpoint audit: 17 cancer-type risk minima, medians, and maxima | Direct cells of Liu et al. 2022 Additional file 2, Table S1; source data covered by the publication's CC0 statement unless separately credited |

The first three files reuse the same risk observations. They are sensitivity
versions, not three independent datasets, and must not be pooled as
independent observations. The fourth file does not supply a matching
stem-cell or division-count vector and is not used as a cross-tissue
risk-versus-divisions fit.

## 2. Original 2015 inputs: source, rights, and exact locators

Original study: Tomasetti C, Vogelstein B. *Variation in cancer risk among
tissues can be explained by the number of stem cell divisions.* Science
2015;347:78–81. [DOI: 10.1126/science.1260825](https://doi.org/10.1126/science.1260825).
The relevant original locator is **Supplementary Table S1, pages 13–15**.

Redistribution relies on a separately licensed republication:
Little MP, Hendry JH, Puskin JS. *Lack of Correlation between Stem-Cell
Proliferation and Radiation- or Smoking-Associated Cancer Risk.* PLOS ONE
2016;11:e0150335.

- [Publication and copyright/Data Availability paragraphs](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0150335).
- [DOI: 10.1371/journal.pone.0150335](https://doi.org/10.1371/journal.pone.0150335).
- [S2 Text downloadable supporting archive](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0150335.s002&type=supplementary).
- Exact archive member: `Tomasetti & Vogelstein Science augmented full.csv`.
- Reuse instrument: [Creative Commons CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/), expressly identified in the publication's copyright paragraph. The Data Availability statement expressly identifies S2 Text as containing the datasets.

Only the small factual input columns and cancer names are reproduced.
The source's radiation estimates, smoking estimates, extra-risk scores,
logged columns, and missing-value sentinels are not imported.

| Repository column | Original Table S1 position | Little S2 exact field |
| --- | --- | --- |
| `risk` | Column 2, lifetime cancer incidence | `Lifetime_cancer_incidence` |
| `stem_cells` | Column 4, stem-cell number (s) | `Number_stem_cells` |
| `divisions_per_year` | Column 5, annual division rate | `Divisions_stem_cell_per_yr` |
| `divisions_per_lifetime` | Column 6, lifetime post-growth divisions (d) | `Divisions_per_stem_cell_lifetime` |
| `total_divisions` | Column 7, reported lifetime total (lscd) | `Divisions_all_stem_cell_lifetime` |

Source order is preserved. The row locators in `tv2015.csv` and
`little2016.csv` count data rows after the header. The corrected file's
locators explicitly count the header where stated. All files also carry a
stable `tissue_id` and source cancer name.

### Independent check against the original table

[Harvard BioNumbers record 111704, version 10](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111704&ver=10)
links an image reproduction of the complete three-page table:

[Original Table S1 scan](https://bionumbers.hms.harvard.edu/files/Lifetime%20cancer%20incidence%20and%20parameters%20related%20to%20the%20normal%20stem%20cells%20that%20are%20precursors%20of%20these%20cancers.pdf).

All three pages were rendered and inspected. All **155 numeric fields**
(31 entries times five fields) were independently transcribed into a
temporary audit file and compared using decimal arithmetic with
`tv2015.csv`; all matched. The original scan provides a check on source
faithfulness, not permission to redistribute an AAAS image or PDF. Its PDF
creation date is not evidence of a Science correction.

The literal later file was also compared directly with the downloaded
CC0 source: all 155 numeric fields in `little2016.csv` matched its five
source columns. Source lexical numerical strings are retained where
possible; scientific notation and decimal notation are compared by value.

### Merge-review independent source spot check (2026-11-18)

The following **12 original TV2015 rows** were checked afresh against the
three-page **printed Table S1 scan linked in section 2**, including all five
numeric fields per row (risk, N, annual divisions, lifetime divisions, D).
The table below reproduces only a short identifier, the risk and N as
cross-check identifiers; the other three fields are checked against the
committed CSV, not recopied here. This is an independent source *spot
check*, not a second 155-field whole-table transcription.

| Original Table S1 row | Risk | N | Check |
| --- | ---: | ---: | --- |
| Acute myeloid leukemia | 0.0041 | 1.35e8 | matched |
| Basal cell carcinoma | 0.3 | 5.82e9 | matched |
| Chronic lymphocytic leukemia | 0.0052 | 1.35e8 | matched |
| Colorectal adenocarcinoma | 0.048 | 2e8 | matched |
| Colorectal adenocarcinoma with FAP | 1 | 2e8 | matched |
| Duodenum adenocarcinoma | 0.0003 | 4e6 | matched |
| Esophageal squamous cell carcinoma | 0.001938 | 8.64e5 | matched |
| Glioblastoma | 0.00219 | 1.35e8 | matched |
| Hepatocellular carcinoma | 0.0071 | 3.01e9 | matched |
| Lung adenocarcinoma (nonsmokers) | 0.0045 | 1.22e9 | matched |
| Ovarian germ cell | 0.000411 | 1.1e7 | matched |
| Thyroid papillary/follicular carcinoma | 0.01026 | 6.5e7 | matched |

Additional version checks on those same source rows: the Little2016
reanalysis retains the original risk, N and D, and the documented
glioblastoma/ovarian maintenance fields become annual 0.025 and lifetime 2;
this agrees with Little et al.'s stated division-by-N then division-by-80
method. The Nunney2020 *derived sensitivity* entries for AML, glioblastoma
and thyroid papillary/follicular recompute respectively to
129734999999, 134999999 and 519999999 by N*(L+1)-1. These three
values are **our arithmetic**, not values claimed as a published Nunney
table. The 2022 source's published supplementary workbook, its exact cell
locators and its data CC0 reuse statement were checked in the original
source audit. This merge review rechecked the publisher's workbook listing,
age endpoint and rights but could not independently inspect the binary XLSX
values; the Liu numeric entries therefore retain the original documented
transcription audit and are not mislabelled as newly verified here.

## 3. The three later division-clock changes

Original Table S1 and Little's republication agree on risk, stem-cell
number, and total divisions for all 31 entries. They differ in exactly
these six fields:

| Tissue | Original annual divisions | Little 2016 annual divisions | Original lifetime post-growth divisions | Little 2016 lifetime divisions |
| --- | ---: | ---: | ---: | ---: |
| Glioblastoma | 0 | 0.025 | 0 | 2 |
| Medulloblastoma | 0 | 0.025 | 0 | 2 |
| Ovarian germ cell | 0 | 0.025 | 0 | 2 |

The zero entries in the original are published values, not absent data.
In all three cases the permitted source inputs satisfy (D/N=2). The
original convention can therefore also be reconstructed from the CC0
columns as (L=D/N-2=0), with annual value (L/80=0), and checked against
the original scan. This calculation is stated in the affected rows'
`reuse_basis` fields. It is a reconstruction of the source's displayed
convention, not an inference of a biological lineage clock or an assertion
that the rounded source formula holds exactly in every row.

`little2016.csv` retains the later values without alteration. These are a
later processing/modeling convention, not a documented numerical erratum
to the Science paper.

The [AAAS Science in the Classroom transcription](https://www.scienceintheclassroom.org/research-papers/does-bad-luck-cause-cancer)
identifies the **23 January 2015 correction** as replacement of a figure
with overlapping labels and a change to the final sentence. It does not
describe a numerical Table S1 amendment. The targeted search did not
locate a formal numerical erratum. This audit does not claim that no such
record could exist elsewhere.

## 4. Developmental-division correction: a derived sensitivity file

Nunney L, Thai K. *Determining cancer risk: the evolutionary multistage
model or total stem cell divisions?* Proceedings of the Royal Society B
2020;287:20202291.
[DOI: 10.1098/rspb.2020.2291](https://doi.org/10.1098/rspb.2020.2291);
[PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7779499/).
The precise correction locator is **section 5, Material and methods,
final paragraph**.

The original bookkeeping formula is

\[
D_{\rm original}=N(L+2)-2.
\]

Nunney and Thai point out that binary growth from one cell to (N) cells
requires (N-1) cell divisions. Each binary division increases the cell
count by one, which proves this identity directly. Adding (NL)
post-growth divisions gives

\[
D_{\rm corrected}=NL+N-1=N(L+1)-1.
\]

`tv2015_nunney2020.csv` is **our own application of that mathematical
correction to the original 2015 (N,L) inputs**. It is not a table
published by Nunney and Thai, not a new biological measurement, and not a
formal Tomasetti–Vogelstein erratum. Only `total_divisions` changes.
All 31 results were checked by exact decimal arithmetic, and the other
124 numeric fields are unchanged from `tv2015.csv`.

The calculation uses the original (L=0) convention for the three
developmental-only rows. Substituting Little's later (L=2) into this
formula would mix conventions. Because the source total (D) is rounded,
the corrected value is computed anew from (N,L); it is not simply the
displayed source total minus (N). Final digits reflect arithmetic on
rounded source inputs, not measurement precision.

**Reuse basis:** the main Nunney article states copyright with all rights
reserved. Its separate Figshare supplement is CC BY 4.0, but that does not
relicense the main article. This file copies no article prose, figure,
table, or numerical input from that article. It contains our elementary
arithmetic on the permitted CC0 inputs and credits the published
mathematical correction. The rights basis does not depend on falsely
calling the main article CC BY.

## 5. Numerical and semantic audit traps

The following checks concern source fidelity, not fitted results:

- The leukemia entries are **acute myeloid leukemia (AML)** and **chronic lymphocytic leukemia (CLL)**. They are not ALL or CML.
- Esophageal stem-cell number is **864000**. Transposing digits changes a covariate.
- Thyroid total divisions are **585000000** for papillary/follicular cancer and **58500000** for medullary cancer. The two entries do not have the same total.
- Popular [unofficial gist transcriptions](https://gist.github.com/lmmx/4ab13755925e068d2999) include leukemia-label, esophageal-number, and thyroid-number discrepancies. No input is taken from that gist. Original Table S1 and the CC0 source agree on the fields retained here, apart from the three explicitly documented division-clock changes.
- Source column names beginning `Ln_` or `ln_` are not reliable evidence of a natural-log transformation. For example, the AML source logged risk matches a base-10 logarithm. All analysis transformations are recomputed from the unlogged factual fields.
- Annual division rate and lifetime maintenance divisions are separate source estimates. Do not impose (L=80b). The liver entries report annual rate **0.9125** and lifetime divisions **88**; both are retained as published.
- The reported total (D) is also rounded. Do not overwrite it merely to force exact equality with (N(L+2)-2).
- (D/N), post-growth lifetime divisions (L), and ancestral lineage depth (L+\log_2N) are different quantities. The CSV names distinguish the first two; neither is silently substituted for the other.
- Colorectal FAP risk is **1**, as published. The CSV retains it. An exact survival-log or inverse model requiring risk below one needs an explicit boundary or censoring treatment, not a fabricated replacement such as an arbitrary value just below one.
- No tissue-specific sampling denominators or confidence intervals are supplied by these input tables. Numeric display precision is not a confidence interval.

Six source-labelled high-risk subgroups are flagged before fitting:
colorectal FAP, colorectal Lynch syndrome, duodenal FAP, HPV16-associated
head and neck, HCV-associated liver, and smokers' lung adenocarcinoma.
The anatomical groups also keep related or overlapping entries together:
AML/CLL, basal-cell/melanoma, colorectal subgroups, duodenum/total small
intestine, glioblastoma/medulloblastoma, head-and-neck subgroups, liver
subgroups, lung subgroups, total/subsite osteosarcomas, both pancreas
entries, and both thyroid entries. Separate esophagus, gallbladder,
ovary, and testis groups complete the 15 groups. Grouping is an analysis
choice recorded in the protocol; it is not a claim that the original
observations are independent samples from 15 biological populations.

## 6. TV2017: what was and was not acquired

Tomasetti C, Li L, Vogelstein B. *Stem cell divisions, somatic mutations,
cancer etiology, and cancer prevention.* Science 2017;355:1330–1334.
[DOI: 10.1126/science.aaf9011](https://doi.org/10.1126/science.aaf9011);
[PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC5852673/).

The paper describes 17 cancer types, 69 countries, and 423 registries.
The supplement captions identify cumulative-risk tables for ages zero
to 75, 80, 85, and 85-plus (Tables S1–S4). Tables S5 and S6 concern the
separate mutation-attribution analysis. The country correlations and the
mutation-source calculation are distinct analyses; reproducing one would
not validate the other.

The following legitimate retrieval routes were checked:

| Source route | Observed result and limitation |
| --- | --- |
| [Science supplementary methods PDF](https://www.science.org/doi/suppl/10.1126/science.aaf9011/suppl_file/tomasetti.sm.pdf) | HTTP 403 during the audit |
| [PMC supplementary methods file](https://pmc.ncbi.nlm.nih.gov/articles/instance/5852673/bin/NIHMS940482-supplement-Suppl_info.pdf) | Download-preparation HTML rather than a PDF |
| [Europe PMC supplementary mirror](https://europepmc.org/articles/PMC5852673/bin/NIHMS940482-supplement-Suppl_info.pdf) | HTTP 403 during the audit |
| [Nunney–Thai Figshare supplement](https://doi.org/10.6084/m9.figshare.13336881.v1) | CC BY 4.0 metadata verified; contains registry lists and summary statistics, but no complete per-tissue input panel |
| [BioNumbers source records](https://bionumbers.hms.harvard.edu/aboutus.aspx) | Citation guidance and selected quoted 2017 parameters were accessible; an explicit database redistribution license and a coherent input panel were not established |

No access control was circumvented. PMC's availability of an author
manuscript does not itself grant permission to mirror an AAAS supplement.
The Nunney Figshare license was checked through
[item metadata](https://api.figshare.com/v2/articles/13336881), within
[collection 5230395](https://api.figshare.com/v2/collections/5230395/articles);
its [DOCX file](https://ndownloader.figshare.com/files/25690221) was inspected.

Nunney and Thai's methods say that their prostate and esophageal
estimates came from TV2017, with the breast stem-cell estimate also from
TV2017 and a separately chosen age-dependent breast division rate.
Thus silently attaching 2015 stem-cell inputs to a later risk table would
create an undocumented hybrid, including an incorrect assumption that
esophageal inputs were unchanged.

**Unresolved acquisition limitation:** this directory does not contain
a coherent licensed TV2017 per-tissue risk-plus-divisions input panel and
does not claim an independent numerical replication of that panel.

## 7. Liu 2022 endpoint audit: exact cells and restricted use

Liu X, Yang J, Li H, et al. *Quantifying substantial carcinogenesis of
genetic and environmental factors from measurement error in the number
of stem cell divisions.* BMC Cancer 2022;22:1194.
[DOI: 10.1186/s12885-022-10219-w](https://doi.org/10.1186/s12885-022-10219-w).

- [Publication, Data Availability, and Rights and permissions](https://link.springer.com/article/10.1186/s12885-022-10219-w).
- [Additional file 1: methods PDF](https://static-content.springer.com/esm/art%3A10.1186%2Fs12885-022-10219-w/MediaObjects/12885_2022_10219_MOESM1_ESM.pdf), **Supplementary Method 1, printed pages 2–3**.
- [Additional file 2: spreadsheet](https://static-content.springer.com/esm/art%3A10.1186%2Fs12885-022-10219-w/MediaObjects/12885_2022_10219_MOESM2_ESM.xlsx), **Table S1**.

The publication licenses article text under CC BY 4.0 and applies the
CC0 waiver to data unless otherwise stated in a credit line. No excluded
credit notice is present for the reproduced Table S1 aggregate cells.
Attribution is retained. These are published aggregate estimates, not
individual records.

The source contains 423 ascending risk values in each of 17 tissue
columns. The committed file copies only the following cells:

| CSV value | Source Table S1 cells | Meaning |
| --- | --- | --- |
| `tissue` | B2:R2 | Source cancer-type headers |
| `risk_min` | B3:R3 | Within-column rank 1 |
| `risk_median` | B214:R214 | Within-column rank 212, the median of 423 |
| `risk_max` | B425:R425 | Within-column rank 423 |

The 17 columns run from head-and-neck squamous cell carcinoma in B
through breast cancer in R. Each row's `source_locator` names its exact
column and the three exact cells. Source XML numeric strings are retained
without rounding away binary-floating-point representation artifacts.
The extraction checked that each full column was ascending and that the
rank labels at the three selected rows were the expected values.

**Each cancer column was independently sorted by the source authors.**
Rows are ranks and do not retain registry or country identity. The 17
minima, medians, or maxima are therefore not matched geographical
populations. No geographical correlation or causal environmental effect
is inferred from this excerpt.

Supplementary Method 1 specifies diagnoses from 2003–2007 in CI5-X,
423 registries in 68 countries, with cumulative risk calculated from
five-year incidence bands covering birth through age 74. This is the
interval `[0,75)` in age notation, not unrestricted lifetime risk. The
different country counts reported by Liu and TV2017 are left as each
source reports them; they are not silently reconciled.

All seven worksheets, shared strings, drawing text, and the methods
supplement were checked for a numeric stem-cell-number, annual-rate,
or total-divisions vector. None was supplied. The publication sends
readers to the TV2015 and TV2017 supplements for those inputs. Therefore
this CSV is an **endpoint/provenance audit only**, not a corrected
TV2017 input table and not a replication fit. Its presence does not
resolve the acquisition limitation above.

## 8. Repository contents and hashes

Publisher PDFs, table scans, images, downloaded archives, spreadsheets,
DOCX supplements, and temporary full-source transcriptions were kept
outside the repository in the session's `source_data/` working directory.
They are not committed. No paper, publisher collection, or full external dataset is mirrored here. Only the small aggregate CSVs and
original documentation/code are included.

The committed verifier computes hashes on content after converting
CRLF and bare CR to LF. The expected normalized SHA-256 values for these
data files at source-audit completion are:

| File | Normalized SHA-256 |
| --- | --- |
| `tv2015.csv` | `40101d378d61c098b5a7e4610861319114644de7f3b25931b1ec3d8380598704` |
| `little2016.csv` | `baee121cc5c08144177d2add71e777821d549f7ee542a7262d1744f1ba80a32d` |
| `tv2015_nunney2020.csv` | `d72da38db78c182aec6659f01b928f00f577e2ca6a2ddb262953a5a67af8e3c6` |
| `liu2022_ranked_risk_audit.csv` | `919a0e48e5455f38cb8d99b85b40f050a591888dd2061f38eb3e0eb36f4550fe` |

Hashes detect local alteration; they do not prove the truth of an
upstream biological estimate or replace the source and visual checks.
