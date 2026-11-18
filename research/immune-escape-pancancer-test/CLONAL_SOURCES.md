# Clonal-fraction source and transcription audit

**Nonclinical, unreviewed research. No treatment advice.** Access date:
2026-11-17 (America/New_York). No fitted response result was inspected during
selection or extraction of this predictor panel.

## Selected cancer-type aggregate panel

[Charoentong et al., Cell Reports (2017)](https://doi.org/10.1016/j.celrep.2016.12.019),
*Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype
Relationships and Predictors of Response to Checkpoint Blockade*, supplies
the [author-hosted supplementary tables](https://icbi.i-med.ac.at/TCIA/SupplementaryTables.pdf).
Table S3 is printed as **Table 3** and contains separate Driver, Passenger,
Clonal and Subclonal columns. The CSV transcribes only the last two numerical
columns and original cancer-category labels, with our own definitions and
provenance fields. It does not redistribute the source table, paper, figure
or database. The clonal column is the predictor; the complementary subclonal
column is retained to make a swapped-column transcription detectable.

This is a **published cancer-type fraction of predicted neoantigens assigned
clonal origin**, not a per-patient median clonal fraction and not a functional
antigen count. The retrieved table does not establish a per-patient weighting
rule or supply a variance; we do not invent either. The primary article says
clonality was inferred from cancer cell fractions and distinguishes clonal
and subclonal neoantigen origin. Single-sample inference can overstate tumour
wide clonality. The fraction is not an estimate of the harmonic per-cell
functional antigen load in PR59.

The table covers nineteen displayed cancer categories because the source
combines colon and rectal cancers into CRC. Seventeen have numeric fractions;
ovarian (OV) and stomach (STAD) rows explicitly print NA. Those missing rows
are preserved with blank numerical fields and original value NA. A TCGA-wide
fraction is not a biomarker-stratified fraction. In particular, BRCA is not
specific to triple-negative breast cancer, CRC is not MSI-stratified, KIRC is
clear-cell renal carcinoma rather than all renal histologies, and SKCM is
cutaneous melanoma rather than uveal melanoma. Joining these predictors to
broader/different response categories requires a visible mapping flag.

## Access and rights: an important source-quality limit

The supplementary PDF and the institutional thesis described below returned
HTTP 404 when directly opened during this run. The [live TCIA site](https://tcia.at/)
reports a maintenance pause. Both search systems nevertheless returned the
complete machine-indexed text of the primary supplementary table, including
every cancer code, column heading and value. All numeric CSV cells were
manually rechecked against this primary-source indexed text on repeat
retrieval; the clonal/subclonal columns add to a complete percentage in every
nonmissing source row. This is **indexed-primary-text verification**, not a
claim that the live PDF was retrieved, rendered or visually checked. The
source could not be independently version-pinned by downloading its bytes.
Every CSV row carries this limitation. Successful numerical tests cannot
remove it. A source-quality sensitivity excluding indexed-only sources
therefore has no clonal predictor panel; it must be reported as unavailable,
not as a negative coefficient result.

The article's original first-page license statement was retrieved in the
indexed text of coauthor [Mirjana Efremova's institutional thesis](https://icbi.i-med.ac.at/theses/Efremova2017.pdf),
which embeds the published article with its title, authors and DOI. It states
[CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/).
The [author TCIA About page](https://www.tcia.at/about), also indexed but
currently returning 404, identifies the paper and links its manuscript,
supplementary figures and supplementary tables. Bibliographic rights metadata
from the [National Diet Library](https://ndlsearch.ndl.go.jp/books/R100000136-I1361981470649694336)
independently lists the same license; it was a rights cross-check, not the
scientific data source. A preprint under a different license exists, but this
analysis does not silently apply the preprint license to the journal version.

Reuse here is limited to unchanged, isolated numerical facts in a new,
attributed noncommercial research compilation, with percentage-to-fraction
unit conversion. No expressive article/table/figure content is redistributed
or adapted. The source's license and attribution remain attached; the
repository's MIT license does not relicense these third-party observations.
The half-width field is half the last displayed digit after unit conversion,
conditional on ordinary rounding. It is not sampling uncertainty, a standard
error, a confidence interval or protection against a source transcription
error.

## Named-anchor audit and why other studies were not substituted

| Source | What was inspected | Decision before fitting |
| --- | --- | --- |
| [McGranahan et al., Science (2016)](https://doi.org/10.1126/science.aaf1490); [primary indexed manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC4984254/) | Main-text neoantigen/clonality descriptions and Figure 1/4 descriptions. Multiregion early lung cases, TCGA lung histologies and treatment-response association are relevant context. | Individual tumour bars cannot be read as cancer-type medians. A small pooled NSCLC heterogeneity summary is not a coherent pan-cancer panel. The inspected manuscript did not establish affirmative redistribution terms; no new numerical rows copied. No patient-level supplement obtained. |
| [Litchfield et al., Cell (2021)](https://doi.org/10.1016/j.cell.2021.01.002); [primary text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7933824/) | Summary, landscape of biomarkers and biomarker-analysis Methods; article-specific CC BY license in primary indexed text. | Clonal TMB response coefficients are not cancer-type clonal neoantigen fractions. No eligible printed pan-cancer fraction panel was located in the inspected text. No fraction derived from individual records or from a ratio of unrelated medians. |
| [Rosenthal et al., Nature (2019), TRACERx](https://doi.org/10.1038/s41586-019-1032-7); [primary text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6954100/) | Abstract, neoantigen depletion/expression Results and main/extended figure descriptions. | Lung-specific evolution and expression findings contextualize measurement failure. The fraction of clonal neoantigens expressed in every region has a different denominator from the fraction of all neoantigens that are clonal. It is not substituted for f. No controlled-access sequences or individual tumour bars used. |
| [Charoentong et al. (2017)](https://doi.org/10.1016/j.celrep.2016.12.019) | Primary indexed Supplementary Table S3, original article text in an author-hosted thesis, author TCIA metadata and license evidence. | Selected as the coherent published aggregate panel. All displayed categories retained, including missing values. Access/version limitations are explicit above and on every CSV row. |

## Search coverage and stopping

Search families included `McGranahan 2016 clonal neoantigens median fraction
cancer types`, `pan cancer clonal neoantigen fraction median`, `Charoentong
Table S3 fractions`, the exact table heading, `Litchfield 2021 clonal fraction
cancer`, and `TRACERx clonal neoantigen median fraction`. Primary article,
supplement, DOI, author-hosted manuscript and institutional copies were
followed. The selected table was found by its subject and coverage, before
response fitting. Publisher and institutional retrieval failures were
retained; no paywall, authentication, robots or data access control was
bypassed. Search-engine text retrieval is distinguished from successful
direct access. We stopped after the named anchors and a coherent aggregate
panel had documented outcomes, with a final check for the panel's exact
values, estimand, rights and missingness. This targeted search cannot establish
that no better published clonal-fraction panel exists elsewhere.
