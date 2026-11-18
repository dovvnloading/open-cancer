# HLA-LOH source and transcription audit

**Nonclinical, unreviewed ecological research. No treatment advice.**

Access date throughout: **2026-11-17, America/New_York**. This audit was
performed after the protocol-only commit and before any model fitting. The
CSV contains only published cancer-type aggregates. Source articles and a
published aggregate figure were inspected; no patient-level source data,
controlled-access data, individual points, source-data workbooks, or
sequencing records were downloaded to obtain these estimates.

## Source selection before fitting

The primary HLA panel is the coherent **Marty Pyke et al. 2022 DASH panel**,
chosen for coverage, identifiable any-HLA-I-LOH definitions and its affirmative
CC BY license. Every cancer-type frequency shown in its Figure 5a is retained.
The source displays only types with sufficient cohort size; omitted types
cannot be reconstructed from its overall sample size. The four percentages
printed in the Results take priority over figure measurements, as specified
by the frozen protocol.

Montesion et al. is the preferred anchor, but an affirmative redistribution
license for its original inspected article/figure was not established.
Lozac'hmeur et al. republishes a smaller, licensed aggregate subset. Its
Montesion column is retained as a separate context panel, with its original
cohort lineage. It is neither an independent cohort nor silently pooled with
DASH. McGranahan's TRACERx lung panel remains a separate context panel.

The CSV's `panel` and `primary` fields identify these roles. Four alternative
Figure 5a readings accompany the corresponding DASH printed-text values.
They are alternative source readings, not additional observations. Any
subsequent analysis combining panels must identify that departure from the
coherent-panel analysis and preserve provenance per mapped predictor.

## DASH: original coherent pan-cancer panel

Rachel Marty Pyke and colleagues. *A machine learning algorithm with
subclonal sensitivity reveals widespread pan-cancer human leukocyte antigen
loss of heterozygosity*. Nature Communications (2022).
[DOI 10.1038/s41467-022-29203-w](https://doi.org/10.1038/s41467-022-29203-w).
[Primary article](https://www.nature.com/articles/s41467-022-29203-w),
[public full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9005524/fullTextXML).

**Definition:** fraction of tumours with at least one loss event affecting an
HLA-A, HLA-B or HLA-C allele, called by DASH. This is an allele-agnostic
any-HLA-I statistic, not the fraction losing a named allele or HLA-A alone.
It does not imply complete absence of antigen presentation.

**Population:** a biobank series profiled using the ImmunoID NeXT platform,
predominantly stages II and III. The Methods subsection on pan-cancer analysis
excludes cases without HLA calls from frequency calculations. This series is
not an ICI-treated ORR cohort and differs from the advanced Foundation
Medicine population and from multi-region TRACERx sampling.

**Printed values:** Results subsection *HLA LOH is pervasive across tumor
types and exhibits signs of immune pressure on tumors*, first paragraph.
The HNSCC, NSCLC-A, NSCLC-SCC and liver percentages are copied as printed and
converted to fractions. Half-widths for these rows record whole-percentage
rounding only, not sampling uncertainty.

**Figure measurements:** Figure 5a, right horizontal aggregate bars. The
[original publisher PNG](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-022-29203-w/MediaObjects/41467_2022_29203_Fig5_HTML.png)
has dimensions 1500 by 2113 pixels. Pixel coordinates start at the upper-left
corner. Fraction-axis tick centres are approximately `(498, 0)`,
`(643, 0.2)`, `(789, 0.4)` and `(935, 0.6)`, where each pair is
`(horizontal pixel, fraction)`. Least-squares calibration gives
`pixel_x = 497.7 + 728.5 * fraction`.

The recorded `pixel_row_y` samples the interior of each grey aggregate bar
away from its dark horizontal confidence-interval line. `pixel_bar_right`
is the midpoint between its last filled pixel and the first background pixel.
Convert with `(pixel_bar_right - 497.7) / 728.5`, then round to two decimal
places. The conservative absolute reading half-width is **0.01**, covering
pixel resolution, anti-aliasing, rounding and visual identification. These
are extraction envelopes, not confidence intervals. No patient points were
measured.

The dark horizontal endpoints are separately retained in `source_ci95_low`
and `source_ci95_high`, with their own reading half-width and pixel
coordinates. The caption identifies these as **95% confidence intervals**
but does not specify their construction. They are source intervals for HLA
prevalence, not ORR intervals, and are not reverse-engineered into binomial
counts or assumed to be Wilson intervals. For printed-text primary rows,
these remain intervals read from the figure; the noted text/figure point
discrepancy is not resolved by attaching an interval. A lower endpoint
coincident with the zero axis is recorded as the plotted zero endpoint.

The digitized bars were checked visually against their labels and magnitudes.
The CSV retains pixel coordinates so a reviewer can repeat this step against
the same primary image. The four printed estimates also retain separate
figure-reading rows. **The Results' NSCLC-SCC percentage exceeds the value
read from its Figure 5a bar.** Both are preserved without choosing by model
fit. The liver text/figure difference is smaller and overlaps the conservative
reading envelope. The source does not explain these differences, and this
audit does not resolve them by inventing denominators.

**Denominators:** no per-type denominator is inferred from a count-bar length
or a rounded percentage. The published Figure 5c caption gives counts for a
different genome-wide-LOH comparison; those are not assumed to be Figure 5a
frequency denominators. Consequently DASH denominator fields remain blank.

**Reuse:** the article's Permissions element explicitly grants
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The original authors
and source are credited on every CSV row. Only isolated numerical aggregates
and the extraction method are included; the paper, figure artwork and raw
database are not redistributed. The repository MIT license does not relicense
these third-party facts or their source material.

## Montesion: licensed aggregate republication, separate context

Ariane Lozac'hmeur and colleagues. *Detecting HLA loss of heterozygosity within
a standard diagnostic sequencing workflow for prognostic and therapeutic
opportunities*. npj Precision Oncology (2024).
[DOI 10.1038/s41698-024-00665-z](https://doi.org/10.1038/s41698-024-00665-z).
[Primary article](https://www.nature.com/articles/s41698-024-00665-z),
[public full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11300791/fullTextXML).

Only the **Montesion et al.** columns of **Table 1, Prevalence of HLA class I
LOH across cancer indications and datasets**, are extracted. Each numeric
row and its published denominator was checked against the article's XML
table. Percentage divided by 100 produces the CSV fraction. Half-width
0.0005 records rounding to one decimal percentage point. No numerator is
reverse-engineered from that percentage and denominator.

The table explicitly identifies its **Tempus** estimates as **HLA-A-only**
in the footnote. Those values are not substituted for the different
Montesion any-HLA-I estimates. The original Montesion source is:

Meagan Montesion and colleagues. *Somatic HLA Class I Loss Is a Widespread
Mechanism of Immune Evasion Which Refines the Use of Tumor Mutational Burden
as a Biomarker of Checkpoint Inhibitor Response*. Cancer Discovery (2021).
[DOI 10.1158/2159-8290.CD-20-0672](https://doi.org/10.1158/2159-8290.CD-20-0672).
[Primary article](https://aacrjournals.org/cancerdiscovery/article/11/2/282/2772/Somatic-HLA-Class-I-Loss-Is-a-Widespread-Mechanism).

The original paper supplies the any-HLA-I definition and Foundation Medicine
cohort lineage, enriched for advanced disease. Its all-rights-reserved notice
was not treated as a permissive source-data license. The numeric data here
come from the explicitly licensed later article, not from an attempted
reconstruction of Montesion's underlying patient database. The republication
does not add independent biological replication.

**Reuse:** Lozac'hmeur's Rights and permissions section grants
[CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/).
Only unchanged isolated numerical facts are transcribed for noncommercial
research, with original authors and source lineage retained. No article,
table artwork or figure is adapted or redistributed. CSV format conversion
and our separate analysis do not purport to relicense the source. This
follows the same factual-excerpt distinction documented in PR63; see the
[Creative Commons FAQ on data and format changes](https://creativecommons.org/faq/).

## McGranahan: lung-only TRACERx context

Nicholas McGranahan and colleagues. *Allele-Specific HLA Loss and Immune
Escape in Lung Cancer Evolution*. Cell (2017).
[DOI 10.1016/j.cell.2017.10.001](https://doi.org/10.1016/j.cell.2017.10.001).
[Primary public article](https://pmc.ncbi.nlm.nih.gov/articles/PMC5720478/),
[public full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5720478/fullTextXML).

Exact aggregate numerators and denominators come from Results, *Prevalence
and Timing of HLA Imbalance and Loss across NSCLC*, paragraph beginning
“In 36/90 (40%) of NSCLCs”, referring to Figures 2A and 2B. Combined NSCLC
and both component histologies are retained. They overlap and must not be
analysed as three independent observations. CSV fractions are ratios of the
published integers; extraction half-width zero means no rounding in those
counts, not no sampling uncertainty. These are early-stage, resected,
multi-region TRACERx tumours, not the ICI outcome cohorts.

**Reuse:** the article identifies the Francis Crick Institute copyright and
an explicit [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license.
Only published aggregate counts and their labelled arithmetic ratios are
included.

## Additional search and excluded material

* **Sewastianik et al. 2025**, *Allele-specific HLA LOH in solid tumors:
  distinct patterns by tumor type and potential prognostic relevance*,
  [DOI 10.1136/jitc-2025-012435](https://doi.org/10.1136/jitc-2025-012435).
  Its [public article XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12421593/fullTextXML)
  establishes CC BY-NC rights and a reanalysis of the Montesion cohort. Its
  main scientific distinction is between loss of any HLA allele and loss of
  a particular allele. The broad graphical panel was not available through
  the inspected public figure routes. Its spreadsheets were not downloaded:
  the inspected descriptions did not establish an aggregate-only schema.
  No rates from that unavailable panel were invented or substituted.
* **Lozac'hmeur Tempus** frequencies are HLA-A-only and therefore excluded
  from the any-HLA primary panel. Other comparison columns were not selected
  opportunistically by cancer type. Its opposite-direction outcome findings
  were already retained in PR63 and reinforce why no uniformly adverse
  HLA-response interpretation is assumed as established fact here.
* **Source-data workbooks and genomic records:** DASH source-data files
  include individual-level measurements; its sequencing data are controlled.
  These were not retrieved. The same restriction applies to Montesion and
  the other cohort records. Public discoverability alone was not treated as
  redistribution permission.
* **Access:** public Europe PMC article XML supplied readable articles when
  some publisher/PMC browser retrievals returned challenges or errors. The
  DASH figure came from its public publisher image endpoint. No access
  control or sign-in wall was circumvented.

## Mapping and interpretation safeguards

Source organ labels remain intact. Kidney, liver and uterine categories do
not automatically establish renal-cell carcinoma, hepatocellular carcinoma
or endometrial histology. Bladder is not automatically every urothelial
primary. Colorectal frequencies are not MSI-stratified; breast frequencies
are not receptor-stratified. The two DASH lung histologies cannot simply be
averaged into an NSCLC predictor without justified aggregate weights.
These mappings require explicit flags and exact-mapping sensitivity analyses.

All panels measure prevalence, not a transition hazard. Different assay
sensitivity, tumour purity, stage, region sampling and eligibility rules can
change the estimated prevalence. A source-checked numerical transcription,
its hash, agreement among agents, or passed software tests cannot establish
the truth of the underlying biological measurement or validate the theory.
