# Aggregate incidence data provenance

**Nonclinical computational research. No medical advice.** These are published,
population-level summary rates, with no individual-level records. The analysis
is unreviewed and does not establish a biological driver count.

## Source and permitted reuse

The sole incidence source is the National Cancer Institute (NCI),
[SEER Cancer Statistics Review, 1975–2017](https://seer.cancer.gov/archive/csr/1975_2017/),
edited by Howlader N et al., Bethesda, Maryland. This archived release is based
on the November 2019 SEER data submission and was released April 15, 2020.
It is deliberately a fixed historical release, not a claim about current
incidence. All source pages and the reuse statement were accessed
**2026-11-17** (America/New_York).

The archive landing page's **Suggested Citation** section expressly places the
report in the public domain and permits copying or reproduction, requesting
source citation. This exact-source permission is the `reuse_terms_url` on every
CSV row. The [NCI reuse policy](https://www.cancer.gov/policies/copyright-reuse)
also permits reuse of NCI text unless a different restriction is stated.
No third-party graphics, logos, or identifiable records are redistributed.
These source rates are not relicensed by the repository's MIT license.

## Selection and transcription

`seer_age_specific_2013_2017.csv` contains **44 rows**: four cancer sites and
eleven consecutive five-year age bands, ages 30–34 through 80–84. Every value
comes from the **Age-Specific SEER Incidence Rates, 2013–2017** subsection,
for **SEER 21 areas**, **all races**, **invasive** disease. Both sexes are used
except for breast, where the female column is used. These are incidence rates
per 100,000 population per year, not mortality rates or the age-adjusted
summary rows that appear above them in the original tables.

| CSV site | Source table | Selected column | HTML table | One-page source PDF |
|---|---|---|---|---|
| `pancreas` | 22.7 | All Races / Both Sexes | [Table 22.7](https://seer.cancer.gov/archive/csr/1975_2017/browse_csr.php?pageSEL=sect_22_table.07&sectionSEL=22) | [PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_22_table.07.pdf) |
| `colon_and_rectum` | 6.10 | All Races / Both Sexes | [Table 6.10](https://seer.cancer.gov/archive/csr/1975_2017/browse_csr.php?pageSEL=sect_06_table.10&sectionSEL=6) | [PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_06_table.10.pdf) |
| `lung_and_bronchus` | 15.10 | All Races / Both Sexes | [Table 15.10](https://seer.cancer.gov/archive/csr/1975_2017/browse_csr.php?pageSEL=sect_15_table.10&sectionSEL=15) | [PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_15_table.10.pdf) |
| `female_breast` | 4.11 | All Races / Females | [Table 4.11](https://seer.cancer.gov/archive/csr/1975_2017/browse_csr.php?pageSEL=sect_04_table.11&sectionSEL=4) | [PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_04_table.11.pdf) |

The CSV preserves the printed decimal strings. The only transformations are
selecting these columns and age bands, giving site names machine-readable
labels, and converting age labels into half-open interval endpoints: `30-34`
means `[30,35)`. No smoothing, imputation, rescaling, or adjustment is applied
to the stored rate column. The open-ended 85+ band is excluded because it has
no finite upper endpoint. Younger bands are excluded consistently across sites;
their occasional suppression is not interpreted as zero incidence.

All 44 selected values were checked against the HTML table and separately
against the corresponding PDF text extraction, including column headings and
the division between incidence and mortality. These are two representations of
one NCI source, not independent epidemiological datasets. The linked official
CSV downloads are recorded for traceability, but this retrieval interface
could not read their response content type; they were not the transcription
source. Source PDFs and HTML were not copied into the repository.

### Compact transcription audit

Each entry below has units of cases per 100,000 population per year and is
printed to one decimal place, matching the selected source column.

| Age | Pancreas | Colon and rectum | Lung and bronchus | Female breast |
|---|---:|---:|---:|---:|
| 30–34 | 0.8 | 6.1 | 1.2 | 29.4 |
| 35–39 | 1.4 | 10.8 | 2.6 | 64.4 |
| 40–44 | 2.8 | 19.4 | 6.3 | 130.7 |
| 45–49 | 5.9 | 33.3 | 16.5 | 197.1 |
| 50–54 | 11.8 | 61.8 | 42.6 | 236.4 |
| 55–59 | 20.6 | 67.8 | 84.6 | 269.0 |
| 60–64 | 32.6 | 88.2 | 134.6 | 342.3 |
| 65–69 | 48.3 | 118.1 | 215.5 | 426.5 |
| 70–74 | 64.3 | 146.2 | 318.8 | 470.2 |
| 75–79 | 82.4 | 188.3 | 390.5 | 458.9 |
| 80–84 | 95.6 | 234.6 | 406.5 | 413.5 |

## Fields and statistical limits

- `age_start` is inclusive and `age_end_exclusive` is exclusive; all endpoints
  are in years. `age_group` preserves the source label.
- `rate_per_100000` is the published age-specific rate. `reported_decimal_places`
  is `1`; it records display precision, not standard error or accuracy.
- `case_count` and `person_years` are deliberately **blank**. These tables do
  not provide the selected age-by-site cell counts and person-years. Neither
  can be recovered uniquely from a rounded rate. Do not treat a blank as zero
  or replace it with an invented effective sample size.
- `source_name`, `source_table`, `source_section`, `source_column`, the three
  source URLs, `access_date`, and `reuse_terms_url` identify each individual row.
- The table footnote suppresses statistics based on fewer than sixteen cases.
  None of the selected cells is suppressed; this does not reveal the actual
  counts.

These are pooled **period** rates, not a followed birth cohort. Age bands can
differ in birth cohort, exposures, detection, screening, competing mortality,
and tumor subtype composition. A constant-rate homogeneous biological model
does not automatically apply to this population summary. The source sites
also need not match sequencing cohorts or histological subtypes.

The reproducible analysis therefore supports shape diagnostics, explicit
assumption-dependent compatibility statements, and observation-model
counterexamples. It does **not** support a sampling confidence interval,
likelihood-ratio significance test, or cancer-specific statistical minimax
number using these rates alone. Any Poisson experiment used elsewhere in this
directory has specified **synthetic** exposure and is labeled accordingly.
Allowing a half-unit in the last reported decimal for a rounding sensitivity
check, if used, is an explicitly stated rounding convention; it is not a
confidence interval.
