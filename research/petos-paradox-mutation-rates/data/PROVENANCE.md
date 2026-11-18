# Aggregate data provenance and observation boundaries

**Nonclinical, unreviewed computational research. No medical advice.**

These small tables contain published species aggregates and explicit source
metadata. They contain no individual human or animal records. Access-date
labels use **2026-11-18, America/New_York**. The protocol was published before
fitting; reading source definitions and transcribing literal published values
preceded analysis. The sources are not an independent validation cohort.

## Vincze cancer mortality table

Source: Vincze et al., *Cancer risk across mammals*, Nature (2022),
[doi:10.1038/s41586-021-04224-5](https://doi.org/10.1038/s41586-021-04224-5).
The article's **Data availability** and **Code availability** sections designate
the authors' public reproducibility repository. Its README explicitly identifies
`data.csv` as containing species-specific cancer mortality, sample sizes, body
mass and remaining adult life expectancy.

The table used here is the openly readable
[publication-linked aggregate CSV](https://github.com/OrsolyaVincze/VinczeEtal2021Nature/blob/d2f8e2dff8e64393ce3512f8b8490507fa9ae803/data.csv).
It was read as text through the GitHub connector. The table is pinned at commit
`d2f8e2dff8e64393ce3512f8b8490507fa9ae803`; its source Git blob is
`d893878432e530a4de3751f95986fe35aca3b128`. `source_row` is the original
numeric identifier in that table, not the line number after extraction.
Each retained row records its binomial and original column locators.

`vincze_mortality.csv` preserves the source decimal strings and integer
counts. Source `NA` is represented by an empty field; it is never a zero.
The following field map is part of the extraction method.

| Local field | Source column | Meaning |
|---|---|---|
| `female_mass_kg`, `male_mass_kg` | `FemaleMeanMass`, `MaleMeanMass` | Mean adult mass by sex, in kilograms. The analysis follows the article's equal-sex averaging convention for the species proxy. |
| `cmr` | `CMR` | Fraction of known-cause deaths for which neoplasia significantly contributed to death. |
| `icm` | `ICM` | Censoring- and delayed-entry-aware cumulative incidence of cancer mortality estimated by the authors. This is not interchangeable with `cmr`. |
| `adult_life_expectancy_days` | `lifeexp` | Expected remaining days after sexual maturity, not total age from birth, maximum longevity, or a mortality-quantile lifespan. |
| `total_animals` | `Sample` | Source animals including those still alive. |
| `total_deaths` | `Deaths` | All recorded deaths, including those without known pathological cause. |
| `cancer_deaths` | `Neoplasia` | Recorded neoplasia cases judged to contribute significantly to death. |
| `known_cause_deaths` | `knownDeaths` | Deaths with recorded pathological cause; the denominator of CMR. |

Per-species CMR confidence intervals are **absent** from this source table.
Both `published_cmr_ci_lower` and `published_cmr_ci_upper` therefore remain
empty, and `ci_status` records their absence. The analysis computes Wilson
intervals from the published counts only after the protocol freeze. These are
derived binomial intervals under an explicit sampling assumption, not intervals
printed by Vincze and not uncertainty for the ICM estimator. They cannot
quantify selection into necropsy, missing causes, zoo conditions, phylogenetic
dependence or mismatch with incidence at a fixed age.

The article's **Rights and permissions** section specifies
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The reuse basis is
the authors' published aggregate table designated by that open-access article,
with attribution and extraction changes retained here. No separate repository
LICENSE was observed; this note does not invent one or extend permission to
independent third-party records. The article explicitly distinguishes its
public aggregate outputs from restricted underlying Species360 data. No
restricted record was requested, accessed or redistributed. The project's
MIT license does not relicense the source aggregates.

## Mutation rates and life-history traits

The mutation and life-history CSV rows retain their own article citation,
versioned table locator and reuse terms. Their route is Cagan's published
somatic-rate estimates as reproduced in the article-linked, openly readable
`Mutation_Data.csv` associated with Compton's separately licensed preprint.
That source is Compton et al., *Cancer Prevalence Across Vertebrates*,
[preprint version 1](https://doi.org/10.21203/rs.3.rs-3117313/v1), published
6 July 2023. Its [public full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10350200/)
explicitly states CC BY 4.0 and links the authors' repository in **Data and
materials availability**. The
[historical aggregate table](https://github.com/zacharycompton/cancerAcrossVertebrates/blob/3b4db10c7c38c63aeebfc9a48ce4bfa376dfd2bc/Mutation_Data.csv)
is pinned at its creation commit
`3b4db10c7c38c63aeebfc9a48ce4bfa376dfd2bc` of 15 December 2022, preceding
that preprint. Its text was checked against the later repository version
`4819e0933358c5fc7eab483a0c05161ea8c074b7` and is unchanged. The data
dictionary confirms grams for adult weight and months for maximum longevity.
The Cagan original estimates retain their original CC BY attribution; an
[institutional deposit](https://doi.org/10.17863/CAM.83740) records that license.
The final Compton Cancer Discovery article has a more restrictive
CC BY-NC-ND license; that final-article license is not used as blanket permission
to republish its cancer-outcome dataset. Only the explicitly documented
mutation/life-trait source route is used, and Cagan retains attribution.

Mean annual SBS and the slope from a zero-intercept regression are different
estimators. Printed confidence limits attached to the latter are not reassigned
to the former. Missing species-level mutation sample sizes or uncertainty stay
missing. A negative regression confidence limit is retained as a source fact;
it is not silently replaced by a positive value for logarithmic calculations.
Genome-wide SBS per year is not an observed per-driver transition intensity.

The reproduced life traits concern maximum longevity and adult weight. Their
units and original values are preserved in their rows. Neither maximum
longevity nor Vincze's remaining adult expectancy is Cagan's mortality-quantile
lifespan. Cross-source comparisons use explicitly labelled sensitivity clocks.
They do not assert that the cohorts, tissues, sex distributions, environments
or age windows match.

## Access limits and excluded substitutes

No browser file download was used and no browser download prompt was accepted.
No binary supplementary spreadsheet or archive was extracted. Openly readable
article text, supplementary text where accessible and publication-linked text
tables supplied the permitted records. Failure to retrieve a supplement was
handled as missingness, not as an invitation to infer plotted coordinates or
obtain individual records elsewhere.

Cagan's original full species-rate and life-history supplementary spreadsheets
were not recovered through a permitted readable-table route. The available
reproduced rows are a selected subset, not the complete original mutation study.
Original mortality-quantile lifespans and some species-level mutation sample
sizes therefore remain unavailable. No claim to replicate Cagan's complete
lifespan analysis follows from the substitute maximum-longevity regression.

Zheng's JITC article links Vincze and Cagan and makes a quantitative-resolution
claim. Its main methods were readable, but the supplemental PDF could not be
read through the available web reader. The analysis consequently audits the
stated mechanistic interpretation; it does not claim exact reproduction of
Zheng's undisclosed-to-this-analysis matched panel or transformations.

Species matching uses the explicit binomials in the published tables. Source
common names and subspecies labels are retained where supplied, with any
species-level collapse documented. Domestic breeds, similarly named relatives,
and unspecified species are not inserted to enlarge the overlap. No row is
removed because its observed cancer outcome is inconvenient.

## What the data do not measure

The mortality panel is an aggregate from managed populations with diagnosed
causes of death. Its estimated mean adult lifetime is not the deterministic
exposure horizon of an initially clean lineage. The mutation measurements
concern normal intestinal crypts, whereas the cancer endpoint pools anatomical
sites. Adult body mass is a proxy for effective susceptible targets. None of
these tables supplies the absolute driver conversion, independent target
count, clonal progression or detection law, competing-death schedule, or a
biological count of rate-limiting steps. Exact computations and content hashes
establish reproducibility of the specified analysis, not those missing bridges.
