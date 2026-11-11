# Dataset-level records

This directory records **specific upstream accessions, studies, or named releases**, rather than portal names. Machine-readable records follow [dataset.schema.json](../../schemas/dataset.schema.json).

A dataset record distinguishes:

- `accession`: upstream identifier and namespace
- `unit_of_observation`: the unit being indexed, not assumed to be a patient
- `cohort_definition`: population/release/filters, including important caveats
- `access`: whether some or all underlying content is protected
- `reuse`: rights to redistribute, **not** synonymous with public web accessibility
- `release.mutable`: whether the upstream resource changes and needs dated provenance and content checksums
- `verified_on`: when the record was checked, not a guarantee it remains live

These are metadata records. They do **not** include raw patient files, hidden authentication tokens, or unlicensed datasets. An accession or portal may contain many datasets with different terms. In the future we'll extend each record with release-specific checksums, persistent DOI citations, cohort ontologies, and validated schema mappings.

## Query offline

```bash
python scripts/search_catalog.py --text lung
python scripts/search_catalog.py --disease DOID:3910 --format json
python scripts/search_catalog.py --source civic
```

To contribute, copy an existing record, verify upstream links and reuse terms, and run `python scripts/validate.py`.

## Refresh accession metadata

The connectors read the GDC TCGA project collection, public cBioPortal studies whose IDs end in `_tcga_pan_can_atlas_2018`, and three named CIViC nightly TSV exports. They use anonymous public endpoints. GDC projects and cBioPortal studies may cover overlapping TCGA participants, and their cohort boundaries can differ. A CIViC record describes an export with a different observation unit.

Run a single lookup from the repository root. The result is a JSON array.

```bash
python scripts/refresh_catalog.py --source gdc --accession TCGA-BRCA
python scripts/refresh_catalog.py --source cbioportal --accession brca_tcga_pan_can_atlas_2018
python scripts/refresh_catalog.py --source civic --accession nightly-AssertionSummaries.tsv
```

To discover everything within this scope, choose a new directory. The command writes one candidate record per accession after all three sources respond successfully.

```bash
python scripts/refresh_catalog.py --output catalog-candidates
```

Review the diff before copying changes into `catalog/datasets/`. Keep existing editorial annotations such as the LUAD disease ontology mapping and the separately indexed GDC EGFR query. The connector output covers the upstream accessions. The filtered research query has its own retrieval script under `results/`. A refresh can reveal removals or changes upstream. Check those individually before retiring any catalog entry.

| Source | Upstream checks and mappings |
| --- | --- |
| [GDC](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/) | Filters on `program.name=TCGA`, requests project name and data categories, and checks the returned total. Disease uses the project name. Modality tags come from the listed categories. The [status service](https://api.gdc.cancer.gov/status) supplies the release label at retrieval. |
| [cBioPortal](https://docs.cbioportal.org/web-api-and-clients/) | Checks the study list against the API's `META` count, then selects public PanCancer Atlas study IDs. The disease label comes from the portal's cancer type and can be broader than the study title. Modality tags reflect positive aggregate sample counts for the mapped assays. Records retain the portal import date and link to the PubMed IDs supplied in the study metadata when available. |
| [CIViC](https://docs.civicdb.org/en/latest/using/data_releases.html) | Verifies the identifier and review columns in the header of the nightly clinical evidence, molecular profile, and assertion exports. No export rows are collected. The project uses [CIViC's published data terms](https://docs.civicdb.org/en/latest/about.html). |

`verified_on` records the UTC day of the lookup. A source's reported release or import date describes the version observed. The live retrieval links can change. Ontology IDs are included only on curated mappings. The connectors do not infer relationships between cohorts or between similar disease labels.
