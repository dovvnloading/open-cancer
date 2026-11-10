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
