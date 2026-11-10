# Source catalog

Each `catalog/sources/<id>.json` is a versioned, human-readable JSON record with a machine-readable schema. Records describe **portals/resources**, not individual harmonized datasets. The next iteration will add dataset-level accession, DOI, version, licensing, cohort, and schema mappings.

Key distinctions:

- `access`: whether public browsing/downloads are possible, or mixed with protected data.
- `reuse.status`: whether a permissive license has been **explicitly verified** for the referenced content, not merely inferred from openness.
- `provenance_url`: documentation used to substantiate the catalog record.
- `reviewed_on`: editorial date the reference was reviewed, **not** a live availability or validity guarantee.

Use `verify-upstream` unless the provider unambiguously documents a permissive license. For mixed-access resources, permissions vary by project, data type, and consent.

To add a source, copy an existing entry, change every field, confirm the official URLs, run the validation script, and submit a PR. Please do not submit copied sensitive data.
