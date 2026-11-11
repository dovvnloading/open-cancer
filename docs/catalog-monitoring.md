# Public catalog freshness and drift monitoring

The [catalog monitor](../.github/workflows/catalog-monitor.yml) runs every Monday at 09:23 UTC and supports manual dispatch. It checks the *metadata* connectors for GDC TCGA projects, cBioPortal public PanCancer Atlas studies, and CIViC nightly export headers. It makes bounded public requests, retries transient HTTP/transport failures no more than three times with finite backoff, and never fetches patient-level rows.

The scan compares **accession presence** against the committed catalog and connector-observed metadata (study names, accessions, disease labels, modalities, releases, citations and retrieval links) against the **previous successful upstream observation** stored in the ledger. Editorial descriptions, ontology mappings, cohort explanations and permission statements are not upstream-native metadata, so differences in those fields are not misrepresented as new upstream drift. It detects additions, removals, unexpected source-count shrinkage, duplicate identifiers, empty upstream responses, malformed responses, and source outages. The custom EGFR-filtered GDC occurrence record is not a GDC project and is excluded from project counts.

For each of the three sources, the monitor also probes the catalogued reuse-terms URL with a bounded HEAD request (GET fallback only on HTTP 405). Failures are reported. When the upstream site supplies ETag or Last-Modified headers, changes between runs trigger a source-rights recheck. This is not a content-level license audit: unchanged headers do not establish unchanged licensing, and the discovery connector necessarily inherits existing editorial reuse fields.

The first successful observation establishes a connector-specific comparison baseline without silently changing the curated records. A source outage retains the last successful observation and last_healthy date for comparison once connectivity resumes. A push changing the monitor, source connector or source definitions on main also triggers one scan; it does not rely only on the weekly schedule.

Every scheduled or triggered run saves a short-lived workflow report artifact and updates one durable GitHub issue entitled **Catalog metadata health ledger (automated)**. That issue records last_checked and last_healthy independently for each source, source failures, source count comparisons, changed fields and link validity. The issue is reopened when review is needed, and new or changed alerts are appended in comments to preserve an audit trail beyond the workflow artifact retention period. The workflow deliberately fails on actionable drift or outages.

**Nothing is automatically published or overwritten.** Maintainers must compare accession changes and source rights against authoritative documentation, then create a normal reviewed PR before updating catalog records. A healthy status means these specific bounded checks found no drift; it does not assert scientific correctness, ongoing redistribution permission, clinical validity, or endpoint uptime at every moment.

Offline tests: \`python -m unittest discover -s tests -p test_catalog_monitor.py -v\`.

Online manual run: \`python scripts/monitor_catalog.py --output catalog-health-report.json\`.


The connector inherits source access/rights descriptions from editorial templates. Therefore monitoring **cannot** infer that source reuse terms have legally changed from their regenerated metadata; it only highlights URL availability and HTTP ETag/Last-Modified changes when present. A human must review official source terms, consent conditions, and source-version changes independently. An upstream field change is not a request to overwrite editorial text automatically.
