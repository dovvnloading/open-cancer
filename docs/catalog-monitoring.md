# Public catalog freshness and drift monitoring

The [catalog monitor](../.github/workflows/catalog-monitor.yml) runs every Monday at 09:23 UTC and supports manual dispatch. It checks the *metadata* connectors for GDC TCGA projects, cBioPortal public PanCancer Atlas studies, and CIViC nightly export headers. It makes bounded public requests, retries transient HTTP/transport failures no more than three times with finite backoff, and never fetches patient-level rows.

The scan compares **accession presence** against the committed catalog and connector-observed metadata (study names, accessions, disease labels, modalities, releases, citations and retrieval links) against the **previous upstream observation** retained in a GitHub Actions run artifact. Editorial descriptions, ontology mappings, cohort explanations and permission statements are not upstream-native metadata, so differences in those fields are not misrepresented as new upstream drift. It detects additions, removals, unexpected source-count shrinkage, duplicate identifiers, empty upstream responses, malformed responses, and source outages. The custom EGFR-filtered GDC occurrence record is not a GDC project and is excluded from project counts.

For each of the three sources, the monitor also probes the catalogued reuse-terms URL with a bounded HEAD request (GET fallback only on HTTP 405). Failures are reported. When the upstream site supplies ETag or Last-Modified headers, changes between runs trigger a source-rights recheck. This is not a content-level license audit: unchanged headers do not establish unchanged licensing, and the discovery connector necessarily inherits existing editorial reuse fields.

The first successful observation establishes a connector-specific comparison baseline without silently changing the curated records. A source outage retains the last successful observation and last_healthy date for comparison once connectivity resumes, as long as a prior report artifact is available. A push changing the monitor, source connector or source definitions on main also triggers one scan; it does not rely only on the weekly schedule.

Every scheduled or triggered run publishes an **Actions run summary** and a 90-day `catalog-health-report` artifact containing the public-source observations, last_checked / last_healthy times, changed fields, source failures and terms-link metadata. The next run retrieves the latest retained artifact as its comparison baseline. An absent/expired previous artifact is explicitly labeled `previous_monitor_artifact_used=false`: historical drift is **not assessed** for the gap. Actual drift/outages cause a **failed Actions run** with the source-specific evidence. No status-only GitHub Issue is created or reopened. A sustained, reproducible problem can be filed separately using the engineering Issue template.

**Nothing is automatically published or overwritten.** Maintainers must compare accession changes and source rights against authoritative documentation, then create a normal reviewed PR before updating catalog records. A healthy status means these specific bounded checks found no drift; it does not assert scientific correctness, ongoing redistribution permission, clinical validity, or endpoint uptime at every moment.

Offline tests: \`python -m unittest discover -s tests -p test_catalog_monitor.py -v\`.

Online manual run: \`python scripts/monitor_catalog.py --output catalog-health-report.json\`.


The connector inherits source access/rights descriptions from editorial templates. Therefore monitoring **cannot** infer that source reuse terms have legally changed from their regenerated metadata; it only highlights URL availability and HTTP ETag/Last-Modified changes when present. A human must review official source terms, consent conditions, and source-version changes independently. An upstream field change is not a request to overwrite editorial text automatically.


## Status storage and artifact retention

Actions artifacts expire after 90 days (or sooner if manually deleted), so
monitoring is **not an immutable archive**. A scheduled scan uses the most recent
retained report from `catalog-monitor.yml`, including on a previous failed run.
The monitor keeps last-good public metadata observations across outages and never
treats missing historical context as proof of unchanged sources. A maintainer
requiring permanent evidence should publish an independently reviewed dated,
rights-safe report or tagged release, not a status-only Issue.

The historical closed **Catalog metadata health ledger (#27)** is retained as a
record of prior automation but is no longer read, written, created or reopened by
the monitor. New reports and errors appear in the
[Actions workflow history](https://github.com/dovvnloading/open-cancer/actions/workflows/catalog-monitor.yml).
The old `scripts/catalog_ledger.cjs` codec and tests remain for inspecting
legacy data, not for driving new issue updates.
