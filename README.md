# Open Cancer

**An open, provenance-first index of cancer research data and reproducible research contributions.**

Open Cancer aims to make public cancer research easier to **discover, understand, verify, and build on** — for people and AI agents. The long-term aspiration is to accelerate research, not to suggest that a catalog alone will cure cancer.

## Open Cancer Explorer

Search by accession, disease, or source. Open a dataset to inspect its unit of observation, cohort, access conditions, reuse terms, and retrieval link. The EGFR lung adenocarcinoma pilot connects its research question to the exact queries, code, and evidence tables behind the published snapshot.

[Open the Explorer (historical public demo)](https://open-cancer-explorer.lagging.chatgpt.site/) · [Explore the EGFR pilot](https://open-cancer-explorer.lagging.chatgpt.site/?view=pilot) · [Download the demo JSON index](https://open-cancer-explorer.lagging.chatgpt.site/api/v1/index.json) · [API and local build guide](docs/explorer-api.md)

**Publication freshness warning:** The public demo was confirmed on November 11, 2026 to be serving an **older Git revision** than `main`. It should not be mistaken for the latest validated research index. Check the [automated live-release status](https://github.com/dovvnloading/open-cancer/issues/29) and [deployment instructions](docs/explorer-releases.md); until a current production deployment passes its smoke test, build the latest Explorer from a pinned source checkout. A passing repository build is not a live deployment.

## What this repository does

1. **Discover:** a curated, machine-readable directory of upstream cancer resources, with access conditions and attribution guidance.
2. **Question:** small, reviewable research tasks that people and agents can work on.
3. **Contribute:** a consistent format for reproducible results with source links, query parameters, methods, checksums, and limitations.
4. **Validate:** automatic structural validation and human scientific review (not interchangeable).

**It is a federated index, not a bulk mirror.** Large biomedical data and protected records should remain with the authorized source. GitHub stores metadata and suitably licensed, non-sensitive derived artifacts when they meet the [data policy](DATA_POLICY.md).

## Explore the initial source catalog

| Source | Focus | Access |
| --- | --- | --- |
| [NCI GDC](catalog/sources/gdc.json) | Harmonized cancer genomics | Mixed open / controlled |
| [cBioPortal](catalog/sources/cbioportal.json) | Cancer studies and molecular data | Public portal; check study reuse |
| [CIViC](catalog/sources/civic.json) | Curated cancer variant evidence | Open; curated data CC0 |
| [ClinicalTrials.gov](catalog/sources/clinicaltrials.json) | Clinical study metadata | Public registry |
| [NCBI GEO](catalog/sources/geo.json) | Functional genomics studies | Public repository |
| [DepMap](catalog/sources/depmap.json) | Cancer cell-line dependencies | Public downloads; check release terms |

These are **source-directory entries**, not replicated datasets. Availability does **not** imply permission to redistribute every source file. Each record links to upstream documentation and distinguishes access from reuse. See [catalog instructions](catalog/README.md).

## Explore specific studies and releases

The [dataset catalog](catalog/datasets/README.md) now includes individually identified public releases and cohorts, beginning with [CIViC nightly clinical evidence](catalog/datasets/civic-nightly-clinical-evidence.json), [cBioPortal TCGA LUAD PanCancer Atlas](catalog/datasets/cbioportal-luad-tcga-pancan-atlas-2018.json), and [GDC TCGA-LUAD](catalog/datasets/gdc-tcga-luad-project.json). Each entry records its accession, retrieval URL, population definition, observation unit, reuse conditions, version mutability and verification date.

```bash
python scripts/search_catalog.py --disease DOID:3910
python scripts/search_catalog.py --source civic --format json
```

## Scientific validation challenge

The first EGFR/lung adenocarcinoma pilot now has an [independent reviewer kit](review_packets/egfr-luad-v1/README.md), a [domain-expert checklist](review_packets/egfr-luad-v1/SCIENTIFIC_CHECKLIST.md), and an [outreach shortlist](review_packets/egfr-luad-v1/REVIEWER_TARGETS.md). A checksum-verified review bundle can be generated with:

~~~bash
python scripts/build_review_packet.py --output /tmp/egfr-review --archive /tmp/egfr-review.zip
~~~

We welcome independent critical cancer-genomics reviews through [issue #3](https://github.com/dovvnloading/open-cancer/issues/3). The computational method has been reproduced by project contributors, but the scientific interpretation is **not independently validated** and the analysis is not suitable for clinical care.

## Start contributing

- **Add or improve a source:** submit `catalog/sources/<id>.json` using [schema](schemas/source.schema.json), citing official documentation.
- **Propose a research task:** add `tasks/<id>.json` using [schema](schemas/task.schema.json).
- **Run an agent or analysis:** pick an [open task](tasks/egfr-luad-evidence-coverage.json), follow [AGENTS.md](AGENTS.md), and submit a reproducible manifest and artifacts. The [result schema](schemas/result.schema.json) is machine-readable.
- **Review a contribution:** use the [scientific review protocol](REVIEW_PROTOCOL.md) and [version-bound review schema](schemas/review.schema.json); independently check methods, sources, and uncertainty rather than treating schema-valid output as medically correct.

Full process: [CONTRIBUTING.md](CONTRIBUTING.md). For featured research, read the [publication governance and required human approval runbook](docs/publication-governance.md). Human review requirements must be enforced by GitHub branch protection; a CODEOWNERS file alone is not an approval gate.

## Validate locally

Requires Python 3.11+:

```bash
python -m pip install -r requirements-dev.txt
python scripts/validate.py
python -m unittest discover -s tests -v
```

Validation runs on pull requests. It checks format, cross-references, local artifact existence, and SHA-256 hashes; it **does not** confirm that biomedical assertions are true.

## Roadmap

- **v0 (this foundation):** source catalog, structured research tasks, agent submission contract, reproducible validation.
- **v0.2:** curated *dataset-level* records (accessions, DOI, release, license, cohort metadata), not just portals.
- **v0.3:** source-specific connectors and reproducible ingestion for permitted, non-sensitive metadata.
- **v0.4:** terminology mapping (e.g., gene identifiers, disease ontologies), search API, provenance and quality dashboards.
- **v1:** community research queue, domain-expert review, benchmarked agents, and versioned derived datasets with auditable evidence.

Following [FAIR data principles](https://data.nih.gov/research/fair-data-principles) means data should be findable, accessible under documented conditions, interoperable, and reusable with provenance — not necessarily universally downloadable.

## Scope and safety

This project is for **research**, not medical advice or clinical decision-making. Never submit patient-identifying or potentially re-identifiable data, credentials, controlled-access data, or copyrighted material without authorization. Agent-generated claims must remain labeled **unreviewed** until independent scientific review. Code in this repository is MIT licensed; third-party datasets retain **their own** permissions and licenses. Read [DATA_POLICY.md](DATA_POLICY.md).
