# Explorer JSON API v1

The Explorer serves a read-only snapshot under `api/v1/` relative to the site's home page. Each endpoint accepts HTTP GET and returns UTF-8 JSON. Download the index for a consistent snapshot across collections. Filter or paginate locally.

| Endpoint | JSON value |
| --- | --- |
| `api/v1/index.json` | Object containing `version` (1), `repository`, `revision` (full Git commit SHA), and all five collections below |
| `api/v1/datasets.json` | Array of accession records with disease, modality, observation unit, cohort, access, reuse, release and verification date |
| `api/v1/sources.json` | Array of source-directory records with publisher and upstream documentation |
| `api/v1/questions.json` | Array of research tasks with question, input source IDs, status and requested deliverables |
| `api/v1/pilots.json` | Array of submitted result manifests with inputs, execution details, artifacts, findings and limitations |
| `api/v1/reviews.json` | Array of version-bound review records, including disputes and revision requests. An empty array means no review records are published. |

For example, from a locally served checkout, look up the exact disease identifier across the accession catalog.

```sh
curl -fsS http://localhost:8000/api/v1/datasets.json |
  jq '.[] | select(.disease.ontology_id == "DOID:3910") | {id, title, unit_of_observation, landing_url}'
```

## Provenance and joins

Dataset `source_id`, question `input_source_ids`, and pilot `inputs[].source_id` point to source `id`. Dataset, source, question and review objects preserve their [repository schema fields](../schemas) and add `record_url`, a GitHub link pinned to the indexed revision. The upstream retrieval and reuse links inside these records retain their original meaning. Entries are ordered by repository path.

A pilot is identified by the pair `(task_id, run_id)`. Its `task_id` points to question `id`. Each pilot adds `manifest_url`, the SHA-256 of the original manifest bytes in `manifest_sha256`, `code_url`, and a commit-pinned `url` per artifact. The EGFR/LUAD pilot embeds the contents of its declared, checksum-validated v1 `coverage.json` in `coverage`. Other runs without that supported artifact have `coverage: null`. Coverage includes the source-specific units, filters, availability and limitations alongside the counts. A missing count keeps its original null value.

Join each review to a pilot on `task_id`, `run_id`, and `reviewed_manifest_sha256 == manifest_sha256`. Keep each review's scope, decision, reproduction outcome, checks, follow-up and rationale together. Read the [review protocol](../REVIEW_PROTOCOL.md) for the distinction between methods reproduction, scientific interpretation assessment and clinical validity. The review collection covers records under `reviews/`. GitHub PR comments and fictional schema examples are outside this collection.

## Build and serve

From a checkout of the revision to publish, run the following commands. The command checks that its source files match the checked-out commit, so each provenance link opens the version shown in the Explorer.

```sh
python -m pip install -r requirements-dev.txt
python scripts/build_explorer.py --output _site
python -m http.server 8000 --bind 127.0.0.1 --directory _site
```

Open the server's root URL for the Explorer or append an endpoint path for JSON. Use a new or empty output directory for each build. The build runs repository validation, including artifact checksums and review-to-manifest bindings, before publishing any files. It reads checked-in metadata and public aggregate snapshots locally. Running it twice from the same clean checkout produces identical JSON. For a fork, set `--repository OWNER/NAME`. A source archive without Git metadata also needs `--revision` with its full source commit SHA.
