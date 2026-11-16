# EGFR/LUAD compute experiments

The pilot runs the existing nine EGFR/LUAD jobs against a pinned public CIViC record selection. For each of three alteration labels (L858R, exon 19 deletion and exon 20 insertion), the investigator considers three questions about terminology, missing denominators and potentially contradictory evidence. A separate critic receives the same records and the investigator's response. All outputs stay local, unreviewed, unpublished and nonclinical. [Workload](../workloads/egfr-luad-compute.json) · [Data policy](../DATA_POLICY.md)

## Launch a real pilot

Run from the repository root on a POSIX system with Python 3.11+ and operator-supplied OpenAI-compatible chat completion servers. Replace the example model IDs and URLs with ones your servers provide. If either server requires authentication, inject its `OPEN_CANCER_INVESTIGATOR_KEY` or `OPEN_CANCER_CRITIC_KEY` into the process environment through your secret manager.

```sh
python -m pip install -r requirements-dev.txt
export OPEN_CANCER_INVESTIGATOR_URL='http://127.0.0.1:8000/v1/chat/completions'
export OPEN_CANCER_INVESTIGATOR_MODEL='your-investigator-model-id'
export OPEN_CANCER_CRITIC_URL='http://127.0.0.1:8001/v1/chat/completions'
export OPEN_CANCER_CRITIC_MODEL='your-critic-model-id'
export OPEN_CANCER_LLM_TIMEOUT=120
export OPEN_CANCER_MAX_TOKENS=2400

# Pin the permitted CIViC projection once. This step makes no inference requests.
python scripts/research_pilot.py prepare --mode real --output .compute/source/evidence.json

# Read the pinned records, then run the nine jobs serially.
python scripts/research_pilot.py run --mode real --snapshot .compute/source/evidence.json \
  --output .compute/real-a
cat .compute/real-a/report.md
```

Use endpoints authorized for these public structured inputs. The pilot sends the job and its records to the selected investigator and critic. The critic also receives the candidate.

Both roles may share a server via `OPEN_CANCER_LLM_URL` and a default model via `OPEN_CANCER_LLM_MODEL`. Role-specific settings take precedence. A shared `OPEN_CANCER_LLM_KEY` is eligible only at the exact configured shared URL. A role with a different URL needs its own key if the endpoint requires one. Use HTTPS for remote inference. Keys belong only in environment variables, outside endpoint URLs, commands, checkpoint files and Git.

For token cost estimates, supply your actual USD prices per million tokens before starting the run. The variables are `OPEN_CANCER_INVESTIGATOR_INPUT_USD_PER_MILLION`, `OPEN_CANCER_INVESTIGATOR_OUTPUT_USD_PER_MILLION`, `OPEN_CANCER_CRITIC_INPUT_USD_PER_MILLION` and `OPEN_CANCER_CRITIC_OUTPUT_USD_PER_MILLION`. The estimate uses reported input/output tokens. Unavailable prices or token counts leave cost unknown. Server fees, local hardware time and energy are outside this estimate.

On the first pass each job has up to two calls, one per role, with 18 calls maximum. A failed investigation prevents its critique. Recover a failed or interrupted run by repeating the `run` command with the same output directory and `--resume`. Resume checks the pinned inputs, configuration and checkpoint integrity. Each failed role can be retried once across the run's lifetime, for an absolute maximum of 36 attempts. There are no automatic retries within a pass.

For an independent rerun, keep the same settings and use the original evidence snapshot in a fresh directory. The comparison reads both runs without calling either endpoint.

```sh
python scripts/research_pilot.py run --mode real --snapshot .compute/real-a/evidence.json \
  --output .compute/real-b
python scripts/research_pilot.py compare .compute/real-a .compute/real-b
```

Omitting `--snapshot` on a fresh real run retrieves CIViC again. Source data and provider deployments may change between runs. Check its SHA-256 and the projected snapshot hash before comparing. A model alias may point to new weights even when its configured name is unchanged. Pin a revision at your server where possible.

## Evidence included in the experiment

The pilot integrates Agent 1’s bounded, anonymous CIViC GraphQL connector from commit `f6d826a`. With `prepare --mode real`, the connector performs one bounded retrieval shared by all nine jobs, using its page, rate and retry limits. The nine-job benchmark applies a further limit of eight evidence items per alteration. [Selectors, pagination and provenance](#graphql-scope-rights-and-provenance) · [GraphQL connector](../scripts/civic_evidence.py) · [Pilot projection](../scripts/pilot_evidence.py)

An existing **native connector bundle** from `compute_fanout.py retrieve` can be passed directly to the benchmark. Supply an operator-trusted file and preserve its reported hash when transferring it. The pilot validates the catalog permission, bundle hash, record hashes, source query, scope and pagination before projection. A supplied bundle must match the requested real or fixture mode. Replaying it makes no source requests. An explicitly unavailable bundle stops this run before inference.

```sh
python scripts/compute_fanout.py retrieve --source civic --output .compute/civic-live.json
python scripts/research_pilot.py run --mode real --snapshot .compute/civic-live.json \
  --output .compute/real-from-connector
```

| Pilot input | Selection and provenance | Interpretation |
| --- | --- | --- |
| CIViC GraphQL | Accepted, unflagged, exact alteration and LUAD label with DOID 3910. The connector excludes deprecated profiles and diseases. It retains categorical evidence annotations, EID, source URL, record fingerprints, page variables and response hashes. | `complete` applies to the bounded server query and local eligibility rules. `truncated` means further pages exist. Each profile’s count then means matches observed so far, including when that count is zero. |
| CIViC nightly TSV, alternate input | On a fresh real retrieval returning `network-error` or `http-error`, the pilot can read the catalog-allowlisted CC0-1.0 nightly export once. It checks source and dataset permissions, pins the raw byte hash and retains provenance for the unsuccessful GraphQL attempt. It requires an unflagged, accepted item with the exact LUAD label, DOID 3910 and exact alteration. | The export lacks the GraphQL deprecated-profile and deprecated-disease flags. Results from the two selectors can differ. A complete or truncated GraphQL retrieval keeps its own observations, even if it yields zero matches. |
| GDC | No matching source connector is invoked. | Somatic mutation occurrences, cases and EGFR-assayed denominators remain unassessed. |
| cBioPortal | No matching source connector is invoked. | Mutation calls, profiled sample denominators and overlaps with GDC remain unassessed. |

For each alteration the selector cycles through groups defined by evidence direction and significance, taking the lowest numeric EID first within each group. Each job receives up to eight observed items and records both its observed match count and included count. Selection excludes compound profiles and unlisted aliases. The model receives a structured projection of a curated annotation. Repeated items can come from the same publication. Distinct directions or significance labels can reflect differences in study context, therapy or source definition. Primary papers and evidence statements were outside the projection. The historical November 10 pilot and its source counts belong to their original snapshot. [CIViC dataset terms](../catalog/datasets/civic-nightly-clinical-evidence.json) · [Historical pilot](../results/egfr-luad-evidence-coverage/pilot-20261110/README.md)

## Reading the outputs

| File under the run directory | What to inspect |
| --- | --- |
| `plan.json` | The nine content-addressed jobs and their question dimensions |
| `evidence.json` | Mode, source scope, connector bundle or TSV byte hash and projected snapshot SHA-256, retained record fields, matching counts and omissions |
| `run-config.json` | Model and prompt configuration, execution limits, fingerprints and operator price settings |
| `records/<job-id>.json` | Checkpoints with investigator observations, assumptions, limitations, counterevidence, critic issues, deterministic evaluations and failures |
| `benchmark.json` | Machine-readable counts, coverage, output fingerprints, timing, reported tokens and cost estimates |
| `report.md` | Technical report with execution status, source/model provenance, measured metrics and limitations |

Citation checks resolve the exact source, record ID and URL in the job's own packet. The deterministic checks compare a structured field assertion with its cited record or count the supplied rows in evidence-item units. Interpretations in free text still require assessment. The critic's unsupported-claim and contradiction counts measure flags raised by the model. The evaluator records invalid references and failed or unassessed assertions even when the critic accepts the candidate.

Coverage reports complete retrievals with no matches, truncated retrievals with no observed matches, and sources lacking retrieval separately. Truncated records can document a counterexample. A critic verdict of `no-internal-issues-found` requires a complete retrieval. An item used by multiple questions counts once in the unique-record coverage total. Compare independent runs for output consistency while accounting for failed jobs and differences in inputs or configuration. For each completed pair, the comparison reports whether the two successful requests returned different model identifiers or system fingerprints for either role. A missing provider value leaves that comparison unknown. Matching outputs measure repeatability. Scientific interpretation, independent validation and clinical assessment remain separate work.

## Reproduce the offline benchmark and validation suites

Install the Python validation dependency once (a local environment or preinstalled copy also works). The following tests and fixture runs use no live source or inference endpoint. Use a POSIX system with Python 3.11+ and Node.

```sh
python -m pip install -r requirements-dev.txt
python scripts/validate.py
python -m unittest discover -s tests -v
node --test tests/*.test.*
python scripts/research_pilot.py run --mode fixture --output .compute/fixture-a
python scripts/research_pilot.py run --mode fixture --output .compute/fixture-b
python scripts/research_pilot.py compare .compute/fixture-a .compute/fixture-b
cat .compute/fixture-a/report.md
```

Use fresh directories for independent fixture runs. The benchmark consumes the connector’s synthetic GraphQL fixture, which has two L858R items, one exon 19 deletion item and no exon 20 insertion item. The native fixture uses invented EIDs. The projected benchmark packet maps these to reserved `fixture.invalid` links. Fixture responses are deterministic, carry explicit synthetic labels and leave token usage and inference cost unmeasured.

## Three-job CIViC annotation workflow

This independent experiment uses one job per alteration to check annotation predicates. The fan-out coordinator enforces its v2 worker protocol for third-party and bundled adapters. Use it to inspect witness citations, counterexamples and abstentions on the original connector bundle. [CIViC-only workload](../workloads/civic-egfr-luad-annotations.json)

From the repository root with Python 3.11+, retrieve once and plan the three supported alteration jobs. The coordinator and adapter use the Python standard library.

```sh
python scripts/compute_fanout.py retrieve --source civic --output .compute/civic-live.json
python scripts/compute_fanout.py plan --spec workloads/civic-egfr-luad-annotations.json \
  --output .compute/civic-jobs.jsonl
```

Inspect the retrieval before using inference. Each row below comes from the retained snapshot and links to its original evidence item. The record's zero-based `page` points into the `pages` array.

```sh
python - <<'PY'
import json
from pathlib import Path
s = json.loads(Path(".compute/civic-live.json").read_text())
print(s["status"], s["failure"], s["retrieved_at"], s["sha256"])
for r in s["records"]:
    print(r["record_id"], r["molecular_profile"], r["evidence_type"], r["evidence_level"],
          r["url"], r["sha256"], s["pages"][r["page"]]["variables"])
PY
```

Set your own inference endpoint and model names. Separate subprocesses handle the two stages, and each can use a different model. `OPEN_CANCER_LLM_MODEL` supplies a fallback for either stage. An optional bearer credential belongs in `OPEN_CANCER_LLM_KEY`. The adapter sends only the selected job's structured records, question, dimensions, retrieval status and time for investigation. Critique adds the validated candidate. It allows HTTP on loopback, requires HTTPS remotely and refuses redirects.

```sh
export OPEN_CANCER_LLM_URL='http://127.0.0.1:8000/v1/chat/completions'
export OPEN_CANCER_INVESTIGATOR_MODEL='your-investigator-model-id'
export OPEN_CANCER_CRITIC_MODEL='your-critic-model-id'
python scripts/compute_fanout.py run --jobs .compute/civic-jobs.jsonl \
  --evidence .compute/civic-live.json --output .compute/civic-run \
  --investigator 'python scripts/openai_compatible_worker.py' \
  --critic 'python scripts/openai_compatible_worker.py' --workers 2 --max-jobs 3
```

Set `OPEN_CANCER_LLM_TIMEOUT` and `OPEN_CANCER_MAX_TOKENS` to bound individual inference requests. Set the coordinator's `--timeout` high enough for each worker process. User-controlled adapters can use any inference backend.

### Reproduce the three-job example offline

These commands use synthetic GraphQL data and deterministic subprocesses. They work without a CIViC connection, model server or credentials. Use fresh output paths the first time.

```sh
python scripts/compute_fanout.py retrieve --source civic --fixture \
  --output .compute/civic-fixture.json
python scripts/compute_fanout.py plan --spec workloads/civic-egfr-luad-annotations.json \
  --output .compute/civic-fixture-jobs.jsonl
python scripts/compute_fanout.py run --jobs .compute/civic-fixture-jobs.jsonl \
  --evidence .compute/civic-fixture.json --output .compute/civic-fixture-run \
  --investigator 'python scripts/demo_compute_worker.py' \
  --critic 'python scripts/demo_compute_worker.py'
python - <<'PY'
import json
from pathlib import Path
for path in sorted(Path(".compute/civic-fixture-run/records").glob("*.json")):
    body = json.loads(path.read_text())["body"]
    print(body["job"]["dimensions"]["alteration"], body["critique"]["verdict"],
          [(c["record_id"], c["relation"]) for c in body["critique"]["evidence"]])
PY
```

| Fixture job | Expected critic verdict | What is exercised |
| --- | --- | --- |
| EGFR L858R | `inconsistent` | The investigator cites the first synthetic item. The critic finds the second item's different direction. |
| EGFR exon 19 deletion | `no-internal-issues-found` | A single scoped item matches the proposed annotation. |
| EGFR exon 20 insertion | `insufficient-evidence` | Neither stage calls a worker when the snapshot contains no matching item. |

The [fixture](../tests/fixtures/civic/egfr_luad_page.json) uses invented IDs in the CIViC wire format. Its date and hashes are deterministic. Every candidate from it is labeled `synthetic: true`, including abstentions. Its generated links are placeholders. The demo adapter accepts only synthetic evidence bundles.

### GraphQL scope, rights and provenance

The connector uses one fixed, anonymous, read-only GraphQL query at `https://civicdb.org/api/graphql`. It requires the checked-in source's exact endpoint, open access, verified-permissive reuse, CC0-1.0 and terms URL before any request and again when loading a snapshot. CIViC is the source of the curated annotations. The reuse decision does not extend to linked publications. See the [source catalog](../catalog/sources/civic.json), [CIViC API reference](https://griffithlab.github.io/civic-v2/), [rate limits](https://docs.civicdb.org/en/latest/api.html) and [licensing](https://docs.civicdb.org/en/latest/about.html).

The server filters for accepted items, EGFR profile names and lung adenocarcinoma. It sorts by evidence ID. The connector then requires the exact LUAD label and DOID 3910, a non-deprecated profile and disease, an unflagged accepted item, and an exact, case-insensitive match to one of the three named profiles. Compound profiles, broader lung diagnoses and alternative spellings can fall outside this definition. One item represents a curated annotation. Multiple items can refer to the same study. The grounded workload can compare annotations but cannot answer patient prevalence, treatment efficacy, or cross-database questions.

| Location in a snapshot | Provenance and meaning |
| --- | --- |
| Snapshot | Source ID, endpoint, CC0 terms, UTC retrieval start, status, exact GraphQL text and bounds. The snapshot SHA-256 covers its canonical JSON excluding its own hash. |
| `pages[i]` | The exact request variables, response byte SHA-256, page size actually returned, server total and continuation flag. Combine `query` and `variables` to reconstruct the request. |
| `records[i]` | Source ID, `EID` identifier, permanent resource URL, page index, exact-scope profile/disease and curated categorical annotations. Its SHA-256 covers that projection excluding its own hash and page index, so it survives page movement. |
| Unavailable snapshot | A fixed failure code, failed request variables, completed page metadata, and zero usable records. Never interpret an unavailable response as a zero-record search. |

The connector persists only the allowlisted projection, page metadata and digests. Narratives, study texts, patients, submitters, credentials, error bodies and arbitrary upstream fields have no destination in the snapshot. Store snapshots and grounded runs under the ignored `.compute/` tree or outside the repository. Retrieve once per investigation and distribute the same snapshot to workers or shards. Hashes check consistency with an operator-trusted snapshot. They do not authenticate a file supplied by someone else. Verify the digest over a trusted channel when transferring a snapshot.

Default retrieval is bounded to four pages of 20 items. `--max-pages` accepts 1–8 and `--page-size` accepts 1–25. A page has a 128 KiB response cap and an eight-second socket timeout. The connector paces pages against the anonymous limit and permits up to three attempts per page for transient failures. It respects `Retry-After` or abandons the request when the required wait exceeds its ten-second backoff budget. It rejects redirects, malformed or partial GraphQL responses, unexpected annotation values and inconsistent pagination. If total counts change or IDs repeat/go backwards, it discards all usable records for that retrieval. Same-count edits during pagination remain possible because the API serves mutable records across separate requests.

- `complete` means the filtered API connection ended normally. Local scope exclusions still apply.
- `truncated` means the configured page budget ended first. Retained records can supply a counterexample. The critic cannot return a passing verdict on that retrieval.
- `unavailable` is written as an inspectable snapshot and `retrieve` exits 1. A subsequent run using it records `insufficient-evidence` for every job, with zero worker attempts. Malformed local inputs, denied permissions and incompatible configuration exit 2. Retry retrieval into a new path, or rerun computation against the previously pinned snapshot.

### Fan-out worker JSON protocol

Grounded requests use outer `protocol: "open-cancer-compute-v2"`. The nested deterministic v1 job and its hash stay intact. Both stages receive `protocol`, `stage`, `job`, and the same `evidence_bundle` over stdin. The `critique` request also includes `candidate`. The bundled model adapter validates the entire bundle, selects only the job's exact profile, and passes those projected records to the model. Grounded runs currently accept the CIViC-only EGFR/LUAD workload with the `alteration` dimension.

| Response | Required fields and checks |
| --- | --- |
| Investigator | `outcome`, `observation`, `evidence`, `limitations`, `synthetic`. A `record-observation` tests whether retrieved items of one `evidence_type` agree on a `field` and `value`. The field is `evidence_level`, `evidence_direction` or `significance`. At least one cited item must match. An abstention uses `outcome: "insufficient-evidence"`, a null observation and empty evidence. |
| Citation in either stage | Exactly `record_id`, `sha256`, `relation`. ID and hash must match a retrieved item in the job. For the selected type, a matching value `supports` the annotation predicate and a different value `contradicts` it. Other evidence types supply `context`. These relations describe the candidate predicate and are distinct from CIViC's own evidence direction. |
| Critic | `verdict`, `issues`, `evidence`. If any scoped item counters the predicate, the critic must report `inconsistent`, identify `annotation-counterexample`, and cite one. A `no-internal-issues-found` verdict requires a supporting citation, complete retrieval and agreement with every relevant returned item. Missing records produce `insufficient-evidence`. |

Both responses use a closed vocabulary for issue and limitation codes. The coordinator validates third-party workers as well as the bundled adapter. Unrecognized IDs, misclassified support, wrong fingerprints, free-text claims and invalid verdicts fail the attempt. A failed job records a bounded error code. Inspect the [validator and exact worker instructions](../scripts/grounded_claims.py) and the [offline protocol tests](../tests/test_grounded_claims.py).

### Existing directory-only jobs (v1)

Omit `--evidence` to use `open-cancer-compute-v1`. Existing deterministic workloads, worker commands, job hashes and JSON field shapes remain supported. The original nine-job smoke test runs as follows.

```sh
python scripts/compute_fanout.py plan --spec workloads/egfr-luad-compute.json \
  --output .compute/jobs.jsonl --max-jobs 100
python scripts/compute_fanout.py run --jobs .compute/jobs.jsonl --output .compute/demo \
  --investigator 'python scripts/demo_compute_worker.py' \
  --critic 'python scripts/demo_compute_worker.py' --workers 4 --max-jobs 100
```

The v1 investigator returns `summary`, `evidence` (source ID, URL, relation), `limitations` and `synthetic`. Only catalog homepages or documentation URLs with the relation `context` are accepted. The v1 critic returns `verdict` and `issues`. A directory-only pass becomes `insufficient-evidence`. An explicit `inconsistent` verdict keeps its issues. Legacy custom workers and resumed checkpoints must meet these citation checks. Use v2 to cite individual items.

### Fan-out checkpoints, sharding and publication

The plan contains content-addressed questions. Use `--max-jobs` to bound expansion and execution, `--workers` to bound subprocess concurrency, and `--retries` to limit repeat attempts. For N shards copy the same plan and snapshot, use `--shard-count N` and distinct `--shard-index` values from 0 to N−1, and assign each shard its own output directory.

`run-config.json` pins the plan, sharding, worker command fingerprint, timeout, retries, and for v2, the validated snapshot and a fingerprint of the adapter's endpoint/model settings. Pin executable and model versions operationally when reproducing a run. `records/<job-id>.json` contains the candidate, critique, attempts, and job/config/evidence hashes. `summary.json` reports completed and failed jobs separately from triage. `--resume` verifies matching configuration, snapshot permissions and checkpoint hashes, rechecks grounded citations, skips completed jobs and retries failed jobs. Refreshing mutable source records requires a new snapshot and run directory.

All outputs retain `review_status: unreviewed`, `publication_status: not-published`, and `clinical_validation_status: not-assessed`. A critic's pass establishes internal agreement within the scoped retrieval. Scientific interpretation, separate reproduction, rights review and any deliberate publication follow the [data policy](../DATA_POLICY.md) and [agent contract](../AGENTS.md). The Explorer does not ingest `.compute/`.
