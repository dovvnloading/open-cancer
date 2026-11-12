# AI compute engine: plan → investigate → critique → triage

Open Cancer is an AI/software investigation project using abundant user-supplied compute. It does not operate a lab or presume academics or clinicians will participate. This repository supplies a **provider-neutral orchestration baseline**, not GPU infrastructure, paid model access or validated biomedical discoveries.

## Runnable local smoke test (synthetic, offline)

From the repository root (Python 3.11+):

```sh
python scripts/compute_fanout.py plan \
  --spec workloads/egfr-luad-compute.json \
  --output .compute/jobs.jsonl --max-jobs 100
python scripts/compute_fanout.py run \
  --jobs .compute/jobs.jsonl --output .compute/demo \
  --investigator 'python scripts/demo_compute_worker.py' \
  --critic 'python scripts/demo_compute_worker.py' \
  --workers 4 --max-jobs 100
cat .compute/demo/summary.json
```

The demo worker **does not call a model, retrieve source records, or produce scientific findings**. All nine records are `synthetic`, not published, unreviewed and not clinically assessed. This smoke test proves the execution protocol, not the biology.

## Scale on your AI hardware

The coordinator is intentionally model- and infrastructure-agnostic. A worker is any executable that reads one JSON object from stdin and prints one JSON object to stdout; it may call local inference engines, a GPU service or a remote endpoint. Credentials belong in secure environment variables/secrets, **never in job JSON, command arguments, logs or Git**. Run separate models/prompts for investigation and criticism when possible. Enforce source API terms and rate limits inside adapters.

A working OpenAI-compatible inference adapter is included as `scripts/openai_compatible_worker.py`. Point it at your inference server with `OPEN_CANCER_LLM_URL` (full `/v1/chat/completions` endpoint), and set `OPEN_CANCER_LLM_MODEL`, optionally `OPEN_CANCER_INVESTIGATOR_MODEL`, `OPEN_CANCER_CRITIC_MODEL`, `OPEN_CANCER_LLM_KEY`, `OPEN_CANCER_LLM_TIMEOUT` and `OPEN_CANCER_MAX_TOKENS`. Plain HTTP is allowed for local loopback services; use HTTPS remotely. Its inputs contain **only checked-in source-directory metadata**, not validated evidence records. It asks the model for source-scoped research leads and forces critiques to `insufficient-evidence` in the absence of records. It cannot establish biological findings. For example:

```sh
export OPEN_CANCER_LLM_URL='http://127.0.0.1:8000/v1/chat/completions'
export OPEN_CANCER_LLM_MODEL='your-local-model-id'
python scripts/compute_fanout.py run --jobs .compute/jobs.jsonl \
  --output .compute/live --investigator 'python scripts/openai_compatible_worker.py' \
  --critic 'python scripts/openai_compatible_worker.py' --workers 8 --max-jobs 100
```

For N machines or partitions, copy the **same immutable jobs.jsonl** to each, use the same worker commands/options, and run each with `--shard-count N --shard-index 0..N-1`, storing each shard in its own output directory. IDs are content-addressed and partitioned by a stable hash, so shards do not overlap. Use `--workers` to limit concurrent subprocesses per machine; `--max-jobs` bounds selected jobs, and `--timeout` and `--retries` limit worker spending. There is no central cluster scheduler or built-in GPU allocator. `--resume` skips completed checkpoints only when the plan and worker configuration match and the checkpoint SHA-256 is intact; failed jobs are retried.

Example for one shard (replace both adapter commands with executables you actually operate):

```sh
python scripts/compute_fanout.py run \
  --jobs .compute/jobs.jsonl --output .compute/shard-2-of-8 \
  --investigator 'python my_workers/investigate.py' \
  --critic 'python my_workers/critique.py' \
  --workers 32 --timeout 180 --retries 1 --max-jobs 5000 \
  --shard-count 8 --shard-index 2 --resume
```

Omit `--resume` for the initial run; pass it when re-running the same configuration. A workload is an allowlisted checked-in `task_id`, associated `source_ids` and a matrix of question dimensions (`axes`). The `plan` command refuses an explosion beyond `--max-jobs`, produces a deterministic JSONL with derived job hashes, and does not consume inference compute.

## Worker JSON protocol, v1

Each subprocess receives one stdin JSON object:

```json
{"protocol":"open-cancer-compute-v1", "stage":"investigate", "job":{"protocol":"open-cancer-compute-v1","id":"...","task_id":"...","question":"...","source_ids":["civic"],"dimensions":{"alteration":"EGFR L858R"}}}
```

The investigator returns exactly:

```json
{"summary":"A source-scoped computational observation of sufficient length.","evidence":[{"source_id":"civic","url":"https://civicdb.org/","relation":"context"}],"limitations":["This URL is contextual; no record-level claim has been established."],"synthetic":false}
```

The critic is called only after a structurally valid investigator result, with the same `protocol`, `job`, `stage: "critique"`, and a `candidate` field containing the investigator response. The critic returns exactly:

```json
{"verdict":"insufficient-evidence","issues":["The context-only URL does not establish the claimed biological relationship."]}
```

Critic verdicts: `no-internal-issues-found`, `inconsistent`, `insufficient-evidence`. The first is **not** truth, peer review or clinical validation. `evidence` requires task-listed source IDs and HTTPS links, but the coordinator cannot certify their truth, licensing, existence or interpretation. Worker responses are size-limited. Worker failures, invalid JSON, unsupported fields and timeouts fail the job closed and remain visible in its record. The coordinator does not execute source text as code or invoke shell expansion.

## Outputs and what they mean

`run-config.json` pins the plan hash, sharding, deadlines/retries and worker command fingerprint. `records/<job-id>.json` has a SHA-256 envelope, the exact job, candidate/critique or failure, and `review_status: unreviewed`, `publication_status: not-published`, `clinical_validation_status: not-assessed`. `summary.json` lists every selected job, completed/failed totals and triage counts. Triage is **workflow prioritization only**; even `candidate-unverified` is not biological or clinical evidence.

None of these files is automatically added to the Explorer or public API. A separate, intentional source/rights/computational verification step is needed before optionally committing non-sensitive derived results using `schemas/result.schema.json`. Running more models can surface or challenge more hypotheses, but consensus alone does not confirm them. Never upload controlled/patient-level input bytes or publish AI-generated treatment recommendations.

## What to build next

Permission-aware record retrieval for the supplied model adapter, source-specific retrieval caches with permission checks, compute/cost/throughput metrics, automated recomputation of *numeric* claims from pinned permitted inputs, and comparison benchmarks belong in engineering PRs—not recruitment or governance tickets. The initial pipeline supplies job orchestration and provenance, not those integrations yet.
