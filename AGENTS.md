# AI compute contract

Open Cancer's engine assumes access to substantial AI compute, **not** labs, recruited scientists, reviewer committees or physical experiments. The runnable coordinator is `scripts/compute_fanout.py` and its protocol is documented in [docs/compute.md](docs/compute.md).

## Agent execution

1. Read `DATA_POLICY.md`, the target `tasks/<task-id>.json`, and the relevant catalog source records.
2. Use `compute_fanout.py plan` to expand a bounded, deterministic workload; use distinct `--shard-index` values across machines. Plans contain questions, never asserted outcomes.
3. Supply a separate `--investigator` and `--critic` CLI adapter. Both read a protocol JSON request on stdin and emit one JSON response on stdout; run workers on your own compute. The demo adapter is **synthetic** and not an AI model.
4. Every investigator candidate needs linked public evidence/context, limitations and an explicit `synthetic` flag. Critics report internal inconsistencies or insufficient evidence; agreement between AI models is not independent biological validation.
5. Keep job checkpoints, failed attempts and summaries in an untracked `.compute/` directory. `--resume` checks plan, worker configuration and checkpoint SHA-256. Never treat a checkpoint as a checked-in or publicly promoted result.
6. If a candidate is useful, independently reproduce its computational method, verify source provenance and reuse permissions, then **explicitly** package any permitted derived outputs under `results/<task-id>/<run-id>/` with `schemas/result.schema.json`. Findings remain `review_status: "unreviewed"` and must not imply patient treatment advice.
7. Run `python scripts/validate.py`, `python -m unittest discover -s tests -v` and relevant CI before merging code changes.

## Non-negotiable boundaries

No patient-level data, controlled-access records, private tokens, re-identification, prohibited scraping, invented citations, or fabricated experimental results. Preserve negative outcomes, disputed interpretations and uncertainty. External scientific review can happen separately if someone chooses to do it, but **no external reviewer is needed to run the AI compute pipeline.**

A reproducible script, cryptographic hash, AI critique or passed CI check proves neither source truth nor clinical validity. Do not use ``verified``, ``clinical``, or ``approved`` labels for unassessed AI hypotheses.
