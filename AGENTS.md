# AI compute contract

Open Cancer's engine assumes access to substantial AI compute, **not** labs, recruited scientists, reviewer committees or physical experiments. Use `scripts/research_pilot.py` for the record-grounded EGFR/LUAD benchmark and `scripts/compute_fanout.py` for generic sharded workloads. [Run commands and protocols](docs/compute.md)

## Agent execution

1. Read `DATA_POLICY.md`, the target `tasks/<task-id>.json`, and the relevant catalog source records.
2. For the nine-job EGFR/LUAD benchmark use `research_pilot.py run --mode real` or `--mode fixture`. Pin the projected public evidence once. A supplied native CIViC bundle or an existing pilot snapshot can be passed with `--snapshot`. Keep source-specific gaps, truncation, and denominators visible. Historical GDC and cBioPortal analyses do not establish that this invocation examined records from those sources.
3. For distributed workloads use `compute_fanout.py plan` with a bounded job count, separate investigator and critic subprocesses, and distinct shard indices. The three-job CIViC annotation workload runs with a snapshot from `compute_fanout.py retrieve --source civic` and `run --evidence`. Omitting evidence selects the directory-only v1 protocol. The demo worker is synthetic. Real adapters use operator-controlled endpoints and environment-only credentials.
4. In grounded fan-out v2, citations must match a job-scoped record ID and SHA-256, with a checked relation to the annotation predicate. Missing evidence requires an explicit abstention. The nine-job benchmark retains free-text observations, assumptions, limitations, citations and potential counterevidence. Both roles see the same job packet. Critic issues and deterministic checks have separate meanings. Preserve negative evidence and explicit synthetic labels. Agreement between AI models is not independent biological validation.
5. Keep job checkpoints, failed attempts and summaries in an untracked `.compute/` directory. `--resume` checks plan, worker configuration and checkpoint SHA-256. Never treat a checkpoint as a checked-in or publicly promoted result.
6. If a candidate is useful, independently reproduce its computational method, verify source provenance and reuse permissions, then **explicitly** package any permitted derived outputs under `results/<task-id>/<run-id>/` with `schemas/result.schema.json`. Findings remain `review_status: "unreviewed"` and must not imply patient treatment advice.
7. Run `python scripts/validate.py`, `python -m unittest discover -s tests -v`, `node --test tests/*.test.*` and relevant CI before merging code changes.

## Non-negotiable boundaries

No patient-level data, controlled-access records, private tokens, re-identification, prohibited scraping, invented citations, or fabricated experimental results. Preserve negative outcomes, disputed interpretations and uncertainty. External scientific review can happen separately if someone chooses to do it, but **no external reviewer is needed to run the AI compute pipeline.**

A reproducible script, cryptographic hash, AI critique or passed CI check proves neither source truth nor clinical validity. Do not use ``verified``, ``clinical``, or ``approved`` labels for unassessed AI hypotheses.
