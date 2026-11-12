# Open Cancer

**Compute-first, AI-driven investigation of public cancer research data.**

Open Cancer is a software project built to use abundant AI compute: fan out bounded investigations over public cancer-data sources, challenge candidate findings with separate AI workers, and retain reproducible evidence and failure records. **It is not a laboratory, clinical service, academic institution, or program that depends on recruiting professors, medical reviewers, volunteers or human research staff.**

The goal is useful computational research leads and tools—not claims of cures, clinical efficacy, or experimentally validated biology. Parallel AI agreement is not scientific proof.

## Run the compute engine

Requires Python 3.11+. The coordinator uses only the standard library; existing repository validation also uses `requirements-dev.txt`.

```bash
# Plan nine source-scoped questions deterministically; no model calls.
python scripts/compute_fanout.py plan \
  --spec workloads/egfr-luad-compute.json \
  --output .compute/jobs.jsonl --max-jobs 100

# Smoke test only: synthetic workers, no API calls or real AI conclusions.
python scripts/compute_fanout.py run \
  --jobs .compute/jobs.jsonl --output .compute/demo \
  --investigator 'python scripts/demo_compute_worker.py' \
  --critic 'python scripts/demo_compute_worker.py' \
  --workers 4 --max-jobs 100
```

**For real AI compute:** use the included [OpenAI-compatible endpoint adapter](scripts/openai_compatible_worker.py) with your own inference server, or replace both commands with your own worker adapters (local GPU models, remote inference, or cluster entrypoints). Example:

```sh
export OPEN_CANCER_LLM_URL='http://127.0.0.1:8000/v1/chat/completions'
export OPEN_CANCER_LLM_MODEL='your-local-model-id'
# Optionally choose a different critic model: OPEN_CANCER_CRITIC_MODEL
python scripts/compute_fanout.py run --jobs .compute/jobs.jsonl \
  --output .compute/live --investigator 'python scripts/openai_compatible_worker.py' \
  --critic 'python scripts/openai_compatible_worker.py' --workers 8 --max-jobs 100
```

The bundled model adapter currently sees **source directory metadata only**, not retrieved source records. It produces unverified research questions/data-quality leads and fails closed on invented record-level evidence; implement permission-aware source retrieval before treating it as a source-level investigation. The investigator accepts a JSON request and returns a candidate with source links; the separate critic attempts to falsify it. Split a stable plan into `--shard-count N` and run each `--shard-index` with its own output directory on different machines. Every completed job remains `unreviewed`, `not-published`, and `not-assessed` for clinical validation, regardless of AI critique outcome. The coordinator does **not** include a model provider or provision hardware. See [compute protocol, sharding, resume and failure handling](docs/compute.md).

## What the project actually does

1. **Search and collect metadata:** link to public cancer datasets, studies and source records without mirroring private or controlled data.
2. **Distribute computational questions:** generate deterministic, content-addressed jobs over selected data/analysis dimensions.
3. **Analyze and challenge:** run an investigator and a separate AI critic per job; explicitly retain contradictions, missing evidence and failed jobs.
4. **Measure rather than assume:** use offline validation, hashes, reproducible scripts and source-level citations. None of these establish medical correctness.
5. **Ship useful software:** a public Explorer, machine-readable catalog, computational workflows and inspectable results.

There is **no requirement** to recruit independent scientists, run physical experiments, or publish a review packet before doing AI research. Historical reviewer materials remain in Git history as optional context; they are not the operating model.

## Explorer and existing data

[Live Explorer](https://dovvnloading.github.io/open-cancer/) · [EGFR/LUAD example](https://dovvnloading.github.io/open-cancer/?view=pilot) · [JSON API](https://dovvnloading.github.io/open-cancer/api/v1/index.json) · [API guide](docs/explorer-api.md)

The source and dataset catalog currently references NCI GDC, cBioPortal, CIViC, ClinicalTrials.gov, GEO and DepMap. It links to upstream records with documented access/reuse caveats rather than hosting protected records. The existing EGFR/LUAD example is a computational snapshot, **not independently validated clinical evidence**. AI fan-out checkpoints do **not** enter the public Explorer or `results/` automatically.

## Code-first contributions

- Improve agent adapters, execution, source connectors, benchmarking, reproducibility, CI or the Explorer.
- Keep GitHub Issues **strictly for actionable code/engineering problems** with reproduction steps and technical acceptance criteria. No recruitment, governance, review committees, open-ended research initiatives or routine status tickets.
- For check-in of an actual computational result, follow [AGENTS.md](AGENTS.md) and [DATA_POLICY.md](DATA_POLICY.md); never publish patient-level/controlled records or AI text as verified medical facts.
- Run `python scripts/validate.py`, `python -m unittest discover -s tests -v`, and `node --test tests/*.test.*` as appropriate.

The codebase includes legacy manifests and publication metadata for backward compatibility. A manifest or API record is not an endorsement. See [CONTRIBUTING.md](CONTRIBUTING.md) and the [compute architecture](docs/compute.md).

## Near-term engineering priorities

- Connect real AI inference backends through the worker JSON protocol.
- Scale bounded job execution across available machines using deterministic sharding and resume.
- Add source-specific, rights-respecting retrieval tools, frozen input snapshots where permitted, and testable computational comparisons.
- Measure throughput, cost, failure rates, duplication, and reproducibility before increasing parallelism.
- Surface *unverified* computational candidates and contrary evidence without auto-promoting clinical or scientific claims.

**Scope:** Public computational research only, not patient-specific diagnosis or treatment. Respect source terms and privacy; see [DATA_POLICY.md](DATA_POLICY.md).
