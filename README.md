# Open Cancer

**Compute-first, AI-driven investigation of public cancer research data.**

Open Cancer is a software project built to use abundant AI compute: fan out bounded investigations over public cancer-data sources, challenge candidate findings with separate AI workers, and retain reproducible evidence and failure records. **It is not a laboratory, clinical service, academic institution, or program that depends on recruiting professors, medical reviewers, volunteers or human research staff.**

The goal is useful computational research leads and tools—not claims of cures, clinical efficacy, or experimentally validated biology. Parallel AI agreement is not scientific proof.

## Run the compute engine

Requires Python 3.11+. The fan-out coordinator uses the standard library. The nine-job benchmark and repository validation also use `requirements-dev.txt`.

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

The research pilot uses this same nine-job workload to investigate pinned public CIViC evidence items. It reads the CIViC GraphQL connector’s record bundles and can use a separately fingerprinted nightly export when a fresh live GraphQL request fails at the network or HTTP layer. Each investigator response and adversarial critique is inspectable beside its citations, deterministic checks and source gaps. GDC and cBioPortal records are outside the pilot connector. Operator-supplied OpenAI-compatible endpoints can use separate URLs, model IDs and credentials for the two roles. A synthetic fixture exercises the full workflow offline.

```sh
python -m pip install -r requirements-dev.txt
python scripts/research_pilot.py run --mode fixture --output .compute/fixture
cat .compute/fixture/report.md
```

The technical report and JSON benchmark include job failures, citation checks, evidence coverage, model configuration, execution time and available usage. Outputs remain unreviewed, unpublished and nonclinical. [Nine-job real pilot launch and offline benchmark](docs/compute.md#launch-a-real-pilot)

The three-job CIViC annotation workflow tests a narrower question. An investigator proposes an annotation predicate for one alteration, and a critic searches the same pinned items for a counterexample. The fan-out coordinator checks citation fingerprints and verdicts. It accepts custom subprocess adapters and supports sharding and resume. Use its native snapshot with either workflow.

```sh
python scripts/compute_fanout.py retrieve --source civic --output .compute/civic-live.json
python scripts/compute_fanout.py plan --spec workloads/civic-egfr-luad-annotations.json \
  --output .compute/civic-jobs.jsonl
```

With inference endpoints configured, launch the three jobs using `compute_fanout.py run --evidence .compute/civic-live.json`, or pass the same bundle to the nine-job benchmark using `research_pilot.py run --mode real --snapshot .compute/civic-live.json`. The guide supplies complete commands for both workflows. [Three-job launch and offline example](docs/compute.md#three-job-civic-annotation-workflow)

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

- Evaluate configured investigator and critic models against the same pinned EGFR/LUAD input.
- Scale bounded job execution across available machines using deterministic sharding and resume.
- Extend record retrieval to additional sources where permissions and observation units allow reproducible comparisons.
- Compare pilot failure rates, output consistency and measured token cost before increasing parallelism.
- Surface *unverified* computational candidates and contrary evidence without auto-promoting clinical or scientific claims.

**Scope:** Public computational research only, not patient-specific diagnosis or treatment. Respect source terms and privacy; see [DATA_POLICY.md](DATA_POLICY.md).
