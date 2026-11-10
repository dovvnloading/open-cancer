# Agent contribution contract (v0)

AI agents can help locate, normalize, and analyze data, but their outputs need the same verification as any research contribution. No agent is automatically a scientific reviewer.

## Required loop

1. Read `README.md`, `DATA_POLICY.md`, and the target `tasks/<task-id>.json`.
2. Use only allowed, accessible upstream resources. Track source ID, source URL, exact query/filter, release if stated, and access date. Honor rate limits.
3. Write reproducible code and small permitted artifacts in `results/<task-id>/<run-id>/`. The run ID should be stable and unique (for example, `agent-20261110-a`).
4. Create `manifest.json` conforming to `schemas/result.schema.json`: input provenance, execution command/environment, artifact paths and SHA-256, findings with record-level evidence, uncertainty, limitations.
5. Mark **every finding** `review_status: "unreviewed"`. Distinguish observations from hypotheses. Include contradicting evidence if discovered; never fabricate missing values, citations, or computational results.
6. Run `python scripts/validate.py` and unit tests. Open a pull request for human and, where appropriate, domain-expert review. **Never auto-merge scientific outputs.**

## Minimal file layout

```text
results/
  egfr-luad-evidence-coverage/
    agent-20261110-a/
      manifest.json
      analyze.py
      records.csv
```

All paths inside `manifest.json` are relative to its directory. The validator checks local files and artifact SHA-256 digests. Commands and environment descriptions should be sufficient for a second contributor to reproduce the result.

## What agents must not do

Do not use patient-level datasets, evade consent or authorization, publish private data or credentials, issue clinical treatment recommendations, or turn absence of a database record into evidence of biological absence.

A valid JSON manifest only demonstrates **structural** compliance. It is not evidence that its citations, code, claims, or methodology are scientifically correct.
