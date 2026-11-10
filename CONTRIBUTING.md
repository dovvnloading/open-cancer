# Contributing to Open Cancer

People in biology, oncology, data engineering, software, clinical research, security, and open science are welcome. Small, auditable submissions are more useful than unsourced bulk dumps.

## Ways to contribute

**1. Improve the source catalog.** Add one JSON file in `catalog/sources/` matching `schemas/source.schema.json`. Provide an official provenance URL, access conditions, source modalities, a review date, and reuse details. If the redistribution license is not confirmed, use `verify-upstream`. Never treat a publicly readable website as a blanket open-data license.

**2. Define a research question.** Add `tasks/<task-id>.json` matching the task schema. Give a narrow research question, upstream sources, concrete deliverables, quality checks, and non-goals. Favor falsifiable questions, negative findings, and analyses that can be independently reproduced.

**3. Submit reproducible results.** Under `results/<task-id>/<run-id>/`, include `manifest.json` matching `schemas/result.schema.json`, executable code, and one or more small, non-sensitive artifacts. Each input needs a source ID, date, URL, and query/filter. Each claimed finding needs record-level evidence (including contradictory evidence where known), limitations, and `review_status: "unreviewed"`. Record code environment and command; include SHA-256 for each artifact. See [AGENTS.md](AGENTS.md).

## Pull request checklist

- [ ] Every source, finding, and derived artifact is traceable to its origin, version/date, and methods.
- [ ] All files are appropriate to redistribute; source data licenses and attribution have been checked.
- [ ] No patient information, sensitive or access-restricted data, secrets, or large raw dumps.
- [ ] Methods can be rerun; uncertainty, missing data, and failed/negative results are disclosed.
- [ ] Validation passes: `python scripts/validate.py` and `python -m unittest discover -s tests -v`.
- [ ] Scientific findings are explicitly **unreviewed**; no clinical recommendations are made.

## Review expectations

Automated schema validation is a minimum quality gate, not scientific peer review. A domain expert should verify study suitability, unit of analysis, selection bias, nomenclature, statistical assumptions, and the cited source records. Prefer multiple independent reviewers for clinical claims. A merged pull request is not an endorsement of clinical validity.

Respect the [data policy](DATA_POLICY.md), upstream attribution, and contributors. If you spot an error, open an issue with a source citation and a reproducible explanation — never include personal health information in issues.
