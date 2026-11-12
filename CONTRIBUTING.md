# Contributing to Open Cancer

People in biology, oncology, data engineering, software, clinical research, security, and open science are welcome. Small, auditable submissions are more useful than unsourced bulk dumps.

## Where to participate

**General introductions and reviewer invitations are not bugs.** See the [community channel guide](docs/community/README.md) and [public reviewer invitation](docs/community/independent-reviewers.md). GitHub Discussions is **currently disabled**, so we maintain the invitation in version-controlled documentation instead of leaving a recruitment issue permanently open. If Discussions is enabled later, maintainers can link to an actual thread.

**GitHub Issues are exclusively for actionable engineering, code, infrastructure, security, technical QA or reproducible research work.** Every issue must name a technical/research problem, evidence or reproduction method, a concrete deliverable, and objective acceptance criteria. The issue chooser offers required [engineering](.github/ISSUE_TEMPLATE/engineering.yml) and [research](.github/ISSUE_TEMPLATE/research.yml) forms rather than blank issues. **Do not create issues** for announcements, recruitment, introductions, general questions, broad discussion, routine monitoring output or status updates. Direct those to the [community guide](docs/community/README.md), [reviewer invitation](docs/community/independent-reviewers.md) or GitHub Actions reports. GitHub Discussions is currently disabled; no Discussion is implied to exist. Maintainers should redirect and close any non-actionable issue as not planned while keeping real research deliverables open. The EGFR/LUAD research workstream is [#3](https://github.com/dovvnloading/open-cancer/issues/3); focused work belongs in [CIViC audit #34](https://github.com/dovvnloading/open-cancer/issues/34), [GDC/cBioPortal cohort audit #35](https://github.com/dovvnloading/open-cancer/issues/35) or [independent expert assessment #36](https://github.com/dovvnloading/open-cancer/issues/36). A comment saying someone may volunteer is not proof of an independent scientific review. Formal findings should be submitted through a review-bound PR.

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

For structured reviews, follow [REVIEW_PROTOCOL.md](REVIEW_PROTOCOL.md) and `schemas/review.schema.json`. Only submit review attestations bound to existing research manifests; reviewer credentials and conflict-of-interest statements need independent human verification. Synthetic examples in `examples/reviews/` are not endorsements.

Automated schema validation is a minimum quality gate, not scientific peer review. A domain expert should verify study suitability, unit of analysis, selection bias, nomenclature, statistical assumptions, and the cited source records. Prefer multiple independent reviewers for clinical claims. A merged pull request is not an endorsement of clinical validity.

Respect the [data policy](DATA_POLICY.md), upstream attribution, and contributors. If you spot an error, open an issue with a source citation and a reproducible explanation — never include personal health information in issues.
