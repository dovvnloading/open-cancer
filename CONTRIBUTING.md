# Contributing to Open Cancer

People in biology, oncology, data engineering, software, clinical research, security, and open science are welcome. Small, auditable submissions are more useful than unsourced bulk dumps.

## Where to participate

**General introductions and reviewer invitations are not bugs.** See the [community channel guide](docs/community/README.md) and [public reviewer invitation](docs/community/independent-reviewers.md). GitHub Discussions is **currently disabled**, so we maintain the invitation in version-controlled documentation instead of leaving a recruitment issue permanently open. If Discussions is enabled later, maintainers can link to an actual thread.

**GitHub Issues are reserved for actionable engineering work:** reproducible
software defects, code/implementation changes, CI, infrastructure, security,
data-pipeline failures and technical production incidents. Each ticket needs
the affected component, reproducible evidence or a precise implementation gap,
and verifiable acceptance criteria. The issue chooser has an
[engineering-only form](.github/ISSUE_TEMPLATE/engineering.yml) and disables
blank issues.

**Do not file GitHub Issues for research coordination, scientific audits
without a specific code defect, independent expert reviews, research
milestones, invitations, general questions, announcements or status reports.**
These belong in [version-controlled research work](docs/research/egfr-luad-open-work.md),
the [reviewer guide](docs/community/independent-reviewers.md), research
result/review PRs, or [community guidance](docs/community/README.md).
Routine checks belong in GitHub Actions. GitHub Discussions is currently
disabled; do not claim a Discussions thread exists.

If scientific work identifies an actual reproducible code defect, open
an **engineering issue about that defect** with reproduction and expected fix,
not an issue to coordinate the research itself. Former coordination issues
#3/#34/#35/#36 are historical, closed trackers; their research has not been
declared complete.

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
