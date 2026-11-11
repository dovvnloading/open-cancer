# Scientific review and correction protocol

Open Cancer keeps a **separation of duties** between automated validation, data engineering, and biomedical interpretation.

1. **Format checks** (JSON Schema, known source IDs, file SHA-256) demonstrate structural validity; they **cannot** verify a clinical or biological claim.
2. **Methods reproduction** means a reviewer independently ran the documented code on matching or declared source versions. It is narrower than scientific validity.
3. **Scientific interpretation assessment** requires an independent human with relevant domain expertise, transparent credentials, a conflict-of-interest statement, and evidence-specific reasoning. A merge or badge **never** means a finding is clinically validated.
4. **Clinical use** is **not assessed** by this repository; no patient-specific treatment, diagnosis, or regulatory endorsement should result from a review record.

## Review procedure

Read the original task, result `manifest.json`, upstream citations, program, raw-response checksums, artifact checksums, source-specific reuse terms, and analysis limitations. Explicitly consider:

- **Source provenance:** original record identifiers, stable versions/access dates, evidence for every assertion, conflicting and negative results.
- **Licensing and privacy:** permissions to access *and* redistribute derived output, authorized consent, non-identifiability, no controlled data.
- **Cohorts and observational units:** tissue/stage/assay, population inclusion criteria, missingness and denominators, one evidence record versus one patient versus one mutation call.
- **Reproducibility:** environment, exact query, determinism, controls and failures, SHA-256 match of outputs.
- **Scientific reasoning:** ontology/gene/variant mappings, confounding, clinical vs preclinical vs observational evidence, absence of causal/efficacy leaps.
- **Uncertainty and limitations:** selection bias, release drift, failed upstream responses, contradiction and limitations acknowledged.

Any component not checked **must** be `not-assessed`, not a silent pass. A missing value is never automatically zero. Reviewers should use `needs-revision` when conclusions exceed available evidence.

## Submit a review

**Only review results that are already present in the repository.** Save a JSON record to

```text
reviews/<task-id>/<run-id>/<review-id>.json
```

matching [the review schema](schemas/review.schema.json). The record points to the **exact** SHA-256 of the result's `manifest.json`, the reviewed finding IDs, reproducibility command and outcomes, decision, typed checks, conflicts and outstanding work. The validator refuses stale review hashes, missing result manifests, mismatched IDs, or invented finding references.

The review schema distinguishes `methods-reproduced` from `scientific-assessment-completed`; agents can identify concerns and reproduce methods but cannot mark a scientific interpretation assessment complete. Credentials are **self-declared** in the record and must be checked by project maintainers separately. Even a submitted `scientific-assessment-completed` record indicates that a human reports an assessment with substantive source/cohort/biological/uncertainty checks; it **does not prove verified identity, independence, clinical correctness, or endorsement**. A disputed or needs-revision assessment must remain visibly contested in the Explorer. A machine-valid review is **not** an identity check, a clinical endorsement, or a peer-reviewed publication.

## Outcomes

- `needs-revision`: missing, incorrect or unclear work; no favorable scientific assertion.
- `methods-reproduced`: evidence that the computation ran as described on declared data; does not verify conclusions.
- `scientific-assessment-completed`: independent human expert **assessed** interpretation with documented caveats; not automatic approval or clinical validity.
- `disputed`: published method or interpretation is contested; preserve both original claim and objections with provenance.

Each review is **append-only in substance**: corrections should create new review records or clearly version old ones with a changelog, rather than silently rewriting conclusions. A reviewer must disclose conflicts and should not review their own scientific contribution as an independent expert. Prefer two independent reviewers for high-stakes clinical claims.

## Synthetic examples

See `examples/reviews/` for a **fictional** methods-reproduced review and a **fictional** needs-revision review. Both are clearly marked `synthetic: true` and do not correspond to real research, people or findings. They validate structure only, and must never be placed in live `reviews/`.

To report a potential problem, open a GitHub issue containing links, record IDs, and a reproducible account of the discrepancy. **Never include patient health information**. A critical safety/privacy concern should be escalated privately to maintainers before posting source data.

This project is for nonclinical research. Expert review does not imply regulatory approval or a care recommendation.
