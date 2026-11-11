# Independent reviewer package: EGFR / lung adenocarcinoma

**Purpose:** support reproducible critique of the Open Cancer public-data pilot. [Public independent reviewer invitation](../../docs/community/independent-reviewers.md) · [focused research review tasks](../../docs/community/README.md). **Status: not independently scientifically or clinically validated.**

[Historical Explorer demo — may be stale](https://open-cancer-explorer.lagging.chatgpt.site/?view=pilot) · [research task](../../tasks/egfr-luad-evidence-coverage.json) · [original run](../../results/egfr-luad-evidence-coverage/pilot-20261110/README.md) · [review protocol](../../REVIEW_PROTOCOL.md) · [issue #3](https://github.com/dovvnloading/open-cancer/issues/3)

## Recorded values (November 10, 2026)

| Source | Value | Exact unit |
| --- | ---: | --- |
| CIViC | 56 | Editor-accepted EGFR-token **evidence items** with exact lung adenocarcinoma label (52 predictive, 2 prognostic, 2 oncogenic) |
| cBioPortal | 86 | EGFR **mutation-call API rows** from TCGA PanCancer Atlas LUAD, not distinct patients |
| GDC | 93 | EGFR **simple somatic mutation occurrence rows** in TCGA-LUAD |
| GDC | 78 | **Distinct GDC case identifiers** among those 93 occurrence rows |

**The measurements have different observation units, study selections, and no harmonized denominators.** They cannot be pooled as prevalence or used to infer clinical efficacy or treatment recommendations. The original research manifest has no asserted scientific findings (`findings: []`).

## Build a reproducible offline review bundle

From a new clone of the repository, using Python 3.12+:

~~~bash
git clone https://github.com/dovvnloading/open-cancer.git
cd open-cancer
python -m pip install -r requirements-dev.txt
python scripts/build_review_packet.py \
  --output /tmp/open-cancer-egfr-review \
  --archive /tmp/open-cancer-egfr-review.zip
cd /tmp/open-cancer-egfr-review
sha256sum --check SHA256SUMS
~~~

The builder verifies the checked-in result manifest and its artifact hashes, recomputes CSV summary totals, and packages **only permitted, non-sensitive derived artifacts**. It generates a targeted `evidence_review_sample.csv` of 12 public CIViC evidence IDs (8 predictive, 2 prognostic, 2 oncogenic), selected deterministically from the original manifest hash. This sample oversamples uncommon categories for quality checks and **is not statistically representative**.

Package contents: the original manifest, coverage summary, CIViC evidence ID table, aggregated cBioPortal mutation table, an independent-review QA sample, analysis entrypoint (for inspection only), guide, checklist, `packet-index.json`, `SHA256SUMS` and source references. It **does not contain patient or case identifiers, controlled-access data, raw upstream mutation records or clinical conclusions**. The bundled `analyze.py` requires a full repository checkout to run.

The manifest SHA-256 identifies the **exact original reviewed snapshot**, even when the public portal's live release changes. The bundle and ZIP are deterministic across identical checkouts. The builder makes **no external network calls**.

## Durable release and provenance status

The 90-day CI packet artifact is for build inspection, **not permanent distribution**. Versioned, checksum-verified public releases are managed by the [reviewer release runbook](../../docs/reviewer-release.md). No release is claimed to exist merely because a workflow or candidate URL exists. A maintainer must publish an explicitly reviewed tag, and set the pinned-release monitor variables after successful verification.

`packet-index.json` now records every original upstream input URL, query/filter, accessed date, response hash and source identifier, plus separate `method_reproduction`, `source_version_drift` and `exact_input_replay` statuses. GDC and cBioPortal accessions and CIViC nightly export names are **not** asserted to be immutable raw-data releases. Original raw source bytes are not bundled, so exact input replay cannot be guaranteed. The original review-bound snapshot and manifest are left untouched.

## Independently rerun the computation

Run from the full repository checkout, with outbound access to the documented public source APIs. Do **not** overwrite the original, review-bound run:

~~~bash
run_id="external-$(date -u +%Y%m%d%H%M%S)"
folder="results/egfr-luad-evidence-coverage/$run_id"
mkdir "$folder"
cp results/egfr-luad-evidence-coverage/pilot-20261110/analyze.py "$folder/analyze.py"
python "$folder/analyze.py" --out-dir "$folder" --run-id "$run_id"
python scripts/validate.py
~~~

Compare `source_content_sha256` from the new run with the original. If upstream versions differ, **report source drift** rather than claiming byte-identical reproduction. If the optional GDC endpoint fails, it must be recorded as unavailable/null, **never zero**.

## How to submit your review

1. Follow [SCIENTIFIC_CHECKLIST.md](SCIENTIFIC_CHECKLIST.md). Inspect each sampled CIViC record at its original source URL; check variant/disease semantics, inclusion errors and exclusions. Negative findings are welcome.
2. Decide which aspects you actually reviewed: metadata, methods or scientific interpretation. Record passed/concern/not-assessed, original evidence links, your reproducibility command, detailed limitations, affiliation/expertise and conflicts.
3. Submit a pull request with one JSON record under `reviews/egfr-luad-evidence-coverage/pilot-20261110/<review-id>.json`, conforming to [the review schema](../../schemas/review.schema.json). The `reviewed_manifest_sha256` must equal the **original manifest's exact byte hash**, as reported in the bundle. See [fictional examples](../../examples/reviews) for shape only; replace every fictitious field and use `synthetic:false` only for real reviews.
4. Do not submit a scientific interpretation assessment as an AI agent. Maintainers must verify self-declared human credentials independently. An expert interpretation assessment is **not clinical approval**; `clinical_validity` must remain `not-assessed`.

## Known scientific review risks

- **DOID:3910 exact-match selection** excludes broader NSCLC and lung cancer evidence labels that may include adenocarcinoma.
- CIViC `Accepted` indicates editor-reviewed **curation**, not proof of biological causality or therapeutic benefit. An EGFR token can occur in a **compound molecular profile**; 56 evidence items are not 56 unique variants, papers or trials.
- The cBioPortal count is mutation-call rows from a selected study/profile and sample list; **no fully profiled denominator** was calculated.
- GDC is queried via public `/ssm_occurrences` by project and gene and case UUIDs are deduplicated in memory only. **Not** a driver-mutation prevalence estimate or a cohort-harmonized match to cBioPortal.
- Public source releases are mutable. The source fingerprints are recorded but complete raw third-party responses are not archived here.
- **Known documentation inconsistency:** the original immutable `coverage.gdc.note` mentions a legacy "aggregate endpoint". The actual recorded retrieval uses `/ssm_occurrences`. Treat the `manifest.inputs[]` URLs and reproducibility code as authoritative; do not silently modify the review-bound snapshot.

**Completion:** at least one genuinely external qualified cancer genomics expert submits a version-bound, evidence-specific review, including flaws/limitations; ideally a second reviews before promoting scientific interpretations. Track its formal acceptance in [scientific assessment #36](https://github.com/dovvnloading/open-cancer/issues/36), not as a permanent recruitment issue. The [public invitation](../../docs/community/independent-reviewers.md) is a request for critique, not an endorsement.
