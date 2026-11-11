# Expert checklist: EGFR/LUAD public-source coverage

**Target:** `egfr-luad-evidence-coverage` / `pilot-20261110`. Check the [original manifest](../../results/egfr-luad-evidence-coverage/pilot-20261110/manifest.json). Mark each domain as **passed / concern / not-assessed** with methods and source links. Do **not** infer treatment recommendations.

## 1. Data integrity and provenance

- [ ] Verify the original manifest SHA-256 and all committed artifact SHA-256 values.
- [ ] Confirm the original manifest includes no endorsed findings and is labeled unreviewed.
- [ ] Run a separate, dated live computation or explain why you could not. Compare source-response hashes and describe any upstream release drift.
- [ ] Check CIViC licensing, the different data access/reuse terms at GDC and cBioPortal, and that no patient, sample or case-level data are published.

## 2. Disease and variant semantics

- [ ] Validate EGFR gene identifiers NCBI Entrez 1956 and Ensembl ENSG00000146648 and distinguish symbols from variant-level mappings.
- [ ] Confirm whether DOID:3910 is an appropriate strict LUAD filter and document broader NSCLC evidence that this method excludes.
- [ ] Follow **every** public record link in `evidence_review_sample.csv`. Inspect disease, molecular profile (including compound profiles), evidence type and reused publications.
- [ ] Explain why CIViC editor `Accepted` evidence is not an independent demonstration of therapeutic benefit.
- [ ] Evaluate cBioPortal study ID, profile/sample-list filters, mutation call annotation, missingness and lack of a sequenced-cohort denominator.
- [ ] Evaluate the GDC SSM occurrence query, consequence and cohort semantics, deduplication, page completeness and lack of a profiled denominator.

## 3. Interpretation and uncertainty

- [ ] Explain why 56 CIViC items, 86 cBioPortal mutation-call rows, 93 GDC SSM occurrence rows and 78 GDC distinct matching case IDs are **not directly comparable**.
- [ ] Document TCGA cohort overlap, different pipelines/assays, ambiguous or missing variants and possible source sampling bias.
- [ ] Identify contradictions/negative evidence and what the current analysis does **not** establish.
- [ ] Inspect the legacy `coverage.gdc.note` wording; record it as a documentation caveat.
- [ ] List specific reproducible changes that would address each scientific concern.

## 4. Public review submission

Choose scope (`metadata` / `methods` / `scientific-interpretation`) and a decision (`needs-revision`, `methods-reproduced`, `scientific-assessment-completed`, `disputed`) as defined in [REVIEW_PROTOCOL.md](../../REVIEW_PROTOCOL.md). For every decision record your independent identity/credentials, conflicts, exact manifest SHA, command/outcome, all six required checks and detailed rationale.

Submit a GitHub PR to `reviews/egfr-luad-evidence-coverage/pilot-20261110/` and link [issue #3](https://github.com/dovvnloading/open-cancer/issues/3). **An agent cannot attest scientific interpretation.** Project-level review is not independent review, and expert assessment is not clinical validation. Mark unassessed items explicitly; accepting critical/negative results is essential.
