# EGFR/LUAD research and scientific review work

**Tracking location:** version-controlled research documentation and review pull requests,
**not GitHub Issues**. This page replaces the former research coordination tickets
[#3](https://github.com/dovvnloading/open-cancer/issues/3),
[#34](https://github.com/dovvnloading/open-cancer/issues/34),
[#35](https://github.com/dovvnloading/open-cancer/issues/35) and
[#36](https://github.com/dovvnloading/open-cancer/issues/36).
Those issues are archived **because they were filed in the wrong tracker, not
because the scientific work was completed**. This page is the current work queue.
GitHub Discussions is not presently enabled; contributions can be proposed in
research result or review pull requests as described below.

**Research task:** [EGFR/LUAD evidence coverage](../../tasks/egfr-luad-evidence-coverage.json)
· **Original run:** [pilot-20261110](../../results/egfr-luad-evidence-coverage/pilot-20261110/README.md)
· **Reviewer kit:** [egfr-luad-v1](../../review_packets/egfr-luad-v1/README.md)
· **Protocol:** [REVIEW_PROTOCOL.md](../../REVIEW_PROTOCOL.md)

**Frozen original manifest SHA-256:**
`9fce9beaaf873b5f0b0a14e338c3b372d8aa5d19c03be11b911ea1eb2b92931e`
(November 10, 2026). Do not rewrite the manifest or artifacts: reviews are
bound to those exact bytes. The original work has internally reproduced
methods but **no independent scientific interpretation approval** and no
clinical-validation status.

The historical observations were **56 editor-accepted CIViC evidence items**,
**86 cBioPortal EGFR mutation-call rows**, **93 GDC simple somatic mutation
occurrence rows**, and **78 distinct GDC case IDs among those occurrences**.
These are non-equivalent units, not unique patients across sources, treatment
outcomes, prevalence estimates or validated biological findings. The original
manifest asserts `findings: []`.

## Research A: CIViC evidence, ontology and variant semantics

**Status: open; not independently assessed.** Historical coordination:
[archived #34](https://github.com/dovvnloading/open-cancer/issues/34).

Research question: Does strict EGFR/LUAD evidence selection correctly
distinguish lung adenocarcinoma from broader NSCLC/lung disease terms, handle
compound molecular profiles, and retain negative/contradictory evidence?

- [ ] Inspect all **12 deterministically sampled** public CIViC records
  against source URLs and primary citations. Record EID, retrieval date,
  source version, record status, inclusion decision and rationale. This
  sample is deliberately **not statistically representative**.
- [ ] Evaluate DOID:3910 exact matching versus broader NSCLC/lung labels,
  EGFR-token false positives, compound profiles, fusion partners,
  categorical labels, repeated publications, duplicate records and missing,
  negative or contradictory evidence.
- [ ] Document false inclusions/exclusions with traceable records, or
  explicitly mark what could not be assessed. Historical excluded
  broader-lung EIDs were not individually preserved.
- [ ] Compare any newly retrieved public-source fingerprints to the
  immutable original; report changed bytes/source release drift rather
  than claiming byte-identical source reproduction.
- [ ] Publish a dated, public-data-only evidence table with source IDs,
  filtering rules, limitations and uncertainties, and a scoped explanation
  of why accepted curation does **not** establish causality, prevalence,
  treatment efficacy or clinical utility.

**Technical corrections**, if demonstrated with tests, belong in an
engineering PR; **scientific observations** belong in a dated non-sensitive
research result, audit report or review PR. Neither replaces independent
scientific expertise. Public preparation and drift reports are described
in [the research preflight](egfr-luad-audit-preflight.md).

## Research B: GDC/cBioPortal units, cohorts and denominators

**Status: open; not independently assessed.** Historical coordination:
[archived #35](https://github.com/dovvnloading/open-cancer/issues/35).

Research question: Are historical EGFR GDC and cBioPortal observations
defined with the right units, cohort and assay coverage, and what is
needed before any cross-cohort comparison?

- [ ] Trace the exact GDC `/ssm_occurrences` project filter
  `TCGA-LUAD`, EGFR Ensembl identifier `ENSG00000146648`,
  pagination/completeness, mutation-consequence assumptions and in-memory
  case-ID deduplication. Distinguish **93 occurrences** from **78 distinct
  matching case IDs** without publishing individual identifiers.
- [ ] Trace cBioPortal study `luad_tcga_pan_can_atlas_2018`,
  mutation profile `luad_tcga_pan_can_atlas_2018_mutations`,
  `_all` sample list and Entrez gene ID `1956`; distinguish
  **86 mutation-call API rows** from unique samples/patients.
- [ ] Record exact source URL, retrieval date, source version/hash,
  field filters, missing assay-specific profiled denominators and access
  restrictions. If a live rerun is unavailable, record the technical or
  access blocker instead of manufacturing a value.
- [ ] Evaluate potential shared TCGA patients, different source pipelines,
  mutation reporting, assay panels, missingness and cohort crosswalks.
  **Never** assume GDC case IDs correspond to cBioPortal samples or use
  their counts as comparable frequencies.
- [ ] Explain the historical `coverage.gdc.note` typo describing an
  aggregate endpoint; the original immutable `manifest.inputs[]` and
  implementation show `/ssm_occurrences`. Do not change the old snapshot.
- [ ] Publish an evidence-linked, rights-safe dated research QA report
  with unresolved denominators recorded as unknown/not assessed (not zero),
  negative observations and uncertainty.

An [engineering provenance enhancement](https://github.com/dovvnloading/open-cancer/pull/46)
has already been merged, with tests and metadata-only audit tooling.
It is **partial data-engineering progress**, not completion of this
scientific assessment or a prevalence calculation.

## Research C: independent human scientific assessment

**Status: open; no qualified independent expert assessment verified.**
Historical coordination:
[archived #36](https://github.com/dovvnloading/open-cancer/issues/36).

- [ ] An appropriately qualified **independent human domain specialist**
  submits a real `synthetic: false` review as a **pull request** to
  `reviews/egfr-luad-evidence-coverage/pilot-20261110/<review-id>.json`,
  conforming to [REVIEW_PROTOCOL.md](../../REVIEW_PROTOCOL.md) and
  [the review schema](../../schemas/review.schema.json).
- [ ] The review pins the **exact original manifest SHA-256 above**.
  Review scope, sources, findings, limitations, conflicts, independence
  and expertise must be explicit.
- [ ] Maintainership verifies claimed credentials, affiliation and
  independence separately and appropriately; an AI-generated document,
  self-declared JSON and a merged PR do not validate credentials.
- [ ] The specialist assesses provenance, CIViC disease/variant
  semantics, GDC/cBioPortal observation units/cohorts, bias,
  reproducibility/drift, uncertainty and conflicting or missing evidence.
  Anything unchecked is marked `not-assessed`.
- [ ] Publish an evidence-bound opinion or decision, including a
  critical, disputed or `needs-revision` outcome if warranted. An
  unfavorable assessment is an acceptable completed **review**, not an
  endorsement. The original result and existing reviews remain unchanged.

## Completion and routing

- [ ] A reproducible, evidence-linked CIViC audit is published with
  limitations and independent interpretation boundaries.
- [ ] A reproducible, source-linked GDC/cBioPortal unit/denominator
  assessment is published (or the exact technical/access blockers are
  documented).
- [ ] A genuinely independent human scientific assessment is submitted,
  version-bound, and verified by maintainers.
- [ ] A nonclinical synthesis notes covered categories, exclusions,
  missingness, conflicts and invalid cross-source comparisons. Clinical
  validity is **not assessed**; no therapeutic recommendation follows.

**Where to put work:** code defects/implementations in engineering Issues
and PRs; scientific analyses, result manifests and review attestations in
the repository's `results/` and `reviews/` directories through PRs;
general questions, invitations and discussion in community documentation
or Discussions once available; machine status in Actions artifacts.

**Privacy:** do not include patient-level, access-restricted or potentially
identifying records, raw genomes, private credentials or unlicensed
third-party material. Respect [DATA_POLICY.md](../../DATA_POLICY.md).

The old linked issue threads are retained as **read-only historical
discussion/provenance**, not live work trackers. Their closure is an
organizational migration only and is **not evidence the research succeeded**.
