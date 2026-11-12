> **Historical invitation only.** Open Cancer now operates as a compute-first AI software project; no recruitment or external reviewer participation is required or assumed. This page is retained for archival context, not as an active project milestone.

# Invitation: independent EGFR/lung adenocarcinoma reviewers

Open Cancer welcomes independent, critical scientific evaluation of its public-data **EGFR/lung adenocarcinoma evidence-coverage pilot**. **This is an invitation, not an unresolved engineering issue, a claim that anyone has accepted, or an endorsement.**

The project has performed automated structural validation and an internal methods rerun, but the original pilot **has not been independently assessed for scientific interpretation or clinical validity**. The goal is to find mistakes, selection bias, limitations, false positives/negatives and missing evidence. Critical, disputed or negative conclusions are as valuable as favorable ones.

## What reviewers can examine

- [Original research task](../../tasks/egfr-luad-evidence-coverage.json), [manifest and non-sensitive outputs](../../results/egfr-luad-evidence-coverage/pilot-20261110/), and [research workstream #3](https://github.com/dovvnloading/open-cancer/issues/3).
- [Reproducible reviewer kit](../../review_packets/egfr-luad-v1/README.md), [scientific checklist](../../review_packets/egfr-luad-v1/SCIENTIFIC_CHECKLIST.md), [review protocol](../../REVIEW_PROTOCOL.md) and [review JSON schema](../../schemas/review.schema.json).
- Original result manifest SHA-256: `9fce9beaaf873b5f0b0a14e338c3b372d8aa5d19c03be11b911ea1eb2b92931e`. A new upstream retrieval may differ from historical inputs because CIViC, GDC and cBioPortal sources are mutable.
- The reviewer archive can be rebuilt deterministically from the repository. The GitHub Actions packet artifact is temporary; **a durable GitHub Release has not yet been verified**. See the [release runbook](../reviewer-release.md) and [issue #17](https://github.com/dovvnloading/open-cancer/issues/17).

The original report records **56 CIViC accepted evidence items**, **86 cBioPortal mutation-call rows**, **93 GDC SSM occurrence rows** and **78 distinct GDC case identifiers among those rows**. These are different record units, from differently selected sources and without harmonized assay denominators. **They are not comparable prevalence estimates or treatment-efficacy evidence.** The original manifest asserts no scientific findings (`findings: []`).

## Focused review opportunities

| Specialty | Work item | Questions that need evidence |
| --- | --- | --- |
| Molecular oncology / cancer variant curation | [CIViC evidence and ontology audit #34](https://github.com/dovvnloading/open-cancer/issues/34) | Exact DOID:3910 filter, broader NSCLC exclusions, compound molecular profiles, acceptance status, repeated publications and contradiction/negative evidence |
| Cancer genomics / bioinformatics | [GDC and cBioPortal cohort audit #35](https://github.com/dovvnloading/open-cancer/issues/35) | SSM occurrence versus distinct case versus mutation call, selected molecular profile/sample lists, missing profiled denominators, cross-cohort comparability |
| Independent qualified human scientific reviewer | [Version-bound expert assessment #36](https://github.com/dovvnloading/open-cancer/issues/36) | Independent, evidence-specific judgment on biological interpretation and limitations, with the possibility of `needs-revision` or `disputed` outcomes |

Reviewing only one relevant section is valuable: do not mark unexamined categories as passed. A reviewer should report source URL/identifier, query or filter, source version and date, observations, uncertainties, limitations, and whether the original input bytes could be recovered. The 12-item CIViC QA sample is intentionally nonrepresentative.

## How to contribute

1. Start with the [reviewer kit](../../review_packets/egfr-luad-v1/README.md) and [checklist](../../review_packets/egfr-luad-v1/SCIENTIFIC_CHECKLIST.md). For a source-specific discrepancy or correction, contribute an evidence-backed audit to issue [#34](https://github.com/dovvnloading/open-cancer/issues/34) or [#35](https://github.com/dovvnloading/open-cancer/issues/35), or a pull request with a new dated, non-sensitive result. Use those issue comments for **specific findings**, not generic introductions.
2. For a formal scientific interpretation assessment, follow the [review protocol](../../REVIEW_PROTOCOL.md) and submit a PR containing a real, `synthetic:false` review under `reviews/egfr-luad-evidence-coverage/pilot-20261110/<review-id>.json`, bound to the **exact original manifest SHA** above. Link the PR to [review deliverable #36](https://github.com/dovvnloading/open-cancer/issues/36).
3. Describe scope, relevant domain expertise, independence and conflicts. Maintainers must independently verify any claimed credentials; a GitHub handle or a schema-valid JSON record is not identity verification. Do not post private identifying documents or credentials in public issues.
4. Publish flaws and negative findings even if the assessment recommends revision. An agent's internal methods rerun, a favorable review, or a code merge cannot establish clinical validity or recommend patient treatment.

Questions about contributing can be raised in a relevant pull request or an actionable issue backed by evidence. **Repository Discussions is currently disabled**, so there is no Discussions thread or contact channel to claim. When a maintainer enables it, this page can link to a dedicated community discussion without turning that invitation back into an open engineering ticket.

This invitation is ongoing; publication of the invitation is not proof of outreach, commitments, review completion, verified scientific credentials or healthcare suitability.
