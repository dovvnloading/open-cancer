# Open Cancer — Operations and Research Integrity Report

**Reporting date:** 2026-11-11 (UTC)  
**Source revision reviewed:** [`402f666863f8c9d9c1a3c1530f8136a71df025ae`](https://github.com/dovvnloading/open-cancer/commit/402f666863f8c9d9c1a3c1530f8136a71df025ae)  
**Report type:** Evidence-linked engineering and publication-readiness assessment; **not** scientific peer review or clinical validation.

> This report is a **point-in-time snapshot** of verified GitHub data. Refer to the live issue and workflow links for changes after the reporting date. A green code/build workflow does not establish that production hosting is current.

## Executive assessment

**Engineering validation passes; public release readiness is blocked.** The most recent checked `main` validation and Explorer build succeeded. The actual externally hosted Explorer continues to serve the older revision `05e3939cab7d791e61c3efe5ef6e89b0a3f651b8`, instead of `402f666863f8c9d9c1a3c1530f8136a71df025ae`, and its public index lacks the `experimental_results` collection. The GitHub Pages deployment job was **skipped**, not passed. [Incident #29](https://github.com/dovvnloading/open-cancer/issues/29) documents the independent live API probe.

**Research publication is intentionally restricted.** New schema-valid result manifests belong to a separate experimental collection unless a SHA-256-bound human editorial publication decision is recorded. However, GitHub reports `main` is **unprotected**; CODEOWNERS is not yet a mandatory approval gate. The existing pilot's `legacy-featured` status is a discoverability decision, **not** a scientific endorsement.

**Catalog health is working.** The successful live run on November 11 checked 33 GDC TCGA projects, 32 cBioPortal PanCancer Atlas studies, and 3 CIViC nightly exports, detecting no monitored accession drift or current endpoint failures. These **68 upstream listing records** are separate from the catalog's additional custom EGFR/GDC query record (69 records total). The source-specific health ledger is kept in [issue #27](https://github.com/dovvnloading/open-cancer/issues/27), which closed after a healthy run.

## Release and QA evidence

| Control | Observed status | Evidence and interpretation |
| --- | --- | --- |
| Repository validation | **PASS** | [CI run 42956840792](https://github.com/dovvnloading/open-cancer/actions/runs/42956840792) at the checked main SHA |
| Explorer static build and local smoke | **PASS** | [Run 42956840818](https://github.com/dovvnloading/open-cancer/actions/runs/42956840818); this is an artifact build, not proof of publication |
| GitHub Pages production deployment | **NOT ACTIVE** | The `deploy` and `smoke-deployed` jobs were skipped in the checked build |
| Public Explorer revision and API completeness | **FAIL** | [Run 42956840828](https://github.com/dovvnloading/open-cancer/actions/runs/42956840828): old deployed SHA and missing experimental collection |
| Public release incident tracking | **ACTIVE — OPEN INCIDENT** | [Issue #29](https://github.com/dovvnloading/open-cancer/issues/29) automatically reopened on stale/missing API data |
| Upstream catalog monitoring | **PASS** | [Run 42955631517](https://github.com/dovvnloading/open-cancer/actions/runs/42955631517), [ledger #27](https://github.com/dovvnloading/open-cancer/issues/27) |
| Durable versioned reviewer-kit release | **NOT PUBLISHED** | No GitHub Releases were listed; reviewer-kit build alone is a time-limited Actions artifact |
| Branch protection and mandatory publication review | **NOT ENFORCED** | The public branch metadata reported `protected: false`; see [governance runbook](../publication-governance.md) |
| Independent cancer-genomics interpretation | **PENDING** | [External reviewer invitation #11](https://github.com/dovvnloading/open-cancer/issues/11); existing version-bound methods review is project-internal/agent-led, not verified external scientific adjudication |

## Scientific pilot: what was and was not established

The original EGFR/lung adenocarcinoma pilot contains source-specific, historical **descriptive observations**, not a clinical finding: **93 GDC simple somatic mutation occurrence rows**, **78 distinct GDC case IDs among those rows**, **86 cBioPortal EGFR mutation-call rows**, and **56 CIViC accepted evidence items selected under a strict EGFR/LUAD filter**. These are different entities across differently selected sources; they must not be pooled, used for cohort prevalence, interpreted as treatment efficacy, or presented as independently scientifically validated. See [original source-coverage snapshot](../../results/egfr-luad-evidence-coverage/pilot-20261110/coverage.json) and [review protocol](../../REVIEW_PROTOCOL.md).

The [reviewer kit](../../review_packets/egfr-luad-v1/README.md) records original manifest SHA-256, source-query identifiers, mutable-source caveats, and deterministic non-sensitive derived artifacts. The original source-response bytes are not archived here: exact upstream input replay may be impossible once source versions change. **A checksum match validates bytes, not scientific correctness or legal reuse permission.**

## Remaining work and acceptance evidence

| Priority | Owner / issue | Next action | Evidence required to close |
| --- | --- | --- | --- |
| **P1** | Hosting owner, [#14](https://github.com/dovvnloading/open-cancer/issues/14) and [#29](https://github.com/dovvnloading/open-cancer/issues/29) | Enable Pages (GitHub Actions source and opt-in variable) **or** republish via the existing host's supported deployment mechanism. Set the public monitor URL to the production host. | Public index reports designated revision; `experimental_results` exists; deployed smoke passes; incident closes and rollback ownership is recorded |
| **P2** | Repository owner plus separate trusted reviewer, [#18](https://github.com/dovvnloading/open-cancer/issues/18) | Enable active branch ruleset/required CI, CODEOWNERS review, author-independent approvals, and bypass controls | An unapproved publication-registry PR is blocked; approved reviewed PR can merge only with green CI |
| **P2** | Release maintainer, [#17](https://github.com/dovvnloading/open-cancer/issues/17) | Inspect non-sensitive reviewer bundle; publish a protected `reviewkit-egfr-luad-v1.*` tag; verify the GitHub Release download and digest; set persistent checksum-monitor variables | Tagged release URL, source commit and original manifest SHA recorded; independent SHA pinned and scheduled link verification passing |
| **Research** | Qualified external cancer-genomics reviewer, [#11](https://github.com/dovvnloading/open-cancer/issues/11) and [#3](https://github.com/dovvnloading/open-cancer/issues/3) | Independently assess original evidence, ontology mapping, cohorts, missingness, interpretation and negative evidence | Credible, version-bound review with verifiable identity, expertise, conflicts and substantive checks; critical reviews remain visible |

## Reporting and change-control conventions

- Report observed CI results, live deployment status, published release artifacts, and scientific review **separately**. A successful build is not a successful deployment; a reproduced calculation is not independent interpretation.
- Include exact Git revision and verified timestamp in future reports. Treat issue comments and status ledgers as evolving snapshots, not a replacement for an append-only research provenance record.
- Maintain a negative and disputed-evidence trail. No clinical endorsement, patient-level redistribution, or self-attested scientific credential verification is implied here.
- Update this report through a normal reviewed pull request; do not silently amend a historical report after publication. Future reports should use new dated files or an automated issue ledger rather than rewriting the original snapshot.

### Evidence links

[GitHub main branch](https://github.com/dovvnloading/open-cancer/tree/main) · [Live Explorer incident #29](https://github.com/dovvnloading/open-cancer/issues/29) · [Catalog health ledger #27](https://github.com/dovvnloading/open-cancer/issues/27) · [Publication governance](../publication-governance.md) · [Explorer release instructions](../explorer-releases.md) · [Reviewer release instructions](../reviewer-release.md)
