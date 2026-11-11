# EGFR/LUAD source audit — reproducible preflight, not expert review

**Original frozen pilot:** [`pilot-20261110`](../../results/egfr-luad-evidence-coverage/pilot-20261110/) · **manifest SHA-256:** `9fce9beaaf873b5f0b0a14e338c3b372d8aa5d19c03be11b911ea1eb2b92931e` · **historical capture date:** November 10, 2026.

**Status:** only **offline structural checks, public evidence-ID text triage, and cross-source arithmetic** have been performed. No individual CIViC evidence page was reopened for this preflight, no new cBioPortal/GDC mutation data were downloaded, and no external qualified human scientific interpretation has been obtained. This is neither an original-source audit nor clinical validation.

## Why the next review needs several independent tracks

The original pilot includes 56 selected CIViC accepted evidence items (52 Predictive, two Prognostic, two Oncogenic), 86 cBioPortal EGFR mutation-call rows, 93 GDC simple somatic mutation occurrence rows and 78 GDC case IDs among those occurrences. These are **different types of observations** and source selections, not cohort-matched patient counts or prevalence. The original manifest asserts zero biological findings (`findings: []`). Responses are mutable public sources identified by their original SHA-256 fingerprints, but complete original raw responses are not archived.

- The official [CIViC evidence model](https://docs.civicdb.org/en/latest/model/evidence.html) describes evidence items as citable statements about molecular profiles that may include one or more variants, including fusions. Its [curation guide](https://docs.civicdb.org/en/latest/curating/evidence.html) explicitly distinguishes moderated acceptance from downstream clinical truth; many evidence items may refer to the same publication.
- The [NCI GDC API](https://gdc.cancer.gov/developers/gdc-application-programming-interface-api) exposes project and mutation metadata with specific query/filter semantics; its `/ssm_occurrences` results are not a sequenced project denominator or a complete EGFR driver-mutation prevalence estimator.
- The [cBioPortal API](https://docs.cbioportal.org/web-api-and-clients/) requires the caller to identify an exact study, molecular profile and sample list. Mutation-call rows should not be equated to distinct patients without explicitly measuring the assayed cohort and duplicates.

## Deterministic files added to the reviewer bundle

Running `python scripts/build_review_packet.py --output OUTPUT_DIRECTORY --archive REVIEW_KIT.zip` from a clean Git checkout generates two additional non-sensitive artifacts, both covered by `packet-index.json` and `SHA256SUMS`:

| File | What it does | Does not establish |
| --- | --- | --- |
| `civic_profile_triage.csv` | Inventories **all 56 selected public CIViC EIDs**, source URLs, disease labels and evidence types. Pure label heuristics flag `AND`/`OR` compound profiles, `::`/fusion partner profiles, and `EGFR Mutation` or `EGFR Rare Mutation` category-level labels | That a flagged item is wrong, that an unflagged item is correct, source status today, evidence significance or independent assessment |
| `cross_source_units.json` | Verifies the stored CIViC ID and cBioPortal aggregate totals, records GDC occurrence/case-count consistency, missing denominators, source response fingerprints and units of observation | Patient cohorts harmonized across sources, biomarker prevalence, clinical efficacy or a new upstream source version |

All CIViC queue rows are marked **`upstream_record_checked=no`** and **`scientific_interpretation_checked=no`**. This avoids turning heuristic text patterns into a false scientific-review badge. Non-flagged entries remain unreviewed. Only original public CIViC evidence IDs and derived aggregate metadata are distributed: no case identifiers, sample IDs, patient details, tokens or raw third-party API responses.

## Remaining work — explicit acceptance gates

**[#34 CIViC evidence and ontology audit](https://github.com/dovvnloading/open-cancer/issues/34):** A specialist must inspect the underlying public CIViC pages and citations, the 12-record nonrepresentative QA sample and any flagged profiles, check DOID:3910 mappings, repeated sources, evidence direction, and false positives. The frozen file **does not include identifiers for 161 excluded EGFR accepted evidence rows under broader lung labels**; a lawful new dated upstream retrieval and explicit exclusion log are required to check false negatives. Any source drift must be recorded, not overwritten.

**[#35 GDC/cBioPortal units and denominators](https://github.com/dovvnloading/open-cancer/issues/35):** A reviewer must validate exact study, gene, profile, sample list, page completeness, source-version hashes, mutation-call/occurrence semantics, consequence filters and denominators. The offline report deliberately records unknown profiled denominators as `null`. No attempt is made to reconcile GDC cases with cBioPortal samples or to estimate frequency.

**[#36 Independent expert interpretation](https://github.com/dovvnloading/open-cancer/issues/36):** A human domain specialist must submit an evidence-bound review of the original manifest using [REVIEW_PROTOCOL.md](../../REVIEW_PROTOCOL.md), with independently checked expertise, independence and conflicts. A `needs-revision` or `disputed` assessment can be an acceptable completed review *deliverable*; it is not project endorsement. No expert review has been performed by this automation.

## Reproduction constraints

The preflight is deterministic from the historical, immutable manifest and saved public derived outputs, and its fields are validated by offline tests. It intentionally **does not fetch new upstream data**, so a green test or identical ZIP does not imply that an original source URL is still live, source reuse rights have been adjudicated, or that scientific conclusions were independently reproduced. All original result files, their checksums, and linked reviews are unchanged.

**Do not close #34, #35 or #36 solely because these machine-generated aids exist.**


## Optional bounded live public metadata drift audit

The [source metadata audit workflow](../../.github/workflows/egfr-source-audit.yml) runs once after a reviewed update to its source code or via manual Actions dispatch. It is not a source of original biological claims. It **retrieves only** the openly published CIViC nightly clinical evidence summary TSV (under CIViC's stated CC0 reuse terms), the public cBioPortal study metadata response and the public GDC TCGA-LUAD project metadata response. It never calls the cBioPortal mutations API, GDC SSM occurrence endpoint, patient/case/sample endpoints or controlled-data downloads.

The [read-only audit script](../../scripts/egfr_live_source_audit.py) checks the original immutable manifest, writes a separate timestamped artifact, and captures current upstream source-response SHA-256 fingerprints:

- `civic-selected-record-drift.csv`: one row per **originally selected** CIViC EID (56), comparing current public TSV fields with the historical selected row, highlighting missing IDs, changed acceptance status, disease/DOID, molecular-profile labels or evidence types. This is **metadata drift**, not expert validation of each evidence statement or medical interpretation.
- `civic-broader-lung-candidates.csv`: current live **accepted EGFR-token records** with broader lung-disease labels omitted by the exact-LUAD rule. These are candidates to inspect, not a finding that the original filter was wrong. The original 161 excluded EIDs were not individually preserved, so today's live count cannot be retroactively treated as the original excluded set.
- `public-source-metadata-summary.json`: original-versus-live selected counts, source SHA fingerprint drift, current cBioPortal study-level sequenced sample count **if reported**, GDC project-level case count **if reported**, both with explicit warnings that neither is the EGFR-assayed denominator, and **no scientific or clinical validation**.
- `SHA256SUMS`: verification of the derived public audit artifact.

The workflow publishes the output as a **90-day Actions artifact**, not into the trusted Explorer API and not as a permanent scientific result. Fetch or schema errors fail the job rather than silently returning zero. To run locally on a clean checkout (with public outbound HTTPS):

~~~bash
python scripts/egfr_live_source_audit.py --output-dir /tmp/egfr-live-public-audit
(cd /tmp/egfr-live-public-audit && sha256sum --check SHA256SUMS)
~~~

Any changed source metadata needs a domain expert to inspect the original CIViC record, cited source and disease/variant context before interpretation. This workflow does not satisfy the human-reviewed acceptance criteria of #34, #35 or #36.


## Follow-up: namespace-equivalent DOIDs are not upstream changes

A first live run on 2026-11-11 reported all 56 historical selected CIViC IDs as changed solely because the original derived CSV recorded `DOID:3910` while the upstream TSV uses `3910`. These strings encode the **same Disease Ontology term**, not a changed cancer type. The audit now normalizes the optional `DOID:` prefix **only for field-change comparison**, retaining the raw original and current values for provenance. Differing DOID numbers and lost strict-LUAD inclusion remain observable. The original source content hash changing still indicates different source **bytes**, not necessarily any changed selected evidence item. An unchanged selected-label comparison cannot establish upstream statement/primary-study correctness or source-wide identity.
