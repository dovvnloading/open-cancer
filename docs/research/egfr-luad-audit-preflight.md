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

**Do not mark these scientific research deliverables complete solely because the machine-generated aids exist. Their archived Issues were closed to clean up tracker routing, not because scientific work was completed.**


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

Any changed source metadata needs a domain expert to inspect the original CIViC record, cited source and disease/variant context before interpretation. This workflow does not satisfy the human-reviewed research acceptance criteria on the [work page](./egfr-luad-open-work.md).


## Follow-up: namespace-equivalent DOIDs are not upstream changes

A first live run on 2026-11-11 reported all 56 historical selected CIViC IDs as changed solely because the original derived CSV recorded `DOID:3910` while the upstream TSV uses `3910`. These strings encode the **same Disease Ontology term**, not a changed cancer type. The audit now normalizes the optional `DOID:` prefix **only for field-change comparison**, retaining the raw original and current values for provenance. Differing DOID numbers and lost strict-LUAD inclusion remain observable. The original source content hash changing still indicates different source **bytes**, not necessarily any changed selected evidence item. An unchanged selected-label comparison cannot establish upstream statement/primary-study correctness or source-wide identity.

## Historical cBioPortal/GDC query provenance and metadata drift

The live audit's `public-source-metadata-summary.json` now includes
`historical_query_provenance`. It is generated **offline** from the original
November 10, 2026 manifest and coverage file, checks that four manifest input
URLs and response fingerprints match the original source selection, and
**fails closed** if the recorded mutation URL, GDC filters, or requested page
shape differ. It does not contact any mutation, case or sample endpoint.

| Original request | Selection and observation unit | What is still unknown |
| --- | --- | --- |
| cBioPortal `luad_tcga_pan_can_atlas_2018_mutations/mutations` | `sampleListId=luad_tcga_pan_can_atlas_2018_all`, `entrezGeneId=1956`, **86 mutation-call rows** | Unique affected samples/patients, EGFR-assayed denominator, gene-panel missingness, mutation-profile coverage |
| GDC `/ssm_occurrences` | `case.project.project_id=TCGA-LUAD` AND `ssm.consequence.transcript.gene.gene_id=ENSG00000146648`; request `case.case_id,ssm.ssm_id`, `from=0`, `size=1000`; **93 occurrence rows**, **78 distinct case IDs in those rows** | Full assayed cohort denominator, variant-consequence inclusion, project/pipeline overlap with cBioPortal, and any new live occurrence counts |

The cBioPortal `_all` list covers the study's samples, whereas the separately
defined `_sequenced` list describes samples with mutation data; even that
list alone does **not** prove every sample was assayed for EGFR across gene
panels. See [cBioPortal case-list and gene-panel formats](https://docs.cbioportal.org/file-formats/)
and the [cBioPortal API guide](https://docs.cbioportal.org/web-api-and-clients/).
The [GDC data-analysis reference](https://docs.gdc.cancer.gov/API/Users_Guide/Data_Analysis/)
defines an SSM *occurrence* as a mutation applied to a case, not the number
of assayed cases. The **historical** code checked that the one requested GDC
page's `pagination.from=0`, `len(hits)=pagination.total`, and
`total<=1000`; the metadata-only audit does **not** repeat that row-level
check and explicitly reports `completeness_reverified_live=false`.
There was **no consequence-type filter** in the original SSM request.

The summary separately provides `metadata_response_byte_identity_with_original`
for CIViC TSV, cBioPortal study metadata, and GDC project metadata. These
booleans mean *only* that the current raw response hash matched the stored
historical SHA-256 at this retrieval; they do **not** establish immutable
upstream releases, scientific agreement, or unchanged unqueried records.
The original November 10 cBioPortal mutation-call and GDC occurrence response
bytes are not archived, and this live metadata audit never pretends to
replay them. Both EGFR-assayed denominators remain `null`, and cross-source
prevalence remains uncomputable. The frozen original `coverage.gdc.note`
mentions an "aggregate endpoint" in error; the version-bound manifest
inputs and analysis code identify the actual `/ssm_occurrences` request.
Do not edit the original checksum-bound pilot to correct that historical
wording.

These assertions have synthetic offline regression tests; an independent
study/profile, sample-list, consequence, coverage and cohort-overlap assessment
by a qualified human specialist remains outstanding under
[issue #35](https://github.com/dovvnloading/open-cancer/issues/35).
