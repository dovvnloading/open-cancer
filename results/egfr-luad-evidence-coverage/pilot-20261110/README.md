# EGFR / lung adenocarcinoma public-data coverage pilot

**Status: unreviewed; descriptive data engineering, not a clinical claim.**

This run attempts five public endpoints covering three sources; an optional GDC somatic-mutation occurrence response is treated as unavailable if the API fails:
- CIViC nightly clinical evidence export (CC0 data). Restrict to accepted evidence, a molecular-profile EGFR token, and DOID:3910 (exact LUAD). Preserve evidence IDs and source links; do not copy claims.
- cBioPortal TCGA PanCancer Atlas LUAD study + mutation calls for NCBI Entrez 1956 (EGFR); aggregate by protein change entirely in memory. **No sample/patient records are saved.**
- NCI GDC TCGA-LUAD project metadata **and a bounded, public SSM occurrence query** for EGFR (ENSG00000146648) and TCGA-LUAD. We count mutation occurrence rows and deduplicate case UUIDs strictly in memory. No case UUIDs or raw mutation records are saved.

Run (Python 3.12+ standard library) from repository root:

```bash
python results/egfr-luad-evidence-coverage/pilot-20261110/analyze.py \
  --out-dir results/egfr-luad-evidence-coverage/pilot-20261110 \
  --run-id pilot-20261110
python scripts/validate.py
```

The generated `manifest.json` records sources, exact query URLs, access date, raw response SHA-256, execution environment, and artifact hashes. The three derived files are `coverage.json`, `civic_evidence_ids.csv`, and `cbio_protein_changes.csv`.

**Interpretation:** CIViC evidence items, cBioPortal mutation-call rows, and GDC SSM occurrence rows/distinct case IDs are *incompatible units*. Do not compare counts as frequencies, clinical efficacy, biological absence, independent patients or cohort-matched observations. CIViC profiles can encode multiple variants. cBioPortal's study does not represent all LUAD patients. Data changes with upstream releases; checksums aid change detection but a mutable nightly file is not a durable archive. Missing or failed mandatory upstream responses cause the run to fail. A GDC mutation-occurrence API failure, warning, or truncated/oversized response is explicitly recorded as `status: unavailable` with null counts (never zero), and its absent source is omitted from input hashes.

**GDC issue #8 addressed:** the old `/analysis/top_cases_counts_by_genes` endpoint intermittently returned HTTP 500 and returned counts with unclear agreement to mutation occurrences. It has been retired from this pilot. Its documented description is cases with a mutation per gene/project, but without independent denominator and mutation-level validation it is not treated as prevalence. The replacement `/ssm_occurrences` query filters by `case.project.project_id` and `ssm.consequence.transcript.gene.gene_id`, requests at most 1,000 records and fails closed on pagination mismatches. It counts **unique case identifiers only among the matching SSM occurrences** and never writes them. Since different projects, variant consequences and assay cohorts remain unharmonized, GDC and cBioPortal numbers cannot be meaningfully compared as frequencies.

**Remaining work:** independent replication with archived releases, sample/assay denominator definitions, nomenclature mapping, and qualified domain-expert scientific review. See issue #3.

Official links:
- [CIViC releases and licensing](https://civicdb.org/releases/licensing)
- [cBioPortal REST API](https://docs.cbioportal.org/web-api-and-clients/)
- [GDC search / retrieval API](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/)
- [Disease Ontology DOID:3910](https://disease-ontology.org/term/DOID:3910/)
