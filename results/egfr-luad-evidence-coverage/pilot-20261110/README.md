# EGFR / lung adenocarcinoma public-data coverage pilot

**Status: unreviewed; descriptive data engineering, not a clinical claim.**

This run queries four public endpoints covering three sources:
- CIViC nightly clinical evidence export (CC0 data). Restrict to accepted evidence, a molecular-profile EGFR token, and DOID:3910 (exact LUAD). Preserve evidence IDs and source links; do not copy claims.
- cBioPortal TCGA PanCancer Atlas LUAD study + mutation calls for NCBI Entrez 1956 (EGFR); aggregate by protein change entirely in memory. **No sample/patient records are saved.**
- NCI GDC TCGA-LUAD **project metadata only**, not GDC mutation evidence.

Run (Python 3.12+ standard library) from repository root:

```bash
python results/egfr-luad-evidence-coverage/pilot-20261110/analyze.py \
  --out-dir results/egfr-luad-evidence-coverage/pilot-20261110 \
  --run-id pilot-20261110
python scripts/validate.py
```

The generated `manifest.json` records sources, exact query URLs, access date, raw response SHA-256, execution environment, and artifact hashes. The three derived files are `coverage.json`, `civic_evidence_ids.csv`, and `cbio_protein_changes.csv`.

**Interpretation:** CIViC evidence items, cBioPortal mutation-call rows and GDC project metadata are *incompatible units*. Do not compare counts as frequencies, clinical efficacy, biological absence, independent patients or cohort-matched observations. CIViC profiles can encode multiple variants. cBioPortal's study does not represent all LUAD patients. Data changes with upstream releases; checksums aid change detection but a mutable nightly file is not a durable archive. Missing or failed upstream responses cause the run to fail, never to output a fabricated zero.

**Remaining work:** add mutation-level GDC analysis for the same clearly defined cohort, independent replication using archived releases, nomenclature mappings, and qualified scientific review. See issue #3.

Official links:
- [CIViC releases and licensing](https://civicdb.org/releases/licensing)
- [cBioPortal REST API](https://docs.cbioportal.org/web-api-and-clients/)
- [GDC search / retrieval API](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/)
- [Disease Ontology DOID:3910](https://disease-ontology.org/term/DOID:3910/)
