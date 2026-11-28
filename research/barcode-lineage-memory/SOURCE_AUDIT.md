# Source audit and reconstruction boundary

This contribution uses public **cell-line artificial-barcode count data**.
It includes no patient records, controlled-access material or full human
sequence files. It is a new retrieval and rerun: see
[RECONSTRUCTION.md](RECONSTRUCTION.md). The previous local files and commit
objects were unavailable. Expected content fingerprints from the recorded
conversation were tested against fresh publisher bytes, not treated as proof
that the absent files had been recovered.

## Primary source and reuse

Whiting and colleagues, [*Quantitative measurement of phenotype dynamics
during cancer drug resistance evolution using genetic barcoding*](https://doi.org/10.1038/s41467-025-59479-7),
Nature Communications **16**, 5282 (2025), published 20 June 2025.
The [publisher article](https://www.nature.com/articles/s41467-025-59479-7)
and [PMC copy](https://pmc.ncbi.nlm.nih.gov/articles/PMC12181405/) identify
**CC BY-NC-ND 4.0** terms. No separate permissive dataset license was located.
Public accessibility does not confer permission to relicense source material.

The committed contribution contains original retrieval and analytical code,
provenance metadata, and small aggregate research calculations. The original
workbook, extracted barcode tables, supplementary PDF and downloaded sequence
prefixes are not redistributed. No publisher figures are adapted. The
repository MIT license does not replace source terms. Publisher-linked
`Barcode_ABC` software is academic-only; it was not executed or incorporated.
The authors' barcode-processing implementation was not copied or executed.

## Fresh workbook retrieval

The [publisher Supplementary Data archive](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-59479-7/MediaObjects/41467_2025_59479_MOESM4_ESM.zip)
is 55,319,940 bytes. This reconstruction retrieved only the byte ranges for
`Supp_Data/Supp_Data_2.xlsx`. Its ZIP local header, member name, compression
method and sizes are validated before extraction. Extraction must match both
the ZIP member CRC32 and the independent workbook SHA256.

| Input property | Freshly matched value |
| --- | --- |
| Workbook size | 32,159,312 bytes |
| Workbook SHA256 | `5c61f93a15134b795e6655ecab18f2d30e7d52daed33168ca0a60406303ee100` |
| ZIP member CRC32 | `7d3fe2bd` |
| HCTbc barcode table | 13,962 rows × 8 integer count columns |
| SW6bc barcode table | 688,439 rows × 8 integer count columns |
| Count column order | DT1_P1, DT1_P2, DT2_P1, DT2_P2, DT3_P1, DT3_P2, DT4_P1, DT4_P2 |

A full-archive SHA256 is **not** claimed: the full archive was not downloaded.
The complete pins, original URLs and exact canonical encodings are in
[data/source_manifest.json](data/source_manifest.json).

`fetch_source.py` verifies nonnegative integral counts, unique 30-base
artificial identifiers, expected headers and row counts. It separately pins
the actual count bytes and barcode identities. Count bytes are encoded as
signed little-endian int64 in C order; barcode strings are source-ordered
ASCII separated by LF including the terminal LF. NumPy container compression,
byte order and memory layout do not change these content hashes. A cache that
keeps its source-SHA metadata but changes values or identifiers is rejected.

| Cell line | Count SHA256 | Barcode SHA256 |
| --- | --- | --- |
| HCTbc | `2213a534ce36abd298d2c587ac6f83f64ccaa2c51c4980d809245388d882fa09` | `2fa829c71b13a0fa8790f6b89bbd735e8b45b294763b9cbaeab9ab7f06a5abea` |
| SW6bc | `790dfdcb71b469afb82f219fe9a133e717dcf202477c445e9f28ac99d3c80585` | `070a2a4fcd1ab01dc1751eb45fae0b1729bc158b0f016b329ddc683cb783448a` |

All four content hashes match the earlier recorded fingerprints after fresh
retrieval and extraction. This establishes input identity at the pinned
content level; it does not establish correctness of the upstream experiment
or recover the earlier analytical code.

## Observation limits relevant to the forecast

The treated cultures share an expanded barcoded parent. P1 is measured after
treatment, at day 67 for HCTbc and day 37 for SW6bc. The authors' analysis
places the HCTbc escape transition before P1. The eight columns in each
selected processed count sheet contain treated P1/P2 observations; they do
not provide full parental POT or untreated CO distributions. This forecast
therefore adjusts for its own P1, not pretreatment abundance or untreated
drift. Other released materials were not exhaustively reprocessed here.

Counts are amplified sequencing reads rather than independent cells. The
published processing has already excluded barcodes; a complete standalone
problem-barcode exclusion list and a raw-to-filtered-count crosswalk were not
recovered in this reconstruction. Thus a dictionary selected only from the
released P1 columns cannot guarantee that the complete upstream processing
was independent of later outcomes.

The barcode-processing discussion in the source Methods explicitly discusses residual
index hopping despite unique dual indexes and exclusion of problem barcodes.
Those observations motivate a technical alternative, not a newly discovered
index-hopping mechanism. The source's diversity metrics emphasize more
abundant lineages; a forecast that depends on rare entries can be fragile
without refuting the source's distinct mechanistic analysis.

The freshly retrieved [Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-59479-7/MediaObjects/41467_2025_59479_MOESM1_ESM.pdf)
SHA256 is `c8402605f8da92bdbc03b7fc3a903c877bda0aaed56f58cb0f2e351e701604cb`.
Supplementary Note 1, section 4.3 (printed page 28, lines 463–469; PDF page
48) selects Hill order `q=2` and explains why emphasizing rare PCR/sequencing
errors can distort diversity statistics. The PDF remains untracked.

## Sequencing provenance

The study deposits raw reads under [ENA PRJEB74046](https://www.ebi.ac.uk/ena/browser/view/PRJEB74046).
The separate bounded first-record audit is described in
[data/sequencing_provenance.json](data/sequencing_provenance.json). It uses
exact sample-title matches restricted to AMPLICON records, then inspects
only the first 65,536 compressed bytes of the relevant R1 file. Only the first
header, source metadata and compressed-prefix SHA256 are retained; sequence
and quality strings are discarded.

The fresh metadata snapshot contains 115 runs. All eight HCTbc DT1–4 P1/P2
names match uniquely; every inspected first header reports instrument
`A00176`, instrument run `418`, flowcell `HW5HLDMXX`, lane `1`. The original
submitted HCTbc DT1_P1 file independently gives the same instrument fields.
The six name-matched SW6bc headers report run `419`, flowcell `HW55GDMXX`,
lane `1`. SW6bc DT1_P2 and DT3_P2 remain unmatched in this snapshot. Thus the
current reconstruction corroborates 14 ENA first headers and one submitted
first-header crosscheck; it does not carry these claims forward solely from
the earlier conversation.

This is **name-based, first-record provenance**. It does not establish the
provenance of every read in a potentially concatenated file, certify sample
identity experimentally, reconstruct the authors' filtering, or show a
raw-count match to the workbook. ENA public release dates are not sequencing
dates. Shared first-record flowcell support establishes an opportunity for
cross-library spillover; it does not estimate contamination or establish that
spillover caused the predictive signal. Unmatched sample names remain
unresolved rather than inferred from neighboring accession numbers.

## Reproduction

From the repository root:

```sh
python -m pip install -r research/barcode-lineage-memory/requirements.txt
python research/barcode-lineage-memory/fetch_source.py
python research/barcode-lineage-memory/fetch_source.py --audit-sequencing
python -m unittest discover -s tests -p test_research_barcode_source.py -v
```

The first command after installation retrieves and prepares the pinned source.
A supplied original workbook can instead be passed with `--workbook`; it must
pass the same checks. Within the checkout, source output is restricted to the
ignored `.compute` directory. The audit writes a new aggregate provenance
report there for comparison with the recorded metadata. It never downloads
complete FASTQ files or reconstructs full control counts.
