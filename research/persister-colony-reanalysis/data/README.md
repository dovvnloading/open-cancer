# Source and licence record

**Nonclinical cell-line experimental aggregates and published model outputs.
Unreviewed research. No patient data. No treatment advice.**

## What is available

The colony file contains literal numerical annotations from the article's
Fig. 3e. The exact printed integers were inspected independently by two AI
readers; this is an internal transcription check, not external review. The
figure's early/late reporting intervals were converted from weeks to days.
The late-positive count is additional to the early-positive count. It is
conditioned on the wells without an early colony, as defined in the
Supplementary Note's fluctuation-estimator derivation.

The population raw data were **not retrieved**. There is no populations.csv.
The selected rows of published_parameters.csv come from Supplementary Tables
1, 2, 4 and 5. They document untreated growth parameters, prior support, the
authors' fitted effective persister loss, and their inferred chronological
phenotypic conversion rates. They are not raw trajectories or independent
colony outcomes. Every row includes its source URL, original file, exact
sheet/cell or figure locator, access date, and licence.

sources.json distinguishes the **advertised, locally unverified** Mendeley
archive digest from hashes of actual downloaded supplementary-file bytes.
manifest.json pins newline-normalized local text contents. A matching hash
establishes file identity, not the truth of an experiment or its interpretation.

## Reuse permission and attribution

The [Mendeley V2 dataset](https://data.mendeley.com/datasets/mvfm7hs9kw/2) and
its file-information panel explicitly identify **CC BY 4.0**. The
[article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9279152/) also explicitly
licenses its content under CC BY 4.0, including the selected supplementary
material in the absence of a contrary credit line. The
[licence](https://creativecommons.org/licenses/by/4.0/) permits redistribution
and adaptation with attribution, a licence link, and an indication of changes.
Both sources were accessed on 2026-11-18. Public accessibility alone was not
treated as permission.

Credit: Mariangela Russo; Simone Pompei; Alberto Sogari; Mattia Corigliano;
Giovanni Crisafulli; Alberto Puliafito; Simona Lamba; Jessica Erriquez; Andrea
Bertotti; Marco Gherardi; Federica Di Nicolantonio; Alberto Bardelli; Marco
Cosentino Lagomarsino. **A modified fluctuation-test framework characterizes
the population dynamics and mutation rate of colorectal cancer persister
cells.** Nature Genetics (2022),
[doi:10.1038/s41588-022-01105-z](https://doi.org/10.1038/s41588-022-01105-z).
Dataset: **Data and code from Russo et al.**, Mendeley Data, V2,
[doi:10.17632/mvfm7hs9kw.2](https://doi.org/10.17632/mvfm7hs9kw.2).

Changes: selected figure annotations and supplementary cells were reorganized
into small CSV tables; weeks were converted to days; the new calculations and
interpretations are this project's work. Original data retain CC BY 4.0 and
are **not relicensed under the repository's MIT licence**. No author
endorsement is implied. No source images, large archives, raw genomes, or
patient records are redistributed.

## Exact source retrieval and independent re-extraction

The article XML and supplementary bundle were obtained from the public
Europe PMC full-text services linked in sources.json. The bundle contains the
publisher's original spreadsheet, supplementary PDF, and figure image.
fetch_sources.py verifies the exact selected member byte hashes. The
API-generated outer ZIP is not pinned, because packaging can change while
the actual member bytes remain identical.

From the repository root, using a scratch directory outside the repository:

```sh
python research/persister-colony-reanalysis/fetch_sources.py --output-dir ../russo-sources --paper-supplement
python research/persister-colony-reanalysis/extract_supplement.py --workbook ../russo-sources/41588_2022_1105_MOESM4_ESM.xlsx
```

The extractor reads the original XLSX XML using only the standard library,
preserves superscript formatting, and checks every selected cell against the
committed CSV. It does not execute the original authors' code or infer values
from figure bar heights. Figure annotations still require visual inspection
of the hash-pinned source image.

For the uncompleted raw-data step:

```sh
python research/persister-colony-reanalysis/fetch_sources.py --output-dir ../russo-raw --download-archive
```

Alternatively provide an official archive obtained manually:

```sh
python research/persister-colony-reanalysis/fetch_sources.py --output-dir ../russo-raw --archive /path/to/Russo_et_al_v2.zip
```

Only a complete checksum match permits source inspection. The helper lists
candidate members; it cannot claim that unknown files are the requested
assays. A future successful raw download still requires mapping file/sheet
locators, checking clone/replicate identities and missingness, cross-checking
the published counts, extracting population measurements, and rerunning the
frozen analysis with any deviations explicitly recorded.

## Source qualifications

The supplementary workbook's interval headings and its caption disagree:
the headings say HDI 3%/97%, while the caption says 2.5%/97.5%. The selected
endpoint cells are preserved with the heading's label and the discrepancy
recorded. They are described only as **published interval endpoints**.
The values in Table 5 are previously inferred phenotypic resistance rates,
not new measurements and not genome-wide mutation rates. The second
subcolumns lack explicit labels; their uncertainty interpretation is linked
to the supplementary derivation and is not used in the new likelihood.

The source's approximate starting-cell population descriptions do not supply
a measured rare-resistant-founder count. No founder absence or exact
population-capacity constraint is inferred from them. Missing raw evidence
is kept missing.
