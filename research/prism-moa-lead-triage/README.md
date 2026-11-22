# Priority review of six PRISM compound–mechanism leads

**Nonclinical, unreviewed research.** None of the six leads warrants a discovery claim. PR #96 recovered an annotated target in the top ten for 8.0% of evaluable compounds, below its frozen 10% requirement. We took its exploratory list as given. The independent drug-response test below failed its prespecified effect-size requirement. [Original screen and calibration](../prism-2024-moa-deconvolution/README.md)

## Priority as of November 22, 2026

| Lead from #96 | Priority | What the prior evidence supports |
| --- | --- | --- |
| Selpercatinib and EXT1, EXT2, SLC35B2, HS2ST1, FGFR2 | **PARTIALLY KNOWN** | Selpercatinib inhibits FGFR2, including in a cellular phosphorylation assay. FGFR-altered lines are preferentially sensitive. The four-gene heparan sulfate (HS) knockout signature was not located. HS can support both FGF and RET signaling. [Kooijman et al.](https://doi.org/10.3389/fonc.2022.953013) · [FDA review](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2020/213246Orig1s000MultidisciplineR.pdf) |
| Zanubrutinib and ERBB2/EGFR | **PUBLISHED** | A prior cell panel found 45-fold higher potency in ERBB2-amplified lines and eightfold higher potency in EGFR-mutant lines. HER2-positive cellular signaling and EGFR biochemical activity also have direct evidence. This covers the biological interpretation. The identical knockout-profile statistic was not located. [Kooijman et al.](https://doi.org/10.3389/fonc.2022.953013) · [Dostálová et al.](https://doi.org/10.1007/s10637-023-01346-7) |
| PBT-1033 (PBT2) and HSCB, LYRM4, IREB2 | **PARTIALLY KNOWN** | PBT2 perturbed Fe–S assembly proteins in *H. pylori*. Human tumor-cell experiments show that copper availability alters its cytotoxicity. The human genetic signature was not located. HSCB and LYRM4 support Fe–S assembly, whereas IREB2 regulates iron homeostasis. [Chen et al.](https://doi.org/10.1128/mbio.01229-26) · [Summers et al.](https://doi.org/10.1039/D0MT00222D) |
| Ruboxistaurin and PPM1D/MDM2 | **PARTIALLY KNOWN** | An overlapping PRISM condition was previously tested against TP53 damaging and hotspot variant status. Neither comparison passed FDR correction (0.342 and 0.980). Functional p53 can explain the paired genetic dependencies. We found no prior drug-specific PPM1D/MDM2 knockout result. [Vermeulen et al., Table S6](https://doi.org/10.1038/s41540-025-00618-7) · [Miller et al.](https://doi.org/10.1182/blood.2023020331) |
| NIM811 and SLC33A1 | **PARTIALLY KNOWN** | A public patent and the investigators’ grant report describe NIM811 and SLC33A1 loss each sensitizing EGFR-mutant lung models to osimertinib. That conditional observation leaves the pan-cancer single-agent claim open. KEAP1/NRF2 status and ER stress can affect the same cells. [Patent, figures 1–2 and 7–9](https://patents.google.com/patent/WO2024220666A2/en) · [Romero et al.](https://doi.org/10.1038/s43018-020-0071-1) |
| DARS1 with eribulin, DM4, soblidotin, mivobulin, batabulin and crolibulin | **PARTIALLY KNOWN** | A breast cancer xenograft shRNA screen identified DARS1 as a paclitaxel sensitivity gene. The six named agents’ baseline DARS1 dependency pattern was not located. Proliferation and the fixed-time assay are plausible common causes. ABCB1-mediated resistance is documented for eribulin. [Sultan et al.](https://doi.org/10.1002/1878-0261.12964) · [Eribulin resistance](https://doi.org/10.1124/molpharm.121.000254) |

Priority concerns the **biological claim**. PUBLISHED requires prior compound-specific pharmacology and a matching cellular phenotype. PARTIALLY KNOWN requires direct precedent for a component, an allied compound in the class, or the same perturbations under a different condition. NOT FOUND means the named drug-specific selectivity claim was absent from the records examined. General pathway knowledge or a change in a signaling marker alone does not establish selective loss of viability. A knockout correlation measures similar patterns of fitness across lines. It does not demonstrate binding to the gene product.

The [six-row triage table](triage.csv) records claims, evidence, citations, confounds and gaps. We also searched the released LinkD correlation table by compound names, aliases and Broad identifiers, and the Vermeulen drug–gene tables by compound and gene. The latter supplied the ruboxistaurin result. [LinkD archive](https://doi.org/10.5281/zenodo.21627794) · [Vermeulen supplementary tables](https://pmc-oa-opendata.s3.amazonaws.com/PMC12738840.1/41540_2025_618_MOESM2_ESM.zip)

The [dated search log](derived/search_log.csv) gives the executed query, source and coverage for each lead. Searches included PubMed, Europe PMC, bioRxiv, drug labels or regulatory reviews, selectivity resources and DepMap/PRISM publications. Native index queries, web-index queries and published-table lookups are identified separately. Search-engine coverage, database updates, name mismatches and unindexed supplements limit every absence claim. Replication flags in #96 apply to its top gene for each drug. They do not validate every gene in the row.

## Independent dose-response test

**FAIL.** The four-gene HS score correlated positively with selpercatinib sensitivity, but the within-lineage correlation of **0.1732** fell below the frozen minimum of **0.20**. The one-sided permutation **p = 0.0492** met the second requirement of p ≤ 0.05. Both were required. The cross-panel predictive claim loses priority at the effect size specified. Published FGFR2 pharmacology remains compatible with this result. [Frozen protocol](PROTOCOL.md) · [Machine-readable result and exclusions](derived/summary.json)

We chose the HS signature because it adds a specific four-gene claim to known kinase pharmacology, had the strongest discovery correlation of the individual leads, and could be challenged in a different drug assay. Another FGFR2 selectivity lookup would revisit an already measured property. The new assay measured 72-hour ATP viability over nine concentrations. The genetic measurements still come from DepMap 24Q4, as in #96. The independent measurement is the drug response in the Oncolines panel. Shared cell models and reused knockout measurements limit generalization. [Kooijman et al., Table S4](https://doi.org/10.3389/fonc.2022.953013.s005)

Of 134 panel models, **101** were evaluable in **19** lineages. We excluded the two known RET-altered models, four unmatched RRIDs, one ambiguous RRID, 25 models lacking lineage or complete dependency measurements, and one singleton lineage. Nineteen IC50s were below the 2,500 nM assay cap. The other 82 responses tied there. Only nine lineages varied in both response and score. Among 9,999 within-lineage permutations, 491 had a statistic at least as large as observed. A modest effect remains possible. Censoring limits resolution, ATP measures viability indirectly, and lineage stratification leaves within-lineage differences in growth, efflux and FGFR drivers. No direct inhibition of an HS enzyme follows from these data.

## Source records and reproduction

The 134-row [panel extract](derived/panel.csv) contains a cell-line name, RRID, DepMap model ID, lineage, source selpercatinib IC50 in nM including its censoring sign, four gene effects and any identifier-matching issue. The script applies the cap at replay. The two third-party sources below permit redistribution under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). This extract selects, matches and reformats records from those sources. The panel contains no PRISM or Sanger measurements. The literature table quotes only published summary statistics.

| Credited source, accessed November 22, 2026 | Original input | Local filename for optional rebuild |
| --- | --- | --- |
| Jeffrey J. Kooijman et al. (2022), Oncolines cell panel metadata | [Supplementary Table S2](https://doi.org/10.3389/fonc.2022.953013.s003) · [Figshare file 37443880](https://ndownloader.figshare.com/files/37443880) | `s2.xlsx` |
| Kooijman et al., Oncolines cell-level proliferation IC50s | [Supplementary Table S4](https://doi.org/10.3389/fonc.2022.953013.s005) · [Figshare file 37443889](https://ndownloader.figshare.com/files/37443889) | `s4.xlsx` |
| Broad Institute Cancer Dependency Map, DepMap Public 24Q4, v1 | [Model metadata](https://ndownloader.figshare.com/files/51065297) · [Release DOI](https://doi.org/10.25452/figshare.plus.27993248.v1) | `Model.csv` |
| Broad Institute Cancer Dependency Map, same release | [CRISPR gene effect](https://ndownloader.figshare.com/files/51064667) | `CRISPRGeneEffect.csv` |

From the repository root, replay the result from the committed extract and run the offline regression tests. Both commands use only the Python standard library.

```sh
python research/prism-moa-lead-triage/analysis.py --out replay
python -m unittest discover -s tests -p test_research_prism_moa_lead_triage.py -v
```

To rebuild the extract as well, place the four pinned files in one local directory and run the following command with openpyxl available. It verifies each upstream file before extracting records. UTF-8 text hashes normalize newlines to LF. XLSX hashes cover the original bytes. Comparison of the resulting panel and summary to the committed files is sufficient to reproduce the study output.

```sh
python research/prism-moa-lead-triage/analysis.py --source-dir /path/to/inputs --out rebuilt
```

The Europe PMC supplement bundle returned a connection-error document. The publisher’s Figshare files were accessible and were used for the analysis. No browser download prompt was encountered.

## Hostile self-review

An independent implementation rebuilt the join from both workbooks and both DepMap files, checked every cell, identifier and exclusion in the extract, and recomputed the ranks and permutations. It reproduced the 33 exclusions, the statistic and 491 exceedances. The protocol-only commit immediately follows the #96 branch tip and predates computation. A search of published supplements also identified the earlier TP53–ruboxistaurin comparison. Its source dataset overlaps #96, and its reference group may include TP53 variants outside the selected class. The narrow p-value crossing, many tied responses and reuse of the genetic data are reasons to retain the failed verdict and keep the biological interpretation exploratory. [Protocol-only commit](https://github.com/dovvnloading/open-cancer/commit/751125824cd49504abf40e083623a84bca4bf9f6)
