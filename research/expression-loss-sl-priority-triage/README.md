# Publication priority and one DNA-loss test of the expression-loss screen

**Nonclinical, unreviewed computational research. No treatment advice.** The published literature predates both headline hypotheses in PR #94. The strongest unlabeled pair, PUDP loss and ITPA dependency, was already an expression-based prediction. Its DNA-loss test here is **inconclusive**. Six matched donors carried deep relative loss, below the frozen minimum of eight.

## What was already reported?

The unit of priority is the claim. For each pair, A ← B means that lower expression of B predicts a stronger fitness effect of knocking out A. PUBLISHED requires prior evidence for the same genes, direction and fitness endpoint, including a signed coefficient in a published prediction table. PARTIALLY KNOWN identifies a specific unsigned predictor, related mechanism, or separately measured links through a shared subtype. NOT FOUND means the search did not identify a prior report of that relationship. None of these categories establishes a causal two-gene interaction.

| Item | Priority | What matters for interpretation |
| --- | --- | --- |
| H1 · GJA1 and metabolic dependencies | PUBLISHED | Dvir et al., Fig. 4, reported weaker metabolic dependencies at high GJA1 and proposed gap-junction rescue by neighbouring cells. PR #94 added drug comparisons and control gene sets. Channel-specific rescue across the panel remains an inference. |
| H2 · ALDH2 and Fanconi/ICL repair | PUBLISHED | Yang et al. identified the dependency and rescued it by restoring ALDH2 in human AML. They also examined solid cancer lines. The original screen's RNAi non-replication remains relevant. |
| ITPA ← PUDP | PARTIALLY KNOWN | PUDP was the top expression feature in a 2023 ITPA model. A public 24Q2 model also ranked DNA copy number at MIR4767, inside the PUDP locus, first. Both feature rankings are unsigned. |
| SRD5A3 ← HTATIP2 | PUBLISHED | Signed model coefficient. HTATIP2/TIP30 can be silenced by promoter methylation. |
| SUCLG1 ← ACOT9 | PUBLISHED | Signed model coefficient. A separate study of SUCLG2 tested a different direction and subunit. The broader pathway prediction failed in #94. |
| HSD17B12 ← LDAH | PUBLISHED | Signed expression coefficient. Differentiation or progression state could contribute. |
| WRAP73 ← HYLS1 | PARTIALLY KNOWN | HYLS1 was the top unsigned predictor. WRAP73 and SPICE1 physically associate. |
| SPICE1 ← HYLS1 | PUBLISHED | Signed expression coefficient. The two genes also have connected roles at centrioles. |
| CCND2 ← RBFOX2 | PARTIALLY KNOWN | Circular transcripts from chicken RBFOX2 regulate CCND2 mRNA. Human CCND2 knockout fitness after loss of linear RBFOX2 remains unresolved. |
| MYCN ← UNC93B1 | NOT FOUND | UNC93B1 appears in a bulk neuroblastoma immune signature. Its ability to mark a cell-line subtype remains uncertain. |
| MYCN ← SLC8B1 | NOT FOUND | MYCN dependency varies with neuroblastoma state. The SLC8B1 link remains open. |
| MYCN ← FUCA1 | PARTIALLY KNOWN | FUCA1 expression falls with regional 1p36 loss in neuroblastoma, which co-occurs with MYCN amplification. This is a strong subtype-proxy candidate. |
| MYCN ← ARRDC2 | NOT FOUND | A mouse neuroblastoma expression observation lacks a fitness contrast. |
| IRF8 ← PSEN2 | NOT FOUND | IRF8 dependency occurs in KMT2A-rearranged AML. Evidence linking low PSEN2 to that subtype is missing. |
| FLVCR1 ← LITAF | NOT FOUND | BCL6 repression and lymphoma-associated LITAF methylation or deletion make subtype a plausible explanation. |
| MYCN ← SHFL | NOT FOUND | MYCN suppresses interferon programs, and SHFL is interferon-inducible. Reverse regulation could generate the association. |
| STAMBP ← VPS35L | PUBLISHED | A 2022 conference abstract reports restored viability and trafficking after VPS35L re-expression in STAMBP-depleted cells. The signed predictor also appears in a journal supplement. |
| MBNL1 ← RRAS2 | PUBLISHED | A small signed expression coefficient occurs in the expression-only model. The multi-omics coefficient is zero. |
| DMXL1 ← CPNE2 | NOT FOUND | A known DMXL1/DMXL2 paralog relationship leaves CPNE2's role unresolved. |
| KRAS ← EVL | PUBLISHED | The KDS30 study reports low EVL in strongly KRAS-dependent mutant lines. The EVL locus is also subject to methylation in colorectal cancer. |

For the claims, DOI-bearing sources, proxy evidence, and exact dated queries for every entry, see [the 20-row triage CSV](triage.csv). Principal sources include [Dvir et al., Fig. 4](https://doi.org/10.1093/genetics/iyac134), [Yang et al., Fig. 4 and Fig. S7](https://doi.org/10.1158/2159-8290.CD-20-1542), [Rosenski et al., Additional file 1, Table 3](https://doi.org/10.1186/s12920-023-01446-6), [Shi et al., Supplementary Tables 1 and 2](https://doi.org/10.1038/s43018-024-00789-y), [Bravo and Kwon, printed pp. 254–255, no abstract DOI located](https://www.sacnas.org/hubfs/AbstractBook2022.pdf), and [Tyc et al., Figs. 2 and 4](https://doi.org/10.1016/j.isci.2023.106082).

## Search coverage and priority limits

Searches on 2026-11-22 included PubMed, Europe PMC, preprints, source-restricted web searches and DepMap publications. Search phrases, gene aliases, platform and access status appear per item. Google Scholar blocked direct results or returned a verification page. Search-engine indexing of Scholar is incomplete. Direct bioRxiv and medRxiv search pages also restricted access. Europe PMC preprint records, publisher-filtered results where indicated, public DOI records and web-indexed preprint pages supplied additional coverage. No browser download was accepted. An inaccessible service contributes no evidence of absence. Supplements, conference abstracts, unindexed tables, and later reports remain potential sources of missed priority.

We matched the 18 unordered symbols in the human SynLethDB 2.0 snapshot (36,746 rows). The single exact symbol hit, EVL/KRAS, points to a perturbation of EVL in KRAS-mutant cells. Direction required checking the primary paper. We also inspected the signed Shi model coefficients, looking up the row EXP.B and column A, and unsigned rankings in Rosenski. A model-wide validation score supplies no individual-feature significance, and validation that reuses DepMap lines supplies no independent cohort. A zero coefficient or a missing row is not evidence against #94. [SynLethDB snapshot](https://github.com/JieZheng-ShanghaiTech/SynLethDB/blob/fbd20be21c73af96850e43831602acf91bc46269/synlethdb_dataset_2020/sl/Human_SL.csv) · [Database paper](https://doi.org/10.1093/database/baac030)

## One prespecified record-level test

PR #94 used RNA expression to define PUDP-off. Here the exposure is DNA copy ratio at PUDP relative to five loci elsewhere on X. The endpoint is the existing ITPA Chronos score, with lower scores indicating greater dependency. We matched on sex and lineage and used one model per donor. The DNA measurement is independent of the RNA definition. The outcome and cell lines overlap discovery. The protocol was committed to the branch before calculating results. [Frozen protocol](PROTOCOL.md) · [Protocol-only commit](https://github.com/dovvnloading/open-cancer/commit/a2cd918dd8f2b1ed0bbebec689f2ccbca588c828)

| Result | Count or estimate |
| --- | ---: |
| Screened models | 1,103 |
| Excluded for missing sex or donor / repeated eligible donor | 51 / 39 |
| Eligible donor representatives | 1,013 |
| Deep loss / retained / intermediate or gain | 6 / 873 / 134 |
| Matched loss / retained donors | 6 / 202, across 6 sex-by-lineage strata |
| Retained minus loss, weighted equally per loss donor | +0.2052 Chronos units |
| PUDP RNA off, loss vs retained | 5/6 vs 7/202 |
| STS deep loss, loss vs retained | 4/6 vs 0/202 |
| PNPLA4 deep loss, loss vs retained | 3/6 vs 0/202 |

**INCONCLUSIVE.** The descriptive contrast was positive. Each stratum contained only one loss donor, two strata had opposite-signed contrasts, and the exposure-count rule stopped inference before any permutations. No p-value was computed. Co-deletion in this panel further limits attribution to PUDP. A public 24Q2 model had already ranked copy number at MIR4767, inside the annotated PUDP interval, ahead of PUDP RNA as an ITPA predictor. The new contribution is a directional, sex- and lineage-matched deep-loss contrast with a frozen inference rule. [Published 24Q2 model](https://doi.org/10.6084/m9.figshare.26955886.v1). [Machine-readable results, diagnostics and input hashes](derived/decisive_test.json)

## Reproduction and data terms

From the repository root, run the offline replay and tests with standard-library Python.

```sh
python research/expression-loss-sl-priority-triage/analysis.py
python -m unittest discover -s tests -p test_research_expression_loss_sl_priority_triage.py -v
```

To rebuild the small extract from the original files, use the versioned Broad release linked below. Pass the unmodified downloads as the first two arguments. The command checks publisher MD5 values, streams the full CN matrix, retains eight genes, maps sex and donor grouping, and exports no source donor IDs. A disposable output path permits comparison against the committed extract.

```sh
python research/expression-loss-sl-priority-triage/analysis.py prepare OmicsCNGene.csv Model.csv reproduced.csv
```

The eight-column DNA projection, sex and eligibility annotations derive from **DepMap, Broad (2024), DepMap Public 24Q4**, v1. The release terms are CC BY 4.0. ModelID is the unit of linkage. The replay reads the existing #94 ITPA effect and PUDP expression extract at its committed precision. Full matrices, patient records and Sanger per-model measurements are excluded from this contribution. [Release DOI](https://doi.org/10.25452/figshare.plus.27993248.v1) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) · [DNA copy ratios](https://ndownloader.figshare.com/files/51065324) · [Model metadata](https://ndownloader.figshare.com/files/51065297) · [Broad's scale documentation](https://forum.depmap.org/t/values-in-omicscngene-csv-24q2-release/3390)

## Hostile self-review

An independent implementation rebuilt the donor set and weighted contrast from the original files, and checked every projected gene value across 1,103 models. It recovered six loss donors, 202 retained donors and the reported contrast. It also confirmed that three loss donors lost both neighbours, one lost STS alone, and two retained both neighbours by the same cutoffs. Even the latter two cannot isolate PUDP from unmeasured alterations. Tests cover invalid inputs, the minimum-exposure stop, both permutation tails on artificial data, portable path keys and normalized hashes. Search review distinguished direct rescue from signed prediction, reverse-direction knockdown, unsigned feature importance and related biology. The causal interpretation of #94's association remains open.

The [published-coefficient lookup](derived/shi2024_published_coefficients.csv) records each of the 18 target/feature cells in both of Shi et al.'s tables. It distinguishes an unavailable target or feature from a coefficient assigned zero. It is a small projection of their CC BY 4.0 supplement. Fig. 3A of the recursive-feature-machine preprint also labels dependencies associated with GJA1 underexpression, including GMPS, ADSL and DHFR. [Shi Table 1](https://pmc-oa-opendata.s3.amazonaws.com/PMC11358024.1/43018_2024_789_MOESM3_ESM.csv) · [Shi Table 2](https://pmc-oa-opendata.s3.amazonaws.com/PMC11358024.1/43018_2024_789_MOESM4_ESM.csv) · [Recursive feature machines, Fig. 3A](https://doi.org/10.1101/2023.12.03.569803)
