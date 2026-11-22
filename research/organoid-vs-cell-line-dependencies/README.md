# Patient-derived organoids vs 2D cell lines: CRISPR dependencies from the same laboratory

**Nonclinical, unreviewed computational research. No treatment advice.** Sanger Cell Model Passports
organoid CRISPR screens (released 2026-07-23) and Project Score cell-line screens (Sanger-screened
columns only). Used under Sanger's research-use licence: only gene-level summary statistics are
committed, never per-model values. No literature priority check has been made.

## Results

| Step (all frozen in `PROTOCOL.md` before running) | Result | Verdict |
| --- | --- | --- |
| **Primary:** colorectal organoids (84) vs colorectal cell lines (29), replicated in oesophageal organoids (59) vs lines (7) | 623 CRC hits, 481 replicated, **but** both directions are dominated by core-essential genes (spliceosome more essential in organoids, ribosome in lines). The organoid−line differences correlate genome-wide between cancer types (r = 0.75), and a multi-copy artefact gene (KRTAP4-11) is a hit | **Platform-confounded.** The difference between screening campaigns (library, guides, protocol) cannot be separated from biology. The hit list is not interpreted gene by gene |
| **Addendum A:** within-platform MSI−MSS contrast (removes between-platform library effects) | WRN positive control passes in **both** platforms (organoids p = 5e-4, lines p = 1e-5). Genome-wide transfer is underpowered: 10 MSI lines, 1 line-specific MSI gene (CEBPZ), Spearman ρ = 0.12 between platforms | Positive control works; transfer unresolved |
| **Addendum B:** is WRN-in-MSI dependency weaker in organoids specifically? | WRN organoid/line ratio 0.46 [0.24, 0.77]; comparator **KRAS** addiction 0.52 [0.34, 0.79]; **BRAF V600E** 0.16 [−0.23, 0.60]; CTNNB1 uninformative (3 mutants per platform) | **General compression, not WRN-specific** (frozen rule) |

## What was learned

1. **Selective oncogene dependencies are about half as strong in colorectal organoids as in
   same-laboratory cell lines** (WRN in MSI, KRAS in KRAS-mutant), even after normalising each screen to
   its own core-essential genes. 60% of MSI organoids vs 90% of MSI lines show a WRN effect < −0.5; 79% vs
   94% of KRAS-mutant models show a KRAS effect < −0.5. Effect sizes measured in cell lines may overstate
   how strongly patient-derived tumour models depend on these targets. Slower growth, subclonality or
   medium could all compress depletion at a fixed time point. This analysis cannot tell them apart.
2. **BRAF V600E dependency nearly disappears in organoids** (difference −0.08 vs −0.54 in lines; 2/8
   organoids vs 4/6 lines below −0.5). Exploratory interpretation: organoid medium supplies EGF, EGFR
   scores as markedly more essential in organoids in the primary comparison (subject to the platform
   caveat), and EGFR feedback is the published route by which BRAF-mutant colorectal cancer escapes
   BRAF dependence (Prahallad et al., *Nature* 2012). This is a hypothesis for a medium-swap experiment;
   it is not a finding.
3. **Method lesson.** Comparing screens across campaigns needs within-platform contrasts. Naive
   organoid-vs-line differences mostly measure the screening campaign.

## Files
`PROTOCOL.md` (primary design plus Addenda A and B, each written before its run); `compare.py`,
`within_platform.py`, `wrn_attenuation.py` (numpy/pandas/scipy; Sanger files not in the repository);
`derived/` gene-level summaries. Test: `tests/test_research_organoid_vs_cell_line_dependencies.py`.
