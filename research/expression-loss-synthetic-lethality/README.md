# Expression-loss synthetic lethality across 1,103 cancer cell lines

**Nonclinical computational research. Unreviewed. No treatment advice.** Record-level
public cell-line data only (DepMap 24Q4, DEMETER2, PRISM 24Q2: CC BY 4.0; Sanger Project
Score: research-use licence, summary statistics only). No patient-level or controlled data.
Literature novelty of the candidates below has **not** yet been assessed; nothing here is a
claim of priority.

This is the programme's first **record-level, genome-wide discovery screen**. Earlier notes
confronted theory with published aggregates. Here the unit of observation is the cell line,
discovery and replication use independent perturbation technologies and laboratories,
and the method is calibrated against known biology before any new candidate is believed.

## Question

When gene **B** is transcriptionally off in a cancer cell line, does the line become
dependent on a different gene **A**? Such B-off → A-dependence pairs are candidate
synthetic-lethal targets for tumours that silence or delete B.

## Results at a glance

| Result | Evidence | Status |
| --- | --- | --- |
| **Method calibration.** A lineage-, screen-strength- and 20-expression-PC-adjusted binary "B off" scan (28.4M tests, Bonferroni) returns 296 pairs. Its top is dominated by **known** loss→dependency biology: FAM50A←FAM50B (r=0.74), RPP25L←RPP25, DDX3X/EIF1AX/ZFX←chrY loss, NAA10←NAA11, DNAJC19←DNAJC15, TTC7A←TTC7B, HBS1L/PELO←FOCAD. | 4/4 screenable benchmark pairs recovered; 38/42 paralog hits and 18/18 chrY hits replicate in independent Sanger CRISPR screens | Method works |
| **v1 failed and was replaced.** A continuous-expression scan without expression PCs gave 29,356 hits dominated by mesenchymal/interferon "hub" genes and missed the benchmark. | `derived/v1_*` | Retained negative result |
| **H1: GJA1 (connexin-43) expression masks nucleotide-synthesis dependencies in pooled CRISPR screens.** GJA1-off lines are more dependent on 15 held-out de novo purine/pyrimidine enzymes (mean r=0.106, p<1e-4; 14/15 positive on adherent lines). The effect is (i) in adherent lines (0.144), not suspension lines (0.021); (ii) absent from four essential **protein-machine** control sets (ribosome, proteasome, spliceosome, Pol II, all ≈0); (iii) replicated in RNAi (0.064) and in independent Sanger CRISPR screens (0.135, p<1e-4); (iv) **attenuated ~5× for drugs inhibiting the same enzymes on the same 625 cell lines** (paired line bootstrap: CRISPR−drug difference 95% CI [0.064, 0.174]; drug effect CI [−0.006, 0.060]); (v) genome-wide, GJA1-protected essential genes are enriched for metabolic enzymes (4.7×, p=1.1e-4), headed by nucleotide, hexosamine, pentose-phosphate, one-carbon, CoA and FMN synthesis. | `HYPOTHESES.md` H1 + Addenda A–C | **Supported**: best explained by gap-junction metabolic cross-feeding from neighbouring cells in pooled screens, i.e. a systematic bias in pooled-screen dependency estimates for small-metabolite pathways |
| **H2: ALDH2 loss → Fanconi anaemia pathway dependency.** 16 held-out FA/ICL-repair genes shift in DepMap CRISPR (mean r=0.085, p<1e-4, 15/16 positive) and in independent Sanger CRISPR (0.111, p<1e-4, only 30 ALDH2-off models); MCM2-7 negative control ≈0 in both. Not seen in RNAi (−0.001). | `HYPOTHESES.md` H2 | **Supported in two independent CRISPR screens; RNAi non-replication recorded.** Mechanistically matches acetaldehyde-induced interstrand crosslinks (shown in mice; Langevin et al. 2011) |
| **H3: ACOT9 loss → succinyl-CoA metabolism.** Seed pairs do not generalise to the set (−0.014, p=0.81). | | **Falsified** |
| **Strict candidate queue.** 51 non-paralog pairs are lineage-robust and replicate in RNAi or in Sanger (≥5 B-off models); 33 are explained by known proxies (NAA11 co-activation with cancer-testis antigens, MSI markers → WRN, 9p21 co-deletion with FOCAD, H1). **18 remain unexplained**, led by ITPA←PUDP (r=0.40; Sanger 0.42), SRD5A3←HTATIP2, HSD17B12←LDAH, WRAP73/SPICE1←HYLS1, STAMBP←VPS35L. | `derived/strict_candidate_queue.csv` | Leads for literature triage and orthogonal testing |

## What this changes

1. **H1 is a measurable bias in the most widely used cancer-target resource.** Pooled
   screens underestimate how much adherent, connexin-43-expressing cells depend on enzymes
   that make small diffusible metabolites, because knockout cells are fed by unedited
   neighbours. Drug treatment removes the donors, and the protection mostly disappears. Target
   prioritisation from pooled CRISPR data (e.g. DHODH, IMPDH2, GART, PPAT) should account for
   GJA1/connexin expression, or use arrayed or drug data. Falsifiable next test: an arrayed
   (non-pooled) knockout or a GJA1-knockout co-culture should remove the GJA1 association.
2. **H2 gives a testable vulnerability hypothesis.** ALDH2-low tumours may carry a higher
   endogenous crosslink load and depend on the FA pathway. Next tests: (a) ALDH2-off lines and
   crosslinker sensitivity in PRISM and GDSC; (b) the germline ALDH2*2 (rs671) genotype of cell
   lines as a natural experiment independent of expression; (c) literature priority check.
3. **The candidate queue is the next discovery front.** It consists of record-level,
   replicated and lineage-robust pairs with known artefact classes removed.

## Files

- `pipeline.py` – discovery screen (numpy/pandas; data not in repo). `annotate.py` – Ensembl
  paralogy (cached in `derived/ensembl_paralogues.json`), chrY and penetrance classes.
  `pathways.py`, `followup.py`, `bootstrap_h1.py`, `enzyme_enrichment.py` – frozen hypothesis
  tests in the order declared in `HYPOTHESES.md`. `queue.py` – strict candidate queue.
- `derived/` – all summary outputs; `depmap24q4_extract_selected_pairs.csv` is a small CC BY 4.0
  DepMap extract (4 gene-effect and 4 expression columns) used by the offline tests.
- `REPORT.md` – methods, every number, limitations and licences.
- Tests: `tests/test_research_expression_loss_synthetic_lethality.py` (stdlib, offline).
