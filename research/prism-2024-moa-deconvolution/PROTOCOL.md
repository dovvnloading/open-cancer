# Protocol (frozen before any drug-gene statistic was computed)

**Nonclinical, unreviewed computational research.** Written 2026-11-22.

## Question
For compounds first screened in PRISM Repurposing 24Q2 (REP.1M and REP.300, 1,519 conditions at 2.5 uM),
which compounds kill the same cancer cell lines as CRISPR knockout of a specific gene, and is that
gene **not** their annotated target? Such pairs are mechanism-of-action (MOA) or off-target hypotheses.

## Data
PRISM Repurposing Public 24Q2 extended primary matrix and compound list (CC BY 4.0); DepMap 24Q4
CRISPR gene effect, expression, Model.csv (CC BY 4.0); DEMETER2 RNAi (CC BY 4.0); Sanger Project Score
(Sanger-screened columns; summary statistics only).

## Statistic
Per compound condition: lines with both PRISM log2 fold-change and DepMap CRISPR data. Both sides are
residualised on intercept, lineage indicators (lineages with <10 lines pooled), and 20 lineage-centred
expression PCs (CRISPR side also on screen strength). Partial Pearson r between residualised drug LFC
and residualised gene effect for every CRISPR gene with <=5% missing values. Positive r = lines killed by
the compound are the lines that depend on the gene.

* Active compound: >= 20 lines with LFC < -1.
* Hit: r above the one-sided Bonferroni threshold over genes for that compound (alpha 0.05 / n_genes).

## Calibration (must pass before any discovery claim)
Over **all** active conditions with a repurposing_target annotation (primary and new screens), the
annotated target must be the top gene, or within the top 10 genes, far more often than expected by
chance (expected top-10 rate ~ 10/n_genes per target). Pass = top-10 hit rate >= 10% of evaluable
annotated compounds **and** >= 20x the chance rate. Report the per-target-class breakdown.

## Discovery set
New-screen (REP.1M/REP.300) active conditions whose top gene is a hit and is not among annotated
targets (or the compound has no annotation). Ranked by r.

## Replication (independent of discovery data)
For each discovery pair: the same partial r using (a) Sanger Project Score fold changes for the gene
(Sanger-screened, QC-passed models mapped to DepMap IDs) and (b) DEMETER2 RNAi for the gene, on the lines
shared with the compound. One-sided, BH 5% within each replication source. A pair is **replicated** if
q <= 0.05 in Sanger or RNAi. Discovery pairs whose gene is a pan-essential "screen strength" gene are
flagged rather than discarded.

## What would count as failure
Calibration failing => no discovery claims. Discovery pairs that do not replicate are reported as
unreplicated. Literature priority is assessed separately before any novelty claim.
