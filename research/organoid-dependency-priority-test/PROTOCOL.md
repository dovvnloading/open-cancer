# Frozen test of the basal EGFR dependency prediction

Nonclinical, unreviewed computational research. No treatment advice. Frozen 2026-11-22 against the open PR #98 head `2b9e7e95ce4115fad49ca6d42badb32b284533e9`. The public files and headers were checked for availability before this commit. No EGFR calls or genotype-by-EGFR summaries were examined. The freeze precedes the analysis code and output in git.

## Question and null

If the proposed EGF–EGFR bypass in claim (c) makes untreated BRAF V600E colorectal organoids more reliant on EGFR, a higher fraction should lose fitness after EGFR knockout than among BRAF-reference organoids in the same campaign. Restrict both groups to models without a reported RAS driver mutation because RAS genotype influences EGFR dependence in this system (Herranz-Ors et al., [Fig. 5](https://doi.org/10.1038/s41586-026-10830-y)). The null is equal probability of an EGFR fitness call in the two groups. This tests the basal prediction. Feedback induced only after BRAF blockade would require a different perturbation.

## Fixed public inputs and unit

Use only these dated [Sanger Cell Model Passports downloads](https://cellmodelpassports.sanger.ac.uk/downloads). Direct URLs start with `https://cog.sanger.ac.uk/cmp/download/`.

| Input | Use |
| --- | --- |
| `Organoids_fitness_scores_20260723.zip`, member `Organoids_fitness_scores_binary_matrix_20260723.tsv` | Only the EGFR binary essentiality row. 1 indicates depletion at the authors' model-specific 5% FDR. #98 used fold changes, with no binary calls or BRAF-to-EGFR genotype contrast. These are different measurements derived from the same screens. |
| `model_list_20261013.csv` | Select organoids classified as Colorectal Carcinoma. Use the anonymous donor key only to avoid pseudoreplication. |
| `model_dataset_availability_20260727.csv` | Require explicit availability of Mutation Sanger Organoid (WGS) in the model itself. |
| `mutations_summary_20260724.zip`, member `mutations_summary_20260724.csv` | Use Sanger calls, gene, protein change and cancer-driver flag. This is a selected variant summary. Absence of a listed variant cannot establish a confidently callable wild-type locus. |

Join on model ID. Select Sanger QC-passed screens. Require a nonmissing donor key, model metadata and organoid WGS coverage. Exclude models with a reported KRAS, NRAS or HRAS cancer driver. The V600E group requires an exact driver call `BRAF p.V600E`. Exclude other or ambiguous BRAF records, including models with both V600E and another reported BRAF variant. The reference group has no reported BRAF variant. If multiple eligible cultures come from one donor, keep the lexicographically first model ID before checking EGFR call availability. Require one unambiguous EGFR row and at most one eligible assay column per model. Exclude missing calls without replacement. Do not impute. Report aggregate attrition and no donor or model IDs.

## One test and decision rule

Use a single two-sided Fisher exact test conditional on the margins (sum probabilities of all tables no more probable than the observed table). Report EGFR-essential and tested counts for both groups, the absolute risk difference in proportions (V600E minus reference), and its descriptive 95% Newcombe hybrid-score interval without continuity correction, using Wilson limits for each group. The interval communicates precision. Only the exact test determines the verdict.

* At p ≤ 0.05, a positive difference supports the basal prediction and a negative difference contradicts it.
* Otherwise the result is inconclusive, including a zero difference. Lack of significance does not establish equivalence.
* Feasibility failure if fewer than five independent evaluable donors remain in either group, if input or outcome coding cannot be established, or if essentiality records cannot be joined unambiguously. Stop without inference if a gate fails. Do not substitute other endpoints, cohorts, scores or thresholds.

This single observational contrast cannot identify EGF as a cause, test inhibitor-induced feedback, establish the strength of BRAF knockout, or calibrate the WRN/KRAS fold-change ratios. Restriction addresses reported RAS driver mutations only. Other pathway alterations, MSI, per-model assay sensitivity and binary thresholding can affect the result.

## Reproduction and output boundary

Read record-level public research-use data locally outside the checkout. Hash the input text members after newline normalization, record source filenames and access date, and release only cohort counts, group-level statistics and original code. No record-level data, identifiers, study-specific example rows, or clinical fields are to be published.
