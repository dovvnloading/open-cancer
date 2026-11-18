# GALAXY serial ctDNA: denominator-gated false-positive occupancy

**Verdict: UNTESTABLE (denominator gate not established). Nonclinical, unreviewed research. No diagnostic, treatment, surveillance, patient-selection, or assay-validation advice.** This is an audit of published aggregate denominators, not a measurement of within-patient error clustering.

This implements **ranked new-discovery job 2** in [research/README.md](../README.md), merged in PR #88, for Nakamura et al. (2026), *npj Precision Oncology* 10:70, [DOI 10.1038/s41698-026-01277-5](https://www.nature.com/articles/s41698-026-01277-5). The question was whether independent exchangeable per-test false-positive calls with a common probability can explain the sample and patient specificity aggregates simultaneously. The article does **not** establish the exact sample-to-person correspondence needed to test that null; the required source gate stops the occupancy calculation.

## Denominator audit and stopping result

The article's Results subsection **“Association of ctDNA status in the post-definitive treatment surveillance window with DFS”**, with Fig. 3a, gives 190 people / 816 plasma samples / 1–8 tests per person; 64 ever-positive people (54 recurred), and 126 serial-negative people (**10 recurred**). From this printed patient-level partition:

| Patient-level ctDNA category | Recurred | No documented recurrence | Total |
| --- | ---: | ---: | ---: |
| Ever-positive | 54 | 10 | 64 |
| Serial-negative | 10 | 116 | 126 |
| Total | 64 | 126 | 190 |

There are **two different sets of 126**, incidentally equal in size: 126 serial-negative people, including 10 recurrences, versus 126 people without documented recurrence, including 10 ever-positive people. Similarly, the 10 recurrences among serial-negatives are **not** the 10 ever-positive nonrecurrences. The patient specificity **116/126** is algebraically consistent with the nonrecurrent column, not evidence that every one of the 126 serial-negatives was recurrence-free. The paper separately prints sample specificity **619/637**, but does not print an explicit list, cross-tab, or aggregate linkage saying that **all 637 denominator samples** were eligible surveillance samples from **those exact 126 nonrecurrent people**.

The article places both specificity numbers under the **same named post-definitive-treatment surveillance window**. This supports a common *nominal* window, **not** the additional identity of sample-level and patient-level eligibility sets. Its Methods separately describe QC/extraction exclusions, unequal follow-up, stopping monitoring at recurrence, different surveillance starts with and without ACT, and a >=7 *post-resection blood draws* requirement for nonrecurrent cohort selection. The 1–8 range is stated for the **190-person surveillance cohort**, not certified for the 637-sample subset. Fig. 3a reports a patient-level survival/recurrence partition, not which samples enter the 619/637 fraction.

Consequently, the arithmetic remainder **637 - 619 = 18** cannot safely be assigned as 18 false-positive tests among those same 126 people, even though this is one conceivable reading of the printed specificity label. Without that bridge, pairing 18 test calls with 10 nonrecurrent ever-positive patients would manufacture an unreported joint denominator. Source gates G2–G4 in [data/denominator_gate.json](data/denominator_gate.json) remain *not documented*. The prespecified rule requires **UNTESTABLE and stop**: **no allocation enumeration, occupancy tail, independent-null rejection, or minimum within-person concentration is reported**. This is an abstention, not evidence for iid testing or for clustering. See [REPORT.md](REPORT.md).

## Prior work and limited novelty screen

The literature already establishes important **mechanisms and repeated-call phenomena**, not the requested exact matched-cohort occupancy test. The [Nakamura et al. article](https://www.nature.com/articles/s41698-026-01277-5) publishes the two specificity fractions and longitudinal prognostic results, but not a denominator crosswalk or false-call concentration by person. Its Methods describe a separately trained methylation classifier and QC, not a quantitative person-level false-positive random-effects estimate. In [Colorectal Cancer Recurrence Prediction Using a Tissue-Free Epigenomic Minimal Residual Disease Assay (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11443202/), another cohort/assay analysis reports six transiently positive nonrecurrent patients with one or two positive samples; those data cannot substitute for this GALAXY cohort's test-count allocation. [Lutfi et al. (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10697744/) describe a methylation-based false-positive case that changed across repeated tests with a noncancer condition; it is a case report, not this cohort's occupancy distribution. The [ESMO Precision Medicine Working Group (2022)](https://www.annalsofoncology.org/article/S0923-7534%2822%2901721-5/fulltext) treats clonal haematopoiesis as a possible source of biological false positives; this does not establish that it caused GALAXY methylation calls. [Diergaarde et al., alpha-CORRECT (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12302971/) expressly define ever-positive surveillance specificity among nonrecurrent people and use follow-up exclusions, illustrating why the target denominator needs a definition; it is a **different** cohort.

This audit builds on [ctdna-paired-calls (PR #87)](../ctdna-paired-calls/) for frozen pre-computation restrictions, nonrejection/untestable distinctions, source provenance, and hostile counterexamples, and [ctdna-surveillance-schedule (PR #58)](../ctdna-surveillance-schedule/) for observation-noise/finite-sampling limits. Neither supplies the missing GALAXY sample mapping. The nearest unsettled testable variant remains **the conditional occupancy envelope after the denominator crosswalk is published**. It was prespecified in [PROTOCOL.md](PROTOCOL.md), not run or passed off as discovered. This is a limited targeted literature screen, not proof of bibliographic priority.

## Licence, source and privacy

**Verified CC BY 4.0** directly in the article's HTML copyright-and-licence statement at the [official PMC-hosted article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12913904/) (2026-11-18; publication 2026-01-19, Nature version of record 2026-02-17); the [Nature article](https://www.nature.com/articles/s41698-026-01277-5) is also indexed as Open Access. The canonical Nature HTML address redirected to an inaccessible publisher identity page in the search tool; the officially hosted PMC HTML and indexed Nature article text supplied the checks. **No article PDF, supplement, source data, images or browser file was downloaded, and no download prompt appeared.** The only retained source facts are sparse attributed aggregate counts with section locators, access date, DOI and licence in [data/printed_counts.csv](data/printed_counts.csv). No source prose, figure, individual sample history, patient record or identifiable data is stored. CC BY terms still require attribution; the repository's MIT licence does not substitute for the article's rights notice. [DATA_POLICY.md](../../DATA_POLICY.md) governs.

## Reproduction and protocol order

The [PROTOCOL.md](PROTOCOL.md) was **committed alone first**, remote commit [`0d39a703`](https://github.com/dovvnloading/open-cancer/commit/0d39a703653540cbd9402f84327dfecdb3c41c10), parent main [`307aa3c3`](https://github.com/dovvnloading/open-cancer/commit/307aa3c336548046c344772f3bedd4f325a1b1a0). The gate and conditional exact polynomial model were frozen before creating the data/code/report. The protocol is unchanged. Because the gate failed, the conditional statistical computation was **not** executed.

```bash
python -B research/ctdna-false-positive-clustering/analysis.py --check
python -B research/ctdna-false-positive-clustering/analysis.py --json
python -B -m unittest discover -s tests -p 'test_research_ctdna_false_positive_clustering.py' -v
```

The audit uses Python standard library only, exact rational printed fractions, line-ending-normalized SHA-256 and POSIX path keys. The independent local workspace ran the 8 dedicated tests on Python 3.13; full-repository validation and Python 3.12/3.14 are delegated to CI where configured. Local `git clone` could not reach github.com due unavailable DNS; the authorized GitHub connection created the branch and commits from the confirmed current main SHA. This does not imply a full local clone was obtained. See the [hostile review and replication limits](REPORT.md).

**Interpretation boundary:** Even with a reconciled denominator and a robust rejection, occupancy concentration cannot identify clonal haematopoiesis, persistent assay background, true subclinical disease that cleared, contamination, follow-up bias or another cause. No clinical or assay recommendation follows.
