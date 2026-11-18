# Source corrections and analysis-plan clarifications

**Nonclinical, unreviewed. No treatment advice.** The original PROTOCOL.md is
unchanged after its local pre-extraction commit `2d46668`. Entries below occurred
during extraction or hostile review, before final analysis freeze and push.
They document changes instead of retrospectively editing the plan.

## Post-extraction scope clarification

The protocol's P6b universal phrase says "improve control" but does not
explicitly settle whether overall survival qualifies or impose one confidence
level on descriptive published intervals. This ambiguity was recognized after
encountering COIN. We preserve **both** interpretations:

- The broad clinical-outcome/any-published-interval reading meets the original
  adverse-interval rule on COIN OS at the source's **80%** confidence level.
  The corresponding `falsified` label describes that limited rule only; it is
  not 95% evidence, a tumour-control result or a PR #61 falsification.
- The explicitly **post-extraction** tumour-control-only sensitivity does not
  meet the falsification rule. Nor does a further sensitivity requiring a
  published 95% tumour-control contrast interval. Both are `untestable` in the
  sense of inconclusive under the allowed three-label vocabulary.

The endpoint ambiguity prevents a unique fully pre-specified P6b verdict.
The narrower readings must not erase the adverse broad result or be described
as wholly pre-specified. No 95% interval is reconstructed from an 80% interval
or p-value, and failure of noninferiority alone is not treated as proof of harm.
Negative point estimates remain in all tables. P1-P5 and the no-pooling rule
are unchanged.

## Transcription and interpretation corrections

| Item | Primary locator | Correction retained |
| --- | --- | --- |
| Zhang 2017 versus later update | [Author-hosted version of record](https://jessica-cunningham.com/wp-content/uploads/2019/11/cunningham_et_al-2017-nature_communications.pdf), clinical results and Figure 5 | Preserve a lower-bound interim median and the comparator printed in that report. Do not substitute the later control median. Same trial, not replication. |
| Zhang 2022 HR | [Article](https://doi.org/10.7554/eLife.76284), clinical-results paragraph and Figure 1 | The reported HR with confidence interval is OS, not radiographic TTP. |
| Zhang 2022 exposure inconsistency | Same article, clinical-results paragraph beginning "All 16 patients" versus Discussion paragraph beginning "Nevertheless" | Results defines dosing rate as 54% of SOC and off-treatment time as 46%; Discussion calls time off 54%. Prefer the explicit Results definition; keep both source locations and flag the internal contradiction. These are rates/time fractions, not lifetime cumulative-dose fractions. |
| Zhang 2022 rate scaling | Same article, equations and Appendix 2 "Optimized parameters" | The fitted growth multiplier multiplies the PSA-slope rates. Conditional inverse growth must include it. An arbitrary carrying-capacity scale cannot identify cell number or mutation supply. |
| INTERIM registration | [Accepted manuscript](https://eprints.whiterose.ac.uk/id/eprint/206919/7/EJC-D-23-02078_R1.pdf), title page and Methods | Use ISRCTN18183156 / EudraCT 2016-005228-27. An initially unverified NCT identifier was removed before integration. |
| COIN intervals | [Article](https://pmc.ncbi.nlm.nih.gov/articles/PMC3159416/), Results preceding Figure 3 | Retain 80% HR CI and median-arm IQRs. Neither is a 95% CI. Use ITT primary comparison; do not select a favourable subgroup. |
| Pressley lung xenograft | [Article](https://doi.org/10.1038/s41598-026-63987-x), Results "Continuous therapy extends survival in NSCLC model" | Negative mean-survival comparison retained; omnibus Kruskal-Wallis result is not a pairwise confidence interval. Absence of density-dependence evidence does not measure zero competition. |
| Reuse | Article-specific license blocks, including official Europe PMC XML where publisher access challenged | PMC availability is not a license. Only explicitly permitted aggregates are transcribed. Unknown/ND rights remain metadata, with positive and negative sources affected. |
| Competition prediction | PR #61 establishment integral and comparison coefficient | The planned P2 is auxiliary. At saturated c, increasing competition reduces A and decreases the advantage. Do not attribute universal increasing-competition benefit to the theorem. |

This audit is not external scientific peer review. Every numerical correction
is backed by a primary locator; agreement between AI workers is not used as
validation. Source data were never changed to make a hypothesis pass.

## Publication transport

The frozen protocol was committed locally as
`2d46668a276568498f2f8d945a882734a8581ed7` at 22:34:19 on 17 November 2026
(America/New_York), before extraction. The command-line push later failed
because Git had no authenticated credential. The authenticated GitHub connector
published the same protocol and audited file contents on a branch from the
then-current main. Its API creates new commit identities and publication
timestamps. The local protocol SHA therefore records the local sequence; it
must not be presented as an independently timestamped remote preregistration.
The protocol's newline-normalized content hash remains fixed in the test.
