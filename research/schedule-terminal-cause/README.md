# Chemotherapy schedules and attributed terminal causes

**Nonclinical preclinical aggregate analysis. Research status: unreviewed.
No treatment, diagnosis, schedule or dose advice.** This is a computation on
published mouse-xenograft summaries, not a new experiment, an independent
biological replication, or clinical guidance.

## Verdict

**The toxicity-attribution contrast fails the frozen robustness rule.**
The designated primary comparison is constant combination MTD minus
source-pooled ping-pong adaptive therapy in Seyedi et al.'s Table S4. Its
printed toxicity-label difference is large and its unadjusted exact tests
reject, but neither passes the prespecified multiplicity correction. The
reclassification analysis separately measures how little attribution change
removes the unadjusted evidence and how many changed labels remove the
direction. Every tandem and vehicle comparison, including null results,
remains in the report below.

This verdict means that the reported composition does **not establish an
attribution-robust toxicity contrast under this declared rule**. It is not
evidence of no toxicity, no schedule association or no survival difference.
The full uncertainty envelope necessarily permits every toxicity label to be
reassigned to cancer; the finite breakpoints and the result at the printed
labels are the informative findings, rather than that limiting observation.
The single-drug adaptive comparisons are **unconditionally untestable** because
their printed cause totals do not reconcile to their allocations. Their
conditional sensitivity results do not repair that failure.

This completes the aggregate question in ranked job 5 of the
[research index](../README.md) and quantifies the comparator warning in
[PR #66's evidence note](../containment-evidence/README.md). It does not test
the mutation-supply boundary in
[PR #61's tumour-scale theorem](../containment-at-tumour-scale/README.md).

## Frozen specification and execution

[PROTOCOL.md](PROTOCOL.md) was committed locally by itself as
`dbf3c54746ba853a81a65f177108531081460978` before inferential code was written
or contrasts, p-values, intervals and fragility thresholds were calculated.
Its content has not been amended. The branch starts from `main` commit
`37c37e2b4cfc711ce938201123ef3085b8868f2f`. Reading the public source and
checking its denominators preceded the freeze; this is transparent
prespecification of a new computation on known results, not blinding or
independent prospective registration of the experiment.

The command-line push lacked GitHub credentials, so publication used the
authenticated GitHub connector. It recreated the commit objects: the
[remote protocol-only commit](https://github.com/dovvnloading/open-cancer/commit/f588ce356c92aba00be33faf3e350af7382587d2)
publishes the identical protocol tree, followed by the analysis commit.
The original local freeze SHA is an audit identifier, **not a remote GitHub
ref**. Its recorded commit time is 2026-11-18 03:46:35 America/New_York;
the remote publication time is later and is not presented as the freeze time.
The unchanged protocol hash, shared tree and original commit metadata make
the correspondence explicit without claiming a remotely timestamped
preregistration.

<details>
<summary>Original local protocol commit object</summary>

The following is the original commit payload (with a final newline). Git's
commit-object hashing reproduces the local freeze identifier. The tree is
also the remote protocol-only commit's tree; the parent is the starting main.
These metadata document local execution history, not an independently trusted
time attestation.

```text
tree 3acce48a734556bbfa18b821d8fd605706f03763
parent 37c37e2b4cfc711ce938201123ef3085b8868f2f
author Codex <codex@openai.com> 1794991595 -0500
committer Codex <codex@openai.com> 1794991595 -0500

research: freeze terminal-cause attribution protocol before inference
```

</details>

The primary estimand is the **composition of attributed terminal events at
the animals' respective observed endpoints**. The direction is MTD minus
ping-pong in the fraction labelled toxicity. The protocol requires that
positive direction and both corrected exact tests to hold for every allowed
attribution state. The same printed-direction retention rule is reported
for secondary pairs. No secondary result can replace a failed primary result.

Run from the repository root, using Python's standard library only:

```bash
python research/schedule-terminal-cause/analysis.py
python research/schedule-terminal-cause/analysis.py --json
python research/schedule-terminal-cause/analysis.py --check
python -m unittest discover -s tests -p test_research_schedule_terminal_cause.py -v
```

The default output reproduces the generated quantitative section below.
The JSON output retains exact rational results and complete sensitivity
summaries and witnesses. `--check` compares the generated section with this
file. The offline test independently recalculates the numerical analysis and
checks the provenance hashes; no source retrieval is part of testing.

## Sources, access and reuse

Seyedi S, Teo R, Foster L, et al. (2024). *Testing Adaptive Therapy Protocols
Using Gemcitabine and Capecitabine in a Preclinical Model of Endocrine-Resistant
Breast Cancer*. **Cancers 16:257**.
[Publisher article](https://www.mdpi.com/2072-6694/16/2/257),
[DOI](https://doi.org/10.3390/cancers16020257),
[public supplement](https://www.mdpi.com/article/10.3390/cancers16020257/s1).
Access date: **2026-11-18 UTC**. Source versions are the publisher's final
article and the PDF identified by the binary hashes in the protocol.

Rights were checked before numerical transcription. The article's closing
copyright paragraph assigns copyright to the authors and grants
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The publisher's [open-access policy](https://www.mdpi.com/openaccess),
"Meaning of Open Access", expressly extends reuse to supplementary material,
while excepting separately credited third-party content. The supplement's
extracted text has no separate licence notice, and visual inspection found
no different licence or third-party credit beside the tables used here.
**The supplement reuse basis is the article's grant plus that explicit
publisher policy; the PDF itself does not print a separate CC notice.**

Changes are selection of printed aggregates, transcription into CSV,
standardization of labels and units, and the declared calculations. The source
data retain their CC BY terms; this repository's MIT licence does not relicense
them. Attribution and locators accompany every CSV row, following
[DATA_POLICY.md](../../DATA_POLICY.md). The PDF, images and individual animal
records are not redistributed.

The supplement was read through ordinary HTTP in the execution environment.
`urllib.request` received the public ZIP, and `zipfile` read its sole member,
`cancers-2752834-supplementary.pdf`. Text extraction and rendered-page inspection
were used together. **No browser file download or browser interaction was
used; no download prompt appeared.** The requested URL redirected to
[the publisher's article-specific supplement endpoint](https://www.mdpi.com/2072-6694/16/2/257/s1).
An independent repeat HTTP fetch returned an identical PDF hash. This is an
integrity check on the retrieved file, not validation of its scientific claims.

| Input | What was transcribed | Source locator |
| --- | --- | --- |
| [counts.csv](data/counts.csv) | Only the printed toxicity, cancer and other counts; unreconciled rows retained verbatim | Tables S2/S3, supplement p.32; S4, pp.32-33; row and individual column names in CSV |
| [allocations.csv](data/allocations.csv) | The stated per-group allocations applied to the listed treatment arms; explicit joins into the source pools | Article Table 1 p.4; section 2.3.1 p.5 |
| [doses.csv](data/doses.csv) | Printed mean cumulative and mean daily doses, for context only | Table S5, supplement p.34 |
| [dose_conflicts.csv](data/dose_conflicts.csv) | Conflicting gemcitabine cumulative-dose means in the article | Article section 3.3.1 p.11 |

The allocation entries apply the article's per-treatment-group and vehicle
sample-size statements to named Table 1 arms; they are not individually
printed sample-size cells. Table 1's defective numeric row labels are avoided
by recording its full treatment-group names. The pooled adaptive rows are
kept intact. The single combination MTD group and single vehicle group are
stored once and reused as shared comparators; they are not replicated animals.

## Denominator and attribution audit

The generated reconciliation table retains both the sum of printed cause
assignments and the allocated size, with the discrepancy visible. The
single-drug adaptive totals are not valid randomized denominators. We neither
choose a cause to subtract nor use Table S1 to manufacture corrected counts.
Only S4 supplies a reconciled primary comparison.

The conditional single-drug calculation asks a narrower hypothetical question:
**if exactly one adaptive cell is one too high, and everything else is
correct, what range results over every possible nonnegative decrement?**
That model excludes other reporting errors and is not evidence for any one
completion. The original values remain unchanged in `counts.csv`. The
conditional comparisons use their own MTD groups, a separate exploratory
Bonferroni family, and cannot rescue the primary verdict.

Table S1's "Reason of euthanizing mouse" column, supplement pp.11-31, and the
cause-definition narrative on p.32 inform the *types* of uncertainty. They
are not an individual-record dataset for this analysis. Relevant reasons
include weight loss, nonspecific weakness, spontaneous death, surgery,
caliper tumour endpoints, and interference with walking. The combination MTD
reasons appear on pp.21-23, the surgery reason on p.29, and the vehicle reasons
on p.31. No additional count of ambiguous cases was reconstructed from those
rows.

The p.32 narrative assigns spontaneous deaths, weight loss and other signs of
ill health to toxicity by assumption. It explicitly asks whether weight loss
might instead be cancer cachexia and says necropsies did not distinguish the
possibilities. Walking-related sacrifices with low viable burden were put in
"something else", though tumour bulk could still cause the endpoint. Tumour
burden triggers were assigned to cancer. The untreated group's toxicity label
is preserved and is a direct warning against interpreting that column as
demonstrated chemotherapy-caused death.

The source narrative refers to cause "Tables S3-S5", whereas their actual
headings are S2-S4. Locators here use the printed headings. The S5-versus-main-
text dose discrepancies are shown rather than silently corrected. No dose
mean is used as a covariate, denominator, observation-time estimate or reason
to alter a terminal label.

## Statistical reference model

Each pair has a conditional exact composition test over toxicity/cancer/other
and a conditional exact toxicity-versus-rest test. The code enumerates all
tables with the observed margins and sums probabilities no larger than the
observed table probability, including exact ties. Integer weights and rational
sums avoid floating-point decisions about inclusion in a tail. Empty columns
remain valid. There is no asymptotic approximation or simulated p-value.

Holm correction is applied to the entire fixed S4 family in every coherent
four-group attribution state. A shared control induces correlated comparisons;
Holm's familywise error bound does not require their independence. The small
sample, broad family and two endpoints make this correction conservative.
Showing unadjusted values makes that cost visible; it does not change the
frozen verdict after seeing the data.

Simultaneous Clopper-Pearson category intervals use the protocol's Bonferroni
allocation. Subtracting arm intervals gives conservative simultaneous
difference envelopes. Those envelopes and the Holm exact tests are different
procedures, so their decisions need not coincide. The analysis assumes an
independent, exchangeable categorical sampling reference model within each
source pool. Intermittent and dose-modulation protocols within a pool could
have different cause probabilities; dependence and heterogeneity are not
identified by these aggregates. Thus exactness is conditional on a reference
model, not a design-based causal guarantee.

### What the sensitivity ranges mean

The **T-only** envelope permits any printed toxicity label to become cancer.
The **T-and-O** envelope additionally permits any other label to become cancer.
All original cancer labels stay fixed, making this an explicitly asymmetric
stress model. The second envelope is deliberately broad: technical and surgery
events are not thereby shown to be cancer deaths. Neither envelope estimates
a misclassification rate, establishes the true cause, or assigns probabilities
to the enumerated states.

Each changed label counts once, with group totals held fixed. The full
four-row state is evaluated before applying Holm. Pair-local breakpoints hold
the other groups fixed and use pair denominators. Global budgets count edits
anywhere and use the full S4 denominator, because changes outside a pair may
alter its adjusted p-value through multiplicity. "At most" budgets include
all smaller changes; no monotonicity of an individual p-value is assumed.

Loss of raw rejection, loss of corrected rejection and loss of the point
direction are separate outputs. A breakpoint of zero means that the stated
criterion already fails at the printed labels. A tied printed direction is
not a positive association. Fractions of animals and of eligible labels answer
different questions and are shown separately. Reclassification envelopes are
not confidence intervals; the interval envelopes add reference-model sampling
uncertainty to each permitted attribution state.

<!-- BEGIN GENERATED RESULTS -->

## Recomputed results

**Primary verdict: FAILS.** The printed-label primary rule already fails; attribution uncertainty is evaluated below. Failure is failure of the frozen robustness rule, not evidence of no toxicity or no association.

All differences are percentage points (pp). Probabilities are displayed to six decimals; calculations and decision comparisons use exact rational p-values. T/C/O means toxicity/cancer/other attributed terminal events.

### Denominator audit

| Source | Group | Printed T/C/O | Printed total | Allocated | Excess |
| --- | --- | --- | --- | --- | --- |
| S2 p.32 | CAP adaptive | (4,5,4) | 13 | 12 | 1 |
| S2 p.32 | CAP MTD | (2,4,0) | 6 | 6 | 0 |
| S3 p.32 | GEM adaptive | (7,1,5) | 13 | 12 | 1 |
| S3 p.32 | GEM MTD | (5,1,0) | 6 | 6 | 0 |
| S4 p.32 | Vehicle | (1,3,0) | 4 | 4 | 0 |
| S4 p.33 | Ping-pong adaptive | (3,6,3) | 12 | 12 | 0 |
| S4 p.33 | Tandem adaptive | (5,4,3) | 12 | 12 | 0 |
| S4 p.33 | Combination MTD | (6,0,0) | 6 | 6 | 0 |

Across the unique printed rows: 72 assignments versus 70 allocated animals. S4 uses 34 unique animals; its MTD and vehicle rows occur once. S2/S3 remain unreconciled.

### S4 compositions and simultaneous reference intervals

Intervals have at least 95% simultaneous coverage for the twelve category probabilities under the frozen exchangeable-binomial reference model; each marginal uses error 1/240.

| Group | n | T/C/O (%) | T interval (%) | C interval (%) | O interval (%) |
| --- | --- | --- | --- | --- | --- |
| Vehicle | 4 | 25.000000/75.000000/0.000000 | [0.052124, 91.781816] | [8.218184, 99.947876] | [0.000000, 78.635650] |
| Ping-pong adaptive | 12 | 25.000000/50.000000/25.000000 | [2.224675, 69.196776] | [12.859021, 87.140979] | [2.224675, 69.196776] |
| Tandem adaptive | 12 | 41.666667/33.333333/25.000000 | [8.476989, 81.980940] | [4.904443, 76.025173] | [2.224675, 69.196776] |
| Combination MTD | 6 | 100.000000/0.000000/0.000000 | [35.737689, 100.000000] | [0.000000, 64.262311] | [0.000000, 64.262311] |

### All fixed S4 comparisons at printed labels

| First minus second | T difference (pp) | C difference (pp) | O difference (pp) | Composition p | Holm p | T vs rest p | Holm p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MTD minus ping-pong (primary) | 75.000000 | -50.000000 | -25.000000 | 0.016484 | 0.181319 | 0.009050 | 0.108597 |
| MTD minus tandem | 58.333333 | -33.333333 | -25.000000 | 0.061409 | 0.429864 | 0.037707 | 0.333333 |
| MTD minus vehicle | 75.000000 | -75.000000 | 0.000000 | 0.033333 | 0.333333 | 0.033333 | 0.333333 |
| Tandem minus ping-pong | 16.666667 | -16.666667 | 0.000000 | 0.869534 | 1.000000 | 0.666846 | 1.000000 |
| Ping-pong minus vehicle | 0.000000 | -25.000000 | 25.000000 | 0.762637 | 1.000000 | 1.000000 | 1.000000 |
| Tandem minus vehicle | 16.666667 | -41.666667 | 25.000000 | 0.446154 | 1.000000 | 1.000000 | 1.000000 |

| First minus second | T difference interval (pp) | C difference interval (pp) | O difference interval (pp) |
| --- | --- | --- | --- |
| MTD minus ping-pong (primary) | [-33.459087, 97.775325] | [-87.140979, 51.403291] | [-69.196776, 62.037637] |
| MTD minus tandem | [-46.243251, 91.523011] | [-76.025173, 59.357869] | [-69.196776, 62.037637] |
| MTD minus vehicle | [-56.044127, 99.947876] | [-99.947876, 56.044127] | [-78.635650, 64.262311] |
| Tandem minus ping-pong | [-60.719787, 79.756265] | [-82.236536, 63.166152] | [-66.972101, 66.972101] |
| Ping-pong minus vehicle | [-89.557141, 69.144652] | [-87.088855, 78.922795] | [-76.410975, 69.196776] |
| Tandem minus vehicle | [-83.304826, 81.928816] | [-95.043433, 67.806989] | [-76.410975, 69.196776] |

The intervals are conservative simultaneous bounds, not inversions of the Holm tests. Shared comparisons are dependent; no independent-replication claim is made.

### Attribution envelope: T-only

336 coherent integer states; 15 eligible labels among 34 animals. The primary verdict is **fails**. These are deterministic sensitivity ranges, not confidence distributions.

| Contrast | T difference range (pp) | C difference range (pp) | O difference range (pp) | Robustness verdict |
| --- | --- | --- | --- | --- |
| MTD minus ping-pong (primary) | [-25.000000, 100.000000] | [-75.000000, 50.000000] | [-25.000000, -25.000000] | fails |
| MTD minus tandem | [-41.666667, 100.000000] | [-75.000000, 66.666667] | [-25.000000, -25.000000] | fails |
| MTD minus vehicle | [-25.000000, 100.000000] | [-100.000000, 25.000000] | [0.000000, 0.000000] | fails |
| Tandem minus ping-pong | [-25.000000, 41.666667] | [-41.666667, 25.000000] | [0.000000, 0.000000] | fails |
| Ping-pong minus vehicle | [-25.000000, 25.000000] | [-50.000000, 0.000000] | [25.000000, 25.000000] | fails |
| Tandem minus vehicle | [-25.000000, 41.666667] | [-66.666667, 0.000000] | [25.000000, 25.000000] | fails |

| Contrast | Composition p range | Composition Holm range | T p range | T Holm range |
| --- | --- | --- | --- | --- |
| MTD minus ping-pong (primary) | [0.000054, 0.760019] | [0.000646, 1.000000] | [0.000054, 1.000000] | [0.000646, 1.000000] |
| MTD minus tandem | [0.000054, 0.760019] | [0.000646, 1.000000] | [0.000054, 1.000000] | [0.000646, 1.000000] |
| MTD minus vehicle | [0.004762, 1.000000] | [0.038095, 1.000000] | [0.004762, 1.000000] | [0.038095, 1.000000] |
| Tandem minus ping-pong | [0.041786, 1.000000] | [0.298137, 1.000000] | [0.037267, 1.000000] | [0.298137, 1.000000] |
| Ping-pong minus vehicle | [0.365385, 0.777473] | [1.000000, 1.000000] | [0.250000, 1.000000] | [1.000000, 1.000000] |
| Tandem minus vehicle | [0.081319, 0.777473] | [0.325275, 1.000000] | [0.244505, 1.000000] | [0.733516, 1.000000] |

Envelope of simultaneous difference intervals (minimum lower bound to maximum upper bound):

| Contrast | T (pp) | C (pp) | O (pp) |
| --- | --- | --- | --- |
| MTD minus ping-pong (primary) | [-69.196776, 100.000000] | [-97.775325, 87.140979] | [-69.196776, 62.037637] |
| MTD minus tandem | [-81.980940, 100.000000] | [-97.775325, 95.095557] | [-69.196776, 62.037637] |
| MTD minus vehicle | [-91.781816, 100.000000] | [-100.000000, 91.781816] | [-78.635650, 64.262311] |
| Tandem minus ping-pong | [-69.196776, 81.980940] | [-92.870882, 84.916304] | [-66.972101, 66.972101] |
| Ping-pong minus vehicle | [-91.781816, 69.196776] | [-87.140979, 89.557141] | [-76.410975, 69.196776] |
| Tandem minus vehicle | [-91.781816, 81.980940] | [-95.095557, 89.557141] | [-76.410975, 69.196776] |

#### Pair-local fragility

Other rows remain at their printed labels. Breakpoint zero means the printed comparison already lacks rejection or strict direction. Fractions below use the pair's animals and eligible labels, respectively. Witness row order is vehicle / ping-pong / tandem / combination MTD, each T,C,O.

| Contrast | Loss criterion | Relabels | Fraction of pair animals | Fraction of eligible pair labels | Witness table |
| --- | --- | --- | --- | --- | --- |
| MTD minus ping-pong (primary) | raw_composition | 1 | 1/18 (5.555556%) | 1/9 (11.111111%) | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus ping-pong (primary) | raw_toxicity | 2 | 2/18 (11.111111%) | 2/9 (22.222222%) | (1,3,0)/(3,6,3)/(5,4,3)/(4,2,0) |
| MTD minus ping-pong (primary) | holm_composition | 0 | 0/18 (0.000000%) | 0/9 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | holm_toxicity | 0 | 0/18 (0.000000%) | 0/9 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | direction | 5 | 5/18 (27.777778%) | 5/9 (55.555556%) | (1,3,0)/(3,6,3)/(5,4,3)/(1,5,0) |
| MTD minus ping-pong (primary) | combined | 0 | 0/18 (0.000000%) | 0/9 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus tandem | raw_composition | 0 | 0/18 (0.000000%) | 0/11 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus tandem | raw_toxicity | 1 | 1/18 (5.555556%) | 1/11 (9.090909%) | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus tandem | holm_composition | 0 | 0/18 (0.000000%) | 0/11 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus tandem | holm_toxicity | 0 | 0/18 (0.000000%) | 0/11 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus tandem | direction | 4 | 4/18 (22.222222%) | 4/11 (36.363636%) | (1,3,0)/(3,6,3)/(5,4,3)/(2,4,0) |
| MTD minus tandem | combined | 0 | 0/18 (0.000000%) | 0/11 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | raw_composition | 1 | 1/10 (10.000000%) | 1/7 (14.285714%) | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus vehicle | raw_toxicity | 1 | 1/10 (10.000000%) | 1/7 (14.285714%) | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus vehicle | holm_composition | 0 | 0/10 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | holm_toxicity | 0 | 0/10 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | direction | 5 | 5/10 (50.000000%) | 5/7 (71.428571%) | (1,3,0)/(3,6,3)/(5,4,3)/(1,5,0) |
| MTD minus vehicle | combined | 0 | 0/10 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | raw_composition | 0 | 0/24 (0.000000%) | 0/8 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | raw_toxicity | 0 | 0/24 (0.000000%) | 0/8 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | holm_composition | 0 | 0/24 (0.000000%) | 0/8 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | holm_toxicity | 0 | 0/24 (0.000000%) | 0/8 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | direction | 2 | 2/24 (8.333333%) | 2/8 (25.000000%) | (1,3,0)/(3,6,3)/(3,6,3)/(6,0,0) |
| Tandem minus ping-pong | combined | 0 | 0/24 (0.000000%) | 0/8 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | raw_composition | 0 | 0/16 (0.000000%) | 0/4 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | raw_toxicity | 0 | 0/16 (0.000000%) | 0/4 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | holm_composition | 0 | 0/16 (0.000000%) | 0/4 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | holm_toxicity | 0 | 0/16 (0.000000%) | 0/4 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | direction | 0 | 0/16 (0.000000%) | 0/4 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | combined | 0 | 0/16 (0.000000%) | 0/4 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | raw_composition | 0 | 0/16 (0.000000%) | 0/6 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | raw_toxicity | 0 | 0/16 (0.000000%) | 0/6 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | holm_composition | 0 | 0/16 (0.000000%) | 0/6 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | holm_toxicity | 0 | 0/16 (0.000000%) | 0/6 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | direction | 2 | 2/16 (12.500000%) | 2/6 (33.333333%) | (1,3,0)/(3,6,3)/(3,6,3)/(6,0,0) |
| Tandem minus vehicle | combined | 0 | 0/16 (0.000000%) | 0/6 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |

#### Primary coherent global fragility

Here edits may affect any row, including changes outside the primary pair that alter the Holm ordering. Fractions use the whole S4 table.

| Loss criterion | Relabels | Fraction of all animals | Fraction of all eligible labels | Witness table |
| --- | --- | --- | --- | --- |
| raw_composition | 1 | 1/34 (2.941176%) | 1/15 (6.666667%) | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| raw_toxicity | 2 | 2/34 (5.882353%) | 2/15 (13.333333%) | (1,3,0)/(3,6,3)/(5,4,3)/(4,2,0) |
| holm_composition | 0 | 0/34 (0.000000%) | 0/15 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| holm_toxicity | 0 | 0/34 (0.000000%) | 0/15 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| direction | 5 | 5/34 (14.705882%) | 5/15 (33.333333%) | (1,3,0)/(3,6,3)/(5,4,3)/(1,5,0) |
| combined | 0 | 0/34 (0.000000%) | 0/15 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |

#### Primary worst case at every cumulative global budget

Each row permits at most k changed labels; extrema can be attained by different tables and are not a synthetic joint table. Complete ranges and attaining witnesses for every metric and every pair at each budget are available with `--json`.

| k | States | k/all animals (%) | k/eligible (%) | T difference range (pp) | Worst composition p | Worst composition Holm | Worst T p | Worst T Holm |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1 | 0.000000 | 0.000000 | [75.000000, 75.000000] | 0.016484 | 0.181319 | 0.009050 | 0.108597 |
| 1 | 5 | 2.941176 | 6.666667 | [58.333333, 83.333333] | 0.085973 | 0.945701 | 0.042986 | 0.515837 |
| 2 | 14 | 5.882353 | 13.333333 | [41.666667, 91.666667] | 0.229261 | 1.000000 | 0.141403 | 1.000000 |
| 3 | 30 | 8.823529 | 20.000000 | [25.000000, 100.000000] | 0.456044 | 1.000000 | 0.344053 | 1.000000 |
| 4 | 54 | 11.764706 | 26.666667 | [8.333333, 100.000000] | 0.636393 | 1.000000 | 1.000000 | 1.000000 |
| 5 | 86 | 14.705882 | 33.333333 | [-8.333333, 100.000000] | 0.636393 | 1.000000 | 1.000000 | 1.000000 |
| 6 | 125 | 17.647059 | 40.000000 | [-25.000000, 100.000000] | 0.760019 | 1.000000 | 1.000000 | 1.000000 |
| 7 | 168 | 20.588235 | 46.666667 | [-25.000000, 100.000000] | 0.760019 | 1.000000 | 1.000000 | 1.000000 |
| 8 | 211 | 23.529412 | 53.333333 | [-25.000000, 100.000000] | 0.760019 | 1.000000 | 1.000000 | 1.000000 |
| 9 | 250 | 26.470588 | 60.000000 | [-25.000000, 100.000000] | 0.760019 | 1.000000 | 1.000000 | 1.000000 |
| 10 | 282 | 29.411765 | 66.666667 | [-25.000000, 100.000000] | 0.760019 | 1.000000 | 1.000000 | 1.000000 |
| 11 | 306 | 32.352941 | 73.333333 | [-25.000000, 100.000000] | 0.760019 | 1.000000 | 1.000000 | 1.000000 |
| 12 | 322 | 35.294118 | 80.000000 | [-25.000000, 100.000000] | 0.760019 | 1.000000 | 1.000000 | 1.000000 |
| 13 | 331 | 38.235294 | 86.666667 | [-25.000000, 100.000000] | 0.760019 | 1.000000 | 1.000000 | 1.000000 |
| 14 | 335 | 41.176471 | 93.333333 | [-25.000000, 100.000000] | 0.760019 | 1.000000 | 1.000000 | 1.000000 |
| 15 | 336 | 44.117647 | 100.000000 | [-25.000000, 100.000000] | 0.760019 | 1.000000 | 1.000000 | 1.000000 |

<details>
<summary>Attaining witnesses for all T-only extrema</summary>

Differences and interval endpoints are pp; p-values are probabilities. Each extremum has its own coherent table. Interval lower/upper extrema are endpoints of reference-model uncertainty, not additional tests.

| Contrast | Metric | Minimum | Minimum witness | Maximum | Maximum witness |
| --- | --- | --- | --- | --- | --- |
| MTD minus ping-pong (primary) | composition_p | 0.000054 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 0.760019 | (1,3,0)/(2,7,3)/(5,4,3)/(1,5,0) |
| MTD minus ping-pong (primary) | toxicity_p | 0.000054 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(2,4,0) |
| MTD minus ping-pong (primary) | toxicity_difference | -25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | toxicity_interval_lower | -69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | -4.481304 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | toxicity_interval_upper | 62.037637 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | cancer_difference | -75.000000 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 50.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus ping-pong (primary) | cancer_interval_lower | -97.775325 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | -51.403291 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus ping-pong (primary) | cancer_interval_upper | 33.459087 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 87.140979 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus ping-pong (primary) | other_difference | -25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | -25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | other_interval_lower | -69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | -69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | other_interval_upper | 62.037637 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 62.037637 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | composition_holm | 0.000646 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(4,2,0) |
| MTD minus ping-pong (primary) | toxicity_holm | 0.000646 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(4,2,0) |
| MTD minus tandem | composition_p | 0.000054 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 0.760019 | (1,3,0)/(3,6,3)/(2,7,3)/(1,5,0) |
| MTD minus tandem | toxicity_p | 0.000054 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(3,3,0) |
| MTD minus tandem | toxicity_difference | -41.666667 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| MTD minus tandem | toxicity_interval_lower | -81.980940 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | -4.481304 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| MTD minus tandem | toxicity_interval_upper | 55.785322 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| MTD minus tandem | cancer_difference | -75.000000 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 66.666667 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus tandem | cancer_interval_lower | -97.775325 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | -40.287484 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus tandem | cancer_interval_upper | 33.459087 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 95.095557 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus tandem | other_difference | -25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | -25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus tandem | other_interval_lower | -69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | -69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus tandem | other_interval_upper | 62.037637 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 62.037637 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus tandem | composition_holm | 0.000646 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus tandem | toxicity_holm | 0.000646 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus vehicle | composition_p | 0.004762 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(2,4,0) |
| MTD minus vehicle | toxicity_p | 0.004762 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(2,4,0) |
| MTD minus vehicle | toxicity_difference | -25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | toxicity_interval_lower | -91.781816 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | -42.897961 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | toxicity_interval_upper | 64.210187 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | cancer_difference | -100.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus vehicle | cancer_interval_lower | -100.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | -64.210187 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus vehicle | cancer_interval_upper | 42.897961 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 91.781816 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus vehicle | other_difference | 0.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 0.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | other_interval_lower | -78.635650 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | -78.635650 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | other_interval_upper | 64.262311 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 64.262311 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | composition_holm | 0.038095 | (0,4,0)/(2,7,3)/(2,7,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus vehicle | toxicity_holm | 0.038095 | (0,4,0)/(2,7,3)/(2,7,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| Tandem minus ping-pong | composition_p | 0.041786 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(4,5,3)/(6,0,0) |
| Tandem minus ping-pong | toxicity_p | 0.037267 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(4,5,3)/(6,0,0) |
| Tandem minus ping-pong | toxicity_difference | -25.000000 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 41.666667 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | toxicity_interval_lower | -69.196776 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | -31.742003 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | toxicity_interval_upper | 37.994318 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 81.980940 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | cancer_difference | -41.666667 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 25.000000 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| Tandem minus ping-pong | cancer_interval_lower | -92.870882 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | -56.337755 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| Tandem minus ping-pong | cancer_interval_upper | 45.221949 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 84.916304 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| Tandem minus ping-pong | other_difference | 0.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 0.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | other_interval_lower | -66.972101 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | -66.972101 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | other_interval_upper | 66.972101 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 66.972101 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | composition_holm | 0.298137 | (0,4,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | toxicity_holm | 0.298137 | (0,4,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | composition_p | 0.365385 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 0.777473 | (1,3,0)/(2,7,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | toxicity_p | 0.250000 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | toxicity_difference | -25.000000 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 25.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | toxicity_interval_lower | -91.781816 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | -76.410975 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | toxicity_interval_upper | 40.166868 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 69.196776 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | cancer_difference | -50.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 0.000000 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | cancer_interval_lower | -87.140979 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | -69.144652 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | cancer_interval_upper | 65.776629 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 89.557141 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | other_difference | 25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | other_interval_lower | -76.410975 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | -76.410975 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | other_interval_upper | 69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | composition_holm | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | toxicity_holm | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | composition_p | 0.081319 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 0.777473 | (1,3,0)/(3,6,3)/(2,7,3)/(6,0,0) |
| Tandem minus vehicle | toxicity_p | 0.244505 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | toxicity_difference | -25.000000 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 41.666667 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | toxicity_interval_lower | -91.781816 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | -70.158660 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | toxicity_interval_upper | 40.166868 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 81.980940 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | cancer_difference | -66.666667 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 0.000000 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| Tandem minus vehicle | cancer_interval_lower | -95.095557 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | -69.144652 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| Tandem minus vehicle | cancer_interval_upper | 54.660823 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 89.557141 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| Tandem minus vehicle | other_difference | 25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | other_interval_lower | -76.410975 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | -76.410975 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | other_interval_upper | 69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | composition_holm | 0.325275 | (0,4,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | toxicity_holm | 0.733516 | (0,4,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |

</details>

### Attribution envelope: T-and-O

5376 coherent integer states; 21 eligible labels among 34 animals. The primary verdict is **fails**. These are deterministic sensitivity ranges, not confidence distributions.

| Contrast | T difference range (pp) | C difference range (pp) | O difference range (pp) | Robustness verdict |
| --- | --- | --- | --- | --- |
| MTD minus ping-pong (primary) | [-25.000000, 100.000000] | [-100.000000, 50.000000] | [-25.000000, 0.000000] | fails |
| MTD minus tandem | [-41.666667, 100.000000] | [-100.000000, 66.666667] | [-25.000000, 0.000000] | fails |
| MTD minus vehicle | [-25.000000, 100.000000] | [-100.000000, 25.000000] | [0.000000, 0.000000] | fails |
| Tandem minus ping-pong | [-25.000000, 41.666667] | [-66.666667, 50.000000] | [-25.000000, 25.000000] | fails |
| Ping-pong minus vehicle | [-25.000000, 25.000000] | [-50.000000, 25.000000] | [0.000000, 25.000000] | fails |
| Tandem minus vehicle | [-25.000000, 41.666667] | [-66.666667, 25.000000] | [0.000000, 25.000000] | fails |

| Contrast | Composition p range | Composition Holm range | T p range | T Holm range |
| --- | --- | --- | --- | --- |
| MTD minus ping-pong (primary) | [0.000054, 1.000000] | [0.000646, 1.000000] | [0.000054, 1.000000] | [0.000646, 1.000000] |
| MTD minus tandem | [0.000054, 1.000000] | [0.000646, 1.000000] | [0.000054, 1.000000] | [0.000646, 1.000000] |
| MTD minus vehicle | [0.004762, 1.000000] | [0.038095, 1.000000] | [0.004762, 1.000000] | [0.038095, 1.000000] |
| Tandem minus ping-pong | [0.001346, 1.000000] | [0.013461, 1.000000] | [0.037267, 1.000000] | [0.260870, 1.000000] |
| Ping-pong minus vehicle | [0.250000, 1.000000] | [1.000000, 1.000000] | [0.250000, 1.000000] | [1.000000, 1.000000] |
| Tandem minus vehicle | [0.081319, 1.000000] | [0.325275, 1.000000] | [0.244505, 1.000000] | [0.733516, 1.000000] |

Envelope of simultaneous difference intervals (minimum lower bound to maximum upper bound):

| Contrast | T (pp) | C (pp) | O (pp) |
| --- | --- | --- | --- |
| MTD minus ping-pong (primary) | [-69.196776, 100.000000] | [-100.000000, 87.140979] | [-69.196776, 64.262311] |
| MTD minus tandem | [-81.980940, 100.000000] | [-100.000000, 95.095557] | [-69.196776, 64.262311] |
| MTD minus vehicle | [-91.781816, 100.000000] | [-100.000000, 91.781816] | [-78.635650, 64.262311] |
| Tandem minus ping-pong | [-69.196776, 81.980940] | [-95.095557, 87.140979] | [-69.196776, 69.196776] |
| Ping-pong minus vehicle | [-91.781816, 69.196776] | [-87.140979, 91.781816] | [-78.635650, 69.196776] |
| Tandem minus vehicle | [-91.781816, 81.980940] | [-95.095557, 91.781816] | [-78.635650, 69.196776] |

#### Pair-local fragility

Other rows remain at their printed labels. Breakpoint zero means the printed comparison already lacks rejection or strict direction. Fractions below use the pair's animals and eligible labels, respectively. Witness row order is vehicle / ping-pong / tandem / combination MTD, each T,C,O.

| Contrast | Loss criterion | Relabels | Fraction of pair animals | Fraction of eligible pair labels | Witness table |
| --- | --- | --- | --- | --- | --- |
| MTD minus ping-pong (primary) | raw_composition | 1 | 1/18 (5.555556%) | 1/12 (8.333333%) | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus ping-pong (primary) | raw_toxicity | 2 | 2/18 (11.111111%) | 2/12 (16.666667%) | (1,3,0)/(3,6,3)/(5,4,3)/(4,2,0) |
| MTD minus ping-pong (primary) | holm_composition | 0 | 0/18 (0.000000%) | 0/12 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | holm_toxicity | 0 | 0/18 (0.000000%) | 0/12 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | direction | 5 | 5/18 (27.777778%) | 5/12 (41.666667%) | (1,3,0)/(3,6,3)/(5,4,3)/(1,5,0) |
| MTD minus ping-pong (primary) | combined | 0 | 0/18 (0.000000%) | 0/12 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus tandem | raw_composition | 0 | 0/18 (0.000000%) | 0/14 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus tandem | raw_toxicity | 1 | 1/18 (5.555556%) | 1/14 (7.142857%) | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus tandem | holm_composition | 0 | 0/18 (0.000000%) | 0/14 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus tandem | holm_toxicity | 0 | 0/18 (0.000000%) | 0/14 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus tandem | direction | 4 | 4/18 (22.222222%) | 4/14 (28.571429%) | (1,3,0)/(3,6,3)/(5,4,3)/(2,4,0) |
| MTD minus tandem | combined | 0 | 0/18 (0.000000%) | 0/14 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | raw_composition | 1 | 1/10 (10.000000%) | 1/7 (14.285714%) | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus vehicle | raw_toxicity | 1 | 1/10 (10.000000%) | 1/7 (14.285714%) | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus vehicle | holm_composition | 0 | 0/10 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | holm_toxicity | 0 | 0/10 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | direction | 5 | 5/10 (50.000000%) | 5/7 (71.428571%) | (1,3,0)/(3,6,3)/(5,4,3)/(1,5,0) |
| MTD minus vehicle | combined | 0 | 0/10 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | raw_composition | 0 | 0/24 (0.000000%) | 0/14 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | raw_toxicity | 0 | 0/24 (0.000000%) | 0/14 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | holm_composition | 0 | 0/24 (0.000000%) | 0/14 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | holm_toxicity | 0 | 0/24 (0.000000%) | 0/14 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | direction | 2 | 2/24 (8.333333%) | 2/14 (14.285714%) | (1,3,0)/(3,6,3)/(3,6,3)/(6,0,0) |
| Tandem minus ping-pong | combined | 0 | 0/24 (0.000000%) | 0/14 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | raw_composition | 0 | 0/16 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | raw_toxicity | 0 | 0/16 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | holm_composition | 0 | 0/16 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | holm_toxicity | 0 | 0/16 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | direction | 0 | 0/16 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | combined | 0 | 0/16 (0.000000%) | 0/7 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | raw_composition | 0 | 0/16 (0.000000%) | 0/9 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | raw_toxicity | 0 | 0/16 (0.000000%) | 0/9 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | holm_composition | 0 | 0/16 (0.000000%) | 0/9 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | holm_toxicity | 0 | 0/16 (0.000000%) | 0/9 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | direction | 2 | 2/16 (12.500000%) | 2/9 (22.222222%) | (1,3,0)/(3,6,3)/(3,6,3)/(6,0,0) |
| Tandem minus vehicle | combined | 0 | 0/16 (0.000000%) | 0/9 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |

#### Primary coherent global fragility

Here edits may affect any row, including changes outside the primary pair that alter the Holm ordering. Fractions use the whole S4 table.

| Loss criterion | Relabels | Fraction of all animals | Fraction of all eligible labels | Witness table |
| --- | --- | --- | --- | --- |
| raw_composition | 1 | 1/34 (2.941176%) | 1/21 (4.761905%) | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| raw_toxicity | 2 | 2/34 (5.882353%) | 2/21 (9.523810%) | (1,3,0)/(3,6,3)/(5,4,3)/(4,2,0) |
| holm_composition | 0 | 0/34 (0.000000%) | 0/21 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| holm_toxicity | 0 | 0/34 (0.000000%) | 0/21 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| direction | 5 | 5/34 (14.705882%) | 5/21 (23.809524%) | (1,3,0)/(3,6,3)/(5,4,3)/(1,5,0) |
| combined | 0 | 0/34 (0.000000%) | 0/21 (0.000000%) | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |

#### Primary worst case at every cumulative global budget

Each row permits at most k changed labels; extrema can be attained by different tables and are not a synthetic joint table. Complete ranges and attaining witnesses for every metric and every pair at each budget are available with `--json`.

| k | States | k/all animals (%) | k/eligible (%) | T difference range (pp) | Worst composition p | Worst composition Holm | Worst T p | Worst T Holm |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1 | 0.000000 | 0.000000 | [75.000000, 75.000000] | 0.016484 | 0.181319 | 0.009050 | 0.108597 |
| 1 | 7 | 2.941176 | 4.761905 | [58.333333, 83.333333] | 0.085973 | 0.945701 | 0.042986 | 0.515837 |
| 2 | 27 | 5.882353 | 9.523810 | [41.666667, 91.666667] | 0.229261 | 1.000000 | 0.141403 | 1.000000 |
| 3 | 77 | 8.823529 | 14.285714 | [25.000000, 100.000000] | 0.456044 | 1.000000 | 0.344053 | 1.000000 |
| 4 | 179 | 11.764706 | 19.047619 | [8.333333, 100.000000] | 0.636393 | 1.000000 | 1.000000 | 1.000000 |
| 5 | 357 | 14.705882 | 23.809524 | [-8.333333, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 6 | 632 | 17.647059 | 28.571429 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 7 | 1015 | 20.588235 | 33.333333 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 8 | 1502 | 23.529412 | 38.095238 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 9 | 2072 | 26.470588 | 42.857143 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 10 | 2688 | 29.411765 | 47.619048 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 11 | 3304 | 32.352941 | 52.380952 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 12 | 3874 | 35.294118 | 57.142857 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 13 | 4361 | 38.235294 | 61.904762 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 14 | 4744 | 41.176471 | 66.666667 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 15 | 5019 | 44.117647 | 71.428571 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 16 | 5197 | 47.058824 | 76.190476 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 17 | 5299 | 50.000000 | 80.952381 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 18 | 5349 | 52.941176 | 85.714286 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 19 | 5369 | 55.882353 | 90.476190 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 20 | 5375 | 58.823529 | 95.238095 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 21 | 5376 | 61.764706 | 100.000000 | [-25.000000, 100.000000] | 1.000000 | 1.000000 | 1.000000 | 1.000000 |

<details>
<summary>Attaining witnesses for all T-and-O extrema</summary>

Differences and interval endpoints are pp; p-values are probabilities. Each extremum has its own coherent table. Interval lower/upper extrema are endpoints of reference-model uncertainty, not additional tests.

| Contrast | Metric | Minimum | Minimum witness | Maximum | Maximum witness |
| --- | --- | --- | --- | --- | --- |
| MTD minus ping-pong (primary) | composition_p | 0.000054 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,7,2)/(5,4,3)/(2,4,0) |
| MTD minus ping-pong (primary) | toxicity_p | 0.000054 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(2,4,0) |
| MTD minus ping-pong (primary) | toxicity_difference | -25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | toxicity_interval_lower | -69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | -4.481304 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | toxicity_interval_upper | 62.037637 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | cancer_difference | -100.000000 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) | 50.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus ping-pong (primary) | cancer_interval_lower | -100.000000 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) | -51.403291 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus ping-pong (primary) | cancer_interval_upper | 4.481304 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) | 87.140979 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus ping-pong (primary) | other_difference | -25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 0.000000 | (1,3,0)/(3,9,0)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | other_interval_lower | -69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | -40.218992 | (1,3,0)/(3,9,0)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | other_interval_upper | 62.037637 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 64.262311 | (1,3,0)/(3,9,0)/(5,4,3)/(6,0,0) |
| MTD minus ping-pong (primary) | composition_holm | 0.000646 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(4,2,0) |
| MTD minus ping-pong (primary) | toxicity_holm | 0.000646 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(4,2,0) |
| MTD minus tandem | composition_p | 0.000054 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,5,2)/(3,3,0) |
| MTD minus tandem | toxicity_p | 0.000054 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(3,3,0) |
| MTD minus tandem | toxicity_difference | -41.666667 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| MTD minus tandem | toxicity_interval_lower | -81.980940 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | -4.481304 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| MTD minus tandem | toxicity_interval_upper | 55.785322 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) |
| MTD minus tandem | cancer_difference | -100.000000 | (1,3,0)/(3,6,3)/(0,12,0)/(6,0,0) | 66.666667 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus tandem | cancer_interval_lower | -100.000000 | (1,3,0)/(3,6,3)/(0,12,0)/(6,0,0) | -40.287484 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus tandem | cancer_interval_upper | 4.481304 | (1,3,0)/(3,6,3)/(0,12,0)/(6,0,0) | 95.095557 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus tandem | other_difference | -25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 0.000000 | (1,3,0)/(3,6,3)/(5,7,0)/(6,0,0) |
| MTD minus tandem | other_interval_lower | -69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | -40.218992 | (1,3,0)/(3,6,3)/(5,7,0)/(6,0,0) |
| MTD minus tandem | other_interval_upper | 62.037637 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 64.262311 | (1,3,0)/(3,6,3)/(5,7,0)/(6,0,0) |
| MTD minus tandem | composition_holm | 0.000646 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus tandem | toxicity_holm | 0.000646 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus vehicle | composition_p | 0.004762 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(2,4,0) |
| MTD minus vehicle | toxicity_p | 0.004762 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(2,4,0) |
| MTD minus vehicle | toxicity_difference | -25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | toxicity_interval_lower | -91.781816 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | -42.897961 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | toxicity_interval_upper | 64.210187 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) | 100.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | cancer_difference | -100.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus vehicle | cancer_interval_lower | -100.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | -64.210187 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus vehicle | cancer_interval_upper | 42.897961 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 91.781816 | (1,3,0)/(3,6,3)/(5,4,3)/(0,6,0) |
| MTD minus vehicle | other_difference | 0.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 0.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | other_interval_lower | -78.635650 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | -78.635650 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | other_interval_upper | 64.262311 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 64.262311 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| MTD minus vehicle | composition_holm | 0.038095 | (0,4,0)/(2,7,3)/(2,7,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| MTD minus vehicle | toxicity_holm | 0.038095 | (0,4,0)/(2,7,3)/(2,7,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(5,1,0) |
| Tandem minus ping-pong | composition_p | 0.001346 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(4,5,3)/(6,0,0) |
| Tandem minus ping-pong | toxicity_p | 0.037267 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(4,5,3)/(6,0,0) |
| Tandem minus ping-pong | toxicity_difference | -25.000000 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 41.666667 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | toxicity_interval_lower | -69.196776 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | -31.742003 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | toxicity_interval_upper | 37.994318 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 81.980940 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | cancer_difference | -66.666667 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) | 50.000000 | (1,3,0)/(3,6,3)/(0,12,0)/(6,0,0) |
| Tandem minus ping-pong | cancer_interval_lower | -95.095557 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) | -27.359972 | (1,3,0)/(3,6,3)/(0,12,0)/(6,0,0) |
| Tandem minus ping-pong | cancer_interval_upper | 16.244165 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) | 87.140979 | (1,3,0)/(3,6,3)/(0,12,0)/(6,0,0) |
| Tandem minus ping-pong | other_difference | -25.000000 | (1,3,0)/(3,6,3)/(5,7,0)/(6,0,0) | 25.000000 | (1,3,0)/(3,9,0)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | other_interval_lower | -69.196776 | (1,3,0)/(3,6,3)/(5,7,0)/(6,0,0) | -37.994318 | (1,3,0)/(3,9,0)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | other_interval_upper | 37.994318 | (1,3,0)/(3,6,3)/(5,7,0)/(6,0,0) | 69.196776 | (1,3,0)/(3,9,0)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | composition_holm | 0.013461 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus ping-pong | toxicity_holm | 0.260870 | (0,4,0)/(0,9,3)/(5,7,0)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | composition_p | 0.250000 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,7,2)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | toxicity_p | 0.250000 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | toxicity_difference | -25.000000 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 25.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | toxicity_interval_lower | -91.781816 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | -76.410975 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | toxicity_interval_upper | 40.166868 | (1,3,0)/(0,9,3)/(5,4,3)/(6,0,0) | 69.196776 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | cancer_difference | -50.000000 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 25.000000 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | cancer_interval_lower | -87.140979 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | -40.166868 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | cancer_interval_upper | 65.776629 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 91.781816 | (1,3,0)/(0,12,0)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | other_difference | 0.000000 | (1,3,0)/(3,9,0)/(5,4,3)/(6,0,0) | 25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | other_interval_lower | -78.635650 | (1,3,0)/(3,9,0)/(5,4,3)/(6,0,0) | -76.410975 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | other_interval_upper | 40.218992 | (1,3,0)/(3,9,0)/(5,4,3)/(6,0,0) | 69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | composition_holm | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Ping-pong minus vehicle | toxicity_holm | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | composition_p | 0.081319 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,7,0)/(6,0,0) |
| Tandem minus vehicle | toxicity_p | 0.244505 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | toxicity_difference | -25.000000 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 41.666667 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | toxicity_interval_lower | -91.781816 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | -70.158660 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | toxicity_interval_upper | 40.166868 | (1,3,0)/(3,6,3)/(0,9,3)/(6,0,0) | 81.980940 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | cancer_difference | -66.666667 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 25.000000 | (1,3,0)/(3,6,3)/(0,12,0)/(6,0,0) |
| Tandem minus vehicle | cancer_interval_lower | -95.095557 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | -40.166868 | (1,3,0)/(3,6,3)/(0,12,0)/(6,0,0) |
| Tandem minus vehicle | cancer_interval_upper | 54.660823 | (0,4,0)/(3,6,3)/(5,4,3)/(6,0,0) | 91.781816 | (1,3,0)/(3,6,3)/(0,12,0)/(6,0,0) |
| Tandem minus vehicle | other_difference | 0.000000 | (1,3,0)/(3,6,3)/(5,7,0)/(6,0,0) | 25.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | other_interval_lower | -78.635650 | (1,3,0)/(3,6,3)/(5,7,0)/(6,0,0) | -76.410975 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | other_interval_upper | 40.218992 | (1,3,0)/(3,6,3)/(5,7,0)/(6,0,0) | 69.196776 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | composition_holm | 0.325275 | (0,4,0)/(0,9,3)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |
| Tandem minus vehicle | toxicity_holm | 0.733516 | (0,4,0)/(1,9,2)/(5,4,3)/(6,0,0) | 1.000000 | (1,3,0)/(3,6,3)/(5,4,3)/(6,0,0) |

</details>

### Single-drug denominator sensitivity

Both unconditional single-drug verdicts are **untestable**. The following bounds assume exactly one excess adaptive category assignment; they are not source corrections. Differences here are adaptive minus its own MTD, unlike the primary MTD-minus-adaptive direction. Bonferroni values use the separate four-test conditional family and cannot rescue the S4 verdict.

#### CAP

All possible one-cell-decrement completions: (3,5,4), (4,4,4), (4,5,3).

| Model | Completion/state scenarios | Unique tables | T difference (pp) | C difference (pp) | O difference (pp) | Composition p | T p | Composition Bonferroni | T Bonferroni |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| as_printed_completion | 3 | 3 | [-8.333333, 0.000000] | [-33.333333, -25.000000] | [25.000000, 33.333333] | [0.303167, 0.557854] | [1.000000, 1.000000] | [1.000000, 1.000000] | [1.000000, 1.000000] |
| toxicity_only | 42 | 30 | [-33.333333, 33.333333] | [-66.666667, 8.333333] | [25.000000, 33.333333] | [0.021978, 0.760019] | [0.098039, 1.000000] | [0.087912, 1.000000] | [0.392157, 1.000000] |
| toxicity_and_other | 195 | 75 | [-33.333333, 33.333333] | [-66.666667, 33.333333] | [0.000000, 33.333333] | [0.021978, 1.000000] | [0.098039, 1.000000] | [0.087912, 1.000000] | [0.392157, 1.000000] |

<details>
<summary>Conditional single-drug attaining witnesses</summary>

| Model | Metric | Minimum | Minimum witness | Maximum | Maximum witness |
| --- | --- | --- | --- | --- | --- |
| as_printed_completion | toxicity_difference | -8.333333 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) | 0.000000 | decrement completion (4,4,4); A/M=(4,4,4)/(2,4,0) |
| as_printed_completion | cancer_difference | -33.333333 | decrement completion (4,4,4); A/M=(4,4,4)/(2,4,0) | -25.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| as_printed_completion | other_difference | 25.000000 | decrement completion (4,5,3); A/M=(4,5,3)/(2,4,0) | 33.333333 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| as_printed_completion | composition_p | 0.303167 | decrement completion (4,4,4); A/M=(4,4,4)/(2,4,0) | 0.557854 | decrement completion (4,5,3); A/M=(4,5,3)/(2,4,0) |
| as_printed_completion | toxicity_p | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| as_printed_completion | composition_bonferroni | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| as_printed_completion | toxicity_bonferroni | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| toxicity_only | toxicity_difference | -33.333333 | decrement completion (3,5,4); A/M=(0,8,4)/(2,4,0) | 33.333333 | decrement completion (4,4,4); A/M=(4,4,4)/(0,6,0) |
| toxicity_only | cancer_difference | -66.666667 | decrement completion (4,4,4); A/M=(4,4,4)/(0,6,0) | 8.333333 | decrement completion (4,5,3); A/M=(0,9,3)/(2,4,0) |
| toxicity_only | other_difference | 25.000000 | decrement completion (4,5,3); A/M=(4,5,3)/(2,4,0) | 33.333333 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| toxicity_only | composition_p | 0.021978 | decrement completion (4,4,4); A/M=(4,4,4)/(0,6,0) | 0.760019 | decrement completion (4,5,3); A/M=(2,7,3)/(1,5,0) |
| toxicity_only | toxicity_p | 0.098039 | decrement completion (3,5,4); A/M=(0,8,4)/(2,4,0) | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| toxicity_only | composition_bonferroni | 0.087912 | decrement completion (4,4,4); A/M=(4,4,4)/(0,6,0) | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| toxicity_only | toxicity_bonferroni | 0.392157 | decrement completion (3,5,4); A/M=(0,8,4)/(2,4,0) | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| toxicity_and_other | toxicity_difference | -33.333333 | decrement completion (3,5,4); A/M=(0,8,4)/(2,4,0) | 33.333333 | decrement completion (4,4,4); A/M=(4,4,4)/(0,6,0) |
| toxicity_and_other | cancer_difference | -66.666667 | decrement completion (4,4,4); A/M=(4,4,4)/(0,6,0) | 33.333333 | decrement completion (3,5,4); A/M=(0,12,0)/(2,4,0) |
| toxicity_and_other | other_difference | 0.000000 | decrement completion (3,5,4); A/M=(3,9,0)/(2,4,0) | 33.333333 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| toxicity_and_other | composition_p | 0.021978 | decrement completion (4,4,4); A/M=(4,4,4)/(0,6,0) | 1.000000 | decrement completion (3,5,4); A/M=(3,7,2)/(2,4,0) |
| toxicity_and_other | toxicity_p | 0.098039 | decrement completion (3,5,4); A/M=(0,8,4)/(2,4,0) | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| toxicity_and_other | composition_bonferroni | 0.087912 | decrement completion (4,4,4); A/M=(4,4,4)/(0,6,0) | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |
| toxicity_and_other | toxicity_bonferroni | 0.392157 | decrement completion (3,5,4); A/M=(0,8,4)/(2,4,0) | 1.000000 | decrement completion (3,5,4); A/M=(3,5,4)/(2,4,0) |

</details>

#### GEM

All possible one-cell-decrement completions: (6,1,5), (7,0,5), (7,1,4).

| Model | Completion/state scenarios | Unique tables | T difference (pp) | C difference (pp) | O difference (pp) | Composition p | T p | Composition Bonferroni | T Bonferroni |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| as_printed_completion | 3 | 3 | [-33.333333, -25.000000] | [-16.666667, -8.333333] | [33.333333, 41.666667] | [0.099925, 0.313833] | [0.315611, 0.600032] | [0.399698, 1.000000] | [1.000000, 1.000000] |
| toxicity_only | 138 | 96 | [-83.333333, 58.333333] | [-100.000000, 50.000000] | [33.333333, 41.666667] | [0.000054, 0.421461] | [0.000700, 1.000000] | [0.000215, 1.000000] | [0.002801, 1.000000] |
| toxicity_and_other | 780 | 288 | [-83.333333, 58.333333] | [-100.000000, 83.333333] | [0.000000, 41.666667] | [0.000054, 1.000000] | [0.000700, 1.000000] | [0.000215, 1.000000] | [0.002801, 1.000000] |

<details>
<summary>Conditional single-drug attaining witnesses</summary>

| Model | Metric | Minimum | Minimum witness | Maximum | Maximum witness |
| --- | --- | --- | --- | --- | --- |
| as_printed_completion | toxicity_difference | -33.333333 | decrement completion (6,1,5); A/M=(6,1,5)/(5,1,0) | -25.000000 | decrement completion (7,0,5); A/M=(7,0,5)/(5,1,0) |
| as_printed_completion | cancer_difference | -16.666667 | decrement completion (7,0,5); A/M=(7,0,5)/(5,1,0) | -8.333333 | decrement completion (6,1,5); A/M=(6,1,5)/(5,1,0) |
| as_printed_completion | other_difference | 33.333333 | decrement completion (7,1,4); A/M=(7,1,4)/(5,1,0) | 41.666667 | decrement completion (6,1,5); A/M=(6,1,5)/(5,1,0) |
| as_printed_completion | composition_p | 0.099925 | decrement completion (7,0,5); A/M=(7,0,5)/(5,1,0) | 0.313833 | decrement completion (7,1,4); A/M=(7,1,4)/(5,1,0) |
| as_printed_completion | toxicity_p | 0.315611 | decrement completion (6,1,5); A/M=(6,1,5)/(5,1,0) | 0.600032 | decrement completion (7,0,5); A/M=(7,0,5)/(5,1,0) |
| as_printed_completion | composition_bonferroni | 0.399698 | decrement completion (7,0,5); A/M=(7,0,5)/(5,1,0) | 1.000000 | decrement completion (7,1,4); A/M=(7,1,4)/(5,1,0) |
| as_printed_completion | toxicity_bonferroni | 1.000000 | decrement completion (6,1,5); A/M=(6,1,5)/(5,1,0) | 1.000000 | decrement completion (6,1,5); A/M=(6,1,5)/(5,1,0) |
| toxicity_only | toxicity_difference | -83.333333 | decrement completion (6,1,5); A/M=(0,7,5)/(5,1,0) | 58.333333 | decrement completion (7,0,5); A/M=(7,0,5)/(0,6,0) |
| toxicity_only | cancer_difference | -100.000000 | decrement completion (7,0,5); A/M=(7,0,5)/(0,6,0) | 50.000000 | decrement completion (7,1,4); A/M=(0,8,4)/(5,1,0) |
| toxicity_only | other_difference | 33.333333 | decrement completion (7,1,4); A/M=(7,1,4)/(5,1,0) | 41.666667 | decrement completion (6,1,5); A/M=(6,1,5)/(5,1,0) |
| toxicity_only | composition_p | 0.000054 | decrement completion (7,0,5); A/M=(7,0,5)/(0,6,0) | 0.421461 | decrement completion (7,1,4); A/M=(6,2,4)/(4,2,0) |
| toxicity_only | toxicity_p | 0.000700 | decrement completion (6,1,5); A/M=(0,7,5)/(5,1,0) | 1.000000 | decrement completion (6,1,5); A/M=(6,1,5)/(3,3,0) |
| toxicity_only | composition_bonferroni | 0.000215 | decrement completion (7,0,5); A/M=(7,0,5)/(0,6,0) | 1.000000 | decrement completion (7,1,4); A/M=(7,1,4)/(5,1,0) |
| toxicity_only | toxicity_bonferroni | 0.002801 | decrement completion (6,1,5); A/M=(0,7,5)/(5,1,0) | 1.000000 | decrement completion (6,1,5); A/M=(6,1,5)/(5,1,0) |
| toxicity_and_other | toxicity_difference | -83.333333 | decrement completion (6,1,5); A/M=(0,7,5)/(5,1,0) | 58.333333 | decrement completion (7,0,5); A/M=(7,0,5)/(0,6,0) |
| toxicity_and_other | cancer_difference | -100.000000 | decrement completion (7,0,5); A/M=(7,0,5)/(0,6,0) | 83.333333 | decrement completion (6,1,5); A/M=(0,12,0)/(5,1,0) |
| toxicity_and_other | other_difference | 0.000000 | decrement completion (6,1,5); A/M=(6,6,0)/(5,1,0) | 41.666667 | decrement completion (6,1,5); A/M=(6,1,5)/(5,1,0) |
| toxicity_and_other | composition_p | 0.000054 | decrement completion (7,0,5); A/M=(7,0,5)/(0,6,0) | 1.000000 | decrement completion (6,1,5); A/M=(6,5,1)/(4,2,0) |
| toxicity_and_other | toxicity_p | 0.000700 | decrement completion (6,1,5); A/M=(0,7,5)/(5,1,0) | 1.000000 | decrement completion (6,1,5); A/M=(6,1,5)/(3,3,0) |
| toxicity_and_other | composition_bonferroni | 0.000215 | decrement completion (7,0,5); A/M=(7,0,5)/(0,6,0) | 1.000000 | decrement completion (6,1,5); A/M=(6,2,4)/(5,1,0) |
| toxicity_and_other | toxicity_bonferroni | 0.002801 | decrement completion (6,1,5); A/M=(0,7,5)/(5,1,0) | 1.000000 | decrement completion (6,1,5); A/M=(6,1,5)/(5,1,0) |

</details>

### Printed dose context and source disagreements

S5 means are not used in any cause test, follow-up reconstruction, or dose-response adjustment. Daily units follow the column heading; the printed cells say mg/kg.

| S5 source row | Drug | Mean cumulative mg/kg | Mean daily mg/kg/day |
| --- | --- | --- | --- |
| Gemcitabine Dose Modulation | GEM | 256 | 5.8 |
| Gemcitabine Intermittent | GEM | 33.3 | 5.5 |
| Gemcitabine Maximum Tolerated Dose | GEM | 783.3 | 12.8 |
| Capecitabine Dose Modulation | CAP | 533.8 | 11.8 |
| Capecitabine Intermittent | CAP | 1226.6 | 26 |
| Capecitabine Maximum Tolerated Dose | CAP | 1373.3 | 28.4 |
| Capecitabine + Gemcitabine Ping-Pong Dose Modulation | GEM | 379.1 | 4.9 |
| Capecitabine + Gemcitabine Ping-Pong Dose Modulation | CAP | 475.5 | 6.1 |
| Capecitabine + Gemcitabine in Tandem Dose Modulation | GEM | 229 | 4.3 |
| Capecitabine + Gemcitabine in Tandem Dose Modulation | CAP | 492.1 | 8.7 |
| Capecitabine + Gemcitabine Ping-Pong Intermittent | GEM | 550 | 6.4 |
| Capecitabine + Gemcitabine Ping-Pong Intermittent | CAP | 660 | 8 |
| Capecitabine + Gemcitabine in Tandem Intermittent | GEM | 341.6 | 9 |
| Capecitabine + Gemcitabine in Tandem Intermittent | CAP | 606.6 | 12.5 |
| Capecitabine + Gemcitabine Maximum Tolerated Dose | GEM | 433.3 | 20.6 |
| Capecitabine + Gemcitabine Maximum Tolerated Dose | CAP | 766.6 | 30.6 |

| Arm | Drug | S5 cumulative mean | Article cumulative mean | Disposition |
| --- | --- | --- | --- | --- |
| gem_dm | GEM | 256 | 256.9 | Unresolved; both retained |
| gem_intermittent | GEM | 33.3 | 333.3 | Unresolved; both retained |

### Reproducibility identifiers

Original local protocol-only freeze commit: `dbf3c54746ba853a81a65f177108531081460978` (publication mapping above). Text hashes below normalize CRLF and CR to LF and preserve all other bytes and final-newline state. Paths use `Path.as_posix()`. These hashes attest content consistency, not source truth.

| Relative path | SHA-256 of newline-normalised content |
| --- | --- |
| `PROTOCOL.md` | `0f228ecb01a7d187b37a57f43d163fb8b4cedb954ee0258c24ccd94363283b90` |
| `analysis.py` | `429decc3d6073b1eace94a9f27a1c45aea6349e4f6a56b8c739d005ba2519501` |
| `data/allocations.csv` | `fb9f938da5d3e61672f6ff8fccc2cab8c278530aff7a5449105c880d97eb9ac4` |
| `data/counts.csv` | `f79fc1a13ebf7b8c36cf480491822c5d523c40b9488f2acf260818e74a97940b` |
| `data/dose_conflicts.csv` | `c9af39a48d9254fc6c63debf8cf48963024dae84aa37b759f4e72c2e9a5ac3c5` |
| `data/doses.csv` | `826de1f5b7c456745922a69644339006239e12b8a24b0c8db96d85cf53686093` |

<!-- END GENERATED RESULTS -->

## Falsification attempts and hostile review

The following checks were chosen to undermine an overly favourable reading.
They are not independent scientific review or validation of the animal model.

| Attempt | What the audit establishes or preserves |
| --- | --- |
| Treat all terminal counts as reconciled animal counts | Fails for the two single-drug adaptive rows. The printed and allocated totals remain separate; conditional bounds never replace them. |
| Interpret the toxicity column as proven drug-caused death | The untreated group and the authors' own unresolved cachexia discussion defeat that interpretation. |
| Treat an unadjusted primary p-value as the declared discovery | The frozen full-family correction is retained; the printed-label result fails that rule. |
| Erase the primary result with the fewest possible label changes | Exhaustive pair-local and coherent global searches supply minimum breakpoints and concrete witnesses, including initially nonrejecting cases. |
| Preserve an apparent benefit by changing the shared control differently for different contrasts | Prohibited by the single four-row attribution state; the same control assignment enters every pair. |
| Obtain extra sample size by splitting pools or duplicating controls | The source pools and allocation joins are fixed; neither split data nor duplicated controls enter the analysis. |
| Promote a conditional single-drug completion to corrected source data | Every possible decrement is retained under a clearly hypothetical error model; the unconditional verdict remains untestable. |
| Resolve source dose disagreements by picking the more plausible number | Both statements are shown, and neither enters the terminal-cause inference. |
| Infer equal follow-up from cumulative and daily dose means | Rejected: a ratio of means is not a mean duration, dose means lack dispersion, and treatment duration depends on outcome. |
| Use a checksum or passing tests to establish biological truth | Rejected: these check source-file identity and computational consistency only. |

An independent source audit checked every transcribed numeric CSV entry
against source text and rendered pages. The analysis and test use separate
exact-probability implementations, check extremizing states and minimum
breakpoints, and exercise normalization of Windows and Linux line endings.
Any execution limitations or deviations from the frozen protocol are recorded
here rather than changing that protocol.

Hostile review corrected the presentation of fragility fractions to preserve
their actual animal and eligible-label denominators, including when no labels
change. A repository-wide test run that overlapped that display edit detected
the difference between its already-loaded renderer and the updated README.
The dedicated suites passed on the final files; the repository-wide suite was
then restarted on the fixed revision. This was generated-output drift during
editing, not a change to the counts, exact tests, bounds or verdict.

**Protocol deviations:** none in analysis specification. Source discrepancies
were known and declared before the freeze. The publication transport recreated
Git commit objects as documented above; it did not change the frozen protocol
or numerical analysis. Execution and validation results are reported with the
pull request.

## Interpretation limits

**Unequal observation times and competing risks prevent interpretation as
equal-horizon cancer mortality, cause-specific survival or hazards, competing-
risk cumulative incidence, or causal attribution.** A smaller terminal
fraction assigned to cancer can occur because another event ends observation
earlier. A larger "other" fraction can reflect later observation, technical
events or euthanasia rules. It is not demonstrated cancer control.

The source's section 2.3.3 censors some locomotion-related and implantation-
surgery events in its survival analysis, subject to the specified burden
criteria. Its terminal-category table therefore cannot be substituted for its
survival event/censoring data. No animal-level times, reconstructed survival
curves, equal-risk windows or competing-risk estimators are available in this
aggregate analysis.

The observed categories do not isolate drug-mediated death from cancer,
frailty, humane endpoints or other causes. Drug combination, dose intensity
and adaptive scheduling differ together. The experiment uses a single
endocrine-resistant breast-cancer cell-line xenograft setting in
immunodeficient mice; it does not identify a general schedule effect.

Neither competition between sensitive and resistant clones nor division-
linked mutation supply is measured by a terminal label. PR #61's theorem
requires a defined population process, establishment-weighted mutation supply,
common finite course, and an uncensored progression-time estimand. These
cause compositions provide none of the missing quantities needed to locate a
real tumour relative to that theorem's boundary. No treatment inference follows.

## Honest novelty assessment

Seyedi et al. already display cause patterns in Figure S7 and Tables S2-S4,
interpret them on supplement p.32, and explicitly discuss the inability to
separate toxicity from cancer. Their article section 4.1 already warns that
the constant combination comparator may exceed tolerable dosing and may
explain much of the multidrug adaptive advantage. Their sections 2.5 and 3.3
analyze survival; section 3.4 relates daily average dose to survival.

This contribution is a **limited repository-specific audit**: an immutable
computation protocol, retained denominator contradictions, exact small-sample
cause-composition comparisons, multiplicity-aware uncertainty, and quantified
attribution fragility using the authors' existing aggregates. No formal exact
cause-composition/reclassification analysis was located in the inspected
article or supplement. That observation supports a narrow extension of the
reported analysis; it does not establish worldwide priority, a new biological
discovery, independent experimental confirmation or treatment efficacy.
