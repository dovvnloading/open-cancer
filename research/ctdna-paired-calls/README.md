# Paired ctDNA calls: a sharp compatibility audit

**Nonclinical, unreviewed research. No treatment, diagnosis, surveillance or
patient-selection advice.** This is an analysis of published aggregate calls,
not an assay validation, biological discovery or clinical decision rule.

This runs [ranked job 3](../README.md) using the DYNAMIC postoperative paired
call categories. The question is how much combined latitude in observation
errors and genuine threshold reversals a persistent-threshold explanation
needs. The answer is a sharp bound on a **sum of unconditional probabilities**,
not an estimate of the assay's error rate. The separate comparable-kernel
restriction is retained even when it is not rejected.

**Outcome:** the zero-latitude special case and the frozen tight total cap
are **falsified under stated caps**, including the adverse missingness
completion. Broader allowances are **compatible with stated latitude**.
The separate noisy-observation restriction is not rejected under any allowed
missingness completion. Actual assay error, biological clearance and the
unconditional no-intervention interpretation remain **untestable**. The
[recomputed results](#recomputed-results) give the exact frontier and each
conditional verdict.

## Protocol, source and availability audit

[PROTOCOL.md](PROTOCOL.md) was frozen locally in commit
`74577826ea6a84fc96f67c9aeafda0e0b7f34838`, whose only added file is that
protocol, at **2026-11-18 04:26:07 America/New_York**. It precedes CSV creation
and numerical computation. Its published protocol-only counterpart is
[`7423b959`](https://github.com/dovvnloading/open-cancer/commit/7423b959f8162310569a733a956b2cfebe948f31).
The branch starts at `main` revision
`620c1df128ce3a30dbaaa190179cb5404048174a`. The source counts were already
disclosed in the task and research index and were inspected before freezing:
this is a pre-computation protocol, **not blinded testing, independent
replication or registered preregistration**. There are no changes to the
frozen model, caps, error allocation or decision rules.

**Publication transport:** shell `git push` lacked credentials. The authorized
GitHub app therefore published a protocol-only commit followed by the
completed analysis. Both the original local freeze and its remote counterpart
have exactly the same Git tree,
`5abb9ca85a74368ad1e158dd19f1b1ed7f8aaea8`, on the same base parent. The
protocol's normalized content hash below is unchanged. The remote commit's
later publication time is **not** claimed to be the pre-computation freeze
time; the original local SHA is recorded here, not presented as a remotely
reachable commit. This delivery note changes no scientific protocol.

The sole empirical source is [Tie, Wang, Lo et al. (2025), DYNAMIC five-year
report](https://pmc.ncbi.nlm.nih.gov/articles/PMC12974608/),
[DOI 10.1038/s41591-025-03579-w](https://doi.org/10.1038/s41591-025-03579-w).
The accessed version is the PMC author manuscript, publicly available there
from 2026-05-01. HTML access and the rights check occurred on **2026-11-18
(America/New_York)**. No source document or source dataset is redistributed.

The small [aggregate CSV](data/aggregates.csv) stores printed counts only.
Every row has a source URL, DOI, section/table locator, access date, population,
unit, version, article rights and reuse basis. Its population field identifies
arm membership; the successful-analysis and full-arm denominators are
separate estimands. The unsuccessful-analysis count is derived, not passed
off as a separately printed datum.

| Source locator | What was checked | Consequence |
| --- | --- | --- |
| Results, **OS, DSS and RFS in the ctDNA-guided population**, first paragraph | Guided-arm and successful-analysis denominators, each explicitly printed positive-call category, and the at-least-one-positive total | The positive-pattern partition is directly observed in print; a positive/negative pair is an observed call reversal. |
| Extended Data Table 1, ctDNA-guided rows, **No. of Patients** column | Printed negative and positive classification counts, cross-checked with the Results total | The CSV retains the printed negative classification. It does not pretend the source directly prints a negative/negative cell. |
| Methods, **Study design and participants**, paragraph beginning with allocation to ctDNA-guided management | Concurrent analysis of both postoperative samples and the definition requiring both calls to be negative | Together with the Results collection statement and the table partition, this supports the **reported paired classification**. The negative/negative cell is a documented inference. |
| Figure 1 caption, reasons for ITT exclusions | Missing later draws, ineligibility and withdrawal before ITT analysis | These all-arm exclusions are outside the guided ITT denominator. They are not added to the unsuccessful analyses. |

**Availability qualification:** collection at both times alone would be
insufficient. The combined reporting cross-check supports using the negative
classification as negative at both occasions. No separate sentence certifying
the successful count as individually audited complete pairs, per-occasion QC
table or individual assay-availability audit is provided. All full-table
calculations below are conditional on the article's reported classification
being complete as defined. `compute_results(...,
paired_classification_supported=False)` refuses full-table inference. If that
reporting interpretation fails, the full-table conclusion is **untestable**;
hidden discordants cannot be assumed absent.

### Rights and retrieval limitations

The [publisher's Rights and permissions
statement](https://www.nature.com/articles/s41591-025-03579-w#rightslink)
retains exclusive rights with Springer Nature or its licensor. **No Creative
Commons licence was established.** PMC availability does not grant a right
to mirror the article. The reuse basis here is the isolated factual counts,
in our own schema and with attribution, under the facts/expression distinction
described by the [U.S. Copyright Office](https://www.copyright.gov/help/faq/faq-protect.html).
This is not a claim of publisher permission, an open article licence or MIT
relicensing of source content. Article prose, expressive table arrangement,
figures and source databases are not copied into the repository.

The initial PMC route intermittently returned a browser-check page; the
article HTML was subsequently readable through its canonical indexed entry.
The direct Nature request ended at an inaccessible identity redirect, so the
rights statement was checked in indexed **official publisher HTML**. No
browser was used, no browser download prompt appeared, and no article PDF,
supplement, source-data file or image was downloaded. A broad text retrieval
incidentally included neighboring patient-table content; none was used,
transcribed, reconstructed or retained. The source audit and CSV use only
the declared aggregate sections. Method-history references below were
checked as publisher metadata, not used as additional empirical data.

## What is being tested

Let `Z4,Z7` be latent binary threshold status and `Y4,Y7` the observed calls.
Within the stated population, define

\[
e_4=P(Y_4\ne Z_4),\qquad e_7=P(Y_7\ne Z_7),\qquad
b=P(Z_4=1,Z_7=0),\qquad q=P(Y_4=1,Y_7=0).
\]

Persistence is `b=0`; positive `b` explicitly relaxes persistence. All rates
are unconditional within the relevant population. They are not conditional
false-positive/false-negative rates, nor fractions conditional on being
positive at the earlier draw. `e4+e7` is the expected number of mismatching
occasions per participant; with overlapping errors it is not the probability
that a participant has any error. Neither occasion has a positive separate
lower bound without an allocation assumption.

An observed reversal requires an error on at least one occasion or a genuine
latent reversal. Therefore

\[
q\le e_4+e_7+b.
\]

For upper caps `E4,E7,B`, the **sharp compatibility frontier** is

\[
\boxed{E_4+E_7+B\ge q},\qquad
\min(e_4+e_7)=\max(0,q-B).
\]

This is sufficient as well as necessary for matching the entire observed
table in the deliberately broad model. Split its positive/negative mass
between three explanations: latent negative/negative with an earlier error,
latent positive/positive with a later error, and a genuine latent reversal
with no error. Keep latent status equal to observed status on every other
cell. Any allocation within the caps constructs a fitting law. These are
observationally indistinguishable witnesses, not fitted mechanisms. Arbitrary
fixed *exact* error triples need additional constraints; the displayed
frontier concerns **upper caps**. Equal per-occasion caps require half of the
remaining total allowance. Unrestricted genuine reversal can absorb the
entire observed reversal mass without any observation error.

The point frontier exactly fits the empirical probability table. The
uncertainty frontier instead uses a one-sided exact binomial lower bound.
A hypothetical cap can fail exact point fitting and still survive the
sampling test. The frozen primary total cap and sensitivity grid are
assumptions, not measured properties of this assay.

### Sampling and missingness

Clopper-Pearson bounds assume independent identically distributed case
pairs within the specified population. They do not require independence of
the two calls within a case. The binomial model is an inferential assumption,
not a conclusion from treatment randomization or proof of transport to a
general patient population.

The frozen two-hypothesis family assigns each of the cap and directional
tests a one-sided significance level of 0.025. This gives a Bonferroni bound
of 0.05 for that family within a stated population. The extra one-sided 95%
bound and all missingness calculations are labelled sensitivities; the note
does not assert simultaneous interval coverage across different populations
or every displayed interval. The nested cap grid is evaluated using the
same decision lower bound. Equality is retained.

Let `x` be observed reversals, `n` successfully classified cases and `N` the
full guided-arm population. Unknown outcomes for `m=N-n` failures are
allocated exhaustively to the four call patterns. The resulting sharp
finite-panel range is

\[
q_{\rm full}\in[x/N,(x+m)/N].
\]

The lower endpoint is sufficient for some completion; the upper endpoint
is sufficient for every completion. These are identification bounds, not
confidence limits. The sampling lower bound `L(x,N)` is conservative over
all completions, and the outer two-sided interval uses `L(x,N)` and
`U(x+m,N)`. These statements require no missing-at-random assumption, but
their confidence interpretation still requires the stated full-arm iid pair
model. Hypothetical completed calls are a sensitivity analysis, not recovered
outcomes. Error probabilities in the successful subset are not silently
transported to the full arm.

### Separate noisy observation restriction

A common nondecreasing detection kernel, stable relevant assay covariates,
nondecreasing **detectable** burden and no intervention-induced change imply

\[
P(Y_7=1)-P(Y_4=1)=P(-/+)-P(+/-)\ge0.
\]

This restriction **requires nondecreasing marginal detection probability**.
Latent persistence with unrestricted occasion-specific errors does not
imply it. Concurrent processing or a shared assay name alone does not prove
stable shedding, molecular input, background DNA or an unchanged kernel.
The exact one-sided test conditions on paired discordants and tests whether
their positive/negative fraction exceeds one half. No independent-samples
test is substituted.

[PR #58](../ctdna-surveillance-schedule/README.md), Section 2 and Theorems 1
and 4, already allow finite-molecule misses. Thus positive/negative calls are
not a falsification of its full noisy model. This audit's broad error-cap
family also allows false positives, so it is not an estimate of that note's
narrower gated molecular parameters. [PR #63](../prediction-confrontation/README.md)
retains untestable schedule and lead-time comparisons; this different
observation channel does not reverse those verdicts.

<!-- BEGIN RECOMPUTED RESULTS -->
### Recomputed results

**Nonclinical, unreviewed. No treatment, diagnosis or surveillance advice.**

All numerical outputs below are calculated from the attributed aggregate CSV. Intervals assume independent identically distributed case pairs in the stated population; calls within a case may be dependent. Rates are unconditional within that population. The two primary tests each use alpha = 0.025, for Bonferroni familywise error at most 0.05 within the stated target population. The full-arm calculations are missingness robustness sensitivities; no additional simultaneous coverage across both populations or all displayed intervals is claimed. Equality at a rejection boundary is retained.

#### Reported paired classification and sharp frontier

The source cross-check supports complete **reported classification**. The negative/negative cell is inferred from the printed negative category and its two-draw definition; an individual assay-availability audit and a separately printed complete-pairs count remain unavailable.

| Week-4 / week-7 category | Count |
| --- | ---: |
| +/+ | 31 |
| +/− | 8 |
| −/+ | 6 |
| −/−, inferred reported category | 246 |
| Successful-analysis denominator | 291 |

| Quantity | Recomputed value |
| --- | --- |
| Exact plug-in total-latitude frontier q | 8/291 = 2.749141% |
| Decision lower frontier, one-sided 97.5% | 1.194198% |
| Marginal descriptive lower frontier, one-sided 95% | 1.375155% |
| Ordinary two-sided 95% exact interval for q | [1.194198%, 5.344834%] |
| Week-4 positive calls / marginal rate | 39; 13/97 = 13.402062% |
| Week-7 positive calls / marginal rate | 37; 37/291 = 12.714777% |
| Week-7 minus week-4 marginal probability | -2/291 = -0.687285% |

For upper caps E4, E7 and B, exact table compatibility holds precisely when E4 + E7 + B >= q. With B = 0, the combined observation-error allowance must reach the plug-in frontier for an exact empirical-table fit. The uncertainty rule rejects only below the decision lower frontier. The sum e4 + e7 is an expected number of occasion mismatches per participant, not necessarily the probability of any mismatch; neither occasion has an identified individual lower bound.

| Hypothetical total cap | Exact binomial upper tail | Exact plug-in table fits? | Population uncertainty verdict |
| --- | ---: | --- | --- |
| 0.000000% | 0 | no | **falsified under stated caps** |
| 0.500000% | 0.00012953562416 | no | **falsified under stated caps** |
| 1.000000% | 0.00968545735207 | no | **falsified under stated caps** |
| 2.000000% | 0.230032685196 | no | **compatible with stated latitude** |
| 3.000000% | 0.646822006956 | yes | **compatible with stated latitude** |
| 5.000000% | 0.979173079558 | yes | **compatible with stated latitude** |

**Primary total cap 1%: falsified under stated caps.** All caps are hypothetical and uncalibrated. The nested sensitivity decisions use the same lower bound. Exact rational tails are retained in the JSON output; the displayed decimals are rounded.

| Hypothetical reversal cap B | Minimum combined error for exact fit | Minimum equal cap per occasion for exact fit | Decision lower combined-error frontier | Decision lower equal-cap frontier |
| --- | ---: | ---: | ---: | ---: |
| 0.000000% | 2.749141% | 1.374570% | 1.194198% | 0.597099% |
| 0.500000% | 2.249141% | 1.124570% | 0.694198% | 0.347099% |
| 1.000000% | 1.749141% | 0.874570% | 0.194198% | 0.097099% |
| 2.000000% | 0.749141% | 0.374570% | 0.000000% | 0.000000% |
| 3.000000% | 0.000000% | 0.000000% | 0.000000% | 0.000000% |

These are max(0, q − B) and half that value, with the same transformations of the decision lower bound. Allowing genuine reversal changes the persistence assumption; unrestricted reversal need not require any observation error. Constructed compatibility does not identify the mechanism.

#### Unsuccessful analyses: all completions, no missing-at-random assumption

The full guided ITT denominator is 294; subtracting 291 successful analyses leaves 3 unsuccessful analyses. There are 20 allocations to (++,+−,−+,−−). These hypothetical outcomes are sensitivity completions, not recovered calls. The full-arm and successful-analysis error probabilities are different estimands.

| Full-arm quantity | Recomputed value |
| --- | --- |
| Sharp finite-panel q range | 4/147, 11/294; [2.721088%, 3.741497%] |
| Completion-specific one-sided 97.5% lower-frontier range | [1.181938%, 1.882281%] |
| Completion-robust outer two-sided 95% interval | [1.181938%, 6.595835%] |

The lower finite-panel endpoint is enough for at least one completion; the upper endpoint allows every completion. They are sharp identification bounds, not confidence intervals. The lower uncertainty endpoint is valid for every completion under the full-arm iid sampling model even when missingness is arbitrary.

| Total cap | Completions rejected | Exact-tail range, rounded | Rejection scope |
| --- | ---: | --- | --- |
| 0.000000% | 20/20 | [0, 0] | every completion |
| 0.500000% | 20/20 | [3.937443501e-07, 0.000138917337304] | every completion |
| 1.000000% | 20/20 | [0.000222984944748, 0.0102568704863] | every completion |
| 2.000000% | 0/20 | [0.0362883905405, 0.238134168838] | no completions |
| 3.000000% | 0/20 | [0.270784585979, 0.657985770168] | no completions |
| 5.000000% | 0/20 | [0.872611735139, 0.980926896667] | no completions |

Every-completion rejection is **falsified under stated caps** despite the missing analyses. Some-completion rejection means the conclusion depends on completion; none means **compatible with stated latitude** by the uncertainty rule for every completion. This is not a missing-at-random estimate.

| Reversal cap B | Combined error needed: some to every completion | Equal per-occasion cap needed: some to every completion |
| --- | --- | --- |
| 0.000000% | [2.721088%, 3.741497%] | [1.360544%, 1.870748%] |
| 0.500000% | [2.221088%, 3.241497%] | [1.110544%, 1.620748%] |
| 1.000000% | [1.721088%, 2.741497%] | [0.860544%, 1.370748%] |
| 2.000000% | [0.721088%, 1.741497%] | [0.360544%, 0.870748%] |
| 3.000000% | [0.000000%, 0.741497%] | [0.000000%, 0.370748%] |

#### Separate comparable-kernel paired-discordant test

This restriction requires nondecreasing **marginal detection probability**, a comparable nondecreasing assay kernel, nondecreasing detectable burden and no intervention-induced change. Persistence plus unrestricted occasion-specific errors alone does not imply it.

| Quantity | Recomputed value |
| --- | --- |
| Discordant cases, +/− versus −/+ | 14; 8 versus 6 |
| Conditional reversal fraction | 4/7 = 57.142857% |
| Two-sided 95% exact interval for conditional reversal fraction | [28.860940%, 82.338891%] |
| Exact one-sided p-value for H0: P(+/−) <= P(−/+) | 1619/4096 = 0.395263672 |
| Reject at 0.025 / descriptive 0.05? | no / no |
| Conditional statistical verdict | **compatible with stated latitude** |
| Minimum p over missing-case completions | 10889/65536 = 0.166152954 |
| Minimum-p attaining allocation(s), (++,+−,−+,−−) | (0, 3, 0, 0) |
| Maximum p over missing-case completions | 44923/65536 = 0.685470581 |
| Maximum-p attaining allocation(s), (++,+−,−+,−−) | (0, 0, 3, 0) |
| Completions rejected at 0.025 | 0/20 |

The noisy restriction's nonrejection is retained. PR #58 already permits finite-molecule misses, so positive/negative pairs alone cannot falsify its noisy model. This tests a necessary marginal restriction, not the full growth model. **Unconditional biological/no-intervention verdict: untestable. Individual intervention dates are not verified.**

#### Source reconciliation and limits

The Figure 1 exclusion categories sum to 14, matching the printed 14 ITT exclusions. Those exclusions are outside the analyzed ITT population and do not replace the unsuccessful-analysis sensitivity panel.

Independent error calibration, an individual availability audit, assay-error/shedding/biological-clearance attribution, recurrence sensitivity, a schedule optimum and treatment implications are all **untestable** from these aggregates. The analysis cannot distinguish assay error, shedding or biological clearance; estimate recurrence sensitivity; validate a schedule optimum; or imply a treatment decision.

#### Reproducibility hashes

SHA-256 of UTF-8 text after CRLF and CR become LF. Path keys use `Path.as_posix()`. These hashes check integrity, not source truth or clinical validity.

| Input | Normalized SHA-256 |
| --- | --- |
| `research/ctdna-paired-calls/PROTOCOL.md` | `6e615dcb5eef9a745d4cbd68e03865455f880d2f56ba382b37e15ec9cbc43e93` |
| `research/ctdna-paired-calls/data/aggregates.csv` | `c96bd1f345900bc5987fe4de3e5be7d957ca96fa75ea53198aaedaec09080310` |
<!-- END RECOMPUTED RESULTS -->

## Hostile self-review and falsification record

This is AI-assisted internal self-review, not external scientific review.
The source transcriptions and analytical frontier were checked separately
from implementation; the test oracle independently recomputes the reported
numbers. Agreement and passed tests do not prove source truth or biology.

| Attempted falsification or overclaim | Check and retained interpretation |
| --- | --- |
| Perfect observations plus strict latent persistence | Observed positive/negative mass contradicts this zero-latitude conjunction. It does not isolate which component fails. |
| The frozen tight total cap | Evaluate the lower uncertainty bound, not merely the empirical proportion; retain the full-arm worst-case completion result separately. |
| A stronger constant/comparable-kernel claim | Use both directions of discordance. Preserve the exact test's nonrejection when it survives, rather than presenting observed reversals as a rejection of all noisy persistence. |
| Failure from only one favorable missingness assumption | Enumerate every allocation of the unsuccessful analyses. No missing-at-random imputation and no selection of a preferred completion. |
| Treat printed negative classification as an audited pair count | Cross-check the definition and reporting, explicitly label the inferred cell, and refuse inference if the gate is false. Individual availability remains unaudited. |
| Confuse ITT exclusions with unsuccessful analyses | Preserve separate all-arm exclusion rows and verify their partition independently. |
| The sum bound is merely loose | Construct latent laws that attain it while preserving all observed cells; test frontier and subfrontier caps. |
| Persistence alone implies nondecreasing observed detection | A synthetic law with persistent latent positives but occasion-dependent errors is a counterexample; the comparable-kernel assumption is essential. |
| Finite-molecule misses cannot produce a reversal | A synthetic shared noisy kernel with nondecreasing burden has positive reversal probability. Tests retain this counterexample to the overclaim. |
| A surviving uncertainty cap fits the exact proportions | Report exact plug-in compatibility and statistical nonrejection in separate columns. |
| Infer no intervention from late result release | The Methods describes both draws being analyzed together and results released at weeks 8–10. **Individual intervention dates are not verified.** The unconditional untreated biological interpretation is untestable. |
| Treat public HTML as an open licence | Check publisher rights; retain only attributed factual counts under the narrow recorded reuse basis. |
| Reproducibility implies validation | Hashes protect these local inputs; CI checks code. Neither establishes independent source truth, calibration, transportability or clinical usefulness. |

No invented assay calibration, intervention record or individual truth label
is used to strengthen a verdict. The result tables retain both failing and
surviving restrictions. If the independence/sampling assumptions are not
adopted, the deterministic finite-panel frontier and completion bounds remain
valid, while the population confidence and significance claims do not follow.

## Novelty and limits

The article already reports the paired-call discordance. This analysis does
not discover it. Union bounds and their sharp mass-allocation constructions,
worst-case missing-data bounds, exact binomial confidence limits and the
paired-discordant test are established mathematics. Method history includes
[Clopper and Pearson (1934)](https://doi.org/10.1093/biomet/26.4.404) and
[McNemar (1947)](https://pubmed.ncbi.nlm.nih.gov/20254758/); their publisher
metadata was checked, and no claim of a new testing method is made.

The contribution is limited to the explicit cap frontier, prespecified
uncertainty calculation, missing-analysis sensitivity and source/assumption
audit for this published observation channel. It is a targeted retrospective
exercise, not a systematic review or an independent experimental replication.
No bibliographic priority or new biological mechanism is established.

These aggregates **cannot distinguish assay error, shedding or biological
clearance; estimate recurrence sensitivity; validate a schedule optimum; or
imply a treatment decision**. There is no individual reference standard for
latent disease, no independently measured error or reversal cap, no verified
individual intervention timing, and no matched alternative schedule.
Recurrence outcomes are not substituted for call truth. A nonrejection is
not equivalence, and a rejection is of the stated assumptions together.
All outputs are nonclinical; no treatment or surveillance advice is given.

## Reproduce and inspect

From the repository root:

```bash
python -B research/ctdna-paired-calls/analysis.py
python -B research/ctdna-paired-calls/analysis.py --json
python -B research/ctdna-paired-calls/analysis.py --check
python -B -m unittest discover -s tests -p test_research_ctdna_paired_calls.py -v
```

The [analysis](analysis.py) and [single test file](../../tests/test_research_ctdna_paired_calls.py)
use only the Python standard library and make no network requests. All
reported derived numbers come from the aggregate CSV and the frozen constants.
The tests recompute probabilities and exact-tail inversions independently,
check every missing-case allocation and the full README result block, and
exercise source gates and mathematical counterexamples. SHA-256 values use
newline-normalised content; all path keys use `Path.as_posix()` so checkout
line endings and Windows separators do not alter identity.

The existing repository validation, unit-test and PR checks are used without
changing workflows, dependencies, indexes, policy files or any other research
directory. Runtime checks and any unresolved execution limits are recorded
in the PR description. The entire new test suite is required to run within
one minute; its timing is an execution check, not a scientific result.
