# Frozen protocol: schedule-associated terminal-cause attribution

**Nonclinical preclinical aggregate analysis; unreviewed. No treatment,
diagnosis, schedule or dose advice.** This is a prospective specification of
this repository computation on already public results, not a blinded analysis,
an independent preregistration, or a new animal experiment.

## Freeze and scope

This file will be committed by itself before writing/running inferential
analysis code or calculating contrasts, p-values, intervals or sensitivity
breakpoints. The first commit containing it is the freeze; the final README
will identify that commit and its newline-normalised SHA-256. Do not amend this
file after the freeze. Any necessary deviation must be identified in README,
with the original analysis retained where executable.

Start from `main` at `37c37e2b4cfc711ce938201123ef3085b8868f2f`.
Add only `research/schedule-terminal-cause/` and
`tests/test_research_schedule_terminal_cause.py`. Read first: ranked job 5 in
`research/README.md`, `DATA_POLICY.md`, PR #66's `containment-evidence` note and
PR #61's `containment-at-tumour-scale` note. Existing files remain untouched.

## Source access, rights and observations made before freeze

Source: Seyedi et al. (2024), *Testing Adaptive Therapy Protocols Using
Gemcitabine and Capecitabine in a Preclinical Model of Endocrine-Resistant
Breast Cancer*, Cancers 16:257,
<https://www.mdpi.com/2072-6694/16/2/257>, DOI
<https://doi.org/10.3390/cancers16020257>.
Supplement: <https://www.mdpi.com/article/10.3390/cancers16020257/s1>.
Access date: 2026-11-18 UTC.

Before numerical transcription, the article's closing copyright paragraph was
checked: copyright 2024 the authors, **CC BY 4.0**
(<https://creativecommons.org/licenses/by/4.0/>).
MDPI's <https://www.mdpi.com/openaccess>, "Meaning of Open Access", explicitly
extends reuse to supplementary material, subject to separately credited
third-party rights. The supplement has no separate licence in its extracted
text and no different notice/third-party credit beside the relevant tables.
Its CC BY 4.0 reuse basis is therefore the article licence plus that explicit
publisher policy, not a licence statement printed inside the supplement.
Attribute the authors, link the sources and licence, and identify transcription
and calculated changes. The repository MIT licence does not relicense these
third-party aggregates. No full PDF, image, individual record or trajectory
will be committed.

The supplement was obtained with ordinary HTTP (`urllib.request`), HTTP 200,
redirecting to <https://www.mdpi.com/2072-6694/16/2/257/s1>, and parsed with
`zipfile`. No browser was used and no download prompt appeared. The ZIP's sole
member is `cancers-2752834-supplementary.pdf`. ZIP SHA-256:
`dc22ed5202643894a9a5e636e6936dba075000b01646a876f1030c89dc151af2`;
PDF SHA-256:
`af5371264cd6f066f0af3200dad68d6c33b1f8a7e30fa24df56f7b61f427067b`.
These are hashes of source binary bytes, distinct from newline-normalised
hashes of repository text. PDF printed and physical page numbers agree.

**Denominators are checked before inference.** Table S2 p.32 prints adaptive
capecitabine `(4,5,4)` in toxicity/cancer/other columns, and Table S3 p.32
prints adaptive gemcitabine `(7,1,5)`. Each totals 13; each combines two
allocated groups of 6, hence 12 animals (article Table 1 and section 2.3.1).
Their MTD rows each total 6. Neither excess count is resolved. Table S4 pp.32-33
prints vehicle `(1,3,0)`, pooled ping-pong adaptive `(3,6,3)`, pooled tandem
adaptive `(5,4,3)`, and combination MTD `(6,0,0)`, matching allocated totals
4, 12, 12 and 6. Preserve those source pools: no splitting dose modulation
from intermittent therapy and no combining ping-pong with tandem.

Read Table S1 pp.11-31 for qualitative attribution reasons, the cause narrative
on p.32, Figure S7 p.10, and Table S5 p.34 for printed dose means. Do not derive
new counts from animal rows or digitize curves. S5's gemcitabine cumulative-dose
means conflict with article section 3.3.1 (256 versus 256.9, and 33.3 versus
333.3 mg/kg). Retain both statements with locators; do not correct S5 or use
dose means to reconstruct durations or adjust the cause analysis. The p.32
narrative's reference to "Tables S3-S5" is inconsistent with the actual cause
table labels S2-S4; use the printed table headings.

## Question, estimand and primary contrast

Does the greater fraction **labelled** toxicity in constant combination MTD
versus source-pooled ping-pong adaptive therapy survive a specified attribution
uncertainty envelope? This is the designated primary contrast, motivated by
PR #66's comparator warning and the source's alternating-drug comparisons.
Tandem and vehicle comparisons are mandatory context, not independent
replications or replacements if the primary result fails.

The estimand is terminal-category composition among the printed animals at
their respective observed endpoints. Use the ordered categories T (toxicity),
C (cancer), O (something else). Show counts, denominators, category proportions
and pairwise percentage-point differences. A difference is first named group
minus second. All six fixed S4 pairs are:

1. MTD minus ping-pong (primary).
2. MTD minus tandem.
3. MTD minus vehicle.
4. Tandem minus ping-pong.
5. Ping-pong minus vehicle.
6. Tandem minus vehicle.

There is one S4 table with one combination MTD row and one vehicle row.
Reusing these rows in correlated contrasts does not duplicate animals or
create independent studies. Do not pool single-drug and combination controls.

## Exact tests and multiplicity-aware uncertainty

For every fixed pair compute a two-sided Fisher-Freeman-Halton conditional
2-by-3 composition test. Condition on row and column margins; enumerate all
feasible first rows x. The table probability is
`product(comb(column_total[j], x[j])) / comb(N, n_first)`.
The two-sided tail sums probabilities no greater than the observed table's
probability, including ties, using exact integers/rationals. Remove empty
columns by the same formula; a one-column table has p=1. No asymptotic
chi-square approximation, Monte Carlo draw or mid-p correction is allowed.

Also compute the two-sided Fisher 2-by-2 T-versus-(C+O) test by the identical
probability-ordering rule. The fixed confirmatory S4 family contains **all
12 tests** (six pairs times two endpoints); apply Holm's step-down correction
at familywise alpha 0.05. Report every raw and adjusted p-value, including
null and adverse comparisons. Overlapping controls do not invalidate Holm's
arbitrary-dependence familywise bound. Do not switch to one-sided testing or
reduce the family after inspecting outcomes.

Report conservative simultaneous 95% intervals for all 12 S4 arm-by-category
probabilities: two-sided Clopper-Pearson binomial bounds at individual error
0.05/12. Invert binomial tails numerically with deterministic bisection;
endpoints for zero/all successes are included. Difference envelopes are
`[lower_first-upper_second, upper_first-lower_second]`. These are conservative
simultaneous reference-model intervals, not Fisher-test inversions; the two
procedures may disagree about exclusion of zero. Exact testing and binomial
coverage assume independent, exchangeable categorical sampling within the
source pools. The pooled protocols need not truly be homogeneous, so this is
an explicit statistical reference model, not a design-based randomization
test or a guarantee under arbitrary animal dependence/heterogeneity.

## Attribution envelopes and fragility

The p.32 narrative assumes that spontaneous death, weight loss and other ill
health indicate toxicity, while acknowledging cancer cachexia and unresolved
necropsies. No printed aggregate identifies an adjudicated subset of T.
Walking-related sacrifice can reflect tumour bulk even with low viable burden.
The vehicle T label remains as printed, even though no chemotherapy was given.

Declare two nested, deterministic envelopes, not misclassification estimates
or confidence distributions:

- **T-only:** in each source row independently choose an integer
  `0 <= t <= printed_T` and transform `(T,C,O)` to `(T-t,C+t,O)`.
- **T-and-O:** also choose `0 <= o <= printed_O` and transform the row to
  `(T-t,C+t+o,O-o)`. This is deliberately broader: not all O reasons are
  demonstrated cancer deaths. It tests dependence on calling them "other".

Keep printed C labels fixed. This asymmetry is an explicit stress model for
the authors' assumptions, not proof those C labels are correct. No records
are deleted, moved between arms, or allocated more than once. The reclassified
count is `k=sum(t+o)`, which counts changed terminal labels once. Enumerate
the complete integer envelope and recompute the 12-test Holm correction for
each coherent four-row state; the shared MTD assignment must be the same in
all simultaneous comparisons. Give envelopes for effects, p-values, adjusted
p-values and simultaneous effect intervals, plus worst-case behaviour at
every cumulative global budget k ("at most k"). Give explicit extremizing
witnesses. Fractions of the full S4 sample and eligible labels use the
appropriate generated denominators, not pseudo-replicated pair totals.

For **each pair** separately, hold the other rows at their printed labels and
find the smallest allowed k causing (a) loss of raw p<=0.05, (b) loss of
Holm-adjusted p<=0.05, for both tests, and (c) loss of the printed strict
toxicity-difference direction (tie or reversal). A contrast already lacking
rejection/direction has breakpoint zero; unreachable changes are labelled
unreachable, never zero. Report witness counts and fractions k/(pair animals)
and k/(eligible pair labels). These pair-local breakpoints differ from global
budgets: edits outside a pair can change its Holm-adjusted p-value. For the
designated primary additionally give the coherent global breakpoint using
the full-table denominator. Under an unrestricted T-to-C envelope a positive
toxicity contrast cannot be guaranteed; the quantified breakpoint, not that
logical observation alone, is the substantive result.

## Single-drug denominator failure and conditional bounds

Retain S2/S3 printed counts and excess totals in the CSV and report. The
unconditional single-drug inferential verdict is **untestable** because the
printed adaptive rows cannot be disjoint counts of the allocated animals.
Do not treat 13 as a randomized sample size and do not select a repaired row.

As an explicitly conditional partial-identification sensitivity, suppose
exactly one printed adaptive entry is one too high and all other entries and
allocations are correct. Enumerate every nonnegative one-cell decrement,
without preferring any cause; this is an assumed error model, not a source
reconciliation, and does not cover every possible reporting error. Compare
each resulting adaptive pool against its own MTD group with the same two exact
tests; report the effect and p-value extrema across all completions and then
across their T-only and T-and-O envelopes. Use Bonferroni across the fixed four
single-drug tests (two drugs times two endpoints). They form a separately
labelled exploratory conditional family and cannot rescue the S4 verdict.
No S1 animal-level reconstruction is permitted.

## Decision rule and forbidden interpretations

The primary verdict is **survives** only if S4 passes source/rights/denominator
gates and, for every state in the full T-and-O envelope, MTD minus ping-pong
has positive T difference and both its tests have Holm-adjusted p<=0.05.
It **fails** if S4 is evaluable but any allowed state violates that rule;
state whether it already fails at the printed labels. It is **untestable**
if source rights, printed counts, or reconciled S4 denominators are unavailable.
Report the T-only verdict and every secondary contrast separately. A secondary
contrast survives only if its printed strict T-difference direction and both
Holm rejections persist throughout the stated envelope; an initially tied
direction fails this rule. Failure
means failure of the declared robustness rule, not absence of an association,
proof of no toxicity, or a causal refutation of treatment benefit.

Unequal observation times and competing terminal events prevent interpretation
as equal-horizon cancer mortality, cause-specific survival/hazards or competing-
risk cumulative incidence. The category table cannot reproduce the article's
survival analysis, whose censoring rules differ. It identifies neither causal
drug toxicity, competition between clones, nor PR #61's mutation-supply
threshold. Dose means cannot identify mediation, equal exposure or follow-up.
The source's own Figure S7/p.32 and article section 4.1 already discuss these
patterns and toxicity limitations; novelty is a limited reproducible aggregate
audit, not a new biological discovery or established bibliographic priority.

## Reproducibility and hostile review

Store only small printed aggregate CSVs, with source URL, table, page/section,
row and column locators, access date and reuse basis. Original analysis code
uses Python stdlib only. Recompute all quantitative report blocks from CSVs.
Hash repository text after CRLF and CR are normalised to LF; retain other
bytes and final-newline state. Build all path keys with `Path.as_posix()`.
The sole new offline test file must recompute the reported analysis, check
exact tails/intervals/breakpoints and normalized hashes, and finish under
60 seconds on supported Python 3.12/3.14 platforms. No network in tests.

Before pushing, independently check every transcribed aggregate against the
printed source, attempt to overturn the result by all declared relabellings,
check the unrepaired denominators, interval assumptions, shared controls and
source-text contradictions, and document corrections/falsification attempts
in README. Passing tests, AI review and hashes establish reproducibility only.
Push the branch, open a PR to main, and inspect actual CI without changing
existing workflow files.
