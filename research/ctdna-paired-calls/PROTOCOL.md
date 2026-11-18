# Frozen protocol: paired postoperative ctDNA calls

**Nonclinical, unreviewed research. No treatment, diagnosis, surveillance or
patient-selection advice.**

Protocol date: **2026-11-18 (America/New_York)**. Base: current `main`,
`620c1df128ce3a30dbaaa190179cb5404048174a`. The first commit containing this
file must contain this file alone, before writing or running the numerical
analysis or transcribing its CSV. Preserve that commit and this file unchanged;
identify any later deviation in the README. This is a pre-computation freeze
for a retrospective, deliberately selected source. The task, research index
and source inspection already disclose the motivating counts. It is neither
blinded testing, independent replication nor registered preregistration.

Read first: [ranked job 3](../README.md),
[PR #58](../ctdna-surveillance-schedule/),
[PR #63](../prediction-confrontation/) and [DATA_POLICY.md](../../DATA_POLICY.md).
No existing file or other research directory will change. The deliverable is
this directory (README, protocol, one small aggregate CSV and analysis code)
and `tests/test_research_ctdna_paired_calls.py` only.

## 1. Question, source and eligibility

How much total latitude in two imperfect observations and genuine threshold
reversal is required to make a persistent-threshold model compatible with
the published DYNAMIC week-4/week-7 calls?

Only new empirical input: Tie, Wang, Lo et al., *Nature Medicine* (2025),
[DYNAMIC five-year report](https://pmc.ncbi.nlm.nih.gov/articles/PMC12974608/),
DOI [10.1038/s41591-025-03579-w](https://doi.org/10.1038/s41591-025-03579-w).
Read article **HTML text only**. No article PDF, source-data file, supplement,
patient table, survival reconstruction or individual figure points. Deny a
browser download prompt if encountered, continue without it, and record the
limitation. Do not fetch files to work around a prompt or access restriction.

Eligible locators are Results, **OS, DSS and RFS in the ctDNA-guided
population**, first paragraph; **Extended Data Table 1**, ctDNA-status rows
and arm denominator; **Figure 1 caption**; and Methods, **Study design and
participants**, especially the paragraph beginning with allocation to
ctDNA-guided management. Pin the access date and article version in README.

The disclosed candidate aggregates are 31 positive/positive, 8
positive/negative, 6 negative/positive, 45 positive at either occasion,
291 successful analyses, 294 guided-arm ITT cases and 246 classified negative.
These are disclosed inputs, not computed results. Independently recheck all
against HTML before committing the CSV. CSV rows retain printed labels,
source URL, DOI, population, units, precise locator, access date, article
rights, reuse basis and interpretation. Store printed facts only, not derived
denominators or a fabricated directly printed negative/negative cell.

**Paired-availability gate.** A statement that samples were collected alone
does not establish successful calls. Cross-check the reported successful
analysis population, concurrent analysis of both draws, explicit definition
of negative at both draws, the positive-pattern partition and the table's
negative classification. If these support complete *reported* paired
classification, analyze the 246 as an explicitly inferred negative/negative
cell at that reporting level. State that no individual assay-availability
audit or separate explicit complete-pairs count was supplied. If the
cross-check fails or reveals single-timepoint classification, do not silently
fill that cell: mark full-table inference untestable and retain only directly
supported discordants with selection caveats. Do not conflate Figure 1's ITT
exclusions for missing draws with unsuccessful analyses inside the ITT arm.

**Rights gate.** Public PMC access is not an open licence. Check the article's
copyright/licence and publisher rights text. The inspected publisher rights
statement retains Springer Nature/licensor exclusive rights; no Creative
Commons licence has been established. Retain only isolated factual numerical
counts in an original schema with attribution, not article prose, source
table arrangement, figures or a source database. The reuse basis is factual
information, as distinguished from protected expression in the
[U.S. Copyright Office FAQ](https://www.copyright.gov/help/faq/faq-protect.html),
not an asserted publisher redistribution grant or MIT relicensing of the
article. If that narrow basis cannot be maintained, retain links and record
an untestable rights exclusion instead of copying protected material.

## 2. Latent model and the sharp cap frontier

For the stated target population, let `Z4,Z7` be latent binary threshold
status at the two occasions and `Y4,Y7` the observed calls. A persistent
threshold has `Z4 <= Z7` almost surely. Define the **unconditional** rates

```
e4 = P(Y4 != Z4)
e7 = P(Y7 != Z7)
b  = P(Z4 = 1, Z7 = 0)
q  = P(Y4 = 1, Y7 = 0)
r  = P(Y4 = 0, Y7 = 1).
```

Here `b=0` means persistence; permitting `b>0` explicitly relaxes it. All
probabilities use the same population denominator. These are not conditional
false-positive rates, conditional false-negative rates, analytical limits of
detection or error rates among week-4 positives. No independent, symmetric,
constant-across-occasion or nondifferential errors are imposed. Dependence
between errors, latent state, burden and the other call is allowed.

On the event `Y4=1,Y7=0`, either occasion is misclassified or there is a
genuine reversal. Therefore

```
{Y4=1,Y7=0} is a subset of {Y4!=Z4} union {Y7!=Z7} union {Z4=1,Z7=0}
q <= e4 + e7 + b.
```

For **upper caps** `E4,E7,B` on these rates, the sharp frontier over all joint
latent/observed laws with the specified observed table is

```
compatible  <=>  E4 + E7 + B >= q
minimum combined error allowance at reversal cap B = max(0, q-B)
minimum equal per-occasion cap at reversal cap B = max(0, q-B)/2.
```

Proof of sufficiency: split the observed `+/-` mass into nonnegative
`a4+a7+ab=q`, each within its corresponding cap. Assign latent `00` on the
first piece (only the week-4 error), `11` on the second (only the week-7
error), and `10` on the third (only reversal). On all other observed cells,
set latent equal to observed. This reproduces the entire observed table and
uses rates `a4,a7,ab`. It proves attainability, not a plausible biological
mechanism or a fitted assay. The inequality alone is **not** a characterization
of arbitrary specified *exact* error-rate triples. For probability caps larger
than needed, no requirement forces those rates to be saturated.

Report the exact rational plug-in frontier, plus the lower-bound frontier
below, as separate objects. Include reversal-cap sections `B` in
`{0, 0.005, 0.01, 0.02, 0.03}` and equal-error-cap sections. These hypothetical
allowances are not independently measured assay performance. In particular,
when reversal is unrestricted the counts need not force observation error.

## 3. Frozen uncertainty, failure rule and verdicts

Use the independent, identically distributed **case-pair** sampling model for
population uncertainty, conditional on the eligible population. Within-case
independence of calls is not required. A finite selected trial panel is not
a probability sample of all patients; population transport is not established.
Do not infer population error calibration from the trial randomization.

For `x` observed `+/-` cases among `n`, the one-sided Clopper-Pearson lower
bound `L_alpha(x,n)` is zero for `x=0`; otherwise it is the unique solution to

```
sum_{j=x}^n choose(n,j) L^j (1-L)^(n-j) = alpha.
```

Freeze two test families: the primary cap restriction at hypothetical
`C0=0.01` and the separate noisy directional restriction in Section 5.
Bonferroni allocates `alpha=0.025` to each, controlling their combined
familywise error at at most 0.05 under their sampling assumptions. Thus the
decision lower bound is **one-sided 97.5%**. Also report the one-sided 95%
lower bound (`alpha=0.05`) as an explicitly marginal descriptive sensitivity.
Report the ordinary two-sided 95% Clopper-Pearson interval for `q`; its lower
endpoint equals the decision bound, but its upper endpoint has a different
purpose. Compute bounds by tail inversion, no asymptotic approximation.

**Failure observable:** reject a hypothetical cap `C` precisely when
`L_0.025(x,n) > C`, equivalently an upper binomial tail at `C` below 0.025
(boundary equality is retained). The primary cap is 0.01. Also display the
prespecified nested sensitivity caps `{0, 0.005, 0.01, 0.02, 0.03, 0.05}`.
All sensitivity decisions come from the same lower bound, so they do not
constitute a search over unrelated tests or additional independent evidence.
Report exact tail probabilities at the listed caps.

Use explicit verdicts:

* **Falsified under stated caps:** the lower bound exceeds the cap under the
  eligible sampling/observation conjunction; not biological falsification.
* **Compatible with stated latitude:** the population restriction is not
  rejected with the stated cap, or a constructed law fits the stated
  probability table. Distinguish sampling nonrejection from *exact plug-in*
  compatibility: a cap below `x/n` can fail exact matching without failing
  the uncertainty rule. Neither type of compatibility validates the model.
* **Untestable:** missing eligibility, independent rate calibration or design
  information prevents the stronger interpretation being tested. Absence of
  evidence does not validate it.

No study-wide biological confirmation follows from these conditional tests.

## 4. Missing analyses without missing-at-random assumptions

Derive `m=N-n` from printed guided-arm and successful-analysis counts. The
`n` successful cases and the full `N` cases are distinct target populations;
rates `e4,e7,b` must be defined separately for each. For failures, hypothetical
completed binary outcomes are a sensitivity device, not recovered calls.

Enumerate every nonnegative integer allocation of `m` to `++,+-,-+,--`.
Without any restriction on why analyses failed, the sharp finite-panel
frontier range is

```
q_full in [x/N, (x+m)/N].
```

Its endpoints are achieved by assigning no failures or all failures to
`+/-`. Report how much error is sufficient for at least one completion and
how much allows every completion. Apply `max(0,q-B)` and the equal-cap
formula to either endpoint where needed. These are identification bounds,
not confidence intervals and not a missing-at-random estimate.

Under the full-arm iid pair model, the conservative lower bound valid for
every completion is `L_0.025(x,N)`; the range of completion-specific lower
bounds is `[L_0.025(x,N), L_0.025(x+m,N)]`. The outer two-sided 95% interval
is `[L_0.025(x,N), U_0.025(x+m,N)]`, where `U` is the upper binomial bound.
For any actual completion its exact interval lies inside this envelope,
so arbitrary missingness is allowed. State that sampling assumptions remain.
Report rejection for every completion, some completions, or none for each
frozen cap. Repeat the directional test under every missing-case allocation.

## 5. Separate comparable-kernel directional test

Freeze the stronger observation model separately: for each case, the same
nondecreasing detection kernel `g(B)` is used at both occasions (including
fixed relevant assay input/background covariates), detectable burden is
nondecreasing, and there is no intervention-induced change between draws.
It implies nondecreasing marginal detection probability, hence

```
P(Y7=1)-P(Y4=1) = r-q >= 0.
H0: q <= r,   H1: q > r.
```

This restriction **requires nondecreasing marginal detection probability**.
Latent persistence with arbitrary occasion-specific errors alone does not
imply it. The same assay name, concurrent processing or increasing tumour
burden alone need not keep input molecules, background DNA, shedding or the
effective detection kernel comparable. The independent-case sampling model
is also not implied by the mathematical kernel restriction.

Condition on discordants `d=x+y`, with `x=+/-` and `y=-/+`. Under an iid
pair law and `H0`, the conditional `+/-` fraction is at most one half. The
one-sided exact paired-discordant (McNemar/binomial) p-value is

```
p = sum_{j=x}^d choose(d,j) / 2^d.
```

Use the frozen alpha 0.025 (also state whether 0.05 would reject). With no
discordants use p=1. Report the rational p-value, discordant fraction and
two-sided 95% exact interval for that fraction. Do not substitute two
independent-proportion tests or discard the `-/+` discordants. Retain the
nonrejection if it survives. Enumerate all full-arm missing-case completions
and report the minimum and maximum p-values and attaining allocations.

Because PR #58 already permits finite-molecule misses, observed `+/-` pairs
alone do not falsify its noisy model. This test concerns a necessary
restriction, not its full growth law or an optimal schedule. Methods says
both draws were processed together and results released later; **individual
intervention dates are not verified**. Report the mathematical test
conditionally, and mark an unconditional biological/no-intervention verdict
untestable when those design controls cannot be confirmed.

## 6. Falsification, novelty and reproducibility commitments

Try to falsify the zero-error/zero-reversal special case, the primary
hypothetical cap, the noisy directional restriction, the robustness to all
missing-case completions and any overclaim about complete availability.
Test the sharpness construction at frontier and subfrontier caps, and
provide counterexamples to deriving directional monotonicity from latent
persistence alone. Include a shared noisy-kernel example with a nonzero
positive/negative probability. Retain contrary and nonrejecting results.

The source already reports the discordance. Union bounds, cap allocation,
worst-case missing-data bounds, exact binomial confidence limits and the
paired-discordant test are established methods, not claimed new mathematics.
The limited contribution is a reproducible quantitative compatibility audit
of this published observation channel beyond PR #63's lead-time summaries.
No new biological finding or bibliographic priority will be asserted.

This cannot distinguish assay error, shedding or biological clearance;
estimate recurrence sensitivity; validate a schedule optimum; or imply a
treatment decision. It supplies no treatment or surveillance advice. No
individual intervention dates, error calibration, truth labels or recurrence
reference standard are available for this question.

Use Python standard library only; no network at analysis or test time.
The single test file must independently recompute every reported derived
number from the CSV and compare the entire generated README result block.
It must complete within 60 seconds, check malformed/absent data gates,
binomial endpoints/tail equations, sharpness and exhaustive missingness.
Hash newline-normalised UTF-8 content (`CRLF` and `CR` converted to `LF`),
with all path keys made using `Path.as_posix()`. Exercise Windows path and
line-ending portability; use syntax/APIs supported by Python 3.12 and 3.14.
Source truth is checked manually against HTML and is not proven by hashes,
reproducibility, AI agreement or CI. Re-read each transcribed count, audit
the final scope, run the existing applicable gates without modifying them,
then push the branch and open a PR with running CI.
