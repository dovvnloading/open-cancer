# Frozen protocol: robustness of the ecDNA rejection-count argument

**Nonclinical, unreviewed computational research. No treatment advice.**

Date: 2026-11-18 (America/New_York). Branch from current main
`9ded69c3586188b173e81fedb7a79abbf4e7e91f`.
This document must be committed before executing any new probability,
frontier, or experiment-design calculation for this contribution. Its content
will then remain unchanged. Record the freeze commit and any subsequent
deviation in README.md; do not amend this file to fit results.

This is a prospective computational specification for a retrospective
sensitivity analysis, not a blinded study or external preregistration. The
three published counts and the independence/Markov endpoints in PRs #63/#65
are already known. The source rows, the ranked third job in research/README.md,
DATA_POLICY.md, AGENTS.md, the E3-X check, and the #55/#60/#65 interpretation
have been read before freezing. Conceptual mathematical criticism before the
freeze is allowed; new numerical analysis is not.

## 1. Question, observations, and exclusions

How much miscalibration of the published nominal 0.01 mother-level tests,
together with dependence between their rejection indicators, suffices to
remove the three-line Holm-adjusted inconsistency with fair segregation?
The target is the robustness of an observation-based argument, not a new
allocation fit, ecDNA dispersion estimate, biological mechanism, or killing
threshold. Failure to reject does not establish fair segregation.

Read only the existing
`research/ecdna-segregation-clustering-data/data/lundberg.csv`, pinned by a
SHA-256 of UTF-8 content with CRLF and CR normalized to LF. Use exactly
`lundberg_chp212_n`, `lundberg_chp212_pct`, `lundberg_lan5_n`,
`lundberg_lan5_pct`, `lundberg_mcix_n`, `lundberg_mcix_pct`, and
`lundberg_alpha`. Require the expected quantities, populations, untreated
condition, DOI, source locator, and provenance fields. Enumerate integers
compatible with a half-percentage-point outward rounding envelope; if
multiple counts were compatible, use the smallest count (largest upper-tail
probability). Expected inputs are CHP-212 2/100, LA-N-5 10/100, MC-IX 1/50.
These are three cell-line panels from one paper, not three new studies.

Source: Lundberg et al. (2008), PLoS ONE, DOI
10.1371/journal.pone.0003099, Results "Statistical distribution of DMs at
mitosis", Figure 3E-G. Reuse the repository's attribution and source-access
date; this run performs no fresh primary-paper retrieval or browser download.
Its tiny numerical facts retain their recorded Creative Commons Attribution
reuse basis (version unspecified), not a new license from this repository.

Exclude colony-inferred divisions, irradiation panels, sister-count Pearson
correlations, medians, and all other sources from these calculations. Do not
reconstruct mother-level observations. Detached-DM exclusion, unspecified
test sidedness/construction, and missing culture/lineage allocation remain
measurement limitations. The nominal level is not verified actual size.

## 2. Count model and sensitivity parameters

For a line with n scored mothers let X_i be the indicator that the source's
whole measurement-and-testing procedure rejects for mother i, and K=sum X_i.
Under fair segregation plus a specified observation pipeline, let
alpha'=p=Pr(X_i=1). Let rho_I=Corr(X_i,X_j), i != j, for exchangeable
indicators within a line. Neither p nor rho_I is estimated by these counts.
In particular rho_I is NOT within-mother ecDNA allocation dispersion, sister
copy-count correlation, or interspecies co-inheritance correlation.

The primary robustness parameters are caps A on actual per-test size and
R on indicator correlation: p <= A and rho_I <= R, with 0 <= A,R <= 1.
Negative feasible correlation is permitted below R; no independence is
inferred from R=0. Allow arbitrary dependence between cell lines. Caps are
shared for the main two-dimensional display but each line may have different
unknown p and rho_I within them. The formula must also accept line-specific
caps without treating them as inferred values.

The primary finite-exchangeable model consists of every count law q_s on
s=0,...,n with the stated moment constraints. Conditional on K=s, a uniformly
chosen s-subset of mothers receives X_i=1. At exact p the moments are

    E K = n p,
    E K(K-1) <= n(n-1)[p^2 + R p(1-p)].

Derive and prove the sharp upper-tail bound over this class. For R>=0 an
independent zero-to-one augmentation should show that maximizing over p<=A
reduces to p=A; prove this rather than substituting a cap without justification.
Use exact rational primal distributions and dual polynomial majorants to
certify any claimed sharpness. If that proof fails, report only a valid
conservative bound and explicitly downgrade any claim of an exact robust
frontier. Do not assume a beta mixing distribution is worst-case exchangeable.

## 3. Required comparison models and bounds

1. **Independent indicators:** K~Binomial(n,p). Compute exact upper tails
   sum_{s=k}^n C(n,s) p^s (1-p)^(n-s). Use p=A for a size cap by a coupling
   argument. This is an additional independence assumption, stronger than
   pairwise zero correlation.
2. **Specified beta-binomial indicators:** for 0<rho_I<1, draw a common
   theta~Beta(p c,(1-p)c), c=(1-rho_I)/rho_I, then independent Bernoulli(theta)
   indicators conditional on theta. The count law is beta-binomial; its
   pairwise indicator correlation is rho_I. Include rho_I=0 as the independent
   limit and rho_I=1 as the all-equal-indicators limit. Compute tails exactly
   for rational parameters or use high-precision calculations checked by
   exact rational endpoint evaluations. This is pointwise model sensitivity,
   not robustness to every dependence distribution. Do not assume the tail
   is monotone in rho_I; retain any loss and reappearance of rejection.
3. **Sharp exchangeable moment envelope:** optimize the count tail at the
   moment caps in section 2. Candidate dual majorants are the constant 1,
   s/k, and (s-j)(s-j-1)/[(k-j)(k-j-1)] for j=0,...,k-2. Prove completeness
   using feasible extremal laws, or use a complete exact finite LP/vertex
   calculation. Account for k=1, p endpoints, integer branch ties, and
   mean nA >= k. Check both normalization and all moment/tail certificates.
4. **Dependence-free envelope:** min(1,nA/k) by Markov, with an explicit
   attaining count law. This needs marginal-size control, but neither
   independence nor exchangeability for validity. Demonstrate that the
   exchangeable envelope reaches it when the correlation cap becomes
   uninformative. Markov non-rejection is not a positive biological finding.

Prove the validity of using these tail bounds as conservative p-values under
their declared null classes. Finite numerical checks supplement the proof;
they are not a proof of an all-parameter assertion.

## 4. Family, decision rule, and meaning of jointly

The only primary multiple-testing family is the three untreated cell lines.
For each model/cap scenario compute all three p-values and the Holm step-down
adjusted values with familywise level delta=0.05. Sort deterministically,
apply multiplicities 3,2,1 and cumulative maxima, and cap at 1. Reject a line
exactly when its adjusted value <= delta; equality remains rejection.
Use exact comparisons for this decision. Prove familywise control under
arbitrary cross-line dependence and the validity of the input tail p-values.

"Joint inconsistency" means rejection of the intersection that all three
line nulls hold: at least one Holm rejection, equivalently
3 min_i p_i <= delta. Report which lines reject, and separately whether all
three reject. It does not mean all three lines reject, nor a pooled fourth
test. Do not add the 13/250 pool, Fisher combination, or PR #63's different
17-comparison family to this primary family. State that this is a targeted
retrospective sensitivity family, not a correction for the repository's
entire research programme or the paper's undisclosed analysis choices.

## 5. Exact frontier and numerical reporting

Define the full two-dimensional rejection regions by the exact tail/bound
formulas and Holm rule, not a visually interpolated grid. For each line and
for the any-rejection intersection test, compute the largest size A still
rejected at each fixed R/rho_I, using monotonicity in A only after proving it.
Give exact rational or algebraic definitions and outward rational brackets
of width <= 10^-12 for nonrational numerical boundary values. This width
describes arithmetic resolution, not sampling uncertainty.

The fixed dependence slices are

    0, 0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02,
    0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.

At nominal A=0.01 obtain the full sharp-cap dependence boundary, with an
exact expression/certificate where possible. For beta-binomial pointwise
sensitivity, isolate every nominal-size crossing on [0,1], including any
re-entry near complete dependence; use a complete polynomial root count
(e.g. rational Sturm isolation) if claiming all crossings. If a computation
is only local or sampled, label it so and do not claim completeness.
Report the independent-size threshold and dependence-free threshold.
Verify any simple local frontier formula over its declared domain rather
than promoting a local active branch to a global expression.

All new derived numerical claims belong in a generated README block or are
explicitly checked against that block. Prefer exact fractions for small
quantities and boundary conditions. The analysis CLI must print its report,
emit structured results, and check README without rewriting it. Hashes cover
the source CSV, frozen protocol, and relevant original code; avoid recursive
self-hashes. Normalize newlines before hashing and form every relative path
key using Path.as_posix(). No timestamp or platform-specific path should make
the report nondeterministic.

## 6. Falsification, experiment, novelty, and delivery

Attempt to falsify: independence inferred from zero pairwise correlation;
monotonic beta-binomial tails in correlation; treating fitted indicator
correlation as physical clustering; fixed Holm ordering; dropping weak lines;
treating a conservative bound as an attained supremum without a witness;
pooling heterogeneous/correlated mothers; using nominal size as calibration;
and claiming that non-rejection proves fair segregation. Preserve failures
and mathematical counterexamples in README, with exact checks where possible.

Translate the frontier into a prospective nonclinical experiment: predeclare
the three-line family, use independent culture/lineage sampling or known
blocks, calibrate the complete per-mother measurement/test pipeline on
separate null controls, and release aggregate conditional maternal totals,
paired-daughter counts, lost-copy outcomes and culture labels. A hypothetical
100-mother-per-line design may report the smallest critical rejection count
under declared size/dependence caps using the conservative first Holm step.
Calibration uncertainty needs simultaneous coverage/error budgeting; a
confidence bound cannot silently be treated as known. Do not invent achieved
experimental precision, power, or a treatment recommendation.

The intended contribution is a repository-specific robustness audit and
reproducible mathematical certificates. Binomial/beta-binomial models,
exchangeable count representations, moment bounds, Markov and Holm are
established ideas. No global priority claim, new biological discovery,
segregation estimate, or empirical validation of #55/#60 is proposed.
Retain #65's opposing evidence and distinction between source fits and
sampling uncertainty. Unconditional biological rejection remains unestablished
if the required calibration/dependence assumptions cannot be verified.

Create only this directory's README.md, this PROTOCOL.md, analysis/verifier
code, and tests/test_research_ecdna_rejection_robustness.py. Modify no existing
file, input data, workflow, index, or other research directory. The new test
must use only the Python standard library, run offline in under 60 seconds,
recompute every reported number, independently challenge probability and
optimization code, and test newline-normalized hashes and POSIX path keys.
Perform hostile mathematical/code review before publication, downgrade every
unproved claim, run the applicable repository gates, push this branch, open
a PR against main, and inspect its actual CI status. Passing CI or AI review
does not establish source truth, scientific priority, or biological validity.
