# Frozen protocol: GALAXY false-positive occupancy denominator audit

**Nonclinical, unreviewed research. No surveillance, treatment, diagnosis, patient selection, or assay-validation advice.**

Source inspected before freeze: Nakamura et al., *npj Precision Oncology* 10:70 (2026), DOI 10.1038/s41698-026-01277-5, article HTML, Results "Association of ctDNA status in the post-definitive treatment surveillance window with DFS", Fig. 3a and Methods "Cohort characteristics and selection" / "Statistical analysis". Access: 2026-11-18. The motivating counts were already published; this is a computational ordering safeguard, not independent or blinded preregistration. Parent must be current main at branch creation. This file is the sole addition in the initial protocol commit and MUST NOT be modified after results are examined.

## 1. Mandatory denominator gate (before model or enumeration)

Define the potential estimand as the 126 *nonrecurrent* persons' 637 eligible plasma tests in the same post-definitive-treatment surveillance window. Four independently necessary source conditions:

G1. The reported 116/126 patient-level specificity has 126 nonrecurrent participants as its denominator, with an ever-positive/serial-negative classification on those participants. Do not equate these persons with the separately printed 126 *serial-negative* participants: 10 of the latter recur.

G2. The sample-level specificity 619/637 explicitly refers to 637 eligible tests taken from **exactly those G1 nonrecurrent persons**, rather than, e.g., all negative assay calls, a different QC/outcome subset, or samples contributed by recurrent people.

G3. Every sample counted in G2 belongs to the same defined surveillance window as G1, and the printed range 1-8 tests per person applies to every G1 participant's *G2-eligible* tests, not just to the 190-person panel.

G4. One patient-level ever-positive in G1 means at least one of their G2-eligible samples was positive; there is no extra confirmation/repeat-positive criterion, different test inclusion, or post-event sample filtering that changes the mapping.

Use only the source's HTML and Fig. 3a caption/visible article text. A standard meaning of the word "specificity", repeated denominators, or algebraically coherent fractions does not by itself certify G2-G4. Missing proof of any gate is **UNTESTABLE**; record what was and was not printed, stop, and do not infer an occupancy p-value or within-person dependence. An explicit contradiction also fails the gate. Never synthesize test histories or fetch supplemental participant records. Data provenance and copyright are audited separately.

## 2. Conditional model: execute ONLY if all G1-G4 are source-supported

Under a successful gate, transcribe printed counts without importing individual-level data. Let N=126 nonrecurrent persons, S=637 eligible tests, m=S-619=18 positive sample calls, and k=N-116=10 people with at least one such call. The overall 190-person / 816-sample / 1-8 range is NOT permission to assign S tests to nonrecurrent participants without G2-G3. Retain all reported source-rounding and QC limitations.

Enumerate the complete, unordered integer allocation set using multiplicities a1,...,a8: each a_j >= 0 integer, sum_j a_j=N, sum_j j*a_j=S. This gives all possible 1-8 test-count vectors up to patient permutation. Iterate a1 through a7 with exact pruning, solve a8, and verify constraints; do not substitute the median, average, a single imputation or random simulation. If empty, classify the allocation restriction as falsified conditionally on G1-G4; never treat emptiness as evidence of biology.

Null H0: independently on each eligible sample, the false-positive indicator is Bernoulli(p) with one common p in [0,1], exchangeable across samples and times, irrespective of person. Conditional on exactly m false-positive calls, every m-subset of S tests is equally likely; p cancels. For each allocation a, the exact occupancy law for K (distinct persons hit) is

    P_a(K=r | M=m) =
      coefficient of x^m y^r in product_{j=1}^8
        [1 + y*((1+x)^j - 1)]^a_j
      divided by binomial(S,m).

Every coefficient and denominator is an integer. Use integer polynomial dynamic programming, exact Fraction probabilities and independent small-case brute-force oracles; never compute binomial probabilities at a plugged-in estimate of p. K=10 is a *conditional* observation, not proof of excess correlation on its own.

Alternative H1: within-person test calls are not independent/exchangeable; repeated false-positive calls may share a person-specific propensity, temporal persistence or other dependence. The test identifies neither a mechanism nor a treatment effect.

## 3. Frozen tails and family-wise error

Two directional tests are in one family: excess concentration (small K; PRIMARY) and excess dispersion (large K; descriptive competing direction). Family-wise alpha=0.05 via Bonferroni 0.025 to each direction, inclusive tails. For EVERY feasible a compute exactly

    L_a = P_a(K <= k | M=m); U_a = P_a(K >= k | M=m).

A direction rejects only if its tail is STRICTLY less than 1/40 for *every* feasible allocation. Equality does not reject. Report the sharp min/max tail envelope, attaining allocation witnesses, number of feasible allocations, and how many allocations would reject, not just a chosen schedule. Requiring rejection for all nuisance allocations is an intersection-union robustness rule: there is no additional multiple-testing penalty for searching allocations, but both directions share the frozen Bonferroni family.

Verdicts: **FALSIFIED** only when the primary small-K direction robustly rejects the independent common-p null after the family rule (or the conditional allocation constraints are empty). **CONSISTENT** when the source gate passes and neither direction rejects for any feasible allocation; this means nonrejection only, not evidence of iid truth. **INCONCLUSIVE** when gate passes but tail conclusions vary by allocation, or another required test precondition remains uncertain after denominators are settled. **UNTESTABLE** on any failed/undocumented G1-G4 gate; this takes precedence over all mathematical verdicts.

## 4. Sharp concentration report, if the gate passes

Given m positive calls occupying k patients, the exact model-free excess beyond the first call is C=m-k; any call allocation has C repeat calls, at least ceil(C/7) multiply-hit patients when no one has >8 tests, and at least C within-person unordered pairs of positive tests. Establish attainability subject to ALL feasible test-count vectors before calling any bound sharp. If independence is robustly rejected, report the sharp feasible lower bound and a concrete *aggregate* multiplicity witness without manufacturing patient sequences; distinguish concentration of realized calls from proof of a correlation parameter. If not rejected, do not assert the calls "must cluster" statistically.

No result can identify clonal haematopoiesis, patient-specific methylation/background, assay contamination, true occult disease subsequently cleared, or any other etiology. Follow-up censoring, QC and selection may change what a "false positive" means. Statistical robustness would not validate the assay or guide surveillance.

## 5. Reproducibility and mandatory stop

Use only offline standard-library text/CSV/JSON, no private data or file downloads. Store sparse attributed printed counts and source locators only. Normalize CRLF and CR to LF before SHA-256; use Path.as_posix() for hash-manifest path keys. Use exact rational arithmetic for fractions and tails; any float display must be checked against declared rounding/error bounds. Test on Python 3.12 and 3.14 where available and with Windows-style path/line-ending fixtures. A gate failure means the implementation returns UNTESTABLE and **must not call allocation enumeration or occupancy code**. A later resolution requires a new, explicitly identified protocol revision/PR, not a silent change to this frozen file.
