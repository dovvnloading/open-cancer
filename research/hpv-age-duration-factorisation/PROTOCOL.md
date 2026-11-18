# Frozen protocol: Swedish HPV age-by-duration factorisation audit

**NONCLINICAL, UNREVIEWED RESEARCH. No vaccination, screening, diagnosis,
treatment or clinical advice.**

Frozen 2026-11-18 (America/New_York), on a branch created from main at
307aa3c336548046c344772f3bedd4f325a1b1a0. This document must be the
ONLY file in the first research commit. Commit it before running any result
calculation or writing the derived data, script, or report. Record the freeze
commit in README. This is a retrospective, source-disclosed protocol freeze,
not independent preregistration or a blind test. No prior repository file
may be edited. Only the new research/hpv-age-duration-factorisation/ files
and tests/test_research_hpv_age_duration_factorisation.py are permitted.

## Literature/novelty gate (survey before computation)

Wu, Deng, Lepp et al., *BMJ* (2026), DOI
10.1136/bmj-2025-087326, 25 February 2026, article HTML:
https://www.bmj.com/content/392/bmj-2025-087326.full.
Use Results > "Invasive cervical cancer by time since vaccination" and Fig 3
with Methods > "Statistical analyses". The article already reports reduced
invasive cervical-cancer incidence, age-at-initiation contrasts, time-since-
vaccination strata and no observed waning. Its reported Methods and Results
do not give a fitted age-at-dose by duration interaction coefficient, a
homogeneity p-value, or a test of the three age-ratio contrasts' equality.
The earlier Swedish cohort (Lei et al., NEJM 2020,
https://pubmed.ncbi.nlm.nih.gov/32997908/, DOI 10.1056/NEJMoa1917338)
already reports age-specific vaccine associations, but not this 7-9 /
10-12 / 13-15-year three-bin homogeneity test. This is a limited literature
screen, not proof of global bibliographic priority. If an existing formal
test of the same ratio constancy is located, STOP the broad novelty claim:
the surviving narrower question is only whether **printed marginal-interval
boxes**, with disclosed rounding, admit a shared ratio; label it as a
reanalysis of a settled question and explain the pivot.

The BMJ article's HTML open-access notice explicitly specifies Creative
Commons Attribution **CC BY 4.0**:
https://creativecommons.org/licenses/by/4.0/.
Attribute authors, journal, DOI and source URL for every sparse numerical
record. Repository MIT does not relicense the underlying data. Do not save
paper images, prose, patient records, model code, supplemental files or other
non-sparse content. Read article HTML only. Do not request or accept browser
file downloads; deny any download prompt and log resulting access limits.
No register-level/patient-level or unvaccinated person-time is available here.

## Target restriction and like-for-like rows

Let E_t be Wu et al.'s fully adjusted vaccinated-versus-unvaccinated
incidence rate ratio (IRR) for vaccination at 10-16 years, and L_t the
corresponding fully adjusted IRR for vaccination at age 17 or later, both
within the SAME time since first dose, t. The necessary multiplicative
factorisation restriction under audit is E_t = q * L_t for one finite
positive q and all three selected bins. Equivalently Q_t = E_t/L_t = q.
This is a descriptive restriction on printed **adjusted relative measures**;
it is not a direct cohort-matched vaccinated-only comparison, biological
constant, causal effect, or a full time-by-age interaction likelihood test.

Freeze exactly three matched pairs from Fig 3 / the Results subsection:
7-9, 10-12 and 13-15 years since first quadrivalent HPV vaccine dose, with
age at first dose 10-16 versus >=17. They use the same figure, fully adjusted
regression definition, endpoint and duration partitions. The published
marginal 95% CIs are for EACH adjusted IRR, not for Q_t. Their dependence and
the adjusted model covariance are unavailable. Retain every pair even if
sparse. Do not use overall age-at-dose IRRs, separate Fig 3 13-18 pooled
estimates, supplementary buffer results, 16-18 with a missing early IRR,
or birth-cohort Table 3 as extra matched cells. Early 1-3 and 4-6 and
16-18 rows have zero vaccinated cases and no published early adjusted IRR;
these are **outside** the primary panel, not hidden zeros in the ratio.

Transcribe precisely the twelve printed estimated/CIs for six strata as
two-decimal strings in an attributed CSV. The user-disclosed six vaccinated
case/person-year pairs (early 5/797583, 11/668369, 8/212315; later
26/264588, 18/195405, 4/87765 in the same bin order) are descriptive only.
They are NOT denominators for the unvaccinated reference and must NEVER be
used to derive/reconstruct adjusted IRRs, standard errors, likelihoods or a
cross-age causal contrast. There are no fabricated source cells.

## Printed precision, exact rational compatibility

Each printed decimal token x with d fractional digits denotes the closed
last-digit rounding band [x-h(x), x+h(x)], where h(x) = 1/(2*10**d).
The printed IRRs and CI endpoints have d=2, so h=0.005, but implement
this from the literal token. Treat both end bounds as inclusive and do not
round before deciding. Validate every point and ordered interval, and that
each point band lies within the widened CI. Nonpositive lower denominator
bounds are not safe to divide; a missing, dash, zero or nonpositive adjusted
IRR triggers the untestable gate rather than epsilon replacement.

For a bin, write early lower/upper CI endpoints a,b and late lower/upper
endpoints c,d (each printed 95% marginal CI). Under the primary rule the
ratio's inclusive feasible *interval box projection* is

    R_t = [(a-h(a))/(d+h(d)), (b+h(b))/(c-h(c))].

The lower denominator must be positive. Use fractions/decimal strings,
with exact rational arithmetic for both endpoints and all comparisons.
This box treats each published CI as if freely selectable independently.
It is a deliberately conservative **compatibility projection** of marginal
information, NOT a 95% CI for Q_t or an estimate of covariance. Display
the printed point ratio E_t/L_t separately, as a descriptive exact fraction
and rounded decimal. Do NOT replace R_t by point-ratio intervals or by
ratios of vaccinated cases/person-years.

Let A = max_t lower(R_t), B = min_t upper(R_t). If A <= B the common-ratio
compatibility set is [A,B] and the **primary rule PASSES** (nonempty).
If A > B it **FAILS** (empty). Record attaining bins; equality passes.
No confidence or statistical rejection is implied by either outcome.

The sharp symmetric *log-ratio* allowance is

    epsilon* = max(0, [ln(A) - ln(B)] / 2).

Equivalently the smallest epsilon >=0 making
[A*exp(-epsilon), B*exp(epsilon)] overlap. It expands EACH bin's
[ln lower(R_t), ln upper(R_t)] by epsilon in both directions around a
common log q. When A>B, the attaining log reference is
(ln(A)+ln(B))/2; smaller epsilon cannot suffice, and that reference
attains the bound. When A<=B, epsilon*=0. Ratios and A/B are exact
fractions; logs and exponentiation are diagnostic floats, verified with
tight documented absolute bounds (1e-11 in log units). Also report
exp(epsilon*) and the point-only allowance separately; do not confuse
the half-gap with a one-sided gap, a p-value, or identified biology.

## Counts, zero events and sparse-cell policy

Provide separately the **exact-count** two-sided 95% Poisson mean
identification interval for each of the six descriptive vaccinated counts,
conditional on its printed person-years as a fixed exposure:
for x>0 lower lambda uniquely solves P_Pois(lambda)(X>=x)=0.025 and
upper uniquely solves P_Pois(lambda)(X<=x)=0.025; for x=0 the interval
is [0,-ln(0.025)]. Divide both endpoints by person-years; display per
100,000 vaccinated person-years. This is exact-tail count-model inversion,
numerically root-solved to 1e-9 mean events or tighter; it is NOT an
exact interval for a fully adjusted IRR. Count uncertainty is NOT added
to the adjusted-IRR 95% CI (no double use of cases). Do not use continuity
corrections or equate crude rates with regression-adjusted IRRs.

For any zero-event or absent age-duration adjusted IRR (including 1-3,
4-6, 16-18 for the early group), the Poisson count interval gives a
finite crude-rate upper limit if exposure is known but **does not identify
an adjusted IRR or Q_t** without the unvaccinated adjustment/reference.
Represent its adjusted ratio identification set as [0,+infinity]/
**untestable**, never Q=0 or an invented point. In the prespecified
out-of-primary zero-bin diagnostic, the three excluded duration bins
remain untestable; this cannot overturn the restricted three-pair verdict.
If any PRIMARY pair has a zero-event missing IRR, fail the primary
data-adequacy gate as UNTESTABLE rather than silently dropping that pair.

## Prespecified sensitivities and multiplicity-aware variant

Run and report all, irrespective of outcome:

1. Primary: 95% printed CI endpoints, widened by one HALF last digit.
2. No-rounding: same printed marginal CIs, exact as printed.
3. Full-unit widening: replace h by 1 full last printed unit (0.01).
4. Point-only: for each row use the printed IRR's own rounding band,
   not its CI; project their ratio bands. This is a stringent **point
   compatibility** diagnostic, not a confidence statement. Compute its
   sharp log allowance.
5. Leave one bin out (three fixed pairs of remaining bins). Report exact
   intersections and allowances, never pick a deletion to change primary.
6. Approximate multiplicity-aware family: **six** underlying IRRs, alpha
   0.05, two-sided Bonferroni target alpha/6 per marginal, and a log-Wald
   approximation only. Let p-,p+ be the point rounding endpoints; a- and
   b+ the source CI lower/upper endpoints widened by half digit. Define
   D=max(ln(p+/a-), ln(b+/p-)) and
   k=NormalDist().inv_cdf(1-0.05/(2*6)) /
     NormalDist().inv_cdf(0.975).
   Give each IRR its approximate inflated range
   [p-*exp(-k*D), p+*exp(k*D)] and project the ratio boxes as above.
   A nonempty intersection is only compatibility under this **explicit
   approximate Bonferroni/Wald sensitivity**. The source does not supply
   log-normality, a six-IRR covariance matrix, or correctly adjusted
   99.17% marginal CIs, so this cannot be advertised as an actual
   simultaneous 95% test. Unknown dependence is not assumed independent.
7. Exact-count diagnostic for all six printed nonzero rows, with zero-bin
   missing-adjusted-IRR rule, as above. No extra crude-rate ratio.

Report whether the **primary verdict** survives no-rounding, full-unit
widening, all bin omissions and the approximate multiplicity-aware variant.
A point-only failure by itself does not falsify the primary uncertainty-box
rule. These sensitivity families were chosen before the computation and
not tuned to results.

## Verdict vocabulary, report and implementation gates

- **CONSISTENT**: complete positive primary input and nonempty primary
  intersection; only consistency with this necessary descriptive condition,
  not truth of factorisation.
- **FALSIFIED**: complete positive input and empty primary intersection
  under the *specified printed-interval descriptive* compatibility rule;
  report epsilon* and any sensitivity that restores compatibility.
- **UNTESTABLE**: primary adjustment/denominator inputs or printed paired
  estimates missing/invalid/nonpositive; no ratio is invented; zero bins
  outside the panel are separately untestable.
- **INCONCLUSIVE**: inferential status without joint covariance and a
  calibrated simultaneous procedure, or when a robustness/inference
  condition cannot establish a directional statistical claim. This
  designation may coexist with a descriptive CONSISTENT or FALSIFIED
  primary label, and is not a backdoor acceptance of the null.

The final report must give exact pair/bound/intersection decisions,
limiting bins, epsilon/log witness, all sensitivities, exact-count
diagnostics, literature screen, CC BY attribution, access limitations,
source checks, a hostile counterargument and specific residual gaps.
State prominently: **THIS DOES NOT TEST WHETHER HPV VACCINATION WORKS OR
WANES; both are already studied/established in the cited literature.**
Programme, birth cohort, prior HPV exposure, background time, unequal
age follow-up and screening are inseparable here from any biological
age-at-dose interaction. Neither passing nor failure proves mechanism.

Implementation only with Python stdlib and offline immutable sparse CSV.
Validate all row keys, provenance, precision and source facts. Store
SHA-256 of PROTOCOL.md, data CSV and analysis.py with CRLF/CR-to-LF
normalised UTF-8; keys are repository-relative Path.as_posix(). Hashes
document exact artifacts, not source truth. Provide independent tests
for formula/endpoints and attaining witness, bad input/zero handling,
precise bounds, multiple sensitivity branches, no-network execution,
Windows path/newline portability and CLI report checks. Python 3.12 and
3.14; <60 seconds; do not add dependencies or edit existing tests.
Hostile review before push; one PR into main; wait for actual CI results.
