# Frozen protocol: incidence shape within published strata

**Nonclinical, unreviewed computational research. No medical advice.**

Protocol date: **2026-11-18, America/New_York**. Branch from current `main`,
commit `88c7e07668887f6e941a3fe7592d53002199669a`. Commit this file alone
before computing any new empirical monotonicity, calibration, age-deletion,
or model-grid result. Preserve this file unchanged thereafter; disclose any
deviation in the README. The final README will identify the protocol commit.
The published All Races breast result in PR #64 and the source tables are
already known. This is a prospective computational specification for public
data, not blinding, independent registration, or independent replication.
Source availability, headers, rights and transcription may be checked before
the freeze; numerical outcome analysis may not.

## Question and predecessors

Run ranked job 4 in [the research programme](../README.md): does PR #64's
rounded-mean shape failure survive the published race strata, and which ages
force relaxation? Read the programme's negative results, [PR #64's note and
code](../multistage-incidence-identifiability/), [PR #70's note and protocol]
(../stem-cell-divisions-risk/), [DATA_POLICY.md](../../DATA_POLICY.md) and
[AGENTS.md](../../AGENTS.md) before analysis. Reuse #64's Theorems 7 and 8,
including their observation assumptions; do not claim new mathematics or
biological driver identification. #70's reproducibility and interpretation
discipline applies, but its tissue-division data are not inputs here.

## Fixed data, strata and exclusions

Use only the NCI **SEER Cancer Statistics Review 1975-2017**, November 2019
submission, released April 15, 2020; **Age-Specific SEER Incidence Rates,
2013-2017**, SEER 21 areas, rates per 100,000. The fixed selections are:

| Role | Site | Table | Printed race columns | Printed sex columns |
|---|---|---|---|---|
| Primary | Invasive female breast | 4.11 | All Races, Whites, Blacks | Females |
| Prespecified contrast | Colon and rectum | 6.10 | All Races, Whites, Blacks | Both Sexes, Males, Females |
| Prespecified contrast | Lung and bronchus | 15.10 | All Races, Whites, Blacks | Both Sexes, Males, Females |
| Prespecified contrast | Pancreas | 22.7 | All Races, Whites, Blacks | Both Sexes, Males, Females |

Use eleven ordered five-year bins: 30-34 through 80-84, operationally
`[30,35),...,[80,85)`. This fixes 30 curves and 330 potential rate cells.
Do not add ages, sites, races or sexes based on outcomes. Exclude all
age-adjusted and mortality sections, younger ages, and the open-ended 85+
bin. Keep source race labels; they are published population categories,
not demonstrated homogeneous biological populations. All Races overlaps the
race-specific columns, and Both Sexes overlaps the sex-specific columns;
do not count these as independent replications or infer mixture weights.

Table URLs share `https://seer.cancer.gov/archive/csr/1975_2017/`:

| Table | PDF suffix | CSV suffix |
|---|---|---|
| 4.11 | `results_single/sect_04_table.11.pdf` | `download_csr_datafile.php/sect_04_table.11.csv` |
| 6.10 | `results_single/sect_06_table.10.pdf` | `download_csr_datafile.php/sect_06_table.10.csv` |
| 15.10 | `results_single/sect_15_table.10.pdf` | `download_csr_datafile.php/sect_15_table.10.csv` |
| 22.7 | `results_single/sect_22_table.07.pdf` | `download_csr_datafile.php/sect_22_table.07.csv` |

Prefer the official CSVs if publicly readable. Otherwise transcribe PDF text,
cross-checking headings and age rows against public HTML/PDF. Never use a
browser file download; abandon any download prompt and record the limitation.
Preserve the printed rate strings. Store only these small selected aggregate
cells, not full source documents. Every row must carry site, sex, race, age,
period, units, display precision, source/version, source URLs, table, section,
column, PDF page, access date, reuse basis and reuse URL. The archive's
Suggested Citation declares the report public domain. Source content is not
relicensed under the repository's MIT license.

Keep case counts and person-years blank: these selected tables do not print
the needed cells. Suppression is missing, never zero. If any required age is
missing, suppressed, nonnumeric or has a nonpositive lower rounding endpoint,
mark the entire affected curve untestable; report the reason and do not fill
it from another release. Duplicates or inconsistent metadata are errors.

## Primary estimand and decision rule

For printed rate `y_j`, use exact decimal/rational values and the closed
rounding interval `l_j=y_j-0.05`, `u_j=y_j+0.05`. This is the same assumed
half-last-unit convention as #64, not a verified rounding algorithm, sampling
confidence interval or measurement-error distribution. Boundary equality is
compatible. Also repeat the monotonicity calculation on exact printed points
(`l_j=u_j=y_j`) as a prespecified sensitivity.

The homogeneous observation model has constant Poisson initiation of
independent identically distributed clones, optionally with an independent
clone-level diagnosis lag. Its direct hazard is a constant times a completion
CDF and is nondecreasing. Nonnegative normalized age weights on ordered,
nonoverlapping bins preserve the order. This necessary restriction does not
assume uniform weights, a particular stage count, or positive net growth.
Published period rates need not satisfy this observation model.

For each eligible curve compute exactly:

1. Every ordered pair `i<j` with `l_i>u_j`, its age labels and minimum
   decline `l_i-u_j` per 100,000. No decline is selected by significance.
2. `R_mono=max(1,max_{i<=j} l_i/u_j)` and
   `rho_mono=sqrt(R_mono)`. Record the exact rational square and all
   maximizing earlier/later witnesses when `R_mono>1`. If `R_mono=1`, say
   no relaxation is forced; do not manufacture a limiting-age pair.
3. The attaining nondecreasing envelope can be checked as
   `z_j=max_{i<=j} l_i/rho`. This is a sharp minimum for **all positive
   nondecreasing mean vectors**, only a necessary lower bound for a specific
   branching model. Allowing multipliers `a_j in [1/rho,rho]` with free
   amplitude quantifies relative shape relaxation, not measured reporting
   bias or an attributable cause.

**Frozen primary verdict:** require all three complete breast curves. If the
All Races curve fails and at least one of Whites/Blacks fails, report
**shape failure persists within strata**, explicitly distinguishing both
strata from only one. If All Races fails but both race-specific curves pass,
report **shape failure disappears within the prespecified strata** (only this
necessary rounded-mean condition). If required inputs are unavailable, or the
All Races benchmark does not reproduce its failure, report **untestable** for
the stated persistence question and explain the discrepancy. Contrasts do not
overrule the primary verdict. Passing is not proof of homogeneous biology.

## Which ages force relaxation: frozen sensitivity checks

For every curve, retain all outcomes, including negative ones:

* Compute the cumulative-prefix monotonicity frontier after every successive
  age bin, and identify the first ending age that forces `rho_mono>1`.
* Remove each age bin in turn and recompute the frontier on the remaining
  ordered bins. Retain all eleven leave-one-age-out results.
* Enumerate minimum-cardinality sets of age bins whose deletion restores
  rounded-interval monotone feasibility; retain every tied minimum set.
  These are diagnostic counterfactual deletions, never the primary panel.
* Compare the full-age outcome with the prespecified prefixes ending at
  65-69, 70-74 and 75-79. Explain whether removing only 80-84 suffices.

These checks localize a constraint; they do not show that an age group is a
cause, an error, or an appropriate clinical exclusion. Do not report
between-stratum significance tests, weighted pooled estimates, or unknown
case counts. Qualitative/rank comparisons of thresholds are descriptive.

## Exact #64 model calibration frontier

Repeat #64's Theorem 7 within each eligible site/race/sex curve. For a fixed
positive bin shape `f`, compute
`L=max_j(l_j/f_j)`, `U=min_j(u_j/f_j)`,
`R_f=max(1,L/U)`, `rho_f=sqrt(R_f)`.
The amplitude interval at any allowance is `[L/rho,rho*U]`. Record all
constraining ages attaining `L` and `U`; these are shape-normalized witnesses
and need not be the monotonicity witnesses or in chronological order.

Freeze the same two limiting families and grid as #64:

* Pure power `f_k(t)=t^(k-1)/(k-1)!`, `k=2,...,10`.
* One net-growth pole
  `f_{k,g}(t)=sum_{n>=0} g^n*t^(n+k-1)/(n+k-1)!`, with the same k values
  and `g=-0.060,-0.055,...,0.060` per year (spacing 0.005, including zero).
  Negative g is net loss. These are the rare-terminal-marking/high-immigration
  limiting shapes, not exact finite-rate fits to an observed biological model.

Use uniform bin means
`fbar_{k,g}=sum_{n>=0} g^n*(b^(n+k)-a^(n+k))/((b-a)*(n+k)!)`.
Use #64's exact rational powers and 70-term-index Taylor enclosure (indices
0 through 70, with its explicit absolute geometric-tail bound) so every grid
candidate and optimum is enclosed. If numerical intervals cannot decide a
winner or constraining age, report all unresolved contenders instead of
claiming unique selection. Compare before decimal rounding. Verify displayed
six-decimal thresholds from both bounds. No continuous or full-MSCE optimum
is claimed. The grid's k is a shape index, never an identified driver count.

Report for every k its pure-power threshold, best grid threshold/g,
constraining ages and uniqueness status. Report overall finite-grid minima
and retain ties. As #64's within-bin-weight sensitivity, also compute the
pure-power frontier allowing each bin shape throughout the endpoint closure
`[a^(k-1)/(k-1)!, b^(k-1)/(k-1)!]`. This replaces L by
`max(l/f_upper)` and U by `min(u/f_lower)` and gives a sharp closure
**infimum**, not guaranteed attainment on a half-open bin. It changes the
shape-family threshold but cannot remove a true ordered-bin decrease.

## Reproduction, falsification and reporting commitments

Confine changes to this directory and the single new file
`tests/test_research_incidence_strata_test.py`. Use standard-library Python
only; no network during analysis/tests. Reuse #64's numerical primitives
without modifying its files, and cross-check the decisive arithmetic with
independent implementations. The new test must run in under 60 seconds,
recompute every reported numerical result from the CSVs, check all generated
README result tables, and check protocol/data hashes. Hash UTF-8 decoded
content after CRLF/CR-to-LF normalization. Build every path key with
`Path.as_posix()` (including any portable-path test equivalent); test LF,
CRLF, CR-only and Windows/POSIX key behavior explicitly.

Before pushing completed work, audit every selected CSV rate against its
source, preferably using the official CSV and PDF text as two representations
of one source. Record failed accesses, discrepancies and remaining limits.
Hostile self-review must examine source headings/columns, age selection,
rounding/equality boundaries, nonadjacent violations, missing data, limiting
ages, interval-vs-sampling distinctions, grid-vs-global optima, portability,
scope and novelty. Run repository validation and relevant CI without changing
existing files. Push the branch and open a PR; report actual CI status.

The falsification ledger must retain whether rounding, race stratification,
sex/site contrasts, terminal-age truncation, individual age deletion,
minimum deletion sets, within-bin weights and the one-growth grid weaken or
remove the corresponding failure. Arbitrary age multipliers can fit every
positive finite shape by construction; this is nonidentification, not evidence
for a mechanism. Without counts/exposures, no p-values or statistical
rejection are justified. Period/cohort composition, observation, detection,
selection and heterogeneity remain possible explanations, not identified
causes. No medical advice, diagnosis or treatment recommendation is supplied.

Novelty is limited to a reproducible, prespecified extension of #64 to these
published strata and its age-constraint audit. The formulas, broad
nonidentifiability warning, and the published incidence patterns are prior
work. No new theorem, first discovery of an age decline, causal mechanism,
biological driver estimate, independent replication or global-priority claim
will be made.
