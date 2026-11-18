# Frozen analysis specification: CH growth and progression

Status: nonclinical, unreviewed mathematical research; no medical advice.

Frozen on 2026-11-17 (America/New_York), 2026-11-18 04:33 UTC,
before extraction of per-gene or VAF-stratified outcome tables. This is an
in-session prospective analysis specification, **not** a registered or fully
blinded study. The literature survey already exposed Fabre's growth-risk
association and headline CHRS, Genovese, and newer-study summaries. These
exposures must be disclosed; they do not constitute independent validation.

## Question and pivot

The association of clone growth with malignant progression is already known.
We ask the narrower question: can *shared extra-hit opportunity*, with measured
growth as the sole gene-varying biological parameter at fixed current clone
size and age, predict the magnitude of subsequent AML or myeloid-neoplasm
risk? Which parts are identifiable using published aggregates alone?

## Quantitative prediction, before table extraction

For a branching clone with n current competent stem cells, net exponential
growth s per year, and common rare successful-progression rate q per cell-year,
the leading cumulative hazard at horizon T is

    H_g(T) = n q A(s_g,T),  A(s,T) = integral_0^T exp(s t) dt.

Thus A=expm1(sT)/s, with A(0,T)=T. The fixed primary horizon is T=10 years;
T=5 and T=20 are sensitivity analyses, never selected on their agreement.
At matched n, age, and q, H_g/H_h=A(s_g,T)/A(s_h,T). Instantaneous hazards
instead have ratio exp((s_g-s_h)t) in the rare-mark limit. A shared nonnegative
latency kernel preserves the bound

    1 <= H_high/H_low <= exp((s_high-s_low)T).

At fixed s and all other parameters, cumulative hazard is proportional to n,
and hence VAF only under the stated diploid, heterozygous, unbiased-output
mapping n=2*N*VAF. Finite branching corrections and assumptions must be proved,
not presumed. Absolute individual risks require additional unobserved rates.

## Eligible aggregates and source hierarchy

1. Fitness: published numerical gene-level growth estimates and intervals,
   preferring Fabre 2022 Supplementary Table 6; retain the published definition
   (including logistic versus exponential growth). Robertson 2022 and Watson
   2020 are prespecified alternative sources if numerical gene-level estimates
   and uncertainty with compatible definitions and redistribution rights are
   accessible. Do not digitize plots or invent missing intervals. Hotspot,
   gene, median, mean, and mechanistic fitness definitions are not interchangeable.
2. AML outcomes: Abelson 2018 / the published aggregates underlying Fabre 2022
   Figure 5a first, Desai 2018 as a separate outcome series. MN: Weeks 2023 as
   a separate series. Include every eligible matched gene in each selected
   table, without selection on sign, magnitude, significance, or fit.
3. VAF-stratified published effect estimates, including their actual bin
   definitions, reference groups, endpoint, adjustment, confidence intervals,
   and cohort denominators. Do not assign arbitrary VAF bin midpoints.
4. Add newer eligible aggregate series located by the literature survey;
   record each addition and its rationale. Do not pool AML with broad MN,
   odds ratios with hazard ratios, or different reference populations.
5. Only small aggregate extracts whose reuse terms allow redistribution.
   Otherwise retain a link and an explicit exclusion reason. No individual
   trajectories, genomes, controlled data, plot-level participant points, or
   source workbooks containing participant records may be committed.

## Analysis and uncertainty

The mechanistic primary target is the conditional cumulative-hazard ratio.
Assess whether each published estimand can identify that target **before**
calling any fit a test. If not, the mechanistic verdict is untestable.

For sufficiently matched conditional cumulative-hazard ratios, compare the
fixed formula directly. For available gene-level OR or Cox-HR summaries,
compute a deliberately labeled, descriptive *transport stress test*, separately
for each endpoint/source: y_g=log(reported ratio_g), x_g=log A(s_g,T),
and y_g=alpha+x_g with slope fixed to one. Estimate alpha by the unweighted
mean of y-x in the training genes. Leave each gene out, estimate alpha on all
other genes, and predict the held-out ratio. Also fit the intercept-only null
on the same training sets. Report log-RMSE, median absolute log error, fold
residuals, and the complete per-gene table. No fitted horizon, fitted exponent,
gene exclusion, or within-sample correlation will substitute for prediction.

Use the reported lower and upper limits to form deterministic uncertainty
envelopes: training alpha ranges from mean(log lower_outcome-log A(upper_s))
to mean(log upper_outcome-log A(lower_s)); propagate the held-out fitness
limits monotonically. These are rectangular sensitivity ranges, **not** joint
95% prediction intervals: unknown covariance, shared controls, multiple
testing, and different studies prevent that interpretation. Use all endpoints
of each interval, even if inconvenient. An exact zero/infinite estimate or
missing limits is reported but excluded from finite-log calculations with a
reason; never add pseudocounts. No p-values assuming independent genes.

Rank deviations by absolute held-out log residual, show their direction, and
discuss the smallest added mechanism (a gene-specific conversion/latency
factor, hidden advanced subclones, turnover, or context) without assigning
causality from these aggregates. If observed and predicted envelopes overlap,
call the *stress test* consistent at that resolution; otherwise falsified
under its explicit transport assumptions. Always report separately whether
the biological model is testable. Failure to outperform the intercept-only
null is a negative predictive result, regardless of interval overlap.

## Required falsification attempts

Check exact versus rare-event formulas; net fitness versus turnover;
logistic versus exponential growth; heterogeneous current size and age;
pre-existing additional hits; phenotype-specific risk; nonidentical reference
groups; uncertain/co-occurring genes; different growth horizons; competition,
death and diagnostic latency. Attempt a counterexample to every necessity or
sufficiency claim. Clearly label unproved claims as conjectures and synthetic
computational checks as non-biological evidence.

## Completion and scope

All changes stay under research/clonal-hematopoiesis-progression/ and the one
new tests/test_research_clonal_hematopoiesis_progression.py. Preserve this
specification unchanged after its initial commit. Every source-derived number
reported in generated results must be recomputable from the small CSVs.
All executable research and its test use only Python's standard library;
the dedicated test must run in under 60 seconds. Normalize CRLF and CR to LF
before hashing file content. Record deviations from this plan in README.md.
