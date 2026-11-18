# Analysis plan frozen before new fitting

**Nonclinical, unreviewed computational research; no treatment advice.**

Date: 2026-11-17 (America/New_York). Base main:
`4f2098f543eb36d2810809af600c4372673794e7`.

This file is committed before this project's numerical fitting or threshold
calculations. It is a prospective computational plan for a retrospective
evidence audit, not a registered or blinded study. The requested prior notes
were read first; their Lundberg rejection percentages, Shibata point fits,
and binomial threshold were already known. Initial source inspection has also
revealed Shibata's fitted maternal-total negative-binomial parameters and
aggregate model-comparison table. These are disclosed prior knowledge, not
held-out validation. Amendments will be dated in README, with their timing
relative to calculations; this file will not be silently revised.

## Inclusion and extraction

Search primary literature and references for ecDNA/double-minute partition
at mitosis, daughter/sister count pairs, clustering, retention, and live
imaging. Include contradictory results, perturbations, and historical work.
Record searched sources and unavailable/excluded measurements; a search
cannot guarantee every publication has been found. Peer-reviewed papers and
preprints must be distinguished. Use only already-published aggregates:
fitted parameters, summary moments, sample sizes, published test results,
or already-tabulated count bins. No plot digitization, individual-cell
records, patient records, controlled data, or inferred raw datasets.

Every numerical source row must identify population/condition, statistic,
units, source/DOI, exact locator, access date, and source-specific reuse
basis. Source materials remain external. Public accessibility alone does
not license redistribution. Copy only tiny permitted aggregate extracts;
source licenses are not replaced by the repository MIT license.

## Estimands and evidence hierarchy

For a mother with pre-replication copies k, the model doubles to n=2k and
draws P~Beta(p*(1/rho-1),(1-p)*(1/rho-1)), J|P~Bin(n,P).
Sisters are complementary. A uniformly selected sister has effective
R=(2p-1)^2+4p*(1-p)*rho, distinct from rho and from an interspecies
correlation coefficient. Report p, rho, and R separately. Nuclear loss,
copy-dependent clustering, count error, selection, and chromosomal controls
are potential violations, not additional estimates of the same parameter.

1. Prefer a likelihood from published aggregate bins or conditional moments
   with denominators/uncertainty. Fit each study/line separately first.
2. If only published maximum-likelihood parameters exist, retain them as
   source fits, not a new refit of inaccessible observations. Transform them
   to R. A source aggregate log-likelihood/AIC comparison may test model
   distinctions but cannot supply an unavailable empirical profile.
3. If source fits also describe the count distribution and the number of
   independent mother divisions is justified, calculate **conditional
   model-based uncertainty** using the source observation likelihood and
   expected Fisher information or parametric repeated-sampling simulation.
   Explicitly disclose reconstruction assumptions, zero floors, sorting,
   boundary parameters, nuisance parameters and normalization. This is not
   resampling observed cells, an exact CI, or a validation of the source fit.
   Do not synthesize observations and present them as published data.
4. If the input summaries do not support (3), report identification bounds,
   rounding sensitivity and scenario uncertainty separately. Do not invent
   standard errors or turn selected clone spread into sampling uncertainty.

Use per-study/per-line inference. Pool across studies only if at least two
studies supply compatible estimands and independent sampling units. Report
heterogeneity descriptively otherwise; no forced hierarchical meta-analysis.
Exclude HSR-only controls from primary ecDNA summaries while retaining them
in the audit; show any disputed clone separately and in sensitivity analysis.

## Baseline tests and uncertainty

Test unbiased binomial p=1/2,rho=0 separately from biased binomial rho=0.
Prefer published aggregate likelihood comparisons and calibrated reported
tests. For Lundberg's rejection counts, enumerate compatible integers from
printed percentage rounding, use exact binomial tails conditional on
calibrated size and independent mothers, Holm adjustment across untreated
lines, and an independence-free Markov sensitivity. Missing calibration or
culture allocation prevents an unconditional empirical rejection claim.
Do not treat paired sisters or thousands of ecDNAs as independent divisions.
Nominal uncertainty level is 95%; if a rectangle is used for a pair of
parameters, use simultaneous coverage (e.g. Bonferroni marginal intervals).
No unreported exact likelihood-ratio p value will be fabricated from MLEs.

## Threshold propagation and numerical verification

Use PR60's linear death law d_k=gamma*k, c=gamma/b and its infinite-type
positive-lineage threshold. The binomial reference is recomputed, not copied
as evidence. Preserve the full (p,rho) kernel: R alone does not determine a
threshold. Report absolute c shifts and relative percentage shifts.

For every retained source fit, compute a numerical threshold with finite
subprocess lower certificates and a valid all-type tail upper bound. If
propagating a parameter interval, prove a uniform envelope or other valid
bound; endpoint interpolation requires a proved monotonicity theorem.
Keep sampling/model uncertainty, source rounding, and numerical error
distinct. The bound c* <= (1-R)/2 is only an upper bound, not the threshold.
Any new proof must include its assumptions, critical/boundary cases, and a
verifier. Unproved general monotonicity remains a conjecture.

## Prediction, falsification, and deliverable checks

Select one sharp future observable using the strongest supported estimate:
prefer a conditional sister-count statistic that distinguishes baseline and
clustered models without requiring unmeasured pharmacological calibration.
Provide a concrete rejection rule with stated independent sampling unit.
An infinite-time threshold is not falsified by finite-time persistence.
Record failed hypotheses, opposing literature, source discrepancies,
nonidentifiability and novelty limitations. No biological discovery claim
will be substituted for evidence that was not located.

Only this new research directory and
`tests/test_research_ecdna_segregation_clustering_data.py` may be added.
The test uses the standard library, runs under 60 seconds, recomputes every
reported derived number from the checked-in CSV inputs, and checks provenance,
likelihood/kernel normalization, uncertainty distinctions and proof
certificates. Existing files, data policy, workflows and research catalogs
remain untouched. Hostile source and mathematical review precedes push.
