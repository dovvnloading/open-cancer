# Frozen protocol: resistance origin from persister-colony observations

**Nonclinical, retrospective computational research. Unreviewed. No treatment
advice.** Protocol version: 2026-11-18. This file is to be committed before
running any fit, scoring calculation, or data-driven model selection. Its
commit must remain an ancestor of the analysis commit. Corrections after that
point must be recorded as deviations in the README, without rewriting this
file. This is a prospective computational freeze within a retrospective
reanalysis, not a registered or outcome-blinded experiment: source inspection
includes the paper's displayed outcomes and the source-file structure.

## Question, target, and scope

Do the DiFi and WiDr experiments in Russo et al. (2022), DOI
[10.1038/s41588-022-01105-z](https://doi.org/10.1038/s41588-022-01105-z), require
stable resistance acquisition after exposure once delayed outgrowth, persister
turnover, and different clone-size/detection kernels are allowed?

The primary estimand is the fraction `theta` of **effective detectable late
clone-event intensity** attributable to stable competence first acquired after
exposure. It is neither a genome-wide mutation rate, a per-division mutation
rate, nor a drug-caused excess. Occupancy (whether a well contains a visible
colony) also does not identify how many clone events occurred in that well.
Drug causation would require a suitable exposure counterfactual, even if event
timing were known.

The user explicitly requested this experiment. At the starting main revision,
`44c7689e6e11458a341dcc0708f5df2ef521c414`, research/README.md instead ranks the
containment/aggressive-treatment control question first. That index will not be
edited. All new tracked files must be inside this directory, except the single
test `tests/test_research_persister_colony_reanalysis.py`. No patient records,
genomes, treatment recommendations, or changes to existing files are permitted.

## Sources, rights, and extraction

Use the public Mendeley **V2** deposit, DOI
[10.17632/mvfm7hs9kw.2](https://doi.org/10.17632/mvfm7hs9kw.2),
`Russo_et_al_v2.zip`, plus the article's Fig. 3e, Supplementary Tables 1, 2, 4
and 5, and the supplementary methods needed to interpret those observations.
The deposit page currently states **CC BY 4.0**. Confirm that attribution and
redistribution terms apply to the extracted files; record the licence link,
authors, version, access date, original-file locators, and changes made.
Retain the repository policy's distinction between public access and reuse.

Download from the official public endpoint, verify the publisher's checksum
where supplied, and record independent file digests. Do not commit the large
archive or images. Commit only the small cell-line experimental aggregates
needed for this analysis, with `source`, `source_file`, sheet/column or other
precise `locator`, and `access_date` on every data row. These are cell-line
cultures, not patient-level observations. If redistribution is disallowed,
commit a fetch-and-verify script and permissible derived summaries instead.
If a source cannot be accessed, record that gap without inventing raw data.

Prefer machine-readable original counts over digitization. Check counts
against Fig. 3e and source analysis code. Do not substitute a simulated result,
published fitted mutation rate, or a figure-derived guess for a raw count.
Population measurements and published parameter tables are calibration or
context only; they must never be appended to the colony target as extra
independent observations. Table 5's inferred rates are a comparison with the
authors, not observations used to fit the new model.

## Observation units and exclusions

Retain the clone and biological-replicate identifiers supplied by the source.
The colony prediction target is the complete early/late observation from an
independent biological experiment. Where the protocol counts wells, define:

* `n`: wells initially assayed;
* `E`: wells with an early visible resistant colony;
* `L`: newly positive late wells among those not already early positive;
* `U = n - E - L`: no recorded resistant colony at the final readout.

The late risk set is `n-E`, not `n`. Preserve actual observation times and
censoring; late visibility is not an acquisition time. Use the source's
clone-specific last readouts rather than silently imposing one common time.
If the archive contradicts this occupancy interpretation, stop the planned
multinomial fit and document a protocol deviation before an appropriate
count-law fit. If raw well IDs are available, aggregate them to this specified
target; do not count technical wells as independent biological experiments.

Exclude patient material, genomic records, unrelated cell lines, molecular
intervention experiments, drug-holiday results not paired with this assay,
simulations, duplicated exports, and already inferred mutation rates from the
primary target. Retain zeros. Exclude a missing/non-numeric outcome only with
an explicit reason; never impute a zero. Do not remove an entire replicate
because it fits poorly. Identify when raw culture identities, death counts,
plating efficiency, or exact first-emergence times are absent.

## Competing families and the observation map

Both families allow separate pre-existing, reversible, and acquired clone
kernels; delayed outgrowth and missed detections; and persister birth/death
turnover. Neither family equates an early/late visibility label with genotype
or resistance-acquisition timing.

* **M0: no post-exposure stable acquisition.** Stable competent clones may
  pre-exist but become visible early or late, or remain invisible. Reversible
  states and turnover may change their readout. The effective acquisition
  intensity is zero.
* **M1: post-exposure stable acquisition allowed.** The same components are
  present, plus events with intensity `eta(t) P(t) dt`, each having its own
  delayed joint readout kernel. M0 is its zero-acquisition boundary.

For endpoint vector arguments `z`, an acquisition component has log-PGF
`integral eta(t) P(t) [K_I(z;t)-1] dt`. Write
`q(t)=1-K_I(0;t)` and `B=integral eta(t) P(t) q(t) dt`. For `B>0`, define the
positive delayed vector-mark kernel

```
K_D*(z) = integral eta(t) P(t) [K_I(z;t)-K_I(0;t)] dt / B.
```

Then its log-PGF is exactly `B[K_D*(z)-1]`. The initial delayed competent pool
and the acquired pool have the same entire observed law if the allowed
pre-existing kernel includes `K_D*`. Keeping the other components fixed
preserves the alias. This is a vector-readout extension of PR #57, Theorem 5.3,
**“Induced clonal amplification”**. For singleton contributions also check
Theorem 5.1, **“Different unknown clone laws”**, equation (21). Invisible-event
and efficiency obstructions in 5.4 and 5.5 must not be overlooked. These
continuous-treatment readouts are not the growth-controlled `0,H,2H` holiday
design required by PR #57's positive reconstruction theorems.

The reduced independent-well observation law is multinomial with
probabilities `(p_E,p_L,p_U)`. An operational Poisson-event representation is

```
a = -log(1-p_E)
b = -log(p_U/(1-p_E))
p_E = 1-exp(-a)
p_L = exp(-a) * (1-exp(-b))
p_U = exp(-a-b)
b_D = (1-theta)*b; b_I = theta*b; 0 <= theta <= 1.
```

Thus timing-origin allocation may change while every observable probability
is fixed. The pre-existing and acquisition kernels need not be equal; their
admissible classes need only overlap at the exhibited alias. This is a
constructive nonidentifiability test within a deliberately broad family, not
a claim that every possible biological M0 is equivalent to every M1. Poisson
event pools are an observation model/rare-event approximation, not an exact
fixed-population genetic branching process. Check the necessary founder stock
against measured starting populations, and report how many delayed competent
events would be needed. Availability of a small stock does not prove that it
biologically exists.

## Fitting, held-out scores, and uncertainty

Fit each clone separately. Hold out each **whole biological replicate** in
turn, including both early and late counts. All other biological replicates
of that clone form its training data. Never split early from late or split
technical wells between training and test. If fewer than two independent
replicates are accessible, report that predictive validation is unavailable
for that clone instead of creating pseudo-replicates.

Primary count fit: a multinomial likelihood and exchangeable Jeffreys prior
`Dirichlet(1/2,1/2,1/2)`. Use the exact Dirichlet-multinomial posterior
predictive probability of the held-out complete count vector as the joint
log score, including its multinomial coefficient. Predictive distributions
are integrated over training-parameter uncertainty, not plug-in likelihoods.
Calculate marginal equal-tail 95% beta-binomial predictive intervals for E
and L. Report every fold's joint score, observed counts, predicted means,
intervals, and coverage, as well as their sum. Score M0 and M1 with the exact
same observable prior under their alias; a prior-dependent Bayes factor on
unidentified latent labels must not be presented as data identification.

An experiment-level random-environment sensitivity caps posterior
concentration at `10`, `100`, or `1000`, preserving the posterior mean before
computing the Dirichlet-multinomial predictive. The uncapped posterior is the
primary model. These caps deliberately explore extra biological variation;
they are not fitted estimates of dispersion. Report all, without choosing a
cap using held-out outcomes. Repeat the primary fit with symmetric prior
pseudocounts `1/4` and `1`. With so few biological replicates, interval
calibration and between-experiment variance cannot be established. Marginal
coverage is descriptive and is not simultaneous coverage. Do not report
well-level standard errors as biological-replicate certainty.

Also report full-data descriptive totals and the distribution of early/late
counts. An apparent M0/M1 score tie must be identified as algebraic equality
of the allowed observation laws, not a low-powered significance test. The
identified set for `theta` is obtained by explicit nonnegative alias
constructions, not by a falsely precise posterior interval.

## Separate dynamics calibration and prespecified sensitivity region

Use raw population measurements from independent growth/single-dose or
dose-response assays. Preserve biological replicate and time labels. Normalize
only within the source assay with a documented baseline. Treat viability/ATP
as a population proxy unless calibrated counts are supplied. Never identify
the complete residual population with known genotyped resistant founders.

For a transparent late-phase calibration, estimate each clone's net log-slope
using positive population observations at days `>=14` when at least two
distinct times exist. Calculate per-biological-replicate slopes where those
identities exist; otherwise report the source aggregation and no invented
replicate interval. Also calculate a `>=10`-day window as sensitivity. If raw
sampling does not support either window, retain that explicit calibration
gap and use the source's Table 4 net-loss interval only as a published
model-dependent sensitivity anchor. Do not treat its posterior as new raw
evidence or claim a re-fit of the full TP model.

Net growth `g` does not determine gross turnover. For nonpositive calibrated
`g`, allow birth `beta in [0,1]` per day and death `delta=beta-g`; for positive
`g`, allow death in that range and birth `beta=delta+g`. Evaluate the endpoints
and midpoints of this interval, as well as the zero-turnover boundary when
admissible. This region is a deliberately broad sensitivity choice, not a
measurement of either rate. Allow detection efficiencies `1/10,1/2,1`, unequal
kernel scales `1/4,1,4`, and arbitrary nonnegative delayed kernels consistent
with the observed early/late bins, including latent mass beyond the final
readout. The broad alias may hold at any calibrated positive persister
exposure; show its nonnegative scaling and do not estimate a mutation-per-
division rate from net decay alone.

Display the latent intensity allocation for `theta=0,1/4,1/2,3/4,1` and verify
that all predictive observables remain unchanged. Consider a restricted
no-delay/no-late-pre-existing model as a falsification control: if its allowed
late probability is zero, any positive late observation falsifies **that
restriction**, not the broad M0. Quantify the delayed pre-existing intensity
needed by M0 and compare to a baseline-population upper capacity; a fit is
not biological validation of that latency distribution. Separately discuss
what externally measured delay bounds, founder absence, birth/death, and
plating efficiency would break the constructed alias.

## Verdict rule, falsification, and interpretation

* **Acquisition required:** only if all admissible zero-acquisition members
  are ruled out by the actual observations/calibration and held-out
  acquisition predictions succeed across the specified nuisance region.
  A fit failure of one restricted null or a favorable likelihood under the
  authors' timing assumptions is insufficient.
* **Not required:** if an admissible zero-acquisition model predicts the
  observations, while the data/design otherwise support identification of
  the origin parameter. This never means acquisition is absent.
* **Not identifiable from these data:** if a nonnegative zero/positive
  acquisition alias survives the actual observed design and calibration.
  Report that verdict even if both observational predictors fit poorly;
  separately report model inadequacy. Do not force an unsupported causal
  yes/no answer.

Attempt to falsify provenance, transcription, endpoint/risk-set semantics,
population-to-count calibration, founder capacity, shared-kernel assumptions,
replicate exchangeability, and predictive coverage. Preserve adverse checks.
Name the relevant PR #57 alias exactly and state which added observations
would exclude it. Identify one sharp follow-up using prospective lineage
tracking plus measured emergence/detection, with randomized exposure controls
if drug causation is the question. Additional endpoint replication alone does
not resolve a structural alias.

Credit Russo et al.'s existing two-step assay, growth/death modeling, timing
simulations, molecular evidence, and rate estimation. Assess novelty as a
limited robustness/identifiability reanalysis unless a broader literature
review supports a stronger claim. This analysis cannot establish absence of
drug-induced mutation, overturn independent molecular experiments, infer
patient effects, or recommend a treatment.

## Reproduction and review gates

Use Python standard library for the committed analysis and offline regression
test. The test must run in less than 60 seconds, reproduce every reported
computed number from committed aggregates, check protocol/data digests on
UTF-8 content after CRLF/CR normalization to LF, and compare relative paths
with `Path.as_posix()`. No network calls in the test. Original binary-source
hashes remain hashes of bytes and must be distinguished from normalized text
digests. A fetch/extract verifier may run separately from the offline test.

Before pushing, perform a hostile source/method/code review and downgrade any
unproven statement. Run the repository's existing validation/test gates
without editing their files. Push the branch, open a PR, and inspect its
actual CI status. Passing checks establish reproducibility properties, not
source truth, scientific priority, or clinical validity.
