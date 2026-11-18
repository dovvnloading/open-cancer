# Persister-colony resistance-origin reanalysis

**Nonclinical computational research. Unreviewed. No treatment advice.**

## Completion status and pre-fit deviation record

**The requested raw Mendeley V2 reanalysis is blocked and incomplete.** The
public dataset's licence and advertised archive checksum were inspected, but
the archive could not be downloaded in this environment. Its bytes have not
been verified, its internal raw-data locators have not been inspected, and no
raw population trajectory has been recalibrated. A successful offline test or
CI run must not be interpreted as completing those missing scientific steps.

The protocol was frozen in local commit `cb76f03` before any empirical fitting. Its
normalized SHA-256 is
`a2a2b44d24977626c5cb395e8832e45378d0db2f9e898c9db23a2cda5b6b35e2`.
The protocol itself remains unchanged. The pre-fit deviation was recorded
in local commit `431881a`. Direct Git push lacked CLI credentials, so publication
uses the authenticated GitHub Git Data API, which creates new server commit
identifiers. [Original local Git objects](data/protocol_freeze.json) preserve
the actual local commit payloads, their tree membership, and the pre-fit
deviation text. The test checks their Git object hashes and protocol content.
These are local chronology records, not an externally trusted timestamp or
independent preregistration. Do not confuse server publication time with the
earlier local freeze.

**Deviation D1, recorded before fitting:** perform a separately labeled,
partial **published-count analysis** using the exact annotated Fig. 3e values
for DiFi B6 and WiDr B7. These are real experimental aggregate observations,
not synthetic data, but they are **not a fetched raw-deposit dataset**. They
are literal printed integers, not pixel-height estimates. The fetched
supplementary workbook supplies context and published model-dependent
sensitivity anchors only. Do not substitute its fitted parameters for raw
population measurements or treat Table 5's inferred rates as count outcomes.
The additional clones cannot enter a raw-data fit without their source data.

The planned whole-biological-replicate scoring, competing observation laws,
priors, sensitivity choices, and verdict rule remain those in
[PROTOCOL.md](PROTOCOL.md). This deviation is a restricted deliverable while
the archive is inaccessible; it does not assert that the user's requested
raw-data job is finished. The eventual PR is to remain draft for that reason.

The acquisition failures were an HTTP 403 / Cloudflare 1010 response at the
documented Mendeley public download and API endpoints, and timeouts without a
saved-file result from supported browser download operations. No access
control was bypassed. [fetch_sources.py](fetch_sources.py) accepts an official
locally downloaded archive and checks its entire advertised digest before
listing candidate data members. A source-access success must precede the
remaining extraction, raw-data cross-check, and population calibration work.

## Pre-fit interpretation constraints

The allowed families can have exactly the same observable occupancy laws:
late visibility may arise from an already competent but delayed lineage, or
from competence acquired after exposure. A fit of that common law cannot
identify the latent allocation. The relevant exact obstruction is the
vector-readout extension of PR #57, Theorem 5.3, **“Induced clonal
amplification”**, with the singleton special case in Theorem 5.1,
**“Different unknown clone laws”**. This is a structural statement conditional
on the broad kernel family, not empirical proof of either mechanism.

The independent-holiday positive identification theorem in PR #57 is not
applicable to repeated readouts during continuous drug exposure. Drug-caused
excess acquisition is a further causal question; Theorem 5.6, **“Drug causation
needs an exposure control”**, remains relevant even if event timing becomes
identifiable.

## Verdict and its practical limit

**Not identifiable from the published count target under the frozen broad
kernel family. The requested raw-data audit remains incomplete.** The
no-acquisition and acquisition families admit the same complete observable
well-occupancy distribution. The latent late-origin fraction can run from
entirely delayed pre-existing competence to entirely acquired competence while
the distribution stays fixed. This is an identified-set statement inside the
specified model, not a confidence interval, a demonstration of absent
acquisition, or a finding that a pre-existing delayed reservoir really existed.

The primary predictive model also has a real weakness: both WiDr late-count
observations fall outside their held-out marginal intervals. The more
dispersed, prespecified sensitivity models broaden those intervals. They were
not selected after seeing held-out performance. With so few biological
experiments, the data cannot establish reliable between-experiment dispersion
or interval calibration. Therefore this analysis does **not** call the
no-acquisition biology validated. A structural origin alias and an imperfect
observation model can coexist.

All computed values, including adverse folds and every sensitivity scenario,
appear in the generated report below. The test regenerates that report from
the committed data. Exact annotated counts and exact selected workbook cells
were checked against the downloaded source image and XLSX. The raw-archive
and population-measurement gaps are separate failures and remain visible.

## What was fitted

The observation is whether a well has become visibly positive, not a count of
all mutated cells or all mutation events. In each biological experiment, the
three categories are early positive, newly late positive, and neither. Late
wells are at risk only if they were not already early positive. Clone identity
is preserved; the biological-replicate labels are local left-to-right labels
for the independent experiments in the published figure, not recovered raw
well identifiers. The actual reporting intervals are retained in the CSV.
Their difference between cell lines prevents a simple pooled timing analysis.
[Figure and methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC9279152/#Fig3)

The observation model uses a Jeffreys Dirichlet prior and an exact
Dirichlet-multinomial posterior predictive. Each fold holds out an entire
biological replicate, including both its early and late observations. Only
the other experiment from the same clone trains that fold. The score is the
log probability of the full held-out category-count vector, including the
multinomial coefficient. The fitted distribution is integrated over parameter
uncertainty; it is not a plug-in fit to the held-out values.

The primary model assumes independent wells conditional on common clone-level
category probabilities. That is not independently established here. Exact
replicate-specific endpoints are unavailable within the figure's reporting
intervals, so different observation durations could contribute to the
apparent biological variation. The concentration-cap sensitivities are
deliberately broader predictive distributions with the same training mean.
They are **not fitted estimates** of biological variation. Prior sensitivity
and dispersion sensitivity are reported separately. Neither source-fitted
growth parameters nor previously inferred mutation rates enter the count
likelihood.

The **prediction target is held out**, but the numerical results cannot be
called a blinded external validation: this is a retrospective reanalysis of
a published figure, with the access-related source deviation recorded before
fitting. No model is tuned on the held-out outcomes. Reuse of each experiment
in the opposite fold means the aggregate scores are not independent
replications of a mechanistic experiment.

## Why origin remains aliased

Let `a` be the effective early detectable-event intensity and let `b_D` and
`b_I` denote delayed pre-existing and post-exposure acquired late-event
intensities. With `b=b_D+b_I`, the model gives

```
p_E = 1 - exp(-a)
p_L = exp(-a) * (1 - exp(-b))
p_U = exp(-a-b).
```

Setting `b_D=(1-theta)b` and `b_I=theta*b` for any admissible `theta` changes
the origin allocation while preserving all category probabilities. It
therefore preserves the full multinomial law and the integrated held-out
predictive law when the same prior is placed on observable probabilities.
A Bayes factor driven only by different priors on unidentified origin labels
would not supply data identification.

This equality also has a joint-readout construction before occupancy
coarsening. For a post-exposure event intensity `eta(t)P(t)` and a possibly
invisible, delayed vector-mark kernel `K_I(z;t)`, define

```
B = integral eta(t) P(t) [1-K_I(0;t)] dt
K_D*(z) = integral eta(t) P(t) [K_I(z;t)-K_I(0;t)] dt / B.
```

For positive `B`, `K_D*` is a valid positive-mark probability generating
function. Direct substitution gives

```
integral eta(t) P(t) [K_I(z;t)-1] dt = B [K_D*(z)-1].
```

The right side can represent lineages already competent at exposure with an
unknown delayed outgrowth/readout kernel. Early and reversible components
can be kept identical. This is precisely the stated vector-readout extension
of **“Induced clonal amplification”**; the singleton reallocation is
**“Different unknown clone laws”**. The tests check an exact rational example
with unequal kernels and invisible events as well as the fitted occupancy
allocations. See the [PR #57 note](../induced-vs-selected-resistance/README.md),
Theorem 5, for the original scalar constructions and their assumptions.

The extension requires the admissible pre-existing delayed-kernel family to
contain the constructed kernel. That is an intentionally broad assumption.
This work neither measures such a kernel nor asserts that all constrained
biological no-acquisition models would fit. A separately established upper
bound on latency, absence of competent founders, or an origin-specific readout
calibration could exclude the construction. Those constraints are not
established by the published occupancy annotations.

The finite-stock issue must not be hidden behind Poisson notation. A literal
fixed cell population does not have an unbounded Poisson number of founders.
At the occupancy level, the endpoint alternatives also admit a finite
construction: assign at most one potential late-visible competent lineage to
each early-negative well with the model's conditional late-positive
probability. It can be competent at exposure and delayed, or first become
competent subsequently, while keeping its visibility distribution the same.
This is an observation-law construction, **not a reconstruction of the actual
pre-treatment expansion process**. No rare-founder stock was measured here,
so the biological capacity check in the requested raw job is still open.

## Turnover, unequal kernels, and causation

All birth/death, reversible-state, extinction, and emergence-delay effects
enter the model's event intensities and readout kernels. The common-law
construction does not require the pre-existing, reversible, and acquisition
kernels to be universally equal; it requires an overlapping observational
class. Positive clone size is largely lost by a binary occupied-well readout.
At a fixed event-level probability of a nonzero observed mark, changing the
positive clone-size kernel leaves occupancy unchanged. The analysis therefore
does not multiply a detection probability by an unconstrained clone-size
scale. Per-cell plating efficiency would require a further clone-size
thinning calculation and is not measured here.

The raw population assays would be the separate calibration target. Their
inaccessibility prevents the planned late-phase slope refit. The published
Table 4 effective persister-loss summaries are used only as conditional
sensitivity anchors. The displayed birth/death pairs have the same net loss;
they are not measurements of gross turnover. The source's death parameter
also aggregates losses such as back-switching to a drug-sensitive state, so
its biological interpretation is more restricted than a directly observed
death count. [Supplementary material](https://doi.org/10.1038/s41588-022-01105-z)

For any positive finite calibrated detectable-persister exposure
`X=integral P(t)q(t)dt`, choosing `eta=b_I/X` realizes a specified positive
acquisition intensity. A changed turnover/readout assumption can change `X`
and the corresponding `eta` while preserving `b_I`. This is an algebraic
parameterization, not a new mutation-rate estimate or independent predictive
validation from dynamics. The alias survives these positive-exposure
normalizations; an empirically justified zero exposure or incompatible
latency bound would matter. Such constraints must be measured, not selected
to force an origin conclusion.

These findings do not establish that treatment caused no mutations, that
resistance was pre-existing, that persister states are not drug induced, or
that the paper's molecular interventions are invalid. The initial persister
fraction `f_0` is also a different quantity from the initial frequency of
stable resistant competence. Even perfectly identifying acquisition **after**
exposure would not by itself identify the drug-caused excess over background
acquisition. That needs a calibrated exposure counterfactual, as explained
by PR #57's **“Drug causation needs an exposure control.”** No patient outcome,
dosing policy, or treatment decision follows from this work.

## Hostile review and falsification record

| Attempt | Result and consequence |
|---|---|
| Verify reuse rights before committing third-party data | Mendeley and the article identify CC BY 4.0. Attribution, licence links, and transformations are recorded in [data/README.md](data/README.md). No raw archive is mirrored. |
| Acquire and check the requested raw deposit | Failed in this environment. The advertised archive hash is not represented as locally verified. The PR remains draft and the raw-data job remains incomplete. |
| Independently check source transcription | Original figure annotations were visually inspected, and all selected workbook cells were re-extracted from hash-pinned XLSX bytes. Local subplot labels were corrected to avoid implying recovered raw IDs. |
| Check denominator, pairing, and unit conversions | Late counts use the early-negative risk set. Entire biological replicates are held out. Week conversion and the different rate-table powers of ten were audited. |
| Check original source uncertainty labels | The workbook/caption interval discrepancy and inferred Table 5 uncertainty labels are preserved rather than silently corrected. |
| Falsify a no-late null | A null that permits neither delayed pre-existing colonies nor acquisition assigns zero probability to positive late counts. It is falsified by the published observations. It is not the full no-acquisition family or a faithful reimplementation of Russo's model. |
| Challenge predictive adequacy | Primary WiDr late intervals miss in both folds. Broader prespecified dispersion sensitivities improve coverage, but neither biological independence nor calibrated uncertainty is established. |
| Challenge identifiable timing | The nonnegative zero/positive-acquisition transformations preserve the entire observation law. More replication of the same endpoint design cannot distinguish these exact aliases. |
| Challenge turnover interpretation | Net loss is compatible with different births, deaths, and state-loss mechanisms. Published anchors are kept separate from the colony likelihood; raw recalibration is missing. |
| Challenge the unequal-kernel implementation | An initial implementation incorrectly scaled occupancy intensity by positive clone size. Internal review caught it; the corrected conversion uses the probability of a nonzero event mark. A regression test guards this correction. |
| Challenge founder capacity and latency | No rare-founder assay or measured maximum delay was available. The construction is explicitly conditional; no existence claim is made for a delayed resistant reservoir. |
| Challenge independence of verification | AI review, source hashes, and passing tests are internal computational checks. They do not constitute external scientific review, biological validation, or source-truth certification. |

## Single sharpest follow-up experiment

Use a **prospective, lineage-linked sister-culture experiment** with continuous
single-cell imaging and randomized drug/control exposure. Before exposure,
barcode and split matched sister lineages, archive baseline material, and
measure stable resistance competence and relevant resistance-associated
genotypes with a documented detection limit. Follow the corresponding
lineages through dormancy, divisions, death, first stable competence, and
visible colony emergence. Independently measure the lag from competence to
the same colony-readout threshold.

The decisive observation is a timestamped transition within the lineage that
later produces the colony, linked to its pretreatment state and its measured
outgrowth delay. Randomized matched exposure supplies the counterfactual for
a drug-caused change in acquisition hazard. Baseline sister sampling is
imperfect: report its sampling/detection error and do not equate a negative
baseline screen with certain absence. Barcodes alone, a late endpoint
genotype, or additional endpoint replicate counts do not timestamp
acquisition and would leave the central alias largely intact. This is a
nonclinical experiment proposal, not an executable clinical recommendation.

## What is new relative to Russo et al.

Russo et al. already developed the two-step fluctuation assay, estimated
population dynamics, incorporated birth/death and stochastic extinction,
modeled drug-effect delay, and simulated early/late emergence under their
model. They also examined resistant-cell fitness variation and provided
molecular/perturbation evidence beyond these endpoint counts. Their
chronological resistance estimates are effective phenotypic conversion
rates, not direct genome-wide mutation measurements. None of those
contributions is new here. [Primary article and supplementary note](https://pmc.ncbi.nlm.nih.gov/articles/PMC9279152/)

The relevant difference is the observation-family audit: permit pre-existing
competent lineages with broader unknown outgrowth delays and unequal
origin-specific kernels, exhibit the complete-law alias explicitly, and
retain whole biological experiments for predictive scoring. Their timing
simulations assume that resistant cells acquire untreated-like growth/death
behavior; successful separation inside that family does not rule out every
broader delayed-readout family. The supplementary rate-ratio posterior also
uses a prior support that excludes zero acquisition, so that posterior alone
is not a test of the broad zero-acquisition hypothesis. These facts do not
refute their independent biological evidence.

**Novelty assessment: limited and preliminary.** This is an application of
already-known compound-Poisson/observation-law reasoning and the repository's
existing alias theorem, plus a small held-out audit of published experimental
counts. It is not the first analysis of these experiments, a new biological
discovery, a complete replication of Russo's estimators, or a completed first
raw-public-data test. No claim of bibliographic priority is made. The new
predictive numbers document the chosen observation model's behavior; the
structural alias was already implied by allowing the broad nuisance family.

## Reproduce and finish the blocked portion

From the repository root:

```sh
python research/persister-colony-reanalysis/analysis.py --format markdown
python research/persister-colony-reanalysis/analysis.py --format json
python -m unittest discover -s tests -p 'test_research_persister_colony_reanalysis.py' -v
```

The test is standard-library-only and offline. It checks the frozen protocol
hash, normalized text digests, portable POSIX path keys, complete regeneration
of the numerical report, independent finite probability identities, source
unit conversions, held-out leakage, and exact aliases. Source-byte checks and
worksheet re-extraction run separately using the commands in
[data/README.md](data/README.md). The archived-source fetch helper never
executes downloaded scientific code.

Finishing the originally requested work requires obtaining the exact official
archive, matching its advertised whole-file hash, mapping actual raw-file
locators, checking culture/replicate and missingness information, comparing
raw counts with the published values, and independently calibrating population
measurements. Only then can the raw-data completion status change. Do not
silently relabel these published annotations or the published model outputs
as raw data. Any changed outcome interpretation or model choices require an
explicit deviation record before refitting.

## Numerical report

The block below is generated verbatim by analysis.py and checked by the
offline test. Its source-context tables retain the original provenance and
uncertainty labels; they are not additional observations for the count fit.

<!-- BEGIN COMPUTED REPORT -->
## Recomputed quantitative report

**Nonclinical, unreviewed; no treatment advice.** These tables are generated by
`analysis.py`. Count probabilities concern well occupancy, not counted genetic
acquisitions.

**Verdict: not identifiable from published count target; requested raw audit incomplete.**

Source status: `partial_published_count_reanalysis_raw_archive_blocked`. Raw archive checksum verified: **false**. Raw population calibration completed: **false**.

Target source labels: `published_figure_annotation_not_raw_deposit`.

### Colony observations

| Cell line / clone | Biological replicates | Wells | Early E | Newly late L | Neither U | Late risk n−E | E/n | L/(n−E) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| DiFi / B6 | 2 | 3840 | 44 | 667 | 3129 | 3796 | 0.011458 | 0.175711 |
| WiDr / B7 | 2 | 3840 | 14 | 667 | 3159 | 3826 | 0.003646 | 0.174334 |

### Held-out predictive scores

Each row holds out a complete biological replicate. E and L intervals are
marginal equal-tail 95% beta-binomial predictive intervals, not simultaneous
intervals or confidence intervals for acquisition. Every joint log score includes
the multinomial coefficient. M0 and M1 share this observable distribution, so
their scores agree algebraically. Caps are prespecified dispersion sensitivities.

#### primary

| Clone | Held out | n; E,L,U | Predicted E | E interval | E covered | Predicted L | L interval | L covered | Joint log score, M0 = M1 |
|---|---|---|---:|---|---|---:|---|---|---:|
| DiFi/B6 | 1 | 1920; 21,322,1577 | 23.482 | [12, 38] | yes | 345.230 | [300, 393] | yes | -7.387315 |
| DiFi/B6 | 2 | 1920; 23,345,1552 | 21.483 | [10, 36] | yes | 322.248 | [278, 369] | yes | -7.458765 |
| WiDr/B7 | 1 | 1920; 6,309,1605 | 8.493 | [2, 18] | yes | 358.220 | [312, 406] | no | -8.572231 |
| WiDr/B7 | 2 | 1920; 8,358,1554 | 6.495 | [1, 15] | yes | 309.258 | [266, 355] | no | -8.768275 |

Joint log-score sum: **-32.186587** for each model; difference **0.000000**. Descriptive marginal coverage: E **4/4**, L **2/4**.

#### cap_10

| Clone | Held out | n; E,L,U | Predicted E | E interval | E covered | Predicted L | L interval | L covered | Joint log score, M0 = M1 |
|---|---|---|---:|---|---|---:|---|---|---:|
| DiFi/B6 | 1 | 1920; 21,322,1577 | 23.482 | [0, 215] | yes | 345.230 | [40, 877] | yes | -11.835682 |
| DiFi/B6 | 2 | 1920; 23,345,1552 | 21.483 | [0, 203] | yes | 322.248 | [32, 846] | yes | -12.078346 |
| WiDr/B7 | 1 | 1920; 6,309,1605 | 8.493 | [0, 100] | yes | 358.220 | [44, 894] | yes | -11.414392 |
| WiDr/B7 | 2 | 1920; 8,358,1554 | 6.495 | [0, 76] | yes | 309.258 | [29, 827] | yes | -12.115409 |

Joint log-score sum: **-47.443828** for each model; difference **0.000000**. Descriptive marginal coverage: E **4/4**, L **4/4**.

#### cap_100

| Clone | Held out | n; E,L,U | Predicted E | E interval | E covered | Predicted L | L interval | L covered | Joint log score, M0 = M1 |
|---|---|---|---:|---|---|---:|---|---|---:|
| DiFi/B6 | 1 | 1920; 21,322,1577 | 23.482 | [1, 80] | yes | 345.230 | [210, 504] | yes | -9.206308 |
| DiFi/B6 | 2 | 1920; 23,345,1552 | 21.483 | [0, 76] | yes | 322.248 | [192, 477] | yes | -9.412300 |
| WiDr/B7 | 1 | 1920; 6,309,1605 | 8.493 | [0, 46] | yes | 358.220 | [221, 519] | yes | -8.723190 |
| WiDr/B7 | 2 | 1920; 8,358,1554 | 6.495 | [0, 40] | yes | 309.258 | [181, 462] | yes | -9.311847 |

Joint log-score sum: **-36.653646** for each model; difference **0.000000**. Descriptive marginal coverage: E **4/4**, L **4/4**.

#### cap_1000

| Clone | Held out | n; E,L,U | Predicted E | E interval | E covered | Predicted L | L interval | L covered | Joint log score, M0 = M1 |
|---|---|---|---:|---|---|---:|---|---|---:|
| DiFi/B6 | 1 | 1920; 21,322,1577 | 23.482 | [10, 42] | yes | 345.230 | [290, 403] | yes | -7.574318 |
| DiFi/B6 | 2 | 1920; 23,345,1552 | 21.483 | [8, 39] | yes | 322.248 | [269, 379] | yes | -7.689077 |
| WiDr/B7 | 1 | 1920; 6,309,1605 | 8.493 | [1, 20] | yes | 358.220 | [303, 417] | yes | -8.158163 |
| WiDr/B7 | 2 | 1920; 8,358,1554 | 6.495 | [0, 17] | yes | 309.258 | [257, 365] | yes | -8.460216 |

Joint log-score sum: **-31.881774** for each model; difference **0.000000**. Descriptive marginal coverage: E **4/4**, L **4/4**.

#### prior_0.25

| Clone | Held out | n; E,L,U | Predicted E | E interval | E covered | Predicted L | L interval | L covered | Joint log score, M0 = M1 |
|---|---|---|---:|---|---|---:|---|---|---:|
| DiFi/B6 | 1 | 1920; 21,322,1577 | 23.241 | [11, 38] | yes | 345.115 | [299, 393] | yes | -7.370546 |
| DiFi/B6 | 2 | 1920; 23,345,1552 | 21.242 | [10, 35] | yes | 322.124 | [278, 368] | yes | -7.478125 |
| WiDr/B7 | 1 | 1920; 6,309,1605 | 8.247 | [2, 17] | yes | 358.110 | [312, 406] | no | -8.526956 |
| WiDr/B7 | 2 | 1920; 8,358,1554 | 6.248 | [1, 14] | yes | 309.129 | [265, 355] | no | -8.824861 |

Joint log-score sum: **-32.200487** for each model; difference **0.000000**. Descriptive marginal coverage: E **4/4**, L **2/4**.

#### prior_1

| Clone | Held out | n; E,L,U | Predicted E | E interval | E covered | Predicted L | L interval | L covered | Joint log score, M0 = M1 |
|---|---|---|---:|---|---|---:|---|---|---:|
| DiFi/B6 | 1 | 1920; 21,322,1577 | 23.963 | [12, 39] | yes | 345.460 | [300, 393] | yes | -7.424607 |
| DiFi/B6 | 2 | 1920; 23,345,1552 | 21.966 | [11, 36] | yes | 322.496 | [278, 369] | yes | -7.424607 |
| WiDr/B7 | 1 | 1920; 6,309,1605 | 8.986 | [2, 19] | yes | 358.440 | [312, 407] | no | -8.672512 |
| WiDr/B7 | 2 | 1920; 8,358,1554 | 6.989 | [1, 16] | yes | 309.516 | [266, 355] | no | -8.672512 |

Joint log-score sum: **-32.194237** for each model; difference **0.000000**. Descriptive marginal coverage: E **4/4**, L **2/4**.

### Explicit origin aliases

The probabilities below are full-data posterior means under the primary prior.
The displayed intensities are deterministic transforms of those means, not
posterior mean mutation rates. All listed allocations leave the observable law
unchanged; numerical equality is checked separately from the analytic identity.

| Clone | Effective a | Effective b | p(E) | p(L) | p(U) |
|---|---:|---:|---:|---:|---:|
| DiFi/B6 | 0.011652 | 0.193338 | 0.011584 | 0.173760 | 0.814656 |
| WiDr/B7 | 0.003782 | 0.191667 | 0.003775 | 0.173760 | 0.822465 |

| Clone | theta | Delayed pre-existing b_D | Acquired b_I |
|---|---:|---:|---:|
| DiFi/B6 | 0.00 | 0.193338 | 0.000000 |
| DiFi/B6 | 0.25 | 0.145004 | 0.048335 |
| DiFi/B6 | 0.50 | 0.096669 | 0.096669 |
| DiFi/B6 | 0.75 | 0.048335 | 0.145004 |
| DiFi/B6 | 1.00 | 0.000000 | 0.193338 |
| WiDr/B7 | 0.00 | 0.191667 | 0.000000 |
| WiDr/B7 | 0.25 | 0.143751 | 0.047917 |
| WiDr/B7 | 0.50 | 0.095834 | 0.095834 |
| WiDr/B7 | 0.75 | 0.047917 | 0.143751 |
| WiDr/B7 | 1.00 | 0.000000 | 0.191667 |

The following conversion assumes the stated probability that a latent
event leaves a nonzero detected mark. It is an event-level probability,
not a measured per-cell plating efficiency. Changing the positive kernel's
relative size scale leaves well occupancy unchanged at fixed nonzero-mark
probability. These are sensitivity assumptions, not measured founder counts.

| Clone | Nonzero-event detection probability | Positive kernel relative scale | Latent late intensity per well |
|---|---:|---:|---:|
| DiFi/B6 | 0.10 | 0.25 | 1.933380 |
| DiFi/B6 | 0.10 | 1.00 | 1.933380 |
| DiFi/B6 | 0.10 | 4.00 | 1.933380 |
| DiFi/B6 | 0.50 | 0.25 | 0.386676 |
| DiFi/B6 | 0.50 | 1.00 | 0.386676 |
| DiFi/B6 | 0.50 | 4.00 | 0.386676 |
| DiFi/B6 | 1.00 | 0.25 | 0.193338 |
| DiFi/B6 | 1.00 | 1.00 | 0.193338 |
| DiFi/B6 | 1.00 | 4.00 | 0.193338 |
| WiDr/B7 | 0.10 | 0.25 | 1.916674 |
| WiDr/B7 | 0.10 | 1.00 | 1.916674 |
| WiDr/B7 | 0.10 | 4.00 | 1.916674 |
| WiDr/B7 | 0.50 | 0.25 | 0.383335 |
| WiDr/B7 | 0.50 | 1.00 | 0.383335 |
| WiDr/B7 | 0.50 | 4.00 | 0.383335 |
| WiDr/B7 | 1.00 | 0.25 | 0.191667 |
| WiDr/B7 | 1.00 | 1.00 | 0.191667 |
| WiDr/B7 | 1.00 | 4.00 | 0.191667 |

### Separate population calibration

No raw-data late-phase population slope could be calibrated.

Calibration gaps:

- no accessible raw population measurements supplied.

### Published parameter context and fallback sensitivity

All source estimates below come from the committed supplementary-table
extraction. They are **published model-dependent context**, never training
observations for the colony fit or a substitute for raw population calibration.
Untreated birth/death means and standard errors, TP-model prior limits,
persister-loss posterior summaries, and resistance-rate estimates retain
their distinct source roles and uncertainty labels.

| Supplementary table | Clone | Parameter | Source statistic | Value | Source units | Exact worksheet cell |
|---:|---|---|---|---:|---|---|
| 1 | DiFi/B6 | `birth_rate` | `mean` | 0.7 | `per_day` | `Supplementary table 1!D4` |
| 1 | DiFi/B6 | `birth_rate` | `standard_error` | 0.02 | `per_day` | `Supplementary table 1!E4` |
| 1 | DiFi/B6 | `death_rate` | `mean` | 0.26 | `per_day` | `Supplementary table 1!D5` |
| 1 | DiFi/B6 | `death_rate` | `standard_error` | 0.02 | `per_day` | `Supplementary table 1!E5` |
| 1 | WiDr/B7 | `birth_rate` | `mean` | 1.02 | `per_day` | `Supplementary table 1!B4` |
| 1 | WiDr/B7 | `birth_rate` | `standard_error` | 0.08 | `per_day` | `Supplementary table 1!C4` |
| 1 | WiDr/B7 | `death_rate` | `mean` | 0.34 | `per_day` | `Supplementary table 1!B5` |
| 1 | WiDr/B7 | `death_rate` | `standard_error` | 0.09 | `per_day` | `Supplementary table 1!C5` |
| 2 | DiFi/B6 | `D_p     (days^−^1)` | `prior_max` | 0.2 | `source_label:D_p     (days^−^1)` | `Supplementary table  2!F32` |
| 2 | DiFi/B6 | `D_p     (days^−^1)` | `prior_min` | 0 | `source_label:D_p     (days^−^1)` | `Supplementary table  2!D32` |
| 2 | DiFi/B6 | `f_0` | `prior_max` | 1 | `source_label:f_0` | `Supplementary table  2!F31` |
| 2 | DiFi/B6 | `f_0` | `prior_min` | 0 | `source_label:f_0` | `Supplementary table  2!D31` |
| 2 | WiDr/B7 | `D_p     (days^−^1)` | `prior_max` | 0.2 | `source_label:D_p     (days^−^1)` | `Supplementary table  2!F9` |
| 2 | WiDr/B7 | `D_p     (days^−^1)` | `prior_min` | 0 | `source_label:D_p     (days^−^1)` | `Supplementary table  2!D9` |
| 2 | WiDr/B7 | `f_0` | `prior_max` | 1 | `source_label:f_0` | `Supplementary table  2!F8` |
| 2 | WiDr/B7 | `f_0` | `prior_min` | 0 | `source_label:f_0` | `Supplementary table  2!D8` |
| 4 | DiFi/B6 | `D_p` | `hdi_lower_source_label_3_percent` | 0.034 | `source_label:D_p` | `Supplementary table  4!H36` |
| 4 | DiFi/B6 | `D_p` | `hdi_upper_source_label_97_percent` | 0.047 | `source_label:D_p` | `Supplementary table  4!J36` |
| 4 | DiFi/B6 | `D_p` | `posterior_mean` | 0.04 | `source_label:D_p` | `Supplementary table  4!E36` |
| 4 | DiFi/B6 | `D_p` | `posterior_sd` | 0.004 | `source_label:D_p` | `Supplementary table  4!G36` |
| 4 | WiDr/B7 | `D_p` | `hdi_lower_source_label_3_percent` | 0.053 | `source_label:D_p` | `Supplementary table  4!H9` |
| 4 | WiDr/B7 | `D_p` | `hdi_upper_source_label_97_percent` | 0.094 | `source_label:D_p` | `Supplementary table  4!J9` |
| 4 | WiDr/B7 | `D_p` | `posterior_mean` | 0.073 | `source_label:D_p` | `Supplementary table  4!E9` |
| 4 | WiDr/B7 | `D_p` | `posterior_sd` | 0.011 | `source_label:D_p` | `Supplementary table  4!F9` |
| 5 | DiFi/B6 | `persister_chronological_resistance_rate` | `estimate` | 3 | `1e-6_per_cell_per_day` | `Supplementary table 5!D5` |
| 5 | DiFi/B6 | `persister_chronological_resistance_rate` | `uncertainty_sd` | 0.3 | `1e-6_per_cell_per_day` | `Supplementary table 5!E5` |
| 5 | DiFi/B6 | `sensitive_chronological_resistance_rate` | `estimate` | 4.5 | `1e-7_per_cell_per_day` | `Supplementary table 5!D4` |
| 5 | DiFi/B6 | `sensitive_chronological_resistance_rate` | `uncertainty_sd` | 0.7 | `1e-7_per_cell_per_day` | `Supplementary table 5!E4` |
| 5 | WiDr/B7 | `persister_chronological_resistance_rate` | `estimate` | 7 | `1e-6_per_cell_per_day` | `Supplementary table 5!B5` |
| 5 | WiDr/B7 | `persister_chronological_resistance_rate` | `uncertainty_sd` | 0.4 | `1e-6_per_cell_per_day` | `Supplementary table 5!C5` |
| 5 | WiDr/B7 | `sensitive_chronological_resistance_rate` | `estimate` | 1.4 | `1e-7_per_cell_per_day` | `Supplementary table 5!B4` |
| 5 | WiDr/B7 | `sensitive_chronological_resistance_rate` | `uncertainty_sd` | 0.4 | `1e-7_per_cell_per_day` | `Supplementary table 5!C4` |

Source interpretation notes:

- Published model-derived estimate; not an independent raw observation and not used to train colony model.
- Uncertainty inferred from Supplementary Note eq13 for sensitive and lines497-499 for persister; workbook lacks subcolumn labels.
- Worksheet headers say HDI 3%/97%; caption says 2.5%/97.5%. Preserve discrepancy; do not silently call these 95% intervals.

#### Comparison with the authors' inferred rates

Unit multipliers are taken from each source row before calculating the
ratio. These are the authors' point estimates and reported uncertainties;
the ratio is descriptive and has no newly inferred uncertainty interval.

| Clone | Sensitive estimate per cell per day | Reported uncertainty | Persister estimate per cell per day | Reported uncertainty | Ratio of point estimates |
|---|---:|---:|---:|---:|---:|
| DiFi/B6 | 4.5e-07 | 7e-08 | 3e-06 | 3e-07 | 6.666667 |
| WiDr/B7 | 1.4e-07 | 4e-08 | 7e-06 | 4e-07 | 50.000000 |

#### Published net-loss anchors; no raw-data recalibration

The source TP-model persister-loss parameter D_p is used only as a
conditional net-loss anchor, g = −D_p. For every published mean and
each reported bound, the frozen birth-rate sensitivity is applied and
death = birth + D_p. These are assumptions consistent with the same net
loss; they do not identify actual persister division or death rates.
The bounds retain their source labels because the worksheet and caption
disagree about interval endpoints. These are not newly estimated intervals.

| Clone | Source D_p statistic | D_p per day | Assumed birth per day | Required death per day | Assumed net per day |
|---|---|---:|---:|---:|---:|
| DiFi/B6 | `hdi_lower_source_label_3_percent` | 0.034000 | 0.000000 | 0.034000 | -0.034000 |
| DiFi/B6 | `hdi_lower_source_label_3_percent` | 0.034000 | 0.500000 | 0.534000 | -0.034000 |
| DiFi/B6 | `hdi_lower_source_label_3_percent` | 0.034000 | 1.000000 | 1.034000 | -0.034000 |
| DiFi/B6 | `hdi_upper_source_label_97_percent` | 0.047000 | 0.000000 | 0.047000 | -0.047000 |
| DiFi/B6 | `hdi_upper_source_label_97_percent` | 0.047000 | 0.500000 | 0.547000 | -0.047000 |
| DiFi/B6 | `hdi_upper_source_label_97_percent` | 0.047000 | 1.000000 | 1.047000 | -0.047000 |
| DiFi/B6 | `posterior_mean` | 0.040000 | 0.000000 | 0.040000 | -0.040000 |
| DiFi/B6 | `posterior_mean` | 0.040000 | 0.500000 | 0.540000 | -0.040000 |
| DiFi/B6 | `posterior_mean` | 0.040000 | 1.000000 | 1.040000 | -0.040000 |
| WiDr/B7 | `hdi_lower_source_label_3_percent` | 0.053000 | 0.000000 | 0.053000 | -0.053000 |
| WiDr/B7 | `hdi_lower_source_label_3_percent` | 0.053000 | 0.500000 | 0.553000 | -0.053000 |
| WiDr/B7 | `hdi_lower_source_label_3_percent` | 0.053000 | 1.000000 | 1.053000 | -0.053000 |
| WiDr/B7 | `hdi_upper_source_label_97_percent` | 0.094000 | 0.000000 | 0.094000 | -0.094000 |
| WiDr/B7 | `hdi_upper_source_label_97_percent` | 0.094000 | 0.500000 | 0.594000 | -0.094000 |
| WiDr/B7 | `hdi_upper_source_label_97_percent` | 0.094000 | 1.000000 | 1.094000 | -0.094000 |
| WiDr/B7 | `posterior_mean` | 0.073000 | 0.000000 | 0.073000 | -0.073000 |
| WiDr/B7 | `posterior_mean` | 0.073000 | 0.500000 | 0.573000 | -0.073000 |
| WiDr/B7 | `posterior_mean` | 0.073000 | 1.000000 | 1.073000 | -0.073000 |

<!-- END COMPUTED REPORT -->
