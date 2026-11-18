# Frozen protocol: MHC-I and cancer mortality transport

**Nonclinical, unreviewed computational research. No medical or veterinary advice.**

## Status and question

This protocol and the derived common-panel inputs are committed before any
model fitting or predictive scoring. Source definitions, licence terms, schema,
missingness, decimal precision, and tree coverage were inspected first. The
published direction of association was already known. This is a retrospective,
outcome-informed research question, not a prospectively collected validation
study. The freeze prevents changing the panel, models, scores, and decision
rule after seeing their predictive results; it does not erase that prior knowledge.

Base: `dovvnloading/open-cancer` main commit
`ff08fb20b9a3c723c11187b5f68d75b5bf07002b`. PR #77 was open and unmerged when
this work began, so the supplied brief is the controlling specification.

Question: does a single `log(MHC-I gene-family size + 1)` term improve
prediction of cancer-attributed deaths in an entirely held-out mammalian
taxonomic order, beyond body mass, adult life expectancy, and diet?
The existing MHC/cancer association is prior work. The proposed contribution
is an incremental out-of-order predictive transport audit, not a discovery
of that association.

## Sources, reuse, and observation unit

Vincze et al. (2025), *Immunological surveillance against cancer across mammals*,
[DOI](https://doi.org/10.1038/s41467-025-65286-x),
[readable article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12644480/), especially
Methods **Cancer mortality risk**, **Confounding effects**, and availability.
The paper-linked repository is
[VinczeEtalNatComm_MHC-cancer](https://github.com/OrsolyaVincze/VinczeEtalNatComm_MHC-cancer),
pinned at `2c8a5c5b5fda6fb34c86b0e5a2b6a58dcc85cc3e`.

The pinned repository has no LICENSE or README and GitHub reports no licence.
The article specifies [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/).
Neither fact is treated as permission to mirror the upstream assets. Commit
only original fetch/analysis code and narrow derived numerical species
aggregates with attribution. Do not commit the source CSV, source tree, R code,
article text, or restricted individual records. The project's MIT licence does
not relicense third-party material. The user explicitly authorized this
aggregate-only fallback. Access/reuse audit date: 2026-11-18 (America/New_York).

Read `SourceData_FINAL.csv` and `phylogeny.tre` as plain text through GitHub's
API/raw endpoint. A future optional fetch reads response text in memory and
verifies pinned content; it does not use a browser file download. No browser
download or prompt was used in preparation. Do not request restricted ZIMS data.

The observational unit is a species aggregate from managed populations.
CMR is cancer-attributed deaths divided by RDI, deaths with known cause.
Adult life expectancy is mean remaining life after sexual maturity, in days;
it is neither maximum longevity nor lifespan from birth. Diet is the source's
binary mammal-prey indicator. MHC is an inferred protein-coding gene-family
size, not a measured expressed or functional antigen repertoire. ICM is a
censoring/truncation-adjusted cumulative mortality estimate, not a simple
count proportion. **Never recover an integer count from ICM times RDI.**

Mortality shares the Species360/ZIMS provenance family of PR #75 (Vincze 2022);
the later extraction does not establish disjoint animals. **This is not
independent mortality validation.** Gene-family size does not identify immune
causation or measure antigen quality, expression, peptide binding, complete
presentation loss, or the functional load in PR #59's theorem. PRs #67, #72,
and #76 motivate demanding a predictive test without supplying those bridges.

## Frozen panel and count reconstruction

Use every source species with nonmissing finite CMR, positive integer RDI,
positive body mass and adult life expectancy, binary diet, and nonnegative
integer MHC-I family size. Match exact source species identifiers; no synonyms,
new sources, outlier screens, or outcome-dependent quality exclusions. Both
models always use precisely the same panel. The sorted identifiers and values
in the first committed `data/species_aggregates.csv` define that panel.

The source has 139 rows. This rule gives 125 species and 16 taxonomic orders,
including 17 recorded zero-CMR species, `Heterocephalus_glaber`, and
`Dolichotis_patagonum`. Twelve rows lack MHC and two further rows lack CMR
(`Bos_gaurus`, `Ovis_orientalis`); an ICM of zero cannot fill a missing CMR.
Record every excluded identifier and reason in `data/exclusions.csv`.
All primary identifiers occur in the pinned tree. No species is removed merely
to make the tree or a sensitivity convenient.

Before fitting, audit each literal decimal CMR token with exact decimal
arithmetic. For non-endpoint token c with decimal exponent e, the admissible
probability interval is c +/- (0.5 * 10^e + 2e-15), clipped to [0,1]. The fixed
2e-15 allowance covers floating-point export at long decimal precision; it is
not estimated from model fit. Admissible integer k satisfy k/n inside this
interval, inclusive of endpoints. Literal 0 and 1 are treated as exact endpoint
proportions, consistent with the source's interpretation of no detected cases;
the displayed integer zero does not mean rounding to the nearest whole probability.

Use nearest integer to n*c (decimal ROUND_HALF_EVEN) as the **nominal** numerator,
but report count_lower/count_upper and reject any row whose nominal numerator
is outside its admissible interval. Do not claim unique reconstruction where
the display precision permits several counts. Here `Eulemur_macaco` permits
2..6 deaths (nominal 4 of 40), and `Myrmecophaga_tridactyla` permits 1..3
(nominal 2 of 20). Exhaust all 15 assignments as a declared sensitivity. All
other panel numerators are unique under this rule. The nominal sum is 697
among 12033 known-cause deaths, with a compatible total range of 694..700.
The largest nearest-integer discrepancy is 1.7e-7 deaths, from a short decimal
for `Lama_guanicoe`; a blanket 1e-8 count tolerance would be inappropriate.
These summaries are a prefit data audit, not predictive results.

The article's Methods summary population counts do not agree with the pinned
file's panel. Record that discrepancy; do not adjust the file to match the text
or claim exact replication of the published panel.

## Paired learning rules

Let n be RDI and k the nominal recovered cancer deaths. Both models are ridge
binomial logistic regressions with an unpenalized intercept:

- Baseline: log(body mass in kg), log(adult life expectancy in days), diet.
- Extended: those same terms plus log(1 + MHC-I family size).

Use natural logarithms. In each training fold, center and scale each continuous
predictor using that fold's unweighted species mean and population standard
deviation. Leave binary diet as 0/1. A zero-variance continuous training column
gets scale 1. All held-out transformations use training parameters. No
interaction, polynomial, data-selected transformation, taxonomic dummy, or
selection of predictors is allowed.

The primary objective is the sum of species binomial negative log likelihoods
(without irrelevant combinatorial constants), plus (lambda/2) times the sum
of squared slopes, with **lambda = 1**. Equivalently the data part is
sum_i [n_i * softplus(eta_i) - k_i * eta_i]. Use a stable damped Newton method
with a convergence check. The half-event correction may initialize the
intercept but must not change the likelihood or labels. A numerical fit
failure is reported, never silently dropped. Predictions may be clipped to
[1e-12, 1-1e-12] only for numerical log scoring; report clipping if it occurs.
No hyperparameter is selected using held-out outcomes.

## Whole-order holdout and estimands

For each source taxonomic order, fit both rules on all other orders and
predict every species in the held-out order. This is leave-one-order-out
cross-validation. A test order contributes no outcomes to training, scaling,
penalty choice, or model selection. Its CMR/RDI are deterministically recovered
by the frozen count rule for scoring only; they do not inform training. Tree
weights are fixed outcome-blind source summaries, as specified below.

Species cancer-death log loss is

    L_i = -[k_i log(p_i) + (n_i-k_i) log(1-p_i)] / n_i.

This is Bernoulli cross-entropy per known-cause death, a proper score including
zero cases; when k=0 it is -log(1-p), not zero. Binomial combinatorial constants
cancel in paired comparisons. Average L_i across species within each test
order, then average the order means equally. This order-balanced loss is the
primary estimand. Define delta = baseline loss - extended loss, in natural-log
units per death; positive favors adding MHC. Report both losses and delta.

Also report the paired loss averaged over all known-cause deaths (RDI weights),
each order's species count and losses/delta, the number of orders with delta>0,
and all held-out species predictions. This separates transport to an equally
weighted observed order from performance dominated by large denominators.
No individual random train/test split is an alternative primary result.

## Order-level uncertainty and validity limits

Delete each order in turn and rerun the **entire** leave-one-order-out pipeline
on the remaining orders, refitting both training models and their transforms.
For G orders, delta on the full panel, and delta_(-g) after deleting order g,
define pseudovalue q_g = G*delta - (G-1)*delta_(-g).
Report the bias-corrected jackknife center q_bar and standard error

    SE = sqrt(sum_g (q_g - q_bar)^2 / [G*(G-1)]).

The approximate 95% interval is q_bar +/- t_(G-1,0.975)*SE. For the frozen
G=16, t_15,0.975 = 2.131449545559323. The primary point estimate remains the
unadjusted full-panel delta; label the interval's bias-corrected center
separately. Also show the full range of delta_(-g). Do not substitute an
individual-death or species bootstrap or a coefficient-significance test.

The jackknife captures order-related sensitivity of both fitting and scoring,
not just variation in fixed predictions. Its population coverage remains
approximate: the orders are few, uneven, selected, and related through deep
phylogeny. It cannot repair ascertainment, cause-of-death misclassification,
sampling selection, or biological proxy mismatch. Do not claim independent
or exhaustive uncertainty from reproducible arithmetic.

## Fixed falsification and sensitivity analyses

Every sensitivity retains the primary holdout/scoring procedure and paired
baseline comparison unless explicitly stated. Report every result; never
choose the favorable one after seeing results. Sensitivities are descriptive
robustness checks, not separately powered confirmations; they receive no
additional positive claim from an unadjusted significance threshold.

1. **Rounding:** rerun both models for all 15 count assignments defined above;
   report all deltas, their range and sign count. The nominal assignment is
   part of this grid. No ICM count substitutions.
2. **Named rodents:** remove naked mole-rat and Patagonian mara together and
   rerun; this is a sensitivity only, never the primary panel.
3. **Recorded zeros:** remove nominal k=0 species and rerun as a deliberately
   outcome-selected, nontransportable sensitivity to the published exclusion
   convention. It cannot rescue the primary claim.
4. **Training weighting:** replace n-weighted likelihood with equally weighted
   species cross-entropies multiplied by the training mean n. Thus effective
   n_i equals the training mean n for every species and k_i/n_i is preserved;
   keep lambda=1. This changes relative ascertainment weights while preserving
   total training likelihood scale and the primary test score.
5. **Penalty:** independently rerun with lambda=0.1 and lambda=10; do not tune.
6. **Phylogeny:** parse the pinned Nexus translation table and branch lengths.
   For each primary species, sum branch length divided by the number of
   primary-panel descendants along its root-to-tip path (fair-proportion
   evolutionary distinctiveness). Store this one derived number per species;
   do not publish the upstream tree or a reconstructable pairwise matrix.
   Multiply n_i and k_i in each training fold by this positive weight normalized
   to mean 1 within that training fold. Keep the primary loss and original
   k,n at scoring. This tests sensitivity to evolutionary representation in
   training. It is **phylogenetic reweighting, not a phylogenetic covariance
   adjustment**, phylogenetic logistic mixed model, or proof of independence.
   Fixed weights use the primary tree panel even in a reduced sensitivity;
   no outcomes enter their calculation. Full phylogenetic covariance modeling
   and uncertainty across alternative trees remain untested. Ignore the root
   stem (zero in this tree); reject nonfinite or negative branch lengths.
7. **Deeper phylogenetic holdout:** set a cut depth equal to half the median
   root-to-tip distance of the primary species, without inspecting outcomes.
   Assign each tip to the first edge crossing this depth. Merge cut groups
   whenever one taxonomic order spans them, so an order cannot straddle folds.
   Sort components by their smallest species identifier to label them. Hold
   out entire resulting components, fit the same two rules, and average losses
   first over species within components, then equally across components. This
   changes the transport target to deeper tree blocks and is not directly the
   primary order-balanced estimand. With at least 4 components and valid fits,
   its delta must also be positive for a positive verdict; with fewer than 4,
  report the low-information descriptive result without using it as a gate.
   Do not retune the cut to obtain a convenient number of blocks. Store only
   the resulting per-species component ID and tree-height audit, not the tree.
8. **Order deletion:** report the uncertainty refits above as falsification
   attempts for dependence on a single influential order.

## Frozen verdict rule

First fail to **untestable** if the panel/count gates cannot be reproduced,
there are fewer than 8 complete orders, or any primary/uncertainty fit fails.
Otherwise apply the following exhaustive rules, without changing thresholds:

A failed required major-sensitivity fit prevents a positive verdict and is
reported as untestable for that robustness requirement. The deeper-tree
minimum-component gate depends on component count alone; invalid fits cannot
turn a component-count-eligible check into an optional one.

- **MHC-I adds portable predictive information** within this evaluation only
  if full-panel delta>0, the jackknife interval lower bound>0, more than half
  the orders improve, every rounding assignment has delta>0, and every major
  sensitivity (rodent removal, equal-species weighting, both penalties,
  phylogenetic reweighting, and deeper tree-block holdout when its gate is met)
  has delta>0. The outcome-selected zero-removal
  result cannot satisfy or veto this criterion.
- **MHC-I does not add portable predictive information** within this evaluation
  only if full-panel delta<=0 and the jackknife interval upper bound<=0.
  This is conditional evidence of non-improvement for the specified learning
  rules and panel, not proof that immunity is irrelevant.
- **Untestable (inconclusive on the available panel)** otherwise. This includes
  intervals spanning zero and failed robustness for an apparently positive
  primary result. It means that the requested directional claim is unresolved,
  even though the predictive experiment was computationally performed.

No positive result supports causation, an independent mortality validation,
functional-repertoire measurement, or clinical/veterinary action. A negative
or inconclusive result does not falsify immune surveillance as a mechanism.

## Reproduction, review, and immutable scope

Implement with Python standard library only. A single new test file under
`tests/` must run offline in under 60 seconds and recompute all reported
analysis numbers from committed derived inputs. Compare the saved structured
results and generated README result block, check count reconstruction and
all species/order holds, and exercise zero loss, held-out-label isolation,
numerical convergence, and the decision rule. Hash UTF-8 text after replacing
CRLF and CR by LF, without trimming text. Construct every file-path dictionary
key using `Path.as_posix()` so hashes and reports agree on Windows and Linux.

Only add this new research directory and
`tests/test_research_mhc_cancer_mortality_transport.py`; do not modify
`research/README.md` or any existing file. Perform hostile review of data,
leakage, uncertainty, novelty, and claims before pushing. Record unresolved
limitations and all protocol deviations. A protocol correction after fitting
must be an explicit dated amendment, never an undisclosed rewrite. Commit
history must retain this prefit protocol separately from fitted results.

Open the requested PR and inspect its actual CI; successful tests establish
computational consistency, not source truth, source reuse rights, or biological
validity. Do not merge the PR as part of this task.
