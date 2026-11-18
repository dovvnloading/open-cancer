# Frozen protocol: species-by-tissue mutation-clock separability

**Nonclinical, unreviewed computational research. No medical advice.**

Protocol date: **2026-11-18, America/New_York**. Branch from current `main`,
commit `09d916cdea7053f25643c6d0bd9c45827436e315`. Commit this file alone
before executing any empirical arithmetic, fitting, scoring, or bound
calculation. Preserve it unchanged thereafter and cite the freeze commit in
the README. Reading source definitions, licenses and printed values is allowed
before the freeze. The public data and the expected tissue differences are
already known; this is a pre-analysis computational specification, not blinded
research, independent registration, or a claim of prospective data collection.

## 1. Question, prior work and scope

Run ranked job 4 of `research/README.md`: quantify the restriction needed to
transport the colon passenger clock across tissues in the bridge discussed in
PR #75 (`research/petos-paradox-mutation-rates/`). Read that note, especially
`data/PROVENANCE.md` and the Ferret row, together with `AGENTS.md` and
`DATA_POLICY.md` before analysis.

The target restriction for positive annual somatic substitution rates is

`r(species, tissue) = species_factor * tissue_factor`.

The mouse/human rate ratio must consequently be identical in blood and colon.
This is a restriction on measured passenger-clock proxies under explicitly
stated measurement assumptions. It does not equate SBS/SNV counts to functional
driver acquisition, rate-limiting oncogenic transitions, or cancer incidence.
Tissue differences are prior work, explicitly discussed by Kapadia and the
specified neuron preprint. The algebra of multiplicative interactions is also
standard. The intended contribution is a reproducible quantitative audit of
this particular bridge, with an honest novelty assessment and no priority
claim.

The peer-reviewed blood/colon comparison is primary. The three-species
colon/neuron comparison is a separately labelled **preprint transport negative
control**. Postmitotic neurons are not cancer targets in this analysis. Neither
panel supports general species-level validation, causal attribution, or a
clinical conclusion.

## 2. Source selection and transcription

Use only small, printed aggregate rates obtained as **article HTML text or
GitHub raw/API text**. No browser file downloads, PDFs, supplementary binary
files, plot digitization, individual records, raw genomes, or controlled-access
data. Deny any download prompt, proceed without it, and record the limitation.
An inaccessible endpoint does not authorize a bypass. Indexed article HTML is
an eligible reading route, but record when direct full-text access failed.

1. **Kapadia et al. (2025), Nature**, DOI
   `10.1038/s41586-025-08625-8`, Results **Somatic mutation accumulation**,
   annual genome-wide SBS estimate and its printed confidence interval for
   cultured mouse HSC/MPP colonies. Preserve the separately printed per-base
   normalization. The primary genotype/strain and culture context must remain
   explicit. In **Fitness landscape of murine clonal haematopoiesis**, retain
   the printed targeted-duplex model's global per-base rate and its
   phylogenetic comparator only for the labelled assay sensitivity below.
2. **Mitchell et al. (2022), Nature**, DOI
   `10.1038/s41586-022-04786-y`, Results **Mutation burden and telomere
   lengths**, precise human HSC/MPP substitutions/cell/year estimate and
   confidence interval. Prefer the precise Results estimate over the rounded
   abstract value. Record that choice; do not count the rounded value as an
   independent observation.
3. **Cagan et al. (2022), Nature**, DOI
   `10.1038/s41586-022-04618-z`,
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC9021023/>, Results
   **Correlation with life-history traits**, opening human/mouse annual
   substitution comparison. These rounded mean per-genome rates are the
   primary colon inputs. Do not borrow regression-line intervals, intervals
   for other species, or intervals for a different rate estimator. No unseen
   supplementary values or graph coordinates may replace them.
4. **Ferret colon:** only the `Ferret` row, columns `Mean mutation rate
   (SBS/genome/year)` and, for an audit-only distinction, `Zero-intercept
   regression slope (MLE)` and its missing interval, in
   <https://github.com/zacharycompton/cancerAcrossVertebrates/blob/3b4db10c7c38c63aeebfc9a48ce4bfa376dfd2bc/Mutation_Data.csv>.
   Use the mean for prediction. Preserve the exact source taxon and the
   missing uncertainty; inspect the source directly as well as PR #75's
   transcription. Do not relicense Compton's other data or extract its cancer
   outcomes. Document the Cagan/Compton preprint reuse route separately.
5. **Caglayan et al., September 2026 preprint, version 1**, 
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC13622217/>, Results accompanying
   **Figs. 1f and 2i-l**, printed human, mouse and ferret glutamatergic-neuron
   annual SNV rates. Preserve approximation marks and printed SEMs. Record
   version, species/strain, genome normalization and assay. Human numeric
   uncertainty absent from the text stays missing even if a figure has a
   shaded interval. Retain exact ferret taxon labels across sources and make
   the species-level mapping an explicit assumption, not a taxonomic proof.

Every row must contain a stable row ID, source citation/URL, version, exact
section/paragraph/figure or row/column locator, access date, original units,
estimator, assay, genotype/strain, uncertainty type, missingness explanation,
reuse basis and terms URL. Record original decimal strings; empty fields mean
missing, never zero. Retain only factual aggregate measurements with a
documented permitted reuse basis; do not infer a CC license from free access.
If necessary inputs or reuse cannot be established through permitted routes,
report that comparison as untestable rather than inventing a value or license.
No source article, figure or full table is redistributed. Repository MIT
licensing does not relicense third-party material.

## 3. Primary statistic and sharp departure bound

Let `B_M, B_H, C_M, C_H` be the positive mouse/human blood/colon rates.
Keep the orientation fixed:

`R_B = B_M/B_H`, `R_C = C_M/C_H`,
`Q = R_B/R_C = B_M*C_H/(B_H*C_M)`.

Separability requires `Q=1`. Report the signed natural-log interaction
`log(Q)` and all three distinctly defined multiplicative quantities:

* `D = exp(abs(log(Q))) = max(Q, 1/Q)`: fold disagreement of the two
  species ratios.
* `delta = D**(1/2)`: the minimum symmetric allowance if **each of the two
  tissue-specific species ratios** may differ from one common ratio by a
  factor in `[1/delta, delta]`.
* `rho = D**(1/4)`: the minimum symmetric **per species-tissue cell**
  allowance in `r_st = a_s*b_t*u_st`, where every `u_st` is in
  `[1/rho, rho]`. This is the principal minimum multiplicative departure bound.

Proof to include: the log interaction is a signed sum of four log residuals,
so its magnitude is at most `4*log(rho)`. Assign the four residuals the
contrast's sign pattern times `log(Q)/4`; subtracting them leaves zero
interaction and hence a separable positive matrix. The bound is attained.
The analogous two-residual argument gives `delta`. These are exact algebraic
bounds for the specified positive aggregate cells, not novel biological laws.

## 4. Intervals, missing uncertainty and verdict rule

For positive input boxes `[L_X,U_X]`, monotonicity gives

`q_low = L_BM*L_CH/(U_BH*U_CM)`,
`q_high = U_BM*U_CH/(L_BH*L_CM)`.

Set `d_min = max(0, log(q_low), -log(q_high))` and report
`D_min=exp(d_min)`, `delta_min=exp(d_min/2)`,
`rho_min=exp(d_min/4)`. The bounds are sharp over the specified box.

Primary box: use only printed confidence limits attached to the primary
blood estimators; hold colon rates at their printed points because matching
intervals are absent. Call this **conditional printed-interval compatibility**.
It is not full uncertainty propagation, a simultaneous confidence region, a
p-value, or a calibrated significance test. Correlations and systematic errors
are unknown. Missing intervals do not assert that a rate is known exactly.

Separate rounding sensitivity: widen each printed blood interval endpoint
outward by half its last printed decimal unit and each colon point by half
its last printed unit. Store the printed resolutions with the input rows.
This is a conventional rounding envelope, not additional sampling confidence.
Do not substitute unrequested, more precise table extrema for the specified
colon paragraph.

Use exactly these interpretation categories:

* **Separable clock consistent** means `1` lies in the declared box for Q;
  it means conditional compatibility, not confirmation.
* **Rejected under the stated measurement assumptions** means `1` is outside
  that box, with a quantified positive departure bound.
* **Untestable** applies when required rate/scale/eligibility information is
  unavailable. The unconditional biological separability claim is also
  untestable if assay, genotype, age-window or normalization interactions are
  unrestricted; distinguish this from the conditional arithmetic verdict.

Failure rejects the **separable passenger-clock proxy under its measurement
assumptions**, not mutation-based cancer models or functional driver clocks.

## 5. Separately labelled assay/normalisation sensitivity

Write `r_observed(st)=k_st*r_latent(st)`. Then
`Q_observed=A*Q_latent`, with
`A=k_MB*k_HC/(k_HB*k_MC)`. The null needs `A=Q_observed`.
Report this break-even measurement interaction and the latent-Q box
`[q_low/K, q_high*K]` for the frozen **hypothetical** bands
`A in [1/K,K]`, `K in {1,2,4,8}`. Report every band, including any rescue.
They are stress scenarios, not measured biases or confidence bounds.

Consistent species-only per-genome/per-base conversions cancel in Q, as do
tissue-only scale factors. Inconsistent assay normalization can create an
interaction. Do not invent genome sizes or call a mixed per-base/per-genome
ratio dimensionless. The primary comparison assumes the published
genome-wide per-cell/per-genome mutation counts are commensurate within each
species comparison. Unknown callable-genome, ploidy, age-window and genotype
differences remain limitations even if an overall conversion cancels.

If both printed Kapadia global per-base rates are available, compute the
dimensionless **targeted-duplex/phylogenetic rate ratio** `s` and its marginal
endpoint ratio box. Retain the source's literal unit wording and explain the
intended common scale; this is not a new per-genome annual observation.
As an explicitly additional bridge, report `Q_assay=Q*s` and the box formed
by multiplying positive endpoint bounds. This asks how the primary contrast
would change if that within-mouse assay ratio multiplied the primary mouse
blood estimate. Dependence between estimates and the model-based targeted
panel prevent calling the result an independent replication or confidence
interval. The matched whole-blood nanorate burden comparison is qualitative
unless its own annual slope is printed; do not infer a slope from age-specific
burdens. Keep genetic backgrounds and selected mutant clones distinct from
species-level passenger-rate estimates; do not pool genotype effects.

## 6. Frozen human/mouse calibration and ferret transfer

This analysis uses the preprint only as a **transport negative control**,
separate from the peer-reviewed primary comparison. Let `N_s` be each
glutamatergic-neuron annual SNV rate and `C_s` the corresponding colon rate.
First report the human/mouse neuron/colon interaction using the same Q and
departure definitions with neurons replacing blood.

Impose the separability-implied unit log slope and fit **only** an equal-weight
human/mouse log offset:

`alpha = (log(N_H/C_H) + log(N_M/C_M))/2`,
`a = exp(alpha)`, `predicted_N_F = C_F*a`.

Report both training ratios `N_H/C_H` and `N_M/C_M`, both signed training
log residuals, and their log RMSE. The offset is a projection even if the
training ratios already disagree; it is not an established calibration.
Do not fit a free slope, choose weights from ferret results, tune the colon
estimator, or use the ferret neuronal outcome in fitting or model selection.

Score the ferret prediction using `observed/predicted`, signed natural-log
error `e=log(N_F/predicted_N_F)`, and absolute fold error `exp(abs(e))`.
Preserve the error's direction. Also show human-only and mouse-only anchor
predictions `C_F*N_H/C_H` and `C_F*N_M/C_M` as predeclared diagnostics of
incompatible calibration; do not select one by its ferret error.

Use printed SEMs only for a separate **plus/minus one-SEM sensitivity**.
For every neuronal row with a numerical SEM use `[rate-SEM,rate+SEM]`;
otherwise hold the point and explicitly label its uncertainty unpropagated.
Do not multiply SEM by a normal quantile, infer degrees of freedom, or report
a confidence level or z-test. Recompute the training-offset/prediction range
by monotone endpoints, and the observed/predicted range by dividing opposite
endpoints. If an endpoint is nonpositive, report the log bound unavailable
instead of clipping. Ferret colon uncertainty stays missing. Check whether
the prediction and observed one-SEM ranges intersect, without a significance
claim. Repeat the human/mouse neuron interaction using this same descriptive
box, keeping all colon inputs fixed.

This is a computational holdout of a public, preselected value, not blinded
or independent biological validation. The ferret mean and regression slope
must remain separate; use the latter only as an audit row, not as an alternate
fit selected after inspecting the neuronal outcome.

## 7. Falsification, implementation and delivery

Retain hostile-review attempts concerning source transcription, rounded
versus precise estimates, missing confidence limits, SEM misuse, source
reuse, genotype/taxon alignment, assay and genome normalization, age slopes
versus average burden/age, shared study/cohort evidence, calibration leakage,
and prior-art duplication. Downgrade unproved claims. An unlimited species by
tissue measurement factor can always rescue the latent restriction; quantify
the required factor instead of presenting assay effects as excluded causes.
Do not turn a rescued sensitivity into confirmation.

Only create files within `research/mutation-clock-separability/` and the one
new test `tests/test_research_mutation_clock_separability.py`. Do not modify
any existing file, workflow, research index or other research directory.
Implement standard-library Python only, offline, with the new test completing
in under 60 seconds. Recompute every reported quantitative result from CSVs
and the explicitly frozen formulas/constants. Generate the numerical README
block and offer a nonwriting `--check` command. No hardcoded empirical result
may substitute for recomputation.

Hash input files after normalizing UTF-8 CRLF and lone CR to LF; build **every
path key with `Path.as_posix()`**. Verify newline variants, scale invariance,
sharp-bound attainment, interval directions, exact/separable synthetic
controls, an assay-rescue case, and ferret-outcome perturbation leaving the
calibration and prediction unchanged. Independently recompute quantitative
outputs in the test rather than testing only cached values against the same
implementation. Include metadata/units/missingness checks. These checks do
not establish source truth, biological validity or external review.

Run repository validation and existing tests as required. Before the final
push, review every retained rate against its original readable source and
review the diff for the strict path scope. Push the branch, open a PR to
`main`, observe the actual CI on the published head, and report its number
and status. All claims and artifacts remain nonclinical and unreviewed.
