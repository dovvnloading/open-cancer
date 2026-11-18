# Does a species mutation-clock factor transport across tissues?

**Nonclinical, unreviewed computational research. No medical advice.**

## Verdict

**Rejected under the stated measurement assumptions:** the mouse/human
passenger-substitution ratio in cultured blood stem/progenitor cells is
incompatible with the specified printed colon ratio, even across the blood
interval endpoints and the separate rounding envelope. This is a deterministic
compatibility result conditional on the colon point inputs and comparable
measurement scales. It is **not** a significance test or an unconditional
biological rejection: colon-rate uncertainty and species-by-tissue assay,
genotype and age-window effects are not bounded by these aggregates.

The frozen human/mouse calibration also overpredicts the ferret neuronal
rate; the separately labelled one-SEM prediction and observation ranges do
not overlap. **A strictly matched ferret colon-to-neuron test is untestable**:
source review found that the ferret intestinal sample was from the **small
intestine**, not colon. The computed ferret score is therefore a proxy transport
diagnostic with an extra, unverified intestinal-site equivalence assumption.
All numerical results and sensitivities are generated below from the CSVs.

Failure concerns the **separable passenger-clock proxy**, under its measurement
assumptions. It does **not** reject mutation-based cancer models or functional
driver clocks. Postmitotic neurons are a transport negative control, not
cancer targets. These sparse, selected species/tissue comparisons do not
validate a general species-level clock or supply clinical guidance.

## Frozen question and audit trail

This note runs [ranked job 4](../README.md) and audits the
colon-to-whole-animal bridge discussed in [PR #75's note](../petos-paradox-mutation-rates/).
The [protocol](PROTOCOL.md) was committed by itself locally as
`3dea146499a321f6fed5e175b2610c9884ac89ca` before any empirical calculation.
The branch starts at then-current `main`,
[`09d916c`](https://github.com/dovvnloading/open-cancer/commit/09d916cdea7053f25643c6d0bd9c45827436e315).
The protocol remains byte-for-byte unchanged after newline normalization.
Ordinary Git push lacked a local GitHub credential, so the connected GitHub
API later published the identical protocol-only tree in
[`15376bd`](https://github.com/dovvnloading/open-cancer/commit/15376bd00a9e1b5f954d53a2cc39f66793a04137).
That publication commit retains the original local commit payload in its
message. Its tree is exactly `f5ddac4c38de1e9672cddf568d5b0160c2d10e96`,
the original freeze tree. Local commitment preceded calculations; **public
publication followed them**. This distinction is explicit rather than
claiming a publicly timestamped pre-analysis registration.
The sources and expected tissue difference were already visible; this is a
pre-analysis computational commitment, not blinded research or independent
registration. No outcome-based exclusion, estimator selection or refitting
was performed.

**Post-freeze source correction, retained rather than concealed.** The task,
inherited audit and protocol call the Ferret input a colon rate. During
hostile source review, after the initial computation, Cagan's **Detection of
somatic mutations across species**, paragraph beginning “A colon sample,”
showed that the ferret was the small-intestine exception. The article motivates
its inclusion using human intestinal-site similarity; that is not a measured
ferret site-conversion factor. [C] The CSV tissue field is corrected, while
the protocol and its `ferret_colon_mean` / `ferret_colon_slope` row IDs remain
traceable. The selected values and frozen scoring formula are unchanged.
The interpretation is downgraded: matched-colon ferret validation is
untestable, and the numerical result is conditional on a substituted-tissue
bridge. The existing index and PR #75 files were not altered.

## Restriction and minimum departure

For positive annual rates, the restriction is

`r(species,tissue) = species_factor * tissue_factor`.

Write `B_M, B_H` for mouse/human blood rates and `C_M, C_H` for the colon
rates. The fixed orientation is

`Q = (B_M/B_H)/(C_M/C_H)`.

Separability requires `Q=1`. Three different measures avoid an ambiguous
claim about a required fold change:

| Measure | Definition | Meaning |
|---|---|---|
| Ratio discrepancy | `D=max(Q,1/Q)` | Fold disagreement between the tissue-specific mouse/human ratios. |
| Common-ratio allowance | `delta=D**(1/2)` | Minimum symmetric allowance around a common species ratio, when each tissue ratio may vary. |
| Per-cell allowance | `rho=D**(1/4)` | Minimum symmetric allowance around a fitted separable matrix, when every species-tissue rate may vary. “Cell” here means an entry in the rate matrix. |

The last quantity is the principal minimum departure bound. If
`r_st=a_s*b_t*u_st` with every residual factor in `[1/rho,rho]`, then
`log(Q)=log(u_MB)-log(u_HB)-log(u_MC)+log(u_HC)` implies
`abs(log(Q)) <= 4*log(rho)`. Equality is attainable: assign signed residuals
`log(Q)/4` according to that contrast. Removing them leaves zero log
interaction, hence a positive separable matrix. The corresponding argument
with two ratio residuals proves the common-ratio allowance. These are sharp
algebraic bounds for this measurement model, not a new biological theorem.

For positive input boxes, the Q endpoints use opposite denominator endpoints.
The closest admissible log interaction to zero is
`d_min=max(0,log(q_low),-log(q_high))`, giving
`rho_min=exp(d_min/4)`. Exact formulas and all prespecified scenarios appear
in [PROTOCOL.md](PROTOCOL.md). Marginal printed confidence intervals do not
form a simultaneous confidence region; the test assigns no p-value.

## Source and measurement decisions

All rate inputs come from article HTML text or the exact GitHub text table.
Every CSV row records source, version, exact locator, access date, units,
estimator, assay, genotype, uncertainty/missingness and reuse basis.

| File | Included evidence and source choice |
|---|---|
| [blood_rates.csv](data/blood_rates.csv) | Kapadia's annual HSC/MPP colony result and Mitchell's precise Results estimate, each with its own printed confidence limits. The rounded Mitchell abstract is not treated as another observation. [K,M] |
| [colon_rates.csv](data/colon_rates.csv) | Cagan's specifically requested rounded human/mouse comparison. The pinned Compton-hosted Ferret mean is the predictor; its separate regression slope is an audit-only row. The ferret tissue is explicitly small intestine. [C,F] |
| [neuron_rates.csv](data/neuron_rates.csv) | The requested preprint version's glutamatergic rates and printed SEMs. Missing numerical human uncertainty remains missing. No graphical intervals are digitized. [N] |
| [assay_rates.csv](data/assay_rates.csv) | The rounded mouse per-base normalization and the distinct targeted-duplex model / phylogenetic comparison. These do not replace the primary per-cell rates. [K] |

The primary blood estimates are age slopes from cultured stem/progenitor
colonies. Cagan's selected species means summarize adjusted burden divided
by age, with its methods specifying individual-level averaging. Treating
these quantities as commensurate annual clocks is an observation assumption;
the source fits, ages, intercept treatment and cohorts are not harmonized.
The studies do not supply a matched multi-tissue assay panel. [K,M,C]

Blood's primary mouse cohort is wild-type C57BL/6J females. The neuron
preprint uses BALB/c and C57BL/6 hybrids; the selected colon paragraph does
not establish a matching strain. Kapadia's other strain, mutant-clone and
perturbation comparisons do not supply replacement genotype-specific
genome-wide annual slopes. No genotype effect is pooled into a species
clock, and its possible contribution cannot be excluded. [K,N]

Genome-wide per-cell/per-genome counts and per-base rates remain distinct in
the CSVs, including a literal-source-unit field. A consistent species-only
genome-size conversion cancels from Q; a consistent tissue-only scale also
cancels. Different callable regions, ploidy conventions, mutation callers,
genotypes or normalizations can still create a species-by-tissue interaction.
The preprint uses its species' own genome sizes and excludes mitochondrial
and sex chromosomes. No unreported genome size is supplied here. [C,K,N]

Kapadia's targeted-duplex paragraph prints an unusual `bp per year` unit.
Its earlier per-base normalization supports the common-scale interpretation,
but the literal wording is retained. The targeted-panel global estimate
depends on a branching model, background synonymous/intronic mutations and
population assumptions. The displayed ratio to the phylogenetic estimate is
a sensitivity, not an independent blood-rate replication. The interval level
and construction are not reiterated in that local paragraph, so those rows
are labelled reported CI rather than assigning a stronger confidence claim.
Matched whole-blood nanorate WGS supplies a burden check without a printed
annual slope in the specified section; it is not converted into a rate. [K]

<!-- BEGIN COMPUTED RESULTS -->

## Recomputed aggregate results

Nonclinical and unreviewed. Values below are generated from the CSVs; display precision is not measurement precision.

Original local protocol freeze commit: `3dea146499a321f6fed5e175b2610c9884ac89ca`.

### Transcribed numerical inputs

`SBS_per_cell_or_genome_year` preserves the source cell/genome distinction; the primary comparison assumes these genome-wide counts are commensurate. The `_interpreted` per-base units preserve the separate unit caveat.

| Input ID | Source tissue | Printed annual rate | Unit | Printed uncertainty | Role |
|---|---|---:|---|---|---|
| `mouse_hsc` | blood_HSC_MPP | 45.3 | `SBS_per_cell_or_genome_year` | 95% CI [42.2, 48.4] | primary |
| `human_hsc` | blood_HSC_MPP | 16.8 | `SBS_per_cell_or_genome_year` | 95% CI [16.5, 17.1] | primary |
| `human_colon` | colon_crypt | approximately 47 | `SBS_per_cell_or_genome_year` | Not printed | primary_and_transfer_training |
| `mouse_colon` | colon_crypt | approximately 796 | `SBS_per_cell_or_genome_year` | Not printed | primary_and_transfer_training |
| `ferret_colon_mean` | small_intestinal_crypt | 496.3995 | `SBS_per_cell_or_genome_year` | Not printed | transfer_predictor |
| `ferret_colon_slope` | small_intestinal_crypt | 496.3765 | `SBS_per_cell_or_genome_year` | Not printed | audit_only_excluded_from_analysis |
| `human_neuron` | cortical_glutamatergic_neuron | 18 | `SNV_per_genome_year` | Not printed | transfer_training |
| `mouse_neuron` | cortical_glutamatergic_neuron | approximately 18.6 | `SNV_per_genome_year` | SEM 5.9 | transfer_training |
| `ferret_neuron` | cortical_glutamatergic_neuron | 24 | `SNV_per_genome_year` | SEM 9.6 | transfer_outcome |
| `mouse_hsc_per_base` | blood_HSC_MPP | 8.3e-9 | `SBS_per_bp_year` | 95% CI [7.7e-9, 8.9e-9] | normalization_audit_only |
| `mouse_phylogeny_per_base` | blood_HSC_MPP | 8.29e-9 | `SBS_per_bp_year_interpreted` | CI [7.73e-9, 8.85e-9] | assay_sensitivity_denominator |
| `mouse_duplex_per_base` | blood_targeted_panel | 11.77e-9 | `SBS_per_bp_year_interpreted` | CI [9.28e-9, 14.94e-9] | assay_sensitivity_numerator |

| Input ID | Source unit wording | Retained uncertainty status |
|---|---|---|
| `mouse_hsc` | SBSs per year in HSC/MPP colonies | printed 95% confidence interval |
| `human_hsc` | substitutions per cell per year | printed 95% confidence interval |
| `human_colon` | substitutions per genome per year | no numerical interval or SEM printed for this rate in selected passage |
| `mouse_colon` | substitutions per genome per year | no numerical interval or SEM printed for this rate in selected passage |
| `ferret_colon_mean` | SBS/genome/year | mean interval absent; source slope CI is NA and not assigned to mean |
| `ferret_colon_slope` | SBS/genome/year; annual zero-intercept slope | Zero-intercept regression slope (95% CI) is literal NA |
| `human_neuron` | SNVs/year; genomewide burden age slope | no numerical interval or SEM printed in selected Results |
| `mouse_neuron` | SNVs/year; genomewide burden age slope | printed SEM; not a confidence interval |
| `ferret_neuron` | SNVs/year; genomewide burden age slope | printed SEM; not a confidence interval |
| `mouse_hsc_per_base` | per base pair per year | printed 95% confidence interval |
| `mouse_phylogeny_per_base` | bp per year (literal source wording) | reported CI; level and construction not reiterated locally |
| `mouse_duplex_per_base` | bp per year (literal source wording) | reported CI; level and construction not reiterated locally |

### Peer-reviewed blood/colon comparison

| Statistic | Value |
|---|---:|
| Mouse/human blood ratio | 2.69643 |
| Mouse/human colon ratio | 16.9362 |
| Q: blood ratio / colon ratio | 0.159211 |
| Signed log interaction | -1.83752 |
| D: fold disagreement | 6.28096 |
| delta: minimum symmetric allowance per tissue-specific species ratio | 2.50619 |
| rho: minimum symmetric allowance per species-tissue cell | 1.58309 |

| Box | Q range | D minimum | delta minimum | rho minimum | Conditional verdict |
|---|---|---:|---:|---:|---|
| Printed blood intervals; colon points fixed | [0.145714, 0.173199] | 5.77369 | 2.40285 | 1.55011 | Rejected under the stated measurement assumptions |
| Outward rounding sensitivity | [0.143483, 0.175866] | 5.68616 | 2.38457 | 1.5442 | Rejected under the stated measurement assumptions |

The printed-interval result is conditional compatibility, not a simultaneous confidence test or a p-value. Missing colon uncertainty is unpropagated. The unconditional biological verdict with unrestricted measurement interactions is **Untestable**.

### Assay/normalisation sensitivity

The latent separable null requires measurement interaction A = 0.159211 at the point estimates; its conditional compatible A range is [0.145714, 0.173199].

| Hypothetical interaction band | Latent Q range | D minimum | delta minimum | rho minimum | Conditional verdict |
|---|---|---:|---:|---:|---|
| A in [1/1, 1] | [0.145714, 0.173199] | 5.77369 | 2.40285 | 1.55011 | Rejected under the stated measurement assumptions |
| A in [1/2, 2] | [0.072857, 0.346399] | 2.88685 | 1.69907 | 1.30348 | Rejected under the stated measurement assumptions |
| A in [1/4, 4] | [0.0364285, 0.692797] | 1.44342 | 1.20143 | 1.0961 | Rejected under the stated measurement assumptions |
| A in [1/8, 8] | [0.0182142, 1.38559] | 1 | 1 | 1 | Separable clock consistent |

These bands are hypothetical stress scenarios. A rescued band establishes compatibility within the allowance, not biological confirmation.

The targeted-duplex/phylogenetic per-base rate ratio s is **1.41978**, with marginal endpoint ratio box **[1.04859, 1.93273]**. Applying it to the primary mouse blood rate is an additional assay bridge, not a matched whole-blood annual observation or independent replication.

| Bridged point statistic | Value |
|---|---:|
| Q | 0.226045 |
| log_Q | -1.48702 |
| D | 4.42389 |
| delta | 2.1033 |
| rho | 1.45028 |

| Bridged box | Q range | D minimum | delta minimum | rho minimum | Conditional verdict |
|---|---|---:|---:|---:|---|
| Q multiplied by s endpoint box | [0.152794, 0.334747] | 2.98733 | 1.72839 | 1.31468 | Rejected under the stated measurement assumptions |

### Preprint neuron transport negative control

Postmitotic neurons are a transport negative control, not cancer targets. These results remain separate from the peer-reviewed blood/colon comparison.

| Human/mouse neuronal comparison | Value |
|---|---:|
| Mouse/human neuron ratio | 1.03333 |
| Mouse/human colon ratio | 16.9362 |
| Q | 0.0610134 |
| log_Q | -2.79666 |
| D | 16.3898 |
| delta | 4.04844 |
| rho | 2.01207 |

| Descriptive box | Q range | D minimum | delta minimum | rho minimum | Box verdict |
|---|---|---:|---:|---:|---|
| Plus/minus one SEM; missing uncertainties fixed | [0.0416597, 0.0803671] | 12.4429 | 3.52745 | 1.87815 | Rejected under the stated measurement assumptions |

| Frozen human/mouse calibration | Value |
|---|---:|
| Human neuron/colon ratio | 0.382979 |
| Mouse neuron/colon ratio | 0.0233668 |
| Equal-weight log offset alpha | -2.35811 |
| Multiplicative offset exp(alpha) | 0.0945992 |
| Human signed training log residual | 1.39833 |
| Mouse signed training log residual | -1.39833 |
| Training log RMSE | 1.39833 |

### Ferret proxy transfer (small-intestine predictor)

Matched-colon ferret transfer verdict: **Untestable**. Postfreeze source review identifies the ferret predictor as small intestine, while the human/mouse training predictors are colon. The frozen numerical rule is retained as a proxy score conditional on an unverified ferret small-intestine=colon bridge. The original `ferret_colon_*` row IDs preserve protocol traceability; they do not assert the source tissue is colon.

| Ferret proxy holdout and frozen anchor diagnostics | Value |
|---|---:|
| Predicted annual neuron SNV/genome rate | 46.959 |
| Observed annual neuron SNV/genome rate | 24 |
| Observed/predicted | 0.511084 |
| Signed natural-log error | -0.67122 |
| Absolute fold error | 1.95662 |
| Human-only anchor prediction, SNV/genome/year | 190.11 |
| Mouse-only anchor prediction, SNV/genome/year | 11.5993 |

| Plus/minus one-SEM transfer sensitivity | Range or result |
|---|---|
| Log offset alpha | [-2.54889, -2.22035] |
| Multiplicative offset | [0.0781687, 0.108571] |
| Predicted ferret SNV/genome/year | [38.8029, 53.8946] |
| Observed ferret SNV/genome/year | [14.4, 33.6] |
| Observed/predicted ratio | [0.267188, 0.865915] |
| Prediction and observation ranges intersect | No |

SEMs were not converted to confidence limits. Neuronal uncertainty left unpropagated: `human_neuron`. Intestinal predictor uncertainty is unpropagated throughout the neuron analysis. The ferret small-intestine mean is used as a colon proxy; the separate regression slope is audit-only. The ferret outcome does not enter calibration or model selection.

### Input integrity

SHA-256 hashes use newline-normalized UTF-8 content and POSIX relative paths. They establish input integrity, not source truth or biological validity.

| Input | SHA-256 |
|---|---|
| `PROTOCOL.md` | `d3dc484ae6ee879c09dcdfc209311fc6ef66ebf2fa49350f1b1cdb8148b18191` |
| `data/blood_rates.csv` | `ead2fbd4b19935e0b7578076165d510f2bb4717685dceb177c71852f77078e21` |
| `data/colon_rates.csv` | `25f05c5684a533dd04484887ed802a57cac6e9c515805fc6296af414c08452a6` |
| `data/neuron_rates.csv` | `d0767ea61b9e65034afe2a32052075b158cac9f11ea281eeaa0874d114e341b1` |
| `data/assay_rates.csv` | `248eaea63165c8460b2613b27a1fe4fa37b046f2d13b6b8984b47e6d7969d6a7` |

<!-- END COMPUTED RESULTS -->

## How to interpret the sensitivities

The primary endpoint box fixes colon estimates only to make the conditional
calculation possible. It does not turn missing colon uncertainty into zero
uncertainty. The rounding envelope accounts for displayed precision under
the conventional nearest-unit assumption; it is not sampling uncertainty.

Let `r_observed(st)=k_st*r_latent(st)` and
`A=k_MB*k_HC/(k_HB*k_MC)`. Then `Q_observed=A*Q_latent`.
The reported break-even value is the interaction a measurement explanation
would require under a separable latent clock, not an estimated laboratory
bias. The fixed hypothetical bands include both rejection and rescue.
Unrestricted `A` can always fit the observed Q, so the unrestricted latent
biological claim remains **untestable**. The targeted-duplex substitution
also adds an unverified transfer assumption and uses dependent estimates;
its endpoint products are sensitivity bounds, not a confidence interval.

For the neuronal control, the unit-slope calibration is

`alpha = (log(N_H/C_H)+log(N_M/C_M))/2`,
`predicted_N_F = C_F*exp(alpha)`.

The human/mouse training ratios already disagree, so the offset is a fixed
least-squares projection, not evidence that training established separability.
The human-only and mouse-only predictions expose this instability; neither
is selected using the ferret outcome. The frozen pooled prediction uses only
human/mouse neuronal outcomes and the ferret intestinal predictor. Changing
the ferret neuronal value cannot change the fitted offset or prediction.

The one-SEM calculations apply the printed SEMs directly. They are descriptive
boxes with unpropagated human-neuron and intestinal-rate uncertainty, not
confidence intervals or a formal rejection test. The public ferret value was
preselected by the task, making this a computational holdout rather than
blinded validation. Its source labels also differ: `Mustela_putorius` in the
table and `Mustela putorius furo` in the preprint. Both are retained; the
species-level mapping and intestinal-site substitution are transport
assumptions, not demonstrated cohort equivalence. [F,N,C]

## Hostile self-review and falsification attempts

| Attempt to defeat or overstate the result | Outcome retained |
|---|---|
| Verify every printed rate against its source rather than trusting inherited transcriptions. | All retained values and uncertainty strings were checked against their specified HTML passages or the pinned GitHub row. This revealed the ferret tissue error; the metadata and interpretation were corrected without changing the frozen protocol or numeric rule. |
| Make the ferret row an exact matched-colon validation. | Fails source eligibility: it is small intestine. The matched test is untestable; the reported proxy score needs an additional site-equivalence assumption. |
| Repair primary separability with the available blood intervals or source rounding. | Neither frozen envelope includes the separable value. The quantitative bound remains conditional on unreported colon uncertainty. |
| Treat missing intervals as evidence of exact measurement. | Rejected. Fixed values are a computational condition; the unconditional uncertainty claim is untestable. |
| Replace cultured blood with a duplex estimate silently. | Rejected. Per-base model estimates are separate; the dimensionless within-mouse assay ratio enters only the predeclared sensitivity. No human duplex slope was supplied by the specified passages. |
| Attribute the interaction uniquely to evolutionary tissue effects. | Rejected. Genotype, cohort, age-window, normalization and assay interactions remain possible. The broader hypothetical measurement band rescues compatibility, and that rescue is reported. |
| Convert SEMs to confidence limits or infer a z-test. | Rejected. The descriptive one-SEM exercise propagates only the printed SEMs and retains missing human and intestinal uncertainty. |
| Improve ferret fit by choosing its intestinal estimator or fitting its neuronal value. | Rejected. The mean, unit slope, equal training weights and score were frozen; the alternate slope is audit-only. Outcome-perturbation tests leave calibration and forecast unchanged. |
| Conceal incompatible training anchors behind a pooled forecast. | Both training residuals and both anchor predictions are reported. A projection can be computed even when its exact model already fails on training data. |
| Claim a new discovery that tissues differ or that colon clocks need not describe the whole animal. | Rejected. Kapadia explicitly states the blood/colon difference; Caglayan explicitly situates neurons after that finding. These are prior observations. [K,N] |
| Use a passed test, hash or AI review as scientific validation. | Rejected. They check source-to-computation consistency and implementation properties; they establish neither source truth, priority nor biological validity. |

## Novelty and implications for the bridge

**Novelty assessment: a reproducible quantitative audit, not a new biological
discovery.** The component tissue contrasts are already published, and the
ratio interaction / multiplicative minimax bound is elementary algebra. The
added work is the frozen application, explicit scale/uncertainty bookkeeping,
ferret outcome holdout, source-eligibility correction and executable negative
evidence. No first-in-literature claim or exhaustive literature search is
asserted. [K,N; prior repository discussion in PR #75]

This audit limits a particular shortcut from intestinal passenger counts to
a whole-animal functional clock. **For PR #75, transporting colon passenger
mutation rates to whole-animal cancer hazards therefore requires an independently
justified species-by-tissue calibration; the colon rate alone cannot serve as a
universal species-wide mutation clock or an observed driver rate.** It does not
resolve Peto's paradox, establish the relevant driver rate, count oncogenic
stages, or quantify cancer risk.
Different tissues can have different species interactions while mutation-based
cancer mechanisms remain operative. The neuron result is useful precisely
as a stress test of transport, with its extra intestinal-site limitation
displayed, and supplies no cancer-target or treatment inference.

## Access, attribution and reuse

Access-date label: **2026-11-18, America/New_York**. Only small printed
aggregates are retained. No patient-level or individual-animal records, raw
genomes, plotted-coordinate extraction, PDFs, supplementary binaries or browser
file downloads were used. No download prompt appeared or was accepted.
Direct access to some PMC/publisher pages returned access challenges; indexed
article HTML supplied the missing readable passages. A publisher browser
landing page also displayed a request-limit notice; no access control was
bypassed. These limitations did not authorize a different numerical source.

Kapadia, Mitchell and Cagan state CC BY terms in their article HTML. The
Compton table route is documented in
[PR #75's provenance](../petos-paradox-mutation-rates/data/PROVENANCE.md):
Cagan's attributed rate columns reproduced in the separately licensed
Compton preprint-linked table. Only the requested Ferret rate fields were
retained; no Compton cancer-outcome data are redistributed.

The Caglayan article-specific license could not be read on the blocked
landing page. The publisher's
[single-DOI metadata](https://api.biorxiv.org/details/biorxiv/10.64898/2026.09.22.753624/na/json)
confirmed the requested title, version/date and `license: cc_by`.
That API was used **only for rights/version verification**; every neuronal
rate remains transcribed from article HTML. The metadata does not specify
the license version, so this note does not invent one. Attribution is retained
and the extraction into labelled aggregate rows is the stated transformation.
No source article text, figures or full source table are redistributed.
The repository MIT license does not relicense these source aggregates.

### Primary sources and exact locators

**[K]** Kapadia CD et al. (2025), *Clonal dynamics and somatic evolution of
haematopoiesis in mouse*, Nature.
[Article HTML](https://pmc.ncbi.nlm.nih.gov/articles/PMC12074984/),
[DOI](https://doi.org/10.1038/s41586-025-08625-8).
Results **Somatic mutation accumulation**, first and final paragraphs;
**Fitness landscape of murine clonal haematopoiesis**, paragraph beginning
“In the targeted gene panel”; Discussion paragraph beginning “Somatic
mutation rates scale inversely”; Methods **Cohort**, **Mutation burden
estimation**. The
[correction](https://doi.org/10.1038/s41586-025-09635-2) concerns an axis
label and does not change these printed rate inputs.

**[M]** Mitchell E et al. (2022), *Clonal dynamics of haematopoiesis across the
human lifespan*, Nature.
[Article HTML](https://pmc.ncbi.nlm.nih.gov/articles/PMC9177428/),
[DOI](https://doi.org/10.1038/s41586-022-04786-y).
Results **Mutation burden and telomere lengths**, opening sentence, Fig. 1b.

**[C]** Cagan A et al. (2022), *Somatic mutation rates scale with lifespan
across mammals*, Nature.
[Article HTML](https://pmc.ncbi.nlm.nih.gov/articles/PMC9021023/),
[DOI](https://doi.org/10.1038/s41586-022-04618-z).
Results **Correlation with life-history traits**, opening human/mouse rate
comparison; **Detection of somatic mutations across species**, ferret
small-intestine exception; Methods **Estimation of mutation rate** and
**Association of mutation rate with life-history traits**.

**[F]** Cagan aggregate fields reproduced in the Compton-hosted
[Mutation_Data.csv](https://github.com/zacharycompton/cancerAcrossVertebrates/blob/3b4db10c7c38c63aeebfc9a48ce4bfa376dfd2bc/Mutation_Data.csv),
commit `3b4db10c7c38c63aeebfc9a48ce4bfa376dfd2bc`, row
`common_name=Ferret`, mean annual mutation rate and separate zero-intercept
slope / missing confidence interval columns. Reuse route:
[Compton et al. preprint version 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10350200/),
[DOI](https://doi.org/10.21203/rs.3.rs-3117313/v1), and
[Cagan's institutional license record](https://doi.org/10.17863/CAM.83740).

**[N]** Caglayan E et al., *Somatic mutations and single-nucleus transcriptomics
reveal uniquely human properties of neuronal aging*, preprint **version 1,
September 24, 2026**, not peer reviewed.
[Article HTML](https://pmc.ncbi.nlm.nih.gov/articles/PMC13622217/),
[DOI](https://doi.org/10.64898/2026.09.22.753624).
Results paragraph beginning “We performed single-nucleus sorting,” Fig. 1f;
**Similar neuronal mutation rates in mammalian species with vastly different
lifespans**, mouse/ferret sentences, Figs. 2i-l; Methods **Genomewide mutation
rate, mutation burden estimates**. Numeric data were readable in indexed
article HTML despite the direct-page challenge.

## Reproduction

From the repository root, with Python's standard library and no network:

```bash
python -B research/mutation-clock-separability/analysis.py
python -B research/mutation-clock-separability/analysis.py --markdown
python -B research/mutation-clock-separability/analysis.py --check
python -B -m unittest discover -s tests -p 'test_research_mutation_clock_separability.py' -v
```

The default command emits the complete numerical result and metadata as JSON;
`--markdown` prints the generated block, and `--check` compares that block to
this README without writing files. The
[new offline test](../../tests/test_research_mutation_clock_separability.py)
independently recomputes the results from the CSVs, checks sharp-bound
attainment and monotonic endpoint choices, exercises compatible and rejected
synthetic controls, and checks assay rescue and ferret outcome leakage.
It also verifies normalized-content hashes under LF, CRLF and lone-CR content
and POSIX path keys on Windows and Linux. The existing repository workflow
discovers the new test; no workflow or existing file is changed.

Hashes identify the local inputs, not their truth. CI and self-review are
implementation checks, not external scientific review. All use remains
nonclinical, unreviewed and subject to [DATA_POLICY.md](../../DATA_POLICY.md).
