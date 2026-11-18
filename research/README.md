# Research programme: independently audited mathematical notes

**Scope and status (17 November 2026).** This index covers only the three
research notes audited here. Each proof concerns its **explicit mathematical
model**, not measured tumour behavior, an approved therapeutic regimen, or a
clinical recommendation. "Proved" denotes a deductive statement inside those
assumptions; it does not certify originality, external peer review, or
biological validity. The original PRs were #52, #53, and #54. The
independent-referee derivations, dependency ledgers, boundary audits, and
limitations are recorded in the respective notes.

## Claim/status index

| Note | Precise central claim and quantifiers | Status | Novelty assessment against primary sources | Standard-library regression/test |
| --- | --- | --- | --- | --- |
| [Clonal-tree identifiability](clonal-tree-identifiability/) (PR #52) | In the exact, labelled, perfect-phylogeny **single malignant founder** model, for each *known target tree* there exists a strictly positive rational designed mixture identifying it: **one** row for a rooted path, **two** for a branching tree. The construction has an entrywise closed-radius certificate. **Arbitrary** finite samples never guarantee uniqueness when `n >= 3`, even generically under full-rank, strictly positive usages. | **Proved (within model)** | **Unverified narrow constructive synthesis; low-to-moderate mathematical originality.** Single-row obstruction, sum condition, order dimension, and ambiguity/counting bounds are prior results, not new. Neither a first identifiability theorem nor an achievable biopsy count. | [tests/test_research_clonal_tree_identifiability.py](../tests/test_research_clonal_tree_identifiability.py) |
| [Adaptive therapy / mutation supply](adaptive-therapy-mutation-supply/) (PR #53) | For the **two-individual-cell** CTMC, ideal containment C versus finite sensitive-clearance M, with initially no R and progression at total burden `>2`, both **uncensored** progression-time means are infinite. For an **independent Exp(1) horizon** and rates `(a,b,d,k,z)=(100,1,3,2,1)`, the sign of capped-time C minus M is that of `18 mu^3 - 185 mu^2 - 517 mu + 400`: C wins for positive `mu` below ~0.63735 or above ~12.44262; M wins strictly between. The two crossing roots are exact algebraic roots; the displayed decimals are approximations. | **Proved (within model)** | **Low, model-specific algebraic counterexample; priority unverified.** Mutation supply versus competitive suppression and objective dependence are established prior themes. The model does **not** identify M with a common bounded-dose MTD optimum. | [tests/test_research_adaptive_therapy_mutation_supply.py](../tests/test_research_adaptive_therapy_mutation_supply.py) |
| [Combination / cross-resistance](combination-cross-resistance/) (PR #54) | For the specified **four-type branching process**, one or finitely many initially sensitive founders, **equal per-drug exposure**, fixed finite A-then-B period count `m`, and `0<w<=1/2, 0<=p<1, c>=1`, A-first raises **extinction by a sufficiently short finite course** over half-dose combination **iff** `c>(1-w)/w`; at equality combination wins at fourth order. At `p=1` sensitive-only starts tie. A nearby finite-course critical curve exists for bounded `c`, and for `w=1/2,m=1` its third-order response to `p` is `-mu*k^2/1440` (computer-assisted exact polynomial certificate). | **Proved (model-specific; sixth-order identity computer-assisted)** | **Narrow, priority-unverified mathematical refinement only.** Four-type branching, cross-resistance, sequencing advantages, and cycling are known. Neither an arbitrary-course ranking nor a first sequencing advantage. | [tests/test_research_combination_cross_resistance.py](../tests/test_research_combination_cross_resistance.py) |

**Verification boundaries.** The first test enumerates all labelled trees through
five vertices and checks selected perturbations independently; that is not an
enumeration of all real datasets. The second compares rational CTMC
resolvents, independent jump chains, root brackets, and fixed-horizon
certificates; it does not demonstrate population-scale mutation rates. The
third compares forward population generators to PGF expansions, verifies exact
rational probability enclosures, and runs a separate generic symbolic
polynomial certificate for its sixth-order identity. No simulation proves
clinical relevance or mathematical novelty. All three tests are Python
standard library only, with bounded deterministic workloads. No patient
records or experimental datasets are included.

## Most promising open follow-ups, ranked by expected scientific impact

The ranking emphasizes changes likely to matter for actual cancer biology or
therapy **if tested**, rather than mathematical elegance or ease of proof.
Every item is an unproved research question.

1. **High — Test containment versus aggressive treatment under one physically
   admissible control model.** Include a shared, bounded dose–response law,
   finite sensitive regrowth, pharmacokinetics, mutation *per division*, and
   measurable resistant-cell competition. Ask whether any ranking reversal
   persists with tumour-relevant population sizes and *pre-specified*
   progression/survival endpoints. Parameterize and falsify with serial
   patient-derived organoids or comparable nonclinical systems before
   considering any therapeutic translation.

2. **High — Measure the joint architecture of collateral sensitivity and
   directly shared resistance.** With genotype-resolved lineage tracing and
   experimentally quantified drug-killing surfaces, estimate the sizes,
   fitness, and supply rates of A-only, B-only, and doubly resistant clones.
   Test whether the proposed short-course order effects survive realistic
   pre-existing resistance, reversible phenotypes, asymmetric sensitivities,
   drug interactions, and toxicity-matched exposures. The current scalar
   `c,p,w` model cannot answer this.

3. **High — Connect clonal identifiability to physically feasible noisy
   sampling.** With copy-number-aware variant frequencies, purity uncertainty,
   missing mutation clusters, read sampling, spatial constraints, and
   realistic biopsy locations, estimate the probability of a *correct*
   tree or defensible equivalence class. Determine whether two accessible
   specimens can ever approach the idealized designed-data bound without
   presupposing the answer.

4. **Medium–high — Establish which stochastic endpoint actually predicts
   treatment benefit.** Compare eventual progression, fixed-time progression
   risk, discounted control time, toxicity-adjusted burden, and competing
   mortality in a shared calibrated model. Determine which policy comparisons
   are robust to observation-horizon laws, which can reverse, and what
   measurement schedule would distinguish those endpoints prospectively.

5. **Medium–high — Quantify spatial and ecological feedback on resistance
   supply.** Extend the finite CTMC and multitype branching models with
   spatial refugia, variable resource competition, persisters, and immune
   pressure. Seek rigorous error controls or validated computational limits
   connecting the small-state theorems to tumour-scale heterogeneous systems.

6. **Medium — Resolve global finite-course switching under defensible drug
   constraints.** For measured PK/PD and stochastic extinction endpoints,
   characterize whether a sequencing advantage persists over *realistic*
   course lengths and initial clonal mixtures; identify any nonmonotone
   schedule thresholds and provide reproducible uncertainty certificates.
   Do not extrapolate the fixed-`m`, `T -> 0` expansion to long treatment
   courses.

## Priority anchors and reproducibility

The literature baseline includes [El-Kebir et al. (2015)](https://doi.org/10.1093/bioinformatics/btv261),
[Pradhan and El-Kebir (2018)](https://www.el-kebir.net/public/preprints/RECOMBCG18_Pradhan_El-Kebir.pdf),
[Marangola et al. (published 2026)](https://doi.org/10.1109/TCBBIO.2026.3726506),
[Hansen, Woods and Read (2017)](https://doi.org/10.1371/journal.pbio.2001110),
[Alvarez and Viossat (2024)](https://doi.org/10.1007/s00285-024-02062-3),
[Katouli and Komarova (2011)](https://doi.org/10.1007/s11538-010-9539-y),
and [Dabi et al. (2026)](https://doi.org/10.1371/journal.pcbi.1014803).
The note-specific bibliographies distinguish the known results from the
unverified-priority local refinements.

Run each linked test individually with `python -m unittest discover -s tests -p
test_research_<name>.py -v` using its actual linked filename (or use
`python -m unittest discover -s tests -p 'test_research_*.py' -v`).
The tests make no network calls and require no nonstandard packages.
