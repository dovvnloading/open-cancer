# Research programme: claims, evidence and next decisions

**Nonclinical, unreviewed research. No treatment, diagnosis, surveillance or
patient-selection advice.** This index covers all **16 research notes merged
through PR #68**, audited against `main` at
[`44c7689`](https://github.com/dovvnloading/open-cancer/commit/44c7689e6e11458a341dcc0708f5df2ef521c414).
It records mathematical results, their observation assumptions, and what the
available data actually tested. No entry establishes a new biological discovery,
clinical validity, external peer review or bibliographic priority.

**Where a discovery is most plausible next:** use observations that can separate
competing mechanisms, beginning with the existing cell-line fluctuation data in
job 1 below. The programme has strong conditional mathematics, but its main
bottleneck is identifying model inputs from observable data. The focused
pan-cancer test did **not** demonstrate improvement over TMB alone; #68's
conditional progression mechanism is likewise untestable with published
aggregates, while its substituted ten-year fitness proxy fails a specified
interval-box rule. Neither finding invalidates the narrower model theorems.

## Negative results and what they teach

These outcomes are part of the research result. **Untestable is an abstention,
not confirmation or an experimental falsification.** Mathematical
counterexamples refute the indicated extension, not the narrower theorem.
The linked notes retain the individual falsification fixtures and review logs.

### Empirical failures, contrary evidence and unavailable tests

| Result retained | What it teaches |
| --- | --- |
| **[Clonal haematopoiesis progression (#68)](clonal-hematopoiesis-progression/): the mechanistic target is untestable; the substituted baseline-fitness proxy fails its frozen ten-year rule.** In a seven-gene leave-one-gene-out proxy using Fabre/Abelson-refit gene Cox HRs, U2AF1 alone has nonoverlapping ten-year envelopes: observed HR **304.052 [134.063, 689.584]**, predicted **8.174 [1.123, 63.617]**. Descriptive log-RMSE is **2.054229** versus **2.276578** for an intercept-only comparator. The five-year proxy also fails for SRSF2; the 20-year per-gene overlap rule is consistent, but the **post-hoc** common-intercept feasibility gap remains **0.344714** (incompatible boxes). | Neither published gene-class baseline coefficients nor unmatched marginal Cox HRs identify the theorem's full clone state or conditional fixed-horizon hazard. A failed CI/HPDI-derived interval-box rule is **not** a calibrated significance test or rejection of established fitness–risk associations. The 20-year overlap is counterevidence to horizon-robust failure, **not** joint calibration. Seven eligible genes are not a patient-level prediction cohort. |
| **[Pan-cancer immune test (#67)](immune-escape-pancancer-test/): no demonstrated predictive gain.** Reciprocal-TMB held-out RMSE is **13.173 percentage points (pp)** versus **12.084 pp** for log-linear TMB (21 rows, 20 parent cancers). The fixed clonal extension worsens matched prediction; its free exponent is **−1.6857**, negative under every parent deletion. Relaxed HLA uncertainty permits either sign. The coherent full panel has nine parents, below its frozen ten-parent gate; deleting the indexed-only clonal source removes clonal/full tests. | The declared improvement criterion failed. The primary reciprocal-versus-log-TMB paired error interval crosses zero, so this is not a calibrated statistical rejection of [#59](immune-escape-neoantigen-clonality/), nor evidence that individual clonality is harmful. TMB, clonal-ID fraction and ordinary HLA loss do not measure functional-load exposure. Historical response-source overlap prevents an independent-replication claim. |
| **[Prediction confrontation (#63)](prediction-confrontation/): 13 of 14 verdicts untestable.** Only the weaker clonal-burden response association is consistent. Opposite nominal HLA-LOH outcome contrasts do not reject after the 14-contrast Holm correction; Ravi's lack of a detected HLA-LOH response association remains visible. ctDNA summaries lack the required schedule comparison or parity-thinned outcome. Iyer's lineage overdispersion lacks founder-conditioned covariance. | No unconditional corrected core prediction is validated or falsified. Ordinary HLA loss is not total presentation loss; variable founders can mimic dependent switching. A lead-time-range conflict with a perfect-periodic/common-window model depends on an unverified imaging-lag bridge. Means, medians and detected-only versus zero-imputed outcomes are different estimands. |
| **[ecDNA calibration (#65)](ecdna-segregation-clustering-data/) and [#63](prediction-confrontation/): new empirical dispersion/threshold calibration is not established.** Lundberg's LA-N-5 rejection count challenges fair allocation only under calibrated tests and independent mothers; dependence-robust bounds do not reject. Lange's eight nonsignificant KS comparisons oppose a blanket rejection but do not establish binomiality. Shibata's point fits lack recovered likelihood/covariance or calibrated dispersion intervals; conservative threshold envelopes include binomial. The source comparison across ecDNA clones and HSR controls favours common dispersion with clone-specific allocation means; HSR controls also have nonzero fitted dispersion. | Separate point estimates do not establish dispersion heterogeneity. Numerical/rounding intervals are not sampling confidence intervals. Folded sister counts can alias bias and dispersion; cross-species correlation is not within-species dispersion. Perturbation arms lack linear-killing calibration; nuclear loss/variable replication break conserved-copy transport. Corrected primary-source HSR exclusions are clones **4, 6 and 7**. |
| **[Containment evidence (#66)](containment-evidence/): all six mechanistic checks untestable**, zero eligible reports/formal mechanistic tests. P6b's broad outcome/any-published-interval reading is contradicted by COIN OS HR **1.084 [1.008, 1.165], an 80% CI**; tumour-control and 95%-tumour-control sensitivities remain untestable. González-Cao and INTERIM have adverse point directions (INTERIM's interval crosses the null); Pressley's lung xenograft comparison is adverse without pairwise uncertainty. Seyedi's single-gemcitabine result is neutral; toxicity complicates combination effects. Bacevic has no monolayer advantage despite spheroid benefit. Restricted-extraction S1320/BFR14 negatives remain linked. | Outcomes do not identify [#61](containment-at-tumour-scale/)'s mutation-supply threshold. P6b's qualified contradiction is neither a 95% rejection nor falsification of #61. Unreached ovarian medians, unequal dosing duration, fitted counterfactuals and overlapping reports cannot become independent mechanistic confirmations. |
| **[Metastatic timing (#62)](metastatic-seeding-timing/): optional broad-witness model fails conditionally.** TRACERx reports 75/90 late-divergence cases called early after single-region downsampling; under independent trials and per-case miss probability at most 1/2, the upper tail is about **4.59 × 10⁻¹¹**. Secondary-seeding aggregates also contradict universal primary-only seeding. Matched two-/three-region predictions remain externally untested. | This is a post-hoc conditional model check, not a migration-date measurement. Divergence is not migration. A one-region summary cannot calibrate two-moment designs; different pool sizes and sampling without replacement require different laws. |
| **[Incidence identifiability (#64)](multistage-incidence-identifiability/): female-breast rounded means fail the homogeneous monotone-mean condition.** The minimum abstract monotonicity calibration factor is **1.066238**; none of four sites identifies a biological transition count under unrestricted age-bin calibration. | This is a descriptive shape failure, not a significance test or rejection of multistage carcinogenesis. Counts/exposures are absent, period rates mix observation effects, and grid indices are not genomic driver counts. |

### Mathematical failures and limits on stronger claims

| Notes | Negative result and implication |
| --- | --- |
| [#52: clonal trees](clonal-tree-identifiability/) | No finite arbitrary sample count guarantees uniqueness for `n>=3`, even with positive usages, distinct frequencies and maximal rank; ambiguity occupies open sets. Dominance order alone does not fix parents. Equality breaks the one-row uniqueness inequality; tied positive columns can exceed the distinct-column bound and a proposed `n!` bound. Root convention changes counts. **Designed mixtures are not a biopsy guarantee.** |
| [#53: mutation supply](adaptive-therapy-mutation-supply/) | Both uncensored means are infinite. A universal one-cutoff rule fails for the capped objective, but the double reversal is not universal: competition, clearance and censoring matter. Eventual progression favours M in the example; fixed and exponential horizons of equal mean can reverse ranking. Correcting exploratory sensitive-cell resurrection removed the original-rate reversal; the retained witness uses different stated rates. A vanishing high-rate cubic coefficient is not a tie; lower coefficients decide. **Specify process, controls and endpoint together.** |
| [#54: combination](combination-cross-resistance/) | The short-course threshold fails globally: the same parameters reverse at longer courses. A zero cubic is not a tie, order-reversed gaps are not exact opposites, and pre-existing resistance can reverse the sensitive-founder result. Cross-resistance changes more than leading gap magnitude. Infinite-switch uniformity and global monotone thresholds are unproved. **Keep period, horizon and founder restrictions.** |
| [#55](ecdna-segregation-therapy/), [#60](ecdna-killing-laws-and-clustering/): ecDNA | Zero-copy dominance can coexist with positive-lineage persistence; normalized distributions do not identify `u`. Clustering has opposite finite-horizon extinction effects for different founders. Entire generation-two mean/Fano dose curves can coincide while thresholds differ. Saturation need not make eradication impossible or have threshold exactly `b`; the specified declining-hazard families persist at every finite intensity. **Low moments and finite survival do not identify asymptotic extinction.** #60 closes #55's linear-hazard gap; interior clustering-threshold monotonicity remains conjectural. |
| [#57: resistance origins](induced-vs-selected-resistance/) | Baseline counts and, over the full unknown-kernel class, two holiday arms have exact aliases. Singleton components, unequal kernels, induced amplification, whole-clone switching, absent reversible components/permanent memory, invisible marks and uncalibrated detection defeat stronger identification, sometimes at every holiday time. Common unknown detection aliases intensities/kernel; arm-specific detection can also alias switching. No uniform finite culture count works near degeneracies; illustrative minima depend on the test class. **Post-exposure acquisition is not drug causation or a molecular mechanism.** |
| [#58: ctDNA](ctdna-surveillance-schedule/) | Finite assay noise or multimodal growth defeats universal geometric optimality; repeat detection can favour overlap. More variance does not universally front-load dates or reduce capture: same-mean gamma and log-normal families differ. Aggregate lead time cannot identify growth, shedding or residual burden. **Specify family, assay and scoring horizon.** |
| [#59: immune escape](immune-escape-neoantigen-clonality/) | Clonal IDs do not identify per-cell loads; equal load multisets can yield equal outcomes with different clonality. Positive successful-null escape prevents almost-sure total eradication at finite load/strength. Exact finite-escape risk is not globally convex; bounded killing defeats the reciprocal law, and persistent low-load founders defeat its uniform harmonic limit. **Source suppression is not total extinction.** |
| [#61](containment-at-tumour-scale/), [#66](containment-evidence/): containment | The large-`N` limit has one crossing, not #53's two. Near-boundary finite-size simulations have opposite estimated signs; leading-coefficient equality is not a finite-`N` tie. Under their respective assumptions, no competition, nonpositive released resistant growth and zero mutation favour M; `lambda↓0` and `N→∞` need not commute. A fully discrete cure atom can make ordinary means infinite, beyond the hybrid/coupling theorem. Hill steepness has no effect at fixed realised killing and admissible controls. A fixed lifetime-budget theorem is unproved. Stronger competition need not improve containment: in #66's saturated `c=L` calculation it reduces establishment and limiting advantage. **Mutation-supply monotonicity is not competition monotonicity.** |
| [#62: timing/sampling](metastatic-seeding-timing/) | Passenger counts, a biopsy endpoint substituted for clonal MRCA, or an uncalibrated clock do not identify migration. Separation fails at `ell=2h`. Positive average visibility does not ensure small robust error; worst-case invisibility is not inevitable in every population. The 17-region illustration is a necessary floor, not a sufficient design. Empirical moments, linked regions and issuance-selected cases invalidate naive guarantees. Genotype support does not count transfers, founding cells or extinct lineages. **Retain lag, sampling, provenance and finite-pilot assumptions.** |
| [#64: incidence](multistage-incidence-identifiability/) | Expansion alone cannot equate different finite transition counts under exact direct observation; unknown clone-level diagnosis delays can. Arbitrarily fast extra converting steps obstruct uniform finite-exposure recovery; any positive calibration allowance permits finite-bin confounding. Mean-count spectra do not classify finite-rate nonlinear hazards; asymmetric precursors and person-level diagnosis lags do not obey the converter-filter identity. **Exact onset identification is not stable inference from adult rounded rates.** |

| [#68: clonal haematopoiesis progression](clonal-hematopoiesis-progression/) | For every integer `n>=1`, positive `q,T` and `c>0`, the same founder count and net fitness `s=0` but `b=d=0` versus `b=d=c` produce different event hazards; hidden advanced stages similarly break fitness-only sufficiency. Exact cumulative hazard is linear in founder count, but risk is nonlinear and replacing stochastic clone exposure by its mean at finite `q` has a nonzero, explicitly bounded error. VAF-to-stem-cell count and a terminal mark-to-diagnosis link require separate assumptions. **Neither marginal Cox HR nor pooled VAF identifies the matched conditional progression law.** |

## Claim/status index

**Status key.** “Proved within model” refers to written assumptions and proof,
including declared computer-assisted components. “Known prior result” credits
established ingredients; “conjecture” is unproved. Empirical **consistent**,
**falsified** and **untestable** apply only to the specified observable/criterion.
#67's **failed predictive criterion** remains distinct from a calibrated
statistical rejection. Priority of candidate local contributions remains **unestablished**.
Symbols are local to their note; identical letters need not mean the same thing.

| Note | Precise central claim and quantifiers | Status as it now stands | Honest novelty level | Test file | Theory–data and dependency links |
| --- | --- | --- | --- | --- | --- |
| [Clonal-tree identifiability](clonal-tree-identifiability/) (#52) | For **each known target** rooted labelled tree in the exact perfect-phylogeny, single-malignant-founder model, some strictly positive rational designed mixture identifies it: **one row for a rooted path, two for branching**, with a closed entrywise perturbation-radius certificate. Arbitrary finite samples cannot guarantee uniqueness for `n>=3`, with open ambiguity sets despite positive usages and maximal rank. | **Proved within model.** | Narrow constructive synthesis. Sum conditions, one-row obstruction, tree-order dimension and counting ingredients are **known prior results**; no first-identifiability or achievable-biopsy-count claim. | [test_research_clonal_tree_identifiability.py](../tests/test_research_clonal_tree_identifiability.py) | [#62](metastatic-seeding-timing/) separately assumes a resolved genealogy; it does not test this designed-mixture theorem. No merged data test. |
| [Adaptive therapy / mutation supply](adaptive-therapy-mutation-supply/) (#53) | Two-individual-cell CTMC, initially `(S,R)=(2,0)`, per-time conversion, ideal containment C versus finite sensitive clearance M, progression at total burden `>2`: both **uncensored** means are infinite. For independent `Exp(1)` censoring and `(a,b,d,k,z)=(100,1,3,2,1)`, for `mu>0`, capped-time C-minus-M has the sign of `18mu^3−185mu^2−517mu+400`. Its positive roots are `alpha≈0.63735`, `beta≈12.44262`: C wins for `0<mu<alpha` or `mu>beta`, M between; ties at the exact roots and zero. | **Proved within model.** | Low, model-specific algebraic counterexample. Mutation supply versus suppression and objective dependence are prior themes. M is not proved a common bounded-control optimum. | [test_research_adaptive_therapy_mutation_supply.py](../tests/test_research_adaptive_therapy_mutation_supply.py) | [#61](containment-at-tumour-scale/) changes scale, mutation, controls and endpoint; [#66](containment-evidence/) confronts that later model. |
| [Combination / cross-resistance](combination-cross-resistance/) (#54) | Specified four-type branching, fixed finite **sensitive-only** founders, positive `b,mu,nu,k`, `d>=0`, equal per-drug exposure, fixed finite A-then-B period count `m`, `0<w<=1/2`, `0<=p<1`, `c>=1`: A-first improves extinction **by every sufficiently short course end** iff `c>(1-w)/w`; combination wins at equality at fourth order. A bounded-`c` local finite-course curve exists; for `w=1/2,m=1` its `T^3` coefficient has `p`-derivative `−mu*k^2/1440`. At `p=1`, sensitive-only starts tie. | **Proved within model; sixth-order identity computer-assisted.** | Narrow local coefficients/threshold refinement. Four-type branching, cross-resistance and sequencing advantages are **known prior results**. | [test_research_combination_cross_resistance.py](../tests/test_research_combination_cross_resistance.py) | No merged empirical test of its joint exposure, founder and killing assumptions. |
| [ecDNA segregation / therapy](ecdna-segregation-therapy/) (#55) | Independent cells, complementary `Bin(2k,1/2)` sisters, no reacquisition, uniformly bounded positive division rates, division-before-death `p_k=u*s^k`, `0<u<=1`, `0<s<1`: **every finite founder multiset containing a positive-copy cell** loses all positive lineages almost surely iff `2u(1−sqrt(1−s))<=1`, including equality. Monomorphic positive founders have pooled **genealogical-generation-two**, all-cell Fano `(3+s)/(2(1+s))`; arbitrary finite mixtures containing a positive-copy founder share the asymptotic pooled Fano limit. | **Proved within model.** Positive-lineage loss, not tumour extinction. | Narrow threshold/critical/Fano package; binomial segregation, branching and neutral moment laws are prior results. | [test_research_ecdna_segregation_therapy.py](../tests/test_research_ecdna_segregation_therapy.py) | Linear gap **closed by [#60](ecdna-killing-laws-and-clustering/)**. [#63](prediction-confrontation/)/[#65](ecdna-segregation-clustering-data/) audit observability/segregation; neither validates this killing law. |
| [Induced versus selected resistance](induced-vs-selected-resistance/) (#57) | Independent Poisson pools; unknown **common positive clone-size law**, finite second moment and non-singleton mass; singleton post-exposure acquisition; quiescent holidays; calibrated readout; positive reversible component and cellular switching `u>=0,v>0`: complete endpoint laws at **`0,H,2H`, known `H>0`**, identify all three intensities, both switching rates and the clone law. Three arms are minimal for a guarantee over the **whole stated class**, not every point; the unrestricted clone law needs the full baseline law. | **Proved within model; empirical application untestable in #63.** One-way-loss two-arm result adds third-moment/tripleton assumptions. | Narrow unknown-common-kernel design synthesis. Fluctuation tests, cumulants, switching and decompounding are known. | [test_research_induced_vs_selected_resistance.py](../tests/test_research_induced_vs_selected_resistance.py) | [#63 R1/R2](prediction-confrontation/) lacks matched holiday laws/conditioned covariance. Job 1 uses a different assay, not this holiday design. |
| [ctDNA surveillance schedule](ctdna-surveillance-schedule/) (#58) | Fixed residual `0<B*<D<C`, perfect threshold, `n>=1` dates in `(0,H]`, `H>0`, separate progression-scoring horizon `T in (0,infinity]`, and continuous positive log-rate density vanishing at both tails, strictly increasing/decreasing about one mode: the unique optimal ordered capture schedule is geometric, ratio `q=b/a`, with `a=log(D/B*)`, `b=log(C/B*)`, placement floor `max(a/H,b/T)`. A separate conditional regret bound handles bounded residual uncertainty and the specified finite-molecule assay; it does not preserve exact geometry. | **Proved within specified models; empirical optimum uncalibrated/untestable in #63.** | Narrow analytic solution/transfer synthesis. Surveillance optimization, highest-density rearrangement and submodular bounds are prior work. | [test_research_ctdna_surveillance_schedule.py](../tests/test_research_ctdna_surveillance_schedule.py) | [#63 C1–C3](prediction-confrontation/) audits capture/lead-time bridges; parity-thinning forecast lacks a matching published outcome. |
| [Immune escape / neoantigen clonality](immune-escape-neoantigen-clonality/) (#59) | Under A1–A8, independent branching, equal antigen weights, **linear unbounded killing**, all rates fixed, `h>0`, null-cell establishment `sigma>0`, actual source exponent `rho=b−d−u*kappa/(ell+kappa)−h<0`, no initial null cells: successful-null escape obeys `k*p_k -> sigma*h/kappa`. For sequences of finite founder populations with **minimum load tending to infinity** and `sum_i 1/k_i -> L` in `[0,infinity]`, total extinction tends to `exp(−sigma*h*L/kappa)`. Equal load multisets imply equal outcomes regardless of antigen identities. | **Proved within model.** Functional-load scaling empirically **untestable** with #63/#67 inputs. | Narrow fixed-positive-escape asymptotic/harmonic limit. Branching/rescue methods and biological clonal-target association are prior work. | [test_research_immune_escape_neoantigen_clonality.py](../tests/test_research_immune_escape_neoantigen_clonality.py) | [#63](prediction-confrontation/) tests bridges; [#67](immune-escape-pancancer-test/) tests an extra ecological surrogate and finds no gain. |
| [ecDNA killing laws / clustering](ecdna-killing-laws-and-clustering/) (#60) | Constant division `b>0`, linear whole-cell killing `d_k=gamma*k`, complementary unbiased binomial sisters, no reacquisition/copy cap: every finite positive-founder population persists with positive probability iff `gamma/b<c*`, extinct almost surely for `gamma/b>=c*`, where **`0.378474994554<c*<0.378474994555`**, with all-type tail certificates. Specified clustered/asymmetric families obey `c*<=(1−R)/2`. Identical full one-copy generation-two mean/Fano dose curves can have different thresholds. | **Proved within model. Interior clustering-threshold monotonicity: conjecture.** | Narrow killing-law classification, certificates and identifiability obstruction; beta-binomial segregation formulas are prior work. | [test_research_ecdna_killing_laws_and_clustering.py](../tests/test_research_ecdna_killing_laws_and_clustering.py) | Extends [#55](ecdna-segregation-therapy/); [#63](prediction-confrontation/) cannot test the threshold; [#65](ecdna-segregation-clustering-data/) transports published kernels conditionally. |
| [Containment at tumour scale](containment-at-tumour-scale/) (#61) | Under A1–A6's **hybrid** (deterministic sensitive bulk; stochastic resistant founders switched to deterministic density at `ceil(N^a)`, `0<a<1/2`; no initial resistance), shared bounded Hill killing, per-division `mu=lambda/N`, fixed `lambda>0`, course `L log N` and common exposure ceiling: C holds sensitive density at `s0`, M uses maximum dose, both stop after the course, and T is first total-burden passage to density `B=s0+eps`, including rebound. Their ordinary means satisfy `E[T_C]/log N -> c`, `E[T_M]/log N -> d+(J−d)exp(−lambda*A)`. Positive competition plus the theorem's admissibility conditions give `d<c<J` and one limiting crossing `log((J−d)/(c−d))/A`. C wins above; a strict limit gives eventual finite-`N` ordering, not a finite-size boundary. | **Proved within hybrid/scaling.** Near-boundary finite-size/biological transport unproved. | Model-specific progression-limit/coupling refinement; ecological trade-off, bounded dosing, branching establishment and thinning are prior work. | [test_research_containment_at_tumour_scale.py](../tests/test_research_containment_at_tumour_scale.py) | Revises [#53](adaptive-therapy-mutation-supply/)'s question. [#66](containment-evidence/) finds no joint aggregate identification of the parameters. |
| [Metastatic seeding timing](metastatic-seeding-timing/) (#62) | In the **conditional witness/clock experiment on witness capture**, calibrated positive calendar mutation rate and known reference/endpoint permit arbitrarily small uniform early/late error for hypotheses separated by `2h` **iff `ell<2h`**, where `ell` bounds ancestor-to-migration lag. For every visibility law `F` on `[0,1]` with population moments `m=E[F]>0,s2=E[F^2]`, sharp miss bounds permit a **maximum** invisible fraction `1−m^2/s2`; finite-pilot/distinct-region laws are separate. Founder-genotype support does not identify all migration events. | **Proved within model.** Useful biological lag/coverage: **conjecture**. Broad-witness submodel **conditionally falsified**; sharper sampling predictions **untested**. | Conditional design synthesis with classical clock, moment and occupancy ingredients; no first timing/sampling warning. | [test_research_metastatic_seeding_timing.py](../tests/test_research_metastatic_seeding_timing.py) | Its own TRACERx/Hessey checks; [#52](clonal-tree-identifiability/) concerns separate genealogy identification. |
| [Prediction confrontation](prediction-confrontation/) (#63) | Aggregate audit with recorded protocol/amendments for five theory notes yields **14 verdicts**: only I2's weaker clonal-burden response association is consistent; **13 are untestable** with located aggregates/controls. No unconditional corrected core prediction is falsified or validated. Exploratory ecDNA rejection-frequency evidence depends on source-test calibration/independence. | **Empirical: consistent (I2 only); untestable (other 13).** | Evidence/eligibility and arithmetic audit; no new biology or independent replication. | [test_research_prediction_confrontation.py](../tests/test_research_prediction_confrontation.py) | Bridges from [#55](ecdna-segregation-therapy/), [#57](induced-vs-selected-resistance/), [#58](ctdna-surveillance-schedule/), [#59](immune-escape-neoantigen-clonality/), [#60](ecdna-killing-laws-and-clustering/); followed by [#65](ecdna-segregation-clustering-data/)/[#67](immune-escape-pancancer-test/). |
| [Multistage incidence identifiability](multistage-incidence-identifiability/) (#64) | Constant positive finite initiation/mutation rates, clean birth, independent clones and exact direct hazard give `k=1+ord_0(h)`, even with intermediate expansion. A nondividing/lossless **converting** precursor exactly aliases a clone-level exponential diagnosis lag. Unbounded fast extra rates prevent uniform finite-Poisson-exposure recovery; with free amplitude and positive finite-bin means, **any `rho>1`** calibration allowance permits a sufficiently fast extra stage with identical means. SEER rounded rates identify no biological `k`; female breast fails the homogeneous monotone-mean condition. | **Proved within model; TSCE equivalence is known prior work. Empirical: descriptive shape falsification; biological count untestable.** | Narrow observation/finite-bin/calibration synthesis. Onset powers, TSCE invariants, spectral symmetries and projective-distance formula have prior foundations. | [test_research_multistage_incidence_identifiability.py](../tests/test_research_multistage_incidence_identifiability.py) | Its own [44-rate SEER provenance](multistage-incidence-identifiability/data/PROVENANCE.md); distinct from #52's tree and #62's migration estimands. |
| [ecDNA segregation / clustering data](ecdna-segregation-clustering-data/) (#65) | Transporting **nine published baseline ecDNA clone MLEs from one COLO320DM system** through #60's full asymmetric beta-binomial linear-killing model gives `R=0.089717–0.148457` and thresholds **11.19%–18.76% below binomial**, conditionally. Bounds cover rounding/numerics, not unknown sampling/model error. Uniform model envelope `L(R_U)<=c*<=(1−R_L)/2` holds for any fixed allocation-mixture law in a justified `R` interval; unordered counts have exact bias/dispersion aliases. | **Model envelope proved. Empirical: conditional source-fit implications; new dispersion inference and empirical eradication-threshold validation untestable.** | Narrow transport/envelope and observation audit. Unequal inheritance/fitted segregation are prior results; no empirical refit or calibrated confidence distribution. | [test_research_ecdna_segregation_clustering_data.py](../tests/test_research_ecdna_segregation_clustering_data.py) | [#55](ecdna-segregation-therapy/)/[#60](ecdna-killing-laws-and-clustering/) theory; revisits [#63](prediction-confrontation/), correcting HSR exclusions and retaining counterevidence. |
| [Containment evidence](containment-evidence/) (#66) | Audit of **69 source/study records, 15 numerically included reports, 179 aggregate rows, 44 contrasts** (not independent-study counts): no report jointly identifies `d,c,J,A,lambda`. **P1/P2/P2a/P3/P4/P5 all untestable**, zero eligible reports/formal mechanistic tests, no pooled effects. P6a is consistent with some benefit; P6b is falsified only under its broad outcome/any-reported-interval reading, by COIN's **80%-CI OS** comparison; narrower post-extraction sensitivities remain untestable. | **Empirical: mostly untestable; limited descriptive consistency and qualified extension falsification.** | Provenance/identifiability audit, not a new efficacy estimate, complete systematic review or mechanistic discovery. | [test_research_containment_evidence.py](../tests/test_research_containment_evidence.py) | Confronts [#61](containment-at-tumour-scale/), with a correction to universal competition-strength benefit; [#53](adaptive-therapy-mutation-supply/) has a different model/endpoint. |
| [Immune escape pan-cancer test](immune-escape-pancancer-test/) (#67) | Extra surrogate `p=[A0(1−L)+A1*L] exp[−B/(T*f)]`, `0<=A1<=A0<=1,B>=0`, uses median TMB `T`, clonal fraction `f`, ordinary-HLA-LOH prevalence `L` and an unvalidated response bridge. Grouped held-out reciprocal-TMB RMSE **13.173 pp vs log-TMB 12.084 pp**, 21 rows/20 parents, fails the improvement criterion. Fixed clonal extension worsens matched prediction; free exponent is negative. Full coherent **nine-parent** panel fails the frozen ten-parent gate; full-model fits are exploratory. | **Empirical: failed predictive criterion; contrary clonal-proxy direction. Full gated test and mechanistic scaling untestable.** | Repository-specific surrogate test with preserved negatives; no demonstrated biomarker, first ecological analysis or independent outcome replication. | [test_research_immune_escape_pancancer_test.py](../tests/test_research_immune_escape_pancancer_test.py) | Tests an extra observation model built on [#59](immune-escape-neoantigen-clonality/), after [#63](prediction-confrontation/); does not refute the theorem. |

| [Clonal haematopoiesis progression](clonal-hematopoiesis-progression/) (#68) | For every integer `n>=0`, finite `b,d,q,T>=0` in the independent marked birth–death model, `S_n(T)=(E exp(-q A_T))^n` with `A_T=∫₀ᵀZ(t)dt`; thus `H_n=-n log E exp(-q A_T)`, not generally `nqK(s,T)` at finite `q`, `s=b−d`. The exact bound is `0<=nqK(s,T)−H_n<=nq² E[A_T²]/2`; for `n→∞,q=κ/n`, `H_n→κK(s,T)` at fixed rates/horizon. For matched initial count and shared nonnegative completion measure, the rare-hit cumulative-hazard ratio at `s_high>=s_low` lies in `[1,exp((s_high−s_low)T)]`. Unknown turnover, advanced states, observation or diagnosis mapping defeat inference from net fitness/VAF alone. | **Proved within specified branching/observation models.** Main conditional biological target **empirically untestable**; seven-gene ten-year substituted baseline-fitness **proxy falsified only by the declared descriptive interval-box rule** (U2AF1). The 20-year per-gene overlap sensitivity is consistent, while a post-hoc common-intercept check remains infeasible. | **Narrow original-to-repository synthesis/audit, global novelty unestablished.** CH growth/AML associations, CHRS, multiple hits and branching machinery are **known prior results**; not a discovered mechanism or clinical risk predictor. | [test_research_clonal_hematopoiesis_progression.py](../tests/test_research_clonal_hematopoiesis_progression.py) | Its own [source/endpoint audit](clonal-hematopoiesis-progression/SOURCE_AUDIT.md) and [frozen prespecification](clonal-hematopoiesis-progression/PRESPEC.md); related but distinct identifiability limits in [#57](induced-vs-selected-resistance/) and [#64](multistage-incidence-identifiability/), and aggregate testing limits in [#63](prediction-confrontation/). |

## Five next jobs, ranked by expected scientific value

This is a prioritization judgment, not a forecast of successful discovery.
Ranking weighs mechanistic discrimination, independent observations, available
measurement, and the chance a negative answer changes the programme. All jobs
are nonclinical and use existing public cell-line data or published aggregates.
Freeze question, estimand, exclusions, sensitivity region and held-out scoring
before fitting; retain null and unavailable outcomes.

1. **Highest: does resistance-origin inference survive realistic outgrowth?**
   Use **Russo et al. 2022, Mendeley Data V2,
   [DOI 10.17632/mvfm7hs9kw.2](https://data.mendeley.com/datasets/mvfm7hs9kw/2)**,
   `Russo_et_al_v2.zip`: DiFi/WiDr early/late colony assays and population
   measurements, with the [paper's Fig. 3e and Supplementary Tables 1, 2, 4,
   5](https://www.nature.com/articles/s41588-022-01105-z).
   Within a prespecified competing-model family, ask whether held-out colony
   counts require increased post-exposure acquisition after allowing delayed
   outgrowth, measured persister turnover and unequal clone-size/readout
   kernels. Hold out entire cultures/biological replicates and keep independent
   dynamics calibration separate from the colony-count prediction target.
   A robust signal or predictive observational alias would sharpen
   [#57](induced-vs-selected-resistance/). This is not its quiescent `0,H,2H`
   experiment and cannot by itself prove drug causation.

2. **High: does the immune surrogate's Merkel-cell residual survive assay
   and viral-status stratification?** Use **JAVELIN Merkel 200 Part B,
   [D'Angelo et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8311489/),
   Fig. 1D, biomarker Results and Supplementary Fig. 3B**, and distinct
   **Part A [2020 report](https://pmc.ncbi.nlm.nih.gov/articles/PMC7239697/),
   response/biomarker Results and Figs. 1/3**. Compute which common assay-scale
   conversions, if any, make the frozen [#67](immune-escape-pancancer-test/)
   TMB-only curve compatible with each cohort's MCPyV-positive/negative ORR
   intervals and TMB medians; compare the compatibility regions. Chalmers
   mutations/Mb (including synonymous alterations) and JAVELIN NSSV/Mb differ:
   this is an observation-bridge sensitivity test, not calibrated external
   validation. Any exploratory viral-status extension must have at most one
   added parameter, fixed on Part A before scoring Part B. Keep first/later-line
   cohorts and TMB-/ORR-evaluable denominators separate; marginal summaries are
   not paired data and overlapping follow-ups are not replications. Median
   substitution remains an assumption; no functional-load or causal inference.

3. **Moderate: how fragile is the ecDNA fair-segregation rejection?** Use
   **[Lundberg et al. 2008](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003099),
   Results “Statistical distribution of DMs at mitosis”, Fig. 3E–G**:
   CHP-212 `2/100`, LA-N-5 `10/100`, MC-IX `1/50` nominal `0.01` test rejections,
   already in [lundberg.csv](ecdna-segregation-clustering-data/data/lundberg.csv).
   Map the joint miscalibration/dependence frontier that removes the
   Holm-adjusted inconsistency, using exact tails under a declared exchangeable
   count model and dependence-robust bounds. This correlation concerns test
   indicators, not ecDNA dispersion, and cannot be estimated from these counts.
   A quantified robustness frontier would extend
   [#65](ecdna-segregation-clustering-data/)'s two dependence endpoints;
   it would not supply a new segregation estimate.

4. **Moderate: does the incidence-model failure persist within published
   demographic strata?** Use **SEER CSR 1975–2017, incidence period 2013–2017,
   [Table 4.11](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_04_table.11.pdf)**,
   age-specific female All Races/White/Black columns at ages 30–84, plus
   **Tables [6.10](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_06_table.10.pdf),
   [15.10](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_15_table.10.pdf),
   [22.7](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_22_table.07.pdf)**
   for prespecified site/sex contrasts. Apply
   [#64](multistage-incidence-identifiability/)'s rounding-aware monotonicity
   and calibration frontier within strata; compare which ages force model
   relaxation. This asks whether pooled shape failure survives the available
   stratification; it cannot identify its cause, supply missing denominators
   or estimate drivers.

5. **Moderate: how robustly do terminal-cause distributions differ across
   containment schedules?** Use **Seyedi et al. 2024,
   [Cancers 16:257](https://www.mdpi.com/2072-6694/16/2/257),
   [supplement](https://www.mdpi.com/article/10.3390/cancers16020257/s1),
   `cancers-2752834-supplementary.pdf`, Tables S2–S4 (pp. 32–33)**: aggregate
   outcomes attributed to toxicity, cancer or other causes across capecitabine,
   gemcitabine and combination schedules. Table S1 records attribution reasons;
   Table S5 supplies dose means. Prespecify exact small-sample cause-composition
   contrasts and reclassification bounds, preserving source-pooled groups and
   counting the shared combination control once. Determine whether the
   toxicity-attribution contrast survives the declared sensitivity range.
   This quantifies [#66](containment-evidence/)'s comparator warning. Unequal
   observation times and competing risks prevent interpretation as equal-horizon
   cancer mortality, cause-specific survival, causal attribution, competition
   or #61's mutation threshold.

The original index's laboratory/new-sampling proposals are not executable
with the present public inputs. A new 100-mother ecDNA experiment, matched
ctDNA schedule outcomes, calibrated metastatic lag/pairwise downsampling
release, or real-tumour `lambda/lambda*` estimation is not an existing-data
job. #60's solved linear threshold is also no longer an open task.

## Verification boundaries and reproducibility

Each linked test checks implementations and stated examples; read it alongside
the note. **Finite checks do not prove all-parameter claims.** #52 independently
enumerates trees through five vertices; #53/#54 compare rational resolvents,
population generators and PGFs, with a generic symbolic certificate for #54's
sixth-order component. #55 uses independent joint-population checks and rational Fano-series remainder
bounds; #60/#65 additionally use independent segregation kernels and explicit
head/tail certificates. A finite copy cap alone cannot certify extinction. #57 checks exact aliases as well as inverses. #58's Decimal
references use high precision, not outward-rounded interval certification.
#59 separates numerical brackets from its all-load proof. #61's stopped-founder
coupling is not a full-discrete progression-time theorem. #62 distinguishes
finite-pilot, finite-pool and conditional-clock guarantees. #64 separates
exact algebra, rounded-rate sensitivity and synthetic exposure calculations.
#68 checks the independent backward equation, stochastic/mean-exposure error
bounds, hidden-state counterexamples, all published gene-proxy LOO predictions,
held-out leakage, sensitivity scenarios and seeded simulations. These are
implementation checks, not validation of a clinical progression mechanism.

Empirical tests in #63/#65/#66/#67/#68 recompute source joins, exclusions,
denominators, uncertainty conventions and generated blocks. They do not
establish source truth, search completeness, undocumented sampling assumptions,
independent replication or biological validity. Respect the frozen protocols
and amendments linked inside the notes. Missing inputs stay missing; rounding
precision, parameter sensitivity, resampling and confidence are distinct.
Source-fit transport and post-hoc checks retain those labels. Unresolved
intervals are not ties; event/search caps do not prove extinction/impossibility.

Each note's bibliography is the priority ledger. Known results include sum
conditions/order-dimension bounds, fluctuation/branching machinery, ecological
trade-offs, surveillance optimization, clonal-target associations,
moment/occupancy bounds and multistage identifiability. Local proofs and
independent algorithms do not establish historical novelty. AI agreement and
passing CI are neither external scientific review nor clinical validation.

From the repository root, run the research-prefixed tests (also includes the
research-pilot infrastructure test):

```bash
python -B -m unittest discover -s tests -p 'test_research_*.py' -v
```

For one note, replace the `-p` pattern with its linked test **basename**.
These dedicated research contributions use Python's standard library and
offline inputs. Individual bounded workloads do not imply the whole programme
finishes within a minute. Check stored empirical reports without rewriting:

```bash
python -B research/prediction-confrontation/analysis.py --check
python -B research/ecdna-segregation-clustering-data/analysis.py --check-readme
python -B research/containment-evidence/analysis.py --check
python -B research/immune-escape-pancancer-test/analysis.py --check
python -B research/multistage-incidence-identifiability/analysis.py
python -B research/clonal-hematopoiesis-progression/analysis.py
python -B research/clonal-hematopoiesis-progression/verify.py
```

Repository gates use the existing dependency in
[requirements-dev.txt](../requirements-dev.txt); the existing
[PR workflow](../.github/workflows/validate.yml) discovers the tests:

```bash
python -B scripts/validate.py
python -B -m unittest discover -s tests -v
node --test tests/*.test.*
```

Follow [DATA_POLICY.md](../DATA_POLICY.md) and [AGENTS.md](../AGENTS.md): no
patient-level/controlled-access records, inferred identities or private
credentials. Retain authoritative sources, versions, locators, reuse terms
and negative evidence. Public availability is not redistribution permission;
repository licensing does not relicense third-party data. The programme
remains **nonclinical and unreviewed**.
