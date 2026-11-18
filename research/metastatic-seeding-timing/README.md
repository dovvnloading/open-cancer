# Metastatic seeding timing: conditional identification and sampling limits

**Nonclinical computational mathematics. No treatment, diagnostic, screening, or
patient-specific advice. Research status: unreviewed.** Proofs concern the models
below; exact arithmetic, simulations, AI critique, and CI do not establish their
biological validity or priority.

## Result and explicit pivot

Multi-region sequencing can date a common ancestor without dating the movement
of a cell. That distinction, sampling-induced apparent clonality, and ambiguity
between a cell cluster and repeated migrations are already established. We do
not claim their discovery or a universal solution of metastatic timing.

The contribution proposed here is a **conditional identification and experiment
design framework**. It joins a sharp separation condition for a calibrated
ancestral clock with sharp coverage/depth bounds given two population sampling
moments, a finite-pilot correction, a distinct-region counterpart, and explicit
limits on founder-clone counting. The central question is:

> Given a separately justified bound on the time from a sampled primary
> ancestor to migration, and only one- and two-region detection summaries, what
> timing error can an experiment actually guarantee?

The principal conclusions are:

* In the bounded-lag clock experiment, early and late hypotheses separated by
  `2h` are distinguishable uniformly exactly when the maximum unobserved lag
  `ell` is smaller than `2h`, provided the clock rate is positive and calibrated.
* If a diagnostic primary witness occurs in a fraction `F` of potential
  regions, two moments `m = E[F]`, `s2 = E[F^2]` give exact best and worst
  possible miss probabilities for every number of regions and every modeled
  depth. They permit a worst-case invisible fraction `1 - m^2/s2`. This is an
  information limit of the two-moment calibration, not an assertion that this
  fraction of real tumours has invisible ancestors.
* Increasing depth cannot compensate for a witness absent from the sampled
  spatial population. An exact budget certificate can therefore say
  **unattainable**, rather than always returning a larger sequencing budget.
* Counting detectable, surviving, genetically distinguishable founder lineages
  is different from counting migration events, founding cells, or extinct
  lineages. We supply recovery laws for the first estimand and counterexamples
  for the others.

These are general-parameter results. The birth-death model has no fixed cell
cap; the coverage bounds hold for every probability distribution on `[0,1]`;
the finite-pool theorem holds for every integer pool size. Small rational
fixtures in the code check formulas and boundary cases; they are not the
domain of the theorems.

**Novelty assessment:** the proposed joint application and certificate structure
may be new; priority has not been established. Poisson marking, classical
two-moment inequalities, Hermite interpolation, occupancy inclusion-exclusion,
and finite-dimensional moment optimization are standard. None is presented as
a newly discovered mathematical technique. The useful extension is the explicit
combination, its attainable/unattainable design boundary, and its qualifications
for metastasis observations. A contribution to biological inference would still
require independent validation of the witness and lag assumptions.

## 1. Literature audit before model selection

Search date: **2026-11-17**, using primary papers, their author manuscripts, and
publisher records. Searches included the named papers, `metastasis seeding time
identifiability`, `private shared mutations sampling`, `metastatic two-moment
sampling`, `Gunnarsson branching process site frequency spectrum`, and recent
2024–2026 studies. This is a targeted audit, not a systematic review or proof
that no matching result exists. Only public aggregate facts are used below;
no patient-level sequence or clinical records were downloaded.

| Prior work | Established result relevant here | Boundary on our claim |
| --- | --- | --- |
| [Hu et al. 2019, Nature Genetics](https://doi.org/10.1038/s41588-019-0423-x) | SCIMET combines tumour evolution, mutation accumulation, and sampling to infer early colorectal dissemination; 19/23 evaluable primary/metastasis pairs were inferred to seed below a primary size of `10^8` cells. | Inferring timing with passengers and growth assumptions is established; those are model-based estimates, not directly observed migration times. |
| [Hu et al. 2020, Nature Genetics](https://doi.org/10.1038/s41588-020-0628-z) | Private clonal counts and growth assumptions estimate timing in 457 samples from 136 patients. The paper explicitly discusses inability to distinguish a multiclonal cluster from separate migrations. | Neither a private-count clock nor cluster/event ambiguity is new here. |
| [Gundem et al. 2015, Nature](https://doi.org/10.1038/nature14347) | Metastasis-to-metastasis spread and transfer of multiple clones in lethal prostate cancer. | A primary-only model is a restricted biological submodel. |
| [Turajlic et al. 2018, Cell / TRACERx Renal](https://doi.org/10.1016/j.cell.2018.03.057) | Distinct routes and evolutionary patterns in 100 patients, with 575 primary and 335 metastatic biopsies. | Anatomical destination, selection, and latency matter. |
| [Reiter et al. 2017, Treeomics](https://doi.org/10.1038/ncomms14114) | Correct tumour phylogenies need not identify anatomical sources or migration order; an early-diverging lineage may migrate late. | The unmarked-edge counterexample in Theorem 2 formalizes a known limitation. |
| [Reiter et al. 2018, Science](https://doi.org/10.1126/science.aat7171) | Functional driver homogeneity across 76 untreated metastases from 20 patients, with a mathematical explanation. | Shared drivers are not an independently calibrated neutral timing clock. |
| [Reiter et al. 2019, Nature Reviews Cancer](https://doi.org/10.1038/s41568-019-0185-x) | Reanalysis of genetic heterogeneity in untreated cancers. | Driver concordance and timing identification are different questions. This review is context, not the source of our technical results. |
| [Heyde, Reiter, Naxerova and Nowak 2019, PNAS](https://doi.org/10.1073/pnas.1819408116) | Immigration/birth-death theory of transferred diversity, polyclonality, and surviving founding cells. | A new generic founder-count or continuous-seeding model is not claimed. |
| [El-Kebir et al. 2018, MACHINA](https://doi.org/10.1038/s41588-018-0106-z) | Parsimonious migration histories from clone trees; multiple compatible anatomical histories. | Parsimony is not proof of the actual migration history. |
| [Sun and Nikolakopoulos 2021, PLoS Computational Biology](https://doi.org/10.1371/journal.pcbi.1008838) | Decomposes measurable private mutations by ancestor detectability and derives distributions; selection can destroy monotonicity with seeding time. | A claim that ordinary private counts always increase with migration time is false. |
| [Al Bakir et al. 2023, TRACERx](https://doi.org/10.1038/s41586-023-05729-x) | 94/126 late and 32/126 early **divergence** cases; single-region downsampling misclassified 75/90 late cases as early. | Divergence is explicitly distinct from cell migration; this is our aggregate sampling benchmark, not a clinical seeding label. |
| [Hong, Shpak and Townsend 2015](https://doi.org/10.1158/0008-5472.CAN-15-1889), [Werner et al. 2017](https://doi.org/10.1038/srep44991), [Opasic et al. 2019](https://doi.org/10.1186/s12885-019-5597-1) | Spatial undersampling and apparent clonality; mathematical multi-region sample calculations. | Broad sampling, geometric miss laws, and the general depth-versus-space warning are established. |
| [Durrett 2013](https://doi.org/10.1214/11-AAP824), [Williams et al. 2016](https://doi.org/10.1038/ng.3489), [Werner et al. 2020](https://doi.org/10.1038/s41467-020-14844-6) | Branching-process passenger spectra and inference of evolutionary/division histories. | We do not claim a new neutral SFS limit or equate calendar time with division count. |
| [Gunnarsson, Leder and Foo 2021](https://doi.org/10.1016/j.tpb.2021.09.004), [Gunnarsson, Leder and Zhang 2025](https://doi.org/10.1016/j.spa.2025.104565) | Exact expected SFS; almost-sure limits at fixed times and population hitting sizes; parameter estimation. | Strong branching-process spectrum theory already exists. |
| [Avanzini and Antal 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6871767/) | Size-dependent seeding with independent metastatic birth-death growth and recurrence-time distributions. | The growth/seeding part of Theorem 1 is background, not a new process. |
| [Angaji et al. 2024, eLife](https://doi.org/10.7554/eLife.95338) | Dense sampling and explicit birth-death pairwise mutation likelihood; Appendix 6 examines lineage-sampling bias in division clocks. | Our mutation clock is explicitly per calendar branch time; no unproved division-clock substitution is used. |
| [Blohmer, Cheek, Hung et al. 2025, Nature Genetics](https://doi.org/10.1038/s41588-025-02078-5) | Polyguanine-based lineage division counts in 505 samples from 37 patients; metastatic divergence is an earliest possible seeding time. | Even a refined division clock is not automatically a migration clock. |
| [Koyyalagunta, Ganesh and Morris, Metient, online 2025 / issue 2026](https://doi.org/10.1038/s41592-025-02924-8) | Multiobjective reconstruction and alternative migration histories in 167 patients / 479 tumours in the final publication. | We do not claim a new migration-history reconstruction algorithm; final-publication denominators supersede differing preprint counts. |
| [Chen et al. 2026](https://doi.org/10.1093/molbev/msag008) | Dense three-dimensional mapping and model-based early continuous dissemination in two colorectal cases. | Dense empirical mapping is valuable but does not establish a universal identification theorem. |
| [Hessey et al. 2026, Nature](https://doi.org/10.1038/s41586-026-10428-4) | In 501 longitudinal samples from 24 patients, 156/258 lesions were inferred seeded by other metastases, 98/258 by the primary, and 4/258 by both. | Applying our primary-only generator as a universal model contradicts published aggregate evidence. |
| [Pregnall et al. 2026](https://doi.org/10.1038/s41698-026-01291-7) | Contemporary application of timing models reports ambiguous timing for 5/15 primary/metastasis pairs. | Clock/model confounding is an active limitation, not a solved empirical problem. |

The abstract clinical question remains open in full generality, while many of
its simpler formulations and impossibility observations are already addressed.
We therefore pivot explicitly to the narrower conditional experiment above.
No result below estimates a seeding date for any published individual.

## 2. Model, observations, and estimands

### Growth and seeding assumptions G1–G6

1. **G1:** `N(t)` is a continuous-time linear birth-death primary, starting from
   an integer `n0 >= 1`, with per-cell birth `bP > dP >= 0` and net growth
   `rP = bP - dP > 0`. Birth replaces one cell by two; death removes one.
   `E[N(t)] = n0 exp(rP t)`. This is exponential growth in expectation, not a
   deterministic replacement for the random population. Extinct primaries are
   retained unless an explicitly stated conditioning event excludes them.
2. **G2:** Conditional on the primary path, seeds form a Poisson process with
   intensity `eta N(t)`, `eta >= 0`. Seeding buds a cell without depleting the
   primary, and chooses a primary cell uniformly. Independent Poisson clocks
   are equivalent to this rule. A model in which migration removes cells and
   changes the observed primary path would need a different conditional law.
3. **G3:** Each seed starts one independent metastatic birth-death process with
   rates `bM, dM >= 0`, independent of the primary, other lesions, and the
   passenger process. In this base model, seeds found separate sites; there is
   no secondary seeding, dormancy, treatment, selection, merging, or migration
   back to the primary. Theorem 6 separately treats multiple founder genotypes
   per observed lesion as an observation-model extension.
4. **G4:** A surviving seeded lineage means at least one descendant at the
   ascertainment time `T`. It need not be radiologically detectable. A known
   size-ascertainment rule can replace survival by its own inclusion probability.
   A statement about surviving seeded lineages is never a count of all human
   metastases or all disseminated cells.
5. **G5:** The first time the primary reaches a specified detection size `ND`
   is `C = inf{t: N(t) >= ND}`. The conditional timing experiment uses an
   independently known value `c` of this threshold-crossing time (or a known
   reference time stipulated by an experiment). **Diagnosis time need not equal
   first detectability.** Neither `C` nor the physical threshold is identified
   by endpoint mutation counts alone. We do not set `C = log(ND)/rP`.
6. **G6:** Only finite trees and finite observation horizons are used. Linear
   birth-death processes are nonexplosive at every finite horizon. All stated
   identities hold for arbitrary finite populations and rates in their stated
   domains, without a small-tumour approximation.

### Passenger and observation assumptions O1–O6

1. **O1:** On each genealogical branch, neutral infinite-sites passengers occur
   as a Poisson process with known rate `theta = E * mu` per unit calendar
   branch time. `mu >= 0` is a calibrated rate per genomic opportunity and
   `E >= 0` is callable genomic exposure. Separate branches and disjoint clock
   loci have independent increments. No recurrence, mutation loss, copy-number
   distortion, or location-dependent clock is allowed. Rates per division
   require another model and are not silently substituted for `mu`.
2. **O2:** The ideal clock experiment observes mutations assigned to disjoint
   edges or disjoint unions of edges of a correctly reconstructed sampled
   genealogy. Conditional statements treat that genealogy and its sampling
   scheme as specified. VAF deconvolution/topology errors are not absorbed into
   Poisson noise. One diploid heterozygous passenger at cancer-cell fraction
   `f` has VAF `f/2` only under the stated purity/copy-number idealization.
3. **O3:** On a chosen metastatic lineage, `a` is its latest common ancestor
   with a sampled primary witness; `s` is migration time; `B >= s` is a known
   terminal time of the counted path. For a single sampled terminal cell,
   `B` is its sampling time. For metastasis-private **clonal bulk** mutations,
   the correct endpoint is the sampled metastatic MRCA, not biopsy time.
   Its calendar age would need separate calibration. The clock observation is
   the path from `a` to `B`, not an arbitrary private mutation total.
4. **O4:** A prespecified diagnostic primary witness, if captured and correctly
   resolved, guarantees `0 <= s - a <= ell`, for a justified constant `ell`.
   This is an additional biological calibration assumption. It does **not**
   follow from exponential growth, neutral passengers, or simply detecting
   any shared mutation. Unobserved/extinct ancestors can violate it.
5. **O5:** Marker loci and their sequencing are independent of the clock loci
   conditional on the history. Anchor selection is prespecified or uses only
   these separate marker data. On capture, the clock panel is correctly called
   and assigned. Thus conditioning on detection does not select high or low
   clock counts. A library artifact, mapping error, wrong ancestry assignment,
   or false marker call violates this idealization.
6. **O6:** In the spatial experiment, an individual has a latent witness-positive
   region fraction `F in [0,1]`. Fresh regions are independent conditional on
   `F`. A positive region has a common known marker VAF `v`; at integer depth
   `D`, its reads are independent Bernoulli(`v`). One mutant read detects the
   witness, with no false positives. Thus `cD = 1 - (1-v)^D`. Higher read
   thresholds are possible by substituting the appropriate binomial tail.
   Cross-depth transfer assumes this common `v`; heterogeneous/correlated VAF
   cannot be inferred from the two moments of `F` alone.

Conditional on a fixed primary/lesion history, `F` describes the sampling
experiment, not a new evolution mechanism. Across histories it can be random.
The sharp moment theorems deliberately allow **all** distributions of `F`
with the prescribed moments. We do not claim their extremizers arise under
every fixed choice of birth-death parameters. Their robustness comes from
enlarging that observation class.

### Which quantities can be counted?

| Estimand | Meaning here | What the experiment can establish |
| --- | --- | --- |
| Observed lesions | Anatomical sites already included in the study | Their recorded number; not the number of unsampled sites |
| Extant founder genotypes | Distinct, diagnostic, pre-seeding genotypes with surviving descendants | Support recovery under Theorem 6's visibility and provenance assumptions |
| Seeding cells / migration events | Individual transfers, including repeated identical genotypes | Not identified by genotype support alone |
| Extinct or never sampled founders | Transfers with no observed descendants | Not recoverable from endpoint genotype support |

## 3. Numbered results with complete proofs

### Theorem 1 — Exact seeding and passenger observation laws

**Assumptions:** G1–G4 and G6; O1–O2 for the mutation statements. No conditioning
on primary survival is required. Let `q(u) = P(Z(u) > 0 | Z(0)=1)` for a
metastatic birth-death process, and let `J(T)` count surviving seed lineages.
Then

$$
q(u)=
\begin{cases}
\displaystyle\frac{b_M-d_M}{b_M-d_M e^{-(b_M-d_M)u}},&b_M\ne d_M,\\[4pt]
\displaystyle\frac{1}{1+b_Mu},&b_M=d_M.
\end{cases}                                                   \tag{1}
$$

This includes `q(u)=1` when both rates are zero and `q(u)=exp(-dM*u)` for pure
death. Conditional on the entire primary path,

$$
J(T)\mid N\sim\operatorname{Pois}(\Lambda_N),\qquad
\Lambda_N=\eta\int_0^T N(t)q(T-t)\,dt.                         \tag{2}
$$

Consequently its unconditional probability generating function is
`E[exp((z-1)*Lambda_N)]`, its mean is
`eta*n0*integral_0^T exp(rP*t)*q(T-t) dt`, and
`Var(J) = E[Lambda_N] + Var(Lambda_N)`. Replacing the random primary by its
mean generally gives a wrong variance and a wrong probability of zero lesions.

For a fixed sampled genealogy `G`, edge `e` has calendar length `le` and
descendant leaf set `Ae`. For each nonempty observed leaf subset `A`, let
`YA` count mutations with exactly that presence pattern. Then

$$
Y_A\mid G \overset{\mathrm{ind}}{\sim}
\operatorname{Pois}\!\left(\theta\sum_{e:A_e=A}l_e\right),
\quad
\mathbb E\!\left[\prod_A z_A^{Y_A}\mid G\right]
=\exp\!\left[\theta\sum_e l_e(z_{A_e}-1)\right].               \tag{3}
$$

In particular, the clock path in O3 has

$$
X\mid a,B\sim\operatorname{Pois}(E\mu(B-a)).                  \tag{4}
$$

Shared mutations on a disjoint root-to-ancestor path of length `a` have law
`Pois(E*mu*a)` and are independent of `X`. Grouping edge patterns by descendant
count gives the corresponding conditional site-frequency spectrum. Shared
summaries that reuse edges are not independent.

**Proof.** For one initial metastatic cell, the first-event decomposition gives
`q' = (bM-dM)q - bM*q^2`, `q(0)=1`: a birth gives two independent descendant
processes and a death gives none. If `bM != dM`, solving the linear equation
for `1/q` gives (1); if rates agree, integrate `q'=-bM*q^2`. These also verify
the stated zero-rate limits. Conditional on `N`, seeding is an inhomogeneous
Poisson process. A seed at `t` survives with probability `q(T-t)`, independently
of all other seeds. Independent marking/thinning gives (2). Iterated
expectation and conditional variance give the displayed marginal formulas;
Tonelli applies because all integrands are nonnegative. The primary mean
follows from `d E[N]/dt = (bP-dP)E[N]`.

On each fixed edge, Poisson increments give independent counts with mean
`theta*le`. Infinite sites and no loss give each mutation the edge's descendant
set. Sums over disjoint groups of edges are independent Poisson variables,
proving (3), and products of their generating functions give the stated PGF.
The chosen path has total length `B-a`, giving (4). No statement replaces a
calendar branch length by a number of divisions. **QED.**

### Theorem 2 — What dates migration, exact confounding, and clock exposure

**Assumptions:** O1–O5, with `mu > 0`; known `B` and known detection reference
`c`; `h > 0` and `0 <= c-h < c+h <= B`. Seeding satisfies `0 <= s <= B`.
The nuisance lag `L=s-a` is any value in `[0,min(ell,s)]`. Thus this is a
conditional marked-lineage experiment, not the full marginal likelihood of a
known-parameter birth-death genealogy. Exposure `E>0` represents independent
genomic clock opportunity on that path.

**(a) Location-blind law.** If the only change is moving an unobserved migration
mark along the same unbranched ancestral edge, keeping its sampled genealogy
and clock lengths fixed, every shared/private pattern in (3) has the same law.
If all admissible placements on that edge form `[a,b]`, this is the sharp
conditional placement set. Early (`s<c`) is certified by `b<c`; late (`s>=c`)
by `a>=c`. An interval containing placements on both sides cannot be resolved
by more location-blind passengers. Endpoints are included as the closure of
admissible placements; conservative certificates use the strict early bound.

**(b) Sharp separated-hypothesis condition.** Consider

$$
H_E:s\le c-h,\qquad H_L:s\ge c+h.                            \tag{5}
$$

If `2h <= ell`, two feasible early/late histories have exactly the same law
of `X` at every exposure. The maximum of their two classification errors is
at least `1/2`, even allowing randomization. If `2h > ell`, put

$$
\lambda_E=\mu(B-c+h),\qquad
\lambda_L=\mu(B-c-h+\ell),\qquad \lambda_E>\lambda_L\ge0.      \tag{6}
$$

For the deterministic rule **early iff `X >= k`**, the exact worst errors are

$$
e_E(E,k)=P_{\operatorname{Pois}(E\lambda_E)}(X<k),\qquad
e_L(E,k)=P_{\operatorname{Pois}(E\lambda_L)}(X\ge k).            \tag{7}
$$

For `0<alpha<1/2`, the exact minimum **integer genomic exposure among these
deterministic threshold rules** is

$$
E_{\min}=\min\{n\ge1:\exists k\in\mathbb Z_{\ge1},\
e_E(n,k)\le\alpha\text{ and }e_L(n,k)\le\alpha\}.               \tag{8}
$$

Equation (8) is an explicit finite search with Poisson-tail certificates, not
a claim that deterministic thresholds dominate randomized tests at a discrete
boundary. Define `Delta = (sqrt(lambdaE)-sqrt(lambdaL))^2`. Uniform error at
most `alpha` is achievable for

$$
E\ge\frac{2\log(1/\alpha)}{\Delta},                          \tag{9}
$$

while **every** test with both errors at most `alpha` requires

$$
E\ge\frac{\log(1/[4\alpha(1-\alpha)])}{\Delta}.               \tag{10}
$$

For `lambdaL=0`, use threshold `k=1`: its errors are
`exp(-E*lambdaE)` and zero. A zero clock rate or zero exposure gives no timing
information. With unknown freely varying `mu`, even separated physical times
can have identical Poisson means.

**(c) Interval version.** Any `1-alpha` confidence interval `[aL,aU]` for `a`
on successful witness capture gives a seed interval
`[aL, min(B,aU+ell)]`, intersected with `[0,B]`. A strict early certificate
requires `aU+ell<c`; a late certificate requires `aL>=c`. Values between the
guard-band hypotheses are deliberately not covered by the binary guarantee.

**Proof.** For (a), (3) depends on edge lengths and descendant patterns, not the
location of an unobserved anatomical transition. An unobserved primary
continuation can be pruned without adding a sampled branch. Each admissible
mark gives the same clock law. Hence no statistic of these passengers can
distinguish the placements. This is conditional equality of observation laws;
different placements may have different probabilities under a fully specified
growth/seeding prior, and extra population-size or longitudinal data can matter.

For (b), if `ell>=2h`, choose `(sE,LE)=(c-h,0)` and
`(sL,LL)=(c+h,2h)`. Both have `a=c-h>=0`, so (4) is identical. If the probability
of calling early under this common law is `p`, the two errors are `1-p` and
`p`; their maximum is at least `1/2`.

If `ell<2h`, the early class has `B-s+L >= B-c+h`, with equality at
`s=c-h,L=0`. The late class has `B-s+L <= B-c-h+ell`, with equality at
`s=c+h,L=ell`. That last lag is feasible since
`ell<2h<=c+h`. Poisson variables increase stochastically with their mean
(couple a larger mean by adding an independent Poisson increment). Thus the
two threshold errors are respectively largest at these attainable boundary
means, proving (7). Their gap is `mu*(2h-ell)>0`. This proves (8).

For completeness, when both means are positive, the likelihood ratio is
`exp[-E*(lambdaE-lambdaL)]*(lambdaE/lambdaL)^X`, an increasing function of `X`.
The test selecting early when this ratio is at least one is a threshold rule.
The Hellinger affinity is obtained by summing the exponential series:

$$
\sum_{j\ge0}\sqrt{P_E(j)P_L(j)}
=\exp[-E(\sqrt{\lambda_E}-\sqrt{\lambda_L})^2/2].              \tag{11}
$$

On the early error region, `PE <= sqrt(PE*PL)`; on the late error region,
`PL <= sqrt(PE*PL)`. Summing separately bounds each error by (11), proving
(9). For any decision rule, the sum of the errors at these two simple
boundary distributions is at least `1-TV(PE,PL)`. Cauchy–Schwarz applied to
`|sqrt(PE)-sqrt(PL)|*(sqrt(PE)+sqrt(PL))` gives
`TV <= sqrt(1-affinity^2)`. Errors each at most `alpha` imply
`1-2alpha <= sqrt(1-exp(-E*Delta))`; rearrangement gives (10). The point-mass
case `lambdaL=0` follows directly and obeys the same affinity bound.

For unknown positive rates, choose any two positive branch lengths `t1,t2`
and set `mu2=mu1*t1/t2`; their Poisson means agree. Finally, on
`a in [aL,aU]`, the inequality `a<=s<=a+ell` proves (c). **QED.**

There is no uniform finite-error classifier across an arbitrarily small
early/late time gap: for `ell>0` the classes overlap, and even with `ell=0`
Poisson laws at times approaching `c` approach each other. Individual times
separated from `c` can still be consistently classified as exposure increases.

### Theorem 3 — Sharp coverage and depth envelope from two moments

**Assumptions:** O6. In a specified target population, `F` is any random variable
on `[0,1]` with **population** moments `m=E[F]`, `s2=E[F^2]`. Thus
`0<=m<=1` and `m^2<=s2<=m`. Let `R>=0` be an integer and `cD in [0,1]` be the
known per-positive-region marker detection probability. This theorem concerns
the probability averaged over that population, not every individual tumour.

The probability of missing every witness is exactly

$$
A_R(D)=\mathbb E[(1-c_DF)^R].                                \tag{12}
$$

For `0<m<1`, define

$$
b=\frac{s_2}{m},\quad p_0=1-\frac{m^2}{s_2},\quad p_b=\frac{m^2}{s_2},
\qquad
a=\frac{m-s_2}{1-m},\quad
w_a=\frac{(1-m)^2}{1-2m+s_2},\quad
w_1=\frac{s_2-m^2}{1-2m+s_2}.                                \tag{13}
$$

Then the exact attainable interval is

$$
\boxed{\quad
L_R(D)=w_a(1-c_Da)^R+w_1(1-c_D)^R
\ \le A_R(D)\le\
U_R(D)=p_0+p_b(1-c_Db)^R .\quad}                             \tag{14}
$$

Both ends are attained simultaneously for every `R,D`, by distributions
`F in {a,1}` with masses `(wa,w1)` and `F in {0,b}` with masses `(p0,pb)`,
respectively. Every intermediate value is attained. For `R=0,1,2`, the two
ends agree and are, respectively, `1`, `1-cD*m`, and
`1-2*cD*m+cD^2*s2`.

If `m=0`, `F=0` and the answer is `1`; if `m=1`, `F=1` and the answer is
`(1-cD)^R`. If `s2=m^2`, `F=m` almost surely; if `s2=m`, `F` is supported
on `{0,1}`. These are exact degeneracies, not divisions by zero. At `R=0`
the empty experiment has miss probability one, including `cD=1,F=1`.

Suppose the population miss probabilities at a pilot depth with `c0>0` are
`A1` for one region and `A2` for both of two independent regions. Then

$$
m=\frac{1-A_1}{c_0},\qquad
s_2=\frac{A_2-2A_1+1}{c_0^2}.                               \tag{15}
$$

The pilot values are feasible if and only if

$$
1-c_0\le A_1\le1,\qquad
A_1^2\le A_2\le(2-c_0)A_1-(1-c_0).                          \tag{16}
$$

**Proof.** Conditional on `F`, each region fails with probability `1-cD*F`;
conditional independence gives (12). For any `F in [0,1]`, nonnegative variance
gives `m^2<=s2`, and `F^2<=F` gives `s2<=m`. Conversely the upper distribution
in (13) has nonnegative masses summing to one, `b in [0,1]`, mean `m` and
second moment `s2`, so these conditions suffice. The lower distribution has
the same properties, as direct substitution verifies.

Put `f(x)=(1-cD*x)^R`. For `R>=3`,
`f'''(x)=-R(R-1)(R-2)cD^3(1-cD*x)^(R-3)<=0`. The quadratic Hermite
interpolant `QU` agrees with `f` at `0` and at `b`, and with its derivative at
`b`. For `x in [0,1]`, the interpolation remainder has the sign of
`f'''(xi)*x*(x-b)^2`, which is nonpositive. Hence `f<=QU`. The quadratic
`QL` agreeing at `a` with multiplicity two and at `1` has remainder with sign
`f'''(xi)*(x-a)^2*(x-1)`, which is nonnegative, so `QL<=f`. Expectations of
quadratics depend only on `1,m,s2`. Evaluating those expectations under the
respective interpolation-node distributions proves (14) and equality at
the endpoints. Mixing the two distributions preserves the moments and
interpolates their expectations, proving the entire attainable interval.
For `R<=2`, expectations depend on these moments directly. Zero variance and
`E[F(1-F)]=0` give the stated degeneracies. Expanding (12) for `R=1,2` gives
(15); substituting the feasible moment inequalities gives (16). **QED.**

**An exact algebra certificate for every integer `R>=3`.** The verifier does
more than check the sign on a grid. Let `h_j(y1,...,yk)` be the sum of all
monomials of total degree `j` in its arguments, with `h_0=1`. With `c=cD`,
the two interpolation remainders factor as

$$
\begin{aligned}
Q_U(x)-f(x)&=c^3x(x-b)^2
 h_{R-3}(1,1-cb,1-cb,1-cx),\\
f(x)-Q_L(x)&=c^3(1-x)(x-a)^2
 h_{R-3}(1-ca,1-ca,1-c,1-cx).                                \tag{17}
\end{aligned}
$$

Every factor and every argument in (17) is nonnegative on `[0,1]`. To prove
the identity without a differentiability heuristic, the third divided
difference of the monomial `y^R` at four arguments is `h_{R-3}`. For distinct
arguments this follows by expanding the partial fractions of
`product_i (1-y_i*t)^(-1)` and comparing the coefficient of `t^(R-3)`;
alternatively induction on the number of divided-difference arguments gives
the same geometric-series identity. Continuity gives repeated arguments.
The affine substitution `y=1-cx` multiplies a third divided difference by
`(-c)^3`. The interpolation remainder is that divided difference times its
three node factors, yielding (17). Thus the positive-monomial representation
and coefficient equality are an exact certificate, including node coincidences.

### Theorem 4 — Exact design frontier and honest finite-pilot guarantees

**Assumptions:** those of Theorem 3; for its ordinary design formulas take
`0<v<1`, `0<m<=1`, positive integers `R,D`, and target spatial miss
`0<beta<1`. Theorem 2 assumptions are additionally required for timing. The
finite-pilot part assumes independent target-population cases within each
pilot panel, a common known `c0>0`, and a new target case independent of the
entire pilot.

**(a) Attainability and conditional minima.** Write `p0,pb,b` as in (13),
with `m=1` interpreted as `p0=0,pb=b=1`. Some finite `R,D` can guarantee
spatial miss at most `beta` uniformly over the moment class **if and only if**

$$
\beta>p_0=1-\frac{m^2}{s_2}.                                \tag{18}
$$

At a fixed positive depth, the exact smallest region count is

$$
R_{\min}(D)=\left\lceil
\frac{\log[(\beta-p_0)/p_b]}{\log(1-bc_D)}
\right\rceil .                                               \tag{19}
$$

This formula applies when (18) holds. The implementation evaluates rational
powers at the candidate and its predecessor, so logarithmic rounding never
determines a certificate. At a fixed positive region count, first reject
`beta<=p0` as unattainable; otherwise put

$$
U_R(\infty)=p_0+p_b(1-b)^R,\quad
q_*=[(\beta-p_0)/p_b]^{1/R},\quad c_*=(1-q_*)/b.               \tag{20}
$$

Finite depth exists precisely when `beta > U_R(infinity)`, equivalently
`c*<1`; its exact minimum is

$$
D_{\min}(R)=\left\lceil\frac{\log(1-c_*)}{\log(1-v)}
\right\rceil .                                               \tag{21}
$$

At `v=0`, detection is impossible. At `v=1`, depth one already gives `cD=1`,
so equality at a limiting bound may be attainable. When `b=cD=1`, a single
region reaches the floor exactly. The code treats these endpoints directly.
`R=0` or `D=0` means no detection; only a target of one is met. Search-budget
exhaustion is reported separately from mathematical unreachability.

**(b) Joint timing guarantee.** Choose a clock rule whose worst conditional
error on witness capture is at most `alpha_mu`, and a spatial design with
miss at most `beta`. Regardless of decisions made when the witness is missed,

$$
P(\text{wrong timing class})\le\min(1,\beta+\alpha_\mu).       \tag{22}
$$

For maximum **class-conditional** error, the moment constraint must hold
separately within each true early/late class. A constraint only for the pooled
population gives a pooled error bound. It gives no guarantee for each tumour.
Independence of the two error events is not required; O5 is required to make
the conditional clock error bound valid. The sum in (22) is sufficient and is
not claimed to be an exact minimax joint risk.

**(c) Finite-pilot correction and explicit pilot sizes.** In `n1` independent
cases record a Bernoulli indicator for one region hitting the witness, with
average `H1hat`. In `n2` independent cases record whether **both** of two
conditionally independent regions hit, with average `H2hat`. Thus the
population expectations are `c0*m` and `c0^2*s2`. The two panels may overlap;
independence between the two averages is unnecessary. For chosen positive
error radii `eps1,eps2`, define

$$
m_L=\max(0,(\widehat H_1-\epsilon_1)/c_0),\quad
s_U=\min(1,(\widehat H_2+\epsilon_2)/c_0^2),\quad
\delta=\min(1,e^{-2n_1\epsilon_1^2}+e^{-2n_2\epsilon_2^2}).     \tag{23}
$$

With probability at least `1-delta`, `m>=mL` and `s2<=sU`. If `mL>1` or
`mL^2>sU`, the rectangle is inconsistent: **abstain**, rather than project it
onto a preferred model. Otherwise, for every `R,D`, its sharp robust upper
miss bound is

$$
\sup_{m\ge m_L,\,s_2\le s_U,\,m^2\le s_2\le m}
 U_R(D;m,s_2)
= U_R(D;m_L,\min(m_L,s_U)).                                 \tag{24}
$$

This one confidence event holds simultaneously for all budgets. A design may
therefore be chosen after seeing the pilot without paying a multiple-budget
penalty. Consider the protocol that issues a label only when it can certify
a fixed target `beta` using (24), and otherwise abstains. Its unconditional
probability of **issuing a wrong label** on the independent new case is at
most `min(1,delta+beta+alpha_mu)`. This does not give the same risk conditional
on issuance or certification; dividing by the probability of issuance gives
a corresponding conditional upper bound. An abstention is not a successful
classification; coverage/abstention rates must be reported separately. To make the pilot
failure bound at most a preassigned `delta0`, sufficient explicit sizes are

$$
n_i\ge\left\lceil\frac{\log(2/\delta_0)}{2\epsilon_i^2}
\right\rceil,\quad i=1,2.                                  \tag{25}
$$

These are Hoeffding sufficient sizes, not optimal pilot sample minima. Pilot
population shifts or dependent cases invalidate (23).

**Proof.** For ordinary finite `D`, `0<cD<1`, and for `m>0` one has
`0<1-b*cD<1`. Hence `U=p0+pb*(1-b*cD)^R` is strictly greater than `p0`
at every finite budget and converges to `p0` as `R` increases. The extremal
distribution in Theorem 3 makes the bound necessary; any fixed positive `D`
and sufficiently large `R` make it sufficient when `beta>p0`. Taking
logarithms of a number in `(0,1)` reverses the appropriate inequality and gives
(19). At fixed `R`, `cD` increases strictly to one. Thus `U` decreases strictly
to (20), and the target is attainable at finite depth iff it lies strictly
above that limit. Solving `cD>=c*`, or `(1-v)^D<=1-c*`, gives (21).
The stated endpoints follow by substitution. Equation (22) follows from
`{wrong} subset {witness missed} union {wrong and witness captured}`.

For (23), a self-contained Bernoulli concentration argument suffices. If a
random variable lies in `[0,1]`, its exponentially tilted variance is at most
`1/4`. Therefore the second derivative of the centered log moment-generating
function is at most `1/4`; integrating from zero twice gives
`log E[exp(t*(Y-EY))] <= t^2/8`. Independence within a panel gives the sum's
bound `n*t^2/8`. Markov's inequality and minimizing over `t=4*eps` give each
one-sided tail `exp(-2*n*eps^2)`. Apply the upper-tail bound to `H1hat` and
the lower-tail bound to `H2hat`, then use the union bound. This yields (23)
and (25) without independence between panels.

To prove (24), for `m>0` put `b=s2/m` and

$$
H(b)=\frac{1-(1-c_Db)^R}{b}
    =c_D\sum_{j=0}^{R-1}(1-c_Db)^j,\qquad U=1-mH(b).          \tag{26}
$$

On `[0,1]`, `H>=0` and `H'<=0`. Thus, at fixed `m`,
`partial U/partial s2=-H'>=0`; at fixed feasible `s2`,
`partial U/partial m=-H+bH'<=0`. If `sU<mL`, these inequalities bound every
feasible pair by `(mL,sU)`, which is feasible when the rectangle is nonempty.
If `sU>=mL`, convexity gives
`(1-cD*x)^R <= 1-x+x*(1-cD)^R`, so
`U <= 1-m+m*(1-cD)^R <= 1-mL+mL*(1-cD)^R`.
The latter is attained by moments `(mL,mL)`. Both cases give (24); `mL=0`
has bound one. If the confidence event holds, an inconsistent rectangle cannot
occur. On its failure use the trivial error bound one; on its success use
(22). This gives the asserted joint guarantee. **QED.**

### Theorem 5 — Sharp correction for distinct-region downsampling

**Assumptions:** a fixed integer pool of `n>=2` potential regions per case;
`U in {0,...,n}` is the number of witness-negative regions. The experiment
chooses `0<=R<=n` **distinct** regions uniformly without replacement. Every
positive region independently detects its marker with probability `cD`.
Pilot summaries are **population** one- and two-distinct-region miss rates
under perfect marker detection:

$$
A_1^{\mathrm{pool}}=E[U/n],\qquad
A_2^{\mathrm{pool}}=E[U(U-1)/(n(n-1))].                       \tag{27}
$$

These are moments of **invisible** regions; they are not the visible-fraction
moments `m,s2` of Theorem 3. Given `U=u`, the exact miss probability is

$$
g_{n,R,c_D}(u)=\sum_{j=0}^R
 \frac{\binom{n-u}{j}\binom{u}{R-j}}{\binom nR}(1-c_D)^j,      \tag{28}
$$

where out-of-range binomial coefficients are zero. At perfect depth this is
`choose(u,R)/choose(n,R)`. The exact sharp bounds given (27) are the minimum
and maximum of

$$
\sum_{u=0}^n p_u g(u),\quad p_u\ge0,\quad
\sum p_u=1,\quad
\sum u p_u=nA_1^{\mathrm{pool}},\quad
\sum u(u-1)p_u=n(n-1)A_2^{\mathrm{pool}}.                     \tag{29}
$$

If the constraint set is nonempty, both extrema have an attaining
distribution on at most **three** integers. Enumerating every support of
size one, two, or three and solving its three moment equations is therefore
an exact general algorithm. For rational inputs it uses rational arithmetic
only; there is no optimization tolerance. Feasibility itself is checked by
this enumeration. The number of candidate supports is `O(n^3)`; this is a
general finite-pool theorem, though the implementation is intended for
region counts rather than cell-sized pools.

**Proof.** Of the `R` selected regions, `j` are positive with the hypergeometric
probability in (28). All of their read experiments miss with probability
`(1-cD)^j`, giving (28). Averaging over `U` gives the linear objective (29),
and the pilot statistics give its constraints. The feasible set is a closed
subset of the probability simplex, so its extrema exist. If an extremizer
has more than three positive masses, its corresponding vectors
`(1,u,u(u-1))` are linearly dependent. Hence a nonzero signed perturbation
of those masses preserves all three constraints. It has both signs because
its coordinates sum to zero. Sufficiently small perturbations in either
direction remain feasible. At an optimum the objective derivative along the
perturbation must be zero; otherwise one direction improves the objective.
Move until at least one mass vanishes, preserving optimality. Repeating
reduces support to at most three. Distinct integer nodes have independent
three-column moment vectors, and one/two-node candidates are checked against
all three equations. Enumeration therefore finds the global extrema and
attaining distributions. Mixtures attain all intermediate objective values.
**QED.**

The fixed-`n` conversion, if a grid-supported visible fraction `F=1-U/n` is
already known, is

$$
A_1^{\mathrm{pool}}=1-m,\qquad
A_2^{\mathrm{pool}}=
\frac{n(1-2m+s_2)-(1-m)}{n-1}.                              \tag{30}
$$

Pooling cases with different `n` without stratification does not preserve
this conversion. A one-region published summary alone does not identify the
two moments in (29). At `n=2,u=1,R=2,cD=1`, (28) is zero, while the
with-replacement probability is `1/4`: a direct counterexample to exchanging
the designs.

### Theorem 6 — Founder-genotype recovery and what remains uncountable

**Assumptions for recovery:** `K` distinct, mutually exclusive, surviving
founder lineages have diagnostic labels known to have existed **before**
seeding; they are retained without recurrence/loss and are distinguishable
from de novo post-seeding subclones. A sample is an independent categorical
cell draw with probabilities `p1,...,pK>=0`, `sum pi<=1` (remaining probability
is background). Its correct clone label is detected with probability
`cD`, independently of identity and other draws. This is a cell-sampling
subexperiment; it does not assert that bulk-region observations are independent
cell labels.

With `M>=0` draws, the exact probability of finding every label is

$$
P(K_{\mathrm{obs}}=K)=
\sum_{A\subseteq\{1,\ldots,K\}}(-1)^{|A|}
\left(1-c_D\sum_{i\in A}p_i\right)^M.                        \tag{31}
$$

If `pi>=pmin>0` for each label, `K<=Kmax` with `Kmax>=1`, `cD>0`, and
`0<beta<1`, then

$$
P(K_{\mathrm{obs}}<K)\le K_{\max}(1-c_Dp_{\min})^M,\qquad
M\ge\left\lceil
\frac{\log(\beta/K_{\max})}{\log(1-c_Dp_{\min})}
\right\rceil                                                \tag{32}
$$

suffices. The endpoint `cD*pmin=1` needs one draw. An empty label set needs
zero draws. The union bound is a
sufficient global sample size, not a claim of exact multi-clone optimality.
Equation (31) supplies exact minima when the full frequency vector is known.

For a spatial lesion-by-founder incidence matrix with at most `Jmax` eligible
nonzero entries, apply Theorem 3 separately to each entry's diagnostic
visibility class. If each has upper miss bound `U_R(D)`, the probability of
missing any entry is at most `Jmax*U_R(D)`. No independence between entries
is needed. For `Jmax>=1`, spatial prevalence at least `fmin>0`, positive-region
VAF at least `vmin>0`, and fixed `D>=1`, a sufficient region count is

$$
R\ge\left\lceil
\frac{\log(\beta/J_{\max})}
{\log[1-f_{\min}(1-(1-v_{\min})^D)]}
\right\rceil.                                               \tag{33}
$$

If the single-region detection probability is one, one region suffices; if it
is zero, a nonempty set cannot reach a target below one. An empty incidence
set needs zero regions. The per-entry miss probability and corresponding
per-entry minimum are sharp;
the simultaneous union bound need not be. With the stated genotype provenance
and no false positives, a fully recovered incidence matrix gives the number
of distinguishable founder lineages per sampled lesion and the number shared
across those lesions.

Without a positive visibility bound, no finite sampling budget uniformly
recovers all extant genotypes. For the non-identifiability statements, allow
clusters and repeated transfers as unrestricted historical nuisance variables,
while fixing the endpoint genotype probabilities. Even perfect recovery of
the matrix does not identify the number of such migration events or genetically
identical seeding cells. Extinct founding lineages (when extinction is allowed)
and unsampled lesions are also unobserved. If a model assumes exactly one
founding event per lesion, its event count is fixed by that assumption, not
inferred from genotype support.

**Proof.** Let `Ei` mean that label `i` is missed. A single draw misses every
label in a subset `A` with probability `1-cD*sum_A pi`. Independence over
draws gives the power in (31). Inclusion-exclusion for the complement of
`union_i Ei` proves (31). Each `P(Ei)=(1-cD*pi)^M`; the union bound and a
prevalence lower bound prove (32). The identical union argument applied to
region detections proves (33) and the incidence-matrix statement.

For the visibility impossibility, move `epsilon` probability mass from
background to an extra diagnostic founder label,
keeping total probability fixed. Couple the observations under its absence and its
presence so that they differ only when the rare label is detected. Their
total variation is at most `M*cD*epsilon`, which tends to zero for every
fixed `M`. Thus a count classifier cannot have uniform maximum error below
`1/2` over the absence/vanishing-frequency alternatives; the probability of
missing the added label tends to one. A label with zero visibility is exactly
unobservable. The analogous spatial construction takes `F` arbitrarily close
to zero, or exactly zero.

Finally, the observation map that records diagnostic genotype presence
forgets the number of transfers of each genotype. One cluster containing
two diagnostic genotypes and two separate transfers of those genotypes
produce the same support. Repeated transfers of an already represented
genotype do not change it. When extinction is allowed (for example `dM>0`),
arbitrarily many extinct transfers do not add a recorded label; neither do
unobserved lesions. Conditional on the same surviving
sampled genealogy, extinct branches contribute no passenger pattern at any
exposure in (3). Consequently these histories cannot be separated by the
stated endpoint experiment. A fully specified stochastic model with sizes,
times, or longitudinal observations may assign different likelihoods to
such histories; that is additional information, not identification by
genotype support. **QED.**

## 4. Numerical consequences, predictions, and external falsification

### 4.1 Reproducible design examples (synthetic parameters)

These are mathematical examples in arbitrary time/exposure units. They are
not parameter estimates or recommended biopsy/sequencing designs.

**Clock example.** Set `B=10`, `c=15/2`, `h=3/2`, `ell=1`, and `mu=1`.
The least-favourable rates are `lambdaE=4`, `lambdaL=2`. At target error
`alpha_mu=1/20`, the exact minimum integer exposure for deterministic
thresholds is **9**, with **early iff `X>=26`**. Its errors are approximately
`0.03453279` (early class) and `0.04460830` (late class). At exposure 8,
the smallest threshold satisfying the late-error target is 24, but the early
error is approximately `0.06096940`; all smaller exposures also fail. The
verifier encloses these tails by rational intervals and certifies every
rejected exposure. More reads of the same already-called sites are not
additional independent clock exposure.

**Coverage example.** Let `m=1/2`, `s2=3/10`, `v=1/2`, and target
`beta=1/5`. At depth four, `cD=15/16`. Four regions have the exact range

$$
[L_4,U_4]=[16667/131072,\ 25847/131072]
\simeq[0.1271591,\ 0.1971970].                              \tag{34}
$$

At this fixed depth three regions give upper miss `1937/8192 > 1/5`, so four
is minimal. At four regions depth three gives `214107/1024000 > 1/5`, so
depth four is minimal. These are **conditional** minima: no claim is made
that the pair minimizes all possible laboratory cost functions.

The same first two moments permit `F=0` with probability `1/6` and `F=3/5`
with probability `5/6`. Consequently a **5%** worst-case miss target is
unattainable from this calibration at any finite region/depth budget. They
also permit `F=2/5` with probability `5/6` and `F=1` with probability `1/6`,
for which misses do vanish as regions increase. The two-moment floor is
therefore uncertainty about the latent distribution, not a proof that the
actual miss rate must remain above `1/6`.

Combining the first two examples by (22) gives a timing error upper bound
below `0.2472` with known moments. It does **not** give a 5% joint timing
guarantee. Using a finite pilot adds its failure allowance as in Theorem 4.
At the moment level the condition for finite-budget spatial error `beta` is
equivalently `s2 < m^2/(1-beta)`. In terms of the coefficient of variation of
positive-mean `F`, it is `CV(F)^2 < beta/(1-beta)`. This is a constraint on
how much heterogeneity the calibration may leave unresolved.

**Distinct-region example.** With `n=8`, invisible-region pilot moments
`A1pool=1/2`, `A2pool=1/4`, `R=4`, and `cD=3/4`, the exact finite-pool range is

$$
[979/7168,\ 1249/7168].                                      \tag{35}
$$

The lower endpoint is attained at `U=0,4,5` with masses `1/10,1/2,2/5`;
the upper at `U=3,4,8` with masses `2/5,1/2,1/10`. The independent simulator
also uses a different `n=6` fixture, so the two implementations are not merely
replaying a single example.

### 4.2 Published aggregate benchmark and an explicitly failed prediction

**Source fact:** [Al Bakir et al. 2023, Figure 2d and accompanying
text](https://doi.org/10.1038/s41586-023-05729-x) report that a single randomly
downsampled primary region changed 75 of 90 late-divergence classifications
to early. Their main series has 94/126 late-divergence cases; the 90-case
downsampling denominator is a different eligible subset and is retained here.
The endpoint is **late divergence being called early**, not physical migration
before first clinical detectability.

**Prediction P1 — optional broad-witness submodel, rejected conditionally.**
Suppose each of these eligible cases has a correctly detectable diagnostic
witness in at least half of its potential primary regions (`F*cD>=1/2`),
false early-divergence calls correspond exactly to witness misses, and
casewise downsampling trials are independent. This predicts an expected
number of misses at most 45 (a fraction at most 45/90) and the following
bound for 75 or more misses:

$$
P(S\ge75)\le\sum_{j=75}^{90}\binom{90}{j}2^{-90}
=\frac{1774990439854189}{38685626227668133590597632}
\simeq4.58824275\times10^{-11}.                              \tag{36}
$$

To prove the bound, couple each independent casewise miss Bernoulli to a
Bernoulli with probability `1/2` using the same independent uniform random
number. The former is never larger; summing gives stochastic domination by
`Bin(90,1/2)`. The rational sum in (36) is evaluated exactly by
`published_aggregate_check()`.

The reported aggregate therefore fails this **optional broad-coverage
submodel** under those additional assumptions. This is a transparent,
post-hoc falsification attempt, not the original study's p-value, a new
empirical discovery, or validation of the lag/clock model. If classification
failure is not a simple witness miss, or the trials are dependent, the tail
calculation is inapplicable. The observed count remains evidence against
casually assuming broad informative spatial coverage. We have not established
which mechanism caused each reported misclassification.

### 4.3 Predictions requiring a matched aggregate sampling design

**Prediction P2 — heterogeneity constrains any resampling extrapolation.**
For fresh independent regions or resampling with replacement, let `Q` be an
individual's one-region miss probability. Then for all integer `R>=1`,

$$
(E Q)^R\le E[Q^R]\le E Q.                                   \tag{37}
$$

The lower bound is Jensen's inequality for `x^R`; the upper is `x^R<=x` on
`[0,1]`. They are sharp, attained by constant `Q` and by `Q in {0,1}`,
respectively. Treating `75/90=5/6` as an exact population value **only for
illustration**, `(5/6)^16≈0.05408789` and `(5/6)^17≈0.04507324`. Thus 17
regions are a **necessary lower floor**, not a uniformly sufficient count
for 5% miss. The compatible mixture `P(Q=1)=5/6` has miss `5/6` forever.
The finite 90-case estimate has uncertainty, so even that floor is not a
confidence-certified design for the actual population.

**Prediction P3 — a sharper third-region interval.** At perfect marker
detection, population one- and two-iid-region miss summaries `A1,A2` imply

$$
\frac{A_2^2}{A_1}\le A_3\le
A_2-\frac{(A_1-A_2)^2}{1-A_1},\quad 0<A_1<1.                \tag{38}
$$

This follows by applying Theorem 3 to `Q=1-F`, or substituting (15) into
(14). For a concrete preregisterable calibration `A1=1/2,A2=3/10`, it
predicts **18–22% three-region misses**, with both endpoints attainable.
The entire curve must obey (14) across depths if the common-VAF assumption
holds. An independently published matched `A3` outside its uncertainty-adjusted
interval would falsify the observation assumptions (or their calibration).
For distinct-region published downsampling, use Theorem 5 at the actual pool
size instead of (37)–(38).

**External status:** P1 has the reported aggregate comparison above. P2–P3
are falsifiable quantitative predictions for matched aggregate downsampling,
but are **not externally tested here**: the cited one-region summary supplies
neither the matched second moment nor pool-size strata. No biological success
is claimed from fitting the same statistic that calibrated a prediction. A
future aggregate release should report its region selection rule, pool sizes,
one/two/three-region miss counts and denominators, VAF/read calling rule, and
casewise dependence. This requests aggregate reporting, not patient records.

### 4.4 Direct counterevidence against universalizing the generator

[Hessey et al. 2026](https://doi.org/10.1038/s41586-026-10428-4) report
metastasis-to-metastasis sources for 156/258 lesions and mixed primary/metastatic
sources for another 4/258; 187/258 lesions had one inferred founding subclone.
Those observations reject applying G3 as a universal description of that
cohort. They do not falsify a theorem explicitly conditional on G3 or turn
founding subclone number into a count of founding cells. This failure is
recorded rather than suppressed by selecting the 98/258 primary-only lesions.

## 5. Falsification and hostile-referee record

The development deliberately tested attractive but invalid shortcuts. These
are counterexamples or failed model extensions, not failed biological
experiments. Review was by independent AI agents plus independent code
implementations; it is not external scientific peer review.

| Attempted claim or attack | Result and correction |
| --- | --- |
| Ordinary metastasis-private counts directly date migration. | Rejected: moving the location mark on a sampled ancestral edge preserves every passenger law. A calibrated lag or additional information is necessary. |
| Use the biopsy date as the endpoint of private **clonal** bulk counts. | Rejected: the relevant endpoint is the sampled metastatic MRCA. The single-cell/path experiment is stated separately. |
| Replace the birth-death primary by its exponential mean. | Rejected as an exact marginal law: the Cox mixture adds `Var(Lambda_N)` and changes zero-count probabilities. The simulation retains a genuine event-generated primary. |
| A gap `2h=ell` is enough for separation. | Rejected at equality: `(s,L)=(c-h,0)` and `(c+h,2h)` give the same anchor. The correct condition is strict. |
| Timing remains identified with an unconstrained clock rate. | Rejected by rate–duration rescaling giving exactly equal Poisson means. |
| Positive mean spatial visibility ensures arbitrarily small robust error. | Rejected: the two-moment extremizer places `1-m^2/s2` mass at zero visibility. |
| A heterogeneity floor is inevitable for every compatible population. | Rejected: the lower-extremizer population can have strictly positive visibility everywhere with the same two moments. The floor describes worst-case uncertainty. |
| The 17-region homogeneous illustration is sufficient for TRACERx. | Rejected: Jensen gives the opposite direction; a compatible heterogeneous population never improves. The downsampling design also matters. |
| Distinct biopsies can be treated as draws with replacement. | Rejected by `n=2,u=1,R=2`: zero miss without replacement versus `1/4` with replacement. Theorem 5 supplies the correction. |
| Replace population moments by empirical moments in an exact guarantee. | Rejected: finite-pilot errors are controlled separately, and inconsistent rectangles cause abstention. |
| The finite-pilot bound is valid conditional on receiving a certificate. | Rejected: selecting certified cases can concentrate pilot failures. The final bound is unconditional for issuing a wrong label on an independent new case. |
| Per-entry union-bound budgets are exact optimal global budgets. | Rejected: the union bound is only sufficient. Exact categorical recovery uses inclusion-exclusion; exact threshold exposure is restricted to the stated deterministic rule family. |
| Any label appearing only in a metastasis is a founder clone. | Rejected: de novo mutations and clone loss require pre-seeding provenance assumptions. |
| Support recovery counts cell clusters or all migration events. | Rejected by repeated identical-genotype transfers and extinct unsampled descendants. |
| Linked/adjacent regions obey the iid moment law automatically. | Rejected: the simulation includes perfectly linked regions, whose miss rate does not follow the iid formula. Spatial design is an assumption. |
| A primary-only migration model describes all published cohorts. | Rejected by the 2026 source-count aggregates in Section 4.4. |
| A correct formula at a few rational inputs proves its general theorem. | Rejected: complete proofs are given above; numerical checks are independent implementation diagnostics and algebra certificates. |
| Reversed Poisson rates mean the models are confounded. | Rejected during API review: reversed unequal rates require the opposite decision direction. The API now rejects this unsupported orientation and reserves `confounded` for equal rates. |

No mathematical claim in the numbered theorems is being left as a
simulation-only conjecture. Two broader conjectures remain explicitly open:

* **Conjecture C1 (biological applicability, unproved):** a suitable untreated
  cohort and prespecified marker protocol can yield a useful independently
  calibrated `ell` smaller than a scientifically relevant separation `2h`,
  together with coverage moments tight enough for a small error target.
  Published improvements with multiregion sampling motivate trying this, but
  neither the cited aggregate counts nor our simulation establish such a lag
  bound. It is not a conclusion of the birth-death model.
* **Conjecture C2 (priority, unestablished):** the exact combined coverage/depth,
  bounded-lag, finite-pilot and finite-pool formulation is not already present
  in the metastasis literature. Targeted searches found no direct match;
  incomplete recall, terminology differences, preprints, and broader moment
  literature can invalidate this conjecture. The separate ingredients are
  acknowledged prior mathematics. Independent specialist review is needed
  before claiming a novel research result.

## 6. Reproduction, exact arithmetic, and test independence

From the repository root, using Python 3.10 or newer and only its standard
library for this contribution:

```bash
python -B research/metastatic-seeding-timing/verifier.py
python -B research/metastatic-seeding-timing/simulation.py --seed 20261117 --replicates 4000
python -B -m unittest discover -s tests -p 'test_research_metastatic_seeding_timing.py' -v
```

The scripts print reports to standard output and require no downloads, patient
data, external credentials, NumPy, or package installation. A maximum search
budget can give `search_exhausted`; that is not a proof of unreachability.
`unreachable` is returned only when a proven limiting obstruction applies.

### Files and roles

| File | Role |
| --- | --- |
| `README.md` | Assumptions, six theorems, complete proofs, literature audit, predictions, limitations, and review record |
| `verifier.py` | Exact rational moment/occupancy/LP calculations, nonnegative polynomial certificates, rigorously enclosed Poisson tails, certified budget searches, aggregate model check |
| `simulation.py` | Separate fixed-seed event/read simulation; does not import the verifier |
| `../../tests/test_research_metastatic_seeding_timing.py` | Independent stdlib oracles, boundary checks, and all stochastic comparisons; discovered automatically by existing CI |

The verifier rejects binary floating-point inputs. Rational values use
`fractions.Fraction`; rational powers and polynomial coefficients are exact.
Transcendental values are **enclosed**, not falsely described as rational.
For rational `x>=0`, let `Sn=sum_{j=0}^n x^j/j!`. Once `n+2>x`, the positive
exponential remainder is bounded by

$$
0\le e^x-S_n\le
\frac{x^{n+1}/(n+1)!}{1-x/(n+2)}=:R_n.                       \tag{39}
$$

Indeed each successive tail ratio is at most `x/(n+2)<1`. Therefore
`1/(Sn+Rn) <= exp(-x) <= 1/Sn`; multiplying by a rational finite Poisson
sum gives rigorous CDF bounds. Subtraction reverses endpoints for upper
tails. Comparisons refine rational tolerances and **raise** if unresolved;
they never choose a side from rounded decimals. Clipping probability bounds
to `[0,1]` is valid. Birth-death survival bounds apply monotone rational
transformations to these exponential enclosures; the critical case is rational.

At each integer clock exposure the verifier finds the smallest threshold
meeting the late-error constraint. Every larger threshold raises early error,
so if that threshold fails, every threshold fails at that exposure. It checks
every smaller exposure before claiming the minimum. Region and depth searches
use exact monotone powers and certify both the selected integer and its
predecessor. The finite-pool algorithm enumerates all supports of size at most
three using exact moment equations. Its exhaustive guarantee is proved in
Theorem 5, not inferred from solver agreement.

The independent test uses high-precision `decimal.Decimal.exp` and direct
Poisson sums to check the rational enclosures; Fraction Gaussian elimination
to reconstruct interpolation and finite-pool moment solutions; literal
enumeration of distinct subsets for hypergeometric probabilities; and a
categorical-state dynamic program to check occupancy rather than reuse
inclusion-exclusion. It checks exposure minima by a separate exhaustive search.
The aggregate null tail is checked by independent integer combinatorial sums.

### Stochastic checks and limits

The simulator generates exponential waiting times for mutations and seeding,
Gillespie events for birth and death, Bernoulli region occupancy and literal
read draws, uniform distinct-region samples, and categorical clone labels.
No fitted parameter is chosen to force agreement. The default fixed seed is
`20261117`; alternative seeds are available. Its source has no import of the
verifier or third-party library.

The default run contains **23** stochastic records across Theorems 1–6. The
primary trace starts with 24 cells and reaches 1,937 through 3,205 birth/death
events; no trajectory is capped or discarded. The lineage-seeding experiment
then conditions on that genuine primary trace, exactly as (2) does. Separate
metastatic simulations retain extinction. Conditional survivor-intensity
integration is a floating-point numerical comparison, with an explicit
refinement discrepancy in the report; it is not an exact-integral certificate.
The code also checks shared/private means, shared-count covariance, confounded
anchor laws, both moment extremizers, an interior law, the finite-pilot
confidence event, finite-pool extrema, and clone recovery.

The tests use tolerances declared from analytic variances, not observed
deviations. A passing fixed-seed comparison is a smoke test for the
implementation, not a probability theorem or empirical cancer evidence.
Non-identifiability is established by equal-law constructions, not by inability
of a simulation to find a difference. Boundary algebra is checked exactly.

Local Python 3.12 runs take roughly two seconds for the new unittest module,
well below the requested 60 seconds; timings are environment-dependent.
The existing repository validation workflow discovers the new module without
any workflow change. Broader verification commands are:

```bash
python -B scripts/validate.py
python -B -m unittest discover -s tests -v
node --test tests/*.test.*
```

Only this directory and the one new test file form the contribution. Existing
research notes, repository indexes, catalogues, policies, requirements,
workflows, and publication registries are outside its change set.

## 7. Limitations that change the conclusions

The largest unresolved step is not numerical precision. It is obtaining a
credible calendar clock, terminal ancestral time, detection reference, and
lag bound from real multi-region data. Without them the exact guarantees
cannot be applied to physical early/late migration. Conditioning on an
observed genealogy intentionally discards size/history likelihood information;
the resulting impossibility theorem is not a claim about every possible
fully specified birth-death inference model.

Primary spatial structure, ongoing selection, variable clock rates, dormancy,
therapy, changes in ploidy/purity, false read calls, variant loss, phylogeny
errors, and metastasis-to-metastasis seeding can change the laws. Population
coverage moments must refer to the actual witness definition and target
population; average coverage cannot certify each patient. Distinct-region
downsampling requires its own sampling law and pool sizes. Founder genotype
provenance cannot be established by calling every newly observed variant a
founder. The finite-pool and coverage calculations do not estimate metastatic
growth or independently reconstruct a migration graph.

No external patient-level validation, calibrated human seeding-time estimate,
clinical performance claim, or established research priority is supplied.
The deliverable is a fully specified conditional mathematical analysis with
reproducible computational checks and an honest record of what failed.
