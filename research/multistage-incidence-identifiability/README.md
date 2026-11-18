# What age-specific incidence can identify about multistage carcinogenesis

**Nonclinical computational research. No medical advice. Unreviewed mathematical
and computational work; no biological or clinical validation is claimed.**

This contribution separates three questions that are often conflated:

1. Does an exactly known, correctly observed hazard determine the number of
   transitions in a specified stochastic model?
2. Does it determine the biological rates or the distinction between mutation
   and diagnosis delays?
3. Can finitely many noisy age bins distinguish the transition count robustly?

The answers are different. With constant positive finite mutation rates, no
premalignant cells at birth, constant Poisson initiation, and immediate
observation, the exact hazard's order of vanishing determines the transition
count **even with expansion at every intermediate stage**. Expansion alone
does not produce exact equality between different counts in that family.
Unknown diagnosis delays can produce exact equality. Arbitrarily fast extra
transitions obstruct uniform estimation from finite exposure. Even arbitrarily
small, positive allowance for age-dependent calibration permits exact
finite-bin equality after adding a sufficiently fast transition.

The potentially new contribution is the explicit **finite-bin resolution and
calibration boundary** in Theorems 5–7, together with an observation-aware
classification in Theorem 4 and reproducible certificates. The building blocks
are classical. Priority for this particular synthesis is **not established**.
This is not a claim to have solved cancer's biological driver count. If a prior
source contains these same results, the contribution should be downgraded to
a reproducible synthesis, without changing the mathematical claims.

## 1. Literature first: the necessary pivot

Search completed **2026-11-17, America/New_York**. Primary papers, their
publisher pages, PubMed records, and available full text were searched using
`multistage clonal expansion`, `identifiability`, `number of stages`, `finite
bins`, `minimax`, `Le Cam`, `observation`, `detection`, `phase type`, and
combinations thereof, including searches for newer work. This was a targeted
research search, not a registered systematic review or a guarantee that every
preprint, mathematical equivalent, or unpublished result was found.

| Prior work | What is already known; implication here |
|---|---|
| Nordling (1953) [1]; Armitage & Doll (1954, 1957) [2,3] | Sequential events and age-incidence power laws, followed by two-stage alternatives. The onset power law is not a new discovery. |
| Moolgavkar–Venzon (1979), Moolgavkar–Knudson (1981) [4,5] | Stochastic two-stage mutation/selection models. The distinction between a stochastic first-event hazard and an expected cell count matters. |
| Little (1995) [6] | General multistage extensions of MVK and effects of parameter changes. Arbitrary numbers of stages are not new. |
| Heidenreich (1996); Heidenreich, Luebeck & Moolgavkar (1997) [7,8] | Exact TSCE identifiable combinations and hazard properties. Theorem 2 below is a transparent rederivation of an established result. |
| Luebeck & Moolgavkar (2002); Meza, Jeon, Moolgavkar & Luebeck (2008) [9,10] | MSCE applications to colorectal and pancreatic incidence, including age regions that are not one power law. Adult log-log slopes cannot simply be read as driver counts. |
| Cox & Huber (2007) [11]; Frank (2007) [12] | Expected-count symmetries, observationally equivalent mechanisms, and the effects of expansion and unequal rates. The spectral permutation symmetry below has this important precedent. |
| Little, Heidenreich & Li (2009) [13] | Parameter redundancy and bounds in broad stochastic carcinogenesis models. Nonidentifiability is established beyond TSCE. |
| **Brouwer, Meza & Eisenberg (online 2016; issue 2017)** [14] | Differential-algebra results for arbitrary pre-initiation count with one final expansion, and multiple expansions through five mutation stages, with a conjecture for arbitrary stage number. This is the closest structural predecessor. |
| **Brouwer, Meza & Eisenberg (2017)** [15] | Pancreatic two-, three-, and four-stage models have only three practically identifiable combinations in their analysis. Rediscovering practical ambiguity would not be novel. |
| Richardson, Anghel & Deng (2021) [16] | Incidence-derived multistage/senescence parameters have already been compared with sequence-derived driver estimates. The comparison itself is not new. |
| Mirzaei et al. (2025) [17] | Tumor size and detection enter MSCE-T. Their general model remains unidentifiable; the identified specification uses observed tumor size and a fixed malignant proliferation rate. Incidence alone is not sufficient. |
| Houlahan et al. (2026) [18]; Bandlamudi et al. (2026) [19] | New incidence-plus-genomic event inference and large driver catalogs add biological information; neither establishes that a single aggregate age curve identifies biological drivers. |

**Explicit pivot.** Fixed-model MSCE nonidentifiability is already substantially
solved. We instead ask: *which information about model order survives a
specified observation operator, and what quantitative obstruction remains with
finite age bins and finite exposure?* We do not claim a new complete
parameter classification for arbitrary finite-rate nonlinear MSCE. We give a
complete finite-rate TSCE classification and a complete arbitrary-stage
classification in a **proved limit**, then prove observation and statistical
results for the broader exact branching family. This scope restriction avoids
silently treating a mean count as an exact stochastic hazard.

## 2. Models, observations, and meaning of a “step”

### 2.1 Exact branching family B

The following assumptions apply to Theorems 1 and 2, and to the biological
processes used in Theorems 5–8 unless a theorem explicitly weakens them.

**B1.** Age starts at known time zero. There are no initiated, premalignant, or
malignant cells initially. Normal cells form a deterministic, nondepleting
reservoir. First mutations create independent stage-1 clones according to a
homogeneous Poisson process of rate `lambda > 0` per person per year.
If `lambda = N mu_0`, neither `N` nor `mu_0` is observed separately.

**B2.** There are `k >= 2` mutation transitions: the initial creation of stage
1, then transitions through stages `1,...,m`, where `m=k-1`, ending in the first
malignant event. An optional `k=1` boundary has direct malignant Poisson
immigration with hazard `lambda` and no intermediate cells.

**B3.** At intermediate stage `i`, a cell gives birth to one additional
same-stage cell at rate `alpha_i >= 0`, dies/is lost at `beta_i >= 0`, and
mutates forward at `mu_i > 0`. All rates are constant and finite. Expansion is
optional at any number of stages. All stages may have positive birth rates.

**B4.** A specified flag `c_i` is either zero or one. For `c_i=0`, a forward
mutation creates a next-stage cell and retains its parent (asymmetric
mutation). For `c_i=1`, it replaces its parent (conversion). At the final
mutation the flag does not affect the first malignant time. These are
different models; a fast asymmetric stage is **not** an exponential delay.

**B5.** Cells and immigrant clones evolve independently; there is no feedback,
competition, mutation-rate variation, frailty, inherited shortcut, or reverse
transition. A cell has finitely many event types at finite per-cell rates, so
the process is nonexplosive on every bounded interval. Clone extinction is
allowed. The first malignant event across all clones defines biological onset.

**B6.** The reference hazard is the cause-specific first-onset hazard for the
homogeneous model, conditional on remaining free of that onset. Competing
death, when relevant, must be independent/noninformative for this hazard
interpretation. The reference observation has no unknown diagnosis delay,
screening effect, birth-cohort effect, or age multiplier.

Here `k` counts **modeled mutation transitions**, not an independently defined
number of rate-limiting bottlenecks or a count of selected sequencing variants.
The fast-stage construction below exposes precisely why those concepts differ.
Calling all modeled transitions “rate limiting” without a time-scale
restriction would make the biological question ill defined.

Let `q_i(t)` be the probability that one stage-`i` cell has produced a first
malignant event by elapsed time `t`; set `q_k(t)=1`. Then

\[
q_i'=\mu_iq_{i+1}+(\alpha_i-\beta_i-c_i\mu_i)q_i
       -\alpha_iq_i^2-(1-c_i)\mu_iq_iq_{i+1},\qquad q_i(0)=0.
\tag{1}
\]

At the last stage this is
`q_m' = mu_m + (alpha_m-beta_m-mu_m) q_m - alpha_m q_m^2`.
These are autonomous backward equations in **elapsed clone time**. We do not
insert calendar-age-dependent rates into them; doing so naively is invalid
in general [25].

### 2.2 Observation models kept distinct

**O0: exact direct observation.** The complete continuous hazard is known
without noise on any nonempty open interval of positive ages, with B1–B6.
Analytic continuation is a uniqueness statement, not a stable numerical
estimator of derivatives at birth.

**OD: clone-level diagnosis lag.** After each immigrant clone's *first*
malignant event, one independent lag `D` is attached to that clone. The first
diagnosed clone is observed. Subsequent malignant events within that clone do
not compete to diagnose earlier. This deliberately explicit observation
model is not a lag applied to the person's already-first cancer, and not
independent lags for all malignant cells. When specified, `D` is a sum of
`d >= 0` independent exponentials with strictly positive finite rates
`ell_1,...,ell_d`, independent of all clone biology.

**OP: finite-bin Poisson experiment.** On disjoint age bins `B_j` in `[0,T]`,
known nonnegative exposure densities `w_j(t)` yield independent counts

\[
Y_j\sim\operatorname{Pois}(m_j),\qquad
m_j=\int_{B_j}w_j(t)h_{\rm obs}(t)\,dt.                 \tag{2}
\]

This is an explicitly assumed aggregate counting experiment. Independent
Poisson counts are not the exact joint law of a fixed birth cohort followed
until one diagnosis; such a cohort instead has survival depletion and a
multinomial counting structure. A known at-risk exposure design or a stated
Poisson approximation is needed. Our lower bounds concern (2), not an
unstated individual-level likelihood.

**OC: age-bin calibration.** For a positive shape vector `f_j`, a free common
amplitude `A>0`, and `rho>=1`, reported mean rates satisfy
`y_j=A a_j f_j`, with `a_j in [1/rho,rho]`. The multipliers represent an
explicit sensitivity allowance for age-shape distortion. They are not
estimated registry completeness, and the allowance is not inferred from the
incidence table. Free overall scale can equivalently normalize relative
thinning probabilities, but is not itself biological evidence of thinning.

The SEER CSV contains **rounded rates only**, not the counts or exposures
needed to instantiate OP. Its analysis uses OC and shape checks, with a
declared uniform within-bin weighting assumption and a separate sensitivity
bound for unknown within-bin weights.

## 3. Numbered results and complete proofs

### Theorem 1 — Exact branching hazard and identifiable onset order

**Assumptions:** B1–B6, O0, `k` a positive finite integer, positive finite
initiation/mutation rates and finite nonnegative birth/death rates. Flags may
differ by stage and all intermediate stages may expand.

For `k>=2`,

\[
S(t)=\exp\{-\lambda\int_0^t q_1(u)\,du\},\quad
h(t)=\lambda q_1(t),\quad
h(t)=\frac{\lambda\prod_{i=1}^{k-1}\mu_i}{(k-1)!}t^{k-1}
       +O(t^k).                                             \tag{3}
\]

Thus `k = 1 + ord_0(h)`. For `k=1`, `h(0)=lambda>0` and the same convention
gives `k=1`. Equality of hazards on any open age interval implies equal `k`.
The product `lambda prod mu_i` is identifiable from the onset derivative;
this statement alone does not classify all biological parameters.

**Proof.** Condition on an event of one stage-`i` cell in an infinitesimal
interval. A birth creates two independent future lineages, with failure
probability `(1-q_i)^2`; death leaves no lineage; asymmetric mutation leaves
two types with failure probability `(1-q_i)(1-q_{i+1})`; conversion leaves
one next type. Expanding and subtracting the no-event probability gives (1).
The final malignant transition has failure probability zero, which is why
`q_k=1`. Its convention cancels from the last equation.

An immigrant at time `s` produces a malignant event by age `t` with probability
`q_1(t-s)`. Independent marking of Poisson arrivals makes the number of
successful immigrant clones by `t` Poisson with mean
`lambda integral_0^t q_1(u)du`. Its zero probability is `S(t)` in (3), and
logarithmic differentiation gives the hazard. This is exact; expected
malignant-cell counts have not been substituted for a hazard.

The polynomial vector field (1) has a locally analytic solution. Probability
bounds `0<=q_i<=1` and finite rates continue it over every finite positive
time interval. At the terminal stage, `q_m(t)=mu_m t+O(t^2)`. Induct downward:
if `q_{i+1}=C t^r+O(t^{r+1})`, write the exact equation as
`q_i'=mu_i q_{i+1}+a_i(t)q_i`, where
`a_i(t)=alpha_i-beta_i-c_i mu_i-alpha_i q_i-(1-c_i)mu_i q_{i+1}`
is bounded and analytic near zero. The integrating factor is `1+O(t)`.
Variation of constants, starting at zero, gives
`q_i=mu_i C t^{r+1}/(r+1)+O(t^{r+2})`. Thus the nonlinear terms cannot
change this coefficient, without assuming the desired order in advance. Equivalently,
recursion of ordinary Taylor coefficients gives the same induction without
assuming a nonzero lower coefficient. Repeating proves (3). The coefficient
is strictly positive, so it cannot cancel. Real-analytic identity on a
connected interval extends equality from an open positive-age interval to
all positive ages and then to the germ at zero. Equal germs have equal
vanishing order. This proves the claim. **QED.**

**Boundaries.** Zero mutation rates give the identically zero hazard and lose
`k`; infinite rates are outside this theorem. Unknown premalignant initial
conditions, a hidden age clock, parallel paths, or unknown observation
filters invalidate its assumptions. A mixture of positive clean-birth paths
has the *smallest* path order at onset, not every path's length. Neither an
adult-age fitted slope nor a rounded table supplies the germ in (3).

### Theorem 2 — Complete finite-rate TSCE equivalence class (established)

**Assumptions:** B1–B6, O0, `k=2`, `lambda,mu>0`, `alpha,beta>=0`, all finite.

Exactly the following three combinations determine the complete hazard:

\[
a=\lambda\mu>0,\qquad b=\alpha-\beta-\mu,\qquad c=\alpha/\lambda\ge0.
\tag{4}
\]

Two feasible TSCE parameter sets have identical hazards **if and only if**
their triples `(a,b,c)` agree. Every identifiable function of the parameters
must therefore be constant on these equivalence classes. This is equivalent
to the established triples involving `lambda mu`, `alpha mu`, and
`alpha-beta-mu` [7,8,14], since `alpha mu=ac`.

The class for a feasible triple consists exactly of all `lambda>0` with

\[
\mu=a/\lambda,\quad \alpha=c\lambda,\quad
\beta=c\lambda-a/\lambda-b\ge0.                         \tag{5}
\]

For `c>0`, the condition is
`lambda >= (b+sqrt(b^2+4ca))/(2c)`. For `c=0`, feasibility requires `b<0`
and `lambda>=a/(-b)`. The no-division boundary is included.

**Proof.** Multiply the last-stage equation by `lambda` to obtain
`h'=a+bh-ch^2`, `h(0)=0`. The triple uniquely determines its solution. In
the other direction,

\[
h'(0)=a,\quad h''(0)=ab,\quad h'''(0)=ab^2-2ca^2,
\]

so `a`, `b=h''/h'`, and
`c=(ab^2-h''')/(2a^2)` are recoverable. Division is legitimate because
`a>0`. Solving the definitions gives (5); multiplying its last inequality
by positive `lambda` gives the stated quadratic condition. For `c=0`, that
condition becomes `-a/lambda-b>=0`, requiring negative `b`. These are all
restrictions, so the description is both necessary and sufficient. **QED.**

### Theorem 3 — Proved many-stage limit and complete spectral invariants

**Assumptions:** Finite `m>=1`; B2–B5 for one clone; any specified conversion
flags; fixed finite positive intermediate mutation rates and finite
nonnegative birth/death rates. Replace the final mutation by an independent
terminal mark of rate `epsilon nu` per terminal cell, where `nu>0`, and
initiate clones at rate `lambda/epsilon`, `lambda>0`, `epsilon>0`. Marking
does not otherwise alter the unmarked process before its first mark. Observe
directly as in O0. Statements of convergence are on a fixed finite `[0,T]`.

Let `Z_m(t)` be the terminal-intermediate population of the **unmarked**
single-clone process and `Y_t=nu integral_0^t Z_m(u)du`. Define

\[
r_i=\alpha_i-\beta_i-c_i\mu_i\ (i<m),\quad
r_m=\alpha_m-\beta_m,\quad
A=\lambda\nu\prod_{i=1}^{m-1}\mu_i>0.
\]

Then `h_epsilon` converges uniformly to a genuine limiting first-event hazard
`h_0`, with explicit bound

\[
0\le h_0(t)-h_\epsilon(t)
\le\frac{\lambda\epsilon}{2}E[Y_T^2],\qquad
\mathcal L h_0(s)=\frac{A}{s\prod_{i=1}^{m}(s-r_i)}.
\tag{6}
\]

The limiting model has `k=m+1`. Its complete hazard identifies **exactly**
`k`, `A`, and the unordered multiset `{r_1,...,r_m}`, with multiplicities.
Any feasible biological parameter sets with these invariants are equivalent;
individual rates, stage order, and the number of stages with positive birth
rate generally are not identifiable. The `k=1` boundary is `h_0=A`.

**Proof.** Conditional on the unmarked clone path, the probability of no mark
by `t` is `exp(-epsilon Y_t)`. Equation (3)'s Poisson argument gives
`h_epsilon=(lambda/epsilon) E[1-exp(-epsilon Y_t)]`. For every `x>=0`,
`0<=x-(1-e^{-x})<=x^2/2`: integrate `0<=1-e^{-u}<=u` from zero to `x`.
Taking expectations proves (6)'s error bound with `h_0=lambda E[Y_t]`.
`Y_t` is nondecreasing, so the same upper bound works for every `t<=T`.

For completeness the moment is finite. Total cell count is dominated by a
Yule process starting at one cell with per-cell rate
`B=max_i(alpha_i+mu_i)` (omit the final mark; an overestimate is harmless).
Deaths and conversions do not increase total count. A Yule process has
`E[N(T)^2]=2e^{2BT}-e^{BT}`, obtained from
`(E N)'=B E N` and `(E N^2)'=2B E N^2+B E N`, with both initial moments
one. Thus
`E[Y_T^2]<=nu^2 T^2(2e^{2BT}-e^{BT})`, also valid at `B=0`.
Uniform convergence of hazards implies convergence of cumulative hazards
and survival functions on `[0,T]`, so this is a proved stochastic limit.

The forward means of one unmarked clone satisfy
`z_1'=r_1 z_1`, `z_1(0)=1`, and
`z_i'=r_i z_i+mu_{i-1}z_{i-1}`, `z_i(0)=0`. Their Laplace transforms,
followed by the terminal time integral, give the rational transform in (6)
for `s>max(0,r_1,...,r_m)`. Net rates `r_i` can be repeated, zero, or negative;
mutation rates remain strictly positive.

The denominator is monic and the numerator is a positive constant, so no
pole can cancel. Equality of two rational transforms forces equality of
their degrees, numerators, and monic denominator polynomials. Unique
factorization gives equality of the multisets including repeated zeros.
Conversely, those quantities make the transforms, and hence the continuous
hazards, equal. The shape is a convolution of nonnegative exponential
kernels and is positive at every `t>0`. Realizing any real `r_i` is possible
by choosing birth/loss rates with their specified difference, and positive
mutation factors with the stated product. No extra hidden restriction
splits an invariant class. This proves completeness for this limit. **QED.**

This is not completeness for finite `epsilon` nonlinear MSCE: higher-order
terms can contain additional information. Neither a finite-epsilon fit nor
the SEER table proves that the limiting regime applies. The spectral
symmetries are closely related to prior expected-count work [11].

### Theorem 4 — Exactly when exponential diagnosis delays hide the count

**Assumptions:** The limiting family in Theorem 3 (including `k=1`), free
`A>0`, OD with unknown finite `d>=0` and strictly positive finite delay
rates, exact observation on an open age interval. Each biological net rate
belongs to a specified allowed set: either all real rates or only
nonnegative rates. There are no other observation effects.

The complete observable invariants are

\[
C=A\prod_{j=1}^d\ell_j,\qquad
V=\{r_1,\ldots,r_{k-1},-\ell_1,\ldots,-\ell_d\},\qquad
\mathcal L h_{\rm obs}(s)=\frac{C}{s\prod_{v\in V}(s-v)}.
\tag{7}
\]

First remove **one** compulsory zero pole for the reservoir. Let `n_+` count
the remaining nonnegative poles, **including zero with multiplicity**, and
let `n_-` count negative poles.

* If all biological net rates are nonnegative, `k=1+n_+` is identifiable,
  as are the biological net-rate multiset, delay-rate multiset, and `A`.
* If signed biological net rates are allowed, the compatible biological
  counts are **exactly every integer** from `1+n_+` through
  `1+n_++n_-`. Choose any submultiset of negative poles to be delays;
  assign the rest to biological stages, and set `A=C/prod ell`.

Thus in the signed family `k` is identifiable if and only if `n_-=0`.
If a separately imposed convention excludes `k=1`, intersect the displayed
set with the allowed counts; it may be empty. Fixing the number of delays
would instead identify `k` from total degree, even when their assignments
remain ambiguous.

**Proof.** The clone-level lag convolves its first-hit-time law with the lag
law; the same Poisson marking argument convolves the hazard with the lag
density. Each exponential contributes the transform `ell/(s+ell)`.
Multiply these transforms with (6) to obtain (7). Its positive constant
numerator prevents cancellation. Monicity and unique factorization identify
`C` and `V` as in Theorem 3. Strictly positive delay rates produce only
negative poles. With nonnegative biological rates, signs therefore separate
the two classes uniquely. With signed biological rates, any negative pole
can be assigned to either class. Every such assignment is feasible by the
rate construction in Theorem 3 and positive rescaling of `A`; every
assignment has the same rational transform. No nonnegative pole can be a
permitted delay, so these are all possibilities. **QED.**

The result counts the order of a specific filter, not the order of arbitrary
tumor growth or detection mechanisms. General phase-type delays can have
numerator zeros and more complicated equivalences. No claim about them is
made here; phase-type representation nonuniqueness has extensive prior work
[28]. The role of detection has substantial precedent [17,26].

### Theorem 5 — Exact mutation/diagnosis equivalence and fast-stage closure

**Assumptions:** A downstream clone has an arbitrary time-homogeneous
first-event distribution `F(t)` (possibly defective through extinction),
with `F(t)=0` for `t<0`. Independent clones arrive at constant finite rate
`lambda>0`; hence `h(t)=lambda F(t)`. This includes B at every finite
parameter value and permits the `k=1` boundary with an immediate event.
Take finite `L>0`. An added precursor is nondividing, has no loss, and its
sole event is **conversion** to one downstream starting cell at rate `L`.

The `(k+1)`-transition model with this precursor has exactly the same
observed hazard as the `k`-transition model with OD lag `Exp(L)`:

\[
h_L(t)=\int_0^t L e^{-Ls}h(t-s)\,ds.                      \tag{8}
\]

For finite `T>=0`,

\[
0\le h_L(t)\le h(t),\qquad
\int_0^T(h-h_L)\,dt
=\frac{h_L(T)}L
=\int_0^T e^{-L(T-u)}h(u)\,du
\le\frac{h(T)(1-e^{-LT})}L\le\frac{h(T)}L.              \tag{9}
\]

The constant one in the last bound cannot be uniformly reduced over this
family and all `LT`: use a constant direct-event hazard and let `LT` grow.
For a downstream member of B, or any hazard with a finite positive leading
monomial, every finite `L` increases the direct-observation onset order by
one. Arbitrary CDFs need not have such an onset order. Exact equality between
the two *direct* biological hazards is not asserted.

**Proof.** A converting precursor waits an independent exponential time and
then starts exactly one downstream clone. Its first-event time is the sum
of that waiting time and a downstream first-event time. OD attaches the
same independent waiting-time law after the clone's first event; sums have
the same distribution. Multiplication of the CDF by `lambda` proves (8).
This reasoning fails for an asymmetric precursor, which can create many
downstream cells, and for lagging only the person's first malignant event.

`F`, and thus `h`, is nondecreasing; extending it by zero to negative ages
shows `h_L(t)=E[h(t-D)]<=h(t)`. Fubini for nonnegative integrands gives
`integral_0^T h_L(t)dt=integral_0^T h(u)(1-e^{-L(T-u)})du`.
Subtract to obtain the middle identity and bound the remaining `h(u)` by
`h(T)`. Formula (8) gives the equality to `h_L(T)/L`. At `T=0` both
integrals vanish. Convolution of the positive leading term in (3) with
`L+O(t)` increases its power by exactly one. **QED.**

### Theorem 6 — Explicit lower bounds for finite age-bin count experiments

**Assumptions:** OP with finitely many disjoint bins in `[0,T]`, finite
Poisson means, and two specified model orders `k` and `k+1`; any estimator
of `k`, including randomized estimators. For the closure specialization use
Theorem 5 and `0<=w_j(t)<=W<infinity` on each bin. For an exposure multiplier
`n`, replace `W` by `nW`. Parameters may vary with `n` when taking a minimax
supremum; the baseline hazard is fixed and nonzero.

Write `m_j` and `m'_j` for the two mean vectors and set

\[
D=\sum_j(\sqrt{m_j}-\sqrt{m'_j})^2,\qquad
d=\sum_j|m_j-m'_j|.
\]

The worst of the two probabilities of estimating the wrong order is at least

\[
\max\left\{\frac{1-\sqrt{1-e^{-D}}}{2},\quad
                  \frac{e^{-d}}2\right\}.                \tag{10}
\]

For the fast-precursor pair, `m'_j<=m_j` and

\[
d\le\frac{Wh(T)}L,\qquad
\max\{P_k(\widehat k\ne k),P_{k+1}(\widehat k\ne k+1)\}
\ge\tfrac12\exp[-Wh(T)/L].                               \tag{11}
\]

Consequently no estimator is uniformly consistent over families permitting
an unbounded extra converting rate. For example, with exposure `n`, choose
`L_n/n -> infinity`; the lower bound tends to `1/2`. This does not preclude
pointwise consistency at a fixed separated finite-rate pair.

**Proof.** For probability laws `P,Q`, let `A` be the event that the
estimator reports `k`. The sum of its two error probabilities is at least
`1-P(A)+Q(A)>=1-TV(P,Q)`, so their maximum is at least half this value.
For scalar Poisson laws the affinity is

\[
\sum_{z=0}^{\infty}\sqrt{e^{-m}m^z/z!\ e^{-m'}(m')^z/z!}
=\exp[-(m+m')/2+\sqrt{mm'}].
\]

Independence multiplies affinities, giving `R=exp(-D/2)` for all bins.
Cauchy–Schwarz applied to
`|p-q|=|sqrt(p)-sqrt(q)|(sqrt(p)+sqrt(q))` yields
`TV<=sqrt(1-R^2)`, including zero means by continuity. This proves the first
bound. Alternatively couple each bin by a common Poisson variable of mean
`min(m,m')` plus independent excess Poisson variables. With probability
`exp(-d)` every excess count is zero. Thus `TV<=1-exp(-d)`, proving the
second bound. This latter coupling need not be optimal, which is why both
bounds are reported.

For (8), nonnegativity of deficits and disjoint bins imply
`d=sum integral_Bj w_j(h-h_L)<=W integral_0^T(h-h_L)`.
Use (9). For every fixed finite exposure the larger class contains all
finite `L`, so taking `L` arbitrarily large proves a minimax lower bound of
at least `1/2`. Under increasing exposure, admissible `L_n/n -> infinity`
gives the same obstruction. Every chosen `L_n` is finite and positive, so
this is not merely substituting an infinite-rate parameter. **QED.**

The hypotheses of (11) express a missing time-scale restriction. They do
not prove that actual cancer has arbitrarily fast biological steps. A lower
bound for a specified bounded-rate subclass would need that restriction and
the observation mechanism in the experiment.

### Theorem 7 — Sharp calibration frontier and exact finite-bin confounding

**Assumptions:** OC, finitely many positive shape entries `f_j`, a free
amplitude `A>0`, and positive reported-mean intervals
`0<l_j<=u_j<infinity`. Intervals may be degenerate exact values. They are
deterministic sets; they are not automatically sampling confidence intervals.

Let `L_f=max_j(l_j/f_j)` and `U_f=min_j(u_j/f_j)`. There exist an amplitude,
admissible age multipliers, and values `y_j in [l_j,u_j]` **if and only if**

\[
\rho^2\ge R_f:=\max\{1,L_f/U_f\}.                       \tag{12}
\]

The sharp minimum is `rho_*=sqrt(R_f)` and the allowed amplitude interval
is `[L_f/rho, rho U_f]`. In log units the minimum distortion is
`delta_*=max(0, (log L_f-log U_f)/2)`. For exact data, the attaining
amplitude at a nontrivial boundary is
`sqrt(min_j(y_j/f_j) max_j(y_j/f_j))`.

If each shape is itself allowed anywhere in a closed positive interval
`[f_j^-,f_j^+]`, replace `L_f` by `max(l_j/f_j^+)` and `U_f` by
`min(u_j/f_j^-)`. For a monotone continuous shape and unknown age weights
on a half-open age bin, endpoint ranges give a sharp **infimum over the
closure**; endpoint attainment by actual age weights is not guaranteed.

For a fixed nonzero exact branching baseline, finite bins of OP with
`m_*:=min_j m_j>0`, bounded exposure density `W`, and any `rho>1`, there is
a finite-rate extra-step model with exactly the same calibrated bin means.
One sufficient choice is

\[
L\ge \frac{Wh(T)}{m_*(1-\rho^{-2})}.                    \tag{13}
\]

The baseline can use multiplier one. The extra-step model uses a rescaled
initiation amplitude and multipliers in `[1/rho,rho]`. Thus a positive
calibration allowance, however small, destroys finite-bin structural
identification over this unbounded-rate family, even at infinite precision
in those bin means. With no finite bound on `rho`, **every** positive shape
and therefore every positive integer `k` is compatible with any positive
finite rate vector after amplitude/multiplier adjustment.

**Proof.** For each bin, the possible rate interval is
`[A f_j/rho,A f_j rho]`. It intersects `[l_j,u_j]` exactly when
`l_j/(rho f_j)<=A<=rho u_j/f_j`. Intersecting all amplitude intervals gives
the displayed interval, nonempty exactly when `L_f<=rho^2 U_f`. Also
`rho>=1`, giving (12). At equality the interval is a singleton if
`L_f>U_f`; if the uncalibrated intervals already intersect, `rho=1`
suffices. This proves necessity, sufficiency, and sharpness without
differentiability or asymptotics.

Allowing `f_j` to vary fills the continuous product interval
`[A f_j^-/rho,A f_j^+ rho]`, giving the interval extension by the same
argument. A continuous strictly monotone function has any interior endpoint
average under a suitable within-bin probability distribution. Its excluded
upper endpoint may only be approached on a half-open bin; taking the
closure gives the stated infimum rather than false guaranteed attainment.

For (13), Theorem 5 bounds each positive deficit by
`Delta=Wh(T)/L`. Hence
`m_{Lj}>=m_j-Delta>=m_j(1-Delta/m_*)`. Under (13),
`1<=m_j/m_{Lj}<=rho^2`. The ratio between the largest and smallest such
ratios is at most `rho^2`, so apply (12) to exact data `m_j` and shape
`m_{Lj}`. Changing `lambda` rescales a branching hazard without changing
its clone law, providing the required amplitude. All means then agree
exactly, so their OP count laws also agree. Finally, every finite positive
ratio vector has finite maximum/minimum; with unrestricted `rho`, (12)
always holds. Pure-power limiting shapes, or any chosen positive exact
branching shapes for each `k`, supply all counts. **QED.**

For exact positive vectors, twice the minimum log distortion is the classical
Hilbert projective distance: `log(max(y/f)/min(y/f))` [29]. Equivalently,
one centers the log ratios in the maximum absolute norm. This formula and
elementary interval intersection are not claimed new. The candidate
contribution is their composition with (8)–(13) for the stated cancer
incidence observation experiment.

The one-bin case is already completely ambiguous with a free amplitude,
even at `rho=1`. We make no claim that `rho=1` identifies `k` from an
arbitrary finite set of bins. Formula (12) is a sharp frontier for a
specified shape, while (13) is a sufficient general construction, not a
claim to be the smallest possible `L`.

### Theorem 8 — A universal monotonicity falsifier and its sharp relaxation

**Assumptions:** Constant Poisson initiation of independent identically
distributed clones as in Theorem 5, optional independent clone-level lag,
and direct bin means with any nonnegative normalized age weights supported
on ordered nonoverlapping bins. To quantify relaxation, use OC and positive
reported-mean intervals `[l_j,u_j]` in bin order.

The underlying hazard is nondecreasing regardless of stage count or
expansion. Therefore every later bin mean is at least every earlier bin
mean. Among **all nondecreasing positive mean vectors**, the sharp smallest
calibration factor is

\[
\rho_{\rm mono}^2=
\max\left\{1,\max_{i\le j}\frac{l_i}{u_j}\right\}.      \tag{14}
\]

This is a necessary compatibility threshold for the branching family, not
a sufficient condition for realizability by a particular MSCE process.

**Proof.** A CDF is nondecreasing, so `h=lambda F` is nondecreasing. A lag
replaces `F` by another CDF. Every point in a later bin lies after every
point in an earlier bin, so integration against any normalized nonnegative
weights preserves the ordering. With calibration, an admissible underlying
mean `z_j` must lie in `[l_j/rho,rho u_j]`. A nondecreasing selection exists
if and only if `l_i/rho<=rho u_j` for all `i<=j`: necessity follows from
`z_i<=z_j`; sufficiency follows by taking
`z_j=max_{i<=j} l_i/rho`. This construction proves (14) including equality.
Not every such vector is a branching-model curve, so sufficiency is only
for the stated shape-constrained class. **QED.**

## 4. Published aggregate incidence: data, method, and results

### 4.1 Provenance and rights

The small CSV contains 44 published summary rates from the **SEER Cancer
Statistics Review 1975–2017**, November 2019 submission, released April 15,
2020 [27]. The chosen incidence period is **2013–2017**, SEER 21 areas, all
races, invasive sites. Ages 30–34 through 80–84 provide eleven complete
five-year bins per site. Both sexes are used except female breast.

| Site | Source table | Column | Source |
|---|---|---|---|
| Pancreas | 22.7 | All races, total | [Table PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_22_table.07.pdf) |
| Colon and rectum | 6.10 | All races, total | [Table PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_06_table.10.pdf) |
| Lung and bronchus | 15.10 | All races, total | [Table PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_15_table.10.pdf) |
| Female breast | 4.11 | All races, females | [Table PDF](https://seer.cancer.gov/archive/csr/1975_2017/results_single/sect_04_table.11.pdf) |

Each CSV row carries source, table, section, column, HTML/PDF/download URLs,
access date **2026-11-17**, and reuse-permission URL. The archive landing
page explicitly puts the report in the public domain and permits
reproduction [27]. This meets `DATA_POLICY.md` for small aggregate
redistribution; the repository license does not relicense source material.
No individual records, linked patient data, genomes, or controlled-access
data are present. `data/PROVENANCE.md` records selection, all 44 printed
values, missing fields, and transcription checks. HTML and PDF text were
cross-checked; they are two presentations of one source, not independent
biological evidence. The inaccessible original CSV response was not
silently claimed as a transcription source.

`case_count` and `person_years` are blank because the selected tables do not
publish them. We do not reverse-engineer an invented denominator from a
rounded rate. No SEER-specific confidence interval, p-value, or numerical
sampling-error lower bound is claimed. The dated archive is a reproducible
historical source, not a statement of current cancer rates.

### 4.2 Reproducible shape/calibration analysis

For display precision of one decimal, the analysis **assumes** a closed
rounding interval of `rate +/- 0.05` per 100,000. This is sensitivity to a
round-to-nearest convention, not a confidence interval or a guarantee about
source rounding. It compares two limiting subfamilies with free amplitude:

\[
f_{k,0}(t)=t^{k-1}/(k-1)!,\qquad
f_{k,g}(t)=\sum_{n=0}^{\infty}
                   \frac{g^n t^{n+k-1}}{(n+k-1)!}.       \tag{15}
\]

The second has one net-growth rate `g` and `k-2` zero net rates. Negative
`g` means net loss, not positive net expansion. The average on `[a,b)` is

\[
\overline f_{k,g}=
\sum_{n=0}^{\infty}
\frac{g^n(b^{n+k}-a^{n+k})}{(b-a)(n+k)!}.                 \tag{16}
\]

`analysis.py` applies (12) using exact rational data, bins, growth-grid
points, and rational Taylor enclosures. The declared grid is `k=2,...,10`
and `g=-0.060,...,0.060` per year in increments of `0.005`. The zero-growth
model is included. These are **finite-grid compatibility minima**, not a
continuous optimum, likelihood fit, stage-count estimate, biological prior,
or optimum over all multistage models. Lower thresholds could occur off
grid or with several nonzero growth rates. Common amplitude absorbs rate
units. Default within-bin averaging is uniform; actual SEER person-time
weights are unavailable. The output separately gives the sharp
endpoint-closure infimum for pure powers with freely varying within-bin
weights. The monotonicity result does not need uniform weights.

Every entry below is computed from the checked-in CSV, and the full per-`k`
grid is in `analysis_results.json`.

| Site | Monotonicity under exact rounded means | Best pure-power grid `k`; minimum `rho` | Best one-growth grid `k`, `g/year`; minimum `rho` | Identified biological `k` |
|---|---|---|---|---|
| Pancreas | Passes necessary condition | 6; 1.333798 | 6, +0.030; 1.245370 | None |
| Colon and rectum | Passes necessary condition | 5; 1.265066 | 5, -0.010; 1.241860 | None |
| Lung and bronchus | Passes necessary condition | 7; 1.576755 | 7, +0.040; 1.482398 | None |
| Female breast | Fails homogeneous monotone-mean condition | 4; 1.578568 | 4, -0.015; 1.504293 | None |

Table thresholds are rounded to six decimals for display and are not exact
feasibility cutoffs. For interpretation, approximately `rho=1.245370` means
the chosen shape needs symmetric relative multipliers approximately in
`[1/1.245370,1.245370]`, after optimizing common amplitude. Feasibility holds
at the unrounded threshold. It does not mean
SEER has that much underreporting. The grid's best `k` is merely a shape
index under those assumptions. Pancreatic `k=7` also has a low grid
threshold, 1.271340; lung `k=8` gives 1.506873. A best grid index is not a
certificate that alternatives are impossible.

Female breast rates are 470.2, 458.9, and 413.5 at ages 70–74, 75–79,
and 80–84, respectively [Table 4.11]. Even the rounding intervals do not
overlap in an order compatible with monotonicity. Formula (14) gives
`rho_mono=1.066238`, the sharp infimum/attained minimum for the abstract
nondecreasing **mean-vector** class. This result is descriptive: random
sample rates can decrease even under a monotone true hazard. Without the
missing denominators we cannot report a significance level. Period rates
also mix cohorts, exposures, screening, competing risks, and subtypes;
failure of this homogeneous observation model is not failure of multistage
carcinogenesis.

**Which cancers have an identifiable count?** In the unrestricted OC
sensitivity class, **none of these four**: Theorem 7 gives exact
observation equivalence for every positive `k`. In the uncalibrated,
constant-initiation model, the breast point summaries fail the necessary
monotone-mean check; the remaining three pass but are not thereby identified.
Under O0 and B, *every* cancer hypothetically satisfying those assumptions
would have an identifiable model order by Theorem 1. These different
answers reflect observation assumptions, not a cancer-specific discovery of
uniqueness or nonuniqueness from eleven rounded values.

### 4.3 Explicit finite-exposure numbers, labeled synthetic

For illustration only, take a specified exposure density bound `W=100000`
persons and an assumed endpoint hazard bound `h(T)<=0.001/year`. These
numbers are **synthetic**, not recovered SEER exposures or fitted hazards.
Theorem 6 gives the following worst-case error lower bounds for `k` versus
one extra converting stage:

| Extra conversion rate `L/year` | Bound `exp(-100/L)/2` |
|---:|---:|
| 100 | 0.18393972 |
| 1000 | 0.45241870 |
| 10000 | 0.49502491 |

These displayed lower bounds are rounded **downward**, preserving their
inequality direction. The independent Poisson simulation uses a different, fully specified
three-bin limiting pair, reports both mean vectors and both bounds, and
compares a likelihood-ratio classifier with the bound. It is not fitted to
the real CSV and is not evidence of actual biological rates.

### 4.4 Comparison with sequencing-based driver estimates

| Study and estimand | Published result | What the incidence analysis can compare |
|---|---|---|
| Martincorena et al. (2017), dN/dS selected coding substitutions [20] | Approximately four per tumor overall, with the abstract reporting fewer than one in thyroid/testicular and more than ten in endometrial/colorectal tumors. | Exome-wide selected substitutions include a different target than rate-limiting initiation transitions. The colorectal shape index 5 is not evidence against this count. |
| Iranzo, Martincorena & Koonin (2018), known-gene driver network [21] | Approximately two drivers on average, with type estimates spanning one to five in the primary 198-gene analysis. | A restricted gene set and a different model change the estimand; this is not the same quantity as exome-wide dN/dS. |
| Anandakrishnan et al. (2019), mutation-count-distribution model [22] | Two to eight inferred hits across 17 types; colon and lung adenocarcinoma estimates are three. | Our site-wide colon/rectum and lung/bronchus tables do not match those histological cohorts. Our grid indices 5 and 7 do not identify alternative driver counts. |
| Tomasetti et al. (2015), incidence ratios plus mutations/exposure strata [23] | Three inferred drivers for lung and colorectal adenocarcinoma under their framework. | Additional strata and molecular information are doing work; this is not single-curve identification. |
| Houlahan et al. (2026), carrier contrasts plus genomic evidence [18] | Inferred one to three deletion events in breast/ovarian tumorigenesis under their framework. | Deletions, inherited susceptibility, and additional observations are not represented by our site-level shape index. |
| Bandlamudi et al. (2026), panel-defined driver alterations [19] | Median four alterations per tumor in a large clinical sequencing cohort. | Observed selected alterations can arise after initiation; this cannot identify our `k`. |

Tomasetti & Vogelstein's **Science 2015** stem-cell-divisions paper [24]
concerns an across-tissue association with lifetime risk. It is distinct from
Tomasetti et al.'s **PNAS 2015** driver paper [23] and does not invert an
age-specific hazard to obtain `k`. The cited sequencing cohorts are not
linked to SEER records here. No patient-level sequencing data are used.

## 5. Verification, falsification attempts, and hostile review

### 5.1 What each implementation checks

`verifier.py` uses only the Python standard library. `Fraction` supplies
exact Taylor jets, products, invariants, inequalities, and interval bounds.
`Decimal` only presents square roots/logarithms/exponentials; Decimal output
is **not called exact arithmetic**. Exp and shape enclosures use exact
rational Taylor sums and geometric tail bounds. For instance, after term
`N` in (16), the absolute tail is at most

\[
\frac{|g|^{N+1} b^{N+1+k}}{(b-a)(N+1+k)!}
\frac1{1-|g|b/(N+k+2)},\qquad |g|b<N+k+2.               \tag{17}
\]

For every term, `b^p-a^p<=b^p` bounds its absolute value by the majorant
`B_n=|g|^n b^(n+k)/((b-a)(n+k)!)`. Starting at `n=N+1`, successive
**majorant** ratios are bounded by `|g|b/(N+k+2)`, giving (17). Actual
bin-integral term ratios need not satisfy that bound. This argument also
encloses alternating series at negative `g`.
The CSV provenance SHA-256 in `analysis_results.json` hashes UTF-8 text
**after universal-newline normalization to LF**. This makes the digest
stable for LF and Windows CRLF checkouts without rewriting the original
SEER values; `test_10` also checks a CR-only input. The digest is a local
content-integrity checksum, not verification of the upstream data.

`analysis.py` uses these enclosures to certify grid comparisons before
formatting any rounded output, and asserts its squared-threshold enclosure
width is below `10^-24`. This is numerical certification, not proof that
the source numbers are true or that the models apply biologically.

`simulation.py` is an independent implementation: it imports neither the
verifier nor the analysis. It uses event-level Gillespie paths, ordinary
double-precision RK4, sampled unmarked cell-time integrals, simplex
integration by Monte Carlo, and Poisson likelihood classification. The
seed is fixed at `19270617`; tolerances are six stated standard errors or
explicit numerical residual bounds. It reports sample sizes and residuals.
There is no adaptive seed selection or claim that passing random examples
proves a theorem. Independent implementation is not independent biological
validation, and agreement between AI contributors is not peer review.

| Result | Exact checks | Independent numerical/simulation checks |
|---|---|---|
| Theorem 1 | 128 rational parameter cases across `k=1,...,8`, mixed mutation conventions and zero birth/death boundaries | Gillespie versus backward RK4 for `k=2,3,4` with expansion at all intermediates; small-age scaling |
| Theorem 2 | Full jet equality for different feasible parameter sets, derivative recovery, no-division boundary | RK4 versus TSCE closed form, transformed biological rates, finite-difference derivatives |
| Theorem 3 | Permutation-invariant polynomial signatures and independent linear Taylor recurrence | Unmarked path-integral rare-mark limit and bound; simplex Monte Carlo versus forward mean ODE |
| Theorem 4 | Exact equal signatures with a negative pole assigned to biology or diagnosis | Every allocation of two negative poles, repeated/zero poles, `k=1` all-delay boundary |
| Theorem 5 | Added converting stage equals the independent scalar delay filter in exact Taylor coefficients; asymmetric-stage counterexample | Gillespie precursor versus clone-level diagnosis lag; integrated deficit and monotonicity |
| Theorem 6 | Rational Poisson affinity partial sums with tail enclosures and synthetic exponential-bound enclosures | Fixed-seed likelihood-ratio decision experiment with stated Poisson means |
| Theorem 7 | Rational feasibility at, below, and above threshold; interval-data and interval-shape boundaries | Independent interval intersection with random positive shapes; limiting close-order checks |
| Theorem 8 | Exact rounded-value inequalities and optimal monotone envelope | Nondecreasing simulated/ODE hazards; independent source-rate reversal audit |

### 5.2 Deliberate attempts to break the claims

The self-review retained the following counterexamples and limitations rather
than explaining them away:

* An extra **asymmetric** fast-mutation precursor does not obey (8). The exact
  jets differ. The theorem explicitly requires a converting precursor.
* Equal mean-count spectra do not establish equality of finite-rate
  nonlinear branching hazards. Theorem 3 uses an explicit stochastic limit
  with an error bound; the general finite-rate MSCE parameter classification
  is not claimed.
* Repeated rates, zero growth, negative growth, zero birth, extinction,
  `k=1`, zero Poisson means, one bin, and `rho=1` were examined. Strictly
  zero mutation/amplitude and zero delay-rate boundaries invalidate selected
  divisions or reduced-order arguments and are excluded where necessary.
* A diagnostic lag on the person's first onset does not generally convolve
  hazards. OD is a marked immigrant-clone construction, with its restrictions
  stated explicitly.
* Unknown within-bin age weights do not necessarily attain endpoint values.
  Endpoint results are labeled closure infima when required.
* An initial implementation ranked grid winners after formatting to six
  decimal places. Hostile code review caught this; selection now uses exact
  rational bounds and formats only afterward.
* Rounding bounds are not statistical uncertainty. Blank denominators stay
  blank, and the synthetic exposure experiment is kept separate.
* Female breast violates the model's monotone-mean condition; the failure is
  reported. A good fit to any of the other sites would not validate the
  biological assumptions or make the best grid index identifiable.

No unsupported mathematical conjecture is presented as a theorem. The open
priority question is explicitly a **novelty hypothesis**: a targeted search
did not locate the same combined finite-bin/certified-calibration statement.
Theorem proofs and code may be assessed independently of that hypothesis.
General nonlinear MSCE global identifiability, continuous global fitting,
and identification with unrestricted diagnosis/mixture mechanisms remain
outside the proved claims.

### 5.3 Reproduction and repository scope

From the repository root, with Python 3.10 or later and no extra packages:

```bash
python research/multistage-incidence-identifiability/verifier.py
python research/multistage-incidence-identifiability/simulation.py
python research/multistage-incidence-identifiability/analysis.py
python -m unittest discover -s tests -p 'test_research_multistage_incidence_identifiability.py' -v
```

To regenerate the checked-in deterministic analysis output:

```bash
python research/multistage-incidence-identifiability/analysis.py --output research/multistage-incidence-identifiability/analysis_results.json
```

All tests are offline and stdlib-only. The single new test file exercises the
exact algebra, independent simulation, data provenance/integrity, and complete
CSV-to-result reproduction within the requested 60-second budget. The
repository's existing validation workflow discovers it without any workflow
change. The existing repository validation and test gates are also run.

Only this new research directory and
`tests/test_research_multistage_incidence_identifiability.py` are added.
No existing research, catalog, results, workflow, policy, or root README is
modified. This is an unreviewed research contribution, not a promotion into
the repository's editorial/publication registry.

## 6. References

Links identify primary sources; access date for the literature review is
2026-11-17. The report summarizes literature without copying its tables,
figures, or full text. The incidence data have separate explicit rights and
per-row attribution in the CSV.

1. Nordling CO (1953). *A new theory on cancer-inducing mechanism.* Br J Cancer 7:68–72. [doi:10.1038/bjc.1953.8](https://doi.org/10.1038/bjc.1953.8).
2. Armitage P, Doll R (1954). *The age distribution of cancer and a multi-stage theory of carcinogenesis.* Br J Cancer 8:1–12. [doi:10.1038/bjc.1954.1](https://doi.org/10.1038/bjc.1954.1).
3. Armitage P, Doll R (1957). *A two-stage theory of carcinogenesis in relation to the age distribution of human cancer.* Br J Cancer 11:161–169. [PubMed 13460138](https://pubmed.ncbi.nlm.nih.gov/13460138/).
4. Moolgavkar SH, Venzon DJ (1979). *Two-event models for carcinogenesis: incidence curves for childhood and adult tumors.* Math Biosci 47:55–77. [doi:10.1016/0025-5564(79)90005-1](https://doi.org/10.1016/0025-5564(79)90005-1).
5. Moolgavkar SH, Knudson AG Jr (1981). *Mutation and cancer: a model for human carcinogenesis.* J Natl Cancer Inst 66:1037–1052. [PubMed 6941039](https://pubmed.ncbi.nlm.nih.gov/6941039/).
6. Little MP (1995). *Are two mutations sufficient to cause cancer? Some generalizations of the two-mutation model of carcinogenesis of Moolgavkar, Venzon, and Knudson, and of the multistage model of Armitage and Doll.* Biometrics 51:1278–1291. [PubMed 8589222](https://pubmed.ncbi.nlm.nih.gov/8589222/).
7. Heidenreich WF (1996). *On the parameters of the clonal expansion model.* Radiat Environ Biophys 35:127–129. [doi:10.1007/BF02434036](https://doi.org/10.1007/BF02434036).
8. Heidenreich WF, Luebeck EG, Moolgavkar SH (1997). *Some properties of the hazard function of the two-mutation clonal expansion model.* Risk Anal 17:391–399. [doi:10.1111/j.1539-6924.1997.tb00878.x](https://doi.org/10.1111/j.1539-6924.1997.tb00878.x).
9. Luebeck EG, Moolgavkar SH (2002). *Multistage carcinogenesis and the incidence of colorectal cancer.* PNAS 99:15095–15100. [doi:10.1073/pnas.222118199](https://doi.org/10.1073/pnas.222118199).
10. Meza R, Jeon J, Moolgavkar SH, Luebeck EG (2008). *Age-specific incidence of cancer: phases, transitions, and biological implications.* PNAS 105:16284–16289. [doi:10.1073/pnas.0801151105](https://doi.org/10.1073/pnas.0801151105).
11. Cox LA Jr, Huber WA (2007). *Symmetry, identifiability, and prediction uncertainties in multistage clonal expansion (MSCE) models of carcinogenesis.* Risk Anal 27:1441–1453. [doi:10.1111/j.1539-6924.2007.00980.x](https://doi.org/10.1111/j.1539-6924.2007.00980.x).
12. Frank SA (2007). *Dynamics of Cancer: Incidence, Inheritance, and Evolution.* Princeton University Press. [NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK1568/).
13. Little MP, Heidenreich WF, Li G (2009). *Parameter identifiability and redundancy in a general class of stochastic carcinogenesis models.* PLOS ONE 4:e8520. [doi:10.1371/journal.pone.0008520](https://doi.org/10.1371/journal.pone.0008520).
14. Brouwer AF, Meza R, Eisenberg MC (2017; online 2016). *A systematic approach to determining the identifiability of multistage carcinogenesis models.* Risk Anal 37:1375–1387. [doi:10.1111/risa.12684](https://doi.org/10.1111/risa.12684), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5472511/).
15. Brouwer AF, Meza R, Eisenberg MC (2017). *Parameter estimation for multistage clonal expansion models from cancer incidence data: a practical identifiability analysis.* PLOS Comput Biol 13:e1005431. [doi:10.1371/journal.pcbi.1005431](https://doi.org/10.1371/journal.pcbi.1005431).
16. Richardson RB, Anghel CV, Deng DS (2021). *Profound synchrony of age-specific incidence rates and tumor suppression for different cancer types as revealed by the multistage-senescence model of carcinogenesis.* Aging (Albany NY) 13:23545–23578. [doi:10.18632/aging.203651](https://doi.org/10.18632/aging.203651).
17. Mirzaei NM et al. (2025). *Modeling early-onset cancer kinetics reveals changes in underlying risk and the impact of population screening.* Cancer Res 85:4558–4570. [doi:10.1158/0008-5472.CAN-25-0186](https://doi.org/10.1158/0008-5472.CAN-25-0186).
18. Houlahan KE et al. (2026). *Quantifying rate-limiting genetic variation in breast and ovarian tumourigenesis.* EBioMedicine 125:106181. [doi:10.1016/j.ebiom.2026.106181](https://doi.org/10.1016/j.ebiom.2026.106181).
19. Bandlamudi C et al. (2026). *Cancer type-specific variation in patterns of driver alterations across 50,000 tumors.* Cancer Cell 44:1045–1062.e6. [doi:10.1016/j.ccell.2026.03.003](https://doi.org/10.1016/j.ccell.2026.03.003).
20. Martincorena I et al. (2017). *Universal patterns of selection in cancer and somatic tissues.* Cell 171:1029–1041.e21. [doi:10.1016/j.cell.2017.09.042](https://doi.org/10.1016/j.cell.2017.09.042).
21. Iranzo J, Martincorena I, Koonin EV (2018). *Cancer-mutation network and the number and specificity of driver mutations.* PNAS 115:E6010–E6019. [doi:10.1073/pnas.1803155115](https://doi.org/10.1073/pnas.1803155115).
22. Anandakrishnan R, Varghese RT, Kinney NA, Garner HR (2019). *Estimating the number of genetic mutations (hits) required for carcinogenesis based on the distribution of somatic mutations.* PLOS Comput Biol 15:e1006881. [doi:10.1371/journal.pcbi.1006881](https://doi.org/10.1371/journal.pcbi.1006881).
23. Tomasetti C et al. (2015). *Only three driver gene mutations are required for the development of lung and colorectal cancers.* PNAS 112:118–123. [doi:10.1073/pnas.1421839112](https://doi.org/10.1073/pnas.1421839112).
24. Tomasetti C, Vogelstein B (2015). *Variation in cancer risk among tissues can be explained by the number of stem cell divisions.* Science 347:78–81. [doi:10.1126/science.1260825](https://doi.org/10.1126/science.1260825).
25. Zheng Q (1998). *To use or not to use? Backward equations in stochastic carcinogenesis models.* Biometrics 54:384–388. [PubMed 9544530](https://pubmed.ncbi.nlm.nih.gov/9544530/).
26. Luebeck EG, Curtius K, Jeon J, Hazelton WD (2013). *Impact of tumor progression on cancer incidence curves.* Cancer Res 73:1086–1096. [doi:10.1158/0008-5472.CAN-12-2198](https://doi.org/10.1158/0008-5472.CAN-12-2198).
27. Howlader N et al., editors (2020). *SEER Cancer Statistics Review, 1975–2017.* National Cancer Institute, November 2019 submission. [Archive and public-domain statement](https://seer.cancer.gov/archive/csr/1975_2017/). Table-specific links and access dates are included on every data row.
28. Lindqvist BH (2023; online 2022). *Phase-type models for competing risks, with emphasis on identifiability issues.* Lifetime Data Anal 29:318–341. [doi:10.1007/s10985-022-09547-7](https://doi.org/10.1007/s10985-022-09547-7).
29. Bushell PJ (1973). *Hilbert's metric and positive contraction mappings in a Banach space.* Arch Rational Mech Anal 52:330–338. [doi:10.1007/BF00247467](https://doi.org/10.1007/BF00247467).
