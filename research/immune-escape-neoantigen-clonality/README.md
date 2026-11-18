# Inherited antigen coverage and escape through loss of presentation

**Nonclinical computational mathematics. Unreviewed research note, 17 November
2026. No treatment advice, patient-level data, or biological validation.**
“Proved” below means a deductive result in the stated stochastic model. Exact
arithmetic and passing tests do not establish clinical validity or historical
priority. All numerical rates and simulated observations are synthetic.

## Result and literature-first pivot

The broad question already has substantial biological and mathematical
answers: clonal targets cover more tumour cells, antigen quality matters,
immune selection changes the neoantigen landscape, and presentation loss can
permit escape. We therefore **pivot from a claimed new explanation of
clonality to a narrower open-to-priority-check question**: what is the exact
large-antigen escape law when individual antigen loss and a common
presentation-loss event remain possible throughout branching growth?

The candidate contribution is Theorem 6. In an infinite antigen-type,
unbounded cell-population model, let each antigen contribute killing rate
`kappa`, let complete loss of the modeled presentation occur at fixed rate
`h > 0`, and let a presentation-null cell establish with probability `sigma`.
Subject to the explicitly stated source-suppression condition,

\[
 \Pr(\text{successful presentation-null escape}\mid k\text{ initial antigens})
 \sim \frac{\sigma h}{\kappa k},\qquad k\longrightarrow\infty.
\]

Both sides of a finite-parameter probability sandwich have this same leading
coefficient. For arbitrarily many starting cells with loads `k_i`, the sharp
scaling variable is the **sum of reciprocal loads**, provided the smallest
load also tends to infinity. This gives a quantitative limit on how increasing
clonal antigen burden suppresses one common escape route. The theorem permits
ongoing antigen acquisition, individual antigen loss, division and death; it
does not assume a tumour of a few cells or truncate the antigen state space.

Two equally important qualifications are proved, not hidden. A count-only
killing model cannot distinguish shared antigen identity from a disjoint
collection of private antigens giving exactly the same per-cell loads
(Theorem 5). Also, a positive rate of successful presentation-null escape
prevents almost-sure total eradication at every finite immune killing strength
(Theorem 1). Neither qualification negates a probability benefit from better
coverage.

### What the literature already establishes

The following primary-source review preceded selection of the theorem. Source
links identify originals; no source datasets or paper text are redistributed.

| Work | Prior result relevant here | Boundary of this note's claim |
| --- | --- | --- |
| [McGranahan et al., Science 2016](https://doi.org/10.1126/science.aaf1490) | Clonal neoantigen enrichment associates with checkpoint-blockade benefit; clonal-reactive T cells and a coverage/dosage interpretation. | The biological clonality explanation is established. |
| [McGranahan et al., Cell 2017](https://doi.org/10.1016/j.cell.2017.10.001) | Allele-specific HLA loss and evidence of selection during lung-cancer evolution. | Ordinary HLA LOH is not loss of every presentation pathway. |
| [Rosenthal et al., Nature 2019](https://doi.org/10.1038/s41586-019-1032-7) | Loss of previously clonal neoantigens, HLA loss and suppressed expression are escape mechanisms. | Both loss mechanisms are biological inputs, not discoveries here. |
| [Lakatos et al., Nature Genetics 2020](https://doi.org/10.1038/s41588-020-0687-1) | Branching growth, inherited antigenicity, antigen-dependent death, random escape, and treatment-driven selection of escaped clones. | The basic branching model architecture overlaps strongly. |
| [Aguadé-Gorgorió and Solé, JRS Interface 2020](https://doi.org/10.1098/rsif.2020.0736) | Neoantigen heterogeneity thresholds separating growth/control and comparison with melanoma cohorts. | This is not the first heterogeneity or immune-control threshold. |
| [George and Levine, J Theor Biol 2018](https://doi.org/10.1016/j.jtbi.2018.09.012); [eLife 2023](https://doi.org/10.7554/eLife.82786) | Stochastic recognition/evasion; antigen loss and replenishment, elimination, escape, and optimal evasion at recognition epochs. | Exact stochastic antigen escape calculations are not new in general. |
| [Łuksza et al., Nature 2017](https://doi.org/10.1038/nature24473); [Nature 2022](https://doi.org/10.1038/s41586-022-04735-9) | Neoantigen fitness based on presentation/recognition and quality-dependent clonal evolution. | Equal-quality additive counts omit established biological structure. |
| [Zapata et al., Nature Genetics 2023](https://doi.org/10.1038/s41588-023-01313-1) | Immune dN/dS and contrasting immune-edited/escaped evolutionary histories associated with treatment response. | Reversible immune evasion must not be equated with permanent presentation loss. |
| [Qi, Vincent and Cao, PLoS Comput Biol 2023](https://doi.org/10.1371/journal.pcbi.1010976) | Antigen-specific effectors, resource competition, inherited antigen gain/loss and MHC-allele loss. | Including both antigen and MHC loss is already modeled. |
| [Chen et al., Communications Biology 2024](https://doi.org/10.1038/s42003-024-06460-7) | Frequency-dependent selection, inherited neoantigens and escape, compared with response data. | Frequency-dependent immune recruitment is an existing alternative to our count-only model. |
| [Wang, Morison and Huang, eLife 2025](https://doi.org/10.7554/eLife.103970.3) | Stochastic cancer–effector coevolution, antigen accumulation, escape and elimination. | Our permanently hidden state differs from susceptibility restored by new antigens. |
| [Sholokhova, Kaveh and Bozic, Nature Communications 2026](https://doi.org/10.1038/s41467-026-71135-2) | Pretreatment branching evolution plus immune-treatment dynamics; strong clonal antigens predict durable response, and nearly complete subclonal coverage can sometimes suffice. | Strong clonal targets and near-complete subclonal coverage are existing explanations; the proposed contribution is a specific probability law. |

The mathematical ingredients also have direct antecedents:

* [Karlin and Tavaré (1982), *Linear birth and death processes with
  killing*](https://doi.org/10.2307/3213507),
  [author-hosted paper](https://www.columbia.edu/cu/simontavare/STpapers-pdf/KT82c.pdf):
  our immigration–death antigen-path transform is a specialization of their
  linear killing framework. The cell branching rate is separate from antigen
  birth in that framework.
* [Martin et al. (2013)](https://doi.org/10.1098/rstb.2012.0088) and
  [Azevedo and Olofsson (2021)](https://doi.org/10.1016/j.mbs.2021.108708):
  evolutionary-rescue probabilities, founder products and rare-rescue Poisson
  laws are established methods.
* [*On the waiting time to escape* (2006)](https://doi.org/10.1239/jap/1143936262):
  successful-mutant timing in branching processes is adjacent prior theory.
* [Gerold's 2019 dissertation](https://dash.harvard.edu/bitstreams/63433f5a-e8a9-43d8-9ee5-2465232d34bc/download),
  *Evolution and Immunity in Cancer and HIV*, includes antigen gain/loss
  extinction recurrences. Search-indexed relevant material was inspected;
  direct retrieval was rate-limited, so a complete dissertation review is
  **not** claimed.
* [Bozic et al. (2013)](https://doi.org/10.7554/eLife.00747) already distinguish
  common cross-resistance mutations from separate resistance paths. A common
  escape route limiting multiple targets is not a new principle.

**Novelty assessment.** Branching fixed points, the first-moment transform,
positive resolvents, convexity of moment sequences, and rare-rescue limits are
classical ingredients. The proposed specific contribution is the fixed,
positive-`h` **nonlinear escape asymptotic with its sharp coefficient**, the
finite-load bounds, and the heterogeneous-founder harmonic limit in Theorem
6. We found no exact antecedent in this targeted primary-literature search.
That is a candidate model-specific contribution, not a guarantee of
publication-level novelty. An equivalent killed-branching or resistance result
under different terminology remains possible. The rest is a rigorous synthesis
and a reproducible verification framework, not claimed new general theory.

## 1. Model, assumptions and outcomes

Every theorem uses the following assumptions unless it explicitly adds or
replaces one. Rates are fixed constants; all times share one arbitrary unit.

**A1 — Finite initial population, unbounded future population.** At time zero
there are `N < infinity` antigen-presenting cells with finite antigen sets
`A_1,...,A_N`, and optionally `m < infinity` presentation-null cells. Set
`k_i = |A_i|`. There is no population or antigen-count cap in the model.

**A2 — Division and death.** A presenting cell of load `k` divides into two
cells of the same type at rate `b > 0` and dies at rate `d + kappa*k`, where
`d > 0` and `kappa > 0`. Both daughters inherit all parental antigen IDs.
Equivalently, a birth adds one copy of the cell. Immune killing is linear in
the number of functional, equally weighted presented antigens.

**A3 — Mutational acquisition and individual loss.** A presenting cell gains
a new, unique antigen ID at rate `u > 0`. Each currently carried antigen is
lost independently at rate `ell > 0`; the total loss transition rate is
`ell*k`. A loss affects that cell and its descendants, not surviving sister
lineages. Lost IDs do not reappear, although new IDs may arise. Thus the count
process has transitions `k -> k+1` at `u` and `k -> k-1` at `ell*k`.
Acquisition is a **per-cell-time** mutation clock, not a per-division model.

**A4 — Common presentation loss.** At rate `h >= 0`, a presenting cell becomes
type `H` and permanently hides **all antigens included in this model**.
It then has birth rate `b_H > 0` and death rate `d_H > 0`, with no return to
presentation. Its descendants inherit this state. Mutation may change hidden
genetic antigens without affecting these rates. Write

\[
 q_H=\min\{1,d_H/b_H\},\qquad \sigma=1-q_H.
\]

For completeness, first-event conditioning on an `H` cell gives
`b_H*q^2-(b_H+d_H)*q+d_H=0`, with roots `1` and `d_H/b_H`.
Finite-event-tree extinction probabilities increase from zero to the smallest
fixed point in `[0,1]`, selecting the displayed `q_H`, including the critical
equality case. This is the standard birth–death extinction calculation.

This is an idealization of loss of the allele/haplotype presenting all modeled
targets, or complete modeled MHC-I/B2M presentation loss. **It is not an
assertion that ordinary HLA LOH hides all tumour antigens.** Residual HLA
presentation, NK killing and other immune pathways are absent.

**A5 — Independent branching.** All clocks and descendants are independent
conditional on their states. There is no shared T-cell budget, spatial
competition, antigen-frequency-dependent recruitment, density dependence,
random common immune environment or feedback from `H` cells.

**A6 — Clonal history and treatment interpretation.** A clonal antigen is an
ID present in every cell at the observation/treatment time, not an
unlosable antigen. Subclonal antigens occupy proper subsets. The same process
with a lower pretreatment `kappa` can generate a growing labelled tree; at a
finite observation or size stopping time, condition on its observed cells and
change `kappa`. The strong Markov/branching property gives the formulas below
conditional on that population. We do not claim its pretreatment distribution
is fitted to any cohort. The simulator includes such labelled growth and
reports extinction rather than silently restarting unsuccessful tumours.

**A7 — Distinct endpoints.** *Source extinction* means eventual absence of
presenting cells. *Total eradication* means the entire population reaches zero.
*Successful H escape* means some `H` family persists forever. One antigen-loss
event is not automatically escape: a cell can reacquire antigens or die.
Loss-enabled persistence concerns the full presenting branching process.
These are infinite-time mathematical events, not RECIST response or durable
clinical benefit, and not mutually exclusive first-route labels.

**A8 — Fixed-parameter limits.** In Theorem 6 all rates are fixed while loads
and founder numbers vary. In Theorem 4 only `h` and the initial population may
vary. Neither theorem is silently uniform in vanishing `kappa`, diverging
`ell`, saturation of killing, or changing antigen quality.

The code additionally accepts several zero-rate boundaries for audit purposes;
their interpretation is stated in Section 8. A1–A8, not that permissive input
validation, define the main theorems.

Define

\[
 a=\ell+\kappa,\quad c=\frac\ell a\in(0,1),\quad
 C=\frac{u\kappa}{a^2},\quad
 r=b-d-\frac{u\kappa}{a},\quad \rho=r-h,
 \quad g(t)=c+(1-c)e^{-at}.
\]

Let `Z(t)` be the number of source cells, `M_k(t)` its expectation starting
from one load-`k` cell **with h=0**, and

\[
 R_k(z)=\int_0^\infty e^{-zt}M_k(t)\,dt\quad(z>r).
\]

In particular `T_k=R_k(0)` is used only when `r<0`. For a nonnegative function
`f` of load, `R f` denotes the baseline expected integrated sum of `f` across
all source cells; `R 1=T`. No clinical unit calibration is implied.

## 2. Theorem 1 — Exact source mean and the correct extinction threshold

**Assumptions:** A1–A8, with one presenting founder of any finite load `k`;
no initial `H` cell when discussing the phase table.

For every `t >= 0`,

\[
 \boxed{M_k(t)=e^{rt+C(1-e^{-at})}g(t)^k},\qquad
 \mathbb E_k Z(t)=e^{-ht}M_k(t).                 \tag{1}
\]

The source becomes extinct almost surely **if and only if** `rho <= 0`.
If `rho>0`, its survival probability is positive for every finite `k`.
For any finite, nonempty initial source population this is the same threshold.
Writing `v=b-d-h`, the threshold can also be expressed as follows:

| Parameters | Source extinction |
| --- | --- |
| `v <= 0` | Almost sure for every positive `kappa`. |
| `0 < v < u` | Almost sure exactly when `kappa >= v*ell/(u-v)`. Equality is included. |
| `v >= u` | Positive source survival probability for every finite `kappa`. |

For total eradication from source-only initial conditions:

| Case | Almost-sure total eradication |
| --- | --- |
| `h=0` | Exactly when `r<=0`. |
| `h>0`, `b_H<=d_H` | Exactly when `rho<=0`. |
| `h>0`, `b_H>d_H` | Never: total survival has positive probability from every finite `k`. |

This last line is why a source threshold must not be advertised as a cure
threshold. Initial clonal burden changes probabilities and times even when it
does not change the infinite-time source phase boundary.

### Proof

**Nonexplosion.** Cell counts are bounded above on finite intervals by a Yule
process with per-cell birth rate `max(b,b_H)`. Its finite-time population and
number of births are finite almost surely. Gains have rate `u` per source cell,
so there are finitely many gains on each finite interval. Starting antigen
sets are finite; division copies only a finite set, and losses remove IDs.
Consequently the state-dependent loss/death rates do not introduce an
explosion. The following path and branching arguments are well defined.

An antigen has an independent lifetime `L~Exp(ell)` along a single lineage.
The probability weight for avoiding its immune killing up to `t` is

\[
 \mathbb E e^{-\kappa\min(L,t)}
 =\int_0^t\ell e^{-(\ell+\kappa)s}\,ds+e^{-(\ell+\kappa)t}
 =g(t).
\]

The `k` initial antigens multiply this weight to `g(t)^k`. Independent Poisson
gains of rate `u` contribute the product functional

\[
 \exp\!\left[u\int_0^t(g(s)-1)ds\right]
 =\exp\!\left[-\frac{u\kappa}{a}t+C(1-e^{-at})\right].
\]

Birth contributes `exp(bt)` to the first moment, baseline death contributes
`exp(-dt)`, and source removal by `H` conversion contributes `exp(-ht)`.
This follows either by conditioning on the first branching event or by the
linear first-moment equation: lineage generator `u` upward, `ell*k` downward,
with potential `b-d-kappa*k-h`. Multiplication gives (1), including `t=0`.
In particular

\[
 c^k e^{\rho t}\le\mathbb E_k Z(t)\le e^C e^{\rho t}.       \tag{2}
\]

If `rho<0`, `P(Z(t)>0)<=E Z(t)->0`, proving extinction. At `rho=0` the means
are uniformly bounded. A cell has probability at least
`epsilon=exp(-b)*(1-exp(-d))>0` of baseline death within one time unit without
division; gain/loss events cannot prevent this event. Conditional on at most
`L` cells at an integer time, the probability all source lineages disappear
within the next unit is at least `epsilon^L`. Repeated visits to that finite
population set on a surviving path have probability zero, by conditioning
successively on the visits. Hence survival implies `Z(n)->infinity`. Fatou's
lemma and the bounded means rule out positive probability of that event.

For `rho>0`, let `m_kj(s)` be the source mean transition kernel, including `h`.
Two distinct cells at time `t` have a unique most recent division at which
their ancestral paths split. Summing over that event gives the factorial
second moment

\[
 \mathbb E_k[Z(t)(Z(t)-1)]
 =2b\int_0^t\sum_j m_{kj}(s)
       [\mathbb E_j Z(t-s)]^2\,ds
 \le\frac{2b e^{3C}}\rho e^{2\rho t}.
\]

The factor two counts the two orders of daughter lineages; all terms are
nonnegative, so summation and integration require no unproved integrability
assumption. Applying Cauchy–Schwarz to `Z(t)*1{Z(t)>0}` and using (2) yields

\[
 \Pr_k(Z(t)>0)\ge
 \frac{c^{2k}}{e^C+2b e^{3C}/\rho}>0.
\]

The events `Z(t)>0` decrease with time because zero is absorbing. Their limit
therefore has positive probability. Solving `rho<=0` gives the displayed
`kappa` threshold, with the stated equality and `v>=u` cases.

If a source family dies, it generates finitely many `H` founders before its
finite extinction time. A finite collection of critical or subcritical
birth–death families dies almost surely. If `h>0` and `sigma>0`, conversion
as the founder's first event followed by `H` establishment already has positive
probability `h*sigma/(b+d+kappa*k+u+ell*k+h)`. This proves the total phase table.
Independent finite founder products give the population claims. **QED.**

## 3. Theorem 2 — Exact probabilities with computable rational enclosures

**Assumptions:** A1–A8. Let `q_k` be eventual total extinction from one
presenting cell. Set

\[
 D_k=b+d+\kappa k+u+\ell k+h,\qquad
 (F x)_k=\frac{d+\kappa k+b x_k^2+u x_{k+1}
                     +\ell k x_{k-1}+h q_H}{D_k},             \tag{3}
\]

where the `k-1` term is omitted at `k=0`. Then `q` is the componentwise
smallest fixed point of `F` in `[0,1]` and

\[
 Q(\mathbf k,m)=q_H^m\prod_{i=1}^N q_{k_i}.                  \tag{4}
\]

Starting `l^(0)=0`, `v^(0)=1`, the exact iterates satisfy
`l^(n)<=q<=v^(n)`. The lower iterates always increase to `q`. The upper
iterates decrease to `q` whenever **either** `rho<=0` **or** `h*sigma>0`.
For rational rates, this gives a terminating, arbitrarily accurate rational
enclosure algorithm throughout the main escape regime, including a supercritical
source with successful `H` escape. When `h*sigma=0` and `rho>0`, only the
least-fixed-point limit is asserted; the upper sequence can stay at one.
The implementation explicitly rejects a promise of a shrinking enclosure
in that exceptional regime.

For initial load `k`, depth `n` uses only types `0,...,k+n`: an event raises
the load by at most one. This is a finite computation on the original infinite
model, **not an antigen cap**. If `rho<0`, an explicit rate certificate is

\[
 v_k^{(n)}-l_k^{(n)}\le\min\{1,W_k/(n+1)\},                 \tag{5}
\]

where

\[
 W_k=1+\frac\ell\kappa+
 \left[2b+u+\frac{\ell(b-d-h)}\kappa\right]R_k(h).          \tag{6}
\]

Rounding every lower/upper update outwards to a grid of denominator `2^p`
preserves the brackets; it adds at most
`2*(2^n-1)/2^p` to (5). The default `p=2n+24` makes this extra error vanish.
The solver stops only when the **achieved exact rational interval width**
meets its requested tolerance. A user-set computation limit raises an error
instead of reporting an uncertified approximation.

### Proof

Conditioning on the founder's first event gives (3): at division both daughters
must become extinct, while at conversion the new `H` family must become
extinct. Monotonicity of `F` gives the bounds and their monotone directions.
The lower iterate is the probability that all source leaves terminate within
`n` event levels and every encountered `H` family eventually dies. On total
extinction the source genealogy is finite, so these events increase to total
extinction. This proves both the limit and minimality among fixed points.

The upper iterate ignores the fate of unresolved source leaves. Its limit is
the probability that no successful `H` family ever appears. If the source
dies almost surely, this is total extinction. If `h*sigma>0`, a surviving
source has `Z(t)>=1` for all time and hence infinite source cell-time.
The counting process of conversions has rate `h*Z(t)`; its usual independent
exponential-clock construction implies infinitely many conversions whenever
that integrated rate diverges. `H` establishment outcomes are independent of
the source tree, since the source is removed regardless of the `H` outcome.
Infinitely many such outcomes with success probability `sigma>0` cannot all
fail. Therefore source survival with no successful `H` has probability zero.
This also shows that, when `h*sigma>0`, `1-q_k` is the probability of successful
`H` escape, even if source survival is itself possible.

For (6), put `L(t)=E sum_cells k`. The mean balance equation is
`(E Z)'=(b-d-h)E Z-kappa*L`. When `rho<0`, integration and (1) give

\[
 \int_0^\infty L(t)dt=
 \frac{1+(b-d-h)R_k(h)}\kappa.
\]

Integrating the total event intensity
`(b+d+u+h)Z+(kappa+ell)sum k` gives (6). An unresolved path after `n`
event levels requires at least `n+1` events in the eventual finite source
tree. Markov's inequality proves (5). The potentially negative coefficient
of `R_k(h)` in (6) is handled with signed interval arithmetic in the code.

On `[0,1]`, each row of `F` has total absolute derivative at most
`(2b+u+ell*k)/D_k<=2`. One outward update costs at most `2^-p` per endpoint;
after `n` updates the endpoint error is at most `(2^n-1)/2^p`. This proves
the rounding bound and convergence of rounded brackets. Equation (4) follows
from conditional independence of the initial cells. **QED.**

## 4. Theorem 3 — Certified occupation formula and complete monotonicity

**Assumptions:** A1–A8 and a fixed real `z>r`. For rational rates and rational
`z`, the following positive series has an exact rational remainder bound.
Write `beta=(z-r)/a>0`. Then

\[
 \boxed{R_k(z)=\frac1a\sum_{j=0}^k {k\choose j}c^{k-j}(1-c)^j
 \sum_{n=0}^\infty
 \frac{C^n}{\prod_{m=0}^n(\beta+j+m)}}.                       \tag{7}
\]

For `v=beta+j`, the inner term `t_n=C^n/prod_(m=0)^n(v+m)` satisfies
`t_(n+1)/t_n=C/(v+n+1)`. If terms through `n=L-1` are retained and
`C/(v+L+1)<1`, the omitted tail lies in

\[
 \left[0,\;\frac{t_L}{1-C/(v+L+1)}\right].                 \tag{8}
\]

For the forward difference `Delta f_k=f_(k+1)-f_k`, every integer `j>=0`
satisfies

\[
 (-1)^j\Delta^j R_k(z)>0,\qquad
 R_k(z)R_{k+2}(z)>R_{k+1}(z)^2.                             \tag{9}
\]

Thus the occupation susceptibility is strictly decreasing, discretely convex,
and strictly log-convex. This is a model-specific application of standard
positive moment-sequence properties, not a new abstract convexity theorem.

### Proof

Substitute `x=exp(-at)` in the definition of `R` and use (1):

\[
 R_k(z)=\frac1a\int_0^1 x^{\beta-1}e^{C(1-x)}
                         [c+(1-c)x]^k dx.
\]

Expand the finite binomial and the nonnegative exponential series. Tonelli
permits exchange of the sum and integral. Repeated integration by parts gives
`integral_0^1 x^(v-1)(1-x)^n dx = n!/prod_(m=0)^n(v+m)` for `v>0`, proving
(7). Successive ratios decrease, so a geometric tail starting at `t_L` gives
(8), including the zero-tail boundary `C=0` used in tests.

Write `w_z(t)=exp((r-z)t+C(1-exp(-at)))>0`. Then

\[
 (-1)^j\Delta^j R_k(z)=
 \int_0^\infty w_z(t)g(t)^k[1-g(t)]^j dt>0.
\]

The integral is finite and `0<g(t)<1` for `t>0`. Cauchy–Schwarz applied to
`g^(k/2)` and `g^((k+2)/2)` in the measure `w_z(t)dt` gives log-convexity.
Equality would require the nonconstant function `g(t)` to be constant almost
everywhere under a positive density, which is impossible. **QED.**

## 5. Theorem 4 — Uniform rare-H error and the large-population transition

**Assumptions:** A1–A8, `r<0`, `sigma>0`, no initial `H`. All rates except `h`
are fixed. Let `S_k` be total source cell-time from one load-`k` founder in
the baseline process with `h=0`, and let `p_k(h)=1-q_k` be successful `H`
escape probability. Put `B=T_0+b*T_0^2`. For every finite `k` and every `h>=0`,

\[
 \mathbb E S_k=T_k,\qquad
 \mathbb E S_k^2\le2B T_k,                                  \tag{10}
\]

\[
 \boxed{0\le\sigma hT_k-p_k(h)\le\sigma h^2 B T_k}.         \tag{11}
\]

This is an absolute error bound, not an unqualified first-order heuristic.
For any finite load vector, set

\[
 A_1=\sum_i T_{k_i},\quad A_2=\sum_i T_{k_i}^2,\quad
 \Lambda=\sigma hA_1,\quad x=\sigma hT_0.
\]

If `x<1`, its total extinction probability obeys

\[
 \boxed{-\Lambda-\frac{\sigma^2h^2 A_2}{2(1-x)}
 \ \le\log Q\le -\Lambda+\sigma h^2 B A_1}.                 \tag{12}
\]

For **any sequence** of finite initial populations and `h_n->0`, with no
restriction on their individual loads, if `Lambda_n->L` in `[0,infinity]`,

\[
 Q_n\longrightarrow e^{-L},\qquad e^{-\infty}:=0.           \tag{13}
\]

Thus vanishing, finite and diverging total exposure give limits one, a
nontrivial probability, and zero. For a fixed empirical average susceptibility
`N^-1 sum T_(k_i)->tau>0`, this is the usual `N*h` transition with coefficient
`sigma*tau`. The general Poisson-rescue principle is prior theory; the
explicit susceptibility and uniform error here make it applicable and
checkable for this antigen model.

### Proof

Work first on the entire baseline genealogy, which is finite almost surely
because `r<0`. Tonelli gives `E S_k=T_k`. To establish the second moment
without assuming it is finite, partition ordered pairs of cell-time points
according to their genealogy. If one point is ancestral to the other, summing
over the earlier point contributes `2 R T`. Otherwise their paths split at a
unique most recent birth; the two daughter orderings, birth intensity `b`, and
independent expected future areas contribute `2b R(T^2)`. The diagonal set of
identical time points has zero two-dimensional measure. All contributions
are nonnegative, so Tonelli gives the extended-real identity

\[
 \mathbb E S_k^2=\{R(2T+2bT^2)\}_k.
\]

By Theorem 3, `T_j<=T_0` for every `j`; hence this expression is at most
`2(T_0+b*T_0^2)(R1)_k=2B T_k<infinity`. This proves (10) without circular
differentiation of an unproved finite second moment.

Place independent Poisson marks of rate `h` on all baseline cell-time, each
with an independent `H` descendant family. The actual source genealogy is
obtained by pruning the baseline descendants after every realized conversion.
The **earliest** potential mark is always realized. Conditional on the
baseline tree, the probability that at least one potential mark exists is
`1-exp(-h*S_k)`; its earliest mark establishes with probability `sigma`.
Conversely every actual successful mark is among the potential successful
marks, whose expected number is `sigma*h*T_k`. Consequently

\[
 \sigma\mathbb E(1-e^{-hS_k})\le p_k(h)\le\sigma hT_k.
\]

It would be incorrect to count *any* successful potential mark as a lower
bound, because an earlier failed conversion may prune its lineage. Using
`0<=y-(1-exp(-y))<=y^2/2` for `y>=0`, followed by (10), proves (11).

For `0<=p<=x<1`, the power series of the logarithm yields

\[
 -p-\frac{p^2}{2(1-x)}\le\log(1-p)\le-p.
\]

Sum this inequality over independent founders and use
`0<=sigma*h*T_(k_i)-p_(k_i)<=sigma*h^2*B*T_(k_i)` to get (12). Since
`A_2<=T_0*A_1`, its lower-error magnitude is at most
`Lambda*x/[2(1-x)]`, and its upper-error magnitude is `Lambda*h*B`.
Both are `o(Lambda)` uniformly over all load vectors as `h->0`. This proves
(13), including `L=0` and `L=infinity`. The case of an empty initial population
or `h=0` has `Q=1` directly. **QED.**

## 6. Theorem 5 — When one clonal antigen exceeds many subclonal targets

**Assumptions:** A1–A8, `N>=1`, and no initial `H` cells. Part (i) holds for all allowed
rates. Parts (ii)–(iv)
add `r<0`, `sigma>0`, compare the same finite number `N` of presenting cells,
and use the rare-`h` error certificate of Theorem 4. All compared antigens
have the same quality and the same individual loss rate.

**(i) Exact identity invariance.** Two initial tumours with the same multiset
of per-cell loads have the same total-extinction and successful-`H`
probabilities. One clonal antigen in every cell and `m` mutually exclusive
private antigens partitioning the cells into nonempty groups therefore
**tie exactly**, however large `m<=N` is. When `2<=m<=N`, their observed clonal
fractions are one and zero respectively. Initial clonal fraction alone does not identify
response even in this model.

**(ii) Coverage theorem at a fixed incidence budget.** Suppose tumour A has
one clonal antigen and no other antigen, so all loads are one. Tumour B has
any number `m` of subclonal antigen IDs with cell fractions `f_1,...,f_m` and

\[
 \sum_{j=1}^m f_j=\frac1N\sum_{i=1}^N k_i\le1.             \tag{14}
\]

Then `sum_i T_(k_i)>=N*T_1`, with strict inequality unless every `k_i=1`.
If strict, A has strictly greater total-eradication probability for all
sufficiently small positive `h`. This is a precise sense in which one
clonal antigen can be worth more than arbitrarily many subclonal antigens:
the condition fixes **antigen incidences per cell**, not just a count of
distinct mutation IDs.

More generally, at any fixed integer total incidence `K=sum k_i`, the
susceptibility sum is minimized when all loads differ by at most one.
If `K/N=s+theta`, `s=floor(K/N)`, `0<=theta<1`, the minimum is
`N*((1-theta)*T_s+theta*T_(s+1))`.

**(iii) An uncovered-fraction certificate.** If a fraction `z` of B's cells
has zero presented antigens and

\[
 z>\frac{T_1}{T_0},                                       \tag{15}
\]

then A has a strictly lower leading escape exposure than B, regardless of
how many antigens the remaining cells carry. This comparison does not require
(14). Equation (7) makes the right side explicitly computable.

**(iv) Finite-`h` rather than asymptotic-only comparison.** For profiles A,B
write `A_j=sum_(i in j) T_(k_i)` and `C_A=sum_(i in A) T_(k_i)^2`.
If `x=sigma*h*T_0<1` and

\[
 \sigma h(A_B-A_A)>
 \sigma h^2 B A_B+
 \frac{\sigma^2h^2 C_A}{2(1-x)},                            \tag{16}
\]

then `Q_A>Q_B` at that particular positive `h`. Every quantity in this
strict sign certificate has a rational enclosure from Theorem 3. No claim
that exact finite-`h` escape risk is globally convex is made.

### Proof

The labelled process maps to the count process because every transition rate
and destination count depends only on `k`, not the identities or the number
of other cells carrying each ID. Couple two populations by matching cells of
the same load, all future birth/gain/loss/death/`H` clocks, and descendant
counts. Their count paths and outcomes agree. This proves (i), including the
fact that a common past ancestor creates no further dependence under A5.

Strict discrete convexity from (9) says successive differences
`T_(k+1)-T_k` strictly increase. If two integer loads satisfy `a>=b+2`, then

\[
 T_a+T_b>T_{a-1}+T_{b+1}.
\]

Moving one incidence from the larger to the smaller load decreases the
susceptibility sum. Repeating finitely many times ends with loads differing
by at most one, proving the stated minimum. If the average is below one,
that minimum is `N*((1-mean)*T_0+mean*T_1)>N*T_1`. At average one equality
requires all loads one. This proves (ii)'s susceptibility assertion.

The zero-load cells alone contribute `z*N*T_0`, so (15) proves the strict
comparison in (iii). For (iv), apply the lower logarithmic bound (12) to A
and the upper bound to B. Inequality (16) makes the former strictly greater
than the latter. If `A_B>A_A` for fixed finite profiles, its left side has
positive order `h` while its right side has order `h^2`, proving the
sufficiently-small-`h` assertions. **QED.**

**Scalable example, not the theorem's scope.** For every even `N`, give half
of B's cells two private antigen IDs each and the other half none. B has
`N` distinct subclonal antigens and the same average load as one clonal
antigen, yet leading susceptibility `N*(T_0+T_2)/2>N*T_1`. Such a pattern can
arise on branches of a labelled tree. The result holds for every even `N`,
and the limit/error results control growing `N`; it is not a two-cell
counterexample. For the synthetic rates in Section 9, `N=20,h=1/100`
already satisfies the finite-`h` certificate (16).

## 7. Theorem 6 — Fixed-H loss produces a reciprocal-load escape law

**Assumptions:** A1–A8; all rates fixed, `h>0`, `sigma>0`, `rho=r-h<0`,
and no initial `H` cells. The baseline `r` may be positive: the requirement
is suppression of the **actual source after conversion**, not necessarily
of the `h=0` comparison process. Define `p_k=1-q_k`.

For every integer `k>=0`,

\[
 \boxed{\sigma h R_k(h+2b)\le p_k\le
              \min\{1,\sigma h R_k(h)\}}.                 \tag{17}
\]

For each fixed `z>r`,

\[
 kR_k(z)\longrightarrow\frac1\kappa;
 \qquad\boxed{kp_k\longrightarrow\frac{\sigma h}{\kappa}}.\tag{18}
\]

Consider any sequence of finite tumour populations with `N_n>=1` and positive initial
loads `k_(n,1),...,k_(n,N_n)` such that `min_i k_(n,i)->infinity`. Set
`H_n=sum_i 1/k_(n,i)`. If `H_n->L` in `[0,infinity]`, then

\[
 \boxed{Q_n\longrightarrow
       \exp\!\left[-\frac{\sigma h}{\kappa}L\right]}.       \tag{19}
\]

For homogeneous clonal burden `k_n`, this is the sharp transition
`N_n/k_n -> 0, L in (0,infinity), infinity`, with extinction limits
`1, exp(-sigma*h*L/kappa), 0`. The same theorem permits heterogeneous initial
trees; large loads enter through their reciprocal sum. A fixed number `m`
of initial `H` founders multiplies these limits by `q_H^m`.

These are fixed-`h`, large-population limits with full antigen gain/loss
dynamics, not a small mutation-rate substitution into a two-type formula.
At finite loads, (7), (17), and (4) give directly computable error brackets.

### Proof

For the lower bound, follow only the original cell until its first division,
death or conversion, leaving all its gain/loss clocks active. A conversion
at time `t` that precedes division and death and subsequently establishes
has density

\[
 \sigma h e^{-(b+d+h)t}
       \mathbb E_k\exp\!\left[-\kappa\int_0^t K(s)ds\right],
\]

where `K` is the immigration–death antigen-count path. Its integral is
`sigma*h*R_k(h+2b)` by (1). The shift is **`2b`**, because (1) already has a
positive first-moment birth factor `exp(bt)`, while this lower event requires
no birth and therefore has factor `exp(-bt)`. This successful route is a
subset of successful escape in the full branching process.

For the upper bound, count all conversions in the actual, pruned source
process. Their expected number is `h*R_k(h)`; each has independent
establishment probability `sigma`. The probability of at least one successful
family is no greater than their expected number. Source suppression makes
this integral finite and all source families extinct, so this is exactly
the complement of total extinction. This proves (17).

To prove (18), write
`R_k(z)=integral_0^infinity f(t)*g(t)^k dt`, with
`f(t)=exp((r-z)t+C(1-exp(-at)))`. Then `f(0)=1`,
`-log(g(t))/t->kappa` as `t->0`, and `integral f(t)dt<infinity` since `z>r`.
For any `epsilon in (0,min(1,kappa))`, choose a sufficiently small `delta>0`
so that for `0<=t<=delta`,

\[
 1-\epsilon\le f(t)\le1+\epsilon,
 \quad e^{-(\kappa+\epsilon)t}\le g(t)
                         \le e^{-(\kappa-\epsilon)t}.
\]

Integrating these elementary exponential bounds gives the small-interval
liminf at least `(1-epsilon)/(kappa+epsilon)` and limsup at most
`(1+epsilon)/(kappa-epsilon)` for `kR_k(z)`. On the remaining interval,

\[
 0\le k\int_\delta^\infty f(t)g(t)^k dt
 \le k g(\delta)^k\int_\delta^\infty f(t)dt\longrightarrow0,
\]

since `g(delta)<1`. Letting `epsilon->0` proves the resolvent limit.
Both values of `z` in (17) satisfy `z>r`, so squeezing proves the nonlinear
escape limit with the sharp coefficient `sigma*h/kappa`.

Let `D=sigma*h/kappa>0`. Since `kp_k->D`, convergence is uniform over the
integer tail `k>=min_i k_(n,i)`. Therefore
`sum_i p_(k_(n,i))=(D+o(1))*H_n`, and `max_i p_(k_(n,i))->0`.
If `H_n` is bounded, `sum_i p_i^2 <= max_i p_i*sum_i p_i ->0`.
The logarithmic series gives `log Q_n=-sum_i p_i+o(1)`, proving (19) for
finite `L`. If `H_n->infinity`, eventually `p_k>=D/(2k)` for every initial
load, so `Q_n<=exp(-D*H_n/2)->0`. Empty populations have probability one
directly. Initial `H` cells contribute the independent factor in (4).
**QED.**

### A quantitative heterogeneity penalty

For large homogeneous load `k` and the same `N`, compare a population with
half its cells at `k/2` and half at `3k/2` (take even integer `k,N`). The
arithmetic average load stays `k`, while the reciprocal-load sum rises from
`N/k` to `4N/(3k)`. Along the limits of Theorem 6,

\[
 Q_{\rm heterogeneous}\longrightarrow
       (Q_{\rm homogeneous,limit})^{4/3}.                  \tag{20}
\]

The power `4/3` is a distributional prediction, not a property of the
fraction of distinct antigen IDs called clonal. A persistent low-load
subpopulation violates the uniform high-load assumption and must retain its
exact product factor from Theorem 2.

## 8. Boundary cases and limits of biological interpretation

The assumptions matter at the boundaries:

* **Zero presentation loss or nonestablishing `H`.** If `h=0` or `sigma=0`,
  successful `H` escape is absent. With a suppressed source, total extinction
  is one. A supercritical source can persist even though successful `H` escape
  is zero; the upper event-tree sequence then need not identify total
  extinction. The code explicitly distinguishes this case.
* **Zero individual antigen loss.** If `ell=0`, (1) still has a meaningful
  limit, but its long-time exponent from load `k` is
  `b-d-u-kappa*k-h`, not simply `r-h`. Types can only increase, so the
  initial antigen floor becomes relevant. The same mean/critical argument
  proves extinction at nonpositive exponent and survival at positive exponent
  (the same-type family alone is then supercritical). The code skips
  zero-weight terms in (7) and requires `z>r-kappa*k` in this boundary.
  Thus taking `ell->0` and taking time to infinity need not commute.
* **Critical `H` families.** `b_H=d_H>0` gives `sigma=0` and almost-sure
  extinction, despite potentially long survival times. A simulation reaching
  its cap for this process is unresolved, not a surviving `H` family.
* **No antigens initially.** `k=0` is permitted in Theorems 1–5 and the
  finite-load bounds of Theorem 6. There is no loss transition out of zero,
  but acquisition remains active. The reciprocal-load population limit
  requires positive loads tending to infinity and cannot include zero-load
  cells by division by zero or by discarding them.
* **Existing escaped cells.** Formula (4) includes them explicitly. In
  particular, `m` pre-existing supercritical `H` cells cap eradication at
  `q_H^m`, even if the source is cleared instantly. Binary “any HLA LOH”
  status does not identify `m` or `q_H`.
* **Mutation timing and antigen quality.** Per-division acquisition, correlated
  deletion of several antigens, allele-specific presentation maps, unequal
  quality, changing birth fitness and changing immune pressure are different
  models. We do not assert (1), the phase boundary or (18) unchanged there.
* **Selection is not mutation induction.** Increasing `kappa` changes survival
  of lineages. Nothing here says immune pressure increases the mutation or
  HLA-loss rate. Escape probabilities are not measurements of those rates.

Two limitations are particularly consequential for biology. First, the model
cannot derive a response coefficient from a clonal-fraction summary alone:
it needs the per-cell load distribution and effective antigen quality, and
Theorem 5 supplies an exact counterexample to identifiability from labels.
Second, the linear, unlimited killing law can fail when effectors or
presentation saturate. If total death is bounded above by `d+K_max`, the
founder can convert before birth or death with probability at least

\[
 \frac{\sigma h}{b+d+K_{\max}+h}>0
\]

regardless of its antigen load, while gains/losses continue. This follows by
integrating `sigma*h*exp(-(b+d+K_max+h)t)` along the founder. Hence a
`1/k` decay is impossible in that bounded-killing model. Theorem 6 is a
falsifiable mechanistic scaling under A2, not a universal law of patients.

## 9. Numerical certificates and independent simulation

The shared CLI example has

\[
 (b,d,u,\ell,\kappa,h,b_H,d_H)
  =(1,1/2,1,1/4,1,1/10,1,1/4).
\]

These are arbitrary synthetic rates, not a calibration. Here
`r=-3/10`, `rho=-2/5`, `sigma=3/4`, and the Theorem 6 coefficient is exactly
`3/40=0.075`. The following decimals display rational brackets computed at
event depth 96; the printed decimal is not itself the exact endpoint.

| Initial load `k` | Total extinction `q_k` (display) | `k*(1-q_k)` (display) | Certified `q_k` interval width, at most |
| --- | --- | --- | --- |
| 0 | 0.7954057017506 | 0 | `6.35e-13` |
| 1 | 0.9103845941743 | 0.0896154058257 | `2.18e-13` |
| 4 | 0.9782751260851 | 0.0868994956598 | `6.06e-15` |
| 16 | 0.9951262190890 | 0.0779804945753 | `1.38e-21` |
| 64 | 0.9988162913607 | 0.0757573529147 | `9.49e-49` |

The approach toward `0.075` is a check of a general theorem, not its proof.
For `h=0` baseline occupation the same parameters give
`T_0=5.6590588081387...`, `T_1=1.8295294040694...`. Equation (15) therefore
has uncovered-fraction threshold `T_1/T_0=0.3232921703231...`. The test also
compares a homogeneous population to zero-load/two-load coverage patterns
with the exact finite-`h` certificate, not an assumed ordering of simulations.

An additional reproducible run of `sample_escape` used 12,000 trials per row,
the same rates, and `hla_threshold=24`; the bound on an `H` misclassification
is `4^-24`. All three runs had zero unresolved trials. Seeds are explicit:

| Starting load multiplicities | Seed | Exact extinction probability (display) | Simulated eradications / 12,000 |
| --- | --- | --- | --- |
| One cell at load 1 | 20261717 | 0.9103845941743 | 10,883 |
| 64 cells, each at load 64 | 20261718 | 0.9269995251636 | 11,093 |
| 32 cells at load 32 and 32 at load 96 | 20261719 | 0.9032315809330 | 10,794 |

The test suite independently repeats these types of comparisons with its own
declared fixed seeds and sampling tolerances. For the lower presentation-loss
rate `h=1/100`, depth-112 certificates give eradication in
`[0.7757015932,0.7757017448]` for 20 load-one cells and
`[0.6493724620,0.6493727462]` for ten zero-load plus ten load-two cells.
These are outward decimal enclosures of rational bounds and satisfy the
strict comparison in (16).

### What is checked independently

The dedicated test is
[`tests/test_research_immune_escape_neoantigen_clonality.py`](../../tests/test_research_immune_escape_neoantigen_clonality.py).
It uses only the Python standard library and deterministic seeds.

| Claim | Exact arithmetic check | Separate stochastic check |
| --- | --- | --- |
| Theorem 1 | Rational generator identities across a parameter grid; critical surface and zero-loss boundary. | Gillespie means below, at and above the source threshold. Finite-time samples do not prove infinite-time survival. |
| Theorem 2 | Dyadic event-tree intervals; known scalar quadratic; an independently assembled **forward population CTMC** with two overflow fates, solved by rational Gaussian elimination. | Source events and explicit `H` birth/death families for single and multiple founders. |
| Theorem 3 | Positive-series remainders; signed resolvent equations; four orders of finite differences; strict log-convexity. | Entire source-genealogy cell-time, including acquired and lost antigens. |
| Theorem 4 | Finite-`h` expansion errors and exponentiated product bounds for decreasing `h`. | Rare-`H` full-population outcomes checked against both exact and expansion bounds. |
| Theorem 5 | Finite-`h` sign certificate at equal incidence budget. | Coverage comparison; identical count paths for shared versus private labels; inherited labels in a growing 64-cell tumour. |
| Theorem 6 | Finite-`k` probability sandwich; scaled coefficient; homogeneous and heterogeneous population products against their limits. | Populations of 64 founders at load 64 or half at 32 and half at 96, with continued mutation/loss. |

`verifier.py` rejects binary floats as certificate inputs. Arithmetic uses
integers and `fractions.Fraction`; `decimal.Decimal` only formats certified
values. The simulator imports no verifier, uses no extinction recurrence,
and does not draw an `H` outcome from a known extinction probability. It
simulates source clocks and `H` birth/death jumps. It stops at the **first**
supercritical `H` family reaching `M` cells, with possible false-persistence
probability at most `(d_H/b_H)^M`. No multiple-family multiplier is needed:
conditional on that first crossing, whole-tumour extinction requires the
currently observed `M`-cell `H` family to die. Other caps return `unresolved`.
The tests include explicit censoring cases to guard against false eradication.

For bounded outcome checks, a Hoeffding/union-bound tolerance covers at most
32 Bernoulli comparisons with family error at most `1e-6`, in addition to
the stated `H` classification error. The unbounded mean checks are explicitly
six-standard-error diagnostics using sample variance, **not** rigorous
distribution-free confidence intervals. Exact certificates and proofs carry
the assertions; fixed-seed simulation is a different implementation check.

## 10. Falsifiable aggregate-cohort predictions and their empirical status

### A quantitative hypothesis, with an explicit observation bridge

For sufficiently large homogeneous functional clonal burden `k`, fixed
effective founder count `N`, fixed presentation-loss hazard and immune
strength, Theorem 6 gives the limiting curve

\[
 Q(k)\simeq e^{-A/k},\qquad A=\sigma hN/\kappa>0.           \tag{21}
\]

Applying this to **clinical response** is an additional, presently unvalidated
hypothesis. Suppose an independently calibrated endpoint bridge has
`s=P(response | model extinction)` and `f=P(response | no model extinction)`,
with `0<=f<s<=1` constant across the comparison strata. Then

\[
 p(k)=f+(s-f)Q(k),\qquad
 \log\!\left[-\log\frac{p(k)-f}{s-f}\right]
       \simeq\log A-\log k.                               \tag{22}
\]

There is also a measurement bridge: `k` counts functional, equally weighted
presented antigens, whereas cohorts often report predicted binders or inferred
clonal mutations. Applying this law to a published burden requires a
calibrated, comparable relationship between that measurement and functional
load across strata. Antigen quality, expression, HLA restriction and
clonality-calling errors can violate it. Neither bridge is inferred from the
branching mathematics.

**Prediction P1:** the transformed clonal-dose slope is `-1`; doubling
functional clonal burden gives

\[
 \frac{p(2k)-f}{s-f}\simeq
       \sqrt{\frac{p(k)-f}{s-f}}.                          \tag{23}
\]

**Prediction P2:** at the same arithmetic mean burden, the specific
half-`k/2`, half-`3k/2` distribution predicts the adjusted response probability
raised to the power `4/3`, as in (20). Testing P2 requires aggregate clone
fraction/load distributions; a distinct-antigen clonal fraction is insufficient.

For the stronger idealized bridge `f=0,s=1`, (21) predicts a particularly
simple quantitative odds relation:

\[
 O(k)=\frac1{e^{A/k}-1},\quad
 \frac{O(2k)}{O(k)}=e^{A/(2k)}+1>2.                         \tag{24}
\]

The equality follows by factoring `e^(A/k)-1`; the test checks that algebra
with exact rationals. It is **not** a statement that RECIST response equals
eradication. If this bridge or common `A` is wrong, (24) need not hold.
Likewise, (21) is a proved scaling limit, not an exact finite-`k` curve.
Finite-cohort claims need a prespecified approximation error derived from
(17), or must be labelled empirical hypotheses.

### How public aggregate data can falsify the hypothesis

Only published group counts and numerical dose bounds are needed; no patient
records are required. Prespecify a histology, treatment and response endpoint,
presentation eligibility, a calibration of `f,s`, and matching assumptions for
the quantities in `A`. For at least three disjoint dose bins collect published
`[L_j,U_j]`, denominator `n_j`, and responders `y_j`, with `0<L_j<=U_j<infinity`.
Do not replace the bounds by an invented or unreported representative dose.
The hypothesis implies every individual probability in bin `j` is between

\[
 p_j^-(A)=f+(s-f)e^{-A/L_j},\qquad
 p_j^+(A)=f+(s-f)e^{-A/U_j}.                                \tag{25}
\]

For any proposed `A`, reject if **any** bin has a lower-tail binomial
probability at `p_j^-(A)` or an upper-tail probability at `p_j^+(A)` below
`alpha/(2J)`. Reject the one-parameter curve if **no** `A>0` passes all bins.
This remains conservative when probabilities vary within a bin: coupling
Bernoulli trials to common uniforms stochastically bounds their sum between
the two endpoint binomials. It assumes independent outcomes across people.
Binomial tails can be evaluated as finite sums of
`comb(n,j)*p**j*(1-p)**(n-j)` with rational probability enclosures from the
verifier. Enlarge (25) first by the justified finite-load and bridge
uncertainty budgets; do not choose those budgets after seeing discrepancies.

For the uninflated curve with fixed `f,s`, equivalently invert the two binomial
tails into acceptance endpoints
`[p_Lj,p_Uj]`, check for intervals wholly outside `[f,s]`, and set
`q_Lj=(p_Lj-f)/(s-f)`, `q_Uj=(p_Uj-f)/(s-f)` clipped to `[0,1]` only after
that check. Every bin permits

\[
 A\in[-L_j\log q_{Uj},\;-U_j\log q_{Lj}],                  \tag{26}
\]

with the usual infinite endpoint at zero, and intersect with **finite `A>0`**.
If `q_Uj=0` or `q_Lj=1`, no such `A` is admissible; infinite or zero singleton
endpoints do not rescue the hypothesis. An empty intersection rejects the
quantitative hypothesis plus its bridge/common-parameter assumptions.
Equation (26) follows by requiring the intervals in (25) to intersect the
observed acceptance interval. Nonzero uncertainty budgets require inversion
of the correspondingly enlarged intervals. Fitting `f,s` freely on the same bins would
weaken or erase this falsifiability, so that is not the proposed protocol.

### What was actually checked against published aggregates

The [McGranahan 2016 main-text/figure-4 summary](https://pmc.ncbi.nlm.nih.gov/articles/PMC4984254/)
reports durable benefit in **12/13** tumours with at least 70 clonal
neoantigens and less than 5% subclonal neoantigens, versus **2/18** in the
complementary group. The crude cross-product odds ratio is exactly
`(12*16)/(1*2)=96`. This is a combined burden/heterogeneity criterion and an
observational association, not an isolated dose experiment or an eradication
endpoint. It is directionally consistent with coverage; it cannot identify
`A`, the reciprocal coefficient, or the shape in (22).

[Ravi et al. 2023, Nature Genetics](https://doi.org/10.1038/s41588-023-01355-5)
reports a clonal-neoantigen response OR of **5.4 (95% CI 2.7–11)** and no
detected pretreatment response association for HLA LOH, HLA allele number,
or B2M loss. The first result is consistent with the direction of coverage;
the latter is contrary evidence against an automatic blanket HLA-loss/response
claim. Its log-transformed predictor summary does not supply eligible dose
bins or the precise neoantigen transformation needed to test (24) directly.

**Empirical status:** the reviewed public aggregate summaries do not provide
verified numerical clonal-dose intervals, response counts, independent bridge
calibration and matching sufficient to execute (25) as a test of the new
coefficient. We did **not** retrieve patient-level data or fit a clinical model.
P1/P2 and the rejection rule are quantitative and testable using an eligible
public aggregate release; their direct cohort test remains **unperformed**.
The existing associations are consistency checks, not validation. A universal
prediction from clonal fraction or binary any-HLA-LOH alone is deliberately
rejected, because neither summary identifies this model's inputs.

## 11. Hostile-referee and falsification log

This log records actual challenges to the claims, including rejected stronger
claims. Separate derivation and critique were performed in independent agent
tasks and checked against distinct code paths. That is computational
cross-checking, not external scientific peer review.

| Attempt | Outcome and resulting restriction |
| --- | --- |
| Search for a first clonality explanation/threshold. | Rejected as prior art: Lakatos, Aguadé-Gorgorió–Solé, Qi–Vincent–Cao, Chen and Sholokhova already address the broad question. Pivoted to the precise fixed-`h` law. |
| Claim the explicit immigration–death transform as new. | Rejected after locating Karlin–Tavaré 1982. It is credited as a classical ingredient. |
| Treat a negative source exponent as guaranteed total eradication. | Rejected by a positive-probability first-event `H` conversion. The source and total phase tables are separated. |
| Extend one-clonal superiority to identical per-cell loads. | Exactly false by count-process coupling; retained as Theorem 5(i). Added the incidence/coverage conditions. |
| Treat initial clonal antigens as permanent. | Rejected: they are inherited but remain individually losable. The zero-loss case has a different initial-load-dependent threshold. |
| Transfer susceptibility convexity to exact escape risk at all parameters. | Falsified by the rationally certified second difference below. Only the proved rare-`h` comparison is claimed. |
| Use any successful baseline-tree mark as a lower bound. | Rejected because failed earlier conversions prune later marks. Theorem 4 uses the earliest mark only. |
| Use `h+b` in the first-cell lower resolvent. | Rejected: removing the first-moment birth factor and imposing no birth requires `h+2b`. |
| Drop the minimum-load condition in the harmonic law. | Rejected: fixed low-load founders retain exact, non-asymptotic factors. |
| Extend reciprocal decay to saturating immune killing. | Analytically false by a positive load-independent escape lower bound (Section 8). |
| Equate ordinary HLA LOH with total presentation loss or clinical response with extinction. | Rejected using the biological source definitions and the explicit observation bridge. Ravi's negative HLA association is retained as contrary evidence to a blanket extension. |
| Use rounded Decimal values as if inside arbitrarily narrow exact intervals. | An initial algebra test failed at decimal rounding precision. Replaced by an exact rational identity; certificate decisions never use formatted decimals. |
| Omit initial-H and singleton-partition boundary quantifiers. | Independent review caught both: Theorem 5 now excludes initial `H`, and its zero-clonal-fraction partition requires at least two groups. |
| Call arbitrary real inputs a rational algorithm. | Restricted the computational guarantee to rational inputs; mathematical formulas retain their stated real-parameter scope. |

The retained finite-`h` convexity counterexample is an infinite-type process,
not a truncated model. For

\[
 (b,d,u,\ell,\kappa,h,b_H,d_H)
 =(10,1/10,10,1/10,1,1/10,1,1/4),
\]

depth-160 event-tree intervals certify

\[
 p_0-2p_1+p_2\in
 [-0.004580404,\;-0.004580403].
\]

The test recomputes the negative sign at depth 176. This is outside the
`r<0` rare-`H` assumptions, so it does not contradict Theorems 3–5; it
disproves an unjustified broader claim. The main contribution remains the
general-parameter probability bounds and scaling limits, not this example.

The targeted novelty search used combinations of “neoantigen extinction
probability loss,” “antigen escape inverse branching probability,” “neoantigen
harmonic,” “antigen loss HLA redundancy mathematical,” “immigration death
linear killing Feynman Kac,” “multitype branching large type killing escape,”
and “evolutionary rescue mutational load reversion,” followed by relevant
primary references. Search date: 17 November 2026. It was not an exhaustive
MathSciNet/zbMATH/dissertation review. Publisher redirects and incomplete
indexing limit negative search conclusions. No unproved extension is promoted
to a theorem: clinical P1/P2 remain conditional hypotheses, and priority
remains explicitly unestablished.

## 12. Reproduce and audit

Run from the repository root with Python 3.10+; the existing CI uses 3.12.
The new verifier, simulator and test require **only the standard library**.

```sh
python -B research/immune-escape-neoantigen-clonality/verifier.py --depth 96 --series 32
python -B research/immune-escape-neoantigen-clonality/simulation.py --trials 20000 --seed 20261117
python -B -m unittest discover -s tests -p 'test_research_immune_escape_neoantigen_clonality.py' -v
```

The dedicated suite has 15 tests and completed in approximately **11–13
seconds** locally, comfortably below 60 seconds. Runtime varies by machine;
the workload is deterministic and bounded. Fixed seeds make regressions
reproducible; changing a seed is not a way to remove an adverse result.

The repository's existing validation gates were also run:

```sh
python scripts/validate.py
python -m unittest discover -s tests -v
node --test tests/*.test.*
```

Local results: structural validation passed; **270 Python tests** passed in
about 29 seconds; **26 JavaScript tests** passed. The existing pull-request
workflow discovers the new Python test automatically. No workflow or
dependency file is changed. CI reports implementation/integrity checks,
not scientific endorsement.

### File and API map

* [`verifier.py`](verifier.py): exact `Parameters`, `mean_bounds`,
  `resolvent_bounds`, `event_tree_bounds`, `solve_extinction`,
  `population_extinction_bounds`, `rare_h_escape_bounds`,
  `finite_h_escape_bounds`, and rational error certificates. Inputs may be
  integers, fractions, decimal literals or rational strings, not binary floats.
  `population_extinction_bounds` covers source founders; multiply its endpoints
  by exact `q_H**m` for specified pre-existing `H` founders.
* [`simulation.py`](simulation.py): independent `Parameters`,
  `sample_escape`, `source_lifetime_samples`, `sample_source_at_time`,
  and labelled `grow_labelled`. Rate names differ explicitly: `u/gain`,
  `ell/loss`, `h/hla`, `bH/b_hla`, `dH/d_hla`. Simulations use ordinary
  floating-point clocks, separately from exact certificate arithmetic.
* [The single new test file](../../tests/test_research_immune_escape_neoantigen_clonality.py):
  exact forward population construction, independent linear solver,
  theorem-specific assertions, falsification cases and fixed-seed comparisons.

Only these three research files and that one test file are added. No existing
research note, root README, catalog, result, script, policy, workflow or
requirement is modified. This note does not claim to solve immune escape in
patients; it provides a precise candidate mathematical refinement, complete
proofs and auditable limits on its interpretation.
