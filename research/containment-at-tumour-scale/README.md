# Containment at tumour scale: a finite-course, rare-mutation theorem

**Nonclinical computational mathematics. Synthetic parameters; no treatment
advice, patient data, experimental validation, or claim of clinical efficacy.**
Research status: **unreviewed**. Proofs below concern the explicitly defined
hybrid model. AI referee checks and passing CI are not external scientific
review. Literature audit: 17 November 2026.

## Result, and what this does not settle

There is no universal advantage of containment. In the model below, both
policies have the **same bounded Hill dose-response law**, mutation occurs
**per sensitive-cell division**, resistant founders undergo stochastic
establishment, and progression is **total burden reaching a fixed multiple
of the initial burden**. Both policies stop dosing after the same finite
course, obey a common cumulative-exposure constraint, and are followed through
subsequent regrowth. The endpoint is the **ordinary, uncensored expected
progression time of this hybrid**, not time to the first mutant.

For initial sensitive population `N*s0`, `R(0)=0`, mutation probability
`mu=lambda/N`, and course length `H_N=L*log(N)`, we prove

$$
\frac{\mathbb E T_C}{\log N}\longrightarrow c,\qquad
\frac{\mathbb E T_M}{\log N}\longrightarrow
d+(J-d)e^{-\lambda A}.                                      \tag{1}
$$

Here `C` maintains a sensitive reservoir by bounded feedback, `M` applies
maximum dose during the course, `d=1/g0` is the resistant expansion coefficient
after competitive release, `c=min(L,1/gC)` is the containment coefficient
(take `c=L` if `gC<=0`), and `J=L*(1+delta/rS)` is the sensitive rebound
coefficient after maximum dose. `A` is an explicit, rigorously computable
establishment-weighted division integral. All definitions and restrictions
appear below.

If competition is positive and the assumptions of Theorem 3 hold, `d<c<J`.
There is **exactly one** boundary in fixed positive mutation supply:

$$
\boxed{\lambda_*=
\frac{\log((J-d)/(c-d))}{A}.}                                \tag{2}
$$

Containment has the larger limiting mean coefficient precisely when
`lambda>lambda_*`; maximum dose has the larger coefficient when
`0<lambda<lambda_*`. A strict coefficient difference proves the same strict
ordering for all sufficiently large `N`. **At equality only the leading
coefficients tie.** It is not a finite-population equality theorem.
With no sensitive-to-resistant competition, or with nonpositive released
resistant growth, maximum dose wins throughout the corresponding regime.
The two reversals in the earlier two-cell note do not occur in this limit.

This is a narrow answer to the first follow-up in [the research index](../README.md),
not a resolution of treatment choice in real tumours. The comparison requires
a course and exposure budget growing as `log(N)`, a deterministic sensitive
bulk, perfect clone measurement, and a specified feedback policy. It neither
optimizes over all admissible controls nor proves the finite-size boundary at
`10^8`--`10^11` cells. Simulations at those sizes are reported separately from
the proven scaling limit. Finite-size corrections are substantial.

## 1. What was already known, and the precise pivot

The [merged mutation-supply note](../adaptive-therapy-mutation-supply/)
(PR #53) is an exact two-cell counterexample with a per-time conversion clock,
singular containment, and an exponentially capped endpoint. It explicitly
does not prove a division-linked, bounded-dose, tumour-scale result.

| Primary source inspected | Established result relevant here | Boundary of applicability |
| --- | --- | --- |
| Hansen, Woods & Read (2017), [DOI](https://doi.org/10.1371/journal.pbio.2001110) | The mutation-supply versus competitive-suppression trade-off, balance thresholds, and cases favouring either strategy. | Aggression is idealized as immediate sensitive removal; the analysis does not supply the common finite-dose hybrid below. |
| Hansen & Read (2020), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7428817/), [DOI](https://doi.org/10.1111/eva.12994) | Exact stochastic birth-death cure/progression calculations. Box 1 already ties mutation to sensitive divisions. Cure versus delay and the chosen objective matter. | Aggression immediately removes sensitive cells. Box 4 describes small-population numerical illustrations with inflated mutation probabilities; its analytical formulas are not restricted to two cells. |
| Viossat & Noble (2021), [DOI](https://doi.org/10.1038/s41559-021-01428-w), [author manuscript](https://openaccess.city.ac.uk/id/eprint/25716/) | General deterministic containment comparisons, including bounded-dose containment, and advantages without a fitness cost. | The main comparison omits ongoing mutation; supplementary numerical mutation investigations are not a rare-mutation scaling theorem. |
| Alvarez & Viossat (2024), [DOI](https://doi.org/10.1007/s00285-024-02062-3), [2022 author preprint](https://arxiv.org/html/2211.17041v1) | Mutation-compatible deterministic comparison under conditions on absolute resistant growth, including Norton-Simon treatment effects. | The accessible preprint assumes positive initial resistance and resistant growth nonincreasing in S. Appendix A already discusses `10^10` cells: large counts alone are not new. The subscription-limited final text was not represented as identical to the preprint. |
| Wang, Scott & Vladimirsky (2024), [full text/DOI](https://doi.org/10.1371/journal.pcbi.1012165) | Stochastic control with objectives depending on cumulative drug/time budgets and probabilities of reaching targets. | Environmental SDE noise; the sensitive/resistant example in section 2.3 explicitly neglects mutation. |

Two additional primary sources prevent broader novelty claims.
[Kuosmanen et al. (2021)](https://doi.org/10.1371/journal.pcbi.1009418)
already combine inhomogeneous mutation arrivals, branching establishment,
Hill pharmacodynamics, and less-aggressive optima; their establishment
approximation uses a constant survival probability.
[Gevertz et al. (2026)](https://doi.org/10.1038/s41540-025-00635-6)
already integrate pharmacokinetics and toxicity with adaptive therapy, without
the division-linked rare-mutation theorem considered here.

**Pivot:** the broad supply/suppression question is established. The question
not established in these audited sources is the exact comparison (1)--(2)
for a common bounded-dose, finite-course hybrid, with a proof connecting
stochastic founders to **macroscopic total-burden progression** and an honest
post-course endpoint. This targeted search cannot establish worldwide
priority. Main-text access to the Nature paper was intermittent; its indexed
primary text and author manuscript metadata were checked. Public full texts,
the explicitly versioned preprint, and publisher metadata are cited rather
than copied into this repository.

## 2. Complete model and assumptions

Every theorem lists its dependencies on the following assumptions. Rates and
parameters are synthetic and are not estimates for a particular cancer.

### A1. Populations, mutation, and the hybrid boundary

`N>=2` is a cell-number scale; `s=S/N` is a deterministic sensitive bulk and
initially `s=s0>0`, `R=0`. A sensitive division occurs at per-cell rate `bS`.
Its added daughter is resistant with probability `mu=lambda/N`; the parent
remains sensitive. Thus sensitive net production uses `(1-mu)*bS`, and mutant
arrivals have rate `mu*bS*N*s=lambda*bS*s`. This is not a conversion clock.
Assume `0<=mu<=1`.

Choose a fixed `0<a<1/2`, and `m_N=ceil(N^a)`. Before the first time `R=m_N`,
R is integer-valued with the following rates, where `r=R/N`:

| Event | Change in R | Rate |
| --- | --- | --- |
| Resistant division | `+1` | `bR*R` |
| Mutant daughter from a sensitive division | `+1` | `lambda*bS*s` |
| Resistant death | `-1` | `(dR+alpha*s+beta*r)*R` |

During this stage,

$$
\dot s=[r_S-\mu b_S-a_{SR}r-k(u)]s,\qquad r_S=b_S-d_S.       \tag{3}
$$

At `R=m_N`, replace R by its density and continue with (3) and

$$
\dot r=(g_0-\alpha s-\beta r)r+\frac{\lambda b_S s}{N},
\qquad g_0=b_R-d_R.                                        \tag{4}
$$

The switch occurs once, with no population jump. Mutation continues after
the switch and after dosing stops. This hybrid is the **model being proved**;
it is not asserted to be an exact finite-N two-type CTMC. The switch scale
goes to infinity while `m_N/N` goes to zero. The simulations use `a=1/4`,
so the initial tumours have `10^8`--`10^11` cells, not 100--563 cells: only
the resistant compartment switches at those smaller counts.

### A2. Demography and competition

`bS,bR>0`, `dS,dR>=0`, `alpha,beta,aSR>=0`, and `rS>0`.
`alpha` measures suppression of R by S, `aSR` reciprocal suppression, and
`beta` R self-competition. One measurable resistance-cost parametrization is
`bR=(1-c_res)*bS`, `0<=c_res<1`; the theorem uses the independently measurable
`bR,dR`. A cost in birth rate alone does not determine net fitness.

The sensitive bulk has no self-limitation in this local, pre-progression
model. Reciprocal competition and resistant self-limitation are retained.
No law beyond the progression boundary is needed. There is no backmutation,
phenotypic switching, spatial structure, initial resistance, immune response,
or direct drug effect on R.

### A3. One bounded dose law and a finite common course

$$
k(u)=k_\infty\frac{u^h}{K_u^h+u^h},\quad
0\le u\le u_{\max},\quad h,K_u,k_\infty,u_{\max}>0,
\quad k_{\max}=k(u_{\max}).                                \tag{5}
$$

This is a quasi-steady pharmacodynamic control model. Bounded dose changes
are allowed; no instantaneous cell removal is allowed. It does **not**
model a separate absorption/clearance compartment or validate a clinical
maximum tolerated dose. Here “maximum dose” means the specified mathematical
upper bound, not a fitted or approved regimen.

Fix `L>0`, `H_N=L*log(N)`, and accumulated exposure

$$
\dot z=u,\quad z(0)=0,\quad z(t)\le B_N^{tox}=u_{\max}H_N.
                                                               \tag{6}
$$

Both policies set `u=0` at and after `H_N`. The exposure constraint is finite
for every N and is binding for M if the full course is completed. It is a
linear cumulative-toxicity proxy, with no recovery. The course and budget
are part of the experimental design. **They grow as `log(N)` in the theorem**;
a fixed, N-independent lifetime budget is not covered. C may leave budget
unused; equal exposure or continued C dosing after M stops would be different
comparisons.

### A4. Policies and a common total-burden endpoint

During the course M uses `u=u_max`. C requests killing

$$
k_C^N(r)=\operatorname{clip}(r_S-\mu b_S-a_{SR}r,0,k_{\max}),
\quad
u_C=K_u\left(\frac{k_C^N}{k_\infty-k_C^N}\right)^{1/h}.       \tag{7}
$$

At the clipping endpoints take `u_C=0` or `u_max`, respectively.
Under the inequalities below, clipping is inactive before progression and
C keeps **sensitive density** `s=s0`. This is a sensitive-reservoir
containment policy; total burden may increase from `s0` to `B=s0+epsilon`.
It is not ideal instantaneous total-burden containment. It requires
measurement of clonal composition and allows no replenishment from an
external sensitive reservoir.

Fix `epsilon>0`, `B=s0+epsilon`. For both policies let

$$
T_X=\inf\{t\ge0:s(t)+r(t)\ge B\},\quad X\in\{C,M\}.         \tag{8}
$$

Continue observing after `H_N`; there is no independent censor clock and no
conditioning on resistant establishment. An R birth can cross the boundary
in the stochastic stage. Initial conditions and B are identical in both
arms. This cell-burden endpoint is not a radiological response criterion.

### A5. Finite-N admissibility and bulk persistence

Write `delta=k_max-rS`. Assume

$$
\delta>0,\qquad r_S>a_{SR}B,\qquad
(\delta+a_{SR}B)L<1.                                      \tag{9}
$$

For finite N also require `r_N=rS-lambda*bS/N>aSR*B`.
The latter holds for all sufficiently large N at each fixed finite lambda.
These restrictions keep surviving S as a large bulk even in the aggressive
arm; this is not an eradication regime.

### A6. Additional assumptions for the nontrivial boundary

For Theorems 2--4 with positive released growth assume

$$
g_0>\beta B,\quad L>1/g_0,\quad
g_C=g_0-\alpha s_0,\quad
\big[g_C>\beta\epsilon\ \text{or}\ g_C\le0\big].             \tag{10}
$$

The excluded interval `0<gC<=beta*epsilon` can allow a sub-progression
resistant equilibrium. It must not be inserted into a formula that assumes
R reaches `epsilon`. Equality `gC=beta*epsilon>0` is excluded too.
The main scaling fixes `0<lambda<infinity` while `N` tends to infinity;
`lambda=0` is treated separately. No uniformity in unbounded lambda is claimed.

## 3. Theorem 1 — Admissibility, finite means, and no subcell reservoir

**Theorem 1.** Under A1--A5, for finite N satisfying the finite-N condition,
both controls are admissible, all event rates are nonnegative, and C has
`s(t)=s0` until `min(T_C,H_N)`. Put `q_N=r_N-aSR*B>0`. For either policy,

$$
T_X\le H_N+
\frac{\log(B/s_0)+(\delta+\lambda b_S/N+a_{SR}B)H_N}{q_N}.
                                                               \tag{11}
$$

Consequently both ordinary expectations are finite. Before progression,
throughout the course, with fixed finite lambda in the displayed limit,

$$
Ns(t)\ge Ns_0\exp[-(\delta+\lambda b_S/N+a_{SR}B)H_N]
=s_0N^{1-(\delta+a_{SR}B)L-o(1)}\longrightarrow\infty.       \tag{12}
$$

**Proof.** Before progression `0<=r<B`. The requested C kill lies in
`[r_N-aSR*B,r_N]`, a positive interval strictly below `k_max`, since
`r_N<=rS<k_max`. The inverse Hill law exists there and is at most `u_max`.
Substitution into (3) gives `s'=0`, also between discrete R events; s does
not jump at those events. M also satisfies the dose bound. Integrating
`0<=u<=u_max` over at most `H_N` proves (6). Birth and death rates are
products or sums of nonnegative quantities. Positive ODE solutions remain
positive; at `r=0` the drift in (4) is nonnegative.

For either policy during the course,
`s'/s >= -(delta+lambda*bS/N+aSR*B)`, yielding (12). After the course,
while total burden is below B, `s'/s>=q_N`; sensitive regrowth alone reaches
B within the additional time in (11). An earlier R-caused crossing can only
shorten this bound. Local event rates are bounded before progression, so
there is no finite-time explosion. The strict exponent in (9) proves the
last limit and rules out using a disappearing subcell S reservoir to obtain
the result. Equation (11), divided by `log(N)`, also gives a deterministic
uniform bound for the expectation limit in Theorem 3. QED.

## 4. Theorem 2 — Establishment with declining competition

This theorem concerns the auxiliary rare-resistance process used to connect
the hybrid to macroscopic progression. It is not itself a progression claim.

**Theorem 2.** Under A1--A6 with `g0>0`, consider the maximum-dose rare
sensitive trajectory `s_M(t)=s0*exp(-delta*t)` continued for all time.
A single resistant founder at time t has birth rate `bR` and time-dependent
death rate `dR+alpha*s_M(v)`, `v>=t`. Its survival probability is

$$
\boxed{p(s)=\left[b_R\int_0^\infty
e^{-g_0v+(\alpha s/\delta)(1-e^{-\delta v})}\,dv\right]^{-1},
\qquad p_M(t)=p(s_M(t)).}                                   \tag{13}
$$

Define

$$
A=\int_0^\infty b_Ss_M(t)p_M(t)\,dt
=\frac{b_S}{\delta}\int_0^{s_0}p(s)\,ds>0.                  \tag{14}
$$

Independent mutant founders arrive at rate `lambda*bS*s_M(t)`.
The probability that **none of these on-course-limit lineages survives** is

$$
q=e^{-\lambda A}.                                         \tag{15}
$$

This is neither cure probability nor permission to stop following S.
An exact positive-series representation is

$$
p(s)=\frac{g_0}{b_R F(s)},\qquad
F(s)=\sum_{k=0}^\infty
\frac{(\alpha s)^k}{\prod_{j=1}^k(g_0+j\delta)}.             \tag{16}
$$

If `alpha=0`, `p=g0/bR` and
`A=bS*s0*g0/(delta*bR)` in closed form. In particular `A` depends on
division rate and demographic turnover, not just net growth.

**Proof.** For a founder born at t, set
`G(v)=exp(integral_t^v [bR-dR-alpha*s_M(w)] dw)` and
`I(v)=integral_t^v bR/G(w) dw`. Its probability generating function is

$$
f(z,v)=1-\frac{G(v)(1-z)}{1+G(v)I(v)(1-z)}.                 \tag{17}
$$

Indeed `G'=gG`, `(GI)'=gGI+bR`; differentiation verifies the forward
equation `f_v=(z-1)(bR*z-d(v))*f_z`, with `f(z,t)=z`.
Thus `P(Z(v)>0)=G/(1+GI)`. Since `G(v)` grows at asymptotic rate `g0>0`,
`I(infinity)` is finite and this probability decreases to `1/I(infinity)`.
Integrating `g(v)=g0-alpha*s_M(t)*exp(-delta*(v-t))` yields (13).

Poisson arrivals can be independently marked by survival because their
branching families do not interact in this auxiliary process. The surviving
founders are an inhomogeneous thinned Poisson process; its integrated
intensity is `lambda*A`. This proves (15). The sensitive-change substitution
`ds=-delta*s*dt` proves (14).

For (16), expand the nonnegative exponential
`exp[(alpha*s/delta)*(1-exp(-delta*v))]` in (13), and interchange sum
and integral by monotone convergence. With `x=exp(-delta*v)`,

$$
\int_0^\infty e^{-g_0v}(1-e^{-\delta v})^k\,dv
=\frac{k!\,\delta^k}{g_0(g_0+\delta)\cdots(g_0+k\delta)}.
$$

This beta-integral identity also follows by k integrations by parts; all
boundary terms vanish except the `k=0` term. Substitution gives (16).
Every denominator is positive. The ratio of successive terms tends to zero,
so the series converges for every finite s. QED.

### A quantitative bridge to the actual stochastic stage

Let `m=m_N`, `H=H_N`. Couple the rare branching/immigration process to the
hybrid using common portions of each event rate, stopping at H, first R=m,
progression, or the first mismatched event. For all sufficiently large N,
`m_N/N<epsilon`, so progression cannot precede this switch during the course.
Then

$$
\begin{split}
P(\text{mismatch in C})&\le D_C=\beta m^2H/N,\\
P(\text{mismatch in M})&\le D_M=D_C+
\frac{s_0(\lambda b_S+\alpha m)(\lambda b_S+a_{SR}m)}{N\delta^2}.
\end{split}                                                \tag{18}
$$

These are genuine finite-N bounds **on this coupling only**, not error bars
for the entire progression-time asymptotic formula.

To prove them, C has exactly `s=s0`, so only additional R self-competition
deaths differ, by at most `beta*m^2/N`. For M, on matched paths,

$$
s_N(t)=s_0e^{-\delta t}
\exp\left[-\frac{\lambda b_S t}{N}
-\frac{a_{SR}}N\int_0^tR(v)dv\right].                     \tag{19}
$$

Using `1-exp(-x)<=x`, its deficit from `s_M` is at most
`s0*exp(-delta*t)*(lambda*bS+aSR*m)*t/N`.
The absolute difference of mutant-arrival and competition-death rates is
at most `(lambda*bS+alpha*m)*abs(s_N-s_M)+beta*m^2/N`.
Integrate and use `integral_0^infinity t*exp(-delta*t)dt=1/delta^2`.
The probability of at least one mismatched Poisson event is at most its
integrated rate. This proves (18), which tends to zero for `a<1/2`.

## 5. Theorem 3 — Macroscopic progression and ordinary means

**Theorem 3.** Under A1--A6, for fixed `0<lambda<infinity`, let

$$
d=1/g_0,\qquad
c=\begin{cases}\min(L,1/g_C),&g_C>\beta\epsilon,\\L,&g_C\le0,
\end{cases}
\quad J=L(1+\delta/r_S),\quad q=e^{-\lambda A}.             \tag{20}
$$

As N tends to infinity,

$$
T_C/\log N\ \xrightarrow{P}\ c,\qquad
T_M/\log N\ \xrightarrow{d}\
\begin{cases}d,&\text{probability }1-q,\\J,&\text{probability }q.
\end{cases}                                               \tag{21}
$$

Their expectations converge to (1). At continuity points h of these limits,

$$
\lim P(T_C>h\log N)=\mathbf1_{h<c},\qquad
\lim P(T_M>h\log N)=
\begin{cases}1,&h<d,\\q,&d<h<J,\\0,&h>J.
\end{cases}                                               \tag{22}
$$

Thus the mean comparison does not establish stochastic dominance. C has
higher control probability in the early interval `d<h<c`; M has the longer
tail in `c<h<J` whenever `alpha>0`. These are mathematical observation
times, not clinical survival claims.

### Complete proof

**Step 1: logarithmic growth of a surviving family.** For (17), `Z(v)/G(v)`
is a nonnegative martingale and
`E[(Z/G)^2]=1/G+2I`. The second moment is bounded when `g0>0`, so the
martingale has a finite limit W in mean square and almost surely. Substituting
`z=exp(-theta/G)` into (17) gives
`E exp(-theta W)=1-theta/(1+I(infinity)*theta)`.
Hence W has mass `1-1/I(infinity)` at zero and, on positivity, an exponential
distribution with mean `I(infinity)`. The extinction probability from (17)
is the same mass at zero. Since extinction implies `W=0`, W is positive
almost surely on nonextinction. It follows that
`log(Z(v))/v -> g0` on survival. The same statement holds with constant
positive `gC`.

This also implies that the first hitting time of a large integer m, on
survival, divided by `log(m)` tends to the inverse growth rate: beyond a
finite random time the population lies between exponentials with rates
`g-eta` and `g+eta`; its earlier finite maximum does not reach unbounded m.
On extinction the lifetime and maximum are finite almost surely.

**Step 2: immigration and the stochastic switch.** M's auxiliary process
has a finite Poisson number of immigrants, with mean
`lambda*bS*s0/delta`. On the event E that any family survives, a finite sum
of such families has logarithmic growth rate `g0`; on its complement all
families die at a finite last extinction time. Therefore its hitting time
of `m_N` divided by `log(N)` tends to `a/g0` on E and is absent for all
large N on the complement. Theorem 2 gives `P(E)=1-q`.

For C with `gC>0`, arrivals continue at rate `lambda*bS*s0`. Each has
survival probability `gC/bR>0`, so the first surviving family arrives at a
finite time almost surely. The aggregate has logarithmic growth rate gC.
For an upper bound, `exp(-gC*t)*(R+lambda*bS*s0/gC)` is a nonnegative,
mean-square-bounded martingale: its quadratic-variation expectation is
bounded by a constant times `integral exp(-gC*t)dt`. The first surviving
family provides a strictly positive lower limit. This proves the claimed
growth rate and switch coefficient `a/gC`.

If `gC<=0`, dominate C's rare process by a critical birth-death process
with the same immigration. That dominating process is a nonnegative
submartingale with mean `lambda*bS*s0*t`. The maximal inequality gives
`P(sup_[0,H_N] R>=m_N)<=lambda*bS*s0*H_N/m_N ->0`.
Thus no switch occurs during the course with probability tending to one.
The bounds (18) transfer all of these stopped conclusions to the hybrid.

**Step 3: from the switch to a macroscopic population.** In C, before
course end, (4) becomes
`r'=r*(gC-beta*r)+lambda*bS*s0/N`. For any fixed sufficiently small
`rho>0`, between `m_N/N` and rho the per-capita growth rate is between
`gC-beta*rho` and `gC+lambda*bS*s0/m_N`. The logarithmic distance is
`(1-a)*log(N)+O(1)`. Its traversal time divided by `log(N)` is squeezed
to `(1-a)/gC` by first taking N to infinity and then rho down to zero.
From rho to epsilon the time is bounded independently of N, since
`gC-beta*epsilon>0`. Adding the stochastic-stage coefficient gives `1/gC`
if this is before course end.

For M on E, the sensitive density at the switch tends to zero by (19),
because its switch time is asymptotic to `a*log(N)/g0`. It continues
decreasing during the course. While `m_N/N<=r<=rho`, its R per-capita
growth lies between `g0-alpha*s(switch)-beta*rho` and
`g0+lambda*bS*s0/m_N`. The identical squeeze gives `(1-a)/g0` for this
stage. From rho to the total-burden crossing takes bounded time because
`g0-beta*B>0`, and S remains vanishingly small. The total coefficient is
`d=1/g0<L`, so progression indeed occurs before dose-off.

If C would reach progression only at or after the course coefficient L,
there is no earlier macroscopic crossing by the preceding squeeze. At H_N,
S is still s0 on every unprogressed path. Theorem 1 then bounds additional
time by `log(B/s0)/q_N=O(1)`. This gives `c=L`, including `1/gC=L` and
the critical/subcritical case.

**Step 4: no established lineage is not cure.** On the auxiliary event
`E^c`, the finite family collection has a finite last extinction time and
finite `integral R(t)dt`, almost surely. For large N it has not hit m_N.
By (18), with probability tending to one the hybrid has `R(H_N)=0` and

$$
s_N(H_N)=s_0N^{-\delta L}(1+o_P(1)).                        \tag{23}
$$

This follows directly from (19), since `lambda*bS*H_N/N` and the integrated
R correction tend to zero. Fix `0<eta<delta*L`. After dose-off, until the
next mutant arrives, S grows at `r_N`. The cumulative **first-arrival
hazard along this mutation-free path**, before s reaches `N^{-eta}`, is at
most `lambda*bS*N^{-eta}/r_N ->0`. The probability of an earlier mutant is
at most that hazard, by `1-exp(-x)<=x`. That density is reached at
`[J-eta/rS]*log(N)+o_P(log(N))`. There is no progression before it with
probability tending to one. Afterwards, regardless of new mutations,
`s'/s>=q_N` before progression, so the remaining time is at most
`[eta*log(N)+log(B)]/q_N`. These bounds squeeze the limiting coefficient
between `J-eta/rS` and `J-eta/rS+eta/(rS-aSR*B)`. Let eta decrease to zero.
This proves coefficient J while allowing all post-course mutations.

**Step 5: pass to ordinary expectations.** The pathwise bound (11) divided
by `log(N)` is uniformly bounded for large N. Consequently the variables
in (21) are uniformly integrable; their means converge. The distribution
limits give (22) at continuity points. QED.

## 6. Theorem 4 — A sharp boundary, and regimes with no reversal

**Theorem 4.** Under A1--A6 with `alpha>0`, fixed `lambda>0`, and (20),
the limit of `(E T_C-E T_M)/log(N)` is

$$
\Delta(\lambda)=(c-d)-(J-d)e^{-\lambda A}.                  \tag{24}
$$

It is strictly increasing, has exactly one positive zero (2), and has the
sign stated at the beginning of this note. A strict sign implies the same
finite-N ordering eventually; no convergence rate for that sign is asserted.

The following no-reversal statements have separate, explicit scopes:

1. **No competition.** Under A1--A6 with `alpha=0`, maximum dose has larger
   limiting expected progression time for every fixed finite `lambda>0`.
2. **No supercritical resistant growth.** Under A1--A5 with `g0<=0`
   (replace A6 entirely), the mean coefficients are C=L and M=J for every
   fixed finite lambda. Thus maximum dose wins throughout this regime.
3. **No mutation.** At `lambda=0`, for all finite N, the exact times are
   `T_C=H_N+log(B/s0)/rS` and
   `T_M=H_N+(delta*H_N+log(B/s0))/rS`; M wins strictly.
4. **No two positive-supply reversals in this scaling.** Under the first
   paragraph's assumptions, two or more sign changes as lambda increases
   are impossible for the limiting means. This says nothing about an
   unproved finite-N phase curve or about a different toxicity/course rule.

**Proof.** Positive competition gives `gC<g0`, so A6 implies
`d<c<=L<J`; (14) gives `A>0`. Subtract (1) to get (24).
Its derivative is `(J-d)*A*exp(-lambda*A)>0`.
Its limit at `lambda down to0` is `c-J<0`, and its limit at infinity
is `c-d>0`. Solving the unique zero gives (2). Positivity of all factors
also fixes the inequality direction. Convergence of the normalized mean
difference implies eventual finite-N strict ordering whenever its limit
is nonzero.

For `alpha=0`, `c=d`, so `Delta=-(J-d)*exp(-lambda*A)<0` for every finite
positive lambda. For `g0<=0`, (11) bounds T by `K*log(N)` for a fixed K
and all sufficiently large N. Until progression or the switch, dominate R
by a critical branching process with immigration `lambda*bS*B` throughout
this entire interval, including after dose-off. The maximal inequality gives
`P(hit m_N)<=lambda*bS*B*K*log(N)/m_N ->0`. On the complementary event,
R stays below m_N and `aSR*integral R/N <=aSR*m_N*K*log(N)/N ->0`.
For M, integrating (3) across the course and rebound therefore yields
`log(s(T)/s0)=r_N*T-k_max*H_N+o_P(1)`, while at progression
`s(T)=B-o_P(1)`. Thus its coefficient is `k_max*L/rS=J`.
C has `s=s0` during the course and reaches B within bounded time
afterwards, so its coefficient is L. This argument does not assume
`R(H_N)=0` in the critical case. Bound (11) justifies expectations again.
Finally, at lambda=0
there is no R; solve the two sensitive exponentials before and after H_N
to get the exact times. Strict monotonicity proves item 4. QED.

**Noncommuting boundary.** When `1/gC<L`, the limit of the C coefficient
as positive lambda decreases to zero is `1/gC`, whereas at exactly zero
it is L. For every fixed positive supply, enough logarithmic time eventually
allows a successful C founder; this need not happen at `10^8` cells when
supply is tiny. Do not insert zero into (1).

### Measurable quantities and dose-response steepness

The criterion uses `N*mu`, the sensitive division rate `bS`, resistance
birth/death rates, competition `alpha`, reciprocal/self-competition feasibility
checks, and the course scale L. They can be estimated separately in a
nonclinical experiment. Four summary labels alone cannot identify the
boundary: net fitness does not identify division-driven mutation supply.

The Hill exponent h enters through
`delta=k_infinity*u_max^h/(K_u^h+u_max^h)-rS` and through the dose needed
for (7). Changing h at fixed `K_u,u_max,k_infinity` can change admissibility,
A and J; there is no universal monotonic “steeper is better” result.
If `k_max` and all demographic parameters are held fixed, and both inverse
controls remain admissible, h changes the dose/exposure parametrization but
does not change the population equations or the comparison. This follows
by substituting the **realized kill rates**, not by treating h as an
independent biological advantage. In the numerical example `u_max=K_u`,
so `k_max=k_infinity/2` for every positive h.

## 7. Reproducible certificates and independent simulation

All code is Python standard library only. It reads no external data, makes
no network requests, and writes no results. From the repository root:

```sh
python research/containment-at-tumour-scale/verifier.py
python research/containment-at-tumour-scale/simulation.py --samples 400 --supplies 0.005 1 --populations 100000000 100000000000
python -m unittest discover -s tests -p test_research_containment_at_tumour_scale.py -v
```

### Exact arithmetic, rather than an unqualified numerical integral

`verifier.py` uses `fractions.Fraction` for scientific certificates;
`decimal.Decimal` only formats bounds. Its Hill calculation accepts positive
integer h for exact rational arithmetic; the theorem permits any real `h>0`.
For a partial sum of (16) through K,
the remaining tail is bounded above by the next term divided by
`1-alpha*s/(g0+(K+2)*delta)`, once that ratio is below one. All terms are
positive, so reciprocation gives rigorous lower and upper bounds for p.

Composite Simpson quadrature of (14) is enclosed, not assumed exact.
Writing `p=1/I`, the integral in (13) gives
`0<=I^(j)/I<=(alpha/delta)^j`. Differentiating `I*p=1` recursively gives
`abs(p^(j)) <= C_j*(alpha/delta)^j*p`, where
`C_j=sum_{k=1}^j binomial(j,k)*C_(j-k)` and
`C_0,...,C_4=1,1,3,13,75`. Since `p<=g0/bR`,

$$
\max_{[0,s_0]}|p^{(4)}|\le75(\alpha/\delta)^4g_0/b_R.
$$

For an even number n of panels the absolute integral error is at most
`s0*(s0/n)^4*max(abs(p''''))/180`. Series intervals and this remainder are
combined with outward rational rounding. Logarithms use the atanh series
with a geometric tail bound; negative exponentials use reciprocals of
positive Taylor-series enclosures. Exact arithmetic and explicit remainders
cover the quoted numerical certificates, not the uncomputed biological
parameter uncertainty.

The dedicated test also integrates the **backward branching Riccati equation**
independently of the series implementation, checks the hazard against
stochastic resistant births/deaths, checks the finite-course zero-mutation
solution, and tests dose bounds, numerical step refinement, scaling slopes,
and the two ranking witnesses. A test does not prove the asymptotic theorem. The dedicated file has 16
tests and completed locally in 22.95 seconds, below the 60-second target;
hardware affects runtime. The standalone exact verifier completed in about
0.5 seconds. No optional packages are required.

### Synthetic parameter set and finite-scale reporting

| Parameter | Value |
| --- | --- |
| `bS,dS,bR,dR` | `1, 0.7, 0.8, 0.6` |
| `alpha,beta,aSR` | `0.12, 0.02, 0.02` |
| `s0,epsilon,B` | `1, 0.2, 1.2` |
| `k_infinity,K_u,u_max,h` | `0.66, 1, 1, 2` |
| `k_max,rS,delta,g0,gC` | `0.33, 0.3, 0.03, 0.2, 0.08` |
| `L,a` | `16, 1/4` |
| `d,c,J` | `5, 12.5, 17.6` |
| Bulk-persistence exponent lower bound | `1-(delta+aSR*B)*L=0.136>0` |

Time is in arbitrary common units. Changing their unit rescales all rates
and course coefficients consistently. These numbers were selected to test
the mathematics, not calibrated to an organoid, xenograft or patient.

`simulation.py` does not import the verifier or sample the limiting mixture.
Before switching it simulates individual R birth/death/mutation events.
Between events, R is fixed and S evolves exponentially, so the total hazard
is `A0+A1*exp(v*t)` and is integrated and inverted directly. After switching,
an independent fourth-order Runge-Kutta integrator evolves (3)--(4), with
event localization at the total-burden boundary. Step-halving checks isolate
ODE error. The simulation includes post-course mutant arrivals and follows
actual regrowth.

The default exact enclosures (256 Simpson panels) are

$$
6.27019085952677871<A<6.27019127344792300,
$$
$$
0.0827397077362023465<\lambda_*<0.0827397131981915061.
$$

At `lambda=0.1`, the limiting C coefficient is exactly `12.5`, and M's
coefficient is enclosed by
`[11.7306901494782831,11.7306904280757858]`. The coupling bounds (18) at
that supply are:

| N | m_N | C mismatch upper bound | M mismatch upper bound |
| --- | --- | --- | --- |
| `10^8` | 100 | 0.000589462 | 0.000871796 |
| `10^9` | 178 | 0.000210111 | 0.000297382 |
| `10^10` | 317 | 0.000074044 | 0.000101335 |
| `10^11` | 563 | 0.000025691 | 0.000034231 |

Displayed upper bounds are rounded upward. They only certify the stopped
branching coupling, as specified in (18).

The following experiment used **400 independent trials per arm**, base seed
`20261117`, `ode_step=0.5`, and `a=1/4`. The CLI increments its seed in
population, supply, C/M order. Reproduce the complete table with:

```sh
python research/containment-at-tumour-scale/simulation.py --samples 400 --seed 20261117 --supplies 0.005 1 --populations 100000000 1000000000 10000000000 100000000000
```

| N | lambda | C mean +/- SE | M mean +/- SE | Larger mean |
| --- | --- | --- | --- | --- |
| `10^8` | 0.005 | 292.9873 +/- 0.5964 | 321.1592 +/- 1.3230 | M |
| `10^9` | 0.005 | 328.4378 +/- 0.8267 | 356.4680 +/- 2.1850 | M |
| `10^10` | 0.005 | 364.8894 +/- 0.8989 | 399.1440 +/- 2.0161 | M |
| `10^11` | 0.005 | 403.0534 +/- 0.7402 | 437.8612 +/- 2.4400 | M |
| `10^8` | 1 | 183.1552 +/- 0.5659 | 105.8332 +/- 0.4222 | C |
| `10^9` | 1 | 213.1352 +/- 0.6452 | 117.5920 +/- 0.4909 | C |
| `10^10` | 1 | 241.7221 +/- 0.6909 | 130.2346 +/- 0.4762 | C |
| `10^11` | 1 | 270.4996 +/- 0.6874 | 142.4217 +/- 0.8813 | C |

SE means the sample standard error across independent trials. Difference
uncertainty uses the square root of the sum of the two variance estimates;
these are not paired trials. The low-supply sign margins exceed 12 such SEs
at every size, and the high-supply margins exceed 100. The simulated M
on-course switch frequencies at `lambda=0.005` were
`0.0200,0.0400,0.0275,0.0300`; the auxiliary asymptotic escape probability
is about `0.03087`. No establishment formula was used to generate these
events. A frequency of one in a finite sample does not establish probability
one; tests use the theoretical Bernoulli variance when checking this limit.

The maximum C dose was about `0.9128709292`; M used the bound 1. Exposure
remained within budget up to floating rounding below `10^-8`. The smallest
observed sensitive population was about 14,454 cells. Equation (12)'s crude
worst-case lower bound is much smaller (about 12 cells at `10^8`); it proves
asymptotic divergence, not a finite-N error certificate for the S fluid model.

**The limiting means are not already accurate finite-N means.** At lambda=1,
the normalized C means move from `9.9429` to `10.6797` across the table,
while the limit is `12.5`; M moves from `5.7453` to `5.6230`, while its limit
is about `5.02384`. The high-supply difference in C mean between `10^8` and
`10^11`, divided by `log(1000)`, is about `12.64`, consistent with its growth
slope while retaining a large startup offset. At lambda=0.005, C's
normalized mean remains about `15.9` rather than `12.5`: the first
successful-founder waiting time is still important. These are recorded
finite-size effects, not hidden discrepancies.

### Attempt to falsify use of the limiting boundary at finite N

With 200 trials per arm, base seed `20262117`, supplies `0.08 0.1`, and
populations `100000000 100000000000`, the following near-boundary results
were obtained:

| N | lambda | C mean +/- SE | M mean +/- SE |
| --- | --- | --- | --- |
| `10^8` | 0.08 | 254.8743 +/- 2.9491 | 237.6194 +/- 7.0089 |
| `10^11` | 0.08 | 355.6719 +/- 3.4383 | 325.1178 +/- 9.8916 |
| `10^8` | 0.1 | 250.2631 +/- 3.0006 | 225.5011 +/- 6.9871 |
| `10^11` | 0.1 | 349.7282 +/- 3.4507 | 317.1807 +/- 10.0955 |

At `lambda=0.08<lambda_*`, the limiting sign favours M, whereas the sample
means favour C by only about 2.27 and 2.92 difference SEs. This is **moderate
evidence of a finite-size reversal**, not a certified finite-N sign. It
rejects treating (2) as an already-converged practical threshold; it does not
contradict a theorem about the limit. The finite-N phase boundary remains
unproved and was not fitted after seeing these outcomes. The test retains
robust far-from-boundary witnesses rather than requiring this noisy result.

## 8. Falsifiable nonclinical prediction

In a serial organoid experiment, use barcoded initially sensitive clones,
estimate `mu` per completed division using an independent fluctuation assay,
measure birth/death turnover separately, and measure resistant invasion
growth against several sensitive densities to estimate `g0` and `alpha`.
Fit the same dose-kill curve in both arms. Prespecify `B`, the finite course,
the common exposure ceiling, and total-burden progression **including after
dose withdrawal**. Do not redefine progression after seeing the trajectories.

The quantitative prediction is that the proportion of maximum-dose cultures
with no expanding resistant lineage by course end has logarithm approximately
`-N*mu*A`, while progression times separate into early resistant expansion
and later sensitive rebound. If the scaling regime applies and competition
is positive, increasing independently measured supply should produce **one**
crossing of the mean advantage, with boundary (2); eliminating competition
should remove the containment-winning region. Serial xenografts could test
the same relationship using clonal measurements and an approved nonclinical
protocol. No new animal experiment was performed or commissioned here.

Estimate rates on separate training cultures and predict a held-out panel
with uncertainty propagated before comparing outcomes. Repeated crossings,
substantial R expansion when measured `g0<=0`, failure of the early/late
mixture, or a persistent disagreement with the independently parameterized
**finite-size simulation**, beyond its sampling and numerical uncertainty,
would falsify this model in that system. A discrepancy from the asymptotic
boundary alone is insufficient at finite size, as the recorded near-boundary
attempt demonstrates. Small organoids cannot
alone validate the `10^8`--`10^11` extrapolation; spatial structure and the
applicability of well-mixed competition must also be tested. Do not vary
mutation rate by a drug that simultaneously changes birth, death or
competition and attribute the result solely to mutation supply.

## 9. Hostile review, limitations, and novelty assessment

The candidate contribution is a **model-specific mathematical refinement**:
the finite-course progression limit, the explicit one-boundary criterion,
the coupling and arithmetic certificates, and independent reproducibility.
Priority is unestablished. The ecological trade-off, branching establishment,
Poisson thinning, Hill dose response, bounded-dose containment, toxicity
constraints and objective dependence are prior work.

### Referee ledger and retained negative results

| Attempt or audit | Disposition |
| --- | --- |
| Claim the broad mutation/competition trade-off or tumour-sized analysis is new | Rejected by the primary-source audit. The claim is restricted to the specified theorem and implementation. |
| Use a fixed observation time with `N*mu` fixed and call establishment progression | Rejected: R stays microscopic on fixed time intervals. The logarithmic bridge is proved instead. |
| Use indefinite dosing to compare ordinary means | Rejected: noescape/cure atoms can make an ordinary mean infinite. Both policies now stop after a stated finite course and all post-course regrowth is followed in the hybrid. |
| Treat a vanishing deterministic S remnant as a reservoir | Excluded by (9) and (12); no sensitive immigration or restoration is present. This does not prove a full discrete S model. |
| Substitute finite-N means for their leading asymptotic coefficients | Rejected by the numerical table and near-boundary falsification attempt. No exact finite-N phase claim remains. |
| Keep `0<gC<=beta*epsilon` in the macroscopic-growth formula | Rejected: the resistant equilibrium can lie below progression. The validator and A6 exclude this region, including equality. |
| Call a first-arrival hazard the expected count after possible feedback | Corrected during referee review. Step 4 explicitly follows the mutation-free path only until the first arrival. |
| Assume critical resistant lineages are extinct at course end | Rejected. Theorem 4 instead uses critical domination over the full pathwise time bound, including rebound. |
| Ignore progression in the finite-N coupling stopping rule | Corrected. The coupling stops at progression; for large N the switch occurs first during dosing. |
| Reverse the phase inequality or claim a tie at finite N | Exact sign witnesses, positivity of A, and the strictly positive derivative in (24) prevent the former; the latter is explicitly unproved. |
| Introduce a universal Hill-steepness benefit | Rejected. Holding the realized kill rates fixed leaves the demographic equations unchanged. |
| Discard failed lineages or turn mutation off after withdrawal | Rejected. The simulator keeps the actual events and the proof bounds post-course first-arrival hazards. |

Independent AI reviewers rederived the PGF, survival series, coupling,
macro-growth squeeze, expectation bound and phase signs, and independently
checked the rational arithmetic against high-precision values and a backward
branching equation. The corrections above were incorporated before pushing.
Their verdict supports the narrow hybrid theorem; it is not external peer
review, a novelty certificate, or biological validation.

Important limitations are mathematical, not merely caveats about calibration:

- The theorem is for the stated hybrid. In a fully discrete S/R process,
  even an extremely small positive cure probability can make the literal
  uncensored progression mean infinite. Vanishing extinction probability
  does not by itself justify convergence of those means. No full-CTMC
  ordinary-mean theorem is claimed.
- The deterministic S bulk is kept large by (9), but stochastic S fluctuations
  and a full discrete-to-hybrid error bound are not established here.
- The mesoscopic coupling bound (18) does not bound discarded demographic
  fluctuations after the switch or give a finite-N error for (1).
- The finite course and its `log(N)` budget scaling are essential. A fixed
  lifetime toxicity budget, equal-exposure policies, adaptive course lengths,
  repeated courses, pharmacokinetic lag and a fixed clinical follow-up
  require new analysis.
- C stabilizes S, not total burden exactly. Instantaneous genotype measurement
  and quasi-steady drug response are assumptions. A delayed/noisy controller
  is not covered, nor is global optimality over all admissible schedules.
- A fixed `N*mu` limit is not the same as increasing N at fixed mutation
  probability. No uniform-in-supply conclusion or high-`mu` extrapolation is
  justified. At the boundary (2), finite-size corrections decide the ranking.
- Pre-existing resistance, heterogeneous costs, persisters, spatial refuges,
  resource feedback and density-dependent division rates can change the
  result. Their absence is not evidence that they are negligible in cancer.

The closest remaining open mathematical variant is a **fully discrete,
finite-PK model with a fixed cumulative budget and a bounded progression
objective**, including pre-existing resistance and spatially measured
competition. The present note does not label that question solved.
