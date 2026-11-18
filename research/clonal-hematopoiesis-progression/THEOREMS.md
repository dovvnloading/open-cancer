## Mathematical model, assumptions, and proofs

**Status: nonclinical, unreviewed mathematical derivations.** A terminal event
below means the specified model endpoint. It is diagnosed AML only if the
completion/detection model is calibrated to that endpoint. No such clinical
calibration is performed here. The branching equations, thinning arguments,
Laplace-transform reasoning, and rare-event limits are classical methods;
writing out their implications for this audit is not a claim to have invented
branching-process cancer theory. Numerical checks check implementations and
counterexamples, not the truth of a biological model.

### Definitions and assumption register

Time is in years. The current age is `a`, the prediction horizon is `T >= 0`,
and `v` is the currently observed VAF. The unmarked clone has per-cell birth
rate `b` and death/differentiation rate `d`; `s = b-d` is its continuous-time
net growth rate. A birth adds one same-type cell, and a death removes one.
The per-cell terminal-mark rate is `q`. The cumulative hazard is
`H(T) = -log S(T)`, risk is `1-S(T)`, and instantaneous event-free hazard is
`h(T) = dH(T)/dT`. These three quantities are different.

The assumptions are deliberately stronger than the available aggregate data.

| Code | Assumption | Where used |
|---|---|---|
| A1 | At baseline there are exactly `n` relevant stem/progenitor cells, a nonnegative integer; their future descendant processes are independent and have identical specified rates. | Theorems 1–3, 5, and 6. |
| A2 | Birth and death form a nonexplosive linear continuous-time branching process with finite nonnegative per-cell rates. Rates are constant unless explicitly allowed to depend on age. There is no density limitation or competition in this model. | Theorems 1–8; Theorem 9 explicitly replaces this assumption. |
| A3 | Conditional on the unmarked clone trajectory, terminal marks are an independent Poisson process of intensity `q Z(t)`; marks do not change the unmarked process used in this construction. There is no prior terminal event. | Theorems 1, 2, 5, and the turnover counterexample. |
| A4 | To equate `n` with `2Nv`, there is a known fixed baseline stem-cell pool `N`, a single heterozygous copy-neutral mutation, unbiased lineage contribution, and negligible measurement error. `2Nv` must be an integer for an exact cell-count interpretation. | Only statements translating cell count into VAF. The stochastic results themselves do not require this mapping. |
| A5 | Delayed pathways are initiated by independent Poisson events from type-0 cells; every initiated pathway evolves independently, does not feed back into type-0 growth, and has the specified possibly defective completion-time distribution. No advanced pathway is already present at baseline. | Theorem 3. |
| A6 | Full multitype branching dynamics, all stage-specific rates, and the current counts at every stage are specified; a transition **replaces** one cell with a next-stage cell. | Theorem 4. |
| A7 | A fitness-only comparison matches current cell count (or VAF **and** `N`), age, outcome definition, horizon, and completion kernel; gene-specific effective hit coefficients are equal unless explicitly retained. | Theorems 7–8 and the conditional quantitative prediction. |
| A8 | A fitted growth coefficient is transportable to the future horizon and the target population. Its meaning is the one required by the chosen exponential or logistic model. | Every empirical application, never established by an algebraic identity. |
| A9 | The numerical rates and counts used in computational checks are declared synthetic scenarios, not inferred patient parameters or outcome observations. | Implementation checks and simulations. |

Copy-number changes, multiple clones, lineage bias, competing mortality,
time-varying environment, hidden intermediate mutations, and ascertainment can
violate these assumptions. They are not silently absorbed into measured
fitness. A gene-level average fitness does not reveal a person's true rate.

### Theorem 1 — Exact marked branching survival and exact clone-count scaling

**Assumptions:** A1–A3, constant finite `b,d,q >= 0`, finite `T >= 0`. A4 is
needed only for the VAF corollary.

Let `Z(0)=1` and `A_T = integral_0^T Z(t) dt`. Then

$$
u(T)=\mathbb E[e^{-qA_T}],\qquad
S_n(T)=u(T)^n,\qquad H_n(T)=-n\log u(T).
\tag{1}
$$

The unique solution in `[0,1]` of the backward equation is

$$
u'=bu^2-(b+d+q)u+d,\qquad u(0)=1.
\tag{2}
$$

For `q > 0`, write `s=b-d`,
`D = sqrt((s-q)^2+4bq) > 0`, and `E=exp(-DT)`. Then

$$
1-u(T)=
\frac{2q(1-E)}{(D-s+q)+(D+s-q)E}.
\tag{3}
$$

For `q=0`, `u=1`. Formula (3) also covers `b=0` and `s=0`. In particular,
when `b=0`,
`u=d/(d+q)+q exp(-(d+q)T)/(d+q)`, with the `d=q=0` case defined as one.
At fixed rates and horizon the cumulative hazard and instantaneous hazard
are exactly linear in `n`. Under A4 they are exactly linear in current
VAF. **Risk itself is not linear:** `P=1-exp(-H_n)`.

**Proof.** Conditional on the unmarked path, the probability of no Poisson
mark is `exp(-qA_T)`, giving the first identity. For `n` initial independent
cells, the probability is a product of `n` identical single-ancestor
probabilities, giving the second and third identities. During an initial
interval `dt`, one initial cell gives birth with probability `b dt+o(dt)`
(two independent future lineages), dies with probability `d dt+o(dt)`
(certain no subsequent event), or is marked with probability `q dt+o(dt)`
(zero no-event probability). Otherwise its future survival is `u(T)`.
Subtracting `u(T)`, dividing by `dt`, and taking the limit yields (2).
The vector field at `u=0` is `d>=0`, and at `u=1` it is `-q<=0`, so the
solution stays in `[0,1]`; the polynomial is locally Lipschitz, so uniqueness
holds there for every finite horizon. Put `y=1-u`; then
`y'=q+(s-q)y-by^2`, `y(0)=0`. Separating this quadratic equation, or
differentiating (3) and checking its initial condition, yields (3). Its
denominator is positive for `q>0`. The scaling assertions follow directly
from (1) and differentiation in `T`. This proves all cases, including `n=0`.

### Theorem 2 — A finite error bound and the exact rare-hit limit

**Assumptions:** those of Theorem 1. Define

$$
K(s,T)=\int_0^T e^{st}dt=
\begin{cases}(e^{sT}-1)/s&s\ne0,\\T&s=0.\end{cases}
\tag{4}
$$

For one initial cell, `E[A_T]=K(s,T)` and

$$
M_2:=\mathbb E[A_T^2]
=K(s,T)^2+
\begin{cases}
\displaystyle\frac{2(b+d)}{s^2}
\{K(2s,T)-T e^{sT}\}&s\ne0,\\[2mm]
\displaystyle\frac{(b+d)T^3}{3}&s=0.
\end{cases}
\tag{5}
$$

For all finite nonnegative parameters,

$$
0\le nqK(s,T)-H_n(T)\le\frac{nq^2M_2}{2}.
\tag{6}
$$

Thus `H_n/(nq) -> K` as `q -> 0` for fixed `n>=1` and `T`. There is also
a finite-risk limit: with `n -> infinity` and `q_n=kappa/n`,

$$
H_n(T)\longrightarrow\kappa K(s,T),\qquad
P_n(T)\longrightarrow1-e^{-\kappa K(s,T)},
\tag{7}
$$

with cumulative-hazard error at most `kappa^2 M_2/(2n)`.
Large `n` **by itself**, with `q` fixed, does not remove the relative
single-ancestor stochastic correction in (6).

For fixed `T`, the exact instantaneous hazard is

$$
h_n(T)=nq\frac{\mathbb E[Z(T)e^{-qA_T}]}{
                         \mathbb E[e^{-qA_T}]},
\quad
\frac{h_n(T)}{nq}\longrightarrow e^{sT}\quad(q\downarrow0).
\tag{8}
$$

Consequently a rare-hit **instantaneous** hazard ratio at a matched time is
`exp((s_g-s_r)T)`, whereas the **cumulative** hazard ratio is `K(s_g,T)/K(s_r,T)`.

**Proof.** The generator applied to `z` gives `s z`, hence
`E Z(t)=e^{st}`. Applied to `z^2`, it gives
`2s z^2+(b+d)z`; solving this linear moment equation gives
`Var Z(t)=(b+d)e^{st}(e^{st}-1)/s`, with continuous limit `(b+d)t`.
For `0<=r<=t`, the branching/Markov property gives
`E[Z(t)|Z(r)]=Z(r)e^{s(t-r)}`, so

$$
\operatorname{Cov}(Z(r),Z(t))=
\begin{cases}
\displaystyle\frac{b+d}{s}e^{st}(e^{sr}-1)&s\ne0,\\
(b+d)r&s=0.
\end{cases}
$$

Integrating the first moment proves (4). Integrating this covariance over
the two triangles of `[0,T]^2` gives
`Var A_T=2(b+d)[K(2s,T)-T exp(sT)]/s^2`; for `s=0`, direct integration gives
`(b+d)T^3/3`. Adding `(E A_T)^2` proves (5); finite moments justify these
integrations. Jensen's inequality gives
`u=E exp(-qA_T) >= exp(-qK)`, hence `H_n<=nqK`. For `x>=0`,
`exp(-x)<=1-x+x^2/2`, so `1-u>=qK-q^2M_2/2`.
Since `-log u>=1-u`, `H_n>=n(qK-q^2M_2/2)`. Combine with nonnegativity
and the upper bound to obtain (6). Dividing (6), or substituting
`q=kappa/n`, proves both limits. Differentiating the path integral's
Laplace transform gives `u'=-q E[Z(T)e^{-qA_T}]` (the derivative holds at
each time by the finite moment bounds and continuity of expectations).
This yields (8). Dominated convergence, using `E Z(T)<infinity`, proves
its limit. The distinction between the two ratios follows from (4) and (8).

The absolute error in the finite-risk approximation
`1-exp(-nqK)` is also at most `nq^2M_2/2`, because `1-exp(-x)` has
derivative at most one for `x>=0`. The additional approximation
`P approximately nqK` requires `nqK` small; rare per-cell hits alone do
not establish that additional approximation.

### Theorem 3 — Age-dependent rates and explicitly delayed additional hits

**Assumptions:** A1, A2 with bounded measurable rates on `[a,a+T]`, and A5.
An initiation at future time `t` has per-cell rate `epsilon c(a+t)`, with
finite `epsilon>=0` and bounded measurable `c>=0`. Its independent
completion delay has a jointly measurable CDF
`F(ell; a+t)` for `ell>=0`; a defective CDF is allowed, representing
pathways that never complete. Every completion counts as the specified
terminal endpoint. There are no unmodelled pre-existing advanced pathways.
The normalized limit below assumes `n>=1`.

Define `q_T(t)=epsilon c(a+t) F(T-t; a+t)`. Then

$$
S_n(T)=\left\{\mathbb E\exp\left[
-\int_0^Tq_T(t)Z(t)dt\right]\right\}^{n}.
\tag{9}
$$

For one ancestor at relative time `r`, its no-completion probability `u(r,T)`
is the unique absolutely continuous solution of the following backward
equation, with the derivative equality interpreted almost everywhere:

$$
-\partial_r u=b(a+r)(u^2-u)+d(a+r)(1-u)-q_T(r)u,
\quad u(T,T)=1.
\tag{10}
$$

Put `m(t;a)=exp(integral_0^t [b(a+r)-d(a+r)]dr)`. As `epsilon -> 0`,

$$
\frac{H_n(T)}{n\epsilon}\longrightarrow
K_F(a,T)=\int_0^T c(a+t)F(T-t;a+t)m(t;a)dt.
\tag{11}
$$

If `W_T=integral_0^T c(a+t)F(T-t;a+t)Z(t)dt`, the finite error is
`0<=n epsilon E W_T-H_n<=n epsilon^2 E W_T^2/2`.

A concrete additional-hit cascade is an initiation followed by `k`
independent exponential waiting times, with specified rates, and then
the endpoint. For equal rates `lambda>0`,

$$
F_k(\ell)=1-e^{-\lambda\ell}
\sum_{j=0}^{k-1}\frac{(\lambda\ell)^j}{j!},\quad k\ge1;
\qquad F_0(\ell)=1.
\tag{12}
$$

This construction includes arbitrarily many further waiting stages and
non-exponential or defective delays. It does not assume that a clinical
AML diagnosis occurs at the instant an extra mutation appears.

**Proof.** Conditional on a type-0 path, initiations form a Poisson random
measure. Mark each initiation independently by whether its cascade
completes before `T`; that event has probability `F(T-t;a+t)`. Independent
Poisson thinning gives the no-completion probability inside the braces
in (9). Independence across ancestors gives its `n`th power. Conditioning
on a first birth, death, or horizon-completing initiation gives the backward
equation (10), with coefficients evaluated at the current absolute age.
After reversing time, its bounded measurable coefficients and uniform
local Lipschitz dependence on `u` on `[0,1]`, together with
inward-pointing endpoint derivatives, give the unique absolutely
continuous solution in `[0,1]`.
The expectation of the linear branching process solves `m'=s(a+t)m`,
giving `m(t;a)`. Apply the Jensen/Taylor proof of Theorem 2 to the
nonnegative random variable `W_T`. Bounded rates and a finite horizon
give its finite second moment by domination by a finite-rate Yule process.
This proves the finite bound and (11). Finally the sum of `k` independent
rate-`lambda` exponential waiting times is at most `ell` exactly when a
rate-`lambda` Poisson process has made at least `k` arrivals by `ell`;
its probability is (12).

Age appears through the future rates, initiation rates, and delay laws.
If they are time homogeneous and the complete current state is supplied,
age itself contributes no extra information. If a prior cascade may
already be in progress, its current stage and residence time are part of
the state; (9) without those components does not apply.

### Theorem 4 — An exact multitype alternative, without collapsing extra hits

**Assumptions:** A6, `m>=1` nonterminal stages numbered `0,...,m-1`, finite
constant rates `b_i,d_i,mu_i>=0`, and known current counts `n_i>=0`. A
current count is a nonnegative integer at every stage. A stage-i birth
produces two stage-i cells; death leaves no cell; a mutation
replaces the cell with one stage-(i+1) cell. Transition out of stage `m-1`
is the first terminal endpoint. Lineages are independent and there is no
density dependence. There are no additional unspecified pathways.

Let `u_i(T)` be the probability of no endpoint from one current stage-i
cell. Put `u_m(T)=0`. Then the exact equations and survival are

$$
u_i'=b_i u_i^2-(b_i+d_i+\mu_i)u_i+d_i+\mu_i u_{i+1},
\quad u_i(0)=1\quad(i<m),
\tag{13}
$$

$$
S(T)=\prod_{i=0}^{m-1}u_i(T)^{n_i}.
\tag{14}
$$

**Proof.** Condition on the first event of one stage-i cell exactly as
in Theorem 1. A mutation now contributes `u_(i+1)`, while a terminal
mutation contributes zero. This yields (13). At the lower face `u_i=0`
the derivative is `d_i+mu_i u_(i+1)>=0`, and at the upper face `u_i=1`
it is `mu_i(u_(i+1)-1)<=0`. The cube `[0,1]^m` is invariant, and the
polynomial vector field is Lipschitz on that compact cube, proving
existence and uniqueness for every finite horizon. Independent initial
lineages give (14).

This theorem is exact for arbitrary numbers of stages and arbitrary
finite nonnegative stage rates. It does not imply that one founder-gene
fitness identifies all these parameters. If mutation instead produces a
next-stage daughter while retaining its parent, its backward contribution
would be `mu_i u_i(u_(i+1)-1)`, a different model. The code implements
the replacement model stated in (13), with explicit RK4 approximation.

### Theorem 5 — Precisely when clone history can be discarded

**Assumptions:** A1–A3 for the first assertion, or A6 for its multitype
version, all required rate parameters known; the complete state at age
`a` is observed. A4 is additionally required if using VAF as a proxy for
the cell count. For strict fitness monotonicity hold `d` and `q>0` fixed,
set `b=d+s>=0`, and assume `n>=1`, `T>0`.

In the time-homogeneous one-stage model, conditional on current `n` and
the specified rates, neither the past growth trajectory, the founding
time, nor chronological age adds information about future progression.
Within a stratum in which **all those rate parameters and `N` are fixed**,
true VAF can replace `n` under A4. This is the limited sense in which
VAF is sufficient. Across strata with distinct known `s` at fixed `d,q`,
progression probability is strictly increasing in `s` at fixed `n`.
Therefore fitness adds information conditional on VAF whenever it has
nondegenerate conditional variation in that restricted model.

In the rare-hit prediction,

$$
\partial_s\log K(s,T)
=\frac{\int_0^T t e^{st}dt}{\int_0^T e^{st}dt}\in(0,T),
\tag{15}
$$

and at `s=0` this derivative equals `T/2`. At time zero the one-stage
instantaneous hazard is `nq`, independent of fitness; future accumulated
opportunity creates the horizon-dependent fitness effect.

**Proof.** The complete birth/death/mark state is Markov. Conditional on
the current state and known generator, the conditional law of the future
does not depend on the past. The formulas of Theorem 1 make this explicit.
For deterministic age-dependent rates, the same assertion holds given
`a` as well as the current state. In the multitype model the required
state is the full vector `(n_0,...,n_(m-1))`, not just its sum. For
monotonicity, couple two processes by using identical births, deaths,
and marks on all shared cells, and add independent births at rate equal
to the difference in `b` for the process with greater birth rate. Extra
cells have their own descendants and marks. Thus every terminal event
in the lower-birth process occurs in the higher-birth process. If the
birth-rate difference is positive, there is positive probability of an
extra birth before `T`, a subsequent mark in that extra lineage, and no
mark in the shared process, establishing strictness. Differentiating
the finite-interval integral in (4) yields (15). The weighted mean of
`t` is strictly between zero and `T`; direct integration at `s=0` gives
`T/2`. Equation (8) at time zero gives `nq`.

This is a conditional Markov statement, not a claim that a measured bulk
VAF is a sufficient statistic for human AML risk or that chronological
age is irrelevant in real cohorts.

### Theorem 6 — Net fitness and VAF do not identify progression in general

**Assumptions:** the class of models allows the finite nonnegative rates
in Theorem 1 or the finite multitype rates in Theorem 4. Current age,
cell count, and the founder-clone net growth rate are specified, but
turnover, completion rates, or intermediate stages need not be known.

There is no function of `(s,n,a,T)` alone that equals progression
probability throughout this class, even when terminal-mark rates are
fixed and known. There is also no function of those summaries alone
that recovers progression when hidden intermediate stages are allowed.

**Proof, turnover counterexample.** Fix arbitrary `n>=1`, `q>0`, and
`T>0`. Model A has `b=d=0`; model B has `b=d=c>0`. Both have net
fitness zero, the same current cell count and age, and the same
terminal-mark rate. In model A, `A_T=T` deterministically, so
`H_A=nqT`. In model B, `E A_T=T` but `A_T` is nonconstant: no demographic
event gives `A_T=T`, while a death before `T` with no birth gives
`A_T<T`, each with positive probability. Strict Jensen's inequality
therefore gives `E exp(-qA_T)>exp(-qT)` and `H_B<nqT`. Their risks
differ despite identical `(s,n,a,T,q)`. This is a family of
counterexamples for every `c,q,T>0`, not a numerical example offered
as a model of actual stem-cell turnover.

**Proof, hidden-stage counterexample.** Let all birth and death rates be
zero and let two successive replacement transitions have the same
arbitrary rate `mu>0`. One current cell already in the second stage
has progression risk `1-exp(-mu T)`. One current cell in the first
stage has risk `1-exp(-mu T)(1+mu T)`, by (12). Their observed founder
clone count, zero growth, and chronological age can be identical.
Their risk difference is `mu T exp(-mu T)>0`. For arbitrary `n>=1`
replace each single-cell survival by its `n`th power; the distinction
persists. This proves nonidentifiability in the larger model class.

If mutation production is proportional to cell division, its intensity
depends on a birth/turnover rate as well as net `b-d`. If successful
completion depends on genotype-specific additional hits or establishment,
its coefficient also varies independently of `s`. These are mechanistic
interpretations of the unidentified parameters, not measured explanations
of a particular gene's residual.

### Theorem 7 — Sharp exposure-ratio bounds and the exact transport condition

**Assumptions:** a fixed finite `T>0`; the rare-initiation cumulative-hazard
limit of Theorem 3; constant future net rates `s_g`; and a finite nonzero
nonnegative completion/opportunity measure `W` on `[0,T]` shared between
the compared genes. Matching age, initial `n`, and initiation coefficient
is required for a fitness-only risk comparison (A7). Define

$$
K_W(s)=\int_{[0,T]}e^{st}\,W(dt).
\tag{16}
$$

If `s_h>=s_l`, then

$$
1\le\frac{K_W(s_h)}{K_W(s_l)}
\le e^{(s_h-s_l)T}.
\tag{17}
$$

Both bounds are sharp over the stated measure class: mass at zero
attains the lower bound and mass at `T` attains the upper bound.
For absolutely continuous weights they are sharp suprema/infima
approached by concentrating density near the endpoints.

More generally, let the leading cumulative hazard for gene `g` be
`L_g=n_g theta_g K_(W_g)(s_g)`, where `theta_g>0` is its effective
initiation coefficient, `n_g>0`, and `W_g` is a finite nonzero nonnegative
completion measure. Thus all ratio denominators below are positive. Relative
to a proposed common measure `W`, all gene ratios are fitness-only,
`L_g/L_r=K_W(s_g)/K_W(s_r)`, **if and only if**

$$
\frac{n_g\theta_g K_{W_g}(s_g)}{K_W(s_g)}
\quad\text{is the same positive constant for every compared gene.}
\tag{18}
$$

If the completion measure really is common, condition (18) reduces
exactly to `n_g theta_g` being common. Thus matching VAF alone is
insufficient if stem-cell pool size or the effective completion
coefficient differs. Compensating nuisance factors can satisfy (18)
without equal mechanisms; equality of all mechanisms is a sufficient,
not a generally necessary, biological claim.

For a genuinely matched cumulative-hazard ratio `R=L_h/L_l`, (17)
implies the residual coefficient bound

$$
\frac{\theta_h}{\theta_l}\ge
R e^{-(s_h-s_l)T}\quad\text{when }n_h=n_l\text{ and }W_h=W_l.
\tag{19}
$$

It is invalid to substitute an unmatched marginal odds ratio or a Cox
hazard ratio into (19).

**Proof.** Let
`P_l(dt)=exp(s_l t) W(dt)/K_W(s_l)`, a probability measure on `[0,T]`.
Then `K_W(s_h)/K_W(s_l)=E_(P_l) exp((s_h-s_l)t)`. The integrand is
between one and `exp((s_h-s_l)T)`, proving (17), and endpoint measures
prove sharpness. Multiplying the ratio identities by their positive
denominators gives (18) in both directions; with common `W` the
integrals cancel. Finally
`R=(theta_h/theta_l) K_W(s_h)/K_W(s_l)` under matched `n`, and the
upper bound in (17) gives (19).

There is also an exact uniqueness result about the **whole response
curve**, which does not hold for a finite set of gene points: for finite
measures `W_1,W_2` on `[0,T]` and a constant `c>0`,

$$
K_{W_1}(s)=cK_{W_2}(s)\text{ on a nonempty open interval of }s
\quad\Longleftrightarrow\quad W_1=cW_2.
\tag{20}
$$

**Proof of (20).** Compact support makes both Laplace transforms entire
functions: their power series converge uniformly on every bounded
complex set, by domination by `exp(|s|T)`. Equality on an open real
interval therefore gives equality everywhere by analytic uniqueness.
Differentiating at zero yields equality of every moment of `W_1` and
`cW_2`. Polynomials are uniformly dense in continuous functions on a
compact interval, so the difference of these finite measures integrates
every continuous function to zero and is the zero measure. The reverse
implication follows immediately from linearity of the integral. A
finite list of gene fitnesses supplies finitely many integral
constraints and does not establish (20)'s open-interval premise.

### Theorem 8 — Why marginal outcomes do not identify the conditional model

**Assumptions:** an individual-level survival law
`S(T|n,g)=exp(-n k_g(T))` with `k_g(T)>=0`, a distribution of nonnegative
clone counts `n`, and a gene group selected independently of future
events given its baseline state. When discussing attenuation, the
distribution of `n` is the **same** for both genes and has finite
positive mean. The hazard identity is asserted at times where `k_g` is
differentiable and `k_g(T)>0`; it also holds at `k_g(T)=0` when the
baseline mean is finite. No proportional-hazards assumption is imposed.

The marginal survival is the Laplace mixture

$$
\bar S_g(T)=\mathbb E_g[e^{-n k_g(T)}],\qquad
\bar h_g(T)=k_g'(T)
\frac{\mathbb E_g[n e^{-n k_g(T)}]}{
      \mathbb E_g[e^{-n k_g(T)}]}.
\tag{21}
$$

Thus a marginal hazard ratio generally depends on the clone-count
distribution and its changing event-free risk set, not just fitness.
Even if the distribution is identical between genes and
`k_h>=k_l>0`,

$$
1\le\frac{-\log\bar S_h}{-\log\bar S_l}
\le\frac{k_h}{k_l}.
\tag{22}
$$

The upper bound is strict when the clone-count distribution is
nondegenerate and `k_h>k_l`. Different baseline distributions can
even reverse an apparent fitness ordering.

**Proof.** Conditional survival averaged over baseline counts gives
the first identity; differentiating its logarithm gives the second,
with integrability at positive `k` supplied by exponential tilting.
Put `f(k)=-log E exp(-kn)`. For `k>0`,
`f'(k)=E_k[n]>=0` and `f''(k)=-Var_k(n)<=0`; the tilted moments exist
because `x^j exp(-kx)` is bounded for `k>0`. Also `f(0)=0`.
Monotonicity gives the lower ratio bound. Concavity and `f(0)=0`
give `f(k_h)/k_h<=f(k_l)/k_l`, hence the upper bound. A nondegenerate
count distribution remains nondegenerate after a positive exponential
tilt, so strict concavity gives strictness. For reversal, choose
deterministic initial counts satisfying `n_h/n_l<k_l/k_h`; then
`n_h k_h<n_l k_l` even though `k_h>k_l`.

For a specified fixed-horizon cumulative hazard `H`, the corresponding
odds are `exp(H)-1`. Therefore an odds ratio is
`(exp(H_g)-1)/(exp(H_r)-1)`, not `H_g/H_r`; they agree only in a
specified joint rare-risk limit. A fitted Cox coefficient is a
model-dependent summary over risk sets; it cannot generally be
replaced by either a terminal-time hazard ratio or a cumulative
hazard ratio. Matching only a broad VAF bin does not match `n` or
its distribution. These distinctions justify an **untestable**
verdict for a conditional mechanistic claim when only unmatched
published marginal effect estimates are available.

### Theorem 9 — Density limitation changes the VAF law

**Assumptions:** this is an explicitly different model replacing A1–A2.
The founder clone's stem-cell fraction `x(t)` is deterministic and obeys
`x'=s x(1-x)` in a fixed pool `N`, with baseline `x in (0,1]` and finite
real `s`. Terminal marks are a Poisson process of intensity `q N x(t)`
with constant `q>=0`, independent of the deterministic path. There are
no unobserved advanced pathways or competing events. A4 is needed to
identify `x=2v`. The value at `x=0` below is a normalized limiting
exposure, not the risk of a clone with no cells.

The integrated opportunity normalized by baseline `n=Nx` is

$$
B(s,x,T)=\int_0^T\frac{x(t)}xdt=
\begin{cases}
\displaystyle\frac{\log[1+x(e^{sT}-1)]}{xs}&s\ne0,\\
T&s=0,
\end{cases}
\quad H(T)=qNxB(s,x,T).
\tag{23}
$$

For `s>=0`, `T<=B<=K(s,T)`; for `s<=0`, `K(s,T)<=B<=T`.
The small-clone limit is `lim_(x->0) B=K`. For `s>0`, `B` decreases
strictly with `x` when `T>0`; cumulative hazard is increasing and
strictly concave in baseline `x` when `q>0` and `T>0`, so the exact linear-VAF
scaling of Theorem 1 no longer holds. For matched `x`,
`s_h>=s_l` still implies

$$
1\le B(s_h,x,T)/B(s_l,x,T)
\le e^{(s_h-s_l)T}.
\tag{24}
$$

**Proof.** Solving the logistic equation gives
`x(t)=x exp(st)/(1-x+x exp(st))`. Its normalized integral is (23),
by differentiating `log(1-x+x exp(st))/s`; continuity gives `s=0`.
For `s>=0`, the normalized fraction is between one and `exp(st)`,
giving `T<=B<=K`; the inequalities reverse for `s<=0`. Dominated
convergence on the finite interval gives the `x->0` limit. For
`s>0` the derivative of the normalized fraction with respect to `x`
is `-exp(st)(exp(st)-1)/(1+x(exp(st)-1))^2`, strictly negative for
`t>0`, proving the first strictness claim. The first derivative of
`H` with respect to `x` is `qN(expm1(sT))/(s(1+x expm1(sT)))>0`;
its second derivative is
`-qN(expm1(sT))^2/(s(1+x expm1(sT))^2)<0` for the stated positive
`s,q,T`. Finally the ratio of the two normalized logistic paths is

$$
e^{(s_h-s_l)t}\frac{1-x+x e^{s_l t}}{1-x+x e^{s_h t}},
$$

which lies between one and `exp((s_h-s_l)t)` because the logistic
solution is increasing in `s` and the second factor is at most one.
The ratio of their integrals is a weighted average of this pointwise
ratio, proving (24).

A slope of `logit(2*VAF)` fits this logistic cell-fraction model, not
an unlimited exponential clone. A slope of `logit(VAF)` would be a
different quantity. An exponential approximation to logistic growth
requires a small clone throughout the horizon; saturation sensitivity
must therefore be reported for larger VAFs. A gene-average slope and
an outcome cohort with unreported baseline VAF distributions cannot
by themselves establish the conditional calibration in (23).
