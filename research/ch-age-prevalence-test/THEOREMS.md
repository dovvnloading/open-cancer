# Conditional mathematics for the CH age-prevalence audit

**Nonclinical, unreviewed research. No medical advice.**

These results concern specified mathematical observation models. They do not
establish a biological mechanism, a calibrated clinical predictor, or historical
priority. The linear branching ingredients are classical. The first result is
an integral reduction of the density used by [Watson et al. 2020,
Science](https://doi.org/10.1126/science.aay9333), with the unchanged published
parameters and precise author-notebook locators in
[data/model_parameters.csv](data/model_parameters.csv). Its density is an
approximation to the biological process. The exact integer branching results
below are separate: they do not prove the approximation or its VAF bridge.

## Notation and observation assumptions

Calendar age is `t>=0`. For the transported Watson density, `A=N*tau>0`,
`mu>=0` is mutation supply per cell per year, and the fitness density has
shape `beta>0`, scale `delta>0` and upper support `S>0`. These DFE parameters
are distinct from a birth or death rate. Write

\[
 Z=\int_0^S e^{-(s/\delta)^\beta}\,ds,\qquad
 w(s)=Z^{-1}e^{-(s/\delta)^\beta},\quad 0\le s\le S.
\]

The source prediction uses `s0=0.0001` as its lower integration limit while
normalizing over `0..S`; the omitted mass is not redistributed. In the VAF
coordinate `0<v<1/2`, set

\[
x(v)=\frac{v}{1-2v},\qquad
\phi(s,t)=\frac{e^{st}-1}{2As},\quad
\phi(0,t)=\frac{t}{2A}.
\]

The transported intensity density, conditional on its population and mutation
assumptions, is

\[
 \rho(v,s;t)=2A\mu\,w(s)
 \frac{\exp[-x(v)/\phi(s,t)]}{v(1-2v)}.
 \tag{1}
\]

Its finite-fraction coordinate is a model convention. Copy neutrality,
heterozygosity, common lineage contribution, constant fitness and supply, and
the absence of competition and survival selection remain transport assumptions.
Perfect detection is an additional observation assumption, not a consequence
of a nominal VAF cutoff or a cohort's mean sequencing depth.

For the separate integer process, `Z_BD(u)` starts with one cell and has
per-cell birth rate `b>=0`, death rate `d>=0`, net rate `s=b-d`, and integer
threshold `m>=1`. A clone's elapsed time `u` is distinct from chronological
age `t`. If an independent Poisson arrival process is used, it starts at birth
with no pre-existing clones. All rates are fixed over time within these results.

## Theorem 1 — Normalization and the E1 reduction of the Watson density

For the positive DFE parameters above,

\[
 Z=\frac{\delta}{\beta}\,
 \gamma\left(\frac1\beta,(S/\delta)^\beta\right),
 \tag{2}
\]

where `gamma(a,z)=integral_0^z y^(a-1) exp(-y) dy` is the lower incomplete
gamma function. For `S>s0`, `t>0` and `0<v_l<v_u<=1/2`, integrating (1) over
the retained fitness support gives

\[
 \Lambda(t;v_l,v_u)=2A\mu\int_{s_0}^{S}w(s)
 \left[E_1\left(\frac{x(v_l)}{\phi(s,t)}\right)
       -E_1\left(\frac{x(v_u)}{\phi(s,t)}\right)\right]ds,
 \tag{3}
\]

with the convention `E1(x(1/2)/phi)=0`, and

\[
 E_1(z)=\int_z^\infty e^{-y}\frac{dy}{y}.
\]

For `S<=s0` the retained fitness domain is empty up to a set of measure zero;
define `Lambda=0`, rather than interpreting the displayed fitness integral
as a reversed signed integral. At `t=0`, also define `Lambda=0`. For fixed
cutoffs the intensity is nondecreasing
in age; it is nonincreasing as the lower VAF cutoff increases.

**Proof.** Substituting `y=(s/delta)^beta` in the integral defining `Z`
gives (2). Also `dx/dv=(1-2v)^(-2)` and hence
`dv/[v(1-2v)]=dx/x`. Substituting `y=x/phi` then gives (3). All integrands
are nonnegative, so the order of integration is justified by Tonelli's
theorem. As `t` increases, `phi` increases; therefore the original integrand
(1) increases pointwise. As the lower cutoff increases, a nonnegative part
of the integration domain is removed. The zero-age limit follows by dominated
convergence, for example by domination with the finite integral at any fixed
positive age. This proves the statements. **QED.**

The reduction is exact **given density (1)**. It does not show that (1) is an
exact finite-population distribution or that its parameters transport between
different variant panels and assays. The pooled ten-gene DFE is its own source
fit; equation (3) does not make it the sum of four separately fitted gene DFEs.

## Theorem 2 — Exact cross-sectional tail of an integer birth–death clone

For `b>0`, `s=b-d !=0` and `u>=0`, define

\[
 q(u)=\frac{b(1-e^{-su})}{b-de^{-su}},\qquad
 a(u)=\frac{s}{b-de^{-su}}.
\]

Then

\[
 \Pr\{Z_{BD}(u)\ge m\}=a(u)q(u)^{m-1}.
 \tag{4}
\]

At criticality `b=d=c>0`,

\[
 \Pr\{Z_{BD}(u)\ge m\}
   =\frac{(cu)^{m-1}}{(1+cu)^m}.
 \tag{5}
\]

Use `0^0=1` only for the `m=1,u=0` boundary. If `b=0`, the probability is
`exp(-d*u)` for `m=1` and zero for `m>1`. This includes the lossless case
`b=d=0`.

**Proof.** The generating function `G(u,z)=E[z^Z_BD(u)]`, for `0<=z<=1`,
satisfies the backward equation

\[
 \partial_uG=bG^2-(b+d)G+d,\qquad G(0,z)=z.
 \tag{6}
\]

Solving this Riccati equation, or substituting into it and using uniqueness,
gives

\[
 G(u,z)=p_0(u)+(1-p_0(u))\frac{(1-q(u))z}{1-q(u)z},
\quad
 p_0(u)=\frac{d(1-e^{-su})}{b-de^{-su}},\quad 1-p_0(u)=a(u).
\]

The geometric-series coefficients are
`P(Z_BD(u)=n)=a(u)(1-q(u))q(u)^(n-1)` for `n>=1`. Summing them from `m`
to infinity proves (4). At `s=0`, the solution of (6) has
`q=cu/(1+cu)` and `a=1/(1+cu)`, which proves (5). If `b=0`, at most the
single original cell remains; its lifetime is exponential with rate `d`.
The stated zero-age and zero-rate limits follow directly. **QED.**

This is a probability of detection **at the observation age**. It is not the
probability of having crossed a threshold at any previous time. The latter
would count clones that subsequently shrink or become extinct.

If a fixed effective cell denominator `N` is posited, the VAF approximation
gives `m=ceil(2*N*v_min)`. That connection additionally requires a heterozygous
copy-neutral variant and unbiased output. An unlimited branching process can
exceed `N`, so this interpretation is not an exact finite-pool theorem. It is
also distinct from the coordinate convention in (1).

## Theorem 3 — An exact integrated threshold kernel and a series error bound

Let

\[
 K_m(t)=\int_0^t\Pr\{Z_{BD}(u)\ge m\}\,du.
\]

For `b>0`, with `q(t)` from Theorem 2 (including its critical definition),

\[
 K_m(t)=\frac1b\left[-\log(1-q(t))-
                   \sum_{k=1}^{m-1}\frac{q(t)^k}{k}\right]
       =\frac1b\sum_{k=m}^{\infty}\frac{q(t)^k}{k}.
 \tag{7}
\]

For `M>=m`, truncating the positive series gives the rigorous remainder bound

\[
 0\le K_m(t)-\frac1b\sum_{k=m}^{M}\frac{q^k}{k}
       \le\frac{q^{M+1}}{b(M+1)(1-q)}.
 \tag{8}
\]

For `b=0,m=1`, the integral is `(1-exp(-d*t))/d` if `d>0`, and `t` if
`d=0`. For `b=0,m>1` it is zero.

**Proof.** Differentiating the expression for `q`, including its continuous
critical limit, gives

\[
 q'(u)=(1-q(u))(b-dq(u)),\qquad
 a(u)=(b-dq(u))/b.
\]

Consequently `a(u)q(u)^(m-1) du=q^(m-1)dq/[b(1-q)]`. Integration from zero
to `q(t)` proves the first expression in (7). Expanding `1/(1-q)` proves the
second. In the remaining sum after `M`, replacing each denominator `k` by
`M+1` and summing the resulting geometric series proves (8). The pure-death
cases are elementary integrals of Theorem 2. **QED.**

For `s=b-d>0` and fixed `m`, (7) also gives

\[
 K_m(t)=\frac{s}{b}t+
        \frac{\log(b/s)-H_{m-1}}b+o(1),\qquad t\longrightarrow\infty,
 \tag{9}
\]

where `H_0=0` and `H_j=sum_(k=1)^j 1/k`. Indeed,
`1-q(t)~(s/b)exp(-s*t)` and the finite sum tends to `H_(m-1)`. This is a
mathematical late-time statement, not a justified extrapolation of unlimited
clones in a finite biological population.

## Theorem 4 — Poisson arrivals, detection, and prevalence

Suppose new single-cell clones arrive independently at a constant rate
`lambda>=0` from age zero, each evolves independently according to Theorem 2,
and no clones exist at birth. The number detectable at age `t` is Poisson with
mean `lambda*K_m(t)`. Thus

\[
 p(t)=1-e^{-\lambda K_m(t)}.
 \tag{10}
\]

More generally, a cohort-specific independent detection rule `D(z) in [0,1]`
replaces `K_m(t)` by `integral_0^t E[D(Z_BD(u))] du`. Independent clone
classes add their intensities before the Poisson zero probability is taken.
In particular, if the transported density (1) is assumed to define a Poisson
point process of clones, its prevalence is `1-exp(-Lambda)` from (3).

**Proof.** An arrival at chronological age `v` is retained at age `t` with
probability `P(Z_BD(t-v)>=m)`. Independent thinning of an inhomogeneously
marked Poisson process leaves a Poisson count with mean
`lambda*integral_0^t P(Z_BD(t-v)>=m) dv=lambda*K_m(t)`. The probability of
at least one clone is one minus the zero probability. Replacing the retention
indicator by a random retention probability `D` gives the generalization.
Independent Poisson counts sum to a Poisson count with summed means. The last
statement applies the same zero-count identity to the separate point-process
assumption. **QED.**

This theorem requires independence at the level of relevant clone arrivals
and observations. Two driver mutations in one clone do not automatically
constitute two independent clones. Overlapping gene detections do not justify
adding person-level prevalences. The observational endpoint must match the
classes and panel represented in the intensity.

## Theorem 5 — Sharp age-band bounds and common-supply feasibility

Let `L(t)>=0` be a fixed finite, continuous, nondecreasing expected-clone
intensity. This can be
`Lambda(t)` from Theorem 1 or `lambda*K_m(t)` from Theorem 4. For a
nonnegative finite common supply multiplier `c`, set

\[
 p_c(t)=1-e^{-cL(t)}.
\]

If an age band contains unknown ages `a<=T<=b`, the sharp possible range of
its mean prevalence is

\[
 p_c(a)\le E[p_c(T)]\le p_c(b).
 \tag{11}
\]

For a reported observed interval `0<=ell<=u<1` and `0<L(a)<=L(b)`, the set
of compatible multipliers is precisely

\[
 \boxed{\quad
 \frac{-\log(1-\ell)}{L(b)}
 \le c\le
 \frac{-\log(1-u)}{L(a)}.\quad}
 \tag{12}
\]

For several bands, the sharp common-`c` set is the intersection of their
intervals. For half-open age bands these assertions concern the closure of
the prediction set; an excluded endpoint can be approached arbitrarily
closely. The midpoint or uniform-age prediction is an additional assumption.

**Proof.** Monotonicity gives (11). A mixture supported at `a` and `b` attains
every mean between the endpoints, proving sharpness. Compatibility with
`[ell,u]` is equivalent to `p_c(b)>=ell` and `p_c(a)<=u`. Solving these two
inequalities gives (12). A multiplier is jointly compatible precisely when
it satisfies every band's inequalities, hence the intersection rule. The
endpoint-mixture argument also proves sufficiency for each band separately;
no shared within-band age distribution was assumed across bands. **QED.**

Boundary conventions are substantive. If `L(a)=0`, the upper multiplier
bound is infinite. If `L(b)=0` and `ell>0`, compatibility is impossible. If
`ell=0`, the lower bound is zero. If `u=1`, there is no finite upper bound.
An interval consisting only of observed probability one is unattainable with
finite intensity and finite `c`. The software returns a reversed interval
for this impossible case. The frozen open-ended-age rule uses the conservative
upper prediction one; it does not invent a maximum lifespan.

Empty intersection is a deterministic incompatibility of the specified
model, age support and observed intervals. Rounding intervals are not sampling
confidence intervals. Marginal source confidence intervals do not become a
joint confidence region by taking their Cartesian product. No significance
level follows from this theorem alone.

## Theorem 6 — Exact supply/thinning alias and limits of fitness-only prediction

**(a) Single-age supply freedom.** For any fixed positive intensity shape
`L(t)>0` and desired `0<p<1`, the unique multiplier giving that prevalence is

\[
 c=\frac{-\log(1-p)}{L(t)}.
 \tag{13}
\]

Consequently a single age-prevalence point cannot distinguish competing
fitness kernels if each is permitted an unrestricted supply multiplier.

**(b) An all-age assay alias.** Suppose detectable clones undergo independent
age-invariant assay thinning with unknown probability `eta in (0,1]`. Then

\[
 p_{\rm obs}(t)=1-e^{-\eta c L(t)}.
 \tag{14}
\]

Any two admissible pairs `(eta,c)` having the same product give the identical
entire age-prevalence curve. For example, replacing them by
`(eta/k,k*c)` for `k>=1` preserves every prediction. Thus even exact
all-age prevalence does not identify the supply multiplier and constant
per-clone sensitivity separately.

**(c) Net birth–death fitness is insufficient.** Fix net fitness `s=0`.
The models `b=d=0` and `b=d=r>0` have the same net fitness. For `m>=2`,
the first gives zero threshold probability at every time, while the second
gives `(r*u)^(m-1)/(1+r*u)^m>0` for every `u>0`. For `m=1` their survival
probabilities are one and `1/(1+r*u)`, respectively. Hence net fitness alone
does not determine the threshold probability or its integrated kernel.

**Proof.** Equation (13) is the unique inverse of the strictly increasing
map `c -> 1-exp(-cL(t))`. Independent thinning multiplies Poisson intensity
by `eta`, proving (14) and its alias. Part (c) follows directly from the
zero-rate and critical cases of Theorem 2. **QED.**

Multiple ages can still constrain a common supply through Theorem 5. The
alias in part (b) does not say that every age-shape restriction is untestable,
nor that a complete noise-free curve under stronger assumptions could never
identify other parameters. A suitably matched published aggregate assay
calibration, such as known-reference-clone detections and denominators at the
relevant VAFs and depths, would constrain `eta`; an independent mutation-supply
estimate would constrain `c`. A mean sequencing depth supplies neither.

## Theorem 7 — Necessary coherence of gene and pooled carrier predictions

Fix a common population, age condition, VAF threshold, variant inclusion rule
and observation rule. Let `A_g` be the event that a person has at least one
detected variant in gene `g` under that rule, and let `A_pool` be the event
that a person has at least one detected variant in a panel containing all
variants included in `A_g`. If `A_g` is a subset of `A_pool`, then

\[
 \Pr(A_g)\le\Pr(A_{\rm pool}).
 \tag{15}
\]

No independence of genes, clones, mutations, or persons is required. If both
predictions are written as Poisson complements with finite nonnegative
intensities, the same condition is equivalent to `Lambda_g<=Lambda_pool`.

**Proof.** Event inclusion implies the pointwise indicator inequality
`1_(A_g)<=1_(A_pool)`. Taking expectations proves (15). Since
`Lambda -> 1-exp(-Lambda)` is strictly increasing on nonnegative finite
intensities, its order is equivalent to the order of those intensities.
**QED.**

This is a necessary condition for interpreting several predictions jointly
as exact carrier probabilities in the same population. A violation rejects
that joint interpretation of the assembled predictions. It does not identify
which fitted forecast is inaccurate, provide sampling uncertainty for either
fit, or reject every fitness estimate or branching model. Different age
supports, populations, panels or detection rules need not define nested
events; the inequality must not be imposed across such unmatched endpoints.

Here the separately fitted gene DFEs and the separately fitted pooled DFE
are retained unchanged. They were not constructed by adding gene intensities,
and the audit does not attribute a claim of an exact coherent joint model to
the source authors. The numerical check compares their candidate forecasts
only where the prediction grid gives the same age and VAF and the mutation
supplies identify the gene as belonging to the pooled panel. It is an
explicitly post-scoring structural audit, separate from empirical cohort
scoring; it does not refit or select a DFE after seeing an outcome.

## Numerical verification and its limits

[model.py](model.py) evaluates `E1` by a convergent series or continued
fraction, the DFE normalization by incomplete-gamma methods, and the remaining
one-dimensional integral over log fitness by successive Gauss–Legendre rules.
This avoids the very small scale parameters and endpoint derivative
singularities that make a naive rectangle sum unreliable. It does not copy
the author's implementation. The frozen normalization and the lower-fitness
cutoff are preserved exactly as specified in the protocol.

Independent checks use the transformed integral

\[
 E_1(x)=\int_0^\infty \exp[-x e^y]dy,
\]

composite Simpson quadrature, a smooth power substitution for DFE
normalization, and a separate incomplete-gamma coordinate for the prevalence
integral. The tail beyond `y=log((x+R)/x)` is at most
`exp(-(x+R))/(x+R)`. Integer-clone checks can independently use the generating
equation (6) or a stopped finite population generator with its lost/overflow
mass retained as an error bound. Formula (8) provides an analytic series
remainder bound. Floating-point arithmetic itself is not outward-rounded
interval arithmetic.

The binomial read illustration is
`P[Binomial(D,v)>=r]` for a declared depth `D` and alternate-read threshold
`r`. It assumes independent reads and omits the other caller filters. It is
not an estimate of an actual cohort's full sensitivity, and it cannot repair
a variant-panel mismatch. Biological validity, source accuracy, source reuse
rights and cohort independence require the separate evidence audit; finite
computational checks prove none of them.
