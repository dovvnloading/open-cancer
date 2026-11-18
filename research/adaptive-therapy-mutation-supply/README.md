# Mutation supply does not imply a universal MTD cutoff

**Nonclinical computational mathematics. Synthetic, two-cell model; no treatment
advice, patient data, experimental validation, or claim of clinical efficacy.**
Research status: **unreviewed**. The proofs and executable checks below concern
the specified mathematical process. They do not establish biological validity
or bibliographic priority.

## Result and explicit pivot

The broad mutation-supply versus competitive-suppression question is already
studied, including stochastic models. It would be incorrect to claim that
adding de novo resistance to containment theory is new. We instead give an
exact, narrowly scoped counterexample to a **single mutation-rate cutoff** for
an explicitly specified progression-time objective.

The model starts with **zero resistant cells**, preserves sensitive extinction,
uses **finite-rate sensitive clearance** in an aggressive-treatment arm, and
defines progression as **total burden strictly exceeding the containment
threshold**. An exact cubic determines the comparison for an entire family of
rates. In a particular instance, containment has longer expected control time
at low and high mutation rates, while the aggressive arm wins only in an
intermediate interval. The two reversals persist under sufficiently small
parameter perturbations.

There are two essential qualifications:

1. The literal, uncensored expected time to progression is **infinite under
   both arms**, because each has positive probability of never progressing.
   Thus there is no strict mutation threshold for that original expectation
   in this model. The finite comparison is
   `E[min(T, H)]`, where `H` is an independent exponential observation horizon.
   It is neither a conditional mean nor a fixed-time-horizon mean.
2. Containment is an **ideal instantaneous population controller**. The
   aggressive arm is a **finite-clearance MTD surrogate**, not a proven
   pharmacological maximum within a common bounded-dose control system.
   We use `C` and `M` for these two precisely defined mathematical policies.

No separate branching-establishment approximation is made: this is a small,
exact interacting birth–death–conversion chain. The model is deliberately
chosen to test a universal mathematical inference, not to estimate a tumour's
mutation threshold.

## 1. Literature audit and novelty assessment

Primary-source audit completed **2026-11-17**. Searches covered containment,
de novo resistance, mutation supply, stochastic establishment, finite killing,
nonmonotonic mutation comparisons, multiple thresholds, and discounted or
restricted progression times. This is a targeted audit, not an exhaustive
systematic review. Paper links are references, not mirrored source data.

| Primary work | What is already known | Boundary relevant here |
| --- | --- | --- |
| Hansen, Woods & Read (2017), [*How to Use a Chemotherapeutic Agent When Resistance to It Threatens the Patient*](https://doi.org/10.1371/journal.pbio.2001110) | Deterministic analysis explicitly balances new resistance arising from sensitive cells against their competitive suppression of resistance; derives balance thresholds and failure regimes. | Aggressive removal is idealized as immediate. The broad trade-off is not our contribution. |
| Hansen & Read (2020), [*Cancer therapy: Attempt cure or manage drug resistance?*](https://doi.org/10.1111/eva.12994) | Birth–death containment with ongoing mutation, exact cure probabilities and progression times, and explicit dependence on how cure and delay are valued. | Aggressive treatment immediately removes sensitive cells. Their illustrative containment chain stops at resistant burden `B`; ours stops at total burden `> B`, and follows the all-resistant state until progression or extinction. |
| Viossat & Noble (2021), [*A theoretical analysis of tumour containment*](https://doi.org/10.1038/s41559-021-01428-w) | General deterministic containment comparisons, including settings without a fitness cost of resistance. | The core general comparison omits mutation; the work is not a license to claim that the whole literature omits it. |
| Alvarez & Viossat (2024), [*Tumor containment: a more general mathematical analysis*](https://doi.org/10.1007/s00285-024-02062-3) | Extends deterministic containment comparisons to mutation-compatible settings and Norton–Simon treatment effects. | The comparison requires net resistant growth to decrease with sensitive abundance. The authors' [2022 preprint](https://arxiv.org/html/2211.17041v1) explicitly assumes positive initial resistance; it is not this zero-origin stochastic problem. The preprint is not represented as textually identical to the final publication. |
| Iwasa, Nowak & Michor (2006), [*Evolution of resistance during clonal expansion*](https://doi.org/10.1534/genetics.105.049791); Foo & Michor (2009), [*Evolution of Resistance to Targeted Anti-Cancer Therapies during Continuous and Pulsed Administration Strategies*](https://doi.org/10.1371/journal.pcbi.1000557); Bozic et al. (2013), [*Evolutionary dynamics of cancer in response to targeted combination therapy*](https://doi.org/10.7554/eLife.00747) | Branching-process resistance theory includes both resistance accumulated before treatment and resistance generated during treatment. In particular, Foo–Michor and Bozic et al. do not restrict resistance to a fixed initial population. | Branching-process resistance and de novo supply are established methods. We do not claim to introduce either. |
| Kuosmanen et al. (2021), [*Drug-induced resistance evolution necessitates less aggressive treatment*](https://doi.org/10.1371/journal.pcbi.1009418) | Inhomogeneous mutation arrivals, branching establishment, extinction probabilities, dose-dependent mutation, and stochastic simulation. | Its establishment approximation assumes a constant establishment probability. Neither Poisson thinning nor the existence of nonmonotonic treatment effects is novel here. |
| Wang, Scott & Vladimirsky (2024), [*Threshold-awareness in adaptive cancer therapy*](https://doi.org/10.1371/journal.pcbi.1012165) | Stochastic optimal-control comparisons depend on the objective, with threshold/budget-aware formulations. | Their threshold is an outcome budget, not this algebraic mutation-rate comparison. Objective dependence itself is not claimed as new. |

Related model context includes Gatenby et al. (2009),
[*Adaptive Therapy*](https://doi.org/10.1158/0008-5472.CAN-08-3658),
Gallaher et al. (2018),
[*Spatial heterogeneity and evolutionary dynamics modulate time to recurrence
in continuous and adaptive cancer therapies*](https://doi.org/10.1158/0008-5472.CAN-17-2649),
and Strobl et al. (2021; online 2020),
[*Turnover modulates the need for a cost of resistance in adaptive therapy*](https://doi.org/10.1158/0008-5472.CAN-20-0806).
These contain much richer biological modeling than the counterexample below.
Martin et al. (1992),
[*Optimal control of tumor size used to maximize survival time when cells are
resistant to chemotherapy*](https://doi.org/10.1016/0025-5564(92)90038-X),
is an earlier optimal-control predecessor to the containment literature.

**Specific candidate contribution.** We did not locate the exact general cubic
reduction below, its high-mutation sign criterion, or the corrected
finite-clearance, zero-initial-resistance example with two ranking reversals
and no sensitive-cell resurrection in the primary sources inspected.
The contribution claimed is this explicit theorem and reproducible
counterexample under its assumptions. The broad trade-off, CTMC backward
equations, discounted resolvents, cure atoms, and objective dependence are
already known. **Priority remains unestablished.** No claim of novelty is based
on AI agreement or on passing tests.

## 2. Model and all assumptions

The following assumptions apply to Theorems 1 and 2 unless explicitly
specialized. They are part of the results, not optional interpretation.

**A1 — Population and origin.** `S` and `R` are integer counts of individual
sensitive and resistant cells. Both policies start at `(S,R)=(2,0)`.
The containment threshold is `B=2`. These are two actual population units in
the model, not two macroscopic tumour compartments or a fitted cell-number
rescaling.

**A2 — Resistance generation.** Each sensitive cell converts irreversibly to
one resistant cell at rate `mu >= 0`, with transition
`(s,r) -> (s-1,r+1)` at rate `mu*s`. This conserves total burden. `mu` is an
inverse-time conversion rate, **not a mutation probability per division**.
It can represent a coarse per-time acquisition of resistance, but its
equivalence to genomic mutation during cell division is not established.
It acts during treatment even when net sensitive reproduction is zero.
There is no backmutation or resistant immigration.

**A3 — Competition and resistant demography.** At total burden `N=1`, each
resistant cell gives birth at rate `a`; at `N=2`, it gives birth at rate `b`,
where `a >= b > 0`. Each resistant cell dies at rate `d > 0`. These rates are
identical under the two policies and are not directly altered by drug.
Conditioned on the current state, all listed events have exponential clocks.
The chain is stopped at progression, so no birth law is needed for `N>=3`.

**A4 — Aggressive policy M.** Sensitive cells have no net reproduction and
are removed at finite rate `k*s`, `k>0`. The conversion clock in A2 remains
active. Resistant births and deaths are uncompensated. This is our
finite-clearance MTD abstraction. No pharmacokinetic or dose-response law
identifies `k` with a measured maximum tolerated dose.

**A5 — Ideal policy C while S is positive.** Total burden is maintained at
two through simultaneous composite transitions. At `(1,1)`, an R birth
removes the surviving S, yielding `(0,2)`; an R death is accompanied by a
division of the surviving S, yielding `(2,0)`. The compensating division is
stipulated sensitive and has no additional mutation draw; A2 is the only
conversion mechanism. There is no other net S birth/death event. Conversion
itself conserves total burden. These instantaneous corrections are a
singular idealization, not a demonstrated bounded-dose treatment.

**A6 — Sensitive extinction is permanent.** Once `S=0`, it stays zero under
both policies. In particular, an R death at `(0,2)` leads to `(0,1)`, **not**
`(1,1)`. Pure-R states have the same dynamics under C and M; total extinction
`(0,0)` is absorbing. C therefore holds total burden at two only while a
sensitive cell remains.

**A7 — Progression endpoint.** `T` is the first completed, post-control jump
with `S+R=3`. It is **not** the first arrival at `R=2`, the first mutant, or
the first established resistant lineage. The simultaneous C compensation in
A5 does not visit a separate intermediate progression state. If progression
never occurs, including after total extinction, set `T=infinity`.

**A8 — Finite objective.** Independently of the population process let
`H_z ~ Exp(z)`, with `z>0`. Define

$$
h_X(\mu,z)=\mathbb E_X[e^{-zT}],\qquad
U_X(\mu,z)=\mathbb E_X[\min(T,H_z)]
          =\frac{1-h_X(\mu,z)}{z},\quad X\in\{C,M\},
$$

where `exp(-z*infinity)=0`. A larger `U` means longer expected control time
for this objective. `H_z` is a mathematical observation horizon, not a
patient's lifetime or a fixed follow-up time. The ordinary expectation
`E[T]` is analyzed separately.

### Complete transition specification

All unlisted transition rates are zero. Identical destinations have their
rates added. A cure/extinction state has value zero for progression
probabilities and Laplace transforms.

| Policy/state | Destination | Rate |
| --- | --- | --- |
| C, `(2,0)` | `(1,1)` | `2*mu` |
| C, `(1,1)` | `(0,2)` | `mu+b` (conversion plus compensated R birth) |
| C, `(1,1)` | `(2,0)` | `d` (R death plus resident S division) |
| M, any `s>0` | `(s-1,r+1)` | `mu*s` |
| M, any `s>0` | `(s-1,r)` | `k*s` |
| M, any `r>0` | `(s,r+1)` or progression if `s+r=2` | `r*a` if `s+r=1`, otherwise `r*b` |
| M, any `r>0` | `(s,r-1)` | `d*r` |
| Either, `(0,1)` | `(0,2)` | `a` |
| Either, `(0,1)` | `(0,0)` | `d` |
| Either, `(0,2)` | progression | `2*b` |
| Either, `(0,2)` | `(0,1)` | `2*d` |

The final four rows specify the shared pure-R boundary and do not add duplicate
events to the M rows. C has five pre-progression states including extinction;
M has six. Both are nonexplosive finite-state chains for every finite `mu`.

## 3. Theorem 1 — General cubic criterion and the uncensored obstruction

**Theorem 1.** Under A1–A8, with fixed finite rates
`a>=b>0`, `d,k,z>0`, put

$$
q=\frac{a}{z+a+d},\qquad
h=\frac{2b}{z+2b+2d-2dq},
$$

$$
K_1=z+k,\quad K_2=z+k+b+d,\quad K_3=z+2k,\quad
D=3z+2b,\quad E=z(z+b+d),
$$

$$
V_1=K_1+\frac b h+(2k+d)q,\qquad
V_0=\frac{K_1b}{h}+kq(K_1+K_2).
$$

For `mu>=0`,

$$
h_C=\frac{2\mu h(\mu+b)}{2\mu^2+D\mu+E},\qquad
h_M=\frac{2\mu h(\mu^2+V_1\mu+V_0)}
 { (\mu+K_1)(\mu+K_2)(2\mu+K_3)}.                         \tag{1}
$$

Define the polynomial

$$
\Psi(\mu)=(\mu+b)(\mu+K_1)(\mu+K_2)(2\mu+K_3)
 -(\mu^2+V_1\mu+V_0)(2\mu^2+D\mu+E).              \tag{2}
$$

The fourth-degree terms cancel, so `Psi` has degree at most three. For every
`mu>0`, the following is a sharp, explicitly computable criterion:

$$
\boxed{\ U_M>U_C\ \Longleftrightarrow\ \Psi(\mu)>0\ } . \tag{3}
$$

Equality and the reverse inequality correspond exactly to `Psi=0` and
`Psi<0`. Its cubic coefficient and the high-mutation comparison are

$$
\ell=4k(1-q)-z
     =\frac{4k(z+d)}{z+a+d}-z,\qquad
\lim_{\mu\to\infty}\mu(U_M-U_C)=\frac{h\ell}{2z}. \tag{4}
$$

Consequently, if `4k(z+d)>z(z+a+d)`, M has larger `U` at all sufficiently
large finite `mu`; if the inequality is reversed, C does. On equality,
the lower-degree terms of (2) decide; no conclusion is obtained by dropping
them. This is a criterion inside this model family, not a rule for patients.

For every finite `mu>0`, **both ordinary expectations `E_C[T]` and `E_M[T]`
are infinite**. At `mu=0`, both policies have `T=infinity` almost surely,
`h_C=h_M=0`, and `U_C=U_M=1/z`.

### Proof of Theorem 1

**Step 1: the objective identity and first-step equations.** For fixed
`T=t`, independence gives

$$
\mathbb E[\min(t,H_z)]=\int_0^t e^{-zu}\,du
 =\frac{1-e^{-zt}}z,
$$

also for `t=infinity` by monotone convergence. Averaging proves A8.
Equivalently, `h_X` is the probability of progression before an independent
rate-`z` censor clock.

At a pre-progression state with outgoing rates `lambda_j` and total `lambda`,
the first-event equation is

$$
(z+\lambda)f=\sum_j\lambda_j f_j,
$$

with progression value one and extinction value zero. This finite linear
system has a unique solution: after moving nonboundary terms to the left,
each diagonal exceeds the sum of absolute off-diagonal entries by at least
`z>0`, so a maximum-absolute-component argument rules out a nonzero vector in
its kernel. Thus solving these equations identifies the actual transform,
not merely a possible algebraic solution.

**Step 2: shared pure-resistant states.** Let `f_1,f_2` be the transforms
starting at `(0,1),(0,2)`. Then

$$
(z+a+d)f_1=af_2,\qquad
(z+2b+2d)f_2=2b+2df_1.
$$

Therefore `f_1=qh`, `f_2=h`. Here `0<q<1`, and the denominator of `h` is
`z+2b+2d(1-q)>0`. In particular, every denominator subsequently used is
strictly positive.

**Step 3: containment.** Write `c_0,c_1` for the transforms at `(2,0),(1,1)`.
A5–A6 give

$$
(z+2\mu)c_0=2\mu c_1,\qquad
(z+\mu+b+d)c_1=(\mu+b)h+dc_0.
$$

Eliminating `c_1` gives

$$
c_0=\frac{2\mu(\mu+b)h}
 {(z+2\mu)(z+\mu+b+d)-2\mu d}.
$$

Expansion of the denominator gives `2mu^2+(3z+2b)mu+z(z+b+d)`, proving
the first part of (1). Notice the permanent pure-R continuation `h` rather
than an absorbing value one at `R=2`.

**Step 4: finite-rate aggressive treatment.** Write `y,t,x` for the
transforms at `(1,0),(1,1),(2,0)`. The direct equations are

$$
y=\frac{\mu qh}{\mu+K_1},\qquad
t=\frac{b+\mu h+kqh+dy}{\mu+K_2},\qquad
x=\frac{2\mu t+2ky}{2\mu+K_3}.                       \tag{5}
$$

The term `b` in the equation for `t` represents immediate **total-burden**
progression after an uncompensated R birth. In C that birth instead reaches
`(0,2)` with transform `h`, which is why the equations differ.
Substitution of `y` and then `t` into `x`, with a common denominator, yields

$$
x=\frac{2\mu h\{\mu^2+
[K_1+b/h+(2k+d)q]\mu+K_1b/h+kq(K_1+K_2)\}}
 {(\mu+K_1)(\mu+K_2)(2\mu+K_3)}.
$$

This is the second part of (1).

**Step 5: sign and high-rate criterion.** Subtract the two transforms:

$$
h_C-h_M=\frac{2\mu h\Psi(\mu)}
 {(2\mu^2+D\mu+E)(\mu+K_1)(\mu+K_2)(2\mu+K_3)}.          \tag{6}
$$

All factors outside `Psi` are positive for `mu>0`.
Since `U_M-U_C=(h_C-h_M)/z`, this proves (3) with the stated inequality
direction. The degree-four coefficients in (2) are both two. The degree-three
coefficient is

$$
K_3+2(K_1+K_2)+2b-D-2V_1
=4k+2b+2d-2b/h-(4k+2d)q
=4k(1-q)-z,
$$

where the last equality uses `2b/h=z+2b+2d-2dq`. The denominator
of (6) has leading term `4mu^5`, proving (4).

**Step 6: ordinary means and zero mutation.** Under M, the probability that
the two sensitive cells are both killed before any mutation is
`[k/(k+mu)]^2>0`. That path ends at extinction, so `P_M(T=infinity)>0`.

Under C and `mu>0`, the mean time `W` to enter `(0,2)` from `(2,0)` is finite.
First-step equations without censoring give

$$
W=\frac{3\mu+b+d}{2\mu(\mu+b)}<\infty.
$$

For example, this follows by writing `w_0=1/(2mu)+w_1` and
`w_1=[1+d w_0]/(mu+b+d)`. The embedded two-state chain has a positive
probability `(mu+b)/(mu+b+d)` of exiting on each visit to `(1,1)`, so exit
occurs almost surely and these expectations are legitimate. From `(0,2)`,
two successive deaths have probability
`[d/(b+d)]*[d/(a+d)]>0`, ending at extinction before progression.
Hence `P_C(T=infinity)>0` also. Any nonnegative extended random variable
with positive mass at infinity has infinite expectation. When `mu=0`, no R
cell can ever arise, establishing the remaining boundary claims. QED.

## 4. Theorem 2 — A bounded mutation interval where M wins

**Theorem 2.** Impose A1–A8 with

$$
(a,b,d,k,z)=(100,1,3,2,1).
$$

These are synthetic inverse-time rates; the competition ratio is deliberately
large. The transform formulas simplify to

$$
h_C=\frac{26\mu(\mu+1)}{21(2\mu^2+5\mu+5)},\qquad
h_M=\frac{\mu(26\mu^2+295\mu+626)}
 {21(\mu+3)(\mu+7)(2\mu+5)}.                      \tag{7}
$$

Let

$$
P(\mu)=18\mu^3-185\mu^2-517\mu+400.
$$

It has exactly two positive roots `alpha < beta`, with certified intervals

$$
0.637350672391423878287347192554
<\alpha<0.637350672391423878287347192555,
$$

$$
12.442615642162591815264324420886
<\beta<12.442615642162591815264324420887.
$$

The exact ranking is

| Mutation rate | Larger expected exponentially capped control time |
| --- | --- |
| `mu=0` | Tie: both equal 1 |
| `0 < mu < alpha` | C |
| `mu=alpha` | Tie |
| `alpha < mu < beta` | M |
| `mu=beta` | Tie |
| `mu > beta` | C |

Thus there is **no finite upper cutoff above which M always wins** for this
objective in this model, and no single-cutoff classification of the two arms.
This is a counterexample to a universal claim, not a claim that every model
has two reversals. The two-reversal property holds on some open neighborhood
of the specified five rates.

### Proof of Theorem 2

The pure-R transforms from Theorem 1 are `f_1=25/42`, `f_2=13/21`.
Substituting the rates into (1) gives (7). Exact cross-multiplication gives

$$
h_C-h_M=
\frac{-\mu P(\mu)}
 {21(\mu+3)(\mu+7)(2\mu+5)(2\mu^2+5\mu+5)}.       \tag{8}
$$

Equivalently, the general comparison polynomial is `Psi=-P/26`.
The denominator is strictly positive for `mu>=0`. At `z=1`,
`U_C-U_M=-(h_C-h_M)`, so its sign for `mu>0` is the sign of `P`.

In descending degree, the nonzero coefficients of `P` have signs `+,-,-,+`.
Descartes' rule gives at most two positive roots, counted with multiplicity.
But `P(0)=400>0`, `P(1)=-284<0`, and `P(13)=1960>0`; the intermediate value
theorem gives one root in `(0,1)` and another in `(1,13)`. These exhaust the
positive roots and are simple. The sign is therefore positive, negative,
positive on the three intervals. The displayed narrower intervals follow
by exact rational substitution of each decimal endpoint, with opposite signs
at the two endpoints. The verifier also bisects in exact arithmetic to width
at most `10^-45`; no numerical root finder is trusted for the root count.

For an open neighborhood statement, the coefficients of `Psi` depend
continuously on the positive rates. At the displayed rates their signs, in
descending order, are `-,+,+,-`, and their values at `mu=1/4,5/2,25` have
strict signs `-,+,-`. All these strict signs persist on a sufficiently small
open neighborhood with `a>b>0` and `d,k,z>0`. Descartes' rule and the two sign
brackets then again give exactly two simple positive roots. No particular
biologically plausible neighborhood is asserted. QED.

### Exact witnesses and boundary checks

Here positive `U_C-U_M` favors C, and negative favors M:

| `mu` | Exact `U_C-U_M` | Decimal |
| --- | --- | --- |
| `1/4` | `16606/4441437` | `0.0037388800066284855` |
| `5/2` | `-101/7524` | `-0.013423710792131845` |
| `25` | `38275/14281344` | `0.002680069886979825` |

As `mu -> infinity`, both transforms tend to `13/21`, and both `U` values
tend to `8/21`. Their difference tends to zero; its sign remains strict
above `beta` for every **finite** mutation rate. Here the high-rate
coefficient is `ell=-9/13`, so
`mu*(U_M-U_C) -> -3/14`, consistent with Theorem 1.

## 5. Corollary 2 — Eventual progression gives a different ordering

Under the rates of Theorem 2, write `p_X=P_X(T<infinity)`. For every `mu>0`,

$$
p_C=\frac{103}{112},\qquad
p_M=\frac{\mu(103\mu^2+1018\mu+1824)}
 {112(\mu+2)^2(\mu+6)},
$$

$$
p_C-p_M=
\frac{3\mu^2+265\mu+618}{28(\mu+2)^2(\mu+6)}>0.    \tag{9}
$$

Thus M has lower eventual progression probability for **every** finite
positive mutation rate, even where C has the larger `U`. Neither comparison
is a clinical recommendation. Both uncensored expected progression times
remain infinite.

**Proof.** At zero censoring rate, the pure-R first-step equations give
`p_1=(100/103)p_2` and `8p_2=2+6p_1`, hence
`p_1=25/28`, `p_2=103/112`. Absorption into progression or extinction occurs
almost surely: the finite pure-R chain has an absorbing exit reachable from
each live state. C reaches the pure-R two-cell state almost surely for
`mu>0`, by Step 6 of Theorem 1; therefore `p_C=p_2` and its extinction
probability is `9/112`.

For M, the direct equations are

$$
y=\frac{\mu p_1}{\mu+2},\qquad
t=\frac{1+\mu p_2+2p_1+3y}{\mu+6},\qquad
x=\frac{2\mu t+4y}{2\mu+4}.
$$

Substitution gives the formula for `p_M`; subtracting it from `103/112`
gives the positive expression (9). At `mu=0` the formulas requiring eventual
exit from the all-sensitive C state do not apply: **both probabilities are
zero**. In particular `p_C(0)=0` but
`lim_{mu down to 0} p_C(mu)=103/112`. This is the difference between an
infinite observation period and an arbitrarily slow positive conversion
clock, not a numerical error. QED.

## 6. Deterministic verification and independent stochastic check

All code is Python standard library only, requires Python 3.11+, reads no
external data, makes no network requests, and writes no result files.
Run from the repository root:

```sh
python research/adaptive-therapy-mutation-supply/verifier.py
python research/adaptive-therapy-mutation-supply/simulation.py --samples 300000
python -m unittest discover -s tests -p test_research_adaptive_therapy_mutation_supply.py -v
```

The standalone verifier uses `fractions.Fraction` for every scientific
calculation. `decimal.Decimal` is used only to display outward-rounded root
intervals; decimal display values are not the certificates. It:

- constructs the general comparison polynomial coefficient by coefficient
  and checks `Psi=-P/26` exactly;
- solves independently assembled occupancy-state backward equations by
  exact Gaussian elimination and compares them with (1) in 96 comparisons
  over four parameter sets and twelve mutation rates, including zero,
  very small rates and very large rates;
- verifies eventual progression probabilities, their strict ordering, and
  the zero-rate exception;
- checks the Descartes sign count and exact sign brackets, then isolates
  each positive root with rational bisection.

The test contains a **third deterministic method**: separately hand-enumerated
embedded jump chains with a censor event. It sums the probability of
progression over successive jumps, without calling the verifier's generator
or Gaussian solver. If `L_n` is the accumulated absorption probability and
`r_n` the remaining live probability, then
`L_n <= h <= L_n+r_n`. This is exact because all uncounted eventual successes
are among the paths still live. The test reduces `r_n` below `10^-25` and
checks that the closed form lies inside that rational interval for both
policies at the three ranking witnesses and at zero.

An additional **exact fixed-horizon falsification** is included in
`fixed_horizon_interval`. At `mu=5/2`, a deterministic horizon of one gives

$$
\mathbb E_C[\min(T,1)]\simeq0.8106969610908228650310043112,
\qquad
\mathbb E_M[\min(T,1)]\simeq0.7520248633528791908957935591.
$$

Both are enclosed by rational intervals whose combined width is less than
`10^-75`. Thus **C wins at fixed horizon one**, although **M wins with an
Exp(1) horizon**, even though the observation horizons have the same mean.
The law of the horizon, not just its mean, matters.

For completeness, the enclosure has the following proof. Uniformize the
population generator at `q >= max(outgoing rate)`, keeping extinction as
a surviving state and progression absorbing. Let `s_j` be the exact
probability of no progression after `j` uniformized jumps and put `x=q*t`.
Conditional on exactly `n` Poisson jumps before a fixed horizon `t`, the
`n+1` expected spacings are all `t/(n+1)`, independently of the embedded
chain. Hence

$$
\mathbb E[\min(T,t)]=e^{-x}\sum_{n=0}^{\infty}
\frac{x^n}{n!}\frac{t}{n+1}\sum_{j=0}^n s_j.
$$

Every conditional mean in this sum lies between zero and `t`. Let `A_K`
be the sum through `K` including those means and
`B_K=sum_{n=0}^K x^n/n!`. For `K+2>x`, the omitted exponential series is
bounded above by

$$
R_K=\frac{x^{K+1}}{(K+1)!}\frac{1}{1-x/(K+2)}.
$$

Positive terms and decreasing successive ratios prove this geometric bound.
It follows that

$$
\frac{A_K}{B_K+R_K}\leq\mathbb E[\min(T,t)]
\leq\frac{A_K+tR_K}{B_K}.
$$

All terms and bounds are computed with `Fraction`; this controls both
arithmetic and truncation error. The recorded check uses `q=103`, `t=1`,
and `K=350`. It establishes this particular endpoint counterexample, not a
fixed-horizon mutation-threshold theorem.

`simulation.py` implements the population events directly and imports none of
the verifier, its matrices, or its formulas. Each trial includes an independent
rate-one censor clock. Sampling the winning exponential clock is an exact
embedded-chain simulation for progression before censoring; no time-step or
tau-leaping approximation is used. To reduce variance, the two arms share
the minimum of each corresponding event rate and have separate residual
clocks. Each marginal retains exactly its stipulated rate. They share an
independent censor time. Trial pairs are independent; uncertainty in a
difference uses the **paired** sample variance, not a false independence
assumption between arms.

Recorded reproducible check: 300,000 pairs at each point, with seeds
`20261117`, `20261118`, `20261119`, respectively:

| `mu` | Exact `U_C-U_M` | Simulation | Estimated paired standard error |
| --- | --- | --- | --- |
| `0.25` | `0.0037388800` | `0.0039866667` | `0.0004943605` |
| `2.5` | `-0.0134237108` | `-0.0138666667` | `0.0007359988` |
| `25` | `0.0026800699` | `0.0025166667` | `0.0002868047` |

The three estimates are each within one estimated standard error of the
exact gap and have the correct sign by more than eight estimated standard
errors in this fixed-seed run. These are simulation diagnostics, not a
probabilistic proof or a confidence statement about patients. The test
requires six-standard-error agreement and a six-standard-error sign margin.
The mathematical proof and exact rational certificates establish the signs.
The 12-test file, including all 900,000 paired trials, runs locally in about
5 seconds, below the requested 60-second bound; hardware may affect timing.

## 7. Falsification attempts and hostile self-review

The following failures and restrictions were retained rather than hidden:

| Attempt or altered regime | Outcome and implication |
| --- | --- |
| Claim the mutation/competition question is new | Rejected by Hansen–Woods–Read 2017, Hansen–Read 2020, and Alvarez–Viossat 2024. The exact theorem above is the narrower candidate contribution. |
| Compare ordinary means and obtain a finite cutoff | Fails: both means are infinite. The objective change to A8 is explicit. Conditioning only one arm on progression would be an invalid comparison. |
| Refill S after S has reached zero | Rejected. An early exploratory model did this after a death from `(0,2)`. It introduced an undeclared sensitive reservoir and produced a different polynomial. It is not the model proved here. Dedicated tests prevent this error. |
| Keep the original exploratory rates after correcting that boundary | With `(a,b,d,k,z)=(11,1,0.1,0.1,0.1)`, the corrected comparison has one reverse-direction crossing, not the initial two. This parameter set is retained among the verifier's exact checks. The final two-crossing example uses different, fully stated rates. |
| Stop as soon as R reaches two | Changes the endpoint and removes the pure-R continuation. Our proof does not apply. The code explicitly checks that conversion into `(0,2)` is not progression. |
| Remove competition: change only `a` from 100 to 1 | Every coefficient of `Psi` becomes positive: M has larger `U` for all `mu>0`. The double reversal is not universal. |
| Use slower censoring: change only `z` to `1/10` | Every coefficient of `Psi` is positive: M wins for all `mu>0`. |
| Use faster censoring: change only `z` to 10 | Every coefficient of `Psi` is negative: C wins for all `mu>0`. |
| Increase sensitive clearance: change only `k` to 20 | Every coefficient of `Psi` is positive: M wins for all `mu>0`. |
| Substitute eventual progression probability for U | M wins throughout, by (9). A larger expected capped time does not imply a larger extinction probability or stochastic dominance. |
| Substitute a fixed deterministic observation time | At `mu=5/2`, fixed horizon one favors C, while Exp(1) favors M. The exact uniformization certificate above reproduces this failure of endpoint substitution. No identical fixed-horizon threshold is asserted. |
| Interpret high mu as a plausible per-division mutation probability | Rejected. It is an unrestricted conversion rate in a two-cell counterexample. No realistic genomic-mutation inference follows. |
| Interpret M as a maximum over the same bounded controls available to C | Not proved. C uses instantaneous paired interventions and regeneration from surviving S; it is a singular mathematical comparator. |

The four all-positive/all-negative coefficient regimes in the table are
checked exactly by the test, as well as by direct occupancy calculations.
They falsify an attempted parameter-independent two-reversal claim.

The main example has `k=2` greater than the mixed-state compensated R-birth
intensity `b=1`, but this **does not** establish pharmacological feasibility
or equivalence of admissible control classes. The assumption remains visible.

All denominators, inequality directions, initial conditions, zero-rate
exceptions, pure-R transitions, and root multiplicities were independently
rederived during adversarial AI review. This is separate computational
checking, **not external scientific peer review**. In particular:

- `h_C-h_M` and `U_C-U_M` have opposite signs; the test includes exact
  rational witnesses to catch a sign reversal.
- No negative birth rate is evaluated at progression: only the two positive
  numbers `a,b` are specified, and the process stops at total three.
- `d>0` and `z>0` are explicit theorem assumptions; taking either to zero
  can change the meaning of the comparison and is not silently substituted.
- The rational root intervals certify location; Descartes' rule and sign
  brackets certify completeness. A plot or sampled grid alone would not.
- Mutation-free containment remains all-sensitive forever. A positive but
  tiny conversion rate is not treated as zero in an infinite-horizon claim.

## 8. Limits and remaining open work

This result says that a particular simple stochastic theory does not support
a universal mutation-only rule. It does not decide whether containment or
MTD is preferable in a tumour. The exact candidate contribution is small:
a tractable model family, a cubic criterion, a sharp high-rate coefficient,
and a certified counterexample to a universal upper mutation cutoff for its
stated objective.

The following are **not proved** and must not be inferred:

- persistence of the two-reversal pattern in a large-population,
  division-linked, rare-mutation limit;
- persistence with both policies derived from one shared, bounded
  pharmacodynamic control system and finite sensitive regrowth;
- a branching-process establishment approximation with a controlled error
  bound for this comparison;
- the same ranking for fixed-horizon means, conditional progression times,
  medians, survival probabilities, or other observation distributions;
- robustness over clinically plausible parameter ranges, spatial structure,
  multiple resistance types, drug toxicity, patient heterogeneity, or
  measured tumour dynamics.

Those are possible subsequent research questions, not conclusions or
implicitly accepted conjectures. No numerical evidence here establishes them.
The result is scoped entirely to the new directory and its dedicated test;
it does not publish a biological finding into the repository's catalog or
results registry.
