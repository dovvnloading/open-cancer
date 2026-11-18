# A sharp short-horizon sequencing threshold with collateral sensitivity and cross-resistance

**Nonclinical computational mathematics. No treatment advice.** These are
statements about an explicitly specified stochastic process, with synthetic,
uncalibrated parameters. This AI-assisted research note has not received external
scientific peer review. A proof about this process, agreement of independent
calculations, and passing software checks do not establish biological validity.

## Result and scope of the contribution

For the four-type branching process below, starting with only sensitive cells,
an exposure-matched course consisting of any **fixed** number of A-then-B periods
has a larger probability of complete extinction **by the end of a sufficiently
short course** precisely when

\[
\boxed{c>\frac{1-w}{w},\qquad 0\le p<1.}
\]

Here `c` is one-way collateral sensitivity, `w` is the fraction of single-mutant
offspring resistant to A, and `p` is the conditional probability that a resistant
offspring acquires resistance to both drugs at once. At equality, combination
strictly wins at fourth order. At `p = 1`, the schedules tie exactly. Theorem 2
also proves a unique exact finite-course threshold near this boundary for short
courses on any fixed bounded collateral-sensitivity interval. Corollary 3 treats
both starting orders and characterizes existence of a locally improving periodic
schedule. Theorem 5 identifies the first nonzero effect of cross-resistance on
the nearby finite-course threshold: for symmetric single-mutant supply and one
A-then-B pair, its third-order coefficient has derivative `-mu*k^2/1440` with
respect to `p`. This is a comparison of thresholds, not a claim that
cross-resistance improves the absolute probability of eradication.

**This is a deliberate narrowing of the original global question.** Sequencing
advantages, cross-resistance, four-type extinction calculations, and
exposure-matched combination comparisons are already in the literature. We did
not establish a closed global dominance classification for arbitrary course
lengths. Indeed, the same parameters produce opposite rankings at two course
lengths below. For arbitrary specified finite courses, `verifier.py` instead
returns rigorous rational probability enclosures and, when separated, a strict
comparison certificate. It never labels overlapping intervals a tie.

The candidate model-specific refinement is the analytic local threshold, its
strict boundary coefficient and period dependence, the accompanying bounded
finite-course threshold theorem, and the leading nonzero cross-resistance
correction to that local curve, certified by an exact polynomial identity. We found no matching result in the primary
sources audited through **17 November 2026**. This is a qualified novelty
assessment, not proof of exhaustive priority.

## 1. Literature audit before model development

The audit used primary papers and author manuscripts, following citations from
the classical resistance literature and searches for combinations of
`branching`, `extinction`, `cross-resistance`, `collateral sensitivity`,
`sequential`, `cyclic`, `equal dose`, `short time`, and `Taylor`. Recent searches
included 2025–2026 papers. No patient records or third-party datasets are used.

| Primary work | What is already known; distinction from this note |
| --- | --- |
| [Komarova & Wodarz (2005), PNAS, doi:10.1073/pnas.0501870102](https://doi.org/10.1073/pnas.0501870102) | Stochastic emergence of resistance and prevention by combinations. Branching-process resistance modeling is foundational prior work. |
| [Iwasa, Nowak & Michor (2006), Genetics, doi:10.1534/genetics.105.049791](https://doi.org/10.1534/genetics.105.049791) | Resistance during clonal expansion and pre-existing resistant lineages. Our sensitive-only initial condition is a significant restriction relative to that motivation. |
| [Katouli & Komarova (2011), *The Worst Drug Rule Revisited*, doi:10.1007/s11538-010-9539-y](https://doi.org/10.1007/s11538-010-9539-y) | The closest mathematical precedent: four-type birth–death processes, direct cross-resistance, exact generating-function extinction equations, cyclic schedules, and finite-duration comparisons. Its average combination is the dose-matched rapid-switch limit. These ingredients and numerical ranking reversals are **not new** here. See Sections 2, 5–7 and [the primary PDF](https://link.springer.com/content/pdf/10.1007/s11538-010-9539-y.pdf). |
| [Bozic et al. (2013), eLife, doi:10.7554/eLife.00747](https://doi.org/10.7554/eLife.00747) | Explicitly includes mutations conferring resistance to multiple drugs at once. Simultaneous therapy outperforms the sequential scenarios studied, with and without cross-resistance. It would be incorrect to describe this paper as covering independent mutations only, or as a universal comparison under our exposure and sensitivity assumptions. |
| [Nichol et al. (2015), PLoS Computational Biology, doi:10.1371/journal.pcbi.1004493](https://doi.org/10.1371/journal.pcbi.1004493) | Sequence-dependent evolutionary steering on fitness landscapes. Noncommutativity of drug sequences is established; the endpoint and process differ from four-type branching extinction. |
| [Chakrabarti & Michor (2017), Cancer Research, doi:10.1158/0008-5472.CAN-16-2871](https://doi.org/10.1158/0008-5472.CAN-16-2871) | Time-dependent branching calculations with pharmacokinetics, interactions, and toxicity constraints. Schedule rankings depend on those assumptions. We omit those mechanisms. |
| [Yoon, Vander Velde, Marusyk & Scott (2018), doi:10.1007/s11538-018-0434-2](https://doi.org/10.1007/s11538-018-0434-2); [Yoon, Krishnan & Scott (2021), doi:10.1007/s00285-021-01671-6](https://doi.org/10.1007/s00285-021-01671-6) | Proved deterministic scheduling results for collaterally sensitive subpopulations, with a shaping stage and rapid cycling; the former also examines stochastic extinction times. These do not supply the joint direct-cross-resistance/CS threshold for our exact extinction endpoint. |
| [Maltas & Wood (2019), PLoS Biology, doi:10.1371/journal.pbio.3000515](https://doi.org/10.1371/journal.pbio.3000515) | Diverse experimental collateral-sensitivity profiles and control strategies. This cautions against treating our scalar, one-way `c` as a generally adequate biological description. |
| [Nyhoegen & Uecker (2023), doi:10.1098/rsif.2022.0793](https://doi.org/10.1098/rsif.2022.0793) | Sequential antibiotic models show the importance of dose response, overlap, and interactions. Their comparisons have different endpoints and pharmacology. |
| [Nyhoegen, Bonhoeffer & Uecker (2024), doi:10.1111/eva.13764](https://doi.org/10.1111/eva.13764) | Exact multitype branching extinction probabilities under a fixed total dose, principally comparing combinations with monotherapy. Exact extinction products and fair-dose constraints are prior work. |
| [Patil, Ahmed, Viossat & Noble (2026), Genetics, doi:10.1093/genetics/iyaf255](https://doi.org/10.1093/genetics/iyaf255) | Online 26 November 2025; February 2026 issue. An evolutionary-rescue approximation and simulations analyze two-strike timing. [Appendix F](https://pmc.ncbi.nlm.nih.gov/articles/PMC13270307/) explicitly gives a four-type, equal-total-toxicity argument for possible sequential superiority. It is not an exact extinction dominance theorem with our two resistance mechanisms. |
| [Molina-Hernández et al. (2026), Physical Review E, doi:10.1103/zty8-yj5p](https://doi.org/10.1103/zty8-yj5p); [author manuscript v4](https://arxiv.org/html/2510.01808v4) | Four-genotype stochastic extinction with collateral sensitivity, density dependence, and switching optimization. The explanatory approximation is explicitly heuristic. The presented model excludes direct sensitive-to-double-mutant transitions; the authors state that adding them leaves the qualitative results unchanged. Four-genotype CS extinction itself is not new. |
| [Dabi et al. (18 September 2026), PLoS Computational Biology, doi:10.1371/journal.pcbi.1014803](https://doi.org/10.1371/journal.pcbi.1014803) | Direct simulation comparisons of sequential extinction therapy and combination, including cross-resistance, already show either ranking. This resolves the broad question of whether sequencing can win; it does not provide the local exact CS criterion proved here. The final publication, not its older preprint, was considered. |

We claim neither the first sequencing advantage nor a new general PGF method.
The open/unproved part left by this note is a global exact phase classification
with arbitrary treatment lengths, initial mixtures, and more realistic
pharmacodynamics. No such conjecture is presented as a theorem.

## 2. Minimal process and equal exposure

Types are sensitive `S`, A-resistant `A`, B-resistant `B`, and doubly resistant
`D`. Cells act independently; the state is in the nonnegative integer lattice.
An offspring event retains the parent and adds one child. Rates are per cell per
unit of an arbitrary mathematical time scale:

| Parent | Child | Rate |
| --- | --- | --- |
| Any type | Same type | `b` |
| S | A | `a = mu (1-p) w` |
| S | B | `r = mu (1-p) (1-w)` |
| S | D | `z = mu p` |
| A or B | D | `nu` |

Thus `b` denotes **clonal** offspring production. Total offspring rates are
`b+mu`, `b+nu`, `b+nu`, and `b`. Equivalently, each sensitive reproductive event
is clonal with probability `b/(b+mu)` and produces a resistant child with
probability `mu/(b+mu)`; conditional on the latter, probability `p` produces
direct double resistance. The cross-resistance parameter redistributes a fixed
total resistant-offspring rate. It is not the rate of adding a new mutation
channel while keeping both single-mutant rates fixed.

Let `u(t)` be the A dose fraction, with B fraction `1-u(t)` and `0 <= u <= 1`.
Death rates are

\[
\delta_S=d+k,\quad \delta_A=d+ck(1-u),\quad
\delta_B=d+ku,\quad \delta_D=d.                         \tag{1}
\]

The factor `c >= 1` represents **one-way** collateral sensitivity: A resistance
raises sensitivity to B. The opposite direction has no extra sensitivity.
Resistance is complete to the named drug; `D` still has natural death. All
births and mutation rates are dose independent. The mutation graph is
irreversible and **decomposable**, not irreducible.

These are modeling assumptions, not measured pharmacology. Linear additive kill
and linear additive toxicity give the instantaneous budget `u+(1-u)=1`. The
comparisons additionally have identical exposure to **each** drug:

\[
\int_0^T u(t)\,dt=\int_0^T(1-u(t))\,dt=T/2.            \tag{2}
\]

`C` is simultaneous half-dose combination, `u = 1/2`. `AB_m` repeats A alone
for `T/(2m)` then B alone for `T/(2m)`, exactly `m` times. `AB_1` is the
two-block sequential course. `BA_m` reverses the order within each period.
The finite positive integer `m` is fixed when taking `T` to zero.

The endpoint is the probability that the population is **exactly zero by time
T**. Zero is absorbing. This is not a thresholded population count, expected
burden, or probability of ultimate extinction after unspecified future therapy.

## 3. Exact extinction probabilities for all finite courses

Let `beta_ij` be the offspring-rate matrix in Section 2 and
`B_i = sum_j beta_ij`. For a constant dose fraction, its backward PGF flow
`Phi_u(t, x)` solves

\[
\dot q_i=F_i^u(q)
=\delta_i(u)(1-q_i)+q_i\sum_j\beta_{ij}(q_j-1),\qquad q(0)=x. \tag{3}
\]

To derive (3), condition on the first event: death contributes `1`, a child of
type `j` contributes `q_i q_j` by independence of the surviving parent and
child, and no event contributes `q_i`. Subtract `q_i`, divide by elapsed time,
and take the limit. Bounded per-cell rates and domination of total births by a
Yule process of rate `max_i B_i` imply nonexplosion on every finite interval.
The PGF interpretation keeps the real solution in `[0,1]^4` and ensures its
uniqueness there. There is no population-size or mutation-order approximation.

For chronological phases `(t_1,u_1),...,(t_L,u_L)`,

\[
q=\Phi_{u_1}(t_1,\Phi_{u_2}(t_2,\ldots\Phi_{u_L}(t_L,0)\ldots)). \tag{4}
\]

The last phase is evaluated first. In particular,

\[
q^C=\Phi_{1/2}(T,0),\quad
q^{AB_1}=\Phi_1(T/2,\Phi_0(T/2,0)),\quad
q^{AB_m}=(\Phi_1(h,\cdot)\circ\Phi_0(h,\cdot))^m(0),\ h=T/(2m). \tag{5}
\]

For a deterministic initial population `n=(n_S,n_A,n_B,n_D)`, the exact
extinction probability is

\[
P(T;n)=\prod_i q_i^{n_i}.                              \tag{6}
\]

Equations (3)–(6) are standard branching machinery; their rigorous rational
evaluation is described in Theorem 4. They apply to arbitrary finite sequences,
including mixed initial populations. The ordering theorems below have a narrower
initial-condition assumption. An ultimate-extinction calculation after a finite
course would require replacing the terminal zero vector by the appropriate
post-course extinction vector; none of our ordering claims transfers to that
different endpoint automatically.

## 4. Theorem 1: sharp short-horizon A-first criterion

**Assumptions.**

1. The independent-cell process is exactly (1)–(3), with no competition,
   back mutation, pharmacokinetic carryover, or drug interaction term.
2. `b > 0`, `mu > 0`, `nu > 0`, `k > 0`, and `d >= 0` are fixed finite rates.
3. `1 <= c < infinity`, `0 < w <= 1/2`, and `0 <= p < 1` are fixed.
4. The initial state is `(N,0,0,0)` for a fixed finite integer `N >= 1`.
5. The endpoint and constraints are exactly (2) and finite-time extinction by
   `T`. Combination means simultaneous half doses.
6. The candidate is `AB_m` for a fixed finite integer `m >= 1`.

Write `s = d+k` (not the type `S`), `a = mu(1-p)w`, `r = mu(1-p)(1-w)`,
and `Delta_m(T) = q_S^{AB_m}(T)-q_S^C(T)`. Then

\[
\Delta_m(T)=C_3T^3+O(T^4),\qquad
C_3=\frac{s k(ac-r)}{16m}.                             \tag{7}
\]

At `ac = r`, its first nonzero coefficient is

\[
\Delta_m(T)=C_4T^4+O(T^5),\qquad
C_4=-\frac{a c k^2s}{384m^3}
\left[(c-1)(8m^2-1)+2m(c+1)\right]<0.                 \tag{8}
\]

Consequently, for every fixed parameter tuple satisfying these assumptions,
there is `T_0 > 0` such that **for every** `0 < T < T_0`,

\[
P_{AB_m}(T;(N,0,0,0))>P_C(T;(N,0,0,0))
\quad\Longleftrightarrow\quad c>(1-w)/w.               \tag{9}
\]

For `c <= (1-w)/w`, the inequality is strictly reversed on that sufficiently
short interval, including equality in the threshold. The threshold is
independent of `p` on `p < 1`; its leading effect is attenuated by `1-p`.

**Degenerate boundaries.** If `p=1`, `mu=0`, or `k=0`, all these schedules tie
exactly from sensitive-only initial states. If `N=0`, extinction probability is
one for every schedule. These are separate statements, not substitutions into a
strict inequality. Positive `b,nu` ensure that the principal theorem genuinely
allows continuing reproduction and second-step resistance; the same expansions
extend to their zero boundaries.

### Complete proof

**Step 1: expand backward in time.** Let `h_i(t;T)` be the extinction probability
by `T` of one type-i cell present at calendar time `t`, so `h(T;T)=0`. Put
`y(t)=h(T-t;T)`; then `y(0)=0` and `y'=F^{u(T-t)}(y)`. Denote its homogeneous terms in
course length by `y_{i,1},y_{i,2},...`. Equivalently rescale time to `[0,1]`
and expand the solution of `dy/dx = T F^{u(T(1-x))}(y)` in `T`.

For `AB_m` define

\[
v(t)=1-2u(T-t),\qquad V(t)=\int_0^t v(x)\,dx,\qquad
J(t)=\int_0^t V(x)\,dx.
\]

`V` is a nonnegative triangular wave: on each backward period it rises with
slope `+1` over its B half, then falls with slope `-1` over its A half. For
combination, `V=0`. Define

\[
\alpha=d+ck/2,\quad \beta=d+k/2,\quad
\lambda_S=s+b+\mu,\quad H=ac-r,\quad
L=ac\alpha-r\beta,\quad M=ac^2+r.                       \tag{10}
\]

The first terms follow directly by integrating the constant part of (3):

\[
y_{S,1}=st,\quad y_{D,1}=dt,\quad
y_{A,1}=\alpha t+(ck/2)V,\quad
y_{B,1}=\beta t-(k/2)V.                                \tag{11}
\]

The degree-two sensitive term is independent of the schedule:
`y_{S,2}=-s lambda_S t^2/2`. Every D term is independent of the schedule.
For the A and B terms, integration of the degree-one part of (3) gives the
following **differences** from combination:

\[
\begin{aligned}
\Delta y_{A,2}&=-\alpha(ck/2)tV-(c^2k^2/8)V^2-(b+\nu)(ck/2)J,\\
\Delta y_{B,2}&=+\beta(k/2)tV-(k^2/8)V^2+(b+\nu)(k/2)J.
\end{aligned}                                        \tag{12}
\]

For example, `y_{A,2}=-y_{A,1}^2/2-(b+nu) int_0^t y_{A,1}`; subtract its
`V=0` version to obtain (12). Quadratic births enter at the next degree, so no
term has been discarded at degree two.

**Step 2: sensitive-cell differences.** Write `W=a y_A+r y_B`. Equations
(11)–(12) imply

\[
\Delta W_1=(kH/2)V,\qquad
\Delta W_2=-(k/2)[L tV+(b+\nu)HJ]-(k^2M/8)V^2.         \tag{13}
\]

In (3) for S, the first difference is therefore

\[
\Delta y_{S,3}'=(skH/2)tV.
\]

At degree four the only differing terms are

\[
\Delta y_{S,4}'=-\lambda_S\Delta y_{S,3}
                 +y_{S,2}\Delta W_1+y_{S,1}\Delta W_2. \tag{14}
\]

The clonal quadratic term contributes no difference at this degree because the
first two S terms agree. The direct S-to-D term contributes none because both
those S terms and the D terms agree. This explains why second-step mutation and
direct cross-resistance do not create extra boundary terms in (8).

Set `I_j=int_0^T t^j V(t) dt` and `I_V=int_0^T t V(t)^2 dt`. Integrating
(13)–(14), including `int_0^T t J(t) dt=(T^2 I_0-I_2)/2`, yields

\[
\begin{aligned}
\Delta y_{S,3}(T)&=(skH/2)I_1,\\
\Delta y_{S,4}(T)&=-\frac{sk}{2}\left[
 \lambda_SH(TI_1-I_2/2)+LI_2+(b+\nu)H(T^2I_0-I_2)/2\right]
 -\frac{sk^2M}{8}I_V.                                 \tag{15}
\end{aligned}
\]

**Step 3: evaluate the schedule moments.** For `m` triangles, each of base
`T/m`, integrate the two linear pieces in each period and sum the period
indices and their squares. The resulting exact moments are

| Moment | Value |
| --- | --- |
| `I_0` | `T^2/(4m)` |
| `I_1` | `T^3/(8m)` |
| `I_2` | `T^4(8m^2-1)/(96m^3)` |
| `I_V` | `T^4/(24m^2)` |
| `T I_1-I_2/2 = (T^2 I_0-I_2)/2` | `T^4(16m^2+1)/(192m^3)` |

For transparency, on period `j=0,...,m-1`, put `ell=T/m` and
`t=j ell+x`. The pieces are `V=x` for `0<=x<=ell/2` and `V=ell-x` otherwise.
This gives `int V=ell^2/4`, `int x V=ell^3/8`,
`int x^2 V=7 ell^4/96`, `int V^2=ell^3/12`, and
`int x V^2=ell^4/24`; substituting the sums of `j` and `j^2` proves every row.

Substituting in (15) proves (7) and also the general fourth coefficient used in
the executable verifier:

\[
C_4=-sk\left[
\frac{H(s+2b+\mu+\nu)(16m^2+1)}{384m^3}
+\frac{L(8m^2-1)}{192m^3}+\frac{kM}{192m^2}\right].     \tag{16}
\]

At `H=0`, `r=ac`, `L=ack(c-1)/2`, and `M=ac(c+1)`, giving (8).
All its factors are positive except its displayed minus sign, because `c>=1`.

**Step 4: justify the remainder and strict signs.** Put

\[
R=\max\{d+ck,b+\mu,b+\nu\}>0.
\]

On normalized backward time the constant coefficient of `F` has magnitude at
most `R`, its original linear coefficient has magnitude at most `2R`, and the
sum of quadratic coefficients in each component is at most `R`. Let `a_n(x)`
be the vector coefficient of `T^n`, with `a_n(0)=0`. Induction in the exact
coefficient recursion gives `||a_n(x)||_infinity <= R^n x^n`: for `n>=2`, the
derivative is bounded by

\[
2R R^{n-1}x^{n-1}+(n-2)R R^{n-1}x^{n-1}
=nR^n x^{n-1}.
\]

There are `n-2` quadratic products because the constant solution term is zero.
For `n=1`, the bound follows by direct integration. The resulting series and
its differentiated series converge for `RT<1`; substitution proves that it is
the unique solution, including across the finitely many switches. Thus after
degree `j` the difference between any two schedules has remainder at most

\[
\frac{2(RT)^{j+1}}{1-RT}.                             \tag{17}
\]

If `C_3 != 0`, its sign is guaranteed whenever

\[
0<T<\min\{1/(2R),\ |C_3|/(4R^4)\}.                  \tag{18}
\]

At the balanced boundary, replace `C_3,R^4` in (18) by `C_4,R^5`. These are
conservative but explicit positive horizons. Equation (7), then (8) at equality,
now proves the one-founder ordering. Since `x -> x^N` is strictly increasing on
`[0,1]` for finite `N>=1`, (6) preserves every strict sign. Finally, if `p=1`
or `mu=0`, only accessible S/D dynamics matter and those are independent of
`u`; if `k=0`, every type's rates are independent of `u`. This proves the
degenerate statements and completes the proof.

## 5. Theorem 2: an exact nearby finite-course boundary

Keep Theorem 1's rates, `w`, and `m` fixed. Choose any finite `C > c_0`, where
`c_0=(1-w)/w >= 1`, and allow `c` to vary in `[1,C]` and `p` in `[0,1)`.
There exists `epsilon > 0`, independent of `p` in that interval, such that for
every `0<T<epsilon` there is a unique critical value
`c_m(T,p)` in `(c_0,C)`, and

\[
P_{AB_m}>P_C\quad\Longleftrightarrow\quad c>c_m(T,p)
\quad (1\le c\le C).                                \tag{19}
\]

Equality of probabilities holds exactly at this critical value. The curve is
real analytic locally in `(T,p)` and extends to `T=0,p=1` after cancellation of
the vanishing gap factor. It satisfies, uniformly in `p`,

\[
c_m(T,p)=c_0+
\frac{c_0 k[(c_0-1)(8m^2-1)+2m(c_0+1)]}{24m^2}\,T
+O(T^2).                                             \tag{20}
\]

In particular, for symmetric single-mutant supply, `w=1/2`,

\[
c_m(T,p)=1+\frac{k}{6m}T+O(T^2).                       \tag{21}
\]

An explicit, conservative choice of the uniform horizon is available. Define

\[
\begin{aligned}
R_*&=\max\{d+(C+1)k,\ b+4\mu,\ b+\nu\},\\
A_*&=\frac{sk\mu w}{16m},\\
B_*&=\frac{\mu w c_0 k^2s}{384m^3}
       [(c_0-1)(8m^2-1)+2m(c_0+1)].
\end{aligned}
\]

Then one may take

\[
\epsilon=\min\left\{\frac1{2R_*},\frac{A_*}{8R_*^4},
\frac{A_*(C-c_0)}{8R_*^4},\frac{B_*}{8R_*^5}\right\}>0. \tag{21a}
\]

This gives an exact local threshold structure, rather than merely the sign of
one expansion. It does **not** assert (19) at arbitrary finite `T`, for unbounded
`c`, or with `m` depending on `T`. Bound (21a) is deliberately very small;
finite-course examples outside it require their separate numerical certificates.

**Proof.** The convergent coefficient construction proves joint analyticity near
`T=0` in the finitely many rates. Coefficients of `T^0,T^1,T^2` of the gap vanish
identically. At `p=1` the entire gap is identically zero. Consequently

\[
G(T,c,p)=\frac{\Delta_m(T,c,p)}{(1-p)T^3}
\]

has an analytic extension through `T=0` and `p=1`. At `T=0`,

\[
G(0,c,p)=\frac{sk\mu[wc-(1-w)]}{16m},\qquad
\partial_cG(0,c,p)=\frac{sk\mu w}{16m}>0.
\]

Compactness of `[1,C] x [0,1]` and continuity give strictly positive
`partial_c G` for all sufficiently small `T`, uniformly on that compact set.
The negative coefficient (8) after dividing by `1-p` gives
`G(T,c_0,p)<0` uniformly for small positive `T`, while
`G(T,C,p)>0`. The intermediate value theorem, strict monotonicity, and the
analytic implicit function theorem prove existence, uniqueness, and (19).
Substitution of `c=c_0+beta T+O(T^2)` into (7)–(8) gives
`beta=-C_4/(partial_c C_3)`, exactly (20). The strictly positive factor
`(1-p)T^3` restores the original gap sign. At `p=1` the actual probabilities
still tie for every `c`; the extended curve does not override that fact.

To verify the explicit bound (21a), complexify `c` and `p` within closed disks
of radius one around their real compact intervals. Then absolute death rates
are at most `d+(C+1)k`, and the sum of absolute offspring coefficients in the S
row is at most `b+4mu`, since `|p|,|1-p|<=2`. The coefficient induction in Step 4
therefore bounds each degree-j gap coefficient `D_j(c,p)` by `2 R_*^j` on these
disks. These coefficients are polynomials. Cauchy's coefficient estimate on
the unit circles bounds both `partial_p D_j` and
`partial_c partial_p D_j` by the same `2 R_*^j`. Because `D_j(c,1)=0`, write

\[
\frac{D_j(c,p)}{1-p}
=-\frac1{1-p}\int_p^1\partial_pD_j(c,v)\,dv.
\]

The same bounds thus hold for this quotient and its c derivative, including
their limits at `p=1`. For `R_*T<=1/2`, the error in `partial_c G` beyond its
constant `A_*` is at most `4 R_*^4 T`. The error in `G(T,C,p)` beyond
`A_*(C-c_0)` has the same bound. At `c=c_0`, the error beyond `-B_*T` is at
most `4 R_*^5 T^2`. Bound (21a) makes these errors at most one half of the
respective positive margins. This proves the asserted derivative and endpoint
signs on the explicit interval, and finishes the proof of Theorem 2.

## 6. Corollary 3: when some periodic order can improve locally

Allow `0<w<1`, retaining the positive rates, `p<1`, `c>=1`, and the
sensitive-only initial state. There exists a **fixed finite** `m` and a choice
of `AB_m` or `BA_m` which strictly beats combination for all sufficiently small
positive `T` **if and only if**

\[
(c,w)\ne(1,1/2).                                     \tag{22}
\]

If `p=1`, no order improves. This corollary is about the specified periodic
family, with `m` chosen before taking the limit; it is not an optimal-control
theorem over all measurable schedules.

**Proof.** Reversing the starting order replaces `V` and `J` by their negatives
and leaves `V^2` unchanged. Thus the cubic changes sign, the first two summands
inside (16) change sign, and its last summand does not. If `H != 0`, choose the
order with a positive cubic. If `H=0` and `c>1`, the B-first quartic is

\[
C_4^{BA}=\frac{ack^2s}{384m^3}
\left[(c-1)(8m^2-1)-2m(c+1)\right].                   \tag{23}
\]

It is positive for some finite integer `m`; for instance, any integer
`m>(c+1)/(c-1)` suffices. The remainder argument from Theorem 1 then applies.
The only remaining case is `c=1,w=1/2`, where `H=0` and both quartics equal
`-a k^2 s/(96m^2)<0` for every fixed finite `m`. This proves both directions.
Mutation bias alone can therefore create a local order advantage even at
`c=1`; the interpretation is not that every such advantage is caused by CS.

## 7. Theorem 4: rational certificates at any finite horizon

For rational nonnegative finite rates and rational piecewise-constant schedules,
the function `extinction_bounds` encloses every component of (4) using only
exact rational arithmetic. It handles any finite course length. Increasing
`order` and `bits`, with the step partition fixed, makes interval widths tend
to zero. For two exposure-matched schedules and any fixed initial state, strict
superiority holds **if and only if** their intervals eventually separate in that
direction under such refinement.

This last statement is a computable strict-comparison criterion, not a closed
global parameter threshold or a new general method. Exact equality need not be
decidable by a finite number of refinements; an overlap is reported as
`unresolved`.

### Proof of the implemented enclosure

On one constant phase let `B=max_i B_i`, `R=max(B,max_i delta_i)`, and `K=2R`.
For an approximate start `z` in the probability cube, write `q=z+v`. The
constant, linear, and quadratic coefficients of `v'=F(z+v)` satisfy

\[
\|F(z)\|_\infty\le R,\qquad
\|DF(z)\|_\infty\le\max_i(\delta_i+2B_i)\le3R,
\qquad \|\hbox{quadratic coefficients}\|\le R.          \tag{24}
\]

For the Jacobian bound, put `U_i=sum_j beta_ij z_j <= B_i`. Its diagonal is
`-delta_i-B_i+U_i+z_i beta_ii`, and its off-diagonal entries are
`z_i beta_ij >= 0`. The triangle inequality bounds the absolute row sum by
`delta_i+B_i-U_i+z_i B_i <= delta_i+2B_i`.
Also `F_i(z)=delta_i(1-z_i)-z_i sum_j beta_ij(1-z_j)` is a difference of
nonnegative quantities bounded by `delta_i` and `B_i`, proving the first bound.

The scalar series solution of `v'=K(1+v)^2,v(0)=0` is
`sum_{n>=1} K^n t^n`. Its coefficient recurrence dominates the absolute vector
recurrence in (24), since `K>=R`, `2K>=3R`, and `K>=R`. Hence the degree-`J`
Taylor remainder for a step `h` is at most

\[
\tau=\frac{(Kh)^{J+1}}{1-Kh}.                         \tag{25}
\]

The program chooses `Kh<=1/4`, evaluates the polynomial exactly with
`fractions.Fraction`, rounds each component down to a multiple of `2^-bits`,
and clips it to `[0,1]`. Rounding changes it by less than `2^-bits`; clipping
cannot increase its distance from the true probability.

For two initial PGF arguments in the cube, the product inequality for a
descendant population gives

\[
\|\Phi_u(h,z)-\Phi_u(h,z')\|_\infty
\le E[\hbox{number of descendants}]\|z-z'\|_\infty
\le e^{Bh}\|z-z'\|_\infty
\le\frac{\|z-z'\|_\infty}{1-Bh}.                     \tag{26}
\]

The first-moment bound follows by ignoring deaths and using maximum offspring
rate `B`. Here `Bh<=1/8`, and `exp(x)<=1/(1-x)` for `0<=x<1` follows by
comparing their power series. If the incoming coordinate error bound is `E`,
the outgoing error is therefore bounded by the exact rational update

\[
E_{new}=E/(1-Bh)+\tau+2^{-bits}.                       \tag{27}
\]

Induction over all reversed phases proves the returned intervals contain the
true probabilities. Zero-duration or zero-rate phases are identities. All
error arithmetic is rational and is never rounded down. For a fixed finite
number of steps, (25) and the rounding error tend to zero as `J,bits` increase;
so do the propagated bounds. Products as in (6) preserve enclosures. A strictly
positive gap eventually exceeds the shrinking total width; separated intervals
conversely prove its sign. This proves Theorem 4.

## 8. Theorem 5: the leading model-specific cross-resistance correction to the threshold

Under Theorem 2's assumptions, specialize to symmetric single-mutant supply
`w=1/2` and one A-then-B pair `m=1`. For any fixed
`0<=p_1<p_2<1`, the **exact** critical curve of Theorem 2 satisfies

\[
\boxed{
c_1(T,p_2)-c_1(T,p_1)
=-\frac{\mu k^2(p_2-p_1)}{1440}T^3+O(T^4).
}                                                     \tag{28}
\]

In particular `c_1(T,p_2)<c_1(T,p_1)` for all sufficiently small positive `T`.
The first two coefficients are independent of `p`:

\[
c_1(T,p)=1+\frac{kT}{6}
+\frac{k(-115b+2d+55k-2\mu+61\nu)}{720}T^2
+\gamma(p)T^3+O(T^4),\qquad
\gamma'(p)=-\frac{\mu k^2}{1440}.                     \tag{29}
\]

This does not say the fixed-duration threshold is exactly a cubic polynomial,
nor that the sign persists for long courses. At `p=1` the actual schedules tie;
an analytic continuation of a normalized threshold is not a unique ranking
threshold at that endpoint. Increasing `p` also changes absolute probabilities,
which are not ordered by (28).

### Proof, with an exact finite algebra certificate

The analytic threshold exists by Theorem 2. Put `h=T/2`, `eta=c-1`,
`a=mu(1-p)/2`, `z=mu p`, and `s=d+k`. The normalized gap has expansion

\[
\frac{\Delta_1(2h,c,p)}{askh^3}
=\frac{\eta}{2}+h f_4(c)+h^2 f_5(c)+h^3 f_6(c)+O(h^4). \tag{30}
\]

Only the following coefficients are required:

\[
\begin{aligned}
L_*&=34b+31d+35k+17\mu+17\nu,\\
f_4(1+\eta)&=-k/6-L_*\eta/24-3k\eta^2/8,\\
f_5(1)&=k(95b+51d+40k+29\mu+8\nu)/120,\\
f_5'(1)&=P_5/120,\qquad f_6(1)=-kP_6/720.
\end{aligned}                                        \tag{31}
\]

Here `P_5,P_6` are the complete polynomials specified by the following table:

| Monomial | Coefficient in `P_5` | Coefficient in `P_6` |
| --- | ---: | ---: |
| `b^2` | 225 | 1291 |
| `b d` | 760 | 2323 |
| `b k` | 855 | 1739 |
| `b mu` | 235 | 779 |
| `b nu` | 215 | 737 |
| `d^2` | 205 | 398 |
| `d k` | 461 | 625 |
| `d mu` | 330 | 668 |
| `d nu` | 280 | 187 |
| `k^2` | 243 | 248 |
| `k mu` | `613/2` | 468 |
| `k mu p` | `-105/2` | -107 |
| `k nu` | 213 | 126 |
| `mu^2` | 65 | 134 |
| `mu nu` | 105 | 83 |
| `nu^2` | 55 | 8 |

For a directly reproducible derivation of every entry, define the polynomial
Lie operator `L_u f = sum_i F_i^u(x) partial_i f`. The exact coefficient of
`h^j` in the gap is

\[
D_j=\left[
\sum_{r=0}^j\frac{L_0^{j-r}L_1^r x_S}{(j-r)!r!}
-\frac{(L_0+L_1)^j x_S}{j!}\right]_{x=0}.              \tag{32}
\]

This follows by expanding the two pullbacks in (5); combination over `2h`
has operator `exp(h(L_0+L_1))` because the vector field is affine in dose.
The order in (32) is essential: chronological A then B means the B pullback
acts on the A-evolved coordinate. Equation (32), differentiation and polynomial
multiplication give (31) and the table. The included
`algebra_certificate.py` implements (32) in a sparse polynomial ring over exact
rational numbers and checks these **generic identities**, without substituting
any parameter values or fitting an interpolant. This is the computer-assisted
finite algebra part of the proof; the identities are fully displayed above.

Its coefficient pruning is exact: a further Lie operation lowers degree in the
four PGF arguments by at most one. With `r` operations left, a term of PGF
degree greater than `r` cannot contribute at `x=0`. No Lie operation
differentiates `eta`, so discarding its degrees above those requested cannot
alter a retained coefficient. Accordingly it retains all of `D_3,D_4`, the
constant and linear eta terms of `D_5`, and the constant eta term of `D_6`.
For `j=3,4` the eta-degree cap is `j`, including every possible power after `j`
Lie operations, each of which increases eta degree by at most one.
For the cross-resistance derivative it works in variables `a,z`, sets
`mu=2a+z`, and applies `mu(partial_z - partial_a/2)`, which differentiates with
respect to `p` at fixed total `mu`. This independently yields from (31)

\[
\partial_p f_4=\partial_p f_5(1)=0,\quad
\partial_p f_5'(1)=-7\mu k/16,\quad
\partial_p f_6(1)=107\mu k^2/720.                       \tag{33}
\]

To finish the proof, insert the exact analytic root
`eta=alpha h+beta h^2+chi h^3+O(h^4)` into (30). At successive orders,

\[
\alpha=k/3,\qquad
\beta=L_*\alpha/12-2f_5(1),\qquad
\chi=L_*\beta/12+3k\alpha^2/4-2\alpha f_5'(1)-2f_6(1).
\]

The first two expressions are independent of `p`. Differentiating the last,
using (33), gives

\[
\partial_p\chi
=-2\left[-\frac{7\mu k}{16}\frac{k}{3}
          +\frac{107\mu k^2}{720}\right]
=-\frac{\mu k^2}{180}.
\]

Since `h^3=T^3/8`, this is precisely `partial_p gamma=-mu k^2/1440`.
The expression for `beta/4` gives the second coefficient in (29). Analyticity
makes the remaining term `O(T^4)` for each fixed pair `p_1,p_2`. Subtraction
therefore proves (28), including the strict sign for sufficiently small `T`.
This proof does not rely on the stochastic experiment or sampled parameter
checks. The test file also checks the sixth-order implication by a separate
population-state generator calculation as an independent falsification check.

## 9. Reproduction and independently checked examples

From the repository root, with ordinary Python 3.10 or later and no extra
dependency for these new files:

```sh
python research/combination-cross-resistance/verifier.py
python research/combination-cross-resistance/algebra_certificate.py
python research/combination-cross-resistance/simulation.py
python -m unittest discover -s tests -p 'test_research_combination_cross_resistance.py' -v
```

For narrower rational enclosures, use `verifier.py --order 24 --bits 96`.
Decimal strings printed by the CLI are rounded displays only. Certificate
decisions use exact rational endpoints. The underlying exact values can be
inspected through `Interval.lower` and `.upper`.

The default rational example has

`b=2/5, d=1/10, k=1, c=2, mu=3/5, p=1/4, w=1/2, nu=1/5`.

| Course / change | Combination probability | Candidate probability | Certified sign and approximate gap |
| --- | ---: | ---: | --- |
| `T=1/2`, sequential `AB_1` | 0.356231204349 | 0.356653738028 | positive, +0.000422533679 |
| `T=1/2`, periodic `AB_2` | 0.356231204349 | 0.356469676695 | positive, +0.000238472346 |
| `T=2`, sequential `AB_1` | 0.661922633865 | 0.659789047833 | negative, -0.002133586031 |
| `T=2`, periodic `AB_2` | 0.661922633865 | 0.661717360758 | negative, -0.000205273106 |
| `c=1,T=1/2`, `AB_1` | 0.354228650934 | 0.354175715499 | negative, -0.000052935435 |

At the default precision, the first gap's exact enclosure is contained in
`(0.00042253367,0.00042253369)`, and the `T=2,AB_1` enclosure is contained in
`(-0.002133587,-0.002133585)`. These are certificates for the displayed finite
courses; they are not claims that those durations lie below Theorem 2's
conservative uniform `epsilon` in (21a).

There is also a certified **cross-resistance-induced ranking reversal at a
fixed finite course length**. Hold the default rates fixed except
`T=1/10`, `c=5084894333/5000000000`, and vary `p`. The exact rational
enclosures for `q_S^(AB_1)-q_S^C`, computed with `order=28,bits=112`, lie in:

| Cross-resistance `p` | Strict gap enclosure, rounded outward | Ranking |
| --- | --- | --- |
| `1/10` | `(-9.502844e-13,-9.502841e-13)` | combination wins |
| `2/5` | `(6.289353e-13,6.289357e-13)` | A-then-B wins |

This synthetic collateral factor was selected by deterministic sign bracketing
to lie near the boundary; it is not a fitted biological parameter. The tiny
effect demonstrates why ordinary floating-point sign inspection is inadequate
here. These signs have exact rational certificates and are explicitly rerun in
`test_certified_finite_course_cross_resistance_reversal`. They are not claimed
to be statistically resolvable by the short Monte Carlo run or to establish a
global monotonic threshold outside Theorem 2's horizon bound.

### Stochastic check

`simulation.py` independently reconstructs population-level birth, death, and
mutant-offspring events and uses the direct Gillespie method. It does not import
the PGF solver or the theorem coefficients. It resamples the exponential clock
at a deterministic phase boundary, using memorylessness, and checks exact zero
population. An event-limit exceedance raises an error and is never counted as
survival or extinction.

The synthetic stress example uses
`b=1/20,d=2,k=1,c=8,mu=6,p=1/10,w=1/2,nu=1/50,T=3/10`.
Its deliberately large mutant-offspring rate makes the small scheduling effect
observable in a short test; it is **not** a plausible mutation-rate estimate for
a tumour. The seed was specified before running, without searching seeds for
agreement. With 150,000 trials per schedule:

| Schedule | Rational-solver midpoint | Seed | Extinctions | Monte Carlo estimate |
| --- | ---: | ---: | ---: | ---: |
| `C` | 0.425327448598 | 20261117 | 63978 | 0.426520000000 |
| `AB_1` | 0.431985219937 | 20261118 | 65139 | 0.434260000000 |
| `AB_2` | 0.428889779322 | 20261119 | 64364 | 0.429093333333 |

Each observed error is below two binomial standard errors. Tests allow six
standard errors to avoid treating one seeded experiment as a proof, and check
the observed sequential-versus-combination direction. The deterministic sign
certificates and analytic proof establish the mathematical result; simulation
is a separate implementation/falsification check. Random-generator details can
change exact counts across Python implementations; the probabilities and
statistical tolerance, not those particular counts, are the test contract.

### Independent deterministic check

The test file separately constructs the **population-state Markov generator**
and multiplies its chronological exponential series through fourth order using
exact fractions. It reconstructs transitions from the event definitions; it
does not call the PGF recurrence to obtain coefficients. With at most four
generator applications, only finitely many population states are reachable,
so these coefficients have no truncation error. This independently checks
both signs, both starting orders, several period counts, nonzero cross mutation,
and the equality boundary against (7), (8), (16), and (23).

The suite also checks per-drug exposure, rational enclosure refinement,
degeneracies, independent-founder products, a single-type birth–death solution,
finite-course reversals, and the sensitive-only restriction. It uses stdlib
only and is intended to finish within 60 seconds on CI.

## 10. Falsification log and hostile referee review

| Attempted stronger claim / likely error | Outcome retained in this note |
| --- | --- |
| “Cross-resistance was absent from Bozic and earlier cycling models.” | False. Corrected before formulating the contribution; Katouli–Komarova and Bozic explicitly include it. |
| “The existence of sequential superiority is new.” | False. The literature audit includes the September 2026 Dabi paper and earlier cyclic calculations. Novelty is restricted to the stated coefficients and local criteria. |
| “`c>(1-w)/w` determines the ranking for any course length.” | Falsified by certified `c=2,w=1/2` examples: `AB_1` wins at `T=1/2` and loses at `T=2`. Global extension is unproved and not claimed. |
| “A zero cubic coefficient means a tie.” | False. The quartic is strictly negative for A-first schedules under Theorem 1's assumptions. At `p=1`, all coefficients vanish for a separate structural reason. |
| “A-first and B-first simply have opposite gaps.” | False beyond the cubic. The `V^2` quartic term has the same sign in both orders. At balance with `c>1`, B-first can win while A-first loses; (23) identifies when. |
| “Cross-resistance only changes the magnitude of a scheduling advantage.” | True only of the leading term here. The fixed-parameter, zero-horizon boundary is unchanged for `p<1`, but Theorem 5 proves a nonzero third-order effect on the nearby finite-course threshold. A global interior cutoff is unproved. |
| “The threshold applies to pre-existing resistant populations.” | False as a general extension. An initially A-resistant lineage responds to schedule timing already at degree two; the solver/test includes an A-founder case where A-first is worse. Theorem 1 requires sensitive-only founders. |
| “Expected population or a Floquet exponent can verify extinction ordering.” | Rejected as insufficient. We use exact PGFs, finite-state generator coefficients, and zero-population simulation. |
| “Rounding or population truncation might create the finite-course sign.” | Fractions arithmetic plus (25)–(27) bounds numerical error without a population cutoff; comparisons require strictly separated intervals. The Monte Carlo event cap raises rather than censors. |
| “The local periodic result is uniform as switching frequency diverges.” | Not proved. Fix `m` before taking `T` to zero. An infinite-switch or jointly scaled limit requires another argument. |
| “A biological mechanism has been validated because agents agree.” | Rejected. Independent derivations and algorithms assess this mathematical model only. No biological or clinical assessment occurred. |

Two separate analytic derivations checked the cubic, boundary quartic, and
period moments. A further review checked the shifted-series majorant,
Lipschitz direction, rounding, and reverse phase composition. The independent
Markov-generator tests protect against shared ODE-expansion mistakes. These are
AI/computational self-reviews, not independent human scientific endorsement.

## 11. Limitations and unresolved questions

The modeled extinction probability can be tiny for a large initial population,
even where its ordering is strict. All rates here are uncalibrated, and the
examples are mathematical witnesses rather than predictions. The process omits
spatial structure, competition, quiescence, changing fitness costs, immune
effects, dose-response saturation, drug interactions, and pharmacokinetics.
One-way scalar collateral sensitivity and fixed irreversible mutation rates are
simplifications. Fully resistant offspring can still die naturally. Equal
linear exposure is only the chosen mathematical fairness constraint.

The local theorem does not solve arbitrary-horizon scheduling, optimizing switch
times, unequal duty fractions, reciprocal CS, pre-treatment clonal expansion,
or ultimate extinction under an infinite schedule. Theorem 4 allows certified
evaluation of specified finite courses but does not convert numerical witnesses
into a universal threshold. Near equality a comparison may remain unresolved
at the requested precision, and exact equality generally needs structural
reasoning. The existence of a globally monotone threshold in `c` or `p` is
**unproved**; no conjectured global extension is used by the tests or proof.

Only original mathematical text, original stdlib code, synthetic parameters,
and source links are included. No patient data or copied source datasets enter
this contribution.


## 12. Independent hostile referee audit (17 November 2026)

**Disposition: Theorem 1's fixed-period local equivalence is confirmed; the
finite-course extension (Theorem 2) and the computer-assisted sixth-order
identity (Theorem 5) retain their stated scopes.** None is an arbitrary-course
or biologically calibrated treatment theorem.

**Independent forward-process route.** Instead of differentiating the backward
PGF, form the continuous-time generator \(Q_u\) on *population-count*
states. One death removes one cell, and an offspring event retains the
parent and adds the child. For one sensitive founder the probability of
extinction after an A-first/B-second course is the \((S,0)\) entry of

\[
e^{(T/(2m))Q_A}e^{(T/(2m))Q_B}\cdots
e^{(T/(2m))Q_A}e^{(T/(2m))Q_B};
\]

for combination it is the corresponding entry of \(e^{TQ_{1/2}}\).
Because the generator is affine in dose, \(Q_{1/2}=(Q_A+Q_B)/2\).
Expand both matrix products as formal Taylor series. For coefficient \(T^j\),
at most \(j\) cell events can contribute; no population truncation is needed
to calculate that coefficient from one initial cell. Enumerating death,
clonal offspring, and all five distinct mutant-offspring channels gives
identical coefficients at orders one and two and the difference

\[
[T^3](P_{AB_m}-P_C)=\frac{(d+k)k(ac-r)}{16m}.
\]

At \(ac=r>0\), the next coefficient reduces to the negative rational
quantity in (8): its bracket consists of a nonnegative
\((c-1)(8m^2-1)\) plus the strictly positive \(2m(c+1)\).
The sign of the first *nonzero* coefficient therefore gives precisely the
strict, fixed-\(m\), sufficiently-short-\(T\) ordering (9). The original
PGF proof supplies the crucial **convergent remainder**; a formal
generator coefficient without that bound would not prove the strict
all-small-\(T\) claim. The independent forward-CTMC test implements this
calculation for multiple rational rate regimes and both phase orders.
The analytic bound, rather than sampled intervals, covers all admissible
parameters.

**Remaining dependency chain.** Theorem 2 divides the gap by
\((1-p)T^3\) before applying compactness and the implicit function theorem;
without that normalization a claim of uniformity near \(p=1\) would be
unjustified. Its explicit bound relies on a coefficient majorant in complex
\(p,c\) disks and requires a fixed finite \(m\), finite upper \(C>c_0\),
and \(p<1\). Theorem 5 uses a **finite symbolic Lie-polynomial
certificate**, not a closed-form elementary hand reduction of all sixth-order
terms. The separate forward-generator test checks rational examples and
the predicted \(p\)-slope, but only the generic polynomial identities
in `algebra_certificate.py` address all algebraic rate values. This is
an explicitly computer-assisted proof component, not a peer-reviewed
symbolic-algebra theorem.

**Quantifier and boundary audit.** At \(ac=r\), A-first loses at order four,
so the condition in (9) must remain strict. At \(p=1\) or \(\mu=0\), the stated identical-schedule conclusion
is guaranteed **from sensitive-only initial populations**; an initially
A-resistant cell can break either equality. If \(k=0\), however, drug
death rates are dose independent, so **all** initial populations tie.
The added regression distinguishes the p=1 initial-state boundary. A dedicated regression now distinguishes
these cases. The parameter \(p\) redistributes mutant offspring toward
the directly doubly resistant type; it is not proof that measured
cross-resistance is beneficial. Equation (28) compares two analytic
thresholds in a specified *short-course* limit; it does not order their
underlying eradication probabilities or apply at \(p=1\).

**Prior art and novelty downgrade.** Four-type cross-resistance, PGF
extinction, cycle-order effects and toxicity-budget comparisons appear
already in Katouli–Komarova (2011); collateral-sensitivity scheduling is
studied by Yoon et al. (2018, 2021); and Dabi et al. (September 2026)
already report simulated extinction-therapy advantages over combination
in the presence of cross-resistance. Molina-Hernández et al. (2026)
also study four-genotype extinction with collateral sensitivity.
Accordingly the broad advantages and machinery are **known prior results**.
The narrow coefficients (7)–(8), finite local threshold, and (28) are
**unverified-priority mathematical refinements**. "Sharp" is justified
only as the exact local *if and only if* with its equality case; "leading"
means the first nonzero term in this model's expansion, not first
discovery in the literature. The numerical finite-horizon enclosures
are witnesses, never a universal global classification. No contradiction
to the narrow theorems was identified.
