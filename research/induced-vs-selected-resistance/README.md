# When endpoint fluctuation counts identify resistance origins

**Nonclinical computational mathematics. Synthetic calculations only. No
treatment advice. Unreviewed research contribution, with complete conditional
proofs and reproducible arithmetic checks; no experimental or clinical
validation is claimed.** Literature search/access date: November 17, 2026.

## Result and scope

The broad premise that resistance-origin inference lacks identifiability
theorems is incorrect. Greene, Gevertz and Sontag (2019), Gunnarsson, Foo and
Leder (2023), and Gevertz et al. (2025) already prove identifiability results.
We explicitly pivot to **minimal finite designs using endpoint resistant-colony
count laws, with an unknown shared clone-size distribution**.

For the model specified below, one exposure level and **three independent
experimental arms**, measured immediately and after holidays of lengths
`H` and `2H`, identify stable pre-exposure clones, stable post-exposure singleton
acquisitions, reversible clones, both directions of switching, and the entire
unknown clone-size law. Means and second factorial cumulants in the three arms,
together with the immediate zero-count probability, identify the finite
parameters. The immediate full count law identifies the clone distribution.
Theorem 3 supplies an explicit global inverse. Three arms are necessary for a
design guaranteed to identify the stated model class; this is a distributional
statement, not a claim that three cultures suffice.

If switching back into tolerance is independently ruled out, **two arms**
suffice when the clone law has a finite third moment and positive probability
of a clone of size at least three (Theorem 2). With an unknown law supported on
sizes one and two, three arms are necessary. With singleton clones, stable
pre-existing and stable post-exposure origins remain confounded at any number
of holidays at one exposure level (Theorem 4).

These are conditional statements. In particular, unequal clone laws, coordinated
switching of a whole clone, uncalibrated detection, or clonal expansion after
induction can destroy identification even with infinitely many holiday arms.
Theorem 5 gives exact observational equivalences for these failures. The result
does not establish genetic versus epigenetic causation from colony counts, or
decide whether any dosing schedule helps a patient.

| Model and information | Sufficient design | Sharp qualification |
|---|---|---|
| Unknown shared clone law, one-way loss of tolerance, finite third moment, `P(J>=3)>0` | `0,H` at one exposure | One arm fails; Theorem 2 |
| Unknown shared law on `{1,2}`, both sizes present, one-way loss | `0,H,2H` | Any two holiday times fail on an open set; Theorem 4 |
| Unknown shared law, finite second moment, `P(J>=2)>0`, two-way switching | `0,H,2H` | No two-arm design identifies the whole class; Theorem 3 |
| Known clone law with positive second factorial moment, one-way loss | `0,H` | Uses calibrated clone moments; Theorem 4 |
| Singleton clones, one-way loss | Three holidays at one exposure plus an immediate arm at a second calibrated exposure | Four scalar count laws are necessary for four interior parameters; Theorem 4 |
| Singleton clones, two-way switching | No finite or infinite dose/holiday design in this family suffices | Stable plateau and tolerant equilibrium confound; Theorem 4 |
| Different unknown stable and reversible clone laws | No fixed-exposure holiday design suffices in general | Stable singleton mass can be reassigned to induction; Theorem 5 |

The positive unknown-kernel reconstruction rows require `c>0` and actual loss
of tolerance; absent-component and permanent-memory boundaries are in 5.5.

“Arm” means a population of independent replicate cultures under one condition.
Structural identifiability assumes its complete probability law is known.
Finite-sample error and power are separate questions (Theorems 6–7). We do not
claim a minimum number of observations for exact recovery of real parameters
or an unrestricted infinite-dimensional clone law.

## 1. Literature survey and novelty assessment

This is a targeted primary-literature survey, not an exhaustive systematic
review. Searches combined fluctuation assay, induced resistance, persister,
identifiability, colony counts, drug holiday, thinning, common clone distribution,
factorial cumulants and minimal experimental design; references and subsequent
papers were followed. The following distinctions determined the pivot.

| Prior work | What it establishes and the boundary of our claim |
|---|---|
| Luria–Delbrück [1], Lea–Coulson [2], Ma–Sandri–Sarkar [3] | Classical jackpot discrimination and mutant-count distributions/recurrences. Neither compound Poisson distributions nor fluctuation-based discrimination is new here. Our finite-growth clone law need not be the ideal infinite-growth LD law. |
| Lang–Murray [4] | Fits an LD distribution convolved with a Poisson postplating component. Adding a Poisson induction term is established practice, not the contribution. |
| Pisco et al. [5] | Quantitative work supports drug-associated non-genetic state transitions in a studied cell system. Selection can act on non-genetic variants, and induction need not be genetic. Thus the three labels are not disjoint molecular categories. |
| Shaffer et al. [6] | Fluctuation experiments, expression measurements and sorting support transient pre-resistance followed by drug-induced reprogramming in studied melanoma populations. This is substantially more biological evidence than counts alone provide. |
| Greene–Gevertz–Sontag [7] | Appendix A13–A17 proves identifiability for a specified two-state ODE using total population derivatives and four controls, including a ramp, from a sensitive initial population under stated positivity assumptions. It is already an identifiability theorem for spontaneous versus induced resistance. |
| Bokes–Singh [8] | Derives an exact distribution and moments for persisters adapting to proliferating resistance; fluctuations can separate adaptation rate and initial abundance. We do not claim first moment-based inference of persister kinetics. |
| Kuosmanen et al. [9] | Models dose-dependent mutation supply and evolutionary rescue. Its control conclusions are model dependent; they do not supply the unknown-clone endpoint design theorem below or justify clinical advice. |
| Russo et al. [10] | A modified two-step fluctuation framework combines early/late colony emergence with population measurements to compare spontaneous and persister-associated mutation rates. This is close experimental antecedence, but a different observation family and model. |
| Saint-Antoine–Grima–Singh [11] | Gives fluctuation-based inference of switching kinetics from state fractions. Closed-form switching inference is established. |
| Gunnarsson–Foo–Leder [12] | Proposition 1 identifies multitype branching division, death and switching parameters from initial slopes of phenotype-resolved means and covariances with isolated initial subpopulations. Our aggregate endpoint colony observations do not supply those measurements. |
| Bokes–Hlubinová–Singh [13] | Reversible transitions change the tail of a fluctuation distribution; their matched asymptotics are not a finite-arm injectivity theorem. We cite the substantive 2023 paper, not only the related 2025 conference abstract. |
| Gevertz et al. [14] | Proves local generic identifiability of a delayed two-state ODE from total count derivatives with a known common delay rate and reports software confirmation of generic global identifiability. Count-only identifiability is therefore not new in general. |
| Browning et al. [15] | Studies practical identifiability of phenotype adaptation from low-count data and model misspecification. Its focus differs from exact equality of endpoint distributions. |
| Wu–Gunnarsson–Foo–Leder [16] | Fits branching models to bulk count dose/time grids and detects induced plasticity. Detailed concentration/timepoint design effects remain a stated future direction; our theorem addresses a narrower assay family. |
| Puig–Valero [17], Gugushvili et al. [18], Potts–Tasche [19] | Factorial-cumulant thinning, compound-Poisson decompounding, and exponential-sequence inversion are established mathematical ingredients. We use and credit them. |
| Whiting et al. [20] | Infers resistance dynamics using lineage tracing and population counts; warns that holiday resensitization can also reflect phenotype-specific growth. Growth control in A4 is therefore consequential. |
| Iyer et al. [21] | Uses time-lapse data and modeling to study pre-/post-drug fate decisions and inherited lineage correlations. This motivates testing independent cellular switching rather than assuming it from a population-level response. |

**Narrow proposed contribution.** We did not locate the exact unknown-common-
kernel reconstruction, the two-versus-three-arm classification, and its explicit
failure strata in these sources. Theorems 2–4 combine established ingredients
into a specific global endpoint-count design theorem. Theorem 5 supplies
general-parameter falsifications of extensions that would be experimentally
tempting. The distribution recurrence, Poisson limit, Neyman–Pearson test and
concentration bound are supporting machinery, not claimed as new mathematics.
The exact priority of the synthesis remains subject to broader specialist
review. Passing the tests is not evidence of biological truth or priority.

## 2. Assay family, operational mechanisms, and assumptions

Use the following assumptions **together**. Each theorem states its additional
restrictions. None is an empirical conclusion of this contribution.

**A1 — Independent cultures and a common preparation.** Each arm uses fresh,
independently grown cultures from the same preparation and the same exposure
pulse unless an exposure contrast is explicitly specified. Cultures are not
technical aliquots counted as independent replicates. There is no shared random
batch multiplier in the mathematical model.

**A2 — Three independent event pools.** At the start of the holiday, the numbers
of stable pre-exposure clone events and reversible clone events are independent
`Poisson(a)` and `Poisson(c)`. Stable post-exposure acquisitions number
`Poisson(b)`, independently. Here `a,b,c>=0` are event intensities per culture,
not per-cell mutation rates. The mechanisms may coexist in the same culture.

**A3 — An unknown shared positive clone law.** Every pre-exposure stable or
reversible event independently leaves `J>=1` countable founder cells with a
common probability mass function `g_k` and probability generating function

\[
G(z)=\sum_{k\ge1}g_k z^k,\qquad G(0)=0,\quad G(1)=1.
\]

The law is otherwise unrestricted. It may have unbounded support. Each stable
post-exposure event leaves exactly one founder; there is no expansion of that
founder before the modeled readout. Independence of marks and event counts is
assumed. Shared `G` is a substantive calibration assumption: stable and
reversible events must have comparable timing/growth/detection kernels. The
theorems identify `G`; they do not establish that the two biological pools have
the same law.

**A4 — Growth-controlled holiday and independent cellular memory.** Initially
all cells in reversible clones are in the tolerant state `T`. During a holiday
each independently follows a two-state continuous-time chain

\[
S\mathrel{\mathop{\rightleftarrows}^{u}_{v}}T,
\qquad u\ge0\text{ for }S\to T,\quad v>0\text{ for }T\to S.
\]

To avoid arrow-label ambiguity, the verbal definitions of `u` and `v` govern.
Both are finite. There are no births, deaths, new stable acquisitions or
additional external sensitive cells entering this modeled pool during the
holiday. Stable pools retain competence. This is an idealized quiescent or
otherwise growth-controlled assay, not arbitrary drug-free regrowth. Founder
normalization alone does not guarantee A4 if descendants have correlated fates.

**A5 — A calibrated readout.** Each stable founder and each reversible founder
in state `T` at readout yields one observed colony; a founder in `S` yields
none. Counting is uncensored and error free in the mathematical observation
model. Subsequent colony initiation neither creates nor erases competence.
An actual colony assay must validate this mapping: a quiescent persister may
not initiate a colony, or a sensitive cell may become tolerant during readout.
Such assays do not automatically satisfy A5. Common detection thinning is
handled separately in Theorem 5; arm-specific detection is not harmless.

**A6 — The design and observations.** A baseline is the immediate arm `h=0`
**after the same pulse**, not an untreated culture. Standard designs use
`h=jH`, with known `H>0`. Only the total resistant-colony count `Y_j` per culture
is used in the theorem; no genotype or phenotype counts are observed. A full
law means all probabilities `P(Y_j=n)`, not a finite observed histogram.

**A7 — Molecular and causal interpretation.** The names `D`, `I`, `P` below
mean stable clonally amplified pre-exposure events, stable unamplified
post-exposure events, and independently reversible clonal events. They are
operational statistical components. A reversible phenotype can be selected;
drug induction can be non-genetic. The timing/mark assumptions do not prove a
mutation mechanism. To interpret `b` as a *drug-caused excess*, background
post-exposure acquisition must be absent, known, or separated with a calibrated
exposure control (Theorem 5.6).

Write

\[
\lambda=u+v,\quad \pi=\frac{u}{u+v},\quad
r=e^{-\lambda H},\quad p_j=\pi+(1-\pi)r^j.
\tag{1}
\]

Thus `0<=pi<1`, `0<r<1`, and `p_0=1`. The one-way submodel has `u=0`, hence
`p_j=q^j`, `q=exp(-vH)`. Complete loss `q=0` is an optional mathematical
endpoint, not a finite transition rate. Define `m_k=E[(J)_k]`, with falling
factorial `(J)_k=J(J-1)...(J-k+1)`.

## 3. Exact laws and the relation to large cultures

### Theorem 1 — Count law, recurrence, and a controlled rare-event limit

**Assumptions:** A1–A6. Moment identities through order `k` additionally require
`E[J^k]<infinity`. The finite-opportunity assertion assumes its explicitly
specified independent-opportunity construction, not an unspecified branching
process.

The exact probability generating function is

\[
F_j(z)=\exp\{a[G(z)-1]+b(z-1)
                 +c[G(1-p_j+p_jz)-1]\},\quad 0\le z\le1.
\tag{2}
\]

Let `L_j=-log F_j(0)` and let the positive log-PGF coefficients be

\[
\ell_{j,k}=a g_k+b\mathbf1_{k=1}
 +c\sum_{l\ge k}g_l {l\choose k}p_j^k(1-p_j)^{l-k},\quad k\ge1.
\tag{3}
\]

Then `L_j=sum_{k>=1} ell_{j,k}` and

\[
P(Y_j=n)=e^{-L_j}w_{j,n},\quad w_{j,0}=1,\quad
w_{j,n}=\frac1n\sum_{k=1}^{n}k\ell_{j,k}w_{j,n-k}.
\tag{4}
\]

The factorial cumulants `kappa_{k,j}=(d^k/dz^k)log F_j(z)|_{z=1}` obey

\[
\mu_j:=\kappa_{1,j}=b+m_1(a+cp_j),\qquad
\kappa_{k,j}=m_k(a+cp_j^k),\quad k\ge2.
\tag{5}
\]

In particular, `kappa_2=Var(Y)-E[Y]`, and

\[
-\log P(Y_0=0)=a+b+c. \tag{6}
\]

For a population-size scaling, set `M=a+b+c>0`. Let there be `N>=M`
independent opportunities, each active with probability `M/N`, and on activity
choose type `D,I,P` with probabilities `a/M,b/M,c/M`; attach the same marks and
holiday dynamics. If its count law is `Q_{N,j}`, then

\[
\|Q_{N,j}-\mathcal L(Y_j)\|_{\rm TV}\le M^2/N.
\tag{7}
\]

For finitely many independent cultures and arms, the bound adds over cultures.
There is no upper bound of a handful of cells in (2), and (7) is a proven
large-opportunity limit. The case `M=0` is identically zero.

**Proof.** Solve `p'(h)=u(1-p(h))-vp(h)`, `p(0)=1`, to obtain (1). Conditional
on a reversible clone size `l`, its competent founder count is
`Binomial(l,p_j)` by independent cellular switching, so its PGF is
`G(1-p_j+p_j z)`. For `N~Poisson(t)` iid marks with PGF `B`, conditioning on
`N` gives `E[B(z)^N]=exp(t(B(z)-1))`. Multiplication of the three independent
PGFs proves (2). Expanding the nonnegative coefficients gives (3); any
zero-size thinned clone contributes zero to the exponent. Since `G(1)=1`,
the constant log coefficient is minus the sum of positive coefficients.
Differentiating `F'_j=F_j sum_{k>=1}k ell_{j,k}z^{k-1}` and equating
coefficients proves (4). Differentiation at one is justified by the stated
finite moments (or monotone convergence of the nonnegative derivative series)
and gives (5). At `j=0`, `G(0)=0`, proving (6).

For (7), `TV(Bernoulli(t),Poisson(t))=t(1-exp(-t))<=t^2` for `0<=t<=1`:
the Bernoulli law exceeds the Poisson law only at one, and by that amount.
Couple the `N` opportunity counts independently and attach identical iid type
and clone marks whenever counts agree. The chance any count differs is at
most `N(M/N)^2`. The sum of the Poisson opportunity counts is `Poisson(M)`;
independent type marking produces the three independent Poisson pools.
Counting and thinning cannot increase total variation. The same coupling and
union bound apply to a product of culture laws. This proves (7). **End proof.**

This construction does not prove a limit for every finite branching model:
overlapping mutant lineages, population feedback, or dependence between
opportunities need separate analysis. The classical infinite-growth LD limit
can lack the moments required below. Finite-growth clone laws or unbounded
laws with the stated moments qualify; the moment-free law (2) still exists
for any positive finite-valued integer mark law.

## 4. Identification theorems

### Theorem 2 — Two arms identify a common unknown clone law under one-way loss

**Assumptions:** A1–A7, `u=0`, `c>0`, `a,b>=0`, `H>0`, `0<q<1`,
`E[J^3]<infinity`, and `P(J>=3)>0`. The optional endpoint `q=0` is also
identified as a retention probability. No clone-size distribution is supplied
to the inverse.

The two laws `F_0,F_1` uniquely identify `(a,b,c,q,G)`. To see precisely which
statistics are sufficient, put

\[
R_k=\frac{\kappa_{k,1}}{\kappa_{k,0}},\ k=2,3,\qquad
d=\frac{1-R_3}{1-R_2}-1.
\]

The inverse is

\[
q=\frac{d+\sqrt{d^2+4d}}2,\qquad
x=\frac{R_2-q^2}{1-q^2},\qquad
V=\frac{\mu_0-\mu_1}{(1-x)(1-q)},
\tag{8}
\]

\[
b=\mu_0-V,\quad A=-\log F_0(0)-b,\quad
a=xA,\quad c=(1-x)A,
\tag{9}
\]

\[
G(z)=1+\frac{\log F_0(z)-b(z-1)}{A}.
\tag{10}
\]

Thus two means, two second factorial cumulants, two third factorial cumulants
and the immediate zero probability identify the finite parameters. The full
immediate law is needed for the unrestricted entire `G`; those seven scalars
alone do not identify an arbitrary infinite-dimensional distribution.

One arm cannot identify this model class. Consequently two is the minimum
number of arm laws for a design guaranteed to identify it.

**Proof.** Set `A=a+c>0`, `x=a/A in [0,1)`. Equation (5) gives
`R_k=x+(1-x)q^k`. Both denominator cumulants are positive because `m_3>0`
implies `m_2>0`. Also `1-R_2=(1-x)(1-q^2)>0`. Hence

\[
\frac{1-R_3}{1-R_2}
=\frac{1-q^3}{1-q^2}=1+\frac{q^2}{1+q}.
\]

The function `q^2/(1+q)` is strictly increasing on `[0,1)`; its derivative
is `q(2+q)/(1+q)^2`, positive for `q>0`, and the endpoint zero is unique.
Its nonnegative inverse solves `q^2-dq-d=0`, yielding (8). The equation for
`R_2` then gives `x`. Moreover `mu_0-mu_1=c m_1(1-q)>0`, so
`V=A m_1`; (5) and (6) prove (9). All divisions used are by strictly
positive quantities. With `A>0`, (2) at baseline gives (10). This recovers
every parameter and every coefficient of `G`, proving injectivity. If a
candidate statistic vector does not yield nonnegative intensities and a valid
PGF in (10), it is outside the model; the algebra is not a projection or a
finite-sample estimator guaranteed to produce admissible fits.

At baseline alone, only `a+c` appears, and `q` is absent. For any one positive
holiday, an alias also exists with arbitrarily large clone support: choose
`g_k>0` for `k=1,...,M`, any `M>=3`, `x in (0,1)` and `b>0`. Put
`nu_k=(a+c)g_k`. Its positive log coefficients are

\[
\ell_j-b\mathbf1_{j=1}=x\nu_j+
 (1-x)\sum_{k=j}^{M}\nu_k{k\choose j}p^j(1-p)^{k-j}.
\tag{11}
\]

For fixed `x,p,b`, this linear map in `nu` is upper triangular with positive
diagonal `x+(1-x)p^j`. Perturb `p` slightly within `(0,1)` and solve the same
system backwards. Its inverse varies continuously, so every `nu_k` remains
positive for a sufficiently small perturbation. Normalizing the new `nu`
gives a different valid kernel and intensities with the *same complete count
law*. A positive holiday length converts each such `p` to a distinct finite
one-way rate. Thus one positive arm also fails. The construction holds for
every `M`, not only a small-cell example. **End proof.**

For a genuine two-way chain, the same calculation identifies the single
occupancy `p(H)` but cannot separate `u` and `v`: infinitely many chains have
the same one-time occupancy. The next theorem removes this ambiguity and
needs only a second moment.

### Theorem 3 — Three arms identify two-way switching using second moments

**Assumptions:** A1–A7, `c>0`, `a,b>=0`, `u>=0`, `v>0`, finite rates,
known `H>0`, `E[J^2]<infinity` and `P(J>=2)>0`. The arms are `0,H,2H` at
the same exposure. In particular, `a=0`, `b=0`, or `u=0` are allowed.

The three laws uniquely identify `(a,b,c,u,v,G)`. Let `W_j=kappa_{2,j}`.
Compute

\[
r=\frac{\mu_1-\mu_2}{\mu_0-\mu_1},\quad
A_2=\frac{r(W_0-W_1)-(W_1-W_2)}{r(1-r)(1-r^2)},
\tag{12}
\]

\[
A_1=\frac{W_0-W_1-A_2(1-r^2)}{1-r},\quad
\pi=\frac{A_1}{A_1+2A_2},\quad p=\pi+(1-\pi)r.
\tag{13}
\]

Then replace `q` by `p` and `R_2` by `W_1/W_0` in (8) to obtain `x,V`,
and use (9)–(10). Finally,

\[
u=\pi\frac{-\log r}{H},\qquad
v=(1-\pi)\frac{-\log r}{H}.
\tag{14}
\]

Three means, three second factorial cumulants and the immediate zero
probability suffice for all finite parameters. Full baseline counts recover
`G`. Three arms are the minimum for a design guaranteed to identify the whole
class, even allowing two arbitrary nonnegative holiday times instead of a
baseline-containing pair.

**Proof of reconstruction.** By (5),

\[
\mu_j=b+(a+c\pi)m_1+c(1-\pi)m_1 r^j.
\]

The first mean drop is strictly positive. The ratio of successive drops is
therefore `r in (0,1)`. Expanding the square in (5) gives

\[
W_j=A_0+A_1r^j+A_2r^{2j},
\quad A_1=2c m_2\pi(1-\pi),\quad
A_2=c m_2(1-\pi)^2>0.
\tag{15}
\]

Write `D_0=W_0-W_1`, `D_1=W_1-W_2`. Then
`D_0=A_1(1-r)+A_2(1-r^2)` and
`D_1=A_1 r(1-r)+A_2 r^2(1-r^2)`.
Subtracting gives
`r D_0-D_1=A_2 r(1-r)(1-r^2)>0`, proving both the sign and the denominator
in (12). Equation (13) follows, with `A_1>=0`, including `A_1=0` when `u=0`.
Now `p in (0,1)`, `W_0>0`, and
`W_1/W_0=x+(1-x)p^2`. Also `mu_0-mu_1=c m_1(1-p)`.
The same positive-denominator argument as in Theorem 2 recovers `x,V,a,b,c,G`.
Equations (1) imply (14). This proves global, not merely local or generic,
injectivity on the stated domain. **End reconstruction proof.**

**Proof of minimality.** A baseline and a single positive holiday determine
at most its occupancy `p in (0,1)`: for every `pi' in [0,p)`, let
`r'=(p-pi')/(1-pi')`. Keeping `a,b,c,G` fixed and choosing rates corresponding
to `(pi',r')` preserves both laws. This argument works for every admitted
kernel and every positive reversible intensity.

For two strictly positive times `t_1<t_2`, there is also a general-kernel
obstruction. Take `a=0`, any `c>0`, `b>=0`, a positive clone law with finite
second moment and `m_2>0`, and an interior pair `0<pi<1,lambda>0`. Write
`p_i=pi+(1-pi)exp(-lambda t_i)`. The Jacobian of the map
`(pi,lambda)->(p_1,p_2)` has determinant

\[
(1-\pi)e^{-\lambda(t_1+t_2)}
 [t_1(e^{\lambda t_2}-1)-t_2(e^{\lambda t_1}-1)]>0.
\]

The inequality follows because `(exp(lambda t)-1)/t` is strictly increasing
for `t>0`: its derivative has numerator `exp(s)(s-1)+1`, positive for `s>0`
since that numerator is zero at zero and has derivative `s exp(s)>0`.
The inverse function theorem therefore allows nearby occupancies
`p_i'=k p_i`, for `k>1` sufficiently close to one, generated by different
valid rates. Set `epsilon=1/k`,

\[
s=1-G(1-\epsilon),\quad c'=cs,\quad
G'(z)=\frac{G(1-\epsilon+\epsilon z)-G(1-\epsilon)}s.
\tag{16}
\]

This is a positive-size binomially thinned clone law, conditioned to be
nonzero. It retains a finite second moment and positive probability of size
at least two. Direct substitution yields
`c'[G'(1-p_i'+p_i'z)-1]=c[G(1-p_i+p_i z)-1]` for both times.
The laws agree exactly. This supplies aliases for arbitrary such kernels and
all `c>0`, not an isolated finite-cell construction. It proves failure of a
two-arm design to identify the entire stated class; it does not assert that
every interior parameter with `a>0` is unidentifiable at two positive times.
Coincident times provide only one arm. Since three specified arms suffice,
the global design minimum is three. **End proof.**

A fully rational infinite-support instance of (16) is useful for checking:
at times `H,2H`, `a=b=0`, `c=1`, `pi=1/4`, `r=1/2`, `G(z)=z/(2-z)` and
`c'=20/21`, `pi'=11/136`, `r'=33/50`, `G'(z)=11z/(21-10z)` give identical
count laws. Their occupancies differ by the factor `11/10`. The immediate
arm distinguishes them.

### Theorem 4 — Exact failure strata and the dose/holiday design boundary

**Assumptions:** A1–A7, except for the explicitly stated kernel restrictions
or exposure extension in each part. Positive intensities in the lower-bound
constructions are interior parameter values; `c>0` and a nonzero holiday are
required in the positive reconstruction statements.

**4.1 — Which pure mechanisms coincide in an ordinary fluctuation assay?**
At `h=0`, pure stable pre-exposure clones with intensity `t` and pure reversible
clones with intensity `t` have exactly the same count distribution, for every
positive clone law and every `t>0`. No statistic or number of baseline
replicates separates them. A pure singleton-induced population is Poisson.
It differs from either nonzero clonal population when `P(J>=2)>0`; a finite
second moment is not needed for this distinction. If `J=1` surely, all three
pure mechanisms at equal intensity coincide. Mixtures at baseline identify
neither the `a:c` split nor the singleton portion of `a+c` versus `b` when
`G` is unknown.

**Proof.** Baseline (2) is
`log F_0=(a+c)(G-1)+b(z-1)`. For a pure clonal source, the log-PGF has a
strictly positive coefficient of `z^k` for some `k>=2` exactly when the clone
law has such a coefficient. A Poisson law has no such coefficient. Analytic
uniqueness on `0<=z<1` gives the assertion without using moments. Moving
stable versus reversible event intensity while preserving `a+c` changes no
law. For the singleton ambiguity, write `nu_k=(a+c)g_k`; only `nu_1+b`
and `nu_k,k>=2` appear. On an interior point with `nu_1,b>0`, a small transfer
between `nu_1` and `b`, followed by renormalizing `G`, preserves the whole
baseline law. **End proof.**

**4.2 — Unknown doubleton kernels require three one-way arms.** Suppose
`G(z)=(1-w)z+wz^2`, `0<w<1`, `a,b,c>0` and a positive finite one-way rate.
For **any two distinct nonnegative holiday times**, there is a continuum of
different nearby parameters giving the same two full laws. Three arms
`0,H,2H` suffice, by Theorem 3 with `u=0`. Thus the requirement of a possible
clone of size at least three in the two-arm theorem cannot simply be dropped.

**Proof.** At times `t_i`, set `p_i=exp(-v t_i)` and
`X=aw`, `Y=cw`, `U=a(1-w)+b`, `V=c(1-w)`. The two positive log coefficients
are

\[
\ell_{2,i}=X+Yp_i^2,\qquad
\ell_{1,i}=U+Vp_i+2Yp_i(1-p_i).
\tag{17}
\]

For any nearby `v'!=v`, its distinct `p'_1,p'_2` permit the following solves:

\[
Y'=\frac{\ell_{2,1}-\ell_{2,2}}{(p'_1)^2-(p'_2)^2},\quad
X'=\ell_{2,1}-Y'(p'_1)^2,
\]

\[
B_i=\ell_{1,i}-2Y'p'_i(1-p'_i),\quad
V'=\frac{B_1-B_2}{p'_1-p'_2},\quad U'=B_1-V'p'_1.
\]

Recover `w'=Y'/(V'+Y')`, `c'=V'+Y'`, `a'=X'/w'`,
`b'=U'-a'(1-w')`. At the unperturbed point these are the original positive
parameters, and all denominators are nonzero. By continuity, sufficiently
small perturbations keep them admissible and `0<w'<1`. They preserve every
positive log coefficient; normalization determines the constant coefficient,
so both entire laws agree. This is an open-parameter family of aliases.
**End proof.**

For exact arithmetic, the baseline/holiday pair

\[
(a,b,c,q,G)=(1,1,1,1/2,(z+z^2)/2)
\]

and

\[
(a',b',c',q',G')=(185/192,4/3,45/64,1/3,(2z+3z^2)/5)
\]

both give `log F_0=2(z-1)+(z^2-1)` and
`log F_1=2(z-1)+(5/8)(z^2-1)`. A third arm separates them. This small
certificate checks the general proof; it is not the main result.

**4.3 — A known kernel removes one nuisance dimension.** In the one-way
model, if `G` is known, `m_2>0` and `E[J^2]<infinity`, two arms identify
`a,b,c,q` even for a doubleton kernel:

\[
q=\frac{m_1}{m_2}\frac{W_0-W_1}{\mu_0-\mu_1}-1,\quad
A=W_0/m_2,\quad c=\frac{\mu_0-\mu_1}{m_1(1-q)},
\quad a=A-c,\quad b=\mu_0-m_1 A.
\tag{18}
\]

**Proof.** The numerator and denominator in the ratio are
`c m_2(1-q^2)` and `c m_1(1-q)`. Cancelling gives `1+q` and the remaining
identities follow from (5). Every denominator is positive. **End proof.**

**4.4 — Singleton clones and a calibrated exposure contrast.** If `G(z)=z`,
all count laws are Poisson with means `a+b+cp(h)`. Even an infinite sequence
of holiday times at one exposure cannot separate `a` and `b`.

For an explicit extension, replace the induced intensity by `b x`, where
`x>=0` is a **known effective exposure covariate**, with all other parameters
fixed across exposure arms. This is not an uncalibrated assertion that
induction is linear in physical drug concentration. Under one-way switching,
use arms `(x,0),(x,H),(x,2H),(x',0)`, with `x'!=x`. If their means are
`s_0,s_1,s_2,t_0`, then

\[
q=\frac{s_1-s_2}{s_0-s_1},\quad c=\frac{s_0-s_1}{1-q},\quad
b=\frac{t_0-s_0}{x'-x},\quad a=s_0-bx-c.
\tag{19}
\]

Four arms are necessary for a design that identifies all interior
`(a,b,c,v)` in this singleton exposure model, even with arbitrary choices
of exposure and holiday times.

**Proof.** The assertions and inverse follow from
`E Y(x,h)=a+bx+c exp(-vh)`. Each entire count law is determined by its one
mean. With at most three arms the observation map is a continuous map from
an open subset of `R^4` into `R^m`, `m<=3`, and cannot be injective.
For completeness, this standard consequence of invariance of domain follows
by padding the output with zero coordinates to obtain a continuous injection
into `R^4`: invariance of domain would make its image open, whereas the
image lies in a proper linear subspace with empty interior. The four-arm
inverse proves sufficiency. This lower bound concerns fixed designs and
structural laws, not adaptive measurements with extra observations.
**End proof.**

For two-way singleton switching, even allowing arbitrary `x` and all holiday
times fails: write

\[
E Y(x,h)=A+bx+C e^{-\lambda h},\quad
A=a+c\pi,\quad C=c(1-\pi).
\tag{20}
\]

For a nearby `pi'` set `c'=C/(1-pi')`, `a'=A-c'pi'` and keep `b,lambda`
fixed. At interior positive parameters these remain admissible while changing
`a,c,u,v`. All means and hence all laws agree. Extra exposure levels alone
do not remove the tolerant-equilibrium ambiguity.

### Theorem 5 — Robustness obstructions and what counts cannot establish

**Assumptions:** Each part is the explicitly stated relaxation of A1–A7.
All other independence and count assumptions are retained. Each construction
holds for every holiday time, not just a particular sampling grid.

**5.1 — Different unknown clone laws.** Allow `G_D` and `G_P` to differ.
When `a>0`, `P(J_D=1)>0`, choose `0<delta<a P(J_D=1)`. Set

\[
a'=a-\delta,\quad b'=b+\delta,\quad
G_D'(z)=\frac{aG_D(z)-\delta z}{a-\delta}.
\tag{21}
\]

This is a valid positive clone law, and
`a'(G_D'-1)+b'(z-1)=a(G_D-1)+b(z-1)`. Keeping the reversible component
fixed gives identical laws for all holidays. Shared `G` is essential to the
unknown-kernel positive result; fitting it well does not verify its biology.

**5.2 — Coordinated switching inside a clone.** If a whole reversible clone
is retained with probability `p(h)` and its cells switch together, the
log-PGF is `[a+c p(h)](G-1)+b(z-1)`. With two-way switching, this depends
on `a,c,pi` only through `a+c pi` and `c(1-pi)`. The transformation after
(20) leaves every full count law unchanged, including for arbitrarily large
or unbounded clone sizes. Independent cell-level switching, rather than
merely reversible population behavior, supplies the `p^2` term in (15).

**5.3 — Induced clonal amplification.** If an induced event has mark law
`G`, its stable contribution is `b(G-1)`. It is interchangeable with the
pre-exposure stable contribution `a(G-1)` at every holiday. More generally,
arbitrary unknown induced mark laws admit overlapping observational classes.
The single-founder acquisition window is therefore an essential assumption
for a timing interpretation.

**5.4 — Unknown common detection efficiency.** Suppose every competent cell
is independently detected with an unknown common `epsilon in (0,1]` in
every arm. Let `s=1-G(1-epsilon)>0` and define `G*` by the zero-truncated
thinning in (16). The observed laws are exactly (2) with

\[
a^*=sa,\qquad c^*=sc,\qquad b^*=\epsilon b,\qquad G^*=G'.
\tag{22}
\]

The holiday occupancies are unchanged, since two independent binomial
thinnings commute. Thus the theorem identifies effective detectable-clone
intensities and switching rates, not uncalibrated biological event rates.
The stable:reversible event-intensity ratio is preserved. To see that this is
a real nonidentifiability, the original model with imperfect detection and
the starred model with perfect detection give the same observations. With
different unknown efficiencies by arm, even retention can be confounded:
in the pure reversible submodel only the products `epsilon_j p_j` enter
`G(1-epsilon_j p_j+epsilon_j p_j z)`.

**5.5 — Absent components, permanent memory, and invisible events.** If
`c=0`, neither switching rate is observable; with unknown `G`, the stable
singleton transfer (21) also confounds `a,b`. If `v=0` and every reversible
cell starts in `T`, all reversible cells stay in `T`, so `a,c` collapse into
their sum. If zero-size marks are allowed, writing
`G=t+(1-t)G_+` shows that `a(G-1)=a(1-t)(G_+-1)`; event intensity and zero
marks cannot be separated. These strata are excluded from the full inverse,
not fixed by additional replication. If `a=c=0`, there is no clone law to
identify. The primary model allows `a=0` or `b=0` separately with `c>0`.

**5.6 — Drug causation needs an exposure control.** Extend `b` to
`b(x)=b_0+eta x` with unknown nonnegative background `b_0` and excess
coefficient `eta`, known exposure covariate `x`, and all other features
calibrated as invariant. At a single exposure only `b_0+eta x` is observable.
Two distinct exposures are necessary and sufficient to separate these two
numbers *once their respective total acquisition intensities are identified*:
`eta=(b(x')-b(x))/(x'-x)` and `b_0=b(x)-eta x`. A sufficient combined design
is the three holiday arms of Theorem 3 at `x`, followed by an immediate
arm at `x'`; that arm gives `b(x')=-log F_{x',0}(0)-(a+c)`. This is a
sufficient four-arm design, not a claim of total-arm optimality among all
joint dose/time designs. If dose changes `a,c,G` or detection, those effects
must also be modeled or independently calibrated. Counting post-exposure
events alone does not establish that exposure caused them.

**Proof of all parts.** Each displayed alternative has nonnegative intensity
and a normalized mark law under its listed strict inequalities. Substituting
(21) cancels the added singleton terms. Whole-clone retention replaces the
thinned PGF by `1-p+pG`, giving 5.2. Addition of independent compound-Poisson
stable sources adds their intensities when their mark laws agree, proving
5.3. For 5.4, substitute `1-epsilon+epsilon z` for `z` in (2), remove the
zero mark by (16), and note
`1-epsilon+epsilon(1-p+p z)=1-epsilon p+epsilon p z`.
This proves both the equivalence and the unchanged switching parameters.
The zero contributions and scalar-factor identities prove 5.5. In 5.6 the
single-exposure exponent contains only `b_0+eta x`; a small transfer preserving
that sum yields aliases at interior values. The two-by-two linear system at
distinct exposures has determinant `x'-x!=0`, giving the stated inverse and
the sufficient combined design. **End proof.**

## 5. Finite cultures: likelihoods, exact minima, and impossibility bounds

The full-count likelihood for independent observations `y_{j,i}` is
`product_{j,i} exp(-L_j) w_{j,y_{j,i}}`, with coefficients (3)–(4).
Likelihood fitting and uncertainty for an arbitrary unknown kernel are not
implemented here. The exact verifier is not presented as a statistically
regularized estimator for real data. The following tests give finite-sample
statements without pretending structural inverses are stable estimators.

### Theorem 6 — An exactly computable minimum for a prespecified test

**Assumptions:** Independent cultures, two simple prespecified hypotheses,
known zero-count probabilities `0<p_0<p_1<1` in a specified arm, type-I target
`0<alpha<1`, and power target `alpha<1-beta<1`. No parameters are estimated
from the cultures used in this test. Restrict the observation to the indicator
of zero colonies in each culture. For the **exact finite enumeration** claim,
`p_0,p_1,alpha,beta` are rational; the likelihood-ratio formulas themselves
hold for real probabilities. Arbitrary real inputs require certified exact
comparisons and are not claimed to have an effective exact decision procedure.

Let `Z` be the number of zero-count cultures among `n` cultures. Write

\[
f_{i,n}(k)={n\choose k}p_i^k(1-p_i)^{n-k},\quad
T_{i,n}(k)=\sum_{l=k}^{n}f_{i,n}(l).
\tag{23}
\]

The most powerful randomized size-`alpha` test rejects for `Z>k` and with
probability `gamma` for `Z=k`, where

\[
T_{0,n}(k+1)\le\alpha<T_{0,n}(k),\qquad
\gamma=\frac{\alpha-T_{0,n}(k+1)}{f_{0,n}(k)}.
\tag{24}
\]

Its power is `T_{1,n}(k+1)+gamma f_{1,n}(k)`. The least `n>=1` for which
this exceeds or equals `1-beta` is the exact minimum among **all randomized
tests of these zero indicators**. Exhaustively checking each smaller `n`
certifies minimality. For deterministic upper-threshold tests, put `gamma=0`
with the largest admissible upper tail and perform the same finite search.
The latter minimum is only over upper-threshold tests; arbitrary deterministic
subsets of outcomes can use a significance budget differently.

**Proof.** The likelihood ratio is
`(p_1/p_0)^Z ((1-p_1)/(1-p_0))^(n-Z)`, strictly increasing in `Z`.
For any test `phi` and the threshold/randomized test `phi*`, choose the
boundary likelihood ratio `l`. Pointwise,
`(phi*-phi)(f_1-l f_0)>=0`; on the boundary this product is zero.
Summing and using `E_0 phi<=E_0 phi*=alpha` gives
`E_1 phi*>=E_1 phi`. The same argument applies to the full binary vector
because its likelihood ratio depends only on `Z`. This proves optimality
including allowed randomization. Within deterministic upper tails, power is
maximized by the least threshold satisfying size. Finally, because `p_1>p_0`,
a threshold between their limiting empirical proportions has both errors
converging to zero, so a finite successful `n` exists for the stated targets.
Enumeration proves the least one in each class; no monotonicity assumption
about deterministic threshold power is needed. **End proof.**

#### Worked exact power certificate with an unbounded clone law

Take `G(z)=z/(2-z)`, the geometric law on positive integers with mean two.
Compare these *simple hypotheses*:

\[
H_D:(a,b,c)=(2\log2,0,0),\qquad
H_P:(a,b,c)=(0,0,2\log2).
\tag{25}
\]

Their entire immediate distributions agree, for an unbounded population
count. At a holiday with `p=1/3`, `G(1-p)=G(2/3)=1/2`. Equation (2) gives
`P_D(Y=0)=1/4` and `P_P(Y=0)=1/2`. For `alpha=1/20`, target power `4/5`:

| Test class | Exact minimum | Rule | Size | Power |
|---|---:|---|---:|---:|
| Deterministic upper threshold on zero counts | 26 holiday cultures | Reject `H_D` if `Z>=11` | `22565567197793/562949953421312` = 0.0400844996 | `1754331/2097152` = 0.8365302086 |
| Randomized Neyman–Pearson test of zero indicators | 23 holiday cultures | Reject if `Z>9`; if `Z=9`, reject with probability `1622209696738/9771486092775` | Exactly 0.05 | `40811361573509/50153065021440` = 0.8137361407 |

Both code paths exhaustively rule out every smaller positive sample size in
the relevant class using `Fraction`. These are not full-count optimal sample
sizes, not estimates of error after fitting nuisance parameters, and not a
universal 23- or 26-culture prescription. Baseline cultures used to check the
assay assumptions are additional to these specified holiday-test counts.

For a matched-mean induced comparator take `H_I: b=4 log 2,a=c=0`.
Its mean equals those in (25) at baseline, while its zero probability is
`1/16`. Thus the same exact program computes baseline zero-indicator
discrimination of this prespecified `I` against `D` or `P`; at holiday their
zero probabilities are `1/16`, `1/4`, `1/2`. The direction of the null must
be chosen explicitly; these pairwise simple-hypothesis calculations are not
a composite three-way classification guarantee.

### Theorem 7 — No uniform finite sample size; a bounded-statistic guarantee

**Assumptions for impossibility:** A1–A6, a finite collection of arms and
finite replicate counts `n_j`, with no lower bound on a reversible intensity
or on the amount of reversion. For the constructive bound, the two compared
count laws and the separating statistic are fixed in advance.

First, compare a model with `c=0` to the same model with an added reversible
component of intensity `c>0`. The joint-law distance is bounded by

\[
\operatorname{TV}(P,Q)\le
1-\exp\left[-c\sum_j n_j\{1-G(1-p_j)\}\right].
\tag{26}
\]

Second, compare intensity `c` of permanent clones with the same intensity of
reversible clones and keep their baseline laws identical. The corresponding
bound is

\[
\operatorname{TV}(P,Q)\le
1-\exp\left[-c\sum_j n_j\{1-G(p_j)\}\right].
\tag{27}
\]

For any randomized size-`alpha` test its power against the respective
alternative is at most `alpha+TV(P,Q)`. Hence no finite uniform number of
cultures guarantees power greater than `alpha` by a fixed amount over the
unrestricted class. More replication does not resolve exactly equal laws.

Conversely, for any two distinct specified arm laws there is a rational
`z in [0,1)` with `delta=|F_1(z)-F_0(z)|>0`. Use `X=z^Y`, interpreting
`0^0=1`, and reject according to which side of the midpoint
`(F_0(z)+F_1(z))/2` the sample mean lies. Both error probabilities are at
most `exp(-n delta^2/2)`. Consequently

\[
n\ge\left\lceil\frac{2}{\delta^2}
\max\left\{\log(1/\alpha),\log(1/\beta)\right\}\right\rceil
\tag{28}
\]

suffices for size at most `alpha` and power at least `1-beta`. This is a
conservative sufficient count, not an optimal minimum. It needs no raw
moment finiteness and extends to a finite prespecified list of alternatives
by a union bound with allocated error budgets.

**Proof.** For (26), couple the stable and induced pools identically. A
reversible event contributes a nonzero mark with probability
`1-G(1-p_j)`. Independent Poisson marking implies no added observed event
in any culture has probability equal to the exponential in (26). On that
event the datasets coincide, proving the coupling bound. For (27), couple
each permanent clone to a reversible clone of the same size, thinning the
latter. They agree if no cell in any such clone is lost. A clone agrees with
probability `E[p_j^J]=G(p_j)`, so Poisson marking gives (27).
The distinction between `G(1-p_j)` and `G(p_j)` is essential. The first
bound tends to zero as `c` tends to zero; the second tends to zero as all
`p_j` tend to one, by dominated convergence. Total variation bounds the
difference of expectations of any `[0,1]`-valued test, proving the power
limitation.

For (28), PGFs are continuous on `[0,1)` and uniquely determine count laws.
If they agreed at every rational there, continuity would make them agree
throughout the interval, and analytic coefficient uniqueness would make the
laws equal. Thus a rational separator exists. The observations `z^Y` lie in
`[0,1]`. The Hoeffding bound for a one-sided deviation of size `delta/2`
is `exp(-2n(delta/2)^2)`, giving both errors and (28). This standard bound
can be obtained by the bounded-variable exponential-moment inequality
`E exp(t(X-EX))<=exp(t^2/8)` and Markov's inequality, minimized at
`t=4(delta/2)`. To justify that inequality, the second derivative of the
log-moment-generating function is the variance under exponential tilting.
The tilted variable still lies in `[0,1]`, so its variance is at most `1/4`
(center at `1/2` and subtract the squared mean displacement). The centered
log-moment-generating function and its first derivative vanish at zero;
integrating the second-derivative bound twice gives `t^2/8` for either sign
of `t`. No assumption on a raw moment of `Y` is used. **End proof.**

The power calculations are for specified separated alternatives. They do not
prove uniform estimation accuracy near `c=0`, `v=0`, `m_2=0`, `m_3=0`, or
very long/short holidays. For example, the denominator in (12) tends to zero
as `r` approaches zero or one. Finite second moments suffice for structural
Theorem 3 but need not give finite variances for empirical variance estimates;
Theorem 2's empirical third cumulant can likewise be extremely unstable.
Bounded PGF statistics and full count likelihoods are distinct from simply
plugging small-sample raw moments into the inverse.

## 6. A worked experimental-design worksheet

This is a nonclinical **assay-design calculation**, not a validated wet-lab
protocol. The unit of observation is a culture's number of competent
single-cell colony founders. A laboratory must first establish that its
readout and holiday actually satisfy A3–A5. Failure of these checks calls for
a different model, not more cultures fitted to this one.

1. Grow independent cultures under a common preparation. Use a short,
   standardized exposure and acquisition window. Verify that late stable
   acquisitions have not already amplified and that the stable and reversible
   clone-size kernels can reasonably be shared. Lineage measurements can
   challenge this assumption; the count-only theorem does not validate it.
2. Randomize cultures into immediate, `H`, and `2H` readout arms at the same
   exposure. Remove the unmodeled external sensitive pool. During the holiday,
   check that cell proliferation, death and new acquisition are negligible
   for the modeled founder population, and that switching is cellular rather
   than coordinated within a clone. Use matched controls for detection and
   state changes during colony outgrowth. Count all cultures, including zeros.
3. Record each complete count, the original culture identifier, exposure,
   holiday duration, plating fraction, and censoring/detection information.
   Treat technical aliquots as dependent unless a joint model is supplied.
   Do not replace individual counts with only their average or discard jackpots.
4. In a pilot, estimate the three means and `Var(Y)-E[Y]`. The inequalities
   `mu_0>mu_1>mu_2`, `0<r<1`, `A_2>0`, `A_1>=0`, `0<=x<1` and
   `a,b>=0,c>0` are necessary model checks, allowing sampling uncertainty.
   Check the full PGFs against (2), not just these seven summaries. An extra
   `3H` arm supplies a direct out-of-design falsification of the predicted
   retention curve. These checks can reject the model; passing does not prove it.
5. Prespecify an effect size and alternative for confirmatory power. For the
   synthetic alternatives (25), the exact holiday-only zero test needs 26
   cultures for a deterministic upper-threshold rule at the stated errors.
   Actual pilot uncertainty and nuisance fitting require their own design
   analysis. If causation by drug is the target, add an exposure control and
   use the conditions in 5.6 before interpreting a post-exposure component as
   a drug-caused excess.

### Arithmetic example for parameter reconstruction

Take the synthetic law

\[
P(J=1)=P(J=2)=1/4,\quad P(J=4)=1/2,
\quad a=3/5,\ b=1/3,\ c=4/5.
\]

Let `pi=1/4`, `r=1/2`. Then `m_1=11/4`, `m_2=13/2`, and occupancies
are `1,5/8,7/16`. The observable summaries are

| Holiday | Mean | Second factorial cumulant |
|---|---:|---:|
| `0` | `251/60` | `91/10` |
| `H` | `403/120` | `949/160` |
| `2H` | `707/240` | `3133/640` |

The immediate negative log-zero probability is `26/15`. Mean drops are
`33/40` and `33/80`, so (12) returns `r=1/2`. Equations (12)–(13) give
`A_2=117/40`, `A_1=39/20`, and `pi=1/4`. Equations (8)–(10) return
`x=3/7`, `V=77/20`, `b=1/3`, `a=3/5`, `c=4/5` and the original kernel.
The hazards are `u=(log 2)/(4H)` and `v=3(log 2)/(4H)`. These numbers are
synthetic and dimensionless until a laboratory chooses and calibrates `H`.

## 7. Verification, falsification attempts, and reproducibility

Only this directory and
`tests/test_research_induced_vs_selected_resistance.py` are added. No dataset,
patient record, external numerical dependency, workflow change or biological
result is needed.

### Exact arithmetic and independent simulation

`verifier.py` uses `fractions.Fraction` for clone thinning, log-PGF
coefficients, factorial cumulants, inverses, alias certificates and binomial
power calculations. It represents count probabilities as
`exp(-L) * rational_weight`; an exponential of a nonzero rational is not
misrepresented as a rational number. Rational alternating Taylor bounds after
range reduction enclose `exp(-L)`. Nonnegative multiplication and repeated
squaring preserve their order. The PMF prefix therefore has certified
intervals, and subtracting its enclosed mass from one gives an explicit tail
interval. The finite coefficient recurrence is separately checked against
expanding the product `product_k exp(ell_k z^k)` and a formal-log inverse.
`Decimal` is used for readable output and an independent high-precision check,
not as a proof of exact equality.

The two-arm rational inverse requires a square rational discriminant for
its exact square-root routine. It succeeds on the generated rational model
instances tested here. Empirical observations need not have such a
discriminant or even lie in the model; use the real formulas in the theorem
and an appropriate constrained statistical analysis for empirical fitting.

`simulation.py` is a separate random engine. It imports no verifier and
evaluates no analytic count PGF. It samples Poisson arrival clocks, independent
clone marks, and independent Bernoulli retention, or directly simulates each
cell's `T->S` and `S->T` exponential clocks. It supports arbitrary finite
positive clone laws and an unbounded geometric law. Its pseudorandom draws
and waiting times use floating-point arithmetic, deliberately independent of
the exact certificate code. The default fixed seed is `20261117`.

The test independently checks coefficient expansions, parameter inverses,
boundary cases, exact aliases, exponential enclosures, finite-test optimality
and event simulations. Simulated PGFs use bounded statistics with explicit
tolerances. Simulation is a consistency/falsification check of the model and
code, never a proof of a universal theorem. No noisy empirical third-moment
agreement is used as a proof or a narrow CI gate.

### Hostile-referee falsification log

| Attempted stronger claim or failure mode | Outcome and correction |
|---|---|
| “No previous structural-identifiability theorem exists.” | Falsified by [7], [12], [14]. The research question was narrowed before the result was developed. |
| “Two endpoints identify both switching rates.” | Falsified by the continuum of `(pi,r)` with fixed `p(H)`. Theorem 2 identifies one-way loss or a single occupancy; Theorem 3 uses a third arm. |
| “Four factorial moments are required for two arms.” | Unnecessary. The third-moment ratio is invertible because the ratio at order zero is fixed at one; (8) reduces the requirement to three. |
| “Two arms suffice for every non-singleton unknown kernel.” | Falsified by the open doubleton alias family (17). The exact rational example is a regression certificate, not the main theorem. |
| “A lower bound only concerns baseline-containing designs.” | Strengthened. Theorem 3 supplies a two-positive-time ambiguity for general kernels; Theorem 4.2 handles arbitrary two times for doubletons. The scope of uniform versus pointwise minimality is explicit. |
| “A shared unknown kernel is merely convenient.” | Falsified by the all-holiday singleton-transfer alias (21). It is a necessary modeling restriction for the stated positive result. |
| “Any reversible population obeys the same inverse.” | Falsified by coordinated whole-clone switching, whose holiday plateau is confounded at all times. |
| “A post-exposure stable component proves drug-induced mutation.” | Falsified as a causal interpretation without control and molecular information. A7 and 5.6 separate timing, drug causation and molecular mechanism. |
| “A long holiday makes the problem automatically easy.” | False as a conditioning claim: the denominator in (12) also vanishes as `r->0`. No optimal holiday duration is asserted. |
| “Unknown plating efficiency is harmless.” | Only a common efficiency preserves occupancy inference; event rates are rescaled and arm-specific efficiency confounds occupancy. |
| “Induced cells can freely grow before readout.” | Falsified by 5.3: if their mark law matches stable clones, their origins become indistinguishable. |
| “An absent reversible component only loses its switching rate.” | Incomplete: unknown stable singleton mass also confounds `a,b`. The full inverse explicitly excludes `c=0`. |
| “Infinite-growth LD laws automatically satisfy the moment assumptions.” | False. The positive theorems require finite stated moments and do not silently exchange an infinite-growth limit with differentiation. |
| “Finite structural summaries mean a few cultures suffice.” | False. Near-boundary coupling bounds rule out a uniform finite guarantee, and high empirical cumulants can be unstable. |
| “26 is the minimum for any possible test.” | False. It is the deterministic upper-threshold minimum for specified zero probabilities; a randomized zero-indicator LR test needs 23. Full-count optimality is not claimed. |
| “Rounded probability calculations prove the reported minimum.” | Avoided. Fractions enumerate every smaller sample size; simulation is checked separately. |
| “Exact finite enumeration works for unspecified arbitrary real inputs.” | Too strong as an algorithmic statement. Theorem 6 restricts exact enumeration to rational probabilities and targets; the real-valued likelihood-ratio formula remains valid. |
| “Valid higher-cumulant ratios alone make an inverse input admissible.” | Falsified by a rising-mean input with decreasing higher cumulants. The verifier explicitly rejects nonpositive mean loss; reconstructed kernels must still satisfy the full model. |

No unproved biological or global model-extension claim is promoted to a
theorem. Remaining research questions include unequal clone kernels with
additional lineage observations, birth/death during holiday, heterogeneous
switching rates, count censoring, and useful composite-hypothesis minimax
sample bounds. These are open directions within this submission, not claims
that the literature has never studied them. No theorem here implies a dosing
benefit or supports an individual treatment decision.

### Rerun from the repository root

```bash
python research/induced-vs-selected-resistance/verifier.py --self-check
python research/induced-vs-selected-resistance/simulation.py --seed 20261117 --cultures 20000
python -m unittest discover -s tests -p test_research_induced_vs_selected_resistance.py -v
```

The new test uses only the standard library and is designed to run below
60 seconds. Repository integration gates remain unchanged:

```bash
python scripts/validate.py
python -m unittest discover -s tests -v
node --test tests/*.test.*
```

The exact checks are finite regression checks of general proofs, not exhaustive
proof assistants. The unit test imports the files by their paths, so no package
or requirements changes are needed. CLI output is explicitly synthetic and
nonclinical. Keep generated console output as a rerun artifact if desired;
the committed source contains the complete reproducible certificate.

## References

Only citations and our original derivations/code are included; no source
datasets or paper text are redistributed. DOI links identify the sources.
Bibliographic and substantive descriptions above were checked against primary
papers or their author-hosted versions.

1. Luria SE, Delbrück M (1943). *Mutations of Bacteria from Virus Sensitivity
   to Virus Resistance*. Genetics 28:491–511.
   <https://doi.org/10.1093/genetics/28.6.491>.
2. Lea DE, Coulson CA (1949). *The distribution of the numbers of mutants in
   bacterial populations*. Journal of Genetics 49:264–285.
   <https://doi.org/10.1007/BF02986080>.
3. Ma WT, Sandri GVH, Sarkar S (1992). *Analysis of the Luria–Delbrück
   distribution using discrete convolution powers*. Journal of Applied
   Probability 29:255–267. <https://doi.org/10.2307/3214564>.
   Related exact-distribution paper: Sarkar S, Ma WT, Sandri GVH (1992),
   Genetica 85:173–179. <https://doi.org/10.1007/BF00120324>.
4. Lang GI, Murray AW (2008). *Estimating the Per-Base-Pair Mutation Rate in
   the Yeast Saccharomyces cerevisiae*. Genetics 178:67–82.
   <https://doi.org/10.1534/genetics.107.071506>;
   [author PDF](https://www.glanglab.com/papers/Genetics2008.pdf).
5. Pisco AO et al. (2013). *Non-Darwinian dynamics in therapy-induced cancer
   drug resistance*. Nature Communications 4:2467.
   <https://doi.org/10.1038/ncomms3467>.
6. Shaffer SM et al. (2017). *Rare cell variability and drug-induced
   reprogramming as a mode of cancer drug resistance*. Nature 546:431–435.
   <https://doi.org/10.1038/nature22794>.
7. Greene JM, Gevertz JL, Sontag ED (2019). *Mathematical Approach to
   Differentiate Spontaneous and Induced Evolution to Drug Resistance During
   Cancer Treatment*. JCO Clinical Cancer Informatics 3:1–20.
   <https://doi.org/10.1200/CCI.18.00087>;
   [author PDF, including identifiability appendix](https://www.sontaglab.org/FTPDIR/2019_greene_gevertz_sontag_asco.pdf).
8. Bokes P, Singh A (2021). *A modified fluctuation test for elucidating drug
   resistance in microbial and cancer cells*. European Journal of Control
   62:130–135. <https://doi.org/10.1016/j.ejcon.2021.06.018>.
9. Kuosmanen T, Cairns J, Noble R, Beerenwinkel N, Mononen T, Mustonen V
   (2021). *Drug-induced resistance evolution necessitates less aggressive
   treatment*. PLOS Computational Biology 17:e1009418.
   <https://doi.org/10.1371/journal.pcbi.1009418>.
10. Russo M et al. (2022). *A modified fluctuation-test framework characterizes
    the population dynamics and mutation rate of colorectal cancer persister
    cells*. Nature Genetics 54:976–984.
    <https://doi.org/10.1038/s41588-022-01105-z>.
11. Saint-Antoine MM, Grima R, Singh A (2022). *Fluctuation-based approaches
    to infer kinetics of cell-state switching*. IEEE CDC, 3878–3883.
    <https://doi.org/10.1109/CDC51059.2022.9992907>.
12. Gunnarsson EB, Foo J, Leder K (2023). *Statistical inference of the rates
    of cell proliferation and phenotypic switching in cancer*. Journal of
    Theoretical Biology 568:111497.
    <https://doi.org/10.1016/j.jtbi.2023.111497>;
    <https://arxiv.org/abs/2306.08096>.
13. Bokes P, Hlubinová A, Singh A (2023). *Reversible Transitions in a
    Fluctuation Assay Modify the Tail of Luria–Delbrück Distribution*. Axioms
    12:249. <https://doi.org/10.3390/axioms12030249>.
14. Gevertz JL, Greene JM, Prosperi S, Comandante-Lou N, Sontag ED (2025).
    *Understanding therapeutic tolerance through a mathematical model of
    drug-induced resistance*. npj Systems Biology and Applications 11:30.
    <https://doi.org/10.1038/s41540-025-00511-3>.
15. Browning AP et al. (2025). *Identifiability of phenotypic adaptation from
    low-cell-count experiments and a stochastic model*. PLOS Computational
    Biology 21:e1013202. <https://doi.org/10.1371/journal.pcbi.1013202>.
16. Wu C, Gunnarsson EB, Foo J, Leder K (2025). *A statistical framework for
    detecting therapy-induced resistance from drug screens*. npj Systems
    Biology and Applications 11:88.
    <https://doi.org/10.1038/s41540-025-00560-8>.
17. Puig P, Valero J (2007). *Characterization of count data distributions
    involving additivity and binomial subsampling*. Bernoulli 13:544–555.
    <https://doi.org/10.3150/07-BEJ6021>;
    <https://arxiv.org/abs/0708.4177>.
18. Gugushvili S, Mariucci E, van der Meulen F (2020). *Decompounding discrete
    distributions: A non-parametric Bayesian approach*. Scandinavian Journal
    of Statistics 47:464–492. <https://doi.org/10.1111/sjos.12413>;
    <https://arxiv.org/abs/1903.11142>.
19. Potts D, Tasche M (2013). *Parameter estimation for nonincreasing
    exponential sums by Prony-like methods*. Linear Algebra and its
    Applications 439:1024–1039. <https://doi.org/10.1016/j.laa.2012.10.036>.
20. Whiting FJH et al. (2025). *Quantitative measurement of phenotype dynamics
    during cancer drug resistance evolution using genetic barcoding*. Nature
    Communications 16:5282. <https://doi.org/10.1038/s41467-025-59479-7>.
21. Iyer A, Alva A, Granada AE, Chakrabarti S (2025). *Inheritable cell-states
    shape drug-persister correlations and population dynamics in cancer
    cells*. PLOS Computational Biology 21:e1013446.
    <https://doi.org/10.1371/journal.pcbi.1013446>.
