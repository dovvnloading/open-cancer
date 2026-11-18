# ctDNA surveillance: growth heterogeneity, geometric coverage, and sampling-error bounds

**Status: unreviewed, nonclinical computational mathematics.** These are
deductive statements about specified mathematical models, not established
tumour biology, externally reviewed research, or treatment or monitoring
advice for any patient. Synthetic simulations and passing CI do not establish
clinical validity, survival benefit, or bibliographic priority.

Literature search and note: **17 November 2026**. Base repository revision:
`d257ac603d5e57a30cbce222446fb799a53a17c0` on `main`.

## 1. Literature first: what is already known, and the explicit pivot

**The broad scheduling question is already established.** In particular,
post-treatment surveillance with false negatives and a detection-delay
objective predates ctDNA, and ctDNA-specific computational scheduling appeared
before this note. We therefore pivot from proposing the first surveillance
optimizer to a narrower question:

> For exponential regrowth, what exact fixed-budget schedule hits a
> growth-linked detectable-to-clinical burden window most often, and how much
> optimality can that schedule lose under random residual burden and finite
> blood sampling?

The candidate contribution is the combination of (i) an exact reduction to
equal-length intervals in **log growth rate**, (ii) global optimal dates for
log-normal and gamma rates, with the correct clinical-horizon boundary,
(iii) different same-mean heterogeneity effects in those two families and an
explicit heterogeneity maximizing the value of another test, and (iv) a
uniform bound transferring the ideal optimum to random-burden, noisy-assay
surveillance. The measure-rearrangement, submodularity, concentration, and
equal-spacing ingredients are classical; their originality is **not** claimed.

### Prior-work ledger

Links identify authoritative publications or authors' manuscripts. No source
PDF, patient record, or cohort-level individual data is redistributed here.
“Abstract” means that the full article was not available in this search;
absence of a theorem from its abstract cannot establish priority.

| Primary source | Established contribution and implication for this note | Access used |
| --- | --- | --- |
| [Zelen, 1993, *Biometrika* 80:279–293](https://doi.org/10.1093/biomet/80.2.279) | General optimal examination spacing, preclinical sojourn distributions, imperfect sensitivity, and links to surveillance and equipment inspection. Equal-spacing conditions are already studied. | Publisher abstract |
| [Parmigiani, 1993, *On Optimal Screening Ages*](https://doi.org/10.1080/01621459.1993.10476315); [1996, *Optimal Scheduling of Fallible Inspections*](https://doi.org/10.1287/opre.44.2.360) | Decision-theoretic ages and continuous-time fallible inspection optimization, with dynamic programming and structural properties. Neither decision-theoretic scheduling nor imperfect inspection is new here. | Publisher abstracts |
| [Lee and Zelen, 1998](https://doi.org/10.1080/01621459.1998.10473788); [Lee, Huang and Zelen, 2004](https://doi.org/10.1191/0962280204sm377ra) | Progressive preclinical/clinical models, schedule sensitivity and risk-dependent examination schedules. Their population disease-entry models differ from the proportional crossing times below. | Primary abstracts |
| [Tsodikov et al., 1995, *Biometrics* 51:437–447](https://pubmed.ncbi.nlm.nih.gov/7662836/) | Explicit post-treatment recurrence surveillance, nonzero false-negative rate, and dynamic programming to minimize expected recurrence-detection delay. This defeats a first-post-treatment-optimization claim. | Primary abstract |
| [Hanin, Tsodikov and Yakovlev, 2001](https://doi.org/10.1016/S0895-7177(01)80023-6) | Stochastic latency, tumour growth and size-dependent detection, with a detection-size distribution objective and fixed-exam numerical schedules. Exponential growth itself is not a new scheduling model. | Publisher abstract; full text unavailable |
| [Ahern, Cheng and Shen, 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3230228/) | Risk-specific fixed-budget nonuniform screening and optimal exam counts within an equal-spacing class. | Full text |
| [Somda et al., online 2014 / print 2016](https://pubmed.ncbi.nlm.nih.gov/24567440/); [Chen and Yu, 2021](https://pubmed.ncbi.nlm.nih.gov/32712953/) | Multistate post-treatment follow-up optimization; the latter includes tailored recurrence schedules and competing death. Those general mechanisms are not claimed as new. | Primary abstracts |
| [Hijazy and Zempléni, 2020](https://doi.org/10.1002/qre.2707) | Random onset, gamma deterioration, threshold symptoms, degradation-dependent test sensitivity, and inspection optimization. | Full text |
| [Avanzini et al., 2020, *Science Advances*](https://pmc.ncbi.nlm.nih.gov/articles/PMC7732186/) | Mechanistic tumour evolution, ctDNA shedding, molecular sampling and sequencing errors; relapse detection size and lead-time predictions depend on testing frequency. The paper already compares monthly and quarterly surveillance. Its stochastic detection sizes are not deterministic assay limits. | Full text |
| [Abbosh et al., 2017, *Nature*](https://doi.org/10.1038/nature22364); [2023, *Nature*](https://doi.org/10.1038/s41586-023-05776-4) | TRACERx longitudinal MRD and molecular lead time. The 2023 median of 119 days (range 0–1137, 63 cases) assigns zero to absent or post-relapse detection; it is not a detected-only median. [Accepted manuscript](https://eprints.whiterose.ac.uk/id/eprint/203072/1/Tracking%20early%20lung%20cancer%20metastatic%20dissemination%20in%20TRACERx%20using%20ctDNA.pdf). | Full primary text/manuscript |
| [Reinert et al., 2019, *JAMA Oncology*](https://pmc.ncbi.nlm.nih.gov/articles/PMC6512280/) | Quarterly longitudinal sampling; 14 of 16 relapses detected in the longitudinal analysis. Detected molecular-to-radiological lead was **mean** 8.7 months, range 0.8–16.5. This supplies the aggregate prediction in Section 9. | Corrected full text |
| [Tie et al., 2022, NEJM DYNAMIC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9701133/) | Randomized ctDNA-guided adjuvant allocation using week-4/week-7 samples. This addresses treatment allocation, not optimal longitudinal surveillance dates. | Full text |
| [Henriksen et al., 2022, *Clinical Cancer Research*](https://pmc.ncbi.nlm.nih.gov/articles/PMC9401484/) | A selected 17-case repeated-positive subgroup had slow and fast fitted exponential ctDNA growth, with reported monthly multiplicative factors 1.25 and 2.43. Their log rates are approximately 0.223 and 0.888/month. This motivates heterogeneity but does not identify a population gamma/log-normal law. | Full text |
| [Henriksen et al., 2024, *Annals of Oncology*](https://doi.org/10.1016/j.annonc.2023.11.009) | Nationwide Danish cohort, including model-based evaluation of sampling frequency and potential lead time. Frequency effects are prior work. | Publisher article and indexed text |
| [Mohammadi et al., January 2025 preprint](https://www.qeios.com/read/5NA8P2/pdf), [medRxiv DOI](https://doi.org/10.1101/2025.01.07.25320131) | Direct ctDNA surveillance optimization for HPV-positive HNSCC: monthly hidden Markov recurrence, constant sensitivity/specificity, constrained POMDP and microsimulation. **Not the first ctDNA scheduling study.** Our question replaces fixed sensitivity with growth-linked thresholds and seeks analytic global guarantees. | Full 21-page preprint, explicitly not treated as a clinical guideline |
| [Mo et al., 2026, phase III FIND](https://pubmed.ncbi.nlm.nih.gov/42525894/) | Randomized quarterly ctDNA plus triggered imaging. Its reported 3.9-month difference is a between-arm difference in median clinical recurrence times, not the within-person molecular lead used by Reinert. Thus it cannot simply calibrate our lead variable. | Primary abstract/indexed publisher text; online 29 July 2026 |
| [Hyndman, 1996](https://robjhyndman.com/papers/hdrgraph.pdf) | Highest-density regions minimize volume for given probability; the paper credits Box and Tiao. The rearrangement step in Theorem 2 is classical. | Author full manuscript |
| [Nemhauser, Wolsey and Fisher, 1978](https://doi.org/10.1007/BF01588971) | Greedy approximation for monotone submodular cardinality optimization. Theorem 1 specializes this established guarantee. | Paper/primary metadata |

Searches combined `ctDNA optimal surveillance schedule`, `growth rate
heterogeneity surveillance`, `geometric inspection intervals`, `lognormal
gamma optimal screening`, `post-treatment surveillance false negatives`, and
`highest density interval coverage`, including recent 2025–2026 work. We found
no matching theorem for the precise model/bound below. This is a documented,
limited search, not an exhaustive systematic review. Related inaccessible
older full texts remain a priority risk. **Novelty is provisional and narrowly
mathematical; clinical meaningfulness is a testable motivation, not a result.**

## 2. Specified stochastic model and two different objectives

All logarithms are natural. Times use one fixed unit. A reference residual
burden is $B_*>0$; actual minimal residual burden is
$B_0=B_*e^Z$. A patient-specific rate $R>0$ and $Z$ are drawn jointly.
The regrowth process is

$$
B(t)=B_*e^{Z+Rt},\qquad t\ge0. \tag{1}
$$

Thus regrowth is stochastic across patients through rate and initial burden,
and deterministic conditional on those variables. This is the random-rate
exponential model allowed by the research question, **not** a birth–death
path approximation claimed without proof. Extinction, dormancy, treatment
changes, multiple lesions and time-varying growth are omitted.

**Shedding and sampling.** Assume a quasi-steady malignant DNA amount
$\kappa B$, fixed background DNA amount $h>0$, and true shed fraction

$$
f(B)=\frac{\kappa B}{h+\kappa B},\quad \kappa>0.
$$

Fix a true-fraction gate $0<\ell<1$. At each blood draw independently
conditional on $(R,Z)$, sample $K\sim\operatorname{Binomial}(m,f(B))$
with integer $m\ge1$. A positive requires both $f(B)\ge\ell$ and
$K\ge\lceil m\ell\rceil$. In particular, there are no false positives
below the true-fraction gate. This is an explicit idealized detection law,
not a claim that real assays implement an observable true-fraction gate.
It keeps a strict analytic LoD separate from false negatives due to sampling:

$$
p_m(B)=\mathbf1_{\{B\ge D\}}
 \sum_{k=\lceil m\ell\rceil}^{m}{m\choose k}
 f(B)^k[1-f(B)]^{m-k},\qquad
D=\frac{h\ell}{\kappa(1-\ell)}. \tag{2}
$$

There is no sequencing-error, clonal-hematopoiesis or false-positive decision
cost in this model. Those require a different objective. “Perfect threshold”
below means replacing (2) by $p_\infty(B)=\mathbf1_{B\ge D}$.
It is not a pointwise binomial limit exactly at the threshold. The fixed-burden
continuous-rate laws in Theorem 2, and the conditional-density assumption in
Theorem 4, give that boundary zero probability at fixed positive test times.
A continuous marginal rate law alone is insufficient under arbitrary
rate–burden dependence; Theorem 1 does not assume such a limit.

Clinical progression occurs on first reaching a fixed burden $C>D$.
Write

$$
a=\log(D/B_*),\quad b=\log(C/B_*),\quad 0<a<b,
\quad A=\frac{a-Z}{R},\quad P=\frac{b-Z}{R}. \tag{3}
$$

Assume $Z<a$ almost surely, so $A>0$ is first ideal detectability and
$P>A$ is progression. A schedule is $S=\{t_1<\cdots<t_n\}\subset(0,H]$,
with fixed integer budget $n$ and $H>0$. It specifies potential draws;
actual draws can stop on a positive or progression. No quantitative-result
adaptive scheduling, minimum visit gap, or mandatory terminal draw is imposed.
Let $D_S$ be the first positive strictly before $P$, or infinity.

An independent **scoring endpoint** $T\in(0,\infty]$ specifies which
progressions count. Define

$$
V_T(S)=\Pr(D_S<\infty,\ P\le T),\quad
p_T=\Pr(P\le T),\quad M_T(S)=p_T-V_T(S). \tag{4}
$$

Maximizing $V_T$ minimizes the absolute missed-progression probability;
if $p_T>0$, it also minimizes conditional miss $1-V_T/p_T$.
For $T=\infty$, all eventually relapsing paths count although testing
stops at $H$; failure to catch a later progression is still a miss for
that endpoint. For $T=H$, only progression within follow-up counts.
These are different estimands. A test at $H$ has no benefit for the
second endpoint, but can benefit the first.

The zero-imputed achieved biological lead and ideal-opportunity loss are

$$
G_T(S)=E[\mathbf1_{P\le T}(P-D_S)_+],
\quad \mathcal L_T(S)=E[\mathbf1_{P\le T}(P-A)]-G_T(S). \tag{5}
$$

Here $(P-\infty)_+=0$: no detection gives zero lead and loses the entire
opportunity $P-A$. If $T=\infty$, assume
$E(1/R)<\infty$ whenever (5) is used; (4) needs no such moment. For a finite
$T$, $0<P-A<P\le T$ on scored paths, so (5) is integrable.
Geometric optimality below concerns (4). Equal-calendar optimality in Theorem
5 concerns (5), under a different residual-burden law. Neither objective
dominates the other in general.

## 3. Theorem 1 — general noisy objective, inserted-test value, and finite-grid bound

**Assumptions.** Equations (1)–(4); any joint law of $(R,Z)$ satisfying
$R>0,Z<a$; conditional independence of test outcomes; fixed assay;
nonnegative weights and an available finite set of candidate draws for the
greedy statement. For the lead formula additionally use the integrability
condition after (5). Set
$q_t(R,Z)=\mathbf1_{t<P}p_m(B(t))$.

Then

$$
V_T(S)=E\left[\mathbf1_{P\le T}
 \left(1-\prod_{t\in S}(1-q_t)\right)\right], \tag{6}
$$

$$
G_T(S)=E\left[\mathbf1_{P\le T}\sum_{i=1}^n
(P-t_i)q_{t_i}\prod_{j<i}(1-q_{t_j})\right], \tag{7}
$$

and the exact value of inserting an available test $u\notin S$, while
retaining all existing dates, is

$$
V_T(S\cup\{u\})-V_T(S)
=E\left[\mathbf1_{P\le T}q_u\prod_{t\in S}(1-q_t)\right]. \tag{8}
$$

The set function $V_T$ is monotone and submodular. On a finite candidate
set with at least $n\ge1$ elements, sequential largest-marginal selection
gives

$$
V_T(S_{\rm greedy})\ge[1-(1-1/n)^n]V_T(S_{\rm grid}^*)
\ge(1-e^{-1})V_T(S_{\rm grid}^*). \tag{9}
$$

This is a grid guarantee, not a guarantee relative to all continuous dates.
For a finite mixture with rational weights and rational $q_t$, (6), (8)
and the grid optimum can be calculated exactly with `Fraction`.

**Proof.** Conditional on the latent path, independence gives the probability
of all failures as the product in (6). First detection at $t_i$ has
probability $q_{t_i}\prod_{j<i}(1-q_{t_j})$, proving (7). Subtracting (6)
with and without $u$ gives (8). Every factor is in $[0,1]$; adding dates
cannot decrease value and can only decrease the marginal of another date.
Taking expectations preserves these inequalities.

For completeness, if $S_j$ is the greedy set after $j$ choices and
$O$ an optimum of size $n$, telescoping and diminishing marginals give
$V_T(O)-V_T(S_j)\le\sum_{u\in O\setminus S_j}\Delta(u\mid S_j)
\le n\max_u\Delta(u\mid S_j)$.
Thus the greedy optimality gap contracts by at most $1-1/n$ at each step.
Starting from value zero and iterating $n$ times proves (9). This is the
classical argument cited in the literature ledger, not a new approximation
theorem. Empty schedules, zero sensitivities and unscored paths satisfy all
formulas directly. **QED.**

## 4. Theorem 2 — exact global geometric schedule and the horizon correction

**Assumptions.** Perfect threshold detection; $Z=0$; fixed $0<a<b$,
$H>0$, $T\in(0,\infty]$, and $n\ge1$. The variable
$X=\log R$ has a continuous strictly positive density $g$ on the real
line, tending to zero at both tails, strictly increasing up to a unique mode
and strictly decreasing afterwards. Define

$$
q=b/a>1,\quad w=\log q,\quad L=nw,\quad
c=\log\max(a/H,b/T), \tag{10}
$$

with $b/\infty=0$. Let $u_L$ be the unique solution of
$g(u_L)=g(u_L+L)$, and let $x_n=\max(c,u_L)$, $r_n=e^{x_n}$.
The unique optimal ordered schedule, up to zero-probability endpoint
conventions, is

$$
\boxed{t_i^*=\frac{a}{r_n q^{n-i}},\quad i=1,\ldots,n},
\qquad \boxed{\frac{t_{i+1}^*}{t_i^*}=q}. \tag{11}
$$

Its absolute detection probability is

$$
\boxed{V_n^*=\Pr(r_n\le R<r_nq^n)}. \tag{12}
$$

For log-normal growth $X\sim N(\mu,\sigma^2)$, $\sigma>0$,

$$
r_n=\max\left(a/H,b/T,e^{\mu-nw/2}\right),\qquad
V_n^*=\Phi\left(\frac{\log r_n+nw-\mu}{\sigma}\right)
-\Phi\left(\frac{\log r_n-\mu}{\sigma}\right). \tag{13}
$$

For $R\sim\operatorname{Gamma}(k,\theta)$, shape $k>0$, scale
$\theta>0$,

$$
\boxed{r_n=\max\left(a/H,b/T,\frac{k\theta\,nw}{q^n-1}\right)},\quad
V_n^*=F_{k,\theta}(r_nq^n)-F_{k,\theta}(r_n). \tag{14}
$$

The gamma formula holds even for $k\le1$: it is the density of **log
rate**, not rate, which must be unimodal. For these positive continuous
families every equally calendar-spaced schedule of $n\ge3$ distinct dates
is strictly suboptimal for (4). This does not say that the first waiting time
and every later gap form one geometric sequence: the inter-test gaps
$t_{i+1}-t_i=(q-1)t_i$ form a geometric sequence, while the initial wait is
$t_1$. For $q<2$, the second waiting interval is shorter than the first.

**Proof.** A perfect test at $t$ detects before progression exactly when

$$
a/R\le t<b/R
\quad\Longleftrightarrow\quad
\log a-\log t\le X<\log b-\log t. \tag{15}
$$

Each test therefore covers an interval of length $w$ in $X$. Since
$t\le H$, its left endpoint is at least $\log(a/H)$; only
$X\ge\log(b/T)$ scores. The scored union is a measurable subset
$E\subset[c,\infty)$ of length at most $L=nw$.

Let the mode of $g$ be $m_g$. On
$[m_g-L,m_g]$, $g(u)-g(u+L)$ is continuous, strictly increasing,
negative at the left endpoint and positive at the right. It therefore has
one zero $u_L$. Outside that range the monotonicity of $g$ precludes
another zero. The interval $I=[x_n,x_n+L]$, with $x_n=\max(c,u_L)$,
contains the highest density values available in $[c,\infty)$.
More explicitly, set $\lambda=g(x_n+L)>0$. On the interior of $I$,
$g>\lambda$; outside $I$ within $[c,\infty)$, $g<\lambda$
almost everywhere. This holds both when endpoints have equal density and
when the lower boundary $c$ binds. Thus

$$
\int_Eg-\int_Ig
=\int_{E\setminus I}(g-\lambda)
-\int_{I\setminus E}(g-\lambda)
+\lambda(|E|-L)\le0. \tag{16}
$$

Equality requires $E=I$ almost everywhere. The dates (11) attain it:
their rate intervals tile $[r_n,r_nq^n)$ without overlaps or gaps, and
$t_n=a/r_n\le H$. Conversely, equality for $n$ intervals each of
length $w$ permits no positive-length overlap or scored-support wastage;
their union must be $I$, so they must tile it. This forces (11).

For a normal log rate, equal-density endpoints are symmetric about $\mu$,
giving (13). For gamma rates the log-rate density is
$g(x)=e^{kx-e^x/\theta}/[\Gamma(k)\theta^k]$, strictly increasing until
$e^x=k\theta$ and decreasing afterwards. Its endpoint equation is
$kL=e^x(e^L-1)/\theta$, proving (14). Three or more equally spaced
positive dates cannot simultaneously have a constant ratio greater than
one; uniqueness therefore proves strict suboptimality. For one or two dates
no such general calendar-spacing comparison is asserted. **QED.**

**Why the endpoint matters.** For eventual relapses $T=\infty$, the floor
is $a/H$. For progression by follow-up $T=H$, it is $b/H$, and the
latest optimal test is at most $H/q$. Substituting $a/H$ for $b/H$
in the second problem is an error. Conditional capture divides (12) by
$\Pr(R\ge b/T)$, not by the number already detected. For infinitely many
tests the ideal attainable value approaches
$\Pr(R\ge\max(a/H,b/T))$; slow rates below $a/H$ never become
detectable while tests remain available.

## 5. Theorem 3 — heterogeneity and the exact value of another test

**Assumptions.** The normal or gamma family, parameters and objective of
Theorem 2; parameters named “fixed” below, including $a,b,H,T$, stay fixed.
The first and second parts apply at all horizons; the explicitly
nonbinding-horizon claims assume the maximum in (13)/(14) selects the
unconstrained expression for the budgets and parameters compared.

1. **Same median log-normal populations.** With $\mu$ fixed, all optimal
   dates are independent of $\sigma$. If the horizon is nonbinding,
   $V_n^*=2\Phi(nw/(2\sigma))-1$, which decreases strictly with
   $\sigma$.
2. **Same arithmetic mean populations.** With mean rate $\bar r$ fixed,
   a log-normal law has $\mu=\log\bar r-\sigma^2/2$. Every optimal date
   is nondecreasing in $\sigma$, strictly increasing and scaling by $e^{\sigma^2/2}$ before
   the horizon binds, and remaining constant after binding. In contrast,
   for gamma laws with $k\theta=\bar r$, the entire optimal schedule is
   independent of shape $k$, hence of coefficient of variation
   $1/\sqrt{k}$, at every horizon. With a nonbinding horizon, gamma
   optimal capture increases strictly with $k$.
3. **Reoptimized additional test.** Put $V_0^*=0$. The exact value is
   $\Delta_n=V_{n+1}^*-V_n^*$, with both terms given by (13) or (14).
   It is positive and strictly decreases with integer $n\ge0$. This
   reoptimizes all dates; it is distinct from insertion value (8).
4. **An explicit peak in the value of another test.** For fixed $\mu$,
   $n\ge1$, and a horizon nonbinding for both budgets,

   $$
   \Delta_n(\sigma)=2\left[\Phi\left(\frac{(n+1)w}{2\sigma}\right)
   -\Phi\left(\frac{nw}{2\sigma}\right)\right]. \tag{17}
   $$

   As a function of heterogeneity, its unique maximum occurs at

   $$
   \boxed{\sigma_*=
    w\sqrt{\frac{2n+1}{8\log(1+1/n)}}}. \tag{18}
   $$

   Thus an additional test has greatest value at an intermediate log-rate
   spread, even though total optimal capture decreases with spread in this
   nonbinding, fixed-median comparison.

**Proof.** In (13) the maximizing lower log-rate endpoint does not depend on
$\sigma$ if $\mu$ is fixed. If the mean is fixed it is
$\max(c,\log\bar r-\sigma^2/2-nw/2)$, a nonincreasing function of
$\sigma$, giving the asserted increasing dates and the plateau. The gamma
lower rate is $\max(e^c,\bar r\,nw/(e^{nw}-1))$, independent of $k$.

To prove the gamma capture comparison, set $Y=R/\bar r$,
$L=nw$, $A_L=L/(e^L-1)<1$, $B_L=e^LA_L>1$.
Define $I(y)=y-1-\log y$. Because $B_L-A_L=\log(B_L/A_L)=L$,
$I(A_L)=I(B_L)$, and $[A_L,B_L]$ is exactly a sublevel set of
$I$. The density of $Y$ is
$C_k y^{-1}e^{-kI(y)}$. Increasing $k$ multiplies it, up to a
normalizing constant, by a strictly decreasing function of $I(y)$.
That multiplier is larger inside this sublevel set than outside, with
strict inequality on sets of positive probability. Comparing its conditional
means inside and outside proves that the normalized probability of the
sublevel set strictly increases. This proves the gamma claim without
equating different distribution families.

For diminishing values, allow a continuous total log-length $L\ge0$
and call the maximizing mass $v(L)$. For normal densities the unconstrained
left endpoint $\mu-L/2$ decreases and the right endpoint increases.
For gamma densities the rate endpoints are
$k\theta L/(e^L-1)$ and $k\theta L/(1-e^{-L})$, respectively
strictly decreasing and increasing (differentiate, using $e^L>1+L$).
Their logarithms bracket the log-rate density mode $\log(k\theta)$.
On binding the lower constraint, the
left endpoint stays at $c$ and the right is $c+L$. The optimal intervals
are nested. Before binding, differentiation of their integral, using equal
endpoint densities, gives $v'(L)=g(u_L+L)$. After binding it gives
$v'(L)=g(c+L)$. The derivatives agree at the transition. These positive
boundary densities strictly decrease with $L$: the right endpoint is above
the log-rate mode, including at and after a binding transition. Thus $v$ is strictly
increasing and strictly concave, including the continuous extension
$v(0)=0$. Integrating $v'$ over consecutive length-$w$ intervals
proves positive, strictly diminishing $\Delta_n$, including $n=0$.

Finally differentiate (17), writing $\phi$ for the standard normal
density:

$$
\Delta_n'(\sigma)=\frac{w}{\sigma^2}
 \left[n\phi\left(\frac{nw}{2\sigma}\right)
 -(n+1)\phi\left(\frac{(n+1)w}{2\sigma}\right)\right].
$$

The bracket is zero exactly when
$(2n+1)w^2/(8\sigma^2)=\log(1+1/n)$; it is positive below that
$\sigma$ and negative above. Both limits of (17) as $\sigma$ tends to
zero or infinity are zero. This proves (18) and uniqueness. **QED.**

**Boundary cautions.** Variance alone does not specify the population law or
what location parameter remains fixed. Nor does greater heterogeneity always
reduce finite-window capture: if the log-normal horizon floor lies above
$\mu$, capture tends to zero as $\sigma\downarrow0$, whereas positive
spread places some mass in the accessible tail. There is no general theorem
here that more heterogeneity requires front-loading.

## 6. Theorem 4 — random residual burden and finite sampling: uniform regret bound

**Assumptions.** Equations (1)–(4), with $0<a<b$, $n\ge1$, and
$|Z|\le\eta<a$ almost surely. Conditional on every $Z$ outside a null
set, $X=\log R$ has a density bounded above by $M<\infty$.
Independence of $X$ and $Z$ is sufficient but not necessary. The
**marginal** density of $X$ satisfies Theorem 2, so its $Z=0$
reference has optimizer $S_0$. For example an independent log-normal has
$M=1/(\sigma\sqrt{2\pi})$; an independent gamma has
$M=k^ke^{-k}/\Gamma(k)$.

Choose $0<\delta<f(C)-\ell$, put

$$
d_\delta=\log\frac{f^{-1}(\ell+\delta)}{D}
=\log\frac{(\ell+\delta)(1-\ell)}{\ell(1-\ell-\delta)},
\quad h_a=\log\frac{a}{a-\eta},\quad h_b=\log\frac{b}{b-\eta},
$$

$$
\epsilon_B=M\{n(h_a+h_b)+\mathbf1_{T<\infty}h_b\},\qquad
\epsilon_N=nM\log\left(1+\frac{d_\delta}{a-\eta}\right)
+e^{-2m\delta^2}. \tag{19}
$$

Let $V_0$ denote reference perfect-threshold value with $Z=0$,
$V_B$ perfect-threshold value with the actual random $Z$, and
$V_N$ the actual noisy value. Uniformly over every feasible $n$-test
schedule,

$$
|V_B(S)-V_0(S)|\le\epsilon_B,\qquad
0\le V_B(S)-V_N(S)\le\epsilon_N. \tag{20}
$$

Consequently, even if an actual noisy optimum is only a supremum,

$$
\boxed{\sup_{|S|=n}V_N(S)-V_N(S_0)
\le\min(1,2\epsilon_B+\epsilon_N)}. \tag{21}
$$

For conditional missed-progression regret divide the uncapped right side by
the **actual** $\Pr(P\le T)>0$, then cap at one. The reference and actual
progression probabilities need not be equal. There is no such denominator
for the absolute-probability bound (21).

**Proof: residual burden.** For fixed $Z=z$, a test covers the log-rate
interval with endpoints $\log(a-z)-\log t$ and
$\log(b-z)-\log t$. Relative to $z=0$, their absolute movements are
at most $h_a$ and $h_b$; the negative perturbation has smaller magnitude
since $\log(1+\eta/a)\le-\log(1-\eta/a)$, and similarly for $b$.
The symmetric difference of two intervals has length at most the sum of
their endpoint movements. For $n$ intervals the union symmetric difference
therefore has length at most $n(h_a+h_b)$.

If $T<\infty$, the scoring half-line also changes from
$X\ge\log b-\log T$ to $X\ge\log(b-z)-\log T$, a movement at most
$h_b$. The symmetric difference of intersections is contained in the
union of these two symmetric differences. Conditional probability is at
most $M$ times their lengths. Integrating over $Z$, while keeping the
same jointly sampled $R$ in both models, proves the first bound in (20).
This argument does not require independence of residual burden and rate.

**Proof: sampling.** The noisy gate can never detect a path that the actual
perfect threshold fails to detect. For a fixed test, the strip
$D\le B(t)<f^{-1}(\ell+\delta)$ corresponds, conditional on $Z=z$,
to a log-rate interval of length
$\log(1+d_\delta/(a-z))\le\log(1+d_\delta/(a-\eta))$.
The probability of at least one strip test is at most the first term of
$\epsilon_N$.

Outside that event, a scored perfect detection implies a scheduled test
before progression with true fraction at least $\ell+\delta$. Choose
the earliest such test as a function of the latent path. Its binomial
failure probability is at most $e^{-2m\delta^2}$ by the bounded-summand
Hoeffding inequality. Failure of every test implies failure of this one.
There is therefore **one**, not $n$, copy of the exponential term.
Ceiling the count threshold causes no extra term: call failure is exactly
$K/m<\ell$. Integrate over latent paths and add the strip bound.

For clarity, the concentration step follows from the exponential Markov
inequality: for independent Bernoulli variables, the log moment-generating
function of each centered summand is at most $s^2/8$. One proof bounds
the second derivative of the tilted Bernoulli log moment-generating function
by $1/4$, with value and derivative zero at the origin, and integrates
twice. Thus the probability that the sample mean lies $\delta$ below
its expectation is at most
$\inf_{s>0}e^{-sm\delta+ms^2/8}=e^{-2m\delta^2}$.

Finally $V_N(S)\le V_B(S)\le V_0(S)+\epsilon_B\le
V_0(S_0)+\epsilon_B$, whereas
$V_N(S_0)\ge V_0(S_0)-\epsilon_B-\epsilon_N$.
Subtract and take the supremum to obtain (21). **QED.**

**Proved scaling regime.** For fixed $a,b,\ell,M,n,H,T$, let
$\eta\to0$, $m\to\infty$, and choose
$\delta_m=\sqrt{\log m/(2m)}$. Eventually it is an admissible margin.
Since $h_a,h_b=O(\eta)$, $d_{\delta_m}=O(\delta_m)$, and
$e^{-2m\delta_m^2}=1/m$, (21) tends to zero at
$O(n\eta+n\sqrt{\log m/m}+1/m)$, with all constants displayed in
(19). No exact finite-noise geometric optimality is inferred. For large
residual uncertainty or sparse molecules the bound can be vacuous; that is
reported rather than hidden.

**A rigorous failure of exact noisy geometry.** Take $a=1,b=2,H=2$,
$Z=0$, $X\sim N(0,10^{-6})$, $m=1$, $\ell=1/2$, and
$T=\infty$. The shed fraction is
$f(B(t))=\operatorname{logistic}(Rt-1)$. Any two dates in geometric
ratio two have disjoint ideal log-rate windows, so at most one can detect
any given path. Before progression its positive probability is below
$\operatorname{logistic}(1)=e/(1+e)<11/15$ (using $e<11/4$).
In contrast, tests at $3/2$ and $8/5$ are both eligible whenever
$|X|\le1/100$. That event has probability at least $99/100$ by
Chebyshev's inequality. Each conditional positive probability is at least
$1/2$, so independent sampling gives total capture at least
$(99/100)(3/4)=297/400>11/15$.
This refutes exact finite-noise geometric optimality even with a strictly
unimodal log-normal rate density. The loose bound (21) correctly does not
promise closeness at $m=1$.

## 7. Theorem 5 — an exact equal-calendar benchmark for lead-time loss

**Assumptions.** Perfect threshold; fixed rate $R=r>0$; random residual
burden $B_0=D e^{-rA}$ with $A\sim\operatorname{Uniform}[0,H]$;
$C=D e^{rW}$ with fixed $W\ge H>0$; scoring endpoint
$T=\infty$; $n\ge1$ dates in $(0,H]$. Thus progression is
$P=A+W$, and a missed path loses the full $W$. The probability-zero
initial endpoint $A=0$ has no effect.

For $t_0=0$,

$$
\mathcal L_\infty(S)
=\frac{1}{2H}\sum_{i=1}^n(t_i-t_{i-1})^2
+\frac{H-t_n}{H}W. \tag{22}
$$

The unique optimum is $t_i=iH/n$, and

$$
\boxed{\mathcal L_n^*=\frac{H}{2n}},\qquad
\boxed{\mathcal L_n^*-\mathcal L_{n+1}^*
=\frac{H}{2n(n+1)}}. \tag{23}
$$

For no tests, loss is $W$; the first test's value is $W-H/2$.
The equally spaced result is a standard integrated waiting-time/convexity
benchmark, not a claimed discovery.

**Proof.** For onset $A\in(t_{i-1},t_i]$, detection occurs at $t_i$
and loses $t_i-A$; the next test is before progression almost surely
because $W\ge H$. Integrating over the uniform onset gives the $i$-th
square in (22). Onsets after the last test are missed and contribute the
remaining term. For fixed $x=t_n$, Cauchy–Schwarz gives
$\sum_i(t_i-t_{i-1})^2\ge x^2/n$, with equality exactly at equal gaps.
The lower bound $x^2/(2nH)+(H-x)W/H$ has derivative
$x/(nH)-W/H<0$ for $x<H$, since $W\ge H$.
Hence $x=H$ uniquely minimizes it, and equality forces all gaps to be
$H/n$. Subtraction proves (23). The boundary $W=H,n=1,x=H$
has zero derivative only at the endpoint and does not destroy uniqueness.
**QED.**

This does not contradict Theorem 2: here onset heterogeneity comes from
residual burden, the rate is fixed, and the objective is lost lead rather
than capture probability. Clipping the loss at $H$ would create a
different optimization problem; no uncounted “free final test” is used.

## 8. Theorem 6 — a falsifiable aggregate lead-time prediction by alternating thinning

**Assumptions.** A fixed original detected cohort has nonnegative leads
$L$, measured from first positive blood draw to an unchanged clinical
recurrence date. Draws occur exactly every $d>0$ time units and positivity
persists after first positive. Retain one of the two alternating draw parities
with probability $1/2$, independently of a case's lead. A retained next
draw is available if it precedes recurrence; otherwise the thinned lead is
zero. Keep the original cohort denominator, including thinned misses. No
newly triggered imaging, changed assay, treatment intervention or follow-up
truncation is introduced.

Then, averaging over the retained parity,

$$
E[L_{\rm thin}\mid L]
=\frac{L+(L-d)_+}{2}
=L-\frac12\min(d,L). \tag{24}
$$

If only mean $\bar L$ and support $l\le L\le u$, $0\le l<u$,
are available, set

$$
Q=\min(d,l)+\frac{\bar L-l}{u-l}
 [\min(d,u)-\min(d,l)].
$$

The sharp mean/support bounds are

$$
\boxed{\bar L-\tfrac12\min(d,\bar L)
\le E L_{\rm thin}\le\bar L-\tfrac12Q}. \tag{25}
$$

For $l=u$, (24) gives the single value. Bounds are sharp over all laws
with that mean and **contained** support; observed endpoints need not have
positive mass in that optimization class.

**Proof.** If the first positive draw is retained, the lead is unchanged. If
it is removed, the next positive draw is $d$ later unless recurrence has
already occurred; zero-imputation gives $(L-d)_+$. This proves (24).
The function $x\mapsto\min(d,x)$ is concave, so it lies above its
chord on $[l,u]$ and Jensen gives
$Q\le E\min(d,L)\le\min(d,\bar L)$. Substituting into (24), with
the inequality directions reversed on subtraction, proves (25).
The endpoint two-point distribution attains the lower expectation of the
minimum, hence the upper thinned-mean bound. A point mass at $\bar L$
attains the other bound. **QED.**

## 9. Quantitative cohort prediction and attempted falsification

Reinert's corrected publication reports a mean of **8.7 months**, range
**0.8–16.5 months**, quarterly blood collection, and persistent positivity
between first post-treatment detection and radiological recurrence. Treating
these rounded values as nominal inputs to (25), with $d=3$, gives the
exact rational statement

$$
\boxed{\frac{36}{5}=7.20\ \le E L_{\rm thin}\le
\frac{6081}{785}=7.746496815\ldots\ \mathrm{months}}. \tag{26}
$$

Equivalently, the predicted mean loss from alternating thinning is between
0.9535031847… and 1.50 months. This is a prediction for the parity-averaged
**same detected cohort**, with new misses assigned zero, not a recommendation
for six-month monitoring and not a prediction of survival or treatment
benefit. It does not recondition on remaining detectable cases.

**Concrete falsification procedure.** In a permitted reanalysis, preserve
clinical recurrence dates, apply each of the two retained parities, compute
zero-imputed lead for every member of the original detected cohort, and
average the two cohort means. Under Theorem 6's exact-cadence and persistence
assumptions, a result outside (26) contradicts the specified observation
model and nominal aggregates. Real missing visits, irregular intervals or
loss of persistent positivity can invalidate its assumptions; record them
before evaluating the prediction. No such reanalysis was performed here,
and no patient records were obtained.

**Rounding sensitivity.** If the three reported one-decimal summaries were
rounded to nearest tenth, use mean in $[8.65,8.75]$ and a conservative
support $[0.75,16.55]$. Then the distribution-free envelope widens to
$[7.15,7.805379747\ldots]$ months. For the upper end insert mean 8.75 into
the chord bound on that wider support; its slope in the mean is positive.
This is an assumption about rounding, not a source-provided uncertainty
interval. Sampling variability and irregular cadence are separate.

**A stronger null fails, but does not identify growth heterogeneity.** A
common biological detectable window $W=\log(C/D)/R$, perfect periodic
testing of spacing $d$ already underway at ideal detectability with no
missing or delayed initial visit, and continuous observation of clinical progression
imply detected leads in an interval of width at most $d$. The nominal
Reinert range has width 15.7 months, exceeding 3. That contradicts the
**conjunction** of those assumptions. It does not isolate its cause:
radiological ascertainment is intermittent; assay misses, threshold
heterogeneity and growth heterogeneity also change the range. If clinical
ascertainment has lag in $[0,I]$, the corresponding width bound becomes
$d+I$; an unknown lag makes the unadjusted rejection biologically
uninformative. Aggregate lead means/medians therefore do not identify
$a,b$, $B_0$, shedding, or a population rate distribution.

The TRACERx 2023 zero-imputed median, Reinert's detected-cohort mean, and
FIND's between-arm difference of clinical recurrence medians are different
estimands. They are recorded in the literature ledger and are deliberately
not pooled. The full model's optimum remains empirically uncalibrated.

## 10. Numerical evidence and reproducibility

Only the three files in this directory and
`tests/test_research_ctdna_surveillance_schedule.py` are added. No workflow,
dependency, existing research note, index, catalog, results directory or
existing test is changed. The repository's existing PR validation workflow
already discovers `tests/test_*.py`.

From the repository root, with Python 3.12 or later:

```bash
python research/ctdna-surveillance-schedule/verifier.py
python research/ctdna-surveillance-schedule/simulation.py
python -m unittest discover -s tests -p 'test_research_ctdna_surveillance_schedule.py' -v
```

All inputs are synthetic mathematical parameters or the explicitly attributed
published aggregate numbers. The programs make no network calls, require no
third-party packages, and print to stdout without producing result files.

**Recorded local execution after referee corrections (Python 3.12):** all
22 dedicated tests passed in 2.073 seconds, below the 60-second requirement;
all 277 repository Python tests passed in 22.404 seconds; schema/reference
validation passed; and all 26 repository JavaScript tests passed. The verifier
and simulation command-line examples also completed successfully. Timings are
environment-dependent. GitHub Actions results are reported on the PR rather
than inferred from these local checks.

### Separation of checks

| Component | What it checks | What it does not establish |
| --- | --- | --- |
| `verifier.py` / `Fraction` | Exact finite-mixture success, binomial tails, insertion values, grid exhaustive/greedy optimization, uniform-onset loss and thinning bounds | Proof for all continuous parameters merely by enumerating examples |
| `verifier.py` / `Decimal` | 70-digit working precision for normal/gamma CDFs, optimum dates/values and explicit perturbation bounds | Outward-rounded interval certification or guaranteed relative tail accuracy |
| `simulation.py` | Direct latent growth trajectories, first-positive/progression logic, individually sampled Bernoulli molecules, uniform onsets and parity thinning | Clinical data, proof, or independence of biological validation |
| New unittest | Independently computed finite-outcome/binomial distributions and interval integration; reference CDFs, parameter families, boundary cases, fixed-seed Monte Carlo and counterexamples | External peer review or novelty certification |

The simulator does **not** import the verifier or evaluate the optimum
formulas. Separate latent-patient and assay random streams preserve the
patient draws when assay settings change. Shared dates use the same assay
outcome across compared schedules; distinct dates have conditionally
independent molecules. Bounded-statistic checks use Hoeffding tolerances
with individual error probability $10^{-10}$ and explicit range widths;
fewer than 100 such comparisons give a union bound of $10^{-8}$ under
independent sampling. The three unbounded log-normal loss checks instead
use Chebyshev's inequality with individual error probability $10^{-3}$ and
the analytic second-moment bound
$(b-a)^2 E[R^{-2}]=(b-a)^2e^{-2\mu+2\sigma^2}$, giving a separate union
bound of $0.003$. No tolerance uses estimated empirical standard errors.
Fixed seeds make these checks
reproducible; stochastic agreement is evidence of implementation consistency,
not proof of the theorems.

The Decimal normal CDF integrates the exponential series with guard digits;
gamma CDFs use a lower-tail series / upper-tail continued fraction with
Stirling-shift normalization. The latter uses exact rational Bernoulli
coefficients. Tests compare against separate `math.erf`, elementary
integer/half-integer gamma identities and independent density quadrature
for general positive shapes, including shapes below one.
For extreme tails cancellation can yield zero at the working precision;
`conditional_value` is then `None` if its numerical denominator is zero.
This is not a claim that the true tail probability is zero. The proofs do
not depend on these transcendental routines.

### Illustrative synthetic values, not fitted cohorts

For $a=1,b=2,n=3,H=12,T=\infty$:

| Rate law | Optimal dates | Optimal absolute capture | Reoptimized fourth-test gain |
| --- | --- | --- | --- |
| $\log R\sim N(0,1)$ | 0.70710678, 1.41421356, 2.82842712 | 0.7015303523 | 0.1328126097 |
| $R\sim\Gamma(2,0.5)$ | 0.84157211, 1.68314421, 3.36628843 | 0.8304052576 | 0.0973193880 |

Those two rows do not hold the same mean rate; they are reproducibility
fixtures, not a variance comparison. The same-mean and same-median
comparisons are separately checked according to Theorem 3.
For $a=1,b=2,H=T=1,n=3,\log R\sim N(0,1)$, the correct optimal dates
are $(1/8,1/4,1/2)$, giving absolute capture 0.2413279779 and conditional
capture 0.9886090948 among progressions by time one. A terminal test at one
would contribute nothing to that endpoint.

## 11. Hostile-referee record and falsification ledger

Proofs and implementations were derived and checked in separate investigator
and critic tasks, then reconciled. This is AI-assisted self-review, not
external scientific review. The following attempted overclaims were rejected
or corrected before submission.

| Attempt or edge case | Outcome and resulting scope |
| --- | --- |
| “First optimized ctDNA surveillance schedule” | Rejected by the 2025 preprint and extensive older surveillance theory. Novelty claim narrowed in Section 1. |
| Use $a/H$ for progressions by $H$ | Incorrect: the scoring event adds floor $b/H$. Theorem 2 and tests handle both endpoints. |
| Assume gamma rate density is unimodal only for shape greater than one | Wrong variable: log-rate density is strictly unimodal for every positive shape. |
| Call geometric dates “all waiting intervals geometric” | Initial wait differs; inter-test gaps form the geometric sequence. Explicitly corrected. |
| Claim more variance always front-loads tests | Refuted by Theorem 3: gamma same-mean dates stay fixed, log-normal same-mean dates move later until constrained. |
| Claim more variance always decreases finite-horizon capture | False when an accessible tail gains mass; monotone capture claims restricted to nonbinding horizons. |
| Extend geometric optimality to arbitrary rate mixtures | False. With $a=1,b=2$, equiprobable rates 1 and 100, dates 1 and 0.01 detect both; any ratio-two pair can detect at most one. Unimodality matters. |
| Extend exact geometry to finite sampling | False even for narrow log-normal rates; Section 6 gives a parameter-explicit proof and tests simulate it. Overlapping windows can provide a second chance. |
| Average perfect optima over random residual burdens | Invalid because the window width and location vary jointly. Theorem 4 supplies a uniform bound instead. |
| Ignore the changing clinical-scoring event under random burden | Corrected by the extra $Mh_b$ term in (19) for finite $T$. |
| Infer zero threshold-boundary probability from a continuous marginal rate | False under arbitrary dependence: $Z=a-Rt_0$ places every path at the threshold at $t_0$. The continuity statement now explicitly uses fixed burden or conditional densities. |
| Multiply the Hoeffding term by every planned test | Valid but unnecessarily loose. Selecting the earliest high-margin test gives the sharper single exponential term. |
| Treat reoptimized extra-test value as insertion value | Different operations, separately expressed in (8) and (17). |
| Count a free final test in the lead-loss benchmark | Rejected. Every date in (22) consumes budget; a missed path loses $W$. |
| Treat Reinert's 8.7 as a median, or TRACERx's 119 days as detected-only | Both incorrect; source endpoint definitions retained. |
| Fit growth heterogeneity from one aggregate lead statistic | Unidentifiable without further assumptions; no fit or clinical calibration claimed. |
| Call fixed-seed simulation proof or successful CI clinical validation | Rejected; numerical and scientific status are distinguished throughout. |
| Trust one numerical integration reference without inspecting failures | Two initial gamma-reference checks failed because a transformed quadrature was inaccurate for shapes above one. Integrating the ordinary density for those shapes, retaining the singularity-removing substitution below one, resolved the discrepancy. The production evaluator was unchanged. |
| Describe empirical-standard-error tolerances as concentration bounds | The mismatch was caught in review. The final tests use explicit Hoeffding ranges and analytic second-moment Chebyshev bounds as documented above. |

**Boundary audit.** No tests give zero capture; $n=1$ has no spacing ratio;
$n=2$ is not covered by the strict arithmetic-spacing comparison; exact
progression-time draws cannot count; touching rate intervals have
zero-probability boundaries under the continuous theorem; horizon equality is
allowed; all $k>0$ are allowed for capture but uncensored expected lead
requires $k>1$ for a gamma rate. The perturbation bound requires
$\eta<a$, a finite conditional density bound, and a margin below clinical
burden. The code rejects invalid inputs rather than silently interpreting
them as another model.

## 12. Novelty assessment, limitations and open questions

| Result | Honest status |
| --- | --- |
| Theorem 1 | Established conditional-independence/submodular machinery, specialized and reproduced for this model. |
| Theorem 2 | Candidate model-specific analytic solution: threshold windows become log-rate coverage with exact horizon-aware placement. Highest-density extremality is prior mathematics. Priority of this application is unestablished. |
| Theorem 3 | Proved consequences and refinements of the explicit solution, including a same-mean family contrast and an interior maximum in added-test value. No claim of first use of exponential tilting or concavity. |
| Theorem 4 | Candidate uniform transfer bound for this particular gated assay and random residual-burden model; fully explicit conditions, no assertion of exact finite-noise scheduling. Priority unestablished. |
| Theorem 5 | Standard equal-spacing benchmark under a specified onset law and a different objective. |
| Theorem 6 | Elementary exact thinning identity plus sharp moment bounds; its application yields an untested, falsifiable published-aggregate prediction. |

The clinically motivated question is whether allocation should follow calendar
time, detectability risk, or growth heterogeneity under a fixed assay budget.
The results show why this cannot be answered by variance, lead time, or assay
LoD alone. They also expose how biological progression and scheduled imaging
are different endpoints. That is potential relevance to experimental and
study design, **not demonstrated clinical usefulness**.

Material limitations include unknown post-treatment residual burden; unknown
correlation of growth, shedding and residual disease; the hard true-fraction
gate; fixed background cfDNA and molecule input; no non-shedding or cured
population; fixed clinical threshold; no competing death or treatment
effects; no harm/cost from false positives; no tumour-location effects;
unrestricted calendar dates; and uncalibrated gamma/log-normal families.
If a cured fraction is independent and all tests on it are negative, it
merely rescales the detection objective; otherwise a new model is required.
Equal burden thresholds across patients are especially consequential for
the ratio $q$. The broad random-burden/noisy problem is not declared
exactly solved by the geometric limiting case.

**Open, unproved questions.** Exact continuous-date optima for finite noisy
assays with overlap, wide or unbounded residual uncertainty, multimodal rate
populations, coupled branching/shedding paths, assay-specific false-positive
costs, and quantitative-result adaptive policies remain outside these
theorems. No proposed rule for these cases is promoted to a theorem or to
patient advice. The present note makes no assertion that a new clinical
surveillance policy is ready for use.
