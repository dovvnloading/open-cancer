# What mutation-rate scaling can resolve about Peto's paradox

**Nonclinical, unreviewed theoretical and computational research. No medical
advice. No biological, clinical, or independent peer-review validation is
claimed.**

## 1. Verdict and contribution

**Biological verdict: not identifiable from these published aggregates.**
Inverse-lifespan mutation scaling does not, on its own, cancel the target-cell
factor in the specified multistage models. Under an explicit proxy mapping of
intestinal passenger mutation rates, adult body mass, adult remaining life
expectancy and zoo cancer mortality, a common-rate-calibration/common-step
model requires substantial additional modulation. That conditional deficit
does not identify extra biological steps, an immune effect, repair efficiency,
or the number of susceptible targets. Missing uncertainty and the observation
mapping prevent a stronger biological verdict.

The analysis deliberately keeps these statements separate:

* **A mathematical result:** equal lifetime functional mutation clocks leave
  cumulative hazard proportional to independent target-cell number, exactly,
  including probability saturation. The required allometric mutation exponent
  depends on the number of steps and on the division-clock convention.
* **A conditional aggregate result:** the cell-aware proxy family predicts
  worse than the intercept and burden-only comparators; an exact finite-rate
  giraffe/naked-mole-rat contrast also requires residual modulation under the
  stated interval boxes. A broad prespecified input stress box can remove that
  requirement for some modeled step counts. All results are retained below.
* **An identification limit:** death fractions are not fixed-age onset risks;
  colonic passenger rates are not whole-animal functional driver clocks; mass
  is not an independent target count. Unrestricted calibration, targets,
  progression or observation factors produce exact endpoint aliases. No
  defensible number of additional biological steps is estimated.

**Novelty assessment:** a reproducible synthesis plus a **candidate narrow
mathematical contribution, priority unestablished**. The potentially useful
addition is the complete finite-rate two-species compatibility classification
with unknown absolute rate/target calibration, its sharp interval modulation
infimum, and the prospectively specified confrontation with accessible
published aggregates. The exact chain and elasticity are already in PR #70;
the interval-intersection geometry is classical and already used in PRs #64
and #70. The broad Peto explanation, rate/step tradeoff, mutation–cancer merge,
and nonidentifiability are prior work. This is not an established new biological
discovery. If the same composed certificate is found in prior literature,
the contribution should be downgraded to an independently reproducible audit.

## 2. Literature first: what was already known

The initial literature survey preceded the protocol freeze and fitting;
source and newer-literature checks continued on **2026-11-18,
America/New_York** without changing the frozen analyses. Searches combined `Peto`, `somatic mutation
rate`, `lifespan`, `multistage`, `Armitage Doll`, `mortality`, `prevalence`,
`identifiability`, `extra steps`, `Cagan`, `Vincze`, `Compton`, and `2025/2026`.
Publisher pages, available full text, supplementary text, institutional copies,
and publication-linked tables were checked. This is not a systematic review
or proof of bibliographic priority. Inaccessible source components are listed
in [data/PROVENANCE.md](data/PROVENANCE.md).

| Primary predecessor | Established result and consequence here |
|---|---|
| Peto et al. (1975), Peto (1977), Armitage–Doll [L1–L3] | Age/exposure and multistage reasoning precede this work. A step count is conditional on the stochastic and observation model. |
| Caulin & Maley (2011); Tollis, Boddy & Maley (2017) [L4,L5] | Evolution can alter targets, rates, tissue architecture and defenses. The expectation of stronger protection in large animals is not novel. |
| Abegglen et al. (2015); Sulak et al. (2016) [L6,L7] | Elephant TP53 expansion and damage responses support particular mechanisms. Gene-copy counts are not measurements of an equal-rate sequential stage count. |
| Nunney (2020) [L8] | Explicitly compares mutation-rate reduction, extra genetic controls, metabolism and immune policing. Its independent unordered-target law differs from the ordered Erlang law at finite rates. Quantifying an extra-step tradeoff is already prior art. |
| Caulin, Graham, Wang & Maley (2015) [L26] | Direct mathematical and cross-species gene analysis already compares reduced mutation/division rates with additional tumor-suppressor requirements. The broad quantitative contrast is not new here. |
| Cagan et al. (2022), especially Supplementary Note 2 [L9] | Establishes the colon-crypt mutation/lifespan association and explicitly treats multistage rate reduction and extra driver requirements. The discussion already separates the lifespan association from residual mass effects. |
| Vincze et al. (2022) [L10] | Reports zoo cancer mortality, life-history traits and relatedness-aware analyses. Mortality risk and cumulative incidence of cancer mortality are distinct source estimators; neither is automatically fixed-horizon first-cancer incidence. |
| Compton et al. (2023 preprint; 2025 Cancer Discovery) [L11,L12] | Already joins mutation rates with comparative cancer prevalence. Neoplasia includes benign disease; malignancy and fatal cancer are different endpoints. This project cannot claim the first merge or association. |
| **Zheng (2025 JITC)** [L13] | Already joins Cagan and Vincze and explicitly claims quantitative resolution from mutation-rate/burden correlations. Correlation alone does not check the target-cell factor or identify functional transition rates. We audit this interpretation; we do not claim exact replication of its inaccessible supplemental analysis. Its correction concerns funding. |
| Butler et al. (2025), the Compton/Butler PNAS exchange, Dujon et al. (2025) [L14–L17] | Later studies disagree about mass associations and emphasize uncertainty and endpoint choice. A positive slope does not itself demonstrate the magnitude predicted by an unsuppressed cell-count model; a zero slope is not assumed as universal empirical truth here. |
| Balmain (2023); Kan & Chen (2024) [L18,L19] | Mutation/tissue/observational caveats and alternative immune-dynamic models are established. This note does not classify every proposed resolution of Peto's paradox. |
| Firsanov et al. (2025) [L20] | Bowhead fibroblasts required fewer experimental oncogenic hits than human fibroblasts, alongside enhanced repair. This opposes a universal claim that larger animals necessarily have more transformation steps. Fibroblast assays are not population-level stage estimates. The correction does not change that conclusion. |
| Kapadia et al. (2025); de Manuel et al. (2026) [L21,L22] | Blood-stem-cell mutation scaling differs from the colon result. Transport of a colon clock to all-cancer mortality is an unsupported empirical bridge, not just an abstract caveat. |
| Matthews et al. (2025); Zhang et al. (2026); Heid et al. (2026) [L23–L25] | Gene-copy associations, cancer-locus conservation, and induced fibroblast mutation assays add comparative evidence but do not supply matched annual functional mutation clocks and onset risks. They do not identify this model's k. |

**Explicit pivot:** the broad resolution was already proposed and partly
quantified. We ask the narrower open question: *what conditional finite-risk
constraint survives unknown absolute calibration, and does the accessible
aggregate panel support that constraint under declared observation bridges?*
Theorems 4–5 answer the mathematical pair question; the empirical audit and
negative identification verdict answer the data question within its limits.

## 3. Frozen protocol, sources and analysis design

[PROTOCOL.md](PROTOCOL.md) was published in commit
[`b98d28e`](https://github.com/dovvnloading/open-cancer/commit/b98d28e5258c09132d197e26529ff4f58d18b8a9)
before the first empirical calculation. The branch starts at then-current
`main`,
[`9ded69c`](https://github.com/dovvnloading/open-cancer/commit/9ded69c3586188b173e81fedb7a79abbf4e7e91f).
This is a public pre-analysis commit, not a blinded study or independent
registration. The literature and printed source values were already visible.
No post-fit exclusion, model selection or protocol amendment was made.

The CSVs preserve published numerical strings and provenance on every row.
The Cagan columns and maximum-longevity traits are read from Compton's
historical, preprint-linked aggregate table; mortality, sex-specific mass and
adult remaining expectancy come from Vincze's publication-linked text table.
The full original Cagan spreadsheet panel was not obtained. Cagan's printed
giraffe/naked-mole-rat comparison is a separate rounded-input sensitivity.
Missing published intervals and sequencing sample sizes are explicit blanks.
Derived Wilson intervals use the published known-cause-death denominators.
The [source audit](data/PROVENANCE.md) describes licenses, exact versioned
locators, missingness and access limits. No browser file downloads, binary
spreadsheet extraction, plotted-coordinate digitization or individual records
were used.

The primary horizon is remaining adult expectancy; the alternative is maximum
longevity. Neither is the mortality-quantile lifespan in Cagan's main analysis,
nor a matched deterministic endpoint. Mean annual passenger rate and
zero-intercept regression slope are analyzed separately; slope confidence
limits are never assigned to the mean. A negative reported slope lower bound
is retained and yields a zero boundary/unbounded log uncertainty when used.
Species labels are joined only through their explicit published binomials;
subspecies labels are retained in source common-name fields.

The fit uses the *rare-lineage model* with free amplitude, log H = a + log M
+ k log(mL), alongside intercept, mass-only, burden-only and free-mass-exponent
comparators. Every preregistered k is shown. The exact pair certificate uses
the finite-rate model and is logically separate from those limiting fits.
LOSO refits use only training outcomes; order-held-out sensitivity keeps
related species together. No phylogenetic causal inference is claimed.
Common time units and amplitude factors do not identify absolute susceptible
targets or driver rates.

## 4. Models and proved results

### Common notation and observational bridge

For species i, m_i is the measured annual *genome-wide passenger* substitution
rate, L_i is a specified deterministic model horizon, and B_i=m_i L_i. They do
not independently measure the transition intensity mu_i. A model imposing
mu_i=c m_i assumes a common positive conversion factor c. Adult body mass
M_i does not measure target-cell number; N_i=nu M_i^a imposes a common target
normalization nu>0 and a specified exponent a. These assumptions must be
displayed whenever the fitted model is interpreted. Neither mortality fraction
nor the average/maximum longevity of a heterogeneous population automatically
equals the model's directly observed onset probability by L_i.

### Theorem 1 — exact chain and two different per-division conventions

Assumptions: (1) k is a positive integer; (2) N is a positive integer number of
initially clean, independent target lineages; allowing real N>0 instead defines
an effective-exposure survival law; (3) one lineage undergoes k sequential
irreversible transitions with independent exponential waits of common annual
intensity mu>0; (4) there is no branching, loss, inherited advanced state,
lineage dependence, reverse transition or competition; (5) the first completed
lineage is directly observed at or before deterministic L>0, without competing
death, diagnostic lag or ascertainment bias.

Let E_{k-1}(x)=sum_{j=0}^{k-1} x^j/j! and psi_k(x)=x-log E_{k-1}(x). Then

    x=mu L,
    R=1-exp(-H),
    H=N psi_k(x),
    psi'_k(x)=x^(k-1)/((k-1)! E_{k-1}(x)).

The function psi_k is continuous, strictly increasing from zero to infinity,
and psi_k(x)~x^k/k! as x decreases to zero. Hence H~N x^k/k! for fixed k and
x→0, even if N varies and H does not vanish. R~N x^k/k! additionally requires
H→0. The elasticity eta_k(x)=x psi'_k(x)/psi_k(x) is in [1,k]; for k>1 and
finite x>0 both inequalities are strict.

Proof. The time to k transitions is Erlang, or the kth event time of a Poisson
process, so one-lineage survival is exp(-x)E_{k-1}(x). Independent survival
multiplies across N lineages, giving H. Differentiation gives psi'. Since
E_{k-1}(x)=1+O(x), integration gives the rare expansion. The derivative psi' is
nondecreasing: its derivative has the sign of
sum_{j=0}^{k-1}(k-1-j)x^j/j!, strictly positive when k>1. Thus
psi(x)=integral_0^x psi'(t)dt<=x psi'(x). Also psi'(t)/t^(k-1) is nonincreasing,
so psi(x)>=psi'(x)integral_0^x(t/x)^(k-1)dt=x psi'(x)/k. These prove the
elasticity bounds, including strictness. Finally, a fixed-degree polynomial
has logarithm O(log x), so psi_k(x)→infinity. QED.

Unequal-rate extension. If step j instead has annual intensity mu_j>0,
write A=L^k product_j(mu_j)/k! and q=P(sum_j Exp(mu_j)<=L). Then

    exp(-L max_j mu_j) A <= q <= A,
    q <= -log(1-q) <= q/(1-q),
    H=-N log(1-q).

The first bound integrates the joint exponential density over the positive
simplex of volume L^k/k!, using sum_j(mu_j t_j)<=L max_j mu_j. The
second follows by integrating 1/(1-z) from zero to q. Therefore
H~N L^k product_j(mu_j)/k! when L max_j(mu_j)→0, with fixed k. These
are inequalities and a limit, not an exact product model at arbitrary rates.

Per-division variants require a further explicit clock convention:

* With independent Poisson divisions at rate d/year and independent success
  probability p in (0,1] at each division, at most one next-stage transition per
  division, thinning gives annual mu=dp and the exact formula above.
* With exactly D integer division opportunities per lineage and independent
  probability p in (0,1) of one next-stage transition at each opportunity,
  survival is sum_{j=0}^{min(k-1,D)} binom(D,j)p^j(1-p)^(D-j). Thus R is one
  minus its Nth power. For D<k, R=0. For fixed D>=k as p→0, the leading
  one-lineage completion probability is binom(D,k)p^k. The approximation
  (Dp)^k/k! requires the additional large-D Poisson regime. Calling a
  continuous division-clock intensity u a literal Bernoulli probability would
  conflate these experiments.

### Theorem 2 — mass scaling and what inverse-lifespan rates cannot cancel

Assumptions: Theorem 1; common k; positive differentiable N(M),mu(M),L(M),
and no unmodeled observation factor. Define a=dlogN/dlogM,
b=dlogL/dlogM, and gamma=-dlogmu/dlogM. At any finite positive H,

    dlogR/dlogM = [H/(exp(H)-1)] [a+eta_k(mu L)(b-gamma)].

Thus local constancy of R requires gamma=b+a/eta_k(mu L), reducing to
gamma=b+a/k in the rare-lineage limit. If mutation is parameterized by a
per-division probability p with annual mu=dp under Poisson divisions d,
and d~M^delta, p~M^-c_p, the rare constant-risk
condition instead is c_p=b+delta+a/k. For unequal annual steps with
mu_j~M^-c_j, the rare hazard exponent is a+kb-sum_j c_j.

Proof. Take logarithms of H=N psi_k(mu L), differentiate, then apply
dlog(1-exp(-H))/dlogH=H/(exp(H)-1)>0. The rare and unequal-rate formulas follow
from Theorem 1. QED.

In particular mu L=x_0 constant gives H=N psi_k(x_0) *exactly*. For a pair
with N_i/N_j=n, R_i=1-(1-R_j)^n, independently of k. Therefore inverse-lifespan
rate scaling does not cancel a growing target count under these assumptions.
Risk saturation can make slopes small near R=1; it never makes two distinct
positive finite N have exactly equal risk at fixed x_0.

### Theorem 3 — exact endpoint equivalence and branching boundary

Assumptions: Theorem 1's model and observation law for each species, observed
0<R_i<1 and positive N_i,L_i; transition rates may differ between species and
are unrestricted positive intensities. Then a complete endpoint equivalence
class is characterized by the scalar

    H_i=N_i psi_{k_i}(mu_i L_i).

For any proposed integer k_i and positive N_i,L_i, exactly one intensity fits:

    mu_i=psi_{k_i}^{-1}(H_i/N_i)/L_i.

Consequently arbitrary endpoint aggregates cannot identify k when rates are
free. Observing passenger m_i does not remove the construction if conversion
c_i=mu_i/m_i is free. If mu and L both are unobserved, their product
is additionally indistinguishable. In the rare product family the identified
combination is only N_i L_i^k product_j(mu_ij)/k!; treating this limiting
invariant as a finite-rate invariant is incorrect.

Proof. The endpoint probability is injective in H. Theorem 1 supplies the
strictly increasing bijection psi and therefore its unique inverse. QED.

Restricted cross-species power-model equivalence. If one additionally imposes
the rare-law model H_i=A M_i^a B_i^k with a single A>0, equality of fitted
endpoint vectors is equivalent to

    (alpha-alpha') + (a-a')logM_i + (k-k')logB_i = 0 for every i,
    alpha=log A.

Thus two admissible parameter vectors are equivalent exactly when their
difference lies in the null space of the matrix with rows (1,logM_i,logB_i);
k must still be a positive integer. A full-column-rank design identifies these model
coefficients; this does not establish their biological meaning. If all B_i
equal B0, only alpha+k logB0 and a are identified (provided masses vary). If
logB_i=u+v logM_i exactly and masses vary, only alpha+ku and a+kv are identified. These
claims follow by taking logarithms and the definition of a matrix null space.
When mu_i=c m_i and N_i=nu M_i^a, A=nu c^k/k!, so even identified A,k do
not separate c and nu. Allowing c_i=c M_i^-s changes the mass coefficient to
a-ks, confounding target scaling with species-dependent mutation calibration.
These are exact equivalences *inside the stated power model*, not exact
finite-rate equivalences. With finite rates and equal B, the identified
common amplitude instead is nu psi_k(cB); all k still alias when nu is free.

Branching extension, with different stated assumptions: a nondepleting normal
reservoir of size N creates independent clones at Poisson rate N mu_0;
one clone has first-completion CDF q_theta(t) under a specified nonexplosive
finite-rate branching model; observation is direct by L. Exact Poisson marking
gives H=N mu_0 integral_0^L q_theta(t)dt. Thus initiation intensity and
integrated clone success are confounded at one endpoint. If mu_0=c_0 m
and *all* clone rates are m times common dimensionless rates, with a common
clone initial type/law and transition structure across species, then
q_theta(t)=Q(mt) by time rescaling and

    H=N c_0 integral_0^(mL) Q(z)dz.

Equal mL still leaves exact proportionality to N, even with expansion. If
birth/loss rates do not share this rescaling, the conclusion does not follow.
This is the existing Poisson marking result from PRs #64/#70, not a new
classification of every multistage branching model.

### Theorem 4 — complete two-species compatibility without absolute calibration

Assumptions: (1) two species obey Theorem 1 with common k; (2) the common
passenger-to-transition calibration c>0 and target normalization nu>0 are
unrestricted, mu_i=c m_i and N_i=nu M_i^a, with fixed a and positive
M_i,m_i,L_i; (3) real effective N is allowed; (4) exact direct onset risks are
strictly between zero and one. Define

    B_i=m_i L_i, b=B_i/B_j,
    n=(M_i/M_j)^a,
    y=(H_i/H_j)/n, H_i=-log(1-R_i).

For k>1 and b!=1, a finite positive common (c,nu) fits both endpoints if and
only if

    min(b,b^k) < y < max(b,b^k).

For k=1, feasibility is exactly y=b. For b=1, feasibility is exactly y=1 for
every k. The closed bounds are sharp limiting bounds, but for k>1,b!=1 their
endpoints require limiting c, not finite c. Two-species feasibility does not
make pairwise conditions jointly sufficient for three or more species: every
pair must share the same c and nu.

Proof. Let x=c B_j. Ratios eliminate nu and give y=psi_k(bx)/psi_k(x).
Integrating the strict elasticity bounds over log x gives the stated strict
interval, with inequalities reversed when b<1. At x→0 the ratio tends to
b^k; at x→infinity it tends to b since psi_k(x)/x→1. Continuity and the
intermediate value theorem give every interior value; no monotonicity of eta
or of the ratio is assumed. Once x is found, take c=x/B_j and
nu=H_j/[M_j^a psi_k(x)]>0. The k=1 and b=1 cases follow identically. QED.

**Synthetic counterexample to joint sufficiency.** Set equal masses,
B=(1,2,3) and H=(1,3,24/5), with any common k>=2. The pairwise hazard ratios
3, 24/5 and 8/5 lie strictly between respectively (2,4), (3,9) and (3/2,9/4),
and thus within the wider-or-equal envelopes for every k>=2. Yet Theorem 1
proves psi_k strictly convex. A common c and nu would require
H_1+H_3>2H_2 at these equally spaced burdens; instead 29/5<6. Consequently
no global fit exists. These are exact mathematical fixture values, not
species observations; the verifier checks their rational inequalities.

For k>1, equal risk, n>1 and b<1, the condition is b^k<1/n<b, or equivalently
k>log(n)/log(1/b) together with nb>1. If b>=1, equal risk is impossible for
all common k. This is a conditional compatibility threshold, not an inferred
biological driver count. If physical constraints bound c or nu, the unrestricted
sufficiency argument must be replaced by their restricted feasible set.

### Theorem 5 — sharp pair interval bound on residual modulation

Assumptions: Theorem 4's law with a multiplicative cumulative-hazard factor
s_i>0, H_i=nu M_i^a s_i psi_k(c B_i). The factors are an effective exposure
parameterization, not automatically clone survival probabilities. The allowed
ratio observations y and b vary *independently* over closed positive finite
intervals [y_-,y_+] and [b_-,b_+]. Their use as a rectangular uncertainty set
is declared and is not a claim of a calibrated joint confidence region.

Put f_k(b)=min(b,b^k), g_k(b)=max(b,b^k), A=f_k(b_-), C=g_k(b_+).
The closure of feasible s_i/s_j is exactly

    [y_-/C, y_+/A].

The sharp infimum of the necessary pairwise fold modulation
F=max(s_i/s_j,s_j/s_i) is

    F_inf=max(1, y_-/C, A/y_+).

In particular, if y_+/A<1, species i must have s_i/s_j<=y_+/A (with strictness
when the limiting boundary is unattainable); the required reduction is at
least A/y_+. If instead y_-/C>1, the direction is increased exposure in i.
An allowance written s_i in [s0/rho,s0*rho] for both species has infimum
rho=sqrt(F_inf). These are different parameter conventions.

Proof. The functions f_k and g_k are increasing for positive b, including
across b=1. The map (b,x)→psi_k(bx)/psi_k(x) is continuous on the connected
domain [b_-,b_+]×(0,infinity), so its image is an interval. Theorem 4's bounds
and limiting corner values show its closure is [A,C]. The independent
positive y interval then gives exactly the stated ratio interval in closure.
Minimize max(s,1/s) on it to obtain F_inf. For two multipliers in a common
symmetric rho interval, their feasible ratio interval is [rho^-2,rho^2],
giving the square root. QED.

For independent input boxes, valid ratio intervals are obtained from
H(R)=-log(1-R), B=mL, and, for fixed a>=0,

    y_-=(H_i^-/H_j^+)(M_j^-/M_i^+)^a,
    y_+=(H_i^+/H_j^-)(M_j^+/M_i^-)^a,
    b_-=B_i^-/B_j^+, b_+=B_i^+/B_j^-.

A source interval including R=0 supplies H^-=0, so any reciprocal bound using
that endpoint is unbounded. Do not clip zero risks to an invented positive
value. A source interval including R=1 similarly gives an infinite upper
hazard. Report one-sided/vacuous bounds instead. If shared input errors impose
covariance, independent boxes are outer sensitivity sets and sharpness applies
to that larger set, not to an unknown smaller joint region.

### Theorem 6 — why an extra-step count needs absolute calibration

Assumptions: Theorem 1; a fixed known x=mu L>0 and N>0; only the integer
step count changes, retaining the same intensity at each old and new stage.
Then psi_{k+d}(x) is strictly decreasing in d>=0 and tends to zero. Thus for
any positive target hazard H_* a finite smallest d>=0 with
N psi_{k+d}(x)<=H_* exists (unless d=0 already satisfies it). This threshold
concerns being no higher than a target, not exact equality to an observed risk.
In the rare limit for fixed k,d,

    psi_{k+d}(x)/psi_k(x) ~ x^d k!/(k+d)!.

Proof. E_{k+d}(x)>E_{k+d-1}(x), so psi decreases strictly; the partial sums
increase to exp(x), hence psi decreases to zero. The ratio limit follows
from Theorem 1. QED.

Conversely, take two species with equal B, target ratio n>1 and equal positive
risk. For *each* fixed base k>=1 and *each* proposed extra count d>=1, there
are positive finite common c and nu, allowing real effective N, fitting that pair when the larger species
uses k+d: the continuous ratio psi_{k+d}(x)/psi_k(x) tends to zero as x→0 and
to one as x→infinity, so it attains 1/n. Choose nu to match the smaller risk.
The common-count d=0 model fails, yet the pair cannot distinguish one extra
step from any other positive extra count when absolute calibrations are free.
Bounds on c, actual target counts, stage-specific intensities and age-specific
incidence could change this conclusion. Unordered target hits have different
finite-rate laws and a different factorial convention; no numerical extra-step
count should be transferred between these models without specifying one.

### Theorem 7 — mortality is a different observation, even with independent death

Assumptions: Theorem 1 defines onset CDF F_C(t)=1-exp(-H_C(t)) for a specified
horizon L. Other-cause death is independent with survival S_D(t). Cancer onset
is instantly fatal (an intentionally restrictive bridge). The probability of
cancer death before other death and before L is

    P_C(L)=integral_0^L S_D(t) dF_C(t)
          =integral_0^L exp[-H_D(t)-H_C(t)] h_C(t)dt.

It generally is not F_C(L). Even restricting other death to Exp(delta) for
delta>0, with delta=0 interpreted as absence of competing death,
P_C(L;delta)=integral_0^L exp(-delta t)dF_C(t) is continuous and strictly
decreasing from F_C(L) at delta=0 to zero as delta→infinity. Hence every
probability in (0,F_C(L)] is compatible with the same onset process under
some such competing-death rate.

Proof. Condition on cancer onset time and independent survival to it. Dominated
convergence gives continuity and the limit; strict monotonicity follows from
positive onset probability at positive times. QED.

This construction treats L as an observation horizon, not as an additionally
fixed mean other-cause lifespan. If a source supplies a lifespan statistic,
one must model what statistic it is and its joint relation with the cancer
population; one summary does not directly supply S_D(t). Noninstantaneous
cancer fatality, prevalence in necropsies and selective postmortem ascertainment
introduce further observation parameters. Thus mortality-only endpoint data
do not inherit the direct-onset pair theorem without an explicitly justified
bridge.

Even a **known mean other-cause lifetime** does not determine the observation.
For the variant with no administrative horizon and independent positive
other-cause lifetime D, cancer-first probability is E[F_C(D)]. Fix E[D]=ell>0.
Deterministic D=ell gives F_C(ell)>0. For 0<epsilon<1, let D equal epsilon*ell
with probability 1-epsilon and ell*(1-epsilon+epsilon^2)/epsilon with
probability epsilon. Its mean is exactly ell; both values are positive, but
E[F_C(D)] <= F_C(epsilon*ell)+epsilon tends to zero. Continuity of the
clean-start CDF at zero proves the claim. This is a mathematical counterexample
to replacing a lifetime distribution by its mean, not a model of zoo survival.
Observed all-cause lifespan, known-cause sampling and cancer fatality introduce
additional constraints that would have to be modeled explicitly.

### Theorem 8 — common-amplitude compatibility in the rare family

**Assumptions:** finitely many species; a fixed positive integer k and the
specified power model H_i=A M_i B_i^k s_i, A>0, with positive finite inputs.
This is an exact model specification or Theorem 1's rare-lineage limit, not
an arbitrary finite-rate equality. Species-specific input sets independently
permit log residual z_i=log H_i-log M_i-k log B_i anywhere in [l_i,u_i].
Permit s_i in [1/rho,rho] about the common amplitude, with rho>=1. The boxes
are deterministic sensitivity sets, not simultaneous confidence regions.

The smallest feasible allowance is

    rho_min=exp(max(0, max_i l_i-min_i u_i)/2).

**Proof.** For a=log A, the equation z_i=a+log s_i permits a exactly in
[l_i-log rho,u_i+log rho]. All equations can hold iff these intervals have
nonempty intersection, equivalently max_i l_i-log rho <= min_i u_i+log rho.
Solving for the nonnegative log rho gives the formula. Closed bounded
intervals attain it. If independent positive input boxes are specified,
monotonicity gives

    l_i=log H_i^- -log M_i^+ -k log(m_i^+ L_i^+),
    u_i=log H_i^+ -log M_i^- -k log(m_i^- L_i^-).

All intermediate z values are attainable by continuity on the connected
input box. These are exact extrema for that rectangle; source covariance
could restrict it. **QED.**

This intersection result is established interval geometry, used in the prior
repository notes. In contrast, Theorems 4–5 retain the exact finite-rate shape
while eliminating absolute calibration for a pair. A failed power-family fit
alone is not a failed exact-chain fit. The exact pair violation supplies a
separate finite-rate obstruction.

### Boundary of the theorem claims

All clean-start results exclude inherited shortcuts, target dependence,
parallel pathways with different counts, changing tissue architectures,
unmodeled growth/loss and observation lags unless specifically introduced.
Exponential stage intensities are effective parameters, not genome-wide
passenger counts. The exact unordered-target alternative has single-lineage
completion (1-exp(-x))^k and a different finite law; its rare factor is x^k,
with the factorial absorbed by an unconstrained amplitude [L8; PR #70].
It is not interchangeable with the ordered chain when comparing extra steps
at a calibrated absolute rate. No theorem here proves a universal biological
driver count or classifies every multistage branching family.

## 5. Recomputed aggregate results


<!-- BEGIN COMPUTED RESULTS -->
### Coverage and observed endpoints

The readable mutation panel contains **9 species**; **6** overlap the mortality table, spanning **4 orders**, **75 cancer-attributed deaths** and **936 known-cause deaths**. These are selected published aggregates, not a random evolutionary sample.

CMR intervals below are **calculated Wilson 95% intervals**, not published intervals. Rates are mean passenger SBS/genome/year; L is adult remaining expectancy; M is the average of published female/male masses.

| Species | Annual SBS | L (years) | M (kg) | Cancer / known deaths | CMR (95% interval) |
|---|---:|---:|---:|---:|---:|
| Colobus guereza | 80.3605 | 12.7401 | 11.8848 | 3/79 | 0.03797 [0.01300, 0.10580] |
| Giraffa camelopardalis | 99.1124 | 12.5077 | 911.597 | 2/163 | 0.01227 [0.00337, 0.04363] |
| Heterocephalus glaber | 92.5653 | 12.9391 | 0.0478751 | 1/43 | 0.02326 [0.00412, 0.12059] |
| Lemur catta | 103.5028 | 14.7759 | 2.86501 | 17/251 | 0.06773 [0.04271, 0.10578] |
| Panthera leo | 159.8666 | 12.8088 | 161.944 | 26/204 | 0.12745 [0.08848, 0.18019] |
| Panthera tigris | 212.0010 | 11.6075 | 121.443 | 26/196 | 0.13265 [0.09216, 0.18727] |

Only **2** matched species have a printed zero-intercept-slope interval; **1** reaches or crosses zero. These intervals do not belong to the mean-rate estimator. Per-species sequencing counts and uncertainty for the lifespan proxies remain missing.

### Frozen prediction comparison

Every k is reported. RMSE is measured in log10 H, with H = −log(1−CMR). All amplitudes and free mass exponents are refit inside each fold. The burden-only comparator deliberately omits the target-cell factor.

| k | Cell-aware LOSO RMSE | Burden-only LOSO RMSE | Free-mass LOSO RMSE | Fitted free mass exponent | Cell-aware leave-order-out RMSE |
|---:|---:|---:|---:|---:|---:|
| 1 | 1.709621 | 0.345595 | 0.533678 | 0.005198 | 1.761197 |
| 2 | 1.757424 | 0.270995 | 0.440579 | -0.032174 | 1.835596 |
| 3 | 1.818334 | 0.281959 | 0.402915 | -0.069546 | 1.925116 |
| 4 | 1.891084 | 0.370976 | 0.435322 | -0.106918 | 2.027754 |
| 5 | 1.974365 | 0.497802 | 0.524980 | -0.144290 | 2.141626 |
| 6 | 2.066905 | 0.640351 | 0.648564 | -0.181662 | 2.265037 |
| 7 | 2.167519 | 0.790161 | 0.790317 | -0.219034 | 2.396515 |
| 8 | 2.275134 | 0.943778 | 0.942073 | -0.256406 | 2.534805 |

Intercept-only LOSO RMSE: **0.466363**; mass-only: **1.676044**.

The next table uses the preregistered reporting index k=3, not a selected biological count. NLL is an equal-species mean of per-known-death Bernoulli log loss. The delete-fold range describes sensitivity of paired squared-error differences from the intercept comparator; it is not a confidence interval.

| Panel | Cell RMSE | Burden RMSE | Free-mass RMSE | Cell NLL | Burden NLL | Cell minus intercept MSE; delete-fold range |
|---|---:|---:|---:|---:|---:|---:|
| primary | 1.818334 | 0.281959 | 0.402915 | 1.358136 | 0.234647 | 3.088845 [1.882570, 3.618648] |
| maximum_longevity | 1.875542 | 0.375877 | 0.517410 | 1.926715 | 0.239875 | 3.300163 [2.045720, 3.906306] |
| slope_estimator | 1.801606 | 0.298635 | 0.444036 | 1.348547 | 0.237625 | 3.028289 [1.890359, 3.537299] |

All fixed-k scores, coefficients, NLLs, fold predictions and sensitivity ranges are in `analysis_results.json`; no favorable k is hidden by this compact display.

### Compatibility and required residual modulation

Rare-family rho is a symmetric per-species hazard allowance about one common amplitude. Exact-pair F is the relative modulation fold between giraffe and naked mole rat, with CMR intervals and both-sex mass endpoints. These are different conventions: pairwise symmetric rho is sqrt(F). Finite-rate boundary values are infima, not necessarily attained minima.

| k | Rare rho, points | Rare rho, CMR intervals | Rare rho, factor-two stress | Exact-pair F, sex-mass/CMR boxes | Exact-pair F, slope-CI subset | Exact-pair F, factor-two stress |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 193.811 | 42.693 | 9.45454 | 1430.21 | 1536.63 | 89.3884 |
| 2 | 197.176 | 43.4343 | 2.40468 | 1430.21 | 1340.84 | 5.78246 |
| 3 | 200.599 | 44.1884 | 1 | 1430.21 | 1170 | 1 |
| 4 | 204.082 | 44.9557 | 1 | 1430.21 | 1020.93 | 1 |
| 5 | 207.626 | 45.7362 | 1 | 1430.21 | 890.847 | 1 |
| 6 | 211.231 | 61.6154 | 1 | 1430.21 | 777.34 | 1 |
| 7 | 254.598 | 88.3184 | 1 | 1430.21 | 678.296 | 1 |
| 8 | 364.937 | 126.594 | 1 | 1430.21 | 591.872 | 1 |

At k=3, the giraffe/naked-mole-rat effective hazard multiplier ratio has closure **[1.0015445e-06, 0.00069919611]**. This bounds a missing model factor under the proxy assumptions; it does not measure an immune mechanism or extra biological stages.

Cagan's separately printed rounded lifespan/rate pair gives burdens **2376** (giraffe) and **2325** (naked mole rat). Using these as horizons with Vincze's mortality/mass boxes remains a cross-source proxy; its k=3 relative-modulation infimum is **1412.12**.

The slope-CI scenario fixes rates lacking intervals at their slope point estimates; it is incomplete uncertainty propagation. The factor-two box is hypothetical. Compatibility inside a box is not confirmation; failure is not a calibrated sampling test. Pairwise success does not establish one common calibration for the entire panel.

### Outliers and lifespan sensitivity

The largest positive cell-aware k=3 LOSO residual is **Heterocephalus glaber**, 3.045178 log10 units (1109.63× observed/predicted H). The largest negative is **Giraffa camelopardalis**, -2.480413 (0.00330816×). These are model residuals, not species-wide rankings of cancer resistance.

In the **9-species maximum-longevity** panel, log annual rate on log longevity has slope **-0.926291**, with leave-species-out slope range **[-1.343776, -0.851907]**. LOSO log10-rate RMSE is **0.205713** for the free slope, **0.130172** for fixed inverse slope, and **0.348049** for intercept-only. This does not reproduce Cagan's original lifespan definition or complete species panel.

### Input integrity

SHA-256 uses newline-normalized UTF-8 content and POSIX relative path keys. Hashes establish local integrity, not truth or biological validity.

| Input | SHA-256 |
|---|---|
| `data/cagan_mutation_rates.csv` | `2cdb1239f394e6341f0d7b2ad11a5101f9a6b36b9d94bb5c2fdd45ba321e67ed` |
| `data/cagan_printed_pair.csv` | `de0b4c1462d9d2838ab17b78ad28aeed8bd7c70f1a359e777636be81706901fd` |
| `data/compton_life_history.csv` | `488e659de57fa22ff5b5101bf2d17dbd4595a80236fb0b06f9a28480e750fe97` |
| `data/vincze_mortality.csv` | `260e6465dd724e9033dfa9876cce56417624c18956e25fed01066ed6fda07a24` |
| `PROTOCOL.md` | `23cdd39e9034fd648cffa470423a5e5ac9425ab2cdef269b4ea8b82e5bcd9a81` |

<!-- END COMPUTED RESULTS -->

## 6. What the results establish, and one falsifiable prediction

The primary comparison supports a **conditional failure of the shared
cell-aware family**, not biological identification of its missing component.
All tested fixed-k cell-aware limiting models predict worse than the simple
intercept comparator. Their poor fit persists under the alternative horizon
and mutation estimator and under holding related species out together. The
burden-only comparator can predict better, but its specification has removed
the very target-cell factor that creates the mechanistic question. Its success
cannot establish that mutation-rate scaling quantitatively resolves Peto's
paradox in a model with more targets in larger animals.

The exact pair calculation closes a possible loophole: the discrepancy is not
merely a misuse of a rare-risk approximation. With primary point clocks,
giraffe has the larger modeled lifetime burden as well as the much larger
mass proxy. Theorem 4 therefore rules out the unmodulated common-k family
for **every positive common k**, under those point-input and observation
assumptions. The uncertainty calculations quantify the residual factor and
preserve their limitations. In particular, changing to slope estimates can
change the burden ordering, and the hypothetical rate-and-lifespan stress
box can restore compatibility for some tested k. Missing uncertainty is not
evidence of exact inputs.

The multiplier s scales cumulative hazard/effective exposure. It is not an
absolute probability of a clone surviving. In a rare family, a change in the
common per-stage driver rate could mimic a multiplier through its kth power;
a change in targets or progression could mimic it directly. In the exact
family these transformations differ and require the actual calibration.
Theorem 6 shows why assigning a numerical extra-step count without that
calibration would be arbitrary. Elephant and bowhead mechanism studies are
valuable evidence [L6,L7,L20], but neither animal is represented by a matched
annual-rate/mortality row in this permitted panel. They cannot be inserted as
synthetic empirical supports for a fitted count.

**Prospective falsifiable prediction.** Before inspecting new cancer outcomes,
select a species pair with independently measured susceptible target counts,
the same tissue and observation design, and equal cumulative **functional**
transition clocks. The shared-k, equal-rate, independent-target model predicts

    H_large/H_small = N_large/N_small,
    R_large = 1-(1-R_small)^(N_large/N_small).

This holds for every common k at finite rates. If clocks are only interval
matched, use Theorem 5's preregistered envelope instead of declaring exact
equality. Direct, age-matched first-onset estimates with an independently
specified competing-death/detection model that fall outside this envelope
falsify the shared mutation-scaling-only specification. Finding similar
hazards despite substantially different independently measured target counts
would challenge it directly. A changed clock, target count or observation
map defeats the premise; it is not confirmation of extra steps. The current
cross-source aggregates have not performed this discriminating test.

Age-specific onset curves, functional mutation targets and calibrated clone
success could then distinguish candidate mechanisms. PR #64 explains why
finite bins, inherited states, observation delays and fast extra transitions
still obstruct biological stage counting. No experimental treatment or
clinical recommendation follows from this prediction.

## 7. Hostile self-review and retained falsification attempts

| Attempt to defeat the conclusion | Outcome or correction retained |
|---|---|
| Find prior quantitative treatments of the broad claim. | Cagan's Supplementary Note 2, Caulin et al., Nunney, Compton and Zheng already cover central themes. Novelty was narrowed before fitting; no first-merge, first-tradeoff or established biological-discovery claim survives. |
| Require a genuinely new theorem rather than renamed algebra. | The exact chain, elasticity and interval-intersection ingredients are explicitly credited. Only the composed two-species classification and sharp interval infimum are candidate refinements; priority remains unestablished. |
| Treat the measured SBS/year as the per-driver rate. | Rejected. A common conversion is an assumption, and free species conversions reproduce every k by Theorem 3. |
| Transport colon rates to all cancer-causing tissues. | Rejected as identified biology. The later blood study supplies direct contrary tissue evidence [L21,L22]. |
| Rescue the homogeneous model by finite-risk saturation. | The exact pair certificate still fails for the stated inputs; the result does not rely only on poor rare-model scores. |
| Infer k from favorable training fit or a positive mutation association. | Rejected. All fixed-k predictions are reported, including worsening high-k scores. Free mass coefficients are model coefficients, not target measurements. |
| Claim endpoint aggregates can never identify any restricted parameter. | Rejected. Theorem 3 gives the design-rank condition under imposed shared power-model assumptions; nonidentifiability requires the explicitly free biological/observation components. |
| Let all stages expand and assume the same conclusion automatically. | Rejected. Only a separately stated common dimensionless branching clock preserves the equal-burden N factor; arbitrary branching kinetics are not classified here. |
| Read pairwise compatibility as joint panel compatibility. | Rejected. Pairwise failure is a valid obstruction, but success need not supply one common calibration for all species. |
| Treat closed pair-envelope boundaries as attained finite fits. | Corrected to closure/infimum language; equality cases and k=1 are separate. |
| Claim partial slope confidence intervals fully quantify rate uncertainty. | Rejected. Their estimand differs from mean rate; missing intervals remain unresolved, and the tiger interval reaches the nonpositive domain. |
| Ignore input uncertainty when claiming required suppression. | Rejected. The prespecified factor-two stress can remove compatibility failures; this limits the conclusion even though the box is not a confidence region. |
| Replace remaining lifespan or death fractions by fixed-age cancer risk without warning. | Rejected. Theorem 7 gives a competing-death construction and a fixed-mean counterexample; the fit is labelled a proxy throughout. |
| Increase apparent precision with the wrong mortality denominator. | Rejected. The known-cause-death denominator was checked against the printed CMR; live animals and unclassified deaths remain separate. |
| Make related species look like independent prediction cases. | LOSO is supplemented by order holdout. No unobserved phylogeny or population sampling design is fabricated. |
| Broaden source coverage through inaccessible spreadsheets or plotted values. | Rejected. The missing original species, sequencing n and lifespan uncertainty are recorded. The maximum-longevity analysis is explicitly not Cagan's original analysis. |
| Infer a universal extra-step mechanism from large-body protection. | Rejected. Bowhead transformation experiments are counterevidence to that universal story; TP53 copies do not count sequential bottlenecks. |
| Accept matching cached output as the sole numeric check. | Separate source transcription checks, independent fit formulas, high-precision survival probabilities, exact discrete recursions, interval-root checks and held-out perturbations corroborate the implementation. |
| Use finite numerical checks as a proof. | Rejected. General claims depend on the written assumptions and proofs. Finite checks only audit their implementation. |
| Ignore underflow or silently combine uncertainty scenarios. | Review caught extreme-log-hazard underflow/overflow and unintended mass widening in slope scenarios. The implementation now uses explicit stable branches and retains mean mass in slope scenarios as frozen. |

All changes remain inside this new note and its single new test file. The
protocol is unchanged after publication. Independent AI implementation checks
are not independent scientific or biological validation. The analysis remains
unreviewed; CI establishes repository/computational properties only.

## 8. Reproduction and verification

From the repository root, with Python and its standard library:

```bash
python research/petos-paradox-mutation-rates/verifier.py
python research/petos-paradox-mutation-rates/analysis.py
python -m unittest discover -s tests -p 'test_research_petos_paradox_mutation_rates.py' -v
```

To regenerate the complete machine-readable result and computed Markdown:

```bash
python research/petos-paradox-mutation-rates/analysis.py \
  --output research/petos-paradox-mutation-rates/analysis_results.json \
  --markdown /tmp/petos_results.md
```

The JSON contains all species predictions and all model/uncertainty scenarios,
not just selected tables. `Infinity` is an explicit string for an unbounded
interval endpoint; no source value is clipped to a convenient positive floor.
The tests recompute the full result tree, check the generated README block,
independently audit fitted quantities and source-denominator arithmetic,
and exercise exact-chain and uncertainty boundaries. Float comparisons allow
small platform differences; display precision is coarser. Every path key for
content hashes is built using `Path.as_posix()`, with line endings normalized
before hashing. The tests also check normalization under alternate newline
encodings. All tests run offline with no added dependencies.

The existing repository CI discovers the new test without workflow edits.
The publication-registry gate is unchanged. Neither it nor a passing test
constitutes approval of the scientific interpretation. The PR reports actual
CI status for its published head.

## 9. Primary references

Access-date label: 2026-11-18. The source table locators and reuse terms are
also carried on every CSV row. This list credits primary studies and authors'
reviews/essays; it is not a claim that every proposed mechanism is established.

**[L1]** Peto R et al. (1975). *Cancer and ageing in mice and men.* British Journal of Cancer. [doi:10.1038/bjc.1975.242](https://doi.org/10.1038/bjc.1975.242).

**[L2]** Peto R (1977). *Epidemiology, multistage models, and short-term mutagenicity tests.* In *Origins of Human Cancer*, Cold Spring Harbor. Historical chapter, discussed by [L4,L8]; no new quantitative data are transcribed from it.

**[L3]** Armitage P, Doll R (1954). *The age distribution of cancer and a multi-stage theory of carcinogenesis.* [doi:10.1038/bjc.1954.1](https://doi.org/10.1038/bjc.1954.1).

**[L4]** Caulin AF, Maley CC (2011). *Peto's Paradox: evolution's prescription for cancer prevention.* [doi:10.1016/j.tree.2011.01.002](https://doi.org/10.1016/j.tree.2011.01.002).

**[L5]** Tollis M, Boddy AM, Maley CC (2017). *Peto's Paradox: how has evolution solved the problem of cancer prevention?* [doi:10.1186/s12915-017-0401-7](https://doi.org/10.1186/s12915-017-0401-7).

**[L6]** Abegglen LM et al. (2015). *Potential mechanisms for cancer resistance in elephants and comparative cellular response to DNA damage in humans.* [doi:10.1001/jama.2015.13134](https://doi.org/10.1001/jama.2015.13134).

**[L7]** Sulak M et al. (2016). *TP53 copy number expansion is associated with the evolution of increased body size and an enhanced DNA damage response in elephants.* [doi:10.7554/eLife.11994](https://doi.org/10.7554/eLife.11994).

**[L8]** Nunney L (2020). *Resolving Peto's paradox: Modeling the potential effects of size-related metabolic changes, and of the evolution of immune policing and cancer suppression.* [doi:10.1111/eva.12993](https://doi.org/10.1111/eva.12993), [primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7428811/).

**[L9]** Cagan A et al. (2022). *Somatic mutation rates scale with lifespan across mammals.* [doi:10.1038/s41586-022-04618-z](https://doi.org/10.1038/s41586-022-04618-z). [Supplementary Note 2, printed pp. 7–8](https://api.repository.cam.ac.uk/server/api/core/bitstreams/09e7cb22-2566-45a8-b91d-036ddd0d6b6f/content), [institutional license record](https://doi.org/10.17863/CAM.83740).

**[L10]** Vincze O et al. (2022; online 2021). *Cancer risk across mammals.* [doi:10.1038/s41586-021-04224-5](https://doi.org/10.1038/s41586-021-04224-5). [Pinned aggregate table](https://github.com/OrsolyaVincze/VinczeEtal2021Nature/blob/d2f8e2dff8e64393ce3512f8b8490507fa9ae803/data.csv).

**[L11]** Compton ZT et al. (2023). *Cancer Prevalence Across Vertebrates.* Research Square preprint version 1. [doi:10.21203/rs.3.rs-3117313/v1](https://doi.org/10.21203/rs.3.rs-3117313/v1). [CC BY full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10350200/), [historical aggregate table](https://github.com/zacharycompton/cancerAcrossVertebrates/blob/3b4db10c7c38c63aeebfc9a48ce4bfa376dfd2bc/Mutation_Data.csv).

**[L12]** Compton ZT et al. (2025; online 2024). *Cancer Prevalence across Vertebrates.* Cancer Discovery. [doi:10.1158/2159-8290.CD-24-0573](https://doi.org/10.1158/2159-8290.CD-24-0573).

**[L13]** Zheng M (2025). *Super-high tumor mutational burden predicts complete remission following immunotherapy: from Peto's paradox to druggable cancer hallmark.* JITC. [doi:10.1136/jitc-2024-010486](https://doi.org/10.1136/jitc-2024-010486), [funding correction](https://doi.org/10.1136/jitc-2024-010486corr1). Only its comparative inference is discussed; no treatment claims are endorsed here.

**[L14]** Butler G et al. (2025). *No evidence for Peto's paradox in terrestrial vertebrates.* PNAS. [doi:10.1073/pnas.2422861122](https://doi.org/10.1073/pnas.2422861122).

**[L15]** Compton et al. (2025). *Paradoxical indeed.* PNAS. [primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12088420/).

**[L16]** Butler et al. (2025). *Reply to Compton et al.: Paradoxical misunderstanding.* PNAS. [doi:10.1073/pnas.2507210122](https://doi.org/10.1073/pnas.2507210122).

**[L17]** Dujon AM, Biro PA, Ujvari B, Thomas F (2025). *Towards a more robust comparative oncology: a Bayesian reanalysis of Peto's paradox and discussion of comparative cancer risk studies in vertebrates.* [doi:10.1098/rsos.250840](https://doi.org/10.1098/rsos.250840).

**[L18]** Balmain A (2023; online 2022). *Peto's paradox revisited: black box vs mechanistic approaches to understanding the roles of mutations and promoting factors in cancer.* [doi:10.1007/s10654-022-00933-x](https://doi.org/10.1007/s10654-022-00933-x).

**[L19]** Kan H, Chen Y (2024). *Revealing endogenous conditions for Peto's paradox via an ordinary differential equation model.* [doi:10.1007/s00285-024-02123-7](https://doi.org/10.1007/s00285-024-02123-7).

**[L20]** Firsanov D et al. (2025). *Evidence for improved DNA repair in the long-lived bowhead whale.* [doi:10.1038/s41586-025-09694-5](https://doi.org/10.1038/s41586-025-09694-5), [title/data-availability correction](https://doi.org/10.1038/s41586-025-09952-6).

**[L21]** Kapadia CD et al. (2025). *Clonal dynamics and somatic evolution of haematopoiesis in mouse.* [doi:10.1038/s41586-025-08625-8](https://doi.org/10.1038/s41586-025-08625-8), [figure-axis correction](https://doi.org/10.1038/s41586-025-09635-2). The correction concerns a duplicated axis label.

**[L22]** de Manuel M, Przeworski M, Spisak N, Stolyarova A (2026). *What sets the mutation rate of a cell type in an animal species?* PLOS Biology, Essay. [doi:10.1371/journal.pbio.3003799](https://doi.org/10.1371/journal.pbio.3003799).

**[L23]** Matthews et al. (2025). *Variable Gene Copy Number in Cancer-Related Pathways Is Associated With Cancer Prevalence Across Mammals.* [doi:10.1093/molbev/msaf056](https://doi.org/10.1093/molbev/msaf056).

**[L24]** Zhang et al. (2026). *Longer-lived mammals exhibit greater evolutionary conservatism and epigenetic modification at cancer-associated loci.* [doi:10.1038/s41467-026-76479-3](https://doi.org/10.1038/s41467-026-76479-3).

**[L25]** Heid et al. (2026). *Mutagen-induced somatic mutation rate in primary mammalian cells in relation to maximum life span.* Geromedicine. [doi:10.70401/geromedicine.2026.0023](https://doi.org/10.70401/geromedicine.2026.0023), [PubMed record](https://pubmed.ncbi.nlm.nih.gov/42245236/).

**[L26]** Caulin AF, Graham TA, Wang LS, Maley CC (2015). *Solutions to Peto's paradox revealed by mathematical modelling and cross-species cancer gene analysis.* [doi:10.1098/rstb.2014.0222](https://doi.org/10.1098/rstb.2014.0222), [primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4581027/).

**Repository predecessors:** [multistage incidence identifiability (PR #64)](../multistage-incidence-identifiability/), [stem-cell divisions and risk (PR #70)](../stem-cell-divisions-risk/), [research claim index](../README.md), [data policy](../../DATA_POLICY.md). These existing files are read and cited, not modified.
