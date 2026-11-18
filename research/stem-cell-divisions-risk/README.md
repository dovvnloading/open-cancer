# What a cross-tissue half-slope can identify about carcinogenesis

**Nonclinical theoretical and computational research. No medical advice.
Unreviewed mathematical work and aggregate-data reanalysis. No biological,
clinical, or independent peer-review validation is claimed.**

The original correlation is reproducible, but its interpretation depends on
which quantities are allowed to differ between tissues. We prove an exact
compatibility result for a specified sequential model: for every proposed
step count, invert each tissue's risk into its required per-division transition
intensity, then calculate the smallest permitted tissue variation that can fit
all observations. The result accounts for the fact that an uncertain stem-cell
number also changes the effective division clock `D/N`. A second result shows
that this uncertainty effect depends on whether the source's
developmental-count formula is imposed instead: it reverses for k>1, and
changes from invariant to decreasing at k=1.

The empirical conclusion is conditional and negative. Homogeneous limiting
models predict these aggregates poorly; adding the tested common expansion
factor does not improve them. The exact finite-rate compatibility calculation
also requires substantial tissue variation, including under the prespecified
uncertainty stress test. With unrestricted tissue rates, every step count is
compatible by construction. No biological step count, amount of clonal
expansion, environmental attributable fraction, or preventable fraction is
identified by this dataset.

**Novelty assessment:** a new reproducible analysis and a **candidate narrow
mathematical refinement; priority unestablished**. The broad question is
already substantially answered in the literature. The potentially new part is
the composition of an exact finite-risk inverse, the coupled observation map,
sharp uncertainty extrema, and the data-version/held-out compatibility audit.
The generic interval-intersection formula is classical and already appears
in this repository's PR #64. It is not claimed as new here.

## 1. Literature first, and the explicit pivot

The targeted search used primary papers, publisher records, public full text
and author-hosted manuscripts. The client-date access label is
**2026-11-17, America/New_York**; GitHub commit timestamps provide the actual
protocol-publication chronology. Search terms included `stem cell divisions`,
`slope`, `multistage`, `bad luck`, `clonal expansion`, `identifiability`,
`Erlang inverse`, `bounded mutation rate`, `measurement error`, and
`compatibility frontier`, including recent work. This was not a systematic
review, a claim of exhaustive recall, or a proof of priority.

| Primary predecessor | Established result and consequence for this contribution |
|---|---|
| Tomasetti & Vogelstein, Science 2015 [L1] | Cross-tissue association from heterogeneous lifetime-risk and stem-cell estimates. A correlation explains variation under a regression; it does not itself specify an intervention. |
| Tomasetti, Li & Vogelstein, Science 2017 [L2] | Extends the association across countries and uses a separate sequencing/epidemiology argument about mutation origins. Country replication does not provide independently measured stem-cell parameters in each country. |
| Wu et al., Nature 2016 [L3] | Shows that correlation can persist under substantial extrinsic contributions. Quantitative lower-envelope attribution additionally assumes a common intrinsic baseline; the nonidentifiability argument does not validate that assumption. |
| Noble, Kaltz & Hochberg 2015; Noble et al. 2016 [L4,L5] | Separates cell number and divisions per cell; examines saturation and anatomical groups; criticizes constant intrinsic risk per division and lower-envelope attribution. These ideas and the known low-risk intestinal outliers are prior art. |
| Little, Hendry & Puskin 2016 [L6] | Tests common-count multistage restrictions and smoking/radiation associations. Failure of its shared assumptions is not rejection of every multistage model. Its permissively reusable supplement supplies the numeric-data route here. |
| Tomasetti & Vogelstein's PLOS reply 2017 [L7] | Explicitly discusses slopes below one, developmental divisions, clonal expansion and heterogeneous driver counts. A new contribution cannot be that a half-slope fails to identify k. |
| Tomasetti et al. and Wu et al., Nature exchange 2017 [L8,L9] | Disputes assumptions behind intrinsic-risk bounds, including false attribution in constructed examples. Neither side's disagreement creates an identified environmental fraction. |
| Nowak & Waclaw 2017 [L10] | Gives one-step and two-step scaling with D and D²/N, and separates initiation from progression. Geometry, target genes, survival, selection and exposure already appear as alternative explanations of scatter. |
| Nunney & Thai 2020 [L11] | Tests total-division and evolutionary multistage models; distinguishes age from tissue slopes and corrects developmental division counting. Its evolutionary interpretation is a competing hypothesis, not a consequence proved here. |
| Zhang et al. 2018 [L26] | Uses an exact independent-target model, the developmental lineage proxy, and the same source-defined exclusions. Those choices and exact forward risk calculations are prior art. Its unordered targets and our ordered waiting stages have the same rare power after amplitude adjustment. |
| Altenberg 2015; Weinberg & Zaykin 2015; Sornette & Favre 2015 [L12–L14] | Square-root scaling, ERS limitations, ecological aggregation and the distinction between explained variance and preventability were already identified. |
| Little & Hendry 2017 [L15] | Fully stochastic stem/transit models include proliferation, loss and mutagen timing. It would be inaccurate to say the critiques uniformly ignored expansion. |
| Armitage–Doll; Ashley 1969; Frank, Iwasa & Nowak 2003 [L16,L18,L25] | Multistage age powers and the importance of division architecture are classical. A modeled transition, a sequenced driver and a rate-limiting bottleneck are different estimands. |
| Lopez-Garcia et al. 2010; McHale & Lander 2014 [L17,L19] | Stem-cell replacement, drift and lineage extinction show why total tissue divisions do not determine independent terminal-lineage histories. We did not locate a distinct Lopez-Garcia proof of the ecological slope rule, and do not invent one. |
| Rozhok, Wahl & DeGregori 2015; Liu et al. 2022 [L20,L21] | Input estimates, selection and measurement error already have substantive critiques and explicit models. Measurement-error sensitivity alone is not novel. |
| Zhang et al. 2024; Blohmer et al. 2025 [L22,L23] | A measured normal-tissue clonality index did not explain organ risk in the former study; the latter develops lineage-division measurements. Neither supplies the missing biological parameters for this aggregate analysis. |

**Pivot.** We ask which *bounded, explicitly stated* rate families remain
compatible with finite-risk aggregates when the observation map and uncertainty
coordinates are specified. Theorem 3 and Corollary 3.1 supply that certificate.
Theorem 4 tests identical-input categories without relying on a cross-tissue
regression. The other theorems make the model-to-slope and equivalence arguments
explicit and checkable; their underlying ideas are established.

Before this work, we read [research/README.md](../README.md),
[PR #64's identifiability note](../multistage-incidence-identifiability/),
[AGENTS.md](../../AGENTS.md) and [DATA_POLICY.md](../../DATA_POLICY.md).
PR #64 concerns complete hazards, observation filters and finite age bins. Here
each tissue supplies one endpoint and unknown biological parameters can differ
between tissues. Its onset-order identifiability theorem therefore does not
identify a biological count in this experiment. Its Theorem 7 is an explicit
predecessor for the final interval-intersection step below.

## 2. Prospective protocol and analysis commitments

[PROTOCOL.md](PROTOCOL.md) was published in commit
[`f5c1f4e`](https://github.com/dovvnloading/open-cancer/commit/f5c1f4e87bfd03c85eb015fa14e72a28d33844ad)
before any regression, empirical inverse-rate calculation or model selection.
The subsequent source audit found the Nunney–Thai counting correction;
[AMENDMENTS.md](AMENDMENTS.md) prospectively added that sensitivity in commit
[`c7d255d`](https://github.com/dovvnloading/open-cancer/commit/c7d255d7bfafe0ac720e1752a4b55058d136ecd3),
also before the first fit. These are publicly inspectable pre-analysis commits,
not a blinded study or independent trial registration. Published summaries and
familiar outliers were known beforehand. No outcome-driven exclusion or change
to the primary model/metric was made after fitting.

The historical replication retains all source categories. The primary panel
excludes the explicitly labeled FAP, Lynch, HPV16, HCV and smoker subgroups,
as specified before fitting. Their removal does not make remaining risks
intrinsic. Anatomical groups hold out overlapping subtypes, bone subsites and
duodenum/small-intestine rows together. Row-level and group-level prediction
scores are both reported; the primary score weights anatomical groups equally.

Fixed k values, the growth grid, training-only parameter selection, both clock
conventions, source variants, deterministic sensitivity boxes and bootstrap
rules are in the frozen protocol. No best-fitting index is called a biological
driver estimate. Growth fits use a proved limiting shape; the exact inverse
certificate uses the finite-rate chain. Their conclusions are kept separate.

## 3. Models and observation boundaries


Use k for the number of modeled irreversible transitions. It is not an identified
number of biological rate-limiting steps or sequenced driver variants. N is a
number/effective number of initially unaltered independent target lineages. v is
a cumulative lineage division clock, and u is a positive transition intensity
per unit of that clock. An exponential transition intensity per unit division
clock is not a literal independent Bernoulli mutation probability at each discrete
division. Neither N nor v is established as biologically observed by a table of
coarse stem-cell aggregates.

In the equal-rate chain all k waits share a time scale by assumption. Outside
that restricted family, strictly positive rates can still contain arbitrarily
fast transitions; calling every transition rate limiting requires additional
time-scale bounds. The aggregates do not supply those bounds. This note counts
model transitions and states where a claim about biological bottlenecks fails.

D=Nv is the independent-lineage model's total clock exposure. If an analysis
sets v=D_published/N, that is an effective equal-allocation surrogate, not a
theorem that published developmental tissue divisions are independent lineage
clocks. Published D=N(s+2)-2 includes development and does not give the
root-to-leaf cell histories or their dependence. A sensitivity v=s+log2(N) is a
balanced-development path proxy, also an assumption, not a correction known true.

Observed R is interpreted by the model as probability of a directly observed
first malignant event by a common, specified endpoint. Source lifetime risks
may instead mix competing mortality, ages, registries and diagnosis mechanisms.
The conversion H=-log(1-R) is always algebraically defined for 0<R<1; it becomes
an integrated biological hazard only under the observation assumptions.


### Theorem 1 — exact no-expansion chain, rare limit, and elasticity

Assumptions:

1. k is a positive integer and 0<N<infinity. For the literal independent-lineage
   process N is a positive integer; extending S_lineage^N to real N>0 defines
   the stated effective-exposure observation model.
2. Each lineage begins at stage zero; it makes k sequential, irreversible
   transitions with independent exponential waits of rate u>0 in division-clock
   units. There is no clonal reproduction, loss, reverse transition, inheritance,
   shared ancestry, variation between lineages, or competing cause.
3. All lineages have the same deterministic v>0; distinct lineages' processes
   are independent. Cancer is the first completed lineage, immediately observed.
4. For the rare limit k is fixed and x=uv tends to zero. A further rare-person
   condition H tends to zero is required only when replacing risk R by H.

Put E_{k-1}(x)=sum_{j=0}^{k-1}x^j/j! and
psi_k(x)=x-log E_{k-1}(x). Then

    R=1-exp(-H), H=N psi_k(uv),
    psi'_k(x)=x^(k-1)/[(k-1)! E_{k-1}(x)] >0,
    psi_k(x) ~ x^k/k!,
    eta_k(x):=x psi'_k(x)/psi_k(x) belongs to [1,k].

For k=1 eta=1 identically; for k>1 and finite x>0 both inequalities are
strict. psi is continuous, strictly increasing from zero to infinity, hence
has a unique continuous increasing inverse on (0,infinity).

Proof. The sum of k equal-rate independent exponential waits is Erlang. Its
survival at v is exp(-uv) E_{k-1}(uv), either from its density or from the
probability of at most k-1 events of a rate-u Poisson clock. Product survival
across N independent lineages yields R. Differentiating gives the displayed
positive derivative, since E_{k-1}-E'_{k-1}=x^(k-1)/(k-1)!.
Integration from zero gives psi(x) ~ x^k/k! because E_{k-1}(x)=1+O(x).
The derivative is nondecreasing: its derivative has sign
(k-1)E_{k-1}-x E'_{k-1}=sum_{j=0}^{k-1}(k-1-j)x^j/j!, nonnegative and
strictly positive for k>1. Thus psi(x)=integral_0^x psi'(t)dt <=x psi'(x),
which proves eta>=1. Also psi'(t)/t^(k-1)=1/[(k-1)!E_{k-1}(t)] is
nonincreasing. Thus for 0<t<x,
psi'(t)>=psi'(x)(t/x)^(k-1). Integrating yields
psi(x)>=x psi'(x)/k, proving eta<=k. Both become strict for k>1 because
E_{k-1} is strictly increasing. Finally psi=x-log E tends to infinity since
log of a fixed-degree positive polynomial is O(log x). QED.

General unequal-rate rare limit. With positive step intensities u_1,...,u_k,
let q be the probability of completing all transitions by v. Writing the
joint density of the independent waiting times over their simplex yields

 exp(-v max_j u_j) [v^k product_j u_j/k!] <= q <= v^k product_j u_j/k!.

The exponent bound is valid since sum_j t_j<=v implies sum_j u_jt_j<=v max u.
The simplex volume is v^k/k!. Also q<=-log(1-q)<=q/(1-q), obtained by
integrating 1/(1-t). Therefore for fixed k as v max u_j ->0,
H=N v^k product_j u_j/k! times (1+o(1)), even when N grows so H need not
vanish. A log-risk approximation additionally requires H->0.
This is a proven limit, not an equality at arbitrary mutation parameters.

#### Corollary 1.1 — Ordered stages and unordered targets share a limiting shape

**Assumptions:** N independent initially clean lineages with clock v as above,
but now each lineage needs k distinct targets, each independently hit at
intensity u per clock unit in any order. The first lineage with all targets
hit causes the observed event. There is no expansion, reverse mutation or loss.
Then a lineage completes with probability `q=(1-exp(-u*v))^k`, so

    R=1-(1-q)^N, H=-N log(1-(1-exp(-uv))^k) ~ N(uv)^k.

For a fixed positive endpoint the unique fitting intensity is

    u=-log(1-[1-exp(-H/N)]^(1/k))/v.

**Proof.** Independent target completion probabilities multiply within a
lineage, and independent lineage survival probabilities multiply across
lineages. Taking logarithms gives H; `1-exp(-x)~x` gives the rare limit.
The displayed inverse follows by taking a kth root and logarithm, each on
its positive domain. Thus an unrestricted tissue-specific rate again fits
every k. In the rare regime, choosing
`u_unordered=u_ordered/(k!)^(1/k)` equates its shape to Theorem 1's ordered
shape. A fitted common amplitude therefore makes their fixed-k limiting
predictions identical. Their exact finite-rate curves generally differ;
no equality at finite rates with a common rate transformation is asserted.
**QED.**

This unordered-target model is closely related to Zhang et al. [L26], whose
per-division probability r corresponds to the continuous intensity
`u=-log(1-r)`. The independent-target result and developmental clock are
credited as prior work. This corollary adds no unregistered empirical fit.


### Theorem 2 — slopes are directional or finite-cloud statements

Assumptions: Theorem 1; differentiable positive N(D),v(D),u(D), with D=Nv;
common k; H>0; no unknown observation multiplier. For the finite-cloud bound,
there are finitely many distinct D values and common u; when D increases both
N and v are nondecreasing. Ties in D may be combined or handled by zero
pair weights. For OLS an intercept is fitted and Var(log D)>0.

The exact local log-risk slope is

 d log R/d log D = sigma(H) [ a + eta_k(uv)(1-a+c) ],
 a=d log N/d log D, c=d log u/d log D,
 sigma(H)=H/(exp(H)-1).

For common u and 0<=a<=1 this lies between sigma(H) and k sigma(H).
For an ordered finite cloud under the assumptions above, its intercept OLS
slope beta_R satisfies

    sigma(max_i H_i) <= beta_R <= k.

Consequently half-slopes at low risks do not follow merely from risk
saturation within this ordered, homogeneous chain class. The result does not
apply to unrestricted clouds where N and v can trade off negatively.

Proof. Differentiating logR=log(1-exp(-H)) with respect to logH gives sigma.
The chain rule for logH=logN+log psi(uv), and logv=logD-logN, gives the
formula. The bound follows from Theorem 1. For two ordered tissues j>i,
write A=log(N_j/N_i)>=0 and B=log(v_j/v_i)>=0. Integrating the elasticity
bound over logv gives
A+B <=log(H_j/H_i)<=A+kB. Since log(D_j/D_i)=A+B, the cumulative-hazard
secant lies in [1,k]. Further, sigma is decreasing: its derivative has the
sign exp(H)(1-H)-1<0, because that numerator is zero at H=0 and has
derivative -H exp(H)<0. Integrating sigma between logH values bounds the
risk secant between sigma(max H) and k. Finally

 beta_R = sum_{i<j}(X_j-X_i)(Y_j-Y_i) / sum_{i<j}(X_j-X_i)^2

with X=logD, Y=logR; this is a convex weighted average of the pairwise
secants, with weights proportional to (X_j-X_i)^2. The identity follows
by expanding both pair sums to n times their centered covariance/variance.
It gives the claimed finite-cloud bound. QED.

Rare-limit OLS decomposition (no sign assumptions). If
H_i=A N_i v_i^k m_i, with A>0, arbitrary m_i>0 absorbing rate products,
expansion/progression etc, and D_i=N_i v_i, then exactly in this specified
power model

 beta_logH = 1+(k-1) Cov(logD,logv)/Var(logD)
                + Cov(logD,logm)/Var(logD).

The exact risk slope is this plus
Cov(logD,log[(1-exp(-H))/H])/Var(logD).
Both follow by taking logs and covariance. A slope one-half is compatible
with any k>=2 if geometry provides covariance ratio -1/[2(k-1)] in the
rare common-m limit. For example positive paths
v proportional D^[-1/(2(k-1))], N proportional D^[1+1/(2(k-1))]
produce H proportional D^(1/2). For k=1 fixed m gives logH slope exactly
one irrespective geometry, but covarying m or appreciable saturation can
reduce it. For every k a suitably varying m can give every real logH slope.
These are algebraic model statements, not claims measured tissue trends.

For the published developmental convention define D*=D+2=N(s+2).
For a model H=A N s^k m with **s>0**, the exact log identity is instead

 logH=logA+log(D+2)+k log s-log(s+2)+logm.

Regressing against logD therefore requires all of these terms. Replacing s
by D/N silently equates different clocks; it must be a declared surrogate.
At s=0 the log identity is undefined; the positive effective or developmental
proxy clocks used in this analysis are separately declared models, not an
extension of log(0).


### Theorem 3 — exact aggregate inverse and coupled-uncertainty frontier

Assumptions:

1. Theorem 1 observation law with specified k, and the effective clock v=D/N.
2. Finitely many tissue aggregates satisfy R_i in (0,1), N_i,D_i>0. Every
   tissue may have its own u_i. k is common only within each candidate class.
3. Interval version: independently for each tissue the permissible triples
   fill the rectangular set [R_i^-,R_i^+] x [N_i^-,N_i^+] x [D_i^-,D_i^+],
   with finite positive endpoints and 0<R_i^-<=R_i^+<1. This is a declared
   deterministic sensitivity set, not a sampling confidence region.
   N,D variation is coupled in v=D/N; treating v as a fourth independent
   interval would be a different experiment. The rectangle is not a claim
   that every perturbation respects source D=N(s+2)-2 with fixed s.
4. A shared unknown baseline u0>0 and multiplicative heterogeneity allowance
   rho>=1 mean u_i is constrained to [u0/rho,u0 rho].

The unique rate fitting a point triple is

    U_k(R,N,D) = (N/D) psi_k^{-1}(H/N), H=-log(1-R).

Thus any fixed k fits any finite positive aggregate vector exactly if all u_i
are free. At fixed R,D, U is nondecreasing in N, strictly increasing for k>1;
it increases strictly with R and decreases strictly with D. Its N elasticity is

    d log U/d log N = 1-1/eta_k(psi_k^{-1}(H/N)),

which belongs to [0,1-1/k]. Therefore the exact feasible rate range over a
rectangle is the full closed interval [ell_i,r_i], where

 ell_i = N_i^-/D_i^+ * psi_k^{-1}(H_i^-/N_i^-),
 r_i   = N_i^+/D_i^- * psi_k^{-1}(H_i^+/N_i^+).

There are choices within all rectangles compatible with one baseline and
allowance rho iff

    rho^2 >= max(1, (max_i ell_i)/(min_i r_i)).

The sharp minimum rho is the square root of that right-hand side; feasible
baselines are exactly [max_i ell_i/rho, rho min_i r_i]. For exact data,
the minimum total rate spread maxu/minu equals max U/min U, and the symmetric
baseline is sqrt(max U * min U). For k=1, U=H/D and N cancels exactly.

Proof. Theorem 1 gives a bijective psi on the positive half-line, so direct
inversion proves uniqueness and unrestricted-k compatibility. Implicit
logarithmic differentiation of psi(x)=H/N gives
(d log x)/(d log N)=-1/eta_k(x); add the prefactor N/D to obtain the stated
elasticity and sign. The signs for R,D follow from the increasing inverse and
positive prefactor. Those monotonicities show the rectangle's minimum and
maximum occur at the displayed corners. The rectangle is connected and U is
continuous; its image in the real line is an interval, so every intermediate
value is attainable. For each tissue [ell_i,r_i] intersects
[u0/rho,u0 rho] iff ell_i/rho<=u0<=rho r_i. Intersect these conditions over
all i. It is nonempty iff maxell/rho<=rho minr. This proves necessity,
sufficiency and sharpness, including endpoints. QED.

**Prior-art boundary.** The interval-intersection end of this proof is
classical and already in that directory. The potentially new cancer-specific
composition is (i) exact finite-risk inverse, (ii) derivative-controlled
N/D coupling, (iii) sharp aggregate rectangle corners, and (iv) empirical
rate-space compatibility across k. Novelty must be qualified accordingly.


#### Corollary 3.1 — Source-consistent uncertainty is a different experiment


For a specified algebraic convention D=N(c+L)-c with c>=0, N>0, L>=0 and
v=c+L-c/N>0, use Theorem 1 and let R,N,L vary in independent closed
rectangles with R strictly between zero and one, N>0, L>=0, and v>0
throughout. The inferred rate

 U(R,N,L)=psi_k^{-1}(H/N)/(c+L-c/N)

increases with R and decreases with N and L. Therefore the exact lower and
upper rate corners are respectively

 U(R^-,N^+,L^+) and U(R^+,N^-,L^-).

Proof. The numerator is strictly decreasing in N by the increasing inverse.
The denominator has derivative c/N^2>=0 in N and derivative one in L.
Both are positive on the rectangle, giving the signs; continuity and
connectedness supply every intermediate value. The R sign is as before. QED.

This is distinct from the R,N,D rectangle theorem: the N derivative changes
because the permissible covariance structure changes. In the source-consistent
experiment N changes both D and v through a prescribed formula; at fixed D it
changes v inversely. Thus an uncertainty frontier is meaningful only after
stating which observables may change independently. TV's c=2 and Nunney--Thai's
corrected c=1 are source conventions, not independently measured lineage trees.

### Theorem 4 — shared-input constraints and a finite-risk rate-ratio bound

Assumptions: Theorem 1, common k; two observed aggregates have identical
N,D (hence v) but R_b>R_a. There is no measurement error in these inputs for
the point statement. Nuisance rates u are allowed to differ.

A common-u model cannot fit both aggregates. The unique rate ratio obeys

    (H_b/H_a)^(1/k) <= u_b/u_a <= H_b/H_a.

For k>1 both pointwise inequalities are strict. For separated risk
intervals 0<R_a^-<=R_a^+<R_b^-<=R_b^+<1, write
H_x^\pm=-log(1-R_x^\pm). A valid uniform (not necessarily sharp) envelope is

    (H_b^-/H_a^+)^(1/k) <= u_b/u_a <= H_b^+/H_a^-.

Using lower H_b and upper H_a supplies **only** a lower bound; it
cannot be used for the upper bound. If the common N,D itself varies
jointly within uncertainty bounds, this envelope remains valid for
each common realization. It does not apply if N,D can differ between
the two observations.

Proof. Equal inputs and u give equal R by Theorem 1, contradiction. Let
x_b=u_bv>x_a=u_av. Integrating 1<=dlogpsi/dlogx<=k gives
log(x_b/x_a)<=log(H_b/H_a)<=k log(x_b/x_a). Rearrange and exponentiate.
Strictness follows from the strict elasticity bounds for k>1. QED.

This quantifies an irreducible within-input inconsistency but does not identify
which omitted mechanism differs. Same-input rows may represent exposures,
hereditary syndromes, histologies or overlapping populations; they are not
independent random tissue replications. Prediction evaluation must group such
rows or explicitly report how duplicates influence leave-one-out error.


### Theorem 5 — expansion equivalence at one endpoint

Assumptions:

1. N,d,T>0 are known constants; normal cells form a nondepleting reservoir.
   Initiations are a homogeneous Poisson process of rate lambda=N d u0,
   u0>0 an initiation intensity per normal division. When D=N d T is used,
   it is the model's reservoir exposure, not an asserted developmental formula.
2. An initiated clone follows any specified time-homogeneous, finite-rate,
   nonexplosive multistage branching model with k-1>=1 remaining transitions,
   strictly positive finite forward rates, arbitrary finite nonnegative birth
   and loss rates, independent clones, and no other observation mechanism.
   All start clean. Their first malignant-time CDF is q_theta(t).
3. The first malignant clone before T gives the observed event. No competing
   mortality, dependence, or unknown lag is inserted. Let J_theta(T)=int_0^T
   q_theta(t)dt; the assumptions ensure 0<J_theta(T)<T.
4. The scalar u0 is unrestricted unless an explicit upper bound is added.

Exactly

    H=N d u0 J_theta(T)=D u0 B_theta(T), B=J/T in (0,1).

Hence, at fixed known D,T the complete single-endpoint observable combination
is w=u0 B_theta(T). Two (u0,theta) pairs are observationally equivalent iff
this product agrees. Every feasible k and any chosen expansion configuration
can fit any 0<R<1 by u0=H/[D B_theta(T)]. In particular free tissue-specific
u0 leaves all these branching families saturated and indistinguishable.
If u0,theta,T are shared across tissues, logH vslogD has exactly unit slope,
and the risk elasticity is sigma(H), independent of k or expansion. If theta
contains rates proportional to tissue d, B is not shared, and that unit-slope
conclusion does not apply.

Proof. An initiation at time s independently succeeds by T with probability
q_theta(T-s). Independent marking of the Poisson process gives a Poisson
number of successful clones of mean lambda J. The probability of no such
clone is exp(-lambda J). Positivity of q at every positive time follows from
positive probability of the finite sequence of forward events with no losses;
q<1 at finite time follows from positive probability of no first event, so
0<J<T. The observable is an injective function 1-exp(-D w) of positive w;
this proves both directions of the equivalence and the fitted amplitude.
Shared w yields logH=logD+logw; differentiation gives the risk statement. QED.

If u0 must be a literal probability per Poisson normal division, impose
0<u0<=1; arbitrary fixed theta then supports only H<=D B_theta. This bound
must not be erased. When H<D, any finite k can still be fitted by choosing
sufficiently fast downstream transition rates, for which B approaches one,
and then setting u0=H/(DB)<1. A proof: a sum of k-1 exponentials with rates
L has expectation (k-1)/L ->0; dominated convergence of its CDF integral
gives B->1. Optional finite birth/loss rates do not invalidate existence if
chosen zero. This fast-stage construction does not show that all modeled
transitions are biologically rate limiting.


### Theorem 6 — a proved clonal-expansion limit and its cross-tissue shape

This is an explicit limit of first-event probabilities. It is not obtained by
silently equating an expected number of malignant cells to a probability.
The construction and rare-event logic have substantial precedent in PR #64;
the cross-tissue clock mapping is the application here.

Assumptions:

1. Fix k>=2, T,d,N>0 and finite strictly positive coefficients
   u_0,...,u_{k-1}. Let epsilon>0 tend to zero. For every epsilon the normal
   reservoir has effective size N_epsilon=N epsilon^(-k), and normal
   divisions occur at rate d per cell. Use its nearest lower integer for a
   literal cell-count model; this does not change the limit. It is
   nondepleting, with independent initiations at rate
   lambda_epsilon=N_epsilon d epsilon u_0.
2. Each initiation starts one stage-1 clone; there are k-1 more forward
   transitions. At intermediate stage j the forward transition intensity
   is epsilon u_j d. A nonterminal forward mutation creates a next-stage
   cell while retaining its parent (asymmetric mutation). The final
   transition is the first malignant mark. Birth adds one same-stage cell.
3. Exactly one of the k-1 intermediate stages has birth intensity g d,
   g>=0 finite; all other birth intensities and all loss intensities are
   zero. Which intermediate stage grows does not affect the limit. All
   clocks and clones are independent, no clone is initially present,
   and first malignancy is directly observed with no competing cause.
4. All parameters except epsilon and N_epsilon remain fixed as epsilon->0;
   convergence is on a fixed finite time interval [0,T]. Intensities are
   continuous division-clock hazards. If an external convention additionally imposes a unit upper bound on the
   mutation intensities, sufficiently small epsilon satisfies that bound;
   this does not turn the continuous-clock construction into a discrete
   Bernoulli-per-division model.

For a fixed actual reservoir N instead of N_epsilon, the same proof below
first gives the uniform asymptotic

    H_epsilon(T)=N epsilon^k A v^k Phi_k(gv)/k! + O(N epsilon^(k+1)).

Thus the limiting shape already follows by normalizing a fixed-reservoir rare
risk. Scaling the reservoir is needed only for a nonzero limiting risk. The
free amplitude used in fitting absorbs epsilon^k. The error constant depends
on the fixed rates, k and horizon; it is not a uniform bound over unbounded
biological rates.

Set A=product_{j=0}^{k-1}u_j, v=dT, and

 Phi_k(z)=k! sum_{n=0}^infinity z^n/(n+k)!
         =k integral_0^1 (1-t)^(k-1) exp(zt) dt.

Then the cumulative first-event hazard converges uniformly, with an O(epsilon)
error on [0,T], to

 H_0(T)= A N v^k/k! Phi_k(gv),  R_0=1-exp[-H_0(T)].

At g=0, Phi=1. At z>0, Phi is increasing and its elasticity
z Phi'_k(z)/Phi_k(z) lies strictly between zero and z. If growth is instead
a shared constant G per chronological year, the factor is Phi_k(GT); at a
shared observation horizon T it is absorbed exactly into the free amplitude.
Thus shared chronological expansion and no expansion are observationally
equivalent in this limiting cross-tissue family at that horizon. Expansion
proportional to tissue division rate creates the nonconstant factor Phi_k(gv).

Proof. Let q_j^epsilon(t) denote the probability one stage-j founder has
produced a malignant mark by clone age t, and put q_k=1. Conditioning on its
first event gives the exact backward equations

 q_j'=epsilon u_j d q_{j+1}+alpha_j q_j-alpha_j q_j^2
                                      -epsilon u_j d q_j q_{j+1}, q_j(0)=0,

where exactly one alpha_j=gd and the others are zero. Probability bounds keep
the solution in [0,1], and finite linear per-cell rates ensure nonexplosion.
Dropping the nonpositive terms and using variation of constants first at the
last stage gives q_{k-1}=O(epsilon) uniformly on [0,T]. Inducting backward
with

 q_j(t)<=epsilon u_j d integral_0^t exp[alpha_j(t-s)] q_{j+1}(s) ds

gives q_j=O(epsilon^(k-j)). Define f_j by

 f_{k-1}'=u_{k-1}d+alpha_{k-1}f_{k-1},
 f_j'=u_j d f_{j+1}+alpha_j f_j  (j<k-1), f_j(0)=0.

At the last stage both omitted nonlinear terms have order epsilon^2.
At stage j<k-1, q_j^2 and epsilon q_j q_{j+1} have order
O(epsilon^(2(k-j))), which is O(epsilon^(k-j+1)) because k-j>=1.
The induction error from epsilon q_{j+1} has the latter order too. A second
variation-of-constants induction therefore establishes, uniformly on [0,T],

 q_j^epsilon=epsilon^(k-j) f_j+O(epsilon^(k-j+1)).

The exact Poisson marking argument of Theorem 5 gives
h_epsilon=lambda_epsilon q_1^epsilon. Thus
h_epsilon=N u_0 d f_1+O(epsilon); integrating gives uniform convergence
of cumulative hazards and survival functions. Rounding N_epsilon to an
integer creates only O(epsilon^k) relative reservoir error, also harmless.

For s>gd, the Laplace transforms of the linear f recursion yield

 Laplace(H_0)(s)= A N d^k / [s^2 product_{j=1}^{k-1}(s-alpha_j)]
               = A N d^k / [s^k(s-gd)].

Expanding 1/(s-gd) and inverting, or convolving the nonnegative exponential
kernels directly, gives H_0=A N d^k sum_{n>=0}(gd)^n T^(n+k)/(n+k)!.
The series converges uniformly on every bounded T interval and equals the
stated Phi expression. Its integral representation follows by expanding
e^{zt} and using the elementary beta integral
int_0^1 t^n(1-t)^(k-1)dt=n!(k-1)!/(n+k)!.
Its logarithmic derivative is the expectation of t under the positive density
proportional to (1-t)^(k-1)e^{zt}, hence lies strictly between zero and one;
multiplying by z proves the elasticity bound. Replacing gd by G in the same
calculation gives Phi_k(GT), independent of tissue when T is common. QED.

Cross-tissue slope implication. In the limiting family H=A N v^k Phi(gv)/k!,
for common A,g and D=Nv, let a=dlogN/dlogD and
zeta(z)=zPhi'(z)/Phi(z). The exact directional risk slope is

 beta_local=sigma(H) {a+[k+zeta(gv)](1-a)}.

For 0<=a<=1, positive division-scaled clonal expansion increases the
cumulative-hazard sensitivity to v; it does not by itself create a subunit
low-risk slope. Tissue-specific growth, negative N-v allocation, altered
mutation rates, progression, or observation effects can change that statement.

Limit-to-data boundary. Actual reservoir counts in the sequence are N_epsilon,
not N. For a common epsilon and k across tissues,
log N_epsilon=log N-k log epsilon and
log D_epsilon=log(Nv)-k log epsilon; these are common shifts absorbed by a
free intercept in cross-tissue slope/shape fitting. The theorem justifies the
limiting functional shape. Fitting that shape to the reported aggregates does
not prove that finite biological parameters are close enough to the limit,
that v=D_published/N is a true lineage clock, or that clonal expansion has been
identified. Do not present the asymptotic shape as an exact finite-rate risk.


### Theorem 7 — extrinsic attribution is observationally confounded

Assumptions: For each tissue i use Theorem 1's **equal-rate** sequential
chain: all k transitions in that tissue have the same positive total
intensity u_i=u_i^I+u_i^E, with nonnegative independent intrinsic and
extrinsic Poisson hazards. The split is constant across steps *within*
a tissue but may differ between tissues. Only (R_i,N_i,D_i) is observed;
there are no cause labels, known interventions, exogenous exposure
constraints or mutation spectra. Arbitrary step-specific intrinsic splits
are not covered by the displayed counterfactual formula.

For every alpha_i in [0,1] (allowing zero intensity by continuous extension), the split
u_i^I=alpha_i u_i, u_i^E=(1-alpha_i)u_i has exactly the same observed law.
If removing extrinsic hazards changes only the latter to zero, its model
counterfactual risk is

    R_i^I=1-exp[-N_i psi_k(alpha_i u_i v_i)].

As alpha_i ranges from zero to one, the corresponding preventable fraction
1-R_i^I/R_i ranges continuously from one to zero, despite identical observed
aggregates. A shared alpha also leaves all baseline aggregates unchanged.
In the rare-person, rare-lineage limit R^I/R approaches alpha^k, showing
why a mutation fraction and a case-preventable fraction are different.

Proof. Superposition of independent exponential/Poisson hazards gives total
transition intensity u_i^I+u_i^E=u_i, so every state-transition law and thus
risk is unchanged. Removing one hazard gives alpha_i u_i and Theorem 1's
counterfactual formula. Continuity and strict monotonicity of psi supply the
range; the rare ratio follows from psi(x)~x^k/k!. QED.

If an external convention requires a common intrinsic intensity b across
all tissues, every b in [0,min_i u_i] still produces the same baseline
observations with u_i^E=u_i-b. The set of counterfactuals is narrower than
independent alpha_i but remains unidentified. Known positive intrinsic bounds
or mechanistic labels would be additional information. Neither an OLS
correlation nor this constructive theorem supplies such information.



## 4. Source data, reuse and reproducibility boundaries

The CSVs contain published tissue-level aggregates and explicitly derived
arithmetic only. No person-level records, patient linkage, genomes or
controlled-access material are included. [data/PROVENANCE.md](data/PROVENANCE.md)
records source files, row/column locators, terms, transcription checks and
version differences. Every data row carries its own provenance and reuse basis.

| CSV | Meaning and permitted reuse |
|---|---|
| `tv2015.csv` | Original Table S1 numeric fields, checked against the printed table. Redistribution uses the CC0-public-domain factual republication in Little et al. S2; the original developmental-only maintenance fields are restored and attributed explicitly. |
| `little2016.csv` | Literal later S2 values. Its changed maintenance-division convention is preserved as a separate version. |
| `tv2015_nunney2020.csv` | Our derived sensitivity using the original parameters and `D=N*(L+1)-1`, following Nunney–Thai's published counting correction. Inputs are CC0; the calculation is our own. Their all-rights-reserved main article is cited, not copied or relicensed. |
| `liu2022_ranked_risk_audit.csv` | A separately labeled endpoint audit of the permissively reusable later matrix related to the international analysis. Its columns are independently sorted risks, not matched country observations. It supplies no division vector. |

The original Science correction notice concerns a figure label and wording;
it is not evidence that Table S1 numbers were formally replaced. Popular
online transcriptions contain thyroid and esophageal errors. This contribution
uses the checked primary table and the permitted republication instead.
Annual division rates and reported lifetime maintenance divisions are retained
separately, because multiplying the former by an assumed lifespan does not
always reproduce the latter. The original D is also rounded; the corrected
dataset recomputes its formula rather than simply subtracting N from a rounded D.

**TV2017 coverage limitation.** A coherent, permissibly reusable international
risk-plus-divisions panel could not be obtained: the publisher supplement
request was blocked, the PMC download did not yield the attachment, and the
accessible later supplements did not include the required division vector.
The later risk-only audit is included, with its truncated age endpoint and
independent column sorting exposed. It is not a refit of the country-level
TV2017 correlations. We did not join old stem-cell inputs to those risks:
the later paper also changes selected inputs, so that would be a new hybrid
dataset. This part of the user's desired empirical coverage remains incomplete.
No access control was bypassed or uncertain data license assumed.

The source risks are heterogeneous summaries, not a random sample of independent
organs. There are no sampling denominators or justified error models for N and
D. The bootstrap describes resampling of the selected groups. The stress box
is a deterministic, hypothetical allowance; it is not a confidence region and
does not enforce the source formula. Corollary 3.1 proves why changing that
uncertainty geometry changes the answer.

<!-- BEGIN RECOMPUTED RESULTS -->
### Published-data replication and panels

All numbers in this block are regenerated from the cited CSVs by `analysis.py`.
Errors and coefficients use base-10 logs. Group-held-out RMSE is in log H;
the primary score gives each anatomical group equal weight.

| Panel | Rows | Groups | log R slope | Correlation |
|---|---:|---:|---:|---:|
| Historical all categories | 31 | — | 0.532638 | 0.803929 |
| Primary | 25 | 15 | 0.436748 | 0.808171 |
| All finite risks | 30 | 15 | 0.503508 | 0.793788 |
| Later reanalysis, primary selection | 25 | 15 | 0.436748 | 0.808171 |
| Corrected development count, primary selection | 25 | 15 | 0.424927 | 0.805922 |

Historical R²: **0.646301**. Conditional group-bootstrap slope intervals:
historical: **[0.313519, 0.636103]** (2000 draws; 0 degenerate draws omitted).
primary: **[0.232630, 0.535653]** (2000 draws; 0 degenerate draws omitted).

These intervals resample this selected tissue collection; they do not quantify
stem-cell measurement error or epidemiological sampling precision.

### Prediction under the primary effective clock

| Model | Training RMSE | Leave-row-out RMSE | Leave-tissue-group-out RMSE | Group CV, row weighted |
|---|---:|---:|---:|---:|
| Descriptive log D | 0.601113 | 0.651007 | 0.647456 | 0.685612 |
| Separable log N, log v | 0.591326 | 0.667872 | 0.692403 | 0.707153 |
| Sequential k=1 | 1.208297 | 1.258642 | 1.254919 | 1.337461 |
| Sequential k=2 | 2.149402 | 2.238960 | 2.241784 | 2.388526 |
| Sequential k=3 | 3.204373 | 3.337888 | 3.389657 | 3.552807 |
| Sequential k=4 | 4.290013 | 4.468764 | 4.579012 | 4.747731 |
| Sequential k=5 | 5.387815 | 5.612307 | 5.784317 | 5.954880 |
| Sequential k=6 | 6.491611 | 6.762095 | 6.997335 | 7.168080 |
| Sequential k=7 | 7.598791 | 7.915407 | 8.214650 | 8.384705 |
| Sequential k=8 | 8.708062 | 9.070898 | 9.434599 | 9.603454 |
| Expansion k=2 | 2.149402 | 2.238960 | 2.241784 | 2.388526 |
| Expansion k=3 | 3.204373 | 3.337888 | 3.389657 | 3.552807 |
| Expansion k=4 | 4.290013 | 4.468764 | 4.579012 | 4.747731 |
| Expansion k=5 | 5.387815 | 5.612307 | 5.784317 | 5.954880 |
| Expansion k=6 | 6.491611 | 6.762095 | 6.997335 | 7.168080 |
| Expansion k=7 | 7.598791 | 7.915407 | 8.214650 | 8.384705 |
| Expansion k=8 | 8.708062 | 9.070898 | 9.434599 | 9.603454 |

Separable coefficients: **b_N=0.505770, b_v=0.358994**.
Geometry covariance ratio Cov(log D,log v)/Var(log D): **0.447708**.

| k | Fixed-k predicted log H slope | Fixed-k predicted log R slope | Exact-chain discrepancy, max log H | Expansion full-fit grid c |
|---:|---:|---:|---:|---:|
| 1 | 1.000000 | 0.975448 | 0.000000 | — |
| 2 | 1.447708 | 1.307569 | 0.000224 | 0.00 |
| 3 | 1.895415 | 1.561173 | 0.017310 | 0.00 |
| 4 | 2.343123 | 1.796491 | 0.163340 | 0.00 |
| 5 | 2.790831 | 2.032943 | 0.643030 | 0.00 |
| 6 | 3.238538 | 2.268747 | 1.571455 | 0.00 |
| 7 | 3.686246 | 2.504312 | 2.861217 | 0.00 |
| 8 | 4.133954 | 2.739604 | 4.365771 | 0.00 |

Maximum primary risk: **0.300000**; its attenuation factor is **0.832242**.
There are **113** pairs violating the joint N/v ordering condition.
The ordered-cloud theorem therefore cannot be used as a bound on this observed cloud.

### Exact rate-heterogeneity certificate

Here rho bounds each tissue intensity within [u0/rho,u0·rho]. It is not an
environmental fraction, a measured mutation-rate spread, or a confidence interval.

| k | Exact point rho | Rare-limit point rho | Stress-box rho | Upper-rate constraint | Lower-rate constraint | Compatible point allowances |
|---:|---:|---:|---:|---|---|---|
| 1 | 122.163202 | 122.163202 | 48.856287 | gallbladder | small_intestine | none |
| 2 | 58.096175 | 58.096128 | 18.371042 | ovarian_germ | small_intestine | 100 |
| 3 | 50.530540 | 50.528898 | 14.775186 | ovarian_germ | small_intestine | 100 |
| 4 | 47.132186 | 47.123349 | 13.253410 | ovarian_germ | small_intestine | 100 |
| 5 | 45.214277 | 45.191207 | 12.421121 | ovarian_germ | small_intestine | 100 |
| 6 | 46.117024 | 46.184780 | 12.471227 | ovarian_germ | colorectal | 100 |
| 7 | 47.129687 | 47.232471 | 12.608637 | ovarian_germ | colorectal | 100 |
| 8 | 47.895104 | 48.033808 | 12.711528 | ovarian_germ | colorectal | 100 |

### Source and clock sensitivity

| Panel / clock | k | Sequential group CV RMSE | Exact point rho | Stress-box rho |
|---|---:|---:|---:|---:|
| All finite / D/N | 1 | 1.180183 | 122.163202 | 48.856287 |
| All finite / D/N | 2 | 2.112594 | 58.096175 | 18.371042 |
| All finite / D/N | 3 | 3.238090 | 50.530540 | 14.775186 |
| All finite / D/N | 4 | 4.411276 | 47.132186 | 13.253410 |
| All finite / D/N | 5 | 5.602271 | 45.214277 | 12.421121 |
| All finite / D/N | 6 | 6.801726 | 46.117024 | 12.471227 |
| All finite / D/N | 7 | 8.005840 | 47.129687 | 12.608637 |
| All finite / D/N | 8 | 9.212785 | 47.895104 | 12.711528 |
| Original / L+log2 N | 1 | 1.080243 | 104.436621 | — |
| Original / L+log2 N | 2 | 1.651110 | 26.257759 | — |
| Original / L+log2 N | 3 | 2.411465 | 17.641538 | — |
| Original / L+log2 N | 4 | 3.230159 | 15.584229 | — |
| Original / L+log2 N | 5 | 4.072156 | 14.473423 | — |
| Original / L+log2 N | 6 | 4.925519 | 14.423579 | — |
| Original / L+log2 N | 7 | 5.785221 | 14.519314 | — |
| Original / L+log2 N | 8 | 6.648804 | 14.590717 | — |
| Later / D/N | 1 | 1.254919 | 122.163202 | 48.856287 |
| Later / D/N | 2 | 2.241784 | 58.096175 | 18.371042 |
| Later / D/N | 3 | 3.389657 | 50.530540 | 14.775186 |
| Later / D/N | 4 | 4.579012 | 47.132186 | 13.253410 |
| Later / D/N | 5 | 5.784317 | 45.214277 | 12.421121 |
| Later / D/N | 6 | 6.997335 | 46.117024 | 12.471227 |
| Later / D/N | 7 | 8.214650 | 47.129687 | 12.608637 |
| Later / D/N | 8 | 9.434599 | 47.895104 | 12.711528 |
| Later / L+log2 N | 1 | 1.079091 | 104.436621 | — |
| Later / L+log2 N | 2 | 1.644235 | 26.257759 | — |
| Later / L+log2 N | 3 | 2.398887 | 17.641538 | — |
| Later / L+log2 N | 4 | 3.212148 | 15.584229 | — |
| Later / L+log2 N | 5 | 4.048853 | 14.473423 | — |
| Later / L+log2 N | 6 | 4.897000 | 14.423579 | — |
| Later / L+log2 N | 7 | 5.751529 | 14.519314 | — |
| Later / L+log2 N | 8 | 6.609966 | 14.590717 | — |
| Corrected development / D/N | 1 | 1.306829 | 124.856154 | 49.939249 |
| Corrected development / D/N | 2 | 2.386528 | 82.146342 | 25.976132 |
| Corrected development / D/N | 3 | 3.619474 | 71.448749 | 20.891694 |
| Corrected development / D/N | 4 | 4.891117 | 66.643573 | 18.739946 |
| Corrected development / D/N | 5 | 6.177606 | 63.931704 | 17.563112 |
| Corrected development / D/N | 6 | 7.471277 | 65.224907 | 17.638489 |
| Corrected development / D/N | 7 | 8.768950 | 66.657152 | 17.832832 |
| Corrected development / D/N | 8 | 10.069079 | 67.739708 | 17.978355 |

### Prespecified outliers

These are group-held-out residuals of the descriptive log-H regression.
The ratio is observed H divided by predicted H; it is a required amplitude
adjustment, not a diagnosed mechanism or a success probability.

| Direction | Tissue / category | Held-out residual | H observed / H predicted |
|---|---|---:|---:|
| Positive | Basal cell carcinoma | 1.194788 | 15.659864 |
| Positive | Thyroid papillary follicular carcinoma | 1.167704 | 14.713099 |
| Positive | Gallbladder non papillary adenocarcinoma | 0.993609 | 9.853929 |
| Negative | Small intestine adenocarcinoma | -1.426125 | 0.037486 |
| Negative | Pancreatic endocrine islet cell carcinoma | -1.111434 | 0.077369 |
| Negative | Duodenum adenocarcinoma | -1.039508 | 0.091305 |

### Identical-input contrasts

All distinct-risk pairs with identical published N and D are retained.
An endpoint printed as one is not clipped to create a finite inverse.

| Lower-risk category | Higher-risk category | Risk ratio | H ratio | Required intensity ratio, k=2 |
|---|---|---:|---:|---:|
| aml | cll | 1.268293 | 1.268993 | 1.126496 |
| colorectal | colorectal_fap | 20.833333 | undefined at R=1 | undefined at R=1 |
| colorectal | colorectal_lynch | 10.416667 | 14.091151 | 3.753895 |
| colorectal_lynch | colorectal_fap | 2.000000 | undefined at R=1 | undefined at R=1 |
| duodenum | duodenum_fap | 116.666667 | 118.739444 | 10.897203 |
| head_neck | head_neck_hpv16 | 5.750000 | 5.949533 | 2.439212 |
| liver | liver_hcv | 10.000000 | 10.335885 | 3.214952 |
| lung_nonsmokers | lung_smokers | 18.000000 | 18.728657 | 4.327675 |

The version audit compares **155** numeric fields and finds **6** changes:

| Tissue | Field | Original | Later reanalysis |
|---|---|---:|---:|
| glioblastoma | divisions_per_year | 0 | 0.025 |
| glioblastoma | divisions_per_lifetime | 0 | 2 |
| medulloblastoma | divisions_per_year | 0 | 0.025 |
| medulloblastoma | divisions_per_lifetime | 0 | 2 |
| ovarian_germ | divisions_per_year | 0 | 0.025 |
| ovarian_germ | divisions_per_lifetime | 0 | 2 |

Annual/lifetime fields exceeding the declared arithmetic check tolerance are preserved:

| Tissue | Annual divisions × 80 | Published lifetime divisions |
|---|---:|---:|
| duodenum | 1920 | 1947 |
| duodenum_fap | 1920 | 1947 |
| liver | 73 | 88 |
| liver_hcv | 73 | 88 |
| osteosarcoma | 5.36 | 5 |
| osteosarcoma_arms | 5.36 | 5 |
| osteosarcoma_head | 5.36 | 5 |
| osteosarcoma_legs | 5.36 | 5 |
| osteosarcoma_pelvis | 5.36 | 5 |
| small_intestine | 2880 | 2920 |

### Later international risk-only endpoint audit

The permissively reusable later table supplies **17** cancer categories.
These are independently ranked registry risks through age 74, not matched countries,
not full lifetime risks, and not a joint risk-versus-divisions replication.

| Cancer category | Registry minimum | Registry median | Registry maximum | Registry entries |
|---|---:|---:|---:|---:|
| head & neck squamous cell carcinoma | 0.002169643 | 0.012133814 | 0.029206721 | 423 |
| esophageal squamous cell carcinoma | 0 | 0.001570111 | 0.15536378 | 423 |
| colorectal adenocarcinoma | 0.001297671 | 0.030898802 | 0.067111297 | 423 |
| hepatocellular carcinoma | 0 | 0.00247951 | 0.028237789 | 423 |
| pancreatic cancer (ductal + endocrine) | 0.000215202 | 0.007753628 | 0.014975921 | 423 |
| lung adenocarcinoma | 0 | 0.01195417 | 0.026455077 | 423 |
| osteosarcoma | 0 | 0.000205915 | 0.001824406 | 423 |
| melanoma | 0 | 0.006615544 | 0.049457625 | 423 |
| ovarian germ cell | 0 | 0.000249753 | 0.001249202 | 423 |
| testicular germ cell | 0.000100616 | 0.003502034 | 0.010162599 | 423 |
| medulloblastoma | 0 | 9.48e-05 | 0.000432123 | 423 |
| thyroid follicular/papillary | 0 | 0.005705464 | 0.040835403 | 423 |
| thyroid medullary | 0 | 0.000141167 | 0.000913568 | 423 |
| chronic lymphocytic leukemia | 0 | 0.002330291 | 0.005817165 | 423 |
| acute myeloid leukemia | 0 | 0.002000186 | 0.007249971 | 423 |
| prostate | 0.001025839 | 0.09924712 | 0.24345741 | 423 |
| breast | 0.007625557 | 0.086177005 | 0.11981572 | 423 |

### Normalized content hashes

| File | SHA-256 after newline normalization |
|---|---|
| AMENDMENTS.md | `c0e1f1edc451de114ffc6fa4f8ffc69983bd6d3e8afca6f28c6e9479edc2cdc8` |
| PROTOCOL.md | `c60a2f2d99106a5e65378f9df1337a43edac79efcb0af1afe8f8f9ec0a1cde4b` |
| data/little2016.csv | `baee121cc5c08144177d2add71e777821d549f7ee542a7262d1744f1ba80a32d` |
| data/liu2022_ranked_risk_audit.csv | `919a0e48e5455f38cb8d99b85b40f050a591888dd2061f38eb3e0eb36f4550fe` |
| data/tv2015.csv | `40101d378d61c098b5a7e4610861319114644de7f3b25931b1ec3d8380598704` |
| data/tv2015_nunney2020.csv | `d72da38db78c182aec6659f01b928f00f577e2ca6a2ddb262953a5a67af8e3c6` |

<!-- END RECOMPUTED RESULTS -->


## 5. Verdicts, mechanisms and falsifiable predictions

| Interpretation | Verdict and exact scope |
|---|---|
| More total stem-cell divisions are associated with greater reported risk across these selected categories. | **Supported descriptively.** The historical association is reproduced. Neither causation nor prevention follows from this regression. |
| The observed half-slope counts rate-limiting transitions. | **Not identified.** Theorems 1–2 separate lineage-time powers, directional derivatives, tissue covariance and finite-risk saturation. Theorem 3 constructs every k with tissue-specific rates. |
| A common-rate, common-k, equal-allocation chain fits the point table. | **Falsified as an exact point specification for each evaluated k.** The exact certificate requires a nonunit allowance; even the prespecified stress box does not rescue the tested modest allowances. This is not a sampling significance test. |
| The homogeneous rare-event family gives good cross-tissue prediction. | **Fails the reported held-out comparison.** The descriptive regression predicts better. High-k exact-chain approximation errors are substantial, so poor limiting-model RMSE alone does not reject the finite-rate family. Its exact point-compatibility failure is a separate result. |
| Common division-scaled positive expansion in the specified family resolves the discrepancy. | **Not supported on the preregistered grid.** The selected grid point is no growth. This does not exclude an off-grid fit, multiple growth stages, loss, niches, heterogeneous rates or observation effects. |
| Common chronological expansion at a common endpoint can be measured from these aggregate fits. | **Exactly confounded in the proved limit.** It changes a common amplitude. Theorem 5 gives a broader finite-rate single-endpoint product equivalence when initiation is free. |
| Residuals or R² identify the extrinsic fraction, intrinsic fraction or preventable fraction. | **Provably confounded in the stated class.** Theorem 7 gives identical baseline observations with different cause splits and counterfactuals. Its counterfactual is model-defined, not a clinical prediction. |
| Choosing the k with the smallest rho identifies k. | **Unsupported.** Each k changes the nonlinear scale on which rate heterogeneity is measured. A minimum across these coordinate systems is not evidence of a driver count. |
| Later source corrections or the developmental proxy remove the need for tissue variation. | **Not supported by the reported sensitivities.** Numerical thresholds change, and those changes are retained. The source revision does not identify the correct lineage history. |

### Outlying tissues and the minimal mechanistic explanation

The prespecified group-held-out diagnostic names basal cell carcinoma,
papillary/follicular thyroid carcinoma and gallbladder cancer as the strongest
positive deviations, and small-intestine adenocarcinoma, pancreatic endocrine
carcinoma and duodenal adenocarcinoma as the strongest negative deviations.
These are residuals against the specified descriptive model, not universal
rankings of tumor biology. Some reproduce known findings in Noble et al. [L4].

The smallest additional *model component* is independent, age-independent
Bernoulli thinning of initiated clones. With underlying clone law `q_0(a)`
and probabilities `p_i,p_ref` in [0,1], write `q_i(a)=p_i q_0(a)` and
`q_ref(a)=p_ref q_0(a)` at every clone age a. With the same initiation-time
weights, Poisson marking gives `H_i/H_ref=p_i/p_ref` by Theorem 5.
These are relative success ratios, so positive residual ratios above one are
legitimate without calling them probabilities. At the level of this free-scale
parameterization, choose `p_ref<=min(1,1/max_i ratio_i)` and
`p_i=ratio_i p_ref`; interpreting those quantities biologically also requires
an independently substantiated underlying clone law and initiation supply.
The descriptive regression does not supply an observed causal reference: using
its tabulated ratio as a relative success factor is a conditional benchmark that requires
independent validation of the reference hazard. This is a minimal
parameterization, not identification of the cause. Different target-cell counts, step numbers,
mutation targets, niche replacement, extinction, exposure and diagnosis could
produce the same fitted factors. Normal-tissue clonality measurements need not
measure this success probability [L22].

**Falsifiable mechanism hypothesis.** Suppose the intestinal deficit is caused
solely by reduced persistence/progression of initiated clones, with the
division estimates, initiating mutation intensities, observation horizon and
comparison baseline independently held fixed. Under the additional
age-independent-thinning hypothesis, a matched lineage-tracing experiment
must show the same relative reduction in success at each clone age, equal to
the negative amplitude factor in the table after accounting for initiation
and validating the reference. An age-varying ratio falsifies that thinning
hypothesis. A test of the more general Theorem 5 factor B instead must match
its initiation-time/clone-age weighting; synchronized q(T) alone is not B.
If independently measured success is near the reference while the risk deficit
and the other assumptions persist, that sole-mechanism hypothesis fails. A
measured difference in initiation instead also falsifies its “solely” clause.
These data have not performed such an experiment; the hypothesis is explicitly
unvalidated. The table supplies the effect to be challenged, not evidence that
the mechanism is true.

**Direct aggregate falsifier.** With common k, common intensity and identical
N,D, the model predicts identical risks. The source subgroup contrasts already
violate that equality. Their intensity-ratio bounds can be tested again with
independently measured inputs and age-matched strata. Re-estimating N,D in each
stratum could defeat the identical-input premise; it would not convert the
original contrast into evidence that a particular parameter was identified.

**Additional discriminating observation.** Repeated age-specific incidence or
lineage-division measurements could constrain what the endpoint loses. Under
the exact clean-birth model, small-clock powers depend on k; PR #64 explains
why observation lags, finite bins and fast transitions still matter. A shared
ratio that persists over several independently measured clocks would be a
stronger restriction than the single endpoint. No such data are supplied here.

## 6. Hostile self-review and retained falsification attempts

| Attempt to defeat the conclusion | Outcome or correction retained |
|---|---|
| Search for prior solutions before proposing novelty. | Broad slope/confounding claims were already solved. The claim was narrowed to the coupled exact certificate and its reproducible application. Generic half-range geometry is credited to prior work and PR #64. |
| Explain the slope by probability saturation alone. | The preregistered one-hit fit remains too steep. The separate ordered-cloud lower bound cannot be applied because the observed N/v order condition fails; that failed premise is reported rather than ignored. |
| Remove source-defined high-risk categories or group related observations. | All-row, finite-risk and primary panels and both prediction schemes are retained. The primary panel does not become a causal intrinsic-risk sample. |
| Change the developmental convention and later maintenance-division fields. | Literal versions and the prospective counting-correction sensitivity are reported. A balanced-development proxy changes thresholds materially. None establishes the actual lineage clock. |
| Let a common intermediate expand. | The tested grid selects its no-growth boundary. The positive-growth class is a rigorous limit, not expected malignant counts substituted for event probabilities. |
| Use exact finite-rate inversion instead of a rare-event approximation. | Independent high-precision inversion reproduces the key certificates. Large high-k limiting-fit discrepancies prevent an unsupported rejection of exact models from those prediction scores. |
| Treat N and D perturbations as independent clocks. | Rejected. Theorem 3 preserves v=D/N; Corollary 3.1 shows the opposite N effect under a source-consistent formula. |
| Interpret a printed probability of one as a finite cumulative hazard. | Rejected. It is retained in historical log-risk replication and excluded from finite-risk inversion without clipping. |
| Apply an interval lower-risk substitution to both sides of a rate-ratio inequality. | Rejected on merge review: the lower-rate bound uses (H_b^-/H_a^+)^(1/k), but the upper-rate bound needs H_b^+/H_a^-; the original wording was ambiguous. |
| Treat a shared intrinsic fraction as a different value at every stage while keeping the equal-rate counterfactual formula. | Rejected: Theorem 7 now states the within-tissue equal-rate and equal-split assumptions explicitly. |
| Treat residual ratios as absolute clone probabilities or synchronized clone survival as the endpoint factor. | Rejected. Relative thinning probabilities and the same clone-age weighting are required. The mechanistic prediction is conditional on an independently validated reference. |
| Use a convenient online transcription or infer missing annual/lifetime values. | Rejected. Printed table checks caught esophageal/thyroid transcription pitfalls and preserved source inconsistencies. |
| Combine the later international sorted risk columns with old division inputs. | Rejected because endpoints, input versions and joint registry identities do not match. Coverage remains explicitly incomplete. |
| Infer an environmental fraction from the successful fit of free rates. | Rejected by the exact constructive equivalences. Saturated fitting is a nonidentifiability certificate, not predictive success. |

The mathematical review checked every implication, interval direction,
endpoint, positivity restriction, common-versus-tissue-specific parameter,
and use of an asymptotic limit. The synthetic verifier includes exact rational
polynomial identities and independent high-precision direct-survival inverses.
A separate review recomputed selected empirical frontiers with that direct
formula. Finite numerical checks corroborate implementations; the general
claims rest on the written proofs. Agreement between AI contributors is not
independent biological evidence or external peer review.

No unsupported biological mechanism is called proved. The novelty hypothesis
remains unestablished. If a prior source contains the same composed certificate,
the contribution becomes an independently reproducible synthesis without
changing the mathematical or numerical statements. The work does not classify
all finite-rate branching models, all lineage architectures, all observation
mechanisms or all possible uncertainty regions.

## 7. How to rerun and what is checked

From the repository root, with Python 3.10 or later and no added dependencies:

```bash
python research/stem-cell-divisions-risk/verifier.py
python research/stem-cell-divisions-risk/analysis.py
python -m unittest discover -s tests -p 'test_research_stem_cell_divisions_risk.py' -v
```

Regenerate the machine-readable analysis and a standalone version of the
computed README block:

```bash
python research/stem-cell-divisions-risk/analysis.py \
  --output research/stem-cell-divisions-risk/analysis_results.json \
  --markdown /tmp/stem_cell_divisions_results.md
```

The tests recompute the complete published summary tree from the CSVs, compare the rendered
empirical block, check source metadata and source-version arithmetic, exercise
independent regression/cross-validation identities, test training-only growth
selection, and check the theorem helpers against independent high-precision
calculations. All tests are offline and standard-library only. They use
tolerances for platform-level floating-point variation; reported display
precision is deliberately coarser. Numerical optimization and comparisons
occur before formatting, and source hashes normalize CRLF and lone CR to LF.
The hashes provide local content integrity, not evidence that source estimates
are biologically correct or licenses valid.

`analyze()` also returns every fold prediction and inferred tissue rate for
the independent checks. The checked-in JSON omits those duplicated intermediate
arrays while retaining every displayed number, all model scores, coefficients,
frontier witnesses, outliers and source audits. The omission changes no fit or
reported result; full intermediate values are reproduced by calling `analyze()`.

Only this new directory and
`tests/test_research_stem_cell_divisions_risk.py` are added. Existing research,
policies, workflows, catalog, results registry and root documentation are
unchanged. The existing CI discovers the new test. Publication-review CI checks
repository policy; it does not confer scientific review or clinical validity.

## 8. Primary references

Access-date label for this targeted review: 2026-11-17. Source-specific data
rights and exact locators are recorded separately in the CSVs. Summaries below
do not reproduce publisher tables or figures.

**[L1]** Tomasetti C, Vogelstein B (2015). *Variation in cancer risk among tissues can be explained by the number of stem cell divisions.* Science. [doi:10.1126/science.1260825](https://doi.org/10.1126/science.1260825).

**[L2]** Tomasetti C, Li L, Vogelstein B (2017). *Stem cell divisions, somatic mutations, cancer etiology, and cancer prevention.* Science. [doi:10.1126/science.aaf9011](https://doi.org/10.1126/science.aaf9011).

**[L3]** Wu S, Powers S, Zhu W, Hannun YA (2016). *Substantial contribution of extrinsic risk factors to cancer development.* Nature. [doi:10.1038/nature16166](https://doi.org/10.1038/nature16166).

**[L4]** Noble R, Kaltz O, Hochberg ME (2015). *Peto's paradox and human cancers.* Philosophical Transactions B. [doi:10.1098/rstb.2015.0104](https://doi.org/10.1098/rstb.2015.0104), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4581036/).

**[L5]** Noble R, Kaltz O, Nunney L, Hochberg ME (2016). *Overestimating the role of environment in cancers.* Cancer Prevention Research. [doi:10.1158/1940-6207.CAPR-16-0126](https://doi.org/10.1158/1940-6207.CAPR-16-0126).

**[L6]** Little MP, Hendry JH, Puskin JS (2016). *Lack of correlation between stem-cell proliferation and radiation- or smoking-associated cancer risk.* PLOS ONE. [doi:10.1371/journal.pone.0150335](https://doi.org/10.1371/journal.pone.0150335), [S2 data](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0150335.s002&type=supplementary).

**[L7]** Tomasetti C, Vogelstein B (2017). *On the slope of the regression between stem cell divisions and cancer risk, and the lack of correlation between stem cell divisions and environmental factors-associated cancer risk.* PLOS ONE. [doi:10.1371/journal.pone.0175535](https://doi.org/10.1371/journal.pone.0175535).

**[L8]** Tomasetti C et al. (2017). *Role of stem-cell divisions in cancer risk.* Nature. [doi:10.1038/nature23302](https://doi.org/10.1038/nature23302), [author-hosted text](https://sites.math.duke.edu/~rtd/reprints/paper209.pdf).

**[L9]** Wu S, Zhu W, Hannun YA (2017). *Wu et al. reply.* Nature. [doi:10.1038/nature23303](https://doi.org/10.1038/nature23303).

**[L10]** Nowak MA, Waclaw B (2017). *Genes, environment, and “bad luck”.* Science. [doi:10.1126/science.aam9746](https://doi.org/10.1126/science.aam9746), [primary magazine mathematical discussion](https://www.sciencemagazinedigital.org/sciencemagazine/24_march_2017?pg=31).

**[L11]** Nunney L, Thai K (2020). *Determining cancer risk: the evolutionary multistage model or total stem cell divisions?* Proceedings B. [doi:10.1098/rspb.2020.2291](https://doi.org/10.1098/rspb.2020.2291), [section 5, final paragraph](https://pmc.ncbi.nlm.nih.gov/articles/PMC7779499/).

**[L12]** Altenberg L (2015). *Statistical problems in a paper on variation in cancer risk among tissues, and new discoveries.* Preprint. [arXiv:1501.04605](https://arxiv.org/abs/1501.04605).

**[L13]** Weinberg CR, Zaykin D (2015). *Is bad luck the main cause of cancer?* JNCI. [doi:10.1093/jnci/djv125](https://doi.org/10.1093/jnci/djv125).

**[L14]** Sornette D, Favre M (2015). *Debunking mathematically the logical fallacy that cancer risk is just “bad luck”.* EPJ Nonlinear Biomedical Physics. [doi:10.1140/epjnbp/s40366-015-0026-0](https://doi.org/10.1140/epjnbp/s40366-015-0026-0).

**[L15]** Little MP, Hendry JH (2017). *Mathematical models of tissue stem and transit target cell divisions and the risk of radiation- or smoking-associated cancer.* PLOS Computational Biology. [doi:10.1371/journal.pcbi.1005391](https://doi.org/10.1371/journal.pcbi.1005391).

**[L16]** Ashley DJB (1969). *The two “hit” and multiple “hit” theories of carcinogenesis.* British Journal of Cancer. [doi:10.1038/bjc.1969.41](https://doi.org/10.1038/bjc.1969.41).

**[L17]** Lopez-Garcia C, Klein AM, Simons BD, Winton DJ (2010). *Intestinal stem cell replacement follows a pattern of neutral drift.* Science. [doi:10.1126/science.1196236](https://doi.org/10.1126/science.1196236).

**[L18]** Frank SA, Iwasa Y, Nowak MA (2003). *Patterns of cell division and the risk of cancer.* Genetics. [doi:10.1093/genetics/163.4.1527](https://doi.org/10.1093/genetics/163.4.1527).

**[L19]** McHale PT, Lander AD (2014). *The protective role of symmetric stem cell division on the accumulation of heritable damage.* PLOS Computational Biology. [doi:10.1371/journal.pcbi.1003802](https://doi.org/10.1371/journal.pcbi.1003802).

**[L20]** Rozhok AI, Wahl GM, DeGregori J (2015). *A critical examination of the “bad luck” explanation of cancer risk.* Cancer Prevention Research. [doi:10.1158/1940-6207.CAPR-15-0229](https://doi.org/10.1158/1940-6207.CAPR-15-0229), [primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4849411/).

**[L21]** Liu et al. (2022). *Quantifying substantial carcinogenesis of genetic and environmental factors from measurement error in the number of stem cell divisions.* BMC Cancer. [doi:10.1186/s12885-022-10219-w](https://doi.org/10.1186/s12885-022-10219-w).

**[L22]** Zhang et al. (2024). *Variation in cancer risk between organs can not be explained by the degree of somatic clonal expansion.* [doi:10.1007/s44307-024-00025-9](https://doi.org/10.1007/s44307-024-00025-9).

**[L23]** Blohmer et al. (2025). *Quantifying cell divisions along evolutionary lineages in cancer.* Nature Genetics. [doi:10.1038/s41588-025-02078-5](https://doi.org/10.1038/s41588-025-02078-5).

**[L24]** Nowak MA, Michor F, Komarova NL, Iwasa Y (2004). *Evolutionary dynamics of tumor suppressor gene inactivation.* PNAS. [doi:10.1073/pnas.0400747101](https://doi.org/10.1073/pnas.0400747101). Distinguishes kinetic regimes; a transition count does not by itself specify the number of effective bottlenecks.

**[L25]** Armitage P, Doll R (1954). *The age distribution of cancer and a multi-stage theory of carcinogenesis.* British Journal of Cancer. [doi:10.1038/bjc.1954.1](https://doi.org/10.1038/bjc.1954.1).

**[L26]** Zhang X, Fröhlich H, Grigoriev D, Vakulenko S, Zimmermann J, Weber AG (2018). *A simple 3-parameter model for cancer incidences.* Scientific Reports. [doi:10.1038/s41598-018-21734-x](https://doi.org/10.1038/s41598-018-21734-x), Methods equation 5 and “Relationship between cancer incidence and stem cell divisions.” Its exact independent-target model cites Calabrese and Shibata (2010) as prior work.
