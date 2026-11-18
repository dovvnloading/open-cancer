# How fragile is the ecDNA fair-segregation rejection?

**Nonclinical, unreviewed computational research. No treatment advice.**
This is an aggregate-only robustness audit of an existing observation-based
argument. It is not a new segregation estimate, biological experiment,
clinical finding, or validation of a treatment model.

## Finding and explicit verdict

**LA-N-5 supplies strong conditional evidence under independent, calibrated
mother-level tests, but substantially weaker evidence under a general
exchangeable dependence model. The published counts do not establish the
calibration and dependence restrictions needed for an unconditional rejection
of fair segregation.** The other two panels do not independently rescue that
claim. The exact surviving region, including simultaneous deterioration of
test size and dependence, is computed below.

The primary result is a **sharp sensitivity frontier**: every stated worst-case
tail bound is attained by an explicitly certified count distribution in the
declared model. Thus crossing this frontier really removes this particular
Holm rejection within that model class. It does not prove that fair segregation
is true, that the study actually has the adversarial dependence, or that a
different test could not detect departure.

The narrower beta-binomial indicator model gives a different frontier,
including a reappearance of rejection near complete dependence. We retain
that result as a falsification of monotone-dependence reasoning. It must not
be interpreted as strong dependence making the evidence robust: a cap on
possible dependence allows every smaller dependence as well.

**These correlations are between different mothers' binary test-rejection
indicators. They are not ecDNA allocation dispersion, sister copy-count
correlation, or physical clustering.** Neither actual null-test size nor
indicator correlation is estimated here. The aggregate observations cannot
provide the independent null calibration and culture replication needed to
estimate those quantities credibly, and they cannot supply a new segregation
parameter or eradication threshold.

## Frozen scope, provenance, and relation to previous work

This executes ranked job 3 in [the research programme](../README.md), extending
the independent and dependence-free endpoints in
[#65](../ecdna-segregation-clustering-data/README.md) and the explicitly
post-extraction E3-X check in
[#63](../prediction-confrontation/README.md). The conditional mathematics in
[#55](../ecdna-segregation-therapy/README.md) and
[#60](../ecdna-killing-laws-and-clustering/README.md) is not refuted by a weaker
observation bridge. Those models still require their stated allocation,
copy-conservation, birth/death, and sampling assumptions.

**What this does and does not say about PR #65:** #65's reported
conditional decrease in the ecDNA-lineage persistence/eradication threshold
came from transporting **Shibata et al.'s published fitted segregation
parameters** through PR #60's specified asymmetric beta-binomial killing
model. That numerical shift was **not derived from Lundberg's mother-level
test-rejection counts**. Weakening the Lundberg rejection argument therefore
does not recalculate or invalidate the Shibata-fit model transport; equally,
those point-fit threshold calculations cannot rescue an unconditional
Lundberg rejection or supply missing test calibration, fit uncertainty,
physical allocation validation, or a clinical conclusion. #65's threshold
shift remains conditional on the Shibata fit and model assumptions, not a
new estimate established by the current three rejection counts.

The [protocol](PROTOCOL.md) was frozen locally in
`acc2dc41f8a6673c17be28a4f80170505cf1fb84`, then preserved unchanged publicly in
[`d58107d841841dfb42577e10f1f6065c11c909bc`](https://github.com/dovvnloading/open-cancer/commit/d58107d841841dfb42577e10f1f6065c11c909bc),
**before executing any new probability, frontier, or design calculation**.
The branch started from current main
[`9ded69c3586188b173e81fedb7a79abbf4e7e91f`](https://github.com/dovvnloading/open-cancer/commit/9ded69c3586188b173e81fedb7a79abbf4e7e91f).
The public commit contains only the protocol. Shell push authentication was
unavailable, so publication uses the connected GitHub Git-data interface;
the local and public freeze contain identical protocol text. This is a
prospective computational specification with already-known counts and
previous endpoints, not blinded or external preregistration. No protocol
content was changed after the freeze.

**Source used:** Lundberg et al. (2008), *PLoS ONE*,
[DOI 10.1371/journal.pone.0003099](https://doi.org/10.1371/journal.pone.0003099),
Results “Statistical distribution of DMs at mitosis,” Fig. 3E–G. This run
reads only the existing
[lundberg.csv](../ecdna-segregation-clustering-data/data/lundberg.csv)
transcription, with its recorded source-access date **2026-11-17**; it does
not claim a fresh inspection of the paper. No browser downloads, source
figures, raw observations, or new data copies were used.

The seven selected IDs are `lundberg_chp212_n`, `lundberg_chp212_pct`,
`lundberg_lan5_n`, `lundberg_lan5_pct`, `lundberg_mcix_n`,
`lundberg_mcix_pct`, and `lundberg_alpha`. Whole-percentage rounding is
handled by enumerating compatible integer numerators, conservatively retaining
the smallest if ambiguous. The source, quantity, population and untreated
condition are checked, so irradiation and colony-inferred comparisons cannot
silently enter the primary family.

The CSV records one anaphase/telophase division with opposing poles per mother,
excluded detached double minutes, unspecified per-mother tail construction
and calibration, and unavailable culture allocation. Its recorded reuse basis
is unchanged small numerical facts with attribution under the article's
Creative Commons Attribution statement, whose version is unspecified. The
repository MIT license does not relicense third-party material. See
[DATA_POLICY.md](../../DATA_POLICY.md). No patient-level, controlled-access,
identifying, or reconstructed individual-cell records are involved.

## Model, family, and meaning of rejection

For each line, let `X_i` indicate that the source procedure reports a
mother-level rejection, and let `K=sum X_i`. Fair segregation itself would
give a daughter count `J | N ~ Binomial(N,1/2)` at a correctly observed
conserved total. That biological statement does **not** automatically specify
the law of the source's reported rejection indicators: counting error,
selection, test construction and common culture effects also matter.

Write `alpha'=p=Pr(X_i=1)` under the fair-segregation-plus-observation null.
The nominal source cutoff is not a measurement of `p`. The primary parameters
are caps `p<=A` and `Corr(X_i,X_j)<=R` for finitely exchangeable mothers within
each line. Negative feasible correlations are permitted below this
nonnegative upper cap. At deterministic `p=0,1`, correlation is undefined;
use the deterministic count laws and limiting moment constraints. The beta
comparison instead specifies one exact nonnegative indicator correlation `r`
and one particular mixing distribution. Neither parameter is #60's physical
allocation dispersion or its effective second moment.

The three untreated lines form the entire primary Holm family, at familywise
level `delta=0.05`. Equality rejects. **Joint inconsistency means at least one
Holm rejection of the three-line intersection null**, equivalently
`3 min(p_i)<=delta`; it does not mean all three lines reject. We also state the
all-rejection condition. Arbitrary cross-line dependence is allowed. There is
no pooled fourth test and no cross-line independence assumption. This targeted
retrospective family is not the different combined family in #63, nor a
correction for every analysis in the paper or repository. The familywise
theorem describes a fixed testing procedure; freezing this sensitivity audit
does not retroactively preregister the original study or E3-X.

## Recomputed results and exact robustness verdict

<!-- BEGIN RECOMPUTED RESULTS -->
### Recomputed source counts and endpoints

All displayed probabilities are rounded; all decisions use exact rational arithmetic. Brackets enclose arithmetic boundary locations, not sampling uncertainty.

| Cell line | Mothers n | Compatible rejected k | Independent tail | Independent Holm | Sharp R=0 tail | Sharp R=0 Holm | Markov tail | Markov Holm |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CHP-212 | 100 | 2 | 0.264238021077 | 0.528476042154 | 0.495000000000 | 0.990000000000 | 0.500000000000 | 1.000000000000 |
| LA-N-5 | 100 | 10 | 0.000000076316 | 0.000000228948 | 0.011000000000 | 0.033000000000 | 0.100000000000 | 0.300000000000 |
| MC-IX | 50 | 1 | 0.394993932862 | 0.528476042154 | 0.500000000000 | 0.990000000000 | 0.500000000000 | 1.000000000000 |

At nominal actual size A = 1/100, the sharp exchangeable-cap family rejects **LA-N-5 only** for **0 <= R <= 17/3267 (approximately 0.005203550658096)**. For larger caps it rejects no line. Equality is included.

| Model / dependence assumption | CHP-212 size frontier | LA-N-5 size frontier | MC-IX size frontier | Joint any-rejection size frontier |
| --- | --- | --- | --- | --- |
| Independent indicators | [0.002431336823, 0.002431336824] | [0.045870556724, 0.045870556725] | [0.001025339867, 0.001025339868] | [0.045870556724, 0.045870556725] |
| Arbitrary dependence (Markov) | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] |

The dependence-free joint size frontier is exactly **1/600**; its displayed bracket is intentionally conservative. These frontiers use the complete three-line Holm family, including both weak panels.

### Complete frontier: certified slices of the exact formulas

Entries are [lower, upper] rational decimal brackets for the largest size cap A still rejected. Their width is 0.000000000001. The equations in Theorems 1–5 define every point between these slices; no linear interpolation is asserted. The all-rejection frontier is the minimum of the three line frontiers; the joint any-rejection frontier is their maximum.

#### Sharp exchangeable cap: rho_I <= R

| R or r | CHP-212 A frontier | LA-N-5 A frontier | MC-IX A frontier | Joint any A frontier |
| ---: | --- | --- | --- | --- |
| 0 | [0.002247332874, 0.002247332875] | [0.014707779495, 0.014707779496] | [0.001000000000, 0.001000000001] | [0.014707779495, 0.014707779496] |
| 1/10000 | [0.002197996478, 0.002197996479] | [0.014549442257, 0.014549442258] | [0.001000000000, 0.001000000001] | [0.014549442257, 0.014549442258] |
| 1/2000 | [0.002011642855, 0.002011642856] | [0.013892150425, 0.013892150426] | [0.001000000000, 0.001000000001] | [0.013892150425, 0.013892150426] |
| 1/1000 | [0.001802988453, 0.001802988454] | [0.012998106644, 0.012998106645] | [0.001000000000, 0.001000000001] | [0.012998106644, 0.012998106645] |
| 1/500 | [0.001460644634, 0.001460644635] | [0.011360147953, 0.011360147954] | [0.001000000000, 0.001000000001] | [0.011360147953, 0.011360147954] |
| 1/200 | [0.000862175156, 0.000862175157] | [0.010080669736, 0.010080669737] | [0.000862175156, 0.000862175157] | [0.010080669736, 0.010080669737] |
| 1/100 | [0.000500000000, 0.000500000001] | [0.008311875308, 0.008311875309] | [0.000500000000, 0.000500000001] | [0.008311875308, 0.008311875309] |
| 1/50 | [0.000500000000, 0.000500000001] | [0.005881022639, 0.005881022640] | [0.000500000000, 0.000500000001] | [0.005881022639, 0.005881022640] |
| 1/20 | [0.000500000000, 0.000500000001] | [0.002873427864, 0.002873427865] | [0.000500000000, 0.000500000001] | [0.002873427864, 0.002873427865] |
| 1/10 | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] |
| 1/4 | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] |
| 1/2 | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] |
| 3/4 | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] |
| 9/10 | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] |
| 99/100 | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] |
| 1 | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] | [0.000500000000, 0.000500000001] | [0.001666666666, 0.001666666667] |

#### Pointwise beta-binomial: rho_I fixed at r

| R or r | CHP-212 A frontier | LA-N-5 A frontier | MC-IX A frontier | Joint any A frontier |
| ---: | --- | --- | --- | --- |
| 0 | [0.002431336823, 0.002431336824] | [0.045870556724, 0.045870556725] | [0.001025339867, 0.001025339868] | [0.045870556724, 0.045870556725] |
| 1/10000 | [0.002393900425, 0.002393900426] | [0.045648469331, 0.045648469332] | [0.001027850069, 0.001027850070] | [0.045648469331, 0.045648469332] |
| 1/2000 | [0.002254442400, 0.002254442401] | [0.044776784160, 0.044776784161] | [0.001037853746, 0.001037853747] | [0.044776784160, 0.044776784161] |
| 1/1000 | [0.002101831521, 0.002101831522] | [0.043722968722, 0.043722968723] | [0.001050276866, 0.001050276867] | [0.043722968722, 0.043722968723] |
| 1/500 | [0.001858771440, 0.001858771441] | [0.041725047661, 0.041725047662] | [0.001074863886, 0.001074863887] | [0.041725047661, 0.041725047662] |
| 1/200 | [0.001451651584, 0.001451651585] | [0.036468486189, 0.036468486190] | [0.001146755705, 0.001146755706] | [0.036468486189, 0.036468486190] |
| 1/100 | [0.001215205854, 0.001215205855] | [0.029550418436, 0.029550418437] | [0.001215205854, 0.001215205855] | [0.029550418436, 0.029550418437] |
| 1/50 | [0.001133488064, 0.001133488065] | [0.020210192641, 0.020210192642] | [0.001133488064, 0.001133488065] | [0.020210192641, 0.020210192642] |
| 1/20 | [0.001307980249, 0.001307980250] | [0.009261772327, 0.009261772328] | [0.001307980249, 0.001307980250] | [0.009261772327, 0.009261772328] |
| 1/10 | [0.001728469974, 0.001728469975] | [0.005832073137, 0.005832073138] | [0.001728469974, 0.001728469975] | [0.005832073137, 0.005832073138] |
| 1/4 | [0.003088124325, 0.003088124326] | [0.005273012177, 0.005273012178] | [0.003088124325, 0.003088124326] | [0.005273012177, 0.005273012178] |
| 1/2 | [0.006034844534, 0.006034844535] | [0.007094533299, 0.007094533300] | [0.006034844534, 0.006034844535] | [0.007094533299, 0.007094533300] |
| 3/4 | [0.010365999186, 0.010365999187] | [0.010365999186, 0.010365999187] | [0.010365999186, 0.010365999187] | [0.010365999186, 0.010365999187] |
| 9/10 | [0.013554096558, 0.013554096559] | [0.013554096558, 0.013554096559] | [0.013554096558, 0.013554096559] | [0.013554096558, 0.013554096559] |
| 99/100 | [0.016304097299, 0.016304097300] | [0.016304097299, 0.016304097300] | [0.016304097299, 0.016304097300] | [0.016304097299, 0.016304097300] |
| 1 | [0.016666666666, 0.016666666667] | [0.016666666666, 0.016666666667] | [0.016666666666, 0.016666666667] | [0.016666666666, 0.016666666667] |

#### Joint trade-offs at nominal or worse size

These points use the closed-form primary frontier proved below. They are simultaneous size/dependence trade-offs, not separate tolerances that may be combined arbitrarily.

| Size cap A | Largest correlation cap R retaining rejection | Decimal R |
| ---: | ---: | ---: |
| 1/100 | 17/3267 | 0.005203550658096 |
| 23/2000 | 2543/1500543 | 0.001694719844750 |
| 3/250 | 109/73359 | 0.001485843591107 |
| 7/500 | 149/341649 | 0.000436120111577 |

### Complete nominal-size beta-binomial crossing audit

The verifier isolates every root of every raw tail minus 0.05/j, j in {1,2,3}. Exact Bernstein sign certificates exclude missing roots; the original rational beta tails certify every bracket endpoint. Only changes in the full Holm rejection set are reported as family transitions.

| Line | Holm divisor j | Polynomial degree | Roots on [0,1] |
| --- | ---: | ---: | ---: |
| CHP-212 | 1 | 99 | 1 |
| CHP-212 | 2 | 99 | 1 |
| CHP-212 | 3 | 99 | 1 |
| LA-N-5 | 1 | 99 | 0 |
| LA-N-5 | 2 | 99 | 2 |
| LA-N-5 | 3 | 99 | 2 |
| MC-IX | 1 | 49 | 1 |
| MC-IX | 2 | 49 | 1 |
| MC-IX | 3 | 49 | 1 |

There are **10 raw crossings** and **2 actual Holm transitions**.

| Correlation root bracket | Rejected immediately below | Rejected at the algebraic boundary | Rejected immediately above |
| --- | --- | --- | --- |
| [0.045806289979, 0.045806289980] | LA-N-5 | LA-N-5 | none |
| [0.728174162473, 0.728174162474] | none | CHP-212, LA-N-5, MC-IX | CHP-212, LA-N-5, MC-IX |

At r=0 the rejection set is **LA-N-5**; at r=1 it is **CHP-212, LA-N-5, MC-IX** and every adjusted value is exactly **3/100**. Thus the apparent reappearance of rejection is real within the restricted beta model. It does not occur for an upper dependence cap and does not make highly dependent data robust evidence.

### Attaining adversarial laws for LA-N-5

Each row is a hypothetical null count distribution, not reconstructed observations. Choose a uniformly random subset of K mothers to receive indicator 1. The laws have mean K=1 and actual per-test size 1/100.

| Cap / witness | Pr(K=0) | Pr(K=1) | Pr(K=10) | Tail | Holm value |
| --- | ---: | ---: | ---: | ---: | ---: |
| R=0 | 99/1000 | 89/100 | 11/1000 | 11/1000 | 33/1000 |
| R=17/3267 | 3/20 | 5/6 | 1/60 | 1/60 | 1/20 |

The R=0 witness has pairwise independent Bernoulli indicators but is not a mutually independent vector. It disproves replacing an uncorrelated-indicator assumption with a binomial count law.

### Prospective design arithmetic (not an experiment performed)

The following critical counts use n=100 mothers per line, A=1/100, and the conservative first Holm step 1/60. The previous count fails that step. Caps here are assumed known; the next example budgets calibration uncertainty separately.

| Indicator assumption | R cap | Minimum rejected mothers | Tail at critical count | Tail one count below |
| --- | ---: | ---: | ---: | ---: |
| independent | 0 | 5 | 0.003432321588 | 0.018374036445 |
| exchangeable | 0 | 9 | 0.013750000000 | 0.017678571429 |
| exchangeable | 1/200 | 10 | 0.016445000000 | 0.020556250000 |
| exchangeable | 1/100 | 12 | 0.014925000000 | 0.017910000000 |
| markov | 1 | 60 | 0.016666666667 | 0.016949152542 |

For an illustrative simultaneous calibration failure budget 1/100 and test-family level 1/25, their sum is 1/20. With A=1/100 and R=1/200, the first-step threshold is 1/75 and the minimum count is **12/100** (tail 0.011212500000; one count below 0.013455000000). This error budget is an illustrative design choice made after the freeze, not achieved calibration or a power calculation.

### Newline-normalized SHA-256 manifest

The source, immutable protocol and original implementation/test files are hashed after UTF-8 CRLF/CR-to-LF normalization. Path keys come from Path.as_posix(). README is excluded to avoid self-reference.

| Path | SHA-256 |
| --- | --- |
| research/ecdna-segregation-clustering-data/data/lundberg.csv | `d6d89d7beb6757db0a93174674c3e0baa38e363f98c5835f66007ad4b917ebea` |
| research/ecdna-rejection-robustness/PROTOCOL.md | `bed446347c21936e5561d712b03cc361b026fc8fb9e96b02d2df27684b414185` |
| research/ecdna-rejection-robustness/analysis.py | `0ceed68520a52f63f573aed5a0bb6aac572f08a6e9ebdacd428dc072ff446ca6` |
| research/ecdna-rejection-robustness/verifier.py | `d5a94717e25f60544eb42f40e90bf2bb6803cd75465c41e71425cd6e9aa81d49` |
| tests/test_research_ecdna_rejection_robustness.py | `ccefefecca5757ee0aec637e7be98109f0cdf20613a9327194b1c099c55ece58` |
<!-- END RECOMPUTED RESULTS -->

## Mathematical guarantees

The numbered theorems prove the bounds and decision rule used above. Exact
finite checks challenge their implementations; they are not substitutes for
these all-parameter proofs. The source's actual satisfaction of the
assumptions remains unverified.

### Theorem 1 — count laws, exact tails, and increasing size

**Assumptions.** Fix integers `n >= 1`, `1 <= k <= n`, and a marginal
indicator probability `p in [0,1]`. Under independence let
`K = sum_{i=1}^n X_i` with iid Bernoulli(p) indicators. Alternatively, for
`0 < r < 1`, let `c=(1-r)/r`, let `Theta ~ Beta(pc,(1-p)c)`, and draw the
indicators conditionally iid Bernoulli(Theta). Treat `p=0,1` by their
deterministic limits. Then the respective upper tails are

$$
T_{\rm ind}(n,k;p)=\sum_{s=k}^{n}\binom ns p^s(1-p)^{n-s},
$$

$$
T_{\rm BB}(n,k;p,r)=\sum_{s=k}^{n}\binom ns
 \frac{(pc)^{\overline s}((1-p)c)^{\overline{n-s}}}
 {c^{\overline n}}.
$$

Here `x^{overline j}` is the rising factorial, with empty product one.
The beta-binomial indicators have mean `p`, and pairwise correlation `r`
when `n>=2` and `0<p<1`.
The limits are `T_BB(p,0)=T_ind(p)` and `T_BB(p,1)=p`.
For every fixed dependence parameter each tail is continuous and
nondecreasing in `p`. Consequently a size cap `p <= A` is handled exactly
by evaluating the stated model at `p=A`.

**Proof.** Independence gives the binomial mass by counting subsets.
Under the mixture, integrate `Theta^s(1-Theta)^{n-s}` against the beta
density; the beta-integral ratio is the displayed rising-factorial ratio.
Conditional expectation and covariance give
`E X_i=p` and `Cov(X_i,X_j)=Var(Theta)=p(1-p)/(c+1)`, hence correlation
`1/(c+1)=r` for interior `p`. At `r downarrow 0`, divide every factor
by `c` to recover the binomial mass. At `r uparrow 1`, the endpoint
masses tend to `1-p,p`; every interior mass has two vanishing numerator
factors and only one vanishing denominator factor, so tends to zero.

For independence, use the same independent uniforms and set
`X_i(p)=1{U_i <= p}`. For interior `r`, the ratio of beta densities
at `p_2 > p_1` is a positive constant times
`(theta/(1-theta))^{c(p_2-p_1)}`, a strictly increasing function.
Since both densities integrate to one, their difference changes sign
once, from negative to positive; integration proves that the larger-p
CDF is everywhere no larger. An inverse-CDF coupling therefore makes
the larger-p mixing variable at least as large, and the same conditional
uniforms make every indicator nondecreasing. The endpoint laws have the
same order. Finite sums or their stated endpoint limits give continuity.
All claimed comparisons follow. This argument makes no monotonicity
claim in `r`. ∎

### Theorem 2 — sharp exchangeable tail at size and correlation caps

**Assumptions.** Fix `n >= 2`, `1 <= k <= n`, and caps `0 <= A,R <= 1`.
Allow every finitely exchangeable Bernoulli vector whose common mean
`p <= A` and, when `0 < p < 1`, pairwise indicator correlation is at most
`R`. Negative feasible correlation is allowed. At `p=0,1` use the
corresponding deterministic law; correlation itself is undefined there.
Write

$$
\mu=nA,\qquad B=n(n-1)\{A^2+RA(1-A)\}.
$$

If `mu >= k`, the sharp tail supremum is one. If `mu < k`, it is

$$
U(n,k;A,R)=\min\left\{\frac\mu k,
\min_{0\le j\le k-2}
 \frac{B-2j\mu+j(j+1)}{(k-j)(k-j-1)}\right\},
$$

with the inner minimum omitted for `k=1`. Every displayed supremum is
attained by an exchangeable indicator law. `U` is continuous and
nondecreasing in each cap. No independence conclusion follows from
`R=0`.

**Proof.** Exchangeability identifies a vector law with its count law:
conditional on `K=s`, all s-subsets are equally likely. Conversely any
probability law on `{0,...,n}` produces an exchangeable vector by this
construction. Its moments satisfy
`E K=np` and `E K(K-1)=n(n-1) E X_1X_2`. Thus the correlation constraint
is precisely the factorial-moment inequality in the protocol.

First reduce the unknown mean to `A`, rather than assuming this reduction.
For `p<A<1`, augment each zero independently with probability
`t=(A-p)/(1-p)`, leaving every one unchanged. The resulting vector `Y`
is exchangeable, has marginal mean `A`, and has count no smaller than
the original count. Conditional independence of augmentation gives
`Cov(Y_i,Y_j)=(1-t)^2 Cov(X_i,X_j)`. For interior `p,A`, its correlation is

$$
\operatorname{Corr}(Y_i,Y_j)
=\operatorname{Corr}(X_i,X_j)\frac{p(1-A)}{A(1-p)}.
$$

The multiplier is in `[0,1]`. Nonnegative correlation decreases, while
negative correlation remains nonpositive; therefore the upper cap `R>=0`
is preserved. At `p=0` use independent Bernoulli(A) augmentation, and
at `A=1` use all ones. Hence the supremum under `p<=A` equals that at
`p=A`.

At fixed `mu`, every integer `s>=0` satisfies
`1{s>=k} <= s/k`. For `j=0,...,k-2`, it also satisfies

$$
\mathbf1_{\{s\ge k\}}\le
Q_j(s):=\frac{(s-j)(s-j-1)}{(k-j)(k-j-1)}.
$$

The numerator is nonnegative on integers and is at least its value at
`s=k` for `s>=k`. Taking expectations and using the positive coefficient
of `E K(K-1)` yields all stated upper bounds. The constant one supplies
the trivial bound.

It remains to prove sharpness. Feasibility is never an issue: a
Binomial(n,A) law has mean `mu` and factorial second moment
`n(n-1)A^2 <= B`. In particular, expectations of every nonnegative
integer polynomial under that feasible law certify the lower moment
inequalities used below.

Suppose `mu<k`. If the Markov bound is the minimum and `k>=2`, comparison
with `Q_0` gives `B >= (k-1)mu`. Put mass `mu/k` at `k` and the remaining
mass at zero. Its factorial moment is `(k-1)mu <= B` and its tail equals
the bound. This same construction works for `k=1`, with factorial moment
zero.

Otherwise choose any minimizing quadratic index `j`. Put
`d=k-j`, `C_j=(B-2j mu+j(j+1))/(d(d-1))`, and assign masses

$$
q_k=C_j,\qquad q_{j+1}=\mu-j-dC_j,\qquad
q_j=j+1-\mu+(d-1)C_j.
$$

They sum to one, have mean `mu`, factorial moment `B`, and tail `C_j`.
Their nonnegativity is essential. `C_j>=0` follows from the nonnegative
integer polynomial `Q_j` and feasibility, while
`C_j<=mu/k<1` follows from minimality. For `j>=1`, direct subtraction
gives

$$
C_{j-1}-C_j=\frac{2q_{j+1}}{d(d+1)},
$$

so minimality proves `q_{j+1}>=0`; for `j=0`, use `C_0<=mu/k` instead.
For `j<=k-3`, subtraction similarly gives

$$
C_{j+1}-C_j=\frac{2q_j}{(d-1)(d-2)}.
$$

At the final index `j=k-2`, where that neighbor is absent, use

$$
d q_j=B-(2k-2)\mu+k(k-1)\ge0.
$$

The inequality follows by taking the expectation of the nonnegative
integer polynomial `(S-(k-1))(S-k)` under a feasible count law and then
replacing its factorial moment by the larger cap `B`. Thus every mass
is valid, including branch ties.

If `mu>=k`, use the law on `floor(mu),ceil(mu)` with mean `mu` (a point
mass when `mu` is an integer). Its tail is one. This law minimizes the
factorial second moment among integer-valued laws of that mean: writing
`m=floor(mu)`, the inequality `(S-m)(S-m-1)>=0` proves the lower bound,
and the adjacent-support law attains it. Its moment is therefore no
larger than that of the feasible binomial law and hence no larger than
`B`. This covers the deterministic endpoints as well.

Every count witness gives the required exchangeable vector by uniform
subset allocation. Since the feasible classes expand with each cap,
their sharp supremum is nondecreasing. For continuity, each candidate
formula is continuous; at `mu=k`, feasibility and nonnegative variance
give `B>=k(k-1)`, so every quadratic candidate is at least one, while
the Markov candidate equals one. Thus the pieces
join continuously. ∎

### Theorem 3 — dependence-free envelope and its sharpness

**Assumptions.** Let `X_i` be arbitrary Bernoulli variables, without any
exchangeability or independence requirement, with `E X_i<=A` for every
mother. Then

$$
\Pr\{K\ge k\}\le M(n,k;A):=\min(1,nA/k).
$$

This is sharp even within the exchangeable subclass. Consequently
`U(n,k;A,1)=M(n,k;A)`.

**Proof.** On `K>=k`, `K>=k 1{K>=k}`; taking expectations gives the
bound. When `nA<k`, put mass `nA/k` at `K=k` and otherwise at zero,
using uniformly chosen subsets. Every marginal is exactly A and the
bound is attained. When `nA>=k`, the fixed count `K=k` gives marginal
`k/n<=A` and tail one. Correlation at most one is no restriction; the
same laws establish equality with Theorem 2 at `R=1`. ∎

### Theorem 4 — valid tail p-values and Holm familywise control

**Assumptions.** For each cell line fix its declared null class and let
`u_i(k)` be a nonincreasing function of the observed count that bounds
`Pr(K_i>=k)` for every law in that class. The exact tails and sharp or
Markov envelopes above have this property. With `p_i=u_i(K_i)`, each
`p_i` is a conservative p-value. For the three prespecified lines,
Holm's step-down rule controls the probability of any false rejection
at delta under arbitrary dependence between lines.

**Proof.** For `0<=t<=1`, the set of counts with `u_i(k)<=t` is either
empty or an upper tail starting at an integer `k_t`. Therefore
`Pr(p_i<=t)=Pr(K_i>=k_t)<=u_i(k_t)<=t`. This holds separately for every
null law and also for the supremum envelopes.

Order the three p-values increasingly, breaking ties deterministically,
and reject successively while `p_(j)<=delta/(4-j)`. Suppose there are
`m_0` true line nulls and some true null is rejected. At the first true
rejection, no true null has yet been removed, so the number of remaining
hypotheses is at least `m_0`. That rejection implies that at least one
true-null p-value is at most `delta/m_0`. The union bound and the
individual p-value inequalities give probability at most
`m_0(delta/m_0)=delta`. No cross-line independence is used. The case of
no true null is vacuous.

For a line i its Holm adjusted value equals

$$
h_i=\min\left(1,
\max_{S\ni i}|S|\min_{j\in S}p_j\right).
$$

Indeed a subset whose least p-value has ordered rank r has at most
`4-r` members. Its contribution is bounded by `(4-r)p_(r)`, and the
entire suffix beginning at r attains that bound whenever it contains i.
This recovers the familiar cumulative maxima over ranks up to i. It
also proves that every adjusted value is coordinatewise nondecreasing
in the raw p-values, regardless of ordering changes. Ties give the same
adjustment to all tied values.

At least one rejection is equivalent to `3 min_i p_i<=delta`. This is
the stated intersection-null test. Rejection of all three additionally
requires the second and third Holm steps. It is not a pooled count test. ∎

### Theorem 5 — complete frontier defined without grid interpolation

**Assumptions.** Use one of the declared tail models, a fixed dependence
parameter (or cap), shared size cap A, and the three-line family above.
Then for each line the exact rejection region in A is a closed interval
`[0,A_i^*]`, where

$$
A_i^*(r)=\sup\{A\in[0,1]:h_i(A,r)\le\delta\}.
$$

The boundary satisfies `h_i(A_i^*,r)=delta`. The any-rejection boundary
is `max_i A_i^*(r)`; the all-rejection boundary is `min_i A_i^*(r)`.
Together with the exact tails and subset formula in Theorem 4, these
equations specify the whole two-dimensional frontier, not only the
reported fixed-dependence slices. Under the sharp-cap model the
boundaries are nonincreasing in the dependence cap. No such direction
is asserted for the pointwise beta-binomial dependence parameter.

**Proof.** Theorems 1–3 give continuous nondecreasing raw p-values in
A; Theorem 4 preserves these properties under Holm adjustment.
At A=0 every observed positive count has p-value zero; at A=1 all tails
and adjusted values equal one. The sublevel set at `0<delta<1` is
therefore a nonempty closed initial interval with an interior endpoint.
Continuity forces equality at its endpoint. Taking the union or
intersection of these three initial intervals gives the any/all rules.
Set inclusion in the dependence cap proves the sharp-cap direction.

For rational dependence slices, binomial tails are polynomials in A;
beta-binomial tails are polynomials in A divided by a positive,
A-independent rising factorial; the sharp envelope is a finite minimum
of rational-coefficient quadratics and a linear function. Holm uses
only finite minima and maxima. The boundaries therefore have exact
algebraic definitions. Rational brackets with exact endpoint decisions
certify their numerical locations. To claim all crossings as dependence
varies in the pointwise beta model requires a complete root-count
certificate in addition; a grid alone supplies none. ∎

### Specialization of Theorem 5: an explicit joint trade-off

For nominal or worse calibration `A>=1/100`, neither CHP-212 nor MC-IX can
reject under the sharp cap model: their raw bounds at the smallest such size
and zero correlation cap already exceed `delta`, and the bounds are
nondecreasing in both caps. Joint rejection therefore holds exactly when
the LA-N-5 raw bound is at most `q=delta/3=1/60`.

Let

$$
A_* = \frac{25+\sqrt{130}}{2475}.
$$

For `1/100<=A<=A_*`, the complete primary frontier is

$$
R\le R_*(A)=
\begin{cases}
\displaystyle\frac{3/2-9900A^2}{9900A(1-A)},
 &1/100\le A\le23/2000,\\[6pt]
\displaystyle\frac{200A-4/5-9900A^2}{9900A(1-A)},
 &23/2000\le A\le A_*.
\end{cases}
$$

For `A>A_*`, no nonnegative correlation cap retains the primary rejection.
At nominal calibration the first branch gives `17/3267`; at zero correlation
cap the second branch gives `A_*`. These are **joint trade-offs**: one cannot
combine the largest separate size and correlation tolerances into a rectangle.
All smaller sizes, including each weak line's own frontier, remain covered
by the complete formulas and certified tables above.

**Derivation.** For LA-N-5 put `mu=100A` and
`B=9900[A^2+RA(1-A)]`. In this range the Markov candidate exceeds `q`.
The sharp tail is at most `q` exactly when

$$
B\le\max_{0\le j\le8}g_j(\mu),\qquad
g_j(\mu)=q(10-j)(9-j)+2j\mu-j(j+1).
$$

Successive differences are
`g_(j+1)-g_j=2[mu-(j+1)-q(9-j)]`, decreasing with `j`.
For `1<=mu<=2`, only indices zero and one can maximize: their tie is
`mu=23/20`, while the next possible tie is `mu=32/15>2`.
Thus `max g_j` is `3/2` followed by `200A-4/5`. Solving the two linear
inequalities in `R` gives the displayed formula. At `R=0` the larger root of
`9900A^2-200A+4/5=0` is `A_*`, which lies below `1/50`.
The same explicit branches give strict non-rejection throughout
`A_*<A<=1/50`; cap monotonicity extends this beyond `1/50`. This last step
does not assume that monotonicity alone rules out an equality plateau.

### Theorem 6 — complete root certificates for the nominal beta frontier

**Assumptions.** Fix rational interior size `p`, a line's integer `n,k`, and
a rational threshold `t`. Cancel the common factor `1-r` from the exact
beta-binomial probability formula before clearing denominators. The resulting
polynomial `f(r)` satisfies

$$
T_{\rm BB}(n,k;p,r)-t
=s\frac{f(r)}{D(r)},\qquad s>0,\qquad
D(r)=\prod_{j=0}^{n-2}(1+jr)>0\quad(0\le r\le1).
$$

Thus its roots and signs are exactly those of the tail comparison, including
the complete-dependence endpoint. Assume `f` is not identically zero;
the zero polynomial is rejected rather than given a finite root certificate.
A finite subdivision of the unit interval
whose nonempty leaves have either zero Bernstein coefficient sign changes,
or one sign change and nonzero opposite-sign endpoint values, supplies a
complete root certificate: zero-change leaves contain no open-interval root;
each one-change leaf contains exactly one simple root. Exact subdivision
endpoint roots are separately recorded. If the leaves do not resolve, this
procedure makes no completeness claim.

**Proof.** In the rising-factorial formula substitute `c=(1-r)/r` and cancel
powers of `r`. The denominator contains one factor `1-r`; every outcome's
numerator has at least one such factor. Its symbolic cancellation leaves
the stated positive product. The explicit endpoints agree with Theorem 1,
so cancellation introduces neither a missing endpoint nor a spurious root.

On any subinterval, rescale to `x in [0,1]` and write its degree-d polynomial
as `sum_i b_i C(d,i) x^i(1-x)^(d-i)`. If the nonzero coefficients all have
one sign (with at least one nonzero coefficient), every Bernstein basis term
is positive in the interior, and their
sum cannot vanish. If there is one sign change, reverse the overall sign
if necessary so that the coefficients have a negative block followed by a
positive block. With `u=x/(1-x)>0`, multiply by `(1+u)^d`. The transformed
polynomial is `sum_i b_i C(d,i)u^i`, with the same coefficient signs.
Let m be the first positive coefficient's index and divide by `u^m`.
Every negative term now has negative exponent and positive derivative;
every positive nonconstant term has positive exponent and positive
derivative. There is at least one strictly positive derivative. The function
therefore increases strictly from a negative limit to a positive limit,
and has exactly one simple zero. This proves the only sign-count cases
actually accepted; no unproved general variation bound is needed.

Exact midpoint de Casteljau subdivision is repeated affine evaluation of
the same polynomial and preserves these facts. The finite leaves partition
the interval; inspecting shared endpoints separately and counting each once
proves completeness. The verifier either completes this certificate or raises
an error at its explicit work limit. For the reported polynomials all roots
are isolated, simple, and in pairwise disjoint brackets.

For Holm, include every equation `p_i(r)=delta/j`, with each line i and
`j=1,2,3`. Between consecutive certified roots all of these comparison signs
are constant, so the closed-testing formula in Theorem 4 makes every Holm
decision constant too, even if raw-p order changes. At a single isolated
root, substitute its exact threshold for the crossing raw value; all other
threshold comparisons remain fixed. This decides equality membership without
pretending a rational bracket midpoint is the algebraic root. The finite
certificates therefore establish every nominal-size family transition. ∎

### Theorem 7 — calibration uncertainty consumes an error budget

**Assumptions.** Separate calibration data produce simultaneous caps for all
tested lines that contain their true null size and dependence parameters
with probability at least `1-epsilon`. Conditional on every calibration
output for which the caps are valid, the test sample satisfies the declared
count null class. Independent calibration/test samples under a stable
measurement pipeline are one sufficient way to justify this conditional
assumption. Apply Holm with familywise level `delta_test` using those caps.
Then the probability of any false line rejection is at most
`epsilon+delta_test`.

**Proof.** Let C denote the event that all calibration caps are valid. On
each such calibration outcome Theorem 4 bounds the conditional false-rejection
probability by `delta_test`. Average over those outcomes and add
`Pr(C^c)<=epsilon` for the remaining outcomes. The union bound gives the
claim. If the calibration/test pipeline changes, or conditioning on the
calibration sample changes the test law beyond the specified class, the
conditional assumption fails and this proof supplies no guarantee. ∎

## What would make a new experiment decisive?

The critical counts in the report specify a concrete prospective decision
rule, not the number of biological cultures needed to estimate correlation
or achieve desired power. No alternative distribution or achieved assay
precision is available for a defensible power claim.

1. **Fix the observational null before scoring.** Measure maternal totals
   and both daughters, check conserved copies and loss, and predeclare the
   exact mother-level test including ties, discreteness, zero handling,
   detection thresholds and exclusions. Randomize and blind scoring. A
   binomial simulation checks test code; it does not calibrate real microscopy
   error, missed double minutes, or selection of mitoses.
2. **Make independence or a small dependence cap credible.** Preselect one
   scored division per independently initiated culture/lineage where feasible,
   and record culture, batch and ancestry. If several mothers share a culture,
   retain that block structure and enough independent cultures for a declared
   dependence analysis. Do not turn sisters, repeated images or related
   mothers into independent replicates. Shared batch effects still need
   controls even when cultures are separately initiated.
3. **Calibrate separately and budget uncertainty.** Use appropriate independent
   null controls for the complete pipeline and establish simultaneous upper
   bounds on actual test size and any retained indicator dependence. Numerical
   calibration of an idealized null alone does not establish the physical
   null-control assumptions. Choose the family/error budget before observing
   biological rejection counts; Theorem 7 explains the cost of uncertain caps.
4. **Release adequate aggregates.** Publish per-culture rejection numerators
   and denominators, conditional maternal-total/paired-daughter count tables,
   zero/loss outcomes, exclusion counts and calibration summaries. This can
   preserve the repository's aggregate-only boundary. Independently replicated
   line/culture assays can test reproducibility without manufacturing raw
   observations from the old percentages.

A repeated LA-N-5 excess that crosses the chosen critical count with those
conditions met would decisively reject the declared fair-allocation-plus-
pipeline null at the chosen error rate. It would still require further
measurements to distinguish biased allocation, correlated copies, loss,
variable replication, or measurement failure, and would not validate the
killing-law models. The present aggregate record cannot retrospectively
establish those experimental conditions.

## Falsification attempts and hostile self-review

| Attempted assertion or shortcut | Outcome retained |
| --- | --- |
| The tiny independent tail is an assumption-free biological rejection. | Rejected. The sharp correlation-cap and Markov results explicitly expose calibration and sampling assumptions. |
| Zero pairwise indicator correlation implies a binomial count. | Refuted by the exact attaining law in the report. It has dependent higher-order structure with pairwise independent Bernoulli indicators. |
| Larger positive correlation always weakens a beta-binomial rejection. | Refuted by complete exact root isolation and the loss/re-entry phases. The primary cap frontier remains monotone because its class expands. |
| A beta-binomial tail is worst-case among exchangeable laws with the same moment restrictions. | Refuted by exact primal witnesses and comparison with the sharp envelope. No minimum across incompatible model classes is used as a robust p-value. |
| A moment inequality alone proves that crossing it removes inconsistency. | Rejected. Every sharp value has a feasible attaining law and a dual polynomial majorant. Independent small-count LP enumeration challenges both. |
| Holm can use a fixed line ordering along the frontier. | Rejected. Actual sorting and an independent closed-testing calculation agree, including ordering changes and ties. |
| The weak panels can be discarded or pooled to improve the result. | Rejected. Every scenario retains the frozen three-line family; no pooled additional experiment is claimed. |
| A numerical grid establishes all beta crossings. | Rejected. All raw threshold polynomials are certified; unresolved root isolation raises an error. |
| A rational midpoint is the algebraic root when deciding equality. | Rejected. Equality uses exact threshold substitution and certified signs for the other comparisons. |
| Monotonicity from equality alone proves no equality plateau beyond the closed-form size root. | Rejected during proof review. The active branches prove strict non-rejection through the next rational endpoint, then monotonicity applies. |
| Calibration confidence bounds may be treated as certain without an error budget. | Rejected. The illustrative design separates calibration failure and testing error explicitly. |
| Passing the new tests proves the source, biology, or historical novelty. | Rejected. Tests establish reproducible arithmetic and finite certificates, with the limits stated throughout. |

Hostile review independently rebuilt the probabilities using Bernoulli
convolution/Pólya-urn recursion and rising factorials, enumerated small-count
LP vertices, checked every primal moment and pointwise dual inequality,
challenged invalid numerical inputs and corrupted certificates, and verified
Bernstein transformations and subdivision identities. It also tested exact
endpoint roots and a repeated-root case that must fail openly when the
isolation algorithm cannot finish. The source selection, every displayed
number and every reported boundary endpoint are recomputed by the new test.

Implementation review caught a dynamic-import/dataclass namespace issue and
corrected it before publication. A separate proof check required explicit
handling of deterministic marginal endpoints, branch ties and equality at
irrational crossings, and excluding the identically zero polynomial from a
finite root certificate. No source count, failed comparison, or frozen protocol
was changed to accommodate the results. This is computational self-review by
AI agents, not independent biological validation or external scientific review.

**Post-freeze choices:** the simultaneous calibration/testing budget in the
report is an illustrative experimental-design choice, explicitly labelled
there. It does not alter the primary family, sensitivity domain, count law,
decision rule, or fixed frontier slices. Bernstein sign subdivision supplies
the protocol's required complete polynomial root count; it is an alternative
to the example Sturm algorithm, not a weaker grid. No other analytical
amendment is needed.

## Honest novelty assessment and remaining limits

The contribution is **original-to-this-repository quantification of the
fragility of a particular published rejection-count argument**: a complete
two-parameter decision definition, a closed-form primary trade-off in the
nominal-or-worse region, attained exchangeable envelopes, and a certified
complete nominal beta comparison. It adds no measurements and does not
discover unequal ecDNA inheritance or a new general statistical method.

Binomial and beta-binomial laws, exchangeable count representations,
integer moment inequalities, Markov's inequality, Holm correction and
Bernstein polynomial subdivision are established mathematical tools. Local
proofs and independent implementations do not establish their historical
novelty. No new literature search or bibliographic-priority claim is made;
global novelty is unestablished. The appropriate claim is a narrow,
reproducible sensitivity audit extending #65's two endpoints.

The scope deliberately preserves the prior note's counterevidence:
[#65's opposing summaries](../ecdna-segregation-clustering-data/README.md#binomial-baseline-challenges-and-opposing-summaries)
report nonsignificant Lange comparisons, which do not establish binomiality,
and fitted nonzero dispersion in chromosomal controls, which warns against
identifying a fit parameter with physical clustering. They are not additional
tests or calibration inputs here. Heterogeneous copy totals, assay selection,
unknown testing choices and culture dependence remain unresolved. The shared
size/correlation caps are sensitivity restrictions, not measured properties
of the paper, and exchangeability itself is an assumption.

An arithmetic enclosure is not a confidence interval; a sharp bound is not
an estimate of realized dependence; loss of rejection is not equivalence or
confirmation. The old counts alone cannot say where the actual study lies
on this frontier. They therefore remain a **conditional challenge**, with
unconditional biological rejection and new segregation inference unestablished.

## Reproduce, verify, and review the PR

From the repository root, offline, using only Python's standard library:

```bash
python -B research/ecdna-rejection-robustness/analysis.py --check
python -B research/ecdna-rejection-robustness/analysis.py --markdown
python -B research/ecdna-rejection-robustness/analysis.py --json
python -B research/ecdna-rejection-robustness/verifier.py
python -B -m unittest discover -s tests -p 'test_research_ecdna_rejection_robustness.py' -v
```

`--check` is read-only; it compares the generated report and normalized hashes
with this README. JSON gives exact fractions, root brackets, primal laws and
root-count certificates. The dedicated test checks the full report and its
under-one-minute workload. Newline-normalized hashing and POSIX-form path
keys make the manifest independent of Windows/Linux checkout line endings.
The README omits its own recursive hash; its numerical block is checked by
exact regeneration instead. File hashes establish file identity, not source
truth or statistical validity.

This contribution adds only `README.md`, the immutable `PROTOCOL.md`,
`analysis.py`, `verifier.py` in this directory and
[the new test](../../tests/test_research_ecdna_rejection_robustness.py).
Existing data, research notes, indexes, workflows and tests remain unchanged.
The existing pull-request workflow discovers the new test and runs repository
validation. The PR records observed CI status; successful CI is not a
scientific-review or clinical-validity claim.
