# ecDNA killing laws, segregation, and an identifiability obstruction

**Nonclinical computational mathematics. No treatment advice.** This is an
unreviewed mathematical contribution about a stipulated branching process.
All numerical experiments are synthetic. No patient data, treatment regimen,
drug efficacy claim, or experimentally validated cancer theory is supplied.
Proofs and executable checks establish model statements; AI self-review and
passing CI do not establish biological validity or scientific priority.

## Results and relation to PR #55

We extend [PR #55](https://github.com/dovvnloading/open-cancer/pull/55) and its
[note](../ecdna-segregation-therapy/README.md), which were read together with
their tests before this work. That contribution supplies the replicated-copy
genealogy, complementary sisters, first-event extinction recursion, the exact
exponential-hazard benchmark, and the previous linear-hazard bracket
`1/3 <= gamma*/b <= 1/2` (its lower endpoint was unresolved). Those ingredients
are reused with attribution.

The principal new result closes that gap through an **exact infinite-state
critical eigenproblem and convergent rational certificates**:

$$
\boxed{\gamma_*(b)=b c_*,\qquad
0.378474994554<c_*<0.378474994555.}
$$

For unbiased replicated binomial segregation and death rate
$d_k=\gamma k$, the entire positive lineage survives with positive
probability exactly when $\gamma<\gamma_*$. **Extinction is certain at
equality.** The decimal endpoints above are exact rationals, certified using
25 and 48 positive types respectively **plus an analytic bound on every
omitted type**. We do not claim an elementary closed form for $c_*$.

Other results are:

* A complete power-law classification for binomial segregation: every
  $\alpha>0$ in $d_k=\gamma k^\alpha$ has a positive finite critical
  value, characterized by the same eigenproblem; $\alpha=0$ has threshold
  $\gamma_*=b$; every $\alpha<0$ permits persistence at every finite
  $\gamma$. The last result also holds for the clustered and biased
  families defined below.
* Saturation in **copy number** does not preclude a finite eradication
  threshold. In contrast, a hazard whose high-copy limit remains strictly
  below $b$ cannot eradicate the binomial model. A general high-copy
  refuge proof makes this distinction precise.
* The infinite-state criterion extends to beta-binomial segregation and
  asymmetric allocation. With effective segregation variance $R$,
  $c_*\leq(1-R)/2$ for linear killing. This proves that sufficiently
  strong clustering/asymmetry makes eradication easier and that its threshold
  tends to zero at complete clumping. It does **not** prove monotonicity
  between every two interior values of $R$.
* A sharp identifiability obstruction: two segregation mechanisms can have
  the **same entire two-generation mean and Fano dose curves, yet different
  eradication thresholds**. A fourth factorial moment, or a specified zero
  fraction difference, distinguishes them. Conditional variance formulas
  themselves are established beta-binomial mathematics, not a novelty claim.

All theorems concern arbitrary finite founder populations and unbounded copy
numbers unless a particular assay design is explicitly stated. An empty or
all-zero founder population is already eradicated.

## 1. Literature audit and model choice

Audit date: **17 November 2026**, using the user's local date. The scope is
primary literature and the inspected versions, not an exhaustive priority
search. We link to sources; no third-party data or figures are redistributed.

| Source | What is reused or constrains the claim |
| --- | --- |
| [PR #55](https://github.com/dovvnloading/open-cancer/pull/55) | The underlying replicated-binomial branching model, first-event recursion, complementary sisters, finite-size visitation extinction argument, and exponential-killing benchmark are reused. Its linear gap is the starting problem. Its Fano formula applies to a different killing law. |
| [Turner et al., Nature 2017](https://doi.org/10.1038/nature21356), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5334176/); [Lange et al., Nature Genetics 2022](https://doi.org/10.1038/s41588-022-01177-x) | ecDNA branching, copy-dependent selection, replicated binomial partition, and neutral moments are prior art. |
| [Shibata et al., Nature Communications 17, 9544 (2026)](https://doi.org/10.1038/s41467-026-76566-5), [primary text](https://pmc.ncbi.nlm.nih.gov/articles/PMC13547322/) | Paired daughter measurements and a **beta-binomial model with allocation bias and intraclass correlation** already exist. Some clones change segregation under stress. We do not claim this parameterization or its variance formula as new. Parameters held fixed during a condition are an assumption to test. |
| [Sankar et al., Nature 649, 152–160 (2026)](https://doi.org/10.1038/s41586-025-09764-8), online 19 November 2025 | Retention elements support chromosome tethering and nuclear transmission. Loss of retention can remove copies from daughter nuclei. Copy-conserving asymmetric allocation below does not model that loss and is not calibrated by this paper. |
| [Scanu, Werner & Huang, bioRxiv v3, 26 July 2025](https://doi.org/10.1101/2024.10.22.619675), [inspected text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12330626/) | A preprint on multiple variants, switching, selection, and de-mixing. The inspected Methods use independent fair binomial allocation, so this version is **not evidence for clustered or biased segregation**. |
| [Kimmel & Axelrod, Genetics 1990](https://doi.org/10.1093/genetics/125.3.633); [Harnevo & Agur, Mathematical Biosciences 1991](https://doi.org/10.1016/0025-5564(91)90094-y) | Amplification, copy-dependent cell loss, and branching extinction are historical prior art. Accessible abstracts/summaries were checked; full historical kernel equivalence is unresolved. |
| [Li, Kibby & Foo, Journal of Mathematical Biology 2023](https://doi.org/10.1007/s00285-023-01992-8) | Amplification-driven escape and extinction conditions exist for other transition laws. The general idea of an amplification threshold is not new. |
| [Marguet & Smadi, SPA 168 (2024), 104262](https://doi.org/10.1016/j.spa.2023.104262), [author preprint](https://arxiv.org/abs/2211.08265); [earlier sharing-law work](https://arxiv.org/abs/2010.16070) | Burden-dependent cell death/division and complementary partition in branching populations are an established mathematical area. |

Searches included the named titles/DOIs and combinations of "ecDNA
beta-binomial extinction", "ecDNA linear killing threshold", "gene
amplification copy-dependent death", "amplification branching extinction",
and "parasite sharing death division rates". Primary indexed passages were
used where publisher/PMC opening encountered a sign-in redirect or CAPTCHA;
those barriers were not bypassed.

**Provisional novelty assessment.** The proposed contribution is the specific
replicated-copy killing-law critical characterization, its critical-extinction
proof, certificates controlling the infinite tail, the classification and
scaling results below, and the threshold identifiability obstruction. The
audit did not locate these results for the same kernel. Finite
Perron–Frobenius theory, monotone branching coupling, beta-binomial moments,
and first-event equations are standard tools. We have not ruled out a
historical theorem implying one of these statements after translation.

## 2. Model and observables

Let $b>0$. A cell with $k\geq1$ copies has independent competing
exponential division and whole-cell death clocks of rates $b$ and
$d_k$, with $0\leq d_k<\infty$ at every finite type. Thus its
division-before-death probability is

$$
p_k=\frac{b}{b+d_k}.
\tag{1}
$$

A dividing mother is replaced by daughters with counts $J$ and
$2k-J$. Exactly $k$ copies replicate to $2k$; sisters are
complementary, never independent samples. Different mothers have independent
split randomness. Copy numbers change only at division. There is no density
dependence, immigration, new ecDNA formation, chromosomal reintegration,
variable replication efficiency, or inherited shared environment. Zero-copy
cells cannot regenerate ecDNA. Their fate is irrelevant to positive
eradication; in all pooled-composition formulas take $d_0=0$, so
$p_0=1$. Even when $\alpha\leq0$, $k^\alpha$ is used only for
$k\geq1$, never at zero.

The three segregation families are:

| Family | Ordered split law $K_k(j)=\Pr(J=j\mid k)$ | Effective $R$ for a uniformly chosen sister |
| --- | --- | --- |
| Binomial | $J\sim\mathrm{Bin}(2k,1/2)$ | $0$ |
| Clustered | Draw $P\sim\mathrm{Beta}(a,a)$ afresh per mother, then $J\mid P\sim\mathrm{Bin}(2k,P)$; $a>0$ | $\rho=1/(2a+1)\in(0,1)$ |
| Asymmetric allocation | $J\sim\mathrm{Bin}(2k,p)$, $0<p<1$, with the other sister complementary | $(2p-1)^2$ |

For the last family a uniformly sampled sister uses the symmetric mixing law
$P\in\{p,1-p\}$, each with probability $1/2$. The **ordered**
daughter has mean $2kp$; the pooled sisters have mean $k$. Both
daughters retain the same subsequent fate law depending only on their own
count. Persistent daughter identities with different future birth/death rates
would be a different, multitype model. We call this asymmetric allocation or
copy-conserving retention, not a model of nuclear retention loss.

The clustered probability mass function is

$$
K_k(j)=\binom{2k}{j}\frac{(a)^{\uparrow j}(a)^{\uparrow(2k-j)}}
 {(2a)^{\uparrow 2k}},
\tag{2}
$$

where rising factorials are used. As $\rho\downarrow0$ the law becomes
binomial; at $\rho=1$ use its weak limit supported on $0,2k$.
The bias endpoints $p=0,1$ also put all copies in one sister.

Write $Z_n(k)$ for birth into genealogical generation $n$,
$Z_n^+=\sum_{k\geq1}Z_n(k)$, and
$\pi_n(k)=\mathbb EZ_n(k)/\mathbb E\sum_jZ_n(j)$. A pooled Fano
factor is the variance/mean under $\pi_n$, including zeros. This is a
ratio of expectations from pooled independent replicate populations, not the
expectation of a within-well ratio or a positive-only gated statistic.

**Eradication** means that the entire positive genealogy eventually becomes
empty; **persistence** means positive probability that it never does. Neither
means tumor elimination. Nonexplosion follows by domination of division counts
by the rate-$b$ Yule population $Y_t$: stopping at $M$ cells gives
$\Pr(\sup_{s\leq t}Y_s\geq M)\leq N_0e^{bt}/M$. Deaths add no
cells and their number is bounded by initial cells plus divisions. Thus even
unbounded $d_k$ causes no finite-time explosion. Finite genealogies finish
in finite time almost surely; infinite genealogies cannot finish in finite
time. Calendar-time and genealogical eradication events agree. Composition at
a fixed calendar time is **not** given by the generational formulas below.

## 3. Exact infinite-state criterion and the linear boundary

### Theorem 1 — the compact-tail branching criterion

Suppose $0<p_k<1$, $p_k\to0$, the split has support in
$0,\ldots,2k$, and $K_k(j)>0$ for every such $j$. This includes
all three interior segregation families with $d_k/b=c k^\alpha$,
$c>0,\alpha>0$. On positive types define

$$
M(k,j)=p_k[K_k(j)+K_k(2k-j)],\quad 1\leq j\leq2k,
\qquad M(k,j)=0\ (j>2k).
\tag{3}
$$

Let $M_N$ be the principal matrix on $1,\ldots,N$, corresponding
to discarding positive daughters above $N$, and let $r_N$ be its
Perron root. Then $r_N\uparrow r\in(0,2)$, and

$$
\boxed{\Pr(\text{persistence})>0\quad\Longleftrightarrow\quad r>1.}
\tag{4}
$$

This is true from every finite initial population containing a positive cell.
At $r=1$, eradication is certain. There are strictly positive vectors
$h\in c_0$, $v\in\ell^1$, with $\|h\|_\infty=1$,
$\sum v_k=1$, satisfying $Mh=rh$, $vM=rv$.

**Proof.** Every finite $M_N$ is irreducible: any type can reach type 1,
and type 1 reaches any type through finitely many steps of at most doubling.
All diagonal entries are positive. Finite Perron–Frobenius theory gives
positive left/right vectors, monotonicity under entrywise comparison, and
$0<M(1,1)\leq r_N\leq2\sup_kp_k<2$. The strict upper inequality
uses $p_k\to0$ and $p_k<1$, hence $\sup_kp_k<1$.

Extend finite matrices and vectors by zero. For bounded sequences,

$$
\|M-M_N\|_\infty\leq
2\sup_{k>\lfloor N/2\rfloor}p_k\longrightarrow0.
\tag{5}
$$

Indeed a row with $k\leq N/2$ has no discarded daughters, and every
remaining row has sum at most $2p_k$. The image of the unit ball under
$M$ has tails bounded by $2p_k$, uniformly tending to zero; its
finitely many head coordinates are bounded. Coordinatewise subsequence
selection plus that tail bound gives compactness in the supremum norm.

Normalize right finite Perron vectors $h_N$ to supremum one. Since
$\|Mh_N-r_Nh_N\|\to0$ and $r_N\geq M(1,1)>0$, compactness
yields a convergent subsequence to $h$ of norm one, with $Mh=rh$.
For left vectors $v_N$, normalized to sum one, the support bound gives

$$
r_N\sum_{j>L}v_N(j)
\leq2\sup_{k>\lfloor L/2\rfloor}p_k.
\tag{6}
$$

Consequently they are tight on the positive integers. A diagonal subsequence
converges in $\ell^1$ to a mass-one $v$, and boundedness of $M$
permits passing to $vM=rv$. Irreducibility propagates any zero coordinate
of either eigenvector to every coordinate, so both are strictly positive.
This constructs the required infinite vectors without assuming an
infinite-dimensional Perron theorem.

Let $s_k$ be the probability of positive survival from type $k$,
with $s_0=0$. Conditioning on the first event, using independent future
subtrees **after conditioning on the complementary split**, gives

$$
s=Ms-Q(s),\qquad
Q_k(s)=p_k\sum_{j=0}^{2k}K_k(j)s_js_{2k-j}.
\tag{7}
$$

Survival requires division, so $s_k\leq p_k$ and $s\in c_0$.
If any $s_k>0$, finite reachability makes all positive coordinates
strictly positive. Pairing (7) with $v$, all sums being absolutely
convergent and nonnegative, gives

$$
(r-1)\langle v,s\rangle=\langle v,Q(s)\rangle>0,
\tag{8}
$$

because $Q_1=p_1K_1(1)s_1^2>0$. This proves extinction for $r\leq1$,
including equality.

If $r>1$, some finite $r_N>1$. Let its positive eigenvector be
$h$, and $H=\max h$. In the finite killed-at-exit process put
$h=0$ outside $1,\ldots,N$. For its quadratic term,
$xy\leq H(x+y)/2$ gives
$Q_k(h)\leq Hr_Nh_k/2$. Choose
$0<\varepsilon\leq\min\{1/H,2(r_N-1)/(Hr_N)\}$.
The extinction map sends $1-\varepsilon h$ to at most
$1-\varepsilon h$. Iteration from extinction-by-generation-zero
probabilities therefore has limit at most this vector, proving actual finite
subprocess survival. Embedding it in the full process and finite reachability
prove (4). The same first-event iteration proves minimality of the extinction
fixed point, exactly as in PR #55. Independent founders finish the proof. ∎

### Theorem 2 — exact characterization, certified bounds, and equality

For binomial segregation and $d_k=\gamma k$, there is a unique
$c_*>0$ with

$$
r(c_*)=1,\qquad
(1+c_*k)h_k=2\sum_{j=1}^{2k}\binom{2k}{j}4^{-k}h_j,
\quad h\in c_0,\ h_k>0,\ \|h\|_\infty=1.
\tag{9}
$$

This eigenproblem characterizes the threshold exactly. Persistence holds for
$c=\gamma/b<c_*$, and extinction for $c\geq c_*$. The rational
interval in the opening section is a rigorous enclosure of its unique value.
Finite lower bounds and explicit infinite-tail upper certificates converge to
the same $c_*$, so arbitrary requested precision is possible in principle.

**Proof of existence and uniqueness.** For $d>c>0$,

$$
\frac cdM_c\leq M_d\leq\frac{1+c}{1+d}M_c.
\tag{10}
$$

Apply finite Perron comparison, then take limits. These inequalities imply
continuity and strict decrease of $r(c)$. At $c=0$, the two-type
matrix is
$\begin{pmatrix}1&1/2\\1/2&3/4\end{pmatrix}$, whose root is
$(7+\sqrt{17})/8>1$. It remains supercritical for small positive
$c$. Conversely $r(c)\leq2/(1+c)\to0$. The intermediate value
theorem gives the unique root. Theorem 1 covers critical extinction; at
$c=0$ a positive division always retains at least one positive daughter,
so survival is certain. A strictly positive bounded solution to (9) at any
$c>0$, paired with $v_c$, forces $r(c)=1$, proving the stated
exact characterization rather than just a necessary equation.

**Finite certificates with no biological copy cap.** Let
$t_N(k)=\sum_{j>N}M_c(k,j)$ for $k\leq N$, and solve exactly

$$
(I-M_{c,N})x=t_N.
\tag{11}
$$

If all $x_k<0$, set $w=-x>0$. Then
$M_{c,N}w-w=t_N\geq0$, with a nonzero entry. Pair with the positive
finite left Perron vector to conclude $r_N>1$: this certifies persistence.

If all $x_k>0$, extend $w_k=x_k$ for $k\leq N$ and
$w_k=1$ for $k>N$. The head equations give $M_cw=w$ there.
For every tail row,

$$
(M_cw)(k)\leq
\frac{2}{1+c(N+1)}\max\{1,\max_{j\leq N}x_j\}.
\tag{12}
$$

When the right side is at most one, $M_cw\leq w$ everywhere and
$\inf w>0$. Therefore the expected positive population is uniformly
bounded by initial weight divided by $\inf w$. Each positive cell has
death probability at least $\delta=c/(1+c)>0$. At population sizes
between 1 and $L$, next-generation extinction probability is at least
$\delta^L$. Summing the decreases of survival probabilities gives
$\sum_n\Pr(1\leq Z_n^+\leq L)\leq\delta^{-L}$. Borel–Cantelli
implies $Z_n^+\to\infty$ on survival, contradicting the bounded mean
by Fatou. Thus (12) certifies extinction in the **unbounded** model.
If (12) is strict, pairing the bounded weight with the strictly positive
infinite left vector gives $r(c)<1$, hence a strict upper bound on
$c_*$.

The verifier uses $c_-=189237497277/500000000000$, $N=25$, and
$c_+=75694998911/200000000000$, $N=48$. It checks exact residuals,
strict signs, finite irreducibility, and (12), whose left side is less than
$1$. Rounded displays of weights are not certificates.

**Why the bounds match.** For each $N\geq2$, the unique root
$c_N$ of $r_N(c_N)=1$ increases to $c_*$: every $c<c_*$
has a finite supercritical subsystem, and no $c\geq c_*$ does.
For fixed $c>c_*$, all finite resolvents in (11) are nonnegative
Neumann series, and $x_N>0$ by irreducibility and positive exits.
Moreover $\|t_N\|_\infty\to0$. If $\|x_N\|\not\to0$,
choose a subsequence with norms bounded below and normalize it to one.
Equation (11), (5), and compactness give a nonzero nonnegative limit
$My=y$. Pairing with $v$ contradicts $r(c)<1$. Thus
$x_N\to0$, and (12) eventually holds strictly. Finite strictly
supercritical Perron witnesses can be approximated by rational vectors with
strict positive residual, giving rational lower certificates. Together these
facts provide matching lower and upper certification algorithms; no stability
of an uncorrected truncation is assumed. ∎

## 4. General killing laws and high-copy refuges

### Theorem 3 — power-law classification and limits

For unbiased binomial segregation, $d_k/b=c k^\alpha$ on positive
types has the following classification, from every positive finite founder:

| Exponent | Exact qualitative statement |
| --- | --- |
| $\alpha>0$ | A unique $0<c_\alpha<\infty$; persistence for $c<c_\alpha$, extinction for $c\geq c_\alpha$. Replace $1+ck$ in (9) by $1+ck^\alpha$. |
| $\alpha=0$ | $c_0=1$, with extinction at equality. |
| $\alpha<0$ | Persistence has positive probability for every finite $c\geq0$: no finite eradication threshold. |

For $\alpha>0$, $c_\alpha$ is nonincreasing in $\alpha$, and

$$
0<c_\alpha\leq2^{-\min(\alpha,1)},\qquad
\lim_{\alpha\downarrow0}c_\alpha=1,\qquad
\boxed{c_\alpha\sim\tfrac12\,2^{-\alpha/2}\quad(\alpha\to\infty).}
\tag{13}
$$

**Proof, increasing laws.** For every $\alpha>0$, $p_k\to0$.
Theorem 1 applies, the two-type block stays supercritical for sufficiently
small $c$, and comparison (10) holds with $k$ replaced by
$k^\alpha\geq1$. This gives a unique positive root and equality
extinction. The upper bound is the $R=0$ case of Theorem 5 below.

Increasing $\alpha$ decreases every division probability. Couple
potential binary-tree nodes using the same uniform fate variable and the same
split whenever both cells divide; the higher-death positive tree is a subtree
of the other. This proves monotonicity. For any fixed $c>0$, put
$w_1=1$, $w_k=\varepsilon$ for $k\geq2$, with
$0<\varepsilon\leq\min(1,2c)$. Then
$(M w)_1=(1+\varepsilon/2)/(1+c)\leq1$, and for $k\geq2$,
$(M w)_k\leq2/(1+c2^\alpha)\leq\varepsilon$ for sufficiently
large $\alpha$. The bounded-weight extinction argument above proves
$c_\alpha\to0$. In fact the following exact two-type calculation
proves the sharper scaling in (13), rather than merely a zero limit.

Put $T=2^\alpha$. The cap-two matrix has rows
$(p_1,p_1/2)$, $(p_2/2,3p_2/4)$, where
$p_1=(1+c)^{-1}$, $p_2=(1+cT)^{-1}$. Its exit vector is
$(0,5p_2/8)$. Equation (11) gives exactly

$$
\Delta=T c^2+c/4-1/4,\qquad
x_1=\frac{5}{16\Delta},\qquad x_2=2cx_1.
\tag{13a}
$$

Take $c=C2^{-\alpha/2}$, for any fixed $C>0$. If $C<1/2$,
then $\Delta\to C^2-1/4<0$; both weights are negative and certify
persistence for all sufficiently large $\alpha$. If $C>1/2$,
both are positive, $x_1$ stays bounded and $x_2\to0$. The exact
tail bound is $2/(1+c3^\alpha)$, tending to zero because
$3/\sqrt2>1$; hence (12) certifies eradication for all sufficiently
large $\alpha$. Squeezing between every $C<1/2$ and every
$C>1/2$ proves (13). This scaling is a theorem for all real exponents
in the limit, not an extrapolation from integer-power simulations.

The remaining claims follow from the next elementary refuge lemma. ∎

### Lemma 3A — a binomial high-copy refuge with a rational certificate

If $\limsup_{k\to\infty}d_k/b<1$, unbiased binomial segregation has
positive persistence probability, even if death at finitely many low types is
arbitrarily large but finite. More explicitly, for integer $L\geq2$,
let $p_{\min}=\min_{L<k<2L}p_k$. A sufficient condition is

$$
2p_{\min}\left(1-\frac5L-\frac{192}{L^3}\right)>1.
\tag{14}
$$

**Proof.** In the finite band $L<k<2L$, use
$h_k=\sin(\pi(k-L)/L)>0$, zero outside. For a full sine on all
integers, the binomial characteristic function gives exactly

$$
\mathbb E\sin(\pi(J-L)/L)
=h_k\cos(\pi/(2L))^{2k}.
$$

The truncated $h$ is at least that sine on reachable types except
$j\geq3L$, where the sine is at most one. Since $k<2L$,
$J\leq4L-2$, so no further positive sine interval can be reached.
For $J\sim\mathrm{Bin}(2k,1/2)$,
$\mathbb Ee^{t(J-k)}=\cosh(t/2)^{2k}\leq e^{kt^2/4}$.
The inequality follows termwise from $(2m)!\geq2^m m!$.
Markov's inequality optimized at $t=2u/k$ gives
$\Pr(J-k\geq u)\leq e^{-u^2/k}$, hence
$\Pr(J\geq3L)\leq e^{-L/2}$.

As $\min h=\sin(\pi/L)\geq2/L$, each band row satisfies

$$
\frac{\mathbb E h_J}{h_k}\geq
\cos(\pi/(2L))^{4L}-\frac L2e^{-L/2}
\geq1-\frac5L-\frac{192}{L^3}.
\tag{15}
$$

The last step uses $\cos x\geq1-x^2/2$, Bernoulli's inequality,
$\pi^2<10$, and $e^{L/2}\geq(L/2)^4/24$. The elementary cosine
bound is positive at the arguments used, so raising it preserves the
inequality. Condition (14) gives a strictly supercritical finite band by
finite Perron comparison, hence actual survival. The right side of (15)
tends to one. A tail division probability strictly above $1/2$
therefore supplies a band satisfying (14). Every finite positive type can
reach that band with positive probability. ∎

**Completion of Theorem 3.** At $\alpha=0$, $c<1$ satisfies the
refuge lemma. For $c\geq1$, dominate all positive offspring by a
zero-or-two process with division probability $1/(1+c)\leq1/2$.
Its extinction probability is one, including the critical case (the least
fixed point of $q=1-p+pq^2$ is one). At $\alpha<0$, the tail
hazard tends to zero for every finite $c$, so Lemma 3A proves the claim.
Finally for any $c<1$ choose a finite supercritical band at
$\alpha=0$. Its entries vary continuously with $\alpha$, so it
stays supercritical for sufficiently small positive $\alpha$.
Together with the upper bound in (13), this proves $c_\alpha\to1$
as $\alpha\downarrow0$. ∎

### Theorem 4 — saturation and declining hazards

1. For any segregation law and $d_k/b=c g(k)$ with
   $\inf_{k\geq1}g(k)=m>0$, every $c\geq1/m$ eradicates,
   including equality. Thus copy-number saturation with positive low-copy
   killing admits a finite eradication parameter.
2. For binomial segregation and $g(k)=k/(H+k)$, $H>0$, every
   $c<1$ persists and every $c\geq H+1$ eradicates. If
   $H>2$, persistence is possible even at $c=1$, and indeed for
   every $c$ satisfying
   $$
   8c^2+(H+2)c-(H+1)(H+2)<0.
   \tag{16}
   $$
3. A binomial hazard family with
   $\limsup_k d_k(D)/b<1$ at **every finite dose parameter** $D$
   never eradicates at a finite $D$. For example,
   $$
   d_k(D)=\Gamma_{\max}\frac{D}{D_0+D}\frac{k}{H+k},
   \quad 0\leq\Gamma_{\max}\leq b,\quad 0\leq D<\infty,\quad D_0,H>0,
   \tag{17}
   $$
   has this property. Here saturation in dose and saturation in copy number
   are distinct assumptions.
4. The no-finite-threshold conclusion for $d_k=\gamma k^\alpha$,
   $\alpha<0$, also holds for every beta-binomial $a>0$, every
   asymmetric $0\leq p\leq1$, and complete clumping.

**Proof.** In (1), $p_k\leq1/(1+cm)\leq1/2$; the critical
zero-or-two domination just used proves extinction. In (2), the tail hazard
is $c$, giving persistence for $c<1$ by Lemma 3A, and
$m=1/(H+1)$ gives the upper bound. In the two-type subsystem,
$p_1=(H+1)/(H+1+c)$, $p_2=(H+2)/(H+2+2c)$, and

$$
M_2=\begin{pmatrix}p_1&p_1/2\\p_2/2&3p_2/4\end{pmatrix}.
$$

The numerator of $\det(I-M_2)$, over a positive denominator, is (16).
When it is negative the two real eigenvalues straddle one, so the positive
Perron root exceeds one. At $c=1$, (16) equals $8-2H-H^2<0$
for $H>2$. This disproves a universal saturation threshold of $b$.
Statement (3) is directly Lemma 3A; (17) has tail strictly below $b$
at every finite $D$, even when $\Gamma_{\max}=b$.

For (4), the fair binomial case is already proved. In a division-only
potential genealogy, follow the larger daughter. Its count never decreases.
For beta or nonfair biased allocation there are $\varepsilon,\eta>0$
and a finite $K$ such that, uniformly for $k\geq K$,
$\Pr(\max(J,2k-J)\geq(1+\varepsilon)k)\geq\eta$.
To see this, choose $\delta>0$ with
$\Pr(|P-1/2|\geq\delta)>0$. Conditional binomial variance is at
most $k/2$; Chebyshev shows that, for large $k$, the realized
fraction is within $\delta/2$ of $P$ with probability bounded
away from zero. Use $\varepsilon=\delta$ and a smaller positive
$\eta$. Uniform random variables couple these increases to independent
Bernoulli($\eta$) successes. The strong law then gives exponential
growth of the followed counts almost surely. For $\alpha<0$,
$\sum_n c k_n^\alpha<\infty$. Conditional on that split sequence,
the independent probability of never being killed along this lineage is
$\prod_n(1+c k_n^\alpha)^{-1}>0$, since
$\log(1+x)\leq x$. Its expectation is strictly positive. Reaching
$K$ initially has positive probability. Complete clumping has
$k_n=2^nk_0$ directly, giving the same product proof. ∎

The exact critical value and critical behavior for a general bounded
saturating hazard **are not claimed**. Its division probabilities do not
vanish at infinity, so Theorem 1 cannot be applied unchanged. Statements
(16)–(17) and the finite/never classification have independent proofs.

## 5. How clustering and asymmetric allocation change the threshold

### Theorem 5 — exact criterion, general shift, and a singular limit

For either interior clustered or asymmetric allocation and $\alpha>0$,
Theorem 1 gives the exact criterion for every $c>0$: replace the kernel
in (9) and $k$ by $k^\alpha$. Define
$c_*=\sup\{c>0:r(c)>1\}$, taking the supremum to be zero if that
set is empty. This is finite. If $c_*>0$, it is the unique solution
$r(c_*)=1$, and extinction holds at equality. If $c_*=0$, all
$c>0$ eradicate, but $c=0$ still gives certain survival.

For a symmetric mixing law of $P$, including both families, put
$R=4\operatorname{Var}(P)$, $a_1=(1-R)/2$. Then for
$0\leq R<1$,

$$
c\geq a_1^{\min(\alpha,1)}\quad\Longrightarrow\quad
\text{eradication},\qquad
c_*\leq a_1^{\min(\alpha,1)}.
\tag{18}
$$

In particular linear killing has $c_*\leq(1-R)/2$. With the certified
binomial lower bound $c_-$, **every**
$R>1-2c_-=0.243050010892$ has a strictly smaller threshold than
binomial segregation. This is a rigorous comparison on a whole parameter
range, not an extrapolation from the numerical table. For every fixed
$\alpha>0$, $c_*\to0$ as $R\to1$.

At complete clumping, or bias $p=0,1$, exactly one positive daughter
has count $2k$. For $\alpha\geq0$, every $c>0$ eradicates,
whereas $c=0$ preserves a positive lineage forever.

**Proof.** The exact criterion follows from Theorem 1. Comparison (10)
holds for $k^\alpha\geq1$, so $r(c)$ is continuous and strictly
decreasing. The upper bound below ensures finiteness. This proves the
statements about a positive root, without incorrectly deducing a finite
supercritical core merely from death-free survival.

At division let $A_k$ be the probability that both sisters are positive.
There is always at least one positive sister, so their expected positive
number is $1+A_k$. We claim

$$
A_k=1-\mathbb E[P^{2k}+(1-P)^{2k}]
\leq\min(k a_1,1).
\tag{19}
$$

For fixed $p$, set $q=1-p$, $f_n=1-p^n-q^n$. For
integers $m,n\geq2$,

$$
f_m+f_n-f_{m+n}
=(1-p^m)(1-p^n)+(1-q^m)(1-q^n)-1
\geq2p^2q^2\geq0.
$$

Here replacing both exponents by two decreases each nonnegative product;
$(1-p^2)^2+(1-q^2)^2-1=2p^2q^2$.
Inducting gives $f_{2k}\leq kf_2$, including equality at $k=1$,
and averaging proves (19), since $\mathbb Ef_2=a_1$.

For $\alpha\geq1$, $A_k/k^\alpha\leq a_1$. For
$0<\alpha\leq1$, write $x=ka_1$; then
$\min(x,1)/k^\alpha=a_1^\alpha\min(x,1)/x^\alpha
\leq a_1^\alpha$. Under (18), the expected positive offspring
$(1+A_k)/(1+ck^\alpha)\leq1$. Positive death probability is
uniformly at least $c/(1+c)>0$, so the finite-size visitation and
Fatou argument in Theorem 2 proves extinction, including equality.
The limit follows by squeezing between zero and (18). For complete clumping
the infinite-lineage probability is
$\prod_{n\geq0}[1+c(2^nk_0)^\alpha]^{-1}=0$ for
$\alpha\geq0,c>0$. At $c=0$ no positive cell dies. ∎

**What is and is not proved about interior parameters.** Exact criterion (4)
and bound (18) hold for every permitted parameter, with no type cap.
Strictly positive thresholds for the examples below have finite survival
witnesses. We do not claim a separate proof here that every nondegenerate
beta parameter has a positive threshold, or that $c_*$ is monotone
between every two interior $R$'s. The latter is a **conjecture**, suggested
by exploratory calculations but unnecessary for any theorem. Distinct
mechanisms at the same $R$ need not have the same threshold.

### Certified segregation comparisons

The verifier supplies rational lower and upper witnesses using (11)–(12)
with the appropriate kernel. Entries below are completed from its reproducible
output; a witness at a finite cutoff includes the all-type tail inequality.

| Segregation law | Exact parameter | Certified interval for $c_*$ | Lower / upper certificate size |
| --- | --- | --- | ---: |
| Binomial | $p=1/2$ | $(0.378474994554,\ 0.378474994555)$ | 25 / 48 |
| Beta-binomial | $\rho=1/41$ | $(0.367,\ 0.368)$ | 8 / 10 |
| Beta-binomial | $\rho=1/4$ | $(0.257,\ 0.258)$ | 14 / 20 |
| Beta-binomial | $\rho=1/3$ | $(0.216,\ 0.217)$ | 16 / 24 |
| Beta-binomial | $\rho=1/2$ | $(0.137,\ 0.138)$ | 32 / 36 |
| Asymmetric binomial | $p=3/4$ | $(0.270,\ 0.271)$ | 20 / 20 |
| Asymmetric binomial | $p=9/10$ | $(0.120,\ 0.122)$ | 28 / 36 |

Every decimal endpoint is an exact rational. All upper witnesses have strict
tail slack, so the interval endpoints are strict. In particular, at the
**same** $c=13/50=0.26$, beta $\rho=1/4$ eradicates (size 14
upper witness), while bias $p=3/4$ persists (size 8 lower witness).
These are infinite-time conclusions proved by certificates; they do not come
from finite simulation survival frequencies.

### A falsified attempted monotonicity proof

An attempted proof asserted that more clustering increased finite-generation
extinction from every founder. Exact arithmetic refutes it. At $c=3/8$,
generation two, compare binomial with beta $a=20$, $\rho=1/41$:

| Founder | Binomial extinction by generation 2 | Clustered extinction by generation 2 | Direction |
| --- | ---: | ---: | --- |
| $k=1$ | approximately 0.455618761 | approximately 0.458760147 | Clustering increases finite-time extinction |
| $k=5$ | approximately 0.790885981 | approximately 0.789199561 | Clustering decreases finite-time extinction |

These are counterexamples to a proposed coupling, not evidence against
Theorem 5 or a proof against monotonicity of the asymptotic threshold.
Selection can favor occasional low-copy daughters of a high-copy mother.
Convex-order claims about segregation alone do not control the nonmonotone
reproductive-value eigenfunction.

## 6. Sharp predictions and threshold identifiability

### Theorem 6 — conditional scaling and dose response

For a uniformly chosen sister, for every mother count $k\geq1$,

$$
\mathbb E[J\mid k]=k,\qquad
\operatorname{Var}(J\mid k)=\frac{1-R}{2}k+Rk^2,
\qquad F_{\mathrm{conditional}}(k)=\frac{1-R}{2}+Rk.
\tag{20}
$$

For any selected mother law of mean $\mu$, variance $V$, daughter
mean and variance are
$\mu_d=\mu$ and
$V_d=(1+R)V+R\mu^2+(1-R)\mu/2$.

For any **fixed finite-support input mother distribution** $\pi$, let
$p_k(c)=(1+cg(k))^{-1}$ with $g$ nonnegative nondecreasing and
$g(0)=0$. Write $\pi_c\propto\pi p(c)$,
$z(k)=\mathbb E[P^{2k}]$, $z(0)=1$. The pooled daughter zero
fraction satisfies

$$
\frac{d f_d^0}{dc}=
-\operatorname{Cov}_{\pi_c}\left(z(k),\frac{g(k)}{1+cg(k)}\right)
\geq0.
\tag{21}
$$

For linear killing, start independent replicates with any positive finite number of identical
**one-copy founders** and pool all generation-two cells, including zeros.
For every $c\geq0$, $0\leq R\leq1$, define

$$
A=2+(3-R)c,\quad D=2+4c+(1+R)c^2,\quad B=3+R+4c.
$$

Then

$$
\boxed{\mu_2=\frac AD,\qquad
F_2=(1+R)\frac BA+\frac{1-R}{2}-\frac AD.}
\tag{22}
$$

The mean strictly decreases and Fano factor strictly increases with $c$,
including in a neighborhood of every positive critical value. These curves
are smooth at the threshold; there is no claimed finite-generation phase
singularity. At $c=0$,
$F_2=1+3R/2+R^2/2$; as $c\to\infty$,
$F_2\to4(1+R)/(3-R)+(1-R)/2$. For binomial this limit is $11/6$,
different from PR #55's exponential-law result.

**Proof.** Conditional total variance for a symmetric $P$ gives
$\mathbb E[2kP(1-P)]+\operatorname{Var}(2kP)
=(1-R)k/2+Rk^2$. Another total-variance decomposition proves the
selected-mother formula. Differentiating the normalized weights $\pi_c$
gives (21). The function $z(k)$ is nonincreasing while
$g(k)/(1+cg(k))$ is nondecreasing. Twice their covariance equals the
expectation of
$(z(K)-z(L))(u(K)-u(L))\leq0$ for independent $K,L\sim\pi_c$.
This proves the sign; it is strict if two positive-mass counts make both
inequalities strict.

For (22), first-generation marginal probabilities are
$w_0=w_2=(1+R)/4$, $w_1=(1-R)/2$. Selecting those mothers by
$1,1/(1+c),1/(1+2c)$ gives mean $A/D$ and raw second moment
$B/D$. The founder's division probability cancels from every normalized
count. Applying the variance formula yields (22). To prove the derivative
signs without relying on a plot, put $C=2+4c+(3-R)c^2$. Direct algebra
gives

$$
\mu_2'=-(1+R)C/D^2<0,\qquad
F_2'=(1+R)\left[\frac C{D^2}-\frac{1-R^2}{A^2}\right]>0.
\tag{23}
$$

Indeed $C\geq D$, and

$$
A^2-(1-R^2)D=(2+2R^2)+(8-4R+4R^2)c
 +(8-7R+2R^2+R^3)c^2>0.
$$

All displayed coefficients are positive on $[0,1]$; for the last,
$8-7R\geq1$. Thus $CA^2>(1-R^2)D^2$. Taking endpoints in
(22) proves the limits. ∎

### Theorem 7 — identical Fano curves can conceal different thresholds

Compare beta clustering with $\rho=1/4$ ($a=3/2$) against
asymmetric allocation with $p=3/4$. They have identical conditional
variance curves (20) for **all** mother counts, and identical entire
$\mu_2(c),F_2(c)$ curves (22). Nevertheless their linear-killing
thresholds differ, as certified by the separated witnesses in Section 5.
Consequently even a noiseless, continuously measured Fano dose curve of this
assay does not identify the eradication threshold across these two families.

An additional statistic for a **uniformly selected or pooled sister** distinguishes
them. For $k\geq2$, write
$(x)_{\underline4}=x(x-1)(x-2)(x-3)$. For general matched
$0<R<1$,

$$
T_4=\frac{16\mathbb E[(J)_{\underline4}\mid k]}{(2k)_{\underline4}}
-1-6R
=\begin{cases}
3R^2/(1+2R),&\text{beta clustering},\\
R^2,&\text{fixed asymmetric allocation}.
\end{cases}
\tag{24}
$$

At $R=1/4$, these values are $1/8$ and $1/16$. Alternatively,
in the generation-two one-copy assay, the pooled zero fractions differ by

$$
f^0_{2,\mathrm{beta}}-f^0_{2,\mathrm{bias}}
=\frac{R^2(1-R^2)(1+c)}{16(1+2R)D}>0.
\tag{25}
$$

For $R=1/4$, this is $5(1+c)/(2048D)$. The distinction is
small enough that measurement noise and independent replicate counts matter.

**Proof.** Both have the same $R$, so Theorem 6 gives the equality
of curves. The exact survival/eradication witnesses in Section 5 prove
different thresholds without floating-point eigenvalue estimates.
Conditional binomial factorial moments give
$\mathbb E[(J)_{\underline4}\mid k]/(2k)_{\underline4}=\mathbb EP^4$.
For symmetric $X=2P-1$,
$16\mathbb EP^4=1+6R+\mathbb EX^4$. Beta moments from (2) yield
$\mathbb EX^4=3R^2/(1+2R)$; for the equally weighted fixed bias,
$X=\pm\sqrt R$, so $\mathbb EX^4=R^2$. Their difference is
$2R^2(1-R)/(1+2R)>0$, proving (24).

The selected first-generation mother laws are identical. Only their
$k=2$ zero-daughter probabilities differ, by
$R^2(1-R)/[8(1+2R)]$. The selected mass at $k=2$ is
$(1+R)(1+c)/(2D)$. Multiplication proves (25). This distinguishes
the two specified families; finitely many moments do not identify an
arbitrary segregation distribution. ∎

### An experimentally testable, nonclinical use of these predictions

Use independent cultured microcolonies or tracked divisions as the sampling
units. Measure maternal counts, paired daughter counts, divisions, deaths,
and censoring. Check $J_1+J_2=2k$ first: retention loss, incomplete
replication, or copy-count error can directly reject the model. Estimate
$p_k$ from observed fates; a linear hazard requires
$p_k^{-1}-1=ck$. Drug concentration is not automatically $c$: the
division rate and death law require independent calibration at each condition.

Conditional daughter Fano versus maternal count should be a line of slope
$R$, with intercept $(1-R)/2$. Test the mean/Fano dose laws (22)
only when the one-copy founder design is actually met. For mixtures, (21)
compares the **same input distribution** across conditions; it is not a claim
about an evolving culture's time trajectory. Once segregation changes with
drug, a total derivative includes its parameter response and need not have
the partial-derivative sign in (23).

Threshold inference additionally requires distinguishing the segregation
family. Measure (24) at $k\geq2$, or the small zero-fraction difference
(25), rather than inferring the threshold from second moments alone. Keep
extinct wells in population averages, include zero-copy cells, and compute
uncertainty across independent wells or independent maternal divisions.
Related sisters are not independent replicates. A finite surviving well does
not disprove eventual extinction, and an extinct well does not disprove a
positive survival probability. Theorems specify an infinite-time model event,
not a finite assay's detection limit.

## 7. Verification, falsification, and hostile self-review

`verifier.py` uses exact `fractions.Fraction` arithmetic. `simulation.py`
independently samples division/death and splits with `random.Random`; it
imports neither the verifier nor its matrices or moment formulas. The new
test independently constructs kernels using Bernoulli bit sequences or a
Pólya urn, enumerates complete small-generation populations, checks rational
certificates and analytic identities, and compares fixed-seed independent
replicate populations with finite-generation predictions. Tests can falsify
an implementation; the infinite-state claims rest on the proofs above.

| Hostile challenge | Resolution or retained limitation |
| --- | --- |
| A finite truncation has a subcritical eigenvalue, therefore the full model dies. | Rejected. Every upper certificate uses the all-type tail inequality (12); convergence of corrected bounds is proved. |
| Growing expectation proves survival. | Rejected. Finite positive Perron witnesses construct an extinction supersolution and actual surviving subprocess. |
| Criticality can be decided by rounded eigenvalues. | Rejected. Equality extinction follows from the survival equation paired with a strictly positive summable left eigenvector. Rational witnesses decide the reported noncritical points. |
| Compactness applies to saturated hazards too. | Rejected: their division probabilities need not tend to zero. Separate bounds and the refuge lemma are used. |
| More clustering always raises extinction probability at every horizon. | Refuted by the two exact opposite-sign founder examples in Section 5. Full threshold monotonicity remains a labelled conjecture. |
| Positive survival without death guarantees a positive threshold for every infinite-type kernel. | Rejected. Complete clumping is a counterexample. The general clustered criterion allows a zero threshold; examples with positive thresholds have explicit finite witnesses. |
| Matching the Fano dose curve identifies the killing threshold. | Refuted by Theorem 7's matched-moment mechanisms and separated exact certificates. |
| Copy-number saturation means eradication is impossible. | Rejected. A positive hazard floor scaled without bound gives finite eradication; a low dose ceiling is a different question. |
| The saturation threshold is always its high-copy division rate. | Refuted for every $H>2$ by (16). |
| Sisters can be sampled independently because their marginals agree. | Rejected. Exact enumeration and simulation preserve their sum; the quadratic extinction/moment terms depend on their joint split. |
| Beta-binomial segregation or its variance law is new to ecDNA. | Rejected after inspecting Shibata et al. The narrower theorem and assay-identifiability contribution is stated instead. |
| Retention-element work validates a copy-conserving bias parameter. | Rejected. Nuclear retention loss violates the copy-conservation assumption and needs another model. |
| Simulated caps can be treated as extinction. | Rejected. Cell/copy resource limits abort with an error and return no censored outcome. |

### Recorded checks and simulation outcomes

The standalone verifier completes **1,589 exact checks**, including the
linear bracket, six segregation brackets, the same-moment counterexample,
general dose identities, and rational refuge/scaling witnesses. The test file
contains **16 tests**, including independent complete-population enumeration,
independent Pascal/urn kernels, and fixed-seed stochastic checks. The initial
complete run took about five seconds. Repository schema validation passed;
the full Python suite, including this file, passed **271 tests in 23.901
seconds**, and all **26 Node tests** passed. The PR records the final remote
CI status. This new test file uses no network or third-party package and stays
below the requested 60-second budget in this workspace.

| Theorem group | Exact/independent numerical check | Stochastic compatibility check |
| --- | --- | --- |
| 1–2: infinite-state boundary | Exact linear systems, independently rebuilt Pascal rows, residual signs, and all-type tail inequalities; a one-type false-extinction shortcut is rejected. | Finite-horizon means, extinction probabilities, and moments on both sides of the linear bracket. |
| 3: general powers and limits | Quadratic certificates; steep-power certificates on both sides of the proved limiting coefficient; rational sine-band bound. | Quadratic and steep-power genealogies; constant killing at its critical parameter; declining killing from a high-copy founder. |
| 4: saturation/refuges | Saturation determinant over several rational parameters, low-copy death floor, and dose-ceiling inequalities. | Saturation with $H=3,c=1$; declining hazards under binomial, beta, asymmetric, and fully clumped segregation. |
| 5: clustering/asymmetry | Six two-sided full-model intervals, universal count inequalities, symmetry, and exact finite-time ordering reversal. | Clustered and biased split experiments, full clumping, and finite populations at the same $c=0.26$. |
| 6–7: falsifiable observables | Exact pooled moments from full population laws, derivative identities, fourth moments, identical Fano curves, and the exact zero-fraction gap. | Independent maternal splits for second/fourth moments and whole replicate populations for dose/Fano predictions. |

No finite stochastic experiment checks the mathematical equality case or
infinite-time survival directly. Those rely on the proof. The simulator is
used to challenge their finite-horizon consequences and model implementation.

The following synthetic values are reproduced by the standalone simulator
with `--trials 10000`; exact predictions use the verifier's independent
finite-generation recursion. The founder is $k=1$, generation is five,
and SE is computed across independent replicate populations.

| Law/kernel | Expected positive cells (exact prediction, rounded) | Simulated mean | SE | Seed |
| --- | ---: | ---: | ---: | ---: |
| Linear, binomial, $c=3/8$ | 1.352533298 | 1.361700000 | 0.023629286 | 2026111701 |
| Linear, binomial, $c=2/5$ | 1.210808544 | 1.217800000 | 0.022062298 | 2026111702 |
| Linear, beta $\rho=1/3,c=1/4$ | 1.154311809 | 1.150400000 | 0.019462679 | 2026111703 |
| Linear, bias $p=3/4,c=3/10$ | 1.123681126 | 1.104800000 | 0.019724619 | 2026111704 |
| Quadratic, binomial, $c=3/10$ | 1.090109264 | 1.095500000 | 0.021375773 | 2026111705 |
| Saturating, binomial, $c=2,H=3$ | 1.033840237 | 1.016300000 | 0.019725677 | 2026111706 |

At $c=0.3$, generation two, the same-variance beta/bias pair has the
same exact Fano factor, approximately **1.491196360**. Independent simulations
give **1.502431537 (SE 0.012651343)** for beta (seed 2026111721) and
**1.487655269 (SE 0.012470202)** for bias (seed 2026111722). Their exact
zero fractions are approximately **0.545902123** and **0.544943986**.
The simulated zero fractions are **0.547008200 (SE 0.002856521)** and
**0.544069463 (SE 0.002858677)**. This sample size is compatible with
the formulas but does not resolve that tiny zero-fraction gap statistically.
Exact arithmetic establishes the gap; power and measurement error would need
attention in an actual assay.

Hostile review independently checked every displayed proof, signs, boundaries,
and the numerical certificates before pushing. It corrected the attribution
of PR #55's lower endpoint, required explicit strict slack for open numerical
intervals, and exposed a typed-cache input-validation issue (now fixed and
regression-tested). An initial test fixture confused three founders' cell
count with their total copy mass; another incorrectly expected every
supercritical point to yield a negative solution at one chosen cutoff. Both
test assumptions were corrected; the verifier correctly returned
`undetermined` when that particular witness did not decide the model. No
failed scientific comparison or unproved monotonicity assertion was hidden.
This is mathematical/computational self-review, not external scientific review.

### Limitations and unresolved questions

The killing law is assumed globally over all copy numbers, including counts
not observed in a finite experiment. No finite dataset establishes such a
global tail law. Copy number alone need not determine expression or fate.
Segregation, birth rates, and death rates may change during adaptation, while
the theorems use time-independent parameters and independent cells.

The general beta/asymmetry threshold is characterized exactly, but no
elementary formula or universal interior monotonicity theorem is claimed.
We do not establish positive critical values for every nondegenerate beta
parameter independently of the spectral criterion. General saturating-law
critical values and equality cases between the stated bounds remain open.
These are explicit limits of the contribution, not simulations promoted to
theorems. No universal calendar-time Fano formula, critical survival-tail
exponent, or finite-time probability of clinical eradication is claimed.

## 8. Reproduce and scope

From the repository root, offline, using the Python standard library:

```bash
python research/ecdna-killing-laws-and-clustering/verifier.py
python research/ecdna-killing-laws-and-clustering/verifier.py --linear-only
python research/ecdna-killing-laws-and-clustering/simulation.py --trials 10000
python -m unittest discover -s tests -p 'test_research_ecdna_killing_laws_and_clustering.py' -v
```

The exact checker has no NumPy/SciPy or other third-party dependencies. Rational
parameters are required; floating-point inputs are rejected. Nonrational
exponents are covered by the proofs, while executable exact checks use
rationally evaluable cases and bounds. Monte Carlo uses floating point and
recorded seeds; extreme hazards or beta parameters can be affected by
floating-point underflow or loss of resolution. Extreme-parameter claims rest
on the proofs and exact checker. It is a separate finite-horizon compatibility check, not an
interval proof or evidence of infinite-time survival.

Repository integration gates, with existing validation dependencies, are:

```bash
python scripts/validate.py
python -m unittest discover -s tests -v
node --test tests/*.test.*
```

This contribution adds only this directory's `README.md`, `verifier.py`,
`simulation.py`, and
`tests/test_research_ecdna_killing_laws_and_clustering.py`. No existing
research note, top-level README, index, policy, catalog, workflow, dependency,
script, result, or existing test is changed.
