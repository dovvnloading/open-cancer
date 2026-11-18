# ecDNA segregation under copy-selective cell killing

**Nonclinical computational mathematics. No treatment advice.** This is an
unreviewed mathematical research contribution, not an experimentally validated
description of a cancer, a drug, or a patient. All numerical examples are
synthetic. The proofs establish statements about the specified model; neither
AI review nor passing tests establishes biological validity or scientific
priority. Literature search date: **17 November 2026**.

## Result and explicit pivot

Binomial ecDNA segregation, branching-process models, copy-dependent selection,
and the neutral constant-mean/linear-variance law already exist. We therefore
pivot from those questions to an explicitly solvable **unbounded-copy model of
whole-cell killing**. A type-$k$ cell divides before dying with probability

$$
p_k=u s^k,\qquad 0<u\leq1,\quad0<s<1.
$$

The exact persistence boundary of the entire ecDNA-positive lineage is

$$
\boxed{\lambda(u,s)=2u(1-\sqrt{1-s})=1.}
$$

Extinction is certain at equality. This holds for **every finite initial
population containing a positive cell, and every finite initial positive copy
number**, without a copy-number cap.
For $u=1$, the boundary is $s=3/4$, or
$\theta=\log(4/3)$ if $s=e^{-\theta}$. These are abstract model parameters,
not pharmacological doses. The theorem covers both arbitrary bounded birth-rate
clocks and an explicit family with increasing copy-dependent birth and death.

A second result gives an exact variance-to-mean time course independent of the
initial copy number for monomorphic founders, and a certified universal Fano
limit for arbitrary finite mixtures. A particularly simple consequence is

$$
\boxed{F_2=\frac{3+s}{2(1+s)},\qquad
s=\frac{3-2F_2}{2F_2-1}.}
$$

Here $F_2$ is measured across **all generation-two cells pooled across
independent replicate populations**, including ecDNA-free cells. It is not the
mean of within-population ratios. The exact benchmark also supplies comparison
certificates for general killing laws: for constant division rate $b$ and
linear killing $d_k=\gamma k$, persistence is proved for $\gamma<b/3$ and
eradication for $\gamma\geq b/2$. The interval between these bounds remains
unresolved here.

**Novelty assessment:** targeted searches did not locate this explicit
unbounded-copy threshold, its critical-case proof, or this therapy-dependent
Fano package for the stated kernel. This is a provisional, narrowly scoped
priority assessment. The offspring model, generating-function techniques,
minimal extinction fixed point, and monotone coupling are established tools.
We do not claim a first ecDNA branching model or a universal cancer threshold.

## 1. Literature audit before model selection

The following sources were checked for overlap or biological limitations. URLs
identify sources rather than redistribute their data. The audit uses accessible full texts, methods, supplements, and abstracts;
historical kernels were not exhaustively checked for equivalence.

| Source | Already known or observed | Consequence here |
| --- | --- | --- |
| [Turner et al., Nature 2017](https://doi.org/10.1038/nature21356), [methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC5334176/) | Discrete-generation branching, copy-dependent division/death, replicated binomial segregation, and a high-copy stopping cutoff. | Our model ingredients are prior art. The exact unbounded-copy theorem is the proposed contribution. |
| [Pichugin, Huang & Werner, 2019 preprint](https://doi.org/10.1101/2019.12.15.876714); [Lange et al., Nature Genetics 2022](https://doi.org/10.1038/s41588-022-01177-x), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9534767/) | Neutral mean and variance dynamics, random segregation, selection, and changing ecDNA-positive fractions. | The neutral variance law below is a consistency check, not a new claim. Relative abundance and extinction of the whole positive lineage must be distinguished. |
| [Kimmel & Axelrod, Genetics 1990](https://doi.org/10.1093/genetics/125.3.633) | Branching models of unstable amplification/double minutes and loss of drug resistance. | Historical amplification theory predates current ecDNA terminology. |
| [Harnevo & Agur, Mathematical Biosciences 1991](https://doi.org/10.1016/0025-5564(91)90094-y), [abstract](https://pubmed.ncbi.nlm.nih.gov/1804437/) | Copy-dependent cell loss, extinction calculations for a single-copy compartment, and dose-sensitive prognosis. | Abstract-level overlap is substantial; the complete historical kernel was not available for an exhaustive equivalence check. |
| [Kimmel & Stivers, Bulletin of Mathematical Biology 1994](https://doi.org/10.1007/BF02460646) | Continuous-time branching walks and explicit expectations for unstable amplification. | Continuous time and exact expectations alone are not new. |
| [Santer & Uecker, Genetics 2020](https://doi.org/10.1534/genetics.119.303012), [Appendix A](https://pmc.ncbi.nlm.nih.gov/articles/PMC7337077/) | Multicopy-plasmid establishment via the first-event extinction equation with correlated daughters. | Theorem 1 is standard branching machinery, included with a complete proof and implementation. |
| [Li, Kibby & Foo, Journal of Mathematical Biology 2023](https://doi.org/10.1007/s00285-023-01992-8) | Amplification-driven resistance, extinction/escape conditions, and recurrence. | We claim only the specific replicated-binomial, exponential-killing solution, not the first amplification extinction criterion. |
| [Nathanson et al., Science 2014](https://doi.org/10.1126/science.1241328) | Dynamic extrachromosomal mutant EGFR regulation under treatment and recovery. | Reacquisition, reintegration, and within-cell copy changes are omitted assumptions, not established impossibilities. |
| [Kim et al., Nature Genetics 2020](https://doi.org/10.1038/s41588-020-0678-2) | Pan-cancer ecDNA prevalence and outcome associations. | Motivates population heterogeneity; does not calibrate our killing law. |
| [Tang et al., Nature 2024](https://doi.org/10.1038/s41586-024-07802-5) | Preclinical preferential killing of ecDNA-containing cancers by enhanced transcription–replication conflict/CHK1 inhibition. | Motivates selective cell death; does not establish $p_k=u s^k$ or a linear copy-dependent hazard. |
| [Montuori et al., Cancer Discovery 2025](https://doi.org/10.1158/2159-8290.CD-24-1738), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12456741/) | Copy-dependent fitness models and high-copy elimination/low-copy tolerance under genotoxic pressure. | Neither preferential high-copy removal nor its adaptive consequences is claimed as new. |
| [Wong, Bailey et al., Cell 2026 review](https://doi.org/10.1016/j.cell.2026.03.011) | Review of ecDNA biology across cancer hallmarks and therapeutic vulnerabilities. | Context, not independent proof of our mathematical or biological assumptions. |
| [Retention elements, Nature 2026; online November 2025](https://doi.org/10.1038/s41586-025-09764-8); [Shibata et al., Nature Communications 2026](https://doi.org/10.1038/s41467-026-76566-5) | Retention mechanisms and active changes in segregation during drug adaptation. | Exact unbiased binomial segregation is a falsifiable restriction, especially under treatment. |
| [Lee et al., 2026 preprint](https://doi.org/10.64898/2026.02.12.705658), [April 30 version](https://pmc.ncbi.nlm.nih.gov/articles/PMC12918801/); [Scanu, Werner & Huang, 2025 preprint](https://doi.org/10.1101/2024.10.22.619675) | Copy-dependent death/density-regulated simulations; multiple ecDNA variants and de-mixing. | These model classes are broader biologically; neither inspected source supplied the exact law proved here. Both are identified as preprints. |

Adjacent mathematical prior art includes
[Bansaye's contaminated-cell branching model](https://doi.org/10.1214/07-AAP465),
[branching within branching](https://doi.org/10.1016/j.spa.2015.12.007),
[burden-dependent parasite birth/death](https://doi.org/10.1016/j.spa.2023.104262),
and [exact intracellular-component partitioning](https://pmc.ncbi.nlm.nih.gov/articles/PMC4550007/).
Generating-function composition and exponential tilting belong to this wider
literature. We have not proved that a differently worded historical theorem
cannot imply the present result.

**Dated trial context:** ecDNA-directed compounds have entered early-phase
trials. A [November 2026 company SEC filing](https://investors.boundlessbio.com/node/7956/html)
reports cessation of POTENTIATE enrollment in January 2026 and that June 2026
BBI-940 exposure findings did not support continued development. Thus “now in
trials” should not be repeated as an unqualified current claim. No clinical
effectiveness claim or drug-specific parameter is used in this work.

Searches included combinations of “ecDNA extinction threshold”, “binomial
segregation therapy”, “gene amplification copy dependent death”, “ecDNA
exponential fitness”, “multiplicative fitness plasmid”, “binomial segregation
extinction 3/4”, and “parasite exponential death branching cells”. Sources and
newer citations were followed across ecDNA, gene-amplification, plasmid, and
host–parasite theory. Search-engine coverage, inaccessible full texts, and
unpublished work prevent an exhaustive priority claim.

## 2. Model, assumptions, and observables

All theorems use the following assumptions unless their statement explicitly
narrows them.

1. The initial state is a **finite deterministic multiset**
   $\mathbf k=(k_1,\ldots,k_N)$, $N\geq0$, $k_i\in\mathbb Z_{\geq0}$.
   There is no uniform upper bound on allowed $k_i$ or subsequent types.
2. A type-$k$ cell has competing, time-independent exponential division and
   death rates $b_k>0$, $d_k\geq0$, both finite at each type. There is a finite
   constant $B$ with $b_k\leq B$ for every $k$. Cells act independently
   conditional on their types. No density dependence, migration, or shared
   environmental randomness is included.
3. At division the mother is replaced by two daughters. Exactly $k$ copies
   replicate to $2k$. Draw $J\sim\operatorname{Bin}(2k,1/2)$; the daughters
   receive **$J$ and $2k-J$**. The sisters are complementary, not independent
   binomial samples. The complete binomial kernel is
   $B_{kj}=\binom{2k}{j}4^{-k}$, $0\leq j\leq2k$.
4. Copies change only at division. A zero-copy cell cannot regenerate ecDNA.
   No copy clustering, biased partition, loss before division, chromosomal
   reintegration, new ecDNA formation, or mutations are included.
5. A **genealogical generation** counts divisions along ancestry. A cell that
   dies contributes no children to the next generation. Its probability of
   dividing before dying is $p_k=b_k/(b_k+d_k)$. Generational results are not
   formulas at fixed elapsed calendar time.
6. For Theorems 2–5, specialize to $p_k=u s^k$, $0<u\leq1$, $0<s<1$.
   This is an assumed log-affine division-before-death law, not an inference
   from a particular drug. Its continuous-time realization is
   $$
   d_k=b_k[(u s^k)^{-1}-1].
   $$
   To make both division and killing copy-dependent, one example is
   $$
   b_k=b\left(1+\eta\frac{k}{h+k}\right),\quad
   b>0,\ h>0,\ \eta>0.
   $$
   This is bounded, strictly increasing division; the displayed $d_k$ is
   also strictly increasing. Constant $b_k=b$ is the $\eta=0$ case.
   With $s=e^{-\theta}$, killing is
   $b_k(u^{-1}e^{\theta k}-1)$, **exponential**, not linear, in copy number
   when $b_k$ is constant.
7. $Z_n(j)$ counts type-$j$ cells at genealogical generation $n$;
   $N_n=\sum_jZ_n(j)$, $Z_n^+=\sum_{j\geq1}Z_n(j)$, and
   $S_0=\sum_i k_i$. Eradication means the entire positive lineage eventually
   becomes empty. This differs from tumor extinction and from a small positive
   fraction. The **pooled cell law** is
   $$
   \pi_n(j)=\frac{\mathbb E Z_n(j)}{\mathbb E N_n}.
   $$
   It describes pooling every cell across many independent replicate
   populations at the same genealogical generation. It is not
   $\mathbb E[Z_n(j)/N_n\mid N_n>0]$. Pooled means and variances below are
   under $\pi_n$. Between-population variances of sums are provided separately.

**Existence and time interpretation.** The number of divisions up to any
finite time is dominated by a Yule process with rate $B$ per cell. For completeness, stop that pure-birth
process when its population first reaches $M$. Its stopped expectation is
at most $N e^{Bt}$ by the integral inequality for the rate $B$ per cell.
Thus the probability of reaching $M$ by time $t$ is at most
$N e^{Bt}/M$, which tends to zero; the process is nonexplosive.
Deaths cannot increase that birth count; their number is at most the initial
cell count plus the number of divisions. Thus the model has no finite-time
event explosion even though $d_k$ can be unbounded over types. Every cell's
total rate is finite and positive. A finite positive genealogy dies out in
finite time almost surely; an infinite positive genealogy cannot complete in
finite time. The extinction event therefore agrees between the continuous-time
process and its genealogical tree. The birth-rate clock affects time to
eradication and calendar-time composition, which are not solved here.

## 3. Numbered theorems and complete proofs

### Theorem 1 — exact extinction and population-sum recursions

Assume 1–5 and 7; $p_k$ may be any such division-before-death probabilities.
Set $q_0^{(n)}=1$ for every $n$, and $q_k^{(0)}=0$ for $k\geq1$. Then

$$
q_k^{(n+1)}=1-p_k+p_k\sum_{j=0}^{2k}B_{kj}
q_j^{(n)}q_{2k-j}^{(n)},\qquad k\geq1.\tag{1}
$$

For every finite initial multiset,

$$
\Pr_{\mathbf k}(Z_n^+=0)=\prod_iq_{k_i}^{(n)},\qquad
\Pr_{\mathbf k}(\text{eradication})=\prod_iq_{k_i},\quad
q_k=\lim_{n\to\infty}q_k^{(n)}.\tag{2}
$$

The vector $q$ is the coordinatewise least fixed point of (1) in $[0,1]$
with $q_0=1$. This is an exact characterization; we do not assert an elementary
closed form for $q_k$ in the persistence regime.

For a real-valued statistic $f(k)$, let $X_n=\sum_j f(j)Z_n(j)$, starting
from one type $k$. Set $m_{0,k}=f(k)$, $v_{0,k}=f(k)^2$. The **raw** first
and second moments satisfy

$$
\begin{split}
m_{n+1,k}&=p_k\sum_jB_{kj}(m_{n,j}+m_{n,2k-j}),\\
v_{n+1,k}&=p_k\sum_jB_{kj}
(v_{n,j}+v_{n,2k-j}+2m_{n,j}m_{n,2k-j}).\tag{3}
\end{split}
$$

For multiple founders, the mean is $\sum_i m_{n,k_i}$ and variance is
$\sum_i(v_{n,k_i}-m_{n,k_i}^2)$. The code implements cell, positive-cell,
zero-cell, copy-mass, and copy-square statistics.

**Proof.** Condition on the first event. Death has no positive descendants;
division has the stated joint daughter types. Given those types, their future
subtrees are independent, which gives the product in (1), not the product of
two separately averaged daughter probabilities. At generation zero a positive
founder is present. The events $\{Z_n^+=0\}$ increase in $n$, because zero
cannot regenerate ecDNA. Their union is eventual eradication; monotone
convergence gives the limit. Each row of (1) is a finite sum, so the limit is a
fixed point. Iteration from the smallest admissible vector stays below any
other fixed point, proving minimality. Independent founders give (2).

For (3), conditional on the split, square the sum of the two future subtree
statistics. Independence gives the cross term $2m_{n,j}m_{n,2k-j}$. Founder
variances add. For a nonempty initial state, at finite $n$, at most $N2^n$ cells and types
up to $2^n\max_i k_i$ are reachable; the empty case is trivial, so all these computations are finite and exact.
This reachable-support bound is not an imposed copy-number cutoff. ∎

### Theorem 2 — exact generating function, means, and variances

Assume 1–7, $N\geq1$, and define

$$
H(z)=\sum_{i=1}^Nz^{k_i},\quad T(z)=\frac{s(1+z)^2}{4},\quad
G_n(z)=\sum_{j\geq0}\mathbb EZ_n(j)z^j.
$$

For every $n\geq0$,

$$
G_n(z)=(2u)^nH(T^{\circ n}(z)).\tag{4}
$$

Write $A_n=T^{\circ n}(1)$, $B_n=T^{\circ n}(0)$,
$D_n=(T^{\circ n})'(1)$, $E_n=(T^{\circ n})''(1)$. Then

$$
\begin{array}{ll}
\mathbb EN_n=(2u)^nH(A_n),&\mathbb EZ_n(0)=(2u)^nH(B_n),\\
\mathbb EZ_n^+=(2u)^n[H(A_n)-H(B_n)],&
\mathbb E\sum_jjZ_n(j)=(2u)^nH'(A_n)D_n.
\end{array}\tag{5}
$$

The pooled copy-number mean $\mu_n$, raw second moment $R_n$, variance
$V_n$, and zero fraction $f_n^0$ are exactly

$$
\begin{split}
\mu_n&=\frac{H'(A_n)D_n}{H(A_n)},\\
R_n&=\frac{H''(A_n)D_n^2+H'(A_n)(E_n+D_n)}{H(A_n)},\\
V_n&=R_n-\mu_n^2,\qquad f_n^0=\frac{H(B_n)}{H(A_n)}.\tag{6}
\end{split}
$$

Positive-only pooled moments use $H(A_n)-H(B_n)$ in place of the denominator
$H(A_n)$. If there are no positive founders those moments are undefined.
All-cell moments are undefined for an empty initial population, whose expected
counts are all zero. Normalized composition is independent of $u$.

**Proof.** One type-$k$ mother contributes the expected polynomial

$$
u s^k\mathbb E(z^J+z^{2k-J})
=2u\left[s\left(\frac{1+z}{2}\right)^2\right]^k.
$$

Summing over mothers gives $G_{n+1}(z)=2uG_n(T(z))$ and induction proves
(4). Evaluate at 1 and 0 for counts. Differentiate once and twice at 1;
$G_n''(1)+G_n'(1)$ is the expected sum of squared copy numbers. Divide by
$G_n(1)>0$ to obtain (6). Conditioning the pooled law on positive copies
changes only the denominator of the raw moments. The finite polynomials make
all differentiations legitimate. For computation,

$$
\begin{split}
A_0=D_0=1,\quad E_0=0;\qquad
A_{n+1}&=s(1+A_n)^2/4,\\
D_{n+1}&=s(1+A_n)D_n/2,\\
E_{n+1}&=sD_n^2/2+s(1+A_n)E_n/2.
\end{split}\tag{7}
$$

The factor $(2u)^n$ cancels only in normalized composition, not in population
counts or extinction probabilities. ∎

### Theorem 3 — sharp eradication threshold, including criticality

Assume 1–7 and at least one positive founder. Set

$$
a=\frac{1-\sqrt{1-s}}{1+\sqrt{1-s}},\quad
r=\frac{2a}{1+a}=1-\sqrt{1-s},\quad
\lambda=2ur=\frac{4ua}{1+a}.
$$

Then $0<a,r<1$, $T(a)=a$, and:

1. The expected positive population satisfies, for all $n\geq0$,
   $$
   \mathbb E Z_n^+\leq3S_0\lambda^n.\tag{8}
   $$
   In fact $\mathbb E Z_n^+\sim\kappa(s)H'(a)\lambda^n$, for a finite
   strictly positive constant $\kappa(s)$ specified in the proof.
2. Eradication has probability 1 **if and only if** $\lambda\leq1$.
   When $\lambda<1$,
   $$
   \Pr(Z_n^+>0)\leq\min(1,3S_0\lambda^n).\tag{9}
   $$
   At criticality extinction is also certain, although (9) supplies no
   vanishing bound there.
3. When $\lambda>1$, define
   $$
   h(0)=0,\quad h(k)=ka^{k-1},\quad
   H_*=(1-a)^{-2},\quad C_*=\frac{us}{(1-sa)^2},\quad
   \varepsilon=\min\{H_*^{-1},(\lambda-1)/C_*\}.
   $$
   The explicit finite-population survival certificate is
   $$
   \Pr_{\mathbf k}(\text{positive lineage survives forever})
   \geq1-\prod_i[1-\varepsilon h(k_i)]>0.\tag{10}
   $$

Equivalently, no $s<1$ permits persistence if $u\leq1/2$. For $u>1/2$,
persistence occurs exactly when

$$
s>\frac{4u-1}{4u^2}.
$$

**Proof of (8) and the asymptotic.** The second fixed point of $T$ is $1/a>1$.
Thus $A_n\downarrow a$ and $B_n\uparrow a$. More quantitatively, put
$c=s/4$, $x_n=A_n-a$, $y_n=a-B_n$, and
$\rho=r+c(1-a)=a(3+a)/(1+a)^2<1$. Direct expansion gives

$$
x_{n+1}=rx_n+cx_n^2,\quad y_{n+1}=ry_n-cy_n^2,
\quad x_n\leq(1-a)\rho^n,\quad y_n\leq ar^n.
$$

Here the bounds follow by induction using $0\leq x_n\leq1-a$ and
$0\leq y_n\leq a$. For $\Delta_n=A_n-B_n$,

$$
\Delta_{n+1}=r\Delta_n\left[1+\frac{x_n-y_n}{2(1+a)}\right].
$$

Use $1+t\leq e^t$, $x_n-y_n\leq x_n$, and
$1-\rho=(1-a)/(1+a)^2$. Since $\Delta_0=1$,

$$
\Delta_n\leq r^n
\exp\left(\frac{1-a}{2(1+a)(1-\rho)}\right)
=r^ne^{(1+a)/2}<3r^n.
$$

For $0\leq B_n\leq A_n\leq1$,
$A_n^k-B_n^k\leq k\Delta_n$. Equation (5) proves (8).
The factors $1+(x_n-y_n)/(2(1+a))$ are strictly positive and their
deviations from 1 are absolutely summable. Their infinite product
$\kappa(s)$ is therefore finite and positive: for all sufficiently large
indices the absolute deviation is at most $1/2$, and
$|\log(1+t)|\leq2|t|$ proves convergence of the tail logarithms. Hence
$\Delta_n/r^n\to\kappa(s)>0$. The polynomial difference quotient
$[H(A_n)-H(B_n)]/\Delta_n\to H'(a)>0$, proving the asymptotic.

**Extinction, including the equality case.** If $\lambda<1$, Markov's
inequality and (8) give (9). To cover equality without silently applying a
finite-type spectral theorem, let $\delta=1-us>0$. Every positive cell has
probability at least $\delta$ of leaving no children. Write
$t_n=\Pr(Z_n^+>0)$. For every integer $L\geq1$, independence implies

$$
t_n-t_{n+1}\geq\delta^L\Pr(1\leq Z_n^+\leq L).
$$

Summation bounds the sum of the right-hand probabilities by $\delta^{-L}$.
The first Borel–Cantelli lemma implies that each finite positive population-size
range is visited only finitely often. On eternal survival, therefore,
$Z_n^+\to\infty$. If that event had positive probability, Fatou's lemma would
give $\liminf_n\mathbb EZ_n^+=\infty$, contradicting the uniform bound
$3S_0$ when $\lambda\leq1$. This proves critical extinction as well.

**Persistence.** Define the positive mean operator
$(Mf)(k)=us^k\sum_jB_{kj}[f(j)+f(2k-j)]$, with $f(0)=0$.
Differentiation of the binomial generating function and $T(a)=a$ give

$$
(Mh)(k)=2us^k k\left(\frac{1+a}{2}\right)^{2k-1}
=\frac{4ua}{1+a}ka^{k-1}=\lambda h(k).
$$

The correlated sister product is exactly

$$
Q_k=us^k\mathbb E[h(J)h(2k-J)]
=us^k a^{2k-2}\frac{k(2k-1)}2.\tag{11}
$$

Indeed $\mathbb E[J(2k-J)]=k(2k-1)/2$. Also
$h(k)\leq\sum_{j\geq1}ja^{j-1}=H_*$. With $x=sa<1$,

$$
\frac{Q_k}{h(k)}=\frac{u(2k-1)x^k}{2a}
\leq\frac{u}{2a}\sum_{j\geq1}(2j-1)x^j
=\frac{us(1+x)}{2(1-x)^2}\leq C_*.
$$

The vector $v_k=1-\varepsilon h(k)$, including $v_0=1$, lies in
$[0,1]$. Substitution into the extinction operator $\mathcal F$ of (1)
gives

$$
\mathcal F(v)_k=1-\varepsilon\lambda h(k)+\varepsilon^2Q_k
\leq1-\varepsilon h(k)=v_k.
$$

Monotonicity and iteration from the least starting vector imply
$q_k\leq v_k<1$ for every positive $k$. Independent founders give (10).
This proves genuine global survival rather than merely growth of an expected
value. Finally solve $2u(1-\sqrt{1-s})>1$, preserving the condition
$u>1/2$ before squaring, to obtain the displayed threshold. ∎

### Theorem 4 — ecDNA-free dominance is compatible with positive survival

Assume 1–7 and $N\geq1$. The pooled zero fraction satisfies
$f_n^0\to1$. If $u>1/2$, then on survival of the **total** genealogical
population,

$$
\frac{Z_n(0)}{N_n}\longrightarrow1\quad\text{almost surely}.\tag{12}
$$

If $u\leq1/2$, the total population dies out almost surely; a post-extinction
fraction is undefined. In particular, when Theorem 3 gives positive ecDNA
survival, (12) still holds on that event. A disappearing proportion is not
proof of eradication.

**Proof.** $A_n,B_n\to a>0$ and $H(a)>0$ prove the pooled statement.
For $u>1/2$, put $m=2u$. By (8), for every $\epsilon>0$,

$$
\sum_n\Pr(Z_n^+/m^n>\epsilon)
\leq3S_0\epsilon^{-1}\sum_nr^n<\infty.
$$

Thus $Z_n^+/m^n\to0$ almost surely. A zero-founder clone is the ordinary
branching process with zero children with probability $1-u$ and two with
probability $u$. If $Y_n$ is its size, $W_n=Y_n/m^n$ is a martingale of
mean 1. With offspring variance $\sigma^2=4u(1-u)$, conditional variance
and summation give

$$
\mathbb EW_n^2=1+\frac{\sigma^2}{m(m-1)}(1-m^{-n}).
$$

By the standard $L^2$-bounded martingale convergence theorem, it converges
in $L^2$ and almost surely to $W\geq0$, with
$\mathbb EW=1$. The branching decomposition implies
$\Pr(W=0)=1-u+u\Pr(W=0)^2$. Its root cannot be 1 because the mean is
positive. Hence $\Pr(W=0)=(1-u)/u$, which also equals the zero clone's extinction
probability by the same least-fixed-point argument applied to the scalar
offspring generating function $f(z)=1-u+uz^2$. Consequently,
$W>0$ on zero-clone survival.

It remains to prove that total survival yields such a zero clone. Choose a
finite integer $K\geq1$ with $2us^{K+1}<1$. A tree restricted to types $>K$,
discarding every daughter of type at most $K$, including zero, has mean offspring at
most $2us^{K+1}<1$ and therefore dies out by a first-moment bound. A surviving
positive tree must have infinitely many nodes of types $1,\ldots,K$:
otherwise the finitely many remaining high-type roots would each generate
finite pruned trees almost surely. Explore the positive tree breadth-first.
At each newly encountered low-type node, there is conditional probability at
least

$$
c=us^K2^{1-2K}\frac{2u-1}{u}>0
$$

of division, an extreme split producing a zero daughter, and survival of that
daughter's zero clone. Distinct newly founded zero clones are disjoint. Conditional on their
founding events, their future offspring variables are independent of all
positive-tree offspring variables. Reveal each such clone with its own trial;
the probability that the first $l$ trials exist and all fail is at most
$(1-c)^l$. No conditioning on future positive survival is used. Infinitely
many trials therefore produce a surviving zero clone almost surely. If the
positive tree dies but the total tree survives, one of the finitely many zero
clones already survives. Initial zero founders are handled in the same way.

A surviving zero clone founded at generation $g$ gives
$\liminf N_n/m^n\geq m^{-g}W>0$. Combine this with
$Z_n^+/m^n\to0$ to prove (12). For $u\leq1/2$, dominate the total population
by the zero-or-two process with division probability $u$; its least fixed
point is 1, including $u=1/2$. ∎

### Theorem 5 — exact Fano collapse, universal limit, and selection response

Assume 1–7 and at least one positive founder. Let $F_n=V_n/\mu_n$.

**(a) Monomorphic founders.** If every founder has the same positive copy count
$k$, then for every $n\geq0$,

$$
F_n=1+E_n/D_n-D_n/A_n.\tag{13}
$$

This is independent of $k$, founder number, and $u$, with
$F_0=0$, $F_1=1/2$, and $F_2=(3+s)/(2(1+s))$.

**(b) Arbitrary finite mixtures.** For every permitted founder multiset,

$$
\mu_n\to0,\qquad
F_n\to C(s)=1+\sum_{j=0}^{\infty}\frac{D_j}{1+A_j},\tag{14}
$$

independently of the initial mixture and $u$. The positive series has the
explicit truncation certificate

$$
0\leq C(s)-\left[1+\sum_{j=0}^{l-1}\frac{D_j}{1+A_j}\right]
\leq\frac{s^l}{1-s}.\tag{15}
$$

In particular $C(s)>3/2$. The verifier encloses the finite series using
outward-rounded exact rationals before adding (15); ordinary rounded decimals
are not treated as certificates.

**(c) One-generation selection response.** For any finite-support mother
distribution $P$, define $P_s(k)=s^kP(k)/\sum_j s^jP(j)$, with mean
$\mu_s$ and variance $V_s$. The pooled daughter distribution obeys

$$
\mu_{\mathrm d}=\mu_s,\qquad
V_{\mathrm d}=V_s+\mu_s/2,\qquad
\frac{d\mu_{\mathrm d}}{d\theta}
=\mu_{\mathrm d}/2-V_{\mathrm d}\leq0,
\quad s=e^{-\theta}.\tag{16}
$$

The derivative compares the same initial mother distribution at different
parameter values. It is not a derivative along an adapting population's time
trajectory.

**Proof.** Divide (6)'s variance by its mean to get

$$
F_n=1+E_n/D_n+D_n\left[
\frac{H''(A_n)}{H'(A_n)}-\frac{H'(A_n)}{H(A_n)}\right].\tag{17}
$$

For $H(z)=Nz^k$, the bracket is $-1/A_n$, giving (13). Substitution of
(7) gives the stated first three values. All denominators are positive when
there is a positive founder. Divide the last equation of (7) by the second:

$$
E_{n+1}/D_{n+1}=E_n/D_n+D_n/(1+A_n).
$$

Since $0<A_n\leq1$, (7) implies $0<D_n\leq s^n$. The series converges
and its remaining terms are bounded by a geometric series, proving (15).
Also $A_n\to a>0$, while $D_n\to0$; all polynomial ratios in (6) and
(17) stay bounded. This proves (14). The first series term is $1/2$ and
the next is strictly positive, proving $C(s)>3/2$.

For (16), selection weights mothers by $us^k$; the constant $u$ cancels
in their conditional law. Each selected mother contributes two daughters.
For a randomly pooled daughter, $\mathbb E[J\mid k]=k$ and
$\operatorname{Var}(J\mid k)=k/2$. Total expectation and total variance give
the first two identities. Differentiate the finite exponential-tilt sum to
obtain $d\mu_s/d\theta=-V_s$, which gives the last identity. ∎

### Theorem 6 — comparison certificates and a linear-hazard bracket

Assume 1–5 and 7 with arbitrary $p_k$. Comparison benchmarks use
$0<u\leq1$, $0<s<1$, as in Theorem 3. The conditions below concern all
positive integers $k$, not just observed types.

* If $p_k\leq u s^k$ for every $k\geq1$, and the benchmark has
  $\lambda(u,s)\leq1$, eradication is certain.
* If $p_k\geq u s^k$ for every $k\geq1$, and $\lambda(u,s)>1$, persistence
  has positive probability for every finite initial state containing a positive
  cell; the lower bound (10) still applies.
* For constant $b_k=b>0$ and $d_k=\gamma k$, $\gamma\geq0$, and an
  initial state containing a positive cell, persistence holds when $\gamma<b/3$, and eradication is certain when
  $\gamma\geq b/2$. No assertion about the gap is made.

**Proof.** Label a potential binary tree in advance, give each node an
independent uniform event variable and binomial partition randomness, and
couple processes at corresponding nodes of the same type. If the lower
division-probability tree divides, the upper tree also divides with the same
daughter split. Induction embeds the lower positive genealogy in the upper.
Theorem 3 then proves both comparisons. Zero-cell behavior cannot affect the
positive tree because reacquisition is prohibited.

For linear killing, write $c=\gamma/b$, so $p_k=(1+ck)^{-1}$. When $c>0$,
the binomial theorem gives $(1+c)^k\geq1+ck$, hence
$p_k\geq[1/(1+c)]^k$. The benchmark with $u=1$, $s=1/(1+c)$ is
supercritical when $c<1/3$. At $c=0$, no cell dies and at least one
positive daughter always remains, so persistence is certain.

For $c\geq1/2$, the mean number of positive daughters from any positive
mother is

$$
\frac{2(1-4^{-k})}{1+ck}
\leq\frac{4(1-4^{-k})}{k+2}\leq1.
$$

The last inequality is equivalent to $2-k\leq4^{1-k}$: it is equality for
$k=1$, and true for every $k\geq2$. Thus $\mathbb EZ_n^+$ is uniformly
bounded by its initial positive count. Every positive cell dies with
probability at least $c/(1+c)>0$. The finite-size visitation and Fatou
argument in Theorem 3 proves extinction, including $c=1/2$. ∎

## 4. Boundary cases and what is not solved

The exact finite-generation formulas also have natural degenerate extensions.
The continuous-time exponential-hazard parameterization itself uses
$u,s>0$.

* An empty population, or an initial population with no positive cells, has
  positive-lineage eradication probability 1 at generation zero. An empty
  population has no copy-number distribution.
* At $s=0$, all positive founders die at their first event in the discrete
  extension. At $u=0$, every cell dies at its first event. The convention
  $0^0=1$ is used only to define the zero type's probability $p_0=u$.
* At $s=1$, $u>0$, and with a nonempty initial state, the pooled mean stays at the initial mean and
  $V_n=V_0+n\mu_0/2$. This follows directly from $A_n=D_n=1$,
  $E_n=n/2$, and reproduces the familiar neutral mechanism. The finite Fano
  limit theorem excludes this boundary. For $u\leq1/2$, total extinction is
  certain. For $u>1/2$, positive survival follows by comparison with some
  $s'<1$ sufficiently close to 1. If $u=s=1$, positive survival is certain
  because a positive division cannot produce two zero daughters.
* At $\lambda=1$, extinction is proved, but a critical survival-tail
  asymptotic is **not** claimed. The explicit exponential error bound only
  decays when $\lambda<1$.
* A zero fraction is undefined after total extinction. Theorem 4 states its
  almost-sure limit only on total survival; (6)'s ratio of expectations is a
  different observable.
* General $p_k$, linear killing within the bracket's gap, calendar-time
  moments/eradication-time distributions, density dependence, multiple ecDNA
  species, therapy schedules, and renewed ecDNA formation remain open in this
  contribution. No numerical pattern in these settings is promoted to a theorem.

## 5. Quantitative, falsifiable experimental predictions

These are conditional tests of a nonclinical model, not proposed patient
interventions. A suitable research design would combine lineage tracking,
copy counts, and cell fates in independent cell-culture microcolonies.
Record each cell when it is born into the specified genealogical generation,
including cells that subsequently die. A later calendar-time census of only
surviving generation-$n$ cells introduces additional survival bias.

**First test the assumptions.** Estimate division-before-death probability by
copy count, distinguishing right-censored observations from death. Check whether
$\log p_k=\log u+k\log s$ is credible over the measured range. Check exact
doubling and the sister constraints $J_1+J_2=2k$,
$\operatorname{Var}(J_1\mid k)=k/2$, and
$\operatorname{Cov}(J_1,J_2\mid k)=-k/2$. Retention loss, clustering, or
stress-dependent segregation can falsify these assumptions before testing the
population predictions. Finite measurements cannot establish a global bound
over unobserved high copy numbers for Theorem 6.

1. **Copy-number-independent Fano time course.** Start replicate sets at two or
   more narrow, independently measured founder copy counts. Pool all cells at
   the same lineage generation, including zeros, and compare the full curve
   (13). For genuinely monomorphic founders the predicted $F_1=0.5$ and
   $F_2=(3+s)/(2(1+s))$ do not depend on initial copy count or founder number.
   Different curves after accounting for initial mixing and measurement noise
   reject the model. Treat whole microcolonies as statistical replicates.
2. **Two-generation parameter inference and boundary.** A measured
   $1<F_2<3/2$ identifies $s=(3-2F_2)/(2F_2-1)$ within this exact model.
   It cannot identify $u$: normalized copy distributions cancel $u$.
   Independently measure baseline fates to estimate $u$, then predict the
   sign of $2u(1-\sqrt{1-s})-1$. For the ideal $u=1$ specialization,
   $F_2\geq15/14$ predicts eventual positive eradication and
   $F_2<15/14$ predicts a nonzero persistence probability. The threshold is
   asymptotic; a surviving finite-duration well does not disprove certain
   eventual extinction, and an extinct well does not disprove positive
   survival probability.
3. **Universal mixed-population Fano plateau.** Under sustained model
   parameters, all-cell pooled $V_n/\mu_n$ approaches the same $C(s)$ for
   distinct initial mixtures. At $s=3/4$, the certified interval routine
   supports $C(s)\approx2.67933549389746$. A stable, reproducible discrepancy
   rejects the stipulated dynamics. Both numerator and denominator eventually
   become small in zero-rich pools, so adequate replicate counts and explicit
   detection-error modeling are necessary; positive-only gating tests a
   different formula.
4. **Selection-response identity.** For a fixed mother distribution across
   conditions, compare finite differences in mean daughter copy count against
   $d\mu_{\mathrm d}/d\theta=\mu_{\mathrm d}/2-V_{\mathrm d}$. Here
   $\theta=-\log s$ must be estimated, not substituted with drug concentration
   without a calibration model. For example, a reliably negative
   $V_{\mathrm d}-\mu_{\mathrm d}/2$ contradicts the assumptions.

For an additional population check, independent replicate wells should give
$\mathbb EZ_{n+1}^+/\mathbb EZ_n^+\to\lambda$, even though their all-cell
zero fraction tends to one. Copy counts alone without fates or generations do
not establish eradication.

## 6. Verification and falsification record

There are three separate levels of evidence: symbolic proofs above, exact
finite arithmetic, and a simulator that implements the offspring experiment
without importing the verifier. Finite checks and simulations can expose
errors; neither proves an infinite-state theorem. The main results are general
theorems and a proven scaling limit, not small-population counterexamples.

| Challenge | Resolution or retained limitation |
| --- | --- |
| “The basic model and neutral moments are new.” | Rejected after Turner/Lange and historical theory were inspected. The target was explicitly narrowed before implementation. |
| Independent binomial draws for both sisters. | Rejected: their sum must be $2k$. Tests enumerate complete population outcomes and check the cross terms in (1), (3), and (11). |
| Inferring persistence from a growing mean. | Rejected. Theorem 3 constructs a genuine extinction supersolution, giving a strictly positive survival lower bound. |
| Inferring critical extinction from a bounded mean alone. | Incomplete argument rejected. The proof additionally establishes escape from every finite positive population-size range on survival. |
| Applying finite-type Perron–Frobenius at an arbitrary type cap. | Avoided. Critical extinction and persistence are proved directly for all positive integers. Finite support in code comes only from a stated finite horizon. |
| Calling the parameterized death hazard linear. | Rejected. It is exponential; a separate proven comparison bracket handles a linear hazard. |
| Calling a zero-rich population eradicated. | Rejected. Theorem 4 proves zero dominance can coexist with positive persistence. |
| Treating pooled moments as the expectation of a within-well ratio. | Rejected. Definitions, formulas, and Monte Carlo estimators use explicit denominators. Extinct wells remain in all population-level statistics. |
| Inferring infinite-time survival from finite simulations. | Rejected. Simulation outputs always report their finite generation; threshold simulations are compatibility checks only. |
| Extending the finite Fano plateau to no therapy. | Rejected. $s=1$ has neutral linear variance growth and is treated separately. |
| Numerical rounding at criticality or in a long series. | Threshold classification uses rational inequalities, including exact equality. Fano intervals use rational outward rounding plus an analytic remainder. Decimal illustrations are labelled non-certified. |
| Assuming binomial segregation remains correct under treatment. | Challenged by recent primary work. This is an explicit experimental failure mode, not a robustness claim. |
| Declaring absolute novelty or clinical relevance from searches and AI review. | Rejected. Priority remains provisional and biology unvalidated. |

The test file independently enumerates complete small-generation population
laws as structural checks, checks arbitrary-copy algebraic identities including
large types, cross-checks analytic derivatives against a forward expected-count
algorithm, tests critical and degenerate boundaries, verifies rational Fano
enclosures, and compares independent fixed-seed trajectories with exact
finite-generation extinction and moment predictions. Monte Carlo tolerances
use independent replicate populations, not related cells. Resource caps abort
with an error; simulations never silently count a censored well as extinct.

### Recorded finite numerical checks

The standalone verifier completed **755 exact checks**. The new test file
completed **17 tests in about two seconds** in the development workspace,
well below its 60-second budget. Tests are deterministic and require no network
or non-stdlib package. The repository's schema validation and complete Python
and Node suites were also run; PR checks provide the final integration gate.

The following synthetic results are reproducible with
`simulation.py --trials 10000`. Every SE uses independent replicate populations.
Rounded displays below do not replace the exact rationals used by the tests.

| Case and observable | Theorem value | Simulation | SE | Seed |
| --- | ---: | ---: | ---: | ---: |
| Founders (0,1,3), n=4, u=0.9, s=0.8: expected positive count | 1.94156202 | 1.92470000 | 0.02455938 | 2718281828 |
| Same case: extinction by generation 4 | 0.48835267 | 0.48990000 | 0.00499923 | 2718281828 |
| Same case: pooled Fano factor | 2.68994982 | 2.67124081 | 0.01987937 | 2718281828 |
| Founder k=1, n=4, u=1, s=0.9: pooled Fano | 2.03742466 | 2.03524582 | 0.01134007 | 2718281859 |
| Founder k=8, otherwise same: pooled Fano | 2.03742466 | 1.97068377 | 0.03695921 | 2718281866 |
| Founders (1,3), n=12, u=1, s=0.75: pooled Fano | 2.67684985 | 2.68213874 | 0.01961729 | 2718281830 |

At the critical parameter s=0.75, u=1, 10,000 simulated one-copy founders
still had a **0.1053 survival frequency at generation 25** (seed 2718281839).
This is compatible with eventual extinction and explicitly prevents equating
finite survival with persistence. The supercritical s=0.8 and subcritical
s=0.7 cases at the same horizon had survival frequencies 0.3054 and 0.0206,
respectively; these illustrate the theorem rather than prove its boundary.

Self-review was performed before opening the PR. Independent AI critiques
checked every displayed proof, with special attention to the critical case,
supersolution direction, and random-versus-pooled fractions. Corrections included
a stopped-process nonexplosion proof, explicit empty/zero-founder and neutral
boundaries, conditional independence of newly founded zero clones, and the
Bansaye DOI. All six stated theorems retained complete proofs; unresolved generalizations
remain explicitly out of scope. This is computational and mathematical self-review,
not independent scientific or clinical validation.

## 7. Reproduce

From the repository root, using only the Python standard library:

```bash
python research/ecdna-segregation-therapy/verifier.py
python research/ecdna-segregation-therapy/verifier.py --illustrate
python research/ecdna-segregation-therapy/simulation.py --trials 10000
python -m unittest discover -s tests -p 'test_research_ecdna_segregation_therapy.py' -v
```

`verifier.py` uses `fractions.Fraction` for exact checks and certified rational
intervals. Its optional `decimal.Decimal` illustration is explicitly not an
interval proof. Exact unrounded iterates have exponentially growing integer
bit lengths, so use the interval routine for long-horizon series rather than
requesting arbitrarily many exact polynomial iterations. `simulation.py` uses
`random.Random` with recorded fixed seeds and direct complementary sister
sampling. Division probabilities and final summaries are floating-point
approximations; extreme probabilities can underflow, so extreme-parameter
claims rest on the exact verifier and proofs. It does not import, call, or fit itself to the verifier.

The repository's required validation commands also apply:

```bash
python scripts/validate.py
python -m unittest discover -s tests -v
node --test tests/*.test.*
```

Only this new research directory and
`tests/test_research_ecdna_segregation_therapy.py` are part of this contribution.
No catalog record, publication registry entry, existing research thread,
workflow, dependency, or top-level documentation is changed. This work remains
nonclinical computational mathematics with unassessed biological hypotheses.
