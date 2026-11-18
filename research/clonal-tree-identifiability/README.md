# How many bulk mixtures can identify a clonal tree?

**Nonclinical computational mathematics. All numerical inputs are synthetic
rational numbers; no patient data are used.** This is an AI-assisted research
note with complete mathematical arguments and independent computational checks,
not an externally reviewed scientific publication or a clinical sampling
recommendation. Literature search cutoff: **2026-11-17**. Drafting and separate
assistant criticism used GPT-6 Astra Pro; executable mathematical evidence is
provided by the two independent standard-library implementations below.

The answer depends on the quantifiers. In the single-founder convention defined
below, every fixed tree can be uniquely identified by **one deliberately designed
sample if it is a path, and two otherwise**. In contrast, for three or more
clones, **no finite number of arbitrary samples guarantees uniqueness**, even
after imposing distinct frequencies, strictly positive clone usages, and maximal
matrix rank. Ambiguity occupies an open set, so removing a measure-zero set of
degenerate inputs does not repair that guarantee.

The unverified candidate mathematical refinement is the universal, rational,
strictly positive two-sample construction, its proof of uniqueness against
*all* alternative trees, and its explicit perturbation radius. Several ingredients and supporting bounds
are already known. The literature review below explains that distinction before
the new arguments.

## 1. Literature first: what is established, and the precise pivot

| Primary source | Relevant established result | Relation to this note |
| --- | --- | --- |
| Jiao et al. (2014), [PhyloSub][jiao], section “Topological constraint rules” | A single frequency vector always admits a chain. The sum and crossing rules can resolve a three-node relationship; two arbitrary samples need not resolve every relationship. | The one-sample obstruction and crossing rule are prior work. Theorem 1 supplies a simultaneous construction for every finite tree. |
| El-Kebir et al. (2015), [AncesTree][ances], Theorems 1–2 | Characterization by ancestry-graph arborescences satisfying the sum condition; NP-completeness of the frequency factorization problem. The paper's VAF convention has an additional factor of one half. | We use cellular prevalence, so the factorization is exactly `F = U B`. Lemma 1 is included for self-containment and is not new. |
| Malikic et al. (2015), [CITUP][citup] | Joint reconstruction across samples by optimization over candidate phylogenies. | An inference method, not an extremal sample-count theorem in the inspected material. |
| El-Kebir et al. (2016), [SPRUCE][spruce]; [full preprint][spruce-preprint], Theorem 3 | Multi-state mixture reconstruction; hardness already with two states and two samples. | Constructing identifying data *from a known target tree* does not solve inference for arbitrary two-sample data. |
| Satas and Raphael (2017), [PASTRI][pastri] | Tree-specific feasible frequency regions, importance sampling, and enumeration of compatible trees. | Feasible-region geometry and enumeration are not new here. |
| Pradhan and El-Kebir (2018), [Theorem 2][pradhan]; Qi, Pradhan and El-Kebir (2019), [Theorems 10–11][qi] | A product-of-indegrees upper bound for distinct columns; the 2019 extension proves counting is #P-complete even for two samples and studies approximation/sampling hardness. | The factorial upper bound in Theorem 3 is an elementary extremal consequence of their bound. Our attainment witnesses and boundary distinctions should not be mistaken for a new counting-hardness result. |
| Sundermann et al. (2021), [SubMARine][submarine] | Partial clone trees and sound inference of common ancestral relations, with experimental unique-solution recovery. | The inspected claims do not establish a universal minimum number of designed samples. |
| Marangola, Sheikholeslami and Bento (2026), [published journal article](https://doi.org/10.1109/TCBBIO.2026.3726506) and [arXiv v2][degeneracy], Theorems 1–3 | Ensemble results on longitudinal constraints and degeneracy, and bounds for dynamic constraints in a restricted family of alternative trees. | Essential recent prior art. This note makes no claim to be the first theoretical or ensemble analysis of ambiguity. |
| Brualdi, Jung and Trotter (1994), [Lemma 2.3][posets], attributing the result to Wolk | Rooted-tree orders have order dimension at most two. | Encoding an ancestor order with two opposite orderings is classical. Additive feasibility and exclusion of skipped parents are the additional issues addressed here. |

**Pivot.** A finite worst-case threshold for unrestricted observations is not the
right open target: simple ambiguous inputs rule it out, and the usual
distinct-column maximum follows from an existing bound. We therefore ask the
closest constructive question:

> For each fixed labelled tree, what is the least positive number of mixtures
> for which some strictly positive usage matrix makes that tree the unique
> explanation among every tree permitted by the model?

We did not find the all-tree constructive upper bound, the explicit robustness
certificate, or the exact one-row criterion below in the inspected primary
sources. This is a **bounded novelty assessment, not proof of bibliographic
priority**. The lower bound, the crossing rule, the two-order ingredient, and the
ancestry-graph counting bound are explicitly credited above. Theorems 1–2 and
Lemma 2 are candidate refinements/syntheses; Theorem 3 is supporting clarification
with a modest novelty claim.

Searches included `"perfect phylogeny mixture" "two samples"`, `"tumor
phylogeny" "sample complexity"`, `"clonal" "two samples suffice"`, `"sum
condition" "unique" tumor tree theorem`, `"phylogeny" "superincreasing"`,
`"tumor phylogeny" "order dimension"`, and combinations of identifiability,
minimum samples, uniqueness, and the named authors. Both broad and targeted
searches were used, with direct inspection of the cited model/theorem sections.
The 2026 paper was inspected in its October 2, 2026 revision, not just its 2025
first version. Absence from these searches is not evidence of absence from all
literature.

## 2. Model, conventions, and quantifiers

All main results use the following assumptions.

1. There are exactly `n >= 1` known, distinct, labelled mutation clusters and
   corresponding malignant clone types. Each clone introduces its own cluster.
   No cluster is missing, merged with another, or inferred from noisy reads.
2. A candidate `T` is any rooted tree on these `n` labels. There is exactly one
   malignant founder, whose label is not supplied in advance. Its mutation is
   inherited by all malignant clones. The normal population is **not** one of
   these `n` vertices. Section 7 treats the other common root convention.
3. Infinite sites / perfect phylogeny holds: each cluster arises once, is
   inherited, and is never lost. Frequencies are exact cellular prevalences,
   without sequencing, copy-number, or purity-estimation error.
4. Samples are unrelated mixtures of these same clone types. There are no
   temporal, spatial, evolutionary-rate, lower-detection-limit, or distributional
   restrictions on usages. In particular, a designed mixture may depend on the
   entire target tree.
5. `B_T[i,j] = 1` exactly when vertex `j` is an ancestor of vertex `i`, including
   equality. The observed `m` by `n` matrix obeys `F = U B_T`, with
   `U[s,i] >= 0` and `sum_i U[s,i] <= 1`. The unused fraction is normal tissue.
   Candidate explanations may have zero usages, even when our constructed true
   usages are strictly positive.
6. Uniqueness means equality of the complete rooted **labelled** tree. Labels
   are not quotiented out by tree isomorphism. Since `B_T` is invertible, each
   candidate tree has at most one usage matrix.

Write `C_T(v)` for the children of `v`, `D_T(v)` for its inclusive descendant
set, and `S(F)` for the set of consistent trees. A *path* here means a rooted
directed chain starting at an endpoint; a tree with a root in the middle of an
undirected path is branching. Mathematical labels below are `1,...,n`; the code
uses `0,...,n-1` and a parent value of `-1` for the malignant root.

Define `d(T)` to be the least **positive** integer `m` for which there exists a
rational `F = U B_T`, with every `U[s,i] > 0` and every row sum `< 1`, satisfying
`S(F) = {T}`. This is an existence/design number, not a guarantee about the first
`m` biopsies. Requiring positive `m` avoids the convention that the only tree on
one vertex can be identified without observing anything.

“Generic” must also be specified. Distinct columns alone are a finite set of
inequalities, not a definition of generic identifiability. Our impossibility
statement concerns the stronger usual question of uniqueness outside a
Lebesgue-null exceptional set in the open usage domain, for a fixed true tree.
An open ambiguity set disproves that property.

### Lemma 1 (standard sum-condition equivalence)

Under assumptions 1–6, a tree `T` is consistent with `F in [0,1]^(m*n)` if and
only if, for every sample and vertex,

$$u_{s,v}=f_{s,v}-\sum_{c\in C_T(v)}f_{s,c}\ \geq\ 0. \tag{1}$$

These are the unique usages, and their row sum equals the founder frequency.

**Proof.** The factorization says
`f[s,v] = sum_{i in D_T(v)} u[s,i]`. The child descendant sets partition
`D_T(v) \ {v}`, giving (1). Conversely, sum (1) over `D_T(v)`; every term except
`f[s,v]` cancels, recovering the factorization. Summing over the whole tree gives
the founder frequency, which is at most one. Finally, a topological labelling
makes `B_T` triangular with diagonal entries one, proving invertibility. This is
the standard result of [El-Kebir et al.][ances].

### Lemma 2 (an exact criterion when the dominance order is known)

Suppose `T` is consistent and the directed graph with edge `a -> v` for distinct
vertices satisfying `f[s,a] >= f[s,v]` in **every** sample is exactly the strict
ancestor relation of `T`. Then `S(F) = {T}` if and only if

$$\text{for every proper ancestor }a\text{ of }v\text{ that is not its parent,}
\quad \exists s:\ f_{s,v}>u_{s,a}. \tag{2}$$

**Proof.** If (2) fails, move the entire subtree rooted at `v` from its parent to
`a`. The usage at `a` decreases by `f[s,v]` and remains nonnegative; the old
parent's usage increases by that amount. No other usage changes. This gives a
different consistent tree with the same root.

Conversely, every edge of a consistent alternative `T'` belongs to the specified
dominance graph. Thus the root is unchanged, and a changed parent can only be a
strict ancestor of the original parent. Among the vertices receiving at least
one new child, choose `a` of minimum depth in `T`. No original child of `a` can
leave: its new parent would be a receiver at smaller depth. Consequently `a`
retains all original children and receives an additional nonchild descendant
`v`. Consistency of `T'` implies `u[s,a] >= f[s,v]` in every sample, contradicting
(2). The dominance-graph hypothesis is essential; Section 9 gives a counterexample
without it.

## 3. Theorem 1: optimal designed sample number, with an explicit radius

**Theorem 1.** Under assumptions 1–6, for every rooted labelled tree `T`,

$$d(T)=\begin{cases}1,&T\text{ is a rooted path},\\2,&T\text{ has a branching vertex}.
\end{cases} \tag{3}$$

In particular, `max_T d(T)` is one for `n = 1,2` and two for `n >= 3`.
The construction below uses rational numbers with `O(n + log n)` bits per entry,
strictly positive usages, and normal fraction one half. If its denominator is
`D`, every real matrix `F'` with

$$\lVert F'-F\rVert_\infty\leq\rho,\qquad
\rho=\frac{1}{4nD}, \tag{4}$$

where the norm is the entrywise maximum `max_(s,v) |F'[s,v]-F[s,v]|`, also
lies in `(0,1)^(m*n)`, has strictly positive usages with row sums below one
for `T`, and has `S(F') = {T}`. Thus these identifying designs have a full
dimensional neighbourhood.

### Construction

Order each vertex's children in any fixed way (the verifier uses increasing
labels), and list the `L` leaves in depth-first order as
`ell_0,...,ell_(L-1)`. Put

$$W=n+1,\qquad S=W(2^L-1)+(n-L),\qquad D=2S.$$

Give each internal vertex unnormalized usage one in both rows. Give leaf
`ell_i` unnormalized usages `W 2^i` and `W 2^(L-1-i)` respectively. Divide each
usage by `D` and obtain `F` by summing usages over descendant sets. Both total
usages equal `S/D = 1/2`. If `L = 1`, the rows are identical; retain just one.

### Complete proof

**Feasibility and positivity.** Every assigned usage is positive and the total
is one half, so Lemma 1 gives a valid factorization. Each subtree contains a
leaf, hence each frequency is at least `W/D`. Every internal vertex has residual
usage exactly `1/D`, and the founder frequency is `1/2`.

**Exactly the intended dominance graph.** If `a` is a proper ancestor of `v`, its
subtree contains that of `v` and includes at least the positive usage at `a`.
Therefore `f[s,a] - f[s,v] >= 1/D` in both rows.

If two vertices are incomparable, their descendant leaf sets are disjoint
intervals in the depth-first leaf order. Suppose all leaf indices of subtree
`A` precede those of subtree `B`, and the smallest index in `B` is `b`. In the
first row, the numerator of `A` is at most

$$W\sum_{i=0}^{b-1}2^i+(n-L)=W(2^b-1)+(n-L),$$

whereas that of `B` is at least `W 2^b`. Their difference is therefore at least
`W-(n-L)=L+1 >= 2`. Reversing leaf weights reverses this strict comparison in the
second row, by the same argument. No incomparable pair can be joined by an edge
of the dominance graph. Hence that graph is exactly the ancestor relation.

**Uniqueness of direct parents.** For every pair in (2), `a` is internal, so
`u[s,a] = 1/D`, while `f[s,v] >= W/D > 1/D`. Lemma 2 excludes every different
tree, including trees that preserve all pairwise dominance inequalities but
skip one or more parents. Equivalently, the shallowest receiver argument gives
a negative candidate usage of at most `-(W-1)/D = -n/D`.

**Paths and the lower bound.** A finite rooted tree has one leaf exactly when
it is a rooted path. For these trees the same proof works with the single
retained row: every pair is comparable and all forbidden parent skips violate
the sum condition. Thus `d(T) <= 1`. For any one-row input, sorting its frequencies
in nonincreasing order produces a consistent rooted chain: successive usages
are adjacent differences and the final usage is the smallest frequency. If
the target tree branches, this chain is a different explanation. Thus no
branching tree has `d(T) = 1`. This obstruction is already in [Jiao et al.][jiao].

**Perturbation bound.** A sum-condition expression contains at most `n` entries
with coefficients `+1` or `-1`; its change is at most
`n rho = 1/(4D)`. True residuals remain at least `3/(4D) > 0`. Every pairwise
frequency difference changes by at most `2 rho < 1/D`, preserving all comparisons
used above and hence the entire dominance graph. For an alternative with edges
in this graph, its original negative residual at a shallowest receiver has
magnitude at least `n/D`; perturbing it by at most `1/(4D)` leaves it negative.
Every other alternative still has a forbidden edge. Finally, original
frequencies are between `W/D` and `1/2`, and `rho < W/D`, `rho < 1/2`, so all
perturbed entries lie strictly between zero and one. The true root remains below
one, which also bounds the perturbed usage row sum. These arguments include
`n = 1`, where there are no competing trees.

**Size of the rationals.** Since `L <= n` and
`D = 2[(n+1)(2^L-1)+n-L]`, the numerators and denominator have
`O(n + log n)` binary digits. All subtree sums have the same denominator. This
is an exact representational bound; the numerical margin can still shrink
exponentially with the number of leaves.

### A five-clone worked example

Take edges `1 -> 2`, `1 -> 3`, `2 -> 4`, `2 -> 5`. The ordered leaves are
`4,5,3`, so `W = 6`, `S = 44`, and `D = 88`. The construction gives

$$U=\frac1{88}\begin{pmatrix}1&1&24&6&12\\1&1&6&24&12\end{pmatrix},
\qquad
F=\frac1{88}\begin{pmatrix}44&19&24&6&12\\44&37&6&24&12\end{pmatrix}.$$

For example, vertices 2 and 3 cross: `19 < 24` but `37 > 6`.
Attaching vertex 4 directly to vertex 1 would add at least `6/88` to its child
sum, although the founder's unused malignant usage is only `1/88`. The full
proof handles all simultaneous reattachments, not just this one. Equation (4)
gives the conservative radius `1/1760` for this example.

## 4. Theorem 2: an exact one-sample uniqueness test

**Theorem 2.** Let `m = 1`, `n >= 1`, and all frequencies lie in `[0,1]`.
For `n = 1` the tree is unique. For `n > 1`, repeated frequencies prevent
uniqueness. Otherwise relabel so that

$$1\geq f_1>f_2>\cdots>f_n\geq0.$$

The input has exactly one consistent tree if and only if

$$f_i<f_{i+1}+f_n\qquad(1\leq i\leq n-2). \tag{5}$$

That tree is the descending chain. For `n = 2` the list of conditions is empty;
two distinct frequencies always give a unique tree.

**Proof.** The descending chain always satisfies the sum condition. If two
frequencies tie, interchanging their labelled positions in a sorted chain gives
a different consistent chain, so uniqueness fails. Assume now all are distinct.
If `f_i >= f_(i+1)+f_n` for an index in (5), start with the descending chain and
reattach the terminal vertex `n` directly to `i`. Its old parent loses a child;
the new parent's child sum is `f_(i+1)+f_n <= f_i`. All other conditions are
unchanged. Equality therefore gives ambiguity, too.

Conversely, any alternative to the sorted chain must branch, since every edge
strictly decreases frequency. Before the first branching vertex its initial
single-child trunk must be `1 -> 2 -> ... -> i`: the root must be the largest
frequency, and a sole child ancestors all remaining vertices, so must be their
largest. At the branch vertex `i`, vertex `i+1` must be a direct child; otherwise
the child on its path would have a smaller frequency while ancestoring it.
There is another child whose frequency is at least `f_n`. The sum condition
would imply `f_i >= f_(i+1)+f_n`, contradicting (5). This proves sufficiency.

No positivity is silently added here. For example, if `f_n = 0` and `n >= 3`,
strict ordering forces every relevant adjacent gap to be positive, so (5)
correctly rules out uniqueness.

## 5. Theorem 3: sharp maximum counts and the failed worst-case threshold

**Theorem 3.** Fix `n >= 1` and `m >= 1`, under assumptions 1–6.

1. Over all matrices in `[0,1]^(m*n)`,
   `max |S(F)| = n^(n-1)` (one when `n = 1`).
2. Over matrices whose `n` columns are pairwise distinct,
   `max |S(F)| = (n-1)!`. This maximum is attained with strictly positive
   usages for every consistent tree, distinct entries in every row, and
   `rank(F) = min(m,n)`. There is also a nonempty open set attaining the maximum.
3. If all feasible inputs are allowed, the least positive number guaranteeing
   uniqueness is one for `n = 1` and does not exist for `n >= 2`.
   If inputs must have pairwise distinct columns, it is one for `n = 1,2` and
   does not exist for `n >= 3`. For `n >= 3`, uniqueness is not an almost-everywhere
   property at any finite `m`, even conditional on any fixed true tree with
   positive usages and normal fraction.

**Proof of 1.** There are `n^(n-2)` unrooted labelled trees for `n >= 2`, by the
Pruefer-code bijection, and `n` choices of root. There cannot be more consistent
trees than this total. The all-zero frequency matrix has all-zero malignant
usages for every tree and attains the total. The one-vertex case is separate and
has one tree. This extremum includes unobserved malignant clones; it is a
mathematical boundary case, not a realistic sequencing observation.

**Upper bound in 2.** Distinct columns make the dominance graph acyclic: a
directed cycle would imply coordinatewise equality of its columns. A feasible
tree supplies a source reaching every vertex; thus the source is unique. In a
topological ordering starting there, vertex `j` has at most `j-1` possible
parents. Each choice of one earlier parent connects recursively to the source,
so the total number of candidate trees is at most
`product_(j=2)^n (j-1) = (n-1)!`; imposing the sum condition only removes choices.
This is the extremal specialization of the bound in
[Pradhan and El-Kebir][pradhan] and [Qi et al.][qi].

**Attainment and maximal rank.** Put

$$D_0=2(2^n-1),\qquad b_j=\frac{2^{n-j}}{D_0},\qquad
\epsilon=\frac1{4nD_0}.$$

For sample `s <= n`, use `b + epsilon e_s`, where `e_s` is the `s`th unit vector.
For `s > n`, use `b + epsilon e_1/s`. This is exactly the zero-based code's
`epsilon/(p+1)` when `p >= n`. Retain the first `m` rows.
In every row, each coordinate strictly exceeds the sum of all later coordinates:
the base excess is `1/D_0`, and the perturbation can decrease it by at most
`epsilon < 1/D_0`. The first coordinate stays below one. Hence **every** choice
of a parent for `j` from `1,...,j-1` satisfies all sum conditions strictly.
There are exactly `(n-1)!` such choices and no other possible edges.

Let `r = min(m,n)`. The leading `r` by `r` block of the first `r` rows is
`epsilon I_r + 1_r (b_1,...,b_r)`. Its determinant is

$$\epsilon^{r-1}\left(\epsilon+\sum_{j=1}^r b_j\right)>0.$$

For completeness, expand the determinant multilinearly in the columns. The
term choosing only identity columns is `epsilon^r`; terms choosing exactly one
multiple of `1_r` sum to `epsilon^(r-1) sum b_j`; any term choosing two such
columns vanishes because they are parallel. Thus the stated rank is exact.
All inequalities are strict, so a sufficiently small ball preserves all these
trees and the distinct order. The same nonzero minor also persists in a small
ball. This proves open-set attainment, including at maximal rank.

**Proof of 3 and the genericity statement.** The all-zero witness already
precludes the unrestricted guarantee for `n >= 2`. With distinct columns,
part 2 gives a maximum of one for `n <= 2` and at least two for `n >= 3`,
independently of `m`. For any fixed true tree, relabel a topological ordering as
`1,...,n`; the strict superincreasing construction includes that tree among all
`(n-1)!` alternatives. Because `F -> F B_T^(-1)` is an invertible linear map,
the open ambiguity ball maps to an open subset of the strictly positive usage
domain with row sums below one. Such a set has positive Lebesgue measure.
It cannot be removed by a measure-zero exceptional set, or by deleting any
finite collection of proper algebraic hypersurfaces. Under a sampling law with
positive density throughout that domain it also has positive probability.
This does not assert anything about distributions deliberately supported only
on an identifying region, or about the infinite-sample limit.

## 6. Exact exhaustive enumeration

The following table is reproduced by `verifier.py`. “Maximum, distinct columns”
counts every candidate tree for the explicit full-rank witnesses in Theorem 3;
the proof supplies the global upper bound. It is **not** a search over all
rational matrices. “Identifying designs checked” means every possible target
tree at that `n` was constructed and tested against every possible candidate.

| n | All rooted labelled trees / unrestricted maximum | Paths needing 1 sample | Branching trees needing 2 | Maximum, distinct columns, m=1 | m=2 | m=3 | Identifying designs checked |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 1 | 0 | 1 | 1 | 1 | 1 |
| 2 | 2 | 2 | 0 | 1 | 1 | 1 | 2 |
| 3 | 9 | 6 | 3 | 2 | 2 | 2 | 9 |
| 4 | 64 | 24 | 40 | 6 | 6 | 6 | 64 |
| 5 | 625 | 120 | 505 | 24 | 24 | 24 | 625 |

There are 701 target trees in total and 394,807 target/candidate pairs across
these five sizes. There are `n!` rooted paths: each permutation lists its labels
from root to terminal leaf once. The complement `n^(n-1)-n!` consists of the
branching trees. The observed patterns suggest (3) and the factorial bound;
the proofs establish them for every `n`, independently of the finite checks.

Two intentionally different implementations check the result:

- The verifier generates labelled trees through Pruefer sequences, then orients
  each tree at every possible root. It evaluates exact sum conditions using
  `fractions.Fraction`, or integer comparisons after clearing denominators.
- The test independently enumerates parent assignments and rejects cycles,
  constructs the ancestry matrices by walking ancestor chains, and inverts them
  by generic rational Gaussian elimination. It checks `U = F B^(-1)` and the
  usage row bounds, without using the verifier's sum-condition oracle.

Clearing rational denominators is exact integer arithmetic, not rounding. No
floating-point tolerance, random trial, NumPy, optimization solver, network
request, or biological dataset is needed. The dedicated test additionally checks:

- The one-row criterion on all 1,023 decreasing rows of lengths one through five
  chosen from `{1/12,2/12,...,11/12}`, plus zero, tie, and equality examples.
- Lemma 2 on all 4,096 two-row, three-column matrices with entries in
  `{1/5,2/5,3/5,4/5}`. Exactly 1,350 `(F,T)` pairs satisfy its hypotheses;
  1,152 are unique and 198 ambiguous. The criterion agrees in every case.
- The robustness certificate against every alternative for every target through
  four vertices, and the normal-root translation for all 145 fixed-root trees
  on two through five total vertices.
- The tied-positive path-partition counts and explicit falsification examples
  in Section 9.

Finite enumeration checks the implementations and examples; it is not a
substitute for the all-`n` proof.

## 7. The normal-root convention changes the counts

Some formulations instead use a fixed normal root `0`, with frequency one, and
`n` malignant vertices whose tree after removing `0` may be a forest. This is a
larger candidate class than assumption 2. The distinction cannot be hidden in a
change of notation.

For this convention the same arguments give:

| Quantity | Single malignant founder, n total malignant vertices | Fixed normal root plus n malignant vertices |
| --- | --- | --- |
| Number of all candidate trees | `n^(n-1)` | `(n+1)^(n-1)` |
| Maximum with distinct observed columns | `(n-1)!` | `n!` |
| First size with open-set worst-case ambiguity | `n = 3` | `n = 2` |
| Designed optimum | 1 for a rooted path, 2 otherwise | 1 for a path including root 0, 2 otherwise |

Here is a proof of the translation. Count labelled trees on `n+1` vertices with
one root fixed, instead of counting all roots on `n` vertices. Topologically
order the distinct malignant columns and prepend the fixed root `0`; edges into
that fixed root are prohibited, even if a malignant column equals its frequency
column of ones. Each malignant vertex now has the additional possible parent `0`, so
the indegree product is `1*2*...*n`. Superincreasing malignant frequencies whose
total is below one attain this bound. For the design theorem, apply the
construction to all `N = n+1` vertices, with normal root `0`, then multiply all
frequencies and usages by two and remove the normal frequency column. The root
frequency is now one, its normal usage is positive, and the same comparisons
and strict sum-condition violations exclude all alternatives. A radius
`2/(4ND)` for the remaining entries works with the root coordinate held fixed.
Every single sample still admits the chain beginning at `0`, proving the lower
bound. The main enumeration table uses the single-founder convention; the test
also checks this translation on small fixed-root trees.

## 8. How to reproduce

From the repository root, using Python 3.10 or later and only the standard library:

```bash
python research/clonal-tree-identifiability/verifier.py
python -m unittest discover -s tests -p test_research_clonal_tree_identifiability.py -v
```

The default verifier enumerates through `n = 5`. Its `--max-n` option changes the
bound. Exhaustive target/candidate verification grows approximately as
`n^(2n-2)`; increasing that bound is optional and is not part of the under-60-second
test requirement. The dedicated test uses a fixed bounded workload.

Repository gates, without modifying existing files:

```bash
python scripts/validate.py
python -m unittest discover -s tests -v
node --test tests/*.test.*
```

The new test is discovered by the existing CI workflow. This directory does not
add dependencies, edit workflows, register a publication, publish clinical
claims, or alter the existing result/data catalog.

Observed local validation used Python 3.12.14: the nine dedicated tests passed
in about 2.3 seconds, all 187 repository Python tests passed, all 26 Node tests
passed, and `scripts/validate.py` passed. Timings depend on the machine; the
fixed test workload is comfortably below the requested 60-second budget.

## 9. Falsification log and hostile self-review

The review used separate mathematical derivations, independent tree generation,
an independently implemented matrix-factorization oracle, and exact boundary
examples. Agreement between AI agents is not external scientific review.

| Tempting claim or proof shortcut | Falsification attempt / outcome |
| --- | --- |
| Enough arbitrary samples must identify a tree. | False. Theorem 3 constructs the maximum ambiguity for every finite `m`, with strictly positive usages. |
| `m >= n` and full rank resolve the problem. | False. The same witnesses have rank `min(m,n)`. In particular `n = m = 3` still admits two trees. |
| Distinct positive values are sufficient genericity. | False. Superincreasing inputs and an open neighbourhood have `(n-1)!` explanations. |
| Two designed samples only work by making some clone absent. | Falsified as an objection: Theorem 1 makes every usage positive and proves an explicit positive radius in the full frequency space. |
| Realizing the correct dominance order alone fixes every parent. | False. For the chain, `(3/4,1/2,1/4)` has the exact total dominance order but also permits the fork `1 -> {2,3}`. Lemma 2 supplies the missing condition. |
| Lemma 2 can omit its dominance-graph hypothesis. | False. For edges `1 -> 2`, `1 -> 3`, `2 -> 4` and row `(9/10,1/2,1/5,3/10)`, the only nonchild ancestor pair is `(1,4)` and `f_4 > u_1`; nevertheless `1 -> 2 -> 4 -> 3` is another explanation. |
| Checking one reattachment proves uniqueness against simultaneous changes. | Inadequate without proof. The shallowest-receiver argument shows that a competing tree retains every original child at a selected receiver and therefore has an explicit sum violation. |
| Equality is harmless in the one-row criterion. | False. `(3/4,1/2,1/4)` lies on equality and admits both chain and fork; `(7/10,1/2,1/4)` admits only the chain. |
| Positive frequencies alone imply the distinct-column factorial maximum. | False. Three equal positive frequencies permit all six labelled chains, exceeding `(3-1)! = 2`. Some usages then vanish. |
| With positive tied columns, `n!` must be the maximum. | Also false: at `n = 9`, a founder of frequency `1/2` and eight other frequencies `1/32` admit **394,353** trees, exceeding `9! = 362,880`. The count is proved just below. This variant is not solved here. |
| The root convention does not affect thresholds. | False. Adding a fixed normal vertex changes `(n-1)!` to `n!` and the first ambiguous size from three to two malignant clones. |
| A novel upper bound contradicts the two-sample NP-hardness result. | False quantifier comparison. The construction takes the target tree as input; arbitrary inverse reconstruction remains a different problem. |

For the nine-node tied example, every nonfounder has at most one child because
two children of frequency `1/32` would exceed its own frequency. Thus the tree
is an unordered collection of ordered paths attached to the founder. All such
collections fit because the founder can support eight children. For exactly
`k` paths, order the eight labels, choose `k-1` of seven cut positions, and divide
by `k!` for the ordering of the paths. The exact total is

$$\sum_{k=1}^{8}\frac{8!}{k!}\binom{7}{k-1}=394353.$$

This is a separate combinatorial construction, not an assertion that all
nine-node trees were exhaustively enumerated. The verifier checks the arithmetic
and checks the same path-partition count against exhaustive enumeration at
smaller sizes.

**Boundary review.** The proofs treat `n = 1`, `n = 2`, root labels other than
the first label, unary vertices, zero frequencies, repeated columns, equality
in sum conditions, normal contamination, and simultaneous parent moves. The
construction does not require binary branching. Rows can be perturbed
independently; identical founder frequencies in the displayed design are not
being passed off as an almost-everywhere assumption. No failed conjecture above
is retained as a theorem.

## 10. Limitations and questions left open

The theorem describes mathematical existence when mixtures can be chosen using
the known target tree. It supplies neither a procedure for selecting physical
biopsies without knowing that tree nor a probabilistic sample complexity under
a specified biological sampling law. Spatial availability, time constraints,
finite sequencing depth, copy-number changes, hidden clones, and uncertain
mutation clustering can alter the problem substantially. No claim about tumour
biology or patient sampling follows from the construction alone.

The perturbation radius certifies feasibility of the same exact mathematical
tree throughout a small frequency ball. It is not a sequencing-error confidence
interval. Its exponential dependence on leaf count can make it impractically
small, despite polynomial bit length. We have not optimized this radius, imposed
a prescribed positive lower bound on every usage, or determined the exact
maximum when frequencies must be positive but repeated columns are allowed.
Those are unproved extensions, not results of this note.

The finite checks reach five malignant vertices for the main exhaustive table.
All larger-size statements rest on the written proofs. The independent test
checks mathematical implementation through a different algorithm; CI success
establishes neither bibliographic novelty nor biological or clinical validity.
The novelty assessment remains provisional and should be revised if an earlier
equivalent construction or theorem is identified.

## 11. Independent referee derivation and dependency/novelty ledger (17 November 2026)

**Disposition: Theorem 1 is proved under assumptions 1–6.** It is a
*designed-data existence* result, not a generic or observational sample-size
guarantee. The independent argument below starts from the frequency constraints
rather than assuming the result of the verifier's search.

**Alternative derivation.** For any candidate rooted tree and row of
frequencies, invert the descendant-sum transform by subtracting the immediately
child-descendant sums. Thus feasibility is equivalent to nonnegative
\(f_v-\sum_{c\text{ child of }v}f_c\), and the root's frequency is the sum
of all inferred malignant usages. This is an already-known sum condition
(El-Kebir et al. 2015).

For the proposed design, multiply each frequency by \(D\). The leaf weights
are \(W2^i\) and \(W2^{L-1-i}\), with unit weights at internal vertices.
If disjoint subtrees \(A,B\) have ordered depth-first leaf blocks, then
the first block's total numerator is at most
\(W(2^j-1)+(n-L)\), while the later block contains at least \(W2^j\).
Thus the later block is larger by at least \(W-(n-L)=L+1\);
reverse weights force the opposite inequality in the other row. Conversely,
if \(a\) is an ancestor of \(v\), the subtree at \(a\) includes at least
one positive usage outside the subtree at \(v\). Therefore the *only*
coordinatewise permitted edges are the actual ancestor pairs. The margins
are at least \(1/D\) for ancestor comparisons and \((L+1)/D\) for
incomparables.

Now suppose another feasible parent map used only those permitted edges.
Choose a vertex \(a\) of smallest original depth that acquires an original
nondirect descendant as a new child. Its old direct children cannot leave,
because their alternative parents would be **shallower receivers**.
Consequently its candidate residual is at most
\(u_a-f_v\le(1-W)/D=-n/D<0\), a contradiction. A path works with one
row because all pairs are comparable. Every single row also admits a sorted
chain with nonnegative consecutive differences, so a branching target
cannot be identified by one row. This independently establishes (3).

**Boundary and perturbation audit.** Changing one entry of \(F\) by at most
\(\rho=1/(4nD)\) changes any residual involving at most \(n\) frequency
entries by no more than \(1/(4D)\), and any pairwise difference by
no more than \(2\rho\). Thus the true positive usages and their row sum
margin persist, all permitted-edge comparisons persist, and the negative
residual at every shallowest receiver remains negative. This includes
the closed-radius endpoints, nonfounder root labels, unary vertices, and
the one-vertex case. The new rational-corner regression checks this full
closed perturbation claim for all nine three-vertex targets independently
of the verifier's sum-condition checker. The all-size assertion rests on
these inequalities, **not** on finite enumeration.

**Dependency and priority ledger.** The sum condition, one-row chain
obstruction, two-order representation of rooted-tree ancestry, and
indegree-product counting bound are established earlier results, not new
theorems of this programme. In particular, Pradhan–El-Kebir (2018) and
Qi et al. (2019) precede the maximum-count discussion, and
Marangola–Sheikholeslami–Bento had a **published** IEEE TCBB article on
PPM degeneracy dated 21 August 2026; its October arXiv revision is not the
sole primary record. After rechecking these sources, there is no support
for claims of the *first* identifiability or ambiguity analysis.
The positive two-row construction and its radius are a **narrow
candidate synthesis with unestablished priority**, not a demonstrated novel
biological finding. The exact one-row criterion is an elementary special
case of sum conditions and should not be advertised as a new inference
method. No counterexample to the stated Theorem 1 was found.

**What the software actually establishes.** The distinct parent-map
enumerator and Gauss–Jordan inversion check designed examples through
five vertices, 4,096 small-grid matrices for the conditional criterion,
and independent positive-margin checks. They do not enumerate arbitrary
real matrices, certify a biological sampling distribution, or establish
bibliographic priority. Real-world bounds under finite reads, copy-number
changes, uncertain purity, and inaccessible physical mixtures remain open.

## References

- Jiao W, Vembu S, Deshwar AG, Stein L, Morris Q (2014). *Inferring clonal
  evolution of tumors from single nucleotide somatic mutations*. BMC
  Bioinformatics 15:35. DOI: [10.1186/1471-2105-15-35][jiao].
- El-Kebir M, Oesper L, Acheson-Field H, Raphael BJ (2015). *Reconstruction of
  clonal trees and tumor composition from multi-sample sequencing data*.
  Bioinformatics 31:i62–i70. DOI: [10.1093/bioinformatics/btv261][ances].
- Malikic S, McPherson AW, Donmez N, Sahinalp SC (2015). *Clonality inference in
  multiple tumor samples using phylogeny*. Bioinformatics 31:1349–1356.
  DOI: [10.1093/bioinformatics/btv003][citup].
- El-Kebir M, Satas G, Oesper L, Raphael BJ (2016). *Inferring the Mutational
  History of a Tumor Using Multi-state Perfect Phylogeny Mixtures*. Cell Systems
  3:43–53. DOI: [10.1016/j.cels.2016.07.004][spruce]. The [2016 preprint][spruce-preprint]
  uses the title *Multi-State Perfect Phylogeny Mixture Deconvolution and
  Applications to Cancer Sequencing*.
- Satas G, Raphael BJ (2017). *Tumor phylogeny inference using tree-constrained
  importance sampling*. Bioinformatics 33:i152–i160.
  DOI: [10.1093/bioinformatics/btx270][pastri].
- Pradhan D, El-Kebir M (2018). *On the Non-uniqueness of Solutions to the Perfect
  Phylogeny Mixture Problem*. [Author-hosted manuscript][pradhan].
- Qi Y, Pradhan D, El-Kebir M (2019). *Implications of non-uniqueness in
  phylogenetic deconvolution of bulk DNA samples of tumors*. Algorithms for
  Molecular Biology 14:19. DOI: [10.1186/s13015-019-0155-6][qi].
- Sundermann LK et al. (2021). *Reconstructing tumor evolutionary histories and
  clone trees in polynomial-time with SubMARine*. PLOS Computational Biology
  17:e1008400. DOI: [10.1371/journal.pcbi.1008400][submarine].
- Marangola J, Sheikholeslami A, Bento J (2026). *Constraints on the Perfect
  Phylogeny Mixture Model and Their Effect on Reducing Degeneracy*.
  DOI: [10.1109/TCBBIO.2026.3726506](https://doi.org/10.1109/TCBBIO.2026.3726506).
  [arXiv:2512.24930v2][degeneracy], revised 2026-10-02.
- Brualdi RA, Jung HC, Trotter WT Jr (1994). *On the poset of all posets on n
  elements*. Discrete Applied Mathematics 50:111–123. [Author-hosted paper][posets].

[jiao]: https://doi.org/10.1186/1471-2105-15-35
[ances]: https://pmc.ncbi.nlm.nih.gov/articles/PMC4542783/
[citup]: https://doi.org/10.1093/bioinformatics/btv003
[spruce]: https://doi.org/10.1016/j.cels.2016.07.004
[spruce-preprint]: https://arxiv.org/abs/1604.02605
[pastri]: https://pmc.ncbi.nlm.nih.gov/articles/PMC5870673/
[pradhan]: https://www.el-kebir.net/public/preprints/RECOMBCG18_Pradhan_El-Kebir.pdf
[qi]: https://doi.org/10.1186/s13015-019-0155-6
[submarine]: https://doi.org/10.1371/journal.pcbi.1008400
[degeneracy]: https://arxiv.org/html/2512.24930v2
[posets]: https://trotter.math.gatech.edu/papers/90.pdf
