# Target-only robustness does not establish predictor integrity

**Original synthetic probability mathematics, nonclinical and unreviewed.**
This note is reconstructed and rerun as described in
[RECONSTRUCTION.md](RECONSTRUCTION.md). It explains a logical distinction in
the barcode analysis; it does not fit biological observations, estimate an
index-hopping rate, or establish historical mathematical priority.

The empirical analysis compares a forecast using a culture's own early
measurements with one that also uses early measurements from other cultures.
Two sensitivity operations act on different arguments of the score:

- Target deletion changes the empirical later distribution, with both
  predictions fixed.
- Predictor perturbation changes the information that produces a forecast.
  The empirical donor-only intervention also refits its parameters.

Passing the first operation cannot certify the second. The construction below
shows this even without refitting: fixed exponents suffice.

## 1. The exact target-deletion problem

Let p be an empirical probability vector, let q0 and q1 be positive fixed
forecasts, and define d_i = log(q1_i/q0_i). The original incremental score is

\[
G=\sum_i p_i d_i.
\]

If mass e_i is removed from each category, with 0 <= e_i <= p_i and
sum(e_i)=epsilon<1, the score on the retained normalized target is

\[
L(e)=\frac{\sum_i(p_i-e_i)d_i}{1-\epsilon}.
\]

For fixed epsilon, maximizing the removed weighted sum minimizes L. Suppose
d_i>d_j, category i is not fully removed, and category j has positive removed
mass. Transferring a small amount of removal from j to i decreases the
retained numerator by that amount times d_i-d_j. At an optimum, removal
therefore exhausts the highest-contrast categories first, with at most one
fractionally removed boundary category, aside from ties.

This exchange argument proves sharpness and constructs an attaining witness.
The optimal retained mean is nonincreasing as epsilon grows: the next removed
contrast is at least the mean of the currently retained distribution. Thus
the same solution is optimal if the constraint permits **at most** epsilon
removed mass. The constant-contrast case can attain the same value at every
budget. Fractional mass is allowed; this is a continuous sensitivity bound,
not a requirement that a physical read be fractionally deleted.

For G>0, a tipping budget solves

\[
\sum_i e_i d_i=G
\]

along that ordered removal path. If every occupied contrast is positive, no
budget below one attains a nonpositive mean. If the minimum occupied contrast
is zero, deleting all positive-contrast mass can attain zero while retaining
the zero-contrast observations. If G<=0, the tipping budget is already zero.

[robustness.py](robustness.py) implements these calculations and can return
the deletion/retention witness for numerical checking. Its default output
contains aggregates only. The exact optimization proof is separate from
ordinary floating-point evaluation: numerical values are not outward-rounded
interval certificates. Extremely degenerate tipping roots that round to one
are marked unresolved at floating-point precision rather than declared
attainable.

## 2. A three-category predictor-contamination construction

Choose parameters

\[
0<w<1,\qquad 0<b<1,\qquad 0<\epsilon<1,
\qquad 0<m<\min\{w,(1-w)/2\}.
\]

Define an own predictor x, a clean donor predictor z, and an observed target y:

\[
x=(w,m,1-w-m),\qquad z=(m,1-2m,m),\qquad y=(1,0,0).
\]

All entries of x and z are strictly positive and each vector sums to one.
Use the fixed geometric forecast with exponents a=1 and b:

\[
Q_i(x,z)=\frac{x_i z_i^b}{\sum_j x_j z_j^b},
\qquad Q_{0i}=x_i.
\]

The clean target probability has the closed form

\[
Q_1(x,z)
=\frac{w m^b}{(1-m)m^b+m(1-2m)^b}
=\frac{w}{1-m+m^{1-b}(1-2m)^b}.
\]

Because 1-b>0, letting m decrease to zero at fixed w and b gives
Q_1(x,z) -> w. The clean incremental target score therefore tends to
log(w/w)=0.

Now contaminate the donor with the outcome category:

\[
z_\epsilon=(1-\epsilon)z+\epsilon y
=((1-\epsilon)m+\epsilon,
  (1-\epsilon)(1-2m),(1-\epsilon)m).
\]

Its total-variation distance from z is

\[
\operatorname{TV}(z_\epsilon,z)
=\epsilon\operatorname{TV}(y,z)
=\epsilon(1-m)<\epsilon.
\]

The three unnormalized forecast weights are

\[
w((1-\epsilon)m+\epsilon)^b,\quad
m((1-\epsilon)(1-2m))^b,\quad
(1-w-m)((1-\epsilon)m)^b.
\]

At fixed positive w, epsilon and b, the first weight tends to w*epsilon^b>0,
and the other two tend to zero. Consequently,

\[
Q_1(x,z_\epsilon)\longrightarrow1,
\qquad
G_\epsilon=\log\frac{Q_1(x,z_\epsilon)}{w}
\longrightarrow-\log w>0.
\]

Since y is a point mass, deleting any fraction eta<1 from the target leaves
the same normalized distribution y. The target-only lower gain equals
G_epsilon for every such eta. Whenever G_epsilon>0, the target-only tipping
mass is unattainable below one.

This proves that an arbitrarily small positive predictor-contamination budget
can produce a finite apparent gain that survives every target-only deletion
budget below one, **provided the probability floor may become sufficiently
small**. It does not assert a discontinuity at a fixed positive floor.

The order of limits matters. First fix w, b, and an arbitrarily small positive
epsilon; then take m sufficiently small. To obtain arbitrarily large gains,
choose w smaller and choose m sufficiently small for that w and epsilon.
Holding m fixed while epsilon tends to zero is a different limit, covered by
the continuity bound below.

## 3. An explicit numerical example

With w=0.2, b=1/4, m=1e-8, and epsilon=1e-4 (0.01% donor mixture mass), the
fresh synthetic calculation gives approximately:

| Quantity | Value |
| --- | ---: |
| Baseline target probability | 0.200000000 |
| Clean augmented target probability | 0.199999802 |
| Contaminated augmented target probability | 0.714295665 |
| Clean incremental score | -0.000000990 nats |
| Contaminated incremental score | +1.272979607 nats |

These are deliberately constructed probability vectors, not estimates of
barcoding error, cancer resistance, or the contamination needed to explain
the empirical result. The point-mass target makes the logical separation
especially transparent; no claim that the experimental target is a point
mass is needed or made.

[predictor_counterexample.py](predictor_counterexample.py) produces a small
declared grid of examples and labels each one `SYNTHETIC`. Its total-variation
field uses the exact expression epsilon*(1-m), avoiding subtraction of two
nearly unit numbers at very small epsilon. Ordinary floating-point values can
round a strictly smaller exact distance to the budget itself.

## 4. The same issue when contamination precedes smoothing

The construction is not restricted to contamination injected after a
smoothing operation. Start with raw profiles

\[
x^{raw}=(w,0,1-w),\quad z^{raw}=(0,1,0),
\quad z^{raw}_\epsilon=(\epsilon,1-\epsilon,0).
\]

Apply the same uniform smoothing to each predictor:

\[
S_\delta(v)=(1-\delta)v+\delta(1/3,1/3,1/3),\qquad 0<\delta<1.
\]

With a=1 and 0<b<1, the clean first and third weights are asymptotic to
w*(delta/3)^b and (1-w)*(delta/3)^b. The second is of order delta, which is
smaller than delta^b. The clean target probability tends to w.

After raw contamination, the first weight instead tends to w*epsilon^b>0;
the other weights vanish. The contaminated target probability tends to one
as delta tends to zero, and the same positive limiting gain follows. The
smoothed donor distance is (1-delta)*epsilon. For every fixed positive delta,
the forecast remains continuous in epsilon. `before_smoothing_example` in the
synthetic script supplies finite numerical examples of this construction.

## 5. A fixed-floor continuity bound

Let x be strictly positive and fix both forecast exponents, let min(z_i)>=m>0, and let
z_epsilon=(1-epsilon)z+epsilon*v for any probability vector v, with
0<=epsilon<1 and b>=0. Then

\[
L=1-\epsilon\le\frac{z_{\epsilon i}}{z_i}
\le1-\epsilon+\frac{\epsilon}{m}=U.
\]

Writing r_i=z_epsilon_i/z_i, normalization of the geometric forecast gives

\[
\frac{Q_i(x,z_\epsilon)}{Q_i(x,z)}
=\frac{r_i^b}{\sum_j Q_j(x,z)r_j^b}.
\]

Both numerator and denominator lie between L^b and U^b. Therefore every
per-category log-ratio has absolute value bounded by b*log(U/L), and for any
target distribution p,

\[
\boxed{\left|\mathbb E_p\log\frac{Q(x,z_\epsilon)}{Q(x,z)}\right|
\le b\log\left(1+\frac{\epsilon}{(1-\epsilon)m}\right).}
\]

For fixed m this bound tends to zero with epsilon. Its dependence on m
explains why a very small uniform floor can give weak protection against
small predictor contamination. The bound concerns fixed x and exponents;
it does not bound refitting the empirical nested procedure, changing its
support, or perturbing both own and donor predictors simultaneously. The
bound is a valid uniform bound, not asserted to be the sharpest bound for
every constrained mixture class.

## 6. Verification and interpretation

The target-deletion tests independently enumerate rational vertices of the
box-plus-equality feasible set using `fractions.Fraction`; they do not repeat
the sorting algorithm. They also independently minimize the mass needed to
erase positive gain. Witness, tie, zero-support, endpoint and invalid-input
checks cover implementation boundaries.

The synthetic tests use 80-digit Decimal arithmetic and fourth roots to
reproduce the numerical forecasts, independently of log-space normalization.
They check point-mass target invariance, multiple small contamination budgets,
the ordered-limit examples, raw-before-smoothing examples, and the fixed-floor
bound at several mixture distributions and exponents. Finite checks support
the implementation; the all-parameter statements follow from the proofs.

Run from the repository root:

```sh
python -B research/barcode-lineage-memory/predictor_counterexample.py \
  --output .compute/barcode-lineage-memory/predictor_counterexamples.json
python -B -m unittest discover -s tests -p 'test_research_barcode_robustness.py' -v
python -B -m unittest discover -s tests -p 'test_research_predictor_contamination.py' -v
```

The empirical donor-only sensitivity is a different operation from adding
the outcome mixture in this construction. The theorem is an explanation of
a possible failure mode, not a fitted explanation of the experimental data.
It neither identifies which experimental reads are erroneous nor rejects the
source study's full mechanistic analysis. No historical priority claim for
general contamination instability, index hopping, or log-score sensitivity
is established by this note.
