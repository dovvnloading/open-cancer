#!/usr/bin/env python3
"""Exact algebra and certified rational enclosures for the README theorems.

Nonclinical computational research; no medical advice. Fraction operations and
the interval routines are exact. Decimal is used only to display square roots,
exponentials, and logarithms, never to decide a certified inequality. Numerical
checks corroborate the proofs; they are not a formal proof assistant.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import permutations
import json
from math import factorial, isqrt
import random


def decimal(x: F | int, precision: int = 60) -> Decimal:
    """A rounded Decimal presentation of a rational, with explicit precision."""
    x = F(x)
    with localcontext() as ctx:
        ctx.prec = precision
        return Decimal(x.numerator) / Decimal(x.denominator)


@dataclass(frozen=True)
class Branching:
    """Poisson initiation; binary birth, death, and forward mutation per stage.

    conversion[i] means a forward mutation removes its parent. At the terminal
    stage the convention has no effect on time to the first malignant event.
    Empty stage tuples describe k=1 (direct malignant immigration).
    """

    initiation: F
    birth: tuple[F, ...]
    death: tuple[F, ...]
    mutation: tuple[F, ...]
    conversion: tuple[bool, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "initiation", F(self.initiation))
        for name in ("birth", "death", "mutation"):
            object.__setattr__(self, name, tuple(F(x) for x in getattr(self, name)))
        m = len(self.mutation)
        if not (len(self.birth) == len(self.death) == len(self.conversion) == m):
            raise ValueError("Stage arrays must have the same length")
        if self.initiation <= 0 or any(x <= 0 for x in self.mutation):
            raise ValueError("Initiation and mutation rates must be positive")
        if any(x < 0 for x in self.birth + self.death):
            raise ValueError("Birth and death rates must be nonnegative")
        if any(c not in (False, True) for c in self.conversion):
            raise ValueError("Mutation convention must be a binary flag")

    @property
    def k(self) -> int:
        return len(self.mutation) + 1


def _product_coefficient(a: list[F], b: list[F], n: int) -> F:
    return sum((a[j] * b[n - j] for j in range(n + 1)), F(0))


def branching_jet(p: Branching, order: int) -> tuple[F, ...]:
    """Ordinary Taylor coefficients h(t)=sum c[n] t**n, exactly (Theorem 1)."""
    if order < 0:
        raise ValueError("order must be nonnegative")
    if p.k == 1:
        return (p.initiation,) + (F(0),) * order
    m = p.k - 1
    q = [[F(0) for _ in range(order + 1)] for _ in range(m)]
    terminal = [F(1)] + [F(0)] * order
    for n in range(order):
        for i in range(m):
            nxt = q[i + 1] if i + 1 < m else terminal
            mu, alpha, beta = p.mutation[i], p.birth[i], p.death[i]
            converting = int(p.conversion[i])
            rhs = mu * nxt[n] + (alpha - beta - converting * mu) * q[i][n]
            rhs -= alpha * _product_coefficient(q[i], q[i], n)
            rhs -= (1 - converting) * mu * _product_coefficient(q[i], nxt, n)
            q[i][n + 1] = rhs / (n + 1)
    return tuple(p.initiation * x for x in q[0])


def tsce_invariants(p: Branching) -> tuple[F, F, F]:
    """The complete (a,b,c) in h'=a+b*h-c*h*h; a>0 (Theorem 2)."""
    if p.k != 2:
        raise ValueError("TSCE requires k=2")
    return (p.initiation * p.mutation[0],
            p.birth[0] - p.death[0] - p.mutation[0],
            p.birth[0] / p.initiation)


def recover_tsce(jet: tuple[F, ...]) -> tuple[F, F, F]:
    if len(jet) < 4 or jet[0] != 0 or jet[1] <= 0:
        raise ValueError("Need a positive clean-birth TSCE jet through order 3")
    a, second, third = jet[1], 2 * jet[2], 6 * jet[3]
    b = second / a
    c = (b * b * a - third) / (2 * a * a)
    return a, b, c


def polynomial_product(a: tuple[F, ...], b: tuple[F, ...]) -> tuple[F, ...]:
    """Ascending powers, including leading zero coefficients."""
    result = [F(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] += ai * bj
    return tuple(result)


def spectral_signature(amplitude: F, growth: tuple[F, ...],
                       delays: tuple[F, ...] = ()) -> tuple[F, tuple[F, ...]]:
    """Reduced C/P(s), with P monic; delay entries are positive rates (Thm 4)."""
    if amplitude <= 0 or any(rate <= 0 for rate in delays):
        raise ValueError("Positive amplitude and delay rates required")
    denominator = (F(0), F(1))  # the reservoir's compulsory zero pole
    numerator = F(amplitude)
    for rate in growth:
        denominator = polynomial_product(denominator, (-rate, F(1)))
    for rate in delays:
        numerator *= rate
        denominator = polynomial_product(denominator, (rate, F(1)))
    return numerator, denominator


def linear_jet(amplitude: F, growth: tuple[F, ...], order: int) -> tuple[F, ...]:
    """Exact cascade Taylor series, independently of polynomial multiplication."""
    if amplitude <= 0 or order < 0:
        raise ValueError("Positive amplitude and nonnegative order required")
    # Unit forcing x0=1, x_i'=x_(i-1)+r_i*x_i, h=A*x_m.
    states = [[F(0)] * (order + 1) for _ in growth]
    forcing = [F(1)] + [F(0)] * order
    if not growth:
        return tuple(amplitude * v for v in forcing)
    for n in range(order):
        for i, rate in enumerate(growth):
            previous = forcing if i == 0 else states[i - 1]
            states[i][n + 1] = (previous[n] + rate * states[i][n]) / (n + 1)
    return tuple(amplitude * v for v in states[-1])


def exponential_delay_jet(jet: tuple[F, ...], rate: F) -> tuple[F, ...]:
    """Solve z'=rate*(h-z), z(0)=0 using exact coefficients (Theorem 5)."""
    if rate <= 0:
        raise ValueError("Positive delay rate required")
    z = [F(0)] * len(jet)
    for n in range(len(jet) - 1):
        z[n + 1] = rate * (jet[n] - z[n]) / (n + 1)
    return tuple(z)


def exp_interval(x: F, terms: int = 120) -> tuple[F, F]:
    """Certified enclosure of exp(x); no floating point or Decimal decisions.

    For x>=0 the omitted positive series is bounded by its first term divided
    by 1-x/(terms+2). Negative x is enclosed by reciprocal positive bounds.
    """
    x = F(x)
    if terms < 0 or abs(x) >= terms + 2:
        raise ValueError("Need terms+2 > abs(x) for geometric tail bound")
    if x < 0:
        lo, hi = exp_interval(-x, terms)
        return 1 / hi, 1 / lo
    term = F(1)
    total = term
    for n in range(1, terms + 1):
        term *= x / n
        total += term
    following = term * x / (terms + 1)
    tail = following / (1 - x / (terms + 2))
    return total, total + tail


def sqrt_interval(x: F, digits: int = 50) -> tuple[F, F]:
    """Certified square-root enclosure using integer arithmetic alone."""
    x = F(x)
    if x < 0 or digits < 0:
        raise ValueError("Nonnegative x and digits required")
    scale = 10 ** digits
    integer = isqrt((x.numerator * scale * scale) // x.denominator)
    lo = F(integer, scale)
    return (lo, lo) if lo * lo == x else (lo, F(integer + 1, scale))


def exp_integral_interval(rate: F, a: F, b: F,
                          terms: int = 160) -> tuple[F, F]:
    if not (0 <= a < b):
        raise ValueError("Require 0 <= a < b")
    if rate == 0:
        return b - a, b - a
    alo, ahi = exp_interval(rate * a, terms)
    blo, bhi = exp_interval(rate * b, terms)
    lo, hi = blo - ahi, bhi - alo
    return (lo / rate, hi / rate) if rate > 0 else (hi / rate, lo / rate)


def power_bin(k: int, a: F, b: F) -> F:
    """Average t**(k-1)/(k-1)! on [a,b); exact, including k=1."""
    a, b = F(a), F(b)
    if k < 1 or not (0 <= a < b):
        raise ValueError("Require k>=1 and 0<=a<b")
    return (b ** k - a ** k) / (factorial(k) * (b - a))


def one_growth_bin_interval(k: int, growth: F, a: F, b: F,
                            terms: int = 100) -> tuple[F, F]:
    """Certified one-nonzero-pole limiting shape average (Theorem 3).

    f=sum_(n>=0) g**n*t**(n+k-1)/(n+k-1)!. Bound the absolute tail
    by (b/(b-a))*b**(k-1)* (abs(g)*b)**(N+1)/(N+1+k)!
    times 1/(1-abs(g)*b/(N+k+2)). Works for negative g as well.
    """
    if k < 2 or not (0 <= a < b) or terms < 0:
        raise ValueError("Require k>=2, 0<=a<b, and terms>=0")
    growth, a, b = F(growth), F(a), F(b)
    if growth == 0:
        v = power_bin(k, a, b)
        return v, v
    ratio = abs(growth) * b / (terms + k + 2)
    if ratio >= 1:
        raise ValueError("Too few terms for tail bound")
    total = sum((growth ** n * (b ** (n + k) - a ** (n + k)) /
                 (factorial(n + k) * (b - a)) for n in range(terms + 1)), F(0))
    first_omitted_bound = (abs(growth) ** (terms + 1) * b ** (terms + 1 + k) /
                           (factorial(terms + 1 + k) * (b - a)))
    tail = first_omitted_bound / (1 - ratio)
    return total - tail, total + tail


def calibration_ratio(lower: tuple[F, ...], upper: tuple[F, ...],
                      shape: tuple[F, ...]) -> F:
    """Exact minimal rho**2 for interval data; rho>=1 (Theorem 7)."""
    if not lower or not (len(lower) == len(upper) == len(shape)):
        raise ValueError("Equal nonempty vectors required")
    if any(l <= 0 or u < l or f <= 0 for l, u, f in zip(lower, upper, shape)):
        raise ValueError("Positive data intervals and shapes required")
    L = max(l / f for l, f in zip(lower, shape))
    U = min(u / f for u, f in zip(upper, shape))
    return max(F(1), L / U)


def calibration_ratio_interval(lower: tuple[F, ...], upper: tuple[F, ...],
                               shape_intervals: tuple[tuple[F, F], ...]
                               ) -> tuple[F, F]:
    """Exact enclosure of rho**2 when each shape has a rational enclosure."""
    if not (len(lower) == len(upper) == len(shape_intervals)) or not lower:
        raise ValueError("Equal nonempty vectors required")
    if any(l <= 0 or u < l for l, u in zip(lower, upper)):
        raise ValueError("Positive data intervals required")
    if any(a <= 0 or b < a for a, b in shape_intervals):
        raise ValueError("Positive shape enclosures required")
    Llo = max(l / hi for l, (lo, hi) in zip(lower, shape_intervals))
    Lhi = max(l / lo for l, (lo, hi) in zip(lower, shape_intervals))
    Ulo = min(u / hi for u, (lo, hi) in zip(upper, shape_intervals))
    Uhi = min(u / lo for u, (lo, hi) in zip(upper, shape_intervals))
    return max(F(1), Llo / Uhi), max(F(1), Lhi / Ulo)


def flexible_shape_calibration(lower: tuple[F, ...], upper: tuple[F, ...],
                               shape_ranges: tuple[tuple[F, F], ...]) -> F:
    """Minimal rho**2 allowing each shape independently anywhere in its range.

    Unlike calibration_ratio_interval, this optimizes over an allowed nuisance
    (e.g. unknown within-bin age weights); it is not a numerical error bound.
    """
    calibration_ratio_interval(lower, upper, shape_ranges)  # validation
    L = max(l / hi for l, (lo, hi) in zip(lower, shape_ranges))
    U = min(u / lo for u, (lo, hi) in zip(upper, shape_ranges))
    return max(F(1), L / U)


def amplitude_interval(lower: tuple[F, ...], upper: tuple[F, ...],
                       shape: tuple[F, ...], rho: F) -> tuple[F, F]:
    """Feasible amplitude interval; empty exactly when returned left > right."""
    calibration_ratio(lower, upper, shape)  # validation
    if rho < 1:
        raise ValueError("rho must be >=1")
    return (max(l / f for l, f in zip(lower, shape)) / rho,
            min(u / f for u, f in zip(upper, shape)) * rho)


def poisson_bounds(means: tuple[F, ...], alternatives: tuple[F, ...]) -> dict[str, str]:
    """Decimal presentations of exact formulas in Theorem 6, not enclosures."""
    if not means or len(means) != len(alternatives):
        raise ValueError("Equal nonempty vectors required")
    if any(x < 0 for x in means + alternatives):
        raise ValueError("Poisson means cannot be negative")
    with localcontext() as ctx:
        ctx.prec = 60
        D = sum(((decimal(x).sqrt() - decimal(y).sqrt()) ** 2
                 for x, y in zip(means, alternatives)), Decimal(0))
        # e^-D may round to 1 at extremely small D; max protects the display.
        affinity_bound = (1 - max(Decimal(0), 1 - (-D).exp()).sqrt()) / 2
        distance = sum((abs(x - y) for x, y in zip(means, alternatives)), F(0))
        coupling_bound = (-decimal(distance)).exp() / 2
        return {"squared_sqrt_mean_distance": str(D),
                "affinity_error_lower_bound": str(affinity_bound),
                "coupling_error_lower_bound": str(coupling_bound)}


def poisson_bounds_exact(means: tuple[F, ...], alternatives: tuple[F, ...],
                         digits: int = 40) -> dict[str, F]:
    """Rationally certified enclosures of Theorem 6's two error bounds.

    These are enclosures of the lower-bound FORMULAS, not an enclosure of the
    actual optimal error (which can be larger). No irrational value is exact.
    """
    if not means or len(means) != len(alternatives):
        raise ValueError("Equal nonempty vectors required")
    if any(x < 0 for x in means + alternatives):
        raise ValueError("Poisson means cannot be negative")
    dlo = dhi = F(0)
    for x, y in zip(means, alternatives):
        lx, ux = sqrt_interval(x, digits)
        ly, uy = sqrt_interval(y, digits)
        lower_abs = max(F(0), lx - uy, ly - ux)
        upper_abs = max(abs(lx - uy), abs(ux - ly))
        dlo += lower_abs * lower_abs
        dhi += upper_abs * upper_abs
    e_lower, _ = exp_interval(-dhi, max(120, int(dhi) + 20))
    _, e_upper = exp_interval(-dlo, max(120, int(dlo) + 20))
    _, root_upper = sqrt_interval(1 - e_lower, digits)
    root_lower, _ = sqrt_interval(max(F(0), 1 - e_upper), digits)
    distance = sum((abs(x - y) for x, y in zip(means, alternatives)), F(0))
    clo, chi = exp_interval(-distance, max(120, int(distance) + 20))
    return {"distance_lower": dlo, "distance_upper": dhi,
            "affinity_bound_lower": max(F(0), (1 - root_upper) / 2),
            "affinity_bound_upper": (1 - root_lower) / 2,
            "coupling_bound_lower": clo / 2, "coupling_bound_upper": chi / 2}


def exact_checks() -> dict[str, int | str]:
    """Finite rational falsification attempts, plus certified analytic examples."""
    rng = random.Random(20261117)
    onset_cases = 0
    for k in range(1, 9):
        for trial in range(16):
            m = k - 1
            p = Branching(F(rng.randrange(1, 8), 7),
                          tuple(F(rng.randrange(0, 5), 11) for _ in range(m)),
                          tuple(F(rng.randrange(0, 5), 13) for _ in range(m)),
                          tuple(F(rng.randrange(1, 6), 17) for _ in range(m)),
                          tuple(bool(rng.randrange(2)) for _ in range(m)))
            jet = branching_jet(p, k + 3)
            leading = p.initiation
            for mu in p.mutation:
                leading *= mu
            assert all(c == 0 for c in jet[:k - 1])
            assert jet[k - 1] == leading / factorial(k - 1)
            onset_cases += 1
    # Whole-curve equivalence is proved in the README; check many exact terms.
    first = Branching(F(2), (F(3),), (F(1),), (F(1),), (False,))
    second = Branching(F(4), (F(6),), (F(9, 2),), (F(1, 2),), (True,))
    assert tsce_invariants(first) == tsce_invariants(second)
    assert branching_jet(first, 18) == branching_jet(second, 18)
    assert recover_tsce(branching_jet(first, 3)) == tsce_invariants(first)
    # No-division boundary and its feasible equivalence class.
    no_birth = Branching(F(2), (F(0),), (F(1),), (F(2),), (False,))
    no_birth_scaled = Branching(F(4), (F(0),), (F(2),), (F(1),), (False,))
    assert branching_jet(no_birth, 16) == branching_jet(no_birth_scaled, 16)
    assert recover_tsce(branching_jet(no_birth, 3)) == tsce_invariants(no_birth)
    growth = (F(0), F(1, 5), F(-1, 3), F(1, 5))
    signatures = {spectral_signature(F(7), tuple(p)) for p in permutations(growth)}
    jets = {linear_jet(F(7), tuple(p), 16) for p in permutations(growth)}
    assert len(signatures) == len(jets) == 1
    # k=3 plus diagnostic lag and k=4 with a subcritical intermediate stage.
    assert spectral_signature(F(2), (F(0), F(1, 5)), (F(3),)) == \
           spectral_signature(F(6), (F(0), F(1, 5), F(-3)))
    # Exact stochastic conversion delay matches its independent scalar filter.
    downstream = Branching(F(2), (F(1, 3), F(1, 2)),
                           (F(1, 4), F(1, 5)), (F(2, 5), F(1, 7)),
                           (False, False))
    extra = Branching(F(2), (F(0),) + downstream.birth,
                      (F(0),) + downstream.death,
                      (F(5),) + downstream.mutation,
                      (True,) + downstream.conversion)
    assert branching_jet(extra, 18) == \
           exponential_delay_jet(branching_jet(downstream, 18), F(5))
    # Mutation with parent retention does not implement an exponential delay.
    wrong = Branching(extra.initiation, extra.birth, extra.death,
                      extra.mutation, (False,) + downstream.conversion)
    assert branching_jet(wrong, 8) != branching_jet(extra, 8)
    for sign in (-1, 0, 1):
        lo, hi = exp_interval(F(sign, 2), 60)
        assert 0 < lo <= hi and hi - lo < F(1, 10**80)
    # Sharp, nontrivial calibration boundary with rational rho=2.
    y, f = (F(2), F(8)), (F(1), F(1))
    assert calibration_ratio(y, y, f) == 4
    assert amplitude_interval(y, y, f, F(2)) == (F(4), F(4))
    left, right = amplitude_interval(y, y, f, F(1999, 1000))
    assert left > right
    # Rounded-interval intersection can be feasible at rho=1.
    assert calibration_ratio((F(1), F(2)), (F(2), F(4)), (F(1), F(2))) == 1
    # Positive kernel and rigorous signed-growth series enclosures.
    for k in (2, 3, 6):
        for g in (F(-1, 20), F(0), F(1, 20)):
            lo, hi = one_growth_bin_interval(k, g, F(30), F(35), 70)
            assert 0 < lo <= hi and hi - lo < F(1, 10**40)
    return {"status": "all finite exact checks passed",
            "onset_parameter_cases": onset_cases,
            "spectral_permutations": 24,
            "seed": 20261117}


if __name__ == "__main__":
    print(json.dumps(exact_checks(), indent=2, sort_keys=True))
