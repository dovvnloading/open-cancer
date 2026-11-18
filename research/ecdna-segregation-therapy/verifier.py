#!/usr/bin/env python3
"""Exact finite checks for the ecDNA segregation/therapy mathematics.

Nonclinical computational mathematics; no treatment advice.  Every arithmetic
operation in the verifier is rational, except the explicitly labelled Decimal
illustration.  Finite-generation calculations retain every reachable copy
number.  They never impose a biological copy-number cutoff.

The analytic recursion and the forward expected-count dynamic program below
are separate implementations.  Neither imports the stochastic simulator.
Numerical checks complement, and do not replace, the proofs in README.md.
"""

from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
import json
from math import comb
from typing import Callable, Iterable


Rational = int | str | Fraction


def _fraction(value: Rational, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, str, Fraction)):
        raise TypeError(f"{name} must be an int, Fraction, or rational string")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be a finite rational number") from exc


def _integer(value: int, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return value


def _parameters(u: Rational, s: Rational) -> tuple[Fraction, Fraction]:
    u, s = _fraction(u, "u"), _fraction(s, "s")
    if not (0 <= u <= 1 and 0 <= s <= 1):
        raise ValueError("require 0 <= u <= 1 and 0 <= s <= 1")
    return u, s


def _initial(initial: Iterable[int]) -> tuple[int, ...]:
    values = tuple(initial)
    for k in values:
        _integer(k, "initial copy number")
    return values


@lru_cache(maxsize=None)
def _kernel(k: int) -> tuple[Fraction, ...]:
    denominator = 1 << (2 * k)
    return tuple(Fraction(comb(2 * k, j), denominator) for j in range(2 * k + 1))


def binomial_kernel(k: int) -> tuple[Fraction, ...]:
    """Law of one daughter; the other daughter is exactly 2*k-j."""
    return _kernel(_integer(k, "k"))


def iterate_jet(s: Rational, z: Rational, generations: int) -> tuple[Fraction, Fraction, Fraction]:
    """Return T_s^n(z), its first derivative, and its second derivative."""
    _, s = _parameters(1, s)
    z = _fraction(z, "z")
    _integer(generations, "generations")
    if not 0 <= z <= 1:
        raise ValueError("require 0 <= z <= 1")
    value, first, second = z, Fraction(1), Fraction(0)
    for _ in range(generations):
        factor = s * (1 + value) / 2
        value, first, second = (
            s * (1 + value) ** 2 / 4,
            factor * first,
            factor * second + s * first**2 / 2,
        )
    return value, first, second


def _polynomial_jet(initial: tuple[int, ...], z: Fraction) -> tuple[Fraction, Fraction, Fraction]:
    value = sum((z**k for k in initial), Fraction(0))
    first = sum((k * z ** (k - 1) for k in initial if k), Fraction(0))
    second = sum((k * (k - 1) * z ** (k - 2) for k in initial if k >= 2), Fraction(0))
    return value, first, second


def _summary(cells: Fraction, zero: Fraction, copies: Fraction, squares: Fraction) -> dict:
    positive = cells - zero
    mean = copies / cells if cells else None
    variance = squares / cells - mean**2 if cells else None
    positive_mean = copies / positive if positive else None
    positive_variance = squares / positive - positive_mean**2 if positive else None
    return {
        "cells": cells,
        "zero": zero,
        "positive": positive,
        "copies": copies,
        "copy_squares": squares,
        "pooled_mean": mean,
        "pooled_variance": variance,
        "positive_mean": positive_mean,
        "positive_variance": positive_variance,
        "pooled_fano": variance / mean if mean else None,
    }


def pgf_summary(initial: Iterable[int], generations: int, u: Rational, s: Rational) -> dict:
    """Exact *pooled expected-count* summaries, not E[within-clone ratios]."""
    initial = _initial(initial)
    u, s = _parameters(u, s)
    _integer(generations, "generations")
    x, first, second = iterate_jet(s, 1, generations)
    y = iterate_jet(s, 0, generations)[0]
    f, fp, fpp = _polynomial_jet(initial, x)
    scale = (2 * u) ** generations
    cells = scale * f
    zero = scale * _polynomial_jet(initial, y)[0]
    copies = scale * fp * first
    squares = scale * (fpp * first**2 + fp * (second + first))
    return _summary(cells, zero, copies, squares)


def forward_expectations(initial: Iterable[int], generations: int, u: Rational, s: Rational) -> dict[int, Fraction]:
    """Independent forward DP over all reachable types at the stated horizon."""
    initial = _initial(initial)
    u, s = _parameters(u, s)
    _integer(generations, "generations")
    current = {k: Fraction(number) for k, number in Counter(initial).items()}
    for _ in range(generations):
        following: dict[int, Fraction] = {}
        for k, expected in current.items():
            factor = 2 * u * s**k * expected
            for j, probability in enumerate(binomial_kernel(k)):
                following[j] = following.get(j, Fraction(0)) + factor * probability
        current = following
    return current


def _generic_probability(division_probability: Callable[[int], Rational], k: int) -> Fraction:
    probability = _fraction(division_probability(k), f"division probability at k={k}")
    if not 0 <= probability <= 1:
        raise ValueError(f"division probability at k={k} must lie in [0, 1]")
    return probability


def generic_forward_expectations(initial: Iterable[int], generations: int, division_probability: Callable[[int], Rational]) -> dict[int, Fraction]:
    """Forward mean counts for any supplied exact per-type division law.

    The callable is evaluated once for each reachable type whose reproduction
    is needed.  Every reachable daughter type is retained; this is a horizon
    bound, not a biological cutoff.
    """
    initial = _initial(initial)
    _integer(generations, "generations")
    if not callable(division_probability):
        raise TypeError("division_probability must be callable")

    @lru_cache(maxsize=None)
    def probability(k: int) -> Fraction:
        return _generic_probability(division_probability, k)

    current = {k: Fraction(number) for k, number in Counter(initial).items()}
    for _ in range(generations):
        following: dict[int, Fraction] = {}
        for k, expected in current.items():
            factor = 2 * probability(k) * expected
            if not factor:
                continue
            for j, daughter_probability in enumerate(binomial_kernel(k)):
                following[j] = following.get(j, Fraction(0)) + factor * daughter_probability
        current = following
    return current


def summary_from_counts(counts: dict[int, Fraction]) -> dict:
    """Summarize an expected-count map from forward_expectations."""
    cells = sum(counts.values(), Fraction(0))
    copies = sum((k * value for k, value in counts.items()), Fraction(0))
    squares = sum((k * k * value for k, value in counts.items()), Fraction(0))
    return _summary(cells, counts.get(0, Fraction(0)), copies, squares)


@lru_cache(maxsize=None)
def _extinction(k: int, generations: int, u: Fraction, s: Fraction) -> Fraction:
    if k == 0:
        return Fraction(1)
    if generations == 0:
        return Fraction(0)
    division = u * s**k
    child_product = sum(
        (probability * _extinction(j, generations - 1, u, s)
         * _extinction(2 * k - j, generations - 1, u, s)
         for j, probability in enumerate(binomial_kernel(k))),
        Fraction(0),
    )
    return 1 - division + division * child_product


def extinction_by_generation(k: int, generations: int, u: Rational, s: Rational) -> Fraction:
    """Probability that one lineage has no positive cells at generation n."""
    _integer(k, "k")
    _integer(generations, "generations")
    u, s = _parameters(u, s)
    return _extinction(k, generations, u, s)


def population_extinction(initial: Iterable[int], generations: int, u: Rational, s: Rational) -> Fraction:
    """Exact finite-generation eradication probability for arbitrary founders."""
    initial = _initial(initial)
    u, s = _parameters(u, s)
    _integer(generations, "generations")
    probability = Fraction(1)
    for k in initial:
        probability *= _extinction(k, generations, u, s)
    return probability


def generic_extinction(initial: Iterable[int], generations: int, division_probability: Callable[[int], Rational]) -> Fraction:
    """Exact finite-generation positive eradication for an arbitrary p_k law.

    q_0 is one because a zero-copy founder cannot produce positive descendants,
    regardless of its own division probability.  The daughter product uses
    complementary counts j and 2*k-j, never independent sister draws.
    """
    initial = _initial(initial)
    _integer(generations, "generations")
    if not callable(division_probability):
        raise TypeError("division_probability must be callable")

    @lru_cache(maxsize=None)
    def probability(k: int) -> Fraction:
        return _generic_probability(division_probability, k)

    @lru_cache(maxsize=None)
    def extinct(k: int, n: int) -> Fraction:
        if k == 0:
            return Fraction(1)
        if n == 0:
            return Fraction(0)
        division = probability(k)
        if not division:
            return Fraction(1)
        offspring = sum(
            (daughter_probability * extinct(j, n - 1) * extinct(2 * k - j, n - 1)
             for j, daughter_probability in enumerate(binomial_kernel(k))),
            Fraction(0),
        )
        return 1 - division + division * offspring

    result = Fraction(1)
    for k in initial:
        result *= extinct(k, generations)
    return result


def _statistic(k: int, statistic: str) -> int:
    if statistic == "cells":
        return 1
    if statistic == "positive":
        return int(k > 0)
    if statistic == "zero":
        return int(k == 0)
    if statistic == "copies":
        return k
    if statistic == "copy_squares":
        return k * k
    raise ValueError("unknown statistic; use cells, positive, zero, copies, or copy_squares")


@lru_cache(maxsize=None)
def _raw_moments(k: int, generations: int, u: Fraction, s: Fraction, statistic: str) -> tuple[Fraction, Fraction]:
    if generations == 0:
        value = Fraction(_statistic(k, statistic))
        return value, value**2
    mean, second = Fraction(0), Fraction(0)
    for j, probability in enumerate(binomial_kernel(k)):
        left_mean, left_second = _raw_moments(j, generations - 1, u, s, statistic)
        right_mean, right_second = _raw_moments(2 * k - j, generations - 1, u, s, statistic)
        mean += probability * (left_mean + right_mean)
        second += probability * (left_second + right_second + 2 * left_mean * right_mean)
    division = u * s**k
    return division * mean, division * second


def lineage_raw_moments(k: int, generations: int, u: Rational, s: Rational, statistic: str = "positive") -> tuple[Fraction, Fraction]:
    """E[X], E[X^2] using correlated sisters and independent future lineages."""
    _integer(k, "k")
    _integer(generations, "generations")
    u, s = _parameters(u, s)
    _statistic(0, statistic)
    return _raw_moments(k, generations, u, s, statistic)


def population_raw_moments(initial: Iterable[int], generations: int, u: Rational, s: Rational, statistic: str = "positive") -> dict:
    """Raw moments of a total population statistic; variance is between clones."""
    initial = _initial(initial)
    u, s = _parameters(u, s)
    _integer(generations, "generations")
    _statistic(0, statistic)
    mean, variance = Fraction(0), Fraction(0)
    for k in initial:
        single_mean, single_second = _raw_moments(k, generations, u, s, statistic)
        mean += single_mean
        variance += single_second - single_mean**2
    return {"mean": mean, "second": variance + mean**2, "variance": variance}


def threshold_certificate(u: Rational, a: Rational) -> dict:
    """A rational parameterization avoids a rounded square root at criticality."""
    u = _fraction(u, "u")
    a = _fraction(a, "a")
    if not 0 < u <= 1 or not 0 < a < 1:
        raise ValueError("require 0 < u <= 1 and 0 < a < 1")
    s = 4 * a / (1 + a) ** 2
    eigenvalue = 4 * u * a / (1 + a)
    h_bound = 1 / (1 - a) ** 2
    quadratic_bound = u * s / (1 - s * a) ** 2
    epsilon = min(1 / h_bound, (eigenvalue - 1) / quadratic_bound) if eigenvalue > 1 else Fraction(0)
    return {"u": u, "a": a, "s": s, "lambda": eigenvalue, "H": h_bound, "C": quadratic_bound, "epsilon": epsilon}


def h_weight(k: int, a: Rational) -> Fraction:
    _integer(k, "k")
    a = _fraction(a, "a")
    if not 0 < a < 1:
        raise ValueError("require 0 < a < 1")
    return k * a ** (k - 1) if k else Fraction(0)


def check_eigen_identity(k: int, u: Rational, a: Rational) -> bool:
    """Check the exact eigenfunction and quadratic bound at any requested k."""
    _integer(k, "k", 1)
    certificate = threshold_certificate(u, a)
    u, s, a = certificate["u"], certificate["s"], certificate["a"]
    expected_weight, quadratic = Fraction(0), Fraction(0)
    for j, probability in enumerate(binomial_kernel(k)):
        left, right = h_weight(j, a), h_weight(2 * k - j, a)
        expected_weight += probability * (left + right)
        quadratic += probability * left * right
    expected_weight *= u * s**k
    quadratic *= u * s**k
    h = h_weight(k, a)
    exact_quadratic = u * s**k * k * (2 * k - 1) * a ** (2 * k - 2) / 2
    return (
        expected_weight == certificate["lambda"] * h
        and quadratic == exact_quadratic
        and 0 < h <= certificate["H"]
        and quadratic <= certificate["C"] * h
        and (not certificate["epsilon"] or certificate["epsilon"] * quadratic <= (certificate["lambda"] - 1) * h)
    )


def survival_lower_bound(initial: Iterable[int], u: Rational, a: Rational) -> Fraction:
    """Certified positive lower bound in the supercritical rational family."""
    initial = _initial(initial)
    certificate = threshold_certificate(u, a)
    upper_extinction = Fraction(1)
    for k in initial:
        upper_extinction *= 1 - certificate["epsilon"] * h_weight(k, certificate["a"])
    return 1 - upper_extinction


def threshold_regime(u: Rational, s: Rational) -> str:
    """Exact threshold classification for the theorem's open domain 0<s<1."""
    u, s = _parameters(u, s)
    if not 0 < s < 1:
        raise ValueError("threshold_regime uses the theorem's domain 0 < s < 1")
    if u <= Fraction(1, 2):
        return "eradication"
    comparison = 4 * u * u * s - (4 * u - 1)
    return "persistence" if comparison > 0 else "critical" if comparison == 0 else "eradication"


def _down(value: Fraction, grid: int) -> Fraction:
    return Fraction(value.numerator * grid // value.denominator, grid)


def _up(value: Fraction, grid: int) -> Fraction:
    return Fraction(-(-value.numerator * grid // value.denominator), grid)


def fano_limit_interval(s: Rational, steps: int = 80, digits: int = 30) -> tuple[Fraction, Fraction]:
    """Certified bounds for 1+sum D_j/(1+A_j), using outward rational rounding.

    A_0=D_0=1, A_(j+1)=s(1+A_j)^2/4 and D_(j+1)=s(1+A_j)D_j/2.
    All finite terms are enclosed.  The omitted series is bounded by
    s**steps/(1-s), because 0<A_j<=1 and 0<D_j<=s**j.
    """
    _, s = _parameters(1, s)
    if not 0 < s < 1:
        raise ValueError("the finite Fano limit requires 0 < s < 1")
    _integer(steps, "steps")
    _integer(digits, "digits", 1)
    grid = 10**digits
    alo = ahi = dlo = dhi = Fraction(1)
    lower = upper = Fraction(0)
    for _ in range(steps):
        lower = _down(lower + dlo / (1 + ahi), grid)
        upper = _up(upper + dhi / (1 + alo), grid)
        alo, ahi, dlo, dhi = (
            _down(s * (1 + alo) ** 2 / 4, grid),
            _up(s * (1 + ahi) ** 2 / 4, grid),
            _down(s * (1 + alo) * dlo / 2, grid),
            _up(s * (1 + ahi) * dhi / 2, grid),
        )
    return 1 + lower, 1 + upper + s**steps / (1 - s)


def decimal_pooled_moments(initial: Iterable[int], generations: int, u: Rational, s: Rational, precision: int = 60) -> dict:
    """NON-CERTIFIED Decimal illustration; use rational intervals for a proof."""
    initial = _initial(initial)
    u, s = _parameters(u, s)
    if not (u > 0 and s > 0):
        raise ValueError("the Decimal illustration requires u > 0 and s > 0; use exact finite-generation APIs for zero endpoints")
    _integer(generations, "generations")
    _integer(precision, "precision", 16)
    with localcontext() as context:
        context.prec = precision
        sd = Decimal(s.numerator) / Decimal(s.denominator)
        ud = Decimal(u.numerator) / Decimal(u.denominator)
        value, first, ratio = Decimal(1), Decimal(1), Decimal(0)
        for _ in range(generations):
            ratio += first / (1 + value)
            value, first = sd * (1 + value) ** 2 / 4, sd * (1 + value) * first / 2
        f = sum((value**k for k in initial), Decimal(0))
        fp = sum((k * value ** (k - 1) for k in initial if k), Decimal(0))
        fpp = sum((k * (k - 1) * value ** (k - 2) for k in initial if k >= 2), Decimal(0))
        if not f:
            return {"pooled_mean": None, "pooled_variance": None, "pooled_fano": None, "certified": False}
        mean = first * fp / f
        fano = 1 + ratio + first * (fpp / fp - fp / f) if fp else None
        return {"cells": (2 * ud) ** generations * f, "pooled_mean": mean, "pooled_variance": mean * fano if fano is not None else Decimal(0), "pooled_fano": fano, "certified": False}


def verify() -> dict:
    """A finite exact verification grid; the README supplies general proofs."""
    checks = Counter()
    for u, s in ((Fraction(1), Fraction(3, 4)), (Fraction(4, 5), Fraction(8, 9)), (Fraction(2, 3), Fraction(1, 2)), (Fraction(0), Fraction(0)), (Fraction(0), Fraction(1)), (Fraction(1), Fraction(0)), (Fraction(1), Fraction(1))):
        for initial in ((), (0,), (1,), (2, 2), (0, 1, 3)):
            for n in range(4):
                analytic = pgf_summary(initial, n, u, s)
                forward = summary_from_counts(forward_expectations(initial, n, u, s))
                assert analytic == forward
                checks["pgf_vs_forward"] += 1
                for statistic in ("cells", "positive", "copies"):
                    raw = population_raw_moments(initial, n, u, s, statistic)
                    assert raw["mean"] == analytic[statistic]
                    assert raw["variance"] >= 0
                    checks["raw_moment_mean"] += 1
                extinction = population_extinction(initial, n, u, s)
                assert 0 <= extinction <= 1
                assert 1 - extinction <= analytic["positive"]
                if n:
                    assert extinction >= population_extinction(initial, n - 1, u, s)
                checks["finite_extinction"] += 1
    for u, a in ((1, Fraction(1, 4)), (1, Fraction(1, 3)), (1, Fraction(1, 2)), (Fraction(3, 4), Fraction(3, 4))):
        for k in (1, 2, 3, 7, 19, 64):
            assert check_eigen_identity(k, u, a)
            checks["eigen_and_quadratic"] += 1
    for u, s, expected in ((1, Fraction(2, 3), "eradication"), (1, Fraction(3, 4), "critical"), (1, Fraction(4, 5), "persistence"), (Fraction(1, 2), Fraction(99, 100), "eradication")):
        assert threshold_regime(u, s) == expected
        checks["threshold_sign"] += 1
    for n in range(1, 5):
        fano = pgf_summary((1,), n, 1, Fraction(3, 4))["pooled_fano"]
        for k in (2, 5, 17):
            assert pgf_summary((k, k), n, Fraction(3, 5), Fraction(3, 4))["pooled_fano"] == fano
            checks["monomorphic_fano"] += 1
    for s in (Fraction(1, 2), Fraction(3, 4), Fraction(8, 9)):
        lower, upper = fano_limit_interval(s)
        # Independent exact short partial sum must lie below the full upper bound.
        exact_partial = Fraction(1)
        for j in range(7):
            value, first, _ = iterate_jet(s, 1, j)
            exact_partial += first / (1 + value)
        assert exact_partial <= upper
        assert lower <= exact_partial + s**7 / (1 - s)
        checks["fano_interval"] += 1
    u, s = Fraction(4, 5), Fraction(8, 9)
    middle = lambda k: u * s**k
    lower_law = lambda k: middle(k) * Fraction(k + 2, k + 3)
    upper_law = lambda k: middle(k) + (1 - middle(k)) / (k + 3)
    for initial in ((1,), (2, 3), (0, 1, 4)):
        for n in range(4):
            exact = population_extinction(initial, n, u, s)
            assert generic_extinction(initial, n, middle) == exact
            assert generic_extinction(initial, n, lower_law) >= exact >= generic_extinction(initial, n, upper_law)
            counts = generic_forward_expectations(initial, n, middle)
            assert summary_from_counts(counts) == pgf_summary(initial, n, u, s)
            checks["generic_comparison"] += 1
    return {"label": "nonclinical computational mathematics; finite checks, not proof", "arithmetic": "exact Fraction; certified outward rational intervals", "checks": dict(checks), "total_checks": sum(checks.values())}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--illustrate", action="store_true", help="also print an explicitly non-certified Decimal long-horizon illustration")
    args = parser.parse_args()
    result = verify()
    if args.illustrate:
        result["decimal_illustration"] = decimal_pooled_moments((1, 3), 100, 1, Fraction(3, 4))
        result["certified_fano_interval_s_3_4"] = fano_limit_interval(Fraction(3, 4))
    print(json.dumps(result, indent=2, default=str, sort_keys=True))


if __name__ == "__main__":
    main()
