#!/usr/bin/env python3
"""Stable mathematical helpers for a nonclinical aggregate research note.

Python standard library only. No network, patient records, clinical predictions,
or medical advice. Finite checks validate implementation; the general proofs are
in README.md. An exponential intensity per division clock is not a Bernoulli
probability per discrete biological division.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys


def _stage(k: int) -> None:
    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise ValueError("k must be a positive integer")


def _positive(value: float, label: str) -> None:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{label} must be finite and positive")


def _log_polynomial(k: int, x: float) -> float:
    """log(sum_{j=0}^{k-1} x**j/j!), evaluated without overflow."""
    if x == 0 or k == 1:
        return 0.0
    lx = math.log(x)
    logs = [j * lx - math.lgamma(j + 1) for j in range(k)]
    largest = max(logs)
    return largest + math.log(math.fsum(math.exp(t - largest) for t in logs))


def _small_log_tail(k: int, x: float) -> float:
    """log P(Poisson(x)>=k) for 0<x<k, with a positive tail series.

    After term n the successive ratios are <= x/(k+n+1)<1. The omitted
    relative tail is bounded by the next term divided by one minus that ratio.
    This is a truncation control, not an enclosure of floating-point roundoff.
    """
    term = 1.0
    terms = [term]
    for n in range(1, 100000):
        term *= x / (k + n)
        terms.append(term)
        total = math.fsum(terms)
        ratio = x / (k + n + 1)
        if term * ratio / (1.0 - ratio) <= 2e-17 * total:
            return -x + k * math.log(x) - math.lgamma(k + 1) + math.log(total)
    raise ArithmeticError("Poisson tail series did not converge")


def log_psi(k: int, x: float) -> float:
    """Natural log of the exact Erlang lineage cumulative hazard.

    Keeping the log avoids underflow when k is large or x is very small.
    """
    _stage(k)
    if not math.isfinite(x) or x < 0:
        raise ValueError("x must be finite and nonnegative")
    if x == 0:
        return -math.inf
    if k == 1:
        return math.log(x)
    if x < k:
        lq = _small_log_tail(k, x)
        if lq < -36:
            # -log(1-q)/q = 1+O(q), below double relative precision here.
            return lq
        return math.log(-math.log1p(-math.exp(lq)))
    # Here the completed-lineage probability is substantial, so subtraction
    # x-log(E) does not suffer the tiny-tail cancellation of the x<k branch.
    return math.log(x - _log_polynomial(k, x))


def psi(k: int, x: float) -> float:
    """Exact-model natural lineage cumulative hazard, numerically evaluated.

    The mathematical function is x-log(sum_{j=0}^{k-1}x**j/j!). This
    implementation uses the positive Poisson tail when direct subtraction
    would lose precision. Values below float range may underflow to zero;
    log_psi retains their logarithm.
    """
    if k == 1:
        _stage(k)
        if not math.isfinite(x) or x < 0:
            raise ValueError("x must be finite and nonnegative")
        return x
    value = log_psi(k, x)
    return math.exp(value) if value != -math.inf else 0.0


def inverse_psi(k: int, h: float) -> float:
    """Unique positive x with psi(k,x)=h, by monotone bisection.

    Bracketing uses the proven psi(x)<=x**k/k! and psi(x)<=x. Bisection
    comparisons are on the log scale. No data-dependent optimization grid is
    used, and no claim of exact floating-point arithmetic is made.
    """
    _stage(k)
    _positive(h, "h")
    if k == 1:
        return h
    target = math.log(h)
    log_rare = (target + math.lgamma(k + 1)) / k
    rare = math.exp(log_rare)
    high = max(rare, h)
    while log_psi(k, high) < target:
        if high > sys.float_info.max / 2:
            high = sys.float_info.max
            if log_psi(k, high) < target:
                raise OverflowError("inverse is outside floating-point range")
            break
        high *= 2.0
    low = 0.0
    for _ in range(90):
        middle = low + (high - low) / 2.0
        if middle == low or middle == high:
            break
        if log_psi(k, middle) < target:
            low = middle
        else:
            high = middle
    return low + (high - low) / 2.0


def cumulative_hazard(risk: float) -> float:
    if not math.isfinite(risk) or not 0 < risk < 1:
        raise ValueError("risk must be strictly between zero and one; do not clip")
    return -math.log1p(-risk)


def rate(k: int, R: float, N: float, D: float) -> float:
    """Exact point inverse U_k(R,N,D), assuming effective clock v=D/N."""
    _positive(N, "N")
    _positive(D, "D")
    h = cumulative_hazard(R) / N
    return (N / D) * inverse_psi(k, h)


def rate_interval(
    k: int, Rlo: float, Rhi: float, Nlo: float, Nhi: float,
    Dlo: float, Dhi: float,
) -> tuple[float, float]:
    """Sharp inverse-rate range over a declared R,N,D rectangle.

    The clock v=D/N remains coupled to both N and D. This is a deterministic
    sensitivity rectangle, not an inferred sampling confidence interval.
    """
    if Rlo > Rhi or Nlo > Nhi or Dlo > Dhi:
        raise ValueError("interval lower endpoints must not exceed upper endpoints")
    return rate(k, Rlo, Nlo, Dhi), rate(k, Rhi, Nhi, Dlo)


def frontier(intervals: list[tuple[float, float]]) -> dict:
    """Sharp common-baseline multiplicative heterogeneity allowance.

    There exist u_i in intervals[i] and u0>0 with u0/rho<=u_i<=u0*rho
    iff rho**2 >= max(1,max(lower)/min(upper)). Generic interval intersection
    is classical; the note's model-specific work constructs these intervals.
    """
    if not intervals:
        raise ValueError("at least one rate interval is required")
    for lower, upper in intervals:
        _positive(lower, "lower rate")
        _positive(upper, "upper rate")
        if lower > upper:
            raise ValueError("rate interval is reversed")
    lower_index = max(range(len(intervals)), key=lambda i: intervals[i][0])
    upper_index = min(range(len(intervals)), key=lambda i: intervals[i][1])
    largest_lower = intervals[lower_index][0]
    smallest_upper = intervals[upper_index][1]
    required_ratio = max(1.0, largest_lower / smallest_upper)
    rho = math.sqrt(required_ratio)
    return {
        "rho": rho,
        "required_rate_ratio": required_ratio,
        "max_lower": largest_lower,
        "min_upper": smallest_upper,
        "lower_witness_index": lower_index,
        "upper_witness_index": upper_index,
        "baseline_lower_at_rho": largest_lower / rho,
        "baseline_upper_at_rho": rho * smallest_upper,
    }


def elasticity(k: int, x: float) -> float:
    """x*psi'(x)/psi(x), with the proven mathematical range [1,k]."""
    _stage(k)
    _positive(x, "x")
    if k == 1:
        return 1.0
    log_derivative = ((k - 1) * math.log(x) - math.lgamma(k)
                      - _log_polynomial(k, x))
    return math.exp(math.log(x) + log_derivative - log_psi(k, x))


def saturation(H: float) -> float:
    """d log(1-exp(-H))/d logH; use 0 at overflow-sized finite H."""
    _positive(H, "H")
    if H > 700:
        return H * math.exp(-H) / (-math.expm1(-H))
    return H / math.expm1(H)


def log_phi(k: int, z: float) -> float:
    """log Phi_k(z), Phi_k=k!*sum_{n>=0} z**n/(n+k)!, z>=0.

    At k>=2 this is the one-growing-intermediate rare-event shape in the
    note. For k=1 and z>0 it is a mathematical kernel, not expansion at an
    intermediate stage of a one-transition model.
    """
    _stage(k)
    if not math.isfinite(z) or z < 0:
        raise ValueError("z must be finite and nonnegative")
    if z == 0:
        return 0.0
    if z < k:
        term = z / (k + 1)
        tail = [term]
        for n in range(1, 100000):
            ratio = z / (k + n + 1)
            total = math.fsum(tail)
            if term * ratio / (1 - ratio) <= 2e-17 * max(total, 1e-300):
                return math.log1p(total)
            term *= ratio
            tail.append(term)
        raise ArithmeticError("Phi series did not converge")
    h = psi(k, z)
    log_tail_probability = math.log(-math.expm1(-h))
    return (math.lgamma(k + 1) + z - k * math.log(z)
            + log_tail_probability)


def normalized_sha256(path: str | Path) -> str:
    """SHA-256 after CRLF and lone-CR normalization to LF; retain all else."""
    payload = Path(path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def _decimal_psi(k: int, x: Decimal) -> Decimal:
    """Independent high-precision direct survivor-polynomial formula."""
    term = Decimal(1)
    polynomial = term
    for j in range(1, k):
        term *= x / j
        polynomial += term
    return x - polynomial.ln()


def _decimal_inverse(k: int, h: Decimal, known_upper: Decimal) -> Decimal:
    """Independent Decimal bisection, used only on synthetic check cases."""
    lo, hi = Decimal(0), known_upper
    assert _decimal_psi(k, hi) >= h
    for _ in range(120):
        mid = (lo + hi) / 2
        if _decimal_psi(k, mid) < h:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _decimal_phi(k: int, z: Decimal) -> Decimal:
    """Independent closed form, evaluated at 100 decimal digits."""
    if not z:
        return Decimal(1)
    polynomial = Decimal(0)
    for j in range(k):
        polynomial += z ** j / math.factorial(j)
    return Decimal(math.factorial(k)) * (z.exp() - polynomial) / z ** k


def _fraction_polynomial_checks() -> int:
    count = 0
    for k in range(1, 13):
        E = [Fraction(1, math.factorial(j)) for j in range(k)]
        derivative = [(j + 1) * E[j + 1] for j in range(k - 1)] + [Fraction(0)]
        difference = [a - b for a, b in zip(E, derivative)]
        assert difference == [Fraction(0)] * (k - 1) + [E[-1]]
        nonnegative = [(k - 1 - j) * E[j] for j in range(k)]
        assert all(c >= 0 for c in nonnegative)
        # Exact multiplication of e^-x by E through order k. All terms before
        # k cancel and the first completion coefficient is 1/k!.
        survivor = []
        for n in range(k + 1):
            survivor.append(sum(
                (E[j] * Fraction((-1) ** (n - j), math.factorial(n - j))
                 for j in range(min(n, k - 1) + 1)), Fraction(0)))
        assert survivor[0] == 1
        assert all(c == 0 for c in survivor[1:k])
        assert survivor[k] == Fraction(-1, math.factorial(k))
        count += 3
    # Independent centered-vs-pairwise OLS identity in exact arithmetic.
    x = list(map(Fraction, [0, 1, 3, 6]))
    y = list(map(Fraction, [2, -1, 4, 10]))
    mx, my = sum(x) / len(x), sum(y) / len(y)
    centered = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sum((a - mx) ** 2 for a in x)
    pairwise = sum((x[j] - x[i]) * (y[j] - y[i]) for i in range(4) for j in range(i + 1, 4)) / sum((x[j] - x[i]) ** 2 for i in range(4) for j in range(i + 1, 4))
    assert centered == pairwise
    # Generic interval-intersection frontier in squared, exact rational form.
    intervals = [(Fraction(1), Fraction(2)), (Fraction(4), Fraction(5)), (Fraction(7), Fraction(9))]
    L = max(a for a, _ in intervals)
    U = min(b for _, b in intervals)
    boundary = L / U
    assert boundary == Fraction(7, 2)
    assert L <= boundary * U
    assert L > (boundary - Fraction(1, 100)) * U
    count += 2
    # Independent forward-mean Taylor recursion for one growing intermediate.
    # Two integrations (terminal marking and immigration) must give the Phi
    # coefficient g**(n-k)/n!, regardless of the expanding stage's position.
    # This checks the limit kernel, not equality with finite-rate branching.
    degree = 16
    for k in range(2, 9):
        m = k - 1
        for growth in (Fraction(0), Fraction(1, 3), Fraction(2)):
            for expanded in range(m):
                coefficients = [[Fraction(0)] * (degree + 1) for _ in range(m)]
                coefficients[0][0] = Fraction(1)
                for p in range(degree):
                    for j in range(m):
                        own = growth * coefficients[j][p] if j == expanded else Fraction(0)
                        incoming = coefficients[j - 1][p] if j else Fraction(0)
                        coefficients[j][p + 1] = (own + incoming) / (p + 1)
                for n in range(2, degree + 1):
                    actual = coefficients[-1][n - 2] / ((n - 1) * n)
                    expected = (growth ** (n - k) / math.factorial(n)
                                if n >= k else Fraction(0))
                    assert actual == expected, (k, growth, expanded, n)
                count += 1
    return count


def exact_checks() -> dict:
    """Deterministic synthetic algebra and numerical checks; no empirical fits."""
    algebra = _fraction_polynomial_checks()
    max_relative_psi_error = 0.0
    max_relative_inverse_error = 0.0
    max_phi_log_error = 0.0
    decimal_cases = 0
    with localcontext() as ctx:
        ctx.prec = 100
        for k in range(1, 9):
            for raw in ("0.001", "0.125", "1", "4", "25"):
                x = Decimal(raw)
                h = _decimal_psi(k, x)
                actual_h = psi(k, float(x))
                relative = abs(actual_h / float(h) - 1)
                max_relative_psi_error = max(max_relative_psi_error, relative)
                assert relative < 2e-13, (k, raw, relative)
                reference_x = _decimal_inverse(k, h, 2 * x)
                actual_x = inverse_psi(k, float(h))
                inverse_error = abs(actual_x / float(reference_x) - 1)
                max_relative_inverse_error = max(max_relative_inverse_error, inverse_error)
                assert inverse_error < 5e-14, (k, raw, inverse_error)
                decimal_cases += 1
        for k in range(1, 9):
            for raw in ("0", "0.000001", "0.25", "2", "20"):
                z = Decimal(raw)
                expected = float(_decimal_phi(k, z).ln())
                discrepancy = abs(log_phi(k, float(z)) - expected)
                max_phi_log_error = max(max_phi_log_error, discrepancy)
                assert discrepancy < 2e-13, (k, raw, discrepancy)

    elasticity_cases = 0
    rectangle_cases = 0
    pair_cases = 0
    for k in range(1, 13):
        for x in (1e-8, 0.01, 0.5, 1.0, 5.0, 50.0):
            e = elasticity(k, x)
            assert 1 - 5e-12 <= e <= k + 5e-12, (k, x, e)
            elasticity_cases += 1
        bounds = rate_interval(k, 0.0001, 0.2, 100, 500, 1000, 3000)
        corners = [rate(k, r, n, d) for r, n, d in itertools.product(
            (0.0001, 0.2), (100, 500), (1000, 3000))]
        assert math.isclose(min(corners), bounds[0], rel_tol=2e-14)
        assert math.isclose(max(corners), bounds[1], rel_tol=2e-14)
        for r, n, d in itertools.product((0.0001, 0.03, 0.2), (100, 250, 500), (1000, 1800, 3000)):
            u = rate(k, r, n, d)
            assert bounds[0] * (1 - 2e-14) <= u <= bounds[1] * (1 + 2e-14)
        rectangle_cases += 1
        ratio = rate(k, 0.2, 100, 1000) / rate(k, 0.0001, 100, 1000)
        h_ratio = cumulative_hazard(0.2) / cumulative_hazard(0.0001)
        assert h_ratio ** (1 / k) * (1 - 2e-13) <= ratio <= h_ratio * (1 + 2e-13)
        pair_cases += 1

    f = frontier([(1, 2), (4, 5), (7, 9)])
    assert math.isclose(f["rho"] ** 2, 3.5, rel_tol=2e-15)
    assert math.isclose(f["baseline_lower_at_rho"], f["baseline_upper_at_rho"], rel_tol=2e-15)
    smaller = f["rho"] * (1 - 1e-8)
    assert 7 / smaller > 2 * smaller
    overlap = frontier([(1, 5), (2, 3)])
    assert overlap["rho"] == 1
    assert overlap["baseline_lower_at_rho"] == 2
    assert overlap["baseline_upper_at_rho"] == 3
    for N in (1, 1e6, 1e12):
        expected = cumulative_hazard(0.03) / 1e8
        assert math.isclose(rate(1, 0.03, N, 1e8), expected, rel_tol=2e-15)
    # Domain errors are substantive: reported R=1 is never silently clipped.
    for bad in (0, 1, -0.1, 1.1):
        try:
            rate(2, bad, 100, 1000)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid risk was accepted")
    # A small positive intensity can make the one-lineage risk astronomically
    # small; log-space inversion remains stable without a float tail value.
    x = inverse_psi(8, 1e-250)
    assert abs(log_psi(8, x) - math.log(1e-250)) < 3e-12
    return {
        "status": "passed; synthetic implementation checks only",
        "fraction_algebra_checks": algebra,
        "independent_decimal_inverse_cases": decimal_cases,
        "elasticity_cases": elasticity_cases,
        "coupled_rectangle_cases": rectangle_cases,
        "same_input_ratio_cases": pair_cases,
        "max_relative_psi_error": max_relative_psi_error,
        "max_relative_inverse_error": max_relative_inverse_error,
        "max_phi_log_error": max_phi_log_error,
    }


if __name__ == "__main__":
    print(json.dumps(exact_checks(), indent=2, sort_keys=True))
