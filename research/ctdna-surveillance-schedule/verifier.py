"""Exact/high-precision checks for a NONCLINICAL ctDNA scheduling model.

Only Python's standard library is used. Fractions certify finite probability
calculations; Decimal calculations use 70 significant digits, not interval
arithmetic. No patient data, fitted clinical parameters, or monitoring advice.
See README.md for assumptions, proofs, provenance and falsification attempts.
"""

from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb
import json


PRECISION = 70
PI = Decimal(
    "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"
)


def decimal(value):
    """Convert without using a binary float as an exact decimal constant."""
    if isinstance(value, Fraction):
        return Decimal(value.numerator) / Decimal(value.denominator)
    return value if isinstance(value, Decimal) else Decimal(str(value))


def _positive_int(value, name, allow_zero=False):
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(name + " must be an integer")
    if value < (0 if allow_zero else 1):
        raise ValueError(name + " is out of range")


def _finite(*values):
    if any(not value.is_finite() for value in values):
        raise ValueError("parameters must be finite")


def binomial_detection(fraction, lod, molecules):
    """Exact gated Bin(m,f) probability, call K/m >= lod, true f >= lod."""
    _positive_int(molecules, "molecules")
    fraction, lod = Fraction(fraction), Fraction(lod)
    if not 0 <= fraction <= 1 or not 0 < lod < 1:
        raise ValueError("require 0 <= fraction <= 1 and 0 < lod < 1")
    if fraction < lod:
        return Fraction(0)
    threshold = (molecules * lod.numerator + lod.denominator - 1) // lod.denominator
    return sum((Fraction(comb(molecules, k)) * fraction**k
                * (1 - fraction)**(molecules - k)
                for k in range(threshold, molecules + 1)), Fraction(0))


def _grid_inputs(weights, rows):
    weights = tuple(map(Fraction, weights))
    rows = tuple(tuple(map(Fraction, row)) for row in rows)
    if not weights or len(weights) != len(rows) or sum(weights) != 1:
        raise ValueError("weights must sum to one and match scenario rows")
    if any(w < 0 for w in weights) or any(len(r) != len(rows[0]) for r in rows):
        raise ValueError("invalid scenario weights or row lengths")
    if any(not 0 <= p <= 1 for row in rows for p in row):
        raise ValueError("probabilities must be in [0,1]")
    return weights, rows


def _selection(selected, width):
    selected = tuple(selected)
    if len(set(selected)) != len(selected) or any(
        isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < width
        for i in selected
    ):
        raise ValueError("selected tests must be distinct valid indices")
    return selected


def _grid_value(weights, rows, selected):
    total = Fraction(0)
    for weight, row in zip(weights, rows):
        failure = Fraction(1)
        for j in selected:
            failure *= 1 - row[j]
        total += weight * (1 - failure)
    return total


def grid_value(weights, rows, selected):
    """Exact unconditional success; zero rows can encode unscored scenarios."""
    weights, rows = _grid_inputs(weights, rows)
    selected = _selection(selected, len(rows[0]))
    return _grid_value(weights, rows, selected)


def grid_marginal(weights, rows, selected, candidate):
    weights, rows = _grid_inputs(weights, rows)
    selected = _selection(selected, len(rows[0]))
    _selection((candidate,), len(rows[0]))
    if candidate in selected:
        return Fraction(0)
    return _grid_value(weights, rows, selected + (candidate,)) - _grid_value(
        weights, rows, selected
    )


def grid_optimum(weights, rows, n):
    """Exhaustive finite-grid optimum, exponential work; no continuum claim."""
    weights, rows = _grid_inputs(weights, rows)
    _positive_int(n, "n", allow_zero=True)
    if n > len(rows[0]):
        raise ValueError("budget exceeds grid size")
    best, value = (), Fraction(-1)
    for selected in combinations(range(len(rows[0])), n):
        candidate = _grid_value(weights, rows, selected)
        if candidate > value:
            best, value = selected, candidate
    return best, value


def grid_greedy(weights, rows, n):
    weights, rows = _grid_inputs(weights, rows)
    _positive_int(n, "n", allow_zero=True)
    if n > len(rows[0]):
        raise ValueError("budget exceeds grid size")
    selected = ()
    for _ in range(n):
        remaining = (j for j in range(len(rows[0])) if j not in selected)
        chosen = max(remaining, key=lambda j: _grid_value(weights, rows, selected + (j,)))
        selected += (chosen,)
    return selected, _grid_value(weights, rows, selected)


def normal_cdf(x):
    """Normal CDF by its integrated exponential series, with guard digits.

    For |x| >= 30 a tail smaller than 1e-197 is rounded to zero/one.
    The routine is high precision, not a certified outward-rounded enclosure.
    """
    x = decimal(x)
    _finite(x)
    if x >= 30:
        return Decimal(1)
    if x <= -30:
        return Decimal(0)
    with localcontext() as ctx:
        # The alternating sum loses O(x^2) decimal digits through cancellation.
        ctx.prec = PRECISION + 40 + int(abs(x) ** 2)
        if x == 0:
            return Decimal("0.5")
        term, total, k = x, x, 0
        threshold = Decimal(10) ** (-(ctx.prec - 10))
        while True:
            term *= -(x * x) / (2 * (k + 1))
            k += 1
            addend = term / (2 * k + 1)
            total += addend
            if abs(addend) < threshold:
                break
            if k > 10000:
                raise ArithmeticError("normal series did not converge")
        result = Decimal("0.5") + total / (2 * PI).sqrt()
        ctx.prec = PRECISION
        return +min(Decimal(1), max(Decimal(0), result))


@lru_cache(maxsize=1)
def _bernoulli_even():
    # Exact Akiyama--Tanigawa recurrence; odd B1 convention is irrelevant here.
    a, result = [], []
    for m in range(65):
        a.append(Fraction(1, m + 1))
        for j in range(m, 0, -1):
            a[j - 1] = j * (a[j - 1] - a[j])
        if m >= 2 and m % 2 == 0:
            result.append(a[0])
    return tuple(result)


@lru_cache(maxsize=128)
def _log_gamma(shape):
    with localcontext() as ctx:
        ctx.prec = PRECISION + 25
        x, correction = shape, Decimal(0)
        while x < 100:
            correction += x.ln()
            x += 1
        value = (x - Decimal("0.5")) * x.ln() - x + (2 * PI).ln() / 2
        for k, coefficient in enumerate(_bernoulli_even(), start=1):
            value += decimal(coefficient) / (2 * k * (2 * k - 1) * x ** (2 * k - 1))
        # At x>=100 the first omitted Stirling term is <1e-90.
        return +(value - correction)


def gamma_cdf(shape, x):
    """Regularized gamma CDF, positive real shape, UNIT scale, Decimal.

    A lower-tail series and an upper-tail continued fraction avoid cancellation.
    The numerical stopping criterion is 1e-85; this is not interval arithmetic.
    """
    shape, x = decimal(shape), decimal(x)
    _finite(shape, x)
    if shape <= 0:
        raise ValueError("shape must be positive")
    if x <= 0:
        return Decimal(0)
    with localcontext() as ctx:
        ctx.prec = PRECISION + 25
        eps = Decimal("1e-85")
        factor = (shape * x.ln() - x - _log_gamma(shape)).exp()
        if x < shape + 1:
            term = total = 1 / shape
            for k in range(1, 100001):
                term *= x / (shape + k)
                total += term
                if abs(term) < abs(total) * eps:
                    break
            else:
                raise ArithmeticError("gamma series did not converge")
            result = total * factor
        else:
            tiny = Decimal("1e-200")
            b = x + 1 - shape
            c, d = 1 / tiny, 1 / b
            h = d
            for k in range(1, 100001):
                an = -Decimal(k) * (Decimal(k) - shape)
                b += 2
                d = an * d + b
                if abs(d) < tiny:
                    d = tiny if d >= 0 else -tiny
                c = b + an / c
                if abs(c) < tiny:
                    c = tiny if c >= 0 else -tiny
                d = 1 / d
                change = d * c
                h *= change
                if abs(change - 1) < eps:
                    break
            else:
                raise ArithmeticError("gamma continued fraction did not converge")
            result = 1 - factor * h
        ctx.prec = PRECISION
        return +min(Decimal(1), max(Decimal(0), result))


def _schedule_parameters(n, a, b, horizon, endpoint):
    _positive_int(n, "n")
    a, b, horizon = map(decimal, (a, b, horizon))
    _finite(a, b, horizon)
    if not 0 < a < b or horizon <= 0:
        raise ValueError("require 0<a<b and positive horizon")
    floor = a / horizon
    if endpoint is not None:
        endpoint = decimal(endpoint)
        _finite(endpoint)
        if endpoint <= 0:
            raise ValueError("clinical endpoint must be positive")
        floor = max(floor, b / endpoint)
    return a, b, horizon, endpoint, floor


def _solution(n, a, b, lower, cdf, endpoint):
    q = b / a
    upper = lower * q**n
    schedule = tuple(a / (lower * q**j) for j in range(n - 1, -1, -1))
    value = cdf(upper) - cdf(lower)
    denominator = Decimal(1) if endpoint is None else 1 - cdf(b / endpoint)
    return {
        "schedule": schedule, "lower_rate": lower, "upper_rate": upper,
        "value": value, "denominator": denominator,
        "conditional_value": value / denominator if denominator else None,
    }


def optimal_lognormal(n, a, b, H, mu, sigma, endpoint=None):
    with localcontext() as ctx:
        ctx.prec = PRECISION
        a, b, _, endpoint, floor = _schedule_parameters(n, a, b, H, endpoint)
        mu, sigma = decimal(mu), decimal(sigma)
        _finite(mu, sigma)
        if sigma <= 0:
            raise ValueError("sigma must be positive")
        width = n * (b / a).ln()
        lower = max(floor, (mu - width / 2).exp())
        return _solution(n, a, b, lower,
                         lambda r: normal_cdf((r.ln() - mu) / sigma), endpoint)


def optimal_gamma(n, a, b, H, shape, scale, endpoint=None):
    with localcontext() as ctx:
        ctx.prec = PRECISION
        a, b, _, endpoint, floor = _schedule_parameters(n, a, b, H, endpoint)
        shape, scale = decimal(shape), decimal(scale)
        _finite(shape, scale)
        if shape <= 0 or scale <= 0:
            raise ValueError("shape and scale must be positive")
        q = b / a
        lower = max(floor, shape * scale * n * q.ln() / (q**n - 1))
        return _solution(n, a, b, lower,
                         lambda r: gamma_cdf(shape, r / scale), endpoint)


def ideal_value(schedule, a, b, cdf, endpoint=None):
    """Direct union integration, independent of the optimum's tiling formula."""
    with localcontext() as ctx:
        ctx.prec = PRECISION
        a, b = decimal(a), decimal(b)
        times = tuple(map(decimal, schedule))
        _finite(a, b, *times)
        if not 0 < a < b or any(t <= 0 for t in times):
            raise ValueError("require 0<a<b and positive times")
        floor = Decimal(0)
        if endpoint is not None:
            endpoint = decimal(endpoint)
            _finite(endpoint)
            if endpoint <= 0:
                raise ValueError("endpoint must be positive")
            floor = b / endpoint
        intervals = sorted((max(floor, a / t), b / t) for t in times
                           if b / t > floor)
        if not intervals:
            return Decimal(0)
        lo, hi = intervals[0]
        value = Decimal(0)
        for left, right in intervals[1:]:
            if left <= hi:
                hi = max(hi, right)
            else:
                value += cdf(hi) - cdf(lo)
                lo, hi = left, right
        return value + cdf(hi) - cdf(lo)


def robustness_bound(n, a, b, eta, density_bound, molecules, lod, delta,
                     finite_endpoint=False):
    """Uniform absolute-probability regret bound; not conditional miss rate.

    Finite endpoint adds the scored clinical-event boundary term. A conditional
    regret bound divides by the actual progression probability, then caps at 1.
    """
    _positive_int(n, "n")
    _positive_int(molecules, "molecules")
    with localcontext() as ctx:
        ctx.prec = PRECISION
        a, b, eta, density_bound, lod, delta = map(
            decimal, (a, b, eta, density_bound, lod, delta)
        )
        _finite(a, b, eta, density_bound, lod, delta)
        if not 0 <= eta < a < b or density_bound <= 0 or not 0 < lod < lod + delta < 1:
            raise ValueError("invalid bound parameters")
        shift = (((lod + delta) * (1 - lod)) / (lod * (1 - lod - delta))).ln()
        if shift >= b - a:
            raise ValueError("margin threshold must be below clinical burden")
        ha, hb = (a / (a - eta)).ln(), (b / (b - eta)).ln()
        burden = density_bound * (n * (ha + hb) + (hb if finite_endpoint else 0))
        noise = (n * density_bound * (1 + shift / (a - eta)).ln()
                 + (-2 * molecules * delta**2).exp())
        return {"burden_error": burden, "noise_error": noise,
                "regret_bound": min(Decimal(1), 2 * burden + noise)}


def uniform_onset_loss(schedule, H, W):
    """Exact Theorem 5 loss, including a missed onset's full window W."""
    H, W = Fraction(H), Fraction(W)
    schedule = tuple(map(Fraction, schedule))
    if H <= 0 or W < H or any(not 0 < t <= H for t in schedule):
        raise ValueError("require H>0, W>=H and times in (0,H]")
    if any(left >= right for left, right in zip(schedule, schedule[1:])):
        raise ValueError("test times must strictly increase")
    previous, total = Fraction(0), Fraction(0)
    for time in schedule:
        total += (time - previous)**2 / (2 * H)
        previous = time
    return total + (H - previous) * W / H


def thinning_bounds(mean, lower, upper, step):
    """Sharp mean-only/support-only bounds after random alternating thinning.

    The original detected cohort is kept and thinned misses receive zero lead.
    The bounds use a support interval; endpoint masses need not be positive.
    """
    mean, lower, upper, step = map(Fraction, (mean, lower, upper, step))
    if not 0 <= lower <= mean <= upper or step <= 0:
        raise ValueError("invalid nonnegative mean/support or thinning step")
    if lower == upper:
        value = mean - min(step, mean) / 2
        return value, value
    chord = (min(step, lower) + (mean - lower)
             * (min(step, upper) - min(step, lower)) / (upper - lower))
    return mean - min(step, mean) / 2, mean - chord / 2


def _demo():
    """Bounded synthetic reproducibility report; printed only, no file writes."""
    with localcontext() as ctx:
        ctx.prec = PRECISION
        report = {"status": "nonclinical computational mathematics; unreviewed"}
        for name, factory in (
            ("lognormal", lambda n: optimal_lognormal(n, 1, 2, 12, 0, 1)),
            ("gamma", lambda n: optimal_gamma(n, 1, 2, 12, 2, "0.5")),
        ):
            optimum, extra = factory(3), factory(4)
            report[name] = {**optimum, "reoptimized_extra_test_value": extra["value"] - optimum["value"]}
        report["finite_endpoint"] = optimal_lognormal(3, 1, 2, 1, 0, 1, endpoint=1)
        report["robustness"] = robustness_bound(
            3, 1, 2, "0.001", "0.4", 1000000, "0.01", "0.002"
        )
        report["reinert_nominal_thinning_bounds_months"] = thinning_bounds(
            Fraction(87, 10), Fraction(4, 5), Fraction(33, 2), 3
        )
        print(json.dumps(report, default=str, indent=2))


if __name__ == "__main__":
    _demo()
