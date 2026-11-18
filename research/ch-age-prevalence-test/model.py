"""Nonclinical, unreviewed CH prevalence mathematics; no medical advice.

Independent standard-library implementation of the equations in PROTOCOL.md.
The transported Watson density is an approximation, whereas the separately
provided integer birth--death identities are exact within their stated model.
No source implementation or participant records are copied here.
"""

from __future__ import annotations

import csv
import math
from collections.abc import Callable, Mapping
from functools import lru_cache
from pathlib import Path


EULER_GAMMA = 0.577215664901532860606512090082402431
FITNESS_LOWER = 0.0001
HERE = Path(__file__).resolve().parent
NUMERIC_FIELDS = (
    "shape", "scale_per_year", "s_max_per_year", "n_tau_years",
    "mu_per_cell_year",
)


def _finite(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _nonnegative(value: float, name: str) -> float:
    value = _finite(value, name)
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _positive(value: float, name: str) -> float:
    value = _finite(value, name)
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _integer(value: int, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def load_parameters(path: Path | str | None = None) -> dict[str, dict]:
    """Read the frozen CSV without fitting or altering a numerical parameter."""
    source = Path(path) if path is not None else HERE / "data" / "model_parameters.csv"
    result = {}
    with source.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            row = dict(row)
            for key in NUMERIC_FIELDS:
                row[key] = _positive(row[key], key)
            key = row["model_id"]
            if key in result:
                raise ValueError(f"duplicate model_id: {key}")
            result[key] = row
    if not result:
        raise ValueError("empty parameter table")
    return result


def e1(x: float) -> float:
    """E1(x)=integral_x^infinity exp(-u)/u du, for x>0.

    A convergent power series is used up to one; above one a modified-Lentz
    continued fraction for the upper incomplete gamma function is used.
    Very small positive results may underflow to zero in binary64.
    """
    x = float(x)
    if x == math.inf:
        return 0.0
    x = _positive(x, "x")
    if x <= 1.0:
        factorial_term = -x
        terms = [factorial_term]
        for k in range(2, 10000):
            factorial_term *= -x / k
            term = factorial_term / k
            terms.append(term)
            if abs(term) <= 2e-17:
                return -EULER_GAMMA - math.log(x) - math.fsum(terms)
        raise ArithmeticError("E1 power series did not converge")
    if x > 745.0:
        return 0.0
    tiny = 1e-300
    denominator = x + 1.0
    c = 1.0 / tiny
    d = 1.0 / denominator
    h = d
    for i in range(1, 10000):
        numerator = -float(i * i)
        denominator += 2.0
        d = numerator * d + denominator
        if abs(d) < tiny:
            d = math.copysign(tiny, d)
        c = denominator + numerator / c
        if abs(c) < tiny:
            c = math.copysign(tiny, c)
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) <= 4e-15:
            return math.exp(-x) * h
    raise ArithmeticError("E1 continued fraction did not converge")


def _e1_log_argument(log_x: float) -> float:
    """Evaluate E1(exp(log_x)) without overflowing or losing tiny arguments."""
    if log_x > math.log(745.0):
        return 0.0
    if log_x < -36.0:
        return -EULER_GAMMA - log_x + math.exp(log_x)
    return e1(math.exp(log_x))


def regularized_gamma_p(a: float, x: float) -> float:
    """Lower regularized incomplete gamma P(a,x), a>0 and x>=0."""
    a = _positive(a, "a")
    x = _nonnegative(x, "x")
    if x == 0:
        return 0.0
    log_prefactor = a * math.log(x) - x - math.lgamma(a)
    if x < a + 1.0:
        term = 1.0 / a
        total = term
        for n in range(1, 10000):
            term *= x / (a + n)
            total += term
            if term <= total * 2e-16:
                return min(1.0, math.exp(log_prefactor) * total)
        raise ArithmeticError("incomplete-gamma series did not converge")
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for n in range(1, 10000):
        an = -n * (n - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = math.copysign(tiny, d)
        c = b + an / c
        if abs(c) < tiny:
            c = math.copysign(tiny, c)
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) <= 4e-15:
            q = math.exp(log_prefactor) * h
            return min(1.0, max(0.0, 1.0 - q))
    raise ArithmeticError("incomplete-gamma continued fraction did not converge")


@lru_cache(maxsize=64)
def _log_normalization(shape: float, scale: float, maximum: float) -> float:
    a = 1.0 / shape
    x = math.exp(shape * math.log(maximum / scale))
    p = regularized_gamma_p(a, x)
    if p <= 0:
        raise ArithmeticError("normalization outside numerical range")
    return math.log(scale / shape) + math.lgamma(a) + math.log(p)


def dfe_normalization(shape: float | Mapping, scale_per_year: float | None = None,
                      s_max_per_year: float | None = None) -> float:
    """Integral_0^smax exp[-(s/scale)^shape] ds (not a rectangle rule).

    Accept either the three numbers or one load_parameters() model dictionary.
    """
    if isinstance(shape, Mapping):
        scale_per_year = shape["scale_per_year"]
        s_max_per_year = shape["s_max_per_year"]
        shape = shape["shape"]
    b = _positive(shape, "shape")
    d = _positive(scale_per_year, "scale_per_year")
    maximum = _positive(s_max_per_year, "s_max_per_year")
    return math.exp(_log_normalization(b, d, maximum))


def dfe_density(params: Mapping, s: float) -> float:
    """Normalized frozen DFE density, including its source support."""
    s = _nonnegative(s, "s")
    maximum = float(params["s_max_per_year"])
    if s > maximum:
        return 0.0
    shape, scale = float(params["shape"]), float(params["scale_per_year"])
    exponent = (s / scale) ** shape if s else 0.0
    return math.exp(-exponent - _log_normalization(shape, scale, maximum))


@lru_cache(maxsize=16)
def _legendre_rule(order: int) -> tuple[tuple[float, float], ...]:
    """Construct Gauss--Legendre nodes/weights independently by Newton's rule."""
    _integer(order, "order", 2)
    points = []
    for i in range((order + 1) // 2):
        z = math.cos(math.pi * (i + 0.75) / (order + 0.5))
        for _ in range(100):
            previous, current = 1.0, z
            for n in range(2, order + 1):
                previous, current = current, ((2*n-1)*z*current - (n-1)*previous)/n
            derivative = order * (z * current - previous) / (z*z - 1.0)
            new_z = z - current / derivative
            if abs(new_z - z) <= 2e-16:
                z = new_z
                break
            z = new_z
        else:
            raise ArithmeticError("Legendre root did not converge")
        # Reevaluate derivative at the final root before calculating its weight.
        previous, current = 1.0, z
        for n in range(2, order + 1):
            previous, current = current, ((2*n-1)*z*current - (n-1)*previous)/n
        derivative = order * (z * current - previous) / (z*z - 1.0)
        weight = 2.0 / ((1.0 - z*z) * derivative*derivative)
        points.append((-z, weight))
        if order % 2 == 0 or i < order // 2:
            points.append((z, weight))
    return tuple(sorted(points))


def gauss_integral(function: Callable[[float], float], lower: float, upper: float,
                   order: int = 64) -> float:
    """One Gauss--Legendre quadrature; used for independent refinement checks."""
    half = (upper - lower) / 2.0
    middle = (upper + lower) / 2.0
    return half * math.fsum(w * function(middle + half*x)
                           for x, w in _legendre_rule(order))


def integrate(function: Callable[[float], float], lower: float, upper: float,
              abs_tol: float = 1e-12, rel_tol: float = 1e-10) -> float:
    """Refined Gauss quadrature, requiring agreement of successive orders.

    This is a checked floating-point numerical approximation, not an
    outward-rounded interval certificate for arbitrary input functions.
    """
    lower, upper = _finite(lower, "lower"), _finite(upper, "upper")
    if lower == upper:
        return 0.0
    old = gauss_integral(function, lower, upper, 32)
    for order in (64, 128, 256):
        new = gauss_integral(function, lower, upper, order)
        if abs(new - old) <= abs_tol + rel_tol * abs(new):
            return new
        old = new
    raise ArithmeticError("quadrature did not agree through order 256")


def simpson_integral(function: Callable[[float], float], lower: float, upper: float,
                     intervals: int = 4096) -> float:
    """Separate composite-Simpson algorithm for numerical cross-checks."""
    _integer(intervals, "intervals", 2)
    if intervals % 2:
        raise ValueError("Simpson intervals must be even")
    step = (upper - lower) / intervals
    odd = math.fsum(function(lower + i*step) for i in range(1, intervals, 2))
    even = math.fsum(function(lower + i*step) for i in range(2, intervals, 2))
    return step * (function(lower) + function(upper) + 4*odd + 2*even) / 3.0


def _log_phi(s: float, age: float, n_tau: float) -> float:
    if s == 0:
        return math.log(age / (2.0 * n_tau))
    st = s * age
    log_increment = (math.log(math.expm1(st)) if st < 50.0
                     else st + math.log1p(-math.exp(-st)))
    return log_increment - math.log(2.0 * n_tau * s)


@lru_cache(maxsize=16384)
def _expected_cached(shape: float, scale: float, maximum: float, n_tau: float,
                     mu: float, age: float, v: float, v_upper: float) -> float:
    if age == 0 or mu == 0 or maximum <= FITNESS_LOWER:
        return 0.0
    log_z = _log_normalization(shape, scale, maximum)
    log_x = math.log(v / (1.0 - 2.0*v))
    log_x_upper = (math.log(v_upper / (1.0 - 2.0*v_upper))
                   if v_upper < 0.5 else None)

    def integrand(log_s: float) -> float:
        s = math.exp(log_s)
        log_phi = _log_phi(s, age, n_tau)
        tail = _e1_log_argument(log_x - log_phi)
        if log_x_upper is not None:
            tail -= _e1_log_argument(log_x_upper - log_phi)
        weight = math.exp(log_s - math.exp(shape * (log_s - math.log(scale))) - log_z)
        return weight * max(0.0, tail)

    factor = 2.0 * n_tau * mu
    integral = integrate(integrand, math.log(FITNESS_LOWER), math.log(maximum),
                         abs_tol=1e-12 / max(1.0, factor), rel_tol=1e-10)
    return factor * integral


def expected_clones(params: Mapping, t: float, v: float = 0.02,
                    n_tau: float | None = None, v_upper: float = 0.5) -> float:
    """Frozen Watson-density detectable-clone intensity (Theorem 1).

    This numerical transport uses the exact frozen DFE and source supply.
    It is not the exact integer birth--death model supplied separately below.
    """
    age = _nonnegative(t, "t")
    v = _positive(v, "v")
    v_upper = _positive(v_upper, "v_upper")
    if not v < v_upper <= 0.5:
        raise ValueError("require 0 < v < v_upper <= 0.5")
    a = _positive(params["n_tau_years"] if n_tau is None else n_tau, "n_tau")
    values = [_positive(params[k], k)
              for k in ("shape", "scale_per_year", "s_max_per_year")]
    mu = _nonnegative(params["mu_per_cell_year"], "mu_per_cell_year")
    return _expected_cached(*values, a, mu, age, v, v_upper)


def prevalence(params: Mapping, t: float, v: float = 0.02,
               n_tau: float | None = None, v_upper: float = 0.5) -> float:
    """At least one clone under the frozen density and Poisson assumptions."""
    return -math.expm1(-expected_clones(params, t, v, n_tau, v_upper))


def band_prediction(params: Mapping, lower_age: float, upper_age: float | None,
                    v: float = 0.02, n_tau: float | None = None,
                    v_upper: float = 0.5) -> tuple[float, float | None, float]:
    """Return (age lower bound, assumed uniform-age mean, age upper bound).

    An open upper age returns (p(lower), None, 1), per frozen protocol.
    """
    lower = _nonnegative(lower_age, "lower_age")
    p_lower = prevalence(params, lower, v, n_tau, v_upper)
    if upper_age is None:
        return p_lower, None, 1.0
    upper = _nonnegative(upper_age, "upper_age")
    if upper < lower:
        raise ValueError("upper_age must be >= lower_age")
    if upper == lower:
        return p_lower, p_lower, p_lower
    mean = integrate(lambda age: prevalence(params, age, v, n_tau, v_upper),
                     lower, upper, abs_tol=1e-11, rel_tol=1e-10) / (upper - lower)
    return p_lower, mean, prevalence(params, upper, v, n_tau, v_upper)


def _bd_logs(birth: float, death: float, t: float) -> tuple[float, float, float]:
    """Return log survival, log q, log(1-q), for birth>0 and t>0."""
    s = birth - death
    if s == 0:
        log_delta = -math.log1p(birth*t)
        q = birth*t / (1.0 + birth*t)
        log_q = (math.log1p(-math.exp(log_delta)) if log_delta < -0.1
                 else math.log(q))
        return log_delta, log_q, log_delta
    if s > 0:
        denominator = s + death * (-math.expm1(-s*t))
        log_survival = math.log(s / denominator)
        log_delta = math.log(s) - s*t - math.log(denominator)
        q = birth * (-math.expm1(-s*t)) / denominator
    else:
        rate = -s
        denominator = rate + birth * (-math.expm1(-rate*t))
        log_survival = math.log(rate) - rate*t - math.log(denominator)
        log_delta = math.log(rate / denominator)
        q = birth * (-math.expm1(-rate*t)) / denominator
    log_q = (math.log1p(-math.exp(log_delta)) if log_delta < -0.1
             else math.log(q))
    return log_survival, log_q, log_delta


def bd_tail(birth: float, death: float, t: float, threshold: int) -> float:
    """Exact P[Z(t)>=threshold] from one ancestor; cross-sectional, not ever."""
    birth, death = _nonnegative(birth, "birth"), _nonnegative(death, "death")
    t = _nonnegative(t, "t")
    m = _integer(threshold, "threshold", 1)
    if t == 0:
        return float(m == 1)
    if birth == 0:
        return math.exp(-death*t) if m == 1 else 0.0
    log_survival, log_q, _ = _bd_logs(birth, death, t)
    return math.exp(log_survival + (m-1)*log_q)


def bd_integrated_tail(birth: float, death: float, t: float, threshold: int) -> float:
    """Exact integrated BD tail evaluated by finite sum or positive series."""
    birth, death = _nonnegative(birth, "birth"), _nonnegative(death, "death")
    t = _nonnegative(t, "t")
    m = _integer(threshold, "threshold", 1)
    if t == 0:
        return 0.0
    if birth == 0:
        if m > 1:
            return 0.0
        return -math.expm1(-death*t) / death if death else t
    _, log_q, log_delta = _bd_logs(birth, death, t)
    if m == 1:
        return -log_delta / birth
    delta = math.exp(log_delta)
    if m * delta < 1.0 and m <= 100000:
        finite = math.fsum(math.exp(k*log_q)/k for k in range(1, m))
        return max(0.0, (-log_delta - finite) / birth)
    power = math.exp(m*log_q)
    if power == 0:
        return 0.0
    q = math.exp(log_q)
    total = power / m
    for k in range(m+1, m+1000001):
        power *= q
        total += power / k
        bound = power*q / ((k+1)*delta) if delta else math.inf
        if bound <= max(1e-15*total, 1e-300):
            return total / birth
    raise ArithmeticError("BD integrated-tail series exceeded its numerical budget")


def binomial_detection(depth: int, v: float, min_alt_reads: int) -> float:
    """Illustrative P[Bin(depth,v)>=min_alt_reads]; not a caller calibration."""
    depth = _integer(depth, "depth")
    minimum = _integer(min_alt_reads, "min_alt_reads")
    v = _nonnegative(v, "v")
    if v > 1:
        raise ValueError("v must be <= 1")
    if minimum == 0:
        return 1.0
    if minimum > depth or v == 0:
        return 0.0
    if v == 1:
        return 1.0
    log_p, log_q = math.log(v), math.log1p(-v)
    normalizer = math.lgamma(depth+1)
    result = math.fsum(math.exp(normalizer - math.lgamma(k+1)
                               - math.lgamma(depth-k+1) + k*log_p
                               + (depth-k)*log_q)
                       for k in range(minimum, depth+1))
    return min(1.0, max(0.0, result))


def multiplier_interval(lambda_at_a: float, lambda_at_b: float,
                         observed_lower: float, observed_upper: float
                         ) -> tuple[float, float]:
    """Sharp nonnegative common-supply multiplier interval (Theorem 5).

    A reversed interval denotes impossibility. Infinite upper bounds are
    intentional; this is a deterministic compatibility statement.
    """
    la = _nonnegative(lambda_at_a, "lambda_at_a")
    lb = _nonnegative(lambda_at_b, "lambda_at_b")
    low = _nonnegative(observed_lower, "observed_lower")
    high = _nonnegative(observed_upper, "observed_upper")
    if la > lb or not low <= high <= 1:
        raise ValueError("require la<=lb and 0<=observed_lower<=observed_upper<=1")
    if low == 1 or (lb == 0 and low > 0):
        return math.inf, -math.inf
    lower = -math.log1p(-low) / lb if lb else 0.0
    upper = (-math.log1p(-high) / la if la and high < 1 else math.inf)
    return lower, upper
