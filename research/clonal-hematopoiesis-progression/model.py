"""Nonclinical branching-process identities and explicitly numerical helpers.

Only the Python standard library is used.  Rates are per year, horizons are
years, and ``s = birth - death``.  A mathematical identity is not biological
validation.  See THEOREMS.md for assumptions and proofs; multistage RK4 is an
approximation to the stated exact backward equations, not a proof.
"""

from __future__ import annotations

import math
from collections.abc import Sequence


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


def _count(n: int) -> int:
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("n must be a nonnegative integer")
    return n


def integrated_exposure(s: float, horizon: float) -> float:
    """A(s,T) = integral_0^T exp(s*t) dt, including s=0 and s<0."""
    s = _finite(s, "s")
    horizon = _nonnegative(horizon, "horizon")
    if s == 0 or horizon == 0:
        return horizon
    return math.expm1(s * horizon) / s


def exposure(s: float, horizon: float) -> float:
    """Short public alias for integrated_exposure."""
    return integrated_exposure(s, horizon)


def logistic_exposure(s: float, fraction: float, horizon: float) -> float:
    """B(s,x,T) = integral_0^T x(t)/x dt for logistic cell fraction x.

    ``fraction`` is the stem-cell fraction x=2*VAF, not VAF.  At x=0 this
    returns the continuous small-clone limit A; at x=1 it returns T.
    """
    s = _finite(s, "s")
    fraction = _nonnegative(fraction, "fraction")
    horizon = _nonnegative(horizon, "horizon")
    if fraction > 1:
        raise ValueError("fraction must be at most one")
    if s == 0 or horizon == 0 or fraction == 1:
        return horizon
    if fraction == 0:
        return integrated_exposure(s, horizon)
    z = s * horizon
    if z < 700:
        log_denom = math.log1p(fraction * math.expm1(z))
    else:
        # log((1-x) + x*exp(z)), without overflowing exp(z).
        log_denom = _logadd(math.log1p(-fraction), math.log(fraction) + z)
    return log_denom / (fraction * s)


def _logadd(x: float, y: float) -> float:
    if x == -math.inf:
        return y
    if y == -math.inf:
        return x
    largest = max(x, y)
    return largest + math.log1p(math.exp(min(x, y) - largest))


def _log_nonnegative(x: float) -> float:
    return math.log(x) if x > 0 else -math.inf


def _one_cell_terms(birth: float, death: float, hit: float,
                    horizon: float) -> tuple[float, float]:
    """Return exact single-cell cumulative hazard and instantaneous hazard."""
    birth = _nonnegative(birth, "birth")
    death = _nonnegative(death, "death")
    hit = _nonnegative(hit, "hit")
    horizon = _nonnegative(horizon, "horizon")
    if hit == 0:
        return 0.0, 0.0
    if horizon == 0:
        return 0.0, hit
    if birth == 0 and death == 0:
        return hit * horizon, hit
    s = birth - death
    delta = math.hypot(s - hit, 2.0 * math.sqrt(birth) * math.sqrt(hit))
    # A=delta-s+hit, B=delta+s-hit, C=delta-s-hit, D=delta+s+hit.
    # Rationalization avoids subtracting almost equal positive quantities.
    aa = (4.0 * birth * hit / (delta + s - hit)
          if s > hit else delta - s + hit)
    bb = (4.0 * birth * hit / aa if aa > 0 else 0.0)
    cc = (4.0 * death * hit / (delta + s + hit)
          if s + hit > 0 else delta - s - hit)
    dd = (4.0 * death * hit / cc if cc > 0 else delta + s + hit)
    z = delta * horizon
    log_den = _logadd(_log_nonnegative(aa), _log_nonnegative(bb) - z)
    log_num = _logadd(_log_nonnegative(cc), _log_nonnegative(dd) - z)
    # For small risks, log1p avoids cancellation in log_den-log_num.
    y = 2.0 * hit * (-math.expm1(-z)) / math.exp(log_den)
    if y < 0.75:
        cumulative = -math.log1p(-y)
    else:
        cumulative = log_den - log_num
    if death == 0:
        instantaneous = hit + birth * (-math.expm1(-cumulative))
    else:
        log_hazard = (math.log(4.0) + 2.0 * math.log(delta)
                      + math.log(hit) - z - log_num - log_den)
        instantaneous = math.exp(log_hazard)
    return cumulative, instantaneous


def cumulative_hazard(n: int, birth: float, death: float, hit: float,
                      horizon: float) -> float:
    """Exact H=-log P(no terminal hit) in the one-stage marked BD model."""
    n = _count(n)
    return n * _one_cell_terms(birth, death, hit, horizon)[0]


def log_survival_one(birth: float, death: float, hit: float,
                     horizon: float) -> float:
    """Log survival for one current cell under the exact marked BD model."""
    return -_one_cell_terms(birth, death, hit, horizon)[0]


def progression_probability(n: int, birth: float, death: float, hit: float,
                            horizon: float) -> float:
    """Exact probability of at least one terminal hit; not diagnosed AML."""
    return -math.expm1(-cumulative_hazard(n, birth, death, hit, horizon))


def instantaneous_hazard(n: int, birth: float, death: float, hit: float,
                         horizon: float) -> float:
    """Exact event-free risk-set hazard, distinct from cumulative hazard."""
    n = _count(n)
    return n * _one_cell_terms(birth, death, hit, horizon)[1]


def exposure_second_moment(birth: float, death: float, horizon: float) -> float:
    """E[(integral_0^T Z(t)dt)^2] for one unmarked ancestor."""
    birth = _nonnegative(birth, "birth")
    death = _nonnegative(death, "death")
    horizon = _nonnegative(horizon, "horizon")
    s = birth - death
    a = integrated_exposure(s, horizon)
    x = s * horizon
    if abs(x) < 0.01:
        # (sinh(x)/x - 1)/x^2 = sum_{k>=0} x^(2k)/(2k+3)!.
        term = 1.0 / 6.0
        series = term
        for k in range(1, 8):
            term *= x * x / ((2 * k + 2) * (2 * k + 3))
            series += term
        variance = 2.0 * (birth + death) * horizon**3 * math.exp(x) * series
    else:
        variance = (2.0 * (birth + death) / s**2
                    * (integrated_exposure(2 * s, horizon)
                       - horizon * math.exp(x)))
    return a * a + variance


def integrated_second_moment(birth: float, death: float, horizon: float) -> float:
    """Short public alias for exposure_second_moment."""
    return exposure_second_moment(birth, death, horizon)


def rare_event_bounds(n: int, birth: float, death: float, hit: float,
                      horizon: float) -> tuple[float, float]:
    """Rigorous lower/upper bounds on exact cumulative hazard (Theorem 2)."""
    n = _count(n)
    birth = _nonnegative(birth, "birth")
    death = _nonnegative(death, "death")
    hit = _nonnegative(hit, "hit")
    horizon = _nonnegative(horizon, "horizon")
    upper = n * hit * integrated_exposure(birth - death, horizon)
    error = n * hit * hit * exposure_second_moment(birth, death, horizon) / 2
    return max(0.0, upper - error), upper


def fitness_ratio_bounds(s_high: float, s_low: float,
                         horizon: float) -> tuple[float, float]:
    """Sharp common-kernel exposure-ratio bounds; no marginal HR claim."""
    s_high = _finite(s_high, "s_high")
    s_low = _finite(s_low, "s_low")
    horizon = _nonnegative(horizon, "horizon")
    if s_high < s_low:
        raise ValueError("s_high must be at least s_low")
    return 1.0, math.exp((s_high - s_low) * horizon)


def erlang_completion(stages: int, rate: float, delay: float) -> float:
    """CDF of the sum of `stages` iid exponential waiting times.

    stages=0 means immediate completion.  This is an optional latency model,
    not a fitted assumption about the number of additional hits in AML.
    """
    stages = _count(stages)
    rate = _nonnegative(rate, "rate")
    delay = _nonnegative(delay, "delay")
    if stages == 0:
        return 1.0
    if rate == 0 or delay == 0:
        return 0.0
    x = rate * delay
    # The lower Poisson tail is evaluated using a positive upper-tail series
    # when 1-exp(-x)*sum would cancel near zero.
    if x < stages + 1:
        log_term = -x + stages * math.log(x) - math.lgamma(stages + 1)
        term = math.exp(log_term)
        total = term
        j = stages
        for _ in range(10000):
            j += 1
            term *= x / j
            total += term
            if term <= max(total, 1e-300) * 1e-16:
                return min(1.0, total)
        raise ArithmeticError("Erlang series did not converge")
    terms = [-x + j * math.log(x) - math.lgamma(j + 1) for j in range(stages)]
    log_sum = -math.inf
    for term in terms:
        log_sum = _logadd(log_sum, term)
    return -math.expm1(log_sum)


def multistage_survival(n_by_stage: Sequence[int], births: Sequence[float],
                        deaths: Sequence[float], mutations: Sequence[float],
                        horizon: float, steps: int = 4096) -> float:
    """RK4 approximation to exact multitype backward survival equations.

    A mutation REPLACES one type-i cell by one type-(i+1) cell; the last
    mutation is a terminal endpoint. Birth produces two same-type cells.
    Current cells of every nonterminal stage must be specified, including
    hidden advanced cells if present. Compare step counts before relying on
    numerical output; this solver is not a symbolic proof.
    """
    n = [_count(value) for value in n_by_stage]
    b = [_nonnegative(value, "birth") for value in births]
    d = [_nonnegative(value, "death") for value in deaths]
    mu = [_nonnegative(value, "mutation") for value in mutations]
    horizon = _nonnegative(horizon, "horizon")
    k = len(n)
    if k == 0 or len(b) != k or len(d) != k or len(mu) != k:
        raise ValueError("all stage arrays must have the same positive length")
    if isinstance(steps, bool) or not isinstance(steps, int) or steps < 1:
        raise ValueError("steps must be a positive integer")
    if horizon == 0 or sum(n) == 0:
        return 1.0
    # Solve completion probabilities y=1-u to retain accuracy near zero.
    y = [0.0] * k
    dt = horizon / steps

    def rhs(z: Sequence[float]) -> list[float]:
        return [(b[i] - d[i] - mu[i]) * z[i] - b[i] * z[i]**2
                + mu[i] * (z[i + 1] if i + 1 < k else 1.0)
                for i in range(k)]

    for _ in range(steps):
        a1 = rhs(y)
        a2 = rhs([y[i] + dt * a1[i] / 2 for i in range(k)])
        a3 = rhs([y[i] + dt * a2[i] / 2 for i in range(k)])
        a4 = rhs([y[i] + dt * a3[i] for i in range(k)])
        y = [y[i] + dt * (a1[i] + 2*a2[i] + 2*a3[i] + a4[i]) / 6
             for i in range(k)]
        if any(not math.isfinite(value) or value < -1e-12 or value > 1 + 1e-12
               for value in y):
            raise ArithmeticError("RK4 left probability domain; increase steps")
    if any(y[i] >= 1 and n[i] > 0 for i in range(k)):
        return 0.0
    log_survival = sum(n[i] * math.log1p(-max(0.0, y[i]))
                       for i in range(k) if n[i] > 0)
    return math.exp(log_survival)
