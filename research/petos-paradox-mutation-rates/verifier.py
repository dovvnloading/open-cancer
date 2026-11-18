"""Nonclinical theorem helpers and independent numerical checks; stdlib only.

The general claims rest on the README proofs. These finite calculations check
their implementation; they establish no biological or clinical validity.
No empirical observations are fitted or read by this module.
"""

from __future__ import annotations

import json
import math
import sys
from decimal import Decimal, localcontext
from fractions import Fraction


def _check_k(k: int) -> None:
    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise ValueError("k must be a positive integer")


def _nonnegative(value: float, name: str, *, finite: bool = False) -> float:
    value = float(value)
    if math.isnan(value) or value < 0 or (finite and not math.isfinite(value)):
        raise ValueError(f"{name} must be nonnegative" + (" and finite" if finite else ""))
    return value


def _logaddexp(a: float, b: float) -> float:
    if a == -math.inf:
        return b
    if b == -math.inf:
        return a
    upper, lower = max(a, b), min(a, b)
    return upper + math.log1p(math.exp(lower - upper))


def _log_polynomial(k: int, x: float) -> float:
    """log(sum_{j=0}^{k-1} x**j/j!) without forming large powers."""
    if k == 1 or x == 0:
        return 0.0
    log_x = math.log(x)
    if x >= k:
        # Scale by the final, largest polynomial term.
        term = total = 1.0
        for j in range(k - 1, 0, -1):
            term *= j / x
            previous = total
            total += term
            if total == previous:
                break
        return (k - 1) * log_x - math.lgamma(k) + math.log(total)
    log_term = log_total = 0.0
    for j in range(1, k):
        log_term += log_x - math.log(j)
        log_total = _logaddexp(log_total, log_term)
    return log_total


def _log_completion_tail(k: int, x: float) -> float:
    """Poisson upper tail in log space, for 0 < x <= k.

    q=exp(-x)*x**k/k! * sum_{j>=0} x**j/((k+1)...(k+j)).
    Every series term decreases in this regime. The sum is positive and
    avoids cancellation in 1-exp(-x)*E_{k-1}(x).
    """
    term = total = 1.0
    j = k
    while True:
        j += 1
        term *= x / j
        previous = total
        total += term
        if total == previous:
            break
    return -x + k * math.log(x) - math.lgamma(k + 1) + math.log(total)


def log_psi(k: int, x: float) -> float:
    """log of the Erlang cumulative hazard, including underflow-sized tails.

    x=0 returns -inf; x=+inf returns +inf. For very small completion
    probability q, log(-log(1-q))=log(q)+O(q); the omitted correction below
    is less than half a double-precision ulp on the logarithmic scale.
    """
    _check_k(k)
    x = _nonnegative(x, "x")
    if x == 0:
        return -math.inf
    if math.isinf(x):
        return math.inf
    if k == 1:
        return math.log(x)
    if x <= k:
        log_q = _log_completion_tail(k, x)
        if log_q < -36.0:
            return log_q
        return math.log(-math.log1p(-math.exp(log_q)))
    return math.log(x - _log_polynomial(k, x))


def psi(k: int, x: float) -> float:
    """Erlang cumulative hazard x-log(E_{k-1}(x)), evaluated stably.

    Values below the floating-point range return zero; log_psi preserves
    their finite logarithms, and risk uses those logarithms directly.
    """
    _check_k(k)
    x = _nonnegative(x, "x")
    if x == 0 or math.isinf(x) or k == 1:
        return x
    if x <= k:
        return -math.log1p(-math.exp(_log_completion_tail(k, x)))
    return x - _log_polynomial(k, x)


def inverse_psi(k: int, h: float) -> float:
    """Unique nonnegative x with psi(k,x)=h, using log-space bisection."""
    _check_k(k)
    h = _nonnegative(h, "h")
    if h == 0 or math.isinf(h) or k == 1:
        return h
    log_h = math.log(h)
    rare_guess = (log_h + math.lgamma(k + 1)) / k
    # psi(x)<=min(x,x**k/k!), giving useful scales at both extremes.
    # Wider initial bounds and direct checks avoid relying on rounded guesses.
    lower = min(log_h, rare_guess) - 1.0
    upper = min(max(log_h, rare_guess) + 1.0, math.log(sys.float_info.max))
    while log_psi(k, math.exp(lower)) > log_h:
        lower -= 1.0
    while log_psi(k, math.exp(upper)) < log_h:
        if upper == math.log(sys.float_info.max):
            return math.inf
        upper = min(upper + 1.0, math.log(sys.float_info.max))
    for _ in range(160):
        midpoint = (lower + upper) / 2.0
        if midpoint == lower or midpoint == upper:
            break
        if log_psi(k, math.exp(midpoint)) < log_h:
            lower = midpoint
        else:
            upper = midpoint
    return math.exp((lower + upper) / 2.0)


def risk(k: int, x: float, n: float) -> float:
    """1-exp[-n*psi(k,x)], preserving tiny-tail times large-n products."""
    n = _nonnegative(n, "n", finite=True)
    if n == 0:
        raise ValueError("n must be positive")
    log_h = math.log(n) + log_psi(k, x)
    if log_h > math.log(746.0):
        return 1.0
    return -math.expm1(-math.exp(log_h))


def elasticity(k: int, x: float) -> float:
    """x*psi'(x)/psi(x); returns its limits k at zero and 1 at infinity."""
    _check_k(k)
    x = _nonnegative(x, "x")
    if x == 0:
        return float(k)
    if k == 1 or math.isinf(x):
        return 1.0
    log_eta = k * math.log(x) - math.lgamma(k) - _log_polynomial(k, x) - log_psi(k, x)
    return math.exp(log_eta)


def _power(value: float, k: int) -> float:
    try:
        return value ** k
    except OverflowError:
        return math.inf


def pair_envelope(k: int, b_low: float, b_high: float) -> tuple[float, float]:
    """Sharp closure bounds on psi(k,b*x)/psi(k,x), x>0, b in a box.

    Zero and +inf are permitted as one-sided limiting box endpoints.
    Endpoints at finite positive b!=1, k>1 need not be attained.
    """
    _check_k(k)
    b_low = _nonnegative(b_low, "b_low", finite=True)
    b_high = _nonnegative(b_high, "b_high")
    if b_high == 0 or b_low > b_high:
        raise ValueError("require 0 <= b_low <= b_high and b_high > 0")
    return min(b_low, _power(b_low, k)), max(b_high, _power(b_high, k))


def pair_modulation(
    k: int, y_low: float, y_high: float, b_low: float, b_high: float
) -> dict[str, float]:
    """Residual hazard multiplier bounds and the minimum fold infima.

    y=(H_i/H_j)/(N_i/N_j). Positive finite independent boxes give the
    proved sharp closure; zero or infinite bounds give conservative one-sided
    extensions. multiplier_* bound s_i/s_j in H_i=N_i*s_i*psi(k,c*B_i).
    fold_infimum is max(s_i/s_j,s_j/s_i); rho_infimum is its square root
    for the distinct convention s_i in [s0/rho,s0*rho]. No factor is
    automatically a biological survival probability.
    """
    y_low = _nonnegative(y_low, "y_low", finite=True)
    y_high = _nonnegative(y_high, "y_high")
    if y_low > y_high:
        raise ValueError("require y_low <= y_high")
    a, c = pair_envelope(k, b_low, b_high)
    multiplier_low = 0.0 if y_low == 0 or math.isinf(c) else y_low / c
    multiplier_high = math.inf if a == 0 or math.isinf(y_high) else y_high / a
    reciprocal_high = math.inf if multiplier_high == 0 else 1.0 / multiplier_high
    fold = max(1.0, multiplier_low, reciprocal_high)
    return {
        "multiplier_low": multiplier_low,
        "multiplier_high": multiplier_high,
        "fold_infimum": fold,
        "rho_infimum": math.sqrt(fold),
    }


def discrete_risk(k: int, p: float, divisions: int, n: float = 1.0) -> float:
    """Exact independent-lineage Bin(divisions,p) completion law.

    At most one sequential next-stage success occurs at a discrete division.
    This is a different clock convention from continuous Erlang intensities.
    """
    _check_k(k)
    p = _nonnegative(p, "p", finite=True)
    n = _nonnegative(n, "n", finite=True)
    if p > 1 or n == 0:
        raise ValueError("require 0 <= p <= 1 and n > 0")
    if isinstance(divisions, bool) or not isinstance(divisions, int) or divisions < 0:
        raise ValueError("divisions must be a nonnegative integer")
    if divisions < k or p == 0:
        return 0.0
    if p == 1:
        return 1.0
    log_p, log_failure = math.log(p), math.log1p(-p)

    def log_probability(j: int) -> float:
        return (
            math.lgamma(divisions + 1) - math.lgamma(j + 1)
            - math.lgamma(divisions - j + 1) + j * log_p
            + (divisions - j) * log_failure
        )

    if divisions * p <= k:
        log_q = -math.inf
        for j in range(k, divisions + 1):
            log_q = _logaddexp(log_q, log_probability(j))
        q = math.exp(log_q)
        if q >= 1:
            return 1.0
        return -math.expm1(n * math.log1p(-q))
    log_survival = -math.inf
    for j in range(k):
        log_survival = _logaddexp(log_survival, log_probability(j))
    return -math.expm1(n * log_survival)


def _decimal_reference(k: int, value: str, precision: int = 180) -> tuple[float, float, float]:
    """Independent direct survival calculation, never the Poisson-tail code."""
    with localcontext() as context:
        context.prec = precision
        x = Decimal(value)
        polynomial = sum((x ** j / Decimal(math.factorial(j)) for j in range(k)), Decimal(0))
        survival = (-x).exp() * polynomial
        h = -survival.ln()
        derivative = x ** (k - 1) / (Decimal(math.factorial(k - 1)) * polynomial)
        eta = x * derivative / h
        return float(h), float(h.ln()), float(eta)


def run_checks() -> dict[str, object]:
    """Deterministic, offline checks of helpers against independent formulas."""
    count = 0

    def check(condition: bool, label: str) -> None:
        nonlocal count
        if not condition:
            raise AssertionError(label)
        count += 1

    def close(a: float, b: float, label: str, *, relative: float = 3e-12, absolute: float = 0.0) -> None:
        check(math.isclose(a, b, rel_tol=relative, abs_tol=absolute), label)

    for k in (1, 2, 3, 8):
        for value in ("1e-12", "0.00001", "0.1", "1", "3", "12", "100"):
            x = float(value)
            h_dec, log_dec, eta_dec = _decimal_reference(k, value)
            close(psi(k, x), h_dec, f"Decimal hazard k={k}, x={value}")
            close(log_psi(k, x), log_dec, f"Decimal log hazard k={k}, x={value}", absolute=3e-12)
            close(elasticity(k, x), eta_dec, f"Decimal elasticity k={k}, x={value}")
            check(1 - 3e-12 <= elasticity(k, x) <= k + 3e-12, "elasticity bounds")
        for h in (1e-250, 1e-30, 1e-6, 0.2, 4.0, 100.0):
            x = inverse_psi(k, h)
            close(psi(k, x), h, f"inverse composition k={k}, h={h}")
        check(psi(k, 0) == 0 and log_psi(k, 0) == -math.inf, "zero clock boundary")
        check(inverse_psi(k, 0) == 0 and risk(k, 0, 10) == 0, "zero hazard boundary")
        check(elasticity(k, 0) == k and elasticity(k, math.inf) == 1, "elasticity limits")

    # psi underflows here, but large n restores a representable risk.
    _, log_dec, _ = _decimal_reference(8, "1e-60", precision=600)
    check(psi(8, 1e-60) == 0, "intentional floating underflow fixture")
    close(log_psi(8, 1e-60), log_dec, "underflow log hazard", absolute=3e-12)
    close(math.log(risk(8, 1e-60, 1e200)), math.log(1e200) + log_dec, "large n restores tail", absolute=3e-12)

    # Exact rational independent state recursion for discrete division trials.
    for k in (1, 2, 3, 5):
        for divisions in (0, 2, 5, 9):
            p = Fraction(1, 4)
            states = [Fraction(1)] + [Fraction(0)] * k
            for _ in range(divisions):
                next_states = [Fraction(0)] * (k + 1)
                for j, probability in enumerate(states):
                    if j == k:
                        next_states[k] += probability
                    else:
                        next_states[j] += probability * (1 - p)
                        next_states[j + 1] += probability * p
                states = next_states
            expected = 1 - (1 - states[k]) ** 3
            close(discrete_risk(k, float(p), divisions, 3), float(expected), "rational discrete recursion", absolute=1e-14)
    check(discrete_risk(3, 1, 2) == 0 and discrete_risk(3, 1, 3) == 1, "discrete division boundaries")

    for k in (2, 3, 8):
        for b in (0.3, 0.9, 1.1, 3.0):
            lower, upper = pair_envelope(k, b, b)
            for x in (0.01, 0.7, 4.0, 100.0):
                ratio = math.exp(log_psi(k, b * x) - log_psi(k, x))
                check(lower < ratio < upper, "strict finite pair envelope")
        r_small = risk(k, 0.4, 2.0)
        r_large = risk(k, 0.4, 20.0)
        close(r_large, -math.expm1(10 * math.log1p(-r_small)), "equal-burden target ratio")
        check(pair_envelope(k, 1, 1) == (1, 1), "equal burden envelope")

    check(pair_envelope(3, 0.5, 2) == (0.125, 8), "box crossing one")
    check(pair_envelope(1, 0.5, 2) == (0.5, 2), "one-stage box")
    modulation = pair_modulation(3, 0.02, 0.03, 0.9, 1.1)
    close(modulation["multiplier_low"], float(Fraction(20, 1331)), "multiplier lower direction")
    close(modulation["multiplier_high"], float(Fraction(30, 729)), "multiplier upper direction")
    close(modulation["fold_infimum"], float(Fraction(729, 30)), "pair fold infimum")
    close(modulation["rho_infimum"] ** 2, modulation["fold_infimum"], "symmetric versus pair allowance")
    check(pair_modulation(3, 0, math.inf, 0, math.inf) == {
        "multiplier_low": 0.0, "multiplier_high": math.inf,
        "fold_infimum": 1.0, "rho_infimum": 1.0,
    }, "unbounded uncertainty is vacuous")
    check(pair_modulation(3, 0, 0, 1, 1)["fold_infimum"] == math.inf, "zero upper risk cannot have positive modulation")

    # Construct equivalent endpoint laws at different k by inverse rates.
    for k in (1, 2, 4, 8):
        x = inverse_psi(k, 0.3 / 17)
        close(risk(k, x, 17), -math.expm1(-0.3), "endpoint rate/count alias")

    # For any proposed positive extra count, an equal-burden pair can be fit
    # when the common absolute clock and target normalization are free.
    for extra in (1, 2, 4):
        lower, upper = -20.0, 20.0
        target = math.log(0.1)
        for _ in range(100):
            middle = (lower + upper) / 2
            x = math.exp(middle)
            difference = log_psi(3 + extra, x) - log_psi(3, x)
            if difference < target:
                lower = middle
            else:
                upper = middle
        x = math.exp((lower + upper) / 2)
        close(psi(3 + extra, x) / psi(3, x), 0.1, "extra-count calibration alias")

    # Rank deficiency in the rare model: log B=u+v log M preserves the
    # combinations alpha+k*u and a+k*v under a different proposed k.
    alpha, a, k, u, v = (Fraction(2), Fraction(1), 3, Fraction(1, 2), Fraction(-1, 4))
    other_k = 5
    other_alpha = alpha + (k - other_k) * u
    other_a = a + (k - other_k) * v
    for log_m in (Fraction(-2), Fraction(0), Fraction(3)):
        log_b = u + v * log_m
        check(alpha + a * log_m + k * log_b == other_alpha + other_a * log_m + other_k * log_b, "rare design null space")

    # Synthetic three-species counterexample: every pair passes its strict
    # envelope, while the shared convex hazard curve cannot pass all points.
    # The README proves the general statement analytically; these exact
    # rational checks audit this displayed witness and convexity numerator.
    burdens = (Fraction(1), Fraction(2), Fraction(3))
    hazards = (Fraction(1), Fraction(3), Fraction(24, 5))
    check(burdens[0] + burdens[2] == 2 * burdens[1], "synthetic equal burden spacing")
    check(hazards[0] + hazards[2] < 2 * hazards[1], "synthetic violation of convexity")
    for k in range(2, 9):
        for i, j in ((0, 1), (0, 2), (1, 2)):
            b, y = burdens[j] / burdens[i], hazards[j] / hazards[i]
            check(b < y < b ** k, "synthetic pairwise compatibility")
        # (k-1)*E(x)-x*E'(x) has these exact coefficients. At k>1
        # its positive constant coefficient and nonnegative others imply
        # psi''(x)>0 at every x>0, as proved in the README.
        coefficients = tuple(Fraction(k - 1 - j, math.factorial(j)) for j in range(k))
        check(coefficients[0] > 0 and all(c >= 0 for c in coefficients), "convexity polynomial coefficients")

    invalid = (
        lambda: psi(0, 1), lambda: psi(2, -1), lambda: psi(2, math.nan),
        lambda: inverse_psi(2, -1), lambda: risk(2, 1, 0),
        lambda: pair_envelope(2, 2, 1), lambda: pair_envelope(2, 0, 0),
        lambda: pair_modulation(2, 2, 1, 1, 1),
    )
    for function in invalid:
        try:
            function()
        except ValueError:
            check(True, "invalid-domain refusal")
        else:
            check(False, "invalid-domain refusal")
    return {"checks": count, "status": "passed", "scope": "synthetic numerical checks; no biological validation"}


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2, sort_keys=True, allow_nan=False))
