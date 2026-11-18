#!/usr/bin/env python3
"""Exact checks for the count-only resistance-origin theorems.

Nonclinical computational mathematics; synthetic inputs only. No treatment
advice. Algebra uses Fraction. Transcendental probability evaluations use
proved rational enclosures, not floating-point equality. Decimal is used only
for human-readable output. This is a verifier, not a fit to biological data.
"""

from fractions import Fraction as F
from math import comb, factorial, isqrt
from decimal import Decimal, localcontext
import argparse
import json


def kernel_checked(kernel):
    """Positive integer clone sizes, probabilities summing exactly to one."""
    if not kernel or any(type(k) is not int or k < 1 for k in kernel):
        raise ValueError("clone sizes must be positive integers")
    result = {k: F(v) for k, v in kernel.items() if F(v) != 0}
    if any(v < 0 for v in result.values()) or sum(result.values()) != 1:
        raise ValueError("clone probabilities must be nonnegative and sum to one")
    return result


def falling(n, k):
    if k < 0 or n < 0:
        raise ValueError("nonnegative integers required")
    return factorial(n) // factorial(n - k) if n >= k else 0


def thin_kernel(kernel, retention):
    """Binomially thin a clone, retaining the zero-size outcome."""
    kernel = kernel_checked(kernel)
    p = F(retention)
    if not 0 <= p <= 1:
        raise ValueError("retention must be in [0, 1]")
    result = {j: F(0) for j in range(max(kernel) + 1)}
    for k, mass in kernel.items():
        for j in range(k + 1):
            result[j] += mass * comb(k, j) * p**j * (1 - p)**(k - j)
    return result


def levy_coefficients(kernel, a, b, c, retention=1, exposure=1):
    """Return positive coefficients l_k in log F(z)=sum l_k(z**k-1)."""
    kernel = kernel_checked(kernel)
    a, b, c, exposure = map(F, (a, b, c, exposure))
    if min(a, b, c, exposure) < 0:
        raise ValueError("intensities and exposure must be nonnegative")
    thinned = thin_kernel(kernel, retention)
    result = {k: a * kernel.get(k, F(0)) + c * thinned[k]
              for k in range(1, max(kernel) + 1)}
    result[1] += b * exposure
    return {k: value for k, value in result.items() if value}


def log_pgf(levy, z):
    z = F(z)
    if not 0 <= z <= 1:
        raise ValueError("z must be in [0, 1]")
    return sum((F(v) * (z**k - 1) for k, v in levy.items()), F(0))


def cumulants(levy, order=3):
    """Factorial cumulants, derivatives of log F at z=1."""
    return tuple(sum((F(v) * falling(k, r) for k, v in levy.items()), F(0))
                 for r in range(1, order + 1))


def reduced_pmf(levy, limit):
    """Exact w_n with P(Y=n)=exp(-sum l_k)*w_n, 0<=n<=limit."""
    if limit < 0 or any(k < 1 or F(v) < 0 for k, v in levy.items()):
        raise ValueError("invalid compound-Poisson coefficients or limit")
    weights = [F(1)]
    for n in range(1, limit + 1):
        weights.append(sum((k * F(v) * weights[n - k]
                            for k, v in levy.items() if k <= n), F(0)) / n)
    return weights


def partition_weights(levy, limit):
    """Independent finite product expansion of exp(sum l_k z**k).

    Unlike reduced_pmf, this enumerates the number of jumps of each size.
    It checks the coefficient recurrence without calling that recurrence.
    """
    result = [F(1)] + [F(0)] * limit
    for k, value in sorted(levy.items()):
        update = [F(0)] * (limit + 1)
        for count in range(limit // k + 1):
            term = F(value)**count / factorial(count)
            for n in range(limit - k * count + 1):
                update[n + k * count] += result[n] * term
        result = update
    return result


def recover_levy_from_weights(weights):
    """Formal logarithm, independent of factorial-moment reconstruction."""
    if not weights or weights[0] != 1:
        raise ValueError("normalized zero coefficient must be one")
    result = {}
    for n in range(1, len(weights)):
        result[n] = F(weights[n]) - sum(
            (k * result[k] * weights[n - k] for k in range(1, n)), F(0)) / n
    return {k: value for k, value in result.items() if value}


def exp_minus_bounds(x, terms=40):
    """Rational lower/upper bounds for exp(-x), x>=0.

    Reduce to 0<=x/2**s<=1. Alternating Taylor sums with an odd/even
    last index bound the exponential; repeated squaring preserves order.
    No Decimal or binary floating-point arithmetic enters the certificate.
    """
    x = F(x)
    if x < 0 or type(terms) is not int or terms < 2:
        raise ValueError("x>=0 and integer terms>=2 required")
    scale = 0
    while x > 1:
        x /= 2
        scale += 1
    even = terms if terms % 2 == 0 else terms + 1
    total, term = F(1), F(1)
    lower = None
    for j in range(1, even + 1):
        term *= -x / j
        total += term
        if j == even - 1:
            lower = total
    upper = total
    for _ in range(scale):
        lower, upper = lower * lower, upper * upper
    return lower, upper


def pmf_bounds(levy, limit, terms=40):
    """Certified prefix intervals and a certified remaining-tail interval."""
    lo, hi = exp_minus_bounds(sum(levy.values(), F(0)), terms)
    weights = reduced_pmf(levy, limit)
    probabilities = [(lo * w, hi * w) for w in weights]
    total_weight = sum(weights, F(0))
    tail = (max(F(0), 1 - hi * total_weight),
            min(F(1), 1 - lo * total_weight))
    return probabilities, tail


def exact_sqrt(value):
    """Square root when the rational is an exact square (verification input)."""
    value = F(value)
    if value < 0:
        raise ValueError("negative square root")
    p, q = isqrt(value.numerator), isqrt(value.denominator)
    if p * p != value.numerator or q * q != value.denominator:
        raise ValueError("inverse is not rational; use the real-valued theorem formula")
    return F(p, q)


def _recover_given_retention(k0, k1, negative_log_zero, p):
    """Finish reconstruction once one holiday's occupancy p is known."""
    k0, k1 = tuple(map(F, k0)), tuple(map(F, k1))
    p = F(p)
    if k0[1] <= 0 or not 0 <= p < 1:
        raise ValueError("positive second factorial cumulant and p<1 required")
    if k0[0] <= k1[0]:
        raise ValueError("positive mean loss is necessary for positive clone sizes")
    x = (k1[1] / k0[1] - p * p) / (1 - p * p)
    if not 0 <= x < 1:
        raise ValueError("no positive reversible component under this model")
    v = (k0[0] - k1[0]) / ((1 - x) * (1 - p))
    b = k0[0] - v
    intensity = F(negative_log_zero) - b
    a, c = x * intensity, (1 - x) * intensity
    if min(a, b) < 0 or c <= 0:
        raise ValueError("observations do not yield admissible intensities")
    return {"a": a, "b": b, "c": c, "p": p}


def recover_two_arm(k0, k1, negative_log_zero):
    """Theorem 2: unknown shared kernel, finite third moment, two arms.

    On exact rational-model inputs the discriminant is an exact square.
    Empirical cumulants need not satisfy this model or yield rational roots.
    """
    k0, k1 = tuple(map(F, k0)), tuple(map(F, k1))
    if len(k0) < 3 or len(k1) < 3 or min(k0[1:3]) <= 0:
        raise ValueError("positive second and third baseline cumulants required")
    r2, r3 = k1[1] / k0[1], k1[2] / k0[2]
    if r2 >= 1:
        raise ValueError("absent loss or absent persisters")
    s = (1 - r3) / (1 - r2) - 1
    p = (s + exact_sqrt(s * s + 4 * s)) / 2
    return _recover_given_retention(k0, k1, negative_log_zero, p)


def recover_three_arm(k0, k1, k2, negative_log_zero):
    """Theorem 3: genuine two-way switching, only finite second moment."""
    k0, k1, k2 = [tuple(map(F, row)) for row in (k0, k1, k2)]
    drop = k0[0] - k1[0]
    if drop <= 0:
        raise ValueError("positive initial mean drop required")
    r = (k1[0] - k2[0]) / drop
    if not 0 < r < 1:
        raise ValueError("finite positive switching rate requires 0<r<1")
    d0, d1 = k0[1] - k1[1], k1[1] - k2[1]
    c2 = (r * d0 - d1) / (r * (1 - r) * (1 - r * r))
    c1 = (d0 - c2 * (1 - r * r)) / (1 - r)
    if c2 <= 0 or c1 < 0:
        raise ValueError("second cumulants contradict two-way common-kernel model")
    pi = c1 / (c1 + 2 * c2)
    p = pi + (1 - pi) * r
    result = _recover_given_retention(k0, k1, negative_log_zero, p)
    result.update({"r": r, "pi": pi})
    return result


def recover_kernel(baseline_levy, a, b, c):
    intensity = F(a) + F(c)
    if intensity <= 0:
        raise ValueError("nonzero clone intensity required")
    kernel = {k: (F(v) - (F(b) if k == 1 else 0)) / intensity
              for k, v in baseline_levy.items()}
    return kernel_checked(kernel)


def recover_singleton_four(means, exposures=(1, 2)):
    """One-way singleton model: arms (x,0),(x,H),(x,2H),(x',0)."""
    u0, u1, u2, other = map(F, means)
    x, y = map(F, exposures)
    if x == y or u0 <= u1:
        raise ValueError("distinct exposures and positive mean loss required")
    q = (u1 - u2) / (u0 - u1)
    if not 0 <= q < 1:
        raise ValueError("invalid one-way retention")
    c = (u0 - u1) / (1 - q)
    b = (other - u0) / (y - x)
    a = u0 - b * x - c
    if min(a, b) < 0:
        raise ValueError("negative intensity")
    return {"a": a, "b": b, "c": c, "q": q}


def doubleton_alias(baseline_levy, holiday_levy, q):
    """Theorem 4's continuum of indistinguishable two-arm doubleton models."""
    q = F(q)
    k0, k1 = cumulants(baseline_levy), cumulants(holiday_levy)
    d = k0[0] - k1[0]
    if not 0 <= q < 1 or d <= 0:
        raise ValueError("positive loss and q<1 required")
    e = (k0[1] - k1[1]) / (2 * d)
    w = e / (1 + q - e)
    intensity = k0[1] / (2 * w)
    c = d / ((1 - q) * (1 + w))
    a, b = intensity - c, k0[0] - intensity * (1 + w)
    if min(a, b, c) < 0 or not 0 < w < 1:
        raise ValueError("candidate lies outside the admissible alias interval")
    return {"a": a, "b": b, "c": c, "q": q,
            "kernel": {1: 1 - w, 2: w}}


def binomial_probabilities(n, p):
    p = F(p)
    if type(n) is not int or n < 0 or not 0 <= p <= 1:
        raise ValueError("invalid binomial parameters")
    return [F(comb(n, k)) * p**k * (1 - p)**(n - k) for k in range(n + 1)]


def best_zero_test(n, p_null, p_alternative, alpha=F(1, 20), randomized=False):
    """Exact increasing-binomial LR test, p_alternative>p_null.

    Randomized: most powerful among ALL tests of the n zero indicators.
    Nonrandomized: best among deterministic upper-threshold tests only.
    Reject for Z>boundary and with probability gamma for Z==boundary.
    """
    p0, p1, alpha = map(F, (p_null, p_alternative, alpha))
    if not 0 < p0 < p1 < 1 or not 0 < alpha < 1:
        raise ValueError("require 0<p_null<p_alternative<1 and 0<alpha<1")
    f0, f1 = binomial_probabilities(n, p0), binomial_probabilities(n, p1)
    size, power = F(0), F(0)
    for k in range(n, -1, -1):
        if size + f0[k] > alpha:
            gamma = (alpha - size) / f0[k] if randomized else F(0)
            return {"n": n, "boundary": k, "gamma": gamma,
                    "size": size + gamma * f0[k],
                    "power": power + gamma * f1[k]}
        size += f0[k]
        power += f1[k]
    raise AssertionError("alpha<1 ensures a boundary")


def minimum_zero_test(p_null, p_alternative, alpha=F(1, 20),
                      target_power=F(4, 5), randomized=False, max_n=1000):
    """Exhaustively certify the least n for this stated class of tests."""
    target_power = F(target_power)
    if not F(alpha) < target_power < 1:
        raise ValueError("target power must exceed alpha and be below one")
    for n in range(1, max_n + 1):
        candidate = best_zero_test(n, p_null, p_alternative, alpha, randomized)
        if candidate["power"] >= target_power:
            return candidate
    raise ValueError("no certified sample size through max_n")


def decimal_string(value, digits=20):
    with localcontext() as ctx:
        ctx.prec = digits
        value = F(value)
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def _jsonable(value):
    if isinstance(value, F):
        return {"exact": str(value), "decimal": decimal_string(value)}
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(v) for v in value]
    return value


def self_check():
    """Non-exhaustive rational instances plus exhaustive finite-power search."""
    kernel = {1: F(1, 4), 2: F(1, 4), 4: F(1, 2)}
    a, b, c = F(3, 5), F(1, 3), F(4, 5)
    base = levy_coefficients(kernel, a, b, c)
    k0 = cumulants(base)
    count = 0
    for q in (F(0), F(1, 10), F(2, 5), F(9, 10)):
        holiday = levy_coefficients(kernel, a, b, c, q)
        recovered = recover_two_arm(k0, cumulants(holiday), sum(base.values()))
        assert recovered == {"a": a, "b": b, "c": c, "p": q}
        assert recover_kernel(base, a, b, c) == kernel
        count += 1
    for pi in (F(0), F(1, 4), F(3, 4)):
        for r in (F(1, 5), F(1, 2), F(4, 5)):
            ks = [cumulants(levy_coefficients(kernel, a, b, c,
                              pi + (1 - pi) * r**j)) for j in range(3)]
            recovered = recover_three_arm(*ks, sum(base.values()))
            assert recovered == {"a": a, "b": b, "c": c,
                                 "p": pi + (1 - pi) * r, "pi": pi, "r": r}
            count += 1
    weights = reduced_pmf(base, 18)
    assert weights == partition_weights(base, 18)
    assert recover_levy_from_weights(weights) == base
    first = levy_coefficients({1: F(1, 2), 2: F(1, 2)}, 1, 1, 1)
    second = levy_coefficients({1: F(1, 2), 2: F(1, 2)}, 1, 1, 1, F(1, 2))
    alias = doubleton_alias(first, second, F(1, 3))
    for original, p in ((first, 1), (second, alias["q"])):
        assert levy_coefficients(alias["kernel"], alias["a"], alias["b"],
                                 alias["c"], p) == original
    deterministic = minimum_zero_test(F(1, 4), F(1, 2))
    randomized = minimum_zero_test(F(1, 4), F(1, 2), randomized=True)
    assert deterministic["n"] == 26 and deterministic["boundary"] == 10
    return {"label": "nonclinical computational mathematics; synthetic inputs",
            "rational_reconstructions": count,
            "coefficient_checks": "recurrence = independent product = formal logarithm",
            "two_arm_doubleton_alias": alias,
            "zero_test_deterministic_threshold_minimum": deterministic,
            "zero_test_randomized_LR_minimum": randomized}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true", help="run exact checks (default)")
    parser.parse_args()
    print(json.dumps(_jsonable(self_check()), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
