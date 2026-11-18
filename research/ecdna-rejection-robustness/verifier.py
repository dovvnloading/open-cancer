"""Independent exact checks for a nonclinical rejection-count sensitivity audit.

The mathematical routines deliberately import no analysis implementation.  The
CLI loads it only to compare independent results.  Finite checks do not replace
the numbered proofs in README.md.  Indicator correlation is not ecDNA allocation
dispersion.  No treatment advice.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb, gcd, lcm, prod
from typing import Iterable, Mapping


F = Fraction


def _integer(value: int, name: str) -> None:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer")


def _probability(value: Fraction, name: str) -> None:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be a Fraction")
    if not 0 <= value <= 1:
        raise ValueError(f"{name} must lie in [0, 1]")


def _count_inputs(n: int, k: int) -> None:
    _integer(n, "n")
    _integer(k, "k")
    if n < 1 or not 1 <= k <= n:
        raise ValueError("require n >= 1 and 1 <= k <= n")


def independent_binomial_tail(n: int, k: int, p: Fraction) -> Fraction:
    """Integer Pascal weights, independently of any probability recurrence."""
    _count_inputs(n, k)
    _probability(p, "p")
    a, b = p.numerator, p.denominator
    return F(sum(comb(n, s) * a**s * (b - a)**(n - s)
                 for s in range(k, n + 1)), b**n)


def independent_beta_pmf(n: int, p: Fraction, rho: Fraction) -> tuple[Fraction, ...]:
    """Rising-factorial beta integral; no neighboring-PMF recurrence."""
    _integer(n, "n")
    if n < 1:
        raise ValueError("n must be positive")
    _probability(p, "p")
    _probability(rho, "rho")
    if p == 0:
        return (F(1),) + (F(0),) * n
    if p == 1:
        return (F(0),) * n + (F(1),)
    if rho == 0:
        return tuple(F(comb(n, s)) * p**s * (1 - p)**(n - s)
                     for s in range(n + 1))
    if rho == 1:
        return (1 - p,) + (F(0),) * (n - 1) + (p,)
    concentration = (1 - rho) / rho
    a, b = p * concentration, (1 - p) * concentration
    ar, br = [F(1)], [F(1)]
    for j in range(n):
        ar.append(ar[-1] * (a + j))
        br.append(br[-1] * (b + j))
    denominator = prod(concentration + j for j in range(n))
    return tuple(comb(n, s) * ar[s] * br[n - s] / denominator
                 for s in range(n + 1))


def independent_beta_tail(n: int, k: int, p: Fraction, rho: Fraction) -> Fraction:
    _count_inputs(n, k)
    return sum(independent_beta_pmf(n, p, rho)[k:], F(0))


def independent_holm(values: Mapping[str, Fraction]) -> dict[str, Fraction]:
    """Closed Bonferroni testing, independent of the step-down implementation."""
    labels = tuple(sorted(values))
    adjusted = {label: F(0) for label in labels}
    for size in range(1, len(labels) + 1):
        for subset in combinations(labels, size):
            local = min(F(1), size * min(values[label] for label in subset))
            for label in subset:
                adjusted[label] = max(adjusted[label], local)
    return adjusted


def verify_moment_law(n: int, k: int, cap: Fraction, rho_cap: Fraction,
                      bound: Fraction,
                      law: Mapping[int, Fraction]) -> None:
    """Check an attaining primal count law, including the upper moment slack."""
    _count_inputs(n, k)
    _probability(cap, "cap")
    _probability(rho_cap, "rho_cap")
    _probability(bound, "bound")
    for s, weight in law.items():
        _integer(s, "support")
        if not 0 <= s <= n or weight < 0:
            raise AssertionError("invalid primal support or negative probability")
    assert sum(law.values(), F(0)) == 1
    assert sum(s * weight for s, weight in law.items()) == n * cap
    moment_cap = n * (n - 1) * (cap**2 + rho_cap * cap * (1 - cap))
    assert sum(s * (s - 1) * weight for s, weight in law.items()) <= moment_cap
    assert sum(weight for s, weight in law.items() if s >= k) == bound


def moment_dual_values(n: int, k: int, cap: Fraction,
                       rho_cap: Fraction) -> tuple[Fraction, ...]:
    """Compute and pointwise check every declared polynomial majorant."""
    _count_inputs(n, k)
    _probability(cap, "cap")
    _probability(rho_cap, "rho_cap")
    mu = n * cap
    second = n * (n - 1) * (cap**2 + rho_cap * cap * (1 - cap))
    values = [F(1), mu / k]
    for j in range(k - 1):
        denominator = (k - j) * (k - j - 1)
        for s in range(n + 1):
            assert F((s - j) * (s - j - 1), denominator) >= int(s >= k)
        values.append((second - 2 * j * mu + j * (j + 1)) / denominator)
    return tuple(values)


def small_count_lp(n: int, k: int, cap: Fraction,
                   rho_cap: Fraction) -> Fraction:
    """Independent complete vertex enumeration; intended for small n.

    With normalization and a fixed first moment, an extreme distribution has
    at most two support points if the second-moment inequality is inactive,
    and at most three if it is active.  All those vertices are enumerated.
    """
    _count_inputs(n, k)
    _probability(cap, "cap")
    _probability(rho_cap, "rho_cap")
    if n > 12:
        raise ValueError("small_count_lp is deliberately bounded to n <= 12")
    mu = n * cap
    factorial_cap = n * (n - 1) * (cap**2 + rho_cap * cap * (1 - cap))
    raw_cap = factorial_cap + mu
    best = F(-1)

    def consider(points: tuple[int, ...], weights: tuple[Fraction, ...]) -> None:
        nonlocal best
        if any(w < 0 for w in weights):
            return
        assert sum(weights, F(0)) == 1
        assert sum(s * w for s, w in zip(points, weights)) == mu
        if sum(s * (s - 1) * w for s, w in zip(points, weights)) > factorial_cap:
            return
        best = max(best, sum((w for s, w in zip(points, weights) if s >= k), F(0)))

    if mu.denominator == 1:
        consider((int(mu),), (F(1),))
    for a, b in combinations(range(n + 1), 2):
        wb = (mu - a) / (b - a)
        consider((a, b), (1 - wb, wb))
    for a, b, c in combinations(range(n + 1), 3):
        wa = (raw_cap - (b + c) * mu + b * c) / ((a - b) * (a - c))
        wb = (raw_cap - (a + c) * mu + a * c) / ((b - a) * (b - c))
        wc = 1 - wa - wb
        consider((a, b, c), (wa, wb, wc))
    assert best >= 0, "the binomial distribution guarantees feasibility"
    return best


def _trim(poly: Iterable[int]) -> tuple[int, ...]:
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result)


def _primitive(poly: Iterable[int]) -> tuple[int, ...]:
    values = _trim(poly)
    divisor = gcd(*values)
    return values if divisor == 0 else tuple(x // divisor for x in values)


def _multiply(a: Iterable[int], b: Iterable[int]) -> tuple[int, ...]:
    a, b = tuple(a), tuple(b)
    answer = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            answer[i + j] += x * y
    return _trim(answer)


def evaluate_polynomial(poly: Iterable[int], x: Fraction) -> Fraction:
    result = F(0)
    for coefficient in reversed(tuple(poly)):
        result = result * x + coefficient
    return result


@dataclass(frozen=True)
class TailPolynomial:
    """tail - threshold = scale * polynomial(rho) / denominator(rho)."""

    coefficients: tuple[int, ...]
    denominator: tuple[int, ...]
    scale: Fraction


@lru_cache(maxsize=None)
def beta_tail_polynomial(n: int, k: int, p: Fraction,
                         threshold: Fraction) -> TailPolynomial:
    """Clear a beta-binomial tail's positive denominator, including rho=1.

    The common factor (1-rho) is canceled symbolically before forming the
    polynomial, preventing an artificial root at complete dependence.
    """
    _count_inputs(n, k)
    _probability(p, "p")
    _probability(threshold, "threshold")
    a, b = p.numerator, p.denominator
    denominator = (1,)
    for j in range(n - 1):
        denominator = _multiply(denominator, (1, j))
    a_prefix, b_prefix = [(1,)], [(1,)]
    for j in range(1, n):
        a_prefix.append(_multiply(a_prefix[-1], (a, b * j - a)))
        b_prefix.append(_multiply(b_prefix[-1], (b - a, b * j - b + a)))
    lower = [0] * n
    for s in range(k):
        if s == 0:
            numerator = tuple((b - a) * x for x in b_prefix[n - 1])
        else:
            product = _multiply(_multiply(a_prefix[s - 1], b_prefix[n - s - 1]), (1, -1))
            factor = comb(n, s) * a * (b - a)
            numerator = tuple(factor * x for x in product)
        for j, coefficient in enumerate(numerator):
            lower[j] += coefficient
    raw = [-threshold.denominator * x for x in lower]
    factor = (threshold.denominator - threshold.numerator) * b**n
    for j, coefficient in enumerate(denominator):
        raw[j] += factor * coefficient
    raw = _trim(raw)
    divisor = gcd(*raw)
    if divisor == 0:
        return TailPolynomial((0,), denominator, F(1))
    return TailPolynomial(tuple(x // divisor for x in raw), denominator,
                          F(divisor, threshold.denominator * b**n))


def _power_to_bernstein(poly: tuple[int, ...]) -> tuple[int, ...]:
    degree = len(poly) - 1
    binomials = tuple(comb(degree, i) for i in range(degree + 1))
    scale = lcm(*binomials)
    result = []
    for j in range(degree + 1):
        result.append(sum(poly[i] * comb(j, i) * (scale // binomials[i])
                          for i in range(j + 1)))
    return _primitive(result)


def _sign_variations(values: Iterable[int]) -> int:
    signs = [1 if x > 0 else -1 for x in values if x]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def _bisect_bernstein(coefficients: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Exact de Casteljau midpoint subdivision, retaining a common scale."""
    degree = len(coefficients) - 1
    row = list(coefficients)
    left, right = [], []
    for depth in range(degree + 1):
        factor = 1 << (degree - depth)
        left.append(row[0] * factor)
        right.append(row[-1] * factor)
        row = [x + y for x, y in zip(row, row[1:])]
    # Do NOT trim trailing zero Bernstein coefficients: they encode endpoints.
    def primitive_bernstein(values: Iterable[int]) -> tuple[int, ...]:
        values = tuple(values)
        divisor = gcd(*values)
        return values if divisor == 0 else tuple(x // divisor for x in values)
    return primitive_bernstein(left), primitive_bernstein(reversed(right))


@dataclass(frozen=True)
class RootIsolation:
    roots: tuple[tuple[Fraction, Fraction], ...]
    zero_variation_leaves: int
    one_variation_leaves: int
    exact_roots: int
    splits: int


def isolate_unit_roots_certificate(poly: Iterable[int], width: Fraction = F(1, 10**13)) -> RootIsolation:
    """Isolate EVERY root on [0,1] using exact Bernstein/Descartes counts.

    A variation count of zero excludes roots in the open interval.  A count
    of one proves exactly one root, counting multiplicity.  Subdivide all
    larger counts until they resolve; fail explicitly at a work limit rather
    than mistake a scan for a complete root count.  Exact endpoint roots are
    retained once.  The result gives arithmetic brackets, not uncertainty.
    """
    poly = _primitive(poly)
    if poly == (0,):
        raise ValueError("the zero polynomial has no finite root isolation")
    if width <= 0:
        raise ValueError("width must be positive")
    if len(poly) == 1:
        return RootIsolation((), 1, 0, 0, 0)
    initial = _power_to_bernstein(poly)
    # Preserve degree when an endpoint happens to be a root.
    if len(initial) != len(poly):
        initial = initial + (0,) * (len(poly) - len(initial))
    roots: set[tuple[Fraction, Fraction]] = set()
    if initial[0] == 0:
        roots.add((F(0), F(0)))
    if initial[-1] == 0:
        roots.add((F(1), F(1)))
    pending = [(F(0), F(1), initial, 0)]
    processed = 0
    zero_leaves = one_leaves = splits = 0
    while pending:
        lo, hi, coefficients, depth = pending.pop()
        processed += 1
        if processed > 20000 or depth > 256:
            raise ArithmeticError("complete root isolation did not resolve within its explicit bound")
        variations = _sign_variations(coefficients)
        if variations == 0:
            zero_leaves += 1
            continue
        if variations == 1 and hi - lo <= width and coefficients[0] and coefficients[-1]:
            assert coefficients[0] * coefficients[-1] < 0
            one_leaves += 1
            roots.add((lo, hi))
            continue
        left, right = _bisect_bernstein(coefficients)
        splits += 1
        middle = (lo + hi) / 2
        if left[-1] == 0:
            roots.add((middle, middle))
        pending.append((middle, hi, right, depth + 1))
        pending.append((lo, middle, left, depth + 1))
    result = tuple(sorted(roots))
    for left, right in zip(result, result[1:]):
        assert left[1] < right[0]
    return RootIsolation(result, zero_leaves, one_leaves,
                         sum(lo == hi for lo, hi in result), splits)


def isolate_unit_roots(poly: Iterable[int], width: Fraction = F(1, 10**13)) -> tuple[tuple[Fraction, Fraction], ...]:
    """Convenience wrapper returning only the complete exact root brackets."""
    return isolate_unit_roots_certificate(poly, width).roots


@lru_cache(maxsize=None)
def _nominal_crossings_cached(observations: tuple[tuple[str, int, int], ...],
                              size: Fraction, level: Fraction,
                              width: Fraction) -> dict:
    raw_roots, polynomial_certificates = [], []
    for label, n, k in observations:
        for divisor in range(1, len(observations) + 1):
            polynomial = beta_tail_polynomial(n, k, size, level / divisor)
            certificate = isolate_unit_roots_certificate(polynomial.coefficients, width)
            polynomial_certificates.append({
                "line": label, "divisor": divisor,
                "degree": len(polynomial.coefficients) - 1,
                "root_count": len(certificate.roots),
                "zero_variation_leaves": certificate.zero_variation_leaves,
                "one_variation_leaves": certificate.one_variation_leaves,
                "exact_roots": certificate.exact_roots,
                "splits": certificate.splits,
            })
            for lo, hi in certificate.roots:
                # These inequalities independently certify every reported
                # nondegenerate endpoint bracket by the original beta law.
                left = independent_beta_tail(n, k, size, lo) - level / divisor
                right = independent_beta_tail(n, k, size, hi) - level / divisor
                assert left * right <= 0
                raw_roots.append({"line": label, "divisor": divisor, "lo": lo, "hi": hi})
    raw_roots.sort(key=lambda item: (item["lo"], item["hi"], item["line"], item["divisor"]))
    # The dataset's distinct crossings must be isolated separately.  If a
    # future dataset has coincident roots, fail instead of silently choosing
    # an ordering or overlooking a possible isolated equality rejection.
    for left, right in zip(raw_roots, raw_roots[1:]):
        if left["hi"] >= right["lo"]:
            raise ArithmeticError("overlapping raw root brackets need joint algebraic classification")

    def rejected(rho: Fraction) -> tuple[str, ...]:
        raw = {label: independent_beta_tail(n, k, size, rho)
               for label, n, k in observations}
        adjusted = independent_holm(raw)
        return tuple(sorted(label for label in raw if adjusted[label] <= level))

    regimes = []
    for index in range(len(raw_roots) + 1):
        lo = F(0) if index == 0 else raw_roots[index - 1]["hi"]
        hi = F(1) if index == len(raw_roots) else raw_roots[index]["lo"]
        regimes.append({"lo": lo, "hi": hi, "rejected": rejected((lo + hi) / 2)})
    transitions = []
    for index, root in enumerate(raw_roots):
        before, after = regimes[index]["rejected"], regimes[index + 1]["rejected"]
        if before != after:
            # All other raw-tail comparisons against every Holm threshold
            # are constant in this isolated bracket.  Set the crossing tail
            # EXACTLY to its algebraic-root value, rather than interpreting
            # a bracket midpoint as the root.  Closed testing then gives the
            # exact rejection set at equality without evaluating irrational
            # numbers or assuming a fixed Holm rank.
            at_root_raw = {label: independent_beta_tail(n, k, size, root["lo"])
                           for label, n, k in observations}
            at_root_raw[root["line"]] = level / root["divisor"]
            at_root_adjusted = independent_holm(at_root_raw)
            at_boundary = tuple(sorted(label for label in at_root_adjusted
                                       if at_root_adjusted[label] <= level))
            assert at_boundary == tuple(sorted(set(before) | set(after)))
            transitions.append({**root, "before": before, "after": after,
                                "at_boundary": at_boundary})
    return {"raw_roots": raw_roots, "polynomial_certificates": polynomial_certificates,
            "transitions": transitions, "regimes": regimes,
            "at_zero": rejected(F(0)), "at_one": rejected(F(1))}


def nominal_crossings(observations: Mapping[str, tuple[int, int]],
                      size: Fraction = F(1, 100), level: Fraction = F(1, 20),
                      width: Fraction = F(1, 10**13)) -> dict:
    """Complete nominal-size beta Holm transitions for a labelled count family.

    observations maps each line name to (number of mothers, rejections).
    The cached return value is read-only by convention; callers must not edit it.
    """
    _probability(size, "size")
    _probability(level, "level")
    rows = tuple(sorted((label, n, k) for label, (n, k) in observations.items()))
    for _, n, k in rows:
        _count_inputs(n, k)
    return _nominal_crossings_cached(rows, size, level, width)


def verify_analysis(analysis) -> None:
    """Bounded independent arithmetic checks, also callable by a test runner."""
    rows = analysis.load_data()
    observations = {row["line"]: (row["n"], row["k"]) for row in rows}
    nominal = rows[0]["nominal"]
    for row in rows:
        n, k = row["n"], row["k"]
        for cap in (F(0), nominal, F(1, 50), F(1, 5), F(1)):
            assert analysis.binomial_tail(n, k, cap) == independent_binomial_tail(n, k, cap)
            for rho in (F(0), F(1, 200), F(1, 100), F(1, 10), F(3, 4), F(1)):
                assert analysis.beta_tail(n, k, cap, rho) == independent_beta_tail(n, k, cap, rho)
                witness = analysis.sharp_certificate(n, k, cap, rho)
                verify_moment_law(n, k, cap, rho, witness["bound"], witness["law"])
                assert witness["bound"] == min(moment_dual_values(n, k, cap, rho))
                c0, c1, c2 = witness["dual"]
                assert c2 >= 0
                assert c0 + c1 * witness["mean"] + c2 * witness["second_cap"] == witness["bound"]
                assert all(c0 + c1 * s + c2 * s * (s - 1) >= int(s >= k)
                           for s in range(n + 1))
        for divisor in (1, 2, 3):
            threshold = F(1, 20 * divisor)
            polynomial = beta_tail_polynomial(n, k, nominal, threshold)
            for rho in (F(0), F(1, 5), F(1, 2), F(1)):
                exact_difference = polynomial.scale * evaluate_polynomial(polynomial.coefficients, rho)
                exact_difference /= evaluate_polynomial(polynomial.denominator, rho)
                assert exact_difference == independent_beta_tail(n, k, nominal, rho) - threshold
    for n in range(1, 7):
        for k in range(1, n + 1):
            for cap in (F(0), F(1, 10), F(1, 2), F(9, 10), F(1)):
                for rho in (F(0), F(1, 10), F(1)):
                    assert analysis.sharp_tail(n, k, cap, rho) == small_count_lp(n, k, cap, rho)
    nominal_crossings(observations, nominal)


def main() -> None:
    import importlib.util
    from pathlib import Path
    import sys

    spec = importlib.util.spec_from_file_location(
        "ecdna_rejection_robustness_analysis_standalone",
        Path(__file__).resolve().with_name("analysis.py"))
    analysis = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = analysis
    spec.loader.exec_module(analysis)
    verify_analysis(analysis)
    print("Independent rational probabilities, primal/dual moment witnesses, finite LPs, "
          "polynomial identities and complete root isolation passed. Nonclinical; unreviewed.")


if __name__ == "__main__":
    main()
