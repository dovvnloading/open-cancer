"""Exact certificates for nonclinical computational mathematics of seeding.

Only Python's standard library is used. All mathematical inputs and calculations
are rational; floats are deliberately rejected. Transcendentals are enclosed by
proved rational Taylor bounds. These checks are not biological validation.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb, factorial
from typing import NamedTuple
import json


F = Fraction
DEFAULT_TOLERANCE = F(1, 10**35)


class RationalInterval(NamedTuple):
    lower: Fraction
    upper: Fraction


def _f(value):
    if isinstance(value, float):
        raise TypeError("use Fraction, integer, or exact decimal string, not float")
    return F(value)


def _integer(value, name, minimum=0):
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _probability(value, name):
    value = _f(value)
    if not 0 <= value <= 1:
        raise ValueError(f"{name} must lie in [0,1]")
    return value


def exp_neg_bounds(x, tolerance=DEFAULT_TOLERANCE):
    """Enclose exp(-x), x >= 0, with a rigorous rational interval.

    For S_n=sum_{j=0}^n x^j/j!, the omitted positive tail is at most
    t_{n+1}/(1-x/(n+2)) whenever n+2>x. Reciprocals reverse endpoints.
    The stopping tolerance bounds the *relative* interval width; no rounded
    floating-point calculation participates in either endpoint.
    """
    return _exp_neg_bounds(_f(x), _f(tolerance))


@lru_cache(maxsize=2048)
def _exp_neg_bounds(x, tolerance):
    if x < 0 or tolerance <= 0:
        raise ValueError("x must be nonnegative and tolerance positive")
    if not x:
        return RationalInterval(F(1), F(1))
    term = total = F(1)
    n = 0
    while True:
        if n + 2 > x:
            remainder = term * x / (n + 1) / (1 - x / (n + 2))
            if remainder <= tolerance * total:
                return RationalInterval(1 / (total + remainder), 1 / total)
        n += 1
        term *= x / n
        total += term


def poisson_pmf_bounds(mean, k, tolerance=DEFAULT_TOLERANCE):
    mean = _f(mean)
    _integer(k, "k")
    if mean < 0:
        raise ValueError("mean must be nonnegative")
    interval = exp_neg_bounds(mean, tolerance)
    multiplier = mean**k / factorial(k)
    return RationalInterval(interval.lower * multiplier, interval.upper * multiplier)


def poisson_cdf_bounds(mean, k, tolerance=DEFAULT_TOLERANCE):
    """Rational enclosure of P(Poisson(mean) <= k), including k < 0."""
    mean = _f(mean)
    if mean < 0 or not isinstance(k, int):
        raise ValueError("mean must be nonnegative and k integral")
    if k < 0:
        return RationalInterval(F(0), F(0))
    term = total = F(1)
    for j in range(1, k + 1):
        term *= mean / j
        total += term
    interval = exp_neg_bounds(mean, tolerance)
    return RationalInterval(
        min(F(1), interval.lower * total),
        min(F(1), interval.upper * total),
    )


def poisson_tail_bounds(mean, k, tolerance=DEFAULT_TOLERANCE):
    """Rational enclosure of P(Poisson(mean) >= k)."""
    interval = poisson_cdf_bounds(mean, k - 1, tolerance)
    return RationalInterval(1 - interval.upper, 1 - interval.lower)


def bd_survival_bounds(birth, death, age, tolerance=DEFAULT_TOLERANCE):
    """T1: P(Z(age)>0 | Z(0)=1), for linear birth-death at any rates."""
    birth, death, age = map(_f, (birth, death, age))
    if min(birth, death, age) < 0:
        raise ValueError("birth, death and age must be nonnegative")
    if birth == death:
        q = 1 / (1 + birth*age)
        return RationalInterval(q, q)
    if birth < death:
        r = death-birth
        e = exp_neg_bounds(r*age, tolerance)
        return RationalInterval(r*e.lower/(death-birth*e.lower),
                                r*e.upper/(death-birth*e.upper))
    e = exp_neg_bounds((birth - death) * age, tolerance)
    r = birth - death
    return RationalInterval(r / (birth - death * e.lower), r / (birth - death * e.upper))


def cox_count_pgf_bounds(integrated_intensity, z, tolerance=DEFAULT_TOLERANCE):
    """T1: conditional Poisson PGF, with the primary path held fixed."""
    integrated_intensity = _f(integrated_intensity)
    z = _probability(z, "z")
    if integrated_intensity < 0:
        raise ValueError("integrated_intensity must be nonnegative")
    return exp_neg_bounds(integrated_intensity * (1 - z), tolerance)


def piecewise_seed_mean(eta, segments):
    """Exact eta*integral N(t)dt for a supplied path; survival thinning = 1.

    segments contains (duration, population). This is the death-free-metastasis
    special case, not an approximation to age-dependent survival thinning.
    """
    eta = _f(eta)
    if eta < 0:
        raise ValueError("eta must be nonnegative")
    total = F(0)
    for duration, population in segments:
        duration = _f(duration)
        _integer(population, "population")
        if duration < 0:
            raise ValueError("duration must be nonnegative")
        total += duration * population
    return eta * total


def mutation_pattern_means(edges, mu, exposure=1):
    """T1: exact Poisson means for descendant-set passenger patterns.

    Each edge is (length, iterable_of_sample_labels). Distinct edges carrying
    the same pattern are aggregated. Labels and lengths describe the *sampled*
    genealogy; geographic migration marks do not enter this conditional law.
    """
    mu, exposure = _f(mu), _f(exposure)
    if mu < 0 or exposure < 0:
        raise ValueError("mu and exposure must be nonnegative")
    means = {}
    for length, labels in edges:
        length = _f(length)
        labels = frozenset(labels)
        if length < 0 or not labels:
            raise ValueError("edge length must be nonnegative and pattern nonempty")
        means[labels] = means.get(labels, F(0)) + exposure * mu * length
    return means


def mutation_pattern_pgf_bounds(pattern_means, leaf_z, tolerance=DEFAULT_TOLERANCE):
    """Joint PGF of per-leaf passenger counts, preserving shared mutations."""
    exponent = F(0)
    for labels, mean in pattern_means.items():
        mean = _f(mean)
        if mean < 0:
            raise ValueError("pattern means must be nonnegative")
        product = F(1)
        for label in labels:
            product *= _probability(leaf_z[label], "leaf PGF argument")
        exponent += mean * (1 - product)
    return exp_neg_bounds(exponent, tolerance)


def anchor_clock_rates(mu, B, c, h, ell):
    """T2: least-favorable clock rates, or an exact overlapping pair.

    Lag is additionally bounded by s, so a=s-lag is a nonnegative anchor.
    ell need not be <=s globally. The confounding witness uses lag 2h at
    s=c+h, which is feasible because c>=h.
    """
    mu, B, c, h, ell = map(_f, (mu, B, c, h, ell))
    if mu <= 0 or h <= 0 or ell < 0 or c - h < 0 or c + h > B:
        raise ValueError("require mu,h>0, ell>=0 and 0<=c-h<c+h<=B")
    early = mu * (B - c + h)
    if 2 * h <= ell:
        return {
            "status": "confounded", "common_rate": early,
            "early_pair": (c - h, F(0)), "late_pair": (c + h, 2 * h),
        }
    return {"status": "separated", "lambda_early": early,
            "lambda_late": mu * (B - c - h + ell)}


def threshold_errors(exposure, lambda_early, lambda_late, k,
                     tolerance=DEFAULT_TOLERANCE):
    """T2: errors for deterministic decision 'early iff X>=k'."""
    exposure, lambda_early, lambda_late = map(_f, (exposure, lambda_early, lambda_late))
    if min(exposure, lambda_early, lambda_late) < 0:
        raise ValueError("exposure and rates must be nonnegative")
    _integer(k, "k")
    return {
        "early": poisson_cdf_bounds(exposure * lambda_early, k - 1, tolerance),
        "late": poisson_tail_bounds(exposure * lambda_late, k, tolerance),
    }


def _certified_compare(evaluator, target):
    """Return (is <= target, interval); refine rather than round a boundary."""
    for digits in (35, 70, 140, 280):
        interval = evaluator(F(1, 10**digits))
        if interval.upper <= target:
            return True, interval
        if interval.lower > target:
            return False, interval
    raise ArithmeticError("comparison unresolved: request a less marginal target")


def _threshold_for_exposure(exposure, lambda_early, lambda_late, alpha):
    mean_late = exposure * lambda_late
    upper = max(1, mean_late.numerator // mean_late.denominator + 1)
    good, _ = _certified_compare(lambda tol: poisson_tail_bounds(mean_late, upper, tol), alpha)
    while not good:
        upper *= 2
        good, _ = _certified_compare(lambda tol: poisson_tail_bounds(mean_late, upper, tol), alpha)
    lower = 0  # P(X>=0)=1 > alpha.
    while upper - lower > 1:
        middle = (lower + upper) // 2
        good, _ = _certified_compare(lambda tol: poisson_tail_bounds(mean_late, middle, tol), alpha)
        if good:
            upper = middle
        else:
            lower = middle
    early_ok, early = _certified_compare(
        lambda tol: poisson_cdf_bounds(exposure * lambda_early, upper - 1, tol), alpha)
    _, late = _certified_compare(lambda tol: poisson_tail_bounds(mean_late, upper, tol), alpha)
    previous_late_ok, previous_late = _certified_compare(
        lambda tol: poisson_tail_bounds(mean_late, upper - 1, tol), alpha)
    assert not previous_late_ok
    return {"exposure": exposure, "threshold": upper, "feasible": early_ok,
            "early": early, "late": late, "previous_threshold_late": previous_late}


def minimal_threshold_exposure(lambda_early, lambda_late, alpha, max_exposure=100):
    """Smallest positive INTEGER exposure among deterministic integer thresholds.

    Every smaller exposure is checked. At each exposure the smallest threshold
    satisfying the late error constraint is certified; any larger threshold
    increases early error. This does not claim randomized minimax optimality.
    """
    lambda_early, lambda_late, alpha = map(_f, (lambda_early, lambda_late, alpha))
    _integer(max_exposure, "max_exposure", 1)
    if lambda_late < 0 or lambda_early < 0 or not 0 < alpha < F(1, 2):
        raise ValueError("rates must be nonnegative and 0<alpha<1/2")
    if lambda_early < lambda_late:
        raise ValueError("require lambda_early >= lambda_late for the early-if-large rule")
    if lambda_early == lambda_late:
        return {"status": "confounded", "rejected_exposures": ()}
    rejected = []
    for exposure in range(1, max_exposure + 1):
        certificate = _threshold_for_exposure(exposure, lambda_early, lambda_late, alpha)
        if certificate["feasible"]:
            return {"status": "attained", **certificate,
                    "rejected_exposures": tuple(rejected)}
        rejected.append(certificate)
    return {"status": "search_exhausted", "max_exposure": max_exposure,
            "rejected_exposures": tuple(rejected)}


def detection_probability(v, depth):
    v = _probability(v, "v")
    _integer(depth, "depth")
    return 1 - (1 - v)**depth


def _moments(m, s):
    m, s = _f(m), _f(s)
    if not 0 <= m <= 1 or not m * m <= s <= m:
        raise ValueError("require 0<=m<=1 and m^2<=s<=m")
    return m, s


def _clean_atoms(atoms):
    result = {}
    for value, weight in atoms:
        if weight:
            result[value] = result.get(value, F(0)) + weight
    return tuple(sorted(result.items()))


def sharp_miss_bounds(m, s, regions, c):
    """T3 sharp E[(1-cF)^R] bounds given E F=m and E F^2=s."""
    m, s = _moments(m, s)
    _integer(regions, "regions")
    c = _probability(c, "c")
    if m == 0 or m == 1:
        atoms = ((m, F(1)),)
        value = (1 - c * m)**regions
        return {"lower": value, "upper": value, "lower_atoms": atoms, "upper_atoms": atoms}
    upper_atoms = _clean_atoms(((F(0), 1 - m*m/s), (s/m, m*m/s)))
    denominator = 1 - 2*m + s
    lower_atoms = _clean_atoms((((m-s)/(1-m), (1-m)**2/denominator),
                                (F(1), (s-m*m)/denominator)))
    evaluate = lambda atoms: sum((w * (1-c*x)**regions for x, w in atoms), F(0))
    return {"lower": evaluate(lower_atoms), "upper": evaluate(upper_atoms),
            "lower_atoms": lower_atoms, "upper_atoms": upper_atoms}


def pilot_moments(a1, a2, c0):
    """Invert MISS summaries A1=P(miss one), A2=P(miss both regions)."""
    a1, a2 = _probability(a1, "a1"), _probability(a2, "a2")
    c0 = _probability(c0, "c0")
    if not c0:
        raise ValueError("a positive pilot detection probability is required")
    return _moments((1-a1)/c0, (a2-2*a1+1)/(c0*c0))


def _trim(poly):
    poly = list(map(_f, poly))
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return tuple(poly) or (F(0),)


def polynomial_add(a, b):
    result = [F(0)] * max(len(a), len(b))
    for i, value in enumerate(a):
        result[i] += value
    for i, value in enumerate(b):
        result[i] += value
    return _trim(result)


def polynomial_multiply(a, b):
    result = [F(0)] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        for j, right in enumerate(b):
            result[i+j] += left * right
    return _trim(result)


def polynomial_evaluate(poly, x):
    x, result = _f(x), F(0)
    for value in reversed(poly):
        result = result*x + value
    return result


def dual_certificate(m, s, regions, c, side="upper"):
    """T3 exact coefficient identity with a nonnegative remainder certificate.

    Coefficients are in ascending powers. For R>=3 the quotient is represented
    as sum_j positive_weights[j]*(1-c*x)^j, so positivity on [0,1]
    follows without numerical root finding or grid checking.
    """
    m, s = _moments(m, s)
    _integer(regions, "regions")
    c = _probability(c, "c")
    if side not in ("upper", "lower"):
        raise ValueError("side must be upper or lower")
    g = tuple(F(comb(regions, j)) * (-c)**j for j in range(regions+1))
    if regions <= 2:
        return {"dual": _trim(g), "factor": (F(0),), "quotient": (F(0),),
                "positive_weights": (), "identity": True}
    # At boundary moments a valid limiting interpolation node is sufficient.
    a = (s/m if m else F(0)) if side == "upper" else ((m-s)/(1-m) if m < 1 else F(1))
    ga = (1-c*a)**regions
    derivative = -regions*c*(1-c*a)**(regions-1)
    if side == "upper":
        if a == 0:
            quadratic = F(regions*(regions-1), 2)*c*c
            linear = -regions*c
        else:
            quadratic = (a*derivative-ga+1)/(a*a)
            linear = derivative-2*a*quadratic
        dual = _trim((F(1), linear, quadratic))
        factor = (F(0), a*a, -2*a, F(1))
        constants = (F(1), 1-c*a, 1-c*a)
        difference = polynomial_add(dual, tuple(-v for v in g))
    else:
        if a == 1:
            quadratic = F(regions*(regions-1), 2)*c*c*(1-c)**(regions-2)
        else:
            quadratic = ((1-c)**regions-ga-derivative*(1-a))/(1-a)**2
        dual = _trim((ga-a*derivative+quadratic*a*a,
                      derivative-2*a*quadratic, quadratic))
        factor = (a*a, -2*a-a*a, 1+2*a, F(-1))
        constants = (1-c*a, 1-c*a, 1-c)
        difference = polynomial_add(g, tuple(-v for v in dual))
    degree = regions-3
    homogeneous = [F(1)] + [F(0)]*degree
    for value in constants:
        for j in range(1, degree+1):
            homogeneous[j] += value*homogeneous[j-1]
    weights = tuple(c**3*homogeneous[degree-j] for j in range(degree+1))
    quotient = (F(0),)
    for j, weight in enumerate(weights):
        term = tuple(weight*comb(j, k)*(-c)**k for k in range(j+1))
        quotient = polynomial_add(quotient, term)
    identity = _trim(difference) == polynomial_multiply(factor, quotient)
    return {"dual": dual, "factor": _trim(factor), "quotient": quotient,
            "positive_weights": weights, "identity": identity, "node": a}


def _budget_search(bound, target, maximum, field):
    if bound(0) <= target:
        return {"status": "attained", field: 0, "upper": bound(0), "previous_upper": None}
    high = 1
    while high < maximum and bound(high) > target:
        high = min(maximum, 2*high)
    if bound(high) > target:
        return {"status": "search_exhausted", "maximum": maximum, "upper": bound(high)}
    low = 0
    while high-low > 1:
        middle = (high+low)//2
        if bound(middle) <= target:
            high = middle
        else:
            low = middle
    return {"status": "attained", field: high, "upper": bound(high),
            "previous_upper": bound(high-1)}


def minimum_regions(m, s, v, depth, alpha, max_regions=10000):
    m, s = _moments(m, s)
    c = detection_probability(v, depth)
    alpha = _probability(alpha, "alpha")
    _integer(max_regions, "max_regions", 1)
    bound = lambda r: sharp_miss_bounds(m, s, r, c)["upper"]
    if alpha == 1:
        return _budget_search(bound, alpha, max_regions, "regions")
    floor = 1-m*m/s if m and c else F(1)
    if alpha < floor or (alpha == floor and (not m or c*s/m != 1)):
        return {"status": "unreachable", "floor": floor}
    return {**_budget_search(bound, alpha, max_regions, "regions"), "floor": floor}


def minimum_depth(m, s, v, regions, alpha, max_depth=10000):
    m, s = _moments(m, s)
    v = _probability(v, "v")
    _integer(regions, "regions")
    _integer(max_depth, "max_depth", 1)
    alpha = _probability(alpha, "alpha")
    bound = lambda d: sharp_miss_bounds(m, s, regions, detection_probability(v, d))["upper"]
    if alpha == 1:
        return _budget_search(bound, alpha, max_depth, "depth")
    floor = sharp_miss_bounds(m, s, regions, F(1) if v else F(0))["upper"]
    if alpha < floor or (alpha == floor and v != 1):
        return {"status": "unreachable", "floor": floor}
    return {**_budget_search(bound, alpha, max_depth, "depth"), "floor": floor}


def pilot_upper_bound(h1, h2, c0, eps1, eps2, regions, c):
    """T4: robust bound from separate HIT pilot summaries and error radii.

    h1 estimates E(one hit)=c0*m; h2 estimates E(both hits)=c0^2*s.
    Inconsistent confidence rectangles cause abstention, never projection.
    """
    h1, h2, c0 = (_probability(x, name) for x, name in ((h1, "h1"), (h2, "h2"), (c0, "c0")))
    eps1, eps2 = _f(eps1), _f(eps2)
    _integer(regions, "regions")
    c = _probability(c, "c")
    if not c0 or eps1 <= 0 or eps2 <= 0:
        raise ValueError("require c0, eps1, eps2 positive")
    m_lower = max(F(0), (h1-eps1)/c0)
    s_upper = min(F(1), (h2+eps2)/(c0*c0))
    if m_lower > 1 or m_lower*m_lower > s_upper:
        return {"status": "inconsistent", "m_lower": m_lower, "s_upper": s_upper}
    s_star = min(m_lower, s_upper)
    upper = sharp_miss_bounds(m_lower, s_star, regions, c)["upper"]
    return {"status": "bounded", "m_lower": m_lower, "s_upper": s_upper,
            "s_star": s_star, "upper": upper}


def pilot_failure_bounds(n1, n2, eps1, eps2):
    """Enclose the Hoeffding UNION UPPER BOUND, not the unknown failure rate.

    The lower endpoint bounds the analytic quantity exp(-2*n1*eps1^2)
    + exp(-2*n2*eps2^2), clipped at one; it is not a lower confidence
    bound on the actual probability that the pilot confidence rectangle fails.
    """
    _integer(n1, "n1", 1)
    _integer(n2, "n2", 1)
    eps1, eps2 = _f(eps1), _f(eps2)
    if eps1 <= 0 or eps2 <= 0:
        raise ValueError("radii must be positive")
    first = exp_neg_bounds(2*n1*eps1*eps1)
    second = exp_neg_bounds(2*n2*eps2*eps2)
    return RationalInterval(min(F(1), first.lower+second.lower),
                            min(F(1), first.upper+second.upper))


def joint_error_bound(miss_upper, clock_error_upper, pilot_error_upper=0):
    """Union bound; independence between spatial and clock errors is not assumed."""
    values = [_probability(x, "error bound") for x in (miss_upper, clock_error_upper, pilot_error_upper)]
    return min(F(1), sum(values, F(0)))


def _choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def finite_pool_miss(n, invisible, regions, c):
    """T5 exact hypergeometric/depth miss probability, sampling without replacement."""
    _integer(n, "n", 2)
    _integer(invisible, "invisible")
    _integer(regions, "regions")
    c = _probability(c, "c")
    if invisible > n or regions > n:
        raise ValueError("invisible and regions must not exceed n")
    return sum((F(_choose(n-invisible, k)*_choose(invisible, regions-k), comb(n, regions))
                * (1-c)**k for k in range(regions+1)), F(0))


def finite_pool_extrema(n, m, s, regions, c):
    """T5 exact moment LP by all basic feasible supports of size at most 3.

    Here m=E(U/n) and s=E[U(U-1)/(n(n-1))] for the number U of invisible
    sites. These are *factorial* moments, so the T3 condition m^2<=s does
    not apply. Infeasibility is detected exactly.
    """
    _integer(n, "n", 2)
    m, s = _probability(m, "m"), _probability(s, "s")
    values = tuple(finite_pool_miss(n, u, regions, c) for u in range(n+1))
    first = n*m
    second = n*(n-1)*s + first
    candidates = []
    for u in range(n+1):
        if first == u and second == u*u:
            candidates.append(((u, F(1)),))
    for u, v in combinations(range(n+1), 2):
        w = (v-first)/(v-u)
        if 0 <= w <= 1 and w*u*u+(1-w)*v*v == second:
            candidates.append(_clean_atoms(((u, w), (v, 1-w))))
    for support in combinations(range(n+1), 3):
        atoms = []
        for i, u in enumerate(support):
            v, w = [value for j, value in enumerate(support) if j != i]
            weight = (second-(v+w)*first+v*w)/((u-v)*(u-w))
            atoms.append((u, weight))
        if all(weight >= 0 for _, weight in atoms):
            candidates.append(_clean_atoms(atoms))
    if not candidates:
        return {"status": "infeasible"}
    candidates = sorted(set(candidates))
    objective = lambda atoms: sum((weight*values[u] for u, weight in atoms), F(0))
    lower_atoms, upper_atoms = min(candidates, key=objective), max(candidates, key=objective)
    return {"status": "attained", "lower": objective(lower_atoms), "upper": objective(upper_atoms),
            "lower_atoms": lower_atoms, "upper_atoms": upper_atoms,
            "feasible_supports": len(candidates)}


def occupancy_all_seen(probabilities, regions):
    """T6 exact IE for independent CATEGORICAL draws, with optional null mass.

    Unlike the union bound, this formula requires that each draw reports at
    most one category and that categorical draws are independent.
    """
    probabilities = tuple(_probability(p, "category probability") for p in probabilities)
    _integer(regions, "regions")
    if sum(probabilities, F(0)) > 1:
        raise ValueError("category probabilities must sum to at most one")
    total = F(0)
    for size in range(len(probabilities)+1):
        for subset in combinations(probabilities, size):
            total += (-1)**size * (1-sum(subset, F(0)))**regions
    return total


def occupancy_miss_bound(entries, f0, c, regions):
    _integer(entries, "entries")
    _integer(regions, "regions")
    f0, c = _probability(f0, "f0"), _probability(c, "c")
    return min(F(1), entries*(1-f0*c)**regions)


def occupancy_budget(entries, f0, v, depth, alpha, max_regions=10000):
    _integer(entries, "entries")
    _integer(max_regions, "max_regions", 1)
    f0 = _probability(f0, "f0")
    c = detection_probability(v, depth)
    alpha = _probability(alpha, "alpha")
    bound = lambda r: occupancy_miss_bound(entries, f0, c, r)
    if entries and alpha < 1 and (not f0*c or (alpha == 0 and f0*c != 1)):
        return {"status": "unreachable", "floor": F(1) if not f0*c else F(0)}
    return _budget_search(bound, alpha, max_regions, "regions")


def binomial_upper_tail(n, k, p):
    """Exact P(Binomial(n,p)>=k), including zero/one and empty-trial cases."""
    _integer(n, "n")
    if not isinstance(k, int) or isinstance(k, bool):
        raise ValueError("k must be an integer")
    p = _probability(p, "p")
    if k <= 0:
        return F(1)
    if k > n:
        return F(0)
    return sum((comb(n, j)*p**j*(1-p)**(n-j) for j in range(k, n+1)), F(0))


def published_aggregate_check():
    """Rerunnable CONDITIONAL model check using Al Bakir et al. (2023).

    Their single-region downsampling misclassified 75/90 late-divergence cases
    as early. Under the additional null that these are independent witness-miss
    indicators, each with probability <=1/2, the upper tail is bounded by the
    binomial tail below. This is neither a published p-value nor validation of
    the mapping from divergence classification to witness failure.

    Treating 75/90 as an exact population miss probability is a separate,
    explicitly idealized calculation: Jensen gives A_R >= (5/6)^R for fresh
    conditionally independent regions. R>=17 is necessary, NOT sufficient,
    for a 5% miss guarantee in that experiment. Finite-pool downsampling and
    sampling uncertainty require separate treatment (Theorems 4 and 5).
    """
    n, observed = 90, 75
    null_p = F(1, 2)
    single_miss, target = F(observed, n), F(1, 20)
    regions = 0
    while single_miss**regions > target:
        regions += 1
    return {
        "n": n,
        "observed": observed,
        "null_p": null_p,
        "model_check_upper_tail": binomial_upper_tail(n, observed, null_p),
        "observed_single_region_miss": single_miss,
        "homogeneous_floor_R5percent": regions,
        "homogeneous_floor": single_miss**regions,
        "previous_homogeneous_floor": single_miss**(regions-1),
        "target": target,
        "source": "https://doi.org/10.1038/s41586-023-05729-x",
        "interpretation": "Conditional null check, not a published p-value or biological validation; 17 regions is necessary, not sufficient.",
    }


def run_checks():
    """Small exact self-checks; independent tests/simulation live separately."""
    certificates = 0
    for m, s in ((F(0), F(0)), (F(1), F(1)), (F(1, 2), F(1, 4)),
                 (F(1, 2), F(1, 3)), (F(1, 2), F(1, 2))):
        for regions in (0, 1, 2, 3, 4, 8, 16):
            for c in (F(0), F(1, 3), F(1)):
                bounds = sharp_miss_bounds(m, s, regions, c)
                assert 0 <= bounds["lower"] <= bounds["upper"] <= 1
                for side in ("upper", "lower"):
                    certificate = dual_certificate(m, s, regions, c, side)
                    assert certificate["identity"]
                    assert all(w >= 0 for w in certificate["positive_weights"])
                    certificates += 1
    exposure = minimal_threshold_exposure(F(4), F(2), F(1, 20), 30)
    assert exposure["status"] == "attained"
    assert all(not item["feasible"] for item in exposure["rejected_exposures"])
    regions = minimum_regions(F(1, 2), F(3, 10), F(1, 2), 4, F(1, 5))
    assert regions["status"] == "attained"
    pool = finite_pool_extrema(8, F(1, 2), F(1, 4), 4, F(3, 4))
    assert pool["status"] == "attained"
    assert occupancy_all_seen((F(1, 3), F(1, 3)), 2) == F(2, 9)
    assert anchor_clock_rates(1, 10, 5, 1, 2)["status"] == "confounded"
    assert bd_survival_bounds(2, 1, 0) == (F(1), F(1))
    aggregate = published_aggregate_check()
    assert aggregate["model_check_upper_tail"] == F(1774990439854189, 38685626227668133590597632)
    assert aggregate["homogeneous_floor_R5percent"] == 17
    assert aggregate["homogeneous_floor"] <= aggregate["target"] < aggregate["previous_homogeneous_floor"]
    return {"status": "passed", "scope": "nonclinical computational mathematics",
            "exact_dual_certificates": certificates,
            "clock_example": {"exposure": exposure["exposure"], "threshold": exposure["threshold"]},
            "region_example": {"regions": regions["regions"], "upper": str(regions["upper"])},
            "finite_pool_example": {"lower": str(pool["lower"]), "upper": str(pool["upper"])},
            "published_aggregate_check": {key: str(value) if isinstance(value, F) else value
                                          for key, value in aggregate.items()},
            "note": "Exact certificates and synthetic checks are not biological validation."}


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2, sort_keys=True))
