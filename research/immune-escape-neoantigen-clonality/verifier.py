#!/usr/bin/env python3
"""Exact certificates for the accompanying nonclinical branching mathematics.

All certificate decisions use integers and fractions.Fraction. Decimal is used
only to display already certified rational intervals. There are no third-party
dependencies, fitted parameters, clinical outcomes, or treatment recommendations.

The event-tree calculation has no antigen cap: depth n for initial load k reads
exactly the necessary finite triangular set of types 0, ..., k+n. Its bounds
come from unresolved event-tree leaves, not a truncated biological model.
"""

import argparse
from collections.abc import Mapping
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from math import comb


def rational(value):
    """Accept exact rational inputs; reject accidental binary floating point."""
    if isinstance(value, (float, bool)):
        raise TypeError("Use an integer, Fraction, Decimal, or rational string")
    return Fraction(value)


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError("%s must be an integer >= %s" % (name, minimum))
    return value


class Parameters:
    """Exact rates; ell=0 is supported with its separate initial-load exponent.

    b,bH,kappa are strictly positive. d,u,ell,h,dH are nonnegative.
    The main loss-and-gain theorem assumes u,ell>0; allowing zero rates here
    supports explicit checks of the degenerate boundaries, not that theorem's
    silent extension to those boundaries.
    """

    names = ("b", "d", "u", "ell", "kappa", "h", "bH", "dH")

    def __init__(self, b=1, d="1/2", u=1, ell="1/4", kappa=1,
                 h="1/10", bH=1, dH="1/4"):
        for name, value in zip(self.names, (b, d, u, ell, kappa, h, bH, dH)):
            value = rational(value)
            if value < 0 or (name in ("b", "bH", "kappa") and value == 0):
                raise ValueError("Invalid nonnegative/positive rate: " + name)
            setattr(self, name, value)

    @property
    def a(self):
        return self.ell + self.kappa

    @property
    def c(self):
        return self.ell / self.a

    @property
    def C(self):
        return self.u * self.kappa / self.a ** 2

    @property
    def r(self):
        return self.b - self.d - self.u * self.kappa / self.a

    @property
    def sigma(self):
        return max(Fraction(0), 1 - self.dH / self.bH)

    def source_exponent(self, k=0):
        """Mean/source threshold exponent before H conversion, for one load."""
        _integer(k, "k")
        return self.r if self.ell > 0 else self.r - self.kappa * k

    def replace(self, **updates):
        if set(updates) - set(self.names):
            raise ValueError("Unknown rate name")
        return Parameters(**{name: updates.get(name, getattr(self, name))
                             for name in self.names})

    def as_dict(self):
        return {name: str(getattr(self, name)) for name in self.names}


def exp_bounds(x, terms=64):
    """Certified rational enclosure of exp(x) by a positive Taylor tail."""
    x = rational(x)
    _integer(terms, "terms", 1)
    if x < 0:
        lower, upper = exp_bounds(-x, terms)
        return 1 / upper, 1 / lower
    term = Fraction(1)
    total = term
    for n in range(1, terms + 1):
        term *= x / n
        total += term
    following = term * x / (terms + 1)
    ratio = x / (terms + 2)
    if ratio >= 1:
        raise ValueError("Increase terms until the exponential tail ratio < 1")
    return total, total + following / (1 - ratio)


def mean_bounds(p, k, t, z=0, terms=64):
    """Enclose exp(-z*t) M_k(t), including source conversion when z=h."""
    _integer(k, "k")
    t, z = rational(t), rational(z)
    if t < 0:
        raise ValueError("t must be nonnegative")
    e_lower, e_upper = exp_bounds(-p.a * t, terms)
    a_lower = p.c + (1 - p.c) * e_lower
    a_upper = p.c + (1 - p.c) * e_upper
    power_lower = (p.r - z) * t + p.C * (1 - e_upper)
    power_upper = (p.r - z) * t + p.C * (1 - e_lower)
    lower = exp_bounds(power_lower, terms)[0] * a_lower ** k
    upper = exp_bounds(power_upper, terms)[1] * a_upper ** k
    return lower, upper


def resolvent_bounds(p, k, z=0, terms=32):
    """Enclose R_k(z)=integral_0^infinity exp(-z*t) M_k(t) dt.

    Each retained term is rational and nonnegative. The remaining terms are
    bounded by their first term and a decreasing-ratio geometric series.
    For ell>0 this requires z>r. When ell=0 it requires z>r-kappa*k;
    zero binomial weights are skipped before testing the denominator.
    """
    _integer(k, "k")
    _integer(terms, "terms", 1)
    z = rational(z)
    if z <= p.source_exponent(k):
        raise ValueError("The occupation resolvent is infinite at this z")
    beta = (z - p.r) / p.a
    lower, upper = Fraction(0), Fraction(0)
    for j in range(k + 1):
        weight = comb(k, j) * p.c ** (k - j) * (1 - p.c) ** j
        if weight == 0:
            continue
        v = beta + j
        if v <= 0:
            raise ValueError("Nonintegrable nonzero binomial component")
        term = 1 / v
        partial = term
        for n in range(1, terms):
            term *= p.C / (v + n)
            partial += term
        following = term * p.C / (v + terms)
        ratio = p.C / (v + terms + 1)
        if ratio >= 1:
            raise ValueError("Increase terms until the resolvent tail ratio < 1")
        lower += weight * partial
        upper += weight * (partial + following / (1 - ratio))
    return lower / p.a, upper / p.a


def _grid_integer(value, denominator, up):
    scaled = value.numerator * denominator
    if up:
        return -((-scaled) // value.denominator)
    return scaled // value.denominator


def event_tree_bounds(p, max_k, depth=64, bits=None):
    """Return q_k intervals, k=0..max_k, from an exact event-tree recursion.

    q is eventual total extinction, not merely extinction of antigen-presenting
    cells. Every iteration is rounded outward to a dyadic grid. The returned
    bounds are valid at every depth, including supercritical and h=0 cases.
    They converge when h*sigma>0 or the source is almost surely extinct.
    If h*sigma=0 and the source is supercritical, the upper iteration can remain
    one; a narrow interval is never asserted without inspecting its width.
    """
    _integer(max_k, "max_k")
    _integer(depth, "depth")
    if bits is None:
        bits = 2 * depth + 24
    _integer(bits, "bits", 1)
    grid = 1 << bits
    size = max_k + depth + 1
    lower, upper = [0] * size, [grid] * size
    for _ in range(depth):
        next_lower, next_upper = [], []
        for k in range(len(lower) - 1):
            death = p.d + p.kappa * k
            loss = p.ell * k
            rate = p.b + death + p.u + loss + p.h
            constant = death + p.h * (1 - p.sigma)
            for old, new, up in ((lower, next_lower, False),
                                 (upper, next_upper, True)):
                value = (p.b * Fraction(old[k] ** 2, grid ** 2)
                         + constant + p.u * Fraction(old[k + 1], grid)
                         + (loss * Fraction(old[k - 1], grid) if k else 0)) / rate
                new.append(_grid_integer(value, grid, up))
        lower, upper = next_lower, next_upper
    return tuple((Fraction(lo, grid), Fraction(hi, grid))
                 for lo, hi in zip(lower, upper))


def extinction_bounds(p, k, depth=64, bits=None):
    """The single-founder component of event_tree_bounds."""
    _integer(k, "k")
    return event_tree_bounds(p, k, depth, bits)[k]


def founder_product_bounds(intervals, counts):
    """Multiply single-founder intervals with exact integer multiplicities."""
    if not isinstance(counts, Mapping):
        raise TypeError("counts must map antigen load to founder multiplicity")
    lower, upper = Fraction(1), Fraction(1)
    for k, number in counts.items():
        _integer(k, "k")
        _integer(number, "founder multiplicity")
        lo, hi = intervals[k]
        if not (0 <= lo <= hi <= 1):
            raise ValueError("Invalid probability interval")
        lower *= lo ** number
        upper *= hi ** number
    return lower, upper


def population_extinction_bounds(p, counts, depth=64, bits=None):
    """Exact product brackets for a finite tumour with specified load counts."""
    if not isinstance(counts, Mapping):
        raise TypeError("counts must map antigen load to founder multiplicity")
    if not counts:
        return Fraction(1), Fraction(1)
    for k, number in counts.items():
        _integer(k, "k")
        _integer(number, "founder multiplicity")
    intervals = event_tree_bounds(p, max(counts), depth, bits)
    return founder_product_bounds(intervals, counts)


def solve_extinction(p, k, tolerance=Fraction(1, 10 ** 8),
                     initial_depth=16, max_depth=256):
    """Refine until an observed exact bracket meets tolerance, or fail loudly."""
    _integer(k, "k")
    tolerance = rational(tolerance)
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    _integer(initial_depth, "initial_depth", 1)
    _integer(max_depth, "max_depth", initial_depth)
    if p.h * p.sigma == 0:
        if p.source_exponent(k) - p.h <= 0:
            return Fraction(1), Fraction(1), 0
        raise ValueError("This bracketing solver cannot shrink a supercritical "
                         "source without successful H escape")
    depth = initial_depth
    while True:
        lower, upper = extinction_bounds(p, k, depth)
        if upper - lower <= tolerance:
            return lower, upper, depth
        if depth == max_depth:
            raise ArithmeticError("Maximum depth reached; certified width is "
                                  + str(upper - lower))
        depth = min(2 * depth, max_depth)


def total_source_event_bounds(p, k, terms=32):
    """Expected total source events, counting gain/loss as events.

    Requires r-h<0 (or the corresponding ell=0 exponent). The possibly
    negative coefficient of R is handled by interval arithmetic, not by
    silently multiplying the upper endpoint.
    """
    lower, upper = resolvent_bounds(p, k, p.h, terms)
    constant = 1 + p.ell / p.kappa
    coefficient = (2 * p.b + p.u
                   + p.ell * (p.b - p.d - p.h) / p.kappa)
    ends = (constant + coefficient * lower,
            constant + coefficient * upper)
    return max(Fraction(0), min(ends)), max(ends)


def event_depth_error_bound(p, k, depth, bits=None, terms=32):
    """A priori gap bound: W/(depth+1) plus outward rounding error."""
    _integer(depth, "depth")
    if bits is None:
        bits = 2 * depth + 24
    _integer(bits, "bits", 1)
    event_upper = total_source_event_bounds(p, k, terms)[1]
    rounding = Fraction(2 * ((1 << depth) - 1), 1 << bits)
    return min(Fraction(1), event_upper / (depth + 1) + rounding)


def rare_h_error_bound(p, k, terms=32):
    """Upper bound on sigma*h*T_k minus the successful-H probability."""
    t0 = resolvent_bounds(p, 0, 0, terms)[1]
    tk = resolvent_bounds(p, k, 0, terms)[1]
    return p.sigma * p.h ** 2 * (t0 + p.b * t0 ** 2) * tk


def rare_h_escape_bounds(p, k, terms=32):
    """Finite-h error enclosure of the rare-H expansion, assuming r<0."""
    t_lower, t_upper = resolvent_bounds(p, k, 0, terms)
    error = rare_h_error_bound(p, k, terms)
    return (max(Fraction(0), p.sigma * p.h * t_lower - error),
            min(Fraction(1), p.sigma * p.h * t_upper))


def finite_h_escape_bounds(p, k, terms=32):
    """First-birth lower bound and all-source-arrivals upper bound.

    sigma*h*R_k(h+2b) <= P(successful H escape) <= sigma*h*R_k(h).
    Source extinction requires r-h<=0; this computable integral bound uses
    the strict condition r-h<0. Under it total extinction is 1 minus H escape.
    """
    lower = resolvent_bounds(p, k, p.h + 2 * p.b, terms)[0]
    upper = resolvent_bounds(p, k, p.h, terms)[1]
    return p.sigma * p.h * lower, min(Fraction(1), p.sigma * p.h * upper)


def eigenfunction_identity(p, k):
    """Return two exact sides of A(c**k)=r*c**k for inspection/testing."""
    _integer(k, "k")
    value = p.c ** k
    gain = p.u * (p.c ** (k + 1) - value)
    loss = p.ell * k * (p.c ** (k - 1) - value) if k else Fraction(0)
    lhs = gain + loss + (p.b - p.d - p.kappa * k) * value
    return lhs, p.r * value


def display(value, digits=18):
    """Decimal text is presentation only and never used in certificates."""
    value = rational(value)
    with localcontext() as context:
        context.prec = digits
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def display_interval(bounds):
    return {"lower": display(bounds[0]), "upper": display(bounds[1]),
            "width": display(bounds[1] - bounds[0])}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", type=int, default=96)
    parser.add_argument("--series", type=int, default=32)
    args = parser.parse_args(argv)
    p = Parameters()
    loads = (0, 1, 4, 16, 64)
    intervals = event_tree_bounds(p, max(loads), args.depth)
    report = {
        "label": "Nonclinical computational mathematics; synthetic rates; "
                 "no treatment advice; no biological validation",
        "rates": p.as_dict(),
        "source_r": str(p.r),
        "source_r_after_H": str(p.r - p.h),
        "H_establishment_probability": str(p.sigma),
        "event_depth": args.depth,
        "series_terms": args.series,
        "single_founder": [],
    }
    for k in loads:
        q = intervals[k]
        report["single_founder"].append({
            "initial_antigens": k,
            "extinction": display_interval(q),
            "occupation_T": display_interval(resolvent_bounds(p, k, 0, args.series)),
            "escape_sandwich": display_interval(finite_h_escape_bounds(p, k, args.series)),
            "scaled_escape_k_p": display_interval((k * (1 - q[1]), k * (1 - q[0]))),
        })
    report["large_load_limit_k_p"] = str(p.sigma * p.h / p.kappa)
    report["three_founders_load_one"] = display_interval(
        founder_product_bounds(intervals, {1: 3}))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
