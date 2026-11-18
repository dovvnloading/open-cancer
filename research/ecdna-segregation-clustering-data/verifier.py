"""Nonclinical exact certificates for copy-selective branching thresholds.

A float eigensolver proposes rational weights. Integer residual comparisons,
including every omitted tail row, decide whether a certificate is accepted.
No eigensolver tolerance is accepted as a proof. Standard library only.
"""
from fractions import Fraction
from functools import lru_cache
from math import comb, gcd, lcm, floor, ceil

F = Fraction


def rational(value):
    """Interpret decimal inputs as the decimal written by the caller."""
    if isinstance(value, bool):
        raise TypeError("boolean is not a model parameter")
    return value if isinstance(value, Fraction) else F(str(value))


def effective_r(p, rho):
    p, rho = rational(p), rational(rho)
    if not 0 <= p <= 1 or not 0 <= rho <= 1:
        raise ValueError("parameters must lie in [0,1]")
    return (2*p-1)**2 + 4*p*(1-p)*rho


def _integer_products(a, b, concentration, n):
    """Return common-denominator exact beta-binomial (or bounding) row."""
    scale = lcm(a.denominator, b.denominator, concentration.denominator)
    aa, bb, ss = (int(x*scale) for x in (a, b, concentration))
    pa, pb, denominator = [1], [1], 1
    for i in range(n):
        pa.append(pa[-1]*(aa+i*scale))
        pb.append(pb[-1]*(bb+i*scale))
        denominator *= ss+i*scale
    numerators = tuple(comb(n, j)*pa[j]*pb[n-j] for j in range(n+1))
    return numerators, denominator


@lru_cache(maxsize=256)
def integer_kernel(k, p, rho):
    """K(j|k), j=0,...,2k, as integers over one denominator."""
    if not isinstance(k, int) or k < 1:
        raise ValueError("k must be a positive integer")
    p, rho = rational(p), rational(rho)
    effective_r(p, rho)  # validation
    n = 2*k
    if rho == 1:
        return ((p.denominator-p.numerator,) + (0,)*(n-1)
                + (p.numerator,)), p.denominator
    if rho == 0 or p in (0, 1):
        a, d = p.numerator, p.denominator
        nums = tuple(comb(n, j)*a**j*(d-a)**(n-j) for j in range(n+1))
        return nums, d**n
    s = 1/rho-1
    nums, den = _integer_products(p*s, (1-p)*s, s, n)
    if sum(nums) != den:
        raise AssertionError("beta-binomial row is not normalized")
    return nums, den


def _matrix(rows, c, cap):
    c = float(c)
    matrix, exits = [], []
    for k, (num, den) in enumerate(rows, 1):
        scale = 1/(1+c*k)
        matrix.append([(num[j]+num[2*k-j])/den*scale if j <= 2*k else 0.0
                       for j in range(1, cap+1)])
        exits.append(sum(num[j]+num[2*k-j] for j in range(cap+1, 2*k+1))/den*scale)
    return matrix, exits


def _solve(matrix, rhs):
    """Float proposal only; the caller must independently check residuals."""
    n = len(rhs)
    a = [list(row)+[y] for row, y in zip(matrix, rhs)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda j: abs(a[j][k]))
        a[k], a[pivot] = a[pivot], a[k]
        if abs(a[k][k]) < 1e-17:
            raise ArithmeticError("singular proposal")
        for i in range(k+1, n):
            mult = a[i][k]/a[k][k]
            if mult:
                for j in range(k+1, n+1):
                    a[i][j] -= mult*a[k][j]
    x = [0.0]*n
    for i in range(n-1, -1, -1):
        x[i] = (a[i][-1]-sum(a[i][j]*x[j] for j in range(i+1, n)))/a[i][i]
    return x


def _finite_root(rows, cap):
    """PF root of I+diag(1/k)(A-I); returns c_N and a weight proposal."""
    a, _ = _matrix(rows, 0, cap)
    b = [[(x-float(i == j))/(i+1)+float(i == j) for j, x in enumerate(row)]
         for i, row in enumerate(a)]
    w = [1.0]*cap
    old = None
    for _ in range(400):
        v = [sum(x*y for x, y in zip(row, w)) for row in b]
        top = max(v)
        if top <= 0:
            return -1.0, w
        w = [x/top for x in v]
        if old is not None and max(abs(x-y) for x, y in zip(w, old)) < 2e-14:
            break
        old = w[:]
    ratios = [sum(x*y for x, y in zip(row, w))/wi for row, wi in zip(b, w)]
    return (min(ratios)+max(ratios))/2-1, w


def _certificate(rows, c, weights, tail=None):
    """Exact integer residual check. A true result is the certificate.

    If tail is None, weights vanish outside the head and every head row
    must have strictly positive M w - w. Otherwise every head row must
    have nonnegative w - M w and all tail rows obey the analytic bound.
    """
    c = rational(c)
    cap = len(rows)
    if len(weights) != cap or any(not isinstance(w, int) or w <= 0 for w in weights):
        return False
    cn, cd = c.numerator, c.denominator
    if c <= 0:
        return False
    for k, ((nums, den), wk) in enumerate(zip(rows, weights), 1):
        lhs = sum((nums[j]+nums[2*k-j])*weights[j-1]
                  for j in range(1, min(2*k, cap)+1))
        if tail is not None:
            lhs += tail*sum(nums[j]+nums[2*k-j] for j in range(cap+1, 2*k+1))
        lhs *= cd
        rhs = (cd+cn*k)*den*wk
        if (tail is None and lhs <= rhs) or (tail is not None and lhs > rhs):
            return False
    if tail is not None:
        if not isinstance(tail, int) or tail <= 0:
            return False
        if 2*cd*max(tail, max(weights)) >= (cd+cn*(cap+1))*tail:
            return False
    return True


def _lower_witness(rows, c, cap, proposal=None):
    if proposal is None:
        _, proposal = _finite_root(rows, cap)
    weights = tuple(max(1, round(w*10**13)) for w in proposal)
    return weights if _certificate(rows, c, weights) else None


def _upper_witness(rows, c, cap):
    matrix, exits = _matrix(rows, c, cap)
    a = [[float(i == j)-x for j, x in enumerate(row)] for i, row in enumerate(matrix)]
    try:
        x = _solve(a, [t+1e-8 for t in exits])
    except ArithmeticError:
        return None
    if min(x) <= 0 or max(x) > 1e5:
        return None
    scale = 10**13
    weights = tuple(max(1, round(w*scale)) for w in x)
    return weights if _certificate(rows, c, weights, scale) else None


@lru_cache(maxsize=128)
def threshold(p, rho, width="0.0001"):
    """Rigorous full infinite-type threshold enclosure for a fitted kernel.

    Includes an exact numerical enclosure only, not sampling uncertainty.
    Decimal arguments are treated as exact; use rectangle_threshold for
    source-rounding or parameter intervals. Runtime is normally <1 second.
    """
    p, rho, width = rational(p), rational(rho), rational(width)
    effective_r(p, rho)
    if width <= 0:
        raise ValueError("positive width required")
    if rho == 1 or p in (0, 1):
        return {"lower": 0.0, "upper": 0.0, "lower_exact": "0", "upper_exact": "0",
                "width": 0.0, "cap": 0, "certificate": "complete-clumping theorem"}
    step = width/4
    for cap in (20, 28, 40, 56, 80):
        rows = tuple(integer_kernel(k, p, rho) for k in range(1, cap+1))
        root, proposal = _finite_root(rows, cap)
        low = max(F(0), floor(root/float(step))*step-step)
        high = (ceil(root/float(step))+1)*step
        lw = _lower_witness(rows, low, cap, proposal) if low > 0 else None
        uw = _upper_witness(rows, high, cap)
        if lw is not None and uw is not None:
            return {"lower": float(low), "upper": float(high),
                    "lower_exact": str(low), "upper_exact": str(high),
                    "width": float(high-low), "cap": cap,
                    "certificate": "exact integer residuals and all-type tail"}
    # Failure to certify is not evidence for a parameter value.
    raise ArithmeticError("requested threshold enclosure not certified; widen width")


def _factor(bias, rho, successes, trials):
    # Successes observed so far in the chosen ordering; trials total so far.
    return (bias*(1-rho)+successes*rho)/(1-rho+trials*rho)


@lru_cache(maxsize=512)
def interval_kernel(k, pl, pu, rl, ru, upper):
    """Outward bound on each PMF entry over an entire parameter rectangle.

    Each Polya-urn factor is monotone in bias and, for fixed bias, monotone
    in rho (derivative sign successes-bias*trials). Its extrema are corners.
    Multiplying positive lower/upper factors gives a valid simultaneous
    entrywise envelope; its entries need not add to one.
    """
    pl, pu, rl, ru = map(rational, (pl, pu, rl, ru))
    if not 0 < pl <= pu < 1 or not 0 <= rl <= ru < 1:
        raise ValueError("interior bias and rho<1 required for finite interval kernel")
    n = 2*k
    extremum = max if upper else min
    values = []
    for j in range(n+1):
        value = F(comb(n, j))
        for t in range(j):
            value *= extremum(_factor(p, r, t, t) for p in (pl, pu) for r in (rl, ru))
        for t in range(n-j):
            value *= extremum(_factor(1-p, r, t, j+t) for p in (pl, pu) for r in (rl, ru))
        values.append(value)
    denominator = lcm(*(v.denominator for v in values))
    return tuple(int(v*denominator) for v in values), denominator


def rectangle_threshold(pl, pu, rl, ru, subdivisions=1, cap=20):
    """Uniform rigorous threshold bounds; no parameter monotonicity assumed.

    Uniformly bound all kernels in a rectangular simultaneous parameter set.
    subdivisions splits each axis into equal cells. Broad rectangles may
    return only the general [0,(1-Rmin)/2] bound; this is honest uncertainty.
    Returned extrema are certificates for all points, not corner guesses.
    """
    pl, pu, rl, ru = map(rational, (pl, pu, rl, ru))
    if not 0 <= pl <= pu <= 1 or not 0 <= rl <= ru <= 1:
        raise ValueError("invalid parameter rectangle")
    if not isinstance(subdivisions, int) or subdivisions < 1:
        raise ValueError("subdivisions must be a positive integer")
    pclosest = min(max(F(1, 2), pl), pu)
    rmin = effective_r(pclosest, rl)
    universal_upper = (1-rmin)/2
    if pl == 0 or pu == 1 or ru == 1:
        return {"lower": 0.0, "upper": float(universal_upper),
                "lower_exact": "0", "upper_exact": str(universal_upper),
                "method": "universal R bound", "cells": 1}
    lower_bounds, upper_bounds = [], []
    step = F(1, 1000)
    for ip in range(subdivisions):
        for ir in range(subdivisions):
            pa, pb = pl+(pu-pl)*ip/subdivisions, pl+(pu-pl)*(ip+1)/subdivisions
            ra, rb = rl+(ru-rl)*ir/subdivisions, rl+(ru-rl)*(ir+1)/subdivisions
            lower_rows = tuple(interval_kernel(k, pa, pb, ra, rb, False) for k in range(1, cap+1))
            upper_rows = tuple(interval_kernel(k, pa, pb, ra, rb, True) for k in range(1, cap+1))
            low_root, w = _finite_root(lower_rows, cap)
            low = max(F(0), floor(low_root/float(step))*step-step)
            if low > 0 and _lower_witness(lower_rows, low, cap, w) is None:
                low = F(0)
            upper_root, _ = _finite_root(upper_rows, cap)
            high = (ceil(upper_root/float(step))+1)*step
            if high <= 0 or high >= universal_upper:
                high = universal_upper
            else:
                while high < universal_upper and _upper_witness(upper_rows, high, cap) is None:
                    high += step
                high = min(high, universal_upper)
            lower_bounds.append(low)
            upper_bounds.append(high)
    lower, upper = min(lower_bounds), max(upper_bounds)
    return {"lower": float(lower), "upper": float(upper),
            "lower_exact": str(lower), "upper_exact": str(upper),
            "method": "uniform Polya-factor rectangles plus all-type tail",
            "cells": subdivisions**2}


def _small_box_rows(pl, pu, rl, ru, cap, upper):
    """Exact positive-factor Pochhammer bounds, efficient for tiny boxes."""
    if rl == 0:
        return tuple(interval_kernel(k, pl, pu, rl, ru, upper) for k in range(1, cap+1))
    sl, su = 1/ru-1, 1/rl-1
    if upper:
        a, b, s = pu*su, (1-pl)*su, sl
    else:
        a, b, s = pl*sl, (1-pu)*sl, su
    return tuple(_integer_products(a, b, s, 2*k) for k in range(1, cap+1))


@lru_cache(maxsize=64)
def rounding_threshold(p, rho, p_unit="0.00000001", rho_unit="0.00000001",
                       numerical_step="0.000001"):
    """Uniform enclosure for nearest-unit source rounding, not a sampling CI.

    p and rho each vary by half their stated published decimal unit. Every
    kernel in this rectangle is included; no threshold monotonicity is used.
    numerical_step controls grid rounding of certified threshold endpoints.
    """
    p, rho, p_unit, rho_unit, step = map(rational, (p, rho, p_unit, rho_unit, numerical_step))
    if p_unit <= 0 or rho_unit <= 0 or step <= 0:
        raise ValueError("positive decimal units and numerical step required")
    pl, pu = p-p_unit/2, p+p_unit/2
    rl, ru = max(F(0), rho-rho_unit/2), rho+rho_unit/2
    if not 0 < pl <= pu < 1 or not 0 <= rl <= ru < 1:
        raise ValueError("rounding rectangle must have interior bias and rho<1")
    for cap in (20, 28, 40, 56, 80):
        lower_rows = _small_box_rows(pl, pu, rl, ru, cap, False)
        upper_rows = _small_box_rows(pl, pu, rl, ru, cap, True)
        lower_root, proposal = _finite_root(lower_rows, cap)
        upper_root, _ = _finite_root(upper_rows, cap)
        lo = max(F(0), floor(lower_root/float(step))*step-step)
        hi = (ceil(upper_root/float(step))+1)*step
        if (lo > 0 and _lower_witness(lower_rows, lo, cap, proposal) is not None
                and _upper_witness(upper_rows, hi, cap) is not None):
            return {"lower": float(lo), "upper": float(hi),
                    "lower_exact": str(lo), "upper_exact": str(hi),
                    "width": float(hi-lo), "cap": cap,
                    "p_lower_exact": str(pl), "p_upper_exact": str(pu),
                    "rho_lower_exact": str(rl), "rho_upper_exact": str(ru),
                    "certificate": "uniform source-rounding rectangle and all-type tail"}
    raise ArithmeticError("source-rounding rectangle not certified; increase numerical step")


def moment_threshold_bounds(rl, ru, precision="0.000000001"):
    """Uniform threshold enclosure from an interval for effective R alone.

    Assumptions: a fixed, copy-conserving exchangeable allocation mixture P,
    complementary sisters, division rate b, and linear death gamma*k. This
    is an identification envelope, not the beta-binomial threshold curve.

    Put X=2P-1 after random sister orientation, R=E X^2, T=E X^4.
    Jensen and |X|<=1 give R^2<=T<=R. For the two-type principal matrix,
    its determinant sign is that of
        (c+R)(1+8c+6R-3T) - (1+R)(1-T).
    This is affine in T. Substituting T=R and T=R^2 gives two quadratics
    increasing in R for c>=0. Below BOTH positive roots at ru, the finite
    matrix is supercritical for every R in [rl,ru]. PR60 gives the upper
    bound (1-rl)/2. Exact bisection bounds the lower endpoint from below.
    """
    rl, ru, precision = map(rational, (rl, ru, precision))
    if not 0 <= rl <= ru <= 1 or precision <= 0:
        raise ValueError("require 0<=rl<=ru<=1 and positive precision")
    r = ru
    coefficients = ((1+11*r, r+4*r*r-1),
                    (1+14*r-3*r*r, 7*r*r-2*r*r*r-1))
    roots = []
    for linear, constant in coefficients:
        if constant >= 0:
            roots.append(F(0))
            continue
        lo, hi = F(0), F(1, 2)
        while hi-lo > precision:
            mid = (lo+hi)/2
            if 8*mid*mid+linear*mid+constant < 0:
                lo = mid
            else:
                hi = mid
        roots.append(lo)
    lower, upper = min(roots), (1-rl)/2
    return {"lower": float(lower), "upper": float(upper),
            "lower_exact": str(lower), "upper_exact": str(upper),
            "lower_numerical_error_bound": float(precision),
            "method": "uniform R-only two-type lower and all-type upper"}


if __name__ == "__main__":
    import json
    print(json.dumps({"binomial": threshold("0.5", "0"),
                      "synthetic_example": threshold("0.45", "0.1")}, indent=2))
