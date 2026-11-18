"""Exact certificates for the accompanying nonclinical mathematical note.

No treatment advice or biological calibration. All scientific calculations use
fractions.Fraction; Decimal is used only for outward-rounded display. The
stochastic simulator is neither imported nor used. The bounds below certify
the explicitly stated limiting model, not an error bound for its progression
time at a particular finite tumour size.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, fields, replace
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
from math import isqrt
from typing import Optional


F = Fraction
Interval = tuple[Fraction, Fraction]


def _rational(value) -> Fraction:
    if isinstance(value, float):
        raise TypeError("Use an integer, string or Fraction, not a binary float")
    return F(value)


@dataclass(frozen=True)
class RationalParameters:
    """Rates and densities; times have the common arbitrary model time unit."""

    b_s: Fraction = F(1)
    d_s: Fraction = F(7, 10)
    b_r: Fraction = F(4, 5)
    d_r: Fraction = F(3, 5)
    alpha: Fraction = F(3, 25)
    beta: Fraction = F(1, 50)
    a_sr: Fraction = F(1, 50)
    s0: Fraction = F(1)
    epsilon: Fraction = F(1, 5)
    L: Fraction = F(16)
    k_max: Fraction = F(33, 100)
    hill_emax: Fraction = F(33, 50)
    hill_ec50: Fraction = F(1)
    u_max: Fraction = F(1)
    hill_h: int = 2

    def __post_init__(self):
        for field in fields(self):
            if field.name != "hill_h":
                object.__setattr__(self, field.name,
                                   _rational(getattr(self, field.name)))
        if (not isinstance(self.hill_h, int)
                or isinstance(self.hill_h, bool) or self.hill_h < 1):
            raise ValueError("hill_h must be a positive integer")

    @property
    def r_s(self):
        return self.b_s - self.d_s

    @property
    def g0(self):
        return self.b_r - self.d_r

    @property
    def delta(self):
        return self.k_max - self.r_s

    @property
    def B(self):
        return self.s0 + self.epsilon

    @property
    def g_c(self):
        return self.g0 - self.alpha * self.s0

    @property
    def J(self):
        return self.L * (1 + self.delta / self.r_s)


DEFAULT = RationalParameters()


def hill_kill(u, p: RationalParameters = DEFAULT) -> Fraction:
    """The same bounded Hill kill law is used by both policies."""
    u = _rational(u)
    if u < 0 or p.hill_ec50 <= 0 or p.hill_emax <= 0:
        raise ValueError("nonnegative dose and positive Hill parameters required")
    return p.hill_emax * u**p.hill_h / (
        p.hill_ec50**p.hill_h + u**p.hill_h)


def validate(p: RationalParameters = DEFAULT) -> tuple[str, ...]:
    """Check the theorem's parameter assumptions, including strict boundaries.

    Mutation probability is lambda/N. Its nonnegativity and being at most one
    are checked by functions taking those arguments. Positive resistant
    growth requires either g_c > beta*epsilon or g_c <= 0; the intermediate
    region is deliberately outside the stated progression-time theorem.
    """
    checks = {
        "positive birth rates": p.b_s > 0 and p.b_r > 0,
        "nonnegative deaths and competition": min(
            p.d_s, p.d_r, p.alpha, p.beta, p.a_sr) >= 0,
        "positive initial density, margin and course": min(
            p.s0, p.epsilon, p.L) > 0,
        "positive sensitive growth and clearance": p.r_s > 0 and p.delta > 0,
        "positive bounded Hill law": min(
            p.hill_emax, p.hill_ec50, p.u_max) > 0,
        "maximum kill matches common dose law": hill_kill(p.u_max, p) == p.k_max,
        "containment kill stays positive": p.r_s > p.a_sr * p.B,
        "sensitive population remains macroscopic":
            (p.delta + p.a_sr * p.B) * p.L < 1,
    }
    if p.g0 > 0:
        checks.update({
            "resistant growth above progression burden": p.g0 > p.beta * p.B,
            "course exceeds resistant logarithmic growth time": p.L > 1 / p.g0,
            "containment endpoint regime":
                p.g_c <= 0 or p.g_c > p.beta * p.epsilon,
        })
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError("Unmet theorem assumptions: " + "; ".join(failed))
    return tuple(checks)


def _ceil(x: Fraction) -> int:
    return -((-x.numerator) // x.denominator)


def _round_interval(lo: Fraction, hi: Fraction, denominator: int) -> Interval:
    """Outward rounding to a rational mesh; never changes enclosure direction."""
    return (F((lo.numerator * denominator) // lo.denominator, denominator),
            F(_ceil(hi * denominator), denominator))


def _positive_tolerance(tolerance) -> Fraction:
    tolerance = _rational(tolerance)
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    return tolerance


def log_interval(x, tolerance=F(1, 10**24)) -> Interval:
    """Enclose log(x) using a positive atanh series and its geometric tail.

    For 1 <= y <= 2, z=(y-1)/(y+1) <= 1/3. After m terms the
    remainder in 2*sum(z**(2*j+1)/(2*j+1)) is at most
    2*z**(2*m+1)/((2*m+1)*(1-z*z)). Range reduction uses powers of two.
    """
    x, tolerance = _rational(x), _positive_tolerance(tolerance)
    if x <= 0:
        raise ValueError("log argument must be positive")
    if x < 1:
        lo, hi = log_interval(1 / x, tolerance)
        return -hi, -lo
    if x == 1:
        return F(0), F(0)
    exponent, reduced = 0, x
    while reduced >= 2:
        reduced /= 2
        exponent += 1
    target = tolerance / (4 * (exponent + 1))

    def series(y):
        z = (y - 1) / (y + 1)
        power, total, m = z, F(0), 0
        while True:
            tail = 2 * power / ((2 * m + 1) * (1 - z * z))
            if tail <= target:
                return total, total + tail
            total += 2 * power / (2 * m + 1)
            power *= z * z
            m += 1

    lo, hi = series(reduced)
    if exponent:
        lo2, hi2 = series(F(2))
        lo, hi = lo + exponent * lo2, hi + exponent * hi2
    return _round_interval(lo, hi, _ceil(8 / tolerance))


def exp_neg_interval(x, tolerance=F(1, 10**24)) -> Interval:
    """Enclose exp(-x), x>=0, by positive Taylor sums and inversion.

    For the partial sum through k, the exponential-series tail is at most
    next_term/(1-y/(k+2)). We reduce x to y<=1 and square interval endpoints;
    rational outward rounding after every squaring controls denominator size.
    The map z->z**2 is 2-Lipschitz on [0,1], controlling accumulated rounding.
    """
    x, tolerance = _rational(x), _positive_tolerance(tolerance)
    if x < 0:
        raise ValueError("exp_neg argument must be nonnegative")
    if x == 0:
        return F(1), F(1)
    squarings, y = 0, x
    while y > 1:
        y /= 2
        squarings += 1
    multiplier = 2**squarings
    denominator = _ceil(16 * multiplier / tolerance)
    target = tolerance / (4 * multiplier)
    total, term, k = F(1), F(1), 0
    while True:
        next_term = term * y / (k + 1)
        tail = next_term / (1 - y / (k + 2))
        lo, hi = 1 / (total + tail), 1 / total
        if hi - lo <= target:
            break
        total += next_term
        term = next_term
        k += 1
    lo, hi = _round_interval(lo, hi, denominator)
    for _ in range(squarings):
        lo, hi = _round_interval(lo * lo, hi * hi, denominator)
    return max(F(0), lo), min(F(1), hi)


def survival_interval(s, p: RationalParameters = DEFAULT,
                      tolerance=F(1, 10**22)) -> Interval:
    """Eventual survival of one R lineage in s(t)=s*exp(-delta*t).

    p_M(s)=(g0/b_r)/F(s), with
    F(s)=sum_{k>=0}(alpha*s)**k/prod_{j=1}^k(g0+j*delta).
    Following the term of index k, all subsequent term ratios are bounded
    by alpha*s/(g0+(k+2)*delta). Once this is below one it bounds the tail.
    """
    s, tolerance = _rational(s), _positive_tolerance(tolerance)
    if s < 0 or p.alpha < 0 or p.b_r <= 0 or p.d_r < 0 or p.delta <= 0:
        raise ValueError("survival formula needs nonnegative density and valid rates")
    if p.g0 <= 0:
        return F(0), F(0)
    p0, z = p.g0 / p.b_r, p.alpha * s
    if not z:
        return p0, p0
    total, term, k = F(1), F(1), 0
    while True:
        next_term = term * z / (p.g0 + (k + 1) * p.delta)
        ratio = z / (p.g0 + (k + 2) * p.delta)
        if ratio < 1:
            tail = next_term / (1 - ratio)
            lo, hi = p0 / (total + tail), p0 / total
            if hi - lo <= tolerance / 2:
                return _round_interval(lo, hi, _ceil(8 / tolerance))
        total += next_term
        term = next_term
        k += 1


def fourth_derivative_bound(p: RationalParameters = DEFAULT) -> Fraction:
    """Bound |p_M''''| by 75*(alpha/delta)^4*(g0/b_r).

    F(s)=E[exp((alpha/delta)*s*X)] for 0<=X<=1, so F^(j)/F<=c^j.
    Differentiating 1/F four times gives coefficients 24,-36,6,8,-1:
    (1/F)''''=(1/F)*(24*a1^4-36*a1^2*a2+6*a2^2+8*a1*a3-a4),
    where aj=F^(j)/F. Their absolute sum is 75 and p_M<=g0/b_r.
    """
    if p.delta <= 0 or p.b_r <= 0 or p.alpha < 0:
        raise ValueError("positive delta and b_r, nonnegative alpha required")
    return 75 * (p.alpha / p.delta)**4 * max(F(0), p.g0 / p.b_r)


def supply_interval(p: RationalParameters = DEFAULT, panels: int = 256,
                    tolerance=F(1, 10**22)) -> Interval:
    """Enclose A=(b_s/delta)*integral_0^s0 p_M(s) ds.

    Composite Simpson integration is accompanied by the rigorous bound
    s0*(s0/panels)^4*sup|p_M''''|/180; endpoint function values are themselves
    rational enclosures. Thus this is a certificate, not convergence evidence.
    """
    if not isinstance(panels, int) or panels < 2 or panels % 2:
        raise ValueError("panels must be a positive even integer >=2")
    tolerance = _positive_tolerance(tolerance)
    if p.s0 <= 0 or p.b_s <= 0 or p.delta <= 0 or p.alpha < 0:
        raise ValueError("positive s0, b_s and delta, nonnegative alpha required")
    if p.g0 <= 0:
        return F(0), F(0)
    if p.alpha == 0:
        value = p.b_s * p.s0 * p.g0 / (p.delta * p.b_r)
        return value, value
    h, low_sum, high_sum = p.s0 / panels, F(0), F(0)
    for i in range(panels + 1):
        lo, hi = survival_interval(i * h, p, tolerance)
        weight = 1 if i in (0, panels) else (4 if i % 2 else 2)
        low_sum += weight * lo
        high_sum += weight * hi
    error = p.s0 * h**4 * fourth_derivative_bound(p) / 180
    factor = p.b_s / p.delta
    return (factor * max(F(0), h * low_sum / 3 - error),
            factor * min(p.s0 * p.g0 / p.b_r, h * high_sum / 3 + error))


def phase_interval(p: RationalParameters = DEFAULT, panels: int = 256,
                   tolerance=F(1, 10**22)) -> Optional[Interval]:
    """Enclose lambda*=log((J-d)/(c-d))/A; None means no finite phase boundary.

    This boundary is for positive lambda and the theorem's normalized ordinary
    expected progression times. Alpha=0 has M strictly better at every finite
    positive lambda. If g0<=0, M also has the larger coefficient, J>L.
    """
    validate(p)
    if p.g0 <= 0 or p.alpha == 0:
        return None
    d = 1 / p.g0
    c = p.L if p.g_c <= 0 else min(p.L, 1 / p.g_c)
    a_lo, a_hi = supply_interval(p, panels, tolerance)
    if a_lo <= 0:
        raise ArithmeticError("Increase Simpson panels to obtain a positive A bound")
    log_lo, log_hi = log_interval((p.J - d) / (c - d), tolerance)
    return log_lo / a_hi, log_hi / a_lo


def asymptotic_mean_intervals(lam, p: RationalParameters = DEFAULT,
                             panels: int = 256,
                             tolerance=F(1, 10**22)) -> dict[str, Interval]:
    """Certified limits of E[T]/log N, under the README theorem assumptions.

    At lambda=0 or g0<=0 the coefficients are C=L, M=J. For lambda>0 and
    g0>0 they are C=min(L,1/g_c) (L when g_c<=0), and
    M=d+(J-d)*exp(-lambda*A), d=1/g0. The zero-lambda boundary is separate.
    """
    validate(p)
    lam = _rational(lam)
    if lam < 0:
        raise ValueError("mutation supply must be nonnegative")
    if lam == 0 or p.g0 <= 0:
        return {"containment": (p.L, p.L), "maximum": (p.J, p.J)}
    d = 1 / p.g0
    c = p.L if p.g_c <= 0 else min(p.L, 1 / p.g_c)
    a_lo, a_hi = supply_interval(p, panels, tolerance)
    exp_lo, _ = exp_neg_interval(lam * a_hi, tolerance)
    _, exp_hi = exp_neg_interval(lam * a_lo, tolerance)
    return {"containment": (c, c),
            "maximum": (d + (p.J - d) * exp_lo,
                        d + (p.J - d) * exp_hi)}


def coupling_bound(N: int, lam, p: RationalParameters = DEFAULT,
                   tolerance=F(1, 10**24)) -> dict[str, Fraction | int]:
    """Finite-N shared-clock mismatch bounds before H=L log N or R=m.

    m=ceil(N**(1/4)); the upper bounds concern the branching approximation
    up to this stopping time only, not the whole progression time. Values
    above one remain valid but uninformative union bounds. The finite-N
    containment kill must remain nonnegative, including its mutation term.
    """
    validate(p)
    if not isinstance(N, int) or isinstance(N, bool) or N < 2:
        raise ValueError("N must be an integer >=2")
    lam = _rational(lam)
    if lam < 0 or lam > N:
        raise ValueError("lambda/N must be a probability in [0,1]")
    if p.r_s - lam * p.b_s / N - p.a_sr * p.B <= 0:
        raise ValueError("N is too small for the stated containment controller")
    root = isqrt(isqrt(N))
    m = root if root**4 == N else root + 1
    _, log_hi = log_interval(N, tolerance)
    H_hi = p.L * log_hi
    containment = p.beta * m * m * H_hi / N
    maximum = containment + (p.s0 * (lam * p.b_s + p.alpha * m)
        * (lam * p.b_s + p.a_sr * m) / (N * p.delta**2))
    return {"m": m, "containment": containment, "maximum": maximum}


def decimal_interval(interval: Interval, precision: int = 18) -> tuple[str, str]:
    """Outward-rounded decimal display; Fraction endpoints remain certificates."""
    output = []
    for value, rounding in zip(interval, (ROUND_FLOOR, ROUND_CEILING)):
        with localcontext() as context:
            context.prec = precision
            context.rounding = rounding
            output.append(str(Decimal(value.numerator) / Decimal(value.denominator)))
    return tuple(output)


def self_check(p: RationalParameters = DEFAULT, panels: int = 256) -> dict:
    """Small executable certificate ledger; independent tests live in tests/."""
    assumptions = validate(p)
    a = supply_interval(p, panels)
    phase = phase_interval(p, panels)
    assert hill_kill(0, p) == 0
    assert hill_kill(p.u_max, p) == p.k_max
    assert log_interval(1) == (F(0), F(0))
    assert exp_neg_interval(0) == (F(1), F(1))
    assert fourth_derivative_bound(p) >= 0
    if phase is not None:
        low_supply, high_supply = phase[0] / 2, 2 * phase[1]
        low = asymptotic_mean_intervals(low_supply, p, panels)
        high = asymptotic_mean_intervals(high_supply, p, panels)
        assert low["maximum"][0] > low["containment"][1]
        assert high["containment"][0] > high["maximum"][1]
    no_competition = replace(p, alpha=F(0))
    assert phase_interval(no_competition, panels) is None
    zero = asymptotic_mean_intervals(0, p, panels)
    assert zero["maximum"][0] > zero["containment"][1]
    return {"assumptions": assumptions, "A": a, "phase": phase}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panels", type=int, default=256)
    parser.add_argument("--lambda", dest="lam", default="0.1")
    args = parser.parse_args()
    ledger = self_check(panels=args.panels)
    print("NONCLINICAL COMPUTATIONAL MATHEMATICS; NO TREATMENT ADVICE")
    print("All parameter assumptions checked:", len(ledger["assumptions"]))
    print("A enclosure:", decimal_interval(ledger["A"]))
    if ledger["phase"] is not None:
        print("lambda* enclosure:", decimal_interval(ledger["phase"]))
    print("E[T]/log N limiting coefficients at lambda =", args.lam)
    for policy, interval in asymptotic_mean_intervals(
            F(args.lam), panels=args.panels).items():
        print(policy, decimal_interval(interval))
    for N in (10**8, 10**9, 10**10, 10**11):
        bounds = coupling_bound(N, F(args.lam))
        print("N =", N, "m =", bounds["m"], "branching mismatch upper bounds:",
              "C", decimal_interval((F(0), bounds["containment"]))[1],
              "M", decimal_interval((F(0), bounds["maximum"]))[1])


if __name__ == "__main__":
    main()
