#!/usr/bin/env python3
"""Nonclinical computational mathematics: exact rational extinction enclosures.

No patient data, treatment recommendations, or biological calibration.  All
arithmetic used for certificates is fractions.Fraction; Decimal is display only.
See README.md for the model, proofs, and the strictly short-horizon theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction


Q = Fraction
Schedule = tuple[tuple[Fraction, Fraction], ...]


@dataclass(frozen=True)
class Model:
    """Rates of clonal births, natural deaths, kill, and mutant offspring.

    Offspring events retain the parent.  b is the CLONAL birth rate; total
    offspring rates are b+mu for S, b+nu for A/B, and b for D.
    """

    b: Fraction = Q(2, 5)
    d: Fraction = Q(1, 10)
    k: Fraction = Q(1)
    c: Fraction = Q(2)
    mu: Fraction = Q(3, 5)
    p: Fraction = Q(1, 4)
    w: Fraction = Q(1, 2)
    nu: Fraction = Q(1, 5)

    def __post_init__(self):
        for name in ("b", "d", "k", "c", "mu", "p", "w", "nu"):
            value = getattr(self, name)
            if not isinstance(value, (int, Fraction)) or isinstance(value, bool):
                raise TypeError("Use integer or Fraction model parameters, not floats")
            object.__setattr__(self, name, Q(value))
        if min(self.b, self.d, self.k, self.mu, self.nu) < 0:
            raise ValueError("Rates must be nonnegative")
        if self.c < 1 or not 0 <= self.p <= 1 or not 0 <= self.w <= 1:
            raise ValueError("Require c >= 1 and p,w in [0,1]")

    def births(self):
        a = self.mu * (1 - self.p) * self.w
        r = self.mu * (1 - self.p) * (1 - self.w)
        z = self.mu * self.p
        return ((self.b, a, r, z), (Q(0), self.b, Q(0), self.nu),
                (Q(0), Q(0), self.b, self.nu), (Q(0), Q(0), Q(0), self.b))

    def deaths(self, u):
        return (self.d + self.k, self.d + self.c * self.k * (1 - u),
                self.d + self.k * u, self.d)


@dataclass(frozen=True)
class Interval:
    lower: Fraction
    upper: Fraction

    def __post_init__(self):
        if self.lower > self.upper:
            raise ValueError("Reversed interval")

    @property
    def width(self):
        return self.upper - self.lower

    @property
    def midpoint(self):
        return (self.lower + self.upper) / 2

    def __sub__(self, other):
        return Interval(self.lower - other.upper, self.upper - other.lower)


def _rational(value, name):
    if not isinstance(value, (int, Fraction)) or isinstance(value, bool):
        raise TypeError(f"{name} must be an integer or Fraction")
    return Q(value)


def combination(duration):
    duration = _rational(duration, "duration")
    if duration < 0:
        raise ValueError("duration must be nonnegative")
    return ((duration, Q(1, 2)),)


def alternating(duration, cycles=1, first="A"):
    duration = _rational(duration, "duration")
    if duration < 0 or not isinstance(cycles, int) or isinstance(cycles, bool) or cycles < 1:
        raise ValueError("Require nonnegative duration and positive integer cycles")
    if first not in ("A", "B"):
        raise ValueError("first must be A or B")
    h = duration / (2 * cycles)
    u = Q(1) if first == "A" else Q(0)
    return tuple((h, value) for _ in range(cycles) for value in (u, 1 - u))


def validate_schedule(schedule):
    checked = []
    for duration, u in schedule:
        duration = _rational(duration, "duration")
        u = _rational(u, "dose fraction")
        if duration < 0 or not 0 <= u <= 1:
            raise ValueError("Invalid duration or dose fraction")
        checked.append((duration, u))
    return tuple(checked)


def exposure(schedule):
    schedule = validate_schedule(schedule)
    return (sum((t * u for t, u in schedule), Q(0)),
            sum((t * (1 - u) for t, u in schedule), Q(0)))


def _taylor_step(start, births, deaths, h, order):
    """Exact polynomial evaluation; the caller supplies its remainder bound."""
    total = tuple(map(sum, births))
    coefficients = [[value] for value in start]
    for degree in range(order):
        next_values = []
        for i in range(4):
            value = deaths[i] if degree == 0 else Q(0)
            value -= (deaths[i] + total[i]) * coefficients[i][degree]
            for j in range(4):
                if births[i][j]:
                    value += births[i][j] * sum(
                        (coefficients[i][v] * coefficients[j][degree - v]
                         for v in range(degree + 1)), Q(0))
            next_values.append(value / (degree + 1))
        for i, value in enumerate(next_values):
            coefficients[i].append(value)
    output = []
    for row in coefficients:
        value = row[-1]
        for coefficient in reversed(row[:-1]):
            value = coefficient + h * value
        output.append(value)
    return output


def extinction_bounds(model, schedule, *, order=18, bits=72):
    """Enclose each one-founder extinction probability by the course's end.

    Chronological segments are integrated in REVERSE order, as required by
    backward generating functions.  No state-space/population truncation occurs.
    Increasing order and bits at fixed step partition makes the bounds converge.
    """
    schedule = validate_schedule(schedule)
    if not isinstance(order, int) or isinstance(order, bool) or order < 1:
        raise ValueError("order must be a positive integer")
    if not isinstance(bits, int) or isinstance(bits, bool) or bits < 1:
        raise ValueError("bits must be a positive integer")
    births = model.births()
    birth_max = max(map(sum, births))
    value = [Q(0)] * 4
    error = Q(0)
    scale = 2 ** bits
    rounding = Q(1, scale)
    for duration, u in reversed(schedule):
        deaths = model.deaths(u)
        rate = max(birth_max, *deaths)
        if not duration or not rate:
            continue
        majorant = 2 * rate
        # h*majorant <= 1/4; the partition does not depend on precision.
        target = 4 * majorant * duration
        steps = max(1, (target.numerator + target.denominator - 1) // target.denominator)
        h = duration / steps
        x = majorant * h
        tail = x ** (order + 1) / (1 - x)
        amplification = 1 / (1 - birth_max * h)
        for _ in range(steps):
            approximate = _taylor_step(value, births, deaths, h, order)
            value = [min(Q(1), max(Q(0), Q((v * scale).__floor__(), scale)))
                     for v in approximate]
            error = amplification * error + tail + rounding
    return tuple(Interval(max(Q(0), v - error), min(Q(1), v + error)) for v in value)


def population_bounds(bounds, initial=(1, 0, 0, 0)):
    if len(initial) != 4 or any(not isinstance(n, int) or isinstance(n, bool) or n < 0
                                for n in initial):
        raise ValueError("initial must contain four nonnegative integers")
    lower = upper = Q(1)
    for interval, count in zip(bounds, initial):
        lower *= interval.lower ** count
        upper *= interval.upper ** count
    return Interval(lower, upper)


def comparison_bounds(model, candidate, reference, *, initial=(1, 0, 0, 0), order=18, bits=72):
    """Positive lower bound certifies strict improvement; overlap is unresolved."""
    candidate = validate_schedule(candidate)
    reference = validate_schedule(reference)
    if exposure(candidate) != exposure(reference):
        raise ValueError("Comparisons must have equal per-drug exposure")
    left = population_bounds(extinction_bounds(model, candidate, order=order, bits=bits), initial)
    right = population_bounds(extinction_bounds(model, reference, order=order, bits=bits), initial)
    return left - right


def local_coefficients(model, cycles=1, first="A"):
    """Exact T^3,T^4 coefficients for one initially sensitive cell.

    These formulas are independently tested using the population-state Markov
    generator, not by differentiating this function or its ODE integrator.
    """
    alternating(Q(1), cycles, first)  # validate the family arguments
    a = model.mu * (1 - model.p) * model.w
    r = model.mu * (1 - model.p) * (1 - model.w)
    h = a * model.c - r
    ell = model.d * h + model.k * (a * model.c ** 2 - r) / 2
    mass = a * model.c ** 2 + r
    ds = model.d + model.k
    orientation = 1 if first == "A" else -1
    cubic = orientation * ds * model.k * h / (16 * cycles)
    odd = (h * (ds + 2 * model.b + model.mu + model.nu) * (16 * cycles ** 2 + 1)
           / (384 * cycles ** 3)
           + ell * (8 * cycles ** 2 - 1) / (192 * cycles ** 3))
    even = model.k * mass / (192 * cycles ** 2)
    quartic = -ds * model.k * (orientation * odd + even)
    return cubic, quartic


def local_witness_horizon(model, cycles=1):
    """Conservative proven horizon for the A-first local sign, not an optimum.

    Returns (T0, sign), with the sign guaranteed for 0<T<T0.  None means a
    degenerate case or a boundary outside the stated theorem's assumptions.
    """
    c3, c4 = local_coefficients(model, cycles)
    rate = max(model.d + model.c * model.k,
               model.b + model.mu, model.b + model.nu)
    if not rate:
        return None
    coefficient, degree = (c3, 3) if c3 else (c4, 4)
    if not coefficient:
        return None
    horizon = min(1 / (2 * rate), abs(coefficient) / (4 * rate ** (degree + 1)))
    return horizon, 1 if coefficient > 0 else -1


def display(value, digits=16):
    """Presentation only: decimal rounding does not enter any certificate."""
    with localcontext() as context:
        context.prec = digits
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=18)
    parser.add_argument("--bits", type=int, default=72)
    args = parser.parse_args()
    print("NONCLINICAL computational mathematics; synthetic rates; no treatment advice.")
    cases = (
        ("short improvement", Model(), Q(1, 2)),
        ("longer-course reversal", Model(), Q(2)),
        ("balanced boundary", Model(c=Q(1)), Q(1, 2)),
        ("Monte Carlo stress test", Model(b=Q(1, 20), d=Q(2), c=Q(8),
                                         mu=Q(6), p=Q(1, 10), nu=Q(1, 50)), Q(3, 10)),
    )
    for name, model, duration in cases:
        reference = extinction_bounds(model, combination(duration), order=args.order, bits=args.bits)[0]
        print(f"\n{name}: T={duration}; combo={display(reference.midpoint)} +/- {display(reference.width / 2, 4)}")
        for cycles in (1, 2):
            result = extinction_bounds(model, alternating(duration, cycles), order=args.order, bits=args.bits)[0]
            gap = result - reference
            verdict = "positive" if gap.lower > 0 else "negative" if gap.upper < 0 else "unresolved"
            print(f"  AB_{cycles}: q={display(result.midpoint)}; "
                  f"gap=[{display(gap.lower)}, {display(gap.upper)}] ({verdict})")


if __name__ == "__main__":
    main()
