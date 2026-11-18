#!/usr/bin/env python3
"""Nonclinical computational mathematics: exact polynomial proof certificate.

No treatment advice.  This module verifies identities in a polynomial ring,
not identities fitted to sampled numerical parameter values.  It uses only
Python integers and fractions.Fraction.  The independent population-state
generator checks live in the test file; this certificate expands PGF flows.

Variables are the four PGF arguments and b,d,k,a,z,nu,delta.  Here delta=c-1,
the two single-mutant offspring rates both equal a, and mu=2*a+z.  Substituting
a=mu*(1-p)/2, z=mu*p recovers the model with symmetric single-mutant supply.
All coefficients below are for half-course length h=T/2.
"""

from fractions import Fraction
from math import factorial
from time import perf_counter


Q = Fraction
_NAMES = ("xS", "xA", "xB", "xD", "b", "d", "k", "a", "z", "nu", "delta")
_INDEX = {name: index for index, name in enumerate(_NAMES)}
_ZERO = (0,) * len(_NAMES)


def _variable(name):
    powers = list(_ZERO)
    powers[_INDEX[name]] = 1
    return {tuple(powers): 1}


def _add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for powers, coefficient in polynomial.items():
            result[powers] = result.get(powers, 0) + coefficient
    return {powers: coefficient for powers, coefficient in result.items() if coefficient}


def _scale(polynomial, scalar):
    return {powers: scalar * coefficient for powers, coefficient in polynomial.items()
            if scalar * coefficient}


def _multiply(left, right):
    result = {}
    for lp, lc in left.items():
        for rp, rc in right.items():
            powers = tuple(a + b for a, b in zip(lp, rp))
            result[powers] = result.get(powers, 0) + lc * rc
    return {powers: coefficient for powers, coefficient in result.items() if coefficient}


def _product(*polynomials):
    result = {_ZERO: 1}
    for polynomial in polynomials:
        result = _multiply(result, polynomial)
    return result


def _derivative(polynomial, name):
    index = _INDEX[name]
    result = {}
    for powers, coefficient in polynomial.items():
        if powers[index]:
            reduced = list(powers)
            reduced[index] -= 1
            result[tuple(reduced)] = coefficient * powers[index]
    return result


def _at_pgf_zero(polynomial):
    return {powers: coefficient for powers, coefficient in polynomial.items()
            if not any(powers[:4])}


def _delta_coefficient(polynomial, degree):
    result = {}
    for powers, coefficient in polynomial.items():
        if powers[-1] == degree:
            reduced = powers[:-1] + (0,)
            result[reduced] = coefficient
    return result


def _lie(polynomial, field, max_pgf_degree, max_delta_degree):
    """Apply F dot grad, discarding terms unable to reach the required constant.

    Each further Lie operation can lower PGF degree by at most one.  Thus with
    r operations left, degrees greater than r cannot contribute at x=0.
    Delta is never differentiated, so its degree can only stay fixed or rise.
    Both truncations preserve every coefficient requested by this certificate.
    """
    result = {}
    for powers, coefficient in polynomial.items():
        for index in range(4):
            if not powers[index]:
                continue
            reduced = list(powers)
            reduced[index] -= 1
            multiplier = coefficient * powers[index]
            for fp, fc in field[index].items():
                degree = tuple(a + b for a, b in zip(reduced, fp))
                if sum(degree[:4]) > max_pgf_degree or degree[-1] > max_delta_degree:
                    continue
                result[degree] = result.get(degree, 0) + multiplier * fc
    return {powers: coefficient for powers, coefficient in result.items() if coefficient}


def _fields():
    v = {name: _variable(name) for name in _NAMES}
    x = [v[name] for name in _NAMES[:4]]
    b, d, k, a, z, nu, delta = (v[name] for name in _NAMES[4:])
    ds = _add(d, k)
    collateral_death = _add(d, k, _multiply(k, delta))
    births = ((b, a, a, z), ({}, b, {}, nu),
              ({}, {}, b, nu), ({}, {}, {}, b))

    def field(deaths):
        output = []
        for index in range(4):
            total = _add(*births[index])
            quadratic = _add(*(_product(x[index], births[index][j], x[j])
                               for j in range(4)))
            output.append(_add(deaths[index],
                               _scale(_multiply(_add(deaths[index], total), x[index]), -1),
                               quadratic))
        return output

    fa = field((ds, d, ds, d))
    fb = field((ds, collateral_death, d, d))
    return v, fa, fb


def _gap_coefficient(order, max_delta_degree, fa, fb):
    """Exact coefficient of h**order in AB(h,h) minus C(2h).

    Pullbacks act in reversed operator order: the S coordinate of
    Phi_A(h,Phi_B(h,0)) is [exp(h L_B) exp(h L_A) xS](0).
    Combination over 2h is exp(h*(L_A+L_B)).
    """
    x = _variable("xS")
    a_powers = [x]
    for j in range(order):
        a_powers.append(_lie(a_powers[-1], fa, order - j - 1, max_delta_degree))
    sequential = {}
    for j in range(order + 1):
        value = a_powers[j]
        for r in range(order - j):
            value = _lie(value, fb, order - j - r - 1, max_delta_degree)
        sequential = _add(sequential, _scale(_at_pgf_zero(value),
                                             Q(1, factorial(j) * factorial(order - j))))
    combined_field = [_add(left, right) for left, right in zip(fa, fb)]
    value = x
    for r in range(order):
        value = _lie(value, combined_field, order - r - 1, max_delta_degree)
    return _add(sequential, _scale(_at_pgf_zero(value), -Q(1, factorial(order))))


def _require_identity(actual, expected, label):
    difference = _add(actual, _scale(expected, -1))
    if difference:
        raise AssertionError(f"Polynomial identity failed: {label}; "
                             f"{len(difference)} nonzero monomials")


def verify():
    """Verify generic identities and return small metadata; raises on failure.

    This is a finite symbolic proof certificate over the rational polynomial
    ring.  No numerical substitutions, floating point, or external CAS occur.
    The return value deliberately contains no machine-dependent timing.
    """
    v, fa, fb = _fields()
    b, d, k, a, z, nu, delta = (v[name] for name in _NAMES[4:])
    ds = _add(d, k)
    mu = _add(_scale(a, 2), z)
    factor = _product(a, ds, k)
    ell = _add(_scale(b, 34), _scale(d, 31), _scale(k, 35),
               _scale(mu, 17), _scale(nu, 17))
    f4 = _add(_scale(k, -Q(1, 6)),
              _scale(_multiply(ell, delta), -Q(1, 24)),
              _scale(_product(k, delta, delta), -Q(3, 8)))
    f50 = _scale(_multiply(k, _add(_scale(b, 95), _scale(d, 51),
                                          _scale(k, 40), _scale(mu, 29),
                                          _scale(nu, 8))), Q(1, 120))

    # P5 and P6 are the complete polynomials, with mu standing for 2*a+z.
    p5 = _add(_scale(_multiply(mu, mu), 65),
              _scale(_multiply(b, mu), 235), _scale(_multiply(d, mu), 330),
              _scale(_multiply(nu, mu), 105),
              _multiply(k, _add(_scale(a, 613), _scale(z, 254))),
              _scale(_multiply(b, b), 225), _scale(_multiply(b, d), 760),
              _scale(_multiply(b, k), 855), _scale(_multiply(b, nu), 215),
              _scale(_multiply(d, d), 205), _scale(_multiply(d, k), 461),
              _scale(_multiply(d, nu), 280), _scale(_multiply(k, k), 243),
              _scale(_multiply(k, nu), 213), _scale(_multiply(nu, nu), 55))
    p6 = _add(_scale(_multiply(mu, mu), 134),
              _scale(_multiply(b, mu), 779), _scale(_multiply(d, mu), 668),
              _scale(_multiply(nu, mu), 83),
              _multiply(k, _add(_scale(a, 936), _scale(z, 361))),
              _scale(_multiply(b, b), 1291), _scale(_multiply(b, d), 2323),
              _scale(_multiply(b, k), 1739), _scale(_multiply(b, nu), 737),
              _scale(_multiply(d, d), 398), _scale(_multiply(d, k), 625),
              _scale(_multiply(d, nu), 187), _scale(_multiply(k, k), 248),
              _scale(_multiply(k, nu), 126), _scale(_multiply(nu, nu), 8))
    f51 = _scale(p5, Q(1, 120))
    f60 = _scale(_multiply(k, p6), -Q(1, 720))

    # Each Lie application raises delta degree by at most one.  Caps equal to
    # the time degree retain ALL collateral powers in these first identities.
    gap3 = _gap_coefficient(3, 3, fa, fb)
    gap4 = _gap_coefficient(4, 4, fa, fb)
    gap5 = _gap_coefficient(5, 1, fa, fb)
    gap6 = _gap_coefficient(6, 0, fa, fb)
    _require_identity(gap3, _scale(_multiply(factor, delta), Q(1, 2)), "cubic gap")
    _require_identity(gap4, _multiply(factor, f4), "complete quartic gap")
    _require_identity(_delta_coefficient(gap5, 0), _multiply(factor, f50),
                      "fifth gap at c=1")
    _require_identity(_delta_coefficient(gap5, 1), _multiply(factor, f51),
                      "collateral derivative of fifth gap at c=1")
    _require_identity(gap6, _multiply(factor, f60), "sixth gap at c=1")

    def redistribution(poly):
        # At fixed mu, da/dp=-mu/2 and dz/dp=mu.  Thus this is d/dp.
        return _multiply(mu, _add(_derivative(poly, "z"),
                                 _scale(_derivative(poly, "a"), -Q(1, 2))))

    _require_identity(redistribution(mu), {}, "mutation rate remains fixed")
    _require_identity(redistribution(f4), {}, "quartic has no p dependence")
    _require_identity(redistribution(f50), {}, "balanced fifth has no p dependence")
    d_f51 = redistribution(f51)
    d_f60 = redistribution(f60)
    _require_identity(d_f51, _scale(_multiply(mu, k), -Q(7, 16)),
                      "p derivative of f5 prime")
    _require_identity(d_f60, _scale(_product(mu, k, k), Q(107, 720)),
                      "p derivative of f6")
    # If c=1+alpha*h+beta*h**2+gamma*h**3+..., alpha=k/3.
    # Differentiate the root equation: gamma_p=-2*(alpha*f51_p+f60_p).
    # Dividing by 8 converts the h**3 coefficient to the T**3 coefficient.
    threshold_slope = _scale(_add(_multiply(_scale(k, Q(1, 3)), d_f51), d_f60),
                             -Q(1, 4))
    _require_identity(threshold_slope, _scale(_product(mu, k, k), -Q(1, 1440)),
                      "third threshold coefficient p slope")
    return {"identities": 11, "maximum_time_degree": 6,
            "arithmetic": "exact rational polynomial identities",
            "threshold_p_slope": "-mu*k**2/1440"}


def main():
    start = perf_counter()
    result = verify()
    print("NONCLINICAL computational mathematics; no treatment advice.")
    for name, value in result.items():
        print(f"{name}: {value}")
    print(f"elapsed_seconds: {perf_counter() - start:.3f}")


if __name__ == "__main__":
    main()
