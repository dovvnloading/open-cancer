"""Exact certificates for nonclinical ecDNA branching mathematics.

No patient data or treatment advice. This module uses only the Python standard
library. Fractions, including rational tail inequalities, are the certificates;
decimal displays are only displays. The README supplies the infinite-type
proofs. A finite matrix alone never certifies eradication of the full model.

``kernel(k)`` is the law of the ordered left daughter J. The right daughter is
always 2*k-J. ``p(k)`` is the division-before-death probability.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from math import comb
import argparse
import json


F = Fraction


def rational(value, name="value"):
    """Accept exact rational inputs; reject implicit conversion of floats."""
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError(name + " must be an int or Fraction")
    return F(value)


def nonnegative_integer(value, name="value"):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(name + " must be an integer")
    if value < 0:
        raise ValueError(name + " must be nonnegative")
    return value


def probability(value, name="probability"):
    value = rational(value, name)
    if not 0 <= value <= 1:
        raise ValueError(name + " must lie in [0,1]")
    return value


@lru_cache(maxsize=None, typed=True)
def binomial_kernel(k, bias=F(1, 2)):
    """Exact Bin(2*k,bias) law; bias labels one daughter, not both sisters."""
    nonnegative_integer(k, "k")
    bias = probability(bias, "bias")
    return tuple(F(comb(2*k, j)) * bias**j * (1-bias)**(2*k-j)
                 for j in range(2*k+1))


@lru_cache(maxsize=None, typed=True)
def beta_binomial_kernel(k, rho, bias=F(1, 2)):
    """Beta-binomial law, with intraclass correlation rho in [0,1].

    For 0<rho<1 and 0<bias<1, draw P~Beta(a,b), then J~Bin(2*k,P),
    where a=bias*(1/rho-1), b=(1-bias)*(1/rho-1). Rho=0 is binomial;
    rho=1 is the weak limit placing all copies in one daughter. Bias=0 or 1
    is deterministic. Every returned entry is an exact Fraction.
    """
    nonnegative_integer(k, "k")
    rho = probability(rho, "rho")
    bias = probability(bias, "bias")
    if rho == 0 or k == 0 or bias in (0, 1):
        return binomial_kernel(k, bias)
    n = 2*k
    if rho == 1:
        return (1-bias,) + (F(0),)*(n-1) + (bias,)
    concentration = 1/rho-1
    a, b = bias*concentration, (1-bias)*concentration
    first = F(1)
    for index in range(n):
        first *= (b+index)/(a+b+index)
    values = [first]
    for j in range(n):
        values.append(values[-1] * F(n-j, j+1) * (a+j)/(b+n-j-1))
    return tuple(values)


def checked_kernel(k, kernel):
    values = tuple(probability(x, "kernel entry") for x in kernel(k))
    if len(values) != 2*k+1 or sum(values) != 1:
        raise ValueError("kernel must be a probability law on 0,...,2*k")
    return values


def division_probability(p, k):
    return probability(p(k), "division probability")


def kernel_moments(k, kernel=binomial_kernel):
    """Ordered-daughter and uniformly pooled-sister moments at fixed mother k."""
    nonnegative_integer(k, "k")
    law = checked_kernel(k, kernel)
    mean = sum(F(j)*value for j, value in enumerate(law))
    second = sum(F(j*j)*value for j, value in enumerate(law))
    pooled_variance = sum(F((j-k)**2)*value for j, value in enumerate(law))
    return {
        "ordered_mean": mean,
        "ordered_variance": second-mean*mean,
        "sister_covariance": -(second-mean*mean),
        "pooled_mean": F(k),
        "pooled_variance": pooled_variance,
        "sister_difference_second": 4*pooled_variance,
        "pooled_zero_fraction": (law[0]+law[-1])/2,
    }


def mean_matrix(cap, p, kernel=binomial_kernel):
    """Positive-type mean matrix for the killed-at-exit finite subprocess.

    Rows are mothers 1,...,cap; columns are positive daughters 1,...,cap.
    Types above cap are discarded only in this mathematical subprocess.
    """
    nonnegative_integer(cap, "cap")
    if cap == 0:
        raise ValueError("cap must be positive")
    result = []
    for k in range(1, cap+1):
        pk = division_probability(p, k)
        law = checked_kernel(k, kernel)
        result.append(tuple(pk*(law[j]+law[2*k-j]) if j <= 2*k else F(0)
                            for j in range(1, cap+1)))
    return tuple(result)


def mean_matrix_and_exit(cap, p, kernel=binomial_kernel):
    """Return the finite matrix and mean number of daughters of type >cap."""
    matrix = mean_matrix(cap, p, kernel)
    exits = []
    for k, row in enumerate(matrix, 1):
        law = checked_kernel(k, kernel)
        total_positive = division_probability(p, k)*(2-law[0]-law[-1])
        exits.append(total_positive-sum(row))
    return matrix, tuple(exits)


def matvec(matrix, vector):
    if any(len(row) != len(vector) for row in matrix):
        raise ValueError("matrix and vector shapes disagree")
    return tuple(sum((x*y for x, y in zip(row, vector)), F(0)) for row in matrix)


def solve_fraction(matrix, rhs):
    """Exact Gaussian elimination with row pivoting, without numerical fitting."""
    n = len(rhs)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("a square matrix matching rhs is required")
    rows = [[rational(x) for x in row]+[rational(y)]
            for row, y in zip(matrix, rhs)]
    for k in range(n):
        if rows[k][k] == 0:
            pivot = next((j for j in range(k+1, n) if rows[j][k]), None)
            if pivot is None:
                raise ValueError("singular matrix")
            rows[k], rows[pivot] = rows[pivot], rows[k]
        pivot_value = rows[k][k]
        for j in range(k+1, n):
            if rows[j][k] == 0:
                continue
            multiplier = rows[j][k]/pivot_value
            for ell in range(k+1, n+1):
                rows[j][ell] -= multiplier*rows[k][ell]
            rows[j][k] = F(0)
    result = [F(0)]*n
    for k in range(n-1, -1, -1):
        result[k] = (rows[k][-1]-sum(rows[k][j]*result[j]
                                      for j in range(k+1, n)))/rows[k][k]
    return tuple(result)


def is_irreducible(matrix):
    """Directed reachability check for a finite nonnegative matrix."""
    n = len(matrix)
    if not n or any(len(row) != n for row in matrix):
        return False
    for transpose in (False, True):
        visited, stack = {0}, [0]
        while stack:
            i = stack.pop()
            for j in range(n):
                value = matrix[j][i] if transpose else matrix[i][j]
                if value > 0 and j not in visited:
                    visited.add(j)
                    stack.append(j)
        if len(visited) != n:
            return False
    return True


def certificate(p, cap, tail_division_bound, kernel=binomial_kernel):
    """Construct an exact finite witness with an explicit infinite-tail bound.

    Caller must prove p(k)<=tail_division_bound for EVERY k>cap. This function
    checks the finite arithmetic but cannot verify an arbitrary callable over
    infinitely many types. The README proves the supplied linear/power bounds.

    Write t for mean exits above cap and solve (I-M_cap)x=t.
    If x>0, set h=x on the head, h=1 on the tail. The checked inequality
    2*tail_division_bound*max(1,max(x))<=1 proves Mh<=h. A checked positive
    uniform death floor then makes this an eradication certificate.
    If x<0 and the head matrix is irreducible, h=-x satisfies M_cap*h>=h,
    strictly somewhere, proving r(M_cap)>1 and finite-subprocess persistence.
    All other outcomes are explicitly undetermined.
    """
    tail_division_bound = probability(tail_division_bound, "tail bound")
    matrix, exits = mean_matrix_and_exit(cap, p, kernel)
    identity_minus = tuple(tuple(F(i == j)-value for j, value in enumerate(row))
                           for i, row in enumerate(matrix))
    try:
        x = solve_fraction(identity_minus, exits)
    except ValueError as error:
        if str(error) != "singular matrix":
            raise
        return {"regime": "undetermined", "reason": "singular finite system",
                "cap": cap}
    # This residual is part of the certificate, not a floating point tolerance.
    if matvec(identity_minus, x) != exits:
        raise AssertionError("exact linear-system residual failed")
    result = {"regime": "undetermined", "cap": cap, "solution": x,
              "exit": exits, "tail_division_bound": tail_division_bound}
    if all(value < 0 for value in x) and any(exits) and is_irreducible(matrix):
        weights = tuple(-value for value in x)
        residual = tuple(a-b for a, b in zip(matvec(matrix, weights), weights))
        if residual != exits or not all(value >= 0 for value in residual):
            raise AssertionError("invalid finite persistence witness")
        result.update(regime="persistence", weights=weights, residual=residual)
    elif all(value > 0 for value in x):
        tail_lhs = 2*tail_division_bound*max(F(1), max(x))
        death_floor = 1-max(tail_division_bound,
                            max(division_probability(p, k) for k in range(1, cap+1)))
        result.update(weights=x, tail_lhs=tail_lhs, death_floor=death_floor)
        if tail_lhs <= 1 and death_floor > 0:
            result["regime"] = "eradication"
    return result


def linear_certificate(c, cap, kernel=binomial_kernel):
    """Certified full-model result for p(k)=1/(1+c*k), c>0."""
    c = rational(c, "c")
    if c <= 0:
        raise ValueError("c must be positive")
    return certificate(lambda k: 1/(1+c*k), cap, 1/(1+c*(cap+1)), kernel)


def certify_linear_threshold():
    """A width-1e-12 rational enclosure, with all-type eradication proof."""
    lower, upper = F("0.378474994554"), F("0.378474994555")
    low = linear_certificate(lower, 25)
    high = linear_certificate(upper, 48)
    if (low["regime"] != "persistence" or high["regime"] != "eradication"
            or high["tail_lhs"] >= F(1, 4)):
        raise AssertionError("linear critical bracket certificate failed")
    return {"lower": lower, "upper": upper, "width": upper-lower,
            "lower_cap": 25, "upper_cap": 48,
            "upper_tail_bound": F(1, 4)}


def certify_segregation_thresholds():
    """Rational full-model intervals and a matched-second-moment obstruction.

    The intervals are deliberately modest in precision: every endpoint is
    proved using a finite witness plus, on the upper side, an infinite-tail
    inequality. They demonstrate distinct thresholds without treating a
    floating point eigensolver or an imposed type cap as proof.
    """
    specifications = (
        ("beta", F(1, 41), F("0.367"), F("0.368"), 8, 10),
        ("beta", F(1, 4), F("0.257"), F("0.258"), 14, 20),
        ("beta", F(1, 3), F("0.216"), F("0.217"), 16, 24),
        ("beta", F(1, 2), F("0.137"), F("0.138"), 32, 36),
        ("bias", F(3, 4), F("0.270"), F("0.271"), 20, 20),
        ("bias", F(9, 10), F("0.120"), F("0.122"), 28, 36),
    )
    cases = []
    for kind, parameter, lower, upper, lower_cap, upper_cap in specifications:
        if kind == "beta":
            kernel = lambda k, value=parameter: beta_binomial_kernel(k, value)
        else:
            kernel = lambda k, value=parameter: binomial_kernel(k, value)
        low = linear_certificate(lower, lower_cap, kernel)
        high = linear_certificate(upper, upper_cap, kernel)
        if (low["regime"] != "persistence" or high["regime"] != "eradication"
                or high["tail_lhs"] >= 1):
            raise AssertionError("segregation threshold certificate failed")
        cases.append({"kind": kind, "parameter": parameter,
                      "lower": lower, "upper": upper,
                      "lower_cap": lower_cap, "upper_cap": upper_cap})
    c = F(13, 50)
    beta = linear_certificate(c, 14, lambda k: beta_binomial_kernel(k, F(1, 4)))
    biased = linear_certificate(c, 8, lambda k: binomial_kernel(k, F(3, 4)))
    if (beta["regime"] != "eradication" or biased["regime"] != "persistence"
            or beta["tail_lhs"] >= 1):
        raise AssertionError("matched-second-moment obstruction failed")
    return {"cases": cases, "same_second_moment": {
        "c": c, "rho": F(1, 4), "bias": F(3, 4),
        "beta_regime": beta["regime"], "bias_regime": biased["regime"],
        "beta_cap": 14, "bias_cap": 8,
    }}


def validate_initial(initial):
    return tuple(nonnegative_integer(k, "initial copy count") for k in initial)


def forward_expectations(initial, generations, p, kernel=binomial_kernel):
    """Exact full reachable-support expected counts at a finite generation."""
    initial = validate_initial(initial)
    nonnegative_integer(generations, "generations")
    counts = {k: F(number) for k, number in Counter(initial).items()}
    for _ in range(generations):
        updated = defaultdict(Fraction)
        for k, count in counts.items():
            pk = division_probability(p, k)
            for j, chance in enumerate(checked_kernel(k, kernel)):
                contribution = count*pk*chance
                if contribution:
                    updated[j] += contribution
                    updated[2*k-j] += contribution
        counts = dict(updated)
    return counts


def population_extinction(initial, generations, p, kernel=binomial_kernel):
    """Exact probability all positive descendants are absent by generation n."""
    initial = validate_initial(initial)
    nonnegative_integer(generations, "generations")

    @lru_cache(maxsize=None)
    def extinction(k, n):
        if k == 0:
            return F(1)
        if n == 0:
            return F(0)
        pk = division_probability(p, k)
        return 1-pk+pk*sum(chance*extinction(j, n-1)*extinction(2*k-j, n-1)
                           for j, chance in enumerate(checked_kernel(k, kernel)))

    answer = F(1)
    for k in initial:
        answer *= extinction(k, generations)
    return answer


def population_raw_moments(initial, generations, p, kernel=binomial_kernel,
                           statistic="positive"):
    """Exact first and raw second moments, retaining correlated sisters."""
    initial = validate_initial(initial)
    nonnegative_integer(generations, "generations")
    statistics = {"cells": lambda k: 1, "positive": lambda k: int(k > 0),
                  "zero": lambda k: int(k == 0), "copies": lambda k: k,
                  "copy_squares": lambda k: k*k}
    if isinstance(statistic, str):
        if statistic not in statistics:
            raise ValueError("unknown statistic")
        function = statistics[statistic]
    elif callable(statistic):
        function = statistic
    else:
        raise TypeError("statistic must be a known name or callable")

    @lru_cache(maxsize=None)
    def moments(k, n):
        if n == 0:
            value = rational(function(k), "statistic")
            return value, value*value
        mean, second = F(0), F(0)
        pk = division_probability(p, k)
        for j, chance in enumerate(checked_kernel(k, kernel)):
            left_mean, left_second = moments(j, n-1)
            right_mean, right_second = moments(2*k-j, n-1)
            mean += pk*chance*(left_mean+right_mean)
            second += pk*chance*(left_second+right_second+2*left_mean*right_mean)
        return mean, second

    components = [moments(k, generations) for k in initial]
    mean = sum((item[0] for item in components), F(0))
    variance = sum((second-first*first for first, second in components), F(0))
    return {"mean": mean, "variance": variance, "second": variance+mean*mean}


def summary_from_counts(counts):
    cells = sum(counts.values(), F(0))
    zero = counts.get(0, F(0))
    copies = sum((k*count for k, count in counts.items()), F(0))
    squares = sum((k*k*count for k, count in counts.items()), F(0))
    mean = copies/cells if cells else None
    variance = squares/cells-mean*mean if cells else None
    return {"cells": cells, "positive": cells-zero, "zero": zero,
            "copies": copies, "copy_squares": squares,
            "pooled_mean": mean, "pooled_variance": variance,
            "pooled_fano": variance/mean if mean else None,
            "pooled_zero_fraction": zero/cells if cells else None}


def two_generation_prediction(c, segregation_variance):
    """Theorem 6: one-copy founders, linear killing, pooled generation two.

    segregation_variance is R=rho for symmetric beta clustering, or
    R=(2*bias-1)**2 for fixed asymmetric allocation. The equality of these
    second-moment curves does not identify the eradication threshold.
    """
    c = rational(c, "c")
    r = probability(segregation_variance, "segregation variance")
    if c < 0:
        raise ValueError("c must be nonnegative")
    a = 2+(3-r)*c
    d = 2+4*c+(1+r)*c*c
    b = 3+r+4*c
    derivative_term = 2+4*c+(3-r)*c*c
    return {
        "mean": a/d,
        "fano": (1+r)*b/a+(1-r)/2-a/d,
        "mean_derivative": -(1+r)*derivative_term/(d*d),
        "fano_derivative": (1+r)*(derivative_term/(d*d)-(1-r*r)/(a*a)),
        "beta_minus_bias_zero_fraction": r*r*(1-r*r)*(1+c)/(16*(1+2*r)*d),
        "beta_fourth_statistic": 3*r*r/(1+2*r),
        "bias_fourth_statistic": r*r,
    }


def refuge_bound(band_start):
    """Exact lower bound on the killed-binomial sine-band row ratios.

    Types L+1,...,2*L-1 and division floor p_min persist if
    2*p_min*refuge_bound(L)>1. The all-L proof is in Lemma 3A.
    """
    nonnegative_integer(band_start, "band_start")
    if band_start < 2:
        raise ValueError("band_start must be at least two")
    return 1-F(5, band_start)-F(192, band_start**3)


def two_type_power_witness(c, two_to_alpha, three_to_alpha):
    """Exact cap-two witness for power killing when these powers are exact.

    Caller supplies 2**alpha and 3**alpha for the SAME positive alpha.
    This finite routine does not certify that relation for arbitrary inputs.
    Tests use positive integer exponents. The analytic limit covers real alpha.
    """
    c = rational(c, "c")
    two = rational(two_to_alpha, "two_to_alpha")
    three = rational(three_to_alpha, "three_to_alpha")
    if c <= 0 or two <= 1 or three <= 1:
        raise ValueError("require c>0 and powers from a positive exponent")
    determinant_numerator = two*c*c+c/4-F(1, 4)
    if determinant_numerator == 0:
        return {"regime": "undetermined", "determinant": F(0)}
    x1 = F(5, 16)/determinant_numerator
    x2 = 2*c*x1
    tail_lhs = 2*max(F(1), x1, x2)/(1+c*three)
    regime = "persistence" if determinant_numerator < 0 else (
        "eradication" if tail_lhs <= 1 else "undetermined")
    return {"regime": regime, "determinant": determinant_numerator,
            "solution": (x1, x2), "tail_lhs": tail_lhs}


def verify():
    """Deterministic exact-arithmetic checks; independent tests add enumeration."""
    checks = 0

    def require(condition):
        nonlocal checks
        if not condition:
            raise AssertionError("exact check " + str(checks+1) + " failed")
        checks += 1

    for k in (0, 1, 2, 5, 13, 31):
        for bias in (F(0), F(1, 5), F(1, 2), F(4, 5), F(1)):
            for rho in (F(0), F(1, 7), F(1, 3), F(3, 4), F(1)):
                kernel = lambda value, r=rho, a=bias: beta_binomial_kernel(value, r, a)
                law = kernel(k)
                require(sum(law) == 1)
                moments = kernel_moments(k, kernel)
                require(moments["ordered_mean"] == 2*k*bias)
                ordered_variance = 2*k*bias*(1-bias)*(1+(2*k-1)*rho)
                require(moments["ordered_variance"] == ordered_variance)
                require(moments["pooled_variance"] == ordered_variance+k*k*(2*bias-1)**2)
                require(moments["sister_covariance"] == -ordered_variance)
                require(sum((j+(2*k-j))*chance for j, chance in enumerate(law)) == 2*k)
    bracket = certify_linear_threshold()
    require(bracket["width"] == F(1, 10**12))
    require(bracket["upper_tail_bound"] < 1)
    segregation = certify_segregation_thresholds()
    for case in segregation["cases"]:
        require(case["lower"] < case["upper"])
    for rho in (F(0), F(1, 3), F(3, 4), F(1)):
        kernel = lambda k, r=rho: beta_binomial_kernel(k, r)
        for c in (F(0), F(1, 4), F(1, 2), F(2)):
            p = lambda k, value=c: 1/(1+value*k)
            for initial in ((), (0,), (1,), (0, 1, 2)):
                for generations in range(3):
                    counts = forward_expectations(initial, generations, p, kernel)
                    summary = summary_from_counts(counts)
                    require(population_raw_moments(initial, generations, p, kernel)["mean"]
                            == summary["positive"])
                    q = population_extinction(initial, generations, p, kernel)
                    require(0 <= q <= 1)
                    if c == 0 and any(initial):
                        require(q == 0)
                    if not any(initial):
                        require(q == 1)
    for r in (F(0), F(1, 4), F(1, 3), F(1, 2), F(1)):
        kernel = lambda k, value=r: beta_binomial_kernel(k, value)
        for c in (F(0), F(1, 10), F(13, 50), F(2, 5), F(3)):
            theory = two_generation_prediction(c, r)
            summary = summary_from_counts(forward_expectations(
                (1,), 2, lambda k, value=c: 1/(1+value*k), kernel))
            require(theory["mean"] == summary["pooled_mean"])
            require(theory["fano"] == summary["pooled_fano"])
            require(theory["mean_derivative"] < 0)
            require(theory["fano_derivative"] > 0)
        a1 = (1-r)/2
        for k in (1, 2, 7, 19, 47):
            law = kernel(k)
            both_positive = 1-law[0]-law[-1]
            require(both_positive <= min(k*a1, 1))
            if r < 1:
                require((1+both_positive)/(1+a1*k) <= 1)
    # The decreasing law c/k can remain lethal at low copies yet has a
    # rigorously supercritical finite high-copy band.
    for c, band in ((F(10), 64), (F(100), 512)):
        p_min = 1/(1+c/(band+1))
        require(2*p_min*refuge_bound(band) > 1)
    # Exact all-integer-power witnesses of the proved alpha->infinity limit.
    for alpha in (8, 12, 20, 40):
        for constant, regime in ((F(2, 5), "persistence"),
                                 (F(3, 5), "eradication")):
            c = constant/F(2**(alpha//2))
            witness = two_type_power_witness(c, 2**alpha, 3**alpha)
            require(witness["regime"] == regime)
            general = certificate(lambda k, value=c, power=alpha:
                                  1/(1+value*k**power), 2,
                                  1/(1+c*3**alpha))
            require(general["regime"] == regime)
            require(general["solution"] == witness["solution"])
    for h in (F(3), F(5), F(20)):
        require(8-(2*h)-h*h < 0)
        require((h+1)/(h+1+(h+1)) == F(1, 2))
    return {"exact_checks": checks, "linear_threshold": bracket,
            "segregation_thresholds": segregation,
            "refuge_example": {"c": 10, "alpha": -1, "band_start": 64,
                                "mean_lower_bound": 2*F(65, 75)*refuge_bound(64)}}


def json_safe(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(item) for item in value]
    return value


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--linear-only", action="store_true",
                        help="run just the rational threshold enclosure")
    args = parser.parse_args()
    result = certify_linear_threshold() if args.linear_only else verify()
    print(json.dumps(json_safe(result), indent=2))
