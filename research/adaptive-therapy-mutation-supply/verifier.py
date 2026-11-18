"""Exact checks for the nonclinical two-cell mutation/containment note.

All model calculations use fractions.Fraction.  Decimal is used only to print
outward-rounded root bounds.  No patient data or treatment recommendations.
See README.md for the assumptions, complete proof, and limits of the model.
"""

from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
import json


F = Fraction
ZERO = F(0)
ONE = F(1)
START = (2, 0)
CURE = (0, 0)
PROGRESSION = None
STATES = {
    "C": ((0, 0), (0, 1), (0, 2), (1, 1), (2, 0)),
    "M": ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0)),
}
# Coefficients throughout are in ascending order.
MAIN_POLYNOMIAL = (F(400), F(-517), F(-185), F(18))


@dataclass(frozen=True)
class Rates:
    """a,b: R birth rates at total 1,2; d: R death; k: S kill; z: censor."""

    a: Fraction = F(100)
    b: Fraction = F(1)
    d: Fraction = F(3)
    k: Fraction = F(2)
    z: Fraction = F(1)

    def __post_init__(self):
        for name in ("a", "b", "d", "k", "z"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
                raise TypeError("rates must be exact integers or Fractions")
            if value <= 0:
                raise ValueError("all five rates must be strictly positive")
            object.__setattr__(self, name, F(value))
        if self.a < self.b:
            raise ValueError("competition requires a >= b")


DEFAULT = Rates()


def mutation_rate(mu):
    if isinstance(mu, bool) or not isinstance(mu, (int, Fraction)):
        raise TypeError("mu must be an exact integer or Fraction")
    mu = F(mu)
    if mu < 0:
        raise ValueError("mu must be nonnegative")
    return mu


def poly_trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p) or (ZERO,)


def poly_add(p, q):
    out = [ZERO] * max(len(p), len(q))
    for i, value in enumerate(p):
        out[i] += value
    for i, value in enumerate(q):
        out[i] += value
    return poly_trim(out)


def poly_multiply(p, q):
    out = [ZERO] * (len(p) + len(q) - 1)
    for i, left in enumerate(p):
        for j, right in enumerate(q):
            out[i + j] += left * right
    return poly_trim(out)


def poly_value(p, x):
    value = ZERO
    for coefficient in reversed(p):
        value = value * x + coefficient
    return value


def sign_variations(p):
    signs = [1 if c > 0 else -1 for c in reversed(p) if c]
    return sum(left != right for left, right in zip(signs, signs[1:]))


def reduced_constants(rates=DEFAULT):
    """Constants in Theorem 1; h is the pure-R two-cell Laplace transform."""
    a, b, d, k, z = rates.a, rates.b, rates.d, rates.k, rates.z
    q = a / (z + a + d)
    h = 2 * b / (z + 2 * b + 2 * d - 2 * d * q)
    aa, bb, cc = z + k, z + k + b + d, z + 2 * k
    dd, ee = 3 * z + 2 * b, z * (z + b + d)
    uu = aa + b / h + (2 * k + d) * q
    vv = aa * b / h + k * q * (aa + bb)
    return q, h, aa, bb, cc, dd, ee, uu, vv


def comparison_polynomial(rates=DEFAULT):
    """Return the exact cubic Psi; its sign equals h_C - h_M for mu > 0."""
    _, _, aa, bb, cc, dd, ee, uu, vv = reduced_constants(rates)
    left = poly_multiply((rates.b, ONE), (aa, ONE))
    left = poly_multiply(left, (bb, ONE))
    left = poly_multiply(left, (cc, F(2)))
    right = poly_multiply((vv, uu, ONE), (ee, dd, F(2)))
    return poly_add(left, tuple(-v for v in right))


def formula_laplace(mu, policy, rates=DEFAULT):
    """Closed-form E[exp(-z*T)]; exp(-z*infinity) means zero."""
    mu = mutation_rate(mu)
    _, h, aa, bb, cc, dd, ee, uu, vv = reduced_constants(rates)
    if policy == "C":
        return 2 * mu * h * (mu + rates.b) / (2 * mu**2 + dd * mu + ee)
    if policy == "M":
        return (2 * mu * h * (mu**2 + uu * mu + vv)
                / ((mu + aa) * (mu + bb) * (2 * mu + cc)))
    raise ValueError("policy must be C or M")


def formula_control_time(mu, policy, rates=DEFAULT):
    """E[min(T,H)] for independent H ~ Exp(z), not the uncensored mean."""
    return (1 - formula_laplace(mu, policy, rates)) / rates.z


def transitions(policy, state, mu, rates=DEFAULT):
    """Build the stopped CTMC directly from occupancy events, not formulas.

    Keys are destination (S,R) tuples or None for total-burden progression.
    Crucially, after S=0 a resistant death NEVER regenerates sensitive cells.
    """
    mu = mutation_rate(mu)
    if policy not in STATES or state not in STATES[policy]:
        raise ValueError("invalid policy or state")
    s, r = state
    out = {}

    def add(destination, rate):
        if rate:
            out[destination] = out.get(destination, ZERO) + rate

    if s:
        add((s - 1, r + 1), mu * s)
        if policy == "M":
            add((s - 1, r), rates.k * s)
    if r:
        birth = r * (rates.a if s + r == 1 else rates.b)
        if policy == "C" and s:
            add((s - 1, r + 1), birth)
            add((s + 1, r - 1), rates.d * r)
        else:
            add(PROGRESSION if s + r == 2 else (s, r + 1), birth)
            add((s, r - 1), rates.d * r)
    return out


def solve_linear(matrix, rhs):
    """Exact Gaussian elimination with pivoting, independent of reduction."""
    n = len(rhs)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("square matrix required")
    rows = [[F(x) for x in row] + [F(y)] for row, y in zip(matrix, rhs)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j]), None)
        if pivot is None:
            raise ValueError("singular system")
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = rows[j][j]
        rows[j] = [x / scale for x in rows[j]]
        for i in range(n):
            if i != j and rows[i][j]:
                scale = rows[i][j]
                rows[i] = [x - scale * y for x, y in zip(rows[i], rows[j])]
    return tuple(row[-1] for row in rows)


def markov_laplace(mu, policy, rates=DEFAULT, *, eventual=False):
    """Solve the occupancy generator independently of the closed forms.

    With eventual=True solve the eventual progression probability (z=0).
    Cure is a boundary with value 0, not an invertible transient state.
    """
    mu = mutation_rate(mu)
    if policy not in STATES:
        raise ValueError("policy must be C or M")
    if eventual and mu == 0:
        return ZERO  # Start has no R and no source of R.
    states = tuple(state for state in STATES[policy] if state != CURE)
    index = {state: i for i, state in enumerate(states)}
    matrix = [[ZERO] * len(states) for _ in states]
    rhs = [ZERO] * len(states)
    discount = ZERO if eventual else rates.z
    for state, i in index.items():
        events = transitions(policy, state, mu, rates)
        matrix[i][i] = discount + sum(events.values(), ZERO)
        for destination, rate in events.items():
            if destination is PROGRESSION:
                rhs[i] += rate
            elif destination != CURE:
                matrix[i][index[destination]] -= rate
    return solve_linear(matrix, rhs)[index[START]]


def fixed_horizon_interval(mu, policy, horizon, terms=350, rates=DEFAULT):
    """Rigorous Fraction enclosure of E[min(T,horizon)] by uniformization.

    This is a deterministic observation horizon; rates.z is not used. Let
    q be the largest outgoing rate, x=q*horizon, and s_j the probability of
    no progression after j steps of the substochastic matrix I+Q/q. Given
    n Poisson clock events before the horizon, the expected interval length
    is horizon/(n+1), so the conditional restricted mean is
    m_n=horizon*sum(s_j for j<=n)/(n+1), in [0,horizon]. Consequently the
    requested mean equals sum(x**n*m_n/n! for n>=0)/exp(x).

    With W and E the numerator and exponential series truncated at terms,
    the omitted exponential tail is at most
    B=(x**(terms+1)/(terms+1)!)/(1-x/(terms+2)). The return interval is
    [W/(E+B), (W+horizon*B)/E], intersected with [0,horizon]. Every operation
    is rational, so this bounds truncation AND arithmetic error exactly.
    The geometric bound requires terms+2 > q*horizon.
    """
    mu = mutation_rate(mu)
    if policy not in STATES:
        raise ValueError("policy must be C or M")
    if isinstance(horizon, bool) or not isinstance(horizon, (int, Fraction)):
        raise TypeError("horizon must be an exact integer or Fraction")
    horizon = F(horizon)
    if horizon < 0:
        raise ValueError("horizon must be nonnegative")
    if isinstance(terms, bool) or not isinstance(terms, int) or terms < 0:
        raise ValueError("terms must be a nonnegative integer")
    if horizon == 0:
        return ZERO, ZERO

    states = STATES[policy]
    index = {state: i for i, state in enumerate(states)}
    events = [transitions(policy, state, mu, rates) for state in states]
    outgoing = [sum(row.values(), ZERO) for row in events]
    q = max(outgoing)
    if q == 0:
        return horizon, horizon
    x = q * horizon
    if terms + 2 <= x:
        raise ValueError("terms + 2 must exceed the uniformization mean q*horizon")

    rows = []
    for i, row in enumerate(events):
        probabilities = {i: ONE - outgoing[i] / q}
        for destination, rate in row.items():
            if destination is not PROGRESSION:
                j = index[destination]
                probabilities[j] = probabilities.get(j, ZERO) + rate / q
        rows.append(tuple((j, p) for j, p in probabilities.items() if p))

    mass = [ZERO] * len(states)
    mass[index[START]] = ONE
    survival_sum = weighted_sum = exponential_sum = ZERO
    weight = ONE
    for n in range(terms + 1):
        survival_sum += sum(mass, ZERO)
        weighted_sum += weight * horizon * survival_sum / (n + 1)
        exponential_sum += weight
        if n < terms:
            following = [ZERO] * len(states)
            for probability, row in zip(mass, rows):
                for j, transition in row:
                    following[j] += probability * transition
            mass = following
        weight *= x / (n + 1)

    tail = weight / (ONE - x / (terms + 2))
    lower = weighted_sum / (exponential_sum + tail)
    upper = (weighted_sum + horizon * tail) / exponential_sum
    return max(ZERO, lower), min(horizon, upper)


def main_eventual_formula(mu, policy):
    """Corollary 2: only for DEFAULT; note the discontinuity at mu=0."""
    mu = mutation_rate(mu)
    if policy not in STATES:
        raise ValueError("policy must be C or M")
    if mu == 0:
        return ZERO
    if policy == "C":
        return F(103, 112)
    if policy == "M":
        return mu * (103 * mu**2 + 1018 * mu + 1824) / (
            112 * (mu + 2)**2 * (mu + 6))
    raise ValueError("policy must be C or M")


def isolate_root(polynomial, lower, upper, places=45):
    """Exact bisection. Uniqueness inside the bracket is proved in README."""
    lower, upper = F(lower), F(upper)
    if places < 1 or lower >= upper:
        raise ValueError("positive precision and ordered endpoints required")
    fl, fu = poly_value(polynomial, lower), poly_value(polynomial, upper)
    if fl * fu >= 0:
        raise ValueError("strict opposite endpoint signs required")
    tolerance = F(1, 10**places)
    while upper - lower > tolerance:
        mid = (lower + upper) / 2
        fm = poly_value(polynomial, mid)
        if fm == 0:
            return mid, mid
        if fl * fm < 0:
            upper = mid
        else:
            lower, fl = mid, fm
    return lower, upper


def decimal_interval(bounds, places=35):
    """Outward-rounded display; the returned exact Fraction bounds certify it."""
    result = []
    for value, rounding in zip(bounds, (ROUND_FLOOR, ROUND_CEILING)):
        with localcontext() as context:
            context.prec = max(80, places + 30)
            context.rounding = rounding
            number = Decimal(value.numerator) / Decimal(value.denominator)
            result.append(str(number.quantize(Decimal(1).scaleb(-places))))
    return result


def verify():
    """Fail loudly on an exact discrepancy; return a reproducibility report."""
    if not __debug__:
        raise RuntimeError("run without -O: exact assertions must remain enabled")
    psi = comparison_polynomial()
    assert psi == tuple(-v / 26 for v in MAIN_POLYNOMIAL)
    assert sign_variations(MAIN_POLYNOMIAL) == 2
    assert poly_value(MAIN_POLYNOMIAL, F(0)) > 0
    assert poly_value(MAIN_POLYNOMIAL, F(1)) < 0
    assert poly_value(MAIN_POLYNOMIAL, F(20)) > 0
    roots = [isolate_root(MAIN_POLYNOMIAL, 0, 1),
             isolate_root(MAIN_POLYNOMIAL, 1, 20)]
    grid = [F(0), F(1, 10**6), F(1, 100), F(1, 10), F(1, 4),
            F(1, 2), F(1), F(5, 2), F(10), F(25), F(100), F(10**6)]
    alternatives = [DEFAULT, Rates(11, 1, F(1, 10), F(1, 10), F(1, 10)),
                    Rates(1, 1, 1, 1, 1), Rates(7, 2, 4, 3, F(2, 5))]
    comparisons = 0
    for rates in alternatives:
        polynomial = comparison_polynomial(rates)
        assert len(polynomial) <= 4
        q = reduced_constants(rates)[0]
        leading = polynomial[3] if len(polynomial) == 4 else ZERO
        assert leading == 4 * rates.k * (1 - q) - rates.z
        for mu in grid:
            hc, hm = (formula_laplace(mu, p, rates) for p in ("C", "M"))
            for policy, closed in (("C", hc), ("M", hm)):
                assert closed == markov_laplace(mu, policy, rates)
                assert 0 <= closed <= 1
                comparisons += 1
            if mu:
                sign = lambda x: (x > 0) - (x < 0)
                assert sign(hc - hm) == sign(poly_value(polynomial, mu))
    for mu in grid:
        for policy in ("C", "M"):
            assert markov_laplace(mu, policy, eventual=True) == main_eventual_formula(mu, policy)
        if mu:
            difference = F(3) * mu**2 + 265 * mu + 618
            difference /= 28 * (mu + 2)**2 * (mu + 6)
            assert main_eventual_formula(mu, "C") - main_eventual_formula(mu, "M") == difference
            assert difference > 0
    samples = []
    for mu in (F(1, 4), F(5, 2), F(25)):
        difference = formula_control_time(mu, "M") - formula_control_time(mu, "C")
        samples.append({"mu": str(mu), "M_minus_C_control_time_exact": str(difference),
                        "M_minus_C_control_time_decimal": float(difference)})
    fixed_c = fixed_horizon_interval(F(5, 2), "C", F(1))
    fixed_m = fixed_horizon_interval(F(5, 2), "M", F(1))
    fixed_gap = (fixed_c[0] - fixed_m[1], fixed_c[1] - fixed_m[0])
    assert fixed_gap[0] > 0  # Opposite to Exp(1) at the same mutation rate.
    assert fixed_gap[1] - fixed_gap[0] < F(1, 10**75)
    return {"scope": "nonclinical computational mathematics; no treatment advice",
            "exact_generator_comparisons": comparisons,
            "polynomial_P_ascending": [int(v) for v in MAIN_POLYNOMIAL],
            "positive_root_intervals": [decimal_interval(bound) for bound in roots],
            "samples": samples,
            "fixed_horizon_1_mu_2.5_C_minus_M_interval": decimal_interval(fixed_gap),
            "ordinary_means": "both infinite: positive no-progression probability",
            "novelty": "specific result not located in targeted audit; priority unestablished"}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
