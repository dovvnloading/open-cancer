#!/usr/bin/env python3
"""Independent fixed-seed numerical falsification attempts; stdlib only.

Nonclinical computational research.  These checks are not mathematical proofs,
biological validation, or medical advice.  This module deliberately imports no
code from verifier.py or analysis.py.  Event-level branching simulations,
ordinary double-precision RK4, simplex Monte Carlo, and Poisson likelihood tests
provide implementations distinct from the exact-arithmetic verifier.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from typing import Callable, Sequence


SEED = 19270617


def _rk4(
    rhs: Callable[[list[float]], list[float]],
    initial: Sequence[float],
    horizon: float,
    steps: int = 1600,
) -> list[float]:
    """Fixed-step fourth-order integration, independently implemented here."""
    y = list(initial)
    dt = horizon / steps
    for _ in range(steps):
        a = rhs(y)
        b = rhs([v + dt * x / 2 for v, x in zip(y, a)])
        c = rhs([v + dt * x / 2 for v, x in zip(y, b)])
        d = rhs([v + dt * x for v, x in zip(y, c)])
        y = [
            v + dt * (aa + 2 * bb + 2 * cc + dd) / 6
            for v, aa, bb, cc, dd in zip(y, a, b, c, d)
        ]
    return y


def clone_probability(
    rates: Sequence[tuple[float, float, float]],
    horizon: float,
    modes: str | Sequence[str] = "asymmetric",
) -> float:
    """Backward ODE success probability of one initially initiated cell.

    Each rates item is (birth, death, mutation).  An asymmetric mutation creates
    a next-type cell and retains its parent.  A converting mutation replaces
    its parent.  Final mutation is the first malignant event in either case.
    """
    n = len(rates)
    mode = [modes] * n if isinstance(modes, str) else list(modes)

    def rhs(q: list[float]) -> list[float]:
        out = []
        for i, (birth, death, mutation) in enumerate(rates):
            following = 1.0 if i + 1 == n else q[i + 1]
            value = (birth - death) * q[i] - birth * q[i] ** 2
            if mode[i] == "converting":
                value += mutation * (following - q[i])
            elif mode[i] == "asymmetric":
                value += mutation * (1 - q[i]) * following
            else:
                raise ValueError("unknown mutation convention")
            out.append(value)
        return out

    return _rk4(rhs, [0.0] * n, horizon)[0]


def _first_mark(
    rng: random.Random,
    rates: Sequence[tuple[float, float, float]],
    horizon: float,
    modes: str | Sequence[str],
) -> bool:
    """Exact Gillespie path up to first malignant mutation or censoring."""
    n = len(rates)
    mode = [modes] * n if isinstance(modes, str) else list(modes)
    cells = [1] + [0] * (n - 1)
    age = 0.0
    while sum(cells):
        weights = [count * sum(rate) for count, rate in zip(cells, rates)]
        total = sum(weights)
        if total == 0:
            return False
        age += rng.expovariate(total)
        if age > horizon:
            return False
        draw = rng.random() * total
        stage = 0
        while draw >= weights[stage]:
            draw -= weights[stage]
            stage += 1
        birth, death, mutation = rates[stage]
        event = rng.random() * (birth + death + mutation)
        if event < birth:
            cells[stage] += 1
        elif event < birth + death:
            cells[stage] -= 1
        elif stage + 1 == n:
            return True
        else:
            cells[stage + 1] += 1
            if mode[stage] == "converting":
                cells[stage] -= 1
    return False


def _terminal_cell_time(
    rng: random.Random,
    rates: Sequence[tuple[float, float, float]],
    horizon: float,
    modes: str,
) -> float:
    """Unmarked branching path integral of final premalignant cell count.

    The terminal mutation rate is ignored, so marks can be attached after the
    path has been drawn.  This is independent of the backward ODE integrator.
    """
    n = len(rates)
    cells = [1] + [0] * (n - 1)
    age = integral = 0.0
    while age < horizon and sum(cells):
        local = [
            (birth, death, 0.0 if i + 1 == n else mutation)
            for i, (birth, death, mutation) in enumerate(rates)
        ]
        weights = [count * sum(rate) for count, rate in zip(cells, local)]
        total = sum(weights)
        if total == 0:
            return integral + cells[-1] * (horizon - age)
        waiting = rng.expovariate(total)
        integral += cells[-1] * min(waiting, horizon - age)
        age += waiting
        if age >= horizon:
            break
        draw = rng.random() * total
        stage = 0
        while draw >= weights[stage]:
            draw -= weights[stage]
            stage += 1
        birth, death, mutation = local[stage]
        event = rng.random() * (birth + death + mutation)
        if event < birth:
            cells[stage] += 1
        elif event < birth + death:
            cells[stage] -= 1
        else:
            cells[stage + 1] += 1
            if modes == "converting":
                cells[stage] -= 1
    return integral


def _mean_hazard(
    net_rates: Sequence[float], amplitude: float, horizon: float
) -> float:
    """Linear chain response with unit transition factors and given product."""
    n = len(net_rates)

    def rhs(state: list[float]) -> list[float]:
        values = [net_rates[0] * state[0]]
        values.extend(
            net_rates[i] * state[i] + state[i - 1] for i in range(1, n)
        )
        values.append(amplitude * state[n - 1])
        return values

    return _rk4(rhs, [1.0] + [0.0] * n, horizon)[-1]


def _delayed_mean_hazard(
    net_rates: Sequence[float],
    amplitude: float,
    delay_rates: Sequence[float],
    horizon: float,
) -> float:
    """Forward biological means followed by independent exponential filters.

    Filters are integrated as additional observation states; negative poles
    are not inserted directly into the biological chain in this implementation.
    """
    n = len(net_rates)

    def rhs(state: list[float]) -> list[float]:
        if n:
            values = [net_rates[0] * state[0]]
            values.extend(
                net_rates[i] * state[i] + state[i - 1] for i in range(1, n)
            )
            values.append(amplitude * state[n - 1])
        else:
            values = [0.0]
        for i, rate in enumerate(delay_rates):
            values.append(rate * (state[n + i] - state[n + i + 1]))
        return values

    initial = ([1.0] + [0.0] * n) if n else [amplitude]
    initial += [0.0] * len(delay_rates)
    return _rk4(rhs, initial, horizon)[-1]


def _sample_summary(values: Sequence[float]) -> tuple[float, float]:
    n = len(values)
    average = math.fsum(values) / n
    variance = math.fsum((x - average) ** 2 for x in values) / (n - 1)
    return average, math.sqrt(variance / n)


def _poisson(rng: random.Random, mean: float) -> int:
    """Knuth sampling; all means in these checks are small."""
    if mean < 0:
        raise ValueError("Poisson means must be nonnegative")
    if mean == 0:
        return 0
    threshold = math.exp(-mean)
    product = 1.0
    count = 0
    while product > threshold:
        product *= rng.random()
        count += 1
    return count - 1


def run_checks(quick: bool = False) -> dict:
    """Run independent checks, raising AssertionError on any failed attempt.

    Six-standard-error Monte Carlo tolerances are specified before sampling.
    Seeds and sample sizes are fixed; a passing check is corroboration only.
    """
    rng = random.Random(SEED)
    paths = 4000 if quick else 16000
    output: dict = {
        "label": "nonclinical computational research; not medical advice",
        "seed": SEED,
        "paths_per_branching_case": paths,
        "event_level_branching": [],
        "onset_order": [],
    }

    # T1: independent event-level branching versus backward probabilities;
    # expansion is positive at every stage, with both mutation conventions.
    for n in (1, 2, 3):
        rates = [(0.20 + 0.04 * i, 0.12 + 0.02 * i, 0.65) for i in range(n)]
        for mode in ("asymmetric", "converting"):
            theoretical = clone_probability(rates, 2.5, mode)
            successes = sum(_first_mark(rng, rates, 2.5, mode) for _ in range(paths))
            observed = successes / paths
            tolerance = 6 * math.sqrt(theoretical * (1 - theoretical) / paths)
            assert abs(observed - theoretical) <= tolerance
            output["event_level_branching"].append(
                {
                    "mutation_steps_k": n + 1,
                    "mutation_convention": mode,
                    "ode_probability": theoretical,
                    "simulated_probability": observed,
                    "six_standard_errors": tolerance,
                }
            )
            # This numerical onset check is separate from rare-event sampling.
            small = 0.002
            leading = math.prod(rate[2] for rate in rates) * small**n / math.factorial(n)
            ratio = clone_probability(rates, small, mode) / leading
            assert abs(ratio - 1) < 0.003
            output["onset_order"].append(
                {"k": n + 1, "mode": mode, "q_over_leading_term": ratio}
            )

    # Established exact TSCE symmetry: q -> c q, lambda -> lambda/c.
    immigration, birth, death, mutation, scale = 0.7, 0.4, 0.2, 0.3, 1.2
    net = birth - death - mutation
    changed = (birth / scale, birth / scale - scale * mutation - net, scale * mutation)
    assert changed[1] >= 0
    residuals = []
    closed_form_residuals = []
    discriminant = math.sqrt(net * net + 4 * birth * mutation)
    for age in (0.25, 1.0, 3.0, 8.0):
        q = clone_probability([(birth, death, mutation)], age)
        transformed = clone_probability([changed], age)
        residuals.append(abs(immigration * q - immigration / scale * transformed))
        decay = math.exp(-discriminant * age)
        formula = 2 * mutation * (1 - decay) / (
            (discriminant - net) + (discriminant + net) * decay
        )
        closed_form_residuals.append(abs(q - formula))
    assert max(residuals) < 1e-10
    assert max(closed_form_residuals) < 1e-9
    dt = 0.001
    hazards = [
        immigration * clone_probability([(birth, death, mutation)], j * dt)
        for j in range(4)
    ]
    first = hazards[1] / dt
    second = (hazards[2] - 2 * hazards[1]) / dt**2
    third = (hazards[3] - 3 * hazards[2] + 3 * hazards[1]) / dt**3
    amplitude = immigration * mutation
    expected_derivatives = (amplitude, amplitude * net, amplitude * (net**2 - 2 * birth * mutation))
    for measured, expected in zip((first, second, third), expected_derivatives):
        assert abs(measured - expected) < 2e-4
    output["tsce_symmetry"] = {
        "max_hazard_residual": max(residuals),
        "max_closed_form_residual": max(closed_form_residuals),
        "forward_difference_derivatives": [first, second, third],
        "expected_derivatives": expected_derivatives,
    }

    # Rare-marking limit from unmarked path integrals.  No limiting formula is
    # used to generate these paths.  Both conventions give the same proof but
    # different intermediate net-growth coefficients.
    output["rare_marking_limit"] = []
    rates = [(0.25, 0.45, 0.6), (0.30, 0.50, 0.8), (0.20, 0.40, 0.7)]
    horizon, immigration = 3.0, 0.6
    for mode in ("asymmetric", "converting"):
        integrals = [_terminal_cell_time(rng, rates, horizon, mode) for _ in range(paths)]
        terminal = rates[-1][2]
        marked_time = [terminal * value for value in integrals]
        estimate, se = _sample_summary(marked_time)
        nets = [
            birth - death - (mutation if mode == "converting" and i + 1 < len(rates) else 0)
            for i, (birth, death, mutation) in enumerate(rates)
        ]
        product = immigration * math.prod(rate[2] for rate in rates)
        limiting = _mean_hazard(nets, product, horizon)
        assert abs(immigration * estimate - limiting) < 6 * immigration * se
        second_moment = math.fsum(x * x for x in marked_time) / paths
        finite = []
        previous = 0.0
        for epsilon in (0.1, 0.01, 0.001):
            value = immigration * math.fsum(-math.expm1(-epsilon * x) / epsilon for x in marked_time) / paths
            upper_gap = immigration * epsilon * second_moment / 2
            gap = immigration * estimate - value
            assert -1e-12 <= gap <= upper_gap + 1e-12
            assert value >= previous
            previous = value
            finite.append({"epsilon": epsilon, "hazard_mc": value, "gap_to_mc_limit": gap, "moment_bound_mc": upper_gap})
        output["rare_marking_limit"].append(
            {"mode": mode, "limiting_ode_hazard": limiting, "limiting_mc_hazard": immigration * estimate, "six_standard_errors": 6 * immigration * se, "finite_epsilon": finite}
        )

    # Exponential convolution as a simplex integral: uniform spacings implement
    # Dirichlet(1,...,1), with volume t^n/n!.  This includes negative, zero,
    # repeated, and positive net rates, independently of the ODE chain.
    output["spectral_limit"] = []
    simplex_paths = 8000 if quick else 32000
    for nets in ((-0.5,), (0.2, -0.4), (0.0, 0.35, -0.2), (0.1, 0.1, -0.3, 0.0)):
        n, age, amplitude = len(nets), 2.0, 0.3
        samples = []
        for _ in range(simplex_paths):
            cut = [0.0] + sorted(rng.random() * age for _ in range(n)) + [age]
            lengths = [cut[i + 1] - cut[i] for i in range(n + 1)]
            samples.append(math.exp(math.fsum(r * length for r, length in zip(nets, lengths))))
        mean, se = _sample_summary(samples)
        factor = amplitude * age**n / math.factorial(n)
        theoretical = _mean_hazard(nets, amplitude, age)
        reversed_value = _mean_hazard(list(reversed(nets)), amplitude, age)
        assert abs(theoretical - factor * mean) < 6 * factor * se + 1e-12
        assert abs(theoretical - reversed_value) < 1e-10
        output["spectral_limit"].append(
            {"net_rates": nets, "hazard_ode": theoretical, "hazard_simplex_mc": factor * mean, "six_standard_errors": 6 * factor * se, "permuted_hazard": reversed_value}
        )

    # Full observation spectra identify biological stages only when their net
    # rates are separated from negative delay poles.  Independently integrate
    # every allocation of two negative poles to biology versus observation.
    poles, observed_numerator, age = [0.25, 0.0, -0.4, -0.9], 0.3, 3.0
    allocation_values = []
    output["delay_pole_assignments"] = []
    for mask in range(4):
        delays = [rate for i, rate in enumerate((0.4, 0.9)) if mask & (1 << i)]
        biological = [0.25, 0.0] + [-rate for i, rate in enumerate((0.4, 0.9)) if not mask & (1 << i)]
        biological_amplitude = observed_numerator / math.prod(delays)
        value = _delayed_mean_hazard(biological, biological_amplitude, delays, age)
        allocation_values.append(value)
        output["delay_pole_assignments"].append(
            {"k": len(biological) + 1, "biological_net_rates": biological, "observation_delay_rates": delays, "biological_amplitude": biological_amplitude, "observed_hazard": value}
        )
    assert max(allocation_values) - min(allocation_values) < 1e-9
    assert abs(allocation_values[0] - _mean_hazard(poles, observed_numerator, age)) < 1e-10
    # k=1 boundary: all poles are delays and the biological hazard is constant.
    no_intermediate = _delayed_mean_hazard([], 0.3 / (0.4 * 0.9), [0.4, 0.9], age)
    two_intermediates = _mean_hazard([-0.4, -0.9], 0.3, age)
    assert abs(no_intermediate - two_intermediates) < 1e-10

    # Fast *converting* precursor, with arbitrary downstream asymmetric
    # branching.  Numerical quadrature independently checks integrated deficit.
    downstream = [(0.35, 0.1, 0.5), (0.25, 0.15, 0.6)]
    horizon, immigration = 3.0, 0.4
    output["fast_precursor"] = []
    grid = 120
    original = [immigration * clone_probability(downstream, horizon * j / grid) for j in range(grid + 1)]
    for speed in (2.0, 10.0, 40.0):
        precursor = [(0.0, 0.0, speed)] + downstream
        modes = ["converting", "asymmetric", "asymmetric"]
        delayed = [immigration * clone_probability(precursor, horizon * j / grid, modes) for j in range(grid + 1)]
        deficits = [a - b for a, b in zip(original, delayed)]
        integral = horizon / grid * (sum(deficits) - (deficits[0] + deficits[-1]) / 2)
        assert min(deficits) > -1e-11
        assert integral <= original[-1] / speed + 1e-7
        probability = delayed[-1] / immigration
        observed = sum(_first_mark(rng, precursor, horizon, modes) for _ in range(paths)) / paths
        diagnosed = 0
        for _ in range(paths):
            lag = rng.expovariate(speed)
            if lag < horizon:
                diagnosed += _first_mark(rng, downstream, horizon - lag, "asymmetric")
        observed_diagnosis = diagnosed / paths
        tolerance = 6 * math.sqrt(probability * (1 - probability) / paths)
        assert abs(observed - probability) <= tolerance
        assert abs(observed_diagnosis - probability) <= tolerance
        assert all(b >= a - 1e-11 for a, b in zip(delayed, delayed[1:]))
        output["fast_precursor"].append(
            {"rate": speed, "integrated_hazard_deficit": integral, "bound_h_T_over_L": original[-1] / speed, "ode_probability": probability, "simulated_precursor_probability": observed, "simulated_clone_diagnosis_probability": observed_diagnosis, "six_standard_errors": tolerance}
        )

    # A fixed-seed positive-rate example checks the sharp bounded-multiplier
    # feasibility criterion independently through intersecting amplitude ranges.
    output["bounded_observation_multiplier"] = []
    for _ in range(8):
        shape = [rng.uniform(0.1, 3) for _ in range(5)]
        observed = [q * rng.uniform(0.6, 1.8) for q in shape]
        ratios = [y / q for y, q in zip(observed, shape)]
        threshold = math.sqrt(max(ratios) / min(ratios))
        central_amplitude = math.sqrt(max(ratios) * min(ratios))
        multipliers = [r / central_amplitude for r in ratios]
        assert min(multipliers) >= 1 / threshold - 1e-12
        assert max(multipliers) <= threshold + 1e-12
        below = (1 + threshold) / 2
        feasible_below = max(r / below for r in ratios) <= min(r * below for r in ratios)
        assert not feasible_below
        output["bounded_observation_multiplier"].append(
            {"minimal_symmetric_factor": threshold, "attaining_amplitude": central_amplitude, "factor_below_threshold_is_feasible": feasible_below}
        )

    # Finite-bin Poisson decision experiment.  A(t)=a*t is a k=2 limiting
    # hazard; the fast-stage alternative is a*(t-(1-exp(-L*t))/L), k=3.
    amplitude, speed, exposure_density = 0.3, 8.0, 20.0
    bins = [(0.0, 1.0), (1.0, 2.0), (2.0, 4.0)]
    base, alternative = [], []
    for left, right in bins:
        integrated = amplitude * (right**2 - left**2) / 2
        deficit = amplitude * ((right - left) / speed + (math.exp(-speed * right) - math.exp(-speed * left)) / speed**2)
        base.append(exposure_density * integrated)
        alternative.append(exposure_density * (integrated - deficit))
    distance = math.fsum((math.sqrt(a) - math.sqrt(b)) ** 2 for a, b in zip(base, alternative))
    hellinger_lower = (1 - math.sqrt(-math.expm1(-distance))) / 2
    coupling_lower = math.exp(-math.fsum(a - b for a, b in zip(base, alternative))) / 2
    trials = 5000 if quick else 25000
    log_ratios = [math.log(a / b) for a, b in zip(base, alternative)]
    offset = math.fsum(a - b for a, b in zip(base, alternative))
    errors = []
    for truth, means in enumerate((base, alternative)):
        count = 0
        for _ in range(trials):
            statistic = math.fsum(_poisson(rng, mean) * weight for mean, weight in zip(means, log_ratios)) - offset
            predict = 0 if statistic >= 0 else 1
            count += predict != truth
        errors.append(count / trials)
    bayes_error = sum(errors) / 2
    # Independent halves have variances <= 1/(4*trials), so this is a
    # conservative six-standard-error bound for their mean.
    tolerance = 6 / math.sqrt(8 * trials)
    assert bayes_error + tolerance >= max(hellinger_lower, coupling_lower)
    output["poisson_decision"] = {
        "bin_means_k2": base,
        "bin_means_k3": alternative,
        "hellinger_lower_bound": hellinger_lower,
        "coupling_lower_bound": coupling_lower,
        "class_conditional_errors": errors,
        "equal_prior_error": bayes_error,
        "six_standard_error_upper_bound": tolerance,
        "trials_per_hypothesis": trials,
    }
    output["passed"] = True
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="use the smaller fixed test sample sizes")
    args = parser.parse_args()
    print(json.dumps(run_checks(quick=args.quick), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
