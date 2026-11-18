"""Independent, offline mathematical checks; nonclinical and unreviewed.

No medical advice. All declared synthetic fixtures are in math_scenarios.csv;
they are mathematical examples, never empirical observations. Numerical checks
support the implementation and do not replace the proofs in THEOREMS.md.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
_SPEC = importlib.util.spec_from_file_location("ch_prevalence_verified_model", HERE / "model.py")
assert _SPEC is not None and _SPEC.loader is not None
model = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(model)


def _number(row: dict, key: str) -> float:
    return float(row[key])


def _integer(row: dict, key: str) -> int:
    result = int(row[key])
    if result < 0:
        raise AssertionError(f"negative integer fixture field: {key}")
    return result


def _assert_close(actual: float, reference: float, tolerance: float, label: str) -> None:
    assert math.isfinite(actual) and math.isfinite(reference), label
    assert abs(actual-reference) <= tolerance, (label, actual, reference, tolerance)


def _simpson(function, lower: float, upper: float, intervals: int) -> float:
    """Verifier-owned composite Simpson rule, separate from production quadrature."""
    assert intervals >= 2 and intervals % 2 == 0
    step = (upper-lower)/intervals
    terms = [function(lower), function(upper)]
    for index in range(1, intervals):
        terms.append((4 if index % 2 else 2) * function(lower+index*step))
    return step * math.fsum(terms) / 3


def _uniformization_bounds(birth: float, death: float, age: float,
                           threshold: int, cap: int) -> tuple[float, float]:
    """Stopped finite-generator brackets with an absorbing overflow bucket.

    Paths are stopped on first exceeding cap, so retained states provide a
    lower bound; all overflow and uncomputed Poisson mass provide an upper
    bound. This computes no generating-function or geometric-tail formula.
    """
    assert cap >= threshold >= 1
    if age == 0 or birth+death == 0:
        probability = float(threshold == 1)
        return probability, probability
    rate = cap*(birth+death)
    poisson_mean = rate*age
    assert poisson_mean < 700, "fixture must avoid underflow of the first Poisson mass"
    distribution = [0.0]*(cap+2)
    distribution[1] = 1.0
    poisson_mass = math.exp(-poisson_mean)
    accumulated = [poisson_mass*x for x in distribution]
    omitted = math.inf
    for count in range(10000):
        # The omitted sum starts at count+1; its subsequent ratios are bounded
        # by poisson_mean/(count+2), giving a positive geometric upper bound.
        if count+2 > poisson_mean:
            next_mass = poisson_mass*poisson_mean/(count+1)
            omitted = next_mass/(1-poisson_mean/(count+2))
            if omitted < 1e-15:
                lower = math.fsum(accumulated[threshold:cap+1])
                upper = lower+accumulated[cap+1]+omitted
                return max(0.0, lower), min(1.0, upper)
        updated = [0.0]*(cap+2)
        updated[0] = distribution[0]
        updated[cap+1] = distribution[cap+1]
        for population in range(1, cap+1):
            mass = distribution[population]
            up = population*birth/rate
            down = population*death/rate
            updated[population-1] += mass*down
            updated[population] += mass*(1-up-down)
            updated[population+1] += mass*up
        distribution = updated
        poisson_mass *= poisson_mean/(count+1)
        for state, mass in enumerate(distribution):
            accumulated[state] += poisson_mass*mass
    raise AssertionError("uniformization exceeded its declared numerical budget")


def _independent_normalization(params: dict, fixture: dict) -> float:
    beta = params["shape"]
    delta = params["scale_per_year"]
    maximum = params["s_max_per_year"]
    power = math.ceil(_number(fixture, "power_multiplier")/beta)

    def integrand(u: float) -> float:
        if u == 0:
            return 0.0
        return (maximum*power*u**(power-1)
                * math.exp(-(maximum/delta)**beta * u**(power*beta)))

    return _simpson(integrand, 0.0, 1.0, _integer(fixture, "intervals"))


def _independent_intensity(params: dict, fixture: dict, normalization: float) -> float:
    """Gamma-coordinate Simpson integral, independent of log-fitness Gauss rule."""
    beta = params["shape"]
    delta = params["scale_per_year"]
    maximum = params["s_max_per_year"]
    age = _number(fixture, "age")
    v = _number(fixture, "v")
    upper_v = _number(fixture, "v_upper")
    a = 1/beta
    lower_y = (model.FITNESS_LOWER/delta)**beta
    upper_y = (maximum/delta)**beta
    log_factor = math.log(delta/beta)-math.log(normalization)

    def integrand(y: float) -> float:
        s = delta*math.exp(a*math.log(y))
        phi = math.expm1(s*age)/(2*params["n_tau_years"]*s)
        e1_difference = model.e1(v/((1-2*v)*phi))
        if upper_v < 0.5:
            e1_difference -= model.e1(upper_v/((1-2*upper_v)*phi))
        return math.exp(log_factor+(a-1)*math.log(y)-y)*e1_difference

    integral = _simpson(integrand, lower_y, upper_y, _integer(fixture, "intervals"))
    return 2*params["n_tau_years"]*params["mu_per_cell_year"]*integral


def _binomial_by_convolution(depth: int, probability: float, minimum: int) -> float:
    distribution = [1.0]
    for _ in range(depth):
        updated = [0.0]*(len(distribution)+1)
        for count, mass in enumerate(distribution):
            updated[count] += mass*(1-probability)
            updated[count+1] += mass*probability
        distribution = updated
    return math.fsum(distribution[minimum:])


def verify_math() -> dict:
    """Recompute deterministic numerical checks; raise AssertionError on failure."""
    with (HERE / "data" / "math_scenarios.csv").open(encoding="utf-8", newline="") as handle:
        fixtures = list(csv.DictReader(handle))
    assert fixtures, "missing mathematical fixtures"
    required_provenance = ("source", "source_url", "locator", "access_date", "reuse_basis")
    assert len({row["fixture_id"] for row in fixtures}) == len(fixtures)
    for row in fixtures:
        assert all(row[key] for key in required_provenance), row["fixture_id"]
        assert row["source"] == "original declared mathematical fixture; not empirical"
    parameters = model.load_parameters()
    normalization_rows = [row for row in fixtures if row["kind"] == "normalization"]
    assert len(normalization_rows) == 1
    normalizations = {
        key: _independent_normalization(params, normalization_rows[0])
        for key, params in parameters.items()
    }
    metrics = {
        "status": "passed",
        "fixture_rows": len(fixtures),
        "parameter_models": len(parameters),
        "ctmc_scenarios": 0,
        "ctmc_max_bracket_width": 0.0,
        "e1_max_relative_error": 0.0,
        "normalization_max_relative_error": 0.0,
        "intensity_max_absolute_error": 0.0,
        "bd_integral_max_absolute_error": 0.0,
        "age_mean_max_absolute_error": 0.0,
        "binomial_max_absolute_error": 0.0,
        "alias_max_absolute_error": 0.0,
        "series_bound_scenarios": 0,
    }
    incompatible_intervals = []
    handled = set()
    for row in fixtures:
        kind = row["kind"]
        ident = row["fixture_id"]
        tolerance = _number(row, "tolerance") if row["tolerance"] else 0.0
        if kind == "bd_ctmc":
            birth, death, age = (_number(row, key) for key in ("birth", "death", "age"))
            threshold, cap = _integer(row, "threshold"), _integer(row, "cap")
            lower, upper = _uniformization_bounds(birth, death, age, threshold, cap)
            exact = model.bd_tail(birth, death, age, threshold)
            assert lower-tolerance <= exact <= upper+tolerance, (ident, lower, exact, upper)
            metrics["ctmc_scenarios"] += 1
            metrics["ctmc_max_bracket_width"] = max(metrics["ctmc_max_bracket_width"], upper-lower)
            integrated = model.bd_integrated_tail(birth, death, age, threshold)
            numerical = _simpson(lambda t: model.bd_tail(birth, death, t, threshold),
                                 0, age, _integer(row, "intervals"))
            error = abs(integrated-numerical)
            _assert_close(integrated, numerical, tolerance, ident)
            metrics["bd_integral_max_absolute_error"] = max(metrics["bd_integral_max_absolute_error"], error)
        elif kind == "e1":
            x, shift = _number(row, "x"), _number(row, "tail_shift")
            numerical = _simpson(lambda y: math.exp(-x*math.exp(y)),
                                 0, math.log((x+shift)/x), _integer(row, "intervals"))
            remainder_bound = math.exp(-(x+shift))/(x+shift)
            exact = model.e1(x)
            _assert_close(exact, numerical, tolerance*exact+remainder_bound, ident)
            metrics["e1_max_relative_error"] = max(metrics["e1_max_relative_error"], abs(exact-numerical)/exact)
        elif kind == "normalization":
            for key, params in parameters.items():
                exact = model.dfe_normalization(params)
                relative = abs(exact-normalizations[key])/exact
                assert relative <= tolerance, (ident, key, relative)
                metrics["normalization_max_relative_error"] = max(metrics["normalization_max_relative_error"], relative)
        elif kind == "intensity":
            for key, params in parameters.items():
                exact = model.expected_clones(params, _number(row, "age"),
                                               _number(row, "v"), v_upper=_number(row, "v_upper"))
                numerical = _independent_intensity(params, row, normalizations[key])
                _assert_close(exact, numerical, tolerance, f"{ident}: {key}")
                metrics["intensity_max_absolute_error"] = max(metrics["intensity_max_absolute_error"], abs(exact-numerical))
                probability = model.prevalence(params, _number(row, "age"),
                                               _number(row, "v"), v_upper=_number(row, "v_upper"))
                _assert_close(probability, 1-math.exp(-exact), tolerance, ident)
        elif kind == "age_band":
            a, b, v = (_number(row, key) for key in ("age_lower", "age_upper", "v"))
            for key, params in parameters.items():
                lower, mean, upper = model.band_prediction(params, a, b, v)
                assert lower <= mean <= upper, (ident, key)
                numerical = _simpson(lambda age: model.prevalence(params, age, v), a, b,
                                     _integer(row, "intervals"))/(b-a)
                _assert_close(mean, numerical, tolerance, f"{ident}: {key}")
                metrics["age_mean_max_absolute_error"] = max(metrics["age_mean_max_absolute_error"], abs(mean-numerical))
        elif kind == "open_age_band":
            for params in parameters.values():
                lower, mean, upper = model.band_prediction(params, _number(row, "age_lower"), None,
                                                           _number(row, "v"))
                assert 0 <= lower <= 1 and mean is None and upper == 1
        elif kind == "bd_series":
            birth, death, age = (_number(row, key) for key in ("birth", "death", "age"))
            threshold, last = _integer(row, "threshold"), _integer(row, "series_last")
            s = birth-death
            q = (birth*math.expm1(s*age)/(birth*math.exp(s*age)-death)
                 if s else birth*age/(1+birth*age))
            partial = math.fsum(q**k/k for k in range(threshold, last+1))/birth
            bound = q**(last+1)/(birth*(last+1)*(1-q))
            exact = model.bd_integrated_tail(birth, death, age, threshold)
            assert partial-tolerance <= exact <= partial+bound+tolerance, ident
            metrics["series_bound_scenarios"] += 1
        elif kind == "critical_limit":
            birth, death, age = (_number(row, key) for key in ("birth", "death", "age"))
            threshold, epsilon = _integer(row, "threshold"), _number(row, "epsilon")
            critical = model.bd_tail(birth, death, age, threshold)
            for direction in (-1, 1):
                nearby = model.bd_tail(birth+direction*epsilon, death, age, threshold)
                _assert_close(nearby, critical, tolerance, ident)
        elif kind == "sharp_mixture":
            la, lb, scale, mix = (_number(row, key) for key in ("lambda_a", "lambda_b", "scale", "mix"))
            pa, pb = -math.expm1(-scale*la), -math.expm1(-scale*lb)
            observed = mix*pa+(1-mix)*pb
            lower, upper = model.multiplier_interval(la, lb, observed, observed)
            assert lower <= scale <= upper, ident
            _assert_close(-math.expm1(-lower*lb), observed, tolerance, ident)
            _assert_close(-math.expm1(-upper*la), observed, tolerance, ident)
        elif kind in ("incompatible_bins", "interval_boundary", "interval_unbounded"):
            values = [_number(row, key) for key in ("lambda_a", "lambda_b", "observed_lower", "observed_upper")]
            lower, upper = model.multiplier_interval(*values)
            if kind == "incompatible_bins":
                incompatible_intervals.append((lower, upper))
            elif kind == "interval_boundary":
                assert lower > upper, ident
            else:
                assert lower == 0 and upper == math.inf, ident
        elif kind == "alias":
            scale, eta, factor = (_number(row, key) for key in ("scale", "eta", "alias_factor"))
            for key in ("lambda_a", "lambda_b"):
                intensity = _number(row, key)
                first = -math.expm1(-eta*scale*intensity)
                second = -math.expm1(-(eta/factor)*(scale*factor)*intensity)
                _assert_close(first, second, tolerance, ident)
                metrics["alias_max_absolute_error"] = max(metrics["alias_max_absolute_error"], abs(first-second))
            # Exercise the production intensity, rather than checking only two
            # algebraically equal products disconnected from the implementation.
            for key, params in parameters.items():
                first_params = dict(params)
                second_params = dict(params)
                first_params["mu_per_cell_year"] *= scale
                second_params["mu_per_cell_year"] *= scale*factor
                first_intensity = model.expected_clones(first_params, _number(row, "age"), _number(row, "v"))
                second_intensity = model.expected_clones(second_params, _number(row, "age"), _number(row, "v"))
                first = -math.expm1(-eta*first_intensity)
                second = -math.expm1(-(eta/factor)*second_intensity)
                _assert_close(first, second, tolerance, f"{ident}: {key}")
                metrics["alias_max_absolute_error"] = max(metrics["alias_max_absolute_error"], abs(first-second))
        elif kind == "poisson":
            a, b = _number(row, "lambda_a"), _number(row, "lambda_b")
            combined = -math.expm1(-(a+b))
            complement = 1-math.exp(-a)*math.exp(-b)
            _assert_close(combined, complement, tolerance, ident)
        elif kind == "binomial":
            depth, v, minimum = _integer(row, "depth"), _number(row, "v"), _integer(row, "min_alt_reads")
            exact = model.binomial_detection(depth, v, minimum)
            independent = _binomial_by_convolution(depth, v, minimum)
            _assert_close(exact, independent, tolerance, ident)
            metrics["binomial_max_absolute_error"] = max(metrics["binomial_max_absolute_error"], abs(exact-independent))
        else:
            raise AssertionError(f"unknown fixture kind: {kind}")
        handled.add(ident)
    assert handled == {row["fixture_id"] for row in fixtures}
    assert len(incompatible_intervals) >= 2
    assert max(item[0] for item in incompatible_intervals) > min(item[1] for item in incompatible_intervals)
    return metrics


if __name__ == "__main__":
    print("Nonclinical, unreviewed mathematical verification; no medical advice.")
    print(json.dumps(verify_math(), indent=2, sort_keys=True))
