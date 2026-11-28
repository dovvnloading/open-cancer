#!/usr/bin/env python3
"""SYNTHETIC probability examples: target checks do not establish predictor integrity.

No biological observations are read or fitted. See PREDICTOR_CONTAMINATION.md
for the limiting argument, a before-smoothing construction, and scope limits.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def _finite(value, name):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _probability(values, name):
    values = [_finite(value, name) for value in values]
    if not values or any(value <= 0 for value in values):
        raise ValueError(f"{name} must contain strictly positive probabilities")
    total = math.fsum(values)
    if not math.isfinite(total):
        raise ValueError(f"{name} total must be finite")
    return [value / total for value in values]


def product_log_prediction(own, donor, a=1.0, b=0.25):
    """Normalized geometric forecast, evaluated without multiplying tiny powers."""
    own, donor = _probability(own, "own"), _probability(donor, "donor")
    if len(own) != len(donor):
        raise ValueError("Own and donor vectors must have matching lengths")
    a, b = _finite(a, "a"), _finite(b, "b")
    if a < 0 or b < 0:
        raise ValueError("Forecast exponents must be nonnegative")
    weights = [a * math.log(x) + b * math.log(z) for x, z in zip(own, donor)]
    if any(not math.isfinite(value) for value in weights):
        raise ValueError("Log weights exceeded floating-point range")
    largest = max(weights)
    log_sum = largest + math.log(math.fsum(math.exp(value-largest) for value in weights))
    return [value-log_sum for value in weights]


def fixed_donor_mixture_score_bound(epsilon, minimum_probability, b=0.25):
    """Uniform score-change bound for fixed own predictor/exponents and mixture.

    If min(z) >= m > 0 and z_e=(1-e)z+e*v with v a probability vector, then
    abs(E_p log(Q(z_e)/Q(z))) <= b*log(1+e/((1-e)*m)). This does not cover
    parameter refitting or simultaneous changes to the own predictor.
    """
    epsilon = _finite(epsilon, "epsilon")
    minimum_probability = _finite(minimum_probability, "minimum_probability")
    b = _finite(b, "b")
    if not 0 <= epsilon < 1:
        raise ValueError("epsilon must satisfy 0 <= epsilon < 1")
    if not 0 < minimum_probability <= 1 or b < 0:
        raise ValueError("A positive probability floor and nonnegative b are required")
    if epsilon == 0 or b == 0:
        return 0.0
    # Evaluate log(1+e/((1-e)*m)) without overflowing its probability ratio.
    log_denominator = math.log1p(-epsilon) + math.log(minimum_probability)
    log_ratio = math.log(epsilon) - log_denominator
    log_one_plus_ratio = (log_ratio + math.log1p(math.exp(-log_ratio))
                          if log_ratio > 0 else math.log1p(math.exp(log_ratio)))
    return b * log_one_plus_ratio


def counterexample(*, w=0.2, m=1e-8, epsilon=1e-4, b=0.25):
    """An explicit synthetic three-category construction, not an estimated assay."""
    w, m = _finite(w, "w"), _finite(m, "m")
    epsilon, b = _finite(epsilon, "epsilon"), _finite(b, "b")
    if not 0 < w < 1:
        raise ValueError("w must satisfy 0 < w < 1")
    if not 0 < m < min(w, (1-w)/2):
        raise ValueError("m must satisfy 0 < m < min(w, (1-w)/2)")
    if not 0 < epsilon < 1 or not 0 < b < 1:
        raise ValueError("epsilon and b must lie strictly between zero and one")
    own = (w, m, 1-w-m)
    donor = (m, 1-2*m, m)
    contaminated = ((1-epsilon)*m+epsilon, (1-epsilon)*(1-2*m), (1-epsilon)*m)
    clean_log = product_log_prediction(own, donor, 1, b)
    contaminated_log = product_log_prediction(own, contaminated, 1, b)
    clean_gain, contaminated_gain = clean_log[0]-math.log(w), contaminated_log[0]-math.log(w)
    return {
        "label": "SYNTHETIC_PROBABILITY_COUNTEREXAMPLE", "a": 1.0, "b": b,
        "w": w, "m": m, "epsilon": epsilon, "own": list(own),
        "clean_donor": list(donor), "contaminated_donor": list(contaminated),
        "target": [1.0, 0.0, 0.0],
        # Exact symbolic value for this mixture, avoiding cancellation when
        # subtracting two nearly unit probabilities at very small epsilon.
        "donor_total_variation": epsilon * (1-m),
        "clean_target_prediction": math.exp(clean_log[0]),
        "contaminated_target_prediction": math.exp(contaminated_log[0]),
        "base_target_prediction": w, "clean_gain": clean_gain,
        "contaminated_gain": contaminated_gain,
        "contaminated_minus_clean_score": contaminated_log[0]-clean_log[0],
        "fixed_donor_mixture_score_bound": fixed_donor_mixture_score_bound(epsilon, m, b),
        "target_only_lower_gain_every_budget_below_one": contaminated_gain,
        "target_only_tipping_status": "unreachable_below_one" if contaminated_gain > 0 else "already_nonpositive",
        "m_to_zero_gain_limit_at_fixed_w_epsilon_b": -math.log(w),
    }


def before_smoothing_example(*, w=0.2, delta=1e-8, epsilon=1e-4, b=0.25):
    """SYNTHETIC contamination in raw probabilities followed by uniform smoothing."""
    w, delta = _finite(w, "w"), _finite(delta, "delta")
    epsilon, b = _finite(epsilon, "epsilon"), _finite(b, "b")
    if not 0 < w < 1 or not 0 < delta < 1 or not 0 < epsilon < 1 or not 0 < b < 1:
        raise ValueError("w, delta, epsilon and b must all lie strictly between zero and one")
    own_raw, donor_raw = (w, 0.0, 1-w), (0.0, 1.0, 0.0)
    dirty_raw = (epsilon, 1-epsilon, 0.0)
    smooth = lambda values: [(1-delta)*value+delta/3 for value in values]
    own, donor, dirty = smooth(own_raw), smooth(donor_raw), smooth(dirty_raw)
    clean_log, dirty_log = product_log_prediction(own, donor, 1, b), product_log_prediction(own, dirty, 1, b)
    return {"label": "SYNTHETIC_BEFORE_SMOOTHING_COUNTEREXAMPLE", "w": w,
            "delta": delta, "epsilon": epsilon, "b": b,
            "base_target_prediction": own[0],
            "clean_target_prediction": math.exp(clean_log[0]),
            "contaminated_target_prediction": math.exp(dirty_log[0]),
            "clean_gain": clean_log[0]-math.log(own[0]),
            "contaminated_gain": dirty_log[0]-math.log(own[0]),
            "smoothed_donor_total_variation": (1-delta)*epsilon}


def default_examples():
    return [counterexample(m=floor) for floor in (1e-3, 1e-5, 1e-8, 1e-12, 1e-20)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"scope": "Original SYNTHETIC probability examples; no empirical or clinical data",
              "examples": default_examples()}
    rendered = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
