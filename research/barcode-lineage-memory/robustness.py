"""Sharp target-only deletion sensitivity for fixed forecast log contrasts.

This is an optimization over probability mass, not independent cells, an assay
error estimate, or a predictor-contamination bound. Arithmetic is floating
point; the proof of sharpness concerns the ideal exact optimization.
"""
from __future__ import annotations

import math
from typing import Iterable


def _inputs(probabilities: Iterable[float], contrasts: Iterable[float]):
    p = [float(value) for value in probabilities]
    d = [float(value) for value in contrasts]
    if not p or len(p) != len(d):
        raise ValueError("Nonempty matching probability and contrast vectors are required")
    if any(not math.isfinite(value) or value < 0 for value in p):
        raise ValueError("Probabilities must be finite and nonnegative")
    if any(not math.isfinite(value) for value in d):
        raise ValueError("Contrasts must be finite")
    total = math.fsum(p)
    if not math.isfinite(total) or total <= 0:
        raise ValueError("Probabilities need finite positive total mass")
    p = [value / total for value in p]
    return p, d


def adverse_deletion_bound(probabilities, contrasts, epsilon, *, return_witness=False):
    """Return the minimum retained mean after deleting at most epsilon mass.

    Write e_i for removed mass. For a fixed mass budget, the remaining numerator
    is sum((p_i-e_i)*d_i), so minimizing it is a fractional-knapsack problem:
    remove the largest d first, splitting the boundary category if necessary.
    The optimal retained mean is nonincreasing in the removal budget, hence
    the optimum for a budget of *at most* epsilon is attained using epsilon.
    Uniform contrasts can attain the same value at every removal budget.

    The optional witness is for numerical checks, not publication of barcode
    measurements. The default result contains aggregate values only.
    """
    p, d = _inputs(probabilities, contrasts)
    epsilon = float(epsilon)
    if not math.isfinite(epsilon) or not 0 <= epsilon < 1:
        raise ValueError("epsilon must be finite with 0 <= epsilon < 1")
    original = math.fsum(mass * value for mass, value in zip(p, d))
    support = [index for index, mass in enumerate(p) if mass > 0]
    constant = all(d[index] == d[support[0]] for index in support)
    if constant and not return_witness:
        lower = d[support[0]]
        return {
            "epsilon": epsilon, "original_gain": original, "lower_gain": lower,
            "removed_mass": epsilon, "retained_mass": 1 - epsilon,
            "retained_weighted_sum": (1 - epsilon) * lower,
            "scope": "fixed_predictions_target_mass_only",
        }
    retained = p.copy()
    deleted = [0.0] * len(p) if return_witness else None
    remaining = epsilon
    for index in sorted(support, key=lambda index: (-d[index], index)):
        if remaining <= 0:
            break
        amount = min(p[index], remaining)
        retained[index] = p[index] - amount
        if deleted is not None:
            deleted[index] = amount
        remaining -= amount
    if remaining > 2e-12:
        raise ArithmeticError("Requested removal mass exceeded normalized support")
    numerator = math.fsum(mass * value for mass, value in zip(retained, d))
    lower = numerator / (1 - epsilon)
    if constant:
        lower = d[support[0]]
        numerator = (1 - epsilon) * lower
    result = {
        "epsilon": epsilon, "original_gain": original, "lower_gain": lower,
        "removed_mass": epsilon, "retained_mass": 1 - epsilon,
        "retained_weighted_sum": numerator,
        "scope": "fixed_predictions_target_mass_only",
    }
    if return_witness:
        result.update(deleted_probabilities=deleted, retained_probabilities=retained)
    return result


def tipping_mass(probabilities, contrasts):
    """First deletion budget attaining a nonpositive retained contrast, if any.

    For positive initial gain, integrate positive contrasts in descending
    order until the removed weighted contrast equals the original gain. If
    the occupied support has strictly positive contrasts, no budget below one
    reaches zero. A zero minimum can attain zero once all positive contrasts
    have been removed. Extremely degenerate floating-point roots that round
    to one are explicitly unresolved rather than reported as attainable.
    """
    p, d = _inputs(probabilities, contrasts)
    original = math.fsum(mass * value for mass, value in zip(p, d))
    support = [index for index, mass in enumerate(p) if mass > 0]
    minimum = min(d[index] for index in support)
    common = {"original_gain": original, "limiting_lower_gain": minimum,
              "scope": "fixed_predictions_target_mass_only"}
    if original <= 0:
        return dict(common, status="already_nonpositive", tipping_mass=0.0)
    if minimum > 0:
        return dict(common, status="unreachable_below_one", tipping_mass=None)
    if minimum == 0:
        root = math.fsum(p[index] for index in support if d[index] > 0)
        if root < 1:
            return dict(common, status="attained", tipping_mass=root)
        return dict(common, status="unresolved_at_float_precision", tipping_mass=None)
    gain_left = original
    removed_masses = []
    for index in sorted(support, key=lambda index: (-d[index], index)):
        if d[index] <= 0:
            break
        contribution = p[index] * d[index]
        if contribution >= gain_left:
            root = math.fsum((*removed_masses, gain_left / d[index]))
            if root < 1:
                return dict(common, status="attained", tipping_mass=max(0.0, root))
            return dict(common, status="unresolved_at_float_precision", tipping_mass=None)
        gain_left = math.fsum((gain_left, -contribution))
        removed_masses.append(p[index])
    return dict(common, status="unresolved_at_float_precision", tipping_mass=None)


# Descriptive aliases retained for users who call the bound by its target scope.
target_deletion_bound = adverse_deletion_bound
target_only_tipping_mass = tipping_mass
