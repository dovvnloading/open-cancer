"""Original log-space forecast implementation for the reconstructed protocol."""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


A_VALUES = (0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
B_VALUES = (0.0, 0.25, 0.5, 1.0, 2.0)
TIE_TOLERANCE = 1e-12


class EmptySampleError(ValueError):
    """A declared intervention empties a required source sample."""


@dataclass
class Observations:
    early: np.ndarray
    late: np.ndarray
    early_dictionary_lineages: int
    retained_mass: np.ndarray
    original_totals: np.ndarray
    retained_totals: np.ndarray
    keep_mask: np.ndarray


@dataclass
class FoldFit:
    base_a: float
    full_a: float
    full_b: float
    base_inner_score: float
    full_inner_score: float
    base_log_prediction: np.ndarray
    full_log_prediction: np.ndarray
    base_score: float
    full_score: float
    gain: float
    untuned_own_score: float
    untuned_other_score: float


def validate_counts(counts):
    raw = np.asarray(counts)
    if raw.ndim != 2 or raw.shape[1] != 8 or not raw.shape[0]:
        raise ValueError("Expected nonempty barcode rows with eight sample columns")
    if raw.dtype.kind not in "iuf" or not np.all(np.isfinite(raw)):
        raise ValueError("Counts must be finite nonnegative integers")
    if np.any(raw < 0) or np.any(raw != np.floor(raw)):
        raise ValueError("Counts must be finite nonnegative integers")
    if np.any(raw >= 2.0**63):
        raise ValueError("Counts must fit signed int64")
    raw = raw.astype(np.int64, copy=False)
    if np.any(raw.sum(axis=0, dtype=np.float64) >= 2.0**63):
        raise ValueError("Sample totals must fit signed int64")
    totals = raw.sum(axis=0, dtype=np.int64)
    if np.any(totals == 0):
        raise EmptySampleError("A required source sample is empty")
    return raw, totals


def validate_floor(floor):
    if not math.isfinite(floor) or not 0 <= floor < 1:
        raise ValueError("Frequency floor must satisfy 0 <= floor < 1")


def collapse_counts(counts, keep_mask):
    """Keep a fixed P1 dictionary and one category for all other barcode rows."""
    return np.concatenate((counts[keep_mask], counts[~keep_mask].sum(axis=0)[None, :]), axis=0)


def prepare_observations(counts, frequency_floor=0.0):
    validate_floor(frequency_floor)
    raw, totals = validate_counts(counts)
    selected = np.where(raw / totals >= frequency_floor, raw, 0)
    selected_totals = selected.sum(axis=0, dtype=np.int64)
    if np.any(selected_totals == 0):
        raise EmptySampleError("Abundance threshold empties a required sample")
    keep = np.any(selected[:, ::2] > 0, axis=1)
    collapsed = collapse_counts(selected, keep)
    probability = (collapsed / selected_totals).T
    return Observations(
        early=probability[::2], late=probability[1::2],
        early_dictionary_lineages=int(keep.sum()),
        retained_mass=(selected_totals / totals).reshape(4, 2),
        original_totals=totals.reshape(4, 2),
        retained_totals=selected_totals.reshape(4, 2), keep_mask=keep,
    )


def unit_rows(values, name="probabilities"):
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or array.shape[0] != 4 or not array.shape[1]:
        raise ValueError(f"{name} must have shape (4,K) with K positive")
    if not np.all(np.isfinite(array)) or np.any(array < 0):
        raise ValueError(f"{name} must be finite and nonnegative")
    totals = array.sum(axis=1)
    if np.any(totals <= 0):
        raise EmptySampleError(f"{name} contains an empty sample")
    return array / totals[:, None]


def smooth_rows(values, delta):
    if not math.isfinite(delta) or not 0 < delta < 1:
        raise ValueError("Smoothing delta must satisfy 0 < delta < 1")
    p = unit_rows(values)
    return (1.0-delta)*p + delta/p.shape[1]


def log_prediction(log_own, log_donor, a, b):
    logits = a*log_own + b*log_donor
    peak = logits.max()
    return logits - (peak + np.log(np.exp(logits-peak).sum()))


def select_candidate(candidates):
    """Select (a,b,score), applying the recorded score and parameter tie rules."""
    if not candidates or any(not math.isfinite(row[2]) for row in candidates):
        raise ValueError("Candidate scores must be finite and nonempty")
    best = max(row[2] for row in candidates)
    tied = [row for row in candidates if best-row[2] <= TIE_TOLERANCE]
    return min(tied, key=lambda row: (row[1], abs(row[0]-1.0), row[0]))


def fit_outer_fold(early, late, target, delta=1e-4, donor_early=None):
    """Tune on other cultures; the held-out late distribution only scores the fit."""
    if target not in range(4):
        raise ValueError("Target must be a culture index 0..3")
    x = smooth_rows(early, delta)
    y = unit_rows(late, "late probabilities")
    donor_x = x if donor_early is None else smooth_rows(donor_early, delta)
    if x.shape != y.shape or x.shape != donor_x.shape:
        raise ValueError("Own, donor and late observations need the same support")
    log_x = np.log(x)
    training = [r for r in range(4) if r != target]
    # An outer target never enters any inner donor mean.
    inner_log_donors = {inner: np.log(donor_x[[r for r in training if r != inner]].mean(axis=0))
                        for inner in training}
    candidates = []
    for a in A_VALUES:
        for b in B_VALUES:
            inner_scores = [float(np.dot(y[r], log_prediction(log_x[r], inner_log_donors[r], a, b)))
                            for r in training]
            candidates.append((a, b, math.fsum(inner_scores)/len(training)))
    base_a, _, base_inner = select_candidate([row for row in candidates if row[1] == 0])
    full_a, full_b, full_inner = select_candidate(candidates)
    outer_donor_log = np.log(donor_x[training].mean(axis=0))
    base_log = log_prediction(log_x[target], outer_donor_log, base_a, 0)
    full_log = log_prediction(log_x[target], outer_donor_log, full_a, full_b)
    return FoldFit(
        base_a=base_a, full_a=full_a, full_b=full_b,
        base_inner_score=base_inner, full_inner_score=full_inner,
        base_log_prediction=base_log, full_log_prediction=full_log,
        base_score=float(np.dot(y[target], base_log)),
        full_score=float(np.dot(y[target], full_log)),
        gain=float(np.dot(y[target], full_log-base_log)),
        untuned_own_score=float(np.dot(y[target], log_x[target])),
        untuned_other_score=float(np.dot(y[target], outer_donor_log)),
    )
