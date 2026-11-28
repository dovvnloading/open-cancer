#!/usr/bin/env python3
"""Recompute declared exploratory interventions and a conditional temporal diagnostic."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

from model import (EmptySampleError, collapse_counts, fit_outer_fold,
                   prepare_observations, validate_counts, validate_floor)
from analysis import fit_scalar_record, read_cache, sha256


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MODES = ("early_only", "late_only", "both", "donor_feature_only")
FLOORS = (1e-5, 1e-4, 1e-3)
WEIGHTS = tuple(i/20 for i in range(21))


def fixed_dictionary_intervention(counts, frequency_floor, mode):
    """Threshold individual barcode entries before collapse, keeping primary support."""
    if mode not in MODES:
        raise ValueError("Unknown intervention mode")
    validate_floor(frequency_floor)
    raw, totals = validate_counts(counts)
    original = prepare_observations(raw)
    thresholded = np.where(raw/totals >= frequency_floor, raw, 0)
    selected = raw.copy()
    if mode in ("early_only", "both"):
        selected[:, ::2] = thresholded[:, ::2]
    if mode in ("late_only", "both"):
        selected[:, 1::2] = thresholded[:, 1::2]
    selected_totals = selected.sum(axis=0, dtype=np.int64)
    if np.any(selected_totals == 0):
        raise EmptySampleError("Intervention empties a required own-early or late sample")
    probability = (collapse_counts(selected, original.keep_mask)/selected_totals).T
    donor_early = None
    donor_mass = selected_totals[::2]/totals[::2]
    if mode == "donor_feature_only":
        donor_counts = thresholded[:, ::2]
        donor_totals = donor_counts.sum(axis=0, dtype=np.int64)
        if np.any(donor_totals == 0):
            raise EmptySampleError("Intervention empties a required donor sample")
        donor_early = (collapse_counts(donor_counts, original.keep_mask)/donor_totals).T
        donor_mass = donor_totals/totals[::2]
    return {"early": probability[::2], "late": probability[1::2],
            "donor_early": donor_early, "retained_mass": (selected_totals/totals).reshape(4, 2),
            "donor_mass": donor_mass,
            "early_dictionary_lineages": original.early_dictionary_lineages}


def ablation(counts, floor, mode, delta=1e-4):
    try:
        data = fixed_dictionary_intervention(counts, floor, mode)
    except EmptySampleError as error:
        return {"mode": mode, "floor": floor, "available": False, "reason": str(error)}
    folds = []
    for target in range(4):
        fit = fit_outer_fold(data["early"], data["late"], target, delta, data["donor_early"])
        record = fit_scalar_record(fit, target)
        record["retained_early_mass"], record["retained_late_mass"] = data["retained_mass"][target].tolist()
        # This field is indexed by SOURCE CULTURE, not the mean across outer donors.
        record["donor_feature_mass_retained"] = float(data["donor_mass"][target])
        folds.append(record)
    return {"mode": mode, "floor": floor, "available": True, "delta": delta,
            "early_dictionary_lineages": data["early_dictionary_lineages"],
            "dictionary_fixed": True, "folds": folds,
            "mean_gain": math.fsum(row["gain"] for row in folds)/4}


def half_parity(barcodes):
    return np.fromiter((hashlib.sha256(str(value).encode("ascii")).digest()[-1] % 2
                        for value in barcodes), dtype=np.int8, count=len(barcodes))


def normalized_profile(values, delta=1e-4):
    total = float(np.sum(values))
    if total <= 0:
        raise EmptySampleError("A selected-half home profile has no observations")
    return (1-delta)*np.asarray(values, dtype=float)/total + delta/len(values)


def temporal_diagnostic(counts, barcodes, delta=1e-4):
    raw, totals = validate_counts(counts)
    identifiers = np.asarray(barcodes)
    if identifiers.ndim != 1 or len(identifiers) != len(raw):
        raise ValueError("Barcodes must match the source rows")
    early = raw[:, ::2]/totals[::2]
    home_results = []
    for home in range(4):
        others = [r for r in range(4) if r != home]
        chosen = (early[:, home] >= 1e-3) & np.all(early[:, others] < 1e-4, axis=1)
        indices = np.flatnonzero(chosen)
        parity = half_parity(identifiers[chosen])
        halves = [indices[parity == half] for half in (0, 1)]
        home_record = {"home_culture": f"DT{home+1}", "selected_barcodes": len(indices),
                       "half_selected_barcodes": [len(idx) for idx in halves]}
        if any(len(idx) < 10 for idx in halves):
            home_record.update(available=False, reason="Fewer than ten selected barcodes in a hash half",
                               receiver_half_results=[])
            home_results.append(home_record)
            continue
        home_record["available"] = True
        receiver_results = []
        for receiver in others:
            profiles = []
            failure = None
            for idx in halves:
                observed = raw[idx, 2*receiver].astype(float)
                if observed.sum() <= 0:
                    failure = "A required receiver half has no observed reads"
                    break
                try:
                    p1 = normalized_profile(raw[idx, 2*home], delta)
                    p2 = normalized_profile(raw[idx, 2*home+1], delta)
                except EmptySampleError as error:
                    failure = str(error)
                    break
                profiles.append((observed/observed.sum(), p1, p2))
            for test_half in (0, 1):
                train_half = 1-test_half
                record = {"receiver_culture": f"DT{receiver+1}", "test_half": test_half,
                          "train_half": train_half,
                          "receiver_test_read_count": int(raw[halves[test_half], 2*receiver].sum()),
                          "receiver_test_occupied_barcodes": int(np.count_nonzero(raw[halves[test_half], 2*receiver])),
                          "receiver_train_read_count": int(raw[halves[train_half], 2*receiver].sum()),
                          "receiver_train_occupied_barcodes": int(np.count_nonzero(raw[halves[train_half], 2*receiver]))}
                if failure:
                    record.update(available=False, reason=failure)
                else:
                    train_target, train_p1, train_p2 = profiles[train_half]
                    train_scores = [float(np.dot(train_target, np.log((1-w)*train_p1+w*train_p2)))
                                    for w in WEIGHTS]
                    best = max(train_scores)
                    selected_weight = next(w for w, value in zip(WEIGHTS, train_scores)
                                           if best-value <= 1e-12)
                    test_target, test_p1, test_p2 = profiles[test_half]
                    base_score = float(np.dot(test_target, np.log(test_p1)))
                    mixture_score = float(np.dot(test_target, np.log((1-selected_weight)*test_p1+selected_weight*test_p2)))
                    record.update(available=True, selected_weight=selected_weight,
                                  base_score=base_score, mixture_score=mixture_score,
                                  gain=mixture_score-base_score)
                receiver_results.append(record)
        available = [row for row in receiver_results if row["available"]]
        home_record["receiver_half_results"] = receiver_results
        home_record["available_receiver_halves"] = len(available)
        # A home is not averaged over a favorable subset when a comparison is absent.
        if len(available) == 6:
            home_record["mean_gain"] = math.fsum(row["gain"] for row in available)/6
        else:
            home_record["available"] = False
            home_record["reason"] = "At least one receiver/half comparison is unavailable"
        home_results.append(home_record)
    available_homes = [home for home in home_results if home["available"]]
    return {"diagnostic_not_prospective_forecast": True,
            "score_unit": "nats per receiver read conditional on the selected barcode half",
            "mixture_weights_are_not_contamination_fractions": True,
            "home_results": home_results, "available_home_cultures": len(available_homes),
            "equal_home_mean_gain": (math.fsum(home["mean_gain"] for home in available_homes)/len(available_homes)
                                     if available_homes else None)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=ROOT/".compute/barcode-source")
    parser.add_argument("--out", type=Path, default=ROOT/".compute/barcode-lineage-memory")
    parser.add_argument("--lines", nargs="+", choices=("HCTbc", "SW6bc"), default=("HCTbc", "SW6bc"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for line in args.lines:
        cache_path = args.cache_dir/f"{line}.npz"
        counts, barcodes = read_cache(cache_path, line)
        report = {"cell_line": line, "review_status": "unreviewed", "exploratory_A2": True,
                  "reconstruction_with_known_results": True, "cache_sha256": sha256(cache_path),
                  "ablations": [ablation(counts, floor, mode) for floor in FLOORS for mode in MODES],
                  "temporal_diagnostic": temporal_diagnostic(counts, barcodes)}
        (args.out/f"{line}_exploratory.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        print(json.dumps({"line": line, "ablations": [{k: value for k, value in row.items()
                             if k in ("mode", "floor", "available", "mean_gain")}
                            for row in report["ablations"]],
                          "temporal_mean_gain": report["temporal_diagnostic"]["equal_home_mean_gain"]}), flush=True)


if __name__ == "__main__":
    main()
