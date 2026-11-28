#!/usr/bin/env python3
"""Recompute aggregate barcode forecasts from authenticated public source counts."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

from model import EmptySampleError, prepare_observations, fit_outer_fold


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TARGET_BUDGETS = (0.0, 0.001, 0.01, 0.05, 0.1)


def read_cache(path, line):
    from fetch_source import read_authenticated_cache
    return read_authenticated_cache(path, line)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fit_scalar_record(fit, target):
    return {
        "target_culture": f"DT{target+1}",
        "base_a": fit.base_a, "full_a": fit.full_a, "full_b": fit.full_b,
        "base_inner_mean_score": fit.base_inner_score,
        "full_inner_mean_score": fit.full_inner_score,
        "base_score": fit.base_score, "full_score": fit.full_score, "gain": fit.gain,
        "untuned_own_score": fit.untuned_own_score,
        "untuned_other_score": fit.untuned_other_score,
    }


def top_category_removal(late, fit):
    top = int(np.argmax(late))
    keep = np.arange(len(late)) != top
    target_mass = float(late[keep].sum())
    if target_mass <= 0:
        return {"available": False, "reason": "The highest category contains every target read"}
    new_target = late[keep]/target_mass
    def renormalize(logp):
        values = logp[keep]
        peak = values.max()
        return values-(peak+np.log(np.exp(values-peak).sum()))
    base, full = renormalize(fit.base_log_prediction), renormalize(fit.full_log_prediction)
    return {"available": True, "removed_target_mass": float(late[top]),
            "gain": float(np.dot(new_target, full-base)),
            "predictions_renormalized_without_refitting": True}


def scenario(counts, frequency_floor=0.0, delta=1e-4, with_bounds=False):
    try:
        obs = prepare_observations(counts, frequency_floor)
    except EmptySampleError as error:
        return {"available": False, "frequency_floor": frequency_floor,
                "delta": delta, "reason": str(error)}
    folds = []
    for target in range(4):
        fit = fit_outer_fold(obs.early, obs.late, target, delta)
        record = fit_scalar_record(fit, target)
        record["retained_early_mass"], record["retained_late_mass"] = obs.retained_mass[target].tolist()
        if with_bounds:
            from robustness import adverse_deletion_bound, tipping_mass
            contrasts = fit.full_log_prediction-fit.base_log_prediction
            record["fixed_predictor_target_only_bounds"] = [
                adverse_deletion_bound(obs.late[target], contrasts, epsilon)
                for epsilon in TARGET_BUDGETS]
            record["target_only_tipping"] = tipping_mass(obs.late[target], contrasts)
            record["highest_late_category_removal"] = top_category_removal(obs.late[target], fit)
        folds.append(record)
    mean_gain = math.fsum(row["gain"] for row in folds)/4
    return {
        "available": True, "frequency_floor": frequency_floor, "delta": delta,
        "early_dictionary_lineages": obs.early_dictionary_lineages,
        "dictionary_categories_including_other": obs.early.shape[1],
        "original_read_totals": obs.original_totals.tolist(),
        "retained_read_totals": obs.retained_totals.tolist(),
        "folds": folds, "mean_gain": mean_gain,
        "geometric_probability_ratio": math.exp(mean_gain),
    }


def decision(primary, required_floor):
    bounds = [next(row["lower_gain"] for row in fold["fixed_predictor_target_only_bounds"]
                   if row["epsilon"] == 0.01) for fold in primary["folds"]]
    conditions = {
        "mean_gain_gt_0_01": primary["mean_gain"] > .01,
        "all_primary_fold_gains_positive": all(row["gain"] > 0 for row in primary["folds"]),
        "required_0_001_abundance_floor_survives": bool(required_floor["available"]
            and required_floor["mean_gain"] > 0
            and all(row["gain"] >= 0 for row in required_floor["folds"])),
        "fixed_predictor_one_percent_target_bound_survives": all(value > 0 for value in bounds),
    }
    return {
        "verdict": "PASS" if all(conditions.values()) else "FAIL",
        "criterion": "Recorded descriptive predictive gate, not a biological significance test",
        "conditions": conditions,
        "mean_fixed_predictor_one_percent_lower_gain": math.fsum(bounds)/4,
        "fixed_predictor_one_percent_lower_gain_range": [min(bounds), max(bounds)],
    }


def analyze_line(counts, line, cache_path):
    primary = scenario(counts, with_bounds=True)
    abundance = [scenario(counts, floor) for floor in (1e-4, 1e-3)]
    smoothing = [scenario(counts, delta=delta) for delta in (1e-6, 1e-3)]
    return {
        "cell_line": line, "review_status": "unreviewed",
        "scope": "Known-result reconstruction; cell-line distribution prediction, not a resistance mechanism",
        "cache_sha256": sha256(cache_path),
        "source_rows": int(counts.shape[0]),
        "execution": {"created_at_utc": datetime.now(timezone.utc).isoformat(),
                      "python": sys.version.split()[0], "numpy": np.__version__},
        "primary": primary, "abundance_sensitivities": abundance,
        "smoothing_sensitivities": smoothing, "decision": decision(primary, abundance[-1]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=ROOT/".compute/barcode-source")
    parser.add_argument("--out", type=Path, default=ROOT/".compute/barcode-lineage-memory")
    parser.add_argument("--lines", nargs="+", choices=("HCTbc", "SW6bc"), default=("HCTbc", "SW6bc"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for line in args.lines:
        cache_path = args.cache_dir/f"{line}.npz"
        counts, _ = read_cache(cache_path, line)
        report = analyze_line(counts, line, cache_path)
        (args.out/f"{line}.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        print(json.dumps({"line": line, "mean_gain": report["primary"]["mean_gain"],
                          "fold_gains": [f["gain"] for f in report["primary"]["folds"]],
                          "decision": report["decision"]}), flush=True)


if __name__ == "__main__":
    main()
