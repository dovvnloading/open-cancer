#!/usr/bin/env python3
"""Independent probability-space reconstruction; never imports forecast code.

The source authentication and the mathematical implementation are explicit here.
Empirical scores are evaluated directly from integer counts, not from the
production preparation, fitting, or scoring functions. All outputs are aggregate.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
A_VALUES = (0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
B_VALUES = (0.0, 0.25, 0.5, 1.0, 2.0)
SOURCE_SHA256 = "5c61f93a15134b795e6655ecab18f2d30e7d52daed33168ca0a60406303ee100"
PINS = {
    "HCTbc": (13962, "2213a534ce36abd298d2c587ac6f83f64ccaa2c51c4980d809245388d882fa09", "2fa829c71b13a0fa8790f6b89bbd735e8b45b294763b9cbaeab9ab7f06a5abea"),
    "SW6bc": (688439, "790dfdcb71b469afb82f219fe9a133e717dcf202477c445e9f28ac99d3c80585", "070a2a4fcd1ab01dc1751eb45fae0b1729bc158b0f016b329ddc683cb783448a"),
}


def read_authenticated_cache(path, line):
    import numpy as np
    with np.load(path, allow_pickle=False) as source:
        if str(source["source_sha256"].item()) != SOURCE_SHA256:
            raise ValueError("wrong workbook label")
        columns = [f"DT{r}_{p}" for r in range(1, 5) for p in ("P1", "P2")]
        if source["columns"].tolist() != columns:
            raise ValueError("wrong sample order")
        counts, barcodes = source["counts"].copy(), source["barcodes"].copy()
    n, chash, bhash = PINS[line]
    if counts.shape != (n, 8) or counts.dtype.kind not in "iu" or np.any(counts < 0):
        raise ValueError("invalid canonical count table")
    if barcodes.shape != (n,) or barcodes.dtype.kind != "U":
        raise ValueError("invalid barcode vector")
    actual = hashlib.sha256(counts.astype("<i8", copy=False).tobytes(order="C")).hexdigest()
    identities = hashlib.sha256(("\n".join(barcodes.tolist()) + "\n").encode("ascii")).hexdigest()
    if actual != chash or identities != bhash:
        raise ValueError("canonical count or barcode fingerprint mismatch")
    return counts, barcodes


def _unit(values):
    total = math.fsum(values)
    if total <= 0 or any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError("invalid probability weights")
    return [v / total for v in values]


def _select(candidates):
    best = max(row[2] for row in candidates)
    return min((row for row in candidates if best - row[2] <= 1e-12),
               key=lambda row: (row[1], abs(row[0] - 1), row[0]))


def reference_nested(early, late, target, delta=1e-4, donor_early=None):
    """Small-array standard-library exhaustive reference, for synthetic tests."""
    k = len(early[0])
    x = [[(1 - delta) * v + delta / k for v in _unit(row)] for row in early]
    source = early if donor_early is None else donor_early
    donor_x = [[(1 - delta) * v + delta / k for v in _unit(row)] for row in source]
    y = [_unit(row) for row in late]
    training = [r for r in range(4) if r != target]

    def forecast(r, donors, a, b):
        z = [math.fsum(donor_x[d][i] for d in donors) / len(donors) for i in range(k)]
        return _unit([x[r][i] ** a * z[i] ** b for i in range(k)])

    def score(r, q):
        return math.fsum(p * math.log(v) for p, v in zip(y[r], q))

    rows = [(a, b, math.fsum(score(r, forecast(r, [d for d in training if d != r], a, b))
                             for r in training) / 3)
            for a, b in itertools.product(A_VALUES, B_VALUES)]
    aa, _, base_inner = _select([row for row in rows if row[1] == 0])
    a, b, full_inner = _select(rows)
    q0, q1 = forecast(target, training, aa, 0), forecast(target, training, a, b)
    return {
        "base_a": aa, "full_a": a, "full_b": b,
        "base_inner_mean_score": base_inner, "full_inner_mean_score": full_inner,
        "base_prediction": q0, "full_prediction": q1,
        "base_score": score(target, q0), "full_score": score(target, q1),
        "gain": math.fsum(p * math.log(v / u) for p, u, v in zip(y[target], q0, q1)),
        "untuned_own_score": score(target, x[target]),
        "untuned_other_score": score(target, forecast(target, training, 0, 1)),
    }


def empirical_reference(counts, frequency_floor=0.0, delta=1e-4, mode="joint_rebuild"):
    """Independent support construction and count-weighted direct-power fit."""
    import numpy as np
    raw = np.asarray(counts)
    totals = raw.sum(axis=0, dtype=np.int64)
    if np.any(totals <= 0):
        return {"available": False}
    chosen = raw.copy()
    rare = raw / totals < frequency_floor
    if mode in ("joint_rebuild", "both"):
        chosen[rare] = 0
    elif mode == "early_only":
        chosen[:, ::2] = np.where(rare[:, ::2], 0, chosen[:, ::2])
    elif mode == "late_only":
        chosen[:, 1::2] = np.where(rare[:, 1::2], 0, chosen[:, 1::2])
    elif mode != "donor_feature_only":
        raise ValueError("unknown intervention")
    support_source = chosen if mode == "joint_rebuild" else raw
    keep = np.any(support_source[:, ::2] > 0, axis=1)
    grouped = np.concatenate([chosen[keep], chosen[~keep].sum(axis=0)[None, :]])
    if np.any(grouped.sum(axis=0) == 0):
        return {"available": False}
    early_counts, late_counts = grouped[:, ::2].T, grouped[:, 1::2].T
    early = early_counts / early_counts.sum(axis=1)[:, None]
    donor_p = early.copy()
    if mode == "donor_feature_only":
        donor_raw = np.where(rare[:, ::2], 0, raw[:, ::2])
        donor_grouped = np.concatenate([donor_raw[keep], donor_raw[~keep].sum(axis=0)[None, :]]).T
        if np.any(donor_grouped.sum(axis=1) == 0):
            return {"available": False}
        donor_p = donor_grouped / donor_grouped.sum(axis=1)[:, None]
    k = early.shape[1]
    x, donor_x = (1 - delta) * early + delta / k, (1 - delta) * donor_p + delta / k

    def forecast(r, donors, a, b):
        donor = sum((donor_x[d] for d in donors)) / len(donors)
        weights = np.power(x[r], a) * np.power(donor, b)
        return weights / weights.sum()

    def score(r, q):
        return float(np.dot(late_counts[r], np.log(q)) / late_counts[r].sum())

    folds = []
    for target in range(4):
        train = [r for r in range(4) if r != target]
        rows = []
        for a, b in itertools.product(A_VALUES, B_VALUES):
            scores = [score(r, forecast(r, [d for d in train if d != r], a, b)) for r in train]
            rows.append((a, b, math.fsum(scores) / 3))
        ba, _, inner0 = _select([row for row in rows if row[1] == 0])
        a, b, inner1 = _select(rows)
        q0, q1 = forecast(target, train, ba, 0), forecast(target, train, a, b)
        gain = float(np.dot(late_counts[target], np.log(q1 / q0)) / late_counts[target].sum())
        folds.append({
            "target_culture": f"DT{target+1}", "base_a": ba, "full_a": a, "full_b": b,
            "base_inner_mean_score": inner0, "full_inner_mean_score": inner1,
            "base_score": score(target, q0), "full_score": score(target, q1), "gain": gain,
            "untuned_own_score": score(target, x[target]),
            "untuned_other_score": score(target, forecast(target, train, 0, 1)),
        })
    return {"available": True, "early_dictionary_lineages": int(keep.sum()),
            "mean_gain": math.fsum(f["gain"] for f in folds) / 4, "folds": folds}


def reference_temporal(counts, barcodes, delta=1e-4):
    """Independent conditional selected-half profile-shape diagnostic."""
    import numpy as np
    early = counts[:, ::2].T
    late = counts[:, 1::2].T
    proportions = early / early.sum(axis=1)[:, None]
    homes = []
    for home in range(4):
        other = [r for r in range(4) if r != home]
        selected = np.flatnonzero((proportions[home] >= 1e-3) &
                                  np.all(proportions[other] < 1e-4, axis=0))
        parity = [hashlib.sha256(str(barcodes[i]).encode("ascii")).digest()[-1] % 2 for i in selected]
        halves = [[int(i) for i, half in zip(selected, parity) if half == h] for h in (0, 1)]
        record = {"home_culture": f"DT{home+1}", "selected_barcodes": len(selected),
                  "half_selected_barcodes": [len(v) for v in halves], "available": min(map(len, halves)) >= 10}
        if not record["available"]:
            homes.append(record)
            continue
        comparisons = []
        for receiver in other:
            for test_half in (0, 1):
                row = {"receiver_culture": f"DT{receiver+1}", "test_half": test_half}
                profiles = []
                for ids in halves:
                    e, l, y = early[home, ids], late[home, ids], early[receiver, ids]
                    if min(e.sum(), l.sum(), y.sum()) == 0:
                        profiles.append(None)
                    else:
                        profiles.append(((1-delta)*e/e.sum()+delta/len(ids),
                                         (1-delta)*l/l.sum()+delta/len(ids), y))
                row["available"] = all(p is not None for p in profiles)
                if row["available"]:
                    e, l, y = profiles[1-test_half]
                    weights = [i / 20 for i in range(21)]
                    values = [float(np.dot(y, np.log((1-w)*e+w*l)) / y.sum()) for w in weights]
                    maximum = max(values)
                    weight = next(w for w, v in zip(weights, values) if maximum-v <= 1e-12)
                    e, l, y = profiles[test_half]
                    gain = float(np.dot(y, np.log(((1-weight)*e+weight*l)/e)) / y.sum())
                    row.update({"selected_weight": weight, "gain": gain,
                                "receiver_test_read_count": int(y.sum()),
                                "receiver_test_occupied_barcodes": int(np.count_nonzero(y))})
                comparisons.append(row)
        available = [r["gain"] for r in comparisons if r["available"]]
        record["receiver_half_results"] = comparisons
        if len(available) == 6:
            record["mean_gain"] = math.fsum(available) / 6
        else:
            record["available"] = False
        homes.append(record)
    means = [r["mean_gain"] for r in homes if r.get("mean_gain") is not None]
    return {"home_results": homes, "available_home_cultures": len(means),
            "equal_home_mean_gain": math.fsum(means)/len(means) if means else None}


def compare_nested(reference, recorded, path="root", tolerance=1e-9):
    """Compare every reference field; production may include extra diagnostics."""
    if isinstance(reference, dict):
        return max([0.0] + [compare_nested(v, recorded[k], path+"."+k, tolerance) for k, v in reference.items()])
    if isinstance(reference, list):
        if len(reference) != len(recorded):
            raise AssertionError(f"{path}: changed list length")
        return max([0.0] + [compare_nested(v, w, f"{path}[{i}]", tolerance)
                           for i, (v, w) in enumerate(zip(reference, recorded))])
    if isinstance(reference, float):
        error = abs(reference-recorded)
        if not math.isfinite(error) or error > tolerance:
            raise AssertionError(f"{path}: {reference} vs {recorded}, error {error}")
        return error
    if reference != recorded:
        raise AssertionError(f"{path}: {reference!r} != {recorded!r}")
    return 0.0


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=ROOT/".compute/barcode-source")
    parser.add_argument("--reports-dir", type=Path, default=ROOT/".compute/barcode-lineage-memory")
    parser.add_argument("--out", type=Path, default=ROOT/".compute/barcode-lineage-memory/independent_verification.json")
    parser.add_argument("--lines", nargs="+", choices=list(PINS), default=list(PINS))
    parser.add_argument("--exploratory-hct", action="store_true")
    args = parser.parse_args()
    checks = []
    for line in args.lines:
        cache = args.cache_dir / f"{line}.npz"
        counts, barcodes = read_authenticated_cache(cache, line)
        path = args.reports_dir / f"{line}.json"
        report = json.loads(path.read_text())
        for floor in (0.0, 1e-4, 1e-3):
            reference = empirical_reference(counts, floor)
            target = report["primary"] if floor == 0 else next(v for v in report["abundance_sensitivities"] if v["frequency_floor"] == floor)
            error = compare_nested(reference, target)
            item = {"cell_line": line, "frequency_floor": floor, "matched": True,
                    "maximum_absolute_numeric_difference": error, "cache_sha256": digest(cache),
                    "report_sha256": digest(path), "reference": reference}
            checks.append(item)
            print(json.dumps({k:v for k,v in item.items() if k not in ("reference", "cache_sha256", "report_sha256")}), flush=True)
    if args.exploratory_hct:
        cache = args.cache_dir / "HCTbc.npz"
        counts, barcodes = read_authenticated_cache(cache, "HCTbc")
        path = args.reports_dir / "HCTbc_exploratory.json"
        report = json.loads(path.read_text())
        for mode in ("early_only", "late_only", "both", "donor_feature_only"):
            reference = empirical_reference(counts, 1e-5, mode=mode)
            target = next(v for v in report["ablations"] if v["mode"] == mode and v["floor"] == 1e-5)
            error = compare_nested(reference, target)
            checks.append({"cell_line": "HCTbc", "exploratory_A2": True, "mode": mode,
                           "floor": 1e-5, "matched": True, "maximum_absolute_numeric_difference": error,
                           "report_sha256": digest(path), "cache_sha256": digest(cache), "reference": reference})
        temporal = reference_temporal(counts, barcodes)
        error = compare_nested(temporal, report["temporal_diagnostic"])
        checks.append({"cell_line": "HCTbc", "exploratory_A2": True, "mode": "temporal",
                       "matched": True, "maximum_absolute_numeric_difference": error,
                       "report_sha256": digest(path), "cache_sha256": digest(cache), "reference": temporal})
        print(json.dumps({"cell_line":"HCTbc", "exploratory_checks":5, "matched":True}), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"method":"independent reconstructed count-weighted probability-space replay", "checks":checks}, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
