#!/usr/bin/env python3
"""Reproduce aggregate-rate compatibility certificates; standard library only.

Nonclinical computational research; no medical advice. This does not estimate
a biological k, supply missing SEER denominators, or conduct a significance
test. Every grid comparison uses certified rational Taylor-tail enclosures.
"""

from __future__ import annotations

import argparse
import csv
from decimal import localcontext
from fractions import Fraction as F
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

from verifier import (calibration_ratio, calibration_ratio_interval, decimal,
                      exp_interval, flexible_shape_calibration,
                      one_growth_bin_interval, power_bin)


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "seer_age_specific_2013_2017.csv"
KS = tuple(range(2, 11))
# A deliberately declared finite grid, NOT continuous optimization or a prior.
GROWTH_GRID = tuple(F(i, 200) for i in range(-12, 13))
SERIES_TERMS = 70


def source_sha256(path: Path = DATA) -> str:
    """SHA-256 of UTF-8 CSV text after universal-newline normalization to LF.

    Git may checkout this public CSV with CRLF under core.autocrlf=true.
    Newline-normalized text hashes identically on Linux, macOS, and Windows.
    This is a content-integrity checksum, not upstream data authentication.
    """
    with path.open("r", encoding="utf-8", newline=None) as stream:
        normalized_text = stream.read()
    return sha256(normalized_text.encode("utf-8")).hexdigest()


def _display(x: F, places: int = 6, sqrt: bool = False) -> str:
    with localcontext() as ctx:
        ctx.prec = 60
        y = decimal(x)
        if sqrt:
            y = y.sqrt()
        return f"{y:.{places}f}"


def _display_lower_bound(x: F, places: int = 8) -> str:
    """Downward rounding by integer division, preserving a certified lower bound."""
    units = x.numerator * 10 ** places // x.denominator
    return f"{units // 10 ** places}.{units % 10 ** places:0{places}d}"


def load_data(path: Path = DATA) -> dict[str, list[dict[str, str]]]:
    groups: dict[str, list[dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    keys = set()
    for row in rows:
        key = row["cancer_site"], row["age_group"]
        if key in keys:
            raise ValueError(f"Duplicate site/age key: {key}")
        keys.add(key)
        if row["case_count"] or row["person_years"]:
            raise ValueError("This analysis expects rates with missing denominators")
        if F(row["rate_per_100000"]) <= 0:
            raise ValueError("Positive unsuppressed rates required")
        if int(row["age_end_exclusive"]) - int(row["age_start"]) != 5:
            raise ValueError("Only five-year bins are used")
        start, end = int(row["age_start"]), int(row["age_end_exclusive"])
        if start < 0 or row["age_group"] != f"{start}-{end - 1}":
            raise ValueError("Age label must match its nonnegative half-open endpoints")
        for name in ("source_url", "source_table", "access_date", "reuse_terms_url"):
            if not row[name]:
                raise ValueError(f"Missing provenance: {name}")
        groups.setdefault(row["cancer_site"], []).append(row)
    for group in groups.values():
        group.sort(key=lambda r: int(r["age_start"]))
        if any(int(a["age_end_exclusive"]) > int(b["age_start"])
               for a, b in zip(group, group[1:])):
            raise ValueError("Age bins must be ordered and nonoverlapping")
        for name in ("behavior", "sex", "race", "registry_set", "period_start",
                     "period_end", "source_table", "source_column", "source_section",
                     "reported_decimal_places"):
            if len({row[name] for row in group}) != 1:
                raise ValueError(f"Cannot pool inconsistent site metadata: {name}")
    return groups


def _rounding_intervals(rows: list[dict[str, str]]) -> tuple[tuple[F, ...], tuple[F, ...]]:
    # Closed +/- half last displayed unit is a sensitivity convention, not a CI.
    values = tuple(F(r["rate_per_100000"]) for r in rows)
    half = tuple(F(1, 2 * 10 ** int(r["reported_decimal_places"])) for r in rows)
    return tuple(y - d for y, d in zip(values, half)), tuple(y + d for y, d in zip(values, half))


def _monotone_certificate(rows: list[dict[str, str]], lower: tuple[F, ...],
                          upper: tuple[F, ...]) -> dict:
    violations = []
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            if lower[i] > upper[j]:
                violations.append({"earlier_age": rows[i]["age_group"],
                                   "later_age": rows[j]["age_group"],
                                   "minimum_decline_per_100000": _display(lower[i] - upper[j], 1)})
    # Exact optimal multiplier for ANY nondecreasing vector, by interval
    # isotonic feasibility: lower_i/rho <= rho*upper_j for every i<=j.
    ratio = max([F(1)] + [lower[i] / upper[j]
                         for i in range(len(rows)) for j in range(i, len(rows))])
    return {"rounding_robust_decreases": violations,
            "homogeneous_monotone_mean_compatible": not violations,
            "minimum_rho_for_any_nondecreasing_mean": _display(ratio, sqrt=True),
            "interpretation": "deterministic rounded-mean compatibility; not a sampling test"}


def analyze() -> dict:
    data = load_data()
    unique_bins = sorted({(F(r["age_start"]), F(r["age_end_exclusive"]))
                          for rows in data.values() for r in rows})
    # Shared across sites; rational enclosures make each grid comparison exact.
    shapes = {(k, g): {age: one_growth_bin_interval(k, g, *age, SERIES_TERMS)
                       for age in unique_bins} for k in KS for g in GROWTH_GRID}
    results: dict[str, dict] = {}
    for site, rows in sorted(data.items()):
        bins = tuple((F(r["age_start"]), F(r["age_end_exclusive"])) for r in rows)
        lower, upper = _rounding_intervals(rows)
        per_k = []
        exact_power_ranks = []
        exact_grid_ranks = []
        for k in KS:
            powers = tuple(power_bin(k, *age) for age in bins)
            pure_ratio = calibration_ratio(lower, upper, powers)
            endpoint_ranges = tuple((a ** (k - 1) / factorial(k - 1),
                                     b ** (k - 1) / factorial(k - 1)) for a, b in bins)
            flex = flexible_shape_calibration(lower, upper, endpoint_ranges)
            candidates = []
            for g in GROWTH_GRID:
                enclosure = tuple(shapes[k, g][age] for age in bins)
                lo, hi = calibration_ratio_interval(lower, upper, enclosure)
                candidates.append((hi, lo, g))
            best_hi, best_lo, best_g = min(candidates)
            # The selected grid point is rigorously smallest if its upper bound
            # is below every competitor's lower bound. Do not silently call a
            # near tie unique; a bounded optimum enclosure still exists.
            contenders = [g for hi, lo, g in candidates if lo <= best_hi]
            optimum_lo = min(lo for hi, lo, g in candidates)
            assert best_hi - optimum_lo < F(1, 10**24)
            exact_power_ranks.append((pure_ratio, k))
            exact_grid_ranks.append((best_hi, optimum_lo, k))
            per_k.append({
                "k": k,
                "pure_power_min_rho": _display(pure_ratio, sqrt=True),
                "pure_power_unknown_within_bin_weights_rho_infimum": _display(flex, sqrt=True),
                "grid_min_rho": _display(best_hi, sqrt=True),
                "grid_growth_per_year": _display(best_g, 3),
                "grid_minimum_certified_unique": len(contenders) == 1,
                "grid_rho_squared_enclosure_width_less_than": "1e-24",
            })
        best_power_k = min(exact_power_ranks)[1]
        global_hi, global_lo, best_grid_k = min(exact_grid_ranks)
        global_contenders = [k for hi, lo, k in exact_grid_ranks if lo <= global_hi]
        best_power = next(r for r in per_k if r["k"] == best_power_k)
        best_grid = next(r for r in per_k if r["k"] == best_grid_k)
        results[site] = {
            "table": rows[0]["source_table"],
            "rows": len(rows),
            "sum_rates_exact": str(sum((F(r["rate_per_100000"]) for r in rows), F(0))),
            "identified_biological_k": None,
            "unrestricted_bin_multiplier_compatible_k": "every positive integer",
            "reason": "Theorem 7: positive bin shapes, free amplitude and unbounded age multipliers",
            "monotonicity": _monotone_certificate(rows, lower, upper),
            "best_pure_power_on_declared_k_grid": best_power["k"],
            "best_pure_power_rho": best_power["pure_power_min_rho"],
            "best_one_growth_on_declared_grid": best_grid["k"],
            "best_one_growth_rho": best_grid["grid_min_rho"],
            "best_one_growth_per_year": best_grid["grid_growth_per_year"],
            "joint_grid_k_minimum_certified_unique": len(global_contenders) == 1,
            "by_k": per_k,
        }
    # Explicitly synthetic exposure certificate valid for ANY downstream clone
    # having h(T)<=1/1000 per year, not estimated from the SEER rates.
    synthetic = []
    for rate in (F(100), F(1000), F(10000)):
        W, H = F(100000), F(1, 1000)
        d_upper = W * H / rate
        elo, ehi = exp_interval(-d_upper, 100)
        synthetic.append({"extra_conversion_rate_per_year": str(rate),
                          "maximum_exposure_density_persons": str(W),
                          "assumed_hazard_upper_per_year": str(H),
                          "coupling_error_lower_bound": _display_lower_bound(elo / 2, 8),
                          "rational_enclosure_width_less_than": "1e-100"})
        assert (ehi - elo) / 2 < F(1, 10**100)
    return {
        "label": "nonclinical computational research; no medical advice; unreviewed",
        "source_sha256": source_sha256(DATA),
        "source_access_date": "2026-11-17",
        "rate_period": "2013-2017",
        "observation_assumptions": {
            "rounding_sensitivity": "closed +/-0.05 per 100000, not a confidence interval",
            "default_within_bin_age_weights": "uniform; source weights are unavailable",
            "unknown_within_bin_weights": "power-family endpoint closure gives a sharp infimum, not guaranteed attainment",
            "multiplier": "a_j in [1/rho,rho] and a free common amplitude",
            "rho_interpretation": "relative age-shape compatibility threshold, not measured reporting bias",
            "limiting_family": "rare-terminal-marking/high-immigration limit; one nonzero net-growth rate",
            "growth_grid": "-0.060 to 0.060 per year inclusive, spacing 0.005",
            "k_grid": list(KS),
            "optimization": "exact finite-grid minimum only; no continuous or full MSCE optimum claimed",
            "person_years_and_counts": "unavailable; no data-specific sampling probabilities",
        },
        "cancers": results,
        "synthetic_poisson_exposure_certificates": synthetic,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional JSON path; default is stdout")
    args = parser.parse_args()
    content = json.dumps(analyze(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(content, encoding="utf-8")
    else:
        print(content, end="")


if __name__ == "__main__":
    main()
