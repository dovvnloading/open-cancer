#!/usr/bin/env python3
"""Frozen-protocol, nonclinical persister-colony occupancy reanalysis.

The numerical observation model is a three-category multinomial. Biological
replicates, including all their early and late observations, are held out as
whole units. A Dirichlet-multinomial posterior predictive integrates sampling
parameter uncertainty. It does not estimate a genetic mutation rate or prove
biological independence of wells. No treatment advice; unreviewed research.

Only the Python standard library is used. This module makes no network calls.
"""

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import date
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROVENANCE = ("source", "source_file", "locator", "access_date")
THETA_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
DETECTION_GRID = (0.1, 0.5, 1.0)
KERNEL_SCALE_GRID = (0.25, 1.0, 4.0)
SCENARIOS = (
    ("primary", 0.5, None),
    ("cap_10", 0.5, 10.0),
    ("cap_100", 0.5, 100.0),
    ("cap_1000", 0.5, 1000.0),
    ("prior_0.25", 0.25, None),
    ("prior_1", 1.0, None),
)


def normalized_text_hash(path):
    """SHA-256 of UTF-8 text after CRLF and CR normalization to LF."""
    content = Path(path).read_bytes().decode("utf-8")
    content = content.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def _integer(value, name):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer, not boolean")
    try:
        number = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be an integer") from error
    if str(value).strip() != str(number) or number < 0:
        raise ValueError(f"{name} must be a nonnegative whole integer")
    return number


def _positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _rows(path, required):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{Path(path).name}: missing columns {sorted(missing)}")
        rows = []
        for line, row in enumerate(reader, 2):
            if None in row or any(row.get(name) is None for name in required):
                raise ValueError(f"{Path(path).name}:{line}: malformed CSV row")
            for name in PROVENANCE:
                if not row[name].strip():
                    raise ValueError(f"{Path(path).name}:{line}: empty {name}")
            date.fromisoformat(row["access_date"])
            rows.append(dict(row))
        return rows


def load_colonies(path):
    """Read one early/late occupancy vector per source biological replicate."""
    fields = (
        "cell_line", "clone", "biological_replicate", "total_wells", "early",
        "late", "early_day", "late_day_min", "late_day_max",
    ) + PROVENANCE
    rows = _rows(path, fields)
    seen = set()
    exposures = defaultdict(set)
    for row in rows:
        for name in ("cell_line", "clone", "biological_replicate"):
            if not row[name].strip():
                raise ValueError(f"colony row requires a nonempty {name}")
        key = (row["cell_line"], row["clone"], row["biological_replicate"])
        if key in seen:
            raise ValueError(f"duplicate biological replicate {key}")
        seen.add(key)
        exposures[key[:2]].add(row.get("exposure", ""))
        for name in ("total_wells", "early", "late"):
            row[name] = _integer(row[name], name)
        n = row["total_wells"]
        if n == 0 or row["early"] + row["late"] > n:
            raise ValueError(f"invalid mutually exclusive well counts for {key}")
        for name in ("early_day", "late_day_min", "late_day_max"):
            row[name] = _positive(row[name], name)
        if not row["early_day"] <= row["late_day_min"] <= row["late_day_max"]:
            raise ValueError(f"inconsistent readout times for {key}")
    if any(len(values) != 1 for values in exposures.values()):
        raise ValueError("different exposures must not be pooled within a clone")
    return sorted(rows, key=lambda row: (
        row["cell_line"], row["clone"], row["biological_replicate"]))


def load_populations(path):
    """Preserve source units; zero values are retained but excluded from log fits."""
    fields = (
        "cell_line", "clone", "biological_replicate", "day", "value",
        "measurement", "assay",
    ) + PROVENANCE
    rows = _rows(path, fields)
    for row in rows:
        for name in ("day", "value"):
            value = float(row[name])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"population {name} must be nonnegative and finite")
            row[name] = value
        for name in ("cell_line", "clone", "measurement", "assay"):
            if not row[name].strip():
                raise ValueError(f"population row requires a nonempty {name}")
    return sorted(rows, key=lambda row: (
        row["cell_line"], row["clone"], row["assay"],
        row["measurement"], row["biological_replicate"], row["day"],
        row["locator"]))


def load_published_parameters(path):
    """Published fitted/context values, never observations for a new raw-data fit."""
    fields = ("source_table", "cell_line", "clone", "parameter", "statistic",
              "value", "units", "raw_cell_value", "rich_text_value", "sheet",
              "cell", "licence", "note") + PROVENANCE
    rows = _rows(path, fields)
    seen = set()
    for row in rows:
        row["source_table"] = _integer(row["source_table"], "source_table")
        row["value"] = float(row["value"])
        if not math.isfinite(row["value"]) or row["value"] < 0:
            raise ValueError("published parameter must be finite and nonnegative")
        key = tuple(row[name] for name in (
            "source_table", "cell_line", "clone", "parameter", "statistic"))
        if key in seen:
            raise ValueError(f"duplicate published parameter {key}")
        seen.add(key)
    return sorted(rows, key=lambda row: (
        row["source_table"], row["cell_line"], row["clone"],
        row["parameter"], row["statistic"]))


def published_parameter_context(rows):
    """Summarize source estimates; do not pass them to the colony likelihood.

    The source TP-model D_p parameter is converted to g=-D_p only for the
    frozen fallback sensitivity. It is not a newly measured net growth rate.
    Source uncertainty labels and worksheet/caption discrepancies are retained.
    """
    lookup = {(row["source_table"], row["cell_line"], row["clone"],
               row["parameter"], row["statistic"]): row for row in rows}
    rate_comparisons = []
    for line, clone in sorted({(row["cell_line"], row["clone"]) for row in rows}):
        values = {}
        for population in ("sensitive", "persister"):
            parameter = f"{population}_chronological_resistance_rate"
            for statistic in ("estimate", "uncertainty_sd"):
                row = lookup.get((5, line, clone, parameter, statistic))
                if row is None:
                    continue
                suffix = "_per_cell_per_day"
                if not row["units"].endswith(suffix):
                    raise ValueError("unsupported published resistance-rate unit")
                # Unit multipliers are read from the committed source row,
                # including the distinct sensitive and persister exponents.
                scale = _positive(row["units"][:-len(suffix)], "rate unit multiplier")
                values[f"{population}_{statistic}_per_cell_per_day"] = row["value"] * scale
        required = {f"{population}_{statistic}_per_cell_per_day"
                    for population in ("sensitive", "persister")
                    for statistic in ("estimate", "uncertainty_sd")}
        if required <= values.keys():
            denominator = values["sensitive_estimate_per_cell_per_day"]
            if denominator <= 0:
                raise ValueError("positive sensitive estimate needed for a published ratio")
            rate_comparisons.append({
                "cell_line": line, "clone": clone, **values,
                "ratio_of_published_point_estimates":
                values["persister_estimate_per_cell_per_day"] / denominator,
                "interpretation": "published model-derived rates, not newly fitted rates",
            })
    anchors = []
    for row in rows:
        if (row["source_table"] != 4 or row["parameter"] != "D_p"
                or row["statistic"] not in (
                    "posterior_mean", "hdi_lower_source_label_3_percent",
                    "hdi_upper_source_label_97_percent")):
            continue
        loss = row["value"]
        anchors.append({
            "cell_line": row["cell_line"], "clone": row["clone"],
            "source_statistic": row["statistic"], "source_locator": row["locator"],
            "published_D_p_per_day": loss, "assumed_net_growth_per_day": -loss,
            "turnover_sensitivity": [{
                "birth_per_day": birth, "death_per_day": birth + loss,
                "net_per_day": -loss,
            } for birth in (0.0, 0.5, 1.0)],
        })
    return {
        "rows": rows, "rate_comparisons": rate_comparisons,
        "fallback_net_loss_anchors": anchors,
        "used_to_fit_colony_counts": False,
        "raw_population_recalibration": False,
        "source_notes": sorted({row["note"] for row in rows if row["note"]}),
    }


def counts_of(row):
    """Mutually exclusive categories; late risk set is n minus early positives."""
    return (row["early"], row["late"],
            row["total_wells"] - row["early"] - row["late"])


def posterior_alpha(training_rows, prior=0.5, cap=None):
    """Training-only posterior, optionally broadened at fixed category means."""
    prior = _positive(prior, "prior")
    alpha = [prior, prior, prior]
    for row in training_rows:
        for index, count in enumerate(counts_of(row)):
            if count < 0:
                raise ValueError("invalid count vector")
            alpha[index] += count
    total = math.fsum(alpha)
    if cap is not None:
        cap = _positive(cap, "concentration cap")
        if total > cap:
            alpha = [value * cap / total for value in alpha]
    return tuple(alpha)


def dirichlet_multinomial_logpmf(counts, alpha):
    """Integrated full-vector log probability, including its multinomial factor."""
    counts = tuple(counts)
    alpha = tuple(_positive(value, "Dirichlet parameter") for value in alpha)
    if len(counts) != len(alpha) or not counts:
        raise ValueError("count and parameter dimensions must agree")
    if any(isinstance(count, bool) or not isinstance(count, int) or count < 0
           for count in counts):
        raise ValueError("counts must be nonnegative integers")
    n = sum(counts)
    concentration = math.fsum(alpha)
    return math.fsum([
        math.lgamma(n + 1),
        -math.fsum(math.lgamma(count + 1) for count in counts),
        math.lgamma(concentration),
        -math.lgamma(concentration + n),
        math.fsum(math.lgamma(value + count) - math.lgamma(value)
                  for value, count in zip(alpha, counts)),
    ])


def beta_binomial_pmf(n, a, b):
    """Complete beta-binomial mass, evaluated in log space without simulation.

    Numerical renormalization removes roundoff from log-gamma evaluations; it
    does not truncate the support. n is the original number of assay wells.
    """
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("n must be a nonnegative integer")
    a, b = _positive(a, "a"), _positive(b, "b")
    constant = (math.lgamma(n + 1) + math.lgamma(a + b)
                - math.lgamma(n + a + b) - math.lgamma(a) - math.lgamma(b))
    logarithms = [constant - math.lgamma(k + 1) - math.lgamma(n - k + 1)
                  + math.lgamma(k + a) + math.lgamma(n - k + b)
                  for k in range(n + 1)]
    maximum = max(logarithms)
    weights = [math.exp(value - maximum) for value in logarithms]
    total = math.fsum(weights)
    return [weight / total for weight in weights]


def beta_binomial_interval(n, a, b, level=0.95):
    """Discrete equal-tail predictive interval, with integrated parameter error."""
    if not 0 < level < 1:
        raise ValueError("interval level must be between zero and one")
    probabilities = beta_binomial_pmf(n, a, b)
    thresholds = ((1.0 - level) / 2.0, (1.0 + level) / 2.0)
    result = []
    for threshold in thresholds:
        cumulative = 0.0
        quantile = n
        for k, probability in enumerate(probabilities):
            cumulative += probability
            if cumulative >= threshold:
                quantile = k
                break
        result.append(quantile)
    return result


def alias_probabilities(a, b_delayed, b_acquired):
    """Operational observable probabilities; origin labels enter only as a sum."""
    if any(not math.isfinite(value) or value < 0
           for value in (a, b_delayed, b_acquired)):
        raise ValueError("event intensities must be nonnegative and finite")
    b = b_delayed + b_acquired
    return (-math.expm1(-a), math.exp(-a) * -math.expm1(-b),
            math.exp(-a - b))


def positive_mark_alias(event_components):
    """Exact finite-mixture version of the protocol's vector-kernel identity.

    Each component is (nonnegative intensity, {integer tuple mark: mass}).
    Zero-vector marks are invisible. Inputs are converted to Fraction; use
    rational inputs for an exact arithmetic certificate. Returns (B, K_D).
    When B is zero the observed component is absent and K_D is unspecified,
    represented by an empty dictionary. This is algebra, not empirical data.
    """
    intensities = defaultdict(Fraction)
    dimensions = set()
    for intensity, kernel in event_components:
        intensity = Fraction(intensity)
        if intensity < 0 or not kernel:
            raise ValueError("nonnegative intensity and a nonempty kernel required")
        kernel = {mark: Fraction(mass) for mark, mass in kernel.items()}
        if sum(kernel.values(), Fraction()) != 1 or any(mass < 0 for mass in kernel.values()):
            raise ValueError("kernel masses must be nonnegative and sum to one")
        for mark, mass in kernel.items():
            if (not isinstance(mark, tuple) or not mark
                    or any(type(value) is not int or value < 0 for value in mark)):
                raise ValueError("marks must be nonempty nonnegative integer tuples")
            dimensions.add(len(mark))
            if any(mark) and mass:
                intensities[mark] += intensity * mass
    if len(dimensions) > 1:
        raise ValueError("all vector marks must have the same dimension")
    total = sum(intensities.values(), Fraction())
    if total == 0:
        return Fraction(), {}
    return total, {mark: value / total for mark, value in sorted(intensities.items())
                   if value > 0}


def _clone_summaries(groups):
    result = []
    for (line, clone), rows in sorted(groups.items()):
        n = sum(row["total_wells"] for row in rows)
        totals = [sum(counts_of(row)[index] for row in rows) for index in range(3)]
        alpha = posterior_alpha(rows)
        probabilities = [value / math.fsum(alpha) for value in alpha]
        p_early, _, p_negative = probabilities
        a = -math.log1p(-p_early)
        b = -math.log(p_negative / (1 - p_early))
        aliases = []
        for theta in THETA_GRID:
            delayed, acquired = (1 - theta) * b, theta * b
            observed = alias_probabilities(a, delayed, acquired)
            aliases.append({
                "theta": theta, "b_delayed": delayed, "b_acquired": acquired,
                "probabilities": list(observed),
                "max_probability_error": max(abs(x - y)
                                             for x, y in zip(observed, probabilities)),
            })
        # For well occupancy, a positive mark's magnitude is irrelevant: each
        # event is either detected or not. Varying a positive clone-size scale
        # therefore leaves this conversion unchanged. Multiplication by scale
        # would incorrectly substitute mean clone size for Pr(nonzero mark).
        # Detection here means event-level nonzero detection, not a measured
        # per-cell plating efficiency. These are conditional sensitivity values.
        required_pool = [{
            "detection_efficiency": detection, "kernel_scale": scale,
            "required_latent_intensity_per_well": b / detection,
        } for detection in DETECTION_GRID for scale in KERNEL_SCALE_GRID]
        late_risk = n - totals[0]
        result.append({
            "cell_line": line, "clone": clone,
            "biological_replicates": len(rows),
            "replicate_ids": [row["biological_replicate"] for row in rows],
            "total_wells": n, "early": totals[0], "late": totals[1],
            "negative": totals[2], "late_risk_wells": late_risk,
            "empirical_early_fraction": totals[0] / n,
            "empirical_late_fraction_of_initial": totals[1] / n,
            "empirical_late_fraction_of_risk": totals[1] / late_risk if late_risk else None,
            "early_by_replicate": [row["early"] for row in rows],
            "late_by_replicate": [row["late"] for row in rows],
            "full_data_posterior_alpha": list(alpha),
            "full_data_posterior_mean_probabilities": probabilities,
            "effective_early_intensity": a, "effective_late_intensity": b,
            "no_delay_no_acquisition_restriction_falsified": totals[1] > 0,
            "aliases": aliases, "latent_intensity_scalings": required_pool,
        })
    return result


def predictive_scenario(groups, label, prior=0.5, cap=None):
    folds = []
    unavailable = []
    for (line, clone), rows in sorted(groups.items()):
        if len(rows) < 2:
            unavailable.append({"cell_line": line, "clone": clone,
                                "reason": "fewer than two biological replicates"})
            continue
        for held_out in rows:
            held_id = held_out["biological_replicate"]
            training = [row for row in rows if row["biological_replicate"] != held_id]
            alpha = posterior_alpha(training, prior=prior, cap=cap)
            concentration = math.fsum(alpha)
            observed = counts_of(held_out)
            n = held_out["total_wells"]
            early_interval = beta_binomial_interval(n, alpha[0], alpha[1] + alpha[2])
            late_interval = beta_binomial_interval(n, alpha[1], alpha[0] + alpha[2])
            score = dirichlet_multinomial_logpmf(observed, alpha)
            folds.append({
                "cell_line": line, "clone": clone, "held_out_replicate": held_id,
                "training_replicates": [row["biological_replicate"] for row in training],
                "counts": list(observed), "total_wells": n, "alpha": list(alpha),
                "predictive_means": [n * value / concentration for value in alpha],
                "early_interval": early_interval, "late_interval": late_interval,
                "early_covered": early_interval[0] <= observed[0] <= early_interval[1],
                "late_covered": late_interval[0] <= observed[1] <= late_interval[1],
                # One observable law, expressed under two latent allocations.
                "log_score_m0": score, "log_score_m1": score,
            })
    total = math.fsum(fold["log_score_m0"] for fold in folds)
    by_clone = []
    for line, clone in sorted(groups):
        selected = [fold for fold in folds
                    if (fold["cell_line"], fold["clone"]) == (line, clone)]
        if selected:
            by_clone.append({
                "cell_line": line, "clone": clone, "n_folds": len(selected),
                "sum_log_score": math.fsum(fold["log_score_m0"] for fold in selected),
                "early_covered": sum(fold["early_covered"] for fold in selected),
                "late_covered": sum(fold["late_covered"] for fold in selected),
            })
    return {
        "label": label, "prior": prior, "cap": cap, "folds": folds,
        "unavailable": unavailable, "n_scored_folds": len(folds),
        "sum_log_score_m0": total, "sum_log_score_m1": total,
        "score_difference": 0.0,
        "early_covered": sum(fold["early_covered"] for fold in folds),
        "late_covered": sum(fold["late_covered"] for fold in folds),
        "by_clone": by_clone,
    }


def log_slope(points):
    """Least-squares slope of log population/proxy against day, with intercept."""
    points = list(points)
    if len(points) < 2 or len({time for time, _ in points}) < 2:
        raise ValueError("two distinct positive-measurement times are required")
    if any(not math.isfinite(time) or not math.isfinite(value) or value <= 0
           for time, value in points):
        raise ValueError("finite times and positive finite measurements required")
    mean_time = math.fsum(time for time, _ in points) / len(points)
    mean_log = math.fsum(math.log(value) for _, value in points) / len(points)
    numerator = math.fsum((time - mean_time) * (math.log(value) - mean_log)
                          for time, value in points)
    denominator = math.fsum((time - mean_time) ** 2 for time, _ in points)
    slope = numerator / denominator
    intercept = mean_log - slope * mean_time
    return slope, intercept


def calibrate_dynamics(rows):
    """Separate source-defined assay/replicate groups; no colony-target leakage."""
    groups = defaultdict(list)
    for row in rows:
        key = (row["cell_line"], row["clone"], row["assay"],
               row["measurement"], row["biological_replicate"])
        groups[key].append(row)
    fits, gaps = [], []
    for key, group in sorted(groups.items()):
        line, clone, assay, measurement, replicate = key
        for minimum_day in (14.0, 10.0):
            positive_by_day = defaultdict(list)
            excluded_nonpositive = 0
            for row in group:
                if row["day"] >= minimum_day:
                    if row["value"] > 0:
                        positive_by_day[row["day"]].append(row["value"])
                    else:
                        excluded_nonpositive += 1
            common = {
                "cell_line": line, "clone": clone, "assay": assay,
                "measurement": measurement, "biological_replicate": replicate,
                "minimum_day": minimum_day,
                "excluded_nonpositive": excluded_nonpositive,
            }
            if len(positive_by_day) < 2:
                gaps.append({**common, "reason": "fewer than two distinct eligible times"})
                continue
            # Multiple source measurements at a time are descriptive technical
            # replicates unless explicitly assigned different biological IDs.
            points = [(day, math.fsum(values) / len(values))
                      for day, values in sorted(positive_by_day.items())]
            g, intercept = log_slope(points)
            turnover = []
            for free_rate in (0.0, 0.5, 1.0):
                birth = free_rate if g <= 0 else free_rate + g
                death = free_rate - g if g <= 0 else free_rate
                turnover.append({"birth_per_day": birth, "death_per_day": death,
                                 "net_per_day": birth - death})
            fits.append({
                **common, "n_source_measurements": sum(map(len, positive_by_day.values())),
                "n_distinct_times": len(points), "days": [point[0] for point in points],
                "time_means": [point[1] for point in points],
                "net_log_slope_per_day": g, "log_intercept": intercept,
                "turnover_sensitivity": turnover,
            })
    if not rows:
        gaps.append({"reason": "no accessible raw population measurements supplied"})
    return {"fits": fits, "gaps": gaps,
            "interpretation": "net log-slope only; birth and death are not identified"}


def analyze(data_dir=None):
    data_dir = Path(data_dir) if data_dir is not None else HERE / "data"
    colonies = load_colonies(data_dir / "colony_counts.csv")
    if not colonies:
        raise ValueError("no colony observations; cannot run a fitted reanalysis")
    population_path = data_dir / "populations.csv"
    populations = load_populations(population_path) if population_path.exists() else []
    parameter_path = data_dir / "published_parameters.csv"
    parameters = load_published_parameters(parameter_path) if parameter_path.exists() else []
    source_path = data_dir / "sources.json"
    metadata = json.loads(source_path.read_text(encoding="utf-8")) if source_path.exists() else {}
    source_status = {
        "analysis_status": metadata.get("analysis_status", "source_status_unspecified"),
        "raw_archive_verified": metadata.get("raw_archive_verified", False),
        "raw_population_calibration_completed":
        metadata.get("raw_population_calibration_completed", False),
        "colony_source_types": sorted({row.get("source_type", "unspecified") for row in colonies}),
    }
    groups = defaultdict(list)
    for row in colonies:
        groups[(row["cell_line"], row["clone"])].append(row)
    hashes = {Path("data/colony_counts.csv").as_posix():
              normalized_text_hash(data_dir / "colony_counts.csv")}
    if population_path.exists():
        hashes[Path("data/populations.csv").as_posix()] = normalized_text_hash(population_path)
    if parameter_path.exists():
        hashes[Path("data/published_parameters.csv").as_posix()] = normalized_text_hash(parameter_path)
    if source_path.exists():
        hashes[Path("data/sources.json").as_posix()] = normalized_text_hash(source_path)
    protocol = HERE / "PROTOCOL.md"
    if protocol.exists():
        hashes[Path("PROTOCOL.md").as_posix()] = normalized_text_hash(protocol)
    return {
        "status": "nonclinical; unreviewed; no treatment advice",
        "estimand": "fraction of effective detectable late event intensity acquired after exposure",
        "primary_unit": "whole source biological replicate",
        "colony_observation_rows": len(colonies),
        "population_observation_rows": len(populations),
        "source_status": source_status,
        "verdict": ("not identifiable from published count target; requested raw audit incomplete"
                    if not source_status["raw_archive_verified"]
                    else "not identifiable within the allowed observation-kernel family"),
        "input_hashes": hashes,
        "settings": {
            "primary_prior": 0.5, "predictive_interval_level": 0.95,
            "theta_grid": list(THETA_GRID), "detection_efficiencies": list(DETECTION_GRID),
            "kernel_scales": list(KERNEL_SCALE_GRID),
            "concentration_caps": [10, 100, 1000],
            "prior_sensitivities": [0.25, 1.0],
        },
        "clone_summaries": _clone_summaries(groups),
        "predictive_scenarios": [predictive_scenario(groups, *scenario)
                                 for scenario in SCENARIOS],
        "dynamics": calibrate_dynamics(populations),
        "published_context": published_parameter_context(parameters),
    }


def _fmt(value, digits=6):
    return "unavailable" if value is None else f"{value:.{digits}f}"


def _interval(values):
    return f"[{values[0]}, {values[1]}]"


def render_markdown(result):
    """Render every displayed computed statistic directly from analyze()."""
    out = [
        "<!-- BEGIN COMPUTED REPORT -->",
        "## Recomputed quantitative report",
        "",
        "**Nonclinical, unreviewed; no treatment advice.** These tables are generated by",
        "`analysis.py`. Count probabilities concern well occupancy, not counted genetic",
        "acquisitions.",
        "",
        f"**Verdict: {result['verdict']}.**",
        "",
        f"Source status: `{result['source_status']['analysis_status']}`. "
        f"Raw archive checksum verified: **{str(result['source_status']['raw_archive_verified']).lower()}**. "
        "Raw population calibration completed: "
        f"**{str(result['source_status']['raw_population_calibration_completed']).lower()}**.",
        "",
        "Target source labels: " + ", ".join(
            f"`{label}`" for label in result["source_status"]["colony_source_types"]) + ".",
        "",
        "### Colony observations",
        "",
        "| Cell line / clone | Biological replicates | Wells | Early E | Newly late L | Neither U | Late risk n−E | E/n | L/(n−E) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in result["clone_summaries"]:
        out.append(f"| {row['cell_line']} / {row['clone']} | {row['biological_replicates']} | "
                   f"{row['total_wells']} | {row['early']} | {row['late']} | {row['negative']} | "
                   f"{row['late_risk_wells']} | {_fmt(row['empirical_early_fraction'])} | "
                   f"{_fmt(row['empirical_late_fraction_of_risk'])} |")
    out.extend([
        "", "### Held-out predictive scores", "",
        "Each row holds out a complete biological replicate. E and L intervals are",
        "marginal equal-tail 95% beta-binomial predictive intervals, not simultaneous",
        "intervals or confidence intervals for acquisition. Every joint log score includes",
        "the multinomial coefficient. M0 and M1 share this observable distribution, so",
        "their scores agree algebraically. Caps are prespecified dispersion sensitivities.",
    ])
    for scenario in result["predictive_scenarios"]:
        out.extend([
            "", f"#### {scenario['label']}", "",
            "| Clone | Held out | n; E,L,U | Predicted E | E interval | E covered | Predicted L | L interval | L covered | Joint log score, M0 = M1 |",
            "|---|---|---|---:|---|---|---:|---|---|---:|",
        ])
        for fold in scenario["folds"]:
            observed = ",".join(map(str, fold["counts"]))
            out.append(f"| {fold['cell_line']}/{fold['clone']} | {fold['held_out_replicate']} | "
                       f"{fold['total_wells']}; {observed} | {_fmt(fold['predictive_means'][0], 3)} | "
                       f"{_interval(fold['early_interval'])} | {'yes' if fold['early_covered'] else 'no'} | "
                       f"{_fmt(fold['predictive_means'][1], 3)} | {_interval(fold['late_interval'])} | "
                       f"{'yes' if fold['late_covered'] else 'no'} | {_fmt(fold['log_score_m0'])} |")
        n = scenario["n_scored_folds"]
        out.extend([
            "", f"Joint log-score sum: **{_fmt(scenario['sum_log_score_m0'])}** for each model; "
            f"difference **{_fmt(scenario['score_difference'])}**. "
            f"Descriptive marginal coverage: E **{scenario['early_covered']}/{n}**, "
            f"L **{scenario['late_covered']}/{n}**.",
        ])
        for gap in scenario["unavailable"]:
            out.extend(["", f"{gap['cell_line']}/{gap['clone']}: predictive validation unavailable "
                        f"({gap['reason']})."])
    out.extend([
        "", "### Explicit origin aliases", "",
        "The probabilities below are full-data posterior means under the primary prior.",
        "The displayed intensities are deterministic transforms of those means, not",
        "posterior mean mutation rates. All listed allocations leave the observable law",
        "unchanged; numerical equality is checked separately from the analytic identity.",
        "",
        "| Clone | Effective a | Effective b | p(E) | p(L) | p(U) |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for row in result["clone_summaries"]:
        p = row["full_data_posterior_mean_probabilities"]
        out.append(f"| {row['cell_line']}/{row['clone']} | {_fmt(row['effective_early_intensity'])} | "
                   f"{_fmt(row['effective_late_intensity'])} | "
                   f"{_fmt(p[0])} | {_fmt(p[1])} | {_fmt(p[2])} |")
    out.extend(["", "| Clone | theta | Delayed pre-existing b_D | Acquired b_I |",
                "|---|---:|---:|---:|"])
    for row in result["clone_summaries"]:
        for alias in row["aliases"]:
            out.append(f"| {row['cell_line']}/{row['clone']} | {_fmt(alias['theta'], 2)} | "
                       f"{_fmt(alias['b_delayed'])} | {_fmt(alias['b_acquired'])} |")
    out.extend([
        "", "The following conversion assumes the stated probability that a latent",
        "event leaves a nonzero detected mark. It is an event-level probability,",
        "not a measured per-cell plating efficiency. Changing the positive kernel's",
        "relative size scale leaves well occupancy unchanged at fixed nonzero-mark",
        "probability. These are sensitivity assumptions, not measured founder counts.",
        "",
        "| Clone | Nonzero-event detection probability | Positive kernel relative scale | Latent late intensity per well |",
        "|---|---:|---:|---:|",
    ])
    for row in result["clone_summaries"]:
        for conversion in row["latent_intensity_scalings"]:
            out.append(f"| {row['cell_line']}/{row['clone']} | "
                       f"{_fmt(conversion['detection_efficiency'], 2)} | "
                       f"{_fmt(conversion['kernel_scale'], 2)} | "
                       f"{_fmt(conversion['required_latent_intensity_per_well'])} |")
    out.extend(["", "### Separate population calibration", ""])
    if result["dynamics"]["fits"]:
        out.extend([
            "These are separate population/proxy measurements. No colony observation",
            "enters the log-slope calibration. Replicate labels are exactly source labels.",
            "",
            "| Clone | Assay / measurement | Replicate | Minimum day | Distinct times | Net log-slope per day |",
            "|---|---|---|---:|---:|---:|",
        ])
        for row in result["dynamics"]["fits"]:
            replicate = row["biological_replicate"] or "source aggregate; unspecified"
            out.append(f"| {row['cell_line']}/{row['clone']} | {row['assay']} / {row['measurement']} | "
                       f"{replicate} | {_fmt(row['minimum_day'], 0)} | {row['n_distinct_times']} | "
                       f"{_fmt(row['net_log_slope_per_day'])} |")
        out.extend(["", "| Clone / assay / replicate / window | Birth per day | Death per day | Net per day |",
                    "|---|---:|---:|---:|"])
        for row in result["dynamics"]["fits"]:
            for rates in row["turnover_sensitivity"]:
                out.append(f"| {row['cell_line']}/{row['clone']} / {row['assay']} / "
                           f"{row['biological_replicate'] or 'unspecified'} / "
                           f"≥{_fmt(row['minimum_day'], 0)} days | "
                           f"{_fmt(rates['birth_per_day'])} | {_fmt(rates['death_per_day'])} | "
                           f"{_fmt(rates['net_per_day'])} |")
    else:
        out.append("No raw-data late-phase population slope could be calibrated.")
    if result["dynamics"]["gaps"]:
        out.extend(["", "Calibration gaps:", ""])
        for gap in result["dynamics"]["gaps"]:
            if "cell_line" in gap:
                out.append(f"- {gap['cell_line']}/{gap['clone']}, {gap['assay']}, "
                           f"replicate {gap['biological_replicate'] or 'unspecified'}, "
                           f"day ≥{_fmt(gap['minimum_day'], 0)}: {gap['reason']}.")
            else:
                out.append(f"- {gap['reason']}.")
    context = result["published_context"]
    out.extend(["", "### Published parameter context and fallback sensitivity", ""])
    if not context["rows"]:
        out.append("No published parameter context was supplied.")
    else:
        out.extend([
            "All source estimates below come from the committed supplementary-table",
            "extraction. They are **published model-dependent context**, never training",
            "observations for the colony fit or a substitute for raw population calibration.",
            "Untreated birth/death means and standard errors, TP-model prior limits,",
            "persister-loss posterior summaries, and resistance-rate estimates retain",
            "their distinct source roles and uncertainty labels.",
            "",
            "| Supplementary table | Clone | Parameter | Source statistic | Value | Source units | Exact worksheet cell |",
            "|---:|---|---|---|---:|---|---|",
        ])
        for row in context["rows"]:
            out.append(f"| {row['source_table']} | {row['cell_line']}/{row['clone']} | "
                       f"`{row['parameter']}` | `{row['statistic']}` | {row['value']:.12g} | "
                       f"`{row['units']}` | `{row['locator']}` |")
        if context["source_notes"]:
            out.extend(["", "Source interpretation notes:", ""])
            out.extend(f"- {note}" for note in context["source_notes"])
        if context["rate_comparisons"]:
            out.extend([
                "", "#### Comparison with the authors' inferred rates", "",
                "Unit multipliers are taken from each source row before calculating the",
                "ratio. These are the authors' point estimates and reported uncertainties;",
                "the ratio is descriptive and has no newly inferred uncertainty interval.",
                "",
                "| Clone | Sensitive estimate per cell per day | Reported uncertainty | Persister estimate per cell per day | Reported uncertainty | Ratio of point estimates |",
                "|---|---:|---:|---:|---:|---:|",
            ])
            for row in context["rate_comparisons"]:
                out.append(f"| {row['cell_line']}/{row['clone']} | "
                           f"{row['sensitive_estimate_per_cell_per_day']:.8g} | "
                           f"{row['sensitive_uncertainty_sd_per_cell_per_day']:.8g} | "
                           f"{row['persister_estimate_per_cell_per_day']:.8g} | "
                           f"{row['persister_uncertainty_sd_per_cell_per_day']:.8g} | "
                           f"{_fmt(row['ratio_of_published_point_estimates'])} |")
        if context["fallback_net_loss_anchors"]:
            out.extend([
                "", "#### Published net-loss anchors; no raw-data recalibration", "",
                "The source TP-model persister-loss parameter D_p is used only as a",
                "conditional net-loss anchor, g = −D_p. For every published mean and",
                "each reported bound, the frozen birth-rate sensitivity is applied and",
                "death = birth + D_p. These are assumptions consistent with the same net",
                "loss; they do not identify actual persister division or death rates.",
                "The bounds retain their source labels because the worksheet and caption",
                "disagree about interval endpoints. These are not newly estimated intervals.",
                "",
                "| Clone | Source D_p statistic | D_p per day | Assumed birth per day | Required death per day | Assumed net per day |",
                "|---|---|---:|---:|---:|---:|",
            ])
            for row in context["fallback_net_loss_anchors"]:
                for rates in row["turnover_sensitivity"]:
                    out.append(f"| {row['cell_line']}/{row['clone']} | "
                               f"`{row['source_statistic']}` | {_fmt(row['published_D_p_per_day'])} | "
                               f"{_fmt(rates['birth_per_day'])} | {_fmt(rates['death_per_day'])} | "
                               f"{_fmt(rates['net_per_day'])} |")
    out.extend(["", "<!-- END COMPUTED REPORT -->"])
    return "\n".join(out) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=HERE / "data")
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    args = parser.parse_args(argv)
    result = analyze(args.data_dir)
    if args.format == "markdown":
        print(render_markdown(result), end="")
    else:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
