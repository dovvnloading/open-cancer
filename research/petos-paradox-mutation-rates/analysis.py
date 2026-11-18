"""Frozen aggregate reanalysis. Nonclinical research; no medical advice.

Only Python's standard library is used. Input tables contain published species
aggregates, never individual records. Mortality-to-onset and mean-lifetime
bridges are deliberately labelled proxies. See PROTOCOL.md before interpreting
a fit, interval box, or exact pair certificate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import statistics
import sys


HERE = Path(__file__).resolve().parent
FREEZE = "b98d28e5258c09132d197e26529ff4f58d18b8a9"
K_VALUES = tuple(range(1, 9))
Z95 = 1.959963984540054
LN10 = math.log(10)
_spec = importlib.util.spec_from_file_location("petos_theorem_verifier", HERE / "verifier.py")
V = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = V
_spec.loader.exec_module(V)


def normalized_bytes(path):
    """Content hashes are invariant to LF, CRLF and lone-CR checkouts."""
    return Path(path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def input_hashes(base=HERE):
    files = sorted((base / "data").glob("*.csv")) + [base / "PROTOCOL.md"]
    return {path.relative_to(base).as_posix(): hashlib.sha256(normalized_bytes(path)).hexdigest()
            for path in files}


def read_csv(name, base=HERE):
    with (base / "data" / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def wilson(events, total, z=Z95):
    """Derived binomial score interval, NOT an interval printed by Vincze."""
    if not isinstance(events, int) or not isinstance(total, int) or not 0 <= events <= total or total <= 0:
        raise ValueError("integer binomial counts required")
    p = events / total
    z2 = z * z
    denom = 1 + z2 / total
    center = (p + z2 / (2 * total)) / denom
    half = z * math.sqrt(p * (1 - p) / total + z2 / (4 * total * total)) / denom
    return max(0.0, center - half), min(1.0, center + half)


def hazard(p):
    if not 0 <= p <= 1:
        raise ValueError("probability outside [0,1]")
    return math.inf if p == 1 else -math.log1p(-p)


def load_data(base=HERE):
    mutation = {r["species"]: r for r in read_csv("cagan_mutation_rates.csv", base)}
    traits = {r["species"]: r for r in read_csv("compton_life_history.csv", base)}
    mortality = read_csv("vincze_mortality.csv", base)
    rows = []
    for source in mortality:
        species = source["species"]
        if species not in mutation or species not in traits:
            continue
        mr, tr = mutation[species], traits[species]
        events, total = int(source["cancer_deaths"]), int(source["known_cause_deaths"])
        risk = events / total
        if not math.isclose(risk, float(source["cmr"]), abs_tol=2e-14, rel_tol=0):
            raise ValueError("published CMR does not match its stated denominator")
        if not 0 < risk < 1:
            raise ValueError("frozen log-hazard fit requires strictly interior CMR")
        female, male = float(source["female_mass_kg"]), float(source["male_mass_kg"])
        lower, upper = wilson(events, total)
        rows.append({
            "species": species, "order": source["order"],
            "m": float(mr["mean_sbs_per_genome_year"]),
            "slope": float(mr["zero_intercept_slope_sbs_per_genome_year"]),
            "slope_ci": ([float(mr["slope_ci95_lower"]), float(mr["slope_ci95_upper"])]
                         if mr["slope_ci95_lower"] else None),
            "M": (female + male) / 2,
            "M_box": [min(female, male), max(female, male)],
            "L": float(source["adult_life_expectancy_days"]) / 365.25,
            "L_max": float(tr["max_longevity_months"]) / 12,
            "events": events, "n": total,
            "total_animals": int(source["total_animals"]),
            "total_deaths": int(source["total_deaths"]),
            "R": risk, "R_box": [lower, upper], "H": hazard(risk),
            "H_box": [hazard(lower), hazard(upper)],
            "ICM": float(source["icm"]) if source["icm"] else None,
        })
    return sorted(rows, key=lambda r: r["species"]), mutation, traits


def panel_version(rows, name):
    answer = []
    for row in rows:
        r = row.copy()
        if name == "maximum_longevity":
            r["L"] = r["L_max"]
        elif name == "slope_estimator":
            r["m"] = r["slope"]
        elif name != "primary":
            raise ValueError("unknown frozen panel")
        r["B"] = r["m"] * r["L"]
        answer.append(r)
    return answer


def ols(x, y):
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("at least two paired observations required")
    xm, ym = statistics.fmean(x), statistics.fmean(y)
    ss = math.fsum((v - xm) ** 2 for v in x)
    if ss <= 1e-24:
        raise ValueError("rank-deficient design")
    b = math.fsum((a - xm) * (v - ym) for a, v in zip(x, y)) / ss
    return ym - b * xm, b


def fit(rows, family, k=0):
    if not rows:
        raise ValueError("empty training set")
    if family == "intercept":
        b, exponent = 0.0, 0
    elif family == "mass":
        b, exponent = 1.0, 0
    elif family in ("burden", "cell", "free_mass") and k in K_VALUES:
        b, exponent = (1.0 if family == "cell" else 0.0), k
    else:
        raise ValueError("unknown model")
    y = [math.log(r["H"]) - exponent * math.log(r["B"]) for r in rows]
    if family == "free_mass":
        a, b = ols([math.log(r["M"]) for r in rows], y)
    else:
        a = statistics.fmean(v - b * math.log(r["M"]) for v, r in zip(y, rows))
    return {"intercept": a, "mass_exponent": b, "burden_exponent": exponent}


def predict(parameters, row):
    return (parameters["intercept"] + parameters["mass_exponent"] * math.log(row["M"])
            + parameters["burden_exponent"] * math.log(row["B"]))


def bernoulli_nll(log_h, risk):
    """Per-known-death log loss; no probability clipping or made-up counts."""
    if math.isnan(log_h) or not 0 <= risk <= 1:
        raise ValueError("invalid log hazard or observed proportion")
    if log_h == -math.inf:
        return 0.0 if risk == 0 else math.inf
    if log_h > math.log(sys.float_info.max):
        return 0.0 if risk == 1 else math.inf
    h = math.exp(log_h)
    if log_h < -36:
        # log(1-exp(-H))=log(H)-H/2+O(H^2); the correction is
        # below log-scale floating precision here, even if H underflows.
        log_risk = log_h
    else:
        log_risk = math.log(-math.expm1(-h)) if h < 0.6931471805599453 else math.log1p(-math.exp(-h))
    return -risk * log_risk + (1 - risk) * h


def risk_from_log_h(log_h):
    if math.isnan(log_h):
        raise ValueError("NaN log hazard")
    return 1.0 if log_h > math.log(746.0) else -math.expm1(-math.exp(log_h))


def cross_validate(rows, family, k=0, grouping="species"):
    predictions = []
    for group in sorted({r[grouping] for r in rows}):
        training = [r for r in rows if r[grouping] != group]
        parameters = fit(training, family, k)
        for r in rows:
            if r[grouping] != group:
                continue
            log_h = predict(parameters, r)
            residual = (math.log(r["H"]) - log_h) / LN10
            predictions.append({
                "species": r["species"], "held_out_group": group,
                "predicted_log10_h": log_h / LN10,
                "predicted_risk": risk_from_log_h(log_h),
                "residual_log10_h": residual,
                "squared_error": residual * residual,
                "nll": bernoulli_nll(log_h, r["R"]),
                "training_mass_exponent": parameters["mass_exponent"],
            })
    predictions.sort(key=lambda p: p["species"])
    return {
        "rmse_log10_h": math.sqrt(statistics.fmean(p["squared_error"] for p in predictions)),
        "mean_nll": statistics.fmean(p["nll"] for p in predictions),
        "predictions": predictions,
    }


def delete_fold_sensitivity(predictions, reference):
    """Paired score differences, NOT a confidence interval or new refit CV."""
    ref = {r["species"]: r for r in reference}
    output = {}
    for key in ("squared_error", "nll"):
        differences = [p[key] - ref[p["species"]][key] for p in predictions]
        deleted = [statistics.fmean(v for j, v in enumerate(differences) if j != i)
                   for i in range(len(differences))]
        output[key] = {"mean": statistics.fmean(differences), "delete_fold_min": min(deleted),
                       "delete_fold_max": max(deleted)}
    return output


def model_specs():
    yield "intercept", "intercept", 0
    yield "mass", "mass", 0
    for family in ("burden", "cell", "free_mass"):
        for k in K_VALUES:
            yield f"{family}_k{k}", family, k


def analyze_panel(rows, include_order=False):
    out = {}
    for name, family, k in model_specs():
        parameters = fit(rows, family, k)
        cv = cross_validate(rows, family, k)
        entry = {"parameters": parameters, "loso": cv}
        if include_order:
            entry["leave_order_out"] = cross_validate(rows, family, k, "order")
        out[name] = entry
    for name in out:
        out[name]["vs_intercept"] = delete_fold_sensitivity(
            out[name]["loso"]["predictions"], out["intercept"]["loso"]["predictions"])
    return out


def rare_certificate(rows, k, scenario):
    intervals = []
    for r in rows:
        mlo = mhi = r["m"]
        llo = lhi = r["L"]
        m_mass_lo = m_mass_hi = r["M"]
        hlo, hhi = ((r["H"], r["H"]) if scenario == "point" else r["H_box"])
        if scenario == "factor_two":
            mlo, mhi = r["m"] / 2, r["m"] * 2
            llo, lhi = r["L"] / 2, r["L"] * 2
            m_mass_lo, m_mass_hi = r["M_box"]
        elif scenario not in ("point", "cmr95"):
            raise ValueError("unknown rare certificate scenario")
        lo = math.log(hlo) - math.log(m_mass_hi) - k * math.log(mhi * lhi)
        hi = math.log(hhi) - math.log(m_mass_lo) - k * math.log(mlo * llo)
        intervals.append((r["species"], lo, hi))
    low_witness = min(intervals, key=lambda v: (-v[1], v[0]))
    high_witness = min(intervals, key=lambda v: (v[2], v[0]))
    gap = max(0.0, low_witness[1] - high_witness[2])
    return {"rho": math.exp(gap / 2), "log_gap": gap,
            "lower_constraint_species": low_witness[0], "upper_constraint_species": high_witness[0],
            "common_amplitude_in_box": gap == 0.0}


def box(row, scenario):
    hlo, hhi = row["H_box"]
    mlo = mhi = row["m"]
    llo = lhi = row["L"]
    nlo = nhi = row["M"]
    if scenario == "point":
        hlo = hhi = row["H"]
    elif scenario == "cmr95":
        pass
    elif scenario in ("sex_mass", "slope_estimator", "slope_ci", "factor_two"):
        if scenario in ("sex_mass", "factor_two"):
            nlo, nhi = row["M_box"]
        if scenario in ("slope_estimator", "slope_ci"):
            mlo = mhi = row["slope"]
        if scenario == "slope_ci" and row["slope_ci"]:
            # Intersect a real-valued slope CI with the physical positive domain.
            # Zero is the closure boundary, NOT an invented positive cutoff.
            mlo = max(0.0, row["slope_ci"][0])
            mhi = row["slope_ci"][1]
        if scenario == "factor_two":
            mlo, mhi = row["m"] / 2, row["m"] * 2
            llo, lhi = row["L"] / 2, row["L"] * 2
    else:
        raise ValueError("unknown pair scenario")
    return {"H": (hlo, hhi), "M": (nlo, nhi), "B": (mlo * llo, mhi * lhi)}


def divide_nonnegative(numerator, denominator):
    return numerator / denominator if denominator else math.inf


def pair_certificate(larger, smaller, k, scenario):
    a, b = box(larger, scenario), box(smaller, scenario)
    ylo = a["H"][0] / b["H"][1] * b["M"][0] / a["M"][1]
    yhi = a["H"][1] / b["H"][0] * b["M"][1] / a["M"][0]
    blo = divide_nonnegative(a["B"][0], b["B"][1])
    bhi = divide_nonnegative(a["B"][1], b["B"][0])
    result = V.pair_modulation(k, ylo, yhi, blo, bhi)
    return {"y_interval": [ylo, yhi], "burden_ratio_interval": [blo, bhi], **result}


def exact_pairs(rows):
    output = []
    for a, b in itertools.combinations(rows, 2):
        larger, smaller = sorted((a, b), key=lambda r: (-r["M"], r["species"]))
        scenarios = {}
        for scenario in ("point", "cmr95", "sex_mass", "slope_estimator", "slope_ci", "factor_two"):
            scenarios[scenario] = {str(k): pair_certificate(larger, smaller, k, scenario) for k in K_VALUES}
        output.append({"larger": larger["species"], "smaller": smaller["species"],
                       "scenarios": scenarios})
    return output


def scaling_analysis(mutation, traits):
    species = sorted(set(mutation) & set(traits))
    x = [math.log(float(traits[s]["max_longevity_months"]) / 12) for s in species]
    y = [math.log(float(mutation[s]["mean_sbs_per_genome_year"])) for s in species]
    a, b = ols(x, y)
    out = {"n_species": len(species), "intercept": a, "slope": b, "models": {}}
    slopes = []
    for mode in ("free", "inverse", "intercept"):
        predictions = []
        for j, name in enumerate(species):
            tx, ty = [v for i, v in enumerate(x) if i != j], [v for i, v in enumerate(y) if i != j]
            if mode == "free":
                fa, fb = ols(tx, ty)
                slopes.append(fb)
            else:
                fb = -1 if mode == "inverse" else 0
                fa = statistics.fmean(yi - fb * xi for xi, yi in zip(tx, ty))
            error = (y[j] - fa - fb * x[j]) / LN10
            predictions.append({"species": name, "predicted_log10_rate": (fa + fb * x[j]) / LN10,
                                "error_log10_rate": error})
        out["models"][mode] = {"rmse_log10_rate": math.sqrt(statistics.fmean(
            p["error_log10_rate"] ** 2 for p in predictions)), "predictions": predictions}
    out["leave_species_out_slope_range"] = [min(slopes), max(slopes)]
    return out


def printed_pair(rows, base=HERE):
    lookup = {r["species"]: r for r in rows}
    for source in read_csv("cagan_printed_pair.csv", base):
        r = lookup[source["species"]].copy()
        r["m"] = float(source["mean_sbs_per_genome_year"])
        r["L"] = float(source["lifespan80_years"])
        lookup[r["species"]] = r
    a, b = lookup["Giraffa_camelopardalis"], lookup["Heterocephalus_glaber"]
    return {"burden_giraffe": a["m"] * a["L"], "burden_naked_mole_rat": b["m"] * b["L"],
            "exact_pair_cmr95": {str(k): pair_certificate(a, b, k, "sex_mass") for k in K_VALUES}}


def analyze(base=HERE):
    rows, mutation, traits = load_data(base)
    primary = panel_version(rows, "primary")
    panels = {name: analyze_panel(panel_version(rows, name), name == "primary")
              for name in ("primary", "maximum_longevity", "slope_estimator")}
    pairs = exact_pairs(primary)
    main_pair = next(p for p in pairs if p["larger"] == "Giraffa_camelopardalis"
                     and p["smaller"] == "Heterocephalus_glaber")
    outliers = []
    for p in panels["primary"]["cell_k3"]["loso"]["predictions"]:
        outliers.append({"species": p["species"], "residual_log10_h": p["residual_log10_h"],
                         "observed_over_predicted_h": 10 ** p["residual_log10_h"]})
    outliers.sort(key=lambda r: (-r["residual_log10_h"], r["species"]))
    return {
        "status": "nonclinical_unreviewed", "biological_verdict": "not_identifiable",
        "protocol_freeze_commit": FREEZE, "input_sha256": input_hashes(base),
        "coverage": {"mutation_rows": len(mutation), "life_history_rows": len(traits),
                     "matched_species": len(rows), "orders": len({r["order"] for r in rows}),
                     "cancer_deaths": sum(r["events"] for r in rows),
                     "known_cause_deaths": sum(r["n"] for r in rows),
                     "mutation_rows_without_mortality": sorted(set(mutation) - {r["species"] for r in rows}),
                     "matched_slope_intervals": sum(r["slope_ci"] is not None for r in rows),
                     "matched_slope_intervals_reaching_zero": sum(
                         r["slope_ci"] is not None and r["slope_ci"][0] <= 0 for r in rows)},
        "observations": primary, "panels": panels, "scaling": scaling_analysis(mutation, traits),
        "rare_certificates": {str(k): {s: rare_certificate(primary, k, s)
                                      for s in ("point", "cmr95", "factor_two")} for k in K_VALUES},
        "exact_pairs": pairs, "main_pair": main_pair,
        "printed_pair_sensitivity": printed_pair(rows, base), "outliers_k3": outliers,
        "exact_worst_pair_by_k": {str(k): max(
            ({"larger": p["larger"], "smaller": p["smaller"],
              "fold_infimum": p["scenarios"]["sex_mass"][str(k)]["fold_infimum"]} for p in pairs),
            key=lambda p: (p["fold_infimum"], p["larger"], p["smaller"])) for k in K_VALUES},
    }


def serializable(value):
    """JSON has no infinity; retain unbounded intervals with explicit strings."""
    if isinstance(value, float) and not math.isfinite(value):
        if math.isnan(value):
            raise ValueError("NaN is not a reportable result")
        return "Infinity" if value > 0 else "-Infinity"
    if isinstance(value, dict):
        return {key: serializable(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(v) for v in value]
    return value


def short_species(name):
    return name.replace("_", " ")


def render_results(r):
    """The README's entire empirical-number block is generated from CSVs."""
    lines = ["<!-- BEGIN COMPUTED RESULTS -->", "### Coverage and observed endpoints", ""]
    c = r["coverage"]
    lines.append(f"The readable mutation panel contains **{c['mutation_rows']} species**; "
                 f"**{c['matched_species']}** overlap the mortality table, spanning "
                 f"**{c['orders']} orders**, **{c['cancer_deaths']} cancer-attributed deaths** and "
                 f"**{c['known_cause_deaths']} known-cause deaths**. These are selected published "
                 "aggregates, not a random evolutionary sample.")
    lines += ["", "CMR intervals below are **calculated Wilson 95% intervals**, not published "
              "intervals. Rates are mean passenger SBS/genome/year; L is adult remaining "
              "expectancy; M is the average of published female/male masses.", "",
              "| Species | Annual SBS | L (years) | M (kg) | Cancer / known deaths | CMR (95% interval) |",
              "|---|---:|---:|---:|---:|---:|"]
    for q in r["observations"]:
        lines.append(f"| {short_species(q['species'])} | {q['m']:.4f} | {q['L']:.4f} | {q['M']:.6g} | "
                     f"{q['events']}/{q['n']} | {q['R']:.5f} [{q['R_box'][0]:.5f}, {q['R_box'][1]:.5f}] |")
    lines += ["", f"Only **{c['matched_slope_intervals']}** matched species have a printed "
              f"zero-intercept-slope interval; **{c['matched_slope_intervals_reaching_zero']}** "
              "reaches or crosses zero. These intervals do not belong to the mean-rate estimator. "
              "Per-species sequencing counts and uncertainty for the lifespan proxies remain missing.",
              "", "### Frozen prediction comparison", "",
              "Every k is reported. RMSE is measured in log10 H, with H = −log(1−CMR). "
              "All amplitudes and free mass exponents are refit inside each fold. "
              "The burden-only comparator deliberately omits the target-cell factor.", "",
              "| k | Cell-aware LOSO RMSE | Burden-only LOSO RMSE | Free-mass LOSO RMSE | "
              "Fitted free mass exponent | Cell-aware leave-order-out RMSE |",
              "|---:|---:|---:|---:|---:|---:|"]
    p = r["panels"]["primary"]
    for k in K_VALUES:
        cell, burden, free = (p[f"{f}_k{k}"] for f in ("cell", "burden", "free_mass"))
        lines.append(f"| {k} | {cell['loso']['rmse_log10_h']:.6f} | "
                     f"{burden['loso']['rmse_log10_h']:.6f} | {free['loso']['rmse_log10_h']:.6f} | "
                     f"{free['parameters']['mass_exponent']:.6f} | {cell['leave_order_out']['rmse_log10_h']:.6f} |")
    lines += ["", f"Intercept-only LOSO RMSE: **{p['intercept']['loso']['rmse_log10_h']:.6f}**; "
              f"mass-only: **{p['mass']['loso']['rmse_log10_h']:.6f}**.", "",
              "The next table uses the preregistered reporting index k=3, not a selected biological count. "
              "NLL is an equal-species mean of per-known-death Bernoulli log loss. "
              "The delete-fold range describes sensitivity of paired squared-error differences "
              "from the intercept comparator; it is not a confidence interval.", "",
              "| Panel | Cell RMSE | Burden RMSE | Free-mass RMSE | Cell NLL | Burden NLL | "
              "Cell minus intercept MSE; delete-fold range |",
              "|---|---:|---:|---:|---:|---:|---:|"]
    for name, pp in r["panels"].items():
        a, b, f = (pp[f"{kind}_k3"] for kind in ("cell", "burden", "free_mass"))
        d = a["vs_intercept"]["squared_error"]
        lines.append(f"| {name} | {a['loso']['rmse_log10_h']:.6f} | {b['loso']['rmse_log10_h']:.6f} | "
                     f"{f['loso']['rmse_log10_h']:.6f} | {a['loso']['mean_nll']:.6f} | "
                     f"{b['loso']['mean_nll']:.6f} | {d['mean']:.6f} "
                     f"[{d['delete_fold_min']:.6f}, {d['delete_fold_max']:.6f}] |")
    lines += ["", "All fixed-k scores, coefficients, NLLs, fold predictions and sensitivity ranges "
              "are in `analysis_results.json`; no favorable k is hidden by this compact display.", "",
              "### Compatibility and required residual modulation", "",
              "Rare-family rho is a symmetric per-species hazard allowance about one common "
              "amplitude. Exact-pair F is the relative modulation fold between giraffe and "
              "naked mole rat, with CMR intervals and both-sex mass endpoints. These are "
              "different conventions: pairwise symmetric rho is sqrt(F). Finite-rate "
              "boundary values are infima, not necessarily attained minima.", "",
              "| k | Rare rho, points | Rare rho, CMR intervals | Rare rho, factor-two stress | "
              "Exact-pair F, sex-mass/CMR boxes | Exact-pair F, slope-CI subset | Exact-pair F, factor-two stress |",
              "|---:|---:|---:|---:|---:|---:|---:|"]
    pair = r["main_pair"]["scenarios"]
    for k in K_VALUES:
        rc = r["rare_certificates"][str(k)]
        lines.append(f"| {k} | {rc['point']['rho']:.6g} | {rc['cmr95']['rho']:.6g} | "
                     f"{rc['factor_two']['rho']:.6g} | {pair['sex_mass'][str(k)]['fold_infimum']:.6g} | "
                     f"{pair['slope_ci'][str(k)]['fold_infimum']:.6g} | "
                     f"{pair['factor_two'][str(k)]['fold_infimum']:.6g} |")
    q = pair["sex_mass"]["3"]
    lines += ["", f"At k=3, the giraffe/naked-mole-rat effective hazard multiplier ratio has "
              f"closure **[{q['multiplier_low']:.8g}, {q['multiplier_high']:.8g}]**. "
              "This bounds a missing model factor under the proxy assumptions; it does not "
              "measure an immune mechanism or extra biological stages."]
    printed = r["printed_pair_sensitivity"]
    lines += ["", f"Cagan's separately printed rounded lifespan/rate pair gives burdens "
              f"**{printed['burden_giraffe']:.0f}** (giraffe) and "
              f"**{printed['burden_naked_mole_rat']:.0f}** (naked mole rat). Using these "
              "as horizons with Vincze's mortality/mass boxes remains a cross-source proxy; "
              f"its k=3 relative-modulation infimum is "
              f"**{printed['exact_pair_cmr95']['3']['fold_infimum']:.6g}**.", "",
              "The slope-CI scenario fixes rates lacking intervals at their slope point estimates; "
              "it is incomplete uncertainty propagation. The factor-two box is hypothetical. "
              "Compatibility inside a box is not confirmation; failure is not a calibrated "
              "sampling test. Pairwise success does not establish one common calibration for "
              "the entire panel.", "", "### Outliers and lifespan sensitivity", ""]
    hi, lo = r["outliers_k3"][0], r["outliers_k3"][-1]
    lines.append(f"The largest positive cell-aware k=3 LOSO residual is "
                 f"**{short_species(hi['species'])}**, {hi['residual_log10_h']:.6f} log10 units "
                 f"({hi['observed_over_predicted_h']:.6g}× observed/predicted H). "
                 f"The largest negative is **{short_species(lo['species'])}**, "
                 f"{lo['residual_log10_h']:.6f} ({lo['observed_over_predicted_h']:.6g}×). "
                 "These are model residuals, not species-wide rankings of cancer resistance.")
    s = r["scaling"]
    lines += ["", f"In the **{s['n_species']}-species maximum-longevity** panel, log annual rate "
              f"on log longevity has slope **{s['slope']:.6f}**, with leave-species-out slope "
              f"range **[{s['leave_species_out_slope_range'][0]:.6f}, "
              f"{s['leave_species_out_slope_range'][1]:.6f}]**. "
              f"LOSO log10-rate RMSE is **{s['models']['free']['rmse_log10_rate']:.6f}** "
              f"for the free slope, **{s['models']['inverse']['rmse_log10_rate']:.6f}** for "
              f"fixed inverse slope, and **{s['models']['intercept']['rmse_log10_rate']:.6f}** "
              "for intercept-only. This does not reproduce Cagan's original lifespan definition "
              "or complete species panel.", "", "### Input integrity", "",
              "SHA-256 uses newline-normalized UTF-8 content and POSIX relative path keys. "
              "Hashes establish local integrity, not truth or biological validity.", "",
              "| Input | SHA-256 |", "|---|---|"]
    for path, digest in r["input_sha256"].items():
        lines.append(f"| `{path}` | `{digest}` |")
    lines += ["", "<!-- END COMPUTED RESULTS -->"]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    result = analyze()
    body = json.dumps(serializable(result), indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(body, encoding="utf-8", newline="\n")
    else:
        print(render_results(result))
    if args.markdown:
        args.markdown.write_text(render_results(result) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
