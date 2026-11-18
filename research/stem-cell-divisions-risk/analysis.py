"""Reproduce the frozen aggregate reanalysis. Nonclinical; no medical advice.

Only the Python standard library is used. There are no network calls. Every
reported empirical table is rendered from the CSVs, never from hand-entered
regression results. Models, folds and sensitivity sets are in PROTOCOL.md.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import random
import statistics
import sys


ROOT = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("stem_divisions_verifier", ROOT / "verifier.py")
V = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = V
_spec.loader.exec_module(V)

LN10 = math.log(10)
BOOTSTRAPS = 2000
SEED = 1260825
COUNTS = tuple(range(1, 9))
GROWTH_GRID = (0, .25, .5, 1, 2, 4, 8)
ALLOWANCES = (1, 2, 5, 10, 100)
NUMERIC_FIELDS = ("risk", "stem_cells", "divisions_per_year",
                  "divisions_per_lifetime", "total_divisions")
PROVENANCE_FIELDS = ("source", "source_url", "source_locator", "source_version",
                     "access_date", "reuse_basis", "reuse_terms_url")
PREREG_COMMIT = "f5c1f4e87bfd03c85eb015fa14e72a28d33844ad"
AMENDMENT_COMMIT = "c7d255d7bfafe0ac720e1752a4b55058d136ecd3"


def relative_posix_key(path, root):
    """Use stable POSIX separators for artifact keys, including on Windows."""
    return path.relative_to(root).as_posix()


def normalized_hash(path):
    """Hash text after universal newline conversion, without trimming content."""
    return hashlib.sha256(Path(path).read_text(encoding="utf-8").encode("utf-8")).hexdigest()


def load_rows(name="tv2015.csv"):
    with (ROOT / "data" / name).open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        if not all(row.get(key, "").strip() for key in PROVENANCE_FIELDS):
            raise ValueError("Missing per-row provenance")
        row["high_risk_subgroup"] = int(row["high_risk_subgroup"])
        for key in NUMERIC_FIELDS:
            row[key] = float(row[key])
        if not 0 < row["risk"] <= 1 or min(row["stem_cells"], row["total_divisions"]) <= 0:
            raise ValueError("Invalid source point")
    if len({row["tissue_id"] for row in rows}) != len(rows):
        raise ValueError("Duplicate tissue identifier")
    return rows


def clock(row, kind="effective"):
    if kind == "effective":
        return row["total_divisions"] / row["stem_cells"]
    if kind == "developmental":
        return row["divisions_per_lifetime"] + math.log2(row["stem_cells"])
    raise ValueError(kind)


def hazard(row):
    if not 0 < row["risk"] < 1:
        raise ValueError("Finite cumulative hazard requires 0<R<1")
    return -math.log1p(-row["risk"])


def centered_cross(xs, ys):
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    return math.fsum((x-mx)*(y-my) for x, y in zip(xs, ys))


def simple_fit(xs, ys):
    xx = centered_cross(xs, xs)
    if xx <= 1e-24:
        raise ValueError("Rank-deficient univariate regression")
    xy, yy = centered_cross(xs, ys), centered_cross(ys, ys)
    slope = xy/xx
    intercept = statistics.fmean(ys)-slope*statistics.fmean(xs)
    residuals = [y-intercept-slope*x for x, y in zip(xs, ys)]
    return {"intercept": intercept, "slope": slope,
            "correlation": xy/math.sqrt(xx*yy) if yy else 0.0,
            "r_squared": xy*xy/(xx*yy) if yy else 0.0,
            "rmse": math.sqrt(statistics.fmean(r*r for r in residuals))}


def separated_fit(rows, kind):
    xs = [math.log10(r["stem_cells"]) for r in rows]
    vs = [math.log10(clock(r, kind)) for r in rows]
    ys = [math.log10(hazard(r)) for r in rows]
    xx, vv, xv = centered_cross(xs, xs), centered_cross(vs, vs), centered_cross(xs, vs)
    xy, vy = centered_cross(xs, ys), centered_cross(vs, ys)
    det = xx*vv-xv*xv
    if det <= 1e-12*max(xx*vv, 1):
        raise ValueError("Rank-deficient separable regression")
    bn, bv = (xy*vv-vy*xv)/det, (vy*xx-xy*xv)/det
    a = statistics.fmean(ys)-bn*statistics.fmean(xs)-bv*statistics.fmean(vs)
    return {"intercept": a, "b_N": bn, "b_v": bv}


def log_shape(row, k, kind="effective", gamma=0.0):
    v = clock(row, kind)
    return (math.log(row["stem_cells"])+k*math.log(v)-math.lgamma(k+1)
            +V.log_phi(k, gamma*v))/LN10


def fit(rows, family, kind="effective", k=None):
    """Training-only fitting, including training-only scaling of growth grid."""
    if family == "descriptive":
        p = simple_fit([math.log10(r["total_divisions"]) for r in rows],
                       [math.log10(hazard(r)) for r in rows])
    elif family == "separable":
        p = separated_fit(rows, kind)
    elif family in ("sequential", "expansion"):
        cs = GROWTH_GRID if family == "expansion" else (0,)
        best = None
        for c in cs:
            gamma = c/max(clock(r, kind) for r in rows)
            residuals = [math.log10(hazard(r))-log_shape(r, k, kind, gamma) for r in rows]
            a = statistics.fmean(residuals)
            sse = math.fsum((r-a)**2 for r in residuals)
            candidate = (sse, c, {"intercept": a, "gamma": gamma, "grid_c": c})
            if best is None or candidate[:2] < best[:2]:
                best = candidate
        p = best[2]
    else:
        raise ValueError(family)
    return {"family": family, "clock": kind, "k": k, **p}


def predict(model, row):
    family = model["family"]
    if family == "descriptive":
        return model["intercept"]+model["slope"]*math.log10(row["total_divisions"])
    if family == "separable":
        return (model["intercept"]+model["b_N"]*math.log10(row["stem_cells"])
                +model["b_v"]*math.log10(clock(row, model["clock"])))
    return model["intercept"]+log_shape(row, model["k"], model["clock"], model["gamma"])


def cross_validate(rows, family, kind="effective", k=None, grouped=True):
    key = "tissue_group" if grouped else "tissue_id"
    folds = sorted({r[key] for r in rows})
    predictions, fold_mse = [], []
    for fold in folds:
        training = [r for r in rows if r[key] != fold]
        testing = [r for r in rows if r[key] == fold]
        model = fit(training, family, kind, k)
        sq = []
        for row in testing:
            estimate = predict(model, row)
            residual = math.log10(hazard(row))-estimate
            sq.append(residual**2)
            predictions.append({"tissue_id": row["tissue_id"], "tissue": row["tissue"],
                                "group": row["tissue_group"], "prediction": estimate,
                                "residual": residual, "amplitude_ratio": 10**residual,
                                "training_grid_c": model.get("grid_c")})
        fold_mse.append(statistics.fmean(sq))
    return {"fold_count": len(folds),
            "group_balanced_rmse": math.sqrt(statistics.fmean(fold_mse)),
            "row_rmse": math.sqrt(statistics.fmean(p["residual"]**2 for p in predictions)),
            "predictions": sorted(predictions, key=lambda x: x["tissue_id"])}


def score_model(rows, family, kind="effective", k=None):
    model = fit(rows, family, kind, k)
    residuals = [math.log10(hazard(r))-predict(model, r) for r in rows]
    result = {"fit": model, "training_rmse": math.sqrt(statistics.fmean(r*r for r in residuals)),
              "leave_group_out": cross_validate(rows, family, kind, k),
              "leave_row_out": cross_validate(rows, family, kind, k, grouped=False)}
    if family == "sequential":
        u = 10**(model["intercept"]/k)
        discrepancies = [abs(math.log10(r["stem_cells"]*V.psi(k, u*clock(r, kind)))
                             -predict(model, r)) for r in rows]
        x = [math.log10(r["total_divisions"]) for r in rows]
        y = [predict(model, r) for r in rows]
        log_r = [math.log10(-math.expm1(-10**z)) for z in y]
        result.update({"implied_equal_rate": u,
                       "max_exact_logH_discrepancy": max(discrepancies),
                       "predicted_logH_slope": simple_fit(x, y)["slope"],
                       "predicted_logR_slope": simple_fit(x, log_r)["slope"]})
    return result


def quantile(values, q):
    values = sorted(values)
    pos = (len(values)-1)*q
    lo, hi = math.floor(pos), math.ceil(pos)
    return values[lo]+(pos-lo)*(values[hi]-values[lo])


def slope_bootstrap(rows):
    rng = random.Random(SEED)
    groups = sorted({r["tissue_group"] for r in rows})
    members = {g: [r for r in rows if r["tissue_group"] == g] for g in groups}
    slopes = []
    skipped = 0
    for _ in range(BOOTSTRAPS):
        sample = [r for _ in groups for r in members[rng.choice(groups)]]
        try:
            slopes.append(simple_fit([math.log10(r["total_divisions"]) for r in sample],
                                     [math.log10(r["risk"]) for r in sample])["slope"])
        except ValueError:
            skipped += 1
    return {"draws": BOOTSTRAPS, "skipped": skipped, "seed": SEED,
            "lower": quantile(slopes, .025), "upper": quantile(slopes, .975)}


def rate_certificate(rows, k, kind="effective", stress=False):
    intervals, points, rare = [], [], []
    for row in rows:
        R, N, D = row["risk"], row["stem_cells"], row["total_divisions"]
        v = clock(row, kind)
        u = V.inverse_psi(k, hazard(row)/N)/v
        points.append(u)
        rare.append(math.exp((math.lgamma(k+1)+math.log(hazard(row)/N))/k)/v)
        if stress:
            if kind != "effective":
                raise ValueError("Stress certificate only preregistered for D/N")
            lower = V.rate(k, R/1.25, N/2, D*2)
            upper = V.rate(k, min(1-1e-12, R*1.25), N*2, D/2)
            intervals.append((lower, upper))
        else:
            intervals.append((u, u))
    max_lower = max(a for a, _ in intervals)
    min_upper = min(b for _, b in intervals)
    rho = math.sqrt(max(1.0, max_lower/min_upper))
    return {"k": k, "rho": rho, "rate_spread": rho*rho,
            "limiting_point_rho": math.sqrt(max(rare)/min(rare)),
            "baseline_interval": [max_lower/rho, rho*min_upper],
            "lower_constraint": [r["tissue_id"] for r, ab in zip(rows, intervals) if ab[0] == max_lower],
            "upper_constraint": [r["tissue_id"] for r, ab in zip(rows, intervals) if ab[1] == min_upper],
            "compatible_allowances": [a for a in ALLOWANCES if a >= rho],
            "max_relative_exact_vs_rare_rate": max(u/a-1 for u, a in zip(points, rare)),
            "rates": [{"tissue_id": r["tissue_id"], "rate": u, "interval": list(ab)}
                      for r, u, ab in zip(rows, points, intervals)]}


def panel(rows, kind="effective", expansion=True, bootstrap=False):
    x = [math.log10(r["total_divisions"]) for r in rows]
    y = [math.log10(r["risk"]) for r in rows]
    v = [math.log10(clock(r, kind)) for r in rows]
    result = {"n": len(rows), "groups": len({r["tissue_group"] for r in rows}),
              "ids": [r["tissue_id"] for r in rows],
              "logR_fit": simple_fit(x, y),
              "geometry_covariance_ratio": centered_cross(x, v)/centered_cross(x, x),
              "descriptive": score_model(rows, "descriptive", kind),
              "separable": score_model(rows, "separable", kind),
              "sequential": [score_model(rows, "sequential", kind, k) for k in COUNTS],
              "exact_frontier": [rate_certificate(rows, k, kind) for k in COUNTS]}
    if expansion:
        result["expansion"] = [score_model(rows, "expansion", kind, k) for k in COUNTS if k >= 2]
    if kind == "effective":
        result["stress_frontier"] = [rate_certificate(rows, k, kind, stress=True) for k in COUNTS]
    if bootstrap:
        result["slope_bootstrap"] = slope_bootstrap(rows)
    ordered = sorted(rows, key=lambda r: r["total_divisions"])
    violations = []
    for a, b in itertools.combinations(ordered, 2):
        if a["total_divisions"] == b["total_divisions"]:
            continue
        if a["stem_cells"] > b["stem_cells"] or clock(a, kind) > clock(b, kind):
            violations.append([a["tissue_id"], b["tissue_id"]])
    hmax = max(hazard(r) for r in rows)
    result["saturation"] = {"max_risk": max(r["risk"] for r in rows),
                            "sigma_at_max_risk": hmax/math.expm1(hmax),
                            "ordering_violations": len(violations),
                            "ordered_cloud_bound_applicable": not violations}
    residuals = result["descriptive"]["leave_group_out"]["predictions"]
    result["outliers_negative"] = sorted(residuals, key=lambda r: (r["residual"], r["tissue_id"]))[:3]
    result["outliers_positive"] = sorted(residuals, key=lambda r: (-r["residual"], r["tissue_id"]))[:3]
    return result


def duplicate_contrasts(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["stem_cells"], row["total_divisions"])].append(row)
    result = []
    for members in groups.values():
        for a, b in itertools.combinations(sorted(members, key=lambda r: (r["risk"], r["tissue_id"])), 2):
            if a["risk"] == b["risk"]:
                continue
            ratios = None
            if b["risk"] < 1:
                ratios = {str(k): V.rate(k, b["risk"], b["stem_cells"], b["total_divisions"])
                          /V.rate(k, a["risk"], a["stem_cells"], a["total_divisions"]) for k in COUNTS}
            result.append({"lower": a["tissue_id"], "higher": b["tissue_id"],
                           "risk_ratio": b["risk"]/a["risk"], "rate_ratios": ratios,
                           "hazard_ratio": hazard(b)/hazard(a) if b["risk"] < 1 else None})
    return sorted(result, key=lambda p: (p["lower"], p["higher"]))


def audit_versions(original, revised):
    if len(original) != len(revised):
        raise ValueError("Source versions have different row counts")
    if {r["tissue_id"] for r in original} != {r["tissue_id"] for r in revised}:
        raise ValueError("Source versions have different tissue identifiers")
    changes = []
    for a, b in zip(original, revised):
        if a["tissue_id"] != b["tissue_id"]:
            raise ValueError("Version row alignment changed")
        for field in NUMERIC_FIELDS:
            if a[field] != b[field]:
                changes.append({"tissue_id": a["tissue_id"], "field": field,
                                "tv2015": a[field], "little2016": b[field]})
    return {"numeric_fields_compared": len(original)*len(NUMERIC_FIELDS), "changes": changes,
            "annual_lifetime_inconsistencies": [
                {"tissue_id": r["tissue_id"], "annual_times_80": r["divisions_per_year"]*80,
                 "source_lifetime": r["divisions_per_lifetime"]}
                for r in original if not math.isclose(r["divisions_per_year"]*80,
                                                      r["divisions_per_lifetime"], rel_tol=.01, abs_tol=.01)]}


def international_endpoint_audit():
    """Later risk-only audit: never join independently ranked columns as people/registries."""
    path = ROOT/"data"/"liu2022_ranked_risk_audit.csv"
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    output = []
    for row in rows:
        if not all(row.get(key, "").strip() for key in PROVENANCE_FIELDS):
            raise ValueError("Missing audit provenance")
        low, middle, high = (float(row[f"risk_{field}"]) for field in ("min", "median", "max"))
        if not 0 <= low <= middle <= high < 1:
            raise ValueError("Invalid ranked risk range")
        output.append({"tissue_id": row["tissue_id"], "tissue": row["tissue"],
                       "risk_min": low, "risk_median": middle, "risk_max": high,
                       "registry_count": int(row["registry_count"]),
                       "age_endpoint": row["age_endpoint"]})
    return {"status": "risk-only endpoint audit; no coherent risk-plus-divisions fit",
            "n": len(output), "rows": output}


def analyze():
    original, revised = load_rows(), load_rows("little2016.csv")
    corrected = load_rows("tv2015_nunney2020.csv")
    primary = [r for r in original if not r["high_risk_subgroup"] and r["risk"] < 1]
    finite = [r for r in original if r["risk"] < 1]
    later = [r for r in revised if not r["high_risk_subgroup"] and r["risk"] < 1]
    corrected_primary = [r for r in corrected if not r["high_risk_subgroup"] and r["risk"] < 1]
    return {"status": "nonclinical; unreviewed; no medical advice",
            "preregistration_commit": PREREG_COMMIT,
            "prospective_amendment_commit": AMENDMENT_COMMIT,
            "hashes": {relative_posix_key(p, ROOT): normalized_hash(p) for p in
                       [ROOT/"PROTOCOL.md", ROOT/"AMENDMENTS.md", *sorted((ROOT/"data").glob("*.csv"))]},
            "historical_all": {"n": len(original),
                **simple_fit([math.log10(r["total_divisions"]) for r in original],
                             [math.log10(r["risk"]) for r in original]),
                "bootstrap": slope_bootstrap(original)},
            "primary": panel(primary, bootstrap=True),
            "all_finite": panel(finite),
            "primary_developmental": panel(primary, "developmental", expansion=False),
            "later_effective": panel(later),
            "later_developmental": panel(later, "developmental", expansion=False),
            "corrected_effective": panel(corrected_primary, expansion=False),
            "duplicate_contrasts": duplicate_contrasts(original),
            "international_endpoint_audit": international_endpoint_audit(),
            "source_version_audit": audit_versions(original, revised)}


def compact(results):
    """Canonical output precision; calculations and model selection are unrounded."""
    if isinstance(results, float):
        if not math.isfinite(results):
            raise ValueError("Nonfinite result")
        return float(format(results, ".12g"))
    if isinstance(results, dict):
        return {k: compact(v) for k, v in results.items()}
    if isinstance(results, (tuple, list)):
        return [compact(v) for v in results]
    return results


def artifact_results(results):
    """Publish all reported numbers without duplicating intermediate fold arrays.

    analyze() retains every prediction and inverse for the independent tests.
    The archived summary omits only arrays named predictions/rates; outliers,
    every model score, coefficients, frontier witnesses and source audit remain.
    """
    if isinstance(results, dict):
        return {key: artifact_results(value) for key, value in results.items()
                if key not in ("predictions", "rates")}
    if isinstance(results, (tuple, list)):
        return [artifact_results(value) for value in results]
    return results


def f(value, places=6):
    return f"{value:.{places}f}"


def render_results(r):
    """Complete empirical report block checked verbatim by the dedicated test."""
    p = r["primary"]
    lines = ["<!-- BEGIN RECOMPUTED RESULTS -->", "### Published-data replication and panels", "",
             "All numbers in this block are regenerated from the cited CSVs by `analysis.py`.",
             "Errors and coefficients use base-10 logs. Group-held-out RMSE is in log H;",
             "the primary score gives each anatomical group equal weight.", "",
             "| Panel | Rows | Groups | log R slope | Correlation |",
             "|---|---:|---:|---:|---:|"]
    hist = r["historical_all"]
    lines.append(f"| Historical all categories | {hist['n']} | — | {f(hist['slope'])} | {f(hist['correlation'])} |")
    for name, label in (("primary", "Primary"), ("all_finite", "All finite risks"),
                        ("later_effective", "Later reanalysis, primary selection"),
                        ("corrected_effective", "Corrected development count, primary selection")):
        a = r[name]
        lines.append(f"| {label} | {a['n']} | {a['groups']} | {f(a['logR_fit']['slope'])} | {f(a['logR_fit']['correlation'])} |")
    lines += ["", f"Historical R²: **{f(hist['r_squared'])}**. Conditional group-bootstrap slope intervals:"]
    for label, a in (("historical", hist["bootstrap"]), ("primary", p["slope_bootstrap"])):
        lines += [f"{label}: **[{f(a['lower'])}, {f(a['upper'])}]** ({a['draws']} draws; {a['skipped']} degenerate draws omitted)."]
    lines += ["", "These intervals resample this selected tissue collection; they do not quantify",
              "stem-cell measurement error or epidemiological sampling precision.", "",
              "### Prediction under the primary effective clock", "",
              "| Model | Training RMSE | Leave-row-out RMSE | Leave-tissue-group-out RMSE | Group CV, row weighted |",
              "|---|---:|---:|---:|---:|"]
    model_rows = [("Descriptive log D", p["descriptive"]), ("Separable log N, log v", p["separable"])]
    model_rows += [(f"Sequential k={i+1}", a) for i, a in enumerate(p["sequential"])]
    model_rows += [(f"Expansion k={a['fit']['k']}", a) for a in p["expansion"]]
    for label, a in model_rows:
        lines.append(f"| {label} | {f(a['training_rmse'])} | {f(a['leave_row_out']['row_rmse'])} | {f(a['leave_group_out']['group_balanced_rmse'])} | {f(a['leave_group_out']['row_rmse'])} |")
    sep = p["separable"]["fit"]
    lines += ["", f"Separable coefficients: **b_N={f(sep['b_N'])}, b_v={f(sep['b_v'])}**.",
              f"Geometry covariance ratio Cov(log D,log v)/Var(log D): **{f(p['geometry_covariance_ratio'])}**.", "",
              "| k | Fixed-k predicted log H slope | Fixed-k predicted log R slope | Exact-chain discrepancy, max log H | Expansion full-fit grid c |",
              "|---:|---:|---:|---:|---:|"]
    for i, a in enumerate(p["sequential"]):
        c = "—" if i == 0 else f(p["expansion"][i-1]["fit"]["grid_c"], 2)
        lines.append(f"| {i+1} | {f(a['predicted_logH_slope'])} | {f(a['predicted_logR_slope'])} | {f(a['max_exact_logH_discrepancy'])} | {c} |")
    sat = p["saturation"]
    lines += ["", f"Maximum primary risk: **{f(sat['max_risk'])}**; its attenuation factor is **{f(sat['sigma_at_max_risk'])}**.",
              f"There are **{sat['ordering_violations']}** pairs violating the joint N/v ordering condition.",
              "The ordered-cloud theorem therefore cannot be used as a bound on this observed cloud.", "",
              "### Exact rate-heterogeneity certificate", "",
              "Here rho bounds each tissue intensity within [u0/rho,u0·rho]. It is not an",
              "environmental fraction, a measured mutation-rate spread, or a confidence interval.", "",
              "| k | Exact point rho | Rare-limit point rho | Stress-box rho | Upper-rate constraint | Lower-rate constraint | Compatible point allowances |",
              "|---:|---:|---:|---:|---|---|---|"]
    for a, b in zip(p["exact_frontier"], p["stress_frontier"]):
        lines.append(f"| {a['k']} | {f(a['rho'])} | {f(a['limiting_point_rho'])} | {f(b['rho'])} | {', '.join(a['lower_constraint'])} | {', '.join(a['upper_constraint'])} | {', '.join(map(str,a['compatible_allowances'])) or 'none'} |")
    lines += ["", "### Source and clock sensitivity", "",
              "| Panel / clock | k | Sequential group CV RMSE | Exact point rho | Stress-box rho |",
              "|---|---:|---:|---:|---:|"]
    for name, label in (("all_finite", "All finite / D/N"),
                        ("primary_developmental", "Original / L+log2 N"),
                        ("later_effective", "Later / D/N"),
                        ("later_developmental", "Later / L+log2 N"),
                        ("corrected_effective", "Corrected development / D/N")):
        for index, (a, b) in enumerate(zip(r[name]["sequential"], r[name]["exact_frontier"])):
            stress = f(r[name]["stress_frontier"][index]["rho"]) if "stress_frontier" in r[name] else "—"
            lines.append(f"| {label} | {b['k']} | {f(a['leave_group_out']['group_balanced_rmse'])} | {f(b['rho'])} | {stress} |")
    lines += ["", "### Prespecified outliers", "",
              "These are group-held-out residuals of the descriptive log-H regression.",
              "The ratio is observed H divided by predicted H; it is a required amplitude",
              "adjustment, not a diagnosed mechanism or a success probability.", "",
              "| Direction | Tissue / category | Held-out residual | H observed / H predicted |",
              "|---|---|---:|---:|"]
    for label, key in (("Positive", "outliers_positive"), ("Negative", "outliers_negative")):
        for a in p[key]:
            lines.append(f"| {label} | {a['tissue']} | {f(a['residual'])} | {f(a['amplitude_ratio'])} |")
    lines += ["", "### Identical-input contrasts", "",
              "All distinct-risk pairs with identical published N and D are retained.",
              "An endpoint printed as one is not clipped to create a finite inverse.", "",
              "| Lower-risk category | Higher-risk category | Risk ratio | H ratio | Required intensity ratio, k=2 |",
              "|---|---|---:|---:|---:|"]
    for a in r["duplicate_contrasts"]:
        hr = "undefined at R=1" if a["hazard_ratio"] is None else f(a["hazard_ratio"])
        ur = "undefined at R=1" if a["rate_ratios"] is None else f(a["rate_ratios"]["2"])
        lines.append(f"| {a['lower']} | {a['higher']} | {f(a['risk_ratio'])} | {hr} | {ur} |")
    audit = r["source_version_audit"]
    lines += ["", f"The version audit compares **{audit['numeric_fields_compared']}** numeric fields and finds **{len(audit['changes'])}** changes:", "",
              "| Tissue | Field | Original | Later reanalysis |", "|---|---|---:|---:|"]
    for a in audit["changes"]:
        lines.append(f"| {a['tissue_id']} | {a['field']} | {a['tv2015']:g} | {a['little2016']:g} |")
    lines += ["", "Annual/lifetime fields exceeding the declared arithmetic check tolerance are preserved:", "",
              "| Tissue | Annual divisions × 80 | Published lifetime divisions |", "|---|---:|---:|"]
    for a in audit["annual_lifetime_inconsistencies"]:
        lines.append(f"| {a['tissue_id']} | {a['annual_times_80']:g} | {a['source_lifetime']:g} |")
    endpoint = r["international_endpoint_audit"]
    lines += ["", "### Later international risk-only endpoint audit", "",
              f"The permissively reusable later table supplies **{endpoint['n']}** cancer categories.",
              "These are independently ranked registry risks through age 74, not matched countries,",
              "not full lifetime risks, and not a joint risk-versus-divisions replication.", "",
              "| Cancer category | Registry minimum | Registry median | Registry maximum | Registry entries |",
              "|---|---:|---:|---:|---:|"]
    for a in endpoint["rows"]:
        lines.append(f"| {a['tissue']} | {a['risk_min']:.8g} | {a['risk_median']:.8g} | {a['risk_max']:.8g} | {a['registry_count']} |")
    lines += ["", "### Normalized content hashes", "", "| File | SHA-256 after newline normalization |", "|---|---|"]
    for name, digest in sorted(r["hashes"].items()):
        lines.append(f"| {name} | `{digest}` |")
    lines += ["", "<!-- END RECOMPUTED RESULTS -->"]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    results = compact(artifact_results(analyze()))
    serialized = json.dumps(results, indent=2, sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(serialized, encoding="utf-8", newline="\n")
    if args.markdown:
        args.markdown.write_text(render_results(results)+"\n", encoding="utf-8", newline="\n")
    if not args.output and not args.markdown:
        print(render_results(results))


if __name__ == "__main__":
    main()
