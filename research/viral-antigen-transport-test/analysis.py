"""NONCLINICAL frozen-parameter aggregate transport audit; Python stdlib only.

No fitting, source downloads, inferred patient records or treatment decisions.
The input medians and ORRs are separate published marginal summaries.
"""

import argparse
import csv
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PR72 = ROOT / "research" / "merkel-viral-antigen-test"
PR67 = ROOT / "research" / "immune-escape-pancancer-test"
BEGIN = "<!-- BEGIN RECOMPUTED RESULTS -->"
END = "<!-- END RECOMPUTED RESULTS -->"
PROTOCOL_SHA256 = "75c9e5a87ad2f3f132efd3f52491ccba5a142494c0f34f7c006ae047bd701ce9"
PINNED = {
    (HERE / "PROTOCOL.md").relative_to(ROOT).as_posix(): PROTOCOL_SHA256,
    (PR72 / "data" / "frozen_curves.csv").relative_to(ROOT).as_posix():
        "1acc4eff8810b2afd59c6e5cac250705f2e4e4f50b1c5f079a0c18308bb893f9",
    (PR72 / "data" / "part_a_calibration.csv").relative_to(ROOT).as_posix():
        "ef597c8dd3ff770848fb06092853e620136a014c1fdd76c9e3873618cbf3af6e",
    (PR67 / "data" / "tmb_orr.csv").relative_to(ROOT).as_posix():
        "b34c1dd8d8d2c403f7251cf02611283d5326cbdc2caf2c72a6afadcf5ddc540b",
    (PR67 / "analysis.py").relative_to(ROOT).as_posix():
        "06a7e15a24fa23b52f4c74b3659c678dd2778a88fdb429b0712e170e4b564b1c",
}
PROVENANCE = ("source_title", "source_doi", "source_url", "locator",
              "access_date", "reuse_basis", "terms_url", "inspection")
REQUIRED = {"hpv_hnscc", "ebv_gastric"}


def canonical_bytes(data):
    return data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def sha256(path):
    return hashlib.sha256(canonical_bytes(Path(path).read_bytes())).hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def number(value):
    return None if value == "" else float(value)


def half_unit(value):
    """Half the last printed decimal unit; preserve source strings in CSVs."""
    return float(Decimal(5).scaleb(Decimal(value).as_tuple().exponent - 1))


def rounded_bounds(value, divisor=1.0, ceiling=math.inf):
    center = float(value) / divisor
    half = half_unit(value) / divisor
    return max(0.0, center - half), min(ceiling, center + half)


def load_frozen():
    for key, expected in PINNED.items():
        if sha256(ROOT / key) != expected:
            raise ValueError("Frozen input changed: " + key)
    curves = {r["name"]: r for r in read_csv(PR72 / "data" / "frozen_curves.csv")}
    reciprocal = curves["reciprocal_holdout"]
    log = curves["log_holdout"]
    constants = {"A": float(reciprocal["A"]), "B": float(reciprocal["B"]),
                 "a": float(log["a"]), "b": float(log["b"])}
    scenarios = []
    for r in read_csv(PR72 / "data" / "part_a_calibration.csv"):
        if r["calibration_set"] != "broad":
            continue
        if float(r["A"]) != constants["A"] or float(r["B"]) != constants["B"]:
            raise ValueError("Offset calibrated for a different curve")
        scenarios.append({"scale_label": r["scale_label"], "c": float(r["scale"]),
                          "target": r["target"], "v": float(r["v"])})
    return constants, scenarios


def probability(load, constants):
    if load < 0 or not math.isfinite(load):
        raise ValueError("Finite nonnegative effective load required")
    return 0.0 if load == 0 else constants["A"] * math.exp(-constants["B"] / load)


def predictions(median, positive, scenario, constants):
    load = scenario["c"] * median
    base = probability(load, constants)
    offset = probability(load + (scenario["v"] if positive else 0), constants)
    log = 0.0 if load == 0 else min(1.0, max(0.0,
        constants["a"] + constants["b"] * math.log(load)))
    return base, offset, log


def interval(row):
    if (not row["ci_low_percent"] or not row["ci_high_percent"]
            or number(row["ci_level_percent"]) != 95):
        return None
    return float(row["ci_low_percent"]) / 100, float(row["ci_high_percent"]) / 100


def expanded_interval(row):
    if interval(row) is None:
        return None
    lower = rounded_bounds(row["ci_low_percent"], 100, 1)[0]
    upper = rounded_bounds(row["ci_high_percent"], 100, 1)[1]
    return lower, upper


def distance(point, bounds):
    return max(bounds[0] - point, 0.0, point - bounds[1])


def range_status(floor, ceiling, bounds):
    if bounds is None:
        return "no_published_95ci"
    lower, upper = bounds
    if upper < floor:
        return "below_floor"
    if lower > ceiling:
        return "above_ceiling"
    if lower == ceiling:
        return "ceiling_boundary_only"
    return "compatible"


def delta_box(base, offset, bounds):
    """Exact offset-minus-baseline squared-error range for fixed predictions."""
    values = [(offset - y) ** 2 - (base - y) ** 2 for y in bounds]
    return min(values), max(values)


def squared_error_box(prediction_bounds, outcome_bounds):
    pl, pu = prediction_bounds
    yl, yu = outcome_bounds
    minimum = max(pl - yu, yl - pu, 0.0) ** 2
    maximum = max((pl - yu) ** 2, (pu - yl) ** 2)
    return minimum, maximum


def median_envelope(point, floor, ceiling):
    return (floor + point) / 2, (point + ceiling) / 2


def overlaps(left, right):
    return max(left[0], right[0]) <= min(left[1], right[1])


def validate_inputs(outcomes, tmb, sources):
    for name, rows in (("outcomes", outcomes), ("tmb", tmb), ("sources", sources)):
        if len({r["id"] for r in rows}) != len(rows):
            raise ValueError("Duplicate IDs in " + name)
        for row in rows:
            for field in PROVENANCE:
                if not row.get(field):
                    raise ValueError(f"Missing {field}: {name}/{row['id']}")
            if row["access_date"] != "2026-11-18":
                raise ValueError("Unexpected access-date snapshot")
            if row.get("nonclinical") != "yes":
                raise ValueError("Missing nonclinical flag")
    outcome_ids = {r["id"] for r in outcomes}
    source_families = {r["source_doi"]: r["family_id"] for r in sources}
    for row in outcomes + tmb:
        if source_families.get(row["source_doi"]) != row["family_id"]:
            raise ValueError("Source and aggregate overlap families disagree: " + row["id"])
    for row in outcomes:
        if row["viral_status"] not in {"positive", "negative", "unknown"}:
            raise ValueError("Unknown status must not be treated as negative")
        if row["role"] not in {"primary", "secondary", "context"}:
            raise ValueError("Unknown report role")
        if row["historical_overlap"] not in {"independent", "possible", "known", "unknown"}:
            raise ValueError("Unknown overlap classification")
        y = number(row["orr_percent"])
        if y is not None and not 0 <= y <= 100:
            raise ValueError("ORR outside percent scale")
        if bool(row["ci_low_percent"]) != bool(row["ci_high_percent"]):
            raise ValueError("Partial confidence interval")
        if row["ci_low_percent"]:
            lower, upper = float(row["ci_low_percent"]), float(row["ci_high_percent"])
            if not (0 <= lower <= upper <= 100 and (y is None or lower <= y <= upper)):
                raise ValueError("Invalid reported interval")
            if not row["ci_level_percent"]:
                raise ValueError("Missing interval level")
        if row["n_orr"]:
            n = int(row["n_orr"])
            if n <= 0:
                raise ValueError("Invalid reported denominator")
            if row["responders"]:
                k = int(row["responders"])
                if not 0 <= k <= n:
                    raise ValueError("Invalid published count")
                if (y is not None and abs(100 * k / n - y) > half_unit(row["orr_percent"]) + 1e-10
                        and not row.get("count_rate_note")):
                    raise ValueError("Published count and rounded ORR disagree")
    links = []
    for row in tmb:
        if row["outcome_id"] and row["outcome_id"] not in outcome_ids:
            raise ValueError("Unknown outcome link")
        if float(row["median"]) < 0:
            raise ValueError("Negative TMB")
        if row["usable"] not in {"yes", "no"}:
            raise ValueError("Unclassified TMB bridge")
        if row["usable"] == "yes":
            if (row["unit"] != "mutations/Mb" or not row["outcome_id"]
                    or row["population_match"] not in {"same", "within_cohort_subset"}):
                raise ValueError("Unsupported median substitution")
            links.append(row["outcome_id"])
    if len(set(links)) != len(links):
        raise ValueError("Multiple selected medians for an outcome")


def evaluate_stratum(row, predictor, scenario, constants):
    positive = row["viral_status"] == "positive"
    floor = probability(scenario["v"], constants) if positive else 0.0
    bounds = interval(row)
    widened = expanded_interval(row)
    result = {"id": row["id"], "scale_label": scenario["scale_label"],
              "target": scenario["target"], "floor": floor,
              "baseline_range": range_status(0.0, constants["A"], bounds),
              "offset_range": range_status(floor, constants["A"], bounds),
              "rounded_offset_range": range_status(floor, constants["A"], widened),
              "rounded_range_gap": None if widened is None else max(
                  floor - widened[1], widened[0] - constants["A"], 0.0),
              "scorable": predictor is not None and row["viral_status"] != "unknown" and bool(row["orr_percent"])
              and row["eligible_monotherapy"] == "yes"}
    if not result["scorable"]:
        return result
    m = float(predictor["median"])
    y = float(row["orr_percent"]) / 100
    base, offset, log = predictions(m, positive, scenario, constants)
    result.update({"baseline": base, "offset": offset, "log": log,
                   "baseline_residual": y - base, "offset_residual": y - offset,
                   "baseline_se": (base - y) ** 2, "offset_se": (offset - y) ** 2,
                   "delta": (offset - y) ** 2 - (base - y) ** 2,
                   "baseline_median_envelope": median_envelope(base, 0.0, constants["A"]),
                   "offset_median_envelope": median_envelope(offset, floor, constants["A"])})
    m_bounds = rounded_bounds(predictor["median"])
    pred_bounds = [predictions(value, positive, scenario, constants) for value in m_bounds]
    base_bounds = pred_bounds[0][0], pred_bounds[1][0]
    offset_bounds = pred_bounds[0][1], pred_bounds[1][1]
    y_bounds = rounded_bounds(row["orr_percent"], 100, 1)
    base_error = squared_error_box(base_bounds, y_bounds)
    offset_error = squared_error_box(offset_bounds, y_bounds)
    result["rounding_delta_bounds"] = ((offset_error[0] - base_error[1],
                                        offset_error[1] - base_error[0])
                                       if positive else (0.0, 0.0))
    result["prediction_bounds"] = {"baseline": base_bounds, "offset": offset_bounds}
    if bounds is not None:
        result.update({"baseline_inside": distance(base, bounds) == 0,
                       "offset_inside": distance(offset, bounds) == 0,
                       "baseline_interval_distance": distance(base, bounds),
                       "offset_interval_distance": distance(offset, bounds),
                       "delta_bounds": delta_box(base, offset, bounds),
                       "rounding_can_include": overlaps(offset_bounds, widened),
                       "baseline_median_compatible": overlaps(
                           result["baseline_median_envelope"], bounds),
                       "offset_median_compatible": overlaps(
                           result["offset_median_envelope"], bounds)})
    return result


def pair_scores(outcomes, strata):
    results = []
    pairs = {}
    for row in outcomes:
        if row["role"] == "primary" and row["pair_id"]:
            pairs.setdefault(row["pair_id"], []).append(row)
    for pair_id, rows in sorted(pairs.items()):
        for field in ("family_id", "cohort_id", "cancer_group"):
            if len({row[field] for row in rows}) != 1:
                raise ValueError("Cross-population pair: " + pair_id)
        statuses = {r["viral_status"] for r in rows}
        record = {"pair_id": pair_id, "family_id": rows[0]["family_id"],
                  "cancer_group": rows[0]["cancer_group"], "point_scorable": False,
                  "full_interval_rule": False, "pass": None,
                  "independent": all(r["historical_overlap"] == "independent" for r in rows)}
        if len(rows) != 2 or statuses != {"positive", "negative"}:
            record["reason"] = "no paired positive/negative ORRs"
        elif not all(strata[r["id"]]["scorable"] for r in rows):
            record["reason"] = "no usable per-Mb median for both statuses"
        else:
            values = [strata[r["id"]] for r in rows]
            record["point_scorable"] = True
            for label in ("baseline", "offset", "log"):
                errors = [v[label] - float(r["orr_percent"]) / 100 for r, v in zip(rows, values)]
                mse = sum(x * x for x in errors) / 2
                record[label] = {"mse": mse, "rmse_pp": math.sqrt(mse) * 100,
                                 "mae_pp": sum(abs(x) for x in errors) * 50}
            record["delta"] = record["offset"]["mse"] - record["baseline"]["mse"]
            record["rounding_delta_bounds"] = tuple(sum(v["rounding_delta_bounds"][i]
                                                        for v in values) / 2 for i in (0, 1))
            if all(interval(r) is not None for r in rows):
                record["full_interval_rule"] = True
                record["pass"] = record["delta"] < 0 and all(v["offset_inside"] for v in values)
                record["delta_bounds"] = tuple(sum(v["delta_bounds"][i] for v in values) / 2
                                                for i in (0, 1))
                for label in ("baseline", "offset"):
                    record[label]["misses"] = sum(not v[label + "_inside"] for v in values)
                    record[label]["interval_mse"] = sum(v[label + "_interval_distance"] ** 2
                                                        for v in values) / 2
                record["reason"] = "pass" if record["pass"] else "failed complete-pair rule"
            else:
                record["reason"] = "point score only; missing published 95% intervals"
        results.append(record)
    return results


def verdict(outcomes, strata, pairs, independent_only=False):
    selected = [r for r in outcomes if r["role"] == "primary"
                and r["cancer_group"] in REQUIRED
                and r["eligible_monotherapy"] == "yes"
                and r["viral_status"] in {"positive", "negative"}
                and (not independent_only or r["historical_overlap"] == "independent")]
    contradictions = [r["id"] for r in selected
                      if strata[r["id"]]["rounded_offset_range"] in {"below_floor", "above_ceiling"}]
    eligible_pairs = [p for p in pairs if p["cancer_group"] in REQUIRED
                      and (not independent_only or p["independent"])]
    failures = [p["pair_id"] for p in eligible_pairs if p["full_interval_rule"] and not p["pass"]]
    if contradictions or failures:
        return {"verdict": "viral offset fails to transport", "range_failures": contradictions,
                "pair_failures": failures}
    passed_groups = {p["cancer_group"] for p in eligible_pairs
                     if p["full_interval_rule"] and p["pass"] and p["independent"]}
    complete_pairs_pass = all(p["full_interval_rule"] and p["pass"]
                              for p in eligible_pairs if p["point_scorable"])
    if passed_groups == REQUIRED and complete_pairs_pass:
        return {"verdict": "viral offset transports", "range_failures": [], "pair_failures": []}
    return {"verdict": "untestable", "range_failures": [], "pair_failures": []}


def balanced_summary(pairs):
    """Empty panels remain unavailable; average within families before across."""
    families = {}
    for pair in pairs:
        if pair["point_scorable"]:
            families.setdefault(pair["family_id"], []).append(pair)
    if not families:
        return None
    family_values = {}
    for name, group in sorted(families.items()):
        family_values[name] = {model: sum(p[model]["mse"] for p in group) / len(group)
                              for model in ("baseline", "offset", "log")}
    summary = {"family_count": len(families), "family_mse": family_values}
    summary["mean_mse"] = {model: sum(v[model] for v in family_values.values()) / len(family_values)
                           for model in ("baseline", "offset", "log")}
    groups = {}
    for family, members in families.items():
        group = members[0]["cancer_group"]
        groups.setdefault(group, []).append(family_values[family])
    summary["cancer_balanced_mse"] = None
    if REQUIRED.issubset(groups):
        summary["cancer_balanced_mse"] = {
            model: sum(sum(v[model] for v in groups[g]) / len(groups[g]) for g in sorted(REQUIRED)) / 2
            for model in ("baseline", "offset", "log")}
    return summary


def analyze():
    constants, scenarios = load_frozen()
    outcomes = read_csv(HERE / "data" / "outcomes.csv")
    tmb = read_csv(HERE / "data" / "tmb.csv")
    sources = read_csv(HERE / "data" / "sources.csv")
    validate_inputs(outcomes, tmb, sources)
    predictors = {r["outcome_id"]: r for r in tmb if r["usable"] == "yes"}
    evaluations = []
    for scenario in scenarios:
        strata = {r["id"]: evaluate_stratum(r, predictors.get(r["id"]), scenario, constants)
                  for r in outcomes}
        pairs = pair_scores(outcomes, strata)
        evaluations.append({"scenario": scenario, "strata": strata, "pairs": pairs,
                            "summary": balanced_summary(pairs),
                            "independent_summary": balanced_summary([p for p in pairs if p["independent"]])})
    primary = next(e for e in evaluations if e["scenario"]["scale_label"] == "1"
                   and e["scenario"]["target"] == "point")
    return {"constants": constants, "outcomes": outcomes, "tmb": tmb, "sources": sources,
            "evaluations": evaluations, "primary": primary,
            "verdict": verdict(outcomes, primary["strata"], primary["pairs"]),
            "independent_verdict": verdict(outcomes, primary["strata"], primary["pairs"], True)}


def fingerprints():
    paths = [ROOT / key for key in PINNED]
    paths.extend(HERE / "data" / name for name in ("outcomes.csv", "tmb.csv", "sources.csv"))
    paths.extend([HERE / "analysis.py", ROOT / "tests" / "test_research_viral_antigen_transport_test.py"])
    return {p.relative_to(ROOT).as_posix(): sha256(p)
            for p in sorted(paths, key=lambda path: path.relative_to(ROOT).as_posix())}


def cell(value):
    return str(value).replace("|", "/").replace("\n", " ")


def table(headers, rows):
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    return "\n".join(lines + ["| " + " | ".join(cell(v) for v in row) + " |" for row in rows])


def percent(value):
    return f"{100 * value:.6f}"


def bounds_text(bounds, factor=1):
    return f"[{factor * bounds[0]:.6f}, {factor * bounds[1]:.6f}]"


def render_report(result):
    c = result["constants"]
    outcomes, tmb = result["outcomes"], result["tmb"]
    primary = result["primary"]
    lines = ["### Verdict and available coverage", "",
             f"**Frozen-rule verdict: {result['verdict']['verdict']}.** "
             f"Strict historical-independence sensitivity: **{result['independent_verdict']['verdict']}**.", ""]
    complete = sum(p["point_scorable"] for p in primary["pairs"])
    full = sum(p["full_interval_rule"] for p in primary["pairs"])
    scorable = sum(v["scorable"] for v in primary["strata"].values())
    lines += [f"The CSVs contain {len(outcomes)} outcome aggregate rows "
              f"({sum(bool(r['orr_percent']) for r in outcomes)} with a response percentage), "
              f"{len(tmb)} published TMB-median rows "
              f"and {len(result['sources'])} source-ledger records. "
              f"There are {complete} complete positive/negative pairs with usable TMB, "
              f"{full} pairs eligible for the full interval rule, and {scorable} individual "
              "strata with a usable ORR–TMB pairing. Rows, follow-ups and overlapping "
              "reports are not independent-study counts.", "",
              "Primary range failures: " + (", ".join(result["verdict"]["range_failures"]) or "none") + ". "
              "Primary complete-pair failures: " + (", ".join(result["verdict"]["pair_failures"]) or "none") + ".", ""]
    lines += ["### Inherited constants and response ranges", "",
              f"A = `{c['A']!r}`, B = `{c['B']!r}`; frozen log-TMB comparator "
              f"a = `{c['a']!r}`, b = `{c['b']!r}`. "
              f"Both reciprocal models have response ceiling {percent(c['A'])}%, "
              "an unattained limit for finite TMB. The baseline lower endpoint is zero.", ""]
    rows = []
    for e in result["evaluations"]:
        s = e["scenario"]
        rows.append([s["scale_label"], s["target"], repr(s["v"]),
                     percent(probability(s["v"], c))])
    lines += [table(["Scale c", "Inherited calibration target", "Offset v", "Viral+ floor %"], rows), "",
              "The CI-low/high calibration targets are inherited training scenarios; "
              "they are not transported predictive intervals.", "",
              "### Published ORR aggregates", "",
              "Blank fields are unavailable or not verified in this aggregate audit; "
              "this does not establish that they were absent from the full source. "
              "no counts are recovered from rounded percentages. Source links and exact "
              "locators are on every CSV row.", ""]
    rows = []
    for r in outcomes:
        ci = (f"[{r['ci_low_percent']}, {r['ci_high_percent']}] ({r['ci_level_percent']}%)"
              if r["ci_low_percent"] else "unavailable")
        rows.append([r["id"], r["role"], r["viral_status"], r["orr_percent"] or "unavailable", ci,
                     r["n_orr"] or "unavailable", r["responders"] or "unavailable"])
    lines += [table(["Outcome ID", "Role", "Viral status", "ORR %", "Published CI", "ORR n", "Responders"], rows), "",
              "Any published count/rate discrepancies are retained, not silently corrected:", ""]
    discrepancies = [[r["id"], r["orr_percent"],
                      f"{100 * int(r['responders']) / int(r['n_orr']):.6f}", r["count_rate_note"]]
                     for r in outcomes if r.get("count_rate_note") and r["responders"] and r["n_orr"]]
    lines += [table(["Outcome", "Printed ORR %", "Printed counts imply %", "Audit note"], discrepancies)
              if discrepancies else "None among rows with both published counts and denominators.", "",
              "### Published TMB medians and eligibility", ""]
    rows = [[r["id"], r["outcome_id"] or "no eligible ORR", r["median"], r["unit"],
             r["n_tmb"] or "unavailable", r["usable"], r["population_match"]] for r in tmb]
    lines += [table(["TMB ID", "Linked ORR", "Median", "Unit", "TMB n", "Usable", "Population match"], rows), "",
              "### Primary model-range checks", "",
              "These checks need no median. Compatibility means only that the published "
              "95% interval intersects the possible response range. It is not a prediction "
              "success. The rounded check expands printed interval endpoints.", ""]
    rows = [[r["id"], primary["strata"][r["id"]]["baseline_range"],
             primary["strata"][r["id"]]["offset_range"],
             primary["strata"][r["id"]]["rounded_offset_range"]] for r in outcomes]
    lines += [table(["Outcome ID", "Baseline", "Offset", "Offset after rounding"], rows), "",
              "### Primary rounding-robust range contradictions", "",
              "The gap is the separation between the widened published interval and "
              "the offset model's closed response range. It is not a paired prediction "
              "error or a p-value. Historical independence is shown separately.", ""]
    failures = set(result["verdict"]["range_failures"])
    rows = []
    for r in outcomes:
        if r["id"] in failures:
            v = primary["strata"][r["id"]]
            rows.append([r["id"], bounds_text(interval(r), 100),
                         bounds_text(expanded_interval(r), 100), percent(v["floor"]),
                         percent(c["A"]), percent(v["rounded_range_gap"]), r["historical_overlap"]])
    lines += [table(["Outcome", "Published CI %", "Widened CI %", "Offset floor %", "Ceiling %", "Gap pp", "Historical overlap"], rows)
              if rows else "None at the primary scenario.", "",
              "### Predictions where a published per-Mb median exists", "",
              "All rows below are single-stratum descriptive calculations. A positive-only "
              "row cannot establish a status-balanced score or a full paired transport pass. "
              "Delta is offset minus baseline squared error in probability-squared units; "
              "negative favors the offset.", ""]
    rows = []
    for e in result["evaluations"]:
        if e["scenario"]["target"] != "point":
            continue
        for r in outcomes:
            v = e["strata"][r["id"]]
            if v["scorable"]:
                rows.append([r["id"], e["scenario"]["scale_label"], percent(v["baseline"]),
                             percent(v["offset"]), percent(v["log"]),
                             percent(v["baseline_residual"]), percent(v["offset_residual"]),
                             f"{v['delta']:.9f}", v.get("baseline_inside", "unavailable"),
                             v.get("offset_inside", "unavailable")])
    lines += [table(["Outcome", "c", "Baseline %", "Offset %", "Log %", "Base residual pp", "Offset residual pp", "Delta SE", "Base in CI", "Offset in CI"], rows)
              if rows else "Unavailable: no eligible ORR–TMB pair.", "",
              "### Paired scoring availability", ""]
    lines += [table(["Pair", "Family", "Point score", "Full interval rule", "Reason"],
                    [[p["pair_id"], p["family_id"], p["point_scorable"], p["full_interval_rule"], p["reason"]]
                     for p in primary["pairs"]]), ""]
    if primary["summary"] is None:
        lines += ["Status-balanced MSE, RMSE, MAE, interval-distance MSE and paired "
                  "uncertainty boxes are **unavailable**. The empty comparison is not "
                  "assigned a zero loss. No cohort- or cancer-balanced comparison is estimated.", ""]
    else:
        rows = []
        for e in result["evaluations"]:
            for p in e["pairs"]:
                if not p["point_scorable"]:
                    continue
                for model in ("baseline", "offset", "log"):
                    m = p[model]
                    rows.append([p["pair_id"], e["scenario"]["scale_label"], e["scenario"]["target"], model,
                                 f"{m['mse']:.9f}", f"{m['rmse_pp']:.6f}", f"{m['mae_pp']:.6f}",
                                 m.get("misses", "unavailable"), f"{m['interval_mse']:.9f}" if "interval_mse" in m else "unavailable"])
        lines += [table(["Pair", "c", "Target", "Model", "MSE", "RMSE pp", "MAE pp", "CI misses", "Interval MSE"], rows), "",
                  "Primary aggregate scores (family first; null means unavailable):", "",
                  "```json", json.dumps(primary["summary"], sort_keys=True, indent=2), "```", "",
                  "Strict historical-independence score summary:", "",
                  "```json", json.dumps(primary["independent_summary"], sort_keys=True, indent=2), "```", ""]
    lines += ["### Inherited uncertainty, rounding and median-substitution checks", "",
              "Every inherited offset is evaluated. The rounding delta bounds below are "
              "conservative interval-arithmetic bounds, allowing printed ORR and median "
              "rounding. They are not sampling confidence intervals. A missing published "
              "ORR interval prevents both interval-box and median-envelope compatibility tests.", ""]
    rows = []
    for e in result["evaluations"]:
        for r in outcomes:
            v = e["strata"][r["id"]]
            if v["scorable"]:
                rows.append([r["id"], e["scenario"]["scale_label"], e["scenario"]["target"],
                             percent(v["offset"]), f"{v['delta']:.9f}", bounds_text(v["rounding_delta_bounds"]),
                             bounds_text(v["delta_bounds"]) if "delta_bounds" in v else "unavailable",
                             v.get("rounding_can_include", "unavailable"),
                             v.get("baseline_median_compatible", "unavailable"),
                             v.get("offset_median_compatible", "unavailable")])
    lines += [table(["Outcome", "c", "Target", "Offset %", "Delta SE", "Rounding delta bounds", "95% ORR-box delta bounds", "Rounding permits CI inclusion", "Base median compatible", "Offset median compatible"], rows)
              if rows else "Unavailable: no eligible ORR–TMB pair.", "",
              "### Population-median envelope endpoints", "",
              "These closure bounds on mean model predictions use the published median "
              "as a population-median scenario. They remain calculable without a published "
              "ORR interval, but they are not confidence intervals; empirical interval "
              "compatibility stays unavailable when that ORR interval is missing.", ""]
    rows = []
    for e in result["evaluations"]:
        for r in outcomes:
            v = e["strata"][r["id"]]
            if v["scorable"]:
                rows.append([r["id"], e["scenario"]["scale_label"], e["scenario"]["target"],
                             bounds_text(v["baseline_median_envelope"], 100),
                             bounds_text(v["offset_median_envelope"], 100)])
    lines += [table(["Outcome", "c", "Target", "Baseline envelope %", "Offset envelope %"], rows)
              if rows else "Unavailable: no eligible ORR–TMB pair.", "",
              "### Range sensitivity across all inherited offsets", "",
              "Lists include primary, secondary and context rows for transparency; only "
              "selected required-cancer primary rows at the identity point scenario "
              "can decide the headline verdict.", ""]
    rows = []
    for e in result["evaluations"]:
        failures = [key for key, v in e["strata"].items()
                    if v["rounded_offset_range"] in {"below_floor", "above_ceiling"}]
        rows.append([e["scenario"]["scale_label"], e["scenario"]["target"], ", ".join(failures) or "none"])
    lines += [table(["c", "Target", "Rounding-robust interval/range contradictions"], rows), "",
              "### Canonical SHA-256 fingerprints", "",
              "UTF-8 bytes with CRLF/CR normalized to LF; keys are built with "
              "`Path.as_posix()` and sorted as canonical POSIX strings.", "",
              table(["Path", "SHA-256"], fingerprints().items())]
    return "\n".join(lines).strip() + "\n"


def report_text(result):
    path = HERE / "README.md"
    text = path.read_text(encoding="utf-8")
    before, remainder = text.split(BEGIN, 1)
    _, after = remainder.split(END, 1)
    return before + BEGIN + "\n\n" + render_report(result) + "\n" + END + after


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    result = analyze()
    if args.json:
        print(json.dumps(result, sort_keys=True, indent=2))
    elif args.write_report:
        updated_text = report_text(result)
        with (HERE / "README.md").open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(updated_text)
    elif args.check:
        if (HERE / "README.md").read_text(encoding="utf-8") != report_text(result):
            raise SystemExit("README numerical report or fingerprints are stale")
        print("NONCLINICAL: frozen inputs, source schemas and generated report agree")
    else:
        print(render_report(result), end="")


if __name__ == "__main__":
    main()
