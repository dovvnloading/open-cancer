"""Recompute the frozen, nonclinical CH age-prevalence audit. No medical advice.

No network, third-party package, target-driven fitness fitting, or source download.
Run without options to verify the committed report; --write regenerates outputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist, mean


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MODEL_SPEC = importlib.util.spec_from_file_location("ch_age_frozen_model", HERE / "model.py")
model = importlib.util.module_from_spec(MODEL_SPEC)
MODEL_SPEC.loader.exec_module(model)
START = "<!-- BEGIN RECOMPUTED RESULTS -->"
END = "<!-- END RECOMPUTED RESULTS -->"
# The integration target is 1e-10 relative with a small absolute tolerance.
# This comparison budget covers binary64/libm and quadrature stopping variation
# across platforms. It is numerical verification tolerance, not uncertainty.
# Integers, metadata, keys, hashes, and the displayed README block remain exact.
REPORT_FLOAT_REL_TOL = 5e-9
REPORT_FLOAT_ABS_TOL = 5e-10


def normalized_bytes(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def digest(path: Path) -> str:
    return hashlib.sha256(normalized_bytes(path)).hexdigest()


def read_csv(name: str) -> list[dict[str, str]]:
    with (HERE / "data" / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def optional_number(row: dict, key: str) -> float | None:
    value = row.get(key, "")
    if value == "":
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"nonfinite source value: {key}")
    return result


def wilson_interval(cases: int, total: int) -> tuple[float, float]:
    """Marginal approximate 95% interval; not a simultaneous model test."""
    if not 0 <= cases <= total or total <= 0:
        raise ValueError("invalid aggregate counts")
    z = NormalDist().inv_cdf(0.975)
    p = cases / total
    denominator = 1 + z*z/total
    centre = (p + z*z/(2*total))/denominator
    half = z*math.sqrt(p*(1-p)/total + z*z/(4*total*total))/denominator
    return max(0.0, centre-half), min(1.0, centre+half)


def observed(row: dict) -> dict:
    k = optional_number(row, "cases")
    n = optional_number(row, "total")
    percent = optional_number(row, "prevalence_percent")
    if k is not None and n is None:
        raise ValueError(f"{row['record_id']}: count requires its own denominator")
    if n is not None:
        if int(n) != n or n <= 0:
            raise ValueError("denominator must be a positive integer")
        n = int(n)
    if k is not None:
        if int(k) != k or not 0 <= k <= n:
            raise ValueError("invalid aggregate count")
        k = int(k)
        point = k/n
    elif percent is not None:
        point = percent/100
    else:
        raise ValueError("missing prevalence")
    lower = optional_number(row, "ci_lower_percent")
    upper = optional_number(row, "ci_upper_percent")
    if (lower is None) != (upper is None):
        raise ValueError("partial printed confidence interval")
    if lower is not None:
        interval = (lower/100, upper/100)
        kind = "published marginal 95% confidence interval"
    elif k is not None:
        interval = wilson_interval(int(k), int(n))
        kind = "computed marginal Wilson 95% interval"
    else:
        radius = 0.5 * 10**(-int(row["percentage_decimals"])) / 100
        interval = (max(0.0, point-radius), min(1.0, point+radius))
        kind = "printed rounding interval; sampling uncertainty unavailable"
    if not 0 <= interval[0] <= point <= interval[1] <= 1:
        raise ValueError(f"invalid observation interval: {row['record_id']}")
    conflict = False
    alternative = None
    if percent is not None and k is not None:
        radius = 0.5 * 10**(-int(row["percentage_decimals"])) / 100
        alternative = percent/100
        conflict = not alternative-radius-1e-15 <= point <= alternative+radius+1e-15
        if conflict:
            # Transparent conservative source-conflict amendment, independent of fit.
            interval = (min(interval[0], max(0.0, alternative-radius)),
                        max(interval[1], min(1.0, alternative+radius)))
            kind += "; hull with conflicting printed percentage rounding interval"
    return dict(point=point, lower=interval[0], upper=interval[1],
                uncertainty=kind, cases=k, total=n, source_conflict=conflict,
                printed_percentage_alternative=alternative)


def binomial_deviance(cases: int, total: int, predicted: float) -> float:
    if not 0 < predicted < 1:
        if (predicted == 0 and cases == 0) or (predicted == 1 and cases == total):
            return 0.0
        return math.inf
    terms = []
    if cases:
        terms.append(cases * math.log(cases / (total*predicted)))
    if cases < total:
        terms.append((total-cases) * math.log((total-cases)/(total*(1-predicted))))
    return 2*math.fsum(terms)


def predict_row(row: dict, params: dict, vaf: float,
                n_tau: float | None = None) -> dict:
    a, b = optional_number(row, "age_lower"), optional_number(row, "age_upper")
    bound_a = a if a is not None else 0.0
    if a is not None and b is not None:
        low, point, high = model.band_prediction(params, a, b, vaf, n_tau)
        midpoint = (a+b)/2
        intensity = model.expected_clones(params, midpoint, vaf, n_tau)
    else:
        low = model.prevalence(params, bound_a, vaf, n_tau)
        high = model.prevalence(params, b, vaf, n_tau) if b is not None else 1.0
        point, midpoint, intensity = None, None, None
    la = model.expected_clones(params, bound_a, vaf, n_tau)
    lb = model.expected_clones(params, b, vaf, n_tau) if b is not None else math.inf
    return dict(lower=low, uniform_age_mean=point, upper=high,
                midpoint_age=midpoint, midpoint_intensity=intensity,
                lower_intensity=la, upper_intensity=lb, vaf=vaf)


def compare(obs: dict, prediction: dict) -> dict:
    low, high = prediction["lower"], prediction["upper"]
    point = prediction["uniform_age_mean"]
    gap = max(0.0, low-obs["upper"], obs["lower"]-high)
    residual = 100*(obs["point"]-point) if point is not None else None
    deviance = (binomial_deviance(int(obs["cases"]), int(obs["total"]), point)
                if obs["cases"] is not None and point is not None else None)
    return dict(signed_error_pp=residual, age_interval_gap_pp=100*gap,
                interval_compatible=gap <= 1e-12, binomial_deviance=deviance)


def supply_interval(prediction: dict, obs: dict) -> tuple[float, float]:
    """Use the shared finite-band theorem, or the stated open-age closure.

    With an unbounded upper age, the frozen conservative upper prediction is
    one. A zero lower multiplier is then an infimum when observed_lower>0,
    rather than an attainable fit at c=0. Such rows are labelled as closures.
    """
    la, lb = prediction["lower_intensity"], prediction["upper_intensity"]
    if lb == math.inf:
        # Exact probability one cannot occur at finite age/intensity and c.
        if obs["lower"] == 1:
            return math.inf, -math.inf
        _, upper = model.multiplier_interval(la, max(la, 1.0), 0.0, obs["upper"])
        return 0.0, upper
    return model.multiplier_interval(la, lb, obs["lower"], obs["upper"])


def score_rows(rows: list[dict], model_parameters: dict, cohorts: dict) -> list[dict]:
    parent_statuses = defaultdict(set)
    for cohort in cohorts.values():
        parent_statuses[cohort["parent_cohort"]].add(cohort["independence"])
    scored = []
    for row in rows:
        cohort = cohorts[row["cohort_id"]]
        endpoint_id = row["endpoint_id"]
        if not endpoint_id:
            raise ValueError(f"{row['record_id']}: endpoint_id is required")
        model_id = row["model_id"]
        obs = observed(row)
        params = model_parameters[model_id]
        vaf = float(cohort["source_vaf"]) if cohort["source_vaf"] else 0.02
        reference = predict_row(row, params, 0.02)
        source = predict_row(row, params, vaf)
        parent_independent = parent_statuses[cohort["parent_cohort"]] == {"independent"}
        finite_age = (optional_number(row, "age_lower") is not None
                      and optional_number(row, "age_upper") is not None)
        strict = (row["analysis_role"] in ("external_age", "external_gene_age")
                  and finite_age
                  and row["endpoint_match"] == "yes" and cohort["assay_calibrated"] == "yes"
                  and parent_independent and not obs["source_conflict"])
        scored.append(dict(record_id=row["record_id"], cohort_id=row["cohort_id"],
                           parent_cohort=cohort["parent_cohort"], gene=row["gene"],
                           model_id=model_id, age_label=row["age_label"],
                           source_status=row.get("source_status", ""),
                           threshold_basis=("printed source cutoff" if cohort["source_vaf"]
                                            else "2% reference only; source has no uniform cutoff"),
                           endpoint_id=endpoint_id, endpoint_match=row["endpoint_match"],
                           independence=cohort["independence"],
                           parent_independence_status=("independent" if parent_independent
                                                       else "not established independent"),
                           analysis_role=row["analysis_role"], role=row["analysis_role"],
                           strict_eligible=strict,
                           observation=obs, reference_prediction=reference,
                           reference_comparison=compare(obs, reference),
                           source_prediction=source, source_comparison=compare(obs, source),
                           supply_bound_kind=("sharp finite-band bounds" if finite_age
                                              else "conservative closure for open or missing age support"),
                           supply_feasible=supply_interval(source, obs)))
    return scored


def summary(rows: list[dict], prediction_key: str = "source_comparison") -> dict:
    by_parent = defaultdict(list)
    for row in rows:
        residual = row[prediction_key]["signed_error_pp"]
        if residual is not None:
            by_parent[row["parent_cohort"]].append(residual)
    mae = mean(mean(abs(e) for e in errors) for errors in by_parent.values()) if by_parent else None
    mse = mean(mean(e*e for e in errors) for errors in by_parent.values()) if by_parent else None
    parents = sorted({r["parent_cohort"] for r in rows})
    return dict(rows=len(rows), point_rows=sum(map(len, by_parent.values())),
                parents=len(parents), parent_ids=parents,
                point_parents=len(by_parent), point_parent_ids=sorted(by_parent), macro_mae_pp=mae,
                macro_rmse_pp=math.sqrt(mse) if mse is not None else None,
                nonoverlap_rows=sum(not row[prediction_key]["interval_compatible"] for row in rows))


def loco(rows: list[dict], held_parent: str) -> dict:
    """Secondary calibration, excluding EVERY row from the held parent."""
    finite_parents = {r["parent_cohort"] for r in rows
                      if r["source_prediction"]["midpoint_intensity"] is not None}
    independence = ("independent parents as audited"
                    if all(r.get("parent_independence_status", r.get("independence")) == "independent"
                           for r in rows)
                    else "non-independent interpretation: unresolved or linked cohort overlap")
    common = dict(held_parent=held_parent, independence_status=independence,
                  finite_parent_ids=sorted(finite_parents))
    if len(finite_parents) < 2:
        return dict(**common, status="untestable: fewer than two external parents with finite-band rows")
    if held_parent not in finite_parents:
        return dict(**common, status="untestable: held parent has no finite-band rows")
    grouped = defaultdict(list)
    exclusions = []
    for row in rows:
        o, p = row["observation"], row["source_prediction"]
        if row["parent_cohort"] == held_parent:
            continue
        if not 0 < o["point"] < 1:
            reason = "endpoint probability 0 or 1 excluded from log transform; no pseudocount"
        elif p["midpoint_intensity"] is None:
            reason = "no finite-band midpoint"
        elif p["midpoint_intensity"] <= 0:
            reason = "nonpositive midpoint intensity excluded from log transform"
        else:
            grouped[row["parent_cohort"]].append(row)
            continue
        exclusions.append(dict(record_id=row["record_id"], reason=reason))
    if not grouped:
        return dict(**common, status="untestable: no usable training parent", training_exclusions=exclusions)
    log_c = mean(mean(math.log(-math.log1p(-r["observation"]["point"]))
                     - math.log(r["source_prediction"]["midpoint_intensity"])
                     for r in group) for group in grouped.values())
    c = math.exp(log_c)
    baseline = mean(mean(r["observation"]["point"] for r in group) for group in grouped.values())
    held = []
    for row in rows:
        intensity = row["source_prediction"]["midpoint_intensity"]
        if row["parent_cohort"] != held_parent or intensity is None:
            continue
        p = -math.expm1(-c*intensity)
        held.append(dict(record_id=row["record_id"], prediction=p,
                         error_pp=100*(row["observation"]["point"]-p),
                         constant_prediction=baseline,
                         constant_error_pp=100*(row["observation"]["point"]-baseline)))
    return dict(**common, status="descriptive; mutation supply refitted on other parents",
                training_parents=sorted(grouped),
                training_ids=sorted(r["record_id"] for g in grouped.values() for r in g),
                training_exclusions=exclusions,
                multiplier=c, constant_baseline=baseline, predictions=held,
                mae_pp=mean(abs(x["error_pp"]) for x in held) if held else None,
                rmse_pp=math.sqrt(mean(x["error_pp"]**2 for x in held)) if held else None,
                constant_mae_pp=mean(abs(x["constant_error_pp"]) for x in held) if held else None,
                constant_rmse_pp=math.sqrt(mean(x["constant_error_pp"]**2 for x in held)) if held else None)


def comparison_groups(rows: list[dict]) -> dict[tuple[str, str, str, str], list[dict]]:
    """Keep alternative outcomes and whole-cohort gene marginals separate."""
    groups = defaultdict(list)
    for row in rows:
        groups[(row["cohort_id"], row["model_id"], row["endpoint_id"], row["analysis_role"])].append(row)
    return groups


def feasibility(rows: list[dict]) -> list[dict]:
    result = []
    for (cohort, model_id, endpoint_id, role), group in sorted(comparison_groups(rows).items()):
        low = max(r["supply_feasible"][0] for r in group)
        high = min(r["supply_feasible"][1] for r in group)
        result.append(dict(cohort_id=cohort, model_id=model_id, endpoint_id=endpoint_id,
                           analysis_role=role, lower=low, upper=high,
                           compatible=low <= high, includes_frozen_supply=low <= 1 <= high,
                           bound_kinds=sorted({r["supply_bound_kind"] for r in group}),
                           row_ids=[r["record_id"] for r in group]))
    return result


def external_strata(rows: list[dict]) -> dict[str, list[dict]]:
    """Uncertain status on any linked release applies to its whole parent."""
    independent = {r["parent_cohort"] for r in rows
                   if r["parent_independence_status"] == "independent"}
    uncertain = {r["parent_cohort"] for r in rows
                 if r["parent_independence_status"] != "independent"}
    independent -= uncertain
    return {"independent": [r for r in rows if r["parent_cohort"] in independent],
            "overlap_uncertain": [r for r in rows if r["parent_cohort"] not in independent]}


def loco_collection(rows: list[dict]) -> dict:
    parents = sorted({r["parent_cohort"] for r in rows
                      if r["source_prediction"]["midpoint_intensity"] is not None})
    return dict(status=("descriptive secondary calibration" if len(parents) >= 2
                        else "untestable: fewer than two parents with finite-band rows"),
                finite_parents=len(parents), finite_parent_ids=parents,
                results=[loco(rows, p) for p in parents] if len(parents) >= 2 else [])


def numeric_clean(value):
    """Stable serialized precision; unbounded intervals are explicit strings."""
    if isinstance(value, float):
        if math.isnan(value):
            raise ValueError("NaN in numerical report")
        if not math.isfinite(value):
            return "infinity" if value > 0 else "minus_infinity"
        return round(value, 12)
    if isinstance(value, dict):
        return {key: numeric_clean(v) for key, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [numeric_clean(v) for v in value]
    return value


def build_report() -> dict:
    rows = read_csv("prevalence.csv")
    cohorts = {r["cohort_id"]: r for r in read_csv("cohorts.csv")}
    parameters = model.load_parameters()
    scored = score_rows(rows, parameters, cohorts)
    external = [r for r in scored if r["role"] == "external_age"]
    development = [r for r in scored if r["role"] == "development"]
    marginal_gene = [r for r in scored if r["role"] == "external_gene_marginal"]
    external_gene_age = [r for r in scored if r["role"] == "external_gene_age"]
    strata = external_strata(external)
    overall_loco = loco_collection(external)
    stratum_loco = {name: loco_collection(group) for name, group in strata.items()}
    strict = [r for r in scored if r["strict_eligible"]]
    strict_loco = loco_collection(strict)
    strict_parents = {r["parent_cohort"] for r in strict}
    strict_has_genes = any(r["model_id"] != "POOLED10" and r["role"] == "external_gene_age"
                           for r in strict)
    if len(strict_parents) < 2 or not strict_has_genes:
        verdict = "untestable"
    elif any(not r["source_comparison"]["interval_compatible"] for r in strict):
        verdict = "fails"
    else:
        verdict = "fitness predicts age-prevalence (compatibility only)"
    cohort_scores = []
    for (cohort_id, model_id, endpoint_id, role), subset in sorted(comparison_groups(scored).items()):
        cohort_scores.append(dict(cohort_id=cohort_id, model_id=model_id,
                                  endpoint_id=endpoint_id, analysis_role=role, **summary(subset)))
    deletions = []
    for parent in sorted({r["parent_cohort"] for r in external}):
        deletions.append(dict(removed_parent=parent,
                             **summary([r for r in external if r["parent_cohort"] != parent])))
    sensitivities = []
    for scenario in read_csv("scenarios.csv"):
        selected = []
        for row in rows:
            if row["analysis_role"] != "external_age":
                continue
            obs = observed(row)
            prediction = predict_row(row, parameters[row["model_id"]],
                                     float(scenario["vaf"]), float(scenario["n_tau"]))
            selected.append(dict(parent_cohort=cohorts[row["cohort_id"]]["parent_cohort"],
                                 source_comparison=compare(obs, prediction)))
        sensitivities.append(dict(scenario_id=scenario["scenario_id"], **summary(selected)))
    detection = []
    for c in cohorts.values():
        if c["depth_for_illustration"] and c["min_alt_reads"]:
            d, r = int(c["depth_for_illustration"]), int(c["min_alt_reads"])
            detection.append(dict(cohort_id=c["cohort_id"], depth=d, min_alt_reads=r,
                                  at_vaf_02=model.binomial_detection(d, .02, r),
                                  at_vaf_05=model.binomial_detection(d, .05, r),
                                  at_vaf_10=model.binomial_detection(d, .10, r),
                                  interpretation=c["depth_illustration_basis"]))
    future = []
    for grid in read_csv("prediction_grid.csv"):
        params = parameters[grid["model_id"]]
        t, v = float(grid["age"]), float(grid["vaf"])
        p = model.prevalence(params, t, v)
        cutoff_p = model.prevalence(params, t, v, v_upper=.4999)
        future.append(dict(model_id=grid["model_id"], age=t, vaf=v, prediction=p,
                           source_cutoff_prediction=cutoff_p, cutoff_difference_pp=100*(p-cutoff_p)))
    pooled = {(x["age"], x["vaf"]): x for x in future if x["model_id"] == "POOLED10"}
    subset_violations = []
    subset_comparisons = 0
    for x in future:
        key = (x["age"], x["vaf"])
        if x["model_id"] == "POOLED10" or key not in pooled:
            continue
        subset_comparisons += 1
        other = pooled[key]["prediction"]
        if x["prediction"] > other:
            subset_violations.append(dict(gene=x["model_id"], age=x["age"], vaf=x["vaf"],
                                          gene_prediction=x["prediction"], pooled_prediction=other,
                                          gap_pp=100*(x["prediction"]-other)))
    comparisons = [(abs(r["source_comparison"]["signed_error_pp"]), r["record_id"])
                   for r in external if r["source_comparison"]["signed_error_pp"] is not None]
    report = dict(
        status="nonclinical; unreviewed; no medical advice", strict_verdict=verdict,
        strict_eligible_rows=len(strict), strict_eligible_parents=len(strict_parents),
        external_substitute_summary=summary(external),
        external_reference_02_summary=summary(external, "reference_comparison"),
        external_independent_summary=summary(strata["independent"]),
        external_overlap_uncertain_summary=summary(strata["overlap_uncertain"]),
        external_independence_note=("A parent is in the independent stratum only if every audited "
                                    "linked release is labelled independent. The other stratum "
                                    "retains unresolved or known overlap; the combined summary "
                                    "is descriptive, not independent validation."),
        strict_summary=summary(strict), external_gene_age_summary=summary(external_gene_age),
        development_summary=summary(development), marginal_gene_summary=summary(marginal_gene),
        descriptive_verdict=("untestable: no external age rows" if not external else
                             "fails frozen interval-compatibility criterion"
                             if any(not r["source_comparison"]["interval_compatible"] for r in external)
                             else "compatible at printed interval/age resolution"),
        scored_rows=scored, cohort_scores=cohort_scores, delete_one_parent=deletions,
        loco_supply_calibration=overall_loco["results"],
        loco_supply_calibration_status=overall_loco["status"],
        loco_supply_calibration_strata=stratum_loco,
        strict_loco_status=strict_loco["status"],
        strict_loco_supply_calibration=strict_loco["results"],
        supply_feasibility=feasibility(scored), sensitivities=sensitivities,
        finite_age_supply_feasibility=feasibility([
            r for r in external if r["source_prediction"]["midpoint_intensity"] is not None]),
        detection_illustrations=detection, future_predictions=future,
        subset_coherence=dict(comparisons=subset_comparisons,
                              violation_count=len(subset_violations), violations=subset_violations,
                              maximum_gap_pp=max((x["gap_pp"] for x in subset_violations), default=0.0),
                              interpretation="separate fitted candidate forecasts; violations preclude a joint carrier model with these subset endpoints"),
        largest_external_age_errors=[record_id for _, record_id in sorted(comparisons, reverse=True)],
        largest_external_age_gaps=[r["record_id"] for r in sorted(
            external, key=lambda r: (-r["source_comparison"]["age_interval_gap_pp"], r["record_id"]))],
        largest_gene_marginal_gaps=[r["record_id"] for r in sorted(
            marginal_gene, key=lambda r: (-r["source_comparison"]["age_interval_gap_pp"], r["record_id"]))],
        normalized_input_sha256={p.relative_to(HERE).as_posix(): digest(p)
                                 for p in sorted((HERE / "data").glob("*.csv"))},
        report_verification=dict(float_relative_tolerance=REPORT_FLOAT_REL_TOL,
                                 float_absolute_tolerance=REPORT_FLOAT_ABS_TOL,
                                 exact="integer values, nonnumeric values, keys, hashes and displayed README",
                                 basis="binary64/libm and numerical-integration variation; not uncertainty"),
        frozen_protocol_sha256=digest(HERE / "PROTOCOL.md"))
    return numeric_clean(report)


def fmt(value, digits: int = 4) -> str:
    return "unavailable" if value is None else f"{value:.{digits}f}" if isinstance(value, (int, float)) else str(value)


def pct(value) -> str:
    return fmt(100*value) if isinstance(value, (int, float)) else fmt(value)


def render_results(r: dict) -> str:
    s = r["external_substitute_summary"]
    lines = [START, "", f"**Strict biological verdict: {r['strict_verdict'].upper()}.**",
             f"Eligible matched rows: **{r['strict_eligible_rows']}**; independent eligible parents: "
             f"**{r['strict_eligible_parents']}**. The descriptive substituted-endpoint result "
             f"**{r['descriptive_verdict']}**.", "",
             f"External substituted endpoints: **{s['rows']} rows**, **{s['point_rows']} finite-band "
             f"point scores**, **{s['parents']} represented parent cohorts**, of which "
             f"**{s['point_parents']} contribute point scores**. With each source's printed threshold "
             f"(or the explicitly labelled reference when missing), equal-parent MAE is "
             f"**{fmt(s['macro_mae_pp'])} percentage points**, RMSE **{fmt(s['macro_rmse_pp'])} pp**; "
             f"**{s['nonoverlap_rows']} rows** have nonoverlapping observation and age-prediction intervals.", "",
             "These errors use an assumed uniform age distribution within finite bands. The intervals "
             "below retain unknown age composition. Rounding intervals do not supply sampling uncertainty. "
             "Model sensitivities do not supply fitted-parameter confidence intervals.", "",
             "### External independence strata", "", r["external_independence_note"], "",
             "| Independence stratum | Rows | Represented parents | Point-score parents | MAE pp | RMSE pp | Nonoverlap rows |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for name, key in (("Independent as audited", "external_independent_summary"),
                      ("Overlap unresolved or linked", "external_overlap_uncertain_summary")):
        x = r[key]
        lines.append(f"| {name} | {x['rows']} | {x['parents']} | {x['point_parents']} | "
                     f"{fmt(x['macro_mae_pp'])} | {fmt(x['macro_rmse_pp'])} | {x['nonoverlap_rows']} |")
    lines += ["",
             "### Complete age-band comparison", "",
             "Predictions use the pooled ten-gene fit for any-CH endpoints. These endpoints are "
             "substitutions for Watson's variant panel. Error = observed minus uniform-age prediction.", "",
             "Rows marked reference use 2% because a uniform source cutoff is unavailable. "
             "Published/Wilson labels denote marginal 95% intervals; rounding-only rows have no "
             "computed sampling interval. Every interval is conditional on reported-age support.", "",
             "| Row | VAF % / basis | Observed % | Observation interval % / kind | Prediction %: age range; uniform mean | Error pp | Gap pp |",
             "|---|---|---:|---|---:|---:|---:|"]
    for row in r["scored_rows"]:
        if row["role"] == "external_gene_marginal":
            continue
        o, p, c = row["observation"], row["source_prediction"], row["source_comparison"]
        kind = ("published 95%" if o["uncertainty"].startswith("published") else
                "Wilson 95%" if o["cases"] is not None else "rounding only")
        basis = "source" if row["threshold_basis"] == "printed source cutoff" else "reference only"
        lines.append(f"| {row['record_id']} ({row['age_label']}) | {100*p['vaf']:g} / {basis} | {pct(o['point'])} | "
                     f"{pct(o['lower'])}–{pct(o['upper'])} / {kind} | {pct(p['lower'])}–{pct(p['upper'])}; "
                     f"{pct(p['uniform_age_mean'])} | {fmt(c['signed_error_pp'])} | {fmt(c['age_interval_gap_pp'])} |")
    lines += ["", "### Cohort scores", "",
              "Alternative endpoints, model panels and whole-cohort gene marginals remain separate.", "",
              "| Cohort / model / endpoint / role | Point rows | MAE pp | RMSE pp | Nonoverlap rows |",
              "|---|---:|---:|---:|---:|"]
    for x in r["cohort_scores"]:
        label = f"{x['cohort_id']} / {x['model_id']} / {x['endpoint_id']} / {x['analysis_role']}"
        lines.append(f"| {label} | {x['point_rows']} | {fmt(x['macro_mae_pp'])} | "
                     f"{fmt(x['macro_rmse_pp'])} | {x['nonoverlap_rows']} |")
    lines += ["", "### Leave-one-parent-out checks", "",
              "With fixed external inputs, deleting a cohort changes the aggregate error, never the "
              "remaining predictions. This is an influence check.", "",
              "| Parent removed | Remaining point rows | MAE pp | RMSE pp |",
              "|---|---:|---:|---:|"]
    for x in r["delete_one_parent"]:
        lines.append(f"| {x['removed_parent']} | {x['point_rows']} | {fmt(x['macro_mae_pp'])} | {fmt(x['macro_rmse_pp'])} |")
    lines += ["", "The secondary calculation changes only the common mutation-supply multiplier using "
              "other parents. All rows sharing the held-out parent are excluded from calibration. "
              "This is a descriptive transport check across unmatched endpoints, not primary validation. "
              "The combined calculation retains a non-independent interpretation whenever any parent "
              "has unresolved or linked overlap.", "",
              f"Combined-calculation status: **{r['loco_supply_calibration_status']}**.", "",
              "| Held-out parent | Independence interpretation | Training-only supply multiplier | Held-out MAE pp | Held-out RMSE pp | Constant-baseline MAE pp | Constant-baseline RMSE pp |",
              "|---|---|---:|---:|---:|---:|---:|"]
    for x in r["loco_supply_calibration"]:
        lines.append(f"| {x['held_parent']} | {x['independence_status']} | {fmt(x.get('multiplier'),6)} | {fmt(x.get('mae_pp'))} | "
                     f"{fmt(x.get('rmse_pp'))} | {fmt(x.get('constant_mae_pp'))} | {fmt(x.get('constant_rmse_pp'))} |")
    lines += ["", "The same calibration is also computed separately within each independence stratum; "
              "a stratum with fewer than two finite-band parents has no LOCO fit.", "",
              "| Independence stratum | Finite-band parents | Status |",
              "|---|---:|---|"]
    for name, collection in r["loco_supply_calibration_strata"].items():
        lines.append(f"| {name} | {collection['finite_parents']} | {collection['status']} |")
    lines += ["", "| Independence stratum / held-out parent | Supply multiplier | MAE pp | RMSE pp | Constant MAE pp | Constant RMSE pp |",
              "|---|---:|---:|---:|---:|---:|"]
    for name, collection in r["loco_supply_calibration_strata"].items():
        for x in collection["results"]:
            lines.append(f"| {name} / {x['held_parent']} | {fmt(x.get('multiplier'),6)} | "
                         f"{fmt(x.get('mae_pp'))} | {fmt(x.get('rmse_pp'))} | "
                         f"{fmt(x.get('constant_mae_pp'))} | {fmt(x.get('constant_rmse_pp'))} |")
    lines += ["", "### Common-supply compatibility", "",
              "The intersection uses every admitted age interval within the same endpoint, model and role. "
              "Finite-band bounds are sharp under their stated assumptions. Open or missing age support "
              "instead gives a conservative closure; a zero lower multiplier can be an unattained infimum. "
              "An empty intersection rules out a common multiplier under the stated kernel/measurement "
              "assumptions. These interval boxes are not a simultaneous confidence test.", "",
              "| Cohort / model / endpoint / role | Lower multiplier | Upper multiplier | Nonempty? | Contains frozen multiplier 1? |",
              "|---|---:|---:|---|---|"]
    for x in r["supply_feasibility"]:
        label = f"{x['cohort_id']} / {x['model_id']} / {x['endpoint_id']} / {x['analysis_role']}"
        lines.append(f"| {label} | {fmt(x['lower'],6)} | {fmt(x['upper'],6)} | "
                     f"{x['compatible']} | {x['includes_frozen_supply']} |")
    lines += ["", "A post-score robustness check repeats the same intersection using only finite "
              "external age bands. It changes no model input or primary score and cannot be used "
              "to select which bands to retain.", "",
              "| Finite-band cohort | Lower multiplier | Upper multiplier | Nonempty? |",
              "|---|---:|---:|---|"]
    for x in r["finite_age_supply_feasibility"]:
        lines.append(f"| {x['cohort_id']} | {fmt(x['lower'],6)} | {fmt(x['upper'],6)} | {x['compatible']} |")
    lines += ["", "### Broad gene marginals", "",
              "These are whole-cohort gene counts or printed percentages over their full age support, not gene-by-decade curves. "
              "Their uniform-age residuals cannot identify the most discordant gene-specific age trend. "
              "Conflicting source percentages are retained in the observation hull.", "",
              "| Gene / source / ages | Observed % | Printed alternative % | Prediction age range %; uniform mean % | Uniform-assumption error pp | Interval gap pp | Source conflict? |",
              "|---|---:|---:|---:|---:|---:|---|"]
    for x in r["scored_rows"]:
        if x["role"] != "external_gene_marginal":
            continue
        o, p = x["observation"], x["source_prediction"]
        lines.append(f"| {x['gene']} / {x['cohort_id']} / {x['age_label']} | {pct(o['point'])} | {pct(o['printed_percentage_alternative'])} | "
                     f"{pct(p['lower'])}–{pct(p['upper'])}; {pct(p['uniform_age_mean'])} | "
                     f"{fmt(x['source_comparison']['signed_error_pp'])} | "
                     f"{fmt(x['source_comparison']['age_interval_gap_pp'])} | {o['source_conflict']} |")
    indexed = {x["record_id"]: x for x in r["scored_rows"]}
    lines += ["", "### Largest conditional deviations", "",
              "Finite-band age errors ranked by absolute error: " + "; ".join(
                  f"{rid}: {fmt(indexed[rid]['source_comparison']['signed_error_pp'])} pp"
                  for rid in r["largest_external_age_errors"][:3]) + ".", "",
              "Largest age-interval gaps, including open bands: " + "; ".join(
                  f"{rid}: {fmt(indexed[rid]['source_comparison']['age_interval_gap_pp'])} pp"
                  for rid in r["largest_external_age_gaps"][:3]) + ".", "",
              "Largest gene-marginal interval gaps: " + "; ".join(
                  f"{indexed[rid]['gene']} ({indexed[rid]['cohort_id']}, {indexed[rid]['age_label']}): "
                  f"{fmt(indexed[rid]['source_comparison']['age_interval_gap_pp'])} pp"
                  for rid in r["largest_gene_marginal_gaps"][:3]) + ". "
              "These identify disagreements of substituted endpoints under their stated bounds. "
              "They do not establish which gene's fitness or age trend is wrong."]
    lines += ["", "### All frozen parameter sensitivities", "",
              "Each scenario applies its common threshold to all external age rows. It is a hypothetical "
              "sensitivity, not a choice of the best cohort-specific cutoff.", "",
              "| Scenario | MAE pp | RMSE pp | Nonoverlap rows |", "|---|---:|---:|---:|"]
    for x in r["sensitivities"]:
        lines.append(f"| {x['scenario_id']} | {fmt(x['macro_mae_pp'])} | {fmt(x['macro_rmse_pp'])} | {x['nonoverlap_rows']} |")
    lines += ["", "### Read-support illustrations", "",
              "These are binomial read-sampling probabilities, not caller sensitivity estimates. "
              "Other filters, variable depths and variant classes remain unmodelled.", "",
              "| Cohort | Illustrative depth / alternate reads | Detection at true VAF 2% | 5% | 10% |",
              "|---|---:|---:|---:|---:|"]
    for x in r["detection_illustrations"]:
        lines.append(f"| {x['cohort_id']} | {x['depth']} / {x['min_alt_reads']} | {pct(x['at_vaf_02'])}% | "
                     f"{pct(x['at_vaf_05'])}% | {pct(x['at_vaf_10'])}% |")
    lines += ["", "### Falsifiable matched-assay predictions", "",
              "Fixed-age prevalence (%) for the source variant panels under perfect detection at true "
              "VAF >= 2%. These are predictions, not observations or clinical probabilities. A future "
              "matched aggregate should supply age support, carrier counts, denominators and caller "
              "calibration. For an age band, use its two endpoint predictions as the age-only bounds. "
              "These are separately fitted candidate forecasts, not a coherent joint carrier model; "
              "the subset-order check immediately below exposes the inconsistency.", "",
              "| Model | Age 40 | Age 50 | Age 60 | Age 70 | Age 80 | Age 90 |",
              "|---|---:|---:|---:|---:|---:|---:|"]
    for name in ("DNMT3A", "TET2", "ASXL1", "TP53", "POOLED10"):
        values = [x for x in r["future_predictions"] if x["model_id"] == name]
        lines.append("| " + name + " | " + " | ".join(pct(x["prediction"]) for x in values) + " |")
    coherence = r["subset_coherence"]
    lines += ["", "### Internal consistency of the separate published fits", "",
              "If a gene's variants are a subset of the pooled endpoint, its carrier prevalence "
              "cannot exceed the prevalence of any pooled-gene clone. This necessary event-inclusion "
              "inequality does not assume independent genes. The unchanged separate DFE fits violate "
              f"it in **{coherence['violation_count']} of {coherence['comparisons']}** grid comparisons; "
              f"the largest excess is **{fmt(coherence['maximum_gap_pp'])} pp**. Therefore these fits "
              "cannot jointly describe one common population/panel as exact probabilities. No fit is "
              "repaired or selected after seeing targets. This structural incompatibility concerns "
              "the assembled candidate predictions; it is not evidence against every individual "
              "fitness estimate or all branching models.", "",
              "| Gene | Age | Gene prevalence % | Pooled prevalence % | Excess pp |",
              "|---|---:|---:|---:|---:|"]
    for x in coherence["violations"]:
        lines.append(f"| {x['gene']} | {x['age']:g} | {pct(x['gene_prediction'])} | "
                     f"{pct(x['pooled_prediction'])} | {fmt(x['gap_pp'])} |")
    maximum_cutoff_difference = max(x["cutoff_difference_pp"] for x in r["future_predictions"])
    lines += ["", f"The maximum difference from the source notebook's VAF upper cutoff over this "
              f"prediction grid is **{fmt(maximum_cutoff_difference,8)} pp**. This checks a numerical "
              "convention; it does not measure biological model error.", "",
              "Report verification allows floating-point differences of at most "
              f"relative **{REPORT_FLOAT_REL_TOL:g}** or absolute **{REPORT_FLOAT_ABS_TOL:g}** "
              "for platform math-library and numerical-integration variation. Integer values, metadata, "
              "keys and input hashes are exact; this displayed results block must reproduce exactly. "
              "These are computational tolerances, not statistical uncertainty.", "", END]
    return "\n".join(lines)


def serialized_report(report: dict) -> str:
    return json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"


def write_manifest() -> None:
    """Record original project payloads; this does not update freeze evidence."""
    payloads = [p for p in HERE.rglob("*") if p.is_file()
                and p.name != "manifest.json" and "__pycache__" not in p.parts
                and p.suffix != ".pyc"]
    payloads.append(ROOT / "tests" / "test_research_ch_age_prevalence_test.py")
    manifest = dict(format_version=1, status="nonclinical; unreviewed; no medical advice",
                    normalization="CRLF and CR to LF before SHA-256; UTF-8 bytes otherwise unchanged",
                    path_policy="Every key constructed using Path.relative_to(ROOT).as_posix()",
                    files={p.relative_to(ROOT).as_posix(): digest(p) for p in sorted(payloads)})
    (HERE / "manifest.json").write_text(serialized_report(manifest), encoding="utf-8", newline="\n")


def assert_reports_close(actual, expected, location: str = "$report") -> None:
    """Compare report structure exactly and floating leaves within one budget.

    Numerical values in CSVs remain protected by exact normalized hashes.
    This function permits no missing/extra keys, reordered arrays, altered
    strings/booleans, or changed integers. The README renderer is checked
    separately without any tolerance.
    """
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            raise AssertionError(f"{location}: report keys differ")
        for key in expected:
            assert_reports_close(actual[key], expected[key], f"{location}[{key!r}]")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise AssertionError(f"{location}: report array differs")
        for index, (a, e) in enumerate(zip(actual, expected)):
            assert_reports_close(a, e, f"{location}[{index}]")
    elif isinstance(expected, float):
        if (isinstance(actual, bool) or not isinstance(actual, (int, float))
                or not math.isfinite(actual) or not math.isfinite(expected)
                or not math.isclose(actual, expected, rel_tol=REPORT_FLOAT_REL_TOL,
                                    abs_tol=REPORT_FLOAT_ABS_TOL)):
            raise AssertionError(f"{location}: numerical result differs: {actual!r} != {expected!r}")
    elif type(actual) is not type(expected) or actual != expected:
        raise AssertionError(f"{location}: exact report value differs: {actual!r} != {expected!r}")


def check_or_write(write: bool = False) -> dict:
    report = build_report()
    expected = serialized_report(report)
    path = HERE / "report.json"
    readme = HERE / "README.md"
    text = readme.read_text(encoding="utf-8")
    block = render_results(report)
    if START not in text or END not in text:
        raise ValueError("README lacks generated-section markers")
    before, rest = text.split(START, 1)
    _, after = rest.split(END, 1)
    full = before + block + after
    if write:
        path.write_text(expected, encoding="utf-8", newline="\n")
        readme.write_text(full, encoding="utf-8", newline="\n")
    else:
        assert_reports_close(json.loads(path.read_text(encoding="utf-8")), report)
        if text != full:
            raise AssertionError("README numerical results differ from recomputation")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="regenerate report, README results and payload manifest")
    args = parser.parse_args()
    report = check_or_write(args.write)
    if args.write:
        write_manifest()
    print(json.dumps({"strict_verdict": report["strict_verdict"],
                      "descriptive_summary": report["external_substitute_summary"]}, indent=2))


if __name__ == "__main__":
    main()
