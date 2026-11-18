"""Recompute the aggregate-only CH transport stress test using the stdlib.

The primary mechanistic verdict is UNTESTABLE: gene-presence Cox coefficients
and baseline gene-class growth terms do not identify a matched-state cumulative
hazard ratio. The separately named proxy stress test can still be falsified.
No clinical prediction, pooled meta-analysis, or patient-level data are used.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import statistics

import model
import simulate

ROOT = Path(__file__).resolve().parent
BEGIN = "<!-- BEGIN RECOMPUTED RESULTS -->"
END = "<!-- END RECOMPUTED RESULTS -->"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def load_inputs(root: Path = ROOT) -> dict:
    return {"fitness": read_csv(root / "data" / "fitness.csv"),
            "outcomes": read_csv(root / "data" / "outcomes.csv"),
            "cohorts": read_csv(root / "data" / "cohorts_outcomes.csv"),
            "fitness_cohorts": read_csv(root / "data" / "cohorts_fitness.csv"),
            "scenarios": read_csv(root / "analysis_scenarios.csv")}


def select_proxy_pairs(data: dict) -> tuple[list[dict], list[dict]]:
    """Complete-case rule fixed by measurement scope, never outcome magnitude.

    A unique gene-class baseline is still only a PROXY for whole-gene fitness.
    Multiple classes are not averaged; sole-hotspot terms are kept separate.
    """
    by_gene: dict[str, list[dict]] = {}
    for fit in data["fitness"]:
        by_gene.setdefault(fit["gene"], []).append(fit)
    pairs, excluded = [], []
    for outcome in data["outcomes"]:
        if outcome["study_id"] != "fabre2022_abelson_refit":
            continue
        matches = by_gene.get(outcome["gene"], [])
        reason = ""
        if not matches:
            reason = "no published numerical fitness term with interval"
        elif len(matches) > 1:
            reason = "multiple mutation classes; no published pooling weights"
        elif matches[0]["estimator_scope"] != "gene_class_baseline_effect":
            reason = "hotspot fitness term versus whole-gene outcome"
        elif not all(outcome.get(k) for k in ("estimate", "lower", "upper")):
            reason = "missing outcome interval"
        elif not all(math.isfinite(float(outcome[k])) and float(outcome[k]) > 0
                     for k in ("estimate", "lower", "upper")):
            reason = "nonpositive or nonfinite ratio; no pseudocounts allowed"
        if reason:
            excluded.append({"gene": outcome["gene"], "reason": reason,
                             "outcome_id": outcome["outcome_id"]})
        else:
            pairs.append({"gene": outcome["gene"], "fitness": matches[0],
                          "outcome": outcome})
    return sorted(pairs, key=lambda p: p["gene"]), sorted(excluded, key=lambda p: p["gene"])


def opportunity(s: float, scenario: dict) -> float:
    horizon = float(scenario["horizon_years"])
    if scenario["model"] == "exponential":
        return model.exposure(s, horizon)
    if scenario["model"] == "logistic":
        return model.logistic_exposure(s, 2 * float(scenario["vaf"]), horizon)
    raise ValueError("unknown model in analysis_scenarios.csv")


def _components(pair: dict, scenario: dict) -> dict:
    fit, outcome = pair["fitness"], pair["outcome"]
    return {
        "x": math.log(opportunity(float(fit["estimate"]), scenario)),
        "xlo": math.log(opportunity(float(fit["lower"]), scenario)),
        "xhi": math.log(opportunity(float(fit["upper"]), scenario)),
        "y": math.log(float(outcome["estimate"])),
        "ylo": math.log(float(outcome["lower"])),
        "yhi": math.log(float(outcome["upper"])),
    }


def calibration(pairs: list[dict], scenario: dict) -> tuple[float, float, float]:
    values = [_components(p, scenario) for p in pairs]
    if not values:
        raise ValueError("calibration requires training genes")
    return (statistics.fmean(v["y"] - v["x"] for v in values),
            statistics.fmean(v["ylo"] - v["xhi"] for v in values),
            statistics.fmean(v["yhi"] - v["xlo"] for v in values))


def loo(pairs: list[dict], scenario: dict) -> dict:
    if len(pairs) < 3:
        raise ValueError("LOO summary requires at least three eligible genes")
    rows = []
    for held in pairs:
        train = [p for p in pairs if p["gene"] != held["gene"]]
        alpha, alo, ahi = calibration(train, scenario)
        value = _components(held, scenario)
        prediction = alpha + value["x"]
        plo, phi = alo + value["xlo"], ahi + value["xhi"]
        null_prediction = statistics.fmean(math.log(float(p["outcome"]["estimate"]))
                                            for p in train)
        overlap = max(value["ylo"], plo) <= min(value["yhi"], phi)
        rows.append({
            "gene": held["gene"],
            "fitness": float(held["fitness"]["estimate"]),
            "fitness_lower": float(held["fitness"]["lower"]),
            "fitness_upper": float(held["fitness"]["upper"]),
            "observed": math.exp(value["y"]),
            "observed_lower": math.exp(value["ylo"]),
            "observed_upper": math.exp(value["yhi"]),
            "prediction": math.exp(prediction),
            "prediction_lower": math.exp(plo), "prediction_upper": math.exp(phi),
            "log_error": value["y"] - prediction,
            "observed_over_predicted": math.exp(value["y"] - prediction),
            "null_log_error": value["y"] - null_prediction,
            "envelopes_overlap": overlap,
            "transport_verdict": "consistent" if overlap else "falsified",
        })
    rmse = math.sqrt(statistics.fmean(r["log_error"]**2 for r in rows))
    null_rmse = math.sqrt(statistics.fmean(r["null_log_error"]**2 for r in rows))
    failed = [r["gene"] for r in rows if not r["envelopes_overlap"]]
    return {"scenario_id": scenario["scenario_id"], "model": scenario["model"],
            "horizon_years": float(scenario["horizon_years"]), "vaf": scenario["vaf"],
            "n_genes": len(rows), "rows": rows, "log_rmse": rmse,
            "null_log_rmse": null_rmse, "rmse_ratio_to_null": rmse / null_rmse,
            "median_absolute_log_error": statistics.median(abs(r["log_error"]) for r in rows),
            "failed_genes": failed,
            "proxy_verdict": "falsified" if failed else "consistent",
            "mechanistic_verdict": "untestable"}


def joint_box_audit(pairs: list[dict], scenario: dict) -> dict:
    """POST HOC hostile-review check, distinct from the frozen LOO criterion.

    A single alpha can pass through all rectangular input intervals iff the
    intersection of [ylo-xhi,yhi-xlo] is nonempty. Intervals are treated as
    chosen sensitivity boxes, not guaranteed simultaneous confidence sets.
    """
    constraints = [(p["gene"], _components(p, scenario)) for p in pairs]
    lower_gene, lower = max(((g, v["ylo"] - v["xhi"]) for g, v in constraints),
                            key=lambda item: item[1])
    upper_gene, upper = min(((g, v["yhi"] - v["xlo"]) for g, v in constraints),
                            key=lambda item: item[1])
    return {"scenario_id": scenario["scenario_id"], "alpha_lower": lower,
            "alpha_upper": upper, "lower_limiting_gene": lower_gene,
            "upper_limiting_gene": upper_gene, "gap": lower - upper,
            "feasible": lower <= upper, "status": "post_hoc_sensitivity_box_check"}


def excluded_gene_sensitivity(data: dict, pairs: list[dict], scenario: dict) -> list[dict]:
    """Never pool class means: show the range of class/hotspot extrapolations.

    These are baseline-term ranges, not uncertainty bounds on total gene
    fitness; omitted site/clone effects can extend beyond them.
    """
    alpha, alo, ahi = calibration(pairs, scenario)
    included = {p["gene"] for p in pairs}
    results = []
    for outcome in data["outcomes"]:
        if outcome["study_id"] != "fabre2022_abelson_refit" or outcome["gene"] in included:
            continue
        fits = [f for f in data["fitness"] if f["gene"] == outcome["gene"]]
        if not fits:
            continue
        centers = [math.exp(alpha) * opportunity(float(f["estimate"]), scenario) for f in fits]
        lows = [math.exp(alo) * opportunity(float(f["lower"]), scenario) for f in fits]
        highs = [math.exp(ahi) * opportunity(float(f["upper"]), scenario) for f in fits]
        observed = float(outcome["estimate"])
        results.append({"gene": outcome["gene"], "classes": [f["mutation_class"] for f in fits],
                        "observed": observed, "point_prediction_min": min(centers),
                        "point_prediction_max": max(centers), "envelope_min": min(lows),
                        "envelope_max": max(highs), "fold_residual_min": observed / max(centers),
                        "fold_residual_max": observed / min(centers)})
    return sorted(results, key=lambda r: r["gene"])


def reference_group_audit(cohort: dict) -> dict:
    counts = {key: int(cohort[key]) for key in
              ("total_n", "total_cases", "carrier_n", "carrier_cases", "sweep_n", "sweep_cases")}
    a = counts["sweep_cases"]
    b = counts["sweep_n"] - a
    all_c = counts["total_cases"] - a
    all_d = counts["total_n"] - counts["sweep_n"] - all_c
    carrier_c = counts["carrier_cases"] - a
    carrier_d = counts["carrier_n"] - counts["sweep_n"] - carrier_c
    if min(a, b, all_c, all_d, carrier_c, carrier_d) <= 0:
        raise ValueError("positive aggregate contingency cells required")
    all_or = a * all_d / (b * all_c)
    carrier_or = a * carrier_d / (b * carrier_c)
    return {**counts, "sweep_controls": b, "non_sweep_cases_all": all_c,
            "non_sweep_controls_all": all_d, "non_sweep_cases_carriers": carrier_c,
            "non_sweep_controls_carriers": carrier_d, "all_participant_or": all_or,
            "within_carriers_or": carrier_or, "comparator_or_factor": all_or / carrier_or,
            "mechanistic_verdict": "untestable"}


def run(root: Path = ROOT, include_simulations: bool = True) -> dict:
    data = load_inputs(root)
    pairs, excluded = select_proxy_pairs(data)
    scenarios = [loo(pairs, s) for s in data["scenarios"]]
    primary_spec = next(s for s in data["scenarios"] if s["scenario_id"] == "primary_10y")
    primary = next(s for s in scenarios if s["scenario_id"] == "primary_10y")
    return {"status": "nonclinical_unreviewed", "mechanistic_verdict": "untestable",
            "proxy_verdict": primary["proxy_verdict"],
            "n_fitness_rows": len(data["fitness"]),
            "n_fitness_genes": len({r["gene"] for r in data["fitness"]}),
            "n_gene_outcomes": sum(r["study_id"] == "fabre2022_abelson_refit" for r in data["outcomes"]),
            "n_outcome_rows": len(data["outcomes"]),
            "fitness_cohorts": data["fitness_cohorts"],
            "vaf_outcome": next(r for r in data["outcomes"] if r["study_id"] == "watson2026_ukb"),
            "primary": primary, "scenarios": scenarios, "excluded": excluded,
            "joint_box_audit": [joint_box_audit(pairs, s) for s in data["scenarios"]],
            "excluded_gene_sensitivity": excluded_gene_sensitivity(data, pairs, primary_spec),
            "reference_group_audit": reference_group_audit(data["cohorts"][0]),
            "simulations": simulate.run_scenarios(root) if include_simulations else []}


def render_report(result: dict) -> str:
    p, ref = result["primary"], result["reference_group_audit"]
    lines = [BEGIN, "## Recomputed quantitative results", "",
             "This entire section is regenerated from the CSVs by `analysis.py`; the test compares it byte for byte.", "",
             f"The extract contains **{result['n_fitness_rows']} fitness terms for {result['n_fitness_genes']} genes**, "
             f"**{result['n_gene_outcomes']} AML gene coefficients**, and one separately sourced VAF-stratified outcome row.",
             "The fitness terms are baseline gene-class effects or explicitly identified hotspot effects, not published pooled gene means. "
             "AML coefficients are from the Fabre reanalysis of Abelson data. The exact refit sample size and coefficient covariance "
             "are not supplied in that aggregate file; neither is reconstructed from participants.", "",
             f"**Conditional fitness-to-progression model: UNTESTABLE. Fixed-shape baseline-proxy transport: {p['proxy_verdict'].upper()}.** "
             "The latter means falsified under the stated interval-box and transport stress-test rule. "
             "CI/HPDI limits are probabilistic intervals, not guaranteed parameter bounds; this is not a formal statistical rejection. "
             "The verdict concerns the compound mapping assumptions, not the branching-process identities or fitness as a useful predictor.", "",
             "### Published growth-cohort denominators", "",
             "| Metric | Value | Unit |", "|---|---:|---|"]
    lines += [f"| {r['metric']} | {r['value']} | {r['unit']} |" for r in result["fitness_cohorts"]]
    lines += ["",
             f"### Primary proxy test: {p['horizon_years']:g}-year opportunity, {p['n_genes']} eligible genes", "",
             "Intervals below combine published fitness HPD limits with reconstructed Wald limits from the published log-HR and SE. "
             "The prediction ranges are rectangular sensitivity envelopes, not calibrated prediction intervals. "
             "Training excludes the displayed gene. HR is a gene-presence Cox coefficient, not a cumulative-hazard ratio at fixed VAF.", "",
             "| Gene | Baseline fitness per year [HPDI] | Observed HR [derived Wald limits] | LOO prediction [envelope] | Observed / predicted | Envelope verdict |",
             "|---|---:|---:|---:|---:|---|"]
    for r in p["rows"]:
        lines.append(f"| {r['gene']} | {r['fitness']:.5f} [{r['fitness_lower']:.5f}, {r['fitness_upper']:.5f}] "
                     f"| {r['observed']:.3f} [{r['observed_lower']:.3f}, {r['observed_upper']:.3f}] "
                     f"| {r['prediction']:.3f} [{r['prediction_lower']:.3f}, {r['prediction_upper']:.3f}] "
                     f"| {r['observed_over_predicted']:.3f} | {r['transport_verdict']} |")
    ranked = sorted(p["rows"], key=lambda r: abs(r["log_error"]), reverse=True)
    lines += ["", f"LOO log-RMSE is **{p['log_rmse']:.6f}**, versus **{p['null_log_rmse']:.6f}** for the "
              f"intercept-only null (ratio **{p['rmse_ratio_to_null']:.6f}**). Median absolute log error is "
              f"**{p['median_absolute_log_error']:.6f}**. These are descriptive gene-level prediction errors, not patient discrimination or calibration.",
              "", "Largest absolute held-out log deviations, in descending order: " +
              "; ".join(f"**{r['gene']}** ({r['observed_over_predicted']:.3f} observed/predicted)" for r in ranked) + ".",
              "", "### Horizon and saturation attempts to rescue the proxy", "",
              "Logistic rows use hypothetical common starting VAFs; no VAF bin midpoint is imputed. The same source genes and LOO rule are retained.", "",
              "| Scenario | Years | Assumed VAF | Log-RMSE | Null log-RMSE | Nonoverlapping genes | Proxy verdict |",
              "|---|---:|---:|---:|---:|---|---|"]
    for s in result["scenarios"]:
        lines.append(f"| {s['scenario_id']} | {s['horizon_years']:g} | {s['vaf'] or 'small-clone limit'} "
                     f"| {s['log_rmse']:.6f} | {s['null_log_rmse']:.6f} "
                     f"| {', '.join(s['failed_genes']) or 'none'} | {s['proxy_verdict']} |")
    long_result = next(s for s in result["scenarios"] if s["scenario_id"] == "long_20y")
    long_u2af1 = next(r for r in long_result["rows"] if r["gene"] == "U2AF1")
    lines += ["", f"**The {long_result['horizon_years']:g}-year sensitivity is consistent under the frozen held-out envelope rule.** "
              f"For U2AF1 its center is {long_u2af1['prediction']:.3f}, envelope "
              f"[{long_u2af1['prediction_lower']:.3f}, {long_u2af1['prediction_upper']:.3f}], "
              f"and observed/predicted remains {long_u2af1['observed_over_predicted']:.3f}. "
              "The primary rejection is therefore not robust to every prespecified horizon under that rule; "
              "overlap does not establish predictive accuracy or joint calibration.",
              "", "### Additional hostile-review check: one common intercept", "",
              "**Post hoc; not the frozen decision rule.** Each gene permits alpha in "
              "`[log(outcome_lower)-log(opportunity_upper), log(outcome_upper)-log(opportunity_lower)]`. "
              "One common alpha exists within all chosen input rectangles exactly when the maximum lower endpoint is no greater "
              "than the minimum upper endpoint. Necessity follows by applying each interval constraint; sufficiency follows by "
              "choosing any alpha in their intersection and a feasible x and y in each rectangle. A positive gap below means "
              "infeasibility within those boxes, not a statistical rejection. This explains why separately calibrated LOO envelopes "
              "can overlap despite no single jointly compatible calibration.", "",
              "| Scenario | Max lower alpha (gene) | Min upper alpha (gene) | Gap | Common alpha feasible? |",
              "|---|---:|---:|---:|---|"]
    for row in result["joint_box_audit"]:
        lines.append(f"| {row['scenario_id']} | {row['alpha_lower']:.6f} ({row['lower_limiting_gene']}) "
                     f"| {row['alpha_upper']:.6f} ({row['upper_limiting_gene']}) | {row['gap']:.6f} "
                     f"| {'yes' if row['feasible'] else 'no'} |")
    lines += ["", "### Complete eligibility accounting", "",
              "Every gene outcome is either in the primary table or below; no removal depends on risk magnitude.", "",
              "| Excluded from primary proxy | Reason |", "|---|---|"]
    lines += [f"| {r['gene']} | {r['reason']} |" for r in result["excluded"]]
    lines += ["", "For excluded genes that have class/hotspot fitness terms, the following exploratory extrapolation keeps all "
              "reported classes and shows the range of their predictions. Calibration uses all primary genes; no excluded gene's "
              "outcome is used to calibrate its own prediction. This range does not bound whole-gene fitness, because unreported "
              "mixture weights and site/clone effects remain unknown.", "",
              "| Gene | Available class/hotspot terms | Observed HR | Point prediction range | Envelope over terms | Observed / predicted range |",
              "|---|---|---:|---:|---:|---:|"]
    for r in result["excluded_gene_sensitivity"]:
        lines.append(f"| {r['gene']} | {', '.join(r['classes'])} | {r['observed']:.3f} "
                     f"| {r['point_prediction_min']:.3f}–{r['point_prediction_max']:.3f} "
                     f"| {r['envelope_min']:.3f}–{r['envelope_max']:.3f} "
                     f"| {r['fold_residual_min']:.3f}–{r['fold_residual_max']:.3f} |")
    vaf = result["vaf_outcome"]
    lines += ["", "### VAF reference-group falsification attempt", "",
              "Watson et al. (Nature, November 2026), Fig. 5f and Methods, report a VAF threshold strictly above "
              f"**{float(vaf['vaf_min']):g}** and a published OR of **{float(vaf['estimate']):g} "
              f"[{float(vaf['lower']):g}, {float(vaf['upper']):g}]** versus all other participants. "
              "Its endpoint is future **myeloid leukaemia (ICD-10 C92)**, not AML alone. "
              "The counts below are published nested cohort aggregates; comparison groups are not independent datasets.", "",
              "| Population | Participants | Future cases |", "|---|---:|---:|",
              f"| All eligible UKB participants | {ref['total_n']} | {ref['total_cases']} |",
              f"| Detected hotspot carriers | {ref['carrier_n']} | {ref['carrier_cases']} |",
              f"| High-VAF sweeps | {ref['sweep_n']} | {ref['sweep_cases']} |", "",
              "| Comparison | Sweep cases / noncases | Comparator cases / noncases | Derived cross-product OR |",
              "|---|---:|---:|---:|",
              f"| Sweep versus all other participants | {ref['sweep_cases']} / {ref['sweep_controls']} "
              f"| {ref['non_sweep_cases_all']} / {ref['non_sweep_controls_all']} | {ref['all_participant_or']:.6f} |",
              f"| Sweep versus other detected carriers | {ref['sweep_cases']} / {ref['sweep_controls']} "
              f"| {ref['non_sweep_cases_carriers']} / {ref['non_sweep_controls_carriers']} | {ref['within_carriers_or']:.6f} |", "",
              f"The comparator change produces a **{ref['comparator_or_factor']:.6f}-fold** difference in the OR without "
              "changing the sweep group. The within-carrier value is a new arithmetic contrast from published counts, "
              "not a reported adjusted effect. No confidence interval is assigned to it here. Neither contrast tests linear "
              "VAF scaling: within-bin sizes, fitness, common follow-up horizon and relevant confounder distributions are unavailable. "
              "The VAF-conditioned mechanistic verdict is **UNTESTABLE**.", "",
              "### Independent synthetic simulation checks", "",
              "Parameters and seeds are in `model_scenarios.csv`; all are synthetic. The check uses a separate Gillespie event simulation. "
              "SE is computed from the exact model probability and the stated replication count, solely for numerical diagnostics.", "",
              "| Scenario | Replicates | Marks | Simulated probability | Exact probability | Error / model SE | Exact cumulative hazard |",
              "|---|---:|---:|---:|---:|---:|---:|"]
    for r in result["simulations"]:
        lines.append(f"| {r['scenario_id']} | {r['replicates']} | {r['events']} | {r['simulation_probability']:.6f} "
                     f"| {r['exact_probability']:.6f} | {r['standardized_error']:.4f} "
                     f"| {r['exact_cumulative_hazard']:.6f} |")
    lines += ["", "The equal-net-fitness supercritical and high-turnover scenarios differ in exact risk. "
              "That is a mathematical counterexample to net-fitness sufficiency at finite marking rate, not an empirical estimate of CH turnover.", END]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--markdown", action="store_true", help="print the generated README section")
    parser.add_argument("--write-readme", action="store_true", help="replace only the marked README result section")
    args = parser.parse_args()
    result = run()
    if args.write_readme:
        path = ROOT / "README.md"
        text = path.read_text(encoding="utf-8")
        before, remaining = text.split(BEGIN, 1)
        _, after = remaining.split(END, 1)
        path.write_text(before + render_report(result) + after, encoding="utf-8", newline="\n")
    elif args.markdown:
        print(render_report(result))
    else:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
