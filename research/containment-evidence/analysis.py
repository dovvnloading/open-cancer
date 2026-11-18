#!/usr/bin/env python3
"""Reproduce a nonclinical aggregate evidence audit (Python standard library).

No network access, individual records, fitted patient trajectories, graph
digitisation, simulation of patients, or clinical recommendations. Eligibility
failures are results: a successful software test is not a biological test.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TABLES = ("studies", "aggregates", "comparisons", "mapping", "gates", "searches")
PROVENANCE = ("source_url", "locator", "access_date", "reuse_basis")
PARAMETERS = ("d", "c", "J", "A", "lambda")
BEGIN = "<!-- BEGIN RECOMPUTED EVIDENCE -->"
END = "<!-- END RECOMPUTED EVIDENCE -->"


def normalized_bytes(path: Path) -> bytes:
    """Hash logical UTF-8 text, treating CRLF, CR and LF identically."""
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(normalized_bytes(path)).hexdigest()


def read_tables(directory: Path = ROOT / "data") -> dict[str, list[dict[str, str]]]:
    result = {}
    for name in TABLES:
        text = normalized_bytes(directory / f"{name}.csv").decode("utf-8-sig")
        result[name] = list(csv.DictReader(io.StringIO(text, newline="")))
    return result


def number(row: dict[str, str], field: str = "value") -> float | None:
    return float(row[field]) if row.get(field, "") else None


def fmt(value: float | None) -> str:
    if value is None:
        return "not estimable"
    return f"{value:.6g}"


def escape(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def markdown_table(headers: list[str], rows: list[list[object]]) -> str:
    return "\n".join([
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
        *["| " + " | ".join(escape(cell) for cell in row) + " |" for row in rows],
    ])


def interval(row: dict[str, str]) -> str:
    low, high = number(row, "lower"), number(row, "upper")
    if low is None or high is None:
        return "not reported"
    return f"[{fmt(low)}, {fmt(high)}] ({row['uncertainty_type']})"


def displayed_value(row: dict[str, str] | None) -> str:
    if row is None:
        return "not reported"
    relation = row.get("relation", "exact")
    if relation == "not_reached":
        return "not reached/undefined"
    prefix = {"at_least": ">=", "greater_than": ">", "less_than": "<", "approximate": "~"}.get(relation, "")
    return prefix + fmt(number(row))


def validate(tables: dict[str, list[dict[str, str]]]) -> None:
    """Reject missing provenance, unsafe reinterpretations, and broken joins."""
    for name, rows in tables.items():
        if not rows:
            raise ValueError(f"Empty table: {name}")
        for index, row in enumerate(rows, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed CSV: {name}:{index}")
            for key in PROVENANCE:
                if not row.get(key):
                    raise ValueError(f"Missing {key}: {name}:{index}")
            if not row["source_url"].startswith("https://"):
                raise ValueError(f"Non-HTTPS source: {name}:{index}")
            if row["access_date"] != "2026-11-17":
                raise ValueError(f"Unexpected extraction date: {name}:{index}")
    studies = {s["study_id"]: s for s in tables["studies"]}
    aggregates = {r["row_id"]: r for r in tables["aggregates"]}
    if len(studies) != len(tables["studies"]) or len(aggregates) != len(tables["aggregates"]):
        raise ValueError("Duplicate study or aggregate ID")
    for row in aggregates.values():
        study = studies[row["study_id"]]
        if study["disposition"] not in {"included", "included_model", "included_single_arm"}:
            raise ValueError(f"Rights/eligibility exclusion contains numbers: {row['row_id']}")
        if "CC BY 4.0" not in row["reuse_basis"]:
            raise ValueError(f"Numeric reuse basis is not explicitly permissive: {row['row_id']}")
        if row["relation"] not in {"exact", "at_least", "greater_than", "less_than", "approximate", "not_reached"}:
            raise ValueError("Unrecognised censoring/rounding relation")
        if row["relation"] == "not_reached" and row["value"]:
            raise ValueError("Not-reached median cannot have a numerical value")
        for key in ("value", "lower", "upper"):
            value = number(row, key)
            if value is not None and (not math.isfinite(value) or (value < 0 and row["unit"] != "correlation")):
                raise ValueError(f"Invalid numeric aggregate: {row['row_id']}:{key}")
        low, high = number(row, "lower"), number(row, "upper")
        if (low is None) != (high is None) or (low is not None and low > high):
            raise ValueError(f"Invalid interval: {row['row_id']}")
        if low is not None and not row["uncertainty_type"]:
            raise ValueError("An interval must retain its published type")
    seen_comparisons = set()
    for comp in tables["comparisons"]:
        if comp["comparison_id"] in seen_comparisons:
            raise ValueError("Duplicate comparison ID")
        seen_comparisons.add(comp["comparison_id"])
        if comp["study_id"] not in studies:
            raise ValueError("Unknown comparison study")
        for key in ("experimental_row", "reference_row", "hr_row"):
            if comp[key]:
                row = aggregates[comp[key]]
                if row["study_id"] != comp["study_id"]:
                    raise ValueError("Comparison combines different study reports")
        if bool(comp["experimental_row"]) != bool(comp["reference_row"]):
            raise ValueError("An arm contrast requires both published summaries")
        if comp["hr_row"] and aggregates[comp["hr_row"]]["statistic"] != "hazard_ratio":
            raise ValueError("HR reference is not a published hazard ratio")
        if comp["experimental_row"]:
            exp, ref = aggregates[comp["experimental_row"]], aggregates[comp["reference_row"]]
            if exp["unit"] != ref["unit"] or exp["statistic"] != ref["statistic"]:
                raise ValueError("Incommensurable arm summaries")
    for study_id in studies:
        mappings = [m for m in tables["mapping"] if m["study_id"] == study_id]
        if sorted(m["parameter"] for m in mappings) != sorted(PARAMETERS):
            raise ValueError(f"Each study needs exactly five parameter mappings: {study_id}")
        gates = [g for g in tables["gates"] if g["study_id"] == study_id]
        if len(gates) != 1:
            raise ValueError(f"Each study needs one gate record: {study_id}")
    for row in tables["mapping"]:
        if row["status"] not in {"measured", "estimated", "conditional", "unidentifiable"}:
            raise ValueError("Unknown parameter status")
        for row_id in filter(None, row["input_ids"].split(";")):
            if row_id not in aggregates:
                raise ValueError("Mapping input missing")
        if row["status"] == "unidentifiable" and row["transformation"]:
            raise ValueError("Unidentifiable parameter cannot have a numerical transformation")


def contrast(comp: dict[str, str], aggregates: dict[str, dict[str, str]]) -> dict:
    result = dict(comp)
    exp = aggregates.get(comp["experimental_row"])
    ref = aggregates.get(comp["reference_row"])
    hr = aggregates.get(comp["hr_row"])
    result.update(experimental=displayed_value(exp), reference=displayed_value(ref),
                  difference=None, ratio=None, bound_type="", hr=None, hr_interval="not reported",
                  experimental_interval=interval(exp) if exp else "not reported",
                  reference_interval=interval(ref) if ref else "not reported",
                  interval_excludes_null=False, direction="not estimable", unit="")
    if exp and ref:
        result["unit"] = exp["unit"]
        if exp["relation"] == "exact" and ref["relation"] == "exact":
            x, y = number(exp), number(ref)
            if x is not None and y is not None:
                result["difference"] = x - y
                result["ratio"] = x / y if y > 0 else None
                result["direction"] = "experimental longer/higher" if x > y else (
                    "experimental shorter/lower" if x < y else "equal point summaries")
        elif exp["relation"] == "at_least" and ref["relation"] == "exact":
            x, y = number(exp), number(ref)
            if x is not None and y is not None and y > 0:
                result["difference"], result["ratio"] = x - y, x / y
                result["bound_type"] = ">="
                result["direction"] = "lower bound only"
    if hr:
        value, low, high = number(hr), number(hr, "lower"), number(hr, "upper")
        result.update(hr=value, hr_interval=interval(hr), unit="HR", hr_level=hr["uncertainty_type"])
        result["interval_excludes_null"] = bool(low is not None and high is not None and (high < 1 or low > 1))
        result["direction"] = "experimental lower hazard" if value < 1 else (
            "experimental higher hazard" if value > 1 else "equal point hazard")
        result["hr_low"], result["hr_high"] = low, high
    result["input_ids"] = ";".join(filter(None, [comp["experimental_row"], comp["reference_row"], comp["hr_row"]]))
    return result


def conditional_parameters(tables: dict[str, list[dict[str, str]]]) -> list[dict]:
    aggregates = {r["row_id"]: r for r in tables["aggregates"]}
    result = []
    for mapping in tables["mapping"]:
        if not mapping["transformation"]:
            continue
        values = [number(aggregates[key]) for key in mapping["input_ids"].split(";")]
        if mapping["transformation"] == "reciprocal_product":
            value = 1.0 / math.prod(values)
        else:
            raise ValueError("Unrecognised conditional transformation")
        result.append({**mapping, "value": value, "interval": "not identifiable from reported uncertainty"})
    return result


def threshold(d: float, c: float, j: float, a: float) -> float:
    """Only meaningful after the study's empirical/model eligibility gates."""
    if not (0 < d < c < j and a > 0):
        raise ValueError("Threshold requires 0 < d < c < J and A > 0")
    return math.log((j - d) / (c - d)) / a


def wilson(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Available for aggregate lineage counts; not applied to censored KM data."""
    if total <= 0 or not (0 <= successes <= total):
        raise ValueError("Invalid binomial counts")
    fraction, z2 = successes / total, z * z
    denominator = 1 + z2 / total
    centre = (fraction + z2 / (2 * total)) / denominator
    half = z * math.sqrt(fraction * (1 - fraction) / total + z2 / (4 * total * total)) / denominator
    return max(0.0, centre - half), min(1.0, centre + half)


def supply_from_zero_fraction(no_establishment: int, total: int) -> dict:
    low, high = wilson(no_establishment, total)
    return {"theta": -math.log(no_establishment / total) if no_establishment else math.inf,
            "lower": -math.log(high), "upper": -math.log(low) if low > 0 else math.inf}


def predictions(tables: dict[str, list[dict[str, str]]], contrasts: list[dict]) -> list[dict]:
    """Apply the frozen eligibility rules; do not force absent biological tests."""
    gate_rows = tables["gates"]
    requirements = {
        "P1": ("supply_series", "mean_burden_endpoint", "model_assumptions"),
        "P2": ("competition_contrast", "measured_supply", "mean_burden_endpoint", "model_assumptions"),
        "P2a": ("zero_competition", "mean_burden_endpoint", "model_assumptions"),
        "P3": ("complete_parameters", "mean_burden_endpoint", "model_assumptions"),
        "P4": ("new_lineage_counts", "measured_supply", "independent_A", "model_assumptions"),
        "P5": ("supply_series", "complete_parameters", "mean_burden_endpoint", "model_assumptions"),
    }
    reasons = {
        "P1": "No independently measured controlled mutation-supply series with the theorem's mean endpoint.",
        "P2": "No fixed-supply measured competition contrast; increasing competition is an auxiliary hypothesis, not a universal consequence of PR #61.",
        "P2a": "No measured zero-competition comparison satisfying the theorem's assumptions; absence of detected density dependence is not alpha=0.",
        "P3": "No report identifies d,c,J,A,lambda jointly; clinical medians and HRs do not identify the uncensored mean coefficients.",
        "P4": "No aggregate count of new established lineages with independent lambda and A; treatment failures are not founder counts.",
        "P5": "No otherwise comparable independently measured mutation-supply series; different cancers are not supply levels.",
    }
    output = []
    for prediction, fields in requirements.items():
        eligible = [r["study_id"] for r in gate_rows if all(r[field] == "yes" for field in fields)]
        # This evidence release contains no eligible mechanistic test. Failing
        # closed protects a later data extension from inheriting a stale verdict.
        if eligible:
            raise ValueError(f"New eligible {prediction} data require the prespecified statistical analysis: {eligible}")
        output.append({"prediction": prediction, "verdict": "untestable", "eligible": len(eligible),
                       "evidence": "data/gates.csv and data/mapping.csv (all study-specific source locators)",
                       "reason": reasons[prediction]})
    observed = [r for r in contrasts if r["role"] == "observed_latest" and r["kind"] == "outcome"]
    benefit = [r for r in observed if (r["hr"] is not None and r["hr"] < 1) or (
        r["hr"] is None and r["difference"] is not None and r["difference"] > 0)]
    adverse_interval = [r for r in observed if r["design"] == "randomized" and r["hr"] is not None
                        and r.get("hr_low") is not None and r["hr_low"] > 1]
    control_adverse = [r for r in adverse_interval if r["endpoint_class"] == "tumour_control"]
    control95_adverse = [r for r in control_adverse if r["hr_level"] == "95% CI"]
    output.append({"prediction": "P6a", "verdict": "consistent" if benefit else "untestable",
                   "eligible": len({r["independence_group"] for r in observed}),
                   "evidence": ";".join(r["comparison_id"] for r in benefit),
                   "reason": "Some observed comparisons favour reduced/adaptive exposure. This is compatibility with possible benefit, not a test of the mutation threshold; no pooled efficacy claim."})
    output.append({"prediction": "P6b", "verdict": "falsified" if adverse_interval else "untestable",
                   "eligible": len({r["independence_group"] for r in observed if r["design"] == "randomized"}),
                   "evidence": ";".join(r["comparison_id"] for r in adverse_interval),
                   "reason": "BROAD outcome reading of the ambiguous universal extension only: COIN OS has an adverse published 80% CI. Thus the original any-reported-interval rule is met. This is a descriptive contradiction at that level, NOT a 95% rejection, proof of harm from noninferiority failure, tumour-control falsification, or PR #61 falsification." if adverse_interval else
                   "No adverse source-reported randomized interval meets the broad-outcome rule."})
    for label, selected, suffix in (("P6b-control", control_adverse, "tumour-control endpoints"),
                                    ("P6b-control95", control95_adverse, "tumour-control endpoints and published 95% intervals")):
        output.append({"prediction": label, "verdict": "falsified" if selected else "untestable",
                       "eligible": len({r["independence_group"] for r in observed if r["design"] == "randomized"
                                        and r["endpoint_class"] == "tumour_control"}),
                       "evidence": ";".join(r["comparison_id"] for r in selected) or "gonzalez_cao_pfs;interim_pfs",
                       "reason": f"Explicit POST-EXTRACTION sensitivity restricted to {suffix}: INTERIM's PFS CI spans 1; Gonzalez-Cao supplies arm-median CIs, not a contrast CI. Original endpoint ambiguity prevents a unique fully pre-specified P6b verdict. Adverse point estimates remain visible."})
    return output


def compute(tables: dict[str, list[dict[str, str]]]) -> dict:
    validate(tables)
    aggregates = {r["row_id"]: r for r in tables["aggregates"]}
    contrasts = [contrast(row, aggregates) for row in tables["comparisons"]]
    return {
        "counts": {"study_records": len(tables["studies"]), "aggregate_rows": len(aggregates),
                   "numeric_reports": len({r["study_id"] for r in aggregates.values()}),
                   "comparison_rows": len(contrasts), "search_records": len(tables["searches"]),
                   "dispositions": dict(sorted(Counter(s["disposition"] for s in tables["studies"]).items()))},
        "contrasts": contrasts, "conditional_parameters": conditional_parameters(tables),
        "predictions": predictions(tables, contrasts),
        "formal_mechanistic_tests": 0, "pooled_effects": 0,
    }


def render(tables: dict[str, list[dict[str, str]]], result: dict) -> str:
    counts = result["counts"]
    parts = [BEGIN, "", "## Recomputed evidence tables", "",
             f"The ledger contains **{counts['study_records']} source/study records**, "
             f"**{counts['numeric_reports']} reports with licensed aggregates**, "
             f"**{counts['aggregate_rows']} aggregate rows**, "
             f"**{counts['comparison_rows']} contrasts**, and "
             f"**{counts['search_records']} search records. These are extraction counts, not independent studies or participants.**",
             "", "Disposition counts: " + "; ".join(f"{key}={value}" for key, value in counts["dispositions"].items()) + ".",
             "", "### Observed outcome contrasts", "",
             "Differences and ratios compare published arm summaries; they are not hazard ratios or uncensored means. "
             "No CI for a difference/ratio is invented. `>=` is a lower bound; `not reached/undefined` is not a numerical median. "
             "Published HRs retain experimental/reference orientation and confidence level. "
             "Earlier reports from the same cohort are explicitly marked. Shared-control rows are dependent.", ""]
    observed_rows, simulated_rows, dose_rows = [], [], []
    for comp in result["contrasts"]:
        bound = comp["bound_type"]
        row = [comp["comparison_id"], comp["endpoint"], comp["experimental"], comp["reference"],
               (bound + fmt(comp["difference"])) if comp["difference"] is not None else "not estimable",
               (bound + fmt(comp["ratio"])) if comp["ratio"] is not None else "not estimable",
               fmt(comp["hr"]) if comp["hr"] is not None else "not reported", comp["hr_interval"],
               comp["role"], f"[{comp['study_id']}]({comp['source_url']}) / `{comp['input_ids']}`"]
        if comp["kind"] == "dose":
            dose_rows.append(row)
        elif comp["role"] == "model_counterfactual":
            simulated_rows.append(row)
        else:
            observed_rows.append(row)
    headers = ["Contrast", "Endpoint / units", "Experimental", "Reference", "Difference", "Ratio", "HR", "HR interval", "Role", "Source / input rows"]
    parts.append(markdown_table(headers, observed_rows))
    arm_intervals = [[c["comparison_id"], c["experimental_interval"], c["reference_interval"]]
                     for c in result["contrasts"] if c["experimental_interval"] != "not reported"
                     or c["reference_interval"] != "not reported"]
    if arm_intervals:
        parts.extend(["", "Published arm-summary intervals (not intervals for the derived ratio or difference):", "",
                      markdown_table(["Contrast", "Experimental uncertainty", "Reference uncertainty"], arm_intervals)])
    parts.extend(["", "### Exposure contrasts", "",
                  "Cumulative amount, mean daily intensity and on-treatment time are different quantities. "
                  "Longer survival changes the accumulation window. Doses below have no published sampling interval unless listed in the inputs.", "",
                  markdown_table(headers, dose_rows), "", "### Published fitted-model counterfactuals", "",
                  "These summaries came from published simulations fitted to clinical cohorts. They are not observed treatment comparisons or new independent clinical replications.", "",
                  markdown_table(headers, simulated_rows), "", "### Other published aggregates", ""])
    used = {key for comp in result["contrasts"] for key in comp["input_ids"].split(";") if key}
    additional = [[r["row_id"], r["measure"], displayed_value(r), r["unit"], interval(r),
                   r["statistic"], f"[{r['study_id']}]({r['source_url']})"]
                  for r in tables["aggregates"] if r["row_id"] not in used]
    parts.append(markdown_table(["Input row", "Measure", "Value", "Unit", "Reported uncertainty", "Statistic", "Source"], additional))
    parts.extend(["", "### Parameter identifiability by report", "",
                  "`conditional` denotes a point transformation inside the source's fitted model, with no transferable uncertainty interval. "
                  "It does not satisfy PR #61's eligibility gates. Full reasons and source locators are in `data/mapping.csv`.", ""])
    mapped = {(r["study_id"], r["parameter"]): r["status"] for r in tables["mapping"]}
    parts.append(markdown_table(["Report", *PARAMETERS], [
        [f"[{s['study_id']}]({s['source_url']})", *[mapped[s["study_id"], p] for p in PARAMETERS]]
        for s in tables["studies"]]))
    if result["conditional_parameters"]:
        parts.extend(["", "Conditional transformations (not tumour threshold estimates):", "",
                      markdown_table(["Report", "Parameter", "Value (days)", "Input rows", "Uncertainty"], [
                          [r["study_id"], r["parameter"], fmt(r["value"]), r["input_ids"], r["interval"]]
                          for r in result["conditional_parameters"]])])
    parts.extend(["", "### Prediction verdicts and declared sensitivities", "",
                  "Eligibility counts are report counts for mechanistic rows and distinct comparison groups for P6. "
                  "No mechanistic slope, threshold, establishment or crossing test passed its data gate; "
                  "**0 formal mechanistic hypothesis tests and 0 pooled effects** were calculated. "
                  "Holm correction therefore has no p-values to adjust. P6 uses descriptive published intervals only.", "",
                  markdown_table(["Prediction", "Verdict", "Eligible", "Evidence rows", "Interpretation"], [
                      [p["prediction"], p["verdict"], p["eligible"], p["evidence"], p["reason"]]
                      for p in result["predictions"]]), "", "### Search coverage and exclusions", "",
                  "Every screened source below is retained regardless of outcome or extraction permission. "
                  "Rights/access exclusions and missing posted results are not negative efficacy outcomes. "
                  "Original dispositions, full reuse basis and exact source locators remain in the CSVs.", "",
                  markdown_table(["Source", "Design / setting", "Disposition", "Limitation / reason"], [
                      [f"[{s['study_id']}: {s['title']}]({s['source_url']})", s["design"] + "; " + s["setting"],
                       s["disposition"], s["limitations"]] for s in tables["studies"]]), "", END])
    return "\n".join(parts)


def update_readme(block: str, check: bool = False) -> None:
    path = ROOT / "README.md"
    text = normalized_bytes(path).decode("utf-8")
    if text.count(BEGIN) != 1 or text.count(END) != 1:
        raise ValueError("README must contain exactly one results block")
    before, rest = text.split(BEGIN)
    _, after = rest.split(END)
    expected = before + block + after
    if check:
        if text != expected:
            raise ValueError("README evidence tables are stale; run analysis.py --write-readme")
    else:
        path.write_text(expected, encoding="utf-8", newline="\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Print complete recomputed result JSON")
    parser.add_argument("--check", action="store_true", help="Check the checked-in README tables")
    parser.add_argument("--write-readme", action="store_true", help="Refresh only this directory's README result block")
    args = parser.parse_args()
    tables = read_tables()
    result = compute(tables)
    block = render(tables, result)
    if args.check or args.write_readme:
        update_readme(block, check=args.check)
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    elif args.check:
        print("Aggregate/provenance checks and README recomputation passed (nonclinical).")
    elif not args.write_readme:
        print(block)


if __name__ == "__main__":
    main()
