#!/usr/bin/env python3
"""Audit published stage-specific AAPCs. Nonclinical; no medical advice.

Only the standard library and the attributed aggregate CSV are used. Exact
decimal comparisons decide interval compatibility; math.log1p evaluates the
reported natural-log summaries. This is not a simultaneous significance test.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
DATA_PATH = ROOT / "data" / "lee_argov_table2.csv"
AGES = ("65-74", "75-84", "85+")
STAGES = ("localized", "regional", "distant", "unknown")
MODES = ("primary", "printed_ci", "full_unit", "point_rounding")
BEGIN = "<!-- BEGIN RECOMPUTED RESULTS -->"
END = "<!-- END RECOMPUTED RESULTS -->"
NUMERIC_FIELDS = ("printed_aapc", "printed_ci_lower", "printed_ci_upper")
PROVENANCE_FIELDS = (
    "source", "doi", "source_url", "table", "source_section", "source_row",
    "source_column", "locator", "period_start", "period_end", "measure",
    "units", "interval_level", "interval_type", "rate_adjustment",
    "dataset_version", "access_date", "access_timezone", "reuse_basis",
    "reuse_url",
)
SOURCE_COLUMNS = dict(zip(AGES, ("65-74 y", "75-84 y", "≥85 y")))
FIXED_METADATA = {
    "doi": "10.1001/jamanetworkopen.2025.16947",
    "source_url": "https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2835519",
    "table": "Table 2",
    "source_section": "Stage at diagnosis",
    "period_start": "2001",
    "period_end": "2019",
    "measure": "AAPC",
    "units": "percent per year",
    "interval_level": "95%",
    "interval_type": "published marginal confidence interval",
    "reuse_url": "https://jamanetwork.com/pages/cc-by-license-permissions",
}


def normalized_sha256(path: Path) -> str:
    """Hash UTF-8 text after CRLF/CR to LF normalization; preserve EOF."""
    content = Path(path).read_text(encoding="utf-8")
    content = content.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def file_hashes(root: Path = ROOT) -> dict[str, str]:
    prefix = Path("research") / "breast-age-period-stage"
    paths = (Path("PROTOCOL.md"), Path("data") / "lee_argov_table2.csv",
             Path("analysis.py"))
    return {(prefix / path).as_posix(): normalized_sha256(root / path)
            for path in paths}


def half_unit(printed: str) -> Decimal:
    """Half of the last printed decimal unit, in AAPC percentage points."""
    if not isinstance(printed, str) or not re.fullmatch(r"-?\d+\.\d+", printed):
        raise ValueError(f"Invalid printed decimal: {printed!r}")
    digits = len(printed.partition(".")[2])
    return Decimal(5).scaleb(-digits - 1)


def percent_interval(cell: dict, mode: str = "primary") -> tuple[Decimal, Decimal]:
    """Build exact percent endpoints; widen before the nonlinear transform."""
    if mode not in MODES:
        raise ValueError(f"Unknown interval mode: {mode}")
    for field in NUMERIC_FIELDS:
        half_unit(cell[field])
    point, lower, upper = (Decimal(cell[field]) for field in NUMERIC_FIELDS)
    if not lower <= point <= upper:
        raise ValueError("Invalid printed interval: estimate outside endpoints")
    if mode == "point_rounding":
        q = half_unit(cell["printed_aapc"])
        low, high = point - q, point + q
    else:
        multiplier = {"primary": 1, "printed_ci": 0, "full_unit": 2}[mode]
        low = lower - multiplier * half_unit(cell["printed_ci_lower"])
        high = upper + multiplier * half_unit(cell["printed_ci_upper"])
    if low <= -100 or high <= -100:
        raise ValueError("Invalid logarithm domain: AAPC endpoint must exceed -100")
    return low, high


def log_trend(percent: Decimal) -> float:
    """Natural-log change per calendar year; do not treat percent as fraction."""
    percent = Decimal(percent)
    if not percent.is_finite() or percent <= -100:
        raise ValueError("Invalid logarithm domain")
    return math.log1p(float(percent / Decimal(100)))


def load_data(path: Path = DATA_PATH) -> dict:
    """Load the complete frozen panel; incomplete inputs stop, never impute."""
    data = {}
    with Path(path).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        expected_fields = {"stage", "age_bin", *NUMERIC_FIELDS, *PROVENANCE_FIELDS}
        if reader.fieldnames is None or set(reader.fieldnames) != expected_fields:
            raise ValueError("Untestable input: missing or unexpected CSV columns")
        if len(reader.fieldnames) != len(expected_fields):
            raise ValueError("Untestable input: duplicate CSV headers")
        for row in reader:
            key = (row["stage"], row["age_bin"])
            if key in data or key[0] not in STAGES or key[1] not in AGES:
                raise ValueError(f"Duplicate or unexpected source cell: {key}")
            if None in row:
                raise ValueError(f"Malformed CSV row: {key}")
            for field in PROVENANCE_FIELDS:
                if not isinstance(row.get(field), str) or not row[field].strip():
                    raise ValueError(f"Missing provenance {field}: {key}")
            for field in NUMERIC_FIELDS:
                if not isinstance(row.get(field), str) or not re.fullmatch(
                    r"-?\d+\.\d{2}", row[field]
                ):
                    raise ValueError(f"Missing/invalid printed decimal {field}: {key}")
            for field, expected in FIXED_METADATA.items():
                if row[field] != expected:
                    raise ValueError(f"Inconsistent source metadata {field}: {key}")
            expected_row = key[0].capitalize()
            if key[0] == "unknown":
                expected_row += " (footnote a)"
            if row["source_row"] != expected_row:
                raise ValueError(f"Inconsistent source row: {key}")
            if row["source_column"] != SOURCE_COLUMNS[key[1]]:
                raise ValueError(f"Inconsistent source age column: {key}")
            expected_locator = (
                "Main article HTML > Table 2 > Stage at diagnosis > "
                f"{expected_row} > {SOURCE_COLUMNS[key[1]]} > AAPC (95% CI), %"
            )
            if row["locator"] != expected_locator:
                raise ValueError(f"Inconsistent source locator: {key}")
            for mode in MODES:
                percent_interval(row, mode)
            data[key] = row
    expected = {(stage, age) for stage in STAGES for age in AGES}
    if set(data) != expected:
        raise ValueError(f"Untestable input: missing source cells {sorted(expected - set(data))}")
    for field in PROVENANCE_FIELDS:
        if field not in {"source_row", "source_column", "locator"}:
            if len({row[field] for row in data.values()}) != 1:
                raise ValueError(f"Inconsistent panel metadata: {field}")
    return data


def summarize_intervals(intervals: dict) -> dict:
    """Exact decimal decision and sharp symmetric log-interval allowance."""
    if not intervals:
        raise ValueError("Empty interval panel")
    bounds = {}
    for age, (low, high) in intervals.items():
        low, high = Decimal(low), Decimal(high)
        if not low.is_finite() or not high.is_finite() or low > high or low <= -100:
            raise ValueError(f"Invalid interval for {age}")
        bounds[age] = (low, high)
    lower = max(pair[0] for pair in bounds.values())
    upper = min(pair[1] for pair in bounds.values())
    compatible = lower <= upper  # Exact decimal comparison; equality passes.
    lower_log, upper_log = log_trend(lower), log_trend(upper)
    gap = 0.0 if compatible else lower_log - upper_log
    return {
        "ages": list(bounds),
        "verdict": "PASS" if compatible else "FAIL",
        "lower_percent": str(lower),
        "upper_percent": str(upper),
        "lower_log": lower_log,
        "upper_log": upper_log,
        "lower_ages": [age for age, pair in bounds.items() if pair[0] == lower],
        "upper_ages": [age for age, pair in bounds.items() if pair[1] == upper],
        "gap": gap,
        "allowance": gap / 2,
        "reference": (lower_log + upper_log) / 2,
    }


def stage_summary(data: dict, stage: str, ages: tuple = AGES,
                  mode: str = "primary") -> dict:
    if stage not in STAGES or not ages or len(set(ages)) != len(ages):
        raise ValueError("Invalid stage or age panel")
    if any(age not in AGES for age in ages):
        raise ValueError("Unknown age group")
    return summarize_intervals({age: percent_interval(data[stage, age], mode)
                                for age in ages})


def analyze(data: dict | None = None) -> dict:
    if data is None:
        data = load_data()
    modes = {mode: {stage: stage_summary(data, stage, mode=mode)
                    for stage in STAGES} for mode in MODES}
    leave_one_out = {
        stage: {omitted: stage_summary(data, stage,
                                      tuple(age for age in AGES if age != omitted))
                for omitted in AGES}
        for stage in STAGES
    }
    cells = []
    for stage in STAGES:
        for age in AGES:
            row = data[stage, age]
            low, high = percent_interval(row)
            point_low, point_high = percent_interval(row, "point_rounding")
            cells.append({
                "stage": stage,
                "age_bin": age,
                **{field: row[field] for field in NUMERIC_FIELDS},
                "point_log": log_trend(Decimal(row["printed_aapc"])),
                "point_lower_log": log_trend(point_low),
                "point_upper_log": log_trend(point_high),
                "ci_lower_percent": str(low),
                "ci_upper_percent": str(high),
                "ci_lower_log": log_trend(low),
                "ci_upper_log": log_trend(high),
            })
    primary = modes["primary"]
    failed = [s for s in STAGES if primary[s]["verdict"] == "FAIL"]
    return {
        "row_count": len(data),
        "cells": cells,
        "modes": modes,
        "leave_one_out": leave_one_out,
        "overall": {
            "verdict": "FAIL" if failed else "PASS",
            "failed_stages": failed,
            "passing_stages": [s for s in STAGES if s not in failed],
            "allowance": max(primary[s]["allowance"] for s in STAGES),
        },
        "hashes": file_hashes(),
    }


def _number(value: float) -> str:
    return f"{value:.12f}"


def _names(values: list) -> str:
    return ", ".join(values) if values else "none"


def render_results(result: dict) -> str:
    """The complete numerical README section; every value comes from the CSV."""
    overall = result["overall"]
    lines = [
        "## Recomputed results",
        "",
        f"Panel: **{result['row_count']} stage-age cells**, with all selected stages and ages retained.",
        f"**Overall descriptive verdict: {overall['verdict']}.** "
        f"Failing stages: **{_names(overall['failed_stages'])}**. "
        f"Passing stages: **{_names(overall['passing_stages'])}**.",
        "",
        "The smallest allowance shared across stages (with a separate reference "
        f"trend per stage) is **{_number(overall['allowance'])} natural-log units per calendar year**.",
        "",
        "### Printed inputs and primary transformed intervals",
        "",
        "AAPCs and printed marginal intervals are percentages per year from "
        "Table 2 [1]. The log estimate and widened intervals are derived. "
        "Rounding is applied to percent endpoints before the logarithm. "
        "The point-rounding intervals, which omit sampling uncertainty, are also "
        "retained in the JSON output.",
        "",
        "| Stage | Age | Printed AAPC | Printed marginal interval | Log estimate | Widened log interval |",
        "|---|---|---:|---|---:|---|",
    ]
    for cell in result["cells"]:
        lines.append(
            f"| {cell['stage']} | {cell['age_bin']} | {cell['printed_aapc']} | "
            f"[{cell['printed_ci_lower']}, {cell['printed_ci_upper']}] | "
            f"{_number(cell['point_log'])} | "
            f"[{_number(cell['ci_lower_log'])}, {_number(cell['ci_upper_log'])}] |"
        )
    lines += [
        "",
        "### Primary compatibility and minimum allowance",
        "",
        "L is the greatest lower log bound; U is the least upper log bound. "
        "For PASS, [L,U] is the common-average compatibility set. For FAIL, "
        "L and U are separating bounds, not a reversed confidence interval. "
        "All log quantities have units per calendar year. Decisions use exact "
        "decimal percent endpoints, not the rounded display.",
        "",
        "| Stage | Verdict | L (attaining ages) | U (attaining ages) | Full positive gap | Symmetric allowance | Attaining reference |",
        "|---|---|---|---|---:|---:|---:|",
    ]
    for stage, summary in result["modes"]["primary"].items():
        lines.append(
            f"| {stage} | **{summary['verdict']}** | "
            f"{_number(summary['lower_log'])} ({_names(summary['lower_ages'])}) | "
            f"{_number(summary['upper_log'])} ({_names(summary['upper_ages'])}) | "
            f"{_number(summary['gap'])} | {_number(summary['allowance'])} | "
            f"{_number(summary['reference'])} |"
        )
    lines += [
        "",
        "The exact percent-scale endpoints below determine the same verdict "
        "because the logarithm is strictly increasing. They are the limiting "
        "widened source endpoints, not new uncertainty estimates.",
        "",
        "| Stage | Greatest lower endpoint (%) | Least upper endpoint (%) |",
        "|---|---:|---:|",
    ]
    for stage, summary in result["modes"]["primary"].items():
        lines.append(f"| {stage} | {summary['lower_percent']} | {summary['upper_percent']} |")
    lines += [
        "",
        "### Frozen rounding and uncertainty challenges",
        "",
        "Each cell shows verdict / symmetric log allowance. Printed CI uses no "
        "rounding widening; full-unit CI uses twice the primary half-unit; "
        "point rounding omits the published sampling uncertainty. These "
        "sensitivities do not replace the primary result.",
        "",
        "| Stage | Printed CI | Primary half-unit CI | Full-unit CI stress | Point rounding only |",
        "|---|---|---|---|---|",
    ]
    for stage in STAGES:
        entries = []
        for mode in ("printed_ci", "primary", "full_unit", "point_rounding"):
            summary = result["modes"][mode][stage]
            entries.append(f"{summary['verdict']} / {_number(summary['allowance'])}")
        lines.append(f"| {stage} | " + " | ".join(entries) + " |")
    lines += [
        "",
        "### Every leave-one-age-out diagnostic",
        "",
        "The omitted age is removed only in this diagnostic. Each row uses the "
        "primary rounding rule on the remaining pair. Passing rows remain visible.",
        "",
        "| Stage | Omitted age | Remaining ages | Verdict | L | U | Symmetric allowance |",
        "|---|---|---|---|---:|---:|---:|",
    ]
    for stage, omissions in result["leave_one_out"].items():
        for omitted, summary in omissions.items():
            lines.append(
                f"| {stage} | {omitted} | {_names(summary['ages'])} | {summary['verdict']} | "
                f"{_number(summary['lower_log'])} | {_number(summary['upper_log'])} | "
                f"{_number(summary['allowance'])} |"
            )
    lines += [
        "",
        "### Reproduction hashes",
        "",
        "SHA-256 of UTF-8 content after CRLF/CR-to-LF normalization, with "
        "repository-relative POSIX path keys. Hashes document exact inputs/code; "
        "they do not prove source truth or scientific validity.",
        "",
        "| Path | SHA-256 |",
        "|---|---|",
    ]
    lines.extend(f"| `{path}` | `{digest}` |" for path, digest in result["hashes"].items())
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("markdown", "json"), default="markdown")
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--check", action="store_true")
    action.add_argument("--write-readme", action="store_true")
    args = parser.parse_args()
    try:
        result = analyze()
        rendered = render_results(result)
        if args.check or args.write_readme:
            path = ROOT / "README.md"
            text = path.read_text(encoding="utf-8")
            if text.count(BEGIN) != 1 or text.count(END) != 1:
                raise ValueError("README needs exactly one generated-block marker pair")
            prefix, rest = text.split(BEGIN, 1)
            current, suffix = rest.split(END, 1)
            replacement = "\n\n" + rendered + "\n"
            if args.check:
                if current != replacement:
                    raise ValueError("README numerical section differs from CSV recomputation")
                print("README numerical section matches the aggregate CSV.")
            else:
                with path.open("w", encoding="utf-8", newline="\n") as stream:
                    stream.write(prefix + BEGIN + replacement + END + suffix)
                print("Recomputed README numerical section.")
        elif args.format == "json":
            print(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False))
        else:
            print(rendered, end="")
    except (ValueError, KeyError) as exc:
        parser.exit(2, f"Untestable input or reproduction error: {exc}\n")


if __name__ == "__main__":
    main()
