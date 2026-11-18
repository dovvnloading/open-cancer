#!/usr/bin/env python3
"""Nonclinical, unreviewed printed-IRR compatibility audit; offline stdlib only.

Do not use vaccinated person-years to estimate the source's adjusted IRRs.
The 95% marginal interval-box result is NOT a simultaneous statistical test.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "wu_figure3.csv"
REPORT = ROOT / "REPORT.md"
BINS = ("7-9", "10-12", "13-15")
AGES = ("10-16", "17+")
SOURCE_URL = "https://www.bmj.com/content/392/bmj-2025-087326.full"
LICENSE_URL = "https://creativecommons.org/licenses/by/4.0/"
BEGIN = "<!-- BEGIN RECOMPUTED RESULTS -->"
END = "<!-- END RECOMPUTED RESULTS -->"
COLUMNS = (
    "duration_bin", "age_at_first_dose", "adjusted_irr", "ci95_lower",
    "ci95_upper", "vaccinated_cases", "vaccinated_person_years", "authors",
    "publication_year", "journal", "doi", "source_url", "source_locator",
    "figure_locator", "measure", "interval_kind", "license", "license_url",
    "access_date",
)
FIXED = {
    "authors": "Wu et al.", "publication_year": "2026", "journal": "BMJ",
    "doi": "10.1136/bmj-2025-087326", "source_url": SOURCE_URL,
    "source_locator": "Results > Invasive cervical cancer by time since vaccination",
    "measure": "fully adjusted vaccinated versus unvaccinated IRR",
    "interval_kind": "published marginal 95% CI", "license": "CC BY 4.0",
    "license_url": LICENSE_URL, "access_date": "2026-11-18",
}
MODES = ("primary", "no_rounding", "full_unit", "point_only")


def normalized_sha256(path: Path) -> str:
    """Hash decoded UTF-8 with LF canonicalisation but unchanged EOF."""
    content = Path(path).read_bytes().decode("utf-8")
    canonical = content.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def file_hashes(root: Path = ROOT) -> dict[str, str]:
    names = (Path("PROTOCOL.md"), Path("data") / "wu_figure3.csv",
             Path("analysis.py"))
    prefix = Path("research") / "hpv-age-duration-factorisation"
    return {(prefix / name).as_posix(): normalized_sha256(root / name)
            for name in names}


def decimal_fraction(token: str) -> Fraction:
    if not isinstance(token, str) or not re.fullmatch(r"\d+\.\d{2}", token):
        raise ValueError(f"UNTESTABLE: malformed/missing printed two-decimal IRR: {token!r}")
    return Fraction(token)


def last_unit(token: str) -> Fraction:
    """The positive exact value of the token's last printed place."""
    decimal_fraction(token)
    return Fraction(1, 10 ** len(token.partition(".")[2]))


def load_data(path: Path = DATA) -> dict[tuple[str, str], dict[str, str]]:
    """Reject missing/extra, zero-event, invented and source-mismatched rows."""
    with Path(path).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != COLUMNS:
            raise ValueError("UNTESTABLE: source columns do not match frozen CSV schema")
        rows = list(reader)
    if len(rows) != 6:
        raise ValueError("UNTESTABLE: expected exactly six printed age-by-bin rows")
    indexed: dict[tuple[str, str], dict[str, str]] = {}
    expected_order = [(b, a) for b in BINS for a in AGES]
    for row in rows:
        if any(row.get(key) != value for key, value in FIXED.items()):
            raise ValueError("UNTESTABLE: source/rights metadata altered or missing")
        if any(not row.get(name) for name in COLUMNS):
            raise ValueError("UNTESTABLE: missing source cell or locator")
        key = (row["duration_bin"], row["age_at_first_dose"])
        if key not in expected_order or key in indexed:
            raise ValueError(f"UNTESTABLE: non-primary or duplicate row {key}")
        label = (f"Fig. 3 > {key[0]} years > vaccinated at "
                 f"{'10-16' if key[1] == '10-16' else '>=17'} years")
        if row["figure_locator"] != label:
            raise ValueError("UNTESTABLE: figure locator mismatch")
        point, low, high = (decimal_fraction(row[name]) for name in
                            ("adjusted_irr", "ci95_lower", "ci95_upper"))
        if not 0 < low <= point <= high:
            raise ValueError("UNTESTABLE: inconsistent/nonpositive marginal interval")
        if not re.fullmatch(r"[0-9]+", row["vaccinated_cases"]) or not re.fullmatch(
                r"[1-9][0-9]*", row["vaccinated_person_years"]):
            raise ValueError("UNTESTABLE: invalid count or person-years")
        if int(row["vaccinated_cases"]) == 0:
            raise ValueError("UNTESTABLE: primary zero-event cell cannot supply a fitted ratio")
        indexed[key] = row
    if list(indexed) != expected_order:
        raise ValueError("UNTESTABLE: incomplete or reordered primary panel")
    return indexed


def irr_range(row: dict[str, str], mode: str = "primary") -> tuple[Fraction, Fraction]:
    if mode not in MODES:
        raise ValueError(f"Unknown mode {mode}")
    point, lower, upper = (decimal_fraction(row[name]) for name in
                           ("adjusted_irr", "ci95_lower", "ci95_upper"))
    if mode == "point_only":
        q = last_unit(row["adjusted_irr"]) / 2
        result = (point - q, point + q)
    else:
        q = (Fraction(0) if mode == "no_rounding" else
             last_unit(row["ci95_lower"]) if mode == "full_unit" else
             last_unit(row["ci95_lower"]) / 2)
        # CI endpoints can in general have different printed precision.
        upper_q = (Fraction(0) if mode == "no_rounding" else
                   last_unit(row["ci95_upper"]) if mode == "full_unit" else
                   last_unit(row["ci95_upper"]) / 2)
        result = (lower - q, upper + upper_q)
    if result[0] <= 0 or result[0] > result[1]:
        raise ValueError("UNTESTABLE: nonpositive denominator/adjusted lower bound")
    return result


def ratio_ranges(rows: dict, mode: str = "primary",
                 exclude: str | None = None) -> dict[str, tuple[Fraction, Fraction]]:
    answer = {}
    for duration in BINS:
        if duration == exclude:
            continue
        e = irr_range(rows[(duration, AGES[0])], mode)
        l = irr_range(rows[(duration, AGES[1])], mode)
        answer[duration] = (e[0] / l[1], e[1] / l[0])
    if not answer:
        raise ValueError("UNTESTABLE: no comparable bin")
    return answer


def summarize(boxes: dict[str, tuple]) -> dict:
    """Works on Fraction boxes and on the approximate floating boxes."""
    left = max(interval[0] for interval in boxes.values())
    right = min(interval[1] for interval in boxes.values())
    if left <= 0 or right <= 0:
        raise ValueError("UNTESTABLE: nonpositive ratio bounds")
    delta = math.log(float(left)) - math.log(float(right))
    allowance = max(0.0, delta / 2)
    center = (math.log(float(left)) + math.log(float(right))) / 2
    return {
        "verdict": "CONSISTENT" if left <= right else "FALSIFIED",
        "lower": left, "upper": right,
        "limiting_lower": [key for key, band in boxes.items() if band[0] == left],
        "limiting_upper": [key for key, band in boxes.items() if band[1] == right],
        "allowance": allowance, "multiplicative_allowance": math.exp(allowance),
        "log_witness": center, "witness": math.exp(center),
        "gap_log": max(0.0, delta),
    }


def bonferroni_irr_range(row: dict[str, str]) -> tuple[float, float]:
    """Explicitly approximate, not published six-way simultaneous confidence."""
    p = float(decimal_fraction(row["adjusted_irr"]))
    a = float(decimal_fraction(row["ci95_lower"]))
    b = float(decimal_fraction(row["ci95_upper"]))
    point_q = float(last_unit(row["adjusted_irr"]) / 2)
    lower_q = float(last_unit(row["ci95_lower"]) / 2)
    upper_q = float(last_unit(row["ci95_upper"]) / 2)
    if min(p-point_q, a-lower_q) <= 0:
        raise ValueError("UNTESTABLE: Bonferroni log interval includes zero")
    radius = max(math.log((p+point_q)/(a-lower_q)),
                 math.log((b+upper_q)/(p-point_q)))
    normal = statistics.NormalDist()
    ratio = normal.inv_cdf(1-0.05/(2*6))/normal.inv_cdf(0.975)
    return ((p-point_q)*math.exp(-ratio*radius),
            (p+point_q)*math.exp(ratio*radius))


def bonferroni_ratio_ranges(rows: dict) -> dict[str, tuple[float, float]]:
    result = {}
    for duration in BINS:
        e = bonferroni_irr_range(rows[(duration, AGES[0])])
        l = bonferroni_irr_range(rows[(duration, AGES[1])])
        result[duration] = (e[0]/l[1], e[1]/l[0])
    return result


def poisson_cdf(n: int, mean: float) -> float:
    """Analytic Poisson cumulative probability for a small integer count."""
    if n < 0:
        return 0.0
    if mean < 0 or not math.isfinite(mean):
        raise ValueError("invalid Poisson mean")
    term = math.exp(-mean)
    total = term
    for j in range(1, n+1):
        term *= mean/j
        total += term
    return min(1.0, max(0.0, total))


def poisson_mean_interval(count: int, alpha: float = 0.05) -> tuple[float, float]:
    """Two-sided exact-tail Poisson (Garwood) mean interval, root-solved."""
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        raise ValueError("invalid case count")
    if not 0 < alpha < 1:
        raise ValueError("invalid alpha")
    tail = alpha/2
    if count == 0:
        return (0.0, -math.log(tail))

    def root(target, cdf_index, increasing):
        lo, hi = 0.0, float(max(2, count))
        criterion = (lambda x: 1-poisson_cdf(cdf_index, x)) if increasing \
                    else (lambda x: poisson_cdf(cdf_index, x))
        while (criterion(hi) < target if increasing else criterion(hi) > target):
            hi *= 2
            if hi > 1_000_000:
                raise ValueError("Poisson inversion failed to bracket")
        for _ in range(100):
            mid = (lo+hi)/2
            if (criterion(mid) < target) == increasing:
                lo = mid
            else:
                hi = mid
        return (lo+hi)/2

    return root(tail, count-1, True), root(tail, count, False)


def exact_count_diagnostics(rows: dict) -> list[dict]:
    result = []
    for b in BINS:
        for a in AGES:
            row = rows[(b, a)]
            count = int(row["vaccinated_cases"])
            years = int(row["vaccinated_person_years"])
            left, right = poisson_mean_interval(count)
            result.append({"bin": b, "age": a, "cases": count, "years": years,
                           "crude_rate_interval_per_100k":
                               [100000*left/years, 100000*right/years]})
    return result


def analyze(rows: dict | None = None, include_hashes: bool = True) -> dict:
    if rows is None:
        rows = load_data()
    summaries = {}
    bands = {}
    for mode in MODES:
        bands[mode] = ratio_ranges(rows, mode)
        summaries[mode] = summarize(bands[mode])
    bands["bonferroni_approx"] = bonferroni_ratio_ranges(rows)
    summaries["bonferroni_approx"] = summarize(bands["bonferroni_approx"])
    omitted = {b: summarize(ratio_ranges(rows, "primary", exclude=b)) for b in BINS}
    # Exact point ratios are reported, but have no inferential status.
    point_ratios = {b: decimal_fraction(rows[(b, AGES[0])]["adjusted_irr"]) /
                       decimal_fraction(rows[(b, AGES[1])]["adjusted_irr"])
                    for b in BINS}
    return {"boxes": bands, "summaries": summaries, "omissions": omitted,
            "point_ratios": point_ratios,
            "counts": exact_count_diagnostics(rows),
            "zero_bin_adjusted_ratio": "UNTESTABLE: [0,+infinity]; no published early adjusted IRR",
            "inferential_verdict": "INCONCLUSIVE: marginal 95% intervals are not simultaneous",
            "hashes": file_hashes() if include_hashes else {}}


def show_fraction(value) -> str:
    if isinstance(value, Fraction):
        return f"{value.numerator}/{value.denominator}"
    return "approximate (floating)"


def fmt(value) -> str:
    return f"{float(value):.12f}"


def describe_summary(name: str, summary: dict) -> str:
    lower, upper = summary["lower"], summary["upper"]
    return (f"| {name} | **{summary['verdict']}** | {fmt(lower)} | "
            f"{fmt(upper)} | {fmt(summary['allowance'])} | "
            f"{','.join(summary['limiting_lower'])} / "
            f"{','.join(summary['limiting_upper'])} |")


def format_results(result: dict) -> str:
    lines = [
        "### Primary exact-rational rounded-interval projection", "",
        "Each ratio is the source's adjusted IRR for dose age 10-16 divided by "
        "its adjusted IRR for dose age 17+ **in the same duration bin**. "
        "Ratios below are projections of **marginal** intervals, not joint "
        "95% confidence intervals or direct vaccinated-only comparisons.", "",
        "| Years since first dose | Printed point ratio | Exact point fraction | "
        "Primary lower | Primary upper | Exact rounded lower / upper |",
        "| --- | ---: | --- | ---: | ---: | --- |",
    ]
    for duration in BINS:
        point = result["point_ratios"][duration]
        lo, hi = result["boxes"]["primary"][duration]
        lines.append(f"| {duration} | {fmt(point)} | {show_fraction(point)} | "
                     f"{fmt(lo)} | {fmt(hi)} | "
                     f"{show_fraction(lo)} / {show_fraction(hi)} |")
    lines += ["", "### Frozen intersection and sensitivity results", "",
              "The two middle bounds are max(lower) and min(upper), not a "
              "valid ordered interval when max(lower) exceeds min(upper).", "",
              "| Prespecified variant | Descriptive verdict | max lower | min upper | "
              "Sharp symmetric log allowance | Limiting lower / upper bins |",
              "| --- | --- | ---: | ---: | ---: | --- |"]
    for key, name in (("primary", "Primary; 0.005-rounded endpoints"),
                      ("no_rounding", "Printed CI; no widening"),
                      ("full_unit", "Full 0.01-unit widening"),
                      ("point_only", "Rounded points only (not CI)"),
                      ("bonferroni_approx", "Approximate six-IRR Bonferroni/Wald")):
        lines.append(describe_summary(name, result["summaries"][key]))
    for b in BINS:
        lines.append(describe_summary(f"Primary minus {b} years", result["omissions"][b]))
    primary = result["summaries"]["primary"]
    point = result["summaries"]["point_only"]
    lines += ["", f"**Primary descriptive verdict: {primary['verdict']}.** "
              f"Common q in [{fmt(primary['lower'])}, {fmt(primary['upper'])}]. "
              f"Exact intersection extrema: {show_fraction(primary['lower'])} "
              f"through {show_fraction(primary['upper'])}. "
              f"Minimum symmetric log allowance epsilon*={fmt(primary['allowance'])}; "
              f"multiplicative allowance exp(epsilon*)="
              f"{fmt(primary['multiplicative_allowance'])}. "
              f"Attaining log midpoint={fmt(primary['log_witness'])} "
              f"(q={fmt(primary['witness'])}).", "",
              f"**Point-only sensitivity: {point['verdict']}**, "
              f"epsilon*={fmt(point['allowance'])}, "
              f"exp(epsilon*)={fmt(point['multiplicative_allowance'])}, "
              f"attaining log reference={fmt(point['log_witness'])}, "
              f"q={fmt(point['witness'])}; "
              "this is a point-rounding diagnostic and does not override the "
              "marginal-CI primary rule.", "",
              "The approximate Bonferroni/Wald variant expands six source "
              "95% marginal log-CI widths using a normal quantile at alpha/6. "
              "Because the actual IRR covariance and underlying interval "
              "method are unavailable, this is **not** an established "
              "95% family-wise test. Statistical inference: **INCONCLUSIVE**.", "",
              "### Descriptive exact-count Poisson identification intervals", "",
              "These invert exact Poisson tails conditional on each **vaccinated** "
              "person-year denominator. They describe crude rates only; "
              "they never reconstruct an adjusted IRR or its standard error.", "",
              "| Years | Dose age | Vaccinated cases / person-years | "
              "95% exact-count crude rate per 100,000 PY |",
              "| --- | --- | ---: | --- |"]
    for c in result["counts"]:
        a, b = c["crude_rate_interval_per_100k"]
        lines.append(f"| {c['bin']} | {c['age']} | {c['cases']} / "
                     f"{c['years']} | [{fmt(a)}, {fmt(b)}] |")
    lines += ["", "For excluded early zero-case cells (1-3, 4-6, "
              "16-18 years), a two-sided exact-count Poisson bound has "
              "mean lower 0 and upper -ln(0.025), conditional on known "
              "exposure. Their **adjusted age-ratio identification set is "
              "[0,+infinity] / UNTESTABLE**, because the early adjusted "
              "IRR and common unvaccinated adjustment/reference are missing. "
              "No continuity corrections or invented point ratios are used.", "",
              "### Normalized reproducibility hashes", "",
              "UTF-8 CRLF and CR are converted to LF before SHA-256; keys use "
              "Path.as_posix(). Hashes document text identity, not source truth.", "",
              "| Input path | SHA-256 |", "| --- | --- |"]
    for path, digest in result["hashes"].items():
        lines.append(f"| {path} | {digest} |")
    return "\n".join(lines)


def json_ready(value):
    if isinstance(value, Fraction):
        return {"numerator": value.numerator, "denominator": value.denominator,
                "decimal": float(value)}
    if isinstance(value, dict):
        return {key: json_ready(v) for key, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_ready(v) for v in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    exclusive = parser.add_mutually_exclusive_group()
    exclusive.add_argument("--check", action="store_true")
    exclusive.add_argument("--write-report", action="store_true")
    exclusive.add_argument("--json", action="store_true")
    args = parser.parse_args()
    try:
        result = analyze()
        block = format_results(result)
        if args.json:
            print(json.dumps(json_ready(result), sort_keys=True, indent=2, allow_nan=False))
        elif args.check or args.write_report:
            source = REPORT.read_text(encoding="utf-8")
            if source.count(BEGIN) != 1 or source.count(END) != 1:
                raise ValueError("REPORT requires one generated block marker pair")
            pre, rest = source.split(BEGIN, 1)
            current, post = rest.split(END, 1)
            if args.check:
                if current.strip() != block.strip():
                    raise ValueError("REPORT computed block is stale")
                print("Nonclinical HPV report matches frozen inputs and script.")
            else:
                REPORT.write_text(pre + BEGIN + "\n\n" + block + "\n\n" + END + post,
                                  encoding="utf-8", newline="\n")
                print("Updated HPV report numerical section.")
        else:
            print(block)
        return 0
    except (ValueError, OSError, KeyError, OverflowError) as error:
        print(f"HPV factorisation UNTESTABLE/reproduction failure: {error}",
              file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
