#!/usr/bin/env python3
"""Recompute the frozen aggregate analysis. Nonclinical; no medical advice.

Only Python's standard library and the attributed Table 1 CSV are used.
Exact conditional Poisson inference means exact-distribution inversion, with
floating-point evaluation (not interval-certified arithmetic).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parents[1]
DATA_PATH = ROOT / "data" / "turbow_table1.csv"
AGES = ("65-69", "70-74", "75-79", "80-84", "85+")
PRIMARY_AGES = AGES[:4]
STAGES = ("overall", "localized", "regional", "distant", "unknown")
NAMED_STAGES = STAGES[1:4]
CONTRASTS = tuple(zip(PRIMARY_AGES, PRIMARY_AGES[1:]))
ALPHA = 0.05
FAMILY_SIZE = len(CONTRASTS) * len(STAGES)
BEGIN = "<!-- BEGIN RECOMPUTED RESULTS -->"
END = "<!-- END RECOMPUTED RESULTS -->"
PROVENANCE_FIELDS = (
    "source", "doi", "source_url", "table", "source_row", "source_column",
    "population_locator", "access_date", "access_timezone", "rate_units",
    "rate_adjustment", "dataset_version", "reuse_basis", "reuse_url",
)


def normalized_sha256(path: Path) -> str:
    """Hash UTF-8 text after universal-newline normalization, preserving EOF."""
    text = Path(path).read_text(encoding="utf-8")
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def file_hashes(root: Path = ROOT) -> dict[str, str]:
    prefix = Path("research") / "breast-incidence-decline-stage"
    paths = (Path("PROTOCOL.md"), Path("data") / "turbow_table1.csv", Path("analysis.py"))
    return {(prefix / path).as_posix(): normalized_sha256(root / path) for path in paths}


def load_data(path: Path = DATA_PATH) -> dict:
    data = {}
    with Path(path).open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            key = (row["age_bin"], row["stage"])
            if key in data or key[0] not in AGES or key[1] not in STAGES:
                raise ValueError(f"Duplicate or unexpected source cell: {key}")
            for field in PROVENANCE_FIELDS:
                if not row.get(field, "").strip():
                    raise ValueError(f"Missing provenance {field}: {key}")
            for field in ("cases", "population"):
                if not re.fullmatch(r"[0-9]+", row[field]):
                    raise ValueError(f"Missing/noninteger {field}: {key}; contrast untestable")
                row[field] = int(row[field])
            if row["population"] <= 0:
                raise ValueError(f"Nonpositive exposure: {key}; contrast untestable")
            for field in ("printed_rate", "printed_ci_lower", "printed_ci_upper"):
                if not re.fullmatch(r"[0-9]+\.[0-9]", row[field]):
                    raise ValueError(f"Invalid printed decimal {field}: {key}")
            if not (Fraction(row["printed_ci_lower"]) <= Fraction(row["printed_rate"])
                    <= Fraction(row["printed_ci_upper"])):
                raise ValueError(f"Invalid printed interval: {key}")
            data[key] = row
    expected = {(age, stage) for age in AGES for stage in STAGES}
    if set(data) != expected:
        raise ValueError(f"Missing source cells: {sorted(expected - set(data))}; untestable")
    for age in AGES:
        if len({data[age, stage]["population"] for stage in STAGES}) != 1:
            raise ValueError(f"Inconsistent population within age {age}")
    return data


def rate(cell: dict) -> Fraction:
    return Fraction(100000 * cell["cases"], cell["population"])


def residuals(data: dict) -> dict[str, int]:
    return {age: data[age, "overall"]["cases"]
            - sum(data[age, stage]["cases"] for stage in STAGES[1:]) for age in AGES}


def _beta_fraction(a: float, b: float, x: float) -> float:
    """Modified Lentz continued fraction for the incomplete beta function."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 10001):
        twice = 2 * m
        aa = m * (b - m) * x / ((qam + twice) * (a + twice))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + twice) * (qap + twice))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 3e-14:
            return h
    raise ArithmeticError("Incomplete-beta continued fraction did not converge")


def beta_cdf(x: float, a: float, b: float) -> float:
    if a <= 0 or b <= 0 or not 0 <= x <= 1:
        raise ValueError("Invalid beta parameters")
    if x == 0 or x == 1:
        return x
    log_front = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                 + a * math.log(x) + b * math.log1p(-x))
    front = math.exp(log_front)
    if x < (a + 1.0) / (a + b + 2.0):
        value = front * _beta_fraction(a, b, x) / a
    else:
        value = 1.0 - front * _beta_fraction(b, a, 1.0 - x) / b
    return min(1.0, max(0.0, value))


def beta_quantile(q: float, a: float, b: float) -> float:
    if not 0 < q < 1 or a <= 0 or b <= 0:
        raise ValueError("Invalid beta quantile parameters")
    if q > 0.5:
        return 1.0 - beta_quantile(1.0 - q, b, a)
    lo, hi = 0.0, 1.0
    for _ in range(70):
        mid = (lo + hi) / 2.0
        if mid == lo or mid == hi:
            break
        if beta_cdf(mid, a, b) < q:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def binomial_tails(k: int, n: int, p: float) -> tuple[float, float]:
    """Return P(X<=k), P(X>=k), each evaluated without a small-tail subtraction."""
    if not 0 <= k <= n or not 0 < p < 1:
        raise ValueError("Invalid binomial parameters")
    lower = 1.0 if k == n else beta_cdf(1.0 - p, n - k, k + 1)
    upper = 1.0 if k == 0 else beta_cdf(p, k, n - k + 1)
    return lower, upper


def conditional_rr(younger: int, older: int, exposure_younger: int,
                   exposure_older: int, alpha: float = ALPHA,
                   family_size: int = FAMILY_SIZE) -> dict:
    for count in (younger, older):
        if not isinstance(count, int) or isinstance(count, bool) or count < 0:
            raise ValueError("Counts must be nonnegative integers")
    if exposure_younger <= 0 or exposure_older <= 0:
        raise ValueError("Exposures must be positive")
    if not 0 < alpha < 1 or not isinstance(family_size, int) or family_size < 1:
        raise ValueError("Invalid error allocation")
    tail = alpha / (2 * family_size)
    low_p = 0.0 if older == 0 else beta_quantile(tail, older, younger + 1)
    high_p = 1.0 if younger == 0 else beta_quantile(1.0 - tail, older + 1, younger)
    one_p = 1.0 if younger == 0 else beta_quantile(1.0 - 2 * tail, older + 1, younger)

    def transform(p):
        return math.inf if p == 1.0 else (exposure_younger / exposure_older) * p / (1 - p)

    lower, upper, one_upper = transform(low_p), transform(high_p), transform(one_p)
    rr = (older * exposure_younger / (younger * exposure_older)
          if younger else (math.inf if older else None))
    p0 = exposure_older / (exposure_younger + exposure_older)
    left, right = binomial_tails(older, younger + older, p0)
    p_two = min(1.0, 2 * min(left, right))
    return {"rr": rr, "lower": lower, "upper": upper,
            "one_sided_upper": one_upper, "p_one_sided": left,
            "p_two_sided": p_two, "adjusted_two_sided_p": min(1.0, family_size * p_two),
            "decision": "decline" if upper < 1 else "increase" if lower > 1 else "unresolved",
            "one_sided_decline": one_upper < 1}


def analyze(data: dict | None = None) -> dict:
    data = load_data() if data is None else data
    rem = residuals(data)
    audit = []
    for age in AGES:
        for stage in STAGES:
            cell = data[age, stage]
            raw, printed = rate(cell), Fraction(cell["printed_rate"])
            audit.append({"age": age, "stage": stage, "cases": cell["cases"],
                          "population": cell["population"], "computed_rate": float(raw),
                          "printed_rate": float(printed), "rate_difference": float(raw - printed),
                          "printed_ci_lower": float(cell["printed_ci_lower"]),
                          "printed_ci_upper": float(cell["printed_ci_upper"]),
                          "rounding_matches": abs(raw - printed) <= Fraction(1, 20),
                          "printed_ci_contains_rate": Fraction(cell["printed_ci_lower"]) <= raw
                          <= Fraction(cell["printed_ci_upper"])})
    accounting = [{"age": age, "overall": data[age, "overall"]["cases"],
                   "stage_sum": sum(data[age, s]["cases"] for s in STAGES[1:]),
                   "residual": rem[age], "unknown_pool": data[age, "unknown"]["cases"] + rem[age]}
                  for age in AGES]
    contrasts, decomposition, allocation = [], [], []
    for young, old in CONTRASTS:
        ey, eo = (data[a, "overall"]["population"] for a in (young, old))
        for stage in STAGES:
            y, o = data[young, stage], data[old, stage]
            contrast = conditional_rr(y["cases"], o["cases"], ey, eo)
            contrast.update(younger=young, older=old, stage=stage,
                            rate_younger=float(rate(y)), rate_older=float(rate(o)),
                            change=float(rate(o) - rate(y)))
            contrasts.append(contrast)
        losses = {stage: rate(data[young, stage]) - rate(data[old, stage]) for stage in STAGES}
        losses["residual"] = Fraction(100000 * rem[young], ey) - Fraction(100000 * rem[old], eo)
        if sum(losses[s] for s in STAGES[1:] + ("residual",)) != losses["overall"]:
            raise ArithmeticError("Stage decomposition failed its exact identity")
        total = losses["overall"]
        shares = {s: float(losses[s] / total) for s in STAGES[1:] + ("residual",)} if total > 0 else None
        pool_y = data[young, "unknown"]["cases"] + rem[young]
        pool_o = data[old, "unknown"]["cases"] + rem[old]
        valid = rem[young] >= 0 and rem[old] >= 0
        share_min = ((losses["localized"] - Fraction(100000 * pool_o, eo)) / total
                     if valid and total > 0 else None)
        share_max = ((losses["localized"] + Fraction(100000 * pool_y, ey)) / total
                     if valid and total > 0 else None)
        robust_total = next(c for c in contrasts if c["younger"] == young and c["stage"] == "overall")["decision"] == "decline"
        decomposition.append({"younger": young, "older": old,
                              "losses": {s: float(v) for s, v in losses.items()}, "shares": shares,
                              "localized_share_min": float(share_min) if share_min is not None else None,
                              "localized_share_max": float(share_max) if share_max is not None else None,
                              "overall_decline_robust": robust_total,
                              "recorded_majority": bool(robust_total and shares and shares["localized"] > 0.5),
                              "allocation_majority": bool(robust_total and share_min is not None and share_min > Fraction(1, 2))})
        for stage in NAMED_STAGES:
            entry = {"younger": young, "older": old, "stage": stage, "allocation_valid": valid}
            if not valid:
                entry.update(lower=None, upper=None, rr_min=None, rr_max=None, decline_survives=None)
            else:
                y, o = data[young, stage]["cases"], data[old, stage]["cases"]
                low = conditional_rr(y + pool_y, o, ey, eo)
                high = conditional_rr(y, o + pool_o, ey, eo)
                entry.update(lower=low["lower"], upper=high["upper"], rr_min=low["rr"],
                             rr_max=high["rr"], decline_survives=high["upper"] < 1,
                             younger_pool=pool_y, older_pool=pool_o)
            allocation.append(entry)
    one_changes = [c for c in contrasts if (c["decision"] == "decline") != c["one_sided_decline"]]
    proportions = [{"age": age,
                    "distant_percent": 100 * data[age, "distant"]["cases"] / data[age, "overall"]["cases"],
                    "unknown_percent": 100 * data[age, "unknown"]["cases"] / data[age, "overall"]["cases"]}
                   for age in AGES]
    return {"reconciliation": audit, "accounting": accounting, "contrasts": contrasts,
            "decomposition": decomposition, "allocation": allocation,
            "one_sided_changes": one_changes, "stage_proportions": proportions,
            "hashes": file_hashes()}


def _table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "| " + " | ".join("---" for _ in headers) + " |"]
                     + ["| " + " | ".join(str(v) for v in row) + " |" for row in rows])


def _edge(row):
    return row["younger"] + " -> " + row["older"]


def render_results(report: dict) -> str:
    """The entire numerical results section, including hashes, comes from this."""
    audit = report["reconciliation"]
    comparisons = report["contrasts"]
    matches = sum(r["rounding_matches"] for r in audit)
    contained = sum(r["printed_ci_contains_rate"] for r in audit)
    largest = max(audit, key=lambda r: abs(r["rate_difference"]))
    sections = ["## Recomputed results",
                f"The CSV contains **{len(audit)} cells** across **{len(AGES)} age groups**. "
                f"**{matches}/{len(audit)}** count/population rates agree with the printed rates "
                f"within half a last digit; **{contained}/{len(audit)}** also lie within their "
                "printed intervals. These checks establish arithmetic compatibility, not source truth. "
                f"The largest absolute rate discrepancy is **{abs(largest['rate_difference']):.6f}** "
                f"per 100,000 ({largest['age']}, {largest['stage']}).",
                "### Every transcribed count/rate reconciliation",
                "Rates are cases per 100,000 women. The source interval is its printed marginal "
                "95% interval, not our simultaneous contrast interval. Gap = count/population rate "
                "minus printed rate. The CSV retains the original decimal strings and provenance."]
    rows = []
    for r in audit:
        rows.append([r["age"], r["stage"], f"{r['cases']:,}", f"{r['population']:,}",
                     f"{r['computed_rate']:.6f}", f"{r['printed_rate']:.1f}",
                     f"[{r['printed_ci_lower']:.1f}, {r['printed_ci_upper']:.1f}]",
                     f"{r['rate_difference']:+.6f}", "yes" if r["rounding_matches"] else "NO"])
    sections.append(_table(["Age", "Endpoint", "Cases", "Population", "Computed rate", "Printed rate",
                            "Printed interval", "Gap", "Rounding match"], rows))
    sections.extend(["### Unassigned residuals",
                     "Residual = overall minus the four printed stage counts. The residual is kept "
                     "separate from unknown stage. The possible-allocation pool is unknown plus "
                     "residual, conditional on these being assignable cases.",
                     _table(["Age", "Overall", "Stage sum", "Residual", "Unknown + residual"],
                            [[r["age"], r["overall"], r["stage_sum"], r["residual"], r["unknown_pool"]]
                             for r in report["accounting"]]),
                     "### Primary simultaneous adjacent-age inference",
                     f"There are **{FAMILY_SIZE}** frozen comparisons. The exact conditional Poisson "
                     f"RR intervals use error **{ALPHA/FAMILY_SIZE:.10f}** per comparison "
                     f"(**{ALPHA/(2*FAMILY_SIZE):.10f}** in each tail), giving at least "
                     f"**{100*(1-ALPHA):.0f}%** simultaneous coverage within this family. "
                     "RR is older/younger; change is older minus younger per 100,000. "
                     "Decisions use unrounded endpoints. 'Unresolved' does not mean equal rates.",
                     _table(["Age contrast", "Endpoint", "Rate change", "RR", "Simultaneous RR interval",
                             "Adjusted two-sided p", "Primary direction"],
                            [[_edge(c), c["stage"], f"{c['change']:+.6f}", f"{c['rr']:.6f}",
                              f"[{c['lower']:.6f}, {c['upper']:.6f}]", f"{c['adjusted_two_sided_p']:.6g}",
                              c["decision"]] for c in comparisons])])
    total_declines = [c for c in comparisons if c["stage"] == "overall" and c["decision"] == "decline"]
    local_declines = [c for c in comparisons if c["stage"] == "localized" and c["decision"] == "decline"]
    advanced_declines = [c for c in comparisons if c["stage"] in ("regional", "distant") and c["decision"] == "decline"]
    sections.append(f"**Primary verdict:** {len(total_declines)} overall declines survive: "
                    + ("; ".join(_edge(c) for c in total_declines) or "none")
                    + f". Localized disease has {len(local_declines)} robust declines on "
                    + ("; ".join(_edge(c) for c in local_declines) or "no edges")
                    + f". Regional/distant disease has {len(advanced_declines)} robust declines "
                    "under this two-sided simultaneous rule. This last finding is limited evidence, "
                    "not proof of no regional or distant decline.")
    sections.extend(["### Signed decomposition of the observed change",
                     "Loss = younger minus older rate. Positive values contribute to decline; negative "
                     "values offset it. Residuals stay explicit. Shares divide by the net overall loss "
                     "only when it is positive. These are arithmetic contributions, not causal fractions "
                     "or confidence bounds.",
                     _table(["Age contrast", "Overall loss", "Localized", "Regional", "Distant", "Unknown", "Residual"],
                            [[_edge(d)] + [f"{d['losses'][s]:+.6f}" for s in STAGES + ("residual",)]
                             for d in report["decomposition"]]),
                     _table(["Declining age contrast", "Localized share", "Regional share", "Distant share",
                             "Unknown share", "Residual share", "Mainly localized (recorded)"],
                            [[_edge(d)] + [f"{100*d['shares'][s]:+.3f}%" for s in STAGES[1:] + ("residual",)]
                             + ["yes" if d["recorded_majority"] else "no"]
                             for d in report["decomposition"] if d["shares"] is not None]),
                     "### Adverse allocation of unknown stage and residuals",
                     "For each named stage separately, assign all older unknown/residual cases to it "
                     "and none younger to maximize the RR; reverse for the minimum. The interval "
                     "envelope also includes the frozen sampling uncertainty. These are conditional "
                     "completion bounds, using the same conservative multiplicity factor in a separate "
                     "sensitivity family. They do not presume all adverse stage allocations occur together.",
                     _table(["Age contrast", "Stage", "Possible point RR range", "Simultaneous RR envelope", "Decline survives"],
                            [[_edge(a), a["stage"], f"[{a['rr_min']:.6f}, {a['rr_max']:.6f}]" if a["allocation_valid"] else "untestable",
                              f"[{a['lower']:.6f}, {a['upper']:.6f}]" if a["allocation_valid"] else "untestable",
                              "yes" if a["decline_survives"] else "no" if a["allocation_valid"] else "untestable"]
                             for a in report["allocation"]]),
                     _table(["Declining age contrast", "Localized share across allocations", "Majority survives allocation"],
                            [[_edge(d), f"[{100*d['localized_share_min']:.3f}%, {100*d['localized_share_max']:.3f}%]"
                              if d["localized_share_min"] is not None else "untestable",
                              "yes" if d["allocation_majority"] else "no"]
                             for d in report["decomposition"] if d["shares"] is not None]),
                     "The share ranges condition on the observed totals and quantify allocation "
                     "uncertainty only. They are **not sampling confidence intervals** for the "
                     "majority claim. The RR envelopes, separately, incorporate the assumed Poisson variability.",
                     "### Prespecified one-sided sensitivity",
                     "These decline-only upper limits spend the full per-comparison error in one "
                     "tail. They do not replace the primary two-sided results. Every primary comparison "
                     "is shown so the sensitivity does not select favorable endpoints.",
                     _table(["Age contrast", "Endpoint", "One-sided simultaneous RR upper", "Decline-only decision"],
                            [[_edge(c), c["stage"], f"{c['one_sided_upper']:.6f}",
                              "decline" if c["one_sided_decline"] else "unresolved"] for c in comparisons]),
                     f"Changed decline decisions: **{len(report['one_sided_changes'])}**. "
                     + ("; ".join(_edge(c) + ", " + c["stage"] for c in report["one_sided_changes"]) or "None.")])
    sections.extend(["### Open-ended age group and source-prose check",
                     "The 85+ cells are reported descriptively in the reconciliation, residual, "
                     "and proportion tables; "
                     "no hypothesis test or finite-bin midpoint is assigned to them. "
                     "Recomputing stage proportions exposes the article's prose/category conflict: "
                     "its sequence described as distant-stage percentages matches the unknown-stage "
                     "count shares after rounding. Here both are calculated from the transcribed counts.",
                     _table(["Age", "Distant / overall", "Unknown / overall"],
                            [[r["age"], f"{r['distant_percent']:.3f}%", f"{r['unknown_percent']:.3f}%"]
                             for r in report["stage_proportions"]]),
                     "### Newline-normalized source and code hashes",
                     "SHA-256 uses UTF-8 after CRLF/CR -> LF normalization; path keys use "
                     "`Path.as_posix()`. The README is excluded to avoid self-reference. "
                     "A hash pins content, not source accuracy or scientific validity.",
                     _table(["Repository path", "SHA-256"], [[k, f"`{v}`"] for k, v in report["hashes"].items()])])
    return "\n\n".join(sections) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Verify the committed generated README block")
    parser.add_argument("--write-readme", action="store_true", help="Update only the generated README block")
    parser.add_argument("--format", choices=("markdown", "json"), default="markdown")
    args = parser.parse_args()
    if args.check and args.write_readme:
        parser.error("Choose --check or --write-readme")
    report = analyze()
    output = render_results(report)
    if args.check or args.write_readme:
        readme = ROOT / "README.md"
        text = readme.read_text(encoding="utf-8")
        if text.count(BEGIN) != 1 or text.count(END) != 1:
            raise ValueError("Expected exactly one pair of generated-result markers")
        before, remaining = text.split(BEGIN)
        _, after = remaining.split(END)
        expected = before + BEGIN + "\n\n" + output + "\n" + END + after
        if args.check:
            if text != expected:
                raise SystemExit("README generated results are stale")
            print("Generated README results and normalized hashes match.")
        else:
            readme.write_text(expected, encoding="utf-8", newline="\n")
    else:
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) if args.format == "json" else output, end="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
