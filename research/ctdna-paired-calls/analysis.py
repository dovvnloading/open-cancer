#!/usr/bin/env python3
"""Nonclinical, offline paired-call compatibility audit; standard library only.

No treatment, diagnosis or surveillance advice. The input consists only of
attributed published aggregate counts. This program does not fetch sources,
write output files, verify individual availability, or calibrate assay errors.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import sys


DIRECTORY = Path(__file__).resolve().parent
REPOSITORY = DIRECTORY.parents[1]
DATA = DIRECTORY / "data" / "aggregates.csv"
PROTOCOL = DIRECTORY / "PROTOCOL.md"
README = DIRECTORY / "README.md"
BEGIN = "<!-- BEGIN RECOMPUTED RESULTS -->"
END = "<!-- END RECOMPUTED RESULTS -->"
FAMILY_ALPHA = 0.05
ALPHA = 0.025
PRIMARY_CAP = Fraction(1, 100)
CAPS = tuple(Fraction(x) for x in ("0", "0.005", "0.01", "0.02", "0.03", "0.05"))
REVERSAL_CAPS = tuple(Fraction(x) for x in ("0", "0.005", "0.01", "0.02", "0.03"))
CELLS = ("++", "+-", "-+", "--")
FIELDS = (
    "id", "count", "population", "unit", "source_url", "doi", "locator",
    "access_date", "source_version", "article_rights", "rights_url",
    "reuse_basis", "reuse_basis_url", "notes",
)
IDS = {
    "guided_itt", "successful_analysis", "both_positive", "week4_only",
    "week7_only", "any_positive", "reported_negative", "itt_excluded",
    "itt_missing_week7", "itt_ineligible", "itt_withdrew",
}
COMPATIBLE = "compatible with stated latitude"
FALSIFIED = "falsified under stated caps"
UNTESTABLE = "untestable"


def normalized_sha256(path: Path) -> str:
    """SHA-256 after UTF-8 validation and CRLF/CR to LF normalization."""
    content = Path(path).read_bytes().decode("utf-8")
    normalized = content.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()


def _validate_counts(counts: dict[str, int]) -> None:
    if set(counts) != IDS:
        raise ValueError("Aggregate IDs differ from the frozen input panel")
    if any(type(value) is not int or value < 0 for value in counts.values()):
        raise ValueError("Counts must be nonnegative integers")
    if not 0 < counts["successful_analysis"] <= counts["guided_itt"]:
        raise ValueError("Successful-analysis denominator must be positive and within ITT")
    if counts["both_positive"] + counts["week4_only"] + counts["week7_only"] != counts["any_positive"]:
        raise ValueError("Printed positive patterns do not partition the positive category")
    if counts["any_positive"] + counts["reported_negative"] != counts["successful_analysis"]:
        raise ValueError("Printed positive/negative categories do not partition successful analyses")
    if sum(counts[key] for key in ("itt_missing_week7", "itt_ineligible", "itt_withdrew")) != counts["itt_excluded"]:
        raise ValueError("Figure-caption exclusions do not reconcile")


def load_counts(path: Path = DATA) -> dict[str, int]:
    """Read only printed integer facts; refuse missing provenance or bad panels."""
    counts: dict[str, int] = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != list(FIELDS):
            raise ValueError("CSV columns differ from the frozen provenance schema")
        for row in reader:
            if None in row or any(row.get(field) is None or not row[field].strip() for field in FIELDS):
                raise ValueError("Every aggregate needs complete provenance")
            identifier = row["id"]
            if identifier in counts:
                raise ValueError("Duplicate aggregate ID: " + identifier)
            value = row["count"]
            if not value.isascii() or not value.isdigit():
                raise ValueError("Counts must be printed nonnegative integers")
            counts[identifier] = int(value)
    _validate_counts(counts)
    return counts


def _validate_binomial(x: int, n: int) -> None:
    if type(x) is not int or type(n) is not int or n < 0 or not 0 <= x <= n:
        raise ValueError("Require integer binomial counts 0 <= x <= n")


def binomial_tail(x: int, n: int, p: float | Fraction) -> float | Fraction:
    """P(Binomial(n,p) >= x); Fraction input gives an exact rational result."""
    if type(n) is not int or n < 0 or type(x) is not int or not 0 <= x <= n + 1:
        raise ValueError("Require integer n >= 0 and 0 <= x <= n+1")
    if not 0 <= p <= 1:
        raise ValueError("Binomial probability must lie in [0,1]")
    exact = isinstance(p, Fraction)
    if x == 0:
        return Fraction(1) if exact else 1.0
    if x == n + 1 or p == 0:
        return Fraction(0) if exact else 0.0
    if p == 1:
        return Fraction(1) if exact else 1.0
    if exact:
        a, denominator = p.numerator, p.denominator
        numerator = sum(
            math.comb(n, j) * a**j * (denominator - a) ** (n - j)
            for j in range(x, n + 1)
        )
        return Fraction(numerator, denominator**n)
    return min(1.0, max(0.0, math.fsum(
        math.comb(n, j) * p**j * (1.0 - p) ** (n - j)
        for j in range(x, n + 1)
    )))


@lru_cache(maxsize=None)
def lower_bound(x: int, n: int, alpha: float = ALPHA) -> float:
    """One-sided exact binomial lower endpoint, inverted by bisection."""
    _validate_binomial(x, n)
    if not 0 < alpha < 1:
        raise ValueError("Tail probability must lie strictly between zero and one")
    if x == 0:
        return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(80):
        middle = (lo + hi) / 2.0
        if binomial_tail(x, n, middle) < alpha:
            lo = middle
        else:
            hi = middle
    return (lo + hi) / 2.0


@lru_cache(maxsize=None)
def cp_interval(x: int, n: int) -> tuple[float, float]:
    """Ordinary two-sided 95% Clopper-Pearson interval."""
    _validate_binomial(x, n)
    return lower_bound(x, n, ALPHA), 1.0 - lower_bound(n - x, n, ALPHA)


def discordant_pvalue(x: int, y: int) -> Fraction:
    """Exact one-sided paired-discordant test of P(+-) <= P(-+)."""
    if type(x) is not int or type(y) is not int or min(x, y) < 0:
        raise ValueError("Discordant counts must be nonnegative integers")
    d = x + y
    return Fraction(sum(math.comb(d, j) for j in range(x, d + 1)), 2**d)


def completions(m: int) -> list[tuple[int, int, int, int]]:
    """Every allocation to (++,+-,-+,--) in deterministic order."""
    if type(m) is not int or m < 0:
        raise ValueError("Missing count must be a nonnegative integer")
    return [
        (a, b, c, m - a - b - c)
        for a in range(m + 1)
        for b in range(m - a + 1)
        for c in range(m - a - b + 1)
    ]


def sharpness_allocation(
    q: Fraction, cap4: Fraction, cap7: Fraction, capb: Fraction
) -> tuple[Fraction, Fraction, Fraction] | None:
    """Allocate the observed reversal mass to disjoint errors/reversals.

    The returned masses receive latent 00,11,10 respectively. Other observed
    cells retain their own latent labels. None means these upper caps cannot fit.
    """
    q, cap4, cap7, capb = (Fraction(value) for value in (q, cap4, cap7, capb))
    if any(value < 0 or value > 1 for value in (q, cap4, cap7, capb)):
        raise ValueError("Mass and probability caps must lie in [0,1]")
    if cap4 + cap7 + capb < q:
        return None
    a4 = min(q, cap4)
    a7 = min(q - a4, cap7)
    return a4, a7, q - a4 - a7


def _discordant_results(x: int, y: int) -> dict:
    p = discordant_pvalue(x, y)
    d = x + y
    fraction = Fraction(x, d) if d else None
    return {
        "n": d,
        "reverse": x,
        "forward": y,
        "fraction": str(fraction) if fraction is not None else None,
        "value": float(fraction) if fraction is not None else None,
        "ci_95": list(cp_interval(x, d)),
        "p_fraction": str(p),
        "p_value": float(p),
        "reject_primary": p < Fraction(1, 40),
        "reject_marginal": p < Fraction(1, 20),
        "verdict": FALSIFIED if p < Fraction(1, 40) else COMPATIBLE,
    }


def compute_results(counts: dict[str, int] | None = None,
                    paired_classification_supported: bool = True) -> dict:
    """Recompute all results; eligibility is a manually audited source gate."""
    counts = load_counts() if counts is None else dict(counts)
    _validate_counts(counts)
    if type(paired_classification_supported) is not bool:
        raise ValueError("The paired-classification source gate must be explicit")
    result = {
        "counts": counts,
        "hashes": {
            path.relative_to(REPOSITORY).as_posix(): normalized_sha256(path)
            for path in (PROTOCOL, DATA)
        },
        "design": {
            "family_alpha": FAMILY_ALPHA,
            "per_test_alpha": ALPHA,
            "primary_cap_fraction": str(PRIMARY_CAP),
            "source_access_date": "2026-11-18",
        },
        "paired_classification_supported": paired_classification_supported,
        "successful": None,
        "missing": None,
        "exclusion_audit": {
            "printed_total": counts["itt_excluded"],
            "parts_total": sum(counts[key] for key in (
                "itt_missing_week7", "itt_ineligible", "itt_withdrew")),
            "matches": True,
        },
        "verdicts": {
            "paired_reporting_gate": COMPATIBLE if paired_classification_supported else UNTESTABLE,
            "individual_availability_audit": UNTESTABLE,
            "unconditional_biological_no_intervention_claim": UNTESTABLE,
            "independently_calibrated_error_or_reversal_rates": UNTESTABLE,
            "assay_error_shedding_or_biological_clearance_attribution": UNTESTABLE,
            "recurrence_sensitivity": UNTESTABLE,
            "schedule_optimum": UNTESTABLE,
            "treatment_decision": UNTESTABLE,
        },
    }
    x, y = counts["week4_only"], counts["week7_only"]
    if not paired_classification_supported:
        result["source_supported_discordants"] = _discordant_results(x, y)
        result["source_supported_discordants"]["verdict"] = UNTESTABLE
        result["source_supported_discordants"]["qualification"] = (
            "Numerical conditional calculation only; ascertainment of all discordants "
            "and a full-pair denominator are not supported."
        )
        return result

    n, total = counts["successful_analysis"], counts["guided_itt"]
    m = total - n
    q = Fraction(x, n)
    lower = lower_bound(x, n, ALPHA)
    table = dict(zip(CELLS, (
        counts["both_positive"], x, y, counts["reported_negative"])))
    n4, n7 = table["++"] + x, table["++"] + y
    cap_rows = []
    for cap in CAPS:
        tail = binomial_tail(x, n, cap)
        reject = tail < Fraction(1, 40)
        cap_rows.append({
            "cap_fraction": str(cap), "cap": float(cap),
            "tail_fraction": str(tail), "tail": float(tail),
            "exact_table_compatible": cap >= q,
            "reject": reject, "verdict": FALSIFIED if reject else COMPATIBLE,
        })
    tradeoff = []
    for b in REVERSAL_CAPS:
        error = max(Fraction(0), q - b)
        lower_error = max(0.0, lower - float(b))
        tradeoff.append({
            "b_fraction": str(b), "b": float(b),
            "combined_error_fraction": str(error), "combined_error": float(error),
            "equal_error_fraction": str(error / 2), "equal_error": float(error / 2),
            "lower_combined_error_975": lower_error,
            "lower_equal_error_975": lower_error / 2,
        })
    successful = {
        "n": n, "table": table,
        "q": {
            "fraction": str(q), "value": float(q),
            "lower_95": lower_bound(x, n, 0.05),
            "lower_975": lower, "ci_95": list(cp_interval(x, n)),
        },
        "marginals": {
            "week4_positive_count": n4, "week7_positive_count": n7,
            "week4_positive_fraction": str(Fraction(n4, n)),
            "week7_positive_fraction": str(Fraction(n7, n)),
            "week4_positive_rate": n4 / n, "week7_positive_rate": n7 / n,
            "difference_fraction": str(Fraction(n7 - n4, n)),
            "difference": (n7 - n4) / n,
        },
        "discordant": _discordant_results(x, y), "caps": cap_rows,
        "tradeoff": tradeoff,
        "primary_verdict": next(row["verdict"] for row in cap_rows
                                if row["cap_fraction"] == str(PRIMARY_CAP)),
    }
    result["successful"] = successful

    completed = []
    for allocation in completions(m):
        nx, ny = x + allocation[1], y + allocation[2]
        probability = Fraction(nx, total)
        p = discordant_pvalue(nx, ny)
        completed.append({
            "allocation": dict(zip(CELLS, allocation)),
            "q_fraction": str(probability), "q": float(probability),
            "lower_975": lower_bound(nx, total, ALPHA),
            "ci_95": list(cp_interval(nx, total)),
            "discordant_p_fraction": str(p), "discordant_p": float(p),
            "discordant_reject": p < Fraction(1, 40),
        })
    missing_caps = []
    for cap in CAPS:
        tails = [binomial_tail(x + row["allocation"]["+-"], total, cap)
                 for row in completed]
        rejects = sum(tail < Fraction(1, 40) for tail in tails)
        scope = "every" if rejects == len(completed) else "some" if rejects else "none"
        missing_caps.append({
            "cap_fraction": str(cap), "cap": float(cap),
            "rejection_scope": scope, "rejecting_completions": rejects,
            "total_completions": len(completed),
            "min_tail_fraction": str(min(tails)), "max_tail_fraction": str(max(tails)),
            "min_tail": float(min(tails)), "max_tail": float(max(tails)),
        })
    q_bounds = (Fraction(x, total), Fraction(x + m, total))
    missing_tradeoff = []
    for b in REVERSAL_CAPS:
        errors = tuple(max(Fraction(0), value - b) for value in q_bounds)
        missing_tradeoff.append({
            "b_fraction": str(b), "b": float(b),
            "combined_error_bounds_fractions": [str(value) for value in errors],
            "combined_error_bounds": [float(value) for value in errors],
            "equal_error_bounds_fractions": [str(value / 2) for value in errors],
            "equal_error_bounds": [float(value / 2) for value in errors],
        })
    p_values = [Fraction(row["discordant_p_fraction"]) for row in completed]
    min_p, max_p = min(p_values), max(p_values)
    result["missing"] = {
        "n_total": total, "m": m,
        "q_bounds": {"fractions": [str(value) for value in q_bounds],
                     "values": [float(value) for value in q_bounds]},
        "lower_975_bounds": [lower_bound(x, total, ALPHA), lower_bound(x + m, total, ALPHA)],
        "ci_95_outer": [lower_bound(x, total, ALPHA), cp_interval(x + m, total)[1]],
        "cap_sensitivity": missing_caps, "tradeoff": missing_tradeoff,
        "completions": completed,
        "directional_extrema": {
            "min_fraction": str(min_p), "max_fraction": str(max_p),
            "min_value": float(min_p), "max_value": float(max_p),
            "min_attainers": [row["allocation"] for row in completed
                              if Fraction(row["discordant_p_fraction"]) == min_p],
            "max_attainers": [row["allocation"] for row in completed
                              if Fraction(row["discordant_p_fraction"]) == max_p],
            "rejecting_completions": sum(row["discordant_reject"] for row in completed),
        },
    }
    result["verdicts"]["primary_cap_successful_population"] = successful["primary_verdict"]
    result["verdicts"]["conditional_noisy_directional_restriction"] = successful["discordant"]["verdict"]
    result["verdicts"]["unrestricted_total_latitude"] = COMPATIBLE
    return result


def _pct(value: float) -> str:
    return f"{100 * value:.6f}%"


def _interval(values: list[float], percentage: bool = True) -> str:
    fmt = _pct if percentage else lambda value: f"{value:.9f}"
    return "[" + ", ".join(fmt(value) for value in values) + "]"


def _allocations(values: list[dict]) -> str:
    return "; ".join("(" + ", ".join(str(row[cell]) for cell in CELLS) + ")" for row in values)


def format_results(result: dict) -> str:
    """Complete generated README block; no writes and no source prose copied."""
    lines = [
        "### Recomputed results",
        "",
        "**Nonclinical, unreviewed. No treatment, diagnosis or surveillance advice.**",
        "",
        "All numerical outputs below are calculated from the attributed aggregate CSV. "
        "Intervals assume independent identically distributed case pairs in the stated "
        "population; calls within a case may be dependent. Rates are unconditional "
        "within that population. The two primary tests each use alpha = 0.025, for "
        "Bonferroni familywise error at most 0.05 within the stated target population. "
        "The full-arm calculations are missingness robustness sensitivities; no "
        "additional simultaneous coverage across both populations or all displayed "
        "intervals is claimed. Equality at a rejection boundary is retained.",
        "",
    ]
    if not result["paired_classification_supported"]:
        d = result["source_supported_discordants"]
        lines.extend([
            "**Paired-availability gate: untestable.** A full paired denominator and "
            "the inferred negative/negative cell are not supported. Full-table and "
            "missing-analysis inference are withheld.", "",
            f"Directly reported discordants: {d['reverse']} +/− and {d['forward']} −/+. "
            f"The purely conditional arithmetic gives p = {d['p_fraction']} "
            f"({d['p_value']:.9f}); unverified directional ascertainment prevents an unconditional test.",
        ])
    else:
        s, missing = result["successful"], result["missing"]
        q, d, marginal = s["q"], s["discordant"], s["marginals"]
        lines.extend([
            "#### Reported paired classification and sharp frontier", "",
            "The source cross-check supports complete **reported classification**. "
            "The negative/negative cell is inferred from the printed negative category "
            "and its two-draw definition; an individual assay-availability audit and a "
            "separately printed complete-pairs count remain unavailable.", "",
            "| Week-4 / week-7 category | Count |", "| --- | ---: |",
            f"| +/+ | {s['table']['++']} |", f"| +/− | {s['table']['+-']} |",
            f"| −/+ | {s['table']['-+']} |", f"| −/−, inferred reported category | {s['table']['--']} |",
            f"| Successful-analysis denominator | {s['n']} |", "",
            "| Quantity | Recomputed value |", "| --- | --- |",
            f"| Exact plug-in total-latitude frontier q | {q['fraction']} = {_pct(q['value'])} |",
            f"| Decision lower frontier, one-sided 97.5% | {_pct(q['lower_975'])} |",
            f"| Marginal descriptive lower frontier, one-sided 95% | {_pct(q['lower_95'])} |",
            f"| Ordinary two-sided 95% exact interval for q | {_interval(q['ci_95'])} |",
            f"| Week-4 positive calls / marginal rate | {marginal['week4_positive_count']}; {marginal['week4_positive_fraction']} = {_pct(marginal['week4_positive_rate'])} |",
            f"| Week-7 positive calls / marginal rate | {marginal['week7_positive_count']}; {marginal['week7_positive_fraction']} = {_pct(marginal['week7_positive_rate'])} |",
            f"| Week-7 minus week-4 marginal probability | {marginal['difference_fraction']} = {_pct(marginal['difference'])} |",
            "",
            "For upper caps E4, E7 and B, exact table compatibility holds precisely when "
            "E4 + E7 + B >= q. With B = 0, the combined observation-error allowance "
            "must reach the plug-in frontier for an exact empirical-table fit. The "
            "uncertainty rule rejects only below the decision lower frontier. The "
            "sum e4 + e7 is an expected number of occasion mismatches per participant, "
            "not necessarily the probability of any mismatch; neither occasion has "
            "an identified individual lower bound.", "",
            "| Hypothetical total cap | Exact binomial upper tail | Exact plug-in table fits? | Population uncertainty verdict |",
            "| --- | ---: | --- | --- |",
        ])
        for row in s["caps"]:
            lines.append(f"| {_pct(row['cap'])} | {row['tail']:.12g} | {'yes' if row['exact_table_compatible'] else 'no'} | **{row['verdict']}** |")
        lines.extend([
            "", f"**Primary total cap 1%: {s['primary_verdict']}.** "
            "All caps are hypothetical and uncalibrated. The nested sensitivity "
            "decisions use the same lower bound. Exact rational tails are retained "
            "in the JSON output; the displayed decimals are rounded.", "",
            "| Hypothetical reversal cap B | Minimum combined error for exact fit | Minimum equal cap per occasion for exact fit | Decision lower combined-error frontier | Decision lower equal-cap frontier |",
            "| --- | ---: | ---: | ---: | ---: |",
        ])
        for row in s["tradeoff"]:
            lines.append(f"| {_pct(row['b'])} | {_pct(row['combined_error'])} | {_pct(row['equal_error'])} | {_pct(row['lower_combined_error_975'])} | {_pct(row['lower_equal_error_975'])} |")
        lines.extend([
            "", "These are max(0, q − B) and half that value, with the same transformations "
            "of the decision lower bound. Allowing genuine reversal changes the "
            "persistence assumption; unrestricted reversal need not require any "
            "observation error. Constructed compatibility does not identify the mechanism.", "",
            "#### Unsuccessful analyses: all completions, no missing-at-random assumption", "",
            f"The full guided ITT denominator is {missing['n_total']}; subtracting "
            f"{s['n']} successful analyses leaves {missing['m']} unsuccessful analyses. "
            f"There are {len(missing['completions'])} allocations to (++,+−,−+,−−). "
            "These hypothetical outcomes are sensitivity completions, not recovered calls. "
            "The full-arm and successful-analysis error probabilities are different estimands.", "",
            "| Full-arm quantity | Recomputed value |", "| --- | --- |",
            f"| Sharp finite-panel q range | {', '.join(missing['q_bounds']['fractions'])}; {_interval(missing['q_bounds']['values'])} |",
            f"| Completion-specific one-sided 97.5% lower-frontier range | {_interval(missing['lower_975_bounds'])} |",
            f"| Completion-robust outer two-sided 95% interval | {_interval(missing['ci_95_outer'])} |", "",
            "The lower finite-panel endpoint is enough for at least one completion; "
            "the upper endpoint allows every completion. They are sharp identification "
            "bounds, not confidence intervals. The lower uncertainty endpoint is valid "
            "for every completion under the full-arm iid sampling model even when "
            "missingness is arbitrary.", "",
            "| Total cap | Completions rejected | Exact-tail range, rounded | Rejection scope |",
            "| --- | ---: | --- | --- |",
        ])
        for row in missing["cap_sensitivity"]:
            scope_label = {"every": "every completion", "some": "some completions", "none": "no completions"}[row["rejection_scope"]]
            lines.append(f"| {_pct(row['cap'])} | {row['rejecting_completions']}/{row['total_completions']} | [{row['min_tail']:.12g}, {row['max_tail']:.12g}] | {scope_label} |")
        lines.extend([
            "", "Every-completion rejection is **falsified under stated caps** despite "
            "the missing analyses. Some-completion rejection means the conclusion "
            "depends on completion; none means **compatible with stated latitude** "
            "by the uncertainty rule for every completion. This is not a missing-at-random estimate.", "",
            "| Reversal cap B | Combined error needed: some to every completion | Equal per-occasion cap needed: some to every completion |",
            "| --- | --- | --- |",
        ])
        for row in missing["tradeoff"]:
            lines.append(f"| {_pct(row['b'])} | {_interval(row['combined_error_bounds'])} | {_interval(row['equal_error_bounds'])} |")
        extrema = missing["directional_extrema"]
        lines.extend([
            "", "#### Separate comparable-kernel paired-discordant test", "",
            "This restriction requires nondecreasing **marginal detection probability**, "
            "a comparable nondecreasing assay kernel, nondecreasing detectable burden "
            "and no intervention-induced change. Persistence plus unrestricted "
            "occasion-specific errors alone does not imply it.", "",
            "| Quantity | Recomputed value |", "| --- | --- |",
            f"| Discordant cases, +/− versus −/+ | {d['n']}; {d['reverse']} versus {d['forward']} |",
            f"| Conditional reversal fraction | {d['fraction']} = {_pct(d['value']) if d['value'] is not None else 'undefined'} |",
            f"| Two-sided 95% exact interval for conditional reversal fraction | {_interval(d['ci_95'])} |",
            f"| Exact one-sided p-value for H0: P(+/−) <= P(−/+) | {d['p_fraction']} = {d['p_value']:.9f} |",
            f"| Reject at 0.025 / descriptive 0.05? | {'yes' if d['reject_primary'] else 'no'} / {'yes' if d['reject_marginal'] else 'no'} |",
            f"| Conditional statistical verdict | **{d['verdict']}** |",
            f"| Minimum p over missing-case completions | {extrema['min_fraction']} = {extrema['min_value']:.9f} |",
            f"| Minimum-p attaining allocation(s), (++,+−,−+,−−) | {_allocations(extrema['min_attainers'])} |",
            f"| Maximum p over missing-case completions | {extrema['max_fraction']} = {extrema['max_value']:.9f} |",
            f"| Maximum-p attaining allocation(s), (++,+−,−+,−−) | {_allocations(extrema['max_attainers'])} |",
            f"| Completions rejected at 0.025 | {extrema['rejecting_completions']}/{len(missing['completions'])} |", "",
            "The noisy restriction's nonrejection is retained. PR #58 already permits "
            "finite-molecule misses, so positive/negative pairs alone cannot falsify "
            "its noisy model. This tests a necessary marginal restriction, not the "
            "full growth model. **Unconditional biological/no-intervention verdict: "
            "untestable. Individual intervention dates are not verified.**", "",
            "#### Source reconciliation and limits", "",
            f"The Figure 1 exclusion categories sum to {result['exclusion_audit']['parts_total']}, "
            f"matching the printed {result['exclusion_audit']['printed_total']} ITT exclusions. "
            "Those exclusions are outside the analyzed ITT population and do not "
            "replace the unsuccessful-analysis sensitivity panel.", "",
            "Independent error calibration, an individual availability audit, "
            "assay-error/shedding/biological-clearance attribution, recurrence "
            "sensitivity, a schedule optimum and treatment implications are all "
            "**untestable** from these aggregates. The analysis cannot distinguish "
            "assay error, shedding or biological clearance; estimate recurrence "
            "sensitivity; validate a schedule optimum; or imply a treatment decision.",
        ])
    lines.extend([
        "", "#### Reproducibility hashes", "",
        "SHA-256 of UTF-8 text after CRLF and CR become LF. Path keys use "
        "`Path.as_posix()`. These hashes check integrity, not source truth or clinical validity.", "",
        "| Input | Normalized SHA-256 |", "| --- | --- |",
    ])
    for path, digest in result["hashes"].items():
        lines.append(f"| `{path}` | `{digest}` |")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--json", action="store_true", help="Print exact-rational and floating results")
    group.add_argument("--check", action="store_true", help="Check the README generated block without writing")
    args = parser.parse_args(argv)
    try:
        result = compute_results()
        if args.json:
            print(json.dumps(result, indent=2, sort_keys=True))
            return 0
        rendered = format_results(result)
        if args.check:
            text = README.read_text(encoding="utf-8")
            if text.count(BEGIN) != 1 or text.count(END) != 1:
                raise ValueError("README must contain exactly one pair of result markers")
            block = text.split(BEGIN, 1)[1].split(END, 1)[0].strip()
            if block != rendered:
                raise ValueError("README result block is stale")
            print("Nonclinical paired-call result block matches the aggregate recomputation.")
        else:
            print(rendered)
        return 0
    except (OSError, ValueError) as error:
        print("Paired-call audit failed: " + str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
