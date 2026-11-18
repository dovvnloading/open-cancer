#!/usr/bin/env python3
"""Frozen, offline incidence-stratum analysis; Python standard library only.

Nonclinical and unreviewed computational research. No medical advice.
Rounding compatibility is not a sampling test, and k is not a driver estimate.
The numerical primitives are inherited unchanged from PR #64. Every empirical
README table is generated here from the selected public NCI aggregate CSVs.
"""

from __future__ import annotations

import argparse
import csv
from datetime import date
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import importlib.util
from itertools import combinations
import json
from math import factorial
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
METHOD = ROOT.parent / "multistage-incidence-identifiability" / "verifier.py"
_SPEC = importlib.util.spec_from_file_location("incidence_strata_pr64_verifier", METHOD)
V = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = V
_SPEC.loader.exec_module(V)

PROTOCOL_COMMIT = "f06822e9f2427cbec8fef16c48589290561b40c3"
LABEL = "Nonclinical, unreviewed computational research. No medical advice."
KS = tuple(range(2, 11))
GROWTH_GRID = tuple(F(i, 200) for i in range(-12, 13))
SERIES_TERMS = 70
AGE_STARTS = tuple(range(30, 85, 5))
AGES = tuple(f"{a}-{a + 4}" for a in AGE_STARTS)
RACES = ("All Races", "Whites", "Blacks")
SEXES = ("Both Sexes", "Males", "Females")
SITES = {
    "female_breast": ("4.11", "sect_04_table.11", ("Females",)),
    "colon_and_rectum": ("6.10", "sect_06_table.10", SEXES),
    "lung_and_bronchus": ("15.10", "sect_15_table.10", SEXES),
    "pancreas": ("22.7", "sect_22_table.07", SEXES),
}
DATA_FILES = tuple(Path("data") / f"{site}.csv" for site in SITES)
SOURCE_BASE = "https://seer.cancer.gov/archive/csr/1975_2017/"
BEGIN = "<!-- BEGIN RECOMPUTED RESULTS -->"
END = "<!-- END RECOMPUTED RESULTS -->"


def path_key(path: Path, root: Path) -> str:
    """Portable relative keys, including when tested with PureWindowsPath."""
    return path.relative_to(root).as_posix()


def source_sha256(path: Path) -> str:
    """Hash decoded UTF-8 with universal-newline normalization, and nothing else."""
    with path.open(encoding="utf-8", newline=None) as stream:
        return sha256(stream.read().encode("utf-8")).hexdigest()


def display(value: F, places: int = 6, *, sqrt: bool = False) -> str:
    with localcontext() as context:
        context.prec = 70
        number = Decimal(value.numerator) / value.denominator
        if sqrt:
            number = number.sqrt()
        return f"{number:.{places}f}"


def group_id(site: str, race: str, sex: str) -> str:
    return "__".join((site, race.lower().replace(" ", "_"),
                     sex.lower().replace(" ", "_")))


def expected_groups() -> dict[str, tuple[str, str, str]]:
    return {group_id(site, race, sex): (site, race, sex)
            for site, (_, _, sexes) in SITES.items()
            for race in RACES for sex in sexes}


def load_data(data_dir: Path = ROOT / "data") -> dict[str, list[dict[str, str]]]:
    """Validate provenance/selection; retain absent rates for explicit abstention."""
    data_dir = Path(data_dir)
    expected = expected_groups()
    groups: dict[str, list[dict[str, str]]] = {key: [] for key in expected}
    seen = set()
    for site, (table, stem, sexes) in SITES.items():
        path = data_dir / f"{site}.csv"
        if not path.exists():
            continue
        with path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed CSV row")
            gid = row.get("group_id")
            if gid not in expected or expected[gid] != (
                    row.get("cancer_site"), row.get("race"), row.get("sex")):
                raise ValueError("Unexpected site/race/sex group")
            if row["cancer_site"] != site or row["sex"] not in sexes:
                raise ValueError("Group is in the wrong site file")
            fixed = {
                "behavior": "invasive", "registry_set": "SEER 21 areas",
                "period_start": "2013", "period_end": "2017",
                "reported_decimal_places": "1", "source_table": table,
                "source_section": "Age-Specific SEER Incidence Rates, 2013-2017",
                "source_column": f"{row['race']} / {row['sex']}",
                "source_page": "1", "rate_units": "cases per 100000 population per year",
                "source_name": "National Cancer Institute, SEER Cancer Statistics Review 1975-2017",
                "source_version": "November 2019 SEER submission; released 2020-04-15",
                "source_url": f"{SOURCE_BASE}download_csr_datafile.php/{stem}.csv",
                "source_pdf_url": f"{SOURCE_BASE}results_single/{stem}.pdf",
                "source_html_url": (f"{SOURCE_BASE}browse_csr.php?pageSEL={stem}"
                                    f"&sectionSEL={int(stem[5:7])}"),
                "reuse_terms_url": SOURCE_BASE,
            }
            for key, value in fixed.items():
                if row.get(key) != value:
                    raise ValueError(f"Inconsistent metadata: {key}")
            for key in ("source_name", "source_version", "source_html_url",
                        "access_date", "reuse_basis"):
                if not row.get(key, "").strip():
                    raise ValueError(f"Missing provenance: {key}")
            if "public domain" not in row["reuse_basis"]:
                raise ValueError("Missing exact-source public-domain reuse basis")
            date.fromisoformat(row["access_date"])
            if row.get("case_count") != "" or row.get("person_years") != "":
                raise ValueError("Selected tables do not supply counts or person-years")
            start, end = int(row["age_start"]), int(row["age_end_exclusive"])
            if (start not in AGE_STARTS or end != start + 5
                    or row["age_group"] != f"{start}-{end - 1}"):
                raise ValueError("Age bin is outside the frozen selection")
            key = gid, start
            if key in seen:
                raise ValueError(f"Duplicate group/age: {key}")
            seen.add(key)
            groups[gid].append(row)
    for rows in groups.values():
        rows.sort(key=lambda row: int(row["age_start"]))
        for key in ("source_name", "source_version", "source_html_url", "access_date"):
            if len({row[key] for row in rows}) > 1:
                raise ValueError(f"Inconsistent group provenance: {key}")
    return groups


def eligibility_reason(rows: list[dict[str, str]]) -> str | None:
    if tuple(row["age_group"] for row in rows) != AGES:
        return "missing required age bins"
    for row in rows:
        raw = row.get("rate_per_100000", "")
        if raw in ("", "-"):
            return f"missing or suppressed rate at {row['age_group']}"
        try:
            value = F(raw)
        except (ValueError, ZeroDivisionError):
            return f"nonnumeric rate at {row['age_group']}"
        if value - F(1, 20) <= 0:
            return f"nonpositive lower rounding endpoint at {row['age_group']}"
        if re.fullmatch(r"\d+\.\d", raw) is None:
            return f"unexpected printed precision at {row['age_group']}"
    return None


def bounds(rows: list[dict[str, str]], half: F = F(1, 20)) -> tuple[tuple[F, ...], tuple[F, ...]]:
    values = tuple(F(row["rate_per_100000"]) for row in rows)
    return tuple(y - half for y in values), tuple(y + half for y in values)


def monotonicity(rows: list[dict[str, str]], half: F = F(1, 20)) -> dict:
    if not rows:
        raise ValueError("Require nonempty ordered rows")
    lower, upper = bounds(rows, half)
    if half < 0 or any(l <= 0 for l in lower):
        raise ValueError("Require positive intervals and nonnegative half width")
    candidates = [(lower[i] / upper[j], i, j)
                  for i in range(len(rows)) for j in range(i, len(rows))]
    ratio = max(F(1), max(item[0] for item in candidates))
    violations = [{
        "earlier_age": rows[i]["age_group"], "later_age": rows[j]["age_group"],
        "earlier_rate": rows[i]["rate_per_100000"],
        "later_rate": rows[j]["rate_per_100000"],
        "minimum_decline": lower[i] - upper[j], "rho_squared": value,
    } for value, i, j in candidates if i < j and value > 1]
    limiting = [{"earlier_age": rows[i]["age_group"],
                 "later_age": rows[j]["age_group"]}
                for value, i, j in candidates if ratio > 1 and value == ratio]
    return {"rho_squared": ratio, "rho": display(ratio, sqrt=True),
            "compatible": ratio == 1, "violations": violations,
            "limiting_pairs": limiting}


def age_sensitivity(rows: list[dict[str, str]]) -> dict:
    prefixes = []
    for i, row in enumerate(rows):
        result = monotonicity(rows[:i + 1])
        prefixes.append({"ending_age": row["age_group"],
                         "rho_squared": result["rho_squared"], "rho": result["rho"]})
    leave_one_out = []
    for i, row in enumerate(rows):
        result = monotonicity(rows[:i] + rows[i + 1:])
        leave_one_out.append({"removed_age": row["age_group"],
                              "rho_squared": result["rho_squared"], "rho": result["rho"]})
    # Restore feasibility iff every incompatible pair loses at least one end.
    # Enumerating vertex covers here differs from the independent test's
    # enumeration of interval-feasible retained subsequences.
    lower, upper = bounds(rows)
    edges = [(i, j) for i in range(len(rows)) for j in range(i + 1, len(rows))
             if lower[i] > upper[j]]
    deletion_sets = []
    for size in range(len(rows)):
        for removed in combinations(range(len(rows)), size):
            indices = set(removed)
            if all(i in indices or j in indices for i, j in edges):
                deletion_sets.append([rows[i]["age_group"] for i in removed])
        if deletion_sets:
            break
    return {"prefixes": prefixes, "leave_one_out": leave_one_out,
            "first_failure_age": next((p["ending_age"] for p in prefixes
                                       if p["rho_squared"] > 1), None),
            "minimum_deletion_sets": deletion_sets}


def exact_frontier(lower: tuple[F, ...], upper: tuple[F, ...],
                   shapes: tuple[F, ...], ages: tuple[str, ...]) -> dict:
    left = [l / f for l, f in zip(lower, shapes)]
    right = [u / f for u, f in zip(upper, shapes)]
    ratio = V.calibration_ratio(lower, upper, shapes)
    return {"rho_squared": ratio, "rho": display(ratio, sqrt=True),
            "lower_ages": [a for a, value in zip(ages, left) if value == max(left)],
            "upper_ages": [a for a, value in zip(ages, right) if value == min(right)]}


def flexible_frontier(lower: tuple[F, ...], upper: tuple[F, ...], k: int) -> dict:
    ranges = tuple((F(a ** (k - 1), factorial(k - 1)),
                    F((a + 5) ** (k - 1), factorial(k - 1))) for a in AGE_STARTS)
    left = [l / hi for l, (_, hi) in zip(lower, ranges)]
    right = [u / lo for u, (lo, _) in zip(upper, ranges)]
    ratio = V.flexible_shape_calibration(lower, upper, ranges)
    return {"rho_squared": ratio, "rho": display(ratio, sqrt=True),
            "lower_ages": [a for a, value in zip(AGES, left) if value == max(left)],
            "upper_ages": [a for a, value in zip(AGES, right) if value == min(right)]}


@lru_cache(maxsize=1)
def shape_grid() -> dict:
    return {(k, g): tuple(V.one_growth_bin_interval(k, g, F(a), F(a + 5), SERIES_TERMS)
                         for a in AGE_STARTS) for k in KS for g in GROWTH_GRID}


def interval_frontier(lower: tuple[F, ...], upper: tuple[F, ...],
                      shapes: tuple[tuple[F, F], ...]) -> dict:
    lo, hi = V.calibration_ratio_interval(lower, upper, shapes)
    left = [(l / b, l / a) for l, (a, b) in zip(lower, shapes)]
    right = [(u / b, u / a) for u, (a, b) in zip(upper, shapes)]
    # Conservative sets include every possible extremizer. A singleton is
    # certified because its lower/upper bound dominates every competitor.
    left_floor, right_ceiling = max(a for a, _ in left), min(b for _, b in right)
    return {"rho_squared_interval": (lo, hi),
            "lower_ages": [age for age, (_, b) in zip(AGES, left) if b >= left_floor],
            "upper_ages": [age for age, (a, _) in zip(AGES, right) if a <= right_ceiling]}


def bounded_rho(interval: tuple[F, F]) -> str:
    a, b = (display(value, sqrt=True) for value in interval)
    if a != b:
        raise ArithmeticError("Numerical enclosure cannot certify six-decimal display")
    return a


def model_frontiers(rows: list[dict[str, str]]) -> dict:
    lower, upper = bounds(rows)
    shapes = shape_grid()
    by_k = []
    all_candidates = []
    for k in KS:
        powers = tuple(V.power_bin(k, F(a), F(a + 5)) for a in AGE_STARTS)
        candidates = []
        for g in GROWTH_GRID:
            result = interval_frontier(lower, upper, shapes[k, g])
            candidates.append(dict(result, g=g, k=k))
        best = min(candidates, key=lambda c: (c["rho_squared_interval"][1],
                                             c["rho_squared_interval"][0], c["g"]))
        hi = best["rho_squared_interval"][1]
        lo = min(c["rho_squared_interval"][0] for c in candidates)
        contenders = [c["g"] for c in candidates if c["rho_squared_interval"][0] <= hi]
        if hi - lo >= F(1, 10**24):
            raise ArithmeticError("Unresolved squared-threshold precision")
        by_k.append({
            "k": k, "power": exact_frontier(lower, upper, powers, AGES),
            "flexible_power": flexible_frontier(lower, upper, k),
            "grid": {"rho_squared_interval": (lo, hi), "rho": bounded_rho((lo, hi)),
                     "selected_g": best["g"], "contenders": contenders,
                     "unique": len(contenders) == 1,
                     "lower_ages": best["lower_ages"], "upper_ages": best["upper_ages"]},
        })
        all_candidates.extend(candidates)
    power_min = min(item["power"]["rho_squared"] for item in by_k)
    best_power_contenders = [item["k"] for item in by_k if item["power"]["rho_squared"] == power_min]
    selected = min(all_candidates, key=lambda c: (c["rho_squared_interval"][1],
                                                  c["rho_squared_interval"][0], c["k"], c["g"]))
    global_hi = selected["rho_squared_interval"][1]
    contenders = [{"k": c["k"], "g": c["g"]} for c in all_candidates
                  if c["rho_squared_interval"][0] <= global_hi]
    return {"by_k": by_k, "best_power_k": best_power_contenders[0],
            "best_power_contenders": best_power_contenders,
            "best_grid_k": selected["k"], "best_grid_contenders": contenders}


def primary_verdict(groups: dict) -> str:
    breast = [groups[group_id("female_breast", race, "Females")] for race in RACES]
    if any(group["status"] != "testable" for group in breast):
        return "untestable"
    if breast[0]["monotonicity"]["compatible"]:
        return "untestable"
    if any(not group["monotonicity"]["compatible"] for group in breast[1:]):
        return "shape failure persists within strata"
    return "shape failure disappears within the prespecified strata"


def analyze(data_dir: Path = ROOT / "data") -> dict:
    data_dir = Path(data_dir)
    data = load_data(data_dir)
    groups = {}
    for gid, (site, race, sex) in expected_groups().items():
        rows = data[gid]
        reason = eligibility_reason(rows)
        result = {"site": site, "race": race, "sex": sex, "n": len(rows),
                  "status": "untestable" if reason else "testable", "reason": reason}
        if reason is None:
            result.update({"monotonicity": monotonicity(rows),
                           "point_monotonicity": monotonicity(rows, F(0)),
                           "age_sensitivity": age_sensitivity(rows)})
            result.update(model_frontiers(rows))
        groups[gid] = result
    hashes = {path_key(ROOT / "PROTOCOL.md", ROOT): source_sha256(ROOT / "PROTOCOL.md")}
    hashes.update({path_key(data_dir / relative.name, data_dir.parent):
                   source_sha256(data_dir / relative.name)
                   for relative in DATA_FILES if (data_dir / relative.name).exists()})
    return {"label": LABEL, "protocol_commit": PROTOCOL_COMMIT, "hashes": hashes,
            "inherited_method": {"path": path_key(METHOD, REPO), "sha256": source_sha256(METHOD)},
            "groups": groups, "verdict": primary_verdict(groups)}


def _name(group: dict) -> str:
    return f"{group['site'].replace('_', ' ')} / {group['race']} / {group['sex']}"


def _ages(ages: list[str]) -> str:
    return ", ".join(ages) or "none"


def _constraints(frontier: dict) -> str:
    return f"{_ages(frontier['lower_ages'])} / {_ages(frontier['upper_ages'])}"


def render_results(results: dict) -> str:
    """All numerical findings are in this mechanically checked README block."""
    groups = results["groups"]
    eligible = [group for group in groups.values() if group["status"] == "testable"]
    failing = [group for group in eligible if not group["monotonicity"]["compatible"]]
    breast = [groups[group_id("female_breast", race, "Females")] for race in RACES]
    breast_failed = [group["race"] for group in breast if group["status"] == "testable"
                     and not group["monotonicity"]["compatible"]]
    lines = [BEGIN, "", "## Recomputed results", "",
             f"**Verdict: {results['verdict']}.**", "",
             f"The CSVs contain **{sum(g['n'] for g in groups.values())} published rates** "
             f"in **{len(groups)} prespecified curves**; **{len(eligible)}** are testable. "
             f"**{len(failing)}** fail the rounded-interval monotone-mean condition.", "",
             "Breast curves failing the condition: **" + ", ".join(breast_failed) + "**.", ""]
    if all(group["status"] == "testable" for group in breast):
        specifics = breast[1:]
        if all(not group["monotonicity"]["compatible"] for group in specifics):
            lines.extend(["The failure persists within **both** printed race-specific breast strata. "
                          "Pooling across those race categories is therefore insufficient, by itself, "
                          "to explain away this descriptive failure.", ""])
    lines.extend([
        "### Primary and contrast monotonicity frontiers", "",
        "`rho` is the minimum for any nondecreasing positive mean vector. The exact "
        "rational `rho^2` is included so a rounded display is never used as a cutoff. "
        "Ages in a witness are earlier -> later. No witness is assigned to a passing curve.", "",
        "| Curve | Rounded rho | Exact rho^2 | Printed-point rho | Limiting age pair(s) | First failing prefix ends |",
        "|---|---:|---|---:|---|---|",
    ])
    for group in groups.values():
        if group["status"] != "testable":
            lines.append(f"| {_name(group)} | untestable | {group['reason']} | - | - | - |")
            continue
        mono = group["monotonicity"]
        witnesses = "; ".join(f"{p['earlier_age']} -> {p['later_age']}" for p in mono["limiting_pairs"]) or "none"
        lines.append(f"| {_name(group)} | {mono['rho']} | {mono['rho_squared']} | "
                     f"{group['point_monotonicity']['rho']} | {witnesses} | "
                     f"{group['age_sensitivity']['first_failure_age'] or 'none'} |")
    lines.extend(["", "### Every rounding-robust decrease", "",
                  "Rates and minimum declines are per 100,000. Every violating earlier/later pair is shown, "
                  "including nonadjacent ages; absence from this table means no violating pair.", "",
                  "| Curve | Earlier age | Printed rate | Later age | Printed rate | Minimum decline |",
                  "|---|---|---:|---|---:|---:|"])
    for group in failing:
        for pair in group["monotonicity"]["violations"]:
            lines.append(f"| {_name(group)} | {pair['earlier_age']} | {pair['earlier_rate']} | "
                         f"{pair['later_age']} | {pair['later_rate']} | {display(pair['minimum_decline'], 1)} |")
    lines.extend(["", "### Age-removal and truncation attempts", "",
                  "The full panel remains primary. These are diagnostic deletions, not a revised fitting panel. "
                  "All tied smallest deletion sets are listed; `{}` denotes the empty set. "
                  "Deleting the oldest bin equals the prefix ending at the preceding bin.", "",
                  "| Curve | Through 65-69 rho | Through 70-74 rho | Through 75-79 rho | Minimum deletion sets restoring compatibility |",
                  "|---|---:|---:|---:|---|"])
    for group in eligible:
        sens = group["age_sensitivity"]
        prefixes = {p["ending_age"]: p["rho"] for p in sens["prefixes"]}
        removals = "; ".join("{" + ", ".join(s) + "}" for s in sens["minimum_deletion_sets"])
        lines.append(f"| {_name(group)} | {prefixes['65-69']} | {prefixes['70-74']} | "
                     f"{prefixes['75-79']} | {removals} |")
    lines.extend(["", "<details>", "<summary>Every prefix and leave-one-age-out frontier</summary>", "",
                  "Vectors follow ages " + ", ".join(AGES) + ". Prefix entries end at the indicated age; "
                  "deletion entries remove only that age bin. All entries are rho.", "",
                  "| Curve | Prefix vector | Leave-one-age-out vector |", "|---|---|---|"])
    for group in eligible:
        sens = group["age_sensitivity"]
        lines.append(f"| {_name(group)} | {', '.join(x['rho'] for x in sens['prefixes'])} | "
                     f"{', '.join(x['rho'] for x in sens['leave_one_out'])} |")
    lines.extend(["", "</details>", "", "### Specified-shape calibration frontiers", "",
                  "These are finite-grid minima under uniform within-bin averaging. `L / U ages` name "
                  "the lower and upper amplitude constraints, in that order, not necessarily chronological order. "
                  "The grid index k is not a biological driver estimate. Negative g is net loss.", "",
                  "| Curve | Pure-power k | Pure rho | Pure L / U ages | One-growth k | g/year | Grid rho | Grid L / U ages | Unique grid minimum? |",
                  "|---|---:|---:|---|---:|---:|---:|---|---|"])
    for group in eligible:
        per_k = {item["k"]: item for item in group["by_k"]}
        power = per_k[group["best_power_k"]]["power"]
        grid = per_k[group["best_grid_k"]]["grid"]
        unique = len(group["best_grid_contenders"]) == 1
        lines.append(f"| {_name(group)} | {group['best_power_k']} | {power['rho']} | {_constraints(power)} | "
                     f"{group['best_grid_k']} | {display(grid['selected_g'], 3)} | {grid['rho']} | "
                     f"{_constraints(grid)} | {'yes' if unique else 'no; see contenders below'} |")
    unique_count = sum(len(group["best_grid_contenders"]) == 1 for group in eligible)
    age_unique_count = sum(len(item["grid"]["lower_ages"]) == len(item["grid"]["upper_ages"]) == 1
                           for group in eligible for item in group["by_k"])
    total_k_rows = sum(len(group["by_k"]) for group in eligible)
    lines.extend(["", f"Certified unique joint grid minima: **{unique_count}/{len(eligible)}**. "
                  f"Per-k selected-grid rows with certified single L and U ages: **{age_unique_count}/{total_k_rows}**. "
                  "Every reported grid threshold has a rational squared-threshold enclosure narrower than "
                  "`1e-24`; both bounds round to the displayed rho. Full rational enclosures are emitted by the CLI.", "",
                  "<details>", "<summary>Complete per-k frontiers, bin-weight sensitivity and constraint ages</summary>", "",
                  "`Weight rho` is the pure-power infimum over freely varying within-bin age weights. "
                  "It is not a confidence bound or a guaranteed attained fit. A list in an age cell "
                  "retains tied/unresolved constraints. g contenders are retained at every k.", "",
                  "| Curve | k | Pure rho | Pure L / U ages | Weight rho | Weight L / U ages | Grid rho | g/year | g contenders | Grid L / U ages |",
                  "|---|---:|---:|---|---:|---|---:|---:|---|---|"])
    for group in eligible:
        for item in group["by_k"]:
            power, flexible, grid = item["power"], item["flexible_power"], item["grid"]
            contenders = ", ".join(display(g, 3) for g in grid["contenders"])
            lines.append(f"| {_name(group)} | {item['k']} | {power['rho']} | {_constraints(power)} | "
                         f"{flexible['rho']} | {_constraints(flexible)} | {grid['rho']} | "
                         f"{display(grid['selected_g'], 3)} | {contenders} | {_constraints(grid)} |")
    lines.extend(["", "Joint-grid contender sets and exact-power ties:", "",
                  "| Curve | Pure-power k contenders | Joint (k, g/year) contenders |", "|---|---|---|"])
    for group in eligible:
        contenders = "; ".join(f"({c['k']}, {display(c['g'], 3)})" for c in group["best_grid_contenders"])
        lines.append(f"| {_name(group)} | {', '.join(map(str, group['best_power_contenders']))} | {contenders} |")
    lines.extend(["", "</details>", "", "### Content integrity", "",
                  "SHA-256 after UTF-8 universal-newline normalization; these establish local content "
                  "integrity, not upstream source truth. Relative path keys use `Path.as_posix()`.", "",
                  "| Input | Normalized SHA-256 |", "|---|---|"])
    for path, digest in sorted(results["hashes"].items()):
        lines.append(f"| `{path}` | `{digest}` |")
    method = results["inherited_method"]
    lines.extend(["", f"Unmodified inherited primitive file `{method['path']}`: `{method['sha256']}`.",
                  "", END])
    return "\n".join(lines)


def json_ready(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="Check the complete generated README block")
    mode.add_argument("--write-readme", action="store_true", help="Replace only the generated README block")
    args = parser.parse_args()
    result = analyze()
    if not (args.check or args.write_readme):
        print(json.dumps(json_ready(result), indent=2, sort_keys=True))
        return
    readme = ROOT / "README.md"
    content = readme.read_text(encoding="utf-8")
    if content.count(BEGIN) != 1 or content.count(END) != 1:
        raise ValueError("Require exactly one results block")
    before, tail = content.split(BEGIN)
    _, after = tail.split(END)
    updated = before + render_results(result) + after
    if args.check:
        if updated != content:
            raise SystemExit("README numerical results differ from recomputation")
        print("All README results and input hashes recomputed successfully. " + LABEL)
    else:
        readme.write_text(updated, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
