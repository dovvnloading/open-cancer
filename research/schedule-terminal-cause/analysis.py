"""Offline nonclinical terminal-attribution analysis of printed aggregates.

No treatment advice. Exact conditional reference tests and deterministic
attribution envelopes are not causal, survival, or mechanistic estimates.
Only Python's standard library is used; this module performs no network I/O.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROTOCOL_COMMIT = "dbf3c54746ba853a81a65f177108531081460978"
PROTOCOL_SHA256 = "0f228ecb01a7d187b37a57f43d163fb8b4cedb954ee0258c24ccd94363283b90"
ALPHA = Fraction(1, 20)
CELL_ALPHA = ALPHA / 12
CAUSES = ("toxicity", "cancer", "other")
S4_GROUPS = ("vehicle", "ping_pong", "tandem", "combo_mtd")
LABELS = {
    "vehicle": "Vehicle", "ping_pong": "Ping-pong adaptive",
    "tandem": "Tandem adaptive", "combo_mtd": "Combination MTD",
    "cap_adaptive": "CAP adaptive", "cap_mtd": "CAP MTD",
    "gem_adaptive": "GEM adaptive", "gem_mtd": "GEM MTD",
}
PAIRS = (
    ("mtd_ping", 3, 1, "MTD minus ping-pong (primary)"),
    ("mtd_tandem", 3, 2, "MTD minus tandem"),
    ("mtd_vehicle", 3, 0, "MTD minus vehicle"),
    ("tandem_ping", 2, 1, "Tandem minus ping-pong"),
    ("ping_vehicle", 1, 0, "Ping-pong minus vehicle"),
    ("tandem_vehicle", 2, 0, "Tandem minus vehicle"),
)
BREAKPOINTS = ("raw_composition", "raw_toxicity", "holm_composition",
               "holm_toxicity", "direction", "combined")
BEGIN = "<!-- BEGIN GENERATED RESULTS -->"
END = "<!-- END GENERATED RESULTS -->"


def normalized_sha256(path):
    """Normalize newline bytes only; preserve all other bytes and final LF."""
    content = Path(path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(content).hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _bounded_rows(columns, total, prefix=()):
    """All first rows with fixed size and nonnegative complementary row."""
    if len(columns) == 1:
        if 0 <= total <= columns[0]:
            yield prefix + (total,)
        return
    lower = max(0, total - sum(columns[1:]))
    for count in range(lower, min(columns[0], total) + 1):
        yield from _bounded_rows(columns[1:], total - count, prefix + (count,))


@lru_cache(maxsize=None)
def exact_conditional(first, second):
    """Probability-ordering Fisher(-Freeman-Halton) exact two-sided tail.

    Conditional multivariate hypergeometric numerators are integers with a
    common denominator. Tail membership and probability ties are exact.
    """
    if len(first) != len(second) or not first:
        raise ValueError("Both rows must have the same positive number of columns")
    if any(type(x) is not int or x < 0 for x in first + second):
        raise ValueError("Counts must be nonnegative integers")
    n_first, n_second = sum(first), sum(second)
    if not n_first or not n_second:
        raise ValueError("Both row totals must be positive")
    columns = tuple(a + b for a, b in zip(first, second))
    observed = math.prod(math.comb(c, a) for c, a in zip(columns, first))
    denominator = math.comb(n_first + n_second, n_first)
    tail = 0
    for row in _bounded_rows(columns, n_first):
        numerator = math.prod(math.comb(c, x) for c, x in zip(columns, row))
        if numerator <= observed:
            tail += numerator
    return Fraction(tail, denominator)


@lru_cache(maxsize=None)
def pair_tests(first, second):
    return (
        exact_conditional(first, second),
        exact_conditional((first[0], sum(first[1:])),
                          (second[0], sum(second[1:]))),
    )


def holm(pvalues):
    """Step-down adjusted values, valid with arbitrary test dependence."""
    pvalues = tuple(pvalues)
    adjusted = [Fraction(0)] * len(pvalues)
    running = Fraction(0)
    for rank, index in enumerate(sorted(range(len(pvalues)), key=pvalues.__getitem__)):
        running = max(running, min(Fraction(1), (len(pvalues) - rank) * pvalues[index]))
        adjusted[index] = running
    return tuple(adjusted)


def _binomial_tail(k, n, probability, upper):
    values = range(k, n + 1) if upper else range(0, k + 1)
    return math.fsum(math.comb(n, j) * probability ** j
                     * (1.0 - probability) ** (n - j) for j in values)


@lru_cache(maxsize=None)
def clopper_pearson(k, n, alpha=CELL_ALPHA):
    """Two-sided exact binomial bounds, deterministic inversion to float precision."""
    if not 0 <= k <= n or n <= 0 or not 0 < alpha < 1:
        raise ValueError("Invalid binomial interval arguments")
    tail = float(alpha / 2)
    lower, upper = 0.0, 1.0
    if k:
        lo, hi = 0.0, 1.0
        for _ in range(64):
            mid = (lo + hi) / 2
            if _binomial_tail(k, n, mid, True) < tail:
                lo = mid
            else:
                hi = mid
        lower = (lo + hi) / 2
    if k < n:
        lo, hi = 0.0, 1.0
        for _ in range(64):
            mid = (lo + hi) / 2
            if _binomial_tail(k, n, mid, False) > tail:
                lo = mid
            else:
                hi = mid
        upper = (lo + hi) / 2
    return lower, upper


@lru_cache(maxsize=None)
def _pair_metrics(first, second):
    n1, n2 = sum(first), sum(second)
    p_comp, p_tox = pair_tests(first, second)
    result = {"composition_p": p_comp, "toxicity_p": p_tox}
    for category, a, b in zip(CAUSES, first, second):
        result[category + "_difference"] = Fraction(a, n1) - Fraction(b, n2)
        lo1, hi1 = clopper_pearson(a, n1)
        lo2, hi2 = clopper_pearson(b, n2)
        result[category + "_interval_lower"] = lo1 - hi2
        result[category + "_interval_upper"] = hi1 - lo2
    return result


def evaluate_s4(state):
    """Evaluate all 12 tests on one coherent table with the shared control once."""
    records, pvalues = [], []
    for pair_id, i, j, label in PAIRS:
        metrics = dict(_pair_metrics(state[i], state[j]))
        pvalues.extend((metrics["composition_p"], metrics["toxicity_p"]))
        records.append({"id": pair_id, "first": S4_GROUPS[i],
                        "second": S4_GROUPS[j], "label": label, "metrics": metrics})
    adjusted = holm(pvalues)
    for index, record in enumerate(records):
        record["metrics"]["composition_holm"] = adjusted[2 * index]
        record["metrics"]["toxicity_holm"] = adjusted[2 * index + 1]
    return records


@lru_cache(maxsize=None)
def row_states(counts, include_other):
    """Each state is (new counts, changed toxicity count, changed other count)."""
    t, c, o = counts
    return tuple(((t - dt, c + dt + do, o - do), dt, do)
                 for dt in range(t + 1)
                 for do in range((o if include_other else 0) + 1))


def _update_ranges(ranges, metrics, witness):
    for metric, value in metrics.items():
        if metric not in ranges:
            ranges[metric] = {"min": value, "max": value,
                              "min_witness": witness, "max_witness": witness}
        else:
            if value < ranges[metric]["min"]:
                ranges[metric]["min"], ranges[metric]["min_witness"] = value, witness
            if value > ranges[metric]["max"]:
                ranges[metric]["max"], ranges[metric]["max_witness"] = value, witness


def _lost_conditions(metrics, baseline_difference):
    direction_lost = baseline_difference * metrics["toxicity_difference"] <= 0
    return {
        "raw_composition": metrics["composition_p"] > ALPHA,
        "raw_toxicity": metrics["toxicity_p"] > ALPHA,
        "holm_composition": metrics["composition_holm"] > ALPHA,
        "holm_toxicity": metrics["toxicity_holm"] > ALPHA,
        "direction": direction_lost,
        "combined": direction_lost or metrics["composition_holm"] > ALPHA
                    or metrics["toxicity_holm"] > ALPHA,
    }


def _update_breakpoints(records, conditions, k, total_n, eligible, witness):
    for name, lost in conditions.items():
        if lost and (records[name] is None or k < records[name]["k"]):
            records[name] = {
                "k": k, "fraction_animals": Fraction(k, total_n),
                "fraction_eligible": Fraction(k, eligible) if eligible else Fraction(0),
                "witness": witness,
            }


def s4_envelope(baseline_state, include_other):
    """Exhaustive statewise multiplicity correction, bounds and local fragility."""
    total_n = sum(map(sum, baseline_state))
    eligible = sum(row[0] + (row[2] if include_other else 0) for row in baseline_state)
    baseline = evaluate_s4(baseline_state)
    directions = {r["id"]: r["metrics"]["toxicity_difference"] for r in baseline}
    ranges = {p[0]: {} for p in PAIRS}
    local = {p[0]: dict.fromkeys(BREAKPOINTS) for p in PAIRS}
    global_breaks = dict.fromkeys(BREAKPOINTS)
    exact_budget = [{"states": 0, "pairs": {p[0]: {} for p in PAIRS}}
                    for _ in range(eligible + 1)]
    state_count = 0
    for choices in itertools.product(*(row_states(row, include_other)
                                       for row in baseline_state)):
        state = tuple(choice[0] for choice in choices)
        moves = tuple(choice[1] + choice[2] for choice in choices)
        k = sum(moves)
        evaluated = evaluate_s4(state)
        exact_budget[k]["states"] += 1
        state_count += 1
        for record, (_, i, j, _) in zip(evaluated, PAIRS):
            pair_id, metrics = record["id"], record["metrics"]
            _update_ranges(ranges[pair_id], metrics, state)
            # Budget envelopes retain every metric and its attained witness.
            _update_ranges(exact_budget[k]["pairs"][pair_id], metrics, state)
            conditions = _lost_conditions(metrics, directions[pair_id])
            if all(moves[index] == 0 for index in range(4) if index not in (i, j)):
                pair_n = sum(baseline_state[i]) + sum(baseline_state[j])
                pair_eligible = sum(baseline_state[index][0]
                                    + (baseline_state[index][2] if include_other else 0)
                                    for index in (i, j))
                _update_breakpoints(local[pair_id], conditions, k, pair_n,
                                    pair_eligible, state)
            if pair_id == PAIRS[0][0]:
                _update_breakpoints(global_breaks, conditions, k, total_n, eligible, state)

    # Convert exact-k ranges to attained at-most-k envelopes. Witnesses are
    # copied with the extrema, so no unattained coordinate-wise table is implied.
    cumulative = {p[0]: {} for p in PAIRS}
    budgets, cumulative_count = [], 0
    for k, bucket in enumerate(exact_budget):
        cumulative_count += bucket["states"]
        for pair_id, bounds in bucket["pairs"].items():
            for metric, bound in bounds.items():
                for side in ("min", "max"):
                    current = cumulative[pair_id].get(metric)
                    value = bound[side]
                    if current is None:
                        cumulative[pair_id][metric] = dict(bound)
                    elif (side == "min" and value < current[side]
                          or side == "max" and value > current[side]):
                        current[side] = value
                        current[side + "_witness"] = bound[side + "_witness"]
        budgets.append({
            "k": k, "states": cumulative_count,
            "fraction_animals": Fraction(k, total_n),
            "fraction_eligible": Fraction(k, eligible) if eligible else Fraction(0),
            "pairs": {pair_id: {metric: dict(bound) for metric, bound in values.items()}
                      for pair_id, values in cumulative.items()},
        })
    verdicts = {pair_id: "survives" if values["combined"] is None else "fails"
                for pair_id, values in local.items()}
    # The global combined breakpoint also covers multiplicity effects caused
    # by edits outside a pair. Local failure implies global failure, but local
    # absence alone would not prove global robustness.
    for pair_id, bounds in ranges.items():
        direction = directions[pair_id]
        same_direction = (bounds["toxicity_difference"]["min"] > 0 if direction > 0
                          else bounds["toxicity_difference"]["max"] < 0
                          if direction < 0 else False)
        verdicts[pair_id] = ("survives" if same_direction
                             and bounds["composition_holm"]["max"] <= ALPHA
                             and bounds["toxicity_holm"]["max"] <= ALPHA else "fails")
    return {
        "state_count": state_count, "eligible_labels": eligible,
        "total_n": total_n, "ranges": ranges, "verdicts": verdicts,
        "pair_local_breakpoints": local,
        "primary_global_breakpoints": global_breaks, "budgets": budgets,
    }


def single_drug_bounds(adaptive, mtd):
    """Bounds conditional on exactly one excess printed adaptive cause count."""
    completions = tuple(tuple(x - (i == category) for i, x in enumerate(adaptive))
                        for category in range(3) if adaptive[category] > 0)
    bounds = {}
    for model, include_other in (("as_printed_completion", None),
                                 ("toxicity_only", False),
                                 ("toxicity_and_other", True)):
        ranges, states = {}, set()
        scenarios = 0
        for completion in completions:
            a_choices = ((completion, 0, 0),) if include_other is None else row_states(completion, include_other)
            m_choices = ((mtd, 0, 0),) if include_other is None else row_states(mtd, include_other)
            for a_choice, m_choice in itertools.product(a_choices, m_choices):
                state = a_choice[0], m_choice[0]
                comp, tox = pair_tests(*state)
                metrics = {category + "_difference": Fraction(state[0][index], sum(state[0]))
                           - Fraction(state[1][index], sum(state[1]))
                           for index, category in enumerate(CAUSES)}
                metrics.update(composition_p=comp, toxicity_p=tox,
                               composition_bonferroni=min(Fraction(1), 4 * comp),
                               toxicity_bonferroni=min(Fraction(1), 4 * tox))
                witness = {"completion": completion, "state": state}
                _update_ranges(ranges, metrics, witness)
                scenarios += 1
                states.add(state)
        bounds[model] = {"ranges": ranges, "completion_scenarios": scenarios,
                         "unique_states": len(states)}
    return {"unconditional_verdict": "untestable", "completions": completions,
            "printed_adaptive": adaptive, "printed_mtd": mtd, "bounds": bounds}


def build_report(base_dir=HERE):
    """Recompute every quantitative report field directly from the small CSVs."""
    base_dir = Path(base_dir)
    counts_rows = read_csv(base_dir / "data" / "counts.csv")
    allocation_rows = read_csv(base_dir / "data" / "allocations.csv")
    dose_rows = read_csv(base_dir / "data" / "doses.csv")
    conflicts = read_csv(base_dir / "data" / "dose_conflicts.csv")
    counts = {row["group_id"]: tuple(int(row[c]) for c in CAUSES) for row in counts_rows}
    if len(counts) != len(counts_rows):
        raise ValueError("Duplicate pooled cause row")
    if len({row["arm_id"] for row in allocation_rows}) != len(allocation_rows):
        raise ValueError("Duplicate allocated arm")
    allocations = {}
    for row in allocation_rows:
        group = row["pooled_group_id"]
        allocations[group] = allocations.get(group, 0) + int(row["animals"])
    denominators = [{"group_id": row["group_id"], "counts": counts[row["group_id"]],
                     "allocated": allocations[row["group_id"]],
                     "printed": sum(counts[row["group_id"]]),
                     "excess": sum(counts[row["group_id"]]) - allocations[row["group_id"]],
                     "table": row["source_table"], "page": row["pdf_page"]}
                    for row in counts_rows]
    state = tuple(counts[group] for group in S4_GROUPS)
    s4_ok = all(sum(counts[group]) == allocations[group] for group in S4_GROUPS)
    if not s4_ok:
        raise ValueError("S4 denominator gate fails: no inferential report may be generated")
    for drug in ("cap", "gem"):
        if sum(counts[drug + "_adaptive"]) - allocations[drug + "_adaptive"] != 1:
            raise ValueError("The frozen exactly-one-excess conditional model does not apply")
        if sum(counts[drug + "_mtd"]) != allocations[drug + "_mtd"]:
            raise ValueError("A single-drug MTD denominator does not reconcile")
    baseline = evaluate_s4(state)
    envelopes = {"toxicity_only": s4_envelope(state, False),
                 "toxicity_and_other": s4_envelope(state, True)}
    arms = [{"group_id": group, "counts": row, "n": sum(row),
             "proportions": tuple(Fraction(x, sum(row)) for x in row),
             "intervals": tuple(clopper_pearson(x, sum(row)) for x in row)}
            for group, row in zip(S4_GROUPS, state)]
    hash_paths = [base_dir / "PROTOCOL.md", base_dir / "analysis.py"]
    hash_paths += sorted((base_dir / "data").glob("*.csv"))
    hashes = {p.relative_to(base_dir).as_posix(): normalized_sha256(p) for p in hash_paths}
    if hashes[Path("PROTOCOL.md").as_posix()] != PROTOCOL_SHA256:
        raise ValueError("Frozen protocol hash mismatch")
    return {
        "protocol_commit": PROTOCOL_COMMIT, "alpha": ALPHA, "cell_alpha": CELL_ALPHA,
        "hashes": hashes, "counts_rows": counts_rows, "allocation_rows": allocation_rows,
        "denominators": denominators, "allocated_total": sum(allocations.values()),
        "printed_total": sum(map(sum, counts.values())),
        "s4": {"state": state, "arms": arms, "baseline": baseline, "envelopes": envelopes,
               "total_n": sum(map(sum, state)),
               "verdict": envelopes["toxicity_and_other"]["verdicts"][PAIRS[0][0]],
               "printed_primary_pass": not _lost_conditions(
                   baseline[0]["metrics"], baseline[0]["metrics"]["toxicity_difference"])["combined"]},
        "singles": {drug: single_drug_bounds(counts[drug + "_adaptive"], counts[drug + "_mtd"])
                    for drug in ("cap", "gem")},
        "doses": dose_rows, "dose_conflicts": conflicts,
    }


def _number(value, percent=False):
    return f"{float(value) * (100 if percent else 1):.6f}"


def _range(bounds, percent=False):
    return f"[{_number(bounds['min'], percent)}, {_number(bounds['max'], percent)}]"


def _counts(row):
    return "(" + ",".join(map(str, row)) + ")"


def _witness(state):
    if isinstance(state, dict):
        return "decrement completion " + _counts(state["completion"]) + "; A/M=" + "/".join(map(_counts, state["state"]))
    return "/".join(map(_counts, state))


def _fraction_label(value, count, denominator):
    # Preserve the actual animal/eligible-label denominator, including at k=0.
    # A reduced Fraction such as 0/1 or 1/9 obscures the sample size here.
    return f"{count}/{denominator} ({100 * float(value):.6f}%)"


def _table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "| " + " | ".join("---" for _ in headers) + " |"]
                     + ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def render_report(report):
    """Generated Markdown; all displayed result numbers derive from build_report."""
    s4, pieces = report["s4"], []
    pieces += ["## Recomputed results", "",
               f"**Primary verdict: {s4['verdict'].upper()}.** The printed-label primary rule "
               + ("passes" if s4["printed_primary_pass"] else "already fails")
               + "; attribution uncertainty is evaluated below. Failure is failure of the frozen "
               "robustness rule, not evidence of no toxicity or no association.", "",
               "All differences are percentage points (pp). Probabilities are displayed to six "
               "decimals; calculations and decision comparisons use exact rational p-values. "
               "T/C/O means toxicity/cancer/other attributed terminal events.", "",
               "### Denominator audit", "",
               _table(["Source", "Group", "Printed T/C/O", "Printed total", "Allocated", "Excess"],
                      [(f"{r['table']} p.{r['page']}", LABELS[r["group_id"]], _counts(r["counts"]),
                        r["printed"], r["allocated"], r["excess"]) for r in report["denominators"]]), "",
               f"Across the unique printed rows: {report['printed_total']} assignments versus "
               f"{report['allocated_total']} allocated animals. S4 uses {s4['total_n']} unique animals; "
               "its MTD and vehicle rows occur once. S2/S3 remain unreconciled.", "",
               "### S4 compositions and simultaneous reference intervals", "",
               "Intervals have at least 95% simultaneous coverage for the twelve category "
               "probabilities under the frozen exchangeable-binomial reference model; "
               "each marginal uses error " + str(report["cell_alpha"]) + ".", "",
               _table(["Group", "n", "T/C/O (%)", "T interval (%)", "C interval (%)", "O interval (%)"],
                      [(LABELS[a["group_id"]], a["n"], "/".join(_number(x, True) for x in a["proportions"]),
                        *[f"[{_number(lo, True)}, {_number(hi, True)}]" for lo, hi in a["intervals"]])
                       for a in s4["arms"]]), "",
               "### All fixed S4 comparisons at printed labels", "",
               _table(["First minus second", "T difference (pp)", "C difference (pp)", "O difference (pp)",
                       "Composition p", "Holm p", "T vs rest p", "Holm p"],
                      [(r["label"], *[_number(r["metrics"][c + "_difference"], True) for c in CAUSES],
                        *[_number(r["metrics"][k]) for k in
                          ("composition_p", "composition_holm", "toxicity_p", "toxicity_holm")])
                       for r in s4["baseline"]]), "",
               _table(["First minus second", "T difference interval (pp)", "C difference interval (pp)",
                       "O difference interval (pp)"],
                      [(r["label"], *[f"[{_number(r['metrics'][c + '_interval_lower'], True)}, "
                                       f"{_number(r['metrics'][c + '_interval_upper'], True)}]" for c in CAUSES])
                       for r in s4["baseline"]]), "",
               "The intervals are conservative simultaneous bounds, not inversions of the "
               "Holm tests. Shared comparisons are dependent; no independent-replication claim is made.", ""]

    labels = {p[0]: p[3] for p in PAIRS}
    for name, envelope in s4["envelopes"].items():
        heading = "T-only" if name == "toxicity_only" else "T-and-O"
        pieces += [f"### Attribution envelope: {heading}", "",
                   f"{envelope['state_count']} coherent integer states; {envelope['eligible_labels']} "
                   f"eligible labels among {envelope['total_n']} animals. The primary verdict is "
                   f"**{envelope['verdicts'][PAIRS[0][0]]}**. These are deterministic sensitivity "
                   "ranges, not confidence distributions.", "",
                   _table(["Contrast", "T difference range (pp)", "C difference range (pp)",
                           "O difference range (pp)", "Robustness verdict"],
                          [(labels[pid], *[_range(bounds[c + "_difference"], True) for c in CAUSES],
                            envelope["verdicts"][pid]) for pid, bounds in envelope["ranges"].items()]), "",
                   _table(["Contrast", "Composition p range", "Composition Holm range", "T p range", "T Holm range"],
                          [(labels[pid], *[_range(bounds[k]) for k in
                            ("composition_p", "composition_holm", "toxicity_p", "toxicity_holm")])
                           for pid, bounds in envelope["ranges"].items()]), "",
                   "Envelope of simultaneous difference intervals (minimum lower bound to maximum upper bound):", "",
                   _table(["Contrast", "T (pp)", "C (pp)", "O (pp)"],
                          [(labels[pid], *[f"[{_number(bounds[c + '_interval_lower']['min'], True)}, "
                                           f"{_number(bounds[c + '_interval_upper']['max'], True)}]"
                                          for c in CAUSES])
                           for pid, bounds in envelope["ranges"].items()]), "",
                   "#### Pair-local fragility", "",
                   "Other rows remain at their printed labels. Breakpoint zero means the printed "
                   "comparison already lacks rejection or strict direction. Fractions below use "
                   "the pair's animals and eligible labels, respectively. Witness row order is "
                   "vehicle / ping-pong / tandem / combination MTD, each T,C,O.", ""]
        break_rows = []
        for pid, records in envelope["pair_local_breakpoints"].items():
            _, first, second, _ = next(pair for pair in PAIRS if pair[0] == pid)
            pair_n = sum(s4["state"][first]) + sum(s4["state"][second])
            pair_eligible = sum(s4["state"][index][0]
                                + (s4["state"][index][2]
                                   if name == "toxicity_and_other" else 0)
                                for index in (first, second))
            for criterion, item in records.items():
                break_rows.append((labels[pid], criterion,
                    "unreachable" if item is None else item["k"],
                    "-" if item is None else _fraction_label(item["fraction_animals"], item["k"], pair_n),
                    "-" if item is None else _fraction_label(item["fraction_eligible"], item["k"], pair_eligible),
                    "-" if item is None else _witness(item["witness"])))
        pieces += [_table(["Contrast", "Loss criterion", "Relabels", "Fraction of pair animals",
                           "Fraction of eligible pair labels", "Witness table"], break_rows), "",
                   "#### Primary coherent global fragility", "",
                   "Here edits may affect any row, including changes outside the primary pair that "
                   "alter the Holm ordering. Fractions use the whole S4 table.", "",
                   _table(["Loss criterion", "Relabels", "Fraction of all animals", "Fraction of all eligible labels", "Witness table"],
                          [(criterion, "unreachable" if item is None else item["k"],
                            "-" if item is None else _fraction_label(item["fraction_animals"], item["k"], envelope["total_n"]),
                            "-" if item is None else _fraction_label(item["fraction_eligible"], item["k"], envelope["eligible_labels"]),
                            "-" if item is None else _witness(item["witness"]))
                           for criterion, item in envelope["primary_global_breakpoints"].items()]), "",
                   "#### Primary worst case at every cumulative global budget", "",
                   "Each row permits at most k changed labels; extrema can be attained by different "
                   "tables and are not a synthetic joint table. Complete ranges and attaining "
                   "witnesses for every metric and every pair at each budget are available with `--json`.", "",
                   _table(["k", "States", "k/all animals (%)", "k/eligible (%)", "T difference range (pp)",
                           "Worst composition p", "Worst composition Holm", "Worst T p", "Worst T Holm"],
                          [(b["k"], b["states"], _number(b["fraction_animals"], True),
                            _number(b["fraction_eligible"], True),
                            _range(b["pairs"][PAIRS[0][0]]["toxicity_difference"], True),
                            *[_number(b["pairs"][PAIRS[0][0]][k]["max"]) for k in
                              ("composition_p", "composition_holm", "toxicity_p", "toxicity_holm")])
                           for b in envelope["budgets"]]), "",
                   "<details>", f"<summary>Attaining witnesses for all {heading} extrema</summary>", "",
                   "Differences and interval endpoints are pp; p-values are probabilities. Each "
                   "extremum has its own coherent table. Interval lower/upper extrema are "
                   "endpoints of reference-model uncertainty, not additional tests.", "",
                   _table(["Contrast", "Metric", "Minimum", "Minimum witness", "Maximum", "Maximum witness"],
                          [(labels[pid], metric, _number(bound["min"], "difference" in metric or "interval" in metric),
                            _witness(bound["min_witness"]),
                            _number(bound["max"], "difference" in metric or "interval" in metric),
                            _witness(bound["max_witness"]))
                           for pid, bounds in envelope["ranges"].items() for metric, bound in bounds.items()]), "",
                   "</details>", ""]

    pieces += ["### Single-drug denominator sensitivity", "",
               "Both unconditional single-drug verdicts are **untestable**. The following bounds "
               "assume exactly one excess adaptive category assignment; they are not source "
               "corrections. Differences here are adaptive minus its own MTD, unlike the primary "
               "MTD-minus-adaptive direction. Bonferroni values use the separate four-test "
               "conditional family and cannot rescue the S4 verdict.", ""]
    for drug, analysis in report["singles"].items():
        pieces += [f"#### {drug.upper()}", "",
                   "All possible one-cell-decrement completions: "
                   + ", ".join(_counts(c) for c in analysis["completions"]) + ".", "",
                   _table(["Model", "Completion/state scenarios", "Unique tables", "T difference (pp)",
                           "C difference (pp)", "O difference (pp)", "Composition p", "T p",
                           "Composition Bonferroni", "T Bonferroni"],
                          [(model, result["completion_scenarios"], result["unique_states"],
                            *[_range(result["ranges"][c + "_difference"], True) for c in CAUSES],
                            *[_range(result["ranges"][m]) for m in
                              ("composition_p", "toxicity_p", "composition_bonferroni", "toxicity_bonferroni")])
                           for model, result in analysis["bounds"].items()]), "",
                   "<details>", "<summary>Conditional single-drug attaining witnesses</summary>", "",
                   _table(["Model", "Metric", "Minimum", "Minimum witness", "Maximum", "Maximum witness"],
                          [(model, metric, _number(bound["min"], "difference" in metric),
                            _witness(bound["min_witness"]), _number(bound["max"], "difference" in metric),
                            _witness(bound["max_witness"]))
                           for model, result in analysis["bounds"].items()
                           for metric, bound in result["ranges"].items()]), "", "</details>", ""]

    pieces += ["### Printed dose context and source disagreements", "",
               "S5 means are not used in any cause test, follow-up reconstruction, or dose-response "
               "adjustment. Daily units follow the column heading; the printed cells say mg/kg.", "",
               _table(["S5 source row", "Drug", "Mean cumulative mg/kg", "Mean daily mg/kg/day"],
                      [(r["source_row"], r["drug"], r["mean_cumulative_mg_kg"], r["mean_daily_mg_kg_day"])
                       for r in report["doses"]]), ""]
    dose_lookup = {(r["arm_id"], r["drug"]): r for r in report["doses"]}
    pieces += [_table(["Arm", "Drug", "S5 cumulative mean", "Article cumulative mean", "Disposition"],
                      [(r["arm_id"], r["drug"], dose_lookup[r["arm_id"], r["drug"]]["mean_cumulative_mg_kg"],
                        r["mean_cumulative_mg_kg"], "Unresolved; both retained") for r in report["dose_conflicts"]]), "",
               "### Reproducibility identifiers", "",
               f"Original local protocol-only freeze commit: `{report['protocol_commit']}` "
               "(publication mapping above). Text hashes below normalize "
               "CRLF and CR to LF and preserve all other bytes and final-newline state. Paths use "
               "`Path.as_posix()`. These hashes attest content consistency, not source truth.", "",
               _table(["Relative path", "SHA-256 of newline-normalised content"],
                      [(f"`{p}`", f"`{digest}`") for p, digest in report["hashes"].items()])]
    return "\n".join(pieces).rstrip() + "\n"


def _json_default(value):
    if isinstance(value, Fraction):
        return {"numerator": value.numerator, "denominator": value.denominator,
                "decimal": float(value)}
    raise TypeError(type(value).__name__)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Print exact report and all attaining witnesses")
    parser.add_argument("--check", action="store_true", help="Check README's generated block without modifying files")
    args = parser.parse_args(argv)
    report = build_report()
    if args.json:
        print(json.dumps(report, default=_json_default, indent=2, sort_keys=True))
    elif args.check:
        readme = (HERE / "README.md").read_text(encoding="utf-8")
        if readme.count(BEGIN) != 1 or readme.count(END) != 1:
            raise SystemExit("README must have exactly one generated-results block")
        actual = readme.split(BEGIN, 1)[1].split(END, 1)[0].strip()
        if actual != render_report(report).strip():
            raise SystemExit("README generated results do not match CSV recomputation")
        print("Generated nonclinical results match CSV recomputation.")
    else:
        print(render_report(report), end="")


if __name__ == "__main__":
    main()
