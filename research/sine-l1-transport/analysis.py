"""Frozen job-5 whole-order SINE-vs-L1 binomial transport audit.

Nonclinical, unreviewed research; not medical or veterinary advice.
Only Python stdlib. Offline analysis uses narrowly derived species aggregates;
source CSV is never fetched or mirrored. See PROTOCOL.md and source_audit.json.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PANEL = HERE / "data" / "species_aggregates.csv"
LAMBDA = 1.0
T95_10 = 2.2281388519649385


class AnalysisError(ValueError):
    """Frozen gate or numerical fit failed: never silently drop a fold."""


def normalized_hash(path: Path) -> str:
    """Hash UTF-8 file contents after canonicalising CRLF and CR to LF."""
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


def load_panel(path: Path = PANEL) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        source = list(csv.DictReader(f))
    if len(source) != 55 or len({r["order"] for r in source}) != 11:
        raise AnalysisError("55-species / 11-order source panel gate failed")
    if [r["species"] for r in source] != sorted(r["species"] for r in source):
        raise AnalysisError("species order is not canonical")
    out = []
    for r in source:
        n, k, mal = (int(r[key]) for key in ("n", "neo_k", "mal_k"))
        nc = tuple(int(x) for x in r["neo_candidates"].split("|"))
        mc = tuple(int(x) for x in r["mal_candidates"].split("|"))
        x = tuple(float(r[key]) for key in ("log1p_l1", "longevity_y", "busco", "log1p_n50", "log1p_sine"))
        if not (n > 0 and 0 <= k <= n and 0 <= mal <= k and k in nc and mal in mc):
            raise AnalysisError("denominator/numerator gate failed: " + r["species"])
        if not (all(0 <= a <= n for a in nc + mc) and 0 <= x[2] <= 100
                and x[1] > 0 and x[3] > 0 and all(map(math.isfinite, x))):
            raise AnalysisError("predictor/count gate failed: " + r["species"])
        if r["ratio_ok"] not in ("0", "1"):
            raise AnalysisError("ratio audit indicator absent")
        out.append({"species": r["species"], "order": r["order"], "n": n,
                    "neo_k": k, "mal_k": mal, "neo_candidates": nc,
                    "mal_candidates": mc, "ratio_ok": bool(int(r["ratio_ok"])), "x": x})
    if len({r["species"] for r in out}) != 55:
        raise AnalysisError("duplicate species identifier")
    if sum(r["n"] for r in out) != 3533 or sum(r["neo_k"] for r in out) != 691:
        raise AnalysisError("pinned denominator/numerator nominal checksum failed")
    if sum(r["mal_k"] for r in out) != 405:
        raise AnalysisError("pinned secondary nominal checksum failed")
    if sum(not r["ratio_ok"] for r in out) != 31:
        raise AnalysisError("secondary provenance/ratio gate changed")
    return out


def softplus(x: float) -> float:
    return x + math.log1p(math.exp(-x)) if x > 0 else math.log1p(math.exp(x))


def sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def cholesky_solve(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Positive definite Newton solve; no numpy/scipy/platform-dependent solver."""
    m = len(rhs)
    l = [[0.0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1):
            v = matrix[i][j] - sum(l[i][k] * l[j][k] for k in range(j))
            if i == j:
                if v <= 0 or not math.isfinite(v):
                    raise AnalysisError("non positive-definite Hessian")
                l[i][j] = math.sqrt(v)
            else:
                l[i][j] = v / l[j][j]
    y = [0.0] * m
    for i in range(m):
        y[i] = (rhs[i] - sum(l[i][j] * y[j] for j in range(i))) / l[i][i]
    out = [0.0] * m
    for i in range(m - 1, -1, -1):
        out[i] = (y[i] - sum(l[j][i] * out[j] for j in range(i + 1, m))) / l[i][i]
    return out


def quantile(values: list[float], proportion: float) -> float:
    s = sorted(values)
    v = (len(s) - 1) * proportion
    i = int(v)
    return s[i] + (s[min(i + 1, len(s) - 1)] - s[i]) * (v - i)


def transform(train: list[dict], held: list[dict], extended: bool, winsor: bool):
    """All feature transforms/cutoffs/centering learned on training species."""
    indexes = range(5 if extended else 4)
    bounds = None
    if winsor:
        vals = [r["x"][3] for r in train]
        bounds = (quantile(vals, 0.1), quantile(vals, 0.9))

    def raw(r):
        z = list(r["x"])
        if bounds is not None:
            z[3] = min(bounds[1], max(bounds[0], z[3]))
        return [z[j] for j in indexes]

    a, b = [raw(r) for r in train], [raw(r) for r in held]
    centers = [sum(x[j] for x in a) / len(a) for j in indexes]
    scales = [math.sqrt(sum((x[j] - centers[j]) ** 2 for x in a) / len(a)) or 1.0
              for j in indexes]
    def standardize(xx):
        return [tuple([1.0] + [(z[j] - centers[j]) / scales[j] for j in indexes]) for z in xx]
    return standardize(a), standardize(b), {"means": centers, "sd": scales,
                                               "n50_winsor_limits": bounds}


def fit(train: list[dict], features: list[tuple], assignment: dict,
        penalty: float = LAMBDA, equal_species: bool = False) -> tuple[list[float], dict]:
    n_total = sum(r["n"] for r in train)
    mean_n = n_total / len(train)
    ns = [mean_n if equal_species else r["n"] for r in train]
    ks = [ns[i] * assignment.get(r["species"], r["neo_k"]) / r["n"]
          for i, r in enumerate(train)]
    total = sum(ns)
    p0 = (sum(ks) + 0.5) / (total + 1.0)
    beta = [math.log(p0 / (1.0 - p0))] + [0.0] * (len(features[0]) - 1)
    tolerance = max(1e-8, 1e-10 * total)
    m = len(beta)

    def evaluate(w, derivatives=False):
        obj = 0.5 * penalty * sum(v * v for v in w[1:])
        grad = [0.0] + [penalty * v for v in w[1:]]
        hess = [[0.0] * m for _ in range(m)]
        for j in range(1, m):
            hess[j][j] = penalty
        for x, n, k in zip(features, ns, ks):
            eta = sum(a * b for a, b in zip(x, w))
            obj += n * softplus(eta) - k * eta
            if derivatives:
                p = sigmoid(eta)
                residual = n * p - k
                curvature = n * p * (1.0 - p)
                for j in range(m):
                    grad[j] += residual * x[j]
                    for h in range(j + 1):
                        hess[j][h] += curvature * x[j] * x[h]
        if derivatives:
            for j in range(m):
                for h in range(j):
                    hess[h][j] = hess[j][h]
        return obj, grad, hess

    steps = backtracks = 0
    for iteration in range(101):
        old, grad, hess = evaluate(beta, True)
        max_gradient = max(map(abs, grad))
        if max_gradient <= tolerance:
            return beta, {"iterations": steps, "line_search_reductions": backtracks,
                          "max_gradient": max_gradient, "tolerance": tolerance,
                          "gradient_tolerance_ratio": max_gradient / tolerance,
                          "effective_n": total, "objective": old, "converged": True}
        if iteration == 100:
            break
        direction = cholesky_solve(hess, grad)
        descent = sum(a * b for a, b in zip(grad, direction))
        if not math.isfinite(descent) or descent <= 0:
            raise AnalysisError("non-descent Newton step")
        for attempt in range(45):
            rate = 0.5 ** attempt
            proposal = [b - rate * d for b, d in zip(beta, direction)]
            current = evaluate(proposal)[0]
            allowance = 1e-12 * max(1.0, abs(old))
            if math.isfinite(current) and current <= old - 1e-4 * rate * descent + allowance:
                beta = proposal
                steps += 1
                backtracks += attempt
                break
        else:
            raise AnalysisError("line search failed")
    raise AnalysisError("Newton fit failed declared convergence tolerance")


def cross_entropy(k: int, n: int, p: float) -> tuple[float, bool]:
    clamped = min(1 - 1e-12, max(1e-12, p))
    return -((k / n) * math.log(clamped) + (1 - k / n) * math.log1p(-clamped)), clamped != p


def score_orders(rows: list[dict], assignment: dict | None = None, *,
                 penalty: float = 1.0, equal_species: bool = False,
                 winsor: bool = False, detailed: bool = False) -> dict:
    assignment = assignment or {}
    orders = sorted({r["order"] for r in rows})
    if len(orders) < 2:
        raise AnalysisError("cannot hold out from fewer than two orders")
    per_order, predictions, diagnostics = [], [], []
    losses0, losses1, weights = [], [], []
    clips = 0
    for group in orders:
        train = [r for r in rows if r["order"] != group]
        held = [r for r in rows if r["order"] == group]
        avg = [0.0, 0.0]
        fitted = []
        for extended in (False, True):
            z_train, z_held, scaling = transform(train, held, extended, winsor)
            beta, diagnostic = fit(train, z_train, assignment, penalty, equal_species)
            fitted.append((beta, z_held))
            if detailed:
                diagnostics.append({"held_out_order": group, "model": "extended" if extended else "baseline",
                                    "train_species": len(train), "test_species": len(held),
                                    "train_mean": scaling["means"], "train_sd": scaling["sd"],
                                    "n50_winsor_limits": scaling["n50_winsor_limits"],
                                    "beta": beta, **diagnostic})
        for i, row in enumerate(held):
            k = assignment.get(row["species"], row["neo_k"])
            p = [sigmoid(sum(a * b for a, b in zip(fitted[j][0], fitted[j][1][i]))) for j in (0, 1)]
            loss0, clipped0 = cross_entropy(k, row["n"], p[0])
            loss1, clipped1 = cross_entropy(k, row["n"], p[1])
            clips += int(clipped0) + int(clipped1)
            avg[0] += loss0 / len(held)
            avg[1] += loss1 / len(held)
            losses0.append(loss0)
            losses1.append(loss1)
            weights.append(row["n"])
            if detailed:
                predictions.append({"species": row["species"], "held_out_order": group,
                                    "n": row["n"], "k": k,
                                    "baseline_p": p[0], "sine_p": p[1],
                                    "baseline_loss": loss0, "sine_loss": loss1})
        per_order.append({"order": group, "species": len(held),
                          "necropsies": sum(r["n"] for r in held),
                          "baseline_loss": avg[0], "sine_loss": avg[1], "gain": avg[0] - avg[1]})
    baseline = sum(x["baseline_loss"] for x in per_order) / len(orders)
    extended = sum(x["sine_loss"] for x in per_order) / len(orders)
    nsum = sum(weights)
    result = {"species": len(rows), "orders": len(orders),
              "baseline_loss": baseline, "sine_loss": extended, "gain": baseline - extended,
              "orders_improving": sum(x["gain"] > 0 for x in per_order),
              "necropsy_weighted_baseline": sum(a * n for a, n in zip(losses0, weights)) / nsum,
              "necropsy_weighted_sine": sum(a * n for a, n in zip(losses1, weights)) / nsum,
              "probability_clips": clips}
    if detailed:
        result.update(per_order=per_order, predictions=predictions, diagnostics=diagnostics)
    return result


def complete_refit(rows: list[dict], assignment: dict, *, detailed: bool = False) -> dict:
    primary = score_orders(rows, assignment, detailed=detailed)
    groups = sorted({r["order"] for r in rows})
    if len(groups) != 11:
        raise AnalysisError("complete primary panel must have eleven orders")
    deletions = []
    for g in groups:
        scored = score_orders([r for r in rows if r["order"] != g], assignment)
        deletions.append({"deleted_order": g, "gain": scored["gain"],
                          "remaining_species": scored["species"], "remaining_orders": scored["orders"]})
    q = [len(groups) * primary["gain"] - (len(groups) - 1) * x["gain"] for x in deletions]
    center = sum(q) / len(q)
    se = math.sqrt(sum((v - center) ** 2 for v in q) / (len(q) * (len(q) - 1)))
    primary["jackknife"] = {"center": center, "se": se,
                             "lower": center - T95_10 * se,
                             "upper": center + T95_10 * se,
                             "deletions": [{**d, "pseudovalue": v} for d, v in zip(deletions, q)]}
    return primary


def all_assignments(rows: list[dict]):
    uncertain = [r for r in rows if len(r["neo_candidates"]) > 1]
    configurations = math.prod(len(r["neo_candidates"]) for r in uncertain)
    if configurations > 4096:
        raise AnalysisError("enumeration cap exceeded")
    for values in itertools.product(*(r["neo_candidates"] for r in uncertain)):
        yield {r["species"]: k for r, k in zip(uncertain, values)}


def sensitivities(rows: list[dict], assignment: dict | None = None) -> dict:
    assignment = assignment or {}
    quality = [r for r in rows if r["x"][2] >= 90 and math.expm1(r["x"][3]) >= 100000]
    qs = len({r["order"] for r in quality})
    result = {"equal_species_training": score_orders(rows, assignment, equal_species=True),
              "ridge_0_1": score_orders(rows, assignment, penalty=0.1),
              "ridge_10": score_orders(rows, assignment, penalty=10.0),
              "training_winsor_n50": score_orders(rows, assignment, winsor=True),
              "assembly_quality": {"species": len(quality), "orders": qs,
                                   "thresholds": "BUSCO >=90; Contig_N50 >=100000 bases"}}
    if qs >= 8:
        result["assembly_quality"]["score"] = score_orders(quality, assignment)
    else:
        result["assembly_quality"]["untestable"] = "fewer than 8 retained orders"
    return result


def verdict(panel: dict, sensitivities_result: dict) -> str:
    j = panel["jackknife"]
    if panel["gain"] <= 0 and j["upper"] < 0:
        return "falsified"
    required = ("equal_species_training", "ridge_0_1", "ridge_10", "training_winsor_n50")
    quality = sensitivities_result["assembly_quality"].get("score")
    if (panel["gain"] > 0 and j["lower"] > 0 and panel["orders_improving"] > 5
            and all(d["gain"] > 0 for d in j["deletions"])
            and quality is not None and quality["gain"] > 0
            and all(sensitivities_result[key]["gain"] > 0 for key in required)):
        return "consistent"
    return "inconclusive"


def run_all(*, exhaustive: bool = True) -> dict:
    rows = load_panel()
    nominal = complete_refit(rows, {}, detailed=True)
    stress = sensitivities(rows)
    if exhaustive:
        variants = []
        for assign in all_assignments(rows):
            # Reuse the already-scored full-panel nominal assignment.
            if all(assign.get(r["species"], r["neo_k"]) == r["neo_k"] for r in rows):
                s, case_stress = nominal, stress
            else:
                s, case_stress = complete_refit(rows, assign), sensitivities(rows, assign)
            variants.append({"counts": assign, "gain": s["gain"],
                             "interval": [s["jackknife"]["lower"], s["jackknife"]["upper"]],
                             "deletion_gains": {d["deleted_order"]: d["gain"] for d in s["jackknife"]["deletions"]},
                             "orders_improving": s["orders_improving"],
                             "sensitivity_gains": {k: case_stress[k]["gain"] for k in ("equal_species_training", "ridge_0_1", "ridge_10", "training_winsor_n50")},
                             "assembly_quality_gain": case_stress["assembly_quality"].get("score", {}).get("gain"),
                             "verdict": verdict(s, case_stress)})
    else:
        variants = [{"counts": {}, "gain": nominal["gain"],
                     "interval": [nominal["jackknife"]["lower"], nominal["jackknife"]["upper"]],
                     "verdict": verdict(nominal, stress)}]
    decisions = set(x["verdict"] for x in variants)
    overall = (decisions.pop() if len(decisions) == 1 else "inconclusive") if exhaustive else "not_final_nominal_only"
    if overall == "consistent":
        # Every assignment must meet all gates; verdict() applies identical stress gains.
        pass
    elif overall == "falsified":
        if any(x["gain"] > 0 or x["interval"][1] >= 0 for x in variants):
            overall = "inconclusive"
    return {"research": "SINE incremental whole-order transport beyond L1",
            "nonclinical_unreviewed": True,
            "source": "VNF1981/Non-LTR-Retrotransposons-Cancer-Project@792af204b4113238276b3ac43e511b62872dd368:Supplementary Data S1.csv",
            "input_hashes": {PANEL.relative_to(HERE).as_posix(): normalized_hash(PANEL),
                             Path(__file__).resolve().relative_to(HERE).as_posix(): normalized_hash(Path(__file__).resolve())},
            "audit": {"species": 55, "orders": 11, "necropsies": 3533,
                      "neo_nominal": 691, "L1_zero": sum(r["x"][0] == 0 for r in rows),
                      "ambiguous_primary_species": [r["species"] for r in rows if len(r["neo_candidates"]) > 1],
                      "compatible_primary_assignments": math.prod(len(r["neo_candidates"]) for r in rows),
                      "secondary_gate": "failed: 31 source Ratio precision mismatches; not scored"},
            "nominal": nominal, "sensitivities": stress,
            "candidate_assignments": variants,
            "candidate_gain_range": [min(v["gain"] for v in variants), max(v["gain"] for v in variants)],
            "candidate_upper_range": [min(v["interval"][1] for v in variants), max(v["interval"][1] for v in variants)],
            "verdict": overall}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nominal-only", action="store_true", help="Fast sanity run; not a final verdict")
    parser.add_argument("--output", type=Path, help="Write JSON report to this path")
    args = parser.parse_args(argv)
    result = run_all(exhaustive=not args.nominal_only)
    content = json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(content, encoding="utf-8", newline="\n")
    else:
        sys.stdout.write(content)


if __name__ == "__main__":
    main()
