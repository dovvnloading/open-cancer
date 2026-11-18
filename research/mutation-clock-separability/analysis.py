#!/usr/bin/env python3
"""Offline aggregate separability audit. Nonclinical; no medical advice.

All empirical inputs come from the small CSVs beside this file. The frozen
PROTOCOL.md defines every calculation. Running this module never writes files.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
CSV_NAMES = ("blood_rates.csv", "colon_rates.csv", "neuron_rates.csv", "assay_rates.csv")
BEGIN = "<!-- BEGIN COMPUTED RESULTS -->"
END = "<!-- END COMPUTED RESULTS -->"
CONSISTENT = "Separable clock consistent"
REJECTED = "Rejected under the stated measurement assumptions"
UNTESTABLE = "Untestable"
PROTOCOL_FREEZE_COMMIT = "3dea146499a321f6fed5e175b2610c9884ac89ca"


def normalized_bytes(data: bytes) -> bytes:
    """Normalize UTF-8 text line endings without changing any other content."""
    return data.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")


def input_hashes(base: Path = HERE) -> dict[str, str]:
    """Hashes use relative POSIX keys on both Windows and Linux."""
    base = Path(base)
    paths = [Path("PROTOCOL.md")] + [Path("data") / name for name in CSV_NAMES]
    return {
        path.as_posix(): hashlib.sha256(normalized_bytes((base / path).read_bytes())).hexdigest()
        for path in paths
    }


def positive(value: float, label: str = "value") -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{label} must be positive and finite")
    return value


def interval(bounds: tuple[float, float], label: str = "interval") -> tuple[float, float]:
    lower, upper = (positive(x, label) for x in bounds)
    if lower > upper:
        raise ValueError(f"{label} endpoints are reversed")
    return lower, upper


def departures(q: float) -> dict:
    """Distinct ratio, common-ratio, and per-cell multiplicative departures."""
    q = positive(q, "Q")
    log_q = math.log(q)
    return {"Q": q, "log_Q": log_q, "D": math.exp(abs(log_q)),
            "delta": math.exp(abs(log_q) / 2), "rho": math.exp(abs(log_q) / 4)}


def ratio_of_ratios(bm: float, bh: float, cm: float, ch: float) -> dict:
    """Mouse/human blood (or neuron) ratio divided by mouse/human colon."""
    bm, bh, cm, ch = (positive(x) for x in (bm, bh, cm, ch))
    rb, rc = bm / bh, cm / ch
    return {"blood_or_neuron_ratio": rb, "colon_ratio": rc, **departures(rb / rc)}


def compatibility(q_low: float, q_high: float) -> dict:
    """Sharp deterministic bounds for a positive Q box; no confidence claim."""
    q_low, q_high = interval((q_low, q_high), "Q box")
    distance = max(0.0, math.log(q_low), -math.log(q_high))
    contains = q_low <= 1 <= q_high
    return {"q_low": q_low, "q_high": q_high, "log_distance_min": distance,
            "D_min": math.exp(distance), "delta_min": math.exp(distance / 2),
            "rho_min": math.exp(distance / 4), "contains_one": contains,
            "verdict": CONSISTENT if contains else REJECTED}


def ratio_box(bm: tuple[float, float], bh: tuple[float, float],
              cm: tuple[float, float], ch: tuple[float, float]) -> dict:
    bm, bh, cm, ch = (interval(x) for x in (bm, bh, cm, ch))
    return compatibility(bm[0] * ch[0] / (bh[1] * cm[1]),
                         bm[1] * ch[1] / (bh[0] * cm[0]))


def fit_calibration(nh: float, nm: float, ch: float, cm: float, cf: float) -> dict:
    """Fit the frozen log offset; cf is the ferret intestinal proxy, not its outcome."""
    nh, nm, ch, cm, cf = (positive(x) for x in (nh, nm, ch, cm, cf))
    rh, rm = nh / ch, nm / cm
    alpha = (math.log(rh) + math.log(rm)) / 2
    residuals = {"human": math.log(rh) - alpha, "mouse": math.log(rm) - alpha}
    multiplier = math.exp(alpha)
    return {"human_ratio": rh, "mouse_ratio": rm, "alpha": alpha,
            "multiplier": multiplier, "training_log_residuals": residuals,
            "training_log_rmse": math.sqrt(sum(x * x for x in residuals.values()) / 2),
            "predicted_ferret": cf * multiplier}


def score_transfer(observed: float, predicted: float) -> dict:
    observed, predicted = positive(observed), positive(predicted)
    ratio = observed / predicted
    error = math.log(ratio)
    return {"observed": observed, "predicted": predicted,
            "observed_over_predicted": ratio, "signed_log_error": error,
            "absolute_fold_error": math.exp(abs(error))}


def load_rates(base: Path = HERE) -> dict[str, dict[str, str]]:
    """Read literal printed values, retaining blanks as missing metadata."""
    rows = {}
    for name in CSV_NAMES:
        with (Path(base) / "data" / name).open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream):
                row_id = row["id"]
                if not row_id or row_id in rows:
                    raise ValueError(f"Missing or duplicate input ID: {row_id!r}")
                positive(row["estimate"], row_id)
                rows[row_id] = row
    return rows


def number(row: dict, field: str) -> float | None:
    return float(row[field]) if row.get(field, "") != "" else None


def printed_bounds(row: dict, rounding: bool = False) -> tuple[float, float]:
    lower, upper = number(row, "ci_lower"), number(row, "ci_upper")
    if (lower is None) != (upper is None):
        raise ValueError(f"Incomplete printed confidence interval: {row['id']}")
    if lower is None:
        lower = upper = float(row["estimate"])
    if rounding:
        resolution = positive(row["printed_resolution"], "printed resolution")
        lower -= resolution / 2
        upper += resolution / 2
    return interval((lower, upper), row["id"])


def sem_bounds(row: dict) -> tuple[float, float] | None:
    """A descriptive one-SEM box; nonpositive endpoints remain unavailable."""
    point, sem = float(row["estimate"]), number(row, "sem")
    if sem is None:
        return point, point
    if not math.isfinite(sem) or sem < 0:
        raise ValueError(f"Invalid SEM: {row['id']}")
    if point - sem <= 0:
        return None
    return point - sem, point + sem


def _input_summary(row: dict) -> dict:
    summary = {field: row.get(field, "") for field in (
        "species", "source_taxon", "tissue", "unit", "original_unit", "estimator",
        "assay", "genotype", "uncertainty_status", "approximate", "role", "source",
        "source_url", "version", "locator", "access_date", "reuse_basis",
        "reuse_terms_url", "notes")}
    for field in ("estimate", "ci_lower", "ci_upper", "ci_level", "sem", "printed_resolution"):
        summary[field] = number(row, field)
        summary[field + "_printed"] = row.get(field, "")
    return summary


def analyze(base: Path = HERE) -> dict:
    """Recompute the complete frozen audit solely from CSV inputs."""
    rows = load_rates(base)
    get = lambda row_id: float(rows[row_id]["estimate"])
    blood_ids = ("mouse_hsc", "human_hsc", "mouse_colon", "human_colon")
    for row_id in (*blood_ids, "ferret_colon_mean", "ferret_colon_slope"):
        if rows[row_id]["unit"] != "SBS_per_cell_or_genome_year":
            raise ValueError(f"Unexpected primary/colon unit: {row_id}")
    neuron_ids = ("human_neuron", "mouse_neuron", "ferret_neuron")
    for row_id in neuron_ids:
        if rows[row_id]["unit"] != "SNV_per_genome_year":
            raise ValueError(f"Unexpected neuron unit: {row_id}")
    for row_id in ("mouse_hsc_per_base", "mouse_phylogeny_per_base", "mouse_duplex_per_base"):
        if rows[row_id]["unit"] not in ("SBS_per_bp_year", "SBS_per_bp_year_interpreted"):
            raise ValueError(f"Unexpected assay unit: {row_id}")

    point = ratio_of_ratios(*(get(row_id) for row_id in blood_ids))
    printed = ratio_box(*(printed_bounds(rows[row_id]) for row_id in blood_ids))
    rounding = ratio_box(*(printed_bounds(rows[row_id], rounding=True) for row_id in blood_ids))
    primary = {"point": point, "printed_interval_box": printed, "rounding_box": rounding,
               "fixed_missing_intervals": [row_id for row_id in blood_ids
                                           if number(rows[row_id], "ci_lower") is None],
               "unconditional_verdict": UNTESTABLE}

    phylo = rows["mouse_phylogeny_per_base"]
    duplex = rows["mouse_duplex_per_base"]
    phylo_bounds, duplex_bounds = printed_bounds(phylo), printed_bounds(duplex)
    assay_ratio = float(duplex["estimate"]) / float(phylo["estimate"])
    sl, su = duplex_bounds[0] / phylo_bounds[1], duplex_bounds[1] / phylo_bounds[0]
    assay = {
        "break_even_point": point["Q"],
        "break_even_box": {"low": printed["q_low"], "high": printed["q_high"]},
        "hypothetical_bands": [{"K": k, **compatibility(printed["q_low"] / k,
                                                        printed["q_high"] * k)}
                               for k in (1, 2, 4, 8)],
        "duplex_phylogenetic_ratio": assay_ratio,
        "duplex_phylogenetic_box": {"low": sl, "high": su},
        "bridged_point": departures(point["Q"] * assay_ratio),
        "bridged_box": compatibility(printed["q_low"] * sl, printed["q_high"] * su),
    }

    nh, nm, nf = (get(row_id) for row_id in neuron_ids)
    ch, cm, cf = (get(row_id) for row_id in ("human_colon", "mouse_colon", "ferret_colon_mean"))
    calibration = fit_calibration(nh, nm, ch, cm, cf)
    neuron = {
        "point": ratio_of_ratios(nm, nh, cm, ch),
        "unpropagated_neuronal_uncertainty": [row_id for row_id in neuron_ids
                                             if number(rows[row_id], "sem") is None],
        "calibration": calibration,
        "ferret_score": score_transfer(nf, calibration["predicted_ferret"]),
        "anchor_predictions": {"human": cf * nh / ch, "mouse": cf * nm / cm},
        "transfer_scope": {
            "predictor_tissue": rows["ferret_colon_mean"]["tissue"],
            "matched_colon_transfer_verdict": UNTESTABLE,
            "assumption": "For numerical proxy scores only: ferret small-intestine rate substitutes "
                          "for unmeasured ferret colon rate; protocol ID names are not tissue claims.",
        },
    }
    hb, mb, fb = (sem_bounds(rows[row_id]) for row_id in neuron_ids)
    if hb is None or mb is None:
        neuron["one_sem_box"] = {"available": False, "reason": "Nonpositive neuronal one-SEM endpoint"}
    else:
        neuron["one_sem_box"] = {"available": True, **ratio_box(mb, hb, (cm, cm), (ch, ch))}
    if hb is None or mb is None or fb is None:
        neuron["one_sem_transfer"] = {"available": False, "reason": "Nonpositive neuronal one-SEM endpoint"}
    else:
        al = (math.log(hb[0] / ch) + math.log(mb[0] / cm)) / 2
        au = (math.log(hb[1] / ch) + math.log(mb[1] / cm)) / 2
        ml, mu = math.exp(al), math.exp(au)
        pl, pu = cf * ml, cf * mu
        neuron["one_sem_transfer"] = {
            "available": True, "alpha_low": al, "alpha_high": au,
            "multiplier_low": ml, "multiplier_high": mu,
            "predicted_low": pl, "predicted_high": pu,
            "observed_low": fb[0], "observed_high": fb[1],
            "ratio_low": fb[0] / pu, "ratio_high": fb[1] / pl,
            "intervals_intersect": max(pl, fb[0]) <= min(pu, fb[1]),
        }
    return {"scope": "Nonclinical, unreviewed computational research; no medical advice",
            "protocol_freeze_commit": PROTOCOL_FREEZE_COMMIT,
            "inputs": {row_id: _input_summary(row) for row_id, row in rows.items()},
            "primary": primary, "assay": assay, "neuron": neuron,
            "hashes": input_hashes(base)}


def _fmt(value: float) -> str:
    return format(value, ".6g")


def _range(lower: float, upper: float) -> str:
    return f"[{_fmt(lower)}, {_fmt(upper)}]"


def _box_row(label: str, box: dict) -> str:
    return (f"| {label} | {_range(box['q_low'], box['q_high'])} | {_fmt(box['D_min'])} | "
            f"{_fmt(box['delta_min'])} | {_fmt(box['rho_min'])} | {box['verdict']} |")


def render_results(result: dict) -> str:
    """Generate all numerical result tables, including source-value auditing."""
    p, a, n = (result[key] for key in ("primary", "assay", "neuron"))
    lines = [BEGIN, "", "## Recomputed aggregate results", "",
             "Nonclinical and unreviewed. Values below are generated from the CSVs; "
             "display precision is not measurement precision.", "",
             f"Original local protocol freeze commit: `{result['protocol_freeze_commit']}`.", "",
             "### Transcribed numerical inputs", "",
             "`SBS_per_cell_or_genome_year` preserves the source cell/genome distinction; "
             "the primary comparison assumes these genome-wide counts are commensurate. "
             "The `_interpreted` per-base units preserve the separate unit caveat.", "",
             "| Input ID | Source tissue | Printed annual rate | Unit | Printed uncertainty | Role |",
             "|---|---|---:|---|---|---|"]
    for row_id, row in result["inputs"].items():
        uncertainty = "Not printed"
        if row["ci_lower"] is not None:
            level = f"{_fmt(row['ci_level'])}% " if row["ci_level"] is not None else ""
            uncertainty = f"{level}CI [{row['ci_lower_printed']}, {row['ci_upper_printed']}]"
        if row["sem"] is not None:
            uncertainty = f"SEM {row['sem_printed']}"
        rate = row["estimate_printed"]
        if row["approximate"].lower() in ("yes", "true", "1"):
            rate = "approximately " + rate
        lines.append(f"| `{row_id}` | {row['tissue']} | {rate} | `{row['unit']}` | {uncertainty} | {row['role']} |")
    lines += ["", "| Input ID | Source unit wording | Retained uncertainty status |",
              "|---|---|---|"]
    for row_id, row in result["inputs"].items():
        lines.append(f"| `{row_id}` | {row['original_unit']} | {row['uncertainty_status']} |")
    lines += ["", "### Peer-reviewed blood/colon comparison", "",
              "| Statistic | Value |", "|---|---:|",
              f"| Mouse/human blood ratio | {_fmt(p['point']['blood_or_neuron_ratio'])} |",
              f"| Mouse/human colon ratio | {_fmt(p['point']['colon_ratio'])} |",
              f"| Q: blood ratio / colon ratio | {_fmt(p['point']['Q'])} |",
              f"| Signed log interaction | {_fmt(p['point']['log_Q'])} |",
              f"| D: fold disagreement | {_fmt(p['point']['D'])} |",
              f"| delta: minimum symmetric allowance per tissue-specific species ratio | {_fmt(p['point']['delta'])} |",
              f"| rho: minimum symmetric allowance per species-tissue cell | {_fmt(p['point']['rho'])} |", "",
              "| Box | Q range | D minimum | delta minimum | rho minimum | Conditional verdict |",
              "|---|---|---:|---:|---:|---|",
              _box_row("Printed blood intervals; colon points fixed", p["printed_interval_box"]),
              _box_row("Outward rounding sensitivity", p["rounding_box"]), "",
              "The printed-interval result is conditional compatibility, not a simultaneous confidence "
              "test or a p-value. Missing colon uncertainty is unpropagated. The unconditional "
              f"biological verdict with unrestricted measurement interactions is **{p['unconditional_verdict']}**.", "",
              "### Assay/normalisation sensitivity", "",
              f"The latent separable null requires measurement interaction A = {_fmt(a['break_even_point'])} "
              f"at the point estimates; its conditional compatible A range is "
              f"{_range(a['break_even_box']['low'], a['break_even_box']['high'])}.", "",
              "| Hypothetical interaction band | Latent Q range | D minimum | delta minimum | rho minimum | Conditional verdict |",
              "|---|---|---:|---:|---:|---|"]
    for box in a["hypothetical_bands"]:
        k = box["K"]
        lines.append(_box_row(f"A in [1/{k}, {k}]", box))
    lines += ["", "These bands are hypothetical stress scenarios. A rescued band establishes compatibility "
              "within the allowance, not biological confirmation.", "",
              f"The targeted-duplex/phylogenetic per-base rate ratio s is "
              f"**{_fmt(a['duplex_phylogenetic_ratio'])}**, with marginal endpoint ratio box "
              f"**{_range(a['duplex_phylogenetic_box']['low'], a['duplex_phylogenetic_box']['high'])}**. "
              "Applying it to the primary mouse blood rate is an additional assay bridge, not a matched "
              "whole-blood annual observation or independent replication.", "",
              "| Bridged point statistic | Value |", "|---|---:|"]
    for key in ("Q", "log_Q", "D", "delta", "rho"):
        lines.append(f"| {key} | {_fmt(a['bridged_point'][key])} |")
    lines += ["", "| Bridged box | Q range | D minimum | delta minimum | rho minimum | Conditional verdict |",
              "|---|---|---:|---:|---:|---|", _box_row("Q multiplied by s endpoint box", a["bridged_box"]), "",
              "### Preprint neuron transport negative control", "",
              "Postmitotic neurons are a transport negative control, not cancer targets. "
              "These results remain separate from the peer-reviewed blood/colon comparison.", "",
              "| Human/mouse neuronal comparison | Value |", "|---|---:|",
              f"| Mouse/human neuron ratio | {_fmt(n['point']['blood_or_neuron_ratio'])} |",
              f"| Mouse/human colon ratio | {_fmt(n['point']['colon_ratio'])} |"]
    for key in ("Q", "log_Q", "D", "delta", "rho"):
        lines.append(f"| {key} | {_fmt(n['point'][key])} |")
    lines.append("")
    if n["one_sem_box"]["available"]:
        lines += ["| Descriptive box | Q range | D minimum | delta minimum | rho minimum | Box verdict |",
                  "|---|---|---:|---:|---:|---|", _box_row("Plus/minus one SEM; missing uncertainties fixed", n["one_sem_box"]), ""]
    else:
        lines += [f"One-SEM log comparison unavailable: {n['one_sem_box']['reason']}.", ""]
    c, f = n["calibration"], n["ferret_score"]
    lines += ["| Frozen human/mouse calibration | Value |", "|---|---:|",
              f"| Human neuron/colon ratio | {_fmt(c['human_ratio'])} |",
              f"| Mouse neuron/colon ratio | {_fmt(c['mouse_ratio'])} |",
              f"| Equal-weight log offset alpha | {_fmt(c['alpha'])} |",
              f"| Multiplicative offset exp(alpha) | {_fmt(c['multiplier'])} |",
              f"| Human signed training log residual | {_fmt(c['training_log_residuals']['human'])} |",
              f"| Mouse signed training log residual | {_fmt(c['training_log_residuals']['mouse'])} |",
              f"| Training log RMSE | {_fmt(c['training_log_rmse'])} |", "",
              "### Ferret proxy transfer (small-intestine predictor)", "",
              f"Matched-colon ferret transfer verdict: **{n['transfer_scope']['matched_colon_transfer_verdict']}**. "
              "Postfreeze source review identifies the ferret predictor as small intestine, while "
              "the human/mouse training predictors are colon. The frozen numerical rule is retained "
              "as a proxy score conditional on an unverified ferret small-intestine=colon bridge. "
              "The original `ferret_colon_*` row IDs preserve protocol traceability; they do not "
              "assert the source tissue is colon.", "",
              "| Ferret proxy holdout and frozen anchor diagnostics | Value |", "|---|---:|",
              f"| Predicted annual neuron SNV/genome rate | {_fmt(f['predicted'])} |",
              f"| Observed annual neuron SNV/genome rate | {_fmt(f['observed'])} |",
              f"| Observed/predicted | {_fmt(f['observed_over_predicted'])} |",
              f"| Signed natural-log error | {_fmt(f['signed_log_error'])} |",
              f"| Absolute fold error | {_fmt(f['absolute_fold_error'])} |",
              f"| Human-only anchor prediction, SNV/genome/year | {_fmt(n['anchor_predictions']['human'])} |",
              f"| Mouse-only anchor prediction, SNV/genome/year | {_fmt(n['anchor_predictions']['mouse'])} |", ""]
    s = n["one_sem_transfer"]
    if s["available"]:
        lines += ["| Plus/minus one-SEM transfer sensitivity | Range or result |", "|---|---|",
                  f"| Log offset alpha | {_range(s['alpha_low'], s['alpha_high'])} |",
                  f"| Multiplicative offset | {_range(s['multiplier_low'], s['multiplier_high'])} |",
                  f"| Predicted ferret SNV/genome/year | {_range(s['predicted_low'], s['predicted_high'])} |",
                  f"| Observed ferret SNV/genome/year | {_range(s['observed_low'], s['observed_high'])} |",
                  f"| Observed/predicted ratio | {_range(s['ratio_low'], s['ratio_high'])} |",
                  f"| Prediction and observation ranges intersect | {'Yes' if s['intervals_intersect'] else 'No'} |", ""]
    else:
        lines += [f"One-SEM transfer log bounds unavailable: {s['reason']}.", ""]
    missing = ", ".join(f"`{row_id}`" for row_id in n["unpropagated_neuronal_uncertainty"])
    lines += ["SEMs were not converted to confidence limits. Neuronal uncertainty left unpropagated: "
              f"{missing or 'none'}. Intestinal predictor uncertainty is unpropagated throughout the neuron analysis. "
              "The ferret small-intestine mean is used as a colon proxy; the separate regression slope is audit-only. "
              "The ferret outcome does not enter calibration or model selection.", "",
              "### Input integrity", "",
              "SHA-256 hashes use newline-normalized UTF-8 content and POSIX relative paths. "
              "They establish input integrity, not source truth or biological validity.", "",
              "| Input | SHA-256 |", "|---|---|"]
    for path, digest in result["hashes"].items():
        lines.append(f"| `{path}` | `{digest}` |")
    lines += ["", END]
    return "\n".join(lines)


def check_readme(base: Path = HERE) -> bool:
    text = (Path(base) / "README.md").read_text(encoding="utf-8")
    if text.count(BEGIN) != 1 or text.count(END) != 1:
        return False
    start, stop = text.index(BEGIN), text.index(END)
    if stop < start:
        return False
    return text[start:stop + len(END)] == render_results(analyze(base))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--markdown", action="store_true", help="print the generated numerical README block")
    group.add_argument("--check", action="store_true", help="verify the numerical README block without writing")
    args = parser.parse_args()
    if args.check:
        ok = check_readme()
        print("Computed README block matches CSV inputs." if ok else "Computed README block differs or is missing.")
        return 0 if ok else 1
    result = analyze()
    print(render_results(result) if args.markdown else json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
