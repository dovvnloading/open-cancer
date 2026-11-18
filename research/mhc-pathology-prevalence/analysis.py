#!/usr/bin/env python3
"""Nonclinical MHC-I pathology-prevalence audit; Python standard library only.

No network is used here. Pinned retrieval and derivation live separately in
``data/fetch_verify.py``. This evaluates endpoint portability of a measured
proxy; it supplies no medical or veterinary advice or immune-causation claim.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import statistics


BASE = Path(__file__).resolve().parent
CONTINUOUS = ("log_body_mass_kg", "log_adult_life_expectancy_days")
MHC = "log1p_mhc_i"
DIET = "mammal_prey"
N = "records"
K = "neoplasia"
SECONDARY_K = "malignant"
T_8 = 2.306004135204166
P_CLIP = 1e-12
MAX_ITERATIONS = 100
RELATIVE_GRADIENT_TOLERANCE = 1e-10
SNAPSHOT_TOLERANCE = 3e-7
STATUS = "nonclinical_unreviewed_research_no_medical_or_veterinary_advice"
NAMED_RODENTS = frozenset(("Heterocephalus_glaber", "Dolichotis_patagonum"))
MAJOR_SENSITIVITIES = (
    "equal_species_training", "penalty_0_1", "penalty_10", "exclude_named_rodents",
)
DIAGNOSTIC_KEYS = frozenset((
    "maximum_gradient_ratio", "maximum_absolute_gradient", "iterations",
    "maximum_fit_iterations", "line_search_halvings",
))
BEGIN = "<!-- BEGIN GENERATED RESULTS -->"
END = "<!-- END GENERATED RESULTS -->"


class FitError(RuntimeError):
    """A failed numerical fit must stop the analysis, never disappear silently."""


def normalize_newlines(text):
    return text.replace("\r\n", "\n").replace("\r", "\n")


def normalized_sha256(path):
    canonical = normalize_newlines(Path(path).read_bytes().decode("utf-8"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def input_hashes(base=BASE):
    base = Path(base)
    paths = [Path("PROTOCOL.md"), Path("analysis.py")]
    paths += [Path("data") / name for name in (
        "species_aggregates.csv", "reconciliation.csv", "endpoint_audit.csv",
        "source_audit.json", "fetch_verify.py",
    )]
    return {path.as_posix(): normalized_sha256(base / path) for path in paths}


def _read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def load_rows(base=BASE):
    """Read only the frozen derived panel; do not alter its membership."""
    rows = _read_csv(Path(base) / "data" / "species_aggregates.csv")
    required = {"species", "order", N, K, SECONDARY_K, DIET, MHC, *CONTINUOUS}
    for row in rows:
        if not required <= set(row):
            raise ValueError("Missing frozen aggregate columns")
        for key in (N, K, SECONDARY_K, DIET):
            row[key] = int(row[key])
        for key in (*CONTINUOUS, MHC):
            row[key] = float(row[key])
            if not math.isfinite(row[key]):
                raise ValueError("Nonfinite predictor: " + row["species"])
        if row[N] <= 0 or not 0 <= row[K] <= row[N]:
            raise ValueError("Invalid denominator-linked neoplasia count")
        if not 0 <= row[SECONDARY_K] <= row[K]:
            raise ValueError("Invalid recorded malignant count")
        if row[DIET] not in (0, 1) or row[MHC] < 0:
            raise ValueError("Invalid diet or MHC-I predictor")
        if not row["species"] or not row["order"]:
            raise ValueError("Empty species/order identifier")
    names = [row["species"] for row in rows]
    if names != sorted(set(names)):
        raise ValueError("Species identifiers must be unique and sorted")
    return rows


def sigmoid(value):
    if value >= 0:
        inverse = math.exp(-value)
        return 1.0 / (1.0 + inverse)
    exponent = math.exp(value)
    return exponent / (1.0 + exponent)


def softplus(value):
    return max(value, 0.0) + math.log1p(math.exp(-abs(value)))


def pathology_loss(k, n, prediction):
    """Per-record cross-entropy, including the k=0 case without pseudocounts."""
    if n <= 0 or not 0 <= k <= n:
        raise ValueError("Invalid pathology count or denominator")
    if not math.isfinite(prediction) or not 0 <= prediction <= 1:
        raise ValueError("Prediction must be a finite probability")
    p = min(1.0 - P_CLIP, max(P_CLIP, prediction))
    return -(k * math.log(p) + (n - k) * math.log1p(-p)) / n


def _cholesky_solve(matrix, vector):
    """Solve a small positive-definite Newton system without dependencies."""
    size = len(vector)
    lower = [[0.0] * size for _ in range(size)]
    for i in range(size):
        for j in range(i + 1):
            residual = matrix[i][j] - sum(
                lower[i][k] * lower[j][k] for k in range(j))
            if i == j:
                if residual <= 0 or not math.isfinite(residual):
                    raise FitError("Newton Hessian is not positive definite")
                lower[i][j] = math.sqrt(residual)
            else:
                lower[i][j] = residual / lower[j][j]
    intermediate = []
    for i in range(size):
        intermediate.append((vector[i] - sum(
            lower[i][j] * intermediate[j] for j in range(i))) / lower[i][i])
    answer = [0.0] * size
    for i in range(size - 1, -1, -1):
        answer[i] = (intermediate[i] - sum(
            lower[j][i] * answer[j] for j in range(i + 1, size))) / lower[i][i]
    return answer


def _feature_vector(row, means, scales, include_mhc):
    vector = [1.0]
    vector += [(row[key] - means[key]) / scales[key] for key in CONTINUOUS]
    vector.append(float(row[DIET]))
    if include_mhc:
        vector.append((row[MHC] - means[MHC]) / scales[MHC])
    return vector


def fit_model(train, include_mhc, penalty=1.0, weighting="counts", outcome=K):
    """Fit the frozen ridge-binomial rule using training species only.

    `species` equalizes effective denominators to the training mean n.
    The score always uses the original denominators and unaltered test labels.
    """
    if not train or penalty <= 0 or not math.isfinite(penalty):
        raise ValueError("Training rows and a positive finite ridge are required")
    if weighting not in ("counts", "species"):
        raise ValueError("Unknown training weighting")
    continuous = CONTINUOUS + ((MHC,) if include_mhc else ())
    means = {key: statistics.fmean(row[key] for row in train) for key in continuous}
    scales = {key: math.sqrt(statistics.fmean(
        (row[key] - means[key]) ** 2 for row in train)) or 1.0
        for key in continuous}
    design = [_feature_vector(row, means, scales, include_mhc) for row in train]
    mean_n = statistics.fmean(row[N] for row in train)
    denominators = []
    cancer = []
    for row in train:
        multiplier = 1.0
        if weighting == "species":
            multiplier = mean_n / row[N]
        denominators.append(row[N] * multiplier)
        cancer.append(row[outcome] * multiplier)
    total_n = math.fsum(denominators)
    total_k = math.fsum(cancer)
    if not 0 < total_k < total_n:
        raise FitError("No finite unpenalized-intercept optimum: constant outcome")
    dimension = len(design[0])
    coefficients = [math.log((total_k + 0.5) / (total_n - total_k + 0.5))]
    coefficients += [0.0] * (dimension - 1)
    gradient_tolerance = max(1e-8, RELATIVE_GRADIENT_TOLERANCE * total_n)

    def objective(beta):
        likelihood = 0.0
        for features, n, k in zip(design, denominators, cancer):
            eta = sum(a * b for a, b in zip(features, beta))
            likelihood += n * softplus(eta) - k * eta
        return likelihood + 0.5 * penalty * sum(v * v for v in beta[1:])

    objective_value = objective(coefficients)
    accepted_steps = 0
    total_halvings = 0
    for iteration in range(MAX_ITERATIONS + 1):
        gradient = [0.0] * dimension
        hessian = [[0.0] * dimension for _ in range(dimension)]
        for features, n, k in zip(design, denominators, cancer):
            eta = sum(a * b for a, b in zip(features, coefficients))
            probability = sigmoid(eta)
            residual = n * probability - k
            curvature = n * probability * (1.0 - probability)
            for a in range(dimension):
                xa = features[a]
                gradient[a] += residual * xa
                curvature_xa = curvature * xa
                for b in range(a + 1):
                    hessian[a][b] += curvature_xa * features[b]
        for a in range(1, dimension):
            gradient[a] += penalty * coefficients[a]
            hessian[a][a] += penalty
        for a in range(dimension):
            for b in range(a):
                hessian[b][a] = hessian[a][b]
        maximum_gradient = max(abs(value) for value in gradient)
        if maximum_gradient <= gradient_tolerance:
            return {
                "include_mhc": bool(include_mhc),
                "outcome": outcome,
                "feature_names": ["intercept", *CONTINUOUS, DIET] +
                                 ([MHC] if include_mhc else []),
                "coefficients": coefficients,
                "means": means,
                "scales": scales,
                "penalty": penalty,
                "weighting": weighting,
                "training_species": len(train),
                "training_records": sum(row[N] for row in train),
                "effective_training_denominator": total_n,
                "converged": True,
                "iterations": accepted_steps,
                "line_search_halvings": total_halvings,
                "objective": objective_value,
                "maximum_absolute_gradient": maximum_gradient,
                "gradient_tolerance": gradient_tolerance,
            }
        if iteration == MAX_ITERATIONS:
            break
        direction = _cholesky_solve(hessian, gradient)
        decrement = sum(g * d for g, d in zip(gradient, direction))
        if not math.isfinite(decrement) or decrement <= 0:
            raise FitError("Invalid Newton descent direction")
        fraction = 1.0
        for halving in range(45):
            candidate = [b - fraction * d for b, d in zip(coefficients, direction)]
            candidate_objective = objective(candidate)
            # The slack is numerical objective roundoff, not a convergence rule.
            roundoff = 1e-12 * max(1.0, abs(objective_value))
            if candidate_objective <= objective_value - 1e-4 * fraction * decrement + roundoff:
                coefficients = candidate
                objective_value = candidate_objective
                accepted_steps += 1
                total_halvings += halving
                break
            fraction *= 0.5
        else:
            raise FitError("Newton line search failed")
    raise FitError("Newton fit failed the gradient convergence criterion")


def predict_one(row, fitted):
    features = _feature_vector(row, fitted["means"], fitted["scales"],
                               fitted["include_mhc"])
    return sigmoid(sum(a * b for a, b in zip(features, fitted["coefficients"])))


def evaluate(rows, penalty=1.0, weighting="counts", outcome=K):
    """Refit the two rules with each complete taxonomic order held out."""
    if outcome not in (K, SECONDARY_K):
        raise ValueError("Only the frozen primary and secondary outcomes are allowed")
    groups = sorted({row["order"] for row in rows})
    if len(groups) < 2:
        raise FitError("At least two complete orders are required")
    predictions, folds, group_scores = [], [], []
    clipping = 0
    for group in groups:
        train = [row for row in rows if row["order"] != group]
        test = [row for row in rows if row["order"] == group]
        baseline = fit_model(train, False, penalty, weighting, outcome)
        extended = fit_model(train, True, penalty, weighting, outcome)
        fold_predictions = []
        for row in test:
            probabilities = (predict_one(row, baseline), predict_one(row, extended))
            clipping += sum(p < P_CLIP or p > 1 - P_CLIP for p in probabilities)
            losses = [pathology_loss(row[outcome], row[N], p) for p in probabilities]
            record = {
                "species": row["species"], "order": group, N: row[N],
                "events": row[outcome], "baseline_probability": probabilities[0],
                "extended_probability": probabilities[1], "baseline_loss": losses[0],
                "extended_loss": losses[1], "delta": losses[0] - losses[1],
            }
            predictions.append(record)
            fold_predictions.append(record)
        means = [statistics.fmean(p[key] for p in fold_predictions)
                 for key in ("baseline_loss", "extended_loss")]
        group_scores.append({
            "group": group, "species": len(test), "records": sum(r[N] for r in test),
            "events": sum(r[outcome] for r in test),
            "zero_case_species": sum(r[outcome] == 0 for r in test),
            "baseline_loss": means[0], "extended_loss": means[1],
            "delta": means[0] - means[1],
        })
        folds.append({
            "heldout_group": group,
            "training_species_ids": [r["species"] for r in train],
            "test_species_ids": [r["species"] for r in test],
            "training_orders": sorted({r["order"] for r in train}),
            "test_orders": [group], "baseline_fit": baseline, "extended_fit": extended,
        })
    total_n = sum(r[N] for r in rows)
    means = [statistics.fmean(g[key] for g in group_scores)
             for key in ("baseline_loss", "extended_loss")]
    record_means = [math.fsum(p[N] * p[key] for p in predictions) / total_n
                    for key in ("baseline_loss", "extended_loss")]
    fits = [fold[key] for fold in folds for key in ("baseline_fit", "extended_fit")]
    return {
        "outcome": outcome, "weighting": weighting, "penalty": penalty,
        "species": len(rows), "orders": len(groups), "groups": len(groups),
        "records": total_n, "events": sum(r[outcome] for r in rows),
        "zero_case_species": sum(r[outcome] == 0 for r in rows),
        "baseline_loss": means[0], "extended_loss": means[1],
        "delta": means[0] - means[1],
        "record_weighted_baseline_loss": record_means[0],
        "record_weighted_extended_loss": record_means[1],
        "record_weighted_delta": record_means[0] - record_means[1],
        "improving_groups": sum(g["delta"] > 0 for g in group_scores),
        "probabilities_clipped_for_scoring": clipping,
        "fit_count": len(fits), "all_fits_converged": all(f["converged"] for f in fits),
        "maximum_fit_iterations": max(f["iterations"] for f in fits),
        "maximum_gradient_ratio": max(f["maximum_absolute_gradient"] /
                                      f["gradient_tolerance"] for f in fits),
        "group_scores": group_scores,
        "predictions": sorted(predictions, key=lambda item: item["species"]),
        "folds": folds,
    }


def _summary(evaluation):
    return {key: value for key, value in evaluation.items()
            if key not in ("predictions", "folds")}


def jackknife_uncertainty(rows, primary, outcome=K):
    """Delete each order and refit every remaining holdout, including scaling."""
    groups = sorted({row["order"] for row in rows})
    count = len(groups)
    if count != 9:
        raise ValueError("The frozen t8 uncertainty requires the nine-order common panel")
    if primary["outcome"] != outcome:
        raise ValueError("Jackknife endpoint and primary estimate disagree")
    deletions = []
    for group in groups:
        retained = [r for r in rows if r["order"] != group]
        result = evaluate(retained, outcome=outcome)
        deletions.append({"deleted_order": group, **_summary(result)})
    pseudovalues = [count * primary["delta"] - (count - 1) * item["delta"]
                    for item in deletions]
    center = statistics.fmean(pseudovalues)
    standard_error = math.sqrt(math.fsum((q - center) ** 2 for q in pseudovalues) /
                               (count * (count - 1)))
    for item, q in zip(deletions, pseudovalues):
        item["pseudovalue"] = q
    return {
        "method": "complete_delete_one_order_refit_jackknife", "outcome": outcome,
        "groups": count, "degrees_of_freedom": count - 1, "t_critical": T_8,
        "full_panel_point_estimate": primary["delta"], "bias_corrected_center": center,
        "standard_error": standard_error,
        "ci_lower": center - T_8 * standard_error,
        "ci_upper": center + T_8 * standard_error,
        "delete_order_delta_min": min(d["delta"] for d in deletions),
        "delete_order_delta_max": max(d["delta"] for d in deletions),
        "deletions": deletions,
    }


def verdict(primary, uncertainty, sensitivities, panel_valid=True):
    """Apply frozen primary-only directional gates; secondary cannot rescue it."""
    groups = [g["group"] for g in primary.get("group_scores", [])]
    deleted = [d["deleted_order"] for d in uncertainty.get("deletions", [])]
    valid = bool(
        panel_valid and primary["orders"] >= 8 and primary["all_fits_converged"]
        and len(groups) == primary["orders"] and len(set(groups)) == len(groups)
        and sorted(deleted) == sorted(groups)
        and all(d["all_fits_converged"] for d in uncertainty["deletions"])
    )
    major_present = all(key in sensitivities for key in MAJOR_SENSITIVITIES)
    major_valid = major_present and all(
        sensitivities[key]["all_fits_converged"] for key in MAJOR_SENSITIVITIES)
    major_positive = major_valid and all(
        sensitivities[key]["delta"] > 0 for key in MAJOR_SENSITIVITIES)
    gates = {
        "valid_panel_and_primary_uncertainty_fits": valid,
        "positive_primary_delta": primary["delta"] > 0,
        "positive_jackknife_lower_bound": uncertainty["ci_lower"] > 0,
        "strict_majority_orders_improve": primary["improving_groups"] > primary["orders"] / 2,
        "all_four_major_sensitivities_valid": major_valid,
        "all_four_major_sensitivities_positive": major_positive,
    }
    nonimprovement = valid and primary["delta"] <= 0 and uncertainty["ci_upper"] < 0
    if not valid:
        code, label = "inconclusive", "Inconclusive (not evaluable: panel or fit validity failure)"
    elif all(gates.values()):
        code, label = "transports", "MHC-I transports to pathology prevalence within this evaluation"
    elif nonimprovement:
        code, label = "does_not_transport", "MHC-I does not transport to pathology prevalence within this evaluation"
    else:
        code, label = "inconclusive", "Inconclusive on the available pathology-prevalence panel"
    return {
        "code": code, "label": label, "positive_claim_gates": gates,
        "nonimprovement_gate": nonimprovement,
        "point_improvement_failure": primary["delta"] <= 0,
        "interval_contains_zero": uncertainty["ci_lower"] <= 0 <= uncertainty["ci_upper"],
        "secondary_can_rescue_primary": False,
    }


def _panel_summary(rows):
    names = {r["species"] for r in rows}
    return {
        "species": len(rows), "orders": len({r["order"] for r in rows}),
        "records": sum(r[N] for r in rows), "neoplasia": sum(r[K] for r in rows),
        "malignant": sum(r[SECONDARY_K] for r in rows),
        "zero_neoplasia_species": sum(r[K] == 0 for r in rows),
        "zero_malignant_species": sum(r[SECONDARY_K] == 0 for r in rows),
        "species_per_order": {g: sum(r["order"] == g for r in rows)
                              for g in sorted({r["order"] for r in rows})},
        "named_rodents_present": sorted(names & NAMED_RODENTS),
        "named_rodents_absent": sorted(NAMED_RODENTS - names),
    }


def prefit_summary(rows, base=BASE):
    """Recount the source union and endpoint audit from committed derived tables."""
    base = Path(base)
    reconciliation = _read_csv(base / "data" / "reconciliation.csv")
    endpoint_audit = _read_csv(base / "data" / "endpoint_audit.csv")
    normalization = lambda name: " ".join(name.replace("_", " ").split()).casefold()
    panel_keys = {normalization(row["species"]) for row in rows}
    ledger_keys = [r["normalized_species"] for r in reconciliation]
    audit_keys = [r["normalized_species"] for r in endpoint_audit]
    if len(set(ledger_keys)) != len(ledger_keys) or len(set(audit_keys)) != len(audit_keys):
        raise ValueError("Duplicate normalized identifier in derived audits")
    vincze = {r["normalized_species"] for r in reconciliation if r["vincze_species"]}
    compton = {r["normalized_species"] for r in reconciliation if r["compton_species"]}
    included = [r for r in reconciliation if r["status"] == "included"]
    if panel_keys != {r["normalized_species"] for r in included} or set(audit_keys) != compton:
        raise ValueError("Reconciliation or endpoint audit does not cover its stated panel")
    if set(ledger_keys) != vincze | compton:
        raise ValueError("Reconciliation does not equal source identifier union")
    counts = {"neoplasia_prevalence_mismatches": 0, "malignant_prevalence_mismatches": 0,
              "malignancyknown_prevalence_mismatches": 0,
              "panel_malignancyknown_prevalence_mismatches": 0}
    for row in endpoint_audit:
        neo_tolerance = Decimal(row["neoplasia_tolerance"])
        malignant_tolerance = Decimal(row["malignancy_tolerance"])
        quantities = [Decimal(row[key]) for key in (
            "neoplasia_residual", "malignant_residual", "malignancyknown_residual")]
        if (not all(q.is_finite() for q in (*quantities, neo_tolerance, malignant_tolerance))
                or neo_tolerance <= 0 or malignant_tolerance <= 0):
            raise ValueError("Invalid endpoint arithmetic audit")
        neo_bad = abs(quantities[0]) > neo_tolerance
        malignant_bad = abs(quantities[1]) > malignant_tolerance
        known_bad = abs(quantities[2]) > malignant_tolerance
        counts["neoplasia_prevalence_mismatches"] += neo_bad
        counts["malignant_prevalence_mismatches"] += malignant_bad
        counts["malignancyknown_prevalence_mismatches"] += known_bad
        counts["panel_malignancyknown_prevalence_mismatches"] += known_bad and row["normalized_species"] in panel_keys
    panel = _panel_summary(rows)
    return {
        "vincze_species": len(vincze), "compton_species": len(compton),
        "normalized_matches": len(vincze & compton), "union_species": len(vincze | compton),
        "primary_species": len(rows), "orders": panel["orders"],
        "species_per_order": panel["species_per_order"],
        **{key: panel[key] for key in ("records", "neoplasia", "malignant",
                                      "zero_neoplasia_species", "zero_malignant_species")},
        "matched_exclusions": [r["vincze_species"] for r in reconciliation
                               if r["normalized_species"] in vincze & compton
                               and r["status"] != "included"],
        "primary_classes": sorted({r["compton_class"] for r in included}),
        "order_label_disagreements": [
            {"species": r["vincze_species"], "vincze_order": r["vincze_order"],
             "compton_order": r["compton_order"]}
            for r in included if r["vincze_order"] != r["compton_order"]],
        "endpoint_audit": {"audited_compton_species": len(endpoint_audit), **counts},
    }


def build_report(base=BASE):
    """Recompute all fits, scores, uncertainties and the verdict offline."""
    base = Path(base)
    rows = load_rows(base)
    source_audit = json.loads((base / "data" / "source_audit.json").read_text(encoding="utf-8"))
    for relative, expected in source_audit.get("derived_file_hashes_lf", {}).items():
        if Path(relative).as_posix() != relative or "\\" in relative:
            raise ValueError("Source audit has a nonportable path key")
        if normalized_sha256(base / "data" / relative) != expected:
            raise ValueError("Derived file differs from source audit: " + relative)
    prefit = prefit_summary(rows, base)
    if prefit != source_audit["prefit_summary"]:
        raise ValueError("Recomputed prefit source audit differs from committed record")
    panel = _panel_summary(rows)
    for key, value in (("species", 42), ("orders", 9), ("records", 2215),
                       ("neoplasia", 312), ("zero_neoplasia_species", 2)):
        if panel[key] != value:
            raise ValueError("Frozen panel gate failed: " + key)
    primary = evaluate(rows)
    uncertainty = jackknife_uncertainty(rows, primary)
    secondary = evaluate(rows, outcome=SECONDARY_K)
    secondary_uncertainty = jackknife_uncertainty(rows, secondary, outcome=SECONDARY_K)
    sensitivities = {
        "equal_species_training": _summary(evaluate(rows, weighting="species")),
        "penalty_0_1": _summary(evaluate(rows, penalty=0.1)),
        "penalty_10": _summary(evaluate(rows, penalty=10.0)),
        "exclude_named_rodents": _summary(evaluate(
            [r for r in rows if r["species"] not in NAMED_RODENTS])),
        "exclude_recorded_zeros": _summary(evaluate([r for r in rows if r[K] > 0])),
    }
    evaluations = [primary, *uncertainty["deletions"], secondary,
                   *secondary_uncertainty["deletions"], *sensitivities.values()]
    return {
        "schema_version": 1, "status": STATUS, "protocol": Path("PROTOCOL.md").as_posix(),
        "analysis_settings": {
            "primary_penalty": 1.0, "score_probability_clip": P_CLIP,
            "relative_gradient_tolerance": RELATIVE_GRADIENT_TOLERANCE,
            "minimum_gradient_tolerance": 1e-8, "maximum_newton_iterations": MAX_ITERATIONS,
            "maximum_halvings_per_step": 44, "uncertainty_t_critical": T_8,
            "snapshot_absolute_and_relative_tolerance": SNAPSHOT_TOLERANCE,
            "minimum_primary_orders": 8,
        },
        "file_hashes_lf": input_hashes(base), "source_audit": source_audit,
        "prefit_summary": prefit,
        "panel": panel, "primary": primary, "uncertainty": uncertainty,
        "secondary": secondary, "secondary_uncertainty": secondary_uncertainty,
        "sensitivities": sensitivities,
        "verdict": verdict(primary, uncertainty, sensitivities),
        "numerical_audit": {
            "complete_evaluations": len(evaluations),
            "fits": sum(e["fit_count"] for e in evaluations),
            "all_fits_converged": all(e["all_fits_converged"] for e in evaluations),
            "maximum_fit_iterations": max(e["maximum_fit_iterations"] for e in evaluations),
            "maximum_gradient_ratio": max(e["maximum_gradient_ratio"] for e in evaluations),
            "probabilities_clipped_for_scoring": sum(
                e["probabilities_clipped_for_scoring"] for e in evaluations),
        },
        "limitations": [
            "The evaluated panel is a mammalian subset of a broader vertebrate source.",
            "This tests endpoint portability of a measured proxy, not PR #59's functional-load theorem.",
            "Pathology prevalence is neither onset risk nor cause-specific death.",
            "Shared zoos and possible animal overlap with PR #79 prevent independent-replication claims.",
            "Ascertainment, necropsy availability, age composition and selection remain limitations.",
            "MHC gene-family size does not measure expression or functional antigen presentation.",
            "The small-order jackknife is approximate and does not remove phylogenetic dependence.",
            "Recorded malignant pathology excludes unknown status from its event numerator; unknown is not benign.",
            "PropMalignant is not used and is never a progression probability.",
            "Neither positive nor negative predictive results establish immune causation.",
            "This bounded endpoint audit does not establish historical novelty or biological validation.",
            "Nonclinical unreviewed research; no medical or veterinary advice.",
        ],
    }


def validate_diagnostics(report):
    """Return violations of declared optimizer bounds, never exact stop values."""
    errors = []

    def visit(value, path):
        if isinstance(value, dict):
            for key in DIAGNOSTIC_KEYS & value.keys():
                number = value[key]
                label = path + "." + key
                if (isinstance(number, bool) or not isinstance(number, (int, float))
                        or not math.isfinite(number) or number < 0):
                    errors.append(label + ": invalid diagnostic")
                    continue
                if key == "maximum_gradient_ratio" and number > 1.0 + 1e-12:
                    errors.append(label + ": exceeds convergence ratio")
                elif key == "maximum_absolute_gradient":
                    tolerance = value.get("gradient_tolerance", float("nan"))
                    if not math.isfinite(tolerance) or tolerance <= 0 or number > tolerance * (1 + 1e-12):
                        errors.append(label + ": exceeds gradient bound")
                elif key in ("iterations", "maximum_fit_iterations"):
                    if not isinstance(number, int) or number > MAX_ITERATIONS:
                        errors.append(label + ": exceeds iteration bound")
                elif key == "line_search_halvings":
                    if not isinstance(number, int) or number > 44 * value.get("iterations", -1):
                        errors.append(label + ": exceeds line-search bound")
            for key, item in value.items():
                visit(item, path + "." + key)
        elif isinstance(value, list):
            for index, item in enumerate(value):
                visit(item, path + "[" + str(index) + "]")
        elif isinstance(value, float) and not math.isfinite(value):
            errors.append(path + ": nonfinite value")

    visit(report, "report")
    return errors


def compare_reports(fresh, saved):
    """Compare scientific values tolerantly and optimizer certificates by bounds."""
    errors = ["fresh " + e for e in validate_diagnostics(fresh)]
    errors += ["saved " + e for e in validate_diagnostics(saved)]

    def compare(left, right, label):
        if isinstance(right, dict):
            if not isinstance(left, dict) or set(left) != set(right):
                errors.append(label + ": dictionary keys differ")
                return
            for key in right:
                if key not in DIAGNOSTIC_KEYS:
                    compare(left[key], right[key], label + "." + key)
        elif isinstance(right, list):
            if not isinstance(left, list) or len(left) != len(right):
                errors.append(label + ": list lengths differ")
                return
            for index, (a, b) in enumerate(zip(left, right)):
                compare(a, b, label + "[" + str(index) + "]")
        elif isinstance(right, float):
            if (isinstance(left, bool) or not isinstance(left, (int, float))
                    or not math.isfinite(left) or not math.isfinite(right)
                    or not math.isclose(left, right, rel_tol=SNAPSHOT_TOLERANCE,
                                        abs_tol=SNAPSHOT_TOLERANCE)):
                errors.append(label + ": scientific numeric value differs")
        elif type(left) is not type(right) or left != right:
            errors.append(label + ": value differs")

    compare(fresh, saved, "report")
    return errors


def render_results(report):
    """Render every reported analysis quantity in one deterministic Markdown block."""
    panel = report["panel"]
    primary, uncertainty = report["primary"], report["uncertainty"]
    secondary, secondary_u = report["secondary"], report["secondary_uncertainty"]
    audit = report["numerical_audit"]
    prefit = report["prefit_summary"]
    endpoint_audit = prefit["endpoint_audit"]
    lines = [
        BEGIN, "### Primary result and scope", "",
        "**Verdict: " + report["verdict"]["label"] + ".**", "",
        f"The pinned sources contain **{prefit['vincze_species']} Vincze species** and",
        f"**{prefit['compton_species']} Compton species**; the normalized union has",
        f"**{prefit['union_species']} identifiers**, with **{prefit['normalized_matches']} matches**.",
        f"The predictor-completeness rule excludes **{len(prefit['matched_exclusions'])} matched species**.", "",
        f"The common panel contains **{panel['species']} mammalian species in {panel['orders']} orders**,",
        f"with **{panel['neoplasia']} denominator-linked neoplasms among {panel['records']} records**.",
        f"It retains **{panel['zero_neoplasia_species']} zero-neoplasm species**. The same panel has",
        f"**{panel['malignant']} recorded malignant cases** and **{panel['zero_malignant_species']} species**",
        "with zero recorded malignancies. It does not cover the other vertebrate classes.", "",
        "### Endpoint arithmetic audit", "",
        f"Across **{endpoint_audit['audited_compton_species']} source species**, printed neoplasia prevalence",
        f"has **{endpoint_audit['neoplasia_prevalence_mismatches']} discrepancies** against the denominator-linked",
        f"neoplasia counts; printed malignancy prevalence has **{endpoint_audit['malignant_prevalence_mismatches']} discrepancies**",
        "against `Malignant / RecordsWithDenominators`, under token-specific rounding allowances.",
        f"The dictionary's `MalignancyKnown` numerator disagrees for **{endpoint_audit['malignancyknown_prevalence_mismatches']} source species**",
        f"and **{endpoint_audit['panel_malignancyknown_prevalence_mismatches']} common-panel species**.",
        "This resolves which recorded count reproduces the published prevalence; it does not",
        "reconstruct individual-record provenance. `PropMalignant` is unused and is never a progression probability.", "",
        "### Paired held-out-order scores", "",
        "Losses are natural-log cross-entropy per denominator-linked record.",
        "The primary target averages species within each order, then orders equally.",
        "Delta is baseline minus extended loss; positive values favor adding MHC-I.", "",
        "| Endpoint / scoring target | Baseline | Baseline + MHC-I | Delta |",
        "|---|---:|---:|---:|",
    ]
    for title, result in (("Primary neoplasia", primary), ("Secondary recorded malignancy", secondary)):
        lines.append(f"| {title}: equal orders | {result['baseline_loss']:.9f} | {result['extended_loss']:.9f} | {result['delta']:+.9f} |")
        lines.append(f"| {title}: record weighted | {result['record_weighted_baseline_loss']:.9f} | {result['record_weighted_extended_loss']:.9f} | {result['record_weighted_delta']:+.9f} |")
    lines += [
        "", f"The primary extension improves **{primary['improving_groups']} of {primary['groups']} orders**.",
        "The complete-refit order jackknife reports the following uncertainty.",
        "Its bias-corrected center is distinct from the full-panel point estimate.", "",
        "| Endpoint | Full-panel delta | Jackknife center | Standard error | Approximate 95% interval |",
        "|---|---:|---:|---:|---|",
    ]
    for title, result, interval in (("Primary neoplasia", primary, uncertainty),
                                   ("Secondary recorded malignancy", secondary, secondary_u)):
        lines.append(f"| {title} | {result['delta']:+.9f} | {interval['bias_corrected_center']:+.9f} | {interval['standard_error']:.9f} | [{interval['ci_lower']:+.9f}, {interval['ci_upper']:+.9f}] |")
    lines += [
        "", f"The t critical value is **{uncertainty['t_critical']:.12f}** with",
        f"**{uncertainty['degrees_of_freedom']} degrees of freedom** for each endpoint.",
        "These are approximate intervals over a small selected set of related orders.",
        "The secondary endpoint is exploratory, has no separate multiplicity-adjusted claim,",
        "and cannot rescue the primary verdict. Unknown malignancy status is not known benign status.", "",
    ]
    for title, result in (("Primary neoplasia", primary), ("Secondary recorded malignancy", secondary)):
        lines += [
            "### " + title + ": every held-out order", "",
            "| Order | Species | Records | Events | Zero species | Baseline | + MHC-I | Delta |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for item in result["group_scores"]:
            lines.append(f"| {item['group']} | {item['species']} | {item['records']} | {item['events']} | {item['zero_case_species']} | {item['baseline_loss']:.9f} | {item['extended_loss']:.9f} | {item['delta']:+.9f} |")
        lines.append("")
    lines += [
        "### Prespecified primary sensitivities and falsification attempts", "",
        "| Sensitivity | Species | Orders | Baseline | + MHC-I | Delta |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    labels = {
        "equal_species_training": "Equal-species training weights",
        "penalty_0_1": "Ridge penalty 0.1",
        "penalty_10": "Ridge penalty 10",
        "exclude_named_rodents": "Remove eligible named rodents",
        "exclude_recorded_zeros": "Remove recorded zeros (outcome-selected)",
    }
    for key, label in labels.items():
        item = report["sensitivities"][key]
        lines.append(f"| {label} | {item['species']} | {item['orders']} | {item['baseline_loss']:.9f} | {item['extended_loss']:.9f} | {item['delta']:+.9f} |")
    lines += [
        "", "Named rodents actually removed: " + ", ".join("`" + name + "`" for name in panel["named_rodents_present"]) + ".",
        "Named rodents already absent from the eligible panel: " + ", ".join("`" + name + "`" for name in panel["named_rodents_absent"]) + ".",
        "Only the first four sensitivities are positive-verdict gates. Zero removal is",
        "deliberately outcome-selected and cannot rescue or veto the primary claim.", "",
    ]
    for title, interval in (("Primary neoplasia", uncertainty), ("Secondary recorded malignancy", secondary_u)):
        lines += [
            "### " + title + ": complete delete-order refits", "",
            f"The delete-order delta range is **[{interval['delete_order_delta_min']:+.9f}, {interval['delete_order_delta_max']:+.9f}]**.",
            "Every row reruns all remaining holdouts, their model fits, and their training transformations.", "",
            "| Deleted order | Remaining species | Remaining orders | Refit delta | Pseudovalue |",
            "|---|---:|---:|---:|---:|",
        ]
        for item in interval["deletions"]:
            lines.append(f"| {item['deleted_order']} | {item['species']} | {item['orders']} | {item['delta']:+.9f} | {item['pseudovalue']:+.9f} |")
        lines.append("")
    lines += ["### Frozen decision gates", "", "| Gate | Outcome |", "|---|---|"]
    for key, value in sorted(report["verdict"]["positive_claim_gates"].items()):
        lines.append("| " + key.replace("_", " ") + " | " + str(value).lower() + " |")
    for key in ("point_improvement_failure", "interval_contains_zero", "nonimprovement_gate"):
        lines.append("| " + key.replace("_", " ") + " | " + str(report["verdict"][key]).lower() + " |")
    lines += [
        "", "A nonpositive point gain fails the point-improvement observable. An interval",
        "containing zero leaves the directional claim inconclusive; an upper bound strictly",
        "below zero contradicts improvement under this frozen procedure.", "",
        "### Numerical reproducibility", "",
        f"Across **{audit['complete_evaluations']} complete evaluations and {audit['fits']} fits**:",
        f"all fits converged **{str(audit['all_fits_converged']).lower()}**;",
        f"**{audit['probabilities_clipped_for_scoring']} probabilities** were clipped for scoring.",
        "Each fit must satisfy the declared maximum-gradient bound",
        "`max(1e-8, 1e-10 * effective_training_denominator)` within 100 Newton steps.",
        "Saved stopping diagnostics are checked against these bounds, not their exact",
        "platform-specific values. Scientific snapshots use absolute/relative tolerance `3e-7`.",
        "All primary and secondary predictions, training-fold parameters, coefficients,",
        "fit certificates, refit summaries, and newline-normalized file hashes are in `data/results.json`.", "",
        "**Nonclinical, unreviewed computational research. No medical or veterinary advice.**", END,
    ]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Write recomputed data/results.json")
    parser.add_argument("--check", action="store_true", help="Recompute and check saved results and README")
    parser.add_argument("--markdown", action="store_true", help="Print generated numerical Markdown block")
    args = parser.parse_args(argv)
    if args.write and args.check:
        parser.error("Choose --write or --check, not both")
    report = build_report()
    diagnostic_errors = validate_diagnostics(report)
    if diagnostic_errors:
        raise ValueError("Invalid numerical diagnostics: " + "; ".join(diagnostic_errors))
    path = BASE / "data" / "results.json"
    if args.write:
        with path.open("w", encoding="utf-8", newline="\n") as stream:
            json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    if args.check:
        saved = json.loads(path.read_text(encoding="utf-8"))
        errors = compare_reports(report, saved)
        readme = normalize_newlines((BASE / "README.md").read_text(encoding="utf-8"))
        if readme.count(BEGIN) != 1 or readme.count(END) != 1:
            errors.append("README must contain one complete generated result block")
        else:
            block = readme[readme.index(BEGIN):readme.index(END) + len(END)] + "\n"
            if block != render_results(saved):
                errors.append("README block differs from archived result rendering")
        if errors:
            raise ValueError("Reproduction failed: " + "; ".join(errors))
        print("Recomputed scientific results and archived README block agree; numerical bounds passed.")
    if args.markdown or not (args.write or args.check):
        print(render_results(report), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
