#!/usr/bin/env python3
"""Frozen, nonclinical MHC-I mortality transport audit; Python stdlib only.

No network operations occur in this module. The source retrieval/derivation
workflow is separate in data/fetch_verify.py. All fits use the committed
per-species aggregates and the learning rules specified in PROTOCOL.md.
This is unreviewed computational research, not medical or veterinary advice.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN, localcontext
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics


BASE = Path(__file__).resolve().parent
CONTINUOUS = ("log_body_mass_kg", "log_adult_life_expectancy_days")
MHC = "log1p_mhc_i"
DIET = "mammal_prey"
N = "known_cause_deaths"
K = "nominal_cancer_deaths"
T_15 = 2.131449545559323
P_CLIP = 1e-12
MAX_ITERATIONS = 100
RELATIVE_GRADIENT_TOLERANCE = 1e-10
STATUS = "nonclinical_unreviewed_research_no_medical_or_veterinary_advice"


class FitError(RuntimeError):
    """A numerical failure must stop the relevant analysis, never drop a fold."""


def normalize_newlines(text):
    return text.replace("\r\n", "\n").replace("\r", "\n")


def normalized_sha256(path):
    content = normalize_newlines(Path(path).read_bytes().decode("utf-8"))
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def input_hashes(base=BASE):
    base = Path(base)
    paths = [Path("PROTOCOL.md"), Path("analysis.py")]
    paths += [Path("data") / name for name in (
        "species_aggregates.csv", "exclusions.csv", "source_audit.json",
        "fetch_verify.py")]
    return {path.as_posix(): normalized_sha256(base / path) for path in paths}


def count_interval(token, n):
    """Independent offline implementation of the frozen decimal count audit."""
    with localcontext() as context:
        context.prec = 60
        c = Decimal(token)
        if not c.is_finite() or not Decimal(0) <= c <= Decimal(1):
            raise ValueError("Invalid CMR probability")
        if n < 1 or int(n) != n:
            raise ValueError("RDI must be a positive integer")
        nominal = int((c * n).to_integral_value(rounding=ROUND_HALF_EVEN))
        if c in (Decimal(0), Decimal(1)):
            return nominal, nominal, nominal
        width = Decimal(1).scaleb(c.as_tuple().exponent) / 2 + Decimal("2e-15")
        lo = int((max(Decimal(0), c - width) * n).to_integral_value(
            rounding=ROUND_CEILING))
        hi = int((min(Decimal(1), c + width) * n).to_integral_value(
            rounding=ROUND_FLOOR))
        if not lo <= nominal <= hi:
            raise ValueError("No rounding-compatible nominal count")
        return nominal, lo, hi


def load_rows(base=BASE):
    """Read and validate derived inputs, without fitting or altering the panel."""
    base = Path(base)
    with (base / "data" / "species_aggregates.csv").open(
            encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    integers = (N, K, "count_lower", "count_upper", DIET)
    floats = CONTINUOUS + (MHC, "phylo_fair_proportion")
    for row in rows:
        for name in integers:
            row[name] = int(row[name])
        for name in floats:
            row[name] = float(row[name])
            if not math.isfinite(row[name]):
                raise ValueError("Nonfinite derived predictor: " + row["species"])
        if row[DIET] not in (0, 1) or row[MHC] < 0:
            raise ValueError("Invalid diet/MHC predictor")
        if row["phylo_fair_proportion"] <= 0:
            raise ValueError("Nonpositive phylogenetic weight")
        recovered = count_interval(row["cmr_decimal"], row[N])
        if recovered != (row[K], row["count_lower"], row["count_upper"]):
            raise ValueError("Count audit mismatch: " + row["species"])
        if not row["species"] or not row["order"]:
            raise ValueError("Missing species/order identifier")
    names = [row["species"] for row in rows]
    if names != sorted(names) or len(names) != len(set(names)):
        raise ValueError("Species must be unique and sorted")
    return rows


def sigmoid(value):
    if value >= 0:
        inverse = math.exp(-value)
        return 1.0 / (1.0 + inverse)
    exponent = math.exp(value)
    return exponent / (1.0 + exponent)


def softplus(value):
    return max(value, 0.0) + math.log1p(math.exp(-abs(value)))


def cancer_death_loss(k, n, prediction):
    """Per-death cross-entropy, including the k=0 case without pseudocounts."""
    if n <= 0 or not 0 <= k <= n:
        raise ValueError("Invalid cancer-death count or denominator")
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


def fit_model(train, include_mhc, penalty=1.0, weighting="counts"):
    """Fit the frozen ridge-binomial rule using training species only.

    `species` equalizes effective denominators to the training mean n.
    `phylogeny` multiplies actual n and k by normalized fair-proportion weights.
    The score always uses the original denominators and unaltered test labels.
    """
    if not train or penalty <= 0 or not math.isfinite(penalty):
        raise ValueError("Training rows and a positive finite ridge are required")
    if weighting not in ("counts", "species", "phylogeny"):
        raise ValueError("Unknown training weighting")
    continuous = CONTINUOUS + ((MHC,) if include_mhc else ())
    means = {key: statistics.fmean(row[key] for row in train) for key in continuous}
    scales = {key: math.sqrt(statistics.fmean(
        (row[key] - means[key]) ** 2 for row in train)) or 1.0
        for key in continuous}
    design = [_feature_vector(row, means, scales, include_mhc) for row in train]
    mean_n = statistics.fmean(row[N] for row in train)
    mean_phylogeny = statistics.fmean(row["phylo_fair_proportion"] for row in train)
    denominators = []
    cancer = []
    for row in train:
        multiplier = 1.0
        if weighting == "species":
            multiplier = mean_n / row[N]
        elif weighting == "phylogeny":
            multiplier = row["phylo_fair_proportion"] / mean_phylogeny
        denominators.append(row[N] * multiplier)
        cancer.append(row[K] * multiplier)
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
                "feature_names": ["intercept", *CONTINUOUS, DIET] +
                                 ([MHC] if include_mhc else []),
                "coefficients": coefficients,
                "means": means,
                "scales": scales,
                "penalty": penalty,
                "weighting": weighting,
                "training_species": len(train),
                "training_known_cause_deaths": sum(row[N] for row in train),
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


def evaluate(rows, penalty=1.0, weighting="counts", group_key="order"):
    """Paired whole-group holdout; no held-out labels reach either fit."""
    labels = sorted({row[group_key] for row in rows})
    if len(labels) < 2:
        raise FitError("At least two held-out groups are required")
    predictions = []
    folds = []
    group_scores = []
    clipping = 0
    for label in labels:
        train = [row for row in rows if row[group_key] != label]
        test = [row for row in rows if row[group_key] == label]
        baseline = fit_model(train, False, penalty, weighting)
        extended = fit_model(train, True, penalty, weighting)
        fold_predictions = []
        for row in test:
            p0, p1 = predict_one(row, baseline), predict_one(row, extended)
            clipping += int(p0 < P_CLIP or p0 > 1 - P_CLIP)
            clipping += int(p1 < P_CLIP or p1 > 1 - P_CLIP)
            loss0 = cancer_death_loss(row[K], row[N], p0)
            loss1 = cancer_death_loss(row[K], row[N], p1)
            record = {
                "species": row["species"], "order": row["order"],
                "heldout_group": label, N: row[N], K: row[K],
                "baseline_probability": p0, "extended_probability": p1,
                "baseline_loss": loss0, "extended_loss": loss1,
                "delta": loss0 - loss1,
            }
            predictions.append(record)
            fold_predictions.append(record)
        baseline_loss = statistics.fmean(x["baseline_loss"] for x in fold_predictions)
        extended_loss = statistics.fmean(x["extended_loss"] for x in fold_predictions)
        group_scores.append({
            "group": label, "species": len(test),
            "orders": sorted({row["order"] for row in test}),
            "known_cause_deaths": sum(row[N] for row in test),
            "nominal_cancer_deaths": sum(row[K] for row in test),
            "zero_case_species": sum(row[K] == 0 for row in test),
            "baseline_loss": baseline_loss, "extended_loss": extended_loss,
            "delta": baseline_loss - extended_loss,
        })
        folds.append({
            "heldout_group": label,
            "training_species_ids": [row["species"] for row in train],
            "test_species_ids": [row["species"] for row in test],
            "training_orders": sorted({row["order"] for row in train}),
            "test_orders": sorted({row["order"] for row in test}),
            "baseline_fit": baseline, "extended_fit": extended,
        })
    total_n = sum(row[N] for row in rows)
    baseline_loss = statistics.fmean(g["baseline_loss"] for g in group_scores)
    extended_loss = statistics.fmean(g["extended_loss"] for g in group_scores)
    death_baseline = math.fsum(p[N] * p["baseline_loss"] for p in predictions) / total_n
    death_extended = math.fsum(p[N] * p["extended_loss"] for p in predictions) / total_n
    diagnostics = [fold[key] for fold in folds for key in ("baseline_fit", "extended_fit")]
    return {
        "group_key": group_key, "weighting": weighting, "penalty": penalty,
        "species": len(rows), "groups": len(labels),
        "orders": len({row["order"] for row in rows}),
        "known_cause_deaths": total_n,
        "nominal_cancer_deaths": sum(row[K] for row in rows),
        "zero_case_species": sum(row[K] == 0 for row in rows),
        "baseline_loss": baseline_loss, "extended_loss": extended_loss,
        "delta": baseline_loss - extended_loss,
        "death_weighted_baseline_loss": death_baseline,
        "death_weighted_extended_loss": death_extended,
        "death_weighted_delta": death_baseline - death_extended,
        "improving_groups": sum(g["delta"] > 0 for g in group_scores),
        "probabilities_clipped_for_scoring": clipping,
        "fit_count": len(diagnostics),
        "maximum_fit_iterations": max(f["iterations"] for f in diagnostics),
        "maximum_gradient_ratio": max(f["maximum_absolute_gradient"] /
                                      f["gradient_tolerance"] for f in diagnostics),
        "all_fits_converged": all(f["converged"] for f in diagnostics),
        "group_scores": group_scores,
        "predictions": sorted(predictions, key=lambda item: item["species"]),
        "folds": folds,
    }


def _summary(evaluation):
    return {key: value for key, value in evaluation.items()
            if key not in ("predictions", "folds")}


def jackknife_uncertainty(rows, primary):
    """Delete a whole order and refit every remaining outer fold."""
    groups = sorted({row["order"] for row in rows})
    if len(groups) != 16:
        raise ValueError("Frozen t15 interval requires the frozen 16-order panel")
    count = len(groups)
    deletions = []
    for group in groups:
        retained = [row for row in rows if row["order"] != group]
        evaluated = evaluate(retained)
        deletions.append({"deleted_order": group, **_summary(evaluated)})
    pseudovalues = [count * primary["delta"] - (count - 1) * item["delta"]
                    for item in deletions]
    center = statistics.fmean(pseudovalues)
    se = math.sqrt(math.fsum((q - center) ** 2 for q in pseudovalues) /
                   (count * (count - 1)))
    for item, pseudovalue in zip(deletions, pseudovalues):
        item["pseudovalue"] = pseudovalue
    return {
        "method": "complete_delete_one_order_refit_jackknife",
        "groups": count, "degrees_of_freedom": count - 1,
        "t_critical": T_15, "full_panel_point_estimate": primary["delta"],
        "bias_corrected_center": center, "standard_error": se,
        "ci_lower": center - T_15 * se, "ci_upper": center + T_15 * se,
        "delete_order_delta_min": min(item["delta"] for item in deletions),
        "delete_order_delta_max": max(item["delta"] for item in deletions),
        "deletions": deletions,
    }


def verdict(primary, uncertainty, rounding, sensitivities, panel_valid=True):
    """Apply the exhaustive frozen directional rule, without fitted thresholds."""
    valid = (panel_valid and primary["orders"] >= 8 and
             primary["all_fits_converged"] and
             len(uncertainty["deletions"]) == primary["orders"] and
             all(item["all_fits_converged"] for item in uncertainty["deletions"]))
    major = ("exclude_named_rodents", "equal_species_training", "penalty_0_1",
             "penalty_10", "phylogeny_weighting")
    deep = sensitivities["deeper_phylogeny_holdout"]
    deeper_gate_applies = deep["groups"] >= 4
    gates = {
        "valid_panel_and_primary_uncertainty_fits": valid,
        "positive_primary_delta": primary["delta"] > 0,
        "positive_jackknife_lower_bound": uncertainty["ci_lower"] > 0,
        "majority_orders_improve": primary["improving_groups"] > primary["groups"] / 2,
        "all_rounding_assignments_positive": (len(rounding) == 15 and all(
            item["delta"] > 0 and item["all_fits_converged"] for item in rounding)),
        "all_major_sensitivities_positive": all(
            sensitivities[key]["delta"] > 0 and sensitivities[key]["all_fits_converged"]
            for key in major),
        "deeper_tree_gate_applies": deeper_gate_applies,
        "deeper_tree_gate_passed": (not deeper_gate_applies or
                                    (deep["delta"] > 0 and deep["all_fits_converged"])),
    }
    positive_keys = [key for key in gates if key != "deeper_tree_gate_applies"]
    if not valid:
        code = "untestable"
        label = "Untestable (panel or primary/uncertainty fit validity failure)"
    elif all(gates[key] for key in positive_keys):
        code = "adds"
        label = "MHC-I adds portable predictive information within this evaluation"
    elif primary["delta"] <= 0 and uncertainty["ci_upper"] <= 0:
        code = "does_not_add"
        label = "MHC-I does not add portable predictive information within this evaluation"
    else:
        code = "untestable"
        label = "Untestable (inconclusive on the available panel)"
    return {"code": code, "label": label, "positive_claim_gates": gates,
            "non_improvement_gate": (valid and primary["delta"] <= 0 and
                                     uncertainty["ci_upper"] <= 0)}


def _panel_summary(rows, exclusions):
    ambiguous = [{"species": row["species"], "lower": row["count_lower"],
                  "nominal": row[K], "upper": row["count_upper"], N: row[N]}
                 for row in rows if row["count_lower"] != row["count_upper"]]
    residuals = [(abs(Decimal(row["cmr_decimal"]) * row[N] - row[K]), row["species"])
                 for row in rows]
    largest = max(residuals)
    return {
        "source_rows": len(rows) + len(exclusions),
        "species": len(rows), "excluded_species": len(exclusions),
        "orders": len({row["order"] for row in rows}),
        "known_cause_deaths": sum(row[N] for row in rows),
        "nominal_cancer_deaths": sum(row[K] for row in rows),
        "count_lower_total": sum(row["count_lower"] for row in rows),
        "count_upper_total": sum(row["count_upper"] for row in rows),
        "zero_case_species": sum(row[K] == 0 for row in rows),
        "unique_count_species": len(rows) - len(ambiguous),
        "ambiguous_count_species": ambiguous,
        "largest_integer_residual": str(largest[0]),
        "largest_integer_residual_species": largest[1],
        "sum_fair_proportion": math.fsum(row["phylo_fair_proportion"] for row in rows),
        "species_per_order": {order: sum(row["order"] == order for row in rows)
                              for order in sorted({row["order"] for row in rows})},
        "species_per_phylo_block": {block: sum(row["phylo_block"] == block for row in rows)
                                    for block in sorted({row["phylo_block"] for row in rows})},
    }


def build_report(base=BASE):
    """Recompute every analysis result from the committed local data, offline."""
    base = Path(base)
    rows = load_rows(base)
    with (base / "data" / "exclusions.csv").open(encoding="utf-8", newline="") as handle:
        exclusions = list(csv.DictReader(handle))
    source_audit = json.loads((base / "data" / "source_audit.json").read_text(encoding="utf-8"))
    for relative, expected in source_audit["derived_file_hashes_lf"].items():
        if normalized_sha256(base / "data" / Path(relative)) != expected:
            raise ValueError("Derived input hash differs from source audit: " + relative)
    panel = _panel_summary(rows, exclusions)
    for field, expected in (("species", 125), ("orders", 16),
                            ("zero_case_species", 17), ("known_cause_deaths", 12033),
                            ("nominal_cancer_deaths", 697), ("source_rows", 139)):
        if panel[field] != expected:
            raise ValueError("Frozen panel gate failed: " + field)
    required = {"Heterocephalus_glaber", "Dolichotis_patagonum"}
    if not required <= {row["species"] for row in rows}:
        raise ValueError("Required primary rodents are absent")
    primary = evaluate(rows)
    uncertainty = jackknife_uncertainty(rows, primary)
    ambiguous = [row for row in rows if row["count_lower"] != row["count_upper"]]
    rounding = []
    for counts in itertools.product(*(range(row["count_lower"], row["count_upper"] + 1)
                                      for row in ambiguous)):
        assignment = dict(zip((row["species"] for row in ambiguous), counts))
        scenario = [dict(row, **{K: assignment.get(row["species"], row[K])}) for row in rows]
        evaluated = evaluate(scenario)
        rounding.append({"assignment": assignment, **_summary(evaluated)})
    if len(rounding) != 15:
        raise ValueError("Frozen rounding grid does not contain 15 scenarios")
    sensitivities = {
        "exclude_named_rodents": _summary(evaluate(
            [row for row in rows if row["species"] not in required])),
        "exclude_recorded_zeros": _summary(evaluate([row for row in rows if row[K] > 0])),
        "equal_species_training": _summary(evaluate(rows, weighting="species")),
        "penalty_0_1": _summary(evaluate(rows, penalty=0.1)),
        "penalty_10": _summary(evaluate(rows, penalty=10.0)),
        "phylogeny_weighting": _summary(evaluate(rows, weighting="phylogeny")),
        "deeper_phylogeny_holdout": _summary(evaluate(rows, group_key="phylo_block")),
    }
    complete_evaluations = [primary, *uncertainty["deletions"], *rounding,
                            *sensitivities.values()]
    return {
        "schema_version": 1, "status": STATUS,
        "protocol": "PROTOCOL.md",
        "analysis_settings": {
            "primary_penalty": 1.0, "score_probability_clip": P_CLIP,
            "relative_gradient_tolerance": RELATIVE_GRADIENT_TOLERANCE,
            "maximum_newton_iterations": MAX_ITERATIONS,
            "uncertainty_t_critical": T_15,
        },
        "file_hashes_lf": input_hashes(base),
        "source_audit": source_audit,
        "panel": panel,
        "primary": primary,
        "uncertainty": uncertainty,
        "rounding": rounding,
        "rounding_summary": {
            "scenarios": len(rounding),
            "minimum_delta": min(item["delta"] for item in rounding),
            "maximum_delta": max(item["delta"] for item in rounding),
            "positive_scenarios": sum(item["delta"] > 0 for item in rounding),
        },
        "sensitivities": sensitivities,
        "verdict": verdict(primary, uncertainty, rounding, sensitivities),
        "numerical_audit": {
            "complete_evaluations": len(complete_evaluations),
            "fits": sum(item["fit_count"] for item in complete_evaluations),
            "all_fits_converged": all(item["all_fits_converged"] for item in complete_evaluations),
            "maximum_fit_iterations": max(item["maximum_fit_iterations"]
                                          for item in complete_evaluations),
            "maximum_gradient_ratio": max(item["maximum_gradient_ratio"]
                                          for item in complete_evaluations),
            "probabilities_clipped_for_scoring": sum(
                item["probabilities_clipped_for_scoring"] for item in complete_evaluations),
        },
        "limitations": [
            "The association was published previously; novelty is incremental out-of-order prediction.",
            "Mortality provenance overlaps PR #75; this is not independent mortality validation.",
            "Gene-family size neither measures functional repertoires nor identifies immune causation.",
            "Nominal counts are rounding-conditional; all compatible count assignments are reported.",
            "The jackknife interval is approximate for selected, related, unevenly sampled orders.",
            "Phylogenetic reweighting is not covariance adjustment; alternative trees remain untested.",
            "Known-cause deaths in managed populations do not represent all mammalian populations.",
            "Cause-of-death error, ascertainment, gene annotation quality and sample selection remain.",
            "Reproducible arithmetic and CI do not establish source truth or biological validity.",
            "Nonclinical unreviewed research; no medical or veterinary advice.",
        ],
    }


def render_results(report):
    """Generate the complete README numerical result block deterministically."""
    panel, primary, uncertainty = report["panel"], report["primary"], report["uncertainty"]
    grid = report["rounding_summary"]
    audit = report["numerical_audit"]
    strongest_order = max(primary["group_scores"], key=lambda item: item["delta"])
    strongest_deleted = next(item for item in uncertainty["deletions"]
                             if item["deleted_order"] == strongest_order["group"])
    minimum_deletion = min(uncertainty["deletions"], key=lambda item: item["delta"])
    if primary["delta"] > 0 and minimum_deletion["delta"] < 0:
        deletion_note = (f"Deleting {minimum_deletion['deleted_order']} and refitting reverses the "
                         f"full-panel improvement: delta **{minimum_deletion['delta']:+.9f}**.")
    else:
        deletion_note = (f"The smallest complete-refit delete-order estimate occurs after deleting "
                         f"{minimum_deletion['deleted_order']}: delta "
                         f"**{minimum_deletion['delta']:+.9f}**.")
    deeper = report["sensitivities"]["deeper_phylogeny_holdout"]
    block_sizes = ", ".join(str(item["species"]) for item in deeper["group_scores"])
    if deeper["groups"] < 4:
        deeper_note = (f"The fixed tree cut yields **{deeper['groups']} blocks**, containing "
                       f"**{block_sizes} species**. This is below the frozen 4-block gate. "
                       "The deeper-holdout result is descriptive and cannot confirm transport.")
    else:
        deeper_note = (f"The fixed tree cut yields **{deeper['groups']} blocks**, containing "
                       f"**{block_sizes} species**. The frozen 4-block gate applies.")
    lines = [
        "<!-- BEGIN GENERATED RESULTS -->",
        "### Primary result",
        "",
        "**Verdict: " + report["verdict"]["label"] + ".**",
        "",
        f"The common panel contains **{panel['species']} species in {panel['orders']} orders**,",
        f"including **{panel['zero_case_species']} recorded zero-case species**. Its",
        f"**{panel['known_cause_deaths']} known-cause deaths** contain **{panel['nominal_cancer_deaths']}",
        f"nominal reconstructed cancer deaths**, with rounding-compatible totals of",
        f"**{panel['count_lower_total']}–{panel['count_upper_total']}**. The",
        f"{panel['source_rows']}-row source excludes {panel['excluded_species']} rows under the frozen rule.",
        f"Count reconstruction is unique for {panel['unique_count_species']} species; the largest",
        f"nearest-integer discrepancy is {panel['largest_integer_residual']} deaths",
        f"for `{panel['largest_integer_residual_species']}`.",
        "",
        "### Paired predictive scores",
        "",
        "All losses are natural-log cross-entropy units per known-cause death.",
        "Positive delta means the MHC extension has lower loss.",
        "",
        "| Estimand | Baseline loss | Baseline + MHC-I loss | Delta |",
        "|---|---:|---:|---:|",
        f"| Equal orders; species averaged within order | {primary['baseline_loss']:.9f} | {primary['extended_loss']:.9f} | {primary['delta']:+.9f} |",
        f"| Death weighted | {primary['death_weighted_baseline_loss']:.9f} | {primary['death_weighted_extended_loss']:.9f} | {primary['death_weighted_delta']:+.9f} |",
        "",
        f"MHC improves {primary['improving_groups']} of {primary['groups']} held-out orders.",
        f"The complete-refit order jackknife has bias-corrected center",
        f"**{uncertainty['bias_corrected_center']:+.9f}**, standard error",
        f"**{uncertainty['standard_error']:.9f}**, and approximate 95% interval",
        f"**[{uncertainty['ci_lower']:+.9f}, {uncertainty['ci_upper']:+.9f}]**",
        f"(t critical {uncertainty['t_critical']:.12f}; {uncertainty['degrees_of_freedom']} degrees of freedom).",
        "The full-panel point estimate and the bias-corrected interval center are different estimators.",
        "The interval is approximate and does not establish independence across orders.",
        "",
        f"The largest held-out-order gain is {strongest_order['group']} "
        f"({strongest_order['species']} species): delta **{strongest_order['delta']:+.9f}**.",
        f"Removing that order and refitting yields delta **{strongest_deleted['delta']:+.9f}**.",
        deletion_note,
        "",
        "### Every held-out order",
        "",
        "| Order | Species | Known-cause deaths | Nominal cancer deaths | Zero species | Baseline | + MHC-I | Delta |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in primary["group_scores"]:
        lines.append(f"| {item['group']} | {item['species']} | {item['known_cause_deaths']} | {item['nominal_cancer_deaths']} | {item['zero_case_species']} | {item['baseline_loss']:.9f} | {item['extended_loss']:.9f} | {item['delta']:+.9f} |")
    lines += [
        "",
        "### Prespecified sensitivities and falsification attempts",
        "",
        "| Sensitivity | Species | Holdout groups | Baseline | + MHC-I | Delta |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    labels = {
        "exclude_named_rodents": "Remove naked mole-rat and Patagonian mara",
        "exclude_recorded_zeros": "Remove recorded zeros (outcome-selected)",
        "equal_species_training": "Equal-species training weights",
        "penalty_0_1": "Ridge penalty 0.1",
        "penalty_10": "Ridge penalty 10",
        "phylogeny_weighting": "Fair-proportion phylogenetic training weights",
        "deeper_phylogeny_holdout": "Deeper tree-block holdout",
    }
    for name, label in labels.items():
        item = report["sensitivities"][name]
        lines.append(f"| {label} | {item['species']} | {item['groups']} | {item['baseline_loss']:.9f} | {item['extended_loss']:.9f} | {item['delta']:+.9f} |")
    lines += [
        "",
        "Zero removal is deliberately outcome-selected and cannot rescue or veto the primary claim.",
        "Tree weighting changes evolutionary representation in training; it is not covariance adjustment.",
        "The deeper-tree result averages tree components equally and therefore changes the transport estimand.",
        "",
        deeper_note,
        "",
        "### Deeper tree-block holdout",
        "",
        "| Block | Included orders | Species | Baseline | + MHC-I | Delta |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for item in report["sensitivities"]["deeper_phylogeny_holdout"]["group_scores"]:
        orders = ", ".join(item["orders"])
        lines.append(f"| {item['group']} | {orders} | {item['species']} | {item['baseline_loss']:.9f} | {item['extended_loss']:.9f} | {item['delta']:+.9f} |")
    lines += [
        "",
        "### Complete rounding grid",
        "",
        f"All {grid['scenarios']} assignments are reported. The delta range is",
        f"**[{grid['minimum_delta']:+.9f}, {grid['maximum_delta']:+.9f}]**;",
        f"{grid['positive_scenarios']} assignments have positive delta.",
        "",
        "| Eulemur macaco cancer deaths | Myrmecophaga tridactyla cancer deaths | Total nominal cancer deaths | Baseline | + MHC-I | Delta |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for item in report["rounding"]:
        assignment = item["assignment"]
        lines.append(f"| {assignment['Eulemur_macaco']} | {assignment['Myrmecophaga_tridactyla']} | {item['nominal_cancer_deaths']} | {item['baseline_loss']:.9f} | {item['extended_loss']:.9f} | {item['delta']:+.9f} |")
    lines += [
        "",
        "### Complete delete-order refits",
        "",
        f"The full-refit delete-order delta range is",
        f"**[{uncertainty['delete_order_delta_min']:+.9f}, {uncertainty['delete_order_delta_max']:+.9f}]**.",
        "Each row reruns fitting, training-only scaling, and all remaining holdouts.",
        "",
        "| Deleted order | Remaining species | Remaining orders | Refit delta | Jackknife pseudovalue |",
        "|---|---:|---:|---:|---:|",
    ]
    for item in uncertainty["deletions"]:
        lines.append(f"| {item['deleted_order']} | {item['species']} | {item['groups']} | {item['delta']:+.9f} | {item['pseudovalue']:+.9f} |")
    lines += [
        "",
        "### Frozen decision gates and numerical audit",
        "",
        "| Positive-claim gate | Outcome |",
        "|---|---|",
    ]
    for gate, passed in sorted(report["verdict"]["positive_claim_gates"].items()):
        lines.append("| " + gate.replace("_", " ") + " | " + str(passed).lower() + " |")
    lines += [
        "",
        f"The non-improvement gate is **{str(report['verdict']['non_improvement_gate']).lower()}**.",
        f"Across {audit['complete_evaluations']} complete evaluations and {audit['fits']} fits,",
        f"all fits converged: **{str(audit['all_fits_converged']).lower()}**.",
        f"Maximum Newton iterations: {audit['maximum_fit_iterations']}; maximum",
        f"gradient/tolerance ratio: {audit['maximum_gradient_ratio']:.9f}.",
        f"Probabilities clipped for scoring: {audit['probabilities_clipped_for_scoring']}.",
        "All primary held-out species probabilities, training-fold transformations, coefficients,",
        "convergence diagnostics, and newline-normalized input/code hashes are in `data/results.json`.",
        "",
        "**Nonclinical, unreviewed computational research. No medical or veterinary advice.**",
        "<!-- END GENERATED RESULTS -->",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Write data/results.json")
    parser.add_argument("--markdown", action="store_true", help="Print the README numerical block")
    args = parser.parse_args(argv)
    report = build_report()
    if args.write:
        path = BASE / "data" / "results.json"
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(report, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
    if args.markdown or not args.write:
        print(render_results(report), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
