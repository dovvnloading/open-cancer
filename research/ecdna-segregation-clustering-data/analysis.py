"""Aggregate-only ecDNA analysis. Nonclinical; no treatment advice.

Source estimates are retained as published estimates. The optional Fisher
calculation is a conditional experimental-design calculation under a normalized
reconstruction, NOT a recovered empirical confidence interval or a refit of
the original observations. All inputs are tiny published aggregates in data/.
Only Python's standard library is used; no network or external data is read.
"""

import argparse
import csv
from fractions import Fraction
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path
from statistics import NormalDist, mean, median


HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
F = Fraction


def load_csv(filename):
    with (DATA / filename).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


@lru_cache(maxsize=1)
def verifier_module():
    spec = importlib.util.spec_from_file_location(
        "clustering_data_verifier", HERE / "verifier.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def effective_r(p, rho):
    """Effective sister-pooled dispersion, including squared allocation bias."""
    return (2*p-1)**2 + 4*p*(1-p)*rho


def binomial_tail(n, count, probability):
    if count <= 0:
        return 1.0
    if count > n:
        return 0.0
    return math.fsum(math.comb(n, j)*probability**j*(1-probability)**(n-j)
                     for j in range(count, n+1))


def compatible_counts(n, printed_percent):
    """Outward whole-percentage rounding; do not manufacture a denominator."""
    center = F(printed_percent)
    return [j for j in range(n+1)
            if center-F(1, 2) <= F(100*j, n) <= center+F(1, 2)]


def wilson(n, count, alpha=0.05):
    z = NormalDist().inv_cdf(1-alpha/2)
    phat = count/n
    denominator = 1+z*z/n
    center = (phat+z*z/(2*n))/denominator
    radius = z*math.sqrt(phat*(1-phat)/n+z*z/(4*n*n))/denominator
    return [max(0.0, center-radius), min(1.0, center+radius)]


def holm(values):
    adjusted = [0.0]*len(values)
    previous = 0.0
    for order, index in enumerate(sorted(range(len(values)), key=values.__getitem__)):
        previous = max(previous, min(1.0, (len(values)-order)*values[index]))
        adjusted[index] = previous
    return adjusted


def solve(matrix, rhs):
    """Small floating-point solver for the conditional information calculation."""
    n = len(rhs)
    a = [list(row)+[value] for row, value in zip(matrix, rhs)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(a[i][k]))
        if abs(a[pivot][k]) < 1e-14:
            raise ArithmeticError("Information is singular; uncertainty is unidentified")
        a[k], a[pivot] = a[pivot], a[k]
        divisor = a[k][k]
        a[k] = [value/divisor for value in a[k]]
        for i in range(n):
            if i != k:
                multiple = a[i][k]
                a[i] = [x-multiple*y for x, y in zip(a[i], a[k])]
    return [row[-1] for row in a]


def nb_positive_rows(size, mean_total, tolerance=1e-11):
    """Positive-truncated NB mass and scores in log(size), log(mean).

    The source fits are inputs to this declared normalized reconstruction.
    Normalizing is NOT a claim that the source's printed likelihood did so.
    """
    q0 = (size/(size+mean_total))**size
    q = q0
    prefix = 0.0
    mass = 0.0
    s0_size = math.log(size/(size+mean_total)) + 1-size/(size+mean_total)
    s0_mean = -size/(size+mean_total)
    rows = []
    for n in range(1, 10001):
        q *= (size+n-1)/n * mean_total/(size+mean_total)
        prefix += 1/(size+n-1)
        probability = q/(1-q0)
        score_size = size*(prefix+math.log(size/(size+mean_total))+1
                           -(size+n)/(size+mean_total)+q0*s0_size/(1-q0))
        score_mean = mean_total*(n/mean_total-(size+n)/(size+mean_total)
                                +q0*s0_mean/(1-q0))
        rows.append((n, probability, score_size, score_mean))
        mass += probability
        if n > mean_total and 1-mass < tolerance:
            return rows, q0, max(0.0, 1-mass)
    raise ArithmeticError("NB tail failed to converge")


def bb_mass_scores(n, p, rho):
    """BB masses and analytic scores in (R,B), B=(2p-1)^2.

    Every source input has 0<p<1/2, 0<rho<1. The squared-bias coordinate
    avoids the first-order loss of identifiability of p near 1/2.
    """
    if not (0 < p < 0.5 and 0 < rho < 1):
        raise ValueError("Information calculation requires an interior source fit")
    concentration = 1/rho-1
    a, b = p*concentration, (1-p)*concentration
    mass = math.prod((b+i)/(concentration+i) for i in range(n))
    probabilities = [mass]
    for j in range(n):
        probabilities.append(probabilities[-1]*(n-j)/(j+1)
                             *(a+j)/(b+n-j-1))
    ha, hb = [0.0], [0.0]
    for i in range(n):
        ha.append(ha[-1]+1/(a+i))
        hb.append(hb[-1]+1/(b+i))
    hs = math.fsum(1/(concentration+i) for i in range(n))
    bias = (2*p-1)**2
    R = effective_r(p, rho)
    scores = []
    for j in range(n+1):
        dp = concentration*(ha[j]-hb[n-j])
        drho = -(p*ha[j]+(1-p)*hb[n-j]-hs)/(rho*rho)
        dR = drho/(1-bias)
        dB = -dp/(4*math.sqrt(bias)) + drho*(R-1)/(1-bias)**2
        scores.append((dR, dB))
    return probabilities, scores


@lru_cache(maxsize=None)
def design_information(p, rho, size, mean_total, tolerance=1e-11):
    """Expected information of a *hypothetical* sorted/floored sister pair.

    Latent total N is positive-truncated NB. Draw J~BB(N,p,rho); observe
    sort(max(J,1),max(N-J,1)). This fixes the printed source's inconsistent
    total in its zero-floor branch, and is explicitly a different normalized
    reconstruction. Nuisance size/mean are jointly estimated. Independence,
    equal observation, constancy of parameters and this model are assumptions.
    """
    nb_rows, q0, tail = nb_positive_rows(size, mean_total, tolerance)
    information = [[0.0]*4 for _ in range(4)]
    total_score = [0.0]*4
    mass_sum = 0.0
    boundary = {}

    def add(probability, derivatives):
        nonlocal mass_sum
        mass_sum += probability
        for i in range(4):
            total_score[i] += derivatives[i]
            for j in range(i+1):
                information[i][j] += derivatives[i]*derivatives[j]/probability

    for n, pn, sn_size, sn_mean in nb_rows:
        masses, scores = bb_mass_scores(n, p, rho)
        for j in range(n//2+1):
            reflected = n-j
            terms = [j] if reflected == j else [j, reflected]
            probability = pn*math.fsum(masses[index] for index in terms)
            derivatives = [pn*math.fsum(masses[index]*scores[index][dim]
                                        for index in terms) for dim in (0, 1)]
            derivatives.extend([probability*sn_size, probability*sn_mean])
            if j <= 1:
                # n=y and j=0 has the same observed pair as n=y+1,j=1.
                y = max(1, reflected)
                item = boundary.setdefault(y, [0.0]*5)
                item[0] += probability
                for index in range(4):
                    item[index+1] += derivatives[index]
            else:
                add(probability, derivatives)
    for item in boundary.values():
        add(item[0], item[1:])
    for i in range(4):
        for j in range(i):
            information[j][i] = information[i][j]
    inverse_first_column = solve(information, [1.0, 0.0, 0.0, 0.0])
    inverse_second_column = solve(information, [0.0, 1.0, 0.0, 0.0])
    if inverse_first_column[0] <= 0:
        raise ArithmeticError("Nonpositive profiled information")
    return {"information": information,
            "r_variance_per_pair": inverse_first_column[0],
            "bias2_variance_per_pair": inverse_second_column[1],
            "zero_total_mass": q0, "included_mass": mass_sum,
            "tail_mass": tail, "max_total": nb_rows[-1][0],
            "score_sum": total_score}


def conditional_design_interval(row, effective_pairs=None, alpha=0.05):
    """Local Fisher sensitivity band; nominal z does not establish coverage.

    The squared-bias nuisance is weakly identified in these fits. Its local
    standard error often exceeds its entire parameter domain. Consequently
    even finite-sample coverage under the reconstruction is not claimed.
    """
    p, rho = float(row["p"]), float(row["rho"])
    information = design_information(p, rho, float(row["nb_size"]), float(row["nb_mean"]))
    pairs = effective_pairs
    if pairs is None:
        pairs = int(row["n_reported"])/(2 if row["id"].startswith("shibata_clone") else 1)
    if pairs <= 0:
        raise ValueError("Positive assumed effective pair count required")
    R = effective_r(p, rho)
    se = math.sqrt(information["r_variance_per_pair"]/pairs)
    bias_se = math.sqrt(information["bias2_variance_per_pair"]/pairs)
    z = NormalDist().inv_cdf(1-alpha/2)
    return {"assumed_independent_pairs": pairs, "r_se": se,
            "r_interval": [max(0.0, R-z*se), min(1.0, R+z*se)],
            "bias2_se": bias_se, "nominal_alpha": alpha,
            "label": "local Fisher sensitivity band; no validated sampling coverage",
            "zero_total_mass": information["zero_total_mass"],
            "omitted_nb_mass": information["tail_mass"]}


def read_lundberg():
    rows = load_csv("lundberg.csv")
    by_population = {}
    for row in rows:
        by_population.setdefault((row["population"], row["condition"]), {})[row["quantity"]] = row
    return rows, by_population


def source_model_comparison():
    rows = load_csv("shibata_models.csv")
    minimum = min(float(row["aic"]) for row in rows)
    denominator = math.fsum(math.exp(-(float(row["aic"])-minimum)/2) for row in rows)
    return [{"model": row["model"], "aic": float(row["aic"]),
             "delta_aic": float(row["aic"])-minimum,
             "akaike_weight": math.exp(-(float(row["aic"])-minimum)/2)/denominator,
             "log_akaike_weight": -(float(row["aic"])-minimum)/2-math.log(denominator),
             "source_column_c": float(row["reported_negative_log_likelihood"])}
            for row in rows]


def settings():
    return {row["key"]: row["value"] for row in load_csv("analysis_design.csv")}


def lundberg_results():
    rows, groups = read_lundberg()
    alpha_row = next(row for row in rows if row["quantity"] == "per_division_alpha")
    alpha = float(alpha_row["value"])
    result = []
    for (population, condition), group in groups.items():
        if condition != "untreated" or "scored_divisions_n" not in group:
            continue
        n = int(group["scored_divisions_n"]["value"])
        counts = compatible_counts(n, group["rejected_divisions_percent"]["value"])
        minimum = min(counts)
        result.append({"population": population, "n": n, "compatible_counts": counts,
                       "rejection_fraction": minimum/n,
                       "wilson_interval": [wilson(n, minimum)[0], wilson(n, max(counts))[1]],
                       "conditional_p": binomial_tail(n, minimum, alpha),
                       "markov_p_bound": min(1.0, n*alpha/minimum)})
    adjusted = holm([row["conditional_p"] for row in result])
    for row, value in zip(result, adjusted):
        row["holm_p"] = value
    total_n = sum(row["n"] for row in result)
    total_count = sum(min(row["compatible_counts"]) for row in result)
    return {"lines": result, "combined_n": total_n,
            "combined_count": total_count,
            "combined_conditional_p": binomial_tail(total_n, total_count, alpha),
            "combined_markov_p_bound": min(1.0, total_n*alpha/total_count)}


def printed_unit(number):
    """Sensitivity unit at the last nonzero decimal in the published scalar.

    A rounding box is a numerical sensitivity scenario, not author-provided
    estimation uncertainty. Trailing display zeroes add no claimed precision.
    """
    digits = number.partition(".")[2].rstrip("0")
    return F(1, 10**len(digits))


def r_rounding_box(row):
    p, rho = F(row["p"]), F(row["rho"])
    dp, dr = printed_unit(row["p"])/2, printed_unit(row["rho"])/2
    candidates = [effective_r(pp, rr) for pp in (p-dp, p+dp)
                  for rr in (rho-dr, rho+dr)]
    if p-dp <= F(1, 2) <= p+dp:
        candidates.extend([rho-dr, rho+dr])
    return min(candidates), max(candidates)


def relative_shift(interval, baseline):
    """Outward relative shift c/c_binomial-1 from two enclosing intervals."""
    return [interval["lower"]/baseline["upper"]-1,
            interval["upper"]/baseline["lower"]-1]


def outward(value, digits, upper=False):
    """Decimal display enclosing an exact endpoint; never nearest rounding."""
    number = value if isinstance(value, F) else F(str(value))
    scale = 10**digits
    scaled = number*scale
    integer = (-((-scaled.numerator)//scaled.denominator) if upper
               else scaled.numerator//scaled.denominator)
    sign = "-" if integer < 0 else ""
    absolute = abs(integer)
    return f"{sign}{absolute//scale}.{absolute%scale:0{digits}d}"


def display_enclosure(interval, digits=6):
    lower = interval.get("lower_exact", str(interval["lower"]))
    upper = interval.get("upper_exact", str(interval["upper"]))
    return f"[{outward(F(lower), digits)}, {outward(F(upper), digits, True)}]"


def display_relative_enclosure(interval, baseline, digits=2):
    lower = F(interval["lower_exact"])/F(baseline["upper_exact"])-1
    upper = F(interval["upper_exact"])/F(baseline["lower_exact"])-1
    return f"[{outward(100*lower, digits)}, {outward(100*upper, digits, True)}]"


def binomial_quantile(n, probability, target):
    cumulative = 0.0
    for k in range(n+1):
        cumulative += math.comb(n, k)*probability**k*(1-probability)**(n-k)
        if cumulative >= target:
            return k
    return n


def experiment_prediction(rows, baseline):
    design = settings()
    selected = next(row for row in rows if row["id"] == design["prediction_source_id"])
    k = int(design["prediction_maternal_k"])
    cutoff = int(design["prediction_minority_cutoff"])
    n_pairs = int(design["prediction_independent_pairs"])
    alpha = float(design["prediction_test_alpha"])
    verifier = verifier_module()

    def probability(p, rho):
        nums, den = verifier.integer_kernel(k, F(p), F(rho))
        return F(sum(num for j, num in enumerate(nums) if min(j, 2*k-j) <= cutoff), den)

    q0 = probability(design["binomial_p"], design["binomial_rho"])
    q1 = probability(selected["p"], selected["rho"])
    critical = next(x for x in range(n_pairs+1)
                    if binomial_tail(n_pairs, x, float(q0)) <= alpha)
    predictive_alpha = float(design["nominal_alpha"])
    band = [binomial_quantile(n_pairs, float(q1), predictive_alpha/2),
            binomial_quantile(n_pairs, float(q1), 1-predictive_alpha/2)]
    R = effective_r(F(selected["p"]), F(selected["rho"]))
    return {"source_id": selected["id"], "maternal_k": k, "total": 2*k,
            "minority_cutoff": cutoff, "independent_pairs": n_pairs,
            "binomial_extreme_probability": float(q0),
            "fitted_extreme_probability": float(q1),
            "binomial_expected_count": n_pairs*float(q0),
            "fitted_expected_count": n_pairs*float(q1),
            "critical_count": critical,
            "null_rejection_probability": binomial_tail(n_pairs, critical, float(q0)),
            "conditional_power": binomial_tail(n_pairs, critical, float(q1)),
            "conditional_prediction_interval": band,
            "lower_falsification_probability": 1-binomial_tail(n_pairs, band[0], float(q1)),
            "upper_falsification_probability": binomial_tail(n_pairs, band[1]+1, float(q1)),
            "conditional_fano": float((1-R)/2+R*k),
            "conditional_variance": float(k*((1-R)/2+R*k)),
            "R": float(R)}


@lru_cache(maxsize=1)
def results():
    """All derived report numbers, recomputed from CSVs; no saved result input."""
    verifier = verifier_module()
    design = settings()
    baseline = verifier.threshold(design["binomial_p"], design["binomial_rho"],
                                  width=design["threshold_width"])
    rows = load_csv("shibata_fits.csv")
    fits, controls = [], []
    for row in rows:
        R = effective_r(F(row["p"]), F(row["rho"]))
        record = {"id": row["id"], "population": row["population"],
                  "condition": row["condition"], "scope": row["scope"],
                  "p": float(row["p"]), "rho": float(row["rho"]), "R": float(R)}
        if row["scope"] == "HSR_control":
            controls.append(record)
            continue
        point = verifier.threshold(row["p"], row["rho"], width=design["threshold_width"])
        rounding = verifier.rounding_threshold(row["p"], row["rho"],
                                              printed_unit(row["p"]), printed_unit(row["rho"]),
                                              design["rounding_numerical_step"])
        info = conditional_design_interval(row, alpha=float(design["nominal_alpha"]))
        # Outward rounding of a *chosen sensitivity band*, not statistical calibration.
        step = int(design["sensitivity_rounding_denominator"])
        rlo = F(math.floor(info["r_interval"][0]*step), step)
        rhi = F(math.ceil(info["r_interval"][1]*step), step)
        propagated = verifier.moment_threshold_bounds(rlo, rhi)
        record.update(point_threshold=point, rounding_threshold=rounding,
                      rounding_R=list(map(float, r_rounding_box(row))),
                      relative_shift=relative_shift(rounding, baseline),
                      local_information=info, outward_sensitivity_R=[float(rlo), float(rhi)],
                      propagated_sensitivity=propagated,
                      propagated_relative_shift=relative_shift(propagated, baseline))
        fits.append(record)
    primary = [row for row in fits if row["scope"] == "ecDNA_baseline"]
    values = [row["R"] for row in primary]
    full_bound = verifier.moment_threshold_bounds(F(0), F(1))
    return {"baseline": baseline, "fits": fits, "controls": controls,
            "primary_summary": {"n_fits": len(primary), "R_min": min(values),
                                "R_max": max(values), "R_median": median(values),
                                "R_mean_descriptive": mean(values),
                                "smallest_downward_percent": -100*max(row["relative_shift"][1] for row in primary),
                                "largest_downward_percent": -100*min(row["relative_shift"][0] for row in primary)},
            "unidentified_R_threshold": full_bound,
            "unidentified_R_relative_shift": relative_shift(full_bound, baseline),
            "source_models": source_model_comparison(),
            "lundberg": lundberg_results(),
            "prediction": experiment_prediction(rows, baseline),
            "inventory": {"literature_records": len(load_csv("literature_audit.csv")),
                          "other_scalar_rows": len(load_csv("other_aggregates.csv")),
                          "lundberg_scalar_rows": len(load_csv("lundberg.csv")),
                          "shibata_fit_rows": len(rows),
                          "shibata_model_rows": len(load_csv("shibata_models.csv"))}}


def render_report():
    """Render every numerical result table directly from CSV-derived results."""
    r = results()
    inventory = r["inventory"]
    base = r["baseline"]
    summary = r["primary_summary"]
    lines = [
        "### Recomputed evidence inventory and calibration",
        "",
        f"The audit contains **{inventory['literature_records']} literature records** (including explicitly labelled duplicates), "
        f"**{inventory['lundberg_scalar_rows']} Lundberg scalar rows**, "
        f"**{inventory['other_scalar_rows']} other-study scalar rows**, "
        f"**{inventory['shibata_fit_rows']} Shibata fitted-parameter rows**, and "
        f"**{inventory['shibata_model_rows']} Shibata model-comparison rows**. "
        "Analysis-design rows are original choices, not published observations.",
        "",
        "Point statistics and probabilities are rounded decimal displays. "
        "Rigorous enclosure endpoints are displayed with outward rounding; JSON output also supplies exact rational endpoints.",
        "",
        f"The independently recomputed binomial threshold is enclosed by "
        f"**[{base['lower']:.8f}, {base['upper']:.8f}]**. "
        "These are exact-certificate endpoints, not a statistical interval.",
        "",
        f"Across the **{summary['n_fits']} ecDNA baseline clone fits**, "
        f"R ranges from **{summary['R_min']:.6f} to {summary['R_max']:.6f}**, "
        f"with descriptive median **{summary['R_median']:.6f}**. "
        f"Their conditional model thresholds are **{summary['smallest_downward_percent']:.2f}%–"
        f"{summary['largest_downward_percent']:.2f}% lower** than binomial. "
        "This is heterogeneity of published point fits in one selected cell-line system, "
        "not a population distribution of true R or a demonstrated treatment effect.",
        "",
        "The threshold intervals below include numerical error and a declared half-last-decimal "
        "rounding box for both source parameters. **They exclude unknown sampling and model error.** "
        "The full asymmetric beta-binomial kernel is used in every row.",
        "",
        "| Source fit | p | rho | Effective R | c* rounding/numerical enclosure | Relative shift (%) |",
        "| --- | ---: | ---: | ---: | --- | --- |",
    ]
    for row in r["fits"]:
        name = row["id"].replace("shibata_", "")
        if row["scope"] != "ecDNA_baseline":
            name += " ("+row["condition"]+")"
        interval = row["rounding_threshold"]
        shift = row["relative_shift"]
        lines.append(f"| {name} | {row['p']:.8f} | {row['rho']:.8f} | {row['R']:.6f} | "
                     f"{display_enclosure(interval)} | "
                     f"{display_relative_enclosure(interval, base, 4)} |")
    lines += [
        "",
        "The perturbation arms are separate conditions from the same study; "
        "their drug names identify source experiments and do not calibrate a drug dose to c. "
        "No averaging across these kernels is interpreted as a mixed-population threshold.",
        "",
        "### Uncertainty that can and cannot be propagated",
        "",
        "The next table is a **post-freeze local-information sensitivity analysis, not empirical "
        "confidence intervals**. The R bands use the nominal normal multiplier from the design CSV "
        "and the normalized hypothetical observation model described below. Source sampling units, "
        "filtering, likelihood normalization, and finite-sample regularity are unresolved. "
        "The mathematical c envelope is uniform over every fixed allocation-mixture law with R "
        "inside the displayed band, using Theorem 1. It does not grant the input band statistical coverage.",
        "",
        "| Source fit | Assumed independent divisions | Local R sensitivity band | Uniform c envelope | Relative shift envelope (%) |",
        "| --- | ---: | --- | --- | --- |",
    ]
    for row in r["fits"]:
        info = row["local_information"]
        band = row["outward_sensitivity_R"]
        envelope = row["propagated_sensitivity"]
        shift = row["propagated_relative_shift"]
        lines.append(f"| {row['id'].replace('shibata_', '')} | {info['assumed_independent_pairs']:.0f} | "
                     f"[{band[0]:.6f}, {band[1]:.6f}] | "
                     f"{display_enclosure(envelope)} | "
                     f"{display_relative_enclosure(envelope, base)} |")
    diagnostic = next(row for row in r["fits"] if row["id"] == "shibata_clone10")
    source = next(row for row in load_csv("shibata_fits.csv") if row["id"] == diagnostic["id"])
    b_est = (2*float(source["p"])-1)**2
    full_shift = r["unidentified_R_relative_shift"]
    lines += [
        "",
        f"A diagnostic of failed local regularity: clone10 has squared bias B={b_est:.6f}, "
        f"but its local B standard error is {diagnostic['local_information']['bias2_se']:.2f}, "
        "larger than the entire admissible domain. This is why even a positive R band is **not "
        "used to reject R=0 empirically**. Every propagated envelope above contains the binomial reference; "
        "the direction of the shift is unresolved by these conservative bounds.",
        "",
        "Without additional validated sampling/measurement assumptions, the aggregates supply no "
        "nontrivial calibrated confidence set for true R. Taking its full admissible range "
        f"[0, 1] gives c in [0, 0.5] and a relative-shift enclosure "
        f"**{display_relative_enclosure(r['unidentified_R_threshold'], base)}%** under the branching assumptions. "
        "This is an outer bound, not a sharp identified set or a confidence interval; "
        "its positive upper endpoint does not prove that an increase is attainable.",
        "",
        "### Published model comparisons and heterogeneity",
        "",
        "All values below are from the source aggregate comparison across the full clone set, "
        "including HSR controls. AIC is a source-model comparison, not a p value or a probability "
        "that a biological mechanism is true. The source's negative-valued column headed "
        "'Negative log-likelihood' behaves as a log-likelihood in its AIC calculation; "
        "the CSV preserves the signed values and the label discrepancy.",
        "",
        "| Model | Parameters | Allocation model | Source AIC | Recomputed delta AIC |",
        "| --- | ---: | --- | ---: | ---: |",
    ]
    source_models = {row["model"]: row for row in load_csv("shibata_models.csv")}
    for row in r["source_models"]:
        original = source_models[row["model"]]
        description = original["kernel"]+"; p "+original["p_structure"]+"; rho "+original["rho_structure"]
        lines.append(f"| {row['model']} | {original['n_parameters']} | {description} | "
                     f"{row['aic']:.3f} | {row['delta_aic']:.3f} |")
    lines += [
        "",
        "The favored source model has clone-specific allocation means and common rho. "
        "The separate clone fits therefore do **not** establish heterogeneous rho. "
        "No compatible second study supplies a new likelihood for the same estimand, "
        "so a cross-study hierarchical dispersion estimate is not fit.",
        "",
        "HSR controls are retained in the audit but excluded from ecDNA calibration:",
        "",
        "| HSR control | Source rho | Transformed R (descriptive only) |",
        "| --- | ---: | ---: |",
    ]
    for row in r["controls"]:
        lines.append(f"| {row['id'].replace('shibata_', '')} | {row['rho']:.8f} | {row['R']:.6f} |")
    lines += [
        "",
        "Nonzero fitted dispersion in these chromosomal controls warns that the fitted parameter "
        "can absorb counting effects or model mismatch; it is not by itself a measurement of physical ecDNA clustering.",
        "",
        "### Binomial-baseline challenges and opposing summaries",
        "",
        "Lundberg's integer counts are derived by enumerating all counts compatible with "
        "the printed whole-percent rounding and published denominators. The exact tails and "
        "Holm adjustment assume calibrated per-division tests and independent sampled mothers; "
        "Wilson intervals have the same sampling qualification. The Markov bound dispenses "
        "with independence but still needs calibrated marginal tests.",
        "",
        "| Line | Divisions | Compatible rejected counts | Conditional upper-tail p | Holm p | Wilson rejection-fraction interval | Dependence-robust Markov p bound |",
        "| --- | ---: | --- | ---: | ---: | --- | ---: |",
    ]
    for row in r["lundberg"]["lines"]:
        wi = row["wilson_interval"]
        lines.append(f"| {row['population']} | {row['n']} | {row['compatible_counts']} | "
                     f"{row['conditional_p']:.8g} | {row['holm_p']:.8g} | "
                     f"[{outward(wi[0], 6)}, {outward(wi[1], 6, True)}] | {row['markov_p_bound']:.6f} |")
    lu = r["lundberg"]
    lines += [
        "",
        f"The secondary pooled count is {lu['combined_count']}/{lu['combined_n']}, "
        f"with conditional p={lu['combined_conditional_p']:.8g} and "
        f"Markov bound {lu['combined_markov_p_bound']:.6f}. "
        "It is not an additional independent experiment or a fourth Holm-family member. "
        "This preserves PR63's verdict: conditional inconsistency, not an unconditional "
        "biological rejection and not an estimate of rho.",
        "",
        "Lange's complete retained Table SI 1 is opposing evidence:",
        "",
        "| Cell line / amplicon | Source KS distance | Source KS p | Source #samples |",
        "| --- | ---: | ---: | ---: |",
    ]
    lange = {}
    for row in load_csv("other_aggregates.csv"):
        if row["study"] == "Lange et al. 2022":
            lange.setdefault(row["population"], {})[row["quantity"]] = row["value"]
    for population, values in lange.items():
        lines.append(f"| {population} | {values['source_ks_d']} | {values['source_ks_p']} | "
                     f"{values['source_sample_size']} |")
    lines += [
        "",
        "These source KS tests do not reject their simulated distribution. Complementary sister "
        "dependence, variable maternal totals and reported main-figure/supplement differences prevent "
        "recovering a calibrated rho confidence interval from them. Failure to reject does not "
        "establish zero overdispersion. No study or source panel is dropped for having an inconvenient sign.",
        "",
        "### Single sharp future prediction",
        "",
    ]
    pr = r["prediction"]
    lines += [
        f"Use untreated **Clone12**, verified maternal pre-replication count **k={pr['maternal_k']}**, "
        f"and conserved post-replication total **N={pr['total']}**. In **{pr['independent_pairs']} independent "
        f"mother divisions**, count H, the number with a minority daughter receiving **at most "
        f"{pr['minority_cutoff']} copies**. Count molecular copies and true zeros, not unresolved hubs. "
        "Eligibility is set from the pre-division mother; nuclear loss is recorded separately and "
        "invalidates the conservation model rather than being silently excluded.",
        "",
        "| Quantity | Unbiased binomial | Source-fit transported prediction |",
        "| --- | ---: | ---: |",
        f"| Probability of an extreme split | {pr['binomial_extreme_probability']:.9f} | {pr['fitted_extreme_probability']:.9f} |",
        f"| Expected H | {pr['binomial_expected_count']:.6f} | {pr['fitted_expected_count']:.6f} |",
        f"| Conditional central 95% predictive interval for H | Not used for the rejection rule | {pr['conditional_prediction_interval']} |",
        "",
        f"Reject the unbiased-binomial experiment model when **H >= {pr['critical_count']}**. "
        f"Its exact null size is **{pr['null_rejection_probability']:.9f}**, and its power under "
        f"the fixed published-fit kernel is **{pr['conditional_power']:.9f}**. "
        f"Conversely, H < {pr['conditional_prediction_interval'][0]} or H > {pr['conditional_prediction_interval'][1]} "
        "falsifies this transported point-kernel prediction at the stated two-tail predictive rule "
        f"(tail probabilities {pr['lower_falsification_probability']:.9f} and "
        f"{pr['upper_falsification_probability']:.9f}). Both models can fail.",
        "",
        f"The same fixed-kernel prediction has conditional daughter Fano **{pr['conditional_fano']:.9f}** "
        f"and variance **{pr['conditional_variance']:.9f}**. "
        "These are alternative summaries of this one experiment, not additional independent tests. "
        "The prediction interval and power condition on the published kernel; they do not include "
        "unknown parameter, assay, clone-transfer or culture-dependence uncertainty.",
    ]
    return "\n".join(lines)+"\n"


def json_safe(value):
    if isinstance(value, F):
        return float(value)
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(v) for v in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", action="store_true", help="conditional model information only")
    parser.add_argument("--json", action="store_true", help="all recomputed results as JSON")
    parser.add_argument("--markdown", action="store_true", help="render numerical report tables")
    parser.add_argument("--check-readme", action="store_true", help="check numerical README block against recomputation")
    args = parser.parse_args()
    rows = load_csv("shibata_fits.csv")
    if args.check_readme:
        content = (HERE / "README.md").read_text(encoding="utf-8")
        found = content.split("<!-- BEGIN RECOMPUTED RESULTS -->\n", 1)[1].split("<!-- END RECOMPUTED RESULTS -->", 1)[0]
        if found != render_report()+"\n":
            raise SystemExit("README numerical block differs from recomputation")
        print("README numerical block matches CSV recomputation; nonclinical, unreviewed.")
    elif args.markdown:
        print(render_report(), end="")
    elif args.design:
        results = []
        for row in rows:
            if row["scope"] == "HSR_control":
                continue
            results.append({"id": row["id"], "R": effective_r(float(row["p"]), float(row["rho"])),
                            **conditional_design_interval(row)})
        print(json.dumps(results, indent=2))
    else:
        print(json.dumps(json_safe(results()), indent=2))


if __name__ == "__main__":
    main()
