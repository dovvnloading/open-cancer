#!/usr/bin/env python3
"""Independent, synthetic simulation of copy-selective ecDNA branching.

Nonclinical computational mathematics; no treatment advice.  No patient data.
This module deliberately does not import the exact-arithmetic verifier.  It
simulates the offspring experiment itself, including complementary sisters.
Generation numbers are genealogical, not elapsed clinical or laboratory time.

Floating-point Monte Carlo is a falsification check, not a mathematical proof.
Division probabilities are evaluated in floating point in the default cases;
extremely small u*s**k can underflow to zero.  Extreme-parameter guarantees
belong to the exact-arithmetic verifier and proofs, not this simulator.
Standard errors use independent replicate populations as the sampling units;
related daughter cells are never treated as independent replicates.
"""

import argparse
import json
import math
import random


DEFAULT_SEED = 2718281828
METRICS = (
    "population", "zero", "positive", "copy_sum", "copy_square_sum",
    "extinct_positive",
)


def _check_inputs(initial, generations, u, s, max_cells):
    initial = tuple(initial)
    if any(isinstance(k, bool) or not isinstance(k, int) or k < 0
           for k in initial):
        raise ValueError("initial copy numbers must be nonnegative integers")
    if isinstance(generations, bool) or not isinstance(generations, int) \
            or generations < 0:
        raise ValueError("generations must be a nonnegative integer")
    if not (0.0 <= u <= 1.0 and 0.0 <= s <= 1.0):
        raise ValueError("u and s must be probabilities in [0, 1]")
    if isinstance(max_cells, bool) or not isinstance(max_cells, int) \
            or max_cells < 1:
        raise ValueError("max_cells must be a positive integer")
    return initial


def simulate_one(initial, generations, u, s, rng, keep_zeros=True,
                 max_cells=1_000_000, *, division_probability=None):
    """Simulate one independent finite population without a type cutoff.

    A k-copy cell divides with probability u*s**k, otherwise it dies.  A
    dividing cell produces J and 2*k-J, where J is a fair-binomial draw made
    with 2*k random bits.  This preserves the exact sister-cell dependence.

    The returned dictionary has ``positive_copies``, ``zero_count`` and
    ``zero_lineages_tracked``.  When keep_zeros=False, the positive process is
    unchanged but zero cells and their descendants are discarded; population
    then means the tracked positive population only.

    At u=1, zero cells double deterministically and are stored as an integer.
    Otherwise each zero cell's Bernoulli experiment is simulated explicitly.
    max_cells limits positive cells held in memory and stochastic zero-cell
    experiments per generation.  An exceeded limit raises RuntimeError; it
    never silently censors a replicate or declares it extinct.

    An optional division_probability(k) callback replaces u*s**k.  It must
    be a pure, time-independent function; probabilities are checked and
    cached by type within a population, including the zero type.  The u and
    s arguments are ignored by this custom law but must remain valid inputs.
    """
    initial = _check_inputs(initial, generations, u, s, max_cells)
    if division_probability is not None and not callable(division_probability):
        raise ValueError("division_probability must be callable or None")
    probability_cache = {}

    def probability(k):
        if division_probability is None:
            return u * s ** k
        if k not in probability_cache:
            value = division_probability(k)
            try:
                valid = math.isfinite(value) and 0.0 <= value <= 1.0
            except (TypeError, ValueError, OverflowError):
                valid = False
            if not valid:
                raise ValueError(
                    f"division_probability({k}) must be finite and in [0, 1]"
                )
            probability_cache[k] = value
        return probability_cache[k]

    p_zero = probability(0)
    positive = [k for k in initial if k]
    zeros = initial.count(0) if keep_zeros else 0
    for generation in range(generations):
        if len(positive) > max_cells \
                or (0.0 < p_zero < 1.0 and zeros > max_cells):
            raise RuntimeError(
                "Simulation resource cap exceeded before generation "
                f"{generation + 1}; no censored result is returned"
            )
        if keep_zeros:
            if p_zero == 1.0:
                new_zeros = 2 * zeros
            elif p_zero == 0.0:
                new_zeros = 0
            else:
                new_zeros = 2 * sum(rng.random() < p_zero for _ in range(zeros))
        else:
            new_zeros = 0
        children = []
        for k in positive:
            if rng.random() >= probability(k):
                continue
            j = rng.getrandbits(2 * k).bit_count()
            for daughter in (j, 2 * k - j):
                if daughter:
                    children.append(daughter)
                elif keep_zeros:
                    new_zeros += 1
            if len(children) > max_cells:
                raise RuntimeError(
                    "Simulation resource cap exceeded during generation "
                    f"{generation + 1}; no censored result is returned"
                )
        positive, zeros = children, new_zeros
    return {
        "positive_copies": tuple(positive),
        "zero_count": zeros,
        "zero_lineages_tracked": bool(keep_zeros),
        "custom_division_law": division_probability is not None,
    }


def outcome_metrics(outcome, weight_base=None):
    """Return per-population observables, optionally the h(k)=k*a**(k-1) sum."""
    positive = outcome["positive_copies"]
    zeros = outcome["zero_count"]
    row = {
        "population": len(positive) + zeros,
        "zero": zeros,
        "positive": len(positive),
        "copy_sum": sum(positive),
        "copy_square_sum": sum(k * k for k in positive),
        "extinct_positive": int(not positive),
    }
    if weight_base is not None:
        if not (0.0 < weight_base <= 1.0):
            raise ValueError("weight_base must be in (0, 1]")
        row["h_sum"] = math.fsum(
            k * weight_base ** (k - 1) for k in positive
        )
    return row


def _ratio_se(numerator, denominator, means, covariances, trials):
    denom = means[denominator]
    if not denom:
        return None
    ratio = means[numerator] / denom
    variance = (covariances[numerator][numerator]
                - 2.0 * ratio * covariances[numerator][denominator]
                + ratio * ratio * covariances[denominator][denominator])
    return math.sqrt(max(0.0, variance) / trials) / denom


def _pooled_summary(means, covariances, trials, denominator):
    """Ratio estimates and replicate-cluster delta-method standard errors."""
    count = means[denominator]
    if count == 0:
        return {"mean": None, "variance": None,
                "mean_standard_error": None, "variance_standard_error": None,
                "fano": None, "fano_standard_error": None}
    mean = means["copy_sum"] / count
    variance = means["copy_square_sum"] / count - mean * mean
    gradients = {
        denominator: (mean * mean - variance) / count,
        "copy_sum": -2.0 * mean / count,
        "copy_square_sum": 1.0 / count,
    }
    variance_of_variance = math.fsum(
        gx * gy * covariances[x][y]
        for x, gx in gradients.items() for y, gy in gradients.items()
    ) / trials
    mass = means["copy_sum"]
    if mass:
        fano = means["copy_square_sum"] / mass - mean
        fano_gradients = {
            denominator: mass / (count * count),
            "copy_sum": -means["copy_square_sum"] / (mass * mass) - 1.0 / count,
            "copy_square_sum": 1.0 / mass,
        }
        variance_of_fano = math.fsum(
            gx * gy * covariances[x][y]
            for x, gx in fano_gradients.items()
            for y, gy in fano_gradients.items()
        ) / trials
        fano_se = math.sqrt(max(0.0, variance_of_fano))
    else:
        fano, fano_se = None, None
    return {
        "mean": mean,
        "variance": max(0.0, variance),
        "mean_standard_error": _ratio_se(
            "copy_sum", denominator, means, covariances, trials
        ),
        "variance_standard_error": math.sqrt(max(0.0, variance_of_variance)),
        "fano": fano,
        "fano_standard_error": fano_se,
    }


def run_ensemble(initial, generations, u, s, trials=10_000,
                 seed=DEFAULT_SEED, keep_zeros=True, weight_base=None,
                 max_cells=1_000_000, *, division_probability=None):
    """Return Monte Carlo observables and uncertainty from independent wells.

    ``means`` are sample means of per-well sums, not normalized cell moments.
    ``pooled`` is computed by pooling every sampled cell across wells; its
    standard errors account for within-well dependence using the covariance
    of whole-well sums.  ``positive_pooled`` additionally gates on k>0.
    Monte Carlo estimates of extinction always refer to the specified finite
    generation; no nonextinction observation is labelled eternal survival.
    A custom division_probability callback has the same meaning and validation
    as in simulate_one; output u and s are null because the callback replaces
    those parameters, and custom_division_law is true.
    """
    initial = _check_inputs(initial, generations, u, s, max_cells)
    if isinstance(trials, bool) or not isinstance(trials, int) or trials < 2:
        raise ValueError("trials must be an integer at least 2")
    if weight_base is not None and not (0.0 < weight_base <= 1.0):
        raise ValueError("weight_base must be in (0, 1]")
    rng = random.Random(seed)
    names = METRICS + (("h_sum",) if weight_base is not None else ())
    sums = dict.fromkeys(names, 0)
    products = {x: dict.fromkeys(names, 0) for x in names}
    for _ in range(trials):
        outcome = simulate_one(
            initial, generations, u, s, rng, keep_zeros, max_cells,
            division_probability=division_probability,
        )
        row = outcome_metrics(outcome, weight_base)
        for x in names:
            sums[x] += row[x]
            for y in names:
                products[x][y] += row[x] * row[y]
    means = {x: sums[x] / trials for x in names}
    covariances = {
        # All non-h observables are integers.  Forming this numerator before
        # division preserves exact cancellation for deterministic copy mass.
        # The optional h statistic and final reported ratios remain floating.
        x: {y: (trials * products[x][y] - sums[x] * sums[y])
            / (trials * (trials - 1))
            for y in names}
        for x in names
    }
    sample_variances = {x: max(0.0, covariances[x][x]) for x in names}
    standard_errors = {
        x: math.sqrt(sample_variances[x] / trials) for x in names
    }
    pooled = _pooled_summary(means, covariances, trials, "population")
    pooled["zero_fraction"] = (
        means["zero"] / means["population"] if means["population"] else None
    )
    pooled["zero_fraction_standard_error"] = _ratio_se(
        "zero", "population", means, covariances, trials
    )
    return {
        "label": "synthetic nonclinical computational mathematics",
        "seed": seed,
        "trials": trials,
        "generations": generations,
        "initial": list(initial),
        "u": u if division_probability is None else None,
        "s": s if division_probability is None else None,
        "custom_division_law": division_probability is not None,
        "zero_lineages_tracked": bool(keep_zeros),
        "weight_base": weight_base,
        "means": means,
        "sample_variances": sample_variances,
        "standard_errors": standard_errors,
        "covariances": covariances,
        "pooled": pooled,
        "positive_pooled": _pooled_summary(
            means, covariances, trials, "positive"
        ),
        "finite_generation_positive_extinction_frequency":
            means["extinct_positive"],
    }


def fixed_seed_suite(trials=3_000):
    """Independent synthetic cases for external comparison with theorems.

    No theoretical answers are imported or used to generate trajectories.
    The threshold cases are finite-horizon illustrations only.  The optional
    h_sum observable lets the separate verifier test its eigenfunction law.
    """
    s = 0.8
    a = (1.0 - math.sqrt(1.0 - s)) / (1.0 + math.sqrt(1.0 - s))
    cases = {
        "finite_population": run_ensemble(
            [0, 1, 3], 4, 0.9, s, trials, DEFAULT_SEED,
            weight_base=a,
        ),
        "one_generation_selection": run_ensemble(
            [0, 1, 1, 2, 4, 8], 1, 0.9, s, trials, DEFAULT_SEED + 1,
        ),
        "therapy_variance_scaling": run_ensemble(
            [1, 3], 12, 1.0, 0.75, trials, DEFAULT_SEED + 2,
        ),
    }
    for index, survival in enumerate((0.7, 0.75, 0.8)):
        cases[f"threshold_s_{survival}"] = run_ensemble(
            [1], 25, 1.0, survival, trials, DEFAULT_SEED + 10 + index,
            keep_zeros=False,
        )
    for k in (1, 8):
        cases[f"monomorphic_fano_k_{k}"] = run_ensemble(
            [k], 4, 1.0, 0.9, trials, DEFAULT_SEED + 30 + k,
        )
    for index, ratio in enumerate((0.25, 0.5)):
        case = run_ensemble(
            [1], 30, 1.0, 1.0, trials, DEFAULT_SEED + 50 + index,
            keep_zeros=False,
            division_probability=lambda k, c=ratio: 1.0 / (1.0 + c * k),
        )
        case["linear_death_to_birth_ratio"] = ratio
        cases[f"linear_hazard_ratio_{ratio}"] = case
    return {
        "status": "synthetic model checks; not biological validation",
        "generation_warning": "genealogical generations, not elapsed time",
        "survival_warning": "finite-horizon survival does not prove persistence",
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=3_000)
    args = parser.parse_args()
    print(json.dumps(fixed_seed_suite(args.trials), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
