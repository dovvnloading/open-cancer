#!/usr/bin/env python3
"""Independent synthetic checks for metastatic-seeding mathematics.

Nonclinical computational mathematics. No treatment advice or patient data.
Only the Python standard library is used. This module deliberately does not
import verifier.py: cell events, immigration, mutations, regions, reads, and
clone labels are generated directly. Agreement is evidence about implementation,
not a proof, a novelty certificate, or biological validation.
"""

import argparse
from fractions import Fraction
from itertools import combinations
import json
import math
import random


DEFAULT_SEED = 20261117


def poisson_waiting(rng, rate, duration=1.0):
    """Count literal exponential interarrival times, not an inverse Poisson CDF."""
    if rate < 0 or duration < 0:
        raise ValueError("rate and duration must be nonnegative")
    if rate == 0 or duration == 0:
        return 0
    count = 0
    elapsed = rng.expovariate(rate)
    while elapsed < duration:
        count += 1
        elapsed += rng.expovariate(rate)
    return count


def bd_trace(rng, birth, death, horizon, initial=1):
    """Gillespie linear birth-death process; retain extinction and every event.

    Returns (piecewise-constant intervals, final size). There is no population
    cap, rejection of extinct runs, or deterministic large-population shortcut.
    The default experiment has thousands of cells; large user-supplied rates
    can take correspondingly long to simulate.
    """
    if min(birth, death, horizon, initial) < 0 or int(initial) != initial:
        raise ValueError("invalid birth-death parameters")
    size = int(initial)
    time = 0.0
    trace = []
    while time < horizon:
        if size == 0 or birth + death == 0:
            trace.append((time, horizon, size))
            break
        next_time = time + rng.expovariate(size * (birth + death))
        if next_time >= horizon:
            trace.append((time, horizon, size))
            break
        trace.append((time, next_time, size))
        size += 1 if rng.random() < birth / (birth + death) else -1
        time = next_time
    return trace, size


def simulate_bd(rng, birth, death, horizon, initial=1):
    """Final population, using the same explicit event generator as bd_trace."""
    return bd_trace(rng, birth, death, horizon, initial)[1]


def seed_trace(rng, trace, seed_rate):
    """Cox immigration along a fixed primary path using exponential hazards."""
    if seed_rate < 0:
        raise ValueError("seed_rate must be nonnegative")
    if seed_rate == 0:
        return []
    times = []
    hazard_to_event = rng.expovariate(1.0)
    for start, stop, size in trace:
        intensity = seed_rate * size
        if intensity == 0:
            continue
        cursor = start
        available = intensity * (stop - cursor)
        while hazard_to_event < available:
            cursor += hazard_to_event / intensity
            times.append(cursor)
            hazard_to_event = rng.expovariate(1.0)
            available = intensity * (stop - cursor)
        hazard_to_event -= available
    return times


def alive_probability(birth, death, age):
    """Independent closed-form oracle for one-cell birth-death survival to age."""
    if age == 0:
        return 1.0
    if birth == death:
        return 1.0 / (1.0 + birth * age)
    growth = birth - death
    return growth / (birth - death * math.exp(-growth * age))


def _simpson(function, start, stop, panels=8):
    """Fixed composite Simpson rule; panels is positive and even."""
    step = (stop - start) / panels
    total = function(start) + function(stop)
    for i in range(1, panels):
        total += (4 if i % 2 else 2) * function(start + i * step)
    return total * step / 3


def poisson_cdf(mean, maximum):
    """Floating-point comparison oracle, independent of the random generator."""
    if maximum < 0:
        return 0.0
    if mean == 0:
        return 1.0
    term = math.exp(-mean)
    total = term
    for count in range(1, maximum + 1):
        term *= mean / count
        total += term
    return min(1.0, total)


def _binomial_cdf(trials, probability, maximum):
    return sum(math.comb(trials, k) * probability ** k
               * (1 - probability) ** (trials - k)
               for k in range(min(maximum, trials) + 1))


def draw_atom(rng, atoms):
    """Draw from (value, mass) pairs; masses must sum to one."""
    draw = rng.random()
    cumulative = 0.0
    for value, mass in atoms:
        cumulative += float(mass)
        if draw < cumulative:
            return value
    # Floating-point conversion can leave cumulative a few ulps below one.
    return atoms[-1][0]


def detection_trial(rng, prevalence, regions, depth, vaf):
    """A region occupancy Bernoulli followed by literal read Bernoullis."""
    for _ in range(regions):
        if rng.random() < prevalence:
            for _ in range(depth):
                if rng.random() < vaf:
                    return True
    return False


def finite_pool_trial(rng, pool_size, invisible, regions, depth, vaf):
    """Sample distinct regions without replacement and simulate their reads."""
    for region in rng.sample(range(pool_size), regions):
        if region >= invisible:
            for _ in range(depth):
                if rng.random() < vaf:
                    return True
    return False


def clone_trial(rng, probabilities, samples, depth, vaf):
    """Categorical cell identity/read draws; unallocated mass is background."""
    probabilities = tuple(float(p) for p in probabilities)
    if any(not math.isfinite(p) or p < 0 or p > 1 for p in probabilities):
        raise ValueError("clone probabilities must be finite and lie in [0,1]")
    total = math.fsum(probabilities)
    if total > 1:
        raise ValueError("clone probabilities must sum to at most one")
    if (not isinstance(samples, int) or isinstance(samples, bool) or samples < 0
            or not isinstance(depth, int) or isinstance(depth, bool) or depth < 0):
        raise ValueError("samples and depth must be nonnegative integers")
    if not math.isfinite(vaf) or not 0 <= vaf <= 1:
        raise ValueError("vaf must lie in [0,1]")
    atoms = list(enumerate(probabilities))
    if total < 1:
        atoms.append((None, 1 - total))
    found = set()
    for _ in range(samples):
        label = draw_atom(rng, atoms)
        if label is None:
            continue
        for _ in range(depth):
            if rng.random() < vaf:
                found.add(label)
                break
    return len(found)


def _record(theorem, name, observations, expected, variance):
    number = len(observations)
    mean = sum(observations) / number
    empirical_variance = sum((x - mean) ** 2 for x in observations) / number
    return {
        "theorem": theorem,
        "name": name,
        "observed": mean,
        "expected": float(expected),
        "trials": number,
        "variance": float(variance),
        "empirical_variance": empirical_variance,
        "standard_error": math.sqrt(float(variance) / number),
    }


def _bernoulli_record(theorem, name, observations, expected):
    expected = float(expected)
    return _record(theorem, name, observations, expected,
                   expected * (1 - expected))


def _atomic_miss(atoms, regions, detection):
    return sum(float(weight) * (1 - detection * float(value)) ** regions
               for value, weight in atoms)


def _finite_pool_miss(pool_size, invisible, regions, detection):
    """Hypergeometric/read oracle, independent of random.sample simulation."""
    visible = pool_size - invisible
    denominator = math.comb(pool_size, regions)
    value = Fraction(0)
    for detected_regions in range(regions + 1):
        if detected_regions <= visible and regions - detected_regions <= invisible:
            ways = (math.comb(visible, detected_regions)
                    * math.comb(invisible, regions - detected_regions))
            value += (Fraction(ways, denominator)
                      * (1 - detection) ** detected_regions)
    return value


def _three_atom_weights(support, first, factorial_second):
    """Solve three finite-pool moment equations by independent elimination."""
    rows = [
        [Fraction(1) for _ in support] + [Fraction(1)],
        [Fraction(x) for x in support] + [first],
        [Fraction(x * (x - 1)) for x in support] + [factorial_second],
    ]
    for column in range(3):
        pivot = next(i for i in range(column, 3) if rows[i][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        factor = rows[column][column]
        rows[column] = [x / factor for x in rows[column]]
        for i in range(3):
            if i != column:
                factor = rows[i][column]
                rows[i] = [x - factor * y
                           for x, y in zip(rows[i], rows[column])]
    return [rows[i][-1] for i in range(3)]


def finite_pool_extrema(pool_size, first, factorial_second, regions, detection):
    """Enumerate all three-atom candidates for the small *oracle* fixture.

    The event simulator accepts any pool size; this enumeration checks a small
    exact example and is not the claimed proof for arbitrary pool sizes.
    """
    candidates = []
    for support in combinations(range(pool_size + 1), 3):
        weights = _three_atom_weights(support, first, factorial_second)
        if min(weights) >= 0:
            atoms = [(x, weight) for x, weight in zip(support, weights) if weight]
            objective = sum(weight * _finite_pool_miss(
                pool_size, x, regions, detection) for x, weight in atoms)
            candidates.append((objective, atoms))
    if not candidates:
        raise ValueError("infeasible finite-pool moments")
    return min(candidates, key=lambda item: item[0]), max(
        candidates, key=lambda item: item[0])


def _clone_recovery(probabilities, samples, detection):
    probability = 0.0
    for count in range(len(probabilities) + 1):
        for selected in combinations(probabilities, count):
            probability += (-1) ** count * (1 - detection * sum(selected)) ** samples
    return probability


def run_simulations(seed=DEFAULT_SEED, replicates=4000):
    """Return fixed-seed, theorem-labeled synthetic observations and oracles.

    Each check provides its known sampling variance, permitting predeclared
    standard-error tolerances without fitting a tolerance to the realized run.
    Birth-death/Cox checks use max(200, replicates//4) complete realizations;
    other checks use replicates. No trajectories are discarded.
    """
    if not isinstance(replicates, int) or replicates < 1:
        raise ValueError("replicates must be a positive integer")
    checks = []
    report = {
        "label": "Synthetic nonclinical computational mathematics; no treatment advice",
        "seed": seed,
        "replicates": replicates,
        "checks": checks,
    }

    # T1: real primary Gillespie trace, explicit Cox arrivals and independent
    # metastatic trajectories. The fixed path is conditioned on, not averaged.
    rng = random.Random(seed + 101)
    horizon = 6.0
    trace, final_size = bd_trace(rng, 0.95, 0.25, horizon, initial=24)
    seed_rate = 0.0015
    metastasis_birth, metastasis_death = 0.65, 0.35
    event_mean = seed_rate * sum(size * (stop - start)
                                 for start, stop, size in trace)
    alive_mean = seed_rate * sum(
        size * _simpson(lambda t: alive_probability(
            metastasis_birth, metastasis_death, horizon - t), start, stop)
        for start, stop, size in trace)
    refined_alive_mean = seed_rate * sum(
        size * _simpson(lambda t: alive_probability(
            metastasis_birth, metastasis_death, horizon - t), start, stop, panels=16)
        for start, stop, size in trace)
    counts, alive_counts = [], []
    process_replicates = max(200, replicates // 4)
    for _ in range(process_replicates):
        times = seed_trace(rng, trace, seed_rate)
        counts.append(len(times))
        alive_counts.append(sum(simulate_bd(
            rng, metastasis_birth, metastasis_death, horizon - time) > 0
                                for time in times))
    checks.append(_record("T1", "cox_event_mean", counts, event_mean, event_mean))
    checks.append(_record("T1", "alive_metastasis_mean", alive_counts,
                          alive_mean, alive_mean))
    checks.append(_record("T1", "alive_metastasis_second_moment",
                          [x * x for x in alive_counts], alive_mean ** 2 + alive_mean,
                          4 * alive_mean ** 3 + 6 * alive_mean ** 2 + alive_mean))
    checks.append(_bernoulli_record("T1", "zero_alive_metastases",
                                   [x == 0 for x in alive_counts], math.exp(-alive_mean)))
    extinction_probability = 1 - alive_probability(0.4, 0.7, 3.0)
    extinct = [simulate_bd(rng, 0.4, 0.7, 3.0) == 0 for _ in range(replicates)]
    checks.append(_bernoulli_record("T1", "extinction_retained", extinct,
                                   extinction_probability))
    report["primary_path"] = {
        "initial_cells": 24, "birth": 0.95, "death": 0.25,
        "horizon": horizon, "final_cells": final_size,
        "population_events": len(trace) - 1, "maximum_cells": max(x[2] for x in trace),
        "seed_rate": seed_rate, "immigration_mean": event_mean,
        "survival_thinned_mean_simpson": alive_mean,
        "survival_mean_quadrature_refinement_difference": abs(alive_mean - refined_alive_mean),
        "extinct_runs_discarded": 0, "population_cap": None,
    }
    shared, private, alternative_private = [], [], []
    for _ in range(replicates):
        shared.append(poisson_waiting(rng, 2.0, 2.0))
        # Last shared ancestor a=2; seed s=4 or s=7; observation B=8.
        private.append(poisson_waiting(rng, 2.0, 2.0)
                       + poisson_waiting(rng, 2.0, 4.0))
        alternative_private.append(poisson_waiting(rng, 2.0, 5.0)
                                   + poisson_waiting(rng, 2.0, 1.0))
    checks.append(_record("T1", "shared_passenger_mean", shared, 4, 4))
    checks.append(_record("T1", "private_before_after_seed_mean", private, 12, 12))
    checks.append(_record("T1", "private_alternative_seed_mean",
                          alternative_private, 12, 12))
    checks.append(_record("T1", "shared_private_covariance",
                          [(x - 4) * (y - 12) for x, y in zip(shared, private)], 0, 48))

    # T2: fixed boundary clock laws, and an explicit shared-law counterexample.
    rng = random.Random(seed + 202)
    early_mean, late_mean, threshold = 35.0, 20.0, 27
    early_error = poisson_cdf(early_mean, threshold - 1)
    late_error = 1 - poisson_cdf(late_mean, threshold - 1)
    early_counts = [poisson_waiting(rng, 1.0, early_mean) for _ in range(replicates)]
    late_counts = [poisson_waiting(rng, 1.0, late_mean) for _ in range(replicates)]
    checks.append(_bernoulli_record("T2", "early_boundary_error",
                                   [x < threshold for x in early_counts], early_error))
    checks.append(_bernoulli_record("T2", "late_boundary_error",
                                   [x >= threshold for x in late_counts], late_error))
    early_coupling_rng = random.Random(seed + 203)
    late_coupling_rng = random.Random(seed + 203)
    early_s, early_lag, late_s, late_lag = 4, 0, 6, 2
    early_a, late_a = early_s - early_lag, late_s - late_lag
    confounded_counts = [poisson_waiting(early_coupling_rng, 1.0, 10 - early_a)
                         for _ in range(replicates)]
    coupled_late_counts = [poisson_waiting(late_coupling_rng, 1.0, 10 - late_a)
                           for _ in range(replicates)]
    checks.append(_record("T2", "confounded_common_clock", confounded_counts, 6, 6))
    report["timing_fixture"] = {
        "exposure": 5, "clock_rate": 1, "B": 10, "c": 5, "h": 2, "ell": 1,
        "threshold": threshold, "early_mean": early_mean, "late_mean": late_mean,
        "confounding": {"early_s": early_s, "early_lag": early_lag, "late_s": late_s,
                        "late_lag": late_lag, "common_a": early_a, "common_mean": 6,
                        "coupled_observations_equal": confounded_counts == coupled_late_counts,
                        "note": "Same a permits the same Poisson path in both histories."},
    }

    # T3: the two sharp atomic laws and a four-atom law with the SAME moments.
    rng = random.Random(seed + 303)
    upper_atoms = [(Fraction(0), Fraction(1, 5)), (Fraction(5, 8), Fraction(4, 5))]
    lower_atoms = [(Fraction(3, 8), Fraction(4, 5)), (Fraction(1), Fraction(1, 5))]
    interior_atoms = [(Fraction(0), Fraction(5, 64)),
                      (Fraction(1, 3), Fraction(27, 64)),
                      (Fraction(2, 3), Fraction(27, 64)),
                      (Fraction(1), Fraction(5, 64))]
    regions, depth, vaf = 8, 5, 0.2
    detection = 1 - (1 - vaf) ** depth
    for name, atoms in [("upper_extremizer", upper_atoms),
                        ("lower_extremizer", lower_atoms),
                        ("four_atom_interior", interior_atoms)]:
        misses = [not detection_trial(rng, float(draw_atom(rng, atoms)),
                                      regions, depth, vaf) for _ in range(replicates)]
        checks.append(_bernoulli_record("T3", name, misses,
                                       _atomic_miss(atoms, regions, detection)))
    linked_misses = []
    for _ in range(replicates):
        occupied = rng.random() < 0.5
        observed = False
        if occupied:
            for _ in range(regions * depth):
                if rng.random() < vaf:
                    observed = True
                    break
        linked_misses.append(not observed)
    linked_expected = 0.5 + 0.5 * (1 - vaf) ** (regions * depth)
    checks.append(_bernoulli_record("T3", "linked_region_counterexample",
                                   linked_misses, linked_expected))
    report["sampling_fixture"] = {
        "first_moment": "1/2", "second_moment": "5/16", "regions": regions,
        "depth": depth, "vaf": vaf, "detection": detection,
        "lower": _atomic_miss(lower_atoms, regions, detection),
        "upper": _atomic_miss(upper_atoms, regions, detection),
        "linked_region_expected": linked_expected,
        "iid_expected_at_fixed_F_half": (1 - detection / 2) ** regions,
        "linked_counterexample_violates": "conditional independence of regions",
    }

    # T4: a literal witness/read stage and clock stage, with deliberately wrong
    # decisions after witness failure. Independence gives an exact sub-bound
    # expectation here; the theorem's union bound does not require independence.
    rng = random.Random(seed + 404)
    witness_miss = _atomic_miss(upper_atoms, regions, detection)
    for name, clock_mean, clock_error, is_early in [
            ("joint_early_error", early_mean, early_error, True),
            ("joint_late_error", late_mean, late_error, False)]:
        mistakes = []
        for _ in range(replicates):
            captured = detection_trial(rng, float(draw_atom(rng, upper_atoms)),
                                       regions, depth, vaf)
            count = poisson_waiting(rng, 1.0, clock_mean)
            clock_mistake = count < threshold if is_early else count >= threshold
            mistakes.append(not captured or clock_mistake)
        expected = witness_miss + (1 - witness_miss) * clock_error
        checks.append(_bernoulli_record("T4", name, mistakes, expected))
    report["joint_bound"] = {
        "witness_upper": witness_miss,
        "worst_clock_error": max(early_error, late_error),
        "union_bound": min(1.0, witness_miss + max(early_error, late_error)),
    }
    # Independent cohorts estimate one hit and a pair of hits. Draw a NEW F
    # per pilot individual, retaining the SAME F across that individual's pair.
    pilot_size, pilot_epsilon, pilot_detection = 50, 0.15, 0.75
    pilot_failures = []
    for _ in range(replicates):
        one_hits, pair_hits = 0, 0
        for _ in range(pilot_size):
            prevalence = float(draw_atom(rng, interior_atoms))
            one_hits += detection_trial(rng, prevalence, 1, 2, 0.5)
            prevalence = float(draw_atom(rng, interior_atoms))
            first_hit = detection_trial(rng, prevalence, 1, 2, 0.5)
            second_hit = detection_trial(rng, prevalence, 1, 2, 0.5)
            pair_hits += first_hit and second_hit
        m_lower = max(0.0, (one_hits / pilot_size - pilot_epsilon) / pilot_detection)
        s_upper = min(1.0, (pair_hits / pilot_size + pilot_epsilon) / pilot_detection ** 2)
        pilot_failures.append(m_lower > 0.5 or s_upper < 5 / 16)
    # The chosen rational fixture has nonintegral cutoffs: failure is H1>=27
    # or H2<=1. Independent cohorts make the union probability explicit.
    first_failure = 1 - _binomial_cdf(pilot_size, pilot_detection / 2, 26)
    second_failure = _binomial_cdf(pilot_size, pilot_detection ** 2 * 5 / 16, 1)
    pilot_failure = first_failure + second_failure - first_failure * second_failure
    checks.append(_bernoulli_record("T4", "finite_pilot_rectangle_failure",
                                   pilot_failures, pilot_failure))
    report["pilot_fixture"] = {
        "cohort_sizes": [pilot_size, pilot_size], "epsilon": pilot_epsilon,
        "pilot_detection": pilot_detection,
        "exact_failure_via_binomial_oracle": pilot_failure,
        "hoeffding_union_upper": min(1.0, 2 * math.exp(-2 * pilot_size * pilot_epsilon ** 2)),
        "observed_failure": sum(pilot_failures) / replicates,
    }

    # T5: finite archived pool. Here U is the INTEGER number of invisible regions.
    rng = random.Random(seed + 505)
    pool_size, pool_regions = 6, 3
    pool_depth, pool_vaf = 2, 0.5
    pool_detection = Fraction(3, 4)
    pool_atoms = [(0, Fraction(1, 4)), (2, Fraction(1, 2)), (5, Fraction(1, 4))]
    lower_pool, upper_pool = finite_pool_extrema(
        pool_size, Fraction(9, 4), Fraction(6), pool_regions, pool_detection)
    for name, atoms in [("finite_pool_lower", lower_pool[1]),
                        ("finite_pool_upper", upper_pool[1]),
                        ("finite_pool_interior", pool_atoms)]:
        misses = [not finite_pool_trial(rng, pool_size, draw_atom(rng, atoms),
                                        pool_regions, pool_depth, pool_vaf)
                  for _ in range(replicates)]
        expected = sum(weight * _finite_pool_miss(
            pool_size, value, pool_regions, pool_detection) for value, weight in atoms)
        checks.append(_bernoulli_record("T5", name, misses, expected))
    report["finite_pool_fixture"] = {
        "pool_size": pool_size, "sampled_regions": pool_regions,
        "pilot_one_region_miss": "3/8", "pilot_two_region_miss": "1/5",
        "future_detection": str(pool_detection),
        "lower": str(lower_pool[0]), "upper": str(upper_pool[0]),
        "lower_atoms": [[x, str(w)] for x, w in lower_pool[1]],
        "upper_atoms": [[x, str(w)] for x, w in upper_pool[1]],
    }

    # T6: sample cells and read their clone labels; extinct or identical seeders
    # do not acquire magically distinguishable labels in this simulation.
    rng = random.Random(seed + 606)
    probabilities, samples = [0.5, 0.3, 0.2], 12
    clone_depth, clone_vaf = 2, 0.5
    clone_detection = 1 - (1 - clone_vaf) ** clone_depth
    recovered = [clone_trial(rng, probabilities, samples, clone_depth, clone_vaf) == 3
                 for _ in range(replicates)]
    exact_recovery = _clone_recovery(probabilities, samples, clone_detection)
    checks.append(_bernoulli_record("T6", "all_distinct_clone_labels_recovered",
                                   recovered, exact_recovery))
    report["clone_fixture"] = {
        "probabilities": probabilities, "cell_samples": samples,
        "read_detection": clone_detection, "exact_recovery": exact_recovery,
        "union_miss_bound": sum((1 - clone_detection * p) ** samples
                                for p in probabilities),
        "indistinguishable_event_counterexample": {
            "event_counts": [1, 2], "endpoint_genotype_sets": [["A"], ["A"]],
            "observation_signature_equal": True,
            "note": "Endpoint labels do not identify repeated identical or extinct seeders.",
            "simulation_is_proof": False,
        },
    }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--replicates", type=int, default=4000)
    arguments = parser.parse_args()
    print(json.dumps(run_simulations(arguments.seed, arguments.replicates),
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
