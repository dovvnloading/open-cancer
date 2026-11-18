"""Independent event simulation of a synthetic containment comparison.

NONCLINICAL COMPUTATIONAL MATHEMATICS. This is not a treatment model or advice.

This file deliberately imports neither the deterministic verifier nor its
transition matrices or closed forms. It implements the occupancy events
directly. Each trial starts both arms at (sensitive, resistant) = (2, 0).

The containment arm instantaneously fills/culls sensitive occupancy while a
sensitive unit remains. It never recreates sensitivity after sensitive
extinction. The aggressive arm kills sensitive units at rate k per unit.
Both arms share resistant birth rates a at total occupancy one and b at total
occupancy two, resistant death rate d, and conversion rate mu per sensitive
unit. Progression is the first unregulated total occupancy of three.

A common independent exponential censor of rate delta is coupled to both
arms. Events of each biological type use a shared rate equal to the minimum
of the two arm rates and an arm-specific residual rate. Each marginal is
therefore the stipulated continuous-time Markov chain. Trial pairs are
independent; within-pair correlation reduces the uncertainty of differences.

Only the event winning each exponential race is sampled: holding times are
unneeded for the Bernoulli endpoint "progression before censor". If that
probability is h, E[min(T, censor)] = (1 - h) / delta. The returned gap is
containment's censored mean minus the aggressive arm's censored mean; a
positive gap favors containment for this mathematical endpoint only.
"""

from __future__ import annotations

import argparse
import json
import math
import random


def simulate_pair(
    mu: float,
    *,
    samples: int = 200_000,
    seed: int = 20261117,
    a: float = 100.0,
    b: float = 1.0,
    d: float = 3.0,
    k: float = 2.0,
    delta: float = 1.0,
) -> dict[str, int | float]:
    """Return paired Monte Carlo estimates and empirical standard errors.

    ``a`` and ``b`` are per-resistant-unit birth rates at total occupancy one
    and two; ``d`` is resistant death, ``k`` aggressive sensitive killing,
    and ``delta`` independent exponential censoring. ``mu`` is a rate, not
    a probability. Standard errors are descriptive Monte Carlo uncertainty,
    not a substitute for the exact deterministic sign certificate.
    """
    if isinstance(samples, bool) or not isinstance(samples, int) or samples < 2:
        raise ValueError("samples must be an integer at least two")
    mu, a, b, d, k, delta = map(float, (mu, a, b, d, k, delta))
    if not all(math.isfinite(value) for value in (mu, a, b, d, k, delta)):
        raise ValueError("all rates must be finite")
    if mu < 0 or d < 0 or k < 0 or b <= 0 or a < b or delta <= 0:
        raise ValueError("require mu,d,k >= 0, a >= b > 0, and delta > 0")

    uniform = random.Random(seed).random
    hits_c = hits_m = discordant = events = 0
    for _ in range(samples):
        sc = sm = 2
        rc = rm = 0
        active_c = active_m = True
        hit_c = hit_m = False

        while active_c or active_m:
            if active_c:
                mutation_c = mu * sc
                birth_c = (a if sc + rc == 1 else b) * rc
                death_c = d * rc
            else:
                mutation_c = birth_c = death_c = 0.0
            if active_m:
                mutation_m = mu * sm
                birth_m = (a if sm + rm == 1 else b) * rm
                death_m = d * rm
                kill_m = k * sm
            else:
                mutation_m = birth_m = death_m = kill_m = 0.0

            mutation_union = max(mutation_c, mutation_m)
            birth_union = max(birth_c, birth_m)
            death_union = max(death_c, death_m)
            total = delta + mutation_union + birth_union + death_union + kill_m
            draw = uniform() * total
            events += 1
            if draw < delta:
                break
            draw -= delta

            if draw < mutation_union:
                # This interval consists of shared and then residual clocks.
                shared = min(mutation_c, mutation_m)
                do_c = draw < shared or mutation_c > mutation_m
                do_m = draw < shared or mutation_m > mutation_c
                if do_c:
                    sc -= 1
                    rc += 1
                if do_m:
                    sm -= 1
                    rm += 1
            else:
                draw -= mutation_union
                if draw < birth_union:
                    shared = min(birth_c, birth_m)
                    do_c = draw < shared or birth_c > birth_m
                    do_m = draw < shared or birth_m > birth_c
                    if do_c:
                        rc += 1
                        if sc > 0:
                            sc -= 1
                        elif rc > 2:
                            hit_c = True
                    if do_m:
                        rm += 1
                        if sm + rm > 2:
                            hit_m = True
                else:
                    draw -= birth_union
                    if draw < death_union:
                        shared = min(death_c, death_m)
                        do_c = draw < shared or death_c > death_m
                        do_m = draw < shared or death_m > death_c
                        if do_c:
                            rc -= 1
                            # Refill only from sensitivity that still exists.
                            if sc > 0:
                                sc += 1
                        if do_m:
                            rm -= 1
                    else:
                        sm -= 1

            active_c = not hit_c and sc + rc > 0
            active_m = not hit_m and sm + rm > 0

        hits_c += hit_c
        hits_m += hit_m
        discordant += hit_c != hit_m

    probability_c = hits_c / samples
    probability_m = hits_m / samples
    probability_gap = probability_c - probability_m
    # Unbiased sample variances of Bernoulli outcomes and paired differences.
    se_c = math.sqrt(probability_c * (1.0 - probability_c) / (samples - 1))
    se_m = math.sqrt(probability_m * (1.0 - probability_m) / (samples - 1))
    variance_gap = (discordant / samples - probability_gap**2) / (samples - 1)
    se_gap = math.sqrt(max(0.0, variance_gap))
    return {
        "samples": samples,
        "seed": seed,
        "mean_c": (1.0 - probability_c) / delta,
        "mean_m": (1.0 - probability_m) / delta,
        "gap": -probability_gap / delta,
        "standard_error_c": se_c / delta,
        "standard_error_m": se_m / delta,
        "gap_standard_error": se_gap / delta,
        "probability_c": probability_c,
        "probability_m": probability_m,
        "probability_gap": probability_gap,
        "probability_gap_standard_error": se_gap,
        "discordant_pairs": discordant,
        "events": events,
    }


def main() -> None:
    """Run the three prespecified synthetic checks without writing files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=200_000)
    parser.add_argument("--seed", type=int, default=20261117)
    parser.add_argument("--mu", type=float, action="append")
    args = parser.parse_args()
    mutation_rates = args.mu if args.mu is not None else [0.25, 2.5, 25.0]
    results = []
    for index, mu in enumerate(mutation_rates):
        result = simulate_pair(mu, samples=args.samples, seed=args.seed + index)
        results.append({"mu": mu, **result})
    print(json.dumps({
        "label": "Nonclinical computational mathematics; synthetic rates; no treatment advice",
        "parameters": {"a": 100, "b": 1, "d": 3, "k": 2, "delta": 1},
        "gap_definition": "E[min(T_containment, censor)] - E[min(T_aggressive, censor)]",
        "results": results,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
