"""Independent Gillespie check for a synthetic, nonclinical branching model.

This module deliberately does not import the deterministic verifier. It samples
the birth, death, and mutant-offspring events directly. It is computational
mathematics, has no patient calibration, and supplies no treatment advice.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
import random
import time


@dataclass(frozen=True)
class SimulationResult:
    """Uncensored counts of extinction by the end of a finite schedule."""

    extinct: int
    trials: int

    @property
    def probability(self) -> float:
        return self.extinct / self.trials


def simulate(
    parameters: dict[str, float],
    schedule: list[tuple[float, float]],
    trials: int,
    seed: int,
    initial: tuple[int, ...] = (1, 0, 0, 0),
    *,
    max_events_per_trial: int = 1_000_000,
) -> SimulationResult:
    """Sample exact continuous-time events, using binary floating-point times.

    Types are S, A-resistant, B-resistant, doubly resistant. Each schedule pair
    is ``(duration, fraction_A)``; the fraction of B is ``1 - fraction_A``.
    Clonal births add one cell. Mutant births retain the parent and add a mutant
    daughter. A rate-cap or numerical failure raises; it is never counted as a
    surviving or extinct trial. Reaching a phase boundary redraws the waiting
    time with the new hazards, as required by exponential memorylessness.
    """
    if isinstance(trials, bool) or not isinstance(trials, int) or trials <= 0:
        raise ValueError("trials must be a positive integer")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")
    if (
        isinstance(max_events_per_trial, bool)
        or not isinstance(max_events_per_trial, int)
        or max_events_per_trial <= 0
    ):
        raise ValueError("max_events_per_trial must be a positive integer")
    if len(initial) != 4 or any(
        isinstance(n, bool) or not isinstance(n, int) or n < 0 for n in initial
    ):
        raise ValueError("initial must contain four nonnegative integer counts")
    keys = ("b", "d", "k", "c", "mu", "p", "w", "nu")
    if set(parameters) != set(keys):
        raise ValueError("parameters must contain exactly " + ", ".join(keys))
    b, d, k, c, mu, p, w, nu = (float(parameters[key]) for key in keys)
    if not all(math.isfinite(x) for x in (b, d, k, c, mu, p, w, nu)):
        raise ValueError("parameters must be finite")
    if min(b, d, k, mu, nu) < 0 or c < 1 or not (0 <= p <= 1 and 0 <= w <= 1):
        raise ValueError("invalid nonnegative rate, sensitivity, or probability")
    if not schedule:
        raise ValueError("schedule must contain at least one phase")

    # Each event is (daughter_type, rate); daughter_type == -1 is a death.
    # These channels express the stochastic population model directly.
    offspring = (
        ((0, b), (1, mu * (1 - p) * w), (2, mu * (1 - p) * (1 - w)), (3, mu * p)),
        ((1, b), (3, nu)),
        ((2, b), (3, nu)),
        ((3, b),),
    )
    phases = []
    for duration, fraction_a in schedule:
        duration, fraction_a = float(duration), float(fraction_a)
        if not math.isfinite(duration) or duration <= 0:
            raise ValueError("phase durations must be finite and positive")
        if not math.isfinite(fraction_a) or not 0 <= fraction_a <= 1:
            raise ValueError("drug fractions must lie in [0, 1]")
        deaths = (d + k, d + c * k * (1 - fraction_a), d + k * fraction_a, d)
        channels = tuple(
            ((-1, deaths[i]),) + offspring[i] for i in range(4)
        )
        hazards = tuple(sum(rate for _, rate in group) for group in channels)
        if not all(math.isfinite(rate) for rate in hazards):
            raise ValueError("rates overflow floating-point arithmetic")
        phases.append((duration, channels, hazards))

    rng = random.Random(seed)
    extinct = 0
    for _ in range(trials):
        counts = list(initial)
        events = 0
        for duration, channels, hazards in phases:
            if not any(counts):
                break
            elapsed = 0.0
            while elapsed < duration:
                weights = tuple(counts[i] * hazards[i] for i in range(4))
                total = sum(weights)
                if not math.isfinite(total):
                    raise RuntimeError("population hazard overflow; no result reported")
                if total == 0:
                    break
                waiting = rng.expovariate(total)
                if waiting >= duration - elapsed:
                    break
                elapsed += waiting
                events += 1
                if events > max_events_per_trial:
                    raise RuntimeError("event cap exceeded; no censored result reported")

                # A single uniform variate chooses a cell type and its event.
                selector = rng.random() * total
                parent = None
                for i, weight in enumerate(weights):
                    if selector < weight:
                        parent = i
                        break
                    selector -= weight
                if parent is None or counts[parent] == 0:
                    raise RuntimeError("event selection lost floating-point resolution")
                selector /= counts[parent]
                daughter = None
                for candidate, rate in channels[parent]:
                    if selector < rate:
                        daughter = candidate
                        break
                    selector -= rate
                if daughter is None:
                    raise RuntimeError("event selection lost floating-point resolution")
                if daughter == -1:
                    counts[parent] -= 1
                else:
                    counts[daughter] += 1
                if not any(counts):
                    break
        extinct += not any(counts)
    return SimulationResult(extinct, trials)


def main() -> None:
    """Run the deliberately synthetic finite-horizon comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=150_000)
    parser.add_argument("--seed", type=int, default=20261117)
    args = parser.parse_args()
    parameters = {"b": 0.05, "d": 2, "k": 1, "c": 8, "mu": 6, "p": 0.1, "w": 0.5, "nu": 0.02}
    schedules = {
        "combination": [(0.3, 0.5)],
        "AB": [(0.15, 1), (0.15, 0)],
        "AB2": [(0.075, 1), (0.075, 0)] * 2,
    }
    output = {"label": "Synthetic nonclinical computational mathematics; no treatment advice", "results": {}}
    started = time.monotonic()
    for offset, (name, schedule) in enumerate(schedules.items()):
        result = simulate(parameters, schedule, args.trials, args.seed + offset)
        output["results"][name] = {
            "extinct": result.extinct,
            "trials": result.trials,
            "probability": result.probability,
            "seed": args.seed + offset,
        }
    output["elapsed_seconds"] = time.monotonic() - started
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
