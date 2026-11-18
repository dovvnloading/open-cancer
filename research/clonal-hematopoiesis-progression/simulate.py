"""Seeded synthetic Gillespie checks of the marked branching model (stdlib).

No empirical parameters are estimated here. A stochastic check tests the
implementation against a separate algorithm; it does not establish a theorem
or validate a biological model.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import random

import model

ROOT = Path(__file__).resolve().parent


def one_run(n: int, birth: float, death: float, hit: float,
            horizon: float, rng: random.Random) -> bool:
    """First terminal mark before horizon, with exact exponential waiting times."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("n must be a nonnegative integer")
    if any(not math.isfinite(x) or x < 0 for x in (birth, death, hit, horizon)):
        raise ValueError("finite nonnegative rates and horizon required")
    total_per_cell = birth + death + hit
    elapsed = 0.0
    while n and total_per_cell:
        elapsed += rng.expovariate(n * total_per_cell)
        if elapsed > horizon:
            return False
        event = rng.random() * total_per_cell
        if event < hit:
            return True
        if event < hit + birth:
            n += 1
        else:
            n -= 1
    return False


def run_scenarios(root: Path = ROOT) -> list[dict]:
    with (root / "model_scenarios.csv").open(encoding="utf-8", newline="") as f:
        scenarios = list(csv.DictReader(f))
    results = []
    for row in scenarios:
        n, reps, seed = (int(row[k]) for k in ("n", "replicates", "seed"))
        birth, death, hit, horizon = (
            float(row[k]) for k in ("birth", "death", "hit", "horizon_years"))
        rng = random.Random(seed)
        events = sum(one_run(n, birth, death, hit, horizon, rng) for _ in range(reps))
        exact = model.progression_probability(n, birth, death, hit, horizon)
        se = math.sqrt(exact * (1 - exact) / reps)
        exact_h = model.cumulative_hazard(n, birth, death, hit, horizon)
        lower, upper = model.rare_event_bounds(n, birth, death, hit, horizon)
        results.append({
            "scenario_id": row["scenario_id"], "replicates": reps,
            "seed": seed, "events": events, "simulation_probability": events / reps,
            "exact_probability": exact, "model_binomial_se": se,
            "standardized_error": (events / reps - exact) / se if se else 0.0,
            "exact_cumulative_hazard": exact_h, "rare_bound_lower": lower,
            "rare_bound_upper": upper, "net_fitness": birth - death,
        })
    return results


if __name__ == "__main__":
    print(json.dumps({"status": "synthetic_nonclinical", "scenarios": run_scenarios()},
                     indent=2, sort_keys=True, allow_nan=False))
