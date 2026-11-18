#!/usr/bin/env python3
"""Independent, fixed-seed simulations for nonclinical surveillance mathematics.

All individuals and assays generated here are synthetic.  This is unreviewed
computational research, not treatment or monitoring advice for any patient.
No theorem evaluator or verifier is imported.  Growth paths, individual sampled
molecules, first positive draws, and progression events are simulated directly.
Only Python's standard library is used; the CLI writes JSON to stdout only.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
import random
from typing import Callable, Iterable, Mapping


class _Moments:
    """Online sample moments; the variance is descriptive, not a confidence bound."""

    def __init__(self) -> None:
        self.n = 0
        self.mean = 0.0
        self.m2 = 0.0

    def add(self, value: float) -> None:
        self.n += 1
        delta = value - self.mean
        self.mean += delta / self.n
        self.m2 += delta * (value - self.mean)

    @property
    def variance(self) -> float:
        return max(0.0, self.m2 / (self.n - 1)) if self.n > 1 else 0.0


def _finite(value: float, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be a finite real number")
    return number


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _schedules(
    schedules: Mapping[str, Iterable[float]],
) -> dict[str, tuple[float, ...]]:
    if not isinstance(schedules, Mapping) or not schedules:
        raise ValueError("schedules must be a nonempty mapping")
    result = {}
    for name, times in schedules.items():
        if not isinstance(name, str) or not name:
            raise ValueError("schedule names must be nonempty strings")
        try:
            sequence = tuple(_finite(t, "test time") for t in times)
        except TypeError as exc:
            raise ValueError("test times must be iterable") from exc
        if any(t < 0 for t in sequence):
            raise ValueError("test times must be nonnegative")
        if any(left >= right for left, right in zip(sequence, sequence[1:])):
            raise ValueError("each schedule must have strictly increasing times")
        result[name] = sequence
    return result


class _Summary:
    def __init__(self, names: Iterable[str]) -> None:
        self.names = tuple(names)
        self.successes = dict.fromkeys(self.names, 0)
        self.leads = {name: _Moments() for name in self.names}
        self.losses = {name: _Moments() for name in self.names}
        self.pairs = {
            (left, right): (_Moments(), _Moments(), _Moments())
            for left in self.names
            for right in self.names
            if left != right
        }
        self.scored = 0
        self.draws = 0

    def add(self, values: dict[str, tuple[int, float, float]], scored: bool) -> None:
        self.draws += 1
        self.scored += int(scored)
        for name, (success, lead, loss) in values.items():
            self.successes[name] += success
            self.leads[name].add(lead)
            self.losses[name].add(loss)
        for (left, right), moments in self.pairs.items():
            for index, item in enumerate(moments):
                item.add(values[left][index] - values[right][index])

    def result(self) -> dict:
        result = {}
        for name in self.names:
            paired = {}
            for other in self.names:
                if name == other:
                    continue
                probability, lead, loss = self.pairs[(name, other)]
                paired[other] = {
                    "probability_difference": probability.mean,
                    "probability_difference_variance": probability.variance,
                    "mean_lead_difference": lead.mean,
                    "lead_difference_variance": lead.variance,
                    "mean_loss_difference": loss.mean,
                    "loss_difference_variance": loss.variance,
                }
            result[name] = {
                "draws": self.draws,
                "scored": self.scored,
                "successes": self.successes[name],
                "probability": self.successes[name] / self.draws,
                "mean_lead": self.leads[name].mean,
                "mean_loss": self.losses[name].mean,
                "lead_variance": self.leads[name].variance,
                "loss_variance": self.losses[name].variance,
                "paired": paired,
            }
        return result


def _path_values(
    schedules: dict[str, tuple[float, ...]],
    positive: set[float],
    onset: float,
    progression: float,
) -> dict[str, tuple[int, float, float]]:
    opportunity = progression - onset
    result = {}
    for name, times in schedules.items():
        first = next((t for t in times if t in positive), None)
        lead = 0.0 if first is None else progression - first
        # min/max only remove last-bit subtraction error; hits occur after onset.
        lead = min(opportunity, max(0.0, lead))
        result[name] = (int(first is not None), lead, max(0.0, opportunity - lead))
    return result


def simulate_growth(
    schedules: dict[str, tuple[float, ...]],
    draws: int,
    seed: int,
    rate_sampler: Callable[[random.Random], float],
    *,
    a: float,
    b: float,
    endpoint: float | None = None,
    burden_sampler: Callable[[random.Random], float] | None = None,
    molecule_count: int | None = None,
    lod: float = 0.01,
) -> dict:
    """Simulate exponential paths B(t)=exp(Z+R*t) and scheduled blood draws.

    ``a`` and ``b`` are log burdens for ideal detectability and progression,
    with a < b.  Z defaults to zero and must be less than a, so the initial
    burden is below ideal detectability and onset is (a-Z)/R.  Each sampled R
    must be finite and strictly positive.  Patient rate/burden draws use one
    RNG; molecule draws use a separate seeded RNG, so assay noise cannot alter
    the sampled patient distribution when schedules or assay settings change.

    With ``molecule_count=None``, eligible blood draws detect perfectly.  With
    m molecules, the true fraction is logistic(logit(lod)+Z+R*t-a).  The draw
    must meet the true fraction gate f>=lod, and K>=ceil(m*lod), where K is
    generated by m independent Bernoulli draws, not an analytic tail formula.
    The cutoff uses lod's decimal-string value to avoid floating ceil artifacts.
    A shared assay outcome is drawn once per time for all compared schedules;
    outcomes at distinct times are independent conditional on the patient.

    A hit must precede progression strictly.  If endpoint is supplied, only
    patients with progression<=endpoint score; others contribute zero to all
    outcomes.  Every returned mean/probability has denominator ``draws``, not
    ``scored`` or detected cases.  Loss equals ideal opportunity minus achieved
    lead.  Each result includes empirical variances and paired differences;
    these alone are not distribution-free confidence intervals for unbounded
    lead times.  Invalid sampler outputs raise rather than being resampled.
    """
    normalized = _schedules(schedules)
    draws = _positive_integer(draws, "draws")
    a, b = _finite(a, "a"), _finite(b, "b")
    if a >= b:
        raise ValueError("a must be less than b")
    if not callable(rate_sampler):
        raise ValueError("rate_sampler must be callable")
    if burden_sampler is not None and not callable(burden_sampler):
        raise ValueError("burden_sampler must be callable or None")
    if endpoint is not None:
        endpoint = _finite(endpoint, "endpoint")
        if endpoint <= 0:
            raise ValueError("endpoint must be positive")
    lod = _finite(lod, "lod")
    if not 0 < lod < 1:
        raise ValueError("lod must be strictly between zero and one")
    threshold = None
    if molecule_count is not None:
        molecule_count = _positive_integer(molecule_count, "molecule_count")
        threshold = math.ceil(Fraction(str(lod)) * molecule_count)
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")
    patient_rng = random.Random(seed)
    assay_rng = random.Random(seed ^ 0xC7D0A5A7)
    all_times = sorted({t for times in normalized.values() for t in times})
    logit_lod = math.log(lod) - math.log1p(-lod)
    summary = _Summary(normalized)
    zero_values = {name: (0, 0.0, 0.0) for name in normalized}
    for _ in range(draws):
        rate = _finite(rate_sampler(patient_rng), "sampled rate")
        if rate <= 0:
            raise ValueError("sampled rates must be strictly positive")
        z = 0.0 if burden_sampler is None else _finite(
            burden_sampler(patient_rng), "sampled log burden"
        )
        if z >= a:
            raise ValueError("initial log burden must be below ideal detectability")
        onset = (a - z) / rate
        progression = (b - z) / rate
        if not math.isfinite(onset) or not math.isfinite(progression):
            raise ValueError("sampled rate/burden produces nonfinite event times")
        if endpoint is not None and progression > endpoint:
            summary.add(zero_values, False)
            continue
        positive = set()
        for time in all_times:
            if time < onset or time >= progression:
                continue
            if molecule_count is None:
                positive.add(time)
                continue
            log_odds = logit_lod + z + rate * time - a
            if log_odds >= 0:
                fraction = 1.0 / (1.0 + math.exp(-log_odds))
            else:
                exp_odds = math.exp(log_odds)
                fraction = exp_odds / (1.0 + exp_odds)
            observed = sum(
                assay_rng.random() < fraction for _ in range(molecule_count)
            )
            if observed >= threshold:
                positive.add(time)
        summary.add(_path_values(normalized, positive, onset, progression), True)
    return summary.result()


def sample_uniform_onsets(
    schedules: dict[str, tuple[float, ...]],
    draws: int,
    seed: int,
    H: float,
    W: float,
) -> dict:
    """Directly sample onset U~Uniform[0,H] and progression U+W.

    Tests detect perfectly between onset and progression, with a strict upper
    endpoint.  All test times must be in [0,H].  No extra terminal test is added:
    callers must include H if their mathematical schedule has a test there.
    Returned moments and paired results have the simulate_growth layout.
    """
    normalized = _schedules(schedules)
    draws = _positive_integer(draws, "draws")
    H, W = _finite(H, "H"), _finite(W, "W")
    if H <= 0 or W <= 0:
        raise ValueError("H and W must be positive")
    if any(t > H for times in normalized.values() for t in times):
        raise ValueError("uniform-onset test times must lie in [0,H]")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")
    rng = random.Random(seed)
    summary = _Summary(normalized)
    all_times = sorted({t for times in normalized.values() for t in times})
    for _ in range(draws):
        onset = rng.random() * H
        progression = onset + W
        if not math.isfinite(progression):
            raise ValueError("H and W produce nonfinite progression times")
        positive = {t for t in all_times if onset <= t < progression}
        summary.add(_path_values(normalized, positive, onset, progression), True)
    return summary.result()


def thinning_simulation(
    leads: Iterable[float], step: float, draws: int, seed: int,
) -> dict:
    """Sample synthetic leads and randomly retain alternate positive visits.

    The supplied finite list defines a uniform empirical distribution; it must
    contain synthetic or otherwise permitted data.  Original first positivity
    is time zero, progression is time L, and the next positive visit is at
    ``step``.  A random parity retains time zero or time step.  Missed cases
    receive zero lead, keeping the original cohort denominator.  This routine
    never reconditions on detection after thinning.
    """
    values = tuple(_finite(x, "lead") for x in leads)
    if not values or any(x < 0 for x in values):
        raise ValueError("leads must be a nonempty sequence of nonnegative values")
    step = _finite(step, "step")
    if step <= 0:
        raise ValueError("step must be positive")
    draws = _positive_integer(draws, "draws")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")
    rng = random.Random(seed)
    original, thinned, reduction = _Moments(), _Moments(), _Moments()
    successes = 0
    for _ in range(draws):
        clinical_time = values[rng.randrange(len(values))]
        first_retained_visit = step if rng.getrandbits(1) else 0.0
        if first_retained_visit < clinical_time:
            successes += 1
            lead = clinical_time - first_retained_visit
        else:
            lead = 0.0
        original.add(clinical_time)
        thinned.add(lead)
        reduction.add(clinical_time - lead)
    return {
        "draws": draws,
        "successes": successes,
        "probability": successes / draws,
        "original_mean": math.fsum(values) / len(values),
        "sampled_original_mean": original.mean,
        "thinned_mean": thinned.mean,
        "mean_reduction": reduction.mean,
        "lead_variance": thinned.variance,
        "reduction_variance": reduction.variance,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=20261117)
    args = parser.parse_args()
    schedules = {
        "equal_synthetic": (1.0, 2.0, 3.0, 4.0),
        "geometric_synthetic": (0.5, 1.0, 2.0, 4.0),
    }
    outputs = {
        "label": "NONCLINICAL SYNTHETIC COMPUTATIONAL MATHEMATICS; UNREVIEWED",
        "lognormal_perfect": simulate_growth(
            schedules, args.draws, args.seed,
            lambda rng: rng.lognormvariate(0.0, 0.7), a=1.0, b=2.0,
            endpoint=4.0,
        ),
        "gamma_noisy": simulate_growth(
            schedules, args.draws, args.seed,
            lambda rng: rng.gammavariate(2.0, 0.5), a=1.0, b=2.0,
            endpoint=4.0, burden_sampler=lambda rng: rng.uniform(-0.2, 0.2),
            molecule_count=1000, lod=0.01,
        ),
        "uniform_onsets": sample_uniform_onsets(
            schedules, args.draws, args.seed, H=4.0, W=8.0,
        ),
        "thinning_synthetic": thinning_simulation(
            (0.8, 2.0, 5.0, 10.0, 16.5), 3.0, args.draws, args.seed,
        ),
    }
    print(json.dumps(outputs, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
