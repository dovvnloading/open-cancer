"""Independent stochastic checks of a nonclinical branching-process model.

All simulated observations are synthetic.  This module does not import the
verifier, use its recurrences, or sample its extinction probabilities.  Source
cells are followed through their actual exponential event clocks.  HLA-hidden
families are followed through their birth/death jump chain.  A population cap
is reported as unresolved; it is never relabelled eradication.

Complete loss of *modeled* presentation is an idealization, not a claim that
ordinary allele-specific HLA loss hides every neoantigen.  No treatment advice.
"""

import argparse
from dataclasses import dataclass
import json
import math
import random


@dataclass(frozen=True)
class Parameters:
    """Per-cell event rates; ``loss * k`` is the total single-antigen loss rate."""

    b: float = 1.0
    d: float = 0.5
    kappa: float = 1.0
    gain: float = 1.0
    loss: float = 0.2
    hla: float = 0.0
    b_hla: float = 1.0
    d_hla: float = 0.5

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0:
                raise ValueError("rates must be finite and nonnegative")
            object.__setattr__(self, name, value)


def _loads(initial_loads):
    values = tuple(initial_loads)
    if any(isinstance(k, bool) or not isinstance(k, int) or k < 0
           for k in values):
        raise ValueError("antigen loads must be nonnegative integers")
    return values


def _positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(name + " must be a positive integer")


def _source_rate(p, k):
    return p.b + p.d + p.kappa * k + p.gain + p.loss * k + p.hla


def _source_event(p, k, draw):
    """Select one event from a uniform draw on the total per-cell rate."""
    if draw < p.b:
        return "birth"
    draw -= p.b
    if draw < p.d + p.kappa * k:
        return "death"
    draw -= p.d + p.kappa * k
    if draw < p.gain:
        return "gain"
    draw -= p.gain
    if draw < p.loss * k:
        return "loss"
    return "hla"


def _hla_family(rng, p, initial, threshold, event_budget):
    """Simulate actual H-cell births/deaths; time cancels in its jump chain.

    Hitting ``threshold`` certifies persistence only up to the explicitly
    reported false-positive probability (d_hla / b_hla)**threshold when the
    H process is supercritical.  No analytic extinction probability is sampled.
    """
    count, used = initial, 0
    while count:
        total = p.b_hla + p.d_hla
        if total == 0:
            # A permanently living cell is persistence, with no classification
            # error.  This degenerate case is outside the main note's theorems.
            return "escape", used
        if count >= threshold:
            if p.b_hla > p.d_hla:
                return "escape", used
            return "unresolved", used
        if used >= event_budget:
            return "unresolved", used
        count += 1 if rng.random() * total < p.b_hla else -1
        used += 1
    return "extinct", used


def _escape_trial(rng, p, initial_loads, initial_hla, threshold, max_events):
    """Depth-first order is valid because distinct descendants never interact."""
    events = 0
    if initial_hla:
        outcome, used = _hla_family(rng, p, initial_hla, threshold, max_events)
        events += used
        if outcome != "extinct":
            return outcome, events
    pending = list(initial_loads)
    while pending:
        k = pending.pop()
        while True:
            if events >= max_events:
                return "unresolved", events
            total = _source_rate(p, k)
            if total == 0:
                return "unresolved", events
            # Independent clocks are simulated even though only the eventual
            # outcome, rather than calendar time, is used in this experiment.
            rng.expovariate(total)
            event = _source_event(p, k, rng.random() * total)
            events += 1
            if event == "birth":
                pending.append(k)
            elif event == "death":
                break
            elif event == "gain":
                k += 1
            elif event == "loss":
                k -= 1
            else:
                outcome, used = _hla_family(
                    rng, p, 1, threshold, max_events - events)
                events += used
                if outcome != "extinct":
                    return outcome, events
                break
    return "eradicated", events


def sample_escape(p, initial_loads, trials, seed, *, initial_hla=0,
                  hla_threshold=32, max_events=200000):
    """Follow source and H families until extinction, H threshold, or censoring.

    The returned frequency bracket assigns every unresolved trial either
    outcome.  It is a *sample* bracket, not a confidence interval.  Sampling
    uncertainty must be added separately.  ``hla_false_positive_bound`` bounds
    the probability of a false H persistence call per experiment.  Stopping at
    the first threshold crossing avoids a multiplicity factor for H families.
    """
    initial_loads = _loads(initial_loads)
    for value, name in ((trials, "trials"), (hla_threshold, "hla_threshold"),
                        (max_events, "max_events")):
        _positive_integer(value, name)
    if (isinstance(initial_hla, bool) or not isinstance(initial_hla, int)
            or initial_hla < 0):
        raise ValueError("initial_hla must be a nonnegative integer")
    rng = random.Random(seed)
    counts = {"eradicated": 0, "escape": 0, "unresolved": 0}
    events = 0
    for _ in range(trials):
        outcome, used = _escape_trial(
            rng, p, initial_loads, initial_hla, hla_threshold, max_events)
        counts[outcome] += 1
        events += used
    error = ((p.d_hla / p.b_hla) ** hla_threshold
             if p.b_hla > p.d_hla else 0.0)
    return dict(counts, trials=trials, seed=seed, events=events,
                hla_threshold=hla_threshold, hla_false_positive_bound=error,
                escape_frequency=counts["escape"] / trials,
                eradication_frequency=counts["eradicated"] / trials,
                escape_sample_bracket=(counts["escape"] / trials,
                    (counts["escape"] + counts["unresolved"]) / trials))


def _source_lifetime_trial(rng, p, initial_loads, max_events):
    pending = list(initial_loads)
    area, events, largest_type = 0.0, 0, max(initial_loads, default=0)
    while pending:
        k = pending.pop()
        while True:
            if events >= max_events:
                return None, area, events, largest_type
            total = _source_rate(p, k)
            if total == 0:
                return None, area, events, largest_type
            area += rng.expovariate(total)
            event = _source_event(p, k, rng.random() * total)
            events += 1
            if event == "birth":
                pending.append(k)
            elif event in ("death", "hla"):
                break
            elif event == "gain":
                k += 1
                largest_type = max(largest_type, k)
            else:
                k -= 1
    return area, area, events, largest_type


def source_lifetime_samples(p, initial_loads, trials, seed, *, max_events=200000):
    """Sample the total source cell-time, following every branching descendant.

    H conversion terminates source occupancy; H cells are not included in the
    occupation integral.  Use hla=0 for the baseline T_k in the rare-H theorem.
    Censored partial areas are reported separately and never used as full areas.
    """
    initial_loads = _loads(initial_loads)
    _positive_integer(trials, "trials")
    _positive_integer(max_events, "max_events")
    rng = random.Random(seed)
    areas, partial_areas = [], []
    events, largest_type = 0, max(initial_loads, default=0)
    for _ in range(trials):
        area, partial, used, largest = _source_lifetime_trial(
            rng, p, initial_loads, max_events)
        areas.append(area)
        partial_areas.append(partial)
        events += used
        largest_type = max(largest_type, largest)
    censored = sum(value is None for value in areas)
    return {"areas": areas, "partial_areas": partial_areas,
            "censored": censored, "trials": trials, "seed": seed,
            "events": events, "largest_type": largest_type,
            "mean": math.fsum(areas) / trials if not censored else None,
            "second_moment": (math.fsum(x * x for x in areas) / trials
                              if not censored else None)}


def _source_at_time_trial(rng, p, initial_loads, horizon, max_events, max_cells):
    counts = {}
    for k in initial_loads:
        counts[k] = counts.get(k, 0) + 1
    now, events, size = 0.0, 0, len(initial_loads)
    while size:
        if events >= max_events or size > max_cells:
            return None, events
        total = math.fsum(number * _source_rate(p, k)
                          for k, number in counts.items())
        if total == 0:
            return size, events
        wait = rng.expovariate(total)
        if now + wait > horizon:
            return size, events
        now += wait
        selector = rng.random() * total
        selected = None
        for k, number in counts.items():
            weight = number * _source_rate(p, k)
            if selector < weight:
                selected = k
                selector /= number
                break
            selector -= weight
        if selected is None:
            # Protect against a floating-point summation endpoint; this event
            # has no effect on the mathematical model and should be exceptional.
            selected = next(reversed(counts))
            selector = rng.random() * _source_rate(p, selected)
        k = selected
        event = _source_event(p, k, selector)
        events += 1
        if event == "birth":
            counts[k] += 1
            size += 1
        else:
            counts[k] -= 1
            if counts[k] == 0:
                del counts[k]
            if event == "gain":
                counts[k + 1] = counts.get(k + 1, 0) + 1
            elif event == "loss":
                counts[k - 1] = counts.get(k - 1, 0) + 1
            else:
                size -= 1
    return 0, events


def sample_source_at_time(p, k, time, trials, seed, *, initial_loads=None,
                          max_events=200000, max_cells=100000):
    """Gillespie source counts at a fixed time, independently of lifetime code."""
    initial_loads = _loads([k] if initial_loads is None else initial_loads)
    if not math.isfinite(time) or time < 0:
        raise ValueError("time must be finite and nonnegative")
    for value, name in ((trials, "trials"), (max_events, "max_events"),
                        (max_cells, "max_cells")):
        _positive_integer(value, name)
    rng = random.Random(seed)
    counts, events = [], 0
    for _ in range(trials):
        count, used = _source_at_time_trial(
            rng, p, initial_loads, time, max_events, max_cells)
        counts.append(count)
        events += used
    censored = sum(value is None for value in counts)
    return {"counts": counts, "censored": censored, "trials": trials,
            "seed": seed, "events": events,
            "mean": math.fsum(counts) / trials if not censored else None,
            "extinct_frequency": sum(value == 0 for value in counts) / trials}


@dataclass(frozen=True)
class LabelledCell:
    """Genetic antigen IDs remain present when presentation is hidden."""

    antigens: frozenset
    hidden: bool = False

    @property
    def presented_antigens(self):
        return frozenset() if self.hidden else self.antigens


def divide_labelled(cell):
    """Immutable antigen sets are copied by value into both descendants."""
    return (LabelledCell(cell.antigens, cell.hidden),
            LabelledCell(cell.antigens, cell.hidden))


def gain_labelled(cell, antigen_id):
    if antigen_id in cell.antigens:
        raise ValueError("a newly acquired antigen must have a new identifier")
    return LabelledCell(cell.antigens | {antigen_id}, cell.hidden)


def lose_labelled(cell, antigen_id):
    if antigen_id not in cell.antigens:
        raise ValueError("cannot lose an absent antigen")
    return LabelledCell(cell.antigens - {antigen_id}, cell.hidden)


def hide_labelled(cell):
    return LabelledCell(cell.antigens, True)


def grow_labelled(p, target_size, seed, *, initial_antigens=(),
                 max_events=100000, max_time=math.inf):
    """Grow one labelled founder; no conditioning by silently restarting.

    A result can be ``target``, ``extinct``, ``event_limit``, ``time_limit`` or
    ``frozen``.  The treatment model conditions on the returned population;
    this sampler makes no assertion about a clinical founder-composition law.
    Every gain uses a unique ID.  Both daughters inherit all parental IDs;
    a later loss only changes the affected daughter's immutable set.
    """
    _positive_integer(target_size, "target_size")
    _positive_integer(max_events, "max_events")
    initial_antigens = frozenset(initial_antigens)
    _loads(initial_antigens)
    if max_time < 0 or math.isnan(max_time):
        raise ValueError("max_time must be nonnegative")
    next_antigen = max(initial_antigens, default=-1) + 1
    cells = [LabelledCell(initial_antigens)]
    rng, now, events = random.Random(seed), 0.0, 0
    status = "target"
    while cells and len(cells) < target_size:
        if events >= max_events:
            status = "event_limit"
            break
        rates = [(p.b_hla + p.d_hla if cell.hidden else
                  _source_rate(p, len(cell.antigens))) for cell in cells]
        total = math.fsum(rates)
        if total == 0:
            status = "frozen"
            break
        wait = rng.expovariate(total)
        if now + wait > max_time:
            status = "time_limit"
            break
        now += wait
        selector = rng.random() * total
        index = len(cells) - 1
        for candidate, rate in enumerate(rates):
            if selector < rate:
                index = candidate
                break
            selector -= rate
        cell = cells[index]
        if cell.hidden:
            event = "birth" if selector < p.b_hla else "death"
        else:
            event = _source_event(p, len(cell.antigens), selector)
        events += 1
        if event == "birth":
            first, second = divide_labelled(cell)
            cells[index] = first
            cells.append(second)
        elif event == "death":
            cells.pop(index)
        elif event == "gain":
            cells[index] = gain_labelled(cell, next_antigen)
            next_antigen += 1
        elif event == "loss":
            antigen = rng.choice(sorted(cell.antigens))
            cells[index] = lose_labelled(cell, antigen)
        else:
            cells[index] = hide_labelled(cell)
    if not cells:
        status = "extinct"
    return {"status": status, "cells": tuple(cells), "time": now,
            "events": events, "next_antigen": next_antigen, "seed": seed}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20261117)
    args = parser.parse_args()
    # Match verifier.py's declared synthetic example exactly, without
    # importing it or sharing its calculations with this simulation.
    p = Parameters(b=1, d=0.5, gain=1, loss=0.25, kappa=1,
                   hla=0.1, b_hla=1, d_hla=0.25)
    summary = {"label": "synthetic nonclinical computational mathematics",
               "parameters": p.__dict__,
               "escape_from_load_1": sample_escape(
                   p, [1], args.trials, args.seed),
               "escape_from_load_8": sample_escape(
                   p, [8], args.trials, args.seed + 1)}
    baseline = Parameters(b=1, d=0.5, gain=1, loss=0.25, kappa=1,
                          hla=0, b_hla=1, d_hla=0.25)
    occupancy = source_lifetime_samples(
        baseline, [1], args.trials, args.seed + 2)
    summary["baseline_source_cell_time"] = {
        key: value for key, value in occupancy.items()
        if key not in ("areas", "partial_areas")}
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
