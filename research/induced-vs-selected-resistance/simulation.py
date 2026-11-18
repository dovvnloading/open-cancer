#!/usr/bin/env python3
"""Independent event simulation for the endpoint-count identifiability model.

NONCLINICAL COMPUTATIONAL MATHEMATICS. All generated data are synthetic.
This file makes no treatment recommendations and supplies no biological validation.

The random engine does not import ``verifier`` or evaluate its count probabilities,
generating function, factorial-cumulant formula, or switching occupancy formula.
It samples Poisson event arrivals, clone marks, and individual cell state changes.
The standard-library RNG and floating-point waiting times make this a numerical
simulation, not the exact-arithmetic part of the accompanying contribution.
"""

import argparse
import json
import math
import random
from collections.abc import Mapping


DEFAULT_CLONE_PMF = {1: 0.25, 2: 0.25, 4: 0.5}


def _nonnegative_finite(value, name):
    """Normalize a numeric simulation input, rejecting infinities and NaNs."""
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite nonnegative number")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite nonnegative number") from exc
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"{name} must be a finite nonnegative number")
    return result


def _probability(value, name, strictly_positive=False):
    value = _nonnegative_finite(value, name)
    if value > 1 or (strictly_positive and value == 0):
        interval = "(0, 1]" if strictly_positive else "[0, 1]"
        raise ValueError(f"{name} must be in {interval}")
    return value


def _poisson(rng, mean):
    """Count unit-rate exponential arrivals in a window of length ``mean``."""
    if mean == 0:
        return 0
    count = 0
    arrival = rng.expovariate(1.0)
    while arrival < mean:
        count += 1
        arrival += rng.expovariate(1.0)
    return count


def _clone_sampler(rng, clone_pmf, geometric_p):
    """Build a mark sampler once, independently of all count-law machinery."""
    if clone_pmf is not None and geometric_p is not None:
        raise ValueError("supply clone_pmf or geometric_p, not both")
    if geometric_p is not None:
        probability = _probability(geometric_p, "geometric_p", True)

        def geometric_mark():
            size = 1
            while rng.random() >= probability:
                size += 1
            return size

        return geometric_mark

    if clone_pmf is None:
        clone_pmf = DEFAULT_CLONE_PMF
    if not isinstance(clone_pmf, Mapping) or not clone_pmf:
        raise ValueError("clone_pmf must be a nonempty positive-size probability mapping")
    if any(isinstance(size, bool) or not isinstance(size, int) or size < 1
           for size in clone_pmf):
        raise ValueError("clone sizes must be positive integers")
    marks = []
    weights = []
    for size, probability in sorted(clone_pmf.items()):
        probability = _probability(probability, "clone probability")
        if probability:
            marks.append(size)
            weights.append(probability)
    total = math.fsum(weights)
    if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("clone probabilities must sum to one")
    # Tiny rounding error in decimal-to-binary conversion is absorbed only here.
    cumulative = []
    running = 0.0
    for weight in weights:
        running += weight / total
        cumulative.append(running)
    cumulative[-1] = 1.0

    def finite_mark():
        draw = rng.random()
        for size, boundary in zip(marks, cumulative):
            if draw < boundary:
                return size
        raise AssertionError("unreachable categorical sampling boundary")

    return finite_mark


def _tolerant_at_end(rng, duration, forward, reverse):
    """Simulate T->S and S->T exponential clocks, starting in T.

    ``reverse`` is the T->S rate; ``forward`` is the S->T rate. Neither
    the stationary fraction nor the closed-form transition probability is used.
    """
    tolerant = True
    elapsed = 0.0
    while elapsed < duration:
        rate = reverse if tolerant else forward
        if rate == 0:
            break
        elapsed += rng.expovariate(rate)
        if elapsed >= duration:
            break
        tolerant = not tolerant
    return tolerant


def sample_counts(
    n,
    seed,
    a,
    b,
    c,
    clone_pmf=None,
    geometric_p=None,
    exposure=1.0,
    retention=1.0,
    holiday_time=None,
    forward_rate=None,
    reverse_rate=None,
):
    """Return independent resistant-colony counts for ``n`` synthetic cultures.

    Each culture has independent Poisson(a), Poisson(c), and Poisson(b*exposure)
    event counts for stable clones, reversible clones, and induced singleton
    founders. Both clone types use the same mark distribution. ``clone_pmf`` is
    a mapping from positive integer clone sizes to probabilities. Alternatively,
    ``geometric_p`` gives a geometric mark on {1, 2, ...} by Bernoulli trials.

    With no ``holiday_time``, each cell in a reversible clone is retained by an
    independent Bernoulli(retention) draw. With ``holiday_time`` specified, both
    switching rates are required and every such cell instead follows direct
    continuous-time T<->S switching clocks. All cells start in T. The retention
    argument must remain its default in that mode, to avoid silently ignoring it.
    Cells neither divide nor die during the modeled holiday. The permanent
    components never switch, and there is no post-induction clone amplification.
    """
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError("n must be a nonnegative integer")
    a = _nonnegative_finite(a, "a")
    b = _nonnegative_finite(b, "b")
    c = _nonnegative_finite(c, "c")
    exposure = _nonnegative_finite(exposure, "exposure")
    retention = _probability(retention, "retention")
    induced_mean = b * exposure
    if not math.isfinite(induced_mean):
        raise ValueError("b * exposure must be finite")
    clock_mode = holiday_time is not None
    if clock_mode:
        if forward_rate is None or reverse_rate is None:
            raise ValueError("holiday_time requires both switching rates")
        if retention != 1:
            raise ValueError("retention must be 1 when switching clocks are requested")
        duration = _nonnegative_finite(holiday_time, "holiday_time")
        forward = _nonnegative_finite(forward_rate, "forward_rate")
        reverse = _nonnegative_finite(reverse_rate, "reverse_rate")
    elif forward_rate is not None or reverse_rate is not None:
        raise ValueError("switching rates require holiday_time")
    rng = random.Random(seed)
    mark = _clone_sampler(rng, clone_pmf, geometric_p)
    counts = []
    for _ in range(n):
        count = _poisson(rng, induced_mean)
        for _ in range(_poisson(rng, a)):
            count += mark()
        for _ in range(_poisson(rng, c)):
            size = mark()
            if clock_mode:
                for _ in range(size):
                    count += _tolerant_at_end(rng, duration, forward, reverse)
            elif retention == 1:
                count += size
            elif retention != 0:
                count += sum(rng.random() < retention for _ in range(size))
        counts.append(count)
    return counts


def summarize(counts):
    """Return descriptive empirical statistics, without confidence claims.

    Factorial cumulants are plug-in statistics and are not unbiased estimators.
    Unlike raw high moments, the displayed empirical PGF statistics are bounded.
    These summaries do not fit parameters or diagnose biological mechanisms.
    """
    counts = list(counts)
    if not counts:
        raise ValueError("at least one count is required")
    if any(isinstance(k, bool) or not isinstance(k, int) or k < 0 for k in counts):
        raise ValueError("counts must be nonnegative integers")
    n = len(counts)
    first = sum(counts) / n
    second = sum(k * (k - 1) for k in counts) / n
    third = sum(k * (k - 1) * (k - 2) for k in counts) / n
    return {
        "n": n,
        "mean": first,
        "factorial_cumulants": [
            first,
            second - first * first,
            third - 3 * second * first + 2 * first * first * first,
        ],
        "empirical_pgf": {
            "0": sum(k == 0 for k in counts) / n,
            "1/2": math.fsum(0.5**k for k in counts) / n,
            "2/3": math.fsum((2 / 3) ** k for k in counts) / n,
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cultures", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20261117)
    args = parser.parse_args(argv)
    if args.cultures < 1:
        parser.error("--cultures must be positive")
    common = {"n": args.cultures, "a": 0.7, "b": 0.4, "c": 0.9}
    summaries = {}
    for arm, retention in enumerate((1.0, 0.5, 0.25)):
        counts = sample_counts(seed=args.seed + arm, retention=retention, **common)
        summaries[f"one_way_arm_{arm}"] = summarize(counts)
    for arm, duration in enumerate((0.0, 1.0, 2.0)):
        counts = sample_counts(
            seed=args.seed + 100 + arm,
            holiday_time=duration,
            forward_rate=0.25,
            reverse_rate=0.75,
            **common,
        )
        summaries[f"two_way_arm_{arm}"] = summarize(counts)
    print(json.dumps({
        "status": "synthetic_nonclinical_unreviewed",
        "seed": args.seed,
        "clone_pmf": DEFAULT_CLONE_PMF,
        "parameters": {"a": common["a"], "b": common["b"], "c": common["c"]},
        "summaries": summaries,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
