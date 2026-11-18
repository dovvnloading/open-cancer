"""Independent fixed-seed ecDNA branching experiments, standard library only.

Nonclinical computational mathematics. No treatment advice or patient data.
This module samples events and complementary sisters directly; it neither
imports the verifier nor uses its matrices, moment formulas, or certificates.
Finite runs cannot prove infinite-time persistence. Resource caps raise errors.
"""

from __future__ import annotations

import argparse
import math
import random


DEFAULT_SEED = 2026111701
METRICS = ("cells", "positive", "zero", "copies", "copy_squares", "extinct_positive")


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(name + " must be an integer")
    if value < minimum:
        raise ValueError(name + " is too small")
    return value


def _probability(value, name):
    result = float(value)
    if not math.isfinite(result) or not 0.0 <= result <= 1.0:
        raise ValueError(name + " must be finite and between zero and one")
    return result


def _kernel(mode, rho, bias):
    if mode not in ("binomial", "beta", "biased"):
        raise ValueError("unknown segregation mode")
    return _probability(rho, "rho"), _probability(bias, "bias")


def sample_split(k, rng, mode="binomial", rho=0.0, bias=0.5):
    """Draw sisters with an exactly conserved total of 2*k copies.

    In beta mode, rho=1/(2*a+1) and each mother draws its own P~Beta(a,a).
    rho=0 is the binomial limit and rho=1 the all-or-nothing endpoint.
    Biased mode labels a daughter with allocation probability ``bias``;
    population summaries contain both sisters and hence symmetrize the bias.
    """
    _integer(k, "copy number")
    rho, bias = _kernel(mode, rho, bias)
    n = 2*k
    if not n:
        return 0, 0
    if mode == "binomial" or (mode == "beta" and rho == 0.0):
        j = rng.getrandbits(n).bit_count()
    elif mode == "beta" and rho == 1.0:
        j = n if rng.random() < 0.5 else 0
    else:
        if mode == "beta":
            a = (1.0-rho)/(2.0*rho)
            probability = rng.betavariate(a, a)
        else:
            probability = bias
        if probability == 0.0:
            j = 0
        elif probability == 1.0:
            j = n
        else:
            j = sum(rng.random() < probability for _ in range(n))
    return j, n-j


def simulate_one(initial, generations, division_probability, rng, *,
                 mode="binomial", rho=0.0, bias=0.5, keep_zeros=True,
                 max_cells=1_000_000, max_copies=1_000_000):
    """Sample one finite genealogy, including extinct populations.

    Birth/death competition is sampled through b/(b+d_k). The continuous-time
    clocks do not change this genealogy, but this is not a calendar-time census.
    A resource cap aborts the run; no capped population is classified extinct.
    Set keep_zeros=False only for positive-lineage observables: resulting pooled
    moments then have a different, positive-only denominator.
    """
    _integer(generations, "generations")
    _integer(max_cells, "max_cells", 1)
    _integer(max_copies, "max_copies", 1)
    _kernel(mode, rho, bias)
    initial = tuple(initial)
    for k in initial:
        _integer(k, "copy number")
    positive = [k for k in initial if k > 0]
    zeros = initial.count(0) if keep_zeros else 0
    p0 = _probability(division_probability(0), "zero division probability") if keep_zeros else 0.0

    def check_size(cells, copies):
        if cells > max_cells:
            raise RuntimeError("cell resource cap reached; simulation is censored, not extinct")
        if copies > max_copies:
            raise RuntimeError("copy resource cap reached; simulation is censored, not extinct")

    check_size(len(positive)+zeros, max(positive, default=0))
    for _ in range(generations):
        if keep_zeros and p0 == 1.0:
            next_zeros = 2*zeros
        elif keep_zeros:
            next_zeros = 2*sum(rng.random() < p0 for _ in range(zeros))
        else:
            next_zeros = 0
        offspring = []
        check_size(next_zeros, 0)
        for k in positive:
            probability = _probability(division_probability(k), "division probability")
            if rng.random() < probability:
                check_size(len(offspring)+next_zeros+2, 2*k)
                for j in sample_split(k, rng, mode=mode, rho=rho, bias=bias):
                    if j:
                        offspring.append(j)
                    elif keep_zeros:
                        next_zeros += 1
        positive, zeros = offspring, next_zeros
        check_size(len(positive)+zeros, max(positive, default=0))
    return {
        "positive_copies": tuple(positive),
        "zero_count": zeros,
        "cells": len(positive)+zeros,
        "positive": len(positive),
        "zero": zeros,
        "copies": sum(positive),
        "copy_squares": sum(k*k for k in positive),
        "extinct_positive": int(not positive),
    }


def run_ensemble(initial, generations, division_probability, *, trials=10000,
                 seed=DEFAULT_SEED, mode="binomial", rho=0.0, bias=0.5,
                 keep_zeros=True, max_cells=1_000_000, max_copies=1_000_000):
    """Estimate means and pooled moments from independent replicate populations.

    Monte Carlo uncertainty uses populations, never related cells, as replicates.
    Ratio standard errors are first-order delta estimates, not exact coverage
    guarantees. An undefined pooled quantity and its SE are returned as None.
    """
    _integer(trials, "trials", 2)
    initial = tuple(initial)
    rng = random.Random(seed)
    dimension = len(METRICS)
    means = [0.0]*dimension
    cross = [[0.0]*dimension for _ in range(dimension)]
    for trial in range(1, trials+1):
        state = simulate_one(
            initial, generations, division_probability, rng, mode=mode,
            rho=rho, bias=bias, keep_zeros=keep_zeros,
            max_cells=max_cells, max_copies=max_copies,
        )
        values = [float(state[name]) for name in METRICS]
        delta = [x-m for x, m in zip(values, means)]
        for i in range(dimension):
            means[i] += delta[i]/trial
        for i in range(dimension):
            for j in range(dimension):
                cross[i][j] += delta[i]*(values[j]-means[j])
    covariance = [[value/(trials-1) for value in row] for row in cross]
    averages = dict(zip(METRICS, means))
    variances = {name: max(0.0, covariance[i][i]) for i, name in enumerate(METRICS)}
    errors = {name: math.sqrt(value/trials) for name, value in variances.items()}

    def delta_se(gradient):
        variance = sum(gradient[i]*covariance[i][j]*gradient[j]
                       for i in range(dimension) for j in range(dimension))
        return math.sqrt(max(0.0, variance)/trials)

    pooled = {}
    for name in ("mean", "variance", "fano", "zero_fraction"):
        pooled[name] = None
        pooled[name+"_standard_error"] = None
    n, _, z, s, q, _ = means
    if n > 0:
        mu = s/n
        pooled["mean"] = mu
        pooled["variance"] = max(0.0, q/n-mu*mu)
        pooled["zero_fraction"] = z/n
        pooled["mean_standard_error"] = delta_se([-s/n**2, 0, 0, 1/n, 0, 0])
        pooled["variance_standard_error"] = delta_se(
            [-q/n**2+2*s*s/n**3, 0, 0, -2*s/n**2, 1/n, 0])
        pooled["zero_fraction_standard_error"] = delta_se([-z/n**2, 0, 1/n, 0, 0, 0])
        if s > 0:
            pooled["fano"] = q/s-s/n
            pooled["fano_standard_error"] = delta_se(
                [s/n**2, 0, 0, -q/s**2-1/n, 1/s, 0])
    return {
        "seed": seed, "trials": trials, "generations": generations,
        "mode": mode, "rho": rho, "bias": bias, "keep_zeros": keep_zeros,
        "means": averages, "sample_variances": variances,
        "standard_errors": errors, "pooled": pooled,
    }


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=10000)
    args = parser.parse_args()
    print("Nonclinical computational mathematics. No treatment advice.")
    print("Finite generations only; none of these runs proves persistence.")
    cases = (
        ("linear lower c=3/8", "binomial", 0.0, 0.5, lambda k: 1/(1+3*k/8)),
        ("linear upper c=2/5", "binomial", 0.0, 0.5, lambda k: 1/(1+2*k/5)),
        ("linear beta rho=1/3 c=1/4", "beta", 1/3, 0.5, lambda k: 1/(1+k/4)),
        ("linear biased p=3/4 c=3/10", "biased", 0.0, 0.75, lambda k: 1/(1+3*k/10)),
        ("quadratic c=3/10", "binomial", 0.0, 0.5, lambda k: 1/(1+3*k*k/10)),
        ("saturating c=2 h=3", "binomial", 0.0, 0.5, lambda k: 1/(1+2*k/(3+k))),
    )
    for index, (label, mode, rho, bias, probability) in enumerate(cases):
        report = run_ensemble((1,), 5, probability, trials=args.trials,
                              seed=DEFAULT_SEED+index, mode=mode, rho=rho, bias=bias)
        print(label, "seed", report["seed"])
        for metric in ("positive", "extinct_positive"):
            print("  {} {:.9f} SE {:.9f}".format(
                metric, report["means"][metric], report["standard_errors"][metric]))
        print("  pooled Fano {:.9f} SE {:.9f}".format(
            report["pooled"]["fano"], report["pooled"]["fano_standard_error"]))
    print("Same-variance pair: one-copy founders, generation two, c=3/10.")
    for index, (mode, rho, bias) in enumerate((("beta", 1/4, 0.5),
                                             ("biased", 0.0, 3/4))):
        report = run_ensemble((1,), 2, lambda k: 1/(1+3*k/10),
                              trials=args.trials, seed=DEFAULT_SEED+20+index,
                              mode=mode, rho=rho, bias=bias)
        print("  {} seed {} Fano {:.9f} SE {:.9f}; zero {:.9f} SE {:.9f}".format(
            mode, report["seed"], report["pooled"]["fano"],
            report["pooled"]["fano_standard_error"],
            report["pooled"]["zero_fraction"],
            report["pooled"]["zero_fraction_standard_error"]))


if __name__ == "__main__":
    _main()
