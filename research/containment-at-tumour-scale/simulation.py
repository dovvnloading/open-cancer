#!/usr/bin/env python3
"""Independent synthetic hybrid-population simulation; no treatment advice.

This program simulates the model specified in the companion README.  It does
not import the verifier, compute an establishment probability, or sample the
limiting mixture asserted by the theorems.  Low-count resistant cells undergo
individual birth/death/immigration events with an exactly integrated waiting
hazard.  Only after a mesoscopic count is reached are they replaced by an ODE.
The sensitive compartment is a deterministic bulk in *both* phases.  This is
not an exact simulation of a fully discrete, two-population tumour.

All parameters and results are synthetic, nonclinical computational
mathematics.  Rates use an arbitrary common time unit.  Only Python's standard
library is used; no network access or result-file writes occur.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
import random
import statistics


@dataclass(frozen=True)
class Parameters:
    """Demographic rates, competition, progression and common dose law."""

    b_s: float = 1.0
    d_s: float = 0.7
    b_r: float = 0.8
    d_r: float = 0.6
    alpha: float = 0.12
    beta: float = 0.02
    a_sr: float = 0.02
    s0: float = 1.0
    epsilon: float = 0.2
    course_scale: float = 16.0
    emax: float = 0.66
    ec50: float = 1.0
    u_max: float = 1.0
    hill: float = 2.0

    @property
    def r_s(self) -> float:
        return self.b_s - self.d_s

    @property
    def threshold(self) -> float:
        return self.s0 + self.epsilon

    @property
    def k_max(self) -> float:
        return hill_kill(self.u_max, self)


@dataclass(frozen=True)
class Trial:
    time: float
    switched: bool
    switched_before_course: bool
    switch_time: float | None
    exposure: float
    max_dose: float
    events: int
    min_sensitive_count: float


def hill_kill(dose: float, p: Parameters) -> float:
    """The same bounded Emax/Hill law is used for both policies."""
    if dose <= 0.0:
        return 0.0
    ratio = (dose / p.ec50) ** p.hill
    return p.emax * ratio / (1.0 + ratio)


def inverse_hill(kill: float, p: Parameters) -> float:
    if kill <= 0.0:
        return 0.0
    if kill >= p.k_max:
        return p.u_max
    return p.ec50 * (kill / (p.emax - kill)) ** (1.0 / p.hill)


def _control(policy: str, r: float, on_course: bool,
             r_sensitive: float, p: Parameters) -> tuple[float, float]:
    if not on_course:
        return 0.0, 0.0
    if policy == "M":
        return p.k_max, p.u_max
    kill = min(p.k_max, max(0.0, r_sensitive - p.a_sr * r))
    return kill, inverse_hill(kill, p)


def _hazard_integral(t: float, constant: float, varying: float,
                     growth: float) -> float:
    if growth == 0.0:
        return (constant + varying) * t
    return constant * t + varying * math.expm1(growth * t) / growth


def _event_wait(exponential: float, constant: float, varying: float,
                growth: float, cap: float) -> float:
    """Invert int_0^t [constant + varying exp(growth*v)] dv.

    ``cap`` is an exact intervening deterministic boundary.  An infinite
    return value means that no stochastic event precedes that boundary.
    Safeguarded Newton iteration controls floating-point root error; there is
    no fixed-step approximation to the stochastic event rates.
    """
    if cap <= 0.0:
        return math.inf
    if math.isfinite(cap):
        available = _hazard_integral(cap, constant, varying, growth)
        if exponential >= available:
            return math.inf
    if constant + varying == 0.0:
        return math.inf
    if growth == 0.0:
        return exponential / (constant + varying)
    if constant == 0.0:
        argument = growth * exponential / varying
        if argument <= -1.0:
            return math.inf
        return math.log1p(argument) / growth

    lo = 0.0
    if growth < 0.0:
        hi = exponential / constant
    else:
        hi = exponential / (constant + varying)
    if math.isfinite(cap):
        hi = min(hi, cap)
    t = min(exponential / (constant + varying), hi)
    for _ in range(40):
        value = _hazard_integral(t, constant, varying, growth) - exponential
        if abs(value) <= 2e-13 * (1.0 + exponential):
            return t
        if value > 0.0:
            hi = t
        else:
            lo = t
        derivative = constant + varying * math.exp(growth * t)
        proposal = t - value / derivative
        t = proposal if lo < proposal < hi else 0.5 * (lo + hi)
    raise ArithmeticError("integrated stochastic hazard did not converge")


def _derivatives(s: float, r: float, policy: str, on_course: bool,
                 r_sensitive: float, immigration: float, p: Parameters
                 ) -> tuple[float, float, float]:
    kill, dose = _control(policy, r, on_course, r_sensitive, p)
    ds = (r_sensitive - p.a_sr * r - kill) * s
    dr = (p.b_r - p.d_r - p.alpha * s - p.beta * r) * r + immigration * s
    return ds, dr, dose


def _rk4(s: float, r: float, step: float, policy: str, on_course: bool,
         r_sensitive: float, immigration: float, p: Parameters
         ) -> tuple[float, float, float]:
    k1 = _derivatives(s, r, policy, on_course, r_sensitive, immigration, p)
    k2 = _derivatives(s + step * k1[0] / 2, r + step * k1[1] / 2,
                      policy, on_course, r_sensitive, immigration, p)
    k3 = _derivatives(s + step * k2[0] / 2, r + step * k2[1] / 2,
                      policy, on_course, r_sensitive, immigration, p)
    k4 = _derivatives(s + step * k3[0], r + step * k3[1],
                      policy, on_course, r_sensitive, immigration, p)
    next_s = s + step * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6
    next_r = r + step * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
    exposure = step * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6
    if next_s <= 0.0 or next_r < 0.0:
        raise ArithmeticError("ODE step lost positivity; reduce ode_step")
    return next_s, next_r, exposure


def trial(policy: str, population: int, supply: float, rng: random.Random,
          p: Parameters = Parameters(), ode_step: float = 0.5,
          switch_exponent: float = 0.25) -> Trial:
    """One complete, uncensored total-burden progression time.

    The course is H=course_scale*log(population); both policies stop at H.
    Unlike an establishment-only calculation, this follows regrowth after H
    when on-course resistant lineages have died out.
    """
    if policy not in ("C", "M"):
        raise ValueError("policy must be C or M")
    if population < 2 or not 0.0 <= supply <= population:
        raise ValueError("require population >= 2 and 0 <= supply <= population")
    if not 0.0 < switch_exponent < 1.0 or ode_step <= 0.0:
        raise ValueError("invalid ODE step or switching exponent")
    if min(p.b_s, p.b_r, p.s0, p.epsilon, p.course_scale,
           p.emax, p.ec50, p.u_max, p.hill) <= 0:
        raise ValueError("birth, scale and pharmacology parameters must be positive")
    if min(p.d_s, p.d_r, p.alpha, p.beta, p.a_sr) < 0:
        raise ValueError("death and competition parameters must be nonnegative")
    mu = supply / population
    r_sensitive = (1.0 - mu) * p.b_s - p.d_s
    if r_sensitive <= p.a_sr * p.threshold:
        raise ValueError("post-course sensitive growth must stay positive before progression")
    horizon = p.course_scale * math.log(population)
    if switch_exponent == 0.25:
        cutoff = math.isqrt(math.isqrt(population))
        cutoff += cutoff ** 4 < population
    else:
        cutoff = math.ceil(population ** switch_exponent)
    t = 0.0
    s = p.s0
    resistant = 0
    exposure = 0.0
    max_dose = 0.0
    events = 0
    min_sensitive = population * s

    # With R fixed, S evolves exponentially and the total jump hazard is
    # A0 + A1 exp(v dt), including competition deaths and mutant births.
    while resistant < cutoff:
        r = resistant / population
        if s + r >= p.threshold:
            return Trial(t, False, False, None, exposure, max_dose,
                         events, min_sensitive)
        on_course = t < horizon
        kill, dose = _control(policy, r, on_course, r_sensitive, p)
        max_dose = max(max_dose, dose)
        growth = r_sensitive - p.a_sr * r - kill
        to_course_end = horizon - t if on_course else math.inf
        to_progression = (math.log((p.threshold - r) / s) / growth
                          if growth > 0.0 else math.inf)
        boundary = min(to_course_end, to_progression)
        constant = (p.b_r + p.d_r + p.beta * r) * resistant
        varying = (p.alpha * resistant + supply * p.b_s) * s
        wait = _event_wait(rng.expovariate(1.0), constant, varying,
                           growth, boundary)
        elapsed = min(wait, boundary)
        if not math.isfinite(elapsed):
            raise ArithmeticError("no event or deterministic progression boundary")
        s *= math.exp(growth * elapsed)
        exposure += dose * elapsed
        t += elapsed
        min_sensitive = min(min_sensitive, population * s)
        if wait >= boundary:
            if to_progression <= to_course_end:
                return Trial(t, False, False, None, exposure, max_dose,
                             events, min_sensitive)
            t = horizon  # Exact boundary: next iteration is dose-off.
            continue
        birth = p.b_r * resistant
        death = (p.d_r + p.alpha * s + p.beta * r) * resistant
        mutation = supply * p.b_s * s
        event = rng.random() * (birth + death + mutation)
        resistant += -1 if birth <= event < birth + death else 1
        events += 1

    switch_time = t
    r = resistant / population
    immigration = supply * p.b_s / population
    # A deterministic bound follows just from S regrowth; it is not a censor.
    lower_growth = r_sensitive - p.a_sr * p.threshold
    upper_time = horizon + (
        math.log(p.threshold / p.s0)
        + (p.k_max - r_sensitive + p.a_sr * p.threshold) * horizon
    ) / lower_growth + 1.0
    while s + r < p.threshold:
        if t > upper_time:
            raise ArithmeticError("trajectory exceeded its sensitive-regrowth bound")
        on_course = t < horizon
        to_course_end = horizon - t if on_course else math.inf
        step = min(ode_step, to_course_end)
        ends_course = on_course and step == to_course_end
        max_dose = max(max_dose, _control(policy, r, on_course, r_sensitive, p)[1])
        next_s, next_r, increment = _rk4(
            s, r, step, policy, on_course, r_sensitive, immigration, p)
        if next_s + next_r >= p.threshold:
            # Find the threshold within the final RK step.  This removes the
            # step-size rounding bias in recorded progression times.
            lo, hi = 0.0, step
            for _ in range(32):
                middle = (lo + hi) / 2.0
                ms, mr, _ = _rk4(s, r, middle, policy, on_course,
                                 r_sensitive, immigration, p)
                if ms + mr >= p.threshold:
                    hi = middle
                else:
                    lo = middle
            step = (lo + hi) / 2.0
            next_s, next_r, increment = _rk4(
                s, r, step, policy, on_course, r_sensitive, immigration, p)
            exposure += increment
            t += step
            min_sensitive = min(min_sensitive, population * next_s)
            break
        t += step
        if ends_course:
            t = horizon
        s, r = next_s, next_r
        exposure += increment
        min_sensitive = min(min_sensitive, population * s)

    return Trial(t, True, switch_time < horizon, switch_time,
                 exposure, max_dose, events, min_sensitive)


def simulate(policy: str, population: int, supply: float, samples: int,
             seed: int, p: Parameters = Parameters(), ode_step: float = 0.5,
             switch_exponent: float = 0.25) -> dict[str, float | int | str]:
    """Independent trials and ordinary sample standard errors, without formulas."""
    if samples < 2:
        raise ValueError("at least two trials are needed to estimate uncertainty")
    rng = random.Random(seed)
    times = []
    switches = 0
    exposures = []
    max_dose = 0.0
    min_sensitive = math.inf
    events = 0
    for _ in range(samples):
        result = trial(policy, population, supply, rng, p, ode_step, switch_exponent)
        times.append(result.time)
        switches += result.switched_before_course
        exposures.append(result.exposure)
        max_dose = max(max_dose, result.max_dose)
        min_sensitive = min(min_sensitive, result.min_sensitive_count)
        events += result.events
    mean = statistics.fmean(times)
    se = statistics.stdev(times) / math.sqrt(samples)
    probability = switches / samples
    return {
        "policy": policy, "population": population, "supply": supply,
        "samples": samples, "seed": seed,
        "mean_time": mean, "se_time": se,
        "mean_time_over_log_n": mean / math.log(population),
        "se_time_over_log_n": se / math.log(population),
        "switch_probability": probability,
        "switch_se": math.sqrt(probability * (1 - probability) / samples),
        "mean_exposure": statistics.fmean(exposures),
        "max_exposure": max(exposures),
        "exposure_budget": p.u_max * p.course_scale * math.log(population),
        "max_dose": max_dose, "min_sensitive_count": min_sensitive,
        "mean_stochastic_events": events / samples,
        "ode_step": ode_step, "switch_exponent": switch_exponent,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261117)
    parser.add_argument("--populations", type=int, nargs="+", default=[10**8, 10**11])
    parser.add_argument("--supplies", type=float, nargs="+", default=[0.02, 1.0])
    parser.add_argument("--ode-step", type=float, default=0.5)
    parser.add_argument("--switch-exponent", type=float, default=0.25)
    args = parser.parse_args()
    print(json.dumps({"status": "synthetic nonclinical computational mathematics; no treatment advice",
                      "parameters": asdict(Parameters())}, sort_keys=True))
    index = 0
    for population in args.populations:
        for supply in args.supplies:
            for policy in ("C", "M"):
                result = simulate(policy, population, supply, args.samples,
                                  args.seed + index, ode_step=args.ode_step,
                                  switch_exponent=args.switch_exponent)
                print(json.dumps(result, sort_keys=True))
                index += 1


if __name__ == "__main__":
    main()
