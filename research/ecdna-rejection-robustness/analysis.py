"""Exact rejection-indicator sensitivity, not segregation inference.

Nonclinical, unreviewed research; no treatment advice. Standard library only.
Reads the existing Lundberg aggregate CSV without changing it. All inferential
comparisons use Fraction; decimal displays are not used to decide rejection.
"""

import argparse
import csv
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys


F = Fraction
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "research" / "ecdna-segregation-clustering-data" / "data" / "lundberg.csv"
PROTOCOL_COMMIT = "d58107d841841dfb42577e10f1f6065c11c909bc"
LOCAL_FREEZE = "acc2dc41f8a6673c17be28a4f80170505cf1fb84"
BASE_COMMIT = "9ded69c3586188b173e81fedb7a79abbf4e7e91f"
LEVEL = F(1, 20)
SLICES = tuple(map(F, ("0", ".0001", ".0005", ".001", ".002", ".005",
                      ".01", ".02", ".05", ".1", ".25", ".5", ".75", ".9",
                      ".99", "1")))
MODELS = ("independent", "beta", "exchangeable", "markov")
BRACKET_DENOMINATOR = 10**12
BEGIN = "<!-- BEGIN RECOMPUTED RESULTS -->"
END = "<!-- END RECOMPUTED RESULTS -->"


def rational(value):
    """Require exact inputs; reject floats and bools instead of silently rounding."""
    if isinstance(value, bool) or not isinstance(value, (int, F, str)):
        raise TypeError("Use an integer, Fraction, or exact decimal/rational string")
    return F(value)


def probability(value):
    value = rational(value)
    if not 0 <= value <= 1:
        raise ValueError("Probability/correlation cap must lie in [0,1]")
    return value


def counts(n, k):
    if type(n) is not int or type(k) is not int or n < 1:
        raise ValueError("n must be a positive integer and k an integer")


def normalized_text(path):
    return Path(path).read_bytes().decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")


def normalized_hash(path):
    return hashlib.sha256(normalized_text(path).encode("utf-8")).hexdigest()


def hashes():
    """No self-referential README hash; all keys are platform-independent paths."""
    paths = (SOURCE, HERE / "PROTOCOL.md", HERE / "analysis.py", HERE / "verifier.py",
             ROOT / "tests" / "test_research_ecdna_rejection_robustness.py")
    return {p.relative_to(ROOT).as_posix(): normalized_hash(p) for p in paths}


@lru_cache(maxsize=1)
def load_data():
    with SOURCE.open(encoding="utf-8", newline="") as handle:
        original = list(csv.DictReader(handle))
    by_id = {row["id"]: row for row in original}
    if len(by_id) != len(original):
        raise ValueError("Duplicate source IDs")
    alpha = by_id["lundberg_alpha"]
    if alpha["quantity"] != "per_division_alpha" or F(alpha["value"]) != F(1, 100):
        raise ValueError("Unexpected source nominal level")
    out = []
    for line, slug in (("CHP-212", "chp212"), ("LA-N-5", "lan5"), ("MC-IX", "mcix")):
        ids = ("lundberg_" + slug + "_n", "lundberg_" + slug + "_pct")
        nrow, prow = (by_id[key] for key in ids)
        if (nrow["quantity"], prow["quantity"]) != (
                "scored_divisions_n", "rejected_divisions_percent"):
            raise ValueError("Wrong source quantities")
        for row in (nrow, prow, alpha):
            if row["condition"] != "untreated" or row["doi"] != "10.1371/journal.pone.0003099":
                raise ValueError("Wrong condition/source")
            for field in ("source_url", "locator", "access_date", "reuse_basis", "notes"):
                if not row[field].strip():
                    raise ValueError("Missing source provenance")
            if "Statistical distribution of DMs at mitosis" not in row["locator"]:
                raise ValueError("Unexpected source locator")
        if nrow["population"] != line or prow["population"] != line:
            raise ValueError("Wrong population")
        n = int(nrow["value"])
        if F(nrow["value"]) != n:
            raise ValueError("Nonintegral denominator")
        printed = F(prow["value"])
        compatible = tuple(k for k in range(n+1)
                           if printed-F(1, 2) <= F(100*k, n) <= printed+F(1, 2))
        if not compatible:
            raise ValueError("No compatible integer count")
        out.append({"line": line, "n": n, "k": min(compatible), "counts": compatible,
                    "printed_percent": printed, "source_ids": ids,
                    "nominal": F(alpha["value"])})
    return tuple(out)


def binomial_tail(n, k, A):
    """Exact integer polynomial, evaluating the shorter of tail and complement."""
    counts(n, k)
    A = probability(A)
    if k <= 0:
        return F(1)
    if k > n:
        return F(0)
    a, d = A.numerator, A.denominator
    b = d-a
    if k <= n//2:
        return 1-F(sum(math.comb(n, j)*a**j*b**(n-j) for j in range(k)), d**n)
    return F(sum(math.comb(n, j)*a**j*b**(n-j) for j in range(k, n+1)), d**n)


def beta_tail(n, k, A, R):
    """Exact beta-binomial indicator tail; R is an exact model parameter here."""
    counts(n, k)
    A, R = probability(A), probability(R)
    if k <= 0:
        return F(1)
    if k > n or A == 0:
        return F(0)
    if A == 1:
        return F(1)
    if R == 0:
        return binomial_tail(n, k, A)
    if R == 1:
        return A
    u, v, h, d = A.numerator, A.denominator, R.numerator, R.denominator
    first, second, step = u*(d-h), (v-u)*(d-h), v*h
    total = first+second
    left, right, denom = [1], [1], 1
    for j in range(n):
        left.append(left[-1]*(first+j*step))
        right.append(right[-1]*(second+j*step))
        denom *= total+j*step
    if k <= n//2:
        return 1-F(sum(math.comb(n, j)*left[j]*right[n-j] for j in range(k)), denom)
    return F(sum(math.comb(n, j)*left[j]*right[n-j] for j in range(k, n+1)), denom)


def markov_tail(n, k, A):
    counts(n, k)
    A = probability(A)
    if k <= 0:
        return F(1)
    if k > n:
        return F(0)
    return min(F(1), n*A/k)


def sharp_certificate(n, k, A, R):
    """Sharp finite-exchangeable upper-cap tail plus a primal/dual certificate.

    The dual polynomial is c0+c1*s+c2*s*(s-1). The law's mean is exactly n*A;
    its factorial second moment is at most the cap. It therefore witnesses
    sharpness even though a mere tail inequality would not establish it.
    """
    counts(n, k)
    A, R = probability(A), probability(R)
    if not 1 <= k <= n:
        raise ValueError("Certificate requires 1 <= k <= n")
    mu = n*A
    B = n*(n-1)*(A*A+R*A*(1-A))
    if mu >= k:
        low = mu.numerator//mu.denominator
        law = {low: 1-(mu-low)}
        if mu != low:
            law[low+1] = mu-low
        bound, dual, branch = F(1), (F(1), F(0), F(0)), "constant"
    else:
        choices = [(mu/k, -1)]
        for j in range(k-1):
            D = (k-j)*(k-j-1)
            choices.append(((B-2*j*mu+j*(j+1))/D, j))
        bound, j = min(choices)
        if j == -1:
            law = {0: 1-bound, k: bound}
            dual, branch = (F(0), F(1, k), F(0)), "markov"
        else:
            D = (k-j)*(k-j-1)
            law = {j: j+1-mu+(k-j-1)*bound,
                   j+1: mu-j-(k-j)*bound, k: bound}
            dual = (F(j*(j+1), D), F(-2*j, D), F(1, D))
            branch = "quadratic_"+str(j)
    law = {s: mass for s, mass in law.items() if mass}
    return {"bound": bound, "mean": mu, "second_cap": B, "law": law,
            "dual": dual, "branch": branch}


def sharp_tail(n, k, A, R):
    counts(n, k)
    A, R = probability(A), probability(R)
    if k <= 0:
        return F(1)
    if k > n:
        return F(0)
    return sharp_certificate(n, k, A, R)["bound"]


def holm(values):
    """Exact step-down adjustment, stable line-name tie order."""
    if not values:
        raise ValueError("Empty Holm family")
    values = {key: probability(value) for key, value in values.items()}
    adjusted, previous = {}, F(0)
    for rank, key in enumerate(sorted(values, key=lambda key: (values[key], key))):
        previous = max(previous, min(F(1), (len(values)-rank)*values[key]))
        adjusted[key] = previous
    return {key: adjusted[key] for key in values}


def scenario(model, A, R=F(0), rows=None):
    if model not in MODELS:
        raise ValueError("Unknown count model")
    A, R = probability(A), probability(R)
    rows = load_data() if rows is None else rows
    tails = {}
    for row in rows:
        n, k = row["n"], row["k"]
        if model == "independent":
            value = binomial_tail(n, k, A)
        elif model == "beta":
            value = beta_tail(n, k, A, R)
        elif model == "exchangeable":
            value = sharp_tail(n, k, A, R)
        else:
            value = markov_tail(n, k, A)
        tails[row["line"]] = value
    adjusted = holm(tails)
    rejected = tuple(key for key in tails if adjusted[key] <= LEVEL)
    return {"raw": tails, "adjusted": adjusted, "rejected": rejected,
            "any": bool(rejected), "all": len(rejected) == len(tails),
            "global": min(F(1), len(tails)*min(tails.values()))}


def decision_value(result, target):
    if target == "any":
        return result["global"]
    if target == "all":
        return max(result["adjusted"].values())
    return result["adjusted"][target]


@lru_cache(maxsize=None)
def _boundary(model, R, target):
    low, high = 0, BRACKET_DENOMINATOR
    while high-low > 1:
        middle = (low+high)//2
        value = decision_value(scenario(model, F(middle, BRACKET_DENOMINATOR), R), target)
        if value <= LEVEL:
            low = middle
        else:
            high = middle
    return F(low, BRACKET_DENOMINATOR), F(high, BRACKET_DENOMINATOR)


def boundary(model, R=F(0), target="any", rows=None):
    """Largest rejected size; exact grid bracket, never a float root decision."""
    R = probability(R)
    if rows is not None:
        raise ValueError("Boundary report uses the frozen three-line family")
    if target not in ("any", "all", *(row["line"] for row in load_data())):
        raise ValueError("Unknown frontier target")
    return _boundary(model, R, target)


@lru_cache(maxsize=1)
def verifier_module():
    spec = importlib.util.spec_from_file_location("ecdna_robustness_verifier", HERE / "verifier.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def robust_rho_candidates(A):
    """All possible sharp-envelope/Holm breakpoints at fixed size, exactly."""
    A = probability(A)
    candidates = {F(0), F(1)}
    if A in (0, 1):
        return tuple(sorted(candidates))
    for row in load_data():
        n, k = row["n"], row["k"]
        mu, base, slope = n*A, n*(n-1)*A*A, n*(n-1)*A*(1-A)
        if slope == 0:
            continue
        for m in range(1, len(load_data())+1):
            t = LEVEL/m
            for j in range(k-1):
                root = (t*(k-j)*(k-j-1)+2*j*mu-j*(j+1)-base)/slope
                if 0 <= root <= 1:
                    candidates.add(root)
    return tuple(sorted(candidates))


def robust_nominal_crossings():
    A = load_data()[0]["nominal"]
    candidates = robust_rho_candidates(A)
    transitions = []
    for index, root in enumerate(candidates[1:-1], start=1):
        left = scenario("exchangeable", A, (candidates[index-1]+root)/2)["rejected"]
        right = scenario("exchangeable", A, (candidates[index+1]+root)/2)["rejected"]
        if left != right:
            transitions.append({"rho": root, "left": left,
                                "at": scenario("exchangeable", A, root)["rejected"],
                                "right": right})
    return tuple(transitions)


def nominal_or_worse_boundary(A):
    """Closed-form sharp joint R frontier for A >= the nominal size.

    None means no rejection even at R=0. Theorem 5's specialization proves
    that only quadratic indices 0 and 1 can be active on this boundary.
    """
    A = probability(A)
    if A < F(1, 100):
        raise ValueError("This specialization requires A >= 0.01")
    if A >= F(1, 50):
        return None
    numerator = (F(3, 2) if A <= F(23, 2000) else 200*A-F(4, 5))-9900*A*A
    if numerator < 0:
        return None
    return numerator/(9900*A*(1-A))


def critical_count(n, A, R, model="exchangeable", family=3, level=LEVEL):
    """Conservative first-Holm-step design, not observed power or a new sample."""
    if family < 1 or type(family) is not int:
        raise ValueError("Positive integer family size required")
    for k in range(1, n+1):
        if model == "independent":
            p = binomial_tail(n, k, A)
        elif model == "exchangeable":
            p = sharp_tail(n, k, A, R)
        elif model == "markov":
            p = markov_tail(n, k, A)
        else:
            raise ValueError("Unsupported experiment bound")
        if p <= level/family:
            return {"n": n, "k": k, "tail": p, "threshold": level/family,
                    "previous_tail": None if k == 1 else (
                        binomial_tail(n, k-1, A) if model == "independent" else
                        sharp_tail(n, k-1, A, R) if model == "exchangeable" else
                        markov_tail(n, k-1, A)), "model": model, "A": A, "R": R}
    return {"n": n, "k": None, "threshold": level/family, "model": model, "A": A, "R": R}


def decimal(value, places=12, rounding=None):
    value = F(value)
    with localcontext() as ctx:
        ctx.prec = max(60, places+25)
        x = Decimal(value.numerator)/Decimal(value.denominator)
        if rounding is not None:
            x = x.quantize(Decimal(1).scaleb(-places), rounding=rounding)
        return format(x, "."+str(places)+"f")


def bracket_text(bracket):
    lo, hi = bracket
    return "["+decimal(lo, 12, ROUND_FLOOR)+", "+decimal(hi, 12, ROUND_CEILING)+"]"


def exact(value):
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def jsonable(value):
    if isinstance(value, F):
        return exact(value)
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    return value


@lru_cache(maxsize=1)
def results():
    rows = load_data()
    targets = tuple(row["line"] for row in rows)+("any",)
    A = rows[0]["nominal"]
    frontiers = {model: tuple({"rho": R,
                              "boundary": {target: boundary(model, R, target) for target in targets}}
                             for R in SLICES)
                 for model in ("exchangeable", "beta")}
    return {"data": rows, "hashes": hashes(), "level": LEVEL, "nominal": A,
            "baseline": {model: scenario(model, A, F(0)) for model in MODELS},
            "independent_boundary": {target: boundary("independent", F(0), target) for target in targets},
            "markov_boundary": {target: boundary("markov", F(0), target) for target in targets},
            "frontiers": frontiers,
            "robust_crossings": robust_nominal_crossings(),
            "beta_crossings": verifier_module().nominal_crossings(
                {row["line"]: (row["n"], row["k"]) for row in rows}),
            "witnesses": {"uncorrelated": sharp_certificate(100, 10, A, F(0)),
                          "nominal_boundary": sharp_certificate(100, 10, A,
                              robust_nominal_crossings()[0]["rho"])},
            "joint_examples": tuple({"A": a, "R_boundary": nominal_or_worse_boundary(a)}
                                    for a in (A, F(23, 2000), F(3, 250), F(7, 500))),
            "experiments": tuple(critical_count(100, A, R, model) for model, R in (
                ("independent", F(0)), ("exchangeable", F(0)),
                ("exchangeable", F(1, 200)), ("exchangeable", F(1, 100)),
                ("markov", F(1)))),
            "calibration_budget": {"epsilon": F(1, 100), "test_level": F(1, 25),
                "design": critical_count(100, A, F(1, 200), level=F(1, 25))}}


def render_report(result=None):
    """Every displayed calculated number is rebuilt from the frozen inputs."""
    r = results() if result is None else result
    lines = ["### Recomputed source counts and endpoints", "",
             "All displayed probabilities are rounded; all decisions use exact rational arithmetic. "
             "Brackets enclose arithmetic boundary locations, not sampling uncertainty.", "",
             "| Cell line | Mothers n | Compatible rejected k | Independent tail | Independent Holm | Sharp R=0 tail | Sharp R=0 Holm | Markov tail | Markov Holm |",
             "| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in r["data"]:
        label = row["line"]
        values = [r["baseline"][model][kind][label]
                  for model in ("independent", "exchangeable", "markov")
                  for kind in ("raw", "adjusted")]
        lines.append(f"| {label} | {row['n']} | {', '.join(map(str, row['counts']))} | "
                     + " | ".join(decimal(v) for v in values)+" |")
    crossing = r["robust_crossings"][0]["rho"]
    lines += ["", "At nominal actual size A = "+exact(r["nominal"])+", "
              "the sharp exchangeable-cap family rejects **LA-N-5 only** for "
              "**0 <= R <= "+exact(crossing)+" (approximately "+decimal(crossing, 15)+")**. "
              "For larger caps it rejects no line. Equality is included.", "",
              "| Model / dependence assumption | CHP-212 size frontier | LA-N-5 size frontier | MC-IX size frontier | Joint any-rejection size frontier |",
              "| --- | --- | --- | --- | --- |"]
    for model, key in (("Independent indicators", "independent_boundary"),
                       ("Arbitrary dependence (Markov)", "markov_boundary")):
        lines.append("| "+model+" | "+" | ".join(bracket_text(r[key][target])
                     for target in ("CHP-212", "LA-N-5", "MC-IX", "any"))+" |")
    lines += ["", "The dependence-free joint size frontier is exactly **1/600**; "
              "its displayed bracket is intentionally conservative. These frontiers use "
              "the complete three-line Holm family, including both weak panels.", "",
              "### Complete frontier: certified slices of the exact formulas", "",
              "Entries are [lower, upper] rational decimal brackets for the largest "
              "size cap A still rejected. Their width is 0.000000000001. The equations "
              "in Theorems 1–5 define every point between these slices; no linear "
              "interpolation is asserted. The all-rejection frontier is the minimum "
              "of the three line frontiers; the joint any-rejection frontier is their maximum."]
    for model, title in (("exchangeable", "Sharp exchangeable cap: rho_I <= R"),
                         ("beta", "Pointwise beta-binomial: rho_I fixed at r")):
        lines += ["", "#### "+title, "",
                  "| R or r | CHP-212 A frontier | LA-N-5 A frontier | MC-IX A frontier | Joint any A frontier |",
                  "| ---: | --- | --- | --- | --- |"]
        for item in r["frontiers"][model]:
            lines.append("| "+exact(item["rho"])+" | "+" | ".join(
                bracket_text(item["boundary"][target])
                for target in ("CHP-212", "LA-N-5", "MC-IX", "any"))+" |")
    lines += ["", "#### Joint trade-offs at nominal or worse size", "",
              "These points use the closed-form primary frontier proved below. "
              "They are simultaneous size/dependence trade-offs, not separate tolerances "
              "that may be combined arbitrarily.", "",
              "| Size cap A | Largest correlation cap R retaining rejection | Decimal R |",
              "| ---: | ---: | ---: |"]
    for item in r["joint_examples"]:
        lines.append("| "+exact(item["A"])+" | "+exact(item["R_boundary"])+" | "
                     +decimal(item["R_boundary"], 15)+" |")
    lines += ["", "### Complete nominal-size beta-binomial crossing audit", "",
              "The verifier isolates every root of every raw tail minus 0.05/j, "
              "j in {1,2,3}. Exact Bernstein sign certificates exclude missing roots; "
              "the original rational beta tails certify every bracket endpoint. "
              "Only changes in the full Holm rejection set are reported as family transitions.", "",
              "| Line | Holm divisor j | Polynomial degree | Roots on [0,1] |",
              "| --- | ---: | ---: | ---: |"]
    beta = r["beta_crossings"]
    for cert in beta["polynomial_certificates"]:
        lines.append(f"| {cert['line']} | {cert['divisor']} | {cert['degree']} | {cert['root_count']} |")
    lines += ["", f"There are **{len(beta['raw_roots'])} raw crossings** and "
              f"**{len(beta['transitions'])} actual Holm transitions**.", "",
              "| Correlation root bracket | Rejected immediately below | Rejected at the algebraic boundary | Rejected immediately above |",
              "| --- | --- | --- | --- |"]
    def names(labels):
        return ", ".join(labels) if labels else "none"
    for item in beta["transitions"]:
        lines.append("| "+bracket_text((item["lo"], item["hi"]))+" | "
                     +names(item["before"])+" | "+names(item["at_boundary"])+" | "
                     +names(item["after"])+" |")
    lines += ["", "At r=0 the rejection set is **"+names(beta["at_zero"])+"**; "
              "at r=1 it is **"+names(beta["at_one"])+"** and every adjusted value "
              "is exactly **3/100**. Thus the apparent reappearance of rejection is "
              "real within the restricted beta model. It does not occur for an upper "
              "dependence cap and does not make highly dependent data robust evidence.", "",
              "### Attaining adversarial laws for LA-N-5", "",
              "Each row is a hypothetical null count distribution, not reconstructed "
              "observations. Choose a uniformly random subset of K mothers to receive "
              "indicator 1. The laws have mean K=1 and actual per-test size 1/100.", "",
              "| Cap / witness | Pr(K=0) | Pr(K=1) | Pr(K=10) | Tail | Holm value |",
              "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for key, label in (("uncorrelated", "R=0"),
                       ("nominal_boundary", "R="+exact(crossing))):
        witness = r["witnesses"][key]
        lines.append("| "+label+" | "+" | ".join(exact(witness["law"].get(k, 0))
                     for k in (0, 1, 10))+" | "+exact(witness["bound"])+" | "
                     +exact(3*witness["bound"])+" |")
    lines += ["", "The R=0 witness has pairwise independent Bernoulli indicators "
              "but is not a mutually independent vector. It disproves replacing "
              "an uncorrelated-indicator assumption with a binomial count law.", "",
              "### Prospective design arithmetic (not an experiment performed)", "",
              "The following critical counts use n=100 mothers per line, A=1/100, "
              "and the conservative first Holm step 1/60. The previous count fails "
              "that step. Caps here are assumed known; the next example budgets "
              "calibration uncertainty separately.", "",
              "| Indicator assumption | R cap | Minimum rejected mothers | Tail at critical count | Tail one count below |",
              "| --- | ---: | ---: | ---: | ---: |"]
    for item in r["experiments"]:
        lines.append("| "+item["model"]+" | "+exact(item["R"])+" | "+str(item["k"])
                     +" | "+decimal(item["tail"])+" | "+decimal(item["previous_tail"])+" |")
    calibration = r["calibration_budget"]
    design = calibration["design"]
    lines += ["", "For an illustrative simultaneous calibration failure budget "
              +exact(calibration["epsilon"])+" and test-family level "
              +exact(calibration["test_level"])+", their sum is "
              +exact(calibration["epsilon"]+calibration["test_level"])+". With "
              "A="+exact(design["A"])+" and R="+exact(design["R"])+", the first-step "
              "threshold is "+exact(design["threshold"])+" and the minimum count "
              "is **"+str(design["k"])+"/"+str(design["n"])+"** (tail "
              +decimal(design["tail"])+"; one count below "
              +decimal(design["previous_tail"])+"). This error budget is an "
              "illustrative design choice made after the freeze, not achieved "
              "calibration or a power calculation.", "",
              "### Newline-normalized SHA-256 manifest", "",
              "The source, immutable protocol and original implementation/test "
              "files are hashed after UTF-8 CRLF/CR-to-LF normalization. Path keys "
              "come from Path.as_posix(). README is excluded to avoid self-reference.", "",
              "| Path | SHA-256 |", "| --- | --- |"]
    for path, digest in r["hashes"].items():
        lines.append("| "+path+" | `"+digest+"` |")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--markdown", action="store_true")
    args = parser.parse_args()
    if args.check:
        text = normalized_text(HERE / "README.md")
        actual = text.split(BEGIN+"\n", 1)[1].split(END, 1)[0]
        if actual != render_report()+"\n":
            raise SystemExit("README results differ from exact recomputation")
        print("Exact report and newline-normalized hashes match; nonclinical, unreviewed.")
    elif args.json:
        print(json.dumps(jsonable(results()), indent=2, sort_keys=True))
    else:
        print(render_report())


if __name__ == "__main__":
    main()
