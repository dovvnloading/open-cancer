#!/usr/bin/env python3
"""Exact, synthetic checks for nonclinical clonal-tree mathematics.

No patient data, floating-point arithmetic, or third-party packages are used.
A tree is a parent tuple on labels 0,...,n-1 with exactly one entry -1.
That root is a malignant founder; the unobserved normal population is not
another tree vertex.  Matrix rows index samples, and B[i][j] is one exactly
when j is an ancestor of i, including i itself.

Run ``python research/clonal-tree-identifiability/verifier.py`` from the
repository root to reproduce the small-n table and counterexample checks.
The table enumerates trees, not a finite grid purported to cover all real
frequency matrices; its maximum claims additionally need the README proof.
"""

from argparse import ArgumentParser
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb, factorial, gcd


def _positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(name + " must be an integer")
    if value < 1:
        raise ValueError(name + " must be positive")
    return value


def _exact(value):
    """Accept exact rational inputs; reject accidental float approximations."""
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("entries must be integers or fractions.Fraction values")
    return Fraction(value)


def _parents(parents):
    result = tuple(parents)
    n = len(result)
    if not n:
        raise ValueError("a tree must contain at least one vertex")
    if any(isinstance(p, bool) or not isinstance(p, int) for p in result):
        raise TypeError("parents must be integer labels, with -1 for the root")
    if any(p < -1 or p >= n for p in result):
        raise ValueError("parent label is outside the tree")
    if result.count(-1) != 1:
        raise ValueError("a tree must have exactly one root")
    for vertex in range(n):
        seen = set()
        cursor = vertex
        while cursor != -1:
            if cursor in seen:
                raise ValueError("parent pointers contain a cycle")
            seen.add(cursor)
            cursor = result[cursor]
    return result


def _matrix(rows, n=None, usage=False):
    result = tuple(tuple(_exact(x) for x in row) for row in rows)
    if not result:
        raise ValueError("at least one sample is required")
    width = len(result[0]) if n is None else n
    if not width or any(len(row) != width for row in result):
        raise ValueError("matrix rows must have the same positive width")
    if any(x < 0 or x > 1 for row in result for x in row):
        raise ValueError("entries must lie in [0, 1]")
    if usage and any(sum(row) > 1 for row in result):
        raise ValueError("usage rows must sum to at most one")
    return result


def rooted_trees(n):
    """Yield every rooted labelled n-vertex tree once using Prüfer codes.

    Prüfer sequences enumerate unrooted trees; orient each one away from
    each possible root.  The singleton case has no Prüfer sequence.
    """
    _positive_integer(n, "n")
    if n == 1:
        yield (-1,)
        return
    for code in product(range(n), repeat=n - 2):
        degree = [1] * n
        for vertex in code:
            degree[vertex] += 1
        adjacency = [[] for _ in range(n)]
        for vertex in code:
            leaf = next(i for i in range(n) if degree[i] == 1)
            adjacency[leaf].append(vertex)
            adjacency[vertex].append(leaf)
            degree[leaf] -= 1
            degree[vertex] -= 1
        a, b = (i for i in range(n) if degree[i] == 1)
        adjacency[a].append(b)
        adjacency[b].append(a)
        for root in range(n):
            parents = [-2] * n
            parents[root] = -1
            pending = [root]
            while pending:
                vertex = pending.pop()
                for child in adjacency[vertex]:
                    if parents[child] == -2:
                        parents[child] = vertex
                        pending.append(child)
            yield tuple(parents)


def children(parents):
    """Return label-sorted child tuples, validating the parent pointers."""
    parents = _parents(parents)
    result = [[] for _ in parents]
    for vertex, parent in enumerate(parents):
        if parent != -1:
            result[parent].append(vertex)
    return tuple(tuple(row) for row in result)


def ancestry_matrix(parents):
    """Return B, with B[i][j] = 1 iff j is an ancestor of i or j == i."""
    parents = _parents(parents)
    result = []
    for vertex in range(len(parents)):
        row = [0] * len(parents)
        cursor = vertex
        while cursor != -1:
            row[cursor] = 1
            cursor = parents[cursor]
        result.append(tuple(row))
    return tuple(result)


def frequencies(parents, usage_rows):
    """Compute F = U B exactly for valid nonnegative usage rows of sum <=1."""
    parents = _parents(parents)
    usage_rows = _matrix(usage_rows, len(parents), usage=True)
    result = []
    for usages in usage_rows:
        row = [Fraction(0) for _ in parents]
        for vertex, usage in enumerate(usages):
            cursor = vertex
            while cursor != -1:
                row[cursor] += usage
                cursor = parents[cursor]
        result.append(tuple(row))
    return tuple(result)


def _integer_rows(rows):
    """Clear a positive common denominator separately in each exact row."""
    result = []
    for row in rows:
        denominator = 1
        for value in row:
            denominator = denominator * value.denominator // gcd(
                denominator, value.denominator
            )
        result.append(tuple(
            value.numerator * (denominator // value.denominator)
            for value in row
        ))
    return tuple(result)


def _fits(child_rows, integer_rows):
    for row in integer_rows:
        for vertex, child_row in enumerate(child_rows):
            total = 0
            for child in child_row:
                total += row[child]
            if total > row[vertex]:
                return False
    return True


@lru_cache(maxsize=None)
def _catalog(n):
    """Cache child indices to avoid rebuilding trees in exhaustive checks."""
    return tuple((parents, children(parents)) for parents in rooted_trees(n))


def consistent(parents, F):
    """Test the exact sum condition; malformed inputs raise an exception.

    For a valid single-root tree the sum of inferred usages telescopes to
    its root frequency, already at most one by input validation.
    """
    parents = _parents(parents)
    rows = _matrix(F, len(parents))
    return _fits(children(parents), _integer_rows(rows))


def consistent_trees(F):
    """Return every consistent parent tuple by exhaustive tree enumeration."""
    rows = _matrix(F)
    integer_rows = _integer_rows(rows)
    return tuple(
        parents
        for parents, child_rows in _catalog(len(rows[0]))
        if _fits(child_rows, integer_rows)
    )


def design_details(parents):
    """Construct identifying rational mixtures for the given target tree.

    Leaves are ordered by DFS with increasing-label children.  With L
    leaves and W=n+1, forward leaf weights are W*2**i and reverse weights
    are W*2**(L-1-i).  Every internal vertex has weight one.  Divide by
    D=2*(W*(2**L-1)+n-L), so every usage row sums to 1/2.

    A path needs only the first row.  A branching tree uses both rows.
    The returned radius is a conservative closed L-infinity perturbation
    bound for F, as proved in the companion README; it is not a statement
    about measurement noise in an actual sequencing experiment.
    """
    parents = _parents(parents)
    child_rows = children(parents)
    leaves = []
    pending = [parents.index(-1)]
    while pending:
        vertex = pending.pop()
        if child_rows[vertex]:
            pending.extend(reversed(child_rows[vertex]))
        else:
            leaves.append(vertex)
    n = len(parents)
    leaf_count = len(leaves)
    weight = n + 1
    total = weight * (2**leaf_count - 1) + n - leaf_count
    denominator = 2 * total
    usage_rows = []
    for reverse in range(1 if leaf_count == 1 else 2):
        numerators = [1] * n
        for position, leaf in enumerate(leaves):
            exponent = leaf_count - 1 - position if reverse else position
            numerators[leaf] = weight * 2**exponent
        usage_rows.append(tuple(Fraction(x, denominator) for x in numerators))
    U = tuple(usage_rows)
    return {
        "F": frequencies(parents, U),
        "U": U,
        "denominator": denominator,
        "radius": Fraction(1, 4 * n * denominator),
        "leaves": tuple(leaves),
    }


def design_samples(parents):
    """Return the one- or two-row identifying frequency matrix."""
    return design_details(parents)["F"]


def max_ambiguity_witness(n, m):
    """Give m distinct rational rows with (n-1)! consistent trees.

    All columns are distinct, every compatible tree has strictly positive
    usages, and the row rank is min(m,n).  Hence repeated samples or rank
    deficiency are unnecessary for maximal ambiguity in this class.
    This witness does not purport to maximize ambiguity when equal
    columns or zero columns are admitted; the zero matrix admits all trees.
    """
    _positive_integer(n, "n")
    _positive_integer(m, "m")
    denominator = 2 * (2**n - 1)
    base = tuple(Fraction(2**(n - 1 - j), denominator) for j in range(n))
    epsilon = Fraction(1, 4 * n * denominator)
    result = []
    for sample in range(m):
        row = list(base)
        if sample < n:
            row[sample] += epsilon
        else:
            row[0] += epsilon / (sample + 1)
        result.append(tuple(row))
    return tuple(result)


def one_sample_unique(row):
    """Apply the exact one-sample uniqueness criterion, including ties.

    After sorting f[0]>...>f[n-1], the answer is true exactly when
    f[i] < f[i+1]+f[n-1] for i=0,...,n-3.  Tied entries give at least
    two labelled chain explanations, except for the singleton case.
    """
    values = _matrix((row,))[0]
    if len(values) == 1:
        return True
    if len(set(values)) != len(values):
        return False
    ordered = sorted(values, reverse=True)
    return all(
        ordered[i] < ordered[i + 1] + ordered[-1]
        for i in range(len(ordered) - 2)
    )


def enumeration_table(max_n=5):
    """Enumerate all target designs and maximal-ambiguity witness trees.

    Each target is checked against the entire candidate tree catalog.
    Counts are computed, not read from a hardcoded expected table.
    The zero matrix supplies the unrestricted maximum witness; the README
    proves the distinct-column maximum that the other witnesses attain.
    """
    _positive_integer(max_n, "max_n")
    result = []
    for n in range(1, max_n + 1):
        catalog = _catalog(n)
        paths = sum(all(len(row) <= 1 for row in kids) for _, kids in catalog)
        unique_designs = 0
        for parents, _ in catalog:
            if consistent_trees(design_samples(parents)) == (parents,):
                unique_designs += 1
        maxima = [len(consistent_trees(max_ambiguity_witness(n, m)))
                  for m in (1, 2, 3)]
        result.append({
            "n": n,
            "total_trees": len(catalog),
            "paths": paths,
            "branching": len(catalog) - paths,
            "max_unrestricted_m1": len(consistent_trees(((Fraction(0),) * n,))),
            "max_separated_m1": maxima[0],
            "max_separated_m2": maxima[1],
            "max_separated_m3": maxima[2],
            "unique_designs": unique_designs,
        })
    return result


def path_partition_count(nonfounders):
    """Count unordered collections of ordered paths on labelled vertices.

    For q nonfounders of equal positive frequency a and a founder of
    frequency strictly greater than a and at least q*a, these are exactly
    the consistent trees.  A nonfounder cannot have two children; the
    founder's children start an unordered family of disjoint ordered paths.
    For k paths, q! orders and binomial(q-1,k-1) cuts count each family k!
    times.  This evaluates that proved formula; it does not enumerate trees
    on q+1 vertices.  The empty family gives the q=0 singleton case.
    """
    if isinstance(nonfounders, bool) or not isinstance(nonfounders, int):
        raise TypeError("nonfounders must be an integer")
    if nonfounders < 0:
        raise ValueError("nonfounders must be nonnegative")
    if nonfounders == 0:
        return 1
    return sum(
        factorial(nonfounders) * comb(nonfounders - 1, k - 1) // factorial(k)
        for k in range(1, nonfounders + 1)
    )


def counterexample_checks():
    """Check explicit synthetic failures of tempting stronger statements."""
    missing_hypothesis_row = (
        Fraction(9, 10), Fraction(1, 2), Fraction(1, 5), Fraction(3, 10)
    )
    missing_hypothesis_F = (missing_hypothesis_row,)
    original_root_usage = (
        missing_hypothesis_row[0]
        - missing_hypothesis_row[1]
        - missing_hypothesis_row[2]
    )
    checks = {
        "a full-rank two-sample matrix can remain ambiguous": (
            len(consistent_trees(max_ambiguity_witness(3, 2))) == 2
        ),
        "three samples for three clones need not imply uniqueness": (
            len(consistent_trees(max_ambiguity_witness(3, 3))) == 2
        ),
        "equality at the one-sample threshold is not unique": (
            len(consistent_trees(((Fraction(9, 10), Fraction(3, 5),
                                   Fraction(3, 10)),))) == 2
        ),
        "distinct single-sample frequencies do not ensure uniqueness": (
            len(consistent_trees(((Fraction(9, 10), Fraction(1, 2),
                                   Fraction(1, 5)),))) == 2
        ),
        "equal columns allow distinct labelled chains": (
            len(consistent_trees(((Fraction(1, 2), Fraction(1, 2)),))) == 2
        ),
        "zero frequencies allow every three-vertex rooted tree": (
            len(consistent_trees(((Fraction(0),) * 3,))) == 9
        ),
        "the ancestor-skip test alone needs its dominance-graph hypothesis": (
            consistent((-1, 0, 0, 1), missing_hypothesis_F)
            and missing_hypothesis_row[3] > original_root_usage
            and consistent((-1, 0, 3, 1), missing_hypothesis_F)
        ),
        "tied-positive path counts match exhaustive enumeration for n=2..5": (
            all(
                len(consistent_trees((
                    (Fraction(1, 2),) + (Fraction(1, 32),) * (n - 1),
                ))) == path_partition_count(n - 1)
                for n in range(2, 6)
            )
        ),
        "the proved tied-positive count is 394353 > 9! (not n=9 enumeration)": (
            path_partition_count(8) == 394353 > factorial(9)
        ),
    }
    for description, passed in checks.items():
        if not passed:
            raise AssertionError(description)
    return tuple(checks)


def main(argv=None):
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=5,
                        help="largest exhaustive size (default: 5)")
    args = parser.parse_args(argv)
    rows = enumeration_table(args.max_n)
    fields = (
        "n", "total_trees", "paths", "branching", "max_unrestricted_m1",
        "max_separated_m1", "max_separated_m2", "max_separated_m3",
        "unique_designs",
    )
    print("Nonclinical computational mathematics; synthetic exact inputs only.")
    print()
    print("| " + " | ".join(fields) + " |")
    print("| " + " | ".join("---" for _ in fields) + " |")
    for row in rows:
        print("| " + " | ".join(str(row[field]) for field in fields) + " |")
        if row["unique_designs"] != row["total_trees"]:
            raise AssertionError("an identifying construction failed")
    print()
    for description in counterexample_checks():
        print("PASS: " + description)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
