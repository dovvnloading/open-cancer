"""Nonclinical computational mathematics; synthetic inputs and no patient data.

The oracle here does not use the research verifier's tree generator, ancestry
matrix, or sum-condition implementation. It enumerates parent maps, constructs
ancestry matrices by walking ancestors, and inverts them by exact Gauss-Jordan
elimination. Denominator clearing only accelerates these exact matrix products.
"""

from fractions import Fraction
from functools import lru_cache
import importlib.util
from itertools import combinations, product
from math import comb, factorial, gcd
from pathlib import Path
import unittest


_VERIFIER_PATH = (
    Path(__file__).resolve().parents[1]
    / "research"
    / "clonal-tree-identifiability"
    / "verifier.py"
)
_SPEC = importlib.util.spec_from_file_location(
    "clonal_tree_identifiability_verifier", _VERIFIER_PATH
)
verifier = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(verifier)


def independent_parent_maps(n):
    """Enumerate root choices and all parent maps, rejecting directed cycles."""
    for root in range(n):
        nonroots = tuple(vertex for vertex in range(n) if vertex != root)
        choices = tuple(
            tuple(parent for parent in range(n) if parent != vertex)
            for vertex in nonroots
        )
        for selected in product(*choices):
            parents = [-1] * n
            for vertex, parent in zip(nonroots, selected):
                parents[vertex] = parent
            valid = True
            for start in nonroots:
                seen = set()
                vertex = start
                while vertex != root:
                    if vertex in seen:
                        valid = False
                        break
                    seen.add(vertex)
                    vertex = parents[vertex]
                if not valid:
                    break
            if valid:
                yield tuple(parents)


def independent_ancestry(parents):
    """Rows are clones; columns are mutation labels on the ancestral path."""
    n = len(parents)
    matrix = []
    for clone in range(n):
        row = [0] * n
        vertex = clone
        while vertex != -1:
            row[vertex] = 1
            vertex = parents[vertex]
        matrix.append(tuple(row))
    return tuple(matrix)


def rational_inverse(matrix):
    """Generic square-matrix inverse, without any tree-specific identity."""
    n = len(matrix)
    augmented = [
        [Fraction(value) for value in row]
        + [Fraction(i == j) for j in range(n)]
        for i, row in enumerate(matrix)
    ]
    for col in range(n):
        pivot = next(
            (row for row in range(col, n) if augmented[row][col]), None
        )
        if pivot is None:
            raise ValueError("singular matrix")
        augmented[col], augmented[pivot] = augmented[pivot], augmented[col]
        scale = augmented[col][col]
        augmented[col] = [value / scale for value in augmented[col]]
        for row in range(n):
            if row != col and augmented[row][col]:
                scale = augmented[row][col]
                augmented[row] = [
                    left - scale * right
                    for left, right in zip(augmented[row], augmented[col])
                ]
    return tuple(tuple(row[n:]) for row in augmented)


@lru_cache(maxsize=None)
def independent_candidates(n):
    """Cache sparse columns of independently computed integral inverses."""
    candidates = []
    for parents in independent_parent_maps(n):
        inverse = rational_inverse(independent_ancestry(parents))
        if any(value.denominator != 1 for row in inverse for value in row):
            raise AssertionError("ancestry inverse unexpectedly nonintegral")
        columns = tuple(
            tuple(
                (row, int(inverse[row][col]))
                for row in range(n)
                if inverse[row][col]
            )
            for col in range(n)
        )
        candidates.append((parents, columns))
    return tuple(candidates)


def integer_rows(frequencies):
    """Represent each rational row exactly by its numerators and denominator."""
    encoded = []
    for row in frequencies:
        row = tuple(Fraction(value) for value in row)
        denominator = 1
        for value in row:
            denominator = (
                denominator * value.denominator
                // gcd(denominator, value.denominator)
            )
        encoded.append((
            tuple(int(value * denominator) for value in row), denominator
        ))
    return tuple(encoded)


def usage_numerators(row, inverse_columns):
    return tuple(
        sum(row[index] * coefficient for index, coefficient in column)
        for column in inverse_columns
    )


def independent_solutions(frequencies):
    """Test F B^-1 >= 0 and row sums <= 1 by integer matrix products."""
    rows = integer_rows(frequencies)
    solutions = []
    for parents, inverse_columns in independent_candidates(len(rows[0][0])):
        for row, denominator in rows:
            usage = usage_numerators(row, inverse_columns)
            if min(usage) < 0 or sum(usage) > denominator:
                break
        else:
            solutions.append(parents)
    return tuple(solutions)


def independent_usages(frequencies, inverse_columns):
    return tuple(
        tuple(Fraction(value, denominator) for value in
              usage_numerators(row, inverse_columns))
        for row, denominator in integer_rows(frequencies)
    )


class ClonalTreeIdentifiabilityTests(unittest.TestCase):
    """Finite, reproducible checks of mathematical claims, not clinical tests."""

    EXPECTED = (
        # n, labelled rooted trees, rooted paths, branching trees, separated max
        (1, 1, 1, 0, 1),
        (2, 2, 2, 0, 1),
        (3, 9, 6, 3, 2),
        (4, 64, 24, 40, 6),
        (5, 625, 120, 505, 24),
    )

    def test_enumeration_table_and_every_designed_tree_independently(self):
        reported = {row["n"]: row for row in verifier.enumeration_table(5)}
        for n, total, paths, branching, separated_max in self.EXPECTED:
            candidates = independent_candidates(n)
            parents_set = {parents for parents, _ in candidates}
            self.assertEqual(len(candidates), total)
            self.assertEqual(len(parents_set), total)
            self.assertEqual(set(verifier.rooted_trees(n)), parents_set)
            path_count = 0
            unique_designs = 0
            for parents, inverse_columns in candidates:
                is_path = all(parents.count(vertex) <= 1 for vertex in range(n))
                path_count += is_path
                details = verifier.design_details(parents)
                frequencies = verifier.design_samples(parents)
                self.assertEqual(details["F"], frequencies)
                self.assertEqual(len(frequencies), 1 if is_path else 2)
                self.assertEqual(
                    verifier.ancestry_matrix(parents),
                    independent_ancestry(parents),
                )
                usage = independent_usages(frequencies, inverse_columns)
                self.assertEqual(usage, details["U"])
                self.assertTrue(all(value > 0 for row in usage for value in row))
                self.assertTrue(all(sum(row) == Fraction(1, 2) for row in usage))
                solutions = independent_solutions(frequencies)
                self.assertEqual(solutions, (parents,), (parents, frequencies))
                unique_designs += 1
            self.assertEqual(path_count, paths)
            self.assertEqual(total - path_count, branching)
            self.assertEqual(unique_designs, total)
            self.assertEqual(reported[n]["total_trees"], total)
            self.assertEqual(reported[n]["paths"], paths)
            self.assertEqual(reported[n]["branching"], branching)
            self.assertEqual(reported[n]["unique_designs"], total)
            for m in (1, 2, 3):
                witness = verifier.max_ambiguity_witness(n, m)
                self.assertEqual(len(witness), m)
                self.assertEqual(len(independent_solutions(witness)), separated_max)
                self.assertEqual(
                    reported[n]["max_separated_m" + str(m)], separated_max
                )
            zeros = (tuple(Fraction(0) for _ in range(n)),)
            self.assertEqual(len(independent_solutions(zeros)), total)

    def test_complete_one_sample_criterion_on_exact_finite_grid(self):
        tested = 0
        for n in range(1, 6):
            for numerators in combinations(range(1, 12), n):
                row = tuple(Fraction(value, 12) for value in reversed(numerators))
                unique = len(independent_solutions((row,))) == 1
                self.assertEqual(verifier.one_sample_unique(row), unique, row)
                tested += 1
        self.assertEqual(tested, 1023)

    def test_sum_boundary_ties_and_repeated_samples_are_counterexamples(self):
        for value in (Fraction(0), Fraction(1)):
            self.assertEqual(independent_solutions(((value,),)), ((-1,),))
            self.assertTrue(verifier.one_sample_unique((value,)))
        self.assertEqual(independent_solutions(((Fraction(1), Fraction(0)),)),
                         ((-1, 0),))
        self.assertTrue(verifier.one_sample_unique((Fraction(1), Fraction(0))))
        self.assertFalse(verifier.one_sample_unique(
            (Fraction(1), Fraction(1, 2), Fraction(0))))
        boundary = (Fraction(3, 4), Fraction(1, 2), Fraction(1, 4))
        self.assertEqual(
            set(independent_solutions((boundary,))),
            {(-1, 0, 1), (-1, 0, 0)},
        )
        self.assertFalse(verifier.one_sample_unique(boundary))
        below = (Fraction(2, 3), Fraction(1, 2), Fraction(1, 4))
        self.assertEqual(independent_solutions((below,)), ((-1, 0, 1),))
        self.assertTrue(verifier.one_sample_unique(below))
        tied = ((Fraction(1, 2),) * 3,)
        self.assertEqual(len(independent_solutions(tied)), 6)
        repeated = ((Fraction(3, 4), Fraction(1, 4), Fraction(1, 8)),) * 7
        self.assertEqual(len(independent_solutions(repeated)), 2)
        self.assertEqual(
            set(verifier.consistent_trees(repeated)),
            set(independent_solutions(repeated)),
        )

    def test_full_rank_positive_usage_does_not_imply_identifiability(self):
        frequencies = tuple(
            tuple(Fraction(value + (i == j), 128)
                  for j, value in enumerate((64, 16, 4)))
            for i in range(3)
        )
        inverse = rational_inverse(frequencies)  # Existence establishes rank 3.
        for i in range(3):
            for j in range(3):
                self.assertEqual(
                    sum(frequencies[i][k] * inverse[k][j] for k in range(3)),
                    int(i == j),
                )
        solutions = independent_solutions(frequencies)
        self.assertEqual(set(solutions), {(-1, 0, 0), (-1, 0, 1)})
        inverse_by_tree = dict(independent_candidates(3))
        for tree in solutions:
            usage = independent_usages(frequencies, inverse_by_tree[tree])
            self.assertTrue(all(value > 0 for row in usage for value in row))
            self.assertTrue(all(sum(row) < 1 for row in usage))
        self.assertEqual(set(verifier.consistent_trees(frequencies)), set(solutions))

    def test_exact_model_boundaries_are_enforced(self):
        # A rounded float must not silently replace an exact boundary value.
        with self.assertRaises(TypeError):
            verifier.consistent((-1, 0), ((0.75, 0.25),))
        with self.assertRaises(ValueError):
            verifier.frequencies((-1, 0), ((Fraction(3, 4), Fraction(3, 4)),))
        with self.assertRaises(ValueError):
            verifier.consistent((-1, 1), ((Fraction(1), Fraction(1, 2)),))
        with self.assertRaises(ValueError):
            verifier.consistent((-1, -1), ((Fraction(1), Fraction(1, 2)),))
        with self.assertRaises(ValueError):
            verifier.consistent((-1,), ())
        with self.assertRaises(ValueError):
            tuple(verifier.rooted_trees(0))

    def test_exact_robustness_certificates_against_all_trees_through_four(self):
        for n in range(1, 5):
            candidates = independent_candidates(n)
            for target, inverse_columns in candidates:
                details = verifier.design_details(target)
                frequencies = details["F"]
                radius = details["radius"]
                self.assertEqual(radius, Fraction(1, 4 * n * details["denominator"]))
                true_usage = independent_usages(frequencies, inverse_columns)
                self.assertGreater(min(value for row in true_usage for value in row),
                                   n * radius)
                self.assertTrue(all(1 - sum(row) > n * radius for row in true_usage))
                for row in frequencies:
                    self.assertTrue(all(radius < value < 1 - radius for value in row))
                    for left, right in combinations(row, 2):
                        self.assertGreater(abs(left - right), 2 * radius)
                for alternative, alternative_inverse in candidates:
                    if alternative == target:
                        continue
                    usage = independent_usages(frequencies, alternative_inverse)
                    worst_violation = max(-value for row in usage for value in row)
                    self.assertGreater(worst_violation, n * radius,
                                       (target, alternative, frequencies))

    def test_fixed_normal_root_translation_independently(self):
        # N includes the normal vertex 0; n = N - 1 counts observed mutations.
        for augmented_n in range(2, 6):
            n = augmented_n - 1
            normal_candidates = []
            for parents, _ in independent_candidates(augmented_n):
                if parents[0] != -1:
                    continue
                augmented_ancestry = independent_ancestry(parents)
                malignant_ancestry = tuple(row[1:] for row in augmented_ancestry[1:])
                inverse = rational_inverse(malignant_ancestry)
                self.assertTrue(all(value.denominator == 1
                                    for row in inverse for value in row))
                inverse_columns = tuple(
                    tuple((i, int(inverse[i][j])) for i in range(n) if inverse[i][j])
                    for j in range(n)
                )
                normal_candidates.append((parents, inverse_columns))
            self.assertEqual(len(normal_candidates), augmented_n ** (n - 1))
            self.assertEqual(
                sum(all(parents.count(v) <= 1 for v in range(augmented_n))
                    for parents, _ in normal_candidates),
                factorial(n),
            )

            def normal_solutions(observed):
                encoded = integer_rows(observed)
                solutions = []
                for parents, inverse_columns in normal_candidates:
                    for row, denominator in encoded:
                        usage = usage_numerators(row, inverse_columns)
                        if min(usage) < 0 or sum(usage) > denominator:
                            break
                    else:
                        solutions.append(parents)
                return tuple(solutions)

            for target, inverse_columns in normal_candidates:
                details = verifier.design_details(target)
                augmented_frequencies = tuple(
                    tuple(2 * value for value in row)
                    for row in verifier.design_samples(target)
                )
                self.assertTrue(all(row[0] == 1 for row in augmented_frequencies))
                observed = tuple(row[1:] for row in augmented_frequencies)
                malignant_usage = independent_usages(observed, inverse_columns)
                self.assertEqual(
                    malignant_usage,
                    tuple(tuple(2 * value for value in row[1:])
                          for row in details["U"]),
                )
                self.assertTrue(all(value > 0
                                    for row in malignant_usage for value in row))
                self.assertTrue(all(sum(row) < 1 for row in malignant_usage))
                self.assertEqual(normal_solutions(observed), (target,))
            for m in (1, 2, 3):
                self.assertEqual(
                    len(normal_solutions(verifier.max_ambiguity_witness(n, m))),
                    factorial(n),
                )

    def test_dominance_criterion_on_all_4096_two_sample_grid_matrices(self):
        candidates = independent_candidates(3)
        ancestry = {tree: independent_ancestry(tree) for tree, _ in candidates}
        closures = {
            tree: frozenset((a, v) for a in range(3) for v in range(3)
                            if a != v and matrix[v][a])
            for tree, matrix in ancestry.items()
        }
        inverse_by_tree = dict(candidates)
        tested = eligible = unique_count = ambiguous_count = 0
        for numerators in product(range(1, 5), repeat=6):
            frequencies = (
                tuple(Fraction(value, 5) for value in numerators[:3]),
                tuple(Fraction(value, 5) for value in numerators[3:]),
            )
            solutions = independent_solutions(frequencies)
            dominance = frozenset(
                (a, v) for a in range(3) for v in range(3)
                if a != v and all(row[a] >= row[v] for row in frequencies)
            )
            for target in solutions:
                if dominance != closures[target]:
                    continue
                eligible += 1
                usage = independent_usages(frequencies, inverse_by_tree[target])
                criterion = all(
                    any(frequencies[s][v] > usage[s][a] for s in range(2))
                    for a, v in closures[target] if target[v] != a
                )
                unique = len(solutions) == 1
                self.assertEqual(criterion, unique, (target, frequencies))
                unique_count += unique
                ambiguous_count += not unique
            tested += 1
        self.assertEqual((tested, eligible, unique_count, ambiguous_count),
                         (4096, 1350, 1152, 198))

    def test_positive_tied_frequencies_exceed_the_naive_factorial_bound(self):
        # The large value is a closed-form count, NOT an n=9 tree enumeration.
        count_eight = sum(
            factorial(8) * comb(7, paths - 1) // factorial(paths)
            for paths in range(1, 9)
        )
        self.assertEqual(count_eight, 394353)
        self.assertGreater(count_eight, factorial(9))
        self.assertEqual(verifier.path_partition_count(8), count_eight)
        self.assertEqual(verifier.path_partition_count(0), 1)
        for n, expected in ((2, 1), (3, 3), (4, 13), (5, 73)):
            row = (Fraction(1, 2),) + (Fraction(1, 32),) * (n - 1)
            self.assertEqual(len(independent_solutions((row,))), expected)
            self.assertEqual(verifier.path_partition_count(n - 1), expected)


if __name__ == "__main__":
    unittest.main()
