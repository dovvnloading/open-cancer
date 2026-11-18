"""Offline, stdlib-only falsification tests for the frozen nonclinical receptor audit."""
from __future__ import annotations

import hashlib
import importlib.util
import math
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path, PureWindowsPath

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'research' / 'breast-receptor-age'
MODULE_PATH = BASE / 'analysis.py'
spec = importlib.util.spec_from_file_location('breast_receptor_age_analysis', MODULE_PATH)
analysis = importlib.util.module_from_spec(spec)
assert spec is not None and spec.loader is not None
spec.loader.exec_module(analysis)


class FrozenReceptorAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = analysis.load_data()
        cls.result = analysis.analyze(cls.data)

    def test_source_html_transcription_all_20_numeric_triplets(self):
        """Independent manual HTML re-read fixture, NOT rate-derived/imputed counts."""
        case_counts = {
            'overall': [37127, 29147, 20575, 13433, 12263],
            'positive': [31317, 24664, 17275, 11117, 9625],
            'negative': [4595, 3475, 2411, 1667, 1533],
            'unknown': [1197, 1002, 887, 644, 1095],
        }
        published = {
            'overall': [('417.7', '413.4', '421.9'), ('458.3', '453.0', '463.6'),
                        ('442.7', '436.7', '448.8'), ('393.5', '386.9', '400.2'),
                        ('295.2', '290.0', '300.5')],
            'positive': [('352.3', '348.4', '356.2'), ('387.8', '383.0', '392.7'),
                         ('371.7', '366.2', '377.3'), ('325.6', '319.6', '331.8'),
                         ('231.7', '227.1', '236.4')],
            'negative': [('51.7', '50.2', '53.2'), ('54.6', '52.8', '56.5'),
                         ('51.9', '49.8', '54.0'), ('48.8', '46.5', '51.2'),
                         ('36.9', '35.1', '38.8')],
            'unknown': [('13.7', '12.9', '14.5'), ('15.8', '14.9', '16.9'),
                        ('19.1', '17.9', '20.4'), ('19.0', '17.6', '20.5'),
                        ('26.6', '25.1', '28.2')],
        }
        populations = [8889353, 6360203, 4647676, 3413802, 4154029]
        self.assertEqual(len(case_counts), 4)
        for name in analysis.CLASSES:
            for i, age in enumerate(analysis.AGES):
                with self.subTest(name=name, age=age):
                    row = self.data[age][name]
                    self.assertEqual(row['cases'], case_counts[name][i])
                    self.assertEqual(row['population'], populations[i])
                    self.assertEqual((f'{row["rate"]:.1f}', f'{row["ci_low"]:.1f}',
                                      f'{row["ci_high"]:.1f}'), published[name][i])

    def test_literal_source_provenance_and_no_her2_substitution(self):
        raw = analysis.DATA.read_text(encoding='utf-8')
        self.assertEqual(raw.count('https://pmc.ncbi.nlm.nih.gov/articles/PMC10846538/#T1'), 20)
        self.assertEqual(raw.count('10.1007/s10549-021-06160-4'), 20)
        self.assertEqual(raw.count('November 2018 submission'), 20)
        self.assertEqual(raw.count('no explicit CC licence shown'), 20)
        self.assertNotIn('HER2,', raw)
        self.assertNotIn('ERBB2,', raw)
        self.assertEqual(set(analysis.CLASSES), {'overall', 'positive', 'negative', 'unknown'})

    def test_residuals_and_integer_proportional_rounding(self):
        expected_residuals = [18, 6, 2, 5, 10]
        for age, residual in zip(analysis.AGES, expected_residuals):
            z = self.data[age]
            self.assertEqual(z['residual'], residual)
            self.assertEqual(z['pool'], z['unknown']['cases'] + residual)
            p, n = z['positive']['cases'], z['negative']['cases']
            expected_k = math.floor(Fraction(z['pool'] * p, p + n) + Fraction(1, 2))
            self.assertEqual(z['prop_k'], expected_k)
            x, y = analysis.completed(z, 'proportional')
            self.assertEqual(x + y, z['overall']['cases'])
            self.assertTrue(p <= x <= p + z['pool'])
            self.assertTrue(n <= y <= n + z['pool'])

    def test_rate_reconciliation_retains_html_conflicts(self):
        mismatches = []
        for age in analysis.AGES:
            for klass in analysis.CLASSES:
                z = self.data[age][klass]
                diff = (100000 * z['cases'] / z['population']) - z['rate']
                if abs(diff) > 0.050000000001:
                    mismatches.append((age, klass))
        self.assertEqual(mismatches, [('65-69', 'unknown'), ('80-84', 'unknown'), ('85+', 'unknown')])
        text = analysis.render_report(self.result)
        self.assertEqual(text.count('**NO**'), len(mismatches))
        self.assertIn('source-known', text)

    def test_binomial_tail_exact_inversions_and_boundaries(self):
        t = analysis.TAIL
        self.assertEqual(t, 0.005)
        self.assertEqual(analysis.cp_interval(0, 10)[0], 0.)
        self.assertEqual(analysis.cp_interval(10, 10)[1], 1.)
        self.assertAlmostEqual(analysis.cp_interval(1, 1)[0], t, places=11)
        self.assertAlmostEqual(analysis.cp_interval(0, 1)[1], 1 - t, places=11)
        lo, hi = analysis.cp_interval(1, 2)
        self.assertAlmostEqual(lo, 1 - math.sqrt(1 - t), places=11)
        self.assertAlmostEqual(hi, math.sqrt(1 - t), places=11)
        for n, x in ((10, 3), (137, 71), (37127, 32377), (12263, 10578)):
            lo, hi = analysis.cp_interval(x, n)
            self.assertLess(lo, x / n)
            self.assertGreater(hi, x / n)
            self.assertAlmostEqual(analysis.regularized_beta(x, n - x + 1, lo), t, delta=1e-8)
            self.assertAlmostEqual(analysis.regularized_beta(x + 1, n - x, hi), 1 - t, delta=1e-8)
        with self.assertRaises(ValueError):
            analysis.cp_interval(11, 10)
        with self.assertRaises(ValueError):
            analysis.cp_interval(0, 0)

    def test_beta_independent_symmetry(self):
        for a, b, x in ((2, 5, .2), (7, 11, .6), (16, 31, .4), (3000, 450, .87)):
            got = analysis.regularized_beta(a, b, x)
            complement = analysis.regularized_beta(b, a, 1 - x)
            self.assertAlmostEqual(got + complement, 1, delta=2e-10)
        self.assertAlmostEqual(analysis.regularized_beta(1, 1, .234), .234, delta=1e-12)

    def test_multiplicity_and_all_ten_pairwise_coverage(self):
        self.assertEqual(len(self.result['pairs']), 10)
        self.assertEqual(analysis.TAIL * 2 * len(analysis.AGES), .05)
        self.assertEqual(sum(p['adjacent'] for p in self.result['pairs']), 4)
        for p in self.result['pairs']:
            self.assertLessEqual(p['envelope'][0], p['envelope'][1])
            for d in p['scenarios'].values():
                self.assertLessEqual(d['lo'], d['q'] + 1e-11)
                self.assertGreaterEqual(d['hi'] + 1e-11, d['q'])

    def test_adversarial_allocation_includes_every_scenario_interval(self):
        # Test stated monotonic bounds rather than brittle exact optimiser decimals.
        for p in self.result['pairs']:
            lo, hi = p['envelope']
            q_points = [d['q'] for d in p['scenarios'].values()]
            lo_b = p['scenarios']['younger_positive_older_negative']['q']
            hi_b = p['scenarios']['younger_negative_older_positive']['q']
            self.assertTrue(all(lo_b - 1e-12 <= q <= hi_b + 1e-12 for q in q_points))
            for d in p['scenarios'].values():
                self.assertLessEqual(lo, d['lo'] + 2e-10)
                self.assertGreaterEqual(hi + 2e-10, d['hi'])
            self.assertEqual(analysis.classify_interval(lo, hi), p['robust'])

    def test_nominal_and_robust_verdict_diverge_without_new_direction_choice(self):
        self.assertEqual(self.result['verdict'], 'INCONCLUSIVE')
        self.assertIn('Minimum completion to restore compatibility: NOT APPLICABLE',
                      analysis.render_report(self.result))
        self.assertIn('does **not** repair source truth.', analysis.render_report(self.result))
        self.assertTrue(all(p['robust'] == 'contains 1' for p in self.result['pairs']))
        nominal_exclusions = {(p['younger'], p['older']) for p in self.result['pairs']
                              if p['scenarios']['proportional']['direction'] != 'contains 1'}
        self.assertEqual(nominal_exclusions, {('70-74', '85+'), ('75-79', '85+')})
        for p in self.result['pairs']:
            d = p['scenarios']['proportional']
            self.assertEqual(analysis.classify_interval(d['lo'], d['hi']), d['direction'])

    def test_rounding_sensitivity_exact_bounds(self):
        for p in self.result['pairs']:
            a = self.data[p['younger']]
            b = self.data[p['older']]
            for half, rng in ((.05, p['rounding_half']), (.10, p['rounding_full'])):
                recomputed = analysis.rounding_q(a, b, half)
                for measured, expected in zip(rng, recomputed):
                    self.assertAlmostEqual(measured, expected, places=13)
                self.assertLess(rng[0], p['printed_q'])
                self.assertGreater(rng[1], p['printed_q'])
            self.assertLessEqual(p['rounding_full'][0], p['rounding_half'][0])
            self.assertGreaterEqual(p['rounding_full'][1], p['rounding_half'][1])

    def test_missing_or_invalid_input_abstains(self):
        with tempfile.TemporaryDirectory() as work:
            path = Path(work) / 'data.csv'
            original = analysis.DATA.read_text(encoding='utf-8')
            path.write_text(original.replace('65-69,overall,37127,', '65-69,overall,1,'), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'UNTESTABLE'):
                analysis.load_data(path)
            path.write_text(original.replace('65-69,positive,31317,', '65-69,positive,x,'), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'Non-numeric'):
                analysis.load_data(path)
            path.write_text(original.replace('65-69,positive,31317,', '65-69,positive,31317,' + '\n65-69,positive,31317,'), encoding='utf-8')
            with self.assertRaises(ValueError):
                analysis.load_data(path)
            path.write_text('\n'.join(original.splitlines()[:-1]) + '\n', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'Missing'):
                analysis.load_data(path)

    def test_newline_normalized_sha_and_windows_posix_paths(self):
        body = 'Receptor + — UTF-8\nsecond line\n'
        expected = hashlib.sha256(body.encode('utf-8')).hexdigest()
        with tempfile.TemporaryDirectory() as work:
            dest = Path(work) / 'case.txt'
            for newline in ('\n', '\r\n', '\r'):
                dest.write_bytes(body.replace('\n', newline).encode('utf-8'))
                self.assertEqual(analysis.sha256_text(dest), expected)
        self.assertEqual(PureWindowsPath('research\\breast-receptor-age\\PROTOCOL.md').as_posix(),
                         'research/breast-receptor-age/PROTOCOL.md')
        self.assertEqual(Path('research', 'breast-receptor-age', 'PROTOCOL.md').as_posix(),
                         'research/breast-receptor-age/PROTOCOL.md')
        self.assertEqual(analysis.key(analysis.PROTOCOL),
                         'research/breast-receptor-age/PROTOCOL.md')

    def test_frozen_hashes_and_exact_generated_markdown(self):
        report = analysis.render_report(self.result)
        readme = analysis.render_readme(self.result)
        self.assertEqual(analysis.normalize_newlines(analysis.REPORT.read_text(encoding='utf-8')), report)
        self.assertEqual(analysis.normalize_newlines(analysis.README.read_text(encoding='utf-8')), readme)
        for path in (analysis.DATA, analysis.PROTOCOL):
            self.assertIn(f'| {analysis.key(path)} | `{analysis.sha256_text(path)}` |', report)
        self.assertEqual(analysis.sha256_text(analysis.DATA),
                         'ad6c7f8888dde1649c4a926f05ab405dd26bd14315775e23cbd6ade8e21b5595')
        self.assertIn('ec33191bec49a945149ab1fc998ed588a41715d3', readme)
        self.assertIn('INCONCLUSIVE', report)
        self.assertIn('UNTESTABLE', report)
        self.assertIn('Not medical advice', report)


if __name__ == '__main__':
    unittest.main()
