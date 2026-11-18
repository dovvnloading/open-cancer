"""Independent, offline checks of NONCLINICAL published-aggregate arithmetic.

Only the small aggregate CSVs are read. The reference calculations below do
not call analysis.py's statistical helpers: they use direct conditional
probabilities, Decimal inversion, score-test inversion and rational convex
bounds. Checking arithmetic/provenance completeness cannot verify a source
transcription, establish sampling assumptions or validate biological claims.
Standard library only; no individual records or network access.
"""

import csv
from decimal import Decimal, localcontext
from fractions import Fraction
import importlib.util
import io
import itertools
import math
from pathlib import Path
import re
import unittest
from unittest.mock import patch


RESEARCH = (Path(__file__).resolve().parents[1] / "research"
            / "prediction-confrontation")
CSV_FILES = ("immune.csv", "hla_response.csv", "ctdna.csv", "ecdna.csv",
             "resistance.csv")
CSV_COLUMNS = ("id", "study", "population", "quantity", "value", "unit",
               "source_url", "doi", "locator", "access_date", "reuse_basis",
               "notes")
F = Fraction
D = Decimal

spec = importlib.util.spec_from_file_location(
    "prediction_confrontation_analysis", RESEARCH / "analysis.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def normal_quantile(probability):
    """Invert the normal CDF independently of statistics.NormalDist."""
    lo, hi = -12.0, 12.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if math.erfc(-mid / math.sqrt(2)) / 2 < probability:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


Z = normal_quantile(.975)


def score_interval(k, n):
    """Wilson bounds by inverting the binomial score inequality."""
    observed = k / n

    def discrepancy(p):
        return (k - n*p)**2 - Z*Z*n*p*(1-p)

    lo, hi = 0.0, observed
    for _ in range(90):
        mid = (lo + hi) / 2
        if discrepancy(mid) > 0:
            lo = mid
        else:
            hi = mid
    lower = 0.0 if k == 0 else (lo + hi) / 2
    lo, hi = observed, 1.0
    for _ in range(90):
        mid = (lo + hi) / 2
        if discrepancy(mid) > 0:
            hi = mid
        else:
            lo = mid
    upper = 1.0 if k == n else (lo + hi) / 2
    return lower, upper


def conditional_weights(a, b, c, d):
    """Condition on the response total, the opposite margin to production."""
    n1, n0, successes = a+b, c+d, a+c
    return {x: math.comb(n1, x)*math.comb(n0, successes-x)
            for x in range(max(0, successes-n0), min(n1, successes)+1)}


def exact_fisher(a, b, c, d):
    weights = conditional_weights(a, b, c, d)
    return float(F(sum(w for w in weights.values() if w <= weights[a]),
                   sum(weights.values())))


def noncentral_tail(weights, observed, odds, upper):
    """Direct Decimal powers, without log-sum-exp or production helpers."""
    values = {x: D(w)*odds**x for x, w in weights.items()}
    numerator = sum(w for x, w in values.items()
                    if (x >= observed if upper else x <= observed))
    return numerator / sum(values.values())


def exact_odds_interval(a, b, c, d):
    weights = conditional_weights(a, b, c, d)
    with localcontext() as context:
        context.prec = 60
        bounds = []
        for upper_tail in (True, False):
            if upper_tail and a == min(weights):
                bounds.append(0.0)
                continue
            if not upper_tail and a == max(weights):
                bounds.append(math.inf)
                continue
            lo, hi = D(0), D(1)
            for _ in range(100):
                probability = noncentral_tail(weights, a, hi, upper_tail)
                if (probability >= D('.025') if upper_tail
                        else probability <= D('.025')):
                    break
                hi *= 2
            else:
                raise AssertionError("Failed to bracket a conditional bound")
            for _ in range(160):
                mid = (lo + hi) / 2
                probability = noncentral_tail(weights, a, mid, upper_tail)
                if (probability < D('.025')) == upper_tail:
                    lo = mid
                else:
                    hi = mid
            bounds.append(float((lo + hi) / 2))
    return tuple(bounds)


def independent_contrast(a, b, c, d):
    p1, p0 = F(a, a+b), F(c, c+d)
    ci1, ci0 = score_interval(a, a+b), score_interval(c, c+d)
    rd = float(p1-p0)
    rd_ci = (max(-1., rd-math.sqrt((float(p1)-ci1[0])**2
                                  + (ci0[1]-float(p0))**2)),
             min(1., rd+math.sqrt((ci1[1]-float(p1))**2
                                  + (float(p0)-ci0[0])**2)))
    rr = float(p1/p0)
    rr_se = math.sqrt(float(F(b, a*(a+b)) + F(d, c*(c+d))))
    return dict(counts=[a, b, c, d], or_=float(F(a*d, b*c)),
                or_ci=exact_odds_interval(a, b, c, d),
                p=exact_fisher(a, b, c, d), risk_loh=float(p1),
                risk_intact=float(p0), risk_loh_ci=ci1, risk_intact_ci=ci0,
                rd=rd, rd_ci=rd_ci, rr=rr,
                rr_ci=(math.exp(math.log(rr)-Z*rr_se),
                       math.exp(math.log(rr)+Z*rr_se)))


def independent_holm(pvalues):
    """Closed-form prefix maxima; valid without independence of tests."""
    keys = sorted(pvalues, key=lambda key: (pvalues[key], key))
    return {key: min(1., max((len(keys)-j)*pvalues[keys[j]]
                             for j in range(rank+1)))
            for rank, key in enumerate(keys)}


def rounding_quantum(text):
    coefficient, _, exponent = text.lower().partition('e')
    places = len(coefficient.partition('.')[2])
    return F(10)**(int(exponent or '0')-places)


def convex_thinning_bounds(mean, lower, upper, cadence):
    """Jensen and endpoint-chord bounds for g(L), in exact arithmetic."""
    def g(lead):
        return (lead + max(F(0), lead-cadence)) / 2
    chord = ((upper-mean)*g(lower) + (mean-lower)*g(upper))/(upper-lower)
    return g(mean), chord


def integer_rounding_envelope(n, percent):
    lo, hi = (percent-F(1, 2))*n/100, (percent+F(1, 2))*n/100
    return list(range(max(0, math.ceil(lo)), min(n, math.floor(hi))+1))


def binomial_tail_by_recurrence(k, n, p):
    """Exact rational probabilities by recurrence, not a combinatorial sum."""
    if p == 0:
        return F(k == 0)
    if p == 1:
        return F(k <= n)
    probability = (1-p)**n
    total = probability if k == 0 else F(0)
    for j in range(1, n+1):
        probability *= F(n-j+1, j)*p/(1-p)
        if j >= k:
            total += probability
    return total


def read_csvs():
    rows, texts, file_counts = {}, {}, {}
    for filename in CSV_FILES:
        text = (RESEARCH / "data" / filename).read_text(encoding="utf-8")
        texts[filename] = text
        reader = csv.DictReader(io.StringIO(text))
        if tuple(reader.fieldnames) != CSV_COLUMNS:
            raise AssertionError(f"Unexpected schema in {filename}")
        records = list(reader)
        file_counts[filename] = len(records)
        for row in records:
            if row['id'] in rows:
                raise AssertionError(f"Duplicate id {row['id']}")
            rows[row['id']] = dict(row, file=filename)
    return rows, texts, file_counts


def recompute(rows, file_counts):
    """Recompute the entire numeric result tree solely from CSV inputs."""
    def f(key):
        return F(rows[key]['value'])

    def v(key):
        return float(f(key))

    def integer(key):
        value = f(key)
        if value.denominator != 1:
            raise AssertionError(f"Expected aggregate integer for {key}")
        return int(value)

    effects, family = [], {}
    for prefix, kind in (("litchfield_clonal", "or"),
                         ("litchfield_total", "or"),
                         ("litchfield_subclonal", "or"),
                         ("litchfield_hla", "or"),
                         ("ravi_clonal", "or"),
                         ("lozac_luad", "hr"), ("lozac_lusc", "hr")):
        item = dict(id=prefix, effect=v(prefix+'_'+kind),
                    ci=(v(prefix+'_low'), v(prefix+'_high')))
        if prefix == 'ravi_clonal':
            se = math.log(item['ci'][1]/item['ci'][0])/(2*Z)
            z = math.log(item['effect'])/se
            item.update(se=se, z=z, p=math.erfc(abs(z)/math.sqrt(2)),
                        source_q=v(prefix+'_q'))
        else:
            item['p'] = v(prefix+'_p')
        if prefix.startswith('lozac_'):
            item['n'] = integer(prefix+'_n')
        if prefix not in ('litchfield_total', 'litchfield_subclonal'):
            family[prefix] = item['p']
        effects.append(item)

    yang = []
    for cohort in ('urothelial', 'melanoma_pd1', 'melanoma_ctla4'):
        for locus in ('a', 'b', 'c'):
            key = f'yang_2022_{cohort}_hla_{locus}'
            counts = [integer(key+'_'+suffix) for suffix in
                      ('response_loh', 'nonresponse_loh', 'response_intact',
                       'nonresponse_intact')]
            item = dict(id=key, **independent_contrast(*counts))
            family[key] = item['p']
            yang.append(item)
    q_upper = float(f('ravi_clonal_q')
                    + rounding_quantum(rows['ravi_clonal_q']['value'])/2)
    sensitivity = dict(family, ravi_clonal=q_upper)
    adjusted, adjusted_bound = (independent_holm(p)
                                for p in (family, sensitivity))
    for item in effects + yang:
        if item['id'] in adjusted:
            item.update(holm=adjusted[item['id']],
                        holm_q_bound=adjusted_bound[item['id']])

    m, a, b, d = [f('reinert_'+name) for name in
                  ('lead_mean', 'lead_min', 'lead_max', 'draw_interval')]
    forecast = convex_thinning_bounds(m, a, b, d)
    h = F(1, 20)
    corners = [convex_thinning_bounds(mm, aa, bb, d)
               for mm, aa, bb in itertools.product((m-h, m+h), (a-h, a+h),
                                                    (b-h, b+h))]
    ctdna = dict(
        mean=float(m), minimum=float(a), maximum=float(b), cadence=float(d),
        forecast=tuple(map(float, forecast)),
        loss=(float(m-forecast[1]), float(m-forecast[0])),
        rounding_envelope=(float(min(x[0] for x in corners)),
                           float(max(x[1] for x in corners))),
        width=float(b-a), width_rounding=(float(b-a-2*h), float(b-a+2*h)),
        recurrence_counts=[integer('reinert_detected_recurrences'),
                           integer('reinert_total_recurrences')],
        fast_log_of_mean=math.log(v('henriksen22_fast_mean')),
        slow_log_of_mean=math.log(v('henriksen22_slow_mean')),
        matched_counts=[integer('henriksen22_matched_early_n'),
                        integer('henriksen22_matched_n')],
        matched_median_iqr=[v('henriksen22_matched_lead_'+s)
                            for s in ('median', 'q1', 'q3')],
        extrapolated_counts=[integer('henriksen24_delay_n'),
                             integer('henriksen24_growth_n')],
        extrapolated_median_iqr=[v('henriksen24_gain_'+s)
                                 for s in ('median', 'q1', 'q3')],
        tracerx_summary=[v('tracerx23_lead_'+s) for s in ('median', 'min', 'max')]
                        + [integer('tracerx23_lead_n')])

    def allocation(key):
        p, rho = f(key+'_p'), f(key+'_rho')
        return dict(p=float(p), rho=float(rho),
                    R=float(rho+(1-rho)*(1-2*p)**2))

    clones, p9, p10 = [], {}, {}
    for i in range(1, 13):
        key = f'shibata_2026_clone{i}'
        item = dict(id=key, clone=i, **allocation(key),
                    ks_p=v(key+'_ks_p'), n=integer(key+'_n'))
        clones.append(item)
        if i not in (4, 6, 7):
            p9[key] = item['ks_p']
        if i not in (4, 6):
            p10[key] = item['ks_p']
    h9, h10 = independent_holm(p9), independent_holm(p10)
    for item in clones:
        item.update(holm9=h9.get(item['id']), holm10=h10.get(item['id']))
    arms = [dict(arm=arm, **allocation('shibata_2026_'+arm)) for arm in
            ('ctrl1_dmso', 'ctrl1_jq1', 'ctrlen2_dmso', 'ctrlen2_jq1')]
    lange = []
    for panel in ('pc3_myc', 'snu16_myc', 'snu16_fgfr2', 'gbm39_egfr',
                  'colo_myc', 'tr14_mycn', 'tr14_cdk4', 'colo320_dm'):
        key = 'lange_2022_'+panel
        lange.append(dict(panel=panel, D=v(key+'_ks_d'), p=v(key+'_ks_p'),
                          n=integer(key+'_n')))

    lundberg, conditional_family = [], {}
    alpha = f('lundberg_2008_alpha')
    for panel in ('chp212', 'la_n5', 'mcix'):
        key = 'lundberg_2008_'+panel
        n, percent = integer(key+'_n'), f(key+'_pct')
        counts = integer_rounding_envelope(n, percent)
        low, high = min(counts), max(counts)
        p = max(binomial_tail_by_recurrence(k, n, alpha) for k in counts)
        item = dict(id=key, panel=panel, n=n, percent=float(percent),
                    counts=counts, rate=(float(F(low, n)), float(F(high, n))),
                    ci=(score_interval(low, n)[0], score_interval(high, n)[1]),
                    p=float(p), excess=(float(F(low, n)-alpha),
                                        float(F(high, n)-alpha)),
                    fold=(float(F(low, n)/alpha), float(F(high, n)/alpha)),
                    markov=float(min(F(1), n*alpha/low)) if low else 1.)
        conditional_family[key] = item['p']
        lundberg.append(item)
    h3 = independent_holm(conditional_family)
    h17 = independent_holm(dict(family, **conditional_family))
    h17_bound = independent_holm(dict(sensitivity, **conditional_family))
    for item in lundberg:
        item.update(holm3=h3[item['id']], holm_combined=h17[item['id']],
                    holm_combined_q_bound=h17_bound[item['id']])

    persister, vmr = f('iyer_persister_fraction'), f('iyer_experimental_vmr')
    resistance = dict(persister=float(persister), vmr=float(vmr),
                      n=integer('iyer_lineage_count'),
                      fixed_size_binomial_vmr=float(1-persister),
                      excess=float(vmr-(1-persister)),
                      ratio=float(vmr/(1-persister)))
    return dict(row_count=len(rows), file_counts=file_counts, effects=effects,
                yang=yang, formal_family=family, formal_holm=adjusted,
                formal_holm_q_bound=adjusted_bound,
                ravi_q_rounding_upper=q_upper, ctdna=ctdna, clones=clones,
                arms=arms, lange=lange, lundberg=lundberg, resistance=resistance)


def numeric_paths(value, path=()):
    """A coverage guard: a new numeric output needs an independent check."""
    if isinstance(value, dict):
        out = {}
        for key, child in value.items():
            out.update(numeric_paths(child, path+(key,)))
        return out
    if isinstance(value, (list, tuple)):
        out = {}
        for index, child in enumerate(value):
            out.update(numeric_paths(child, path+(index,)))
        return out
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return {path: value}
    return {}


def display(value):
    return 'excluded' if value is None else format(value, '.6g')


def display_interval(values):
    return '['+', '.join(display(value) for value in values)+']'


class PredictionConfrontationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows, cls.csv_texts, cls.counts = read_csvs()
        cls.expected = recompute(cls.rows, cls.counts)
        cls.result = analysis.calculate()

    def assertClose(self, actual, expected):
        self.assertTrue(math.isclose(actual, expected, rel_tol=3e-12,
                                     abs_tol=3e-14), (actual, expected))

    def test_every_numeric_output_is_independently_recomputed(self):
        actual = numeric_paths(self.result)
        expected = numeric_paths(self.expected)
        self.assertEqual(set(actual), set(expected),
                         "Every new numeric result requires an independent calculation")
        for path, value in expected.items():
            with self.subTest(path=path):
                self.assertClose(actual[path], value)

    def test_csv_schema_provenance_and_source_units(self):
        self.assertEqual(analysis.load_data(), self.rows)
        self.assertEqual(set(self.counts), set(CSV_FILES))
        for key, row in self.rows.items():
            with self.subTest(row=key):
                self.assertEqual(set(row), set(CSV_COLUMNS) | {'file'})
                self.assertTrue(all(row[column].strip() for column in CSV_COLUMNS))
                self.assertTrue(math.isfinite(float(row['value'])))
                self.assertRegex(row['source_url'], r'^https://[^/]+/')
                self.assertRegex(row['doi'], r'^10\.\d{4,9}/\S+$')
                self.assertEqual(row['access_date'], '2026-11-17')
                self.assertTrue('CC BY' in row['reuse_basis']
                                or 'Creative Commons Attribution' in row['reuse_basis'])
                self.assertNotIn('patient_id', row)
                self.assertNotIn('cell_id', row)
                if row['file'] == 'hla_response.csv':
                    self.assertEqual(F(row['value']).denominator, 1)
                    self.assertGreaterEqual(int(row['value']), 0)
                if row['file'] == 'immune.csv' and key.endswith('_p'):
                    self.assertEqual(row['unit'], 'p_value')
        self.assertEqual(self.rows['ravi_clonal_q']['unit'], 'BH_adjusted_p_value')
        self.assertIn('95 percent CI', self.rows['ravi_clonal_low']['notes'])
        self.assertIn('not CI', self.rows['henriksen22_matched_lead_q1']['notes'])
        self.assertIn('no CI', self.rows['shibata_2026_clone1_p']['notes'])

    def test_source_effects_and_formal_family_are_not_silently_changed(self):
        actual = {item['id']: item for item in self.result['effects']}
        for expected in self.expected['effects']:
            item = actual[expected['id']]
            self.assertEqual(item['ci'], expected['ci'])
            self.assertEqual(item['effect'], expected['effect'])
            if item['id'] == 'ravi_clonal':
                self.assertIn('approximate', item['p_kind'])
                self.assertGreater(self.result['ravi_q_rounding_upper'], item['source_q'])
                self.assertNotEqual(item['p'], item['source_q'])
            else:
                self.assertEqual(item['p_kind'], 'source reported')
                self.assertEqual(item['p'], expected['p'])
        family = self.result['formal_family']
        self.assertEqual(len(family), 14)
        self.assertNotIn('litchfield_total', family)
        self.assertNotIn('litchfield_subclonal', family)
        self.assertEqual(sum(key.startswith('yang_') for key in family), 9)
        self.assertTrue(all(not key.startswith(('shibata_', 'lange_', 'lundberg_'))
                            for key in family))
        # The opposite nominal associations remain in the family and outputs.
        self.assertLess(actual['lozac_luad']['ci'][1], 1)
        self.assertGreater(actual['lozac_luad']['holm'], .05)
        urothelial = next(x for x in self.result['yang']
                          if x['id'] == 'yang_2022_urothelial_hla_a')
        self.assertGreater(urothelial['or_ci'][0], 1)
        self.assertLess(urothelial['p'], .05)
        self.assertGreater(urothelial['holm'], .05)

    def test_exact_or_interval_tails_all_nine_panels(self):
        for item in self.result['yang']:
            a, b, c, d = item['counts']
            weights = conditional_weights(a, b, c, d)
            self.assertEqual(item['p'], exact_fisher(a, b, c, d))
            with self.subTest(panel=item['id']), localcontext() as context:
                context.prec = 55
                for bound, upper in zip(item['or_ci'], (True, False)):
                    probability = noncentral_tail(weights, a, D(str(bound)), upper)
                    self.assertClose(float(probability), .025)
        self.assertEqual(analysis.conditional_or_interval(0, 4, 4, 0)[0], 0)
        self.assertTrue(math.isinf(analysis.conditional_or_interval(4, 0, 0, 4)[1]))

    def test_score_and_risk_intervals_are_on_their_declared_scales(self):
        for item in self.result['yang']:
            for p, bounds in ((item['risk_loh'], item['risk_loh_ci']),
                              (item['risk_intact'], item['risk_intact_ci'])):
                self.assertLessEqual(0, bounds[0])
                self.assertLessEqual(bounds[0], p)
                self.assertLessEqual(p, bounds[1])
                self.assertLessEqual(bounds[1], 1)
            self.assertLessEqual(-1, item['rd_ci'][0])
            self.assertLessEqual(item['rd_ci'][1], 1)
            self.assertGreater(item['rr_ci'][0], 0)
        for k, n in ((0, 10), (10, 10), (1, 2), (9, 17)):
            for actual, expected in zip(analysis.wilson(k, n), score_interval(k, n)):
                self.assertClose(actual, expected)

    def test_holm_keeps_dependent_contrasts_and_uses_joint_families(self):
        family = self.expected['formal_family']
        conditional = {x['id']: x['p'] for x in self.expected['lundberg']}
        self.assertEqual(len(dict(family, **conditional)), 17)
        self.assertEqual(analysis.holm(family), independent_holm(family))
        duplicate = {'same_panel_a': .02, 'same_panel_b': .02, 'third': .9}
        self.assertEqual(analysis.holm(duplicate), independent_holm(duplicate))
        self.assertEqual(analysis.holm(dict(reversed(list(duplicate.items())))),
                         analysis.holm(duplicate))
        for key, value in family.items():
            self.assertGreaterEqual(self.result['formal_holm'][key], value)
            self.assertGreaterEqual(self.result['formal_holm_q_bound'][key],
                                    self.result['formal_holm'][key])
        lan5 = next(x for x in self.result['lundberg'] if x['panel'] == 'la_n5')
        self.assertLess(lan5['holm_combined'], .05)
        self.assertGreater(lan5['markov'], .05)
        self.assertEqual(self.result['verdicts']['E3-X'],
                         'untestable with current public aggregates')

    def test_ctdna_bounds_include_decimal_rounding_without_ci_claims(self):
        c = self.result['ctdna']
        self.assertLess(c['rounding_envelope'][0], c['forecast'][0])
        self.assertGreater(c['rounding_envelope'][1], c['forecast'][1])
        # Exercise both other locations of the kink without production algebra.
        for m, a, b, d in ((1, 0, 2, 3), (5, 4, 6, 3), (2, 0, 4, 3)):
            expected = convex_thinning_bounds(F(m), F(a), F(b), F(d))
            for actual, target in zip(analysis.thinning_bounds(m, a, b, d), expected):
                self.assertClose(actual, float(target))
        self.assertGreater(c['width_rounding'][0], c['cadence'])
        for prediction in ('C1', 'C2', 'C3'):
            self.assertEqual(self.result['verdicts'][prediction],
                             'untestable with current public aggregates')
        self.assertIn('Detected subgroup',
                      self.rows['henriksen22_matched_lead_median']['notes'])
        self.assertIn('delayed subgroup', self.rows['henriksen24_gain_median']['notes'])
        self.assertEqual(c['tracerx_summary'][1], 0)

    def test_ecdna_classification_and_source_tests_remain_separate(self):
        for item in self.result['clones']:
            clone = item['clone']
            if clone in (4, 6):
                self.assertEqual(item['scope'], 'HSR control')
                self.assertIsNone(item['holm9'])
                self.assertIsNone(item['holm10'])
            elif clone == 7:
                self.assertEqual(item['scope'], 'classification conflict')
                self.assertIsNone(item['holm9'])
                self.assertIsNotNone(item['holm10'])
                self.assertIn('conflict', self.rows[item['id']+'_p']['notes'])
            else:
                self.assertEqual(item['scope'], 'ecDNA')
                self.assertIsNotNone(item['holm9'])
            self.assertLessEqual(0, item['R'])
            self.assertLessEqual(item['R'], 1)
        self.assertEqual(sum(x['holm9'] is not None for x in self.result['clones']), 9)
        self.assertEqual(sum(x['holm10'] is not None for x in self.result['clones']), 10)
        self.assertEqual(len(self.result['lange']), 8)
        self.assertTrue(all('CI' in self.rows['shibata_2026_'+x['arm']+'_p']['notes']
                            for x in self.result['arms']))

    def test_lundberg_rounding_exact_tails_and_dependence_bound(self):
        nominal = F(self.rows['lundberg_2008_alpha']['value'])
        for item in self.result['lundberg']:
            expected_counts = integer_rounding_envelope(item['n'], F(str(item['percent'])))
            self.assertEqual(item['counts'], expected_counts)
            for k in expected_counts:
                exact = binomial_tail_by_recurrence(k, item['n'], nominal)
                self.assertEqual(analysis.binomial_upper(k, item['n'], nominal), float(exact))
                self.assertLessEqual(float(exact), item['markov'])
        # At an unspecified half-percent tie, retain both outer boundaries.
        self.assertEqual(analysis.rounded_counts(200, 10), [19, 20, 21])
        self.assertEqual(analysis.rounded_counts(3, 2), [])
        self.assertEqual(analysis.binomial_upper(0, 5, 0), 1)
        self.assertEqual(analysis.binomial_upper(1, 5, 0), 0)
        self.assertEqual(analysis.binomial_upper(5, 5, 1), 1)

    def test_eligibility_missing_inputs_do_not_turn_into_zeros(self):
        expected = {'I1', 'I2-H', 'C1', 'C2', 'C3', 'E1', 'E2', 'E3',
                    'E4', 'R1', 'R2'}
        audit = self.result['eligibility']
        self.assertEqual(set(audit), expected)
        self.assertEqual(set(self.result['verdicts']), expected | {'I2', 'I3', 'E3-X'})
        for prediction, entry in audit.items():
            with self.subTest(prediction=prediction):
                self.assertTrue(entry['missing'])
                self.assertTrue(entry['reason'])
                self.assertTrue(all(anchor in self.rows for anchor in entry['anchors']))
                self.assertEqual(entry['verdict'], 'untestable with current public aggregates')
                damaged = dict(self.rows)
                del damaged[entry['anchors'][0]]
                with self.assertRaises(ValueError):
                    analysis.eligibility(damaged)
        self.assertIn('one_copy_founders', audit['E4']['missing'])
        self.assertIn('exact_grid_adherence', audit['C1']['missing'])
        self.assertIn('continuous_recurrence_or_actual_imaging_lag_bound',
                      audit['C3']['missing'])
        self.assertIn('simultaneous_uncertainty', audit['R1']['missing'])
        self.assertIn('uncertainty_or_coefficient_covariance', audit['I2-H']['missing'])
        self.assertTrue(self.result['verdicts']['I2'].startswith('consistent'))
        self.assertIn('nominal counterevidence', self.result['verdicts']['I3'])

    def test_verdict_directions_and_primary_family_with_synthetic_aggregates(self):
        # These altered aggregate scalars exercise the decision rule only;
        # they are not additional observations attributed to the publications.
        effects = {x['id']: x for x in self.expected['effects']}
        adjusted = self.expected['formal_holm']
        i2_counter = any(effects[key]['effect'] < 1 and adjusted[key] < .05
                         for key in ('litchfield_clonal', 'ravi_clonal'))
        i3_counter = any(effects[key]['effect'] < 1 and adjusted[key] < .05
                         for key in ('lozac_luad', 'lozac_lusc'))
        i3_counter |= (effects['litchfield_hla']['effect'] > 1
                       and adjusted['litchfield_hla'] < .05)
        i3_counter |= any(x['or_'] > 1 and adjusted[x['id']] < .05
                          for x in self.expected['yang'])
        for key, counter in (('I2', i2_counter), ('I3', i3_counter)):
            expected_prefix = 'inconsistent (falsified)' if counter else 'consistent'
            self.assertTrue(self.result['verdicts'][key].startswith(expected_prefix))

        def altered(changes):
            rows = {key: dict(row) for key, row in self.rows.items()}
            for key, value in changes.items():
                rows[key]['value'] = value
            return analysis.calculate(rows)

        # A coarse q-bound is sensitivity analysis, not a replacement for the
        # frozen primary family. This example makes the distinction observable.
        result = altered({'ravi_clonal_or': '0.2', 'ravi_clonal_low': '0.1',
                          'ravi_clonal_high': '0.4', 'ravi_clonal_q': '0.9'})
        self.assertLess(result['formal_holm']['ravi_clonal'], .05)
        self.assertGreater(result['formal_holm_q_bound']['ravi_clonal'], .05)
        self.assertEqual(result['verdicts']['I2'], 'inconsistent (falsified)')
        result = altered({'lozac_luad_p': '0.000000001'})
        self.assertEqual(result['verdicts']['I3'], 'inconsistent (falsified)')
        # HR > 1 means worse survival, whereas response OR > 1 means better
        # response. A strong adverse HR must not be called contrary evidence.
        result = altered({'lozac_luad_hr': '2', 'lozac_luad_low': '1.2',
                          'lozac_luad_high': '3', 'lozac_luad_p': '0.000000001'})
        self.assertTrue(result['verdicts']['I3'].startswith('consistent'))

    def test_loader_rejects_missing_provenance_and_nonaggregate_schema(self):
        def check_bad(filename, mutate, columns=CSV_COLUMNS):
            records = list(csv.DictReader(io.StringIO(self.csv_texts[filename])))
            mutate(records)
            handle = io.StringIO()
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(records)
            supplied = dict(self.csv_texts, **{filename: handle.getvalue()})

            def fake_open(path, *args, **kwargs):
                return io.StringIO(supplied[path.name])

            with patch.object(Path, 'open', fake_open), self.assertRaises(ValueError):
                analysis.load_data(Path('unused-in-memory-directory'))

        for field, value in (('locator', ''), ('value', 'NaN'),
                             ('source_url', 'http://example.invalid/')):
            check_bad('immune.csv', lambda rows, f=field, v=value: rows[0].update({f: v}))
        check_bad('immune.csv', lambda rows: rows.append(dict(rows[0])))
        check_bad('immune.csv', lambda rows: None, CSV_COLUMNS+('patient_id',))

    def test_numerical_endpoint_and_parameter_guards(self):
        invalid = (
            (analysis.wilson, (-1, 10)), (analysis.wilson, (1, 0)),
            (analysis.wilson, (11, 10)),
            (analysis.hypergeom_weights, (-1, 2, 3, 4)),
            (analysis.hypergeom_weights, (1.5, 2, 3, 4)),
            (analysis.conditional_or_interval, (1, 2, 3, 4, 0)),
            (analysis.conditional_or_interval, (1, 2, 3, 4, 1)),
            (analysis.contrast, (0, 1, 1, 1)),
            (analysis.log_or_p, (1, 2, 3)),
            (analysis.thinning_bounds, (2, 0, 1, 3)),
            (analysis.thinning_bounds, (1, 0, 2, 0)),
            (analysis.rounded_thinning_bounds, (2, 0, 4, 3)),
            (analysis.rounded_counts, (0, 10)),
            (analysis.rounded_counts, (3.5, 10)),
            (analysis.binomial_upper, (6, 5, .1)),
            (analysis.binomial_upper, (1, 5, 1.1)),
            (analysis.effective_R, (-.1, .1)),
            (analysis.effective_R, (.5, 1.1)),
        )
        for function, args in invalid:
            with self.subTest(function=function.__name__, args=args):
                with self.assertRaises(ValueError):
                    function(*args)
        for p in (-.1, 1.1, float('nan')):
            with self.assertRaises(ValueError):
                analysis.holm({'invalid': p})

    def test_generated_readme_and_each_numeric_table_match_csv_recomputation(self):
        text = (RESEARCH / 'README.md').read_text(encoding='utf-8')
        begin, end = '<!-- BEGIN RECOMPUTED RESULTS -->', '<!-- END RECOMPUTED RESULTS -->'
        self.assertEqual(text.count(begin), 1)
        self.assertEqual(text.count(end), 1)
        block = text.split(begin, 1)[1].split(end, 1)[0].strip()
        self.assertEqual(block, analysis.markdown(self.result).strip())
        analysis.check_readme(self.result)
        self.assertIn('Nonclinical', text)
        table_rows = [[cell.strip() for cell in line.strip().strip('|').split('|')]
                      for line in block.splitlines() if line.startswith('|')]

        for filename, count in self.counts.items():
            self.assertIn([f'`{filename}`', str(count)], table_rows)
        for actual, expected in zip(self.result['effects'], self.expected['effects']):
            row = next(r for r in table_rows if r[0] == actual['label'])
            self.assertEqual(row[1:], [actual['kind'].upper()+' '+display(expected['effect']),
                                      display_interval(expected['ci']), display(expected['p']),
                                      actual['p_kind'], display(expected.get('holm')),
                                      display(expected.get('holm_q_bound'))])
        for actual, expected in zip(self.result['yang'], self.expected['yang']):
            label = actual['label']+' / '+actual['locus']
            a, b, c, d = expected['counts']
            first = next(r for r in table_rows if r[0] == label and len(r) == 7)
            self.assertEqual(first[1:], [f'{a}, {b}', f'{c}, {d}', display(expected['or_']),
                                        display_interval(expected['or_ci']), display(expected['p']),
                                        display(expected['holm'])])
            second = next(r for r in table_rows if r[0] == label and len(r) == 5)
            self.assertEqual(second[1:], [
                display(expected['risk_loh'])+' '+display_interval(expected['risk_loh_ci']),
                display(expected['risk_intact'])+' '+display_interval(expected['risk_intact_ci']),
                display(expected['rd'])+' '+display_interval(expected['rd_ci']),
                display(expected['rr'])+' '+display_interval(expected['rr_ci'])])
        c = self.expected['ctdna']
        ctdna_cells = {
            'Reinert original': [c['mean'], c['minimum'], c['maximum'], c['cadence']],
            'Both-parity thinning': c['forecast'], 'Forecast loss': c['loss'],
            'Outward decimal': c['rounding_envelope'],
            'Original lead range': [c['width'], *c['width_rounding']],
            'Reinert detected': c['recurrence_counts'],
            'Henriksen 2022 CT-date': [*c['matched_counts'], *c['matched_median_iqr']],
            'Log of Henriksen': [c['slow_log_of_mean'], c['fast_log_of_mean']],
            'Henriksen 2024 frequency': [*c['extrapolated_counts'], *c['extrapolated_median_iqr']],
            'TRACERx median': c['tracerx_summary'],
        }
        for prefix, values in ctdna_cells.items():
            cell = next(r[1] for r in table_rows if r[0].startswith(prefix))
            tokens = re.findall(r'(?<![A-Za-z])[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?', cell)
            self.assertEqual([F(t) for t in tokens], [F(display(v)) for v in values])
        for actual, expected in zip(self.result['clones'], self.expected['clones']):
            row = next(r for r in table_rows if r[0] == str(expected['clone']) and len(r) == 9)
            self.assertEqual(row[1:], [display(expected['p']), display(expected['rho']),
                                      display(expected['R']), display(expected['ks_p']),
                                      str(expected['n']), actual['scope'],
                                      display(expected['holm9']), display(expected['holm10'])])
        classified = [r['R'] for r in self.expected['clones'] if r['clone'] not in (4, 6, 7)]
        self.assertIn(display_interval((min(classified), max(classified))), block)
        for arm in self.expected['arms']:
            self.assertIn([arm['arm'], display(arm['p']), display(arm['rho']), display(arm['R'])],
                          table_rows)
        for panel in self.expected['lange']:
            self.assertIn([panel['panel'], display(panel['D']), display(panel['p']), str(panel['n'])],
                          table_rows)
        for item in self.expected['lundberg']:
            row = next(r for r in table_rows if r[0] == item['panel'] and len(r) == 9)
            self.assertEqual(row[1:], [display(item['percent'])+'% / '+str(item['n']),
                                      ', '.join(map(str, item['counts'])),
                                      display(item['rate'][0])+' '+display_interval(item['ci']),
                                      display(item['excess'][0])+' / '+display(item['fold'][0]),
                                      display(item['p']), display(item['holm3']),
                                      display(item['holm_combined']), display(item['markov'])])
        ravi = next(x for x in self.expected['effects'] if x['id'] == 'ravi_clonal')
        for phrase, number in (('Ravi source BH q = ', ravi['source_q']),
                               ('rounding-conditional upper bound = ', self.expected['ravi_q_rounding_upper']),
                               ('reconstructed log-OR SE = ', ravi['se']), ('z = ', ravi['z'])):
            self.assertIn(phrase+display(number), block)
        r = self.expected['resistance']
        for phrase, number in (('approximately p = ', r['persister']), ('VMR = ', r['vmr']),
                               ('1-p = ', r['fixed_size_binomial_vmr']),
                               ('excess = ', r['excess']), ('ratio = ', r['ratio'])):
            self.assertIn(phrase+display(number), block)
        self.assertIn(f"with {r['n']} retained lineages", block)
        # Synchronization must fail if an otherwise valid numeric result is stale.
        changed = dict(self.result, row_count=self.result['row_count']+1)
        changed['file_counts'] = dict(self.result['file_counts'], **{'immune.csv': -1})
        with self.assertRaises(ValueError):
            analysis.check_readme(changed)


if __name__ == '__main__':
    unittest.main()
