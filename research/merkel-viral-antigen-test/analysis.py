#!/usr/bin/env python3
"""Offline NONCLINICAL aggregate sensitivity analysis; standard library only.

The immutable protocol precedes computation. The calibration command reads
Part A only; the report command checks that committed calibration before using
Part B. Median plug-in predictions are ecological, not individual estimates.
"""

import argparse
import csv
from decimal import Decimal
from functools import lru_cache
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PR67 = ROOT / 'research' / 'immune-escape-pancancer-test'
PROTOCOL_COMMIT = '09ecbb18e6b159d8cf74d30eac82059dc54e472e'
BASE_COMMIT = '16c7596637a60836056593c1fb22b1f5c9d07501'
CALIBRATION_SHA256 = 'ef597c8dd3ff770848fb06092853e620136a014c1fdd76c9e3873618cbf3af6e'
SCALES = (('0.5', 0.5), ('1', 1.0), ('4/3', 4/3), ('2', 2.0), ('5', 5.0))
BANDS = {'counting_assumption': (1.0, 2.0),
         'wide_assumption': (0.5, 5.0)}
CURVE_LABELS = {
    'reciprocal_holdout': 'Reciprocal, Merkel held out (primary)',
    'reciprocal_all': 'Reciprocal, full panel',
    'log_holdout': 'Log TMB, Merkel held out',
    'log_all': 'Log TMB, full panel',
}
BEGIN = '<!-- BEGIN RECOMPUTED RESULTS -->'
END = '<!-- END RECOMPUTED RESULTS -->'


def normalized_bytes(path):
    """Canonical text bytes, independent of git's checkout newline convention."""
    return Path(path).read_bytes().decode('utf-8').replace('\r\n', '\n').replace(
        '\r', '\n').encode('utf-8')


def digest(path):
    return hashlib.sha256(normalized_bytes(path)).hexdigest()


def read_csv(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def load_inputs():
    rows = read_csv(HERE / 'data' / 'part_a.csv') + read_csv(HERE / 'data' / 'part_b.csv')
    result = {}
    required = ('source_title', 'doi', 'source_url', 'data_url', 'locator',
                'access_date', 'reuse_basis', 'license_url', 'inspection', 'notes')
    for row in rows:
        if row['id'] in result:
            raise ValueError('Duplicate aggregate identifier')
        if any(not row[field].strip() for field in required):
            raise ValueError('Missing row provenance: ' + row['id'])
        for key in ('value', 'ci_low', 'ci_high', 'n', 'responders'):
            if row[key] and not math.isfinite(float(row[key])):
                raise ValueError('Nonfinite source aggregate')
        if row['metric'] == 'orr_percent' and row['value']:
            lo, mid, hi = (float(row[k]) for k in ('ci_low', 'value', 'ci_high'))
            if not 0 <= lo <= mid <= hi <= 100 or not row['n']:
                raise ValueError('Invalid published ORR row')
        if row['metric'] == 'tmb_median' and (float(row['value']) <= 0 or not row['n']):
            raise ValueError('Invalid TMB row')
        result[row['id']] = row
    return result


def replay_curves():
    """Recreate exactly the pinned PR67 fits, using only its TMB/ORR CSV.

    Reading the raw CSV avoids unrelated clonal/HLA data dependencies. Their
    fields do not enter these original fits. This is a deterministic replay,
    never a JAVELIN refit or a new model-selection step.
    """
    spec = importlib.util.spec_from_file_location('merkel_original_pr67', PR67 / 'analysis.py')
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    rows = []
    for raw in read_csv(PR67 / 'data' / 'tmb_orr.csv'):
        if raw['tmb']:
            rows.append(dict(raw, tmb=float(raw['tmb']), orr=float(raw['orr'])))
    result = {}
    for suffix, subset in (('holdout', [r for r in rows if r['parent_type'] != 'merkel']),
                           ('all', rows)):
        for kind, model in (('reciprocal', 'reciprocal_tmb'), ('log', 'log_tmb')):
            fit = original.fit(subset, model)
            curve = {'kind': kind, 'n': len(subset),
                     'parents': len({r['parent_type'] for r in subset})}
            if kind == 'reciprocal':
                curve.update(A=fit['A0'], B=fit['B'])
            else:
                curve.update(a=fit['intercept'], b=fit['slope'])
            result[kind + '_' + suffix] = curve
    return result


@lru_cache(maxsize=1)
def frozen_curves():
    """Validate the original fits, then use their recorded full-precision values.

    Platform libm rounding can move the last digits of a flat least-squares
    optimum. Fixed snapshots keep the same frozen curve on every platform;
    a bounded numerical comparison still independently recomputes each fit.
    Neither replay nor snapshots use JAVELIN data.
    """
    replay = replay_curves()
    result = {}
    for raw in read_csv(HERE / 'data' / 'frozen_curves.csv'):
        name = raw['name']
        if name in result or name not in replay:
            raise ValueError('Unexpected/duplicate frozen-curve snapshot')
        if (raw['pr67_code_sha256'] != digest(PR67 / 'analysis.py') or
                raw['pr67_csv_sha256'] != digest(PR67 / 'data' / 'tmb_orr.csv')):
            raise ValueError('Frozen original source fingerprint mismatch')
        curve = {'kind': raw['kind'], 'n': int(raw['n']), 'parents': int(raw['parents'])}
        keys = ('A', 'B') if curve['kind'] == 'reciprocal' else ('a', 'b')
        curve.update({key: float(raw[key]) for key in keys})
        if any(curve[k] != replay[name][k] for k in ('kind', 'n', 'parents')):
            raise ValueError('Frozen curve identity mismatch')
        if any(not math.isclose(curve[k], replay[name][k], rel_tol=0, abs_tol=2e-6) for k in keys):
            raise ValueError('Original curve replay exceeds numerical tolerance')
        result[name] = curve
    if set(result) != set(replay):
        raise ValueError('Missing frozen curve')
    return result


def predict(curve, x):
    if x < 0 or not math.isfinite(x):
        raise ValueError('Load must be finite and nonnegative')
    if curve['kind'] == 'reciprocal':
        if curve['B'] == 0:
            return curve['A']
        return curve['A'] * math.exp(-curve['B']/x) if x else 0.0
    if curve['b'] == 0:
        return min(1.0, max(0.0, curve['a']))
    if not x:
        return 0.0 if curve['b'] > 0 else 1.0
    return min(1.0, max(0.0, curve['a'] + curve['b']*math.log(x)))


def intersect(*regions):
    """Closed finite endpoints; zero/infinity denote excluded domain limits."""
    if any(region is None for region in regions):
        return None
    lo = max((r[0] for r in regions), default=0.0)
    hi = min((r[1] for r in regions), default=math.inf)
    if lo > hi or hi <= 0 or math.isinf(lo):
        return None
    return lo, hi


def compatibility(curve, m, lower, upper):
    """All finite positive scales whose median plug-in is in [lower, upper]."""
    if not math.isfinite(m) or m <= 0 or lower > upper:
        raise ValueError('Invalid median/response interval')
    if upper < 0 or lower > 1:
        return None
    if curve['kind'] == 'reciprocal':
        a, b = curve['A'], curve['B']
        if b == 0 or a == 0:
            return (0.0, math.inf) if lower <= a <= upper else None
        if a < 0 or b < 0:
            raise ValueError('Curve is not monotone in the declared family')
        if upper <= 0 or lower >= a:
            return None  # neither 0 nor A is attained at positive finite x
        lo = 0.0 if lower <= 0 else b / math.log(a/lower) / m
        hi = math.inf if upper >= a else b / math.log(a/upper) / m
    else:
        a, b = curve['a'], curve['b']
        if b == 0:
            p = min(1.0, max(0.0, a))
            return (0.0, math.inf) if lower <= p <= upper else None
        if b < 0:
            raise ValueError('Decreasing log curve not in the frozen family')
        lo = 0.0 if lower <= 0 else math.exp((lower-a)/b)/m
        hi = math.inf if upper >= 1 else math.exp((upper-a)/b)/m
    return intersect((lo, hi))


def half_unit(value):
    """Half the last printed unit, including explicitly retained trailing zeros."""
    return float(Decimal(1).scaleb(Decimal(value).as_tuple().exponent)/2)


def rounded_compatibility(curve, tmb_row, orr_row):
    m = float(tmb_row['value'])
    dm = half_unit(tmb_row['value'])
    lower = max(0.0, (float(orr_row['ci_low']) - half_unit(orr_row['ci_low']))/100)
    upper = min(1.0, (float(orr_row['ci_high']) + half_unit(orr_row['ci_high']))/100)
    if m <= dm:
        raise ValueError('Rounded TMB crosses zero')
    # Union over possible medians: the low-c boundary uses the largest m;
    # the high-c boundary uses the smallest m. Monotonicity makes it exact.
    largest = compatibility(curve, m+dm, lower, upper)
    smallest = compatibility(curve, m-dm, lower, upper)
    return None if largest is None else (largest[0], smallest[1])


def median_compatibility(curve, m, lower, upper):
    if curve['kind'] != 'reciprocal' or curve['B'] <= 0:
        raise ValueError('Median envelope requires a positive-B reciprocal curve')
    transformed_lower = max(0.0, 2*lower-curve['A'])
    transformed_upper = min(curve['A'], 2*upper)
    if transformed_lower > transformed_upper:
        return None
    return compatibility(curve, m, transformed_lower, transformed_upper)


def viral_calibration(curve, m, target, scale):
    """Only one additional parameter; the scale is an externally fixed scenario."""
    if not (m > 0 and scale > 0 and 0 <= target <= 1):
        raise ValueError('Invalid calibration target')
    baseline = predict(curve, scale*m)
    if target <= baseline:
        return {'v': 0.0, 'status': 'zero_boundary', 'fitted_orr': baseline}
    if target >= curve['A']:
        return {'v': None, 'status': 'unattained_ceiling', 'fitted_orr': None}
    v = curve['B']/math.log(curve['A']/target) - scale*m
    return {'v': v, 'status': 'finite', 'fitted_orr': predict(curve, scale*m+v)}


def calibrate_part_a(rows=None, curves=None):
    """Intentionally has no Part B read, even when called with a combined dict."""
    if rows is None:
        rows = {r['id']: r for r in read_csv(HERE / 'data' / 'part_a.csv')}
    curves = frozen_curves() if curves is None else curves
    curve = curves['reciprocal_holdout']
    tmb = rows['a_tmb_positive']
    m = float(tmb['value'])
    result = []
    for analysis_set, orr_id in (('broad', 'a_orr_positive'),
                                 ('tmb_subset', 'a_subset_orr_positive')):
        orr = rows[orr_id]
        for scale_label, scale in SCALES:
            for target, field in (('point', 'value'), ('ci_low', 'ci_low'),
                                  ('ci_high', 'ci_high')):
                y = float(orr[field])/100
                row = {'scale_label': scale_label, 'scale': scale,
                       'calibration_set': analysis_set, 'target': target,
                       'A': curve['A'], 'B': curve['B'],
                       'source_orr_id': orr_id, 'source_tmb_id': tmb['id'],
                       'target_orr': y, 'median': m}
                row.update(viral_calibration(curve, m, y, scale))
                result.append(row)
    return result


CALIBRATION_FIELDS = (
    'scale_label', 'scale', 'calibration_set', 'target', 'A', 'B',
    'source_orr_id', 'source_tmb_id', 'target_orr', 'median', 'v', 'status',
    'fitted_orr', 'n_orr', 'n_tmb', 'source_title', 'source_url', 'locator',
    'access_date', 'reuse_basis', 'license_url', 'part_a_sha256',
    'pr67_code_sha256', 'pr67_csv_sha256', 'protocol_sha256')


def csv_text(rows, fields):
    handle = io.StringIO(newline='')
    writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
    writer.writeheader()
    for row in rows:
        writer.writerow({key: '' if row[key] is None else row[key] for key in fields})
    return handle.getvalue()


def calibration_text():
    source = {r['id']: r for r in read_csv(HERE / 'data' / 'part_a.csv')}
    rows = []
    for fit in calibrate_part_a(source):
        orr, tmb = source[fit['source_orr_id']], source[fit['source_tmb_id']]
        row = dict(fit, n_orr=orr['n'], n_tmb=tmb['n'],
                   source_title='Derived Part A calibration; ' + orr['source_title'],
                   source_url=orr['source_url'],
                   locator=orr['locator'] + '; ' + tmb['locator'] + '; PROTOCOL.md viral formula',
                   access_date=orr['access_date'],
                   reuse_basis='Derived nonclinical calculation from attributed CC BY 4.0 aggregate facts; not a published antigen measurement',
                   license_url=orr['license_url'],
                   part_a_sha256=digest(HERE / 'data' / 'part_a.csv'),
                   pr67_code_sha256=digest(PR67 / 'analysis.py'),
                   pr67_csv_sha256=digest(PR67 / 'data' / 'tmb_orr.csv'),
                   protocol_sha256=digest(HERE / 'PROTOCOL.md'))
        rows.append(row)
    return csv_text(rows, CALIBRATION_FIELDS)


def check_calibration():
    path = HERE / 'data' / 'part_a_calibration.csv'
    if digest(path) != CALIBRATION_SHA256:
        raise ValueError('Immutable Part A calibration changed')
    actual = read_csv(path)
    expected = list(csv.DictReader(io.StringIO(calibration_text())))
    if len(actual) != len(expected):
        raise ValueError('Part A calibration row count mismatch')
    numerical = {'scale', 'A', 'B', 'target_orr', 'median', 'v', 'fitted_orr'}
    for observed, recomputed in zip(actual, expected):
        for key in CALIBRATION_FIELDS:
            first, second = observed[key], recomputed[key]
            if key in numerical and first and second:
                same = math.isclose(float(first), float(second), rel_tol=5e-12, abs_tol=5e-12)
            else:
                same = first == second
            if not same:
                raise ValueError('Part A calibration replay mismatch: ' + key)
    return True


def recorded_calibration():
    """Use the committed pre-Part-B values after their recomputation check."""
    keys = CALIBRATION_FIELDS[:13]
    numerical = {'scale', 'A', 'B', 'target_orr', 'median', 'v', 'fitted_orr'}
    return [{key: (float(row[key]) if row[key] else None) if key in numerical else row[key]
             for key in keys}
            for row in read_csv(HERE / 'data' / 'part_a_calibration.csv')]


def manifest():
    paths = [PR67 / 'analysis.py', PR67 / 'data' / 'tmb_orr.csv',
             HERE / 'PROTOCOL.md', HERE / 'analysis.py',
             ROOT / 'tests' / 'test_research_merkel_viral_antigen_test.py']
    paths += [HERE / 'data' / name for name in
              ('part_a.csv', 'part_b.csv', 'part_a_calibration.csv', 'frozen_curves.csv')]
    keyed = {path.relative_to(ROOT).as_posix(): path for path in paths}
    return {key: digest(keyed[key]) for key in sorted(keyed)}


def check_manifest():
    expected = {row['path']: row['sha256'] for row in read_csv(HERE / 'data' / 'fingerprints.csv')}
    if manifest() != expected:
        raise ValueError('Fingerprint mismatch (canonical LF text / POSIX path keys)')
    return True


def response(row):
    return tuple(float(row[key])/100 for key in ('value', 'ci_low', 'ci_high'))


def region_panel(curve, inputs, mode='printed'):
    result = {}
    for cohort in ('A', 'B'):
        prefix = cohort.lower()
        groups = {}
        for status in ('positive', 'negative'):
            tmb = inputs[f'{prefix}_tmb_{status}']
            source = f'{prefix}_subset_orr_{status}' if mode == 'subset_a' and cohort == 'A' else f'{prefix}_orr_{status}'
            orr = inputs[source]
            _, lower, upper = response(orr)
            if mode == 'rounded':
                region = rounded_compatibility(curve, tmb, orr)
            elif mode == 'median_envelope':
                region = median_compatibility(curve, float(tmb['value']), lower, upper)
            else:
                region = compatibility(curve, float(tmb['value']), lower, upper)
            groups[status] = region
        groups['joint'] = intersect(groups['positive'], groups['negative'])
        result[cohort] = groups
    result['joint'] = intersect(result['A']['joint'], result['B']['joint'])
    result['bands'] = {key: intersect(result['joint'], band) for key, band in BANDS.items()}
    return result


def score(predictions, observations):
    """Equal-status descriptive point loss, not a patient-level likelihood."""
    errors = [y-p for p, (y, _, _) in zip(predictions, observations)]
    distances = [max(lo-p, 0.0, p-hi) for p, (_, lo, hi) in zip(predictions, observations)]
    n = len(errors)
    mse = math.fsum(e*e for e in errors)/n
    return {'mse': mse, 'rmse': math.sqrt(mse),
            'mae': math.fsum(abs(e) for e in errors)/n,
            'misses': sum(d > 0 for d in distances),
            'interval_mse': math.fsum(d*d for d in distances)/n,
            'residuals': errors, 'distances': distances}


def build_report(inputs=None, curves=None):
    check_calibration()
    inputs = load_inputs() if inputs is None else inputs
    curves = frozen_curves() if curves is None else curves
    regions = {key: {mode: region_panel(curve, inputs, mode)
                     for mode in ('printed', 'rounded', 'subset_a')}
               for key, curve in curves.items()}
    primary = curves['reciprocal_holdout']
    median = region_panel(primary, inputs, 'median_envelope')
    observations_b = [response(inputs['b_orr_'+status]) for status in ('positive', 'negative')]
    baseline = []
    for name, curve in curves.items():
        for label, scale in SCALES:
            predictions = {cohort: [predict(curve, scale*float(inputs[f'{cohort.lower()}_tmb_{status}']['value']))
                                    for status in ('positive', 'negative')]
                           for cohort in ('A', 'B')}
            baseline.append({'curve': name, 'scale_label': label, 'scale': scale,
                             'predictions': predictions, 'score': score(predictions['B'], observations_b)})
    calibrated = recorded_calibration()
    extensions = []
    for fit in calibrated:
        if fit['target'] != 'point':
            continue
        base = next(row for row in baseline if row['curve'] == 'reciprocal_holdout'
                    and row['scale_label'] == fit['scale_label'])
        if fit['v'] is None:
            extensions.append(dict(fit, predictions=None, score=None, repair=False))
            continue
        predictions = [predict(primary, fit['scale']*float(inputs['b_tmb_'+status]['value'])
                               + (fit['v'] if status == 'positive' else 0.0))
                       for status in ('positive', 'negative')]
        loss = score(predictions, observations_b)
        endpoint_predictions = {}
        for endpoint in ('ci_low', 'ci_high'):
            edge = next(r for r in calibrated if r['scale_label'] == fit['scale_label']
                        and r['calibration_set'] == fit['calibration_set'] and r['target'] == endpoint)
            endpoint_predictions[endpoint] = (None if edge['v'] is None else
                predict(primary, fit['scale']*float(inputs['b_tmb_positive']['value']) + edge['v']))
        extensions.append(dict(fit, predictions=predictions, score=loss,
                               delta_mse=loss['mse']-base['score']['mse'],
                               repair=loss['mse'] < base['score']['mse'] and loss['misses'] == 0,
                               endpoint_predictions=endpoint_predictions))
    merkel = next(r for r in read_csv(PR67 / 'data' / 'tmb_orr.csv') if r['parent_type'] == 'merkel')
    historic = {'tmb': float(merkel['tmb']), 'orr': float(merkel['orr']),
                'prediction': predict(primary, float(merkel['tmb']))}
    historic['residual'] = historic['orr'] - historic['prediction']
    return {'inputs': inputs, 'curves': curves, 'regions': regions,
            'median_envelope': median, 'calibration': calibrated,
            'baseline': baseline, 'extensions': extensions, 'historic': historic}


def fmt_region(region):
    if region is None:
        return 'empty'
    lo, hi = region
    return ('(0, ' if lo == 0 else f'[{lo:.3f}, ') + ('infinity)' if math.isinf(hi) else f'{hi:.3f}]')


def fmt_pct(value):
    if value is None:
        return 'unattained'
    return f'{100*value:.3f}' if value >= 0.000005 else f'{100*value:.3e}'


def render_report(report):
    """All empirical numeric tables in the README are derived here."""
    inputs, curves = report['inputs'], report['curves']
    lines = ['### Published marginal inputs', '',
             'ORR is percent [published 95% CI]; TMB is median NSSV/Mb. '
             'Adjacent cells are an ecological pairing, not individual records. '
             'Unknown responder counts remain unreported.', '',
             '| Cohort / viral stratum | ORR % [95% CI] | ORR denominator | Published responders | TMB median | TMB denominator |',
             '| --- | --- | ---: | ---: | ---: | ---: |']
    for co in ('a', 'b'):
        for status in ('all', 'positive', 'negative'):
            r, t = inputs[f'{co}_orr_{status}'], inputs[f'{co}_tmb_{status}']
            lines.append(f'| {co.upper()} / {status} | {r["value"]} [{r["ci_low"]}, {r["ci_high"]}] | '
                         f'{r["n"]} | {r["responders"] or "unreported"} | {t["value"]} | {t["n"]} |')
    r = inputs['b_orr_unknown']
    lines.append(f'| B / unknown | {r["value"]} [{r["ci_low"]}, {r["ci_high"]}] | {r["n"]} | {r["responders"]} | missing | missing |')
    lines += ['', 'Additional published margins remain separate:', '',
              '| Part A context | ORR % [95% CI] | Denominator | Published responders |',
              '| --- | --- | ---: | ---: |']
    for key in ('a_orr_prior_one', 'a_orr_prior_multiple', 'a_subset_orr_all',
                'a_subset_orr_positive', 'a_subset_orr_negative'):
        r = inputs[key]
        lines.append(f'| {r["population"]} | {r["value"]} [{r["ci_low"]}, {r["ci_high"]}] | {r["n"]} | {r["responders"] or "unreported"} |')
    lines += ['',
              f'Part A explicitly reports {inputs["a_count_viral_unevaluable"]["value"]} participants '
              f'with unevaluable viral status and {inputs["a_count_bor_unevaluable"]["value"]} '
              'with unevaluable best response; these categories must not be equated.', '',
              '| Cohort / endpoint | Source cutoff |', '| --- | --- |',
              f'| A broad ORR | {inputs["a_orr_all"]["data_cutoff"]} |',
              f'| A TMB and TMB-subset ORR | {inputs["a_tmb_all"]["data_cutoff"]} |',
              f'| B ORR and TMB | {inputs["b_orr_all"]["data_cutoff"]} |', '',
              f'Transcribed source rows: Part A {sum(r["cohort"] == "A" for r in inputs.values())}; '
              f'Part B {sum(r["cohort"] == "B" for r in inputs.values())}. '
              'The Part B count includes the explicitly missing viral-specific TMB-subset ORR rows.', '',
              '### Exact frozen-curve replay', '',
              '| Curve | Training rows / parents | A or a | B or b |',
              '| --- | --- | ---: | ---: |']
    for name in CURVE_LABELS:
        curve = curves[name]
        first, second = ('A', 'B') if curve['kind'] == 'reciprocal' else ('a', 'b')
        lines.append(f'| {CURVE_LABELS[name]} | {curve["n"]} / {curve["parents"]} | {curve[first]:.17g} | {curve[second]:.17g} |')
    h = report['historic']
    lines += ['', f'Replayed historical Merkel input: median TMB {h["tmb"]:.2f}, '
              f'ORR {100*h["orr"]:.3f}%; held-out reciprocal prediction '
              f'{100*h["prediction"]:.3f}%, residual {100*h["residual"]:+.3f} pp. '
              'This is the original mixed-source motivating residual, not a JAVELIN viral contrast.', '',
              '### Common-scale compatibility', '',
              'All intervals below are sets of finite positive conversion factors c. '
              'They use printed median substitution and published ORR intervals. '
              'They are not calibrated confidence regions. Infinity denotes an open limit.', '',
              '| Curve | Cohort | MCPyV+ scales | MCPyV- scales | Both statuses |',
              '| --- | --- | --- | --- | --- |']
    for name in CURVE_LABELS:
        panel = report['regions'][name]['printed']
        for co in ('A', 'B'):
            p = panel[co]
            lines.append(f'| {CURVE_LABELS[name]} | {co} | {fmt_region(p["positive"])} | '
                         f'{fmt_region(p["negative"])} | {fmt_region(p["joint"])} |')
    lines += ['', '| Curve | Both cohorts / both statuses | Within assumed [1, 2] | Within assumed [0.5, 5] |',
              '| --- | --- | --- | --- |']
    for name in CURVE_LABELS:
        p = report['regions'][name]['printed']
        lines.append(f'| {CURVE_LABELS[name]} | {fmt_region(p["joint"])} | '
                     f'{fmt_region(p["bands"]["counting_assumption"])} | {fmt_region(p["bands"]["wide_assumption"])} |')
    lines += ['', 'Rounding expands medians and interval endpoints by half the last printed unit; '
              'it does not account for sampling error in the medians:', '',
              '| Curve | Rounded A intersection | Rounded B intersection | Rounded common intersection |',
              '| --- | --- | --- | --- |']
    for name in CURVE_LABELS:
        p = report['regions'][name]['rounded']
        lines.append(f'| {CURVE_LABELS[name]} | {fmt_region(p["A"]["joint"])} | '
                     f'{fmt_region(p["B"]["joint"])} | {fmt_region(p["joint"])} |')
    lines += ['', 'Substituting the published Part A TMB-subset ORRs (secondary) gives:', '',
              '| Curve | A subset MCPyV+ | A subset MCPyV- | A subset intersection | With broad B intersection |',
              '| --- | --- | --- | --- | --- |']
    for name in CURVE_LABELS:
        p = report['regions'][name]['subset_a']
        lines.append(f'| {CURVE_LABELS[name]} | {fmt_region(p["A"]["positive"])} | '
                     f'{fmt_region(p["A"]["negative"])} | {fmt_region(p["A"]["joint"])} | {fmt_region(p["joint"])} |')
    lines += ['', '**The matched Part A-to-Part B test is untestable:** '
              'viral-specific TMB-subset ORRs in Part B are unpublished in the inspected sources. '
              'The last column above still uses broad B and therefore retains the denominator mismatch.', '',
              '### Relaxing median substitution', '',
              'For the primary reciprocal curve, replace f(c*m) by the conservative mean-response '
              'closure [f(c*m)/2, (A+f(c*m))/2]. This assumes the TMB median represents the ORR '
              'population; it does not identify its distribution or remove selection bias.', '',
              '| Cohort | MCPyV+ scales | MCPyV- scales | Both statuses |',
              '| --- | --- | --- | --- |']
    p = report['median_envelope']
    for co in ('A', 'B'):
        lines.append(f'| {co} | {fmt_region(p[co]["positive"])} | {fmt_region(p[co]["negative"])} | {fmt_region(p[co]["joint"])} |')
    lines += ['', f'Common envelope-compatible scales: **{fmt_region(p["joint"])}**. '
              f'Assumed-band intersections: [1, 2] -> {fmt_region(p["bands"]["counting_assumption"])}; '
              f'[0.5, 5] -> {fmt_region(p["bands"]["wide_assumption"])}.', '',
              '### Predictions of every frozen curve at every fixed scale', '',
              'Percent predictions; Part B RMSE is status-balanced, in percentage points. '
              'Misses count predictions outside the published marginal intervals. '
              'Very small predictions use scientific notation to avoid implying an exact zero.', '',
              '| Curve | c | A+ % | A- % | B+ % | B- % | B RMSE pp | B misses |',
              '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for name in CURVE_LABELS:
        for row in (r for r in report['baseline'] if r['curve'] == name):
            a, b = row['predictions']['A'], row['predictions']['B']
            lines.append(f'| {CURVE_LABELS[name]} | {row["scale_label"]} | {fmt_pct(a[0])} | {fmt_pct(a[1])} | '
                         f'{fmt_pct(b[0])} | {fmt_pct(b[1])} | {100*row["score"]["rmse"]:.3f} | {row["score"]["misses"]} |')
    lines += ['', '### Part-A-fixed viral extension scored on Part B', '',
              'Primary calibration uses broad Part A MCPyV+ ORR. Only v is fitted; '
              'the scale is fixed by the scenario. Virus-negative predictions remain unchanged. '
              'Residuals are observed minus predicted, in pp.', '',
              '| c | Frozen v | B+ predicted % | B- predicted % | B+ residual pp | B- residual pp | B+ interval distance pp | B- interval distance pp |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for row in (r for r in report['extensions'] if r['calibration_set'] == 'broad'):
        if row['score'] is None:
            lines.append(f'| {row["scale_label"]} | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable | unavailable |')
            continue
        s = row['score']
        lines.append(f'| {row["scale_label"]} | {row["v"]:.6f} | {fmt_pct(row["predictions"][0])} | '
                     f'{fmt_pct(row["predictions"][1])} | {100*s["residuals"][0]:+.3f} | {100*s["residuals"][1]:+.3f} | '
                     f'{100*s["distances"][0]:.3f} | {100*s["distances"][1]:.3f} |')
    lines += ['', 'MSE and interval MSE are in squared percentage points; RMSE/MAE are pp. '
              'Repair requires smaller MSE and no interval misses. These are descriptive aggregate '
              'scores, not patient-level prediction scores or significance tests.', '',
              '| c | Baseline MSE | Extension MSE | Change MSE | Baseline RMSE | Extension RMSE | Extension MAE | Extension interval MSE | Misses | Repair rule |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |']
    for row in (r for r in report['extensions'] if r['calibration_set'] == 'broad'):
        if row['score'] is None:
            continue
        base = next(r['score'] for r in report['baseline'] if r['curve'] == 'reciprocal_holdout' and r['scale_label'] == row['scale_label'])
        s = row['score']
        lines.append(f'| {row["scale_label"]} | {10000*base["mse"]:.3f} | {10000*s["mse"]:.3f} | '
                     f'{10000*row["delta_mse"]:+.3f} | {100*base["rmse"]:.3f} | {100*s["rmse"]:.3f} | '
                     f'{100*s["mae"]:.3f} | {10000*s["interval_mse"]:.3f} | {s["misses"]} | {"passes" if row["repair"] else "fails"} |')
    lines += ['', 'Required calibration sensitivities (all scored against broad B, not a missing '
              'matched B outcome). The CI-endpoint columns are conditional Part A calibration '
              'scenarios for predicted B+; they are not predictive confidence intervals.', '',
              '| Calibration | c | v at point ORR | B+ % | B RMSE pp | B misses | B+ % from A lower CI | B+ % from A upper CI |',
              '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for row in report['extensions']:
        if row['score'] is None:
            continue
        edge = row['endpoint_predictions']
        lines.append(f'| {row["calibration_set"]} | {row["scale_label"]} | {row["v"]:.6f} | '
                     f'{fmt_pct(row["predictions"][0])} | {100*row["score"]["rmse"]:.3f} | {row["score"]["misses"]} | '
                     f'{fmt_pct(edge["ci_low"])} | {fmt_pct(edge["ci_high"])} |')
    lines += ['', '### Canonical input and implementation fingerprints', '',
              'SHA-256 uses UTF-8 text after CRLF/CR normalization to LF. Path keys are '
              'repository-relative `Path.as_posix()` values. README and the manifest itself '
              'are excluded to avoid a circular fingerprint; the test regenerates this entire block.', '',
              '| Path | SHA-256 |', '| --- | --- |']
    for path, sha in manifest().items():
        lines.append(f'| `{path}` | `{sha}` |')
    return '\n'.join(lines) + '\n'


def check_readme(report=None):
    text = normalized_bytes(HERE / 'README.md').decode('utf-8')
    actual = text.split(BEGIN + '\n', 1)[1].split(END, 1)[0]
    expected = '\n' + render_report(build_report() if report is None else report) + '\n'
    if actual != expected:
        raise ValueError('README numeric block does not match the CSV calculations')
    return True


def json_safe(value):
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and math.isinf(value):
        return 'infinity'
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--freeze-part-a', action='store_true',
                        help='Write calibration once, reading only Part A; commit before any Part B scoring')
    parser.add_argument('--write-report', action='store_true',
                        help='Regenerate fingerprints and the README numeric block; never alter calibration')
    parser.add_argument('--check', action='store_true', help='Recompute and check all stored artifacts without writing')
    parser.add_argument('--json', action='store_true', help='Print the recomputed report as strict JSON')
    args = parser.parse_args()
    if args.freeze_part_a:
        target = HERE / 'data' / 'part_a_calibration.csv'
        if target.exists():
            raise SystemExit('Refusing to overwrite frozen Part A calibration')
        target.write_text(calibration_text(), encoding='utf-8', newline='\n')
        print('Part A calibration written; commit it before scoring Part B.')
        return
    report = build_report()
    if args.write_report:
        fingerprints = [{'path': path, 'sha256': sha} for path, sha in manifest().items()]
        (HERE / 'data' / 'fingerprints.csv').write_text(csv_text(fingerprints, ('path', 'sha256')),
                                                      encoding='utf-8', newline='\n')
        path = HERE / 'README.md'
        text = normalized_bytes(path).decode('utf-8')
        before = text.split(BEGIN, 1)[0]
        after = text.split(END, 1)[1]
        path.write_text(before + BEGIN + '\n\n' + render_report(report) + '\n' + END + after,
                        encoding='utf-8', newline='\n')
        print('Recomputed README and fingerprints; Part A calibration unchanged.')
    if args.check or not (args.write_report or args.json):
        check_manifest()
        check_readme(report)
        print('All canonical fingerprints, Part A calibration and README calculations match.')
    if args.json:
        print(json.dumps(json_safe(report), indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
