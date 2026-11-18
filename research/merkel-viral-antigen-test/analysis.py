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


@lru_cache(maxsize=1)
def frozen_curves():
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
    return compatibility(curve, m, max(0.0, 2*lower-curve['A']),
                         min(curve['A'], 2*upper))


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
    actual = normalized_bytes(HERE / 'data' / 'part_a_calibration.csv').decode('utf-8')
    if actual != calibration_text():
        raise ValueError('Part A calibration changed; do not score with an altered freeze')
    return True


def manifest():
    paths = [PR67 / 'analysis.py', PR67 / 'data' / 'tmb_orr.csv',
             HERE / 'PROTOCOL.md', HERE / 'analysis.py',
             ROOT / 'tests' / 'test_research_merkel_viral_antigen_test.py']
    paths += [HERE / 'data' / name for name in
              ('part_a.csv', 'part_b.csv', 'part_a_calibration.csv')]
    return {path.relative_to(ROOT).as_posix(): digest(path) for path in sorted(paths)}


def check_manifest():
    expected = {row['path']: row['sha256'] for row in read_csv(HERE / 'data' / 'fingerprints.csv')}
    if manifest() != expected:
        raise ValueError('Fingerprint mismatch (canonical LF text / POSIX path keys)')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--freeze-part-a', action='store_true',
                        help='Write calibration once, reading only Part A; commit before any Part B scoring')
    args = parser.parse_args()
    if args.freeze_part_a:
        target = HERE / 'data' / 'part_a_calibration.csv'
        if target.exists():
            raise SystemExit('Refusing to overwrite frozen Part A calibration')
        target.write_text(calibration_text(), encoding='utf-8', newline='\n')
        print('Part A calibration written; commit it before scoring Part B.')


if __name__ == '__main__':
    main()
