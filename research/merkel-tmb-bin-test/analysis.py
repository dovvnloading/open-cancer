#!/usr/bin/env python3
"""NONCLINICAL offline sharp-bin compatibility audit; standard library only.

No JAVELIN parameter is fitted. Original curve snapshots stay fixed. Intervals
retain attainable endpoints; distance infima alone never decide compatibility.
"""

import argparse
import csv
from dataclasses import dataclass
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = HERE / 'data' / 'tmb_bins.csv'
CURVES_PATH = ROOT / 'research' / 'merkel-viral-antigen-test' / 'data' / 'frozen_curves.csv'
BASE_COMMIT = '37c37e2b4cfc711ce938201123ef3085b8868f2f'
LOCAL_PROTOCOL_COMMIT = '6d6d4088c11099c73fd23b79c24ffe0c5c781b79'
PUBLISHED_PROTOCOL_COMMIT = '892df395d5b3fa0c293e76c03d8312f79832dcf3'
BEGIN = '<!-- BEGIN RECOMPUTED RESULTS -->'
END = '<!-- END RECOMPUTED RESULTS -->'
SCALES = (('0.5', 0.5), ('1', 1.0), ('4/3', 4 / 3), ('2', 2.0), ('5', 5.0))
CURVE_LABELS = {
    'reciprocal_holdout': 'Reciprocal holdout (primary)',
    'reciprocal_all': 'Reciprocal full panel',
    'log_holdout': 'Clipped log holdout',
    'log_all': 'Clipped log full panel',
}
FROZEN_HASHES = {
    (Path('research') / 'merkel-tmb-bin-test' / 'PROTOCOL.md').as_posix():
        '4789e3562692bf563d0eed4754e519b5f5ae45424092516e074117fd7a1ae1d6',
    (Path('research') / 'merkel-viral-antigen-test' / 'PROTOCOL.md').as_posix():
        '22923506eff606d81284a476a56813ec0773d73d238a5dcf633067344e711ea1',
    (Path('research') / 'merkel-viral-antigen-test' / 'analysis.py').as_posix():
        '39c6d7be8b81a02c8b6ee5c9ffacff37df30e9d705a6745039424af4cebcde31',
    (Path('research') / 'merkel-viral-antigen-test' / 'data' / 'frozen_curves.csv').as_posix():
        '1acc4eff8810b2afd59c6e5cac250705f2e4e4f50b1c5f079a0c18308bb893f9',
    (Path('research') / 'immune-escape-pancancer-test' / 'analysis.py').as_posix():
        '06a7e15a24fa23b52f4c74b3659c678dd2778a88fdb429b0712e170e4b564b1c',
    (Path('research') / 'immune-escape-pancancer-test' / 'data' / 'tmb_orr.csv').as_posix():
        'b34c1dd8d8d2c403f7251cf02611283d5326cbdc2caf2c72a6afadcf5ddc540b',
}


@dataclass(frozen=True)
class Interval:
    """Nonempty interval; infinite endpoints are never included."""

    lo: float
    hi: float
    left_closed: bool
    right_closed: bool

    def __post_init__(self):
        if (math.isnan(self.lo) or math.isnan(self.hi) or self.lo > self.hi or
                (self.lo == self.hi and not (self.left_closed and self.right_closed)) or
                (not math.isfinite(self.lo) and self.left_closed) or
                (not math.isfinite(self.hi) and self.right_closed)):
            raise ValueError('Invalid or empty interval')

    def contains(self, value):
        if not math.isfinite(value):
            return False
        return ((value > self.lo or (value == self.lo and self.left_closed)) and
                (value < self.hi or (value == self.hi and self.right_closed)))

    def as_dict(self):
        return {'lo': self.lo if math.isfinite(self.lo) else '-infinity',
                'hi': self.hi if math.isfinite(self.hi) else 'infinity',
                'left_closed': self.left_closed, 'right_closed': self.right_closed}


POSITIVE = Interval(0.0, math.inf, False, False)
BANDS = {'counting_assumption': Interval(1.0, 2.0, True, True),
         'wide_assumption': Interval(0.5, 5.0, True, True)}


def normalized_bytes(path):
    return Path(path).read_bytes().decode('utf-8').replace('\r\n', '\n').replace(
        '\r', '\n').encode('utf-8')


def digest(path):
    return hashlib.sha256(normalized_bytes(path)).hexdigest()


def read_csv(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def file_hashes():
    paths = [ROOT / key for key in FROZEN_HASHES]
    paths += [DATA, HERE / 'analysis.py',
              ROOT / 'tests' / 'test_research_merkel_tmb_bin_test.py']
    return {path.relative_to(ROOT).as_posix(): digest(path)
            for path in sorted(set(paths), key=lambda p: p.as_posix())}


def load_inputs(path=DATA):
    rows = read_csv(path)
    required = ('id', 'part', 'bin', 'tmb_rule', 'cutoff_nssv_mb',
                'orr_percent', 'ci_low_percent', 'ci_high_percent',
                'tmb_evaluable_cohort_n', 'population', 'endpoint', 'source_title',
                'doi', 'source_url', 'html_url', 'section', 'paragraph',
                'denominator_locator', 'access_date', 'access_timezone',
                'reuse_basis', 'license_url', 'inspection', 'limitations', 'notes')
    seen = set()
    totals = {}
    for row in rows:
        if any(not row.get(key, '').strip() for key in required):
            raise ValueError('Missing aggregate value or row provenance')
        key = (row['part'], row['bin'])
        if key in seen or key not in {(p, b) for p in ('A', 'B') for b in ('low', 'high')}:
            raise ValueError('Duplicate or unexpected cohort/bin')
        seen.add(key)
        if row['id'] != '_'.join(key):
            raise ValueError('Unexpected bin identifier')
        low, point, high = [float(row[k]) for k in
                            ('ci_low_percent', 'orr_percent', 'ci_high_percent')]
        if not 0 <= low <= point <= high <= 100:
            raise ValueError('Invalid published percentage interval')
        if (float(row['cutoff_nssv_mb']) != 2 or
                row['tmb_rule'] != ('<2' if row['bin'] == 'low' else '>=2')):
            raise ValueError('Unexpected TMB bin definition')
        n = int(row['tmb_evaluable_cohort_n'])
        if n <= 0 or (row['part'] in totals and totals[row['part']] != n):
            raise ValueError('Inconsistent TMB-evaluable cohort denominator')
        totals[row['part']] = n
        # The selected prose supplies neither field. Never recover them from ORRs.
        if row.get('bin_n') != '' or row.get('responders') != '':
            raise ValueError('Unprinted bin denominator or responder count')
    if len(seen) != 4:
        raise ValueError('Both bins in both cohorts are required')
    return sorted(rows, key=lambda r: (r['part'], r['bin'] == 'high'))


def load_curves():
    for key, expected in FROZEN_HASHES.items():
        if digest(ROOT / key) != expected:
            raise ValueError('Frozen source fingerprint mismatch: ' + key)
    found = {}
    for row in read_csv(CURVES_PATH):
        name = row['name']
        if name in found or name not in CURVE_LABELS:
            raise ValueError('Unexpected frozen curve')
        kind = row['kind']
        if kind not in ('reciprocal', 'log'):
            raise ValueError('Unknown frozen curve family')
        curve = {'kind': kind, 'n': int(row['n']), 'parents': int(row['parents'])}
        keys = ('A', 'B') if kind == 'reciprocal' else ('a', 'b')
        curve.update({k: float(row[k]) for k in keys})
        if any(not math.isfinite(curve[k]) for k in keys):
            raise ValueError('Nonfinite coefficient')
        if kind == 'reciprocal' and not (0 < curve['A'] <= 1 and curve['B'] > 0):
            raise ValueError('Reciprocal curve is outside the frozen family')
        if kind == 'log' and curve['b'] <= 0:
            raise ValueError('Log curve is not increasing')
        found[name] = curve
    if set(found) != set(CURVE_LABELS):
        raise ValueError('Missing frozen curve')
    return {name: found[name] for name in CURVE_LABELS}


def intersect(*regions):
    if not regions:
        return POSITIVE
    if any(region is None for region in regions):
        return None
    lo = max(region.lo for region in regions)
    hi = min(region.hi for region in regions)
    left = all(r.left_closed for r in regions if r.lo == lo)
    right = all(r.right_closed for r in regions if r.hi == hi)
    if lo > hi or (lo == hi and not (left and right)):
        return None
    return Interval(lo, hi, left, right)


def predict(curve, x):
    if x < 0 or not math.isfinite(x):
        raise ValueError('Curve input must be finite and nonnegative')
    if x == 0:
        return 0.0
    if curve['kind'] == 'reciprocal':
        return curve['A'] * math.exp(-curve['B'] / x)
    return min(1.0, max(0.0, curve['a'] + curve['b'] * math.log(x)))


def attainable(curve, bin_name, scale, cutoff=2.0):
    """Exact mathematical response ranges, with numerical boundary values.

    The inverse scale calculation, rather than epsilon-expanded forward values,
    decides compatibility. Finite-grid forward values are ordinary floats.
    """
    if not math.isfinite(scale) or scale <= 0 or not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError('Positive finite scale and cutoff required')
    if bin_name not in ('low', 'high'):
        raise ValueError('Unknown bin')
    if curve['kind'] == 'reciprocal':
        s = predict(curve, scale * cutoff)
        if s <= 0 or s >= curve['A']:
            raise ArithmeticError('Forward range exceeds floating-point resolution')
        return (Interval(0.0, s, True, False) if bin_name == 'low'
                else Interval(s, curve['A'], True, False))
    c0 = math.exp(-curve['a'] / curve['b']) / cutoff
    c1 = math.exp((1 - curve['a']) / curve['b']) / cutoff
    # Explicit thresholds prevent a rounded log value from changing attainment.
    if scale <= c0:
        s = 0.0
    elif scale >= c1:
        s = 1.0
    else:
        s = predict(curve, scale * cutoff)
    if bin_name == 'high':
        return Interval(s, 1.0, True, True)
    if scale <= c0:
        return Interval(0.0, 0.0, True, True)
    return Interval(0.0, s, True, scale > c1)


def scale_region(curve, bin_name, lower, upper, cutoff=2.0):
    """Sharp scale set for arbitrary finite-T distributions, no optimizer."""
    if not 0 <= lower <= upper <= 1 or not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError('Invalid target interval or cutoff')
    if bin_name not in ('low', 'high'):
        raise ValueError('Unknown bin')
    if curve['kind'] == 'reciprocal':
        a, b = curve['A'], curve['B']
        if lower >= a:
            return None  # Full response ceiling, not merely the threshold bound.
        if bin_name == 'low':
            return (POSITIVE if lower == 0 else
                    Interval(b / math.log(a / lower) / cutoff, math.inf, False, False))
        if upper <= 0:
            return None
        return (POSITIVE if upper >= a else
                Interval(0.0, b / math.log(a / upper) / cutoff, False, True))
    a, b = curve['a'], curve['b']
    if bin_name == 'low':
        return (POSITIVE if lower == 0 else
                Interval(math.exp((lower - a) / b) / cutoff, math.inf, False, False))
    return (POSITIVE if upper == 1 else
            Interval(0.0, math.exp((upper - a) / b) / cutoff, False, True))


def half_unit(value):
    """Half the last printed decimal unit; trailing zeros are significant."""
    return Decimal(1).scaleb(Decimal(value).as_tuple().exponent) / 2


def target(row, mode='interval'):
    if mode == 'point':
        value = float(Decimal(row['orr_percent']) / 100)
        return value, value
    low, high = (Decimal(row[k]) for k in ('ci_low_percent', 'ci_high_percent'))
    if mode == 'rounded':
        low = max(Decimal(0), low - half_unit(row['ci_low_percent']))
        high = min(Decimal(100), high + half_unit(row['ci_high_percent']))
    elif mode != 'interval':
        raise ValueError('Unknown scoring mode')
    return float(low / 100), float(high / 100)


def distance_infimum(value, region):
    """Distance to the closure; zero does not imply an attainable value."""
    return max(region.lo - value, value - region.hi, 0.0)


def witness(curve, row, scale):
    """One finite point-mass existence certificate, not an observed record.

    No witness is retained in the source CSV and none estimates a distribution.
    Choosing a feasible mean here cannot change any scale set or score.
    """
    low, high = target(row)
    cutoff = float(row['cutoff_nssv_mb'])
    scales = scale_region(curve, row['bin'], low, high, cutoff)
    # A rounded forward evaluation must never close an analytic open endpoint.
    if scales is None or not scales.contains(scale):
        return None
    region = attainable(curve, row['bin'], scale, cutoff)
    overlap = intersect(region, Interval(low, high, True, True))
    if overlap is None:
        raise ArithmeticError('Witness overlap exceeds floating-point resolution')
    mean = (overlap.lo + overlap.hi) / 2
    if row['bin'] == 'high' and mean == region.lo:
        tmb = cutoff
    elif mean == 0:
        tmb = 0.0
    elif curve['kind'] == 'reciprocal':
        tmb = curve['B'] / math.log(curve['A'] / mean) / scale
    else:
        tmb = math.exp((mean - curve['a']) / curve['b']) / scale
        if row['bin'] == 'high':
            tmb = max(cutoff, tmb)
    in_bin = (0 <= tmb < cutoff if row['bin'] == 'low' else cutoff <= tmb)
    if not math.isfinite(tmb) or not in_bin:
        raise ArithmeticError('Finite witness exceeds floating-point resolution')
    return {'mean': mean, 'tmb': tmb}


def build_report(rows=None, curves=None):
    rows = load_inputs() if rows is None else rows
    curves = load_curves() if curves is None else curves
    modes = {}
    for mode in ('interval', 'point', 'rounded'):
        models = {}
        for name, curve in curves.items():
            bins = {row['id']: scale_region(curve, row['bin'], *target(row, mode),
                                           float(row['cutoff_nssv_mb'])) for row in rows}
            cohorts = {part: intersect(*(bins[r['id']] for r in rows if r['part'] == part))
                       for part in ('A', 'B')}
            joint = intersect(*cohorts.values())
            models[name] = {'bins': bins, 'cohorts': cohorts, 'joint': joint,
                            'bands': {band: intersect(joint, limits)
                                      for band, limits in BANDS.items()}}
        modes[mode] = models
    grid = {}
    for name, curve in curves.items():
        grid[name] = []
        for label, scale in SCALES:
            ranges = {r['id']: attainable(curve, r['bin'], scale,
                                         float(r['cutoff_nssv_mb'])) for r in rows}
            passes = {mode: sum(region is not None and region.contains(scale)
                                for region in modes[mode][name]['bins'].values())
                      for mode in ('interval', 'point')}
            gaps = {r['id']: 100 * distance_infimum(target(r, 'point')[0], ranges[r['id']])
                    for r in rows}
            grid[name].append({'scale_label': label, 'scale': scale, 'ranges': ranges,
                               'interval_passes': passes['interval'],
                               'point_passes': passes['point'], 'point_gaps_pp': gaps,
                               'rms_gap_pp': math.sqrt(sum(v * v for v in gaps.values()) / len(rows)),
                               'witnesses': {r['id']: witness(curve, r, scale) for r in rows}})
    return {'inputs': rows, 'curves': curves, 'modes': modes, 'grid': grid,
            'hashes': file_hashes()}


def fmt(region, percent=False):
    if region is None:
        return 'empty'
    multiplier = 100 if percent else 1
    lo = '-infinity' if not math.isfinite(region.lo) else f'{multiplier * region.lo:.6f}'
    hi = 'infinity' if not math.isfinite(region.hi) else f'{multiplier * region.hi:.6f}'
    return ('[' if region.left_closed else '(') + lo + ', ' + hi + (
        ']' if region.right_closed else ')')


def render_report(report):
    rows, curves, modes = report['inputs'], report['curves'], report['modes']
    lines = ['### Transcribed outcome-matched aggregates', '',
             'All percentages and intervals below are printed source values. The n column is the '
             '**total TMB-evaluable cohort**, repeated for context; it is not either bin denominator. '
             'Bin denominators and responder counts are not printed in the selected prose and remain blank in the CSV.', '',
             '| Part | TMB bin (NSSV/Mb) | ORR % | Printed marginal 95% interval % | TMB-evaluable cohort n |',
             '| --- | --- | ---: | --- | ---: |']
    for row in rows:
        lines.append(f"| {row['part']} | {row['tmb_rule']} | {row['orr_percent']} | "
                     f"[{row['ci_low_percent']}, {row['ci_high_percent']}] | {row['tmb_evaluable_cohort_n']} |")
    lines += ['', '### Frozen curve identity and full response range', '',
              'Coefficients are the unchanged full-precision PR #72 CSV values. Its original fit replay '
              'is checked with absolute coefficient tolerance 2e-6; optimizer diagnostics are bounded, '
              'not compared with exact saved outputs. These curves are never fitted to the new bins.', '',
              '| Curve | Coefficients | Training rows / parents | Response range for finite nonnegative input |',
              '| --- | --- | --- | --- |']
    for name, curve in curves.items():
        keys = ('A', 'B') if curve['kind'] == 'reciprocal' else ('a', 'b')
        coefficients = ', '.join(f'{key}={curve[key]!r}' for key in keys)
        response_range = f"[0, {100 * curve['A']:.6f}%)" if curve['kind'] == 'reciprocal' else '[0%, 100%]'
        lines.append(f"| {CURVE_LABELS[name]} | {coefficients} | {curve['n']} / {curve['parents']} | {response_range} |")
    lines += ['', '### Primary printed-interval scale regions', '',
              'All sets concern finite positive c. Parentheses are scientifically meaningful: a low-bin '
              'mean cannot attain the threshold response while the curve is increasing. Values are displayed '
              'to six decimals; endpoint flags apply to the unrounded analytic values. These sets do not '
              'have simultaneous confidence coverage.', '',
              '| Curve | Part | Low-bin scales | High-bin scales | Both bins |',
              '| --- | --- | --- | --- | --- |']
    for name in curves:
        result = modes['interval'][name]
        for part in ('A', 'B'):
            lines.append(f"| {CURVE_LABELS[name]} | {part} | {fmt(result['bins'][part + '_low'])} | "
                         f"{fmt(result['bins'][part + '_high'])} | {fmt(result['cohorts'][part])} |")
    lines += ['', '| Curve | Both cohorts, all c>0 | In assumed [1,2] | In assumed [0.5,5] | Primary rule |',
              '| --- | --- | --- | --- | --- |']
    for name in curves:
        result = modes['interval'][name]
        verdict = 'survives with stated latitude' if result['joint'] is not None else 'bridge falsified under this rule'
        lines.append(f"| {CURVE_LABELS[name]} | {fmt(result['joint'])} | "
                     f"{fmt(result['bands']['counting_assumption'])} | "
                     f"{fmt(result['bands']['wide_assumption'])} | {verdict} |")
    for mode, title, explanation in (
            ('point', 'Exact printed point-ORR sensitivity',
             'Points are treated as exact means only in this secondary diagnostic. A point failure is '
             'not a population rejection when the printed-interval rule survives.'),
            ('rounded', 'Source-rounding sensitivity',
             'Only interval endpoints are widened by half their last printed unit (0.05 percentage '
             'points here). Thresholds and coefficients are unchanged. This is not added sampling uncertainty.')):
        lines += ['', '### ' + title, '', explanation, '',
                  '| Curve | Part A | Part B | Both cohorts | In [1,2] | In [0.5,5] |',
                  '| --- | --- | --- | --- | --- | --- |']
        for name in curves:
            result = modes[mode][name]
            lines.append(f"| {CURVE_LABELS[name]} | {fmt(result['cohorts']['A'])} | "
                         f"{fmt(result['cohorts']['B'])} | {fmt(result['joint'])} | "
                         f"{fmt(result['bands']['counting_assumption'])} | {fmt(result['bands']['wide_assumption'])} |")
        lines += ['', 'Per-bin sets (including ceiling failures):', '',
                  '| Curve | A low | A high | B low | B high |',
                  '| --- | --- | --- | --- | --- |']
        for name in curves:
            values = ' | '.join(fmt(modes[mode][name]['bins'][r['id']]) for r in rows)
            lines.append(f'| {CURVE_LABELS[name]} | {values} |')
    lines += ['', '### Scale-independent reciprocal ceiling deficits', '',
              'Positive deficits are infima at infinite load: the finite-T ceiling is excluded. '
              'Zero deficits here do not decide low-bin feasibility.', '',
              '| Reciprocal curve | A high point above ceiling (pp) | B high point above ceiling (pp) |',
              '| --- | ---: | ---: |']
    for name, curve in curves.items():
        if curve['kind'] != 'reciprocal':
            continue
        deficits = [100 * max(0, target(row, 'point')[0] - curve['A'])
                    for row in rows if row['bin'] == 'high']
        lines.append(f'| {CURVE_LABELS[name]} | {deficits[0]:.6f} | {deficits[1]:.6f} |')
    lines += ['', '### Every inherited fixed scale: attainable means and feasibility', '',
              'Ranges are percentages, identical across cohorts at a given scale because the cutoff is '
              'shared. Counts are numbers of compatible bins out of the four published bins. Each bin '
              'may use a different admissible distribution.', '',
              '| Curve | c | Low mean range % | High mean range % | Interval-compatible bins | Point-compatible bins |',
              '| --- | ---: | --- | --- | ---: | ---: |']
    for name in curves:
        for result in report['grid'][name]:
            lines.append(f"| {CURVE_LABELS[name]} | {result['scale_label']} | "
                         f"{fmt(result['ranges']['A_low'], True)} | {fmt(result['ranges']['A_high'], True)} | "
                         f"{result['interval_passes']} | {result['point_passes']} |")
    lines += ['', '### Point-distance infima at every fixed scale', '',
              'Distances in percentage points grant arbitrary within-bin distributions. RMS is the root '
              'mean of the four squared distances. These are lower bounds on aggregate discrepancy, '
              '**not fitted prediction errors**. An excluded endpoint can have zero distance but fail '
              'exact compatibility, which is evaluated separately above.', '',
              '| Curve | c | A low gap | A high gap | B low gap | B high gap | RMS gap |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for name in curves:
        for result in report['grid'][name]:
            gaps = ' | '.join(f"{result['point_gaps_pp'][r['id']]:.6f}" for r in rows)
            lines.append(f"| {CURVE_LABELS[name]} | {result['scale_label']} | {gaps} | {result['rms_gap_pp']:.6f} |")
    lines += ['', '### Canonical input and implementation fingerprints', '',
              'SHA-256 of UTF-8 after CRLF/CR normalization to LF. Every path key is built with '
              '`Path.as_posix()`. README is excluded from its own fingerprint. Original source files '
              'are read without modification; the protocol hash is pinned in code.', '',
              '| Path | SHA-256 |', '| --- | --- |']
    for path, value in report['hashes'].items():
        lines.append(f'| `{path}` | `{value}` |')
    return '\n'.join(lines) + '\n'


def jsonable(value):
    if isinstance(value, Interval):
        return value.as_dict()
    if isinstance(value, dict):
        return {key: jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true', help='Check README generated results without writes')
    group.add_argument('--write', action='store_true', help='Update only the generated README block')
    group.add_argument('--json', action='store_true', help='Print all sets, scores and mathematical witnesses')
    args = parser.parse_args()
    report = build_report()
    if args.json:
        print(json.dumps(jsonable(report), indent=2, allow_nan=False))
        return
    block = render_report(report)
    if args.check or args.write:
        path = HERE / 'README.md'
        old = normalized_bytes(path).decode('utf-8')
        if old.count(BEGIN) != 1 or old.count(END) != 1:
            raise SystemExit('Expected one generated block')
        before, tail = old.split(BEGIN)
        unused, after = tail.split(END)
        new = before + BEGIN + '\n\n' + block + '\n' + END + after
        if args.write:
            with path.open('w', encoding='utf-8', newline='\n') as handle:
                handle.write(new)
            print('Updated the generated README block only.')
        elif old != new:
            raise SystemExit('README differs from recomputed source-CSV results')
        else:
            print('README matches all recomputed source-CSV results and normalized fingerprints.')
    else:
        print(block, end='')


if __name__ == '__main__':
    main()
