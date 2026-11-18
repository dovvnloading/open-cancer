#!/usr/bin/env python3
"""Offline, aggregate-only, NONCLINICAL ecological model comparison.

Standard library only. See PROTOCOL.md for the fixed predictions and README.md
for the observation assumptions. Numerical reproducibility is not biological
validation. No network, individual records, or reconstruction of individuals.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import statistics

HERE = Path(__file__).resolve().parent
BEGIN = '<!-- BEGIN RECOMPUTED RESULTS -->'
END = '<!-- END RECOMPUTED RESULTS -->'
SEED = 20261117
PARAMETER_REPLICATES = 1000
ERROR_REPLICATES = 10000
BASELINES = ('mean', 'linear_tmb', 'log_tmb', 'reciprocal_tmb')
LABELS = {
    'mean': 'Mean only', 'linear_tmb': 'Linear TMB',
    'log_tmb': 'Log-linear TMB', 'reciprocal_tmb': 'Reciprocal TMB',
    'reciprocal_clonal': 'Reciprocal clonal',
    'reciprocal_hla': 'Reciprocal TMB + HLA',
    'reciprocal_full': 'Full reciprocal',
    'free_hla': 'HLA sign unrestricted',
    'free_full': 'Full HLA sign unrestricted',
    'power_tmb': 'Free TMB exponent',
    'power_clonal': 'Free clonal exponent',
}

# Source-to-target categorical decisions frozen before empirical fitting.
# A join is an ecological measurement assumption, not a patient-level merge.
DASH_MAP = {
    'urothelial': 'dash2022_bladder', 'breast': 'dash2022_breast',
    'cervical': 'dash2022_cervical', 'head_neck': 'dash2022_hnsc',
    'renal_cell': 'dash2022_kidney', 'hepatocellular': 'dash2022_liver',
    'ovarian': 'dash2022_ovarian', 'pancreatic': 'dash2022_pancreatic',
    'skin_melanoma': 'dash2022_melanoma', 'endometrial': 'dash2022_uterine',
}
MONTESION_MAP = {
    'head_neck': 'montesion2021_republished_hnsc',
    'breast': 'montesion2021_republished_breast',
    'pancreatic': 'montesion2021_republished_pancreatic',
    'prostate': 'montesion2021_republished_prostate',
    'ovarian': 'montesion2021_republished_ovarian',
    'nsclc': 'montesion2021_republished_nsclc',
    'mesothelioma': 'montesion2021_republished_mesothelioma',
}
APPROXIMATE_CLONAL = {'urothelial', 'cervical', 'head_neck', 'renal_cell',
                      'pancreatic', 'prostate', 'endometrial'}
APPROXIMATE_HLA = {'urothelial', 'cervical', 'head_neck', 'renal_cell',
                   'hepatocellular', 'ovarian', 'pancreatic', 'skin_melanoma',
                   'endometrial', 'mesothelioma'}


def normalized_bytes(path):
    """Every file checksum uses LF-normalized UTF-8 content, including CSVs."""
    return Path(path).read_bytes().replace(b'\r\n', b'\n').replace(b'\r', b'\n')


def digest(path):
    return hashlib.sha256(normalized_bytes(path)).hexdigest()


def read_csv(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle))
    ids = set()
    for row in rows:
        if not row.get('id') or row['id'] in ids:
            raise ValueError(f'Missing or duplicate id in {path}')
        ids.add(row['id'])
        if 'locator' not in row and 'source_locator' in row:
            row['locator'] = row['source_locator']
        if 'source_url' not in row and 'orr_source_url' in row:
            for field in ('source_url', 'locator', 'access_date', 'reuse_basis'):
                row[field] = row['orr_'+field]
            if row.get('tmb'):
                for field in ('source_url', 'locator', 'access_date', 'reuse_basis'):
                    if not row.get('tmb_'+field, '').strip():
                        raise ValueError(f'Missing TMB {field}: {row["id"]}')
        for field in ('source_url', 'locator', 'access_date', 'reuse_basis'):
            if not row.get(field, '').strip():
                raise ValueError(f'Missing {field}: {row["id"]}')
        if not row['source_url'].startswith('https://'):
            raise ValueError(f'Invalid source URL: {row["id"]}')
    return rows


def clip(x, lo=0.0, hi=1.0):
    return max(lo, min(hi, x))


def quantile(values, p):
    values = sorted(values)
    if not values:
        return None
    pos = (len(values) - 1) * p
    low = int(pos)
    high = min(len(values) - 1, low + 1)
    return values[low] + (pos - low) * (values[high] - values[low])


def interval(values):
    return [quantile(values, .025), quantile(values, .975)]


def grouped(rows):
    result = {}
    for row in rows:
        result.setdefault(row['parent_type'], []).append(row)
    return result


def group_count(rows):
    return len(grouped(rows))


def _bounded_ceiling(z, outcomes, loh=None, adverse=True):
    """Exact bounded linear LS for fixed B, including all boundary faces.

    The adverse model is the triangle 0 <= A1 <= A0 <= 1; its sign-relaxed
    diagnostic is the square [0,1]^2. No numerical optimizer for this step.
    """
    if loh is None:
        zz = math.fsum(x*x for x in z)
        a = clip(math.fsum(x*y for x, y in zip(z, outcomes))/zz) if zz else 0.
        fitted = [a*x for x in z]
        loss = math.fsum((y-p)**2 for y, p in zip(outcomes, fitted))
        return loss, a, a
    u = [x*(1-l) for x, l in zip(z, loh)]
    v = [x*l for x, l in zip(z, loh)]
    uu = math.fsum(x*x for x in u)
    vv = math.fsum(x*x for x in v)
    uv = math.fsum(x*y for x, y in zip(u, v))
    uy = math.fsum(x*y for x, y in zip(u, outcomes))
    vy = math.fsum(x*y for x, y in zip(v, outcomes))
    yy = math.fsum(y*y for y in outcomes)
    det = uu*vv-uv*uv
    candidates = [(0., 0.), (1., 0.), (1., 1.)]
    if det > 1e-13 * max(uu*vv, 1e-300):
        a0, a1 = (uy*vv-vy*uv)/det, (vy*uu-uy*uv)/det
        if 0 <= a0 <= 1 and 0 <= a1 <= 1 and (not adverse or a1 <= a0):
            candidates.append((a0, a1))
    # A1=0; A0=1; the diagonal A0=A1.
    candidates.append((clip(uy/uu) if uu else 0., 0.))
    candidates.append((1., clip((vy-uv)/vv) if vv else 0.))
    diagonal = uu+2*uv+vv
    a = clip((uy+vy)/diagonal) if diagonal else 0.
    candidates.append((a, a))
    if not adverse:
        candidates.extend([(0., 1.), (0., clip(vy/vv) if vv else 0.),
                           (clip((uy-uv)/uu) if uu else 0., 1.)])
    def loss(pair):
        a0, a1 = pair
        return yy-2*(a0*uy+a1*vy)+a0*a0*uu+2*a0*a1*uv+a1*a1*vv
    a0, a1 = min(candidates, key=loss)
    # Direct residuals avoid cancellation near an exact fit.
    error = math.fsum((y-a0*x-a1*w)**2 for y, x, w in zip(outcomes, u, v))
    return error, a0, a1


def _golden(function, left, right, iterations=54):
    ratio = (math.sqrt(5)-1)/2
    c, d = right-ratio*(right-left), left+ratio*(right-left)
    fc, fd = function(c), function(d)
    for _ in range(iterations):
        if fc[0] < fd[0]:
            right, d, fd = d, c, fc
            c = right-ratio*(right-left)
            fc = function(c)
        else:
            left, c, fc = c, d, fd
            d = left+ratio*(right-left)
            fd = function(d)
    return min(((fc, c), (fd, d), (function(left), left),
                (function(right), right)), key=lambda item: item[0][0])


def reciprocal_fit(loads, outcomes, loh=None, adverse=True, dense=False):
    if not loads or len(loads) != len(outcomes):
        raise ValueError('Need paired nonempty aggregate loads and outcomes')
    if any(x <= 0 or not math.isfinite(x) for x in loads):
        raise ValueError('Functional-load proxy must be positive')
    if any(not 0 <= y <= 1 for y in outcomes):
        raise ValueError('ORR must be a probability')
    if loh is not None and (len(loh) != len(loads) or any(not 0 <= x <= 1 for x in loh)):
        raise ValueError('LOH must be a paired probability')
    # B=0 is explicitly allowed. At B>=50*max(load), predictions <= exp(-50).
    # This is the zero-prediction limit, not an attractive hidden high-B minimum.
    scale = statistics.median(loads)
    maximum = 50*max(loads)
    ngrid = 192 if dense else 48
    loglo, loghi = math.log(scale*1e-7), math.log(maximum)
    grid = [0.] + [math.exp(loglo+(loghi-loglo)*i/ngrid) for i in range(ngrid+1)]
    def objective(b):
        return _bounded_ceiling([math.exp(-b/x) for x in loads], outcomes, loh, adverse)
    results = [objective(b) for b in grid]
    best = min(zip(results, grid), key=lambda pair: pair[0][0])
    # Refine every grid-local minimum, not just a chosen starting value.
    for i in range(1, len(grid)-1):
        if (results[i][0] <= results[i-1][0] and results[i][0] <= results[i+1][0]
                and (results[i][0] < results[i-1][0] or results[i][0] < results[i+1][0])):
            candidate = _golden(objective, grid[i-1], grid[i+1])
            if candidate[0][0] < best[0][0]:
                best = candidate
    (loss, a0, a1), b = best
    singular = (max(loads)-min(loads) <= 1e-10*max(loads) or not any(outcomes)
                or (loh is not None and max(loh)-min(loh) <= 1e-10))
    return {'A0': a0, 'A1': a1, 'B': b, 'sse': loss, 'singular': singular,
            'eta': 1-a1/a0 if a0 > 1e-14 else None}


def fit(rows, model, dense=False, clonal_with_hla=None):
    if not rows:
        raise ValueError('Cannot fit an empty panel')
    y = [r['orr'] for r in rows]
    if model in ('mean', 'linear_tmb', 'log_tmb'):
        x = [0. if model == 'mean' else
             (math.log(r['tmb']) if model == 'log_tmb' else r['tmb']) for r in rows]
        mx, my = statistics.mean(x), statistics.mean(y)
        xx = math.fsum((v-mx)**2 for v in x)
        slope = math.fsum((a-mx)*(b-my) for a, b in zip(x, y))/xx if xx > 1e-20 else 0.
        return {'model': model, 'intercept': my-slope*mx, 'slope': slope,
                'singular': model != 'mean' and xx <= 1e-20}
    if model == 'power_tmb':
        def profile(q):
            result = reciprocal_fit([r['tmb']**q for r in rows], y)
            return result['sse'], result
        grid = [.05+3.95*i/24 for i in range(25)]
        profiles = [profile(q) for q in grid]
        best = min(zip(profiles, grid), key=lambda item: item[0][0])
        for i in range(1, len(grid)-1):
            if profiles[i][0] <= profiles[i-1][0] and profiles[i][0] <= profiles[i+1][0]:
                value = _golden(profile, grid[i-1], grid[i+1], 28)
                if value[0][0] < best[0][0]:
                    best = value
        result = dict(best[0][1], q=best[1], model=model)
        return result
    if model == 'power_clonal':
        if clonal_with_hla is None:
            clonal_with_hla = all(r.get('loh_fraction') is not None for r in rows)
        if clonal_with_hla and any(r.get('loh_fraction') is None for r in rows):
            raise ValueError('Fixed HLA diagnostic design has missing HLA')
        loh = [r['loh_fraction'] for r in rows] if clonal_with_hla else None
        def profile(theta):
            result = reciprocal_fit([r['tmb']*r['fraction']**theta for r in rows], y, loh)
            return result['sse'], result
        grid = [-2+5*i/20 for i in range(21)]
        profiles = [profile(q) for q in grid]
        best = min(zip(profiles, grid), key=lambda item: item[0][0])
        for i in range(1, len(grid)-1):
            if profiles[i][0] <= profiles[i-1][0] and profiles[i][0] <= profiles[i+1][0]:
                value = _golden(profile, grid[i-1], grid[i+1], 28)
                if value[0][0] < best[0][0]:
                    best = value
        return dict(best[0][1], theta=best[1], with_hla=loh is not None, model=model)
    has_f = model in ('reciprocal_clonal', 'reciprocal_full', 'free_full')
    has_l = model in ('reciprocal_hla', 'reciprocal_full', 'free_hla', 'free_full')
    if model not in LABELS:
        raise ValueError(f'Unknown model {model}')
    if has_f and any(r.get('fraction') is None for r in rows):
        raise ValueError('Missing clonal fraction cannot be imputed')
    if has_l and any(r.get('loh_fraction') is None for r in rows):
        raise ValueError('Missing HLA-LOH cannot be imputed')
    loads = [r['tmb']*(r['fraction'] if has_f else 1.) for r in rows]
    loh = [r['loh_fraction'] for r in rows] if has_l else None
    result = reciprocal_fit(loads, y, loh, not model.startswith('free_'), dense=dense)
    if has_l:
        result['ceiling_difference'] = result['A1']-result['A0']
    return dict(result, model=model)


def predict(parameters, row):
    model = parameters['model']
    if model in ('mean', 'linear_tmb', 'log_tmb'):
        x = 0. if model == 'mean' else (math.log(row['tmb']) if model == 'log_tmb' else row['tmb'])
        return clip(parameters['intercept']+parameters['slope']*x)
    if model == 'power_tmb':
        load = row['tmb']**parameters['q']
    elif model == 'power_clonal':
        load = row['tmb']*row['fraction']**parameters['theta']
    else:
        load = row['tmb']*(row['fraction'] if model in ('reciprocal_clonal', 'reciprocal_full', 'free_full') else 1.)
    has_l = model in ('reciprocal_hla', 'reciprocal_full', 'free_hla', 'free_full') or parameters.get('with_hla', False)
    l = row['loh_fraction'] if has_l else 0.
    ceiling = parameters['A0']*(1-l)+parameters['A1']*l
    return ceiling*math.exp(-parameters['B']/load)


def scores(rows, predictions):
    if len(rows) != len(predictions) or not rows:
        raise ValueError('Need one prediction per nonempty aggregate row')
    errors = {}
    for row, predicted in zip(rows, predictions):
        errors.setdefault(row['parent_type'], []).append(row['orr']-predicted)
    squared = {g: statistics.mean(e*e for e in values) for g, values in errors.items()}
    absolute = {g: statistics.mean(abs(e) for e in values) for g, values in errors.items()}
    bias = statistics.mean(statistics.mean(e) for e in errors.values())
    mse = statistics.mean(squared.values())
    return {'mse': mse, 'rmse': math.sqrt(mse), 'mae': statistics.mean(absolute.values()),
            'mean_residual': bias, 'group_mse': squared,
            'row_rmse': math.sqrt(statistics.mean((r['orr']-p)**2 for r, p in zip(rows, predictions)))}


def cross_validate(rows, model):
    groups = grouped(rows)
    if len(groups) < 3:
        raise ValueError('Need at least three parent groups for held-out prediction')
    predictions, folds = {}, {}
    for group, test in groups.items():
        train = [r for r in rows if r['parent_type'] != group]
        parameters = fit(train, model)
        folds[group] = parameters
        for row in test:
            predictions[row['id']] = predict(parameters, row)
    ordered = [predictions[r['id']] for r in rows]
    return dict(scores(rows, ordered), predictions=predictions, folds=folds)


def evaluate(rows, model):
    parameters = fit(rows, model)
    fitted = [predict(parameters, r) for r in rows]
    sse = math.fsum((r['orr']-p)**2 for r, p in zip(rows, fitted))
    mean = statistics.mean(r['orr'] for r in rows)
    tss = math.fsum((r['orr']-mean)**2 for r in rows)
    return {'parameters': parameters, 'training_sse': sse,
            'training_r2': 1-sse/tss if tss else None,
            'training_predictions': dict(zip((r['id'] for r in rows), fitted)),
            'cv': cross_validate(rows, model)}


def paired_error_interval(first, second, replicates=ERROR_REPLICATES):
    keys = sorted(first['group_mse'])
    if keys != sorted(second['group_mse']):
        raise ValueError('Prediction comparisons require identical parent groups')
    if ('predictions' in first and 'predictions' in second
            and set(first['predictions']) != set(second['predictions'])):
        raise ValueError('Prediction comparisons require identical aggregate rows')
    delta = [first['group_mse'][k]-second['group_mse'][k] for k in keys]
    rng = random.Random(SEED)
    boot = [statistics.mean(rng.choices(delta, k=len(delta))) for _ in range(replicates)]
    return {'delta_mse': statistics.mean(delta), 'interval': interval(boot),
            'replicates': replicates, 'groups': len(keys)}


def parameter_bootstrap(rows, model, replicates=PARAMETER_REPLICATES):
    groups = grouped(rows)
    keys = sorted(groups)
    rng = random.Random(SEED)
    values = {}
    skipped = 0
    boundaries = 0
    for _ in range(replicates):
        draws = rng.choices(keys, k=len(keys))
        if len(set(draws)) < 3:
            skipped += 1
            continue
        sample = [row for key in draws for row in groups[key]]
        result = fit(sample, model)
        if result.get('singular', False):
            skipped += 1
            continue
        for key in ('A0', 'A1', 'B', 'eta', 'ceiling_difference', 'intercept', 'slope'):
            if key in result and result[key] is not None:
                values.setdefault(key, []).append(result[key])
        at_boundary = result.get('B', 1) < 1e-8 or result.get('A0', .5) < 1e-8 or result.get('A0', .5) > 1-1e-8
        if model in ('reciprocal_hla', 'reciprocal_full', 'free_hla', 'free_full'):
            at_boundary = at_boundary or result['A1'] < 1e-8 or result['A1'] > 1-1e-8
            if not model.startswith('free_'):
                at_boundary = at_boundary or abs(result['A0']-result['A1']) < 1e-8
        if at_boundary:
            boundaries += 1
    return {'replicates': replicates, 'skipped': skipped, 'boundary_fits': boundaries,
            'usable': {key: len(value) for key, value in values.items()},
            'intervals': {key: interval(value) for key, value in values.items()}}


def available_panels(rows):
    clonal = [r for r in rows if r.get('fraction') is not None]
    hla = [r for r in rows if r.get('loh_fraction') is not None]
    full = [r for r in clonal if r.get('loh_fraction') is not None]
    specs = [('all', rows, 3, 'reciprocal_tmb'),
             ('clonal', clonal, 5, 'reciprocal_clonal'),
             ('hla', hla, 10, 'reciprocal_hla'),
             ('full', full, 10, 'reciprocal_full')]
    return [{'name': name, 'rows': panel, 'gate': gate, 'model': model,
             'eligible': group_count(panel) >= gate} for name, panel, gate, model in specs]


def analyze_rows(rows, bootstrap=True, diagnostics=True):
    if len({r['id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate aggregate row')
    panels = []
    for spec in available_panels(rows):
        output = {k: v for k, v in spec.items() if k != 'rows'}
        output['row_ids'] = [r['id'] for r in spec['rows']]
        output['n_rows'] = len(spec['rows'])
        output['n_groups'] = group_count(spec['rows'])
        if spec['eligible']:
            models = list(BASELINES)
            if spec['name'] == 'full':
                models.extend(['reciprocal_clonal', 'reciprocal_hla'])
            if spec['model'] not in models:
                models.append(spec['model'])
            output['models'] = {model: evaluate(spec['rows'], model) for model in models}
            theory = output['models'][spec['model']]['cv']
            output['comparisons'] = {
                model: paired_error_interval(theory, output['models'][model]['cv'])
                for model in ('log_tmb', 'linear_tmb', 'reciprocal_tmb')
                if model != spec['model']}
            if spec['name'] == 'full':
                output['comparisons'].update({
                    model: paired_error_interval(theory, output['models'][model]['cv'])
                    for model in ('reciprocal_clonal', 'reciprocal_hla')})
            if bootstrap:
                output['parameter_uncertainty'] = parameter_bootstrap(spec['rows'], spec['model'])
            if diagnostics and spec['name'] in ('hla', 'full'):
                diagnostic = 'free_full' if spec['name'] == 'full' else 'free_hla'
                output['sign_diagnostic'] = evaluate(spec['rows'], diagnostic)
                if bootstrap:
                    output['sign_uncertainty'] = parameter_bootstrap(spec['rows'], diagnostic)
        panels.append(output)
    primary = next((p for p in panels if p['name'] == 'full' and p['eligible']), panels[0])
    summary = {'panels': panels, 'primary_panel': primary['name']}
    if diagnostics:
        power = fit(rows, 'power_tmb')
        summary['power_diagnostic'] = {'parameters': power,
            'jackknife_q': {g: fit([r for r in rows if r['parent_type'] != g], 'power_tmb')['q'] for g in sorted(grouped(rows))}}
        clonal = [r for r in rows if r.get('fraction') is not None]
        if group_count(clonal) >= 10:
            fixed_hla = all(r.get('loh_fraction') is not None for r in clonal)
            power = fit(clonal, 'power_clonal', clonal_with_hla=fixed_hla)
            summary['clonal_power_diagnostic'] = {'parameters': power,
                'jackknife_theta': {g: fit([r for r in clonal if r['parent_type'] != g],
                    'power_clonal', clonal_with_hla=fixed_hla)['theta'] for g in sorted(grouped(clonal))}}
    return summary


def source_tables(directory=HERE / 'data'):
    directory = Path(directory)
    return {path.name: read_csv(path) for path in sorted(directory.glob('*.csv'))}


def load_rows(directory=HERE / 'data', mixed_hla=False, figure_values=False):
    tables = source_tables(directory)
    clonal = {r['tcga_code']: r for r in tables['clonal_fractions.csv']}
    hla = {r['id']: r for r in tables['hla_loh.csv']}
    result, excluded = [], []
    for raw in tables['tmb_orr.csv']:
        if not raw['tmb']:
            excluded.append(raw)
            continue
        row = dict(raw)
        for field in ('tmb', 'orr', 'tmb_half_width', 'orr_half_width'):
            row[field] = float(raw[field])
            if not math.isfinite(row[field]):
                raise ValueError('Nonfinite source value')
        if row['tmb'] <= 0 or not 0 <= row['orr'] <= 1:
            raise ValueError('Invalid TMB or response probability')
        source_f = clonal.get(raw['tcga_code'])
        row['fraction'] = float(source_f['fraction']) if source_f and source_f['fraction'] else None
        row['fraction_half_width'] = float(source_f['fraction_half_width']) if row['fraction'] is not None else None
        row['clonal_source_id'] = source_f['id'] if source_f else None
        row['clonal_mapping_quality'] = 'approximate' if row['id'] in APPROXIMATE_CLONAL else 'exact'
        if row['fraction'] is not None and not 0 < row['fraction'] <= 1:
            raise ValueError('Invalid clonal fraction')
        key = MONTESION_MAP.get(row['id']) if mixed_hla else None
        key = key or DASH_MAP.get(row['id'])
        if figure_values and key and key+'_figure' in hla:
            key += '_figure'
        source_l = hla.get(key)
        row['loh_fraction'] = float(source_l['loh_fraction']) if source_l else None
        row['loh_half_width'] = float(source_l['loh_half_width']) if source_l else None
        row['hla_source_id'] = key
        row['hla_mapping_quality'] = 'approximate' if row['id'] in APPROXIMATE_HLA else 'exact'
        for field in ('source_ci95_low', 'source_ci95_high', 'source_ci_reading_half_width'):
            row[field] = float(source_l[field]) if source_l and source_l.get(field) else None
        if row['loh_fraction'] is not None and not 0 <= row['loh_fraction'] <= 1:
            raise ValueError('Invalid HLA-LOH fraction')
        result.append(row)
    return result, excluded, tables


def panel_analysis(rows, name, model, bootstrap=True, status='exploratory'):
    """Below-gate and mixed-source additions remain explicitly exploratory."""
    models = list(BASELINES)
    if model == 'reciprocal_full':
        models += ['reciprocal_clonal', 'reciprocal_hla']
    if model not in models:
        models.append(model)
    output = {'name': name, 'model': model, 'status': status,
              'row_ids': [r['id'] for r in rows], 'n_rows': len(rows),
              'n_groups': group_count(rows),
              'models': {m: evaluate(rows, m) for m in models}}
    theory = output['models'][model]['cv']
    output['comparisons'] = {m: paired_error_interval(theory, output['models'][m]['cv'])
                             for m in models if m not in ('mean', model)}
    if bootstrap:
        output['parameter_uncertainty'] = parameter_bootstrap(rows, model)
    if model in ('reciprocal_full', 'reciprocal_hla'):
        free = 'free_full' if model == 'reciprocal_full' else 'free_hla'
        output['sign_diagnostic'] = evaluate(rows, free)
        if bootstrap:
            output['sign_uncertainty'] = parameter_bootstrap(rows, free)
    return output


def comparison_summary(rows, model):
    a, b = evaluate(rows, model), evaluate(rows, 'log_tmb')
    return {'n_rows': len(rows), 'n_groups': group_count(rows),
            'theory_rmse': a['cv']['rmse'], 'log_rmse': b['cv']['rmse'],
            'delta_mse': a['cv']['mse']-b['cv']['mse'],
            'A0': a['parameters']['A0'], 'A1': a['parameters']['A1'],
            'B': a['parameters']['B'], 'eta': a['parameters']['eta']}


def sensitivity_analysis(rows, model, envelope_draws=100):
    """Prespecified deterministic extraction and case-influence sensitivities.

    Envelope ranges are not uncertainty intervals or exhaustive worst cases.
    HLA source-CI scenarios do not convert marginal intervals to a joint CI.
    """
    answer = {}
    exact = [r for r in rows if r['mapping_quality'] == 'exact'
             and (model not in ('reciprocal_full', 'reciprocal_clonal') or r['clonal_mapping_quality'] == 'exact')
             and (model not in ('reciprocal_full', 'reciprocal_hla') or r['hla_mapping_quality'] == 'exact')]
    answer['exact_mapping'] = comparison_summary(exact, model) if group_count(exact) >= 3 else {'n_rows': len(exact), 'n_groups': group_count(exact), 'unavailable': True}
    answer['single_parent_deletions'] = {g: comparison_summary([r for r in rows if r['parent_type'] != g], model)
                                         for g in sorted(grouped(rows))}
    rng = random.Random(SEED+1)
    comparisons = []
    fields = [('tmb', 'tmb_half_width'), ('orr', 'orr_half_width')]
    if model in ('reciprocal_full', 'reciprocal_clonal'):
        fields.append(('fraction', 'fraction_half_width'))
    if model in ('reciprocal_full', 'reciprocal_hla'):
        fields.append(('loh_fraction', 'loh_half_width'))
    for _ in range(envelope_draws):
        perturbed = []
        for row in rows:
            copy = dict(row)
            for field, width in fields:
                copy[field] += rng.uniform(-1, 1)*row[width]
                if field != 'tmb':
                    copy[field] = clip(copy[field], 1e-12 if field == 'fraction' else 0., 1.)
            perturbed.append(copy)
        comparisons.append(comparison_summary(perturbed, model)['delta_mse'])
    answer['rounding_envelope'] = {'draws': envelope_draws,
                                   'delta_mse_range': [min(comparisons), max(comparisons)]}
    answer['biomarker_selection'] = {
        'known_selected': sum(r['biomarker_selected'] == 'yes' for r in rows),
        'unknown': sum(r['biomarker_selected'] == 'unknown' for r in rows),
        'status': 'unavailable: trial selection is not resolved for the pooled source rows'}
    if any(r.get('source_ci95_low') is not None for r in rows) and model in ('reciprocal_full', 'reciprocal_hla'):
        scenarios = {}
        for label, field in (('all_lower', 'source_ci95_low'), ('all_upper', 'source_ci95_high')):
            perturbed = [dict(r, loh_fraction=r[field] if r.get(field) is not None else r['loh_fraction']) for r in rows]
            scenarios[label] = comparison_summary(perturbed, model)
        answer['hla_source_interval_scenarios'] = scenarios
    return answer


def build_report(directory=HERE / 'data'):
    rows, excluded, tables = load_rows(directory)
    report = analyze_rows(rows)
    report.update(rows=rows, source_rows=tables, excluded_rows=excluded,
                  hashes={p.name: digest(p) for p in sorted(Path(directory).glob('*.csv'))},
                  protocol_sha256=digest(HERE / 'PROTOCOL.md'))
    full = [r for r in rows if r['fraction'] is not None and r['loh_fraction'] is not None]
    exploratory = []
    if 5 <= group_count(full) < 10:
        exploratory.append(panel_analysis(full, 'full_below_gate', 'reciprocal_full'))
    mixed, _, _ = load_rows(directory, mixed_hla=True)
    mixed_full = [r for r in mixed if r['fraction'] is not None and r['loh_fraction'] is not None]
    if group_count(mixed_full) >= 5:
        panel = panel_analysis(mixed_full, 'full_mixed_sources', 'reciprocal_full')
        panel['rows'] = mixed_full
        # This joint diagnostic is exploratory; keep its design fixed on deletion.
        panel['clonal_power'] = fit(mixed_full, 'power_clonal', clonal_with_hla=True)
        panel['clonal_power_jackknife'] = {g: fit([r for r in mixed_full if r['parent_type'] != g],
            'power_clonal', clonal_with_hla=True)['theta'] for g in sorted(grouped(mixed_full))}
        exploratory.append(panel)
    report['exploratory_panels'] = exploratory
    primary = next(p for p in report['panels'] if p['name'] == report['primary_panel'])
    selected = [r for r in rows if r['id'] in primary['row_ids']]
    report['sensitivities'] = {'primary': sensitivity_analysis(selected, primary['model'])}
    if full:
        report['sensitivities']['coherent_full'] = sensitivity_analysis(full, 'reciprocal_full')
        figure_rows, _, _ = load_rows(directory, figure_values=True)
        figure_full = [r for r in figure_rows if r['fraction'] is not None and r['loh_fraction'] is not None]
        report['sensitivities']['figure_instead_of_printed'] = comparison_summary(figure_full, 'reciprocal_full')
    report['source_quality'] = {'clonal_without_indexed_only': sum(
        bool(r['fraction']) and 'indexed' not in r['extraction_method']
        for r in tables['clonal_fractions.csv']),
                              'verdict': 'No clonal/full test remains when the indexed-only source is excluded.'}
    # Check the fitted minima with a denser grid independently of CV selection.
    grid_check = {}
    for panel in report['panels']+exploratory:
        if 'models' not in panel:
            continue
        panel_rows = panel.get('rows') or [r for r in rows if r['id'] in panel['row_ids']]
        for model, output in panel['models'].items():
            if not model.startswith('reciprocal'):
                continue
            dense = fit(panel_rows, model, dense=True)
            difference = output['parameters']['sse']-dense['sse']
            if difference > 1e-9:
                raise ArithmeticError('Coarse profile missed a better minimum')
            grid_check[panel['name']+'/'+model] = difference
    report['grid_check_max_difference'] = max(grid_check.values(), default=0.)
    return report


def _number(value, places=4):
    return 'NA' if value is None else f'{value:.{places}f}'


def _interval(values, multiplier=1., places=4):
    return '['+', '.join(_number(v*multiplier, places) for v in values)+']'


def ranked_outliers(rows, predictions, limit=5):
    groups = grouped(rows)
    ranked = sorted(groups, key=lambda g: (-statistics.mean(
        abs(r['orr']-predictions[r['id']]) for r in groups[g]), g))[:limit]
    return [(g, groups[g]) for g in ranked]


def render_report(report):
    """All displayed empirical numbers are generated from the CSV report."""
    rows = report['rows']
    panels = [p for p in report['panels'] if p['eligible']]+report['exploratory_panels']
    primary = next(p for p in panels if p['name'] == report['primary_panel'])
    lines = ['### Coverage and preregistered gates', '',
             '| Panel | Rows | Parent cancers | Required parents | Status |',
             '| --- | ---: | ---: | ---: | --- |']
    for p in report['panels']:
        lines.append(f'| {p["name"]} | {p["n_rows"]} | {p["n_groups"]} | {p["gate"]} | '+
                     ('eligible' if p['eligible'] else '**unavailable under frozen gate**')+' |')
    for p in report['exploratory_panels']:
        lines.append(f'| {p["name"]} | {p["n_rows"]} | {p["n_groups"]} | unchanged | **exploratory** |')
    source = report['source_rows']
    total = sum(len(v) for v in source.values())
    lines += ['', f'The CSVs retain {total} source rows: '+
              ', '.join(f'{name}: {len(value)}' for name, value in sorted(source.items()))+'.',
              f'TMB is missing for {len(report["excluded_rows"])} original ORR rows. '
              f'The source ORR table reports {sum(int(r["orr_source_n"]) for r in source["tmb_orr.csv"])} '
              'trial participants across its rows; this is a transcribed table total, not a verified unique-person count.',
              '', '### Source inputs and coherent-panel joins', '',
              'ORR and predictor fractions below are percentages; TMB is mutations/Mb. '
              'A dash means missing. HLA is the coherent DASH panel. Source-specific '
              'reading/rounding envelopes and provenance are in the CSVs.', '',
              '| Cancer category | TMB | ORR % | Clonal % | HLA-LOH % | TMB category match |',
              '| --- | ---: | ---: | ---: | ---: | --- |']
    lookup = {r['id']: r for r in rows}
    for raw in source['tmb_orr.csv']:
        r = lookup.get(raw['id'])
        def display(field, multiplier=1, places=2):
            return _number(r[field]*multiplier, places) if r and r.get(field) is not None else '—'
        lines.append(f'| {raw["cancer_type"]} | {display("tmb")} | {100*float(raw["orr"]):.1f} | '
                     f'{display("fraction",100)} | {display("loh_fraction",100)} | {raw["mapping_quality"]} |')
    lines += ['', '### Held-out prediction', '',
              'Each model is refitted after holding out an entire parent cancer. '
              'RMSE and MAE give equal weight to parents. Training R² uses row-level '
              'squared residuals and is descriptive. Models within a panel use identical rows; '
              'comparisons across different panels are not performance improvements.', '',
              '| Panel | Model | Training R² | Held-out RMSE (pp) | Held-out MAE (pp) | Mean residual (pp) |',
              '| --- | --- | ---: | ---: | ---: | ---: |']
    for panel in panels:
        for model, output in panel['models'].items():
            cv = output['cv']
            lines.append(f'| {panel["name"]} | {LABELS[model]} | {_number(output["training_r2"])} | '
                         f'{100*cv["rmse"]:.3f} | {100*cv["mae"]:.3f} | {100*cv["mean_residual"]:.3f} |')
    lines += ['', '### Paired prediction-error differences', '',
              'ΔMSE = designated theory model minus comparator, in squared percentage '
              'points; negative favors theory. Intervals are conditional descriptive '
              f'percentile intervals from {ERROR_REPLICATES} parent resamples, seed {SEED}. '
              'They are not a formal significance test: folds share training data, '
              'cancers are not a random sample, and these are multiple descriptive comparisons.', '',
              '| Panel | Comparator | ΔMSE (pp²) | Descriptive 95% interval (pp²) |',
              '| --- | --- | ---: | --- |']
    for p in panels:
        for model, result in p['comparisons'].items():
            lines.append(f'| {p["name"]} | {LABELS[model]} | {10000*result["delta_mse"]:.3f} | '
                         f'{_interval(result["interval"],10000,3)} |')
    comparison = primary['comparisons']['log_tmb']
    beats_both = all(primary['comparisons'][m]['delta_mse'] < 0 for m in ('log_tmb','linear_tmb'))
    criterion = beats_both and comparison['interval'][1] < 0
    lines += ['', f'**Primary descriptive improvement criterion: {"met" if criterion else "not met"}.** '
              f'The available primary panel is `{report["primary_panel"]}`; '
              'the preregistered full model remains below its complete-case gate.', '',
              '### Fitted composite parameters and uncertainty', '',
              'A₀ and A₁ are fitted endpoint ceilings, B is an uncalibrated composite. '
              'η is a population-mixture penalty, not the presentation-loss hazard. '
              f'Intervals use {PARAMETER_REPLICATES} cancer-group bootstrap draws, seed {SEED}; '
              'they condition on the source values and do not include between-study '
              'measurement error or cure/response calibration.', '',
              '| Panel | Parameter | Estimate | Descriptive 95% interval | Usable draws |',
              '| --- | --- | ---: | --- | ---: |']
    for p in panels:
        param = p['models'][p['model']]['parameters']
        boot = p['parameter_uncertainty']
        keys = ('A0', 'B', 'A1', 'eta') if p['model'] in ('reciprocal_hla','reciprocal_full') else ('A0','B')
        for key in keys:
            lines.append(f'| {p["name"]} | {key} | {_number(param[key])} | '
                         f'{_interval(boot["intervals"][key])} | {boot["usable"][key]} |')
    lines += ['', '| Panel | Skipped singular draws | Fits on any active constraint boundary |',
              '| --- | ---: | ---: |']
    for p in panels:
        b = p['parameter_uncertainty']
        lines.append(f'| {p["name"]} | {b["skipped"]} | {b["boundary_fits"]} |')
    lines += ['', '### HLA sign and exponent diagnostics', '',
              'The HLA diagnostic removes A₁ ≤ A₀ while keeping both ceilings in '
              '[0,1]. Positive A₁ − A₀ (negative η) would reverse the adverse bridge. '
              'The relaxed intervals, rather than an enforced nonnegative sign, show uncertainty.', '',
              '| Panel | A₁ − A₀ | Descriptive 95% interval | η | Descriptive 95% interval | Held-out RMSE (pp) |',
              '| --- | ---: | --- | ---: | --- | ---: |']
    for p in panels:
        if 'sign_diagnostic' not in p:
            continue
        d, b = p['sign_diagnostic'], p['sign_uncertainty']
        v = d['parameters']
        lines.append(f'| {p["name"]} | {_number(v["ceiling_difference"])} | '
                     f'{_interval(b["intervals"]["ceiling_difference"])} | {_number(v["eta"])} | '
                     f'{_interval(b["intervals"]["eta"])} | {100*d["cv"]["rmse"]:.3f} |')
    lines += ['', '| Shape diagnostic | Fitted exponent | Range after each parent deletion | Boundary deletions |',
              '| --- | ---: | --- | ---: |']
    power = report['power_diagnostic']
    values = list(power['jackknife_q'].values())
    lines.append(f'| TMB exponent q; specified q = 1 | {power["parameters"]["q"]:.4f} | '
                 f'{_interval([min(values),max(values)])} | {sum(v<=.050001 or v>=3.999999 for v in values)} |')
    if 'clonal_power_diagnostic' in report:
        power = report['clonal_power_diagnostic']
        values = list(power['jackknife_theta'].values())
        lines.append(f'| Clonal exponent θ without HLA; specified θ = 1 | {power["parameters"]["theta"]:.4f} | '
                     f'{_interval([min(values),max(values)])} | {sum(v<=-1.999999 or v>=2.999999 for v in values)} |')
    for panel in report['exploratory_panels']:
        if 'clonal_power' in panel:
            values = list(panel['clonal_power_jackknife'].values())
            lines.append(f'| Joint θ with mixed-source HLA, exploratory | {panel["clonal_power"]["theta"]:.4f} | '
                         f'{_interval([min(values),max(values)])} | {sum(v<=-1.999999 or v>=2.999999 for v in values)} |')
    lines += ['', 'Deletion ranges are influence diagnostics, **not confidence intervals**. '
              'A fitted boundary is not a precise magnitude estimate. The negative clonal '
              'exponent is descriptive contrary evidence to the positive ecological '
              'coverage proxy, not a causal claim that clonality is harmful.', '',
              '### Every primary held-out prediction', '',
              '| Cancer category | Observed ORR % | Reciprocal % | Log-linear % | Linear % | Reciprocal residual (pp) |',
              '| --- | ---: | ---: | ---: | ---: | ---: |']
    models = primary['models']
    for r in rows:
        if r['id'] not in primary['row_ids']:
            continue
        predictions = [models[m]['cv']['predictions'][r['id']] for m in (primary['model'],'log_tmb','linear_tmb')]
        lines.append(f'| {r["cancer_type"]} | {100*r["orr"]:.3f} | '
                     +' | '.join(f'{100*p:.3f}' for p in predictions)+
                     f' | {100*(r["orr"]-predictions[0]):+.3f} |')
    lines += ['', '### Largest parent-cancer outliers', '',
              'Rank by parent mean absolute held-out residual; preserve each row sign '
              'within a parent. Positive means the observed response exceeds prediction.', '',
              '| Model | Rank | Parent / constituent category | Observed % | Predicted % | Residual (pp) |',
              '| --- | ---: | --- | ---: | ---: | ---: |']
    for model in (primary['model'], 'log_tmb'):
        predicted = models[model]['cv']['predictions']
        for rank, (parent, members) in enumerate(ranked_outliers(rows, predicted), 1):
            for r in members:
                p = predicted[r['id']]
                lines.append(f'| {LABELS[model]} | {rank} | {parent}: {r["cancer_type"]} | '
                             f'{100*r["orr"]:.3f} | {100*p:.3f} | {100*(r["orr"]-p):+.3f} |')
    lines += ['', '### Sensitivities and failed falsification attempts', '',
              'All differences below compare the designated theory with log-linear '
              'TMB on identical rows. Scenario ranges are neither CIs nor exhaustive bounds.', '',
              '| Panel / sensitivity | Parent cancers | Theory RMSE (pp) | Log-linear RMSE (pp) | ΔMSE (pp²) |',
              '| --- | ---: | ---: | ---: | ---: |']
    for name in ('primary','coherent_full'):
        s = report['sensitivities'].get(name)
        if not s:
            continue
        exact = s['exact_mapping']
        if exact.get('unavailable'):
            lines.append(f'| {name}: exact mappings only | {exact["n_groups"]} | unavailable | unavailable | unavailable |')
        else:
            lines.append(f'| {name}: exact mappings only | {exact["n_groups"]} | {100*exact["theory_rmse"]:.3f} | '
                         f'{100*exact["log_rmse"]:.3f} | {10000*exact["delta_mse"]:.3f} |')
        for label, scenario in s.get('hla_source_interval_scenarios', {}).items():
            lines.append(f'| {name}: HLA source intervals {label} | {scenario["n_groups"]} | '
                         f'{100*scenario["theory_rmse"]:.3f} | {100*scenario["log_rmse"]:.3f} | {10000*scenario["delta_mse"]:.3f} |')
    fig = report['sensitivities'].get('figure_instead_of_printed')
    if fig:
        lines.append(f'| coherent_full: figure instead of printed HLA | {fig["n_groups"]} | {100*fig["theory_rmse"]:.3f} | '
                     f'{100*fig["log_rmse"]:.3f} | {10000*fig["delta_mse"]:.3f} |')
    lines += ['', '| Panel | Extraction scenarios | ΔMSE scenario range (pp²) | Parent-deletion ΔMSE range (pp²) | Deletions favoring theory |',
              '| --- | ---: | --- | --- | ---: |']
    for name in ('primary','coherent_full'):
        s = report['sensitivities'].get(name)
        if not s:
            continue
        values = [v['delta_mse'] for v in s['single_parent_deletions'].values()]
        env = s['rounding_envelope']
        lines.append(f'| {name} | {env["draws"]} | {_interval(env["delta_mse_range"],10000,3)} | '
                     f'{_interval([min(values),max(values)],10000,3)} | {sum(v<0 for v in values)}/{len(values)} |')
    lines += ['', f'Excluding indexed-only clonal evidence leaves {report["source_quality"]["clonal_without_indexed_only"]} '
              'usable source clonal categories: clonal/full tests are unavailable. '
              f'Trial selection is unknown for {sum(r["biomarker_selected"] == "unknown" for r in rows)} '
              'included ORR rows, so no unselected-only sensitivity is claimed.', '',
              '### Normalized provenance fingerprints', '',
              '| File | SHA-256 of LF-normalized content |', '| --- | --- |']
    for filename, checksum in report['hashes'].items():
        lines.append(f'| data/{filename} | `{checksum}` |')
    lines.append(f'| PROTOCOL.md | `{report["protocol_sha256"]}` |')
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Recompute and compare the README numeric block')
    parser.add_argument('--write', action='store_true', help='Replace the README numeric block from CSVs')
    parser.add_argument('--json', action='store_true', help='Print every recomputed output as JSON')
    args = parser.parse_args()
    report = build_report()
    block = render_report(report)
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    elif not args.check and not args.write:
        print(block)
    if args.check or args.write:
        path = HERE / 'README.md'
        text = path.read_text(encoding='utf-8')
        before, rest = text.split(BEGIN, 1)
        old, after = rest.split(END, 1)
        expected = '\n\n'+block+'\n\n'
        if args.write:
            path.write_text(before+BEGIN+expected+END+after, encoding='utf-8')
        elif old != expected:
            raise SystemExit('README recomputed results differ: run analysis.py --write')
        else:
            print('Aggregate arithmetic and README results match (nonclinical).')


if __name__ == '__main__':
    main()
