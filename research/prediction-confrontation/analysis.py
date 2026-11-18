#!/usr/bin/env python3
"""Published aggregate confrontation; standard library, offline, nonclinical.

No individual observations, patient reconstruction, source fetching or model
fitting to individual records. Read README.md and PROTOCOL.md before interpreting
outputs. Eligibility is a documented manual source audit, not a statistical test.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
from statistics import NormalDist

HERE = Path(__file__).resolve().parent
FILES = ('immune.csv', 'hla_response.csv', 'ctdna.csv', 'ecdna.csv', 'resistance.csv')
COLUMNS = ('id', 'study', 'population', 'quantity', 'value', 'unit', 'source_url',
           'doi', 'locator', 'access_date', 'reuse_basis', 'notes')
Z95 = NormalDist().inv_cdf(.975)
ALPHA = .05
BEGIN = '<!-- BEGIN RECOMPUTED RESULTS -->'
END = '<!-- END RECOMPUTED RESULTS -->'
UNTESTABLE = 'untestable with current public aggregates'

# These are source-audit annotations, not measurements manufactured from absence.
# Anchors identify the closest retained aggregates. Missing fields cannot be
# replaced by these proxies, and are deliberately NOT assigned numerical zeros.
AUDIT = {
    'I1': dict(anchors=('litchfield_clonal_or', 'ravi_clonal_or'),
               missing=('functional_load_bins', 'founder_count', 'independent_rates',
                        'escape_or_calibrated_response_bridge', 'finite_load_error'),
               reason='Per-SD/log-feature response ORs are not calibrated escape versus N/k.'),
    'I2-H': dict(anchors=('litchfield_subclonal_or', 'litchfield_clonal_or'),
                 missing=('joint_heterogeneity_contrast_at_fixed_burden', 'uncertainty_or_coefficient_covariance'),
                 reason='Separate clonal/subclonal ORs cannot identify a matched heterogeneity effect.'),
    'C1': dict(anchors=('reinert_lead_mean', 'henriksen22_matched_lead_median'),
               missing=('exact_grid_adherence', 'observed_both_parity_zero_imputed_thinned_mean'),
               reason='The CT-date restriction reports a different statistic and design.'),
    'C2': dict(anchors=('henriksen22_fast_mean', 'henriksen24_gain_median'),
               missing=('independent_thresholds_and_residual', 'matched_schedule_capture_counts',
                        'unimodal_log_growth_and_horizon_conditions'),
               reason='Data-defined growth groups and model-extrapolated gains do not test an optimum.'),
    'C3': dict(anchors=('reinert_lead_min', 'reinert_lead_max', 'reinert_draw_interval'),
               missing=('common_detection_window', 'continuous_recurrence_or_actual_imaging_lag_bound'),
               reason='The observed recurrence endpoint is intermittent imaging, not continuous threshold crossing.'),
    'E1': dict(anchors=('lange_2022_pc3_myc_ks_p', 'shibata_2026_clone1_rho'),
               missing=('generation_two_all_cell_moments', 'replicate_uncertainty',
                        'independent_survival_s', 'monomorphic_founders'),
               reason='Sister-fraction fits and calendar-time populations are not the generation-two estimand.'),
    'E2': dict(anchors=('shibata_2026_ctrl1_jq1_p', 'shibata_2026_ctrl1_jq1_rho'),
               missing=('independent_copy_specific_death_division', 'founder_state',
                        'kernel_calibration', 'finite_generation_surviving_well_counts'),
               reason='Drug concentration and fitted segregation parameters do not calibrate gamma/b.'),
    'E3': dict(anchors=('shibata_2026_clone3_ks_p', 'lundberg_2008_la_n5_pct'),
               missing=('conditional_moments_with_uncertainty_or_calibrated_test',
                        'independent_mother_sampling_or_dependence_model'),
               reason='Fit parameters lack uncertainty; published tests have calibration/dependence or classification gaps.'),
    'E4': dict(anchors=('shibata_2026_ctrl1_dmso_rho', 'shibata_2026_ctrl1_jq1_rho'),
               missing=('one_copy_founders', 'generation_two_all_cell_moments',
                        'independent_c_and_R', 'replicate_uncertainty'),
               reason='One-copy founder/generation design and independent dose calibration are absent.'),
    'R1': dict(anchors=('iyer_experimental_vmr', 'iyer_lineage_count'),
               missing=('matched_0_H_2H_3H_count_moments', 'simultaneous_uncertainty',
                        'growth_and_readout_controls', 'shared_clone_law'),
               reason='Across-lineage VMR is not a quartet of independent holiday-culture count laws.'),
    'R2': dict(anchors=('iyer_persister_fraction', 'iyer_experimental_vmr'),
               missing=('state_and_founder_conditioned_covariance', 'uncertainty'),
               reason='Variable founder counts and shared initial states can explain unconditioned excess dispersion.'),
}


def load_data(directory=HERE / 'data'):
    rows = {}
    for filename in FILES:
        with (Path(directory) / filename).open(newline='', encoding='utf-8') as handle:
            reader = csv.DictReader(handle)
            if tuple(reader.fieldnames or ()) != COLUMNS:
                raise ValueError(f'Unexpected columns in {filename}')
            for row in reader:
                if any(not row[c] or not row[c].strip() for c in COLUMNS):
                    raise ValueError(f'Missing value/provenance in {filename}')
                key = row['id']
                if key in rows:
                    raise ValueError(f'Duplicate id: {key}')
                if not math.isfinite(float(row['value'])):
                    raise ValueError(f'Nonfinite value: {key}')
                if not row['source_url'].startswith('https://'):
                    raise ValueError(f'Missing primary URL: {key}')
                rows[key] = dict(row, file=filename)
    return rows


def wilson(k, n, z=Z95):
    if n <= 0 or not 0 <= k <= n:
        raise ValueError('Need 0 <= k <= n and n > 0')
    p = k / n
    den = 1 + z*z/n
    mid = (p + z*z/(2*n)) / den
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return max(0., mid-half), min(1., mid+half)


def hypergeom_weights(a, b, c, d):
    if any(x < 0 or int(x) != x for x in (a, b, c, d)):
        raise ValueError('Contingency counts must be nonnegative integers')
    n1, successes, total = a+b, a+c, a+b+c+d
    if not 0 < n1 < total:
        raise ValueError('Both groups must be nonempty')
    support = range(max(0, n1-(total-successes)), min(n1, successes)+1)
    weights = {x: math.comb(successes, x)*math.comb(total-successes, n1-x)
               for x in support}
    return weights, math.comb(total, n1)


def fisher_two_sided(a, b, c, d):
    weights, den = hypergeom_weights(a, b, c, d)
    return float(Fraction(sum(w for w in weights.values() if w <= weights[a]), den))


def conditional_or_interval(a, b, c, d, alpha=ALPHA):
    """Central exact interval at fixed margins; not a Wald interval.

    Lower solves P_theta(X>=a)=alpha/2; upper solves P_theta(X<=a)=alpha/2.
    The reported point estimate elsewhere is the crude cross-product OR.
    """
    if not 0 < alpha < 1:
        raise ValueError('Need 0 < alpha < 1')
    weights, _ = hypergeom_weights(a, b, c, d)
    xs = list(weights)
    logw = [math.log(weights[x]) for x in xs]

    def tail(logtheta, upper):
        logits = [w+x*logtheta for w, x in zip(logw, xs)]
        mx = max(logits)
        ws = [math.exp(v-mx) for v in logits]
        return sum(w for x, w in zip(xs, ws) if (x >= a if upper else x <= a))/sum(ws)

    def solve(upper):
        lo, hi = -80., 80.
        for _ in range(110):
            mid = (lo+hi)/2
            value = tail(mid, upper)
            if (value < alpha/2) == upper:
                lo = mid
            else:
                hi = mid
        return math.exp((lo+hi)/2)

    low = 0. if a == min(xs) else solve(True)
    high = math.inf if a == max(xs) else solve(False)
    return low, high


def contrast(a, b, c, d):
    if min(a, b, c, d) <= 0:
        raise ValueError('This risk-ratio calculation requires positive cells; no pseudocounts')
    p1, p0 = a/(a+b), c/(c+d)
    ci1, ci0 = wilson(a, a+b), wilson(c, c+d)
    rd = p1-p0
    rd_ci = (max(-1., rd-math.hypot(p1-ci1[0], ci0[1]-p0)),
             min(1., rd+math.hypot(ci1[1]-p1, p0-ci0[0])))
    rr = p1/p0
    rrse = math.sqrt(1/a-1/(a+b)+1/c-1/(c+d))
    return dict(counts=[a,b,c,d], or_=a*d/(b*c), or_ci=conditional_or_interval(a,b,c,d),
                p=fisher_two_sided(a,b,c,d), risk_loh=p1, risk_intact=p0,
                risk_loh_ci=ci1, risk_intact_ci=ci0, rd=rd, rd_ci=rd_ci,
                rr=rr, rr_ci=(rr*math.exp(-Z95*rrse), rr*math.exp(Z95*rrse)))


def holm(pvalues):
    if any(not 0 <= p <= 1 for p in pvalues.values()):
        raise ValueError('Invalid probability')
    order = sorted(pvalues, key=lambda key: (pvalues[key], key))
    out, previous = {}, 0.
    for rank, key in enumerate(order):
        previous = max(previous, min(1., (len(order)-rank)*pvalues[key]))
        out[key] = previous
    return out


def log_or_p(effect, low, high):
    if not 0 < low < effect < high:
        raise ValueError('Need positive, ordered source CI')
    se = (math.log(high)-math.log(low))/(2*Z95)
    return dict(se=se, z=math.log(effect)/se,
                p=math.erfc(abs(math.log(effect)/se)/math.sqrt(2)))


def thinning_bounds(mean, lower, upper, cadence):
    """Sharp plug-in mean/support bounds, conditional on the theorem's design."""
    if not 0 <= lower < upper or not lower <= mean <= upper or cadence <= 0:
        raise ValueError('Invalid support/mean/cadence')
    trunc_lower = ((upper-mean)*min(lower,cadence)
                   +(mean-lower)*min(upper,cadence))/(upper-lower)
    trunc_upper = min(mean, cadence)
    return mean-trunc_upper/2, mean-trunc_lower/2


def rounded_thinning_bounds(mean, lower, upper, cadence, half_unit=.05):
    # For this data rectangle a<d<m<b everywhere. Lower bound m-d/2
    # is increasing in m. Upper bound is increasing in m,b, decreasing
    # in a; corner evaluation therefore gives the global outer envelope.
    if not 0 <= lower-half_unit < lower+half_unit < cadence < mean-half_unit:
        raise ValueError('Unsupported rounding rectangle')
    if not mean+half_unit < upper-half_unit:
        raise ValueError('Overlapping mean/support rounding intervals')
    vals = [thinning_bounds(m,a,b,cadence) for m,a,b in itertools.product(
        (mean-half_unit, mean+half_unit),
        (lower-half_unit, lower+half_unit),
        (upper-half_unit, upper+half_unit))]
    return min(v[0] for v in vals), max(v[1] for v in vals)


def rounded_counts(n, percent, half_unit=Fraction(1,2)):
    if n <= 0 or int(n) != n:
        raise ValueError('An explicit integer denominator is required')
    p = Fraction(str(percent))
    return [k for k in range(n+1) if p-half_unit <= Fraction(100*k,n) <= p+half_unit]


def binomial_upper(k, n, p):
    if not 0 <= k <= n:
        raise ValueError('Invalid binomial count')
    p = Fraction(str(p))
    if not 0 <= p <= 1:
        raise ValueError('Invalid binomial probability')
    return float(sum((Fraction(math.comb(n,j))*p**j*(1-p)**(n-j)
                      for j in range(k,n+1)), Fraction()))


def effective_R(p, rho):
    if not 0 <= p <= 1 or not 0 <= rho <= 1:
        raise ValueError('Allocation parameters outside the declared family')
    return (2*p-1)**2 + 4*p*(1-p)*rho


def eligibility(rows):
    out = {}
    for key, item in AUDIT.items():
        if any(anchor not in rows for anchor in item['anchors']):
            raise ValueError(f'Missing evidence anchor for {key}')
        out[key] = dict(item, verdict=UNTESTABLE)
    return out


def calculate(rows=None):
    rows = load_data() if rows is None else rows
    def v(key): return float(rows[key]['value'])
    def integer(key):
        x = v(key)
        if x != int(x): raise ValueError(f'Non-integer count: {key}')
        return int(x)

    effects = []
    specs = (
        ('litchfield_clonal','Clonal TMB; Litchfield','or','I2',True),
        ('litchfield_total','Total TMB; Litchfield','or','context',False),
        ('litchfield_subclonal','Subclonal TMB; Litchfield','or','context',False),
        ('litchfield_hla','HLA LOH; Litchfield','or','I3',True),
        ('ravi_clonal','Clonal neoantigens; Ravi','or','I2',True),
        ('lozac_luad','HLA-A LOH; Lozac’hmeur LUAD','hr','I3',True),
        ('lozac_lusc','HLA-A LOH; Lozac’hmeur LUSC','hr','I3',True),
    )
    family = {}
    for prefix,label,kind,prediction,formal in specs:
        entry = dict(id=prefix,label=label,kind=kind,prediction=prediction,
                     effect=v(prefix+'_'+kind),ci=(v(prefix+'_low'),v(prefix+'_high')))
        if prefix == 'ravi_clonal':
            entry.update(log_or_p(entry['effect'],*entry['ci']))
            entry.update(p_kind='approximate from rounded CI', source_q=v(prefix+'_q'))
        else:
            entry.update(p=v(prefix+'_p'),p_kind='source reported')
        if prefix.startswith('lozac'): entry['n'] = integer(prefix+'_n')
        if formal: family[prefix] = entry['p']
        effects.append(entry)

    yang = []
    for cohort,label in (('urothelial','Urothelial / anti-PDL1'),
                         ('melanoma_pd1','Melanoma / anti-PD1'),
                         ('melanoma_ctla4','Melanoma / anti-CTLA4')):
        for locus in ('a','b','c'):
            key = f'yang_2022_{cohort}_hla_{locus}'
            counts = [integer(key+'_'+suffix) for suffix in
                      ('response_loh','nonresponse_loh','response_intact','nonresponse_intact')]
            record = dict(id=key,label=label,locus=locus.upper(),**contrast(*counts))
            family[key] = record['p']
            yang.append(record)
    adjusted = holm(family)
    # A printed BH-adjusted p bounds the raw p only before rounding. This
    # sensitivity assumes rounding to the nearest last displayed decimal.
    q_decimal = Decimal(rows['ravi_clonal_q']['value'])
    q_upper = float(q_decimal + Decimal(5).scaleb(q_decimal.as_tuple().exponent-1))
    conservative = dict(family, ravi_clonal=q_upper)
    adjusted_bound = holm(conservative)
    for record in effects + yang:
        key = record['id']
        if key in adjusted:
            record['holm'] = adjusted[key]
            record['holm_q_bound'] = adjusted_bound[key]

    m,a,b,d = (v('reinert_'+name) for name in
               ('lead_mean','lead_min','lead_max','draw_interval'))
    th = thinning_bounds(m,a,b,d)
    ctdna = dict(mean=m,minimum=a,maximum=b,cadence=d,
                 forecast=th,loss=(m-th[1],m-th[0]),
                 rounding_envelope=rounded_thinning_bounds(m,a,b,d),
                 width=b-a, width_rounding=(b-a-.1,b-a+.1),
                 recurrence_counts=[integer('reinert_detected_recurrences'),integer('reinert_total_recurrences')],
                 fast_log_of_mean=math.log(v('henriksen22_fast_mean')),
                 slow_log_of_mean=math.log(v('henriksen22_slow_mean')),
                 matched_counts=[integer('henriksen22_matched_early_n'),integer('henriksen22_matched_n')],
                 matched_median_iqr=[v('henriksen22_matched_lead_median'),v('henriksen22_matched_lead_q1'),v('henriksen22_matched_lead_q3')],
                 extrapolated_counts=[integer('henriksen24_delay_n'),integer('henriksen24_growth_n')],
                 extrapolated_median_iqr=[v('henriksen24_gain_median'),v('henriksen24_gain_q1'),v('henriksen24_gain_q3')],
                 tracerx_summary=[v('tracerx23_lead_median'),v('tracerx23_lead_min'),v('tracerx23_lead_max'),integer('tracerx23_lead_n')])

    clones, clone_p, clone_p_ambiguous = [], {}, {}
    for i in range(1,13):
        key = f'shibata_2026_clone{i}'
        p,rho = v(key+'_p'),v(key+'_rho')
        scope = 'HSR control' if i in (4,6) else ('classification conflict' if i==7 else 'ecDNA')
        entry = dict(id=key,clone=i,p=p,rho=rho,R=effective_R(p,rho),
                     ks_p=v(key+'_ks_p'),n=integer(key+'_n'),scope=scope)
        clones.append(entry)
        if i not in (4,6,7): clone_p[key] = entry['ks_p']
        if i not in (4,6): clone_p_ambiguous[key] = entry['ks_p']
    adj9,adj10 = holm(clone_p),holm(clone_p_ambiguous)
    for record in clones:
        record['holm9'] = adj9.get(record['id'])
        record['holm10'] = adj10.get(record['id'])
    arms = []
    for arm in ('ctrl1_dmso','ctrl1_jq1','ctrlen2_dmso','ctrlen2_jq1'):
        key = 'shibata_2026_'+arm
        p,rho=v(key+'_p'),v(key+'_rho')
        arms.append(dict(arm=arm,p=p,rho=rho,R=effective_R(p,rho)))
    lange = []
    for panel in ('pc3_myc','snu16_myc','snu16_fgfr2','gbm39_egfr','colo_myc','tr14_mycn','tr14_cdk4','colo320_dm'):
        key = 'lange_2022_'+panel
        lange.append(dict(panel=panel,D=v(key+'_ks_d'),p=v(key+'_ks_p'),n=integer(key+'_n')))

    lundberg, lp = [], {}
    nominal = v('lundberg_2008_alpha')
    for panel in ('chp212','la_n5','mcix'):
        key = 'lundberg_2008_'+panel
        n,pct = integer(key+'_n'),v(key+'_pct')
        counts = rounded_counts(n,pct)
        if not counts: raise ValueError('No count compatible with published percentage')
        tails = [binomial_upper(k,n,nominal) for k in counts]
        # Preserve the full set; use conservative extrema if rounding is nonunique.
        ci = (wilson(min(counts),n)[0],wilson(max(counts),n)[1])
        k=min(counts)
        item = dict(id=key,panel=panel,n=n,percent=pct,counts=counts,
                    rate=(min(counts)/n,max(counts)/n),ci=ci,p=max(tails),
                    excess=(min(counts)/n-nominal,max(counts)/n-nominal),
                    fold=(min(counts)/n/nominal,max(counts)/n/nominal),
                    markov=min(1.,n*nominal/k) if k else 1.)
        lundberg.append(item)
        lp[key]=item['p']
    lh=holm(lp)
    combined=holm(dict(family,**lp))
    combined_bound=holm(dict(conservative,**lp))
    for record in lundberg:
        key=record['id']
        record.update(holm3=lh[key],holm_combined=combined[key],
                      holm_combined_q_bound=combined_bound[key])

    persister,vmr=v('iyer_persister_fraction'),v('iyer_experimental_vmr')
    fixed=1-persister
    resistance=dict(persister=persister,vmr=vmr,n=integer('iyer_lineage_count'),
                    fixed_size_binomial_vmr=fixed,excess=vmr-fixed,ratio=vmr/fixed)
    verdicts = {key: value['verdict'] for key,value in eligibility(rows).items()}
    # Holm non-rejection does not establish directional consistency for the
    # universal HLA-LOH bridge, especially with opposite signs and missing
    # functional-presentation calibration.
    i2 = [r for r in effects if r['prediction']=='I2']
    i3 = [r for r in effects if r['prediction']=='I3'] + yang
    i2_counter = any(r['effect']<1 and r['holm']<ALPHA for r in i2)
    i3_counter = any((r.get('effect',r.get('or_')) > 1 if r.get('kind','or')=='or'
                      else r['effect'] < 1) and r['holm']<ALPHA for r in i3)
    verdicts['I2'] = 'inconsistent (falsified)' if i2_counter else 'consistent (clonal-burden association only)'
    verdicts['I3'] = ('inconsistent (falsified)' if i3_counter
                      else 'untestable with current public aggregates (conflicting directions; no corrected rejection)')
    verdicts['E3-X'] = UNTESTABLE
    return dict(row_count=len(rows),file_counts={f:sum(r['file']==f for r in rows.values()) for f in FILES},
                effects=effects,yang=yang,formal_family=family,formal_holm=adjusted,
                formal_holm_q_bound=adjusted_bound,ravi_q_rounding_upper=q_upper,
                ctdna=ctdna,clones=clones,arms=arms,lange=lange,
                lundberg=lundberg,resistance=resistance,eligibility=eligibility(rows),verdicts=verdicts)


def fmt(x):
    if x is None: return 'excluded'
    return f'{x:.6g}'


def interval(x): return '['+', '.join(fmt(v) for v in x)+']'


def markdown(result):
    lines = ['### Recomputed results', '',
             'Sampling intervals below are marginal 95% intervals; model/forecast bounds, ranges, IQRs and rounding envelopes are not confidence intervals. '
             'Holm uses all 14 specified immune contrasts; they are not 14 independent studies. '
             'Source regression/KS test inputs are absent from the retained aggregates, so their p-values are transcribed; all subsequent transformations are checked.', '',
             '| CSV | Published scalar rows |', '| --- | ---: |']
    for f,n in result['file_counts'].items(): lines.append(f'| `{f}` | {n} |')
    lines += ['', '#### Immune associations', '',
              'OR > 1 favors response; HR > 1 means greater death hazard. These are different endpoints. '
              'Litchfield ORs are per within-study SD, including the binary HLA feature; Ravi uses a log-transformed feature with no calibrated doubling interpretation. '
              'The q-bound column uses an upper bound on the source BH-adjusted value conditional on nearest-last-digit rounding; it is not a claimed raw p.', '',
              '| Contrast | Effect | Source CI | p used | p provenance | Holm 14 | Holm, Ravi q bound |',
              '| --- | ---: | --- | ---: | --- | ---: | ---: |']
    for r in result['effects']:
        lines.append(f"| {r['label']} | {r['kind'].upper()} {fmt(r['effect'])} | {interval(r['ci'])} | {fmt(r['p'])} | {r['p_kind']} | {fmt(r.get('holm'))} | {fmt(r.get('holm_q_bound'))} |")
    ravi=next(r for r in result['effects'] if r['id']=='ravi_clonal')
    lines += ['',f"Ravi source BH q = {fmt(ravi['source_q'])}; rounding-conditional upper bound = {fmt(result['ravi_q_rounding_upper'])}; reconstructed log-OR SE = {fmt(ravi['se'])}, z = {fmt(ravi['z'])}. Total/subclonal TMB are contextual comparators, excluded from formal inference; no difference between their correlated ORs is tested.", '',
              '| Yang cohort / HLA locus | LOH R, NR | Intact R, NR | Crude response OR | Exact conditional OR CI | Fisher two-sided p | Holm 14 |',
              '| --- | --- | --- | ---: | --- | ---: | ---: |']
    for r in result['yang']:
        a,b,c,d=r['counts']
        lines.append(f"| {r['label']} / {r['locus']} | {a}, {b} | {c}, {d} | {fmt(r['or_'])} | {interval(r['or_ci'])} | {fmt(r['p'])} | {fmt(r['holm'])} |")
    lines += ['', '| Yang cohort / locus | LOH response risk (Wilson CI) | Intact response risk (Wilson CI) | Risk difference (Newcombe CI) | Risk ratio (log-Wald CI) |',
              '| --- | --- | --- | --- | --- |']
    for r in result['yang']:
        lines.append(f"| {r['label']} / {r['locus']} | {fmt(r['risk_loh'])} {interval(r['risk_loh_ci'])} | {fmt(r['risk_intact'])} {interval(r['risk_intact_ci'])} | {fmt(r['rd'])} {interval(r['rd_ci'])} | {fmt(r['rr'])} {interval(r['rr_ci'])} |")
    c=result['ctdna']
    lines += ['', '#### ctDNA: forecasts and endpoint audit', '',
              '| Quantity | Recomputed value | Interpretation |','| --- | --- | --- |',
              f"| Reinert original detected-only mean / range / planned cadence | {fmt(c['mean'])} / {interval((c['minimum'],c['maximum']))} / {fmt(c['cadence'])} months | Range is not a CI; exact grid not documented |",
              f"| Both-parity thinning mean forecast | {interval(c['forecast'])} months | Conditional forecast, no observed target aggregate |",
              f"| Forecast loss of lead | {interval(c['loss'])} months | Sharp plug-in bounds, not a CI |",
              f"| Outward decimal-rounding forecast envelope | {interval(c['rounding_envelope'])} months | Assumes nearest tenth, not sampling uncertainty |",
              f"| Original lead range width | {fmt(c['width'])} months; rounding {interval(c['width_rounding'])} | Cannot reject C3 without its observation assumptions |",
              f"| Reinert detected / recurrent | {c['recurrence_counts'][0]} / {c['recurrence_counts'][1]} | Descriptive; no cross-cohort rate comparison |",
              f"| Henriksen 2022 CT-date restriction | {c['matched_counts'][0]} / {c['matched_counts'][1]} earlier; among earlier detections median {fmt(c['matched_median_iqr'][0])}, IQR {interval(c['matched_median_iqr'][1:])} months | Not alternate-grid, both-parity, zero-imputed mean |",
              f"| Log of Henriksen 2022 group mean multiplier | slow {fmt(c['slow_log_of_mean'])}, fast {fmt(c['fast_log_of_mean'])} per month | Not mean log-growth; data-defined groups, no bimodality test |",
              f"| Henriksen 2024 frequency extrapolation | {c['extrapolated_counts'][0]} / {c['extrapolated_counts'][1]} delayed; delayed subgroup median gain {fmt(c['extrapolated_median_iqr'][0])}, IQR {interval(c['extrapolated_median_iqr'][1:])} months | Model-derived, not an independently observed monthly schedule |",
              f"| TRACERx median / range / recurrent n | {fmt(c['tracerx_summary'][0])} / {interval(c['tracerx_summary'][1:3])} days / {c['tracerx_summary'][3]} | Zeros included for absent/post-relapse detection; median not mean |", '',
              '#### ecDNA: fit transformations and excluded sensitivity tests', '',
              'Shibata R values below are transformations of fitted parameters with no uncertainty; n is printed in the source, with fractions plotted per daughter, not asserted independent mothers. '
              'KS adjustments are sensitivity arithmetic only. Including the disputed clone cannot resolve its classification or test calibration.', '',
              '| Clone | Source p allocation | Source rho | Derived R | Source KS p | n | Scope | Holm 9 | Holm 10, ambiguous 7 included |',
              '| --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: |']
    for r in result['clones']:
        lines.append(f"| {r['clone']} | {fmt(r['p'])} | {fmt(r['rho'])} | {fmt(r['R'])} | {fmt(r['ks_p'])} | {r['n']} | {r['scope']} | {fmt(r['holm9'])} | {fmt(r['holm10'])} |")
    er=[r['R'] for r in result['clones'] if r['scope']=='ecDNA']
    lines += ['',f"The eligible-by-classification fitted R range is {interval((min(er),max(er)))}; this is not a confidence interval.", '',
              '| Shibata arm ID | p allocation | rho | Derived R |','| --- | ---: | ---: | ---: |']
    for r in result['arms']: lines.append(f"| {r['arm']} | {fmt(r['p'])} | {fmt(r['rho'])} | {fmt(r['R'])} |")
    lines += ['', 'Lange supplementary Table SI1 is retained only for the source-version audit; conflicting main/SI representations and missing conditional moments prevent formal inference.', '',
              '| SI1 panel | Source KS D | Source KS p | Source n |','| --- | ---: | ---: | ---: |']
    for r in result['lange']: lines.append(f"| {r['panel']} | {fmt(r['D'])} | {fmt(r['p'])} | {r['n']} |")
    lines += ['', '#### E3-X: exploratory after extraction, conditional on valid tests and independent mitoses', '',
              'The percentage-to-count conversion uses the published denominator and every compatible integer numerator. '
              'The Wilson intervals and binomial tail require independent sampling. Markov removes independence but still requires marginal test size at most the source cutoff. '
              'Excess/fold refer to the rejection-rate ceiling, not a biological allocation effect.', '',
              '| Cell line | Published % / n | Compatible k | Rate, Wilson CI | Excess / fold vs nominal | Conditional p | Holm 3 | Conditional Holm combined 17 | Markov upper bound |',
              '| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: |']
    for r in result['lundberg']:
        lines.append(f"| {r['panel']} | {fmt(r['percent'])}% / {r['n']} | {', '.join(map(str,r['counts']))} | {fmt(r['rate'][0])} {interval(r['ci'])} | {fmt(r['excess'][0])} / {fmt(r['fold'][0])} | {fmt(r['p'])} | {fmt(r['holm3'])} | {fmt(r['holm_combined'])} | {fmt(r['markov'])} |")
    r=result['resistance']
    lines += ['', '#### Resistance: contextual dispersion, no holiday-law test', '',
              f"Iyer reports approximately p = {fmt(r['persister'])}, VMR = {fmt(r['vmr'])}, with {r['n']} retained lineages. "
              f"The fixed-founder independent-binomial comparator is 1-p = {fmt(r['fixed_size_binomial_vmr'])}; "
              f"excess = {fmt(r['excess'])}, ratio = {fmt(r['ratio'])}. "
              'No CI or p-value is calculable for the required conditioned contrast from these aggregates. Variable founder counts alone change this comparator.', '',
              '#### Verdicts under the frozen rules', '', '| ID | Verdict |', '| --- | --- |']
    for key in ('I1','I2','I2-H','I3','C1','C2','C3','E1','E2','E3','E4','R1','R2','E3-X'):
        lines.append(f"| {key} | **{result['verdicts'][key]}** |")
    lines += ['', 'I3 remains unresolved: nominal opposite HLA-LOH contrasts do not survive Holm correction, '
              'and the public aggregates cannot establish the universal adverse bridge. '
              'E3-X gives a conditional challenge, not an unconditional empirical falsification.']
    return '\n'.join(lines)+'\n'


def check_readme(result=None):
    text=(HERE/'README.md').read_text(encoding='utf-8')
    if text.count(BEGIN)!=1 or text.count(END)!=1:
        raise ValueError('Missing/duplicate generated-result markers')
    actual=text.split(BEGIN,1)[1].split(END,1)[0].strip()
    expected=markdown(calculate() if result is None else result).strip()
    if actual!=expected:
        raise ValueError('README result block is stale; regenerate from the CSVs')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json',action='store_true',help='Print all computed results as JSON')
    parser.add_argument('--check',action='store_true',help='Fail if README result tables differ')
    args=parser.parse_args()
    result=calculate()
    if args.check:
        check_readme(result)
        print('Aggregate arithmetic and README synchronization passed; not biological validation.')
    elif args.json:
        print(json.dumps(result,indent=2,allow_nan=False))
    else:
        print(markdown(result),end='')


if __name__=='__main__':
    main()
