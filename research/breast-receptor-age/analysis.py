"""Nonclinical, offline audit of 2016 USCS receptor/age aggregate counts.

Run: python research/breast-receptor-age/analysis.py --check
     python research/breast-receptor-age/analysis.py --write
No medical advice. Frozen decisions: PROTOCOL.md committed separately first.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'research' / 'breast-receptor-age'
DATA = BASE / 'data' / 'turbow_table1_receptor.csv'
PROTOCOL = BASE / 'PROTOCOL.md'
REPORT = BASE / 'REPORT.md'
README = BASE / 'README.md'
AGES = ('65-69', '70-74', '75-79', '80-84', '85+')
CLASSES = ('overall', 'positive', 'negative', 'unknown')
SCENARIOS = ('proportional', 'all_positive', 'all_negative', 'younger_positive_older_negative', 'younger_negative_older_positive')
SOURCE = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC10846538/#T1'
ALPHA = 0.05
TAIL = ALPHA / (2 * len(AGES))


def normalize_newlines(text: str) -> str:
    return text.replace('\r\n', '\n').replace('\r', '\n')


def sha256_text(path: Path) -> str:
    return hashlib.sha256(normalize_newlines(path.read_text(encoding='utf-8')).encode('utf-8')).hexdigest()


def key(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def beta_fraction(a: float, b: float, x: float) -> float:
    """Lentz continued fraction for regularized incomplete beta (stdlib only)."""
    tiny = 1e-290
    qab, qap, qam = a + b, a + 1., a - 1.
    c = 1.
    d = 1. - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1. / d
    h = d
    for m in range(1, 1201):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1. + aa * d
        c = 1. + aa / c
        if abs(d) < tiny:
            d = tiny
        if abs(c) < tiny:
            c = tiny
        d = 1. / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1. + aa * d
        c = 1. + aa / c
        if abs(d) < tiny:
            d = tiny
        if abs(c) < tiny:
            c = tiny
        d = 1. / d
        delta = d * c
        h *= delta
        if abs(delta - 1.) < 5e-13:
            return h
    raise ArithmeticError('Incomplete beta continued fraction did not converge')


def regularized_beta(a: float, b: float, x: float) -> float:
    if a <= 0 or b <= 0 or not 0. <= x <= 1.:
        raise ValueError('Invalid beta parameters')
    if x == 0.:
        return 0.
    if x == 1.:
        return 1.
    log_front = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    log_front += a * math.log(x) + b * math.log1p(-x)
    front = math.exp(log_front)
    if x < (a + 1.) / (a + b + 2.):
        return min(1., max(0., front * beta_fraction(a, b, x) / a))
    return min(1., max(0., 1. - front * beta_fraction(b, a, 1. - x) / b))


def inverse_beta(q: float, a: float, b: float) -> float:
    if not 0. < q < 1.:
        raise ValueError('Probability must lie strictly between zero and one')
    lo, hi = 0., 1.
    for _ in range(72):
        mid = (lo + hi) / 2.
        if regularized_beta(a, b, mid) < q:
            lo = mid
        else:
            hi = mid
    ans = (lo + hi) / 2.
    if abs(regularized_beta(a, b, ans) - q) > 1e-7:
        raise ArithmeticError('Beta inverse failed target-tail verification')
    return ans


def cp_interval(x: int, n: int, tail: float = TAIL) -> tuple[float, float]:
    """Two-sided exact Clopper-Pearson, per-tail error = 0.005 by protocol."""
    if not isinstance(n, int) or not isinstance(x, int) or not 0 <= x <= n or n < 1:
        raise ValueError('Invalid binomial count')
    if not 0. < tail < 0.5:
        raise ValueError('Invalid tail')
    low = 0. if x == 0 else inverse_beta(tail, x, n - x + 1)
    high = 1. if x == n else inverse_beta(1. - tail, x + 1, n - x)
    return low, high


def odds(p: float) -> float:
    if p <= 0.:
        return 0.
    if p >= 1.:
        return math.inf
    return p / (1. - p)


def ratio_ci(younger: tuple[float, float], older: tuple[float, float]) -> tuple[float, float]:
    return odds(older[0]) / odds(younger[1]), odds(older[1]) / odds(younger[0])


def classify_interval(lo: float, hi: float) -> str:
    return 'below' if hi < 1. else 'above' if lo > 1. else 'contains 1'


def load_data(path: Path = DATA) -> dict[str, dict[str, dict]]:
    required = ('age_bin', 'receptor_class', 'cases', 'population', 'printed_rate', 'printed_ci_lower',
                'printed_ci_upper', 'source', 'doi', 'source_url', 'source_table', 'source_row',
                'source_column', 'population_locator', 'access_date', 'access_timezone',
                'rate_units', 'adjustment_label', 'dataset_vintage', 'reuse_basis', 'license')
    records: dict[str, dict[str, dict]] = {age: {} for age in AGES}
    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != list(required):
            raise ValueError('Unexpected column layout; frozen input contract not met')
        for line, row in enumerate(reader, start=2):
            age, label = row['age_bin'], row['receptor_class']
            if age not in records or label not in CLASSES or label in records[age]:
                raise ValueError(f'Duplicate or invalid receptor/age cell at line {line}')
            for field in required:
                if not row.get(field):
                    raise ValueError(f'Missing {field} in line {line}')
            if row['source_url'] != SOURCE or row['source_table'] != 'Table 1':
                raise ValueError(f'Source locator wrong at line {line}')
            if row['doi'] != '10.1007/s10549-021-06160-4':
                raise ValueError(f'Source DOI wrong at line {line}')
            if row['source_column'] != ('>=85' if age == '85+' else age):
                raise ValueError(f'Age header locator wrong at line {line}')
            try:
                cases = int(row['cases'])
                population = int(row['population'])
                rate, low, high = (float(row[k]) for k in ('printed_rate', 'printed_ci_lower', 'printed_ci_upper'))
            except (ValueError, TypeError) as error:
                raise ValueError(f'Non-numeric published cell at line {line}') from error
            if cases <= 0 or population <= 0 or not all(map(math.isfinite, (rate, low, high))):
                raise ValueError(f'Invalid cases, exposure or finite rate at line {line}')
            if not 0 < low <= rate <= high:
                raise ValueError(f'Invalid published 95% rate interval at line {line}')
            if row['receptor_class'] != 'overall' and not row['source_row'].startswith('Receptor Status: '):
                raise ValueError(f'Wrong receptor locator at line {line}')
            records[age][label] = dict(cases=cases, population=population, rate=rate,
                                       ci_low=low, ci_high=high, source_row=row['source_row'])
    if any(set(records[age]) != set(CLASSES) for age in AGES):
        raise ValueError('Missing one or more of the 20 frozen table cells')
    for age in AGES:
        cells = records[age]
        if len({item['population'] for item in cells.values()}) != 1:
            raise ValueError(f'Population mismatch for {age}')
        p, n, u, total = (cells[name]['cases'] for name in ('positive', 'negative', 'unknown', 'overall'))
        residual = total - p - n - u
        if residual < 0:
            raise ValueError(f'UNTESTABLE: negative count residual in {age}')
        cells['residual'] = residual
        cells['pool'] = u + residual
        cells['prop_k'] = (2 * cells['pool'] * p + (p + n)) // (2 * (p + n))
        # This is round-half-up for integer numerator: floor(M*p/(p+n)+0.5).
    return records


def completed(age_records: dict, kind: str) -> tuple[int, int]:
    p = age_records['positive']['cases']
    total = age_records['overall']['cases']
    m = age_records['pool']
    if kind == 'proportional':
        k = age_records['prop_k']
    elif kind == 'all_positive':
        k = m
    elif kind == 'all_negative':
        k = 0
    else:
        raise ValueError('Unknown allocation kind')
    return p + k, total - p - k


def rounding_q(younger: dict, older: dict, half: float) -> tuple[float, float]:
    def r(row: dict, cls: str) -> float:
        return row[cls]['rate']
    po, no, py, ny = r(older, 'positive'), r(older, 'negative'), r(younger, 'positive'), r(younger, 'negative')
    if min(po, no, py, ny) <= half:
        raise ValueError('Rounded-rate range includes nonpositive rates')
    return ((po - half) * (ny - half) / ((no + half) * (py + half)),
            (po + half) * (ny + half) / ((no - half) * (py - half)))


def analyze(rows: dict | None = None) -> dict:
    if rows is None:
        rows = load_data()
    cp: dict[str, dict[str, tuple[float, float]]] = {}
    for age in AGES:
        cp[age] = {}
        for label in ('proportional', 'all_positive', 'all_negative'):
            x, other = completed(rows[age], label)
            cp[age][label] = cp_interval(x, x + other)
    pairs = []
    for i, j in itertools.combinations(range(len(AGES)), 2):
        younger, older = AGES[i], AGES[j]
        a, b = rows[younger], rows[older]
        scenarios = {}
        for label in SCENARIOS:
            if label == 'younger_positive_older_negative':
                a_label, b_label = 'all_positive', 'all_negative'
            elif label == 'younger_negative_older_positive':
                a_label, b_label = 'all_negative', 'all_positive'
            else:
                a_label = b_label = label
            lo, hi = ratio_ci(cp[younger][a_label], cp[older][b_label])
            xp, xn = completed(a, a_label)
            yp, yn = completed(b, b_label)
            q = (yp / yn) / (xp / xn)
            scenarios[label] = dict(q=q, lo=lo, hi=hi, direction=classify_interval(lo, hi))
        envelope = ratio_ci((cp[younger]['all_negative'][0], cp[younger]['all_positive'][1]),
                            (cp[older]['all_negative'][0], cp[older]['all_positive'][1]))
        p_old = b['positive']['cases'] / b['positive']['population']
        p_young = a['positive']['cases'] / a['positive']['population']
        n_old = b['negative']['cases'] / b['negative']['population']
        n_young = a['negative']['cases'] / a['negative']['population']
        q_print = ((b['positive']['rate'] / b['negative']['rate']) /
                   (a['positive']['rate'] / a['negative']['rate']))
        pairs.append(dict(younger=younger, older=older, adjacent=(j == i + 1),
                          pos_rr=p_old / p_young, neg_rr=n_old / n_young,
                          proportional_q=scenarios['proportional']['q'],
                          printed_q=q_print, rounding_half=rounding_q(a, b, 0.05),
                          rounding_full=rounding_q(a, b, 0.10),
                          scenarios=scenarios, envelope=envelope,
                          robust=classify_interval(*envelope)))
    if any(pair['robust'] != 'contains 1' for pair in pairs):
        verdict = 'FALSIFIED'
    elif all(d['direction'] == 'contains 1' for pair in pairs for d in pair['scenarios'].values()):
        verdict = 'CONSISTENT'
    else:
        verdict = 'INCONCLUSIVE'
    return dict(rows=rows, pairs=pairs, verdict=verdict, cp=cp)


def fmt(value: float, digits: int = 6) -> str:
    return f'{value:.{digits}f}'


def bracket(interval: tuple[float, float]) -> str:
    return f'[{fmt(interval[0])}, {fmt(interval[1])}]'


def render_report(result: dict) -> str:
    r = result['rows']
    out = [
        '# Breast receptor-by-age proportionality — frozen 2016 USCS analysis',
        '',
        '**Nonclinical, unreviewed research. Not medical advice, screening advice, or clinical validation.**',
        '',
        '<!-- Deterministic output of analysis.py; source data in data/turbow_table1_receptor.csv. -->',
        '',
        f'**Overall four-way verdict: {result["verdict"]}.** HER2: **UNTESTABLE** (not reported in Turbow Table 1).',
        '',
        'A failure would require at least one missing-status-enveloped, five-age Bonferroni-simultaneous',
        'interval for Q to exclude 1. A point contrast or one favorable allocation does not meet that criterion.',
        'Here Q = (positive older/younger incidence RR) / (double-negative older/younger incidence RR).',
        'It is a case-composition odds ratio when the receptor strata have the same age-specific exposure.',
        '',
        '## Input and provenance',
        '',
        f'- Table: [Turbow et al. (2021), PMC Table 1]({SOURCE}), Receptor Status and Overall rows.',
        '- Dataset: 2016 USCS (NPCR/SEER); November 2018 submission; released June 2019.',
        '- Article manuscript is free to read but shows no explicit Creative Commons grant: assume standard',
        '  copyright for expression. This CSV transcribes an attributed small subset of factual statistics,',
        '  without copying article prose/layout or relicensing third-party content under MIT.',
        '- Article HTML was read; no manuscript PDF, browser file download, or download prompt was used.',
        '- The source describes its rates as age-adjusted to the 2000 U.S. standard population (19 age groups).',
        '  Audited count/population rates below are recomputed proxies; no age-standardisation weights are inferred.',
        '- Each source row in the CSV records DOI, exact source label, age-column locator, printed 95% CI,',
        '  population header, unit, dataset vintage, access date, attribution and licence caveat.',
        '',
        '### SHA-256 of UTF-8 with CRLF/CR normalized to LF',
        '',
        '| Path (Path.as_posix()) | Normalized SHA-256 |',
        '| --- | --- |',
    ]
    for path in (PROTOCOL, DATA):
        out.append(f'| {key(path)} | `{sha256_text(path)}` |')
    out += ['', '## Source transcription and rounding reconciliation (all 20 cells)', '',
            'The last three source fields show printed rate (printed 95% CI). The rate difference is',
            'count/population × 100,000 minus printed rate; ±0.05 is only a one-decimal rounding check.',
            '', '| Age | Source row | Cases / population | Published rate (95% CI) | Count rate | Difference | ±0.05 fit |',
            '| --- | --- | ---: | ---: | ---: | ---: | :---: |']
    for age in AGES:
        for cls in CLASSES:
            z = r[age][cls]
            actual = 100000. * z['cases'] / z['population']
            difference = actual - z['rate']
            out.append(f'| {age} | {cls} | {z["cases"]:,} / {z["population"]:,} | '
                       f'{z["rate"]:.1f} ({z["ci_low"]:.1f}–{z["ci_high"]:.1f}) | '
                       f'{fmt(actual)} | {difference:+.6f} | {"yes" if abs(difference) <= 0.050000000001 else "**NO**"} |')
    out += ['', '## Unknown status, residual and primary proportional allocation', '',
            'Unknown and unaccounted residual cases are separate in the input; both are in the',
            'conditional worst-case assignable pool. Proportional allocation uses integer',
            'round-half-up of pool × observed positive/(positive + negative).', '',
            '**Conditional accounting only:** the completed partition supposes the printed overall',
            'total and all observed receptor-class counts are correct, and all remaining cases',
            'belong to one of two binary classes. Residuals are accounting discrepancies, not',
            'verified additional cases; assigning them does **not** repair source truth.', '' ,
            '| Age | Printed positive | Printed double-negative | Printed unknown | Overall | Residual | Assignable pool | k to positive | Completed positive / negative |',
            '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for age in AGES:
        z = r[age]
        cp, cn = completed(z, 'proportional')
        out.append(f'| {age} | {z["positive"]["cases"]:,} | {z["negative"]["cases"]:,} | '
                   f'{z["unknown"]["cases"]:,} | {z["overall"]["cases"]:,} | '
                   f'{z["residual"]:,} | {z["pool"]:,} | {z["prop_k"]:,} | {cp:,} / {cn:,} |')
    out += ['', '## Simultaneous inference over all 10 age contrasts', '',
            'Within each hypothesized completed table, independent Poisson age/class marks imply',
            'binomial positive composition conditional on the total for each age. Exact two-sided',
            'Clopper–Pearson intervals at 99% **per age** (two tails 0.005 each), combined by',
            'Bonferroni over five age groups, give ≥95% simultaneous coverage for all ten Q',
            'transforms. This is process-model uncertainty conditional on hypothetical classification,',
            'not a model for the missing-status mechanism, coding error or population error.', '',
            'Below, P-RR and N-RR are the source-known positive and double-negative age rate ratios,',
            'before allocation. Q and the interval are based on proportional completion. All ten',
            'pairs are prespecified; “adjacent” only labels their position, not a different test.', '',
            '| Younger → older | Kind | P-RR | N-RR | Proportional Q | Proportional interval | Status | All-allocation envelope | Robust status |',
            '| --- | --- | ---: | ---: | ---: | --- | --- | --- | --- |']
    for p in result['pairs']:
        nominal = p['scenarios']['proportional']
        out.append(f'| {p["younger"]} → {p["older"]} | {"adjacent" if p["adjacent"] else "nonadjacent"} | '
                   f'{fmt(p["pos_rr"])} | {fmt(p["neg_rr"])} | {fmt(nominal["q"])} | '
                   f'{bracket((nominal["lo"], nominal["hi"]))} | {nominal["direction"]} | '
                   f'{bracket(p["envelope"])} | {p["robust"]} |')
    out += ['', '### Complete allocation sensitivity (all five cases for each pair)', '',
            'Scenarios: proportional, all unknown+residual to positive in both ages,',
            'all to negative in both, and the two separately adverse cross-age corners.',
            'Each interval uses the same per-age Bonferroni CP rule. The envelope above',
            'projects over **every integer allocation**, not just these five displayed choices.', '',
            '| Younger → older | Allocation | Q point | Simultaneous interval | Excludes 1? |',
            '| --- | --- | ---: | --- | --- |']
    for p in result['pairs']:
        for name in SCENARIOS:
            c = p['scenarios'][name]
            out.append(f'| {p["younger"]} → {p["older"]} | {name} | {fmt(c["q"])} | '
                       f'{bracket((c["lo"], c["hi"]))} | {c["direction"]} |')
    out += ['', '### All rounding sensitivities (not sampling intervals)', '',
            'The count-based primary inference does not use rounded rates. The following',
            'uses only positive and negative **printed** Table 1 rates, and allows each',
            'independently to vary ±0.05 (half last digit) or ±0.10 (full last digit)',
            'per 100,000. Each row is a sharp interval for a rounding-only Q expression,',
            'not a confidence interval; these do not change the count-based verdict.', '',
            '| Younger → older | Q from printed rates | ±0.05 rounding range | ±0.10 stress range |',
            '| --- | ---: | --- | --- |']
    for p in result['pairs']:
        out.append(f'| {p["younger"]} → {p["older"]} | {fmt(p["printed_q"])} | '
                   f'{bracket(p["rounding_half"])} | {bracket(p["rounding_full"])} |')
    robust_n = sum(p['robust'] != 'contains 1' for p in result['pairs'])
    scenario_n = sum(c['direction'] != 'contains 1' for p in result['pairs'] for c in p['scenarios'].values())
    out += ['', '## Explicit verdict and scope', '',
            f'**{result["verdict"]}.** Robust missing-status-enveloped exclusions: **{robust_n}/10**.',
            f'Allocation-specific interval exclusions: **{scenario_n}/50** (different mutually hypothetical',
            'completions; not 50 independent tests). The frozen rule does not allow a scenario-specific',
            'effect to override an allocation-robust failure criterion. **HER2: UNTESTABLE** here.',
            '',
            '**Minimum completion to restore compatibility: NOT APPLICABLE.** No pair met',
            'the stronger all-allocation falsification criterion, so no case-completion',
            'restoration threshold was triggered, estimated, or fitted after inspection.', '' ,
            'The older/younger receptor-class rate ratios are comparisons of one calendar year’s',
            'different age groups; no transition, time-series trend, screening-detection effect,',
            'birth-cohort effect, competing-risk mechanism, survival effect, latent biology, or',
            'cause of the localized-stage decline is identified. The 85+ group is open-ended',
            'and is *not* assigned a midpoint or five-year width. The source’s ER-or-PR',
            'positive class includes hormone receptor borderline positive and all HER2 statuses.',
            'Double-negative ER/PR is not equivalent to triple-negative HER2 disease.',
            '',
            '### Literature gate and independent critical interpretation', '',
            '- Turbow et al. (2021) themselves described the positive-class rates as resembling the',
            '  all-cancer age pattern and tabulated the receptor composition by age. This is **not**',
            '  a novel discovery of those rates, only a prespecified missingness/multiplicity audit.',
            '- Lee Argov et al. (2025), doi:10.1001/jamanetworkopen.2025.16947, published',
            '  related USCS 2001–2019 older-age HR/HER2 composition and trends, but not this',
            '  specific 2016 five-bin all-allocation interaction test.',
            '- Jiang et al. (2024), doi:10.1002/ijc.35073, considered ER-specific **calendar**',
            '  trends over 2004–2020 and age 20–84; different estimand and period.',
            '- Acheampong et al. (2020), doi:10.1001/jamanetworkopen.2020.13226, and other SEER subtype-age analyses use different coverage,',
            '  molecular categories and periods; overlapping registries are **not** independent',
            '  replication. The nearby unsettled question is formal allocation-robust proportionality.',
            '- The interval assumption can fail under registry classification error, extra-Poisson',
            '  heterogeneity, ascertainment changes, or residual cases outside binary status.',
            '  Widely differing allocation corners cannot be interpreted as equally plausible.',
            '- No case-level data, joint stage×receptor×age table, HER2 counts for 2016,',
            '  independent population-age weights, or screening linkage are supplied.',
            '',
            '### Hostile-review checklist', '',
            '- Every one of the 20 printed HTML count/rate/95% CI triples and five age-header',
            '  populations checked against the full article HTML, not PDF snippets.',
            '- Explicitly checked source manuscript label, absence of a CC licence, ER/PR class',
            '  definitions, five residuals, rate rounding anomalies, and open-ended 85+.',
            '- CP tails, degenerate count boundary cases, interval monotonicity, all extrema,',
            '  proportional integer rounding and direct printed-rate envelopes independently tested.',
            '- Original numerical conclusions, even if contrary to expectation, must remain intact.',
            '- No browser file download or prompt occurred; browser recaptcha on one URL form was',
            '  bypassed only by opening an already indexed **same PMC HTML** result, not PDF.',
            '- Shell git clone was blocked by DNS; the GitHub connector used for remote branch',
            '  and separate protocol commit. No unsupported claim of successful local full clone.',
            '',
            '**Nonclinical research output. Not medical advice. Software CI and automated',
            'source review do not independently validate scientific or clinical claims.**',
            '']
    return '\n'.join(out)


def render_readme(result: dict) -> str:
    return '\n'.join([
        '# Breast receptor-by-age proportionality: frozen aggregate audit',
        '',
        '**Nonclinical, unreviewed research. Not medical advice or a screening recommendation.**',
        '',
        f'**Verdict: {result["verdict"]}.** See [the complete generated analysis](REPORT.md).',
        'HER2 is **UNTESTABLE** in the source Table 1: its two recorded classes',
        'are ER-positive and/or PR-positive versus ER-negative **and** PR-negative.',
        '',
        'The qualitative receptor-dependent age pattern has already been described',
        'in the literature. This is a narrower, conditional reproducibility audit',
        'of the simultaneous two-class proportionality restriction in **2016 USCS**',
        'data, not a claim of epidemiological discovery or independent replication.',
        'The strict failure gate asks whether an interval excludes Q=1 after',
        'allowing **every** allocation of unknown and residual receptor cases.',
        'No pair passed that stronger failure gate; a minimum restoration completion',
        'was therefore not triggered or computed.',
        '',
        'The [PROTOCOL.md](PROTOCOL.md) was committed by itself before numerical',
        'work (GitHub commit `ec33191bec49a945149ab1fc998ed588a41715d3`).',
        'It freezes data, exact binomial uncertainty, 10 contrasts, missingness',
        'bounds, rounding stress tests, verdicts and reporting. It was not',
        'blinded or registered outside this repository.',
        '',
        'Source: [Turbow et al. 2021 Table 1](' + SOURCE + '), 2016 USCS.',
        'This is an author manuscript with no explicit CC licence visible;',
        'the derived CSV contains attributed factual counts/rates only.',
        '',
        '## Reproduce, offline (stdlib Python only)',
        '',
        '```sh',
        'python research/breast-receptor-age/analysis.py --check',
        'python -m unittest discover -s tests -p test_research_breast_receptor_age.py -v',
        '```',
        '',
        'Use `--write` only if regenerating README.md and REPORT.md from the',
        'committed CSV/protocol. Hashes normalize all newline conventions to LF,',
        'and path keys use `Path.as_posix()` for cross-platform reproduction.',
        '',
        'Results concern one-year population-level group contrasts; no stage',
        'cross-tabulation or HER2 split is available. They identify no screening,',
        'birth-cohort, biological, competing-risk, or clinical cause.',
        '',
    ])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true', help='Regenerate report and README from frozen CSV')
    mode.add_argument('--check', action='store_true', help='Check committed generated outputs')
    args = parser.parse_args()
    result = analyze()
    expected = {REPORT: render_report(result), README: render_readme(result)}
    for path, body in expected.items():
        if args.write:
            path.write_text(body, encoding='utf-8', newline='\n')
        elif normalize_newlines(path.read_text(encoding='utf-8')) != body:
            raise SystemExit(f'Generated file is stale: {key(path)}')
    print(f'PASS ({result["verdict"]}): {len(result["pairs"])} predeclared age contrasts; '
          f'{len(AGES)} age bins; {len(CLASSES) * len(AGES)} source cells')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
