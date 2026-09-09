"""Sex-stratified adult body proportions from the NHANES body-measures file.

    python scripts/index_anthropometry.py

`BMX_J.parquet` -- standing height, weight, upper leg length, upper arm length,
waist, hip and arm circumference for 8,704 participants -- has been decoded and
sitting on disk since the population index was built, and nothing has ever read
it.  `scripts/index_nhanes.py` loads the table and then only uses the 23
laboratory and vital bindings; stature and segment lengths are in `tables` and
never leave it.  So this repository has had measured human body proportions all
along and has been describing its stature bounds as engineering guesses.

What this produces is a **population statistic, not a body**.  It says how tall
adult men and women are and how their limb lengths sit against their stature.
It does not scale anything, and `docs/BODY_PARAMETERS.md` sets out why the
current isotropic scaler cannot consume the ratios even though it can consume
the stature.

Everything is survey-weighted with `WTMEC2YR`, because NHANES oversamples and an
unweighted mean is a mean over the sample design rather than over the country.
Both are reported so the difference is visible rather than asserted.
"""
from pathlib import Path
import json, sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SOURCE = 'data/derived/population/nhanes-2017-2018'
OUT = 'data/derived/anthropometry/nhanes-2017-2018'

#: NHANES body-measures columns, with the units the codebook gives them.
MEASURES = {
    'BMXHT': ('standing_height_cm', 'cm'),
    'BMXWT': ('weight_kg', 'kg'),
    'BMXBMI': ('bmi_kg_m2', 'kg/m2'),
    'BMXLEG': ('upper_leg_length_cm', 'cm'),
    'BMXARML': ('upper_arm_length_cm', 'cm'),
    'BMXARMC': ('arm_circumference_cm', 'cm'),
    'BMXWAIST': ('waist_circumference_cm', 'cm'),
    'BMXHIP': ('hip_circumference_cm', 'cm'),
}

#: Dimensionless proportions, which are what a *shape* parameter would need.
#: Each is computed per participant and then summarised, never as a ratio of
#: two group means -- those are not the same number.
RATIOS = {
    'upper_leg_over_stature': ('BMXLEG', 'BMXHT'),
    'upper_arm_over_stature': ('BMXARML', 'BMXHT'),
    'waist_over_hip': ('BMXWAIST', 'BMXHIP'),
    'waist_over_stature': ('BMXWAIST', 'BMXHT'),
    'hip_over_stature': ('BMXHIP', 'BMXHT'),
}

#: A case whose answer is known before the code runs. Published survey-weighted
#: means for adults 20 and over, from the CDC/NCHS anthropometric reference
#: report, which pools the 2015-2016 and 2017-2018 cycles. This file uses the
#: 2017-2018 cycle alone, so agreement is expected to be close but not exact;
#: the tolerance is set to admit one cycle's worth of difference and nothing
#: like a weighting error, which would show up as several centimetres.
KNOWN_ANSWER = {
    'male': {'standing_height_cm': (175.3, 1.5), 'weight_kg': (90.6, 3.0)},
    'female': {'standing_height_cm': (161.3, 1.5), 'weight_kg': (77.5, 3.0)},
}

AGE_RANGE = (20, 80)


def weighted_stats(values, weights):
    """Survey-weighted mean, sd and percentiles over the finite pairs."""
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[keep], weights[keep]
    if len(values) < 30:
        return None
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    total = weights.sum()
    mean = float(values @ weights / total)
    variance = float(((values - mean) ** 2) @ weights / total)
    cumulative = (np.cumsum(weights) - 0.5 * weights) / total
    percentiles = {'p%g' % p: float(np.interp(p / 100, cumulative, values))
                   for p in (0.5, 2.5, 5, 25, 50, 75, 95, 97.5, 99.5)}
    return {'n': int(len(values)), 'weighted_mean': mean,
            'weighted_sd': float(np.sqrt(variance)),
            'unweighted_mean': float(values.mean()),
            'percentiles': percentiles}


def cohens_d(a, b):
    """Weighted standardised difference, male minus female."""
    if a is None or b is None:
        return None
    pooled = np.sqrt((a['weighted_sd'] ** 2 + b['weighted_sd'] ** 2) / 2)
    return float((a['weighted_mean'] - b['weighted_mean']) / pooled) if pooled else None


def build(root):
    root = Path(root)
    measures = pd.read_parquet(root / SOURCE / 'BMX_J.parquet')
    demographics = pd.read_parquet(root / SOURCE / 'DEMO_J.parquet')
    frame = measures.merge(
        demographics[['SEQN', 'RIAGENDR', 'RIDAGEYR', 'WTMEC2YR']], on='SEQN')
    adults = frame[(frame.RIDAGEYR >= AGE_RANGE[0]) & (frame.RIDAGEYR < AGE_RANGE[1])].copy()

    for name, (numerator, denominator) in RATIOS.items():
        adults[name] = adults[numerator] / adults[denominator]

    groups = {'male': adults[adults.RIAGENDR == 1],
              'female': adults[adults.RIAGENDR == 2],
              'all_adults': adults}
    report = {'schema': 'ihm.nhanes-anthropometry.v1',
              'source': {'cycle': 'NHANES 2017-2018 (J)',
                         'tables': ['BMX_J.parquet', 'DEMO_J.parquet'],
                         'origin': 'CDC/NCHS; staged by scripts/index_nhanes.py',
                         'weight': 'WTMEC2YR', 'age_years': list(AGE_RANGE)},
              'groups': {}}

    for label, group in groups.items():
        entry = {'participants': int(len(group)), 'measures': {}, 'ratios': {}}
        for column, (name, unit) in MEASURES.items():
            stats = weighted_stats(group[column], group.WTMEC2YR)
            if stats:
                entry['measures'][name] = dict(stats, unit=unit, column=column)
        for name in RATIOS:
            stats = weighted_stats(group[name], group.WTMEC2YR)
            if stats:
                entry['ratios'][name] = dict(stats, unit='dimensionless')
        report['groups'][label] = entry

    male, female = report['groups']['male'], report['groups']['female']
    report['sex_difference'] = {
        'note': ('Male minus female. d is the weighted standardised difference; '
                 'the point of reporting it beside the means is that a large '
                 'absolute difference in a measure can be a small one in shape, '
                 'and the ratios are where that shows.'),
        'measures': {name: {'male': male['measures'][name]['weighted_mean'],
                            'female': female['measures'][name]['weighted_mean'],
                            'difference': (male['measures'][name]['weighted_mean']
                                           - female['measures'][name]['weighted_mean']),
                            'd': cohens_d(male['measures'][name], female['measures'][name])}
                     for name in male['measures'] if name in female['measures']},
        'ratios': {name: {'male': male['ratios'][name]['weighted_mean'],
                          'female': female['ratios'][name]['weighted_mean'],
                          'relative_difference': (male['ratios'][name]['weighted_mean']
                                                  / female['ratios'][name]['weighted_mean'] - 1),
                          'd': cohens_d(male['ratios'][name], female['ratios'][name])}
                   for name in male['ratios'] if name in female['ratios']}}

    checks = []
    for sex, expectations in KNOWN_ANSWER.items():
        for name, (expected, tolerance) in expectations.items():
            got = report['groups'][sex]['measures'][name]['weighted_mean']
            checks.append({'check': '%s %s' % (sex, name), 'published': expected,
                           'measured': got, 'tolerance': tolerance,
                           'passed': abs(got - expected) <= tolerance})
    report['known_answer_checks'] = checks
    report['passed'] = all(check['passed'] for check in checks)
    if not report['passed']:
        raise ValueError('Weighted means disagree with the published reference: '
                         + json.dumps([c for c in checks if not c['passed']]))

    report['limitations'] = [
        'These are population statistics for US adults 20-79 in one NHANES '
        'cycle. They describe a population; the body in this repository is one '
        'specimen and is not a member of it.',
        'Standard errors are not computed. NHANES has a complex design '
        '(SDMVPSU/SDMVSTRA) and a correct variance needs Taylor linearisation '
        'or replicate weights; the weighted sd reported here is the spread of '
        'the population, NOT the uncertainty of the mean.',
        'BMXLEG is upper leg length measured from the inguinal crease to the '
        'proximal tibia and BMXARML is acromion to olecranon. Neither is a '
        'segment length in the OpenSim sense, and neither can be substituted '
        'for one without a correspondence this repository does not have.',
        'Waist and hip circumference are soft-tissue girths. They are the best '
        'shape proxies in this file, and they are still not pelvic bone width: '
        'no catalogued source in this repository measures bi-iliac breadth, '
        'pelvic inlet shape or Q-angle in either sex.',
    ]
    output = root / OUT
    output.mkdir(parents=True, exist_ok=True)
    (output / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    report = build(ROOT)
    checks = report['known_answer_checks']
    print('known-answer checks against the published CDC/NCHS reference')
    for check in checks:
        print('  %-28s published %7.1f   measured %7.1f   %s'
              % (check['check'], check['published'], check['measured'],
                 'ok' if check['passed'] else 'FAIL'))
    print('\nsurvey-weighted, adults %d-%d' % AGE_RANGE)
    print('  %-26s %10s %10s %10s %8s' % ('', 'male', 'female', 'difference', 'd'))
    for name, row in report['sex_difference']['measures'].items():
        print('  %-26s %10.2f %10.2f %10.2f %8.2f'
              % (name, row['male'], row['female'], row['difference'], row['d']))
    print('  proportions')
    for name, row in report['sex_difference']['ratios'].items():
        print('  %-26s %10.4f %10.4f %9.2f%% %8.2f'
              % (name, row['male'], row['female'],
                 100 * row['relative_difference'], row['d']))
    print('\nwritten to %s/summary.json' % OUT)


if __name__ == '__main__':
    main()
