"""Is geometric similarity true of actual people?  Measure it, do not assume it.

    .venv/bin/python scripts/measure_stature_allometry.py

`ihm/body_scaling.py` scales the body by one isotropic linear factor, and every
exponent in it follows from that assumption.  The assumption makes a *testable*
prediction about a population: if a taller person were a scaled-up shorter one,
mass would go as stature cubed and BMI would go as stature.

`BMX_J.parquet` -- the NHANES body measures `scripts/index_anthropometry.py`
already reads for stature and weight percentiles -- can answer that directly, and
the answer is not 3.

This is the same file, the same survey weights and the same adult window as
`index_anthropometry.py`, deliberately: a disagreement between the two would then
be about the regression and not about the sample.  The check with an answer known
in advance is that file's, reused -- the survey-weighted mean stature must land on
the published CDC/NCHS value -- because a weighting error would move the slope and
would otherwise be invisible in a slope.

What it does NOT establish.  This is a between-person association measured across
individuals at one time, not a within-person law, and it is confounded by
adiposity: heavier-for-height people are in the sample and pull the slope. That
confound is not a defect for the question actually being asked, which is "someone
asks for a 2.03 m body -- what should it weigh", because that too is a
between-person question.
"""
from pathlib import Path
import json, sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SOURCE = 'data/derived/population/nhanes-2017-2018'
OUT = 'data/derived/anthropometry/stature-allometry'
AGE_RANGE = (20, 80)

#: The same check with a known answer that scripts/index_anthropometry.py uses.
#: Published survey-weighted means for adults 20 and over, CDC/NCHS.
KNOWN_MEAN_HEIGHT_CM = {'male': (175.3, 1.5), 'female': (161.3, 1.5)}

#: What geometric similarity predicts, before anything is measured.
PREDICTED = {'weight_kg': 3.0, 'bmi_kg_m2': 1.0,
             'waist_circumference_cm': 1.0, 'hip_circumference_cm': 1.0,
             'arm_circumference_cm': 1.0, 'upper_leg_length_cm': 1.0,
             'upper_arm_length_cm': 1.0}

COLUMN = {'weight_kg': 'BMXWT', 'bmi_kg_m2': 'BMXBMI',
          'waist_circumference_cm': 'BMXWAIST', 'hip_circumference_cm': 'BMXHIP',
          'arm_circumference_cm': 'BMXARMC', 'upper_leg_length_cm': 'BMXLEG',
          'upper_arm_length_cm': 'BMXARML'}


def weighted_log_slope(height_cm, value, weight):
    """Survey-weighted slope of log(value) on log(height), with its standard error.

    A log-log slope IS the allometric exponent: value ~ height**b.  The standard
    error is the weighted sandwich form, because the survey weights are not
    counts and treating them as counts would understate it by roughly the design
    effect.
    """
    keep = (np.isfinite(height_cm) & np.isfinite(value) & np.isfinite(weight)
            & (weight > 0) & (height_cm > 0) & (value > 0))
    x, y, w = np.log(height_cm[keep]), np.log(value[keep]), weight[keep]
    if len(x) < 100:
        return None
    p = w / w.sum()
    xb, yb = float(p @ x), float(p @ y)
    sxx = float(p @ (x - xb) ** 2)
    slope = float(p @ ((x - xb) * (y - yb))) / sxx
    residual = y - yb - slope * (x - xb)
    # Weighted heteroskedasticity-consistent variance of the slope.
    variance = float((p ** 2) @ ((x - xb) ** 2 * residual ** 2)) / sxx ** 2
    return dict(n=int(keep.sum()), exponent=slope,
                standard_error=float(np.sqrt(variance)),
                unweighted_exponent=float(np.polyfit(x, y, 1)[0]),
                mean_log_height=xb)


def build(root):
    root = Path(root)
    measures = pd.read_parquet(root / SOURCE / 'BMX_J.parquet')
    demographics = pd.read_parquet(root / SOURCE / 'DEMO_J.parquet')
    frame = measures.merge(
        demographics[['SEQN', 'RIAGENDR', 'RIDAGEYR', 'WTMEC2YR']], on='SEQN')
    adults = frame[(frame.RIDAGEYR >= AGE_RANGE[0]) & (frame.RIDAGEYR < AGE_RANGE[1])]

    groups = {'male': adults[adults.RIAGENDR == 1],
              'female': adults[adults.RIAGENDR == 2],
              'all_adults': adults}

    report = {'schema': 'ihm.stature-allometry.v1',
              'source': {'cycle': 'NHANES 2017-2018 (J)',
                         'tables': ['BMX_J.parquet', 'DEMO_J.parquet'],
                         'weight': 'WTMEC2YR', 'age_years': list(AGE_RANGE)},
              'question': ('Geometric similarity says a taller adult is a scaled-up '
                           'shorter one. Under that assumption mass goes as '
                           'stature**3 and every circumference and segment length '
                           'as stature**1. These are the measured exponents.'),
              'known_answer_check': {}, 'groups': {}}

    for label in ('male', 'female'):
        g = groups[label]
        keep = np.isfinite(g.BMXHT) & np.isfinite(g.WTMEC2YR) & (g.WTMEC2YR > 0)
        h, w = g.BMXHT[keep].to_numpy(), g.WTMEC2YR[keep].to_numpy()
        mean = float(h @ w / w.sum())
        want, tolerance = KNOWN_MEAN_HEIGHT_CM[label]
        report['known_answer_check'][label] = dict(
            weighted_mean_height_cm=mean, published_cm=want, tolerance_cm=tolerance,
            passed=abs(mean - want) <= tolerance,
            note='a weighting error would bias the slope and is invisible in a slope')

    for label, g in groups.items():
        entry = {'participants': int(len(g)), 'exponents': {}}
        for name, column in COLUMN.items():
            fit = weighted_log_slope(g.BMXHT.to_numpy(), g[column].to_numpy(),
                                     g.WTMEC2YR.to_numpy())
            if fit is None:
                continue
            predicted = PREDICTED[name]
            entry['exponents'][name] = dict(
                fit, predicted_by_geometric_similarity=predicted,
                deviation=fit['exponent'] - predicted,
                z_against_prediction=((fit['exponent'] - predicted)
                                      / fit['standard_error']))
            entry['exponents'][name]['consistent_with_isometry'] = bool(
                abs(entry['exponents'][name]['z_against_prediction']) < 2.0)
        report['groups'][label] = entry

    pooled = report['groups']['all_adults']['exponents']['weight_kg']['exponent']
    male = report['groups']['male']['exponents']['weight_kg']['exponent']
    female = report['groups']['female']['exponents']['weight_kg']['exponent']
    report['mixture_check'] = dict(
        pooled=pooled, male=male, female=female,
        pooled_lies_between=bool(min(male, female) <= pooled <= max(male, female)),
        note=('Pooling two groups with different means can manufacture a slope that '
              'is outside both. If the pooled exponent lies between the two '
              'within-sex ones, the pooled number is not that artefact.'))

    from ihm.body_parameters import MECHANICAL_STATURE_M, MECHANICAL_TARGET_MASS_KG
    for stature in (1.40, 1.60, 1.80, 2.03, 2.05):
        s = stature / MECHANICAL_STATURE_M
        report.setdefault('consequence_kg', []).append(dict(
            stature_m=stature, stature_scale=s,
            geometric_similarity_kg=MECHANICAL_TARGET_MASS_KG * s ** 3,
            population_exponent_kg=MECHANICAL_TARGET_MASS_KG * s ** pooled,
            relative_difference=(s ** 3 / s ** pooled - 1.0)))

    report['disposition'] = (
        'Reported, not applied. mass_kg is an independent knob in '
        'ihm/body_parameters.py, and ihm/native/model_scaling.py defaults '
        'mass_scale to s**3 because that is what geometric similarity says. '
        'This file is the measurement that says geometric similarity is wrong '
        'about it, so that a caller who takes the default takes it knowingly.')
    return report


def main():
    report = build(ROOT)
    out = ROOT / OUT
    out.mkdir(parents=True, exist_ok=True)
    (out / 'allometry.json').write_text(json.dumps(report, indent=2) + '\n')
    for label in ('male', 'female', 'all_adults'):
        row = report['groups'][label]['exponents']
        print('%-11s n=%5d  mass ~ stature**%.3f (se %.3f, isometry says 3.00, z=%+.1f)'
              % (label, row['weight_kg']['n'], row['weight_kg']['exponent'],
                 row['weight_kg']['standard_error'],
                 row['weight_kg']['z_against_prediction']))
        print('%-11s          BMI  ~ stature**%.3f (isometry says 1.00, z=%+.1f)'
              % ('', row['bmi_kg_m2']['exponent'], row['bmi_kg_m2']['z_against_prediction']))
    print(json.dumps(report['mixture_check'], indent=2))
    print(json.dumps(report['known_answer_check'], indent=2))
    print(json.dumps(report['consequence_kg'], indent=2))
    print('wrote', out / 'allometry.json')


if __name__ == '__main__':
    main()
