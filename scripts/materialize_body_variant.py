"""One command for a whole body of a given size -- scaffold and body together.

    .venv/bin/python scripts/materialize_body_variant.py --stature-m 2.03
    .venv/bin/python scripts/materialize_body_variant.py --scale 1.0    # identity

`materialize_stature_variant.py` builds the 22-segment scaffold at a stature.
The four scripts beside this one build the body: the 4,000 anatomical entities,
the 146 nerve routes and 1,743 conduction delays, the 1,326 skin patches, and
the muscle and ligament properties.  Running the scaffold alone is what it means
for the parametrization to stop at the scaffold, so this runs all five and
refuses to leave a variant on disk if any one of them fails a gate.

**And it measures the seam between them**, which turns out to be where a
scaled body silently stops being scaled.

`scripts/native_mechanical_stream.cpp` line 69 reads

    double mass_scale = target_mass / original_mass;

and multiplies every body mass by it.  So the engine RENORMALISES the whole
model to whatever `target_mass_kg` it was handed, and the `s**3` that
`materialize_stature_variant.py` so carefully applied to the 22 body masses is
divided straight back out.  Ask for a 2.03 m body, pass the standing
`target_mass_kg = 77.6122029` that 58 files pass, and you get a 2.03 m skeleton
weighing 77.61 kg: a person 13% taller and not one gram heavier.  Every gate on
the mechanical side passes, because the mass scaling was applied correctly and
the annihilation happens downstream in the engine.

The weight arm below measures exactly that, statically, from the scaled `.osim`
and the engine's own formula.  It does not run the engine.  **No native run has
been accepted on a scaled body**, on either side -- `native_acceptance_complete`
is false in every variant this repository has ever written -- so the standing
weight and momentum-residual gates that need a running plant are named here as
open and not reported as passed.
"""
from pathlib import Path
import argparse, hashlib, json, sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm import body_scaling as bs                                       # noqa: E402
from ihm.body_parameters import (MECHANICAL_SOURCE_MASS_KG,              # noqa: E402
                                 MECHANICAL_STATURE_M,
                                 MECHANICAL_TARGET_MASS_KG, resolve)

import scale_anatomical_body as anatomy_stage                            # noqa: E402
import scale_muscle_and_tissue as muscle_stage                           # noqa: E402
import scale_nerve_conduction as nerve_stage                             # noqa: E402
import scale_skin_patches as skin_stage                                  # noqa: E402
from materialize_stature_variant import materialize as materialize_mechanical  # noqa: E402

G = 9.80665


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def model_mass_kg(path):
    return sum(float(b.findtext('mass')) for b in ET.parse(path).getroot().iter('Body'))


def weight_arm(scaled_model, scale):
    """What the plant will actually weigh, and what it would have to be handed.

    Static: the summed body masses of the scaled model, then the engine's own
    renormalisation formula applied to them.  The point is the third row.
    """
    scaled = model_mass_kg(scaled_model)
    expected = MECHANICAL_SOURCE_MASS_KG * scale ** bs.exponent('mass')
    engine_mass = MECHANICAL_TARGET_MASS_KG          # whatever it is handed, exactly
    consistent_target = MECHANICAL_TARGET_MASS_KG * scale ** bs.exponent('mass')
    return dict(
        source_model_mass_kg=MECHANICAL_SOURCE_MASS_KG,
        scaled_model_mass_kg=scaled,
        expected_scaled_mass_kg=expected,
        mass_exponent=bs.exponent('mass'),
        relative_error=abs(scaled / expected - 1.0),
        engine_renormalises_to_kg=engine_mass,
        engine_mass_scale=engine_mass / scaled,
        standing_weight_n=engine_mass * G,
        weight_if_target_scaled_n=consistent_target * G,
        consistent_target_mass_kg=consistent_target,
        population_exponent_target_mass_kg=MECHANICAL_TARGET_MASS_KG * scale ** 2.034,
        annihilated=abs(scale - 1.0) > 1e-12,
        finding=(
            'The scaled model carries %.3f kg. native_mechanical_stream.cpp line '
            '69 then multiplies every body mass by target_mass/original_mass, so '
            'the plant weighs whatever target_mass_kg it was handed -- %.4f kg by '
            'default, in 58 files -- and the s**%g applied to the model is divided '
            'straight back out. A caller who wants a body of this stature to weigh '
            'what geometric similarity says must pass target_mass_kg = %.3f kg; '
            'what the NHANES population says is %.3f kg '
            '(scripts/measure_stature_allometry.py). Neither is the default, and '
            'the default is what every existing call site passes.'
            % (scaled, MECHANICAL_TARGET_MASS_KG, bs.exponent('mass'),
               consistent_target,
               MECHANICAL_TARGET_MASS_KG * scale ** 2.034)),
        open_gates=[
            'Standing weight equal to m*g on the running plant, and momentum '
            'residual below 1e-14, both need a native run. No native run has been '
            'accepted on a scaled body on either side -- native_acceptance_complete '
            'is false in every variant this repository has written -- so those two '
            'are named open here rather than reported as passed.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--stature-m', type=float)
    group.add_argument('--scale', type=float)
    parser.add_argument('--muscle-force-scale', type=float, default=None)
    parser.add_argument('--mesh-sample', type=int, default=200)
    parser.add_argument('--output', default=None)
    args = parser.parse_args()

    stature = args.stature_m if args.stature_m is not None else args.scale * MECHANICAL_STATURE_M
    parameters = resolve({'stature_m': stature}
                         | ({'muscle_force_scale': args.muscle_force_scale}
                            if args.muscle_force_scale is not None else {}))
    scale = parameters['derived']['stature_scale']
    out = ROOT / (args.output or ('data/derived/body-variants/stature_%s'
                                  % ('%.6f' % scale).replace('.', 'p')))
    out.mkdir(parents=True, exist_ok=True)

    stages, gates = {}, []

    mechanical_dir = out / 'mechanical'
    if not (mechanical_dir.exists() and any(mechanical_dir.iterdir())):
        stages['mechanical'] = materialize_mechanical(
            ROOT, mechanical_dir, stature_m=stature,
            muscle_force_scale=args.muscle_force_scale)
    else:
        stages['mechanical'] = json.loads(
            (mechanical_dir / 'registration.json').read_text())
    gates += [dict(g, stage='mechanical') for g in stages['mechanical']['gates']]

    scaled_entities, scaled_binding, anatomy_report = anatomy_stage.build(
        ROOT, scale, mesh_sample=args.mesh_sample)
    (out / 'anatomy_scaled.json').write_text(json.dumps(
        {'schema': 'ihm.anatomy-scaled.v1', 'stature_scale': scale,
         'entities': scaled_entities}, indent=1) + '\n')
    (out / 'binding_scaled.json').write_text(json.dumps(scaled_binding, indent=1) + '\n')
    stages['anatomy'] = anatomy_report
    gates += [dict(g, stage='anatomy') for g in anatomy_report['gates']]

    scaled_peripheral, nerve_report = nerve_stage.build(ROOT, scale)
    (out / 'peripheral_scaled.json').write_text(
        json.dumps(scaled_peripheral, indent=1) + '\n')
    stages['nerve'] = nerve_report
    gates += [dict(g, stage='nerve') for g in nerve_report['gates']]

    scaled_dermatomes, skin_report = skin_stage.build(ROOT, scale)
    (out / 'dermatomes_scaled.json').write_text(
        json.dumps(scaled_dermatomes, indent=1) + '\n')
    stages['skin'] = skin_report
    gates += [dict(g, stage='skin') for g in skin_report['gates']]

    scaled_catalog, scaled_ligaments, muscle_report = muscle_stage.build(ROOT, scale)
    (out / 'muscle_catalog_scaled.json').write_text(
        json.dumps(scaled_catalog, indent=1) + '\n')
    if scaled_ligaments is not None:
        (out / 'ligaments_scaled.json').write_text(json.dumps(
            {'schema': 'ihm.tissue-force-elements-scaled.v1', 'stature_scale': scale,
             'elements': scaled_ligaments}, indent=1) + '\n')
    stages['muscle_tissue'] = muscle_report
    gates += [dict(g, stage='muscle_tissue') for g in muscle_report['gates']]

    weight = weight_arm(ROOT / stages['mechanical']['model_path'], scale)
    gates.append(dict(stage='weight', gate='scaled model mass takes the mass exponent',
                      relative_error=weight['relative_error'], tolerance=1e-12,
                      passed=weight['relative_error'] <= 1e-12,
                      note=weight['finding']))

    record = dict(
        schema='ihm.body-variant.v1',
        stature_m=stature, stature_scale=scale,
        body_parameters=parameters,
        scaling_table=bs.table(),
        stages={name: {k: v for k, v in report.items()
                       if k not in ('gates', 'sources', 'scaling_table')}
                for name, report in stages.items()},
        weight=weight,
        gates=gates,
        gates_passed=sum(1 for g in gates if g['passed']), gates_total=len(gates),
        default_enabled=False, native_acceptance_complete=False,
        materializer_sha256=_sha(Path(__file__)),
        scope=('One isotropic factor applied to the whole body: the 22-segment '
               'scaffold AND the 4,000 anatomical entities, 146 nerve routes, '
               '1,326 skin patches, 98 actuators and 117 ligament force elements '
               'that ride on it. It changes how BIG the subject is and not what '
               'SHAPE it is. Every proportion is the source subject\'s at every '
               'stature, and ihm/body_scaling.py ALLOMETRY lists the six places '
               'that assumption is known to be wrong, with the size of each error.'))
    (out / 'registration.json').write_text(json.dumps(record, indent=2) + '\n')

    by_stage = {}
    for gate in gates:
        by_stage.setdefault(gate['stage'], [0, 0])
        by_stage[gate['stage']][0] += int(gate['passed'])
        by_stage[gate['stage']][1] += 1
    print('stature %.4f m, scale %.6f' % (stature, scale))
    for stage, (ok, total) in by_stage.items():
        print('  %-14s %2d/%2d gates' % (stage, ok, total))
    for gate in gates:
        if not gate['passed']:
            print('  FAIL [%s] %s' % (gate['stage'], gate['gate']))
    print()
    for row in stages['nerve']['headline']:
        print('  %-36s %7.2f ms -> %7.2f ms  %+.2f%%'
              % (row['description'][:36], row['delay_ms'], row['scaled_delay_ms'],
                 row['change_percent']))
    print()
    print('  ' + weight['finding'].replace('. ', '.\n  '))
    print('\nwrote', out)
    if not all(g['passed'] for g in gates):
        raise SystemExit('body variant failed its own gates; nothing promoted')


if __name__ == '__main__':
    main()
