"""What a proportional change costs, if the fitted muscle paths are not refitted.

    python scripts/measure_anisotropic_pelvis_error.py

The stature knob is isotropic, and isotropy is not a simplification that was
chosen for convenience -- it is the condition under which the fitted path
polynomials can be scaled at all.  Every musculotendon length in this model is a
polynomial in joint angles, and one factor multiplied through it is exactly
right when the whole body scales by that factor and has no defensible value when
different parts scale differently.

The proportional parameters a sex knob would need are anisotropic.  NHANES gives
women a hip circumference 12.6% larger relative to stature than men's
(`scripts/index_anthropometry.py`), and widening a pelvis is a change in one
direction of one body.

So: widen the pelvis medio-laterally, hold everything else, and measure how far
each muscle's path length actually moves.  That number is the error a caller
would incur by making the change and leaving the polynomials alone -- which,
since there is no refitting tool in this repository, is the only thing a caller
could currently do.  It is reported per muscle rather than as a single figure,
because the point is that it is not uniform: it is large across the hip and zero
below the knee, and no single corrective factor exists.

Refitting is what would be needed, and it needs OpenSim: `FunctionBasedPath`
coefficients come from `PolynomialPathFitter` sampling the real `GeometryPath`
with its wrap objects over the coordinate ranges.  The OpenSim Python bindings
are not installed in this environment, and the native adapter that is built
exposes the stream protocol, not the fitter.
"""
from pathlib import Path
import argparse, json, sys
import xml.etree.ElementTree as ET

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.native.model_scaling import path_point_polyline_lengths  # noqa: E402

BASE = 'data/derived/mechanics/whole_body_lumbar_current'
#: Medio-lateral is the model's z axis; the pelvis frame is aligned with ground
#: at the default pose, so widening z widens the pelvis.
AXIS = 2


def widen_pelvis(model_bytes, factor, axis=AXIS):
    """Scale every pelvis-frame coordinate along one axis, and nothing else.

    Touched: the pelvis body's own mass centre, its wrap objects, the markers
    and muscle path points attached to it, and the parent offset frames of the
    joints whose parent is the pelvis -- which is what carries the hip joint
    centres outward.  Nothing distal moves, so the femurs are the same femurs at
    a wider stance width.
    """
    root = ET.fromstring(model_bytes)
    touched = 0

    def scale(node):
        nonlocal touched
        values = [float(x) for x in node.text.split()]
        values[axis] *= factor
        node.text = ' '.join(repr(v) for v in values)
        touched += 1

    for body in root.iter('Body'):
        if body.get('name') != 'pelvis':
            continue
        scale(body.find('mass_center'))
        for wrap in body.iter('WrapObjectSet'):
            for node in wrap.iter('translation'):
                scale(node)
    for element in root.iter():
        parent = element.findtext('socket_parent_frame') or element.findtext('socket_parent')
        if parent not in ('/bodyset/pelvis', 'pelvis'):
            continue
        for tag in ('location', 'translation'):
            node = element.find(tag)
            if node is not None and node.text and node.text.split():
                scale(node)
    return ET.tostring(root, encoding='utf-8', xml_declaration=True), touched


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--factor', type=float, default=None,
                        help='medio-lateral pelvis factor; the default is read '
                             'from the NHANES summary as the female/male ratio '
                             'of hip circumference over stature')
    parser.add_argument('--out', default=None)
    args = parser.parse_args()

    if args.factor is None:
        summary = json.loads(
            (ROOT / 'data/derived/anthropometry/nhanes-2017-2018/summary.json').read_text())
        ratio = summary['sex_difference']['ratios']['hip_over_stature']
        args.factor = ratio['female'] / ratio['male']

    registration = json.loads((ROOT / BASE / 'registration.json').read_text())
    model = (ROOT / registration['model_path']).read_bytes()
    widened, touched = widen_pelvis(model, args.factor)

    work = ROOT / 'data/derived/anisotropic-pelvis-probe'
    work.mkdir(parents=True, exist_ok=True)
    (work / 'widened.osim').write_bytes(widened)

    before = path_point_polyline_lengths(ROOT / registration['model_path'])
    after = path_point_polyline_lengths(work / 'widened.osim')
    shifts = {name: after[name] / before[name] - 1 for name in before}
    moved = {name: value for name, value in shifts.items() if abs(value) > 1e-9}

    ranked = sorted(moved.items(), key=lambda kv: -abs(kv[1]))
    report = {
        'schema': 'ihm.anisotropic-pelvis-error.v1',
        'factor': args.factor,
        'factor_basis': ('NHANES 2017-2018 female/male ratio of hip circumference '
                         'over stature, from data/derived/anthropometry/'
                         'nhanes-2017-2018/summary.json'),
        'elements_moved': touched,
        'muscles': len(before),
        'muscles_whose_path_changed': len(moved),
        'muscles_unaffected': len(before) - len(moved),
        'largest_change': ranked[0][1] if ranked else 0.0,
        'median_change_among_affected': float(np.median([abs(v) for v in moved.values()]))
                                        if moved else 0.0,
        'per_muscle': [{'muscle': name, 'relative_change': value} for name, value in ranked],
        'conclusion': (
            'A single uniform coefficient factor cannot represent this. %d of %d '
            'muscle paths change and %d do not, and among those that do the '
            'change spans %.2f%% to %.2f%%. Any one number applied to the '
            'polynomials is wrong for every muscle but at most one. The fitted '
            'paths would have to be REFITTED, not rescaled.'
            % (len(moved), len(before), len(before) - len(moved),
               100 * min(abs(v) for v in moved.values()) if moved else 0,
               100 * max(abs(v) for v in moved.values()) if moved else 0)),
        'caveat': (
            'Measured on the straight-line path-point polyline, which ignores '
            'wrap objects. Where a hip muscle wraps, the true change differs '
            'from this; the polyline agrees with the fitted polynomial to within '
            '5.1% at the reference pose, so this is the right order of magnitude '
            'and not a substitute for a refit.'),
    }
    (work / 'report.json').write_text(json.dumps(report, indent=2) + '\n')

    print('pelvis widened medio-laterally by %.4f; %d XML elements moved'
          % (args.factor, touched))
    print('%d of %d muscle paths changed length; %d did not'
          % (len(moved), len(before), len(before) - len(moved)))
    print('\nlargest changes')
    for name, value in ranked[:12]:
        print('  %-18s %+7.3f%%' % (name, 100 * value))
    print('\n%s' % report['conclusion'])
    print('written to data/derived/anisotropic-pelvis-probe/report.json')


if __name__ == '__main__':
    main()
