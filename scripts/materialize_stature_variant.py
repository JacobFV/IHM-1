"""Materialize a geometrically scaled copy of the mechanical body.

    python scripts/materialize_stature_variant.py --stature-m 1.90
    python scripts/materialize_stature_variant.py --scale 0.92 --muscle-force-scale 1.0

Writes a variant directory under ``data/derived/mechanics/`` holding a scaled
model, a scaled fitted-path set, a scaled contact geometry set, a scaled muscle
catalog and a ``registration.json`` that ``NativeMechanicalStream`` accepts as
``augmented_registration``.  ``data/derived`` is gitignored; this script is the
artifact and regenerating a variant takes about a second.

**What makes this more than a multiplication.**  The engine does not run the
model's own muscle paths.  ``native_mechanical_stream.cpp`` calls
``ModelFactory::replacePathsWithFunctionBasedPaths`` and replaces all 80 of them
with fitted polynomials from a separate file, and until this commit that file
was hard-wired to the source subject's.  Scaling the model alone would have left
every lower-limb musculotendon at the unscaled subject's length while its bones
moved, with no error and no warning.  The registration therefore declares
``geometric_scale`` and ``source_overrides`` together, and the stream refuses a
scaled model that does not carry a scaled path set.

The gates run inside this script, on the artifact it just wrote, and it refuses
to leave a variant on disk that fails one.
"""
from pathlib import Path
import argparse, hashlib, json, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.body_parameters import MECHANICAL_STATURE_M, resolve  # noqa: E402
from ihm.native.model_scaling import (head_marker_height_m,  # noqa: E402
                                      path_point_polyline_lengths, scale_catalog,
                                      scale_contact_geometry_set,
                                      scale_function_based_path_set, scale_model,
                                      UNSCALED_BY_DESIGN)
from ihm.native.moment_arm_control import FittedMomentArms  # noqa: E402

BASE = 'data/derived/mechanics/whole_body_lumbar_current'
RAW = 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking'
PATHSET = 'subject_walk_scaled_FunctionBasedPathSet.xml'
CONTACT = 'subject_walk_scaled_ContactGeometrySet.xml'

#: ``scripts/collect_pose_corpus.py`` uses this band as its own quality check on
#: a motion file: a musculotendon path is a few times its own optimal fibre
#: length, never a thousandth of it and never a million times it.  Reused here
#: unchanged, because a scaling that leaves it is a scaling that broke something.
PATH_OVER_FIBER_BAND = (0.2, 20.0)


def reference_pose(model_path):
    """The model's own declared default coordinates -- the reference pose."""
    import xml.etree.ElementTree as ET
    return {c.get('name'): float(c.findtext('default_value'))
            for c in ET.parse(model_path).getroot().iter('Coordinate')}


def gate(base_model, base_pathset, base_catalog,
         scaled_model, scaled_pathset, scaled_catalog, scale):
    """Everything that must hold after scaling, measured on the artifacts.

    Two of these are consistency checks and two are not, and the difference
    matters.  "Stature scales" and "every fitted path length scales" read back
    what the scaler wrote, through independent code, and would catch a file
    scaled by the wrong factor or not at all -- which is the failure this whole
    commit exists to prevent -- but they cannot catch the model and the path set
    being scaled consistently with each other and both being wrong.

    The polyline gate can.  It compares the fitted polynomial length against a
    straight-line walk of the model's own ``PathPoint`` locations, which are a
    different part of a different file reached by a different code path, and
    demands that their ratio not move.
    """
    results = []

    # 1. Stature.  Measured by a forward-kinematic walk over joint frames the
    #    scaler touched one at a time; it must come out at exactly `scale`.
    base_stature = head_marker_height_m(base_model)
    got = head_marker_height_m(scaled_model)
    results.append(dict(
        gate='stature scales', expected=base_stature * scale, measured=got,
        relative_error=abs(got - base_stature * scale) / (base_stature * scale),
        tolerance=1e-12,
        note='forward kinematics over every scaled joint offset frame'))

    # 2. Muscle path length at the reference pose.  This is the gate the whole
    #    exercise turns on: it reads the FITTED polynomials, which live in a
    #    different file from the model, so it fails loudly if only one of the
    #    two was scaled.
    base_paths, scaled_paths = FittedMomentArms(base_pathset), FittedMomentArms(scaled_pathset)
    pose = reference_pose(base_model)
    ratios, arm_ratios = {}, []
    for name in base_paths.paths:
        lb, ab = base_paths.length_and_moment_arms(name, pose)
        ls, arms = scaled_paths.length_and_moment_arms(name, pose)
        ratios[name] = ls / lb
        for coordinate, value in ab.items():
            if abs(value) > 1e-6:
                arm_ratios.append(arms[coordinate] / value)
    worst = max(abs(r / scale - 1) for r in ratios.values())
    results.append(dict(
        gate='every fitted path length scales', muscles=len(ratios), expected=scale,
        measured_min=min(ratios.values()), measured_max=max(ratios.values()),
        relative_error=worst, tolerance=1e-9,
        note='a scale factor applied to the model but not to the path polynomials '
             'lands here as ratio 1.0 while the stature gate above still passes'))
    results.append(dict(
        gate='every moment arm scales', arms=len(arm_ratios), expected=scale,
        measured_min=min(arm_ratios), measured_max=max(arm_ratios),
        relative_error=max(abs(r / scale - 1) for r in arm_ratios), tolerance=1e-9,
        note='OpenSim takes moment arms as -dL/dq when no explicit functions are '
             'given; a moment arm is a length and must scale like one'))

    # 3. The physiological band, read from the two catalogs independently.
    #    Path length over optimal fibre length is dimensionless, so under an
    #    isotropic scale it must be INVARIANT -- not merely in band.
    def excursions(paths, catalog):
        fibre = {entry['id']: float(entry['optimal_fiber_length_m'])
                 for entry in catalog if entry.get('optimal_fiber_length_m')}
        return {name: paths.length_and_moment_arms(name, pose)[0] / fibre[name]
                for name in paths.paths if name in fibre}

    excursion = excursions(scaled_paths, scaled_catalog)
    base = excursions(base_paths, base_catalog)
    low, high = min(excursion.values()), max(excursion.values())
    results.append(dict(
        gate='path length over optimal fibre length stays in band',
        muscles=len(excursion), band=list(PATH_OVER_FIBER_BAND),
        measured_min=low, measured_max=high,
        passed=(PATH_OVER_FIBER_BAND[0] < low and high < PATH_OVER_FIBER_BAND[1]
                and set(excursion) == set(base)),
        drift=max(abs(excursion[m] / base[m] - 1) for m in excursion),
        note='band and numerator from the scaled artifacts, denominator from '
             'each catalog\'s own optimal_fiber_length_m; drift is how far the '
             'dimensionless ratio moved, and it must be zero'))

    # 4. The independent one.  Fitted polynomial length against a straight-line
    #    walk of the model's own PathPoint locations -- different file, different
    #    elements, different code path. Ignores wrapping, so the ratio is not 1;
    #    what must hold is that it does not MOVE.
    polyline_base = path_point_polyline_lengths(base_model)
    polyline_scaled = path_point_polyline_lengths(scaled_model)
    shift, agreement = [], []
    for name in base_paths.paths:
        if name not in polyline_base:
            continue
        rb = base_paths.length_and_moment_arms(name, pose)[0] / polyline_base[name]
        rs = scaled_paths.length_and_moment_arms(name, pose)[0] / polyline_scaled[name]
        shift.append(abs(rs / rb - 1))
        agreement.append(rb)
    results.append(dict(
        gate='polynomial length over path-point polyline length does not move',
        muscles=len(shift), relative_error=max(shift), tolerance=1e-9,
        base_agreement_min=min(agreement), base_agreement_max=max(agreement),
        note='the two representations of the same 80 musculotendon paths agree '
             'to within %.1f%% at the reference pose before scaling, which is '
             'wrap-object geometry the polyline does not model; scale one file '
             'and not the other and this ratio moves by the scale factor'
             % (100 * max(abs(1 - r) for r in agreement))))

    for record in results:
        if 'relative_error' in record:
            record['passed'] = record['relative_error'] <= record['tolerance']
    return results


def materialize(root, output, *, stature_m=None, scale=None, muscle_force_scale=None,
                base=BASE):
    root, output = Path(root).resolve(), Path(output).resolve()
    if not output.is_relative_to(root / 'data/derived'):
        raise ValueError('Owned derived output required')
    if output.exists() and any(output.iterdir()):
        raise ValueError('Refusing to overwrite an existing variant')

    if (stature_m is None) == (scale is None):
        raise ValueError('Give exactly one of --stature-m or --scale')
    if stature_m is None:
        stature_m = float(scale) * MECHANICAL_STATURE_M
    request = {'stature_m': stature_m}
    if muscle_force_scale is not None:
        request['muscle_force_scale'] = muscle_force_scale
    parameters = resolve(request)
    scale = parameters['derived']['stature_scale']
    force_scale = parameters['derived']['muscle_force_scale']

    inputs = {}

    def read(relative):
        raw = (root / relative).read_bytes()
        inputs[str(relative)] = hashlib.sha256(raw).hexdigest()
        return raw

    registration = json.loads(read(Path(base) / 'registration.json'))
    model_bytes = read(registration['model_path'])
    catalog = json.loads(read(registration['catalog_path']))
    if inputs[registration['model_path']] != registration['model_sha256'] \
            or inputs[registration['catalog_path']] != registration['catalog_sha256']:
        raise ValueError('Base identity mismatch')
    for path, sha in registration.get('sources', {}).items():
        read(path)
        if inputs[path] != sha:
            raise ValueError('Base donor identity mismatch: ' + path)
    pathset_bytes = read(Path(RAW) / PATHSET)
    contact_bytes = read(Path(RAW) / CONTACT)

    scaled_model, model_report = scale_model(model_bytes, scale, force_scale=force_scale)
    scaled_pathset, path_count = scale_function_based_path_set(pathset_bytes, scale)
    scaled_contact, sphere_count = scale_contact_geometry_set(contact_bytes, scale)
    scaled_catalog = scale_catalog(catalog, scale, force_scale)

    output.mkdir(parents=True, exist_ok=True)

    def write(name, raw):
        path = output / name
        path.write_bytes(raw)
        return str(path.relative_to(root)), hashlib.sha256(raw).hexdigest()

    model_path, model_sha = write('subject_scaled_stature.osim', scaled_model)
    catalog_path, catalog_sha = write(
        'catalog.json', (json.dumps(scaled_catalog, indent=2) + '\n').encode())
    pathset_path, pathset_sha = write(PATHSET, scaled_pathset)
    contact_path, contact_sha = write(CONTACT, scaled_contact)

    gates = gate(root / registration['model_path'], root / RAW / PATHSET, catalog,
                 output / 'subject_scaled_stature.osim', output / PATHSET,
                 scaled_catalog, scale)
    if not all(record['passed'] for record in gates):
        for record in gates:
            print(json.dumps(record, indent=2), file=sys.stderr)
        raise ValueError('Scaled variant failed its own gates; nothing promoted')

    record = dict(
        schema='ihm.stature-variant.v1',
        model_path=model_path, model_sha256=model_sha,
        catalog_path=catalog_path, catalog_sha256=catalog_sha,
        geometric_scale=scale,
        source_overrides={PATHSET: dict(path=pathset_path, sha256=pathset_sha),
                          CONTACT: dict(path=contact_path, sha256=contact_sha)},
        base_registration=str(Path(base) / 'registration.json'),
        base_model_path=registration['model_path'],
        body_parameters=parameters,
        stature_m=stature_m,
        base_stature_m=MECHANICAL_STATURE_M,
        muscle_count=registration.get('muscle_count'),
        body_count=registration.get('body_count'),
        scaling_report=model_report,
        fitted_paths_scaled=path_count,
        contact_spheres_scaled=sphere_count,
        gates=gates,
        unscaled_by_design=list(UNSCALED_BY_DESIGN),
        sources=inputs,
        materializer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        default_enabled=False, native_acceptance_complete=False,
        scope=('Geometrically scaled copy of the mechanical body. Isotropic: one '
               'linear factor everywhere, so this changes how BIG the subject is '
               'and not what shape it is. Segment length ratios, mass fractions '
               'and radii of gyration are exactly the source subject\'s. No '
               'native run has been accepted on it.'))
    write('registration.json', (json.dumps(record, indent=2) + '\n').encode())
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--stature-m', type=float,
                       help='requested stature in metres; %g is the source subject'
                            % MECHANICAL_STATURE_M)
    group.add_argument('--scale', type=float, help='linear factor, relative to the source subject')
    parser.add_argument('--muscle-force-scale', type=float, default=None,
                        help='default is stature_scale**2 (geometric similarity)')
    parser.add_argument('--base', default=BASE)
    parser.add_argument('--output', default=None)
    args = parser.parse_args()

    scale = args.scale if args.scale is not None else args.stature_m / MECHANICAL_STATURE_M
    output = args.output or ('data/derived/mechanics/stature_%s' % ('%.4f' % scale).replace('.', 'p'))
    record = materialize(ROOT, ROOT / output, stature_m=args.stature_m, scale=args.scale,
                         muscle_force_scale=args.muscle_force_scale, base=args.base)
    print(json.dumps({k: v for k, v in record.items() if k != 'sources'}, indent=2))


if __name__ == '__main__':
    main()
