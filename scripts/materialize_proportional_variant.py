"""Materialize a body whose SHAPE differs from the source subject's, not only its size.

    python scripts/materialize_proportional_variant.py --sex female
    python scripts/materialize_proportional_variant.py --pelvis-breadth-ratio 1.10

``materialize_stature_variant.py`` applies one isotropic factor.  That is exact,
and it is also the reason it cannot change proportions: the fitted muscle path
polynomials can only be carried across an isotropic change, because path length
is homogeneous of degree one in the geometry only when the whole geometry moves
together.  A wider pelvis moves 56 of 98 paths and leaves 42 alone.

So this script refits.  ``scripts/native_polynomial_path_fit.cpp`` drives
OpenSim's ``PolynomialPathFitter`` over the shaped model's real ``GeometryPath``s,
wrap objects included, and the result is gated against the shaped model's own
paths on a held-out trajectory -- with the refit noise floor, measured by
``verify_path_refitting.py``, as the baseline that says whether a difference is
real.  A refit takes two to six minutes.

**Order of operations, and why it is that order.**  Every proportional parameter
is defined RELATIVE TO STATURE, because that is how the anthropometry measures
it.  So: apply the shape factors, measure what the shaped body's stature and
proportions actually became, solve the fixed point (the legs are 48% of stature,
so shortening them shortens stature and moves the very denominator the request
was written against), then apply ONE isotropic factor to land stature exactly on
the request.  Nothing here assumes a factor equals a ratio.

**What this does not do.**  It does not make a female body.  ``--sex female``
gives male anatomy at female-typical proportions and the registration says so in
those words.  The entity set is unchanged: 0 female-specific entities, 0 mammary
glands in either sex, and the male genital tract still bound to the pelvis.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.body_parameters import (  # noqa: E402
    LEG_SEGMENT_BODIES, MECHANICAL_STATURE_M, resolve)
from ihm.native import path_refitting as pf  # noqa: E402
from ihm.native.anisotropic_scaling import scale_model_anisotropic  # noqa: E402
from ihm.native.body_measurements import compare, measure  # noqa: E402
from ihm.native.model_scaling import (  # noqa: E402
    UNSCALED_BY_DESIGN, path_point_polyline_lengths, scale_contact_geometry_set,
    scale_model)
from ihm.native.moment_arm_control import FittedMomentArms  # noqa: E402

BASE = 'data/derived/mechanics/whole_body_lumbar_current'
RAW = 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking'
PATHSET = 'subject_walk_scaled_FunctionBasedPathSet.xml'
CONTACT = 'subject_walk_scaled_ContactGeometrySet.xml'
COORDINATES = 'coordinates.sto'
PATH_OVER_FIBER_BAND = (0.2, 20.0)

#: Measured by ``scripts/verify_path_refitting.py`` on this model: two
#: independent refits differ by this much in path length RMS, because
#: ``LatinHypercubeDesign`` seeds its shuffle from ``std::random_device``.
#: Nothing measured on a refit means anything unless it clears this.
REFIT_NOISE_FLOOR_M = 3.12e-5

#: Bodies whose contact spheres ride a segment that ``leg_length_ratio`` scales.
CONTACT_BODIES = ('calcn_r', 'calcn_l', 'toes_r', 'toes_l')


def solve_shape_factors(model_bytes, requested, *, iterations=60, tolerance=1e-13):
    """Raw per-body factors that land the MEASURED ratios on the request.

    Fixed point, on the model, not on an assumption.  ``leg_length_ratio``'s
    denominator is stature and the legs are most of stature, so the naive
    identity "raw factor equals ratio" is wrong by about 1.2% for the female
    preset -- small, and exactly the size of the effect being asked for, which
    is why it is solved rather than absorbed.

    Returns ``(factors, trace)``; ``trace`` carries every iterate so the
    convergence is visible rather than asserted.
    """
    base = _measure_bytes(model_bytes)
    targets = {
        'hip_separation_over_stature': base['hip_separation_over_stature'] * requested['pelvis'],
        'acromial_over_stature': (base['acromial_separation_m'] / base['stature_m']
                                  * requested['shoulder']),
        'leg_over_stature': base['leg_over_stature'] * requested['leg'],
    }
    raw = {'pelvis': requested['pelvis'], 'shoulder': requested['shoulder'],
           'leg': requested['leg']}
    trace = []
    for _ in range(iterations):
        got = _measure_bytes(_apply_shape(model_bytes, raw)[0])
        achieved = {
            'hip_separation_over_stature': got['hip_separation_over_stature'],
            'acromial_over_stature': got['acromial_separation_m'] / got['stature_m'],
            'leg_over_stature': got['leg_over_stature'],
        }
        residual = {k: achieved[k] / targets[k] - 1 for k in targets}
        trace.append({'raw': dict(raw), 'achieved': achieved, 'residual': residual})
        if max(abs(v) for v in residual.values()) < tolerance:
            break
        raw['pelvis'] /= 1 + residual['hip_separation_over_stature']
        raw['shoulder'] /= 1 + residual['acromial_over_stature']
        raw['leg'] /= 1 + residual['leg_over_stature']
    return raw, {'targets': targets, 'base': {k: base[k] for k in
                                              ('hip_separation_over_stature',
                                               'leg_over_stature', 'stature_m')},
                 'iterations': trace}


def _apply_shape(model_bytes, raw):
    factors = {}
    if raw['pelvis'] != 1.0:
        factors['pelvis'] = (1.0, 1.0, raw['pelvis'])
    if raw['shoulder'] != 1.0:
        factors['torso'] = (1.0, 1.0, raw['shoulder'])
    if raw['leg'] != 1.0:
        for body in LEG_SEGMENT_BODIES:
            factors[body] = (raw['leg'],) * 3
    if not factors:
        return model_bytes, {'factors': {}, 'anisotropic_scale': False}
    return scale_model_anisotropic(model_bytes, factors)


def _measure_bytes(model_bytes):
    scratch = ROOT / 'data/derived/.proportional-measure.osim'
    scratch.parent.mkdir(parents=True, exist_ok=True)
    scratch.write_bytes(model_bytes)
    try:
        return measure(scratch)
    finally:
        scratch.unlink(missing_ok=True)


def per_muscle_fibre_scale(base_model, variant_model, work):
    """OpenSim's own rule: fibre and tendon lengths follow the muscle's own path.

    ``Muscle::extendPostScale`` scales ``optimal_fiber_length`` and
    ``tendon_slack_length`` by the ratio of the scaled to the unscaled
    musculotendon length.  Under an isotropic change that ratio is one number;
    under a shape change it is 98 different numbers, and using a single one
    would put every muscle but at most one on the wrong operating point of its
    force-length curve.

    Measured from each model's own ``GeometryPath`` at its own default pose --
    the real path with its wrap objects, not a polynomial approximation of it --
    so every muscle gets a factor, including the 18 the shipped 80-path set
    never covered.
    """
    def at_default(model, name):
        table = work / (name + '_pose.sto')
        default_pose_table(model, table)
        series = pf.sample(model, table, work / (name + '_pose.csv'),
                           log=work / (name + '_pose.log'))
        return {key[0].rsplit('/', 1)[-1]: values[0]
                for key, values in series.items() if key[1] == 'length'}

    before = at_default(base_model, 'fibre_base')
    after = at_default(variant_model, 'fibre_variant')
    shared = sorted(set(before) & set(after))
    if not shared:
        raise ValueError('No muscle is present in both models')
    return {name: after[name] / before[name] for name in shared}


def default_pose_table(model_path, destination):
    """A one-row coordinate trajectory at a model's declared default pose."""
    root = ET.parse(model_path).getroot()
    labels, values = ['time'], [0.0]
    for joint in root.find('.//JointSet/objects'):
        for coordinate in joint.iter('Coordinate'):
            labels.append('/jointset/%s/%s/value' % (joint.get('name'), coordinate.get('name')))
            values.append(float(coordinate.findtext('default_value')))
    pf.write_sto(destination, labels, [values])
    return len(labels) - 1


def scale_catalog_per_muscle(catalog, fibre_scale, global_scale, force_scale):
    """Per-muscle length scaling; a muscle with no refitted path is named, not guessed."""
    out, missing = [], []
    for entry in catalog:
        entry = dict(entry)
        factor = fibre_scale.get(entry['id'])
        if factor is None:
            factor = global_scale
            missing.append(entry['id'])
        for key in ('optimal_fiber_length_m', 'tendon_slack_length_m'):
            if entry.get(key) is not None:
                entry[key] = float(entry[key]) * factor
        if entry.get('max_isometric_force_n') is not None:
            entry['max_isometric_force_n'] = float(entry['max_isometric_force_n']) * force_scale
        if entry.get('path_points'):
            entry['path_points'] = [point for point in entry['path_points']]
        out.append(entry)
    return out, sorted(missing)


def materialize(output, request, *, threads=6, keep_work=False):
    output = Path(output).resolve()
    if not output.is_relative_to(ROOT / 'data/derived'):
        raise ValueError('Owned derived output required')
    if output.exists() and any(output.iterdir()):
        raise ValueError('Refusing to overwrite an existing variant')

    parameters = resolve(request)
    p = parameters['parameters']
    requested = {'pelvis': p['pelvis_breadth_ratio'],
                 'shoulder': p['shoulder_breadth_ratio'],
                 'leg': p['leg_length_ratio']}
    if all(v == 1.0 for v in requested.values()):
        raise ValueError('No proportional change requested; use '
                         'materialize_stature_variant.py for a pure size change')

    inputs = {}

    def read(relative):
        raw = (ROOT / relative).read_bytes()
        inputs[str(relative)] = hashlib.sha256(raw).hexdigest()
        return raw

    registration = json.loads(read(Path(BASE) / 'registration.json'))
    model_bytes = read(registration['model_path'])
    catalog = json.loads(read(registration['catalog_path']))
    contact_bytes = read(Path(RAW) / CONTACT)
    read(Path(RAW) / PATHSET)
    for path, sha in registration.get('sources', {}).items():
        read(path)
        if inputs[path] != sha:
            raise ValueError('Base donor identity mismatch: ' + path)

    base_measurements = _measure_bytes(model_bytes)

    # --- 1. shape ---------------------------------------------------------
    raw, solver_trace = solve_shape_factors(model_bytes, requested)
    shaped_bytes, shape_report = _apply_shape(model_bytes, raw)
    shaped = _measure_bytes(shaped_bytes)

    # --- 2. size ----------------------------------------------------------
    global_scale = p['stature_m'] / shaped['stature_m']
    force_scale = parameters['derived']['muscle_force_scale']
    final_bytes, size_report = scale_model(shaped_bytes, global_scale, force_scale=force_scale)

    output.mkdir(parents=True, exist_ok=True)
    model_path = output / 'subject_proportional.osim'
    model_path.write_bytes(final_bytes)
    final = measure(model_path)

    # Contact spheres ride calcn and toes, which the leg factor scales, and then
    # the global factor scales everything.  Both, or the feet are the wrong size.
    scaled_contact, spheres = scale_contact_geometry_set(
        contact_bytes, raw['leg'] * global_scale)
    (output / CONTACT).write_bytes(scaled_contact)

    # --- 3. refit ---------------------------------------------------------
    work = output / 'refit'
    fitted = pf.fit(model_path, ROOT / RAW / COORDINATES, work,
                    threads=threads, log=output / 'refit.log')
    shutil.copyfile(fitted, output / PATHSET)

    # --- 4. gates ---------------------------------------------------------
    gates, evidence = gate(model_path, output / PATHSET, ROOT / registration['model_path'],
                           ROOT / RAW / PATHSET, base_measurements, final,
                           requested, catalog, work, threads)

    fibre_scale = per_muscle_fibre_scale(
        ROOT / registration['model_path'], model_path, work)
    scaled_catalog, without_paths = scale_catalog_per_muscle(
        catalog, fibre_scale, raw['leg'] * global_scale, force_scale)
    (output / 'catalog.json').write_text(json.dumps(scaled_catalog, indent=2) + '\n')

    if not keep_work:
        shutil.rmtree(work, ignore_errors=True)

    def sha(name):
        return hashlib.sha256((output / name).read_bytes()).hexdigest()

    record = dict(
        schema='ihm.proportional-variant.v1',
        model_path=str((output / 'subject_proportional.osim').relative_to(ROOT)),
        model_sha256=sha('subject_proportional.osim'),
        catalog_path=str((output / 'catalog.json').relative_to(ROOT)),
        catalog_sha256=sha('catalog.json'),
        geometric_scale=global_scale,
        anisotropic_scale=True,
        raw_shape_factors=raw,
        requested_ratios=requested,
        shape_solver=solver_trace,
        source_overrides={
            PATHSET: dict(path=str((output / PATHSET).relative_to(ROOT)), sha256=sha(PATHSET)),
            CONTACT: dict(path=str((output / CONTACT).relative_to(ROOT)), sha256=sha(CONTACT))},
        base_registration=str(Path(BASE) / 'registration.json'),
        base_model_path=registration['model_path'],
        body_parameters=parameters,
        base_measurements=base_measurements,
        shaped_measurements=shaped,
        final_measurements=final,
        measured_change=compare(base_measurements, final),
        shape_report={k: v for k, v in shape_report.items() if k != 'wraps'},
        wrap_scaling=shape_report.get('wraps', []),
        size_report=size_report,
        contact_spheres_scaled=spheres,
        fitted_paths=len(FittedMomentArms(output / PATHSET).paths),
        muscles_without_a_refitted_path=without_paths,
        paths_the_shipped_set_never_covered=sorted(
            set(FittedMomentArms(output / PATHSET).paths)
            - set(FittedMomentArms(ROOT / RAW / PATHSET).paths)),
        fitted_path_coverage_note=(
            'The refit covers every PathActuator in the model, which is 98 here '
            'against the shipped set\'s 80. That is a change with two sides: the '
            '18 arm and trunk muscles the shipped set omitted now have fitted '
            'moment arms, so JointPosturalController can actuate them, AND they '
            'lose their explicit GeometryPath in the engine, which is the '
            'geometry a retinaculum or a fascia would need to constrain. Stated '
            'rather than chosen silently.'),
        per_muscle_fibre_scale=fibre_scale,
        gates=gates,
        gate_evidence=evidence,
        unscaled_by_design=list(UNSCALED_BY_DESIGN),
        sources=inputs,
        materializer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        default_enabled=False, native_acceptance_complete=False,
        scope=parameters['derived']['sex_realisation']['honest_summary'],
        honesty=('This variant differs from the source subject in SIZE and in '
                 'PROPORTION. Its anatomical entity set is unchanged and male: 0 '
                 'female-specific entities, 0 mammary glands in either sex, and '
                 'the male genital tract still bound to the pelvis. No native '
                 'run has been accepted on it.'))
    (output / 'registration.json').write_text(json.dumps(record, indent=2) + '\n')
    if not all(g['passed'] for g in gates):
        raise ValueError('Proportional variant failed its own gates; registration '
                         'written for inspection but the variant is not usable')
    return record


def gate(model_path, pathset, base_model, base_pathset, base_measurements,
         final, requested, catalog, work, threads):
    """Everything that must hold, measured on the artifacts that were written."""
    results, evidence = [], {}

    # 1-3. The proportions the request asked for, as known answers.
    stature_target = final['stature_m']
    for name, key, factor in (
            ('pelvic breadth over stature', 'hip_separation_over_stature', requested['pelvis']),
            ('leg length over stature', 'leg_over_stature', requested['leg'])):
        expected = base_measurements[key] * factor
        got = final[key]
        results.append(dict(gate=name, expected=expected, measured=got,
                            relative_error=abs(got / expected - 1), tolerance=1e-9,
                            note='solved on the model, not assumed from the factor'))
    results.append(dict(
        gate='shoulder breadth over stature',
        expected=(base_measurements['acromial_separation_m'] / base_measurements['stature_m']
                  * requested['shoulder']),
        measured=final['acromial_separation_m'] / final['stature_m'],
        relative_error=abs((final['acromial_separation_m'] / final['stature_m'])
                           / (base_measurements['acromial_separation_m']
                              / base_measurements['stature_m'] * requested['shoulder']) - 1),
        tolerance=1e-9, note='1.0 in every sex preset: no measured value exists'))

    # 4. The refit tracks THIS body's own muscle paths, on held-out frames.
    held_out = work / 'gate_held_out.sto'
    frames = pf.held_out_trajectory(ROOT / RAW / COORDINATES, held_out, stride=50, offset=5)
    truth = pf.sample(model_path, held_out, work / 'gate_truth.csv')
    refit = pf.sample(model_path, held_out, work / 'gate_refit.csv', pathset=pathset)
    tracking = pf.compare(refit, truth, 'length')
    arms = pf.compare(refit, truth, 'moment_arm')
    results.append(dict(
        gate='refitted paths track this body\'s GeometryPaths',
        expected=0.0, measured=tracking['rms_m'],
        relative_error=tracking['rms_m'], tolerance=1e-3,
        note=('RMS metres over %d held-out frames against the variant\'s own '
              'GeometryPaths, wrap objects included. The fitter\'s own printed '
              'RMS is a training error and is not used. Worst muscle %s at '
              '%.3g m; moment arm RMS %.3g m.'
              % (frames, tracking['worst_actuator'], tracking['max_actuator_rms_m'],
                 arms['rms_m']))))

    # 5. The refit was NECESSARY. Compare against what a coefficient scale would
    #    have produced: the base path set times the same global factor. If that
    #    were within the refit noise floor, this whole mechanism would be
    #    unnecessary, and the gate would be saying so.
    scaled_base = pf.sample(base_model, held_out, work / 'gate_base_truth.csv')
    coefficient_scaled = {key: [v * (final['stature_m'] / base_measurements['stature_m'])
                                for v in values]
                          for key, values in scaled_base.items()}
    necessity = pf.compare(truth, coefficient_scaled, 'length')
    results.append(dict(
        gate='a coefficient scale could not have done this',
        expected='> %g m (the refit noise floor)' % REFIT_NOISE_FLOOR_M,
        measured=necessity['rms_m'],
        passed=necessity['rms_m'] > 10 * REFIT_NOISE_FLOOR_M,
        note=('how far this body\'s true muscle paths sit from the isotropically '
              'scaled source subject\'s, which is what a coefficient scale of the '
              'shipped polynomials would have produced. %.0fx the refit noise '
              'floor; worst muscle %s at %.3g m. If this number were small the '
              'proportional change would not have reached the muscles at all.'
              % (necessity['rms_m'] / REFIT_NOISE_FLOOR_M, necessity['worst_actuator'],
                 necessity['max_actuator_rms_m']))))

    # 6. Physiology, from the catalog's own fibre lengths.
    pose = {c.get('name'): float(c.findtext('default_value'))
            for c in ET.parse(model_path).getroot().iter('Coordinate')}
    paths = FittedMomentArms(pathset)
    fibre = {e['id']: float(e['optimal_fiber_length_m']) for e in catalog
             if e.get('optimal_fiber_length_m')}
    excursion = {n: paths.length_and_moment_arms(n, pose)[0] / fibre[n]
                 for n in paths.paths if n in fibre}
    low, high = min(excursion.values()), max(excursion.values())
    results.append(dict(
        gate='path length over optimal fibre length stays in band',
        band=list(PATH_OVER_FIBER_BAND), measured_min=low, measured_max=high,
        passed=PATH_OVER_FIBER_BAND[0] < low < high < PATH_OVER_FIBER_BAND[1],
        note=('necessary and nowhere near sufficient: verify_stature_scaling.py '
              'showed this band passing on a 13% mis-scaled body. It is here to '
              'catch a fit that produced nonsense, not to certify one that did '
              'not. Denominator is the UNSCALED catalog, so this also carries '
              'the per-muscle fibre rescale that follows.')))

    # 7. Polyline: an entirely Python-side representation of the same paths.
    polyline = path_point_polyline_lengths(model_path)
    ratios = {n: paths.length_and_moment_arms(n, pose)[0] / polyline[n]
              for n in paths.paths if n in polyline}
    base_paths = FittedMomentArms(base_pathset)
    base_polyline = path_point_polyline_lengths(base_model)
    base_pose = {c.get('name'): float(c.findtext('default_value'))
                 for c in ET.parse(base_model).getroot().iter('Coordinate')}
    base_ratios = {n: base_paths.length_and_moment_arms(n, base_pose)[0] / base_polyline[n]
                   for n in base_paths.paths if n in base_polyline}
    common = sorted(set(ratios) & set(base_ratios))
    worst = max(common, key=lambda n: abs(ratios[n] - 1))
    results.append(dict(
        gate='polynomial length over path-point polyline stays near unity',
        muscles=len(common),
        measured_min=min(ratios[n] for n in common),
        measured_max=max(ratios[n] for n in common),
        worst_muscle=worst,
        shipped_worst=max(abs(base_ratios[n] - 1) for n in common),
        passed=max(abs(ratios[n] - 1) for n in common)
               <= max(abs(base_ratios[n] - 1) for n in common) + 1e-9,
        note=('computed in Python by forward kinematics over the .osim\'s own '
              'PathPoint elements -- different file section, different code from '
              'anything OpenSim evaluated. It ignores wrapping, so the ratio is '
              'not 1; the gate is that this refit is no further from the '
              'polyline than the shipped set is on the base body.')))

    for record in results:
        if 'passed' not in record:
            record['passed'] = record['relative_error'] <= record['tolerance']
    evidence['held_out_frames'] = frames
    evidence['refit_tracking_worst'] = tracking['per_actuator']
    return results, evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sex', choices=('male', 'female'), default=None)
    parser.add_argument('--stature-m', type=float, default=None)
    parser.add_argument('--mass-kg', type=float, default=None)
    parser.add_argument('--pelvis-breadth-ratio', type=float, default=None)
    parser.add_argument('--shoulder-breadth-ratio', type=float, default=None)
    parser.add_argument('--leg-length-ratio', type=float, default=None)
    parser.add_argument('--threads', type=int, default=6)
    parser.add_argument('--keep-work', action='store_true')
    parser.add_argument('--output', default=None)
    args = parser.parse_args()

    request = {}
    for flag, name in (('sex', 'sex'), ('stature_m', 'stature_m'), ('mass_kg', 'mass_kg'),
                       ('pelvis_breadth_ratio', 'pelvis_breadth_ratio'),
                       ('shoulder_breadth_ratio', 'shoulder_breadth_ratio'),
                       ('leg_length_ratio', 'leg_length_ratio')):
        value = getattr(args, flag)
        if value is not None:
            request[name] = value
    output = args.output or ('data/derived/mechanics/proportional_%s'
                             % (args.sex or 'custom'))
    record = materialize(ROOT / output, request, threads=args.threads,
                         keep_work=args.keep_work)
    print(json.dumps({k: v for k, v in record.items()
                      if k not in ('sources', 'per_muscle_fibre_scale', 'wrap_scaling',
                                   'gate_evidence', 'shape_solver')}, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
