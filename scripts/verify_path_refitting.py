"""Gate the path refitter, on a body whose answer is known before the code runs.

    python scripts/verify_path_refitting.py [--work DIR] [--threads N] [--reuse]

Refitting the muscle path polynomials is the keystone for every ANISOTROPIC body
parameter -- pelvic breadth, shoulder:hip, limb proportion -- because a single
coefficient factor is only exact under an isotropic scale.  ``docs/BODY_PARAMETERS.md``
called the refitter the smallest item on the blocker list.  It was not missing:
``libosimActuators.so`` in ``data/runtime/opensim/install`` exports
``OpenSim::PolynomialPathFitter``.  What was missing was a driver and a gate.

The gate has to answer a specific question -- *is a refitted path set as good a
representation of this model's muscles as the shipped one?* -- and it has to be
able to fail.  Five arms:

``truth``
    Refit the UNCHANGED model.  Compare the refit's lengths and moment arms
    against the model's own ``GeometryPath``s, on a HELD-OUT trajectory, and
    compare the shipped path set against the same truth.  The shipped set is the
    baseline; the refit must not be materially worse.

``noise``
    Two independent refits of the same model.  The fitter is not deterministic
    (``LatinHypercubeDesign::computeRandomHypercube`` seeds ``std::mt19937`` from
    ``std::random_device``), so its run-to-run spread is a real floor, and no
    later claim of the form "widening the pelvis moved this muscle" means
    anything unless it clears this number.  Reproducing the shipped coefficients
    bit for bit is impossible for the same reason, and this arm is what replaces
    that demand.

``isotropic``
    The known answer.  Scale the model by 1.10 isotropically.  Path length is
    homogeneous of degree one in the geometry at fixed pose, so the true lengths
    of the scaled model are EXACTLY 1.10x the unscaled ones -- checked first, on
    the ``GeometryPath``s themselves.  Then refit the scaled model from scratch
    and require the refit to land on 1.10x the base refit within the noise floor.
    The analytic route and the numerical route share no code; agreeing is the
    strongest statement available that the refit is measuring the body it was
    given.

``polyline``
    The independent-representation check this repository already trusts, pointed
    at the right question.  A straight-line walk of the model's own ``PathPoint``
    locations, computed in Python by forward kinematics over the ``.osim`` -- a
    different file section reached by different code from anything OpenSim
    evaluated.  It ignores wrapping, so its ratio to a true path length is not 1.
    What is checked is that the refit and the model's own ``GeometryPath``s put
    that ratio in the SAME place.  An earlier draft of this arm demanded instead
    that refitting not move the ratio away from the SHIPPED set's, which quietly
    assumed the shipped set was right; the truth arm shows it is 16.7x further
    from the model's own paths than the refit is.

``sabotage``
    A refit of a DIFFERENT body -- the pelvis widened by the measured female
    hip/stature ratio -- scored against the ORIGINAL body's truth.  This is the
    mistake the whole mechanism exists to prevent, injected on purpose, and the
    arm passes only if the gate notices.  A gate that has never been seen to
    fail is not a gate.
"""
from pathlib import Path
import argparse
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.native import path_refitting as pf            # noqa: E402
from ihm.native.model_scaling import (                 # noqa: E402
    default_pose_frames, path_point_polyline_lengths, scale_model)
from ihm.native.moment_arm_control import FittedMomentArms  # noqa: E402
from ihm.native.anisotropic_scaling import scale_model_anisotropic  # noqa: E402

EXAMPLE = ROOT / 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking'
MODEL = EXAMPLE / 'subject_walk_scaled.osim'
COORDINATES = EXAMPLE / 'coordinates.sto'
SHIPPED = EXAMPLE / 'subject_walk_scaled_FunctionBasedPathSet.xml'

#: Measured female/male ratio of hip circumference over stature, NHANES
#: 2017-2018, survey-weighted over 4,883 adults.  See
#: data/derived/anthropometry/nhanes-2017-2018/summary.json.
PELVIS_WIDTH_RATIO = 1.1265
ISOTROPIC_SCALE = 1.10

#: Held-out evaluation: rows 5, 55, 105, ... of the reference trajectory.  The
#: fits centre their sample clouds on rows 0, 10, 20, ...; this offset picks
#: nominal frames none of them sat on.
EVAL_STRIDE, EVAL_OFFSET = 50, 5


def reference_pose(model_path):
    """The model's own default coordinate values, as the fitted paths want them."""
    import xml.etree.ElementTree as ET
    return {c.get('name'): float(c.findtext('default_value'))
            for c in ET.parse(model_path).getroot().iter('Coordinate')}


def default_pose_table(model_path, destination):
    """A one-row coordinate trajectory at the model's declared default pose."""
    import xml.etree.ElementTree as ET
    root = ET.parse(model_path).getroot()
    labels, values = ['time'], [0.0]
    for joint in root.find('.//JointSet/objects'):
        for coordinate in joint.iter('Coordinate'):
            labels.append('/jointset/%s/%s/value' % (joint.get('name'), coordinate.get('name')))
            values.append(float(coordinate.findtext('default_value')))
    pf.write_sto(destination, labels, [values])
    return len(labels) - 1


def run(work, threads, reuse):
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    summary = {'schema': 'ihm.path-refitting-verification.v1', 'arms': {},
               'inputs': {'model': str(MODEL.relative_to(ROOT)),
                          'shipped_pathset': str(SHIPPED.relative_to(ROOT)),
                          'trajectory': str(COORDINATES.relative_to(ROOT))}}

    # ---- the held-out evaluation trajectory ------------------------------
    evaluation = work / 'held_out.sto'
    frames = pf.held_out_trajectory(COORDINATES, evaluation,
                                    stride=EVAL_STRIDE, offset=EVAL_OFFSET)
    summary['inputs']['held_out_frames'] = frames

    def do_fit(name, model_path, settings=None):
        out = work / name
        produced = sorted(out.glob('*_FunctionBasedPathSet.xml')) if out.exists() else []
        if reuse and len(produced) == 1:
            return produced[0]
        started = time.time()
        result = pf.fit(model_path, COORDINATES, out, settings=settings,
                        threads=threads, log=work / (name + '.log'))
        summary.setdefault('fit_seconds', {})[name] = round(time.time() - started, 1)
        return result

    def do_sample(name, model_path, pathset=None):
        out = work / (name + '.csv')
        if reuse and out.exists():
            series = {}
            with open(out) as handle:
                next(handle)
                for line in handle:
                    _row, _t, actuator, quantity, coordinate, value = line.rstrip('\n').split(',')
                    series.setdefault((actuator, quantity, coordinate), []).append(float(value))
            return series
        return pf.sample(model_path, evaluation, out, pathset=pathset,
                         log=work / (name + '.samplelog'))

    # ---- bodies ----------------------------------------------------------
    scaled_model = work / 'subject_iso110.osim'
    if not (reuse and scaled_model.exists()):
        blob, _report = scale_model(MODEL.read_bytes(), ISOTROPIC_SCALE)
        scaled_model.write_bytes(blob)

    wide_model = work / 'subject_widepelvis.osim'
    if not (reuse and wide_model.exists()):
        blob, _report = scale_model_anisotropic(
            MODEL.read_bytes(), {'pelvis': (1.0, 1.0, PELVIS_WIDTH_RATIO)})
        wide_model.write_bytes(blob)

    # ---- fits ------------------------------------------------------------
    fit_a = do_fit('base-a', MODEL)
    fit_b = do_fit('base-b', MODEL)
    fit_iso = do_fit('iso110', scaled_model)
    fit_wide = do_fit('widepelvis', wide_model)

    # ---- samples ---------------------------------------------------------
    truth = do_sample('truth-geometry', MODEL)
    shipped = do_sample('shipped', MODEL, SHIPPED)
    refit_a = do_sample('refit-a', MODEL, fit_a)
    refit_b = do_sample('refit-b', MODEL, fit_b)
    truth_iso = do_sample('truth-iso110', scaled_model)
    refit_iso = do_sample('refit-iso110', scaled_model, fit_iso)
    truth_wide = do_sample('truth-widepelvis', wide_model)
    refit_wide_on_wide = do_sample('refit-wide-on-wide', wide_model, fit_wide)
    refit_wide_on_base = do_sample('refit-wide-on-base', MODEL, fit_wide)

    # ---- arm: truth ------------------------------------------------------
    shipped_error = pf.compare(shipped, truth, 'length')
    refit_error = pf.compare(refit_a, truth, 'length')
    shipped_arms = pf.compare(shipped, truth, 'moment_arm')
    refit_arms = pf.compare(refit_a, truth, 'moment_arm')
    summary['arms']['truth'] = {
        'note': ('refit vs the model\'s own GeometryPaths on %d held-out frames, '
                 'with the shipped path set scored against the same truth as the '
                 'baseline. A raw RMS on its own would say nothing.' % frames),
        'shipped_length_rms_m': shipped_error['rms_m'],
        'refit_length_rms_m': refit_error['rms_m'],
        'shipped_moment_arm_rms_m': shipped_arms['rms_m'],
        'refit_moment_arm_rms_m': refit_arms['rms_m'],
        'refit_over_shipped_length': refit_error['rms_m'] / shipped_error['rms_m'],
        'refit_over_shipped_moment_arm': refit_arms['rms_m'] / shipped_arms['rms_m'],
        'refit_worst_actuator': refit_error['worst_actuator'],
        'refit_worst_actuator_rms_m': refit_error['max_actuator_rms_m'],
        'passed': (refit_error['rms_m'] <= 1.5 * shipped_error['rms_m']
                   and refit_arms['rms_m'] <= 1.5 * shipped_arms['rms_m'])}

    # ---- arm: noise ------------------------------------------------------
    noise_length = pf.compare(refit_a, refit_b, 'length')
    noise_arms = pf.compare(refit_a, refit_b, 'moment_arm')
    summary['arms']['noise'] = {
        'note': ('two independent refits of the SAME model. The fitter draws its '
                 'Latin hypercube from std::random_device, so this is the floor '
                 'under every refit-based claim, and the reason a bitwise '
                 'reproduction of the shipped coefficients is not achievable.'),
        'length_rms_m': noise_length['rms_m'],
        'length_max_actuator_rms_m': noise_length['max_actuator_rms_m'],
        'moment_arm_rms_m': noise_arms['rms_m'],
        'moment_arm_max_actuator_rms_m': noise_arms['max_actuator_rms_m'],
        'passed': noise_length['rms_m'] < 1e-3}

    # ---- arm: isotropic (the answer is known in advance) -----------------
    exact = {key: [v * ISOTROPIC_SCALE for v in values] for key, values in truth.items()}
    geometry_check = pf.compare(truth_iso, exact, 'length')
    predicted = {key: [v * ISOTROPIC_SCALE for v in values] for key, values in refit_a.items()}
    refit_check = pf.compare(refit_iso, predicted, 'length')
    tracks_truth = pf.compare(refit_iso, truth_iso, 'length')
    summary['arms']['isotropic'] = {
        'note': ('scale 1.10 isotropic. True path lengths must be exactly 1.10x, '
                 'and an independent refit of the scaled body must land on 1.10x '
                 'the base refit within the noise floor.'),
        'scale': ISOTROPIC_SCALE,
        'geometrypath_vs_exact_scaling_rms_m': geometry_check['rms_m'],
        'refit_vs_exact_scaling_rms_m': refit_check['rms_m'],
        'refit_vs_scaled_truth_rms_m': tracks_truth['rms_m'],
        'noise_floor_scaled_m': ISOTROPIC_SCALE * noise_length['rms_m'],
        'refit_over_noise_floor': refit_check['rms_m'] / (ISOTROPIC_SCALE * noise_length['rms_m']),
        'geometry_tolerance_m': 1e-5,
        'geometry_tolerance_basis': (
            'Not 0. The sampler calls Model::assemble, and this model carries 4 '
            'CoordinateCouplerConstraints; the assembler solves them to its own '
            'tolerance, and the two models solve them fractionally differently. '
            'The residual is 1.4e-7 m and its worst muscle is gaslat_l, which '
            'crosses the coupled walker knee -- so it is the assembler, not the '
            'scaling. That scale_model is EXACT is established without a solver '
            'by verify_stature_scaling.py, whose identity arm reproduces every '
            'path length to relative error 0.0 straight from the XML.'),
        'passed': (geometry_check['rms_m'] < 1e-5
                   and refit_check['rms_m'] <= 2.0 * ISOTROPIC_SCALE * noise_length['rms_m'])}

    # ---- arm: polyline ---------------------------------------------------
    # The polyline is computed in Python, by forward kinematics over the .osim's
    # own PathPoint elements -- a different file section reached by different
    # code from anything OpenSim evaluated above. It ignores wrap objects, so
    # its ratio to a true path length is not 1. What it can do is say whether
    # the refit sits where the model's own GeometryPaths sit: r_refit and
    # r_truth must agree, and a systematic error anywhere in the OpenSim-side
    # pipeline would separate them.
    polyline = path_point_polyline_lengths(MODEL, default_pose_frames(MODEL))
    pose_table = work / 'default_pose.sto'
    default_pose_table(MODEL, pose_table)

    def at_default_pose(name, pathset=None):
        out = work / (name + '.csv')
        if not (reuse and out.exists()):
            pf.sample(MODEL, pose_table, out, pathset=pathset,
                      log=work / (name + '.samplelog'))
        lengths = {}
        with open(out) as handle:
            next(handle)
            for line in handle:
                _row, _t, actuator, quantity, _c, value = line.rstrip('\n').split(',')
                if quantity == 'length':
                    lengths[actuator.rsplit('/', 1)[-1]] = float(value)
        return lengths

    truth_pose = at_default_pose('pose-truth')
    refit_pose = at_default_pose('pose-refit-a', fit_a)
    shipped_pose = at_default_pose('pose-shipped', SHIPPED)
    common = sorted(set(polyline) & set(refit_pose) & set(shipped_pose))
    ratio = lambda d: {n: d[n] / polyline[n] for n in common}  # noqa: E731
    r_truth, r_refit, r_shipped = ratio(truth_pose), ratio(refit_pose), ratio(shipped_pose)
    refit_drift = {n: r_refit[n] / r_truth[n] - 1 for n in common}
    shipped_drift = {n: r_shipped[n] / r_truth[n] - 1 for n in common}
    worst = max(refit_drift, key=lambda n: abs(refit_drift[n]))
    worst_shipped = max(shipped_drift, key=lambda n: abs(shipped_drift[n]))
    summary['arms']['polyline'] = {
        'note': ('polynomial length over a Python-side path-point polyline at '
                 'the default pose, for the refit and for the GeometryPath '
                 'truth. The ratio is not 1 -- the polyline does not wrap -- and '
                 'that is not what is checked; what is checked is that the refit '
                 'and the truth put it in the same place. The shipped path set '
                 'is scored the same way as the baseline.'),
        'muscles': len(common),
        'truth_mean_ratio': sum(r_truth[n] for n in common) / len(common),
        'refit_mean_ratio': sum(r_refit[n] for n in common) / len(common),
        'shipped_mean_ratio': sum(r_shipped[n] for n in common) / len(common),
        'refit_worst_muscle': worst,
        'refit_worst_drift_from_truth': refit_drift[worst],
        'shipped_worst_muscle': worst_shipped,
        'shipped_worst_drift_from_truth': shipped_drift[worst_shipped],
        'wrap_gap_truth_vs_polyline': max(abs(r_truth[n] - 1) for n in common),
        'passed': (abs(refit_drift[worst]) < 0.02
                   and abs(refit_drift[worst]) <= abs(shipped_drift[worst_shipped]))}

    # ---- arm: sabotage ---------------------------------------------------
    wrong_body = pf.compare(refit_wide_on_base, truth, 'length')
    right_body = pf.compare(refit_wide_on_wide, truth_wide, 'length')
    geometry_moved = pf.compare(truth_wide, truth, 'length')
    summary['arms']['sabotage'] = {
        'note': ('a refit of the WIDE-PELVIS body scored against the ORIGINAL '
                 "body's truth. The gate has to notice; it also has to still pass "
                 'that same refit against its own body.'),
        'pelvis_medial_lateral_factor': PELVIS_WIDTH_RATIO,
        'geometrypath_length_change_rms_m': geometry_moved['rms_m'],
        'geometrypath_worst_actuator': geometry_moved['worst_actuator'],
        'geometrypath_worst_change_m': geometry_moved['max_actuator_rms_m'],
        'wide_refit_against_original_truth_rms_m': wrong_body['rms_m'],
        'wide_refit_against_its_own_truth_rms_m': right_body['rms_m'],
        'noise_floor_m': noise_length['rms_m'],
        'ratio_wrong_to_right': wrong_body['rms_m'] / right_body['rms_m'],
        'passed': (wrong_body['rms_m'] > 10 * noise_length['rms_m']
                   and right_body['rms_m'] <= 1.5 * refit_error['rms_m'])}

    summary['passed'] = all(arm['passed'] for arm in summary['arms'].values())
    summary['findings'] = [
        'PolynomialPathFitter is NOT deterministic: LatinHypercubeDesign seeds '
        'std::mt19937 from std::random_device on every call. Two refits of the '
        'same model differ by %.3g m RMS in path length. The shipped coefficients '
        'therefore cannot be reproduced bit for bit by anyone, including the '
        'people who produced them, and a gate that demanded it would be testing '
        'the wrong property.' % noise_length['rms_m'],
        'The fitter\'s own printed RMS is a TRAINING error, measured on the '
        'samples it fitted. Every number here is measured on %d held-out frames '
        'against the model\'s own GeometryPaths.' % frames,
        'Widening the pelvis medio-laterally by %.4f moves the true GeometryPath '
        'lengths by %.3g m RMS, worst muscle %s at %.3g m. That is %.0fx the '
        'refit noise floor, so the change is measurable and not a fitting '
        'artefact.' % (PELVIS_WIDTH_RATIO, geometry_moved['rms_m'],
                       geometry_moved['worst_actuator'],
                       geometry_moved['max_actuator_rms_m'],
                       geometry_moved['rms_m'] / noise_length['rms_m'])]
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--work', default=str(ROOT / 'data/derived/path-refitting/verify'))
    parser.add_argument('--threads', type=int, default=6)
    parser.add_argument('--reuse', action='store_true',
                        help='reuse fits and samples already in the work directory')
    args = parser.parse_args()
    summary = run(args.work, args.threads, args.reuse)
    print(json.dumps(summary, indent=2))
    return 0 if summary['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
