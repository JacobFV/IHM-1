"""Drive and read back OpenSim's PolynomialPathFitter.

The engine substitutes fitted polynomial paths for the model's ``GeometryPath``s
(``ModelFactory::replacePathsWithFunctionBasedPaths``), so any change to the
skeleton that is not a single isotropic factor has to be followed by a REFIT of
those polynomials.  ``scripts/native_polynomial_path_fit.cpp`` is the driver;
this module launches it, reads its output, and provides the comparison that
turns a refit into a gated result rather than a file.

Two facts about the fitter that shape every gate built on it:

1.  **It is not deterministic.**  ``LatinHypercubeDesign::computeRandomHypercube``
    seeds ``std::mt19937`` from ``std::random_device`` on every call, so two runs
    of the same fit on the same model draw different sample sets and produce
    different coefficients.  Reproducing the shipped coefficients bit for bit is
    therefore impossible by construction, and any gate that demanded it would be
    testing the wrong thing.  What can be reproduced -- and is -- is the fitted
    paths' *predictions*.
2.  **Its own summary is a training error.**  The RMS it prints is measured on
    the samples it fitted.  Everything here evaluates on a held-out trajectory
    instead, against the model's own ``GeometryPath`` as ground truth.
"""
from pathlib import Path
import math
import subprocess

ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / 'data/runtime/opensim/native_polynomial_path_fit'

#: The settings that produced the shipped ``subject_walk_scaled_FunctionBasedPathSet.xml``.
#: They are the upstream example's (``OpenSim/Examples/PolynomialPathFitter``),
#: and the shipped file's uniform order 5 with mostly-zero coefficients is the
#: fingerprint of exactly this configuration: stepwise regression at order 5.
UPSTREAM_SETTINGS = dict(
    row_stride=10, max_order=5, min_order=2, samples_per_frame=10, stepwise=1,
    path_length_tolerance=1e-3, moment_arm_tolerance=1e-3,
    global_bounds=(-30.0, 30.0),
    bounds={'/jointset/hip_r/hip_flexion_r': (-15.0, 15.0),
            '/jointset/hip_l/hip_flexion_l': (-15.0, 15.0)},
    fill_sine={'/jointset/mtp_r/mtp_angle_r/value': (0.5, 10.0, 0.0),
               '/jointset/mtp_l/mtp_angle_l/value': (0.5, 10.0, 0.0)},
)


def _require_binary():
    if not BINARY.exists():
        raise RuntimeError('%s is not built; run scripts/build_native_path_fitter.py'
                           % BINARY.relative_to(ROOT))


def read_sto(path):
    """(column labels, list of rows) from an OpenSim Storage file.

    Refuses ``inDegrees=yes`` rather than converting: every consumer here works
    in radians, and a silent unit change is the shape of half this programme's
    withdrawn results.
    """
    lines = Path(path).read_text().splitlines()
    header_end = next(i for i, line in enumerate(lines) if line.strip() == 'endheader')
    header = dict(line.split('=', 1) for line in lines[:header_end] if '=' in line)
    if header.get('inDegrees', 'no').strip() != 'no':
        raise ValueError('%s is in degrees; this module works in radians' % path)
    labels = lines[header_end + 1].split('\t')
    rows = [[float(x) for x in line.split('\t')]
            for line in lines[header_end + 2:] if line.strip()]
    if any(len(row) != len(labels) for row in rows):
        raise ValueError('Ragged storage file %s' % path)
    return labels, rows


def write_sto(path, labels, rows, name='Coordinates'):
    text = ['%s' % name, 'version=3', 'nRows=%d' % len(rows),
            'nColumns=%d' % len(labels), 'inDegrees=no', 'endheader',
            '\t'.join(labels)]
    text += ['\t'.join(repr(v) for v in row) for row in rows]
    Path(path).write_text('\n'.join(text) + '\n')


def held_out_trajectory(source, destination, *, stride, offset):
    """Rows ``offset``, ``offset+stride``, ... of a coordinate trajectory.

    The fits below use ``offset=0``.  Evaluating on a different offset gives
    nominal frames the fit never centred a sample cloud on, which is the closest
    thing to a held-out set this trajectory affords.
    """
    labels, rows = read_sto(source)
    picked = rows[offset::stride]
    if not picked:
        raise ValueError('Held-out selection is empty')
    write_sto(destination, labels, picked)
    return len(picked)


def fit(model, coordinates, output_dir, *, settings=None, threads=6, log=None,
        extra_constant=()):
    """Run PolynomialPathFitter; return the path to the fitted path set."""
    _require_binary()
    settings = dict(UPSTREAM_SETTINGS if settings is None else settings)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    lo, hi = settings['global_bounds']
    argv = [str(BINARY), 'fit',
            '--model', str(model), '--coordinates', str(coordinates),
            '--output-dir', str(output_dir),
            '--row-stride', str(settings['row_stride']),
            '--max-order', str(settings['max_order']),
            '--min-order', str(settings['min_order']),
            '--samples-per-frame', str(settings['samples_per_frame']),
            '--stepwise', str(settings['stepwise']),
            '--path-length-tolerance', repr(settings['path_length_tolerance']),
            '--moment-arm-tolerance', repr(settings['moment_arm_tolerance']),
            '--global-bounds', '%r,%r' % (lo, hi),
            '--threads', str(threads)]
    for coordinate, (lo, hi) in settings.get('bounds', {}).items():
        argv += ['--bounds', '%s,%r,%r' % (coordinate, lo, hi)]
    for column, (amplitude, frequency, phase) in settings.get('fill_sine', {}).items():
        argv += ['--fill-sine', '%s,%r,%r,%r' % (column, amplitude, frequency, phase)]
    for column, value in extra_constant:
        argv += ['--fill-constant', '%s,%r' % (column, value)]
    handle = open(log, 'w') if log else subprocess.DEVNULL
    try:
        subprocess.run(argv, check=True, stdout=handle, stderr=subprocess.STDOUT)
    finally:
        if log:
            handle.close()
    produced = sorted(output_dir.glob('*_FunctionBasedPathSet.xml'))
    if len(produced) != 1:
        raise RuntimeError('Expected one fitted path set in %s, found %d'
                           % (output_dir, len(produced)))
    return produced[0]


def sample(model, coordinates, output, *, pathset=None, moment_arm_coordinates=None,
           log=None):
    """Evaluate lengths and moment arms; return ``{(actuator, quantity, coord): [values]}``.

    With ``pathset`` the fitted polynomials are substituted into the model
    first, so this reads the fitted representation through OpenSim's own
    evaluator rather than through a Python reimplementation of it.  Without it,
    the model's own ``GeometryPath``s answer -- wrap objects included -- and that
    is the ground truth every comparison here is against.
    """
    _require_binary()
    argv = [str(BINARY), 'sample', '--model', str(model),
            '--coordinates', str(coordinates), '--output', str(output)]
    if pathset:
        argv += ['--pathset', str(pathset)]
    if moment_arm_coordinates:
        argv += ['--moment-arm-coordinates', ','.join(moment_arm_coordinates)]
    handle = open(log, 'w') if log else subprocess.DEVNULL
    try:
        subprocess.run(argv, check=True, stdout=handle, stderr=subprocess.STDOUT)
    finally:
        if log:
            handle.close()
    series = {}
    with open(output) as f:
        next(f)
        for line in f:
            row, _time, actuator, quantity, coordinate, value = line.rstrip('\n').split(',')
            series.setdefault((actuator, quantity, coordinate), []).append(float(value))
    if not series:
        raise RuntimeError('Sampler produced no rows')
    return series


def compare(a, b, quantity):
    """Per-actuator RMS difference between two sample dictionaries.

    Returns ``{'per_actuator': {name: rms}, 'rms': overall, 'max': worst}`` in
    metres.  Keys present in only one of the two are an error, not a skip.
    """
    keys = [k for k in a if k[1] == quantity]
    missing = [k for k in keys if k not in b]
    if missing:
        raise ValueError('Sample sets disagree on %d %s keys, e.g. %r'
                         % (len(missing), quantity, missing[0]))
    per_actuator, total, count = {}, 0.0, 0
    for key in keys:
        x, y = a[key], b[key]
        if len(x) != len(y):
            raise ValueError('Sample sets disagree in length for %r' % (key,))
        squares = [(u - v) ** 2 for u, v in zip(x, y)]
        if not all(math.isfinite(s) for s in squares):
            raise ValueError('Nonfinite sample difference for %r' % (key,))
        per_actuator.setdefault(key[0], []).extend(squares)
        total += sum(squares)
        count += len(squares)
    if not count:
        raise ValueError('Nothing to compare for quantity %r' % quantity)
    rms = {name: math.sqrt(sum(s) / len(s)) for name, s in per_actuator.items()}
    return {'per_actuator': rms, 'rms_m': math.sqrt(total / count),
            'max_actuator_rms_m': max(rms.values()),
            'worst_actuator': max(rms, key=rms.get), 'samples': count}
