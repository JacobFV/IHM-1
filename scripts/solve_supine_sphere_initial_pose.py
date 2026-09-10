"""Find a lower-acceleration supine sphere-supported pose using native statics.

This is initialization only, not a balance controller or an equilibrium claim.
Every free mobility acceleration remains in the objective and final report.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

import numpy as np
from scipy.optimize import least_squares, brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG


def solve(max_nfev=80, seed_path=None):
    output = Path(tempfile.mkdtemp(prefix='supine-sphere-pose-', dir=ROOT/'data/derived'))
    started = time.monotonic()
    protocol_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    stream = NativeMechanicalStream(ROOT, output/'native', environment='supine',
        target_mass_kg=MECHANICAL_TARGET_MASS_KG,
        augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json')
    def write(name, data):
        (output/name).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    def query(q):
        args = ['evaluate_static_pose', str(len(q))]
        for name, value in q.items():
            args.extend([name, str(float(value))])
        return stream._request(' '.join(args))
    try:
        before = stream.snapshot()
        initial = query({'pelvis_tx': before['coordinates']['pelvis_tx']['value']})
        write('initial.json', initial)
        coordinates = {n: c['value'] for n, c in initial['coordinates'].items() if c['independent']}
        if seed_path:
            seed = json.loads(Path(seed_path).read_text())
            coordinates.update(seed['coordinates'])
        # Plane translations and rotation about the support normal are gauges.
        gauges = ('pelvis_ty', 'pelvis_tz', 'pelvis_list')
        names = [n for n in coordinates if n not in gauges]
        model = ET.parse(output/'native/inputs/subject_walk_scaled.osim')
        source_bounds = {c.attrib['name']: list(map(float, c.findtext('range').split()))
                         for c in model.findall('.//Coordinate')}
        posture_bounds = {n: [max(source_bounds[n][0], -.15), min(source_bounds[n][1], .15)]
                          for n in ('pelvis_tilt','pelvis_list','pelvis_rotation',
                                    'lumbar_extension','lumbar_bending','lumbar_rotation')}
        effective_bounds = {**source_bounds, **posture_bounds}
        coordinates.update({n: float(np.clip(coordinates[n], *b)) for n, b in posture_bounds.items()})
        bounds = np.array([effective_bounds[n] for n in names]).T
        if seed_path is None:
            from build_supine_initial_state import passive_neutral_seed
            neutral = passive_neutral_seed(dict(seed_coordinates=coordinates, coordinate_bounds=source_bounds),
                (output/'native/inputs/subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml').read_bytes())
            coordinates.update({n: v for n, v in neutral['seed_coordinates'].items() if n.startswith('elbow_')})
            write('passive_neutral_seed.json', neutral)
        weight = initial['mass_kg']*np.linalg.norm(initial['gravity_m_s2'])
        # Shift along the support normal only as far as required for total normal force balance.
        # Subsequent all-coordinate solve retains native root accelerations too.
        def force_error(tx):
            return query({**coordinates, 'pelvis_tx': tx})['contact_force_n'][0]-weight
        tx = coordinates['pelvis_tx']
        lower = max(source_bounds['pelvis_tx'][0], tx-.5)
        upper = min(source_bounds['pelvis_tx'][1], tx+.5)
        if force_error(lower)*force_error(upper) < 0:
            coordinates['pelvis_tx'] = float(brentq(force_error, lower, upper, xtol=1e-12))
        else:
            raise ValueError('No bounded vertical support-force root for seed')
        write('preload.json', query(coordinates))
        best = None
        evaluations = 0
        def residual(x):
            nonlocal best, evaluations
            q = {**coordinates, **dict(zip(names, map(float, x)))}
            native = query(q)
            a = np.asarray(native['udot'])
            # Numerical values correspond to explicit 1 rad/s² and 1 m/s² scales.
            r = np.r_[a, 100*np.asarray(native['com_acceleration_m_s2'])]
            cost = float(r@r)
            evaluations += 1
            if best is None or cost < best['cost']:
                best = dict(cost=cost, coordinates=q, native=native, evaluation=evaluations)
                write('best_candidate.json', best)
            return r
        x0 = np.clip(np.array([coordinates[n] for n in names]), bounds[0]+1e-10, bounds[1]-1e-10)
        result = least_squares(residual, x0, bounds=bounds, max_nfev=max_nfev,
            diff_step=1e-4, x_scale='jac', ftol=1e-10, xtol=1e-10, gtol=1e-8)
        native = best['native']
        after = stream._request('observe')
        unchanged = all(before[k] == after[k] for k in ('time_s','coordinates'))
        artifact = dict(schema='ihm.native-initial-pose.v1', environment='supine',
            target_mass_kg=MECHANICAL_TARGET_MASS_KG, coordinates=best['coordinates'],
            protocol_sha256=protocol_sha256, posture_bounds_rad=posture_bounds,
            accepted_equilibrium=False, native_static_residual=native,
            model_sha256=hashlib.sha256((output/'native/inputs/subject_walk_scaled.osim').read_bytes()).hexdigest(),
            scope='Lower acceleration native sphere-supported initialization; no controller, constraint or force-law change; forward verification required.')
        write('initial_pose.json', artifact)
        report = dict(output=str(output), optimizer_success=bool(result.success),
            objective='All mobility accelerations divided by 1 rad/s² or 1 m/s² plus 100 times COM acceleration divided by 1 m/s²; quadratic residual',
            protocol_sha256=protocol_sha256, posture_bounds_rad=posture_bounds,
            supported=bool(abs(native['contact_force_n'][0]-weight)/weight < .01),
            optimizer_message=result.message, evaluations=evaluations, wall_s=time.monotonic()-started,
            continuing_state_unchanged=unchanged, accepted_equilibrium=False,
            initial_max_abs_udot=max(map(abs, initial['udot'])),
            final_max_abs_udot=max(map(abs,native['udot'])),
            initial_acceleration_norm=float(np.linalg.norm(initial['udot'])),
            final_acceleration_norm=float(np.linalg.norm(native['udot'])),
            weight_n=weight, final_contact_force_n=native['contact_force_n'],
            final_com_acceleration_m_s2=native['com_acceleration_m_s2'],
            largest_accelerations=sorted(zip(native['mobility_coordinate_names'], native['udot']), key=lambda x:abs(x[1]),reverse=True)[:10])
        write('report.json', report)
        print(json.dumps(report, indent=2))
    finally:
        stream.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-nfev', type=int, default=80)
    parser.add_argument('--seed', type=Path)
    args = parser.parse_args()
    solve(args.max_nfev, args.seed)
