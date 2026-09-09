"""Find a lower-acceleration supine sphere-supported pose using native statics.

This is initialization only, not a balance controller or an equilibrium claim.
Every free mobility acceleration remains in the objective and final report.
The six donor lumbar activation values are additional bounded search variables.
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
from scipy.optimize import least_squares, brentq, lsq_linear

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream


def solve(max_nfev=80, seed_path=None, surface=None, bed=None):
    activation_names = [f"gait2392_{muscle}_{side}" for side in ("r", "l") for muscle in ("ercspn", "intobl", "extobl")]
    output = Path(tempfile.mkdtemp(prefix='supine-lumbar-activation-', dir=ROOT/'data/derived'))
    started = time.monotonic()
    protocol_bytes=Path(__file__).read_bytes()
    protocol_sha256 = hashlib.sha256(protocol_bytes).hexdigest()
    (output/'solver.py').write_bytes(protocol_bytes)
    stream = NativeMechanicalStream(ROOT, output/'native', environment='supine',
        target_mass_kg=77.6122029,
        augmented_registration='data/derived/mechanics/whole_body_lumbar_current/registration.json',
        surface_contact_manifest=surface, bed_material=bed)
    def write(name, data):
        (output/name).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    def query(q, activations=None):
        args = ['evaluate_static_pose', str(len(q))]
        for name, value in q.items():
            args.extend([name, str(float(value))])
        if activations:
            args.append(str(len(activations)))
            for name, value in activations.items():
                args.extend([name, str(float(value))])
        return stream._request(' '.join(args))
    try:
        before = stream.snapshot()
        initial = query({'pelvis_tx': before['coordinates']['pelvis_tx']['value']})
        write('initial.json', initial)
        coordinates = {n: c['value'] for n, c in initial['coordinates'].items() if c['independent']}
        seed={}
        if seed_path:
            seed = json.loads(Path(seed_path).read_text())
            coordinates.update(seed['coordinates'])
            if seed.get('activation_selection_path'):
                selection=json.loads((ROOT/seed['activation_selection_path']).read_text())
                if set(selection['actuators'])!=set(seed['activations']) or any(n not in before['muscles'] for n in selection['actuators']):
                    raise ValueError('Selected activation identities differ from native model/selection receipt')
                activation_names=list(selection['actuators'])
        write('seed_native.json',query(coordinates,seed.get('activations')))
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
        bounds = np.array([effective_bounds[n] for n in names]+[[.01,1.] for _ in activation_names]).T
        if seed_path is None:
            from build_supine_initial_state import passive_neutral_seed
            neutral = passive_neutral_seed(dict(seed_coordinates=coordinates, coordinate_bounds=source_bounds),
                (output/'native/inputs/subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml').read_bytes())
            coordinates.update({n: v for n, v in neutral['seed_coordinates'].items() if n.startswith('elbow_')})
            write('passive_neutral_seed.json', neutral)
        weight = initial['mass_kg']*np.linalg.norm(initial['gravity_m_s2'])
        body_length_m=1.71
        if surface is not None:
            support_manifest=json.loads((ROOT/surface).read_text())
            with np.load(ROOT/support_manifest['arrays_path']) as arrays:
                body_length_m=float(np.ptp(arrays['reference_points_source_m'][:,1]))
        # Shift along the support normal only as far as required for total normal force balance.
        # Subsequent all-coordinate solve retains native root accelerations too.
        preload_trials=[]
        def force_error(tx):
            write('preload_last_request.json',dict(pelvis_tx=float(tx),status='pending'))
            native=query({**coordinates, 'pelvis_tx': tx},seed.get('activations'))
            preload_trials.append(dict(pelvis_tx=float(tx),contact_force_n=native['contact_force_n'],maximum_penetration_m=native['maximum_penetration_m'],surface_foundation=native.get('surface_foundation')))
            write('preload_trials.json',preload_trials)
            write('preload_last_request.json',dict(pelvis_tx=float(tx),status='completed'))
            return native['contact_force_n'][0]-weight
        tx = coordinates['pelvis_tx']
        if surface is not None:
            origin=force_error(tx);previous=tx
            direction=-1 if origin<0 else 1
            lower=upper=tx
            for count in range(1,401):
                trial=tx+direction*.0005*count
                if not source_bounds['pelvis_tx'][0]<=trial<=source_bounds['pelvis_tx'][1]:break
                try:value=force_error(trial)
                except ValueError as error:
                    if 'domain' not in str(error):raise
                    write('preload_domain_rejection.json',dict(error=str(error),rejected_pelvis_tx=trial,retained_original_seed_pelvis_tx=tx,scope='Rejected material-domain trial; refine original admissible seed without preload force balance, no constitutive clamp.'))
                    break
                if value*origin<=0:
                    lower,upper=sorted((previous,trial));break
                previous=trial
        else:
            lower = max(source_bounds['pelvis_tx'][0], tx-.5)
            upper = min(source_bounds['pelvis_tx'][1], tx+.5)
        if lower!=upper and force_error(lower)*force_error(upper) < 0:
            coordinates['pelvis_tx'] = float(brentq(force_error, lower, upper, xtol=1e-12))
        elif surface is None:
            raise ValueError('No bounded vertical support-force root for seed')
        write('preload.json', query(coordinates))
        best = None
        evaluations = 0
        def residual(x):
            nonlocal best, evaluations
            if time.monotonic()-started > (90 if surface is not None else 170):
                raise TimeoutError('Bounded diagnostic wall cap')
            q = {**coordinates, **dict(zip(names, map(float, x[:len(names)])))}
            activations = dict(zip(activation_names, map(float, x[len(names):])))
            native = query(q, activations)
            a = np.asarray(native['udot'])
            # Numerical values correspond to explicit 1 rad/s² and 1 m/s² scales.
            root_indices = [native['mobility_coordinate_names'].index(n) for n in ('pelvis_tilt','pelvis_list','pelvis_rotation')]
            root_translation_indices=[native['mobility_coordinate_names'].index(n) for n in ('pelvis_tx','pelvis_ty','pelvis_tz')]
            tree=np.asarray(native['tree_zero_acceleration_residual_mobility_force'])
            r = np.r_[a, 1000*tree[root_translation_indices]/weight, 1000*tree[root_indices]/(weight*body_length_m)]
            cost = float(r@r)
            evaluations += 1
            if best is None or cost < best['cost']:
                best = dict(cost=cost, coordinates={n:native['coordinates'][n]['value'] for n in q}, requested_coordinates=q, activations=activations, native=native, evaluation=evaluations)
                write('best_candidate.json', best)
            return r
        x0 = np.clip(np.array([coordinates[n] for n in names]+[seed.get('activations',{}).get(n,.01) for n in activation_names]), bounds[0]+1e-10, bounds[1]-1e-10)
        def jacobian(x):
            base=residual(x)
            jac=np.zeros((len(base),len(x)))
            for i in range(len(x)):
                h=1e-5 if i<len(names) else 1e-4
                if x[i]+h>bounds[1,i]:h=-h
                shifted=x.copy();shifted[i]+=h
                jac[:,i]=(residual(shifted)-base)/h
            singular=np.linalg.svd(jac,compute_uv=False)
            write('jacobian_diagnostic.json',dict(evaluation=evaluations,coordinate_step=1e-5,activation_step=1e-4,
                columns=names+activation_names,shape=list(jac.shape),singular_values=singular.tolist(),
                relative_rank_tolerance=1e-8,numerical_rank=int(np.sum(singular>singular[0]*1e-8)),
                condition_number=None if singular[-1]==0 else float(singular[0]/singular[-1]),jacobian=jac.tolist()))
            return jac
        try:
            if surface is None:
                result = least_squares(residual, x0, bounds=bounds, max_nfev=max_nfev,
                    jac=jacobian, x_scale='jac', ftol=1e-10, xtol=1e-10, gtol=1e-8)
            else:
                from types import SimpleNamespace
                x=x0.copy();rejected=[];iterations=[]
                coordinate_radius=.03 if seed.get('support_geometry') else .01
                translation_radius=.001 if seed.get('support_geometry') else .0005
                result=SimpleNamespace(success=False,message='Bounded domain-aware local iteration cap')
                for iteration in range(min(max_nfev,20)):
                    base=residual(x);jac=jacobian(x)
                    radius=np.array([translation_radius if n=='pelvis_tx' else coordinate_radius for n in names]+[.02]*len(activation_names))
                    step=lsq_linear(jac,-base,bounds=(np.maximum(-radius,bounds[0]-x),np.minimum(radius,bounds[1]-x)),tol=1e-9,max_iter=100).x
                    accepted=False
                    for fraction in (1.,.5,.25,.125,.0625,.03125):
                        trial=x+fraction*step
                        try:value=residual(trial)
                        except ValueError as error:
                            if 'domain' not in str(error):raise
                            rejected.append(dict(iteration=iteration,fraction=fraction,error=str(error),values=trial.tolist()))
                            write('rejected_domain_trials.json',rejected)
                            continue
                        if float(value@value)<float(base@base):
                            x=trial;accepted=True;break
                    iterations.append(dict(iteration=iteration,accepted=accepted,prior_cost=float(base@base),best_cost=best['cost'],coordinate_radius=coordinate_radius,translation_radius=translation_radius,accepted_fraction=fraction if accepted else None))
                    if accepted and fraction==1.:
                        coordinate_radius=min(.03,coordinate_radius*1.5);translation_radius=min(.005,translation_radius*1.5)
                    else:
                        coordinate_radius*=.5;translation_radius*=.5
                    write('local_iterations.json',iterations)
                    if not accepted and coordinate_radius<1e-5:
                        result.message='No decreasing admissible local step after radius reductions';break
        except (TimeoutError,ValueError) as error:
            from types import SimpleNamespace
            result = SimpleNamespace(success=False, message=str(error))
        native = best['native']
        after = stream._request('observe')
        unchanged = all(before[k] == after[k] for k in ('time_s','coordinates'))
        artifact = dict(schema='ihm.native-initial-pose.v1', environment='supine',
            target_mass_kg=77.6122029, coordinates=best['coordinates'], requested_coordinates=best['requested_coordinates'], activations=best['activations'],
            activation_selection_path=seed.get('activation_selection_path'),
            protocol_sha256=protocol_sha256, posture_bounds_rad=posture_bounds, body_length_m=body_length_m,
            support_geometry=None if surface is None else dict(manifest_path=str(surface),manifest_sha256=hashlib.sha256((ROOT/surface).read_bytes()).hexdigest(),bed_material=bed),
            accepted_equilibrium=False, native_static_residual=native,
            model_sha256=hashlib.sha256((output/'native/inputs/subject_walk_scaled.osim').read_bytes()).hexdigest(),
            scope='Lower acceleration native initialization with explicit support identity; selected tonic activations are inverse-statics parameters, not learned cortical control; no geometry or force-law change; forward verification required.')
        write('initial_pose.json', artifact)
        report = dict(output=str(output), optimizer_success=bool(result.success), activations=best['activations'],
            objective='All mobility accelerations divided by 1 rad/s² or 1 m/s² plus 1000 times actual tree root forces/(mg) and 1000 times actual tree root moments/(mg*body length); quadratic residual. Pelvis angular acceleration is not used as an external-wrench proxy',
            protocol_sha256=protocol_sha256, posture_bounds_rad=posture_bounds, body_length_m=body_length_m,
            support_geometry=None if surface is None else dict(manifest_path=str(surface),manifest_sha256=hashlib.sha256((ROOT/surface).read_bytes()).hexdigest(),bed_material=bed),
            supported=bool(abs(native['contact_force_n'][0]-weight)/weight < .01),
            optimizer_message=result.message, evaluations=evaluations, wall_s=time.monotonic()-started,
            continuing_state_unchanged=unchanged, accepted_equilibrium=False,
            initial_max_abs_udot=max(map(abs, initial['udot'])),
            final_max_abs_udot=max(map(abs,native['udot'])),
            initial_acceleration_norm=float(np.linalg.norm(initial['udot'])),
            final_acceleration_norm=float(np.linalg.norm(native['udot'])),
            weight_n=weight, final_contact_force_n=native['contact_force_n'],
            final_com_acceleration_m_s2=native['com_acceleration_m_s2'],
            root_tree_generalized_forces={n:native['tree_zero_acceleration_residual_mobility_force'][i] for i,n in enumerate(native['mobility_coordinate_names']) if n.startswith('pelvis_')},
            largest_accelerations=sorted(zip(native['mobility_coordinate_names'], native['udot']), key=lambda x:abs(x[1]),reverse=True)[:10])
        write('report.json', report)
        print(json.dumps(report, indent=2))
    except Exception as error:
        write('failure.json',dict(error_type=type(error).__name__,error=str(error),wall_s=time.monotonic()-started,accepted_equilibrium=False))
        raise
    finally:
        stream.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-nfev', type=int, default=80)
    parser.add_argument('--seed', type=Path, default=ROOT/'data/derived/supine-sphere-pose-7i99jsg8/initial_pose.json')
    parser.add_argument('--surface',type=Path)
    parser.add_argument('--bed',choices=('MM','HM'))
    args = parser.parse_args()
    solve(args.max_nfev, args.seed, args.surface, args.bed)
