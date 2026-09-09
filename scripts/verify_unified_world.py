#!/usr/bin/env python3
"""Retained sustained native brain/body/world runs and matched kernel ablations.

Survival is checked separately from cortical contribution. No walking success is
inferred from a moving centroid or a nonzero excitation.
"""
import argparse
import cProfile
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time

import numpy as np

SCRIPT_BYTES = Path(__file__).read_bytes()
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.embodied import EmbodiedRuntime


def run(root, output, seconds, environment, scene, kind, sever=False, no_cord=False, target_rad=.12, profile=False, objects=(), pelvis_push=None):
    body = None
    started = time.monotonic()
    records = []
    profiler=cProfile.Profile() if profile else None
    result = {'controller': {'kind': kind, 'sever': sever, 'no_cord': no_cord},
              'environment': environment, 'requested_seconds': seconds,
              'passed': False, 'scope': 'Sustained native integration; not walking validation'}
    if kind in ('implicit_ankle_primitive','implicit_cortical_ankle'):
        result['controller']['target_rad'] = target_rad
    output.mkdir()
    try:
        body = EmbodiedRuntime.from_workspace(root, output / 'runtime', environment=environment,
                environment_selection={'scene': scene, 'objects': list(objects)} if scene else None,
                controller=result['controller'])
        initial = body.snapshot()
        owner = body.environment_dynamics
        from ihm.assembly.controller_selection import STANCE_KINDS
        stance = kind in STANCE_KINDS
        if stance:
            from ihm.native.moment_arm_control import native_center_of_mass
            initial_native=body.plant.native.snapshot()
            initial_com,_=native_center_of_mass(initial_native)
            result['pelvis_push_source_n']=pelvis_push
        previous_skin = None if owner is None else owner.skin_points(initial['mechanics']['entities'])
        if owner is not None:
            assert owner.native_frame_bound
            np.testing.assert_allclose(owner.world_frame.canonical_gravity_m_s2,
                body.plant.registration.basis @ np.asarray(body.plant.registration.reference['gravity_m_s2']), atol=1e-12)
            result['world_frame'] = owner.world_frame.metadata()
            result['world_diagnostic_scope'] = 'State sampled each 20 ms; constraint work is last 1 ms substep only, not cumulative work'

        print(json.dumps({'event':'initialized','output':str(output)}),flush=True)
        if profiler:profiler.enable()
        for _ in range(round(seconds / .02)):
            inputs={'seconds': .02}
            if pelvis_push is not None and 1.-1e-9 <= body.time_s < 1.1-1e-9:
                native=body.plant.native.snapshot();registration=body.plant.registration
                point=body.plant.native.body_point(body='pelvis',station_m=native['bodies']['pelvis']['mass_center_local_m'])['point_source_m']
                force=np.array([pelvis_push[0],0.,pelvis_push[1]])
                ident=registration.groups['pelvis']['canonical_bones'][0]
                port={'id':ident,'point_m':(registration.basis@point+registration.maps['pelvis'][:3,3]).tolist(),
                      'force_n':(registration.basis@force).tolist()}
                mapped=registration.force(ident,port['point_m'],port['force_n'],native)
                assert mapped['body']=='pelvis'
                np.testing.assert_allclose(mapped['point_m'],point,atol=1e-12)
                np.testing.assert_allclose(mapped['force_n'],force,atol=1e-12)
                inputs['forces']=[port]
            frame = body.step(inputs)
            t = frame['time_s']
            assert abs(frame['mechanics']['time_s'] - t) < 1e-8
            world = frame.get('environment_state')
            if world: assert abs(world['time_s'] - t) < 1e-8
            neural = frame['neural']
            assert abs(neural['time_s'] - t) < 1e-8
            commands = neural['motor_excitations']
            assert commands and all(np.isfinite(v) and 0 <= v <= 1 for v in commands.values())
            row = {'time_s': t, 'motor_excitations': commands, 'arc_max': neural.get('arc_max', {}),
                   'contact_count': (world or {}).get('contact_count', 0),
                   'controller': neural.get('controller'),
                   'joints': frame['mechanics'].get('joints', {}),
                   'motor_primitive': neural.get('motor_primitive'),
                   'lqr_stance': neural.get('lqr_stance'),
                   'cortical_stance': neural.get('cortical_stance'),
                   'physiology_elapsed_s': frame['physiology'].get('elapsed_s')}
            if owner is not None:
                skin = owner.skin_points(frame['mechanics']['entities'])
                row['skin_speed_over_tick_max_m_s'] = float(np.linalg.norm(skin-previous_skin,axis=1).max()/.02)
                previous_skin = skin
                row['world_objects'] = {mesh.id: {
                    'max_speed_m_s':float(np.linalg.norm(mesh.v,axis=1).max()),
                    'material_stretch':mesh.stretch_state,
                    'constraint_last_substep':mesh.constraint_state,
                    'bending':mesh.bending_state,
                    'triangle_self_contact':mesh.triangle_contact_state,
                } for mesh in owner.soft}
                row['rigid_objects'] = {prop.id:{'position_m':prop.x.tolist(),'velocity_m_s':prop.v.tolist(),
                    'contact_dissipation_j':prop.contact_dissipation_j} for prop in owner.rigid}
            if stance:
                native=body.plant.native.snapshot();com,velocity=native_center_of_mass(native)
                row['com_horizontal_displacement_m']=float(np.linalg.norm((com-initial_com)[[0,2]]))
                row['com_speed_m_s']=float(np.linalg.norm(velocity))
                row['pelvis_height_m']=native['coordinates']['pelvis_ty']['value']
                row['fallen']=row['pelvis_height_m']<.6
            records.append(row)
            if stance and row['fallen']:break
            if len(records)%50==0:print(json.dumps({'event':'progress','time_s':t,'output':str(output)}),flush=True)
        result.update(passed=True, duration_s=frame['time_s'], steps=len(records))
        if stance:
            result['stance']={'displacement_basis':'horizontal source-ground x/z; velocity is full3D', 'fallen':records[-1]['fallen'],
                'peak_com_displacement_m':max(r['com_horizontal_displacement_m'] for r in records),
                'final_com_displacement_m':records[-1]['com_horizontal_displacement_m'],
                'final_com_speed_m_s':records[-1]['com_speed_m_s'],
                'completed_requested_duration':len(records)==round(seconds/.02)}
    except Exception as error:
        result.update(error=f'{type(error).__name__}: {error}', steps=len(records),
                      duration_s=records[-1]['time_s'] if records else 0.)
    finally:
        if profiler:
            profiler.disable();profiler.dump_stats(str(output/'steps.prof'))
        if body is not None:
            try: body.close()
            except Exception as error: result.update(passed=False, cleanup_error=str(error))
        result['wall_s'] = time.monotonic() - started
        (output / 'trace.json').write_text(json.dumps(records, indent=2, allow_nan=False) + '\n')
        (output / 'report.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result, records


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--seconds', type=float, default=5.)
    ap.add_argument('--environment', choices=['supine', 'upright', 'free'], default='supine')
    ap.add_argument('--scene', default='bedroom')
    ap.add_argument('--arms', nargs='+', choices=['regional', 'full', 'sever', 'no-cord'],
                    default=['full', 'sever', 'no-cord'])
    ap.add_argument('--objects', nargs='*', default=[], help='Additional catalogue objects in both matched arms')
    ap.add_argument('--output', type=Path)
    ap.add_argument('--profile',action='store_true',help='Profile accepted-step work, excluding native initialization')
    ap.add_argument('--controller', choices=['implicit', 'implicit_curriculum16', 'implicit_ankle_primitive', 'implicit_cortical_ankle','engineering_stance','implicit_cortical_stance','implicit_curriculum16_stance'], default='implicit')
    ap.add_argument('--pelvis-push',nargs=2,type=float,metavar=('X_N','Z_N'),help='100 ms native-ground COM push starting at 1 second; stance modes only')
    ap.add_argument('--target-rad', type=float, default=.12)
    args = ap.parse_args()
    if args.controller=='engineering_stance' and (args.environment!='upright' or args.arms!=['full']):
        ap.error('Engineered stance requires --environment upright --arms full')
    from ihm.assembly.controller_selection import STANCE_KINDS,CORTICAL_STANCE_KINDS
    if args.pelvis_push is not None and (args.controller not in STANCE_KINDS or args.environment!='upright' or not all(np.isfinite(x) and abs(x)<=100 for x in args.pelvis_push)):
        ap.error('Pelvis push requires upright stance controller and finite components within 100 N')
    if not np.isfinite(args.seconds) or not .02 <= args.seconds <= 120 or abs(args.seconds / .02 - round(args.seconds / .02)) > 1e-8:
        ap.error('--seconds must be a multiple of 20 ms in [.02,120]')
    output = args.output or Path(tempfile.mkdtemp(prefix='unified-world-', dir=ROOT/'data/derived'))
    output.mkdir(exist_ok=True, parents=True)
    (output/'executed_verifier.py').write_bytes(SCRIPT_BYTES)
    reports, traces = {}, {}
    for arm in args.arms:
        reports[arm], traces[arm] = run(ROOT, output/arm, args.seconds, args.environment,
              args.scene or None, 'regional' if arm == 'regional' else args.controller,
              sever=arm == 'sever', no_cord=arm == 'no-cord', target_rad=args.target_rad,profile=args.profile,objects=args.objects,pelvis_push=args.pelvis_push)
        print(json.dumps({'arm': arm, **reports[arm]}), flush=True)
    comparison = {'cortical_control_demonstrated': False,
                  'reason': 'Trajectory difference alone is not task performance'}
    if traces.get('full') and traces.get('sever'):
        count = min(len(traces['full']), len(traces['sever']))
        ids = sorted(traces['full'][0]['motor_excitations'])
        arrays = {arm: np.array([[row['motor_excitations'][m] for m in ids]
                                for row in traces[arm][:count]]) for arm in ('full', 'sever')}
        comparison['matched_steps'] = count
        comparison['command_rms_difference'] = float(np.sqrt(np.mean((arrays['full'] - arrays['sever'])**2)))
        comparison['temporal_sd_of_command_mean'] = {arm: float(value.mean(axis=1).std()) for arm, value in arrays.items()}
    if args.controller in ('implicit_ankle_primitive','implicit_cortical_ankle'):
        comparison['ankle_target_mse_rad2'] = {arm: float(np.mean([(row['motor_primitive']['angle_rad'] - row['motor_primitive']['target_rad'])**2 for row in rows])) for arm, rows in traces.items() if rows and all(row['motor_primitive'] for row in rows)}
        comparison['motor_owner'] = ('Trained eight-site materialization; larger cortex observed with inactive motor output' if args.controller == 'implicit_ankle_primitive' else 'Trained 128-site IBM cortical E/I dynamics')
        mse = comparison['ankle_target_mse_rad2']
        if all(arm in mse and reports[arm]['passed'] for arm in ('full','sever')):
            improved = mse['full'] < mse['sever'] - 1e-8
            comparison['task_improvement_demonstrated'] = improved
            comparison['relative_mse_reduction'] = 1 - mse['full']/mse['sever'] if mse['sever'] else None
            comparison['cortical_control_demonstrated'] = improved and args.controller == 'implicit_cortical_ankle'
            comparison['reason'] = 'Ankle tracking only, for this target/initial state and duration; not walking or general motor competence'
    if args.controller in CORTICAL_STANCE_KINDS:
        comparison['stance']={arm:r.get('stance') for arm,r in reports.items()}
        comparison['motor_owner']='Persistent 1024-site trained IBM E/I cortex; privileged native state feedback'
        comparison['kernel_lineage']=next((rows[0]['controller'].get('kernel_identity') for rows in traces.values() if rows), None)
        if args.pelvis_push is not None and all(reports.get(arm,{}).get('passed') and reports[arm].get('stance') for arm in ('full','sever')):
            full=reports['full']['stance'];sever=reports['sever']['stance']
            comparison['cortical_control_demonstrated']=bool(full['completed_requested_duration'] and not full['fallen'] and
                full['final_com_displacement_m']<.01 and full['final_com_speed_m_s']<.001 and sever['fallen'])
            comparison['reason']='Bounded push recovery on this fixed body/model/physiology only; no walking or general brain competence claim'
    report = {'passed': all(r['passed'] for r in reports.values()), 'arms': reports,
              'comparison': comparison, 'script_sha256': hashlib.sha256(SCRIPT_BYTES).hexdigest()}
    (output/'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'output': str(output), **report}, indent=2))
    if not report['passed']: raise SystemExit(1)


if __name__ == '__main__':
    main()
