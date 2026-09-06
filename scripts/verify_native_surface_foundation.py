"""Explicitly gated two-millisecond controlled native surface-contact fixture."""
from pathlib import Path
import argparse,json,signal,sys,tempfile,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.supine_contact import foundation


def run(manifest_path):
    from ihm.native.mechanical_stream import NativeMechanicalStream
    manifest_path=manifest_path.resolve();manifest=json.loads(manifest_path.read_text())
    data=np.load(ROOT/manifest['arrays_path']);names=manifest['bodies'];owners=data['body_indices']
    output=Path(tempfile.mkdtemp(prefix='native-surface-foundation-',dir=ROOT/'data/derived'));started=time.monotonic();stream=None
    report=dict(passed=False,dt_s=.002,wall_budget_s=60,scope='Controlled contact implementation fixture; no equilibrium or mattress calibration acceptance')
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    def deadline(signum,stack):
        if stream is not None and stream.process.poll() is None:stream.process.kill()
        raise TimeoutError('Surface contact fixture wall cap')
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,60)
    try:
        for variant in ('sphere','surface'):
            stream=NativeMechanicalStream(ROOT,output/variant,environment='supine',target_mass_kg=77.6122029,
                 augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json',
                 surface_contact_manifest=manifest_path if variant=='surface' else None)
            initial=stream.snapshot();write(variant+'_initial.json',initial)
            if variant=='surface':
                assert initial['contact_model']=='retained_skin_foundation'
                assert all(c['geometry_type']=='retained_skin_foundation' for c in initial['contacts'])
                checkpoint=stream.checkpoint()
                evaluated=stream._request('evaluate_static_pose 1 pelvis_tx -0.001')
                reference=np.array(data['reference_points_source_m']);reference[:,0]-=.001
                origins=np.array([initial['bodies'][name]['transform_ground'] for name in names])[:,:3,3];origins[:,0]-=.001
                expected=foundation(reference,np.zeros_like(reference),data['area_m2'],owners,origins,manifest['plane_source_x_m'],manifest['material'])
                assert np.allclose(evaluated['contact_force_n'],expected['body_forces_n'].sum(0),rtol=1e-8,atol=1e-8)
                assert np.isclose(evaluated['surface_foundation']['elastic_energy_j'],expected['elastic_energy_j'],rtol=1e-8,atol=1e-10)
                write('static_translation.json',evaluated)
                stream.restore(checkpoint);stream.release(checkpoint)
            frame=stream.advance(.002);write(variant+'_advanced.json',frame)
            assert np.linalg.norm(frame['momentum_balance_residual_n'])<1e-5
            if variant=='surface':
                points=np.zeros_like(data['stations_local_m']);velocities=np.zeros_like(points);origins=[]
                for index,name in enumerate(names):
                    body=frame['bodies'][name];transform=np.array(body['transform_ground']);selected=owners==index
                    offset=data['stations_local_m'][selected]@transform[:3,:3].T
                    points[selected]=offset+transform[:3,3]
                    velocities[selected]=np.array(body['origin_velocity_m_s'])+np.cross(body['angular_velocity_rad_s'],offset)
                    origins.append(transform[:3,3])
                expected=foundation(points,velocities,data['area_m2'],owners,np.array(origins),manifest['plane_source_x_m'],manifest['material'])
                assert np.allclose(frame['contact_force_n'],expected['body_forces_n'].sum(0),rtol=1e-7,atol=1e-8)
                assert np.allclose(frame['surface_foundation']['bed_moment_about_source_origin_nm'],expected['bed_moment_about_source_origin_nm'],rtol=1e-7,atol=1e-8)
                assert frame['surface_foundation']['dissipative_power_w']<=1e-8
                report['surface_active_points']=expected['contacting_points'];report['surface_support_n']=frame['contact_force_n']
                report['surface_indentation_m']=expected['maximum_penetration_m']
            else:report['sphere_support_n']=frame['contact_force_n']
            stream.close();stream=None
        report['passed']=True
    except Exception as error:report['error']=f'{type(error).__name__}: {error}'
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if stream is not None:stream.close()
        report['wall_s']=time.monotonic()-started;write('report.json',report);print(json.dumps({**report,'output':str(output)},indent=2))
    return 0 if report['passed'] else 1


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-native',action='store_true');parser.add_argument('--manifest',type=Path);args=parser.parse_args()
    if not args.run_native or args.manifest is None:raise SystemExit('Coordinated --run-native slot and --manifest required')
    raise SystemExit(run(args.manifest))
