"""Held-tone native body + reciprocal bedroom contact acceptance diagnostic."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import tempfile
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.assembly.articulated import ArticulatedBodyPlant
from ihm.assembly.environment_dynamics import EnvironmentDynamics,prepare_cloth


def run(registration,pose_path,duration=.2,surface=None,bed=None):
    output=Path(tempfile.mkdtemp(prefix='supine-tonic-world-',dir=ROOT/'data/derived'))
    (output/'verification_source.py').write_bytes(Path(__file__).read_bytes())
    (output/'environment_dynamics_source.py').write_bytes((ROOT/'ihm/assembly/environment_dynamics.py').read_bytes())
    artifact=json.loads(pose_path.read_text());started=time.monotonic();rows=[];body=None
    report={'passed':False,'registration':registration,'pose_path':str(pose_path),
            'pose_sha256':hashlib.sha256(pose_path.read_bytes()).hexdigest(),'duration_requested_s':duration,
            'surface_contact_manifest':surface,'bed_material':bed,
            'scope':'Real finite-mass native body and reciprocal bedroom contact with held native muscle excitation defaults; finite-state and momentum diagnostic only, no energy-passivity, brain, standing, or equilibrium claim'}
    try:
        body=ArticulatedBodyPlant(ROOT,output/'body',environment='supine',target_mass_kg=artifact['target_mass_kg'],
            augmented_registration=registration,initial_pose=artifact['coordinates'],surface_contact_manifest=surface,bed_material=bed)
        frame=body.snapshot();environment=EnvironmentDynamics(ROOT,'supine',{'scene':'bedroom','objects':[]},body.registration)
        points=environment.skin_points(frame['entities'])
        for mesh in environment.soft:
            if mesh.kind=='cloth':prepare_cloth(mesh,points,environment.base_plane)
        initial_cloth={mesh.id:mesh.x.copy() for mesh in environment.soft}
        initial_native=body.native.snapshot()
        report['initial_native_coordinates']=initial_native['coordinates']
        report['initial_excitations']={name:m['excitation'] for name,m in initial_native['muscles'].items()}
        report['initial_native_energy']={name:initial_native[name] for name in ('kinetic_energy_j','potential_energy_j','metabolic_reference')}
        for step in range(round(duration/.005)):
            ports=environment.advance(.005,frame['entities'])
            frame=body.advance(.005,ports)
            native=body.native.snapshot()
            points_next=environment.skin_points(frame['entities'])
            rows.append({'time_s':frame['time_s'],'environment_time_s':environment.time_s,
                         'maximum_skin_step_m':float(np.max(np.linalg.norm(points_next-points,axis=1))),
                         'maximum_cloth_speed_m_s':max(float(np.max(np.linalg.norm(mesh.v,axis=1))) for mesh in environment.soft),
                         'maximum_cloth_displacement_m':max(float(np.max(np.linalg.norm(mesh.x-initial_cloth[mesh.id],axis=1))) for mesh in environment.soft),
                         'kinetic_energy_j':native['kinetic_energy_j'],
                         'momentum_residual_n':float(np.linalg.norm(native['momentum_balance_residual_n'])),
                         'contact_ports':len(ports)})
            points=points_next
            assert abs(frame['time_s']-environment.time_s)<1e-9
            assert rows[-1]['momentum_residual_n']<1e-5
            if step%10==0:print(json.dumps(rows[-1]),flush=True)
        report['final_native_coordinates']=native['coordinates']
        report['final_native_surface_foundation']=native.get('surface_foundation')
        report['final_absolute_coordinate_displacement']={name:abs(value['value']-initial_native['coordinates'][name]['value']) for name,value in native['coordinates'].items()}
        report['passed']=True
    except Exception as error:
        report['error']=type(error).__name__+': '+str(error)
    finally:
        if body is not None:body.close()
        report.update(rows=rows,wall_s=time.monotonic()-started)
        (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:v for k,v in report.items() if k not in ('rows','initial_excitations','initial_native_coordinates','final_native_coordinates','final_native_surface_foundation','final_absolute_coordinate_displacement')}|{'output':str(output)},indent=2))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('registration');parser.add_argument('pose',type=Path)
    parser.add_argument('--duration',type=float,default=.2);parser.add_argument('--surface');parser.add_argument('--bed')
    args=parser.parse_args()
    if args.duration<=0 or abs(args.duration/.005-round(args.duration/.005))>1e-8:parser.error('duration must be positive multiple of .005')
    report=run(args.registration,args.pose,args.duration,args.surface,args.bed)
    if not report['passed']:raise SystemExit(1)
