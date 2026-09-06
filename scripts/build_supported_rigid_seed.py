"""Bounded source-only force/moment-balanced rigid pose seed on measured bed.

This does not solve native joint/muscle residuals or install a resting reference.
"""
from pathlib import Path
import argparse,hashlib,json,signal,sys,tempfile,time
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.bed_compression import load_bed,series_response,maximum_approach


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--material',choices=('MM','HM'),default='MM');args=parser.parse_args()
    manifest_path=ROOT/'data/derived/supine-surface-contact-exmzq9pq/manifest.json';manifest=json.loads(manifest_path.read_text())
    data=np.load(ROOT/manifest['arrays_path']);reference_path=ROOT/'data/derived/supine-support-5ma720yd/initial_native.json';native=json.loads(reference_path.read_text())
    points=data['reference_points_source_m'];areas=data['area_m2'];bed=load_bed(ROOT,args.material);skin=manifest['material']
    center=np.zeros(3);mass=0.
    for body in native['bodies'].values():
        transform=np.array(body['transform_ground']);center+=body['mass_kg']*(transform[:3,:3]@body['mass_center_local_m']+transform[:3,3]);mass+=body['mass_kg']
    center/=mass;weight=mass*9.81;length=float(np.ptp(points[:,1]));limit=maximum_approach(skin,bed)
    output=Path(tempfile.mkdtemp(prefix='supported-rigid-seed-',dir=ROOT/'data/derived'));started=time.monotonic();evaluations=[];best=None
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    def deadline(signum,stack):raise TimeoutError('Source rigid-seed search reached60s wall cap')
    def evaluate(values):
        nonlocal best
        if len(evaluations)>=200:raise RuntimeError('Source rigid-seed search reached200 actual evaluations')
        fraction,pitch,roll=values;c,s=np.cos(pitch),np.sin(pitch);a,b=np.cos(roll),np.sin(roll)
        rotation=np.array([[c,-s,0],[s,c,0],[0,0,1]])@np.array([[a,0,b],[0,1,0],[-b,0,a]])
        world=(points-center)@rotation.T+center
        translation=world[:,0].min()-manifest['plane_source_x_m']+fraction*limit
        world[:,0]-=translation;approach=np.maximum(0,manifest['plane_source_x_m']-world[:,0])
        response=series_response(approach,skin,bed);force=areas*response['pressure_pa'];total=float(force.sum())
        displaced_center=center+[-translation,0,0]
        moment=np.cross(world-displaced_center,np.column_stack((force,np.zeros_like(force),np.zeros_like(force)))).sum(0)
        residual=np.array([(total-weight)/weight,moment[1]/(weight*length),moment[2]/(weight*length)])
        score=float(np.linalg.norm(residual));row=dict(evaluation=len(evaluations)+1,wall_s=time.monotonic()-started,residual_norm=score,
              support_n=total,moment_about_com_nm=moment.tolist(),pitch_rad=float(pitch),roll_rad=float(roll),downward_translation_m=float(translation))
        evaluations.append(row)
        if best is None or score<best['residual_norm']:
            root_origin=np.array(native['bodies']['pelvis']['transform_ground'])[:3,3]
            shifted_root=rotation@(root_origin-center)+center+[-translation,0,0]
            coordinates={name:c['value'] for name,c in native['coordinates'].items() if not name.endswith('_beta')}
            coordinates.update(pelvis_tilt=float(pitch),pelvis_rotation=float(roll),pelvis_list=0.,pelvis_tx=float(shifted_root[0]),pelvis_ty=float(shifted_root[1]),pelvis_tz=float(shifted_root[2]))
            best={**row,'seed_coordinates':coordinates,'maximum_skin_indentation_m':float(response['skin_indentation_m'].max()),
                 'maximum_bed_indentation_m':float(response['bed_indentation_m'].max()),'contacting_points':int(np.sum(force>0)),
                 'source_com_m':center.tolist(),'source_to_candidate_rotation':rotation.tolist()}
        return residual
    report=dict(passed=False,native_run=False,accepted_equilibrium=False,material=args.material,maximum_evaluations=200,maximum_wall_s=60,
                scope='Net force/two support moments of whole-body rigid seed only; native joint and muscle equilibrium and forward support remain unverified',
                search_domain='Global pitch/long-axis roll limited to±0.2rad to retain posterior-envelope scope; source pelvis bounds are wider')
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,60)
    try:
        result=least_squares(evaluate,[.7,0.,0.],bounds=([0,-.2,-.2],[1-1e-9,.2,.2]),max_nfev=200,ftol=1e-10,xtol=1e-10,gtol=1e-10)
        report.update(passed=best['residual_norm']<1e-7,status='rigid_support_constraints_solved' if best['residual_norm']<1e-7 else 'unresolved',optimizer_message=str(result.message))
    except Exception as error:report.update(status='stopped',error=f'{type(error).__name__}: {error}')
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        report.update(wall_s=time.monotonic()-started,evaluations=len(evaluations),best=best,
          source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (manifest_path,reference_path,Path(__file__).resolve(),ROOT/'ihm/assembly/bed_compression.py')})
        write('report.json',report);write('evaluations.json',evaluations);print(json.dumps({**report,'output':str(output)},indent=2))


if __name__=='__main__':main()
