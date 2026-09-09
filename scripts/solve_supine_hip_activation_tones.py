"""Fixed-pose native hip/lumbar inverse statics; no learned brain controller."""
from pathlib import Path
import json,sys,tempfile,time,hashlib,argparse
from types import SimpleNamespace
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream

def solve(seed_path,seconds=100):
    seed_path=seed_path.resolve();seed=json.loads(seed_path.read_text());q=seed['coordinates'];support=seed['support_geometry']
    output=Path(tempfile.mkdtemp(prefix='supine-hip-tones-',dir=ROOT/'data/derived'));started=time.monotonic();protocol=Path(__file__).read_bytes();(output/'solver.py').write_bytes(protocol)
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    stream=NativeMechanicalStream(ROOT,output/'native',environment='supine',target_mass_kg=77.6122029,augmented_registration='data/derived/mechanics/whole_body_lumbar_current/registration.json',initial_pose=q,surface_contact_manifest=support['manifest_path'],bed_material=support['bed_material'])
    try:
        initial=stream.snapshot();write('initialized.json',initial)
        hip_coordinates=[f'hip_{axis}_{side}' for side in ('r','l') for axis in ('rotation','adduction','flexion')]
        arms=stream.moment_arms(muscles=list(initial['muscles']),coordinates=hip_coordinates);write('native_moment_arms.json',arms)
        selection=[];selected=set(seed['activations'])
        for coordinate in hip_coordinates:
            for sign in (-1,1):
                candidates=[]
                for muscle,values in arms['moment_arms_m'].items():
                    arm=values[coordinate]
                    if sign*arm>.01:
                        maximum_force=initial['muscles'][muscle]['max_isometric_force_n']
                        candidates.append(dict(muscle=muscle,moment_arm_m=arm,max_isometric_force_n=maximum_force,maximum_isometric_torque_nm=maximum_force*abs(arm)))
                candidates.sort(key=lambda item:item['maximum_isometric_torque_nm'],reverse=True)
                chosen=candidates[:2];selected.update(item['muscle'] for item in chosen)
                selection.append(dict(coordinate=coordinate,sign=sign,candidates=candidates,selected=chosen))
        names=sorted(selected);write('selection.json',dict(criteria='At actual fixed candidate pose, abs native moment arm >0.01m; top2 Fmax*abs(moment arm) per sign per hip axis, both legs, plus existing six lumbar actuators.',actuators=names,groups=selection))
        evaluations=0;best=None
        def evaluate(x):
            nonlocal best,evaluations
            if time.monotonic()-started>seconds:raise TimeoutError('Fixed-pose activation diagnostic wall cap')
            activations=dict(zip(names,map(float,x)));args=['evaluate_static_pose',str(len(q))]
            for n,v in q.items():args.extend([n,str(v)])
            args.append(str(len(activations)))
            for n,v in activations.items():args.extend([n,str(v)])
            native=stream._request(' '.join(args));a=np.asarray(native['udot']);root=[native['mobility_coordinate_names'].index(n) for n in ('pelvis_tilt','pelvis_list','pelvis_rotation')]
            residual=np.r_[a,100*np.asarray(native['com_acceleration_m_s2']),20*a[root]];cost=float(residual@residual);evaluations+=1
            if best is None or cost<best['cost']:
                best=dict(coordinates=q,activations=activations,cost=cost,native=native,evaluation=evaluations);write('best_candidate.json',best)
            return residual
        def jacobian(x):
            r=evaluate(x);j=np.zeros((len(r),len(x)))
            for i in range(len(x)):
                h=1e-4 if x[i]+1e-4<=1 else -1e-4;y=x.copy();y[i]+=h;j[:,i]=(evaluate(y)-r)/h
            s=np.linalg.svd(j,compute_uv=False);write('jacobian_diagnostic.json',dict(columns=names,singular_values=s.tolist(),rank_relative1e_8=int(sum(s>s[0]*1e-8)),step=1e-4))
            return j
        x0=np.array([seed['activations'].get(n,initial['muscles'][n]['activation']) for n in names]);x0=np.clip(x0,.0100000001,.9999999999)
        base=evaluate(x0);baseline=json.loads(json.dumps(best));write('baseline.json',baseline)
        try:result=least_squares(evaluate,x0,jac=jacobian,bounds=(.01,1.),x_scale='jac',max_nfev=40,ftol=1e-10,gtol=1e-8,xtol=1e-10)
        except TimeoutError as error:result=SimpleNamespace(success=False,message=str(error))
        native=best['native'];observed=stream._request('observe')
        artifact={**seed,'coordinates':q,'activations':best['activations'],'native_static_residual':native,'protocol_sha256':hashlib.sha256(protocol).hexdigest(),'activation_selection_path':str((output/'selection.json').relative_to(ROOT)),'accepted_equilibrium':False,'scope':'Fixed-pose native selected hip/lumbar activation inverse statics; tonic candidate, not learned cortical control; geometry unchanged.'}
        write('initial_pose.json',artifact)
        report=dict(output=str(output),optimizer_success=bool(result.success),optimizer_message=str(result.message),wall_s=time.monotonic()-started,evaluations=evaluations,selected_count=len(names),actuators=names,
            continuing_state_unchanged=all(initial[k]==observed[k] for k in ('coordinates','time_s','muscles')),coordinates_held_fixed=True,
            baseline_max_abs_udot=max(map(abs,baseline['native']['udot'])),final_max_abs_udot=max(map(abs,native['udot'])),baseline_udot_norm=float(np.linalg.norm(baseline['native']['udot'])),final_udot_norm=float(np.linalg.norm(native['udot'])),baseline_cost=baseline['cost'],final_cost=best['cost'],contact_force_n=native['contact_force_n'],com_acceleration_m_s2=native['com_acceleration_m_s2'],accepted_equilibrium=False,
            largest_accelerations=sorted(zip(native['mobility_coordinate_names'],native['udot']),key=lambda item:abs(item[1]),reverse=True)[:12],activations=best['activations'])
        write('report.json',report);print(json.dumps(report,indent=2))
    finally:stream.close()
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--seed',type=Path,default=ROOT/'data/derived/supine-lumbar-activation-caxd5_y_/initial_pose.json');parser.add_argument('--seconds',type=float,default=100);args=parser.parse_args();solve(args.seed,args.seconds)
