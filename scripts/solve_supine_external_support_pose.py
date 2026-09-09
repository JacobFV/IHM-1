"""Balance actual external root force/moment using only three rigid pose coordinates."""
from pathlib import Path
import json,sys,tempfile,time,argparse,hashlib
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import lsq_linear
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream

def solve(seed_path,seconds=80):
    seed_path=seed_path.resolve();seed=json.loads(seed_path.read_text());q=seed['coordinates'];activations=seed.get('activations',{});support=seed['support_geometry']
    output=Path(tempfile.mkdtemp(prefix='supine-external-support-',dir=ROOT/'data/derived'));started=time.monotonic();protocol=Path(__file__).read_bytes();(output/'solver.py').write_bytes(protocol)
    def write(name,data):(output/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    manifest=json.loads((ROOT/support['manifest_path']).read_text())
    with np.load(ROOT/manifest['arrays_path']) as arrays:length=float(np.ptp(arrays['reference_points_source_m'][:,1]))
    stream=NativeMechanicalStream(ROOT,output/'native',environment='supine',target_mass_kg=77.6122029,augmented_registration='data/derived/mechanics/whole_body_lumbar_current/registration.json',surface_contact_manifest=support['manifest_path'],bed_material=support['bed_material'])
    names=['pelvis_tx','pelvis_tilt','pelvis_rotation'];ranges={c.get('name'):list(map(float,c.findtext('range').split())) for c in ET.parse(output/'native/inputs/subject_walk_scaled.osim').findall('.//Coordinate')}
    bounds=np.array([ranges[n] if n=='pelvis_tx' else [max(-.15,ranges[n][0]),min(.15,ranges[n][1])] for n in names]).T
    records=[];rejected=[];best=None
    def evaluate(x):
        nonlocal best
        if time.monotonic()-started>seconds:raise TimeoutError('Rigid external-support wall cap')
        candidate={**q,**dict(zip(names,map(float,x)))};args=['evaluate_static_pose',str(len(candidate))]
        for n,v in candidate.items():args.extend([n,str(v)])
        args.append(str(len(activations)))
        for n,v in activations.items():args.extend([n,str(v)])
        native=stream._request(' '.join(args));mapping=dict(zip(native['mobility_coordinate_names'],range(len(native['udot']))));force=np.asarray(native['tree_zero_acceleration_residual_mobility_force']);weight=native['mass_kg']*np.linalg.norm(native['gravity_m_s2'])
        r=np.array([force[mapping[n]]/(weight if n=='pelvis_tx' else weight*length) for n in names]);all_names=['pelvis_tx','pelvis_ty','pelvis_tz','pelvis_tilt','pelvis_list','pelvis_rotation'];all_r=[force[mapping[n]]/(weight if n in all_names[:3] else weight*length) for n in all_names]
        entry=dict(coordinates={n:native['coordinates'][n]['value'] for n in candidate},requested_coordinates=candidate,max_assembly_q_adjustment=max(abs(native['coordinates'][n]['value']-v) for n,v in candidate.items()),activations=activations,native=native,cost=float(r@r),support_residual=r.tolist(),all_root_residual=all_r,root_coordinate_order=all_names,weight_n=float(weight),body_length_m=length)
        records.append(dict(evaluation=len(records)+1,values=x.tolist(),residual=r.tolist(),all_root_residual=all_r))
        write('evaluations.json',records)
        if best is None or entry['cost']<best['cost']:best=entry;write('best_candidate.json',best)
        return entry
    try:
        before=stream.snapshot();x=np.clip(np.array([q[n] for n in names]),bounds[0]+1e-11,bounds[1]-1e-11);first=evaluate(x);write('initial.json',first)
        radius=np.array([.005,.03,.03]);iterations=[];message='Iteration cap'
        for iteration in range(25):
            entry=evaluate(x);r=np.asarray(entry['support_residual'])
            if max(map(abs,entry['all_root_residual']))<1e-7:message='All external root force/moment residuals converged';break
            j=np.zeros((3,3))
            for i in range(3):
                h=1e-5 if x[i]+1e-5<=bounds[1,i] else -1e-5;y=x.copy();y[i]+=h
                column=None
                for multiplier in (1.,-1.,.1,-.1,.01,-.01,.001,-.001):
                    delta=h*multiplier;y=x.copy();y[i]+=delta
                    if not bounds[0,i]<=y[i]<=bounds[1,i]:continue
                    try:column=(np.asarray(evaluate(y)['support_residual'])-r)/delta;break
                    except ValueError as error:
                        if 'domain' not in str(error):raise
                        rejected.append(dict(iteration=iteration,kind='derivative',coordinate=names[i],delta=delta,error=str(error)));write('rejected_domain_trials.json',rejected)
                if column is None:raise ValueError('No admissible coordinate derivative inside retained material domain')
                j[:,i]=column
            step=lsq_linear(j,-r,bounds=(np.maximum(-radius,bounds[0]-x),np.minimum(radius,bounds[1]-x)),tol=1e-12,max_iter=100).x
            accepted=False
            for fraction in (1.,.5,.25,.125,.0625,.03125):
                trial=x+fraction*step
                try:value=evaluate(trial)
                except ValueError as error:
                    if 'domain' not in str(error):raise
                    rejected.append(dict(iteration=iteration,fraction=fraction,values=trial.tolist(),error=str(error)));write('rejected_domain_trials.json',rejected);continue
                if value['cost']<entry['cost']:x=trial;accepted=True;break
            iterations.append(dict(iteration=iteration,accepted=accepted,fraction=fraction,cost=best['cost'],radius=radius.tolist()));write('iterations.json',iterations)
            if accepted and fraction==1:radius=np.minimum(radius*1.5,[.01,.03,.03])
            else:radius*=.5
            if not accepted and max(radius)<1e-5:message='No decreasing admissible rigid step';break
    except (TimeoutError,ValueError) as error:message=str(error)
    finally:
        if best is not None:
            observed=stream._request('observe');stream.close();artifact={**seed,'coordinates':best['coordinates'],'native_static_residual':best['native'],'external_support_residual':best['all_root_residual'],'external_support_order':best['root_coordinate_order'],'body_length_m':length,'protocol_sha256':hashlib.sha256(protocol).hexdigest(),'accepted_equilibrium':False,'scope':'Three rigid root coordinates adjusted to true external force/moment balance; all internal q and tones fixed; joint equilibrium and forward stability remain separate.'}
            write('initial_pose.json',artifact)
            report=dict(output=str(output),message=message,wall_s=time.monotonic()-started,evaluations=len(records),domain_rejections=len(rejected),body_length_m=length,initial_root_residual=first['all_root_residual'],final_root_residual=best['all_root_residual'],root_order=best['root_coordinate_order'],external_wrench_balanced=bool(max(map(abs,best['all_root_residual']))<1e-7),accepted_equilibrium=False,max_abs_udot=max(map(abs,best['native']['udot'])),udot_norm=float(np.linalg.norm(best['native']['udot'])),contact_force_n=best['native']['contact_force_n'],coordinates={n:best['coordinates'][n] for n in names},continuing_state_unchanged=all(before[k]==observed[k] for k in ('time_s','coordinates','muscles')))
            write('report.json',report);print(json.dumps(report,indent=2))
        stream.close()
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--seed',type=Path,default=ROOT/'data/derived/supine-hip-tones-w3_krvzv/initial_pose.json');p.add_argument('--seconds',type=float,default=80);a=p.parse_args();solve(a.seed,a.seconds)
