"""Refine joint statics while preserving actual external root wrench balance."""
from pathlib import Path
import argparse,hashlib,json,sys,tempfile,time
import xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream
from bounded_static_root import constrained_local_step


def run(seed_path,seconds=90,distal_only=False,joint_block=None,objective="torque",activation_only=False,retain_selection=False):
    block=joint_block.split(",") if joint_block else (["mtp_angle_l","mtp_angle_r","ankle_angle_l","ankle_angle_r"] if distal_only else [])
    distal_only=bool(block)
    seed_path=seed_path.resolve();seed=json.loads(seed_path.read_text());q=dict(sorted(seed['coordinates'].items()));support=seed['support_geometry']
    output=Path(tempfile.mkdtemp(prefix='supine-supported-tones-',dir=ROOT/'data/derived'));started=time.monotonic();protocol=Path(__file__).read_bytes();(output/'solver.py').write_bytes(protocol)
    def write(name,data):(output/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    stream=NativeMechanicalStream(ROOT,output/'native',environment='supine',target_mass_kg=77.6122029,augmented_registration='data/derived/mechanics/whole_body_lumbar_current/registration.json',initial_pose=q,surface_contact_manifest=support['manifest_path'],bed_material=support['bed_material'])
    best=None;records=[];rejections=[];message='Iteration cap';first=None
    try:
        before=stream.snapshot();write('initialized.json',before)
        names=[n for n in q if n not in ('pelvis_ty','pelvis_tz','pelvis_list')]
        if distal_only:names=[n for n in names if n in ['pelvis_tx','pelvis_tilt','pelvis_rotation']+block]
        if activation_only:names=['pelvis_tx','pelvis_tilt','pelvis_rotation']
        internal=[n for n in q if not n.startswith('pelvis_')]
        arms=stream.moment_arms(muscles=list(before['muscles']),coordinates=internal);write('native_moment_arms.json',arms)
        selected={f'gait2392_{m}_{side}' for side in ('r','l') for m in ('ercspn','intobl','extobl')};groups=[]
        for coordinate in [n for n in internal if n.startswith('hip_')]:
            for sign in (-1,1):
                candidates=[dict(muscle=m,moment_arm_m=row[coordinate],capacity_nm=abs(row[coordinate])*before['muscles'][m]['max_isometric_force_n']) for m,row in arms['moment_arms_m'].items() if row[coordinate]*sign>.01]
                candidates.sort(key=lambda value:value['capacity_nm'],reverse=True);chosen=candidates[:2];selected.update(value['muscle'] for value in chosen);groups.append(dict(coordinate=coordinate,sign=sign,selected=chosen))
        fixed_activations=dict(seed['activations']) if distal_only else {}
        if distal_only:selected=set(fixed_activations)
        if activation_only or retain_selection:selected=set(seed['activations'])
        tonic_names=sorted(selected);activation_names=[] if distal_only and not activation_only else tonic_names
        write('selection.json',dict(actuators=tonic_names,groups=groups,held_fixed=distal_only,retained_selection_source=seed.get('activation_selection_path') if distal_only else None,criteria='In distal-only mode retain seed tones exactly; otherwise Top2 native Fmax*abs(moment arm) per hip axis/sign at this actually initialized balanced pose, abs arm>0.01m; plus six lumbar.'))
        ranges={c.get('name'):list(map(float,c.findtext('range').split())) for c in ET.parse(output/'native/inputs/subject_walk_scaled.osim').findall('.//Coordinate')}
        for n in ('pelvis_tilt','pelvis_rotation','lumbar_extension','lumbar_bending','lumbar_rotation'):ranges[n]=[max(-.15,ranges[n][0]),min(.15,ranges[n][1])]
        bounds=np.array([ranges[n] for n in names]+[[.01,1.] for _ in activation_names]).T
        capacities={}
        for n in internal:
            all_pos=sum(max(0.,row[n])*before['muscles'][m]['max_isometric_force_n'] for m,row in arms['moment_arms_m'].items());all_neg=sum(max(0.,-row[n])*before['muscles'][m]['max_isometric_force_n'] for m,row in arms['moment_arms_m'].items())
            sel_pos=sum(max(0.,arms['moment_arms_m'][m][n])*before['muscles'][m]['max_isometric_force_n'] for m in selected);sel_neg=sum(max(0.,-arms['moment_arms_m'][m][n])*before['muscles'][m]['max_isometric_force_n'] for m in selected)
            capacities[n]=dict(all_positive_nm=all_pos,all_negative_nm=all_neg,selected_positive_nm=sel_pos,selected_negative_nm=sel_neg,normalization_nm=max(1.,all_pos,all_neg),uses_engineering_1nm_floor=max(all_pos,all_neg)<1.)
        write('torque_normalization.json',dict(reference_pose=q,coordinates=capacities,basis='Fixed reference normalization by max positive/negative sum Fmax*native moment arm over all98 modeled muscles, not claimed calibrated strength. Selected-tone directional capacities reported separately.1Nm explicit engineering floor where modeled torque capacity is smaller; dependent knee beta uses associated independent knee scale.'))
        manifest=json.loads((ROOT/support['manifest_path']).read_text())
        with np.load(ROOT/manifest['arrays_path']) as arrays:length=float(np.ptp(arrays['reference_points_source_m'][:,1]))
        weight=77.6122029*9.81;root_names=['pelvis_tx','pelvis_ty','pelvis_tz','pelvis_tilt','pelvis_list','pelvis_rotation'];root_variables=[names.index(n) for n in ('pelvis_tx','pelvis_tilt','pelvis_rotation')]
        def evaluate(x):
            nonlocal best
            if time.monotonic()-started>seconds:raise TimeoutError('Supported joint-force refinement wall cap')
            requested=dict(sorted({**q,**dict(zip(names,map(float,x[:len(names)])))}.items()));activations={**fixed_activations,**dict(zip(activation_names,map(float,x[len(names):])))};args=['evaluate_static_pose',str(len(requested))]
            for n,v in requested.items():args.extend([n,str(v)])
            args.append(str(len(activations)))
            for n,v in activations.items():args.extend([n,str(v)])
            native=stream._request(' '.join(args));mobility_names=native['mobility_coordinate_names'];mapping=dict(zip(mobility_names,range(len(mobility_names))))
            scales=[]
            for n,rotational in zip(mobility_names,native['mobility_rotational']):
                if n.startswith('pelvis_'):scales.append(weight*length if rotational else weight)
                else:scales.append(capacities[n.replace('_beta','')]['normalization_nm'])
            force=np.asarray(native['constrained_zero_acceleration_residual_mobility_force'])/scales
            force_names=mobility_names
            if distal_only:
                force_names=block;force=force[[mapping[n] for n in force_names]]
            if objective=='acceleration':
                force_names=mobility_names;force=np.asarray(native['udot'])/np.array([9.81/length if rotational else 9.81 for rotational in native['mobility_rotational']])
            tree=native['tree_zero_acceleration_residual_mobility_force'];all_root=np.array([tree[mapping[n]]/(weight if i<3 else weight*length) for i,n in enumerate(root_names)]);s=all_root[[0,3,5]]
            entry=dict(coordinates={n:native['coordinates'][n]['value'] for n in requested},requested_coordinates=requested,activations=activations,native=native,normalized_joint_force=force.tolist(),objective_mobility_names=force_names,objective=objective,all_root_residual=all_root.tolist(),support=s.tolist(),cost=float(force@force),feasible=bool(max(abs(all_root))<=1e-7))
            records.append(dict(evaluation=len(records)+1,cost=entry['cost'],maximum_root_residual=float(max(abs(all_root))),maximum_abs_udot=max(map(abs,native['udot'])),feasible=entry['feasible']));write('evaluations.json',records)
            if entry['feasible'] and (best is None or entry['cost']<best['cost']):best=entry;write('best_candidate.json',best)
            return entry
        def combined(entry):return np.r_[entry['normalized_joint_force'],entry['support']]
        x=np.array([q[n] for n in names]+[seed.get('activations',{}).get(n,before['muscles'][n]['activation']) for n in activation_names]);x=np.clip(x,bounds[0]+1e-10,bounds[1]-1e-10)
        first=evaluate(x);write('initial.json',first)
        if not first['feasible']:raise ValueError('Initial candidate lacks actual external wrench balance')
        iterations=[];coordinate_radius=.03;translation_radius=.001
        for iteration in range(20):
            base=evaluate(x);base_residual=combined(base);rows=len(base['normalized_joint_force']);jac=np.zeros((len(base_residual),len(x)))
            for i in range(len(x)):
                h=1e-5 if i<len(names) else 1e-4
                if x[i]+h>bounds[1,i]:h=-h
                derivative=None
                for multiplier in (1.,-1.,.1,-.1):
                    y=x.copy();y[i]+=h*multiplier
                    if not bounds[0,i]<=y[i]<=bounds[1,i]:continue
                    try:derivative=(combined(evaluate(y))-base_residual)/(h*multiplier);break
                    except ValueError as error:
                        if 'domain' not in str(error):raise
                        rejections.append(dict(kind='derivative',iteration=iteration,column=i,error=str(error)));write('rejected_domain_trials.json',rejections)
                if derivative is None:raise ValueError('No admissible finite-difference direction')
                jac[:,i]=derivative
            # Exact internal-force invariance, independently measured in the fixed-q audit.
            maximum_raw_activation_root_derivative=float(np.max(np.abs(jac[rows:,len(names):]))) if activation_names else 0.;jac[rows:,len(names):]=0.
            radius=np.array([translation_radius if n=='pelvis_tx' else coordinate_radius for n in names]+[.02]*len(activation_names));scale=radius/.03
            step_scaled,diagnostic=constrained_local_step(jac[:rows]*scale,np.asarray(base['normalized_joint_force']),x/scale,np.array([bounds[0]/scale,bounds[1]/scale]).T,jac[rows:]*scale,np.asarray(base['support']),radius=.03)
            direction=step_scaled*scale;root_block=jac[rows:,root_variables];accepted=False
            for fraction in (1.,.5,.25,.125,.0625,.03125):
                trial=x+fraction*direction
                try:
                    for correction in range(6):
                        value=evaluate(trial)
                        if value['feasible']:
                            if value['cost']<base['cost']:
                                x=trial;accepted=True
                            break
                        adjustment=np.linalg.solve(root_block,-np.asarray(value['support']));next_trial=trial.copy();next_trial[root_variables]+=adjustment
                        if np.any(next_trial<np.maximum(bounds[0],x-radius)-1e-12) or np.any(next_trial>np.minimum(bounds[1],x+radius)+1e-12):break
                        trial=next_trial
                except ValueError as error:
                    if 'domain' not in str(error):raise
                    rejections.append(dict(kind='trial',iteration=iteration,fraction=fraction,error=str(error)));write('rejected_domain_trials.json',rejections)
                if accepted:break
            diagnostic['predicted_normalized_force_norm']=diagnostic.pop('predicted_acceleration_norm')
            iterations.append(dict(iteration=iteration,accepted=accepted,fraction=fraction,cost=best['cost'],coordinate_radius=coordinate_radius,translation_radius=translation_radius,maximum_raw_activation_root_derivative=maximum_raw_activation_root_derivative,linear_step=diagnostic));write('iterations.json',iterations)
            if accepted and fraction==1:coordinate_radius=min(.03,coordinate_radius*1.5);translation_radius=min(.005,translation_radius*1.5)
            else:coordinate_radius*=.5;translation_radius*=.5
            if not accepted and coordinate_radius<1e-5:message='No feasible decreasing step after radius reductions';break
    except (TimeoutError,ValueError) as error:message=str(error)
    finally:
        if best is not None:
            observed=stream._request('observe');native=best['native'];artifact={**seed,'coordinates':dict(sorted(best['coordinates'].items())),'requested_coordinates':best['requested_coordinates'],'activations':best['activations'],'native_static_residual':native,'activation_selection_path':str((output/'selection.json').relative_to(ROOT)),'external_support_residual':best['all_root_residual'],'external_support_order':root_names,'protocol_sha256':hashlib.sha256(protocol).hexdigest(),'accepted_equilibrium':False,'objective':objective,'scope':'For acceleration objective: full udot normalized by gravity (translations) and gravity/body length (rotations), an engineering aggregate. For torque objective: strength-normalized joint force refinement with actual external root wrench enforced by linear equality and nonlinear rigid-root correction. No geometry/domain/guard change; tonic inverse statics, not learned brain control.'}
            write('initial_pose.json',artifact)
            report=dict(output=str(output),message=message,wall_s=time.monotonic()-started,evaluations=len(records),domain_rejections=len(rejections),selected_actuators=tonic_names,distal_only=distal_only,objective_mobility_names=best['objective_mobility_names'],objective=objective,activation_only=activation_only,initial_cost=first['cost'],final_cost=best['cost'],external_wrench_balanced=best['feasible'],root_residual_normalized=best['all_root_residual'],root_order=root_names,accepted_equilibrium=False,max_abs_udot=max(map(abs,native['udot'])),udot_norm=float(np.linalg.norm(native['udot'])),contact_force_n=native['contact_force_n'],continuing_state_unchanged=all(before[k]==observed[k] for k in ('time_s','coordinates','muscles')),largest_accelerations=sorted(zip(native['mobility_coordinate_names'],native['udot']),key=lambda v:abs(v[1]),reverse=True)[:10])
            write('report.json',report);print(json.dumps(report,indent=2))
        stream.close()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--seed',type=Path,default=ROOT/'data/derived/supine-external-support-aiui5g4g/initial_pose.json');p.add_argument('--seconds',type=float,default=90);p.add_argument('--distal-only',action='store_true');p.add_argument('--joint-block');p.add_argument('--objective',choices=['torque','acceleration'],default='torque');p.add_argument('--activation-only',action='store_true');p.add_argument('--retain-selection',action='store_true');a=p.parse_args();run(a.seed,a.seconds,a.distal_only,a.joint_block,a.objective,a.activation_only,a.retain_selection)
