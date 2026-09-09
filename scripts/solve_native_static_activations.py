#!/usr/bin/env python3
"""Finite-difference native static activation identification at a fixed pose.

Optimizes leg/lumbar force residuals only. Root and unactuated arm residuals are
reported separately and cannot be silently claimed as a full static equilibrium.
"""
from pathlib import Path
import json,sys,tempfile,time
import numpy as np
from scipy.optimize import lsq_linear
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(tempfile.mkdtemp(prefix='native_static_activation_solve_',dir=ROOT/'data/research/locomotion_control'))
native=NativeMechanicalStream(ROOT,OUT/'native',environment='upright',target_mass_kg=70,augmented_registration='data/derived/mechanics/initial_stance86/registration.json')
try:
    initial=native.snapshot();names=list(initial['muscles']);q={'pelvis_ty':initial['coordinates']['pelvis_ty']['value']}
    a=np.array([initial['muscles'][n]['activation'] for n in names]);history=[];evaluations=[]
    def evaluate(values,label):
        start=time.monotonic();result=native.evaluate_static_pose(q,activations=dict(zip(names,map(float,values))))
        all_residual=dict(zip(result['mobility_coordinate_names'],result['constrained_zero_acceleration_residual_mobility_force']))
        selected=[n for n in result['mobility_coordinate_names'] if n.startswith(('hip_','knee_angle','ankle_','lumbar_')) and result['coordinates'][n]['independent']]
        residual=np.array([all_residual[n] for n in selected])
        path=OUT/f'evaluation_{len(evaluations):03d}.json';path.write_text(json.dumps(result,allow_nan=False))
        evaluations.append({'label':label,'selected_residual_norm_nm':float(np.linalg.norm(residual)),'path':str(path.relative_to(ROOT)),'wall_seconds':time.monotonic()-start})
        return residual,result,selected,all_residual
    residual,result,selected,all_residual=evaluate(a,'initial')
    for iteration in range(3):
        jacobian=np.zeros((len(residual),len(names)))
        for j,name in enumerate(names):
            probe=a.copy();delta=.02 if a[j]<.97 else -.02;probe[j]+=delta
            jacobian[:,j]=(evaluate(probe,f'iteration{iteration}_jacobian_{name}')[0]-residual)/delta
        penalty=.2
        matrix=np.vstack((jacobian,penalty*np.eye(len(names))))
        rhs=np.r_[-residual,np.zeros(len(names))]
        step=lsq_linear(matrix,rhs,bounds=(.01-a,1.-a),tol=1e-8,max_iter=300).x
        proposal=np.clip(a+step,.01,1)
        candidate,candidate_result,_,candidate_all=evaluate(proposal,f'iteration{iteration}_candidate')
        before=float(np.linalg.norm(residual));after=float(np.linalg.norm(candidate));accepted=after<before
        history.append({'iteration':iteration,'before_nm':before,'after_nm':after,'accepted':bool(accepted)})
        if accepted:a,residual,result,all_residual=proposal,candidate,candidate_result,candidate_all
        print(json.dumps(history[-1]),flush=True)
        if np.linalg.norm(residual)<.01:break
    report={'schema':'ihm.native-static-muscle-activation-fit.v1','registration':'data/derived/mechanics/initial_stance86/registration.json','activations':dict(zip(names,map(float,a))),'selected_coordinates':selected,'selected_residual_norm_nm':float(np.linalg.norm(residual)),'selected_residual_nm':{n:all_residual[n] for n in selected},'excluded_root_and_arm_residual':{n:v for n,v in all_residual.items() if n not in selected},'all_coordinate_acceleration':{n:r['acceleration'] for n,r in result['coordinates'].items()},'history':history,'evaluations':evaluations,'physical_time_advanced_s':0,'full_static_equilibrium_claimed':False,'sustained_stance_demonstrated':False,'initial_source_activations':{n:initial['muscles'][n]['activation'] for n in names},'candidate_state_variables':result['state_variables']}
    (OUT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'report':str(OUT/'report.json'),'residual_nm':report['selected_residual_norm_nm'],'excluded':report['excluded_root_and_arm_residual']},indent=2),flush=True)
finally:native.close()
