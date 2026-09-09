#!/usr/bin/env python3
"""Opt-in native resting86 activation refit with a0 >= .03 reserve.

The pose is fixed to the retained refined resting stance. Increased baseline
coactivation is an engineering controller prior, not learned recruitment.
"""
from pathlib import Path
import hashlib,json,re,sys,tempfile
import numpy as np
from scipy.optimize import lsq_linear
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.assembly.embodied import _prepare_mechanical_registration
ROOT=Path(__file__).resolve().parents[1]
BASE='data/derived/mechanics/resting_stance86_refined/registration.json'
OUT=ROOT/'data/derived/mechanics/coactive_stance86'
sha=lambda data:hashlib.sha256(data).hexdigest()

def main():
    if OUT.exists():raise ValueError('Refusing to overwrite coactivation candidate')
    _,base,_=_prepare_mechanical_registration(ROOT,BASE)
    q=dict(base['initial_pose_only_coordinates']);work=Path(tempfile.mkdtemp(prefix='coactive_stance86_',dir=ROOT/'data/research/locomotion_control'))
    native=NativeMechanicalStream(ROOT,work/'native',environment='upright',target_mass_kg=70,augmented_registration=BASE,initial_pose=q)
    initial=native.snapshot();names=list(initial['muscles']);a=np.array([max(.03,initial['muscles'][n]['activation']) for n in names]);history=[];evaluations=[]
    def evaluate(values,label):
        result=native.evaluate_static_pose(q,activations=dict(zip(names,map(float,values))))
        all_residual=dict(zip(result['mobility_coordinate_names'],result['constrained_zero_acceleration_residual_mobility_force']))
        selected=[n for n in result['mobility_coordinate_names'] if n.startswith(('hip_','knee_angle','ankle_','lumbar_')) and result['coordinates'][n]['independent']]
        residual=np.array([all_residual[n] for n in selected]);path=work/f'evaluation_{len(evaluations):03d}.json';path.write_text(json.dumps(result,allow_nan=False));evaluations.append({'label':label,'selected_residual_norm_nm':float(np.linalg.norm(residual)),'path':str(path.relative_to(ROOT))})
        return residual,result,selected,all_residual
    try:
        residual,result,selected,all_residual=evaluate(a,'coactive_floor_initial')
        for iteration in range(4):
            jacobian=np.zeros((len(residual),len(names)))
            for j,name in enumerate(names):
                probe=a.copy();delta=.02 if a[j]<.97 else -.02;probe[j]+=delta
                jacobian[:,j]=(evaluate(probe,f'iteration{iteration}_jacobian_{name}')[0]-residual)/delta
            matrix=np.vstack((jacobian,.1*np.eye(len(names))));rhs=np.r_[-residual,np.zeros(len(names))]
            step=lsq_linear(matrix,rhs,bounds=(.03-a,1.-a),tol=1e-9,max_iter=300).x
            proposal=np.clip(a+step,.03,1);candidate,candidate_result,_,candidate_all=evaluate(proposal,f'iteration{iteration}_candidate')
            before=float(np.linalg.norm(residual));after=float(np.linalg.norm(candidate));accepted=after<before
            history.append({'iteration':iteration,'before_nm':before,'after_nm':after,'accepted':bool(accepted)})
            if accepted:a,residual,result,all_residual=proposal,candidate,candidate_result,candidate_all
            print(json.dumps(history[-1]),flush=True)
            if np.linalg.norm(residual)<.0001:break
        if np.linalg.norm(residual)>.01:raise RuntimeError('Coactive selected-force solve did not reach.01Nm')
        report={'schema':'ihm.native-coactive-stance86-static.v1','base_registration':BASE,'coordinates':q,'activations':dict(zip(names,map(float,a))),'activation_lower_bound':.03,'native_activation_floor':.01,'minimum_downward_activation_reserve':float(a.min()-.01),'selected_coordinates':selected,'selected_residual_norm_nm':float(np.linalg.norm(residual)),'all_independent_residual':{n:v for n,v in all_residual.items() if result['coordinates'][n]['independent']},'excluded_root_arm_and_toe_residual':{n:v for n,v in all_residual.items() if n not in selected},'all_coordinate_acceleration':{n:c['acceleration'] for n,c in result['coordinates'].items()},'history':history,'evaluations':evaluations,'initial_pose_unchanged':True,'physical_time_advanced_s':0,'external_support':False,'prescribed_motion':False,'sustained_stance_demonstrated':False,'full_static_equilibrium_claimed':False,'coactivation_basis':'Engineering bidirectional activation reserve for feedback, not measured recruitment or brain training.'}
        reportpath=work/'report.json';reportpath.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    finally:native.close()
    model=(ROOT/base['model_path']).read_bytes()
    for name,value in report['activations'].items():
        pattern=rb'(<(?P<kind>Millard2012EquilibriumMuscle|Thelen2003Muscle) name="'+re.escape(name.encode())+rb'">)(?P<body>.*?)(</(?P=kind)>)'
        def replace(m):
            body=m['body'];v=format(value,'.17g').encode()
            if b'<default_activation>' in body:body=re.sub(rb'(<default_activation>).*?(</default_activation>)',lambda n:n[1]+v+n[2],body,count=1,flags=re.S)
            else:body=b'\n<default_activation>'+v+b'</default_activation>\n'+body
            return m[1]+body+m[4]
        model,count=re.subn(pattern,replace,model,count=1,flags=re.S)
        if count!=1:raise ValueError('Missing muscle')
    OUT.mkdir(parents=True);modelpath=OUT/'model.osim';modelpath.write_bytes(model)
    registration=dict(base);registration['sources']=dict(base['sources']);registration['sources'].update({BASE:sha((ROOT/BASE).read_bytes()),base['model_path']:base['model_sha256'],str(reportpath.relative_to(ROOT)):sha(reportpath.read_bytes()),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__).read_bytes())})
    registration.update(model_path=str(modelpath.relative_to(ROOT)),model_sha256=sha(model),equilibrium_excitations=report['activations'],static_coactivation_solve=str(reportpath.relative_to(ROOT)),native_acceptance_complete=False)
    regpath=OUT/'registration.json';regpath.write_text(json.dumps(registration,indent=2)+'\n');_prepare_mechanical_registration(ROOT,str(regpath.relative_to(ROOT)))
    verify=NativeMechanicalStream(ROOT,OUT/'native_initialization',environment='upright',target_mass_kg=70,augmented_registration=str(regpath.relative_to(ROOT)),initial_pose=q)
    try:
        state=verify.snapshot();actual=verify.evaluate_static_pose(q);r=dict(zip(actual['mobility_coordinate_names'],actual['constrained_zero_acceleration_residual_mobility_force']));error=max(abs(state['muscles'][n]['activation']-v) for n,v in report['activations'].items())
        if error>1e-12 or min(m['activation'] for m in state['muscles'].values())<.03-1e-12:raise ValueError('Native activation initialization mismatch')
        acceptance={'activation_default_max_error':error,'minimum_actual_activation':min(m['activation'] for m in state['muscles'].values()),'selected_residual_norm_nm':float(np.linalg.norm([r[n] for n in selected])),'all_independent_residual':{n:v for n,v in r.items() if actual['coordinates'][n]['independent']},'maximum_abs_coordinate_acceleration':max(abs(c['acceleration']) for c in actual['coordinates'].values()),'physical_time_advanced_s':0,'sustained_stance_demonstrated':False}
        (OUT/'initial_snapshot.json').write_text(json.dumps(state,allow_nan=False));acceptancepath=OUT/'native_initial_acceptance.json';acceptancepath.write_text(json.dumps(acceptance,indent=2)+'\n')
    finally:verify.close()
    registration['sources'][str(acceptancepath.relative_to(ROOT))]=sha(acceptancepath.read_bytes());registration['coactive_initial_acceptance']=str(acceptancepath.relative_to(ROOT));regpath.write_text(json.dumps(registration,indent=2)+'\n')
    print(json.dumps({'registration':str(regpath.relative_to(ROOT)),'report':str(reportpath),**acceptance},indent=2),flush=True)
if __name__=='__main__':main()
