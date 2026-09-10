#!/usr/bin/env python3
"""Solve initial passive arm/toe pose and muscle activation force balance.

All queries operate on copied native States. The output changes only INITIAL
coordinate/activation defaults; no support actuator or motion constraint exists.
"""
from pathlib import Path
import argparse,hashlib,json,re,sys,tempfile
import numpy as np
from scipy.optimize import least_squares,lsq_linear
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.assembly.embodied import _prepare_mechanical_registration
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG
ROOT=Path(__file__).resolve().parents[1]
BASE='data/derived/mechanics/patient_rightlift_stance98/registration.json'
SEED='data/derived/mechanics/patient_rightlift_stance98/registration.json'
OUT=ROOT/'data/derived/mechanics/patient_rightlift5mm_stance98'
TARGET_MASS_KG=MECHANICAL_TARGET_MASS_KG
sha=lambda data:hashlib.sha256(data).hexdigest()

def main():
    if OUT.exists():raise ValueError('Refusing to overwrite candidate')
    work=Path(tempfile.mkdtemp(prefix='resting_stance98_',dir=ROOT/'data/research/locomotion_control'))
    _,base,_=_prepare_mechanical_registration(ROOT,BASE)
    _,seed,_=_prepare_mechanical_registration(ROOT,SEED)
    FIT=seed['static_resting_solve']
    native=NativeMechanicalStream(ROOT,work/'native',environment='upright',target_mass_kg=TARGET_MASS_KG,augmented_registration=BASE,initial_pose=base.get('initial_pose'))
    initial_native=native.snapshot()
    fit=json.loads((ROOT/FIT).read_text());names=list(native.state['muscles']);activation=np.array([fit['activations'].get(n,.01) for n in names]);evaluations=[]
    def coordinates(x):
        a,h,listing,hip_common=x[:4]
        q={'pelvis_tilt':-float(a),'pelvis_ty':float(h),'pelvis_list':float(listing),'pelvis_tz':0.}
        for i,(side,sign) in enumerate([('r',1),('l',-1)]):
            af,aa,ar,el,ps,mtp=x[4+i*6:10+i*6]
            q.update({f'ankle_angle_{side}':float(a),f'hip_adduction_{side}':sign*(float(hip_common)-float(listing)),f'arm_flex_{side}':float(af),f'arm_add_{side}':float(aa),f'arm_rot_{side}':float(ar),f'elbow_flex_{side}':float(el),f'pro_sup_{side}':float(ps),f'mtp_angle_{side}':float(mtp)})
        q.update(hip_flexion_r=.1115, knee_angle_r=.223, ankle_angle_r=float(a)+.1115)
        return q
    def evaluate(q,a,label):
        result=native.evaluate_static_pose(q,activations=dict(zip(names,map(float,a))))
        residual=dict(zip(result['mobility_coordinate_names'],result['constrained_zero_acceleration_residual_mobility_force']))
        path=work/f'query_{len(evaluations):04d}.json';path.write_text(json.dumps(result,allow_nan=False));evaluations.append({'label':label,'path':str(path.relative_to(ROOT))})
        return result,residual
    def pose_residual(x):
        result,r=evaluate(coordinates(x),activation,'pose')
        # Bilateral passive force equations plus root support force/moments. Their
        # stated scaling balances force(N) and torque(Nm) numerical objectives.
        values=[r[n+'_'+side] for side in ('r','l') for n in ('arm_flex','arm_add','arm_rot','elbow_flex','pro_sup','mtp_angle')]
        load={side:sum(c['force_n'][1] for c in result['contacts'] if c['name'].startswith('contact') and c['name'].endswith('_'+side)) for side in ('r','l')}
        return np.r_[np.array(values)/10,r['pelvis_tilt']/10,r['pelvis_ty']/100,r['pelvis_list']/10,(x[3]-(reference['hip_adduction_r']+reference['pelvis_list']))*10]
    reference=fit['coordinates']
    x=np.array([-reference['pelvis_tilt'],reference['pelvis_ty'],reference['pelvis_list'],reference['hip_adduction_r']+reference['pelvis_list']]+[reference[n+'_'+side] for side in ('r','l') for n in ('arm_flex','arm_add','arm_rot','elbow_flex','pro_sup','mtp_angle')])
    bounds=([-.3,.95,-.2,-.25]+[-1,-1,-.5,.01,0,-.5]*2,[.5,1.3,.1,.25]+[1,.1,.5,2.5,2.08,.5]*2)
    history=[]
    try:
        for outer in range(3):
            solved=least_squares(pose_residual,x,bounds=bounds,x_scale='jac',diff_step=1e-4,max_nfev=80,xtol=1e-10,ftol=1e-10,gtol=1e-9)
            x=solved.x;q=coordinates(x)
            result,residual=evaluate(q,activation,'post_pose')
            print(json.dumps({'stage':'pose','outer':outer,'x':x.tolist(),'scaled_norm':float(np.linalg.norm(solved.fun)),'root_ty':residual['pelvis_ty'],'elbow':residual['elbow_flex_r'],'mtp':residual['mtp_angle_r']}),flush=True)
            selected=[n for n,c in result['coordinates'].items() if c['independent'] and n.startswith(('hip_','knee_','ankle_','lumbar_'))];r=np.array([residual[n] for n in selected])
            for iteration in range(3):
                jacobian=np.zeros((len(r),len(names)))
                for j,name in enumerate(names):
                    probe=activation.copy();delta=.02 if activation[j]<.97 else -.02;probe[j]+=delta
                    _,response=evaluate(q,probe,'activation_jacobian_'+name);jacobian[:,j]=(np.array([response[n] for n in selected])-r)/delta
                step=lsq_linear(np.vstack((jacobian,.1*np.eye(len(names)))),np.r_[-r,np.zeros(len(names))],bounds=(.01-activation,1.-activation),tol=1e-8,max_iter=300).x
                proposal=np.clip(activation+step,.01,1);candidate,newresidual=evaluate(q,proposal,'activation_candidate');newr=np.array([newresidual[n] for n in selected])
                before=float(np.linalg.norm(r));after=float(np.linalg.norm(newr));history.append({'outer':outer,'iteration':iteration,'before_nm':before,'after_nm':after})
                if after<before:activation,r,result,residual=proposal,newr,candidate,newresidual
                print(json.dumps(history[-1]),flush=True)
                if np.linalg.norm(r)<.001:break
        result,residual=evaluate(q,activation,'final')
        independent={n:v for n,v in residual.items() if result['coordinates'][n]['independent']}
        report={'schema':'ihm.native-patient-rightlift-static.v1','target_left_load_fraction':1.0,'weight_transfer_demonstrated':False,'target_mass_kg':TARGET_MASS_KG,'source_registration':BASE,'seed_registration':SEED,'coordinates':q,'activations':dict(zip(names,map(float,activation))),'all_independent_residual':independent,'all_independent_residual_norm_mixed_units':float(np.linalg.norm(list(independent.values()))),'all_coordinate_acceleration':{n:c['acceleration'] for n,c in result['coordinates'].items()},'selected_residual_norm_nm':float(np.linalg.norm([residual[n] for n in selected])),'history':history,'evaluations':evaluations,'candidate_state_variables':result['state_variables'],'physical_time_advanced_s':0,'external_support':False,'prescribed_motion':False,'sustained_stance_demonstrated':False,'full_static_equilibrium_claimed':False}
        reportpath=work/'report.json';reportpath.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    finally:native.close()
    model=(ROOT/base['model_path']).read_bytes()
    # Canonical source coordinate defaults and material registration remain exact.
    # The separately declared initial_pose is applied by native startup only.
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
    registration=dict(base);registration['sources']=dict(base['sources']);registration['sources'].update({BASE:sha((ROOT/BASE).read_bytes()),base['model_path']:base['model_sha256'],FIT:sha((ROOT/FIT).read_bytes()),SEED:sha((ROOT/SEED).read_bytes()),str(reportpath.relative_to(ROOT)):sha(reportpath.read_bytes()),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__).read_bytes())})
    registration.update(model_path=str(modelpath.relative_to(ROOT)),model_sha256=sha(model),initial_pose=q,initial_pose_only_coordinates=q,target_mass_kg=TARGET_MASS_KG,equilibrium_excitations=report['activations'],static_resting_solve=str(reportpath.relative_to(ROOT)),native_acceptance_complete=False)
    regpath=OUT/'registration.json';regpath.write_text(json.dumps(registration,indent=2)+'\n');_prepare_mechanical_registration(ROOT,str(regpath.relative_to(ROOT)))
    provisional=NativeMechanicalStream(ROOT,OUT/'native_before_foot_recentering',environment='upright',target_mass_kg=TARGET_MASS_KG,augmented_registration=str(regpath.relative_to(ROOT)),initial_pose=q)
    try:
        current=provisional.snapshot()
        def midpoint_z(state):return float(np.mean([c['center_m'][2] for c in state['contacts'] if c['name'].startswith('contact') and c['name'].endswith(('_r','_l'))]))
        displacement=midpoint_z(initial_native)-midpoint_z(current)
        q['pelvis_tz']+=displacement
        report['initial_foot_midpoint_recenter_z_m']=displacement
        report['coordinates']=q
    finally:provisional.close()
    reportpath.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    registration['sources'][str(reportpath.relative_to(ROOT))]=sha(reportpath.read_bytes())
    registration['initial_pose']=q;registration['initial_pose_only_coordinates']=q
    regpath.write_text(json.dumps(registration,indent=2)+'\n')
    verify=NativeMechanicalStream(ROOT,OUT/'native_initialization',environment='upright',target_mass_kg=TARGET_MASS_KG,augmented_registration=str(regpath.relative_to(ROOT)),initial_pose=q)
    try:
        state=verify.snapshot();actual=verify.evaluate_static_pose(q)
        actual_residual=dict(zip(actual['mobility_coordinate_names'],actual['constrained_zero_acceleration_residual_mobility_force']))
        error=max(abs(state['muscles'][n]['activation']-v) for n,v in report['activations'].items())
        if error>1e-12:raise ValueError('Actual native initial activations differ')
        foot_force={side:sum(c['force_n'][1] for c in state['contacts'] if c['name'].startswith('contact') and c['name'].endswith('_'+side)) for side in ('r','l')}
        acceptance={'right_minimum_contact_clearance_m':min(c['center_m'][1]-c['radius_m'] for c in state['contacts'] if c['name'].startswith('contact') and c['name'].endswith('_r')),'foot_vertical_force_n':foot_force,'actual_left_load_fraction':foot_force['l']/sum(foot_force.values()),'schema':'ihm.patient-rightlift-initial-acceptance.v1','target_left_load_fraction':1.0,'weight_transfer_demonstrated':False,'target_mass_kg':TARGET_MASS_KG,'actual_native_mass_kg':state['mass_kg'],'activation_default_max_error':error,'all_independent_residual':{n:v for n,v in actual_residual.items() if actual['coordinates'][n]['independent']},'maximum_absolute_coordinate_acceleration':max(abs(c['acceleration']) for c in actual['coordinates'].values()),'physical_time_advanced_s':0,'canonical_source_coordinate_defaults_unchanged':True,'initial_pose':q,'sustained_stance_demonstrated':False}
        (OUT/'initial_snapshot.json').write_text(json.dumps(state,allow_nan=False));acceptancepath=OUT/'native_initial_acceptance.json';acceptancepath.write_text(json.dumps(acceptance,indent=2)+'\n')
        (OUT/'initial_pose.json').write_text(json.dumps(q,indent=2)+'\n');(OUT/'equilibrium_excitations.json').write_text(json.dumps(report['activations'],indent=2)+'\n')
    finally:verify.close()
    registration['sources'][str(acceptancepath.relative_to(ROOT))]=sha(acceptancepath.read_bytes());registration['native_initial_acceptance']=str(acceptancepath.relative_to(ROOT));regpath.write_text(json.dumps(registration,indent=2)+'\n')
    _prepare_mechanical_registration(ROOT,str(regpath.relative_to(ROOT)))
    print(json.dumps({'registration':str(regpath.relative_to(ROOT)),'report':str(reportpath),'independent_norm':report['all_independent_residual_norm_mixed_units'],**acceptance},indent=2),flush=True)
if __name__=='__main__':main()
