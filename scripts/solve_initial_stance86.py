#!/usr/bin/env python3
"""Solve initial pose geometry/contact using fresh native equilibrium snapshots.

Only default coordinates change; no prescribed trajectory, support actuator,
continuing-state teleportation, or dynamic-stability claim is introduced.
"""
from pathlib import Path
import hashlib,json,re,sys,tempfile
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.moment_arm_control import native_center_of_mass,foot_support_center
from ihm.assembly.embodied import _prepare_mechanical_registration
ROOT=Path(__file__).resolve().parents[1]
BASE='data/derived/mechanics/leg_lumbar86/registration.json'
OUT=ROOT/'data/derived/mechanics/initial_stance86'
sha=lambda raw:hashlib.sha256(raw).hexdigest()

def main():
    if OUT.exists():raise ValueError('Refusing to replace retained initial pose')
    _,base,catalog=_prepare_mechanical_registration(ROOT,BASE)
    raw=(ROOT/base['model_path']).read_bytes()
    scratch=Path(tempfile.mkdtemp(prefix='initial_stance86_solve_',dir=ROOT/'data/research/locomotion_control'))
    evaluations=[]
    def materialize(directory,angles):
        a,height=map(float,angles)
        coordinates={'pelvis_tilt':-a,'ankle_angle_r':a,'ankle_angle_l':a,'pelvis_ty':height}
        model=raw
        for name,value in coordinates.items():
            pattern=rb'(<Coordinate name="'+name.encode()+rb'">.*?<default_value>)(.*?)(</default_value>)'
            model,count=re.subn(pattern,lambda m:m[1]+format(value,'.17g').encode()+m[3],model,count=1,flags=re.S)
            if count!=1:raise ValueError('Coordinate default not uniquely found')
        directory.mkdir(parents=True)
        modelpath=directory/'model.osim';modelpath.write_bytes(model)
        reg=dict(base);reg['sources']=dict(base['sources']);reg['sources'][BASE]=sha((ROOT/BASE).read_bytes());reg['sources'][base['model_path']]=base['model_sha256']
        reg['model_path']=str(modelpath.relative_to(ROOT));reg['model_sha256']=sha(model);reg['initial_pose_only_coordinates']=coordinates
        reg['native_acceptance_complete']=False
        path=directory/'registration.json';path.write_text(json.dumps(reg,indent=2)+'\n')
        return path,reg
    def evaluate(x):
        directory=scratch/f'eval_{len(evaluations):03d}'
        path,reg=materialize(directory,x)
        native=NativeMechanicalStream(ROOT,directory/'native',environment='upright',target_mass_kg=70,augmented_registration=str(path.relative_to(ROOT)))
        try:state=native.snapshot()
        finally:native.close()
        com,velocity=native_center_of_mass(state);support=foot_support_center(state)
        feet=[c for c in state['contacts'] if c.get('body_frame') in ('calcn_r','calcn_l','toes_r','toes_l')]
        footforce=sum(c['force_n'][1] for c in feet)
        if support is None:support=np.mean([c['center_m'] for c in feet],axis=0)
        residual=[(com[0]-support[0])/.1,(footforce-70*9.81)/(70*9.81)]
        record={'parameters':list(map(float,x)),'coordinates':reg['initial_pose_only_coordinates'],'com_ground_m':com.tolist(),'support_center_ground_m':support.tolist(),'foot_vertical_force_n':footforce,'residual':residual,'snapshot_path':str((directory/'initial_snapshot.json').relative_to(ROOT))}
        (directory/'initial_snapshot.json').write_text(json.dumps(state,allow_nan=False))
        evaluations.append(record)
        print(json.dumps({'evaluation':len(evaluations),'x':record['parameters'],'residual':residual,'force':footforce}),flush=True)
        return residual
    result=least_squares(evaluate,[.12,1.035],bounds=([-.1,.95],[.4,1.3]),diff_step=1e-4,xtol=1e-8,ftol=1e-8,gtol=1e-7,max_nfev=35)
    final_residual=evaluate(result.x)
    if max(abs(v) for v in final_residual)>.005:raise RuntimeError('Static initial support solve did not converge')
    path,registration=materialize(OUT,result.x)
    receipt={'schema':'ihm.initial-stance86-static-solve.v1','solver':'scipy least_squares on isolated fresh native initialized states','parameters':result.x.tolist(),'residual':final_residual,'target_mass_kg':70,'target_foot_vertical_force_n':70*9.81,'termination':str(result.message),'evaluations':evaluations,'physical_time_advanced_s':0,'prescribed_coordinates':False,'external_support':False,'initial_muscle_state':'Native default activations and native equilibrateMuscles at changed pose; no dynamic equilibrium assertion','dynamic_stability_demonstrated':False,'limits':['COM and vertical load balance do not ensure joint torque equilibrium.','86-muscle search plant retains inactive articulated arms;98-muscle transfer requires separate validation.']}
    receipt_path=OUT/'static_solve.json';receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
    registration['sources'][str(receipt_path.relative_to(ROOT))]=sha(receipt_path.read_bytes())
    registration['sources'][str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__).read_bytes())
    registration['static_solve_receipt']=str(receipt_path.relative_to(ROOT))
    path.write_text(json.dumps(registration,indent=2)+'\n')
    _prepare_mechanical_registration(ROOT,str(path.relative_to(ROOT)))
    print(json.dumps({'registration':str(path.relative_to(ROOT)),'final':evaluations[-1],'preflight_passed':True}),flush=True)
if __name__=='__main__':main()
