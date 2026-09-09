#!/usr/bin/env python3
"""Export exact native equilibrium state in an existing LQR chart, read-only."""
from pathlib import Path
import argparse,hashlib,json,sys,tempfile
import xml.etree.ElementTree as ET
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.stance_lqr import NativeStanceLQR
ROOT=Path(__file__).resolve().parents[1]
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()

def export(registration,policy_path,output,snapshot_path):
    path=ROOT/registration;reg=json.loads(path.read_text());out=Path(output);out.mkdir(parents=True,exist_ok=False)
    with np.load(policy_path,allow_pickle=False) as policy:
        state_names=policy['state_names'].tolist();muscle_names=policy['muscle_names'].tolist();source_model_sha256=str(policy['model_sha256'].item());mass=float(policy['target_mass_kg'].item())
    if abs(mass-reg['target_mass_kg'])>1e-9:raise ValueError('Target/policy mass mismatch')
    source_reg=json.loads((ROOT/'data/derived/mechanics/patient_stance98/registration.json').read_text())
    if source_reg['model_sha256']!=source_model_sha256:raise ValueError('Policy source model must match patient source')
    def canonical(node):return(node.tag,tuple(sorted(node.attrib.items())),(node.text or '').strip(),tuple(canonical(c) for c in node if c.tag!='default_activation'))
    if canonical(ET.parse(ROOT/source_reg['model_path']).getroot())!=canonical(ET.parse(ROOT/reg['model_path']).getroot()):raise ValueError('Target model changes physical parameters beyond initial activation defaults')
    state=json.loads(Path(snapshot_path).read_text())
    policy=NativeStanceLQR(policy_path)
    if set(muscle_names)!=set(state['muscles']):raise ValueError('Exact muscle catalog required')
    x=policy.state_vector(state);u=np.array([state['muscles'][n]['activation'] for n in muscle_names])
    if not np.isfinite(x).all() or not np.isfinite(u).all():raise ValueError('Nonfinite native target')
    np.savez(out/'target.npz',x_target=x,u_target=u,x0=x,u0=u,state_names=np.array(state_names),muscle_names=np.array(muscle_names),target_mass_kg=mass,source_model_sha256=source_model_sha256,target_model_sha256=reg['model_sha256'])
    q=reg['initial_pose']
    metadata={'schema':'ihm.native-equilibrium-target.v1','target_registration':registration,'target_registration_sha256':sha(path),'policy_source_path':str(Path(policy_path).resolve().relative_to(ROOT)),'policy_sha256':sha(Path(policy_path)),'source_model_sha256':source_model_sha256,'target_model_sha256':reg['model_sha256'],'same_physical_model_except_initial_activation_defaults':True,'target_mass_kg':mass,'state_count':len(state_names),'muscle_count':len(muscle_names),'state_path_order_matches_policy':True,'state_values_from':'Retained actual native initialization snapshot after muscle equilibration and final initial pose; charted with NativeStanceLQR.state_vector, includes actual joint values/speeds and all muscle activations/fiber lengths','snapshot_path':str(Path(snapshot_path).resolve().relative_to(ROOT)),'snapshot_sha256':sha(Path(snapshot_path)),'initial_pose':q,'target_sha256':sha(out/'target.npz'),'physical_time_advanced_s':0,'weight_transfer_demonstrated':False,'prescribed_motion':False,'dynamic_transition_validation_required':True,'new_native_process_started':False}
    (out/'manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--registration',required=True);p.add_argument('--policy',required=True);p.add_argument('--output',required=True);p.add_argument('--snapshot',required=True);args=p.parse_args();export(args.registration,Path(args.policy),Path(args.output),Path(args.snapshot))
