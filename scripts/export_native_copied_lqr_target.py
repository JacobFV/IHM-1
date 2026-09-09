#!/usr/bin/env python3
"""Export a copied-equilibrium policy reference with a verified baseline chart."""
from pathlib import Path
import argparse,hashlib,json
import xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--registration',required=True);p.add_argument('--policy',required=True);p.add_argument('--baseline-policy',default='data/models/engineering_stance_v1/linearization.npz');p.add_argument('--copied-state',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    rp=ROOT/a.registration;pp=ROOT/a.policy;bp=ROOT/a.baseline_policy;cp=ROOT/a.copied_state;out=ROOT/a.output
    reg=json.loads(rp.read_text());copied=json.loads(cp.read_text());source_reg=json.loads((ROOT/'data/derived/mechanics/patient_stance98/registration.json').read_text())
    if out.exists():raise ValueError('Refusing to overwrite target')
    with np.load(pp,allow_pickle=False) as policy,np.load(bp,allow_pickle=False) as baseline:
        names=policy['state_names'].copy();muscles=policy['muscle_names'].copy();x=policy['x0'].copy();u=policy['u0'].copy();mass=float(policy['target_mass_kg'].item());model=str(policy['model_sha256'].item());source_model=str(baseline['model_sha256'].item())
        if not np.array_equal(names,baseline['state_names']) or not np.array_equal(muscles,baseline['muscle_names']):raise ValueError('Baseline and target charts differ')
        if mass!=float(baseline['target_mass_kg'].item()) or mass!=reg['target_mass_kg'] or mass!=copied['target_mass_kg']:raise ValueError('Mass mismatch')
    if model!=reg['model_sha256'] or model!=copied['source_model_sha256']:raise ValueError('Target model identity mismatch')
    if source_model!=source_reg['model_sha256']:raise ValueError('Baseline model identity mismatch')
    def canonical(n):return(n.tag,tuple(sorted(n.attrib.items())),(n.text or '').strip(),tuple(canonical(c) for c in n if c.tag!='default_activation'))
    if canonical(ET.parse(ROOT/source_reg['model_path']).getroot())!=canonical(ET.parse(ROOT/reg['model_path']).getroot()):raise ValueError('Physical model changed')
    expected_x=np.array([copied['state_variables'][n] for n in names]);expected_u=np.array([copied['equilibrium_excitations'][n] for n in muscles])
    if not np.array_equal(x,expected_x) or not np.array_equal(u,expected_u):raise ValueError('Policy reference does not exactly match copied equilibrium')
    if not np.isfinite(x).all() or not np.isfinite(u).all() or np.any(u<0) or np.any(u>1):raise ValueError('Invalid target')
    out.mkdir(parents=True)
    np.savez(out/'target.npz',x0=x,u0=u,x_target=x,u_target=u,state_names=names,muscle_names=muscles,target_mass_kg=mass,source_model_sha256=source_model,target_model_sha256=model)
    with np.load(out/'target.npz',allow_pickle=False) as check:
        if not np.array_equal(check['x0'],x) or not np.array_equal(check['u0'],u):raise ValueError('Export changed arrays')
    metadata={'schema':'ihm.native-copied-equilibrium-lqr-target.v1','target_registration':a.registration,'target_registration_sha256':sha(rp),'copied_state_path':a.copied_state,'copied_state_sha256':sha(cp),'policy_source_path':a.policy,'policy_sha256':sha(pp),'baseline_policy_path':a.baseline_policy,'baseline_policy_sha256':sha(bp),'source_model_sha256':source_model,'target_model_sha256':model,'same_physical_model_except_initial_activation_defaults':True,'state_values_from':'Exact target-policy x0/u0, verified against copied equilibrated native state including all fibers; NOT untouched initializer snapshot','state_count':len(names),'muscle_count':len(muscles),'target_mass_kg':mass,'state_path_order_matches_baseline_and_target_policy':True,'x0_max_absolute_error_against_target_policy':0.,'u0_max_absolute_error_against_target_policy':0.,'x0_max_absolute_error_against_copied_state':0.,'u0_max_absolute_error_against_copied_state':0.,'target_sha256':sha(out/'target.npz'),'physical_time_advanced_s':0,'new_native_process_started':False,'actual_dynamic_swing_demonstrated':False}
    (out/'manifest.json').write_text(json.dumps(metadata,indent=2)+'\n');print(json.dumps(metadata,indent=2))
if __name__=='__main__':main()
