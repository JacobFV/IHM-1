#!/usr/bin/env python3
"""Export retained copied equilibrium, translated to its declared initial anchor.

The original query already equilibrated all fiber states. A horizontal global
translation on the retained flat ground changes neither gravity nor contact gaps.
No native dynamics or runtime state is changed by this export.
"""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--registration',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    rp=ROOT/a.registration;reg=json.loads(rp.read_text());reportpath=ROOT/reg['static_resting_solve'];report=json.loads(reportpath.read_text());states=dict(report['candidate_state_variables'])
    shifts={}
    for coordinate in ('pelvis_tx','pelvis_tz'):
        paths=[n for n in states if n.endswith('/'+coordinate+'/value')]
        if len(paths)!=1:raise ValueError('Missing root translation state')
        n=paths[0];old=states[n];states[n]=reg['initial_pose'].get(coordinate,old);shifts[n]={'copied_query_value':old,'anchored_value':states[n],'translation_m':states[n]-old}
    out=ROOT/a.output
    if out.exists():raise ValueError('Refusing to overwrite evidence')
    payload={'schema':'ihm.native-copied-equilibrium-state.v1','state_variables':states,'source_registration':a.registration,'source_registration_sha256':sha(rp),'source_model_path':reg['model_path'],'source_model_sha256':reg['model_sha256'],'target_mass_kg':reg['target_mass_kg'],'source_static_query_report':str(reportpath.relative_to(ROOT)),'source_static_query_report_sha256':sha(reportpath),'state_semantics':'Copied native State after muscle equilibration; NOT untouched initializer snapshot','horizontal_translation_only':shifts,'flat_ground_translation_invariance_assumed':True,'equilibrium_excitations':reg['equilibrium_excitations'],'copied_query_all_independent_residual':report['all_independent_residual'],'physical_time_advanced_s':0,'dynamic_transition_demonstrated':False}
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(payload,indent=2)+'\n');print(json.dumps({'output':str(out.relative_to(ROOT)),'state_count':len(states),'translations':shifts},indent=2))
if __name__=='__main__':main()
