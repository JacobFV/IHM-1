#!/usr/bin/env python3
"""Recenter an INITIAL pose on one physical stance foot; never advances time."""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.assembly.embodied import _prepare_mechanical_registration
ROOT=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--registration',required=True);p.add_argument('--anchor-registration',required=True);p.add_argument('--output',required=True);p.add_argument('--side',choices=('l','r'),default='l');a=p.parse_args()
    source=ROOT/a.registration;anchor=ROOT/a.anchor_registration;out=ROOT/a.output
    if out.exists():raise ValueError('Refusing to overwrite artifact')
    reg=json.loads(source.read_text());anchorreg=json.loads(anchor.read_text())
    if reg['target_mass_kg']!=anchorreg['target_mass_kg']:raise ValueError('Mass mismatch')
    current=json.loads((source.parent/'initial_snapshot.json').read_text());reference=json.loads((anchor.parent/'initial_snapshot.json').read_text())
    def centroid(s):
        contacts=[c['center_m'] for c in s['contacts'] if c['name'].startswith('contact') and c['name'].endswith('_'+a.side)]
        if len(contacts)!=6:raise ValueError('Expected six source foot contacts')
        return np.mean(contacts,axis=0)
    delta=centroid(reference)-centroid(current);delta[1]=0
    q=dict(reg['initial_pose']);q['pelvis_tx']=q.get('pelvis_tx',0)+float(delta[0]);q['pelvis_tz']=q.get('pelvis_tz',0)+float(delta[2])
    reg.update(initial_pose=q,initial_pose_only_coordinates=q)
    reg['sources']=dict(reg['sources']);reg['sources'].update({a.registration:sha(source),a.anchor_registration:sha(anchor),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))})
    out.mkdir(parents=True);rp=out/'registration.json';rp.write_text(json.dumps(reg,indent=2)+'\n')
    native=NativeMechanicalStream(ROOT,out/'native_initialization',environment='upright',target_mass_kg=reg['target_mass_kg'],augmented_registration=str(rp.relative_to(ROOT)),initial_pose=q)
    try:
        state=native.snapshot();r=native.evaluate_static_pose(q)
        error=centroid(state)[[0,2]]-centroid(reference)[[0,2]]
        if np.max(np.abs(error))>1e-10:raise ValueError('Stance ground station not preserved')
        forces={s:sum(c['force_n'][1] for c in state['contacts'] if c['name'].startswith('contact') and c['name'].endswith('_'+s)) for s in ('r','l')}
        report={'schema':'ihm.initial-stance-foot-anchor.v1','source_registration':a.registration,'anchor_registration':a.anchor_registration,'stance_side':a.side,'initial_horizontal_translation_m':delta.tolist(),'stance_contact_centroid_m':centroid(state).tolist(),'anchor_contact_centroid_m':centroid(reference).tolist(),'horizontal_centroid_error_m':error.tolist(),'right_minimum_contact_clearance_m':min(c['center_m'][1]-c['radius_m'] for c in state['contacts'] if c['name'].startswith('contact') and c['name'].endswith('_r')),'foot_vertical_force_n':forces,'maximum_absolute_coordinate_acceleration':max(abs(c['acceleration']) for c in r['coordinates'].values()),'physical_time_advanced_s':0,'prescribed_motion':False,'dynamic_transition_demonstrated':False,'initial_pose':q}
        (out/'initial_snapshot.json').write_text(json.dumps(state));(out/'initial_pose.json').write_text(json.dumps(q,indent=2)+'\n');(out/'equilibrium_excitations.json').write_text(json.dumps(reg['equilibrium_excitations'],indent=2)+'\n')
        reportpath=out/'native_initial_acceptance.json';reportpath.write_text(json.dumps(report,indent=2)+'\n');reg['native_initial_acceptance']=str(reportpath.relative_to(ROOT));reg['sources'][str(reportpath.relative_to(ROOT))]=sha(reportpath);rp.write_text(json.dumps(reg,indent=2)+'\n');_prepare_mechanical_registration(ROOT,str(rp.relative_to(ROOT)))
        print(json.dumps(report,indent=2))
    finally:native.close()
if __name__=='__main__':main()
