#!/usr/bin/env python3
"""Transfer solved initial defaults to98 muscles and verify native static contact.

This creates an opt-in initialization artifact; it does not promote runtime
selection defaults or claim controller transfer/dynamic stability.
"""
from pathlib import Path
import hashlib,json,re,sys
import xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.embodied import _prepare_mechanical_registration
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.moment_arm_control import native_center_of_mass,foot_support_center
ROOT=Path(__file__).resolve().parents[1]
BASE='data/derived/mechanics/whole_body_lumbar_current/registration.json'
SOLVED='data/derived/mechanics/initial_stance86/registration.json'
OUT=ROOT/'data/derived/mechanics/initial_stance98'
sha=lambda raw:hashlib.sha256(raw).hexdigest()

def main():
    if OUT.exists():raise ValueError('Refusing to replace retained candidate')
    _,base,_=_prepare_mechanical_registration(ROOT,BASE)
    _,solved,_=_prepare_mechanical_registration(ROOT,SOLVED)
    receipt=json.loads((ROOT/solved['static_solve_receipt']).read_text())
    if max(abs(v) for v in receipt['residual'])>.005:raise ValueError('Source initial solve not acceptable')
    coordinates=solved['initial_pose_only_coordinates'];raw=(ROOT/base['model_path']).read_bytes();model=raw
    for name,value in coordinates.items():
        pattern=rb'(<Coordinate name="'+name.encode()+rb'">.*?<default_value>)(.*?)(</default_value>)'
        model,count=re.subn(pattern,lambda m:m[1]+format(value,'.17g').encode()+m[3],model,count=1,flags=re.S)
        if count!=1:raise ValueError('Coordinate default not found')
    # Ensure transfer has the same source bodies/mass geometry and every joint.
    transferred=ET.fromstring(model);reference=ET.parse(ROOT/solved['model_path'])
    for xpath in ('.//JointSet',):
        if ET.tostring(transferred.find(xpath))!=ET.tostring(reference.find(xpath)):raise ValueError('Initial joint geometry differs')
    for b in transferred.findall('.//BodySet/objects/Body'):
        ref=next(x for x in reference.findall('.//BodySet/objects/Body') if x.get('name')==b.get('name'))
        for field in ('mass','mass_center','inertia'):
            if b.findtext(field)!=ref.findtext(field):raise ValueError('Effective initial mass geometry differs')
    OUT.mkdir(parents=True)
    path=OUT/'model.osim';path.write_bytes(model)
    registration=dict(base);registration['sources']=dict(base['sources'])
    registration['sources'].update({BASE:sha((ROOT/BASE).read_bytes()),SOLVED:sha((ROOT/SOLVED).read_bytes()),base['model_path']:base['model_sha256'],solved['static_solve_receipt']:sha((ROOT/solved['static_solve_receipt']).read_bytes()),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__).read_bytes())})
    registration['model_path']=str(path.relative_to(ROOT));registration['model_sha256']=sha(model);registration['initial_pose_only_coordinates']=coordinates
    registration['native_acceptance_complete']=False
    regpath=OUT/'registration.json';regpath.write_text(json.dumps(registration,indent=2)+'\n')
    native=NativeMechanicalStream(ROOT,OUT/'native_initialization',environment='upright',target_mass_kg=70,augmented_registration=str(regpath.relative_to(ROOT)))
    try:state=native.snapshot()
    finally:native.close()
    com,velocity=native_center_of_mass(state);support=foot_support_center(state)
    force=sum(c['force_n'][1] for c in state['contacts'] if c.get('body_frame') in ('calcn_r','calcn_l','toes_r','toes_l'))
    if support is None or abs(com[0]-support[0])>5e-4 or abs(force-686.7)>.5 or len(state['muscles'])!=98:raise ValueError('Native98 initial support mismatch')
    snapshot=OUT/'initial_snapshot.json';snapshot.write_text(json.dumps(state,allow_nan=False))
    result={'schema':'ihm.initial-stance98-transfer.v1','source_solve':SOLVED,'only_default_coordinate_changes':coordinates,'muscles':98,'body_mass_kg':state['mass_kg'],'initial_com_ground_m':com.tolist(),'initial_support_center_ground_m':support.tolist(),'initial_foot_vertical_force_n':force,'native_initialization_accepted':True,'native_muscle_equilibrium':True,'physical_time_advanced_s':0,'dynamic_stability_demonstrated':False,'external_support':False,'prescribed_motion':False,'default_runtime_selection_changed':False,'initial_snapshot_sha256':sha(snapshot.read_bytes())}
    report=OUT/'initial_acceptance.json';report.write_text(json.dumps(result,indent=2)+'\n')
    registration['sources'][str(report.relative_to(ROOT))]=sha(report.read_bytes());registration['initial_acceptance_receipt']=str(report.relative_to(ROOT))
    regpath.write_text(json.dumps(registration,indent=2)+'\n')
    _prepare_mechanical_registration(ROOT,str(regpath.relative_to(ROOT)))
    print(json.dumps({'registration':str(regpath.relative_to(ROOT)),**result},indent=2))
if __name__=='__main__':main()
