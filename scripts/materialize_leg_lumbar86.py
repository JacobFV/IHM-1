#!/usr/bin/env python3
"""Opt-in 86-muscle controller-search plant: source80 plus retained lumbar6.

Retains all source bodies/inertias and joints. Arm26 muscles/wraps are not added;
this is an explicitly different muscle authority from the 98-muscle model.
"""
from pathlib import Path
import hashlib,json,sys
import xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.embodied import _prepare_mechanical_registration
ROOT=Path(__file__).resolve().parents[1]
REG='data/derived/mechanics/whole_body_lumbar_current/registration.json'
SOURCE='data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled.osim'
OUT=ROOT/'data/derived/mechanics/leg_lumbar86'
sha=lambda raw:hashlib.sha256(raw).hexdigest()

def main():
    if OUT.exists():raise ValueError('Refusing to replace existing derived candidate')
    _,lumbar,rows=_prepare_mechanical_registration(ROOT,REG)
    raw=(ROOT/SOURCE).read_bytes()
    if sha(raw)!=lumbar['sources'][SOURCE]:raise ValueError('Source80 identity mismatch')
    original=ET.fromstring(raw)
    reference=ET.parse(ROOT/lumbar['model_path'])
    a=next(x for x in original.findall('.//JointSet/objects/*') if x.get('name')=='back')
    b=next(x for x in reference.findall('.//JointSet/objects/*') if x.get('name')=='back')
    if ET.tostring(a)!=ET.tostring(b):raise ValueError('Lumbar target frame differs from source80')
    muscles=[m for m in original.findall('.//ForceSet/objects/*') if 'Muscle' in m.tag]
    if len(muscles)!=80 or [m.get('name') for m in muscles]!=[row['id'] for row in rows[:80]]:raise ValueError('Original80 catalog mismatch')
    insert=(ROOT/lumbar['insert_path']).read_bytes()
    if sha(insert)!=lumbar['insert_sha256']:raise ValueError('Lumbar insertion changed')
    start=raw.index(b'<ForceSet');stop=raw.index(b'</ForceSet>',start);end=raw.rindex(b'</objects>',start,stop)
    model=raw[:end]+insert+raw[end:]
    parsed=ET.fromstring(model)
    names=[m.get('name') for m in parsed.findall('.//ForceSet/objects/*') if 'Muscle' in m.tag]
    catalog=rows[:80]+rows[92:]
    if names!=[r['id'] for r in catalog] or len(names)!=86:raise ValueError('Candidate muscle catalog mismatch')
    for xpath in ('.//BodySet','.//JointSet'):
        if ET.tostring(original.find(xpath))!=ET.tostring(parsed.find(xpath)):raise ValueError('Source body or joint changed')
    OUT.mkdir(parents=True)
    def write(name,data):
        path=OUT/name;path.write_bytes(data);return str(path.relative_to(ROOT)),sha(data)
    mp,mh=write('subject_leg_lumbar86.osim',model)
    cp,ch=write('catalog.json',(json.dumps(catalog,indent=2)+'\n').encode())
    ip,ih=write('insert.xml',insert)
    sources={SOURCE:sha(raw),REG:sha((ROOT/REG).read_bytes()),lumbar['model_path']:lumbar['model_sha256'],lumbar['catalog_path']:lumbar['catalog_sha256'],lumbar['insert_path']:lumbar['insert_sha256'],str(Path(__file__).relative_to(ROOT)):sha(Path(__file__).read_bytes())}
    for row in catalog:
        source=row['source_path'];data=(ROOT/source).read_bytes()
        if sha(data)!=row['source_sha256']:raise ValueError('Catalog donor identity mismatch')
        sources[source]=sha(data)
    provenance={'schema':'ihm.leg-lumbar86-materialization.v1','source80_exact_bytes_preserved_outside_insertion':True,'insertion_matches_retained_lumbar6_exact_bytes':True,'back_joint_matches_lumbar_registration_exactly':True,'body_and_joint_sets_unchanged':True,'full_source_body_inertia_retained':True,'native_dynamic_acceptance':False,'arm26_muscles_added':False,'scope':'Opt-in engineering stance-controller search plant with80 leg muscles and6 transferred lumbar muscles; no active arm muscle authority or walking evidence.'}
    pp,ph=write('provenance.json',(json.dumps(provenance,indent=2)+'\n').encode());sources[pp]=ph
    registration={'schema':'ihm.lumbar-muscle-variant.v1','model_path':mp,'model_sha256':mh,'catalog_path':cp,'catalog_sha256':ch,'insert_path':ip,'insert_sha256':ih,'base_model_path':SOURCE,'sources':sources,'muscle_count':86,'body_count':len(parsed.findall('.//BodySet/objects/Body')),'default_enabled':False,'native_acceptance_complete':False,'registration':lumbar['registration'],'limits':['Experimental controller-search plant;98-muscle transfer must be evaluated independently.','Arm26 muscles and wrapping geometry are not present; retained arm bodies remain inertial and passively articulated.','Lumbar donor parameters are transferred priors, not subject calibration.'],'provenance_path':pp}
    rp,_=write('registration.json',(json.dumps(registration,indent=2)+'\n').encode())
    _,_,validated=_prepare_mechanical_registration(ROOT,rp)
    print(json.dumps({'registration':rp,'muscles':len(validated),'bodies':registration['body_count'],'preflight_passed':True},indent=2))
if __name__=='__main__':main()
