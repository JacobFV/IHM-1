"""Opt-in exact-byte extension of the retained 92-muscle model; no native job."""
from pathlib import Path
import argparse,hashlib,json
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]

def materialize(root,output):
    root=Path(root).resolve();output=Path(output).resolve()
    if not output.is_relative_to(root/'data/derived'):raise ValueError('Owned derived output required')
    if output.exists() and any(output.iterdir()):raise ValueError('Refusing to overwrite materialized variant')
    inputs={}
    def read(p):
        raw=(root/p).read_bytes();inputs[p]=hashlib.sha256(raw).hexdigest();return raw
    base_reg=json.loads(read('data/derived/mechanics/whole_body_arm26_v2/registration.json'))
    base=read(base_reg['model_path']);catalog=json.loads(read(base_reg['catalog_path']))
    if inputs[base_reg['model_path']]!=base_reg['model_sha256'] or inputs[base_reg['catalog_path']]!=base_reg['catalog_sha256']:raise ValueError('Base identity mismatch')
    for p,h in base_reg['sources'].items():
        read(p)
        if inputs[p]!=h:raise ValueError('Base donor identity mismatch')
    audit=json.loads(read('data/research/lumbar_shoulder_coverage/audit.json'))
    fragment=read('data/research/lumbar_shoulder_coverage/lumbar_candidate_forces.xml')
    if hashlib.sha256(fragment).hexdigest()!=audit['candidate_xml_sha256']:raise ValueError('Fragment identity mismatch')
    for p,h in audit['sources'].items():
        read(p)
        if inputs[p]!=h:raise ValueError('Audit source identity mismatch')
    forces=ET.fromstring(fragment).findall('.//ForceSet/objects/*')
    original=ET.fromstring(base);old=[m for m in original.findall('.//ForceSet/objects/*') if 'Muscle' in m.tag]
    if len(old)!=92 or [m.get('name') for m in old]!=[m['id'] for m in catalog]:raise ValueError('Base muscle catalog order differs')
    if len(forces)!=6:raise ValueError('Expected exact six lumbar muscles')
    for f,row in zip(forces,audit['lumbar_candidate']):
        if f.get('name')!=row['candidate_name']:raise ValueError('Candidate order differs')
        side=row['source_name'][-1];hemi='lh' if side=='r' else 'rh'
        catalog.append(dict(id=f.get('name'),side=side,body_group='lumbar_trunk',attachment_bodies=['pelvis','torso'],
            sensory_region=f'brain-{hemi}-postcentral',motor_region=f'brain-{hemi}-precentral',
            max_isometric_force_n=row['parameters']['max_isometric_force'],optimal_fiber_length_m=row['parameters']['optimal_fiber_length'],
            tendon_slack_length_m=row['parameters']['tendon_slack_length'],source_path='data/raw/anatomy/opensim-models/source/Models/Gait2392_Simbody/gait2392_thelen2003muscle.osim',
            source_muscle_name=row['source_name'],source_muscle_law=row['type'],path_points=row['stations'],
            assignment_basis='Contralateral regional cortical engineering prior only; no measured trunk recruitment or excitation policy.',
            parameter_basis='Unscaled retained donor parameters; unit-scale homologous back-frame station transfer, anatomical correspondence unvalidated.',
            native_control_ready=False,default_excitation_assignment=None))
        catalog[-1]['source_sha256']=inputs[catalog[-1]['source_path']]
    # Insert immediately before the closing ForceSet objects tag. Everything
    # already present, including comments, whitespace and floating literals, survives.
    start=base.index(b'<ForceSet');stop=base.index(b'</ForceSet>',start);end=base.rindex(b'</objects>',start,stop)
    insert=b'\n'+b'\n'.join(ET.tostring(f,encoding='utf-8') for f in forces)+b'\n'
    model=base[:end]+insert+base[end:]
    output.mkdir(parents=True,exist_ok=True)
    def write(name,raw):
        p=output/name;p.write_bytes(raw);return str(p.relative_to(root)),hashlib.sha256(raw).hexdigest()
    mp,mh=write('subject_with_lumbar.osim',model);cp,ch=write('catalog.json',(json.dumps(catalog,indent=2)+'\n').encode());ip,ih=write('insert.xml',insert)
    r=dict(schema='ihm.lumbar-muscle-variant.v1',model_path=mp,model_sha256=mh,catalog_path=cp,catalog_sha256=ch,
           insert_path=ip,insert_sha256=ih,base_model_path=base_reg['model_path'],sources=inputs,
           muscle_count=98,body_count=22,default_enabled=False,native_acceptance_complete=False,
           registration=audit['registration'],limits=audit['limits'],materializer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    write('registration.json',(json.dumps(r,indent=2)+'\n').encode());return r
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='data/derived/mechanics/whole_body_lumbar_v1');a=p.parse_args()
    print(json.dumps(materialize(ROOT,ROOT/a.output),indent=2))
