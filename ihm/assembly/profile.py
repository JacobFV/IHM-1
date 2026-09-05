"""One explicit generic demographic profile shared by anatomy and physiology."""
from pathlib import Path
import json,hashlib
import xml.etree.ElementTree as ET

def build_profile(root):
    root=Path(root);raw=root/'data/derived/anatomy/bodyparts3d_index.json'
    skin=next(e for e in json.loads(raw.read_text())['meshes'] if e['element_id']=='FJ2810')
    height=(skin['bounds_in_source_coordinates'][1][2]-skin['bounds_in_source_coordinates'][0][2])*.001
    directory=root/'data/runtime/physiology/biogears-build/runtime/patients';source=directory/'StandardMale.xml'
    tree=ET.parse(source)
    for node in tree.getroot():
        name=node.tag.split('}')[-1]
        if name=='Name':node.text='IHMGenericMale'
        elif name=='Height':node.set('unit','m');node.set('value',repr(height))
    target=directory/'IHMGenericMale.xml';tree.write(target,encoding='UTF-8',xml_declaration=True)
    profile=dict(id='ihm-generic-male-v1',native_patient='IHMGenericMale',sex='male',height_m=height,mass_kg=77.1107029,age_years=44,body_fat_fraction=.21,
        reference_posture='supine reference condition; acquired anatomical coordinates retained; no bed-contact equilibrium',
        parameters={'height_m':{'basis':'canonical skin extent','source':skin['source_path'],'source_sha256':skin['sha256']},
                    'other_demographics':{'basis':'inherited generic source-model prior','source':str(source.relative_to(root)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}},
        native_patient_path=str(target.relative_to(root)),native_patient_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        limitations=['Generic synthesized reference profile, not a measured patient.','Changing native height triggers source-engine initialization; it does not independently calibrate organ geometry or material properties.','Native environment retains source posture limitations; mechanical bed-rest support is separate.'])
    out=root/'data/derived/canonical';out.mkdir(parents=True,exist_ok=True);(out/'profile.json').write_text(json.dumps(profile,indent=2)+'\n')
    return profile
