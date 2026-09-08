"""One explicit generic demographic profile shared by anatomy and physiology."""
from pathlib import Path
import json,hashlib
import xml.etree.ElementTree as ET

# Body mass is measured on this specimen's own geometry, not inherited. The
# BioGears StandardMale 77.1107029 kg constant cannot be reconciled with the
# acquired envelope: the composed fill weighs 31.7461 kg over the 0.032117 m3
# void, a mean 988.4 kg/m3, and reaching 77.1107029 kg would need the void to
# weigh 1185.8 kg/m3 -- denser than every soft tissue in ICRU-44 -- which by the
# Siri relation is a negative fat fraction, i.e. not a body. The two free
# choices in that ledger move the answer by 0.024 kg (8 mm vs 10 mm partition)
# and 0.534 kg (residual density swept over the whole [1000, 1060] band). Only
# the height in the inherited profile ever came from the atlas.
MASS_LEDGER='data/derived/interstitial-composition-prior-v1/ledger.json'
MASS_KG=70.7713

def build_profile(root):
    root=Path(root)
    # The declared mass is gated on the artifact it came from before anything
    # else is read, so the number can never drift away from its own evidence.
    ledger_path=root/MASS_LEDGER;ledger=json.loads(ledger_path.read_text())
    composed=ledger['verdict']['composed_total_body_mass_kg']
    if abs(composed-MASS_KG)>5e-5:
        raise ValueError('Declared body mass disagrees with its own composition ledger: %r vs %r'%(MASS_KG,composed))
    raw=root/'data/derived/anatomy/bodyparts3d_index.json'
    skin=next(e for e in json.loads(raw.read_text())['meshes'] if e['element_id']=='FJ2810')
    height=(skin['bounds_in_source_coordinates'][1][2]-skin['bounds_in_source_coordinates'][0][2])*.001
    directory=root/'data/runtime/physiology/biogears-build/runtime/patients';source=directory/'StandardMale.xml'
    tree=ET.parse(source)
    for node in tree.getroot():
        name=node.tag.split('}')[-1]
        if name=='Name':node.text='IHMGenericMale'
        elif name=='Height':node.set('unit','m');node.set('value',repr(height))
    target=directory/'IHMGenericMale.xml';tree.write(target,encoding='UTF-8',xml_declaration=True)
    profile=dict(id='ihm-generic-male-v1',native_patient='IHMGenericMale',sex='male',height_m=height,mass_kg=MASS_KG,age_years=44,body_fat_fraction=.21,
        reference_posture='supine reference condition; acquired anatomical coordinates retained; no bed-contact equilibrium',
        parameters={'height_m':{'basis':'canonical skin extent','source':skin['source_path'],'source_sha256':skin['sha256']},
                    'mass_kg':{'basis':'composed over this specimen\'s own measured interior: the 8 mm voxel partition of the acquired envelope, filled at sourced per-constituent densities','source':MASS_LEDGER,'source_sha256':hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
                               'composed_total_body_mass_kg':composed,'declared_rounded_to_kg':1e-4,'implied_bmi':composed/height**2,
                               'supersedes':{'value_kg':77.1107029,'basis':'inherited BioGears StandardMale constant','reason':'unreachable in the acquired envelope; requires a 1185.8 kg/m3 void, denser than every ICRU-44 soft tissue, and a negative Siri fat fraction'},
                               'partition_sensitivity_kg':ledger['verdict']['partition_sensitivity_kg'],'residual_density_sensitivity_kg':ledger['verdict']['residual_density_sensitivity_kg'],
                               'calibrated':False,'measured_patient':False},
                    'other_demographics':{'basis':'inherited generic source-model prior','source':str(source.relative_to(root)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}},
        native_patient_path=str(target.relative_to(root)),native_patient_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        limitations=['Generic synthesized reference profile, not a measured patient.','Body mass is composed from the acquired geometry and sourced constituent densities, not weighed; its sensitivity to the partition and to the residual density band is recorded in parameters.mass_kg.','Changing native height triggers source-engine initialization; it does not independently calibrate organ geometry or material properties.','Native environment retains source posture limitations; mechanical bed-rest support is separate.'])
    out=root/'data/derived/canonical';out.mkdir(parents=True,exist_ok=True);(out/'profile.json').write_text(json.dumps(profile,indent=2)+'\n')
    return profile
