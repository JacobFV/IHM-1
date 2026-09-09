"""Materialize diagnostic native-selected activation defaults, preserving other XML bytes."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import tempfile
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(pose_path):
    pose_path=pose_path.resolve(); pose=json.loads(pose_path.read_text())
    source=ROOT/'data/derived/mechanics/whole_body_lumbar_current/registration.json'
    registration=json.loads(source.read_text());model_path=ROOT/registration['model_path']
    if pose['model_sha256']!=sha(model_path):raise ValueError('Solved pose model differs from source registration')
    selection_path=(ROOT/pose['activation_selection_path']).resolve()
    if not selection_path.is_relative_to(ROOT):raise ValueError('Owned native activation selection receipt required')
    selection=json.loads(selection_path.read_text());expected=set(selection['actuators'])
    lumbar={f'gait2392_{muscle}_{side}' for side in ('r','l') for muscle in ('ercspn','intobl','extobl')}
    if not lumbar.issubset(expected):raise ValueError('Six lumbar tones must remain in selected set')
    activations=pose['activations']
    if set(activations)!=expected or any(not .01<=v<=1 for v in activations.values()):raise ValueError('Exactly the native-selected bounded activation identities required')
    original=model_path.read_text();updated=original
    for name,value in activations.items():
        pattern=rf'(<(?:Thelen2003Muscle|Millard2012EquilibriumMuscle) name="{re.escape(name)}">)(.*?)(</(?:Thelen2003Muscle|Millard2012EquilibriumMuscle)>)'
        def replace(match):
            body=match[2]
            if '<default_activation>' in body:
                body,n=re.subn(r'<default_activation>[^<]*</default_activation>',f'<default_activation>{value:.17g}</default_activation>',body)
                if n!=1:raise ValueError('Ambiguous default activation')
            else:body+='\n        <default_activation>'+format(value,'.17g')+'</default_activation>\n      '
            return match[1]+body+match[3]
        updated,n=re.subn(pattern,replace,updated,flags=re.S)
        if n!=1:raise ValueError('Missing/ambiguous selected muscle')
    # Compare parsed trees after removing only the explicitly changed property.
    old_tree=ET.fromstring(original);new_tree=ET.fromstring(updated)
    for tree in (old_tree,new_tree):
        for muscle in tree.iter():
            if muscle.tag in ('Thelen2003Muscle','Millard2012EquilibriumMuscle') and muscle.get('name') in expected:
                for child in list(muscle):
                    if child.tag=='default_activation':muscle.remove(child)
        for element in tree.iter():
            element.text=(element.text or '').strip();element.tail=(element.tail or '').strip()
    if ET.tostring(old_tree)!=ET.tostring(new_tree):raise AssertionError('Non-activation model change')
    output=Path(tempfile.mkdtemp(prefix='supine-selected-defaults-',dir=ROOT/'data/derived'))
    out_model=output/'subject_with_lumbar.osim';out_model.write_text(updated)
    registration.update(model_path=str(out_model.relative_to(ROOT)),model_sha256=sha(out_model),default_enabled=False,native_acceptance_complete=False,
        activation_initialization=dict(activations=activations,pose_path=str(pose_path.relative_to(ROOT)),pose_sha256=sha(pose_path),scope='Diagnostic native-selected hip/lumbar defaults only; held excitation and coupled runtime require independent validation.'))
    registration['limits']=[line for line in registration.get('limits',[]) if not line.startswith('No excitation assignment')] + ['Hip/lumbar activation defaults selected by native static optimization; this is a candidate initialization, not a trained controller or equilibrium certification.']
    registration['sources'].update({str(source.relative_to(ROOT)):sha(source),str(model_path.relative_to(ROOT)):sha(model_path),str(pose_path.relative_to(ROOT)):sha(pose_path),str(Path(__file__).resolve().relative_to(ROOT)):sha(__file__)})
    registration['sources'][str(selection_path.relative_to(ROOT))]=sha(selection_path)
    arms=selection_path.parent/'native_moment_arms.json'
    registration['sources'][str(arms.relative_to(ROOT))]=sha(arms)
    (output/'registration.json').write_text(json.dumps(registration,indent=2)+'\n')
    variant_pose={**pose,'coordinates':dict(sorted(pose['coordinates'].items())),'source_model_sha256':pose['model_sha256'],'model_sha256':sha(out_model),'registration_path':str((output/'registration.json').relative_to(ROOT))}
    (output/'initial_pose.json').write_text(json.dumps(variant_pose,indent=2)+'\n')
    print(json.dumps(dict(output=str(output),registration=str(output/'registration.json'),initial_pose=str(output/'initial_pose.json')),indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('pose',type=Path);args=parser.parse_args();main(args.pose)
