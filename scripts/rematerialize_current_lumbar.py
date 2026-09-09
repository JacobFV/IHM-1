#!/usr/bin/env python3
"""Re-run lumbar source audit and retain a mechanically byte-identical variant.

Canonical anatomy/mechanics are coverage metadata in the audit, not inputs to the
six donor muscle path transformations. Rebuilding and comparing the insertion
fragment proves this for the current inputs rather than merely changing hashes.
"""
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.audit_lumbar_shoulder_coverage import build
from ihm.assembly.embodied import _prepare_mechanical_registration
ROOT=Path(__file__).resolve().parents[1]
OLD='data/derived/lumbar-muscle-native-lb45uirs/variant/registration.json'
OUT=ROOT/'data/derived/mechanics/whole_body_lumbar_current'
sha=lambda raw:hashlib.sha256(raw).hexdigest()

def main():
    if OUT.exists():raise ValueError('Choose fresh output; existing materialization retained')
    original=json.loads((ROOT/OLD).read_text())
    old_audit='data/research/lumbar_shoulder_coverage/audit.json'
    old_fragment='data/research/lumbar_shoulder_coverage/lumbar_candidate_forces.xml'
    coverage={'data/derived/canonical/anatomy.json','data/derived/canonical/mechanics.json'}
    changed={}
    for path,digest in original['sources'].items():
        current=sha((ROOT/path).read_bytes())
        if current!=digest:
            if path not in coverage:raise ValueError('Mechanical donor or retained audit changed: '+path)
            changed[path]={'previous_sha256':digest,'current_sha256':current}
    report,fragment=build()
    if fragment!=(ROOT/old_fragment).read_bytes():raise ValueError('Recomputed donor insertion changed mechanically')
    retained_audit=json.loads((ROOT/old_audit).read_text())
    for key in ('lumbar_candidate','registration','arm_paths','fitted_paths'):
        if report[key]!=retained_audit[key]:raise ValueError('Recomputed mechanical audit changed: '+key)
    for field,digest in [('model_path','model_sha256'),('catalog_path','catalog_sha256'),('insert_path','insert_sha256')]:
        if sha((ROOT/original[field]).read_bytes())!=original[digest]:raise ValueError('Retained output changed: '+field)
    OUT.mkdir(parents=True)
    def write(name,raw):
        path=OUT/name;path.write_bytes(raw);return str(path.relative_to(ROOT)),sha(raw)
    audit_path,audit_sha=write('audit.json',(json.dumps(report,indent=2)+'\n').encode())
    fragment_path,fragment_sha=write('lumbar_candidate_forces.xml',fragment)
    registration=dict(original);registration['sources']=dict(original['sources'])
    for field,digest,name in [('model_path','model_sha256','subject_with_lumbar.osim'),('catalog_path','catalog_sha256','catalog.json'),('insert_path','insert_sha256','insert.xml')]:
        registration[field],registration[digest]=write(name,(ROOT/original[field]).read_bytes())
    for path in (old_audit,old_fragment,*coverage):registration['sources'].pop(path,None)
    registration['sources'].update(report['sources'])
    registration['sources'].update({audit_path:audit_sha,fragment_path:fragment_sha,OLD:sha((ROOT/OLD).read_bytes()),'scripts/audit_lumbar_shoulder_coverage.py':sha((ROOT/'scripts/audit_lumbar_shoulder_coverage.py').read_bytes()),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__).read_bytes())})
    provenance={'schema':'ihm.lumbar-source-rematerialization.v1','previous_registration':OLD,'changed_coverage_sources':changed,'current_audit_recomputed':True,'donor_fragment_byte_identical':True,'mechanical_audit_fields_equal':['lumbar_candidate','registration','arm_paths','fitted_paths'],'model_catalog_and_insertion_byte_identical':True,'new_dynamic_acceptance_claimed':False,'explanation':'Canonical files affect coverage metadata only; current audit rebuild reproduces exactly the prior donor path fragment and mechanical registration.'}
    receipt_path,receipt_sha=write('rematerialization.json',(json.dumps(provenance,indent=2)+'\n').encode())
    registration['sources'][receipt_path]=receipt_sha
    registration['rematerialization_receipt']=receipt_path
    path,_=write('registration.json',(json.dumps(registration,indent=2)+'\n').encode())
    frozen,validated,rows=_prepare_mechanical_registration(ROOT,path)
    print(json.dumps({'registration':path,'preflight_passed':True,'muscles':len(rows),'model_sha256':validated['model_sha256'],'changed_metadata':list(changed)},indent=2))
if __name__=='__main__':main()
