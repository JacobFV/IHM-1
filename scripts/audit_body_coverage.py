"""Record artifact and execution coverage without confusing it with validation."""
from pathlib import Path
import argparse
import collections
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.body import write_json
from ihm.app.experiments import read_experiment
from ihm.human import ImplicitHuman


def audit(root=ROOT):
    root = Path(root)
    directory = root/'data/derived/canonical'
    artifacts = {}
    for path in sorted(directory.glob('*.json')):
        raw = path.read_bytes()
        artifacts[path.stem] = {'path': str(path.relative_to(root)),
                                'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw)}
    anatomy = json.loads((directory/'anatomy.json').read_text())
    trajectory = json.loads((directory/'trajectory.json').read_text()) if (directory/'trajectory.json').exists() else {}
    frames = trajectory.get('frames', [])
    moving = sorted(set().union(*(set(f.get('entities', {})) for f in frames)))
    manifest = json.loads((directory/'body.json').read_text())
    stale = []
    for kind in ['sources', 'runtime_sources']:
        for name, record in manifest.get(kind, {}).items():
            path = root/record['path']
            if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
                stale.append(f'{kind}:{name}')
    if trajectory.get('runtime_sources') != manifest.get('runtime_sources'):
        stale.append('trajectory:runtime_sources')
    for name, record in manifest.get('sources', {}).items():
        if trajectory.get('sources', {}).get(name) != record:
            stale.append('trajectory:source:'+name)
    channels = set().union(*(set(f.get('physiology', {})) for f in frames))
    entries = []
    def add(id, owner, path, evidence, field=None, channel=None, limitation=''):
        integrated = bool(frames) and (all(field in f for f in frames) if field else channel in channels if channel else False)
        entries.append(dict(id=id, state_owner=owner, implementation_path=path,
            implemented=(root/path).exists(), recorded_integration=integrated,
            recorded_evidence=evidence if integrated else None,
            numerical_verification='see per-experiment reports; presence is not a passing test',
            empirical_constraint='source-model evidence only; inspect parameter cards',
            empirically_validated_whole_body=False, limitation=limitation))
    add('coarse_mechanics','canonical mechanics','ihm/assembly/mechanics.py','frames[].mechanical_audit',field='mechanical_audit',limitation='Affine solids and constrained bone orientations; no local tissue contact or bed equilibrium.')
    add('thoracic_kinematics','thoracic mechanics','ihm/assembly/respiration.py','frames[].respiration',field='respiration',limitation='Three-mode reduction; prescribed native gas in replay mode, no full pleural contact.')
    add('somatic_peripheral','somatic peripheral','ihm/assembly/peripheral.py','frames[].peripheral',field='peripheral',limitation='Named inferred pathways; target inventory is not complete anatomical innervation.')
    add('reduced_brain','IBM-derived reduction','ihm/assembly/brain.py','frames[].brain',field='brain',limitation='Selected neural laws and 80 populations do not establish full IBM runtime reuse.')
    for id, channel in [('cardiovascular','MeanArterialPressure(mmHg)'),('respiratory','RespirationRate(1/min)'),('blood','OxygenSaturation'),('renal','UrineProductionRate(mL/min)'),('thermal','CoreTemperature(degC)')]:
        add(id,'BioGears','scripts/native_biogears_rest.cpp',channel,channel=channel,limitation='Native compartment-level physiology; geometry does not establish spatial state resolution.')
    for id, path, owner, experiment in [
        ('microvascular','ihm/assembly/vascular.py','regional passive pressure/flow solve','details'),
        ('lymph_interstitial','ihm/assembly/skin_transport.py','isolated regional fluid/albumin reservoirs','skin-transport'),
        ('body_hair','ihm/assembly/hair.py','material attachments transported by computed skin field','details'),
        ('skin_electrical','ihm/assembly/skin_bioelectric.py','regional non-neural capacitor states','skin-electric'),
        ('reference_contact','ihm/assembly/mechanics_backend.py','regional quasistatic tetrahedral solve','forearm-touch'),
        ('ibm_materializer','ihm/brain/ibm_backend.py','pinned selected IBM transfers','forearm-touch')]:
        add(id,owner,path,'frames[].'+id,field=id,limitation='Separate regional/source evidence is available; it does not establish global state coupling or complete anatomical coverage.')
        try:
            data=read_experiment(root,experiment)
            entries[-1]['regional_evidence']={'experiment':experiment,'source_receipts_current':True,
                'limitations':data.get('limitations',[]),'global_native_storage_connected':False}
        except (OSError,ValueError,KeyError) as error:
            entries[-1]['regional_evidence']={'experiment':experiment,'source_receipts_current':False,'error':str(error)}
    for id in ['endocrine','gastrointestinal_hepatic','immune_hematologic','reproductive','special_senses','autonomic_visceral','connective_adipose_marrow']:
        add(id,'BioGears where implemented; detailed spatial ownership unresolved','scripts/native_biogears_rest.cpp',None,limitation='Requires a system-specific state/perturbation audit; native code presence alone is insufficient.')
    references={}
    for kind in ('compression','contact'):
        path=root/'data/derived/mechanics-reference'/kind/'benchmark.json'
        if path.exists():
            raw=path.read_bytes();data=json.loads(raw)
            references[kind]={'path':str(path.relative_to(root)),'sha256':hashlib.sha256(raw).hexdigest(),
                'recorded_status':data.get('status'),'scope':data.get('scope'),
                'interpretation':'Native benchmark receipt; this inventory does not rerun or revalidate the reference solver.'}
    return dict(schema_version=2, created_utc=datetime.now(timezone.utc).isoformat(),
        revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
        workspace_changes=subprocess.check_output(['git','status','--short'],cwd=root,text=True).splitlines(),
        artifacts=artifacts, anatomy={'representation_count':len(anatomy['entities']), 'roles':dict(collections.Counter(e['role'] for e in anatomy['entities']))},
        execution={'frames':len(frames),'moving_entity_count':len(moving),'moving_entity_ids':moving,'physiology_channels':sorted(channels), 'stale_dependencies':stale},
        capabilities=entries,reference_mechanics=references,microstructure_evidence=ImplicitHuman.open(root).microstructure_evidence(),
        interpretation='Canonical recorded integration is established only by frame fields/channels. Regional artifacts, source acquisitions and native reference benchmarks have separate receipts. This audit does not assert calibrated mechanisms or empirical validity.')


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'data/derived/audits/execution-coverage.json');a=p.parse_args()
    result=audit();write_json(a.output,result)
    print(json.dumps({'output':str(a.output),'execution':{k:v for k,v in result['execution'].items() if k!='moving_entity_ids'}},indent=2))
