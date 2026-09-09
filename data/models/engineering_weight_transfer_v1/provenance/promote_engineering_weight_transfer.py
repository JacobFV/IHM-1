#!/usr/bin/env python3
"""Freeze the demonstrated scheduled weight-transfer reference program.

Runtime remains the existing baseline stance model. Endpoint model XML is retained
only to verify that its physical definition differs solely in activation defaults.
"""
import argparse,gzip,hashlib,io,json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def sha(raw):return hashlib.sha256(raw).hexdigest()
def canonical(node):
    return (node.tag,tuple(sorted(node.attrib.items())),(node.text or '').strip(),tuple(canonical(child) for child in node if child.tag!='default_activation'))
def promote(output='data/models/engineering_weight_transfer_v1'):
    out=(ROOT/output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('Fresh owned bundle required')
    baseline=ROOT/'data/models/engineering_stance_v1';source=ROOT/'data/research/locomotion_control/patient_weight_transfer80_scheduled'
    targetdir=ROOT/'data/derived/mechanics/patient_left80_stance98/lqr_target';endpoint=ROOT/'data/research/locomotion_control/linearization_bkcuq9xn/discrete_margin/linearization.npz'
    receipt_raw=(source/'report.json').read_bytes();receipt=json.loads(receipt_raw);summary=json.loads((source/'balance_summary.json').read_bytes())
    if not receipt['completed_horizon'] or receipt['seconds_requested']<10 or receipt['error'] is not None or receipt['fraction']!=1 or not receipt['gain_scheduled']:raise ValueError('Successful10s full scheduled transfer required')
    if not .78<=summary['final_left_load_fraction']<=.82 or summary['final_com_speed_m_s']>.01 or summary['max_clipped_count']!=0:raise ValueError('Observed transfer outside declared acceptance')
    raws={name:(source/name).read_bytes() for name in receipt['source_sha256']}
    for name,digest in receipt['source_sha256'].items():
        if sha(raws[name])!=digest:raise ValueError('Recorded acceptance source changed: '+name)
    if raws['policy.npz']!=(baseline/'linearization.npz').read_bytes() or raws['target.npz']!=(targetdir/'target.npz').read_bytes() or raws['target_policy.npz']!=endpoint.read_bytes():raise ValueError('Program differs from native accepted artifacts')
    regraw=(baseline/'registration.json').read_bytes();reg=json.loads(regraw);base_model=(ROOT/reg['model_path']).read_bytes()
    targetmeta_raw=(targetdir/'manifest.json').read_bytes();targetmeta=json.loads(targetmeta_raw)
    targetreg_raw=(ROOT/targetmeta['target_registration']).read_bytes();targetreg=json.loads(targetreg_raw);target_model=(ROOT/targetreg['model_path']).read_bytes()
    if sha(base_model)!=reg['model_sha256'] or sha(target_model)!=targetreg['model_sha256'] or sha(targetreg_raw)!=targetmeta['target_registration_sha256']:raise ValueError('Endpoint model provenance differs')
    a,b=ET.fromstring(base_model),ET.fromstring(target_model)
    if canonical(a)!=canonical(b):raise ValueError('Endpoint changes physical model beyond activation defaults')
    activations=lambda node:{item.attrib['name']:item.findtext('default_activation') for item in node.findall('.//ForceSet/objects/*') if 'Muscle' in item.tag}
    before,after=activations(a),activations(b)
    changed={key:{'baseline':before[key],'target':after[key]} for key in before if before[key]!=after[key]}
    with np.load(io.BytesIO(raws['policy.npz']),allow_pickle=False) as d:policy={k:d[k].copy() for k in d.files}
    with np.load(io.BytesIO(raws['target.npz']),allow_pickle=False) as d:target={k:d[k].copy() for k in d.files}
    with np.load(io.BytesIO(raws['target_policy.npz']),allow_pickle=False) as d:gain={k:d[k].copy() for k in d.files}
    for key in ('state_names','muscle_names'):
        if not np.array_equal(policy[key],target[key]) or not np.array_equal(policy[key],gain[key]):raise ValueError('State or muscle order differs')
    if any(float(d['target_mass_kg'])!=reg['target_mass_kg'] for d in (policy,target,gain)):raise ValueError('Patient masses differ')
    if str(target['source_model_sha256'].item())!=reg['model_sha256'] or str(gain['model_sha256'].item())!=targetreg['model_sha256']:raise ValueError('Program native model identity mismatch')
    if any(not np.array_equal(x,y) for x,y in [(target['x0'],target['x_target']),(target['u0'],target['u_target']),(target['x_target'],gain['x0']),(target['u_target'],gain['u0'])]):raise ValueError('Gain and target equilibrium mismatch')
    out.mkdir(parents=True);files={};origins={}
    def write(name,raw,origin=None,encoding=None):
        path=out/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw);relative=str(path.relative_to(ROOT));files[relative]=sha(raw)
        if origin is not None:origins[str(origin)]={'bundle_path':relative,'sha256':sha(raw),**({'encoding':encoding} if encoding else {})}
        return relative
    write('target.npz',raws['target.npz'],str((targetdir/'target.npz').relative_to(ROOT)))
    write('target_policy.npz',raws['target_policy.npz'],str(endpoint.relative_to(ROOT)))
    write('provenance/baseline_policy.npz',raws['policy.npz'],'data/models/engineering_stance_v1/linearization.npz')
    write('provenance/baseline_registration.json',regraw,'data/models/engineering_stance_v1/registration.json')
    write('provenance/baseline_manifest.json',(baseline/'manifest.json').read_bytes(),'data/models/engineering_stance_v1/manifest.json')
    write('provenance/target_model.osim',target_model,targetreg['model_path'])
    write('provenance/target_registration.json',targetreg_raw,targetmeta['target_registration'])
    write('provenance/target_manifest.json',targetmeta_raw,str((targetdir/'manifest.json').relative_to(ROOT)))
    write('provenance/target_snapshot.json.gz',gzip.compress((ROOT/targetmeta['snapshot_path']).read_bytes(),mtime=0),targetmeta['snapshot_path'],'gzip')
    write('acceptance/report.json',receipt_raw,str((source/'report.json').relative_to(ROOT)))
    write('acceptance/balance_summary.json',(source/'balance_summary.json').read_bytes(),str((source/'balance_summary.json').relative_to(ROOT)))
    trajectory=(source/'trajectory.json').read_bytes();write('acceptance/trajectory.json.gz',gzip.compress(trajectory,mtime=0),str((source/'trajectory.json').relative_to(ROOT)),'gzip')
    write('provenance/evaluator.py',raws['evaluator.py'],str((source/'evaluator.py').relative_to(ROOT)))
    write('provenance/promote_engineering_weight_transfer.py',Path(__file__).read_bytes(),'scripts/promote_engineering_weight_transfer.py')
    for name in ('execution.json',):write('provenance/native_'+name,(source/'plant'/name).read_bytes(),str((source/'plant'/name).relative_to(ROOT)))
    proof={'schema':'ihm.weight-transfer-physical-model-equivalence.v1','baseline_model_sha256':sha(base_model),'target_model_sha256':sha(target_model),
        'canonical_tree_sha256_without_default_activation':sha(json.dumps(canonical(a),separators=(',',':')).encode()),
        'same_except_default_activation':True,'changed_default_activations':changed,'comparison':'Exact parsedXML tags/attributes/text/childorder, ignoring only default_activation elements; not blanket parameter equivalence assumption',
        'runtime_model_substitution':False}
    write('provenance/model_equivalence.json',(json.dumps(proof,indent=2)+'\n').encode())
    program={'schema':'ihm.engineering-weight-transfer-program.v1','sampling_interval_s':.01,'start_s':1.,'transition_s':receipt['transition_s'],'fraction':1.,
        'blend':'10*p^3-15*p^4+6*p^5 withp=clip((t-start)/transition,0,1)',
        'reference':'Use baseline policy x0/u0, interpolate toward target.npz x_target/u_target; add derivative-of-blend times position delta to joint speed references',
        'target_alias_warning':'target.npz x0/u0 are aliases for ENDPOINT target state/excitation; never use them as baseline x0/u0',
        'gain':'Linear blend of baseline and endpoint K; muscle excitation clipped[0,1]',
        'target_path':str((out/'target.npz').relative_to(ROOT)),'target_policy_path':str((out/'target_policy.npz').relative_to(ROOT)),
        'baseline_binding':{'bundle_path':'data/models/engineering_stance_v1','manifest_sha256':sha((baseline/'manifest.json').read_bytes()),'registration_sha256':sha(regraw),
            'policy_sha256':sha(raws['policy.npz']),'native_model_sha256':sha(base_model),'target_mass_kg':reg['target_mass_kg']},
        'runtime_model_substitution':False,'external_root_forces':False,'prescribed_motion':False,'walking_demonstrated':False,'brain_trained':False,
        'scope':'Single demonstrated left-weight-transfer trajectory; bothfeet supported, no footstep or gait. Endpoint gain interpolation has no general nonlinear guarantee.'}
    write('program.json',(json.dumps(program,indent=2)+'\n').encode())
    manifest={'schema':'ihm.engineering-weight-transfer-bundle.v1','files':files,'origins':origins,'baseline_binding':program['baseline_binding'],
        'accepted_program':program,'acceptance':summary,'trajectory_uncompressed_sha256':sha(trajectory),
        'archival_path_scope':'Historical paths inside byte-retained provenance identify origins only. Runtime program uses bundlefiles plus digest-bound baseline stance bundle. Target model is proofonly; it mustneverreplace baseline runtime model.'}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for path,digest in files.items():
        if sha((ROOT/path).read_bytes())!=digest:raise AssertionError('Bundle filehash mismatch')
    print(json.dumps({'bundle':str(out),'files':len(files),'changed_activation_defaults':len(changed),'accepted_left_load_fraction':summary['final_left_load_fraction']}))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',default='data/models/engineering_weight_transfer_v1');a=p.parse_args();promote(a.output)
