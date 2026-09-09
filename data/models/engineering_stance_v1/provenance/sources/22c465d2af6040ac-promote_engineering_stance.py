#!/usr/bin/env python3
"""Freeze a patient-matched stance bundle after corrected native paired acceptance.

Never overwrites a bundle. Retained archival receipts preserve original bytes and
origin paths; only promoted registration runtime paths are rewritten.
"""
import argparse,gzip,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def sha(raw):return hashlib.sha256(raw).hexdigest()
def promote(artifact,registration,feedback,held,output='data/models/engineering_stance_v1'):
    artifact,registration,feedback,held=map(lambda p:(ROOT/p).resolve() if not Path(p).is_absolute() else Path(p).resolve(),(artifact,registration,feedback,held))
    out=(ROOT/output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('Fresh owned bundle path required')
    regraw=registration.read_bytes();reg=json.loads(regraw);npz=artifact.read_bytes()
    with np.load(artifact,allow_pickle=False) as d:
        if str(d['model_sha256'].item())!=reg['model_sha256'] or abs(float(d['dt_s'])-.01)>1e-12 or abs(float(d['target_mass_kg'])-reg['target_mass_kg'])>1e-9:raise ValueError('Policy and registered patient identity mismatch')
    reports={arm:json.loads((directory/'report.json').read_bytes()) for arm,directory in [('feedback',feedback),('held',held)]}
    full=reports['feedback']
    if not full['completed_horizon'] or full['seconds_requested']<10 or full['error'] is not None:raise ValueError('Completed error-free10s feedback acceptance required')
    if full['final_com_speed_m_s']>.05 or np.linalg.norm(full['final_com_displacement_m'])>.05:raise ValueError('Feedback has not recovered within declared COM acceptance bounds')
    if full['pelvis_height_m']<.6 or abs(full['pelvis_tilt_rad'])>1 or abs(full['pelvis_list_rad'])>1:raise ValueError('Feedback acceptance contains a fall')
    for arm,report in reports.items():
        if report['arm']!=arm or report['target_mass_kg']!=reg['target_mass_kg']:raise ValueError('Paired acceptance arm/mass mismatch')
        if report['perturbation'].get('point_basis')!='current pelvis COM in source ground; recomputed each step':raise ValueError('Corrected COM force application required')
        if sha(npz) not in report['source_sha256'].values():raise ValueError('Acceptance did not use promoted artifact')
    if reports['feedback']['perturbation']!=reports['held']['perturbation']:raise ValueError('Perturbations are not matched')
    out.mkdir(parents=True);evidence=out/'provenance';evidence.mkdir();origins={};files={}
    def write(relative,raw,origin=None):
        path=out/relative;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
        rel=str(path.relative_to(ROOT));files[rel]=sha(raw)
        if origin is not None:origins[str(origin)]={'bundle_path':rel,'sha256':sha(raw)}
        return rel
    def retain(path,digest=None,fallback=None):
        original=str(path);resolved=Path(path) if Path(path).is_absolute() else ROOT/path
        raw=resolved.read_bytes() if resolved.is_file() else None
        if raw is None or digest is not None and sha(raw)!=digest:
            if fallback is not None and fallback.is_file():raw=fallback.read_bytes()
            if raw is None or digest is not None and sha(raw)!=digest:raise ValueError('Evidence source missing/changed: '+original)
        prior=origins.get(original)
        if prior and prior['sha256']==sha(raw):return prior['bundle_path']
        return write('provenance/sources/'+sha(raw)[:16]+'-'+resolved.name,raw,original)
    promoted=dict(reg)
    # Preserve every authoritative native source check, changing only its path.
    promoted['sources']={}
    canonical={}
    for origin,digest in reg['sources'].items():
        mapped=retain(origin,digest);promoted['sources'][mapped]=digest
        if origin.startswith('data/derived/canonical/'):
            promoted['sources'][origin]=digest;canonical[origin]=digest
    for field,name in [('model_path','model.osim'),('catalog_path','catalog.json'),('insert_path','insert.xml')]:
        raw=(ROOT/reg[field]).read_bytes();digestfield=field.replace('_path','_sha256')
        if sha(raw)!=reg[digestfield]:raise ValueError('Registered artifact source changed')
        promoted[field]=write(name,raw,reg[field])
    # Additional registration references are runtime-resolvable bundle evidence.
    for field in ('base_model_path','rematerialization_receipt','static_resting_solve','native_initial_acceptance','canonical_source_audit'):
        if field in reg:promoted[field]=retain(reg[field])
    write('linearization.npz',npz,str(artifact.relative_to(ROOT)))
    write('initial_pose.json',(json.dumps(reg['initial_pose'],indent=2)+'\n').encode())
    write('equilibrium_excitations.json',(json.dumps(reg['equilibrium_excitations'],indent=2)+'\n').encode())
    write('provenance/original_registration.json',regraw,str(registration.relative_to(ROOT)))
    acceptance={}
    for arm,directory in [('feedback',feedback),('held',held)]:
        report=reports[arm];retained_sources={}
        for origin,digest in report['source_sha256'].items():
            mapped=retain(origin,digest,directory/Path(origin).name);retained_sources[mapped]=digest
        reportpath=write('acceptance/'+arm+'/report.json',(directory/'report.json').read_bytes(),str((directory/'report.json').relative_to(ROOT)))
        trajectory=(directory/'trajectory.json').read_bytes()
        compressed=write('acceptance/'+arm+'/trajectory.json.gz',gzip.compress(trajectory,mtime=0),str((directory/'trajectory.json').relative_to(ROOT)))
        acceptance[arm]={'report':reportpath,'trajectory_gzip':compressed,'uncompressed_trajectory_sha256':sha(trajectory),'retained_source_sha256':retained_sources}
    for extra in (artifact.parent/'report.json',feedback/'paired_acceptance.json',Path(__file__),ROOT/'ihm/native/stance_lqr_controller.py',ROOT/'ihm/native/stance_lqr.py'):
        if extra.is_file():retain(str(extra.relative_to(ROOT)))
    promoted['default_enabled']=False
    promoted['native_acceptance_complete']=True
    promoted['promoted_bundle_scope']='10s patient-matched engineering stance COM-push acceptance; no walking or cortical motor claim'
    promoted['canonical_source_bindings']=canonical
    promoted['sources'].update({path:digest for path,digest in files.items() if path!=str((out/'registration.json').relative_to(ROOT))})
    write('registration.json',(json.dumps(promoted,indent=2)+'\n').encode())
    manifest={'schema':'ihm.engineering-stance-bundle.v1','target_mass_kg':reg['target_mass_kg'],'native_model_sha256':reg['model_sha256'],
        'sampling_interval_s':.01,'files':files,'origins':origins,'acceptance':acceptance,'canonical_source_bindings':canonical,
        'acceptance_limits':{'minimum_feedback_horizon_s':10.,'maximum_final_com_speed_m_s':.05,'maximum_final_com_displacement_m':.05},
        'controller':'engineering discrete native LQR; ill-conditioned DARE residual about0.415, not an optimality certificate','walking_demonstrated':False,'brain_trained':False,
        'archival_paths':'Paths inside byte-preserved original evidence receipts identify historical origins. Use origins mapping for retained copies; promoted registration runtime paths resolve within bundle except explicit live canonical bindings.'}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    # Verify every promoted source check with the same filesystem rule as native.
    for path,digest in promoted['sources'].items():
        if sha((ROOT/path).read_bytes())!=digest:raise AssertionError('Promoted source verification failed')
    print(json.dumps({'bundle':str(out),'files':len(files),'origins':len(origins),'canonical_bindings':canonical}));return manifest
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--artifact',required=True);p.add_argument('--registration',required=True);p.add_argument('--feedback',required=True);p.add_argument('--held',required=True);p.add_argument('--output',default='data/models/engineering_stance_v1');a=p.parse_args();promote(a.artifact,a.registration,a.feedback,a.held,a.output)
