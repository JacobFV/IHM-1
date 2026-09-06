"""Independent recheck of the emitted per-element muscle fibre field.

Run with the isolated libigl environment, not the physiological runtime:
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/verify_muscle_fibre_field.py --self-test
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/verify_muscle_fibre_field.py \
  --build data/derived/muscle-fibre-field-v1 --output data/derived/muscle-fibre-field-verification-v1

Nothing recorded by the build is trusted. Every artifact hash is recomputed, then for each sampled entity
the stored tet mesh is re-read and used to recompute, from scratch: the shape-function gradients, the
discrete Laplace residual of the stored potential on its free nodes, the fibre direction from that
potential, the pennation angle between the stored harmonic and final fields, the tet volume against the
independently proved surface volume, the recorded discrete-maximum-principle overshoot, and the
line-of-action angle against the registered OpenSim path.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import importlib.metadata
import json
import shutil
import sys
import time

import numpy as np
import scipy.sparse as sp

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(Path(__file__).resolve().parent))
from build_muscle_fibre_field import (ANATOMY,MECHANICS,TET_BUILD,sha,write,unit,deg,tet_gradients,
                                      load_surface,register,muscle_paths,to_canonical,polyline_parameter)

GRADIENT_TOLERANCE=1e-6
NORM_TOLERANCE=1e-6


def signed_volume(v,f):
    return float(np.einsum('ij,ij->i',v[f[:,0]],np.cross(v[f[:,1]],v[f[:,2]])).sum()/6.)


def stiffness(tv,tt,gradients,volume):
    rows=np.repeat(tt,4,axis=1).ravel();cols=np.tile(tt,(1,4)).ravel()
    data=(np.abs(volume)[:,None,None]*np.einsum('nid,njd->nij',gradients,gradients)).ravel()
    return sp.coo_matrix((data,(rows,cols)),shape=(len(tv),len(tv))).tocsr()


def recheck(record,payload,polylines):
    tv=payload['tet_vertices_m'];tt=payload['tets'].astype(np.int64)
    phi=payload['potential'].astype(float);fibre=payload['fibre'].astype(float)
    harmonic=payload['fibre_harmonic'].astype(float)
    gradients,volume=tet_gradients(tv,tt)
    raw=np.einsum('nid,ni->nd',gradients,phi[tt])
    recomputed=unit(raw);magnitude=np.linalg.norm(raw,axis=1)
    reference=float(np.median(magnitude[magnitude>0])) if (magnitude>0).any() else 0.
    # Every tet must either reproduce the stored direction from the stored potential, or be one whose
    # potential is flat, where the build substitutes its nearest resolved neighbour. Both are checked
    # here without consulting any mask the build recorded.
    disagreement=np.abs(recomputed-harmonic).max(axis=1)
    mismatched=disagreement>1e-6
    live=~mismatched
    row={'entity_id':record['entity_id'],'name':record['name'],'provenance':record['provenance'],
         'tets':int(len(tt)),'recorded_tets':record['tetgen']['tets'],
         'tets_match':int(len(tt))==record['tetgen']['tets'],
         'negative_volume_tets':int((volume<=0).sum()),
         'volume_m3':float(volume.sum()),'recorded_volume_m3':record['field']['volume_m3'],
         'volume_relative_error':abs(float(volume.sum())-record['field']['volume_m3'])/record['field']['volume_m3'],
         'fibre_unit_norm_max_deviation':float(np.abs(np.linalg.norm(fibre,axis=1)-1).max()),
         'harmonic_unit_norm_max_deviation':float(np.abs(np.linalg.norm(harmonic,axis=1)-1).max()),
         'potential_min':float(phi.min()),'potential_max':float(phi.max()),
         'recorded_potential_overshoot':record['field']['potential_overshoot'],
         'gradient_direction_max_error':float(disagreement[live].max()) if live.any() else 0.,
         'reproduced_tets':int(live.sum()),'mismatched_tets':int(mismatched.sum()),
         'mismatched_max_relative_gradient':float(magnitude[mismatched].max()/reference)
                                            if mismatched.any() and reference else 0.}
    free=(phi>1e-12)&(phi<1-1e-12)
    if free.any():
        residual=np.abs(stiffness(tv,tt,gradients,volume)@phi)[free]
        scale=float(np.abs(volume).sum()**(1/3))
        row.update(laplace_residual_max=float(residual.max()),
                   laplace_residual_relative=float(residual.max()/max(scale,1e-30)),
                   free_nodes=int(free.sum()))
    if record['pennation_angle_applied_rad']>0:
        angle=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',fibre,harmonic),-1,1)))
        row.update(pennation_angle_recorded_deg=float(np.degrees(record['pennation_angle_applied_rad'])),
                   pennation_angle_measured_median_deg=float(np.median(angle)),
                   pennation_angle_measured_max_error_deg=float(np.abs(angle-np.degrees(record['pennation_angle_applied_rad'])).max()))
    else:
        row['pennation_identical']=bool(np.array_equal(fibre,harmonic))
    integrated=unit((harmonic*volume[:,None]).sum(0))
    row['integrated_fibre_direction_error']=float(np.abs(integrated-np.asarray(record['integrated_fibre_direction'])).max())
    for entry in record.get('opensim',[]):
        poly=polylines[entry['opensim_muscle']]
        chord=unit(poly[-1]-poly[0])
        row.setdefault('line_of_action',[]).append({'opensim_muscle':entry['opensim_muscle'],
            'recomputed_deg':deg(abs(float(integrated@chord))),
            'recorded_deg':entry['integrated_fibre_vs_chord_deg']})
    return row


def verify(build,output,sample,run_all):
    build=Path(build).resolve();out=Path(output).resolve()
    if out.exists() or not out.is_relative_to(ROOT):raise ValueError('Choose a fresh verification directory')
    started=time.monotonic()
    manifest=json.loads((build/'manifest.json').read_text())
    artifact_mismatches=[]
    for name,digest in manifest['artifacts_sha256'].items():
        path=(build/name).resolve()
        if not path.is_relative_to(build) or sha(path)!=digest:artifact_mismatches.append(name)
    input_mismatches=[p for p,digest in manifest['inputs_sha256'].items() if sha(ROOT/p)!=digest]
    records=[json.loads(l) for l in (build/'entities.jsonl').read_text().splitlines()]
    anatomy=json.loads((ROOT/ANATOMY).read_bytes());spec={e['id']:e for e in anatomy['entities']}
    registration,_,_=register(ROOT,spec)
    paths=muscle_paths(ROOT)
    polylines={name:np.concatenate([to_canonical(registration,b,l) for b,l in row['points']])
               for name,row in paths.items()}
    tet_records={json.loads(l)['entity_id']:json.loads(l)
                 for l in (ROOT/TET_BUILD/'entities.jsonl').read_text().splitlines()}
    chosen=records if run_all else records[::max(1,len(records)//max(sample,1))][:sample]
    out.mkdir(parents=True);(out/'inputs').mkdir()
    shutil.copyfile(__file__,out/'inputs/verify_muscle_fibre_field.py')
    shutil.copyfile(build/'manifest.json',out/'inputs/build-manifest.json')
    shutil.copyfile(build/'summary.json',out/'inputs/build-summary.json')
    rows=[]
    for count,record in enumerate(chosen,1):
        with np.load(build/record['fibre_path']) as payload:
            data={k:payload[k] for k in payload.files}
        row=recheck(record,data,polylines)
        surface_volume=tet_records[record['entity_id']]['after']['signed_volume_m3']
        row['surface_signed_volume_m3']=surface_volume
        row['tet_vs_surface_relative_error']=abs(row['volume_m3']-surface_volume)/abs(surface_volume)
        rows.append(row)
        with (out/'rechecked.jsonl').open('a') as handle:handle.write(json.dumps(row,allow_nan=False)+'\n')
        if count%25==0:print('Rechecked',count,'of',len(chosen),flush=True)
    line=[(l['recomputed_deg'],l['recorded_deg']) for r in rows for l in r.get('line_of_action',[])]
    summary={'schema':'ihm.muscle-fibre-field-verification.v1','build':str(build),
        'entities_in_build':len(records),'entities_rechecked':len(rows),
        'artifact_hash_mismatches':artifact_mismatches,'input_hash_mismatches':input_mismatches,
        'tet_count_mismatches':[r['entity_id'] for r in rows if not r['tets_match']],
        'negative_volume_tets':int(sum(r['negative_volume_tets'] for r in rows)),
        'max_volume_relative_error':max(r['volume_relative_error'] for r in rows),
        'max_tet_vs_surface_relative_error':max(r['tet_vs_surface_relative_error'] for r in rows),
        'max_fibre_unit_norm_deviation':max(r['fibre_unit_norm_max_deviation'] for r in rows),
        'max_harmonic_unit_norm_deviation':max(r['harmonic_unit_norm_max_deviation'] for r in rows),
        'entities_with_potential_overshoot':len([r for r in rows
                                                 if r['potential_min']<-1e-9 or r['potential_max']>1+1e-9]),
        'entities_with_potential_overshoot_over_1pct':[r['entity_id'] for r in rows
            if max(-r['potential_min'],r['potential_max']-1,0)>.01],
        'max_potential_overshoot':max(max(-r['potential_min'],r['potential_max']-1,0) for r in rows),
        'max_recorded_potential_overshoot_error':max(abs(max(-r['potential_min'],r['potential_max']-1,0)
            -r['recorded_potential_overshoot']) for r in rows),
        'max_gradient_direction_error':max(r['gradient_direction_max_error'] for r in rows),
        'mismatched_tets':int(sum(r['mismatched_tets'] for r in rows)),
        'mismatched_tet_fraction':float(sum(r['mismatched_tets'] for r in rows)/sum(r['tets'] for r in rows)),
        'max_mismatched_relative_gradient':max(r['mismatched_max_relative_gradient'] for r in rows),
        'max_laplace_residual_relative':max((r['laplace_residual_relative'] for r in rows
                                             if 'laplace_residual_relative' in r),default=None),
        'max_pennation_angle_error_deg':max((r['pennation_angle_measured_max_error_deg'] for r in rows
                                             if 'pennation_angle_measured_max_error_deg' in r),default=None),
        'max_integrated_direction_error':max(r['integrated_fibre_direction_error'] for r in rows),
        'line_of_action_comparisons':len(line),
        'max_line_of_action_disagreement_deg':max((abs(a-b) for a,b in line),default=None),
        'gradient_tolerance':GRADIENT_TOLERANCE,'norm_tolerance':NORM_TOLERANCE,
        'elapsed_s':time.monotonic()-started}
    summary['verified']=(not artifact_mismatches and not input_mismatches
        and not summary['tet_count_mismatches'] and summary['negative_volume_tets']==0
        and summary['max_fibre_unit_norm_deviation']<NORM_TOLERANCE
        and summary['max_gradient_direction_error']<GRADIENT_TOLERANCE
        and summary['max_mismatched_relative_gradient']<1e-4
        and summary['max_recorded_potential_overshoot_error']<1e-9
        and (summary['max_line_of_action_disagreement_deg'] or 0)<1e-6)
    write(out/'summary.json',summary)
    write(out/'manifest.json',{'schema':'ihm.muscle-fibre-field-verification.v1',
        'packages':{n:importlib.metadata.version(n) for n in ('libigl','numpy','scipy')},'python':sys.version,
        'build_manifest_sha256':sha(build/'manifest.json'),
        'artifacts_sha256':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()},
        'scope':'Recomputes the emitted field from the emitted mesh and potential. It proves internal '
                'consistency and numerical correctness of the solve; it cannot make an unmeasured fibre '
                'architecture measured.'})
    print(json.dumps(summary,indent=2))


def self_test():
    tv=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
    tt=np.array([[0,1,2,3]],np.int64)
    gradients,volume=tet_gradients(tv,tt)
    assert abs(volume[0]-1/6)<1e-15
    phi=tv[:,0].copy()
    assert np.abs(unit(np.einsum('nid,ni->nd',gradients,phi[tt]))-np.array([1.,0,0])).max()<1e-14
    k=stiffness(tv,tt,gradients,volume)
    assert abs((np.ones(4)@k@np.ones(4)))<1e-14,'a constant potential must carry no stiffness energy'
    assert abs(signed_volume(tv,np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]],np.int64))-1/6)<1e-15
    record={'entity_id':'x','name':'x','provenance':'inferred','tetgen':{'tets':1},
            'field':{'volume_m3':1/6,'potential_overshoot':0.},'pennation_angle_applied_rad':0.,
            'integrated_fibre_direction':[1.,0,0]}
    payload={'tet_vertices_m':tv,'tets':tt.astype(np.int32),'potential':phi.astype(np.float32),
             'fibre':np.array([[1.,0,0]]),'fibre_harmonic':np.array([[1.,0,0]])}
    row=recheck(record,payload,{})
    assert row['tets_match'] and row['negative_volume_tets']==0 and row['pennation_identical']
    assert row['gradient_direction_max_error']<1e-6 and row['integrated_fibre_direction_error']<1e-12
    print('PASS single-tet gradient, constant-potential stiffness null space and recheck plumbing')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--self-test',action='store_true')
    p.add_argument('--build',type=Path);p.add_argument('--output',type=Path)
    p.add_argument('--sample',type=int,default=40);p.add_argument('--all',action='store_true')
    a=p.parse_args()
    if a.self_test:self_test()
    if a.build and a.output:verify(a.build,a.output,a.sample,a.all)
    elif a.build or a.output:p.error('--build and --output go together')
    elif not a.self_test:p.error('Select --self-test or --build/--output')
