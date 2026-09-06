"""Independent proof that repaired muscular surfaces tetrahedralise.

Run with the isolated libigl environment, not the physiological runtime:
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/verify_muscle_tet_ready_surfaces.py --self-test
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/verify_muscle_tet_ready_surfaces.py \
  --build data/derived/muscle-tet-ready-v1 --output data/derived/muscle-tet-ready-verification-v1 --sample 24

Every recorded diagnostic is recomputed here from the stored geometry rather than trusted, then a
stratified sample is handed to TetGen twice: once as the losslessly welded canonical surface (the
control) and once as the repaired surface. TetGen writes to stdout from native code, so the file
descriptor is captured and the message stored verbatim.
"""
from pathlib import Path
import argparse
import ctypes
import gzip
import hashlib
import importlib.metadata
import json
import os
import shutil
import sys
import tempfile
import time

import numpy as np
import igl
import igl.copyleft.cgal as cgal
from igl.copyleft import tetgen

LIBC=ctypes.CDLL(None)
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(Path(__file__).resolve().parent))
from build_muscle_tet_ready_surfaces import (sha,write,diagnose,signed_volume,weld_exact,
                                             drop_repeated_index,drop_zero_area,dedupe_faces,orient_outward)


def load_surface(path):
    payload=json.loads(gzip.decompress(Path(path).read_bytes()))
    v=np.ascontiguousarray(np.asarray(payload['positions'],float).reshape(-1,3))
    f=np.ascontiguousarray(np.asarray(payload['indices'],np.int64).reshape(-1,3))
    return payload,v,f


def welded_control(v,f):
    v,f=weld_exact(v,f);f,_=drop_repeated_index(f);f,_=drop_zero_area(v,f)
    return np.ascontiguousarray(v),np.ascontiguousarray(f)


def tetrahedralize(v,f,flags):
    """Run TetGen in a forked child with its native stdout and stderr captured verbatim.

    TetGen aborts the whole process on some inputs instead of raising, so the call is isolated;
    a child killed by a signal is recorded as a failure rather than losing the run."""
    handle,name=tempfile.mkstemp(suffix='.tetgen');os.close(handle)
    reader,writer=os.pipe();started=time.monotonic();outcome={'flags':flags}
    pid=os.fork()
    if pid==0:
        try:
            os.close(reader);sink=os.open(name,os.O_WRONLY)
            sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None);os.dup2(sink,1);os.dup2(sink,2)
            child={}
            try:
                result=tetgen.tetrahedralize(v,f,flags=flags)
                child.update(status=int(result[-1]),tet_vertices=int(len(result[0])),tets=int(len(result[1])))
                if child['status']==0 and len(result[1]):
                    tv=np.ascontiguousarray(np.asarray(result[0],float))
                    tt=np.ascontiguousarray(np.asarray(result[1],np.int64))
                    child['tet_volume_m3']=float(np.abs(np.asarray(igl.volume(tv,tt),float)).sum())
            except BaseException as error:
                child.update(status=None,exception=type(error).__name__,message=str(error))
            sys.stdout.flush();sys.stderr.flush();LIBC.fflush(None)
            os.write(writer,json.dumps(child,allow_nan=False).encode());os.close(writer)
        finally:
            os._exit(0)
    os.close(writer);chunks=[]
    while True:
        chunk=os.read(reader,65536)
        if not chunk:break
        chunks.append(chunk)
    os.close(reader);_,wait_status=os.waitpid(pid,0)
    if chunks:outcome.update(json.loads(b''.join(chunks).decode()))
    if os.WIFSIGNALED(wait_status):
        outcome.update(status=None,exception='ProcessAborted',
                       message='TetGen terminated the process with signal %d'%os.WTERMSIG(wait_status))
    elif not chunks:
        outcome.update(status=None,exception='NoResult',message='TetGen child exited without a result')
    outcome['tetgen_stdout']=Path(name).read_text(errors='replace').strip().splitlines()
    Path(name).unlink()
    # A run can return status 0 having meshed almost nothing: the right inferior lung
    # lobe yields 9 tets over 4.77e-12 m3 against a 6.66e-4 m3 surface. Gate on the
    # volume the tets actually account for. The ratio is computed here, against this
    # call's own input, so the gate cannot be defeated by a caller that never sets it.
    outcome['volume_gate']=VOLUME_GATE
    surface=abs(signed_volume(v,f))
    meshed=outcome.get('tet_volume_m3')
    error=(meshed-surface)/surface if meshed is not None and surface>0 else None
    outcome['tet_volume_vs_surface_relative_error']=error
    outcome['succeeded']=(outcome.get('status')==0 and outcome.get('tets',0)>0
                          and error is not None and abs(error)<=VOLUME_GATE)
    outcome['seconds']=time.monotonic()-started
    return outcome


VOLUME_GATE=.05


def strata(records,count):
    """Face-count terciles, plus every entity whose repair needed the self-union step and every
    entity the prior readiness screen blocked, up to the requested budget."""
    ordered=sorted(records,key=lambda r:(r['after']['faces'],r['entity_id']))
    n=len(ordered);picked=[]
    for lo,hi in ((0,n//3),(n//3,2*n//3),(2*n//3,n)):
        band=ordered[lo:hi]
        if not band:continue
        take=max(1,count//3)
        step=max(1,len(band)//take)
        picked+= band[::step][:take]
    prioritised=[r for r in ordered if r['prior_blocking_reasons'] or
                 any(s['stage']=='extract_outer_manifold' for s in r['operations'])]
    for r in prioritised:
        if len(picked)>=count+len(prioritised):break
        picked.append(r)
    seen={};[seen.setdefault(r['entity_id'],r) for r in picked]
    return sorted(seen.values(),key=lambda r:r['after']['faces'])


def self_test():
    v=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]])
    f=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]],np.int64)
    good=tetrahedralize(v,f,'pYq1.414')
    assert good['succeeded'] and good['tets']>0 and good['status']==0
    assert abs(good['tet_volume_m3']-signed_volume(v,f))<1e-12
    assert not any('KABOOOM' in line for line in good['tetgen_stdout'])
    crossing_v=np.array([[-1.,-1,0],[1.,-1,0],[0,1.,0],[0,-.5,-1],[0,-.5,1],[0,.5,0]])
    bad=tetrahedralize(crossing_v,np.array([[0,1,2],[3,4,5]],np.int64),'pYq1.414')
    assert not bad['succeeded'] and bad['tetgen_stdout'],'a failing run must retain its native message'
    open_shell=tetrahedralize(v,f[:3],'pYq1.414')
    assert not open_shell['succeeded']
    assert abs(good['tet_volume_vs_surface_relative_error'])<1e-12,'a sound mesh accounts for its volume'
    original_volume=igl.volume
    try:
        igl.volume=lambda tv,tt:np.asarray(original_volume(tv,tt),float)*1e-9
        starved=tetrahedralize(v,f,'pYq1.414')
    finally:igl.volume=original_volume
    assert starved['status']==0 and starved['tets']>0,'the starved fixture must still exit cleanly'
    assert not starved['succeeded'],'a zero-exit run that meshed no volume is not a success'
    original=tetgen.tetrahedralize
    try:
        tetgen.tetrahedralize=lambda *args,**kwargs:os.abort()
        aborted=tetrahedralize(v,f,'pYq1.414')
    finally:tetgen.tetrahedralize=original
    assert not aborted['succeeded'] and aborted['exception']=='ProcessAborted','a native abort must not end the run'
    assert tetrahedralize(v,f,'pYq1.414')['succeeded'],'the parent must survive an aborted child'
    seam=welded_control(v[f].reshape(-1,3),np.arange(12,dtype=np.int64).reshape(-1,3))
    assert len(seam[0])==4 and len(seam[1])==4
    records=[dict(entity_id='e%03d'%i,after=dict(faces=10*i),prior_blocking_reasons=[],operations=[]) for i in range(30)]
    chosen=strata(records,9)
    assert 9<=len(chosen)<=12 and len({r['entity_id'] for r in chosen})==len(chosen)
    assert chosen[0]['after']['faces']<chosen[-1]['after']['faces']
    print('PASS captured-stdout TetGen success, failure and stratification checks')


def verify(build,output,sample,flags,run_all):
    build=Path(build).resolve();out=Path(output).resolve()
    if out.exists():raise ValueError('Choose a fresh verification directory')
    manifest=json.loads((build/'manifest.json').read_text())
    for name,digest in manifest['artifacts_sha256'].items():
        path=(build/name).resolve()
        if not path.is_relative_to(build) or sha(path)!=digest:raise ValueError('Build artifact changed: '+name)
    for name,digest in manifest['source_geometry_sha256'].items():
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT) or sha(path)!=digest:raise ValueError('Canonical source changed: '+name)
    if sha(ROOT/'data/derived/canonical/anatomy.json')!=manifest['anatomy_sha256']:
        raise ValueError('Canonical anatomical identity changed')
    packages={n:importlib.metadata.version(n) for n in ('libigl','numpy')}
    if packages['libigl']!='2.6.2':raise ValueError('This verification pins libigl 2.6.2')
    records=[json.loads(l) for l in (build/'entities.jsonl').read_text().splitlines()]
    out.mkdir(parents=True);(out/'inputs').mkdir()
    shutil.copyfile(__file__,out/'inputs/verify_muscle_tet_ready_surfaces.py')
    shutil.copyfile(build/'manifest.json',out/'inputs/build-manifest.json')
    shutil.copyfile(build/'summary.json',out/'inputs/build-summary.json')
    started=time.monotonic();rediagnosed=[];mismatches=[]
    for record in records:
        _,v,f=load_surface(build/record['output_path'])
        current=diagnose(v,f)
        differing={k:(record['after'][k],current[k]) for k in current
                   if k=='signed_volume_m3' and abs(current[k]-record['after'][k])>1e-18
                   or k!='signed_volume_m3' and current[k]!=record['after'][k]}
        if differing:mismatches.append(dict(entity_id=record['entity_id'],differing=differing))
        rediagnosed.append(dict(entity_id=record['entity_id'],**current))
        with (out/'rediagnosed.jsonl').open('a') as handle:
            handle.write(json.dumps(rediagnosed[-1],allow_nan=False)+'\n')
    if mismatches:write(out/'diagnostic-mismatches.json',mismatches)
    chosen=records if run_all else strata(records,sample)
    by_id={r['entity_id']:r for r in records};trials=[]
    for record in chosen:
        payload,rv,rf=load_surface(build/record['output_path'])
        source=json.loads(gzip.decompress((ROOT/record['source_path']).read_bytes()))
        cv,cf=welded_control(np.asarray(source['positions'],float).reshape(-1,3),
                             np.asarray(source['indices'],np.int64).reshape(-1,3))
        control=tetrahedralize(cv,cf,flags);repaired=tetrahedralize(rv,rf,flags)
        surface_volume=signed_volume(rv,rf)
        trial=dict(entity_id=record['entity_id'],name=record['name'],
            faces_welded_control=int(len(cf)),faces_repaired=int(len(rf)),
            prior_blocking_reasons=record['prior_blocking_reasons'],
            prior_topological_candidate=record['prior_topological_candidate'],
            self_intersecting_pairs_before=record['before']['self_intersecting_face_pairs'],
            nonmanifold_vertices_after=record['after']['nonmanifold_vertices'],
            nonmanifold_edges_after=record['after']['nonmanifold_edges'],
            face_components_after=record['after']['face_components'],
            signed_volume_before_m3=record['before']['signed_volume_m3'],
            signed_volume_after_m3=surface_volume,
            signed_volume_relative_drift=record['signed_volume_relative_drift'],
            tet_volume_vs_surface_relative_error=(repaired['tet_volume_m3']-surface_volume)/surface_volume
                if repaired.get('tet_volume_m3') and surface_volume else None,
            welded_control=control,repaired=repaired)
        trials.append(trial)
        with (out/'tetgen-trials.jsonl').open('a') as handle:
            handle.write(json.dumps(trial,allow_nan=False)+'\n')
        print('%-22s control=%-5s repaired=%-5s tets=%s'%(record['entity_id'],control['succeeded'],
              repaired['succeeded'],repaired.get('tets')),flush=True)
    successes=[t for t in trials if t['repaired']['succeeded']]
    errors=[abs(t['tet_volume_vs_surface_relative_error']) for t in successes
            if t['tet_volume_vs_surface_relative_error'] is not None]
    drifts=[abs(r['signed_volume_relative_drift']) for r in records if r['signed_volume_relative_drift'] is not None]
    summary=dict(schema='ihm.muscle-tet-ready-verification.v1',build=str(build),
        entities_rediagnosed=len(rediagnosed),diagnostic_mismatches=len(mismatches),
        recorded_diagnostics_reproduced=not mismatches,
        tetgen_flags=flags,tetgen_trials=len(trials),
        repaired_tetrahedralized=len(successes),welded_control_tetrahedralized=sum(t['welded_control']['succeeded'] for t in trials),
        repaired_success_rate=len(successes)/len(trials) if trials else 0.0,
        trials_with_prior_blockers=sum(bool(t['prior_blocking_reasons']) for t in trials),
        prior_blocked_tetrahedralized=sum(t['repaired']['succeeded'] and bool(t['prior_blocking_reasons']) for t in trials),
        trials_with_nonmanifold_vertices=sum(t['nonmanifold_vertices_after']>0 for t in trials),
        nonmanifold_vertex_trials_tetrahedralized=sum(t['repaired']['succeeded'] and t['nonmanifold_vertices_after']>0 for t in trials),
        total_tets=sum(t['repaired'].get('tets',0) for t in successes),
        median_tets=float(np.median([t['repaired']['tets'] for t in successes])) if successes else 0.0,
        max_absolute_tet_vs_surface_volume_error=max(errors,default=0.0),
        max_absolute_repair_volume_drift=max(drifts,default=0.0),
        median_absolute_repair_volume_drift=float(np.median(drifts)) if drifts else 0.0,
        entities_with_volume_drift_over_1pct=[r['entity_id'] for r in records
            if r['signed_volume_relative_drift'] is not None and abs(r['signed_volume_relative_drift'])>.01],
        repaired_process_aborts=sum(t['repaired'].get('exception')=='ProcessAborted' for t in trials),
        welded_control_process_aborts=sum(t['welded_control'].get('exception')=='ProcessAborted' for t in trials),
        failures=[dict(entity_id=t['entity_id'],name=t['name'],faces=t['faces_repaired'],
                       nonmanifold_vertices=t['nonmanifold_vertices_after'],
                       nonmanifold_edges_after=t.get('nonmanifold_edges_after'),
                       exception=t['repaired'].get('exception'),message=t['repaired'].get('message'),
                       tetgen_stdout=t['repaired']['tetgen_stdout'][-12:])
                  for t in trials if not t['repaired']['succeeded']],
        elapsed_s=time.monotonic()-started)
    write(out/'summary.json',summary)
    artifacts={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
    write(out/'manifest.json',dict(schema='ihm.muscle-tet-ready-verification.v1',packages=packages,python=sys.version,
        build_manifest_sha256=sha(build/'manifest.json'),anatomy_sha256=manifest['anatomy_sha256'],
        native_cgal_sha256=sha(Path(cgal.pyigl_copyleft_cgal.__file__)),
        native_tetgen_sha256=sha(Path(tetgen.pyigl_copyleft_tetgen.__file__)),
        artifacts_sha256=artifacts,limitations=[
            'A successful tetrahedralisation proves the surface is a valid closed PLC for TetGen; it does not '
            'establish anatomical correctness, inter-entity disjointness or material ownership.',
            'The welded control is the losslessly welded canonical surface with repeated-index and exactly '
            'collinear faces dropped, so it isolates the contribution of the later repair steps only.',
            'Tet volume is compared against the repaired surface integral, not against any tissue mass ledger.']))
    print(json.dumps({k:v for k,v in summary.items() if k!='failures'},indent=2))
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--self-test',action='store_true')
    p.add_argument('--build',type=Path,default=ROOT/'data/derived/muscle-tet-ready-v1')
    p.add_argument('--output',type=Path);p.add_argument('--sample',type=int,default=24)
    p.add_argument('--all',action='store_true');p.add_argument('--flags',default='pYq1.414')
    a=p.parse_args()
    if a.self_test:self_test()
    if a.output:verify(a.build,a.output,a.sample,a.flags,a.all)
    if not a.self_test and not a.output:p.error('Select --self-test or --output')
