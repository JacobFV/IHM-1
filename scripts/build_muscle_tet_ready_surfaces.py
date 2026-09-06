"""Conservative repair of canonical muscular surfaces into tetrahedralisation-ready shells.

Run with the isolated libigl environment, not the physiological runtime:
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_muscle_tet_ready_surfaces.py --self-test
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_muscle_tet_ready_surfaces.py --output data/derived/muscle-tet-ready-v1

Repair order, each step recorded per entity in operations[]:
  1 weld_exact                  numerically identical coordinate triples merged; no point moves
  2 drop_repeated_index_faces   faces reusing a vertex; zero volume, edge parity preserved
  3 drop_zero_area_faces        exactly collinear faces; zero volume, may expose a degenerate slit
  4 dedupe_coincident_faces     same-winding copies collapse to one; opposite-winding pairs are
                                zero-volume fins and both go, which preserves edge parity
  5 orient_consistent_outward   per-patch bfs_orient then flip patches with negative signed volume
  6 resolve_self_intersections  CGAL remesh_self_intersections(stitch_all) exact arrangement, then 1-5 again
  7 extract_outer_manifold      CGAL self-union (winding-number outer shell) when 6 leaves nonmanifold edges
No step fills a hole or moves an existing coordinate. An output that still carries boundary edges is
flagged hole_filling_required and left open.
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
import igl
import igl.copyleft.cgal as cgal

ROOT=Path(__file__).resolve().parents[1]
EMPTY_V=np.zeros((0,3));EMPTY_F=np.zeros((0,3),np.int64)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def signed_volume(v,f):
    if not len(f):return 0.0
    return float(np.einsum('ij,ij->i',v[f[:,0]],np.cross(v[f[:,1]],v[f[:,2]])).sum()/6.0)


def edge_incidence(f):
    if not len(f):return 0,0
    e=np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1)
    _,count=np.unique(e,axis=0,return_counts=True)
    return int((count==1).sum()),int((count>2).sum())


def winding_sign(f):
    p=np.argsort(f,axis=1)
    return np.where((p[:,1]-p[:,0])*(p[:,2]-p[:,0])*(p[:,2]-p[:,1])>0,1,-1)


def weld_exact(v,f):
    w,inverse=np.unique(v,axis=0,return_inverse=True)
    return w,np.ascontiguousarray(inverse[f.ravel()].reshape(-1,3),dtype=np.int64)


def drop_repeated_index(f):
    keep=(f[:,0]!=f[:,1])&(f[:,1]!=f[:,2])&(f[:,0]!=f[:,2])
    return np.ascontiguousarray(f[keep]),int((~keep).sum())


def drop_zero_area(v,f):
    n=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]])
    keep=np.any(n!=0,axis=1)
    return np.ascontiguousarray(f[keep]),int((~keep).sum())


def dedupe_faces(f):
    """Collapse coincident triangles. Removing both members of an opposite-winding pair
    changes every shared edge's incidence by two, so closed surfaces stay closed."""
    if not len(f):return f,0,0
    key=np.sort(f,axis=1);order=np.lexsort(key.T[::-1]);ordered=key[order]
    bounds=np.r_[0,1+np.flatnonzero(np.any(ordered[1:]!=ordered[:-1],axis=1)),len(ordered)]
    sign=winding_sign(f);keep=np.ones(len(f),bool);same=0;fins=0
    for lo,hi in zip(bounds[:-1],bounds[1:]):
        if hi-lo==1:continue
        rows=order[lo:hi];s=sign[rows];pos=rows[s>0];neg=rows[s<0]
        paired=min(len(pos),len(neg));fins+=2*paired
        same+=len(rows)-2*paired-(1 if len(rows)>2*paired else 0)
        keep[rows]=False
        if len(pos)>paired:keep[pos[paired]]=True
        elif len(neg)>paired:keep[neg[paired]]=True
    return np.ascontiguousarray(f[keep]),same,fins


def orient_outward(v,f):
    if not len(f):return f,0,0,0
    oriented,patch=igl.bfs_orient(f)
    oriented=np.ascontiguousarray(np.asarray(oriented,np.int64));patch=np.asarray(patch).ravel()
    if not np.array_equal(np.sort(oriented,axis=1),np.sort(f,axis=1)):
        raise ValueError('bfs_orient must only reverse windings')
    reversed_faces=int((oriented!=f).any(axis=1).sum());flipped=0
    for index in range(int(patch.max())+1):
        mask=patch==index
        if signed_volume(v,oriented[mask])<0:
            oriented[mask]=oriented[mask][:,[0,2,1]];flipped+=1
    return oriented,int(patch.max())+1,reversed_faces,flipped


def self_intersecting_pairs(v,f):
    if not len(f):return np.zeros((0,2),np.int64)
    result=cgal.remesh_self_intersections(v,f,detect_only=True,first_only=False,cutoff=2_000_000)
    pairs=np.asarray(result[2],np.int64).reshape(-1,2)
    return np.unique(np.sort(pairs,axis=1),axis=0) if len(pairs) else pairs


def diagnose(v,f,intersections=True):
    boundary,nonmanifold=edge_incidence(f)
    used=np.unique(f) if len(f) else np.zeros(0,np.int64)
    links=np.asarray(igl.is_vertex_manifold(f),bool) if len(f) else np.zeros(0,bool)
    components=int(igl.facet_components(f)[0]) if len(f) else 0
    key=np.sort(f,axis=1) if len(f) else np.zeros((0,3),np.int64)
    _,count=np.unique(key,axis=0,return_counts=True) if len(f) else (None,np.zeros(0,int))
    normals=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]) if len(f) else np.zeros((0,3))
    oriented=np.asarray(igl.bfs_orient(f)[0],np.int64) if len(f) else f
    record=dict(vertices=int(len(v)),faces=int(len(f)),
        repeated_index_faces=int(((f[:,0]==f[:,1])|(f[:,1]==f[:,2])|(f[:,0]==f[:,2])).sum()) if len(f) else 0,
        zero_area_faces=int(np.all(normals==0,axis=1).sum()),
        duplicate_face_copies=int((count-1)[count>1].sum()),
        boundary_edges=boundary,nonmanifold_edges=nonmanifold,
        nonmanifold_vertices=int((~links[used]).sum()) if len(f) else 0,
        face_components=components,
        orientation_consistent=bool(len(f)) and bool(np.array_equal(oriented,f)),
        signed_volume_m3=signed_volume(v,f))
    record['closed']=boundary==0 and nonmanifold==0
    if intersections:
        pairs=self_intersecting_pairs(v,f)
        record['self_intersecting_face_pairs']=int(len(pairs))
        record['self_intersection_free']=not len(pairs)
    return record


def normalize(v,f,log,tag):
    """Steps 1-5 as one auditable block."""
    before=len(f)
    v,f=weld_exact(v,f)
    f,repeated=drop_repeated_index(f)
    f,zero=drop_zero_area(v,f)
    f,same,fins=dedupe_faces(f)
    f,patches,reversed_faces,flipped=orient_outward(v,f)
    log.append(dict(stage=tag,input_faces=before,welded_vertices=int(len(v)),
        repeated_index_faces_dropped=repeated,zero_area_faces_dropped=zero,
        duplicate_same_winding_faces_dropped=same,opposite_winding_fin_faces_dropped=fins,
        orientable_patches=patches,faces_rewound=reversed_faces,patches_flipped_outward=flipped,
        output_faces=int(len(f))))
    return np.ascontiguousarray(v),np.ascontiguousarray(f)


def repair(v,f):
    log=[];v,f=normalize(v,f,log,'weld_and_clean')
    pairs=self_intersecting_pairs(v,f)
    if len(pairs):
        base_v,base_f=v,f
        result=cgal.remesh_self_intersections(v,f,detect_only=False,first_only=False,
                                              stitch_all=True,cutoff=2_000_000)
        rv=np.ascontiguousarray(np.asarray(result[0],float));rf=np.ascontiguousarray(np.asarray(result[1],np.int64))
        log.append(dict(stage='resolve_self_intersections',input_face_pairs=int(len(pairs)),
                        arrangement_vertices=int(len(rv)),arrangement_faces=int(len(rf))))
        v,f=normalize(rv,rf,log,'weld_and_clean_after_arrangement')
        boundary,nonmanifold=edge_incidence(f)
        residual=len(self_intersecting_pairs(v,f))
        # The exact arrangement removes the crossings but keeps every overlapped sheet, so the
        # surface stays nonmanifold. Only the winding-number outer shell of the pre-arrangement
        # mesh is a single closed boundary; it is taken from step 5, which the boolean accepts.
        if nonmanifold or residual:
            try:
                uv,uf,_=cgal.mesh_boolean(base_v,base_f,EMPTY_V,EMPTY_F,type_str='union')
                uv=np.ascontiguousarray(np.asarray(uv,float));uf=np.ascontiguousarray(np.asarray(uf,np.int64))
                log.append(dict(stage='extract_outer_manifold',source='weld_and_clean',
                    reason='nonmanifold_edges' if nonmanifold else 'residual_self_intersections',
                    arrangement_nonmanifold_edges=nonmanifold,arrangement_boundary_edges=boundary,
                    arrangement_self_intersecting_pairs=residual,
                    input_faces=int(len(base_f)),output_faces=int(len(uf))))
                v,f=normalize(uv,uf,log,'weld_and_clean_after_union')
            except RuntimeError as error:
                log.append(dict(stage='extract_outer_manifold',source='weld_and_clean',failed=str(error),
                    arrangement_nonmanifold_edges=nonmanifold,arrangement_boundary_edges=boundary,
                    arrangement_self_intersecting_pairs=residual))
    used=np.unique(f) if len(f) else np.zeros(0,np.int64)
    if len(used)!=len(v):
        remap=np.full(len(v),-1,np.int64);remap[used]=np.arange(len(used))
        v=np.ascontiguousarray(v[used]);f=np.ascontiguousarray(remap[f.ravel()].reshape(-1,3))
        log.append(dict(stage='drop_unreferenced_vertices',removed=int(len(remap)-len(used))))
    return v,f,log


def self_test():
    tet_v=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]])
    tet_f=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]],np.int64)
    base=diagnose(tet_v,tet_f)
    assert base['closed'] and base['self_intersection_free'] and base['orientation_consistent']
    assert base['signed_volume_m3']>0,'reference tetrahedron must already be outward'
    seam_v=tet_v[tet_f].reshape(-1,3);seam_f=np.arange(12,dtype=np.int64).reshape(-1,3)
    v,f,log=repair(seam_v,seam_f)
    assert len(v)==4 and len(f)==4 and log[0]['welded_vertices']==4
    assert abs(signed_volume(v,f)-base['signed_volume_m3'])<1e-15,'welding must not change volume'
    inverted=repair(tet_v,tet_f[:,[0,2,1]])
    assert inverted[2][0]['patches_flipped_outward']==1 and signed_volume(*inverted[:2])>0
    dup=repair(tet_v,np.r_[tet_f,tet_f[:1],tet_f[:1]])
    assert dup[2][0]['duplicate_same_winding_faces_dropped']==2 and len(dup[1])==4
    fin=repair(tet_v,np.r_[tet_f,tet_f[:1],tet_f[:1][:,[0,2,1]]])
    assert fin[2][0]['opposite_winding_fin_faces_dropped']==2 and len(fin[1])==4
    assert diagnose(*fin[:2])['closed']
    repeated=repair(tet_v,np.r_[tet_f,np.array([[0,0,1]],np.int64)])
    assert repeated[2][0]['repeated_index_faces_dropped']==1 and diagnose(*repeated[:2])['closed']
    collinear_v=np.r_[tet_v,[[2.,0,0]]]
    collinear=repair(collinear_v,np.r_[tet_f,np.array([[0,1,4]],np.int64)])
    assert collinear[2][0]['zero_area_faces_dropped']==1 and len(collinear[1])==4
    cross_v=np.array([[-1.,-1,0],[1.,-1,0],[0,1.,0],[0,-.5,-1],[0,-.5,1],[0,.5,0]])
    assert len(self_intersecting_pairs(cross_v,np.array([[0,1,2],[3,4,5]],np.int64)))==1
    shift=np.array([.2,.2,.2])
    twin_v=np.r_[tet_v,tet_v+shift];twin_f=np.r_[tet_f,tet_f+4]
    assert not diagnose(twin_v,twin_f)['self_intersection_free']
    rv,rf,rlog=repair(twin_v,twin_f)
    after=diagnose(rv,rf)
    assert after['self_intersection_free'] and after['closed'] and after['nonmanifold_vertices']==0
    assert any(step['stage']=='resolve_self_intersections' for step in rlog)
    assert after['signed_volume_m3']<2*base['signed_volume_m3'],'self-union must remove the double-counted lens'
    disjoint=repair(np.r_[tet_v,tet_v+9],np.r_[tet_f,tet_f+4])
    assert diagnose(*disjoint[:2])['face_components']==2
    open_v,open_f=tet_v,tet_f[:3]
    assert diagnose(open_v,open_f)['boundary_edges']==3
    print('PASS weld/degeneracy/duplicate/fin/orientation/self-union repair and diagnostic checks')


def build(output,limit,order):
    out=Path(output).resolve()
    if out.exists():raise ValueError('Choose a fresh output directory')
    anatomy_path=ROOT/'data/derived/canonical/anatomy.json'
    readiness=ROOT/'data/derived/surface-volume-readiness-v1'
    anatomy=json.loads(anatomy_path.read_text())
    prior={r['entity_id']:r for r in map(json.loads,(readiness/'entities.jsonl').read_text().splitlines())}
    entities=[e for e in anatomy['entities'] if e['system']=='muscular']
    if not entities:raise ValueError('No muscular entities found')
    if order=='faces':entities.sort(key=lambda e:e['source_face_count'])
    if limit:entities=entities[:limit]
    packages={n:importlib.metadata.version(n) for n in ('libigl','numpy')}
    if packages['libigl']!='2.6.2':raise ValueError('This build pins libigl 2.6.2')
    out.mkdir(parents=True);(out/'geometry').mkdir();(out/'inputs').mkdir()
    shutil.copyfile(__file__,out/'inputs/build_muscle_tet_ready_surfaces.py')
    shutil.copyfile(readiness/'manifest.json',out/'inputs/surface-volume-readiness-manifest.json')
    native=Path(cgal.pyigl_copyleft_cgal.__file__);shutil.copyfile(native,out/'inputs'/native.name)
    sources={};records=[];started=time.monotonic()
    try:
        for entity in entities:
            ident=entity['id'];reference=entity['reference_geometry'];path=ROOT/reference['path']
            digest=sha(path)
            if digest!=reference['sha256']:raise ValueError('Canonical geometry changed: '+ident)
            sources[reference['path']]=digest
            raw=json.loads(gzip.decompress(path.read_bytes()))
            v=np.asarray(raw['positions'],float).reshape(-1,3)
            f=np.ascontiguousarray(np.asarray(raw['indices'],np.int64).reshape(-1,3))
            if raw.get('display_decimation'):raise ValueError('Decimated display geometry rejected: '+ident)
            if not np.isfinite(v).all() or f.min()<0 or f.max()>=len(v):
                raise ValueError('Invalid source geometry: '+ident)
            raw=diagnose(v,f,intersections=False)
            before=diagnose(*weld_exact(v,f))
            rv,rf,log=repair(v,f)
            after=diagnose(rv,rf)
            drift=after['signed_volume_m3']-before['signed_volume_m3']
            payload=dict(schema='ihm.muscle-tet-ready-surface.v1',entity_id=ident,name=entity['name'],
                         units='m',frame=reference['frame'],representation='triangular_surface',
                         source_path=reference['path'],source_sha256=digest,
                         positions=[float(x) for x in rv.ravel()],indices=[int(i) for i in rf.ravel()])
            geometry=out/'geometry'/(ident+'.json.gz')
            geometry.write_bytes(gzip.compress(json.dumps(payload,allow_nan=False).encode(),mtime=0))
            record=dict(entity_id=ident,name=entity['name'],source_path=reference['path'],source_sha256=digest,
                prior_blocking_reasons=prior.get(ident,{}).get('blocking_reasons'),
                prior_topological_candidate=prior.get(ident,{}).get('topological_candidate'),
                raw=raw,before=before,after=after,operations=log,
                output_path=str(geometry.relative_to(out)),output_sha256=sha(geometry),
                signed_volume_drift_m3=drift,
                signed_volume_relative_drift=drift/before['signed_volume_m3'] if before['signed_volume_m3'] else None,
                hole_filling_required=after['boundary_edges']>0,
                vertex_links_manifold=after['nonmanifold_vertices']==0,
                tet_ready=after['closed'] and after['self_intersection_free'] and
                          after['orientation_consistent'] and after['signed_volume_m3']>0)
            records.append(record)
            with (out/'entities.jsonl').open('a') as handle:handle.write(json.dumps(record,allow_nan=False)+'\n')
            if len(records)%25==0:print('Repaired',len(records),'of',len(entities),flush=True)
        def tally(key,side):return sum(bool(r[side][key]) for r in records)
        summary=dict(schema='ihm.muscle-tet-ready-build.v1',entities=len(records),
            stage_meaning=dict(raw='canonical rendering topology, seam-duplicated vertices, no intersection test',
                before='after lossless exact-coordinate welding only; this is the defect baseline',
                after='after the full documented repair order'),
            raw=dict(vertices=sum(r['raw']['vertices'] for r in records),
                closed=sum(r['raw']['closed'] for r in records),
                total_boundary_edges=sum(r['raw']['boundary_edges'] for r in records),
                total_nonmanifold_edges=sum(r['raw']['nonmanifold_edges'] for r in records)),
            faces_in=sum(r['before']['faces'] for r in records),faces_out=sum(r['after']['faces'] for r in records),
            before=dict(closed=tally('closed','before'),self_intersection_free=tally('self_intersection_free','before'),
                orientation_consistent=tally('orientation_consistent','before'),
                with_duplicate_faces=sum(r['before']['duplicate_face_copies']>0 for r in records),
                with_zero_area_faces=sum(r['before']['zero_area_faces']>0 for r in records),
                with_repeated_index_faces=sum(r['before']['repeated_index_faces']>0 for r in records),
                with_nonmanifold_edges=sum(r['before']['nonmanifold_edges']>0 for r in records),
                with_nonmanifold_vertices=sum(r['before']['nonmanifold_vertices']>0 for r in records),
                with_boundary_edges=sum(r['before']['boundary_edges']>0 for r in records),
                multi_component=sum(r['before']['face_components']>1 for r in records),
                total_duplicate_face_copies=sum(r['before']['duplicate_face_copies'] for r in records),
                total_zero_area_faces=sum(r['before']['zero_area_faces'] for r in records),
                total_repeated_index_faces=sum(r['before']['repeated_index_faces'] for r in records),
                total_nonmanifold_edges=sum(r['before']['nonmanifold_edges'] for r in records),
                total_boundary_edges=sum(r['before']['boundary_edges'] for r in records),
                total_self_intersecting_face_pairs=sum(r['before']['self_intersecting_face_pairs'] for r in records)),
            after=dict(closed=tally('closed','after'),self_intersection_free=tally('self_intersection_free','after'),
                orientation_consistent=tally('orientation_consistent','after'),
                with_duplicate_faces=sum(r['after']['duplicate_face_copies']>0 for r in records),
                with_zero_area_faces=sum(r['after']['zero_area_faces']>0 for r in records),
                with_repeated_index_faces=sum(r['after']['repeated_index_faces']>0 for r in records),
                with_nonmanifold_edges=sum(r['after']['nonmanifold_edges']>0 for r in records),
                with_nonmanifold_vertices=sum(r['after']['nonmanifold_vertices']>0 for r in records),
                with_boundary_edges=sum(r['after']['boundary_edges']>0 for r in records),
                multi_component=sum(r['after']['face_components']>1 for r in records),
                total_self_intersecting_face_pairs=sum(r['after']['self_intersecting_face_pairs'] for r in records)),
            tet_ready=sum(r['tet_ready'] for r in records),
            tet_ready_with_manifold_vertex_links=sum(r['tet_ready'] and r['vertex_links_manifold'] for r in records),
            pinched_vertex_links=[r['entity_id'] for r in records if not r['vertex_links_manifold']],
            hole_filling_required=[r['entity_id'] for r in records if r['hole_filling_required']],
            self_union_applied=[r['entity_id'] for r in records
                                if any(s['stage']=='extract_outer_manifold' for s in r['operations'])],
            arrangement_applied=sum(any(s['stage']=='resolve_self_intersections' for s in r['operations']) for r in records),
            max_absolute_relative_volume_drift=max((abs(r['signed_volume_relative_drift']) for r in records
                                                    if r['signed_volume_relative_drift'] is not None),default=0.0),
            elapsed_s=time.monotonic()-started)
        write(out/'summary.json',summary)
        artifacts={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
        write(out/'manifest.json',dict(schema='ihm.muscle-tet-ready-build.v1',packages=packages,python=sys.version,
            native_module_sha256=sha(native),anatomy_sha256=sha(anatomy_path),
            readiness_manifest_sha256=sha(readiness/'manifest.json'),
            source_geometry_sha256=sources,artifacts_sha256=artifacts,
            method='exact weld; drop repeated-index and exactly-collinear faces; collapse coincident faces '
                   '(opposite-winding pairs removed as zero-volume fins); bfs_orient patches flipped to positive '
                   'signed volume; CGAL remesh_self_intersections(stitch_all) arrangement; CGAL self-union outer '
                   'shell only where the arrangement leaves nonmanifold edges',
            canonical_geometry_modified=False,limitations=[
                'The self-union step discards material enclosed twice by overlapping lobes of one surface; '
                'its volume change is recorded per entity and is not a modelling decision.',
                'No hole is filled and no existing coordinate is moved; entities still carrying boundary edges '
                'are flagged hole_filling_required and are not certified.',
                'Closed, oriented, intersection-free shells are geometric candidates only: shell nesting, '
                'inter-entity overlap, exclusive material ownership and tissue occupancy remain unresolved.',
                'Tetrahedralisation itself is proved separately by scripts/verify_muscle_tet_ready_surfaces.py.']))
        print(json.dumps({k:v for k,v in summary.items()
                          if k not in ('hole_filling_required','self_union_applied','pinched_vertex_links')},indent=2))
    except BaseException as error:
        write(out/'failure.json',dict(type=type(error).__name__,message=str(error),completed=len(records)))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--self-test',action='store_true');p.add_argument('--output',type=Path)
    p.add_argument('--limit',type=int,default=0);p.add_argument('--order',choices=('anatomy','faces'),default='anatomy')
    a=p.parse_args()
    if a.self_test:self_test()
    if a.output:build(a.output,a.limit,a.order)
    if not a.self_test and not a.output:p.error('Select --self-test or --output')
