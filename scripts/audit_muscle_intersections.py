"""Read-only CGAL intersection and vertex-link audit of retained muscle bulk.

Run with the isolated libigl environment, not the physiological runtime:
data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/audit_muscle_intersections.py --self-test
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


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def inspect(vertices,triangles):
    x=np.asarray(vertices,dtype=np.float64)
    t=np.asarray(triangles,dtype=np.int64)
    if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all():
        raise ValueError('Finite 3D vertices required')
    if t.ndim!=2 or t.shape[1]!=3 or not len(t) or t.min()<0 or t.max()>=len(x):
        raise ValueError('Nonempty valid triangle indices required')
    points=x[t]
    if np.any(np.all(np.cross(points[:,1]-points[:,0],points[:,2]-points[:,0])==0,axis=1)):
        raise ValueError('Degenerate triangles are not passed to CGAL')
    # Numeric coordinate equality changes neither vertex positions nor faces;
    # the source coordinate array and original triangle rows are retained.
    welded,inverse=np.unique(x,axis=0,return_inverse=True)
    faces=np.ascontiguousarray(inverse[t],dtype=np.int64)
    used=np.unique(faces)
    links=np.asarray(igl.is_vertex_manifold(faces),dtype=bool)
    invalid=used[~links[used]]
    result=cgal.remesh_self_intersections(welded,faces,detect_only=True,
                                         first_only=False,cutoff=2_000_000)
    pairs=np.asarray(result[2],dtype=np.int64).reshape(-1,2)
    if len(pairs):
        pairs=np.unique(np.sort(pairs,axis=1),axis=0)
        if pairs.min()<0 or pairs.max()>=len(t) or np.any(pairs[:,0]==pairs[:,1]):
            raise ValueError('Invalid intersection face identities')
    return dict(intersecting_face_pairs=pairs,nonmanifold_welded_vertices=invalid,
                source_to_welded_vertex=inverse,welded_vertices_m=welded,
                self_intersection_free=not len(pairs),vertex_links_manifold=not len(invalid))


def self_test():
    x=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.]])
    f=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]])
    clean=inspect(x,f)
    assert clean['self_intersection_free'] and clean['vertex_links_manifold']
    # Shared ordinary edges are not intersections; duplicate raw seam vertices
    # must not turn a proper tetrahedron into four disconnected triangles.
    seam=inspect(x[f].reshape(-1,3),np.arange(12).reshape(-1,3))
    assert seam['self_intersection_free'] and seam['vertex_links_manifold']
    crossing=np.array([[-1.,-1.,0],[1.,-1.,0],[0,1.,0],
                       [0,-.5,-1],[0,-.5,1],[0,.5,0]])
    pair=np.array([[0,1,2],[3,4,5]])
    assert np.array_equal(inspect(crossing,pair)['intersecting_face_pairs'],[[0,1]])
    coplanar=np.array([[0.,0,0],[2.,0,0],[0,2.,0],[.2,.2,0],[1,.2,0],[.2,1,0]])
    assert not inspect(coplanar,pair)['self_intersection_free']
    # Two disjoint solid interiors meeting at a single vertex have a bad link.
    pinched_x=np.r_[x,-x[1:]];pinched_f=np.r_[f,np.where(f==0,0,f+3)]
    assert not inspect(pinched_x,pinched_f)['vertex_links_manifold']
    # Nesting is deliberately NOT inferred from no intersections.
    nested=inspect(np.r_[x,x*.1+.1],np.r_[f,f+4])
    assert nested['self_intersection_free'] and nested['vertex_links_manifold']
    bad=x.copy();bad[3]=bad[0]
    try:inspect(bad,f)
    except ValueError:pass
    else:raise AssertionError('Degenerate geometry reached native intersection code')
    print('PASS native CGAL crossing/coplanar, seam/adjacency, vertex-link, nesting distinction and degeneracy checks')


def validate_collection(directory):
    directory=Path(directory).resolve();m=json.loads((directory/'manifest.json').read_text())
    for name,digest in m['artifacts_sha256'].items():
        path=(directory/name).resolve()
        if not path.is_relative_to(directory) or sha(path)!=digest:
            raise ValueError('Derived muscle input changed: '+name)
    for name,digest in m.get('source_files_sha256',{}).items():
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT) or sha(path)!=digest:
            raise ValueError('Underlying muscle source changed: '+name)
    return m


def build(output,parent,residual):
    out=Path(output).resolve();parent=Path(parent).resolve();residual=Path(residual).resolve()
    if out.exists():raise ValueError('Choose a fresh audit directory')
    pm=validate_collection(parent);rm=validate_collection(residual)
    if rm['parent_manifest_sha256']!=sha(parent/'manifest.json'):
        raise ValueError('Residual and parent partitions do not match')
    if sha(ROOT/'data/derived/canonical/anatomy.json')!=pm['anatomy_sha256']:
        raise ValueError('Canonical anatomical identity changed')
    out.mkdir(parents=True);(out/'inputs').mkdir();(out/'entities').mkdir()
    shutil.copyfile(__file__,out/'inputs/audit_muscle_intersections.py')
    for name,path in [('parent-manifest.json',parent/'manifest.json'),
                      ('residual-manifest.json',residual/'manifest.json'),
                      ('entities.jsonl',parent/'entities.jsonl'),('sources.jsonl',parent/'sources.jsonl')]:
        shutil.copyfile(path,out/'inputs'/name)
    native=Path(cgal.pyigl_copyleft_cgal.__file__)
    shutil.copyfile(native,out/'inputs'/native.name)
    packages={n:importlib.metadata.version(n) for n in ('libigl','numpy','scipy')}
    if packages['libigl']!='2.6.2':raise ValueError('This audit pins libigl 2.6.2')
    sources={r['source_sha256']:r for r in map(json.loads,(parent/'sources.jsonl').read_text().splitlines())}
    records=[];started=time.monotonic()
    try:
        for entity in map(json.loads,(parent/'entities.jsonl').read_text().splitlines()):
            ident=entity['entity_id'];source=sources[entity['source_record_sha256']]
            rp=residual/'entities'/ident/'report.json'
            if rp.exists():
                previous=json.loads(rp.read_text());candidate=previous['bulk_candidate']
                geometry=rp.with_name('bulk.json.gz');mapping=rp.with_name('mapping.json')
            else:
                candidate=source['bulk_candidate'];previous=source
                geometry=parent/source['representations']['bulk']['path'];mapping=parent/source['mapping']
            record=dict(entity_id=ident,name=entity['name'],prior_bulk_candidate=candidate,
                        input_geometry_path=str(geometry.relative_to(ROOT)),input_geometry_sha256=sha(geometry),
                        input_mapping_path=str(mapping.relative_to(ROOT)),input_mapping_sha256=sha(mapping),
                        geometric_screen_passed=False,canonical_handoff_applied=False,mass_assigned_kg=None)
            dest=out/'entities'/ident;dest.mkdir()
            shutil.copyfile(geometry,dest/'input-bulk.json.gz');shutil.copyfile(mapping,dest/'input-mapping.json')
            write(dest/'prior-report.json',previous)
            if candidate:
                g=json.loads(gzip.decompress(geometry.read_bytes()))
                x=np.asarray(g['positions'],float).reshape(-1,3);t=np.asarray(g['indices'],np.int64).reshape(-1,3)
                result=inspect(x,t);original=np.asarray(json.loads(mapping.read_text())['bulk_source_face_indices'],np.int64)
                if len(original)!=len(t):raise ValueError('Face provenance length mismatch')
                np.savez_compressed(dest/'intersection-evidence.npz',
                    intersecting_bulk_face_pairs=result['intersecting_face_pairs'],
                    intersecting_original_face_pairs=original[result['intersecting_face_pairs']],
                    nonmanifold_welded_vertices=result['nonmanifold_welded_vertices'],
                    source_to_welded_vertex=result['source_to_welded_vertex'],
                    welded_vertices_m=result['welded_vertices_m'])
                record.update(self_intersection_free=result['self_intersection_free'],
                    vertex_links_manifold=result['vertex_links_manifold'],
                    intersection_pair_count=len(result['intersecting_face_pairs']),
                    nonmanifold_vertex_count=len(result['nonmanifold_welded_vertices']),
                    geometric_screen_passed=result['self_intersection_free'] and result['vertex_links_manifold'])
            else:record['skipped_reason']='Prior topology/degeneracy gate failed; no native query performed'
            write(dest/'report.json',record);records.append(record)
            with (out/'records.jsonl').open('a') as f:f.write(json.dumps(record,allow_nan=False)+'\n')
            if len(records)%25==0:print('Audited',len(records),'muscles',flush=True)
        summary=dict(entities=len(records),prior_bulk_candidates=sum(r['prior_bulk_candidate'] for r in records),
                     screened_candidates=sum(r['geometric_screen_passed'] for r in records),
                     self_intersecting_entities=sum(r.get('intersection_pair_count',0)>0 for r in records),
                     nonmanifold_vertex_entities=sum(r.get('nonmanifold_vertex_count',0)>0 for r in records),
                     elapsed_s=time.monotonic()-started)
        write(out/'summary.json',summary)
        artifacts={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}
        write(out/'manifest.json',dict(schema='ihm.muscle-intersection-audit.v1',packages=packages,
              python=sys.version,native_module_sha256=sha(native),parent_manifest_sha256=sha(parent/'manifest.json'),
              residual_manifest_sha256=sha(residual/'manifest.json'),artifacts_sha256=artifacts,
              method='libigl 2.6.2 CGAL remesh_self_intersections detect_only; exact-coordinate seam welding; vertex-link manifold test',
              source_geometry_modified=False,limitations=[
                'Intersection-free surfaces do not establish shell nesting, outward/cavity convention or biological tissue occupancy.',
                'Cross-muscle and all-body overlaps, registration, material parameters and exclusive mass ownership remain unresolved.',
                'The library performs detection only; no remeshing, filling, normal changes or new anatomical surfaces are applied.']))
        print(json.dumps(summary,indent=2))
    except BaseException as e:
        write(out/'failure.json',dict(type=type(e).__name__,message=str(e),completed_entities=len(records)))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--self-test',action='store_true');p.add_argument('--output',type=Path)
    p.add_argument('--parent',type=Path,default=ROOT/'data/derived/muscle-dimensional-decomposition-v2')
    p.add_argument('--residual',type=Path,default=ROOT/'data/derived/remaining-muscle-bulk-v1');a=p.parse_args()
    if a.self_test:self_test()
    if a.output:build(a.output,a.parent,a.residual)
    if not a.self_test and not a.output:p.error('Select --self-test or --output')
