"""How far muscle surfaces can be decimated before volume, shape or meshability degrade.

Element count is the cost of the conforming route, so this measures what decimation buys and
what it costs: for a stratified sample of muscles at several reduction targets it records the
surviving topology, the divergence volume and area error, the sampled symmetric surface
deviation against the undecimated surface, the exact self-intersection count, and whether
TetGen still produces a conforming mesh.

Two environments are needed, so the run is staged through a work directory:
  .venv/bin/python scripts/measure_muscle_decimation_tolerance.py --phase decimate
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/measure_muscle_decimation_tolerance.py --phase measure
  .venv/bin/python scripts/measure_muscle_decimation_tolerance.py --phase report

Not measured here: mechanical response of a decimated mesh, decimation of anything but muscle,
and whether a decimated set is a valid global PLC.
"""
from pathlib import Path
import argparse,gzip,hashlib,importlib.util,json,sys,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
MUSCLE=ROOT/'data/derived/muscle-tet-ready-v1'
REDUCTIONS=(.5,.75,.9,.95,.98)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def surface(path):
    d=json.loads(gzip.decompress(Path(path).read_bytes()))
    return (np.ascontiguousarray(np.asarray(d['positions'],float).reshape(-1,3)),
            np.ascontiguousarray(np.asarray(d['indices'],np.int64).reshape(-1,3)))

def clean(V,F):
    """Exact-coordinate weld then drop repeated-index and zero-area faces; no coordinate is moved."""
    V=np.asarray(V,float);F=np.asarray(F,np.int64)
    u,inverse=np.unique(V,axis=0,return_inverse=True);F=np.ravel(inverse)[F]
    keep=(F[:,0]!=F[:,1])&(F[:,1]!=F[:,2])&(F[:,0]!=F[:,2]);dropped=int((~keep).sum());F=F[keep]
    A,B,C=u[F[:,0]],u[F[:,1]],u[F[:,2]]
    area=np.linalg.norm(np.cross(B-A,C-A),axis=1)/2;keep=area>0;dropped+=int((~keep).sum())
    return np.ascontiguousarray(u),np.ascontiguousarray(F[keep]),dropped

def topology(V,F):
    e=np.concatenate((F[:,[0,1]],F[:,[1,2]],F[:,[2,0]]))
    unique,inverse,counts=np.unique(np.sort(e,axis=1),axis=0,return_inverse=True,return_counts=True)
    inverse=np.ravel(inverse)
    net=np.rint(np.bincount(inverse,weights=np.where(e[:,0]<e[:,1],1.,-1.),minlength=len(unique))).astype(int)
    A,B,C=V[F[:,0]],V[F[:,1]],V[F[:,2]]
    return {'vertices':int(len(V)),'faces':int(len(F)),
        'boundary_edges':int((counts==1).sum()),'nonmanifold_edges':int((counts>2).sum()),
        'closed':bool((counts==1).sum()==0),'manifold':bool((counts==2).all()),
        'orientation_conflicts':int((net!=0).sum()),
        'signed_volume_m3':float(np.einsum('ij,ij->i',A,np.cross(B,C)).sum()/6),
        'surface_area_m2':float(np.linalg.norm(np.cross(B-A,C-A),axis=1).sum()/2)}

def decimate(work,sample,reductions):
    import fast_simplification
    work=Path(work);work.mkdir(parents=True,exist_ok=True)
    entities=[json.loads(l) for l in (MUSCLE/'entities.jsonl').read_text().splitlines()]
    volume={}
    for e in entities:
        V,F,_=clean(*surface(MUSCLE/e['output_path']));volume[e['entity_id']]=(abs(topology(V,F)['signed_volume_m3']),e)
    order=sorted(volume,key=lambda k:volume[k][0])
    picks=[order[i] for i in np.unique(np.linspace(0,len(order)-1,sample).astype(int))]
    rows=[];began=time.perf_counter()
    for eid in picks:
        e=volume[eid][1];V,F=surface(MUSCLE/e['output_path']);base=topology(V,F)
        np.savez(work/f'{eid}__base.npz',V=V,F=F)
        row={'entity_id':eid,'name':e['name'],'base':base,'levels':[]}
        for r in reductions:
            t=time.perf_counter()
            dv,df=fast_simplification.simplify(V.astype(np.float32),F.astype(np.int32),r)
            seconds=time.perf_counter()-t
            dv,df,dropped=clean(dv.astype(np.float64),df.astype(np.int64))
            top=topology(dv,df) if len(df) else None
            np.savez(work/f'{eid}__{r:g}.npz',V=dv,F=df)
            row['levels'].append({'target_reduction':r,'seconds':seconds,'topology':top,'degenerate_faces_dropped':dropped,
                'face_reduction_achieved':float(1-len(df)/base['faces']),
                'volume_relative_error':(float(abs(top['signed_volume_m3'])/abs(base['signed_volume_m3'])-1) if top else None),
                'area_relative_error':(float(top['surface_area_m2']/base['surface_area_m2']-1) if top else None)})
        rows.append(row)
    (work/'decimated.json').write_text(json.dumps({'schema':'ihm.muscle-decimation.decimate.v1',
        'reductions':list(reductions),'entities':rows,'seconds':time.perf_counter()-began,
        'library':'fast_simplification','python':sys.version},indent=2,allow_nan=False)+'\n')
    print(json.dumps({'entities':len(rows),'reductions':list(reductions),'seconds':time.perf_counter()-began},indent=2))

def measure(work):
    import igl,igl.copyleft.cgal as cgal
    sys.path.insert(0,str(ROOT/'scripts'))
    spec=importlib.util.spec_from_file_location('plc',ROOT/'scripts/probe_conforming_complement_plc.py')
    plc=importlib.util.module_from_spec(spec);argv=sys.argv;sys.argv=['x'];spec.loader.exec_module(plc);sys.argv=argv
    work=Path(work);payload=json.loads((work/'decimated.json').read_text());began=time.perf_counter()
    for row in payload['entities']:
        d=np.load(work/f"{row['entity_id']}__base.npz");V=np.ascontiguousarray(d['V']);F=np.ascontiguousarray(d['F'])
        qs=plc_sample(V,F)
        for level in row['levels']:
            e=np.load(work/f"{row['entity_id']}__{level['target_reduction']:g}.npz")
            dv=np.ascontiguousarray(e['V']);df=np.ascontiguousarray(e['F'])
            if not len(df):level['measurement']=None;continue
            qd=plc_sample(dv,df)
            d1=np.sqrt(igl.point_mesh_squared_distance(qs,dv,df)[0])
            d2=np.sqrt(igl.point_mesh_squared_distance(qd,V,F)[0])
            si=int(len(cgal.remesh_self_intersections(dv,df,True,False,False,False,100000)[2]))
            trial=plc.tetrahedralize(dv,df,'pY')
            level['measurement']={'symmetric_hausdorff_m':float(max(d1.max(),d2.max())),
                'symmetric_mean_m':float((d1.mean()+d2.mean())/2),
                'p95_original_to_decimated_m':float(np.percentile(d1,95)),
                'self_intersecting_face_pairs':si,
                'tetgen_succeeded':bool(trial.get('succeeded')),'tetgen_tets':trial.get('tets'),
                'tetgen_exception':trial.get('exception')}
    payload['measure_seconds']=time.perf_counter()-began
    (work/'measured.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'entities':len(payload['entities']),'seconds':payload['measure_seconds']},indent=2))

def plc_sample(V,F,points=20000):
    """Deterministic area-proportional barycentric samples; enough to bound deviation, not a certified Hausdorff."""
    A,B,C=V[F[:,0]],V[F[:,1]],V[F[:,2]]
    area=np.linalg.norm(np.cross(B-A,C-A),axis=1)/2
    n=np.maximum(1,np.rint(points*area/max(area.sum(),1e-30)).astype(int))
    rng=np.random.default_rng(0);out=[V]
    for k in np.unique(n):
        idx=np.flatnonzero(n==k)
        u=rng.random((len(idx),k));v=rng.random((len(idx),k));flip=u+v>1
        u=np.where(flip,1-u,u);v=np.where(flip,1-v,v)
        out.append((A[idx][:,None,:]+u[...,None]*(B-A)[idx][:,None,:]+v[...,None]*(C-A)[idx][:,None,:]).reshape(-1,3))
    return np.ascontiguousarray(np.concatenate(out))

def report(work,out,force,volume_tolerance,deviation_tolerance):
    work=Path(work);out=Path(out)
    if out.exists() and any(out.iterdir()) and not force:raise ValueError('Choose a fresh output directory or pass --force')
    payload=json.loads((work/'measured.json').read_text())
    levels={}
    for row in payload['entities']:
        for level in row['levels']:
            k=f"{level['target_reduction']:g}";s=levels.setdefault(k,{'target_reduction':level['target_reduction'],'entities':0,
                'closed':0,'manifold':0,'self_intersection_free':0,'tetgen_succeeded':0,'within_volume_tolerance':0,
                'within_deviation_tolerance':0,'within_both':0,'volume_errors':[],'hausdorff':[],'mean_dev':[],
                'faces_in':0,'faces_out':0,'seconds':0.})
            s['entities']+=1;s['faces_in']+=row['base']['faces'];s['seconds']+=level['seconds']
            m=level.get('measurement');t=level.get('topology')
            if t:
                s['faces_out']+=t['faces'];s['closed']+=int(t['closed'] and t['orientation_conflicts']==0)
                s['manifold']+=int(t['manifold'] and t['closed'] and t['orientation_conflicts']==0)
                s['volume_errors'].append(abs(level['volume_relative_error']))
            if m:
                s['self_intersection_free']+=int(m['self_intersecting_face_pairs']==0)
                s['tetgen_succeeded']+=int(m['tetgen_succeeded'])
                s['hausdorff'].append(m['symmetric_hausdorff_m']);s['mean_dev'].append(m['symmetric_mean_m'])
                v=abs(level['volume_relative_error'])<=volume_tolerance
                dd=m['symmetric_hausdorff_m']<=deviation_tolerance
                s['within_volume_tolerance']+=int(v);s['within_deviation_tolerance']+=int(dd)
                s['within_both']+=int(v and dd and t['closed'] and t['manifold'] and m['self_intersecting_face_pairs']==0 and m['tetgen_succeeded'])
    for s in levels.values():
        s['face_reduction_achieved']=float(1-s['faces_out']/s['faces_in'])
        for key,src in (('median_abs_volume_error','volume_errors'),('median_symmetric_hausdorff_m','hausdorff'),
                        ('median_symmetric_mean_m','mean_dev')):
            s[key]=float(np.median(s[src])) if s[src] else None
        s['p95_symmetric_hausdorff_m']=float(np.percentile(s['hausdorff'],95)) if s['hausdorff'] else None
        s['max_abs_volume_error']=float(max(s['volume_errors'])) if s['volume_errors'] else None
        for k in ('volume_errors','hausdorff','mean_dev'):s.pop(k)
    passing=[s for s in levels.values() if s['within_both']==s['entities']]
    summary={'schema':'ihm.muscle-decimation-tolerance.v1','entities':len(payload['entities']),
        'tolerance':{'abs_volume_relative_error':volume_tolerance,'symmetric_hausdorff_m':deviation_tolerance,
            'also_required':'closed, edge-manifold, consistently oriented, zero exact self-intersections, TetGen pY succeeds'},
        'levels':dict(sorted(levels.items(),key=lambda kv:kv[1]['target_reduction'])),
        'largest_reduction_meeting_tolerance_for_every_sampled_muscle':
            (max(s['target_reduction'] for s in passing) if passing else None),
        'library':'fast_simplification 0.2.0','per_entity':payload['entities'],
        'decimate_seconds':payload['seconds'],'measure_seconds':payload['measure_seconds'],
        'unverified':['surface deviation is a sampled bound, not a certified Hausdorff distance',
            'decimated surfaces are not checked for mutual intersection with neighbouring structures',
            'no mechanical solve is run on a decimated mesh',
            'only muscle surfaces are decimated; bone, vessel and organ surfaces are untouched']}
    out.mkdir(parents=True,exist_ok=True)
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (out/'manifest.json').write_text(json.dumps({'schema':'ihm.muscle-decimation-tolerance-manifest.v1',
        'inputs_sha256':{'data/derived/muscle-tet-ready-v1/manifest.json':sha(MUSCLE/'manifest.json'),
            'scripts/measure_muscle_decimation_tolerance.py':sha(__file__)},
        'artifacts_sha256':{'summary.json':sha(out/'summary.json')},'canonical_assets_modified':False},indent=2)+'\n')
    printable={k:v for k,v in summary.items() if k!='per_entity'}
    print(json.dumps(printable,indent=2,allow_nan=False))
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--phase',choices=('decimate','measure','report'),required=True)
    p.add_argument('--work',default=str(ROOT/'data/derived/muscle-decimation-tolerance-v1/work'))
    p.add_argument('--out',default=str(ROOT/'data/derived/muscle-decimation-tolerance-v1'))
    p.add_argument('--sample',type=int,default=20);p.add_argument('--force',action='store_true')
    p.add_argument('--volume-tolerance',type=float,default=.01);p.add_argument('--deviation-tolerance',type=float,default=.001)
    a=p.parse_args()
    if a.phase=='decimate':decimate(a.work,a.sample,REDUCTIONS)
    elif a.phase=='measure':measure(a.work)
    else:report(a.work,a.out,a.force,a.volume_tolerance,a.deviation_tolerance)
