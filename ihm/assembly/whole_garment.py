"""Exact retained garment meshes as source-verified finite-mass elastic owners."""
from collections import defaultdict
from pathlib import Path
import hashlib,json
import numpy as np
from .clothing import Cloth


def audit_pattern(positions_m,triangles):
    x=np.asarray(positions_m,float);f=np.asarray(triangles)
    if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all() or f.ndim!=2 or f.shape[1]!=3 or not len(f) or f.dtype.kind not in 'iu' or f.min()<0 or f.max()>=len(x):raise ValueError('Invalid indexed pattern')
    area=.5*np.linalg.norm(np.cross(x[f[:,1]]-x[f[:,0]],x[f[:,2]]-x[f[:,0]]),axis=1)
    if np.any(area<=1e-16):raise ValueError('Degenerate pattern face')
    if len(np.unique(np.sort(f,axis=1),axis=0))!=len(f):raise ValueError('Duplicate pattern face')
    incidences=defaultdict(list);adj=[[] for _ in f]
    for i,t in enumerate(f):
        for a,b in zip(t,np.roll(t,-1)):incidences[tuple(sorted((int(a),int(b))))].append((i,int(a),int(b)))
    boundary=defaultdict(list)
    for edge,uses in incidences.items():
        if len(uses)>2:raise ValueError('Nonmanifold pattern edge')
        if len(uses)==2:
            a,b=uses
            if a[1:]!=b[1:][::-1]:raise ValueError('Inconsistent pattern winding')
            adj[a[0]].append(b[0]);adj[b[0]].append(a[0])
        else:
            a,b=edge;boundary[a].append(b);boundary[b].append(a)
    if any(len(v)!=2 for v in boundary.values()):raise ValueError('Ambiguous open boundary vertex')
    unseen=set(range(len(f)));components=0
    while unseen:
        stack=[unseen.pop()];components+=1
        while stack:
            for v in adj[stack.pop()]:
                if v in unseen:unseen.remove(v);stack.append(v)
    loops=[];unseen=set(boundary)
    while unseen:
        first=min(unseen);loop=[];previous=None;current=first
        while True:
            if current not in unseen:raise ValueError('Boundary revisits a different loop')
            unseen.remove(current);loop.append(current)
            nxt=min(boundary[current]) if previous is None else next(v for v in boundary[current] if v!=previous)
            previous,current=current,nxt
            if current==first:break
        loops.append(loop)
    if len(np.unique(f))!=len(x):raise ValueError('Unused pattern nodes')
    return {'vertices':len(x),'triangles':len(f),'edges':len(incidences),'face_components':components,
            'area_m2':float(area.sum()),'boundary_loops':loops,'interior_shared_edges':sum(len(v)==2 for v in incidences.values()),
            'euler_characteristic':int(len(x)-len(incidences)+len(f)),
            'self_intersections_checked':False,'seam_semantics':'Shared indices enforce exact mesh connectivity; stitched seam stiffness and strength are not identified.'}


def materialize_garments(root,*,areal_density_kg_m2,edge_stiffness_n_m):
    """No hidden density/stiffness defaults and no body support/contact assumptions."""
    root=Path(root);path=root/'data/derived/clothing/garments.json';raw=path.read_bytes();data=json.loads(raw)
    source=data['source'];relative=Path(source['path'])
    if relative.is_absolute() or '..' in relative.parts:raise ValueError('Unsafe garment source path')
    if hashlib.sha256((root/relative).read_bytes()).hexdigest()!=source['geometry_sha256']:raise ValueError('Canonical skin source hash mismatch')
    if hashlib.sha256((root/'app/src/clothing.js').read_bytes()).hexdigest()!=data['constructor_sha256']:raise ValueError('Garment constructor hash mismatch; materialize a fresh source artifact explicitly')
    bodies={}
    for g in data['garments']:
        name=g['id']
        if name not in ('shirt','shorts') or name in bodies:raise ValueError('Unexpected or duplicate garment owner')
        if g['source']['id']!=source['id'] or g['source']['geometry_sha256']!=source['geometry_sha256']:raise ValueError('Per-garment source identity mismatch')
        x=np.asarray(g['positions'],float).reshape(-1,3);tri=np.asarray(g['indices']).reshape(-1,3);audit=audit_pattern(x,tri)
        if audit['face_components']!=1 or len(audit['boundary_loops'])!={'shirt':4,'shorts':3}[name]:raise ValueError('Unexpected garment topology')
        loops=audit['boundary_loops'];ordered=sorted(range(len(loops)),key=lambda i:np.mean(x[loops[i],1]))
        labels={}
        if name=='shirt':
            labels[ordered[0]]='lower_hem';labels[ordered[-1]]='neck_opening'
            for i in ordered[1:-1]:labels[i]='left_arm_opening' if np.mean(x[loops[i],0])>0 else 'right_arm_opening'
        else:
            labels[ordered[-1]]='waistband'
            for i in ordered[:-1]:labels[i]='left_leg_hem' if np.mean(x[loops[i],0])>0 else 'right_leg_hem'
        audit['opening_labels']={str(k):v for k,v in labels.items()};audit['opening_classification_evidence']='Engineering pattern/frame labels; not measured skin/body contact'
        body=Cloth(x,tri,areal_density_kg_m2=areal_density_kg_m2,edge_stiffness_n_m=edge_stiffness_n_m)
        body.material_owner_ids=('garment-'+name,);body.pattern_audit=audit;body.source_node_indices=np.arange(len(x));body.source_face_indices=np.arange(len(tri));bodies[name]=body
    if set(bodies)!={'shirt','shorts'}:raise ValueError('Whole garment set missing')
    return bodies,{'garments_sha256':hashlib.sha256(raw).hexdigest(),'skin':source,'constructor_sha256':data['constructor_sha256'],
                   'areal_density_kg_m2':float(areal_density_kg_m2),'edge_stiffness_n_m':float(edge_stiffness_n_m),'parameter_status':'Caller-supplied engineering parameters; no calibrated fabric claim'}
