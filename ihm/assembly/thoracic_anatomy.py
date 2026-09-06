"""Source-bound reduced thoracic geometry and inertia, without passive laws."""
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from .cervical_inertia import _physical


def points(value):
    a=np.asarray(value,dtype=float)
    if a.ndim!=2 or a.shape[1]!=3 or not len(a) or len(a)>250000 or not np.isfinite(a).all():
        raise ValueError('Bounded finite 3D points required')
    return a


def surface_mass_prior(vertices,faces,mass_kg):
    v=points(vertices);f=np.asarray(faces,dtype=int)
    if f.ndim!=2 or f.shape[1]!=3 or f.size==0 or f.min()<0 or f.max()>=len(v) or len(f)>500000:
        raise ValueError('Valid bounded triangle indices required')
    mass=float(mass_kg)
    if not np.isfinite(mass) or mass<=0:raise ValueError('Positive mass allocation required')
    # Exact triangle-lamina integration. No closed volume or uniform bulk density
    # is inferred from open anatomical surfaces; total mass is supplied explicitly.
    origin=v.mean(0);area=0.;first=np.zeros(3);second=np.zeros((3,3));degenerate=0
    for begin in range(0,len(f),20000):
        tri=v[f[begin:begin+20000]]-origin
        a,b,c=np.moveaxis(tri,1,0);weights=np.linalg.norm(np.cross(b-a,c-a),axis=1)/2
        degenerate+=int(np.sum(weights<=1e-20));s=a+b+c
        area+=weights.sum();first+=np.einsum('i,ij->j',weights,s)/3
        moment=np.einsum('i,ij,ik->jk',weights,s,s)
        for t in (a,b,c):moment+=np.einsum('i,ij,ik->jk',weights,t,t)
        second+=moment/12
    if area<=1e-15:raise ValueError('Zero-area material surface')
    center=first/area;cov=second/area-np.outer(center,center)
    moment=mass*cov;inertia=np.trace(moment)*np.eye(3)-moment
    _physical(inertia)
    return {'mass_kg':mass,'center_m':(center+origin).tolist(),'inertia_kg_m2':inertia.tolist(),
        'source_area_m2':float(area),'degenerate_faces':degenerate,'areal_mass_kg_m2':mass/float(area),
        'basis':'Exact uniform triangular-lamina moments normalized to explicit source proxy mass; no closed-volume or measured shell-thickness claim'}


def nearest_node_pair(source,target):
    a=points(source);b=points(target);dist,index=cKDTree(b).query(a)
    i=int(np.argmin(dist));j=int(index[i])
    return {'source_vertex':i,'target_vertex':j,'source_point_m':a[i].tolist(),'target_point_m':b[j].tolist(),
        'gap_m':float(dist[i]),'basis':'Closest retained source vertices; not continuous-surface distance or measured enthesis'}


def hinge_candidate(rib,vertebra):
    rib=points(rib);vertebra=points(vertebra)
    distances,_=cKDTree(vertebra).query(rib)
    # Nearby posterior rib patch identifies a reproducible geometric axis proxy.
    count=max(12,int(np.ceil(.05*len(rib))))
    selected=np.argsort(distances,kind='stable')[:count];patch=rib[selected]
    origin=patch.mean(0);_,s,vh=np.linalg.svd(patch-origin,full_matrices=False)
    axis=vh[0]
    # Deterministic orientation only; positive q is not labeled inspiration.
    if axis[np.argmax(abs(axis))]<0:axis=-axis
    if s[0]<1e-8:raise ValueError('Degenerate posterior rib patch')
    return {'origin_m':origin.tolist(),'axis':axis.tolist(),'source_vertices':selected.tolist(),
        'posterior_patch_fraction':.05,'patch_singular_values_m':s.tolist(),
        'first_to_second_axis_ratio':float(s[0]/max(s[1],1e-15)),
        'max_patch_distance_to_vertebra_nodes_m':float(distances[selected].max()),
        'nearest_pair':nearest_node_pair(rib,vertebra),
        'basis':'Geometric PCA of nearest5percent rib vertices to matched vertebra; inferred hinge, not anatomical costovertebral/costotransverse calibration'}


def hinge_motion(vertices,origin,axis,angle_rad):
    v=points(vertices);o=np.asarray(origin,float);a=np.asarray(axis,float)
    if o.shape!=(3,) or a.shape!=(3,) or not np.isfinite([o,a]).all() or not np.isclose(np.linalg.norm(a),1,rtol=0,atol=1e-10) or not np.isfinite(angle_rad):
        raise ValueError('Finite proper hinge required')
    rotated=(v-o)@Rotation.from_rotvec(a*angle_rad).as_matrix().T
    return o+rotated,np.cross(a,rotated)


def descent_weights(vertices,anchor_indices):
    v=points(vertices);indices=np.asarray(anchor_indices,int)
    if indices.ndim!=1 or len(indices)<2 or indices.min()<0 or indices.max()>=len(v):raise ValueError('Valid diaphragm anchors required')
    distance=cKDTree(v[indices]).query(v)[0];maximum=float(distance.max())
    if maximum<=1e-12:raise ValueError('No free diaphragm span')
    return (distance/maximum)**2


def diaphragm_motion(vertices,weights,axis,descent_m):
    v=points(vertices);w=np.asarray(weights,float);a=np.asarray(axis,float)
    if w.shape!=(len(v),) or a.shape!=(3,) or not np.isfinite(w).all() or not np.isfinite(a).all() or (w<0).any() or (w>1).any() or not np.isclose(np.linalg.norm(a),1,atol=1e-10) or not np.isfinite(descent_m):
        raise ValueError('Invalid fixed-anchor diaphragm material map')
    jacobian=w[:,None]*a
    return v+descent_m*jacobian,jacobian


def validate_manifest(path):
    """Check frozen source bindings and independent composed mass ledger."""
    import gzip,hashlib,json
    from pathlib import Path
    from .cervical_inertia import combine_bodies
    path=Path(path);payload=json.loads(path.read_bytes());folder=path.parent
    if payload['native_activation_allowed'] is not False:raise ValueError('Anatomy recipe cannot authorize native mechanics')
    vertices={};faces={}
    for ident,row in payload['entities'].items():
        raw=(folder/row['source_geometry']['retained_copy']).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=row['source_geometry']['sha256']:raise ValueError('Source geometry hash mismatch')
        geometry=json.loads(gzip.decompress(raw));v=points(np.asarray(geometry['positions']).reshape(-1,3));f=np.asarray(geometry['indices']).reshape(-1,3)
        transform=np.asarray(row['source_geometry_to_torso']);vertices[ident]=v@transform[:3,:3].T+transform[:3,3];faces[ident]=f
    pairs=list(payload['attachment_candidates'])+list(payload['diaphragm_attachment_candidates'])
    pairs += [h['attachment'] for h in payload['rib_hinges'].values()]
    pairs += [p[k] for p in payload['intercostal_candidates'] for k in ('upper_candidate','lower_candidate')]
    for pair in pairs:
        for endpoint in ('a','b'):
            b=pair[endpoint];f=faces[b['entity_id']];v=vertices[b['entity_id']]
            bary=np.asarray(b['barycentric']);position=bary@v[f[b['face_index']]]
            if not np.isclose(bary.sum(),1) or (bary<0).any() or not np.allclose(position,b['point_in_torso_m'],rtol=0,atol=1e-12) or not np.allclose(v[b['vertex_index']],position,rtol=0,atol=1e-12):
                raise ValueError('Attachment is not bound to declared source face/node')
        gap=np.linalg.norm(np.asarray(pair['a']['point_in_torso_m'])-pair['b']['point_in_torso_m'])
        if not np.isclose(gap,pair['gap_m'],rtol=0,atol=1e-12):raise ValueError('Attachment gap mismatch')
    raw=(folder/payload['diaphragm_mode']['weights_file']).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=payload['diaphragm_mode']['weights_sha256']:raise ValueError('Diaphragm map identity mismatch')
    mode=json.loads(gzip.decompress(raw));w=np.asarray(mode['weights'])
    if not np.array_equal(w[mode['anchor_vertices']],np.zeros(len(mode['anchor_vertices']))):raise ValueError('Diaphragm anchors move in fixed-anchor mode')
    prior_entry=payload['inputs'][0];repo=Path(__file__).resolve().parents[2];prior_raw=(repo/prior_entry['path']).read_bytes()
    if hashlib.sha256(prior_raw).hexdigest()!=prior_entry['sha256']:raise ValueError('Prior partition identity mismatch')
    prior=json.loads(prior_raw)['residual_torso']
    combined=combine_bodies([payload['replace_reduced_torso_with']]+[r['torso_frame_mass_prior'] for r in payload['entities'].values() if r['kind']!='posterior_support'])
    for key in ('mass_kg','center_m','inertia_kg_m2'):
        if not np.allclose(combined[key],prior[key],rtol=0,atol=1e-10):raise ValueError('Thoracic mass debit violates prior torso ledger')
    return {'source_entities':len(vertices),'verified_attachment_pairs':len(pairs),'conserved_parent_mass_kg':combined['mass_kg']}
