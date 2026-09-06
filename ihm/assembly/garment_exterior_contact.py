"""Exact exterior source selection and fail-closed garment tunneling guard.

Shared MovingSurfaceContact remains unchanged for other clients. This adapter
rejects events needing a contact schedule rather than pretending to resolve CCD.
"""
import copy
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from .skin_layers import physical_skin_support
from .garment_surface_contact import MovingSurfaceContact
from .garment_body_ports import SourceSurface
from .garment_swept_vertex_face import swept_vertex_face


def exterior_source(root,skin_reference,pairings,garment_positions_m):
    root=Path(root);relative=Path(skin_reference['path'])
    if relative.is_absolute() or '..' in relative.parts:raise ValueError('Unsafe skin source path')
    raw=(root/relative).read_bytes();evidence=(root/'data/research/engineered_skin_territories/materialization.json').read_bytes()
    reference=dict(path=str(relative),sha256=skin_reference['geometry_sha256'],units='m',frame='bodyparts3d-display-m')
    support=physical_skin_support(reference,raw,evidence)
    data=json.loads(gzip.decompress(raw));x=np.asarray(data['positions'],float).reshape(-1,3);tri=np.asarray(data['indices'],int).reshape(-1,3)
    face_ids=np.asarray(json.loads(evidence)['contact_eligible_triangle_ids'],int)
    if len(np.unique(face_ids))!=len(face_ids):raise ValueError('Duplicate exterior face identity')
    # This is a NEW candidate support mapping. Preserve the complete prior
    # record; source nodes are unchanged, inferred tether locations may change.
    result=copy.deepcopy(pairings);surface=SourceSurface(x,tri[face_ids])
    for name,rows in result['garments'].items():
        for category in ('support_pairings','contact_samples'):
            for i,previous in enumerate(rows[category]):
                point=np.asarray(garment_positions_m[name][previous['garment_node']],float)
                current=surface.closest(point);local=current.pop('triangle_index')
                rows[category][i]={**previous,**current,'triangle_index':int(face_ids[local]),
                    'previous_raw_source_pairing':copy.deepcopy(previous),
                    'mapping_status':'new_exterior_only_engineering_candidate_not_original_tether',
                    'exterior_local_face_index':int(local)}
    result['schema']='ihm.garment-exterior-pairings.v1'
    result['mapping_identity_sha256']=hashlib.sha256(json.dumps(result,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return x,tri[face_ids],face_ids,result,support


class StrictGarmentSurfaceContact(MovingSurfaceContact):
    """Transactional guard; any unresolved event prevents garment acceptance."""
    def resolve(self,positions_m,velocities_m_s,mass_kg,dt_s):
        x=np.asarray(positions_m,float);v=np.asarray(velocities_m_s,float)
        if x.ndim!=2 or x.shape[1]!=3 or v.shape!=x.shape or not np.isfinite(x).all() or not np.isfinite(v).all():raise ValueError('Finite compatible cloth endpoints required')
        if not np.isfinite(dt_s) or dt_s<=0 or dt_s>self.duration or not 0<=self.fraction<=1:raise ValueError('Invalid swept interval')
        lower=self.fraction-dt_s/self.duration
        if lower< -1e-12:raise ValueError('Swept interval precedes source endpoints')
        a=x-dt_s*v # Cloth.step's actual kick/drift trial; linear within substep.
        start=self.start+max(0.,lower)*self.delta;end=self.start+self.fraction*self.delta
        for lo in range(0,len(x),128):
            hi=min(lo+128,len(x));mid=.5*(a[lo:hi]+x[lo:hi])
            radius=self.radius+self.motion_bound+.5*np.linalg.norm(x[lo:hi]-a[lo:hi],axis=1)+self.distance
            sizes=self.tree.query_ball_point(mid,radius,return_length=True)
            if int(sizes.sum())>self.max_pairs:raise MemoryError('Swept candidate budget exceeded; no contact acceptance')
            queries=self.tree.query_ball_point(mid,radius,return_sorted=True)
            for offset,candidates in enumerate(queries):
                node=lo+offset
                for face in candidates:
                    ids=self.triangles[face]
                    # Conservative swept AABB exclusion before polynomial work.
                    floor=np.minimum(start[ids].min(0),end[ids].min(0))-1e-12
                    ceiling=np.maximum(start[ids].max(0),end[ids].max(0))+1e-12
                    if np.any(np.maximum(a[node],x[node])<floor) or np.any(np.minimum(a[node],x[node])>ceiling):continue
                    event=swept_vertex_face(a[node],x[node],start[ids],end[ids])
                    if event['status']!='clear':
                        source_ids=getattr(self,'source_face_indices',None)
                        source_face=None if source_ids is None else int(source_ids[face])
                        raise ValueError('Unaccepted swept garment contact: '+json.dumps(dict(node=int(node),local_face=int(face),source_face=source_face,**event)))
        result=super().resolve(x,v,mass_kg,dt_s)
        if result['unresolved_edge_contacts']:raise ValueError('Unaccepted nearest-edge garment contact')
        result['acceptance_scope']='No detected swept interior node-face event; edge/vertex/self-contact and finite thickness not certified'
        return result
