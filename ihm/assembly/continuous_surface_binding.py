"""Shared graph-regularized skin attachment and work-conjugate force scatter.

This is a geometric linear-blend embedding, not a constitutive skin model.
Weights are solved once on retained rest topology, never from a displayed pose.
"""
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import cg
from scipy.spatial import cKDTree
from .surface_binding import SegmentSurfaceBinding, SKIN_ASSET, _digest

DEFAULT_ASSET = 'data/derived/canonical/continuous_surface_binding.json.gz'


def materialize(root, registration, *, screening_length_m=.08):
    if not np.isfinite(screening_length_m) or screening_length_m<=0: raise ValueError('Positive finite surface screening length required')
    source_bytes=Path(__file__).read_bytes()
    root=Path(root); hard=SegmentSurfaceBinding.from_root(root,registration)
    segments=hard.manifest()['segments']; sources=hard.manifest()['source_identity']
    sources['continuous_materialization_implementation_sha256']=hashlib.sha256(source_bytes).hexdigest()
    geometry=json.loads(gzip.decompress((root/SKIN_ASSET).read_bytes()))
    original=np.asarray(geometry['positions'],float).reshape(-1,3)
    # Exact coordinate duplicates share one material support. No proximity weld
    # may connect a touching hand to a hip or merge the two nested skin shells.
    points,inverse=np.unique(original,axis=0,return_inverse=True)
    triangles=inverse[np.asarray(geometry['indices']).reshape(-1,3)]
    edges=np.unique(np.sort(np.concatenate([triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]]]),axis=1),axis=0)
    edges=edges[edges[:,0]!=edges[:,1]]
    length=np.linalg.norm(points[edges[:,0]]-points[edges[:,1]],axis=1)
    conductance=1/np.maximum(length,.001)**2
    adjacency=coo_matrix((conductance,(edges[:,0],edges[:,1])),shape=(len(points),len(points))).tocsr()
    adjacency=adjacency+adjacency.T
    count,labels=connected_components(adjacency)
    distances=[]
    for segment in segments:
        bone_points=[]
        for bone in registration.groups[segment['id']]['canonical_bones']:
            relative='data/derived/canonical/geometry/'+bone+'.json.gz'
            raw=(root/relative).read_bytes();sources[relative]=hashlib.sha256(raw).hexdigest()
            bone_points.extend(json.loads(gzip.decompress(raw))['positions'])
        distances.append(cKDTree(np.asarray(bone_points).reshape(-1,3)).query(points)[0])
    seeds=np.argmin(np.stack(distances,axis=1),axis=1)
    coefficient=float(screening_length_m)**2/6
    degree=np.asarray(adjacency.sum(axis=1)).ravel()
    operator=diags(1+coefficient*degree)-coefficient*adjacency
    preconditioner=diags(1/operator.diagonal())
    weights=np.empty((len(points),len(segments)))
    residuals=[]
    for index in range(len(segments)):
        rhs=(seeds==index).astype(float)
        solution,status=cg(operator,rhs,M=preconditioner,rtol=1e-8,atol=1e-12,maxiter=5000)
        if status: raise RuntimeError('Surface graph solve failed: '+str(status))
        residuals.append(float(np.linalg.norm(operator@solution-rhs)))
        weights[:,index]=np.maximum(solution,0)
    # Tiny numerical tails are removed before both physics and display consume
    # the same quantized sidecar, with an explicit bound on the alteration.
    weights[weights<1e-8]=0
    weights/=weights.sum(axis=1)[:,None]
    weights=weights[inverse].astype(np.float32).astype(float)
    weights/=weights.sum(axis=1)[:,None]
    payload=dict(schema='ihm.continuous-surface-binding.v1', segments=segments,
        surface_entity_ids=['body-bp3d-FJ2810'], reference_positions_m=original.tolist(),
        weights=weights.tolist(), source_identity=sources,
        materialization_source={'path':'ihm/assembly/continuous_surface_binding.py',
            'sha256':sources['continuous_materialization_implementation_sha256'],
            'text':source_bytes.decode()},
        algorithm=dict(seed='nearest_retained_bone_mesh_vertex',tie_break='lexicographic_native_segment_id',
            graph='retained skin triangle edges, exact-coordinate duplicate weld only',
            screening_length_m=float(screening_length_m), minimum_edge_length_m=.001,
            equation='(I + screening_length_m^2 / 6 * weighted_graph_laplacian) weights = one_hot_bone_seed',
            weight_tail_cutoff=1e-8, quantization='float32 then row normalization in float64',
            connected_components=count, largest_component_vertices=sorted(np.bincount(labels).tolist(),reverse=True)[:10],
            duplicate_vertices_welded=len(original)-len(points), maximum_solver_residual=max(residuals)),
        scope='Approximate rest-topology graph-regularized linear blend attachment shared by skin display and world contact. No skin FEM, material stiffness calibration, volume preservation, respiration deformation, or joint-continuity certification.')
    payload['binding_identity']=_digest(payload)
    return payload


class ContinuousSurfaceBinding:
    def __init__(self,payload):
        if not isinstance(payload,dict) or payload.get('schema')!='ihm.continuous-surface-binding.v1':
            raise ValueError('Unknown continuous surface binding schema')
        content=deepcopy(payload); identity=content.pop('binding_identity')
        if _digest(content)!=identity: raise ValueError('Surface binding identity mismatch')
        source=payload.get('materialization_source')
        if source is not None:
            if not isinstance(source,dict) or not isinstance(source.get('text'),str):
                raise ValueError('Invalid retained materialization source')
            digest=hashlib.sha256(source['text'].encode()).hexdigest()
            if source.get('sha256')!=digest or payload['source_identity'].get('continuous_materialization_implementation_sha256')!=digest:
                raise ValueError('Materialization source identity mismatch')
        self._payload=deepcopy(payload);self.identity=identity
        self._asset_bytes=json.dumps(payload,separators=(',', ':'),allow_nan=False).encode()
        self._asset_sha256=hashlib.sha256(self._asset_bytes).hexdigest()
        self.segments=deepcopy(payload['segments']);self.rest=np.asarray(payload['reference_positions_m'],float)
        self.weights=np.asarray(payload['weights'],float)
        if not isinstance(self.segments,list) or not self.segments:
            raise ValueError('Surface binding needs segment supports')
        for segment in self.segments:
            if not isinstance(segment,dict) or any(not isinstance(segment.get(k),str) or not segment[k] for k in ('id','bone_id')):
                raise ValueError('Invalid surface segment identity')
            center=np.asarray(segment.get('reference_centroid_m'),float)
            if center.shape!=(3,) or not np.isfinite(center).all():raise ValueError('Invalid surface reference centroid')
        for key in ('id','bone_id'):
            if len({s[key] for s in self.segments})!=len(self.segments):raise ValueError('Duplicate surface segment support')
        self._validate_binding((self.rest,self.weights))
        if not len(self.rest):raise ValueError('Surface binding needs rest vertices')
        self.rest.setflags(write=False);self.weights.setflags(write=False)
        self._tree=None

    @classmethod
    def from_root(cls,root,registration,*,asset=DEFAULT_ASSET):
        root=Path(root).resolve();path=(root/asset).resolve()
        if not path.is_relative_to(root):raise ValueError('Surface binding asset escaped workspace')
        payload=json.loads(gzip.decompress(path.read_bytes()))
        # Registration coordinates and bone grouping must match the retained
        # support even when native mass, initial pose or controller differs.
        expected=SegmentSurfaceBinding(registration).manifest()['segments']
        if payload['segments']!=expected: raise ValueError('Surface binding registration supports changed')
        sources=payload['source_identity']
        if not {SKIN_ASSET,'data/derived/canonical/mechanics.json'}<=sources.keys():
            raise ValueError('Surface binding geometry provenance missing')
        for relative,digest in sources.items():
            if relative.startswith('data/'):
                source=(root/relative).resolve()
                if not source.is_relative_to(root):raise ValueError('Surface source escaped workspace')
                if hashlib.sha256(source.read_bytes()).hexdigest()!=digest:
                    raise ValueError('Surface binding geometry source changed: '+relative)
        return cls(payload)

    def manifest(self):
        return dict(schema=self._payload['schema'],binding_identity=self.identity,
            coordinate_frame='canonical_current_world',reference_coordinate_frame='canonical_rest',
            rule='graph_regularized_linear_blend',segments=deepcopy(self.segments),
            surface_entity_ids=deepcopy(self._payload['surface_entity_ids']),
            weights_url='/api/surface-binding/'+self._asset_sha256,weights_sha256=self._asset_sha256,
            vertex_count=len(self.rest),source_identity=deepcopy(self._payload['source_identity']),
            materialization_source_retained='materialization_source' in self._payload,
            algorithm=deepcopy(self._payload['algorithm']),scope=self._payload['scope'])

    def asset_records(self): return {self._asset_sha256:self._asset_bytes}

    def sidecar(self): return deepcopy(self._payload)

    def frame(self,entities):
        result={}
        for segment in self.segments:
            center,rotation=self._pose(entities[segment['bone_id']])
            result[segment['id']]={'centroid_m':center.tolist(),'rotation_matrix':rotation.tolist()}
        return result

    @staticmethod
    def _pose(pose):
        center=np.asarray(pose['centroid_m'],float);rotation=np.asarray(pose['rotation_matrix'],float)
        if (center.shape!=(3,) or rotation.shape!=(3,3) or not np.isfinite(center).all()
                or not np.isfinite(rotation).all() or not np.allclose(rotation.T@rotation,np.eye(3),atol=1e-8,rtol=0)
                or np.linalg.det(rotation)<=0):raise ValueError('Surface pose requires finite centroid and proper rotation')
        return center,rotation

    def _validate_binding(self,binding):
        points,weights=map(lambda value:np.asarray(value,float),binding)
        if points.ndim!=2 or points.shape[1]!=3 or not np.isfinite(points).all():
            raise ValueError('Invalid surface rest vertices')
        if (weights.shape!=(len(points),len(self.segments)) or not np.isfinite(weights).all()
                or np.any(weights<0) or not np.allclose(weights.sum(axis=1),1,atol=1e-12,rtol=0)):
            raise ValueError('Invalid continuous surface weights')
        return points,weights

    def bind(self,rest_vertices):
        """Transfer existing skin weights to other surfaces by nearest rest vertex.

        Deterministic exact-distance ties choose the smallest source vertex index.
        This is a declared extension prior, not independent garment dynamics.
        """
        points=np.asarray(rest_vertices,float)
        if points.ndim!=2 or points.shape[1]!=3 or not np.isfinite(points).all(): raise ValueError('Invalid rest surface vertices')
        if self._tree is None:self._tree=cKDTree(self.rest)
        distances,indices=self._tree.query(points)
        for i,(point,distance) in enumerate(zip(points,distances)):
            candidates=self._tree.query_ball_point(point,np.nextafter(distance,np.inf))
            if len(candidates)>1: indices[i]=min(candidates,key=lambda j:(float(np.linalg.norm(self.rest[j]-point)),j))
        return points.copy(), self.weights[indices].copy()

    def sampled_binding(self,indices):
        indices=np.asarray(indices)
        if indices.ndim!=1 or indices.dtype.kind not in 'iu' or np.any(indices<0) or np.any(indices>=len(self.rest)):
            raise ValueError('Invalid surface sample indices')
        return self.rest[indices].copy(),self.weights[indices].copy()

    def stations(self,binding,transforms):
        points,_=self._validate_binding(binding)
        stations=[]
        for segment in self.segments:
            center,rotation=self._pose(transforms[segment['id']])
            stations.append(center+(points-np.asarray(segment['reference_centroid_m']))@rotation.T)
        return np.stack(stations,axis=1)

    def project(self,binding,transforms):
        return np.sum(self.stations(binding,transforms)*binding[1][:,:,None],axis=1)

    def scatter_forces(self,binding,transforms,forces):
        """Exact transpose of the material-point velocity map.

        Apply w_j F at x_j, not at the blended point. Then total force, moment
        about every origin, and virtual power equal the surface port values.
        """
        forces=np.asarray(forces,float);points,weights=binding
        if forces.shape!=points.shape or not np.isfinite(forces).all():raise ValueError('Invalid material surface forces')
        stations=self.stations(binding,transforms);ports=[]
        for i in np.flatnonzero(np.linalg.norm(forces,axis=1)>0):
            for j in np.flatnonzero(weights[i]>0):
                ports.append(dict(id=self.segments[j]['bone_id'],point_m=stations[i,j].tolist(),force_n=(weights[i,j]*forces[i]).tolist()))
        return ports
