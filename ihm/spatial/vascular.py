"""Finite-element audit of original tetrahedral CFD states, without unit inference.

Velocity and pressure are interpreted as continuous, piecewise-linear nodal
fields on the supplied tetrahedra. This is a postprocessing audit, not another
flow solver or a certificate of physical/clinical validity.
"""
import numpy as np
from scipy.spatial import cKDTree


def _points(value):
    p=np.asarray(value,dtype=float)
    if p.ndim!=2 or p.shape[1]!=3 or not len(p) or not np.isfinite(p).all():
        raise ValueError('Expected nonempty finite Nx3 source coordinates')
    return p


def _connectivity(value,width,count):
    raw=np.asarray(value)
    if raw.ndim!=2 or raw.shape[1]!=width or not len(raw) or not np.issubdtype(raw.dtype,np.integer):
        raise ValueError('Expected nonempty integer connectivity')
    f=raw.astype(np.int64,copy=False)
    if f.min()<0 or f.max()>=count or np.any(np.diff(np.sort(f,axis=1),axis=1)==0):
        raise ValueError('Invalid or repeated connectivity index')
    return f


def triangle_geometry(points,faces):
    p=_points(points);f=_connectivity(faces,3,len(p));v=p[f]
    vectors=np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0])/2
    areas=np.linalg.norm(vectors,axis=1)
    if not np.isfinite(vectors).all() or np.any(areas<=0):
        raise ValueError('Degenerate or overflowing boundary triangle')
    return v.mean(axis=1),vectors,areas


def surface_integrals(faces,area_vectors,velocity,pressure=None):
    """Exact P1 face flow/pressure integral and area-weighted normal-velocity RMS.

    Signed flow uses outward area vectors. Sum of absolute face-mean flow is
    explicitly a facewise diagnostic, not the exact integral of absolute flux.
    """
    u=_points(velocity);f=_connectivity(faces,3,len(u));av=np.asarray(area_vectors,float)
    if av.shape!=(len(f),3) or not np.isfinite(av).all():raise ValueError('Invalid face area vectors')
    areas=np.linalg.norm(av,axis=1)
    if (areas<=0).any() or not np.isfinite(areas).all():raise ValueError('Invalid boundary area')
    uv=u[f];flux=np.einsum('fi,fi->f',uv.mean(axis=1),av)
    normal=np.einsum('fvi,fi->fv',uv,av/areas[:,None])
    square_mean=(normal.sum(axis=1)**2+(normal**2).sum(axis=1))/12
    total_area=float(areas.sum())
    result={'area':total_area,'outward_flux':float(flux.sum()),
            'sum_absolute_face_mean_flux':float(np.abs(flux).sum()),
            'normal_velocity_area_rms':float(np.sqrt(np.sum(areas*square_mean)/total_area)),
            'maximum_nodal_speed':float(np.linalg.norm(uv,axis=2).max())}
    if pressure is not None:
        pressure=np.asarray(pressure,float)
        if pressure.shape!=(len(u),) or not np.isfinite(pressure).all():raise ValueError('Invalid nodal pressure')
        result['mean_pressure']=float(np.sum(pressure[f].mean(axis=1)*areas)/total_area)
    if not all(np.isfinite(v) for v in result.values()):raise ValueError('Nonfinite boundary integral')
    return result


class TetrahedralAuditMesh:
    """Cache outward boundary geometry and independent P1 volume derivatives."""
    def __init__(self,points,tetrahedra):
        self.points=_points(points).copy();self.tets=_connectivity(tetrahedra,4,len(self.points)).copy()
        p=self.points[self.tets];edges=p[:,1:]-p[:,:1]
        determinant=np.linalg.det(edges)
        scale=np.max(np.linalg.norm(edges,axis=2),axis=1)**3
        degenerate=(np.abs(determinant)<=128*np.finfo(float).eps*scale)|~np.isfinite(determinant)
        if degenerate.any():raise ValueError(f'{int(degenerate.sum())} degenerate or overflowing tetrahedra')
        self.volumes=np.abs(determinant)/6
        self.gradients=np.empty((len(self.tets),4,3))
        self.gradients[:,1:]=np.linalg.inv(edges).transpose(0,2,1)
        self.gradients[:,0]=-self.gradients[:,1:].sum(axis=1)
        template=np.array([[1,2,3],[0,3,2],[0,1,3],[0,2,1]])
        faces=self.tets[:,template].reshape(-1,3)
        parent=np.repeat(np.arange(len(self.tets)),4);opposite=np.tile(np.arange(4),len(self.tets))
        centers,vectors,areas=triangle_geometry(self.points,faces)
        interior_vector=self.points[self.tets[parent,opposite]]-centers
        reverse=np.einsum('ij,ij->i',vectors,interior_vector)>0
        faces[reverse]=faces[reverse][:,[0,2,1]];vectors[reverse]*=-1
        _,first,inverse,counts=np.unique(np.sort(faces,axis=1),axis=0,return_index=True,return_inverse=True,return_counts=True)
        if (counts>2).any():raise ValueError('Nonmanifold tetrahedral face shared by more than two cells')
        # Adjacent cells must lie on opposite sides, including when input winding differs.
        summed=np.zeros((len(first),3));np.add.at(summed,inverse,vectors)
        internal=counts==2
        if np.any(np.linalg.norm(summed[internal],axis=1)>1e-10*areas[first[internal]]):
            raise ValueError('Overlapping or inconsistently connected tetrahedra')
        boundary=first[counts==1]
        if not len(boundary):raise ValueError('Tetrahedral domain has no exterior boundary')
        self.boundary_faces=faces[boundary];self.boundary_parent=parent[boundary];self.boundary_opposite=opposite[boundary]
        self.face_centers=centers[boundary];self.area_vectors=vectors[boundary]
        self.face_areas=areas[boundary]
        self._face_lookup={tuple(sorted(face)):i for i,face in enumerate(self.boundary_faces)}
        self._tree=cKDTree(self.points)
        self.metadata={'nodes':len(self.points),'tetrahedra':len(self.tets),'boundary_triangles':len(boundary),
                       'interior_triangles':int(internal.sum()),'degenerate_tetrahedra':0,
                       'negative_input_winding_tetrahedra':int((determinant<0).sum()),
                       'minimum_tetrahedron_volume':float(self.volumes.min()),'total_volume':float(self.volumes.sum()),
                       'boundary_area':float(self.face_areas.sum()),
                       'coordinate_bounds':{'min':self.points.min(0).tolist(),'max':self.points.max(0).tolist()}}

    def match_surface(self,points,faces,global_node_ids=None,tolerance=None):
        """Require source-point correspondence AND exact unordered mesh triangles.

        No nearest-triangle substitution or cross-specimen registration occurs.
        Returned partial matches must not be treated as a complete cap integral.
        """
        p=_points(points);f=_connectivity(faces,3,len(p))
        if tolerance is None:tolerance=float(np.linalg.norm(np.ptp(self.points,axis=0))*5e-7)
        if not np.isfinite(tolerance) or tolerance<0:raise ValueError('Invalid source-point matching tolerance')
        distances,nearest=self._tree.query(p,k=1);method='nearest original point with distance and uniqueness checks'
        if global_node_ids is not None:
            ids=np.asarray(global_node_ids)
            if ids.shape!=(len(p),) or not np.isfinite(ids).all() or not np.equal(ids,np.round(ids)).all() or (ids<1).any() or (ids>len(self.points)).any():
                raise ValueError('Invalid one-based GlobalNodeID')
            nearest=ids.astype(np.int64)-1;distances=np.linalg.norm(self.points[nearest]-p,axis=1)
            method='one-based original GlobalNodeID, verified against original coordinates'
        valid=distances<=tolerance
        unique=len(np.unique(nearest))==len(nearest)
        matched=[];unmatched=0
        for face in f:
            index=self._face_lookup.get(tuple(sorted(nearest[face]))) if valid[face].all() else None
            if index is None:unmatched+=1
            else:matched.append(index)
        duplicate=len(matched)-len(set(matched))
        complete=bool(valid.all() and unique and unmatched==0 and duplicate==0)
        return {'complete':complete,'method':method,'input_points':len(p),'input_faces':len(f),
                'matched_faces':len(matched),'unmatched_faces':unmatched,'duplicate_faces':duplicate,
                'point_mapping_unique':unique,'unmatched_points':int((~valid).sum()),
                'maximum_point_distance':float(distances.max()),'point_tolerance':tolerance,
                'boundary_face_indices':matched}

    def integrate(self,velocity,pressure=None,face_indices=None):
        idx=np.arange(len(self.boundary_faces)) if face_indices is None else np.asarray(face_indices)
        if idx.ndim!=1 or not len(idx) or not np.issubdtype(idx.dtype,np.integer) or idx.min()<0 or idx.max()>=len(self.boundary_faces) or len(np.unique(idx))!=len(idx):
            raise ValueError('Invalid or duplicate boundary face selection')
        return surface_integrals(self.boundary_faces[idx],self.area_vectors[idx],velocity,pressure)

    def audit(self,velocity,pressure=None):
        u=_points(velocity)
        if u.shape!=self.points.shape:raise ValueError('Velocity must match original mesh nodes')
        result=self.integrate(u,pressure)
        divergence=np.einsum('tni,tni->t',u[self.tets],self.gradients)
        integral=float(divergence@self.volumes)
        absolute=float(np.abs(divergence)@self.volumes)
        residual=result['outward_flux']-integral
        scale=max(absolute,result['sum_absolute_face_mean_flux'],np.finfo(float).tiny)
        result.update(volume_divergence_integral=integral,
                      absolute_volume_divergence_integral=absolute,
                      divergence_volume_rms=float(np.sqrt((divergence**2)@self.volumes/self.volumes.sum())),
                      divergence_theorem_residual=residual,divergence_theorem_relative_residual=abs(residual)/scale)
        if not all(np.isfinite(v) for v in result.values()):raise ValueError('Nonfinite volume audit')
        return result
