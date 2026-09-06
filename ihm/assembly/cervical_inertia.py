"""Geometry-conditioned inertial priors; no native model or geometry mutation.

Convex envelopes are explicit approximations, never repaired source volumes.
All body ledgers use mass, COM, and a COM inertia tensor in one common frame.
"""
from collections import Counter
import numpy as np
from scipy.spatial import ConvexHull,cKDTree


def _topology(faces):
    directed=Counter((int(a),int(b)) for f in faces for a,b in zip(f,np.roll(f,-1)))
    undirected=Counter()
    for (a,b),count in directed.items():undirected[tuple(sorted((a,b)))]+=count
    return {'faces':len(faces),'unique_edges':len(undirected),
            'boundary_edges':sum(n==1 for n in undirected.values()),
            'nonmanifold_edges':sum(n>2 for n in undirected.values()),
            'closed_oriented_edge_manifold':bool(len(faces) and all(n==2 for n in undirected.values())
                and all(directed[(a,b)]==directed[(b,a)] for a,b in directed)),
            'scope':'Indexed edge incidence/orientation only; does not prove source absence of self-intersection'}


def _physical(inertia):
    value=np.asarray(inertia,dtype=float)
    if value.shape!=(3,3) or not np.isfinite(value).all() or not np.allclose(value,value.T,atol=1e-12):
        raise ValueError('Invalid physical inertia tensor')
    second=np.trace(value)/2*np.eye(3)-value
    eigen=np.linalg.eigvalsh(second)
    if np.linalg.eigvalsh(value).min()<=0 or eigen.min()<-1e-12:
        raise ValueError('Impossible physical inertia; no clipping or native fallback allowed')
    return eigen


def _body(body):
    mass=float(body['mass_kg']);center=np.asarray(body['center_m'],dtype=float);inertia=np.asarray(body['inertia_kg_m2'],dtype=float)
    if not np.isfinite(mass) or mass<=0 or center.shape!=(3,) or not np.isfinite(center).all():raise ValueError('Positive finite mass and center required')
    _physical(inertia)
    return mass,center,inertia


def convex_geometry_prior(vertices,faces):
    vertices=np.asarray(vertices,dtype=float);faces=np.asarray(faces,dtype=int).reshape(-1,3)
    if vertices.ndim!=2 or vertices.shape[1]!=3 or len(vertices)<4 or len(vertices)>250000 or not np.isfinite(vertices).all():raise ValueError('Bounded finite 3D vertices required')
    if faces.size and (faces.min()<0 or faces.max()>=len(vertices)):raise ValueError('Invalid surface indices')
    # Translate before QHull/integration to avoid cancellation at large origins.
    origin=vertices.mean(axis=0);local=vertices-origin
    hull=ConvexHull(local);triangles=hull.simplices.copy()
    a,b,c=np.moveaxis(local[triangles],1,0)
    wrong=np.einsum('ij,ij->i',np.cross(b-a,c-a),hull.equations[:,:3])<0
    triangles[wrong]=triangles[wrong][:,[0,2,1]]
    a,b,c=np.moveaxis(local[triangles],1,0)
    volumes=np.einsum('ij,ij->i',a,np.cross(b,c))/6
    if (volumes<=0).any():raise ValueError('Convex hull does not form outward positive tetrahedra')
    volume=volumes.sum();sums=a+b+c
    first=np.einsum('i,ij->j',volumes,sums)/4
    second=np.einsum('i,ij,ik->jk',volumes,sums,sums)
    for point in (a,b,c):second+=np.einsum('i,ij,ik->jk',volumes,point,point)
    second/=20
    mean=first/volume;covariance=second/volume-np.outer(mean,mean)
    if not np.isclose(volume,hull.volume,rtol=1e-10):raise ValueError('Hull integration disagrees with QHull volume')
    topology=_topology(triangles)
    if not topology['closed_oriented_edge_manifold']:raise ValueError('Hull is not a closed oriented edge manifold')
    def depths(points):
        result=[]
        for start in range(0,len(points),256):
            # Facet normals are normalized by QHull; points lie inside hull.
            result.extend(np.maximum(0,-np.max(points[start:start+256]@hull.equations[:,:3].T+hull.equations[:,3],axis=1)))
        return np.asarray(result)
    vertex_depth=depths(local)
    source_centers=local[faces].mean(axis=1) if len(faces) else np.empty((0,3))
    source_depth=depths(source_centers)
    # Finite samples: distance to nearest source vertex bounds distance to mesh
    # at each hull-facet center; it is not a global Hausdorff certificate.
    reverse_distance=cKDTree(local).query(local[triangles].mean(axis=1))[0]
    lower=vertices.min(axis=0);upper=vertices.max(axis=0);box_cov=np.diag((upper-lower)**2/12)
    return {'basis':'Uniform solid convex hull; canal/concavity filling is an explicit geometry prior',
            'volume_m3':float(volume),'centroid_m':(origin+mean).tolist(),'covariance_m2':covariance.tolist(),
            'bounds_min_m':lower.tolist(),'bounds_max_m':upper.tolist(),
            'source_topology':_topology(faces),'hull_topology':topology,
            'geometry_error':{'hull_to_aabb_volume_ratio':float(volume/np.prod(upper-lower)),
                'source_vertex_max_hull_depth_m':float(vertex_depth.max()),
                'source_face_centroid_max_hull_depth_m':float(source_depth.max()) if len(source_depth) else None,
                'hull_face_centroid_max_nearest_source_vertex_m':float(reverse_distance.max()),
                'aabb_vs_hull_covariance_relative_frobenius':float(np.linalg.norm(box_cov-covariance)/np.linalg.norm(covariance)),
                'scope':'Sampled convexification discrepancy and AABB sensitivity; not a certified surface Hausdorff bound or bone-volume accuracy'},
            'source_vertices':len(vertices),'hull_vertices':len(hull.vertices)}


def segment_prior(geometry,mass_kg,*,bone_density_kg_m3=1900.,soft_density_kg_m3=1000.):
    volume=float(geometry['volume_m3']);covariance=np.asarray(geometry['covariance_m2']);mass=float(mass_kg)
    if not all(np.isfinite(x) and x>0 for x in (volume,mass,bone_density_kg_m3,soft_density_kg_m3)):raise ValueError('Positive finite mass/density required')
    bone=volume*bone_density_kg_m3;soft=mass-bone
    if soft<0:raise ValueError('Convex bone proxy exceeds donor segment mass; reidentify the model')
    scale=(1+soft/(soft_density_kg_m3*volume))**(1/3)
    # Dilation about the convex-volume centroid: volume scales s^3 and
    # unnormalized second moments scale s^5. The shell excludes the core.
    second=(bone_density_kg_m3*volume+soft_density_kg_m3*volume*(scale**5-1))*covariance
    inertia=np.trace(second)*np.eye(3)-second
    eigen=_physical(inertia);center=np.asarray(geometry['centroid_m'])
    return {'mass_kg':mass,'center_m':center.tolist(),'inertia_kg_m2':inertia.tolist(),
            'second_moment_eigenvalues_kg_m2':eigen.tolist(),
            'basis':'Disjoint bone-proxy core plus homothetic soft shell within this segment; uniform source-informed density priors',
            'components':{'bone_proxy_mass_kg':bone,'soft_proxy_mass_kg':soft,'bone_density_kg_m3':bone_density_kg_m3,
                'soft_density_kg_m3':soft_density_kg_m3,'soft_shell_scale':scale,
                'envelope_bounds_min_m':(center+scale*(np.asarray(geometry['bounds_min_m'])-center)).tolist(),
                'envelope_bounds_max_m':(center+scale*(np.asarray(geometry['bounds_max_m'])-center)).tolist(),
                'scope':'Bone hull fills canals/concavities; shell is a lumped distribution prior, not a measured exclusive tissue compartment. Adjacent segment envelopes may overlap.'}}


def homogeneous_prior(geometry,mass_kg):
    covariance=np.asarray(geometry['covariance_m2']);second=float(mass_kg)*covariance
    body={'mass_kg':float(mass_kg),'center_m':geometry['centroid_m'],
          'inertia_kg_m2':(np.trace(second)*np.eye(3)-second).tolist(),
          'basis':'Uniform composite mass inside convex anatomical envelope; head includes intracranial/soft mass, not bone density',
          'effective_density_kg_m3':float(mass_kg)/geometry['volume_m3']}
    _body(body);return body


def transform_body(body,transform):
    transform=np.asarray(transform,dtype=float)
    if transform.shape!=(4,4) or not np.isfinite(transform).all() or not np.allclose(transform[3],[0,0,0,1]):raise ValueError('A proper rigid transform is required')
    rotation=transform[:3,:3]
    if not np.allclose(rotation.T@rotation,np.eye(3),atol=1e-10) or not np.isclose(np.linalg.det(rotation),1,atol=1e-10):raise ValueError('A proper rigid transform is required')
    mass,center,inertia=_body(body)
    return {**body,'mass_kg':mass,'center_m':(rotation@center+transform[:3,3]).tolist(),'inertia_kg_m2':(rotation@inertia@rotation.T).tolist()}


def _parallel(center):return np.dot(center,center)*np.eye(3)-np.outer(center,center)


def combine_bodies(bodies):
    rows=[_body(b) for b in bodies];mass=sum(m for m,_,_ in rows)
    if not rows:raise ValueError('At least one body required')
    center=sum(m*c for m,c,_ in rows)/mass
    inertia=sum(i+m*_parallel(c-center) for m,c,i in rows)
    return {'mass_kg':mass,'center_m':center.tolist(),'inertia_kg_m2':inertia.tolist()}


def partition_body(parent,children):
    mass,center,inertia=_body(parent);rows=[_body(b) for b in children]
    residual_mass=mass-sum(m for m,_,_ in rows)
    if residual_mass<=0:raise ValueError('No positive residual torso mass')
    residual_center=(mass*center-sum((m*c for m,c,_ in rows),start=np.zeros(3)))/residual_mass
    residual_inertia=inertia+mass*_parallel(center)-sum((i+m*_parallel(c) for m,c,i in rows),start=np.zeros((3,3)))-residual_mass*_parallel(residual_center)
    eigen=_physical(residual_inertia)
    result={'mass_kg':residual_mass,'center_m':residual_center.tolist(),'inertia_kg_m2':residual_inertia.tolist(),
            'second_moment_eigenvalues_kg_m2':eigen.tolist()}
    reconstructed=combine_bodies(children+[result])
    result['reconstruction_residual']={'mass_kg':reconstructed['mass_kg']-mass,
        'first_moment_kg_m':(mass*(np.asarray(reconstructed['center_m'])-center)).tolist(),
        'inertia_frobenius_kg_m2':float(np.linalg.norm(np.asarray(reconstructed['inertia_kg_m2'])-inertia))}
    return result
