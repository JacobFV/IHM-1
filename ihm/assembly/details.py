"""Material-attached skin detail and passive microcirculation primitives (SI)."""
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve
from scipy.sparse.csgraph import connected_components


def physical_mesh(vertices, faces, entity_ids=None):
    """Weld exact coordinates only; optional per-vertex entity IDs prevent cross-entity welds."""
    import trimesh
    v=np.asarray(vertices,float); f=np.asarray(faces,int)
    if not np.isfinite(v).all(): raise ValueError('Nonfinite vertices')
    if entity_ids is None: _,first,inverse=np.unique(v,axis=0,return_index=True,return_inverse=True)
    else:
        ids=np.asarray(entity_ids)
        if ids.shape != (len(v),): raise ValueError('entity_ids must be per vertex')
        _,codes=np.unique(ids,return_inverse=True)
        _,first,inverse=np.unique(np.column_stack((codes,v)),axis=0,return_index=True,return_inverse=True)
    return trimesh.Trimesh(v[first],inverse[f],process=False)


def solve_network(positions, edges, radius_m, boundary_pressures_pa, viscosity_pa_s=None):
    """Sparse steady Dirichlet hydraulic solve. Default viscosity is a coefficient PRIOR.

    Diameter-dependent apparent viscosity priors are intentionally not an audited
    Pries law: 1.8 mPa s above 30 um, 2.4 at 12–30 um, 4.0 below 12 um.
    No claim of red-cell, phase separation, compliance, or native circulation coupling.
    """
    p=np.asarray(positions,float); e=np.asarray(edges,int); r=np.asarray(radius_m,float)
    if e.ndim!=2 or e.shape[1]!=2 or r.shape!=(len(e),) or np.any(r<=0): raise ValueError('Invalid edges/radii')
    if not np.isfinite(p).all() or not np.isfinite(r).all(): raise ValueError('Nonfinite geometry')
    length=np.linalg.norm(p[e[:,1]]-p[e[:,0]],axis=1)
    if np.any(length<=0): raise ValueError('Zero length vessel')
    mu=np.where(2*r<12e-6,.004,np.where(2*r<30e-6,.0024,.0018)) if viscosity_pa_s is None else np.broadcast_to(viscosity_pa_s,r.shape).copy()
    if np.any(mu<=0) or not np.isfinite(mu).all(): raise ValueError('Invalid viscosity')
    resistance=8*mu*length/(np.pi*r**4);conductance=1/resistance
    incidence=coo_matrix((np.tile([1.,-1.],len(e)),(e.ravel(),np.repeat(np.arange(len(e)),2))),shape=(len(p),len(e))).tocsr()
    lap=incidence@diags(conductance)@incidence.T
    fixed=np.array(sorted(boundary_pressures_pa),int)
    if not len(fixed) or np.any(fixed<0) or np.any(fixed>=len(p)): raise ValueError('Invalid boundaries')
    nc,labels=connected_components(lap,directed=False)
    if len(np.unique(labels[fixed]))!=nc: raise ValueError('Every component needs pressure boundary')
    pressure=np.zeros(len(p));pressure[fixed]=[boundary_pressures_pa[int(i)] for i in fixed]
    if not np.isfinite(pressure).all(): raise ValueError('Nonfinite pressure')
    internal=np.setdiff1d(np.arange(len(p)),fixed)
    pressure[internal]=spsolve(lap[internal][:,internal],-lap[internal][:,fixed]@pressure[fixed])
    q=conductance*(incidence.T@pressure); net=incidence@q
    return dict(pressure_pa=pressure,flow_m3_s=q,resistance_pa_s_m3=resistance,viscosity_pa_s=mu,
                maximum_internal_residual_m3_s=float(np.max(np.abs(net[internal]),initial=0)),
                boundary_flow_m3_s=net[fixed],boundary_nodes=fixed,
                dissipation_w=float(np.sum(q*q*resistance)))


def microvascular_unit(center, extent_m, capillaries, seed=0):
    """Paired binary supply/return trees linked by exactly N capillaries; inferred geometry."""
    if capillaries<1 or extent_m<=0: raise ValueError('Positive extent and count required')
    rng=np.random.default_rng(seed); center=np.asarray(center,float)
    p=[[-.7*extent_m,0,0],[.7*extent_m,0,0]];edges=[];r=[];kind=[]
    leaves=np.column_stack((np.zeros(capillaries),rng.uniform(-.45,.45,capillaries)*extent_m,rng.uniform(-.2,.2,capillaries)*extent_m))
    def branch(indices, supply, ret, depth=0):
        if len(indices)==1:
            edges.append((supply,ret));r.append(4e-6);kind.append(2);return
        groups=np.array_split(indices,2)
        for group in groups:
            loc=leaves[group].mean(axis=0);x=.5*extent_m/(depth+1)
            a=len(p);p.append([loc[0]-x,loc[1],loc[2]])
            b=len(p);p.append([loc[0]+x,loc[1],loc[2]])
            radius=8e-6*len(group)**(1/3)
            edges.extend([(supply,a),(b,ret)]);r.extend([radius,radius*1.2]);kind.extend([0,1]);branch(group,a,b,depth+1)
    branch(np.arange(capillaries),0,1)
    return np.asarray(p)+center,np.asarray(edges,int),np.asarray(r),np.asarray(kind)


def sample_hair(vertices, faces, density_per_cm2, seed=0):
    """Poisson follicles, area-uniform on each source triangle. No display thinning here."""
    v=np.asarray(vertices,float);f=np.asarray(faces,int);tri=v[f]
    density=np.broadcast_to(density_per_cm2,(len(f),))
    if np.any(density<0) or not np.isfinite(density).all(): raise ValueError('Invalid density')
    cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(cross,axis=1)/2
    rng=np.random.default_rng(seed);counts=rng.poisson(area*density*1e4)
    face=np.repeat(np.arange(len(f)),counts);uv=rng.random((len(face),2));a=np.sqrt(uv[:,0]);bary=np.column_stack((1-a,a*(1-uv[:,1]),a*uv[:,1]))
    roots=np.einsum('ni,nij->nj',bary,tri[face]);normals=cross[face]/np.linalg.norm(cross[face],axis=1)[:,None]
    return dict(ids=np.arange(len(face),dtype=np.int64),face_index=face,barycentric=bary,roots_m=roots,
                normals=normals,radius_m=rng.uniform(8e-6,16e-6,len(face)),length_m=rng.uniform(.0005,.002,len(face)))


def shaft_mesh(roots, directions, lengths, radii, sides=4):
    """Physical radius polygonal open shaft tubes; no visual radius amplification."""
    roots=np.asarray(roots);n=np.asarray(directions);n=n/np.linalg.norm(n,axis=1)[:,None]
    ref=np.tile([0.,1.,0.],(len(n),1));ref[np.abs(n[:,1])>.9]=[1.,0.,0.]
    a=np.cross(n,ref);a/=np.linalg.norm(a,axis=1)[:,None];b=np.cross(n,a)
    angles=np.arange(sides)*2*np.pi/sides;offset=np.asarray(radii)[:,None,None]*(np.cos(angles)[None,:,None]*a[:,None,:]+np.sin(angles)[None,:,None]*b[:,None,:])
    vertices=np.concatenate((roots[:,None,:]+offset,(roots+n*np.asarray(lengths)[:,None])[:,None,:]+offset),axis=1)
    face=[]
    for j in range(sides):face.extend([[j,(j+1)%sides,j+sides],[(j+1)%sides,(j+1)%sides+sides,j+sides]])
    indices=(np.asarray(face)[None,:,:]+(np.arange(len(n))*2*sides)[:,None,None]).reshape(-1,3)
    return dict(positions=vertices.reshape(-1).tolist(),indices=indices.reshape(-1).tolist())
