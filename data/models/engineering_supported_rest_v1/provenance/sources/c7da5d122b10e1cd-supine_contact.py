"""Provenance-bound posterior skin quadrature and finite-strain foundation law.

No native integration; no additional inertia. The retained skin material values
are generic priors, not measured mattress or patient contact properties.
"""
import numpy as np


def posterior_envelope(vertices, faces, spacing_m=.005, maximum_cells=100000):
    """Intersect source triangles with a +X ray grid and retain the first surface.

Each point is on a retained face. Raster area is reference projected area;
perimeter-cell area error and between-sample extrema are discretization limits.
"""
    vertices=np.asarray(vertices,float);faces=np.asarray(faces,int)
    if vertices.ndim!=2 or vertices.shape[1]!=3 or not np.isfinite(vertices).all():raise ValueError('Finite surface vertices required')
    if faces.ndim!=2 or faces.shape[1]!=3 or len(faces)>300000:raise ValueError('Bounded triangular surface required')
    if not np.isfinite(spacing_m) or not .002<=spacing_m<=.02:raise ValueError('Quadrature spacing outside declared bounds')
    lower=vertices[:,1:].min(0);upper=vertices[:,1:].max(0)
    shape=np.ceil((upper-lower)/spacing_m).astype(int)
    if np.prod(shape)>maximum_cells:raise ValueError('Posterior grid exceeds cell budget')
    depth=np.full(tuple(shape),np.inf);owner=np.full(tuple(shape),-1,dtype=int)
    for index,triangle in enumerate(vertices[faces]):
        yz=triangle[:,1:];a,b=yz[1]-yz[0],yz[2]-yz[0]
        determinant=a[0]*b[1]-a[1]*b[0]
        if abs(determinant)<1e-16:continue
        start=np.maximum(0,np.ceil((yz.min(0)-lower)/spacing_m-.5).astype(int))
        stop=np.minimum(shape-1,np.floor((yz.max(0)-lower)/spacing_m-.5).astype(int))
        if (stop<start).any():continue
        yy,zz=np.meshgrid(np.arange(start[0],stop[0]+1),np.arange(start[1],stop[1]+1),indexing='ij')
        y=lower[0]+(yy+.5)*spacing_m-yz[0,0];z=lower[1]+(zz+.5)*spacing_m-yz[0,1]
        u=(y*b[1]-z*b[0])/determinant;v=(a[0]*z-a[1]*y)/determinant
        x=triangle[0,0]+u*(triangle[1,0]-triangle[0,0])+v*(triangle[2,0]-triangle[0,0])
        inside=(u>=-1e-10)&(v>=-1e-10)&(u+v<=1+1e-10)&(x<depth[yy,zz])
        depth[yy[inside],zz[inside]]=x[inside];owner[yy[inside],zz[inside]]=index
    yy,zz=np.where(owner>=0)
    points=np.column_stack((depth[yy,zz],lower[0]+(yy+.5)*spacing_m,lower[1]+(zz+.5)*spacing_m))
    if not len(points):raise ValueError('No posterior skin intersections')
    return dict(points_source_m=points,face_indices=owner[yy,zz],area_m2=np.full(len(points),spacing_m**2),
                raster_shape=shape,minimum_full_surface_x_m=float(vertices[:,0].min()),
                minimum_sampled_surface_x_m=float(points[:,0].min()))


def foundation(points,velocities,areas,body_indices,body_origins,plane_x_m,material):
    """Normal confined neo-Hookean columns, with dissipative contact friction.

The face geometry is rigidly attached; column energy models local compression.
A rigid stationary plane receives the equal/opposite resultant and moment.
"""
    points=np.asarray(points,float);velocities=np.asarray(velocities,float);areas=np.asarray(areas,float)
    body_indices=np.asarray(body_indices,int);body_origins=np.asarray(body_origins,float)
    if points.shape!=velocities.shape or points.ndim!=2 or points.shape[1]!=3 or areas.shape!=(len(points),):raise ValueError('Contact array shape mismatch')
    if not all(np.isfinite(a).all() for a in (points,velocities,areas,body_origins)) or (areas<=0).any():raise ValueError('Finite positive quadrature required')
    h=float(material['total_layer_thickness_m']);mu=float(material['shear_modulus_pa']);lam=float(material['lame_lambda_pa'])
    minimum=float(material['minimum_thickness_ratio'])
    if h<=0 or mu<=0 or lam<0 or not 0<minimum<1:raise ValueError('Invalid retained material/domain')
    penetration=np.maximum(0.,plane_x_m-points[:,0]);stretch=1-penetration/h
    if (stretch<minimum-1e-12).any():raise ValueError('Skin foundation compression exceeds declared finite-strain domain')
    log=np.log(stretch)
    pressure=-(mu*(stretch-1/stretch)+lam*log/stretch)
    elastic=areas*pressure
    normal=np.maximum(0.,elastic*(1-material['dissipation_s_m']*velocities[:,0]))
    tangent=velocities[:,1:];speed=np.linalg.norm(tangent,axis=1)
    friction_coefficient=material['dynamic_friction']*np.tanh(speed/material['transition_velocity_m_s'])+material['viscous_friction']*speed
    forces=np.zeros_like(points);forces[:,0]=normal
    forces[:,1:]=-normal[:,None]*friction_coefficient[:,None]*tangent/np.maximum(speed[:,None],1e-30)
    stored=areas*h*(.5*mu*(stretch**2-1)-mu*log+.5*lam*log**2)
    body_forces=np.zeros_like(body_origins);body_moments=np.zeros_like(body_origins)
    np.add.at(body_forces,body_indices,forces)
    np.add.at(body_moments,body_indices,np.cross(points-body_origins[body_indices],forces))
    bed_force=-forces.sum(0);bed_moment=-np.cross(points,forces).sum(0)
    dissipative=forces.copy();dissipative[:,0]-=elastic
    dissipation_power=float(np.sum(dissipative*velocities))
    if dissipation_power>1e-8:raise ValueError('Contact dissipative law created energy')
    return dict(point_forces_n=forces,body_forces_n=body_forces,body_moments_nm=body_moments,
                bed_force_n=bed_force,bed_moment_about_source_origin_nm=bed_moment,
                force_balance_residual_n=body_forces.sum(0)+bed_force,
                moment_balance_residual_nm=(body_moments+np.cross(body_origins,body_forces)).sum(0)+bed_moment,
                elastic_energy_j=float(stored.sum()),power_to_body_w=float(np.sum(forces*velocities)),
                dissipative_power_w=dissipation_power,maximum_penetration_m=float(penetration.max(initial=0)),
                contacting_points=int(np.sum(penetration>0)))
