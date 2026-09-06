import numpy as np
from analyze_arm26_wrap_contacts import ellipsoid,AXES
basis=np.array([[1.,-1,0],[1,1,-2.]])
basis/=np.linalg.norm(basis,axis=1)[:,None]
t=np.linspace(0,1,20);points=AXES*(np.cos(t)[:,None]*basis[0]+np.sin(t)[:,None]*basis[1])
a,b=points[0],points[-1];ta=AXES*basis[1];tb=AXES*(-np.sin(1)*basis[0]+np.cos(1)*basis[1]);ta/=np.linalg.norm(ta);tb/=np.linalg.norm(tb)
r={'p1':a-.02*ta,'p2':b+.02*tb,'r1':a,'r2':b,'c1_unscaled':points[9],'wrap_points':points,'wrap_length':np.linalg.norm(np.diff(points,axis=0),axis=1).sum()}
x=ellipsoid(r)
assert x['plane_residual_m']<1e-12 and max(x['endpoint_surface_tangency_dot'])<1e-12
assert x['contact_projection_displacement_max_m']<1e-12 and x['quadrature_32_64_difference_m']<1e-12
assert x['projected_arc_minus_native_stored_m']>0 and x['planar_section_geodesic_curvature_max_per_m']>1
print('PASS synthetic oblique ellipsoid diagnostics: exact contacts/plane/tangency, nonzero geodesic curvature, chord deficit, quadrature agreement')
