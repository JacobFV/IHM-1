"""Small analytic work/geodesic fixtures; no native or anatomical geometry load."""
import json
import numpy as np
from audit_arm26_wrap_geometry import ROOT,audit,cylinder_branch,plane_section_geodesic_curvature
r=audit();assert r==json.loads((ROOT/'data/research/arm26_wrap_geometry/audit.json').read_text())
assert len(r['paths'])==4 and all(p['wrap_count']==1 for p in r['paths'])
errors=[]
for sense in (-1,1):
    for dz in (-.1,0,.08):
        points=np.array([[.04,.004,-.03],[-.035,.006,-.03+dz]])
        result=cylinder_branch(*points,.012,sense)
        contacts=result['contacts'];gradient=result['endpoint_gradients']
        # Independent reconstructed two straight segments + helical arc.
        np.testing.assert_allclose(result['length'],np.linalg.norm(points[0]-contacts[0])+result['wrap_length']+np.linalg.norm(points[1]-contacts[1]),atol=1e-15)
        for i in range(2):
            for j in range(3):
                step=np.zeros((2,3));step[i,j]=1e-7
                fd=(cylinder_branch(*(points+step),.012,sense)['length']-cylinder_branch(*(points-step),.012,sense)['length'])/2e-7
                errors.append(abs(fd-gradient[i,j]))
        # Axial endpoint stationarity: no residual virtual work in wrap sliding.
        helix_slope=(contacts[1,2]-contacts[0,2])/result['wrap_length']
        assert abs(-gradient[0,2]-helix_slope)<1e-14
        assert abs(-gradient[1,2]+helix_slope)<1e-14
assert max(errors)<1e-8
axes=np.array([.027559229367297551,.02204738349383804,.02204738349383804]);oblique=np.array([[1,-1,0],[1,1,-2.]])
oblique/=np.linalg.norm(oblique,axis=1)[:,None]
generic=plane_section_geodesic_curvature(axes,oblique,.7)
principal=plane_section_geodesic_curvature(axes,np.eye(3)[:2],.7)
assert generic>1 and principal<1e-12
print(json.dumps({'passed':True,'max_exact_cylinder_endpoint_gradient_error':max(errors),'generic_plane_section_geodesic_curvature_per_m':generic,'principal_plane_section_geodesic_curvature_per_m':principal,'scope':'Synthetic geometric identities only; no corrected native model acceptance'}))
