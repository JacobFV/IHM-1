"""Behavioral checks for the independent source CFD numerical audit."""
from pathlib import Path
import json
import numpy as np
from ihm.spatial.vascular import TetrahedralAuditMesh, triangle_geometry, surface_integrals

# A reversed tetra winding must not invert physical outward normals.
p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
t=np.array([[0,1,2,3],[0,2,1,4]])
a=TetrahedralAuditMesh(p,t)
b=TetrahedralAuditMesh(p,t[:,[0,2,1,3]])
v=p@np.diag([2.,3.,-1.])+[4.,-2.,3.]
for mesh in [a,b]:
 assert len(mesh.boundary_faces)==6 and np.isclose(mesh.volumes.sum(),1/3)
 result=mesh.audit(v,np.full(len(p),7.))
 assert np.isclose(result['outward_flux'],4/3,atol=1e-13)
 assert np.isclose(result['volume_divergence_integral'],4/3,atol=1e-13)
 assert abs(result['divergence_theorem_residual'])<1e-13
 assert np.isclose(result['mean_pressure'],7)
 assert np.all(np.einsum('ij,ij->i',mesh.area_vectors,p[mesh.tets[mesh.boundary_parent,mesh.boundary_opposite]]-mesh.face_centers)<0)
assert abs(a.audit(np.tile([1.,2.,3.],(5,1)))['outward_flux'])<1e-13

# Unequal triangle areas must weight scalar pressure, not average faces equally.
q=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[2,0,1],[0,2,1]])
f=np.array([[0,1,2],[3,4,5]])
centers,vectors,areas=triangle_geometry(q,f)
r=surface_integrals(f,vectors,np.tile([0.,0.,2.],(6,1)),np.array([2,2,2,8,8,8]))
assert np.isclose(r['outward_flux'],5) and np.isclose(r['mean_pressure'],6.8)
assert np.isclose(r['normal_velocity_area_rms'],2)

# Degenerate/nonmanifold tetrahedra and corrupt fields must be rejected.
for points,tets in [(p,[[0,1,2,2]]),(p,[[0,1,2,3]]*3),(p,[[0,1,2,99]]),(np.zeros((4,3)),[[0,1,2,3]])]:
 try:TetrahedralAuditMesh(points,tets)
 except ValueError:pass
 else:raise AssertionError('invalid tetrahedral mesh accepted')
try:a.audit(np.full((5,3),np.nan))
except ValueError:pass
else:raise AssertionError('nonfinite field accepted')

# A shuffled cap must match actual triangles; a non-boundary triangle must not.
faces=a.boundary_faces[:2];nodes=np.unique(faces)[::-1];lookup={v:i for i,v in enumerate(nodes)}
local=np.array([[lookup[v] for v in face] for face in faces])
m=a.match_surface(p[nodes],local,global_node_ids=nodes+1,tolerance=1e-12)
assert m['complete'] and m['matched_faces']==2 and m['maximum_point_distance']==0
m=a.match_surface(p,np.array([[0,1,2]]),global_node_ids=np.arange(1,6),tolerance=1e-12)
assert not m['complete'] and m['unmatched_faces']==1
print('verified: affine-tet divergence theorem, outward winding, pressure-area weights, exact normal RMS, corruption rejection, cap connectivity matching')

root=Path(__file__).resolve().parents[1];report=root/'data/derived/vascular/audit.json';series=root/'data/derived/vascular/cap-flow.json'
if report.exists() and series.exists():
 data=json.loads(report.read_text());flow=json.loads(series.read_text())
 assert data['mesh']['nodes']==153082 and data['mesh']['tetrahedra']==738037
 assert data['mesh']['degenerate_tetrahedra']==0
 assert data['conservation']['maximum_divergence_theorem_relative_residual']<1e-10
 assert len(flow['time_s'])==len(flow['step'])>0 and np.all(np.diff(flow['time_s'])>0)
 assert not data['units']['absolute_si_confirmed']
 for cap in flow['caps']:
  if cap['matching_complete']:
   assert len(cap['outward_flow'])==len(flow['time_s']) and np.isfinite(cap['outward_flow']).all()
 print('verified: actual source mesh, decoded source clock, finite cap trajectories and discrete conservation identity')
