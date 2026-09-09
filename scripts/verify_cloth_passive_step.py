"""Experimental implicit cloth acceptance, work and atomic rejection tests."""
import unittest
from types import SimpleNamespace
import numpy as np
from ihm.assembly.cloth_passive_step import ClothPassiveStepper,PassiveStepRejected


def mesh(points,edges=(),lengths=None):
 x=np.asarray(points,float);a=np.array([e[0] for e in edges],int);b=np.array([e[1] for e in edges],int)
 return SimpleNamespace(x=x,v=np.zeros_like(x),fixed=np.zeros(len(x),bool),mass=1.,a=a,b=b,
  edges=np.column_stack((a,b,np.ones(len(a)))),length=np.linalg.norm(x[b]-x[a],axis=1) if lengths is None else np.asarray(lengths,float),stiffness=np.ones(len(a))*45)


class Tests(unittest.TestCase):
 def test_free_fall_consistent_kinematics_and_energy(self):
  m=mesh([[0,0,1]]);s=ClothPassiveStepper(m);x=m.x.copy()
  r=s.step(m,.01)
  np.testing.assert_allclose(m.v,(m.x-x)/.01)
  np.testing.assert_allclose(m.v,[[0,0,-.0981]],atol=1e-7)
  self.assertLessEqual(r['energy_excess_j'],1e-8)
  np.testing.assert_allclose(r['momentum_residual_ns'],0,atol=1e-6)

 def test_stationary_contact_dissipates_and_tracks_reaction(self):
  m=mesh([[0,0,0]]);m.v[0,2]=-1
  r=ClothPassiveStepper(m).step(m,.01,planes=[{'node':0,'normal':[0,0,1],'offset_m':0}])
  self.assertGreaterEqual(m.x[0,2],-1e-7);self.assertLessEqual(r['energy_excess_j'],1e-8)
  self.assertLess(r['plane_reactions'][0]['impulse_ns'][2],-1)
  np.testing.assert_allclose(r['momentum_residual_ns'],0,atol=1e-5)

 def test_moving_contact_work_and_consistent_motion(self):
  m=mesh([[0,0,0]])
  r=ClothPassiveStepper(m).step(m,.01,gravity=[0,0,0],planes=[{'node':0,'normal':[0,0,1],'offset_m':0,'velocity_m_s':[0,0,.5]}])
  np.testing.assert_allclose(m.v,[[0,0,.5]],atol=1e-4)
  self.assertGreater(r['prescribed_boundary_work_j'],0)
  self.assertLessEqual(r['energy_excess_j'],1e-8)
  np.testing.assert_allclose(r['momentum_residual_ns'],0,atol=1e-5)

 def test_spring_release_preserves_momentum_and_dissipates(self):
  m=mesh([[0,0,0],[1.1,0,0]],[(0,1)],lengths=[1.])
  s=ClothPassiveStepper(m)
  for _ in range(10):
   r=s.step(m,.005,gravity=[0,0,0]);self.assertLessEqual(r['energy_excess_j'],1e-8)
   np.testing.assert_allclose(m.v.sum(axis=0),0,atol=1e-6)

 def test_initial_infeasibility_rejected_atomically(self):
  m=mesh([[0,0,0],[2,0,0]],[(0,1)],lengths=[1.]);before=m.x.copy()
  with self.assertRaises(PassiveStepRejected):ClothPassiveStepper(m).step(m,.01)
  np.testing.assert_array_equal(m.x,before);np.testing.assert_array_equal(m.v,0)

 def test_nonconvex_stationary_point_with_energy_creation_is_rejected(self):
  # A compressed triangular spring network can converge to an implicit
  # stationary point that creates energy. Solver convergence alone is insufficient.
  x=[[-.004023660657323963,-.09693797494702056,.10325537880331428],
     [.0115875810844058,.14323847701548728,.08567532486474887],
     [-.013048738737548983,.08274959864925857,-.06246175594314021]]
  m=mesh(x,[(0,1),(1,2),(0,2)],lengths=[1.,1.,1.])
  m.v[:]=[[-.8584621878220764,-1.2278016556855142,-4.577373147261053],
           [-3.8503055215453625,1.4179226397338889,-1.2837034379499594],
           [1.0209454804015294,-3.690753783402222,.7383759468367054]]
  before_x=m.x.copy();before_v=m.v.copy()
  with self.assertRaises(PassiveStepRejected) as rejected:
   ClothPassiveStepper(m).step(m,.05039960623421782,gravity=[0,0,0])
  self.assertTrue(rejected.exception.diagnostics['solver_success'])
  self.assertGreater(rejected.exception.diagnostics['energy_excess_j'],.6)
  np.testing.assert_array_equal(m.x,before_x);np.testing.assert_array_equal(m.v,before_v)

 def test_four_overlapping_tangents_do_not_become_four_equalities(self):
  m=mesh([[0,0,1.]]);m.v[0,2]=-20
  normals=np.array([[0,0,1],[.1,0,1],[-.1,.1,1],[0,-.1,1]],float)
  normals/=np.linalg.norm(normals,axis=1)[:,None]
  planes=[{'node':0,'normal':n.tolist(),'offset_m':height} for n,height in zip(normals,[.3,.31,.32,.29])]
  r=ClothPassiveStepper(m).step(m,.05,gravity=[0,0,0],planes=planes)
  self.assertEqual(r['solver_message'],'Sparse active-set Newton converged')
  self.assertLessEqual(r['constraint_violation_m'],1e-7)
  np.testing.assert_allclose(r['momentum_residual_ns'],0,atol=1e-6)

 def test_releasing_contact_has_no_spurious_normal_impulse(self):
  m=mesh([[0,0,0]]);m.v[0,2]=1
  r=ClothPassiveStepper(m).step(m,.01,gravity=[0,0,0],planes=[{'node':0,'normal':[0,0,1],'offset_m':0}])
  np.testing.assert_allclose(m.v,[[0,0,1]],atol=1e-12)
  np.testing.assert_array_equal(r['plane_reactions'][0]['impulse_ns'],0)

 def test_compressed_spring_can_release_passively(self):
  m=mesh([[0,0,0],[.3,0,0]],[(0,1)],lengths=[1.]);m.v[:]=[[0,1,0],[0,-1,0]]
  r=ClothPassiveStepper(m).step(m,.02,gravity=[0,0,0])
  self.assertLess(r['energy_excess_j'],0)
  self.assertGreater(np.linalg.norm(m.x[1]-m.x[0]),.3)
  np.testing.assert_allclose(m.v.sum(axis=0),0,atol=1e-9)

 def test_external_force_work_and_sparse_full_mesh_support(self):
  from ihm.assembly.environment_dynamics import SpringMesh
  m=SpringMesh('test','cloth',[-.42,-.9,0],[.42,-.04,0],1.5)
  s=ClothPassiveStepper(m);planes=[{'node':i,'normal':[0,0,1],'offset_m':0} for i in range(len(m.x))]
  r=s.step(m,.005,planes=planes)
  # A fully supported resting sheet must gain no energy at all. The cloth now
  # also carries a hinge bending term, so this sum is round-off rather than a
  # bit-exact zero; the bound is still twelve orders below the 1e-8 acceptance.
  self.assertLessEqual(abs(r['energy_excess_j']),1e-20)
  self.assertEqual(r['bending']['folds_over_60_deg'],0)
  force=np.zeros_like(m.x);force[-1,2]=.5
  for _ in range(3):
   r=s.step(m,.005,planes=planes,forces_n=force)
   self.assertLessEqual(r['energy_excess_j'],1e-8)
   self.assertLessEqual(r['max_extension_ratio'],1.120001)
   np.testing.assert_allclose(r['momentum_residual_ns'],0,atol=1e-6)
  self.assertGreater(r['external_force_work_j'],0)
  self.assertGreater(m.x[-1,2],0)

 def test_overlapping_contact_and_material_rows_no_longer_break_the_solve(self):
  # Two anchored nodes exactly at the extension ceiling. The six anchor
  # equalities already fix both positions, so the active material-limit row is
  # linearly dependent and the KKT matrix is exactly singular. That used to
  # reject the step and advise a smaller timestep, which cannot repair a rank
  # deficiency. The redundant set is consistent, so the step must be accepted.
  m=mesh([[0,0,0],[1.12,0,0]],[(0,1)],lengths=[1.]);m.fixed[:]=True
  before=m.x.copy()
  r=ClothPassiveStepper(m).step(m,.01,gravity=[0,0,0])
  self.assertEqual(r['solver_message'],'Sparse active-set Newton converged')
  self.assertTrue(any('redundant' in row['note'] for row in r['redundant_constraint_solves']))
  self.assertLessEqual(r['constraint_violation_m'],1e-7)
  self.assertLessEqual(r['energy_excess_j'],1e-8)
  np.testing.assert_allclose(m.x,before,atol=1e-12)
  np.testing.assert_allclose(m.v,0,atol=1e-12)

 def test_duplicated_contact_planes_share_one_total_reaction(self):
  # Sampled skin spheres routinely produce several identical tangents on one
  # vertex. Splitting one non-unique multiplier between them must not change
  # the total reaction the support receives.
  def impulse(count):
   m=mesh([[0,0,1.]]);m.v[0,2]=-4
   planes=[{'node':0,'normal':[0,0,1],'offset_m':.5} for _ in range(count)]
   r=ClothPassiveStepper(m).step(m,.05,gravity=[0,0,0],planes=planes)
   np.testing.assert_allclose(r['momentum_residual_ns'],0,atol=1e-6)
   self.assertLessEqual(r['energy_excess_j'],1e-8)
   return np.sum([p['impulse_ns'] for p in r['plane_reactions']],axis=0),m.x[0,2]
  single,height=impulse(1)
  for count in (2,5):
   total,repeat=impulse(count)
   np.testing.assert_allclose(total,single,atol=1e-9)
   self.assertAlmostEqual(repeat,height,places=12)

 def test_dynamically_infeasible_overlapping_planes_reject_atomically(self):
  # Two prescribed planes that are compatible now and conflict by the end of
  # the step. Tolerating redundancy must not tolerate contradiction.
  m=mesh([[0,0,0.]]);before=m.x.copy()
  planes=[{'node':0,'normal':[0,0,1],'offset_m':0.,'velocity_m_s':[0,0,1.]},
          {'node':0,'normal':[0,0,-1],'offset_m':-.005,'velocity_m_s':[0,0,-1.]}]
  with self.assertRaises(PassiveStepRejected):
   ClothPassiveStepper(m).step(m,.05,gravity=[0,0,0],planes=planes)
  np.testing.assert_array_equal(m.x,before);np.testing.assert_array_equal(m.v,0)

 def test_outward_motion_cannot_break_material_limit(self):
  m=mesh([[0,0,0],[1.11,0,0]],[(0,1)],lengths=[1.]);m.v[:]=[[-2,0,0],[2,0,0]]
  r=ClothPassiveStepper(m).step(m,.01,gravity=[0,0,0])
  self.assertLessEqual(r['max_extension_ratio'],1.1200001)
  self.assertLessEqual(r['energy_excess_j'],1e-8)


if __name__=='__main__':unittest.main()
