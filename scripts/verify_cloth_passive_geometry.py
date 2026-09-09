"""Experimental nonlinear geometry adapter acceptance and atomicity checks."""
import unittest
from types import SimpleNamespace
import numpy as np
from ihm.assembly.cloth_passive_step import ClothPassiveStepper, PassiveStepRejected
from ihm.assembly.cloth_passive_geometry import ClothPassiveGeometry


def mesh(points):
    x = np.array(points, dtype=float)
    return SimpleNamespace(x=x, v=np.zeros_like(x), fixed=np.zeros(len(x),bool), mass=1.,
        a=np.array([],int), b=np.array([],int), edges=np.empty((0,3)), length=np.array([]), stiffness=np.array([]))


def sphere(center=(0,0,0), radius=1., velocity=(0,0,0), owner=73):
    return dict(center_m=center, radius_m=radius, velocity_m_s=velocity, owner=owner)


class Tests(unittest.TestCase):
    def test_resting_sphere_reaction_retains_skin_owner(self):
        m = mesh([[0,0,1]])
        r = ClothPassiveGeometry(ClothPassiveStepper(m)).step(m,.01,spheres=[sphere()])
        self.assertGreaterEqual(m.x[0,2],1.-1e-7)
        self.assertEqual(r['plane_reactions'][0]['owner'],73)
        self.assertLess(r['plane_reactions'][0]['impulse_ns'][2],0)
        self.assertLessEqual(r['energy_excess_j'],1e-8)
        np.testing.assert_allclose(r['momentum_residual_ns'],0,atol=1e-7)

    def test_moving_sphere_boundary_work(self):
        m = mesh([[0,0,1]])
        r = ClothPassiveGeometry(ClothPassiveStepper(m)).step(m,.01,spheres=[sphere(velocity=[0,0,.5])],gravity=[0,0,0])
        np.testing.assert_allclose(m.x,[[0,0,1.005]],atol=1e-7)
        self.assertGreater(r['prescribed_boundary_work_j'],0)
        self.assertLessEqual(r['energy_excess_j'],1e-8)

    def test_sweep_detects_tunneling_with_feasible_endpoints(self):
        m = mesh([[-2,0,0]]); m.v[0,0] = 400.
        r = ClothPassiveGeometry(ClothPassiveStepper(m)).step(m,.01,spheres=[sphere(radius=.5)],gravity=[0,0,0])
        self.assertLessEqual(m.x[0,0],-.5+1e-7)
        self.assertGreaterEqual(r['geometry']['minimum_swept_sphere_gap_m'],-1e-7)

    def test_all_spheres_not_nearest_only(self):
        m = mesh([[0,0,1]])
        objects = [sphere(center=[.05,0,.9],radius=.05,owner=9), sphere()]
        r = ClothPassiveGeometry(ClothPassiveStepper(m)).step(m,.01,spheres=objects)
        self.assertGreaterEqual(m.x[0,2],1.-1e-7)
        self.assertIn(73,[p['owner'] for p in r['plane_reactions']])

    def test_retry_uses_original_state_not_preceding_trial(self):
        m = mesh([[0,0,1.1]]); initial=m.x.copy()
        r = ClothPassiveGeometry(ClothPassiveStepper(m)).step(m,.1,spheres=[sphere()],gravity=[0,0,0],forces_n=[[0,0,-100]])
        self.assertEqual(r['geometry']['outer_iterations'],2)
        np.testing.assert_allclose(m.v,(m.x-initial)/.1,atol=1e-10)
        self.assertGreaterEqual(m.x[0,2],1.-1e-7)

    def test_outer_limit_rejects_atomically(self):
        m = mesh([[0,0,1.1]]); initial=m.x.copy()
        with self.assertRaises(PassiveStepRejected):
            ClothPassiveGeometry(ClothPassiveStepper(m),max_outer_iterations=1).step(m,.1,spheres=[sphere()],gravity=[0,0,0],forces_n=[[0,0,-100]])
        np.testing.assert_array_equal(m.x,initial); np.testing.assert_array_equal(m.v,0)

    def test_initial_overlap_rejected_atomically(self):
        m=mesh([[0,0,.5]]); initial=m.x.copy()
        with self.assertRaises(PassiveStepRejected):
            ClothPassiveGeometry(ClothPassiveStepper(m)).step(m,.01,spheres=[sphere()])
        np.testing.assert_array_equal(m.x,initial)

    def test_bounded_box_top_and_outside_floor(self):
        m=mesh([[0,0,1],[2,0,.1]])
        r=ClothPassiveGeometry(ClothPassiveStepper(m)).step(m,.01,
            boxes=[dict(min_m=[-1,-1,0],max_m=[1,1,1],owner='mattress')],floor_height_m=0)
        self.assertGreaterEqual(m.x[0,2],1.-1e-7)
        self.assertLess(m.x[1,2],.1)
        self.assertIn('mattress',[p['owner'] for p in r['plane_reactions']])

    def test_box_sweep_both_endpoints_outside(self):
        m=mesh([[-2,0,1.1]]);m.v[0]=[400,0,-100]
        r=ClothPassiveGeometry(ClothPassiveStepper(m)).step(m,.01,gravity=[0,0,0],
            boxes=[dict(min_m=[-1,-1,0],max_m=[1,1,1])])
        self.assertGreaterEqual(m.x[0,2],1.-1e-7)

    def test_box_entry_below_top_rejected_without_side_model(self):
        m=mesh([[-2,0,.5]]);m.v[0,0]=400;initial=m.x.copy();initial_v=m.v.copy()
        with self.assertRaises(PassiveStepRejected):
            ClothPassiveGeometry(ClothPassiveStepper(m)).step(m,.01,gravity=[0,0,0],
                boxes=[dict(min_m=[-1,-1,0],max_m=[1,1,1])])
        np.testing.assert_array_equal(m.x,initial);np.testing.assert_array_equal(m.v,initial_v)

    def test_adaptive_moving_spheres_match_manual_substeps(self):
        class SmallSteps(ClothPassiveStepper):
            def step(self, mesh, dt_s, **kwargs):
                if dt_s > .005000001:
                    raise PassiveStepRejected('fixture requests refinement')
                return super().step(mesh, dt_s, **kwargs)
        m=mesh([[0,0,1]])
        r=ClothPassiveGeometry(SmallSteps(m)).advance(m,.02,min_dt_s=.0025,
            spheres=[sphere(velocity=[0,0,.5])],gravity=[0,0,0])
        manual=mesh([[0,0,1]]); adapter=ClothPassiveGeometry(ClothPassiveStepper(manual))
        results=[]
        for i in range(4):
            results.append(adapter.step(manual,.005,spheres=[sphere(center=[0,0,i*.005*.5],velocity=[0,0,.5])],gravity=[0,0,0]))
        np.testing.assert_allclose(m.x,manual.x,atol=1e-10)
        np.testing.assert_allclose(m.v,manual.v,atol=1e-10)
        self.assertEqual(r['accepted_substeps'],4);self.assertEqual(r['finest_dt_s'],.005)
        self.assertEqual(r['rejected_attempts'],3)
        self.assertTrue(all(p['owner']==73 for p in r['plane_reactions']))
        self.assertAlmostEqual(r['prescribed_boundary_work_j'],sum(a['prescribed_boundary_work_j'] for a in results))
        self.assertAlmostEqual(r['energy_excess_j'],r['energy_after_j']-r['energy_before_j']-r['prescribed_boundary_work_j']-r['external_force_work_j'])
        np.testing.assert_allclose(r['total_support_impulse_ns'],np.sum([p['impulse_ns'] for a in results for p in a['plane_reactions']],axis=0),atol=1e-10)

    def test_adaptive_failure_after_accepted_half_rolls_back_whole_interval(self):
        class LaterFailure(ClothPassiveStepper):
            def step(self, mesh, dt_s, **kwargs):
                if dt_s > .005000001 or mesh.x[0,2] < .99999:
                    raise PassiveStepRejected('fixture rejects later motion')
                return super().step(mesh,dt_s,**kwargs)
        m=mesh([[0,0,1]]); initial=m.x.copy(); original_v=m.v.copy()
        with self.assertRaises(PassiveStepRejected):
            ClothPassiveGeometry(LaterFailure(m)).advance(m,.01,min_dt_s=.00125)
        np.testing.assert_array_equal(m.x,initial);np.testing.assert_array_equal(m.v,original_v)

    def test_adaptive_internal_kinematics_are_not_outer_average(self):
        class SmallSteps(ClothPassiveStepper):
            def step(self,mesh,dt_s,**kwargs):
                if dt_s > .005000001: raise PassiveStepRejected('refine')
                return super().step(mesh,dt_s,**kwargs)
        m=mesh([[0,0,1]]);initial=m.x.copy()
        r=ClothPassiveGeometry(SmallSteps(m)).advance(m,.02,forces_n=[[0,0,1]])
        self.assertEqual(r['accepted_substeps'],4)
        self.assertGreater(np.linalg.norm(m.v-(m.x-initial)/.02),.01)
        self.assertLessEqual(r['energy_excess_j'],1e-8)
        np.testing.assert_allclose(r['momentum_residual_ns'],0,atol=1e-7)

    def test_upright_y_axis_floor(self):
        m=mesh([[0,0,0]])
        ClothPassiveGeometry(ClothPassiveStepper(m),up_axis=1).step(m,.01,gravity=[0,-9.81,0],floor_height_m=0)
        np.testing.assert_allclose(m.x,0,atol=1e-7)


if __name__ == '__main__':
    unittest.main()
