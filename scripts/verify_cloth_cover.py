"""Geometric and causal checks for explicitly one-sided supine cover contact."""
import unittest
from types import SimpleNamespace
import numpy as np
from ihm.assembly.cloth_cover import SupineClothCover


def mesh(x, v=None):
    x=np.asarray(x,float)
    return SimpleNamespace(x=x, v=np.zeros_like(x) if v is None else np.asarray(v,float),
        mass=.003, fixed=np.zeros(len(x),bool))


class Tests(unittest.TestCase):
    def test_upper_envelope_owner_and_no_coverage(self):
        cover=SupineClothCover([[0,0,0],[0,0,.02]],radius_m=.05)
        depth,normal,owner=cover.query([[0,0,0],[1,0,0]])
        np.testing.assert_allclose(depth,[.07,0]);np.testing.assert_array_equal(owner,[1,-1])
        np.testing.assert_allclose(normal,[[0,0,1],[0,0,1]])
        tie=SupineClothCover([[0,0,0],[0,0,0]],radius_m=.05)
        self.assertEqual(tie.query([[0,0,0]])[2][0],0)

    def test_normal_matches_height_gradient(self):
        cover=SupineClothCover([[0,0,0]],radius_m=.05)
        point=np.array([[.02,.01,-1.]])
        depth,normal,_=cover.query(point)
        gradient=[]
        for axis in (0,1):
            plus=point.copy();minus=point.copy();plus[0,axis]+=1e-7;minus[0,axis]-=1e-7
            gradient.append((cover.query(plus)[0][0]-cover.query(minus)[0][0])/2e-7)
        expected=np.array([-gradient[0],-gradient[1],1.]);expected/=np.linalg.norm(expected)
        np.testing.assert_allclose(normal[0],expected,atol=1e-8)

    def test_moving_boundary_velocity_impulse_work_and_reciprocal_momentum(self):
        cover=SupineClothCover([[0,0,.01]],[[0,0,.5]])
        cloth=mesh([[0,0,.048]])
        impulse=cover.project(cloth,.001)
        np.testing.assert_allclose(cloth.x,[[0,0,.058]])
        np.testing.assert_allclose(cloth.v,[[0,0,.5]])
        np.testing.assert_allclose(cloth.mass*cloth.v.sum(axis=0)+impulse.sum(axis=0),0)
        state=cloth.skin_barrier_state
        self.assertAlmostEqual(state['contact_impulse_dissipation_j'],.5*cloth.mass*.5**2)
        self.assertAlmostEqual(state['gravity_projection_energy_change_j'],cloth.mass*.01*9.81)

    def test_stationary_oblique_surface_only_removes_normal_approach(self):
        cover=SupineClothCover([[0,0,0]],radius_m=.05)
        cloth=mesh([[.03,0,.01]],[[-2,3,-1]])
        before=cloth.v.copy();_,normal,_=cover.query(cloth.x)
        impulse=cover.project(cloth,.001,friction=0)
        self.assertAlmostEqual(float(cloth.v[0]@normal[0]),0)
        np.testing.assert_allclose(cloth.v-before, np.maximum(0,-before@normal[0])[:,None]*normal)
        np.testing.assert_allclose(cloth.mass*(cloth.v-before).sum(axis=0)+impulse.sum(axis=0),0,atol=1e-14)
        self.assertLess(cloth.skin_barrier_state['kinetic_energy_change_j'],0)

    def test_stationary_projection_creates_no_velocity_and_fixed_nodes_stay(self):
        cover=SupineClothCover([[0,0,0]])
        cloth=mesh([[0,0,0],[0,0,0]]);cloth.fixed[1]=True
        impulse=cover.project(cloth,.001)
        np.testing.assert_array_equal(cloth.v,0);np.testing.assert_array_equal(impulse,0)
        np.testing.assert_array_equal(cloth.x[1],0)
        self.assertEqual(cloth.skin_barrier_state['vertex_count'],1)

    def test_coulomb_budget_and_moving_boundary_dissipation(self):
        cover=SupineClothCover([[0,0,0]],[[1,0,.5]])
        cloth=mesh([[0,0,.01]],[[3,0,-1.5]])
        before=cloth.v.copy()
        impulse=cover.project(cloth,.001,friction=.45)
        # Relative normal speed -2 gives 2m/s normal correction and .9m/s
        # tangential correction: mu times normal impulse, never overshooting.
        np.testing.assert_allclose(cloth.v,[[2.1,0,.5]])
        np.testing.assert_allclose(cloth.mass*(cloth.v-before).sum(axis=0)+impulse.sum(axis=0),0,atol=1e-14)
        self.assertGreater(cloth.skin_barrier_state['contact_impulse_dissipation_j'],0)
        cloth=mesh([[0,0,.01]],[[1.1,0,-1.5]])
        cover.project(cloth,.001,friction=.45)
        np.testing.assert_allclose(cloth.v,[[1,0,.5]])

    def test_rim_and_empty_topology_are_finite(self):
        cover=SupineClothCover([[0,0,0]],radius_m=.05)
        depth,normal,owner=cover.query([[.05,0,-.01]])
        self.assertTrue(np.isfinite(normal).all());self.assertAlmostEqual(np.linalg.norm(normal),1)
        empty=SupineClothCover(np.empty((0,3)))
        depth,normal,owner=empty.query([[0,0,0]])
        np.testing.assert_array_equal(depth,0);np.testing.assert_array_equal(owner,-1)

    def test_does_not_mutate_cached_source_arrays(self):
        points=np.array([[0.,0,0]]);velocity=np.array([[0.,0,.5]])
        cover=SupineClothCover(points,velocity);points[:]=10;velocity[:]=10
        np.testing.assert_array_equal(cover.points,[[0,0,0]])
        np.testing.assert_array_equal(cover.velocity,[[0,0,.5]])


if __name__=='__main__':unittest.main()
