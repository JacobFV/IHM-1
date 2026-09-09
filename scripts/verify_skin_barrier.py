"""Passive sampled skin response; position repair is tracked separately."""
import unittest
from types import SimpleNamespace
import numpy as np
from scipy.spatial import cKDTree
from ihm.assembly.environment_dynamics import project_skin_barrier


def mesh(x, v=None):
    x = np.array(x, float)
    return SimpleNamespace(x=x, rest=x.copy(), v=np.zeros_like(x) if v is None else np.array(v, float), mass=1.5/437, fixed=np.zeros(len(x), bool))


class Tests(unittest.TestCase):
    def test_moving_boundary_does_not_amplify_position_error(self):
        cloth = mesh([[.048,0,0]])
        points = np.array([[.010,0,0]])
        impulse = project_skin_barrier(cloth,cKDTree(points),points,.001,[[.5,0,0]])
        np.testing.assert_allclose(cloth.v, [[.5,0,0]])
        np.testing.assert_allclose(cloth.x, [[.058,0,0]])
        np.testing.assert_allclose(impulse.sum(axis=0)+cloth.mass*cloth.v.sum(axis=0),0,atol=1e-14)
        self.assertGreater(cloth.skin_barrier_state['normal_impulse_dissipation_j'],0)
        self.assertAlmostEqual(cloth.skin_barrier_state['max_projection_m'],.01)

    def test_stationary_projection_cannot_create_velocity(self):
        cloth = mesh([[.038,0,0]])
        points = np.zeros((1,3))
        impulse = project_skin_barrier(cloth,cKDTree(points),points,.001)
        np.testing.assert_array_equal(cloth.v,0)
        np.testing.assert_array_equal(impulse,0)
        self.assertEqual(cloth.skin_barrier_state['kinetic_energy_change_j'],0)

    def test_stationary_collision_dissipates_only_normal_motion(self):
        cloth = mesh([[.038,0,0]], [[-2,3,0]])
        points = np.zeros((1,3))
        initial = cloth.mass*cloth.v.sum(axis=0)
        impulse = project_skin_barrier(cloth,cKDTree(points),points,.001)
        np.testing.assert_allclose(cloth.v,[[0,3,0]])
        np.testing.assert_allclose(cloth.mass*cloth.v.sum(axis=0)+impulse.sum(axis=0),initial)
        self.assertLess(cloth.skin_barrier_state['kinetic_energy_change_j'],0)
        self.assertAlmostEqual(cloth.skin_barrier_state['normal_impulse_dissipation_j'],2*cloth.mass)

    def test_separating_and_fixed_particles_are_unchanged(self):
        cloth = mesh([[.038,0,0],[-.038,0,0]], [[2,0,0],[0,0,0]])
        cloth.fixed[1] = True
        points = np.zeros((1,3))
        impulse = project_skin_barrier(cloth,cKDTree(points),points,.001)
        np.testing.assert_allclose(cloth.v,[[2,0,0],[0,0,0]])
        np.testing.assert_allclose(cloth.x[1],[-.038,0,0])
        np.testing.assert_array_equal(impulse,0)

    def test_spring_projection_energy_is_reported(self):
        cloth = mesh([[.038,0,0],[.2,0,0]])
        cloth.a=np.array([0]);cloth.b=np.array([1]);cloth.length=np.array([.162]);cloth.stiffness=np.array([45.])
        points = np.zeros((1,3))
        project_skin_barrier(cloth,cKDTree(points),points,.001)
        self.assertAlmostEqual(cloth.skin_barrier_state['spring_energy_change_j'],.5*45*.01**2)
        self.assertEqual(cloth.skin_barrier_state['kinetic_energy_change_j'],0)

    def test_coincident_vertex_remains_finite(self):
        cloth=mesh([[0,0,0]])
        points=np.zeros((1,3))
        project_skin_barrier(cloth,cKDTree(points),points,.001)
        self.assertAlmostEqual(np.linalg.norm(cloth.x),.048)
        np.testing.assert_array_equal(cloth.v,0)


if __name__=='__main__':unittest.main()
