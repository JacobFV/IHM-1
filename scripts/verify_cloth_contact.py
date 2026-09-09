"""Particle self-contact checks, not a claim of continuous cloth collision."""
import unittest
from types import SimpleNamespace
import numpy as np
from ihm.assembly.cloth_contact import ClothSelfContact
from ihm.assembly.environment_dynamics import SpringMesh


def particles(x, v=None, mass=1., edges=(), fixed=None):
    x = np.array(x, dtype=float)
    return SimpleNamespace(x=x, rest=x.copy(),
        v=np.zeros_like(x) if v is None else np.array(v, dtype=float), mass=mass,
        fixed=np.zeros(len(x), bool) if fixed is None else np.array(fixed, bool),
        a=np.array([e[0] for e in edges], int), b=np.array([e[1] for e in edges], int))


class Tests(unittest.TestCase):
    def test_free_pair_preserves_momentum_and_com_without_kinetic_injection(self):
        mesh = particles([[0, 0, 0], [.005, 0, 0]], [[1, 2, 0], [-2, 3, 0]], mass=np.array([2., 3.]))
        momentum = np.sum(mesh.mass[:, None] * mesh.v, axis=0)
        com = np.sum(mesh.mass[:, None] * mesh.x, axis=0)
        result = ClothSelfContact(mesh).resolve(mesh)
        self.assertEqual(result['pair_count'], 1)
        self.assertLess(result['kinetic_energy_change_j'], 0)
        self.assertAlmostEqual(np.linalg.norm(mesh.x[1]-mesh.x[0]), .012)
        np.testing.assert_allclose(np.sum(mesh.mass[:, None]*mesh.v, axis=0), momentum)
        np.testing.assert_allclose(np.sum(mesh.mass[:, None]*mesh.x, axis=0), com)
        np.testing.assert_array_equal(mesh.v[:, 1], [2., 3.])

    def test_projection_does_not_create_motion(self):
        mesh = particles([[0, 0, 0], [0, 0, 0]])
        result = ClothSelfContact(mesh).resolve(mesh)
        self.assertEqual(result['kinetic_energy_change_j'], 0)
        np.testing.assert_array_equal(mesh.v, 0)
        self.assertAlmostEqual(np.linalg.norm(mesh.x[1]-mesh.x[0]), .012)

    def test_neighbors_and_two_hop_neighbors_excluded(self):
        mesh = particles([[0, 0, 0], [.001, 0, 0], [.002, 0, 0]], edges=[(0, 1), (1, 2)])
        self.assertEqual(ClothSelfContact(mesh).resolve(mesh)['pair_count'], 0)
        np.testing.assert_array_equal(mesh.x, mesh.rest)

    def test_fixed_support_reaction_and_separating_velocity(self):
        mesh = particles([[0, 0, 0], [.005, 0, 0]], [[0, 0, 0], [-2, 0, 0]], fixed=[True, False])
        result = ClothSelfContact(mesh).resolve(mesh)
        np.testing.assert_array_equal(mesh.x[0], mesh.rest[0])
        np.testing.assert_allclose(result['support_impulse_ns'], [-2., 0, 0])
        np.testing.assert_array_equal(mesh.v, 0)
        mesh.x[1, 0] = .005
        mesh.v[1, 0] = 2
        self.assertEqual(ClothSelfContact(mesh).resolve(mesh)['kinetic_energy_change_j'], 0)
        self.assertEqual(mesh.v[1, 0], 2)

    def test_folded_real_mesh_and_deterministic_replay(self):
        mesh = SpringMesh('test', 'cloth', [0, 0, 0], [1.8, 2.2, 0], 1.5)
        # Fold right half onto left with a sub-thickness layer spacing.
        right = mesh.x[:, 0] > .9
        mesh.x[right, 0] = 1.8 - mesh.x[right, 0]
        mesh.x[right, 2] = .004
        mesh.v[right, 2] = -.1
        contact = ClothSelfContact(mesh)
        x, v = mesh.x.copy(), mesh.v.copy()
        result = contact.resolve(mesh)
        self.assertGreater(result['pair_count'], 100)
        self.assertLessEqual(result['kinetic_energy_change_j'], 1e-12)
        expected_x, expected_v = mesh.x.copy(), mesh.v.copy()
        mesh.x[:], mesh.v[:] = x, v
        self.assertEqual(contact.resolve(mesh), result)
        np.testing.assert_array_equal(mesh.x, expected_x)
        np.testing.assert_array_equal(mesh.v, expected_v)


if __name__ == '__main__':
    unittest.main()
