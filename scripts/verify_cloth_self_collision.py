"""Triangle-level and continuous self-collision: momentum, energy, penetration.

Every assertion here is a conservation or geometry statement. Motion is never
accepted as evidence: the checks are that the impulse adds no kinetic energy,
carries no net force or torque, and that the step ends with no primitive pair
closer than the declared sheet thickness.
"""
import unittest
from types import SimpleNamespace
import numpy as np
from ihm.assembly.cloth_self_collision import ClothSelfCollision
from ihm.assembly.environment_dynamics import SpringMesh


def sheet(x, v=None, mass=1., fixed=None):
    x = np.asarray(x, float)
    return SimpleNamespace(x=x, v=np.zeros_like(x) if v is None else np.asarray(v, float),
                           mass=mass if np.ndim(mass) else float(mass),
                           fixed=np.zeros(len(x), bool) if fixed is None else np.asarray(fixed, bool))


def ledger(test, mesh, before_v, before_x, result, mass):
    momentum = np.sum(mass[:, None]*(mesh.v-before_v), axis=0)
    np.testing.assert_allclose(momentum, 0, atol=1e-14)
    test.assertLessEqual(result['impulse_torque_residual_nms'], 1e-12)
    test.assertLessEqual(result['kinetic_energy_change_j'], 0.)


class Tests(unittest.TestCase):
    def test_vertex_through_triangle_is_caught_by_the_swept_test(self):
        start = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0], [.25, .25, .4]])
        mesh = sheet(start.copy(), [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, -.8]], mass=np.ones(4))
        solver = ClothSelfCollision([[0, 1, 2]], 4, thickness_m=.05, iterations=8)
        before_v = mesh.v.copy()
        mesh.x[3, 2] = -.4
        result = solver.resolve(mesh, previous_positions=start)
        self.assertEqual(result['swept_contacts'], 1)
        ledger(self, mesh, before_v, start, result, np.ones(4))
        # It must come back out on the side it entered from, not merely stop.
        self.assertGreater(mesh.x[3, 2], mesh.x[:3, 2].mean())
        self.assertGreaterEqual(result['residual_min_distance_m'], .05-1e-9)

    def test_without_the_previous_positions_tunnelling_is_reported_undetected(self):
        start = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0], [.25, .25, .4]])
        mesh = sheet(start.copy(), [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, -.8]], mass=np.ones(4))
        mesh.x[3, 2] = -.4
        result = ClothSelfCollision([[0, 1, 2]], 4, thickness_m=.05).resolve(mesh)
        self.assertFalse(result['continuous_test_applied'])
        self.assertEqual(result['swept_contacts'], 0)
        self.assertLess(mesh.x[3, 2], 0)

    def test_edge_edge_crossing_separates_and_conserves_momentum(self):
        # Two non-adjacent edges of two separate triangles crossing near their
        # middles: no vertex is near another vertex, so a particle guard sees
        # nothing at all here.
        points = np.array([[-1., 0, 0], [1., 0, 0], [0, -2., .3],
                           [0, -1., .002], [0, 1., .002], [0, 0, 2.]])
        mesh = sheet(points.copy(), mass=np.ones(6))
        mesh.v[3] = [0, 0, .5]; mesh.v[0] = [0, 0, -.5]
        solver = ClothSelfCollision([[0, 1, 2], [3, 4, 5]], 6, thickness_m=.02, iterations=20)
        self.assertLess(solver.minimum_distance(mesh.x), .02)
        before_v = mesh.v.copy(); before_x = mesh.x.copy()
        result = solver.resolve(mesh)
        self.assertGreater(result['vertex_triangle_and_edge_edge_contacts'], 0)
        ledger(self, mesh, before_v, before_x, result, np.ones(6))
        self.assertGreaterEqual(result['residual_min_distance_m'], .02-1e-9)

    def test_a_primitive_is_never_tested_against_one_it_shares_a_vertex_with(self):
        # A vertex against a triangle it belongs to, or two edges meeting at a
        # vertex, are always within any thickness. Testing them would jam the
        # sheet solid, so the exclusion is structural rather than tolerance based.
        points = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0], [-1., .2, 0], [-.2, -1., 0]])
        solver = ClothSelfCollision([[0, 1, 2], [0, 3, 4]], 5, thickness_m=2.)
        mesh = sheet(points.copy(), mass=np.ones(5))
        vertex_face, edge_edge = solver._pairs(mesh.x, mesh.x)
        faces = solver.triangles[vertex_face[:, 1]]
        self.assertFalse(np.any(faces == vertex_face[:, 0][:, None]))
        first = solver.edges[edge_edge[:, 0]]; second = solver.edges[edge_edge[:, 1]]
        self.assertFalse(np.any(first[:, :, None] == second[:, None, :]))
        single = ClothSelfCollision([[0, 1, 2]], 3, thickness_m=2.)
        alone = sheet(points[:3].copy(), mass=np.ones(3))
        before = alone.x.copy()
        self.assertEqual(single.minimum_distance(alone.x), float('inf'))
        self.assertEqual(single.resolve(alone)['vertex_triangle_and_edge_edge_contacts'], 0)
        np.testing.assert_array_equal(alone.x, before)

    def test_a_hinge_cannot_fold_below_the_declared_sheet_thickness(self):
        # Two triangles of one hinge folded flat onto each other. Their opposite
        # edges share no vertex, so this pair is real contact: fabric of finite
        # thickness cannot close a crease completely.
        points = np.array([[0., 0, 0], [1., 0, 0], [.5, .5, 0], [.5, .4, .0005]])
        mesh = sheet(points.copy(), mass=np.ones(4))
        mesh.v[3] = [0, 0, -.6]
        solver = ClothSelfCollision([[0, 1, 2], [1, 0, 3]], 4, thickness_m=.02, iterations=30)
        self.assertLess(solver.minimum_distance(mesh.x), .02)
        before_v = mesh.v.copy(); before_x = mesh.x.copy()
        result = solver.resolve(mesh)
        ledger(self, mesh, before_v, before_x, result, np.ones(4))
        self.assertGreaterEqual(solver.minimum_distance(mesh.x), .02-1e-9)

    def test_separating_contact_receives_no_impulse(self):
        points = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0], [.25, .25, .004]])
        mesh = sheet(points.copy(), mass=np.ones(4))
        mesh.v[3] = [0, 0, 3.]
        solver = ClothSelfCollision([[0, 1, 2]], 4, thickness_m=.02, iterations=8)
        before_v = mesh.v.copy()
        result = solver.resolve(mesh)
        self.assertEqual(result['kinetic_energy_change_j'], 0.)
        np.testing.assert_array_equal(mesh.v, before_v)
        self.assertGreater(result['max_overlap_m'], 0.)

    def test_fixed_vertices_act_as_supports_and_receive_the_reaction(self):
        points = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0], [.25, .25, .004]])
        mesh = sheet(points.copy(), mass=np.ones(4), fixed=[True, True, True, False])
        mesh.v[3] = [0, 0, -2.]
        solver = ClothSelfCollision([[0, 1, 2]], 4, thickness_m=.02, iterations=8)
        before = np.sum(mesh.v, axis=0)
        result = solver.resolve(mesh)
        np.testing.assert_array_equal(mesh.x[:3], points[:3])
        np.testing.assert_allclose(np.sum(mesh.v, axis=0)+np.asarray(result['fixed_support_impulse_ns']),
                                   before, atol=1e-14)
        self.assertLess(result['kinetic_energy_change_j'], 0.)

    def test_folded_production_blanket_separates_without_energy_or_momentum_error(self):
        mesh = SpringMesh('blanket', 'cloth', [-.42, -.88, 0], [.42, -.04, 0], 1.5)
        grid = mesh.x.reshape(19, 23, 3)
        grid[10:, :, 0] = 2*grid[10, 0, 0]-grid[10:, :, 0]
        grid[10:, :, 2] += .003
        mesh.v[:] = 0.
        mesh.v.reshape(19, 23, 3)[10:, :, 2] = -.4
        mass = np.full(len(mesh.x), mesh.mass)
        solver = ClothSelfCollision(np.asarray(mesh.indices, int).reshape(-1, 3), len(mesh.x),
                                    thickness_m=.006, iterations=40)
        self.assertLess(solver.minimum_distance(mesh.x), .006)
        before_v = mesh.v.copy(); before_x = mesh.x.copy()
        result = solver.resolve(mesh)
        self.assertGreater(result['vertex_triangle_and_edge_edge_contacts'], 500)
        ledger(self, mesh, before_v, before_x, result, mass)
        self.assertGreaterEqual(result['residual_min_distance_m'], .00599)
        self.assertGreaterEqual(solver.minimum_distance(mesh.x), .00599)

    def test_replay_is_deterministic(self):
        mesh = SpringMesh('blanket', 'cloth', [-.42, -.88, 0], [.42, -.04, 0], 1.5)
        grid = mesh.x.reshape(19, 23, 3)
        grid[10:, :, 0] = 2*grid[10, 0, 0]-grid[10:, :, 0]
        grid[10:, :, 2] += .003
        mesh.v.reshape(19, 23, 3)[10:, :, 2] = -.4
        solver = ClothSelfCollision(np.asarray(mesh.indices, int).reshape(-1, 3), len(mesh.x),
                                    thickness_m=.006, iterations=6)
        x, v = mesh.x.copy(), mesh.v.copy()
        first = solver.resolve(mesh)
        expected_x, expected_v = mesh.x.copy(), mesh.v.copy()
        mesh.x[:], mesh.v[:] = x, v
        self.assertEqual(solver.resolve(mesh), first)
        np.testing.assert_array_equal(mesh.x, expected_x)
        np.testing.assert_array_equal(mesh.v, expected_v)

    def test_invalid_topology_and_state_are_refused(self):
        with self.assertRaises(ValueError):
            ClothSelfCollision([[0, 1, 5]], 3)
        with self.assertRaises(ValueError):
            ClothSelfCollision([[0, 1, 2]], 3, thickness_m=0.)
        with self.assertRaises(ValueError):
            ClothSelfCollision([[0, 1, 2]], 3, iterations=0)
        solver = ClothSelfCollision([[0, 1, 2]], 3)
        mesh = sheet([[0., 0, 0], [1, 0, 0], [0, 1, 0]], mass=np.ones(3))
        mesh.mass = np.array([1., -1., 1.])
        with self.assertRaises(ValueError):
            solver.resolve(mesh)


if __name__ == '__main__':
    unittest.main()
