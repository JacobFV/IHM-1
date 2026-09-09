"""Hinge bending checks: conservative forces, invariances, and the fold mode.

These assert energy and force properties, not appearance. The decisive one is
``test_isometric_fold_costs_no_spring_energy_but_does_cost_bending_energy``:
it shows the defect this term repairs. A crease that preserves every edge
length is free under distance springs alone, which is why sharp folds appeared
in the simulated blanket rather than only in the render.
"""
import unittest
import numpy as np
from ihm.assembly.cloth_bending import ClothBending, hinges_from_triangles
from ihm.assembly.environment_dynamics import SpringMesh


def flap(height=0.):
    """Two triangles sharing edge (0, 1); vertex 3 lifted by ``height``."""
    return np.array([[0., 0, 0], [1., 0, 0], [.5, -.6, 0], [.5, .4, height]])


FLAP_TRIANGLES = [[0, 1, 2], [1, 0, 3]]


class Tests(unittest.TestCase):
    def setUp(self):
        self.rest = flap()
        self.model = ClothBending(self.rest, FLAP_TRIANGLES, rigidity_n_m=.0025)
        self.mesh = SpringMesh('blanket', 'cloth', [-.42, -.88, 0], [.42, -.04, 0], 1.5)
        self.faces = np.asarray(self.mesh.indices, int).reshape(-1, 3)

    def test_rest_state_is_a_stress_free_minimum(self):
        self.assertLess(self.model.energy(self.rest), 1e-24)
        np.testing.assert_allclose(self.model.forces(self.rest), 0, atol=1e-12)
        self.assertLess(self.model.rest_stencil_residual_m, 1e-12)

    def test_bending_never_resists_in_plane_stretch_or_shear(self):
        # Any affine image of a planar rest flap stays planar, so the hinge term
        # must contribute exactly nothing to stretch, shear or scale. A bending
        # model that leaked into those would double count the spring stiffness.
        transform = np.array([[1.7, .4, 0], [-.3, .6, 0], [0, 0, 1.]])
        stretched = self.rest@transform.T+np.array([3., -2., 5.])
        self.assertLess(self.model.energy(stretched), 1e-20)
        np.testing.assert_allclose(self.model.forces(stretched), 0, atol=1e-10)

    def test_rigid_motion_invariance_and_zero_net_force_and_torque(self):
        angle = .7
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1.]])
        folded = flap(.3)
        moved = folded@rotation.T+np.array([1., 2., 3.])
        self.assertAlmostEqual(self.model.energy(moved), self.model.energy(folded), places=12)
        force = self.model.forces(folded)
        np.testing.assert_allclose(force.sum(axis=0), 0, atol=1e-12)
        np.testing.assert_allclose(np.cross(folded, force).sum(axis=0), 0, atol=1e-12)

    def test_forces_are_the_exact_negative_energy_gradient(self):
        state = flap(.25)
        state[2, 1] += .05
        analytic = self.model.forces(state)
        numeric = np.zeros_like(state)
        step = 1e-6
        for node in range(len(state)):
            for axis in range(3):
                plus = state.copy(); plus[node, axis] += step
                minus = state.copy(); minus[node, axis] -= step
                numeric[node, axis] = -(self.model.energy(plus)-self.model.energy(minus))/(2*step)
        np.testing.assert_allclose(analytic, numeric, atol=1e-7)

    def test_hessian_is_constant_symmetric_and_positive_semidefinite(self):
        dense = self.model.hessian().toarray()
        np.testing.assert_allclose(dense, dense.T, atol=1e-12)
        self.assertGreaterEqual(np.min(np.linalg.eigvalsh(dense)), -1e-10)
        # A constant Hessian handed to an implicit solve cannot create energy.
        np.testing.assert_allclose(self.model.hessian().toarray(), dense)
        state = flap(.3)
        step = 1e-6
        for node in range(len(state)):
            for axis in range(3):
                plus = state.copy(); plus[node, axis] += step
                minus = state.copy(); minus[node, axis] -= step
                column = -(self.model.forces(plus)-self.model.forces(minus)).ravel()/(2*step)
                np.testing.assert_allclose(dense[:, 3*node+axis], column, atol=1e-6)

    def test_energy_matches_the_discrete_shell_hinge_law(self):
        # 3 * B * |edge|^2 * theta^2 / (A0 + A1) for one hinge, to leading order.
        rigidity = .0025
        for angle in (.002, .02):
            state = flap()
            state[3] = [.5, .4*np.cos(angle), .4*np.sin(angle)]
            expected = 3*rigidity*1.**2*angle**2/(.5*1.*.6+.5*1.*.4)
            self.assertAlmostEqual(self.model.energy(state)/expected, 1., places=3)
            measured = np.radians(self.model.dihedral_angles_deg(state)[0])
            self.assertAlmostEqual(measured, angle, places=9)

    def _folded(self, angle):
        """Rotate half the production sheet about a grid line by ``angle``."""
        mesh = self.mesh
        grid = mesh.rest.reshape(19, 23, 3).copy()
        hinge = grid[9, 0, 0]
        offset = grid[10:, :, 0]-hinge
        grid[10:, :, 0] = hinge+offset*np.cos(angle)
        grid[10:, :, 2] = grid[10:, :, 2]+offset*np.sin(angle)
        return grid.reshape(-1, 3)

    def _spring_energy(self, x, material):
        mesh = self.mesh
        length = np.linalg.norm(x[mesh.b]-x[mesh.a], axis=1)
        select = (np.asarray(mesh.edges)[:, 2] >= .99) == material
        return float(.5*np.sum(mesh.stiffness[select]*(length[select]-mesh.length[select])**2))

    def test_isometric_fold_leaves_every_material_edge_and_its_energy_untouched(self):
        # This is the defect the hinge term repairs. Folding half the sheet is a
        # rigid motion of that half, so no yarn edge changes length and the
        # material spring energy is exactly zero: the crease was free.
        mesh = self.mesh
        folded = self._folded(np.pi/2)
        material = np.asarray(mesh.edges)[:, 2] >= .99
        lengths = np.linalg.norm(folded[mesh.b]-folded[mesh.a], axis=1)
        np.testing.assert_allclose(lengths[material], mesh.length[material], atol=1e-15)
        self.assertLess(self._spring_energy(folded, True), 1e-25)
        self.assertGreater(mesh.bending.energy(folded), .05)
        self.assertGreater(mesh.bending.dihedral_angles_deg(folded).max(), 89.)
        self.assertGreater(mesh.bending.state(folded)['folds_over_60_deg'], 20)

    def test_two_apart_springs_answer_a_crease_quartically_not_quadratically(self):
        # The only pre-existing resistance to a crease came from the two-apart
        # springs, whose chord shortens as the square of the fold angle, so
        # their energy falls as its fourth power: at small creases they supply
        # almost nothing. A bending law must fall as the square.
        mesh = self.mesh
        angle = .4
        coarse = self._folded(angle)
        fine = self._folded(angle/2)
        two_apart = self._spring_energy(coarse, False)/self._spring_energy(fine, False)
        hinge = mesh.bending.energy(coarse)/mesh.bending.energy(fine)
        self.assertAlmostEqual(two_apart, 16., delta=1.)
        self.assertAlmostEqual(hinge, 4., delta=.05)
        self.assertGreater(mesh.bending.energy(fine), 4*self._spring_energy(fine, False))

    def test_smooth_drape_stays_much_cheaper_than_a_crease(self):
        mesh = self.mesh
        model = ClothBending(mesh.rest, self.faces, rigidity_n_m=.0025)
        grid = mesh.rest.reshape(19, 23, 3)
        index = np.arange(19)[:, None]*np.ones((1, 23))
        smooth = mesh.rest.copy()
        smooth[:, 2] += (.01*np.cos(np.pi*index/18)).reshape(-1)
        zigzag = mesh.rest.copy()
        zigzag[:, 2] += (.01*(-1.)**index).reshape(-1)
        self.assertLess(model.energy(smooth), 1e-4)
        self.assertGreater(model.energy(zigzag), 100*model.energy(smooth))

    def test_nonplanar_rest_and_broken_topology_are_refused(self):
        with self.assertRaises(ValueError):
            ClothBending(flap(.1), FLAP_TRIANGLES)
        with self.assertRaises(ValueError):
            ClothBending(self.rest, [[0, 1, 2], [1, 0, 3]], rigidity_n_m=-1.)
        with self.assertRaises(ValueError):
            ClothBending(np.array([[0., 0, 0], [1., 0, 0], [2., 0, 0], [3., 0, 0]]),
                         [[0, 1, 2], [1, 0, 3]])
        with self.assertRaises(ValueError):
            hinges_from_triangles([[0, 0, 1]])

    def test_boundary_and_nonmanifold_edges_carry_no_hinge(self):
        self.assertEqual(len(hinges_from_triangles([[0, 1, 2]])), 0)
        self.assertEqual(len(hinges_from_triangles([[0, 1, 2], [1, 0, 3], [0, 1, 4]])), 0)
        self.assertEqual(len(self.model.flaps), 1)

    def test_production_blanket_resolves_its_own_drape(self):
        # A bending length below the vertex spacing would leave creases the mesh
        # cannot represent, so the reported value must not fall under it.
        mesh = self.mesh
        state = mesh.bending.state(mesh.rest, mesh.areal_density_kg_m2)
        spacing = float(np.median(mesh.length[np.asarray(mesh.edges)[:, 2] >= .99]))
        self.assertGreaterEqual(state['bending_length_m'], spacing)
        self.assertEqual(state['folds_over_60_deg'], 0)
        self.assertEqual(state['hinge_count'], 1148)


if __name__ == '__main__':
    unittest.main()
