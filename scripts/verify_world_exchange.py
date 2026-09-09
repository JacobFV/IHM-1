"""Subcycling preserves physical time, material stations and interval work."""
import unittest
import numpy as np
from ihm.assembly.world_exchange import advance_body_world


class Tests(unittest.TestCase):
    def test_material_force_impulse_and_work_quadrature(self):
        def state(t):
            return {'time_s': t, 'positive_muscle_work_j': .015,
                    'entities': {'body': {'centroid_m': [t, 0., 0.],
                                         'rotation_matrix': np.eye(3).tolist()}}}
        class Plant:
            def __init__(self): self.t = 0.; self.ports = []
            def advance(self, dt, forces, actuation):
                self.ports.append(forces); self.t += dt
                return state(self.t)
        class World:
            def __init__(self): self.t = 0.; self.seen = []
            def advance(self, dt, entities, **kwargs):
                self.seen.append(entities['body']['centroid_m'][0]); self.t += dt
                return []
        plant, world = Plant(), World()
        command = {'id': 'body', 'point_m': [1., 0., 0.], 'force_n': [2., 0., 0.]}
        out, loads, audit = advance_body_world(plant, world, .02, state(0.), [command], {})
        self.assertEqual(audit['substeps'], 4)
        self.assertEqual(plant.t, world.t)
        self.assertAlmostEqual(out['time_s'], .02)
        self.assertAlmostEqual(out['positive_muscle_work_j'], .06)
        np.testing.assert_allclose([ports[0]['point_m'][0] for ports in plant.ports], [1., 1.005, 1.010, 1.015])
        np.testing.assert_allclose(world.seen, [0., .005, .010, .015])
        np.testing.assert_allclose(np.sum([p['force_n'] for p in loads], axis=0)*.02, [.04, 0., 0.])
        self.assertEqual(command['point_m'], [1., 0., 0.])

    def test_respiratory_virtual_work_uses_each_rotating_material_pose(self):
        from ihm.assembly.embodied_respiration import EmbodiedRespiration
        from scripts.verify_embodied_respiration import fixture, entities
        respiration = EmbodiedRespiration(fixture(), 3000.)
        material = np.array([.04, .02, .03])
        force = np.array([3., -2., 5.])

        def state(t):
            angle = t / .005 * np.pi / 6
            rotation = np.array([[np.cos(angle), -np.sin(angle), 0.],
                                 [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
            pose = entities(rotation)
            for entity in pose.values():
                entity['centroid_m'][0] += 2*t
                entity['translation_m'][0] += 2*t
            return {'time_s': t, 'positive_muscle_work_j': 0., 'entities': pose}

        class Plant:
            t = 0.
            def advance(self, dt, forces, actuation):
                self.t += dt
                return state(self.t)

        class World:
            def advance(self, dt, entities): return []

        initial = state(0.)
        command = {'id': 'lung', 'force_n': force.tolist(),
                   'point_m': respiration.point_position('lung', material, 3000., initial['entities']).tolist()}
        endpoint, historical_loads, audit = advance_body_world(
            Plant(), World(), .02, initial, [command], {},
            respiratory_projector=lambda ports, pose: respiration.project_load(ports, pose, 3000.))
        load = audit['respiratory_load']
        # Independent finite-difference virtual work, at the material station
        # and articulation used by each force sample.
        expected = []
        for t in (0., .005, .01, .015):
            pose = state(t)['entities']
            epsilon_ml = .001
            jacobian = (respiration.point_position('lung', material, 3000.+epsilon_ml, pose)
                        - respiration.point_position('lung', material, 3000.-epsilon_ml, pose)) / (2*epsilon_ml*1e-6)
            expected.append(-force @ jacobian)
        self.assertAlmostEqual(load['external_pressure_pa'], float(np.mean(expected)), places=6)
        stale = respiration.project_load(historical_loads, endpoint['entities'], 3000.)
        self.assertGreater(abs(stale['external_pressure_pa']-load['external_pressure_pa']), .1)
        self.assertAlmostEqual(sum(p['generalized_force_pa'] for p in load['force_ports']),
                               -load['external_pressure_pa'])
        for port in load['force_ports']:
            self.assertAlmostEqual(np.dot(port['force_n'], port['point_jacobian_m_per_m3']),
                                   port['generalized_force_pa'])
        np.testing.assert_allclose([p['sample_time_s'] for p in load['force_ports']], [0., .005, .01, .015])
        self.assertFalse(load['passive_recoil_added'])
        self.assertEqual(load['limitations'], stale['limitations'])


if __name__ == '__main__': unittest.main()
