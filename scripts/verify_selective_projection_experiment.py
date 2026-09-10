"""Retained-state prototype verification; --run-native adds actual native replay.

The default fixture exercises real canonical projection and production frame
assembly, but its deterministic mechanical transition is not physics evidence.
"""
from copy import deepcopy
from pathlib import Path
import argparse
import hashlib
import json
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.articulated import ArticulatedBodyPlant, CanonicalRegistration
from ihm.assembly.selective_projection import SelectiveProjectionPlant
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG


def reference():
    native = json.loads((ROOT/'data/derived/native-stream-smoke-n6e0pvqi/supine/smoke.json').read_bytes())['initial']
    # This older retained geometry predates metabolic/mass-transfer channels.
    # These declared fixture values exercise frame ownership, not physiology.
    native.update(mass_transfer={'enabled': False}, metabolic_reference={'id': 'fixture'},
                  signed_active_fiber_work_j=0., muscle_heat_energy_j=0.,
                  muscle_metabolic_energy_j=0., signed_active_fiber_power_w=0.)
    return native


class FixtureNative:
    def __init__(self, initial):
        self.state = deepcopy(initial)
        self.history = []
        self.fail = False
        self.closed = False

    def snapshot(self): return deepcopy(self.state)
    def checkpoint(self): return deepcopy((self.state, self.history))
    def restore(self, token): self.state, self.history = deepcopy(token)
    def release(self, token): pass

    def advance(self, dt, forces=(), actuation=None):
        self.history.append(deepcopy((dt, forces, actuation)))
        self.state['time_s'] += dt
        angle = dt * .3
        transform = np.eye(4)
        transform[:3, :3] = [[np.cos(angle), -np.sin(angle), 0.],
                             [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]]
        transform[:3, 3] = [dt*.1, dt*.05, -dt*.02]
        for body in self.state['bodies'].values():
            body['transform_ground'] = (transform @ body['transform_ground']).tolist()
        for field, rate in (('positive_active_fiber_work_j', 3.), ('signed_active_fiber_work_j', 2.),
                            ('muscle_heat_energy_j', 4.), ('muscle_metabolic_energy_j', 6.)):
            self.state[field] += dt*rate
        if self.fail: raise RuntimeError('Injected failure after native mutation')
        return self.snapshot()

    def transfer_mass(self, **payload):
        self.state['mass_kg'] += payload['mass_kg']
        self.state['mass_transfer'] = {'enabled': True, 'increment_kg': payload['mass_kg']}
        return self.snapshot()


def fixture_plant(kind, registration, initial):
    plant = object.__new__(kind)
    plant.registration = registration
    plant.native = FixtureNative(initial)
    plant.garments = None
    plant.state = plant._project(initial, 0.)
    return plant


class Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.initial = reference()
        cls.registration = CanonicalRegistration(
            json.loads((ROOT/'data/derived/canonical/mechanics.json').read_bytes()), cls.initial)
        cls.ids = {group['canonical_bones'][0] for group in cls.registration.groups.values()}
        cls.force_id = next(k for k in cls.registration.ids if k not in cls.registration.named)
        cls.ids.add(cls.force_id)  # Include a spatially routed soft-tissue force.

    def pair(self):
        return tuple(fixture_plant(k, self.registration, self.initial)
                     for k in (ArticulatedBodyPlant, SelectiveProjectionPlant))

    def command(self, state):
        return [{'id': self.force_id, 'point_m': state['entities'][self.force_id]['centroid_m'],
                 'force_n': [2., -3., 1.]}]

    def test_exact_full_frames_force_routes_work_and_checkpoint_replay(self):
        full, selective = self.pair()
        baseline_checkpoint, selective_checkpoint = full.checkpoint(), selective.checkpoint()
        outputs = []
        for i in range(8):
            forces = self.command(full.snapshot())
            a = full.advance(.005, forces, {'tibant_r': .2})
            b = selective.advance_observation(.005, forces, {'tibant_r': .2}, entity_ids=self.ids)
            expected = deepcopy(a)
            expected['entities'] = {k: v for k, v in a['entities'].items() if k in self.ids}
            self.assertEqual(expected, b)
            self.assertEqual(full.native.snapshot(), selective.native.snapshot())
            self.assertEqual(full.native.history, selective.native.history)
            if (i+1) % 4 == 0:
                self.assertEqual(a, selective.snapshot())
                self.assertEqual(a, selective.snapshot())
            outputs.append(b)
        full.restore(baseline_checkpoint)
        selective.restore(selective_checkpoint)
        for expected in outputs:
            forces = self.command(full.snapshot())
            full.advance(.005, forces, {'tibant_r': .2})
            self.assertEqual(expected, selective.advance_observation(.005, forces, {'tibant_r': .2}, entity_ids=self.ids))

    def test_pending_snapshot_and_returned_observations_are_isolated(self):
        full, selective = self.pair()
        expected = full.advance(.005)
        observation = selective.advance_observation(.005, entity_ids=self.ids)
        self.assertIsNone(selective._full_frame)
        observation['entities'].clear()
        observation['native_bodies'].clear()
        # Mutating the owner's native endpoint after capture must not change
        # an already accepted pending projection; restore it after this probe.
        token = selective.native.checkpoint()
        selective.native.state['bodies'].clear()
        self.assertEqual(expected, selective.snapshot())
        selective.native.restore(token)
        exposed = selective.snapshot()
        exposed['entities'].clear()
        exposed['muscles'].clear()
        self.assertEqual(expected, selective.snapshot())

    def test_late_failure_restores_native_and_pending_projection(self):
        full, selective = self.pair()
        full.advance(.005)
        selective.advance_observation(.005, entity_ids=self.ids)
        for mode in ('native', 'projection'):
            self.assertIsNone(selective._full_frame)
            native = selective.native.snapshot()
            selective.native.fail = mode == 'native'
            context = patch.object(selective, '_project_observation', side_effect=RuntimeError('Projection failed'))
            if mode == 'projection': context.start()
            try:
                with self.assertRaises(RuntimeError):
                    selective.advance_observation(.005, entity_ids=self.ids)
            finally:
                selective.native.fail = False
                if mode == 'projection': context.stop()
            self.assertIsNone(selective._full_frame)
            self.assertEqual(native, selective.native.snapshot())
        self.assertEqual(full.snapshot(), selective.snapshot())

    def test_third_substep_outer_rollback_and_mass_refresh(self):
        full, selective = self.pair()
        saved = selective.checkpoint()
        selective.advance_observation(.005, entity_ids=self.ids)
        selective.advance_observation(.005, entity_ids=self.ids)
        selective.native.fail = True
        with self.assertRaises(RuntimeError): selective.advance_observation(.005, entity_ids=self.ids)
        selective.native.fail = False
        selective.restore(saved)
        self.assertEqual(full.snapshot(), selective.snapshot())
        full.advance(.005)
        selective.advance_observation(.005, entity_ids=self.ids)
        self.assertEqual(full.transfer_mass(mass_kg=.01), selective.transfer_mass(mass_kg=.01))
        self.assertEqual(full.snapshot(), selective.snapshot())
        self.assertEqual(full.advance(.005), selective.advance(.005))

    def test_unknown_ids_rejected_before_native_mutation(self):
        _, selective = self.pair()
        before = selective.snapshot()
        with self.assertRaises(ValueError): selective.advance_observation(.005, entity_ids={'typo'})
        self.assertEqual(before, selective.snapshot())
        self.assertFalse(selective.native.history)

    def test_garment_owner_restored_after_projection_failure(self):
        _, selective = self.pair()

        class Garment:
            def __init__(self): self.time = 0.
            def checkpoint(self): return self.time
            def restore(self, value): self.time = value
            def frame(self): return {'fixture_time_s': self.time}
            def advance(self, native, dt, mapped, actuation):
                self.time += dt
                return native.advance(dt, mapped, actuation), {}

        selective.garments = Garment()
        selective.advance_observation(.005, entity_ids=self.ids)
        before = selective.snapshot()
        with patch.object(selective, '_project_observation', side_effect=RuntimeError('Late projection failure')):
            with self.assertRaises(RuntimeError):
                selective.advance_observation(.005, entity_ids=self.ids)
        self.assertEqual(selective.garments.time, .005)
        self.assertEqual(selective.native.snapshot()['time_s'], .005)
        self.assertEqual(before, selective.snapshot())

    def test_control_world_exchange_preserves_exact_public_frame_and_loads(self):
        from ihm.assembly.control_exchange import advance_feedback_exchange

        class Feedback:
            control_interval_s = .01
            def __init__(self): self.time = 0.
            def step(self, dt, observation, **inputs):
                assert observation['time_s'] == self.time
                assert observation['muscles'] and observation['joints']
                self.time += dt
                return {'time_s': self.time, 'motor_excitations': {'tibant_r': .1}}

        world_id = next(iter(self.registration.named))

        class World:
            ids = [world_id]
            def __init__(self): self.time = 0.
            def advance(self, dt, entities):
                self.time += dt
                return [{'id': world_id, 'point_m': entities[world_id]['centroid_m'],
                         'force_n': [1., 2., 3.]}]

        def projector(forces, entities):
            ports = []
            for f in forces:
                jacobian = np.array(entities[f['id']]['rotation_matrix'])[:, 0]
                ports.append({**f, 'generalized_force_pa': float(np.dot(f['force_n'], jacobian))})
            return {'external_pressure_pa': -sum(p['generalized_force_pa'] for p in ports),
                    'force_ports': ports, 'ignored_nonrespiratory_ids': []}

        for with_world in (False, True):
            full, selective = self.pair()
            state = full.snapshot()
            force = self.command(state)
            with patch.object(selective, 'advance_observation', wraps=selective.advance_observation) as narrow:
                expected = advance_feedback_exchange(full, Feedback(), World() if with_world else None,
                    .02, state, force, {}, respiratory_projector=projector)
                actual = advance_feedback_exchange(selective, Feedback(), World() if with_world else None,
                    .02, state, force, {}, respiratory_projector=projector)
            self.assertEqual(expected, actual)
            self.assertEqual(full.snapshot(), selective.snapshot())
            self.assertEqual(full.native.history, selective.native.history)
            self.assertEqual(narrow.call_count, 3 if with_world else 1)
            for call in narrow.call_args_list:
                self.assertEqual(set(call.kwargs['entity_ids']),
                                 {self.force_id, world_id} if with_world else {self.force_id})
            self.assertEqual(len(actual[1]['entities']), len(self.registration.ids))


def native_acceptance():
    """One actual owner, matched checkpoint branches; no concurrent native runs."""
    output = Path(tempfile.mkdtemp(prefix='selective-projection-', dir=ROOT/'data/derived'))
    # Retain source bytes before owner creation, never hash a possibly edited
    # observer at the end of its run. The native owner retains binary identity.
    sources = {Path(__file__).resolve()}
    for name, module in tuple(sys.modules.items()):
        if name.startswith('ihm.') and getattr(module, '__file__', None):
            path = Path(module.__file__).resolve()
            if path.is_relative_to(ROOT) and path.suffix == '.py': sources.add(path)
    hashes = {}
    for path in sorted(sources):
        raw = path.read_bytes()
        relative = path.relative_to(ROOT)
        destination = output/'sources'/relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(raw)
        hashes[str(relative)] = hashlib.sha256(raw).hexdigest()
    (output/'source_manifest.json').write_text(json.dumps(hashes, indent=2)+'\n')
    plant = SelectiveProjectionPlant(ROOT, output/'plant', environment='supine',
                                    target_mass_kg=MECHANICAL_TARGET_MASS_KG,
                                    augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json')
    try:
        initial = plant.snapshot()
        ids = {group['canonical_bones'][0] for group in plant.registration.groups.values()}
        force_id = plant.registration.groups['hand_r']['canonical_bones'][0]
        saved = plant.checkpoint()
        traces = {}
        timings = {}
        for mode in ('full', 'selective'):
            plant.restore(saved)
            traces[mode] = []
            started = time.perf_counter()
            state = initial
            for i in range(8):
                forces = [{'id': force_id, 'point_m': state['entities'][force_id]['centroid_m'],
                           'force_n': [2., 0., 0.]}]
                if mode == 'full' or (i+1) % 4 == 0:
                    state = plant.advance(.005, forces, {'tibant_r': .1})
                else: state = plant.advance_observation(.005, forces, {'tibant_r': .1}, entity_ids=ids)
                narrow = deepcopy(state)
                narrow['entities'] = {k: v for k, v in narrow['entities'].items() if k in ids}
                traces[mode].append({'observation': narrow})
                if (i+1) % 4 == 0: traces[mode][-1]['public_frame'] = plant.snapshot()
            timings[mode] = time.perf_counter()-started
        assert traces['full'] == traces['selective'], 'Actual native replay differs'
        plant.release(saved)
        report = {'passed': True, 'exact_native_frames': True, 'steps': 8,
                  'subset_entities': len(ids), 'full_entities': len(initial['entities']),
                  'wall_s': timings, 'speedup': timings['full']/timings['selective'],
                  'source_sha256_at_start': hashes,
                  'scope': 'Short native checkpoint replay; not sustained world or physiology acceptance'}
        (output/'traces.json').write_text(json.dumps(traces, allow_nan=False)+'\n')
        (output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps({'output': str(output), **report}, indent=2))
    finally: plant.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-native', action='store_true')
    args = parser.parse_args()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests))
    if not result.wasSuccessful(): raise SystemExit(1)
    if args.run_native: native_acceptance()
