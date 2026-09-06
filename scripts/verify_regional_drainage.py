#!/usr/bin/env python3
"""Tiny area/force/ownership boundary fixtures; no native solver run."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ihm.assembly.regional_drainage import RegionalDrainageBoundary
from verify_body_exchange import snapshot


def cells():
    return [{'id': str(i), 'native_owner': 'Skin.extracellular', 'area_m2': area,
             'outward_normal': [0, 0, 1], 'source_identity': {'mesh_sha256': 'fixture', 'triangle': i}}
            for i, area in enumerate([.1, .3])]


def loads():
    return [{'cell_id': r['id'], 'area_m2': r['area_m2'], 'traction_pa': [0, 0, -20]} for r in cells()]


def model():
    return RegionalDrainageBoundary(cells(), {'Skin.extracellular': {
        'complete_native_owner_surface': True, 'source_identity': 'explicit-fixture-scope'}})


def native():
    s = snapshot()
    s['values'].update({'tissue.compression.Skin.requested_pa': 0., 'tissue.compression.Skin.applied_pa': 0.})
    return s


class Checks(unittest.TestCase):
    def test_uniform_area_and_force_conserved(self):
        source, load = native(), loads()
        before = deepcopy((source, load))
        report = model().map(load, source)
        command, = report['native_command_candidates']
        self.assertEqual(command['pressure_pa'], 20.)
        self.assertEqual(command['normal_load_residual_n'], 0.)
        self.assertEqual(report['owner_audits'][0]['normal_load_n'], 8.)
        self.assertEqual(report['owner_audits'][0]['force_n'], [0., 0., -8.])
        self.assertFalse(report['native_commands_sent'])
        self.assertEqual((source, load), before)

    def test_nonuniform_is_not_averaged_into_command(self):
        load = loads(); load[0]['traction_pa'][2] = -40
        report = model().map(load, native())
        self.assertFalse(report['native_command_candidates'])
        self.assertIn('nonuniform_regional_pressure_has_no_native_port', report['unresolved'][0]['reasons'])
        self.assertEqual(report['owner_audits'][0]['normal_load_n'], 10.)

    def test_missing_scope_partial_area_or_port_is_unresolved(self):
        cases = [(RegionalDrainageBoundary(cells()), loads(), native()),
                 (model(), loads()[:1], native()), (model(), loads(), snapshot())]
        for m, load, source in cases:
            report = m.map(load, source)
            self.assertFalse(report['native_command_candidates'])
            self.assertTrue(report['unresolved'])

    def test_duplicate_ownership_and_area_mismatch_fail(self):
        source = cells(); source[1]['source_identity'] = source[0]['source_identity']
        with self.assertRaises(ValueError): RegionalDrainageBoundary(source)
        load = loads(); load[0]['area_m2'] *= 2
        with self.assertRaises(ValueError): model().map(load, native())
        with self.assertRaises(ValueError): model().map(loads()+loads(), native())
        source = cells(); source[0]['outward_normal'] = [0, 0, 2]
        with self.assertRaises(ValueError): RegionalDrainageBoundary(source)

    def test_shear_tension_and_excess_pressure_have_no_command(self):
        for traction in [[1, 0, -20], [0, 0, 20], [0, 0, -5001]]:
            load = loads()
            for r in load: r['traction_pa'] = traction
            self.assertFalse(model().map(load, native())['native_command_candidates'])

    def test_other_native_owner_is_not_skin(self):
        source = cells()
        for r in source: r['native_owner'] = 'Muscle.extracellular'
        m = RegionalDrainageBoundary(source, {'Muscle.extracellular': {
            'complete_native_owner_surface': True, 'source_identity': 'fixture'}})
        report = m.map(loads(), native())
        self.assertFalse(report['native_command_candidates'])
        self.assertIn('owner_has_no_audited_external_pressure_adapter', report['unresolved'][0]['reasons'])


if __name__ == '__main__': unittest.main()
