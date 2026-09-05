#!/usr/bin/env python3
"""Behavioral checks for the reduced thoracic mechanics; no biological validation."""
import copy
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import numpy as np

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))


def verify_native(payload):
    """Exercise every native 50-Hz interval through both mechanical solvers."""
    from ihm.assembly.respiration import BodyRespiration
    from ihm.assembly.mechanics import BodyMechanics
    path = BASE/'data/derived/canonical/native_baseline_v1/body_compartments.csv'
    rows = list(csv.DictReader(path.open()))
    reference = rows.pop(0)
    anatomy = {e['id']: e for e in json.loads((BASE/'data/derived/canonical/anatomy.json').read_text())['entities']}
    columns = {b['entity_id']: ('Left' if 'left' in anatomy[b['entity_id']]['name'].lower() else 'Right')+'LungPulmonaryGasVolume(mL)' for b in payload['bindings'] if b['kind'] == 'lung'}
    body = BodyRespiration.from_dict(payload)
    mechanics = BodyMechanics.from_dict(json.loads((BASE/'data/derived/canonical/mechanics.json').read_text()))
    previous, max_motion, max_error = 0., 0., 0.
    excursions = []
    started = time.perf_counter()
    for index, row in enumerate(rows):
        current = float(row['Time(s)'])
        dt, previous = current-previous, current
        assert abs(dt-.02) < 1e-8, 'Native respiratory verification requires its 50-Hz source clock'
        delta = sum(float(row[c])-float(reference[c]) for c in ['LeftLungPulmonaryGasVolume(mL)', 'RightLungPulmonaryGasVolume(mL)'])*1e-6
        state = body.step(dt, {'volume_change_m3': delta, 'lung_volume_ratios': {id: float(row[c])/float(reference[c]) for id, c in columns.items()}})
        mechanical = mechanics.step(dt, state['mechanics_drivers'])
        assert abs(state['time_s']-current) < 1e-8 and abs(mechanical['time_s']-current) < 1e-8
        assert abs(state['audit']['volume_constraint_residual_m3']) < 1e-12
        assert mechanical['audit']['internal_force_residual_n'] < 1e-6
        max_error = max(max_error, mechanical['audit']['dirichlet_position_residual_m'])
        max_motion = max(max_motion, float(np.max(np.linalg.norm(mechanics.x-mechanics.x0, axis=1))))
        excursions.append(state['diaphragm_descent_m'])
        if (index+1) % 250 == 0:
            print(f'Respiratory native verification {current:.2f}s', flush=True)
    assert previous >= 30.-1e-8 and max_motion < .05 and max_error < 1e-12
    assert max(excursions)-min(excursions) > .005
    return {'samples': len(rows), 'native_duration_s': previous, 'wall_time_s': time.perf_counter()-started, 'csv_sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'maximum_body_translation_m': max_motion, 'maximum_dirichlet_residual_m': max_error, 'diaphragm_peak_to_peak_m': max(excursions)-min(excursions), 'respiratory_energy_residual_j': state['audit']['energy_balance_residual_j']}


def main():
    assert importlib.util.find_spec('ihm.assembly.respiration'), 'Reduced respiratory mechanics is not implemented'
    from ihm.assembly.respiration import BodyRespiration, build, deform_skin
    payload = build(BASE)
    body = BodyRespiration.from_dict(payload)
    zero = body.step(.02, {})
    assert np.linalg.norm(zero['displacement_m']) < 1e-14, 'Unforced chest moved'
    assert abs(zero['audit']['elastic_energy_j']) < 1e-14
    # A sustained imposed pressure must expand the cavity and settle to C*p.
    for _ in range(250):
        pressured = body.step(.02, {'pressure_pa': 98.0665})
    assert pressured['volume_change_m3'] > 0
    assert abs(pressured['volume_change_m3'] - .0002) < 2e-7
    assert np.min(pressured['displacement_m']) > 0
    # Native volume ownership: a volume increase must move actual chest structures.
    body = BodyRespiration.from_dict(payload)
    for i in range(1, 101):
        state = body.step(.02, {'volume_change_m3': .0005 * min(i / 50, 1)})
        assert abs(state['audit']['volume_constraint_residual_m3']) < 1e-13
        assert state['audit']['force_balance_residual_n'] < 1e-8
    assert state['chest_dimensions_m']['transverse'] > payload['reference_dimensions_m']['transverse']
    assert state['chest_dimensions_m']['anteroposterior'] > payload['reference_dimensions_m']['anteroposterior']
    assert .002 < state['diaphragm_descent_m'] < .03
    diaphragm = next(e for e in payload['bindings'] if e['kind'] == 'diaphragm')
    ds = state['entities'][diaphragm['entity_id']]
    assert ds['translation_m'][1] < 0
    assert abs(np.linalg.det(ds['deformation_gradient']) - 1) < 1e-12
    assert state['audit']['accumulated_dissipation_j'] >= 0
    assert abs(state['audit']['energy_balance_residual_j']) < 1e-10
    # Skin must move regionally, without moving head/feet or accumulating drift.
    center = np.array(payload['skin_field']['center_m'])
    points = np.array([center + [0, 0, payload['reference_dimensions_m']['anteroposterior']/2], [0, .8, 0], [0, -.8, 0]])
    moved = deform_skin(points, state['skin_field'])
    assert moved[0, 2] > points[0, 2]
    assert np.array_equal(moved[1:], points[1:])
    assert np.array_equal(points, np.array([center + [0, 0, payload['reference_dimensions_m']['anteroposterior']/2], [0, .8, 0], [0, -.8, 0]]))
    # Lung volume remains native-owned even as anisotropy follows the chest.
    lung = next(e for e in payload['bindings'] if e['kind'] == 'lung')
    state = body.step(.02, {'volume_change_m3': .0005, 'lung_volume_ratios': {lung['entity_id']: 1.2}})
    assert abs(np.linalg.det(state['entities'][lung['entity_id']]['deformation_gradient']) - 1.2) < 1e-12
    # Local external movement redistributes the constrained chest, without gas creation.
    loaded = BodyRespiration.from_dict(payload)
    for _ in range(100):
        load = loaded.step(.02, {'volume_change_m3': 0., 'generalized_forces_n': [0, 1., 0]})
    assert load['displacement_m'][1] > 0 and min(load['displacement_m']) < 0
    assert abs(load['volume_change_m3']) < 1e-13
    # Timestep refinement catches integration bias rather than a fixed display curve.
    ends = []
    for h in [.02, .01, .005]:
        b = BodyRespiration.from_dict(payload)
        b.max_substep = h
        for _ in range(round(.4 / h)):
            s = b.step(h, {'pressure_pa': 98.0665})
        ends.append(np.array(s['displacement_m']))
    assert np.linalg.norm(ends[1]-ends[2]) < np.linalg.norm(ends[0]-ends[2])
    # Global origin shifts must not alter the mechanics or relative shape.
    shifted = copy.deepcopy(payload)
    shift = np.array([2., -3., 5.])
    for bind in shifted['bindings']:
        bind['centroid_m'] = (np.array(bind['centroid_m']) + shift).tolist()
    shifted['skin_field']['center_m'] = (np.array(shifted['skin_field']['center_m'])+shift).tolist()
    for key in ['min', 'max']:
        shifted['skin_field']['bounds_m'][key] = (np.array(shifted['skin_field']['bounds_m'][key])+shift).tolist()
    a = BodyRespiration.from_dict(payload).step(.02, {'pressure_pa': 20.})
    b = BodyRespiration.from_dict(shifted).step(.02, {'pressure_pa': 20.})
    assert np.allclose(a['displacement_m'], b['displacement_m'], atol=1e-14)
    assert np.allclose(deform_skin(points+shift,b['skin_field'])-shift, deform_skin(points,a['skin_field']), atol=1e-14)
    # Dirichlet thorax boundaries must transmit forces to unconstrained neighbors.
    from ihm.assembly.mechanics import BodyMechanics
    mechanics_payload = json.loads((BASE/'data/derived/canonical/mechanics.json').read_text())
    mechanics = BodyMechanics.from_dict(mechanics_payload)
    sternum = next(e['entity_id'] for e in payload['bindings'] if e['kind'] == 'sternum')
    prescribed = {'prescribed_translations_m': {sternum: [0., 0., .001]}}
    mechanical = mechanics.step(.02, prescribed)
    assert abs(mechanical['entities'][sternum]['translation_m'][2]-.001) < 1e-12
    neighbors = {e['b'] if e['a'] == sternum else e['a'] for e in mechanics_payload['links'] if sternum in [e['a'],e['b']]}
    assert max(np.linalg.norm(mechanical['entities'][id]['translation_m']) for id in neighbors) > 1e-8
    assert np.linalg.norm(mechanical['prescribed_reactions_n'][sternum]) > 0
    assert mechanical['audit']['accumulated_prescribed_boundary_work_j'] > 0
    assert mechanical['audit']['dirichlet_position_residual_m'] < 1e-12
    # Anisotropic prescribed lung F must retain native volume and account for energy once.
    mechanical = mechanics.step(.02, {'deformation_gradients': {lung['entity_id']: [[1.2,0,0],[0,1,0],[0,0,1]]}})
    assert np.allclose(mechanical['entities'][lung['entity_id']]['deformation_gradient'], np.diag([1.2,1,1]))
    assert mechanical['audit']['prescribed_affine_boundary_work_j'] > 0
    mx, mt = mechanics.x.copy(), mechanics.time
    for invalid in [{'prescribed_translations_m': {sternum: [0., float('nan'), 0.]}}, {'deformation_gradients': {sternum: np.eye(3).tolist()}}, {'deformation_gradients': {lung['entity_id']: np.diag([-1,1,1]).tolist()}}, {'deformation_gradients': {lung['entity_id']: np.eye(3).tolist()}, 'volume_ratios': {lung['entity_id']: 1.}}]:
        try:
            mechanics.step(.02, invalid)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid Dirichlet/affine boundary accepted')
        assert mechanics.time == mt and np.array_equal(mechanics.x, mx)
    # Rejected drivers may not mutate the shared simulation clock or state.
    previous = (body.time, body.q.copy(), body.v.copy(), body.work)
    for bad in [{'pressure_pa': float('nan')}, {'pressure_pa': 1., 'volume_change_m3': 0.}, {'volume_change_m3': True}, {'generalized_forces_n': [0, 1]}, {'volume_change_m3': 10.}, {'pressure_pa': 1e12}]:
        try:
            body.step(.02, bad)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid respiratory forcing accepted')
        assert body.time == previous[0] and np.array_equal(body.q, previous[1]) and np.array_equal(body.v, previous[2]) and body.work == previous[3]
    report = {'status': 'pass', 'bound_entities': len(payload['bindings']), 'pressure_compliance_response_m3': pressured['volume_change_m3'], 'half_litre_diaphragm_descent_m': state['diaphragm_descent_m'], 'constrained_volume_residual_m3': state['audit']['volume_constraint_residual_m3'], 'energy_balance_residual_j': state['audit']['energy_balance_residual_j'], 'refinement_difference_m': float(np.linalg.norm(ends[1]-ends[2])), 'regional_skin_and_volume_conservation': True, 'empirical_validation': False}
    if '--native' in sys.argv:
        report['native_trajectory'] = verify_native(payload)
    out = BASE/'data/derived/canonical/respiration_verification.json'
    out.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
