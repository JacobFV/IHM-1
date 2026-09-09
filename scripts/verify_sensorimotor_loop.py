#!/usr/bin/env python3
"""Mechanical causality and IBM stretch-delay checks; no biological validation."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.mechanical_peripheral import MechanicalPeripheral


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ibm-root', type=Path, default=ROOT.parent / 'IBM-1')
    args = parser.parse_args()
    sys.path.insert(0, str(args.ibm_root))
    from ibm.processes.cord import SegmentalCord

    def plant():
        return MechanicalPeripheral.from_directory(ROOT / 'data/derived/canonical')

    active, passive = plant(), plant()
    # Rejected peripheral commands must not leave the mechanical clock ahead.
    before = active.mechanics.x.copy()
    try:
        active.step(.01, brain_state={'motor_commands': {'proprio:wrong': .5}})
    except ValueError:
        pass
    else:
        raise AssertionError('Unknown motor ID accepted')
    assert active.mechanics.time == active.peripheral.time_s == 0.
    np.testing.assert_array_equal(active.mechanics.x, before)
    # Force alone must not masquerade as spindle afference.
    from ihm.assembly.peripheral import BodyPeripheral
    isolated = BodyPeripheral.from_dict(active.peripheral.data)
    muscle = next(iter(active.ids))
    force_only = isolated.step(.01, mechanical_state={'muscle_forces_n': {muscle: 100.}})
    assert max(force_only['spindle_rates_hz'].values()) == 0.
    assert force_only['tendon_rates_hz'][muscle] > 0.
    ids = sorted(active.ids)
    peak = 0.
    for _ in range(20):
        out = active.step(.01, brain_state={'motor_commands': dict.fromkeys(ids, .5)})
        rest = passive.step(.01)
        peak = max(peak, max(out['spindle_rates_hz'].values()))
        assert max(rest['spindle_rates_hz'].values()) < 1e-9
    assert peak > .01, 'Activation failed to create length feedback'
    assert any(np.linalg.norm(e['translation_m']) > 1e-8
               for e in out['mechanical_state']['entities'].values())
    assert any(np.linalg.norm(np.asarray(e['deformation_gradient']) - np.eye(3)) > 1e-5
               for e in out['mechanical_state']['entities'].values())

    # A real mechanics boundary perturbation, compared with a fresh unperturbed
    # plant and cord. Both start with zero descending drive.
    perturbed, control = plant(), plant()
    c1, c0 = SegmentalCord(ids, dt=.001), SegmentalCord(ids, dt=.001)
    target = next(b for b in perturbed.peripheral.bindings.values()
                  if b['muscle_id'].endswith('recfem_r'))
    axis = np.asarray(target['fiber_axis'])
    gradient = np.eye(3) + .05 * np.outer(axis, axis)
    j = ids.index(target['muscle_id'])
    first_input = first_arc = first_alpha = first_activation = None
    command1 = command0 = None
    for i in range(65):
        drivers = {'deformation_gradients': {target['canonical_entity_id']: gradient.tolist()}}
        o1 = perturbed.step(.001, drivers=drivers, brain_state=command1)
        o0 = control.step(.001, brain_state=command0)
        st1 = np.array([o1['spindle_rates_hz'][m] / 100 for m in ids])
        st0 = np.array([o0['spindle_rates_hz'][m] / 100 for m in ids])
        r1, r0 = c1.step(np.zeros(len(ids)), stretch=st1), c0.step(np.zeros(len(ids)), stretch=st0)
        command1 = {'motor_commands': dict(zip(ids, map(float, r1['alpha'])))}
        command0 = {'motor_commands': dict(zip(ids, map(float, r0['alpha'])))}
        t = (i + 1) * .001
        if first_input is None and st1[j] - st0[j] > 1e-8: first_input = t
        if first_arc is None and r1['stretch'][j] - r0['stretch'][j] > 1e-8: first_arc = t
        if first_alpha is None and r1['alpha'][j] - r0['alpha'][j] > 1e-8: first_alpha = t
        if first_activation is None and o1['motor_activations'].get(ids[j], 0) - o0['motor_activations'].get(ids[j], 0) > 1e-8:
            first_activation = t
        if i < 30: assert abs(r1['stretch'][j] - r0['stretch'][j]) < 1e-8
    assert first_input is not None and first_arc is not None and first_alpha is not None
    assert abs(first_arc - first_input - .030) < 1e-9
    assert first_alpha == first_arc
    assert first_activation is not None and first_activation > first_alpha
    print(json.dumps({'activation_induced_spindle_peak_hz': peak,
                      'first_spindle_s': first_input, 'first_stretch_arc_s': first_arc,
                      'first_alpha_s': first_alpha, 'first_muscle_activation_s': first_activation,
                      'cord_delay_s': first_arc - first_input,
                      'biological_validation': False}, indent=2))


if __name__ == '__main__':
    main()
