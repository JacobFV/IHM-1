"""Executable checks of attached RC dynamics, source parity and charge closure."""
from pathlib import Path
from dataclasses import replace
import numpy as np
from scipy.linalg import expm


def verify(root):
    from ihm.assembly.skin_bioelectric import build_skin_electric, electrical_system
    from ihm.materialize.skin import SkinPatch, skin_model
    patch = SkinPatch.line(wound=False)
    a, b, capacitance = electrical_system(patch)
    native = skin_model(patch)
    names = [f'skin.cell_{i}.membrane_voltage' for i in range(9)] + [f'skin.surface_{i}.potential' for i in range(9)]
    indices = [native.components.index(n) for n in names]
    x = np.linspace(-.07, -.01, 18)
    native.mean[indices] = x
    assert np.max(abs((native.a @ native.mean + native.b)[indices] - (a @ x + b))) < 1e-12
    augmented = np.zeros((19, 19)); augmented[:18, :18] = a; augmented[:18, 18] = b
    expected = (expm(augmented*.2) @ np.r_[x, 1.])[:18]
    native.advance(.2)
    assert np.max(abs(native.mean[indices]-expected)) < 1e-12
    # Unequal capacitances expose erroneous symmetric voltage-rate exchange.
    heterogeneous = replace(patch, cells=tuple(replace(c, capacitance_f=c.capacitance_f*(i+1)) for i,c in enumerate(patch.cells)))
    disconnected = replace(heterogeneous, gap_edges=(), surface_edges=())
    ha, _, hc = electrical_system(heterogeneous)
    da, _, _ = electrical_system(disconnected)
    assert abs(float(hc @ ((ha-da) @ x))) < 1e-21
    for current in ([0.], [float('nan')]*9):
        try: electrical_system(patch, current)
        except ValueError: pass
        else: raise AssertionError('invalid electrode accepted')
    artifact = build_skin_electric(root)
    ex = artifact['experiments']
    baseline = ex['baseline']['frames']
    assert np.max(abs(np.array(baseline[-1]['apical_voltage_V']) + .03)) < 1e-12
    assert np.max(abs(np.array(baseline[-1]['membrane_voltage_V']) - baseline[0]['membrane_voltage_V'])) < 1e-12
    wound = ex['wound_shunt']['frames'][-1]
    assert min(wound['edge_field_V_m']) < -1 < 1 < max(wound['edge_field_V_m'])
    assert np.allclose(wound['membrane_voltage_V'], baseline[-1]['membrane_voltage_V'], atol=1e-12)
    assert ex['electrode']['frames'][-1]['apical_voltage_V'] != baseline[-1]['apical_voltage_V']
    assert abs(ex['membrane_perturbation']['frames'][0]['membrane_voltage_V'][0] - baseline[0]['membrane_voltage_V'][0]) > .009
    for experiment in ex.values():
        assert experiment['charge_audit']['max_node_residual_c'] < 1e-18
        assert experiment['charge_audit']['max_internal_pair_residual_a'] < 1e-22
        assert abs(experiment['charge_audit']['global_residual_c']) < 1e-18
    fine = build_skin_electric(root, dt_s=.025)
    for name in ex:
        lhs = ex[name]['frames'][-1]
        rhs = fine['experiments'][name]['frames'][-1]
        for key in ('membrane_voltage_V', 'apical_voltage_V', 'edge_field_V_m'):
            assert np.max(abs(np.array(lhs[key])-rhs[key])) < 1e-10
    geometry = artifact['geometry']
    assert np.asarray(geometry['surface_positions_m']).shape == (9, 3)
    assert np.linalg.norm(np.asarray(geometry['surface_positions_m'])[4]-artifact['anchor']['origin_m']) < 1e-12
    for kwargs in ({'dt_s': 0}, {'dt_s': float('nan')}, {'dt_s': .03}):
        try: build_skin_electric(root, **kwargs)
        except ValueError: pass
        else: raise AssertionError('invalid clock accepted')
    return {'source_generator_parity': True, 'charge_closure': True, 'equilibrium': True,
            'barrier_shunt': True, 'electrode': True, 'membrane_perturbation': True,
            'time_refinement_max_tolerance_V_per_m': 1e-10}


if __name__ == '__main__':
    print(verify(Path(__file__).resolve().parents[1]))
