"""skin network checks against conservation and an analytic RC solution."""
import numpy as np
from ihm import body
from ihm.materialize.skin import SkinPatch, skin_model, electric_field

patch = SkinPatch.line(n=7)
r = body(skin=patch)
m = skin_model(patch, 'skin-subject')
assert 'blood.pressure' not in m.components
assert 'skin.cell_0.membrane_voltage' in m.components
assert r.topologies['skin_gap_junction'].edges != r.topologies['skin_extracellular'].edges
# Exercise the assembled process contributions with unequal capacitances.
from dataclasses import replace
from ihm.materialize.skin import Cell, SurfaceSite
heterogeneous = replace(patch, cells=tuple(replace(c, capacitance_f=c.capacitance_f*(i+1))
                                          for i, c in enumerate(patch.cells)))
registered = body(skin=heterogeneous)
voltage = np.linspace(-.08, -.03, len(heterogeneous.cells))
state = {f'skin.cell_{i}.membrane_voltage': v for i, v in enumerate(voltage)}
charge_rate = 0.
for process in registered.processes.values():
    if process.topology == 'skin_gap_junction':
        cell_index = int(process.output.split('.')[1].split('_')[1])
        rate = sum(w*state[c] for c, w in zip(process.inputs, process.weights)) + process.bias
        charge_rate += heterogeneous.cells[cell_index].capacitance_f*rate
assert abs(charge_rate) < 1e-24
# A single membrane with fixed ion reservoirs has an analytic RC relaxation.
cell = Cell((0., 0., -1e-4))
single = SkinPatch((cell,), (SurfaceSite((0., 0., 0.)),), (), ())
rc = skin_model(single, view='membrane')
i = rc.components.index('skin.cell_0.membrane_voltage')
resting = rc.mean[i]; rc.mean[i] += .01
rc.advance(.1)
rate = (cell.sodium_conductance_s+cell.potassium_conductance_s+cell.chloride_conductance_s)/cell.capacitance_f
assert np.isclose(rc.mean[i], resting+.01*np.exp(-rate*.1), atol=1e-12)
before = m.mean.copy()
m.advance(5)
assert np.isfinite(m.mean).all() and np.linalg.eigvalsh(m.cov).min() > -1e-9
field = electric_field(m, patch)
assert len(field['V_per_m']) == len(patch.surface_edges)
# Uniform intact TEP has no lateral field. A shunt creates opposing gradients.
intact = skin_model(SkinPatch.line(n=7, wound=False), 's')
intact.advance(5)
assert np.max(np.abs(electric_field(intact, SkinPatch.line(n=7, wound=False))['V_per_m'])) < 1e-8
assert min(field['V_per_m']) < 0 < max(field['V_per_m'])
assert not np.array_equal(before, m.mean)
print('verified: heterogeneous skin graphs, finite covariance, intact symmetry, wound field direction')
