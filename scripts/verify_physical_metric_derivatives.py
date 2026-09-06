"""Retained native directional prediction and shoulder force accounting fixtures."""
from diagnose_physical_metric_derivatives import diagnose
r=diagnose()
assert r['product_rule_corrected_vs_full_relative_error']<2e-6
for trial in r['trials']:
    assert trial['full_error_norm']<.11*trial['frozen_error_norm']
for direction in r['recomputed_full_jacobian_directions']:
    assert direction['diagnostic']['predicted_acceleration_cost']<r['initial_cost']
    assert direction['diagnostic']['maximum_linear_support_error']<1e-8
    assert direction['active_manifold_stationarity_norm']<1e-6
    assert not direction['source_bound_active']
for shoulder in r['shoulder_force_attribution']:
    assert abs(shoulder['closure_error_nm'])<1e-10
    assert len(shoulder['muscles'])==3
print('PASS: omitted mass derivative attribution, improved full-Jacobian predictions, QP tangent stationarity and shoulder force closure')
