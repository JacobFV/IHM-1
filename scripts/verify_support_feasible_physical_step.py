"""Exact support prevents the previous inverse-stiffness merit regression."""
import numpy as np
from support_feasible_physical_step import physical_step,feasible_trial,native_metric
# Unconstrained descent would trade support for a much smaller limb residual.
a=np.array([0.,10.]);J=np.array([[100.,0.],[100.,1.]])
d,diag=physical_step(a,J,[0.,0.],[[-1.,1.],[-1.,1.]],[0.],[[100.,0.]],.03)
assert abs(d[0])<1e-10 and abs(d[1]+.03)<1e-8
current=dict(udot=[0.,10.]);candidate=dict(udot=[0.,1.],support_constraints=[.18,0.,0.],gauge_residual=[0.,0.,0.])
assert not feasible_trial(current,candidate,49.)['accepted']
candidate['support_constraints']=[0.,0.,0.]
assert feasible_trial(current,candidate,49.)['accepted']
response=dict(mass_matrix=[[2.,0.],[0.,3.]],inverse_mass_matrix=[[.5,0.],[0.,1/3]],constrained_residual=[-2.,-3.],tree_residual=[3.,-8.],udot=[1.,1.])
assert np.allclose(native_metric(response),[1.,1.])
print('PASS: hard support blocks tradeoff, nonlinear infeasibility rejected despite lower cost, native metric sign verified')
from diagnose_support_feasible_physical_step import diagnose
retained=diagnose()
assert retained['anchor_q_difference']==0 and retained['anchor_force_difference']==0
assert all(row['diagnostic']['maximum_linear_support_error']<1e-8 and row['diagnostic']['predicted_acceleration_cost']<retained['initial_acceleration_cost'] for row in retained['rows'])
print('PASS: all four retained98 bounded QPs preserve linear support and predict lower physical cost')
