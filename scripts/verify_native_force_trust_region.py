"""Analytic units, bound and model-quality fixtures for the force-root method."""
import numpy as np
from native_force_trust_region import sensitivity_scales,bounded_step,assess_trial
J=np.array([[100.,2.],[.2,.01],[0.,3.]])
bounds=np.array([[-1.,1.],[-2.,2.]]);q=np.zeros(2);r=np.array([3.,.1,1.])
c,s,_=sensitivity_scales(J,bounds)
d,info=bounded_step(J,r,q,bounds,c,s,.01)
a=assess_trial(r,r+J@d,J,d,s)
assert a['accepted'] and abs(a['reduction_ratio']-1)<1e-12
# Changing force units (N -> mN) must not change the requested coordinate step.
units=np.array([1000.,1.,1000.]);c2,s2,_=sensitivity_scales(J*units[:,None],bounds)
d2,_=bounded_step(J*units[:,None],r*units,q,bounds,c2,s2,.01)
assert np.allclose(d,d2,atol=1e-12)
q=np.array([1.,-2.]);d,_=bounded_step(J,r,q,bounds,c,s,.01)
assert np.all(q+d<=bounds[:,1]+1e-12) and np.all(q+d>=bounds[:,0]-1e-12)
assert not assess_trial(r,2*r,J,d,s)['accepted']
# A KKT bound optimum cannot masquerade as a physical root.
d,_=bounded_step(np.eye(1),[-1.],[1.],[[0.,1.]],[1.],[1.],.01)
assert abs(d[0])<1e-10 and abs(-1.+d[0])>.9
print('PASS: force-unit invariance, exact source bounds, poor prediction rejection and bound residual retained')
from solve_native_force_supine import physical_gate
from copy import deepcopy
names=['pelvis_tx','pelvis_ty','pelvis_tz','pelvis_tilt','pelvis_list','pelvis_rotation']
zero=dict(mobility_names=names,constrained_residual=[0.]*6,udot=[0.]*6,constraint_position_error=0.,constraint_velocity_error=0.,constraint_acceleration_error=0.)
assert physical_gate(zero,760.,1.7)['passed']
for key,value in [('udot',[0.,0.,0.,0.,0.,1.01e-4]),('constrained_residual',[.077,0.,0.,0.,0.,0.]),('constraint_position_error',1.01e-5)]:
    invalid=deepcopy(zero);invalid[key]=value;assert not physical_gate(invalid,760.,1.7)['passed']
print('PASS: original acceleration, normalized support and constraint limits independently enforced')
