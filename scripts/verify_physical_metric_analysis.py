"""Sign, ideal-reaction, and physical units fixtures for native metric audit."""
import numpy as np
from verify_native_physical_metric import audit_metric
mass=np.array([[2.,.1],[.1,3.]]);inverse=np.linalg.inv(mass);a=np.array([1.,1.]);r=-mass@a
response=dict(mass_matrix=mass.tolist(),inverse_mass_matrix=inverse.tolist(),udot=a.tolist(),constrained_residual=r.tolist(),tree_residual=(r+np.array([5.,-5.])).tolist())
audit=audit_metric(response)
assert audit['passed'] and audit['tree_without_constraint_acceleration_error']>1.
wrong=dict(response,constrained_residual=(-r).tolist())
assert not audit_metric(wrong)['passed']
print('PASS: constrained sign identity accepted, wrong sign rejected, omitted ideal reaction exposed')
