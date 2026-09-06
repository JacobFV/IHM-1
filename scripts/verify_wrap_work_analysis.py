"""Analytic conservative path and deliberately inconsistent force fixtures."""
from verify_native_wrap_virtual_work import directional
class Port:
    def muscle(self,m):return {'corrected_passive_energy_j':m['passive_energy_j']}
def state(q,bias=0):
    length=.3+.02*q+.01*q*q
    return dict(independent_names=['arm_flex_r'],independent_q=[q],muscles=[dict(name='arm26_BIClong_r',moment_arms_m=[-(.02+.02*q)+bias],length_m=length,passive_energy_j=10*length,active_effective_energy_j=0,tendon_force_n=10,fiber_length_m=.1,minimum_fiber_length_m=.01,internal_stiffness_n_m=100)])
for h in (1e-3,1e-4,1e-5,-1e-3,-1e-4,-1e-5):
    row=directional(state(.2),state(.2+h),'arm_flex_r',Port())['muscles'][0]
    assert abs(row['trapezoid_virtual_work_error_m_per_rad'])<1e-10
    assert abs(row['effective_energy_length_work_error_nm'])<1e-9
    bad=directional(state(.2,.0002),state(.2+h,.0002),'arm_flex_r',Port())['muscles'][0]
    assert abs(bad['trapezoid_virtual_work_error_m_per_rad']-.0002)<1e-10
print('PASS: sided trapezoid work converges for conservative path and preserves fixed moment-arm discrepancy')
