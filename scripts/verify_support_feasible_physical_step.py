"""Exact support prevents the previous inverse-stiffness merit regression."""
import json
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

# Recorded support-physical-root-xgpiykt1 accepted chain: (max|udot|, sum udot^2, ratio).
# feasible_trial reads only these two invariants of udot, so a 33-component vector carrying
# them exactly replays the native decision without the native solver.
CHAIN=[(5.688100756161002,106.00234988834035,None),(6.079168405925805,85.90940672941792,.8813567701581103),
    (5.942333259077085,83.52096395570562,.5374220790331374),(6.004820636009027,82.24470525999271,.544957027115058),
    (5.8635569757992245,82.00206726451077,.12314061743444306),(5.947271812104625,80.8057378103373,.9871369259241405),
    (5.9038453021011765,80.69980083287598,.45107601088990573)]
SUPPORTED=dict(support_constraints=[0.,0.,0.],gauge_residual=[0.,0.,0.])

def replay(maximum,total,n=33):
    rest=(total-maximum**2)/(n-1)
    assert 0<=rest<=maximum**2
    return [maximum]+[-rest**.5]*(n-1)

def superseded(current,candidate,predicted_reduction):
    # Pre-fix rule: sum-of-squares descent only.
    actual=.5*(np.sum(np.asarray(current['udot'])**2)-np.sum(np.asarray(candidate['udot'])**2))
    ratio=actual/predicted_reduction if predicted_reduction>0 else None
    return bool(max(map(abs,candidate['support_constraints']+candidate['gauge_residual']))<=1e-4
        and predicted_reduction>0 and actual>0 and ratio>=.1)

rows=[]
for i in range(1,len(CHAIN)):
    (m0,s0,_),(m1,s1,ratio)=CHAIN[i-1],CHAIN[i]
    current=dict(udot=replay(m0,s0));candidate=dict(udot=replay(m1,s1),**SUPPORTED)
    predicted=.5*(s0-s1)/ratio
    old=superseded(current,candidate,predicted);new=feasible_trial(current,candidate,predicted)
    assert abs(new['actual_reduction']-.5*(s0-s1))<1e-9 and abs(new['reduction_ratio']-ratio)<1e-9
    assert old,'recorded step was accepted natively'
    assert new['accepted']==(m1<=m0),'new rule tracks the gated maximum'
    rows.append(dict(step=i,maximum=m0,candidate_maximum=m1,sum_of_squares=s0,candidate_sum_of_squares=s1,
        sum_of_squares_fell=s1<s0,maximum_rose=m1>m0,superseded_accepted=old,accepted=new['accepted']))
regressions=[r for r in rows if r['maximum_rose']]
assert len(regressions)==3 and all(r['sum_of_squares_fell'] and r['superseded_accepted'] and not r['accepted'] for r in regressions)
improving=[r for r in rows if not r['maximum_rose']]
assert len(improving)==3 and all(r['accepted'] for r in improving)
assert all(r['sum_of_squares_fell'] for r in rows) and CHAIN[-1][0]>CHAIN[0][0]
print('REPLAY '+json.dumps(rows))
print('PASS: superseded sum-of-squares rule accepts all %d recorded xgpiykt1 steps including the %d that raise max|udot| (%.9f -> %.9f overall); lexicographic rule rejects exactly those and keeps the %d genuinely improving steps'
    %(len(rows),len(regressions),CHAIN[0][0],CHAIN[-1][0],len(improving)))

# Plateau: maximum held by a stuck coordinate, sum of squares still falling, stays acceptable.
plateau=feasible_trial(dict(udot=[5.,4.,4.]),dict(udot=[5.,4.,1.],**SUPPORTED),7.5)
assert plateau['accepted'] and plateau['maximum_reduction']==0
# Strict maximum regression is rejected even at an arbitrarily favourable predicted ratio.
assert not feasible_trial(dict(udot=[5.,5.,5.]),dict(udot=[6.,0.,0.],**SUPPORTED),19.5)['accepted']
# Genuine descent on both quantities still accepted.
assert feasible_trial(dict(udot=[5.,5.,5.]),dict(udot=[1.,1.,1.],**SUPPORTED),36.)['accepted']
print('PASS: plateau tie-break retained, strict max regression rejected, joint descent accepted')
