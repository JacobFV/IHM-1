"""Physical conservation, hydrostatics, charge and two-way interface checks."""
import numpy as np
from ihm.coupling import ExchangeState,Flux,hydrostatic_pressure,electrodiffusive_flux,advance_ions,HydraulicNetwork

s=ExchangeState(('vascular','interstitium','lymph'),{'volume':np.array([5.,12.,2.])*1e-3},{'volume':'m3'})
r=s.transfer([Flux('filtration','volume',0,1,2,'mL/min'),Flux('lymph-uptake','volume',1,2,2,'mL/min'),Flux('return','volume',2,0,2,'mL/min')],60)
assert np.allclose(r.amounts['volume'],s.amounts['volume'])
assert abs(r.audit[-1]['internal_balance']['volume'])<1e-15
try:s.transfer([Flux('wrong','volume',0,1,2,'kg/s')],1)
except ValueError:pass
else:raise AssertionError('dimension mismatch accepted')
try:s.transfer([Flux('overdraw','volume',0,1,99,'m3/s')],1)
except ValueError:pass
else:raise AssertionError('overdraw accepted')
assert np.allclose(s.amounts['volume'],[.005,.012,.002])
p=hydrostatic_pressure([[0,0,0],[0,-1,0]],10000,[0,0,0],1050,[0,-9.80665,0])
assert np.allclose(p,[10000,10000+1050*9.80665])
F=96485.33212;R=8.314462618;T=310.;dv=.02
cj=100*np.exp(-F*dv/(R*T))
assert abs(electrodiffusive_flux(100,cj,0,dv,1,1e-12,T))<1e-20
n=np.array([[.1,.1],[.2,.2]])
new,audit=advance_ions(n,[.001,.002],[0,0],[(0,1,1e-6)],[1,-1],1e6,T)
assert (new>=0).all() and np.allclose(new.sum(0),n.sum(0))
assert np.max(np.abs(audit['mole_balance']))<1e-14
h=HydraulicNetwork([.001,.002],[100,0],[1e-5,2e-5],[(0,1,1e-5)])
v,flow=h.advance(10)
assert np.isclose(v.sum(),.003) and v[0]<.001 and flow[0]>0
print('verified: fluid mass/dimension/positivity, hydrostatics, electrochemical equilibrium, ionic conservation, two-way hydraulic exchange')
for invalid in [lambda:electrodiffusive_flux(1,1,0,1e308,1,1e-12,310),lambda:hydrostatic_pressure([[1e308,0,0]],0,[0,0,0],1000,[10,0,0]),lambda:advance_ions(n,[.001,.002],[0,0],[(True,0,1e-6)],[1,-1],1,T)]:
    try:invalid()
    except ValueError:pass
    else:raise AssertionError('nonfinite derived physics or boolean node accepted')
