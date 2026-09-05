"""Conserved exchange contracts and dimensional physical interface solvers.

Coefficients are supplied explicitly by source models or fitted artifacts. No
universal body interaction coefficients are introduced by these operators.
"""
from dataclasses import dataclass,field
import copy
import math
import numpy as np

RATE_UNITS={'m3/s':('m3',1.),'mL/s':('m3',1e-6),'mL/min':('m3',1e-6/60),'L/min':('m3',1e-3/60),
            'kg/s':('kg',1.),'g/s':('kg',1e-3),'mol/s':('mol',1.),'mmol/s':('mol',1e-3),
            'A':('C',1.),'C/s':('C',1.),'W':('J',1.),'J/s':('J',1.)}
FARADAY=96485.33212;GAS=8.314462618

def finite(x):return np.isfinite(x).all()
def positive_time(dt):
    if isinstance(dt,bool) or not np.isscalar(dt) or not math.isfinite(dt) or dt<0:raise ValueError('duration must be finite and nonnegative')

@dataclass(frozen=True)
class Flux:
    id:str
    quantity:str
    source:int|None
    target:int|None
    rate:float
    unit:str
    reference:str=''

@dataclass
class ExchangeState:
    compartments:tuple
    amounts:dict
    units:dict
    audit:list=field(default_factory=list)
    def __post_init__(self):
        self.compartments=tuple(self.compartments)
        if not self.compartments or len(set(self.compartments))!=len(self.compartments) or set(self.amounts)!=set(self.units):raise ValueError('invalid compartments/quantity units')
        self.amounts={k:np.array(v,float,copy=True) for k,v in self.amounts.items()}
        for k,v in self.amounts.items():
            if v.shape!=(len(self.compartments),) or not finite(v) or self.units[k] not in ('m3','kg','mol','C','J'):raise ValueError('invalid conserved state')
            if self.units[k]!='C' and (v<0).any():raise ValueError('negative extensive amount')
    def transfer(self,fluxes,seconds):
        """Atomic simultaneous transfers; internal edges cancel exactly in ledger.

        None is an explicitly external boundary. Negative rates reverse direction.
        The same physical exchange must appear once, identified by a unique id.
        """
        positive_time(seconds);fluxes=list(fluxes)
        if len({f.id for f in fluxes})!=len(fluxes):raise ValueError('duplicate physical exchange IDs')
        delta={k:np.zeros(len(self.compartments)) for k in self.amounts};external={k:0. for k in self.amounts};internal={k:0. for k in self.amounts}
        for f in fluxes:
            if not f.id or f.quantity not in self.amounts or f.unit not in RATE_UNITS or not math.isfinite(f.rate):raise ValueError('invalid flux')
            unit,scale=RATE_UNITS[f.unit]
            if unit!=self.units[f.quantity]:raise ValueError('flux dimension does not match conserved quantity')
            if f.source==f.target:raise ValueError('exchange must connect different compartments')
            for node in (f.source,f.target):
                if node is not None and (isinstance(node,bool) or not isinstance(node,(int,np.integer)) or not 0<=node<len(self.compartments)):raise ValueError('unknown exchange port')
            amount=f.rate*scale*seconds
            if f.source is not None:delta[f.quantity][f.source]-=amount
            else:external[f.quantity]+=amount
            if f.target is not None:delta[f.quantity][f.target]+=amount
            else:external[f.quantity]-=amount
        updated={k:self.amounts[k]+delta[k] for k in self.amounts}
        for k,v in updated.items():
            if not finite(v) or (self.units[k]!='C' and (v<0).any()):raise ValueError('exchange overdrafts a compartment; shorten step or solve implicitly')
            internal[k]=float(delta[k].sum()-external[k])
            scale=max(float(np.abs(delta[k]).sum()),np.finfo(float).tiny)
            if abs(internal[k])>1e-12*scale:raise ValueError('internal exchange conservation failure')
        result=ExchangeState(self.compartments,updated,self.units.copy(),copy.deepcopy(self.audit))
        result.audit.append({'seconds':seconds,'exchange_ids':[f.id for f in fluxes],'internal_balance':internal,'external_amount':external})
        return result

def hydrostatic_pressure(points_m,reference_pressure_pa,reference_position_m,density_kg_m3,gravity_m_s2):
    p=np.asarray(points_m,float);r=np.asarray(reference_position_m,float);g=np.asarray(gravity_m_s2,float)
    if p.ndim!=2 or p.shape[1]!=3 or r.shape!=(3,) or g.shape!=(3,) or not all(finite(v) for v in (p,r,g)):
        raise ValueError('hydrostatics needs finite three-dimensional SI coordinates')
    if not math.isfinite(reference_pressure_pa) or not math.isfinite(density_kg_m3) or density_kg_m3<=0:raise ValueError('invalid pressure or density')
    with np.errstate(over='raise',invalid='raise'):
        try:result=reference_pressure_pa+density_kg_m3*((p-r)@g)
        except FloatingPointError as e:raise ValueError('hydrostatic pressure overflows supplied scales') from e
    if not finite(result):raise ValueError('nonfinite hydrostatic pressure')
    return result

def bernoulli(x):
    x=np.asarray(x,float);out=np.empty_like(x);small=np.abs(x)<1e-5;high=x>50;low=x<-50;mid=~(small|high|low)
    out[small]=1-x[small]/2+x[small]**2/12
    out[high]=x[high]*np.exp(-x[high]);out[low]=-x[low];out[mid]=x[mid]/np.expm1(x[mid])
    return out

def electrodiffusive_flux(ci,cj,phi_i,phi_j,valence,conductance_m3_s,temperature_k):
    """Scharfetter–Gummel mol/s; conductance is D*area/length, not electric siemens."""
    values=np.array([ci,cj,phi_i,phi_j,valence,conductance_m3_s,temperature_k],float)
    if not finite(values) or ci<0 or cj<0 or not valence or conductance_m3_s<0 or temperature_k<=0:raise ValueError('invalid electrodiffusive parameters')
    psi=valence*FARADAY*(phi_j-phi_i)/(GAS*temperature_k)
    if not math.isfinite(psi):raise ValueError('electrochemical potential overflows supplied scales')
    with np.errstate(over='raise',invalid='raise'):
        try:result=float(conductance_m3_s*(ci*bernoulli(psi)-cj*bernoulli(-psi)))
        except FloatingPointError as e:raise ValueError('electrodiffusive flux overflows supplied scales') from e
    if not math.isfinite(result):raise ValueError('nonfinite ionic flux')
    return result

def advance_ions(amounts_mol,volumes_m3,potentials_v,edges,valences,seconds,temperature_k):
    """Implicit finite-volume transport at fixed potentials; exact species conservation.

    Edges (i,j,D*area/length) define source-specific geometric permeability.
    Potentials are imposed boundary data; electrode/ATP energetics are not solved.
    """
    positive_time(seconds);n=np.asarray(amounts_mol,float);v=np.asarray(volumes_m3,float);phi=np.asarray(potentials_v,float);z=np.asarray(valences,float)
    if n.ndim!=2 or v.shape!=(len(n),) or phi.shape!=v.shape or z.shape!=(n.shape[1],):raise ValueError('invalid ionic dimensions')
    if not all(finite(x) for x in (n,v,phi,z)) or (n<0).any() or (v<=0).any() or (z==0).any() or not math.isfinite(temperature_k) or temperature_k<=0:raise ValueError('invalid ionic state')
    edges=list(edges);out=np.empty_like(n)
    for species,charge in enumerate(z):
        a=np.zeros((len(v),len(v)))
        for i,j,g in edges:
            if isinstance(i,(bool,np.bool_)) or isinstance(j,(bool,np.bool_)) or not isinstance(i,(int,np.integer)) or not isinstance(j,(int,np.integer)) or not 0<=i<len(v) or not 0<=j<len(v) or i==j or not math.isfinite(g) or g<0:raise ValueError('invalid ionic edge')
            psi=charge*FARADAY*(phi[j]-phi[i])/(GAS*temperature_k)
            if not math.isfinite(psi):raise ValueError('electrochemical potential overflows supplied scales')
            kij=g*float(bernoulli(psi))/v[i];kji=g*float(bernoulli(-psi))/v[j]
            a[i,i]-=kij;a[j,i]+=kij;a[j,j]-=kji;a[i,j]+=kji
        out[:,species]=np.linalg.solve(np.eye(len(v))-seconds*a,n[:,species])
    if not finite(out) or (out<0).any():raise ValueError('ionic solve failed positivity')
    balance=out.sum(0)-n.sum(0)
    if not np.allclose(balance,0,atol=1e-11*max(n.sum(),np.finfo(float).tiny),rtol=0):raise ValueError('ionic conservation tolerance exceeded')
    return out,{'mole_balance':balance.tolist(),'charge_balance_C':float(FARADAY*balance@z),'field_status':'imposed potentials; no independent energetic closure'}

class HydraulicNetwork:
    """Closed compliant compartments connected by reciprocal pressure/flow ports.

    volumes[m3], pressure[Pa], compliance[m3/Pa], edge conductance[m3/(Pa*s)].
    Backward Euler solves the coupled pressure response, not sequential overwrites.
    """
    def __init__(self,volumes_m3,pressure_pa,compliance_m3_pa,edges):
        self.v0=np.asarray(volumes_m3,float).copy();self.p0=np.asarray(pressure_pa,float).copy();self.c=np.asarray(compliance_m3_pa,float).copy();self.v=self.v0.copy();self.edges=tuple(edges)
        if self.v.ndim!=1 or not len(self.v) or self.p0.shape!=self.v.shape or self.c.shape!=self.v.shape or not all(finite(x) for x in (self.v,self.p0,self.c)) or (self.v<=0).any() or (self.c<=0).any():raise ValueError('invalid hydraulic state')
        self.l=np.zeros((len(self.v),len(self.v)))
        for i,j,g in self.edges:
            if isinstance(i,(bool,np.bool_)) or isinstance(j,(bool,np.bool_)) or not isinstance(i,(int,np.integer)) or not isinstance(j,(int,np.integer)) or not 0<=i<len(self.v) or not 0<=j<len(self.v) or i==j or not math.isfinite(g) or g<0:raise ValueError('invalid hydraulic port')
            self.l[i,i]+=g;self.l[j,j]+=g;self.l[i,j]-=g;self.l[j,i]-=g
    def advance(self,seconds):
        positive_time(seconds)
        offset=self.p0-self.v0/self.c
        v=np.linalg.solve(np.eye(len(self.v))+seconds*self.l/self.c[None,:],self.v-seconds*self.l@offset)
        if not finite(v) or (v<=0).any():raise ValueError('hydraulic model left positive-volume domain')
        if not np.isclose(v.sum(),self.v.sum(),rtol=1e-10,atol=0):raise ValueError('hydraulic conservation failure')
        pressure=offset+v/self.c;flow=np.array([g*(pressure[i]-pressure[j]) for i,j,g in self.edges])
        self.v=v
        return v.copy(),flow
