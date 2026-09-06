"""Source-derived Thelen effective potential fixture, not native stored energy."""
import csv,hashlib,json,math
from pathlib import Path
import xml.etree.ElementTree as ET
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1]


class ThelenBranch:
    def __init__(self,parameters,activation):
        self.p=parameters;self.activation=activation
        if parameters['pennation_angle_at_optimal']!=0:raise ValueError('This source fixture is restricted to zero-pennation lumbar donors')
    def fiber_force(self,l):
        p=self.p;x=l/p['optimal_fiber_length'];active=math.exp(-(x-1)**2/p['KshapeActive'])
        passive=math.expm1(p['KshapePassive']*(x-1)/p['FmaxMuscleStrain'])/math.expm1(p['KshapePassive']) if x>1 else 0
        return p['max_isometric_force']*(self.activation*active+passive)
    def fiber_effective_potential(self,l):
        p=self.p;x=l/p['optimal_fiber_length'];ka=p['KshapeActive'];kp=p['KshapePassive'];e=p['FmaxMuscleStrain'];d=max(0.,x-1)
        active=self.activation*math.sqrt(math.pi*ka)/2*math.erf((x-1)/math.sqrt(ka))
        passive=(e/kp*math.expm1(kp*d/e)-d)/math.expm1(kp)
        return p['max_isometric_force']*p['optimal_fiber_length']*(active+passive)
    def tendon(self,l):
        p=self.p;x=l/p['tendon_slack_length']-1;e0=p['FmaxTendonStrain'];etoe=99*e0*math.exp(3)/(166*math.exp(3)-67);k=.67/(e0-etoe);c=.33/math.expm1(3)
        if x<=0:f=energy=0.
        elif x<=etoe:f=c*math.expm1(3*x/etoe);energy=c*(etoe/3*math.expm1(3*x/etoe)-x)
        else:
            d=x-etoe;f=.33+k*d;energy=c*(etoe/3*math.expm1(3)-etoe)+.33*d+.5*k*d*d
        return p['max_isometric_force']*f,p['max_isometric_force']*p['tendon_slack_length']*energy
    def equilibrate(self,L):
        ofl=self.p['optimal_fiber_length']
        # Bounded fixture branch only; do not imitate native lower-bound policy.
        lo=.2*ofl;hi=min(2*ofl,L)
        l=brentq(lambda l:self.fiber_force(l)-self.tendon(L-l)[0],lo,hi,xtol=1e-15)
        force,Ut=self.tendon(L-l);h=1e-7
        stiffness=((self.fiber_force(l+h)-self.tendon(L-l-h)[0])-(self.fiber_force(l-h)-self.tendon(L-l+h)[0]))/(2*h)
        if stiffness<=0:raise ValueError('Unstable internal equilibrium branch')
        return dict(force_n=force,effective_potential_j=Ut+self.fiber_effective_potential(l),fiber_length_m=l,internal_stiffness_n_m=stiffness)


def audit():
    model=ROOT/'data/derived/lumbar-muscle-native-lb45uirs/variant/subject_with_lumbar.osim';observations=ROOT/'data/derived/lumbar-muscle-native-lb45uirs/observations.csv';root=ET.parse(model).getroot();models={}
    keys=('pennation_angle_at_optimal','optimal_fiber_length','tendon_slack_length','max_isometric_force','KshapeActive','KshapePassive','FmaxMuscleStrain','FmaxTendonStrain')
    for m in root.iter('Thelen2003Muscle'):
        if m.get('name').startswith('gait2392_'):models[m.get('name')]=ThelenBranch({k:float(m.findtext(k)) for k in keys},.05)
    result=[]
    for row in csv.DictReader(observations.open()):
        branch=models[row['muscle']];L=float(row['length']);state=branch.equilibrate(L);h=1e-7
        derivative=(branch.equilibrate(L+h)['effective_potential_j']-branch.equilibrate(L-h)['effective_potential_j'])/(2*h)
        result.append(dict(muscle=row['muscle'],pose=int(row['pose']),axis=row['axis'],native_force_n=float(row['tendon_force']),native_initialization_tolerance_n=1e-8*branch.p['max_isometric_force'],source_branch_force_n=state['force_n'],force_error_n=state['force_n']-float(row['tendon_force']),energy_length_derivative_error_n=derivative-state['force_n'],effective_torque_derivative_error_nm=(derivative-state['force_n'])*float(row['moment_arm']),internal_stiffness_n_m=state['internal_stiffness_n_m']))
    force_error=max(abs(r['force_error_n']) for r in result);work_error=max(abs(r['energy_length_derivative_error_n']) for r in result)
    sources=[model,observations,Path(__file__),ROOT/'data/raw/mechanics/opensim-core/OpenSim/Actuators/Thelen2003Muscle.cpp']
    return dict(scope='Offline reproduction of six zero-pennation lumbar Thelen fixed-activation static branches at retained native poses; active effective potential is not physiological stored energy',native_runs=0,observations=len(result),maximum_native_force_error_n=force_error,maximum_energy_length_derivative_error_n=work_error,minimum_internal_stiffness_n_m=min(r['internal_stiffness_n_m'] for r in result),rows=result,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})

if __name__=='__main__':print(json.dumps(audit(),indent=2,allow_nan=False))
