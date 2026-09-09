"""Engineering spindle feedback policy for the native muscle-driven plant.

Excitations only: this controller never changes coordinates or applies external
support. It is not the IBM cortical kernel and must be reported separately.
"""
from dataclasses import dataclass, asdict
import math

@dataclass(frozen=True)
class PosturalConfig:
    baseline: float = .03
    length_gain: float = 2.
    velocity_gain: float = .1
    def __post_init__(self):
        for name, lo, hi in [('baseline',0,1),('length_gain',0,100),('velocity_gain',0,10)]:
            value=getattr(self,name)
            if isinstance(value,bool) or not math.isfinite(value) or not lo<=value<=hi:
                raise ValueError('Invalid '+name)

class PosturalController:
    """Length-error servo with normalized native fiber velocity damping."""
    def __init__(self, reference, config=PosturalConfig(), muscle_baselines=None):
        self.config=config
        self.reference={k:float(v['fiber_length_m']) for k,v in reference['muscles'].items()}
        self.baselines=dict(muscle_baselines or {})
        if set(self.baselines)-set(self.reference):raise ValueError('Unknown native muscle baseline')
        if any(not math.isfinite(v) or not 0<=v<=1 for v in self.baselines.values()):raise ValueError('Invalid muscle baseline')
    def commands(self,state):
        if set(state['muscles'])!=set(self.reference):raise ValueError('Native muscle catalog changed')
        commands={}
        for name,m in state['muscles'].items():
            optimal=float(m['optimal_fiber_length_m'])
            if optimal<=0:raise ValueError('Positive optimal fiber length required')
            length=(m['fiber_length_m']-self.reference[name])/optimal
            velocity=m['fiber_velocity_m_s']/optimal
            u=self.baselines.get(name,self.config.baseline)+self.config.length_gain*length+self.config.velocity_gain*velocity
            if not math.isfinite(u):raise ValueError('Nonfinite postural control')
            commands[name]=max(0.,min(1.,u))
        return commands
    def identity(self):
        return {'schema':'ihm.engineering-postural-policy.v1','config':asdict(self.config),
                'muscle_baselines':self.baselines,'reference_fiber_length_m':self.reference,
                'brain_trained':False,'external_support':False,'prescribed_motion':False}

def posture_metrics(initial,final):
    """Endpoints alone are insufficient for a walking claim; retain trajectories."""
    q=final['coordinates'];q0=initial['coordinates']
    return {'elapsed_s':final['time_s']-initial['time_s'],
            'pelvis_height_m':q['pelvis_ty']['value'],
            'forward_displacement_m':q['pelvis_tx']['value']-q0['pelvis_tx']['value'],
            'lateral_displacement_m':q['pelvis_tz']['value']-q0['pelvis_tz']['value'],
            'pelvis_tilt_rad':q['pelvis_tilt']['value'],
            'pelvis_list_rad':q['pelvis_list']['value']}
