"""Atomic mechanical owner with opaque rollback and an explicit pressure port.

No native physiological state is touched. Checkpoints roll back this mechanical
state and its work ledger only. Callers must serialize access to one owner.
"""
import copy,uuid
import numpy as np
from .thoracic_free_dynamics import pose_rates
from .thoracic_trajectory import rotation_increment


class ThoracicCompositeOwner:
    def __init__(self,metric,state=None):
        self.metric=metric;self.owner_id=uuid.uuid4().hex;self.revision=0;self._checkpoints={}
        if state is None:state={'q':np.zeros(32),'velocity':np.zeros(38),'rotation':np.eye(3),'translation_m':np.zeros(3),'time_s':0.,'pressure_work_J':0.}
        self._state=self._validate(state)

    def _validate(self,state):
        q=self.metric.coordinates(state['q']).copy();u=np.asarray(state['velocity'],float).copy();r=np.asarray(state['rotation'],float).copy();t=np.asarray(state['translation_m'],float).copy()
        if u.shape!=(38,) or not np.isfinite(u).all() or any(abs(u[6+k])>1e-14 for k in self.metric.thorax.locked):raise ValueError('Invalid composite speeds')
        pose_rates(r,u[:32])
        if t.shape!=(3,) or not np.isfinite(t).all():raise ValueError('Finite parent world translation required')
        time=float(state['time_s']);work=float(state['pressure_work_J'])
        if not np.isfinite([time,work]).all() or time<0:raise ValueError('Invalid mechanical time/work ledger')
        return {'q':q,'velocity':u,'rotation':r,'translation_m':t,'time_s':time,'pressure_work_J':work}

    @property
    def state(self):return copy.deepcopy(self._state)

    def checkpoint(self):
        """Opaque instance-bound token; callers cannot edit saved physical state."""
        if len(self._checkpoints)>=8:raise ValueError('Release a checkpoint before creating more than eight')
        token=self.owner_id+':'+uuid.uuid4().hex;self._checkpoints[token]=self.state;return token

    def restore(self,token):
        if token not in self._checkpoints:raise ValueError('Unknown or foreign mechanical checkpoint')
        self._state=copy.deepcopy(self._checkpoints[token]);self.revision+=1
        return self.state

    def release_checkpoint(self,token):
        if token not in self._checkpoints:raise ValueError('Unknown or foreign mechanical checkpoint')
        del self._checkpoints[token]

    def pressure_port(self,pressure_pa):
        return self.metric.pressure_port(self._state['q'],self._state['velocity'],pressure_pa)

    def advance(self,dt,pressure_pa=0.):
        """Atomic second-order mechanical step at constant outward pressure.

The pressure work ledger uses p*(V_end-V_start), the exact work of this
constant geometric pressure. Midpoint p*Vdot*dt is reported separately so its
quadrature error is visible. No native-flow sign mapping is assumed.
        """
        if isinstance(dt,(bool,np.bool_)) or not np.isscalar(dt) or not np.isfinite(dt) or dt<=0 or not np.isfinite(pressure_pa):raise ValueError('Finite positive dt and finite pressure required')
        start=self.state;q=start['q'];u=start['velocity'];r=start['rotation'];t=start['translation_m']
        initial=self.metric.evaluate(q,u,pressure_pa);half_q=q+.5*dt*u[6:];half_u=u+.5*dt*initial['velocity_derivative'];half_r=r@rotation_increment(.5*dt*u[3:6])
        middle=self.metric.evaluate(half_q,half_u,pressure_pa)
        candidate={**start,'q':q+dt*half_u[6:],'velocity':u+dt*middle['velocity_derivative'],
                   'rotation':r@rotation_increment(dt*half_u[3:6]),'translation_m':t+dt*(half_r@half_u[:3]),'time_s':start['time_s']+dt}
        candidate=self._validate(candidate);end_port=self.metric.pressure_port(candidate['q'],candidate['velocity'],pressure_pa)
        work=float(pressure_pa*(end_port['geometric_volume_m3']-initial['pressure_port']['geometric_volume_m3']))
        candidate['pressure_work_J']+=work
        if not np.isfinite(candidate['pressure_work_J']):raise ValueError('Mechanical work ledger overflow')
        self._state=candidate;self.revision+=1
        return {'state':self.state,'mechanical_pressure_work_J':work,'midpoint_pressure_work_J':float(dt*middle['pressure_port']['mechanical_power_W']),
                'work_quadrature_residual_J':float(work-dt*middle['pressure_port']['mechanical_power_W']),
                'native_state_advanced':False,'revision':self.revision}
