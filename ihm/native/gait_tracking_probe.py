"""Small finite donor joint-reference offsets for native muscle-only LQR probes.

Only the controller's desired observations change. No plant coordinates, root
references, muscle seeds, residual actuators, or ground forces are assigned.
"""
import math
import numpy as np
from .gait_reference import GaitReference


class SlowJointReference:
    def __init__(self, reference=None, *, amplitude=.05, start_s=1., duration_s=5., ramp_s=.5):
        self.reference=reference or GaitReference.load()
        if any(n.startswith('pelvis_') for n in self.reference.coordinate_names):raise ValueError('Root coordinates prohibited')
        if not all(math.isfinite(v) for v in (amplitude,start_s,duration_s,ramp_s)) or not 0<=amplitude<=.1 or start_s<0 or duration_s<=0 or not 0<ramp_s<=duration_s:raise ValueError('Invalid bounded joint-reference probe')
        self.amplitude=amplitude;self.start_s=start_s;self.duration_s=duration_s;self.ramp_s=ramp_s
        self.names=tuple(self.reference.coordinate_names)
        self.source_start=float(self.reference.coordinate_data[0,0]);self.source_end=float(self.reference.coordinate_data[-1,0])
        self.q0=self.reference.coordinate_data[0,1:].copy()

    def sample(self,time_s):
        time_s=float(time_s)
        if not math.isfinite(time_s) or time_s<0:raise ValueError('Finite nonnegative time required')
        elapsed=max(0.,min(self.duration_s,time_s-self.start_s));phase=elapsed/self.duration_s
        source_time=self.source_start+phase*(self.source_end-self.source_start)
        q=np.array(list(self.reference.joint_targets(source_time).values()))
        data=self.reference.coordinate_data
        index=min(max(int(np.searchsorted(data[:,0],source_time,side='right'))-1,0),len(data)-2)
        q_slope=(data[index+1,1:]-data[index,1:])/(data[index+1,0]-data[index,0])
        source_rate=(self.source_end-self.source_start)/self.duration_s if self.start_s<time_s<self.start_s+self.duration_s else 0.
        ramp=min(1.,elapsed/self.ramp_s)
        gain=ramp*ramp*(3-2*ramp)
        gain_rate=6*ramp*(1-ramp)/self.ramp_s if 0<elapsed<self.ramp_s else 0.
        offset=self.amplitude*gain*(q-self.q0)
        speed=self.amplitude*(gain_rate*(q-self.q0)+gain*q_slope*source_rate)
        return {'joint_offsets_rad':dict(zip(self.names,map(float,offset))),'joint_target_speed_offsets_rad_s':dict(zip(self.names,map(float,speed))),'recording_phase':phase,'source_time_s':source_time,'ramp':gain}

    def controller_observation(self,snapshot,time_s):
        """Return a separate observation with desired offsets subtracted."""
        target=self.sample(time_s)
        if set(self.names)-set(snapshot['coordinates']):raise ValueError('Native joint catalog lacks donor targets')
        observation=dict(snapshot);observation['coordinates']={n:dict(v) for n,v in snapshot['coordinates'].items()}
        for name,offset in target['joint_offsets_rad'].items():
            observation['coordinates'][name]['value']-=offset
            observation['coordinates'][name]['speed']-=target['joint_target_speed_offsets_rad_s'][name]
        return observation,target
