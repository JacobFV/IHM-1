"""Research-only conservative scalar restraint at one explicitly calibrated pose.

No default human spring, multi-DOF potential, damping, or native registration.
The aggregate periarticular measurement includes labrum and ligamentum teres.
"""
import json,math
from dataclasses import dataclass


@dataclass(frozen=True)
class FixedPoseRestraint:
    center: float
    half_slack: float
    toe: float
    stiffness_per_degree: float
    endpoint_torque: float = 5.

    def __post_init__(self):
        values=(self.center,self.half_slack,self.toe,self.stiffness_per_degree,self.endpoint_torque)
        if not all(math.isfinite(x) for x in values) or self.half_slack<0 or min(self.toe,self.stiffness_per_degree,self.endpoint_torque)<=0:
            raise ValueError('Finite positive calibration required')
        if self.exponent<=1:raise ValueError('Calibration does not support a zero-stiffness continuous toe')

    @property
    def exponent(self):
        return self.stiffness_per_degree*180/math.pi*self.toe/self.endpoint_torque

    @classmethod
    def from_file(cls,path):
        payload=json.loads(path.read_bytes())
        if payload['schema']!='ihm.hip-capsule-fixed-pose.v1' or payload['native_activation_allowed'] is not False:
            raise ValueError('Expected non-native research calibration')
        if payload['calibrated_pose_degrees']!={'flexion':0.,'adduction':0.}:
            raise ValueError('Only retained neutral slice supported')
        p=payload['nominal']
        return cls(math.radians(p['mid_slack_internal_degrees']),math.radians(p['slack_width_degrees']/2),
                   math.radians(p['slack_to_5Nm_degrees']),p['stiffness_at_5Nm_Nm_per_degree'])

    def evaluate(self,internal_rotation_rad,*,flexion_rad,adduction_rad,compression_N=110.,load_medial_angle_degrees=20.):
        """Return (potential J, restoring torque Nm, positive tangent Nm/rad).

        Flexion/abduction are fixed boundary conditions, not unmodeled forces.
        Reject other poses and rotations requiring >5 Nm rather than extrapolate.
        """
        if not all(math.isfinite(q) for q in (internal_rotation_rad,flexion_rad,adduction_rad,compression_N,load_medial_angle_degrees)):
            raise ValueError('Finite coordinates required')
        if abs(compression_N-110.)>1e-12 or abs(load_medial_angle_degrees-20.)>1e-12:
            raise ValueError('Outside retained loading condition; no load-dependence calibration')
        if abs(flexion_rad)>1e-12 or abs(adduction_rad)>1e-12:
            raise ValueError('Outside calibrated fixed flexion/abduction slice; no multi-DOF potential identified')
        displacement=internal_rotation_rad-self.center
        stretch=max(0.,abs(displacement)-self.half_slack)
        if stretch>self.toe+1e-12:raise ValueError('Beyond measured 5 Nm domain; extrapolation forbidden')
        if stretch<=1e-15:return 0.,0.,0.
        ratio=min(stretch/self.toe,1.);p=self.exponent
        torque=self.endpoint_torque*ratio**p
        energy=self.endpoint_torque*self.toe/(p+1)*ratio**(p+1)
        tangent=self.endpoint_torque*p/self.toe*ratio**(p-1)
        return energy,-math.copysign(torque,displacement),tangent
