"""Muscle-only engineering joint feedback using retained fitted path derivatives.

This is neither a learned brain nor evidence of walking. Allocation approximates
strength by isometric force times a Gaussian active force-length curve; it ignores
pennation, tendon compliance, passive force, force-velocity and activation lag.
"""
from pathlib import Path
import xml.etree.ElementTree as ET
import math
import numpy as np


def polynomial_exponents(dimension, order):
    """OpenSim ordering: last coordinate varies fastest, total degree <= order."""
    if dimension < 1 or order < 0:
        raise ValueError('Positive dimension and nonnegative order required')
    if dimension == 1:
        return [(i,) for i in range(order + 1)]
    return [(i,) + tail for i in range(order + 1)
            for tail in polynomial_exponents(dimension - 1, order - i)]


class FittedMomentArms:
    def __init__(self, path):
        self.paths = {}
        for item in ET.parse(path).findall('.//FunctionBasedPath'):
            name = item.attrib['name'].split('/')[-1]
            coordinates = tuple(x.split('/')[-1] for x in item.findtext('coordinate_paths').split())
            function = item.find('length_function/MultivariatePolynomialFunction')
            if function is None:
                raise ValueError('Only explicit polynomial path functions are supported')
            coefficients = np.array([float(x) for x in function.findtext('coefficients').split()])
            exponents = np.array(polynomial_exponents(int(function.findtext('dimension')), int(function.findtext('order'))))
            if len(coefficients) != len(exponents) or exponents.shape[1] != len(coordinates) or not np.isfinite(coefficients).all():
                raise ValueError('Malformed fitted polynomial')
            if name in self.paths:
                raise ValueError('Duplicate fitted muscle')
            self.paths[name] = (coordinates, coefficients, exponents)
        if not self.paths:
            raise ValueError('No fitted muscle paths')

    def length_and_moment_arms(self, name, coordinates):
        names, coefficients, exponents = self.paths[name]
        q = np.array([coordinates[n]['value'] if isinstance(coordinates[n], dict) else coordinates[n] for n in names], dtype=float)
        if not np.isfinite(q).all():
            raise ValueError('Nonfinite coordinates')
        length = float(coefficients @ np.prod(q[None, :] ** exponents, axis=1))
        arms = {}
        for j, coordinate in enumerate(names):
            active = exponents[:, j] > 0
            powers = exponents[active].copy()
            multipliers = powers[:, j].copy()
            powers[:, j] -= 1
            arms[coordinate] = -float((coefficients[active] * multipliers) @ np.prod(q[None, :] ** powers, axis=1))
        return length, arms


def allocate_excitation(matrix, torques, baseline=0.02, regularization=1.0):
    """Bounded least squares; torque residual stays observable, never overridden."""
    matrix, torques = np.asarray(matrix, dtype=float), np.asarray(torques, dtype=float)
    if matrix.ndim != 2 or torques.shape != (matrix.shape[0],) or not np.isfinite(matrix).all() or not np.isfinite(torques).all():
        raise ValueError('Finite compatible torque matrix required')
    n = matrix.shape[1]
    baseline = np.broadcast_to(np.asarray(baseline,dtype=float),(n,)).copy()
    if not np.isfinite(baseline).all() or np.any(baseline<0) or np.any(baseline>1) or not math.isfinite(regularization) or regularization <= 0:
        raise ValueError('Invalid allocation regularization/baseline')
    augmented = np.vstack((matrix, math.sqrt(regularization) * np.eye(n)))
    rhs = np.r_[torques, math.sqrt(regularization) * baseline]
    try:
        from scipy.optimize import lsq_linear
        result = lsq_linear(augmented, rhs, bounds=(0, 1), method='bvls', tol=1e-7, max_iter=100)
        if not result.success:
            raise RuntimeError('Muscle allocation failed: ' + result.message)
        # BVLS may leave machine-roundoff excursions at an active bound.
        # Reject materially infeasible results, then canonicalize the wire values.
        if not np.isfinite(result.x).all() or np.any(result.x < -1e-10) or np.any(result.x > 1 + 1e-10):
            raise RuntimeError('Muscle allocation returned infeasible excitation')
        return np.clip(result.x, 0., 1.)
    except ImportError:
        # Projected gradient minimizes the same strictly convex objective.
        hessian = augmented.T @ augmented
        gradient_rhs = augmented.T @ rhs
        step = 1 / np.linalg.norm(hessian, ord=2)
        value = baseline.copy()
        for _ in range(3000):
            updated = np.clip(value - step * (hessian @ value - gradient_rhs), 0, 1)
            if np.max(np.abs(updated - value)) < 1e-8:
                return updated
            value = updated
        return value


def native_center_of_mass(state):
    """Ground-frame COM and velocity from effective native body mass properties."""
    total = 0.; position = np.zeros(3); velocity = np.zeros(3)
    for body in state['bodies'].values():
        mass = float(body['mass_kg'])
        transform = np.asarray(body['transform_ground'], dtype=float)
        local = np.asarray(body['mass_center_local_m'], dtype=float)
        omega = np.asarray(body['angular_velocity_rad_s'], dtype=float)
        origin_velocity = np.asarray(body['origin_velocity_m_s'], dtype=float)
        if mass < 0 or not math.isfinite(mass) or transform.shape != (4,4) or any(x.shape != (3,) for x in (local,omega,origin_velocity)) or not all(np.isfinite(x).all() for x in (transform,local,omega,origin_velocity)):
            raise ValueError('Invalid native COM input')
        offset = transform[:3,:3] @ local
        position += mass * (transform[:3,3] + offset)
        velocity += mass * (origin_velocity + np.cross(omega,offset))
        total += mass
    if total <= 0:
        raise ValueError('Positive native body mass required')
    return position / total, velocity / total


def foot_support_center(state):
    """Vertical-load weighted foot contact center, or None without foot support."""
    points=[]; weights=[]
    for contact in state.get('contacts', []):
        if contact.get('body_frame') in ('calcn_r','calcn_l','toes_r','toes_l'):
            load=max(0.,float(contact['force_n'][1]))
            if load>0:
                points.append(contact['center_m']); weights.append(load)
    return np.average(points,axis=0,weights=weights) if weights else None


class JointPosturalController:
    """Joint PD torque requests realized exclusively through muscle excitations.

    Pelvis feedback modifies the hip and ankle torque requests. Gains and signs
    are engineering policy choices requiring closed-loop evaluation. Unmodeled
    muscles retain baseline excitation and no global-root actuator is introduced.
    """
    def __init__(self, reference, path=None, kp=100.0, kd=15.0,
                 baseline=0.02, pelvis_gain=60.0, pelvis_damping=10.0,
                 regularization=1.0, targets=None, feedforward_torques=None,
                 native_moment_arms=None, com_position_gain=0., com_velocity_gain=0., com_target_x_m=None,
                 com_lateral_position_gain=0., com_lateral_velocity_gain=0., com_target_z_m=None,
                 equilibrium_excitations=None):
        path = path or Path(__file__).resolve().parents[2] / 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled_FunctionBasedPathSet.xml'
        self.geometry = FittedMomentArms(path)
        self.reference_time_s = reference.get('time_s')
        self.reference_muscles = {n:dict(m) for n,m in reference['muscles'].items()}
        self.muscles = tuple(reference['muscles'])
        self.fitted = tuple(n for n in self.muscles if n in self.geometry.paths)
        self.joints = tuple(sorted({q for n in self.fitted for q in self.geometry.paths[n][0]}))
        self.reference_coordinates = {q: float(v['value']) for q,v in reference['coordinates'].items()}
        self.targets = {q: self.reference_coordinates[q] for q in self.joints}
        self.native_moment_arms = {}
        self.native_moment_arms_time_s = None
        self.controlled = self.fitted
        if native_moment_arms is not None:
            self.update_native_moment_arms(native_moment_arms,time_s=self.reference_time_s)
        if set(targets or {}) - set(self.joints):
            raise ValueError('Target coordinate lacks fitted muscle actuation')
        self.targets.update(targets or {})
        self.feedforward = dict(feedforward_torques or {})
        if set(self.feedforward) - set(self.joints):
            raise ValueError('Feedforward coordinate lacks fitted muscle actuation')
        self.kp, self.kd, self.baseline = kp, kd, baseline
        self.pelvis_gain, self.pelvis_damping = pelvis_gain, pelvis_damping
        self.regularization = regularization
        self.com_position_gain, self.com_velocity_gain = com_position_gain, com_velocity_gain
        self.com_lateral_position_gain, self.com_lateral_velocity_gain = com_lateral_position_gain, com_lateral_velocity_gain
        self.com_enabled = any((com_position_gain,com_velocity_gain,com_lateral_position_gain,com_lateral_velocity_gain))
        if not all(math.isfinite(x) and x >= 0 for x in (com_position_gain,com_velocity_gain,com_lateral_position_gain,com_lateral_velocity_gain)):
            raise ValueError('Nonnegative finite COM gains required')
        if com_target_x_m is not None and (isinstance(com_target_x_m,bool) or not math.isfinite(com_target_x_m)):
            raise ValueError('Finite optional COM target required')
        if com_target_z_m is not None and (isinstance(com_target_z_m,bool) or not math.isfinite(com_target_z_m)):
            raise ValueError('Finite optional lateral COM target required')
        initial_com = native_center_of_mass(reference)[0] if self.com_enabled else None
        self.initial_com_x = float(initial_com[0]) if self.com_enabled else None
        self.initial_com_z = float(initial_com[2]) if self.com_enabled else None
        self.com_target_z_m = self.initial_com_z if com_target_z_m is None else float(com_target_z_m)
        self.com_target_x_m = self.initial_com_x if com_target_x_m is None else float(com_target_x_m)
        support = foot_support_center(reference)
        self.initial_support_center = support.tolist() if support is not None else None
        if not all(math.isfinite(x) and x >= 0 for x in (kp,kd,baseline,pelvis_gain,pelvis_damping)) or baseline > 1:
            raise ValueError('Invalid joint control gains')
        if not all(math.isfinite(x) for x in [*self.targets.values(), *self.feedforward.values()]):
            raise ValueError('Nonfinite target or feedforward torque')
        self.equilibrium_excitations = None if equilibrium_excitations is None else dict(equilibrium_excitations)
        if self.equilibrium_excitations is not None:
            if set(self.equilibrium_excitations)!=set(self.muscles) or any(isinstance(v,bool) or not math.isfinite(v) or not 0<=v<=1 for v in self.equilibrium_excitations.values()):
                raise ValueError('Equilibrium excitations require the exact native muscle catalog and finite[0,1] values')
        self.initial_moment_arms = {n:self.geometry.length_and_moment_arms(n,reference['coordinates'])[1] for n in self.fitted}
        for n,arms in self.native_moment_arms.items():self.initial_moment_arms.setdefault(n,{}).update(arms)
        self.last_allocation = {}

    def update_native_moment_arms(self, moment_arms, *, time_s=None):
        """Refresh muscle→coordinate→signed meter arms measured by native model.

        Caller must refresh with body motion; last supplied geometry is held until
        the next update. Missing previously registered muscles/joints is rejected
        to prevent silently dropping actuation. Native values override fitted ones.
        """
        if time_s is not None and (isinstance(time_s,bool) or not math.isfinite(time_s) or time_s<0):
            raise ValueError('Finite nonnegative native moment-arm timestamp required')
        if not isinstance(moment_arms, dict) or set(moment_arms)-set(self.muscles):
            raise ValueError('Native moment arms contain unknown muscles')
        updated = {}
        for name, row in moment_arms.items():
            if not isinstance(row, dict) or not row:
                raise ValueError('Native moment arm row must contain coordinates')
            updated[name] = {}
            for coordinate, value in row.items():
                if coordinate not in self.reference_coordinates or coordinate.startswith('pelvis_'):
                    raise ValueError('Unknown or global-root native moment arm coordinate')
                if isinstance(value, bool) or not math.isfinite(value):
                    raise ValueError('Finite native moment arm required')
                updated[name][coordinate] = float(value)
        for name, previous in self.native_moment_arms.items():
            if name not in updated or set(previous)-set(updated[name]):
                raise ValueError('Native moment arm update removed registered actuation')
        self.native_moment_arms = updated
        self.native_moment_arms_time_s = None if time_s is None else float(time_s)
        if hasattr(self,'initial_moment_arms') and time_s is not None and time_s==self.reference_time_s:
            for n,arms in updated.items():self.initial_moment_arms.setdefault(n,{}).update(arms)
        self.controlled = tuple(n for n in self.muscles if n in self.geometry.paths or n in updated)
        self.joints = tuple(sorted({q for n in self.fitted for q in self.geometry.paths[n][0]} | {q for row in updated.values() for q in row}))
        for q in self.joints:
            self.targets.setdefault(q, self.reference_coordinates[q])

    def commands(self, state):
        if set(state['muscles']) != set(self.muscles):
            raise ValueError('Native muscle catalog changed')
        q = state['coordinates']
        requested = {n: self.kp * (self.targets[n] - q[n]['value']) - self.kd * q[n]['speed'] + self.feedforward.get(n, 0.) for n in self.joints}
        pitch = self.pelvis_gain*(q['pelvis_tilt']['value']-self.reference_coordinates['pelvis_tilt']) + self.pelvis_damping*q['pelvis_tilt']['speed']
        roll = self.pelvis_gain*(q['pelvis_list']['value']-self.reference_coordinates['pelvis_list']) + self.pelvis_damping*q['pelvis_list']['speed']
        com_torque = 0.; lateral_torque = 0.
        com_diagnostics = {}
        if self.com_enabled:
            com, com_velocity = native_center_of_mass(state)
            com_torque = -self.com_position_gain*(com[0]-self.com_target_x_m)-self.com_velocity_gain*com_velocity[0]
            # Native matched-pulse evidence: +right/-left hip torque moves COM
            # toward +z even while pelvis list decreases. Translational COM
            # feedback therefore has the opposite sign from list feedback.
            lateral_torque = -self.com_lateral_position_gain*(com[2]-self.com_target_z_m)-self.com_lateral_velocity_gain*com_velocity[2]
            support = foot_support_center(state)
            com_diagnostics = {'position_ground_m':com.tolist(),'velocity_ground_m_s':com_velocity.tolist(),'target_x_ground_m':self.com_target_x_m,'ankle_correction_total_nm':float(com_torque),'target_z_ground_m':self.com_target_z_m,'hip_lateral_correction_total_nm':float(lateral_torque),'support_center_ground_m':None if support is None else support.tolist()}
        # +pelvis tilt is +z/backward lean. +ankle torque acts -z on the
        # planted tibia and corrects toward +x; its vestibular sign is positive.
        for side, sign in [('r', 1), ('l', -1)]:
            for name, correction in [(f'hip_flexion_{side}', pitch), (f'ankle_angle_{side}', .5*pitch + .5*com_torque), (f'hip_adduction_{side}', sign*(roll+.5*lateral_torque))]:
                if name in requested:
                    requested[name] += correction
        matrix = np.zeros((len(self.joints), len(self.controlled)))
        rows = {n: i for i, n in enumerate(self.joints)}
        current_moment_arms = {}
        for j, name in enumerate(self.controlled):
            muscle = state['muscles'][name]
            optimum = float(muscle['optimal_fiber_length_m'])
            strength = float(muscle['max_isometric_force_n'])
            if not math.isfinite(optimum) or optimum <= 0 or not math.isfinite(strength) or strength <= 0:
                raise ValueError('Positive muscle strength and fiber optimum required')
            normalized = float(muscle['fiber_length_m']) / optimum
            strength *= math.exp(-((normalized-1)/.45)**2)
            arms = self.geometry.length_and_moment_arms(name, q)[1] if name in self.geometry.paths else {}
            arms.update(self.native_moment_arms.get(name, {}))
            current_moment_arms[name] = dict(arms)
            for coordinate, arm in arms.items():
                matrix[rows[coordinate], j] = strength * arm
        torque = np.array([requested[n] for n in self.joints])
        reference_activation = np.array([self.equilibrium_excitations[n] for n in self.controlled]) if self.equilibrium_excitations is not None else None
        equilibrium_torque = matrix @ reference_activation if reference_activation is not None else np.zeros(len(torque))
        activation = allocate_excitation(matrix, torque+equilibrium_torque, self.baseline if reference_activation is None else reference_activation, self.regularization)
        commands = dict.fromkeys(self.muscles, self.baseline) if self.equilibrium_excitations is None else dict(self.equilibrium_excitations)
        commands.update(zip(self.controlled, map(float, activation)))
        achieved = matrix @ activation
        self.last_allocation = {'com_feedback':com_diagnostics, 'requested_torque_nm':requested, 'estimated_torque_nm':dict(zip(self.joints,map(float,achieved-equilibrium_torque))),'estimated_total_torque_nm':dict(zip(self.joints,map(float,achieved))), 'residual_norm_nm':float(np.linalg.norm(achieved-equilibrium_torque-torque))}
        def measured_torques(muscles, moment_arms):
            missing_forces=[n for n in self.controlled if 'tendon_force_n' not in muscles[n]]
            missing_geometry=[n for n in self.controlled if n not in moment_arms or set(current_moment_arms[n])-set(moment_arms[n])]
            if missing_forces or missing_geometry:return None,missing_forces,missing_geometry
            torques=dict.fromkeys(self.joints,0.)
            for n in self.controlled:
                force=float(muscles[n]['tendon_force_n'])
                if not math.isfinite(force):raise ValueError('Nonfinite measured tendon force')
                for joint,arm in moment_arms[n].items():torques[joint]+=force*arm
            return torques,[],[]
        actual,missing_forces,missing_geometry=measured_torques(state['muscles'],current_moment_arms)
        initial_actual,initial_missing_forces,initial_missing_geometry=measured_torques(self.reference_muscles,self.initial_moment_arms)
        actual_delta=None if actual is None or initial_actual is None else {n:actual[n]-initial_actual[n] for n in self.joints}
        state_time=state.get('time_s')
        age=None if state_time is None or self.native_moment_arms_time_s is None else state_time-self.native_moment_arms_time_s
        self.last_allocation.update({
            'requested_delta_torque_nm':dict(requested),
            'estimated_commanded_delta_torque_nm':dict(zip(self.joints,map(float,achieved-equilibrium_torque))),
            'estimated_commanded_total_torque_nm':dict(zip(self.joints,map(float,achieved))),
            'actual_muscle_torque_nm':actual,
            'initial_reference_muscle_torque_nm':initial_actual,
            'actual_delta_torque_nm':actual_delta,
            'actual_activation':{n:state['muscles'][n].get('activation') for n in self.muscles},
            'commanded_excitations':dict(commands),
            'excitations_at_lower_bound':[n for n,u in commands.items() if u<=1e-10],
            'excitations_at_upper_bound':[n for n,u in commands.items() if u>=1-1e-10],
            'torque_observation_clock':{'state_time_s':state_time,'reference_time_s':self.reference_time_s,'fitted_moment_arm_time_s':state_time,'native_moment_arm_time_s':self.native_moment_arms_time_s,'native_moment_arm_age_s':age},
            'torque_observation_coverage':{'missing_current_tendon_force':missing_forces,'missing_current_moment_arms':missing_geometry,'missing_reference_tendon_force':initial_missing_forces,'missing_reference_moment_arms':initial_missing_geometry,'excluded_muscles_without_moment_arms':sorted(set(self.muscles)-set(self.controlled))},
            'torque_observation_basis':'Measured native tendon force times current fitted derivative or last supplied native moment arm. Muscle contributions only; excludes joint passive forces, gravity and contact. Actual is the observed pre-command state, not a prediction of the issued command.'})
        return commands

    def identity(self):
        return {'schema':'ihm.engineering-joint-posture.v1', 'kp':self.kp, 'kd':self.kd,
                'pelvis_gain':self.pelvis_gain,'pelvis_damping':self.pelvis_damping,
                'baseline':self.baseline,'regularization':self.regularization,
                'com_lateral_position_gain_nm_per_m':self.com_lateral_position_gain,'com_lateral_velocity_gain_nm_s_per_m':self.com_lateral_velocity_gain,
                'com_target_z_m':self.com_target_z_m,'initial_com_z_ground_m':self.initial_com_z,
                'com_position_gain_nm_per_m':self.com_position_gain,'com_velocity_gain_nm_s_per_m':self.com_velocity_gain,
                'initial_com_x_ground_m':self.initial_com_x,'com_target_x_m':self.com_target_x_m,
                'initial_support_center_m':self.initial_support_center,'initial_support_center_ground_m':self.initial_support_center,
                'equilibrium_excitations':self.equilibrium_excitations,
                'targets_rad':self.targets, 'feedforward_torques_nm':self.feedforward,
                'fitted_muscles':len(self.fitted), 'native_moment_arm_muscles':sorted(self.native_moment_arms),
                'native_moment_arm_refresh':'caller refreshes with body motion; last geometry held between updates',
                'unmodeled_muscles':sorted(set(self.muscles)-set(self.controlled)),
                'brain_trained':False,'external_support':False,'prescribed_motion':False,
                'walking_demonstrated':False,'strength_model':'isometric Gaussian force-length approximation; no pennation, force-velocity, passive force or activation lag'}
