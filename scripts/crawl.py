#!/usr/bin/env python3
"""Prone locomotion on the native plant: settle the body on its belly, then crawl.

Why prone.  Bipedal gait on this model is a balance problem that the deflated
stance LQR only barely solves: one genuine step, then a fall.  Prone locomotion
is not a balance problem at all -- the support polygon is the whole ventral
surface, the centre of mass is already inside it, and the body cannot fall
because it is already down.  That makes it the cheap developmental precursor the
programme actually asked for, and it isolates *propulsion* from *balance*.

What drives the body.  Muscle excitation in [0,1] on the 98 source muscles, plus
the source model's declared CoordinateActuator torque ports for the lumbar and
both arms.  Those ports are torque actuators, NOT muscles: the model has no
shoulder musculature to innervate, and nothing here may be reported as
muscle-driven.  No coordinate is prescribed, no external force is applied to any
body, and no motion constraint exists.

Contact.  ``environment='upright'`` gives the source foot contacts plus one
inertia-inscribed sphere per non-foot body against the floor, so a body lying on
its front is supported by its own trunk and limb proxies.  Those spheres are
engineering proxies from segment inertia, not an anatomical skin reconstruction,
and the report says so.
"""
from __future__ import annotations

import argparse, json, math, os, shutil, sys, time, traceback, uuid
from pathlib import Path

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream  # noqa: E402

BUNDLE = ROOT / 'data/models/engineering_stance_v1'
REGISTRATION = 'data/models/engineering_stance_v1/registration.json'
TARGET_MASS_KG = 77.6122029
DT = 0.01
WORK = ROOT / 'data/derived/crawl-work'

# The five declared torque ports per arm, and the three lumbar ones.
ARM_PORT = {'arm_flex': 'shoulder_flex', 'arm_add': 'shoulder_add',
            'arm_rot': 'shoulder_rot', 'elbow_flex': 'elbow_flex',
            'pro_sup': 'pro_sup'}
LUMBAR_PORT = {'lumbar_extension': 'lumbar_ext', 'lumbar_bending': 'lumbar_bend',
               'lumbar_rotation': 'lumbar_rot'}

# Every port gets the SAME shaping, applied to its own effective inertia:
# kp = I w^2 and kd = 2 z I w, with one w and one z for the whole body.  The
# numbers below were measured, not chosen.
#
# The failure being avoided is at BOTH ends.  Too soft and the port lets its
# coordinate run: the first crawl attempt put the right elbow at 10.6 rad/s
# while its own target moved at 4.5 rad/s, which dragged arm26_BIClong_r to
# 38.1 optimal fibre lengths per second against the muscle's maximum of 10 --
# 3.8x outside the force-velocity curve's domain, where the equilibrium muscle
# inversion is ill-conditioned and the Simbody error controller collapses its
# step.  Too stiff and the 100 Hz command update itself is the instability: a
# sweep over (w, z) at 0.6 s of settling counted torque commands driven to their
# +/-1 ceiling, and z = 2.1 chattered the arm chain to a peak coordinate speed of
# 47.9 rad/s at w = 14 and 110.3 at w = 20, saturating 392 and 500 commands out
# of 780, against 1-8 for every z <= 1.0.  Over-damping a zero-order-held port is
# not a safe direction; it is the other cliff.
#
# w = 10 rad/s, z = 0.9 sits between them, and puts kd*dt/I at 0.18 -- an order
# of magnitude inside the kd*dt/I < 2 stability limit for a held damper.
# Effective inertias are read off the source segment masses and inertias.
PORT_INERTIA_KG_M2 = {
    'arm_flex': 0.35, 'arm_add': 0.35,     # whole arm about the shoulder
    'arm_rot': 0.02,                       # humerus about its own long axis
    'elbow_flex': 0.05,                    # forearm + hand about the elbow
    'pro_sup': 0.0015,                     # radius + hand about the forearm axis
    'lumbar_extension': 1.00, 'lumbar_bending': 1.00, 'lumbar_rotation': 1.00,
}
PORT_NATURAL_FREQUENCY_RAD_S = 10.0
PORT_DAMPING_RATIO = 0.9


def _gains(inertia):
    w, z = PORT_NATURAL_FREQUENCY_RAD_S, PORT_DAMPING_RATIO
    return inertia * w * w, 2 * z * inertia * w


ARM_GAINS = {k: _gains(v) for k, v in PORT_INERTIA_KG_M2.items() if k in ARM_PORT}
LUMBAR_GAINS = {k: _gains(v) for k, v in PORT_INERTIA_KG_M2.items() if k in LUMBAR_PORT}
OPTIMAL_FORCE_NM = 50.0

MUSCLE_GROUPS = {
    'hipflex':  ('iliacus', 'psoas', 'recfem', 'sart', 'tfl'),
    'hipext':   ('glmax1', 'glmax2', 'glmax3', 'bflh', 'semimem', 'semiten', 'addmagIsch'),
    'kneeext':  ('vasint', 'vaslat', 'vasmed', 'recfem'),
    'kneeflex': ('bfsh', 'bflh', 'semimem', 'semiten', 'grac', 'gaslat', 'gasmed'),
    'plantar':  ('soleus', 'gaslat', 'gasmed', 'tibpost', 'fhl', 'fdl', 'perlong', 'perbrev'),
    'dorsi':    ('tibant', 'edl', 'ehl'),
    'abduct':   ('glmed1', 'glmed2', 'glmed3', 'glmin1', 'glmin2', 'glmin3', 'tfl'),
    'adduct':   ('addbrev', 'addlong', 'addmagDist', 'addmagMid', 'addmagProx', 'grac'),
}


def pd_commands(state, targets, gains, ports):
    out = {}
    coords = state['coordinates']
    for coordinate, port in ports.items():
        kp, kd = gains[coordinate]
        c = coords[coordinate]
        torque = kp * (targets[coordinate] - c['value']) - kd * c['speed']
        out[port] = float(np.clip(torque / OPTIMAL_FORCE_NM, -1.0, 1.0))
    return out


def limb_pd(state, targets):
    """PD hold/drive on every declared torque port, arms per side plus lumbar."""
    out = {}
    coords = state['coordinates']
    for stem, port in ARM_PORT.items():
        kp, kd = ARM_GAINS[stem]
        for side in 'rl':
            name = stem + '_' + side
            c = coords[name]
            torque = kp * (targets[name] - c['value']) - kd * c['speed']
            out[port + '_' + side] = float(np.clip(torque / OPTIMAL_FORCE_NM, -1.0, 1.0))
    for coordinate, port in LUMBAR_PORT.items():
        kp, kd = LUMBAR_GAINS[coordinate]
        c = coords[coordinate]
        torque = kp * (targets[coordinate] - c['value']) - kd * c['speed']
        out[port] = float(np.clip(torque / OPTIMAL_FORCE_NM, -1.0, 1.0))
    return out


def open_stream(out_dir, pose):
    return NativeMechanicalStream(ROOT, out_dir, environment='upright',
                                  target_mass_kg=TARGET_MASS_KG, initial_pose=pose,
                                  augmented_registration=REGISTRATION)


def contact_summary(state):
    rows = []
    for c in state['contacts']:
        f = math.sqrt(sum(v * v for v in c['force_n']))
        if f > 1.0:
            rows.append({'name': c['name'], 'force_n': f, 'centre_y_m': c['center_m'][1]})
    rows.sort(key=lambda r: -r['force_n'])
    return rows


def probe(pose, seconds, out_dir, excitation=0.02):
    """Drop the body in a candidate pose and report where it comes to rest."""
    shutil.rmtree(out_dir, ignore_errors=True)
    native = open_stream(out_dir, pose)
    try:
        state = native.snapshot()
        muscles = sorted(state['muscles'])
        targets = {k: float(pose.get(k, 0.0)) for k in
                   [s + '_' + t for s in ARM_PORT for t in 'rl'] + list(LUMBAR_PORT)}
        track = []
        t0 = time.time()
        for i in range(int(round(seconds / DT))):
            state = native.advance(DT, actuation={m: excitation for m in muscles},
                                   coordinate_actuation=limb_pd(state, targets))
            if i % 10 == 0:
                c = state['coordinates']
                track.append({'time_s': state['time_s'],
                              'pelvis_ty': c['pelvis_ty']['value'],
                              'pelvis_tx': c['pelvis_tx']['value'],
                              'pelvis_tilt': c['pelvis_tilt']['value'],
                              'pelvis_list': c['pelvis_list']['value'],
                              'contact_force_n': math.sqrt(sum(v * v for v in state['contact_force_n'])),
                              'kinetic_energy_j': state['kinetic_energy_j']})
        return {'track': track, 'wall_s': time.time() - t0,
                'final_contacts': contact_summary(state),
                'final_coordinates': {k: v['value'] for k, v in state['coordinates'].items()},
                'bodies_on_floor': sorted({c['name'] for c in contact_summary(state)})}
    finally:
        native.close()


# ------------------------------------------------------------------- controller
# The prone pose the body settles into from a face-down drop.  Only the pelvis is
# set: every other coordinate is left at the model default, because ``assemble``
# rejects a pose that violates the source coupler constraints and a hand-written
# whole-body pose does.
PRONE_POSE = {'pelvis_tilt': -1.5708, 'pelvis_ty': 0.25}

BOUNDS = {
    'period':      (0.60, 2.60),   # one full crawl cycle, s
    'settle_s':    (0.20, 0.80),   # let the drop damp out before driving
    'ramp_s':      (0.15, 0.80),   # fade the pattern in; a step change in 98
                                   # excitations at settle_s whipped edl_r to
                                   # 16.7 optimal fibre lengths per second
    'e_base':      (0.01, 0.25),   # tonic excitation on every muscle
    'a_hipflex':   (0.00, 1.00),   # draw the knee up under the body
    'a_hipext':    (0.00, 1.00),   # drive it back down and push
    'a_kneeflex':  (0.00, 1.00),
    'a_kneeext':   (0.00, 1.00),
    'a_plantar':   (0.00, 1.00),   # toe push-off against the floor
    'a_dorsi':     (0.00, 0.80),   # lift the toe clear on recovery
    'a_abduct':    (0.00, 0.80),   # frog-kick the knee out to the side
    'a_adduct':    (0.00, 0.80),
    'psi_knee':    (-3.14, 3.14),  # knee phase relative to the hip
    'psi_ankle':   (-3.14, 3.14),
    'arm_flex_mid': (-0.40, 1.60), # shoulder hold about which the arm sweeps, rad
    'arm_flex_amp': (0.00, 1.20),  # reach forward / pull back
    'arm_add_mid': (-0.60, 0.60),
    'arm_add_amp': (0.00, 0.60),
    'elbow_mid':   (0.10, 2.00),
    'elbow_amp':   (0.00, 1.00),
    'psi_arm':     (-3.14, 3.14),  # arm phase relative to the SAME-side leg
    'psi_elbow':   (-3.14, 3.14),
    'lumbar_bend_amp': (0.00, 0.50),   # lateral trunk undulation, rad
    'lumbar_ext_mid':  (-0.40, 0.40),
    'lumbar_ext_amp':  (0.00, 0.40),
    'psi_lumbar':  (-3.14, 3.14),
}
ORDER = tuple(BOUNDS)

SEED = {
    'period': 1.40, 'settle_s': 0.40, 'ramp_s': 0.40, 'e_base': 0.05,
    # Timid on the ankle: the seed's original 0.40 plantar against 0.20 dorsi
    # whipped edl_l to 15.1 optimal fibre lengths per second within 0.85 s, on an
    # unloaded foot that has nothing to push against while the body is prone.
    'a_hipflex': 0.35, 'a_hipext': 0.55, 'a_kneeflex': 0.30, 'a_kneeext': 0.30,
    'a_plantar': 0.15, 'a_dorsi': 0.05, 'a_abduct': 0.30, 'a_adduct': 0.20,
    'psi_knee': 0.60, 'psi_ankle': -1.20,
    'arm_flex_mid': 0.50, 'arm_flex_amp': 0.60,
    'arm_add_mid': 0.00, 'arm_add_amp': 0.20,
    'elbow_mid': 0.90, 'elbow_amp': 0.40,
    'psi_arm': 3.14, 'psi_elbow': 0.00,
    'lumbar_bend_amp': 0.15, 'lumbar_ext_mid': 0.00, 'lumbar_ext_amp': 0.10,
    'psi_lumbar': 0.00,
}


def clamp(p):
    return {k: float(np.clip(p[k], *BOUNDS[k])) for k in ORDER}


class CrawlPattern:
    """One open-loop crawl cycle: muscle excitation for the legs, torque ports
    for the arms and lumbar.  Contralateral by construction -- the right arm is
    phase-locked to the left leg through ``psi_arm``.

    Nothing here is feedback.  A prone body cannot fall, so the whole point of
    the balance regulator that bipedal gait needs is absent, and an open-loop
    pattern is the honest minimum: whatever travel it produces is produced by
    the muscles and the floor, not by a controller correcting the plant.
    """

    def __init__(self, params, muscle_names):
        self.q = clamp(params)
        self.muscle_names = list(muscle_names)
        self.index = {n: i for i, n in enumerate(self.muscle_names)}
        self.group = {}
        for group, stems in MUSCLE_GROUPS.items():
            for side in 'rl':
                self.group[(group, side)] = [self.index[s + '_' + side]
                                             for s in stems if s + '_' + side in self.index]

    def phase(self, t_s, side):
        q = self.q
        drive = max(0.0, t_s - q['settle_s'])
        return 2 * math.pi * drive / q['period'] + (0.0 if side == 'r' else math.pi)

    def ramp(self, t_s):
        q = self.q
        return float(np.clip((t_s - q['settle_s']) / max(q['ramp_s'], DT), 0.0, 1.0))

    def excitation(self, t_s):
        q = self.q
        u = np.full(len(self.muscle_names), q['e_base'])
        gain = self.ramp(t_s)
        if gain <= 0.0:
            return u
        for side in 'rl':
            phi = self.phase(t_s, side)
            half = lambda x: max(0.0, math.sin(x))            # noqa: E731
            other = lambda x: max(0.0, -math.sin(x))          # noqa: E731
            add = {
                ('hipflex', side):  q['a_hipflex'] * half(phi),
                ('hipext', side):   q['a_hipext'] * other(phi),
                ('kneeflex', side): q['a_kneeflex'] * half(phi + q['psi_knee']),
                ('kneeext', side):  q['a_kneeext'] * other(phi + q['psi_knee']),
                ('dorsi', side):    q['a_dorsi'] * half(phi + q['psi_ankle']),
                ('plantar', side):  q['a_plantar'] * other(phi + q['psi_ankle']),
                ('abduct', side):   q['a_abduct'] * half(phi),
                ('adduct', side):   q['a_adduct'] * other(phi),
            }
            for key, value in add.items():
                u[self.group[key]] += gain * value
        return np.clip(u, 0.0, 1.0)

    def limb_targets(self, t_s):
        q = self.q
        gain = self.ramp(t_s)
        targets = {}
        for side in 'rl':
            # contralateral: the right arm follows the LEFT leg's phase
            phi = self.phase(t_s, 'l' if side == 'r' else 'r') + q['psi_arm']
            targets['arm_flex_' + side] = q['arm_flex_mid'] + q['arm_flex_amp'] * math.sin(phi)
            targets['arm_add_' + side] = ((q['arm_add_mid'] + q['arm_add_amp'] * math.sin(phi))
                                          * (1.0 if side == 'r' else -1.0))
            targets['arm_rot_' + side] = 0.0
            targets['elbow_flex_' + side] = float(np.clip(
                q['elbow_mid'] + q['elbow_amp'] * math.sin(phi + q['psi_elbow']), 0.05, 2.55))
            targets['pro_sup_' + side] = 0.3
        lum = self.phase(t_s, 'r') + q['psi_lumbar']
        targets['lumbar_bending'] = q['lumbar_bend_amp'] * math.sin(lum)
        targets['lumbar_extension'] = q['lumbar_ext_mid'] + q['lumbar_ext_amp'] * math.sin(lum)
        targets['lumbar_rotation'] = 0.0
        # fade from the settled hold pose into the pattern
        for k in targets:
            rest = 0.3 if k.startswith('pro_sup') else 0.0
            targets[k] = rest + gain * (targets[k] - rest)
        return targets


# ---------------------------------------------------------------------- rollout
def plant_summary(state, command, muscle_names):
    """What the plant looked like entering an advance the integrator struggled on."""
    coords = state['coordinates']
    muscles = state['muscles']
    speeds = sorted(((abs(v['speed']), k) for k, v in coords.items()), reverse=True)[:10]
    fibres = sorted(((abs(m['fiber_velocity_m_s']) / max(m['optimal_fiber_length_m'], 1e-9), n)
                     for n, m in muscles.items()), reverse=True)[:10]
    short = sorted(((m['fiber_length_m'] / max(m['optimal_fiber_length_m'], 1e-9), n)
                    for n, m in muscles.items()))[:10]
    contacts = sorted(({'name': c['name'],
                        'force_n': math.sqrt(sum(v * v for v in c['force_n'])),
                        'penetration_proxy_m': c['radius_m'] - c['center_m'][1]}
                       for c in state['contacts']), key=lambda r: -r['force_n'])[:10]
    return {
        'fastest_coordinates': [{'coordinate': k, 'abs_speed': v} for v, k in speeds],
        'fastest_fibres_optimal_lengths_per_s': [{'muscle': n, 'normalized_speed': v}
                                                 for v, n in fibres],
        'shortest_fibres_normalized': [{'muscle': n, 'normalized_length': v} for v, n in short],
        'loaded_contacts': contacts,
        'commanded_at_ceiling': int(sum(1 for v in command if v >= 0.999)),
        'coordinate_actuator_commands': {k: v['command']
                                         for k, v in state['coordinate_actuators'].items()},
    }


# Both Thelen2003 and Millard2012 muscles in this model carry the default
# maximum contraction velocity of 10 optimal fibre lengths per second.  Outside
# that the force-velocity curve is being extrapolated, the equilibrium muscle's
# inversion is ill-conditioned, and the error-controlled integrator collapses its
# step: every native stall measured on this plant, in gait and in prone, has this
# signature.  Stopping here converts an unbounded wall-clock stall into a fast
# termination that names the muscle -- the plant has left the muscle model's
# validity domain, and a result computed past this point is not about the model
# anyone declared.
MAXIMUM_CONTRACTION_VELOCITY_OFL_S = 10.0
# The guard is set at 1.5x the model's maximum, not at it.  Measured: the plant
# crosses 10 briefly and harmlessly on an ankle whip (11.8 ofl/s on edl_r with
# every advance still costing 0.35 s or less), while every stall observed on this
# plant was already past 24 by the time the first slow advance appeared.  A guard
# AT the maximum ends almost every run on a transient and measures nothing; a
# guard well past the observed stall threshold measures nothing either.  Because
# 10 is nonetheless where the force-velocity curve stops being defined, every
# report also carries how long the run spent outside it, so a run that ends
# normally still says how far outside the declared muscle model it went.
FIBRE_VELOCITY_GUARD_OFL_S = 15.0


def maximum_normalized_fibre_velocity(state):
    worst, culprit = 0.0, None
    for name, m in state['muscles'].items():
        v = abs(m['fiber_velocity_m_s']) / max(m['optimal_fiber_length_m'], 1e-9)
        if v > worst:
            worst, culprit = v, name
    return worst, culprit


def diverged(state):
    coords = state['coordinates']
    if coords['pelvis_ty']['value'] > 1.20:
        return 'pelvis_above_1.20m'
    for name, c in coords.items():
        if abs(c['speed']) > 25.0:
            return 'joint_speed_beyond_25_' + name
    worst, culprit = maximum_normalized_fibre_velocity(state)
    if worst > FIBRE_VELOCITY_GUARD_OFL_S:
        return 'fibre_velocity_beyond_%.0f_ofl_s_%s' % (FIBRE_VELOCITY_GUARD_OFL_S, culprit)
    return None


def rollout(params, muscle_names, horizon_s=8.0, record=False, work_dir=None,
            wall_budget_s=600.0, trace_path=None, slow_s=2.0):
    params = clamp(params)
    out = Path(work_dir or (WORK / ('run-' + uuid.uuid4().hex)))
    shutil.rmtree(out, ignore_errors=True)
    pattern = CrawlPattern(params, muscle_names)
    frames, native = [], None
    failure, stop_reason = None, 'horizon'
    tx0 = tz0 = None
    tx_last = tz_last = 0.0
    tx_max = -1e9
    lateral = 0.0
    t_end = 0.0
    started = time.time()
    advance_wall_max = 0.0
    fibre_worst, fibre_culprit = 0.0, None
    outside_domain_s = 0.0
    trace = None if trace_path is None else Path(trace_path).open('w')
    try:
        native = open_stream(out, PRONE_POSE)
        state = native.snapshot()
        for _ in range(int(round(horizon_s / DT))):
            t_s = state['time_s']
            coords = state['coordinates']
            tx, tz = coords['pelvis_tx']['value'], coords['pelvis_tz']['value']
            if tx0 is None:
                tx0, tz0 = tx, tz
            tx_last, tz_last = tx - tx0, tz - tz0
            tx_max = max(tx_max, tx_last)
            lateral = max(lateral, abs(tz_last))
            worst, culprit = maximum_normalized_fibre_velocity(state)
            if worst > fibre_worst:
                fibre_worst, fibre_culprit = worst, culprit
            if worst > MAXIMUM_CONTRACTION_VELOCITY_OFL_S:
                outside_domain_s += DT
            u = pattern.excitation(t_s)
            targets = pattern.limb_targets(t_s)
            if record:
                frames.append({
                    'time_s': t_s,
                    'joints': {k: {'value': v['value'], 'speed': v['speed'], 'unit': v['unit']}
                               for k, v in coords.items()},
                    'contacts': [{'name': c['name'],
                                  'force_n': math.sqrt(sum(v * v for v in c['force_n'])),
                                  'center_m': c['center_m']}
                                 for c in state['contacts']
                                 if math.sqrt(sum(v * v for v in c['force_n'])) > 1.0],
                    'motor_excitations': {n: float(v) for n, v in zip(muscle_names, u)},
                    'coordinate_actuator_commands': {
                        k: v['command'] for k, v in state['coordinate_actuators'].items()},
                })
            reason = diverged(state)
            if reason is not None:
                stop_reason, t_end = reason, t_s
                break
            if time.time() - started > wall_budget_s:
                stop_reason, t_end = 'wall_budget', t_s
                break
            step_started = time.time()
            entering = state
            try:
                state = native.advance(DT, actuation=dict(zip(muscle_names, map(float, u))),
                                       coordinate_actuation=limb_pd(state, targets))
            except BaseException as exc:
                if trace is not None:
                    trace.write(json.dumps({'time_s': t_s, 'wall_s': time.time() - step_started,
                                            'failed': True,
                                            'exception': type(exc).__name__ + ': ' + str(exc)[:200],
                                            'entering_state': plant_summary(entering, u, muscle_names)}) + '\n')
                    trace.flush()
                raise
            wall = time.time() - step_started
            advance_wall_max = max(advance_wall_max, wall)
            if trace is not None:
                row = {'time_s': t_s, 'wall_s': wall, 'failed': False}
                row['max_normalized_fibre_velocity'] = worst
                if wall >= slow_s:
                    row['entering_state'] = plant_summary(entering, u, muscle_names)
                trace.write(json.dumps(row) + '\n')
                trace.flush()
            t_end = state['time_s']
    except Exception as exc:
        failure = type(exc).__name__ + ': ' + str(exc)[:200]
        stop_reason = 'native_failure'
    finally:
        if trace is not None:
            trace.close()
        if native is not None:
            try:
                native.close()
            except Exception:
                pass
        if not record:
            shutil.rmtree(out, ignore_errors=True)
    report = {
        'params': params,
        'horizon_s': horizon_s,
        'simulated_s': t_end,
        'wall_s': time.time() - started,
        'slowest_advance_wall_s': advance_wall_max,
        'peak_normalized_fibre_velocity_ofl_s': fibre_worst,
        'peak_normalized_fibre_velocity_muscle': fibre_culprit,
        'muscle_maximum_contraction_velocity_ofl_s': MAXIMUM_CONTRACTION_VELOCITY_OFL_S,
        'seconds_outside_muscle_force_velocity_domain': outside_domain_s,
        'fraction_outside_muscle_force_velocity_domain': (
            None if t_end <= 0 else outside_domain_s / t_end),
        'completed_horizon': stop_reason == 'horizon',
        'stop_reason': stop_reason,
        'failure': failure,
        'pelvis_forward_travel_m': None if tx0 is None else tx_last,
        'pelvis_peak_forward_travel_m': None if tx0 is None else tx_max,
        'pelvis_lateral_excursion_m': None if tx0 is None else lateral,
        'mean_forward_velocity_m_s': None if tx0 is None or t_end <= 0 else tx_last / t_end,
        'posture': 'prone',
        'drive': ('muscle excitation on 98 source muscles plus the declared '
                  'CoordinateActuator torque ports for the lumbar and both arms; '
                  'no prescribed coordinate, no external force, no motion constraint'),
        'contact_basis': ('source foot contacts plus one inertia-inscribed sphere per '
                          'non-foot body against the floor: engineering proxies, not an '
                          'anatomical skin surface'),
    }
    if record:
        report['output_dir'] = str(out.relative_to(ROOT))
    return report, frames


def score(report):
    if report['pelvis_forward_travel_m'] is None:
        return -1e6
    s = 100.0 * report['pelvis_forward_travel_m']
    s -= 20.0 * max(0.0, report['pelvis_lateral_excursion_m'] - 0.10)
    s += 2.0 * report['simulated_s']
    s -= 40.0 * max(0.0, report['horizon_s'] - report['simulated_s'])
    return s


# ----------------------------------------------------------------------- search
_MUSCLES = None


def muscle_names():
    """Read the muscle ordering from the plant itself, once per process."""
    global _MUSCLES
    if _MUSCLES is None:
        out = WORK / ('names-' + uuid.uuid4().hex)
        shutil.rmtree(out, ignore_errors=True)
        native = open_stream(out, PRONE_POSE)
        try:
            _MUSCLES = sorted(native.snapshot()['muscles'])
        finally:
            native.close()
            shutil.rmtree(out, ignore_errors=True)
    return _MUSCLES


def _worker(job):
    index, params, horizon = job
    try:
        report, _ = rollout(params, muscle_names(), horizon_s=horizon)
    except Exception:
        return index, {'params': clamp(params), 'stop_reason': 'worker_exception',
                       'failure': traceback.format_exc()[-300:], 'simulated_s': 0.0,
                       'horizon_s': horizon, 'completed_horizon': False,
                       'pelvis_forward_travel_m': None}, -1e6
    return index, report, score(report)


def _kill_orphan_engines():
    import signal, subprocess
    try:
        listing = subprocess.run(['ps', '-eo', 'pid,args'], capture_output=True, text=True).stdout
    except Exception:
        return
    for line in listing.splitlines():
        if 'native_mechanical_stream' in line and 'crawl-work' in line:
            try:
                os.kill(int(line.split()[0]), signal.SIGKILL)
            except Exception:
                pass


def search(generations, population, horizon, workers, sigma0, out_path, seed=0, initial=None,
           generation_timeout_s=1200.0):
    import multiprocessing as mp
    WORK.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    span = np.array([BOUNDS[k][1] - BOUNDS[k][0] for k in ORDER])
    lo = np.array([BOUNDS[k][0] for k in ORDER])
    hi = np.array([BOUNDS[k][1] for k in ORDER])
    incumbent = np.array([(initial or SEED)[k] for k in ORDER], float)
    best_report, best_score, history = None, -1e9, []
    sigma = sigma0
    started = time.time()
    ctx = mp.get_context('spawn')
    pool = ctx.Pool(workers)
    try:
        for generation in range(generations):
            jobs = []
            for i in range(population):
                if generation == 0 and i == 0:
                    candidate = incumbent.copy()
                else:
                    candidate = incumbent + rng.normal(0, sigma, len(ORDER)) * span
                jobs.append((i, dict(zip(ORDER, np.clip(candidate, lo, hi))), horizon))
            try:
                results = pool.map_async(_worker, jobs).get(timeout=generation_timeout_s)
            except Exception as exc:
                print(json.dumps({'generation': generation, 'pool_reset': str(exc)[:120]}),
                      flush=True)
                pool.terminate()
                pool.join()
                _kill_orphan_engines()
                pool = ctx.Pool(workers)
                continue
            results.sort(key=lambda r: -r[2])
            top = results[0]
            if top[2] > best_score:
                best_score, best_report = top[2], top[1]
                incumbent = np.array([best_report['params'][k] for k in ORDER])
                sigma = sigma0
            else:
                sigma = max(sigma * 0.9, sigma0 * 0.4)
            row = {'generation': generation, 'wall_s': time.time() - started, 'sigma': sigma,
                   'generation_best_score': top[2], 'best_score': best_score,
                   'best_travel_m': best_report['pelvis_forward_travel_m'],
                   'best_mean_velocity_m_s': best_report['mean_forward_velocity_m_s'],
                   'best_simulated_s': best_report['simulated_s'],
                   'best_stop': best_report['stop_reason']}
            history.append(row)
            print(json.dumps(row), flush=True)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(json.dumps(
                {'best_score': best_score, 'best_report': best_report, 'history': history},
                indent=2) + '\n')
    finally:
        pool.terminate()
        pool.join()
        _kill_orphan_engines()
    return best_report, best_score, history


def write_best(params, horizon, destination):
    destination.mkdir(parents=True, exist_ok=True)
    work = WORK / ('best-' + uuid.uuid4().hex)
    report, frames = rollout(params, muscle_names(), horizon_s=horizon, record=True,
                             work_dir=work, wall_budget_s=1e9)
    trajectory = {
        'schema': 'ihm.crawl-trajectory.v1',
        'dt_s': DT,
        'target_mass_kg': TARGET_MASS_KG,
        'muscle_names': muscle_names(),
        'initial_pose': PRONE_POSE,
        'basis': ('Prone locomotion.  Muscle excitation on the 98 source muscles plus the '
                  'source model\'s declared CoordinateActuator torque ports for the lumbar '
                  'and both arms.  Those ports are torque actuators, NOT muscles: the model '
                  'carries no shoulder musculature, so this motion is not wholly '
                  'muscle-driven and must not be reported as such.  No prescribed '
                  'coordinate, no external force on any body, no motion constraint; joint '
                  'values are integrated native OpenSim/Simbody output.'),
        'frames': frames,
    }
    (destination / 'trajectory.json').write_text(json.dumps(trajectory) + '\n')
    report['frames'] = len(frames)
    report['trajectory_path'] = str((destination / 'trajectory.json').relative_to(ROOT))
    (destination / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    shutil.rmtree(work, ignore_errors=True)
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', choices=('probe', 'single', 'search', 'best'), default='single')
    ap.add_argument('--horizon', type=float, default=8.0)
    ap.add_argument('--generations', type=int, default=12)
    ap.add_argument('--population', type=int, default=10)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--sigma', type=float, default=0.10)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--params')
    ap.add_argument('--search-out', default='data/derived/crawl-search/search.json')
    ap.add_argument('--out', default='data/derived/crawl-best')
    ap.add_argument('--trace', help='write a per-advance wall-time trace here')
    a = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    params = dict(SEED)
    if a.params:
        loaded = json.loads(Path(a.params).read_text())
        if 'best_report' in loaded:
            loaded = loaded['best_report']['params']
        elif 'params' in loaded:
            loaded = loaded['params']
        params.update(loaded)
    if a.mode == 'probe':
        out = WORK / 'probe'
        result = probe(PRONE_POSE, a.horizon, out)
        shutil.rmtree(out, ignore_errors=True)
        print(json.dumps(result, indent=1))
        return
    if a.mode == 'search':
        best, best_score, _ = search(a.generations, a.population, a.horizon, a.workers,
                                     a.sigma, ROOT / a.search_out, seed=a.seed, initial=params)
        print(json.dumps({'best_score': best_score,
                          'travel_m': best['pelvis_forward_travel_m'],
                          'mean_velocity_m_s': best['mean_forward_velocity_m_s']}, indent=2))
        return
    if a.mode == 'best':
        report = write_best(params, a.horizon, ROOT / a.out)
    else:
        report, _ = rollout(params, muscle_names(), horizon_s=a.horizon, wall_budget_s=1e9,
                            trace_path=a.trace)
    print(json.dumps({k: v for k, v in report.items() if k != 'params'}, indent=2))


if __name__ == '__main__':
    main()
