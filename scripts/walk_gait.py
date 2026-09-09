#!/usr/bin/env python3
"""Contact-triggered finite-state gait controller over the engineered LQR stance policy.

Everything this script commands is muscle excitation in [0,1]. No coordinate is
prescribed, no external force is applied to the root, and no motion constraint
exists: the native OpenSim/Simbody plant integrates whatever the muscles do.

Design
------
* The engineered stance feedback ``u = u0 - K (x - x0)`` from
  ``data/models/engineering_stance_v1/linearization.npz`` is the balance core.
* ``K`` is redesigned with the ground-plane symmetry deflated out, exactly as
  ``scripts/design_native_deflated_lqr.py`` does.  ``docs/NATIVE_LANDING_GROUND_PLANE_SYMMETRY.md``
  section 3 shows that a gain with feedback on absolute ``pelvis_tx``/``pelvis_tz``/
  ``pelvis_rotation`` cannot walk, because walking is the act of increasing
  ``pelvis_tx`` and that error is not reducible by any muscle.  Deflation is a
  precondition, not a tuning choice.
* A finite state machine (settle -> double support -> swing -> reach -> ...)
  moves the LQR *reference* ``x0 -> x0 + dx(phase)`` on a handful of joint
  coordinates and one velocity (``pelvis_tx/speed``), and adds per-phase muscle
  group biases on top of the feedback command.  Phase exits are driven by
  measured per-foot contact load, not by a clock, except for guard timeouts.
* The reference offset is low-pass filtered so no phase entry is a step change.

Sign conventions were measured on this model with ``evaluate_static_pose``
(zero-time-advance queries), not assumed:
  +x anterior, +y up, +z to the subject's right;
  hip_flexion +   -> foot forward/up          (swing the leg forward)
  knee_angle +    -> flexion, foot back/up
  ankle_angle +   -> dorsiflexion, foot up    (so - is plantarflexion / push-off)
  hip_adduction + -> that foot toward midline, i.e. pelvis toward that foot
  pelvis_tilt +   -> posterior (backward) lean
"""
from __future__ import annotations

import argparse, json, math, os, shutil, sys, time, traceback, uuid
from pathlib import Path

# One BLAS thread per process: the search runs many workers and a threaded BLAS
# oversubscribes the machine badly enough to dominate the wall time.
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_v, '1')

import numpy as np
from scipy.linalg import expm, solve_discrete_are

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream  # noqa: E402

BUNDLE = ROOT / 'data/models/engineering_stance_v1'
REGISTRATION = 'data/models/engineering_stance_v1/registration.json'
TARGET_MASS_KG = 77.6122029
DT = 0.01
WORK = ROOT / 'data/derived/gait-work'

# ------------------------------------------------------------------- arm hold
# The source model declares CoordinateActuator torque ports for the lumbar and
# for both arms (shoulder flexion/adduction/rotation, elbow, forearm pronation)
# and, until the native engine grew a port for them, NO controller was connected,
# so their controls were identically zero for the whole run.  The arm chain was a
# passive rag doll on a walking body.  Measured entering the advance that stalls
# the engine: pro_sup_r at 11.0 rad/s -- a coordinate the source passive force set
# does not damp at all -- and arm26_BIClong_l at 3.05 m/s of fibre velocity
# against a 1.28 m/s maximum contraction velocity.  Outside the force-velocity
# curve's domain the Simbody error controller collapses its step, and one 10 ms
# advance costs minutes.
#
# These are torque actuators, not muscles.  They are declared scaffolding of the
# same status as the engine's inertia-inscribed fall-support spheres: the model
# has no shoulder musculature to innervate, so nothing here can be reported as
# muscle-driven motion.  Gains are sized to each segment's own inertia for an
# undamped natural frequency near 20 rad/s at critical damping, so a command
# updated at 100 Hz stays well inside the stable region: a single shared gain
# would be 380 rad/s on forearm pronation and would make the stiffness worse
# than the flailing.
ARM_PORT = {'arm_flex': 'shoulder_flex', 'arm_add': 'shoulder_add',
            'arm_rot': 'shoulder_rot', 'elbow_flex': 'elbow_flex',
            'pro_sup': 'pro_sup'}
ARM_HOLD_GAINS = {                     # coordinate stem: (kp N.m/rad, kd N.m.s/rad)
    'arm_flex':   (140.0, 14.0),       # whole arm about the shoulder, ~0.35 kg.m2
    'arm_add':    (140.0, 14.0),
    'arm_rot':    (8.0, 0.8),          # humerus about its own long axis, ~0.02
    'elbow_flex': (20.0, 2.0),         # forearm + hand about the elbow, ~0.05
    'pro_sup':    (0.6, 0.06),         # radius + hand about the forearm axis, ~0.0015
}
COORDINATE_ACTUATOR_OPTIMAL_FORCE_NM = 50.0   # every port in the source model


def arm_hold_targets(pose):
    """Hold each arm coordinate where the registered initial pose puts it."""
    targets = {}
    for stem in ARM_PORT:
        for side in 'rl':
            name = stem + '_' + side
            targets[name] = float(pose.get(name, 0.0))
    return targets


def arm_hold_commands(state, targets):
    """Critically damped PD on the declared arm torque ports, in [-1, 1]."""
    out = {}
    coords = state['coordinates']
    for stem, port in ARM_PORT.items():
        kp, kd = ARM_HOLD_GAINS[stem]
        for side in 'rl':
            name = stem + '_' + side
            c = coords[name]
            torque = kp * (targets[name] - c['value']) - kd * c['speed']
            out[port + '_' + side] = float(np.clip(
                torque / COORDINATE_ACTUATOR_OPTIMAL_FORCE_NM, -1.0, 1.0))
    return out


# ---------------------------------------------------------------- muscle groups
GROUPS = {
    'hipflex':  ('iliacus', 'psoas', 'recfem', 'sart', 'tfl'),
    'hipext':   ('glmax1', 'glmax2', 'glmax3', 'bflh', 'semimem', 'semiten', 'addmagIsch'),
    'kneeext':  ('vasint', 'vaslat', 'vasmed', 'recfem'),
    'kneeflex': ('bfsh', 'bflh', 'semimem', 'semiten', 'grac', 'gaslat', 'gasmed'),
    'plantar':  ('soleus', 'gaslat', 'gasmed', 'tibpost', 'fhl', 'fdl', 'perlong', 'perbrev'),
    'dorsi':    ('tibant', 'edl', 'ehl'),
    'abduct':   ('glmed1', 'glmed2', 'glmed3', 'glmin1', 'glmin2', 'glmin3', 'tfl'),
    'adduct':   ('addbrev', 'addlong', 'addmagDist', 'addmagMid', 'addmagProx', 'grac'),
}

# ------------------------------------------------------------------- parameters
# name: (low, high).  Everything the search is allowed to move.
BOUNDS = {
    'v_forward':    (0.00, 0.70),   # pelvis_tx/speed reference, m/s
    't_ds':         (0.05, 0.45),   # max double-support duration, s
    't_swing':      (0.15, 0.60),   # swing duration, s
    't_reach_max':  (0.15, 0.70),   # touchdown guard timeout, s
    'a_hip':        (0.15, 1.00),   # swing hip flexion reference offset, rad
    'a_knee':       (0.10, 1.30),   # swing knee flexion offset, rad
    'a_ankle':      (-0.10, 0.45),  # swing ankle dorsiflexion offset, rad
    'a_hip_land':   (0.00, 0.70),   # hip flexion held at touchdown, rad
    'a_knee_land':  (-0.10, 0.45),  # knee flexion held at touchdown, rad
    'a_ankle_land': (-0.10, 0.40),  # ankle dorsiflexion at touchdown, rad
    'a_stance_knee': (0.00, 0.60),  # THE HYPOTHESIS: stance knee flexion, rad
    'a_stance_ankle': (0.00, 0.35), # stance ankle dorsiflexion (lowers pelvis), rad
    'a_push':       (0.00, 0.50),   # stance-ankle plantarflexion late in swing, rad
    'a_push_ds':    (0.00, 0.60),   # trailing-ankle plantarflexion in double support, rad
    'a_add':        (0.00, 0.30),   # stance hip adduction (lateral weight shift), rad
    'a_swing_abd':  (0.00, 0.35),   # swing hip ABduction -> wider step, bigger base, rad
    'a_list':       (0.00, 0.25),   # pelvis roll toward the stance foot, rad
    'a_bend':       (-0.35, 0.35),  # lumbar bending, signed toward the stance side, rad
    'a_stance_hipext': (0.00, 1.00),# stance hip extension -> forward propulsion, rad
    'a_tilt':       (-0.20, 0.15),  # pelvis tilt reference offset, rad
    'b_hipflex':    (0.00, 1.00),   # swing hip flexor excitation bias
    'b_dorsi':      (0.00, 0.80),   # swing dorsiflexor bias
    'b_kneeflex':   (0.00, 0.80),   # swing knee flexor bias (foot clearance)
    'b_kneeext':    (0.00, 1.50),   # swing knee extensor bias during reach (extend to land)
    'b_push':       (0.00, 1.00),   # stance/trailing plantarflexor bias at push-off
    'b_vasti':      (-0.40, 0.40),  # stance knee extensor bias (negative = let it flex)
    'b_abduct':     (0.00, 1.00),   # stance hip abductor bias (frontal-plane support)
    'b_hipext':     (0.00, 0.80),   # stance hip extensor bias (propulsion)
    'gain':         (0.40, 1.60),   # LQR regulation scale about the standing pose
    'ff_gain':      (0.005, 0.15),  # weight on the +K dx feedforward synergy (never 0:
                                    # at 0 every a_* reference parameter is inert)
    'swing_reg':    (0.50, 1.00),   # how much the swing leg regulates itself to x0
                                    # (floored at 0.5: an ablation showed that
                                    # weakening the regulator, not the muscle
                                    # biases, is what makes the plant diverge)
    'swing_relax':  (0.00, 1.00),   # scale on the swing leg's own baseline u0 for its
                                    # antigravity muscles: at 1.0 the leg keeps its
                                    # stance excitation and never unloads the foot
    'stance_boost': (1.00, 3.00),   # scale on the STANCE leg's baseline u0 for the same
                                    # muscles.  u0 is a DOUBLE-support equilibrium, so
                                    # each leg's baseline carries about half the body.
                                    # In single support one leg must carry all of it.
    'cross_gate':   (0.00, 1.00),   # how much swing-leg error the other muscles see
    'tau':          (0.02, 0.30),   # reference low-pass time constant, s
    'load_off':     (0.02, 0.12),   # fraction of body weight below which a foot is "off"
    'load_on':      (0.05, 0.45),   # fraction above which touchdown is accepted
    'q_tx_speed':   (0.0, 40000.0), # extra LQR stage cost on forward velocity.  THE
                                    # transport knob: 0 gives 2 mm of travel, 1e4 gives
                                    # 120 mm, and it trades against stability.
}
ORDER = tuple(BOUNDS)

# The search only moves these.  With 34 knobs and a population of 14 a
# (1+lambda) search spends every generation on directions that do not matter;
# the rest are pinned at the values below and left alone.
FREE = ('v_forward', 't_ds', 't_swing', 't_reach_max', 'a_hip', 'a_knee',
        'a_stance_knee', 'a_add', 'a_swing_abd', 'b_hipflex', 'b_abduct', 'b_dorsi',
        'b_push', 'b_kneeflex', 'ff_gain', 'swing_reg', 'gain', 'load_on',
        'a_tilt', 'b_hipext', 'a_stance_ankle', 'swing_relax', 'stance_boost',
        # propulsion: trailing-limb push-off is the main source of forward COM
        # velocity, and freezing these left the search with no axis on the
        # quantity that was actually failing.
        'a_push', 'a_push_ds', 'a_stance_hipext',
        # landing geometry: what the second step has to land into.
        'a_hip_land', 'a_knee_land', 'a_ankle_land', 'a_ankle',
        'b_kneeext', 'b_vasti', 'tau', 'load_off', 'a_bend',
        'q_tx_speed')

# Deliberately timid: with max|K| = 156, large offsets and biases saturate the
# command vector and the plant flails within a second.  The search grows the
# amplitudes from here.
SEED_PARAMS = {
    'v_forward': 0.20, 't_ds': 0.20, 't_swing': 0.30, 't_reach_max': 0.30,
    'a_hip': 0.35, 'a_knee': 0.40, 'a_ankle': 0.12, 'a_hip_land': 0.18,
    'a_knee_land': 0.05, 'a_ankle_land': 0.05, 'a_stance_knee': 0.15,
    'a_stance_ankle': 0.08, 'a_push': 0.08, 'a_push_ds': 0.10, 'a_add': 0.08,
    'a_swing_abd': 0.06, 'a_list': 0.00, 'a_bend': 0.00, 'a_stance_hipext': 0.08,
    'a_tilt': -0.03,
    'b_hipflex': 0.12, 'b_dorsi': 0.08, 'b_kneeflex': 0.08, 'b_kneeext': 0.08,
    'b_push': 0.10, 'b_vasti': -0.04, 'b_abduct': 0.10, 'b_hipext': 0.06,
    'gain': 1.0, 'ff_gain': 0.05, 'swing_reg': 0.85, 'swing_relax': 0.40,
    'stance_boost': 1.6, 'cross_gate': 1.0, 'tau': 0.10,
    'load_off': 0.12, 'load_on': 0.16, 'q_tx_speed': 1500.0,
}


def clamp_params(p):
    return {k: float(np.clip(p[k], *BOUNDS[k])) for k in ORDER}


# ------------------------------------------------------------------ LQR policy
def deflation_basis(Ad, Bd, tolerance=1e-8, input_tolerance=1e-5, gap=100.0, maximum=8):
    """Unit-eigenvalue subspace of Ad that Bd cannot reach (see design_native_deflated_lqr.py)."""
    U, s, Vt = np.linalg.svd(Ad - np.eye(Ad.shape[0]))
    order = np.argsort(s)
    relative = s[order] / s[order[-1]]
    reach = np.linalg.norm(U[:, order].T @ Bd, axis=1) / np.linalg.norm(U[:, order].T, axis=1)
    k = 0
    while k < min(maximum, len(s) - 1) and relative[k] <= tolerance and reach[k] <= input_tolerance:
        k += 1
    if not k:
        raise ValueError('No unreachable unit-eigenvalue subspace found')
    if relative[k] < gap * max(relative[k - 1], np.finfo(float).eps):
        raise ValueError('Unit-eigenvalue subspace not separated from the spectrum')
    keep = order[:k]
    return Vt[keep].T.copy(), U[:, keep].T.copy()


GAIN_CACHE = WORK / 'deflated_gain.npz'


def design_gain(artifact=BUNDLE / 'linearization.npz', dt=DT, q_tx_speed=0.0):
    """Solve the deflated discrete LQR.  Expensive (255-state DARE); cached.

    ``q_tx_speed`` adds stage cost on ``pelvis_tx/speed`` and nothing else.  The
    shipped design is solved about a standing equilibrium whose forward velocity
    is exactly zero and weights that state at 0.1, which leaves the gain column
    on forward velocity at L2 3.12 against a median column of 1.10.  The
    controller therefore has almost no authority over the one quantity walking
    consists of: measured, a commanded 0.319 m/s produced 0.00127 m/s, short by
    a factor of 252.  Adding cost here is the difference between 2 mm and 120 mm
    of travel; it also costs stability, and that trade is the result.
    """
    with np.load(artifact, allow_pickle=False) as d:
        data = {k: d[k] for k in d.files}
    A, B = data['A'], data['B']
    n, m = B.shape
    transition = expm(np.block([[A, B], [np.zeros((m, n + m))]]) * dt)
    Ad, Bd = transition[:n, :n], transition[:n, n:]
    Q = data['Q'] * dt
    if q_tx_speed:
        names = [str(x) for x in data['state_names']]
        Q[names.index('/jointset/ground_pelvis/pelvis_tx/speed'),
          names.index('/jointset/ground_pelvis/pelvis_tx/speed')] += float(q_tx_speed) * dt
    R = np.diag(np.diag(data['R'] * dt)
                * (.05 / np.maximum(data['u0'] - data['minimum_activation'], .005)) ** 2)
    N, _ = deflation_basis(Ad, Bd)
    k = N.shape[1]
    Qr, _ = np.linalg.qr(np.hstack([N, np.eye(n)]))
    T = Qr[:, :n].copy()
    T[:, :k] = N
    Tinv = np.linalg.inv(T)
    S, C = Tinv[k:, :], T[:, k:]
    A22, B2, Q22 = S @ Ad @ C, S @ Bd, C.T @ Q @ C
    P = solve_discrete_are(A22, B2, Q22, R)
    K2 = np.linalg.solve(R + B2.T @ P @ B2, B2.T @ P @ A22)
    K = K2 @ S
    full = np.abs(np.linalg.eigvals(Ad - Bd @ K))
    return {
        'K': K,
        'deflated_dimension': np.array(k),
        'deflation_residual': np.array(float(np.abs(K @ N).max())),
        'deflated_states': np.array([str(data['state_names'][i])
                                     for i in np.argsort(-np.abs(N).sum(axis=1))[:4]]),
        'closed_loop_radius_excluding_deflated': np.array(float(np.sort(full)[-(k + 1)])),
        'gain_max': np.array(float(np.abs(K).max())),
        'q_tx_speed': np.array(float(q_tx_speed)),
        'forward_velocity_gain_norm': np.array(float(np.linalg.norm(
            K[:, [str(x) for x in data['state_names']].index(
                '/jointset/ground_pelvis/pelvis_tx/speed')]))),
    }


class Policy:
    """Deflated discrete LQR built from the engineered stance linearization."""

    def __init__(self, artifact=BUNDLE / 'linearization.npz', dt=DT, cache=GAIN_CACHE):
        with np.load(artifact, allow_pickle=False) as d:
            data = {k: d[k] for k in d.files}
        self.state_names = [str(x) for x in data['state_names']]
        self.muscle_names = [str(x) for x in data['muscle_names']]
        self.x0 = data['x0'].copy()
        self.u0 = data['u0'].copy()
        self.model_sha256 = str(data['model_sha256'].item())
        design = None
        if cache is not None and Path(cache).exists():
            try:
                with np.load(cache, allow_pickle=False) as d:
                    design = {k: d[k] for k in d.files}
            except Exception:
                design = None
        if design is None:
            design = design_gain(artifact, dt)
            if cache is not None:
                Path(cache).parent.mkdir(parents=True, exist_ok=True)
                tmp = Path(cache).with_name('tmp-' + uuid.uuid4().hex + '.npz')
                np.savez(tmp, **design)
                os.replace(tmp, cache)
        self._artifact = artifact
        self._dt = dt
        self._gain_cache = {0.0: design}
        self.K = design['K']
        self.deflated_dimension = int(design['deflated_dimension'])
        self.deflation_residual = float(design['deflation_residual'])
        self.deflated_states = [str(x) for x in design['deflated_states']]
        self.closed_loop_radius_excluding_deflated = float(design['closed_loop_radius_excluding_deflated'])
        self.gain_max = float(design['gain_max'])

        # index tables
        self.state_index = {name: i for i, name in enumerate(self.state_names)}
        self.muscle_index = {name: i for i, name in enumerate(self.muscle_names)}
        self.coord_value_index, self.coord_speed_index = {}, {}
        for name, i in self.state_index.items():
            parts = name.strip('/').split('/')
            if parts[0] == 'jointset':
                (self.coord_value_index if parts[-1] == 'value' else self.coord_speed_index)[parts[-2]] = i
        # per-side masks
        self.side_state_mask = {}
        for side in 'rl':
            mask = np.zeros(len(self.state_names))
            for name, i in self.state_index.items():
                parts = name.strip('/').split('/')
                token = parts[-2]
                if token.endswith('_' + side) and (parts[0] == 'forceset' or token.split('_')[0] in
                                                   ('hip', 'knee', 'ankle', 'mtp', 'subtalar')):
                    mask[i] = 1.0
            self.side_state_mask[side] = mask
        self.side_muscle_mask = {
            side: np.array([1.0 if n.endswith('_' + side) else 0.0 for n in self.muscle_names])
            for side in 'rl'}
        self.group_index = {}
        for group, stems in GROUPS.items():
            for side in 'rl':
                self.group_index[(group, side)] = np.array(
                    [self.muscle_index[s + '_' + side] for s in stems if s + '_' + side in self.muscle_index],
                    dtype=int)

        # antigravity muscles per side: what has to relax for a foot to unload
        self._relax_index = {}
        for side in 'rl':
            names = set()
            for group in ('plantar', 'kneeext', 'hipext', 'adduct'):
                names.update(int(i) for i in self.group_index[(group, side)])
            self._relax_index[side] = np.array(sorted(names), dtype=int)
        self._relax_cache = {}

        # parsed once; state_vector runs every 10 ms
        self._plan = []
        for path in self.state_names:
            parts = path.strip('/').split('/')
            name, variable = parts[-2], parts[-1]
            if parts[0] == 'jointset':
                self._plan.append((0, name, variable))
            else:
                self._plan.append((1, name, 'fiber_length_m' if variable == 'fiber_length' else variable))
        self._buffer = np.empty(len(self.state_names))

    def gain_for(self, q_tx_speed):
        """Deflated gain for one forward-velocity stage cost.  Cached per value."""
        w = float(q_tx_speed)
        if w <= 0:
            return self.K, self._gain_cache[0.0]
        # quantise to one significant figure so a search over a continuous knob
        # does not resolve a 255-state DARE for every candidate
        exponent = math.floor(math.log10(w))
        w = round(w / 10 ** exponent, 1) * 10 ** exponent
        hit = self._gain_cache.get(w)
        if hit is None:
            hit = design_gain(self._artifact, self._dt, q_tx_speed=w)
            if len(self._gain_cache) > 24:
                self._gain_cache = {0.0: self._gain_cache[0.0]}
            self._gain_cache[w] = hit
        return hit['K'], hit

    def antigravity_scale(self, side, scale):
        """Multiplier on u0 for one side's weight-bearing muscles."""
        key = (side, scale)
        mask = self._relax_cache.get(key)
        if mask is None:
            mask = np.ones(len(self.muscle_names))
            mask[self._relax_index[side]] = scale
            if len(self._relax_cache) > 8:
                self._relax_cache = {}
            self._relax_cache[key] = mask
        return mask

    def state_vector(self, snapshot):
        values = self._buffer
        coords, muscles = snapshot['coordinates'], snapshot['muscles']
        for i, (kind, name, variable) in enumerate(self._plan):
            values[i] = (coords if kind == 0 else muscles)[name][variable]
        if not np.isfinite(values).all():
            raise ValueError('Nonfinite native feedback')
        return values


# ------------------------------------------------------------- gait controller
SWING_ARTIFACT = (ROOT / 'data/research/locomotion_control/linearization_j0b9dci6'
                  / 'deflated_margin/linearization.npz')


def load_single_support(policy, artifact=SWING_ARTIFACT):
    """The repo's own right-swing single-support equilibrium and deflated gain.

    ``docs/NATIVE_LANDING_GROUND_PLANE_SYMMETRY.md`` section 4 designs this gain
    on the plant with the RIGHT foot airborne; on that plant it has a closed-loop
    spectral radius of 0.9786 excluding the symmetry modes, where the standing
    gain has no guarantee at all.  State and muscle orderings are verified equal
    to the stance artifact, and the target mass is identical; the model digest
    differs only because the waypoint's model file bakes different
    ``default_activation`` values, which is an initial condition, not dynamics.

    There is NO mirrored left-swing waypoint in the repository, so this can only
    be applied to a right swing.
    """
    if not artifact.exists():
        return None
    with np.load(artifact, allow_pickle=False) as d:
        data = {k: d[k] for k in d.files}
    if [str(x) for x in data['state_names']] != policy.state_names:
        raise ValueError('single-support artifact state ordering differs')
    if [str(x) for x in data['muscle_names']] != policy.muscle_names:
        raise ValueError('single-support artifact muscle ordering differs')
    if abs(float(data['target_mass_kg']) - TARGET_MASS_KG) > 1e-9:
        raise ValueError('single-support artifact identification mass differs')
    return {'K': data['K'], 'x0': data['x0'].copy(), 'u0': data['u0'].copy(),
            'artifact': str(artifact.relative_to(ROOT)),
            'gain_max': float(np.abs(data['K']).max())}


class GaitController:
    """settle -> ds -> swing -> reach -> ds(other side) -> ...

    ``ds`` and ``reach`` end on measured foot load; ``swing`` ends on its clock.
    """

    def __init__(self, policy: Policy, params, max_steps=None):
        self.p = policy
        self.q = clamp_params(params)
        # After max_steps completed cycles the machine stops driving and hands
        # the plant back to the plain stance regulator.  Used to separate "the
        # step itself is fatal" from "starting the next step is fatal".
        self.max_steps = max_steps
        self.dx = np.zeros(len(policy.state_names))
        self.phase = 'settle'
        self.phase_t = 0.0
        self.swing = 'r'          # side that will swing next / is swinging
        self.steps = []           # completed steps
        self.events = []
        self.settle_s = 0.40
        self.bias_now = np.zeros(len(policy.muscle_names))
        self.lift_x = self.lift_y = 0.0
        self.lift_t = 0.0
        self.min_load, self.max_clear = 1.0, 0.0
        self.single_support = None      # set by rollout when enabled
        self.ss_blend = 0.0
        self.ss_frames = 0

    # -- reference offset for the current phase -----------------------------
    def target_offset(self, load):
        q = self.q
        p = self.p
        dx = np.zeros(len(p.state_names))
        s, t = self.swing, ('l' if self.swing == 'r' else 'r')  # swing, stance

        def put(coord, value):
            i = p.coord_value_index.get(coord)
            if i is not None:
                dx[i] += value

        def putspeed(coord, value):
            i = p.coord_speed_index.get(coord)
            if i is not None:
                dx[i] += value

        if self.phase in ('settle', 'hold'):
            return dx
        # forward velocity reference is what actually asks the body to travel
        putspeed('pelvis_tx', q['v_forward'])
        put('pelvis_tilt', q['a_tilt'])
        # lateral weight shift: adduct the stance hip and roll/bend toward it so
        # the pelvis rides over the foot that is about to carry everything.
        sign = 1.0 if t == 'r' else -1.0   # +z is the subject's right
        put('hip_adduction_' + t, q['a_add'])
        put('pelvis_list', sign * q['a_list'])
        put('lumbar_bending', sign * q['a_bend'])
        if self.phase == 'ds':
            # the trailing leg (the one about to swing) pushes off
            put('ankle_angle_' + s, -q['a_push_ds'])
            return dx
        # stance leg lowers the pelvis (the touchdown hypothesis) and extends the
        # hip to carry the body over the foot
        put('knee_angle_' + t, q['a_stance_knee'])
        put('ankle_angle_' + t, q['a_stance_ankle'])
        put('hip_flexion_' + t, -q['a_stance_hipext'])
        if self.phase == 'swing':
            put('hip_flexion_' + s, q['a_hip'])
            put('knee_angle_' + s, q['a_knee'])
            put('ankle_angle_' + s, q['a_ankle'])
            put('hip_adduction_' + s, -q['a_swing_abd'])
            put('ankle_angle_' + t, -q['a_push'])
        else:  # reach: extend the leg out and down, wait for contact
            put('hip_flexion_' + s, q['a_hip_land'])
            put('knee_angle_' + s, q['a_knee_land'])
            put('ankle_angle_' + s, q['a_ankle_land'])
            put('hip_adduction_' + s, -q['a_swing_abd'])
            put('ankle_angle_' + t, -q['a_push'])
        return dx

    def bias(self):
        p, q = self.p, self.q
        b = np.zeros(len(p.muscle_names))
        if self.phase in ('settle', 'hold'):
            return b
        s, t = self.swing, ('l' if self.swing == 'r' else 'r')
        b[p.group_index[('abduct', t)]] += q['b_abduct']
        b[p.group_index[('kneeext', t)]] += q['b_vasti']
        b[p.group_index[('hipext', t)]] += q['b_hipext']
        if self.phase == 'ds':
            b[p.group_index[('plantar', s)]] += q['b_push']    # trailing push-off
        elif self.phase == 'swing':
            b[p.group_index[('hipflex', s)]] += q['b_hipflex']
            b[p.group_index[('dorsi', s)]] += q['b_dorsi']
            b[p.group_index[('kneeflex', s)]] += q['b_kneeflex']
            b[p.group_index[('plantar', t)]] += q['b_push']
        elif self.phase == 'reach':
            b[p.group_index[('dorsi', s)]] += q['b_dorsi'] * 0.5
            b[p.group_index[('kneeext', s)]] += q['b_kneeext']   # extend to reach the floor
            b[p.group_index[('plantar', t)]] += q['b_push']
        return b

    # -- transitions --------------------------------------------------------
    # A step is only counted when the foot actually left the ground and came
    # down somewhere new.  A foot that merely unweights to 13% of body weight
    # and reloads is a weight shift, and the first search run found exactly that
    # (3 mm of "step"), so the thresholds below are FIXED, not searchable.
    AIRBORNE_LOAD = 0.02      # fraction of body weight
    MIN_CLEARANCE_M = 0.010   # rise of the foot contact centre above lift-off
    MIN_ADVANCE_M = 0.030     # forward travel of the foot over the step

    def advance_phase(self, t_s, load, centres):
        q = self.q
        s = self.swing
        moved = False
        if self.phase in ('swing', 'reach'):
            self.min_load = min(self.min_load, load[s])
            self.max_clear = max(self.max_clear, centres[s][1] - self.lift_y)
        if self.phase == 'hold':
            return
        if self.max_steps is not None and len(self.steps) >= self.max_steps:
            self.phase, self.phase_t = 'hold', 0.0
            self.events.append({'time_s': t_s, 'phase': 'hold', 'swing': s})
            return
        if self.phase == 'settle':
            if t_s >= self.settle_s:
                self.phase, self.phase_t, moved = 'ds', 0.0, True
        elif self.phase == 'ds':
            if load[s] <= q['load_off'] or self.phase_t >= q['t_ds']:
                self.lift_x, self.lift_y = centres[s][0], centres[s][1]
                self.lift_t = t_s
                self.min_load, self.max_clear = load[s], 0.0
                self.phase, self.phase_t, moved = 'swing', 0.0, True
        elif self.phase == 'swing':
            if self.phase_t >= q['t_swing']:
                self.phase, self.phase_t, moved = 'reach', 0.0, True
        elif self.phase == 'reach':
            landed = load[s] >= q['load_on']
            if landed or self.phase_t >= q['t_reach_max']:
                advance = centres[s][0] - self.lift_x
                airborne = self.min_load <= self.AIRBORNE_LOAD
                self.steps.append({
                    'side': s, 'lift_time_s': self.lift_t, 'land_time_s': t_s,
                    'foot_advance_m': advance,
                    'peak_clearance_m': self.max_clear,
                    'minimum_load_fraction': self.min_load,
                    'landed_on_contact': bool(landed),
                    'airborne': bool(airborne),
                    'genuine': bool(landed and airborne
                                    and self.max_clear >= self.MIN_CLEARANCE_M
                                    and advance >= self.MIN_ADVANCE_M),
                })
                self.swing = 'l' if s == 'r' else 'r'
                self.phase, self.phase_t, moved = 'ds', 0.0, True
        if moved:
            self.events.append({'time_s': t_s, 'phase': self.phase, 'swing': self.swing})

    def commands(self, t_s, x, load, centres):
        self.advance_phase(t_s, load, centres)
        p, q = self.p, self.q
        target = self.target_offset(load)
        alpha = min(DT / max(q['tau'], DT), 1.0)
        # both the reference offset and the muscle bias are low-passed: a phase
        # entry must not be a step change in the command vector.
        self.dx += (target - self.dx) * alpha
        self.bias_now += (self.bias() - self.bias_now) * alpha
        # u = u0 - K (x - x0 - dx) splits exactly into a *regulation* term
        # -K (x - x0) and a *feedforward* term +K dx.  They are scaled
        # separately: max|K| is 156, so a 0.5 rad reference offset fed through
        # the gain at unit weight saturates most of the 98 commands and the
        # plant flails (measured: 40/98 clipped, pelvis_tilt ran to +0.6 rad).
        # The regulator is what balances; the feedforward only has to nudge the
        # right synergy, so it gets its own small weight.
        # The forward-velocity reference has to be tracked at full loop gain.
        # Routing it through ff_gain (~0.08) while the regulator damped the
        # actual pelvis_tx velocity at full gain made it inert: an ablation with
        # v_forward = 0 moved the pelvis the same 6 mm and fell at the same time.
        dx_velocity = np.zeros_like(self.dx)
        i_tx = p.coord_speed_index.get('pelvis_tx')
        if i_tx is not None:
            dx_velocity[i_tx] = self.dx[i_tx]
        dx_position = self.dx - dx_velocity
        # Blend toward the single-support equilibrium and gain while the right
        # foot is actually unloading.  The repo's own finding is that a handoff
        # is only safe once the plant has already become the target plant, so the
        # blend is driven by measured foot load, not by a clock.
        K, _ = p.gain_for(q['q_tx_speed'])
        x0, u0 = p.x0, p.u0
        if self.single_support is not None and self.swing == 'r':
            want = 1.0 - min(max(load['r'] / 0.5, 0.0), 1.0) if self.phase in ('swing', 'reach') else 0.0
            self.ss_blend += (want - self.ss_blend) * alpha
            b = self.ss_blend
            if b > 1e-6:
                self.ss_frames += 1
                ss = self.single_support
                K = (1.0 - b) * K + b * ss['K']
                x0 = (1.0 - b) * p.x0 + b * ss['x0']
                u0 = (1.0 - b) * p.u0 + b * ss['u0']
        error = x - x0 - dx_velocity
        swing_mask = p.side_state_mask[self.swing]
        # the swing leg regulates itself toward the standing pose only weakly
        # (or it will refuse to leave it); everything else sees the swing leg's
        # deviation attenuated so a commanded swing does not corrupt balance.
        swing_reg = 1.0 if self.phase == 'hold' else q['swing_reg']
        error_swing = error * (1.0 - swing_mask * (1.0 - swing_reg))
        error_other = error * (1.0 - swing_mask * (1.0 - q['cross_gate']))
        feedforward = q['ff_gain'] * (K @ dx_position)
        # The swing leg's own antigravity muscles hold it at its stance
        # excitation, so the foot keeps pressing on the ground and never
        # unloads.  Scale that baseline down while the leg is meant to be off.
        u0_swing, u0_other = u0, u0
        if self.phase in ('swing', 'reach') and self.phase != 'hold':
            stance = 'l' if self.swing == 'r' else 'r'
            u0_swing = u0 * p.antigravity_scale(self.swing, q['swing_relax'])
            u0_other = u0 * p.antigravity_scale(stance, q['stance_boost'])
        u_swing = u0_swing - q['gain'] * (K @ error_swing) + feedforward
        u_other = u0_other - q['gain'] * (K @ error_other) + feedforward
        m = p.side_muscle_mask[self.swing]
        raw = m * u_swing + (1.0 - m) * u_other + self.bias_now
        u = np.clip(raw, 0.0, 1.0)
        self.phase_t += DT
        return u, int(np.count_nonzero(raw != u))


# ------------------------------------------------------------------- rollout
def foot_loads(state):
    total = TARGET_MASS_KG * 9.81
    f = state['foot_contact_force_n']
    return {'r': f['r'] / total, 'l': f['l'] / total}


def foot_centres(state):
    out = {}
    for side in 'rl':
        pts = [c['center_m'] for c in state['contacts']
               if c['name'].startswith('contact') and c['name'].endswith('_' + side)]
        out[side] = [sum(p[i] for p in pts) / len(pts) for i in range(3)]
    return out


def fall_load(state):
    return sum(math.sqrt(sum(v * v for v in c['force_n'])) for c in state['contacts']
               if c['name'].startswith('fall_support'))


def diverged(state, tx0, tz0):
    """Stop before the state is violent enough to stall the native integrator.

    A flailing plant drives the Simbody error controller to tens of thousands of
    internal steps for one 10 ms advance, which costs minutes of wall time per
    rollout and returns nothing useful.  These guards are all falls or worse.
    """
    coords = state['coordinates']
    if fall_load(state) > 5.0:
        return 'fall_support_contact'
    if coords['pelvis_ty']['value'] < 0.80:
        return 'pelvis_below_0.80m'
    if abs(coords['pelvis_tilt']['value']) > 0.55:
        return 'pelvis_tilt_beyond_0.55rad'
    if abs(coords['pelvis_list']['value']) > 0.45:
        return 'pelvis_list_beyond_0.45rad'
    if tz0 is not None and abs(coords['pelvis_tz']['value'] - tz0) > 0.30:
        return 'lateral_excursion_beyond_0.30m'
    for name in ('pelvis_tx', 'pelvis_ty', 'pelvis_tz'):
        if abs(coords[name]['speed']) > 4.0:
            return 'pelvis_speed_beyond_4m_s'
    for name, c in coords.items():
        if abs(c['speed']) > 25.0:
            return 'joint_speed_beyond_25_' + name
    return None


def rollout(params, policy, horizon_s=5.0, record=False, work_dir=None, wall_budget_s=300.0,
            max_steps=None, single_support=False, arm_hold=True):
    """Integrate the plant under the controller.  Returns a report dict."""
    params = clamp_params(params)
    out = Path(work_dir or (WORK / ('run-' + uuid.uuid4().hex)))
    if out.exists():
        shutil.rmtree(out)
    pose = json.loads((BUNDLE / 'initial_pose.json').read_text())
    arm_targets = arm_hold_targets(pose)
    controller = GaitController(policy, params, max_steps=max_steps)
    if single_support:
        controller.single_support = load_single_support(policy)
    frames, native = [], None
    fell, failure, clipped_max = False, None, 0
    stop_reason = 'horizon'
    started_wall = time.time()
    tx0 = ty0 = tz0 = None
    tx_max = -1e9
    tx_speed_max = 0.0
    tz_min, tz_max = 1e9, -1e9
    com_h_min, com_h_max = 1e9, -1e9
    t_end = 0.0
    try:
        native = NativeMechanicalStream(ROOT, out, environment='upright',
                                        target_mass_kg=TARGET_MASS_KG, initial_pose=pose,
                                        augmented_registration=REGISTRATION)
        state = native.snapshot()
        n_steps = int(round(horizon_s / DT))
        for _ in range(n_steps):
            t_s = state['time_s']
            coords = state['coordinates']
            tx, ty, tz = coords['pelvis_tx']['value'], coords['pelvis_ty']['value'], coords['pelvis_tz']['value']
            if tx0 is None:
                tx0, ty0, tz0 = tx, ty, tz
            tx_max = max(tx_max, tx - tx0)
            tx_speed_max = max(tx_speed_max, coords['pelvis_tx']['speed'])
            tz_min, tz_max = min(tz_min, tz), max(tz_max, tz)
            com_h = state['potential_energy_j'] / (state['mass_kg'] * 9.81)
            com_h_min, com_h_max = min(com_h_min, com_h), max(com_h_max, com_h)
            load = foot_loads(state)
            centres = foot_centres(state)
            x = policy.state_vector(state)
            u, clipped = controller.commands(t_s, x, load, centres)
            clipped_max = max(clipped_max, clipped)
            if record:
                frames.append({
                    'time_s': t_s,
                    'joints': {k: {'value': v['value'], 'speed': v['speed'], 'unit': v['unit']}
                               for k, v in coords.items()},
                    'contacts': {
                        'foot_load_fraction': load,
                        'foot_contact_force_n': dict(state['foot_contact_force_n']),
                        'foot_centre_m': centres,
                        'fall_support_force_n': fall_load(state),
                    },
                    'motor_excitations': {n: float(v) for n, v in zip(policy.muscle_names, u)},
                    'phase': controller.phase,
                    'swing_side': controller.swing,
                })
            reason = diverged(state, tx0, tz0)
            if reason is not None:
                fell, stop_reason, t_end = True, reason, t_s
                break
            if time.time() - started_wall > wall_budget_s:
                stop_reason, t_end = 'wall_budget', t_s
                break
            state = native.advance(
                DT, actuation=dict(zip(policy.muscle_names, map(float, u))),
                coordinate_actuation=(arm_hold_commands(state, arm_targets)
                                      if arm_hold else None))
            t_end = state['time_s']
    except Exception as exc:  # native rejection, timeout, integrator blow-up
        failure = type(exc).__name__ + ': ' + str(exc)[:200]
        stop_reason = 'native_failure'
        fell = True
    finally:
        if native is not None:
            try:
                native.close()
            except Exception:
                pass
        if not record:
            shutil.rmtree(out, ignore_errors=True)

    steps = controller.steps
    contact_steps = [s for s in steps if s['landed_on_contact']]
    airborne_steps = [s for s in steps if s['airborne'] and s['landed_on_contact']]
    genuine_steps = [s for s in steps if s['genuine']]

    def longest_alternating(selected):
        run, best, expect = 0, 0, None
        for s in steps:
            if s in selected and (expect is None or s['side'] == expect):
                run += 1
                best = max(best, run)
                expect = 'l' if s['side'] == 'r' else 'r'
            else:
                run, expect = 0, None
        return best

    report = {
        'params': params,
        'horizon_s': horizon_s,
        'simulated_s': t_end,
        'wall_s': time.time() - started_wall,
        'fell': fell,
        'completed_horizon': stop_reason == 'horizon',
        'stop_reason': stop_reason,
        'failure': failure,
        'steps_attempted': len(steps),
        'steps_on_contact': len(contact_steps),
        'steps_airborne': len(airborne_steps),
        'steps_genuine': len(genuine_steps),
        'consecutive_alternating_steps': longest_alternating(genuine_steps),
        'consecutive_alternating_airborne': longest_alternating(airborne_steps),
        'consecutive_alternating_contact': longest_alternating(contact_steps),
        'step_criterion': {
            'airborne_load_fraction_below': GaitController.AIRBORNE_LOAD,
            'minimum_clearance_m': GaitController.MIN_CLEARANCE_M,
            'minimum_forward_advance_m': GaitController.MIN_ADVANCE_M,
            'note': ('"genuine" requires the swing foot to unload below 2% of body '
                     'weight, rise at least 10 mm, advance at least 30 mm and land '
                     'on measured contact.  Everything weaker is reported '
                     'separately and is a weight shift, not a step.'),
        },
        'step_records': steps,
        'phase_events': controller.events,
        'pelvis_forward_travel_m': None if tx0 is None else tx_max,
        'pelvis_lateral_excursion_m': None if tx0 is None else tz_max - tz_min,
        'com_height_range_m': None if tx0 is None else [com_h_min, com_h_max],
        'com_height_excursion_m': None if tx0 is None else com_h_max - com_h_min,
        'max_clipped_muscles': clipped_max,
        'arm_hold': ({'ports': sorted(ARM_PORT[k] + '_' + s for k in ARM_PORT for s in 'rl'),
                      'gains_nm_per_rad_and_nms_per_rad': ARM_HOLD_GAINS,
                      'basis': ('declared source CoordinateActuator torque ports held at the '
                                'registered initial arm pose; torque actuators, not muscles')}
                     if arm_hold else None),
        'commanded_forward_velocity_m_s': params['v_forward'],
        'achieved_mean_forward_velocity_m_s': (
            None if tx0 is None or t_end <= 0 else tx_max / t_end),
        'forward_velocity_tracking_ratio': (
            None if tx0 is None or t_end <= 0 or params['v_forward'] <= 0
            else (tx_max / t_end) / params['v_forward']),
        'peak_forward_velocity_m_s': tx_speed_max,
        'max_steps_before_hold': max_steps,
        'single_support_gain': None if controller.single_support is None else {
            'artifact': controller.single_support['artifact'],
            'gain_max': controller.single_support['gain_max'],
            'frames_blended': controller.ss_frames,
            'right_swing_only': True,
        },
    }
    if record:
        report['output_dir'] = str(out.relative_to(ROOT))
    return report, frames


def score(report):
    if report['pelvis_forward_travel_m'] is None:
        return -1e6
    s = 250.0 * report['consecutive_alternating_steps']       # genuine only
    s += 60.0 * report['steps_genuine']
    s += 40.0 * report['consecutive_alternating_airborne']    # partial credit
    s += 15.0 * report['steps_airborne']
    s += 300.0 * report['pelvis_forward_travel_m']
    # shaping toward an actual step: without these the search happily settles on
    # a forward lean that never lifts a foot, because lifting is destabilising
    # before it is useful and a step is far away in parameter space.
    s += 1000.0 * sum(max(0.0, min(r['peak_clearance_m'], 0.10))
                      for r in report['step_records'])
    s += 250.0 * sum(max(0.0, 0.25 - r['minimum_load_fraction'])
                     for r in report['step_records'])
    s += 10.0 * report['simulated_s']
    # the known failure mode of this plant is a frontal-plane collapse off the
    # single-support base (docs/NATIVE_LANDING_GROUND_PLANE_SYMMETRY.md s1), so
    # sideways travel beyond a normal step width is charged for directly.
    s -= 200.0 * max(0.0, (report['pelvis_lateral_excursion_m'] or 0.0) - 0.12)
    # Only a run that reached its horizon counts as survived.  A run cut off by
    # the wall budget is not a success: violent states are what make the native
    # integrator slow, so "ran out of wall clock" correlates with falling.
    if report.get('completed_horizon'):
        s += 60.0
    # Falling has to cost more than one step is worth, or the search converges on
    # a single lunge that scores a step and then collapses -- which is exactly
    # what it did before this term was added.
    s -= 120.0 * max(0.0, report['horizon_s'] - report['simulated_s'])
    return s


# --------------------------------------------------------------------- search
_POLICY = None


def _worker(job):
    global _POLICY
    if _POLICY is None:
        _POLICY = Policy()
    index, params, horizon = job
    try:
        report, _ = rollout(params, _POLICY, horizon_s=horizon)
    except Exception:
        return index, {'params': clamp_params(params), 'fell': True, 'wall_s': 0.0,
                       'stop_reason': 'worker_exception', 'completed_horizon': False,
                       'failure': traceback.format_exc()[-300:], 'simulated_s': 0.0,
                       'steps_attempted': 0, 'steps_on_contact': 0, 'steps_airborne': 0, 'steps_genuine': 0, 'consecutive_alternating_airborne': 0, 'consecutive_alternating_contact': 0,
                       'consecutive_alternating_steps': 0, 'step_records': [],
                       'phase_events': [], 'pelvis_forward_travel_m': None,
                       'pelvis_lateral_excursion_m': None, 'com_height_range_m': None,
                       'com_height_excursion_m': None, 'max_clipped_muscles': 0,
                       'horizon_s': horizon}, -1e6
    return index, report, score(report)


def _kill_orphan_engines():
    """Native engines outlive a killed worker; they hold cores and never exit."""
    import signal, subprocess
    try:
        out = subprocess.run(['ps', '-eo', 'pid,args'], capture_output=True, text=True).stdout
    except Exception:
        return
    for line in out.splitlines():
        if 'native_mechanical_stream' in line and 'gait-work' in line:
            try:
                os.kill(int(line.split()[0]), signal.SIGKILL)
            except Exception:
                pass


def search(generations, population, horizon, workers, sigma0, out_path, seed=0, initial=None,
           generation_timeout_s=600.0):
    import multiprocessing as mp
    WORK.mkdir(parents=True, exist_ok=True)
    Policy()                      # populate the gain cache once, not per worker
    rng = np.random.default_rng(seed)
    span = np.array([BOUNDS[k][1] - BOUNDS[k][0] for k in ORDER])
    free_mask = np.array([1.0 if k in FREE else 0.0 for k in ORDER])
    incumbent = np.array([(initial or SEED_PARAMS)[k] for k in ORDER], float)
    best_report = None
    best_score = -1e9
    history = []
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
                    step = rng.normal(0, sigma, len(ORDER)) * span * free_mask
                    candidate = incumbent + step
                lo = np.array([BOUNDS[k][0] for k in ORDER])
                hi = np.array([BOUNDS[k][1] for k in ORDER])
                candidate = np.clip(candidate, lo, hi)
                jobs.append((i, dict(zip(ORDER, candidate)), horizon))
            # A worker that dies (the native engine is memory hungry and this
            # machine is shared) makes Pool.map hang forever waiting for a
            # result that will never arrive.  Bound it, and rebuild the pool.
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
            gen_best = top[2]
            if gen_best > best_score:
                best_score, best_report = gen_best, top[1]
                incumbent = np.array([best_report['params'][k] for k in ORDER])
                sigma = sigma0
            else:
                sigma = max(sigma * 0.92, sigma0 * 0.45)
            row = {
                'generation': generation, 'wall_s': time.time() - started,
                'sigma': sigma, 'generation_best_score': gen_best, 'best_score': best_score,
                'best_consecutive_steps': best_report['consecutive_alternating_steps'],
                'best_steps_genuine': best_report['steps_genuine'],
                'best_steps_airborne': best_report['steps_airborne'],
                'best_alt_airborne': best_report['consecutive_alternating_airborne'],
                'best_travel_m': best_report['pelvis_forward_travel_m'],
                'best_fell': best_report['fell'],
                'best_completed': best_report.get('completed_horizon'),
                'best_stop': best_report.get('stop_reason'),
            }
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


# ----------------------------------------------------------------------- main
def write_best(params, policy, horizon, destination, arm_hold=True):
    destination.mkdir(parents=True, exist_ok=True)
    work = WORK / ('best-' + uuid.uuid4().hex)
    report, frames = rollout(params, policy, horizon_s=horizon, record=True, work_dir=work,
                             arm_hold=arm_hold, wall_budget_s=1e9)
    trajectory = {
        'schema': 'ihm.gait-trajectory.v1',
        'dt_s': DT,
        'model_sha256': policy.model_sha256,
        'target_mass_kg': TARGET_MASS_KG,
        'muscle_names': policy.muscle_names,
        'basis': ('Muscle excitation only.  No prescribed coordinate, no external root force, '
                  'no motion constraint.  Joint values are integrated native OpenSim/Simbody output.'),
        'frames': frames,
    }
    (destination / 'trajectory.json').write_text(json.dumps(trajectory) + '\n')
    report['controller'] = {
        'kind': 'contact_triggered_fsm_over_deflated_lqr',
        'lqr_source': 'data/models/engineering_stance_v1/linearization.npz',
        'deflated_dimension': policy.deflated_dimension,
        'deflated_dominant_states': policy.deflated_states,
        'residual_gain_on_symmetry': policy.deflation_residual,
        'closed_loop_radius_excluding_deflated': policy.closed_loop_radius_excluding_deflated,
        'discrete_gain_max': policy.gain_max,
        'forward_velocity_stage_cost': report['params']['q_tx_speed'],
        'forward_velocity_gain_norm': float(
            policy.gain_for(report['params']['q_tx_speed'])[1]['forward_velocity_gain_norm']),
        'forward_velocity_gain_norm_as_shipped': 3.12,
        'prescribed_motion': False,
        'external_root_forces': False,
        'muscle_excitation_only': True,
    }
    report['frames'] = len(frames)
    report['trajectory_path'] = str((destination / 'trajectory.json').relative_to(ROOT))
    (destination / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    shutil.rmtree(work, ignore_errors=True)
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', choices=('single', 'search', 'best'), default='single')
    ap.add_argument('--horizon', type=float, default=5.0)
    ap.add_argument('--generations', type=int, default=8)
    ap.add_argument('--population', type=int, default=16)
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--sigma', type=float, default=0.12)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--params', help='JSON file with a parameter dict (or a search output)')
    ap.add_argument('--search-out', default='data/derived/gait-search/search.json')
    ap.add_argument('--out', default='data/derived/gait-best')
    ap.add_argument('--no-arm-hold', action='store_true',
                    help='leave the declared arm torque ports at zero, as before')
    a = ap.parse_args()

    WORK.mkdir(parents=True, exist_ok=True)
    params = dict(SEED_PARAMS)
    if a.params:
        loaded = json.loads(Path(a.params).read_text())
        if 'best_report' in loaded:
            loaded = loaded['best_report']['params']
        elif 'params' in loaded:
            loaded = loaded['params']
        params.update(loaded)

    if a.mode == 'search':
        best, best_score, _ = search(a.generations, a.population, a.horizon, a.workers,
                                     a.sigma, ROOT / a.search_out, seed=a.seed, initial=params)
        print(json.dumps({'best_score': best_score,
                          'consecutive_steps': best['consecutive_alternating_steps'],
                          'travel_m': best['pelvis_forward_travel_m']}, indent=2))
        return
    policy = Policy()
    if a.mode == 'best':
        report = write_best(params, policy, a.horizon, ROOT / a.out,
                            arm_hold=not a.no_arm_hold)
    else:
        report, _ = rollout(params, policy, horizon_s=a.horizon,
                            arm_hold=not a.no_arm_hold, wall_budget_s=1e9)
    printable = {k: v for k, v in report.items() if k not in ('params', 'phase_events')}
    print(json.dumps(printable, indent=2))


if __name__ == '__main__':
    main()
