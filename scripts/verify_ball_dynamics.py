#!/usr/bin/env python3
"""Broadened rigid-sphere contact checks: normal, oblique, rolling, prop and body.

SCOPE. Every "analytic" comparison in this module is a closed-form solution of the
*stated* rigid-impulse contact model (instantaneous Newton restitution plus a
clamped Coulomb tangential impulse on a uniform solid sphere). Agreement shows the
solver implements the model it declares, and that time-of-impact splitting does not
leak energy or momentum. It is NOT a measurement against a real rubber ball: the
coefficients (restitution .75, friction .35) are uncalibrated engineering values.

Cases with no closed form -- ball against the compliant penalty colliders used for
furniture, cloth and body skin -- are checked for the conservation and containment
properties the implementation must have (equal/opposite exchange, no energy gain,
no tunnelling). Those are consistency checks, explicitly not physics validation.

Run: PYTHONPATH=.:scripts .venv/bin/python -m unittest scripts.verify_ball_dynamics
Receipt: PYTHONPATH=. .venv/bin/python scripts/verify_ball_dynamics.py
"""
import json
import unittest
from pathlib import Path

import numpy as np

from ihm.assembly.environment_dynamics import EnvironmentDynamics, RigidProp
from ihm.assembly.rigid_contact import BOUNCY_BALL, sphere_plane_step

ROOT = Path(__file__).resolve().parents[1]
G = 9.81
# Uniform solid sphere: I = 2/5 m r^2, so the contact-point effective mass for a
# tangential impulse is m / (1 + m r^2 / I) = m / 3.5, i.e. the classical 2/7 law.
SLIP_FACTOR = 3.5

RECORD = {'bounds_m': {'min': [-.1] * 3, 'max': [.1] * 3}, 'mass_kg': .4, 'radius_m': .1}
ENTITY = {'bone': {'centroid_m': [0, 0, 0], 'rotation_matrix': np.eye(3).tolist()}}


def fixture(points, *, gravity=(0, 0, -G), rigid=(), static=(), axis=2, base_plane=-1.):
    """Held-canonical EnvironmentDynamics with no catalogue or native registration."""
    owner = EnvironmentDynamics.__new__(EnvironmentDynamics)
    owner.ids = ['bone'] * len(points)
    owner.offsets = np.array(points, float)
    owner.previous = None
    owner.soft = []
    owner.rigid = list(rigid)
    owner.static = list(static)
    owner.object_static = []
    owner.gravity = np.array(gravity, float)
    owner.axis = axis
    owner.base_plane = base_plane
    owner.object_plane = base_plane
    owner.time_s = 0.
    owner.contacts = []
    owner.last_impulse = np.zeros(3)
    owner.selection = {}
    owner.placements = []
    owner.sources = {}
    return owner


def ball(offset=(0, 0, 0), *, mass=.4, radius=.1):
    record = dict(RECORD, mass_kg=mass, radius_m=radius,
                  bounds_m={'min': [-radius] * 3, 'max': [radius] * 3})
    return RigidProp('ball', record, np.array(offset, float))


def energy(prop, gravity, axis):
    inertia = .4 * prop.mass * prop.radius ** 2
    return (.5 * prop.mass * prop.v @ prop.v + .5 * inertia * prop.omega @ prop.omega
            - prop.mass * gravity[axis] * prop.x[axis])


def impact_step(v, omega, *, e=.75, mu=.35, radius=.1, mass=.4, axis=2, plane=0.,
                threshold=.05, dt=1e-4, gravity=0.):
    """One isolated impact at the plane, in free flight, and the resulting state.

    gravity defaults to zero so the impulse law is measured on its own. With
    gravity on, the post-impact velocity also carries -g*(dt - t_impact) of free
    flight, which is what makes a naive rebound-ratio reading fall short of e.
    """
    x = np.zeros(3)
    x[axis] = plane + radius
    v = np.array(v, float)
    omega = np.array(omega, float)
    acceleration = np.zeros(3)
    acceleration[axis] = -gravity
    loss, impulse, repair = sphere_plane_step(x, v, omega, dt, acceleration, radius, mass,
                                              axis, plane, e, mu, threshold)
    return x, v, omega, loss, impulse, repair


def analytic_impact(vt, vn, omega_t, *, e=.75, mu=.35, radius=.1, mass=.4):
    """Closed-form post-impact tangential velocity/spin for the stated impulse law.

    The contact point sits at lever -r*n below the centre, so its tangential speed
    is vt - omega_t*r for a spin omega_t about the axis completing the right-handed
    triad. vn>0 is the closing normal speed.
    Returns (tangential centre velocity, spin, sticking flag).
    """
    inertia = .4 * mass * radius ** 2
    slip = vt - omega_t * radius            # contact-point tangential speed
    jn = (1 + e) * mass * vn
    stick = mass * abs(slip) / SLIP_FACTOR
    sticking = stick <= mu * jn
    jt = -np.sign(slip) * min(stick, mu * jn)
    return vt + jt / mass, omega_t - jt * radius / inertia, sticking


class Analytic(unittest.TestCase):
    """Closed-form solutions of the declared rigid-impulse model."""

    def test_normal_restitution_ratio_is_exact_for_a_range_of_drops(self):
        # The rebound *speed* ratio is sampling-free and must equal e exactly.
        for height in (.05, .25, 1., 2., 5.):
            with self.subTest(height=height):
                approach = np.sqrt(2 * G * height)
                _, v, omega, loss, impulse, _ = impact_step([0, 0, -approach], [0, 0, 0])
                self.assertAlmostEqual(v[2] / approach, BOUNCY_BALL['restitution'], places=12)
                np.testing.assert_allclose(omega, 0, atol=1e-15)
                # Impulse and dissipation match the closed-form values.
                self.assertAlmostEqual(impulse[2], (1 + .75) * .4 * approach, places=12)
                self.assertAlmostEqual(loss, .5 * .4 * approach ** 2 * (1 - .75 ** 2), places=12)

    def test_rebound_ratio_shortfall_under_gravity_is_exactly_the_free_flight_term(self):
        """Naive readings undershoot e by g*dt/approach, not by a solver error."""
        dt = 1e-4
        for height in (.05, 1., 5.):
            with self.subTest(height=height):
                approach = np.sqrt(2 * G * height)
                _, v, _, _, _, _ = impact_step([0, 0, -approach], [0, 0, 0], gravity=G, dt=dt)
                self.assertAlmostEqual(v[2], .75 * approach - G * dt, places=12)

    def test_drop_apex_matches_e_squared_height_when_sampled_finely(self):
        """h' = e^2 h. Coarse tick sampling is the only source of shortfall."""
        rows = []
        for height in (.25, 1., 2.):
            prop = ball((0, 0, -1. + .1 + height))
            owner = fixture([[0, 50, 0]], rigid=[prop])
            start = prop.x[2] - (owner.base_plane + prop.radius)
            fine, seen, previous = [], False, 0.
            for _ in range(4000):
                owner.advance(.001, ENTITY)   # 1 ms sampling of the apex
                if prop.v[2] > 0 and previous < 0:
                    seen = True
                previous = prop.v[2]
                if seen and prop.v[2] >= 0:
                    fine.append(prop.x[2] - (owner.base_plane + prop.radius))
            apex = max(fine)
            rows.append({'drop_m': start, 'apex_m': apex, 'analytic_m': .75 ** 2 * start,
                         'relative_error': apex / (.75 ** 2 * start) - 1})
            self.assertTrue(seen)
            self.assertAlmostEqual(apex / start, .75 ** 2, delta=2e-4)
        self.rows = rows

    def test_oblique_impact_matches_closed_form_in_both_friction_regimes(self):
        approach = 3.
        # Shallow tangential speed -> Coulomb budget suffices -> contact sticks,
        # so the centre keeps 5/7 of its tangential speed (classical 2/7 loss).
        for tangential, expect_stick in ((.5, True), (12., False)):
            with self.subTest(tangential=tangential):
                _, v, omega, loss, impulse, _ = impact_step(
                    [tangential, 0, -approach], [0, 0, 0])
                vt, spin, sticking = analytic_impact(tangential, approach, 0.)
                self.assertIs(sticking, expect_stick)
                self.assertAlmostEqual(v[0], vt, places=12)
                self.assertAlmostEqual(v[2], .75 * approach, places=12)
                self.assertAlmostEqual(omega[1], spin, places=12)
                np.testing.assert_allclose([v[1], omega[0], omega[2]], 0, atol=1e-15)
                self.assertGreaterEqual(loss, 0.)
                if expect_stick:
                    # Sticking means zero contact-point slip afterwards.
                    self.assertAlmostEqual(v[0] - omega[1] * .1, 0., places=12)
                    self.assertAlmostEqual(v[0], 5 / 7 * tangential, places=12)
                else:
                    self.assertAlmostEqual(v[0], tangential - .35 * 1.75 * approach, places=12)

    def test_oblique_reflection_angle_follows_the_impulse_law(self):
        approach, tangential = 4., 20.       # deliberately sliding
        _, v, _, _, _, _ = impact_step([tangential, 0, -approach], [0, 0, 0])
        incoming = np.arctan2(tangential, approach)
        outgoing = np.arctan2(v[0], v[2])
        analytic = np.arctan2(tangential - .35 * 1.75 * approach, .75 * approach)
        self.assertAlmostEqual(outgoing, analytic, places=12)
        # Not specular: restitution scales the normal component by e=.75 while
        # sliding friction scales the tangential one by .878 here, so this impact
        # rebounds SHALLOWER than it arrived. The direction of that inequality
        # depends on the incidence, so compare the two ratios rather than guess.
        self.assertGreater(outgoing, incoming)
        self.assertGreater((tangential - .35 * 1.75 * approach) / tangential, .75)

    def test_steep_impact_rebounds_steeper_than_it_arrived(self):
        approach, tangential = 4., 1.       # sticking regime: tangential keeps 5/7
        _, v, _, _, _, _ = impact_step([tangential, 0, -approach], [0, 0, 0])
        self.assertLess(np.arctan2(v[0], v[2]), np.arctan2(tangential, approach))
        self.assertLess(5 / 7, .75)         # the reason, stated explicitly

    def test_spin_can_reverse_the_tangential_direction(self):
        approach, tangential, spin = 3., 1., -60.
        _, v, omega, _, _, _ = impact_step([tangential, 0, -approach], [0, spin, 0])
        vt, _, sticking = analytic_impact(tangential, approach, spin)
        self.assertFalse(sticking)           # Coulomb budget saturates
        self.assertAlmostEqual(v[0], vt, places=12)
        self.assertLess(v[0], 0.)            # leaves the surface travelling backwards

    def test_sliding_sphere_reaches_five_sevenths_rolling(self):
        """Classical result: a sliding sphere rolls at 5/7 v0 after 2 v0 / (7 mu g)."""
        speed = 4.
        prop = ball((0, 0, -.9))             # already resting on the plane
        prop.v[:] = [speed, 0, 0]
        owner = fixture([[0, 50, 0]], rigid=[prop])
        analytic_time = 2 * speed / (7 * .35 * G)
        rolled_at, samples = None, []
        for tick in range(2000):
            owner.advance(.001, ENTITY)
            slip = prop.v[0] - prop.omega[1] * prop.radius
            samples.append({'t': (tick + 1) * .001, 'v_x': float(prop.v[0]), 'slip': float(slip)})
            if rolled_at is None and abs(slip) < 1e-9:
                rolled_at = (tick + 1) * .001
        self.assertIsNotNone(rolled_at)
        self.assertAlmostEqual(prop.v[0], 5 / 7 * speed, delta=2e-3)
        self.assertAlmostEqual(rolled_at, analytic_time, delta=1.5e-3)
        # Once rolling, an ideal sphere neither speeds up nor slows down.
        self.assertAlmostEqual(prop.v[0], samples[-1]['v_x'], delta=1e-9)
        self.samples = samples

    def test_slip_decays_at_the_analytic_seven_halves_mu_g_rate(self):
        speed = 4.
        prop = ball((0, 0, -.9))
        prop.v[:] = [speed, 0, 0]
        owner = fixture([[0, 50, 0]], rigid=[prop])
        rate = SLIP_FACTOR * .35 * G
        for tick in range(1, 60):            # sample well inside the sliding phase
            owner.advance(.001, ENTITY)
            slip = prop.v[0] - prop.omega[1] * prop.radius
            self.assertAlmostEqual(slip, speed - rate * tick * .001, places=9)

    def test_rolling_ball_keeps_rolling_and_does_not_creep_through_the_plane(self):
        prop = ball((0, 0, -.9))
        prop.v[:] = [2., 0, 0]
        prop.omega[:] = [0, 2. / .1, 0]      # already rolling without slipping
        owner = fixture([[0, 50, 0]], rigid=[prop])
        for _ in range(3000):
            owner.advance(.001, ENTITY)
        self.assertAlmostEqual(prop.v[0], 2., delta=1e-9)
        self.assertAlmostEqual(prop.x[2], owner.base_plane + prop.radius, places=12)
        self.assertEqual(prop.contact_projection_m, 0.)


class Conservation(unittest.TestCase):
    """Ledger properties that must hold whether or not a closed form exists."""

    def test_energy_ledger_closes_through_oblique_bounces_and_roll_out(self):
        prop = ball((0, 0, .4))
        prop.v[:] = [2.5, 1.5, 0]
        owner = fixture([[0, 50, 0]], rigid=[prop])
        initial = energy(prop, owner.gravity, owner.axis)
        worst = 0.
        for _ in range(6000):
            owner.advance(.001, ENTITY)
            worst = max(worst, abs(energy(prop, owner.gravity, owner.axis)
                                   + prop.contact_dissipation_j - initial))
        self.assertLess(worst, 1e-9)
        self.assertGreater(prop.contact_dissipation_j, 0.)
        # It must come to rest on the plane, not sink into it or hover.
        self.assertAlmostEqual(prop.x[2], owner.base_plane + prop.radius, places=9)
        self.worst = worst

    def test_contact_never_increases_kinetic_energy_over_random_impacts(self):
        rng = np.random.default_rng(20260908)
        inertia = .4 * .4 * .1 ** 2
        for _ in range(400):
            v = rng.normal(scale=4., size=3)
            v[2] = -abs(v[2]) - .2
            omega = rng.normal(scale=30., size=3)
            before = .5 * .4 * v @ v + .5 * inertia * omega @ omega
            _, out_v, out_omega, loss, impulse, _ = impact_step(v.copy(), omega.copy())
            after = .5 * .4 * out_v @ out_v + .5 * inertia * out_omega @ out_omega
            self.assertLessEqual(after, before + 1e-12)
            self.assertAlmostEqual(loss, before - after, places=12)
            # Linear momentum change equals the reported contact impulse exactly.
            np.testing.assert_allclose(.4 * (out_v - v), impulse, atol=1e-12)
            # A unilateral plane can only push, and Coulomb must stay in its cone.
            self.assertGreaterEqual(impulse[2], 0.)
            self.assertLessEqual(np.linalg.norm(impulse[:2]), .35 * impulse[2] + 1e-12)

    def test_gravity_axis_choice_does_not_change_the_bounce(self):
        """Supine (axis 2) and upright (axis 1) must give the identical trajectory."""
        histories = []
        for axis in (2, 1):
            gravity = np.zeros(3)
            gravity[axis] = -G
            offset = np.zeros(3)
            offset[axis] = -1. + .1 + 1.
            prop = ball(offset)
            prop.v[(axis + 1) % 3] = 1.5     # some tangential motion as well
            owner = fixture([[0, 50, 0]], gravity=gravity, rigid=[prop], axis=axis)
            history = []
            for _ in range(1200):
                owner.advance(.001, ENTITY)
                history.append((float(prop.x[axis]), float(prop.v[axis]),
                                float(prop.contact_dissipation_j)))
            histories.append(history)
        np.testing.assert_allclose(histories[0], histories[1], atol=0., rtol=0.)


class Penalty(unittest.TestCase):
    """Compliant-collider cases with no closed form: consistency only."""

    @staticmethod
    def _offset_collision(substep):
        """Off-axis prop/prop strike in free fall; returns momentum bookkeeping."""
        left = RigidProp('a', dict(RECORD, mass_kg=.4, radius_m=.1), np.array([-.15, 0, 0.]))
        right = RigidProp('b', dict(RECORD, mass_kg=.6, radius_m=.1), np.array([.15, .02, 0.]))
        left.v[:] = [3., 0, 0]
        owner = fixture([[0, 50, 0]], gravity=(0, 0, 0), rigid=[left, right])
        props = (left, right)
        linear = lambda: sum(p.mass * p.v for p in props)
        angular = lambda: sum(p.mass * np.cross(p.x, p.v)
                              + .4 * p.mass * p.radius ** 2 * p.omega for p in props)
        kinetic = lambda: sum(.5 * p.mass * p.v @ p.v
                              + .5 * .4 * p.mass * p.radius ** 2 * p.omega @ p.omega
                              for p in props)
        before = linear(), angular(), kinetic()
        for _ in range(round(.4 / substep)):
            owner.advance(substep, ENTITY)
        return props, before, (linear(), angular(), kinetic())

    def test_ball_against_ball_conserves_linear_momentum_exactly(self):
        """No ground here, so plane friction cannot hide a broken exchange."""
        (left, right), before, after = self._offset_collision(.001)
        np.testing.assert_allclose(after[0], before[0], atol=1e-12)
        self.assertGreater(right.v[0], 0.)                 # struck ball moves off
        self.assertLess(left.v[0], 3.)
        self.assertGreater(abs(right.omega[2]), 1e-6)      # off-axis strike spins it
        self.assertLess(after[2], before[2] + 1e-12)       # compliant, so lossy

    def test_prop_angular_momentum_drift_is_first_order_integrator_truncation(self):
        """FINDING, stated plainly: total angular momentum is NOT conserved.

        Contact forces and torques are reciprocal, but the explicit update reads
        them at the start of a substep while the centres move within it. An
        off-axis 3 m/s strike therefore leaves a spurious net angular momentum of
        about 2.5% of the angular momentum actually exchanged, at the production
        1 ms environment substep. It converges first order in the substep -
        halving the substep halves it - so it is truncation, not a broken
        exchange, but it does not vanish at the substep the runtime actually uses.
        """
        drifts, exchanged = [], None
        for substep in (.001, .0005, .00025, .000125):
            props, before, after = self._offset_collision(substep)
            np.testing.assert_allclose(after[0], before[0], atol=1e-12)
            drifts.append(float(np.linalg.norm(after[1] - before[1])))
            if exchanged is None:
                exchanged = float(np.linalg.norm(sum(
                    .4 * p.mass * p.radius ** 2 * p.omega for p in props)))
        ratios = [a / b for a, b in zip(drifts, drifts[1:])]
        for ratio in ratios:
            self.assertAlmostEqual(ratio, 2., delta=.4)     # first order in h
        self.assertLess(drifts[0] / exchanged, .03)
        self.assertGreater(drifts[0] / exchanged, .01)      # pin the actual size
        self.drifts, self.exchanged = drifts, exchanged

    def test_ball_landing_on_furniture_penetrates_the_soft_penalty_collider(self):
        """CHARACTERIZATION, not validation: measure how soft the box contact is.

        Furniture uses a 1500 N/m penalty spring sampled at six probe points. That
        is far softer than the analytic ground plane, so a fast ball sinks visibly
        into a table before being pushed out. The number is recorded, not blessed.
        """
        top = -.3
        table = (np.array([-.5, -.5, -1.]), np.array([.5, .5, top]), 'table')
        results = {}
        for drop in (.05, .2, .6):
            prop = ball((0, 0, top + .1 + drop))
            owner = fixture([[0, 50, 0]], rigid=[prop], static=[table])
            lowest = prop.x[2]
            for _ in range(6000):
                owner.advance(.001, ENTITY)
                lowest = min(lowest, prop.x[2])
            penetration = max(0., top + prop.radius - lowest)
            results[drop] = penetration
            # Containment is the property that must hold: it never passes through.
            self.assertGreater(lowest, table[0][2] + prop.radius)
            self.assertGreater(prop.x[2], top)
            self.assertGreater(prop.x[2], owner.base_plane + prop.radius + .01)
            # Peak penetration follows the linear-spring scaling v*sqrt(m/k).
            expected = np.sqrt(2 * G * drop) * np.sqrt(.4 / 1500.)
            self.assertAlmostEqual(penetration, expected, delta=.4 * expected)
        self.assertGreater(results[.6], results[.05])
        self.penetration_m = results

    def test_ball_beside_furniture_falls_past_it_to_the_floor(self):
        table = (np.array([-.5, -.5, -1.]), np.array([.5, .5, -.3]), 'table')
        prop = ball((1.2, 0, .3))
        owner = fixture([[0, 50, 0]], rigid=[prop], static=[table])
        for _ in range(12000):               # long enough for the bounces to die
            owner.advance(.001, ENTITY)
        self.assertAlmostEqual(prop.x[2], owner.base_plane + prop.radius, delta=1e-6)
        np.testing.assert_allclose(prop.v, 0, atol=1e-9)

    def test_ball_against_body_skin_transfers_equal_and_opposite_impulse(self):
        """The body's returned force ports must mirror what the ball receives."""
        prop = ball((0, 0, -.5))
        prop.v[:] = [0, 0, -2.]
        skin = [[0, 0, -.62], [.05, 0, -.62], [-.05, 0, -.62], [0, .05, -.62], [0, -.05, -.62]]
        owner = fixture(skin, rigid=[prop])
        exchanged, ports_seen = np.zeros(3), 0
        for _ in range(120):
            before = prop.mass * prop.v.copy()
            ports = owner.advance(.001, ENTITY)
            gravity_impulse = prop.mass * owner.gravity * .001
            delta = prop.mass * prop.v - before - gravity_impulse
            body = np.sum([p['force_n'] for p in ports], axis=0) * .001 if ports else np.zeros(3)
            if ports:
                ports_seen += 1
                # Equal and opposite to machine precision, substep by substep.
                np.testing.assert_allclose(delta, -body, atol=1e-12)
            exchanged += body
        self.assertGreater(ports_seen, 0)
        self.assertLess(exchanged[2], 0.)                  # ball pushes the skin down
        self.assertGreater(prop.v[2], -2.)                 # and is slowed by it


def main():
    """Write a receipt with explicit scope alongside the assertions."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite(loader.loadTestsFromTestCase(case)
                               for case in (Analytic, Conservation, Penalty))
    result = unittest.TextTestRunner(verbosity=2).run(suite)

    measured = {}
    for height in (.25, 1., 2., 5.):
        approach = np.sqrt(2 * G * height)
        _, v, _, _, _, _ = impact_step([0, 0, -approach], [0, 0, 0])
        measured[f'drop_{height}m'] = {
            'approach_m_s': approach, 'rebound_m_s': float(v[2]),
            'analytic_rebound_m_s': .75 * approach,
            'rebound_speed_ratio': float(v[2] / approach), 'analytic_ratio': .75,
            'apex_m': float(v[2] ** 2 / (2 * G)), 'analytic_apex_m': .75 ** 2 * height}
    oblique = {}
    for tangential in (.5, 3., 12., 20.):
        _, v, omega, _, _, _ = impact_step([tangential, 0, -3.], [0, 0, 0])
        vt, spin, sticking = analytic_impact(tangential, 3., 0.)
        oblique[f'tangential_{tangential}'] = {
            'sticking': bool(sticking), 'solver_vt_m_s': float(v[0]), 'analytic_vt_m_s': float(vt),
            'solver_spin_rad_s': float(omega[1]), 'analytic_spin_rad_s': float(spin),
            'abs_error_m_s': abs(float(v[0]) - float(vt))}

    # Rolling transition, measured rather than asserted.
    speed = 4.
    prop = ball((0, 0, -.9))
    prop.v[:] = [speed, 0, 0]
    owner = fixture([[0, 50, 0]], rigid=[prop])
    rolled_at = None
    for tick in range(2000):
        owner.advance(.001, ENTITY)
        if rolled_at is None and abs(prop.v[0] - prop.omega[1] * prop.radius) < 1e-9:
            rolled_at = (tick + 1) * .001
    rolling = {'initial_speed_m_s': speed, 'final_speed_m_s': float(prop.v[0]),
               'analytic_final_speed_m_s': 5 / 7 * speed,
               'roll_onset_s': rolled_at, 'analytic_roll_onset_s': 2 * speed / (7 * .35 * G),
               'analytic_slip_decay_rate_m_s2': SLIP_FACTOR * .35 * G,
               'basis': 'Classical uniform-sphere sliding-to-rolling transition'}

    # Penalty-collider characterization: how far a ball sinks into furniture.
    furniture = {}
    for drop in (.05, .2, .6):
        top = -.3
        table = (np.array([-.5, -.5, -1.]), np.array([.5, .5, top]), 'table')
        prop = ball((0, 0, top + .1 + drop))
        owner = fixture([[0, 50, 0]], rigid=[prop], static=[table])
        lowest = prop.x[2]
        for _ in range(6000):
            owner.advance(.001, ENTITY)
            lowest = min(lowest, prop.x[2])
        furniture[f'drop_{drop}m'] = {
            'impact_speed_m_s': float(np.sqrt(2 * G * drop)),
            'peak_penetration_m': float(max(0., top + prop.radius - lowest)),
            'linear_spring_estimate_m': float(np.sqrt(2 * G * drop) * np.sqrt(.4 / 1500.)),
            'rested_above_top': bool(prop.x[2] > top)}

    # Off-axis prop/prop strike: linear momentum exact, angular momentum not.
    drifts = []
    for substep in (.001, .0005, .00025, .000125):
        props, before, after = Penalty._offset_collision(substep)
        drifts.append({'substep_s': substep,
                       'linear_residual_kg_m_s': float(np.linalg.norm(after[0] - before[0])),
                       'angular_residual_kg_m2_s': float(np.linalg.norm(after[1] - before[1])),
                       'exchanged_spin_kg_m2_s': float(np.linalg.norm(sum(
                           .4 * p.mass * p.radius ** 2 * p.omega for p in props)))})
    for row in drifts:
        row['angular_residual_fraction'] = row['angular_residual_kg_m2_s'] / row['exchanged_spin_kg_m2_s']

    report = {
        'schema': 'ihm.ball-dynamics-validation.v1',
        'passed': result.wasSuccessful(),
        'tests_run': result.testsRun,
        'failures': [str(f[0]) for f in result.failures + result.errors],
        'material': dict(BOUNCY_BALL),
        'analytic_scope': (
            'Closed-form solutions of the DECLARED rigid-impulse model (Newton restitution '
            'plus clamped Coulomb tangential impulse, uniform solid sphere I=2/5 m r^2). '
            'Agreement demonstrates the solver implements its stated model exactly; it is '
            'NOT a measurement against a physical ball. Restitution .75 and friction .35 '
            'are uncalibrated engineering values.'),
        'penalty_scope': (
            'Ball against furniture, other props, cloth and body skin uses a compliant '
            'penalty law with velocity-regularized friction and, for spheres against '
            'boxes, only six axis-aligned surface probe points. No closed form exists and '
            'none is claimed: those cases are checked for equal/opposite exchange, '
            'non-increasing energy and absence of tunnelling only.'),
        'known_limitations': [
            'Sphere-vs-AABB furniture contact samples six probe points, so an edge or '
            'corner strike is approximated and a box thinner than the probe spacing '
            'could be missed.',
            'Only the ground plane uses the analytic time-of-impact solver; every other '
            'contact is compliant penalty, so a ball does not bounce elastically off '
            'skin, cloth or furniture.',
            'Restitution is velocity-independent above the .05 m/s bounce threshold; no '
            'measured coefficient, rolling resistance, spin decay or air drag is modelled.',
            'Furniture penalty stiffness is 1500 N/m, which is soft: a .4 kg ball '
            'arriving at 3.4 m/s sinks about 47 mm into a table before being pushed '
            'back out. See furniture_penetration.',
            'Total angular momentum is conserved only to first order in the environment '
            'substep. At the production 1 ms substep an off-axis prop/prop strike '
            'leaves a spurious net angular momentum of about 2.5% of the angular '
            'momentum exchanged. Linear momentum is exact to machine precision. '
            'See prop_momentum_residual.'],
        'normal_impacts': measured,
        'oblique_impacts': oblique,
        'rolling': rolling,
        'furniture_penetration': furniture,
        'prop_momentum_residual': drifts,
        'not_covered': [
            'Ball against the real 98-muscle native body is exercised only as '
            'integration survival in the unified world run; no analytic prediction '
            'exists for it and none is claimed here.',
            'Ball against cloth and against the volumetric pillow lattice.',
            'Repeated edge/corner strikes on furniture, and stacked props.'],
    }
    output = ROOT / 'data/derived/ball-dynamics-v1'
    output.mkdir(parents=True, exist_ok=True)
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'output': str(output), 'passed': report['passed'],
                      'tests_run': report['tests_run']}, indent=2))
    raise SystemExit(0 if report['passed'] else 1)


if __name__ == '__main__':
    main()
