#!/usr/bin/env python3
"""One measured gravity vector and one support plane, shared by body and world.

The native body's coordinates are not axis-aligned with the authored environment
templates: the registered canonical basis is tilted about 4.5 degrees away from
them. Before the common-frame correction the world integrated props and cloth
under an authored axis-aligned gravity while the body fell along its own measured
one, so "down" differed between them by that angle and the bed/floor planes did
not coincide. This verifier confirms, on the ACTUAL native body in both supported
environments, that a single measured gravity vector and a single support plane now
cross the boundary, and reports the residual numerically.

SCOPE. This is a coordinate-consistency proof. It shows the body and the world
agree about which way is down and where the support surface is, and that the map
between their frames is a proper rigid transform that preserves norms, angles,
power and momentum. It says nothing about whether either model is biomechanically
calibrated, and it is not a contact or balance validation.

Run: PYTHONPATH=. .venv/bin/python scripts/verify_common_frame.py
"""
import json
import tempfile
import time
from pathlib import Path

import numpy as np

from ihm.assembly.articulated import ArticulatedBodyPlant
from ihm.assembly.environment_dynamics import EnvironmentDynamics
from ihm.assembly.world_frame import GravityAlignedWorldFrame
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG

ROOT = Path(__file__).resolve().parents[1]
REGISTRATION = 'data/derived/mechanics/whole_body_lumbar_current/registration.json'
MASS_KG = MECHANICAL_TARGET_MASS_KG
# The gravity direction each environment template is authored around. The whole
# point of the correction is that the body does NOT share these axes natively.
AUTHORED = {'supine': ([0., 0., -1.], 2), 'upright': ([0., -1., 0.], 1)}
SCENES = {'supine': 'bedroom', 'upright': 'play-floor'}


def angle_between(a, b):
    a = np.asarray(a, float) / np.linalg.norm(a)
    b = np.asarray(b, float) / np.linalg.norm(b)
    return float(np.arccos(np.clip(a @ b, -1., 1.)))


def check(environment, output):
    authored, axis = AUTHORED[environment]
    plant = ArticulatedBodyPlant(ROOT, output / 'plant', environment=environment,
                                 target_mass_kg=MASS_KG, augmented_registration=REGISTRATION)
    try:
        state = plant.snapshot()
        native = plant.native.snapshot()
        registration = plant.registration
        owner = EnvironmentDynamics(ROOT, environment,
                                    {'scene': SCENES[environment], 'objects': []}, registration)
        frame = owner.world_frame
        report = {'environment': environment, 'scene': SCENES[environment],
                  'native_frame_bound': bool(owner.native_frame_bound)}
        assert owner.native_frame_bound

        # ---- A. One measured gravity vector -------------------------------
        source_gravity = np.asarray(native['gravity_m_s2'], float)
        canonical_gravity = registration.basis @ source_gravity
        # The body reports its own gravity in canonical coordinates; the frame
        # must have been built from that same measured vector, not an authored one.
        np.testing.assert_allclose(state['gravity_m_s2'], canonical_gravity, atol=1e-12)
        np.testing.assert_allclose(frame.canonical_gravity_m_s2, canonical_gravity, atol=1e-12)
        np.testing.assert_allclose(frame.source_gravity_m_s2, source_gravity, atol=1e-12)
        # What the world actually integrates props and cloth with.
        np.testing.assert_allclose(owner.gravity, frame.world_gravity_m_s2, atol=0., rtol=0.)
        # A rotation cannot change the magnitude of the measured field.
        magnitudes = [float(np.linalg.norm(v)) for v in
                      (source_gravity, canonical_gravity, frame.world_gravity_m_s2)]
        assert max(magnitudes) - min(magnitudes) < 1e-12, magnitudes
        # In the world frame gravity is exactly along the authored axis.
        np.testing.assert_allclose(frame.world_gravity_m_s2,
                                   np.array(authored) * magnitudes[0], atol=1e-12)
        off_axis = np.delete(frame.world_gravity_m_s2, axis)
        assert float(np.abs(off_axis).max()) < 1e-13, off_axis

        # The mismatch the correction removes, and what is left of it.
        uncorrected = angle_between(canonical_gravity, authored)
        corrected = angle_between(frame.world_gravity_m_s2, authored)
        report['gravity'] = {
            'source_m_s2': source_gravity.tolist(),
            'canonical_m_s2': canonical_gravity.tolist(),
            'world_m_s2': frame.world_gravity_m_s2.tolist(),
            'magnitude_m_s2': magnitudes,
            'magnitude_spread_m_s2': max(magnitudes) - min(magnitudes),
            'uncorrected_mismatch_rad': uncorrected,
            'uncorrected_mismatch_deg': np.degrees(uncorrected),
            'residual_after_correction_rad': corrected,
            'basis': ('Angle between the body\'s measured gravity in canonical '
                      'coordinates and the axis the environment template is authored '
                      'around. The first is the error a naive axis-aligned world would '
                      'carry; the second is what remains once the shared frame maps it.')}
        assert corrected < 1e-12, corrected
        assert .07 < uncorrected < .09, uncorrected      # the ~4.5 degrees on record

        # ---- B. One support plane ------------------------------------------
        support = frame.support
        assert support is not None
        world_normal = np.asarray(support['world_normal'], float)
        # The support opposes gravity in every frame, not just in the template.
        np.testing.assert_allclose(world_normal, -np.asarray(authored), atol=1e-12)
        np.testing.assert_allclose(np.asarray(support['source_normal'], float),
                                   -source_gravity / magnitudes[0], atol=1e-10)
        np.testing.assert_allclose(np.asarray(support['canonical_normal'], float),
                                   -canonical_gravity / magnitudes[0], atol=1e-10)
        # The world's own template height is where the body's plane lands.
        world_point = np.zeros(3)
        world_point[axis] = owner.base_plane
        canonical_point = frame.points_to_canonical(world_point)
        plane = {'template_height_m': owner.base_plane,
                 'object_floor_height_m': float(owner.object_plane),
                 'object_floor_is_template_plane':
                     bool(abs(owner.object_plane - owner.base_plane) < 1e-9)}

        if environment == 'supine':
            # The body publishes its own ideal plane; it must be the same plane.
            body_plane = state['body_environment']['plane']
            normal = np.asarray(body_plane['normal'], float)
            offset = float(abs((canonical_point - np.asarray(body_plane['point_m'], float)) @ normal))
            plane.update(source='body_environment.plane (native supine ideal plane)',
                         signed_offset_m=offset,
                         normal_angle_rad=angle_between(normal, -canonical_gravity))
            assert offset < 1e-10, offset
            np.testing.assert_allclose(frame.vectors_to_canonical(np.array(authored) * -1.),
                                       normal, atol=1e-12)
        else:
            # Upright has no published plane: measure it from the native contact
            # spheres the body is actually standing on. A loaded sphere must push
            # along +source Y, and its lowest point must sit at source Y = 0.
            loaded = [c for c in native['contacts']
                      if 'center_m' in c and np.linalg.norm(c['force_n']) > 1.]
            assert loaded, 'no loaded native contact spheres to measure the floor from'
            directions = [angle_between(c['force_n'], [0., 1., 0.]) for c in loaded]
            bottoms = [c['center_m'][1] - c['radius_m'] for c in loaded]
            weight = MASS_KG * magnitudes[0]
            total = float(sum(c['force_n'][1] for c in loaded))
            spread = float(max(bottoms) - min(bottoms))
            plane.update(source='native loaded contact spheres (measured, not declared)',
                         loaded_sphere_count=len(loaded),
                         max_force_direction_error_rad=max(directions),
                         lowest_loaded_sphere_source_y_m=min(bottoms),
                         highest_loaded_sphere_source_y_m=max(bottoms),
                         coplanarity_spread_m=spread,
                         total_normal_force_n=total,
                         total_normal_force_body_weights=total / weight,
                         initial_pose_penetration_m=float(-min(bottoms)),
                         measurement_basis=(
                             'The loaded spheres are coplanar to within '
                             f'{spread * 1e3:.1f} mm, which is the frame evidence: they '
                             'rest on one plane whose normal is the source floor normal '
                             '+Y, opposing the measured gravity. Their common depth '
                             'below source Y=0 is compliant-contact compression of the '
                             'INITIAL POSE at t=0, which is not in contact equilibrium; '
                             'that transient is a native initial-condition property, '
                             'not a frame disagreement.'))
            # Every loaded contact pushes the body up along the source floor normal.
            assert max(directions) < 1e-6, max(directions)
            # They lie on one plane; that is what "a shared support plane" means here.
            assert spread < 5e-3, spread
            # All of them are below the declared plane, i.e. loaded by penetration.
            assert max(bottoms) < 0., max(bottoms)
            assert min(bottoms) > -.2, min(bottoms)
            np.testing.assert_allclose(support['source_point_m'], [0., 0., 0.], atol=1e-12)
            np.testing.assert_allclose(support['source_normal'], [0., 1., 0.], atol=1e-12)
        report['support_plane'] = plane

        # ---- C. The map is a proper rigid transform ------------------------
        rotation = frame.world_to_canonical[:3, :3]
        np.testing.assert_allclose(rotation.T @ rotation, np.eye(3), atol=1e-12)
        determinant = float(np.linalg.det(rotation))
        assert abs(determinant - 1.) < 1e-10, determinant
        skin = owner.skin_points(state['entities'])           # already world-frame
        back = frame.points_to_canonical(skin)
        round_trip = float(np.abs(frame.points_to_world(back) - skin).max())
        assert round_trip < 1e-12, round_trip
        rng = np.random.default_rng(908)
        force, velocity = rng.normal(size=(2, 64, 3))
        power = float(np.abs(np.sum(force * velocity, axis=1)
                             - np.sum(frame.vectors_to_world(force)
                                      * frame.vectors_to_world(velocity), axis=1)).max())
        assert power < 1e-12, power
        distances = lambda p: np.linalg.norm(p[1:] - p[:-1], axis=1)
        length = float(np.abs(distances(skin) - distances(back)).max())
        assert length < 1e-12, length
        report['transform'] = {'determinant': determinant,
                               'point_round_trip_m': round_trip,
                               'preserved_pair_distance_error_m': length,
                               'force_port_power_error_w': power,
                               'skin_sample_count': int(len(skin))}

        # ---- D. Live exchange keeps using that one frame -------------------
        ports = owner.advance(.005, state['entities'])
        impulse = np.sum([p['force_n'] for p in ports], axis=0) * .005 if ports else np.zeros(3)
        np.testing.assert_allclose(impulse, owner.last_impulse, atol=1e-9)
        exported = owner.frame()
        np.testing.assert_allclose(exported['world_frame']['world_gravity_m_s2'],
                                   frame.world_gravity_m_s2, atol=0., rtol=0.)
        assert exported['native_frame_bound'] is True
        assert exported['body_force_coordinate_frame'] == 'canonical'
        assert exported['object_coordinate_frame'] == 'gravity_aligned_world'
        report['exchange'] = {'body_force_ports': len(ports),
                              'body_impulse_ns': owner.last_impulse.tolist(),
                              'exported_frames': [exported['body_force_coordinate_frame'],
                                                  exported['object_coordinate_frame']]}
        return report
    finally:
        plant.close()


def synthetic_regression():
    """A deliberately tilted registration must still land on one shared frame."""
    theta = .12
    # Tilt about X, so the registered "up" genuinely leaves the template axis.
    tilt = np.array([[1., 0., 0.], [0., np.cos(theta), -np.sin(theta)],
                     [0., np.sin(theta), np.cos(theta)]])
    transform = np.eye(4)
    transform[:3, :3] = tilt
    transform[:3, 3] = [.4, -.2, .7]
    frame = GravityAlignedWorldFrame.from_native(
        transform, {'gravity_m_s2': [0., -9.81, 0.]}, environment='upright',
        template_support_height_m=-.96)
    np.testing.assert_allclose(frame.world_gravity_m_s2, [0., -9.81, 0.], atol=1e-12)
    # The tilt must live entirely in the map, never in the magnitude.
    assert abs(np.linalg.norm(frame.canonical_gravity_m_s2) - 9.81) < 1e-12
    assert angle_between(frame.canonical_gravity_m_s2, [0., -1., 0.]) > .1
    np.testing.assert_allclose(frame.points_to_world(transform[:3, 3])[1], -.96, atol=1e-12)
    # A gravity that does not oppose the declared support must be refused.
    for bad in ([0., 9.81, 0.], [9.81, 0., 0.]):
        try:
            GravityAlignedWorldFrame(transform, source_to_world_rotation=np.eye(3),
                                     source_gravity_m_s2=bad,
                                     source_support_point_m=[0., 0., 0.],
                                     source_support_normal=[0., 1., 0.],
                                     world_support_point_m=[0., -.96, 0.])
        except ValueError:
            continue
        raise AssertionError('support normal inconsistent with gravity was accepted')
    return {'tilt_rad': theta, 'refused_inconsistent_gravity': True}


def main():
    started = time.monotonic()
    output = Path(tempfile.mkdtemp(prefix='common-frame-', dir=ROOT / 'data/derived'))
    report = {'schema': 'ihm.common-frame-consistency.v1', 'passed': False,
              'scope': ('Coordinate consistency between the actual native body and the '
                        'server-owned world, in both supported environments. Proves one '
                        'measured gravity vector and one support plane cross the '
                        'boundary and that the map is a proper rigid transform. Not a '
                        'contact, balance or biomechanical validation.'),
              'environments': {}}
    try:
        report['synthetic_regression'] = synthetic_regression()
        for environment in ('supine', 'upright'):
            report['environments'][environment] = check(environment, output / environment)
            print(json.dumps({'event': 'checked', 'environment': environment}), flush=True)
        report['passed'] = True
        report['summary'] = {
            'uncorrected_mismatch_deg': {e: v['gravity']['uncorrected_mismatch_deg']
                                         for e, v in report['environments'].items()},
            'residual_after_correction_rad': {e: v['gravity']['residual_after_correction_rad']
                                              for e, v in report['environments'].items()},
            'gravity_magnitude_spread_m_s2': {e: v['gravity']['magnitude_spread_m_s2']
                                              for e, v in report['environments'].items()}}
    except BaseException as error:
        report['error'] = f'{type(error).__name__}: {error}'
        raise
    finally:
        report['wall_s'] = time.monotonic() - started
        (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps({'output': str(output), 'passed': report['passed'],
                          **report.get('summary', {})}, indent=2), flush=True)
    raise SystemExit(0 if report['passed'] else 1)


if __name__ == '__main__':
    main()
