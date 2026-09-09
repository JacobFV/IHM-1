#!/usr/bin/env python3
"""Coordinate invariance checks independent of an OpenSim installation."""
import copy
import numpy as np
from ihm.assembly.world_frame import GravityAlignedWorldFrame, SUPINE_SOURCE_TO_WORLD


def main():
    theta = .079
    tilt = np.array([[1., 0., 0.], [0., np.cos(theta), -np.sin(theta)], [0., np.sin(theta), np.cos(theta)]])
    c = np.eye(4)
    c[:3, :3] = tilt @ SUPINE_SOURCE_TO_WORLD
    c[:3, 3] = [.11, -.03, .19]
    before = c.copy()
    snapshot = dict(gravity_m_s2=[-9.81, 0., 0.], support_plane_source_x_m=-.137)
    frame = GravityAlignedWorldFrame.from_native(c, snapshot, environment='supine')
    q = frame.world_to_canonical[:3, :3]
    assert np.isclose(np.linalg.det(q), 1.)
    np.testing.assert_allclose(frame.world_gravity_m_s2, [0., 0., -9.81])
    canonical_plane = c[:3, :3] @ [-.137, 0., 0.] + c[:3, 3]
    np.testing.assert_allclose(frame.points_to_world(canonical_plane)[2], -.24, atol=1e-14)
    mapped_template = frame.points_to_canonical([0., 0., -.24])
    np.testing.assert_allclose(np.dot(mapped_template - canonical_plane, c[:3, 0]), 0., atol=1e-14)
    shift = frame.world_to_canonical[:3, 3]
    np.testing.assert_allclose(shift - np.dot(shift, c[:3, 0]) * c[:3, 0], 0., atol=1e-14)
    # Arbitrary plane anchor tangential coordinates cannot relocate the bed.
    alternate = GravityAlignedWorldFrame(c, source_to_world_rotation=SUPINE_SOURCE_TO_WORLD,
        source_gravity_m_s2=[-9.81,0.,0.], source_support_point_m=[-.137,40.,-50.],
        source_support_normal=[1.,0.,0.], world_support_point_m=[0.,0.,-.24])
    np.testing.assert_allclose(alternate.world_to_canonical, frame.world_to_canonical, atol=1e-14)
    np.testing.assert_allclose(frame.vectors_to_world(c[:3, 0]), [0., 0., 1.], atol=1e-14)
    rng = np.random.default_rng(71)
    x, v, f = rng.normal(size=(3, 23, 3))
    np.testing.assert_allclose(frame.points_to_canonical(frame.points_to_world(x)), x, atol=1e-14)
    np.testing.assert_allclose(frame.vectors_to_canonical(frame.vectors_to_world(v)), v, atol=1e-14)
    np.testing.assert_allclose(np.sum(f*v, axis=1), np.sum(frame.vectors_to_world(f)*frame.vectors_to_world(v), axis=1), atol=1e-14)
    np.testing.assert_allclose(np.sum(frame.vectors_to_world(f), axis=0), frame.vectors_to_world(np.sum(f, axis=0)), atol=1e-14)
    # Angular momentum about one shared physical origin also transforms rigidly.
    origin = [.3, -.5, .6]
    angular = np.cross(x-origin, v)
    world_angular = np.cross(frame.points_to_world(x)-frame.points_to_world(origin), frame.vectors_to_world(v))
    np.testing.assert_allclose(world_angular, frame.vectors_to_world(angular), atol=1e-14)
    port = dict(id='skin', point_m=x[0].tolist(), force_n=f[0].tolist(), moment_nm=[.1,.2,.3], velocity_m_s=v[0].tolist(), angular_velocity_rad_s=[.4,.2,.1])
    original = copy.deepcopy(port)
    wp = frame.force_port_to_world(port)
    restored = frame.force_port_to_canonical(wp)
    for key in port:
        if key != 'id': np.testing.assert_allclose(restored[key], port[key], atol=1e-14)
    power = lambda p: np.dot(p['force_n'],p['velocity_m_s']) + np.dot(p['moment_nm'],p['angular_velocity_rad_s'])
    np.testing.assert_allclose(power(wp), power(port), atol=1e-14)
    assert original == port
    np.testing.assert_array_equal(c, before)
    np.testing.assert_allclose(frame.rotations_to_canonical(frame.rotations_to_world(tilt)), tilt, atol=1e-14)
    tensor = np.diag([1., 2., 3.])
    np.testing.assert_allclose(frame.tensors_to_canonical(frame.tensors_to_world(tensor)), tensor, atol=1e-14)
    upright = GravityAlignedWorldFrame.from_native(c, dict(gravity_m_s2=[0.,-9.81,0.]), environment='upright')
    np.testing.assert_allclose(upright.points_to_world(c[:3,3])[1], -.96, atol=1e-14)
    np.testing.assert_allclose(np.dot(upright.points_to_canonical([0.,-.96,0.])-c[:3,3], c[:3,1]), 0., atol=1e-14)
    upshift = upright.world_to_canonical[:3,3]
    np.testing.assert_allclose(upshift - np.dot(upshift,c[:3,1])*c[:3,1], 0., atol=1e-14)
    np.testing.assert_allclose(upright.vectors_to_world(c[:3,1]), [0.,1.,0.], atol=1e-14)
    free = GravityAlignedWorldFrame.from_native(c, dict(gravity_m_s2=[0.,0.,0.]), environment='free')
    np.testing.assert_allclose(free.world_to_canonical, np.eye(4), atol=1e-14)
    for bad in (np.diag([-1.,1.,1.,1.]), np.diag([2.,1.,1.,1.])):
        try: GravityAlignedWorldFrame.from_native(bad, snapshot, environment='supine')
        except ValueError: pass
        else: raise AssertionError('improper registration accepted')
    try: GravityAlignedWorldFrame.from_native(c, dict(snapshot,gravity_m_s2=[0.,-9.81,0.]), environment='supine')
    except ValueError: pass
    else: raise AssertionError('gravity/support mismatch accepted')
    metadata = frame.metadata()
    metadata['support']['source_point_m'][0] = 100.
    assert frame.metadata()['support']['source_point_m'][0] == -.137
    print('PASS: support/gravity alignment; point, vector, tensor and rotation roundtrip; force-port power, linear/angular momentum; immutable registration; upright and free modes; invalid frames rejected')


if __name__ == '__main__':
    main()
