"""Rigid coordinate bridge; canonical material registration is never modified.

Environment templates retain their own gravity-aligned coordinates. All world
points/loads crossing the native boundary must use this same immutable map.
Moments in force ports are free couples at the accompanying application point,
not moments about either frame's origin.
"""
from copy import deepcopy
import numpy as np


SUPINE_SOURCE_TO_WORLD = np.array([[0., 0., -1.], [0., 1., 0.], [1., 0., 0.]])


def _array(value, shape, name):
    result = np.array(value, dtype=float, copy=True)
    if result.shape != shape or not np.isfinite(result).all():
        raise ValueError(f'{name} must be finite with shape {shape}')
    return result


def _rotation(value, name):
    result = _array(value, (3, 3), name)
    if not np.allclose(result.T @ result, np.eye(3), atol=1e-10, rtol=0) or not np.isclose(np.linalg.det(result), 1., atol=1e-10, rtol=0):
        raise ValueError(f'{name} must be a proper orthonormal rotation')
    return result


class GravityAlignedWorldFrame:
    def __init__(self, global_map, *, source_to_world_rotation, source_gravity_m_s2,
                 source_support_point_m=None, source_support_normal=None,
                 world_support_point_m=None, support_identity='explicit'):
        c = _array(global_map, (4, 4), 'global_map')
        if not np.allclose(c[3], [0, 0, 0, 1], atol=1e-12, rtol=0):
            raise ValueError('global_map must be affine rigid')
        r = _rotation(c[:3, :3], 'canonical registration')
        p = _rotation(source_to_world_rotation, 'source_to_world_rotation')
        gravity = _array(source_gravity_m_s2, (3,), 'source_gravity_m_s2')
        supplied = [x is not None for x in (source_support_point_m, source_support_normal, world_support_point_m)]
        if any(supplied) and not all(supplied):
            raise ValueError('support point, normal and template point must be supplied together')
        q = r @ p.T
        shift = np.zeros(3)
        self.support = None
        if all(supplied):
            point = _array(source_support_point_m, (3,), 'source_support_point_m')
            normal = _array(source_support_normal, (3,), 'source_support_normal')
            wp = _array(world_support_point_m, (3,), 'world_support_point_m')
            if not np.isclose(np.linalg.norm(normal), 1., atol=1e-10):
                raise ValueError('support normal must be unit length')
            if np.linalg.norm(gravity) == 0 or not np.allclose(normal, -gravity / np.linalg.norm(gravity), atol=1e-10, rtol=0):
                raise ValueError('support normal must oppose gravity')
            # A plane fixes one translational degree of freedom. Preserve the
            # canonical tangential origin rather than treating its arbitrary
            # source plane point as the template's bed/floor center.
            canonical_normal = r @ normal
            displacement = r @ point + c[:3, 3] - q @ wp
            shift = canonical_normal * np.dot(canonical_normal, displacement)
            self.support = dict(identity=support_identity,
                                alignment_basis='minimum-normal-translation; tangential origin preserved',
                                source_point_m=point.tolist(),
                                source_normal=normal.tolist(), world_point_m=wp.tolist(),
                                world_normal=(p @ normal).tolist(),
                                canonical_point_m=(r @ point + c[:3, 3]).tolist(),
                                canonical_normal=(r @ normal).tolist())
        self.world_to_canonical = np.eye(4)
        self.world_to_canonical[:3, :3] = q
        self.world_to_canonical[:3, 3] = shift
        self.canonical_to_world = np.linalg.inv(self.world_to_canonical)
        self.source_to_canonical = c
        self.source_gravity_m_s2 = gravity
        self.world_gravity_m_s2 = p @ gravity
        self.canonical_gravity_m_s2 = r @ gravity
        for a in (self.world_to_canonical, self.canonical_to_world, self.source_to_canonical,
                  self.source_gravity_m_s2, self.world_gravity_m_s2, self.canonical_gravity_m_s2):
            a.setflags(write=False)

    @classmethod
    def identity(cls):
        """Explicit canonical-coordinate fixture with zero gravity and no support."""
        return cls(np.eye(4), source_to_world_rotation=np.eye(3),
                   source_gravity_m_s2=[0., 0., 0.], support_identity='explicit identity fixture')

    @classmethod
    def from_native(cls, registration, native_state, *, environment, template_support_height_m=None,
                    upright_support_point_source_m=None, upright_support_normal_source=None):
        """Use native gravity and support metadata; upright uses the source floor.

        The legacy native snapshot's support_plane_source_x_m is NOT the upright
        floor position. The retained source ContactHalfSpace named floor has
        location zero and orientation [0, 0, -pi/2]: source Y=0, upward +Y.
        """
        c = registration.global_map if hasattr(registration, 'global_map') else registration
        gravity = _array(native_state['gravity_m_s2'], (3,), 'native gravity')
        if environment == 'free':
            if np.linalg.norm(gravity) != 0:
                raise ValueError('free environment requires zero native gravity')
            return cls(c, source_to_world_rotation=np.asarray(c)[:3, :3], source_gravity_m_s2=gravity,
                       support_identity='free: no support')
        if environment == 'supine':
            return cls(c, source_to_world_rotation=SUPINE_SOURCE_TO_WORLD,
                       source_gravity_m_s2=gravity,
                       source_support_point_m=[native_state['support_plane_source_x_m'], 0., 0.],
                       source_support_normal=[1., 0., 0.],
                       world_support_point_m=[0., 0., -.24 if template_support_height_m is None else template_support_height_m],
                       support_identity='native supine support_plane_source_x_m')
        if environment == 'upright':
            if upright_support_point_source_m is None:
                upright_support_point_source_m = [0., 0., 0.]
            if upright_support_normal_source is None:
                upright_support_normal_source = [0., 1., 0.]
            return cls(c, source_to_world_rotation=np.eye(3), source_gravity_m_s2=gravity,
                       source_support_point_m=upright_support_point_source_m,
                       source_support_normal=upright_support_normal_source,
                       world_support_point_m=[0., -.96 if template_support_height_m is None else template_support_height_m, 0.],
                       support_identity='native source ContactHalfSpace floor: source Y=0 unless explicitly overridden')
        raise ValueError('unknown native environment: ' + str(environment))

    def _transform(self, values, to_world, point=False):
        x = np.asarray(values, dtype=float)
        if x.ndim == 0 or x.shape[-1] != 3 or not np.isfinite(x).all():
            raise ValueError('coordinates must be finite (..., 3) arrays')
        m = self.canonical_to_world if to_world else self.world_to_canonical
        return x @ m[:3, :3].T + (m[:3, 3] if point else 0.)

    def points_to_world(self, points): return self._transform(points, True, True)
    def points_to_canonical(self, points): return self._transform(points, False, True)
    def vectors_to_world(self, vectors): return self._transform(vectors, True)
    def vectors_to_canonical(self, vectors): return self._transform(vectors, False)

    def rotations_to_world(self, rotation):
        return self.canonical_to_world[:3, :3] @ _rotation(rotation, 'rotation')

    def rotations_to_canonical(self, rotation):
        return self.world_to_canonical[:3, :3] @ _rotation(rotation, 'rotation')

    def tensors_to_world(self, tensor):
        r = self.canonical_to_world[:3, :3]
        return r @ _array(tensor, (3, 3), 'tensor') @ r.T

    def tensors_to_canonical(self, tensor):
        r = self.world_to_canonical[:3, :3]
        return r @ _array(tensor, (3, 3), 'tensor') @ r.T

    def _port(self, port, to_world):
        result = deepcopy(port)
        result['point_m'] = self._transform(port['point_m'], to_world, True).tolist()
        result['force_n'] = self._transform(port['force_n'], to_world).tolist()
        for key in ('moment_nm', 'torque_nm', 'velocity_m_s', 'angular_velocity_rad_s'):
            if key in port:
                result[key] = self._transform(port[key], to_world).tolist()
        return result

    def force_port_to_world(self, port): return self._port(port, True)
    def force_port_to_canonical(self, port): return self._port(port, False)

    def metadata(self):
        return dict(world_to_canonical=self.world_to_canonical.tolist(),
                    canonical_to_world=self.canonical_to_world.tolist(),
                    source_to_canonical=self.source_to_canonical.tolist(),
                    source_gravity_m_s2=self.source_gravity_m_s2.tolist(),
                    world_gravity_m_s2=self.world_gravity_m_s2.tolist(),
                    canonical_gravity_m_s2=self.canonical_gravity_m_s2.tolist(),
                    support=deepcopy(self.support),
                    basis='Proper rigid coordinate change; no anatomy, mass or material registration change')
