"""Per-body, per-axis geometric scaling of the mechanical model.

``ihm.native.model_scaling`` scales the whole body by one factor.  That is what
``stature_m`` is, and it is exact: path length is homogeneous of degree one in
the geometry at fixed pose, so the fitted path polynomials can be scaled by a
single multiplication.  It changes how BIG the subject is and not what SHAPE it
is.

Shape is what this module changes, and it is where the sex-dimorphic
proportions live.  NHANES gives women a hip circumference 12.6% larger relative
to stature than men's (``d = -0.91``) and a waist:hip ratio 8% smaller
(``d = 1.07``), while limb-to-trunk proportions differ by under 3%.  Those are
per-body, per-axis changes.

**Nothing here may be combined with a coefficient scale of the fitted paths.**
An anisotropic change moves 56 of 98 muscle paths by between 0.01% and 6.23% and
leaves 42 alone; a single factor is wrong for every muscle but at most one.
Every model this module produces must have its paths REFITTED --
``ihm.native.path_refitting`` -- and ``anisotropic_scale`` is written into the
report so a caller that skips the refit is making a recorded mistake rather than
a silent one.

What is exact here, and what is not
-----------------------------------
Exact: every point attached to a scaled body -- path points, markers, mass
centre, the offset frames that carry child joints -- moves componentwise in that
body's own frame.  Mass and the full inertia tensor follow analytically from the
second moments, and reduce to ``s**3`` and ``s**5`` under isotropy.

Not exact: **wrap object shape**.  Every wrap object in this model is a
``WrapCylinder``, most of them rotated in their body's frame.  A rotated
cylinder under an anisotropic scale is an elliptic cylinder, and OpenSim has no
such wrap.  ``wrap_shape='axis_aware'`` scales its length by the stretch along
its own axis and its radius by the geometric mean of the two principal
stretches perpendicular to that axis -- which preserves cross-sectional area and
is exactly right under isotropy.  ``wrap_shape='translate_only'`` moves wrap
objects without resizing them.  The difference between the two policies is
measurable and is measured; it is not assumed away.
"""
import math
import xml.etree.ElementTree as ET

import numpy as np

from .model_scaling import _scale_function_output

WRAP_SHAPE_POLICIES = ('axis_aware', 'translate_only')


def _vec(text):
    values = [float(x) for x in (text or '').split()]
    if not values:
        raise ValueError('Empty vector')
    return np.array(values, float)


def _set(node, values):
    node.text = ' '.join(repr(float(v)) for v in values)


def _body_rotation(xyz):
    """OpenSim body-fixed XYZ Euler angles, as used by WrapObject."""
    a, b, c = xyz
    rx = np.array([[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]])
    ry = np.array([[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]])
    rz = np.array([[math.cos(c), -math.sin(c), 0], [math.sin(c), math.cos(c), 0], [0, 0, 1]])
    return rx @ ry @ rz


def scale_inertia(mass, inertia, factors):
    """Exact inertia of a body scaled by ``factors`` at constant density.

    ``inertia`` is OpenSim's (Ixx Iyy Izz Ixy Ixz Iyz) about the mass centre.
    The diagonal recovers the three second moments ``Jx = integral x^2 dm`` etc.,
    each of which scales by ``volume * factor**2`` along its own axis; the
    products scale by ``volume * fi * fj``.  Under isotropy this returns
    ``s**5`` times the input, which is what ``model_scaling`` does by hand.

    A body whose declared inertia violates the triangle inequality has no
    consistent mass distribution and raises rather than being scaled into a
    plausible-looking wrong answer.
    """
    ixx, iyy, izz, ixy, ixz, iyz = [float(v) for v in inertia]
    fx, fy, fz = [float(f) for f in factors]
    volume = fx * fy * fz
    jx = 0.5 * (iyy + izz - ixx)
    jy = 0.5 * (ixx + izz - iyy)
    jz = 0.5 * (ixx + iyy - izz)
    if min(jx, jy, jz) < -1e-12 * max(abs(ixx), abs(iyy), abs(izz), 1e-30):
        raise ValueError('Inertia (%r) violates the triangle inequality; it has no '
                         'consistent mass distribution to scale' % (inertia,))
    jx, jy, jz = volume * fx * fx * jx, volume * fy * fy * jy, volume * fz * fz * jz
    return (jy + jz, jx + jz, jx + jy,
            volume * fx * fy * ixy, volume * fx * fz * ixz, volume * fy * fz * iyz)


def scale_model_anisotropic(model_bytes, factors, *, wrap_shape='axis_aware',
                            scale_mass=True):
    """Scale named bodies by per-axis factors; return (bytes, report).

    ``factors`` maps body name -> ``(fx, fy, fz)`` in that body's own frame.
    Bodies not named are untouched, so the femur stays the source subject's
    femur while the pelvis widens -- which is exactly the pelvic-breadth
    parameter, and is why the hip joint centres separate.

    ``scale_mass`` also scales the body's mass and inertia.  The native engine
    renormalises total mass to ``target_mass_kg`` afterwards, so what this
    actually changes is the mass DISTRIBUTION between segments -- which is the
    point, and is the one thing the existing ``mass_kg`` knob cannot do.
    """
    if wrap_shape not in WRAP_SHAPE_POLICIES:
        raise ValueError('wrap_shape must be one of %r' % (WRAP_SHAPE_POLICIES,))
    factors = {name: np.array([float(f) for f in triple], float)
               for name, triple in factors.items()}
    for name, triple in factors.items():
        if triple.shape != (3,) or not np.isfinite(triple).all() or (triple <= 0).any():
            raise ValueError('Body %s needs three positive finite factors' % name)

    root = ET.fromstring(model_bytes)
    bodies = {body.get('name'): body for body in root.iter('Body')}
    unknown = set(factors) - set(bodies)
    if unknown:
        raise ValueError('Model has no bodies named %r' % sorted(unknown))
    if not factors:
        raise ValueError('No bodies to scale; use model_scaling.scale_model for isotropy')

    report = {'factors': {k: v.tolist() for k, v in factors.items()},
              'wrap_shape_policy': wrap_shape, 'elements': {}, 'wraps': [],
              'mass_kg': {}, 'anisotropic_scale': True,
              'requires_path_refit': True}

    def bump(key, n=1):
        report['elements'][key] = report['elements'].get(key, 0) + n

    # 1. Everything attached directly to a scaled body's frame.
    for element in root.iter():
        socket = (element.findtext('socket_parent_frame')
                  or element.findtext('socket_parent'))
        if not socket:
            continue
        body = socket.rsplit('/', 1)[-1]
        if socket.rsplit('/', 1)[0] != '/bodyset' or body not in factors:
            continue
        for tag in ('location', 'translation'):
            node = element.find(tag)
            if node is not None and node.text and node.text.split():
                _set(node, _vec(node.text) * factors[body])
                bump(element.tag + '.' + tag)

    # 2. The bodies themselves: mass centre, attached display geometry, mass,
    #    inertia, and the wrap objects they carry.
    for name, factor in factors.items():
        body = bodies[name]
        centre = body.find('mass_center')
        if centre is not None and centre.text and centre.text.split():
            _set(centre, _vec(centre.text) * factor)
            bump('Body.mass_center')
        for node in body.iter('scale_factors'):
            if node.text and node.text.split():
                _set(node, _vec(node.text) * factor)
                bump('Geometry.scale_factors')
        if scale_mass:
            mass = float(body.findtext('mass'))
            volume = float(np.prod(factor))
            body.find('mass').text = repr(mass * volume)
            inertia = body.find('inertia')
            _set(inertia, scale_inertia(mass, _vec(inertia.text), factor))
            report['mass_kg'][name] = {'before': mass, 'after': mass * volume}
            bump('Body.mass')

        wraps = body.find('WrapObjectSet/objects')
        for wrap in (wraps if wraps is not None else []):
            translation = wrap.find('translation')
            if translation is not None and translation.text and translation.text.split():
                _set(translation, _vec(translation.text) * factor)
            entry = {'body': name, 'wrap': wrap.get('name'), 'type': wrap.tag,
                     'length_factor': 1.0, 'radius_factor': 1.0}
            if wrap_shape == 'axis_aware':
                rotation = _body_rotation(_vec(wrap.findtext('xyz_body_rotation') or '0 0 0'))
                if wrap.tag == 'WrapCylinder':
                    # OpenSim's WrapCylinder axis is its local z.
                    axis = rotation[:, 2]
                    perpendicular = (rotation[:, 0], rotation[:, 1])
                elif wrap.tag == 'WrapSphere':
                    axis, perpendicular = None, (rotation[:, 0], rotation[:, 1])
                else:
                    raise ValueError('Wrap type %s has no declared anisotropic '
                                     'scaling rule' % wrap.tag)
                stretch = lambda d: float(np.linalg.norm(factor * d))  # noqa: E731
                radius_factor = math.sqrt(stretch(perpendicular[0]) * stretch(perpendicular[1]))
                node = wrap.find('radius')
                if node is not None:
                    node.text = repr(float(node.text) * radius_factor)
                entry['radius_factor'] = radius_factor
                if axis is not None:
                    length_factor = stretch(axis)
                    node = wrap.find('length')
                    if node is not None:
                        node.text = repr(float(node.text) * length_factor)
                    entry['length_factor'] = length_factor
            report['wraps'].append(entry)
            bump('WrapObject')

    # 3. CustomJoint translation transform functions.  These output a LENGTH in
    #    the joint's parent offset frame -- the walker knee's three translations
    #    are the tibiofemoral roll-glide -- so they move with a scaled parent.
    #    Componentwise scaling is only defined when the offset frame is aligned
    #    with the body it sits in; ``walker_knee_r``'s is rotated by
    #    (-1.64, 1.45, 1.57) rad, so an anisotropic femur raises here rather than
    #    silently producing a knee that translates along the wrong axes.
    for joint in root.iter('CustomJoint'):
        transform = joint.find('SpatialTransform')
        if transform is None:
            continue
        socket = joint.findtext('socket_parent_frame') or ''
        frames = {frame.get('name'): frame for frame in (joint.find('frames') or [])}
        frame = frames.get(socket.rsplit('/', 1)[-1])
        if frame is None:
            continue
        parent = (frame.findtext('socket_parent') or '')
        body = parent.rsplit('/', 1)[-1]
        if parent.rsplit('/', 1)[0] != '/bodyset' or body not in factors:
            continue
        factor = factors[body]
        isotropic = float(factor.max() - factor.min()) < 1e-12
        aligned = float(np.abs(_vec(frame.findtext('orientation') or '0 0 0')).max()) < 1e-12
        for axis in transform:
            if not (axis.get('name') or '').startswith('translation'):
                continue
            direction = _vec(axis.findtext('axis'))
            if isotropic:
                component = float(factor[0])
            else:
                nonzero = np.flatnonzero(np.abs(direction) > 1e-12)
                if not aligned or nonzero.size != 1:
                    raise ValueError(
                        'Joint %s takes its translations in a frame that is '
                        'rotated in body %s, or along a non-principal axis %r. '
                        'An anisotropic scale of that body has no componentwise '
                        'meaning here.' % (joint.get('name'), body, direction.tolist()))
                component = float(factor[nonzero[0]])
            for child in axis:
                if child.tag in ('coordinates', 'axis'):
                    continue
                _scale_function_output(child, component)
                bump('translation_axis')

    report['limitation'] = (
        'Wrap objects: every one in this model is a rotated WrapCylinder, and a '
        'rotated cylinder under an anisotropic scale is an elliptic cylinder, '
        'which OpenSim cannot express. Policy %r was applied and the per-wrap '
        'factors are listed. Muscle paths are affected only where a wrap is '
        'active, and the size of that effect is measured by re-running with '
        "policy 'translate_only'." % wrap_shape)
    return ET.tostring(root, encoding='utf-8', xml_declaration=True), report
