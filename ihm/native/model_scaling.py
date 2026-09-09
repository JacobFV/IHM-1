"""Measure the mechanical body, and scale it geometrically.

Two things live here.

**Measurement.**  ``head_marker_height_m`` needs a forward-kinematic pose of the
whole model at its default coordinates, which is why a small FK walk of the
joint tree is in this file rather than borrowed.  The walk evaluates each
joint's transform functions at the declared ``default_value`` of their
coordinates rather than assuming the joint sits at identity: the first draft did
assume it, and ``patellofemoral_r`` -- whose ``rotation1`` is 1.14e-3 rad at a
zero knee angle -- is the case that says why not.

**Scaling.**  ``scale_model`` applies one isotropic linear factor to the model.
The trap it exists to avoid is stated in ``docs/BODY_PARAMETERS.md`` and is worth
repeating here: the muscle paths this model actually runs on are NOT the
``GeometryPath`` point sets in the ``.osim``.  ``native_mechanical_stream.cpp``
calls ``ModelFactory::replacePathsWithFunctionBasedPaths`` and substitutes 80
fitted multivariate polynomials from a *separate file*.  Scale the model and not
that file and every muscle in the lower limb keeps the unscaled subject's
length while pulling on a resized skeleton -- silently, with no error and no
warning.  ``scale_function_based_path_set`` is therefore not optional, and
``scaled_source_files`` returns the model and the path set together so a caller
cannot take one without the other.
"""
from pathlib import Path
import math
import xml.etree.ElementTree as ET

import numpy as np

#: Markers AddBiomechanics places on the floor plane.  Their common height is
#: the ground the subject stood on.
GROUND_MARKERS = ('R.HeelGround', 'L.HeelGround', 'R.MT5Ground', 'L.MT5Ground',
                  'R.ToeGround', 'L.ToeGround')
HEAD_MARKER = 'Head'


def _vec3(text):
    v = np.array([float(x) for x in (text or '0 0 0').split()], float)
    if v.shape != (3,) or not np.isfinite(v).all():
        raise ValueError('Three finite components required')
    return v


def _rotation(xyz):
    """OpenSim body-fixed XYZ Euler angles."""
    a, b, c = xyz
    rx = np.array([[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]])
    ry = np.array([[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]])
    rz = np.array([[math.cos(c), -math.sin(c), 0], [math.sin(c), math.cos(c), 0], [0, 0, 1]])
    return rx @ ry @ rz


def _function_value(node, q):
    """Value of an OpenSim scalar joint-transform function at argument ``q``.

    Only the forms this model actually uses are supported, and an unknown form
    raises rather than defaulting to zero -- a silent zero here would turn an
    unmodelled joint offset into a stature error.
    """
    tag = node.tag
    if tag == 'MultiplierFunction':
        inner = [c for c in node if c.tag == 'function']
        if len(inner) != 1 or len(inner[0]) != 1:
            raise ValueError('MultiplierFunction needs exactly one inner function')
        return float(node.findtext('scale', '1')) * _function_value(inner[0][0], q)
    if tag == 'Constant':
        return float(node.findtext('value'))
    if tag in ('PolynomialFunction', 'LinearFunction'):
        # OpenSim writes these highest-degree first.
        return float(np.polyval([float(c) for c in node.findtext('coefficients').split()], q))
    if tag == 'SimmSpline':
        x = [float(v) for v in node.findtext('x').split()]
        y = [float(v) for v in node.findtext('y').split()]
        if len(x) != len(y) or not x:
            raise ValueError('Malformed SimmSpline')
        return float(np.interp(q, x, y))
    raise ValueError('Unsupported joint transform function: ' + tag)


def _axis_value(axis, defaults):
    """(value, axis vector) of one TransformAxis at the model's default pose."""
    names = (axis.findtext('coordinates') or '').split()
    if len(names) > 1:
        raise ValueError('Multi-coordinate transform axes are not supported')
    q = defaults.get(names[0], 0.0) if names else 0.0
    functions = [c for c in axis if c.tag not in ('coordinates', 'axis')]
    if len(functions) > 1:
        raise ValueError('More than one function on a transform axis')
    value = _function_value(functions[0], q) if functions else 0.0
    return value, _vec3(axis.findtext('axis'))


def _frame(name, frames):
    """Resolve a socket path to (body, translation, rotation) in that body."""
    local = name.split('/')[-1]
    if local in frames:
        translation, orientation, parent = frames[local]
        body = 'ground' if parent == '/ground' else parent.split('/')[-1]
        return body, translation, _rotation(orientation)
    return ('ground' if name == '/ground' else local), np.zeros(3), np.eye(3)


def _joint_transform(joint, defaults):
    """(rotation, translation) from a joint's parent frame to its child frame."""
    tag = joint.tag
    if tag == 'WeldJoint':
        return np.eye(3), np.zeros(3)
    if tag == 'PinJoint':
        names = [c.get('name') for c in joint.iter('Coordinate')]
        angle = defaults.get(names[0], 0.0) if names else 0.0
        return _rotation((0.0, 0.0, angle)), np.zeros(3)
    if tag != 'CustomJoint':
        raise ValueError('Unsupported joint type: ' + tag)
    transform = joint.find('SpatialTransform')
    rotation, translation = np.eye(3), np.zeros(3)
    for axis in transform:
        name = axis.get('name') or ''
        value, direction = _axis_value(axis, defaults)
        if name.startswith('rotation'):
            # Rotations compose in declaration order about their own body-fixed
            # axes, which is what SimTK's SpatialTransform does.
            magnitude = np.linalg.norm(direction)
            if magnitude == 0:
                raise ValueError('Zero rotation axis')
            unit = direction / magnitude
            skew = np.array([[0, -unit[2], unit[1]], [unit[2], 0, -unit[0]], [-unit[1], unit[0], 0]])
            rotation = rotation @ (np.eye(3) + math.sin(value) * skew
                                   + (1 - math.cos(value)) * (skew @ skew))
        elif name.startswith('translation'):
            translation = translation + value * direction
        else:
            raise ValueError('Unknown transform axis ' + name)
    return rotation, translation


def default_pose_frames(model_path):
    """World transform of every body at the model's default coordinates.

    Returns ``{body: (origin_m, rotation)}``.
    """
    root = ET.parse(model_path).getroot()
    defaults = {c.get('name'): float(c.findtext('default_value'))
                for c in root.iter('Coordinate')}
    edges = {}
    for joint in root.find('.//JointSet/objects'):
        frames = {}
        container = joint.find('frames')
        if container is not None:
            for offset in container:
                frames[offset.get('name')] = (_vec3(offset.findtext('translation')),
                                              _vec3(offset.findtext('orientation')),
                                              offset.findtext('socket_parent'))
        rj, tj = _joint_transform(joint, defaults)
        parent_body, tp, rp = _frame(joint.findtext('socket_parent_frame'), frames)
        child_body, tc, rc = _frame(joint.findtext('socket_child_frame'), frames)
        # parent body -> parent offset -> joint -> child offset -> child body
        rotation = rp @ rj @ rc.T
        edges[child_body] = (parent_body, tp + rp @ tj - rotation @ tc, rotation)

    pose = {'ground': (np.zeros(3), np.eye(3))}

    def place(body, seen=()):
        if body in pose:
            return pose[body]
        if body in seen:
            raise ValueError('Cycle in the joint tree at ' + body)
        parent, translation, rotation = edges[body]
        origin, frame = place(parent, seen + (body,))
        pose[body] = (origin + frame @ translation, frame @ rotation)
        return pose[body]

    for body in list(edges):
        place(body)
    return pose


def marker_positions(model_path, pose=None):
    root = ET.parse(model_path).getroot()
    pose = pose or default_pose_frames(model_path)
    out = {}
    for marker in root.iter('Marker'):
        body = marker.findtext('socket_parent_frame').split('/')[-1]
        origin, frame = pose[body]
        out[marker.get('name')] = origin + frame @ _vec3(marker.findtext('location'))
    return out


def head_marker_height_m(model_path):
    """Stature proxy: floor-marker plane to Head marker, at the default pose.

    Why a marker pair and not a mesh bounding box.  The model's display geometry
    is a *generic* mesh plus per-body ``scale_factors``; its bounding box is
    therefore a statement about the generic donor, refracted through nine
    anisotropic scale triples, and the head mesh in particular rides ``torso``
    with no neck joint.  The ``*Ground`` markers, by contrast, were placed by
    AddBiomechanics on the plane the subject actually stood on, and ``Head`` on
    the subject's own head.  Both are properties of this subject.

    Both are also markers, so this is a proxy and not a standing height: the
    ``Head`` marker sits above the scalp.  Carried onto the anatomical atlas by
    the fitted registration scale it overshoots the atlas skin extent by 0.66%,
    which bounds the error from above.
    """
    markers = marker_positions(model_path)
    heights = [markers[name][1] for name in GROUND_MARKERS if name in markers]
    if len(heights) < 4:
        raise ValueError('Model does not carry the AddBiomechanics ground markers')
    # 1 um. The six markers agree to 5.5e-9 m in the source model; the check is
    # here to catch a model whose feet are not both on the floor at its default
    # pose, which would make the measurement below meaningless.
    if max(heights) - min(heights) > 1e-6:
        raise ValueError('Ground markers are not coplanar in height: %r' % (heights,))
    if HEAD_MARKER not in markers:
        raise ValueError('Model does not carry a Head marker')
    return float(markers[HEAD_MARKER][1] - heights[0])


# ---------------------------------------------------------------------------
# Scaling
# ---------------------------------------------------------------------------

#: Every element whose text is a length in metres and which therefore takes one
#: factor of ``scale``.  Read as: tag -> (component count, where it appears).
_LENGTH_TAGS = {
    'translation': 3,      # PhysicalOffsetFrame, WrapObject
    'location': 3,         # PathPoint, Marker, ContactSphere, ContactGeometry
    'mass_center': 3,      # Body
    'radius': 1,           # WrapCylinder, WrapSphere, ContactSphere
    'length': 1,           # WrapCylinder
    'dimensions': 3,       # WrapEllipsoid
    'optimal_fiber_length': 1,
    'tendon_slack_length': 1,
    'scale_factors': 3,    # Mesh, FrameGeometry -- display scale of the donor mesh
}

#: Coordinates whose value is a length rather than an angle.
_TRANSLATIONAL_COORDINATES = ('pelvis_tx', 'pelvis_ty', 'pelvis_tz')


def _scale_text(node, factor):
    values = [float(x) for x in node.text.split()]
    node.text = ' '.join(repr(v * factor) for v in values)
    return len(values)


def _scale_function_output(node, factor):
    """Multiply the OUTPUT of a scalar joint-transform function by ``factor``."""
    tag = node.tag
    if tag == 'MultiplierFunction':
        node.find('scale').text = repr(float(node.findtext('scale', '1')) * factor)
        return 1
    if tag == 'Constant':
        node.find('value').text = repr(float(node.findtext('value')) * factor)
        return 1
    if tag in ('PolynomialFunction', 'LinearFunction'):
        return _scale_text(node.find('coefficients'), factor)
    if tag == 'SimmSpline':
        return _scale_text(node.find('y'), factor)
    raise ValueError('Unsupported joint transform function: ' + tag)


def scale_model(model_bytes, scale, *, mass_scale=None, force_scale=None):
    """Return (scaled model bytes, report).

    ``scale`` is one isotropic linear factor.  Isotropy is the assumption doing
    the work: it is what makes every musculotendon path length scale by exactly
    ``scale`` at every pose, which is what lets the fitted path polynomials be
    scaled by a single multiplication.  A per-segment scale would break that and
    would need the polynomials refitted, not rescaled.

    ``mass_scale`` defaults to ``scale**3`` and ``force_scale`` to ``scale**2``
    (geometric similarity).  Body inertia takes ``mass_scale * scale**2``, so
    the mass/inertia ratio stays dimensionally right no matter what the native
    engine's own uniform ``target_mass_kg`` renormalisation then does to both.
    """
    scale = float(scale)
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError('Positive finite scale required')
    mass_scale = scale ** 3 if mass_scale is None else float(mass_scale)
    force_scale = scale ** 2 if force_scale is None else float(force_scale)

    root = ET.fromstring(model_bytes)
    counts = {}

    def bump(key, n=1):
        counts[key] = counts.get(key, 0) + n

    # 1. Plain length-valued elements.
    for tag, components in _LENGTH_TAGS.items():
        for node in root.iter(tag):
            if node.text is None or not node.text.split():
                continue
            n = _scale_text(node, scale)
            if n != components:
                raise ValueError('<%s> has %d components, expected %d' % (tag, n, components))
            bump(tag)

    # 2. Mass and inertia.  Inertia is (Ixx Iyy Izz Ixy Ixz Iyz) about the COM.
    for body in root.iter('Body'):
        body.find('mass').text = repr(float(body.findtext('mass')) * mass_scale)
        _scale_text(body.find('inertia'), mass_scale * scale ** 2)
        bump('Body')

    # 3. Joint transform axes whose output is a translation, not a rotation.
    #    A CustomJoint's knee translations are lengths and must move with the
    #    skeleton; its rotations are angles and must not.
    for joint in root.iter('CustomJoint'):
        transform = joint.find('SpatialTransform')
        if transform is None:
            continue
        for axis in transform:
            if not (axis.get('name') or '').startswith('translation'):
                continue
            for child in axis:
                if child.tag in ('coordinates', 'axis'):
                    continue
                _scale_function_output(child, scale)
                bump('translation_axis')

    # 4. Translational coordinates: default value and clamp range are lengths.
    for coordinate in root.iter('Coordinate'):
        if coordinate.get('name') not in _TRANSLATIONAL_COORDINATES:
            continue
        for tag in ('default_value', 'range'):
            node = coordinate.find(tag)
            if node is not None and node.text and node.text.split():
                _scale_text(node, scale)
        bump('translational_coordinate')

    # 5. Muscle strength.  Separated from everything above because it is a
    #    modelling choice rather than a geometric consequence.
    for node in root.iter('max_isometric_force'):
        node.text = repr(float(node.text) * force_scale)
        bump('max_isometric_force')

    # 6. Reserve actuators. All 13 in this model act on ROTATIONAL coordinates,
    #    so their optimal_force is a moment and takes force x length. One acting
    #    on a translational coordinate would be a force and must not take the
    #    extra factor, so the type is checked rather than assumed.
    for actuator in root.iter('CoordinateActuator'):
        coordinate = (actuator.findtext('coordinate') or '').split('/')[-1]
        if coordinate in _TRANSLATIONAL_COORDINATES:
            raise ValueError(
                'CoordinateActuator %s drives the translational coordinate %s; '
                'its optimal_force is a force, not a moment, and the moment '
                'scaling below would be wrong for it'
                % (actuator.get('name'), coordinate))
        node = actuator.find('optimal_force')
        node.text = repr(float(node.text) * force_scale * scale)
        bump('optimal_force')

    report = dict(scale=scale, mass_scale=mass_scale, force_scale=force_scale,
                  elements=counts,
                  unscaled=UNSCALED_BY_DESIGN)
    return ET.tostring(root, encoding='utf-8', xml_declaration=True), report


#: Things this scaling deliberately leaves alone, and why.  Printed into every
#: variant's registration so no caller has to go looking.
UNSCALED_BY_DESIGN = (
    'Coordinate ranges of ROTATIONAL coordinates: angles are dimensionless '
    'under an isotropic scale.',
    'ExpressionBasedCoordinateForceSet (passive shoulder/elbow/hip stops): the '
    'expressions give a moment in N m as a function of angle. Under geometric '
    'similarity these should go as scale**3, but they are authored range-of-'
    'motion stops rather than measured tissue properties, and no evidence in '
    'this repository fixes their scaling. They are left at the source values, '
    'so a scaled body has the SOURCE subject\'s passive joint stops.',
    'ContactForceSet stiffness/dissipation: Hunt-Crossley parameters in '
    'N/m^(3/2); their scaling is a contact-model question, not a geometry one, '
    'and they are left alone.',
    'Muscle curve shape parameters (ActiveForceLengthCurve, ForceVelocityCurve, '
    'FiberForceLengthCurve, TendonForceLengthCurve): all are in NORMALISED '
    'fibre-length and force units and are correctly invariant.',
    'pennation_angle_at_optimal: an angle.',
    'max_contraction_velocity: expressed in optimal fibre lengths per second, '
    'so it is already relative and correctly invariant.',
)


def scale_function_based_path_set(pathset_bytes, scale):
    """Scale the fitted musculotendon path polynomials by ``scale``.

    This is exact, and the reason it is exact is worth stating.  Each path's
    length is a polynomial in JOINT COORDINATES, and all 80 of them in this
    model take only rotational coordinates as arguments (verified below).  Under
    an isotropic geometric scale with the pose held fixed, every point on every
    path moves to ``scale`` times its position, so the path length at every pose
    is ``scale`` times what it was: ``L'(q) = scale * L(q)`` identically in
    ``q``.  Multiplying every coefficient by ``scale`` realises exactly that.

    The moment arms follow for free.  OpenSim's ``FunctionBasedPath`` takes them
    as ``-dL/dq`` when no explicit ``moment_arm_functions`` are given -- none are
    here -- so they scale by ``scale`` too, which is right: a moment arm is a
    length.

    If any path took a TRANSLATIONAL coordinate as an argument this would be
    wrong, because the argument itself would carry length units and the
    polynomial would need its variable rescaled as well as its value.  The
    function raises in that case rather than producing a plausible wrong number.
    """
    scale = float(scale)
    root = ET.fromstring(pathset_bytes)
    paths = 0
    for path in root.iter('FunctionBasedPath'):
        for coordinate in path.findtext('coordinate_paths').split():
            if coordinate.split('/')[-1] in _TRANSLATIONAL_COORDINATES:
                raise ValueError(
                    'Path %s is a function of the translational coordinate %s; '
                    'a uniform coefficient scale is not valid for it'
                    % (path.get('name'), coordinate))
        function = path.find('length_function/MultivariatePolynomialFunction')
        if function is None:
            raise ValueError('Only explicit polynomial path functions are supported')
        if path.find('moment_arm_functions') is not None:
            raise ValueError('Explicit moment arm functions would need scaling too')
        _scale_text(function.find('coefficients'), scale)
        paths += 1
    if not paths:
        raise ValueError('No fitted paths found')
    return ET.tostring(root, encoding='utf-8', xml_declaration=True), paths


def scale_contact_geometry_set(bytes_, scale):
    """Scale the foot contact spheres.  Locations and radii, nothing else."""
    root = ET.fromstring(bytes_)
    spheres = 0
    for sphere in root.iter('ContactSphere'):
        _scale_text(sphere.find('location'), scale)
        sphere.find('radius').text = repr(float(sphere.findtext('radius')) * scale)
        spheres += 1
    if not spheres:
        raise ValueError('No contact spheres found')
    return ET.tostring(root, encoding='utf-8', xml_declaration=True), spheres


def scale_catalog(catalog, scale, force_scale):
    """Scale the muscle catalog that rides alongside an augmented model."""
    out = []
    for entry in catalog:
        entry = dict(entry)
        for key in ('optimal_fiber_length_m', 'tendon_slack_length_m'):
            if key in entry and entry[key] is not None:
                entry[key] = float(entry[key]) * scale
        if entry.get('max_isometric_force_n') is not None:
            entry['max_isometric_force_n'] = float(entry['max_isometric_force_n']) * force_scale
        if entry.get('path_points'):
            entry['path_points'] = [
                dict(point, location=[c * scale for c in point['location']])
                if isinstance(point, dict) and 'location' in point else point
                for point in entry['path_points']]
        out.append(entry)
    return out


def path_point_polyline_lengths(model_path, pose=None):
    """Straight-line length through each muscle's declared path points, in ground.

    This is deliberately a *different* quantity from the fitted polynomial
    length, computed from a *different* part of the XML: the ``PathPoint``
    locations in the model's own ``GeometryPath``, placed by the same forward
    kinematics as the stature measurement.  It ignores wrapping, so it
    underestimates the true path wherever a ``PathWrap`` is active, and it is
    useless as an absolute check.

    Its use is the ratio.  ``polynomial_length / polyline_length`` is a number
    per muscle that depends on both files, and under a correct isotropic scale
    it must not move at all -- numerator and denominator both scale.  If either
    file were scaled and the other were not, this ratio moves by exactly the
    scale factor, and it moves whichever of the two was missed.  That is the
    check the fitted-path gate on its own cannot make, because the fitted path
    set and the model can be scaled consistently with each other and still both
    be wrong together.
    """
    root = ET.parse(model_path).getroot()
    pose = pose or default_pose_frames(model_path)
    lengths = {}
    for force in root.iter('ForceSet'):
        for muscle in force.find('objects'):
            points = muscle.find('GeometryPath/PathPointSet/objects')
            if points is None:
                continue
            world = []
            for point in points:
                body = point.findtext('socket_parent_frame').split('/')[-1]
                origin, frame = pose[body]
                world.append(origin + frame @ _vec3(point.findtext('location')))
            if len(world) < 2:
                raise ValueError('Muscle %s has fewer than two path points' % muscle.get('name'))
            lengths[muscle.get('name')] = float(
                sum(np.linalg.norm(b - a) for a, b in zip(world, world[1:])))
    if not lengths:
        raise ValueError('No muscle geometry paths found')
    return lengths
