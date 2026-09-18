"""Pose every bound anatomical entity from the mechanical body's state -- a per-frame call.

WHY.  The engine integrates 22 OpenSim bodies.  The real body is ~4,000 atlas entities in
`data/derived/canonical/anatomy.json`.  `data/derived/anatomy-segment-binding/binding.json`
(`scripts/bind_anatomy_to_segments.py`) says which segment carries each entity, and
`scripts/render_anatomical_motion.py` used it to render and export trajectories after the fact.
Nothing a running simulation could call turned the state it had just integrated into the pose
of the anatomy.  This module is that call:

    poser = AnatomyPoser.from_workspace(root)            # once
    pose = poser.pose(bodies)                              # every frame: {body: 4x4 ground}
    pose = poser.pose_from_native(native_frame)            # a NativeMechanicalStream frame
    pose = poser.pose_from_coordinates(q)                  # a motion file's coordinates
    x = pose.transform_vertices(entity_id, rest_vertices)  # any entity's surface, posed
    skin = poser.skin_vertices(pose)                       # the skin, linear-blended

`pose.rotation[i]`, `pose.translation[i]` carry entity `pose.entity_ids[i]` from its atlas rest
position to its current one: `x_now = R x_rest + t`, metres, frame `bodyparts3d-display-m`.

GEOMETRY.  Let A be the similarity carrying OpenSim ground into the atlas frame (binding.json,
scale 0.963) and T_ref[s] segment s at the registered reference pose.  In the default
`pivot='opensim'` mode an entity on s moves by

    M[s] = A . T_cur[s] . T_ref[s]^-1 . A^-1

which is the identity at the reference pose, rigid (the scale cancels), and follows the
simulated segment exactly.  The joint it rotates about is OpenSim's joint centre carried by A,
which misses the atlas's own joint by 38 mm median (binding report), and a rigid binding opens a
gap there as the joint moves.  `pivot='anatomical'` keeps every segment's ORIENTATION exactly
and re-seats each joint's rotation at the point where the two segments' bone surfaces are
closest -- IBM-1's site exporter's fix for the knee it had put 5 cm below the kneecap
(`docs/EMBODIMENT.md` s.9, `docs/LOG.md` 17 Sep) -- at the price of no longer following the
simulated segment's translation exactly.  Which opens less is measured, not assumed:
`scripts/verify_anatomy_pose.py`.

FRAMES, asserted at the boundary, never inferred from a bounding box (the extended atlas's
display export is a normalised +-1 box, 15.7% from these metres, and "the boxes agree" hid it).
  * load: anatomy.json and binding.json both declare `bodyparts3d-display-m` in metres; the
    binding was fitted on THIS anatomy.json and THIS model (sha256); A's scale is the declared
    0.963 registration (a binding fitted to the normalised box would read ~1.11); the skin mesh
    measures the anatomical stature 1.7195 m (a +-1 box reads ~2.0).
  * every call: exactly the 22 bodies; each a proper rigid 4x4; and every rotation-only joint's
    parent and child offset frames coincide to 0.1 mm -- an invariant of the model that fails
    for millimetres, for a scaled frame and for the atlas frame (see `joint_residuals_m`).
  * output: `transform_vertices` refuses any frame but `bodyparts3d-display-m`.

WHAT IS NOT POSED, AND SAYS SO.  `pose.unbound` lists every anatomy.json entity with no rigid
pose and why, on every frame; nothing is silently left at rest.  The skin is posed by the
existing graph-regularised linear blend (`continuous_surface_binding.json.gz`) through
`skin_vertices`, and its three layers have no geometry of their own and follow it.

LEFT/RIGHT.  binding.json assigns 11 of 1,342 mirrored pairs to unmirrored segments (left
tibialis anterior on the tibia, right on the calcaneus).  An asymmetric binding tore the site
figure's skin at the waist.  With `symmetric=True` (default) each such pair takes the
assignment of its higher-coherence side, mirrored; ties go to the smaller mean surface distance.
Every override is listed in `poser.symmetry_overrides`.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

FRAME = 'bodyparts3d-display-m'
OPENSIM_FRAME = 'opensim-ground-m'
REFUSED_FRAMES = ('z-anatomy-display-normalized',)
BINDING = 'data/derived/anatomy-segment-binding/binding.json'
ANATOMY = 'data/derived/canonical/anatomy.json'
MODEL = 'data/models/engineering_stance_v1/model.osim'
SKIN_BINDING = 'data/derived/canonical/continuous_surface_binding.json.gz'
SKIN_ID = 'body-bp3d-FJ2810'
SKIN_LAYERS = ('body-skin-epidermis', 'body-skin-dermis', 'body-skin-hypodermis')
JOINT_TOLERANCE_M = 1e-4
REGISTRATION_SCALE = 0.963          # ihm.body_parameters.ANATOMY_REGISTRATION_SCALE
REGISTRATION_SCALE_TOLERANCE = 0.005
PIVOT_PAIR_FRACTION = 0.04


class FrameError(ValueError):
    """An input or output is not in the frame this module's contract names."""


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------------------------
# forward kinematics of model.osim, pure python.  Checked against Simbody's own transforms in
# native frames by scripts/verify_anatomy_pose.py.
# ---------------------------------------------------------------------------------------------
def _floats(node, default=(0.0, 0.0, 0.0)):
    if node is None or not (node.text or '').strip():
        return np.array(default, float)
    return np.array([float(v) for v in node.text.split()], float)


def _rot_xyz(a, b, c):
    ca, sa, cb, sb, cc, sc = math.cos(a), math.sin(a), math.cos(b), math.sin(b), math.cos(c), math.sin(c)
    rx = np.array([[1, 0, 0], [0, ca, -sa], [0, sa, ca]], float)
    ry = np.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]], float)
    rz = np.array([[cc, -sc, 0], [sc, cc, 0], [0, 0, 1]], float)
    return rx @ ry @ rz


def _axis_rot(axis, angle):
    axis = np.asarray(axis, float)
    n = np.linalg.norm(axis)
    if n < 1e-12:
        return np.eye(3)
    k = axis / n
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]], float)
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * (K @ K)


def _xform(R, t):
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = t
    return T


class SimmSpline:
    """OpenSim's SimmSpline: the Forsythe-Malcolm-Moler cubic, end slopes from the cubic through
    the first and last four knots, linear extrapolation outside.  A natural spline is not the
    same curve and moved the patella by up to 6.5 mm against Simbody."""

    def __init__(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        n = len(x)
        self.x, self.y = x, y
        b, c, d = np.zeros(n), np.zeros(n), np.zeros(n)
        if n == 2:
            b[:] = (y[1] - y[0]) / (x[1] - x[0])
        elif n >= 3:
            nm1 = n - 1
            d[0] = x[1] - x[0]
            c[1] = (y[1] - y[0]) / d[0]
            for i in range(1, nm1):
                d[i] = x[i + 1] - x[i]
                b[i] = 2.0 * (d[i - 1] + d[i])
                c[i + 1] = (y[i + 1] - y[i]) / d[i]
                c[i] = c[i + 1] - c[i]
            b[0] = -d[0]
            b[nm1] = -d[n - 2]
            c[0] = c[nm1] = 0.0
            if n > 3:
                d1 = c[2] / (x[3] - x[1]) - c[1] / (x[2] - x[0])
                d2 = c[nm1 - 1] / (x[nm1] - x[n - 3]) - c[n - 3] / (x[nm1 - 1] - x[n - 4])
                c[0] = d1 * d[0] * d[0] / (x[3] - x[0])
                c[nm1] = -d2 * d[n - 2] * d[n - 2] / (x[nm1] - x[n - 4])
            for i in range(1, n):
                t = d[i - 1] / b[i - 1]
                b[i] -= t * d[i - 1]
                c[i] -= t * c[i - 1]
            c[nm1] /= b[nm1]
            for i in range(nm1 - 1, -1, -1):
                c[i] = (c[i] - d[i] * c[i + 1]) / b[i]
            b[nm1] = (y[nm1] - y[n - 2]) / d[n - 2] + d[n - 2] * (c[n - 2] + 2.0 * c[nm1])
            for i in range(nm1):
                b[i] = (y[i + 1] - y[i]) / d[i] - d[i] * (c[i + 1] + 2.0 * c[i])
                d[i] = (c[i + 1] - c[i]) / d[i]
                c[i] *= 3.0
            c[nm1] *= 3.0
            d[nm1] = d[n - 2]
        self.b, self.c, self.d = b, c, d

    def __call__(self, q):
        x, y, n = self.x, self.y, len(self.x)
        if n == 1:
            return float(y[0])
        if q < x[0]:
            return float(y[0] + (q - x[0]) * self.b[0])
        if q > x[-1]:
            return float(y[-1] + (q - x[-1]) * self.b[-1])
        k = min(int(np.searchsorted(x, q, side='right')) - 1, n - 2)
        dx = q - x[k]
        return float(y[k] + dx * (self.b[k] + dx * (self.c[k] + dx * self.d[k])))


def _function(node):
    if node is None:
        return lambda q: 0.0
    tag, kids = node.tag, list(node)
    if tag in ('function', 'coupled_coordinates_function') or tag.endswith('_function'):
        return _function(kids[0]) if kids else (lambda q: 0.0)
    if tag == 'Constant':
        v = float(node.findtext('value', '0'))
        return lambda q, v=v: v
    if tag == 'LinearFunction':
        c = _floats(node.find('coefficients'), (1.0, 0.0))
        return lambda q, c=c: c[0] * q + c[1]
    if tag == 'PolynomialFunction':
        c = _floats(node.find('coefficients'), (0.0,))
        return lambda q, c=c: float(np.polyval(c, q))
    if tag == 'MultiplierFunction':
        inner = node.find('function')
        sub = _function(list(inner)[0]) if inner is not None and len(inner) else (lambda q: 0.0)
        s = float(node.findtext('scale', '1'))
        return lambda q, sub=sub, s=s: sub(q) * s
    if tag == 'SimmSpline':
        return SimmSpline(_floats(node.find('x'), (0.0,)), _floats(node.find('y'), (0.0,)))
    if tag == 'PiecewiseLinearFunction':
        x, y = _floats(node.find('x'), (0.0,)), _floats(node.find('y'), (0.0,))
        return lambda q, x=x, y=y: float(np.interp(q, x, y))
    raise ValueError(f'unsupported OpenSim function {tag!r}; add it rather than guess')


class OsimKinematics:
    """Bodies, joints, coordinates and coupler constraints of an .osim, and its FK."""

    JOINT_TYPES = ('CustomJoint', 'PinJoint', 'WeldJoint')

    def __init__(self, path: Path):
        model = ET.parse(path).getroot().find('Model')
        self.bodies = [b.get('name') for b in model.iter('Body')]
        self.coordinates, self.joints = {}, []
        for jt in self.JOINT_TYPES:
            for j in model.iter(jt):
                self.joints.append(self._joint(j, jt))
        for other in ('BallJoint', 'FreeJoint', 'SliderJoint', 'UniversalJoint', 'PlanarJoint',
                      'GimbalJoint', 'EllipsoidJoint'):
            if next(model.iter(other), None) is not None:
                raise ValueError(f'{other} in {path.name}: not implemented here')
        self.couplers = []
        for c in model.iter('CoordinateCouplerConstraint'):
            self.couplers.append(dict(
                independent=(c.findtext('independent_coordinate_names') or '').split(),
                dependent=(c.findtext('dependent_coordinate_name') or '').strip(),
                f=_function(c.find('coupled_coordinates_function'))))
        placed, order, pending = {'ground'}, [], list(self.joints)
        while pending:
            rest = [j for j in pending if j['parent'] not in placed]
            for j in pending:
                if j['parent'] in placed:
                    order.append(j)
                    placed.add(j['child'])
            if len(rest) == len(pending):
                raise ValueError('disconnected joints: ' + ', '.join(j['name'] for j in rest))
            pending = rest
        self.order = order
        self.parent = {j['child']: j['parent'] for j in order}
        self.joint_of = {j['child']: j for j in order}

    def _joint(self, j, jt):
        frames = {}
        fr = j.find('frames')
        if fr is not None:
            for f in fr.iter('PhysicalOffsetFrame'):
                frames[f.get('name')] = ((f.findtext('socket_parent') or '').split('/')[-1],
                                         _xform(_rot_xyz(*_floats(f.find('orientation'))),
                                                _floats(f.find('translation'))))

        def resolve(sock):
            short = (j.findtext(sock) or '').strip().split('/')[-1]
            return frames.get(short, (short, np.eye(4)))

        parent, Xpf = resolve('socket_parent_frame')
        child, Xcf = resolve('socket_child_frame')
        coords = []
        cs = j.find('coordinates')
        if cs is not None:
            for c in cs.iter('Coordinate'):
                name = c.get('name')
                coords.append(name)
                rng = _floats(c.find('range'), (-np.inf, np.inf))
                self.coordinates[name] = dict(default=float(c.findtext('default_value', '0')),
                                              range=(float(rng[0]), float(rng[1])), joint=j.get('name'))
        axes, translates = [], False
        st = j.find('SpatialTransform')
        if st is not None:
            for ax in st.iter('TransformAxis'):
                cname = (ax.findtext('coordinates') or '').strip()
                fnode = next((c for c in ax if c.tag in ('LinearFunction', 'Constant', 'PolynomialFunction',
                                                          'MultiplierFunction', 'SimmSpline',
                                                          'PiecewiseLinearFunction')), None)
                fn = _function(fnode) if fnode is not None else (lambda q: 0.0)
                name = ax.get('name')
                if name.startswith('translation') and cname:
                    translates = True
                axes.append((name, cname.split()[0] if cname else None,
                             _floats(ax.find('axis'), (1.0, 0.0, 0.0)), fn))
        return dict(name=j.get('name'), type=jt, parent=parent, child=child, Xpf=Xpf, Xcf=Xcf,
                    Xcf_inv=np.linalg.inv(Xcf), coords=coords, axes=axes, translates=translates)

    def complete(self, q: dict, fill: dict | None = None) -> dict:
        """every coordinate a value: from q, else from `fill`, else an error; couplers applied."""
        out = dict(q)
        missing = [c for c in self.coordinates if c not in out]
        coupled = {c['dependent'] for c in self.couplers}
        for c in missing:
            if c in coupled:
                continue
            if fill is None or c not in fill:
                raise KeyError(f'coordinate {c!r} has no value and no fill')
            out[c] = float(fill[c])
        for c in self.couplers:
            out[c['dependent']] = c['f'](out[c['independent'][0]])
        unknown = set(out) - set(self.coordinates)
        if unknown:
            raise KeyError(f'coordinates not in the model: {sorted(unknown)}')
        return out

    def _motion(self, j, q):
        if j['type'] == 'WeldJoint':
            return np.eye(4)
        if j['type'] == 'PinJoint':
            return _xform(_axis_rot((0, 0, 1), q[j['coords'][0]]), np.zeros(3))
        R, t = np.eye(3), np.zeros(3)
        for name, cname, axis, fn in j['axes']:
            val = fn(q[cname] if cname else 0.0)
            if name.startswith('rotation'):
                R = R @ _axis_rot(axis, val)
            else:
                t = t + axis * val
        return _xform(R, t)

    def forward(self, q: dict) -> dict:
        """complete coordinates -> {body: 4x4 body-to-ground}, OpenSim ground, metres."""
        T = {'ground': np.eye(4)}
        for j in self.order:
            T[j['child']] = T[j['parent']] @ j['Xpf'] @ self._motion(j, q) @ j['Xcf_inv']
        del T['ground']
        return T

    def joint_residuals_m(self, T: dict) -> dict:
        """for every joint with no translational DOF, the distance between its parent-side and
        child-side frame origins.  Zero in any consistent state of this model, in metres; a
        state in millimetres, in the atlas frame or in a rescaled frame breaks it."""
        out = {}
        for j in self.order:
            if j['translates'] or j['parent'] == 'ground':
                continue
            a = (T[j['parent']] @ j['Xpf'])[:3, 3]
            b = (T[j['child']] @ j['Xcf'])[:3, 3]
            out[j['name']] = float(np.linalg.norm(a - b))
        return out

    def descendants(self, body: str) -> set:
        out, frontier = {body}, [body]
        while frontier:
            b = frontier.pop()
            for j in self.order:
                if j['parent'] == b and j['child'] not in out:
                    out.add(j['child'])
                    frontier.append(j['child'])
        return out


# ---------------------------------------------------------------------------------------------
# the pose
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class AnatomyPose:
    frame: str
    pivot: str
    entity_ids: tuple
    segment_index: np.ndarray            # (N,) index into `segments`
    segments: tuple
    rotation: np.ndarray                 # (N,3,3)
    translation: np.ndarray              # (N,3): x_now = R x_rest + t
    centroid_m: np.ndarray               # (N,3) posed centroids
    segment_motion: np.ndarray           # (22,4,4) atlas-frame rigid motion per segment
    unbound: dict = field(default_factory=dict)
    index: dict = field(default_factory=dict)

    def transform_vertices(self, entity_id: str, vertices, frame: str = FRAME) -> np.ndarray:
        """rest-pose vertices of `entity_id` (in `frame`) -> posed, same frame."""
        if frame != FRAME:
            raise FrameError(f'vertices declared in {frame!r}; this pose is in {FRAME!r}. '
                             + ('The extended atlas display export is a normalised +-1 box, 15.7% '
                                'from these metres; carry it into metres first.'
                                if frame in REFUSED_FRAMES else ''))
        if entity_id in self.unbound:
            raise KeyError(f'{entity_id} is not rigidly posed: {self.unbound[entity_id]}')
        i = self.index[entity_id]
        v = np.asarray(vertices, float)
        return v @ self.rotation[i].T + self.translation[i]


class AnatomyPoser:
    def __init__(self, root: Path, binding: dict, anatomy: dict, kin: OsimKinematics, *,
                 pivot: str = 'opensim', symmetric: bool = True, check_sources: bool = True):
        if pivot not in ('opensim', 'anatomical'):
            raise ValueError(f'pivot must be opensim or anatomical, not {pivot!r}')
        self.root, self.kin, self.pivot = Path(root), kin, pivot
        # ---- frames, at the boundary
        for label, fr in (('anatomy.json', anatomy['frame']['id']), ('binding.json', binding['frame'])):
            if fr != FRAME:
                raise FrameError(f'{label} declares {fr!r}, not {FRAME!r}')
        if anatomy['frame'].get('units') != 'm':
            raise FrameError(f'anatomy.json units {anatomy["frame"].get("units")!r}, not metres')
        A = np.asarray(binding['similarity_atlas_from_opensim_ground'], float)
        scale = float(np.cbrt(np.linalg.det(A[:3, :3])))
        if abs(scale - REGISTRATION_SCALE) > REGISTRATION_SCALE_TOLERANCE:
            raise FrameError(f'binding similarity scale {scale:.4f} is not the declared '
                             f'{REGISTRATION_SCALE} OpenSim->atlas registration; a binding fitted '
                             'to the normalised display box reads ~1.11')
        self.A, self.Ainv, self.scale = A, np.linalg.inv(A), scale
        self.segments = tuple(binding['segments'])
        if set(self.segments) != set(kin.bodies):
            raise ValueError('binding segments are not the model bodies')
        self.seg_index = {s: i for i, s in enumerate(self.segments)}
        ref_q = kin.complete(binding['reference_pose_rad'])
        self.reference_pose = ref_q
        T_ref = kin.forward(ref_q)
        self.T_ref_inv = np.stack([np.linalg.inv(T_ref[s]) for s in self.segments])

        # ---- which entity rides which segment
        ents = {e['id']: e for e in anatomy['entities']}
        rows = binding['entities']
        self.unbound = {}
        for eid in ents:
            if eid not in rows:
                self.unbound[eid] = 'absent from binding.json (binding older than anatomy.json)'
        seg = {k: r['segment'] for k, r in rows.items() if r['segment'] is not None and k in ents}
        self.symmetry_overrides = self._symmetrise(rows, seg) if symmetric else []
        for k, r in rows.items():
            if r['segment'] is None:
                self.unbound[k] = ('skin: posed by the continuous linear blend, AnatomyPoser.skin_vertices'
                                   if k == SKIN_ID else
                                   'skin layer: no geometry of its own; follows the skin blend'
                                   if k in SKIN_LAYERS else
                                   f'excluded by the binding ({r.get("basis")}): a whole-body '
                                   f'{r.get("role", "structure")} no single rigid segment carries')
        self.entity_ids = tuple(sorted(seg))
        self.index = {k: i for i, k in enumerate(self.entity_ids)}
        self.segment_of = np.array([self.seg_index[seg[k]] for k in self.entity_ids])
        self.rest_centroid = np.array([ents[k]['centroid_m'] for k in self.entity_ids], float)
        bc = np.array([binding['centroids_m'][k] for k in self.entity_ids], float)
        if np.abs(bc - self.rest_centroid).max() > 1e-12:
            raise FrameError('binding centroids are not this anatomy.json\'s centroids')
        self.coverage = dict(anatomy_entities=len(ents), rigid=len(self.entity_ids),
                             skin_blend=int(SKIN_ID in ents), skin_layers_follow=sum(k in ents for k in SKIN_LAYERS),
                             not_posed=sorted(k for k in self.unbound if k != SKIN_ID and k not in SKIN_LAYERS))
        self._binding, self._anatomy_entities = binding, ents
        self._skin = None
        self._pivots = None
        if check_sources:
            self._check_sources(binding)

    # ---- construction ---------------------------------------------------------------------
    @classmethod
    def from_workspace(cls, root, **kw) -> 'AnatomyPoser':
        root = Path(root)
        binding = json.loads((root / BINDING).read_text())
        anatomy = json.loads((root / ANATOMY).read_text())
        return cls(root, binding, anatomy, OsimKinematics(root / MODEL), **kw)

    def _check_sources(self, binding):
        prov = binding['provenance']
        for key, rel in (('anatomy', ANATOMY), ('model', MODEL)):
            got = _sha256(self.root / rel)
            if prov[key]['sha256'] != got:
                raise FrameError(f'binding.json was fitted on a different {rel} '
                                 f'({prov[key]["sha256"][:12]} vs {got[:12]}); rebuild it')
        skin = self._anatomy_entities.get(SKIN_ID)
        if skin is not None:
            with gzip.open(self.root / skin['reference_geometry']['path']) as f:
                y = np.asarray(json.load(f)['positions'], float)[1::3]
            from ihm.body_parameters import ANATOMICAL_STATURE_M
            if abs((y.max() - y.min()) - ANATOMICAL_STATURE_M) > 1e-6:
                raise FrameError(f'skin mesh is {y.max() - y.min():.4f} tall, not the anatomical '
                                 f'{ANATOMICAL_STATURE_M} m: this geometry is not in metres')

    def _symmetrise(self, rows, seg):
        def mirror_name(n):
            return re.sub(r'\bleft\b', 'right', n) if re.search(r'\bleft\b', n) else None

        def mirror_seg(s):
            return s[:-2] + ('_l' if s.endswith('_r') else '_r') if s.endswith(('_r', '_l')) else s
        by_name: dict = {}
        for k, r in rows.items():
            by_name.setdefault(r['name'], []).append(k)
        out = []
        for k, r in sorted(rows.items()):
            m = mirror_name(r['name'])
            if m is None or len(by_name.get(m, [])) != 1 or len(by_name[r['name']]) != 1:
                continue
            j = by_name[m][0]
            if k not in seg or j not in seg or mirror_seg(seg[k]) == seg[j]:
                continue
            rk, rj = rows[k], rows[j]
            keep_left = (rk['coherence'], -rk['mean_distance_m']) >= (rj['coherence'], -rj['mean_distance_m'])
            if keep_left:
                old, seg[j] = seg[j], mirror_seg(seg[k])
                out.append(dict(changed=j, name=rj['name'], from_segment=old, to_segment=seg[j],
                                mirror_of=k, coherence_kept=rk['coherence'], coherence_dropped=rj['coherence']))
            else:
                old, seg[k] = seg[k], mirror_seg(seg[j])
                out.append(dict(changed=k, name=rk['name'], from_segment=old, to_segment=seg[k],
                                mirror_of=j, coherence_kept=rj['coherence'], coherence_dropped=rk['coherence']))
        return out

    # ---- the per-frame call ---------------------------------------------------------------
    def check_bodies(self, bodies: dict) -> dict:
        """assert the input is 22 proper rigid transforms in OpenSim ground metres."""
        if set(bodies) != set(self.segments):
            raise FrameError(f'expected the {len(self.segments)} model bodies, got '
                             f'missing {sorted(set(self.segments) - set(bodies))} '
                             f'extra {sorted(set(bodies) - set(self.segments))}')
        T = {}
        for b, m in bodies.items():
            m = np.asarray(m, float)
            if m.shape != (4, 4) or not np.isfinite(m).all() or not np.allclose(m[3], [0, 0, 0, 1]):
                raise FrameError(f'{b}: not a homogeneous 4x4')
            R = m[:3, :3]
            if not np.allclose(R.T @ R, np.eye(3), atol=1e-6) or np.linalg.det(R) <= 0:
                raise FrameError(f'{b}: rotation is not proper orthonormal (a scaled frame?)')
            T[b] = m
        res = self.kin.joint_residuals_m(T)
        worst = max(res, key=res.get)
        if res[worst] > JOINT_TOLERANCE_M:
            raise FrameError(f'joint {worst}: parent and child frames {1000 * res[worst]:.2f} mm apart; '
                             f'these transforms are not {OPENSIM_FRAME} states of this model '
                             '(millimetres? the atlas frame? a rescaled body?)')
        return T

    def segment_motions(self, bodies: dict) -> np.ndarray:
        T = self.check_bodies(bodies)
        cur = np.stack([T[s] for s in self.segments])
        M = self.A[None] @ cur @ self.T_ref_inv @ self.Ainv[None]
        if self.pivot == 'anatomical':
            M = self._reseat(M)
        return M

    def pose(self, bodies: dict) -> AnatomyPose:
        M = self.segment_motions(bodies)
        R = M[self.segment_of, :3, :3]
        t = M[self.segment_of, :3, 3]
        c = np.einsum('nij,nj->ni', R, self.rest_centroid) + t
        return AnatomyPose(FRAME, self.pivot, self.entity_ids, self.segment_of, self.segments,
                           R, t, c, M, dict(self.unbound), self.index)

    def pose_from_native(self, frame: dict) -> AnatomyPose:
        return self.pose({b: r['transform_ground'] for b, r in frame['bodies'].items()})

    def pose_from_coordinates(self, q: dict, fill: str | dict | None = 'reference') -> AnatomyPose:
        """a motion file's coordinates; absent ones take the registered reference pose by default
        (held at zero the pelvis would sit at the ground origin), or `fill='default'`, or a dict,
        or None to require all 33."""
        fill_map = (self.reference_pose if fill == 'reference' else
                    {k: v['default'] for k, v in self.kin.coordinates.items()} if fill == 'default' else fill)
        return self.pose(self.kin.forward(self.kin.complete(q, fill_map)))

    # ---- the skin ---------------------------------------------------------------------------
    def skin_vertices(self, pose: AnatomyPose) -> np.ndarray:
        """the FJ2810 skin, every vertex, by its graph-regularised linear blend over the same
        segment motions (`continuous_surface_binding.json.gz`)."""
        if self._skin is None:
            p = json.loads(gzip.decompress((self.root / SKIN_BINDING).read_bytes()))
            order = [self.seg_index[s['id']] for s in p['segments']]
            self._skin = (np.asarray(p['reference_positions_m'], float), np.asarray(p['weights'], float),
                          np.array(order))
        rest, w, order = self._skin
        M = pose.segment_motion[order]
        out = np.zeros_like(rest)
        for k in range(len(order)):
            nz = w[:, k] > 0
            out[nz] += w[nz, k, None] * (rest[nz] @ M[k, :3, :3].T + M[k, :3, 3])
        return out

    # ---- anatomical pivots ------------------------------------------------------------------
    def pivots(self) -> dict:
        """joint -> point where the parent's and child's bone surfaces are closest, atlas rest."""
        if self._pivots is None:
            from scipy.spatial import cKDTree
            bones = {}
            for s in self.segments:
                pts = []
                for bid in self._binding['segment_named_bones'][s]:
                    e = self._anatomy_entities[bid]
                    with gzip.open(self.root / e['reference_geometry']['path']) as f:
                        pts.append(np.asarray(json.load(f)['positions'], float).reshape(-1, 3))
                P = np.concatenate(pts)
                bones[s] = P[:: max(1, len(P) // 40000)]
            piv = {}
            for j in self.kin.order:
                if j['parent'] == 'ground':
                    continue
                a, b = bones[j['parent']], bones[j['child']]
                d, idx = cKDTree(b).query(a)
                k = max(8, int(PIVOT_PAIR_FRACTION * len(a)))
                sel = np.argsort(d)[:k]
                piv[j['name']] = dict(point_m=(0.5 * (a[sel] + b[idx[sel]])).mean(0),
                                      pair_distance_m=float(np.median(d[sel])),
                                      parent=j['parent'], child=j['child'])
            self._pivots = piv
        return self._pivots

    def _reseat(self, M):
        piv = self.pivots()
        out = M.copy()
        for j in self.kin.order:
            if j['parent'] == 'ground':
                continue
            p, c = self.seg_index[j['parent']], self.seg_index[j['child']]
            rel = np.linalg.inv(M[p]) @ M[c]
            q = piv[j['name']]['point_m']
            R = rel[:3, :3]
            seat = np.eye(4)
            seat[:3, :3] = R
            seat[:3, 3] = q - R @ q
            out[c] = out[p] @ seat
        return out
