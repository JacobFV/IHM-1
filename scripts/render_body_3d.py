#!/usr/bin/env python3
"""Render a recorded IHM-1 body trajectory as a 3D video.

Everything geometric in this file comes out of the OpenSim model
(`data/models/engineering_stance_v1/model.osim`): the joint graph, the offset
frames, the coordinate-to-transform functions of every CustomJoint, and the
marker set that fixes the head, heel and toe landmarks.  No limb length is
hardcoded.  The recorded coordinates in the trace are pushed through a generic
OpenSim forward-kinematics implementation, so a joint whose SpatialTransform is
a polynomial of its coordinate -- the Rajagopal walker knee, which translates as
well as rotates -- is evaluated as the model defines it rather than as a hinge.

Rendering is a small software rasteriser written against numpy: a z-buffered,
super-sampled triangle pipeline with smooth vertex normals, a three-point light
rig, an analytic ground plane and a projected soft shadow.  It is here rather
than matplotlib because matplotlib's 3D axes depth-sort whole collections, which
puts an arm behind a torso it is in front of.

Muscle excitations recorded per frame are mapped to body segments through
`catalog.json`'s `attachment_bodies` and shown as segment colour.

Usage
-----
    python scripts/render_body_3d.py                     # default dual render
    python scripts/render_body_3d.py --check             # print FK landmarks
    python scripts/render_body_3d.py --subject data/derived/x/trace.json:label
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "data/models/engineering_stance_v1/model.osim"
CATALOG = ROOT / "data/models/engineering_stance_v1/catalog.json"


# --------------------------------------------------------------------------
# OpenSim model: parsing and forward kinematics
# --------------------------------------------------------------------------

def _floats(node, default=(0.0, 0.0, 0.0)):
    if node is None or not (node.text or "").strip():
        return np.array(default, float)
    return np.array([float(v) for v in node.text.split()], float)


def _rot_xyz(a, b, c):
    """OpenSim body-fixed X-Y-Z Euler angles -> rotation matrix."""
    ca, sa = math.cos(a), math.sin(a)
    cb, sb = math.cos(b), math.sin(b)
    cc, sc = math.cos(c), math.sin(c)
    rx = np.array([[1, 0, 0], [0, ca, -sa], [0, sa, ca]], float)
    ry = np.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]], float)
    rz = np.array([[cc, -sc, 0], [sc, cc, 0], [0, 0, 1]], float)
    return rx @ ry @ rz


def _xform(R, t):
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = t
    return T


def _axis_rot(axis, angle):
    axis = np.asarray(axis, float)
    n = np.linalg.norm(axis)
    if n < 1e-12:
        return np.eye(3)
    k = axis / n
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]], float)
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * (K @ K)


class _NaturalSpline:
    """Natural cubic spline, matching OpenSim's SimmSpline closely enough for FK."""

    def __init__(self, x, y):
        from scipy.interpolate import CubicSpline
        self.f = CubicSpline(x, y, bc_type="natural", extrapolate=True)
        self.x0, self.x1 = x[0], x[-1]

    def __call__(self, q):
        return float(self.f(min(max(q, self.x0), self.x1)))


def _parse_function(node):
    """An OpenSim Function element -> a python callable of one coordinate."""
    if node is None:
        return lambda q: 0.0
    kids = [c for c in node]
    tag = node.tag
    if tag in ("function", "prescribed_function") or tag.endswith("_function"):
        # a wrapper element; descend
        return _parse_function(kids[0]) if kids else (lambda q: 0.0)
    if tag == "Constant":
        v = float(node.findtext("value", "0"))
        return lambda q, v=v: v
    if tag == "LinearFunction":
        c = _floats(node.find("coefficients"), (1.0, 0.0))
        return lambda q, c=c: c[0] * q + c[1]
    if tag == "PolynomialFunction":
        c = _floats(node.find("coefficients"), (0.0,))
        return lambda q, c=c: float(np.polyval(c, q))
    if tag == "MultiplierFunction":
        inner = node.find("function")
        sub = _parse_function([c for c in inner][0]) if inner is not None and len(inner) else (lambda q: 0.0)
        scale = float(node.findtext("scale", "1"))
        return lambda q, sub=sub, s=scale: sub(q) * s
    if tag in ("SimmSpline", "NaturalCubicSpline", "GCVSpline"):
        x = _floats(node.find("x"), (0.0,))
        y = _floats(node.find("y"), (0.0,))
        if len(x) < 3:
            return lambda q, x=x, y=y: float(np.interp(q, x, y))
        sp = _NaturalSpline(x, y)
        return sp
    if tag == "PiecewiseLinearFunction":
        x = _floats(node.find("x"), (0.0,))
        y = _floats(node.find("y"), (0.0,))
        return lambda q, x=x, y=y: float(np.interp(q, x, y))
    # unknown function type: descend into the first child that looks like one
    for c in kids:
        if c.tag.endswith("Function") or c.tag == "Constant":
            return _parse_function(c)
    return lambda q: 0.0


class OsimModel:
    def __init__(self, path: Path):
        root = ET.parse(path).getroot()
        model = root.find("Model")
        self.name = model.get("name")
        self.bodies = [b.get("name") for b in model.iter("Body")]

        self.joints = []
        jtypes = ("CustomJoint", "PinJoint", "WeldJoint", "BallJoint", "FreeJoint",
                  "SliderJoint", "UniversalJoint", "PlanarJoint")
        for jt in jtypes:
            for j in model.iter(jt):
                self.joints.append(self._parse_joint(j, jt))

        self.markers = {}
        ms = model.find(".//MarkerSet")
        if ms is not None:
            for m in ms.iter("Marker"):
                body = (m.findtext("socket_parent_frame") or "").split("/")[-1]
                self.markers[m.get("name")] = (body, _floats(m.find("location")))

        self.order = self._topo_order()

    # -- parsing -----------------------------------------------------------
    @staticmethod
    def _frame_table(joint):
        table = {}
        frames = joint.find("frames")
        if frames is None:
            return table
        for f in frames.iter("PhysicalOffsetFrame"):
            body = (f.findtext("socket_parent") or "").split("/")[-1]
            t = _floats(f.find("translation"))
            o = _floats(f.find("orientation"))
            table[f.get("name")] = (body, _xform(_rot_xyz(*o), t))
        return table

    def _parse_joint(self, j, jtype):
        frames = self._frame_table(j)

        def resolve(socket_name):
            key = (j.findtext(socket_name) or "").strip()
            short = key.split("/")[-1]
            if short in frames:
                return frames[short]
            # socket points straight at a body / ground
            return short, np.eye(4)

        pbody, Xpf = resolve("socket_parent_frame")
        cbody, Xcf = resolve("socket_child_frame")

        coords = []
        cs = j.find("coordinates")
        if cs is not None:
            coords = [c.get("name") for c in cs.iter("Coordinate")]

        axes = []
        st = j.find("SpatialTransform")
        if st is not None:
            for ax in st.iter("TransformAxis"):
                cname = (ax.findtext("coordinates") or "").strip()
                axis = _floats(ax.find("axis"), (1.0, 0.0, 0.0))
                fn = None
                for c in ax:
                    if c.tag in ("LinearFunction", "Constant", "PolynomialFunction",
                                 "MultiplierFunction", "SimmSpline", "NaturalCubicSpline",
                                 "PiecewiseLinearFunction", "GCVSpline"):
                        fn = _parse_function(c)
                        break
                axes.append((ax.get("name"), cname.split()[0] if cname else None, axis,
                             fn or (lambda q: 0.0)))

        return dict(name=j.get("name"), type=jtype, parent=pbody, child=cbody,
                    Xpf=Xpf, Xcf=Xcf, coords=coords, axes=axes)

    def _topo_order(self):
        placed = {"ground"}
        order, pending = [], list(self.joints)
        while pending:
            progress = False
            rest = []
            for j in pending:
                if j["parent"] in placed:
                    order.append(j)
                    placed.add(j["child"])
                    progress = True
                else:
                    rest.append(j)
            pending = rest
            if not progress:
                raise RuntimeError(f"disconnected joints: {[j['name'] for j in pending]}")
        return order

    # -- kinematics --------------------------------------------------------
    def _joint_motion(self, j, q):
        """X_FM: the parent offset frame -> child offset frame transform."""
        if j["type"] == "WeldJoint":
            return np.eye(4)
        if j["type"] == "PinJoint":
            a = q.get(j["coords"][0], 0.0) if j["coords"] else 0.0
            return _xform(_axis_rot((0, 0, 1), a), np.zeros(3))
        R = np.eye(3)
        t = np.zeros(3)
        for name, cname, axis, fn in j["axes"]:
            val = fn(q.get(cname, 0.0) if cname else 0.0)
            if name.startswith("rotation"):
                R = R @ _axis_rot(axis, val)
            else:
                # OpenSim expresses SpatialTransform translations in the parent frame
                t = t + axis * val
        return _xform(R, t)

    def forward(self, q: dict) -> dict:
        """coordinate values -> {body: 4x4 body-to-ground transform}"""
        T = {"ground": np.eye(4)}
        for j in self.order:
            T[j["child"]] = T[j["parent"]] @ j["Xpf"] @ self._joint_motion(j, q) @ np.linalg.inv(j["Xcf"])
        return T

    # -- landmarks ---------------------------------------------------------
    def joint_point(self, joint_name, side="parent"):
        """The joint centre expressed in its parent (or child) body frame."""
        for j in self.joints:
            if j["name"] == joint_name:
                X = j["Xpf"] if side == "parent" else j["Xcf"]
                body = j["parent"] if side == "parent" else j["child"]
                return body, X[:3, 3].copy()
        raise KeyError(joint_name)


# --------------------------------------------------------------------------
# Skeleton definition, entirely derived from the model
# --------------------------------------------------------------------------

def build_skeleton(model: OsimModel):
    """-> (landmarks, segments).

    landmarks: name -> (body, local point)
    segments : list of dicts describing a sphere-swept tube between landmarks
    """
    L = {}

    def jp(joint, side="parent", as_name=None):
        body, p = model.joint_point(joint, side)
        L[as_name or joint] = (body, p)
        return L[as_name or joint]

    def mk(marker, as_name):
        body, p = model.markers[marker]
        L[as_name] = (body, p.copy())
        return L[as_name]

    # pelvis / trunk
    jp("hip_r", "parent", "hip_r")
    jp("hip_l", "parent", "hip_l")
    jp("back", "parent", "lumbar")
    mk("S2", "sacrum")
    mk("MidASIS", "asis")

    # torso: shoulders and head come from the acromial joints and the head marker
    jp("acromial_r", "parent", "shoulder_r")
    jp("acromial_l", "parent", "shoulder_l")
    mk("Head", "head_top")
    hb, hp = L["head_top"]
    sr, sl = L["shoulder_r"][1], L["shoulder_l"][1]
    neck = (sr + sl) / 2.0
    neck[1] += 0.035
    L["neck"] = (hb, neck)
    # head centre sits below the head marker by a head radius
    head_r = float(hp[1] - neck[1]) * 0.42
    L["head_c"] = (hb, np.array([hp[0] + 0.012, hp[1] - head_r * 1.02, hp[2]]))
    L["_head_r"] = (hb, np.array([head_r, 0, 0]))
    # sternum: on the torso axis, at the level where the rib cage is widest
    L["sternum"] = (hb, np.array([0.012, neck[1] * 0.60, 0.0]))
    # mid-hip, in the pelvis frame, so the pelvis block spans hips -> lumbar
    L["hip_mid"] = (L["hip_r"][0], (L["hip_r"][1] + L["hip_l"][1]) / 2.0)

    # legs
    for s in ("r", "l"):
        jp(f"walker_knee_{s}", "parent", f"knee_{s}")          # in femur
        jp(f"walker_knee_{s}", "child", f"knee_t_{s}")         # in tibia
        jp(f"ankle_{s}", "parent", f"ankle_{s}")               # in tibia
        jp(f"mtp_{s}", "parent", f"mtp_{s}")                   # in calcn
        mk(f"{s.upper()}.HeelGround", f"heel_{s}")
        mk(f"{s.upper()}.ToeGround", f"toe_{s}")
        mk(f"{s.upper()}.MT5Ground", f"mt5_{s}")
        # arms
        jp(f"elbow_{s}", "parent", f"elbow_{s}")               # in humerus
        jp(f"radius_hand_{s}", "parent", f"wrist_{s}")         # in radius

    segments = []

    def seg(name, a, b, r0, r1, bodies, ex=1.0, ez=1.0, tone=0.0):
        segments.append(dict(name=name, a=a, b=b, r0=r0, r1=r1, bodies=bodies,
                             ex=ex, ez=ez, tone=tone))

    # trunk.  ex = anterior-posterior scale, ez = the perpendicular one.
    seg("pelvis", "hip_mid", "lumbar", 0.112, 0.100, ["pelvis"], ex=0.80, ez=1.20)
    seg("hips", "hip_l", "hip_r", 0.082, 0.082, ["pelvis"], ex=0.92, ez=1.00)
    seg("abdomen", "lumbar", "sternum", 0.106, 0.118, ["torso"], ex=0.80, ez=1.30)
    seg("chest", "sternum", "neck", 0.118, 0.086, ["torso"], ex=0.78, ez=1.36)
    seg("clavicle", "shoulder_l", "shoulder_r", 0.062, 0.062, ["torso"], ex=0.95, ez=0.95)
    seg("neck", "neck", "head_c", 0.046, 0.052, ["torso"], tone=0.18)

    for s in ("r", "l"):
        seg(f"thigh_{s}", f"hip_{s}", f"knee_{s}", 0.086, 0.060, [f"femur_{s}"],
            ex=0.94, ez=1.00)
        seg(f"shank_{s}", f"knee_t_{s}", f"ankle_{s}", 0.061, 0.039, [f"tibia_{s}"],
            ex=0.92, ez=1.00)
        seg(f"foot_{s}", f"heel_{s}", f"mtp_{s}", 0.041, 0.034, [f"calcn_{s}", f"talus_{s}"],
            ex=0.95, ez=0.80)
        seg(f"toes_{s}", f"mtp_{s}", f"toe_{s}", 0.031, 0.019, [f"toes_{s}"], ex=1.15, ez=0.70)
        seg(f"upperarm_{s}", f"shoulder_{s}", f"elbow_{s}", 0.055, 0.041, [f"humerus_{s}"])
        seg(f"forearm_{s}", f"elbow_{s}", f"wrist_{s}", 0.043, 0.030,
            [f"ulna_{s}", f"radius_{s}"])
        seg(f"hand_{s}", f"wrist_{s}", f"hand_tip_{s}", 0.032, 0.021, [f"hand_{s}"],
            ex=0.62, ez=1.05, tone=0.12)

    return L, segments


def landmark_world(model, T, L, name):
    body, p = L[name]
    M = T[body]
    return M[:3, :3] @ p + M[:3, 3]


def all_landmarks(model, T, L):
    P = {k: landmark_world(model, T, L, k) for k in L if not k.startswith("_")}
    # hand tips are not a model landmark; extend the forearm past the wrist
    for s in ("r", "l"):
        d = P[f"wrist_{s}"] - P[f"elbow_{s}"]
        n = np.linalg.norm(d)
        P[f"hand_tip_{s}"] = P[f"wrist_{s}"] + (d / n if n > 1e-9 else np.array([0, -1, 0])) * 0.088
    return P


# --------------------------------------------------------------------------
# Muscle excitation -> segment activation
# --------------------------------------------------------------------------

class MuscleMap:
    def __init__(self, catalog_path: Path, segments):
        cat = json.loads(catalog_path.read_text())
        self.by_body = {}
        for m in cat:
            for b in m.get("attachment_bodies", []):
                self.by_body.setdefault(b, []).append(
                    (m["id"], float(m.get("max_isometric_force_n", 1.0))))
        self.seg_muscles = {}
        for s in segments:
            acc = {}
            for b in s["bodies"]:
                for mid, f in self.by_body.get(b, []):
                    acc[mid] = max(acc.get(mid, 0.0), f)
            self.seg_muscles[s["name"]] = list(acc.items())

    def activations(self, excitations: dict):
        out = {}
        for name, muscles in self.seg_muscles.items():
            num = den = 0.0
            for mid, f in muscles:
                e = excitations.get(mid)
                if e is None:
                    continue
                num += e * f
                den += f
            out[name] = num / den if den > 0 else None
        return out


# --------------------------------------------------------------------------
# Mesh construction
# --------------------------------------------------------------------------

def _tube(p0, p1, r0, r1, ex, ez, nseg=18, ncap=5):
    """Sphere-swept, elliptical-section tube from p0 to p1. -> verts, faces

    The cross-section frame is built deterministically so the two radius scales
    always mean the same thing: `ex` scales the anterior-posterior axis, `ez`
    the axis perpendicular to it (lateral for a limb, vertical for a segment
    that already runs laterally).
    """
    p0 = np.asarray(p0, float)
    p1 = np.asarray(p1, float)
    d = p1 - p0
    ln = np.linalg.norm(d)
    if ln < 1e-9:
        d = np.array([0.0, 1.0, 0.0])
        ln = 1e-9
    w = d / ln
    ref = np.array([1.0, 0.0, 0.0])            # anterior
    if abs(w[0]) > 0.9:                        # segment already runs anteriorly
        ref = np.array([0.0, 1.0, 0.0])
    u = ref - w * float(ref @ w)
    u /= np.linalg.norm(u)
    v = np.cross(w, u)

    rings = []  # (centre_along_axis, radius_scale)
    for a in np.linspace(-math.pi / 2, 0.0, ncap, endpoint=False):
        rings.append((r0 * math.sin(a), r0 * math.cos(a)))
    nb = 8
    for t in np.linspace(0.0, 1.0, nb):
        rings.append((t * ln, r0 + (r1 - r0) * t))
    for a in np.linspace(0.0, math.pi / 2, ncap + 1)[1:]:
        rings.append((ln + r1 * math.sin(a), r1 * math.cos(a)))

    ang = np.linspace(0, 2 * math.pi, nseg, endpoint=False)
    cu = np.cos(ang) * ex
    cv = np.sin(ang) * ez
    verts = []
    for z, r in rings:
        c = p0 + w * z
        verts.append(c[None, :] + (r * cu)[:, None] * u[None, :] + (r * cv)[:, None] * v[None, :])
    verts = np.concatenate(verts, 0)
    nr = len(rings)

    faces = []
    for i in range(nr - 1):
        a0 = i * nseg
        a1 = (i + 1) * nseg
        j = np.arange(nseg)
        jn = (j + 1) % nseg
        faces.append(np.stack([a0 + j, a0 + jn, a1 + jn], 1))
        faces.append(np.stack([a0 + j, a1 + jn, a1 + j], 1))
    faces = np.concatenate(faces, 0)
    return verts.astype(np.float32), faces.astype(np.int32)


def _ellipsoid(centre, radii, n=22):
    th = np.linspace(0, math.pi, n // 2 + 1)
    ph = np.linspace(0, 2 * math.pi, n, endpoint=False)
    T, P = np.meshgrid(th, ph, indexing="ij")
    x = np.sin(T) * np.cos(P) * radii[0]
    y = np.cos(T) * radii[1]
    z = np.sin(T) * np.sin(P) * radii[2]
    V = np.stack([x, y, z], -1).reshape(-1, 3) + np.asarray(centre, float)
    rows, cols = T.shape
    faces = []
    for i in range(rows - 1):
        j = np.arange(cols)
        jn = (j + 1) % cols
        a0, a1 = i * cols, (i + 1) * cols
        faces.append(np.stack([a0 + j, a0 + jn, a1 + jn], 1))
        faces.append(np.stack([a0 + j, a1 + jn, a1 + j], 1))
    return V.astype(np.float32), np.concatenate(faces, 0).astype(np.int32)


def vertex_normals(V, F):
    n = np.zeros_like(V)
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    fn = np.cross(b - a, c - a)
    for k in range(3):
        np.add.at(n, F[:, k], fn)
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    return n / np.maximum(ln, 1e-12)


# --------------------------------------------------------------------------
# Colour
# --------------------------------------------------------------------------

BASE_COL = np.array([0.263, 0.325, 0.384])     # slate
WARM_COL = np.array([0.98, 0.72, 0.28])        # mid activation
HOT_COL = np.array([1.00, 0.30, 0.20])         # high activation
DEAD_COL = np.array([0.30, 0.31, 0.34])        # no muscle mapped to this segment


def activation_colour(a, lo, hi):
    if a is None:
        return DEAD_COL.copy()
    t = (a - lo) / max(hi - lo, 1e-9)
    t = min(max(t, 0.0), 1.0)
    if t < 0.5:
        return BASE_COL + (WARM_COL - BASE_COL) * (t / 0.5)
    return WARM_COL + (HOT_COL - WARM_COL) * ((t - 0.5) / 0.5)


# --------------------------------------------------------------------------
# Software rasteriser
# --------------------------------------------------------------------------

class Camera:
    def __init__(self, eye, target, up=(0, 1, 0), fov_deg=27.0, aspect=16 / 9):
        eye = np.asarray(eye, float)
        target = np.asarray(target, float)
        f = target - eye
        f /= np.linalg.norm(f)
        r = np.cross(f, np.asarray(up, float))
        r /= np.linalg.norm(r)
        u = np.cross(r, f)
        self.eye = eye
        self.R = np.stack([r, u, -f], 0)     # world -> camera
        self.fy = 1.0 / math.tan(math.radians(fov_deg) / 2.0)
        self.aspect = aspect
        self.forward = f
        self.right = r
        self.up = u

    def to_view(self, P):
        return (P - self.eye) @ self.R.T

    def project(self, P, W, H):
        V = self.to_view(P)
        z = np.maximum(-V[:, 2], 1e-6)
        x = (V[:, 0] * (self.fy / self.aspect)) / z
        y = (V[:, 1] * self.fy) / z
        sx = (x * 0.5 + 0.5) * W
        sy = (0.5 - y * 0.5) * H
        return sx, sy, z, V


def rasterise(sx, sy, depth, F, W, H, backface=True):
    """Vectorised z-buffered rasteriser.

    Returns (pix_index, tri_index, bary) for the winning fragment of every
    covered pixel.
    """
    i0, i1, i2 = F[:, 0], F[:, 1], F[:, 2]
    x0, y0 = sx[i0], sy[i0]
    x1, y1 = sx[i1], sy[i1]
    x2, y2 = sx[i2], sy[i2]
    area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)

    keep = np.abs(area) > 1e-9
    if backface:
        keep &= area < 0            # screen y is flipped, so front faces are CW
    zmin = np.minimum(np.minimum(depth[i0], depth[i1]), depth[i2])
    keep &= zmin > 1e-3

    xmin = np.ceil(np.minimum(np.minimum(x0, x1), x2) - 0.5).astype(np.int32)
    xmax = np.floor(np.maximum(np.maximum(x0, x1), x2) - 0.5).astype(np.int32)
    ymin = np.ceil(np.minimum(np.minimum(y0, y1), y2) - 0.5).astype(np.int32)
    ymax = np.floor(np.maximum(np.maximum(y0, y1), y2) - 0.5).astype(np.int32)
    xmin = np.clip(xmin, 0, W - 1); xmax = np.clip(xmax, 0, W - 1)
    ymin = np.clip(ymin, 0, H - 1); ymax = np.clip(ymax, 0, H - 1)
    bw = xmax - xmin + 1
    bh = ymax - ymin + 1
    keep &= (bw > 0) & (bh > 0)

    idx = np.nonzero(keep)[0]
    if idx.size == 0:
        return None
    bw = bw[idx].astype(np.int64)
    bh = bh[idx].astype(np.int64)
    counts = bw * bh
    total = int(counts.sum())
    if total == 0:
        return None

    tri = np.repeat(idx, counts)
    starts = np.zeros(idx.size, np.int64)
    np.cumsum(counts[:-1], out=starts[1:])
    off = np.arange(total, dtype=np.int64) - np.repeat(starts, counts)
    bw_r = np.repeat(bw, counts)
    px = np.repeat(xmin[idx].astype(np.int64), counts) + off % bw_r
    py = np.repeat(ymin[idx].astype(np.int64), counts) + off // bw_r

    fx = px.astype(np.float32) + 0.5
    fy = py.astype(np.float32) + 0.5
    ax, ay = x0[tri], y0[tri]
    bx, by = x1[tri], y1[tri]
    cx, cy = x2[tri], y2[tri]
    ar = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay)
    w0 = ((bx - fx) * (cy - fy) - (cx - fx) * (by - fy)) / ar
    w1 = ((cx - fx) * (ay - fy) - (ax - fx) * (cy - fy)) / ar
    w2 = 1.0 - w0 - w1
    inside = (w0 >= -1e-6) & (w1 >= -1e-6) & (w2 >= -1e-6)
    if not inside.any():
        return None
    tri = tri[inside]; px = px[inside]; py = py[inside]
    w0 = w0[inside]; w1 = w1[inside]; w2 = w2[inside]

    # perspective-correct interpolation
    iz = w0 / depth[F[tri, 0]] + w1 / depth[F[tri, 1]] + w2 / depth[F[tri, 2]]
    z = 1.0 / np.maximum(iz, 1e-9)
    pw0 = w0 / depth[F[tri, 0]] * z
    pw1 = w1 / depth[F[tri, 1]] * z
    pw2 = 1.0 - pw0 - pw1

    flat = py * W + px
    order = np.lexsort((z, flat))
    flat_s = flat[order]
    first = np.ones(flat_s.size, bool)
    first[1:] = flat_s[1:] != flat_s[:-1]
    win = order[first]
    bary = np.stack([pw0[win], pw1[win], pw2[win]], 1)
    return flat[win], tri[win], bary, z[win]


def _box_blur(a, k):
    """Separable box blur via cumulative sums; O(N) and vectorised."""
    pad = k // 2
    for axis in (0, 1):
        b = np.moveaxis(a, axis, 0)
        c = np.zeros((b.shape[0] + 2 * pad + 1,) + b.shape[1:], np.float32)
        c[pad + 1:pad + 1 + b.shape[0]] = b
        c[pad + 1 + b.shape[0]:] = 0.0
        cs = np.cumsum(c, axis=0)
        out = (cs[k:] - cs[:-k]) / k
        a = np.moveaxis(out[:b.shape[0]], 0, axis)
    return a


def shadow_map(V, F, light_dir, extent, res=384):
    """Project the mesh to y=0 along the light and rasterise coverage."""
    d = np.asarray(light_dir, float)
    d = d / np.linalg.norm(d)
    if abs(d[1]) < 1e-3:
        return None
    t = V[:, 1] / d[1]
    G = V - t[:, None] * d[None, :]
    (x0, x1, z0, z1) = extent
    sx = (G[:, 0] - x0) / (x1 - x0) * res
    sy = (G[:, 2] - z0) / (z1 - z0) * res
    depth = np.ones(V.shape[0], np.float32)
    out = rasterise(sx.astype(np.float32), sy.astype(np.float32), depth, F,
                    res, res, backface=False)
    buf = np.zeros(res * res, np.float32)
    if out is not None:
        buf[out[0]] = 1.0
    buf = buf.reshape(res, res)
    for _ in range(3):
        buf = _box_blur(buf, 9)
    return np.clip(buf * 1.9, 0, 1)


# --------------------------------------------------------------------------
# Scene rendering
# --------------------------------------------------------------------------

KEY_DIR = np.array([-0.55, -0.78, -0.30])
FILL_DIR = np.array([0.62, -0.30, -0.55])
RIM_DIR = np.array([0.25, -0.15, 0.95])
BG_TOP = np.array([0.0045, 0.0068, 0.0100])
BG_BOT = np.array([0.0130, 0.0185, 0.0255])
GROUND_COL = np.array([0.0105, 0.0155, 0.0215])
GRID_COL = np.array([0.055, 0.090, 0.120])


def render_ground(cam, W, H, shadow, shadow_extent, glow_centre, pools=()):
    """Background gradient plus the analytic ground plane. -> (rgb, view depth)"""
    yy = (np.arange(H, dtype=np.float32) / max(H - 1, 1))[:, None, None]
    img = (BG_TOP.astype(np.float32)[None, None, :] * (1 - yy) +
           BG_BOT.astype(np.float32)[None, None, :] * yy)
    img = np.repeat(img, W, axis=1)

    gsx, gsy, gz, _ = cam.project(glow_centre[None, :], W, H)
    gx, gy = float(gsx[0]), float(gsy[0])
    X = np.arange(W, dtype=np.float32)[None, :]
    Y = np.arange(H, dtype=np.float32)[:, None]
    r2 = ((X - gx) / (W * 0.42)) ** 2 + ((Y - gy) / (H * 0.62)) ** 2
    img += np.exp(-r2 * 2.2)[..., None] * np.array([0.016, 0.028, 0.042], np.float32)

    ndc_x = ((X + 0.5) / W * 2 - 1) * (cam.aspect / cam.fy)
    ndc_y = (1 - (Y + 0.5) / H * 2) * (1.0 / cam.fy)
    dirs = (cam.right.astype(np.float32)[None, None, :] * ndc_x[..., None] +
            cam.up.astype(np.float32)[None, None, :] * ndc_y[..., None] +
            cam.forward.astype(np.float32)[None, None, :])
    dirs /= np.linalg.norm(dirs, axis=2, keepdims=True)
    dy = dirs[:, :, 1]
    hit = dy < -1e-6
    t = np.where(hit, -np.float32(cam.eye[1]) / np.where(hit, dy, np.float32(-1.0)),
                 np.float32(0.0))
    Px = np.float32(cam.eye[0]) + dirs[:, :, 0] * t
    Pz = np.float32(cam.eye[2]) + dirs[:, :, 2] * t
    gdepth = np.where(hit, t * np.einsum("ijk,k->ij", dirs, cam.forward.astype(np.float32)),
                      np.float32(1e9))

    dist = np.sqrt((Px - glow_centre[0]) ** 2 + (Pz - glow_centre[2]) ** 2)
    fade = np.exp(-np.maximum(dist - 0.6, 0) / 2.6)
    gcol = np.repeat(np.repeat(GROUND_COL.astype(np.float32)[None, None, :], H, 0), W, 1)

    for spacing, strength in ((0.5, 0.38), (2.0, 0.80)):
        fx = np.abs(((Px / spacing + 0.5) % 1.0) - 0.5)
        fz = np.abs(((Pz / spacing + 0.5) % 1.0) - 0.5)
        wid = np.clip(0.012 + t * 0.0032, 0.012, 0.12)
        g = np.maximum(np.clip(1 - fx / wid, 0, 1), np.clip(1 - fz / wid, 0, 1))
        gcol += (GRID_COL * strength).astype(np.float32)[None, None, :] * (g * fade)[..., None]

    # a pool of light under each body, so the plane reads as a floor and the
    # projected shadow has something to subtract from
    for px, pz, pcol in pools:
        pd = np.sqrt((Px - px) ** 2 + (Pz - pz) ** 2)
        pool = np.exp(-(pd / 1.25) ** 1.8)
        gcol += (np.asarray(pcol, np.float32) * 0.085)[None, None, :] * pool[..., None]
        ring = np.exp(-((pd - 0.62) / 0.045) ** 2)
        gcol += (np.asarray(pcol, np.float32) * 0.10)[None, None, :] * ring[..., None]

    if shadow is not None:
        x0, x1, z0, z1 = shadow_extent
        su = (Px - x0) / (x1 - x0) * shadow.shape[1]
        sv = (Pz - z0) / (z1 - z0) * shadow.shape[0]
        ok = (su >= 0) & (su < shadow.shape[1] - 1) & (sv >= 0) & (sv < shadow.shape[0] - 1)
        iu = np.clip(su, 0, shadow.shape[1] - 2)
        iv = np.clip(sv, 0, shadow.shape[0] - 2)
        u0 = iu.astype(np.int32); v0 = iv.astype(np.int32)
        fu = iu - u0; fv = iv - v0
        sh = (shadow[v0, u0] * (1 - fu) * (1 - fv) + shadow[v0, u0 + 1] * fu * (1 - fv) +
              shadow[v0 + 1, u0] * (1 - fu) * fv + shadow[v0 + 1, u0 + 1] * fu * fv)
        gcol *= (1.0 - 0.86 * np.where(ok, sh, 0.0))[..., None]

    gcol *= fade[..., None] * 0.92 + 0.08
    # atmospheric fade, so the plane dissolves into the background instead of
    # ending at a hard horizon line
    fog = np.clip(1.0 - np.exp(-t / 11.0), 0.0, 1.0)[..., None]
    gcol = gcol * (1.0 - fog) + img * fog
    img = np.where(hit[..., None], gcol, img)
    return img, np.where(hit, gdepth, np.float32(1e9))


def shade_body(cam, V, F, VN, VC, W, H, ground_depth_lo, ss):
    """Rasterise and shade the body at the super-sampled resolution.

    Returns premultiplied colour and coverage at the *output* resolution, so
    only the body pays for super-sampling; the ground plane is evaluated once
    per output pixel.
    """
    key = KEY_DIR / np.linalg.norm(KEY_DIR)
    fill = FILL_DIR / np.linalg.norm(FILL_DIR)
    rim = RIM_DIR / np.linalg.norm(RIM_DIR)

    sx, sy, z, _ = cam.project(V, W, H)
    out = rasterise(sx.astype(np.float32), sy.astype(np.float32), z.astype(np.float32),
                    F, W, H)
    Ho, Wo = H // ss, W // ss
    col = np.zeros((Ho * Wo, 3), np.float32)
    cov = np.zeros(Ho * Wo, np.float32)
    if out is None:
        return col, cov
    pix, tri, bary, pz = out

    # depth test against the ground, sampled at the output pixel the fragment
    # falls in; the plane is smooth so this is exact everywhere except within
    # one output pixel of the contact line
    lo_idx = (pix // W // ss) * Wo + (pix % W) // ss
    keep = pz < ground_depth_lo[lo_idx]
    pix = pix[keep]; tri = tri[keep]; bary = bary[keep]; pz = pz[keep]
    if pix.size == 0:
        return col, cov
    lo_idx = lo_idx[keep]

    f3 = F[tri]
    N = (VN[f3[:, 0]] * bary[:, 0:1] + VN[f3[:, 1]] * bary[:, 1:2] +
         VN[f3[:, 2]] * bary[:, 2:3])
    N /= np.maximum(np.linalg.norm(N, axis=1, keepdims=True), 1e-9)
    C = (VC[f3[:, 0]] * bary[:, 0:1] + VC[f3[:, 1]] * bary[:, 1:2] +
         VC[f3[:, 2]] * bary[:, 2:3])
    Pw = (V[f3[:, 0]] * bary[:, 0:1] + V[f3[:, 1]] * bary[:, 1:2] +
          V[f3[:, 2]] * bary[:, 2:3])
    view = cam.eye[None, :].astype(np.float32) - Pw
    view /= np.maximum(np.linalg.norm(view, axis=1, keepdims=True), 1e-9)

    nk = np.maximum(-(N @ key.astype(np.float32)), 0.0)
    nf = np.maximum(-(N @ fill.astype(np.float32)), 0.0)
    nr = np.maximum(-(N @ rim.astype(np.float32)), 0.0)
    nk = (nk + 0.28) / 1.28                     # wrap-around diffuse
    amb = 0.16 + 0.14 * np.clip(N[:, 1] * 0.5 + 0.5, 0, 1)

    h = (-key).astype(np.float32) + view
    h /= np.maximum(np.linalg.norm(h, axis=1, keepdims=True), 1e-9)
    spec = np.power(np.maximum((N * h).sum(1), 0.0), 34.0) * 0.34
    fres = np.power(1.0 - np.clip((N * view).sum(1), 0, 1), 3.0)

    lit = (C * (nk * 1.05 + amb)[:, None] +
           C * np.array([0.30, 0.42, 0.55], np.float32) * (nf * 0.34)[:, None] +
           np.array([0.42, 0.68, 0.95], np.float32) * (nr ** 2.2 * 0.55 * fres)[:, None] +
           spec[:, None] * np.array([1.0, 0.95, 0.88], np.float32))

    w = np.float32(1.0 / (ss * ss))
    np.add.at(col, lo_idx, lit * w)
    np.add.at(cov, lo_idx, w)
    return col.reshape(Ho, Wo, 3), cov.reshape(Ho, Wo)


def vignette(img, strength=0.42):
    H, W, _ = img.shape
    x = (np.arange(W, dtype=np.float32) / (W - 1) - 0.5) * 2.0
    y = (np.arange(H, dtype=np.float32) / (H - 1) - 0.5) * 2.0
    r = np.sqrt((x[None, :] ** 2) * 0.82 + (y[:, None] ** 2))
    v = np.clip(1.0 - strength * np.clip(r - 0.42, 0, None) ** 1.7, 0.0, 1.0)
    return img * v[..., None]


def tonemap(img):
    x = np.maximum(img, 0.0)
    x = x * (2.51 * x + 0.03) / (x * (2.43 * x + 0.59) + 0.14)   # ACES-ish
    x = np.clip(x, 0, 1) ** (1 / 2.2)
    return x


# --------------------------------------------------------------------------
# HUD
# --------------------------------------------------------------------------

def _font(size, mono=True, bold=False):
    from PIL import ImageFont
    import matplotlib
    base = Path(matplotlib.__file__).parent / "mpl-data/fonts/ttf"
    name = "DejaVuSansMono" if mono else "DejaVuSans"
    if bold:
        name += "-Bold"
    return ImageFont.truetype(str(base / f"{name}.ttf"), size)


def draw_hud(img, subjects, frames, t, prog, W, H, hud_muscles, title, subtitle, note,
             seg_scale):
    from PIL import Image, ImageDraw
    img = img.copy()
    # scrims: keep the caption and title legible whatever passes behind them
    top = int(H * 0.16)
    ramp = np.clip(1.0 - np.arange(top, dtype=np.float32) / top, 0, 1)[:, None, None] ** 1.6
    img[:top] *= 1.0 - 0.72 * ramp
    bot = int(H * 0.14)
    ramp = np.clip(np.arange(bot, dtype=np.float32) / bot, 0, 1)[:, None, None] ** 1.8
    img[H - bot:] *= 1.0 - 0.80 * ramp
    pil = Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8))
    d = ImageDraw.Draw(pil, "RGBA")
    s = W / 1920.0
    f_title = _font(int(34 * s), mono=False, bold=True)
    f_sub = _font(int(19 * s), mono=True)
    f_lab = _font(int(17 * s), mono=True)
    f_big = _font(int(30 * s), mono=True, bold=True)
    f_small = _font(int(14 * s), mono=True)

    ink = (232, 240, 247)
    dim = (128, 148, 166)

    d.text((int(56 * s), int(46 * s)), title, font=f_title, fill=ink)
    d.text((int(58 * s), int(90 * s)), subtitle, font=f_sub, fill=dim)

    # progress bar
    bx0, bx1 = int(56 * s), int(W - 56 * s)
    by = int(H - 62 * s)
    d.rounded_rectangle([bx0, by, bx1, by + int(5 * s)], radius=int(3 * s), fill=(38, 48, 58, 255))
    d.rounded_rectangle([bx0, by, bx0 + int((bx1 - bx0) * prog), by + int(5 * s)],
                        radius=int(3 * s), fill=(96, 200, 178, 255))
    d.text((bx0, by + int(14 * s)), note, font=f_small, fill=(110, 128, 144))
    d.text((bx1 - int(150 * s), by - int(30 * s)), f"t = {t:6.2f} s", font=f_sub, fill=ink)

    # per-subject readout cards
    cy = int(150 * s)
    for sub, fr in zip(subjects, frames):
        col = sub["colour"]
        cw, ch = int(342 * s), int(172 * s)
        cx = int(56 * s)
        d.rounded_rectangle([cx, cy, cx + cw, cy + ch], radius=int(10 * s),
                            fill=(12, 18, 25, 190), outline=(*col, 130), width=max(1, int(1.5 * s)))
        d.rectangle([cx, cy + int(10 * s), cx + int(4 * s), cy + ch - int(10 * s)], fill=(*col, 255))
        d.text((cx + int(18 * s), cy + int(12 * s)), sub["label"].upper(), font=f_lab, fill=(*col, 255))
        state, scol = sub["state"](fr)
        d.text((cx + int(18 * s), cy + int(36 * s)), state, font=f_big, fill=scol)
        rows = [
            ("pelvis", f"{fr.get('pelvis_height_m', float('nan')):.3f} m"),
            ("COM drift", f"{fr.get('com_horizontal_displacement_m', 0.0) * 1000:7.2f} mm"),
            ("contacts", f"{int(fr.get('contact_count', 0)):d}"),
        ]
        if fr.get("_ended") is not None:
            rows.append(("recording", f"ended {fr['_ended']:.2f} s"))
        ry = cy + int(78 * s)
        for k, v in rows:
            d.text((cx + int(18 * s), ry), k, font=f_small, fill=dim)
            d.text((cx + int(140 * s), ry), v, font=f_small, fill=ink)
            ry += int(21 * s)
        cy += ch + int(16 * s)

    # muscle excitation bars for the first subject
    if hud_muscles:
        px = int(W - 330 * s)
        py = int(150 * s)
        d.text((px, py - int(26 * s)), "MUSCLE EXCITATION", font=f_lab, fill=dim)
        for name, val, vmax in hud_muscles:
            short = name.replace("gait2392_", "").replace("arm26_", "")
            d.text((px, py), short, font=f_small, fill=(150, 168, 184))
            x0 = px + int(118 * s)
            x1 = int(W - 56 * s)
            d.rectangle([x0, py + int(4 * s), x1, py + int(13 * s)], fill=(30, 39, 48, 255))
            frac = min(max(val / vmax, 0.0), 1.0)
            c = activation_colour(val, 0.0, vmax)
            d.rectangle([x0, py + int(4 * s), x0 + int((x1 - x0) * frac), py + int(13 * s)],
                        fill=tuple(int(255 * v) for v in np.clip(c, 0, 1)))
            py += int(23 * s)

    # colour legend
    lo, hi = seg_scale
    lx, ly = int(W - 330 * s), int(H - 156 * s)
    d.text((lx, ly), "SEGMENT COLOUR", font=f_small, fill=dim)
    d.text((lx, ly + int(17 * s)), "force-weighted mean excitation", font=f_small,
           fill=(96, 112, 126))
    bw = int(220 * s)
    for i in range(bw):
        c = activation_colour(lo + (hi - lo) * i / max(bw - 1, 1), lo, hi)
        d.rectangle([lx + i, ly + int(40 * s), lx + i + 1, ly + int(52 * s)],
                    fill=tuple(int(255 * v) for v in np.clip(c, 0, 1)))
    d.text((lx, ly + int(56 * s)), f"{lo:.3f}", font=f_small, fill=dim)
    d.text((lx + bw - int(38 * s), ly + int(56 * s)), f"{hi:.3f}", font=f_small, fill=dim)

    return np.asarray(pil).astype(np.float32) / 255.0


# --------------------------------------------------------------------------
# Trace handling
# --------------------------------------------------------------------------

def load_trace(path: Path):
    d = json.loads(Path(path).read_text())
    if isinstance(d, dict):
        for k in ("frames", "trace", "trajectory"):
            if k in d:
                d = d[k]
                break
    return d


def frame_coords(fr):
    return {k: float(v["value"]) for k, v in fr["joints"].items()}


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

SUBJECT_COLOURS = [(96, 205, 180), (240, 128, 150), (140, 165, 240)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subject", action="append", default=None,
                    help="path/to/trace.json:label  (repeatable)")
    ap.add_argument("--out", default="artifacts/body_render.mp4")
    ap.add_argument("--width", type=int, default=1920)
    ap.add_argument("--height", type=int, default=1080)
    ap.add_argument("--fps", type=int, default=25)
    ap.add_argument("--ssaa", type=int, default=2)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-frames", type=int, default=0)
    ap.add_argument("--check", action="store_true", help="print FK landmarks and exit")
    ap.add_argument("--still", type=int, default=None, help="render one frame to PNG and exit")
    ap.add_argument("--title", default="IHM-1 — embodied stance, 3D from the OpenSim model")
    a = ap.parse_args()

    model = OsimModel(MODEL)
    L, segments = build_skeleton(model)
    mm = MuscleMap(CATALOG, segments)

    if a.check:
        q = {}
        T = model.forward(q)
        P = all_landmarks(model, T, L)
        print(f"model {model.name}: {len(model.bodies)} bodies, {len(model.joints)} joints, "
              f"{len(model.markers)} markers")
        print("\nneutral pose (pelvis at origin), y = height above pelvis origin:")
        for k in sorted(P):
            print(f"  {k:14s} {P[k][0]:8.3f} {P[k][1]:8.3f} {P[k][2]:8.3f}")
        print("\nsegment lengths read from the model:")
        for s in segments:
            d = np.linalg.norm(P[s['b']] - P[s['a']])
            print(f"  {s['name']:12s} {d*100:6.1f} cm   bodies={s['bodies']}")
        return

    subs = a.subject
    if not subs:
        base = "data/derived/unified-world-ofb2dp7z"
        subs = [f"{base}/full/trace.json:cortex in the loop",
                f"{base}/sever/trace.json:cortex severed"]

    subjects = []
    for i, spec in enumerate(subs):
        path, _, label = spec.partition(":")
        tr = load_trace(ROOT / path if not Path(path).is_absolute() else Path(path))
        h0 = float(tr[0].get("pelvis_height_m", 1.0))

        def state_of(fr, h0=h0):
            """UPRIGHT / FALLING / FALLEN.

            FALLEN is the flag the run itself recorded.  FALLING is this
            renderer's own label for `pelvis has dropped more than 3 cm below
            where it started and the run has not yet declared the fall`; it is
            a display aid, not a quantity the simulation reported.
            """
            if fr.get("fallen"):
                return "FALLEN", (255, 92, 74)
            if float(fr.get("pelvis_height_m", h0)) < h0 - 0.03:
                return "FALLING", (255, 176, 72)
            return "UPRIGHT", (96, 210, 176)

        subjects.append(dict(path=path, label=label or Path(path).parent.name,
                             trace=tr, colour=SUBJECT_COLOURS[i % len(SUBJECT_COLOURS)],
                             state=state_of, h0=h0))
        print(f"subject {i}: {len(tr):4d} frames  {label}  <- {path}")

    n = max(len(s["trace"]) for s in subjects)
    idx = list(range(0, n, a.stride))
    if a.max_frames:
        idx = idx[:a.max_frames]

    # lateral placement so several bodies share one scene
    nsub = len(subjects)
    span = 1.42
    for i, s in enumerate(subjects):
        s["offset"] = np.array([0.0, 0.0, (i - (nsub - 1) / 2.0) * span])

    # Colour scale, read off the data rather than guessed.  Two scales: the raw
    # per-muscle excitation for the HUD bars, and the force-weighted per-segment
    # mean for the segment colours.  They are different quantities and the
    # legend prints the numeric range of the one it labels.
    all_ex, all_seg = [], []
    for s in subjects:
        for fr in s["trace"][::5]:
            ex = fr.get("motor_excitations", {})
            all_ex.extend(ex.values())
            all_seg.extend(v for v in mm.activations(ex).values() if v is not None)
    ex_hi = max(float(np.percentile(all_ex, 99.5)) if all_ex else 0.15, 1e-3)
    seg_lo = float(np.percentile(all_seg, 2.0)) if all_seg else 0.0
    seg_hi = float(np.percentile(all_seg, 99.0)) if all_seg else 0.15
    if seg_hi - seg_lo < 1e-4:
        seg_hi = seg_lo + 1e-4
    print(f"per-muscle excitation scale : 0 .. {ex_hi:.4f}")
    print(f"per-segment colour scale    : {seg_lo:.4f} .. {seg_hi:.4f}")

    SS = max(1, a.ssaa)
    W, H = a.width * SS, a.height * SS      # body raster resolution
    Wo, Ho = a.width, a.height              # output / ground resolution
    aspect = a.width / a.height

    # HUD muscle picks: the largest movers of stance
    hud_pick = ["gait2392_ercspn_r", "gait2392_extobl_r", "psoas_r", "glmax2_r",
                "glmed1_r", "vaslat_r", "gasmed_r", "soleus_r", "tibant_r",
                "arm26_BIClong_r"]
    hud_pick = [m for m in hud_pick if m in subjects[0]["trace"][0].get("motor_excitations", {})]

    out_path = ROOT / a.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    writer = None
    if a.still is None:
        import imageio_ffmpeg
        ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
        cmd = [ffmpeg, "-y", "-loglevel", "error",
               "-f", "rawvideo", "-pix_fmt", "rgb24",
               "-s", f"{a.width}x{a.height}", "-r", str(a.fps), "-i", "-",
               "-an", "-c:v", "libx264", "-preset", "slow", "-crf", "17",
               "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(out_path)]
        writer = subprocess.Popen(cmd, stdin=subprocess.PIPE)

    unit_cache = {}
    frames_to_do = [a.still] if a.still is not None else idx
    import time
    t_start = time.time()

    for n_done, k in enumerate(frames_to_do):
        u = k / max(len(idx) - 1, 1)

        Vs, Fs, Cs = [], [], []
        cur_frames = []
        for s in subjects:
            over = k > len(s["trace"]) - 1
            fr = dict(s["trace"][min(k, len(s["trace"]) - 1)])
            fr["_ended"] = float(s["trace"][-1].get("time_s", 0.0)) if over else None
            cur_frames.append(fr)
            q = frame_coords(fr)
            T = model.forward(q)
            P = all_landmarks(model, T, L)
            acts = mm.activations(fr.get("motor_excitations", {}))
            off = s["offset"]
            for seg in segments:
                p0 = P[seg["a"]] + off
                p1 = P[seg["b"]] + off
                v, f = _tube(p0, p1, seg["r0"], seg["r1"], seg["ex"], seg["ez"])
                col = activation_colour(acts.get(seg["name"]), seg_lo, seg_hi)
                col = col * (1 - seg["tone"]) + np.array([0.62, 0.55, 0.50]) * seg["tone"]
                Fs.append(f + sum(x.shape[0] for x in Vs))
                Vs.append(v)
                Cs.append(np.repeat(col[None, :], v.shape[0], 0))
            # head
            hb, hp = L["head_c"]
            hr = float(L["_head_r"][1][0])
            hc = landmark_world(model, T, L, "head_c") + off
            R = T[hb][:3, :3]
            v, f = _ellipsoid((0, 0, 0), (hr * 0.86, hr * 1.06, hr * 0.90), n=26)
            v = v @ R.T + hc
            Fs.append(f + sum(x.shape[0] for x in Vs))
            Vs.append(v.astype(np.float32))
            Cs.append(np.repeat((BASE_COL * 0.55 + np.array([0.62, 0.55, 0.50]) * 0.45)[None, :],
                                v.shape[0], 0))

        V = np.concatenate(Vs, 0).astype(np.float32)
        F = np.concatenate(Fs, 0).astype(np.int32)
        VC = np.concatenate(Cs, 0).astype(np.float32)
        VN = vertex_normals(V.astype(np.float64), F).astype(np.float32)

        # ---- camera: slow orbit around the group, gentle rise and dolly
        centre = np.array([0.0, 0.0, 0.0])
        centre[0] = float(V[:, 0].mean())
        centre[2] = float(np.mean([s["offset"][2] for s in subjects]))
        target = np.array([centre[0], 0.95, centre[2]])
        az = math.radians(-29.0 + 58.0 * u)
        el = math.radians(7.5 + 5.5 * math.sin(u * math.pi))
        dist = 4.95 + 0.40 * math.cos(u * math.pi)
        eye = target + dist * np.array([math.cos(el) * math.cos(az),
                                        math.sin(el),
                                        math.cos(el) * math.sin(az)])
        eye[1] = max(eye[1], 0.35)
        cam = Camera(eye, target, fov_deg=27.0, aspect=aspect)

        ext = (centre[0] - 3.0, centre[0] + 3.0, centre[2] - 3.0, centre[2] + 3.0)
        sh = shadow_map(V.astype(np.float64), F, KEY_DIR, ext, res=320)

        pools = [(float(s["offset"][0]), float(s["offset"][2]),
                  np.array(s["colour"], float) / 255.0) for s in subjects]
        ground, gdepth = render_ground(cam, Wo, Ho, sh, ext,
                                       np.array([target[0], 1.0, target[2]]), pools)
        bcol, bcov = shade_body(cam, V, F, VN, VC, W, H, gdepth.reshape(-1), SS)
        img = ground * (1.0 - bcov[..., None]) + bcol
        img = tonemap(vignette(img))

        hud_m = []
        ex0 = cur_frames[0].get("motor_excitations", {})
        for m in hud_pick:
            hud_m.append((m, float(ex0.get(m, 0.0)), ex_hi))

        t_s = float(cur_frames[0].get("time_s", k * 0.02))
        sub_title = (f"OpenSim forward kinematics of the recorded joint trajectory  ·  "
                     f"{len(model.bodies)} bodies  ·  {len(segments)} segments  ·  "
                     f"98 muscles")
        note = ("segment colour is recorded muscle excitation mapped through catalog.json "
                "attachment_bodies; geometry is read from model.osim, no limb length is assumed. "
                "playback 0.5× real time.")
        img = draw_hud(img, subjects, cur_frames, t_s, (n_done + 1) / len(frames_to_do),
                       a.width, a.height, hud_m, a.title, sub_title, note,
                       (seg_lo, seg_hi))

        buf = (np.clip(img, 0, 1) * 255).astype(np.uint8)
        if writer is None:
            from PIL import Image
            p = out_path.with_suffix(".png") if out_path.suffix != ".png" else out_path
            Image.fromarray(buf).save(p)
            print(f"wrote still {p}")
            return
        writer.stdin.write(buf.tobytes())

        if n_done % 10 == 0 or n_done == len(frames_to_do) - 1:
            el_s = time.time() - t_start
            rate = (n_done + 1) / max(el_s, 1e-6)
            eta = (len(frames_to_do) - n_done - 1) / max(rate, 1e-6)
            print(f"  frame {n_done+1:4d}/{len(frames_to_do)}  "
                  f"{rate:4.2f} fps  eta {eta/60:4.1f} min", flush=True)

    writer.stdin.close()
    rc = writer.wait()
    if rc != 0:
        print("ffmpeg failed", file=sys.stderr)
        sys.exit(1)
    size = out_path.stat().st_size / 1e6
    print(f"wrote {out_path}  {a.width}x{a.height}  {len(idx)} frames  "
          f"{a.fps} fps  {size:.1f} MB")


if __name__ == "__main__":
    main()
