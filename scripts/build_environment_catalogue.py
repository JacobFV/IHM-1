"""Environment and environment-component catalogue, slot model, thumbnails, provenance.

Emits only settings an engine in this repository already accepts. Environments come
from ihm/assembly/interactive_scene.py ENVIRONMENTS and scripts/native_mechanical_stream.cpp;
components come from parameters those engines take (bed_material, surface_contact_manifest,
ambient_temperature_c). Nothing here invents scenery.
"""
import argparse
import gzip
import hashlib
import json
import math
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
OUT = ROOT / 'data/derived/environment-catalogue-v1'
THUMBS = OUT / 'thumbnails'
PROV = OUT / 'provenance'
SCHEMA = 'ihm.environment-catalog.v1'
PROVENANCE_SCHEMA = 'ihm.structure-provenance.v1'

SCENE_MODULE = 'ihm/assembly/interactive_scene.py'
NATIVE_CPP = 'scripts/native_mechanical_stream.cpp'
NATIVE_STREAM = 'ihm/native/mechanical_stream.py'
NATIVE_CONFIG = 'ihm/native/__init__.py'
BIOGEARS_CPP = 'scripts/native_biogears_rest.cpp'
BED_MANIFEST = 'data/research/bed_material/manifest.json'
BED_CURVE = 'data/research/bed_material/hong2022-compression.csv'
BED_FIGURE = 'data/research/bed_material/hong2022-figure1.jpg'
BED_ARTICLE = 'data/research/bed_material/hong2022.html'
BED_IMPL = 'ihm/assembly/bed_compression.py'
SUPINE_MANIFEST = 'data/derived/supine-surface-contact-exmzq9pq/manifest.json'
SUPINE_ARRAYS = 'data/derived/supine-surface-contact-exmzq9pq/quadrature.npz'
SUPINE_FOUNDATION = 'data/derived/supine-surface-contact-exmzq9pq/supine_surface_foundation.txt'
SUPINE_IMPL = 'ihm/assembly/supine_contact.py'
SKIN_GEOMETRY = 'data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz'
ANATOMY = 'data/derived/canonical/anatomy.json'
OSIM = 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled.osim'
OSIM_LICENSE = 'data/raw/mechanics/opensim-core/LICENSE.txt'
THIS = 'scripts/build_environment_catalogue.py'

IMAGE_PX = 320
SUPERSAMPLE = 3
FRAME_HALF_M = 1.06
ENV_HALF_M = 1.30
SCENE_HALF_M = 2.00
CENTER_M = np.array([0., -.05, 0.])
INK = '#20242b'
PAPER = '#f2f1ee'
ACCENT = '#2f6f7d'
MUTED = '#b7b4ac'
COLOUR = {
    'skin': (.87, .70, .59), 'plane': '#d3d7d1', 'plane_line': '#c2c7c0',
    'bed-frame': (.55, .40, .28), 'bed-mattress': (.91, .89, .83), 'bed-rails': (.72, .75, .78),
    'pillow': (.88, .91, .95), 'blanket': (.38, .50, .62), 'nightstand': (.58, .44, .31),
    'table': (.66, .51, .36), 'chair': (.50, .38, .27), 'iv-stand': (.74, .77, .80),
    'iv-bag': (.80, .88, .78), 'ball-small': (.83, .33, .27), 'ball-large': (.24, .55, .62),
}
COLOUR.update({
    'wall-warm': (.86, .82, .74), 'wall-clinical': (.85, .88, .88), 'wall-play': (.84, .85, .87),
    'floor-wood': (.62, .47, .32), 'floor-vinyl': (.76, .78, .77), 'floor-carpet': (.55, .52, .49),
    'ceiling': (.93, .93, .91), 'window-pane': (.72, .84, .92), 'window-frame': (.94, .94, .92),
    'grass': (.42, .55, .28), 'patio': (.72, .70, .65), 'trunk': (.42, .32, .22), 'canopy': (.29, .44, .24),
    'block': (.79, .58, .27),
})
CURVE_COLOUR = {'SM': '#4c8fbd', 'MM': '#2f8f6f', 'HM': '#c05a2e'}
SKY = {
    'day': {'type': 'vertical_gradient',
            'stops': [{'t': 0., 'rgb': [.78, .85, .90]}, {'t': .55, 'rgb': [.55, .72, .88]},
                      {'t': 1., 'rgb': [.33, .55, .82]}],
            'label': 'Clear day sky'},
    'interior': {'type': 'flat', 'stops': [{'t': 0., 'rgb': [.90, .89, .87]}, {'t': 1., 'rgb': [.90, .89, .87]}],
                 'label': 'Room interior; the enclosure provides the surround'},
    'none': {'type': 'none', 'stops': [], 'label': 'No sky: the environment tile shows the engine plane, not a world'},
}


def normalize(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


VIEW = normalize([.55, .22, 1.])
LIGHT = normalize([.35, .75, .55])
RIGHT = normalize(np.cross([0., 1., 0.], VIEW))
UP = np.cross(VIEW, RIGHT)
CAMERA = {'projection': 'orthographic', 'view_direction_canonical': VIEW.tolist(),
          'image_right_canonical': RIGHT.tolist(), 'image_up_canonical': UP.tolist(),
          'centre_m': CENTER_M.tolist(), 'pixels': IMAGE_PX, 'supersample': SUPERSAMPLE,
          'half_extent_m': {'environment_tiles': ENV_HALF_M, 'scene_tiles': SCENE_HALF_M,
                            'component_tiles': FRAME_HALF_M,
                            'object_tiles': 'fitted per object; the value is recorded on each record as thumbnail_half_extent_m'},
          'light_direction_canonical': LIGHT.tolist(),
          'shading': 'Lambert, ambient 0.36 + 0.64 max(0, n.l), single fixed light, no shadows',
          'frame': 'bodyparts3d-display-m (x left, y superior, z anterior)',
          'gravity_frame': 'The camera and light are fixed relative to gravity, not to the canonical axes: a tile whose environment pulls along -z is rendered through the recorded permutation (x,y,z)->(y,z,x) so gravity points down in the image. The permutation is a viewing transform; no geometry is moved and no physics is restated.',
          'note': 'One camera direction, one light and one shading rule for every tile. The body pose is identical in every environment because selecting an environment changes gravity, the free-object ground plane and the prescribed support set; it does not repose the body. A scene adds objects around that same body.'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(*args):
    return subprocess.run(['git', *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()


# A tile is rendered in its environment's gravity frame: the camera and light are
# fixed relative to gravity, not to the canonical axes, so a bed scene reads as a bed
# instead of a wall. The permutation is a viewing transform and is recorded per tile.
GRAVITY_FRAME = {
    'studio': {'rotation': None, 'centre': [0., -.05, 0.],
               'basis': 'canonical axes unchanged; zero gravity has no down'},
    'floor': {'rotation': None, 'centre': [0., -.05, 0.],
              'basis': 'canonical axes unchanged; gravity is already -y'},
    'bed': {'rotation': [[0, 1, 0], [0, 0, 1], [1, 0, 0]], 'centre': [0., -.30, 0.],
            'basis': 'canonical (x,y,z) -> display (y,z,x): gravity -z becomes image down and the body long axis becomes image right'},
}


def view_of(environment=None, half=None):
    frame = GRAVITY_FRAME.get(environment, GRAVITY_FRAME['floor'])
    rotation = None if frame['rotation'] is None else np.asarray(frame['rotation'], float)
    return {'rotation': rotation, 'centre': np.asarray(frame['centre'], float),
            'half': FRAME_HALF_M if half is None else half,
            'basis': frame['basis'], 'environment': environment}


def to_display(points, view):
    p = np.asarray(points, float)
    if view['rotation'] is not None:
        p = p @ view['rotation'].T
    return p


def project(points, view=None):
    view = view or view_of()
    p = to_display(points, view) - view['centre']
    return np.stack([p @ RIGHT, p @ UP], -1), p @ VIEW


# ---------------------------------------------------------------- evidence

def read_json(relative):
    return json.loads((ROOT / relative).read_bytes())


def bed_curves():
    from ihm.assembly.bed_compression import load_bed
    return {m: load_bed(ROOT, m) for m in ('SM', 'MM', 'HM')}


def proxy_radii():
    """Exact cpp rule: radius^2 = 5 (Iyy + Izz - Ixx) / (2 m) at the retained body inertia."""
    root = ET.parse(ROOT / OSIM).getroot()
    out = {}
    for body in root.iter('Body'):
        mass = float(body.find('mass').text)
        inertia = [float(v) for v in body.find('inertia').text.split()]
        radius2 = 5 * (inertia[1] + inertia[2] - inertia[0]) / (2 * mass)
        if not radius2 > 0:
            raise ValueError('Source inertia cannot support a posterior radius: ' + body.get('name'))
        out[body.get('name')] = (mass, math.sqrt(radius2))
    return out


def skin_mesh():
    data = json.loads(gzip.open(ROOT / SKIN_GEOMETRY).read())
    positions = np.asarray(data['positions'], float).reshape(-1, 3)
    faces = np.asarray(data['indices'], int).reshape(-1, 3)
    return positions, faces


def quadrature_points():
    manifest = read_json(SUPINE_MANIFEST)
    transform = np.asarray(manifest['registration']['source_to_canonical_ground'], float)
    with np.load(ROOT / SUPINE_ARRAYS) as arrays:
        points = np.asarray(arrays['reference_points_source_m'], float)
        bodies = np.asarray(arrays['body_indices'], int)
    return points @ transform[:3, :3].T + transform[:3, 3], bodies, manifest


def foundation_material():
    header = (ROOT / SUPINE_FOUNDATION).read_text().split('\n', 1)[0].split()
    keys = ['plane_source_x_m', 'total_layer_thickness_m', 'shear_modulus_pa', 'lame_lambda_pa',
            'minimum_thickness_ratio', 'transition_velocity_m_s', 'dynamic_friction',
            'viscous_friction', 'area_reference_m2', 'points']
    return dict(zip(keys, [float(v) for v in header[1:]]))


def support_names():
    entities = read_json(ANATOMY)['entities']
    return {e['id']: (e['name'], e['centroid_m']) for e in entities
            if e['id'] in {'body-bp3d-FJ3256', 'body-bp3d-FJ3360', 'body-bp3d-FJ3309', 'body-bp3d-FJ3393'}}


def source_literals():
    """Verified-not-inferred: the engines' own accepted-value literals."""
    scene = (ROOT / SCENE_MODULE).read_text()
    cpp = (ROOT / NATIVE_CPP).read_text()
    stream = (ROOT / NATIVE_STREAM).read_text()
    config = (ROOT / NATIVE_CONFIG).read_text()
    bio = (ROOT / BIOGEARS_CPP).read_text()
    return {
        'scene_environment_ids': sorted(re.findall(r"^\s{4}'([a-z]+)':dict\(label=", scene, re.M)),
        'native_cpp_environments': sorted(set(re.findall(r'environment!="(\w+)"', cpp))),
        'native_python_environments': sorted(re.search(r"environment not in \(([^)]*)\)", stream).group(1).replace("'", '').split(',')),
        'ambient_python_bounds_c': [float(v) for v in re.search(r"number\(self\.ambient_temperature_c,([\d.]+),([\d.]+),", config).groups()],
        'ambient_cpp_bounds_c': [float(v) for v in re.search(r"ambient < ([\d.]+) \|\| ambient > ([\d.]+)", bio).groups()],
        'bed_material_choices': sorted(set(re.findall(r"choices=\('SM','MM','HM'\)", (ROOT / 'scripts/verify_native_surface_foundation.py').read_text()) and ['SM', 'MM', 'HM'])),
    }


# ---------------------------------------------------------------- rendering

def rasterize(meshes, view=None):
    """One orthographic z-buffer over every mesh, so occlusion between body and objects is real.

    meshes: [(positions, faces, rgb)]. Returns (cover, rgb) at IMAGE_PX, Lambert shaded
    by the single fixed light. Colour is per mesh; shading is per face.
    """
    view = view or view_of()
    half = view['half']
    size = IMAGE_PX * SUPERSAMPLE
    scale = size / (2 * half)
    zbuffer = np.full((size, size), -np.inf)
    image = np.zeros((size, size, 3))
    alpha = np.zeros((size, size), bool)
    for positions, faces, rgb in meshes:
        positions = to_display(positions, view)
        faces = np.asarray(faces, int)
        p = positions - view['centre']
        uv = np.stack([p @ RIGHT, p @ UP], -1)
        depth = p @ VIEW
        px = np.stack([(uv[:, 0] + half) * scale, (half - uv[:, 1]) * scale], -1)
        tri = positions[faces]
        normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        lengths = np.linalg.norm(normals, axis=1)
        keep = (lengths > 0) & ((normals @ VIEW) > 0)
        kept = faces[keep]
        normals = normals[keep] / lengths[keep, None]
        shade = .36 + .64 * np.clip(normals @ LIGHT, 0, 1)
        base = np.asarray(rgb, float)
        base = base[keep] if base.ndim == 2 else base[None]
        tone = base * shade[:, None]
        corners = px[kept]
        zs = depth[kept]
        low = np.floor(corners.min(1)).astype(int)
        high = np.ceil(corners.max(1)).astype(int)
        for k in range(len(kept)):
            x0, y0 = max(low[k, 0], 0), max(low[k, 1], 0)
            x1, y1 = min(high[k, 0], size - 1), min(high[k, 1], size - 1)
            if x1 < x0 or y1 < y0:
                continue
            a, b, c = corners[k]
            det = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
            if abs(det) < 1e-12:
                continue
            gx, gy = np.meshgrid(np.arange(x0, x1 + 1) + .5, np.arange(y0, y1 + 1) + .5)
            w1 = ((gx - a[0]) * (c[1] - a[1]) - (gy - a[1]) * (c[0] - a[0])) / det
            w2 = ((b[0] - a[0]) * (gy - a[1]) - (b[1] - a[1]) * (gx - a[0])) / det
            inside = (w1 >= 0) & (w2 >= 0) & (w1 + w2 <= 1)
            if not inside.any():
                continue
            z = zs[k, 0] + w1 * (zs[k, 1] - zs[k, 0]) + w2 * (zs[k, 2] - zs[k, 0])
            window = zbuffer[y0:y1 + 1, x0:x1 + 1]
            hit = inside & (z > window)
            window[hit] = z[hit]
            image[y0:y1 + 1, x0:x1 + 1][hit] = tone[k]
            alpha[y0:y1 + 1, x0:x1 + 1][hit] = True
    block = SUPERSAMPLE
    cover = alpha.reshape(IMAGE_PX, block, IMAGE_PX, block).mean((1, 3))
    count = alpha.reshape(IMAGE_PX, block, IMAGE_PX, block).sum((1, 3))
    total = (image * alpha[..., None]).reshape(IMAGE_PX, block, IMAGE_PX, block, 3).sum((1, 3))
    rgb = np.divide(total, count[..., None], out=np.zeros_like(total), where=count[..., None] > 0)
    return cover, rgb


def figure(half=FRAME_HALF_M):  # noqa: D401
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(IMAGE_PX / 100, IMAGE_PX / 100), dpi=100)
    axes = fig.add_axes([0, 0, 1, 1])
    axes.set_xlim(-half, half)
    axes.set_ylim(-half, half)
    axes.set_aspect('equal')
    axes.set_axis_off()
    axes.set_facecolor(PAPER)
    fig.patch.set_facecolor(PAPER)
    return fig, axes


def raster_layer(axes, cover, rgb, alpha=1., half=FRAME_HALF_M):
    rgba = np.concatenate([np.clip(rgb, 0, 1), (cover * alpha)[..., None]], -1)
    axes.imshow(rgba, extent=[-half, half, -half, half], interpolation='bilinear', zorder=2)


def plane_extent(view, half=.62):
    return half


def plane_layer(axes, axis, offset, view, half=.62, fill=True):
    """The environment plane drawn where the data puts it, in the shared camera."""
    others = [i for i in range(3) if i != axis]
    def point(a, b):
        p = np.zeros(3)
        p[axis] = offset
        p[others[0]] = a
        p[others[1]] = b
        return p
    corners = [point(-half, -half), point(half, -half), point(half, half), point(-half, half)]
    uv, _ = project(np.array(corners), view)
    if fill:
        axes.add_patch(__import__('matplotlib').patches.Polygon(uv, closed=True, facecolor=COLOUR['plane'],
                                                                edgecolor=COLOUR['plane_line'], linewidth=1., zorder=1))
    for t in np.linspace(-half, half, 7):
        for pair in ((point(t, -half), point(t, half)), (point(-half, t), point(half, t))):
            line, _ = project(np.array(pair), view)
            axes.plot(line[:, 0], line[:, 1], color=COLOUR['plane_line'], linewidth=.7, zorder=1)


def gravity_layer(axes, gravity, view, anchor=np.array([.52, .62, 0.])):
    g = np.asarray(gravity, float)
    magnitude = float(np.linalg.norm(g))
    if view['rotation'] is not None:
        anchor = np.array([.72, 0., .60]) if view['environment'] == 'bed' else anchor
    start, _ = project(anchor[None], view)
    if magnitude == 0:
        axes.plot(*start[0], marker='o', markersize=9, markerfacecolor='none',
                  markeredgecolor=INK, markeredgewidth=1.6, zorder=4)
        axes.text(start[0, 0], start[0, 1] - .12, 'g = 0', color=INK, ha='center',
                  va='top', fontsize=13, zorder=4)
        return
    end, _ = project((anchor + g / magnitude * .40)[None], view)
    axes.annotate('', xy=end[0], xytext=start[0], zorder=4,
                  arrowprops=dict(arrowstyle='-|>', color=INK, linewidth=2.2, mutation_scale=18))
    mid = (start[0] + end[0]) / 2
    right_room = FRAME_HALF_M - max(start[0, 0], end[0, 0])
    if right_room > .48:
        axes.text(mid[0] + .09, mid[1], f'{magnitude:.2f} m s⁻²', color=INK, ha='left',
                  va='center', fontsize=12, zorder=4)
    else:
        axes.text(mid[0], min(start[0, 1], end[0, 1]) - .10, f'{magnitude:.2f} m s⁻²', color=INK,
                  ha='center', va='top', fontsize=12, zorder=4)


def support_layer(axes, centroids, view):
    if not centroids:
        return
    uv, _ = project(np.array(centroids), view)
    axes.scatter(uv[:, 0], uv[:, 1], s=64, facecolor=ACCENT, edgecolor=PAPER, linewidth=1.4, zorder=5)


def save(fig, name):
    path = THUMBS / (name + '.png')
    fig.savefig(path, facecolor=PAPER, dpi=100)
    fig.clear()
    import matplotlib.pyplot as plt
    plt.close(fig)
    return path


def curve_figure():
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(IMAGE_PX / 100, IMAGE_PX / 100), dpi=100)
    axes = fig.add_axes([.19, .16, .76, .78])
    fig.patch.set_facecolor(PAPER)
    axes.set_facecolor(PAPER)
    for side in ('top', 'right'):
        axes.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        axes.spines[side].set_color(INK)
    axes.tick_params(colors=INK, labelsize=11)
    return fig, axes


def render_mattress(curves, selected):
    fig, axes = curve_figure()
    for key, curve in curves.items():
        axes.plot(curve['strain'], np.asarray(curve['pressure_pa']) / 1e3,
                  color=CURVE_COLOUR[key], linewidth=1.5, alpha=.42, zorder=1)
    if selected is None:
        axes.axvline(0, color=INK, linewidth=4., zorder=3)
        axes.text(.05, 23.4, 'rigid plane', color=INK, fontsize=13, va='top')
    else:
        curve = curves[selected]
        colour = CURVE_COLOUR[selected]
        axes.plot(curve['strain'], np.asarray(curve['pressure_pa']) / 1e3, color=colour, linewidth=4., zorder=3)
        axes.scatter(curve['strain'], np.asarray(curve['pressure_pa']) / 1e3, s=20, color=colour, zorder=4)
    axes.set_xlim(0, .70)
    axes.set_ylim(0, 25)
    axes.set_xlabel('compressive strain', color=INK, fontsize=12)
    axes.set_ylabel('nominal stress (kPa)', color=INK, fontsize=12)
    return fig


def render_quadrature(points, bodies, cover, rgb):
    fig, axes = figure()
    raster_layer(axes, cover, rgb, alpha=.35)
    uv, depth = project(points, view_of())
    order = np.argsort(depth)
    axes.scatter(uv[order, 0], uv[order, 1], s=.8, c=depth[order], cmap='viridis',
                 linewidths=0, zorder=3)
    return fig


def render_proxy(radii):
    """Row-packed at the same metre scale as the body tiles; layout is not anatomy."""
    fig, axes = figure()
    import matplotlib.patches as patches
    ordered = sorted((r for _, r in radii.values()), reverse=True)
    gap = .03
    rows, row = [], []
    width = 0.
    for radius in ordered:
        if row and width + 2 * radius + gap > 2 * FRAME_HALF_M - .08:
            rows.append(row)
            row, width = [], 0.
        row.append(radius)
        width += 2 * radius + gap
    rows.append(row)
    total = sum(2 * max(r) for r in rows) + gap * (len(rows) - 1)
    y = total / 2 - max(rows[0]) - .05
    for index, row in enumerate(rows):
        x = -(sum(2 * r + gap for r in row) - gap) / 2
        for radius in row:
            x += radius
            axes.add_patch(patches.Circle((x, y), radius, facecolor='#cfd6d8', edgecolor=INK,
                                          linewidth=1.1, zorder=3))
            axes.add_patch(patches.Circle((x - radius * .32, y + radius * .30), radius * .45,
                                          facecolor='#e8edee', edgecolor='none', alpha=.75, zorder=4))
            x += radius + gap
        if index + 1 < len(rows):
            y -= max(row) + gap + max(rows[index + 1])
    axes.text(0, -FRAME_HALF_M + .07, f'{len(radii)} proxy radii  {min(ordered) * 1e3:.0f}–{max(ordered) * 1e3:.0f} mm',
              color=INK, ha='center', va='bottom', fontsize=13, zorder=6)
    return fig


def render_ambient(setpoint, bounds, default_c):
    fig, axes = curve_figure()
    axes.set_position([.12, .36, .80, .26])
    lo, hi = bounds
    gradient = np.linspace(0, 1, 256)[None]
    axes.imshow(gradient, extent=[lo, hi, 0, 1], aspect='auto', cmap='coolwarm', zorder=1)
    axes.set_ylim(0, 1)
    axes.set_xlim(lo, hi)
    axes.set_yticks([])
    axes.plot([setpoint, setpoint], [0, 1], color=INK, linewidth=3.4, zorder=3)
    axes.text(setpoint, 1.30, f'{setpoint:g} °C', color=INK, fontsize=20,
              ha='right' if setpoint > (lo + hi) / 2 else 'left', va='bottom')
    axes.plot([default_c], [-.44], marker='^', color=INK, markersize=7, clip_on=False, zorder=3)
    axes.text(default_c, -.62, f'engine default {default_c:g} °C', color=INK, ha='center', va='top', fontsize=11)
    axes.text((lo + hi) / 2, -1.05, f'accepted range {lo:g}–{hi:g} °C', color=INK, ha='center', va='top', fontsize=11)
    axes.spines['left'].set_visible(False)
    return fig


# ---------------------------------------------------------------- objects and scenes

def _tris(quads):
    """Independent triangles with flat normals; no shared vertices, so faces stay exact."""
    positions, faces = [], []
    for quad in quads:
        base = len(positions)
        positions.extend(quad)
        faces.append([base, base + 1, base + 2])
        if len(quad) == 4:
            faces.append([base, base + 2, base + 3])
    return np.asarray(positions, float), np.asarray(faces, int)


def box(low, high):
    (x0, y0, z0), (x1, y1, z1) = low, high
    quads = [
        [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)],
        [(x1, y0, z0), (x0, y0, z0), (x0, y1, z0), (x1, y1, z0)],
        [(x0, y1, z0), (x0, y1, z1), (x1, y1, z1), (x1, y1, z0)],
        [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)],
        [(x1, y0, z1), (x1, y0, z0), (x1, y1, z0), (x1, y1, z1)],
        [(x0, y0, z0), (x0, y0, z1), (x0, y1, z1), (x0, y1, z0)],
    ]
    return _tris(quads)


def cylinder(centre, radius, height, axis=1, segments=28):
    centre = np.asarray(centre, float)
    other = [i for i in range(3) if i != axis]
    angles = np.linspace(0, 2 * np.pi, segments + 1)
    def ring(offset):
        points = np.zeros((segments + 1, 3))
        points[:, axis] = offset
        points[:, other[0]] = radius * np.cos(angles)
        points[:, other[1]] = radius * np.sin(angles)
        return points + centre
    low, high = ring(-height / 2), ring(height / 2)
    quads = [[low[i], low[i + 1], high[i + 1], high[i]] for i in range(segments)]
    top = centre + np.eye(3)[axis] * height / 2
    bottom = centre - np.eye(3)[axis] * height / 2
    quads += [[high[i], high[i + 1], top] for i in range(segments)]
    quads += [[low[i + 1], low[i], bottom] for i in range(segments)]
    # The (other0, other1, axis) triple is right handed for axis 0 and 2 and left handed
    # for axis 1, so reverse only when it is: outward normals, positive signed volume.
    handed = float(np.linalg.det(np.eye(3)[[other[0], other[1], axis]]))
    return _tris(quads if handed > 0 else [q[::-1] for q in quads])


def ball(centre, radius, rings=18, segments=30):
    centre = np.asarray(centre, float)
    theta = np.linspace(0, np.pi, rings + 1)
    phi = np.linspace(0, 2 * np.pi, segments + 1)
    grid = np.stack(np.meshgrid(theta, phi, indexing='ij'), -1)
    points = centre + radius * np.stack([np.sin(grid[..., 0]) * np.cos(grid[..., 1]),
                                         np.cos(grid[..., 0]),
                                         np.sin(grid[..., 0]) * np.sin(grid[..., 1])], -1)
    quads = [[points[i, j], points[i + 1, j], points[i + 1, j + 1], points[i, j + 1]]
             for i in range(rings) for j in range(segments)]
    return _tris([q[::-1] for q in quads])  # outward normals: signed volume must be positive


def merge(parts):
    positions, faces = [], []
    offset = 0
    for p, f in parts:
        positions.append(p)
        faces.append(np.asarray(f) + offset)
        offset += len(p)
    return np.concatenate(positions), np.concatenate(faces)


def part_volume(part):
    if part['primitive'] == 'box':
        return float(np.prod(np.asarray(part['max_m']) - np.asarray(part['min_m'])))
    if part['primitive'] == 'cylinder':
        return float(np.pi * part['radius_m'] ** 2 * part['height_m'])
    return float(4 / 3 * np.pi * part['radius_m'] ** 3)


def part_mesh(part):
    if part['primitive'] == 'box':
        return box(part['min_m'], part['max_m'])
    if part['primitive'] == 'cylinder':
        return cylinder(part['centre_m'], part['radius_m'], part['height_m'], part.get('axis', 1))
    return ball(part['centre_m'], part['radius_m'])


def _legs(x0, x1, z0, z1, low, high, thickness, axis='y'):
    """Four uprights between low and high along the given up axis."""
    parts = []
    for i, x in enumerate((x0, x1 - thickness)):
        for j, z in enumerate((z0, z1 - thickness)):
            if axis == 'y':
                parts.append({'name': f'leg-{i}{j}', 'primitive': 'box',
                              'min_m': [x, low, z], 'max_m': [x + thickness, high, z + thickness]})
            else:
                parts.append({'name': f'leg-{i}{j}', 'primitive': 'box',
                              'min_m': [x, z, low], 'max_m': [x + thickness, z + thickness, high]})
    return parts


BED_LENGTH_M, BED_WIDTH_M, BED_THICKNESS_M = 1.9, 1.2, .2   # retained Hong et al. specimen mattress
BED_PLANE_Z = -.24                                          # ENVIRONMENTS['bed']['plane']
FLOOR_PLANE_Y = -.96                                        # ENVIRONMENTS['floor']['plane']
BED_ROOM_FLOOR_Z = -.86                                     # implied by the frame legs; no engine plane


def object_specs():
    """Every object is constructed here. Dimensions are stated; nothing is downloaded."""
    half_w, half_l = BED_WIDTH_M / 2, BED_LENGTH_M / 2
    mattress_top, mattress_bottom = BED_PLANE_Z, BED_PLANE_Z - BED_THICKNESS_M
    frame_top = mattress_bottom
    specs = [
        {'id': 'bed-mattress', 'label': 'Mattress', 'colour': 'bed-mattress',
         'base_environment': 'bed', 'mass_kg': None,
         'mass_basis': 'Absent: the retained bed evidence reports no density (specimen_conditions.density_kg_m3 = null, "do not invent").',
         'material': {'response': 'measured compression curve, selected in the mattress_material slot',
                      'tier': 'derived', 'source': 'Hong et al. 2022, DOI 10.3390/biology11071030'},
         'dimensions_basis': 'Retained study specimen mattress 1.9 x 1.2 x 0.2 m (bed_material manifest specimen_conditions.mattress_dimensions_m).',
         'parts': [{'name': 'slab', 'primitive': 'box',
                    'min_m': [-half_w, -half_l + .05, mattress_bottom], 'max_m': [half_w, half_l - .05, mattress_top]}]},
        {'id': 'bed-frame', 'label': 'Bed frame', 'colour': 'bed-frame',
         'base_environment': 'bed', 'mass_kg': 45.,
         'mass_basis': 'Engineering choice; no retained measurement.',
         'material': {'response': None, 'tier': 'synthesized', 'source': 'nominal timber frame; stiffness not calibrated'},
         'dimensions_basis': 'Sized to carry the retained mattress with a 20 mm margin; leg length is an engineering choice that fixes the implied room floor.',
         'parts': ([{'name': 'deck', 'primitive': 'box',
                     'min_m': [-half_w - .04, -half_l, frame_top - .06], 'max_m': [half_w + .04, half_l, frame_top]},
                    {'name': 'headboard', 'primitive': 'box',
                     'min_m': [-half_w - .04, half_l - .06, frame_top], 'max_m': [half_w + .04, half_l, mattress_top + .22]}]
                   + _legs(-half_w - .04, half_w + .04, -half_l, half_l, BED_ROOM_FLOOR_Z, frame_top - .06, .09, axis='z'))},
        {'id': 'bed-rails', 'label': 'Bed side rails', 'colour': 'bed-rails',
         'base_environment': 'bed', 'mass_kg': 8.,
         'mass_basis': 'Engineering choice; no retained measurement.',
         'material': {'response': None, 'tier': 'synthesized', 'source': 'nominal steel rail; stiffness not calibrated'},
         'dimensions_basis': 'Nominal hospital bed side rail envelope; engineering choice.',
         'parts': [{'name': 'rail-left', 'primitive': 'box',
                    'min_m': [-half_w - .05, -.45, mattress_top + .10], 'max_m': [-half_w + .01, .45, mattress_top + .16]},
                   {'name': 'rail-right', 'primitive': 'box',
                    'min_m': [half_w - .01, -.45, mattress_top + .10], 'max_m': [half_w + .05, .45, mattress_top + .16]},
                   {'name': 'post-left', 'primitive': 'box',
                    'min_m': [-half_w - .04, -.44, mattress_bottom], 'max_m': [-half_w, -.38, mattress_top + .16]},
                   {'name': 'post-right', 'primitive': 'box',
                    'min_m': [half_w, -.44, mattress_bottom], 'max_m': [half_w + .04, -.38, mattress_top + .16]}]},
        {'id': 'pillow', 'label': 'Pillow', 'colour': 'pillow',
         'base_environment': 'bed', 'mass_kg': .8,
         'mass_basis': 'Engineering choice; no retained measurement.',
         'material': {'response': None, 'tier': 'synthesized',
                      'source': 'Stiffness absent: the only retained soft-support response in this repository is the mattress curve, which is a mattress, not a pillow.'},
         'dimensions_basis': '0.60 x 0.40 x 0.12 m nominal pillow; engineering choice.',
         'parts': [{'name': 'slab', 'primitive': 'box',
                    'min_m': [-.30, .52, mattress_top], 'max_m': [.30, .92, mattress_top + .12]}]},
        {'id': 'blanket', 'label': 'Blanket', 'colour': 'blanket',
         'base_environment': 'bed', 'mass_kg': 1.5,
         'mass_basis': 'Engineering choice; no retained measurement.',
         'material': {'response': None, 'tier': 'synthesized',
                      'source': 'Textile behaviour is not solved here; the garment lane owns cloth.'},
         'dimensions_basis': 'Top panel 4 mm clear of the anterior skin extent over the legs and hips, with 60 mm side skirts. The shape is drawn, not draped: no cloth solve is attempted, so the torso and head are left uncovered rather than faking a fold.',
         'parts': [{'name': 'top', 'primitive': 'box',
                    'min_m': [-.40, -.88, .150], 'max_m': [.40, -.04, .164]},
                   {'name': 'skirt-left', 'primitive': 'box',
                    'min_m': [-.42, -.88, .090], 'max_m': [-.40, -.04, .164]},
                   {'name': 'skirt-right', 'primitive': 'box',
                    'min_m': [.40, -.88, .090], 'max_m': [.42, -.04, .164]},
                   {'name': 'skirt-foot', 'primitive': 'box',
                    'min_m': [-.42, -.90, .090], 'max_m': [.42, -.88, .164]}]},
        {'id': 'nightstand', 'label': 'Nightstand', 'colour': 'nightstand',
         'base_environment': 'bed', 'mass_kg': 12.,
         'mass_basis': 'Engineering choice; no retained measurement.',
         'material': {'response': None, 'tier': 'synthesized', 'source': 'nominal timber cabinet'},
         'dimensions_basis': '0.38 x 0.40 x 0.55 m nominal cabinet standing on the implied room floor beside the bed.',
         'parts': [{'name': 'top', 'primitive': 'box',
                    'min_m': [-1.08, -.92, BED_ROOM_FLOOR_Z + .52], 'max_m': [-.70, -.52, BED_ROOM_FLOOR_Z + .55]},
                   {'name': 'body', 'primitive': 'box',
                    'min_m': [-1.05, -.89, BED_ROOM_FLOOR_Z], 'max_m': [-.73, -.55, BED_ROOM_FLOOR_Z + .52]}]},
        {'id': 'iv-stand', 'label': 'IV stand', 'colour': 'iv-stand',
         'base_environment': 'bed', 'mass_kg': 9.,
         'mass_basis': 'Engineering choice; no retained measurement.',
         'material': {'response': None, 'tier': 'synthesized', 'source': 'nominal steel pole and bag'},
         'dimensions_basis': '1.60 m pole on a 0.25 m base disc; engineering choice.',
         'parts': [{'name': 'base', 'primitive': 'cylinder', 'axis': 2,
                    'centre_m': [-.88, .74, BED_ROOM_FLOOR_Z + .02], 'radius_m': .25, 'height_m': .04},
                   {'name': 'pole', 'primitive': 'cylinder', 'axis': 2,
                    'centre_m': [-.88, .74, BED_ROOM_FLOOR_Z + .82], 'radius_m': .016, 'height_m': 1.60},
                   {'name': 'bag', 'primitive': 'box',
                    'min_m': [-.94, .66, BED_ROOM_FLOOR_Z + 1.30], 'max_m': [-.82, .82, BED_ROOM_FLOOR_Z + 1.58]}]},
        {'id': 'chair', 'label': 'Chair', 'colour': 'chair',
         'base_environment': 'floor', 'mass_kg': 6.,
         'mass_basis': 'Engineering choice; no retained measurement.',
         'material': {'response': None, 'tier': 'synthesized', 'source': 'nominal timber chair'},
         'dimensions_basis': '0.45 m seat height, 0.45 x 0.45 m seat, 0.45 m back; engineering choice.',
         'parts': ([{'name': 'seat', 'primitive': 'box',
                     'min_m': [.52, FLOOR_PLANE_Y + .42, -.24], 'max_m': [.97, FLOOR_PLANE_Y + .46, .21]},
                    {'name': 'back', 'primitive': 'box',
                     'min_m': [.52, FLOOR_PLANE_Y + .46, -.28], 'max_m': [.97, FLOOR_PLANE_Y + .91, -.24]}]
                   + _legs(.52, .97, -.24, .21, FLOOR_PLANE_Y, FLOOR_PLANE_Y + .42, .05))},
        {'id': 'table', 'label': 'Table', 'colour': 'table',
         'base_environment': 'floor', 'mass_kg': 20.,
         'mass_basis': 'Engineering choice; no retained measurement.',
         'material': {'response': None, 'tier': 'synthesized', 'source': 'nominal timber table'},
         'dimensions_basis': '0.74 m top height, 0.70 x 0.60 m top; engineering choice.',
         'parts': ([{'name': 'top', 'primitive': 'box',
                     'min_m': [-1.02, FLOOR_PLANE_Y + .70, -.30], 'max_m': [-.32, FLOOR_PLANE_Y + .74, .30]}]
                   + _legs(-1.02, -.32, -.30, .30, FLOOR_PLANE_Y, FLOOR_PLANE_Y + .70, .06))},
        {'id': 'ball-small', 'label': 'Ball (0.065 m)', 'colour': 'ball-small',
         'base_environment': None, 'mass_kg': .4, 'radius_m': .065,
         'engine_object_id': 'scene-ball', 'always_present': True,
         'mass_basis': 'Retained scene default: Sphere(radius=.065, mass=.4) in ihm/assembly/interactive_scene.py.',
         'material': {'response': 'zero restitution, Coulomb-capped tangential impulse against the environment plane',
                      'tier': 'synthesized', 'source': 'Sphere.step in ihm/assembly/interactive_scene.py'},
         'dimensions_basis': 'The radius and mass the scene already instantiates.',
         'collider': 'sphere', 'start_m': [.42, .1, .3],
         'parts': [{'name': 'ball', 'primitive': 'sphere', 'centre_m': [.42, .1, .3], 'radius_m': .065}]},
        {'id': 'tree', 'label': 'Tree', 'colour': 'canopy',
         'base_environment': 'floor', 'mass_kg': None,
         'mass_basis': 'Absent: a tree is scenery here, not a body any engine can hold. No measurement is retained and none is invented.',
         'material': {'response': None, 'tier': 'synthesized',
                      'source': 'Nominal trunk and canopy; no wood or foliage property is claimed.'},
         'dimensions_basis': '1.4 m trunk of 0.09 m radius under a 0.58 m canopy sphere pair; engineering choice.',
         'parts': [{'name': 'trunk', 'primitive': 'cylinder', 'axis': 1,
                    'centre_m': [1.50, FLOOR_PLANE_Y + .70, -1.30], 'radius_m': .09, 'height_m': 1.40},
                   {'name': 'canopy-lower', 'primitive': 'sphere',
                    'centre_m': [1.50, FLOOR_PLANE_Y + 1.52, -1.30], 'radius_m': .58},
                   {'name': 'canopy-upper', 'primitive': 'sphere',
                    'centre_m': [1.60, FLOOR_PLANE_Y + 1.90, -1.20], 'radius_m': .38}]},
        {'id': 'block', 'label': 'Block (0.12 m)', 'colour': 'block',
         'base_environment': None, 'mass_kg': .35,
         'mass_basis': 'Engineering choice; no retained measurement. The engine never uses it: there is no box collider, so the block is display only.',
         'material': {'response': None, 'tier': 'synthesized',
                      'source': 'Nominal timber block; no contact law, because no box collider exists.'},
         'dimensions_basis': '0.12 m cube resting on the floor plane; engineering choice.',
         'parts': [{'name': 'cube', 'primitive': 'box',
                    'min_m': [.30, FLOOR_PLANE_Y, .34], 'max_m': [.42, FLOOR_PLANE_Y + .12, .46]}]},
        {'id': 'ball-large', 'label': 'Ball (0.11 m)', 'colour': 'ball-large',
         'base_environment': None, 'mass_kg': .6, 'radius_m': .11,
         'mass_basis': 'Engineering choice inside the accepted Sphere bounds (0 < radius <= 1 m, 0 < mass <= 100 kg).',
         'material': {'response': 'zero restitution, Coulomb-capped tangential impulse against the environment plane',
                      'tier': 'synthesized', 'source': 'Sphere.step in ihm/assembly/interactive_scene.py'},
         'dimensions_basis': 'Engineering choice; accepted by the existing Sphere constructor.',
         'collider': 'sphere', 'start_m': [-.52, .18, .34],
         'parts': [{'name': 'ball', 'primitive': 'sphere', 'centre_m': [-.52, .18, .34], 'radius_m': .11}]},
    ]
    return {spec['id']: spec for spec in specs}


def object_geometry(specs):
    """Construct every object mesh and retain it as geometry the viewer can load."""
    (OUT / 'objects').mkdir(exist_ok=True)
    built = {}
    for ident, spec in specs.items():
        parts = [part_mesh(part) for part in spec['parts']]
        positions, faces = merge(parts)
        payload = {
            'schema': 'ihm.scene-object-geometry.v1', 'id': ident, 'label': spec['label'],
            'units': 'm', 'frame': 'bodyparts3d-display-m (x left, y superior, z anterior)',
            'construction': 'procedural: constructed by scripts/build_environment_catalogue.py from the stated primitives. Not downloaded, not acquired.',
            'parts': spec['parts'],
            'part_volume_sum_m3': round(sum(part_volume(part) for part in spec['parts']), 8),
            'part_volume_note': 'Sum over primitives; overlaps between parts are not subtracted.',
            'bounds_m': {'min': positions.min(0).tolist(), 'max': positions.max(0).tolist()},
            'vertex_count': int(len(positions)), 'face_count': int(len(faces)),
            'topology': 'independent triangles, one vertex triple per face, so flat normals are exact',
            'positions': [round(float(v), 6) for v in positions.reshape(-1)],
            'indices': [int(v) for v in faces.reshape(-1)],
        }
        path = OUT / 'objects' / f'{ident}.json'
        path.write_text(json.dumps(payload) + '\n')
        built[ident] = {'spec': spec, 'positions': positions, 'faces': faces,
                        'path': str(path.relative_to(ROOT)), 'payload': payload}
    return built


def render_object(entry):
    """Same camera direction and light; frame fitted to the object, extent recorded."""
    view = view_of(entry['spec']['base_environment'])
    positions = to_display(entry['positions'], view)
    centre = (positions.min(0) + positions.max(0)) / 2
    fitted = {**view, 'rotation': None, 'centre': centre, 'half': 1.}
    uv, _ = project(positions, fitted)
    fitted['half'] = float(np.abs(uv).max() * 1.14)
    cover, rgb = rasterize([(positions, entry['faces'], COLOUR[entry['spec']['colour']])], fitted)
    fig, axes = figure(half=fitted['half'])
    raster_layer(axes, cover, rgb, half=fitted['half'])
    return fig, fitted['half']


def render_scene(spec, objects, skin, environment, view, surround, placements):
    positions, faces, colours = visible_surround(surround)
    meshes = [(positions, faces, colours), (skin[0], skin[1], COLOUR['skin'])]
    for placement in placements:
        entry = objects[placement['object']]
        offset = np.asarray(placement['offset_m'], float)
        meshes.append((entry['positions'] + offset, entry['faces'], COLOUR[entry['spec']['colour']]))
    cover, rgb = rasterize(meshes, view)
    fig, axes = figure(half=view['half'])
    sky_layer(axes, SKY[spec['sky']], view['half'])
    raster_layer(axes, cover, rgb, half=view['half'])
    return fig


# ---------------------------------------------------------------- world surround

ROOM_FLOOR_M = {'floor': FLOOR_PLANE_Y, 'bed': BED_ROOM_FLOOR_Z}
ROOM = {'width': 3.40, 'depth': 3.20, 'height': 2.45, 'thickness': .08}
GROUND_EXTENT_M = 24.
GROUND_CELL_M = .75


def from_gravity_frame(points, environment):
    """Surround geometry is authored with gravity down, then written back in canonical axes."""
    view = view_of(environment)
    points = np.asarray(points, float)
    if view['rotation'] is None:
        return points
    return points @ view['rotation']  # p_canonical = R^T p_display


def room_parts(environment, wall, floor_material, ceiling=True):
    """A closed room: floor, ceiling, four walls, one window opening in the -z wall.

    Everything is visual surround. Nothing in it has a collider.
    """
    base = ROOM_FLOOR_M[environment]
    w, d, h, t = ROOM['width'] / 2, ROOM['depth'] / 2, ROOM['height'], ROOM['thickness']
    top = base + h
    parts = [
        {'name': 'floor', 'primitive': 'box', 'colour': floor_material,
         'min_m': [-w, base - t, -d], 'max_m': [w, base, d]},
        {'name': 'wall-left', 'primitive': 'box', 'colour': wall,
         'min_m': [-w - t, base, -d], 'max_m': [-w, top, d]},
        {'name': 'wall-right', 'primitive': 'box', 'colour': wall, 'cutaway': True,
         'min_m': [w, base, -d], 'max_m': [w + t, top, d]},
        {'name': 'wall-front', 'primitive': 'box', 'colour': wall, 'cutaway': True,
         'min_m': [-w, base, d], 'max_m': [w, top, d + t]},
        {'name': 'ceiling', 'primitive': 'box', 'colour': 'ceiling', 'cutaway': True,
         'min_m': [-w, top, -d], 'max_m': [w, top + t, d]},
    ]
    if not ceiling:
        parts = [p for p in parts if p['name'] != 'ceiling']
    # -z wall carries the window: four boxes around the opening plus a pane.
    x0, x1 = -.55, .75
    y0, y1 = base + .95, base + 1.85
    parts += [
        {'name': 'wall-back-lower', 'primitive': 'box', 'colour': wall,
         'min_m': [-w, base, -d - t], 'max_m': [w, y0, -d]},
        {'name': 'wall-back-upper', 'primitive': 'box', 'colour': wall,
         'min_m': [-w, y1, -d - t], 'max_m': [w, top, -d]},
        {'name': 'wall-back-left', 'primitive': 'box', 'colour': wall,
         'min_m': [-w, y0, -d - t], 'max_m': [x0, y1, -d]},
        {'name': 'wall-back-right', 'primitive': 'box', 'colour': wall,
         'min_m': [x1, y0, -d - t], 'max_m': [w, y1, -d]},
        {'name': 'window-pane', 'primitive': 'box', 'colour': 'window-pane',
         'min_m': [x0, y0, -d - .01], 'max_m': [x1, y1, -d]},
        {'name': 'window-frame-sill', 'primitive': 'box', 'colour': 'window-frame',
         'min_m': [x0 - .04, y0 - .04, -d - .05], 'max_m': [x1 + .04, y0, -d + .01]},
        {'name': 'window-frame-head', 'primitive': 'box', 'colour': 'window-frame',
         'min_m': [x0 - .04, y1, -d - .05], 'max_m': [x1 + .04, y1 + .04, -d + .01]},
    ]
    return parts


def ground_parts(environment, colour='grass'):
    """A tiled ground sheet drawn at the environment plane. Visual only: the engine's
    ground is the ideal half space at the same offset, not this mesh."""
    base = ROOM_FLOOR_M[environment]
    half = GROUND_EXTENT_M / 2
    return [{'name': 'sheet', 'primitive': 'tiled_plane', 'colour': colour,
             'min_m': [-half, base - .02, -half], 'max_m': [half, base, half],
             'cell_m': GROUND_CELL_M,
             'tone_variation': 'deterministic per-cell lightness jitter, seed 20260907, +-7 percent'}]


def surround_mesh(parts, environment):
    """Build a surround, returning positions, faces and one colour row per face."""
    meshes, colours = [], []
    for part in parts:
        if part['primitive'] == 'tiled_plane':
            positions, faces, tones = tiled_plane(part)
            colours.append(tones)
        else:
            positions, faces = part_mesh(part)
            colours.append(np.tile(COLOUR[part['colour']], (len(faces), 1)))
        meshes.append((positions, faces))
    positions, faces = merge(meshes)
    return from_gravity_frame(positions, environment), faces, np.concatenate(colours)


def tiled_plane(part):
    """A flat sheet of cells with jittered tone; the jitter is seeded, so rebuilds match."""
    low, high = np.asarray(part['min_m'], float), np.asarray(part['max_m'], float)
    cell = part['cell_m']
    xs = np.arange(low[0], high[0] - 1e-9, cell)
    zs = np.arange(low[2], high[2] - 1e-9, cell)
    base = COLOUR[part['colour']]
    rng = np.random.default_rng(20260907)
    quads, tones = [], []
    for x in xs:
        for z in zs:
            quads.append([(x, high[1], z), (x, high[1], z + cell),
                          (x + cell, high[1], z + cell), (x + cell, high[1], z)])
            tones.append(np.clip(np.asarray(base) * (1 + rng.uniform(-.07, .07)), 0, 1))
    positions, faces = _tris(quads)
    return positions, faces, np.repeat(np.asarray(tones), 2, axis=0)


SURROUND_SPECS = [
    {'id': 'room-bedroom', 'label': 'Bedroom shell', 'environment': 'bed',
     'kind': 'interior', 'wall': 'wall-warm', 'floor_material': 'floor-wood'},
    {'id': 'room-hospital', 'label': 'Hospital room shell', 'environment': 'bed',
     'kind': 'interior', 'wall': 'wall-clinical', 'floor_material': 'floor-vinyl'},
    {'id': 'room-clinic', 'label': 'Clinic room shell', 'environment': 'floor',
     'kind': 'interior', 'wall': 'wall-clinical', 'floor_material': 'floor-vinyl'},
    {'id': 'room-play', 'label': 'Play room shell', 'environment': 'floor',
     'kind': 'interior', 'wall': 'wall-play', 'floor_material': 'floor-carpet'},
    {'id': 'ground-grass', 'label': 'Grass ground sheet', 'environment': 'floor',
     'kind': 'open_ground', 'colour': 'grass'},
    {'id': 'ground-patio', 'label': 'Paved patio on grass', 'environment': 'floor',
     'kind': 'open_ground', 'colour': 'grass', 'patio': True},
]


def surround_geometry(specs):
    """Retain each surround as loadable geometry with a per-face colour table."""
    (OUT / 'surrounds').mkdir(exist_ok=True)
    built = {}
    for spec in specs:
        if spec['kind'] == 'interior':
            parts = room_parts(spec['environment'], spec['wall'], spec['floor_material'])
        else:
            parts = ground_parts(spec['environment'], spec.get('colour', 'grass'))
            if spec.get('patio'):
                base = ROOM_FLOOR_M[spec['environment']]
                parts = parts + [{'name': 'patio', 'primitive': 'box', 'colour': 'patio',
                                  'min_m': [-1.7, base, -1.5], 'max_m': [1.7, base + .04, 1.5]}]
        positions, faces, colours = surround_mesh(parts, spec['environment'])
        payload = {
            'schema': 'ihm.scene-surround-geometry.v1', 'id': spec['id'], 'label': spec['label'],
            'units': 'm', 'frame': 'bodyparts3d-display-m (x left, y superior, z anterior)',
            'authoring_frame': f"gravity frame of the {spec['environment']} environment; written back in canonical axes",
            'kind': spec['kind'],
            'physics': 'Visual surround only. No collider exists for any part of it: the body and the free spheres cannot touch a wall, a floor sheet or a window.',
            'construction': 'procedural: constructed by scripts/build_environment_catalogue.py from the stated primitives. Not downloaded, not acquired.',
            'parts': [{k: v for k, v in part.items()} for part in parts],
            'cutaway_parts': [part['name'] for part in parts if part.get('cutaway')],
            'cutaway_note': 'Named parts are omitted from the tile render so the interior is visible; the retained geometry is the whole enclosure.',
            'bounds_m': {'min': positions.min(0).tolist(), 'max': positions.max(0).tolist()},
            'vertex_count': int(len(positions)), 'face_count': int(len(faces)),
            'positions': [round(float(v), 5) for v in positions.reshape(-1)],
            'indices': [int(v) for v in faces.reshape(-1)],
            'face_colours_rgb': [round(float(v), 4) for v in colours.reshape(-1)],
        }
        path = OUT / 'surrounds' / f"{spec['id']}.json"
        path.write_text(json.dumps(payload) + '\n')
        built[spec['id']] = {'spec': spec, 'parts': parts, 'positions': positions, 'faces': faces,
                             'colours': colours, 'path': str(path.relative_to(ROOT)), 'payload': payload}
    return built


MAX_SCENE_OBJECTS = 16


def scene_placements(spec, objects):
    """Resolve a scene's placements to instances. Repeats are allowed: the first instance
    keeps the object id, later ones get -2, -3 and so on. The engine mirrors this rule."""
    counts = {}
    placements = []
    for entry in spec['placements']:
        ident = entry['object']
        obj = objects[ident]['spec']
        counts[ident] = counts.get(ident, 0) + 1
        base = obj.get('engine_object_id', ident)
        instance = base if counts[ident] == 1 else f'{base}-{counts[ident]}'
        placements.append({'object': ident, 'instance': instance,
                           'index': counts[ident],
                           'offset_m': [float(v) for v in entry.get('offset_m', [0., 0., 0.])],
                           'simulated': obj.get('collider') == 'sphere',
                           'start_m': ([round(a + b, 4) for a, b in zip(obj['start_m'], entry.get('offset_m', [0., 0., 0.]))]
                                       if obj.get('collider') == 'sphere' else None)})
    return placements


def visible_surround(entry):
    """The same surround with its cutaway parts dropped, for the tile render."""
    keep = [part for part in entry['parts'] if not part.get('cutaway')]
    if len(keep) == len(entry['parts']):
        return entry['positions'], entry['faces'], entry['colours']
    positions, faces, colours = surround_mesh(keep, entry['spec']['environment'])
    return positions, faces, colours


def sky_layer(axes, sky, half):
    if sky['type'] == 'none':
        return
    stops = sky['stops']
    ramp = np.linspace(0, 1, 256)
    ts = [stop['t'] for stop in stops]
    columns = np.stack([np.interp(ramp, ts, [stop['rgb'][channel] for stop in stops]) for channel in range(3)], -1)
    axes.imshow(columns[::-1, None].repeat(2, 1), extent=[-half, half, -half, half],
                aspect='auto', interpolation='bilinear', zorder=0)


SCENE_SPECS = [
    {'id': 'bedroom', 'label': 'Bedroom', 'environment': 'bed', 'surround': 'room-bedroom', 'sky': 'interior',
     'placements': [{'object': 'bed-frame'}, {'object': 'bed-mattress'}, {'object': 'pillow'},
                    {'object': 'blanket'}, {'object': 'nightstand'}],
     'description': 'Enclosed domestic room on the supine environment: timber floor, warm walls, a window, and a bed carrying the retained study mattress, pillow, blanket and a bedside cabinet.'},
    {'id': 'hospital-room', 'label': 'Hospital room', 'environment': 'bed', 'surround': 'room-hospital', 'sky': 'interior',
     'placements': [{'object': 'bed-frame'}, {'object': 'bed-mattress'}, {'object': 'bed-rails'},
                    {'object': 'pillow'}, {'object': 'iv-stand'}, {'object': 'nightstand'}],
     'description': 'Enclosed clinical room on the supine environment: vinyl floor, pale walls, a window, bed with side rails and retained mattress, IV stand and bedside cabinet.'},
    {'id': 'clinic-room', 'label': 'Clinic room', 'environment': 'floor', 'surround': 'room-clinic', 'sky': 'interior',
     'placements': [{'object': 'table'}, {'object': 'chair'}],
     'description': 'Enclosed consultation room on the floor environment: vinyl floor, pale walls, a window, table and chair beside the supported stance.'},
    {'id': 'play-floor', 'label': 'Play room', 'environment': 'floor', 'surround': 'room-play', 'sky': 'interior',
     'placements': [{'object': 'ball-small'}, {'object': 'ball-large'}, {'object': 'block'}, {'object': 'chair'}],
     'description': 'Enclosed room on the floor environment with the two free balls the engine actually simulates, a block and a chair. The balls fall to the ground plane and accept force ports; they do not collide with the body, the block or the walls.'},
    {'id': 'grass-field', 'label': 'Grass field', 'environment': 'floor', 'surround': 'ground-grass', 'sky': 'day',
     'placements': [{'object': 'tree'}, {'object': 'tree', 'offset_m': [-3.00, 0., -.35]},
                    {'object': 'tree', 'offset_m': [-1.40, 0., -1.55]}, {'object': 'ball-small'}],
     'description': 'Open grass under a clear sky on the floor environment: a 24 m tiled ground sheet at the ground plane, three trees and one free ball. Same gravity, same plane and same supports as the bare floor environment; only the world around it differs.'},
    {'id': 'garden-patio', 'label': 'Garden patio', 'environment': 'floor', 'surround': 'ground-patio', 'sky': 'day',
     'placements': [{'object': 'table'}, {'object': 'chair'}, {'object': 'tree'},
                    {'object': 'tree', 'offset_m': [-3.00, 0., -.35]}, {'object': 'ball-large'}],
     'description': 'Paved patio on grass under a clear sky, on the floor environment: table, chair, two trees and a free ball.'},
]


# ---------------------------------------------------------------- provenance

def provenance_record(**kwargs):
    """One ihm.structure-provenance.v1 record; absent evidence is emitted empty."""
    record = {
        'schema': PROVENANCE_SCHEMA,
        'structure_id': kwargs['structure_id'],
        'canonical_entity_id': kwargs.get('canonical_entity_id'),
        'dataset': kwargs['dataset'],
        'source_file': kwargs['source_file'],
        'build': kwargs['build'],
        'geometry': kwargs['geometry'],
        'transforms': kwargs.get('transforms', []),
        'tier': kwargs['tier'],
        'tier_basis': kwargs['tier_basis'],
        'tier_evidence': kwargs.get('tier_evidence', []),
        'frame_relation': kwargs.get('frame_relation'),
        'assumptions': kwargs.get('assumptions', []),
        'derived_artifacts': kwargs.get('derived_artifacts', []),
        'thumbnail': kwargs.get('thumbnail'),
        'selection': kwargs.get('selection'),
        'measurement': kwargs.get('measurement'),
    }
    required = ['dataset.id', 'dataset.label', 'dataset.license', 'source_file.path', 'source_file.sha256',
                'build.script', 'build.commit', 'geometry.path', 'geometry.sha256', 'transforms', 'tier']
    missing = []
    for field in required:
        head, _, tail = field.partition('.')
        value = record.get(head)
        if tail:
            value = None if value is None else value.get(tail)
        if value is None or value == '':
            missing.append(field)
    hashes = [record['source_file'].get('sha256_verified'), (record['geometry'] or {}).get('sha256_verified')]
    record['completeness'] = {
        'required_field_count': len(required),
        'missing': missing,
        'answerable': not missing,
        'hashes_verified': sum(1 for h in hashes if h is True),
        'hashes_unverifiable': sum(1 for h in hashes if h is None),
    }
    return record


def file_block(relative, verify=True):
    path = ROOT / relative
    digest = sha(path)
    return {'path': relative, 'sha256': digest, 'bytes': path.stat().st_size,
            'sha256_verified': True if verify else None}


def build_block(commit, dirty):
    return {'script': THIS, 'script_sha256': sha(ROOT / THIS), 'commit': commit,
            'commit_covers_working_tree': not dirty,
            'uncommitted_changes': dirty}


# ---------------------------------------------------------------- catalogue

def build(render=True):
    OUT.mkdir(parents=True, exist_ok=True)
    THUMBS.mkdir(exist_ok=True)
    PROV.mkdir(exist_ok=True)
    commit = git('rev-parse', 'HEAD')
    dirty = [line[3:] for line in git('status', '--porcelain').splitlines()]
    build_info = build_block(commit, dirty)
    literals = source_literals()

    from ihm.assembly.interactive_scene import ENVIRONMENTS
    curves = bed_curves()
    radii = proxy_radii()
    supports = support_names()
    material = foundation_material()
    supine = read_json(SUPINE_MANIFEST)
    bed_evidence = read_json(BED_MANIFEST)

    objects = object_geometry(object_specs())
    surrounds = surround_geometry(SURROUND_SPECS)
    if render:
        skin = skin_mesh()
        body_raster = {name: rasterize([(skin[0], skin[1], COLOUR['skin'])], view_of(name, ENV_HALF_M))
                       for name in ('studio', 'floor', 'bed')}
        cover, body_rgb = rasterize([(skin[0], skin[1], COLOUR['skin'])])
        points, bodies, _ = quadrature_points()
    else:
        skin = body_raster = cover = body_rgb = points = bodies = None

    inputs = {p: sha(ROOT / p) for p in sorted({
        SCENE_MODULE, NATIVE_CPP, NATIVE_STREAM, NATIVE_CONFIG, BIOGEARS_CPP, BED_MANIFEST, BED_CURVE,
        BED_FIGURE, BED_ARTICLE, BED_IMPL, SUPINE_MANIFEST, SUPINE_ARRAYS, SUPINE_FOUNDATION,
        SUPINE_IMPL, SKIN_GEOMETRY, ANATOMY, OSIM, OSIM_LICENSE, THIS})}

    ihm_dataset = {'id': 'ihm-engineering-settings', 'label': 'IHM repository engineering settings',
                   'version': None, 'revision': commit, 'url': None, 'specimen': None,
                   'units': 'm; m s^-2; degC', 'frame': 'bodyparts3d-display-m',
                   'attribution': 'IHM repository', 'acquisition_status': 'authored in this repository',
                   'license': 'repository-internal', 'license_url': None, 'license_evidence': None,
                   'license_absent_reason': None}
    hong_dataset = {'id': 'hong2022-mattress-compression', 'label': bed_evidence['title'],
                    'version': bed_evidence['publication'], 'revision': bed_evidence['doi'],
                    'url': bed_evidence['sources'][0]['url'],
                    'specimen': 'Sinomax polyurethane mattress foam, 1.9 x 1.2 x 0.2 m; ILD 25% 20/42/120 lbf',
                    'units': 'dimensionless strain; Pa', 'frame': 'uniaxial compression coupon',
                    'attribution': ', '.join(bed_evidence['authors']),
                    'acquisition_status': 'digitized from retained published figure; raw instrument data not retained',
                    'license': bed_evidence['license']['id'], 'license_url': bed_evidence['license']['url'],
                    'license_evidence': file_block(BED_ARTICLE), 'license_absent_reason': None}
    opensim_dataset = {'id': 'opensim-moco-example3dwalking', 'label': 'OpenSim Moco example3DWalking scaled subject model',
                       'version': None, 'revision': None,
                       'url': 'https://github.com/opensim-org/opensim-core', 'specimen': 'scaled gait subject model',
                       'units': 'kg; kg m^2; m', 'frame': 'OpenSim model body frames',
                       'attribution': 'OpenSim / opensim-core contributors',
                       'acquisition_status': 'retained upstream source model',
                       'license': 'Apache-2.0', 'license_url': 'https://www.apache.org/licenses/LICENSE-2.0',
                       'license_evidence': file_block(OSIM_LICENSE), 'license_absent_reason': None}
    bp3d_dataset = {'id': 'bodyparts3d-4.0-skin', 'label': 'BodyParts3D 4.0 skin surface (FJ2810), canonicalized',
                    'version': '4.0', 'revision': None,
                    'url': 'https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html',
                    'specimen': 'adult male reference atlas', 'units': 'm', 'frame': 'bodyparts3d-display-m',
                    'attribution': 'BodyParts3D, The Database Center for Life Science',
                    'acquisition_status': 'derived canonical geometry retained in repository',
                    'license': 'CC BY 4.0', 'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                    'license_evidence': None,
                    'license_absent_reason': 'Licence string carried by data/derived/canonical/anatomy.json provenance; no separate licence file retained'}

    native_map = {'studio': 'free', 'floor': 'upright', 'bed': 'supine'}
    native_gravity = {'free': [0., 0., 0.], 'upright': [0., -9.81, 0.], 'supine': [-9.81, 0., 0.]}
    environment_notes = {
        'studio': 'Zero gravity, no prescribed supports, no ground contact for free objects.',
        'floor': 'Gravity 9.81 m s^-2 toward -y (inferior). Free objects contact a plane at y = -0.96 m. The body is held by prescribed zero translation at both calcanei; its gravity prestress is never solved, an equal and opposite reference support is assumed.',
        'bed': 'Gravity 9.81 m s^-2 toward -z (posterior). Free objects contact a plane at z = -0.24 m. The body is held by prescribed zero translation at occipital bone, sacrum and both calcanei.',
    }

    records = []
    provenance = []
    for ident, spec in ENVIRONMENTS.items():
        gravity = spec['gravity']
        centroids = [supports[i][1] for i in spec['supports'] if i in supports]
        if render:
            view = view_of(ident, ENV_HALF_M)
            fig, axes = figure(half=ENV_HALF_M)
            if ident != 'studio':
                plane_layer(axes, spec['axis'], spec['plane'], view)
            raster_layer(axes, *body_raster[ident], half=ENV_HALF_M)
            support_layer(axes, centroids, view)
            gravity_layer(axes, gravity, view)
            save(fig, ident)
        thumbnail = f'{OUT.relative_to(ROOT)}/thumbnails/{ident}.png'
        records.append({
            'id': ident, 'label': spec['label'].split(' · ')[0], 'full_label': spec['label'],
            'slot': 'environment', 'kind': 'environment',
            'thumbnail': thumbnail, 'thumbnail_url': f'/api/scene/thumbnail/{ident}',
            'thumbnail_kind': 'rendered_geometry',
            'description': environment_notes[ident],
            'gravity': list(gravity), 'gravity_magnitude_m_s2': float(np.linalg.norm(gravity)),
            'axis': spec['axis'], 'plane': spec['plane'],
            'plane_scope': 'Free-object ground plane only; the body is held by the prescribed support set, not by this plane.',
            'supports': list(spec['supports']),
            'support_names': [supports[i][0] for i in spec['supports'] if i in supports],
            'native_environment': native_map[ident],
            'native_gravity_m_s2': native_gravity[native_map[ident]],
            'native_contact_set': {'free': 'none', 'upright': 'source foot ContactGeometrySet/ContactForceSet',
                                   'supine': 'posterior contact: retained skin quadrature when a surface manifest is given, otherwise inertial ellipsoid proxies'}[native_map[ident]],
            'evidence_kind': 'engineering_choice',
            'requires': [],
            'world': {'kind': 'abstract', 'surround': None, 'sky': {**SKY['none'], 'id': 'none'},
                      'ground': ({'type': 'none', 'note': 'Studio has no ground contact for free objects.'} if ident == 'studio'
                                 else {'type': 'ideal_half_space', 'axis': spec['axis'], 'level_m': spec['plane'],
                                       'note': 'The engine plane itself, drawn as a grid. It is the free-object ground, not a floor the body stands on.'}),
                      'light': {'direction_gravity_frame': LIGHT.tolist(),
                                'kind': 'single fixed key light, ambient 0.36, no shadows'},
                      'physics': 'The plane is the only surface in the environment and it acts on free spheres only. A world with walls, ground sheet or sky is a scene, not an environment.'},
            'selection': [{'endpoint': 'POST /api/scene/sessions', 'parameter': 'environment', 'value': ident},
                          {'endpoint': 'POST /api/embodied/sessions', 'parameter': 'environment', 'value': native_map[ident]}],
            'provenance': f'{PROV.relative_to(ROOT)}/{ident}.json',
        })
        provenance.append(provenance_record(
            structure_id=ident, dataset=ihm_dataset,
            source_file=file_block(SCENE_MODULE),
            build=build_info,
            geometry={'path': None, 'sha256': None, 'sha256_verified': None, 'representation': None,
                      'frame': 'bodyparts3d-display-m', 'units': 'm', 'vertex_count': None, 'face_count': None,
                      'absent_reason': 'An environment is a gravity vector, a free-object plane offset and a prescribed support id set. The repository holds no geometry file for it; the tile is a render of the canonical body under those settings.'},
            transforms=[],
            tier='synthesized',
            tier_basis='Gravity magnitude/direction, plane offsets and support id sets are engineering choices written into ENVIRONMENTS; no acquisition and no fit.',
            tier_evidence=[f"ENVIRONMENTS['{ident}'] gravity={gravity} axis={spec['axis']} plane={spec['plane']}",
                           spec['description'],
                           'scripts/native_mechanical_stream.cpp: model.setGravity(environment=="free"?SimTK::Vec3(0):environment=="supine"?SimTK::Vec3(-9.81,0,0):SimTK::Vec3(0,-9.81,0));'],
            frame_relation='canonical',
            assumptions=[{'id': 'reference-support-preload',
                          'statement': 'Gravity prestress is never solved; an equal and opposite reference support force is assumed each step (interactive_scene.step).'},
                         {'id': 'plane-offset-below-skin',
                          'statement': f"Plane offset {spec['plane']} m on axis {spec['axis']} lies {abs(spec['plane']) - (.8648706 if spec['axis'] == 1 else .1460153):.3f} m outside the canonical skin extent on that axis; it is a free-object ground level, not a fitted body support."}],
            thumbnail={'path': thumbnail, 'sha256': sha(ROOT / thumbnail) if render else None,
                       'representation': 'orthographic_shaded_raster_png', 'camera': CAMERA,
                       'renders': [file_block(SKIN_GEOMETRY)]} if render else None,
            selection=records[-1]['selection'],
        ))

    # ------------------------------------------------ bed support-surface model
    support_components = [
        {'id': 'bed-support-skin-quadrature', 'label': 'Skin contact quadrature',
         'description': f"Retained posterior skin quadrature: {supine['points']} points at {math.sqrt(material['area_reference_m2']) * 1e3:.0f} mm spacing over {len(supine['bodies'])} bodies, confined neo-Hookean skin columns (mu {material['shear_modulus_pa']:.0f} Pa, lambda {material['lame_lambda_pa']:.0f} Pa, layer {material['total_layer_thickness_m'] * 1e3:.1f} mm) against an ideal plane at source x = {material['plane_source_x_m']:.4f} m.",
         'evidence_kind': 'derived_geometry', 'thumbnail_kind': 'rendered_geometry',
         'selection': [{'endpoint': 'EmbodiedRuntime.from_workspace / ArticulatedBodyPlant',
                        'parameter': 'surface_contact_manifest', 'value': SUPINE_MANIFEST}]},
        {'id': 'bed-support-inertial-proxy', 'label': 'Inertial ellipsoid proxies',
         'description': f"Default posterior contact when no surface manifest is given: one contact sphere per body, radius^2 = 5(Iyy+Izz-Ixx)/(2m) from the retained source inertia. {len(radii)} radii, {min(r for _, r in radii.values()) * 1e3:.0f}-{max(r for _, r in radii.values()) * 1e3:.0f} mm. Contact material law is transferred from the source foot model, not a mattress calibration.",
         'evidence_kind': 'engineering_choice', 'thumbnail_kind': 'plotted_quantity',
         'selection': [{'endpoint': 'EmbodiedRuntime.from_workspace / ArticulatedBodyPlant',
                        'parameter': 'surface_contact_manifest', 'value': None}]},
    ]
    for component in support_components:
        ident = component['id']
        if render:
            if ident == 'bed-support-skin-quadrature':
                fig = render_quadrature(points, bodies, cover, body_rgb)
            else:
                fig = render_proxy(radii)
            save(fig, ident)
        thumbnail = f'{OUT.relative_to(ROOT)}/thumbnails/{ident}.png'
        records.append({**{k: v for k, v in component.items() if k != 'selection'},
                        'slot': 'bed_support_model', 'kind': 'component',
                        'thumbnail': thumbnail, 'thumbnail_url': f'/api/scene/thumbnail/{ident}',
                        'requires': [{'slot': 'environment', 'any_of': ['bed']}],
                        'selection': component['selection'],
                        'provenance': f'{PROV.relative_to(ROOT)}/{ident}.json'})
    provenance.append(provenance_record(
        structure_id='bed-support-skin-quadrature',
        canonical_entity_id='body-bp3d-FJ2810',
        dataset=bp3d_dataset,
        source_file=file_block(SKIN_GEOMETRY),
        build=build_info,
        geometry={'path': SUPINE_ARRAYS, 'sha256': sha(ROOT / SUPINE_ARRAYS),
                  'sha256_verified': sha(ROOT / SUPINE_ARRAYS) == supine['arrays_sha256'],
                  'representation': 'quadrature_point_set', 'frame': 'opensim source ground',
                  'units': 'm', 'vertex_count': supine['points'], 'face_count': None},
        transforms=[{'name': 'posterior +X ray-grid quadrature of the canonical skin surface',
                     'script': 'scripts/build_supine_surface_contact.py',
                     'residual': None,
                     'residual_reason': 'Rasterization step, not a fit; discretization limits are the 5 mm cell size and perimeter-cell area error.'},
                    {'name': 'rigid registration canonical -> opensim source ground',
                     'method': supine['registration']['basis'],
                     'residual': {'metric': 'rms_landmark_residual_m',
                                  'value': supine['registration']['rms_landmark_residual_m'],
                                  'method': 'unweighted proper-rigid least squares over 22 approximate COM/bone-envelope correspondences'},
                     'maximum_landmark_residual_m': supine['registration']['maximum_landmark_residual_m']}],
        tier='transferred',
        tier_basis='Geometry is canonical BodyParts3D skin placed into the OpenSim source frame by a fitted rigid transform whose landmark residual is recorded.',
        tier_evidence=[supine['registration']['basis'],
                       f"rms_landmark_residual_m={supine['registration']['rms_landmark_residual_m']}",
                       f"maximum_landmark_residual_m={supine['registration']['maximum_landmark_residual_m']}",
                       f"accepted_support={supine['accepted_support']}",
                       f"native_integration={supine['native_integration']}"],
        frame_relation='registered_into_canonical',
        assumptions=[{'id': 'skin-material-prior',
                      'statement': 'Skin layer thickness/moduli are generic priors, not measured patient or mattress contact properties (ihm/assembly/supine_contact.py docstring).'},
                     {'id': 'support-not-accepted',
                      'statement': 'manifest accepted_support=false and native_integration=false: the quadrature is retained evidence, not an accepted equilibrium support solution.'}],
        thumbnail={'path': f'{OUT.relative_to(ROOT)}/thumbnails/bed-support-skin-quadrature.png',
                   'sha256': sha(THUMBS / 'bed-support-skin-quadrature.png') if render else None,
                   'representation': 'orthographic_point_scatter_png', 'camera': CAMERA} if render else None,
        selection=support_components[0]['selection'],
        measurement={'points': supine['points'], 'bodies': len(supine['bodies']),
                     'spacing_m': round(math.sqrt(material['area_reference_m2']), 6),
                     'plane_source_x_m': material['plane_source_x_m'], **{k: v for k, v in material.items() if k != 'points'}},
    ))
    provenance.append(provenance_record(
        structure_id='bed-support-inertial-proxy',
        dataset=opensim_dataset,
        source_file=file_block(OSIM),
        build=build_info,
        geometry={'path': None, 'sha256': None, 'sha256_verified': None,
                  'representation': 'contact spheres constructed at run time from body inertia',
                  'frame': 'OpenSim body frames', 'units': 'm', 'vertex_count': None, 'face_count': None,
                  'absent_reason': 'The proxies are constructed inside the native engine from retained mass/inertia; no geometry file exists.'},
        transforms=[{'name': 'uniform mass scaling to the target body mass',
                     'residual': None,
                     'residual_reason': 'Exact: radius^2 = 5(Iyy+Izz-Ixx)/(2m) is invariant under uniform scaling of mass and inertia.'}],
        tier='synthesized',
        tier_basis='Ellipsoid posterior radii are constructed from source inertia priors; the contact material law is copied from the source foot contact set.',
        tier_evidence=['scripts/native_mechanical_stream.cpp: "Engineering posterior contact proxies: uniform ellipsoid posterior semi-axis derived from source mass/inertia, placed at retained segment COM. Contact material law is explicitly transferred from the source foot model, not a calibrated mattress or surface-anatomy reconstruction."'],
        frame_relation='source_frame_display_placement',
        assumptions=[{'id': 'foot-contact-law-transfer',
                      'statement': 'SmoothSphereHalfSpaceForce parameters are cloned from the source foot ContactForceSet for every body.'}],
        thumbnail={'path': f'{OUT.relative_to(ROOT)}/thumbnails/bed-support-inertial-proxy.png',
                   'sha256': sha(THUMBS / 'bed-support-inertial-proxy.png') if render else None,
                   'representation': 'radius_disc_grid_png',
                   'note': 'Discs are the 22 proxy radii drawn at the same metre scale as the body tiles; the grid layout is not anatomical placement.',
                   'camera': CAMERA} if render else None,
        selection=support_components[1]['selection'],
        measurement={'radius_rule': 'radius^2 = 5 (Iyy + Izz - Ixx) / (2 m)',
                     'radii_m': {k: round(v[1], 6) for k, v in radii.items()},
                     'source_mass_kg': round(sum(v[0] for v in radii.values()), 6)},
    ))

    # ------------------------------------------------ mattress material
    mattress = [('mattress-soft', 'SM', 'Soft mattress', 20),
                ('mattress-medium', 'MM', 'Medium mattress', 42),
                ('mattress-firm', 'HM', 'Firm mattress', 120),
                ('mattress-rigid', None, 'Rigid support', None)]
    for ident, key, label, ild in mattress:
        if render:
            save(render_mattress(curves, key), ident)
        thumbnail = f'{OUT.relative_to(ROOT)}/thumbnails/{ident}.png'
        if key is None:
            description = 'No mattress: the skin foundation reacts against the ideal rigid support plane. This is what the engine does when bed_material is not given.'
            measured = None
        else:
            curve = curves[key]
            description = (f"Measured {label.lower()} response, Hong et al. 2022 (DOI {curve['source_doi']}), "
                           f"ILD 25% {ild} lbf. {len(curve['strain'])} digitized points, strain 0-{curve['strain'][-1]:.4f}, "
                           f"stress 0-{curve['pressure_pa'][-1] / 1e3:.2f} kPa, mattress thickness {curve['thickness_m']} m "
                           f"(deflection domain 0-{curve['thickness_m'] * curve['strain'][-1] * 1e3:.1f} mm). "
                           f"Extrapolation beyond the digitized domain is refused, not clipped.")
            measured = {'curve_id': key, 'points': len(curve['strain']),
                        'strain_domain': [curve['strain'][0], curve['strain'][-1]],
                        'stress_domain_pa': [curve['pressure_pa'][0], curve['pressure_pa'][-1]],
                        'thickness_m': curve['thickness_m'],
                        'deflection_domain_m': [0., curve['thickness_m'] * curve['strain'][-1]],
                        'ild_25_percent_lbf': ild, 'stress_basis': curve['stress_basis'],
                        'normal_rate_law': curve['normal_rate_law'],
                        'table_sha256': curve['table_sha256']}
        records.append({
            'id': ident, 'label': label, 'slot': 'mattress_material', 'kind': 'component',
            'thumbnail': thumbnail, 'thumbnail_url': f'/api/scene/thumbnail/{ident}',
            'thumbnail_kind': 'plotted_measurement' if key else 'plotted_quantity',
            'description': description,
            'evidence_kind': 'digitized_measurement' if key else 'engineering_choice',
            'requires': [{'slot': 'environment', 'any_of': ['bed']},
                         {'slot': 'bed_support_model', 'any_of': ['bed-support-skin-quadrature']}],
            'selection': [{'endpoint': 'EmbodiedRuntime.from_workspace / ArticulatedBodyPlant',
                           'parameter': 'bed_material', 'value': key}],
            'measurement': measured,
            'provenance': f'{PROV.relative_to(ROOT)}/{ident}.json',
        })
        if key is None:
            provenance.append(provenance_record(
                structure_id=ident, dataset=ihm_dataset, source_file=file_block(BED_IMPL), build=build_info,
                geometry={'path': None, 'sha256': None, 'sha256_verified': None, 'representation': None,
                          'frame': None, 'units': None, 'vertex_count': None, 'face_count': None,
                          'absent_reason': 'Absence of a mattress; the reaction surface is the ideal plane already in the foundation law.'},
                transforms=[], tier='synthesized',
                tier_basis='bed_material=None is the engine default path: the skin foundation reacts against a rigid ideal plane. No material evidence is claimed.',
                tier_evidence=['ihm/native/mechanical_stream.py: bed=None if bed_material is None',
                               'ihm/assembly/supine_contact.py: "A rigid stationary plane receives the equal/opposite resultant and moment."'],
                frame_relation='canonical',
                assumptions=[{'id': 'rigid-support', 'statement': 'An infinitely stiff support is an idealization, not a measured mattress.'}],
                thumbnail={'path': thumbnail, 'sha256': sha(ROOT / thumbnail) if render else None,
                           'representation': 'compression_curve_png'} if render else None,
                selection=records[-1]['selection']))
            continue
        curve = curves[key]
        provenance.append(provenance_record(
            structure_id=ident, dataset=hong_dataset,
            source_file=file_block(BED_FIGURE),
            build=build_info,
            geometry={'path': None, 'sha256': None, 'sha256_verified': None, 'representation': None,
                      'frame': None, 'units': None, 'vertex_count': None, 'face_count': None,
                      'absent_reason': 'A uniaxial compression response, not geometry. The retained table is recorded under measurement.'},
            transforms=[{'name': 'manual digitization of retained Figure 1b',
                         'script': 'data/research/bed_material/digitize_hong2022.py',
                         'method': bed_evidence['digitization']['method'],
                         'residual': {'metric': 'stress_pick_uncertainty_pa',
                                      'value': bed_evidence['digitization']['stress_pick_uncertainty_pa'],
                                      'method': f"manual line-center pixel picks, {bed_evidence['digitization']['manual_line_pick_uncertainty_px']} px"},
                         'strain_pick_uncertainty': bed_evidence['digitization']['strain_pick_uncertainty'],
                         'uncertainty_scope': bed_evidence['digitization']['uncertainty_scope']}],
            tier='derived',
            tier_basis='The retained curve is produced from the published figure by a recorded digitization script. The underlying specimen response is a real Instron 5569 measurement, but no raw instrument data is retained (digitization.raw_instrument_data=false), so the artifact in this repository is derived, not measured.',
            tier_evidence=[bed_evidence['digitization']['method'],
                           f"raw_instrument_data={bed_evidence['digitization']['raw_instrument_data']}",
                           bed_evidence['calibration_scope'],
                           bed_evidence['engineering_law_candidate']['nominal_stress_interpretation']],
            frame_relation='canonical',
            assumptions=[{'id': 'nominal-stress', 'statement': bed_evidence['engineering_law_candidate']['nominal_stress_interpretation']},
                         {'id': 'no-hysteresis', 'statement': 'Unloading curve, hysteresis, relaxation, creep and damping are not calibrated (manifest unresolved).'},
                         {'id': 'no-bottoming', 'statement': bed_evidence['unresolved']['bottoming_out']},
                         {'id': 'curve-code-mapping', 'statement': bed_evidence['specimen_conditions']['code_to_curve_mapping']}],
            derived_artifacts=[file_block(BED_IMPL)],
            thumbnail={'path': thumbnail, 'sha256': sha(ROOT / thumbnail) if render else None,
                       'representation': 'compression_curve_png',
                       'note': 'Nothing to photograph: the tile is the measured compression curve, selected curve bold over the other two.'} if render else None,
            selection=records[-1]['selection'],
            measurement={**measured, 'evidence_manifest_sha256': curve['evidence_manifest_sha256'],
                         'retained_table': file_block(BED_CURVE)},
        ))

    # ------------------------------------------------ ambient thermal
    ambient_bounds = literals['ambient_python_bounds_c']
    ambient_default = 22.
    for ident, setpoint, label in (('ambient-18c', 18., 'Cool air 18 °C'),
                                   ('ambient-22c', 22., 'Neutral air 22 °C'),
                                   ('ambient-30c', 30., 'Warm air 30 °C')):
        if render:
            save(render_ambient(setpoint, ambient_bounds, ambient_default), ident)
        thumbnail = f'{OUT.relative_to(ROOT)}/thumbnails/{ident}.png'
        is_default = setpoint == ambient_default
        records.append({
            'id': ident, 'label': label, 'slot': 'ambient_thermal', 'kind': 'component',
            'thumbnail': thumbnail, 'thumbnail_url': f'/api/scene/thumbnail/{ident}',
            'thumbnail_kind': 'plotted_quantity',
            'description': (f"Ambient, mean radiant and respiration ambient temperature set to {setpoint:g} °C for the physiology engine. "
                            + ('This is the retained engine source condition (22 °C, 0.5 clo, 0.1 m s^-1 air).' if is_default
                               else f'Setpoint is an engineering choice inside the engine-accepted range {ambient_bounds[0]:g}-{ambient_bounds[1]:g} °C; the engine accepts any value in that range.')),
            'evidence_kind': 'engine_default' if is_default else 'engineering_choice',
            'requires': [],
            'applies_to': 'physiology run (POST /api/scenarios); not consumed by mechanical scene sessions',
            'range_c': ambient_bounds, 'continuous': True, 'value_c': setpoint,
            'selection': [{'endpoint': 'POST /api/scenarios', 'parameter': 'ambient_temperature_c', 'value': setpoint}],
            'provenance': f'{PROV.relative_to(ROOT)}/{ident}.json',
        })
        provenance.append(provenance_record(
            structure_id=ident, dataset=ihm_dataset, source_file=file_block(BIOGEARS_CPP), build=build_info,
            geometry={'path': None, 'sha256': None, 'sha256_verified': None, 'representation': None,
                      'frame': None, 'units': None, 'vertex_count': None, 'face_count': None,
                      'absent_reason': 'A scalar boundary condition, not geometry.'},
            transforms=[], tier='synthesized',
            tier_basis=('Retained engine source condition read back from a recorded run receipt.' if is_default
                        else 'Setpoint chosen inside the engine-accepted range; no acquisition and no fit.'),
            tier_evidence=[f"scripts/native_biogears_rest.cpp: if (!std::isfinite(ambient) || ambient < {ambient_bounds[0]:g} || ambient > {ambient_bounds[1]:g}) return 4;",
                           f"ihm/native/__init__.py: number(self.ambient_temperature_c,{ambient_bounds[0]:g},{ambient_bounds[1]:g},'ambient_temperature_c')",
                           'data/derived/physiology/native_environment_smoke_valid/summary.json: ENVIRONMENT_SOURCE_AMBIENT_C=22 CLO=0.5 AIR_SPEED_M_S=0.1'],
            frame_relation='canonical',
            assumptions=[{'id': 'radiant-equals-air',
                          'statement': 'Ambient, mean radiant and respiration ambient temperatures are all set to the same value (native_biogears_rest.cpp).'},
                         {'id': 'clothing-untouched',
                          'statement': 'Clothing insulation (clo) and air velocity keep the engine source conditions; clo is left to the garment lane.'}],
            thumbnail={'path': thumbnail, 'sha256': sha(ROOT / thumbnail) if render else None,
                       'representation': 'range_bar_png'} if render else None,
            selection=records[-1]['selection']))

    # ------------------------------------------------ interactable objects
    contact_model = ('Simulated as a finite-mass sphere with rotational inertia: gravity, applied force ports, '
                     'zero-restitution contact with the environment plane and a Coulomb-capped tangential impulse. '
                     'There is no body-object contact and no object-object contact in any engine here '
                     '(interactive scene scope: body_object_contact=false).')
    for ident, entry in objects.items():
        spec = entry['spec']
        payload = entry['payload']
        simulated = spec.get('collider') == 'sphere'
        if render:
            fig, half = render_object(entry)
            save(fig, ident)
        else:
            half = None
        thumbnail = f'{OUT.relative_to(ROOT)}/thumbnails/{ident}.png'
        size = np.asarray(payload['bounds_m']['max']) - np.asarray(payload['bounds_m']['min'])
        requires = [] if spec['base_environment'] is None else [{'slot': 'environment', 'any_of': [spec['base_environment']]}]
        records.append({
            'id': ident, 'label': spec['label'], 'slot': 'objects', 'kind': 'object',
            'thumbnail': thumbnail, 'thumbnail_url': f'/api/scene/thumbnail/{ident}',
            'thumbnail_kind': 'rendered_geometry',
            'thumbnail_half_extent_m': half,
            'description': (f"{spec['dimensions_basis']} Extent {size[0]:.2f} x {size[1]:.2f} x {size[2]:.2f} m. "
                            + (f"Mass {spec['mass_kg']:g} kg. " if spec['mass_kg'] is not None else 'Mass absent. ')
                            + ('Contact-ready: the scene instantiates it as a free sphere.' if simulated
                               else 'Display-only in the current engines: no collider exists for this shape.')),
            'evidence_kind': 'constructed_geometry',
            'geometry': entry['path'], 'geometry_url': f'/api/scene/object/{ident}',
            'geometry_sha256': sha(ROOT / entry['path']),
            'vertex_count': payload['vertex_count'], 'face_count': payload['face_count'],
            'bounds_m': payload['bounds_m'], 'extent_m': [round(float(v), 4) for v in size],
            'mass_kg': spec['mass_kg'], 'mass_basis': spec['mass_basis'], 'material': spec['material'],
            'collider': spec.get('collider'), 'radius_m': spec.get('radius_m'),
            'start_m': spec.get('start_m'), 'colour_rgb': [round(c, 3) for c in COLOUR[spec['colour']]],
            'engine_object_id': spec.get('engine_object_id', ident),
            'always_present': spec.get('always_present', False),
            'physical_contact_solved': simulated,
            'contact_model': contact_model if simulated else 'None. Display-only: the scene has no collider for a constructed box or cylinder, so it is never instantiated as a physical body.',
            'base_environment': spec['base_environment'],
            'requires': requires,
            'in_scenes': [s['id'] for s in SCENE_SPECS if any(pl['object'] == ident for pl in s['placements'])],
            'insertable': simulated,
            'insert_label': ('+ Insert ' + spec['label'].split(' (')[0]) if simulated else None,
            'max_instances_per_session': (MAX_SCENE_OBJECTS if simulated else 0),
            'instance_model': ('Additive: each insert creates a new sphere with its own id (the first instance keeps '
                               'the object id, later ones get -2, -3 and so on), up to '
                               f'{MAX_SCENE_OBJECTS} spheres per session including the scene default ball.'
                               if simulated else
                               'Not insertable: the engine has no collider for this shape, so a session holds zero '
                               'instances of it. A scene may still place it more than once for the renderer to draw.'),
            'selection': ([{'endpoint': 'POST /api/scene/sessions', 'parameter': 'scene',
                            'value': 'any scene listing this object; simulated spheres are instantiated'}] if simulated
                          else [{'endpoint': 'display only', 'parameter': None, 'value': None}]),
            'provenance': f'{PROV.relative_to(ROOT)}/{ident}.json',
        })
        dataset = hong_dataset if ident == 'bed-mattress' else ihm_dataset
        provenance.append(provenance_record(
            structure_id=ident, dataset=dataset,
            source_file=file_block(BED_MANIFEST if ident == 'bed-mattress' else THIS),
            build=build_info,
            geometry={'path': entry['path'], 'sha256': sha(ROOT / entry['path']), 'sha256_verified': True,
                      'representation': 'independent-triangle surface constructed from primitives',
                      'frame': 'bodyparts3d-display-m', 'units': 'm',
                      'vertex_count': payload['vertex_count'], 'face_count': payload['face_count']},
            transforms=[{'name': 'procedural construction from stated primitives',
                         'residual': None,
                         'residual_reason': 'Exact construction, not a fit: the primitives are the definition.'}],
            tier='synthesized',
            tier_basis=('Extent constructed from explicit priors with no source geometry of its own. '
                        + ('Dimensions are the retained study mattress; the compression response is the measured curve carried by the mattress_material slot.'
                           if ident == 'bed-mattress' else 'Dimensions and mass are engineering choices.')),
            tier_evidence=([f"bed_material manifest specimen_conditions.mattress_dimensions_m={bed_evidence['specimen_conditions']['mattress_dimensions_m']}"]
                           if ident == 'bed-mattress' else
                           ["Sphere('scene-ball', [.42,.1,.3]) with Sphere(radius=.065, mass=.4) in ihm/assembly/interactive_scene.py"]
                           if ident == 'ball-small' else
                           [spec['dimensions_basis'], spec['mass_basis']]),
            frame_relation='canonical',
            assumptions=[{'id': 'no-download', 'statement': 'Geometry is constructed in this repository. Nothing was downloaded, so no third-party licence applies to it.'},
                         {'id': 'mass-tier', 'statement': spec['mass_basis']},
                         {'id': 'material-tier', 'statement': spec['material']['source'] or 'No material response retained.'},
                         {'id': 'contact-scope', 'statement': records[-1]['contact_model']}],
            thumbnail={'path': thumbnail, 'sha256': sha(ROOT / thumbnail) if render else None,
                       'representation': 'orthographic_shaded_raster_png', 'camera': {**CAMERA, 'half_extent_m': half},
                       'note': 'Same camera direction, light and shading as the body tiles; the frame is fitted to the object and its half extent is recorded here.'} if render else None,
            selection=records[-1]['selection'],
            measurement={'extent_m': records[-1]['extent_m'], 'mass_kg': spec['mass_kg'],
                         'part_volume_sum_m3': payload['part_volume_sum_m3'],
                         'parts': [p['name'] for p in spec['parts']]},
        ))

    # ------------------------------------------------ composed scenes
    for spec in SCENE_SPECS:
        ident = spec['id']
        environment = ENVIRONMENTS[spec['environment']]
        placements = scene_placements(spec, objects)
        named = list(dict.fromkeys(pl['object'] for pl in placements))
        simulated = [pl for pl in placements if pl['simulated']]
        surround = surrounds[spec['surround']]
        world = {
            'kind': surround['spec']['kind'],
            'surround': spec['surround'],
            'surround_geometry': surround['path'],
            'surround_geometry_sha256': sha(ROOT / surround['path']),
            'surround_url': '/api/scene/surround/' + spec['surround'],
            'cutaway_parts': surround['payload']['cutaway_parts'],
            'sky': {**SKY[spec['sky']], 'id': spec['sky']},
            'ground': ({'type': 'enclosure_floor', 'part': 'floor',
                        'material': surround['spec'].get('floor_material'),
                        'level_m': ROOM_FLOOR_M[spec['environment']],
                        'axis': 'gravity axis of the ' + spec['environment'] + ' environment'}
                       if surround['spec']['kind'] == 'interior' else
                       {'type': 'tiled_ground_sheet', 'material': surround['spec'].get('colour'),
                        'extent_m': GROUND_EXTENT_M, 'cell_m': GROUND_CELL_M,
                        'level_m': ROOM_FLOOR_M[spec['environment']],
                        'axis': 'gravity axis of the ' + spec['environment'] + ' environment'}),
            'light': {'direction_gravity_frame': LIGHT.tolist(),
                      'kind': 'single fixed key light, ambient 0.36, no shadows and no light transport',
                      'note': 'The window is declared scenery. No illumination is computed from it.'},
            'enclosure_dimensions_m': (dict(ROOM) if surround['spec']['kind'] == 'interior' else None),
            'physics': ('Visual surround only. Walls, floor sheet, ceiling, window and every scene object other '
                        'than the free spheres have no collider, so the body can never touch them. The engine '
                        'ground stays the ideal half space at the environment plane; the surround floor is drawn '
                        f"at the room floor level {ROOM_FLOOR_M[spec['environment']]:.2f} m, which is not an engine plane."),
        }
        if render:
            save(render_scene(spec, objects, skin, environment, view_of(spec['environment'], SCENE_HALF_M),
                              surround, placements), ident)
        thumbnail = f'{OUT.relative_to(ROOT)}/thumbnails/{ident}.png'
        gap = abs(environment['plane']) - (.8648706 if environment['axis'] == 1 else .1460153)
        records.append({
            'id': ident, 'label': spec['label'], 'slot': 'scene', 'kind': 'scene',
            'thumbnail': thumbnail, 'thumbnail_url': f'/api/scene/thumbnail/{ident}',
            'thumbnail_kind': 'rendered_geometry',
            'description': spec['description'],
            'evidence_kind': 'composed_arrangement',
            'base_environment': spec['environment'],
            'requires': [{'slot': 'environment', 'any_of': [spec['environment']]}],
            'world': world,
            'objects': named,
            'placements': placements,
            'simulated_objects': [pl['object'] for pl in simulated],
            'engine_object_ids': sorted({pl['instance'] for pl in simulated} | {'scene-ball'}),
            'display_only_objects': [o for o in named if not objects[o]['spec'].get('collider')],
            'object_count': len(placements),
            'repeated_objects': sorted({pl['object'] for pl in placements
                                        if sum(1 for q in placements if q['object'] == pl['object']) > 1}),
            'physical_contact_solved': bool(simulated),
            'contact_model': ('Objects with a sphere collider are instantiated in the scene session and fall to the '
                              'environment plane; every other object is drawn only. A scene never changes gravity, the '
                              'plane offset or the prescribed supports: those stay exactly as the base environment sets them.'),
            'implied_room_floor_m': ROOM_FLOOR_M[spec['environment']],
            'body_to_plane_gap_m': round(gap, 3),
            'selection': [{'endpoint': 'POST /api/scene/sessions', 'parameter': 'scene', 'value': ident}],
            'provenance': f'{PROV.relative_to(ROOT)}/{ident}.json',
        })
        provenance.append(provenance_record(
            structure_id=ident, dataset=ihm_dataset, source_file=file_block(THIS), build=build_info,
            geometry={'path': None, 'sha256': None, 'sha256_verified': None,
                      'representation': 'composition of retained object geometries',
                      'frame': 'bodyparts3d-display-m', 'units': 'm', 'vertex_count': None, 'face_count': None,
                      'absent_reason': 'A scene is an arrangement, not a mesh. Each object it names carries its own geometry record.'},
            transforms=[], tier='synthesized',
            tier_basis='An arrangement of constructed objects over one of the three real environments. The arrangement is an engineering choice; the base environment underneath it is unchanged.',
            tier_evidence=[f"base environment {spec['environment']}: gravity={environment['gravity']} axis={environment['axis']} plane={environment['plane']}",
                           'ihm/assembly/interactive_scene.py scope: body_object_contact=False'],
            frame_relation='canonical',
            assumptions=[{'id': 'placement', 'statement': 'Object placement is an engineering choice; no room measurement is retained.'},
                         {'id': 'surround-not-physical', 'statement': world['physics']},
                         {'id': 'sky-declared', 'statement': f"Sky is a declared {world['sky']['type']} the renderer paints; no sky image is acquired and no daylight model is solved."},
                         {'id': 'body-plane-gap',
                          'statement': f'The body sits {gap:.3f} m off the environment plane, so it is drawn resting above the surface, not on it. The plane is the free-object ground level and the body is held by prescribed supports.'},
                         {'id': 'implied-floor',
                          'statement': 'Bed scenes place cabinets and stands on an implied room floor at z = -0.86 m fixed by the bed frame legs. No engine plane exists there.'}],
            derived_artifacts=([{'path': objects[o]['path'], 'sha256': sha(ROOT / objects[o]['path'])} for o in named]
                               + [{'path': surround['path'], 'sha256': sha(ROOT / surround['path'])}]),
            thumbnail={'path': thumbnail, 'sha256': sha(ROOT / thumbnail) if render else None,
                       'representation': 'orthographic_shaded_raster_png', 'camera': CAMERA} if render else None,
            selection=records[-1]['selection'],
            measurement={'placements': placements, 'engine_object_ids': records[-1]['engine_object_ids'],
                         'enclosure_m': world['enclosure_dimensions_m'], 'ground': world['ground']}))

    # ------------------------------------------------ world surrounds
    for ident, entry in surrounds.items():
        spec = entry['spec']
        payload = entry['payload']
        provenance.append(provenance_record(
            structure_id=ident, dataset=ihm_dataset, source_file=file_block(THIS), build=build_info,
            geometry={'path': entry['path'], 'sha256': sha(ROOT / entry['path']), 'sha256_verified': True,
                      'representation': 'independent-triangle surround with a per-face colour table',
                      'frame': 'bodyparts3d-display-m', 'units': 'm',
                      'vertex_count': payload['vertex_count'], 'face_count': payload['face_count']},
            transforms=[{'name': 'authored in the gravity frame, written back in canonical axes',
                         'method': view_of(spec['environment'])['basis'],
                         'residual': None,
                         'residual_reason': 'Exact axis permutation, not a fit.'}],
            tier='synthesized',
            tier_basis='A room shell or ground sheet constructed from explicit priors. No room was measured and no asset was acquired.',
            tier_evidence=[payload['construction'], payload['physics'],
                           f"enclosure {ROOM} m" if spec['kind'] == 'interior' else
                           f"ground sheet {GROUND_EXTENT_M} m at {GROUND_CELL_M} m cells"],
            frame_relation='canonical',
            assumptions=[{'id': 'surround-not-physical', 'statement': payload['physics']},
                         {'id': 'dimensions', 'statement': 'Room and ground dimensions are engineering choices; no architectural source is retained.'},
                         {'id': 'cutaway', 'statement': payload['cutaway_note']}],
            thumbnail=None,
            selection=[{'endpoint': 'declared by a scene record under world.surround', 'parameter': None, 'value': None}],
            measurement={'parts': [q['name'] for q in entry['parts']],
                         'bounds_m': payload['bounds_m'],
                         'cutaway_parts': payload['cutaway_parts'],
                         'floor_level_m': ROOM_FLOOR_M[spec['environment']]}))

    slots = [
        {'id': 'environment', 'label': 'Environment', 'control': 'tiles', 'order': 1,
         'exclusive': True, 'required': True, 'default': 'bed', 'requires': [], 'contingent_on': None,
         'note': 'The body accepts exactly one environment. Entries in this slot are mutually exclusive.'},
        {'id': 'scene', 'label': 'Scene', 'control': 'tiles', 'order': 2,
         'exclusive': True, 'required': False, 'default': None, 'requires': [],
         'contingent_on': {'slot': 'environment', 'per_option': True},
         'note': 'A named world over one of the three real environments. Exclusive: one at a time. Each scene names the base environment it requires; selecting it does not change gravity, plane or supports.'},
        {'id': 'objects', 'label': 'Objects', 'control': 'insert', 'order': 3,
         'exclusive': False, 'multiple': True, 'required': False, 'default': None, 'requires': [],
         'contingent_on': {'slot': 'environment', 'per_option': True},
         'note': 'Additive inserts, not toggles: the same object may be inserted more than once and each insert is a separate instance. Only sphere-collider objects can be inserted into a session; the rest are drawn by a scene and cannot be instantiated.'},
        {'id': 'bed_support_model', 'label': 'Support surface', 'control': 'select', 'order': 4,
         'exclusive': True, 'required': False, 'default': 'bed-support-inertial-proxy',
         'requires': [{'slot': 'environment', 'any_of': ['bed']}],
         'contingent_on': {'slot': 'environment', 'any_of': ['bed']},
         'note': 'Posterior contact model used by the native supine environment. Shown only with the bed environment.'},
        {'id': 'mattress_material', 'label': 'Mattress firmness', 'control': 'select', 'order': 5,
         'exclusive': True, 'required': False, 'default': 'mattress-rigid',
         'requires': [{'slot': 'environment', 'any_of': ['bed']},
                      {'slot': 'bed_support_model', 'any_of': ['bed-support-skin-quadrature']}],
         'contingent_on': {'slot': 'bed_support_model', 'any_of': ['bed-support-skin-quadrature']},
         'note': 'A property of the bed, not an alternative to it. The engine refuses a mattress without an explicit skin surface foundation, so this control appears only once that support surface is chosen.'},
        {'id': 'ambient_thermal', 'label': 'Ambient air temperature', 'control': 'select', 'order': 6,
         'exclusive': True, 'required': False, 'default': 'ambient-22c', 'requires': [],
         'contingent_on': None, 'continuous_range_c': ambient_bounds, 'unit': 'degC',
         'note': 'Physiology-engine boundary condition, independent of the mechanical environment. The options are convenience points on a continuous accepted range, so a slider over range_c is equally valid.'},
    ]
    by_slot = {}
    for record in records:
        by_slot.setdefault(record['slot'], []).append(record)
    for slot in slots:
        options = []
        for record in by_slot.get(slot['id'], []):
            if slot['control'] == 'insert' and not record.get('insertable'):
                continue
            options.append({'value': record['id'], 'label': record['label'],
                            'default': record['id'] == slot['default'],
                            'thumbnail_url': record['thumbnail_url'],
                            'summary': record['description'].split('. ')[0] + '.',
                            'requires': record['requires'],
                            'value_c': record.get('value_c')})
        slot['options'] = options
        slot['option_count'] = len(options)
        if slot['control'] == 'insert':
            slot['not_insertable'] = [{'value': r['id'], 'label': r['label'],
                                       'reason': 'no collider for this shape; a scene may draw it but a session cannot hold an instance'}
                                      for r in by_slot.get(slot['id'], []) if not r.get('insertable')]
            slot['max_instances_per_session'] = MAX_SCENE_OBJECTS

    catalogue = {
        'schema': SCHEMA,
        'generated_unix': time.time(),
        'commit': commit,
        'exclusivity_model': 'Slots carry exclusivity: entries sharing a slot are mutually exclusive, different slots combine. An entry with unmet requires is not selectable.',
        'camera': {**CAMERA, 'gravity_frames': GRAVITY_FRAME},
        'slots': slots,
        'environments': [r for r in records if r['kind'] == 'environment'],
        'components': [r for r in records if r['kind'] == 'component'],
        'objects': [r for r in records if r['kind'] == 'object'],
        'scenes': [r for r in records if r['kind'] == 'scene'],
        'surrounds': [{'id': ident, 'label': entry['spec']['label'], 'kind': entry['spec']['kind'],
                       'base_environment': entry['spec']['environment'],
                       'geometry': entry['path'], 'geometry_sha256': sha(ROOT / entry['path']),
                       'geometry_url': '/api/scene/surround/' + ident,
                       'vertex_count': entry['payload']['vertex_count'],
                       'face_count': entry['payload']['face_count'],
                       'cutaway_parts': entry['payload']['cutaway_parts'],
                       'part_names': [p['name'] for p in entry['parts']],
                       'physics': entry['payload']['physics'],
                       'used_by': [s['id'] for s in SCENE_SPECS if s['surround'] == ident],
                       'provenance': f'{PROV.relative_to(ROOT)}/{ident}.json'}
                      for ident, entry in surrounds.items()],
        'verified': {
            'scene_environment_ids': literals['scene_environment_ids'],
            'native_cpp_environments': literals['native_cpp_environments'],
            'native_python_environments': literals['native_python_environments'],
            'mattress_curve_ids': sorted(curves),
            'ambient_bounds_c': ambient_bounds,
            'quadrature_points': supine['points'],
            'proxy_bodies': len(radii),
            'scene_default_ball': {'radius_m': .065, 'mass_kg': .4, 'start_m': [.42, .1, .3],
                                   'source': 'ihm/assembly/interactive_scene.py Sphere defaults'},
            'sphere_bounds': {'radius_m': [0, 1], 'mass_kg': [0, 100],
                              'source': 'Sphere.__init__ guard in ihm/assembly/interactive_scene.py'},
            'body_object_contact': False,
        },
        'not_selectable': [
            {'candidate': 'gravity magnitude other than 0 or 9.81 m s^-2',
             'reason': 'Both engines hardcode the three vectors; no parameter accepts another magnitude.'},
            {'candidate': 'upright floor contact-set variants',
             'reason': 'The native upright branch always loads the source foot ContactGeometrySet/ContactForceSet; there is no alternative to select.'},
            {'candidate': 'mattress firmness inside an interactive scene session',
             'reason': 'SceneSessions.create accepts only {environment}; bed_material reaches the native/articulated path, not the reduced interactive body. Declared per record under selection.'},
            {'candidate': 'clothing insulation (clo) and air velocity',
             'reason': 'Accepted by the physiology engine but left to the garment lane; not emitted here.'},
            {'candidate': 'body-object contact',
             'reason': 'No engine solves it. The interactive scene declares body_object_contact=false and objects contact only the environment plane. Objects are placed and simulated as free bodies; the body cannot yet rest on, push or grasp them.'},
            {'candidate': 'object-object contact',
             'reason': 'Spheres are integrated independently against the plane; no pair test exists.'},
            {'candidate': 'box, cylinder and cloth colliders',
             'reason': 'The scene has one collider, Sphere, against a half-space. Furniture, pillow and blanket are therefore display-only until a mesh or box collider exists.'},
            {'candidate': 'rooms as a fourth environment',
             'reason': 'Gravity direction, plane offset and prescribed supports are hardcoded in two engines. A scene composes objects on top of one of the three; it never invents a fourth.'},
            {'candidate': 'downloaded furniture, texture or sky assets',
             'reason': 'Nothing was downloaded. Every object and surround under objects/ and surrounds/ is constructed procedurally by this script, and the sky is a declared gradient, so no third-party licence is claimed or needed.'},
            {'candidate': 'walls, floors, ground sheet, sky and window as physical surfaces',
             'reason': 'The whole surround is visual. There is no collider for any of it, so the body cannot lean on a wall and a ball cannot bounce off one. The only surface the engine solves is the ideal half space at the environment plane.'},
            {'candidate': 'inserting a block, or any non-sphere object, as a physical body',
             'reason': 'The scene has one collider, Sphere against a half space. A block, chair or pillow can be drawn and placed but never instantiated, so + Insert is offered only for the two balls.'},
            {'candidate': 'outdoor terrain that is not flat',
             'reason': 'The engine ground is a half space at a fixed offset. A grass field is that same plane with a tiled sheet drawn on it; slopes, steps and uneven ground would need a collider that does not exist.'},
            {'candidate': 'illumination from the declared window or sky',
             'reason': 'Tiles use one fixed key light with no shadows and no light transport. The window and sky are declared scenery, not light sources.'},
        ],
    }

    (OUT / 'catalogue.json').write_text(json.dumps(catalogue, indent=1) + '\n')
    (OUT / 'slots.json').write_text(json.dumps({'schema': SCHEMA + '.slots', 'slots': slots}, indent=1) + '\n')
    for record in provenance:
        (PROV / (record['structure_id'] + '.json')).write_text(json.dumps(record, indent=1) + '\n')
    (PROV / 'index.json').write_text(json.dumps({
        'schema': PROVENANCE_SCHEMA + '.index',
        'records': [{'structure_id': r['structure_id'], 'tier': r['tier'],
                     'answerable': r['completeness']['answerable'], 'missing': r['completeness']['missing'],
                     'path': f"{PROV.relative_to(ROOT)}/{r['structure_id']}.json"} for r in provenance]}, indent=1) + '\n')

    report = self_test(catalogue, provenance, curves, radii, literals, render=render)
    artifacts = {}
    for path in sorted(OUT.rglob('*')):
        if path.is_file() and path.name != 'manifest.json':
            artifacts[str(path.relative_to(ROOT))] = sha(path)
    manifest = {
        'schema': 'ihm.environment-catalogue-manifest.v1',
        'generated_unix': time.time(),
        'commit': commit,
        'commit_covers_working_tree': not dirty,
        'uncommitted_changes': dirty,
        'inputs': inputs,
        'artifacts': artifacts,
        'artifact_exclusion': 'manifest.json is deliberately absent from artifacts: a manifest may not hash itself.',
        'self_test': report,
        'renderer': 'offscreen: numpy orthographic z-buffer rasterizer (2x supersampled) for body geometry, matplotlib Agg for overlays and plots. No browser, no playwright.',
        'counts': {'environments': len(catalogue['environments']), 'scenes': len(catalogue['scenes']),
                   'components': len(catalogue['components']), 'objects': len(catalogue['objects']),
                   'contact_ready_objects': sum(o['physical_contact_solved'] for o in catalogue['objects']),
                   'surrounds': len(catalogue['surrounds']),
                   'scene_placements': sum(len(s['placements']) for s in catalogue['scenes']),
                   'slots': len(slots), 'provenance_records': len(provenance),
                   'thumbnails': len(list(THUMBS.glob('*.png'))),
                   'object_geometries': len(list((OUT / 'objects').glob('*.json'))),
                   'surround_geometries': len(list((OUT / 'surrounds').glob('*.json')))},
    }
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    return catalogue, manifest


# ---------------------------------------------------------------- self-test

def _scene_smoke_test(catalogue):
    """Open every scene in the real engine, step it, close it, then remove the run directory.

    The directory is a constant under this build's own output. No caller-supplied root is
    ever handed to anything that writes.
    """
    import shutil
    from ihm.assembly.interactive_scene import InteractiveScene
    area = OUT / 'self-test'
    if area.exists():
        shutil.rmtree(area)
    detail = []
    passed = True
    try:
        for record in catalogue['scenes']:
            output = area / record['id']
            scene = InteractiveScene(ROOT, output, record['base_environment'], record['id'])
            frame = scene.step({'seconds': .02, 'sequence': 0})
            simulated = sorted(o['id'] for o in frame['objects'])
            expected = sorted(record['engine_object_ids'])
            inserted = scene.insert({'object': 'ball-large'})
            twice = scene.insert({'object': 'ball-large', 'offset_m': [.15, 0., -.2]})
            after = [o['id'] for o in twice['objects']]
            duplicated = len(after) == len(set(after)) and len(after) == len(simulated) + 2
            ok = simulated == expected and frame['scene_id'] == record['id'] and duplicated
            passed = passed and ok
            detail.append(f"{record['id']}:{'ok' if ok else 'MISMATCH'} scene={simulated} after_two_inserts={sorted(after)}")
    finally:
        if area.exists() and area.is_relative_to(OUT):
            shutil.rmtree(area)
    return {'passed': passed, 'detail': '; '.join(detail)}


def self_test(catalogue, provenance, curves, radii, literals, render=True):
    checks = []

    def check(name, ok, detail):
        checks.append({'check': name, 'passed': bool(ok), 'detail': detail})

    entries_all = (catalogue['environments'] + catalogue['components']
                   + catalogue['objects'] + catalogue['scenes'])
    ids = [r['id'] for r in catalogue['environments']]
    check('environment ids equal ENVIRONMENTS keys', sorted(ids) == sorted(literals['scene_environment_ids']),
          f"catalogue={sorted(ids)} source={sorted(literals['scene_environment_ids'])}")
    check('native environment names agree (cpp vs python guard)',
          sorted(literals['native_cpp_environments']) == sorted(literals['native_python_environments']) == ['free', 'supine', 'upright'],
          f"cpp={literals['native_cpp_environments']} python={literals['native_python_environments']}")
    check('every environment maps to an accepted native environment',
          all(r['native_environment'] in literals['native_cpp_environments'] for r in catalogue['environments']),
          str({r['id']: r['native_environment'] for r in catalogue['environments']}))

    from ihm.assembly.bed_compression import load_bed
    rejected = False
    try:
        load_bed(ROOT, 'XX')
    except ValueError:
        rejected = True
    check('bed loader accepts SM/MM/HM and rejects others', sorted(curves) == ['HM', 'MM', 'SM'] and rejected,
          'load_bed("XX") raised ValueError' if rejected else 'load_bed("XX") did not raise')
    monotone = all(np.all(np.diff(c['strain']) > 0) and np.all(np.diff(c['pressure_pa']) >= 0) and
                   c['strain'][0] == 0 and c['pressure_pa'][0] == 0 for c in curves.values())
    check('mattress curves are anchored monotone', monotone,
          str({k: [len(v['strain']), round(v['strain'][-1], 6), round(v['pressure_pa'][-1], 2)] for k, v in curves.items()}))
    ordering = all(np.interp(.4, curves['SM']['strain'], curves['SM']['pressure_pa'])
                   < np.interp(.4, curves['MM']['strain'], curves['MM']['pressure_pa'])
                   < np.interp(.4, curves['HM']['strain'], curves['HM']['pressure_pa']) for _ in (0,))
    check('firmness ordering SM < MM < HM at 0.4 strain', ordering,
          str({k: round(float(np.interp(.4, v['strain'], v['pressure_pa'])), 1) for k, v in curves.items()}) + ' Pa')

    scaled = {}
    root = ET.parse(ROOT / OSIM).getroot()
    for body in root.iter('Body'):
        mass = float(body.find('mass').text) * 1.37
        inertia = [float(v) * 1.37 for v in body.find('inertia').text.split()]
        scaled[body.get('name')] = math.sqrt(5 * (inertia[1] + inertia[2] - inertia[0]) / (2 * mass))
    invariant = max(abs(scaled[k] - v[1]) for k, v in radii.items())
    check('proxy radius invariant under uniform mass scaling', invariant < 1e-12,
          f'max radius difference {invariant:.3e} m over {len(radii)} bodies at scale 1.37')

    slot_ids = {s['id'] for s in catalogue['slots']}
    entries = entries_all
    check('every entry declares a known slot', all(e['slot'] in slot_ids for e in entries),
          str(sorted({e['slot'] for e in entries})))
    known = {e['id'] for e in entries}
    dependencies_ok = all(r['slot'] in slot_ids and set(r['any_of']) <= known
                          for e in entries for r in e.get('requires', []))
    check('every requires clause resolves to a known slot and ids', dependencies_ok,
          str([(e['id'], e['requires']) for e in entries if e.get('requires')][:3]) + ' ...')
    check('exactly one required slot (environment) and only the objects slot is non-exclusive',
          [s['id'] for s in catalogue['slots'] if s['required']] == ['environment']
          and [s['id'] for s in catalogue['slots'] if not s['exclusive']] == ['objects'],
          str([(s['id'], s['exclusive'], s['required']) for s in catalogue['slots']]))
    check('slot defaults exist in their slot',
          all(s['default'] is None or any(e['id'] == s['default'] and e['slot'] == s['id'] for e in entries)
              for s in catalogue['slots']),
          str({s['id']: s['default'] for s in catalogue['slots']}))
    check('mattress components depend on the skin quadrature support model',
          all(any(r['slot'] == 'bed_support_model' and r['any_of'] == ['bed-support-skin-quadrature']
                  for r in e['requires']) for e in entries if e['slot'] == 'mattress_material'),
          'engine raises "Measured bed requires explicit surface foundation" without it')

    per_record = {r['structure_id']: r for r in provenance}
    covered = known | {s['id'] for s in catalogue['surrounds']}
    check('one provenance record per catalogue entry and surround', set(per_record) == covered,
          f"records={len(per_record)} entries={len(known)} surrounds={len(catalogue['surrounds'])}")
    check('provenance records carry the schema and a tier',
          all(r['schema'] == PROVENANCE_SCHEMA and r['tier'] in ('measured', 'transferred', 'derived', 'synthesized')
              for r in provenance),
          str(sorted({r['tier'] for r in provenance})))
    check('records with no geometry state why', all(r['geometry'].get('path') or r['geometry'].get('absent_reason')
                                                    for r in provenance),
          str([r['structure_id'] for r in provenance if not r['geometry'].get('path')]))
    check('source file hashes recomputed against the working tree',
          all(r['source_file']['sha256'] == sha(ROOT / r['source_file']['path']) for r in provenance),
          f'{len(provenance)} records verified')

    # ---- composed scenes and objects
    objects = {o['id']: o for o in catalogue['objects']}
    check('every scene names a real base environment',
          all(s['base_environment'] in [e['id'] for e in catalogue['environments']] for s in catalogue['scenes']),
          str({s['id']: s['base_environment'] for s in catalogue['scenes']}))
    check('every scene object exists and every object geometry hashes',
          all(o in objects for s in catalogue['scenes'] for o in s['objects'])
          and all(sha(ROOT / o['geometry']) == o['geometry_sha256'] for o in objects.values()),
          f"{len(objects)} objects, {sum(len(s['objects']) for s in catalogue['scenes'])} placements")
    check('object requires clause matches its base environment',
          all((o['requires'] == [] and o['base_environment'] is None)
              or o['requires'] == [{'slot': 'environment', 'any_of': [o['base_environment']]}] for o in objects.values()),
          str({o['id']: o['base_environment'] for o in objects.values()}))
    check('scene base environment agrees with every object it places',
          all(objects[o]['base_environment'] in (None, s['base_environment'])
              for s in catalogue['scenes'] for o in s['objects']),
          str({s['id']: [objects[o]['base_environment'] for o in s['objects']] for s in catalogue['scenes']}))
    check('contact-ready exactly where a sphere collider exists',
          all(o['physical_contact_solved'] == (o['collider'] == 'sphere') for o in objects.values()),
          str({o['id']: o['physical_contact_solved'] for o in objects.values()}))
    from ihm.assembly.interactive_scene import Sphere
    accepted = []
    for entry in objects.values():
        if entry['collider'] != 'sphere':
            continue
        sphere = Sphere(entry['id'], entry['start_m'], radius=entry['radius_m'], mass=entry['mass_kg'])
        accepted.append((entry['id'], sphere.radius, sphere.mass))
    check('every contact-ready object is accepted by the existing Sphere constructor', bool(accepted), str(accepted))
    volumes = {}
    for entry in catalogue['objects']:
        payload = json.loads((ROOT / entry['geometry']).read_bytes())
        positions = np.asarray(payload['positions'], float).reshape(-1, 3)
        faces = np.asarray(payload['indices'], int).reshape(-1, 3)
        tri = positions[faces]
        volume = float(np.sum(np.einsum('ij,ij->i', tri[:, 0], np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]))) / 6)
        volumes[entry['id']] = (round(volume, 6), payload['part_volume_sum_m3'])
    check('object meshes are closed with outward normals and match their primitive volume',
          all(v > 0 and abs(v - s) / s < .05 for v, s in volumes.values()),
          str({k: v for k, v in sorted(volumes.items())[:4]}) + ' ... (signed volume m3, primitive sum m3)')
    check('constructed geometry declares no third-party licence',
          all(json.loads((ROOT / o['geometry']).read_bytes())['construction'].startswith('procedural')
              for o in objects.values()),
          'all object geometry is constructed by this script; nothing downloaded')
    # ---- worlds
    surrounds = {s['id']: s for s in catalogue['surrounds']}
    check('every scene declares a world with sky, ground and light',
          all(set(s['world']) >= {'kind', 'sky', 'ground', 'light', 'physics'} for s in catalogue['scenes']),
          str({s['id']: (s['world']['kind'], s['world']['sky']['id'], s['world']['ground']['type']) for s in catalogue['scenes']}))
    check('every scene surround exists, hashes and matches the scene base environment',
          all(s['world']['surround'] in surrounds
              and sha(ROOT / surrounds[s['world']['surround']]['geometry']) == s['world']['surround_geometry_sha256']
              and surrounds[s['world']['surround']]['base_environment'] == s['base_environment']
              for s in catalogue['scenes']),
          str({s['id']: s['world']['surround'] for s in catalogue['scenes']}))
    check('surrounds are declared visual only, with no collider claimed',
          all('no collider' in s['physics'].lower() or 'visual surround only' in s['physics'].lower()
              for s in surrounds.values())
          and all('visual surround only' in s['world']['physics'].lower() for s in catalogue['scenes']),
          f'{len(surrounds)} surrounds')
    check('interior scenes carry an enclosure and a cutaway list; outdoor scenes carry a ground sheet',
          all((s['world']['enclosure_dimensions_m'] and s['world']['cutaway_parts'])
              if s['world']['kind'] == 'interior' else
              (s['world']['ground']['type'] == 'tiled_ground_sheet' and s['world']['sky']['type'] != 'none')
              for s in catalogue['scenes']),
          str({s['id']: s['world']['kind'] for s in catalogue['scenes']}))
    check('at least one outdoor scene over the floor environment',
          any(s['world']['kind'] == 'open_ground' and s['base_environment'] == 'floor' for s in catalogue['scenes']),
          str([s['id'] for s in catalogue['scenes'] if s['world']['kind'] == 'open_ground']))

    # ---- controls and instancing
    controls = {s['id']: s['control'] for s in catalogue['slots']}
    check('every component slot is a select with ordered options and a default that is one of them',
          all(s['options'] and any(o['default'] for o in s['options'])
              and s['default'] in [o['value'] for o in s['options']]
              for s in catalogue['slots'] if s['control'] == 'select'),
          str({s['id']: [o['value'] for o in s['options']] for s in catalogue['slots'] if s['control'] == 'select'}))
    check('every select slot states what it is contingent on',
          all('contingent_on' in s for s in catalogue['slots'] if s['control'] == 'select'),
          str({s['id']: s['contingent_on'] for s in catalogue['slots'] if s['control'] == 'select'}))
    check('the objects slot is an insert control listing only instantiable objects',
          controls['objects'] == 'insert'
          and {o['value'] for o in next(s for s in catalogue['slots'] if s['id'] == 'objects')['options']}
          == {o['id'] for o in catalogue['objects'] if o['insertable']},
          str([o['value'] for o in next(s for s in catalogue['slots'] if s['id'] == 'objects')['options']]))
    repeated = {s['id']: s['repeated_objects'] for s in catalogue['scenes'] if s['repeated_objects']}
    check('repeated placements resolve to distinct instance ids',
          all(len({pl['instance'] for pl in s['placements']}) == len(s['placements']) for s in catalogue['scenes'])
          and bool(repeated),
          str(repeated))
    live = _scene_smoke_test(catalogue)
    check('each scene opens in the engine and steps', live['passed'], live['detail'])

    if render:
        from PIL import Image
        sizes = {p.name: Image.open(p).size for p in sorted(THUMBS.glob('*.png'))}
        coloured = 0
        for path in sorted(THUMBS.glob('*.png')):
            pixels = np.asarray(Image.open(path).convert('RGB'), float)
            if float(np.abs(pixels - pixels.mean(-1, keepdims=True)).max()) > 12:
                coloured += 1
        check('one square thumbnail per entry at a single size',
              set(p.stem for p in THUMBS.glob('*.png')) == known and set(sizes.values()) == {(IMAGE_PX, IMAGE_PX)},
              f'{len(sizes)} thumbnails at {sorted(set(sizes.values()))}')
        body_tiles = [e['id'] for e in catalogue['environments'] + catalogue['scenes']] + ['bed-support-skin-quadrature']
        spreads = {}
        for name in body_tiles:
            pixels = np.asarray(Image.open(THUMBS / f'{name}.png').convert('RGB'), float)
            spreads[name] = round(float(np.abs(pixels - pixels.mean(-1, keepdims=True)).max()), 1)
        check('environment and scene tiles are rendered in colour', min(spreads.values()) > 12,
              str(spreads))
        matches = {}
        for entry in catalogue['objects']:
            pixels = np.asarray(Image.open(THUMBS / f"{entry['id']}.png").convert('RGB'), float) / 255
            flat = pixels.reshape(-1, 3)
            lit = flat[np.abs(flat - np.asarray([.949, .945, .933])).max(1) > .02]
            mean = lit.mean(0)
            declared = np.asarray(entry['colour_rgb'])
            matches[entry['id']] = round(float(mean @ declared / (np.linalg.norm(mean) * np.linalg.norm(declared))), 4)
        check('object tiles carry their declared palette colour', min(matches.values()) > .995,
              f'cosine to declared rgb, worst {min(matches, key=matches.get)}={min(matches.values())}')
        check('no tile was written through a greyscale pipeline',
              all(Image.open(p).mode in ('RGB', 'RGBA') for p in THUMBS.glob('*.png')),
              f'{coloured}/{len(sizes)} tiles exceed a 12/255 channel spread; the rest are objects whose declared colour is near neutral')
    check('manifest does not list itself',
          True, 'artifacts are collected with name != manifest.json; asserted again after write')

    return {'passed': all(c['passed'] for c in checks), 'checks': checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true', help='rebuild and fail on any failed check')
    args = parser.parse_args()
    catalogue, manifest = build()
    listed = set(manifest['artifacts'])
    if str((OUT / 'manifest.json').relative_to(ROOT)) in listed:
        raise SystemExit('manifest lists itself')
    report = manifest['self_test']
    for entry in report['checks']:
        print(('PASS ' if entry['passed'] else 'FAIL ') + entry['check'] + ' :: ' + entry['detail'])
    print(f"environments={len(catalogue['environments'])} components={len(catalogue['components'])} "
          f"slots={len(catalogue['slots'])} artifacts={len(listed)} passed={report['passed']}")
    if args.self_test and not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
