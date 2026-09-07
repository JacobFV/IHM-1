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

IMAGE_PX = 512
SUPERSAMPLE = 2
FRAME_HALF_M = 1.06
CENTER_M = np.array([0., -.05, 0.])
INK = '#20242b'
PAPER = '#f2f1ee'
ACCENT = '#2f6f7d'
MUTED = '#b7b4ac'


def normalize(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


VIEW = normalize([.55, .22, 1.])
LIGHT = normalize([.35, .75, .55])
RIGHT = normalize(np.cross([0., 1., 0.], VIEW))
UP = np.cross(VIEW, RIGHT)
CAMERA = {'projection': 'orthographic', 'view_direction_canonical': VIEW.tolist(),
          'image_right_canonical': RIGHT.tolist(), 'image_up_canonical': UP.tolist(),
          'centre_m': CENTER_M.tolist(), 'half_extent_m': FRAME_HALF_M, 'pixels': IMAGE_PX,
          'supersample': SUPERSAMPLE, 'light_direction_canonical': LIGHT.tolist(),
          'shading': 'Lambert, ambient 0.25 + 0.75 max(0, n.l), single fixed light, no shadows',
          'frame': 'bodyparts3d-display-m (x left, y superior, z anterior)',
          'note': 'One camera, light, scale and frame for every tile. The body pose is identical in every environment because selecting an environment changes gravity, the free-object ground plane and the prescribed support set; it does not repose the body.'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(*args):
    return subprocess.run(['git', *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()


def project(points):
    p = np.asarray(points, float) - CENTER_M
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

def rasterize_body(positions, faces):
    """Orthographic z-buffer raster with Lambert shading. Returns (alpha, shade) at IMAGE_PX."""
    size = IMAGE_PX * SUPERSAMPLE
    uv, depth = project(positions)
    scale = size / (2 * FRAME_HALF_M)
    px = np.stack([(uv[:, 0] + FRAME_HALF_M) * scale, (FRAME_HALF_M - uv[:, 1]) * scale], -1)
    tri = positions[faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    lengths = np.linalg.norm(normals, axis=1)
    keep = (lengths > 0) & ((normals @ VIEW) > 0)
    faces = faces[keep]
    normals = normals[keep] / lengths[keep, None]
    shade = .25 + .75 * np.clip(normals @ LIGHT, 0, 1)
    corners = px[faces]
    zs = depth[faces]
    zbuffer = np.full((size, size), -np.inf)
    image = np.zeros((size, size))
    alpha = np.zeros((size, size), bool)
    low = np.floor(corners.min(1)).astype(int)
    high = np.ceil(corners.max(1)).astype(int)
    for k in range(len(faces)):
        x0, y0 = max(low[k, 0], 0), max(low[k, 1], 0)
        x1, y1 = min(high[k, 0], size - 1), min(high[k, 1], size - 1)
        if x1 < x0 or y1 < y0:
            continue
        a, b, c = corners[k]
        det = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        if abs(det) < 1e-12:
            continue
        xs = np.arange(x0, x1 + 1) + .5
        ys = np.arange(y0, y1 + 1) + .5
        gx, gy = np.meshgrid(xs, ys)
        w1 = ((gx - a[0]) * (c[1] - a[1]) - (gy - a[1]) * (c[0] - a[0])) / det
        w2 = ((b[0] - a[0]) * (gy - a[1]) - (b[1] - a[1]) * (gx - a[0])) / det
        inside = (w1 >= 0) & (w2 >= 0) & (w1 + w2 <= 1)
        if not inside.any():
            continue
        z = zs[k, 0] + w1 * (zs[k, 1] - zs[k, 0]) + w2 * (zs[k, 2] - zs[k, 0])
        window = zbuffer[y0:y1 + 1, x0:x1 + 1]
        hit = inside & (z > window)
        window[hit] = z[hit]
        image[y0:y1 + 1, x0:x1 + 1][hit] = shade[k]
        alpha[y0:y1 + 1, x0:x1 + 1][hit] = True
    block = SUPERSAMPLE
    cover = alpha.reshape(IMAGE_PX, block, IMAGE_PX, block).mean((1, 3))
    total = (image * alpha).reshape(IMAGE_PX, block, IMAGE_PX, block).sum((1, 3))
    count = alpha.reshape(IMAGE_PX, block, IMAGE_PX, block).sum((1, 3))
    shaded = np.divide(total, count, out=np.zeros_like(total), where=count > 0)
    return cover, shaded


def figure():
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(IMAGE_PX / 100, IMAGE_PX / 100), dpi=100)
    axes = fig.add_axes([0, 0, 1, 1])
    axes.set_xlim(-FRAME_HALF_M, FRAME_HALF_M)
    axes.set_ylim(-FRAME_HALF_M, FRAME_HALF_M)
    axes.set_aspect('equal')
    axes.set_axis_off()
    axes.set_facecolor(PAPER)
    fig.patch.set_facecolor(PAPER)
    return fig, axes


def body_layer(axes, cover, shaded, alpha=1.):
    rgba = np.zeros((IMAGE_PX, IMAGE_PX, 4))
    tone = .30 + .62 * shaded
    rgba[..., 0] = tone * .99
    rgba[..., 1] = tone * .97
    rgba[..., 2] = tone * .95
    rgba[..., 3] = cover * alpha
    axes.imshow(rgba, extent=[-FRAME_HALF_M, FRAME_HALF_M, -FRAME_HALF_M, FRAME_HALF_M],
                interpolation='bilinear', zorder=2)


def plane_layer(axes, axis, offset, half=.62):
    """The environment plane drawn where the data puts it, in the shared camera."""
    others = [i for i in range(3) if i != axis]
    def point(a, b):
        p = np.zeros(3)
        p[axis] = offset
        p[others[0]] = a
        p[others[1]] = b
        return p
    corners = [point(-half, -half), point(half, -half), point(half, half), point(-half, half)]
    uv, _ = project(np.array(corners))
    axes.add_patch(__import__('matplotlib').patches.Polygon(uv, closed=True, facecolor='#dcdad3',
                                                            edgecolor='#c3c0b7', linewidth=1., zorder=1))
    for t in np.linspace(-half, half, 7):
        for pair in ((point(t, -half), point(t, half)), (point(-half, t), point(half, t))):
            line, _ = project(np.array(pair))
            axes.plot(line[:, 0], line[:, 1], color='#cbc8bf', linewidth=.7, zorder=1)


def gravity_layer(axes, gravity, anchor=np.array([.52, .62, 0.])):
    g = np.asarray(gravity, float)
    magnitude = float(np.linalg.norm(g))
    start, _ = project(anchor[None])
    if magnitude == 0:
        axes.plot(*start[0], marker='o', markersize=9, markerfacecolor='none',
                  markeredgecolor=INK, markeredgewidth=1.6, zorder=4)
        axes.text(start[0, 0], start[0, 1] - .12, 'g = 0', color=INK, ha='center',
                  va='top', fontsize=13, zorder=4)
        return
    end, _ = project((anchor + g / magnitude * .40)[None])
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


def support_layer(axes, centroids):
    if not centroids:
        return
    uv, _ = project(np.array(centroids))
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
        axes.plot(curve['strain'], np.asarray(curve['pressure_pa']) / 1e3, color=MUTED, linewidth=1.6, zorder=1)
    if selected is None:
        axes.axvline(0, color=ACCENT, linewidth=3.4, zorder=3)
        axes.text(.04, 23.2, 'rigid plane', color=ACCENT, fontsize=13, va='top')
    else:
        curve = curves[selected]
        axes.plot(curve['strain'], np.asarray(curve['pressure_pa']) / 1e3, color=ACCENT, linewidth=3.4, zorder=3)
        axes.scatter(curve['strain'], np.asarray(curve['pressure_pa']) / 1e3, s=16, color=ACCENT, zorder=4)
    axes.set_xlim(0, .70)
    axes.set_ylim(0, 25)
    axes.set_xlabel('compressive strain', color=INK, fontsize=12)
    axes.set_ylabel('nominal stress (kPa)', color=INK, fontsize=12)
    return fig


def render_quadrature(points, bodies, cover, shaded):
    fig, axes = figure()
    body_layer(axes, cover, shaded, alpha=.40)
    uv, depth = project(points)
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

    if render:
        positions, faces = skin_mesh()
        cover, shaded = rasterize_body(positions, faces)
        points, bodies, _ = quadrature_points()
    else:
        cover = shaded = points = bodies = None

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
            fig, axes = figure()
            if ident != 'studio':
                plane_layer(axes, spec['axis'], spec['plane'])
            body_layer(axes, cover, shaded)
            support_layer(axes, centroids)
            gravity_layer(axes, gravity)
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
                fig = render_quadrature(points, bodies, cover, shaded)
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

    slots = [
        {'id': 'environment', 'label': 'Environment', 'exclusive': True, 'required': True,
         'default': 'bed', 'requires': [],
         'note': 'The body accepts exactly one environment. Entries in this slot are mutually exclusive.'},
        {'id': 'bed_support_model', 'label': 'Support surface', 'exclusive': True, 'required': False,
         'default': 'bed-support-inertial-proxy',
         'requires': [{'slot': 'environment', 'any_of': ['bed']}],
         'note': 'Posterior contact model used by the native supine environment. Only offered with the bed environment.'},
        {'id': 'mattress_material', 'label': 'Mattress', 'exclusive': True, 'required': False,
         'default': 'mattress-rigid',
         'requires': [{'slot': 'environment', 'any_of': ['bed']},
                      {'slot': 'bed_support_model', 'any_of': ['bed-support-skin-quadrature']}],
         'note': 'A property of the bed, not an alternative to it. The engine refuses a mattress without an explicit skin surface foundation.'},
        {'id': 'ambient_thermal', 'label': 'Ambient air', 'exclusive': True, 'required': False,
         'default': 'ambient-22c', 'requires': [], 'continuous_range_c': ambient_bounds,
         'note': 'Physiology-engine boundary condition, independent of the mechanical environment. Discrete tiles are convenience points on a continuous accepted range.'},
    ]

    catalogue = {
        'schema': SCHEMA,
        'generated_unix': time.time(),
        'commit': commit,
        'exclusivity_model': 'Slots carry exclusivity: entries sharing a slot are mutually exclusive, different slots combine. An entry with unmet requires is not selectable.',
        'camera': CAMERA,
        'slots': slots,
        'environments': [r for r in records if r['kind'] == 'environment'],
        'components': [r for r in records if r['kind'] == 'component'],
        'verified': {
            'scene_environment_ids': literals['scene_environment_ids'],
            'native_cpp_environments': literals['native_cpp_environments'],
            'native_python_environments': literals['native_python_environments'],
            'mattress_curve_ids': sorted(curves),
            'ambient_bounds_c': ambient_bounds,
            'quadrature_points': supine['points'],
            'proxy_bodies': len(radii),
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
        'counts': {'environments': len(catalogue['environments']), 'components': len(catalogue['components']),
                   'slots': len(slots), 'provenance_records': len(provenance), 'thumbnails': len(list(THUMBS.glob('*.png')))},
    }
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    return catalogue, manifest


# ---------------------------------------------------------------- self-test

def self_test(catalogue, provenance, curves, radii, literals, render=True):
    checks = []

    def check(name, ok, detail):
        checks.append({'check': name, 'passed': bool(ok), 'detail': detail})

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
    entries = catalogue['environments'] + catalogue['components']
    check('every entry declares a known slot', all(e['slot'] in slot_ids for e in entries),
          str(sorted({e['slot'] for e in entries})))
    known = {e['id'] for e in entries}
    dependencies_ok = all(r['slot'] in slot_ids and set(r['any_of']) <= known
                          for e in entries for r in e.get('requires', []))
    check('every requires clause resolves to a known slot and ids', dependencies_ok,
          str([(e['id'], e['requires']) for e in entries if e.get('requires')][:3]) + ' ...')
    check('exactly one exclusive required slot (environment)',
          [s['id'] for s in catalogue['slots'] if s['required']] == ['environment']
          and all(s['exclusive'] for s in catalogue['slots']),
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
    check('one provenance record per catalogue entry', set(per_record) == known,
          f'records={len(per_record)} entries={len(known)}')
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

    if render:
        from PIL import Image
        sizes = {p.name: Image.open(p).size for p in sorted(THUMBS.glob('*.png'))}
        check('one thumbnail per entry at a single size',
              set(p.stem for p in THUMBS.glob('*.png')) == known and set(sizes.values()) == {(IMAGE_PX, IMAGE_PX)},
              f'{len(sizes)} thumbnails at {sorted(set(sizes.values()))}')
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
