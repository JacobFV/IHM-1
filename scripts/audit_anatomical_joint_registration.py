"""Bounded registration audit: metadata/XML only, no native launch or mesh loads.

AABB distances are lower bounds to surface distance; envelope endpoints are
screening proxies, never measured articular landmarks. CustomJoint offset-frame
separation can be prescribed kinematics and is not a joint constraint residual.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import xml.etree.ElementTree as ET

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.articulated import CanonicalRegistration

MECHANICS = 'data/derived/canonical/mechanics.json'
SNAPSHOT = 'data/derived/native-stream-smoke-n6e0pvqi/supine/smoke.json'
MODEL = 'data/derived/native-stream-smoke-n6e0pvqi/supine/inputs/subject_walk_scaled.osim'
EXECUTION = 'data/derived/native-stream-smoke-n6e0pvqi/supine/execution.json'
MAX_INPUT_BYTES = 16 * 1024 * 1024


def read_bounded(path):
    if path.stat().st_size > MAX_INPUT_BYTES:
        raise ValueError(f'Metadata exceeds 16 MiB ceiling: {path}')
    return path.read_bytes()


def point(transform, local):
    return transform[:3, :3] @ local + transform[:3, 3]


def box_distance(p, group):
    return float(np.linalg.norm(np.maximum(np.maximum(
        group['bounds_min_m'] - p, p - group['bounds_max_m']), 0)))


def audit(root=ROOT):
    paths = (MECHANICS, SNAPSHOT, MODEL, EXECUTION)
    raw = {p: read_bounded(root / p) for p in paths}
    payload = json.loads(raw[MECHANICS])
    retained = json.loads(raw[SNAPSHOT])
    reg = CanonicalRegistration(payload, retained['initial'])
    rows = []
    for joint in ET.fromstring(raw[MODEL]).findall('.//JointSet/objects/*'):
        frames = {f.get('name'): f for f in joint.findall('frames/PhysicalOffsetFrame')}
        selected = [frames[joint.findtext('socket_' + side + '_frame')]
                    for side in ('parent', 'child')]
        bodies = [f.findtext('socket_parent').split('/')[-1] for f in selected]
        if 'ground' in bodies:
            continue
        locals_ = [np.fromstring(f.findtext('translation'), sep=' ') for f in selected]
        row = {'joint': joint.get('name'), 'type': joint.tag, 'bodies': bodies}
        for phase in ('initial', 'step'):
            state = retained[phase]
            source = [point(np.asarray(state['bodies'][b]['transform_ground']), p)
                      for b, p in zip(bodies, locals_)]
            world = [point(reg.global_map, p) for p in source]
            row[phase + '_offset_origin_separation_m'] = float(np.linalg.norm(source[0] - source[1]))
            row[phase + '_common_world_distance_error_m'] = abs(
                float(np.linalg.norm(world[0] - world[1])) - row[phase + '_offset_origin_separation_m'])
            if phase == 'initial':
                row['canonical_offset_origins_m'] = [p.tolist() for p in world]
                row['offset_origin_to_own_bone_aabb_lower_bound_m'] = [
                    box_distance(p, reg.groups[b]) for b, p in zip(bodies, world)]
        rows.append(row)
    endpoints = []
    for side in ('r', 'l'):
        for name, parent, child in [('knee', 'femur', 'tibia'), ('ankle', 'tibia', 'talus'), ('elbow', 'humerus', 'ulna')]:
            a, b = (reg.groups[x + '_' + side] for x in (parent, child))
            pa, pb = (g['anchor_m'].copy() for g in (a, b))
            pa[1], pb[1] = a['bounds_min_m'][1], b['bounds_max_m'][1]
            joint = next(r for r in rows if r['joint'] == ('walker_knee' if name == 'knee' else name) + '_' + side)
            origins = np.asarray(joint['canonical_offset_origins_m'])
            endpoints.append({'joint': name + '_' + side,
                'basis': 'distal/proximal Y-envelope face centers; not articular landmarks',
                'envelope_endpoint_separation_m': float(np.linalg.norm(pa - pb)),
                'signed_y_gap_m': float(pa[1] - pb[1]),
                'offset_origin_to_envelope_endpoint_m': [float(np.linalg.norm(origins[0] - pa)), float(np.linalg.norm(origins[1] - pb))]})
    rigid_rows = [r for r in rows if r['type'] in ('PinJoint', 'WeldJoint')]
    checks = {
        'retained_source_identity': hashlib.sha256(raw[MODEL]).hexdigest() == json.loads(raw[EXECUTION])['source_sha256']['subject_walk_scaled.osim'],
        'one_common_map': all(np.array_equal(m, reg.global_map) for m in reg.maps.values()),
        'proper_rotation': bool(np.allclose(reg.basis.T @ reg.basis, np.eye(3), atol=1e-12) and abs(np.linalg.det(reg.basis) - 1) < 1e-12),
        'common_world_preserves_joint_distances': all(r[p + '_common_world_distance_error_m'] < 1e-12 for r in rows for p in ('initial', 'step')),
        'pin_weld_offset_origins_close': all(r[p + '_offset_origin_separation_m'] < 1e-8 for r in rigid_rows for p in ('initial', 'step')),
    }
    return {'schema': 'ihm.anatomical-registration-audit.v1',
        'scope': 'Retained initial/step metadata only; no surface continuity or anatomical calibration claim',
        'input_sha256': {p: hashlib.sha256(v).hexdigest() for p, v in raw.items()},
        'input_bytes': {p: len(v) for p, v in raw.items()},
        'canonical_registration_source_files': payload['source_files'],
        'global_fit': reg.global_fit, 'joint_offset_frames': rows,
        'envelope_endpoint_proxies': endpoints, 'invariants': checks,
        'mass_ownership': 'Audit adds no bodies or mass; canonical material mass is not independently integrated'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, help='Optional fresh JSON output; defaults to stdout')
    args = parser.parse_args()
    result = audit(args.root)
    serialized = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        with args.output.open('x') as stream:
            stream.write(serialized)
    else:
        print(serialized, end='')
    if not all(result['invariants'].values()):
        raise SystemExit('Registration invariant failure; see audit output')


if __name__ == '__main__':
    main()
