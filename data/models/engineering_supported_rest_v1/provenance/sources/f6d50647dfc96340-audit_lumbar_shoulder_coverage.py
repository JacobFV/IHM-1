#!/usr/bin/env python3
"""Bounded retained-XML/metadata audit. No OpenSim import or mesh decoding."""
import copy
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TARGET = 'data/derived/mechanics/whole_body_arm26_v2/subject_with_arms.osim'
DONOR = 'data/raw/anatomy/opensim-models/source/Models/Gait2392_Simbody/gait2392_thelen2003muscle.osim'
FITTED = 'data/derived/locomotion/moco3d-torque-v1/subject_walk_scaled_FunctionBasedPathSet.xml'
OUT = ROOT / 'data/research/lumbar_shoulder_coverage'
NAMES = [f'{m}_{s}' for s in ('r', 'l') for m in ('ercspn', 'intobl', 'extobl')]

def vec(s):
    return np.array([float(v) for v in s.split()])

def joint(root, name):
    return next(x for x in root.findall('.//JointSet/objects/*') if x.get('name') == name)

def frames(j):
    return {f.findtext('socket_parent').split('/')[-1]: vec(f.findtext('translation'))
            for f in j.findall('./frames/PhysicalOffsetFrame')}

def build():
    hashes = {}
    def read(path, parser):
        raw = (ROOT / path).read_bytes()
        if len(raw) > 32 * 1024 * 1024:
            raise ValueError('metadata size budget exceeded')
        hashes[path] = hashlib.sha256(raw).hexdigest()
        return parser(raw)
    target, donor, fitted = [read(p, ET.fromstring) for p in (TARGET, DONOR, FITTED)]
    anatomy = read('data/derived/canonical/anatomy.json', json.loads)
    mechanics = read('data/derived/canonical/mechanics.json', json.loads)
    muscles = target.findall('.//ForceSet/objects/*')
    muscles = [m for m in muscles if 'Muscle' in m.tag]
    arm = []
    for m in muscles:
        if m.get('name').startswith('arm26_'):
            points = [{'body': p.findtext('socket_parent_frame').split('/')[-1],
                       'station_m': vec(p.findtext('location')).tolist()}
                      for p in m.findall('.//PathPoint')]
            arm.append({'name': m.get('name'), 'points': points,
                        'crosses_shoulder': any(p['body'] == 'torso' for p in points),
                        'wraps': [ET.tostring(w, encoding='unicode') for w in m.findall('.//PathWrap')]})
    fitted_rows = [{'name': p.get('name'),
                    'coordinates': (p.findtext('coordinate_paths') or '').split()}
                   for p in fitted.findall('.//FunctionBasedPath')]
    # Explicit rigid registration about homologous back-joint frames. No scale
    # inferred from body COMs or target anthropometry. Both frames must align.
    dj, tj = joint(donor, 'back'), joint(target, 'back')
    for j in (dj, tj):
        for f in j.findall('./frames/PhysicalOffsetFrame'):
            assert np.allclose(vec(f.findtext('orientation')), 0, atol=1e-14)
    da, ta = frames(dj), frames(tj)
    axes = [vec(a.findtext('axis')) for a in tj.findall('./SpatialTransform/TransformAxis')[:3]]
    assert [a.findtext('coordinates') for a in dj.findall('./SpatialTransform/TransformAxis')[:3]] == [a.findtext('coordinates') for a in tj.findall('./SpatialTransform/TransformAxis')[:3]]
    for a, b in zip(axes, [vec(a.findtext('axis')) for a in dj.findall('./SpatialTransform/TransformAxis')[:3]]):
        assert np.array_equal(a, b)
    fragment = ET.Element('OpenSimDocument', Version=target.get('Version'))
    forces = ET.SubElement(fragment, 'ForceSet', name='gait2392_lumbar_candidate')
    objects = ET.SubElement(forces, 'objects')
    rows = []
    for name in NAMES:
        source = next(m for m in donor.findall('.//ForceSet/objects/*') if m.get('name') == name)
        candidate = copy.deepcopy(source)
        candidate.set('name', 'gait2392_' + name)
        candidate.find('GeometryPath').set('name', 'path')
        path = candidate.findall('.//PathPoint')
        assert len(path) == 2 and not candidate.findall('.//PathWrap')
        stations = []
        for p in path:
            body = p.findtext('socket_parent_frame').split('/')[-1]
            old = vec(p.findtext('location'))
            new = old - da[body] + ta[body]
            p.find('location').text = ' '.join(format(v, '.17g') for v in new)
            stations.append({'body': body, 'donor_station_m': old.tolist(), 'target_station_m': new.tolist()})
        assert [p['body'] for p in stations] == ['pelvis', 'torso']
        a = np.array(stations[0]['target_station_m']) - ta['pelvis']
        b = np.array(stations[1]['target_station_m']) - ta['torso']
        line = a - b
        torque_per_N = np.cross(b, line / np.linalg.norm(line))
        moment_arms = np.array([axis @ torque_per_N for axis in axes])
        params = {n: float(source.findtext(n)) for n in ('max_isometric_force', 'optimal_fiber_length', 'tendon_slack_length', 'pennation_angle_at_optimal', 'max_contraction_velocity', 'activation_time_constant', 'deactivation_time_constant')}
        rows.append({'source_name': name, 'candidate_name': candidate.get('name'),
                     'type': source.tag, 'stations': stations, 'parameters': params,
                     'neutral_path_length_m': float(np.linalg.norm(line)),
                     'neutral_moment_arms_m': moment_arms.tolist(),
                     'neutral_max_isometric_torque_Nm': (moment_arms * params['max_isometric_force']).tolist()})
        objects.append(candidate)
    selected = []
    terms = ('deltoid', 'pectoralis', 'latissimus', 'supraspinatus', 'infraspinatus', 'subscapularis', 'teres ', 'external oblique', 'internal oblique', 'erector', 'multifidus', 'longissimus thoracis', 'iliocostalis lumborum', 'iliocostalis thoracis', 'spinalis thoracis', 'quadratus lumborum', 'rectus abdominis')
    by_entity = {}
    for m in mechanics['muscles']:
        by_entity.setdefault(m.get('canonical_entity_id'), []).append(m)
    for e in anatomy['entities']:
        if e['id'] in by_entity and any(t in e['name'].lower() for t in terms):
            selected.append({'id': e['id'], 'name': e['name'],
                             'reference_geometry': e.get('reference_geometry'),
                             'mechanics_rows': by_entity[e['id']]})
    report = {'schema': 'ihm.lumbar-shoulder-coverage.v1', 'sources': hashes,
              'native_enabled': False, 'source_registration_validated_anatomically': False,
              'native_muscle_count': len(muscles), 'fitted_paths': fitted_rows,
              'arm_paths': arm, 'lumbar_candidate': rows,
              'registration': {'method': 'joint-frame rigid station transfer, unit scale, matched rotation axes',
                               'donor_joint': 'back', 'target_joint': 'back',
                               'donor_frames_m': {k:v.tolist() for k,v in da.items()},
                               'target_frames_m': {k:v.tolist() for k,v in ta.items()},
                               'coordinate_order': ['lumbar_extension', 'lumbar_bending', 'lumbar_rotation'],
                               'neutral_moment_arm_rank': int(np.linalg.matrix_rank(np.array([r['neutral_moment_arms_m'] for r in rows]).T)),
                               'world_transform_change': False, 'mass_added_kg': 0.0},
              'canonical_candidates': selected,
              'limits': ['Donor strengths and dynamics are transferred model priors, not subject calibration.',
                         'XML is a ForceSet insertion fragment; no native loading or runtime activation certified.',
                         'No excitation assignment, reserve removal, body creation, or inertia change.',
                         'Canonical inferred insertion paths are retained evidence, not certified muscle routing.',
                         'Shoulder crossing is topological; wrapping-dependent moment arms and feasible torque cone require native validation.']}
    ET.indent(fragment)
    xml = ET.tostring(fragment, encoding='utf-8', xml_declaration=True)
    report['candidate_xml_sha256'] = hashlib.sha256(xml).hexdigest()
    return report, xml

if __name__ == '__main__':
    report, xml = build()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    (OUT / 'lumbar_candidate_forces.xml').write_bytes(xml)
    print(json.dumps({'native_muscles': report['native_muscle_count'], 'fitted_paths': len(report['fitted_paths']), 'shoulder_crossing_paths': sum(r['crosses_shoulder'] for r in report['arm_paths']), 'lumbar_candidate_count': len(report['lumbar_candidate']), 'canonical_candidates': len(report['canonical_candidates']), 'neutral_lumbar_rank': report['registration']['neutral_moment_arm_rank']}))
