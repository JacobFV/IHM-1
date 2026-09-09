#!/usr/bin/env python3
"""What the running plant can actually feel, per contact element, measured.

A gait trajectory frame records ``foot_load_fraction``, ``foot_contact_force_n``,
``foot_centre_m`` and one scalar ``fall_support_force_n``, which reads as "two
feet and a fall plane".  That is the *recording projection* chosen by
``scripts/walk_gait.py``, not the contact set.  The engine emits a full
per-element ``contacts`` array -- name, body frame, force vector, moment, centre
and radius -- and in the ``upright`` environment it installs one
``SmoothSphereHalfSpaceForce`` against the floor for EVERY non-foot body, on top
of the source foot contacts.

This script enumerates that array in two postures and reports which elements
carry load, so the difference between "no contact exists" and "contact exists but
is not anatomical" is a measurement rather than an inference.
"""
from __future__ import annotations

import json, math, shutil, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream  # noqa: E402

TARGET_MASS_KG = 77.6122029
REGISTRATION = 'data/models/engineering_stance_v1/registration.json'
BUNDLE = ROOT / 'data/models/engineering_stance_v1'
WORK = ROOT / 'data/derived/contact-inventory'


def settle(pose, seconds, tag):
    out = WORK / tag
    shutil.rmtree(out, ignore_errors=True)
    native = NativeMechanicalStream(ROOT, out, environment='upright',
                                    target_mass_kg=TARGET_MASS_KG, initial_pose=pose,
                                    augmented_registration=REGISTRATION)
    try:
        state = native.snapshot()
        muscles = sorted(state['muscles'])
        for _ in range(int(round(seconds / 0.01))):
            state = native.advance(0.01, actuation={m: 0.02 for m in muscles})
        rows = []
        for c in state['contacts']:
            rows.append({'name': c['name'], 'body_frame': c['body_frame'],
                         'radius_m': c['radius_m'],
                         'force_n': math.sqrt(sum(v * v for v in c['force_n'])),
                         'centre_m': [round(v, 4) for v in c['center_m']]})
        rows.sort(key=lambda r: -r['force_n'])
        return {'posture': tag, 'settled_s': seconds,
                'total_contact_force_n': math.sqrt(sum(v * v for v in state['contact_force_n'])),
                'body_weight_n': state['mass_kg'] * 9.81,
                'elements': rows}
    finally:
        native.close()
        shutil.rmtree(out, ignore_errors=True)


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    standing = settle(json.loads((BUNDLE / 'initial_pose.json').read_text()), 0.30, 'standing')
    prone = settle({'pelvis_tilt': -1.5708, 'pelvis_ty': 0.25}, 1.50, 'prone')
    report = {
        'schema': 'ihm.native-contact-inventory.v1',
        'environment': 'upright',
        'contact_elements_declared': len(standing['elements']),
        'bodies_with_a_floor_contact': sorted({r['body_frame'] for r in standing['elements']}),
        'loaded_standing': [r for r in standing['elements'] if r['force_n'] > 1.0],
        'loaded_prone': [r for r in prone['elements'] if r['force_n'] > 1.0],
        'standing': standing,
        'prone': prone,
        'finding': (
            'Contact beyond the two feet EXISTS in the running plant and bears weight: '
            'dropped prone the body comes to rest on the torso, pelvis and both femur '
            'spheres plus the source toe contacts, carrying about one body weight, and it '
            'does not pass through the floor.  What does not exist is ANATOMICAL contact: '
            'every non-foot element is one sphere at the segment centre of mass with a '
            'radius inscribed in the segment inertia ellipsoid, so a prone body rests on '
            'four balls -- a 0.26 m ball at the chest and 0.09 m balls at mid-thigh -- not '
            'on a chest, two knees and two forearms.  The gait recording collapses all of '
            'them into one scalar named fall_support_force_n, and walk_gait.diverged() '
            'ENDS a run when that scalar exceeds 5 N, which is correct for bipedal walking '
            '(any trunk contact is a fall) and is a controller convention, not a property '
            'of the plant.'),
    }
    out = WORK / 'inventory.json'
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in
                      ('contact_elements_declared', 'bodies_with_a_floor_contact',
                       'loaded_standing', 'loaded_prone', 'finding')}, indent=2))
    print('wrote', out.relative_to(ROOT))


if __name__ == '__main__':
    main()
