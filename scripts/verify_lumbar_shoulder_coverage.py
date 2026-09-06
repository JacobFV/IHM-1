#!/usr/bin/env python3
"""Small independent geometric work/registration checks; no native process."""
import hashlib
import json
import xml.etree.ElementTree as ET
import numpy as np
from audit_lumbar_shoulder_coverage import build, OUT

r, xml = build()
assert json.loads((OUT / 'audit.json').read_text()) == r
assert (OUT / 'lumbar_candidate_forces.xml').read_bytes() == xml
assert hashlib.sha256(xml).hexdigest() == r['candidate_xml_sha256']
assert r['native_muscle_count'] == 92 and len(r['fitted_paths']) == 80
assert all(p['coordinates'] for p in r['fitted_paths'])
assert not any('lumbar_' in c or 'arm_' in c for p in r['fitted_paths'] for c in p['coordinates'])
assert sum(p['crosses_shoulder'] for p in r['arm_paths']) == 6
assert len(ET.fromstring(xml).findall('.//Thelen2003Muscle')) == 6
axes = np.array([[0.,0.,1.], [1.,0.,0.], [0.,1.,0.]])
for m in r['lumbar_candidate']:
    points = m['stations']
    target_frames = r['registration']['target_frames_m']
    donor_frames = r['registration']['donor_frames_m']
    a, b = [np.array(p['target_station_m']) - target_frames[p['body']] for p in points]
    for p in points:
        np.testing.assert_allclose(np.array(p['target_station_m']) - target_frames[p['body']], np.array(p['donor_station_m']) - donor_frames[p['body']], atol=1e-15)
    arms = []
    # Independent central finite difference of length: generalized tension is -dL/dq.
    for axis in axes:
        def length(q):
            rotated = b*np.cos(q) + np.cross(axis,b)*np.sin(q) + axis*(axis@b)*(1-np.cos(q))
            return np.linalg.norm(a-rotated)
        arms.append(-(length(1e-6)-length(-1e-6))/2e-6)
    np.testing.assert_allclose(arms, m['neutral_moment_arms_m'], atol=5e-11, rtol=0)
    force = (a-b)/np.linalg.norm(a-b)
    np.testing.assert_allclose(np.cross(b,force)+np.cross(a,-force),0,atol=1e-16)
assert r['registration']['neutral_moment_arm_rank'] == 3
assert not r['native_enabled'] and r['registration']['mass_added_kg'] == 0
print('PASS: deterministic retained-source extraction, six shoulder crossings, lumbar registration, finite-difference virtual work, internal torque balance, rank and disabled activation')
