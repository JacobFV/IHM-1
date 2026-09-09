#!/usr/bin/env python3
"""Source-backed reference math checks; not native gait validation."""
from pathlib import Path
import copy,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.gait_tracking_probe import SlowJointReference
r=SlowJointReference()
assert len(r.names)==23 and not any(n.startswith('pelvis_') for n in r.names)
for t in (0.,.9,1.):
    s=r.sample(t)
    assert all(v==0 for v in s['joint_offsets_rad'].values())
    assert all(v==0 for v in s['joint_target_speed_offsets_rad_s'].values())
for t in (1.123,2.123,5.123):
    eps=1e-6;a=r.sample(t-eps);b=r.sample(t+eps);v=r.sample(t)['joint_target_speed_offsets_rad_s']
    assert max(abs((b['joint_offsets_rad'][n]-a['joint_offsets_rad'][n])/(2*eps)-v[n]) for n in r.names)<1e-8
snapshot={'coordinates':{n:{'value':.1,'speed':.2} for n in (*r.names,'pelvis_tilt')},'muscles':{}}
before=copy.deepcopy(snapshot);observed,target=r.controller_observation(snapshot,2.123)
assert snapshot==before and observed['coordinates']['pelvis_tilt']==snapshot['coordinates']['pelvis_tilt']
assert observed['coordinates']!=snapshot['coordinates']
assert r.sample(8.)['recording_phase']==1 and all(v==0 for v in r.sample(8.)['joint_target_speed_offsets_rad_s'].values())
for kwargs in ({'amplitude':.2},{'ramp_s':0},{'duration_s':-1}):
    try:SlowJointReference(**kwargs)
    except ValueError:pass
    else:raise AssertionError('Unbounded reference parameters accepted')
print('PASS: source-backed23 non-root targets, smooth start, matching target speed, unchanged actual snapshot/root, finite endpoint hold, parameter bounds.')
