#!/usr/bin/env python3
"""Native copied-state static activation protocol and exact replay checks."""
from pathlib import Path
import json,sys,tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(tempfile.mkdtemp(prefix='static_activation_verify_',dir=ROOT/'data/research/locomotion_control'))
report={'schema':'ihm.native-static-activation-check.v1','environments':{}}
for environment in ('upright','supine','free'):
    native=NativeMechanicalStream(ROOT,OUT/environment,environment=environment,target_mass_kg=70,augmented_registration='data/derived/mechanics/initial_stance86/registration.json')
    try:
        state=native.snapshot();saved=native.checkpoint();commands={n:.03 for n in state['muscles']}
        expected=native.advance(.01,actuation=commands)
        native.restore(saved)
        coordinates={'pelvis_ty':state['coordinates']['pelvis_ty']['value']}
        old=native.evaluate_static_pose(coordinates)
        updated=native.evaluate_static_pose(coordinates,activations={'soleus_r':.35,'gait2392_ercspn_r':.2})
        assert updated['activation_overrides']['soleus_r']['actual']==.35
        assert updated['activation_overrides']['gait2392_ercspn_r']['actual']==.2
        assert updated['dynamic_residual_identity_norm']<1e-7
        assert native.snapshot()==state | {'kind':'restored'}
        replay=native.advance(.01,actuation=commands)
        assert replay==expected,'Candidate query changed continuing trajectory'
        native.release(saved)
        report['environments'][environment]={'old_coordinate_only_protocol':True,'activation_override_verified':True,'exact_native_trajectory_replay':True,'dynamic_residual_identity_norm':updated['dynamic_residual_identity_norm'],'maximum_penetration_m':updated['maximum_penetration_m']}
    finally:native.close()
(OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'path':str(OUT/'report.json'),**report},indent=2))
