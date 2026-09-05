#!/usr/bin/env python3
"""Primary measurement lookup and actual Coulomb contact regression."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ihm.calibration.textile import textile_friction, measurement_conditions
from ihm.assembly.clothing import PlaneContact


def main():
    conditions = measurement_conditions()
    estimates = {}
    for region, expected in [('chest', (.90,.19,.79,.20)), ('dorsal_forearm', (.31,.07,.25,.05))]:
        estimate = textile_friction(region, fabric='temel_2022_cotton_single_jersey', conditions=conditions)
        assert tuple(estimate[k] for k in ['static_mean','static_sd','dynamic_mean','dynamic_sd']) == expected
        assert estimate['participants_n'] == 10
        contact = PlaneContact((0,1,0), 0, estimate['static_mean'], estimate['dynamic_mean'])
        result = contact.resolve([[0,0,0]], [[.02,-.04,0]], [.05], .001)
        expected_v = 0 if region == 'chest' else .01
        assert np.isclose(result['velocities_m_s'][0,0], expected_v, atol=1e-15)
        assert result['sticking_count'] == (region == 'chest')
        assert np.allclose(result['impulses_ns'] + result['body_reaction_impulses_ns'], 0)
        estimates[region] = {'estimate': estimate, 'tangential_velocity_after_m_s': float(result['velocities_m_s'][0,0]), 'sticking_count':result['sticking_count']}
    for region,fabric,cond in [('genital_skin','temel_2022_cotton_single_jersey',conditions),('chest','polyester',conditions),('chest','temel_2022_cotton_single_jersey',{}),('chest','temel_2022_cotton_single_jersey',{**conditions,'relative_humidity_percent':90}),('chest','temel_2022_cotton_single_jersey',{**conditions,'direction':'inferior'}),('chest','temel_2022_cotton_single_jersey',{**conditions,'wet':True})]:
        try: textile_friction(region, fabric=fabric, conditions=cond)
        except ValueError: pass
        else: raise AssertionError('unsupported request accepted')
    conditions['temperature_c'] = -100
    assert measurement_conditions()['temperature_c'] == 25
    estimates['chest']['estimate']['static_mean'] = -1
    assert textile_friction('chest', fabric='temel_2022_cotton_single_jersey', conditions=measurement_conditions())['static_mean'] == .90
    root=Path(__file__).resolve().parents[1]
    card=json.loads((root/'data/sources/human-skin-textile-friction.json').read_text())
    for artifact in card['artifacts']:
        p=root/artifact['path']
        assert p.exists(), f'acquired source artifact absent: {p}'
        assert hashlib.sha256(p.read_bytes()).hexdigest()==artifact['sha256']
    out=root/'data/derived/audits/textile';out.mkdir(parents=True,exist_ok=True)
    # Re-fetch after the deliberate mutation test.
    estimates['chest']['estimate']=textile_friction('chest',fabric='temel_2022_cotton_single_jersey',conditions=measurement_conditions())
    (out/'verification.json').write_text(json.dumps({'passed':True,'contact_impulse_ratio':.5,'conditions':measurement_conditions(),'results':estimates,'limitation':'Contact transfer demonstration, not garment validation or impact friction measurement.'},indent=2)+'\n')
    print('PASS: two measured cohorts, unsupported requests, defensive copies, source hashes, actual chest stick / dorsal forearm slide')

if __name__=='__main__': main()
