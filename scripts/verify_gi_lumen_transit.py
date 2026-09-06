#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.assembly.gi_lumen_transit import transit

def main():
    # Snapshot-based synchronous updates prevent same-tick multi-hop transit.
    water=np.array([10.,0.,0.]);solute=np.array([[2.,3.,0.],[0.,0.,0.],[0.,0.,0.]])
    r=transit(water,solute,[10.,10.,10.],provenance='synthetic requested volumes')
    np.testing.assert_array_equal(r['water_ml'],[0,10,0]);np.testing.assert_array_equal(r['solute_mol'],[[0,0,0],[2,3,0],[0,0,0]])
    assert r['fecal_water_ml']==0
    initial=solute.sum(axis=0);totalwater=water.sum();fw=0.;fs=np.zeros(3)
    for _ in range(3):
        r=transit(water,solute,[100]*3,provenance='synthetic emptying upper bound')
        water=r['water_ml'];solute=r['solute_mol'];fw+=r['fecal_water_ml'];fs+=r['fecal_solute_mol']
        assert np.all(water>=0) and np.all(solute>=0)
        np.testing.assert_allclose(solute.sum(axis=0)+fs,initial)
        assert water.sum()+fw==totalwater
    assert fw==10 and water.sum()==0
    # Dry retained solids are not dissolved aqueous payload and do not teleport.
    r=transit([0,1,0],[[4,0],[1,2],[0,0]],[1,.25,0],provenance='dry and partial fixture')
    np.testing.assert_array_equal(r['solute_mol'][0],[4,0])
    np.testing.assert_allclose(r['solute_mol'][1],[.75,1.5])
    np.testing.assert_allclose(r['solute_mol'][2],[.25,.5])
    for bad in ([1,-1,0],[1,float('nan'),0]):
        try:transit([1,1,1],np.zeros((3,2)),bad,provenance='invalid')
        except ValueError:pass
        else:raise AssertionError('invalid request accepted')
    print('PASS: finite water/species conservation, full/tail/zero transit, explicit fecal sink, no same-tick multi-hop, dry retention, invalid input rejection')

if __name__=='__main__':main()
