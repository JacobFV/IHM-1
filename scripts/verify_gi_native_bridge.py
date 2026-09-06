#!/usr/bin/env python3
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.assembly.gi_native_bridge import prepare, validate_commit
from ihm.assembly.gi_colon_perfusion import response


def main():
    snap=dict(epoch=7,owner_ids=['native:SI','native:Colon','native:Rectum','native:LIV','external:Feces'],
              water_ml=[10,20,0,100,0],ion_mol=[[.002,.001,.002],[.003,.001,.003],[0,0,0],[.02,.01,.02],[0,0,0]],
              phase='after_native_SI_absorption',current_ledger_owner='explicit:unresolved_colon_current')
    condition=response(150,flow_ml_hr=850,protocol='intact_human_colon_isotonic_NaCl_mannitol')
    tx=prepare(snap,[10,20,10],1,condition)
    validate_commit(snap,tx)
    assert not tx['native_commit_ready']  # Complete electrical paths are still missing.
    d=np.array(tx['water_delta_ml']);n=np.array(tx['ion_delta_mol'])
    assert abs(d.sum())<1e-12
    np.testing.assert_allclose(n.sum(axis=0),0,atol=1e-18)
    assert (np.array(snap['water_ml'])+d>=0).all()
    assert (np.array(snap['ion_mol'])+n>=0).all()
    assert d[4]==0 # Initially empty rectum cannot immediately pass incoming colon packet.
    assert tx['colon_uncompensated_charge_mol']!=0
    # Existing intestinal absorption is not applied a second time.
    assert d[0]==-10
    for key,value in [('epoch',8),('phase','before_native_SI_absorption')]:
        bad=dict(snap);bad[key]=value
        try:validate_commit(bad,tx)
        except ValueError:pass
        else:raise AssertionError('stale/wrong-phase snapshot accepted')
    bad=dict(snap);bad['owner_ids']=['same']*5
    try:prepare(bad,[0]*3,0,condition)
    except ValueError:pass
    else:raise AssertionError('aliased stores accepted')
    bad=dict(snap);bad['current_ledger_owner']=''
    try:prepare(bad,[0]*3,1,condition)
    except ValueError:pass
    else:raise AssertionError('unowned current accepted')
    for key,value in [('epoch',True),('epoch',-1),('current_ledger_owner',7),('current_ledger_owner','   ')]:
        bad=dict(snap);bad[key]=value
        try:prepare(bad,[0]*3,0,condition)
        except ValueError:pass
        else:raise AssertionError('invalid epoch/current owner accepted')
    print('PASS: bridge shared donor conservation, snapshot stale guard, distinct owners, no double SI absorption, explicit charge gap')

if __name__=='__main__':main()
