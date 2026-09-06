#!/usr/bin/env python3
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.assembly.gi_colon_perfusion import fit,response,exchange,rectal_isotonic_observation

def main():
    rectum=rectal_isotonic_observation();assert rectum['water_absorption_ml_min']==0 and rectum['uncertainty'] is None
    f=fit(); protocol=f['protocol']
    low=response(0,flow_ml_hr=850,protocol=protocol);high=response(150,flow_ml_hr=850,protocol=protocol)
    assert low['rates']['water_ml_min']<0 and high['rates']['water_ml_min']>0
    assert high['rates']['k_meq_hr']<0
    assert abs(high['rates']['water_ml_min']-2.8)<1.0
    assert abs(high['rates']['na_meq_hr']-32)<10
    try:response(220,flow_ml_hr=850,protocol=protocol)
    except ValueError:pass
    else:raise AssertionError('hypertonic/outside protocol silently accepted')
    assert response(220,flow_ml_hr=850,protocol=protocol,transfer_prior=True)['transfer_prior']
    for c in [low,high]:
        w=np.array([1.,2.]);n=np.array([[1e-5,1e-5,1e-5],[1e-5,1e-5,1e-5]])
        r=exchange(w,n,60,c)
        assert (r['water_ml']>=0).all() and (r['ion_mol']>=0).all()
        np.testing.assert_allclose(r['water_ml'].sum(),w.sum())
        np.testing.assert_allclose(r['ion_mol'].sum(axis=0),n.sum(axis=0))
        assert r['applied_time_fraction']<1
    r=exchange([0,0],np.zeros((2,3)),60,high);assert r['applied_time_fraction']==0
    assert len(f['leave_one_concentration_out_error'])==6
    path=Path('data/research/gi_colon/conditional_fit.json');path.write_text(json.dumps(f,indent=2)+'\n')
    print('PASS: condition-specific fit, held-concentration errors, absorption/secretion, finite signed donors, charge residual and protocol boundary')
if __name__=='__main__':main()
