"""Source identity, dimensional audit, mass balance and independent solver convergence."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.native.reproductive import load_source, run_reproductive, build_reproductive


def main():
    root=Path(__file__).resolve().parents[1];source,audit=load_source(root)
    assert audit['states_count']==4 and audit['constants_count']==20 and audit['time_unit']=='day'
    y,c=source.initConsts();np.testing.assert_allclose(y,[467,40,0,150])
    for t in [0.,2.,7.,10.]:
        algebraic=source.computeAlgebraic(c,np.array(y)[:,None],np.array([t]))[:,0]
        rates=np.array(source.computeRates(t,y,c))
        # Internal release cancels when pool mass and circulating mass are combined.
        np.testing.assert_allclose(rates[0]+c[8]*rates[1],algebraic[10]-c[8]*algebraic[0],rtol=1e-12)
        np.testing.assert_allclose(rates[2]+c[15]*rates[3],algebraic[11]-c[15]*algebraic[2],rtol=1e-12)
    a=run_reproductive(root,method='BDF');b=run_reproductive(root,method='DOP853')
    np.testing.assert_allclose(a['state_values'],b['state_values'],rtol=2e-6,atol=2e-5)
    assert np.min(a['state_values'])>=-1e-9
    lh=next(c for c in a['channels'] if c['id']=='LH');assert lh['si_unit']=='kg/m^3' and lh['si_scale']==1e-6
    inhibin=next(c for c in a['channels'] if c['id']=='Ih');assert inhibin['si_unit'] is None
    np.testing.assert_allclose(a['time_s'][-1],864000)
    for bad in [0,31,np.nan]:
        try:run_reproductive(root,duration_days=bad)
        except ValueError:pass
        else:raise AssertionError('invalid source horizon accepted')
    print('reproductive: source hashes, XML state/constant/unit audit, internal release mass cancellation, BDF/DOP853 agreement, positivity and time basis PASS')

if __name__=='__main__':main()
