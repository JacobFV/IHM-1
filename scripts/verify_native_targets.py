"""Analytic unit, target parser and source window binding checks."""
from pathlib import Path
import sys,tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.calibration.native_targets import parse_target,convert,compare_run


def main():
    np.testing.assert_allclose(convert(125,'mL/min','L/day'),180)
    np.testing.assert_allclose(convert(5.6,'L/min','mL/min'),5600)
    np.testing.assert_allclose(convert(37,'degC','K'),310.15)
    np.testing.assert_allclose(convert(1790,'kcal/day','W'),1790*4184/86400)
    assert parse_target('[60,90]')==dict(kind='range',low=60.,high=90.)
    for text in ['7.2, [7.35,7.4]','[12,20], [13,19]','NaN','[9,2]','180 ± 10']:assert parse_target(text) is None
    try:convert(1,'mL','mmHg')
    except ValueError:pass
    else:raise AssertionError('incompatible quantities accepted')
    with tempfile.TemporaryDirectory() as d:
        root=Path(d);path=root/'native.csv';t=np.arange(1,121)
        data=np.column_stack([t,np.where(t<=60,0,125),np.full(120,36),np.full(120,80),np.full(120,3000)])
        np.savetxt(path,data,delimiter=',',header='#Time(s),GlomerularFiltrationRate(mL/min),CoreTemperature(degC),MeanArterialPressure(mmHg),TotalLungVolume(mL)',comments='')
        targets=[dict(name=n,units=u,reference_value=v) for n,u,v in [('GlomerularFiltrationRate','L/day','180'),('CoreTemperature','degC','37'),('MeanArterialPressure','mmHg','[70,105]'),('TotalLungVolume','L','3.2'),('BloodPH',None,'7.4')]]
        run=compare_run(path,targets,root);rows=run['rows']
        assert run['window']['samples']==60 and run['window']['start_s']==61
        assert rows[0]['absolute_difference']==0
        assert rows[1]['signed_difference']==-1 and rows[1]['relative_difference'] is None
        assert rows[2]['signed_distance_to_reference_range']==0
        assert rows[3]['status']=='undersampled_oscillatory_quantity' and rows[4]['status']=='unmapped_quantity'
    print('Native targets: exact bindings, final-minute window, GFR/flow/energy/temperature conversions, range parsing and ambiguous reference rejection PASS')

if __name__=='__main__':main()
