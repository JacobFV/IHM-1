"""Source CSF equations, mass balance, reference regime and unit-checked pressure port."""
from pathlib import Path
import sys,json,hashlib,tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.native.csf import CSFModel, PressureInput, classify_source_response, run_csf


def main():
    assert classify_source_response(b'<html><title>Just a moment...</title></html>','pdf')=='access_challenge'
    assert classify_source_response(b'not a PDF','pdf')=='wrong_content_type'
    m=CSFModel();d=m.evaluate([9.5,.15],100,0)
    np.testing.assert_allclose(d['Va_ml'],.15*90.5)
    np.testing.assert_allclose(d['Ra_mmHg_s_per_ml'],49100/90.5**2)
    assert abs(d['q_ml_s']-12.5)<.03
    assert abs(d['qf_ml_s']-d['qo_ml_s'])<.0002
    assert abs(d['capillary_balance_ml_s'])<1e-12
    assert abs(d['intracranial_balance_ml_s'])<1e-12
    np.testing.assert_allclose(m.regulated_compliance(0),.15)
    eps=1e-8
    np.testing.assert_allclose((m.regulated_compliance(eps)-m.regulated_compliance(-eps))/(2*eps),-1.5,rtol=1e-6)
    np.testing.assert_allclose(m.regulated_compliance(-1e5),.15+.75/2)
    np.testing.assert_allclose(m.regulated_compliance(1e5),.15-.075/2)
    p=PressureInput.from_samples([0,1,2],[100,90,95],unit='mmHg')
    q=PressureInput.from_samples([0,1,2],np.array([100,90,95])*133.322387415,unit='Pa')
    np.testing.assert_allclose(p.evaluate(.5),[95,-10]);np.testing.assert_allclose(p.evaluate(.5),q.evaluate(.5))
    for t in [-1,3]:
        try:p.evaluate(t)
        except ValueError:pass
        else:raise AssertionError('pressure extrapolation accepted')
    try:PressureInput.from_samples([0,1],[100,90],unit='mL')
    except ValueError:pass
    else:raise AssertionError('wrong pressure unit accepted')
    root=Path(__file__).resolve().parents[1]
    a=run_csf(root,duration_s=300,method='DOP853');b=run_csf(root,duration_s=300,method='Radau')
    np.testing.assert_allclose(a['state_values'],b['state_values'],rtol=2e-7,atol=2e-8)
    assert max(abs(np.array(a['channels'][0]['values'])-9.5))<.2
    assert max(m.linearization()['eigenvalues_real'])<0
    with tempfile.TemporaryDirectory() as directory:
        path=Path(directory)/'map.csv'
        path.write_text('#Time(s),MeanArterialPressure(mmHg)\n10,100\n11,99\n12,98\n')
        port=PressureInput.from_native_map(path)
        assert port.provenance['sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
        assert port.provenance['source_start_s']==10
        np.testing.assert_allclose(port.evaluate(.5),[99.5,-1])
        driven=run_csf(root,pressure_input=port,duration_s=2,samples=21)
        driven2=run_csf(root,pressure_input=port,duration_s=2,samples=21,method='Radau')
        np.testing.assert_allclose(driven['state_values'],driven2['state_values'],rtol=2e-7,atol=2e-8)
        assert driven['validation']['max_intracranial_balance_ml_s']<1e-12
        path.write_text('#Time(s),MeanArterialPressure(mL)\n10,100\n11,99\n')
        try:PressureInput.from_native_map(path)
        except ValueError:pass
        else:raise AssertionError('dimensionally invalid native pressure accepted')
    print('CSF: primary coefficient targets, nonlinear vascular/CSF conservation, sigmoid slope/saturation, SI pressure port and solver agreement PASS')

if __name__=='__main__':main()
