"""Analytic and out-of-sample numerical checks: python scripts/verify_temporal.py."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.temporal import finite_laplace, spectral_estimate, cross_spectrum, resolvent, fit_predictor, ReducedPredictor


def main():
    t = np.linspace(0, 10, 10001)
    s = np.array([0, .2+2j, 1+4j])
    expected = -np.expm1(-(s+.7)*t[-1])/(s+.7)
    np.testing.assert_allclose(finite_laplace(t, np.exp(-.7*t), s), expected, rtol=3e-6)
    # Translation of the absolute source clock must not alter elapsed-time transform.
    np.testing.assert_allclose(finite_laplace(t+91, np.exp(-.7*t), s), expected, rtol=3e-6)
    fs=100.; t=np.arange(4000)/fs; x=3*np.sin(2*np.pi*5*t)
    p=spectral_estimate(x, fs, nperseg=1000)
    assert p['frequency_hz'][np.argmax(p['psd'])] == 5
    np.testing.assert_allclose(np.sum(p['psd'])*(fs/1000), 4.5, rtol=1e-10)
    c=cross_spectrum(x, 2*x, fs, nperseg=1000)
    np.testing.assert_allclose(c['coherence'][50],1,atol=1e-12)
    assert cross_spectrum(x,np.zeros_like(x),fs,nperseg=1000)['coherence'][50] is None
    assert all(v is None for v in cross_spectrum(np.full(1000,.1),np.full(1000,.1),100,nperseg=100)['coherence'])
    np.testing.assert_allclose(resolvent([[-2]],[[3]],[[4]], [0,1j])[:,0,0],12/(2+np.array([0,1j])))
    theta=.13; a=.998*np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    y=np.empty((500,2)); y[0]=[1,0]
    for i in range(1,len(y)): y[i]=a@y[i-1]
    m=fit_predictor(y,.1,rank=2,train_fraction=.7,ridge=1e-10)
    assert max(abs(np.linalg.eigvals(m.transition))) <= 1
    assert np.mean(m.diagnostics['holdout_rmse'])<.025
    restored=ReducedPredictor.from_dict(json.loads(json.dumps(m.to_dict())))
    np.testing.assert_allclose(m.forecast(y[349],150),restored.forecast(y[349],150))
    changed=y.copy(); changed[350:]+=100
    other=fit_predictor(changed,.1,rank=2,train_fraction=.7,ridge=1e-10)
    np.testing.assert_allclose(m.transition,other.transition)
    np.testing.assert_allclose(m.center,other.center)
    # A change of physical unit must not erase a resolved oscillation or change
    # the standardized predictor. These amplitudes remain far above underflow.
    for amplitude in [1e-13, 1e3]:
        scaled=fit_predictor(y*amplitude,.1,rank=2,train_fraction=.7,ridge=1e-10)
        assert scaled.diagnostics['constant_training_variables']==[]
        np.testing.assert_allclose(scaled.forecast(y[349]*amplitude,150)/amplitude,
                                   m.forecast(y[349],150),rtol=1e-8,atol=1e-10)
    deficient=fit_predictor(np.column_stack([y[:,0],2*y[:,0],np.ones(500)]),.1)
    assert deficient.diagnostics['data_rank'] == 1
    for args in [(x,0),([1,np.nan],10)]:
        try: spectral_estimate(*args)
        except ValueError: pass
        else: raise AssertionError('invalid input accepted')
    try: finite_laplace([0,0,1],[1,2,3],[1])
    except ValueError: pass
    else: raise AssertionError('duplicate times accepted')
    from ihm.temporal.atlas import _read, _run
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as directory:
        root=Path(directory);path=root/'units.csv';clock=np.arange(400)/10
        wave=np.sin(2*np.pi*.5*clock)
        np.savetxt(path,np.column_stack([clock,wave,wave*1e-13,np.zeros(400),np.full(400,.1)]),
                   delimiter=',',header='Time(s),Large(V),Small(V),Zero(V),Constant(V)',comments='')
        run=_run(path,root,'analytic_units','analytic_fixture')
        assert [v['constant'] for v in run['variables']]==[False,False,True,True]
        np.testing.assert_allclose([v['dominant_frequency_hz'] for v in run['variables'][:2]],[.5,.5])
        assert [v['dominant_frequency_hz'] for v in run['variables'][2:]]==[None,None]
        assert run['predictor']['model']['diagnostics']['constant_training_variables']==[2,3]
        json.dumps(run,allow_nan=False)
    with TemporaryDirectory() as directory:
        path=Path(directory)/'rounded.csv'
        path.write_text('Time [s],RESP\n100.00,1\n100.01,2\n100.02,3\n100.02,4\n')
        _, clock, _, cadence = _read(path, sampling_rate_hz=125)
        np.testing.assert_allclose(clock, 100+np.arange(4)/125)
        assert cadence == .008
    with TemporaryDirectory() as directory:
        root=Path(directory);path=root/'slow.csv';t=np.arange(3600.)
        values=np.column_stack([t,np.sin(2*np.pi*t/300),np.sin(2*np.pi*t*.8),np.sin(2*np.pi*t*.7)])
        np.savetxt(path,values,delimiter=',',header='Time(s),CoreTemperature(degC),ArterialPressure(mmHg),TotalLungVolume(mL)',comments='')
        run=_run(path,root,'analytic_slow','analytic_fixture',slow=True)
        assert len(run['variables'])==1 and len(run['omitted_undersampled_channels'])==2
        assert run['window_duration_s']==900
        np.testing.assert_allclose(run['variables'][0]['dominant_frequency_hz'],1/300)
        assert run['laplace']['sigma_per_s']==[0.,1/3600,1/600,1/60]
        json.dumps(run,allow_nan=False)
    print('temporal: analytic Laplace, clock origin, PSD power, coherence, resolvent, stable forecast, unit scaling, serialization, no leakage, deficient rank and validation PASS')

if __name__=='__main__': main()
