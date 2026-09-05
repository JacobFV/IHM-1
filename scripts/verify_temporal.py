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
    deficient=fit_predictor(np.column_stack([y[:,0],2*y[:,0],np.ones(500)]),.1)
    assert deficient.diagnostics['data_rank'] == 1
    for args in [(x,0),([1,np.nan],10)]:
        try: spectral_estimate(*args)
        except ValueError: pass
        else: raise AssertionError('invalid input accepted')
    try: finite_laplace([0,0,1],[1,2,3],[1])
    except ValueError: pass
    else: raise AssertionError('duplicate times accepted')
    from ihm.temporal.atlas import _read
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as directory:
        path=Path(directory)/'rounded.csv'
        path.write_text('Time [s],RESP\n100.00,1\n100.01,2\n100.02,3\n100.02,4\n')
        _, clock, _, cadence = _read(path, sampling_rate_hz=125)
        np.testing.assert_allclose(clock, 100+np.arange(4)/125)
        assert cadence == .008
    print('temporal: analytic Laplace, clock origin, PSD power, coherence, resolvent, stable forecast, serialization, no leakage, deficient rank and validation PASS')

if __name__=='__main__': main()
