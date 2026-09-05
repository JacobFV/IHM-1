"""Independent analytical and actual-donor causal checks; run from repository root."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy import signal
from types import SimpleNamespace
from ihm.brain.ibm_backend import IBMBackend


def main():
    import ihm.brain.causal as causal
    backend = IBMBackend()
    sites = np.array([[0., 0., 0.]])
    receipts = []
    for kind in ('rapid', 'slow', 'twitch', 'activation'):
        for semantics in ('direct', 'drift'):
            model = causal.CausalIBM(backend, kind=kind, sites_m=sites, semantics=semantics)
            assert model.audit['donor_relative_error'] < 1e-9
            assert max(model.audit['poles_real_per_s']) < 0
            # Independent steady sinusoidal response of the actual ZOH matrices.
            # The finite hold introduces O(omega*dt) phase error against continuous H.
            omega = np.geomspace(.1, 1000., 31)
            dt = 1e-7
            ad, bd, c, d, _ = signal.cont2discrete(
                (model._A, model._B, model._C, model._D), dt)
            sampled = np.array([(c @ np.linalg.solve(np.exp(1j*w*dt)*np.eye(len(ad))-ad, bd)+d).item()
                                for w in omega])
            donor = backend.transfer(kind, SimpleNamespace(omega=omega, k=len(omega)))
            if semantics == 'drift':
                donor = donor/(1+1j*omega)
            np.testing.assert_allclose(sampled, donor, rtol=7e-5, atol=1e-9)
            receipts.append(model.audit)
    def twitch(**kwargs):
        return causal.CausalIBM(backend, kind='twitch', sites_m=sites,
                               parameters={'contraction_time_s': .05}, **kwargs)
    model = twitch(delay_s=.037)
    assert model.advance([1.], .02)[0] == 0
    np.testing.assert_allclose(model.advance([1.], .017), 0., atol=1e-28)
    value = model.advance([1.], .05)[0]
    np.testing.assert_allclose(value, 1 - 2 / np.e, atol=1e-12)
    saved = model.checkpoint()
    continuation = model.advance([0.], .1)
    model.restore(saved)
    np.testing.assert_array_equal(model.advance([0.], .1), continuation)
    restart = twitch(delay_s=.037)
    restart.restore(json.loads(json.dumps(saved)))
    np.testing.assert_array_equal(restart.advance([0.], .1), continuation)
    intact = restart.checkpoint()
    corrupt = json.loads(json.dumps(intact))
    corrupt['state'][0][0] = float('nan')
    try:
        restart.restore(corrupt)
    except ValueError:
        pass
    else:
        raise AssertionError('nonfinite checkpoint accepted')
    assert restart.checkpoint() == intact
    corrupt = json.loads(json.dumps(intact))
    corrupt['queue'] = [(intact['time_s'] + 100., [1.])]
    try:
        restart.restore(corrupt)
    except ValueError:
        pass
    else:
        raise AssertionError('impossible future queue timestamp accepted')
    assert restart.checkpoint() == intact
    coarse, fine = twitch(), twitch()
    a = coarse.advance([2.], .05)
    fine.advance([2.], .025)
    np.testing.assert_allclose(fine.advance([2.], .025), a, atol=1e-12)
    # Exact rectangular pulse equals the difference of two analytic step responses.
    pulse = twitch()
    pulse.advance([1.], .02)
    y = pulse.advance([0.], .08)[0]
    step = lambda t: 1 - np.exp(-t/.05) * (1+t/.05)
    np.testing.assert_allclose(y, step(.1)-step(.08), atol=1e-12)
    rapid = causal.CausalIBM(backend, kind='rapid', sites_m=sites)
    assert abs(rapid.advance([1.], 2.)[0]) < 1e-10
    assert rapid.advance([0.], .01)[0] < 0  # signed release, no rectification
    for operation in (
        lambda: causal.CausalIBM(backend, sites_m=np.empty((0, 3))),
        lambda: causal.CausalIBM(backend, sites_m=sites, parameters={'tau_adapt_s': -1}),
        lambda: twitch().restore(saved),  # delay is part of checkpoint identity
        lambda: causal.CausalIBM(backend, kind='twitch', sites_m=[[1., 0., 0.]],
                                delay_s=.037).restore(saved),
        lambda: causal.CausalIBM(backend, kind='twitch', sites_m=sites,
                    parameters={'contraction_time_s': .06}, delay_s=.037).restore(saved),
        lambda: twitch().advance([float('nan')], .01),
    ):
        try:
            operation()
        except ValueError:
            pass
        else:
            raise AssertionError('invalid operation accepted')
    print(json.dumps({'status': 'pass', 'receipts': receipts}, indent=2))


if __name__ == '__main__':
    main()
