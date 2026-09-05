"""numerical and integration checks; synthetic evidence is not biological validation."""
import tempfile
from pathlib import Path
import numpy as np
from ihm import body, materialize, Request
from ihm.forge.ingest import Evidence, SourceCard, load
from ihm.materialize.model import Model
from ihm.forge.fit import fit_affine, evaluate


def main():
    registry = body()
    m = materialize(registry, Request(('metabolic.glucose',), 'alice'))
    i = m.components.index('metabolic.glucose')
    old_mean, old_var = m.mean[i], m.cov[i, i]
    e = Evidence('lab', 'one', 'alice', 0, 'metabolic.glucose', 7., .25, 'measured')
    m.assimilate([e])
    expected = (old_mean / old_var + 7 / .25) / (1 / old_var + 1 / .25)
    assert np.isclose(m.mean[i], expected)
    assert np.isclose(m.cov[i, i], 1 / (1 / old_var + 1 / .25))
    for bad in [e, Evidence('lab', 'two', 'bob', 0, e.component, 7., .25, 'measured')]:
        before = m.mean.copy()
        try:
            m.assimilate([bad])
        except ValueError:
            pass
        else:
            raise AssertionError('duplicate or subject mismatch accepted')
        assert np.array_equal(m.mean, before)
    assert 'endocrine.insulin' in m.components
    assert 'mechanical.force' not in m.components
    try:
        materialize(registry, Request(('metabolic.glucose',), 'alice', max_states=1))
    except ValueError:
        pass
    else:
        raise AssertionError('budget ignored')
    result = m.forecast([0, 1, 10])
    assert np.all(np.isfinite(result['mean']))
    assert np.linalg.eigvalsh(m.cov).min() > -1e-9
    m.advance(2, {'metabolic.glucose': 8.})
    assert m.mean[i] == 8 and m.cov[i, i] == 0
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / 'model.npz'
        m.save(path)
        restored = Model.load(path)
        assert np.allclose(restored.forecast([3])['mean'], m.forecast([3])['mean'])
        try:
            restored.assimilate([e])
        except ValueError:
            pass
        else:
            raise AssertionError('restored model forgot evidence')
        csv = Path(d) / 'labs.csv'
        csv.write_text('id,subject,time,component,value,unit,variance\na,alice,0,metabolic.glucose,90,mg/dL,324\n')
        card = SourceCard('lab', 'measured', ('metabolic.glucose',), 'local')
        ev = load(csv, card, registry)
        assert ev[0].reference == 'local' and len(ev[0].data_sha256) == 64
        source_model = materialize(registry, Request(('metabolic.glucose',), 'alice'))
        source_model.assimilate(ev); source_model.save(path)
        assert Model.load(path).provenance['sources']['lab']['card']['reference'] == 'local'
        assert np.isclose(ev[0].value, 90 / 18.01559)
        assert np.isclose(ev[0].variance, 324 / 18.01559**2)
    rng = np.random.default_rng(12)
    x = rng.normal(size=(100, 2)); y = 2*x[:, 0] - 3*x[:, 1] + 4
    fitted = fit_affine(x, y, ('a', 'b'), 'c', 'fixture', ('train',))
    heldout_x = rng.normal(size=(40, 2))
    heldout_y = 2*heldout_x[:, 0] - 3*heldout_x[:, 1] + 4
    assert evaluate(fitted, heldout_x, heldout_y, ('heldout',))['rmse'] < 1e-8
    try:
        evaluate(fitted, x, y, ('train',))
    except ValueError:
        pass
    else:
        raise AssertionError('training subject leakage')
    print('verified: conditioning, units, isolation, deduplication, closure, budget, forecast, clamps, export, fit, leakage')


if __name__ == '__main__':
    main()

# Negative variance must be rejected even at small physical scales.
from ihm.runtime.gaussian import covariance
try:
    covariance([[-1e-10]], 1)
except ValueError:
    pass
else:
    raise AssertionError('negative small-scale variance accepted')
