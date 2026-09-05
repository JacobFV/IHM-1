#!/usr/bin/env python3
"""Executable counterexamples for the bounded, real IBM backend."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np


def main():
    import importlib.util
    assert importlib.util.find_spec('ihm.brain') is not None, 'real IBM backend is missing'
    from ihm.brain.ibm_backend import IBMBackend
    from ihm.brain.port_mapping import pressure_to_indentation_um
    backend = IBMBackend()
    assert backend.identity['files'], 'artifact has no pinned package bytes'
    np.testing.assert_allclose(pressure_to_indentation_um([100], stiffness_pa_per_m=1e6), [100])
    try:
        pressure_to_indentation_um([100], stiffness_pa_per_m=0)
    except ValueError:
        pass
    else:
        raise AssertionError('zero stiffness accepted')
    try:
        backend.materialize_tactile([])
    except ValueError:
        pass
    else:
        raise AssertionError('empty materialization accepted')
    model = backend.materialize_tactile([[0, 0, 0], [.01, 0, 0]])
    assert model.layout['transduction.mechanoreceptor'].n_sites == 2
    assert model.implementations['transduction'].name == 'mechanoreceptor_rapid'
    assert backend.materialize_tactile([[0,0,0]], kind='slow').implementations['transduction'].name == 'mechanoreceptor_slow'
    pulse = np.zeros((2, 512)); pulse[:, 80:120] = 1
    errors = {}
    for kind in ('rapid', 'slow', 'twitch', 'activation'):
        result = backend.run_window(pulse, kind=kind, delay_s=.01)
        basis = backend.basis(512, .002)
        transfer = backend.transfer(kind, basis)
        # Independent frequency-domain solution of donor's drift equation.
        expected = basis.synthesize(basis.analyze(pulse) * transfer *
                                    np.exp(-1j*basis.omega*.01) / (1+1j*basis.omega))
        np.testing.assert_allclose(result.values, expected, atol=1e-10, rtol=1e-8)
        assert result.residual < 1e-8
        pure = basis.synthesize(basis.analyze(pulse)*transfer)
        assert np.max(np.abs(pure-result.values)) > 1e-5, 'drift mistaken for direct transfer'
        errors[kind] = float(np.max(np.abs(result.values-expected)))
    assert abs(backend.transfer('rapid', basis)[0]) == 0
    assert backend.transfer('slow', basis)[0].real > 0
    # Independent source-function parity at donor prior medians.
    w = 1j*basis.omega
    np.testing.assert_allclose(backend.transfer('rapid', basis),
                              (w*.02/(1+w*.02))/(1+w*.002), atol=1e-12)
    np.testing.assert_allclose(backend.transfer('slow', basis),
                              (.4+.6*w*2/(1+w*2))/(1+w*.01), atol=1e-12)
    try:
        backend.transfer('rapid', basis, {'threshold_um': 2})
    except ValueError:
        pass
    else:
        raise AssertionError('unused donor prior accepted as implemented parameter')
    # Pinned package bytes, not a commit label, determine source identity.
    import tempfile, shutil, json
    with tempfile.TemporaryDirectory() as directory:
        artifact = Path(directory)
        shutil.copytree(backend.artifact_dir, artifact, dirs_exist_ok=True)
        path = artifact/'source/ibm/runtime/step.py'
        path.write_bytes(path.read_bytes()+b'\n# tamper\n')
        try:
            IBMBackend(artifact)
        except ValueError as error:
            assert 'hash mismatch' in str(error)
        else:
            raise AssertionError('changed package source accepted')
    with tempfile.TemporaryDirectory() as directory:
        artifact=Path(directory)
        shutil.copytree(backend.artifact_dir,artifact,dirs_exist_ok=True)
        shadow=artifact/'source/ibm/frames/__init__.py'
        shadow.parent.mkdir(exist_ok=True)
        shadow.write_text('raise RuntimeError("unmanifested package executed")\n')
        try:IBMBackend(artifact)
        except ValueError as error:assert 'unmanifested' in str(error)
        else:raise AssertionError('unmanifested shadow package accepted')
    print({'status': 'pass', 'artifact': backend.identity['package_sha256'], 'errors': errors,
           'semantics': 'periodic single-window drift; no causal streaming claim'})


if __name__ == '__main__':
    main()
