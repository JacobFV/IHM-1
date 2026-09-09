"""Verify retained gait references cannot smuggle root force or silent channels."""
from pathlib import Path
import sys
import tempfile
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ihm.native.gait_reference import GaitReference, _read


def must_raise(function):
    try:
        function()
    except ValueError:
        return
    raise AssertionError('Expected strict rejection')


def main():
    reference = GaitReference.load()
    targets = reference.joint_targets(.5)
    assert len(targets) == 23 and not any('pelvis' in name for name in targets)
    assert all(not any(token in name for token in ('reserve', 'corrector')) for name in targets)
    assert not any(name.endswith('_beta') for name in targets)
    excluded = reference.provenance['excluded_dependent_coordinates']
    assert set(excluded) == {'knee_angle_r_beta', 'knee_angle_l_beta'}
    assert all(item['maximum_error_vs_primary_degrees'] < 1e-8 for item in excluded.values())
    seed = reference.muscle_seed(.6)
    assert sum(value is not None for value in seed.values()) == 72
    assert {name for name, value in seed.items() if value is None} == {
        f'addmag{part}_{side}' for part in ('Dist', 'Isch', 'Mid', 'Prox') for side in ('r', 'l')}
    assert reference.muscle_seed(.6, ['new_lumbar_r']) == {'new_lumbar_r': None}
    for name, source in reference.provenance['muscle_mapping'].items():
        assert name[-2:] == source[-2:]
    must_raise(lambda: reference.muscle_seed(.6, ['proprio:soleus_r']))
    must_raise(lambda: reference.joint_targets(.449))
    must_raise(lambda: reference.muscle_seed(2))
    must_raise(lambda: reference.sample_phase(1.01))
    must_raise(lambda: reference.sample_phase(float('nan')))
    # Independent source clocks must not masquerade as synchronized trials.
    start = reference.sample_phase(0)
    assert start['coordinate_source_time_s'] == .45
    assert start['muscle_source_time_s'] == .53
    assert reference.provenance['periodic'] is False
    assert reference.provenance['synchronized_trials'] is False
    midpoint = (reference.coordinate_data[3, 0] + reference.coordinate_data[4, 0]) / 2
    np.testing.assert_allclose(list(reference.joint_targets(midpoint).values()),
        (reference.coordinate_data[3, 1:] + reference.coordinate_data[4, 1:]) / 2)
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'bad.sto'
        path.write_text('inDegrees=no\nendheader\ntime q\n0 1\n0 2\n1 3\n')
        must_raise(lambda: _read(path))
        path.write_text('inDegrees=yes\nendheader\ntime q\n0 1\n1 3\n')
        must_raise(lambda: _read(path))
    print('PASS: joint-only references; 72 mapped/8 explicit missing muscles; strict clocks, units, mapping, interpolation, duplicate validation')


if __name__ == '__main__':
    main()
