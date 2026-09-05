"""Lossless selection of renderer observables from a full body trajectory."""
from copy import deepcopy


def display_trajectory(trajectory, *, source_sha256):
    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    result = {k: deepcopy(v) for k, v in trajectory.items() if k != 'frames'}
    result['projection'] = {
        'kind': 'exact_observable_selection', 'source_sha256': source_sha256,
        'precision_reduction': False, 'temporal_resampling': False,
        'identity_convention': 'Omitted rotation and deformation matrices are exactly identity.',
        'omitted': 'Internal stress, forces, compartments and neural state; retained in full trajectory.',
    }
    result['frames'] = []
    for frame in trajectory['frames']:
        selected = {'time_s': frame['time_s'], 'physiology': frame.get('physiology', {}),
                    'entities': {id: {k: v for k, v in state.items() if k in
                        ('translation_m', 'rotation_matrix', 'deformation_gradient') and
                        (k == 'translation_m' or v != identity)}
                        for id, state in frame['entities'].items()}}
        if 'skin_field' in frame.get('respiration', {}):
            selected['respiration'] = {'skin_field': frame['respiration']['skin_field']}
        result['frames'].append(selected)
    return result
