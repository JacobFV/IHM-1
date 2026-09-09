"""Verify actual checkpoint adapter replay, validation and30ms cord latency."""
from copy import deepcopy
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import torch
from ihm.assembly.ibm_controller import IBMImplicitController
from ihm.assembly.sensorimotor_catalog import native_muscle_catalog


def main():
    torch.set_num_threads(1)
    c = IBMImplicitController.from_root(ROOT, muscle_catalog=native_muscle_catalog(ROOT), sites=128)
    initial = c.checkpoint()
    observation = {'time_s':0., 'foot_contact_force_n':{'r':0.,'l':0.},
        'muscles':{m:{'sensor_basis':'test perturbation on native observation ports',
                     'optimal_fiber_length_m':1.,'max_isometric_force_n':100.,
                     'fiber_length_m':1.1,'tendon_force_n':0.} for m in c.muscles}}
    assert c.step(.03, observation)['arc_max']['stretch'] == 0.
    observation['time_s'] = .03
    saved = c.checkpoint()
    first = c.step(.001, observation)
    assert first['arc_max']['stretch'] > 0.
    endpoint = c.checkpoint()
    c.restore(saved)
    assert first == c.step(.001, observation)
    assert endpoint == c.checkpoint()
    c.restore(saved)
    blocked = c.step(.001, observation, sensory_blocks=c.muscles)
    assert blocked['arc_max']['stretch'] == 0.
    c.restore(initial)
    observation['time_s'] = 0.
    for kwargs in ({'additional_sensory_inputs_hz':{'typo':1.}}, {'motor_blocks':['typo']}):
        try: c.step(.02, observation, **kwargs)
        except ValueError: pass
        else: raise AssertionError('Malformed input was accepted')
        assert c.checkpoint() == initial
    bad = deepcopy(initial)
    bad['cortical_state'][0][0][0] = float('nan')
    try: c.restore(bad)
    except ValueError: pass
    else: raise AssertionError('Nonfinite checkpoint accepted')
    assert c.checkpoint() == initial
    blocked = c.step(.02, observation, motor_blocks=c.muscles)
    assert all(v == 0. for v in blocked['motor_excitations'].values())
    assert any(v > 0. for v in blocked['cortical_commands'].values())
    assert set(blocked['motor_excitations']) == set(c.muscles)
    c.restore(initial)
    observation['time_s'] = 0.
    normal = c.step(.02, observation)
    c.restore(initial)
    hypoxic = c.step(.02, observation, physiology={'oxygen_saturation':.3})
    assert normal['cortical_commands'] != hypoxic['cortical_commands']
    c.restore(initial)
    warm = c.step(.02, observation, physiology={'core_temperature_C':40.})
    assert normal['cortical_commands'] != warm['cortical_commands']
    malformed = c.checkpoint()
    del malformed['delays']['stretch']
    intact = c.checkpoint()
    try: c.restore(malformed)
    except ValueError: pass
    else: raise AssertionError('Missing delay history accepted')
    assert c.checkpoint() == intact
    c.restore(initial)
    for m in observation['muscles'].values(): m['tendon_force_n'] = 20.
    result = c.step(.04, observation)
    assert all(result['arc_max'][name] > 0. for name in ('stretch','autogenic','renshaw','reciprocal'))
    print(json.dumps({'passed':True, 'checks':['exact replay','30ms delay','sensory flush',
        'strict inputs','checkpoint atomicity','motor block','oxygen and temperature causality','all four named arcs'],
        'source_checkpoint_sha256':c.identity['checkpoint_sha256'],
        'muscles':len(c.muscles),'spinal_mapped':len(c.muscles)-len(c.cord.unmapped),
        'stretch_at_30ms':first['arc_max']['stretch'],
        'scope':'Adapter native-observation fixture; actual world verification is separate'},indent=2))

if __name__ == '__main__': main()
