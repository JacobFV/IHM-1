"""cycle generation, fluid conservation and supine cardiopulmonary mechanics."""
import numpy as np
from scipy.signal import find_peaks
from ihm.materialize.cardiopulmonary import cardiopulmonary_model

m = cardiopulmonary_model(subject='supine-example')
result = m.forecast(np.linspace(0, 30, 3001))
x = result['state']; signals = result['signals']; time = result['time']
assert np.all(np.isfinite(x)) and np.min(x[:, :8]) > 0
assert np.max(np.abs(x[:, :8].sum(axis=1)-x[0, :8].sum())) < 1e-5
assert np.max(np.abs(signals['net_circulatory_flow_mL_s'])) < 1e-8
late = time > 15
beats, _ = find_peaks(signals['aortic_flow_mL_s'][late], prominence=10, distance=50)
breaths, _ = find_peaks(signals['lung_volume_L'][late], prominence=.1, distance=200)
assert 14 <= len(beats) <= 19, len(beats)
assert 2 <= len(breaths) <= 4, len(breaths)
assert np.ptp(signals['systemic_arterial_pressure_mmHg'][late]) > 10
assert .2 < np.ptp(signals['lung_volume_L'][late]) < 1.5
assert np.min(signals['airflow_L_s'][late]) < 0 < np.max(signals['airflow_L_s'][late])
assert result['posture'] == 'supine' and result['validated_biology'] is False
print('verified: heartbeat, inspiration/expiration, pulsatile pressure, closed blood-volume conservation, supine posture')

from ihm.processes.cardiopulmonary import Parameters
for kwargs in ({'left_ventricular_emax_mmHg_mL': .01}, {'right_ventricular_emax_mmHg_mL': .01}):
    try:
        Parameters(**kwargs)
    except ValueError:
        pass
    else:
        raise AssertionError('inverted elastance range accepted')

assert np.min(signals['flow_mL_s'][:, [0,1,4,5]]) >= 0
custom = cardiopulmonary_model(parameters=Parameters(heart_rate_bpm=90, respiratory_rate_per_min=15))
custom_result = custom.forecast(np.linspace(0, 12, 1201))
assert np.isclose(custom_result['state'][-1, 9], 18)
assert np.isclose(custom_result['state'][-1, 10], 3)
import tempfile
from pathlib import Path
from ihm.materialize.cardiopulmonary import CardiopulmonaryModel
with tempfile.TemporaryDirectory() as directory:
    file = Path(directory)/'cycle.npz'
    m.advance(2); m.save(file)
    restored = CardiopulmonaryModel.load(file)
    assert np.allclose(m.forecast([3])['state'], restored.forecast([3])['state'])
print('verified: custom pacing, nonnegative valves, parameter validation, nonlinear save/load continuation')
