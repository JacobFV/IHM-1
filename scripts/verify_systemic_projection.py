"""Check actual native-to-thorax projection, energy, cadence and spectrum limits."""
import json
from pathlib import Path
import tempfile
import numpy as np
from ihm.assembly.systemic_projection import project_systemic
from ihm.app.experiments import read_experiment

root=Path(__file__).resolve().parents[1]
rest=read_experiment(root,'systemic-rest')
apnea=read_experiment(root,'systemic-apnea')
for data in (rest,apnea):
    assert data['has_body_projection'] and len(data['frames'])==1801
    assert data['projection_audit']['maximum_energy_balance_residual_j']<1e-9
    initial=data['frames'][0]['systemic_values']['lung_volume_ml']
    for frame in data['frames'][1:]:
        assert abs(frame['respiration']['time_s']-frame['time_s'])<1e-8
        assert abs(frame['respiration']['volume_change_m3']-(frame['systemic_values']['lung_volume_ml']-initial)*1e-6)<1e-12
    assert max(data['spectra']['frequency_hz'])<=.5/data['clock']['sample_interval_s']+1e-10
steady=[f['respiration']['diaphragm_descent_m'] for f in apnea['frames'] if 70<=f['time_s']<=85]
recovered=[f['respiration']['diaphragm_descent_m'] for f in apnea['frames'] if 150<=f['time_s']<=180]
assert np.ptp(steady)<1e-6 and np.ptp(recovered)>.001
# A deliberately sparse *observational* reduction must never invent chest cycles.
source=next(root/p for p in rest['runtime_sources'] if p.endswith('/rest/systemic.json'))
coarse=json.loads(source.read_text())
coarse['frames']=coarse['frames'][::300]
coarse['configuration']['sample_interval_s']=30
temp=Path(tempfile.mkdtemp(prefix='projection-',dir=root/'data/derived/audits'))/'coarse.json'
temp.write_text(json.dumps(coarse))
projected=project_systemic(root,temp)
assert not projected['has_body_projection']
assert all(not frame['entities'] and 'respiration' not in frame for frame in projected['frames'])
assert max(projected['spectra']['frequency_hz'])<=1/60+1e-12
print(json.dumps({'passed':True,'apnea_diaphragm_excursion_m':float(np.ptp(steady)),
                 'recovered_diaphragm_excursion_m':float(np.ptp(recovered)),
                 'energy_residual_j':apnea['projection_audit']['maximum_energy_balance_residual_j']},indent=2))
