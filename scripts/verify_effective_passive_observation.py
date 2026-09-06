"""Derivative, raw-value preservation and retained observation fixtures."""
import json, math
from pathlib import Path
from effective_passive_energy_observation import PassiveEnergyObservation, fiber_primitive, ROOT

lo=.1; lf=.14; fiso=900.; k=4.; strain=.6; h=1e-7
assert fiber_primitive(lf,lo,fiso,k,strain)==0
slope=(fiber_primitive((lf+h)/lo,lo,fiso,k,strain)-fiber_primitive((lf-h)/lo,lo,fiso,k,strain))/(2*h)
assert math.isclose(slope,fiso*math.expm1(k*(lf/lo-1)/strain)/math.expm1(k),rel_tol=1e-8)
base=json.loads((ROOT/'data/derived/native-effective-potential-_livswt7/base.json').read_text())
raw=json.dumps(base,sort_keys=True)
port=PassiveEnergyObservation(ROOT/'data/derived/effective-potential-build-3a9juno_/inputs/subject_walk_scaled.osim')
result=port.observe(base)
assert len(result['muscles'])==98 and result['source_passive_correction_j']>0
assert json.dumps(base,sort_keys=True)==raw
assert math.isclose(result['corrected_passive_energy_j'],result['raw_native_passive_energy_j']+result['source_passive_correction_j'],abs_tol=1e-12)
for m,r in zip(base['muscles'],result['muscles']):
    if m['type']=='Millard2012EquilibriumMuscle':assert r['source_passive_correction_j']==0
print('PASS: source derivative, all98 raw observations preserved, correction separately exposed')
