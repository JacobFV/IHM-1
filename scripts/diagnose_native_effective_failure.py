"""Offline attribution of retained effective-potential gate failures."""
from pathlib import Path
import hashlib,json,math
import xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def fiber_primitive(normalized_length,lo,fiso,k,e):
    d=max(0.,normalized_length-1)
    return fiso*lo*(e/k*math.expm1(k*d/e)-d)/math.expm1(k)


def audit():
    folder=ROOT/'data/derived/native-effective-potential-_livswt7';records=[json.loads(s) for s in (folder/'observations.jsonl').read_text().splitlines()];base=records[0]['response'];names=base['independent_names'];lookup={r['label']:r['response'] for r in records};pairs=[(lookup[n+'_1'],lookup[n+'_2']) for n in names]
    model=ROOT/'data/derived/effective-potential-build-3a9juno_/inputs/subject_walk_scaled.osim';root=ET.parse(model).getroot();parameters={}
    for m in root.iter('Thelen2003Muscle'):parameters[m.get('name')]={k:float(m.findtext(k)) for k in ('optimal_fiber_length','max_isometric_force','KshapePassive','FmaxMuscleStrain')}
    def correction(m):
        if m['name'] not in parameters:return 0.
        p=parameters[m['name']];lo=p['optimal_fiber_length'];args=(lo,p['max_isometric_force'],p['KshapePassive'],p['FmaxMuscleStrain']);lf=m['fiber_length_m']
        return fiber_primitive(lf/lo,*args)-fiber_primitive(lf,*args)
    D=np.array([2*(np.array(a['independent_q'])-base['independent_q'])-.5*(np.array(b['independent_q'])-base['independent_q']) for a,b in pairs])
    def grad(get):return np.linalg.solve(D,np.array([2*(get(a)-get(base))-.5*(get(b)-get(base)) for a,b in pairs]))
    corrected_gradient=grad(lambda r:r['effective_energy_j']+sum(correction(m) for m in r['muscles']));original=json.loads((folder/'report.json').read_text());physical=np.array(original['physical_gradient']);rows=[]
    for i,m in enumerate(base['muscles']):
        g=grad(lambda r:r['muscles'][i]['passive_energy_j']+r['muscles'][i]['active_effective_energy_j']+correction(r['muscles'][i]));ma=np.array(m['moment_arms_m']);length=grad(lambda r:r['muscles'][i]['length_m'])
        rows.append(dict(name=m['name'],added_passive_fiber_energy_j=correction(m),corrected_energy_virtual_work_error=float(np.max(np.abs(g+m['tendon_force_n']*ma))),corrected_energy_length_work_error=float(np.max(np.abs(g-m['tendon_force_n']*length))),length_moment_arm_error_m=float(np.max(np.abs(length+ma)))))
    paths=[model,folder/'observations.jsonl',folder/'report.json',ROOT/'data/raw/mechanics/opensim-core/OpenSim/Actuators/Thelen2003Muscle.cpp',Path(__file__)]
    return dict(scope='Offline source-derived numerical merit correction only; original failed native receipt unchanged; no physiological library or ledger edit',
        defect='Thelen calcMusclePotentialEnergyInfo passes mli.fiberLength to calcfpefisoPE expecting normalized length',maximum_original_gradient_error=original['maximum_gradient_error'],maximum_corrected_gradient_error=float(np.max(np.abs(corrected_gradient-physical))),corrected_gradient_errors=dict(zip(names,(corrected_gradient-physical).tolist())),muscles=rows,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':print(json.dumps(audit(),indent=2,allow_nan=False))
