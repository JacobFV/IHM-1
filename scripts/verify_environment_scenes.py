"""Sustained environment integration against a fixed canonical body fixture."""
import json,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from ihm.assembly.environment_dynamics import EnvironmentDynamics
from ihm.assembly.articulated import CanonicalRegistration
root=Path(__file__).resolve().parents[1]
payload=json.loads((root/'data/derived/canonical/mechanics.json').read_bytes())
specs={e['id']:e for e in payload['entities']};groups={}
for name,row in payload['registration'].items():
    bones=[specs[i] for i in row.get('canonical_bones',[]) if i in specs]
    if bones:groups[name]={'canonical_bones':[e['id'] for e in bones],
        'bounds_min_m':np.min([e['bounds_m']['min'] for e in bones],axis=0),'bounds_max_m':np.max([e['bounds_m']['max'] for e in bones],axis=0)}
registration=SimpleNamespace(specs=specs,groups=groups)
registration._ranking=lambda point:CanonicalRegistration._ranking(registration,point)
entities={i:{'centroid_m':e['centroid_m'],'rotation_matrix':np.eye(3).tolist()} for i,e in specs.items()}
catalog=json.loads((root/'data/derived/environment-catalogue-v1/catalogue.json').read_bytes());results=[]
for scene in catalog['scenes']:
    start=time.monotonic();base={'bed':'supine','floor':'upright'}[scene['base_environment']]
    owner=EnvironmentDynamics(root,base,{'scene':scene['id'],'objects':[]},registration)
    peak=0
    for _ in range(150):
        ports=owner.advance(.02,entities);peak=max(peak,len(ports))
    row={'scene':scene['id'],'duration_s':owner.time_s,'peak_contact_count':peak,'wall_s':time.monotonic()-start,'passed':True}
    print(json.dumps(row),flush=True);results.append(row)
Path('test-results/environment-upgrade').mkdir(parents=True,exist_ok=True)
Path('test-results/environment-upgrade/scene-stability.json').write_text(json.dumps({'scope':'3 s fixed-body fixture; excludes articulated feedback stability','results':results},indent=2))
