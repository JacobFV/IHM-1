"""Offline scaling design from actual98-muscle/seam-omitted native samples."""
from pathlib import Path
import json,hashlib
import xml.etree.ElementTree as ET
import numpy as np
from verify_native_effective_potential import pullback
from native_force_trust_region import sensitivity_scales,bounded_step
ROOT=Path(__file__).resolve().parents[1]
GAUGES=('pelvis_ty','pelvis_tz','pelvis_list')

def model_info(path):
    model=ET.parse(path).getroot();bounds={c.get('name'):list(map(float,c.findtext('range').split())) for c in model.iter('Coordinate')};couplings={}
    for c in model.findall('.//ConstraintSet/objects/*'):
        coefficients=list(map(float,c.findtext('.//LinearFunction/coefficients').split()))
        if c.tag!='CoordinateCouplerConstraint' or len(coefficients)!=2 or coefficients[1]!=0:raise ValueError('Unsupported coupling')
        couplings[c.findtext('dependent_coordinate_name')]=(c.findtext('independent_coordinate_names'),coefficients[0])
    return bounds,couplings

def force(response,couplings):return pullback(response['mobility_names'],response['independent_names'],couplings,response['tree_residual'])

def prepare():
    folder=ROOT/'data/derived/native-effective-potential-_livswt7';log=folder/'observations.jsonl';model=ROOT/'data/derived/effective-potential-build-3a9juno_/inputs/subject_walk_scaled.osim'
    records=[json.loads(s)['response'] for s in log.read_text().splitlines()];base=records[0];names=base['independent_names'];bounds,couplings=model_info(model)
    f=force(base,couplings);D=[];Y=[]
    for i in range(len(names)):
        a,b=records[1+2*i:3+2*i]
        D.append(2*(np.array(a['independent_q'])-base['independent_q'])-.5*(np.array(b['independent_q'])-base['independent_q']))
        Y.append(2*(force(a,couplings)-f)-.5*(force(b,couplings)-f))
    J=np.linalg.solve(D,Y).T;indices=[i for i,n in enumerate(names) if n not in GAUGES];free=[names[i] for i in indices];B=np.array([bounds[n] for n in free]);J=J[:,indices]
    c,s,description=sensitivity_scales(J,B)
    q=np.array(base['independent_q'])[indices];d,diag=bounded_step(J,f,q,B,c,s,.01)
    paths=[log,model,Path(__file__),ROOT/'scripts/native_force_trust_region.py']
    return dict(scope='Offline scaling only from actual98-muscle/seam-omitted observations; no92 response reuse; no new equilibrium claim',
        independent_names=names,free_names=free,gauges=list(GAUGES),coordinate_units=['m' if n in ('pelvis_tx','pelvis_ty','pelvis_tz') else 'rad' for n in free],
        residual_units=['N' if n in ('pelvis_tx','pelvis_ty','pelvis_tz') else 'Nm' for n in names],
        coordinate_scale=c.tolist(),residual_scale=s.tolist(),basis=description,initial_scaled_cost=float(.5*np.sum((f/s)**2)),
        local_step=dict(zip(free,d.tolist())),diagnostic=diag,
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':print(json.dumps(prepare(),indent=2,allow_nan=False))
