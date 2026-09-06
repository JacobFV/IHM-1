"""Offline new physical QP using retained actual98 force Jacobian and native MInv."""
from pathlib import Path
import json,hashlib
import numpy as np
from prepare_native_force_scaling import model_info,GAUGES,ROOT
from support_feasible_physical_step import physical_step,native_metric
from solve_native_force_supine import physical_gate

def diagnose():
    folder=ROOT/'data/derived/native-effective-potential-_livswt7'
    records=[json.loads(s)['response'] for s in (folder/'observations.jsonl').read_text().splitlines()];base=records[0]
    metriclog=ROOT/'data/derived/native-physical-metric-04jdqqbn/observations.jsonl';metric=json.loads(metriclog.read_text().splitlines()[0])['response']
    qerror=float(np.max(np.abs(np.array(base['independent_q'])-metric['independent_q'])))
    forceerror=float(np.max(np.abs(np.array(base['constrained_residual'])-metric['constrained_residual'])))
    if qerror>1e-8 or forceerror>1e-6:raise ValueError('Retained metric/Jacobian anchor differs')
    names=base['independent_names'];free=[n for n in names if n not in GAUGES];idx=[names.index(n) for n in free]
    bounds,_=model_info(ROOT/'data/derived/effective-potential-build-_jnt2x88/inputs/subject_walk_scaled.osim');D=[];Y=[];Ya=[]
    for i in range(len(names)):
        a,b=records[1+2*i:3+2*i]
        D.append(2*(np.array(a['independent_q'])-base['independent_q'])-.5*(np.array(b['independent_q'])-base['independent_q']))
        Y.append(2*(np.array(a['constrained_residual'])-base['constrained_residual'])-.5*(np.array(b['constrained_residual'])-base['constrained_residual']))
        Ya.append(2*(np.array(a['udot'])-base['udot'])-.5*(np.array(b['udot'])-base['udot']))
    Jr=np.linalg.solve(D,Y).T[:,idx];Ja=np.linalg.solve(D,Ya).T[:,idx]
    inverse=np.array(metric['inverse_mass_matrix']);B=-inverse@Jr;a=native_metric(metric)
    case=json.load(open(ROOT/'data/derived/lumbar-supine-reference-1qex2x9i/manifest.json'));surface=json.load(open(ROOT/case['surface_manifest']))
    height=float(np.ptp(np.load(ROOT/surface['arrays_path'])['reference_points_source_m'][:,1]));weight=case['mass_kg']*9.81
    gate=physical_gate(metric,weight,height);rootnames=('pelvis_tx','pelvis_tilt','pelvis_rotation');rootindices=[base['mobility_names'].index(n) for n in rootnames]
    A=Jr[rootindices]/np.array([weight,weight*height,weight*height])[:,None]
    rows=[]
    for radius in (.03,.015,.0075,.00375):
        d,diag=physical_step(a,B,np.array(base['independent_q'])[idx],[bounds[n] for n in free],gate['support_constraints'],A,radius)
        rows.append(dict(radius=radius,step=dict(zip(free,d.tolist())),diagnostic=diag,full_acceleration_jacobian_predicted_cost=float(.5*np.sum((a+Ja@d)**2))))
    return dict(anchor_q_difference=qerror,anchor_force_difference=forceerror,native_inverse_identity_passed=True,initial_acceleration_cost=float(.5*a@a),
        frozen_mass_vs_full_acceleration_jacobian_relative_difference=float(np.linalg.norm(B-Ja)/np.linalg.norm(Ja)),rows=rows,
        scope='Offline bounded linear predictions only, including neglected configuration derivative of MInv; no native candidate acceptance',
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (metriclog,folder/'observations.jsonl',Path(__file__))})
if __name__=='__main__':print(json.dumps(diagnose(),indent=2,allow_nan=False))
