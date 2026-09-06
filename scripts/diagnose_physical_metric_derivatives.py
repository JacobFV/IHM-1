"""Held-data attribution of frozen inverse-mass derivative error and shoulder forces."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.linalg import null_space
from support_feasible_physical_step import physical_step
from prepare_native_force_scaling import model_info,GAUGES,ROOT
from solve_support_physical_metric import physical_gate


def diagnose():
    folder=ROOT/'data/derived/support-physical-root-xgpiykt1';log=folder/'observations.jsonl'
    records=[json.loads(s) for s in log.read_text().splitlines()];lookup={r['label']:r['response'] for r in records}
    base=lookup['trial:4:0:correction:0'];names=base['independent_names'];free=[n for n in names if n not in GAUGES];idx=[names.index(n) for n in free];x=np.array(base['independent_q'])[idx]
    a=np.array(base['udot']);r=np.array(base['constrained_residual']);inverse=np.array(base['inverse_mass_matrix']);D=[];Yr=[];Ya=[];Ym=[]
    for n in free:
        v=lookup['jacobian:5:'+n];D.append(np.array(v['independent_q'])[idx]-x);Yr.append(np.array(v['constrained_residual'])-r);Ya.append(np.array(v['udot'])-a);Ym.append((np.array(v['inverse_mass_matrix'])-inverse).ravel())
    Jr=np.linalg.solve(D,Yr).T;Ja=np.linalg.solve(D,Ya).T
    derivatives=np.linalg.solve(D,Ym).reshape(len(free),len(a),len(a))
    frozen=-inverse@Jr;missing=-np.stack([m@r for m in derivatives],axis=1);corrected=frozen+missing
    trials=[]
    for label,v in lookup.items():
        if not label.startswith('trial:5:'):continue
        d=np.array(v['independent_q'])[idx]-x;actual=np.array(v['udot'])
        trials.append(dict(label=label,maximum_actual_step=float(np.max(abs(d))),actual_cost=float(.5*actual@actual),
            frozen_prediction_cost=float(.5*np.sum((a+frozen@d)**2)),full_prediction_cost=float(.5*np.sum((a+Ja@d)**2)),
            frozen_error_norm=float(np.linalg.norm(actual-a-frozen@d)),full_error_norm=float(np.linalg.norm(actual-a-Ja@d)),
            omitted_mass_term_norm=float(np.linalg.norm(missing@d)),finite_product_rule_error_norm=float(np.linalg.norm((corrected-Ja)@d))))
    case=json.load(open(ROOT/'data/derived/lumbar-supine-reference-1qex2x9i/manifest.json'));surface=json.load(open(ROOT/case['surface_manifest']));height=float(np.ptp(np.load(ROOT/surface['arrays_path'])['reference_points_source_m'][:,1]));weight=case['mass_kg']*9.81
    gate=physical_gate(base,weight,height);rootidx=[base['mobility_names'].index(n) for n in ('pelvis_tx','pelvis_tilt','pelvis_rotation')];A=Jr[rootidx]/np.array([weight,weight*height,weight*height])[:,None]
    bounds,_=model_info(ROOT/'data/derived/effective-potential-build-_jnt2x88/inputs/subject_walk_scaled.osim');bounds=np.array([bounds[n] for n in free]);directions=[]
    for radius in (.00375,.001875,.0009375,.00046875,.000234375):
        d,diag=physical_step(a,Ja,x,bounds,gate['support_constraints'],A,radius)
        active=diag['active_source_or_box_columns'];c=np.array(diag['coordinate_preconditioner']);scale=max(1.,np.linalg.norm(a));C=Ja*c/scale;b=a/scale;z=d/c
        constraints=np.vstack([A*c,np.eye(len(free))[active]]) if active else A*c
        tangent=null_space(constraints);gradient=C.T@(b+C@z)
        stationarity=float(np.linalg.norm(tangent.T@gradient))
        directions.append(dict(radius=radius,diagnostic=diag,active_manifold_stationarity_norm=stationarity,source_bound_active=[free[i] for i in active if min(abs(x[i]+d[i]-bounds[i]))<1e-8],step=dict(zip(free,d.tolist()))))
    final=json.load(open(folder/'last_accepted.json'))['response'];shoulders=[]
    for n in (stem+side for side in ('r','l') for stem in ('arm_flex_','arm_add_','arm_rot_')):
        i=final['independent_names'].index(n);k=final['mobility_names'].index(n)
        muscles=[dict(name=m['name'],force_nm=m['tendon_force_n']*m['moment_arms_m'][i],moment_arm_m=m['moment_arms_m'][i],activation=m['activation']) for m in final['muscles'] if abs(m['moment_arms_m'][i])>1e-10]
        parts={part:final[part+'_generalized_forces'][k] for part in ('gravity','contact','joint')};parts['muscle']=sum(m['force_nm'] for m in muscles)
        shoulders.append(dict(coordinate=n,native_tree_residual_nm=final['tree_residual'][k],native_acceleration_rad_s2=final['udot'][k],forces_nm=parts,muscles=muscles,closure_error_nm=final['tree_residual'][k]+sum(parts.values())))
    return dict(initial_cost=float(.5*a@a),actual_difference_condition=float(np.linalg.cond(D)),
        frozen_vs_full_jacobian_relative_error=float(np.linalg.norm(frozen-Ja)/np.linalg.norm(Ja)),
        product_rule_corrected_vs_full_relative_error=float(np.linalg.norm(corrected-Ja)/np.linalg.norm(Ja)),
        trials=trials,recomputed_full_jacobian_directions=directions,shoulder_force_attribution=shoulders,
        interpretation='Directional frozen-mass derivative error explains failed predictions; no causal claim that limited muscle coverage or wrap mismatch explains stagnation. New QP predictions remain unverified native candidates.',
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (log,folder/'last_jacobian.json',folder/'last_accepted.json',Path(__file__))})
if __name__=='__main__':print(json.dumps(diagnose(),indent=2,allow_nan=False))
