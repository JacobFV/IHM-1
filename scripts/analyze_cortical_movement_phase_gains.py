#!/usr/bin/env python3
"""Audit slow endpoint gains under pure40ms delay; no interpolation/native claim."""
import argparse,hashlib,json,sys,time,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.linalg import solve_discrete_are,expm


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(node):
 return (node.tag,tuple(sorted(node.attrib.items())),(node.text or '').strip(),tuple(canonical(child) for child in node if child.tag!='default_activation'))


def spectrum(Ad,Bd,K,ticks=4):
 n,m=Bd.shape
 if ticks:
  F=np.zeros((n+m*ticks,n+m*ticks));F[:n,:n]=Ad;F[:n,n:n+m]=Bd
  for i in range(ticks-1):F[n+i*m:n+(i+1)*m,n+(i+1)*m:n+(i+2)*m]=np.eye(m)
  F[-m:,:n]=-K
 else:F=Ad-Bd@K
 eig=np.linalg.eigvals(F);mag=abs(eig)
 return {'spectral_radius':float(max(mag)),'unstable_modes':int(np.count_nonzero(mag>1+1e-8)),
  'numerically_neutral_modes':int(np.count_nonzero(abs(mag-1)<=1e-8)),
  'largest_non_neutral_magnitude':float(max(mag[abs(mag-1)>1e-8],default=0)),'bounded_with_neutral_tolerance':bool(max(mag)<=1+1e-8)}


def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',default='data/runtime/motor-learning/cortical-movement-phase-gains-20260908');args=ap.parse_args()
 out=ROOT/args.output;out.mkdir(parents=True,exist_ok=False);bundle=ROOT/'data/models/engineering_weight_transfer_v1'
 manifest=json.loads((bundle/'manifest.json').read_text());equivalence=json.loads((bundle/'provenance/model_equivalence.json').read_text())
 required=[bundle/'provenance/baseline_policy.npz',bundle/'target_policy.npz',bundle/'target.npz',bundle/'provenance/model_equivalence.json',
  bundle/'provenance/baseline_registration.json',bundle/'provenance/target_registration.json',bundle/'provenance/target_model.osim']
 for path in required:
  assert manifest['files'][str(path.relative_to(ROOT))]==sha(path), 'Promoted file digest differs'
 baseline_model=ROOT/'data/models/engineering_stance_v1/model.osim';target_model=bundle/'provenance/target_model.osim'
 assert sha(baseline_model)==equivalence['baseline_model_sha256'] and sha(target_model)==equivalence['target_model_sha256']
 exact_equivalence=canonical(ET.parse(baseline_model).getroot())==canonical(ET.parse(target_model).getroot())
 assert exact_equivalence and equivalence['same_except_default_activation']
 data={}
 for phase,path in [('baseline',required[0]),('target',required[1])]:
  with np.load(path,allow_pickle=False) as z:data[phase]={k:z[k].copy() for k in z.files}
 with np.load(bundle/'target.npz',allow_pickle=False) as z:target={k:z[k].copy() for k in z.files}
 assert np.array_equal(target['x_target'],data['target']['x0']) and np.array_equal(target['u_target'],data['target']['u0'])
 assert np.array_equal(data['baseline']['state_names'],data['target']['state_names']) and np.array_equal(data['baseline']['muscle_names'],data['target']['muscle_names'])
 report={'schema':'ihm.cortical-movement-phase-gains.v1','source_sha256':sha(__file__),'bundle_manifest_sha256':sha(bundle/'manifest.json'),
  'physical_model_equivalence':{'exact_xml_equal_ignoring_only_default_activation':exact_equivalence,
   'baseline_model_sha256':sha(baseline_model),'target_model_sha256':sha(target_model),
   'note':'Distinct artifact identities and distinct operating points; same physical XML except default activation. Baseline runtime model remains required; target model is provenance only.'},
  'target_reference_exactly_matches_target_policy':True,'dt_s':.01,'pure_delay_ticks':4,'phases':{},
  'scope':['All258states and98muscles retained; original Q_discrete unchanged, R_discrete scaled per endpoint',
   'Offline teacher gain construction only; no production controller or live LQR motor output',
   'Pure40ms input delay is not the actual neural transfer function',
   'Independent endpoint spectral checks do not prove interpolation, trajectory tracking, nonlinear stability, or movement',
   'Numerically neutral modes near1 prevent a strict asymptotic-stability claim']}
 started=time.monotonic();retained=out/'sources';retained.mkdir()
 for path in required:
  (retained/path.name).write_bytes(path.read_bytes())
 (retained/'manifest.json').write_bytes((bundle/'manifest.json').read_bytes())
 for phase,path in [('baseline',required[0]),('target',required[1])]:
  d=data[phase];Ad,Bd=d['Ad'],d['Bd'];n,m=Bd.shape
  assert (n,m)==(258,98) and float(d['dt_s'])==.01 and abs(float(d['target_mass_kg'])-77.6122029)<1e-9
  assert str(d['model_sha256'].item())==equivalence[phase+'_model_sha256']
  transition=expm(np.block([[d['A'],d['B']],[np.zeros((m,n+m))]])*.01)
  assert np.allclose(transition[:n,:n],Ad,atol=1e-12,rtol=1e-12) and np.allclose(transition[:n,n:],Bd,atol=1e-12,rtol=1e-12)
  reserve=d['u0']-d['minimum_activation'];rows=[]
  item={'source_artifact':str(path.relative_to(ROOT)),'source_sha256':sha(path),'model_sha256':str(d['model_sha256'].item()),
   'state_count':n,'muscle_count':m,'target_mass_kg':float(d['target_mass_kg']),
   'minimum_excitation_reserve':float(reserve.min()),'reserve_limiting_muscle':str(d['muscle_names'][np.argmin(reserve)]),
   'reserve_below_1e_minus5_count':int(np.count_nonzero(reserve<1e-5)),
   'original_40ms':spectrum(Ad,Bd,d['K']),'original_undelayed':spectrum(Ad,Bd,d['K'],0),'candidates':rows}
  for scale in (300.,1000.,3000.,10000.,100000.,1000000.):
   Q=d['Q_discrete'];R=d['R_discrete']*scale
   try:
    P=solve_discrete_are(Ad,Bd,Q,R);K=np.linalg.solve(R+Bd.T@P@Bd,Bd.T@P@Ad)
    delay=spectrum(Ad,Bd,K);candidate=out/f'{phase}_R{int(scale)}.npz'
    np.savez_compressed(candidate,**dict(d,K=K,Q_discrete=Q,R_discrete=R,design_delay_ticks=np.array(4),
     phase_name=np.array(phase),source_policy_sha256=np.array(sha(path)),R_scale=np.array(scale)))
    terms=(Ad.T@P@Ad,-P,-Ad.T@P@Bd@K,Q)
    residual=np.linalg.norm(sum(terms))
    row={'R_scale':scale,'delayed_40ms':delay,'undelayed':spectrum(Ad,Bd,K,0),'artifact':str(candidate.relative_to(ROOT)),
     'sha256':sha(candidate),'dare_relative_residual':float(residual/np.linalg.norm(Q)),
     'dare_term_relative_residual':float(residual/sum(np.linalg.norm(t) for t in terms)),
     'riccati_condition_number':float(np.linalg.cond(P)),
     'optimality_note':'Large Q-relative residual means optimality is not established; spectral values measure returned K directly regardless of Riccati conditioning.'}
   except Exception as exc:row={'R_scale':scale,'error':type(exc).__name__+': '+str(exc)}
   rows.append(row);print(json.dumps({'phase':phase,**row}),flush=True)
  stable=[r for r in rows if r.get('delayed_40ms',{}).get('bounded_with_neutral_tolerance')]
  item['selected']=min(stable,key=lambda r:abs(np.log(r['R_scale']/1000))) if stable else None
  report['phases'][phase]=item;(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 report['both_endpoints_have_40ms_bounded_candidate']=all(v['selected'] is not None for v in report['phases'].values())
 report['wall_seconds']=time.monotonic()-started
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');(out/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
 print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
