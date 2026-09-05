"""Fail-closed optimized-vs-forward acceptance for native 3D gait artifacts."""
import argparse,hashlib,json
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

def read_sto(path):
 lines=path.read_text().splitlines();start=lines.index('endheader');meta=dict(s.split('=',1) for s in lines[:start] if '=' in s)
 names=lines[start+1].split();values=np.loadtxt(lines[start+2:])
 if values.ndim!=2 or values.shape[1]!=len(names) or not np.isfinite(values).all() or not np.all(np.diff(values[:,0])>0):raise ValueError('Malformed/nonfinite/native time trajectory')
 if len(set(names))!=len(names):raise ValueError('Duplicate native columns')
 return names,values,meta

def main():
 p=argparse.ArgumentParser();p.add_argument('artifact',type=Path);args=p.parse_args();out=args.artifact.resolve();execution=json.loads((out/'execution.json').read_text());manifest=json.loads((out/'manifest.json').read_text());checks={};metrics={}
 checks['native_process_success']=execution['native_success']
 checks['retained_output_hashes_match']=bool(execution.get('outputs')) and all((out/name).is_file() and hashlib.sha256((out/name).read_bytes()).hexdigest()==digest for name,digest in execution.get('outputs',{}).items())
 checks['retained_executable_hash_matches']=hashlib.sha256((out/'moco3d').read_bytes()).hexdigest()==manifest['executable_sha256']
 model_path=out/'optimized_plant.osim'
 checks['optimized_plant_present']=model_path.exists()
 if model_path.exists():
  model=ET.parse(model_path)
  mass=sum(float(b.findtext('mass')) for b in model.iter('Body'))
  contacts=list(model.iter('SmoothSphereHalfSpaceForce'))
  checks['source_mass_and_3d_contacts']=abs(mass-85.26984854173146)<1e-8 and len(contacts)==12 and len(list(model.iter('Body')))==22
  checks['no_measured_external_loads']=not list(model.iter('ExternalLoads'))
  checks['no_prescribed_coordinates']=all(c.findtext('prescribed','false')=='false' for c in model.iter('Coordinate'))
  checks['no_pelvis_actuators']=not any('pelvis_' in a.findtext('coordinate','') for a in model.iter('CoordinateActuator'))
  metrics['source_mass_kg']=mass
  metrics['actuator_inventory']=[{'name':a.attrib.get('name'),'coordinate':a.findtext('coordinate'),'optimal_force':a.findtext('optimal_force')} for a in model.iter('CoordinateActuator')]
 stage='muscle_driven_tracking' if manifest['mode']=='muscle' else 'torque_driven_tracking'
 solution=out/f'example3DWalking_{stage}_solution.sto';forward=out/'forward_validation.sto'
 checks['solution_and_forward_present']=solution.exists() and forward.exists()
 if checks['solution_and_forward_present']:
  names,a,meta=read_sto(solution);other,b,_=read_sto(forward)
  checks['solver_declared_success']=meta.get('success','').lower() in ['true','1']
  checks['source_interval']=abs(a[0,0]-.48)<1e-9 and abs(a[-1,0]-1.61)<1e-9 and abs(b[0,0]-.48)<1e-9 and abs(b[-1,0]-1.61)<.0011
  translations=[];angles=[];errors={}
  for j,key in enumerate(names):
   if not key.endswith('/value'):continue
   if key not in other:raise ValueError('Missing forward coordinate '+key)
   delta=np.interp(b[:,0],a[:,0],a[:,j])-b[:,other.index(key)];error=float(max(abs(delta)));errors[key]=error
   (translations if any(x in key for x in ['pelvis_tx','pelvis_ty','pelvis_tz']) else angles).append(error)
  metrics['coordinate_maximum_absolute_error']=errors
  checks['translation_forward_tolerance_5cm']=bool(translations) and max(translations)<.05
  checks['angular_forward_tolerance_point1rad']=bool(angles) and max(angles)<.1
  grf=out/f'example3DWalking_{stage}_ground_reactions.sto';forward_grf=out/'forward_ground_reactions.sto'
  checks['ground_force_outputs']=grf.exists() and forward_grf.exists()
  if checks['ground_force_outputs']:
   gn,ga,_=read_sto(grf);fn,fa,_=read_sto(forward_grf)
   vertical=[j for j,n in enumerate(gn) if n.endswith('_vy')];checks['nonzero_vertical_support']=bool(vertical) and np.max(ga[:,vertical])>1
   metrics['ground_force_columns']=gn;metrics['peak_vertical_ground_force_n']=float(np.max(ga[:,vertical])) if vertical else None
  muscles=[j for j,n in enumerate(names) if n.endswith('/activation')]
  checks['expected_muscle_states']=len(muscles)==(80 if manifest['mode']=='muscle' else 0)
  if muscles:
   checks['activation_bounds']=np.min(a[:,muscles])>=-1e-6 and np.max(a[:,muscles])<=1+1e-6
   forward_muscles=[other.index(names[j]) for j in muscles]
   checks['forward_activation_bounds']=np.min(b[:,forward_muscles])>=-1e-6 and np.max(b[:,forward_muscles])<=1+1e-6
  control_columns=[j for j,n in enumerate(names) if n.startswith('/forceset/') and n.count('/')==2]
  checks['expected_control_count']=len(control_columns)==(93 if manifest['mode']=='muscle' else 25)
  if manifest['mode']=='muscle':
   muscle_control_indices=[names.index('/forceset/'+m.attrib['name']) for m in model.iter('DeGrooteFregly2016Muscle')]
   checks['muscle_excitation_bounds']=len(muscle_control_indices)==80 and np.min(a[:,muscle_control_indices])>=-1e-6 and np.max(a[:,muscle_control_indices])<=1+1e-6
  checks['control_bounds']=bool(control_columns) and np.min(a[:,control_columns])>=-1-1e-6 and np.max(a[:,control_columns])<=1+1e-6
  forcepath=out/'optimized_muscle_and_actuator_outputs.sto';checks['physical_actuation_telemetry']=forcepath.exists()
  if forcepath.exists():
   on,oa,_=read_sto(forcepath);metrics['actuation_output_columns']=len(on)-1
   for suffix,label,lower in [('fiber_length','fiber_length_positive',0),('tendon_force','tendon_force_nonnegative',-1e-6)]:
    idx=[j for j,n in enumerate(on) if n.endswith(suffix)]
    if manifest['mode']=='muscle':checks[label]=len(idx)==80 and np.min(oa[:,idx])>=lower
 checks={k:bool(v) for k,v in checks.items()}
 report={'passed':all(checks.values()),'checks':checks,'metrics':metrics,'manifest_sha256':hashlib.sha256((out/'manifest.json').read_bytes()).hexdigest(),'execution_sha256':hashlib.sha256((out/'execution.json').read_bytes()).hexdigest(),'limitation':'Open-loop source tracking/forward validation, not neural walking or canonical subject registration.'}
 (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
 if not report['passed']:raise SystemExit(1)
if __name__=='__main__':main()
