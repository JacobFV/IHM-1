"""Forward native-only stability of an explicitly initialized supported body."""
from pathlib import Path
import argparse,json,sys,tempfile,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream

def run(registration,pose_path,duration):
 out=Path(tempfile.mkdtemp(prefix='supine-native-equilibrium-',dir=ROOT/'data/derived'));artifact=json.loads(pose_path.read_text());start=time.monotonic();rows=[];native=None
 report={'passed':False,'registration':registration,'pose_path':str(pose_path),'duration_s':duration,'scope':'Native MM supported body only, no cloth, brain or physiology; held explicit initialized muscle excitation defaults','thresholds':{'maximum_coordinate_change':1e-5,'maximum_kinetic_energy_j':1e-6,'maximum_momentum_residual_n':1e-5}}
 (out/'verification_source.py').write_bytes(Path(__file__).read_bytes())
 try:
  native=NativeMechanicalStream(ROOT,out/'native',target_mass_kg=artifact['target_mass_kg'],environment='supine',augmented_registration=registration,initial_pose=artifact['coordinates'],surface_contact_manifest='data/derived/supine-surface-contact-5jqy1juo/manifest.json',bed_material='MM')
  first=native.snapshot();report['initial_state']=first
  for _ in range(round(duration/.02)):
   state=native.advance(.02)
   rows.append({'time_s':state['time_s'],'maximum_coordinate_change':max(abs(v['value']-first['coordinates'][k]['value']) for k,v in state['coordinates'].items()),'kinetic_energy_j':state['kinetic_energy_j'],'momentum_residual_n':sum(x*x for x in state['momentum_balance_residual_n'])**.5})
  report['maximum_coordinate_change']=max(x['maximum_coordinate_change'] for x in rows);report['maximum_kinetic_energy_j']=max(x['kinetic_energy_j'] for x in rows);report['maximum_momentum_residual_n']=max(x['momentum_residual_n'] for x in rows)
  report['final_state']=state;report['passed']=all(report[k]<=v for k,v in report['thresholds'].items())
 except Exception as error:report['error']=type(error).__name__+': '+str(error)
 finally:
  if native is not None:native.close()
  report.update(rows=rows,wall_s=time.monotonic()-start);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('rows','initial_state','final_state')}|{'output':str(out)},indent=2))
 return report
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('registration');parser.add_argument('pose',type=Path);parser.add_argument('--duration',type=float,default=1);args=parser.parse_args();raise SystemExit(0 if run(args.registration,args.pose,args.duration)['passed'] else 1)
