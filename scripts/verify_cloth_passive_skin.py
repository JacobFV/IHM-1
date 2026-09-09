"""Experimental passive cloth on held canonical skin in the native world frame.

No production solver changes, brain, physiology, or live native body. Uses the
retained proper gravity transform from the actual common-frame integration proof.
"""
import argparse,gzip,hashlib,json,tempfile,time
from pathlib import Path
import numpy as np
from ihm.assembly.environment_dynamics import SpringMesh
from ihm.assembly.cloth_passive_step import ClothPassiveStepper
from ihm.assembly.cloth_passive_geometry import ClothPassiveGeometry
ROOT=Path(__file__).resolve().parents[1]

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--seconds',type=float,default=.5);args=ap.parse_args()
 output=Path(tempfile.mkdtemp(prefix='passive-cloth-skin-',dir=ROOT/'data/derived'));start=time.monotonic();rows=[]
 geometry_path=ROOT/'data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz'
 points=np.asarray(json.loads(gzip.decompress(geometry_path.read_bytes()))['positions']).reshape(-1,3)
 _,indices=np.unique(np.floor(points/.045).astype(int),axis=0,return_index=True);points=points[indices]
 frame_path=ROOT/'data/derived/common-world-native-907cp3gk/report.json';frame=json.loads(frame_path.read_text())['world_frame'];matrix=np.asarray(frame['canonical_to_world'])
 points=points@matrix[:3,:3].T+matrix[:3,3]
 mesh=SpringMesh('blanket-experimental','cloth',[-.42,-.88,0],[.42,-.04,0],1.5)
 mesh.x[:,2]=points[:,2].max()+.052;mesh.rest=mesh.x.copy()
 stepper=ClothPassiveStepper(mesh);solver=ClothPassiveGeometry(stepper)
 spheres=[{'center_m':point.tolist(),'radius_m':.048,'owner':i} for i,point in enumerate(points)]
 sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [geometry_path,frame_path,Path(__file__).resolve(),ROOT/'ihm/assembly/cloth_passive_step.py',ROOT/'ihm/assembly/cloth_passive_geometry.py']}
 for name,expected in sources.items():
  source=ROOT/name;target=output/'inputs'/name;target.parent.mkdir(parents=True,exist_ok=True);raw=source.read_bytes();assert hashlib.sha256(raw).hexdigest()==expected;target.write_bytes(raw)
 report={'passed':False,'scope':'Experimental437-node cloth, actual sampled canonical skin held fixed in native-aligned frame; no native dynamics, brain or physiology','source_sha256':sources}
 try:
  for i in range(round(args.seconds/.005)):
   result=solver.advance(mesh,.005,spheres=spheres,floor_height_m=-.86,boxes=[{'min_m':[-.6,-.95,-.44],'max_m':[.6,.95,-.24],'owner':'mattress-top'}])
   assert result['energy_excess_j']<=1e-8
   assert np.linalg.norm(result['momentum_residual_ns'])<=1e-6, 'Unbalanced cloth/support momentum'
   assert result['stationarity_residual']<=1e-9, 'Implicit force balance did not converge'
   rows.append({'time_s':(i+1)*.005,'max_speed_m_s':float(np.linalg.norm(mesh.v,axis=1).max()),'diagnostics':result})
   if (i+1)%20==0:print(json.dumps({'output':str(output),'time_s':rows[-1]['time_s'],'speed':rows[-1]['max_speed_m_s']}),flush=True)
  report.update(passed=True)
 except Exception as error:
  report.update(error=f'{type(error).__name__}: {error}',failure_diagnostics=getattr(error,'diagnostics',None))
 finally:
  report.update(accepted_seconds=len(rows)*.005,wall_s=time.monotonic()-start,rows=rows,final_positions=mesh.x.tolist())
  (output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k not in ('rows','final_positions')},indent=2),flush=True);print(str(output),flush=True)
 if not report['passed']:raise SystemExit(1)
if __name__=='__main__':main()
