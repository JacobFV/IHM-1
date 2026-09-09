"""Native/world step-refinement experiment from identical checkpoints.

Reports sensitivity, not a convergence certification or physiology validation.
Run from the IHM root with PYTHONPATH=. .venv/bin/python.
"""
import json,tempfile,time
from pathlib import Path
import numpy as np
from ihm.assembly.articulated import ArticulatedBodyPlant
from ihm.assembly.environment_dynamics import EnvironmentDynamics
root=Path(__file__).resolve().parents[1];output=Path(tempfile.mkdtemp(prefix='native-environment-dt-',dir=root/'data/derived'));plant=None;checkpoint=None;start=time.monotonic();results=[]
try:
 plant=ArticulatedBodyPlant(root,output/'plant',environment='supine');state=plant.snapshot();owner=EnvironmentDynamics(root,'supine',{'scene':'bedroom','objects':[]},plant.registration)
 checkpoint=plant.checkpoint();env_checkpoint=owner.checkpoint();initial_skin=owner.skin_points(state['entities']);cloth=next(m for m in owner.soft if m.kind=='cloth');initial_cloth=cloth.x.copy()
 for dt in (.02,.005):
  plant.restore(checkpoint);owner.restore(env_checkpoint);state=plant.snapshot();cloth=next(m for m in owner.soft if m.kind=='cloth');rows=[]
  for _ in range(round(.4/dt)):
   before=owner.skin_points(state['entities']);ports=owner.advance(dt,state['entities']);state=plant.advance(dt,forces=ports);skin=owner.skin_points(state['entities'])
   rows.append({'t':state['time_s'],'cloth_speed_max_m_s':float(np.max(np.linalg.norm(cloth.v,axis=1))),'cloth_displacement_max_m':float(np.max(np.linalg.norm(cloth.x-initial_cloth,axis=1))),'skin_speed_max_m_s':float(np.max(np.linalg.norm(skin-before,axis=1)))/dt,'skin_displacement_max_m':float(np.max(np.linalg.norm(skin-initial_skin,axis=1))),'native_external_work_j':state['audit']['external_work_j'],'force_max_n':max([np.linalg.norm(p['force_n']) for p in ports],default=0.)})
  r={'dt_s':dt,'duration_s':state['time_s'],'peak_cloth_speed_m_s':max(r['cloth_speed_max_m_s'] for r in rows),'peak_skin_speed_m_s':max(r['skin_speed_max_m_s'] for r in rows),'final':rows[-1],'rows':rows};results.append(r);print(json.dumps({k:v for k,v in r.items() if k!='rows'}),flush=True)
 report={'passed':True,'results':results,'wall_s':time.monotonic()-start,'scope':'Same native and environment initial checkpoints, .4s passive raw mechanics, varying body/environment exchange20ms vs5ms; no physiology/brain'}
except BaseException as e:
 report={'passed':False,'error':str(e),'results':results};raise
finally:
 (output/'report.json').write_text(json.dumps(report,indent=2));print(str(output),flush=True)
 if plant:
  if checkpoint:plant.release(checkpoint)
  plant.close()
