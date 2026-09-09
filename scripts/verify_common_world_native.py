"""Actual native body, loose cover, and ball in a common gravity/support frame.
Run with PYTHONPATH=. .venv/bin/python scripts/verify_common_world_native.py.
Retains 0.1 s mechanics evidence; does not instantiate brain or physiology.
"""
import json,tempfile,time
from pathlib import Path
import numpy as np
from ihm.assembly.articulated import ArticulatedBodyPlant
from ihm.assembly.environment_dynamics import EnvironmentDynamics
root=Path(__file__).resolve().parents[1];output=Path(tempfile.mkdtemp(prefix='common-world-native-',dir=root/'data/derived'));plant=None;start=time.monotonic();rows=[];report={}
try:
 plant=ArticulatedBodyPlant(root,output/'plant',environment='supine',target_mass_kg=77.6122029,augmented_registration='data/derived/mechanics/whole_body_lumbar_current/registration.json')
 state=plant.snapshot();owner=EnvironmentDynamics(root,'supine',{'scene':'bedroom','objects':['ball-small']},plant.registration);bridge=owner.world_frame
 expected_g=plant.registration.basis@np.asarray(plant.registration.reference['gravity_m_s2']);np.testing.assert_allclose(bridge.canonical_gravity_m_s2,expected_g,atol=1e-12)
 plane=state['body_environment']['plane'];world_plane=bridge.points_to_canonical([0,0,owner.base_plane]);normal=np.asarray(plane['normal']);plane_error=float(abs((world_plane-plane['point_m'])@normal));assert plane_error<1e-10
 np.testing.assert_allclose(bridge.vectors_to_canonical([0,0,1]),normal,atol=1e-12)
 ball=next(p for p in owner.rigid if p.radius);ball.x[:]=[.9,0,owner.object_plane+ball.radius+.01];cloth=next(p for p in owner.soft if p.kind=='cloth')
 for tick in range(20):
  raw=[{'id':ball.id,'force_n':bridge.vectors_to_canonical([1,0,0]).tolist(),'point_m':bridge.points_to_canonical(ball.x).tolist()}]
  bound=owner.validate_object_forces(raw);np.testing.assert_allclose(bound[0]['force_n'],[1,0,0],atol=1e-12)
  ports=owner.advance(.005,state['entities'],bound);np.testing.assert_allclose(np.sum([p['force_n'] for p in ports],axis=0)*.005,owner.last_impulse,atol=1e-9)
  state=plant.advance(.005,forces=ports)
  if tick%4==3:rows.append({'t':state['time_s'],'ball_world_z_m':float(ball.x[2]),'ball_world_vz_m_s':float(ball.v[2]),'ball_canonical_position_m':bridge.points_to_canonical(ball.x).tolist(),'cloth_max_extension':cloth.stretch_state.get('max_extension_ratio'),'cloth_constraint_converged':cloth.stretch_state.get('converged'),'cloth_max_speed_m_s':float(np.linalg.norm(cloth.v,axis=1).max()),'body_contact_count':len(ports)})
 assert any(r['ball_world_vz_m_s']>0 for r in rows)
 report={'passed':True,'plane_error_m':plane_error,'canonical_gravity_m_s2':expected_g.tolist(),'world_frame':bridge.metadata(),'rows':rows,'wall_s':time.monotonic()-start,'scope':'Actual 98-muscle native body, loose cloth, and ball with shared gravity and bed plane; 5 ms exchange, 0.1 s; no brain/physiology or calibrated equilibrium claim'}
except BaseException as e:
 report={'passed':False,'error':str(e),'rows':rows,'wall_s':time.monotonic()-start};raise
finally:
 (output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps({'output':str(output),**report},indent=2),flush=True)
 if plant:plant.close()
