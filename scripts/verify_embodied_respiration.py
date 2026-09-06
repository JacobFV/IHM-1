"""Tiny geometry/virtual-work tests; no native runtime or browser."""
from pathlib import Path
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.embodied_respiration import EmbodiedRespiration

def fixture():
 return {'volume_jacobian_m2':[.1,.2,.3],'stiffness_n_m':np.diag([100,200,300]).tolist(),'reference_dimensions_m':{'transverse':.3,'anteroposterior':.2,'height':.2},'bindings':[{'entity_id':n,'kind':kind,'centroid_m':[0,0,0],'translation_basis':np.eye(3).tolist()} for n,kind in [('rib','rib'),('lung','lung'),('diaphragm','diaphragm')]],'skin_field':{'entity_ids':['skin'],'center_m':[0,0,0],'bounds_m':{'min':[-.3,-.3,-.2],'max':[.3,.3,.2]},'thorax_y_m':[-.1,.1],'reference_radii_m':[.1,.1],'basis':'thoracic-smoothstep-v1'}}
def entities(rotation=np.eye(3)):
 return {n:{'centroid_m':[.2,.1,-.3],'translation_m':[.2,.1,-.3],'rotation_matrix':rotation.tolist(),'deformation_gradient':np.eye(3).tolist()} for n in ['rib','lung','diaphragm','skin']}
class Tests(unittest.TestCase):
 def test_held_canonical_bindings_include_posterior_supports(self):
  import json
  payload=json.loads((ROOT/'data/derived/canonical/respiration.json').read_text());e={b['entity_id']:{'centroid_m':b['centroid_m'],'translation_m':[0,0,0]} for b in payload['bindings']}
  for ident in payload['skin_field']['entity_ids']:e[ident]={'centroid_m':[0,0,0],'translation_m':[0,0,0]}
  g=EmbodiedRespiration(payload,3000).geometry(3300,e)
  for binding in payload['bindings']:
   if binding['kind']=='lung':self.assertAlmostEqual(np.linalg.det(g['entities'][binding['entity_id']]['deformation_gradient']),1.1)
   if binding['kind']=='posterior_support':self.assertEqual(g['entities'][binding['entity_id']]['translation_m'],[0,0,0])
 def test_volume_and_material_ownership(self):
  p=EmbodiedRespiration(fixture(),3000);g=p.geometry(3300,entities(),time_s=1)
  self.assertAlmostEqual(np.linalg.det(g['entities']['lung']['deformation_gradient']),1.1)
  self.assertAlmostEqual(np.linalg.det(g['entities']['diaphragm']['deformation_gradient']),1.)
  self.assertFalse(g['independent_mass_or_recoil']);self.assertAlmostEqual(np.dot([.1,.2,.3],g['displacement_m']),.0003)
 def test_rotated_rib_force_conjugacy(self):
  r=np.array([[0,-1,0],[1,0,0],[0,0,1]]);p=EmbodiedRespiration(fixture(),3000);f=np.array([2,-3,4]);x=np.array([.02,.03,.01]);e=entities(r);point=p.point_position('rib',x,3200,e);load=p.project_load([{'id':'rib','point_m':point.tolist(),'force_n':f.tolist()}],e,3200)
  self.assertTrue(np.allclose(load['force_ports'][0]['point_jacobian_m_per_m3'],r@p.b));self.assertAlmostEqual(load['external_pressure_pa'],-f@(r@p.b))
 def test_offcenter_shape_and_skin_derivative(self):
  p=EmbodiedRespiration(fixture(),3000);r=np.array([[1,0,0],[0,0,-1],[0,1,0]]);e=entities(r);f=np.array([3,-2,5]);x=np.array([.04,.02,.03]);volume=3200.;h=.001
  for name in ['lung','diaphragm','skin']:
   point=p.point_position(name,x,volume,e);load=p.project_load([{'id':name,'point_m':point.tolist(),'force_n':f.tolist()}],e,volume);j=np.array(load['force_ports'][0]['point_jacobian_m_per_m3']);fd=(p.point_position(name,x,volume+h,e)-p.point_position(name,x,volume-h,e))/(2*h*1e-6)
   self.assertTrue(np.allclose(j,fd,rtol=1e-6,atol=1e-7),(name,j,fd));self.assertAlmostEqual(-load['external_pressure_pa'],f@j)
 def test_resultant_wrench_scope(self):
  p=EmbodiedRespiration(fixture(),3000);e=entities();f={'id':'rib','point_m':[.2,.1,-.3],'force_n':[1,0,0],'moment_nm':[0,2,0]}
  load=p.project_load([f],e,3100);self.assertEqual(load['force_ports'][0]['angular_jacobian_rad_per_m3'],[0,0,0])
  f['id']='lung'
  with self.assertRaises(ValueError):p.project_load([f],e,3100)
 def test_no_double_articulation_and_invalid_volume(self):
  p=EmbodiedRespiration(fixture(),3000);e=entities();g=p.geometry(3000,e,time_s=0);self.assertTrue(np.allclose(g['entities']['rib']['translation_m'],e['rib']['translation_m']))
  with self.assertRaises(ValueError):p.geometry(-1,e,time_s=0)
if __name__=='__main__':unittest.main()
