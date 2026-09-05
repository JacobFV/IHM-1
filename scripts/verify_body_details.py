from pathlib import Path
import sys,unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
class DetailTests(unittest.TestCase):
 def test_flow_conserves_at_branch(self):
  from ihm.assembly.details import solve_network
  p=np.array([[0,0,0],[1,0,0],[2,1,0],[2,-1,0]],float);e=np.array([[0,1],[1,2],[1,3]])
  state=solve_network(p,e,np.array([.01,.01,.01]),{0:100.,2:0.,3:0.})
  self.assertLess(state['maximum_internal_residual_m3_s'],1e-15)
  self.assertAlmostEqual(state['flow_m3_s'][0],sum(state['flow_m3_s'][1:]),places=15)
  self.assertGreater(state['dissipation_w'],0)
 def test_paired_tree_has_return_and_capillaries(self):
  from ihm.assembly.details import microvascular_unit,solve_network
  p,e,r,k=microvascular_unit(np.zeros(3),.01,32,seed=9)
  self.assertEqual(sum(k==2),32);self.assertTrue(np.all(r>0))
  self.assertLessEqual(np.max(np.linalg.norm(p,axis=1)),.01)
  s=solve_network(p,e,r,{0:100.,1:0.})
  self.assertLess(s['maximum_internal_residual_m3_s'],1e-17)
  self.assertTrue(np.all(s['flow_m3_s']>0))
 def test_skin_samples_keep_surface_attachment(self):
  from ihm.assembly.details import sample_hair
  vertices=np.array([[0,0,0],[.01,0,0],[0,.01,0]],float);faces=np.array([[0,1,2]])
  h=sample_hair(vertices,faces,np.array([100.]),seed=1)
  self.assertGreater(len(h['roots_m']),20)
  self.assertLess(len(h['roots_m']),80)
  self.assertTrue(np.allclose(h['roots_m'][:,2],0))
  self.assertTrue(np.all(h['roots_m'][:,:2].sum(1)<=.01))
  self.assertTrue(np.allclose(h['barycentric'].sum(1),1))
  self.assertTrue(np.all(h['radius_m']>0))
 def test_welded_surface_ignores_attribute_seams(self):
  from ihm.assembly.details import physical_mesh
  v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float);f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
  split=v[f].reshape(-1,3);m=physical_mesh(split,np.arange(12).reshape(-1,3))
  self.assertTrue(m.is_watertight);self.assertAlmostEqual(m.volume,1/6)
class InputIntegrityTests(unittest.TestCase):
 def test_fractional_and_invalid_mesh_indices_are_rejected(self):
  from ihm.assembly.details import physical_mesh,sample_hair
  v=np.eye(3)
  for f in [[[0.5,1,2]],[[-1,1,2]],[[0,1,3]],[[0,float('nan'),2]],[[0,1]]]:
   with self.subTest(faces=f):
    with self.assertRaises(ValueError):physical_mesh(v,f)
    with self.assertRaises(ValueError):sample_hair(v,f,[10])
 def test_pinned_skin_hash_frame_and_units(self):
  import tempfile,gzip,json,hashlib
  from unittest.mock import patch
  from scripts import build_body_details as builder
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);path=root/'skin.json.gz';path.write_bytes(gzip.compress(json.dumps({'positions':[],'indices':[]}).encode()))
   ref={'path':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'frame':'bodyparts3d-display-m','units':'m'}
   with patch.object(builder,'ROOT',root):
    self.assertEqual(builder.load_skin_geometry({'reference_geometry':ref})[0],path)
    for key,value in [('sha256','0'*64),('frame','another-frame'),('units','mm')]:
     with self.subTest(key=key),self.assertRaises(ValueError):builder.load_skin_geometry({'reference_geometry':{**ref,key:value}})
    path.write_bytes(gzip.compress(b'{"positions":[1],"indices":[]}'))
    with self.assertRaises(ValueError):builder.load_skin_geometry({'reference_geometry':ref})
 def test_append_publishes_current_geometry_before_manifest(self):
  import tempfile,json,hashlib
  from unittest.mock import patch
  from scripts import build_body_details as builder
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);app=root/'data/derived/app';app.mkdir(parents=True)
   manifest=app/'manifest.json';manifest.write_text(json.dumps({'structures':[{'id':'unrelated'}]}))
   source=root/'fresh.json.gz';source.write_bytes(b'new geometry bytes')
   s={'id':'detail','geometry_path':source.name,'geometry_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
   with patch.object(builder,'ROOT',root):
    builder.append_display([s]);published=app/'geometry/detail.json.gz'
    self.assertEqual(published.read_bytes(),source.read_bytes())
    published.write_bytes(b'stale');builder.append_display([s]);self.assertEqual(published.read_bytes(),source.read_bytes())
    self.assertEqual([e['id'] for e in json.loads(manifest.read_text())['structures']],['unrelated','detail'])
    before=manifest.read_bytes();source.write_bytes(b'corrupt')
    with self.assertRaises(ValueError):builder.append_display([s])
    self.assertEqual(manifest.read_bytes(),before)
if __name__=='__main__':unittest.main()
