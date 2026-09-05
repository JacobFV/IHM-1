"""Exact supported default-pose tests and real-source closure checks."""
import sys, unittest, xml.etree.ElementTree as E
import json,gzip,hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.spatial.opensim import offset_transform, evaluate_function, solve_frames, load_model
class PoseTests(unittest.TestCase):
    def test_rotation_and_chain(self):
        t=offset_transform([1,2,3],[0,0,np.pi/2]);np.testing.assert_allclose(t@np.array([1,0,0,1]),[1,3,3,1],atol=1e-12)
        joints=[('ground','a',offset_transform([1,0,0],[0,0,np.pi/2]),np.eye(4),np.eye(4)),('a','b',offset_transform([0,2,0],[0,0,0]),np.eye(4),offset_transform([0,1,0],[0,0,0]))]
        f=solve_frames(['a','b'],joints);np.testing.assert_allclose(f['b'][:3,3],[0,0,0],atol=1e-12)
        with self.assertRaises(ValueError):solve_frames(['a'],[])
    def test_functions(self):
        self.assertEqual(evaluate_function(E.fromstring('<LinearFunction><coefficients>2 3</coefficients></LinearFunction>'),[4]),11)
        f=E.fromstring('<SimmSpline><x>-1 0 1</x><y>2 3 4</y></SimmSpline>')
        self.assertEqual(evaluate_function(f,[0]),3)
        with self.assertRaises(ValueError):evaluate_function(f,[.1])
    def test_actual_pose(self):
        model=load_model(Path('data/raw/anatomy/opensim-models/source/Models/Rajagopal/Rajagopal2016.osim'))
        self.assertEqual(len(model['frames']),23)
        self.assertEqual(len(model['muscles']),80)
        self.assertEqual(sum(len(m['points']) for m in model['muscles']),288)
        for t in model['frames'].values():
            np.testing.assert_allclose(t[:3,:3].T@t[:3,:3],np.eye(3),atol=1e-12)
        np.testing.assert_allclose(model['frames']['pelvis'][:3,3],[0,.94,0])
class DisplayTests(unittest.TestCase):
    def test_assets(self):
        path=Path('data/derived/app/manifest.json')
        if not path.exists():self.skipTest('Build app manifest first')
        manifest=json.loads(path.read_text());models=[x for x in manifest['models'] if x['id']=='opensim-rajagopal']
        if not models:self.skipTest('Build OpenSim display first')
        self.assertEqual(len(models),1)
        rows=[x for x in manifest['structures'] if x['model_id']=='opensim-rajagopal']
        self.assertEqual(len(rows),161);self.assertEqual(len(set(x['id'] for x in rows)),161)
        for row in rows:
            p=path.parent/'geometry'/(row['id']+'.json.gz')
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),row['geometry_sha256'])
            data=json.loads(gzip.decompress(p.read_bytes()));points=np.asarray(data['positions']).reshape(-1,3)
            self.assertTrue(np.isfinite(points).all());self.assertEqual(data['units'],'m')
            if row['kind']=='mesh':
                indices=np.asarray(data['indices']);self.assertGreaterEqual(indices.min(),0);self.assertLess(indices.max(),len(points))
            else:
                self.assertEqual(len(points)%2,0)
                if data['wrap_solved']:
                    self.assertGreater(data['native_mechanics']['length_m'],0)
                    self.assertEqual(data['native_mechanics']['engine_variant'],'wrap_8_0.0005_cache')
if __name__=='__main__':unittest.main()
