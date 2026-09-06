"""Light materialization checks; --run-native requires the shared heavy slot."""
from pathlib import Path
import json,sys,tempfile,unittest
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from materialize_lumbar_muscle_variant import materialize
class MaterializationTests(unittest.TestCase):
    def test_preserves_model_and_appends_exact_six(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'data/derived') as d:
            r=materialize(ROOT,Path(d)); model=(ROOT/r['model_path']).read_bytes()
            old=(ROOT/r['base_model_path']).read_bytes(); insert=(ROOT/r['insert_path']).read_bytes()
            self.assertEqual(model.replace(insert,b'',1),old)
            root=ET.fromstring(model);self.assertEqual(len(root.findall('.//BodySet/objects/Body')),22)
            muscles=[m for m in root.findall('.//ForceSet/objects/*') if 'Muscle' in m.tag]
            cat=json.loads((ROOT/r['catalog_path']).read_text())
            self.assertEqual([m.get('name') for m in muscles],[m['id'] for m in cat])
            self.assertEqual(len(cat),98)
            self.assertFalse(r['default_enabled']);self.assertFalse(r['native_acceptance_complete'])
            with self.assertRaises(ValueError):materialize(ROOT,Path(d))
if __name__=='__main__':unittest.main()
