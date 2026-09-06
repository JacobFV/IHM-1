import tempfile,unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from materialize_arm26_wrap_probe import ROOT,materialize
class Tests(unittest.TestCase):
    def test_exact_scope_and_overwrite(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'data/derived') as d:
            r=materialize(Path(d));self.assertEqual(len(r['changed_wraps']),4)
            self.assertNotEqual(r['model_sha256']['observed'],r['model_sha256']['candidate'])
            candidate=ET.parse(Path(d)/'candidate.osim').getroot();source=ET.parse(ROOT/r['source']).getroot()
            a=source.findall('.//ForceSet/objects/*');b=candidate.findall('.//ForceSet/objects/*');self.assertEqual(len(a),len(b))
            for old,new in zip(a,b):
                if old.get('name','').startswith('arm26_BRA_'):
                    refs=new.findall('.//wrap_object');self.assertEqual(len(refs),1);self.assertTrue(refs[0].text.endswith('__bra_exact'))
                else:self.assertEqual(ET.tostring(old),ET.tostring(new))
            self.assertEqual(len(candidate.findall('.//ExactArm26Cylinder')),2)
            with self.assertRaises(ValueError):materialize(Path(d))
if __name__=='__main__':unittest.main()
