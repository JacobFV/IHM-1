import tempfile,unittest,xml.etree.ElementTree as ET
from pathlib import Path
from materialize_stationary_wrap_probe import ROOT,materialize
class Tests(unittest.TestCase):
    def test_private_four_and_exact_force_parameters(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'data/derived') as d:
            r=materialize(Path(d));original=ET.parse(ROOT/r['source']).getroot();candidate=ET.parse(Path(d)/'candidate.osim').getroot()
            self.assertEqual(len(candidate.findall('.//ExactArm26Cylinder')),2);self.assertEqual(len(candidate.findall('.//StationaryArm26Ellipsoid')),2)
            for old,new in zip(original.findall('.//ForceSet/objects/*'),candidate.findall('.//ForceSet/objects/*')):
                if old.get('name','').startswith(('arm26_BRA_','arm26_BIClong_')):
                    for a,b in zip(old.findall('.//wrap_object'),new.findall('.//wrap_object')):b.text=a.text
                self.assertEqual(ET.tostring(old),ET.tostring(new))
            for tag in ('BodySet','JointSet'):
                # BodySet adds only massless wrap objects; physical body parameters stay exact.
                if tag=='JointSet':self.assertEqual(ET.tostring(original.find('.//'+tag)),ET.tostring(candidate.find('.//'+tag)))
            for a,b in zip(original.findall('.//BodySet/objects/Body'),candidate.findall('.//BodySet/objects/Body')):
                for field in ('mass','mass_center','inertia'):self.assertEqual(a.findtext(field),b.findtext(field))
if __name__=='__main__':unittest.main()
