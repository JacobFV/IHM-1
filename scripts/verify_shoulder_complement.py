#!/usr/bin/env python3
"""Source-only shoulder dependency closure; no OpenSim or native initialization."""
from pathlib import Path
import sys,unittest,xml.etree.ElementTree as E
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from materialize_shoulder_complement import generate,SELECTED
class Checks(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.report,cls.forces,cls.dependencies=generate(ROOT)
 def test_source_parameters_and_frames(self):
  r=self.report;self.assertEqual(len(r['muscles']),15);self.assertEqual(r['donor']['muscle_count'],50)
  self.assertEqual(r['muscles'][0]['parameters']['max_isometric_force']['value'],1218.9)
  self.assertEqual(r['muscles'][0]['parameters']['optimal_fiber_length']['unit'],'m')
  self.assertEqual(set(r['attachment_body_names']),{'humerus','scapula','clavicle','thorax'})
  self.assertFalse(r['native_installation']);self.assertFalse(r['mirror_equals_official_release_verified'])
 def test_xml_exact_subtree_and_dependency_closure(self):
  src=E.parse(ROOT/self.report['donor']['path']).getroot();dst=E.fromstring(self.forces)
  self.assertEqual([x.get('name') for x in dst.findall('./ForceSet/objects/*')],list(SELECTED))
  for name in SELECTED:
   a=src.find(f".//Millard2012EquilibriumMuscle[@name='{name}']");b=dst.find(f".//Millard2012EquilibriumMuscle[@name='{name}']")
   self.assertEqual(E.tostring(a),E.tostring(b))
  dep=E.fromstring(self.dependencies);wraps={x.get('name') for x in dep.findall('.//WrapObjectSet/objects/*')}
  self.assertTrue(set(self.report['referenced_wrap_names'])<=wraps)
  self.assertEqual(len(dep.findall('.//CoordinateCouplerConstraint')),13)
  self.assertEqual(len(dep.findall('.//Body')),12)
  self.assertEqual(len(dep.findall('.//CustomJoint')),12)
 def test_no_silent_coordinate_or_route_substitution(self):
  r=self.report;self.assertGreater(r['path_point_type_counts']['MovingPathPoint'],0);self.assertGreater(r['path_point_type_counts']['ConditionalPathPoint'],0)
  self.assertEqual(r['registration']['body_transforms'],None);self.assertEqual(r['registration']['coordinate_mapping'],None)
  self.assertEqual(r['registration']['added_native_mass_kg'],0)
  self.assertEqual(r['unrepresented_girdle_muscles'],['trapezius','serratus anterior','rhomboids','levator scapulae','pectoralis minor'])
 def test_reproducible(self):self.assertEqual(generate(ROOT),(self.report,self.forces,self.dependencies))
if __name__=='__main__':unittest.main()
