"""Physical shell area must not count both sides of the source skin asset."""
import copy
import gzip
import hashlib
import json
import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.skin_layers import physical_skin_support

class Checks(unittest.TestCase):
    def setUp(self):
        self.raw=gzip.compress(json.dumps({'positions':[0,0,0,1,0,0,1,1,0,0,1,0,0,0,.01,1,0,.01,1,1,.01,0,1,.01], 'indices':[0,1,2,0,2,3,4,6,5,4,7,6]}).encode())
        self.ref={'path':'skin.json.gz','sha256':hashlib.sha256(self.raw).hexdigest(),'units':'m','frame':'bodyparts3d-display-m'}
        self.evidence={'schema':'engineered_skin_territories_v1','source_receipts':[{'entity_id':'skin','path':self.ref['path'],'sha256':self.ref['sha256']}], 'contact_eligible_triangle_ids':[0,1], 'surface_diagnostic':{'selected_component_id':0,'face_component_ids':[0,0,1,1],'selected_area_m2':1.,'selection_evidence':'inferred_largest_component_positive_signed_volume_integral','physical_surface_exclusivity_validated':False,'self_intersections_tested':False,'signed_integral_is_closed_volume':False}}
    def run_support(self,e=None,raw=None):
        return physical_skin_support(self.ref,self.raw if raw is None else raw,json.dumps(self.evidence if e is None else e).encode())
    def test_two_shells(self):
        support=self.run_support()
        self.assertEqual(support['area_m2'],1.)
        self.assertEqual(support['raw_source_area_m2'],2.)
        self.assertEqual(support['selected_triangle_count'],2)
        self.assertFalse(support['physical_surface_exclusivity_validated'])
        self.assertAlmostEqual(support['area_m2']*.0066,.0066)
    def test_reject_changed_source(self):
        with self.assertRaises(ValueError):self.run_support(raw=self.raw+b'changed')
    def test_reject_mismatched_mask_area_component(self):
        for key,value in [('contact_eligible_triangle_ids',[0,0]),('contact_eligible_triangle_ids',[0,2]),('contact_eligible_triangle_ids',[True,1])]:
            e=copy.deepcopy(self.evidence);e[key]=value
            with self.assertRaises(ValueError):self.run_support(e)
        e=copy.deepcopy(self.evidence);e['surface_diagnostic']['selected_area_m2']=2.
        with self.assertRaises(ValueError):self.run_support(e)
    def test_reject_nonselected_area_overflow(self):
        geometry=json.loads(gzip.decompress(self.raw))
        geometry['positions'][15]=1e200;geometry['positions'][18]=1e200
        geometry['positions'][19]=1e200;geometry['positions'][22]=1e200
        raw=gzip.compress(json.dumps(geometry).encode())
        self.ref['sha256']=hashlib.sha256(raw).hexdigest()
        self.evidence['source_receipts'][0]['sha256']=self.ref['sha256']
        with self.assertRaises(ValueError):self.run_support(raw=raw)
    def test_reject_stale_receipt(self):
        e=copy.deepcopy(self.evidence);e['source_receipts'][0]['sha256']='0'*64
        with self.assertRaises(ValueError):self.run_support(e)

if __name__=='__main__':unittest.main()
