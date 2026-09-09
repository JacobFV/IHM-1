"""Snapshot data ownership, graph structure, and public-boundary checks; no native process."""
from copy import deepcopy
import pickle
import unittest
from collections import OrderedDict
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.snapshot_data import clone_snapshot_data
from ihm.assembly.articulated import ArticulatedBodyPlant
from ihm.assembly.embodied_respiration import EmbodiedRespiration
from scripts.verify_embodied_respiration import fixture,entities
from scripts import verify_embodied_runtime as runtime_fixture

class Tests(unittest.TestCase):
    def test_graph_types_aliases_cycles_and_arrays_are_isolated(self):
        shared=[{'x':1}];array=np.arange(12,dtype=np.float32).reshape(3,4)
        source=OrderedDict(shared=shared,alias=shared,array=array,array_alias=array,tuple=(1,2),set={3,4},scalar=np.float32(.5))
        source['cycle']=source
        result=clone_snapshot_data(source)
        self.assertIsInstance(result,OrderedDict);self.assertIs(result['cycle'],result)
        self.assertIs(result['shared'],result['alias']);self.assertIs(result['array'],result['array_alias'])
        self.assertIsInstance(result['tuple'],tuple);self.assertIsInstance(result['set'],set)
        self.assertEqual(type(result['scalar']),type(source['scalar']))
        np.testing.assert_array_equal(result['array'],array);self.assertEqual(result['array'].dtype,array.dtype)
        self.assertFalse(np.shares_memory(result['array'],array))
        result['shared'][0]['x']=8;result['array'][0,0]=9
        self.assertEqual(source['shared'][0]['x'],1);self.assertEqual(source['array'][0,0],0)
        source['array'][1,1]=20;self.assertEqual(result['array'][1,1],5)

    def test_serialized_input_is_copied_as_bytes_not_deserialized(self):
        encoded=pickle.dumps({'would_be_decoded':True},protocol=5)
        self.assertEqual(clone_snapshot_data(encoded),encoded)
        self.assertIsInstance(clone_snapshot_data(encoded),bytes)

    def test_articulated_snapshots_checkpoints_and_restore_are_independent(self):
        class Native:
            def checkpoint(self):return {'native_tick':0}
            def restore(self,value):self.restored=deepcopy(value)
        plant=ArticulatedBodyPlant.__new__(ArticulatedBodyPlant)
        plant.state={'time_s':0.,'entities':{'bone':{'centroid_m':[1.,2.,3.]}}};plant.native=Native();plant.garments=None
        public=plant.snapshot();checkpoint=plant.checkpoint();public['entities']['bone']['centroid_m'][0]=100
        self.assertEqual(plant.state['entities']['bone']['centroid_m'][0],1.)
        plant.state['entities']['bone']['centroid_m'][0]=200;plant.restore(checkpoint)
        self.assertEqual(plant.state,checkpoint['frame']);self.assertEqual(plant.native.restored,{'native_tick':0})
        checkpoint['frame']['entities']['bone']['centroid_m'][0]=300
        self.assertEqual(plant.state['entities']['bone']['centroid_m'][0],1.)

    def test_embodied_frame_and_snapshot_remain_independent(self):
        body=runtime_fixture.Tests().body();result=body.step({});expected=deepcopy(body.frame)
        self.assertEqual(result,expected)
        result['mechanics']['entities']['injected']={'x':1};result['physiology']['values']['oxygen_saturation']=0
        self.assertEqual(body.frame,expected)
        public=body.snapshot();public['physiology']['values']['oxygen_saturation']=.1
        self.assertEqual(body.frame,expected)

    def test_respiration_geometry_and_skin_do_not_alias_owner_or_input(self):
        payload=fixture();owner=EmbodiedRespiration(payload,3000);mechanical=entities();original=deepcopy(mechanical)
        result=owner.geometry(3300,mechanical);second=owner.geometry(3300,mechanical)
        self.assertEqual(result,second);self.assertEqual(mechanical,original)
        result['entities']['skin']['centroid_m'][0]=100
        result['skin_field']['entity_transforms']['skin']['centroid_m'][0]=200
        result['skin_field']['bounds_m']['min'][0]=-100
        self.assertEqual(owner.geometry(3300,mechanical),second);self.assertEqual(mechanical,original)
        self.assertNotEqual(result['entities']['skin']['centroid_m'],result['skin_field']['entity_transforms']['skin']['centroid_m'])

if __name__=='__main__':unittest.main()
