"""Pressure command capability receipts reflect the actual owner; no native launch."""
import unittest
from ihm.assembly.embodied import EmbodiedRuntime
from verify_embodied_runtime import Plant,Neural,Native,Load,Exchange

class Tests(unittest.TestCase):
    def body(self,native):return EmbodiedRuntime(Plant(),Neural(),native,Exchange(),Load())

    def test_unknown_owner_never_claims_pressure_support(self):
        body=self.body(Native())
        self.assertEqual(body.snapshot()['input_capabilities']['skin_pressure'],{'whole_skin':False,'regional_ids':[]})
        self.assertEqual(body.step({})['input_capabilities'],body.snapshot()['input_capabilities'])

    def test_whole_and_regional_owners_use_actual_methods_and_ids(self):
        class Whole(Native):
            def skin_compression(self,value):pass
        class Regional(Whole):
            supports_whole_skin_compression=False
            regions=('configured_left','residual')
            def regional_skin_pressure(self,region,value):pass
        for native,expected in ((Whole(),{'whole_skin':True,'regional_ids':[]}),
                                (Regional(),{'whole_skin':False,'regional_ids':['configured_left','residual']})):
            with self.subTest(owner=type(native).__name__):
                body=self.body(native)
                self.assertEqual(body.snapshot()['input_capabilities']['skin_pressure'],expected)
                self.assertEqual(body.step({})['input_capabilities']['skin_pressure'],expected)
                self.assertEqual(body.snapshot()['input_capabilities']['skin_pressure'],expected)

    def test_malformed_owner_metadata_fails_closed(self):
        native=Native();native.regions=('bad','bad');native.regional_skin_pressure=lambda *a:None
        native.skin_compression=lambda *a:None;native.supports_whole_skin_compression=False
        self.assertEqual(self.body(native).snapshot()['input_capabilities']['skin_pressure'],{'whole_skin':False,'regional_ids':[]})

if __name__=='__main__':unittest.main()
