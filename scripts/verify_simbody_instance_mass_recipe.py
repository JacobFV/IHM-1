"""Source-only recipe/immutability acceptance; never compiles native code."""
from pathlib import Path
import hashlib,json,shutil,tempfile,unittest
from patch_simbody_instance_mass import ROOT,RECIPE,patch
from build_simbody_instance_mass_variant import validate
class RecipeTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name);self.identity=json.loads(RECIPE.read_text())['source_sha256']
        self.raw=ROOT/'data/raw/mechanics/simbody/Simbody/src'
        for name in self.identity:shutil.copy2(self.raw/name,self.base/name)
    def tearDown(self):self.temp.cleanup()
    def test_exact_patch_changes_each_required_source(self):
        self.assertEqual(set(patch(self.base)),set(self.identity))
        for name,digest in self.identity.items():
            self.assertNotEqual(hashlib.sha256((self.base/name).read_bytes()).hexdigest(),digest)
            self.assertEqual(hashlib.sha256((self.raw/name).read_bytes()).hexdigest(),digest)
        self.assertIn('return iv.bodyMassProperties[getMyMobilizedBodyIndex()]',(self.base/'MobilizedBodyImpl.h').read_text())
        self.assertIn('getMk_G(pc).getMass() *',(self.base/'RigidBodyNode.cpp').read_text())
        self.assertNotIn('eps/getMass()',(self.base/'RigidBodyNode_LoneParticle.cpp').read_text())
    def test_unpinned_source_rejected_before_patch(self):
        p=self.base/'RigidBodyNode.cpp';p.write_text(p.read_text()+'\n')
        with self.assertRaises(ValueError):patch(self.base)
        self.assertEqual(hashlib.sha256((self.base/'MobilizedBodyImpl.h').read_bytes()).hexdigest(),self.identity['MobilizedBodyImpl.h'])
    def test_frozen_overlay_and_donor_guard(self):
        (self.base/'src').mkdir();(self.base/'objects').mkdir();(self.base/'src'/'x').write_text('source');(self.base/'objects'/'x.o').write_text('object')
        digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
        m={'inputs':{},'source_overlay_sha256':{'x':digest(self.base/'src'/'x')},'donor_objects':{'x.o':digest(self.base/'objects'/'x.o')},'built_objects':{}}
        validate(self.base,m);(self.base/'objects'/'x.o').write_text('mutated')
        with self.assertRaises(ValueError):validate(self.base,m)
if __name__=='__main__':unittest.main()
