"""Kernel path ablation and learning regression, no native installation needed."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import torch
from ihm.native.motor_learning import KernelAnklePolicy
class Tests(unittest.TestCase):
    def test_sever(self):
        torch.manual_seed(5)
        p=KernelAnklePolicy(torch.randn(8,16));x=torch.tensor([[.2,.1],[-.2,-.1],[.1,-1.]])
        self.assertTrue(torch.equal(p(x,sever=True),torch.zeros(3,2)))
        self.assertGreater(torch.count_nonzero(p(x)),0)
        self.assertEqual(list(dict(p.named_parameters())),['embed'])
    def test_learning(self):
        torch.manual_seed(17);torch.set_num_threads(1)
        p=KernelAnklePolicy(torch.randn(8,16)*.02);x=torch.rand(128,2)*.4-.2
        d=2*x[:,0]-.1*x[:,1];y=torch.stack((d.clamp(0,1),(-d).clamp(0,1)),-1)
        before=p.embed.detach().clone();initial=float((p(x)-y).square().mean().detach())
        opt=torch.optim.Adam(p.parameters(),lr=.003)
        for _ in range(100):
            loss=(p(x)-y).square().mean();opt.zero_grad();loss.backward();opt.step()
        self.assertLess(float((p(x)-y).square().mean().detach()),initial*.1)
        self.assertFalse(torch.equal(before,p.embed))
    def test_primitive_blocks_and_checkpoint(self):
        from ihm.native.motor_learning_controller import IBMAnklePrimitiveController
        from types import SimpleNamespace
        class Inner:
            muscles=('tibant_r','soleus_r');sever=False
            def __init__(self):self.time_s=0.
            def checkpoint(self):return {'time_s':self.time_s}
            def restore(self,c):self.time_s=c['time_s']
            def step(self,dt,observation,**kwargs):
                self.time_s+=dt
                return {'brain':{'oxygen_perfusion_availability':1.},'cortical_commands':{},'arc_max':{'stretch':.2}}
        p=IBMAnklePrimitiveController();p.inner=Inner();p.policy=KernelAnklePolicy(torch.randn(8,16));p.target_rad=.12
        p.artifact_sha256='test';p.model_sha256='test-model';p.controller_metadata={};p.excitations={n:0. for n in p.inner.muscles}
        # Execute the real projection method so native-to-runtime key drift is tested.
        from ihm.assembly.articulated import ArticulatedBodyPlant
        import numpy as np
        native={key:0. for key in ('mass_kg','time_s','signed_active_fiber_work_j','muscle_heat_energy_j','signed_active_fiber_power_w','positive_active_fiber_work_j','contact_force_n','momentum_balance_residual_n','constraint_position_error','constraint_velocity_error','external_work_j','external_power_w')}
        native.update(mass_transfer={},muscles=dict(p.excitations),foot_contact_force_n={},metabolic_reference={},
            coordinates={'ankle_angle_r':{'value':0.,'speed':0.}},bodies={},environment='free',contacts=[],gravity_m_s2=[0.,-9.81,0.])
        plant=ArticulatedBodyPlant.__new__(ArticulatedBodyPlant)
        plant.registration=SimpleNamespace(project=lambda n:{},contact_wrenches=lambda n:[],cutaneous_contacts=lambda n:None,basis=np.eye(3));plant.garments=None
        observation=plant._project(native,0.)
        self.assertIn('joints',observation);self.assertNotIn('coordinates',observation)
        saved=p.checkpoint()
        result=p.step(.01,observation,sensory_blocks=['tibant_r'])
        self.assertEqual(result['motor_excitations'],{'tibant_r':0.,'soleus_r':0.})
        self.assertEqual(result['arc_max'],{'stretch':0.})
        self.assertEqual(result['inactive_cord_arc_max'],{'stretch':.2})
        p.restore(saved);self.assertEqual(p.time_s,0.)
        with self.assertRaises(ValueError):p.step(.01,dict(observation,coordinates={'ankle_angle_r':{'value':.3,'speed':0.}}))
        self.assertEqual(p.time_s,0.)
        result=p.step(.01,observation,motor_blocks=['tibant_r','soleus_r'])
        self.assertEqual(result['motor_excitations'],{'tibant_r':0.,'soleus_r':0.})
    def test_retains_loaded_artifact_and_sources(self):
        from ihm.native.motor_learning_controller import IBMAnklePrimitiveController
        from types import SimpleNamespace
        from unittest.mock import patch
        import tempfile,hashlib,json
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);artifact=folder/'model.pt'
            torch.save({'schema':'ihm.ibm-native-ankle-kernel.v1','embed':torch.randn(8,16),'provenance':{}},artifact)
            original=artifact.read_bytes()
            inner=SimpleNamespace(identity={},controller_metadata={},excitations={'tibant_r':0.,'soleus_r':0.},retain_sources=lambda output:None)
            with patch('ihm.assembly.ibm_controller.IBMImplicitController.from_root',return_value=inner):
                p=IBMAnklePrimitiveController.from_root(folder,muscle_catalog=[],artifact_path=artifact)
            artifact.write_bytes(b'changed after load')
            retained=folder/'retained';retained.mkdir();p.retain_sources(retained)
            saved=retained/'ankle-motor-primitive'
            self.assertEqual((saved/'motor_kernel.pt').read_bytes(),original)
            self.assertEqual(hashlib.sha256(original).hexdigest(),p.identity['artifact_sha256'])
            self.assertEqual(hashlib.sha256((saved/'motor_learning.py').read_bytes()).hexdigest(),p.identity['policy_sha256'])
            self.assertEqual(hashlib.sha256((saved/'motor_learning_controller.py').read_bytes()).hexdigest(),p.identity['adapter_sha256'])
if __name__=='__main__':unittest.main()
