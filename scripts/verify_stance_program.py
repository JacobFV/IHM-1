"""Compare reference scheduler against immutable accepted evaluator expressions."""
import ast,json,sys,unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.stance_program import StanceReferenceProgram
ROOT=Path(__file__).resolve().parents[1]
class Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.program=StanceReferenceProgram.from_bundle(ROOT)
        bundle=ROOT/'data/models/engineering_weight_transfer_v1'
        with np.load(bundle/'provenance/baseline_policy.npz') as d:cls.base={k:d[k].copy() for k in d.files}
        with np.load(bundle/'target.npz') as d:cls.target={k:d[k].copy() for k in d.files}
        with np.load(bundle/'target_policy.npz') as d:cls.goal={k:d[k].copy() for k in d.files}
        tree=ast.parse((bundle/'provenance/evaluator.py').read_text())
        run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run')
        loop=next(n for n in ast.walk(run) if isinstance(n,ast.While))
        statements=[]
        for statement in loop.body:
            if isinstance(statement,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='state' for t in statement.targets):break
            statements.append(statement)
        cls.accepted_code=compile(ast.fix_missing_locations(ast.Module(body=statements,type_ignores=[])),'immutable_accepted_evaluator','exec')
    def test_matches_accepted_commands_across_phases(self):
        base,target,goal=self.base,self.target,self.goal
        fixed=base['x0']+np.linspace(-.00001,.00001,len(base['x0']))
        policy=SimpleNamespace(x0=base['x0'],u0=base['u0'],K=base['K'],state_names=base['state_names'].tolist(),muscle_names=base['muscle_names'].tolist(),state_vector=lambda state:fixed)
        for t in (0.,.99,1.,1.01,2.,3.5,5.99,6.,10.,1e6):
            namespace={'np':np,'stream':SimpleNamespace(state={'time_s':t}),'transition_s':5.,'fraction':1.,'policy':policy,
                'delta':target['x_target']-base['x0'],'du':target['u_target']-base['u0'],'goal_gain':goal['K'],
                'values':{path.rsplit('/',1)[0]:i for i,path in enumerate(policy.state_names) if path.endswith('/value')}}
            exec(self.accepted_code,namespace)
            observed=self.program.reference_at(t)
            np.testing.assert_array_equal(observed['state_reference'],namespace['reference'])
            np.testing.assert_array_equal(observed['gain'],namespace['gain'])
            actual=observed['equilibrium_excitations']-observed['gain']@(fixed-observed['state_reference'])
            np.testing.assert_array_equal(actual,namespace['raw'])
            self.assertTrue(np.isfinite(actual).all());json.dumps(observed['telemetry'],allow_nan=False)
    def test_time_checkpoint_and_endpoint(self):
        p=self.program;p.restore({'schema':'ihm.stance-reference-program-state.v1','model_sha256':p.model_sha256,'time_s':0.})
        saved=p.checkpoint();first=p.advance(.01);self.assertEqual(first['telemetry']['elapsed_s'],0.);self.assertEqual(p.time_s,.01)
        p.restore(saved);self.assertEqual(p.time_s,0.)
        for dt in (.02,float('nan'),True):
            with self.assertRaises(ValueError):p.advance(dt)
        with self.assertRaises(ValueError):p.restore(dict(saved,model_sha256='wrong'))
        endpoint=p.reference_at(6.);later=p.reference_at(1000.)
        np.testing.assert_array_equal(endpoint['state_reference'],later['state_reference'])
        self.assertEqual(later['telemetry']['blend_rate_per_s'],0.);self.assertEqual(later['telemetry']['phase'],'endpoint_hold')
        np.testing.assert_allclose(p.reference_at(0.)['state_reference'],self.base['x0'],rtol=0,atol=0)
        self.assertGreater(np.max(np.abs(p.reference_at(0.)['state_reference']-self.target['x0'])),.1)
if __name__=='__main__':unittest.main()
