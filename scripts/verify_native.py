#!/usr/bin/env python3
"""Contract checks plus optional real-engine experiments (--engine)."""
import json, tempfile, unittest, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native import NativeConfig, Intervention, load_trajectory, run_native, summarize
class Contracts(unittest.TestCase):
    def test_bounds(self):
        for kw in ({'seconds':0},{'seconds':float('nan')},{'seconds':True},{'sample_hz':13},{'patient':'../StandardMale'},{'seconds':0.03}):
            with self.assertRaises(ValueError): NativeConfig(**kw)
        for kw in ({'unknown':2},{'interventions':[{'time_s':0,'kind':'unknown','value':1}]},{'interventions':[{'time_s':61,'kind':'exercise','value':.2}]}):
            with self.assertRaises(ValueError): NativeConfig.from_dict(kw)
        with self.assertRaises(ValueError): Intervention(1,'exercise',.8)
    def test_engine_error_despite_finite_csv(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'native_multisystem.csv').write_text('#Time(s),X\n0.02,0\n')
            (p/'runner_stdout.log').write_text('ENGINE_VERSION=test\nFATAL Unknown Data Request : X\n')
            with self.assertRaisesRegex(ValueError,'Engine logged'):summarize(p)
    def test_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.csv';p.write_text('Time(s),X(mmHg)\n0.02,3\n0.04,4\n')
            t=load_trajectory(p);self.assertEqual(t['values']['X(mmHg)'],[3.,4.])
            p.write_text('Time(s),X\n0.02,nan\n')
            with self.assertRaises(ValueError):load_trajectory(p)
            p.write_text('Time(s),X\n0.04,1\n0.02,2\n')
            with self.assertRaises(ValueError):load_trajectory(p)
def verify_engine_clock():
    from ihm.native import BASE
    state=BASE/'data/derived/physiology/native_baseline_v2/states/native_stabilized.xml'
    with tempfile.TemporaryDirectory(prefix='native-clock-') as directory:
        for hz in (1,2,5,10,25,50):
            config=NativeConfig(seconds=2.04,state_path=str(state),sample_hz=hz,
                interventions=(Intervention(.4,'exercise',.1),Intervention(1.4,'exercise',0)))
            out=Path(directory)/str(hz);summary=run_native(config,out)
            assert summary['rows']==102//(50//hz)
            log=(out/'runner_stdout.log').read_text()
            assert 'FINAL_TIME_S=2.04' in log
            assert 'ACTION_TIME_S=0.4' in log and 'ACTION_TIME_S=1.4' in log
    print('Native integer scheduling and all sample cadences PASS')
if __name__=='__main__':
    if '--engine-clock' in sys.argv:
        sys.argv.remove('--engine-clock');verify_engine_clock()
    unittest.main()
