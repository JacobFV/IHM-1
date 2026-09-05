"""Canonical API job version pinning and optional fresh native execution."""
from pathlib import Path
import sys,json,tempfile,threading,unittest,argparse
from unittest.mock import patch
from urllib.request import Request,urlopen
from urllib.error import HTTPError
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.app import Jobs,create_server

class CanonicalJobTests(unittest.TestCase):
    def test_missing_body_returns_structured_error(self):
        with tempfile.TemporaryDirectory() as d:
            server=create_server(d,port=0);thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
            try:
                request=Request(f'http://127.0.0.1:{server.server_port}/api/body/scenarios',data=b'{"seconds": 1}',headers={'Content-Type':'application/json'})
                with self.assertRaises(HTTPError) as caught:urlopen(request)
                self.assertEqual(caught.exception.code,503)
                self.assertIn('artifact',json.load(caught.exception)['error'])
            finally:server.shutdown();server.server_close();thread.join()

    def test_queue_pins_sources_before_native_execution(self):
        with tempfile.TemporaryDirectory() as d:
            jobs=Jobs(ROOT);jobs.directory=Path(d);jobs.runs=[]
            with patch.object(jobs.pool,'submit') as enqueue:
                job=jobs.submit({'seconds':1,'engine_variant':'saturation_bounds_heatflux'},canonical=True)
                self.assertEqual(job['config']['patient'],'IHMGenericMale')
                self.assertIn('anatomy',job['canonical_sources'])
                fn,queued,config=enqueue.call_args.args
                queued['canonical_sources']={}
                with patch('ihm.native.run_native') as native:
                    fn(queued,config);native.assert_not_called()
                self.assertEqual(queued['status'],'failed');self.assertIn('changed',queued['error'])
            jobs.pool.shutdown()

    def test_snapshot_hashes_follow_consumed_bytes(self):
        from ihm.assembly.body import read_native,digest
        p=ROOT/'data/derived/canonical/native_baseline_v1'
        if not p.exists():self.skipTest('No native execution')
        native=read_native(p)
        self.assertEqual(native['input_hashes'],{name:digest(p/name) for name in ['native_multisystem.csv','body_compartments.csv','summary.json']})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--native',action='store_true');args=parser.parse_args()
    result=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(CanonicalJobTests))
    if not result.wasSuccessful():raise SystemExit(1)
    if args.native:
        import time
        jobs=Jobs(ROOT)
        job=jobs.submit({'seconds':2,'engine_variant':'saturation_bounds_heatflux'},canonical=True)
        print('Fresh canonical native job:',job['id'],flush=True)
        deadline=time.monotonic()+240
        while time.monotonic()<deadline:
            row=next(j for j in jobs.list()['runs'] if j['id']==job['id'])
            if row['status'] in ('failed','completed'):break
            time.sleep(.5)
        assert row['status']=='completed',row
        trajectory=json.loads((jobs.directory/job['id']/'body-trajectory.json').read_text())
        assert trajectory['clock']['end_s']==2 and trajectory['clock']['native_samples']==100
        from ihm.assembly.body import CanonicalBody
        result=CanonicalBody.from_workspace(ROOT).simulate(jobs.directory/job['id']/'output',ROOT/'artifacts/body-30hz.json',output_hz=30)
        assert len(result['frames'])==61,len(result['frames'])
        assert result['frames'][0]['time_s']==.02
        print('Fresh native execution, pinned canonical job, physiology/compartment receipt and anchored 30 Hz materialization PASS',flush=True)
        jobs.pool.shutdown()
