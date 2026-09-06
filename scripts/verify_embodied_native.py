"""Explicit opt-in tiny integrated native run; retains every source/command."""
from pathlib import Path
import argparse,json,tempfile,time,resource
from ihm.assembly.embodied import EmbodiedRuntime


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',action='store_true');args=p.parse_args()
    if not args.run:raise SystemExit('Use --run only in the coordinated single native-job window')
    root=Path(__file__).resolve().parents[1]
    parent=Path(tempfile.mkdtemp(prefix='embodied-native-',dir=root/'data/derived/audits'))
    start=time.monotonic();body=None
    try:
        body=EmbodiedRuntime.from_workspace(root,parent/'body',environment='supine')
        initial=body.snapshot();frames=[]
        for _ in range(3):frames.append(body.step({}))
        assert frames[-1]['time_s']==.06
        assert all(abs(f['physiology']['elapsed_s']-f['time_s'])<1e-10 and abs(f['mechanics']['time_s']-f['time_s'])<1e-10 for f in frames)
        assert len(frames[-1]['entities'])==len(initial['entities'])
        assert frames[-1]['physiology']['values']['lung_volume_ml']!=initial['physiology']['values']['lung_volume_ml']
        record={'passed':True,'directory':str(parent.relative_to(root)),'duration_s':.06,
            'entities':len(initial['entities']),'native_muscle_effectors':len(initial['mechanics']['muscles']),
            'wall_s':time.monotonic()-start,'parent_peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'scope':'Three shared-clock intervals; no equilibrium, gait, long-horizon or biological calibration acceptance'}
        (parent/'verification.json').write_text(json.dumps(record,indent=2)+'\n')
        print(json.dumps(record,indent=2))
    except BaseException as error:
        (parent/'failure.json').write_text(json.dumps({'error':str(error),'wall_s':time.monotonic()-start},indent=2)+'\n')
        raise
    finally:
        if body:body.close()

if __name__=='__main__':main()
