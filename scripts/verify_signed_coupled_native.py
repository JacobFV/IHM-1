"""Bounded whole-engine signed-command smoke; no mechanical equivalence claim."""
from pathlib import Path
import argparse,hashlib,json,tempfile,time,resource
from ihm.native.session import SessionConfig,Meal
from ihm.native.coupled_session import SignedCoupledNativeSession

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',action='store_true');a=p.parse_args()
    if not a.run:raise SystemExit('Explicit --run and coordinated heavy slot required')
    root=Path(__file__).resolve().parents[1];out=Path(tempfile.mkdtemp(prefix='signed-coupled-',dir=root/'data/derived/audits'));start=time.monotonic()
    source=json.loads((root/'data/derived/systemic/exertion_v3/exercise/native/manifest.json').read_text());config=SessionConfig(state_path=source['configuration']['state_path'],engine_variant='whole_body_integrity_signed_muscle_v2',horizon_s=.14)
    ref=hashlib.sha256(b'explicit boundary protocol fixture; no mechanical reference calibration').hexdigest();s=None
    try:
        s=SignedCoupledNativeSession(config,out/'disabled');disabled=s.step(.02);s.close();s=None
        s=SignedCoupledNativeSession(config,out/'signed');zero=s.signed_step(ref,0,0,0)
        residual={k:(disabled['values'][k],v) for k,v in zero['values'].items() if not k.startswith('coupling.') and disabled['values'][k]!=v}
        assert not residual,dict(list(residual.items())[:5])
        frames=[zero]
        for m,h,w in ((.01,.01,0),(-.01,-.01,0),(.01,.02,-.01),(0,0,0)):
            frame=s.signed_step(ref,m,h,w);values=frame['values'];frames.append(frame)
            assert values['coupling.muscle_heat_count']==values['coupling.muscle_tissue_count']==1
            assert abs(values['coupling.effective_heat_w']-values['coupling.native_heat_w']-h)<1e-10
            assert values['coupling.muscle_unmet_kcal']<=1e-12
        ack=s.meal(Meal(name='protocol_water',water_ml=10))
        assert ack['pending_meal'] is True
        frame=s.signed_step(ref,0,0,0);frames.append(frame)
        assert frame['pending_meal'] is False and abs(frame['elapsed_s']-.12)<1e-10
        s.close();s=None
        record={'passed':True,'zero_nonboundary_ports_exact':True,'signed_steps':len(frames),'elapsed_s':frame['elapsed_s'],
            'wall_s':time.monotonic()-start,'parent_peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'scope':'Whole-engine signed boundaries, matched zero, return, eccentric heat and native meal consumption; not a mechanical resting reference or long-run validation','output':str(out.relative_to(root))}
        (out/'verification.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
    except BaseException as error:
        (out/'failure.json').write_text(json.dumps({'error':str(error),'wall_s':time.monotonic()-start},indent=2)+'\n');raise
    finally:
        if s:s.close(graceful=False)
if __name__=='__main__':main()
