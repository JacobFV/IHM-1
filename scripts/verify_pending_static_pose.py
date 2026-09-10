"""One gated pending static-pose replay, with exact rejection receipt."""
from pathlib import Path
import argparse,hashlib,json,signal,sys,tempfile,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG

def run(path):
    from ihm.native.mechanical_stream import NativeMechanicalStream
    pending=json.loads(path.read_text());output=Path(tempfile.mkdtemp(prefix='pending-static-pose-',dir=ROOT/'data/derived'))
    stream=None;started=time.monotonic();report=dict(physical_time_advanced_s=0,pending_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),accepted_equilibrium=False)
    def timeout(signum,frame):
        if stream is not None and stream.process.poll() is None:stream.process.kill()
        raise TimeoutError('Pending static replay exceeded15s')
    old=signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,15)
    try:
        stream=NativeMechanicalStream(ROOT,output/'native',environment='supine',target_mass_kg=MECHANICAL_TARGET_MASS_KG,
            augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json',
            surface_contact_manifest='data/derived/supine-surface-contact-exmzq9pq/manifest.json',bed_material='MM')
        before=stream.snapshot();q=pending['coordinates'];command=['evaluate_static_pose',str(len(q))]
        for name,value in q.items():command.extend([name,str(value)])
        (output/'request.json').write_text(json.dumps(pending,indent=2)+'\n')
        try:
            observed=stream._request(' '.join(command));(output/'observed.json').write_text(json.dumps(observed,indent=2)+'\n')
            report['classification']='native_static_response'
        except ValueError as error:
            report.update(classification='native_rejection',native_failure_text=str(error))
            if isinstance(error,json.JSONDecodeError):report.update(classification='malformed_protocol',response_payload=error.doc[:65536])
        live=stream._request('observe');report['continuing_state_unchanged']=all(before[k]==live[k] for k in before if k!='kind')
        report['protocol_passed']=report['classification']!='malformed_protocol' and report['continuing_state_unchanged']
    except Exception as error:report.update(classification='fixture_failure',failure=f'{type(error).__name__}: {error}')
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if stream is not None:stream.close()
        report['wall_s']=time.monotonic()-started;(output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({**report,'output':str(output)},indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-native',action='store_true');parser.add_argument('pending',type=Path);args=parser.parse_args()
    if not args.run_native:raise SystemExit('Explicit coordinated --run-native required')
    run(args.pending.resolve())
