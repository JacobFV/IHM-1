"""Source-only frozen adapter identity and tamper-rejection fixtures."""
from pathlib import Path
import argparse,copy,json,signal,sys,tempfile,time
import numpy as np
from frozen_static_stream import validate,load


def main(path):
    path=Path(path);record,origin=validate(path);klass=load(path)
    assert callable(klass) and len(record['original_build_validation_map'])==len(origin['files'])
    temporary=path.parent/'tampered-fixture-manifest.json'
    try:
        for key in ['adapter_original.py','adapter_frozen.py',next(k for k in record['files'] if k.startswith('libraries/'))]:
            changed=copy.deepcopy(record);changed['files'][key]='0'*64;temporary.write_text(json.dumps(changed))
            try:validate(temporary)
            except ValueError:pass
            else:raise AssertionError('Tampered frozen record accepted: '+key)
        changed=copy.deepcopy(record);changed['original_build_validation_map'].pop(next(iter(changed['original_build_validation_map'])))
        temporary.write_text(json.dumps(changed))
        try:validate(temporary)
        except ValueError:pass
        else:raise AssertionError('Missing compiled dependency accepted')
    finally:temporary.unlink(missing_ok=True)
    print(json.dumps({'passed':True,'native_run':False,'compiled_dependencies_verified':len(origin['files']),'copied_files':len(record['files']),'tamper_checks':4}))

def native_check(path,candidate_path):
    root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
    expected=json.loads(candidate_path.read_text());klass=load(path)
    output=Path(tempfile.mkdtemp(prefix='frozen-static-acceptance-',dir=root/'data/derived'))
    started=time.monotonic();stream=None;report=dict(passed=False,physical_time_advanced_s=0,accepted_equilibrium=False)
    def deadline(signum,frame):
        if stream is not None and stream.process.poll() is None:stream.process.kill()
        raise TimeoutError('Frozen static attestation acceptance exceeded15s')
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,15)
    try:
        stream=klass(root,output/'native',environment='supine',target_mass_kg=77.6122029,
            augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json',
            surface_contact_manifest='data/derived/supine-surface-contact-exmzq9pq/manifest.json',bed_material='MM')
        before=stream.snapshot();q=expected['coordinates'];command=['evaluate_static_pose',str(len(q))]
        for name,value in q.items():command.extend([name,str(value)])
        observed=stream._request(' '.join(command));(output/'observed.json').write_text(json.dumps(observed,indent=2)+'\n')
        report['maximum_numeric_differences']={key:float(np.max(np.abs(np.asarray(observed[key])-np.asarray(expected['native'][key]))))
            for key in ('q','u','udot','mass_times_udot','constrained_zero_acceleration_residual_mobility_force','contact_force_n')}
        assert all(value<=1e-9 for value in report['maximum_numeric_differences'].values()),report
        live=stream._request('observe');assert all(before[k]==live[k] for k in before if k!='kind')
        report.update(passed=True,continuing_state_unchanged=True)
    except Exception as error:report['error']=f'{type(error).__name__}: {error}'
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if stream is not None:stream.close()
        validate(path);report['archive_intact_after_run']=True;report['wall_s']=time.monotonic()-started
        (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({**report,'output':str(output)},indent=2))
    return report['passed']


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('manifest',type=Path);parser.add_argument('--run-native',action='store_true');parser.add_argument('--candidate',type=Path);args=parser.parse_args()
    if args.run_native:
        if args.candidate is None:raise SystemExit('Explicit candidate required')
        raise SystemExit(0 if native_check(args.manifest.resolve(),args.candidate.resolve()) else 1)
    main(args.manifest)
