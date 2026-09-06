"""Mass-metric fixtures and optionally one bounded native identity observation."""
from pathlib import Path
import argparse,json,signal,sys,tempfile,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'scripts'))
from solve_constrained_supine_pose import dynamic_metric


def fixtures():
    matrix=np.array([[2.,.4],[.4,.3]]);acceleration=np.array([3.,-4.]);force=matrix@acceleration
    native=dict(constrained_zero_acceleration_residual_mobility_force=(-force).tolist(),udot=acceleration.tolist(),
         mass_times_udot=force.tolist(),mass_kg=2.,gravity_m_s2=[-10.,0,0])
    assert np.isclose(dynamic_metric(native),float(force@np.linalg.solve(matrix,force))/200)
    altered={**native,'mass_times_udot':(force+1).tolist()}
    try:dynamic_metric(altered)
    except ValueError:pass
    else:raise AssertionError('Broken r=-M udot identity accepted')
    negative={**native,'mass_times_udot':(-acceleration).tolist(),'constrained_zero_acceleration_residual_mobility_force':acceleration.tolist()}
    try:dynamic_metric(negative)
    except ValueError:pass
    else:raise AssertionError('Negative mass metric clamped/accepted')
    return dict(passed=True,native_run=False,checks=['primal-dual SPD identity','constraint/sign identity rejection','negative cost rejection without floor'])


def native_fixture():
    from ihm.native.mechanical_stream import NativeMechanicalStream
    parent=Path(tempfile.mkdtemp(prefix='native-static-metric-',dir=ROOT/'data/derived'));stream=None;start=time.monotonic()
    report=dict(passed=False,physical_time_advanced_s=0)
    def deadline(signum,stack):
        if stream is not None and stream.process.poll() is None:stream.process.kill()
        raise TimeoutError('Native metric identity fixture wall cap')
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,15)
    try:
        stream=NativeMechanicalStream(ROOT,parent/'native',environment='supine',target_mass_kg=77.6122029,
            augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json',
            surface_contact_manifest='data/derived/supine-surface-contact-exmzq9pq/manifest.json',bed_material='MM')
        before=stream.snapshot();seed=json.loads((ROOT/'data/derived/supported-rigid-seed-ljyewh9y/report.json').read_text())['best']['seed_coordinates']
        command=['evaluate_static_pose',str(len(seed))]
        for name,value in seed.items():command.extend([name,str(value)])
        observed=stream._request(' '.join(command));(parent/'observed.json').write_text(json.dumps(observed,indent=2)+'\n')
        report['dynamic_metric']=dynamic_metric(observed);report['identity_residual_norm']=observed['dynamic_residual_identity_norm']
        live=stream._request('observe');assert all(before[k]==live[k] for k in before if k!='kind')
        report.update(passed=True,continuing_state_unchanged=True)
    except Exception as error:report['error']=f'{type(error).__name__}: {error}'
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if stream is not None:stream.close()
        report['wall_s']=time.monotonic()-start;(parent/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({**report,'output':str(parent)},indent=2))
    return report['passed']


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-native',action='store_true');args=parser.parse_args()
    if args.run_native:raise SystemExit(0 if native_fixture() else 1)
    print(json.dumps(fixtures(),indent=2))
