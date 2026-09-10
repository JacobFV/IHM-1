"""Bounded current-pose Jacobian and assembly-history numerical compatibility."""
from pathlib import Path
import argparse,hashlib,json,signal,sys,tempfile,time,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG
from frozen_static_stream import load,validate

ASSEMBLY_ACCURACY=1e-10


def compare_group(observed,reference,jacobian,rows,unit):
    # Each independently assembled coordinate may differ by accuracy in its
    # native coordinate unit: the pair envelope is twice that accuracy.
    rows=np.asarray(rows,bool);difference=np.asarray(observed)-np.asarray(reference)
    envelope=2*ASSEMBLY_ACCURACY*np.max(np.sum(np.abs(jacobian[rows,:]),axis=1))
    discrepancy=float(np.max(np.abs(difference[rows])))
    return dict(unit=unit,maximum_difference=discrepancy,jacobian_assembly_envelope=float(envelope),passed=bool(discrepancy<=envelope))


def run(manifest,candidate_path,history_path):
    expected=json.loads(candidate_path.read_text());history=json.loads(history_path.read_text())
    names=list(expected['coordinates']);requested=np.array([expected['coordinates'][n] for n in names])
    model=ET.parse(ROOT/'data/derived/constrained-supine-xojr73kr/native/inputs/subject_walk_scaled.osim').getroot()
    bounds={e.attrib['name']:tuple(map(float,e.findtext('range').split())) for e in model.iter('Coordinate') if e.find('range') is not None}
    archive,origin=validate(manifest)
    archived_cpp=(manifest.parent/'sources/native_mechanical_stream.cpp').read_text()
    assert 'model.set_assembly_accuracy(1e-10)' in archived_cpp
    output=Path(tempfile.mkdtemp(prefix='frozen-equivalence-investigation-',dir=ROOT/'data/derived'))
    started=time.monotonic();stream=None;count=0;records=[]
    report=dict(passed=False,accepted_equilibrium=False,physical_time_advanced_s=0,maximum_native_evaluations=40,maximum_wall_s=15,
        assembly_accuracy=ASSEMBLY_ACCURACY,comparison_basis='Grouped physical-unit infinity norms, bounded by 2*assembly_accuracy times measured current-pose Jacobian row absolute sums; no change to static or forward thresholds',
        limitation='Operational first-order assembly-error envelope, not a proof of bitwise equality or global nonlinear error bound')
    def write(name,value):(output/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    def deadline(signum,frame):
        if stream is not None and stream.process.poll() is None:stream.process.kill()
        raise TimeoutError('Frozen equivalence investigation exceeded15s')
    def evaluate(values,stage):
        nonlocal count
        if count>=40:raise RuntimeError('Frozen equivalence reached40native evaluations')
        count+=1;command=['evaluate_static_pose',str(len(names))]
        for name,value in zip(names,values):command.extend([name,str(float(value))])
        write('pending_request.json',dict(stage=stage,evaluation=count,coordinates=dict(zip(names,map(float,values)))))
        observed=stream._request(' '.join(command));records.append(dict(stage=stage,requested=values.tolist(),observed=observed))
        with (output/'observations.jsonl').open('a') as destination:destination.write(json.dumps(records[-1],allow_nan=False)+'\n')
        return observed
    old=signal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,15)
    try:
        stream=load(manifest)(ROOT,output/'native',environment='supine',target_mass_kg=MECHANICAL_TARGET_MASS_KG,
            augmented_registration='data/derived/mechanics/whole_body_arm26_v2/registration.json',
            surface_contact_manifest='data/derived/supine-surface-contact-exmzq9pq/manifest.json',bed_material='MM')
        live_before=stream.snapshot();base=evaluate(requested,'fresh requested pose');samples=[]
        for i,name in enumerate(names):
            lo,hi=bounds[name];step=1e-5 if requested[i]+1e-5<=hi else -1e-5
            shifted=requested.copy();shifted[i]+=step;samples.append(evaluate(shifted,'current Jacobian '+name))
        comparisons={'cached_multi_evaluation_pose':expected['native'], 'after_jacobian_history':evaluate(requested,'repeat after Jacobian')}
        evaluate(np.array([history['coordinates'][n] for n in names]),'different supported history')
        comparisons['after_different_supported_history']=evaluate(requested,'repeat after different supported history')
        backward={}
        for name in ('pelvis_tx','pelvis_tilt','pelvis_rotation'):
            shifted=requested.copy();shifted[names.index(name)]-=1e-5;backward[name]=evaluate(shifted,'backward root probe '+name)
        comparisons['after_backward_history']=evaluate(requested,'repeat after backward probes')
        actual_q=lambda native:np.array([native['coordinates'][n]['value'] for n in names])
        D=np.column_stack([actual_q(sample)-actual_q(base) for sample in samples]);report['actual_coordinate_difference_condition']=float(np.linalg.cond(D))
        if not np.isfinite(np.linalg.cond(D)) or np.linalg.cond(D)>10:raise ValueError('Current-pose perturbation matrix is ill-conditioned')
        rotational=np.asarray(base['mobility_rotational'],bool);independent_rotational=np.array([base['coordinates'][n]['rotational'] for n in names],bool)
        fields=('udot','constrained_zero_acceleration_residual_mobility_force','mass_times_udot','contact_force_n')
        jacobians={key:np.column_stack([np.asarray(sample[key])-np.asarray(base[key]) for sample in samples])@np.linalg.inv(D) for key in fields}
        report['comparisons']={}
        for label,observed in comparisons.items():
            dq=actual_q(observed)-actual_q(base);entry={'coordinates':{}}
            for mask,unit in [(independent_rotational,'rad'),(~independent_rotational,'m')]:
                discrepancy=float(np.max(np.abs(dq[mask])));entry['coordinates'][unit]=dict(maximum_difference=discrepancy,pair_assembly_envelope=2*ASSEMBLY_ACCURACY,passed=discrepancy<=2*ASSEMBLY_ACCURACY)
            for key in fields:
                if key=='contact_force_n':groups=[(np.ones(3,bool),'N')]
                elif key=='udot':groups=[(rotational,'rad/s^2'),(~rotational,'m/s^2')]
                else:groups=[(rotational,'N*m'),(~rotational,'N')]
                entry[key]=[compare_group(observed[key],base[key],jacobians[key],mask,unit) for mask,unit in groups]
                entry[key+'_linear_prediction_error_max']=float(np.max(np.abs(np.asarray(observed[key])-np.asarray(base[key])-jacobians[key]@dq)))
            report['comparisons'][label]=entry
        report['passed']=all(all(x['passed'] for x in entry['coordinates'].values()) and all(x['passed'] for key in fields for x in entry[key]) for entry in report['comparisons'].values())
        live=stream._request('observe');report['continuing_state_unchanged']=all(live_before[k]==live[k] for k in live_before if k!='kind')
        report['passed']=report['passed'] and report['continuing_state_unchanged']
        write('current_jacobians.json',{key:value.tolist() for key,value in jacobians.items()})
    except Exception as error:report['error']=f'{type(error).__name__}: {error}'
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,old)
        if stream is not None:stream.close()
        validate(manifest);report.update(archive_intact=True,evaluations=count,wall_s=time.monotonic()-started,
            source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (manifest,candidate_path,history_path,Path(__file__).resolve())})
        write('report.json',report);print(json.dumps({key:value for key,value in report.items() if key!='comparisons'}|{'output':str(output)},indent=2))
    return report['passed']


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-native',action='store_true');parser.add_argument('--manifest',type=Path,required=True);parser.add_argument('--candidate',type=Path,required=True);parser.add_argument('--history',type=Path,required=True);args=parser.parse_args()
    if not args.run_native:raise SystemExit('Explicit bounded native grant required')
    raise SystemExit(0 if run(args.manifest.resolve(),args.candidate.resolve(),args.history.resolve()) else 1)
