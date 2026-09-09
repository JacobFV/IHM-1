"""Execute source-attached TA mechanics → stretch feedback → motor force experiments."""
from copy import deepcopy
from dataclasses import asdict
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import urllib.request
import numpy as np
from ihm.assembly.mechanics import BodyMechanics
from ihm.assembly.reflexes import BodyReflex, ReflexParameters, GEYER_HERR_URL

MECHANICS_STATE=('x','v','omega','r','deformation','time','work','dissipation','affine_work',
    '_soft_energy','prescribed_work','_soft_previous','prescribed_reactions','orientation_reactions')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def subdomain(root):
    path=Path(root)/'data/derived/canonical/mechanics.json'
    source=json.loads(path.read_text())
    for source_path,digest in source['source_files'].items():
        if sha(Path(root)/source_path)!=digest:raise ValueError('Canonical mechanics source changed: '+source_path)
    peripheral=json.loads((Path(root)/'data/derived/canonical/peripheral.json').read_text())
    muscle=deepcopy(next(m for m in source['muscles'] if m.get('source_name')=='tibant_r'))
    binding=next(m for m in peripheral['muscle_bindings'] if m['muscle_id']==muscle['id'])
    nerve=binding['nerve_id']
    muscle['neural_conduction']={'delays_s':binding['delays_s']}
    endpoint=muscle['anchors'][-1]['entity_id']
    support=deepcopy(next(link for link in source['links'] if link['kind']=='inferred_skeletal_support' and endpoint in (link['a'],link['b'])))
    ids={a['entity_id'] for a in muscle['anchors']}|{support['a'],support['b'],muscle['canonical_entity_id']}
    entities=[deepcopy(e) for e in source['entities'] if e['id'] in ids]
    if any(e['constitutive']!='rigid' for e in entities if e['id']!=muscle['canonical_entity_id']):raise ValueError('Expected rigid bone attachments')
    payload=dict(model_id='canonical-right-TA-isolated-subdomain',entities=entities,links=[support],muscles=[muscle])
    fixed=sorted(ids-{endpoint})
    points=muscle['anchors'];direction=np.asarray(points[-1]['point_m'])-points[-2]['point_m'];direction/=np.linalg.norm(direction)
    return payload,muscle,nerve,fixed,endpoint,direction


def observe(mechanics,muscle):
    return dict(time_s=mechanics.time,entities={key:dict(translation_m=(mechanics.x[i]-mechanics.x0[i]).tolist(),
        rotation_matrix=mechanics.r[i].tolist(),deformation_gradient=mechanics.deformation[i].tolist()) for i,key in enumerate(mechanics.ids)},
        muscle_forces_n={muscle['id']:0.})


# Keep the 1 mm time-refinement bound fixed while route-derived neural timing
# changes. Delay causality is checked against the actual Ia/alpha loop below.
def run_trial(payload,muscle,nerve,fixed,endpoint,direction,*,dt=.00005,duration=.4,mode='intact',resume_test=False):
    params=ReflexParameters.from_binding(muscle['neural_conduction'],extra_delay_s=.04 if mode=='delayed' else 0.)
    mechanics=BodyMechanics.from_dict(deepcopy(payload))
    controller=BodyReflex(muscle,{e['id']:e['centroid_m'] for e in payload['entities']},params,nerve)
    mechanics.max_substep=dt
    state=observe(mechanics,muscle);activation=0.;frames=[];external_work=0.;active_work=0.
    checkpoint=None;continued=[];replayed=[]
    n=round(duration/dt);start=round(.1/dt);stop=round(.25/dt);checkpoint_tick=round(.2/dt)
    def one(i):
        nonlocal state,activation,external_work,active_work
        time=i*dt;force=20. if start<=i<stop and mode!='no_load' else 0.
        # Action-switch contrasts share prior dynamics; delayed mode changes a model parameter from t=0.
        perturbed=i>=start
        switches=dict(motor_block=perturbed and (mode=='nerve_block' or (mode=='nerve_release' and i<stop)),sensory_block=perturbed and mode=='sensory_block',
            controller_enabled=not(perturbed and mode=='removed'),descending_gain=.5 if perturbed and mode=='descending_half' else 1.)
        neural=controller.step(dt,state,**switches)
        applied=activation
        prior_work=mechanics.work;prior_prescribed_work=mechanics.prescribed_work
        state=mechanics.step(dt,{'activation':{muscle['id']:applied},'external_forces_n':{endpoint:(force*direction).tolist()},
            'prescribed_translations_m':{key:[0.,0.,0.] for key in fixed}})
        increment=dt*float(np.dot(force*direction,mechanics.v[mechanics.index[endpoint]]))
        external_work+=increment
        active_work+=mechanics.work-prior_work-increment-(mechanics.prescribed_work-prior_prescribed_work)
        activation=neural['activation']
        sensed=controller.observe(state)
        return dict(time_s=(i+1)*dt,external_force_n=force,external_work_j=external_work,active_mechanical_work_j=active_work,activation=activation,applied_activation=applied,stimulation=neural['stimulation'],
            path_length_m=sensed['path_length_m'],normalized_fiber_length=sensed['normalized_fiber_length'],
            delayed_length=neural['delayed_normalized_length'],force_n=state['muscle_forces_n'][muscle['id']],
            endpoint_translation_m=state['entities'][endpoint]['translation_m'],
            **{key:state['audit'][key] for key in ('internal_force_residual_n','internal_torque_residual_nm','elastic_energy_j','kinetic_energy_j',
                'accumulated_active_external_work_j','accumulated_dissipation_j','energy_balance_residual_j')})
    for i in range(n):
        frame=one(i);frames.append(frame)
        if resume_test and i+1==checkpoint_tick:
            checkpoint=dict(controller=controller.checkpoint(),mechanics={k:deepcopy(getattr(mechanics,k)) for k in MECHANICS_STATE if hasattr(mechanics,k)},
                state=deepcopy(state),activation=activation,external_work=external_work,active_work=active_work)
        elif checkpoint is not None:continued.append(frame)
    if resume_test:
        controller.restore(checkpoint['controller'])
        for key,value in checkpoint['mechanics'].items():setattr(mechanics,key,deepcopy(value))
        state=deepcopy(checkpoint['state']);activation=checkpoint['activation'];external_work=checkpoint['external_work'];active_work=checkpoint['active_work']
        for i in range(checkpoint_tick,n):replayed.append(one(i))
    return dict(mode=mode,dt_s=dt,duration_s=duration,parameters=asdict(params),frames=frames,
        checkpoint_equal=continued==replayed if resume_test else None,
        final_controller_checkpoint=controller.checkpoint())


def run_suite(root,output_dir):
    root=Path(root).resolve();parent=Path(output_dir);parent.mkdir(parents=True,exist_ok=True)
    out=Path(tempfile.mkdtemp(prefix='run-',dir=parent))
    paths=['ihm/assembly/reflexes.py','ihm/assembly/mechanics.py','scripts/build_neuromechanical_experiments.py',
        'data/derived/canonical/mechanics.json','data/derived/canonical/peripheral.json']
    paper=root/'data/derived/audits/neuromechanical/sources/geyer_herr_2010.pdf'
    if not paper.exists():
        paper.parent.mkdir(parents=True,exist_ok=True)
        with urllib.request.urlopen(GEYER_HERR_URL,timeout=60) as response:paper.write_bytes(response.read())
    if sha(paper)!='71619c24b0c7f2f4a0c2eca2149e6b1a17542487a82a558f3bff140327e7d381':
        raise ValueError('Published controller source digest differs from inspected paper')
    paths.append(str(paper.relative_to(root)))
    receipts={p:sha(root/p) for p in paths}
    payload,muscle,nerve,fixed,endpoint,direction=subdomain(root)
    modes=('intact','no_load','nerve_block','nerve_release','sensory_block','removed','descending_half','delayed')
    trials={mode:run_trial(payload,muscle,nerve,fixed,endpoint,direction,mode=mode,resume_test=mode=='intact') for mode in modes}
    refined=run_trial(payload,muscle,nerve,fixed,endpoint,direction,dt=.000025)
    for name,trial in {**trials,'refined':refined}.items():
        (out/(name+'.json')).write_text(json.dumps(trial,indent=2,allow_nan=False)+'\n')
    def difference(a,b,key):
        return max(abs(u[key]-v[key]) for u,v in zip(trials[a]['frames'],trials[b]['frames']) if u['time_s']>.1)
    checks={mode+'_finite':all(np.isfinite([f['path_length_m'],f['force_n'],f['activation'],f['kinetic_energy_j']]).all() and f['path_length_m']>0 and 0<=f['activation']<=1 for f in trial['frames']) for mode,trial in trials.items()}
    checks['checkpoint_exact']=trials['intact']['checkpoint_equal']
    for mode in ('nerve_block','nerve_release','sensory_block','removed','descending_half','delayed'):
        checks[mode+'_changes_motion']=difference('intact',mode,'path_length_m')>1e-5
    for mode in ('nerve_block','nerve_release','sensory_block','removed','descending_half'):
        checks[mode+'_matched_before_intervention']=all(a==b for a,b in zip(trials['intact']['frames'],trials[mode]['frames']) if a['time_s']<=.1)
    changes=[a['time_s'] for a,b in zip(trials['intact']['frames'],trials['no_load']['frames']) if abs(a['activation']-b['activation'])>1e-10]
    activation_onset=changes[0] if changes else None
    checks['load_to_motor_delay_respected']=activation_onset is not None and activation_onset>=.1+trials['intact']['parameters']['loop_delay_s']-1e-10
    checks['external_load_changes_motor_activation']=difference('intact','no_load','activation')>1e-4
    checks['nerve_release_restores_motor']=any(f['activation']>.05 for f in trials['nerve_release']['frames'] if f['time_s']>.3)
    checks['external_load_changes_sensory']=difference('intact','no_load','normalized_fiber_length')>1e-4
    checks['internal_force_balance']=max(f['internal_force_residual_n'] for t in trials.values() for f in t['frames'])<1e-8
    checks['internal_torque_balance']=max(f['internal_torque_residual_nm'] for t in trials.values() for f in t['frames'])<1e-8
    refinement=max(abs(a['path_length_m']-b['path_length_m']) for a,b in zip(trials['intact']['frames'],refined['frames'][1::2]))
    checks['time_refinement_path_under_1mm']=refinement<.001
    maxdiff={mode:difference('intact',mode,'path_length_m') for mode in modes if mode!='intact'}
    report=dict(schema='ihm.neuromechanical.v1',output=str(out),checks=checks,passed=all(checks.values()),
        motor_response_onset_s=activation_onset,maximum_path_difference_m=maxdiff,time_refinement_maximum_path_error_m=refinement,
        muscle=muscle,nerve_id=nerve,retained_subdomain=payload,fixed_entities=fixed,free_entity=endpoint,
        external_force_direction=direction.tolist(),runtime_sources=receipts,controller_source=GEYER_HERR_URL,
        experiment=dict(load_n=20,start_s=.1,end_s=.25,duration_s=.4,exchange_step_s=trials['intact']['dt_s'],feedback_exchange_latency_s=trials['intact']['dt_s'],refined_step_s=refined['dt_s']),
        ownership='Independent extracted mechanics experiment; not additive mass in a simultaneous canonical body solve',
        limitations=['TA fiber length is a reference-calibrated stiff-tendon path proxy, not native CE.',
            'Published feedback coefficients transferred to canonical muscle and inferred support; not human stretch-reflex calibration.',
            'Three rigid bone entities and the retained TA tissue carrier, with fixed proximal supports; free calcaneus translation is not articulated ankle dorsiflexion.',
            'Spinal Ia/alpha route delays plus a 1 ms synaptic prior; no gamma controller, identified cortical recruitment law or IBM motor circuit.',
            'Mechanical active/external work is recorded; metabolic cost or BioGears exercise intensity is not inferred from mechanical work.'])
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='data/derived/audits/neuromechanical/experiments')
    args=parser.parse_args();root=Path(__file__).resolve().parents[1];report=run_suite(root,args.output)
    print(json.dumps({'output':report['output'],'passed':report['passed'],'checks':report['checks'],'maximum_path_difference_m':report['maximum_path_difference_m'],'time_refinement_maximum_path_error_m':report['time_refinement_maximum_path_error_m']},indent=2))
    raise SystemExit(0 if report['passed'] else 1)
