"""Strict stance widths, baseline ablation and persistent recorded-state replay."""
import argparse,json,sys
from pathlib import Path
from copy import deepcopy
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import torch,numpy as np
from ihm.native.cortical_stance import load_cortical_stance,PersistentStanceCommands
from ihm.native.stance_lqr import NativeStanceLQR


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact',default='data/runtime/motor-learning/cortical-stance-pilot-20260908/cortical_stance.pt')
    parser.add_argument('--trajectory',default='data/research/locomotion_control/lqr98_discrete_push/trajectory.json')
    parser.add_argument('--teacher',default='data/research/locomotion_control/linearization_3q5uh3o3/discrete/linearization.npz')
    args=parser.parse_args();torch.set_num_threads(1)
    teacher=NativeStanceLQR(args.teacher,dt_s=.01)
    policy,artifact=load_cortical_stance(args.artifact,model_sha256=teacher.model_sha256,dt_s=.01)
    assert policy.state_names==teacher.state_names and policy.muscle_names==teacher.muscle_names
    frames=[f for f in json.loads(Path(args.trajectory).read_text()) if f['time_s']<=2.+1e-9]
    vector=policy.snapshot_vector(frames[0]);initial=policy.state()
    for bad in (vector[:,:-1],torch.full_like(vector,float('nan'))):
        try:policy.advance(bad,initial)
        except ValueError:pass
        else:raise AssertionError('Wrong afferent width or nonfinite feedback accepted')
    bad=deepcopy(frames[0]);bad['muscles'].pop(policy.muscle_names[0])
    try:policy.snapshot_vector(bad)
    except ValueError:pass
    else:raise AssertionError('Missing native muscle accepted')
    try:load_cortical_stance(args.artifact,model_sha256='wrong',dt_s=.01)
    except ValueError:pass
    else:raise AssertionError('Wrong native model identity accepted')
    try:load_cortical_stance(args.artifact,model_sha256=teacher.model_sha256,dt_s=.02)
    except ValueError:pass
    else:raise AssertionError('Wrong sample interval accepted')
    replay={};target=[]
    for frame in frames:
        commands,_=teacher.commands(frame)
        target.append([max(.01,commands[m]) for m in policy.muscle_names])
    target=np.array(target)
    trained_embedding=policy.dyn.embed.detach().clone()
    for arm in ('full','sever','embedding_reset'):
        if arm=='embedding_reset':
            with torch.no_grad():policy.dyn.embed.copy_(artifact['initial_state_dict']['dyn.embed'])
        adapter=PersistentStanceCommands(policy);commands=[]
        for frame in frames:
            output=adapter.commands(frame,sever=arm=='sever')
            assert set(output)==set(policy.muscle_names)
            commands.append([output[m] for m in policy.muscle_names])
        replay[arm]=np.array(commands)
    with torch.no_grad():policy.dyn.embed.copy_(trained_embedding)
    assert np.allclose(replay['sever'],policy.u0.numpy()[None],rtol=0,atol=1e-8)
    # Measure cortical response kinetics separately from endpoint imitation MSE.
    # Both branches start from identical persistent rest history; one receives a
    # real retained perturbed native vector and the other remains at x0.
    strongest=int(np.argmax(np.linalg.norm(target-policy.u0.numpy()[None],axis=1)))
    perturbed=policy.snapshot_vector(frames[strongest]);reference=policy.x0[None]
    with torch.no_grad():
        weights=policy.dyn.edge_weights().detach()
        _,warm=policy.advance(reference,policy.state(),ticks=500,weights=weights)
        response_state=tuple(x.clone() for x in warm);reference_state=tuple(x.clone() for x in warm)
        response_norm=[]
        for tick in range(100):
            response,response_state=policy.advance(perturbed,response_state,ticks=1,weights=weights)
            baseline,reference_state=policy.advance(reference,reference_state,ticks=1,weights=weights)
            response_norm.append(float(torch.linalg.vector_norm(response-baseline)))
    plateau=response_norm[-1]
    kinetics={str(fraction):next(((i+1)*.001 for i,value in enumerate(response_norm) if value>=fraction*plateau),None)
              for fraction in (.1,.5,.9)} if plateau>1e-10 else {}
    nominal_drift=0.
    if policy.reference_normalization:
        nominal_state=policy.state()
        with torch.no_grad():
            for _ in range(4):
                nominal,nominal_state=policy.advance(policy.x0[None],nominal_state,ticks=250,weights=weights)
                nominal_drift=max(nominal_drift,float((nominal-policy.u0).abs().max()))
        assert nominal_drift==0., 'Persistent zero-error decoder drifted from tonic baseline'
    report={'passed_contract_checks':True,'muscles':len(policy.muscle_names),'states':len(policy.state_names),
        'model_sha256':teacher.model_sha256,'artifact_sha256':artifact['artifact_sha256'],
        'recorded_state_persistent_replay_mse':{k:float(np.mean((v-target)**2)) for k,v in replay.items()},
        'reference_normalization':policy.reference_normalization,'zero_error_1s_max_command_drift':nominal_drift,
        'step_response_s_at_fraction_of_100ms_amplitude':kinetics,
        'step_response_100ms_command_norm':plateau,'perturbed_snapshot_time_s':frames[strongest]['time_s'],
        'sever_retains_exact_baseline':True,'native_stance_demonstrated':False,
        'scope':'Open-loop replay of actual native snapshots; cortical history persists; mechanical commands not applied',
        'comparison_limit':'Recorded-state replay is not closed-loop native control; embedding_reset retains trained decoder and restores original materialized association embedding'}
    path=Path(args.artifact).with_name('verification.json');path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
