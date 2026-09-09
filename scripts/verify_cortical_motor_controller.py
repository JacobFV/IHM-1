"""Verify trained actual E/I controller protocol, ablations and transactional state."""
from pathlib import Path
from copy import deepcopy
import json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import torch
from ihm.native.cortical_motor_controller import IBMCorticalAnkleController
from ihm.assembly.sensorimotor_catalog import native_muscle_catalog


def main():
    torch.set_num_threads(1)
    catalog=native_muscle_catalog(ROOT)
    c=IBMCorticalAnkleController.from_root(ROOT,muscle_catalog=catalog)
    observation={'time_s':0.,'foot_contact_force_n':{'r':0.,'l':0.},
        'joints':{'ankle_angle_r':{'value':0.,'speed':0.}},
        'muscles':{m:{'sensor_basis':'verification fixture on actual native port schema',
            'optimal_fiber_length_m':1.,'max_isometric_force_n':100.,
            'fiber_length_m':1.1,'tendon_force_n':20.} for m in c.muscles}}
    initial=c.checkpoint();first=c.step(.02,observation)
    c.restore(initial);assert first==c.step(.02,observation)
    observation['time_s']=.02;second=c.step(.02,observation)
    assert all(v>0 for v in second['arc_max'].values())
    assert second['controller']['cortical_motor_output_active']
    # Block after afference has entered every delay queue. Different future
    # positions/velocities must produce identical commands from the same cortical
    # history, while all afferent-driven cord arcs are flushed immediately.
    delayed=c.checkpoint();observation['time_s']=.04
    forward=deepcopy(observation);backward=deepcopy(observation)
    forward['joints']['ankle_angle_r']={'value':.2,'speed':.5}
    backward['joints']['ankle_angle_r']={'value':-.2,'speed':-.5}
    blocked_forward=c.step(.02,forward,sensory_blocks=c.muscles)
    blocked_endpoint=c.checkpoint()
    for name in ('stretch','autogenic','reciprocal'):
        assert blocked_forward['arc_max'][name]==0., name+' retained queued afference'
        assert all(not any(sample) for sample in blocked_endpoint['inner']['delays'][name])
    c.restore(delayed)
    blocked_backward=c.step(.02,backward,sensory_blocks=c.muscles)
    assert blocked_forward['cortical_commands']==blocked_backward['cortical_commands']
    assert blocked_forward['motor_excitations']==blocked_backward['motor_excitations']
    assert blocked_endpoint==c.checkpoint()
    # Previously integrated cortical memory and Renshaw recurrence intentionally
    # persist across a sensory block. The block is not a brain reset.
    assert blocked_forward['arc_max']['renshaw']>0.
    c.restore(delayed)
    unblocked=c.step(.02,forward)
    assert unblocked['cortical_commands']!=blocked_forward['cortical_commands']
    valid=c.checkpoint()
    malformed=[]
    bad=deepcopy(valid);bad['excitations']['tibant_r']=.987;malformed.append(bad)
    bad=deepcopy(valid);del bad['excitations'];malformed.append(bad)
    bad=deepcopy(valid);bad['excitations']=[];malformed.append(bad)
    bad=deepcopy(valid);bad['inner']['cortical_state'][0][0][0]=float('nan');malformed.append(bad)
    bad=deepcopy(valid);bad['inner']['cortical_state'][0]=[];malformed.append(bad)
    bad=deepcopy(valid);del bad['inner']['delays']['stretch'];malformed.append(bad)
    bad=deepcopy(valid);bad['inner']['delays']['reciprocal'][0]=[];malformed.append(bad)
    bad=deepcopy(valid);bad['inner']['time_s']+=.0005;malformed.append(bad)
    bad=deepcopy(valid);bad['inner']['alpha'][0]=1.01;malformed.append(bad)
    for bad in malformed:
        try:c.restore(bad)
        except ValueError:pass
        else:raise AssertionError('Malformed trained checkpoint accepted')
        assert c.checkpoint()==valid, 'Failed restore mutated accepted state'
    # Late failures after E/I and cord updates must also rewind every owner.
    c.restore(delayed)
    original_advance=c.policy.advance
    calls=[0]
    def fail_after_cortical_progress(*args,**kwargs):
        calls[0]+=1
        if calls[0]==3:raise RuntimeError('injected late integrator failure')
        return original_advance(*args,**kwargs)
    c.policy.advance=fail_after_cortical_progress
    try:
        try:c.step(.02,forward)
        except RuntimeError:pass
        else:raise AssertionError('Injected failure missing')
        assert c.checkpoint()==delayed
    finally:c.policy.advance=original_advance
    c.restore(initial);observation['time_s']=0.
    for argument in ({'descending':{'wrong':1.}},{'additional_sensory_inputs_hz':{'wrong':1.}},
                     {'physiology':{'oxygen_saturation':float('nan')}}):
        try:c.step(.02,observation,**argument)
        except ValueError:pass
        else:raise AssertionError('Invalid input accepted')
        assert c.checkpoint()==initial
    blocked=c.step(.02,observation,motor_blocks=c.muscles)
    assert not any(blocked['motor_excitations'].values())
    c.restore(initial)
    ordinary=c.step(.1,observation)
    c.restore(initial)
    hypoxic=c.step(.1,observation,physiology={'oxygen_saturation':.2})
    assert ordinary['cortical_commands']!=hypoxic['cortical_commands']
    sever=IBMCorticalAnkleController.from_root(ROOT,muscle_catalog=catalog,sever=True,no_cord=True)
    disconnected=sever.step(.1,observation)
    assert not any(disconnected['motor_excitations'].values())
    print(json.dumps({'passed':True,'checks':['persistent replay','all named cord arcs','strict inputs',
        'motor blocks','sensory block flush and input isolation','atomic malformed restore','late-step rollback','physiology causality','severed cortical motor silence'],
        'trained_tibant_command':ordinary['cortical_commands']['tibant_r'],
        'artifact_sha256':c.artifact['artifact_sha256'],'scope':'Native-port fixtures; native body tracking is a separate evaluation'},indent=2))

if __name__=='__main__':main()
