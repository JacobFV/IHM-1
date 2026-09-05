"""Causality, delay, block, checkpoint and real mechanical-feedback regressions."""
from copy import deepcopy
import json
import math
from pathlib import Path
from ihm.assembly.reflexes import BodyReflex, ReflexParameters

ROOT=Path(__file__).resolve().parents[1]

def fixture():
    return dict(id='test-muscle',max_isometric_force_n=100,optimal_fiber_length_m=.1,
        rest_path_length_m=.2,pennation_angle_rad=0,anchors=[
            dict(entity_id='a',point_m=[0,0,0]),dict(entity_id='b',point_m=[.2,0,0])])

def observation(time=0,length=.2):
    return dict(time_s=time,entities={
        'a':dict(translation_m=[0,0,0],rotation_matrix=[[1,0,0],[0,1,0],[0,0,1]],deformation_gradient=[[1,0,0],[0,1,0],[0,0,1]]),
        'b':dict(translation_m=[length-.2,0,0],rotation_matrix=[[1,0,0],[0,1,0],[0,0,1]],deformation_gradient=[[1,0,0],[0,1,0],[0,0,1]])},muscle_forces_n={'test-muscle':0})

def main():
    p=ReflexParameters(length_gain=2,length_offset=1,prestimulation=0,minimum_stimulation=0,
        loop_delay_s=.02,activation_tau_s=.01)
    body=BodyReflex(fixture(),{'a':[0,0,0],'b':[.2,0,0]},p)
    # A 10% fiber extension yields .2 stimulation only after the 20 ms loop delay.
    for i in range(20):
        out=body.step(.001,observation(i*.001,.21))
        assert abs(out['activation'])<1e-14, 'motor activity preceded sensory arrival'
    out=body.step(.001,observation(.02,.21))
    assert math.isclose(out['stimulation'],.2,abs_tol=1e-12)
    assert math.isclose(out['activation'],.2*(1-math.exp(-.1)),abs_tol=1e-12)
    assert out['sensor']['path_length_m']==.21
    # Block removes in-flight sensory samples and decays existing activation.
    previous=out['activation']
    out=body.step(.001,observation(.021,.21),motor_block=True)
    assert math.isclose(out['activation'],previous*math.exp(-.1),abs_tol=1e-12)
    out=body.step(.001,observation(.022,.21),sensory_block=True)
    assert out['stimulation']==0 and out['pending_samples']==0
    for i in range(20):
        out=body.step(.001,observation(.023+i*.001,.21))
        assert out['stimulation']==0 or i==19
    out=body.step(.001,observation(.043,.21),controller_enabled=False,descending_drive=.1)
    assert out['stimulation']==.1, 'controller removal must preserve explicit descending input'
    checkpoint=body.checkpoint()
    a=body.step(.001,observation(.044,.21),descending_gain=.5)
    body.restore(checkpoint)
    b=body.step(.001,observation(.044,.21),descending_gain=.5)
    assert a==b and body.checkpoint()!=checkpoint
    assert math.isclose(a['stimulation'],.1,abs_tol=1e-12), 'descending gain did not scale sensory feedback'
    before=body.checkpoint()
    for kwargs in ({'descending_drive':float('nan')},{'motor_block':1},{'descending_gain':-1}):
        try:body.step(.001,observation(.045,.21),**kwargs)
        except ValueError:pass
        else:raise AssertionError('invalid input accepted')
        assert body.checkpoint()==before, 'invalid input mutated controller'
    bad=observation(.045,.21);del bad['entities']['b']['rotation_matrix']
    try:body.step(.001,bad)
    except ValueError:pass
    else:raise AssertionError('missing mechanical feedback accepted')
    assert body.checkpoint()==before
    bad=deepcopy(checkpoint);bad['model_sha256']='wrong'
    try:body.restore(bad)
    except ValueError:pass
    else:raise AssertionError('cross-model checkpoint accepted')
    assert body.checkpoint()==before
    # Full attachment transforms preserve length under a common rigid motion.
    rigid=BodyReflex(fixture(),{'a':[0,0,0],'b':[.2,0,0]},p)
    rotated=observation();rotated['entities']['a']['rotation_matrix']=[[0,-1,0],[1,0,0],[0,0,1]]
    rotated['entities']['b']['rotation_matrix']=[[0,-1,0],[1,0,0],[0,0,1]]
    rotated['entities']['a']['translation_m']=[1,2,3]
    rotated['entities']['b']['translation_m']=[.8,2.2,3]
    assert math.isclose(rigid.observe(rotated)['path_length_m'],.2,abs_tol=1e-12)
    # Constant held sensory input gives the same activation after time partitioning.
    whole=BodyReflex(fixture(),{'a':[0,0,0],'b':[.2,0,0]},p)
    split=BodyReflex(fixture(),{'a':[0,0,0],'b':[.2,0,0]},p)
    a=whole.step(.1,observation(0,.21))
    for i in range(100):b=split.step(.001,observation(i*.001,.21))
    assert math.isclose(a['activation'],b['activation'],abs_tol=1e-12)
    # Unlike the unit sensor fixture, this trial invokes the existing mechanics solver.
    from scripts.build_neuromechanical_experiments import run_suite
    report=run_suite(ROOT,ROOT/'data/derived/audits/neuromechanical/verification')
    assert report['passed'],report['checks']
    print(json.dumps({'passed':True,'checks':report['checks']},indent=2))

if __name__=='__main__':main()
