"""Bounded computed 1-DOF test plant: software causality, not human validation."""
import json
import math
from pathlib import Path
from ihm.assembly.sensorimotor import SensorimotorController
from scripts.verify_sensorimotor import observation
ROOT=Path(__file__).resolve().parents[1]

def run(dt=.002, block=None):
    c=SensorimotorController.from_root(ROOT)
    q=.01;v=0.;activation=0.;frames=[];command=0.
    # Explicit engineering test plant, independent of anatomical mass/calibration.
    mass=1.;stiffness=50.;damping=5.;maximum_force=10.;tau=.01
    for i in range(round(.3/dt)):
        obs=observation(c.time_s,1+q/.1,0.)
        obs['muscles']['tibant_r']['tendon_force_n']=maximum_force*activation
        obs['muscles']['tibant_r']['max_isometric_force_n']=maximum_force
        # Sensor perturbation is actual integrated position, force actual activation.
        kwargs={}
        if block and c.time_s<.15:kwargs[block]=['tibant_r']
        out=c.step(dt,obs,**kwargs)
        activation=command+(activation-command)*math.exp(-dt/tau)
        force=-maximum_force*activation
        v+=dt*(force-stiffness*q-damping*v)/mass;q+=dt*v
        command=out['motor_excitations']['tibant_r']
        if i%max(1,round(.01/dt))==0:
            frames.append({'time_s':c.time_s,'position_m':q,'velocity_m_s':v,'activation':activation,'active_force_n':force,
                           'excitation':out['motor_excitations']['tibant_r'], 'brain_gain':out['descending_gain']['tibant_r']})
    return {'final_position_m':q,'final_velocity_m_s':v,'frames':frames,'model_sha256':c.model_sha256}

def main():
    cases={k:run(block=b) for k,b in [('intact',None),('sensory_block_release','sensory_blocks'),('motor_block_release','motor_blocks')]}
    fine=run(.001)
    early=lambda k:cases[k]['frames'][10]['active_force_n']
    assert abs(early('intact'))>abs(early('sensory_block_release'))>abs(early('motor_block_release'))
    assert cases['motor_block_release']['frames'][-1]['activation']>0
    error=abs(cases['intact']['final_position_m']-fine['final_position_m'])
    assert error<.002
    report={'schema':'ihm.sensorimotor-feedback-probe.v1','passed':True,'cases':cases,
            'dt_2ms_vs_1ms_final_position_difference_m':error,'source':'explicit synthetic 1kg linear test plant',
            'scope':'software feedback, blocks, release and temporal refinement; not native human or gait validation'}
    out=ROOT/'data/derived/audits/sensorimotor_feedback_v1';out.mkdir(parents=True,exist_ok=True)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':True,'early_active_forces_n':{k:early(k) for k in cases},'position_refinement_difference_m':error}))
if __name__=='__main__':main()
