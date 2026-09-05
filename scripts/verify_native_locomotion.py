"""Behavioral acceptance for the native articulated/contact forward plant."""
from pathlib import Path
import sys, tempfile, json
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np

def verify():
    try:
        from ihm.native.locomotion import LocomotionConfig, run_locomotion
    except ImportError:
        raise AssertionError('Native time-integrated locomotion API is missing') from None
    for kwargs in ({'seconds':float('nan')}, {'excitation':True}, {'friction_scale':-1}, {'maximum_step_s':0}):
        try: LocomotionConfig(**kwargs)
        except ValueError: pass
        else: raise AssertionError(f'Unsafe configuration accepted: {kwargs}')
    root=Path(tempfile.mkdtemp(prefix='locomotion-',dir=Path('data/derived')))
    results={}
    for name,kwargs in [('baseline',{}),('velocity',{'pelvis_speed_delta_m_s':.05}),('contact',{'stiffness_scale':.5}),('friction',{'friction_scale':.25}),('excitation',{'excitation':.1}),('feedback',{'reflex_gain':.5}),('refined',{'maximum_step_s':.0005,'accuracy':1e-7}),('finest',{'maximum_step_s':.00025,'accuracy':1e-8})]:
        data=run_locomotion(LocomotionConfig(seconds=.1,**kwargs),root/name)
        assert data.get('termination',{}).get('completed') is True, 'Explicit native termination receipt missing'
        frames=data['frames']; assert len(frames)==21
        assert abs(frames[-1]['time_s']-.1)<1e-10
        assert data['model']['muscle_count']==80 and data['model']['body_count']==22
        assert data['model']['prescribed_coordinate_count']==0
        assert data['model']['contact_count']==12
        assert abs(data['model']['mass_kg']-85.26984854173146)<1e-8
        assert max(f['constraint_position_error'] for f in frames)<1e-5
        assert max(f['constraint_velocity_error'] for f in frames)<1e-4
        assert np.linalg.norm(np.array(frames[-1]['q'])-frames[0]['q'])>1e-3
        assert max(np.linalg.norm(f['contact_force_N']) for f in frames)>1
        assert all('momentum_balance_residual_N' in f for f in frames), 'Native force/momentum audit missing'
        assert max(np.linalg.norm(f['momentum_balance_residual_N']) for f in frames)<1e-5
        assert max(np.linalg.norm(c['force_pair_residual_N']) for f in frames for c in f['contacts'])<1e-8
        assert all(.009<=m['activation']<=1.001 for f in frames for m in f['muscles'])
        results[name]=data
    b=results['baseline']['frames'][-1]
    metrics={}
    for name in ('velocity','contact','friction','excitation','feedback','refined','finest'):
        metrics[name+'_q_difference']=float(np.linalg.norm(np.array(results[name]['frames'][-1]['q'])-b['q']))
        if name not in ('refined','finest'):assert metrics[name+'_q_difference']>1e-6
    assert metrics['refined_q_difference']<1e-3
    metrics['last_refinement_q_difference']=float(np.linalg.norm(np.array(results['refined']['frames'][-1]['q'])-results['finest']['frames'][-1]['q']))
    assert metrics['last_refinement_q_difference']<metrics['refined_q_difference']
    fall=run_locomotion(LocomotionConfig(seconds=1),root/'fall')
    assert fall['termination']['completed'] is False
    assert fall['termination']['reason']=='upright_plant_envelope_exceeded'
    assert fall['frames'][-1]['time_s']<1
    metrics['fall_time_s']=fall['frames'][-1]['time_s']
    metrics.update(output_dir=str(root), status='passed', horizon_s=.1)
    (root/'verification.json').write_text(json.dumps(metrics,indent=2)+'\n')
    print(json.dumps(metrics,indent=2))
    return metrics
if __name__=='__main__':verify()
