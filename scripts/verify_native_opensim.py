"""Validate native mechanics arguments and original engine numerical results."""
import sys,unittest,json
import numpy as np
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.opensim_backend import OpenSimConfig
class Contracts(unittest.TestCase):
    def test_bounds(self):
        for kw in ({'activation':float('nan')},{'activation':2},{'delta_rad':.5},{'coordinate':'../x'},{'delta_rad':True}):
            with self.assertRaises(ValueError):OpenSimConfig(**kw)
def verify_engine(variant="upstream"):
    from ihm.native.opensim_backend import BASE,run_opensim
    from ihm.spatial.opensim import load_model
    root=BASE/'data/derived/opensim'/('native' if variant=='upstream' else 'native_corrected');runs={}
    configs={'baseline':OpenSimConfig(),'plus':OpenSimConfig(delta_rad=.1001),'minus':OpenSimConfig(delta_rad=.0999),'knee_bend':OpenSimConfig(delta_rad=.1),'activation':OpenSimConfig(activation=.1)}
    for name,config in configs.items():
        config=OpenSimConfig(coordinate=config.coordinate,delta_rad=config.delta_rad,activation=config.activation,engine_variant=variant)
        out=root/name
        runs[name]=json.loads((out/'summary.json').read_text()) if (out/'summary.json').exists() else run_opensim(config,out)
    baseline=runs['baseline'];muscles={m['name']:m for m in baseline['muscles']}
    assert len(muscles)==80 and baseline['equilibrated']
    xml=load_model(BASE/'data/raw/anatomy/opensim-models/source/Models/Rajagopal/Rajagopal2016.osim')
    transform_error=max(float(np.max(abs(np.asarray(t)-xml['frames'][name]))) for name,t in baseline['body_transforms_ground'].items())
    assert transform_error<1e-9
    plus={m['name']:m for m in runs['plus']['muscles']};minus={m['name']:m for m in runs['minus']['muscles']}
    bend={m['name']:m for m in runs['knee_bend']['muscles']}
    errors={n:abs(m['moment_arms_m']['knee_angle_r']+(plus[n]['length_m']-minus[n]['length_m'])/.0002) for n,m in bend.items()}
    equilibrium={n:abs(m['tendon_force_N']-m['fiber_force_along_tendon_N']) for n,m in muscles.items()}
    report={'source_simulation':True,'experimental_validation':False,'body_transform_max_error':transform_error,'moment_arm_finite_difference_errors_m':errors,'equilibrium_force_residual_N':equilibrium,'active_wrap_segments':sum(m['active_wrap_segments'] for m in muscles.values()),'muscle_count':80}
    (root/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if not isinstance(v,dict)},indent=2))
    print('max moment arm difference',max(errors.values()),'max equilibrium residual',max(equilibrium.values()))
    assert report['active_wrap_segments']>0
    failed={n:e for n,e in errors.items() if e>=.002}
    if variant=='upstream':assert set(failed)=={'gasmed_r','gaslat_r'}, 'Unexpected native moment-arm discrepancies'
    else:
        assert max(errors.values())<1e-4
        assert max(m['max_wrap_location_cache_error_m'] for m in runs['knee_bend']['muscles'])<1e-10
    report['moment_arm_consistency_passed']=not failed
    report['unverified_moment_arms']=failed
    report['independent_checks_passed']=True
    (root/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    assert max(equilibrium.values())<.1
    assert any(abs(m['length_m']-muscles[m['name']]['length_m'])>.001 for m in runs['knee_bend']['muscles'])
    assert sum(m['tendon_force_N'] for m in runs['activation']['muscles'])>sum(m['tendon_force_N'] for m in muscles.values())
    for run in runs.values():
        for muscle in run['muscles']:
            muscle['moment_arm_check']={'coordinate':'knee_angle_r','at_coordinate_rad':.1,'absolute_error_m':errors[muscle['name']],'verified_at_test_pose':muscle['name'] not in failed}
    for name,run in runs.items():(root/name/'summary.json').write_text(json.dumps(run,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':
    if '--engine' in sys.argv:sys.argv.remove('--engine');verify_engine()
    if '--corrected' in sys.argv:sys.argv.remove('--corrected');verify_engine('wrap_8_0.0005_cache')
    unittest.main()
