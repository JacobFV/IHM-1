#!/usr/bin/env python3
"""Paired native baseline, exercise recovery and hemorrhage/saline verification."""
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native import *
base=BASE/'data/derived/physiology'
state=base/'native_baseline_v2/states/native_stabilized.xml'
for name,events in [('native_exercise_low',(Intervention(30,'exercise',.05),Intervention(60,'exercise',0))),('native_fluid_recovery',(Intervention(30,'hemorrhage',60),Intervention(90,'hemorrhage',0),Intervention(90,'saline',60),Intervention(150,'saline',0)))]:
    out=base/name
    if not (out/'execution.json').exists():run_native(NativeConfig(seconds=180,state_path=str(state),interventions=events),out)
    print(name,load_trajectory(out)['values']['BloodVolume(mL)'][-1],flush=True)
# Comparisons use the reloaded control, preserving the reload discrepancy separately.
control_dir=base/'native_state_baseline'
if not (control_dir/'execution.json').exists():run_native(NativeConfig(seconds=180,state_path=str(state)),control_dir)
control=load_trajectory(control_dir);fresh=load_trajectory(base/'native_baseline_v2')
exercise=load_trajectory(base/'native_exercise_low');fluid=load_trajectory(base/'native_fluid_recovery')
def delta(run,col,lo,hi):
    values=[x-y for t,x,y in zip(control['time_s'],run['values'][col],control['values'][col]) if lo<t<=hi]
    return sum(values)/len(values)
report={'evidence_kind':'source_simulation','experimental_validation':False,
 'baseline_relative_blood_volume_range':(max(control['values']['BloodVolume(mL)'])-min(control['values']['BloodVolume(mL)']))/control['values']['BloodVolume(mL)'][0],
 'state_reload_max_absolute_error':{c:max(abs(x-y) for x,y in zip(control['values'][c],fresh['values'][c])) for c in control['values']},
 'exercise_hr_delta_bpm':delta(exercise,'HeartRate(1/min)',55,60),
 'exercise_power_delta_W':delta(exercise,'TotalMetabolicRate(W)',55,60),
 'recovery_power_delta_W':delta(exercise,'TotalMetabolicRate(W)',170,180),
 'hemorrhage_volume_delta_mL':delta(fluid,'BloodVolume(mL)',85,90),
 'after_saline_volume_delta_mL':delta(fluid,'BloodVolume(mL)',145,150),
 'scope':'Blood-volume stability and perturbation accounting only, not proof of full whole-body mass conservation. Hemorrhage rate is pressure dependent. Recovery is partial.',
 'failed_experiment':{'path':'native_exercise_recovery','exercise_intensity':.3,'failure':'irreversible hypercapnia at 130.32 s; exit 3 rejected; preserved'}}
assert report['baseline_relative_blood_volume_range']<.002
assert report['state_reload_max_absolute_error']['BloodVolume(mL)']<.1
assert report['exercise_hr_delta_bpm']>1
assert report['exercise_power_delta_W']>.5
assert abs(report['recovery_power_delta_W'])<.2*report['exercise_power_delta_W']
assert report['hemorrhage_volume_delta_mL']<-30
assert 45<report['after_saline_volume_delta_mL']-report['hemorrhage_volume_delta_mL']<70
report['checks_passed']=True
(base/'native_experiment_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
