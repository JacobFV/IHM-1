"""Validate sampled systemic protocols and their causal contrasts."""
from ihm.assembly.systemic import SystemicConfig, protocol_events, verify_contrasts
from copy import deepcopy
import json
from pathlib import Path
import tempfile
from ihm.native import _sha
from ihm.assembly.systemic_projection import accepted_systemic_sources


def main():
    for kw in ({'protocol': 'unknown'}, {'seconds': 1.03}, {'sample_interval_s': .03},
               {'seconds': 10, 'sample_interval_s': 3}, {'seconds': 0}, {'seconds': True}):
        try:
            SystemicConfig(**kw)
        except ValueError:
            pass
        else:
            raise AssertionError(f'invalid configuration accepted: {kw}')
    events = protocol_events(SystemicConfig(protocol='meal_exercise', seconds=21600, sample_interval_s=30))
    assert events[0]['kind']=='meal' and events[1]['kind']=='exercise'
    assert events[1]['time_s']==1800 and events[2]['time_s']==2400
    assert not verify_contrasts({})['passed']
    fixture={'checks':{'local_store_nonnegativity_passed':True},'native_manifest':{
        key:'shared' for key in ('state_sha256','library_sha256','executable_sha256','patient_identity_input_sha256')},
        'runtime_sources':{'code':'same'},'actions':[],
        'frames':[{'time_s':t,'values':{'lung_volume_ml':3000.,'arterial_co2_mmhg':40.}} for t in (0,30,60,180)]}
    fixture['native_manifest'].update(dependency_sha256={'cdm':'same'},native_step_s=.02)
    apnea=deepcopy(fixture);apnea['actions']=[{'time_s':30,'kind':'apnea','value':1}]
    for frame in apnea['frames'][2:]:
        frame['values']={'lung_volume_ml':2800.,'arterial_co2_mmhg':45.}
    assert not verify_contrasts({'apnea':apnea})['passed']
    assert verify_contrasts({'apnea':apnea,'rest':fixture})['passed']
    bad=deepcopy(apnea);bad['native_manifest']['state_sha256']='other'
    try:verify_contrasts({'apnea':bad,'rest':fixture})
    except ValueError:pass
    else:raise AssertionError('mismatched-state comparison accepted')
    bad=deepcopy(apnea);bad['frames'][0]['values']['lung_volume_ml']=2000
    try:verify_contrasts({'apnea':bad,'rest':fixture})
    except ValueError:pass
    else:raise AssertionError('initial-only difference accepted')
    # A published pass must bind the exact compared bytes, and its acceptance
    # must be reproducible. A boolean in a neighboring JSON file is insufficient.
    with tempfile.TemporaryDirectory() as directory:
        root=Path(directory); inputs={}
        for name,data in [('rest',fixture),('apnea',apnea)]:
            data=deepcopy(data);data['configuration']={'protocol':name}
            path=root/'experiment'/name/'systemic.json';path.parent.mkdir(parents=True)
            path.write_text(json.dumps(data))
            inputs[name]={'path':str(path.relative_to(root)),'sha256':_sha(path)}
        report_path=root/'experiment/contrasts.json'
        report_path.write_text(json.dumps({'passed':True,'inputs':inputs}))
        source=root/inputs['apnea']['path']
        assert len(accepted_systemic_sources(root,source))==3
        source.write_text(source.read_text()+' ')
        try:accepted_systemic_sources(root,source)
        except ValueError:pass
        else:raise AssertionError('changed accepted experiment published')
        inputs['apnea']['sha256']=_sha(source)
        report_path.write_text(json.dumps({'passed':False,'inputs':inputs}))
        try:accepted_systemic_sources(root,source)
        except ValueError:pass
        else:raise AssertionError('failed acceptance published')
    print('PASS systemic protocol validation')


if __name__ == '__main__':
    main()
