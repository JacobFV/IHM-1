"""Validate sampled systemic protocols and their causal contrasts."""
from ihm.assembly.systemic import SystemicConfig, protocol_events, verify_contrasts
from copy import deepcopy
import json
from pathlib import Path
import tempfile
from unittest.mock import patch
from ihm.native import _sha
from ihm.assembly.systemic_projection import accepted_systemic_sources,require_generic_thermal_domain,require_generic_homeostasis
from ihm.assembly.systemic_evidence import freeze_sources


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
        'configuration':{'protocol':'rest','seconds':180,'sample_interval_s':30},
        'frames':[{'time_s':t,'values':{'lung_volume_ml':3000.,'arterial_co2_mmhg':40.,'core_temperature_c':37.,
                  'arterial_ph':7.4,'Aorta.Glucose.concentration_mg_per_dl':85.}} for t in range(0,181,30)]}
    for frame in fixture['frames']:
        frame['values'].update({key:1. for key in ['stomach_'+x for x in ('carbohydrate_g','protein_g','fat_g','sodium_g','calcium_mg','water_ml')]+['liver_glycogen_g','muscle_glycogen_g','stored_protein_g','stored_fat_g']})
    fixture['native_manifest'].update(dependency_sha256={'cdm':'same'},native_step_s=.02)
    apnea=deepcopy(fixture);apnea['configuration']['protocol']='apnea'
    apnea['actions']=protocol_events(SystemicConfig(protocol='apnea',seconds=180,sample_interval_s=30))
    for frame in apnea['frames'][2:]:
        frame['values'].update(lung_volume_ml=2800.,arterial_co2_mmhg=45.)
    assert not verify_contrasts({'apnea':apnea})['passed']
    assert verify_contrasts({'apnea':apnea,'rest':fixture})['passed']
    for label,mutate in [('negative hidden by cached pass',lambda d:d['frames'][-1]['values'].update(stomach_water_ml=-1000)),
                         ('wrong action kind',lambda d:d['actions'][0].update(kind='exercise')),
                         ('missing required store',lambda d:d['frames'][-1]['values'].pop('stomach_water_ml')),
                         ('invalid clock',lambda d:d['frames'][-1].update(time_s=181))]:
        bad=deepcopy(apnea);mutate(bad)
        try:verify_contrasts({'apnea':bad,'rest':fixture})
        except ValueError:pass
        else:raise AssertionError(label)
    bad=deepcopy(apnea);bad['native_manifest']['state_sha256']='other'
    try:verify_contrasts({'apnea':bad,'rest':fixture})
    except ValueError:pass
    else:raise AssertionError('mismatched-state comparison accepted')
    bad=deepcopy(apnea);bad['frames'][0]['values']['lung_volume_ml']=2000
    try:verify_contrasts({'apnea':bad,'rest':fixture})
    except ValueError:pass
    else:raise AssertionError('initial-only difference accepted')
    for temperature in (24.,39.,None,float('nan')):
        cold=deepcopy(fixture);cold['frames'][-1]['values']['core_temperature_c']=temperature
        try:require_generic_thermal_domain(cold)
        except ValueError:pass
        else:raise AssertionError('out-of-domain generic physiology published')
    require_generic_homeostasis(fixture)
    for key,values in [('arterial_ph',(7.29,7.50,None,float('nan'))),
                       ('Aorta.Glucose.concentration_mg_per_dl',(32.9,69.99,None,float('inf')))]:
        for value in values:
            bad=deepcopy(fixture);bad['frames'][-1]['values'][key]=value
            try:require_generic_homeostasis(bad)
            except ValueError:pass
            else:raise AssertionError('Abnormal or absent homeostatic observation accepted: '+key)
    # Apnea deliberately perturbs respiratory acid-base state; it is not a
    # normal-rest claim. Its glucose and finite-pH requirements still apply.
    perturbation=deepcopy(apnea);perturbation['frames'][-1]['values']['arterial_ph']=7.48
    require_generic_homeostasis(perturbation)
    perturbation['frames'][-1]['values']['arterial_ph']=None
    try:require_generic_homeostasis(perturbation)
    except ValueError:pass
    else:raise AssertionError('Apnea accepted without observed pH')
    # A published pass must bind the exact compared bytes, and its acceptance
    # must be reproducible. A boolean in a neighboring JSON file is insufficient.
    # Native environment identity has its own real-archive rejection suite;
    # isolate paired acceptance here from ELF/resource acquisition.
    with tempfile.TemporaryDirectory() as directory, patch('ihm.assembly.systemic_projection.resolve_native_environment',return_value={}), patch('ihm.assembly.systemic_projection.resolve_systemic_execution',return_value={}):
        root=Path(directory); inputs={}
        solver=root/'solver.py';solver.write_text('retained numerical source')
        for name,data in [('rest',fixture),('apnea',apnea)]:
            data=deepcopy(data)
            data['runtime_sources']={'solver.py':_sha(solver)}
            path=root/'experiment'/name/'systemic.json';path.parent.mkdir(parents=True)
            path.write_text(json.dumps(data))
            freeze_sources(root,path.parent,data['runtime_sources'])
            inputs[name]={'path':str(path.relative_to(root)),'sha256':_sha(path)}
        report_path=root/'experiment/contrasts.json'
        report_path.write_text(json.dumps({'passed':True,'inputs':inputs}))
        source=root/inputs['apnea']['path']
        resolved=accepted_systemic_sources(root,source)
        control_copy=root/'experiment/rest/inputs/solver.py'
        assert str(control_copy.relative_to(root)) in resolved
        original=control_copy.read_bytes();control_copy.write_text('changed control solver')
        try:accepted_systemic_sources(root,source)
        except ValueError:pass
        else:raise AssertionError('corrupted control evidence accepted')
        control_copy.write_bytes(original)
        extra=deepcopy(data);extra['configuration']['protocol']='hydration'
        extra['native_manifest']['library_sha256']='unrelated'
        extra_path=root/'experiment/hydration/systemic.json';extra_path.parent.mkdir()
        extra_path.write_text(json.dumps(extra))
        extra_inputs={**inputs,'hydration':{'path':str(extra_path.relative_to(root)),'sha256':_sha(extra_path)}}
        report_path.write_text(json.dumps({'passed':True,'inputs':extra_inputs}))
        try:accepted_systemic_sources(root,extra_path)
        except ValueError:pass
        else:raise AssertionError('unpaired extra experiment published')
        extra_inputs['../../escaped']=extra_inputs.pop('hydration')
        report_path.write_text(json.dumps({'passed':True,'inputs':extra_inputs}))
        try:accepted_systemic_sources(root,extra_path)
        except ValueError:pass
        else:raise AssertionError('unsafe protocol accepted')
        report_path.write_text(json.dumps({'passed':True,'inputs':inputs}))
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
