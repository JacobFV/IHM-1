"""Validate coupling inputs, then optionally execute short source/feedback parity."""
import argparse
import json
from pathlib import Path
import tempfile
from ihm.native.session import NativeSession,SessionConfig
from ihm.native.coupled_session import CoupledNativeSession


def contracts():
    session=object.__new__(CoupledNativeSession);calls=[]
    command=['/held/engine','patient','/held/state','50']
    bounded=session.launch_command(command)
    assert bounded[:3]==['/usr/bin/prlimit','--as=4294967296','--'] and bounded[-4:]==command
    session._command=lambda *args:calls.append(args)
    for value in [True,float('nan'),float('inf'),5001,-5001]:
        try:session.respiratory_load(value)
        except ValueError:pass
        else:raise AssertionError('Invalid pressure accepted')
    assert not calls
    session.respiratory_load(100);session.respiratory_load(0)
    assert calls==[('respiratory_load',100),('respiratory_load',0)]


def native():
    root=Path(__file__).resolve().parents[1]
    state=root/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml'
    output=Path(tempfile.mkdtemp(prefix='native-mechanical-feedback-',dir=root/'data/derived/audits'))
    config=SessionConfig(state_path=state,engine_variant='whole_body_integrity_depletion',horizon_s=1)
    traces={}
    for label,cls,load in [('reference',NativeSession,0),('coupled_zero',CoupledNativeSession,0),('loaded',CoupledNativeSession,98.0665)]:
        with cls(config,output/label) as session:
            initial=session.snapshot()
            if cls is CoupledNativeSession:session.respiratory_load(load)
            frames=[session.step(.02) for _ in range(25)]
            if cls is CoupledNativeSession:session.respiratory_load(0)
            frames.extend(session.step(.02) for _ in range(25))
            traces[label]={'initial':initial,'frames':frames}
    shared=set(traces['reference']['frames'][0]['values'])
    for a,b in zip(traces['reference']['frames'],traces['coupled_zero']['frames']):
        for key in shared:
            assert a['values'][key]==b['values'][key],(key,a['values'][key],b['values'][key])
    pressure_error=max(abs(f['values']['coupling.applied_driver_pa']-f['values']['coupling.generated_driver_pa']-98.0665) for f in traces['loaded']['frames'][:25])
    assert pressure_error<1e-10
    assert all(f['values']['coupling.external_pressure_pa']==0 for f in traces['loaded']['frames'][25:])
    delta=max(abs(a['values']['lung_volume_ml']-b['values']['lung_volume_ml']) for a,b in zip(traces['loaded']['frames'],traces['coupled_zero']['frames']))
    assert delta>1e-5
    result={'passed':True,'native_ticks_per_condition':50,'zero_boundary_exact_source_port_parity':True,
        'pressure_residual_pa':pressure_error,'maximum_volume_change_ml':delta,
        'scope':'One-second native external-load feedback; no long-horizon or physiological calibration claim',
        'directory':str(output.relative_to(root))}
    (output/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--native',action='store_true');args=p.parse_args()
    contracts()
    if args.native:native()
    else:print('Coupled native input validation passed; no native process started.')
