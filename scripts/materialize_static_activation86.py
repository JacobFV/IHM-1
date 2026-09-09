#!/usr/bin/env python3
"""Opt-in native86 initial muscle activation candidate from static force fitting."""
from pathlib import Path
import hashlib,json,re,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.assembly.embodied import _prepare_mechanical_registration
from ihm.native.mechanical_stream import NativeMechanicalStream
ROOT=Path(__file__).resolve().parents[1]
BASE='data/derived/mechanics/initial_stance86/registration.json'
REPORT='data/research/locomotion_control/native_static_activation_solve_ah5yuur3/report.json'
OUT=ROOT/'data/derived/mechanics/initial_stance86_activation'
sha=lambda raw:hashlib.sha256(raw).hexdigest()

def main():
    if OUT.exists():raise ValueError('Refusing to overwrite activation candidate')
    _,base,rows=_prepare_mechanical_registration(ROOT,BASE)
    fit=json.loads((ROOT/REPORT).read_text());activations=fit['activations']
    if set(activations)!={r['id'] for r in rows}:raise ValueError('Exact native activation catalog required')
    raw=(ROOT/base['model_path']).read_bytes();model=raw
    for name,value in activations.items():
        pattern=rb'(<(?P<kind>Millard2012EquilibriumMuscle|Thelen2003Muscle) name="'+re.escape(name.encode())+rb'">)(?P<body>.*?)(</(?P=kind)>)'
        def replace(m):
            body=m['body'];v=format(value,'.17g').encode()
            if b'<default_activation>' in body:body=re.sub(rb'(<default_activation>).*?(</default_activation>)',lambda n:n[1]+v+n[2],body,count=1,flags=re.S)
            else:body=b'\n<default_activation>'+v+b'</default_activation>\n'+body
            return m[1]+body+m[4]
        model,count=re.subn(pattern,replace,model,count=1,flags=re.S)
        if count!=1:raise ValueError('Expected exactly one muscle')
    OUT.mkdir(parents=True);path=OUT/'model.osim';path.write_bytes(model)
    registration=dict(base);registration['sources']=dict(base['sources']);registration['sources'].update({BASE:sha((ROOT/BASE).read_bytes()),REPORT:sha((ROOT/REPORT).read_bytes()),base['model_path']:base['model_sha256'],str(Path(__file__).relative_to(ROOT)):sha(Path(__file__).read_bytes())})
    registration['model_path']=str(path.relative_to(ROOT));registration['model_sha256']=sha(model);registration['equilibrium_excitations']=activations
    registration['native_acceptance_complete']=False;registration['static_activation_fit_report']=REPORT
    regpath=OUT/'registration.json';regpath.write_text(json.dumps(registration,indent=2)+'\n')
    _prepare_mechanical_registration(ROOT,str(regpath.relative_to(ROOT)))
    native=NativeMechanicalStream(ROOT,OUT/'native_initialization',environment='upright',target_mass_kg=70,augmented_registration=str(regpath.relative_to(ROOT)))
    try:
        state=native.snapshot()
        error=max(abs(state['muscles'][n]['activation']-v) for n,v in activations.items())
        if error>1e-12:raise ValueError('Native initialized activation differs')
        static=native.evaluate_static_pose({'pelvis_ty':state['coordinates']['pelvis_ty']['value']})
        residual=dict(zip(static['mobility_coordinate_names'],static['constrained_zero_acceleration_residual_mobility_force']))
        norm=sum(residual[n]**2 for n in fit['selected_coordinates'])**.5
        receipt={'schema':'ihm.static-activation86-initialization.v1','activation_default_max_error':error,'selected_residual_norm_nm':norm,'selected_coordinates':fit['selected_coordinates'],'excluded_residual':{n:v for n,v in residual.items() if n not in fit['selected_coordinates']},'physical_time_advanced_s':0,'full_static_equilibrium_claimed':False,'sustained_stance_demonstrated':False}
        (OUT/'initial_snapshot.json').write_text(json.dumps(state,allow_nan=False));report=OUT/'initial_acceptance.json';report.write_text(json.dumps(receipt,indent=2)+'\n')
    finally:native.close()
    registration['sources'][str(report.relative_to(ROOT))]=sha(report.read_bytes());registration['initial_activation_acceptance']=str(report.relative_to(ROOT));regpath.write_text(json.dumps(registration,indent=2)+'\n')
    print(json.dumps({'registration':str(regpath.relative_to(ROOT)),**receipt},indent=2))
if __name__=='__main__':main()
