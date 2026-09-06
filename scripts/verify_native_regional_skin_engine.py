#!/usr/bin/env python3
"""Standalone actual-native state-load/split acceptance with bounded engine steps.

Default prepares artifacts only. --run requires the root's heavy execution slot.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import tempfile
import time

from verify_native_regional_skin import ROOT,RUNTIME,SOURCE,sha
VARIANT=RUNTIME/'variants/whole_body_integrity_regional_skin_graph_v2'
LIBRARY_PIN='bc91cbab829c04bfa2df7490433af5bee752e7514332ac715e9f0bd11a17650a'
MANIFEST_PIN='52de403e3dab2682f774755b7bf1370711bc5ab218244c581460b84a82aad130'
ENGINE=SOURCE/'projects/biogears/libBiogears/src/engine/Controller/BioGearsEngine.cpp'
ENGINE_PIN='5c29a09c536a876d569069625bbe4ef2c9e375e3ba12c48c64b94ceed49cbc95'


def step_source(text):
    signature='bool BioGearsEngine::AdvanceModelTime(bool appendDataTrack)'
    start=text.index(signature);end=text.index('\n//-------------------------------------------------------------------------------',start)
    method=text[start:end].replace(signature,'bool AdvanceModelTime(bool appendDataTrack=false) override')
    for old,new in [('  PreProcess();','  PreProcess();\n  if(regional)regional->after_preprocess();'),
                    ('  PostProcess();','  PostProcess();\n  if(regional)regional->after_postprocess();')]:
        if method.count(old)!=1:raise ValueError('Unexpected native lifecycle anchor')
        method=method.replace(old,new)
    return method+'\n'


def cap(memory,cpu):
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(memory,memory));resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu))
    return limits


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',action='store_true');args=parser.parse_args()
    assert sha(ENGINE)==ENGINE_PIN
    assert sha(VARIANT/'manifest.json')==MANIFEST_PIN and sha(VARIANT/'libbiogears.so.8.0.0')==LIBRARY_PIN
    vm=json.loads((VARIANT/'manifest.json').read_text())
    header=ROOT/'scripts/native_regional_skin.h'
    assert sha(header)==vm['regional_skin_header_sha256'], 'Built variant regional header changed'
    assert sha(ROOT/'scripts/native_regional_species.h')==vm['regional_species_header_sha256'], 'Built sweat helper changed'
    source_receipt=ROOT/'data/derived/systemic/exertion_v3/exercise/native/manifest.json'
    state=Path(json.loads(source_receipt.read_text())['configuration']['state_path'])
    if not state.is_file():raise ValueError('Missing retained stabilized state')
    out=Path(tempfile.mkdtemp(prefix='native-regional-engine-',dir=ROOT/'data/derived/audits'))
    files=[header,ROOT/'scripts/native_regional_species.h',ROOT/'scripts/native_regional_sweat_fixture.cpp',ROOT/'scripts/native_regional_skin_engine_probe.cpp',ROOT/'scripts/native_body_ports.h',ROOT/'scripts/native_tissue_ports.h']
    for p in files:(out/p.name).write_bytes(p.read_bytes())
    generated=out/'native_regional_engine_step.inc';generated.write_text(step_source(ENGINE.read_text()))
    preparation={'variant':VARIANT.name,'variant_manifest_sha256':MANIFEST_PIN,'library_sha256':LIBRARY_PIN,
        'native_engine_source_sha256':ENGINE_PIN,'generated_step_sha256':sha(generated),'native_state_path':str(state),'native_state_sha256':sha(state),
        'state_reference_manifest_sha256':sha(source_receipt),'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},
        'native_modes':['baseline','zero','load'],'steps_per_mode':12,'native_step_s':.02,
        'scope':'Actual native LoadState, not fresh stabilization; full-engine lifecycle with aggregate source chemical laws retained'}
    (out/'preparation.json').write_text(json.dumps(preparation,indent=2)+'\n');print(out,flush=True)
    if not args.run:return
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';binary=out/'probe'
    command=['c++','-std=c++20','-O0','-g1','-DBIOGEARS_THROW_READONLY_EXCEPTIONS',str(out/'native_regional_skin_engine_probe.cpp')]
    for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:
        command+=['-I',str(p)]
    command+=[str(VARIANT/'libbiogears.so.8.0.0'),'-L',str(lib),'-lbiogears_cdm','-o',str(binary)]
    env={**os.environ,'LD_LIBRARY_PATH':str(VARIANT)+':'+str(lib),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    started=time.monotonic()
    def run(cmd,label,cwd,memory,cpu,wall):
        with (out/(label+'.stdout')).open('w') as stdout,(out/(label+'.stderr')).open('w') as stderr:
            result=subprocess.run(['/usr/bin/time','-v','-o',str(out/(label+'.resources'))]+cmd,cwd=cwd,env=env,stdout=stdout,stderr=stderr,preexec_fn=cap(memory,cpu),timeout=wall)
        if result.returncode:raise RuntimeError(f'{label} failed ({result.returncode}); retained {out}')
    try:
        run(command,'compile',out,1024**3,60,90)
        linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(out/'ldd.txt').write_text(linkage)
        assert str(VARIANT/'libbiogears.so.8.0.0') in linkage
        sweat_command=list(command)
        sweat_command[sweat_command.index(str(out/'native_regional_skin_engine_probe.cpp'))]=str(out/'native_regional_sweat_fixture.cpp')
        sweat_command[sweat_command.index('-o')+1]=str(out/'sweat_fixture')
        run(sweat_command,'sweat_compile',out,1024**3,60,90)
        sweat_work=out/'sweat_fixture_work';sweat_work.mkdir()
        for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:
            (sweat_work/name).symlink_to(build/'runtime'/name)
        run([str(out/'sweat_fixture')],'sweat_fixture',sweat_work,1024**3,30,40)
        sweat_result=json.loads(next(line[7:] for line in (out/'sweat_fixture.stdout').read_text().splitlines() if line.startswith('RESULT ')))
        assert sweat_result['passed']
        frames={}
        for mode in preparation['native_modes']:
            work=out/mode;work.mkdir()
            for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:
                (work/name).symlink_to(build/'runtime'/name)
            run([str(binary),mode,str(state)],mode,work,2*1024**3,30,40)
            frames[mode]=[json.loads(line[len('RESULT '):]) for line in (out/(mode+'.stdout')).read_text().splitlines() if line.startswith('RESULT ')]
            assert len(frames[mode])==13, 'Incomplete native frames: '+mode
        checks={};max_relative=max_absolute=0.;unavailable=[];parity_failures=[]
        for base,zero in zip(frames['baseline'],frames['zero']):
            checks['native_clock_'+str(base['tick'])]=abs(base['time_s']-zero['time_s'])<1e-12
            for key,a in base['values'].items():
                if key.startswith('acceptance.'):continue
                b=zero['values'].get(key)
                if a is None or b is None:
                    if a is not None and b is None:unavailable.append(key)
                    continue
                error=abs(a-b);relative=error/max(abs(a),abs(b),1.)
                max_absolute=max(max_absolute,error);max_relative=max(max_relative,relative)
                if error>1e-8+1e-8*max(abs(a),abs(b)):
                    parity_failures.append({'tick':base['tick'],'key':key,'baseline':a,'regional':b,'absolute_error':error})
        checks['zero_load_common_native_ports']=not parity_failures
        max_fluid=max_species=max_sweat=0.
        for mode,records in frames.items():
            for record in records:
                for key,value in record['values'].items():
                    if key.endswith(('fluid_step_residual_ml','volume_ownership_residual_ml','install_volume_residual_ml')):
                        assert value is not None
                        max_fluid=max(max_fluid,abs(value))
                    if key.endswith(('ownership_residual_ug','install_mass_residual_ug')):
                        assert value is not None
                        max_species=max(max_species,abs(value))
                    if key.endswith('paired_mass_residual_mg'):
                        assert value is not None
                        max_sweat=max(max_sweat,abs(value))
        checks['native_fluid_incidence_and_ownership']=max_fluid<1e-8
        checks['species_ownership']=max_species<1e-6
        checks['native_sweat_paired_mass']=max_sweat<1e-8
        loaded=frames['load'][6]['values'];unloaded=frames['load'][8]['values']
        checks['local_load_applied']=abs(loaded['acceptance.region.region_a.external_pressure_pa']-133.322387415)<1e-9
        checks['residual_unloaded']=loaded['acceptance.region.residual.external_pressure_pa']==0
        checks['unload_applied']=unloaded['acceptance.region.region_a.external_pressure_pa']==0
        zero=frames['zero'][6]['values']
        dp=loaded['acceptance.region.region_a.pressure_mmhg']-zero['acceptance.region.region_a.pressure_mmhg']
        dq=loaded['acceptance.region.region_a.lymph_flow_ml_per_s']-zero['acceptance.region.region_a.lymph_flow_ml_per_s']
        checks['local_native_pressure_and_flow_response']=abs(dp)>1e-6 and abs(dq)>1e-10
        report={**preparation,'passed':all(checks.values()),'checks':checks,'compile_command':command,'binary_sha256':sha(binary),'sweat_fixture':sweat_result,
            'wall_s':time.monotonic()-started,'maximum_common_port_absolute_error':max_absolute,'maximum_common_port_scaled_error':max_relative,
            'parity_failures':parity_failures,'unavailable_original_active_circuit_ports':sorted(set(unavailable)),
            'maximum_native_fluid_residual_ml':max_fluid,'maximum_species_ownership_residual_ug':max_species,
            'maximum_native_sweat_paired_residual_mg':max_sweat,
            'local_pressure_change_mmhg':dp,'local_lymph_flow_change_ml_per_s':dq,
            'regional_chemical_independence_claimed':False,'whole_body_species_closure_claimed':False,
            'resources':{name:(out/(name+'.resources')).read_text() for name in ['compile','sweat_compile','sweat_fixture','baseline','zero','load']}}
        (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps({k:v for k,v in report.items() if k not in ['resources','compile_command','parity_failures']},indent=2))
        if not report['passed']:raise RuntimeError('Native regional engine acceptance failed; see retained verification.json')
    except BaseException as error:
        (out/'failure.json').write_text(json.dumps({'error':str(error),'wall_s':time.monotonic()-started},indent=2)+'\n');raise


if __name__=='__main__':main()
