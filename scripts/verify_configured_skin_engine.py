#!/usr/bin/env python3
"""Isolated named-region native full-engine acceptance; default source preparation."""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import tempfile
from prepare_configured_regional_skin import ROOT,PARENT,PARENT_MANIFEST_SHA,sha
from verify_native_regional_skin_engine import ENGINE,ENGINE_PIN,step_source
from ihm.assembly.regional_skin_configuration import validate_configuration

RUNTIME=ROOT/'data/runtime/physiology';SOURCE=ROOT/'data/raw/physiology/biogears'
CONFIGURED=ROOT/'data/research/configured_regional_skin/prepared_v1'
ACCEPTED=ROOT/'data/derived/audits/native-regional-engine-cu9dn9or'
STAGE=ROOT/'data/research/configured_skin_engine/prepared_v3'


def prepare():
    config=json.loads((CONFIGURED/'configuration.json').read_text());validate_configuration(config)
    configured=json.loads((CONFIGURED/'preparation.json').read_text())
    old=json.loads((ACCEPTED/'preparation.json').read_text())
    if sha(PARENT/'manifest.json')!=PARENT_MANIFEST_SHA or sha(ENGINE)!=ENGINE_PIN:
        raise ValueError('Changed accepted native source/parent')
    state=Path(old['native_state_path'])
    if sha(state)!=old['native_state_sha256']:raise ValueError('Changed held native initial state')
    inputs={'native_configured_regional_skin.h':CONFIGURED/'native_configured_regional_skin.h',
            'native_configured_skin_engine_probe.cpp':ROOT/'scripts/native_configured_skin_engine_probe.cpp'}
    if sha(inputs['native_configured_regional_skin.h'])!=configured['configured_header_sha256']:raise ValueError('Changed configured header')
    for name in ['native_regional_species.h','native_body_ports.h','native_tissue_ports.h']:
        path=ACCEPTED/name
        if sha(path)!=old['source_sha256']['scripts/'+name]:raise ValueError('Changed frozen accepted port/helper')
        inputs[name]=path
    blobs={name:path.read_bytes() for name,path in inputs.items()}
    generated=step_source(ENGINE.read_text())
    # Exactly the held engine method plus two regional hooks: no extra chemistry.
    if __import__('hashlib').sha256(generated.encode()).hexdigest()!=old['generated_step_sha256']:
        raise ValueError('Native engine lifecycle changed')
    blobs['native_regional_engine_step.inc']=generated.encode()
    receipt={'schema':'configured_skin_engine_preparation_v1','configuration_sha256':config['configuration_sha256'],
             'region_names':[r['name'] for r in config['regions']],
             'fractions':[r['initial_fraction'] for r in config['regions']],
             'parent_variant':PARENT.name,'parent_manifest_sha256':PARENT_MANIFEST_SHA,
             'parent_library_sha256':configured['parent_library_sha256'],
             'native_state_path':str(state),'native_state_sha256':sha(state),
             'native_engine_source_sha256':ENGINE_PIN,'accepted_fullengine_preparation_sha256':sha(ACCEPTED/'preparation.json'),
             'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs.values()},
             'prepared_sha256':{name:__import__('hashlib').sha256(raw).hexdigest() for name,raw in blobs.items()},
             'allowed_missing_original_ports':json.loads((ACCEPTED/'verification.json').read_text())['unavailable_original_active_circuit_ports'],
             'steps_per_mode':12,'dt_s':.02,'modes':['baseline','zero','load'],
             'native_executed':False,'native_activation_allowed':False,
             'scope':'Held native LoadState plus full circuit/solute lifecycle; no restart equivalence or fresh stabilization claim'}
    if STAGE.exists():
        for name,raw in blobs.items():
            if (STAGE/name).read_bytes()!=raw:raise ValueError('Frozen configured engine source changed')
        if json.loads((STAGE/'preparation.json').read_text())!=receipt:raise ValueError('Frozen preparation changed')
    else:
        STAGE.mkdir(parents=True)
        for name,raw in blobs.items():(STAGE/name).write_bytes(raw)
        (STAGE/'preparation.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def assess(frames,receipt):
    names=receipt['region_names'];fractions=receipt['fractions'];failures=[];unavailable=set()
    maxima={'common_port_scaled_error':0.,'mass_scaled_parity':0.,'fluid_residual_ml':0.,'species_ownership_ug':0.,'sweat_residual_mg':0.,'law_scaled_residual':0.}
    def check(condition,label):
        if not condition:failures.append(label)
    def near(a,b,label,group,tol=1e-8):
        if a is None or b is None or not math.isfinite(a) or not math.isfinite(b):check(False,label);return
        scaled=abs(a-b)/max(1.,abs(a),abs(b));maxima[group]=max(maxima[group],scaled);check(scaled<tol,label)
    for mode,records in frames.items():
        check(len(records)==13,'frame_count.'+mode)
        for tick,record in enumerate(records):
            check(record['tick']==tick and record['configuration_sha256']==receipt['configuration_sha256'],'identity.'+mode+'.'+str(tick))
            near(record['elapsed_s'],tick*.02,'clock.'+mode+'.'+str(tick),'common_port_scaled_error',1e-12)
            v=record['values']
            check(v.get('acceptance.liquid_owners.external_boundary_count')==1,'declared_Ambient_boundary_count')
            check(v.get('acceptance.boundary.Ambient.volume_is_positive_infinity')==1,'declared_Ambient_volume')
            for key,value in v.items():
                if key.endswith(('fluid_step_residual_ml','volume_ownership_residual_ml','install_volume_residual_ml')):
                    check(value is not None and abs(value)<1e-8,key);maxima['fluid_residual_ml']=max(maxima['fluid_residual_ml'],abs(value or 0))
                if key.endswith(('ownership_residual_ug','install_mass_residual_ug')):
                    check(value is not None and abs(value)<1e-6,key);maxima['species_ownership_ug']=max(maxima['species_ownership_ug'],abs(value or 0))
                if key.endswith('paired_mass_residual_mg'):
                    check(value is not None and abs(value)<1e-8,key);maxima['sweat_residual_mg']=max(maxima['sweat_residual_mg'],abs(value or 0))
            if mode!='baseline':
                for i,name in enumerate(names):near(v.get('acceptance.region.'+name+'.initial_fraction'),fractions[i],'fraction.'+name,'law_scaled_residual',1e-15)
                for key,value in v.items():
                    if not key.startswith('acceptance.original_path.'):continue
                    suffix=key[len('acceptance.original_path.'):]
                    regional=[v.get('acceptance.region.'+name+'.path.'+suffix) for name in names]
                    if suffix.endswith('resistance_mmhg_s_per_ml'):
                        for i,r in enumerate(regional):near(None if r is None else r*fractions[i],value,'resistance.'+suffix,'law_scaled_residual')
                    else:
                        near(None if None in regional else math.fsum(regional),value,'extensive_law.'+suffix,'law_scaled_residual')
    for b,z in zip(frames['baseline'],frames['zero']):
        bv,zv=b['values'],z['values']
        check(zv['acceptance.liquid_owners.count']==bv['acceptance.liquid_owners.count']+2,'owning_leaf_count')
        for key,a in bv.items():
            if key.startswith('acceptance.boundary.'):
                near(a,zv.get(key),key,'common_port_scaled_error');continue
            if key.startswith('acceptance.'):
                if not (key.startswith(('acceptance.skin.species.','acceptance.liquid_owners.species.')) and key.endswith('.mass_ug')):continue
                near(a,zv.get(key),key,'mass_scaled_parity');continue
            value=zv.get(key)
            if a is None or value is None:
                if a is not None and value is None:unavailable.add(key)
                continue
            near(a,value,key,'common_port_scaled_error')
        near(bv['acceptance.liquid_owners.volume_ml'],zv['acceptance.liquid_owners.volume_ml'],'all_liquid_volume','common_port_scaled_error')
    check(unavailable.issubset(set(receipt['allowed_missing_original_ports'])),'unexpected_missing_native_ports')
    loaded=frames['load'][6]['values'];zero=frames['zero'][6]['values'];released=frames['load'][8]['values']
    key='acceptance.region.'+names[0]+'.'
    near(loaded[key+'external_pressure_pa'],133.322387415,'load_applied','law_scaled_residual')
    check(released[key+'external_pressure_pa']==0,'release_applied')
    check(all(loaded['acceptance.region.'+name+'.external_pressure_pa']==0 for name in names[1:]),'other_regions_unloaded')
    dp=loaded[key+'pressure_mmhg']-zero[key+'pressure_mmhg'];dq=loaded[key+'lymph_flow_ml_per_s']-zero[key+'lymph_flow_ml_per_s']
    check(abs(dp)>1e-6 and abs(dq)>1e-10,'local_pressure_and_drainage_response')
    species=sorted(k[len('acceptance.skin.species.'):-len('.mass_ug')] for k in loaded if k.startswith('acceptance.skin.species.') and k.endswith('.mass_ug'))
    check(len(species)==28 and 'Albumin' in species,'complete_skin_species')
    for mode in ['zero','load']:
        for record in frames[mode]:
            for name in names:
                for sub in species:
                    value=record['values'].get('acceptance.region.'+name+'.'+sub+'.mass_ug')
                    check(value is not None and math.isfinite(value) and value>=0,'complete_region_mass.'+name+'.'+sub)
    return {'passed':not failures,'failures':failures,'maxima':maxima,'species_names':species,
            'unavailable_original_active_circuit_ports':sorted(unavailable),'local_pressure_change_mmhg':dp,
            'local_lymph_flow_change_ml_per_s':dq,'whole_body_species_closure_claimed':False,
            'regional_chemical_independence_claimed':False,'exact_restart_equivalence_claimed':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',action='store_true');args=parser.parse_args();receipt=prepare()
    if not args.run:print(json.dumps({'prepared':True,'native_executed':False,'configuration_sha256':receipt['configuration_sha256']}));return
    if sha(PARENT/'libbiogears.so.8.0.0')!=receipt['parent_library_sha256']:raise ValueError('Changed parent library')
    out=Path(tempfile.mkdtemp(prefix='configured-skin-engine-',dir=ROOT/'data/derived/audits'))
    for p in STAGE.iterdir():shutil.copyfile(p,out/p.name)
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';binary=out/'probe'
    cmd=['c++','-std=c++20','-O0','-g1','-DBIOGEARS_THROW_READONLY_EXCEPTIONS',str(out/'native_configured_skin_engine_probe.cpp')]
    for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:
        cmd+=['-I',str(p)]
    cmd+=[str(PARENT/'libbiogears.so.8.0.0'),'-L',str(lib),'-lbiogears_cdm','-o',str(binary)]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':str(PARENT)+':'+str(lib)}
    def run(command,label,cwd,memory,wall):
        def limits():resource.setrlimit(resource.RLIMIT_AS,(memory,memory));resource.setrlimit(resource.RLIMIT_CPU,(60,60))
        with (out/(label+'.stdout')).open('w') as stdout,(out/(label+'.stderr')).open('w') as stderr:
            result=subprocess.run(['/usr/bin/time','-v','-o',str(out/(label+'.resources'))]+command,cwd=cwd,env=env,stdout=stdout,stderr=stderr,preexec_fn=limits,timeout=wall)
        if result.returncode:raise RuntimeError(label+' failed; retained '+str(out))
    print(str(out),flush=True)
    try:
        run(cmd,'compile',out,1024**3,75)
        linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(out/'ldd.txt').write_text(linkage)
        if str(PARENT/'libbiogears.so.8.0.0') not in linkage:raise ValueError('Wrong native lineage')
        frames={}
        for mode in receipt['modes']:
            work=out/mode;work.mkdir()
            for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:(work/name).symlink_to(build/'runtime'/name)
            run([str(binary),mode,receipt['native_state_path']],mode,work,2*1024**3,10)
            frames[mode]=[json.loads(line[7:]) for line in (out/(mode+'.stdout')).read_text().splitlines() if line.startswith('RESULT ')]
        report=assess(frames,receipt);report.update(preparation=receipt,native_executed=True,compile_command=cmd,binary_sha256=sha(binary),production_activation_allowed=False)
        (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');prepare()
        print(json.dumps({k:v for k,v in report.items() if k not in ['preparation','compile_command']}))
        if not report['passed']:raise RuntimeError('Configured full-engine acceptance failed')
    except BaseException as error:
        (out/'failure.json').write_text(json.dumps({'error':str(error)},indent=2)+'\n');raise


if __name__=='__main__':main()
