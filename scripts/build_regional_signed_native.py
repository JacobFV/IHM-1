#!/usr/bin/env python3
"""Prepare an immutable regional signed adapter; compile only with explicit --run."""
from pathlib import Path
import argparse,hashlib,json,os,resource,shlex,shutil,subprocess,tempfile,time
BASE=Path(__file__).resolve().parents[1]
SOURCE=BASE/'data/raw/physiology/biogears';RUNTIME=BASE/'data/runtime/physiology'
VARIANT=RUNTIME/'variants/whole_body_integrity_regional_skin_graph_v2'
MANIFEST_PIN='52de403e3dab2682f774755b7bf1370711bc5ab218244c581460b84a82aad130'
LIBRARY_PIN='bc91cbab829c04bfa2df7490433af5bee752e7514332ac715e9f0bd11a17650a'
REVISION='3f16a5fa1dade9c511b88d923606fa51cc35e95d'
ENGINE=SOURCE/'projects/biogears/libBiogears/src/engine/Controller/BioGearsEngine.cpp'
ENGINE_PIN='5c29a09c536a876d569069625bbe4ef2c9e375e3ba12c48c64b94ceed49cbc95'
ADAPTER_SOURCES=('native_biogears_regional_signed.cpp','native_regional_coupled_engine.h','native_body_ports.h',
                 'native_tissue_ports.h','native_regional_skin.h','native_regional_species.h','native_signed_muscle_port.h')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def step_source(text):
    signature='bool BioGearsEngine::AdvanceModelTime(bool appendDataTrack)'
    if text.count(signature)!=1:raise ValueError('Native step signature changed')
    start=text.index(signature);end=text.index('\n//-------------------------------------------------------------------------------',start)
    method=text[start:end].replace(signature,'bool AdvanceModelTime(bool appendDataTrack=false) override')
    edits=[('  PreProcess();','  begin_coupled_step();\n  PreProcess();\n  coupled_after_preprocess();'),
           ('  PostProcess();','  PostProcess();\n  coupled_after_postprocess();')]
    for old,new in edits:
        if method.count(old)!=1:raise ValueError('Native lifecycle anchor changed')
        method=method.replace(old,new)
    # All native guards, Process, events, clocks and tracking must remain byte-identical.
    restored=method
    for old,new in reversed(edits):restored=restored.replace(new,old)
    restored=restored.replace('bool AdvanceModelTime(bool appendDataTrack=false) override',signature)
    if restored!=text[start:end]:raise ValueError('Generated native lifecycle changed beyond hooks')
    return method+'\n'

def verify_wire(binary,out):
    directory=out/'protocol';directory.mkdir()
    runtime=RUNTIME/'biogears-build/runtime';lib=RUNTIME/'biogears-build/outputs/Release/lib'
    for name in ('patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml'):
        (directory/name).symlink_to(runtime/name,target_is_directory=(runtime/name).is_dir())
    (directory/'states').mkdir()
    source_receipt=BASE/'data/derived/systemic/exertion_v3/exercise/native/manifest.json'
    state=Path(json.loads(source_receipt.read_text())['configuration']['state_path']);frozen_state=directory/'input-state.xml'
    shutil.copyfile(state,frozen_state)
    ref=hashlib.sha256(b'regional protocol fixture, no physiologic reference calibration').hexdigest()
    operations=[f'signed_step {ref} 0 0 0','regional_skin_pressure region_a 133.322387415',
        f'signed_step {ref} 0 0 0','regional_skin_pressure region_a 0',f'signed_step {ref} 0 0 0',
        'skin_compression 0','save','step 1','regional_skin_pressure unknown 0','regional_skin_pressure region_b -1',
        'respiratory_load 50',f'signed_step {ref} .01 .01 0','respiratory_load 0',f'signed_step {ref} 0 0 0','quit']
    commands=''.join(f'{i} {op}\n' for i,op in enumerate(operations,1));(directory/'commands.txt').write_text(commands)
    launch=['prlimit','--as=2147483648','--','nice','-n','10',str(binary),'StandardMale',str(frozen_state),'5']
    env={**os.environ,'LD_LIBRARY_PATH':str(VARIANT)+':'+str(lib),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    resolved=subprocess.run(['ldd',str(binary)],env=env,text=True,capture_output=True,check=True).stdout
    (directory/'resolved-libraries.txt').write_text(resolved)
    if str(VARIANT/'libbiogears.so.8.0.0') not in resolved:raise ValueError('Wrong regional library resolved')
    started=time.monotonic();proc=subprocess.run(launch,cwd=directory,env=env,input=commands,text=True,capture_output=True,timeout=45)
    (directory/'stdout.log').write_text(proc.stdout);(directory/'stderr.log').write_text(proc.stderr);proc.check_returncode()
    frames=[json.loads(line.split('\t',1)[1]) for line in proc.stdout.splitlines() if line.startswith('IHM\t')]
    assert len(frames)==len(operations)+1
    by_seq={frame['sequence']:frame for frame in frames};assert by_seq[0]['status']=='ready'
    assert [frame['sequence'] for frame in frames]==list(range(16))
    for seq in (6,7,8,9,10):
        assert by_seq[seq]['status']=='rejected'
        assert by_seq[seq]['elapsed_s']==by_seq[5]['elapsed_s'] and by_seq[seq]['values']==by_seq[5]['values']
    assert by_seq[15]['status']=='closed'
    step_seq=(1,3,5,12,14)
    for tick,seq in enumerate(step_seq,1):
        frame=by_seq[seq];values=frame['values'];assert frame['status']=='ok' and abs(frame['elapsed_s']-.02*tick)<1e-12
        assert values['coupling.muscle_heat_count']==values['coupling.muscle_tissue_count']==1
        assert values['tissue.regional_skin.completed_native_steps']==tick
        assert abs(values['tissue.regional_skin.aggregate.volume_ownership_residual_ml'])<1e-8
        assert abs(values['tissue.regional_skin.aggregate.fluid_step_residual_ml'])<1e-8
        for key,value in values.items():
            if key.startswith('tissue.regional_skin.') and key.endswith(('ownership_residual_ug','paired_mass_residual_mg')):assert abs(value)<1e-7,key
        for region,fraction in [('region_a',.2),('region_b',.3),('residual',.5)]:
            prefix='tissue.regional_skin.'+region+'.'
            assert values[prefix+'fraction']==fraction
            assert abs(values[prefix+'fluid_step_residual_ml'])<1e-8
    assert by_seq[2]['elapsed_s']==by_seq[1]['elapsed_s']
    assert by_seq[2]['values']['tissue.regional_skin.region_a.requested_pressure_pa']==133.322387415
    assert by_seq[2]['values']['tissue.regional_skin.region_a.applied_pressure_pa']==0
    assert by_seq[3]['values']['tissue.regional_skin.region_a.applied_pressure_pa']==133.322387415
    assert by_seq[5]['values']['tissue.regional_skin.region_a.applied_pressure_pa']==0
    assert by_seq[3]['values']['tissue.regional_skin.region_a.pressure_mmhg']-by_seq[1]['values']['tissue.regional_skin.region_a.pressure_mmhg']>.5
    assert by_seq[12]['values']['coupling.external_pressure_pa']==50
    assert by_seq[14]['values']['coupling.external_pressure_pa']==0
    report={'passed':True,'native_steps':5,'wire_commands':len(operations),'state_sha256':sha(frozen_state),
        'state_source_manifest_sha256':sha(source_receipt),'library_sha256':LIBRARY_PIN,'launch':launch,
        'wall_s':time.monotonic()-started,'frames':frames,
        'scope':'Regional topology lifecycle + signed/respiratory wire integration; engineering fractions, no anatomic calibration or long-run validation'}
    (directory/'verification.json').write_text(json.dumps(report,indent=2)+'\n');return report

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',action='store_true');parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if args.verify and not args.run:raise ValueError('--verify requires explicit --run and coordinated native slot')
    if sha(VARIANT/'manifest.json')!=MANIFEST_PIN or sha(VARIANT/'libbiogears.so.8.0.0')!=LIBRARY_PIN:raise ValueError('Regional variant changed')
    if sha(ENGINE)!=ENGINE_PIN:raise ValueError('Held native lifecycle source changed')
    manifest=json.loads((VARIANT/'manifest.json').read_text())
    for name,expected in manifest['header_receipts'].items():
        if sha(BASE/'scripts'/name)!=expected or sha(VARIANT/name)!=expected:raise ValueError('Regional/signed ABI header changed: '+name)
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears';objects=shlex.split((VARIANT/'objects.rsp').read_text())
    if len(objects)!=len(set(objects)) or set(objects)!=set(manifest['object_sha256']):raise ValueError('Regional object inventory mismatch')
    for obj in objects:
        if sha(cwd/obj)!=manifest['object_sha256'][obj]:raise ValueError('Regional object changed: '+obj)
    sources=[BASE/'scripts'/name for name in ADAPTER_SOURCES]+[ENGINE]
    frozen={p:p.read_bytes() for p in sources};out=Path(tempfile.mkdtemp(prefix='regional-signed-adapter-',dir=RUNTIME))
    for path,raw in frozen.items():(out/path.name).write_bytes(raw)
    generated=out/'native_regional_coupled_step.inc';generated.write_text(step_source(ENGINE.read_text()))
    record={'source_revision':REVISION,'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(raw).hexdigest() for p,raw in frozen.items()},
        'variant':VARIANT.name,'variant_manifest_sha256':MANIFEST_PIN,'library_sha256':LIBRARY_PIN,'generated_step_sha256':sha(generated),
        'builder_sha256':sha(Path(__file__)),'engineering_fractions':{'region_a':.2,'region_b':.3,'residual':.5},'retained_build':str(out.relative_to(BASE))}
    (out/'preparation.json').write_text(json.dumps(record,indent=2)+'\n')
    if not args.run:print(json.dumps({'prepared':True,'output':str(out)}));return
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';binary=out/'native_biogears_regional_signed'
    command=['prlimit','--as=4294967296','--','nice','-n','10','c++','-std=c++20','-O1','-DBIOGEARS_THROW_READONLY_EXCEPTIONS',str(out/'native_biogears_regional_signed.cpp')]
    for path in (SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release'):
        command+=['-I',str(path)]
    command += [str(VARIANT/'libbiogears.so.8.0.0'),'-L',str(lib),f'-Wl,-rpath,{VARIANT}:{lib}','-lbiogears_cdm','-o',str(binary)]
    started=time.monotonic()
    with (out/'compile.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
    if any(p.read_bytes()!=raw for p,raw in frozen.items()) or sha(VARIANT/'manifest.json')!=MANIFEST_PIN:raise ValueError('Build inputs changed')
    record.update(executable_sha256=sha(binary),command=command,wall_s=time.monotonic()-started,peak_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    if args.verify:record['verification']= {k:v for k,v in verify_wire(binary,out).items() if k!='frames'}
    (out/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    target=RUNTIME/binary.name;temporary=target.with_suffix('.building');shutil.copy2(binary,temporary);temporary.replace(target)
    target.with_suffix('.manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'passed':True,'verified':args.verify,'output':str(out),'wall_s':record['wall_s'],'peak_child_rss_kib':record['peak_child_rss_kib']}))
if __name__=='__main__':main()
