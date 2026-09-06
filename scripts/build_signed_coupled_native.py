"""Build the signed thin adapter without replacing existing native executables."""
from pathlib import Path
import argparse,hashlib,json,os,resource,shlex,shutil,subprocess,tempfile,time
BASE=Path(__file__).resolve().parents[1]
SOURCE=BASE/'data/raw/physiology/biogears';RUNTIME=BASE/'data/runtime/physiology'
SOURCE_REVISION='3f16a5fa1dade9c511b88d923606fa51cc35e95d'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def intake_source_receipt(variant):
    """Bind the actual compiled consumer and its native payload/action ownership."""
    manifest=json.loads((variant/'manifest.json').read_text())
    if sha(variant/'libbiogears.so.8.0.0')!=manifest['library_sha256']:raise ValueError('Intake runtime library changed')
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    objects=shlex.split((variant/'objects.rsp').read_text())
    if len(objects)!=len(set(objects)) or set(objects)!=set(manifest['object_sha256']):raise ValueError('Intake variant object inventory changed')
    records={}
    for name,relative in [('Gastrointestinal','engine/Systems/Gastrointestinal.cpp'),
                          ('SEPatientActionCollection','cdm/scenario/SEPatientActionCollection.cpp'),
                          ('SENutrition','cdm/patient/SENutrition.cpp')]:
        matches=[o for o in objects if o.endswith('/'+name+'.cpp.o')]
        if len(matches)!=1:raise ValueError('Ambiguous native intake source object')
        obj=matches[0];resolved=cwd/obj
        if sha(resolved)!=manifest['object_sha256'][obj]:raise ValueError('Compiled intake owner changed')
        original=SOURCE/'projects/biogears/libBiogears/src'/relative
        committed=subprocess.check_output(['git','-C',str(SOURCE),'show',f'{SOURCE_REVISION}:{original.relative_to(SOURCE)}'])
        if original.read_bytes()!=committed:raise ValueError('Held native intake owner source changed')
        source=resolved.with_suffix('') if 'variants' in resolved.parts else original
        if source!=original and sha(source) not in (source.parent/'manifest.json').read_text():raise ValueError('Inherited intake source not bound by manifest')
        records[name]={'object_path':str(resolved),'object_sha256':sha(resolved),'source_path':str(source.relative_to(BASE)),'source_sha256':sha(source)}
        if name=='Gastrointestinal':
            text=source.read_text();start=text.index('void Gastrointestinal::PreProcess()');end=text.index('    DigestNutrient();',start)
            block=text[start:end]
            if hashlib.sha256(block.encode()).hexdigest()!='64ce7582c64ff38e9ab0161b2f7d7f156a9f3a605edd06dc73ed21a1855d6e5b':raise ValueError('Audited Active GI intake consumer changed')
            if text.count('RemoveConsumeNutrients()')!=1:raise ValueError('Unexpected native GI action-removal path')
            records[name]['active_consumer_sha256']=hashlib.sha256(block.encode()).hexdigest()
    # The adapter does not cancel or reload actions during an advance. A new
    # engine-side remover would invalidate inference from disappearance.
    engine=SOURCE/'projects/biogears/libBiogears/src/engine'
    removers=[str(p.relative_to(engine)) for p in engine.rglob('*.cpp') if b'RemoveConsumeNutrients(' in p.read_bytes()]
    if removers!=['Systems/Gastrointestinal.cpp']:raise ValueError('Intake action no longer has the audited sole engine consumer')
    return {'variant':variant.name,'variant_manifest_sha256':sha(variant/'manifest.json'),'library_sha256':manifest['library_sha256'],'owners':records,
            'scope':'Audited Active GI increment/remove boundary; no cancellation, native state reload or nutrition-file command during adapter advance'}

def main():
    from ihm.native.coupled_session import SignedCoupledNativeSession
    parser=argparse.ArgumentParser();parser.add_argument('--variant',default='whole_body_integrity_signed_muscle_v2');args=parser.parse_args()
    variant=RUNTIME/'variants'/args.variant;manifest_path=variant/'manifest.json';manifest_raw=manifest_path.read_bytes();variant_manifest=json.loads(manifest_raw)
    if sha(variant/'libbiogears.so.8.0.0')!=variant_manifest['library_sha256']:raise ValueError('Signed library changed')
    if sha(BASE/'scripts/native_signed_muscle_port.h')!=variant_manifest['header_sha256']:raise ValueError('Signed ABI header changed')
    intake_guards={name:intake_source_receipt(RUNTIME/'variants'/name) for name in ('whole_body_integrity_signed_muscle_v2','whole_body_integrity_gi_absorption')}
    sources=[BASE/'scripts'/name for name in SignedCoupledNativeSession.adapter_sources]
    sources.append(SOURCE/'projects/biogears/libBiogears/src/engine/Controller/BioGearsEngine.cpp')
    frozen={p:p.read_bytes() for p in sources}
    out=Path(tempfile.mkdtemp(prefix='signed-adapter-',dir=RUNTIME));build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';sysroot=RUNTIME/'sysroot'
    for p,raw in frozen.items():(out/p.name).write_bytes(raw)
    command=['prlimit','--as=4294967296','--','nice','-n','10','c++','-std=c++20','-O1',str(out/'native_biogears_signed.cpp')]
    for p in (SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',sysroot/'usr/include',sysroot/'usr/include/eigen3',build/'projects/biogears/generated/Release'):
        command+=['-I',str(p)]
    binary=out/'native_biogears_signed'
    command+=['-L',str(variant),'-L',str(lib),f'-Wl,-rpath,{variant}:{lib}','-l:libbiogears.so.8.0.0','-lbiogears_cdm','-o',str(binary)]
    started=time.monotonic()
    with (out/'compile.log').open('w') as log:
        subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),timeout=120)
    if any(p.read_bytes()!=raw for p,raw in frozen.items()) or manifest_path.read_bytes()!=manifest_raw:raise ValueError('Signed build inputs changed')
    manifest={'executable_sha256':sha(binary),'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(raw).hexdigest() for p,raw in frozen.items()},
        'intake_consumer_guards':intake_guards,'source_revision':SOURCE_REVISION,'command':command,'variant':args.variant,'variant_manifest_sha256':hashlib.sha256(manifest_raw).hexdigest(),
        'library_sha256':variant_manifest['library_sha256'],'wall_s':time.monotonic()-started,
        'peak_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'retained_build':str(out.relative_to(BASE))}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    target=RUNTIME/SignedCoupledNativeSession.executable_name;temporary=target.with_suffix('.building');shutil.copy2(binary,temporary);temporary.replace(target)
    target.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'passed':True,'retained_build':str(out),'wall_s':manifest['wall_s'],'peak_child_rss_kib':manifest['peak_child_rss_kib']}))

if __name__=='__main__':main()
