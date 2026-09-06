#!/usr/bin/env python3
"""Two bounded builds with isolated corrected methods; actual schema objects."""
import os,json,subprocess,tempfile,resource,argparse,shutil,shlex
from pathlib import Path
from prepare_gi_serialization_repair import prepare,BASE
from build_biogears_shared_donor_variant import ROOT,RUNTIME,SOURCE,sha
from verify_gi_shared_donor import limits

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--version',choices=['original','corrected'],required=True);parser.add_argument('--objects-manifest',type=Path);args=parser.parse_args()
    out=Path(tempfile.mkdtemp(prefix='native-gi-codec-',dir=ROOT/'data/derived/audits'));print(out,flush=True);repair=prepare(out)
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';variant=RUNTIME/'variants/whole_body_integrity_gi_absorption'
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    probe=ROOT/'scripts/native_gi_serialization_probe.cpp';cmd=['c++','-std=c++20','-O0',str(probe)]
    includes=build/'projects/biogears/libBiogears/CMakeFiles/libbiogears.dir/includes_CXX.rsp';cmd+=['@'+str(includes)]
    inventory=json.loads((variant/'manifest.json').read_text())['object_sha256'];object_cwd=build/'projects/biogears/libBiogears'
    objects=[];replaced={}
    replacements={}
    if args.version=='corrected':
        if not args.objects_manifest:raise ValueError('Compiled full-TU manifest required')
        replacements=json.loads(args.objects_manifest.read_text())['objects']
        assert len(replacements)==3
    for name,digest in inventory.items():
        path=(object_cwd/name).resolve();assert sha(path)==digest
        if name in replacements:
            rec=replacements[name];path=Path(rec['object']);assert sha(path)==rec['sha256'];assert sha(Path(rec['patched_source']))==rec['patched_source_sha256'];replaced[name]=rec
        objects.append(str(path))
    assert len(replaced)==len(replacements)
    rsp=out/'objects.rsp';rsp.write_text(' '.join(shlex.quote(p) for p in objects))
    cmd+=['@'+str(rsp),str(lib/'libbiogears_common_st.a'),str(lib/'libbiogears_cdm.so.8.0.0'),'-ldl',str(RUNTIME/'sysroot/usr/lib/aarch64-linux-gnu/libxerces-c-3.2.so')]
    def bounded():
        limits();resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    def run(args,name):
        with (out/(name+'.stdout')).open('w') as stdout,(out/(name+'.stderr')).open('w') as stderr:
            p=subprocess.run(['/usr/bin/time','-v','-o',str(out/(name+'.resources'))]+args,cwd=out,env=env,stdout=stdout,stderr=stderr,timeout=90,preexec_fn=bounded)
        return p.returncode
    commands={}
    for version in [args.version]:
        command=cmd+['-o',str(out/version)];commands[version]=command
        if run(command,version+'-compile'):raise RuntimeError(f'compile failed {out}')
    multi=subprocess.check_output(['c++','-print-multiarch'],text=True).strip();env['LD_LIBRARY_PATH']=f'{variant}:{lib}:{RUNTIME}/sysroot/usr/lib/{multi}'
    for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:(out/name).symlink_to(build/'runtime'/name)
    codes={}
    for version in [args.version]:
        for mode in ['roundtrip','nested','cleanup']:codes[version+'/'+mode]=run([str(out/version),mode],version+'-'+mode)
    passed=all((codes[args.version+'/'+m]!=0 if args.version=='original' else codes[args.version+'/'+m]==0) for m in ['roundtrip','nested','cleanup'])
    report=dict(binary_sha256=sha(out/args.version),passed=passed,meaning="Original failure reproduction" if args.version=="original" else "Corrected roundtrip verification",version=args.version,returncodes=codes,scope='Actual original/corrected native schema codecs; FixtureGI only overrides SetUp to avoid patient preparation. Schema-object roundtrip, not whole-engine XML reload.',commands=commands,native_object_sha256=inventory,replaced_full_translation_units=replaced,source_hashes=json.loads((out/'source_hashes.json').read_text()),probe_sha256=sha(probe),repair_sha256=sha(repair),include_response_sha256=sha(includes),export_header_sha256=sha(BASE.parent/'libCDM/include/biogears/cdm-exports.h'),library_sha256=sha(variant/'libbiogears.so.8.0.0'),resources={p.name:p.read_text() for p in out.glob('*.resources')})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'passed':passed,'returncodes':codes}))
    if not passed:raise RuntimeError('regression failed')
if __name__=='__main__':main()
