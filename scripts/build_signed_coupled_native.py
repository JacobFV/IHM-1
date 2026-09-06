"""Build the signed thin adapter without replacing existing native executables."""
from pathlib import Path
import argparse,hashlib,json,os,resource,shutil,subprocess,tempfile,time
from ihm.native import SOURCE,RUNTIME,BASE,SOURCE_REVISION
from ihm.native.coupled_session import SignedCoupledNativeSession

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--variant',default='whole_body_integrity_signed_muscle_v2');args=parser.parse_args()
    variant=RUNTIME/'variants'/args.variant;manifest_path=variant/'manifest.json';manifest_raw=manifest_path.read_bytes();variant_manifest=json.loads(manifest_raw)
    if sha(variant/'libbiogears.so.8.0.0')!=variant_manifest['library_sha256']:raise ValueError('Signed library changed')
    if sha(BASE/'scripts/native_signed_muscle_port.h')!=variant_manifest['header_sha256']:raise ValueError('Signed ABI header changed')
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
        'source_revision':SOURCE_REVISION,'command':command,'variant':args.variant,'variant_manifest_sha256':hashlib.sha256(manifest_raw).hexdigest(),
        'library_sha256':variant_manifest['library_sha256'],'wall_s':time.monotonic()-started,
        'peak_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'retained_build':str(out.relative_to(BASE))}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    target=RUNTIME/SignedCoupledNativeSession.executable_name;temporary=target.with_suffix('.building');shutil.copy2(binary,temporary);temporary.replace(target)
    target.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'passed':True,'retained_build':str(out),'wall_s':manifest['wall_s'],'peak_child_rss_kib':manifest['peak_child_rss_kib']}))

if __name__=='__main__':main()
