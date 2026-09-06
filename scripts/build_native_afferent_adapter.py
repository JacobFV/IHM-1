"""Prepare isolated observation-only adapters; --run requires the heavy queue."""
from pathlib import Path
import argparse,hashlib,json,os,resource,shutil,subprocess,tempfile,time
ROOT=Path(__file__).resolve().parents[1];RUNTIME=ROOT/'data/runtime/physiology'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def prepare(regional=False,run=False):
    base='native_biogears_regional_signed' if regional else 'native_biogears_signed'
    name='native_biogears_afferent_regional' if regional else 'native_biogears_afferent_signed'
    base_path=RUNTIME/(base+'.manifest.json');raw=base_path.read_bytes();record=json.loads(raw)
    if sha(RUNTIME/base)!=record['executable_sha256']:raise ValueError('Parent adapter executable changed')
    frozen={ROOT/p:(ROOT/p).read_bytes() for p in record['source_sha256']}
    if any(hashlib.sha256(v).hexdigest()!=record['source_sha256'][str(p.relative_to(ROOT))] for p,v in frozen.items()):raise ValueError('Parent adapter source changed')
    source=ROOT/'data/raw/physiology/biogears/projects/biogears/libBiogears'
    for path in (ROOT/'scripts/native_nervous_afferents.h',ROOT/'scripts'/f'{name}.cpp',source/'include/biogears/engine/Systems/Nervous.h',source/'src/engine/Systems/Nervous.cpp'):
        frozen[path]=path.read_bytes()
    out=Path(tempfile.mkdtemp(prefix='native-afferent-adapter-',dir=RUNTIME))
    retained=ROOT/record['retained_build']
    for path in retained.iterdir():
        if path.suffix in ('.h','.cpp','.inc'):shutil.copy2(path,out/path.name)
    for path,content in frozen.items():(out/path.name).write_bytes(content)
    command=list(record['command']);old_binary=command[command.index('-o')+1];command[command.index('-o')+1]=str(out/name)
    matches=[i for i,item in enumerate(command) if item.endswith('/'+base+'.cpp')]
    if len(matches)!=1:raise ValueError('Parent adapter compile source is ambiguous')
    command[matches[0]]=str(out/(name+'.cpp'))
    library=RUNTIME/'variants'/record['variant']/'libbiogears.so.8.0.0'
    if sha(library)!=record['library_sha256']:raise ValueError('Native observer library differs from parent adapter')
    result={**record,'parent_adapter_manifest_sha256':hashlib.sha256(raw).hexdigest(),
        'parent_adapter_executable_sha256':record['executable_sha256'],'observer_scope':'Source-internal readonly native Nervous observer; no layout cast, state writer or receptor equation',
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(v).hexdigest() for p,v in frozen.items()},
        'observer_header_sha256':sha(ROOT/'scripts/native_nervous_afferents.h'),'nervous_header_sha256':sha(source/'include/biogears/engine/Systems/Nervous.h'),
        'builder_sha256':sha(Path(__file__)),'retained_build':str(out.relative_to(ROOT)),'command':command}
    result.pop('executable_sha256',None)
    (out/'preparation.json').write_text(json.dumps(result,indent=2)+'\n')
    if not run:return {'prepared':True,'output':str(out)}
    started=time.monotonic()
    with (out/'compile.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120,
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1'))
    if base_path.read_bytes()!=raw or sha(library)!=record['library_sha256'] or any(p.read_bytes()!=v for p,v in frozen.items()):raise ValueError('Native observer build inputs changed')
    result.update(executable_sha256=sha(out/name),wall_s=time.monotonic()-started,peak_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    (out/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    target=RUNTIME/name;temporary=target.with_suffix('.building');shutil.copy2(out/name,temporary);temporary.replace(target)
    target.with_suffix('.manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    return {'passed':True,'output':str(out),'wall_s':result['wall_s']}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--regional',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args();print(json.dumps(prepare(a.regional,a.run)))
