"""Explicit one-TU isolated experiment build; no adapters/default pointers changed."""
from pathlib import Path
import argparse,hashlib,json,os,resource,shlex,shutil,subprocess,time
from patch_signed_cardiovascular_reader import ROOT,SOURCE,SOURCE_SHA,prepare
RUNTIME=ROOT/'data/runtime/physiology'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def build(overlay,parent_id,kind,variant_id):
    if kind not in ('reader_only','vascular_prior'):raise ValueError('Explicit experiment kind required')
    parent=RUNTIME/'variants'/parent_id;pm=json.loads((parent/'manifest.json').read_text());cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
    manifest=json.loads((overlay/'manifest.json').read_text())
    if sha(parent/'libbiogears.so.8.0.0')!=pm['library_sha256'] or sha(SOURCE)!=SOURCE_SHA:raise ValueError('Pinned parent changed')
    for name,digest in manifest['files'][kind].items():
        if sha(overlay/kind/name)!=digest:raise ValueError('Frozen overlay changed')
    objects=shlex.split((parent/'objects.rsp').read_text())
    if len(objects)!=len(set(objects)) or set(objects)!=set(pm['object_sha256']) or any(sha(cwd/o)!=pm['object_sha256'][o] for o in objects):raise ValueError('Parent object inventory changed')
    cardio=[o for o in objects if o.endswith('/Cardiovascular.cpp.o')]
    if len(cardio)!=1 or (cwd/cardio[0]).with_suffix('').resolve()!=SOURCE.resolve():raise ValueError('Unexpected inherited cardiovascular source')
    variant=RUNTIME/'variants'/variant_id
    if variant.exists():raise ValueError('Fresh experiment name required')
    variant.mkdir();frozen={parent/'manifest.json':sha(parent/'manifest.json'),SOURCE:SOURCE_SHA,Path(__file__):sha(__file__),overlay/'manifest.json':sha(overlay/'manifest.json')}
    for name in manifest['files'][kind]:shutil.copyfile(overlay/kind/name,variant/name);frozen[overlay/kind/name]=sha(overlay/kind/name)
    source=variant/'Cardiovascular.cpp';obj=variant/'Cardiovascular.cpp.o'
    cmd=['prlimit','--as=4294967296','--','nice','-n','10','c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O1','-DNDEBUG','-c',str(source),'-o',str(obj)]
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1');started=time.monotonic()
    with (variant/'compile.log').open('w') as log:subprocess.run(cmd,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
    objects[objects.index(cardio[0])]=str(obj);response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    link=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text());lib=variant/'libbiogears.so.8.0.0';link[link.index('-o')+1]=str(lib);link=['@'+str(response) if x=='@CMakeFiles/libbiogears.dir/objects1.rsp' else x for x in link]
    link=['prlimit','--as=4294967296','--','nice','-n','10',*link]
    with (variant/'link.log').open('w') as log:subprocess.run(link,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
    if any(sha(p)!=digest for p,digest in frozen.items()):raise ValueError('Input changed; unpublished experiment retained')
    record={'schema':'ihm.signed-cardiovascular-experiment.v1','variant':variant_id,'experiment':kind,'parent_variant':parent_id,'parent_manifest_sha256':sha(parent/'manifest.json'),'parent_library_sha256':pm['library_sha256'],'library_sha256':sha(lib),'header_sha256':sha(variant/'native_signed_muscle_port.h'),'source_sha256':sha(source),'object_sha256':{o:sha(cwd/o) for o in objects},'compile_command':cmd,'link_command':link,'wall_s':time.monotonic()-started,'maximum_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'scope':'Isolated experimental library only; no adapter/default selection or calibrated vascular claim'}
    (variant/'manifest.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k!='object_sha256'},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--overlay',type=Path);p.add_argument('--parent',default='whole_body_integrity_gi_absorption');p.add_argument('--kind',choices=('reader_only','vascular_prior'));p.add_argument('--variant');p.add_argument('--build',action='store_true');a=p.parse_args()
    if not a.build:print(prepare(a.overlay))
    else:
        if not a.overlay or not a.kind or not a.variant:p.error('Explicit overlay/kind/fresh variant required; coordinated heavy slot must be granted')
        build(a.overlay.resolve(),a.parent,a.kind,a.variant)
