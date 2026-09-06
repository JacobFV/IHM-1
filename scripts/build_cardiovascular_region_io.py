"""One isolated Circuit IO TU; preserve all inherited variant objects and libraries."""
from pathlib import Path
import argparse,json,os,resource,shlex,shutil,subprocess,time
from patch_cardiovascular_region_io import ROOT,IO,IO_SHA,sha
RUNTIME=ROOT/'data/runtime/physiology'
def build(overlay,parent_id,variant_id,reuse=None):
 parent=RUNTIME/'variants'/parent_id;pm=json.loads((parent/'manifest.json').read_text());om=json.loads((overlay/'manifest.json').read_text());cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears'
 if 'cdm_object_sha256' in pm:raise ValueError('Sleep-schema branch composition requires its explicit paired CDM build recipe; not silently supported')
 if sha(parent/'libbiogears.so.8.0.0')!=pm['library_sha256'] or sha(IO)!=IO_SHA:raise ValueError('Changed parent/IO source')
 for name,digest in om['files'].items():
  if sha(overlay/name)!=digest:raise ValueError('Changed IO overlay')
 objects=shlex.split((parent/'objects.rsp').read_text())
 if len(objects)!=len(set(objects)) or set(objects)!=set(pm['object_sha256']) or any(sha(cwd/o)!=pm['object_sha256'][o] for o in objects):raise ValueError('Changed inherited object inventory')
 slots=[o for o in objects if o.endswith('/io/cdm/Circuit.cpp.o')]
 if len(slots)!=1:raise ValueError('One original Circuit IO object required')
 out=RUNTIME/'variants'/variant_id
 if out.exists():raise ValueError('Fresh variant name required')
 out.mkdir();frozen={parent/'manifest.json':sha(parent/'manifest.json'),overlay/'manifest.json':sha(overlay/'manifest.json'),Path(__file__).resolve():sha(__file__),parent/'libbiogears.so.8.0.0':pm['library_sha256'],cwd/'../../../outputs/Release/lib/libbiogears_cdm.so.8.0.0':sha(cwd/'../../../outputs/Release/lib/libbiogears_cdm.so.8.0.0')}
 frozen.update({cwd/o:pm['object_sha256'][o] for o in objects})
 shutil.copyfile(overlay/'Circuit.cpp',out/'Circuit.cpp');frozen[out/'Circuit.cpp']=sha(out/'Circuit.cpp')
 for name in ('Circuit.h','Property.h'):
  original=IO.parent/name
  if sha(original)!=sha(overlay/name):raise ValueError('Changed pinned IO header')
  frozen[original]=sha(original)
 cardio=[o for o in objects if o.endswith('/Cardiovascular.cpp.o')];assert len(cardio)==1;cardio_source=(cwd/cardio[0]).with_suffix('')
 shutil.copyfile(cardio_source,out/'Cardiovascular.cpp');frozen[cardio_source]=sha(cardio_source)
 for name in ('native_signed_muscle_port.h','native_signed_vascular_prior.h'):
  source=cardio_source.parent/name
  if source.exists():shutil.copyfile(source,out/name);frozen[out/name]=sha(out/name)
 source=out/'Circuit.cpp';obj=out/'Circuit.cpp.o';command=None;started=time.monotonic();env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
 if reuse:
  prior=json.loads((reuse/'manifest.json').read_text())
  if prior['io_source_sha256']!=sha(source) or sha(reuse/'Circuit.cpp.o')!=prior['io_object_sha256'] or prior['paired_cdm_library_sha256']!=sha(cwd/'../../../outputs/Release/lib/libbiogears_cdm.so.8.0.0'):raise ValueError('Changed/mismatched reusable IO object')
  shutil.copyfile(reuse/'Circuit.cpp.o',obj);frozen[reuse/'Circuit.cpp.o']=sha(obj)
 else:
  command=['prlimit','--as=4294967296','--','nice','-n','10','c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-I',str(IO.parent),'-std=gnu++20','-fPIC','-O1','-DNDEBUG','-c',str(source),'-o',str(obj)]
  with (out/'compile.log').open('w') as log:subprocess.run(command,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=30)
 objects[objects.index(slots[0])]=str(obj);response=out/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects));link=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text());lib=out/'libbiogears.so.8.0.0';link[link.index('-o')+1]=str(lib);link=['@'+str(response) if x=='@CMakeFiles/libbiogears.dir/objects1.rsp' else x for x in link];link=['prlimit','--as=4294967296','--','nice','-n','10',*link]
 with (out/'link.log').open('w') as log:subprocess.run(link,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=30)
 if any(sha(p)!=digest for p,digest in frozen.items()):raise ValueError('Build inputs changed; no publication')
 record={'schema':'ihm.cardiovascular-region-io-variant.v1','variant':variant_id,'parent_variant':parent_id,'parent_manifest_sha256':sha(parent/'manifest.json'),'parent_library_sha256':pm['library_sha256'],'library_sha256':sha(lib),'paired_cdm_library_sha256':sha(cwd/'../../../outputs/Release/lib/libbiogears_cdm.so.8.0.0'),'header_sha256':pm['header_sha256'],'source_sha256':sha(out/'Cardiovascular.cpp'),'io_source_sha256':sha(source),'io_object_sha256':sha(obj),'reuse_io_from':None if reuse is None else str(reuse),'object_sha256':{o:sha(cwd/o) for o in objects},'compile_command':command,'link_command':link,'wall_s':time.monotonic()-started,'maximum_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'scope':'Isolated existing-field serialization correction; no automatic legacy metadata migration or default promotion'}
 (out/'manifest.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:v for k,v in record.items() if k not in ('object_sha256','link_command')},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--overlay',required=True,type=Path);p.add_argument('--parent',default='whole_body_integrity_gi_absorption');p.add_argument('--variant',required=True);p.add_argument('--reuse-io-from',type=Path);p.add_argument('--build',action='store_true');a=p.parse_args()
 if not a.build:raise SystemExit('Explicit coordinated --build required')
 build(a.overlay.resolve(),a.parent,a.variant,None if a.reuse_io_from is None else a.reuse_io_from.resolve())
