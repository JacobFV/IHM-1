"""Bounded isolated existing-schema metadata roundtrip; no patient advance."""
from pathlib import Path
import argparse,json,os,resource,shutil,subprocess,tempfile,time
from patch_cardiovascular_region_io import ROOT,IO,sha

def main():
 p=argparse.ArgumentParser();p.add_argument('--variant',required=True);p.add_argument('--run',action='store_true');a=p.parse_args()
 if not a.run:raise SystemExit('Explicit coordinated native slot required')
 started=time.monotonic();runtime=ROOT/'data/runtime/physiology';variant=runtime/'variants'/a.variant;manifest=json.loads((variant/'manifest.json').read_text());library=variant/'libbiogears.so.8.0.0'
 assert sha(library)==manifest['library_sha256']
 out=Path(tempfile.mkdtemp(prefix='cardiovascular-region-io-native-',dir=ROOT/'data/derived/audits'));source=ROOT/'scripts/native_cardiovascular_region_io_roundtrip.cpp';shutil.copyfile(source,out/source.name)
 donor=ROOT/'data/raw/physiology/biogears';build=runtime/'biogears-build';libs=build/'outputs/Release/lib';sysroot=runtime/'sysroot';binary=out/'probe'
 frozen={source:sha(source),library:sha(library),variant/'manifest.json':sha(variant/'manifest.json'),IO.parent/'Circuit.h':sha(IO.parent/'Circuit.h'),Path(__file__).resolve():sha(__file__)}
 command=['prlimit','--as=4294967296','--','nice','-n','10','c++','-std=c++20','-O1',str(out/source.name)]
 for path in (IO.parent,IO.parent.parent.parent,donor/'projects/biogears/libCDM/include',donor/'projects/biogears/libBiogears/include',donor/'projects/biogears-common/include',sysroot/'usr/include',sysroot/'usr/include/eigen3',build/'projects/biogears/generated/Release'):command+=['-I',str(path)]
 object_root=build/'projects/biogears/libBiogears';property_slots=[name for name in manifest['object_sha256'] if name.endswith('/io/cdm/Property.cpp.o')];assert len(property_slots)==1
 objects=[variant/'Circuit.cpp.o',object_root/property_slots[0]]
 assert sha(objects[0])==manifest['io_object_sha256'] and sha(objects[1])==manifest['object_sha256'][property_slots[0]]
 frozen.update({obj:sha(obj) for obj in objects})
 command += [str(obj) for obj in objects]
 command+=['-L',str(variant),'-L',str(libs),f'-Wl,-rpath,{variant}:{libs}','-l:libbiogears.so.8.0.0','-lbiogears_cdm','-L',str(sysroot/'usr/lib/aarch64-linux-gnu'),'-lxerces-c','-o',str(binary)]
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',LD_LIBRARY_PATH=f'{variant}:{libs}:{sysroot}/usr/lib/aarch64-linux-gnu')
 (out/'command.json').write_text(json.dumps(command,indent=2))
 with (out/'compile.log').open('w') as log:subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=30)
 linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);assert str(library) in linkage and 'not found' not in linkage;(out/'ldd.txt').write_text(linkage)
 run=subprocess.run(['prlimit','--as=4294967296','--','nice','-n','10',str(binary)],env=env,capture_output=True,text=True,timeout=5);(out/'stdout.log').write_text(run.stdout);(out/'stderr.log').write_text(run.stderr);run.check_returncode();result=json.loads(run.stdout)
 assert result['passed'] and all(sha(path)==digest for path,digest in frozen.items())
 result.update(wall_s=time.monotonic()-started,maximum_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,source_sha256={str(path):digest for path,digest in frozen.items()},executable_sha256=sha(binary),scope='Exact corrected Circuit object and inherited Property object linked into isolated probe for private IO calls; existing Circuit DTO roundtrip, no patient advance or history recovery')
 (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'output':str(out),**result},indent=2))
if __name__=='__main__':main()
