"""Enable the held Moco source backend in an isolated object/library variant."""
from pathlib import Path
import concurrent.futures,hashlib,json,os,re,shlex,subprocess,time,urllib.request,zipfile
ROOT=Path(__file__).resolve().parents[1];RUNTIME=ROOT/'data/runtime/moco3d';SOURCE=ROOT/'data/raw/mechanics/opensim-core';OPENSIM=ROOT/'data/runtime/opensim'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if os.uname().machine!='aarch64':raise RuntimeError('This pinned native artifact is for aarch64; acquire a matching primary CasADi distribution for another host')
 RUNTIME.mkdir(parents=True,exist_ok=True);raw=ROOT/'data/raw/mechanics/casadi-3.6.5';raw.mkdir(exist_ok=True)
 name='casadi-3.6.5-cp311-none-manylinux2014_aarch64.whl';archive=raw/name
 if not (raw/'wheel_receipt.json').exists():
  metadata=json.load(urllib.request.urlopen('https://pypi.org/pypi/casadi/3.6.5/json'));entry=next(f for f in metadata['urls'] if f['filename']==name)
  urllib.request.urlretrieve(entry['url'],archive)
  (raw/'wheel_receipt.json').write_text(json.dumps({'artifact':entry,'source':'https://pypi.org/pypi/casadi/3.6.5/json'},indent=2)+'\n')
 receipt=json.loads((raw/'wheel_receipt.json').read_text());assert sha(archive)==receipt['artifact']['digests']['sha256']
 casadi=RUNTIME/'casadi-wheel/casadi'
 if not casadi.exists():
  with zipfile.ZipFile(archive) as z:z.extractall(RUNTIME/'casadi-wheel')
 with zipfile.ZipFile(archive) as z:
  for entry in z.infolist():
   if not entry.is_dir() and (RUNTIME/'casadi-wheel'/entry.filename).read_bytes()!=z.read(entry):raise RuntimeError('Extracted CasADi artifact changed: '+entry.filename)
 revision=subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip();assert revision=='86b30588374650fbaf012a345a836a64f6855522'
 subprocess.run(['git','-C',str(SOURCE),'diff','--exit-code','HEAD','--','OpenSim/Moco'],check=True,stdout=subprocess.DEVNULL)
 cwd=OPENSIM/'opensim-build/OpenSim/Moco';flags=(cwd/'CMakeFiles/osimMoco.dir/flags.make').read_text()
 args=[]
 for key in ['CXX_DEFINES','CXX_INCLUDES','CXX_FLAGS']:args+=shlex.split(re.search(r'^'+key+r' = (.*)$',flags,re.M).group(1))
 args+=['-DOPENSIM_WITH_CASADI','-I',str(casadi/'include')]
 backend=SOURCE/'OpenSim/Moco/MocoCasADiSolver'
 files=[backend/'MocoCasADiSolver.cpp']+sorted(backend.glob('CasOC*.cpp'))+[backend/'MocoCasOCProblem.cpp']
 objects=RUNTIME/'objects';objects.mkdir(exist_ok=True);logs=RUNTIME/('build-'+str(time.time_ns()));logs.mkdir()
 commands=[]
 def compile_file(p):
  command=['c++',*args,'-c',str(p),'-o',str(objects/(p.name+'.o'))];commands.append(command)
  with (logs/(p.name+'.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(compile_file,files))
 link=shlex.split((cwd/'CMakeFiles/osimMoco.dir/link.txt').read_text());old='CMakeFiles/osimMoco.dir/MocoCasADiSolver/MocoCasADiSolver.cpp.o';assert link.count(old)==1
 link[link.index(old)]=str(objects/'MocoCasADiSolver.cpp.o');link[link.index('-o')+1]=str(RUNTIME/'libosimMoco.so')
 link += [str(objects/(p.name+'.o')) for p in files if p.name!='MocoCasADiSolver.cpp']+[str(casadi/'libcasadi.so'),'-Wl,-rpath,'+str(casadi)]
 with (logs/'link.log').open('w') as log:subprocess.run(link,cwd=cwd,stdout=log,stderr=subprocess.STDOUT,check=True)
 all_objects=[Path(x) if Path(x).is_absolute() else cwd/x for x in link if x.endswith('.o')]
 record={'architecture':os.uname().machine,'source_revision':revision,'source_files_sha256':{str(p.relative_to(SOURCE)):sha(p) for p in (SOURCE/'OpenSim/Moco').rglob('*') if p.is_file()},'wheel_receipt_sha256':sha(raw/'wheel_receipt.json'),'wheel_sha256':sha(archive),'compile_commands':commands,'link_command':link,'build_cwd':str(cwd),'object_sha256':{str(p):sha(p) for p in all_objects},'library_sha256':sha(RUNTIME/'libosimMoco.so'),'original_library_sha256':sha(OPENSIM/'install/opensim/lib/libosimMoco.so'),'logs':str(logs),'scope':'Backend-enabled isolated Moco library. Existing OpenSim and SCONE binaries/libraries unchanged.'}
 (RUNTIME/'build_manifest.json').write_text(json.dumps(record,indent=2)+'\n');print(RUNTIME/'libosimMoco.so')
if __name__=='__main__':main()
