"""Source-compiled constitutive probe, with immutable execution inputs."""
from pathlib import Path
import argparse,json,os,subprocess,sys,tempfile
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from scripts.build_biogears_thermal_boundary_variant import RUNTIME,sha
from ihm.assembly.systemic_evidence import freeze_sources
from ihm.assembly.native_environment_evidence import freeze_native_environment

def run(variant):
    if variant not in ('whole_body_integrity_skin_perfusion','whole_body_integrity_sweat_evaporation','whole_body_integrity_evaporation_humidity'):raise ValueError('Explicit audit variant required')
    library=RUNTIME/'variants'/variant/'libbiogears.so.8.0.0';vm=json.loads((library.parent/'manifest.json').read_text())
    if sha(library)!=vm['library_sha256']:raise ValueError('Variant library changed')
    out=Path(tempfile.mkdtemp(prefix='native-evaporation-'+variant.removeprefix('whole_body_integrity_')+'-',dir=BASE/'data/derived/audits'))
    donor=BASE/'data/raw/physiology/biogears';build=RUNTIME/'biogears-build';sr=RUNTIME/'sysroot';lib=build/'outputs/Release/lib'
    original=BASE/'scripts/native_evaporation_probe.cpp';source=out/original.name;source.write_bytes(original.read_bytes());exe=out/'probe'
    command=['c++','-std=c++20','-O2',str(source)]
    for p in (donor/'projects/biogears/libBiogears/include',donor/'projects/biogears-common/include',sr/'usr/include',sr/'usr/include/eigen3',build/'projects/biogears/generated/Release'):command+=['-I',str(p)]
    command+=['-L',str(lib),f'-Wl,-rpath,{lib}','-lbiogears','-lbiogears_cdm','-o',str(exe)]
    (out/'compile.json').write_text(json.dumps(command,indent=2)+'\n');subprocess.run(command,check=True)
    for name in ('patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml'):
        target=build/'runtime'/name;(out/name).symlink_to(target,target_is_directory=target.is_dir())
    state=BASE/'data/derived/audits/thermal-perfusion-fresh-hour/states/native_stabilized.xml'
    sources={str(p.relative_to(BASE)):sha(p) for p in (Path(__file__),original,source,exe,library,library.parent/'manifest.json',state)}
    freeze_sources(BASE,out,sources)
    env=os.environ.copy();env['LD_LIBRARY_PATH']=str(library.parent)+os.pathsep+env.get('LD_LIBRARY_PATH','')
    linkage=subprocess.check_output(['ldd',str(exe)],env=env,text=True);deps={}
    for line in linkage.splitlines():
        if 'not found' in line:raise ValueError('Unresolved dependency')
        if '=>' in line:
            p=Path(line.split('=>',1)[1].split(' (',1)[0].strip())
            if p.is_file():deps[str(p.resolve())]=sha(p)
    if deps.get(str(library.resolve()))!=sha(library):raise ValueError('Incorrect native variant linkage')
    args=[str(exe),str(state),'evaporation.csv']
    (out/'manifest.json').write_text(json.dumps(dict(schema='ihm.native-session.v1',adapter_kind='audit_constitutive_probe',command=args,dependency_sha256=deps,executable_sha256=sha(exe),library_sha256=sha(library),source_hashes=sources),indent=2)+'\n')
    (out/'receipts.jsonl').write_text(json.dumps({'command':'CONSTITUTIVE_PROBE','argv':args})+'\n')
    with (out/'stdout.log').open('w') as log:r=subprocess.run(args,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT)
    (out/'execution.json').write_text(json.dumps(dict(returncode=r.returncode,command=args,source_hashes=sources),indent=2)+'\n')
    if r.returncode:raise RuntimeError('Native probe failed; outputs retained at '+str(out))
    freeze_native_environment(BASE,out)
    print(json.dumps({'output_dir':str(out),'status':'executed','variant':variant},indent=2));return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('variant');run(p.parse_args().variant)
