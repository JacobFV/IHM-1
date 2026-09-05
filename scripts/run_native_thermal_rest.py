"""Retained, uninterrupted long-rest diagnostic; shared runner is not modified."""
from pathlib import Path
import argparse,difflib,json,math,os,subprocess,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts.build_biogears_thermal_boundary_variant import BASE,RUNTIME,sha
from ihm.assembly.systemic_evidence import freeze_sources
from ihm.assembly.native_environment_evidence import freeze_native_environment

def run(output,variant,seconds=21600,state=None):
    if variant not in ('whole_body_integrity_thermal_boundary_v2','whole_body_integrity_skin_perfusion','whole_body_integrity_sweat_evaporation','whole_body_integrity_evaporation_humidity'):raise ValueError('Explicit thermal audit variant required')
    if isinstance(seconds,bool) or not isinstance(seconds,(int,float)) or not math.isfinite(seconds) or not 0<seconds<=21600 or abs(seconds*50-round(seconds*50))>1e-8:raise ValueError('Bounded native clock required')
    out=Path(output).resolve()
    if not out.is_relative_to(BASE) or out.exists():raise ValueError('Fresh in-workspace diagnostic directory required')
    if state is not None and not Path(state).resolve().is_file():raise ValueError('State unavailable')
    library=RUNTIME/'variants'/variant/'libbiogears.so.8.0.0';vm=json.loads((library.parent/'manifest.json').read_text())
    if sha(library)!=vm['library_sha256']:raise ValueError('Variant library changed')
    original=BASE/'scripts/native_biogears_rest.cpp';before=original.read_text();old='seconds <= 3600'
    if before.count(old)!=1:raise ValueError('Bounded source runner changed')
    after=before.replace(old,'seconds <= 21600')
    out.mkdir();source=out/'long_rest_runner.cpp';source.write_text(after)
    (out/'runner_horizon.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='native_biogears_rest.cpp',tofile='audit/long_rest_runner.cpp')))
    exe=out/'long_rest_runner';donor=BASE/'data/raw/physiology/biogears';build=RUNTIME/'biogears-build';sysroot=RUNTIME/'sysroot';lib=build/'outputs/Release/lib'
    command=['c++','-std=c++20','-O2',str(source)]
    for p in (donor/'projects/biogears/libBiogears/include',donor/'projects/biogears-common/include',sysroot/'usr/include',sysroot/'usr/include/eigen3',build/'projects/biogears/generated/Release'):command+=['-I',str(p)]
    command+=['-L',str(lib),f'-Wl,-rpath,{lib}','-lbiogears','-lbiogears_cdm','-o',str(exe)];subprocess.run(command,check=True)
    for name in ('patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml'):
        target=build/'runtime'/name;(out/name).symlink_to(target,target_is_directory=target.is_dir())
    (out/'states').mkdir();(out/'timeline.tsv').write_text('')
    sources={str(p.relative_to(BASE)):sha(p) for p in (Path(__file__),original,source,exe,library,library.parent/'manifest.json')}
    if state:sources[str(Path(state).resolve().relative_to(BASE))]=sha(state)
    freeze_sources(BASE,out,sources)
    env=os.environ.copy();env['LD_LIBRARY_PATH']=str(library.parent)+os.pathsep+env.get('LD_LIBRARY_PATH','')
    linkage=subprocess.check_output(['ldd',str(exe)],env=env,text=True);dependencies={}
    for line in linkage.splitlines():
        if 'not found' in line:raise ValueError('Unresolved native library')
        if '=>' in line:
            path=Path(line.split('=>',1)[1].split(' (',1)[0].strip())
            if path.is_file():dependencies[str(path.resolve())]=sha(path)
    if dependencies.get(str(library.resolve()))!=sha(library):raise ValueError('Incorrect native library resolution')
    args=[str(exe),str(seconds),'IHMGenericMale',str(Path(state).resolve()) if state else '-','1','timeline.tsv','-','-','-']
    config={'variant':variant,'seconds':seconds,'state_path':str(state) if state else None,'ambient_temperature_c':22,'clothing_clo':.5,
        'kind':'batch thermal audit, not NativeSession protocol or accepted systemic experiment','source_hashes':sources}
    (out/'configuration.json').write_text(json.dumps(config,indent=2)+'\n')
    # Generic execution-input schema for the independent environment archiver.
    # Adapter kind explicitly distinguishes this batch executable from the stream.
    (out/'manifest.json').write_text(json.dumps({'schema':'ihm.native-session.v1','adapter_kind':'audit_batch_rest','command':args,
        'executable_sha256':sha(exe),'dependency_sha256':dependencies,'library_sha256':sha(library),'configuration':config},indent=2)+'\n')
    (out/'receipts.jsonl').write_text(json.dumps({'command':'EXECUTE_BATCH_REST','before_ticks':0,'argv':args})+'\n')
    start=time.monotonic()
    with (out/'runner_stdout.log').open('w') as log:r=subprocess.run(args,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT)
    (out/'execution.json').write_text(json.dumps({'returncode':r.returncode,'wall_seconds':time.monotonic()-start,'command':args,'source_hashes':sources},indent=2)+'\n')
    if r.returncode:raise RuntimeError('Thermal rest diagnostic failed; outputs retained')
    freeze_native_environment(BASE,out)
    print(json.dumps({'output':str(out),'returncode':r.returncode,'wall_seconds':time.monotonic()-start},indent=2),flush=True)
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--variant',required=True);p.add_argument('--seconds',type=float,default=21600);p.add_argument('--state',type=Path);a=p.parse_args();run(a.output,a.variant,a.seconds,a.state)
