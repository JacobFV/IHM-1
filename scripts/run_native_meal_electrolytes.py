"""Bounded same-state nutrition counterion probes; no model-equation changes."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json,os,subprocess,sys,tempfile,time
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from scripts.audit_long_horizon_thermal import sha
from ihm.assembly.systemic_evidence import freeze_sources
from ihm.assembly.native_environment_evidence import freeze_native_environment,materialize_native_resources

def run():
    runtime=BASE/'data/runtime/physiology';build=runtime/'biogears-build';lib=build/'outputs/Release/lib';donor=BASE/'data/raw/physiology/biogears';sr=runtime/'sysroot'
    variant='whole_body_integrity_evaporation_humidity';library=runtime/'variants'/variant/'libbiogears.so.8.0.0'
    assert sha(library)=='ab485812968b17e16915b04a6e5144bb1e1ec7f5bf78213936c3334d087c0a5a'
    state=BASE/'data/derived/audits/thermal-final-initial-state-ctkr8bap/inputs/data/derived/audits/thermal-evaporation-humidity-fresh-hour/states/native_stabilized.xml'
    out=Path(tempfile.mkdtemp(prefix='meal-electrolyte-native-',dir=BASE/'data/derived/audits'));print(out,flush=True)
    for name in ('native_meal_electrolyte_probe.cpp','native_body_ports.h'):(out/name).write_bytes((BASE/'scripts'/name).read_bytes())
    exe=out/'probe';command=['c++','-std=c++20','-O2',str(out/'native_meal_electrolyte_probe.cpp')]
    for p in (donor/'projects/biogears/libBiogears/include',donor/'projects/biogears-common/include',sr/'usr/include',sr/'usr/include/eigen3',build/'projects/biogears/generated/Release'):command+=['-I',str(p)]
    command+=['-L',str(lib),f'-Wl,-rpath,{lib}','-lbiogears','-lbiogears_cdm','-o',str(exe)];subprocess.run(command,check=True)
    (out/'compile.json').write_text(json.dumps(command,indent=2)+'\n')
    env={**os.environ,'LD_LIBRARY_PATH':str(library.parent)+':'+str(lib),'OPENBLAS_NUM_THREADS':'1'}
    linkage=subprocess.check_output(['ldd',str(exe)],env=env,text=True);assert 'not found' not in linkage
    deps={str(Path(s.split('=>')[1].split(' (')[0].strip()).resolve()):sha(Path(s.split('=>')[1].split(' (')[0].strip())) for s in linkage.splitlines() if '=>' in s}
    assert deps[str(library.resolve())]==sha(library)
    def case(mode):
        work=out/mode;work.mkdir();initial=work/'initial_state.xml';initial.write_bytes(state.read_bytes());assert sha(initial)==sha(state)
        names=('patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml')
        for name in names:
            target=build/'runtime'/name;(work/name).symlink_to(target,target_is_directory=target.is_dir())
        args=[str(exe),str(initial),mode,'600']
        sources={str(p.relative_to(BASE)):sha(p) for p in (Path(__file__),out/'native_meal_electrolyte_probe.cpp',out/'native_body_ports.h',exe,initial,state,library,library.parent/'manifest.json')}
        freeze_sources(BASE,work,sources)
        (work/'manifest.json').write_text(json.dumps(dict(schema='ihm.native-session.v1',adapter_kind='nutrition_counterion_audit',command=args,dependency_sha256=deps,executable_sha256=sha(exe),library_sha256=sha(library),configuration=dict(mode=mode,seconds=600,state_sha256=sha(initial),variant=variant),source_hashes=sources),indent=2)+'\n')
        freeze_native_environment(BASE,work,before_start=True);resources=materialize_native_resources(BASE,work)
        for name in names:
            (work/name).unlink();(work/name).symlink_to(resources/name,target_is_directory=(resources/name).is_dir())
        (work/'receipts.jsonl').write_text(json.dumps({'command':'RUN_NUTRITION_COUNTERION_AUDIT','argv':args})+'\n')
        start=time.monotonic()
        with (work/'stdout.log').open('w') as log:result=subprocess.run(args,cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=900)
        receipt=dict(mode=mode,returncode=result.returncode,wall_seconds=time.monotonic()-start,source_hashes=sources)
        (work/'execution.json').write_text(json.dumps(receipt,indent=2)+'\n')
        if result.returncode:raise RuntimeError('Native counterion probe failed: '+str(work))
        print(json.dumps({'mode':mode,'status':'executed','output_dir':str(work)}),flush=True);return receipt
    with ThreadPoolExecutor(max_workers=5) as pool:results=list(pool.map(case,('rest','water_only','na_only','nutrients_only','na_water')))
    (out/'execution-group.json').write_text(json.dumps(results,indent=2)+'\n');return out

if __name__=='__main__':run()
