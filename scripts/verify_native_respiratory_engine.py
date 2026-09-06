#!/usr/bin/env python3
"""Prepare matched whole-engine respiratory observations; --run needs native slot."""
from pathlib import Path
import argparse,hashlib,json,os,resource,shlex,shutil,subprocess,tempfile,time
import xml.etree.ElementTree as ET
from verify_native_respiratory_work import prepare as isolated_sources
ROOT=Path(__file__).resolve().parents[1];SOURCE=ROOT/'data/raw/physiology/biogears';RUNTIME=ROOT/'data/runtime/physiology'
VARIANT=RUNTIME/'variants/whole_body_integrity_gi_absorption'
LIBRARY_PIN='9792d857c47a5907f571a03495fe9f4f1144f114afd72c0451869e1a7049588b'
MANIFEST_PIN='de1aa254b868b4b51e3fbd4f90323e09370a409957d83ac7552a09bb0f3e7125'
ENGINE=SOURCE/'projects/biogears/libBiogears/src/engine/Controller/BioGearsEngine.cpp'
ENGINE_PIN='5c29a09c536a876d569069625bbe4ef2c9e375e3ba12c48c64b94ceed49cbc95'
FILES=['native_respiratory_engine_probe.h','native_respiratory_engine_probe.cpp','native_respiratory_work.h','native_body_ports.h','native_signed_muscle_port.h']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def step_source(text):
    signature='bool BioGearsEngine::AdvanceModelTime(bool appendDataTrack)'
    assert text.count(signature)==1
    start=text.index(signature);end=text.index('\n//-------------------------------------------------------------------------------',start)
    original=text[start:end];method=original.replace(signature,'bool probe_advance(bool appendDataTrack=false)')
    edits=[('  PreProcess();','  before_preprocess();\n  PreProcess();\n  after_preprocess();'),('  Process();','  Process();\n  after_process();')]
    for old,new in edits:
        assert method.count(old)==1,'Unexpected native lifecycle anchor'
        method=method.replace(old,new)
    restored=method
    for old,new in reversed(edits):restored=restored.replace(new,old)
    assert restored.replace('bool probe_advance(bool appendDataTrack=false)',signature)==original
    return method+'\n'
def prepare():
    pins=isolated_sources();assert sha(ENGINE)==ENGINE_PIN
    assert sha(VARIANT/'manifest.json')==MANIFEST_PIN and sha(VARIANT/'libbiogears.so.8.0.0')==LIBRARY_PIN
    for name in FILES:
        p=ROOT/'scripts'/name;assert p.is_file(),'Missing full-engine respiratory probe '+name;pins[str(p.relative_to(ROOT))]=sha(p)
    method=step_source(ENGINE.read_text())
    for anchor in ['  PreProcess();','  Process();']:
        try:step_source(ENGINE.read_text().replace(anchor,anchor+'\n'+anchor,1))
        except AssertionError:pass
        else:raise AssertionError('Duplicate lifecycle anchor accepted')
    comparator_tests()
    return pins,method
def state_tree(path):
    root=ET.parse(path).getroot()
    # Only this collection is reordered: native serializer iterates unordered_set<SESubstance*>.
    entries=[x for x in root if x.tag.rsplit('}',1)[-1]=='ActiveSubstance']
    names=[next(x.text for x in row if x.tag.rsplit('}',1)[-1]=='Name') for row in entries]
    assert len(names)==len(set(names)), 'Duplicate native active substance identity'
    ordered=iter(row for _,row in sorted(zip(names,entries)))
    children=[next(ordered) if x.tag.rsplit('}',1)[-1]=='ActiveSubstance' else x for x in root]
    def content(text):return text if text and text.strip() else ''
    def shape(e,kids=None):return (e.tag,tuple(sorted(e.attrib.items())),(content(e.text),content(e.tail)),tuple(shape(c) for c in (list(e) if kids is None else kids)))
    return shape(root,children)
def compare_state(a,b):
    x,y=state_tree(a),state_tree(b);differences=[]
    def visit(x,y,path):
        if x[:3]!=y[:3]:differences.append({'path':path,'control':x[:3],'observer':y[:3]})
        if len(x[3])!=len(y[3]):differences.append({'path':path,'control_child_count':len(x[3]),'observer_child_count':len(y[3])});return
        for i,(c,d) in enumerate(zip(x[3],y[3])):
            if c!=d:visit(c,d,path+'/'+c[0].rsplit('}',1)[-1]+f'[{i}]')
    visit(x,y,'/BioGearsState')
    return {'raw_bytes_equal':a.read_bytes()==b.read_bytes(),'exact_state_values_equal':not differences,'differences':differences,'comparison':'All XML tags, attributes and text exact; only root named ActiveSubstance siblings sorted because source collection is unordered. No numerical tolerance or removed state fields.'}
def comparator_tests():
    with tempfile.TemporaryDirectory(prefix='respiratory-state-comparison-') as directory:
        a,b=Path(directory)/'a.xml',Path(directory)/'b.xml'
        one='<ActiveSubstance><Name>A</Name><Mass value="1"/></ActiveSubstance>'
        two='<ActiveSubstance><Name>B</Name><Mass value="2"/></ActiveSubstance>'
        original='<BioGearsState>'+one+two+'<System><ReactionTime value="0"/></System></BioGearsState>'
        a.write_text(original);b.write_text(original.replace(one+two,two+one))
        assert compare_state(a,b)['exact_state_values_equal'] and not compare_state(a,b)['raw_bytes_equal']
        b.write_text(original.replace('ReactionTime value="0"','ReactionTime value="1e-313"'))
        assert not compare_state(a,b)['exact_state_values_equal']
        b.write_text(original.replace(one+two,one+one))
        try:state_tree(b)
        except AssertionError:pass
        else:raise AssertionError('Duplicate substance identity accepted')
        a.write_text('<BioGearsState><First/><Second/></BioGearsState>');b.write_text('<BioGearsState><Second/><First/></BioGearsState>')
        assert not compare_state(a,b)['exact_state_values_equal']
        b.write_text('<BioGearsState><First/>unexpected tail<Second/></BioGearsState>')
        assert not compare_state(a,b)['exact_state_values_equal']
def comparison_sources():
    pins={}
    for name in ['include/biogears/cdm/substance/SESubstanceManager.h','src/io/biogears/BioGears.cpp','src/io/biogears/BioGearsPhysiology.cpp']:
        p=SOURCE/'projects/biogears/libBiogears'/name
        assert p.read_bytes()==subprocess.check_output(['git','-C',str(SOURCE),'show','3f16a5fa1dade9c511b88d923606fa51cc35e95d:'+str(p.relative_to(SOURCE))])
        pins[str(p.relative_to(ROOT))]=sha(p)
    return pins
def limits(memory,cpu):
    def cap():resource.setrlimit(resource.RLIMIT_AS,(memory,memory));resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu))
    return cap
def execute(pins,method,retained=None):
    vm=json.loads((VARIANT/'manifest.json').read_text());build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib'
    cwd=build/'projects/biogears/libBiogears';objects=shlex.split((VARIANT/'objects.rsp').read_text())
    assert len(objects)==len(set(objects)) and set(objects)==set(vm['object_sha256'])
    for name in objects:assert sha(cwd/name)==vm['object_sha256'][name],'Native object changed'
    source_manifest=ROOT/'data/derived/systemic/exertion_v3/exercise/native/manifest.json'
    state=Path(json.loads(source_manifest.read_text())['configuration']['state_path'])
    assert sha(state)=='cba7ffb523728d86e8f07b2a30b215b040676c30580e8accda722d175bd687b3'
    out=retained.resolve() if retained else Path(tempfile.mkdtemp(prefix='native-respiratory-engine-',dir=ROOT/'data/derived/audits'));print(out,flush=True)
    for name in FILES:
        if retained:assert (out/name).read_bytes()==(ROOT/'scripts'/name).read_bytes(),'Retained compiled source changed'
        else:(out/name).write_bytes((ROOT/'scripts'/name).read_bytes())
    frozen=out/'input.xml'
    if retained:
        assert (out/'native_respiratory_engine_step.inc').read_text()==method and sha(frozen)==sha(state)
    else:
        (out/'native_respiratory_engine_step.inc').write_text(method);shutil.copyfile(state,frozen)
    report={'source_sha256':pins,'variant':VARIANT.name,'variant_manifest_sha256':MANIFEST_PIN,'library_sha256':LIBRARY_PIN,'engine_source_sha256':ENGINE_PIN,'generated_step_sha256':sha(out/'native_respiratory_engine_step.inc'),'state_sha256':sha(frozen),'state_source_manifest_sha256':sha(source_manifest),'builder_sha256':sha(Path(__file__)),'modes':['baseline','observer','pulse_control','pulse_observer'],'steps_per_mode':8}
    if retained:
        original_report=json.loads((out/'preparation.json').read_text());assert original_report['source_sha256']==pins;report=original_report
        report['analysis_builder_sha256']=sha(Path(__file__))
    else:(out/'preparation.json').write_text(json.dumps(report,indent=2)+'\n')
    command=['c++','-std=c++20','-O0','-g1','-DBIOGEARS_THROW_READONLY_EXCEPTIONS',str(out/'native_respiratory_engine_probe.cpp')]
    for p in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',RUNTIME/'sysroot/usr/include',RUNTIME/'sysroot/usr/include/eigen3',build/'projects/biogears/generated/Release']:
        command+=['-I',str(p)]
    binary=out/'probe';command+=[str(VARIANT/'libbiogears.so.8.0.0'),'-L',str(lib),'-lbiogears_cdm','-o',str(binary)]
    env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':str(VARIANT)+':'+str(lib)}
    def run(cmd,label,work,memory,cpu,wall):
        with (out/(label+'.stdout')).open('w') as stdout,(out/(label+'.stderr')).open('w') as stderr:
            result=subprocess.run(['/usr/bin/time','-v','-o',str(out/(label+'.resources')),'nice','-n','10']+cmd,cwd=work,env=env,stdout=stdout,stderr=stderr,preexec_fn=limits(memory,cpu),timeout=wall)
        if result.returncode:raise RuntimeError(f'{label} failed ({result.returncode}); retained {out}')
    started=time.monotonic()
    if not retained:run(command,'compile',out,1024**3,60,70)
    linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(out/'ldd.txt').write_text(linkage);assert str(VARIANT/'libbiogears.so.8.0.0') in linkage
    frames={};states={}
    for mode in report['modes']:
        work=out/mode
        if not retained:
            work.mkdir()
            for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:(work/name).symlink_to(build/'runtime'/name)
            run([str(binary),mode,str(frozen),str(work/'final.xml')],mode,work,2*1024**3,20,30)
        frames[mode]=[json.loads(line[7:]) for line in (out/(mode+'.stdout')).read_text().splitlines() if line.startswith('RESULT ')]
        assert len(frames[mode])==9 and [r['tick'] for r in frames[mode]]==list(range(9))
        states[mode]=sha(work/'final.xml')
    report['state_comparison_source_sha256']=comparison_sources()
    state_comparisons={}
    for control,observed in [('baseline','observer'),('pulse_control','pulse_observer')]:
        for a,b in zip(frames[control],frames[observed]):assert a['time_s']==b['time_s'] and a['values']==b['values'],'Observer changed native ports'
        state_comparisons[control+'__'+observed]=compare_state(out/control/'final.xml',out/observed/'final.xml')
    max_kcl=max_lung=max_constitutive=max_partition=0.
    for mode in ['observer','pulse_observer']:
        cumulative=0
        for row in frames[mode][1:]:
            w=row['work'];assert w and abs(row['time_s']-frames[mode][0]['time_s']-.02*row['tick'])<1e-8
            assert w['external_pressure_pa']==(10 if mode=='pulse_observer' and row['tick'] in (3,4) else 0)
            cumulative+=w['source_work_j'];assert cumulative==w['cumulative_source_work_j']
            assert abs(w['source_pressure_residual_pa'])<1e-8
            max_kcl=max(max_kcl,abs(w['source_kcl_m3']));max_lung=max(max_lung,abs(w['source_minus_lung_stroke_m3']))
            max_partition=max(max_partition,abs(w['source_work_j']-w['generated_work_j']-w['external_work_j']))
            for side in ['left','right']:max_constitutive=max(max_constitutive,abs(w[side+'.constitutive_residual_m3']))
    assert max_kcl<1e-10 and max_constitutive<1e-10 and max_partition<1e-10
    pulse=frames['pulse_observer'];base=frames['observer']
    assert pulse[2]['values']==base[2]['values']
    delta_volume=pulse[4]['values']['lung_volume_ml']-base[4]['values']['lung_volume_ml']
    delta_flow=pulse[3]['values']['airway_flow_l_per_s']-base[3]['values']['airway_flow_l_per_s']
    assert abs(delta_volume)>1e-6 and abs(delta_flow)>1e-6,'External pressure did not perturb native respiration'
    assert pulse[5]['work']['external_pressure_pa']==0 and pulse[5]['work']['external_work_j']==0
    assert prepare()==(pins,method) and sha(frozen)==report['state_sha256']
    accepted=all(x['exact_state_values_equal'] for x in state_comparisons.values())
    report.update(passed=accepted,respiratory_observation_checks_passed=True,state_comparisons=state_comparisons,analysis_only=retained is not None,wall_s=time.monotonic()-started,frames=frames,final_state_sha256=states,exact_serialized_state_observer_parity=all(states[a]==states[b] for a,b in [('baseline','observer'),('pulse_control','pulse_observer')]),max_kcl_m3=max_kcl,max_constitutive_m3=max_constitutive,max_partition_j=max_partition,max_source_minus_lung_stroke_m3=max_lung,pulse_lung_delta_vs_control_ml=delta_volume,pulse_airway_flow_delta_vs_control_l_s=delta_flow,compile_command=command,executable_sha256=sha(binary),resources={name:(out/(name+'.resources')).read_text() for name in ['compile']+report['modes']},scope='Four matched frozen native human-state loads,8steps each; observational sourcework and10Pa externalpulse only. No physical thorax, actuator energetics or production adapter mutation.')
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k in ['passed','wall_s','exact_serialized_state_observer_parity','max_kcl_m3','max_constitutive_m3','max_partition_j','max_source_minus_lung_stroke_m3','pulse_lung_delta_vs_control_ml','pulse_airway_flow_delta_vs_control_l_s']}))
    if not accepted:raise RuntimeError('Exact native state parity failed; all differences retained in verification.json')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',action='store_true');p.add_argument('--analyze-retained',type=Path);args=p.parse_args();pins,method=prepare()
    if args.run and args.analyze_retained:raise ValueError('Choose new native run or offline retained analysis')
    if args.run or args.analyze_retained:execute(pins,method,args.analyze_retained)
    else:print(json.dumps({'prepared':True,'compiled':False,'source_sha256':pins,'generated_step_sha256':hashlib.sha256(method.encode()).hexdigest()}))
