"""Fresh isolated native observer/kinematics build; compilation explicitly gated."""
from pathlib import Path
import argparse,json,os,shlex,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.source_skin import file_sha256
POSE=r'''
   if(command=="hair_pose"){
    std::string body;int count;in>>body>>count;if(!in||count<1||count>8)throw std::runtime_error("bounded hair pose coordinates required");
    SimTK::State candidate=state;
    for(int i=0;i<count;i++){std::string name;double value;in>>name>>value;const auto& c=model.getCoordinateSet().get(name);if(!in||!std::isfinite(value)||c.isDependent(state)||c.getLocked(state)||value<c.getRangeMin()||value>c.getRangeMax())throw std::runtime_error("invalid hair pose coordinate");c.setValue(candidate,value,false);}
    std::string extra;if(in>>extra)throw std::runtime_error("trailing hair pose data");model.assemble(candidate);model.realizePosition(candidate);
    std::cout<<"@IHM {\"kind\":\"hair_pose\",\"body\":";str(std::cout,body);std::cout<<",\"transform_ground\":";matrix(std::cout,model.getBodySet().get(body).getTransformInGround(candidate));std::cout<<"}"<<std::endl;continue;
   }
'''

def prepare(output,compile=False):
    output=Path(output).resolve()
    if output.exists() or not output.is_relative_to(ROOT):raise ValueError('Fresh owned observer build required')
    old=ROOT/json.loads((ROOT/'data/runtime/mechanical-stream/latest.json').read_bytes())['build'];old_manifest=json.loads((old/'manifest.json').read_bytes());output.mkdir(parents=True);dependencies={}
    for path in [old/'native_mechanical_stream.cpp',*old.glob('*.h')]:
        raw=path.read_bytes();digest=file_sha256(raw)
        if old_manifest['files'].get(str(path.relative_to(ROOT)))!=digest:raise ValueError('Frozen native source/header hash mismatch')
        dependencies[str(path.relative_to(ROOT))]=digest;(output/path.name).write_bytes(raw)
    source=output/'native_mechanical_stream.cpp';text=source.read_text();anchor='   if(command=="close")break;';field='   o<<",\\"optimal_fiber_length_m\\":";'
    if text.count(anchor)!=1 or text.count(field)!=1:raise ValueError('Native observer anchors changed')
    text=text.replace(anchor,anchor+POSE).replace(field,'   o<<",\\"passive_energy_j\\":";num(o,m.getMusclePotentialEnergy(state));\n'+field);source.write_text(text)
    runtime=ROOT/'data/runtime/opensim';flags=['-std=c++20','-O0','-DSWIG_PYTHON']
    for part in ('install/opensim/include','install/opensim/include/OpenSim','install/simbody/include/simbody'):flags+=['-isystem',str(runtime/part)]
    previous=shlex.split((runtime/'dynamics-adapter/build/CMakeFiles/native_opensim_dynamics.dir/link.txt').read_text());libraries=[v for v in previous if v.endswith('.so') or '.so.' in v or v.startswith('-Wl,') or v.startswith('-l')]
    exe=output/'native_mechanical_stream';command=['prlimit','--as=3221225472','--','nice','-n','10','taskset','-c','0','c++',*flags,str(source),'-o',str(exe),*libraries,'-ldl']
    report={'schema':'ihm.hair-native-observer-build.v1','complete':False,'source_dependencies':dependencies,'command':command,'changes':['Read-only raw muscle passive potential emission','Read-only copied-state bounded articulated pose transform query'],'force_laws_unchanged':True}
    (output/'command.json').write_text(json.dumps(command,indent=2)+'\n')
    if compile:
        started=time.monotonic()
        with (output/'compile.log').open('w') as log:result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),timeout=45)
        if result.returncode:raise RuntimeError('Observer compile failed: '+str(output))
        report.update(complete=True,compile_wall_s=time.monotonic()-started)
    report['files']={str(p.relative_to(ROOT)):file_sha256(p) for p in [source,*output.glob('*.h'),*([exe] if exe.exists() else [])]};(output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');return report

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);parser.add_argument('--compile',action='store_true');args=parser.parse_args();r=prepare(args.output,args.compile);print(json.dumps({'complete':r['complete'],'force_laws_unchanged':r['force_laws_unchanged']}))
