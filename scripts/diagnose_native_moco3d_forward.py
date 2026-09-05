"""Repeat source ODE rollout with explicit accuracy; never changes gait acceptance."""
import argparse,json,os,subprocess,time
from pathlib import Path
from run_native_moco3d import ROOT,RUNTIME,OPENSIM,sha
from verify_native_moco3d import read_sto
import numpy as np

def main():
 p=argparse.ArgumentParser();p.add_argument('candidate',type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--accuracy',type=float,default=1e-8);a=p.parse_args()
 if not 1e-10<=a.accuracy<=1e-5:raise ValueError('Accuracy outside bounded diagnostic range')
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);candidate=a.candidate.resolve();manifest=json.loads((candidate/'manifest.json').read_text());model=candidate/'optimized_plant.osim';solution=candidate/('example3DWalking_'+('muscle' if manifest['mode']=='muscle' else 'torque')+'_driven_tracking_solution.sto')
 cpp=out/'forward.cpp';cpp.write_text('''#include <OpenSim/OpenSim.h>
#include <OpenSim/Moco/osimMoco.h>
#include <iostream>
using namespace OpenSim;
int main(int argc,char**argv) { try {
 Model model(argv[1]); auto state=model.initSystem(); MocoTrajectory solution(argv[2]);
 for(const auto& a:model.getComponentList<ScalarActuator>())
 std::cout<<"ACTUATOR "<<a.getAbsolutePathString()<<" override="<<a.isActuationOverridden(state)<<std::endl;
 auto forward=simulateTrajectoryWithTimeStepping(solution,model,std::stod(argv[3]));forward.write("forward_validation.sto");
 STOFileAdapter::write(analyzeMocoTrajectory<double>(model,forward,{".*actuation",".*activation",".*tendon_force"}),"forward_actuation.sto");
 return 0;}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 10;}}
''')
 cmd=manifest['compile_command'].copy();cmd[cmd.index(str(candidate/'executed_example.cpp'))]=str(cpp);cmd[-1]=str(out/'forward')
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':':'.join(map(str,[RUNTIME,RUNTIME/'casadi-wheel/casadi',OPENSIM/'install/opensim/lib',OPENSIM/'install/simbody/lib',OPENSIM/'sysroot/usr/lib/aarch64-linux-gnu']))}
 with (out/'compile.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
 deps=subprocess.check_output(['ldd',str(out/'forward')],env=env,text=True);(out/'ldd.txt').write_text(deps)
 record={'candidate':str(candidate),'model_sha256':sha(model),'solution_sha256':sha(solution),'accuracy':a.accuracy,'compile_command':cmd,'executable_sha256':sha(out/'forward'),'source_sha256':sha(cpp),'dependency_sha256':{str(Path(s.split('=>')[1].split(' (')[0].strip()).resolve()):sha(Path(s.split('=>')[1].split(' (')[0].strip())) for s in deps.splitlines() if '=>' in s}}
 (out/'manifest.json').write_text(json.dumps(record,indent=2));start=time.monotonic()
 with (out/'runner.log').open('w') as f:r=subprocess.run([str(out/'forward'),str(model),str(solution),str(a.accuracy)],cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=180)
 report={'exit_code':r.returncode,'wall_s':time.monotonic()-start}
 if r.returncode==0:
  names,x,meta=read_sto(solution);other,y,_=read_sto(out/'forward_validation.sto');states=[n for n in names if n.endswith(('/value','/speed','/activation'))];controls=[n for n in names if n.startswith('/forceset/') and n.count('/')==2]
  report['maximum_initial_state_error']=max(abs(x[0,names.index(n)]-y[0,other.index(n)]) for n in states)
  report['maximum_control_interpolation_error']=max(float(np.max(abs(np.interp(y[:,0],x[:,0],x[:,names.index(n)])-y[:,other.index(n)]))) for n in controls)
  report['coordinate_errors']={n:float(np.max(abs(np.interp(y[:,0],x[:,0],x[:,names.index(n)])-y[:,other.index(n)]))) for n in states if n.endswith('/value')}
  report['all_actuator_overrides_disabled']='override=1' not in (out/'runner.log').read_text()
 (out/'diagnosis.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
