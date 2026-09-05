"""Run the primary 3D contact gait example with explicit solve budgets and receipts."""
import argparse,difflib,hashlib,json,os,shlex,shutil,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];RUNTIME=ROOT/'data/runtime/moco3d';OPENSIM=ROOT/'data/runtime/opensim';SOURCE=ROOT/'data/raw/mechanics/opensim-core';EXAMPLE=SOURCE/'OpenSim/Examples/Moco/example3DWalking'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--mode',choices=['preflight','torque','muscle'],default='preflight');p.add_argument('--max-iterations',type=int,default=200);p.add_argument('--timeout-s',type=int,default=900);p.add_argument('--guess');p.add_argument('--mesh-interval',type=float,default=.02);p.add_argument('--continuation-guess');p.add_argument('--original-backend',action='store_true');args=p.parse_args()
 if not 1<=args.max_iterations<=1000 or not 1<=args.timeout_s<=3600:p.error('bounded iterations/time required')
 if args.mesh_interval not in [.02,.01]:p.error('Only source mesh or declared half-mesh refinement allowed')
 if args.continuation_guess and args.mode!='torque':p.error('Continuation currently limited to identical torque problem')
 card=json.loads((ROOT/'data/sources/moco-3d-contact-gait.json').read_text())
 for name,digest in card['source_files_sha256'].items():
  if sha(EXAMPLE/name)!=digest:raise RuntimeError('Pinned primary source changed: '+name)
 out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=True)
 if any(out.iterdir()):raise ValueError('Fresh artifact directory required')
 for f in EXAMPLE.iterdir():
  if f.is_file() and f.suffix not in ['.cpp','.txt']:shutil.copyfile(f,out/f.name)
 if args.guess:shutil.copyfile(args.guess,out/'example3DWalking_torque_driven_tracking_solution.sto')
 if args.continuation_guess:shutil.copyfile(args.continuation_guess,out/'continuation_guess.sto')
 before=(EXAMPLE/'example3DWalking.cpp').read_text();after=before.replace('#include <format>','#include <format>\n#include <iostream>\n#include <stdexcept>')
 after=after.replace('int main() {','int main(int argc, char** argv) {\n    std::cout << "IHM_BACKEND_AVAILABLE " << MocoCasADiSolver::isAvailable() << std::endl;\n    if (!MocoCasADiSolver::isAvailable()) return 7;\n    std::string mode = argc > 1 ? argv[1] : "preflight";')
 after=after.replace('    runTrackingStudy(model, false);','''    auto& state = model.initSystem();
    std::cout << "IHM_MODEL mass_kg=" << model.getTotalMass(state) << " bodies=" << model.getBodySet().getSize() << " muscles=" << model.getMuscles().getSize() << " coordinates=" << model.getCoordinateSet().getSize() << " contacts=" << contactForceSet.getSize() << std::endl;
    model.print("prepared_contact_model.osim");
    if (mode == "preflight") return 0;
    if (mode == "torque") runTrackingStudy(model, false);''')
 after=after.replace('    runTrackingStudy(model, true);','    if (mode == "muscle") runTrackingStudy(model, true);')
 after=after.replace('    solver.set_optim_convergence_tolerance(1e-2);',f'    solver.set_parallel(4);\n    solver.set_output_interval(25);\n    solver.set_optim_max_iterations({args.max_iterations});\n    solver.set_optim_convergence_tolerance(1e-2);')
 after=after.replace('track.set_mesh_interval(0.02);',f'track.set_mesh_interval({args.mesh_interval});')
 if args.continuation_guess:after=after.replace('    solver.setGuess(guess);','    solver.setGuessFile("continuation_guess.sto");')
 after=after.replace('    MocoSolution solution = study.solve();','''    study.print("executed_study.omoco");
    Model prepared = modelProcessor.process(); prepared.initSystem(); prepared.print("optimized_plant.osim");
    MocoSolution solution = study.solve();
    const bool converged = solution.success();
    if (!converged) solution.unseal();
    std::cout << "IHM_SOLVE success=" << converged << " status=" << solution.getStatus() << " iterations=" << solution.getNumIterations() << std::endl;
    if (!converged) { solution.write("failed_iterate.sto"); throw std::runtime_error("Moco solve did not converge; failed iterate retained"); }''')
 after=after.replace('    study.visualize(solution);','''    auto forces = analyzeMocoTrajectory<double>(modelSolution, solution,
        {".*activation", ".*fiber_length", ".*fiber_velocity", ".*tendon_force", ".*actuation"});
    STOFileAdapter::write(forces, "optimized_muscle_and_actuator_outputs.sto");
    // Independent ODE rollout from optimized initial state and controls.
    auto forward = simulateTrajectoryWithTimeStepping(solution, modelSolution, 1e-6);
    forward.write("forward_validation.sto");
    auto forwardForces = analyzeMocoTrajectory<double>(modelSolution, forward,
        {".*activation", ".*fiber_length", ".*fiber_velocity", ".*tendon_force", ".*actuation"});
    STOFileAdapter::write(forwardForces, "forward_muscle_and_actuator_outputs.sto");
    auto forwardGRF = createExternalLoadsTableForGait(modelSolution, forward, contactForcesRight, contactForcesLeft);
    STOFileAdapter::write(forwardGRF, "forward_ground_reactions.sto");''')
 after=after.replace('int main(int argc, char** argv) {','int main(int argc, char** argv) { try {',1)
 after=after.replace('    return EXIT_SUCCESS;\n}', '    return EXIT_SUCCESS;\n} catch (const std::exception& e) { std::cerr << "IHM_NATIVE_FAILURE " << e.what() << std::endl; return 10; }\n}',1)
 cpp=out/'executed_example.cpp';cpp.write_text(after);(out/'source_changes.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='primary/example3DWalking.cpp',tofile='executed_example.cpp')))
 binary=out/'moco3d';selected=OPENSIM/'install/opensim/lib' if args.original_backend else RUNTIME
 command=['c++','-std=c++20','-O2',str(cpp),'-I',str(OPENSIM/'install/opensim/include'),'-I',str(OPENSIM/'install/simbody/include/simbody'),'-L',str(selected),'-L',str(OPENSIM/'install/opensim/lib'),'-L',str(OPENSIM/'install/simbody/lib')]
 command+=['-I',str(SOURCE/'Vendors/lepton/include')]
 command+=['-l'+x for x in ['osimMoco','osimTools','osimAnalyses','osimActuators','osimSimulation','osimCommon','SimTKsimbody','SimTKmath','SimTKcommon']]
 command+=['-Wl,-rpath-link,'+str(RUNTIME/'casadi-wheel/casadi'),'-Wl,-rpath-link,'+str(OPENSIM/'sysroot/usr/lib/aarch64-linux-gnu'),'-o',str(binary)]
 with (out/'compile.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','LD_LIBRARY_PATH':':'.join(map(str,[selected,RUNTIME/'casadi-wheel/casadi',OPENSIM/'install/opensim/lib',OPENSIM/'install/simbody/lib',OPENSIM/'sysroot/usr/lib/aarch64-linux-gnu']))}
 linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True);(out/'ldd.txt').write_text(linkage);assert 'not found' not in linkage
 deps={str(Path(s.split('=>')[1].split(' (')[0].strip()).resolve()):sha(Path(s.split('=>')[1].split(' (')[0].strip())) for s in linkage.splitlines() if '=>' in s}
 record={'runner_sha256':sha(Path(__file__)), 'executed_input_sha256':{f.name:sha(out/f.name) for f in EXAMPLE.iterdir() if f.is_file() and f.suffix not in ['.cpp','.txt']},'casadi_plugin_sha256':{f.name:sha(f) for f in (RUNTIME/'casadi-wheel/casadi').glob('*.so*') if f.is_file()},'mesh_interval_s':args.mesh_interval,'continuation_guess_sha256':sha(Path(args.continuation_guess)) if args.continuation_guess else None,'mode':args.mode,'max_iterations':args.max_iterations,'timeout_s':args.timeout_s,'source_revision':subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip(),'primary_file_sha256':{f.name:sha(f) for f in EXAMPLE.iterdir() if f.is_file()},'executed_source_sha256':sha(cpp),'executable_sha256':sha(binary),'compile_command':command,'dependency_sha256':deps,'model_scope':'Source 3D contact gait subject; canonical-body registration not performed','guess_sha256':sha(Path(args.guess)) if args.guess else None}
 (out/'manifest.json').write_text(json.dumps(record,indent=2)+'\n');start=time.monotonic()
 with (out/'runner.log').open('w') as log:
  try:r=subprocess.run([str(binary),args.mode],cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout_s);status={'exit_code':r.returncode,'timed_out':False}
  except subprocess.TimeoutExpired:status={'exit_code':None,'timed_out':True}
 status.update(elapsed_wall_s=time.monotonic()-start,outputs={f.name:sha(f) for f in out.iterdir() if f.is_file() and f.suffix in ['.sto','.osim','.omoco']})
 status['source_unchanged']=all(sha(EXAMPLE/f)==expected for f,expected in record['primary_file_sha256'].items())
 status['inputs_unchanged']=all(sha(out/name)==expected for name,expected in record['executed_input_sha256'].items())
 status['native_success']=status['exit_code']==0 and not status['timed_out'] and status['inputs_unchanged'] and status['source_unchanged']
 (out/'execution.json').write_text(json.dumps(status,indent=2)+'\n');print(out);print(json.dumps(status,indent=2))
 if not status['native_success']:raise SystemExit(1)
if __name__=='__main__':main()
