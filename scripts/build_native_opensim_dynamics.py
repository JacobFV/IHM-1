"""Build the contact forward runner against the retained native OpenSim install."""
from pathlib import Path
import subprocess,json,hashlib
BASE=Path(__file__).resolve().parents[1]
ROOT=BASE/'data/runtime/opensim'

def build():
    project=ROOT/'dynamics-adapter';project.mkdir(parents=True,exist_ok=True)
    (project/'CMakeLists.txt').write_text('cmake_minimum_required(VERSION 3.15)\nproject(NativeDynamics LANGUAGES CXX)\nfind_package(OpenSim REQUIRED)\nadd_executable(native_opensim_dynamics "'+str(BASE/'scripts/native_opensim_dynamics.cpp')+'")\ntarget_link_libraries(native_opensim_dynamics ${OpenSim_LIBRARIES})\ntarget_include_directories(native_opensim_dynamics PRIVATE ${OpenSim_INCLUDE_DIRS})\nset_target_properties(native_opensim_dynamics PROPERTIES CXX_STANDARD 17 RUNTIME_OUTPUT_DIRECTORY "'+str(ROOT)+'")\n')
    commands=[['cmake','-S',str(project),'-B',str(project/'build'),f'-DCMAKE_PREFIX_PATH={ROOT}/install/opensim;{ROOT}/install/simbody',f'-DCMAKE_EXE_LINKER_FLAGS=-Wl,-rpath,{ROOT}/sysroot/usr/lib/aarch64-linux-gnu'],['cmake','--build',str(project/'build'),'-j','2']]
    with (ROOT/'dynamics-build.log').open('w') as log:
        for cmd in commands:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    files=[BASE/'scripts/native_opensim_dynamics.cpp',ROOT/'native_opensim_dynamics']
    for part in ('opensim','simbody'):files+=list((ROOT/f'install/{part}/lib').glob('lib*.so'))
    manifest={'files':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},'source_build':json.loads((ROOT/'build_manifest.json').read_text())['sources']}
    (ROOT/'dynamics_build_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(ROOT/'native_opensim_dynamics')
if __name__=='__main__':build()
