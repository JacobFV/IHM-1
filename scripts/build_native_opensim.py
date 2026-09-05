"""Build official OpenSim/Simbody locally for the host, without system installs."""
from pathlib import Path
import subprocess,os,argparse,json,hashlib
BASE=Path(__file__).resolve().parents[1];RAW=BASE/'data/raw/mechanics';ROOT=BASE/'data/runtime/opensim'
PINS={'opensim-core':('https://github.com/opensim-org/opensim-core.git','86b30588374650fbaf012a345a836a64f6855522'),'simbody':('https://github.com/simbody/simbody.git','59c6e7b89b3bdf266a2f3d54c611599d964205f7')}
def run(args,name,cwd=BASE):
    print(name,flush=True)
    with (ROOT/(name+'.log')).open('w') as f:subprocess.run([str(a) for a in args],cwd=cwd,stdout=f,stderr=subprocess.STDOUT,check=True)
def build(stage):
    ROOT.mkdir(parents=True,exist_ok=True);RAW.mkdir(parents=True,exist_ok=True)
    if stage in ('all','dependencies'):
        for name,(url,pin) in PINS.items():
            source=RAW/name
            if not source.exists():
                run(['git','clone','--depth','1',url,source],name+'-clone')
                run(['git','-C',source,'fetch','--depth','1','origin',pin],name+'-fetch')
                run(['git','-C',source,'checkout','--detach','FETCH_HEAD'],name+'-checkout')
            if subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()!=pin:raise ValueError('Unexpected source identity '+name)
        debs=ROOT/'debs';debs.mkdir(exist_ok=True);sysroot=ROOT/'sysroot';sysroot.mkdir(exist_ok=True)
        if len(list(debs.glob('*.deb')))<5:run(['apt-get','download','liblapack-dev','liblapack3','libblas-dev','libblas3','libgfortran5'],'dependencies',debs)
        for p in debs.glob('*.deb'):run(['dpkg-deb','-x',p,sysroot],'extract-'+p.stem)
        arch=subprocess.check_output(['dpkg-architecture','-qDEB_HOST_MULTIARCH'],text=True).strip();libs=sysroot/'usr/lib'/arch
        run(['cmake','-S',RAW/'simbody','-B',ROOT/'simbody-build','-DCMAKE_BUILD_TYPE=Release',f'-DCMAKE_INSTALL_PREFIX={ROOT}/install/simbody','-DBUILD_VISUALIZER=OFF','-DBUILD_EXAMPLES=OFF','-DBUILD_TESTING=OFF',f'-DBUILD_USING_OTHER_LAPACK={libs}/lapack/liblapack.so;{libs}/blas/libblas.so'],'simbody-configure')
        run(['cmake','--build',ROOT/'simbody-build','--target','install','-j','16'],'simbody-build')
    if stage in ('all','engine'):
        run(['cmake','-S',RAW/'opensim-core','-B',ROOT/'opensim-build','-DCMAKE_BUILD_TYPE=Release',f'-DCMAKE_INSTALL_PREFIX={ROOT}/install/opensim',f'-DSIMBODY_HOME={ROOT}/install/simbody','-DBUILD_TESTING=OFF','-DBUILD_API_ONLY=ON','-DBUILD_JAVA_WRAPPING=OFF','-DBUILD_PYTHON_WRAPPING=OFF','-DOPENSIM_WITH_CASADI=OFF','-DOPENSIM_COPY_DEPENDENCIES=OFF','-DOPENSIM_C3D_PARSER=None'],'opensim-configure')
        run(['cmake','--build',ROOT/'opensim-build','--target','install','-j','16'],'opensim-build')
    if stage in ('all','adapter'):
        (ROOT/'install/opensim/bin').mkdir(exist_ok=True) # API-only install omits bin but OpenSimConfig validates it.
        project=ROOT/'adapter';project.mkdir(exist_ok=True)
        (project/'CMakeLists.txt').write_text('cmake_minimum_required(VERSION 3.15)\nproject(NativeOpenSim LANGUAGES CXX)\nfind_package(OpenSim REQUIRED)\nadd_executable(native_opensim_export "'+str(BASE/'scripts/native_opensim_export.cpp')+'")\ntarget_link_libraries(native_opensim_export ${OpenSim_LIBRARIES})\ntarget_include_directories(native_opensim_export PRIVATE ${OpenSim_INCLUDE_DIRS})\nset_target_properties(native_opensim_export PROPERTIES CXX_STANDARD 17 RUNTIME_OUTPUT_DIRECTORY "'+str(ROOT)+'")\n')
        run(['cmake','-S',project,'-B',project/'build',f'-DCMAKE_PREFIX_PATH={ROOT}/install/opensim;{ROOT}/install/simbody',f'-DCMAKE_EXE_LINKER_FLAGS=-Wl,-rpath,{ROOT}/sysroot/usr/lib/aarch64-linux-gnu'],'adapter-configure')
        run(['cmake','--build',project/'build','-j','4'],'adapter-build')
    manifest={'sources':{name:{'url':url,'revision':pin} for name,(url,pin) in PINS.items()},'architecture':os.uname().machine,'files':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for pattern in ['debs/*.deb','install/opensim/lib/lib*.so','install/simbody/lib/lib*.so','native_opensim_export'] for p in ROOT.glob(pattern) if p.is_file()}}
    (ROOT/'build_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['all','dependencies','engine','adapter'],default='all');build(p.parse_args().stage)
