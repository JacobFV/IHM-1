#!/usr/bin/env python3
"""Build the collected BioGears source natively without system package installation."""
from pathlib import Path
import os,subprocess
base=Path(__file__).resolve().parents[1]
source=base/'data/raw/physiology/biogears'
root=base/'data/runtime/physiology'; deps=root/'debs'; sysroot=root/'sysroot'; build=root/'biogears-build'
for p in [deps,sysroot,build]:p.mkdir(parents=True,exist_ok=True)
def run(args,log,cwd=base,env=None):
    print('Running:', ' '.join(map(str,args)),flush=True)
    with (root/log).open('w') as f:subprocess.run(list(map(str,args)),cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
run(['apt-get','download','libeigen3-dev','libxerces-c-dev','libxerces-c3.2t64','liblog4cpp5-dev','liblog4cpp5v5','xsdcxx'],'dependencies.log',deps)
for p in deps.glob('*.deb'):subprocess.run(['dpkg-deb','-x',str(p),str(sysroot)],check=True)
arch=subprocess.check_output(['dpkg-architecture','-qDEB_HOST_MULTIARCH'],text=True).strip()
if not (sysroot/'lib').exists():(sysroot/'lib').symlink_to('usr/lib/'+arch)
env=dict(os.environ,LD_LIBRARY_PATH=str(sysroot/'lib'))
run(['cmake','-S',source,'-B',build,'-DCMAKE_BUILD_TYPE=Release',f'-DCMAKE_PREFIX_PATH={sysroot}/usr',f'-DBiogears_EXTERNAL={sysroot}',f'-DCodeSynthesis_EXECUTABLE={sysroot}/usr/bin/xsdcxx','-DBiogears_BUILD_DOCUMENTATION=OFF','-DBiogears_BUILD_IO_LIBRARY=OFF'],'configure-reproducible.log',env=env)
generated=build/'projects/biogears/generated'
for name in ['cdm','biogears']:(generated/'Release/biogears/schema'/name).mkdir(parents=True,exist_ok=True)
if not (generated/'biogears').exists():(generated/'biogears').symlink_to('Release/biogears')
# GCC 13 reveals an upstream missing <iterator> include. Supply it only to
# DataTrack.cpp through the generated recipe, leaving archived source intact.
recipe=build/'projects/biogears/libBiogears/CMakeFiles/libbiogears.dir/build.make'
text=recipe.read_text();text='\n'.join(line.replace('/usr/bin/c++ ','/usr/bin/c++ -include iterator ') if ' -c ' in line and 'DataTrack.cpp' in line else line for line in text.split('\n'));recipe.write_text(text)
run(['cmake','--build',build,'--target','libbiogears','prepare_runtime_dir','-j','16'],'build-reproducible.log')
lib=build/'outputs/Release/lib'
run(['c++','-std=c++20','-O2',base/'scripts/native_biogears_rest.cpp',
     '-I',source/'projects/biogears/libBiogears/include','-I',source/'projects/biogears-common/include',
     '-I',sysroot/'usr/include','-I',sysroot/'usr/include/eigen3','-I',generated/'Release',
     '-L',lib,f'-Wl,-rpath,{lib}','-lbiogears','-lbiogears_cdm','-o',root/'native_biogears_rest'],'adapter-reproducible.log')
print(root/'native_biogears_rest')
