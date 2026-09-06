"""Build only the optional thin coupled adapter; original executables unchanged."""
from pathlib import Path
import hashlib
import json
import subprocess
from ihm.native import SOURCE,SOURCE_REVISION,RUNTIME,BASE
from ihm.native.coupled_session import CoupledNativeSession

def main():
    actual=subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip()
    if actual!=SOURCE_REVISION:raise ValueError('Native source revision changed')
    build=RUNTIME/'biogears-build';lib=build/'outputs/Release/lib';sysroot=RUNTIME/'sysroot'
    executable=RUNTIME/CoupledNativeSession.executable_name
    command=['c++','-std=c++20','-O1',str(BASE/'scripts/native_biogears_coupled.cpp')]
    for path in [SOURCE/'projects/biogears/libBiogears/include',SOURCE/'projects/biogears-common/include',sysroot/'usr/include',sysroot/'usr/include/eigen3',build/'projects/biogears/generated/Release']:
        command+=['-I',str(path)]
    temporary=executable.with_suffix('.building')
    command+=['-L',str(lib),f'-Wl,-rpath,{lib}','-lbiogears','-lbiogears_cdm','-o',str(temporary)]
    subprocess.run(command,check=True)
    temporary.replace(executable)
    inputs=[BASE/'scripts'/name for name in CoupledNativeSession.adapter_sources]
    inputs.append(SOURCE/'projects/biogears/libBiogears/src/engine/Controller/BioGearsEngine.cpp')
    manifest={'executable_sha256':hashlib.sha256(executable.read_bytes()).hexdigest(),
        'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'source_revision':actual,'command':command,'original_adapters_modified':False}
    executable.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Built optional coupled native adapter')

if __name__=='__main__':main()
