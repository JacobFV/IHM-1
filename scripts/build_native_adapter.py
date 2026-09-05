"""Recompile the thin adapter against the already acquired native engine."""
from pathlib import Path
import subprocess
import hashlib
import json
root=Path(__file__).resolve().parents[1];runtime=root/'data/runtime/physiology';source=root/'data/raw/physiology/biogears'
build=runtime/'biogears-build';sysroot=runtime/'sysroot';lib=build/'outputs/Release/lib'
for name in ('native_biogears_rest', 'native_biogears_stream'):
    command=['c++','-std=c++20','-O2',str(root/f'scripts/{name}.cpp')]
    for p in [source/'projects/biogears/libBiogears/include',source/'projects/biogears-common/include',sysroot/'usr/include',sysroot/'usr/include/eigen3',build/'projects/biogears/generated/Release']:command+=['-I',str(p)]
    temporary=runtime/f'{name}.building'
    command+=['-L',str(lib),f'-Wl,-rpath,{lib}','-lbiogears','-lbiogears_cdm','-o',str(temporary)]
    subprocess.run(command,check=True)
    temporary.replace(runtime/name)
    inputs=[root/f'scripts/{name}.cpp']
    if name=='native_biogears_stream':inputs.append(root/'scripts/native_body_ports.h')
    manifest={'executable_sha256':hashlib.sha256((runtime/name).read_bytes()).hexdigest(),
              'source_sha256':{str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
              'command':command}
    (runtime/f'{name}.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Built {name}')
