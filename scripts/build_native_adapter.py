"""Recompile the thin adapter against the already acquired native engine."""
from pathlib import Path
import subprocess
root=Path(__file__).resolve().parents[1];runtime=root/'data/runtime/physiology';source=root/'data/raw/physiology/biogears'
build=runtime/'biogears-build';sysroot=runtime/'sysroot';lib=build/'outputs/Release/lib'
command=['c++','-std=c++20','-O2',str(root/'scripts/native_biogears_rest.cpp')]
for p in [source/'projects/biogears/libBiogears/include',source/'projects/biogears-common/include',sysroot/'usr/include',sysroot/'usr/include/eigen3',build/'projects/biogears/generated/Release']:command+=['-I',str(p)]
command+=['-L',str(lib),f'-Wl,-rpath,{lib}','-lbiogears','-lbiogears_cdm','-o',str(runtime/'native_biogears_rest')]
subprocess.run(command,check=True)
print('Built native adapter with named compartment telemetry')
