#!/usr/bin/env python3
"""Compile only the thin adapter against the already built upstream libraries."""
from pathlib import Path
import subprocess
base=Path(__file__).resolve().parents[1];source=base/'data/raw/physiology/biogears';root=base/'data/runtime/physiology';build=root/'biogears-build';lib=build/'outputs/Release/lib';sysroot=root/'sysroot'
subprocess.run(['c++','-std=c++20','-O2',str(base/'scripts/native_biogears_rest.cpp'),'-I',str(source/'projects/biogears/libBiogears/include'),'-I',str(source/'projects/biogears-common/include'),'-I',str(sysroot/'usr/include'),'-I',str(sysroot/'usr/include/eigen3'),'-I',str(build/'projects/biogears/generated/Release'),'-L',str(lib),f'-Wl,-rpath,{lib}','-lbiogears','-lbiogears_cdm','-o',str(root/'native_biogears_rest')],check=True)
