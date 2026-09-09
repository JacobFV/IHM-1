"""Build the native PolynomialPathFitter driver against the local OpenSim install.

`scripts/build_native_opensim.py` already builds and installs OpenSim and Simbody
under ``data/runtime/opensim/install``; ``libosimActuators.so`` in that install
exports ``OpenSim::PolynomialPathFitter`` (97 symbols).  The fitter was never
unavailable -- only undriven.  This script adds the driver, in its own CMake
project so it does not disturb the two adapters already built there.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE / 'data/runtime/opensim'
INSTALL = ROOT / 'install'
PROJECT = ROOT / 'path-fitter'
BINARY = ROOT / 'native_polynomial_path_fit'
SOURCE = BASE / 'scripts/native_polynomial_path_fit.cpp'

CMAKE = '''cmake_minimum_required(VERSION 3.15)
project(NativePathFitter LANGUAGES CXX)
find_package(OpenSim REQUIRED)
add_executable(native_polynomial_path_fit "{source}")
target_link_libraries(native_polynomial_path_fit ${{OpenSim_LIBRARIES}})
target_include_directories(native_polynomial_path_fit PRIVATE ${{OpenSim_INCLUDE_DIRS}})
set_target_properties(native_polynomial_path_fit PROPERTIES CXX_STANDARD 17
    RUNTIME_OUTPUT_DIRECTORY "{out}")
'''


def run(args, name, cwd=BASE):
    print(name, flush=True)
    log = ROOT / (name + '.log')
    with log.open('w') as handle:
        subprocess.run([str(a) for a in args], cwd=cwd, stdout=handle,
                       stderr=subprocess.STDOUT, check=True)


def build(jobs=6):
    if not (INSTALL / 'opensim/lib/libosimActuators.so').exists():
        raise SystemExit('OpenSim is not installed; run scripts/build_native_opensim.py first')
    symbols = subprocess.check_output(
        ['nm', '-D', '--defined-only', str(INSTALL / 'opensim/lib/libosimActuators.so')],
        text=True)
    if 'PolynomialPathFitter' not in symbols:
        raise SystemExit('The installed libosimActuators.so exports no PolynomialPathFitter')

    PROJECT.mkdir(parents=True, exist_ok=True)
    (PROJECT / 'CMakeLists.txt').write_text(CMAKE.format(source=SOURCE, out=ROOT))
    run(['cmake', '-S', PROJECT, '-B', PROJECT / 'build',
         '-DCMAKE_BUILD_TYPE=Release',
         f'-DCMAKE_PREFIX_PATH={INSTALL}/opensim;{INSTALL}/simbody'],
        'path-fitter-configure')
    run(['cmake', '--build', PROJECT / 'build', '-j', str(jobs)], 'path-fitter-build')
    if not BINARY.exists():
        raise SystemExit('Build reported success but produced no binary')

    manifest = {
        'binary': str(BINARY.relative_to(BASE)),
        'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'binary_sha256': hashlib.sha256(BINARY.read_bytes()).hexdigest(),
        'opensim_install': str(INSTALL.relative_to(BASE)),
    }
    (ROOT / 'path_fitter_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--jobs', type=int, default=6)
    build(parser.parse_args().jobs)
