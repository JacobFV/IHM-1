"""Build an isolated numerical-convergence experiment; original engine stays intact."""
from pathlib import Path
import subprocess,shlex,json,hashlib,argparse,difflib
BASE=Path(__file__).resolve().parents[1];ROOT=BASE/'data/runtime/opensim';SOURCE=BASE/'data/raw/mechanics/opensim-core/OpenSim/Simulation/Model/GeometryPath.cpp'
def build(iterations,tolerance,cache_fix=False):
    tag=f'wrap_{iterations}_{tolerance:g}'+('_cache' if cache_fix else '');out=ROOT/'variants'/tag;out.mkdir(parents=True,exist_ok=True)
    original=SOURCE.read_text();patched=original.replace('const int maxIterations = wrapSetSize < 2 ? 1 : 8;',f'const int maxIterations = wrapSetSize < 2 ? 1 : {iterations};').replace('std::abs(length - last_length) < 0.0005',f'std::abs(length - last_length) < {tolerance}')
    if original==patched and not cache_fix:raise ValueError('Expected upstream patch anchors absent')
    cpp=out/'GeometryPath.cpp';cpp.write_text(patched)
    obj=out/'GeometryPath.cpp.o';build=ROOT/'opensim-build/OpenSim/Simulation';metadata=build/'CMakeFiles/osimSimulation.dir'
    flags={line.split(' = ',1)[0]:shlex.split(line.split(' = ',1)[1]) for line in (metadata/'flags.make').read_text().splitlines() if ' = ' in line}
    command=['/usr/bin/c++','-I',str(SOURCE.parent),*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],*flags['CXX_FLAGS'],'-c',str(cpp),'-o',str(obj)]
    with (out/'compile.log').open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
    link=shlex.split((metadata/'link.txt').read_text());link[link.index('-o')+1]=str(out/'libosimSimulation.so')
    anchor='CMakeFiles/osimSimulation.dir/Model/GeometryPath.cpp.o';link[link.index(anchor)]=str(obj)
    patch_sources={}
    if cache_fix:
        wrap_source=SOURCE.parent.parent/'Wrap/PathWrapPoint.cpp';wrap_cpp=out/'PathWrapPoint.cpp';wrap_obj=out/'PathWrapPoint.cpp.o'
        anchor='    markCacheVariableValid(s, _location);'
        replacement=anchor+'\n    markCacheVariableInvalid(s, "location");\n    markCacheVariableInvalid(s, "velocity");\n    markCacheVariableInvalid(s, "acceleration");'
        contents=wrap_source.read_text()
        if contents.count(anchor)!=1:raise ValueError('Cache patch anchor differs')
        wrap_cpp.write_text(contents.replace(anchor,replacement))
        (out/'cache_fix.patch').write_text(''.join(difflib.unified_diff(contents.splitlines(True),wrap_cpp.read_text().splitlines(True),fromfile='a/OpenSim/Simulation/Wrap/PathWrapPoint.cpp',tofile='b/OpenSim/Simulation/Wrap/PathWrapPoint.cpp')))
        compile_wrap=['/usr/bin/c++','-I',str(wrap_source.parent),*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],*flags['CXX_FLAGS'],'-c',str(wrap_cpp),'-o',str(wrap_obj)]
        with (out/'compile-cache.log').open('w') as f:subprocess.run(compile_wrap,stdout=f,stderr=subprocess.STDOUT,check=True)
        link[link.index('CMakeFiles/osimSimulation.dir/Wrap/PathWrapPoint.cpp.o')]=str(wrap_obj)
        patch_sources={'PathWrapPoint_original_sha256':hashlib.sha256(wrap_source.read_bytes()).hexdigest(),'PathWrapPoint_patched_sha256':hashlib.sha256(wrap_cpp.read_bytes()).hexdigest()}
    with (out/'link.log').open('w') as f:subprocess.run(link,cwd=build,stdout=f,stderr=subprocess.STDOUT,check=True)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps({'variant':tag,'source_sha256':sha(SOURCE),'patched_sha256':sha(cpp),'library_sha256':sha(out/'libosimSimulation.so'),'cache_invalidation_fix':cache_fix,**patch_sources,'max_iterations':iterations,'length_convergence_tolerance_m':tolerance,'changes':'Separate cache-invalidation correction and explicitly reported multiwrap iteration/tolerance settings; raw source and original installation unchanged'},indent=2)+'\n')
    print(out)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--iterations',type=int,default=200);p.add_argument('--tolerance',type=float,default=1e-10);p.add_argument('--cache-fix',action='store_true');a=p.parse_args();build(a.iterations,a.tolerance,a.cache_fix)
