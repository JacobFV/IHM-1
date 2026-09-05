#!/usr/bin/env python3
"""Build isolated original-equation variant removing a proven out-of-bounds access."""
from pathlib import Path
import argparse,difflib,hashlib,json,shlex,subprocess
BASE=Path(__file__).resolve().parents[1]
RUNTIME=BASE/'data/runtime/physiology'
def build(heatflux=False,thermal_units=False):
    heatflux=heatflux or thermal_units
    source=BASE/'data/raw/physiology/biogears';revision=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    if revision!='3f16a5fa1dade9c511b88d923606fa51cc35e95d':raise RuntimeError('Unexpected source revision')
    original=source/'projects/biogears/libBiogears/src/engine/Systems/Saturation.cpp'
    text=original.read_text();bad='    if (std::abs(x(3)) < approxZero) {\n      x(3) = 0.0;\n    }\n'
    if text.count(bad)!=1:raise RuntimeError('Expected exact known defect once')
    corrected=text.replace(bad,'');variant=RUNTIME/'variants'/('saturation_bounds_heatflux_thermal_units' if thermal_units else 'saturation_bounds_heatflux' if heatflux else 'saturation_bounds');variant.mkdir(parents=True,exist_ok=True)
    patched=variant/'Saturation.cpp';patched.write_text(corrected)
    patch=variant/'saturation_bounds.patch';patch.write_text(''.join(difflib.unified_diff(text.splitlines(True),corrected.splitlines(True),fromfile='original/Saturation.cpp',tofile='corrected/Saturation.cpp')))
    cwd=RUNTIME/'biogears-build/projects/biogears/libBiogears';obj=variant/'Saturation.cpp.o'
    subprocess.run(['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(patched),'-o',str(obj)],cwd=cwd,check=True)
    extra_manifest={}
    if heatflux:
        environment_original=source/'projects/biogears/libBiogears/src/engine/Systems/Environment.cpp'
        environment_text=environment_original.read_text();bad_flux='GetEvaporativeHeatLoss().SetValue(dTotalHeatLoss_W, PowerUnit::W);'
        if environment_text.count(bad_flux)!=1:raise RuntimeError('Expected exact heat telemetry defect')
        environment_corrected=environment_text.replace(bad_flux,'GetEvaporativeHeatLoss().SetValue(eHeatLoss_W, PowerUnit::W);')
        environment_telemetry=environment_corrected
        if thermal_units:
            for coefficient in ('Radiative','Convective'):
                wrong=f'dSurfaceArea_m2 / d{coefficient}HeatTransferCoefficient_WPerM2_K'
                if environment_corrected.count(wrong)!=1:raise RuntimeError('Expected exact resistance unit defect')
                environment_corrected=environment_corrected.replace(wrong,f'1.0 / (dSurfaceArea_m2 * d{coefficient}HeatTransferCoefficient_WPerM2_K)')
            thermal_patch=variant/'thermal_units.patch';thermal_patch.write_text(''.join(difflib.unified_diff(environment_telemetry.splitlines(True),environment_corrected.splitlines(True),fromfile='telemetry_corrected/Environment.cpp',tofile='thermal_units_corrected/Environment.cpp')))
        environment_file=variant/'Environment.cpp';environment_file.write_text(environment_corrected)
        environment_patch=variant/'evaporation_telemetry.patch';environment_patch.write_text(''.join(difflib.unified_diff(environment_text.splitlines(True),environment_telemetry.splitlines(True),fromfile='original/Environment.cpp',tofile='corrected/Environment.cpp')))
        subprocess.run(['c++','-DBIOGEARS_COMMON_BUILD_STATIC','-DBIOGEARS_THROW_NAN_EXCEPTIONS','-DBIOGEARS_THROW_READONLY_EXCEPTIONS','-Dbiogears_EXPORTS','@CMakeFiles/libbiogears.dir/includes_CXX.rsp','-std=gnu++20','-fPIC','-O3','-DNDEBUG','-c',str(environment_file),'-o',str(variant/'Environment.cpp.o')],cwd=cwd,check=True)
        sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
        extra_manifest={'environment_original_sha256':sha(environment_original),'environment_patched_sha256':sha(environment_file),'environment_patch_sha256':sha(environment_patch),'environment_scope':'Correct evaporation output accumulator; no consumer in physiological equations.'}
        if thermal_units:extra_manifest.update({'thermal_patch_sha256':sha(thermal_patch),'thermal_scope':'Correct convective/radiative resistance dimensions: A/h to 1/(A*h). Original empirical clothing multipliers retained. Physics variant, not clinically calibrated.'})
    objects=shlex.split((cwd/'CMakeFiles/libbiogears.dir/objects1.rsp').read_text());objects=[str(obj) if o.endswith('/Saturation.cpp.o') else o for o in objects]
    if heatflux:objects=[str(variant/'Environment.cpp.o') if o.endswith('/src/engine/Systems/Environment.cpp.o') else o for o in objects]
    response=variant/'objects.rsp';response.write_text(' '.join(shlex.quote(o) for o in objects))
    command=shlex.split((cwd/'CMakeFiles/libbiogears.dir/link.txt').read_text());lib=variant/'libbiogears.so.8.0.0';command[command.index('-o')+1]=str(lib)
    command=[('@'+str(response)) if c=='@CMakeFiles/libbiogears.dir/objects1.rsp' else c for c in command]
    subprocess.run(command,cwd=cwd,check=True)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest={'variant':variant.name,**extra_manifest,'source_revision':revision,'original_source_sha256':sha(original),'patched_source_sha256':sha(patched),'patch_sha256':sha(patch),'library_sha256':sha(lib),'scope':'Remove erroneous x(3) read/write on a three-element vector. Existing x(0..2) handling and all physiological equations retained.','reproduction':'StandardFemale initialization aborts upstream; debug Eigen assertion in Saturation.cpp:536 confirms index3 out of bounds.'}
    (variant/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--heatflux',action='store_true');parser.add_argument('--thermal-units',action='store_true');args=parser.parse_args();build(args.heatflux,args.thermal_units)
