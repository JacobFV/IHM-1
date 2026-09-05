"""Build the pinned open predictive-motion engine against held OpenSim4."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import shutil
import difflib
import time

ROOT=Path(__file__).resolve().parents[1]
raw=ROOT/'data/raw/mechanics/scone-core'
source=raw/'source'
runtime=ROOT/'data/runtime/scone'
runtime.mkdir(parents=True,exist_ok=True)
receipt=json.loads((raw/'acquisition.json').read_text())
if subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()!=receipt['revision']:
    raise RuntimeError('Unexpected SCONE source')
for name,info in receipt['files'].items():
    if hashlib.sha256((source/name).read_bytes()).hexdigest()!=info['sha256']:
        raise RuntimeError('SCONE source changed: '+name)
subprocess.run(['git','-C',str(source),'submodule','update','--init','--depth','1','submodules/xo','submodules/spot'],check=True)
submodules={}
for name in ('xo','spot'):
    path=source/'submodules'/name
    revision=subprocess.check_output(['git','-C',str(path),'rev-parse','HEAD'],text=True).strip()
    expected=subprocess.check_output(['git','-C',str(source),'ls-tree','HEAD',f'submodules/{name}'],text=True).split()[2]
    if revision!=expected:raise RuntimeError('SCONE submodule revision differs')
    if subprocess.check_output(['git','-C',str(path),'status','--porcelain'],text=True).strip():raise RuntimeError('Modified SCONE submodule '+name)
    tracked=subprocess.check_output(['git','-C',str(path),'ls-files','-z']).decode().split('\0')
    submodules[name]={'revision':revision,'files':{n:hashlib.sha256((path/n).read_bytes()).hexdigest() for n in tracked if n and (path/n).is_file()}}
(raw/'submodule_receipts.json').write_text(json.dumps(submodules,indent=2)+'\n')
# Preserve the acquired donor. Linux/OpenSim4 API compatibility is an explicit
# build-source variant, with exact patches and hashes retained beside the build.
variant=runtime/'source-opensim4'
# Refresh all donor bytes so stale edits cannot become an undocumented variant.
if variant.exists():
    unexpected=[f for f in variant.rglob('*') if f.is_file() and not (source/f.relative_to(variant)).is_file()]
    if unexpected:raise RuntimeError('Unexpected files in build source variant: '+str(unexpected[:3]))
shutil.copytree(source,variant,ignore=shutil.ignore_patterns('.git'),dirs_exist_ok=True)
patches={}
def patch(relative, transform):
    original=(source/relative).read_text()
    modified=transform(original)
    (variant/relative).write_text(modified)
    patches[relative]={'original_sha256':hashlib.sha256(original.encode()).hexdigest(),
                      'patched_sha256':hashlib.sha256(modified.encode()).hexdigest(),
                      'diff':''.join(difflib.unified_diff(original.splitlines(True),modified.splitlines(True),fromfile='original/'+relative,tofile='variant/'+relative))}
for path in (source/'src/sconelib/sconeopensim4').glob('*.cpp'):
    relative=str(path.relative_to(source))
    def compatibility(text, name=path.name):
        text='#include <OpenSim/OpenSim.h>\n'+text
        if name=='ModelOpenSim4.cpp':
            text=text.replace('cg_osim->getBody().getName()', 'cg_osim->getFrame().findBaseFrame().getName()')
            text=text.replace('cg_osim->getLocation()', 'cg_osim->get_location()')
            text=text.replace('cg_osim->getOrientation()', 'cg_osim->get_orientation()')
            text=text.replace('auto& name = cg_osim->getName();', 'auto& name = cg_osim->getName();\n\t\t\tif (&cg_osim->getFrame() != &cg_osim->getFrame().findBaseFrame()) SCONE_THROW("Offset contact frames require explicit transform mapping");')
            # This source requires fixed_control_step_size and uses its own
            # TimeStepper. Model::RequestTermination is checked in that loop.
            text=text.replace('m_pOsimManager->halt(); // needed when using an OpenSim::Manager', '// Fixed-step TimeStepper loop observes Model::RequestTermination; obsolete Manager path is disabled.')
        return text
    patch(relative,compatibility)
patch('submodules/xo/xo/system/version.h',lambda text:text.replace('#include "xo/string/string_type.h"','#include "xo/string/string_type.h"\n#include "xo/string/string_cast.h"').replace('string str() const { }','string str() const { return xo::to_str(*this); }'))
patch('src/sconelib/scone/core/system_tools.cpp',lambda text:'#include <cstdlib>\n'+text.replace(
    'return xo::get_config_dir() / "SCONE";',
    'if (const char* directory = std::getenv("IHM_SCONE_CONFIG_DIR")) return path(directory);\n\t\treturn xo::get_config_dir() / "SCONE";'))
for name in ('ValueArg','MultiArg'):
    patch(f'contrib/tclap-1.2.1/include/tclap/{name}.h',lambda text,n=name:text.replace(f'{n}<T>(const {n}<T>&',f'{n}(const {n}<T>&'))
(runtime/'compatibility_patches.json').write_text(json.dumps(patches,indent=2)+'\n')
opensim=ROOT/'data/runtime/opensim'
command=['cmake','-S',str(variant),'-B',str(runtime/'build-opensim4'),'-DCMAKE_BUILD_TYPE=Release',
         '-DSCONE_OPENSIM_4=ON','-DSCONE_OPENSIM_3=OFF','-DSCONE_LUA=OFF','-DSCONE_HYFYDY=OFF',
         '-DSCONE_PYTHON=OFF','-DXO_TEST_ENABLED=OFF','-DCMAKE_POLICY_VERSION_MINIMUM=3.5',
         f'-DCMAKE_PREFIX_PATH={opensim}/install/opensim;{opensim}/install/simbody',
         f'-DCMAKE_EXE_LINKER_FLAGS=-Wl,-rpath,{opensim}/sysroot/usr/lib/aarch64-linux-gnu',
         f'-DCMAKE_INSTALL_PREFIX={runtime}/install']
commands=[command,['cmake','--build',str(runtime/'build-opensim4'),'--target','sconecmd','-j','6']]
if (runtime/'build-opensim4.log').exists():
    shutil.copyfile(runtime/'build-opensim4.log',runtime/f'build-attempt-{time.time_ns()}.log')
with (runtime/'build-opensim4.log').open('w') as log:
    for command in commands:
        run=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        if run.returncode:raise RuntimeError(f'Build failed; inspect {runtime}/build-opensim4.log')
binary=runtime/'build-opensim4/bin/sconecmd'
(runtime/'build-opensim4/.sconeroot').write_text(str(variant))
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
link_env={**os.environ,'LD_LIBRARY_PATH':':'.join(map(str,[runtime/'build-opensim4/bin',opensim/'install/opensim/lib',opensim/'install/simbody/lib',opensim/'sysroot/usr/lib/aarch64-linux-gnu']))}
linkage=subprocess.check_output(['ldd',str(binary)],env=link_env,text=True)
dependencies={}
for line in linkage.splitlines():
    if 'not found' in line:raise RuntimeError('Missing SCONE dependency')
    if '=>' in line:
        path=Path(line.split('=>',1)[1].split(' (',1)[0].strip())
        if path.is_file():dependencies[str(path.resolve())]=sha(path)
record={'variant_file_sha256':{str(f.relative_to(variant)):sha(f) for f in sorted(variant.rglob('*')) if f.is_file()},
        'dependency_sha256':dependencies,'source_revision':receipt['revision'],'source_receipt_sha256':hashlib.sha256((raw/'acquisition.json').read_bytes()).hexdigest(),
        'submodule_receipt_sha256':hashlib.sha256((raw/'submodule_receipts.json').read_bytes()).hexdigest(),
        'commands':commands,'executable_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
        'compatibility_patches_sha256':hashlib.sha256((runtime/'compatibility_patches.json').read_bytes()).hexdigest(),
        'opensim_build_manifest_sha256':hashlib.sha256((opensim/'build_manifest.json').read_bytes()).hexdigest()}
(runtime/'build_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
print(binary)
