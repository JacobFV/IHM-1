"""Evaluate acquired reflex gait parameters in the actual native source engine."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess

p=argparse.ArgumentParser()
p.add_argument('--output',required=True)
p.add_argument('--controller-off',action='store_true')
p.add_argument('--allow-early-termination',action='store_true')
p.add_argument('--seconds',type=float,default=10)
p.add_argument('--accuracy',type=float,default=.002)
args=p.parse_args()
if not math.isfinite(args.seconds) or not 0<args.seconds<=30:p.error('seconds must be in (0,30]')
if not math.isfinite(args.accuracy) or not 1e-7<=args.accuracy<=.01:p.error('accuracy must be 1e-7..0.01')
root=Path(__file__).resolve().parents[1]
runtime=root/'data/runtime/scone'
source=root/'data/raw/mechanics/scone-core/source'
binary=runtime/'build-opensim4/bin/sconecmd'
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
build=json.loads((runtime/'build_manifest.json').read_text())
if sha(binary)!=build['executable_sha256']:raise RuntimeError('Native SCONE executable changed')
for relative, expected in build['variant_file_sha256'].items():
    if sha(runtime/'source-opensim4'/relative)!=expected:raise RuntimeError('Build source changed: '+relative)
for relative,key in [('acquisition.json','source_receipt_sha256'),('submodule_receipts.json','submodule_receipt_sha256')]:
    if sha(source.parent/relative)!=build[key]:raise RuntimeError('Source receipt changed')
acquisition=json.loads((source.parent/'acquisition.json').read_text())
for relative,info in acquisition['files'].items():
    if sha(source/relative)!=info['sha256']:raise RuntimeError('Acquired source changed: '+relative)
submodules=json.loads((source.parent/'submodule_receipts.json').read_text())
for name,info in submodules.items():
    for relative,expected in info['files'].items():
        if sha(source/'submodules'/name/relative)!=expected:raise RuntimeError('Acquired submodule changed')
if sha(runtime/'compatibility_patches.json')!=build['compatibility_patches_sha256']:raise RuntimeError('Compatibility receipt changed')
out=Path(args.output).resolve()
if out.exists() and any(out.iterdir()):raise ValueError('Choose a fresh output directory')
out.mkdir(parents=True,exist_ok=True)
scenario=source/'scenarios/Examples2/Gait - H0914 - OpenSim4.scone'
parameters=source/'scenarios/Examples2/data/ResultH0914Gait10.par'
(out/'data').symlink_to(source/'scenarios/Examples2/data',target_is_directory=True)
text=scenario.read_text().replace('max_duration = 10',f'max_duration = {args.seconds}').replace('integration_accuracy = 0.002',f'integration_accuracy = {args.accuracy}')
if args.controller_off:
    controller=(source/'scenarios/Examples2/data/ControllerGH2010v9.scone').read_text()
    if controller.count('GaitStateController {')!=1:raise RuntimeError('Controller source layout changed')
    controller=controller.replace('GaitStateController {','GaitStateController {\n stop_time = 1',1)
    (out/'controller_off.scone').write_text(controller)
    text=text.replace('<< data/ControllerGH2010v9.scone >>','<< controller_off.scone >>')
(out/'evaluation.scone').write_text(text)
shutil.copyfile(parameters,out/'evaluation.par')
settings=out/'settings';settings.mkdir()
(settings/'scone-settings.zml').write_text('folders {\nresults = "'+str(out/'results')+'"\nscenarios = "'+str(out)+'"\n}\ndata { frequency = 100\nbody = 1\nmuscle_detail = 1\ncontact = 1\npower = 1\n}\noptimizer { max_threads = 1\nevaluator = 0\n}\n')
opensim=root/'data/runtime/opensim'
env={**os.environ,'OPENBLAS_NUM_THREADS':'1','IHM_SCONE_CONFIG_DIR':str(settings),
     'LD_LIBRARY_PATH':':'.join(map(str,[runtime/'build-opensim4/bin',opensim/'install/opensim/lib',opensim/'install/simbody/lib',opensim/'sysroot/usr/lib/aarch64-linux-gnu']))}
command=[str(binary),'-e',str(out/'evaluation.par'),'-r',str(out/'trajectory.sto')]
linkage=subprocess.check_output(['ldd',str(binary)],env=env,text=True)
dependencies={}
for line in linkage.splitlines():
    if 'not found' in line:raise RuntimeError('Missing native SCONE dependency')
    if '=>' in line:
        path=Path(line.split('=>',1)[1].split(' (',1)[0].strip())
        if path.is_file():dependencies[str(path.resolve())]=sha(path)
if dependencies!=build['dependency_sha256']:raise RuntimeError('Native dependency identity changed')
receipt={'controller_stop_time_s':1 if args.controller_off else None,'controller_off':args.controller_off,'source_revision':build['source_revision'],'build_manifest_sha256':sha(runtime/'build_manifest.json'),
         'command':command,'source_scenario_sha256':sha(scenario),'source_parameters_sha256':sha(parameters),
         'executed_scenario_sha256':sha(out/'evaluation.scone'),'dependency_sha256':dependencies,
         'scope':'Source planar nine-DOF, fourteen-muscle reflex gait evaluation; not canonical 3D walking or neural subject calibration',
         'requested_seconds':args.seconds,'integration_accuracy':args.accuracy,
         'executed_file_sha256':{str(f.relative_to(out)):sha(f) for f in out.rglob('*') if f.is_file() and not f.is_symlink()}}
(out/'manifest.json').write_text(json.dumps(receipt,indent=2)+'\n')
with (out/'runner.log').open('w') as log:
    try:result=subprocess.run(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600)
    except subprocess.TimeoutExpired:
        (out/'execution.json').write_text(json.dumps({'exit_code':None,'status':'timeout','requested_seconds':args.seconds})+'\n')
        raise
log=(out/'runner.log').read_text()
matches=re.findall(r'simulation time\s*=\s*([\d.eE+-]+)',log)
report={'exit_code':result.returncode,'native_simulation_time_s':float(matches[-1]) if matches else None,
        'requested_seconds':args.seconds,'outputs':[str(f.relative_to(out)) for f in out.glob('trajectory*')]}
report['completed_requested_duration']=report['native_simulation_time_s'] is not None and abs(report['native_simulation_time_s']-args.seconds)<=.011
report['output_sha256']={str(f.relative_to(out)):sha(f) for f in out.glob('trajectory*') if f.is_file()}
report['input_integrity_after_run']=all(sha(out/name)==expected for name,expected in receipt['executed_file_sha256'].items()) and all(sha(Path(name))==expected for name,expected in dependencies.items()) and sha(binary)==build['executable_sha256']
(out/'execution.json').write_text(json.dumps(report,indent=2)+'\n')
if not report['input_integrity_after_run']:raise RuntimeError('Run inputs changed during evaluation')
print(json.dumps(report,indent=2))
if result.returncode:raise SystemExit(result.returncode)
if not report['outputs'] or report['native_simulation_time_s'] is None:raise RuntimeError('Native run produced no trajectory/time')
if abs(report['native_simulation_time_s']-args.seconds)>.011 and not args.allow_early_termination:
    raise RuntimeError('Native engine terminated before requested duration; retained execution is incomplete')
