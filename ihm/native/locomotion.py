"""Native articulated muscle/contact plant; source subject, not canonical twin."""
from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib, json, math, os, subprocess

BASE = Path(__file__).resolve().parents[2]
RUNTIME = BASE / 'data/runtime/opensim'
SOURCE = BASE / 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking'

@dataclass(frozen=True)
class LocomotionConfig:
    seconds: float = .1
    sample_hz: float = 200
    excitation: float = .03
    reflex_gain: float = 0
    pelvis_speed_delta_m_s: float = 0
    stiffness_scale: float = 1
    friction_scale: float = 1
    accuracy: float = 1e-6
    maximum_step_s: float = .001

    def __post_init__(self):
        bounds = {'seconds':(.005,10), 'sample_hz':(20,1000), 'excitation':(.01,.5),
                  'reflex_gain':(0,2), 'pelvis_speed_delta_m_s':(-.5,.5),
                  'stiffness_scale':(.1,10), 'friction_scale':(0,2),
                  'accuracy':(1e-10,1e-4), 'maximum_step_s':(1e-5,.005)}
        for name,(low,high) in bounds.items():
            value=getattr(self,name)
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not low<=value<=high:
                raise ValueError(f'Invalid finite {name}')
        if abs(self.seconds*self.sample_hz-round(self.seconds*self.sample_hz))>1e-9:
            raise ValueError('Horizon must align to output sampling clock')

def _sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run_locomotion(config: LocomotionConfig, output_dir):
    if not isinstance(config,LocomotionConfig): raise ValueError('Expected LocomotionConfig')
    exe=RUNTIME/'native_opensim_dynamics'
    manifest=json.loads((RUNTIME/'dynamics_build_manifest.json').read_text())
    for relative,digest in manifest['files'].items():
        if _sha(BASE/relative)!=digest: raise ValueError(f'Stale native dynamics build: {relative}')
    out=Path(output_dir).resolve()
    if out.exists(): raise ValueError('Native dynamics requires a fresh output directory')
    out.mkdir(parents=True)
    libdirs=[RUNTIME/'install/opensim/lib',RUNTIME/'install/simbody/lib']
    libdirs += [RUNTIME/'sysroot/usr/lib/aarch64-linux-gnu'/x for x in ('lapack','blas','')]
    env=dict(os.environ,LD_LIBRARY_PATH=':'.join(map(str,libdirs)),OPENBLAS_NUM_THREADS='1')
    command=[str(exe),str(SOURCE),str(out),*[str(v) for v in asdict(config).values()]]
    execution={'config':asdict(config),'command':command,'build':manifest,
               'input_sha256':{p.name:_sha(p) for p in SOURCE.glob('*') if p.is_file()}}
    (out/'execution.json').write_text(json.dumps(execution,indent=2)+'\n')
    with (out/'engine.log').open('w') as log:
        result=subprocess.run(command,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
    execution['returncode']=result.returncode
    (out/'execution.json').write_text(json.dumps(execution,indent=2)+'\n')
    if result.returncode: raise RuntimeError(f'Native locomotion failed: {out}/engine.log')
    def reject(value): raise ValueError(f'Nonfinite native output: {value}')
    data={'model':json.loads((out/'model.json').read_text(),parse_constant=reject),
          'termination':json.loads((out/'termination.json').read_text(),parse_constant=reject),
          'frames':[json.loads(line,parse_constant=reject) for line in (out/'frames.jsonl').read_text().splitlines()],
          'execution':execution,'experimental_validation':False}
    if data['termination']['completed'] and len(data['frames'])!=round(config.seconds*config.sample_hz)+1:
        raise ValueError('Incomplete native dynamics trajectory')
    if not data['frames'] or abs(data['frames'][-1]['time_s']-data['termination']['final_time_s'])>1e-10:
        raise ValueError('Inconsistent native dynamics termination receipt')
    return data
