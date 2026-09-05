"""Source-engine OpenSim geometry/force export, independent of BioGears physiology."""
from dataclasses import dataclass,asdict
from pathlib import Path
import math,re,json,subprocess,hashlib,os
BASE=Path(__file__).resolve().parents[2];RUNTIME=BASE/'data/runtime/opensim'
@dataclass(frozen=True)
class OpenSimConfig:
    coordinate:str='knee_angle_r'
    delta_rad:float=0
    activation:float=.05
    engine_variant:str='upstream'
    def __post_init__(self):
        if self.engine_variant!='upstream' and not re.fullmatch(r'wrap_[0-9]+_[0-9.e+-]+(?:_cache)?',self.engine_variant):raise ValueError('Invalid engine variant')
        if not isinstance(self.coordinate,str) or not re.fullmatch('[A-Za-z0-9_]+',self.coordinate):raise ValueError('Invalid coordinate')
        for x,lo,hi in [(self.delta_rad,-.15,.15),(self.activation,.01,.5)]:
            if isinstance(x,bool) or not isinstance(x,(float,int)) or not math.isfinite(x) or not lo<=x<=hi:raise ValueError('Unsupported finite mechanics input')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run_opensim(config,output_dir):
    if not isinstance(config,OpenSimConfig):raise ValueError('Expected OpenSimConfig')
    out=Path(output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'execution.json').exists():raise ValueError('Choose fresh output directory')
    exe=RUNTIME/'native_opensim_export'
    if not exe.is_file():raise RuntimeError('Build native OpenSim first')
    source=BASE/'data/raw/anatomy/opensim-models/source/Models/Rajagopal/Rajagopal2016.osim'
    env=dict(os.environ);libdirs=[RUNTIME/'install/opensim/lib',RUNTIME/'install/simbody/lib',RUNTIME/'sysroot/usr/lib/aarch64-linux-gnu/lapack',RUNTIME/'sysroot/usr/lib/aarch64-linux-gnu/blas',RUNTIME/'sysroot/usr/lib/aarch64-linux-gnu']
    variant=RUNTIME/'variants'/config.engine_variant
    simulation_library=RUNTIME/'install/opensim/lib/libosimSimulation.so'
    if config.engine_variant!='upstream':
        if not (variant/'manifest.json').is_file():raise ValueError('Experimental variant unavailable')
        simulation_library=variant/'libosimSimulation.so';libdirs.insert(0,variant)
    env['LD_LIBRARY_PATH']=':'.join(str(p) for p in libdirs)
    command=[str(exe),str(source),str(out/'mechanics.json'),config.coordinate,str(config.delta_rad),str(config.activation)]
    with (out/'engine.log').open('w') as f:result=subprocess.run(command,env=env,cwd=out,stdout=f,stderr=subprocess.STDOUT,timeout=600)
    execution={'returncode':result.returncode,'config':asdict(config),'command':command,'model_sha256':sha(source),'executable_sha256':sha(exe),'adapter_sha256':sha(BASE/'scripts/native_opensim_export.cpp'),'engine_source_revision':subprocess.check_output(['git','-C',str(BASE/'data/raw/mechanics/opensim-core'),'rev-parse','HEAD'],text=True).strip(),'simbody_source_revision':subprocess.check_output(['git','-C',str(BASE/'data/raw/mechanics/simbody'),'rev-parse','HEAD'],text=True).strip()}
    execution['simulation_library_sha256']=sha(simulation_library)
    execution['engine_variant']=config.engine_variant
    (out/'execution.json').write_text(json.dumps(execution,indent=2)+'\n')
    if result.returncode:raise RuntimeError(f'OpenSim failed; see {out}/engine.log')
    data=json.loads((out/'mechanics.json').read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Nonfinite native mechanics')))
    data['execution']=execution;data['experimental_validation']=False
    (out/'summary.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    return data
