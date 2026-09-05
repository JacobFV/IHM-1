"""Validated, source-preserving local BioGears predictor. No replacement physiology."""
from dataclasses import dataclass, asdict
from pathlib import Path
import csv, hashlib, json, math, os, re, subprocess, time
import xml.etree.ElementTree as ET
from .bindings import variable_bindings
BASE=Path(__file__).resolve().parents[2]
RUNTIME=BASE/'data/runtime/physiology'
SOURCE=BASE/'data/raw/physiology/biogears'
SOURCE_REVISION='3f16a5fa1dade9c511b88d923606fa51cc35e95d'
def number(x,lo,hi,name):
    if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or not lo<=x<=hi: raise ValueError(f'{name} must be finite in [{lo}, {hi}]')
def available_patients():
    return sorted(p.stem for p in (RUNTIME/'biogears-build/runtime/patients').glob('*.xml'))
@dataclass(frozen=True)
class Intervention:
    time_s: float
    kind: str
    value: float
    def __post_init__(self):
        number(self.time_s,0,3600,'time_s')
        if self.kind not in ('exercise','hemorrhage','saline'):raise ValueError('Unsupported intervention')
        number(self.value,0,.5 if self.kind=='exercise' else 100,'value')
        if abs(self.time_s*50-round(self.time_s*50))>1e-7:raise ValueError('Intervention must align to 0.02 s solver step')
@dataclass(frozen=True)
class NativeConfig:
    seconds: float=60
    patient: str='StandardMale'
    state_path: str|None=None
    interventions: tuple=()
    sample_hz: int=50
    engine_variant: str="upstream"
    ambient_temperature_c: float|None=None
    clothing_clo: float|None=None
    def __post_init__(self):
        number(self.seconds,.02,3600,'seconds')
        if self.engine_variant not in ('upstream','saturation_bounds','saturation_bounds_heatflux','saturation_bounds_heatflux_thermal_units'):raise ValueError('Unknown engine variant')
        if self.ambient_temperature_c is not None:number(self.ambient_temperature_c,10,35,'ambient_temperature_c')
        if self.clothing_clo is not None:number(self.clothing_clo,0,3,'clothing_clo')
        if abs(self.seconds*50-round(self.seconds*50))>1e-7:raise ValueError('Duration must align to solver step')
        if isinstance(self.sample_hz,bool) or self.sample_hz not in (1,2,5,10,25,50):raise ValueError('sample_hz must divide 50')
        if self.seconds*self.sample_hz<1-1e-9:raise ValueError('Duration must include at least one sample')
        if not isinstance(self.patient,str) or not re.fullmatch(r'[A-Za-z0-9_]+',self.patient) or self.patient not in available_patients():raise ValueError('Unknown upstream patient')
        if self.state_path is not None and (not isinstance(self.state_path,(str,Path)) or not Path(self.state_path).is_file()):raise ValueError('State file unavailable')
        if not isinstance(self.interventions,(list,tuple)) or len(self.interventions)>100:raise ValueError('Invalid intervention list')
        last=-1;seen=set()
        for event in self.interventions:
            if not isinstance(event,Intervention):raise ValueError('Expected Intervention')
            if event.time_s<last or event.time_s>=self.seconds:raise ValueError('Timeline must be ordered within duration')
            if (event.time_s,event.kind) in seen:raise ValueError('Duplicate simultaneous action')
            seen.add((event.time_s,event.kind));last=event.time_s
        object.__setattr__(self,'interventions',tuple(self.interventions))
    @classmethod
    def from_dict(cls,data):
        if not isinstance(data,dict) or set(data)-{'seconds','patient','state_path','interventions','sample_hz','ambient_temperature_c','clothing_clo','engine_variant'}:raise ValueError('Unknown native configuration fields')
        try:
            data=dict(data)
            if 'interventions' in data:data['interventions']=tuple(Intervention(**e) for e in data['interventions'])
            return cls(**data)
        except (TypeError,KeyError) as e:raise ValueError('Invalid native configuration') from e

def load_trajectory(path):
    path=Path(path)
    if path.is_dir():path=path/'native_multisystem.csv'
    with path.open() as f:
        reader=csv.reader(f)
        try:header=next(reader);rows=[[float(v) for v in r] for r in reader if r]
        except (StopIteration,ValueError) as e:raise ValueError('Malformed native CSV') from e
    if not header or header[0] not in ('Time(s)', '#Time(s)') or len(set(header))!=len(header) or not rows or any(len(r)!=len(header) for r in rows):raise ValueError('Empty, duplicate or ragged native CSV')
    if not all(math.isfinite(v) for r in rows for v in r):raise ValueError('Nonfinite native trajectory')
    times=[r[0] for r in rows]
    if any(b<=a for a,b in zip(times,times[1:])):raise ValueError('Non-increasing time')
    return {'columns':header,'time_s':times,'values':{h:[r[j] for r in rows] for j,h in enumerate(header) if j}}

def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def state_nutrition_context(path):
    """Read original serialized stomach inputs without inferring a fasted state."""
    path=Path(path)
    if not path.is_file():return None
    tree=ET.parse(path)
    stomach=next((node for node in tree.iter() if node.tag.split('}')[-1]=='StomachContents'),None)
    if stomach is None:return None
    return {node.tag.split('}')[-1]:{'value':float(node.attrib['value']),'unit':node.attrib.get('unit')} for node in stomach if 'value' in node.attrib}

def summarize(output_dir,config=None):
    out=Path(output_dir);trajectory=load_trajectory(out);times=trajectory['time_s'];steps=[b-a for a,b in zip(times,times[1:])]
    log='\n'.join(p.read_text(errors='replace') for p in [out/'runner_stdout.log',out/'native_engine.log'] if p.exists())
    if re.search(r'FATAL|Unknown Data Request|\[ERROR\]|<ERROR>|entered irreversible state',log,re.I):raise ValueError('Engine logged fatal/error/unknown request; inspect preserved logs')
    version=re.search(r'ENGINE_VERSION=(.*)',log)
    if not version:raise ValueError('Missing engine version')
    if config:
        expected=1/config.sample_hz
        if steps and max(abs(s-expected) for s in steps)>1e-5:raise ValueError('Unexpected temporal sampling')
        start=float(re.search(r'STABILIZED_TIME_S=(.*)',log)[1])
        final=float(re.search(r'FINAL_TIME_S=(.*)',log)[1])
        if abs(final-start-config.seconds)>1e-4:raise ValueError('Truncated native execution')
        count=round(config.seconds*50)//round(50/config.sample_hz)
        if len(times)!=count or abs(times[0]-start-expected)>1e-4 or abs(times[-1]-start-count*expected)>1e-4:raise ValueError('Truncated or misaligned native trajectory')
    stats={h:{'initial':v[0],'final':v[-1],'minimum':min(v),'maximum':max(v)} for h,v in trajectory['values'].items()}
    summary={'engine':'BioGears','source_revision':subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip(),'compiled_engine_version':version[1],'rows':len(times),'columns':len(trajectory['columns']),'output_columns':trajectory['columns'],'time_start_s':times[0],'time_end_s':times[-1],'sample_interval_min_s':min(steps) if steps else None,'sample_interval_max_s':max(steps) if steps else None,'all_finite':True,'summary':stats,'csv_sha256':_sha(out/'native_multisystem.csv'),'executable_sha256':_sha(RUNTIME/'native_biogears_rest'),'adapter_sha256':_sha(BASE/'scripts/native_biogears_rest.cpp'),'experimental_validation':False,'uncertainty':{'parameter_covariance':None,'status':'unquantified; deterministic source simulation, not calibrated human evidence'},'limitations':['No explicit posture action found in upstream source.','Original upstream equations; no independent clinical validation.','A stopped intervention initiates recovery but does not guarantee return to baseline.','Short upstream stabilization is not evidence of hourly thermal or fluid steady state.','Environment radiation uses an upstream standing-area factor; a supine hydrostatic adjunct does not change this thermal assumption.','EvaporativeHeatLoss repeats convective heat loss due to an upstream telemetry assignment defect; excluded from energy accounting.'],'configuration':asdict(config) if config else None,'states':{p.name:_sha(p) for p in (out/'states').glob('*.xml')}}
    compartment_file=out/'body_compartments.csv'
    if compartment_file.exists():
        with compartment_file.open() as f:
            reader=csv.reader(f);header=next(reader,[]);rows=[[float(v) for v in r] for r in reader]
        if not header or header[0]!='Time(s)' or len(set(header))!=len(header) or any(len(r)!=len(header) for r in rows) or not all(math.isfinite(v) for r in rows for v in r):raise ValueError('Malformed compartment telemetry')
        if len(rows)!=len(times)+1 or any(abs(row[0]-t)>1e-8 for row,t in zip(rows[1:],times)):raise ValueError('Compartment telemetry clock mismatch')
        summary['compartment_telemetry']={'csv_sha256':_sha(compartment_file),'rows':len(rows),'columns':header,'receipt_basis':'validated and hashed during native run summarization'}
    summary['environment_context_lines']=[line for line in log.splitlines() if line.startswith('ENVIRONMENT_')]
    variant=config.engine_variant if config else 'upstream'
    summary['engine_variant']=variant
    summary['variant_manifest']=json.loads((RUNTIME/'variants'/variant/'manifest.json').read_text()) if variant!='upstream' else None
    summary['nutrition_context']={'initial_stomach_contents':state_nutrition_context(out/'states/native_stabilized.xml'),'final_stomach_contents':state_nutrition_context(out/'states/native_final.xml'),'interpretation':'Default stomach water and electrolytes are model inputs; resting does not imply fasting or fixed blood volume.'}
    summary['input_state_sha256']=_sha(config.state_path) if config and config.state_path else None
    summary['patient_sha256']=_sha(RUNTIME/'biogears-build/runtime/patients'/f'{config.patient}.xml') if config else None
    summary['variable_bindings']=variable_bindings(trajectory['columns'][1:],summary['source_revision'])
    if variant in ('saturation_bounds_heatflux','saturation_bounds_heatflux_thermal_units'):
        summary['limitations']=[text for text in summary['limitations'] if not text.startswith('EvaporativeHeatLoss repeats')]
        summary['limitations'].append('Evaporation telemetry corrected in explicit source variant; skin heat loss is an aggregate, not an additional independent heat pathway.')
        binding=summary['variable_bindings'].get('EvaporativeHeatLoss(W)')
        if binding:
            binding['usable_for_energy_accounting']=True
            binding['source_diagnostic_defect']='Corrected by isolated evaporation_telemetry.patch; original physiological equations retained.'
    if variant=='saturation_bounds_heatflux_thermal_units':
        summary['limitations']=[text for text in summary['limitations'] if not text.startswith('Original upstream equations;')]
        summary['limitations'].append('Experimental thermal dimensions correction changes source heat-transfer equations; empirical clothing factors retained, no clinical calibration.')
    if summary['source_revision'] != SOURCE_REVISION:raise ValueError('Unexpected upstream source revision')
    (out/'summary.json').write_text(json.dumps(summary,indent=2,default=str)+'\n');return summary

def run_native(config,output_dir):
    if isinstance(config,dict):config=NativeConfig.from_dict(config)
    if not isinstance(config,NativeConfig):raise ValueError('Expected NativeConfig')
    exe=RUNTIME/'native_biogears_rest'
    if not exe.is_file():raise RuntimeError('Native backend unavailable; run scripts/build_native_biogears.py')
    out=Path(output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'execution.json').exists():raise ValueError('Output already contains an execution; choose a fresh directory')
    runtime=RUNTIME/'biogears-build/runtime'
    for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:
        p=out/name
        if not p.exists():p.symlink_to(runtime/name,target_is_directory=(runtime/name).is_dir())
    (out/'states').mkdir(exist_ok=True)
    (out/'configuration.json').write_text(json.dumps(asdict(config),indent=2,default=str)+'\n')
    (out/'timeline.tsv').write_text(''.join(f'{e.time_s} {e.kind} {e.value}\n' for e in config.interventions))
    state=str(Path(config.state_path).resolve()) if config.state_path else '-'
    command=[str(exe),str(config.seconds),config.patient,state,str(config.sample_hz),'timeline.tsv',str(config.ambient_temperature_c) if config.ambient_temperature_c is not None else '-',str(config.clothing_clo) if config.clothing_clo is not None else '-']
    environment=os.environ.copy()
    if config.engine_variant!='upstream':
        variant=RUNTIME/'variants'/config.engine_variant
        manifest=json.loads((variant/'manifest.json').read_text())
        if _sha(variant/'libbiogears.so.8.0.0')!=manifest['library_sha256']:raise ValueError('Variant library integrity mismatch')
        environment['LD_LIBRARY_PATH']=str(variant)+os.pathsep+environment.get('LD_LIBRARY_PATH','')
    start=time.monotonic()
    try:
        with (out/'runner_stdout.log').open('w') as f:result=subprocess.run(command,cwd=out,stdout=f,stderr=subprocess.STDOUT,timeout=1800,env=environment)
    except subprocess.TimeoutExpired:
        (out/'execution.json').write_text(json.dumps({'exit_code':None,'timed_out':True,'wall_seconds':time.monotonic()-start,'command':command},indent=2)+'\n')
        raise RuntimeError(f'Native engine timed out; see {out}/runner_stdout.log') from None
    (out/'execution.json').write_text(json.dumps({'exit_code':result.returncode,'wall_seconds':time.monotonic()-start,'command':command},indent=2)+'\n')
    if result.returncode:raise RuntimeError(f'Native engine failed ({result.returncode}); see {out}/runner_stdout.log')
    return summarize(out,config)
