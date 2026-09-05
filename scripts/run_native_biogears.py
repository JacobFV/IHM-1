#!/usr/bin/env python3
"""Run the native upstream whole-body engine and summarize unmodified tracker CSV."""
from pathlib import Path
import argparse,csv,hashlib,json,math,re,subprocess,time
base=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=60);parser.add_argument('--summarize-only',action='store_true');parser.add_argument('--output-dir',type=Path,default=base/'data/derived/physiology/biogears_native_run');args=parser.parse_args()
root=base/'data/runtime/physiology';source=base/'data/raw/physiology/biogears';out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
exe=root/'native_biogears_rest';runtime=root/'biogears-build/runtime'
for name in ['patients','substances','environments','nutrition','config','ecg','xsd','UCEDefs.conf','BioGearsConfiguration.xml']:
    p=out/name
    if not p.exists():p.symlink_to(runtime/name,target_is_directory=(runtime/name).is_dir())
if not args.summarize_only:
    start=time.monotonic()
    with (out/'runner_stdout.log').open('w') as f:
        result=subprocess.run([str(exe),str(args.seconds)],cwd=out,stdout=f,stderr=subprocess.STDOUT)
    (out/'execution.json').write_text(json.dumps(dict(exit_code=result.returncode,wall_seconds=time.monotonic()-start,command=[str(exe),str(args.seconds)]),indent=2)+'\n')
    if result.returncode:raise SystemExit(f'Native engine failed ({result.returncode}); see {out}/runner_stdout.log')
path=out/'native_multisystem.csv'
with path.open() as f:
    reader=csv.reader(f);header=next(reader);rows=[[float(v) for v in r] for r in reader if r]
if not rows or any(len(r)!=len(header) for r in rows):raise ValueError('Empty or ragged native CSV')
if not all(math.isfinite(v) for r in rows for v in r):raise ValueError('Nonfinite native trajectory')
times=[r[0] for r in rows];steps=[b-a for a,b in zip(times,times[1:])]
if not all(s>0 for s in steps):raise ValueError('Non-increasing time')
stats={h:dict(initial=rows[0][j],final=rows[-1][j],minimum=min(r[j] for r in rows),maximum=max(r[j] for r in rows)) for j,h in enumerate(header)}
log=(out/'runner_stdout.log').read_text(errors='replace')
if '<FATAL>' in log or 'Unknown Data Request' in log:raise ValueError('Engine logged fatal error despite its exit status')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
summary=dict(engine='BioGears',source_revision=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip(),
    compiled_engine_version=re.search(r'ENGINE_VERSION=(.*)',log)[1],rows=len(rows),columns=len(header),output_columns=header,
    csv_sha256=sha(path),executable_sha256=sha(exe),adapter_sha256=sha(base/'scripts/native_biogears_rest.cpp'),
    time_start_s=times[0],time_end_s=times[-1],sample_interval_min_s=min(steps),sample_interval_max_s=max(steps),
    all_finite=True,summary=stats,
    stabilization_log_lines=[line for line in log.splitlines() if any(s in line for s in ['stabilized','Stabilized','STABILIZED_TIME','steady state','Steady State','Convergence took','Engine is now Active','Resting Stabilization','Secondary Stabilization'])],
    experimental_validation=False,
    limitations=['Initialization/resting execution only; no independent experimental validation.',
                 'StandardMale upstream baseline, no explicit supine-posture action configured.',
                 'Original integrated engine equations preserved; not attached to IHM reduced runtime.',
                 'Full SHA is the source version authority; shallow-clone build reports fallback 8.0.0-Source.',
                 'No subject-specific parameter fit, sensitivity analysis or physiological bounds certification.'])
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k not in ['summary','stabilization_log_lines']},indent=2))
