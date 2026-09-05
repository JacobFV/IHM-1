"""Run continuing native physiology and retain all actions and source receipts."""
import argparse
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from ihm.assembly.systemic import SystemicConfig, run_systemic, verify_contrasts
from ihm.native import _sha,RUNTIME

p=argparse.ArgumentParser()
p.add_argument('--protocols', nargs='+', default=['hydration','meal'])
p.add_argument('--seconds', type=float, default=21600)
p.add_argument('--sample-interval', type=float, default=30)
p.add_argument('--variant', default='whole_body_integrity_depletion')
p.add_argument('--state',help='Explicit retained initial state shared by every condition')
p.add_argument('--output', required=True)
p.add_argument('--workers', type=int, default=2)
args=p.parse_args()
if not 1<=args.workers<=6:p.error('--workers must be 1..6')
if len(set(args.protocols))!=len(args.protocols):p.error('duplicate protocols')
root=Path(__file__).resolve().parents[1]
shared=SystemicConfig(seconds=args.seconds,sample_interval_s=args.sample_interval,engine_variant=args.variant,
                      **({'state_path':args.state} if args.state else {}))
state_path=Path(shared.state_path).resolve()
if not state_path.is_relative_to(root):p.error('Retain the initial state inside this workspace')
out=Path(args.output).resolve()
if not out.is_relative_to(root):p.error('Experiment groups must be retained inside this workspace')
if out.exists() and any(out.iterdir()):
    p.error('Choose a fresh experiment group directory')
out.mkdir(parents=True,exist_ok=True)
# Queued workers must not attribute an already-imported numerical module to
# newer on-disk source bytes while an interactive development session edits it.
source_identity={name:_sha(root/name) for name in ('ihm/assembly/systemic.py','ihm/native/session.py',
                                                'ihm/assembly/systemic_evidence.py','ihm/assembly/native_environment_evidence.py',
                                                'ihm/native/__init__.py','scripts/native_body_ports.h',
                                                'scripts/native_biogears_stream.cpp')}
library=(RUNTIME/'biogears-build/outputs/Release/lib' if args.variant=='upstream' else RUNTIME/'variants'/args.variant)/'libbiogears.so.8.0.0'
for path in (state_path,library,RUNTIME/'native_biogears_stream'):
    source_identity[str(path.relative_to(root))]=_sha(path)
(out/'group-configuration.json').write_text(json.dumps(dict(options=vars(args),
    source_identity=source_identity),indent=2)+'\n')
def run(protocol):
    if any(_sha(root/name)!=expected for name,expected in source_identity.items()):
        raise RuntimeError('Executing source changed before queued protocol start; restart in a fresh group')
    print('Running '+protocol, flush=True)
    result=run_systemic(root,out/protocol,SystemicConfig(protocol=protocol,
        seconds=args.seconds,sample_interval_s=args.sample_interval,engine_variant=args.variant,state_path=state_path))
    print(json.dumps({'protocol':protocol,'checks':result['checks']}),flush=True)
    return protocol,result
results={};failures={}
with ThreadPoolExecutor(max_workers=args.workers) as pool:
    futures={pool.submit(run,protocol):protocol for protocol in args.protocols}
    for future in as_completed(futures):
        protocol=futures[future]
        try:
            _,result=future.result()
            results[protocol]=result
        except Exception as error:
            failures[protocol]=dict(error_type=type(error).__name__,error=str(error),
                                    retained_directory=str((out/protocol).relative_to(root)) if (out/protocol).is_dir() else None)
            (out/'failures.json').write_text(json.dumps(failures,indent=2)+'\n')
            print(json.dumps({'protocol':protocol,'failed':failures[protocol]}),flush=True)
try:
    report=verify_contrasts(results)
except Exception as error:
    report=dict(passed=False,checks={},causal_contrast_count=0,
                acceptance_error=dict(error_type=type(error).__name__,error=str(error)))
report['failed_protocols']=failures
report['requested_protocols']=args.protocols
report['passed']=report['passed'] and not failures
report['inputs']={protocol:dict(path=str((out/protocol/'systemic.json').relative_to(root)),
    sha256=_sha(out/protocol/'systemic.json')) for protocol in results}
(out/'contrasts.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
raise SystemExit(0 if report['passed'] else 1)
