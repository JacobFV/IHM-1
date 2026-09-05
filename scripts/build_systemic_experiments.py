"""Run continuing native physiology and retain all actions and source receipts."""
import argparse
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from ihm.assembly.systemic import SystemicConfig, run_systemic, verify_contrasts

p=argparse.ArgumentParser()
p.add_argument('--protocols', nargs='+', default=['hydration','meal'])
p.add_argument('--seconds', type=float, default=21600)
p.add_argument('--sample-interval', type=float, default=30)
p.add_argument('--variant', default='whole_body_integrity_renal')
p.add_argument('--output', required=True)
p.add_argument('--workers', type=int, default=2)
args=p.parse_args()
if not 1<=args.workers<=6:p.error('--workers must be 1..6')
if len(set(args.protocols))!=len(args.protocols):p.error('duplicate protocols')
root=Path(__file__).resolve().parents[1]
out=Path(args.output).resolve()
def run(protocol):
    print('Running '+protocol, flush=True)
    result=run_systemic(root,out/protocol,SystemicConfig(protocol=protocol,
        seconds=args.seconds,sample_interval_s=args.sample_interval,engine_variant=args.variant))
    print(json.dumps({'protocol':protocol,'checks':result['checks']}),flush=True)
    return protocol,result
results={}
with ThreadPoolExecutor(max_workers=args.workers) as pool:
    for future in as_completed([pool.submit(run,protocol) for protocol in args.protocols]):
        protocol,result=future.result()
        results[protocol]=result
report=verify_contrasts(results)
(out/'contrasts.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
raise SystemExit(0 if report['passed'] else 1)
