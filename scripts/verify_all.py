"""Run numerical/data/API checks and persist their actual outputs and return codes."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import time

root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--native',action='store_true');p.add_argument('--app',action='store_true');args=p.parse_args()
names=['verify.py','verify_evidence.py','verify_skin.py','verify_fluids.py','verify_cardiopulmonary.py','verify_population.py','verify_collected_data.py','verify_spatial.py','verify_coupling.py','verify_calibration.py','verify_temporal.py','verify_native.py','verify_native_circuits.py','verify_circuit_predictor.py','verify_opensim.py','verify_native_opensim.py','verify_reproductive.py','verify_vascular_flow.py','verify_human.py','verify_app.py']
commands=[[sys.executable,'scripts/'+name] for name in names]
if args.native:commands += [[sys.executable,'scripts/verify_native.py','--engine-clock'],[sys.executable,'scripts/verify_native_opensim.py','--corrected']]
if args.app:commands += [['npm','--prefix','app','test'],['npm','--prefix','app','run','build'],['npm','--prefix','app','run','test:browser']]
results=[];out=root/'artifacts/verification';out.mkdir(parents=True,exist_ok=True)
for i,command in enumerate(commands):
 start=time.time();run=subprocess.run(command,cwd=root,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 log=out/f'{i:02d}.log';log.write_text(run.stdout)
 result=dict(command=command,returncode=run.returncode,elapsed_s=time.time()-start,log=str(log.relative_to(root)))
 results.append(result);print(('PASS' if run.returncode==0 else 'FAIL')+' '+' '.join(command),flush=True)
 (out/'report.json').write_text(json.dumps(dict(all_passed=all(x['returncode']==0 for x in results),results=results),indent=2)+'\n')
sys.exit(int(any(x['returncode'] for x in results)))
