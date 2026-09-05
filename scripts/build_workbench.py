"""Rebuild portable views from the acquired corpus in dependency order.

This does not redownload datasets or rerun patient simulations. The optional
--native-adapter flag recompiles the thin adapter against existing libraries.
Use source-specific collectors/builders to acquire or build native engines.
"""
from pathlib import Path
import subprocess
import sys
import json
import argparse

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--native-adapter',action='store_true',help='Recompile the thin native adapter before canonical profile generation')
p.add_argument('--plan',action='store_true',help='Print the ordered command list without executing it')
args=p.parse_args()
root=Path(__file__).resolve().parents[1]
commands=[
 [sys.executable,'scripts/build_spatial_atlas.py'],
 [sys.executable,'scripts/build_opensim_display.py'],
 [sys.executable,'scripts/build_extended_anatomy.py','--append-manifest','data/derived/app/manifest.json'],
 [sys.executable,'scripts/build_lymph_network.py','--append'],
 [sys.executable,'scripts/audit_anatomy_coverage.py'],
 [sys.executable,'scripts/verify_anatomy_fidelity.py'],
 [str(root/'.venv-physiology/bin/python'),'scripts/export_betse_tissue.py'],
 [sys.executable,'scripts/calibrate_skin.py'],
 [sys.executable,'-c','from ihm.coupling.circuit import build_skin_circuit; build_skin_circuit(".")'],
 [sys.executable,'scripts/build_temporal_atlas.py'],
 [sys.executable,'scripts/build_coupling_atlas.py'],
 [sys.executable,'scripts/build_system_coverage.py'],
 [sys.executable,'scripts/build_evidence_catalog.py'],
 [sys.executable,'scripts/audit_native_targets.py'],
 ['npm','--prefix','app','run','build']]
commands += [[sys.executable,'scripts/build_canonical_anatomy.py','--append']]
if args.native_adapter:commands += [[sys.executable,'scripts/build_native_adapter.py']]
commands += [
 [sys.executable,'scripts/build_body_profile.py'],
 [sys.executable,'scripts/build_body_brain.py'],
 [sys.executable,'scripts/build_body_mechanics.py'],
 [sys.executable,'scripts/build_canonical_body.py'],
 [sys.executable,'scripts/verify_canonical_anatomy.py'],
 [sys.executable,'scripts/verify_body_brain.py'],
 [sys.executable,'scripts/verify_body_mechanics.py'],
 [sys.executable,'scripts/verify_body_certainty.py'],
 [sys.executable,'scripts/verify_canonical_body.py'],
 [sys.executable,'-m','ihm','integrated','--output','artifacts/integrated-human.json']]
if args.plan:
 print(json.dumps(commands,indent=2));raise SystemExit(0)
logs=root/'artifacts/build';logs.mkdir(parents=True,exist_ok=True)
for i,command in enumerate(commands):
 print('BUILD '+' '.join(command),flush=True)
 with (logs/f'{i:02d}.log').open('w') as log:
  result=subprocess.run(command,cwd=root,stdout=log,stderr=subprocess.STDOUT)
 if result.returncode:raise SystemExit(f'Build failed; inspect {logs}/{i:02d}.log')
print('Built and checked canonical body, source-derived workbench views and artifacts/integrated-human.json')
