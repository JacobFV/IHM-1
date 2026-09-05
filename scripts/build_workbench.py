"""Rebuild portable views from the acquired corpus in dependency order.

This does not redownload datasets, compile native engines, or rerun patient
simulations. Use the documented source-specific collectors/builders for those.
"""
from pathlib import Path
import subprocess
import sys
import json
root=Path(__file__).resolve().parents[1]
commands=[
 [sys.executable,'scripts/build_spatial_atlas.py'],
 [sys.executable,'scripts/build_opensim_display.py'],
 [str(root/'.venv-physiology/bin/python'),'scripts/export_betse_tissue.py'],
 [sys.executable,'scripts/calibrate_skin.py'],
 [sys.executable,'-c','from ihm.coupling.circuit import build_skin_circuit; build_skin_circuit(".")'],
 [sys.executable,'scripts/build_temporal_atlas.py'],
 [sys.executable,'scripts/build_coupling_atlas.py'],
 [sys.executable,'scripts/build_system_coverage.py'],
 [sys.executable,'scripts/build_evidence_catalog.py'],
 [sys.executable,'scripts/audit_native_targets.py'],
 ['npm','--prefix','app','run','build'],
 [sys.executable,'-m','ihm','integrated','--output','artifacts/integrated-human.json']]
logs=root/'artifacts/build';logs.mkdir(parents=True,exist_ok=True)
for i,command in enumerate(commands):
 print('BUILD '+' '.join(command),flush=True)
 with (logs/f'{i:02d}.log').open('w') as log:
  result=subprocess.run(command,cwd=root,stdout=log,stderr=subprocess.STDOUT)
 if result.returncode:raise SystemExit(f'Build failed; inspect {logs}/{i:02d}.log')
print('Built source-derived workbench views and artifacts/integrated-human.json')
