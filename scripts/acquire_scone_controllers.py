"""Acquire a pinned predictive-locomotion controller source, with byte receipts."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
REVISION='39c4b7a52a90dff838ebe7a1cb8afe58af1aaae1'
URL='https://github.com/tgeijten/scone-core.git'
out=ROOT/'data/raw/mechanics/scone-core'
source=out/'source'
if not source.exists():
    source.mkdir(parents=True)
    subprocess.run(['git','init',str(source)],check=True)
    subprocess.run(['git','-C',str(source),'remote','add','origin',URL],check=True)
    subprocess.run(['git','-C',str(source),'fetch','--depth','1','origin',REVISION],check=True)
    subprocess.run(['git','-C',str(source),'checkout','--detach',REVISION],check=True)
actual=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
if actual!=REVISION:raise RuntimeError('Existing controller source revision differs')
if subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip():
    raise RuntimeError('Controller source is modified')
tracked=subprocess.check_output(['git','-C',str(source),'ls-files','-z']).decode().split('\0')
files={}
for name in tracked:
    path=source/name
    if name and path.is_file():
        files[name]={'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
record={'source_id':'scone-predictive-locomotion','repository':URL,'revision':REVISION,
        'retrieved_utc':datetime.now(timezone.utc).isoformat(),'files':files,
        'submodules_acquired':False,'source_only':True,'built':False,
        'role':'Published controller implementations and source muscle/contact model variants; not a generic-body gait calibration',
        'license_path':'source/LICENSE','third_party_notice_path':'source/THIRD_PARTY_NOTICES.md'}
receipt_path=out/'acquisition.json'
if receipt_path.exists():
    previous=json.loads(receipt_path.read_text())
    if previous['revision']!=actual or previous['files']!=files:raise RuntimeError('Existing acquisition receipt disagrees with donor bytes')
else:receipt_path.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'revision':actual,'files':len(files),'bytes':sum(f['bytes'] for f in files.values())},indent=2))
