"""Acquire exact retained primary-source bytes or fail visibly on publisher changes."""
from pathlib import Path
import hashlib,json,urllib.request
ROOT=Path(__file__).resolve().parents[1]
def acquire(root=ROOT):
 root=Path(root);card=json.loads((root/'data/sources/human-hair-mechanics.json').read_text());directory=root/'data/raw/hair';directory.mkdir(parents=True,exist_ok=True)
 for name,receipt in card['acquisition_receipts'].items():
  path=directory/name
  if path.exists():payload=path.read_bytes()
  else:
   request=urllib.request.Request(receipt['url'],headers={'User-Agent':'IHM research source acquisition'})
   with urllib.request.urlopen(request,timeout=60) as response:payload=response.read()
  if len(payload)!=receipt['bytes'] or hashlib.sha256(payload).hexdigest()!=receipt['sha256']:raise ValueError('Publisher bytes changed: '+name)
  if not path.exists():path.write_bytes(payload)
  print(name,receipt['sha256'])
if __name__=='__main__':acquire()
