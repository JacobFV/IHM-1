"""Small fail-closed hair receipt and frozen-public-materialization checks."""
from pathlib import Path
import hashlib,json,sys,tempfile
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.app.experiments import read_experiment
from ihm.human import ImplicitHuman

def verify():
 with tempfile.TemporaryDirectory(prefix='ihm-hair-receipts-') as tmp:
  root=Path(tmp);directory=root/'data/derived/hair/elastic_v3';directory.mkdir(parents=True);(root/'app/src').mkdir(parents=True);module=root/'app/src/hair_dynamics.js';module.write_text('fixture');geometry=root/'geometry.json';geometry.write_text('{}');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
  manifest=directory/'manifest_fragment.json';manifest.write_text(json.dumps({'source_hashes':{'app/src/hair_dynamics.js':sha(module)},'structures':[{'geometry_path':'geometry.json','geometry_sha256':sha(geometry)}]}))
  body=ImplicitHuman.open(root);assert len(body.materialize('hair-strands')['structures'])==1
  geometry.write_text('changed')
  try:read_experiment(root,'hair-strands')
  except ValueError:pass
  else:raise AssertionError('Changed geometry accepted')
  geometry.write_text('{}');module.write_text('changed')
  for call in [lambda:read_experiment(root,'hair-strands'),lambda:body.materialize('hair-strands')]:
   try:call()
   except ValueError:pass
   else:raise AssertionError('Changed source accepted')
 print('PASS: valid materialization, changed geometry, changed source, frozen public source')
if __name__=='__main__':verify()
