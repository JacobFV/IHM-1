"""Verify detached allocation, preserved support, and hair equivalence receipts."""
import argparse,json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from prepare_skin_layer_migration import digest
from ihm.assembly.contact_dynamics import DynamicTetrahedra


def verify(root,epoch):
 root=Path(root).resolve();epoch=Path(epoch).resolve();r=json.loads((epoch/'acceptance.json').read_bytes())
 for relative,expected in r['outputs'].items():
  p=(epoch/relative).resolve()
  if not p.is_relative_to(epoch) or digest(p.read_bytes())!=expected:raise ValueError('Consumer output changed')
 for domain in r['material_domains']:
  p=root/domain['output'];m=json.loads((p/'manifest.json').read_bytes())
  with np.load(root/domain['source']/'pelvic-domain.npz',allow_pickle=False) as old,np.load(p/'pelvic-domain.npz',allow_pickle=False) as new:
   for key in domain['unchanged_arrays']:
    if not np.array_equal(old[key],new[key]):raise ValueError('Retained domain geometry changed')
   b=DynamicTetrahedra(new['vertices_m'],new['tetrahedra'],mu_pa=new['mu_pa'],lambda_pa=new['lambda_pa'],density_kg_m3=new['density_kg_m3'],fixed_nodes=new['fixed_nodes'])
   if not np.isclose(b.mass_kg.sum(),m['mass_claim_kg'],atol=1e-12,rtol=1e-10):raise ValueError('Mass claim differs')
   if b.max_explicit_dt_s!=m['max_explicit_dt_s']:raise ValueError('Stability limit differs')
  if m['mass_activation_contract']['canonical_runtime_handoff_applied'] is not False:raise ValueError('Unexpected domain activation')
 support=json.loads((epoch/'support/manifest.json').read_bytes());parent=json.loads((root/'data/derived/lumbar-supine-reference-1qex2x9i/contact/manifest.json').read_bytes())
 for key in ('arrays_path','native_input_path'):
  if (root/support[key]).read_bytes()!=(root/parent[key]).read_bytes():raise ValueError('Frozen support changed')
 if support['accepted_support'] or support['native_integration']:raise ValueError('Unsupported support acceptance')
 hair=json.loads((epoch/'hair-equivalence.json').read_bytes())
 for path,entry in hair['assets'].items():
  if digest((root/path).read_bytes())!=entry['sha256']:raise ValueError('Hair asset changed')
 if hair['roots_regenerated'] or hair['native_mass_modified']:raise ValueError('Unexpected hair mutation claim')
 return {'domains':len(r['material_domains']),'hair_assets':len(hair['assets']),'support_bytes_equal':True,'native_executed':False}

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('epoch',type=Path);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args();print(json.dumps(verify(a.root,a.epoch)))
