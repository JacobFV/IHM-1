"""Small numerical and ownership tests; no native run or body-sized array."""
from pathlib import Path
import sys,unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.body_exchange import NativeTissueExchange

def snapshot():
 v={}
 for pool in ('Skin.vascular','Skin.extracellular','Skin.intracellular','Lymph','VenaCava'):
  p='tissue.'+pool;v[p+'.volume_ml']=100.;v[p+'.pressure_mmhg']=2.
  for s in ('Albumin','Glucose','Oxygen','CarbonDioxide','Sodium','Potassium','Chloride'):
   v[p+'.'+s+'.mass_g']=.1;v[p+'.'+s+'.concentration_g_per_l']=1.
 for path,q in [('SkinE1ToSkinE2',.3),('SkinE3ToSkinI',.1),('SkinE3ToSkinL1',.2),('LymphToVenaCava',.2),('SkinSweating',0.)]:v['tissue.path.'+path+'.flow_ml_per_s']=q
 return {'time_s':1.,'elapsed_s':1.,'values':v}
class Checks(unittest.TestCase):
 def test_partition_exact_no_extra_native_store(self):
  x=NativeTissueExchange(organs=('Skin',),partitions={'Skin':[{'id':'patch','fraction':.2}]});s=snapshot();d=x.observe(s)
  for pool in ('vascular','extracellular','intracellular'):
   rows=[r for r in d['partitions'] if r['owner']=='Skin.'+pool]
   self.assertAlmostEqual(sum(r['volume_ml'] for r in rows),100.)
   self.assertAlmostEqual(sum(r['mass_g']['Albumin'] for r in rows),.1)
  self.assertEqual(s,snapshot());self.assertEqual(len(d['native_compartments']),5)
 def test_fluid_incidence_cancels_every_internal_transfer(self):
  d=NativeTissueExchange(organs=('Skin',)).observe(snapshot());self.assertAlmostEqual(sum(d['internal_volume_rate_ml_per_s'].values()),0.)
  self.assertAlmostEqual(d['internal_volume_rate_ml_per_s']['Skin.extracellular'],0.)
 def test_missing_and_negative_fail_closed(self):
  x=NativeTissueExchange(organs=('Skin',));s=snapshot();del s['values']['tissue.Skin.extracellular.volume_ml']
  with self.assertRaises(ValueError):x.observe(s)
  s=snapshot();s['values']['tissue.Skin.extracellular.Glucose.mass_g']=-1e-9
  with self.assertRaises(ValueError):x.observe(s)
 def test_absent_species_stays_missing(self):
  s=snapshot();del s['values']['tissue.Skin.intracellular.Albumin.mass_g']
  d=NativeTissueExchange(organs=('Skin',)).observe(s)
  self.assertIsNone(d['native_compartments']['Skin.intracellular']['mass_g']['Albumin'])
  self.assertIsNone(next(r for r in d['partitions'] if r['owner']=='Skin.intracellular')['mass_g']['Albumin'])
 def test_overlapping_partition_budget_rejected(self):
  with self.assertRaises(ValueError):NativeTissueExchange(organs=('Skin',),partitions={'Skin':[{'id':'a','fraction':.8},{'id':'b','fraction':.8}]})
 def test_no_false_albumin_flux_or_whole_body_closure(self):
  d=NativeTissueExchange(organs=('Skin',)).observe(snapshot());self.assertFalse(d['whole_body_mass_closure_claimed']);self.assertEqual(d['solute_fluxes'],[])
def retained(directory):
 import json,math,hashlib,tempfile
 from ihm.assembly.skin_bioelectric import native_skin_ionic_boundary
 directory=Path(directory).resolve();first=last=None
 for line in (directory/'receipts.jsonl').read_text().splitlines():
  record=json.loads(line)
  if 'acknowledgment' in record:
   last=record['acknowledgment']
   if first is None:first=last
 identity=json.loads((directory/'manifest.json').read_text())
 model=NativeTissueExchange.from_workspace(ROOT,first,identity);view=model.observe(last);networks=model.project_networks(view);surfaces=model.project_anatomical_supports(view)
 max_volume=max_mass=0.
 for owner,quantity in view['native_compartments'].items():
  if owner in ('Lymph','VenaCava'):continue
  regions=[r for r in view['partitions'] if r['owner']==owner]
  max_volume=max(max_volume,abs(math.fsum(r['volume_ml'] for r in regions)-quantity['volume_ml']))
  for sub,mass in quantity['mass_g'].items():
   if mass is not None:max_mass=max(max_mass,abs(math.fsum(r['mass_g'][sub] for r in regions)-mass))
 assert max_volume<1e-10 and max_mass<1e-12
 for network in networks:
  parent=next(r for r in view['partitions'] if r['id']==network['id']+'.vascular');assert abs(math.fsum(network['edge_volume_ml'])-parent['volume_ml'])<1e-12
 ionic=native_skin_ionic_boundary(view)
 out=Path(tempfile.mkdtemp(prefix='body-exchange-',dir=ROOT/'data/derived/audits'))
 report={'passed':True,'native_record':str(directory.relative_to(ROOT)),'native_receipts_sha256':hashlib.sha256((directory/'receipts.jsonl').read_bytes()).hexdigest(),'source_hashes':model.source_hashes,'time_s':view['time_s'],'native_owners':len(view['native_compartments']),'runtime_partitions':len(view['partitions']),'on_demand_source_surface_partitions':len(surfaces),'source_fine_networks':len(networks),'max_partition_volume_residual_ml':max_volume,'max_partition_mass_residual_g':max_mass,'selected_internal_flow_incidence_residual_ml_per_s':math.fsum(view['internal_volume_rate_ml_per_s'].values()),'maximum_native_mass_concentration_residual_g':max(abs(r) for c in view['native_compartments'].values() for r in c['mass_concentration_residual_g'].values() if r is not None),'skin_ionic_boundary':ionic,'unlocalized_organs':[o for o,rs in model.partitions.items() if any(r['evidence_kind']=='unlocalized_native_owner' for r in rs)],'scope':'Read-only actual native snapshot materialization; no new native run, compression validation or integrated solute ledger claim'}
 (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print(out)
if __name__=='__main__':
 if len(sys.argv)==3 and sys.argv[1]=='--retained':retained(sys.argv[2])
 else:unittest.main()
