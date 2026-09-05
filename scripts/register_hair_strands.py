"""Register verified hair assets, replacing the old rigid display only."""
from pathlib import Path
import hashlib,json,shutil,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.app.experiments import read_experiment

def register(root=ROOT):
 root=Path(root);fragment=read_experiment(root,'hair-strands');manifest_path=root/'data/derived/app/manifest.json';manifest=json.loads(manifest_path.read_bytes());old=[s for s in manifest['structures'] if s['id']=='body-detail-hair']
 replaced={'body-detail-hair',*(s['id'] for s in fragment['structures'])};manifest['structures']=[s for s in manifest['structures'] if s['id'] not in replaced]
 for row in fragment['structures']:
  row=dict(row);region='scalp' if 'scalp' in row['id'] else 'body';population=fragment['populations'][region]
  row.update(population_count=population['population_count'],display_count=population['render_fiber_count'],simulated_guides=population['simulated_guide_count'],display_density_fraction=population['display_density_fraction'],physical_radius_multiplier=1,maximum_update_hz=10,material_source='10.3390/molecules25092143; 10.3390/cosmetics6020024',material_condition='Human scalp reference; 5.7 GPa bending is a representative example, not cohort mean. Body-vellus transfer unvalidated.',geometry_certainty=population['coverage'],scalp_geometry_cohort='120 Arabian adults, 60 men/60 women, age 18–60; 10.2147/CCID.S394045',collision_model='none',experiment_url='/api/body/experiments/hair-strands')
  target=root/'data/derived/app/geometry'/(row['id']+'.json.gz');shutil.copyfile(root/row['geometry_path'],target)
  assert hashlib.sha256(target.read_bytes()).hexdigest()==row['geometry_sha256'];manifest['structures'].append(row)
 manifest['hair_strands']={'experiment_url':'/api/body/experiments/hair-strands','replaced_display_ids':['body-detail-hair'],'full_population_dynamics':False,'retained_old_display':old or manifest.get('hair_strands',{}).get('retained_old_display',[])}
 temp=manifest_path.with_suffix('.tmp');temp.write_text(json.dumps(manifest,separators=(',',':'))+'\n');temp.replace(manifest_path)
 print(json.dumps({'registered':[r['id'] for r in fragment['structures']],'replaced_display_ids':['body-detail-hair'],'source_population_preserved':True}))
if __name__=='__main__':register()
