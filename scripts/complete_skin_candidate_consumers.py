#!/usr/bin/env python3
"""Rebind bounded retained consumer artifacts to a separate skin metadata epoch."""
from copy import deepcopy
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from prepare_skin_layer_migration import ANATOMY,MECHANICS,encode,digest
from ihm.assembly.contact_dynamics import DynamicTetrahedra
from ihm.assembly.material_domains import MaterialOwnership

HAIR='data/derived/hair/elastic_v3'
SUPPORT='data/derived/lumbar-supine-reference-1qex2x9i/contact/manifest.json'


def read_checked(root,path,expected=None):
 p=(root/path).resolve()
 if not p.is_relative_to(root):raise ValueError('Source outside root')
 raw=p.read_bytes()
 if expected is not None and digest(raw)!=expected:raise ValueError('Source digest differs: '+path)
 return raw


def geometry_plan(root,anatomy):
 refs={}
 for e in anatomy['entities']:
  r=e['reference_geometry'];p=r['path']
  if p in refs and refs[p]['sha256']!=r['sha256']:raise ValueError('Conflicting geometry identity')
  source=(root/p).resolve()
  if not source.is_relative_to(root) or not source.is_file():raise ValueError('Invalid geometry path')
  refs[p]={'sha256':r['sha256'],'bytes':source.stat().st_size}
 return {'schema':'ihm.viewer-geometry-staging-plan.v1','files':refs,'file_count':len(refs),
  'total_bytes':sum(v['bytes'] for v in refs.values()),'contents_hashed':False,'copied':False,
  'method':'Explicit byte-budgeted hash pass then independent reflink/copy and destination hash verification; no shared writable hardlinks'}


def material_candidate(root,stage,out,source):
 manifest_raw=read_checked(root,source+'/manifest.json');manifest=json.loads(manifest_raw)
 registry=json.loads(read_checked(root,source+'/registry.json',manifest['artifacts']['registry.json']))
 npz=read_checked(root,source+'/pelvic-domain.npz',manifest['artifacts']['pelvic-domain.npz'])
 if len(npz)>4*1024*1024 or manifest['tetrahedra']>60000:raise ValueError('Domain exceeds bounded retained-array budget')
 import io
 with np.load(io.BytesIO(npz),allow_pickle=False) as data:arrays={k:data[k].copy() for k in data.files}
 mechanics=json.loads((stage/MECHANICS).read_bytes());anatomy=json.loads((stage/ANATOMY).read_bytes());specs={e['id']:e for e in mechanics['entities']};entities={e['id']:e for e in anatomy['entities']}
 regions=deepcopy(manifest['material_regions']);ratios=[]
 for region in regions:
  ident=region['source_id'];spec=specs[ident];surface=next(s for s in manifest['source_surfaces'] if s['id']==ident)
  if surface['source_sha256']!=entities[ident]['reference_geometry']['sha256']:raise ValueError('Retained tetra source geometry differs')
  read_checked(root,surface['path'],surface['source_sha256'])
  if region['source_material']!=spec['material']:raise ValueError('Constitutive material changed')
  selected=arrays['material_index']==regions.index(region)
  if not np.allclose(arrays['density_kg_m3'][selected],region['effective_inertial_density_kg_m3'],rtol=1e-12,atol=0):raise ValueError('Old density array differs')
  ratio=spec['mass_kg']/region['allocated_mass_kg'];ratios.append(ratio)
  region['allocated_mass_kg']=spec['mass_kg'];region['effective_inertial_density_kg_m3']=spec['mass_kg']/region['volume_m3']
  arrays['density_kg_m3'][selected]=region['effective_inertial_density_kg_m3']
 body=DynamicTetrahedra(arrays['vertices_m'],arrays['tetrahedra'],mu_pa=arrays['mu_pa'],lambda_pa=arrays['lambda_pa'],density_kg_m3=arrays['density_kg_m3'],fixed_nodes=arrays['fixed_nodes'])
 for index,region in enumerate(regions):
  actual_volume=body.region.volumes[arrays['material_index']==index].sum()
  if not np.isclose(actual_volume,region['volume_m3'],rtol=1e-10,atol=1e-15):raise ValueError('Per-owner retained tetra volume differs')
 if not np.array_equal(body.triangles,arrays['boundary_triangles']):raise ValueError('Retained boundary topology changed')
 ownership=MaterialOwnership({r['source_id']:r['allocated_mass_kg'] for r in regions});ownership.claim(manifest['domain_id'],[r['source_id'] for r in regions])
 claimed=ownership.mass_kg(manifest['domain_id'])
 if not np.isclose(body.mass_kg.sum(),claimed,rtol=1e-10,atol=1e-12):raise ValueError('Nodal mass does not match exclusive allocation')
 for row in registry['entities']:
  spec=specs[row['entity_id']]
  if row['source_geometry']!=entities[row['entity_id']]['reference_geometry'] or row['material_prior']!=spec['material']:raise ValueError('Registry geometry/material differs')
  row['mass_kg']=spec['mass_kg'];row['mass_role']=spec['mass_role']
 registry['single_body_mass_kg']=sum(r['mass_kg'] for r in registry['entities'])
 out.mkdir(parents=True);np.savez_compressed(out/'pelvic-domain.npz',**arrays);(out/'registry.json').write_bytes(encode(registry))
 result=deepcopy(manifest);result.update(material_regions=regions,mass_kg=float(body.mass_kg.sum()),mass_claim_kg=claimed,max_explicit_dt_s=body.max_explicit_dt_s)
 result['sources']={str((stage/p).relative_to(root)):digest((stage/p).read_bytes()) for p in (ANATOMY,MECHANICS)}
 result['parent_materialization']={'path':source+'/manifest.json','sha256':digest(manifest_raw)}
 result['migration_scope']='Retained tetra topology and constitutive priors; density/nodal allocation/stability limit recomputed. Detached only, no native debit.'
 result['artifacts']={name:digest((out/name).read_bytes()) for name in ('pelvic-domain.npz','registry.json')}
 (out/'manifest.json').write_bytes(encode(result))
 return {'source':source,'output':str(out.relative_to(root)),'tetrahedra':len(arrays['tetrahedra']),
  'old_mass_kg':manifest['mass_kg'],'new_mass_kg':result['mass_kg'],'density_ratios':ratios,
  'old_max_explicit_dt_s':manifest['max_explicit_dt_s'],'new_max_explicit_dt_s':result['max_explicit_dt_s'],
  'unchanged_arrays':[k for k in arrays if k!='density_kg_m3'],'native_mass_modified':False}


def hair_equivalence(root,stage):
 raw=read_checked(root,HAIR+'/manifest_fragment.json');old=json.loads(raw);anatomy=json.loads((stage/ANATOMY).read_bytes());skin=next(e for e in anatomy['entities'] if e['role']=='skin')
 if skin['reference_geometry']['sha256']!=old['source_geometry_sha256']:raise ValueError('Hair skin geometry differs')
 read_checked(root,old['source_geometry_path'],old['source_geometry_sha256'])
 # Preserve all root/population/render assets; hash rather than regenerate.
 assets={}
 for path in sorted((root/HAIR).iterdir()):
  if path.suffix in ('.npz','.gz'):
   assets[str(path.relative_to(root))]={'sha256':digest(path.read_bytes()),'bytes':path.stat().st_size}
 for structure in old['structures']:read_checked(root,structure['geometry_path'],structure['geometry_sha256'])
 for path,sha in old['source_hashes'].items():
  if path not in (ANATOMY,'scripts/build_hair_strands.py','app/src/hair_dynamics.js'):read_checked(root,path,sha)
 return {'schema':'ihm.hair-anatomy-equivalence.v1','parent_manifest':{'path':HAIR+'/manifest_fragment.json','sha256':digest(raw)},
  'historical_anatomy_sha256':old['anatomy_sha256'],'candidate_anatomy_sha256':digest((stage/ANATOMY).read_bytes()),
  'source_geometry_sha256':old['source_geometry_sha256'],'assets':assets,'roots_regenerated':False,'native_mass_modified':False,
  'acceptance':'Same source geometry and unchanged retained populations/root assets; historical hair metadata retained. Runtime/build source evolution is outside this geometric equivalence.'}


def support_candidate(root,stage,out):
 raw=read_checked(root,SUPPORT);old=json.loads(raw)
 arrays=read_checked(root,old['arrays_path'],old['arrays_sha256']);native=read_checked(root,old['native_input_path'],old['native_input_sha256'])
 specs={e['id']:e for e in json.loads((stage/MECHANICS).read_bytes())['entities']}
 layers=deepcopy(old['layers'])
 for layer in layers:
  spec=specs[layer['id']]
  if layer['thickness_m']!=spec['shell']['thickness_m'] or layer['young_modulus']!=spec['material']['young_modulus'] or layer['poisson_ratio']!=spec['material']['poisson_ratio']:raise ValueError('Support physical law differs')
  layer['thickness_basis']={'method':'explicit_shell_thickness','shell':deepcopy(spec['shell']),'physical_surface_support':deepcopy(spec['physical_surface_support'])}
 result=deepcopy(old);out.mkdir(parents=True)
 for name,content in [('quadrature.npz',arrays),('supine_surface_foundation.txt',native)]:(out/name).write_bytes(content)
 result.update(layers=layers,accepted_support=False,native_integration=False,
  arrays_path=str((out/'quadrature.npz').relative_to(root)),native_input_path=str((out/'supine_surface_foundation.txt').relative_to(root)))
 result['source_files']={SUPPORT:digest(raw),str((stage/MECHANICS).relative_to(root)):digest((stage/MECHANICS).read_bytes())}
 result['migration_scope']='Only explicit thickness/provenance metadata changed. Exact held quadrature and native text preserved, including seam omission. No new equilibrium claim.'
 (out/'manifest.json').write_bytes(encode(result))
 return {'output':str(out.relative_to(root)),'parent_manifest_sha256':digest(raw),'native_input_equal':True,'quadrature_equal':True,
  'omitted_unresolved_contact':old.get('omitted_unresolved_contact'),'accepted_support':False,'native_integration':False}


def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--epoch',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 root=a.root.resolve();stage=a.epoch.resolve()/'root';out=a.output.resolve()
 if out.exists() or not out.is_relative_to(root/'data/derived') or out.is_relative_to(root/'data/derived/canonical'):raise ValueError('Fresh isolated derived epoch required')
 from verify_skin_candidate_stage import verify
 verify(a.epoch);out.mkdir(parents=True)
 material=[material_candidate(root,stage,out/('pelvis-'+spacing), 'data/derived/material-domains/pelvis-'+spacing) for spacing in ('0.004m','0.002m')]
 hair=hair_equivalence(root,stage);(out/'hair-equivalence.json').write_bytes(encode(hair))
 support=support_candidate(root,stage,out/'support')
 plan=geometry_plan(root,json.loads((stage/ANATOMY).read_bytes()));(out/'viewer-geometry-plan.json').write_bytes(encode(plan))
 report={'schema':'ihm.skin-consumer-acceptance.v1','stage_epoch':str(a.epoch.resolve().relative_to(root)),
  'material_domains':material,'hair_equivalence':'hair-equivalence.json','support':support,'viewer_inventory_bytes':plan['total_bytes'],
  'outputs':{str(path.relative_to(out)):digest(path.read_bytes()) for path in out.rglob('*') if path.is_file()},
  'implementation_sha256':digest(Path(__file__).read_bytes()),'native_executed':False,'published':False}
 (out/'acceptance.json').write_bytes(encode(report));print(json.dumps({'output':str(out),'material_domains':material,'support':support,'viewer_bytes':plan['total_bytes']}))

if __name__=='__main__':main()
