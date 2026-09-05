"""Materialize scalp and body strands without changing canonical/app manifests."""
from pathlib import Path
import argparse,gzip,hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.build_body_details import load_skin_geometry,write_npz
from ihm.assembly.details import sample_hair
from ihm.assembly.hair import display_indices

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def tube(centers,offsets,radii):
 p=[];indices=[]
 for s,(a,b) in enumerate(zip(offsets[:-1],offsets[1:])):
  tangent=centers[a+1]-centers[a];tangent/=np.linalg.norm(tangent);ref=np.array([1.,0,0]) if abs(tangent[1])>.9 else np.array([0.,1,0]);u=np.cross(tangent,ref);u/=np.linalg.norm(u);w=np.cross(tangent,u)
  for j in range(a,b):
   for k in range(4):p.append(centers[j]+radii[s]*(np.cos(k*np.pi/2)*u+np.sin(k*np.pi/2)*w))
   if j<b-1:
    for k in range(4):indices.extend([[j*4+k,j*4+(k+1)%4,(j+1)*4+k],[j*4+(k+1)%4,(j+1)*4+(k+1)%4,(j+1)*4+k]])
 return np.asarray(p).ravel().tolist(),np.asarray(indices).ravel().tolist()

def build(out,scalp_count=64,body_count=32,scalp_render_count=2048,body_render_count=512):
 out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
 if any(out.iterdir()):raise ValueError('Use a fresh output directory')
 reference_path=ROOT/'data/measurements/hair/strand_reference.json';reference=json.loads(reference_path.read_text());regions=reference['scalp_geometry']['regions'];rho=reference['density']['value_kg_m3']
 anatomy_path=ROOT/'data/derived/canonical/anatomy.json';anatomy=json.loads(anatomy_path.read_text());skin=next(e for e in anatomy['entities'] if e['id']=='body-bp3d-FJ2810');source,g=load_skin_geometry(skin);v=np.asarray(g['positions']).reshape(-1,3);f=np.asarray(g['indices']).reshape(-1,3);tri=v[f];center=tri.mean(axis=1);cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);norm=np.linalg.norm(cross,axis=1);normal=cross/np.maximum(norm[:,None],1e-30)
 # Atlas-frame scalp coverage prior, not a measured scalp segmentation.
 def scalp_mask(p):return (abs(p[...,0])<.105)&((p[...,1]>.79)|((p[...,1]>.72)&(p[...,2]<.015)))
 radial=center-np.array([0,.755,0]);radial/=np.maximum(np.linalg.norm(radial,axis=1)[:,None],1e-30);mask=np.all(scalp_mask(tri),axis=1)&(np.einsum('ij,ij->i',radial,normal)>.45)
 region=np.where(center[:,2]>.045,'frontal',np.where(center[:,2]<-.025,'occipital','vertex'));densities={k:r['density_per_cm2'] for k,r in regions.items()};diameters={k:r['diameter_um']*1e-6 for k,r in regions.items()};density=np.array([densities[x] for x in region])*mask
 scalp=sample_hair(v,f,density,seed=20260906);scalp['region']=region[scalp['face_index']];scalp['radius_m']=np.array([diameters[x]/2 for x in scalp['region']]);scalp['length_m']=np.full(len(scalp['ids']),.03);write_npz(out/'scalp_population.npz',scalp)
 body_source=ROOT/'data/derived/canonical/hair_samples.npz';body=dict(np.load(body_source,allow_pickle=False));structures=[];summaries={}
 if str(body['source_geometry_sha256'].item())!=sha(source) or str(body['source_entity_id'].item())!=skin['id']:raise ValueError('Body root source mismatch')
 def materialize(data,limit,nodes):
  selected=display_indices(data,limit);roots=data['roots_m'][selected];directions=data['normals'][selected];lengths=data['length_m'][selected];radii=data['radius_m'][selected];offsets=np.arange(len(selected)+1)*nodes;centers=(roots[:,None,:]+directions[:,None,:]*lengths[:,None,None]*np.linspace(0,1,nodes)[None,:,None]).reshape(-1,3)
  triangles=tri[data['face_index'][selected]]
  if not np.allclose(np.einsum('ni,nij->nj',data['barycentric'][selected],triangles),roots,rtol=0,atol=1e-9):raise ValueError('Root barycentric registration mismatch')
  attachment={'kind':'ElasticHairMaterialPoint','skin_entity_id':skin['id'],'source_geometry_sha256':sha(source),'sample_ids':data['ids'][selected].tolist(),'barycentric':data['barycentric'][selected].tolist(),'reference_triangles_m':triangles.ravel().tolist(),'root_rule':'triangle barycentric position; face normal sets clamped first segment; free nodes dynamically simulated'}
  strands={'centerlines_m':centers.ravel().tolist(),'strand_offsets':offsets.tolist(),'radius_m':radii.tolist(),'tensile_modulus_pa':reference['tension']['value_pa'],'bending_modulus_pa':reference['bending']['reference_value_pa'],'density_kg_m3':rho,'mechanical_scope':'Mixed human scalp reference measurements transferred to a circular short straight beam; body-vellus material transfer unvalidated','small_deflection_slope_limit':.3}
  return selected,roots,lengths,radii,offsets,centers,attachment,strands
 for name,data,limit,render_limit,nodes in [('scalp',scalp,scalp_count,scalp_render_count,7),('body',body,body_count,body_render_count,3)]:
  if not 1<=limit<=512 or not limit<=render_limit<=4096:raise ValueError('Invalid guide/render budget')
  selected,roots,lengths,radii,offsets,centers,attachment,strands=materialize(data,limit,nodes)
  rs,rr,rl,rad,ro,rc,ra,rstrands=materialize(data,render_limit,nodes);positions,indices=tube(rc,ro,rad)
  # Nearest guides are restricted to the same retained anatomical region.
  interpolation=[]
  for j,point in enumerate(rr):
   candidates=np.flatnonzero(data['region'][selected]==data['region'][rs[j]])
   if not len(candidates):raise ValueError('Guide budget missed a region')
   distance=np.linalg.norm(roots[candidates]-point,axis=1);order=np.argsort(distance)[:3];ids=candidates[order];distance=distance[order]
   weights=np.array([1.]+[0.]*(len(ids)-1)) if distance[0]<1e-12 else (1/distance)/(1/distance).sum()
   interpolation.append([[int(i),float(w)] for i,w in zip(ids,weights)])
  sid=f'body-dynamic-{name}-hair';geometry={'positions':positions,'indices':indices,'strands':strands,'attachment':attachment,'render_strands':rstrands,'render_attachment':ra,'guide_interpolation':interpolation};path=out/(sid+'.json.gz');path.write_bytes(gzip.compress(json.dumps(geometry,separators=(',',':'),allow_nan=False).encode(),mtime=0))
  structures.append({'id':sid,'name':f'{name.title()} hair — physical elastic strands','model_id':'ihm-body','system':'hair','kind':'mesh','color':'#382b20','default_visible':True,'evidence_kind':'measured_material_reference_with_atlas_root_prior','canonical_entity_id':skin['id'],'geometry_url':'/api/geometry/'+sid,'geometry_path':str(path.relative_to(ROOT)),'geometry_sha256':sha(path)})
  summaries[name]={'population_count':len(data['ids']),'simulated_guide_count':len(selected),'render_fiber_count':len(rs),'display_density_fraction':len(rs)/len(data['ids']),'node_count':len(centers),'radius_range_m':[float(radii.min()),float(radii.max())],'length_range_m':[float(lengths.min()),float(lengths.max())],'physical_radius_multiplier':1,'coverage':'scalp atlas-coordinate/outward-face mask prior' if name=='scalp' else 'held regional body density and vellus morphology priors','length_evidence':'30mm haircut scenario, not a measured individual' if name=='scalp' else 'held .5–2mm scenario prior','total_simulated_guide_mass_kg':float(np.sum(rho*np.pi*radii**2*lengths))}
 report={'structures':structures,'populations':summaries,'source_geometry_sha256':sha(source),'source_geometry_path':str(source.relative_to(ROOT)),'source_hashes':{str(p.relative_to(ROOT)):sha(p) for p in [source,body_source,anatomy_path,reference_path,ROOT/'data/sources/human-hair-mechanics.json',ROOT/'app/src/hair_dynamics.js',Path(__file__)]},'body_population_sha256':sha(body_source),'anatomy_sha256':sha(anatomy_path),'builder_sha256':sha(Path(__file__)),'source_card':'data/sources/human-hair-mechanics.json','source_card_sha256':sha(ROOT/'data/sources/human-hair-mechanics.json'),'measurement_sha256':sha(reference_path),'maximum_update_hz':10,'render_motion':'three nearest same-region guide displacements; not independent rendered-fiber dynamics','scalp_mask_area_m2':float(np.sum(norm[mask]/2)),'limitations':['Independent strand motion; no hair-hair collisions, torsion, growth or fluid drag law.','Root triangles are actual source surface; scalp/body region masks and transfer to canonical subject remain explicit priors.','Displayed sample does not thicken fibers or represent full population dynamics.']}
 (out/'manifest_fragment.json').write_text(json.dumps(report,indent=2)+'\n');print(out);print(json.dumps(summaries,indent=2));return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',default='data/derived/hair/elastic_v3');p.add_argument('--scalp-count',type=int,default=64);p.add_argument('--body-count',type=int,default=32);a=p.parse_args();build(a.output,a.scalp_count,a.body_count)
