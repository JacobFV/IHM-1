"""Bounded conditional cutaneous vessels attached to exact eligible source faces.

A source-conditioned geometry/flow view, never a second native blood store.
The paired plexus template, depth, rheology and pressure inputs are engineering.
"""
import gzip
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from .vascular_priors import solve_passive

TERRITORY_PATH='data/research/engineered_skin_territories/materialization.json'
TERRITORY_SHA='de669f96f6db5b150dd467f149464677e19f22ec2e7be5a77d895d873b3fef1a'


def default_scenario():
 """Explicit engineering scenario; these constants are not source measurements."""
 return {'id':'forearm_baseline_transfer_to_lower_leg_v1',
         'epidermal_thickness_m':.0001,'loop_apex_depth_m':.00015,
         'superficial_plexus_depth_m':.0003,'deep_plexus_depth_m':.0014,
         'dermal_bottom_depth_m':.0016,'arteriole_radius_m':20e-6,
         'venule_radius_m':25e-6,'viscosity_pa_s':.003,
         'arterial_pressure_pa':4500.,'venous_pressure_pa':1500.,
         'allow_forearm_to_lower_leg_transfer':True}


def _read(root,path,sha=None):
 b=(root/path).read_bytes()
 if sha is not None and hashlib.sha256(b).hexdigest()!=sha:raise ValueError('Changed source bytes: '+path)
 return b


def materialize_skin_patch(root,*,territory_id,scenario,source_triangle_id=None,
                           patch_area_mm2=1.,seed=0):
 """Materialize 1..56 papillary loops on one original exterior-proxy triangle.

 Whole source mesh is loaded without simplification. The triangle footprint is a
 centroid-scaled source face, so every projected centerline stays on that face.
 No physical tissue containment, packing or macro perfusion binding is inferred.
 """
 root=Path(root)
 if type(seed) is not int or not 0<=seed<2**32:raise ValueError('Seed must be an unsigned 32-bit integer')
 if isinstance(patch_area_mm2,bool) or not isinstance(patch_area_mm2,(int,float)) or not math.isfinite(patch_area_mm2) or not .15<=patch_area_mm2<=8:raise ValueError('Patch area must be .15..8 mm2')
 expected=default_scenario()
 if not isinstance(scenario,dict) or set(scenario)!=set(expected) or scenario['id']!=expected['id'] or scenario['allow_forearm_to_lower_leg_transfer'] is not True:raise ValueError('Explicit supported scenario and anatomical-site transfer required')
 s=dict(scenario)
 for k in expected.keys()-{'id','allow_forearm_to_lower_leg_transfer'}:
  if isinstance(s[k],bool) or not isinstance(s[k],(int,float)) or not math.isfinite(s[k]):raise ValueError('Finite scenario numbers required')
 if not 0<s['epidermal_thickness_m']<s['loop_apex_depth_m']<s['superficial_plexus_depth_m']<s['deep_plexus_depth_m']<s['dermal_bottom_depth_m']<=.005:raise ValueError('Ordered bounded skin layer depths required')
 if not 1e-6<=s['arteriole_radius_m']<=100e-6 or not 1e-6<=s['venule_radius_m']<=100e-6 or not .0005<=s['viscosity_pa_s']<=.01:raise ValueError('Radius or viscosity outside engineering bounds')
 if not 0<=s['venous_pressure_pa']<s['arterial_pressure_pa']<=20000:raise ValueError('Ordered pressure scenario required')
 regpath='data/sources/skin_microvascular_priors.json';regraw=_read(root,regpath);priors=json.loads(regraw)
 for receipt in priors['primary_receipts']:_read(root,receipt['path'],receipt['sha256'])
 material=json.loads(_read(root,TERRITORY_PATH,TERRITORY_SHA))
 regions=[r for r in material['regions'] if r['id']==territory_id]
 if len(regions)!=1:raise ValueError('Exact engineered territory ID required')
 receipt=material['source_receipts'][0];raw=_read(root,receipt['path'],receipt['sha256'])
 mesh=json.loads(gzip.decompress(raw) if receipt['path'].endswith('.gz') else raw)
 xyz=np.asarray(mesh['positions'],float).reshape(-1,3);faces=np.asarray(mesh['indices'],int).reshape(-1,3)
 ids=regions[0]['triangle_ids'];eligible=set(material['contact_eligible_triangle_ids'])
 if source_triangle_id is not None and (type(source_triangle_id) is not int or source_triangle_id not in ids or source_triangle_id not in eligible):raise ValueError('Face is not eligible in selected territory')
 area=patch_area_mm2*1e-6
 chosen=None
 for fid in ids if source_triangle_id is None else [source_triangle_id]:
  if fid not in eligible:continue
  tri=xyz[faces[fid]];cross=np.cross(tri[1]-tri[0],tri[2]-tri[0]);twice=float(np.linalg.norm(cross))
  if twice<=0:continue
  scale=math.sqrt(area/(twice/2))
  # A margin around the projected footprint protects finite vessel radius at face boundaries.
  maxside=max(np.linalg.norm(tri[(i+1)%3]-tri[i]) for i in range(3))
  margin=(1-scale)*twice/maxside/3
  if scale<=.8 and margin>max(s['arteriole_radius_m'],s['venule_radius_m'],9.59e-6/2):chosen=(fid,tri,cross/twice,scale,margin);break
 if chosen is None:raise ValueError('No source triangle supports requested bounded footprint and radius margin')
 fid,tri,normal,scale,margin=chosen
 prior=priors['human_confocal'];density=prior['density']['value'];radius=prior['lumen_diameter']['value']*.5e-6
 count=round(density*patch_area_mm2)
 if not 1<=count<=56:raise ValueError('Loop count outside 1..56')
 bary=[];depths=[];nodes=[];edges=[];rows=[]
 def surface(b):return (1-scale)/3+scale*np.asarray(b,float)
 def node(b,depth):
  bc=surface(b);bary.append(bc.tolist());depths.append(depth);nodes.append((bc@tri-depth*normal).tolist());return len(nodes)-1
 def edge(a,b,r,kind,polyline=None):
  pts=np.asarray([nodes[a],nodes[b]] if polyline is None else polyline)
  length=float(np.linalg.norm(np.diff(pts,axis=0),axis=1).sum())
  if length<=0:raise ValueError('Degenerate generated vessel')
  edges.append([a,b]);rows.append({'id':f'skin:{territory_id}:face:{fid}:seed:{seed}:edge:{len(edges)-1}',
   'source_edge_identity':'conditional_template_not_individual_measured_vessel','vessel_class':kind,
   'endpoints_m':[nodes[a],nodes[b]],'centerline_samples_m':pts.tolist(),
   'radius_m':r,'length_m':length,'resistance_pa_s_per_m3':8*s['viscosity_pa_s']*length/(math.pi*r**4)})
 deep=s['deep_plexus_depth_m'];sup=s['superficial_plexus_depth_m']
 aroot=node([.7,.2,.1],deep);ar=node([.5,.3,.2],deep);a=node([.45,.35,.2],sup)
 v=node([.2,.35,.45],sup);vr=node([.2,.3,.5],deep);vroot=node([.1,.2,.7],deep)
 edge(aroot,ar,s['arteriole_radius_m'],'deep_arterial_plexus')
 edge(ar,a,s['arteriole_radius_m'],'ascending_arteriole')
 edge(v,vr,s['venule_radius_m'],'descending_venule')
 edge(vr,vroot,s['venule_radius_m'],'deep_venous_plexus')
 rng=np.random.default_rng(seed)
 for i in range(count):
  center=.15+.55*rng.dirichlet([2,2,2]);ba=center+np.array([.03,-.03,0]);bv=center-np.array([.03,-.03,0])
  ca=node(ba,sup);cv=node(bv,sup)
  edge(a,ca,s['arteriole_radius_m'],'superficial_arterial_plexus')
  ts=np.linspace(0,1,33);path=[]
  for t in ts:
   d=sup-(sup-s['loop_apex_depth_m'])*math.sin(math.pi*t)
   path.append(surface(ba*(1-t)+bv*t)@tri-d*normal)
  edge(ca,cv,radius,'papillary_capillary_loop',path)
  edge(cv,v,s['venule_radius_m'],'superficial_venous_plexus')
 # Finite radii must remain within the declared nominal dermal slab.
 for row in rows:
  d=-(np.asarray(row['centerline_samples_m'])-tri[0])@normal
  if np.any(d-row['radius_m']<s['epidermal_thickness_m']) or np.any(d+row['radius_m']>s['dermal_bottom_depth_m']):raise ValueError('Vessel radius crosses nominal dermal boundary')
 hydraulic=solve_passive(np.asarray(edges,int),[r['resistance_pa_s_per_m3'] for r in rows],{aroot:s['arterial_pressure_pa'],vroot:s['venous_pressure_pa']})
 for i,((a,b),row) in enumerate(zip(edges,rows)):
  row['endpoint_pressure_pa']=[hydraulic['pressure_pa'][a],hydraulic['pressure_pa'][b]];row['flow_m3_per_s']=hydraulic['flow_m3_per_s'][i]
 identity=hashlib.sha256(json.dumps({'face':fid,'territory_sha':TERRITORY_SHA,'priors_sha':hashlib.sha256(regraw).hexdigest(),'scenario':s,'seed':seed,'area':area},sort_keys=True).encode()).hexdigest()
 for row in rows:row['id']=identity+':'+row['id']
 return {'schema':'ihm.skin-microvascular-patch.v1','id':identity,'seed':seed,'scenario':s,
  'registration':{'entity_id':'body-bp3d-FJ2810','frame':'bodyparts3d-display-m','source_geometry':receipt,
   'source_triangle_id':fid,'source_triangle_positions_m':tri.tolist(),'territory_id':territory_id,
   'territory_path':TERRITORY_PATH,'territory_sha256':TERRITORY_SHA,'exterior_eligible':True,
   'exterior_basis':'largest positive-integral open component engineering prior',
   'outward_normal_prior':normal.tolist(),'minimum_footprint_face_margin_m':margin,
   'physical_tissue_containment_verified':False,'local_depth_basis':'inward original face normal; nominal scenario layers'},
  'source_conditioning':priors,'source_registry_sha256':hashlib.sha256(regraw).hexdigest(),
  'nodes_m':nodes,'node_surface_barycentric':bary,'node_depth_m':depths,'edges':edges,'zoom_edges':rows,
  'hydraulics':hydraulic,'constraints':{'patch_area_m2':area,'loop_count':count,'target_loop_density_per_mm2':density,'realized_loop_density_per_mm2':count/patch_area_mm2,
   'capillary_diameter_m':2*radius,'radius_distribution':'constant source mean; reported SD is not sampled as vessel variation',
   'deep_plexus_and_superficial_connectivity':'engineering paired star plexus and papillary loop template',
   'packing_or_vessel_intersection_validated':False},
  'ownership':{'prospective_native_owner':'Skin','native_boundary_binding':None,'macro_vessel_binding':None,
   'binding_ready':False,'additional_native_volume_m3':0.,'independent_blood_store':False,
   'geometry_lumen_volume_m3':sum(math.pi*r['radius_m']**2*r['length_m'] for r in rows),
   'geometry_volume_role':'diagnostic local geometry only, never added to native Skin blood volume'},
  'lymphatics_materialized':False,'exchange_model':None,'solute_or_charge_model':None,
  'native_commands':[],'read_only':True,'limits':{'maximum_loops':56,'maximum_area_mm2':8,'whole_source_mesh_preserved':True}}
