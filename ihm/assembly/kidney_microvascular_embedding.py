"""Canonical kidney containment and polyline-preserving read-only zoom views."""
import gzip,hashlib,json,math
from pathlib import Path
import numpy as np
from .kidney_microvascular_patch import materialize_kidney_patch
from .muscle_microvascular_embedding import vector,point_surface_certificate,clip_segment_box
from .vascular_priors import positive


def _polyline_views(points,lower,upper,resolution):
 lengths=np.linalg.norm(np.diff(points,axis=0),axis=1);total=float(lengths.sum());offset=np.r_[0,np.cumsum(lengths)]
 runs=[];current=[];parameters=[]
 def finish():
  if len(current)>1:runs.append((list(current),list(parameters)))
  current.clear();parameters.clear()
 for i,(a,b) in enumerate(zip(points[:-1],points[1:])):
  interval=(0.,1.) if lower is None else clip_segment_box(a,b,lower,upper)
  if interval is None or interval[0]==interval[1]:finish();continue
  x,y=interval;p=a+x*(b-a);q=a+y*(b-a);t0=(offset[i]+x*lengths[i])/total;t1=(offset[i]+y*lengths[i])/total
  if current and np.linalg.norm(np.asarray(current[-1])-p)>1e-12:finish()
  if not current:current.append(p.tolist());parameters.append(float(t0))
  n=max(1,int(math.ceil(float(np.linalg.norm(q-p))/resolution)))
  if n>20000:raise ValueError('Kidney zoom sampling exceeds limit')
  for k in range(1,n+1):current.append((p+(q-p)*k/n).tolist());parameters.append(float(t0+(t1-t0)*k/n))
 finish();return runs


def materialize_registered_kidney(root, *, entity_id,position_m,resolution_m,seed,scenario,zoom_bounds_m=None):
 root=Path(root).resolve();position=vector(position_m);resolution=positive(resolution_m,'resolution');p=root/'data/derived/canonical/anatomy.json';raw=p.read_bytes();anatomy=json.loads(raw)
 entities=[x for x in anatomy['entities'] if x['id']==entity_id]
 if len(entities)!=1 or entities[0]['name'] not in ('left kidney','right kidney'):raise ValueError('Canonical kidney required')
 entity=entities[0];source=entity['reference_geometry'];p=(root/source['path']).resolve()
 if not p.is_relative_to(root):raise ValueError('Escaping source geometry')
 data=p.read_bytes()
 if hashlib.sha256(data).hexdigest()!=source['sha256'] or source['units']!='m':raise ValueError('Changed or non-SI kidney geometry')
 mesh=json.loads(gzip.decompress(data));v=np.array(mesh['positions']).reshape(-1,3);f=np.array(mesh['indices']).reshape(-1,3)
 certificate=point_surface_certificate(position,v,f)
 if not certificate['inside']:raise ValueError('Position not certified inside canonical kidney')
 owner='LeftKidney' if entity['name']=='left kidney' else 'RightKidney';graph=materialize_kidney_patch(root,owner=owner,seed=seed,scenario=scenario)
 axis=vector(entity['principal_axis']);axis/=np.linalg.norm(axis);second=np.cross(axis,np.eye(3)[np.argmin(np.abs(axis))]);second/=np.linalg.norm(second);rotation=np.column_stack((axis,second,np.cross(axis,second)))
 radius=max(graph['radius_m']);sphere=max(float(np.linalg.norm(points,axis=1).max()) for points in graph['edge_polylines_m'])+radius
 margin=certificate['clearance_m']-sphere-certificate['numerical_margin_m']
 if margin<=0:raise ValueError('Kidney patch polylines/radii exceed certified surface clearance')
 graph['positions_m']=(np.asarray(graph['positions_m'])@rotation.T+position).tolist()
 graph['edge_polylines_m']=[(np.asarray(x)@rotation.T+position).tolist() for x in graph['edge_polylines_m']]
 graph['body_registered']=True;graph['coordinate_frame']=source['frame']
 graph['flow_solution']={'pressure_pa':graph['physical_node_pressure_pa'],'flow_m3_per_s':graph['physical_edge_flow_m3_per_s'],'max_internal_residual_m3_per_s':graph['hydraulics']['max_internal_residual_m3_per_s']}
 identity={'entity_id':entity_id,'position_m':position.tolist(),'seed':seed,'scenario':scenario,'geometry_sha256':source['sha256'],'anatomy_sha256':hashlib.sha256(raw).hexdigest(),'registry_sha256':graph['source_registry_sha256'],'generator_sha256':graph['generator_sha256'],'embedding_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 patch_id='kidney-microvascular-'+hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:24]
 lower=upper=None
 if zoom_bounds_m is not None:
  lower,upper=map(vector,zoom_bounds_m)
  if np.any(lower>=upper):raise ValueError('Invalid zoom bounds')
 rows=[];total_samples=0
 for i,(a,b) in enumerate(graph['edges']):
  points=np.asarray(graph['edge_polylines_m'][i]);runs=_polyline_views(points,lower,upper,resolution)
  pa,pb=graph['physical_node_pressure_pa'][a],graph['physical_node_pressure_pa'][b]
  for fragment,(samples,parameters) in enumerate(runs):
   total_samples+=len(samples)
   if total_samples>20000 or len(rows)>=256:raise ValueError('Kidney zoom fragment/sample limit exceeded')
   rows.append({'edge_id':patch_id+':'+str(i)+':'+str(fragment),'source_edge_id':patch_id+':'+str(i),'edge_index':i,'fragment_index':fragment,'source_nodes':[a,b],
    'source_t':[parameters[0],parameters[-1]],'endpoints_m':[samples[0],samples[-1]],'centerline_samples_m':samples,'sample_source_t':parameters,
    'radius_m':graph['radius_m'][i],'edge_kind':graph['edge_kind'][i],'endpoint_pressure_pa':[pa+(pb-pa)*parameters[0],pa+(pb-pa)*parameters[-1]],
    'flow_m3_per_s':graph['physical_edge_flow_m3_per_s'][i],'clipped_view_only':parameters[0]>0 or parameters[-1]<1,'independent_store':False})
 certificate.update(enclosing_capsule_ball_radius_m=sphere,clearance_margin_m=margin,whole_patch_inside_authored_surface=True,cortical_location_verified=False)
 return {'schema':'ihm.registered_kidney_microvascular_patch.v1','patch_id':patch_id,'request':identity,'graph':graph,'zoom_edges':rows,'resolution_m':resolution,'zoom_bounds_m':zoom_bounds_m,
  'native_owner':owner,'independent_store':False,'native_volume_allocated_ml':None,'macro_boundary_correspondence':None,'containment':certificate,
  'registration':{'entity_id':entity_id,'frame':source['frame'],'source_geometry':source,'anatomy_sha256':hashlib.sha256(raw).hexdigest(),'translation_m':position.tolist(),'local_to_body_rotation':rotation.tolist(),
   'evidence':'Whole canonical kidney containment only; cortical layer and measured afferent/native vessel registration unresolved; atlas principal-axis orientation inferred'},
  'boundary_scope':'Four declared diagnostic pressure/exchange boundaries; no native perfusion/filtration ownership',
  'ownership':'Read-only native kidney view; geometric demand is not another store. Zoom fragments refer to original arc-length segments.'}
