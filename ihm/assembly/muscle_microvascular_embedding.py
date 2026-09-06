"""On-demand registered vascular patches and geometric zoom, never new stores."""
import gzip
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from .vascular_priors import materialize_muscle, positive


def vector(value):
 v=np.asarray(value,dtype=float)
 if v.shape!=(3,) or not np.isfinite(v).all():raise ValueError('Finite 3-vector required')
 return v


def clip_segment_box(start,end,lower,upper):
 """Closed AABB slab intersection in original segment parameters; no raster mask."""
 a,b,lo,hi=map(vector,(start,end,lower,upper))
 if np.any(lo>=hi):raise ValueError('Positive box extent required')
 t0,t1=0.,1.;delta=b-a
 for k in range(3):
  if delta[k]==0:
   if a[k]<lo[k] or a[k]>hi[k]:return None
  else:
   x,y=sorted(((lo[k]-a[k])/delta[k],(hi[k]-a[k])/delta[k]));t0=max(t0,x);t1=min(t1,y)
   if t0>t1:return None
 return float(t0),float(t1)


def point_surface_certificate(point,vertices,faces):
 """Numerical winding and analytic point/triangle distance on exact-coordinate weld.

 Closure is tested after exact coordinate equality only. This establishes a
 local authored-surface interior, not exclusive biological tissue occupancy.
 """
 q=vector(point);v=np.asarray(vertices,float);raw=np.asarray(faces)
 if v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all() or raw.ndim!=2 or raw.shape[1]!=3 or raw.dtype.kind not in 'iu' or not len(raw):raise ValueError('Finite triangular surface required')
 f=raw.astype(int)
 if f.min()<0 or f.max()>=len(v):raise ValueError('Invalid triangle index')
 unique,index=np.unique(v,axis=0,return_inverse=True);wf=index[f]
 directed=np.concatenate((wf[:,[0,1]],wf[:,[1,2]],wf[:,[2,0]]));undirected=np.sort(directed,axis=1)
 _,inverse,count=np.unique(undirected,axis=0,return_inverse=True,return_counts=True)
 signs=np.where(directed[:,0]<directed[:,1],1,-1);balance=np.bincount(inverse,weights=signs)
 if np.any(count!=2) or np.any(balance!=0) or np.any(directed[:,0]==directed[:,1]):raise ValueError('Surface not closed and consistently wound after exact weld')
 triangles=v[f];a,b,c=triangles[:,0],triangles[:,1],triangles[:,2];ab=b-a;ac=c-a;n=np.cross(ab,ac);norm2=np.einsum('ij,ij->i',n,n)
 if np.any(norm2<=0):raise ValueError('Degenerate triangles have unresolved containment')
 aq=q-a;signed=np.einsum('ij,ij->i',aq,n);projection=q-signed[:,None]*n/norm2[:,None]
 ap=projection-a;d00=np.einsum('ij,ij->i',ab,ab);d01=np.einsum('ij,ij->i',ab,ac);d11=np.einsum('ij,ij->i',ac,ac)
 d20=np.einsum('ij,ij->i',ap,ab);d21=np.einsum('ij,ij->i',ap,ac);den=d00*d11-d01*d01
 if np.any(den<=0):raise ValueError('Ill-conditioned triangle containment')
 u=(d11*d20-d01*d21)/den;w=(d00*d21-d01*d20)/den
 within=(u>=0)&(w>=0)&(u+w<=1);distance=np.where(within,signed*signed/norm2,np.inf)
 for x,y in ((a,b),(b,c),(c,a)):
  edge=y-x;t=np.clip(np.einsum('ij,ij->i',q-x,edge)/np.einsum('ij,ij->i',edge,edge),0,1)
  difference=q-(x+t[:,None]*edge);distance=np.minimum(distance,np.einsum('ij,ij->i',difference,difference))
 clearance=float(math.sqrt(float(distance.min())))
 ra=a-q;rb=b-q;rc=c-q;la=np.linalg.norm(ra,axis=1);lb=np.linalg.norm(rb,axis=1);lc=np.linalg.norm(rc,axis=1)
 numerator=np.einsum('ij,ij->i',ra,np.cross(rb,rc));denominator=la*lb*lc+np.einsum('ij,ij->i',ra,rb)*lc+np.einsum('ij,ij->i',rb,rc)*la+np.einsum('ij,ij->i',rc,ra)*lb
 winding=float(math.fsum(2*np.arctan2(numerator,denominator))/(4*math.pi))
 tol=1e-9
 return {'winding_number':winding,'clearance_m':clearance,'inside':clearance>tol and abs(abs(winding)-1)<1e-7,
         'numerical_margin_m':tol,'method':'Exact-coordinate welded closed surface; solid-angle winding and point-triangle distance',
         'physical_tissue_occupancy_verified':False,'self_intersections_globally_verified':False}


def materialize_patch(root, *, entity_id, position_m, resolution_m, seed=0, cohort='young_men',
                      section_mm=(.2,.2), extent_mm=.5, hydraulic_scenario=None, zoom_bounds_m=None):
 """Embed one stable patch; changing zoom resolution never regenerates its physics.

 The atlas principal axis is an explicit orientation inference. A clearance
 ball enclosing all centerlines plus largest radius certifies no surface
 crossing. Boundary-near requests fail; physiological graphs are not clipped.
 """
 root=Path(root).resolve();position=vector(position_m);resolution=positive(resolution_m,'zoom sampling resolution')
 anatomy_path=root/'data/derived/canonical/anatomy.json';anatomy_bytes=anatomy_path.read_bytes();anatomy=json.loads(anatomy_bytes)
 matches=[e for e in anatomy['entities'] if e['id']==entity_id]
 if len(matches)!=1 or matches[0].get('role')!='muscle':raise ValueError('Canonical muscle entity required')
 entity=matches[0]
 if 'vastus lateralis' not in entity['name'].lower():raise ValueError('This prior is registered only to vastus lateralis; other muscles require explicit transfer evidence')
 source=entity['reference_geometry'];path=(root/source['path']).resolve()
 if not path.is_relative_to(root):raise ValueError('Escaping geometry path')
 raw=path.read_bytes()
 if hashlib.sha256(raw).hexdigest()!=source['sha256'] or source.get('units')!='m':raise ValueError('Changed or non-SI canonical geometry')
 surface=json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw);vertices=np.asarray(surface['positions'],float).reshape(-1,3);faces=np.asarray(surface['indices'],int).reshape(-1,3)
 certificate=point_surface_certificate(position,vertices,faces)
 if not certificate['inside']:raise ValueError('Requested position is not certified inside canonical muscle')
 graph=materialize_muscle(root,cohort=cohort,section_mm=section_mm,extent_mm=extent_mm,seed=seed,hydraulic_scenario=hydraulic_scenario)
 local=np.asarray(graph['positions_m']);axis=vector(entity['principal_axis']);axis/=np.linalg.norm(axis)
 reference=np.eye(3)[int(np.argmin(np.abs(axis)))];second=np.cross(axis,reference);second/=np.linalg.norm(second);third=np.cross(axis,second);rotation=np.column_stack((axis,second,third))
 max_radius=0. if graph['radius_m'] is None else max(graph['radius_m']);sphere=float(np.linalg.norm(local,axis=1).max())+max_radius
 margin=certificate['clearance_m']-sphere-certificate['numerical_margin_m']
 if margin<=0:raise ValueError('Entire patch capsule is not certified inside; request a smaller or deeper patch')
 world=local@rotation.T+position;graph['positions_m']=world.tolist();graph['body_registered']=True;graph['coordinate_frame']=source['frame'];graph['region']=entity['name']+'; conditional synthetic patch'
 request={'entity_id':entity_id,'position_m':position.tolist(),'seed':seed,'cohort':cohort,'section_mm':list(section_mm),'extent_mm':extent_mm,'hydraulic_scenario':hydraulic_scenario,'source_geometry_sha256':source['sha256'],'registry_sha256':graph['registry_sha256'],'generator_sha256':graph['generator_sha256'],'topology_generator_sha256':graph['topology_generator_sha256'],'anatomy_sha256':hashlib.sha256(anatomy_bytes).hexdigest(),'embedding_module_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 patch_id='muscle-microvascular-'+hashlib.sha256(json.dumps(request,sort_keys=True).encode()).hexdigest()[:24]
 certificate.update(enclosing_capsule_ball_radius_m=sphere,clearance_margin_m=margin,whole_patch_inside_authored_surface=True,radius_containment='realized cylinder radius conservatively bounded by sphere' if graph['radius_m'] is not None else 'centerlines only; radius unresolved')
 bounds=None
 if zoom_bounds_m is not None:
  if len(zoom_bounds_m)!=2:raise ValueError('Two zoom box corners required')
  bounds=[vector(x) for x in zoom_bounds_m]
  if np.any(bounds[0]>=bounds[1]):raise ValueError('Positive zoom box required')
 zoom=[];sample_count=0;pressures=None if graph['flow_solution'] is None else graph['flow_solution']['pressure_pa']
 for i,(a,b) in enumerate(graph['edges']):
  interval=(0.,1.) if bounds is None else clip_segment_box(world[a],world[b],*bounds)
  if interval is None or interval[0]==interval[1]:continue
  t0,t1=interval;p0=world[a]+t0*(world[b]-world[a]);p1=world[a]+t1*(world[b]-world[a]);length=float(np.linalg.norm(p1-p0))
  needed=length/resolution
  if not math.isfinite(needed) or needed>100000:raise ValueError('Zoom sampling exceeds bounded size')
  segments=max(1,int(math.ceil(needed)));sample_count+=segments+1
  if sample_count>100000:raise ValueError('Zoom sampling exceeds bounded 100000 vertices')
  pressure=None if pressures is None else [pressures[a]+t*(pressures[b]-pressures[a]) for t in interval]
  zoom.append({'edge_id':patch_id+':'+str(i),'edge_index':i,'source_nodes':[a,b],'source_t':[t0,t1],
   'endpoints_m':[p0.tolist(),p1.tolist()],'radius_m':None if graph['radius_m'] is None else graph['radius_m'][i],
   'endpoint_pressure_pa':pressure,'flow_m3_per_s':None if graph['flow_solution'] is None else graph['flow_solution']['flow_m3_per_s'][i],
   'centerline_samples_m':np.linspace(p0,p1,segments+1).tolist(),'clipped_view_only':t0>0 or t1<1,'independent_store':False})
 return {'schema':'ihm.registered_muscle_microvascular_patch.v1','patch_id':patch_id,'request':request,'graph':graph,'zoom_edges':zoom,
  'resolution_m':resolution,'zoom_bounds_m':None if bounds is None else [x.tolist() for x in bounds],
  'native_owner':'Muscle','independent_store':False,'native_volume_allocated_ml':None,'macro_boundary_correspondence':None,
  'registration':{'entity_id':entity_id,'source_geometry':source,'anatomy_sha256':hashlib.sha256(anatomy_bytes).hexdigest(),
   'embedding_module_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'frame':source['frame'],'local_to_body_rotation':rotation.tolist(),
   'translation_m':position.tolist(),'evidence':'Synthetic patch positioned in retained canonical authored surface; principal-axis direction is inferred, not measured fibers',
   'source_specimen':entity.get('provenance',{}).get('source',{}).get('specimen'),'source_license':entity.get('provenance',{}).get('license')},
  'containment':certificate,'ownership':'On-demand native-owner view; geometric volume is unallocated demand, never added to native blood. Zoom fragments are views of original edges.',
  'boundary_scope':'Supply/return graph nodes carry only declared scenario pressure; macro-vessel and native terminal correspondence unresolved'}
