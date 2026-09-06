"""Kidney-specific conditional serial vascular beds with explicit exchange ports."""
import hashlib,json,math
from pathlib import Path
import numpy as np
from .vascular_priors import positive,solve_passive


def load_kidney_priors(root):
 root=Path(root).resolve();path=root/'data/sources/kidney_microvascular_priors.json';registry=json.loads(path.read_text())
 for name,row in registry['source_receipts'].items():
  if 'sha256' not in row:continue
  p=root/'data/research/kidney_microvascular'/name;b=p.read_bytes()
  if len(b)!=row['bytes'] or hashlib.sha256(b).hexdigest()!=row['sha256']:raise ValueError('Changed kidney primary evidence')
 return registry,hashlib.sha256(path.read_bytes()).hexdigest()


def _length(points):return float(np.linalg.norm(np.diff(points,axis=0),axis=1).sum())


def _tortuous(a,b,target):
 direction=b-a;axis=np.eye(3)[np.argmin(np.abs(direction))];normal=np.cross(direction,axis);normal/=np.linalg.norm(normal)
 t=np.linspace(0,1,129);base=a[None,:]+t[:,None]*direction;wave=np.sin(32*math.pi*t);wave[[0,-1]]=0.
 def curve(amplitude):return base+amplitude*wave[:,None]*normal
 lo,hi=0.,25e-6
 if _length(curve(hi))<target:raise ValueError('Length target cannot fit bounded synthetic capillary')
 for _ in range(55):
  mid=(lo+hi)/2
  if _length(curve(mid))<target:lo=mid
  else:hi=mid
 return curve((lo+hi)/2)


def materialize_kidney_patch(root, *, owner='RightKidney', seed=0, scenario=None):
 """Source-conditioned chambers/lobular loops, then efferent/peritubular mesh.

 Scenario must declare viscosity, unmeasured PTC radius and four pressure
 boundaries. Exchange conductances may be zero; nonzero values are purely
 hydraulic engineering choices without oncotic/solute/charge calibration.
 """
 if owner not in ('LeftKidney','RightKidney'):raise ValueError('Kidney native owner required')
 if type(seed) is not int or seed<0:raise ValueError('Nonnegative integer seed required')
 required={'viscosity_pa_s','peritubular_radius_m','arterial_pressure_pa','venous_pressure_pa','bowman_pressure_pa','interstitial_pressure_pa','glomerular_hydraulic_conductance_m3_s_pa','peritubular_hydraulic_conductance_m3_s_pa'}
 if isinstance(scenario,dict) and 'peritubular_radius_mode' in scenario:required=(required-{'peritubular_radius_m'})|{'peritubular_radius_mode'}
 if not isinstance(scenario,dict) or set(scenario)!=required:raise ValueError('Complete explicit kidney hydraulic scenario required')
 priors,registry_hash=load_kidney_priors(root)
 values={k:positive(v,k) for k,v in scenario.items() if not k.endswith('conductance_m3_s_pa') and k!='peritubular_radius_mode'}
 radius_mode=scenario.get('peritubular_radius_mode','caller_engineering_radius')
 if radius_mode!='caller_engineering_radius':
  groups={'human_control_area_equivalent':'control','human_injury_area_equivalent':'injured'}
  if radius_mode not in groups:raise ValueError('Unknown human peritubular radius mode')
  area=priors['peritubular']['two_d_stats'][groups[radius_mode]]['luminal_area_um2']['mean']
  values['peritubular_radius_m']=math.sqrt(area/math.pi)*1e-6
 for name in required-set(values)-{'peritubular_radius_mode'}:
  x=scenario[name]
  if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x<0:raise ValueError('Nonnegative hydraulic conductance required')
  values[name]=float(x)
 rng=np.random.default_rng(seed)
 pos=[[-160e-6,-40e-6,0],[-60e-6,-40e-6,0],[-60e-6,40e-6,0],[-160e-6,40e-6,0],[-260e-6,40e-6,0],[-500e-6,40e-6,0],[-80e-6,-40e-6,0],[-80e-6,40e-6,0]]
 edges=[];kinds=[];radii=[];lobule_ids=[]
 def edge(a,b,kind,radius,lobule=None):edges.append([a,b]);kinds.append(kind);radii.append(radius);lobule_ids.append(lobule)
 def measured(kind):return priors['diameters_um'][kind]['mean']*.5e-6
 edge(0,6,'afferent_arteriole',measured('afferent_arteriole'));edge(6,1,'afferent_chamber',measured('afferent_chamber'))
 edge(2,7,'efferent_chamber',measured('efferent_chamber'));edge(7,3,'efferent_arteriole',measured('efferent_arteriole'));edge(3,4,'efferent_arteriole',measured('efferent_arteriole'))
 lobules=[];filtration_nodes=[]
 for l in range(7):
  theta=2*math.pi*l/7;center=np.array([45*math.cos(theta),45*math.sin(theta),0])*1e-6
  radial=np.array([math.cos(theta),math.sin(theta),0]);vertical=np.array([0,0,1.]);phase=rng.uniform(0,2*math.pi);nodes=[]
  for i in range(8):
   angle=2*math.pi*i/8+phase;nodes.append(len(pos));pos.append((center+22e-6*(math.cos(angle)*radial+math.sin(angle)*vertical)).tolist())
  lobules.append(nodes);filtration_nodes.extend(nodes)
  edge(1,nodes[0],'afferent_conduit',measured('afferent_conduit'),l)
  for i in range(8):edge(nodes[i],nodes[(i+1)%8],'filtration_capillary',priors['filtration_capillary_radius_um']['value']*1e-6,l)
  for a,b in [(1,5),(2,6)]:edge(nodes[a],nodes[b],'filtration_capillary',priors['filtration_capillary_radius_um']['value']*1e-6,l)
  for i in [3,7]:edge(nodes[i],2,'efferent_first_order',measured('efferent_first_order'),l)
 for l in range(7):edge(lobules[l][2],lobules[(l+1)%7][5],'filtration_capillary',priors['filtration_capillary_radius_um']['value']*1e-6,None)
 ptc=[]
 for row in range(3):
  for col in range(4):ptc.append(len(pos));pos.append([-260e-6-row*80e-6,-30e-6+col*40e-6,10e-6*math.sin(row+col)])
 r_ptc=values['peritubular_radius_m'];edge(4,ptc[0],'peritubular_capillary',r_ptc)
 for row in range(3):
  for col in range(4):
   i=row*4+col
   if row<2:edge(ptc[i],ptc[i+4],'peritubular_capillary',r_ptc)
   if col<3:edge(ptc[i],ptc[i+1],'peritubular_capillary',r_ptc)
 edge(ptc[-1],5,'peritubular_capillary',r_ptc)
 pos=np.asarray(pos);edge_array=np.asarray(edges);mask=np.array(kinds)=='filtration_capillary';chords=np.linalg.norm(pos[edge_array[:,1]]-pos[edge_array[:,0]],axis=1)
 target=priors['total_filtration_capillary_length_m']['value'];scale=target/float(chords[mask].sum())
 if scale<=1:raise ValueError('Source length target shorter than generated chords')
 polylines=[_tortuous(pos[a],pos[b],chords[i]*scale) if mask[i] else np.array([pos[a],pos[b]]) for i,(a,b) in enumerate(edges)]
 length=np.array([_length(x) for x in polylines]);radius=np.asarray(radii);resistance=8*values['viscosity_pa_s']*length/(math.pi*radius**4)
 if not np.isfinite(resistance).all() or np.any(resistance<=0):raise ValueError('Invalid hydraulic resistance')
 # Exchange ports are pressure reservoirs in the diagnostic, not vessel geometry.
 solver_edges=list(edges);solver_resistance=resistance.tolist();boundaries={0:values['arterial_pressure_pa'],5:values['venous_pressure_pa']};port_ids={'arterial':0,'venous':5,'bowman':None,'interstitial':None};exchange_indices={};next_node=len(pos)
 for name,nodes,coefficient in [('bowman',filtration_nodes,values['glomerular_hydraulic_conductance_m3_s_pa']),('interstitial',ptc,values['peritubular_hydraulic_conductance_m3_s_pa'])]:
  indices=[]
  if coefficient>0:
   port=next_node;next_node+=1;port_ids[name]=port;boundaries[port]=values[name+'_pressure_pa']
   for n in nodes:indices.append(len(solver_edges));solver_edges.append([n,port]);solver_resistance.append(len(nodes)/coefficient)
  exchange_indices[name]=indices
 solution=solve_passive(solver_edges,solver_resistance,boundaries);volume=math.pi*radius**2*length
 return {'schema':'ihm.conditional_kidney_microvascular.v1','native_owner':owner,'independent_store':False,'native_volume_allocated_ml':None,'native_feedback_enabled':False,
  'species':'Homo sapiens','body_registered':False,'coordinate_frame':'synthetic kidney patch local meters','seed':seed,
  'positions_m':pos.tolist(),'edges':edges,'edge_kind':kinds,'edge_lobule':lobule_ids,'edge_polylines_m':[x.tolist() for x in polylines],'edge_length_m':length.tolist(),'radius_m':radii,
  'edge_resistance_pa_s_per_m3':resistance.tolist(),'geometric_lumen_volume_m3':float(math.fsum(volume)),
  'hydraulics':solution,'physical_node_pressure_pa':solution['pressure_pa'][:len(pos)],'physical_edge_flow_m3_per_s':solution['flow_m3_per_s'][:len(edges)],'boundary_nodes':port_ids,
  'exchange_flows_m3_s':{name:math.fsum(solution['flow_m3_per_s'][i] for i in indices) for name,indices in exchange_indices.items()},
  'exchange_boundaries':{'bowman':{'pressure_pa':values['bowman_pressure_pa'],'conductance_m3_s_pa':values['glomerular_hydraulic_conductance_m3_s_pa']},'interstitial':{'pressure_pa':values['interstitial_pressure_pa'],'conductance_m3_s_pa':values['peritubular_hydraulic_conductance_m3_s_pa']}},
  'solute_fluxes':None,'charge_fluxes':None,'oncotic_model':None,'filtration_fraction_calibrated':False,
  'peritubular_radius_conditioning':{'mode':radius_mode,'realized_radius_m':r_ptc,'measured_radius_distribution':False},
  'conditioning':{'lobules':7,'afferent_conduits':7,'efferent_first_order_vessels':14,'afferent_diameter_um':priors['diameters_um']['afferent_arteriole']['mean'],'priors':priors,'scenario':values},
  'realized':{'filtration_capillary_length_m':float(math.fsum(length[mask])),'cycle_rank':len(edges)-len(pos)+1,'filtration_capillary_edges':int(mask.sum()),'physical_nodes':len(pos),'physical_edges':len(edges)},
  'source_registry_sha256':registry_hash,'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  'topology_evidence':'Human chamber/conduit and cyclic lobular organization conditions a new inferred graph; individual connectivity/positions and PTC lattice are engineering synthesis, not reconstructed donor vessels',
  'unmet_constraints':['measured individual topology','glomerular diameter and no-overlap packing','human PTC radius/3D geometry','measured filtration-capillary radius distribution','native terminal registration','oncotic and solute/charge transport'],
  'exchange_interpretation':'Optional purely hydraulic shunts to declared external pressures; not calibrated filtration, tubular reabsorption, or GFR. No new native pool or native mass transfer is created.'}
