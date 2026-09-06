"""Evidence-conditioned bounded graphs; no new native blood or flow owner."""
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from .details import microvascular_unit


def load_registry(root):
 root=Path(root).resolve();path=root/'data/sources/organ_microvascular_priors.json'
 registry=json.loads(path.read_text())
 for row in registry['sources']:
  source=(root/row['path']).resolve()
  if not source.is_relative_to(root) or hashlib.sha256(source.read_bytes()).hexdigest()!=row['sha256']:
   raise ValueError('Changed or escaping vascular prior source')
 return registry


def positive(value,label):
 if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:
  raise ValueError('Positive finite '+label+' required')
 return float(value)


def materialize_muscle(root, *, cohort='young_men', section_mm=(.2,.2), extent_mm=.5,
                       seed=0, density_per_mm2=None):
 """Match an observed section-density target, leaving 3D topology inferred.

 One capillary per transverse-plane crossing. Mean target by default; no
 probability distribution is inferred from mean/SD alone. Integer rounding
 produces a reported residual, never a claim of exact continuous density.
 """
 if isinstance(seed,bool) or not isinstance(seed,int) or seed<0:raise ValueError('Nonnegative integer seed required')
 if len(section_mm)!=2:raise ValueError('Two transverse dimensions required')
 width,height=(positive(x,'section dimension') for x in section_mm)
 extent=positive(extent_mm,'longitudinal extent')
 registry=load_registry(root)
 if cohort not in registry['organs']['Muscle']['cohorts']:raise ValueError('Unknown human muscle cohort')
 prior=registry['organs']['Muscle']['cohorts'][cohort]
 density=positive(prior['density']['mean'] if density_per_mm2 is None else density_per_mm2,'density')
 area=positive(width*height,'section area');expected=positive(density*area,'expected crossing count')
 count=int(math.floor(expected+.5))
 if not 1<=count<=128:raise ValueError('Bounded fixture requires 1..128 capillary crossings')
 # Reuse topology only. Discard its arbitrary radius coefficients entirely.
 positions,edges,unused_radii,kinds=microvascular_unit([0.,0.,0.],extent*1e-3,count,seed)
 # Generator's random transverse box is .9*extent by .4*extent.
 positions[:,1]*=width/(.9*extent);positions[:,2]*=height/(.4*extent)
 crossing=(positions[edges[:,0],0]*positions[edges[:,1],0]<0)
 if int(crossing.sum())!=count or not np.all(kinds[crossing]==2):raise ValueError('Unexpected transverse topology')
 achieved=count/area
 path=Path(root)/'data/sources/organ_microvascular_priors.json'
 return {'schema':'ihm.conditional_muscle_microvascular.v1','evidence_kind':'conditional_synthetic_graph',
  'species':'Homo sapiens','region':'vastus lateralis prior; donor-local synthetic fixture','cohort':cohort,'seed':seed,
  'registry_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'source_receipts':registry['sources'],
  'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  'topology_generator_sha256':hashlib.sha256(Path(__file__).with_name('details.py').read_bytes()).hexdigest(),
  'positions_m':positions.tolist(),'edges':edges.tolist(),'edge_kind':kinds.tolist(),'radius_m':None,'flow_solution':None,
  'native_owner':'Muscle','independent_store':False,'native_flow_owned':False,'native_feedback_enabled':False,
  'ownership':'Geometry only; native volume, mass and perfusion remain native-owned and unallocated here',
  'conditioning':{'section_mm':[width,height],'extent_mm':extent,'target_kind':'cohort_mean' if density_per_mm2 is None else 'caller_scenario_override','cohort_summary':prior},
  'density_constraint':{'target_per_mm2':density,'achieved_per_mm2':achieved,'residual_per_mm2':achieved-density,'area_mm2':area,'plane_intersections':count,'plane':'x=0','integer_rounding_bound_per_mm2':.5/area},
  'topology_assumption':'Inherited paired binary supply/return trees; connectivity inferred, not measured human branching',
  'unmet_constraints':['domain_area_distribution','capillary_per_fiber','log_domain_sd','3D_radius','tortuosity','permeability','registered_fiber_alignment','native_terminal_correspondence'],
  'body_registered':False,'population_distribution_sampled':False}


def solve_passive(edges, resistance_pa_s_per_m3, pressure_boundaries_pa):
 """Optional scenario-only Kirchhoff solve using explicit edge resistances.

 This is a diagnostic boundary-value solve, not native perfusion or feedback.
 No viscosity, radius or pressure coefficients are supplied by this module.
 """
 raw=np.asarray(edges)
 if raw.ndim!=2 or raw.shape[1]!=2 or not len(raw) or raw.dtype.kind not in 'iu':raise ValueError('Integer edges required')
 edges=raw.astype(int)
 if edges.min()<0 or np.any(edges[:,0]==edges[:,1]):raise ValueError('Invalid edge endpoints')
 n=int(edges.max())+1
 if n>1024 or len(edges)>2048:raise ValueError('Diagnostic graph exceeds bounded size')
 r=np.asarray(resistance_pa_s_per_m3,float)
 if r.shape!=(len(edges),) or not np.isfinite(r).all() or np.any(r<=0):raise ValueError('Positive finite edge resistance required')
 bounds=dict(pressure_boundaries_pa)
 if not bounds or any(isinstance(k,bool) or not isinstance(k,int) or not 0<=k<n or isinstance(v,bool) or not math.isfinite(v) for k,v in bounds.items()):raise ValueError('Finite indexed pressure boundaries required')
 # Every node must be connected to a boundary; reject floating components.
 adj=[[] for _ in range(n)]
 for a,b in edges:adj[a].append(b);adj[b].append(a)
 seen=set(bounds);pending=list(bounds)
 while pending:
  for v in adj[pending.pop()]:
   if v not in seen:seen.add(v);pending.append(v)
 if len(seen)!=n:raise ValueError('Unanchored component')
 with np.errstate(over='ignore',under='ignore',divide='ignore'):g=1/r
 if not np.isfinite(g).all() or np.any(g<=0):raise ValueError('Unrepresentable conductance')
 scale=float(g.max());matrix=np.zeros((n,n))
 for (a,b),c in zip(edges,g/scale):matrix[a,a]+=c;matrix[b,b]+=c;matrix[a,b]-=c;matrix[b,a]-=c
 pressure=np.zeros(n);fixed=np.array(sorted(bounds),int);free=np.array([x for x in range(n) if x not in bounds],int)
 pressure[fixed]=[bounds[int(x)] for x in fixed]
 if len(free):pressure[free]=np.linalg.solve(matrix[np.ix_(free,free)],-matrix[np.ix_(free,fixed)]@pressure[fixed])
 flow=g*(pressure[edges[:,0]]-pressure[edges[:,1]])
 residual=np.zeros(n)
 for (a,b),q in zip(edges,flow):residual[a]+=q;residual[b]-=q
 return {'pressure_pa':pressure.tolist(),'flow_m3_per_s':flow.tolist(),
  'max_internal_residual_m3_per_s':float(np.max(np.abs(residual[free]))) if len(free) else 0.,
  'boundary_outflow_m3_per_s':{str(int(i)):float(residual[i]) for i in fixed},
  'evidence_kind':'caller_resistance_pressure_scenario','native_flow_owned':False,'native_feedback_enabled':False}
