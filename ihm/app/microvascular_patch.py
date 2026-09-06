"""Read-only bounded vascular detail service; no native actor or quantity mutation."""
from copy import deepcopy
import json
import math
from pathlib import Path
from threading import Lock
from ihm.assembly.muscle_microvascular_embedding import materialize_patch
from ihm.assembly.vascular_priors import load_registry
from ihm.assembly.kidney_microvascular_embedding import materialize_registered_kidney

_QUERY_LOCK=Lock()
LIMITS={'max_request_bytes':32768,'max_capillaries':64,'max_edges':256,'max_zoom_samples':20000,'max_response_bytes':2097152,
        'min_resolution_m':5e-6,'max_resolution_m':1e-3,'max_extent_mm':2.,'max_parallel_queries':1,'max_kidney_filtration_edges':96}
REGIONS={'body-bp3d-FJ1442':('right vastus lateralis',[-.14,-.23,0.]),
         'body-bp3d-FJ1442M':('left vastus lateralis',[.141294,-.23,0.]),
         'body-bp3d-FJ3145':('left kidney',[.075,.25,-.02]),'body-bp3d-FJ3147':('right kidney',[-.07,.23,-.012])}


class PatchRequestError(ValueError):
 def __init__(self,status,code,message):super().__init__(message);self.status=status;self.code=code
 def as_dict(self):return {'error':str(self),'code':self.code}


def _number(value,label,lo,hi):
 if isinstance(value,bool) or not isinstance(value,(float,int)) or not math.isfinite(value) or not lo<=value<=hi:
  raise PatchRequestError(400,'invalid_request',label+' outside supported finite range')
 return float(value)


def _vector(value,label):
 if not isinstance(value,(list,tuple)) or len(value)!=3:raise PatchRequestError(400,'invalid_request',label+' needs three coordinates')
 return [_number(x,label,-10.,10.) for x in value]


def _scenario(scenario_id):
 if scenario_id in ('human_kidney_control_v1','human_kidney_injury_v1'):
  return {'id':scenario_id,'label':'Human glomerular chambers and '+('preserved' if 'control' in scenario_id else 'injured')+' biopsy-region PTC area; declared hydraulic exchange',
   'hydraulic_scenario':{'viscosity_pa_s':.003,'peritubular_radius_mode':'human_control_area_equivalent' if 'control' in scenario_id else 'human_injury_area_equivalent','arterial_pressure_pa':12000.,'venous_pressure_pa':1000.,'bowman_pressure_pa':2500.,'interstitial_pressure_pa':500.,'glomerular_hydraulic_conductance_m3_s_pa':1e-17,'peritubular_hydraulic_conductance_m3_s_pa':1e-18},
   'interpretation':'Human chamber diameter means and PTC mean-area transform condition synthetic serial beds. Filtration radius comes from a published model assumption; topology, viscosity, pressures and hydraulic-only exchange are declared engineering choices, not native GFR/perfusion.'}
 ids={'human_quadriceps_1988_corrected_v1':'shrinkage_corrected','human_quadriceps_1988_uncorrected_v1':'uncorrected'}
 if scenario_id not in ids:raise PatchRequestError(400,'invalid_request','Unknown declared radius scenario')
 return {'id':scenario_id,'label':'Human quadriceps diameter means with '+ids[scenario_id]+' preparation; declared hydraulic assumptions',
  'hydraulic_scenario':{'radius_mode':'human_quadriceps_1988_equal_area','diameter_state':ids[scenario_id],
   'capillary_radius_cv':.2,'supply_radius_multiplier':2.,'return_radius_multiplier':2.4,
   'viscosity_pa_s':.003,'pressure_boundaries_pa':{0:4000.,1:1000.}},
  'interpretation':'Diameter means are human histological evidence. CV/lognormal shape, supply/return radius multipliers, viscosity and pressure boundaries are engineering choices; source does not identify their distribution or native correspondence.'}


class MicrovascularPatchService:
 """One nonblocking query at a time across service instances in this process.

 Every request rechecks source hashes. There is no cached stale evidence,
 native simulation, actor acquisition, output file or independently owned pool.
 """
 def __init__(self,root):self.root=Path(root).resolve()

 def materialize(self,data):
  allowed={'entity_id','position_m','resolution_m','section_mm','extent_mm','seed','cohort','scenario_id','zoom_bounds_m'}
  if not isinstance(data,dict) or set(data)-allowed:raise PatchRequestError(400,'invalid_request','Unknown vascular query fields or invalid JSON object')
  try:request_size=len(json.dumps(data,allow_nan=False).encode())
  except (TypeError,ValueError,OverflowError):raise PatchRequestError(400,'invalid_request','Query must be finite JSON') from None
  if request_size>LIMITS['max_request_bytes']:raise PatchRequestError(413,'request_limit','Vascular query exceeds request byte limit')
  entity=data.get('entity_id')
  if not isinstance(entity,str) or entity not in REGIONS:raise PatchRequestError(422,'unsupported_region','Select left/right vastus lateralis or kidney')
  kidney=entity in ('body-bp3d-FJ3145','body-bp3d-FJ3147')
  if kidney and any(k in data for k in ('cohort','section_mm','extent_mm')):raise PatchRequestError(400,'invalid_request','Muscle density/extent conditioning is not accepted for kidney')
  name,default_position=REGIONS[entity];position=_vector(data.get('position_m',default_position),'position_m')
  resolution=_number(data.get('resolution_m',50e-6),'resolution_m',LIMITS['min_resolution_m'],LIMITS['max_resolution_m'])
  section=data.get('section_mm',[.2,.2])
  if not isinstance(section,(list,tuple)) or len(section)!=2:raise PatchRequestError(400,'invalid_request','section_mm needs two dimensions')
  section=[_number(x,'section_mm',.01,1.) for x in section]
  extent=_number(data.get('extent_mm',.5),'extent_mm',.01,LIMITS['max_extent_mm'])
  seed=data.get('seed',19)
  if type(seed) is not int or not 0<=seed<=2**32-1:raise PatchRequestError(400,'invalid_request','seed must be a uint32 integer')
  cohort=data.get('cohort','young_men')
  if cohort not in ('young_men','young_women','older_men','older_women'):raise PatchRequestError(400,'invalid_request','Unknown human muscle cohort')
  scenario_id=data.get('scenario_id','human_kidney_control_v1' if kidney else 'human_quadriceps_1988_corrected_v1')
  if not isinstance(scenario_id,str):raise PatchRequestError(400,'invalid_request','scenario_id must be a string')
  if kidney!=scenario_id.startswith('human_kidney_'):raise PatchRequestError(400,'invalid_request','Scenario does not belong to selected organ')
  scenario=_scenario(scenario_id);zoom=data.get('zoom_bounds_m')
  if zoom is not None:
   if not isinstance(zoom,(list,tuple)) or len(zoom)!=2:raise PatchRequestError(400,'invalid_request','zoom_bounds_m needs two corners')
   zoom=[_vector(x,'zoom_bounds_m') for x in zoom]
   if any(a>=b for a,b in zip(*zoom)):raise PatchRequestError(400,'invalid_request','Zoom extents must be positive')
  if not _QUERY_LOCK.acquire(blocking=False):raise PatchRequestError(503,'busy','Another bounded vascular query is running; retry after it completes')
  try:
   if kidney:
    result=materialize_registered_kidney(self.root,entity_id=entity,position_m=position,resolution_m=resolution,seed=seed,scenario=deepcopy(scenario['hydraulic_scenario']),zoom_bounds_m=zoom)
    if result['graph']['realized']['filtration_capillary_edges']>LIMITS['max_kidney_filtration_edges']:raise PatchRequestError(413,'geometry_limit','Kidney filtration edge budget exceeded')
   else:
    registry=load_registry(self.root);density=registry['organs']['Muscle']['cohorts'][cohort]['density']['mean']
    count=int(math.floor(density*section[0]*section[1]+.5));edge_bound=max(1,5*count-4)
    diagonal=math.sqrt((1.4*extent)**2+sum(x*x for x in section))*1e-3
    sample_bound=edge_bound*(math.ceil(diagonal/resolution)+1)
    if not 1<=count<=LIMITS['max_capillaries'] or edge_bound>LIMITS['max_edges'] or sample_bound>LIMITS['max_zoom_samples']:
     raise PatchRequestError(413,'geometry_limit','Requested patch or conservative sampling bound exceeds local geometry limits')
    result=materialize_patch(self.root,entity_id=entity,position_m=position,resolution_m=resolution,seed=seed,cohort=cohort,
        section_mm=section,extent_mm=extent,hydraulic_scenario=deepcopy(scenario['hydraulic_scenario']),zoom_bounds_m=zoom)
   zoom_edges=result['zoom_edges'];edges=result['graph']['edges']
   if len(edges)>LIMITS['max_edges'] or len(zoom_edges)>LIMITS['max_edges'] or sum(len(e.get('centerline_samples_m',[])) for e in zoom_edges)>LIMITS['max_zoom_samples']:
    raise PatchRequestError(413,'geometry_limit','Materialized result exceeds geometry budget')
   response={'schema':'ihm.microvascular-patch-response.v1','patch':result,'scenario':scenario,'limits':dict(LIMITS),
     'selection':{'entity_id':entity,'name':name,'position_m':position,'position_basis':'caller-selected canonical point' if 'position_m' in data else 'explicit retained-atlas fixture point; inside and clearance reverified per request'},
     'read_only':True,'native_commands_issued':0,'response_units':{'position':'m','radius':'m','pressure':'Pa','flow':'m3/s'}}
   if len(json.dumps(response,allow_nan=False).encode())>LIMITS['max_response_bytes']:
    raise PatchRequestError(413,'response_limit','Serialized vascular result exceeds response budget')
   return response
  except PatchRequestError:raise
  except (KeyError,TypeError):raise PatchRequestError(409,'evidence_schema','Retained vascular evidence or materializer result has an unsupported schema') from None
  except (OSError,json.JSONDecodeError):raise PatchRequestError(503,'evidence_unavailable','Required retained vascular evidence is unavailable or unreadable') from None
  except ValueError as error:
   message=str(error)
   if any(x in message.lower() for x in ('changed','escaping','source','non-si')):raise PatchRequestError(409,'evidence_integrity','Retained vascular source identity or units failed verification') from None
   if any(x in message.lower() for x in ('limit exceeded','sampling exceeds')):raise PatchRequestError(413,'geometry_limit',message) from None
   raise PatchRequestError(422,'patch_not_supported',message) from None
  finally:_QUERY_LOCK.release()
