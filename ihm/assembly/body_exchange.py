"""Native-owned tissue exchange and conservative spatial views.

There is no Python blood, lymph or solute integration here. The native solver
owns oncotic/osmotic feedback, transport, metabolism and drainage. Fine regions
partition an observed native quantity; they never add a second reservoir.
"""
from copy import deepcopy
import hashlib,json,math
from pathlib import Path

ORGANS=('Fat','Bone','Brain','Gut','LeftKidney','RightKidney','Liver','LeftLung','RightLung','Muscle','Myocardium','Skin','Spleen')
SUBSTANCES=('Albumin','Glucose','Oxygen','CarbonDioxide','Sodium','Potassium','Chloride')
POOLS=('vascular','extracellular','intracellular')

def number(value,label,nonnegative=False):
 if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or (nonnegative and value<0):raise ValueError('Invalid native quantity: '+label)
 return float(value)

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

class NativeTissueExchange:
 def __init__(self,organs=ORGANS,partitions=None,native_identity=None,source_hashes=None):
  self.organs=tuple(organs)
  if not self.organs or len(set(self.organs))!=len(self.organs) or any(o not in ORGANS for o in self.organs):raise ValueError('Unknown or duplicate native organ')
  self.partitions=deepcopy(partitions or {});self.native_identity=deepcopy(native_identity or {});self.source_hashes=dict(source_hashes or {})
  if any(o not in self.organs for o in self.partitions):raise ValueError('Partition has no selected native owner')
  for organ,rows in self.partitions.items():
   if len({r['id'] for r in rows})!=len(rows) or any(r['id']=='unresolved_remainder' for r in rows):raise ValueError('Duplicate or reserved partition ID')
   weights=[number(r['fraction'],'partition fraction',True) for r in rows]
   if math.fsum(weights)>1.:raise ValueError('Spatial partition exceeds native owner')

 @classmethod
 def from_workspace(cls,root,reference_snapshot,native_identity,organs=ORGANS):
  """Allocate synthetic cutaneous graph lumen against an explicit native reference.

  Lumen/reference vascular volume is a geometric allocation assumption, not a
  measured mapping of terminal capillaries into native lumped blood regions.
  The same fraction is transferred to interstitial/intracellular spatial views.
  """
  root=Path(root).resolve();micro_path=root/'data/derived/canonical/microvascular.json';anatomy_path=root/'data/derived/canonical/anatomy.json'
  micro=json.loads(micro_path.read_text());anatomy=json.loads(anatomy_path.read_text());entities={e['id']:e for e in anatomy['entities']};rows=[]
  reference_volume=number(reference_snapshot['values'].get('tissue.Skin.vascular.volume_ml'),'reference Skin vascular volume',True)
  if reference_volume<=0:raise ValueError('Positive reference vascular volume needed for allocation')
  hashes={str(p.relative_to(root)):digest(p) for p in [micro_path,anatomy_path,Path(__file__),root/'scripts/native_tissue_ports.h',root/'scripts/native_tissue_compression.h']}
  for unit in micro['units']:
   if 'forearm skin' not in unit['territory']:raise ValueError('Unregistered microvascular territory')
   attachment=unit['material_attachment'];entity=entities[attachment['entity_id']]
   geometry=entity.get('reference_geometry',{});path=geometry.get('path') or entity.get('geometry_path');expected=geometry.get('sha256') or entity.get('geometry_sha256')
   if not path or not expected:raise ValueError('Missing canonical surface geometry receipt')
   source=(root/path).resolve()
   if not source.is_relative_to(root) or digest(source)!=expected:raise ValueError('Changed microvascular attachment geometry')
   hashes[str(source.relative_to(root))]=expected
   positions=unit['positions_m'];edges=unit['edges'];radii=unit['radius_m']
   if len(edges)!=len(radii):raise ValueError('Malformed microvascular graph')
   volumes=[]
   for (a,b),radius in zip(edges,radii):
    r=number(radius,'lumen radius',True)
    if r<=0 or not 0<=a<len(positions) or not 0<=b<len(positions):raise ValueError('Invalid microvascular edge')
    length=math.dist(positions[a],positions[b]);number(length,'edge length',True)
    if length<=0:raise ValueError('Zero-length microvascular edge')
    volumes.append(math.pi*r*r*length)
   lumen=math.fsum(volumes)
   rows.append({'id':unit['id'],'fraction':lumen*1e6/reference_volume,'territory':unit['territory'],'material_attachment':deepcopy(attachment),'source_geometry_sha256':expected,'geometry':{'positions_m':deepcopy(positions),'edges':deepcopy(edges),'radius_m':list(radii),'edge_volume_fractions':[v/lumen for v in volumes]},'evidence_kind':'synthetic_graph_with_inferred_native_territory_allocation','reference_lumen_volume_m3':lumen,'reference_native_vascular_volume_ml':reference_volume,'reference_time_s':reference_snapshot['time_s']})
  identity=deepcopy(native_identity)
  if not identity:raise ValueError('Explicit native identity receipt is required')
  from .body_microstructure import anatomical_supports
  raw_supports=anatomical_supports(anatomy,organs)
  supports={organ:[{'id':'native-anatomical-support-'+organ,'fraction':1.,'territory':organ,'material_attachment':None,'anatomical_supports':regions,'evidence_kind':'grouped_source_anatomy_allocation' if any(r['material_attachment'] for r in regions) else 'unlocalized_native_owner','allocation_measure':'source surface-area subdivisions expanded only on demand'}] for organ,regions in raw_supports.items()}
  if 'Skin' in organs:
   fraction=math.fsum(r['fraction'] for r in rows)
   if fraction>1:raise ValueError('Fine skin graph exceeds native reference vascular volume')
   for region in supports['Skin']:region['fraction']*=1.-fraction
   supports['Skin']=rows+supports['Skin']
  hashes['ihm/assembly/body_microstructure.py']=digest(root/'ihm/assembly/body_microstructure.py')
  return cls(organs,supports,identity,hashes)

 def observe(self,snapshot):
  time=number(snapshot.get('time_s'),'time',True);values=snapshot.get('values',{});compartments={}
  owners=[o+'.'+p for o in self.organs for p in POOLS]+['Lymph','VenaCava']
  for owner in owners:
   prefix='tissue.'+owner;volume=number(values.get(prefix+'.volume_ml'),prefix+'.volume_ml',True);mass={};concentration={};consistency={};ionic={};partial={}
   for sub in SUBSTANCES:
    key=prefix+'.'+sub;m=values.get(key+'.mass_g');mass[sub]=None if m is None else number(m,key+'.mass_g',True)
    c=values.get(key+'.concentration_g_per_l');concentration[sub]=None if c is None else number(c,key+'.concentration_g_per_l',True)
    consistency[sub]=None if c is None or mass[sub] is None else mass[sub]-c*volume/1000.
    if sub in ('Sodium','Potassium','Chloride'):ionic[sub]=values.get(key+'.molarity_mmol_per_l')
    if sub in ('Oxygen','CarbonDioxide'):partial[sub]=values.get(key+'.partial_pressure_mmhg')
   pressure=values.get(prefix+'.pressure_mmhg')
   if pressure is not None:number(pressure,prefix+'.pressure_mmhg')
   for label,v in [*ionic.items(),*partial.items()]:
    if v is not None:number(v,label,True)
   compartments[owner]={'native_owner':owner,'volume_ml':volume,'pressure_mmhg':pressure,'mass_g':mass,'concentration_g_per_l':concentration,'mass_concentration_residual_g':consistency,'ionic_molarity_mmol_per_l':ionic,'gas_partial_pressure_mmhg':partial}
  transfers=[];rates={o:0. for o in owners}
  def transfer(path,source,target):
   q=number(values.get('tissue.path.'+path+'.flow_ml_per_s'),path+' flow');transfers.append({'native_path':path,'source':source,'target':target,'flow_ml_per_s':q,'kind':'observed_native_fluid_flow'})
   if source in rates:rates[source]-=q
   if target in rates:rates[target]+=q
  for organ in self.organs:
   transfer(organ+'E1To'+organ+'E2',organ+'.vascular',organ+'.extracellular')
   transfer(organ+'E3To'+organ+'I',organ+'.extracellular',organ+'.intracellular')
   transfer(organ+'E3To'+organ+'L1',organ+'.extracellular','Lymph')
  transfer('LymphToVenaCava','Lymph','VenaCava')
  # This is a selected incidence ledger, not dV/dt of the whole native system:
  # vascular perfusion, unselected organs, sweat, renal output and GI input exist.
  regions=[]
  for organ in self.organs:
   rows=self.partitions.get(organ,[])
   for pool in POOLS:
    owner=organ+'.'+pool;parent=compartments[owner];used_volume=[];used_mass={s:[] for s in SUBSTANCES}
    for row in rows:
     fraction=row['fraction'];volume=parent['volume_ml']*fraction;mass={s:None if parent['mass_g'][s] is None else parent['mass_g'][s]*fraction for s in SUBSTANCES};used_volume.append(volume)
     for sub in SUBSTANCES:
      if mass[sub] is not None:used_mass[sub].append(mass[sub])
     regions.append({'id':row['id']+'.'+pool,'owner':owner,'fraction':fraction,'volume_ml':volume,'mass_g':mass,'independent_store':False,'material_attachment':deepcopy(row.get('material_attachment')),'evidence_kind':row.get('evidence_kind'),'allocation_measure':row.get('allocation_measure','reference geometric lumen fraction')})
    regions.append({'id':organ+'.unresolved_remainder.'+pool,'owner':owner,'fraction':1.-math.fsum(r['fraction'] for r in rows),'volume_ml':parent['volume_ml']-math.fsum(used_volume),'mass_g':{s:None if parent['mass_g'][s] is None else parent['mass_g'][s]-math.fsum(used_mass[s]) for s in SUBSTANCES},'independent_store':False,'material_attachment':None})
  return {'schema':'native_tissue_exchange_v1','time_s':time,'native_identity':deepcopy(self.native_identity),'source_hashes':dict(self.source_hashes),'native_compartments':compartments,'native_circuit':{k:v for k,v in values.items() if k.startswith(('tissue.path.','tissue.node.','tissue.compression.'))},'partitions':regions,'fluid_transfers':transfers,'internal_volume_rate_ml_per_s':rates,'solute_fluxes':[],'whole_body_mass_closure_claimed':False,'thermal_boundary_c':{'core':values.get('tissue.core_temperature_c'),'skin':values.get('tissue.skin_temperature_c')},'feedback_owner':'Native Tissue oncotic/osmotic laws, cardiovascular fluid solve and Diffusion albumin/solute transport','limitations':['Spatial regions are conservative views, not new blood/lymph stores.','Snapshot flow incidence is not integrated mass transfer or complete dV/dt.','Gas masses are named native free-substance pools, not total hemoglobin-bound oxygen/carbon.','Native Skin is one lumped compartment; source graph pressures and terminal vessel correspondence are unresolved.','No solute flux is reconstructed from post-step concentration: native reactions, capping and transport order matter.']}

 def project_networks(self,observation):
  """Return synthetic graph states as subdivisions of existing observed partitions."""
  owners={r['id']:r for r in observation['partitions']};networks=[]
  for organ,rows in self.partitions.items():
   for row in rows:
    if 'geometry' not in row:continue
    parent=owners[row['id']+'.vascular'];geometry=deepcopy(row['geometry']);weights=geometry['edge_volume_fractions'];volumes=[parent['volume_ml']*w for w in weights];volumes[-1]=parent['volume_ml']-math.fsum(volumes[:-1]);mass={}
    for sub,total in parent['mass_g'].items():
     mass[sub]=None if total is None else [total*w for w in weights]
     if mass[sub] is not None:mass[sub][-1]=total-math.fsum(mass[sub][:-1])
    networks.append({'id':row['id'],'owner':parent['owner'],'time_s':observation['time_s'],'geometry':geometry,'edge_volume_ml':volumes,'edge_mass_g':mass,'independent_store':False,'pressure_resolution':'only native parent pressure observed; no inferred microvascular gradient','native_pressure_mmhg':observation['native_compartments'][parent['owner']]['pressure_mmhg']})
  return networks

 def project_anatomical_supports(self,observation):
  """Expand source-surface quantities only when the caller requests spatial detail."""
  owners={r['id']:r for r in observation['partitions']};result=[]
  for organ,rows in self.partitions.items():
   for row in rows:
    supports=row.get('anatomical_supports')
    if not supports:continue
    for pool in POOLS:
     parent=owners[row['id']+'.'+pool];used=[];used_mass={s:[] for s in SUBSTANCES}
     for i,support in enumerate(supports):
      last=i==len(supports)-1;fraction=support['fraction'];volume=parent['volume_ml']-math.fsum(used) if last else parent['volume_ml']*fraction;mass={s:None if m is None else m-math.fsum(used_mass[s]) if last else m*fraction for s,m in parent['mass_g'].items()};used.append(volume)
      for sub,m in mass.items():
       if m is not None:used_mass[sub].append(m)
      result.append({'id':support['id']+'.'+pool,'owner':parent['owner'],'volume_ml':volume,'mass_g':mass,'independent_store':False,'material_attachment':deepcopy(support['material_attachment']),'source_geometry':deepcopy(support.get('source_geometry')),'evidence_kind':support['evidence_kind'],'allocation_measure':support['allocation_measure']})
  return result
