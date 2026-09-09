#!/usr/bin/env python3
"""Build named somatic priors without changing canonical anatomy identity."""
from pathlib import Path
import gzip
import hashlib
import json
import shutil
import numpy as np
from enrich_peripheral_routes import enrich
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/derived/canonical'
SOURCES={
 'touch':{'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6666033/','finding':'Human fingertip afferent recordings relate force stimuli to discharge and timing; present pressure gain is an uncalibrated reduction.'},
 'thermal':{'url':'https://pubmed.ncbi.nlm.nih.gov/3225599/','finding':'Human warming/cooling conduction estimates: means 0.5 and 2.1 m/s; transfer to all patches is a prior.'},
 'spindle':{'url':'https://pubmed.ncbi.nlm.nih.gov/2141632/','finding':'Human finger extensor spindle afferents respond to stretch; reduction omits fusimotor drive and detailed dynamics.'},
 'median':{'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC7580294/','finding':'Cadaver study maps median motor branches; supports named forearm relationships and documents variation.'},
 'arm':{'url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4079987/','finding':'Cadaver study of musculocutaneous/anterior-arm innervation with anatomical variation.'},
 'leg':{'url':'https://pubmed.ncbi.nlm.nih.gov/1766452/','finding':'Human stimulation maps lower-limb segmental motor supply; roots overlap broadly, vary and may be asymmetric.'},
 'femoral':{'url':'https://pubmed.ncbi.nlm.nih.gov/2254163/','finding':'Cadaver dissections trace vastus medialis branches to femoral/lumbar supply; coarse segment grouping here is inferred.'}}
# Explicit anatomical name rules, not nearest-surface innervation. Split/ambiguous
# muscles remain unsupported unless the source model identifies the compartment.
RULES=[
 ('musculocutaneous','cervical',('biceps brachii','brachialis','coracobrachialis')),
 ('axillary','cervical',('deltoid','teres minor')),
 ('suprascapular','cervical',('supraspinatus','infraspinatus')),
 ('radial','cervical',('triceps brachii','brachioradialis','anconeus','supinator','extensor carpi','extensor digitorum','extensor digiti minimi','extensor indicis','extensor pollicis','abductor pollicis longus')),
 ('median','cervical',('pronator','flexor carpi radialis','palmaris longus','flexor digitorum superficialis','flexor pollicis longus','abductor pollicis brevis','opponens pollicis')),
 ('ulnar','cervical',('flexor carpi ulnaris','adductor pollicis','interossei of hand','abductor digiti minimi of hand','flexor digiti minimi brevis of hand','opponens digiti minimi of hand')),
 ('femoral','lumbar',('rectus femoris','vastus','sartorius','iliacus')),
 ('obturator','lumbar',('adductor brevis','adductor longus','gracilis','obturator externus')),
 ('superior_gluteal','sacral',('gluteus medius','gluteus minimus','tensor fasciae latae')),
 ('inferior_gluteal','sacral',('gluteus maximus',)),
 ('sciatic_tibial','sacral',('long head of biceps femoris','semimembranosus','semitendinosus')),
 ('sciatic_fibular','sacral',('short head of biceps femoris',)),
 ('deep_fibular','sacral',('tibialis anterior','extensor digitorum longus','extensor hallucis','fibularis tertius')),
 ('superficial_fibular','sacral',('fibularis longus','fibularis brevis')),
 ('tibial','sacral',('gastrocnemius','soleus','plantaris','popliteus','tibialis posterior','flexor digitorum longus','flexor hallucis longus')),
 ('medial_plantar','sacral',('abductor hallucis','flexor digitorum brevis','first lumbrical of foot')),
 ('lateral_plantar','sacral',('adductor hallucis','flexor accessorius','plantar interosseous','abductor digiti minimi of foot','flexor digiti minimi brevis of foot','opponens digiti minimi of foot','second lumbrical of foot','third lumbrical of foot','fourth lumbrical of foot')),
 ('oculomotor','cranial',('medial rectus','superior rectus','inferior rectus','inferior oblique','levator palpebrae')),
 ('trochlear','cranial',('superior oblique',)),('abducens','cranial',('lateral rectus',)),
 ('accessory','cervical',('trapezius','sternocleidomastoid')),
 ('hypoglossal','cranial',('genioglossus','hyoglossus')),
 ('facial','cranial',('platysma','stylohyoid')),
]
# Region centers are anatomical display priors, projected to the real skin mesh.
PATCHES=[('palm','median','cervical',[.31,.00,.035]),('forearm','musculocutaneous','cervical',[.28,.21,.035]),
 ('upper_arm','axillary','cervical',[.22,.43,.04]),('thigh','lateral_femoral_cutaneous','lumbar',[.12,-.20,.04]),
 ('calf','sural','sacral',[.10,-.58,-.06]),('foot','medial_plantar','sacral',[.11,-.77,.06]),
 ('face','trigeminal','cranial',[.055,.71,.07]),('trunk','intercostal','thoracic',[.12,.34,.11])]

def build():
 anatomy=json.loads((OUT/'anatomy.json').read_text());mechanics=json.loads((OUT/'mechanics.json').read_text());brain=json.loads((OUT/'brain.json').read_text())
 entities={e['id']:e for e in anatomy['entities']};bn={n['id']:n for n in brain['nodes']}
 skin=next(e for e in anatomy['entities'] if e['role']=='skin')
 geo=json.loads(gzip.decompress((ROOT/skin['reference_geometry']['path']).read_bytes()))
 vertices=np.asarray(geo['positions'],float).reshape(-1,3)
 relays=[];nerves={};lines=[];bindings=[];unsupported=[];patches=[]
 ys={'cranial':.65,'cervical':.53,'thoracic':.30,'lumbar':.04,'sacral':-.08}
 def nerve(side,name,level):
  nid=f'peripheral-nerve-{side}-{name}';rid=f'peripheral-relay-{side}-{level}'
  if nid not in nerves:
   nerves[nid]={'id':nid,'name':f'{side} {name.replace("_"," ")} nerve','side':side,'relay_id':rid,'evidence_kind':'named_anatomical_prior','geometry_kind':'schematic_route','measured_axon_geometry':False}
  return nid,rid
 for side,sign in [('left',1),('right',-1)]:
  for level,y in ys.items():relays.append({'id':f'peripheral-relay-{side}-{level}','name':f'{side} {level} somatic relay','side':side,'position_m':[sign*.012,y,-.04],'kind':'cranial_sensory_motor_relay' if level=='cranial' else 'spinal_segment_group','evidence_kind':'regional_group_prior'})
 rp={r['id']:r['position_m'] for r in relays}
 for m in mechanics['muscles']:
  e=entities[m['canonical_entity_id']];name=e['name'].lower()
  side='left' if 'left ' in name else 'right' if 'right ' in name else None
  if side is None:unsupported.append({'muscle_id':m['id'],'reason':'No explicit laterality; no geometric innervation guess'});continue
  # Long digital extensors belong to the leg, before generic hand rule.
  candidates=sorted(RULES,key=lambda r:0 if r[0]=='deep_fibular' else 1)
  rule=next((r for r in candidates if any(token in name for token in r[2])),None)
  if 'longus colli' in name:rule=None  # neck muscle is not an extraocular oblique
  if 'adductor magnus' in name:
   rule=('sciatic_tibial','sacral',()) if 'Isch' in m.get('source_name','') else ('obturator','lumbar',()) if m.get('source_name') else None
  if not rule:unsupported.append({'muscle_id':m['id'],'reason':'Named innervation rule absent or mixed/ambiguous compartment'});continue
  nid,rid=nerve(side,*rule[:2]);pos=e['centroid_m'];length=max(.03,float(np.linalg.norm(np.array(pos)-rp[rid]))*1.15)
  target=f'brain-{"rh" if side=="left" else "lh"}-postcentral';motor=target.replace('postcentral','precentral')
  binding={'muscle_id':m['id'],'canonical_entity_id':e['id'],'name':e['name'],'side':side,'nerve_id':nid,'relay_id':rid,'brain_target_id':target,'brain_motor_id':motor,'path_length_m':length,'motor_delay_s':length/50.+.012,'afferent_delay_s':length/60.+.012,'rest_path_length_m':m['rest_path_length_m'],'max_isometric_force_n':m['max_isometric_force_n'],'fiber_axis':next(x['fiber_axis'] for x in mechanics['entities'] if x['id']==e['id']),'anchors':m['anchors'],'evidence_kind':'named_anatomical_prior','sensorimotor_region':'upper_limb' if 'biceps brachii' in name else 'lower_limb' if 'tibialis anterior' in name else None,'somatic_readout_gain':.025 if ('biceps brachii' in name or 'tibialis anterior' in name) else 0.}
  bindings.append(binding)
  lines.append({'id':f'peripheral-route-{m["id"]}','nerve_id':nid,'entity_ids':[e['id']],'points_m':[bn[motor]['position_m'],rp[rid],pos],'kind':'motor_and_proprioceptive','evidence_kind':'schematic_anatomical_prior'})
 for side,sign in [('left',1),('right',-1)]:
  for name,nn,level,center in PATCHES:
   center=np.array(center)*[sign,1,1];pos=vertices[np.argmin(np.linalg.norm(vertices-center,axis=1))].tolist();nid,rid=nerve(side,nn,level)
   length=max(.03,float(np.linalg.norm(np.array(pos)-rp[rid]))*1.15);target=f'brain-{"rh" if side=="left" else "lh"}-postcentral'
   p={'id':f'peripheral-skin-{side}-{name}','name':f'{side} {name.replace("_"," ")} receptor patch','side':side,'body_entity_id':skin['id'],'position_m':pos,'radius_m':.018,'nerve_id':nid,'relay_id':rid,'brain_target_id':target,'path_length_m':length,'evidence_kind':'regional_cutaneous_prior','surface_anchor':'nearest skin vertex to authored regional landmark; not measured receptor distribution','sensorimotor_region':'upper_limb' if name in ['palm','forearm','upper_arm'] else 'lower_limb' if name in ['thigh','calf','foot'] else None,'modalities':['pressure_pa','temperature_C','stretch_fraction']}
   patches.append(p);lines.append({'id':f'peripheral-route-{p["id"]}','nerve_id':nid,'entity_ids':[skin['id']],'points_m':[pos,rp[rid],bn[target]['position_m']],'kind':'cutaneous_afferent','evidence_kind':'schematic_anatomical_prior'})
 preserved=OUT/'peripheral-sources';preserved.mkdir(exist_ok=True);receipts=[]
 for rel in ['ibm/processes/transduction.py','ibm/processes/effector.py','ibm/runtime/step.py']:
  src=ROOT.parent/'IBM-1'/rel;dst=preserved/('ibm-'+Path(rel).parent.name+'-'+Path(rel).name);shutil.copyfile(src,dst)
  receipts.append({'source_path':str(src),'preserved_path':str(dst.relative_to(ROOT)),'sha256':hashlib.sha256(dst.read_bytes()).hexdigest(),'role':'inspected_reference_not_executed'})
 data={'schema_version':1,'id':'ihm-body-peripheral','frame':anatomy['frame'],'units':{'position':'m','time':'s','activity':'Hz','pressure':'Pa','temperature':'C','activation':'1'},'relays':relays,'nerves':list(nerves.values()),'receptor_patches':patches,'muscle_bindings':bindings,'unsupported_muscles':unsupported,'sources':SOURCES,'source_receipts':receipts,'parameters':{'pressure_gain_hz_pa':.005,'stretch_gain_hz':200.,'temperature_gain_hz_C':8.,'baseline_skin_temperature_C':32.,'max_receptor_rate_hz':200.,'receptor_tau_s':.02,'activation_tau_s':.03,'tactile_velocity_m_s':50.,'warm_velocity_m_s':.5,'cold_velocity_m_s':2.1,'central_afferent_delay_s':.012},'parameter_scope':'Uncalibrated illustrative transfer priors; velocity anchors differ by fiber class. Rate is excess evoked activity above omitted spontaneous baseline.','biological_validation':False,'ibm_reuse_scope':'Existing BodyBrain executes preserved IBM Wilson-Cowan/shunting functions. Peripheral runtime is a new causal time-domain reduction informed by IBM receptor/effector declarations; IBM spectral runtime, full multimodal materialization, uncertainty, inference and learned motor policies are not integrated.','counts':{'nerves':len(nerves),'relays':len(relays),'receptor_patches':len(patches),'muscle_bindings':len(bindings),'unsupported_muscles':len(unsupported)},'limitations':['Schematic centerlines are inferred visualization; do not represent dissected or measured nerves.','Coarse contralateral postcentral routing; modality-specific thalamic nuclei and decussation sites unresolved.','Named muscle bindings are generic anatomical priors; mixed innervation and absent name rules remain unsupported.','No autonomic controller or respiration override; somatic commands only.']}
 enrich(data,lines,ROOT)
 relays=data['relays']
 for fn,payload in [('peripheral.json',data),('peripheral_display.json',{'schema_version':1,'frame':anatomy['frame'],'lines':lines,'nodes':relays+patches,'evidence_kind':'schematic_anatomical_prior'})]:
  (OUT/fn).write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
 print(json.dumps(data['counts']))
if __name__=='__main__':build()
