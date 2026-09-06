"""Explicit anatomical support priors for native tissue-owner partitions.

Surface-area weights are display allocation proxies, not measured perfused
volumes, vessel correspondences, or an assertion that atlas meshes are disjoint.
"""
import math
from copy import deepcopy


def anatomical_supports(anatomy,organs):
 result={}
 for organ in organs:
  selected=[]
  for entity in anatomy['entities']:
   name=entity['name'].lower();role=entity.get('role');match=False
   if organ=='Skin':match=role=='skin'
   elif organ=='Bone':match=role=='rigid_bone'
   elif organ=='Muscle':match=role=='muscle'
   elif organ=='Brain':match=role=='nerve' and any(x in name for x in ('cerebell','cerebral hemisphere','thalamus','caudate nucleus','pons','medulla oblongata','midbrain','corpus callosum'))
   elif organ=='Gut':match=role=='soft_organ' and any(x in name for x in ('stomach','colon','duodenum','ileum','jejunum','cecum','rectum')) and not any(x in name for x in ('meso','duct'))
   elif organ=='Liver':match=role=='soft_organ' and ('hepatovenous segment' in name or name=='caudate lobe of liver')
   elif organ=='Myocardium':match=name in ('wall of ventricle','wall of left atrium','wall of right atrium')
   elif organ in ('LeftKidney','RightKidney'):match=name==('left kidney' if organ=='LeftKidney' else 'right kidney')
   elif organ in ('LeftLung','RightLung'):match=role=='soft_organ' and ('left lung' if organ=='LeftLung' else 'right lung') in name
   elif organ=='Spleen':match=name=='spleen'
   if match:
    area=entity.get('surface_area_m2');geometry=entity.get('reference_geometry')
    if not isinstance(area,(int,float)) or not math.isfinite(area) or area<=0 or not geometry:continue
    selected.append((entity,area))
  total=math.fsum(a for _,a in selected);rows=[]
  for entity,area in selected:
   rows.append({'id':'native-support-'+entity['id'],'fraction':area/total,'territory':organ,'material_attachment':{'entity_id':entity['id'],'kind':'whole_source_surface_support'},'source_geometry':deepcopy(entity['reference_geometry']),'evidence_kind':'anatomical_surface_with_engineering_quantity_allocation','allocation_measure':'normalized source surface area, not perfused tissue volume','source_name':entity['name']})
  if rows:
   rows[-1]['fraction']=1.-math.fsum(r['fraction'] for r in rows[:-1])
  else:
   rows=[{'id':'native-unlocalized-'+organ,'fraction':1.,'territory':organ,'material_attachment':None,'evidence_kind':'native_compartment_without_registered_anatomical_support','allocation_measure':'entire native owner retained unlocalized; no fabricated fine network'}]
  result[organ]=rows
 return result
