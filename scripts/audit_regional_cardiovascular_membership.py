"""Source-bound membership audit; no region guesses, native runs, or shared edits."""
import hashlib,json,re,xml.etree.ElementTree as ET
from pathlib import Path
from patch_cardiovascular_region_io import ROOT,topology
HEADER=ROOT/'scripts/native_regional_skin.h'
HEADER_SHA='1948643357d0d0acdfab82582af886e6c727c7938d56a65965d8674cc157dee6'
PARENT=ROOT/'data/runtime/physiology/variants/whole_body_integrity_regional_skin_graph_v2'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def audit():
 if sha(HEADER)!=HEADER_SHA:raise ValueError('Changed regional implementation; repeat ownership audit')
 header=HEADER.read_text();block=header.split('static inline const std::array<std::string,9> path_names{',1)[1].split('};',1)[0];names=re.findall(r'"([^"]+)"',block);assert len(names)==9
 assigned,_=topology();reference=json.loads((ROOT/'data/derived/systemic/exertion_v3/exercise/native/manifest.json').read_text());state=Path(reference['configuration']['state_path']);assert sha(state)==reference['state_sha256'];tree=ET.fromstring(state.read_bytes());paths={}
 for node in tree.iter():
  if node.tag.split('}')[-1]=='FluidPath':
   row={x.tag.split('}')[-1]:x.text for x in node};paths[row['Name']]=row
 rows=[]
 for name in names:
  p=paths[name];rows.append({'source_path':name,'source_node':p['SourceNode'],'target_node':p['TargetNode'],'source_constructor_region':assigned.get(name),'resistance_baseline':p.get('ResistanceBaseline'),'has_resistance_baseline':'ResistanceBaseline' in p,'regional_paths':['IHM_'+r+'_'+name for r in ('region_a','region_b','residual')],'direct_routine_vascular_tone_member':False})
 assert not (set(names)&set(assigned))
 assert {r['source_path'] for r in rows if r['has_resistance_baseline']}=={'SkinE1ToSkinE2','SkinL1ToSkinL2'}
 retained={name:assigned[name] for name in ('Aorta1ToSkin1','Skin1ToSkin2')};assert all(r['region']=='Extrasplanchnic' for r in retained.values())
 assert all(name not in names for name in retained)
 return {'schema':'ihm.regional-cardiovascular-membership.v1','header_sha256':HEADER_SHA,'regional_parent_manifest_sha256':sha(PARENT/'manifest.json'),'reference_state_sha256':sha(state),'cloned_paths':rows,'retained_vascular_tone_paths':retained,'new_external_pressure_paths':['IHM_'+r+'_ExternalPressureDrive' for r in ('region_a','region_b','residual')],'cloned_resistance_path_count':6,'cloned_direct_vascular_members':0,'source_conclusion':'No missing region label on these nine source laws: vascular Skin resistors remain unsplit. Filtration and lymph resistance are distinct controls. Do not assign Extrasplanchnic to tissue/lymph clones.','native_response_verified':False}
if __name__=='__main__':print(json.dumps(audit(),indent=2))
