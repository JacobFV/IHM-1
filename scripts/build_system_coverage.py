"""Join declared physiology, executed native mechanisms and actual observation coverage."""
import json
from collections import Counter
from pathlib import Path
from ihm import body
from ihm.forge.acquisition import sha256

# Group membership is an explicit application ontology, not a geometric registration.
DOMAINS=[
('cardiovascular','Heart and vascular circulation',['cardiac','vascular','venous','capillary'],['Cardiovascular'],['Aorta','VenaCava','Heart','Myocardium'],['arterial','venous','cardiac']),
('blood','Blood composition and chemistry',['blood'],['BloodChemistry'],[],[]),
('respiratory','Respiration and gas exchange',['respiratory','airway','alveolar'],['Respiratory'],['Lung','Bronchi','Alveoli','Pleural'],['respiratory']),
('renal','Renal and urinary',['renal','urinary'],['Renal'],['Kidney','Nephron','Glomer','Tubules','Bladder','Ureter'],['urinary']),
('interstitial','Interstitial fluid and cell exchange',['interstitial','extracellular'],['Tissue'],['Extracellular','Intracellular'],[]),
('lymphatic','Lymph transport and nodes',['lymphatic','lymph_node'],['Tissue'],['Lymph'],['lymphatic']),
('integumentary','Skin, glands and follicles',['integumentary','dermal','sweat','follicular'],['Tissue','Energy'],['Skin'],['integumentary']),
('bioelectric','Cellular and epithelial bioelectricity',[],[],[],[]),
('neural','Central, peripheral and autonomic nervous',['neural','effector'],['Nervous'],['Brain'],['nervous']),
('endocrine','Endocrine and pancreatic regulation',['endocrine','pancreatic'],['Endocrine'],[],['endocrine']),
('metabolic','Metabolism and adipose',['metabolic','adipose'],['Energy','Tissue'],['Fat'],[]),
('digestive','Digestive transport and absorption',['digestive'],['Gastrointestinal'],['Gut','Intestine','Splanchnic','Stomach'],['digestive']),
('hepatic','Hepatic processing',['hepatic'],['Hepatic'],['Liver','Portal'],[]),
('immune','Immune and inflammation',['immune'],['BloodChemistry'],['Spleen'],[]),
('hematopoietic','Hematopoiesis and spleen',['hematopoietic','splenic'],['BloodChemistry'],['Spleen','Bone'],[]),
('musculoskeletal','Bone, muscle, joints and tendons',['bone','joint','tendon','mechanical','structural'],['Tissue'],['Bone','Muscle'],['skeletal','muscular','connective']),
('sensory','Sensory systems',['sensory'],['Nervous'],[],['sensory']),
('csf','Cerebrospinal fluid',['csf'],[],[],[]),
('reproductive','Reproductive and pregnancy-specific',['reproductive','uterine','placental'],[],[],['reproductive']),
('thermal','Thermoregulation and environment',['thermal'],['Energy','Environment'],['Internal','External','Clothing'],[]),
('devices','Devices and pharmacology',['device'],['Drug','AnesthesiaMachine','Inhaler','ElectroCardioGram'],[],[])]


def build(root=Path('.')):
    root=Path(root);out=root/'data/derived';registry=body()
    graph_path=out/'native-circuits/graph.json';g=json.loads(graph_path.read_text())
    manifest=json.loads((out/'app/manifest.json').read_text())
    prior=json.loads((out/'population/nhanes-2017-2018/joint-population-prior.json').read_text())
    components=registry.describe();measured=set(prior['components'])
    rows=[]
    for id,name,prefixes,native,compnames,geometry in DOMAINS:
        physical=[c for c in components if c['id'].split('.')[0] in prefixes]
        systems=[s['type'] for s in g['systems'] if any(s['type']=='BioGears'+key+'SystemData' or s['type']=='BioGears'+key+'Data' for key in native)]
        comps=[c['id'] for c in g['compartments'] if any(key.lower() in c['name'].lower() for key in compnames)]
        matched=[s for s in manifest['structures'] if s['system'] in geometry]
        meshes=[s['id'] for s in matched]
        geometry_by_family=dict(Counter(s['model_id'] for s in matched))
        observations=[c['id'] for c in physical if c['id'] in measured]
        limitations=[]
        if id in ('immune','hematopoietic','sensory','lymphatic'):limitations.append('Native presence is partial coverage; it does not establish all minor mechanisms or cell populations.')
        if id=='lymphatic':limitations.append('Native lumped lymph drainage is executed; subject-matched whole-body lymphatic geometry/flow is unavailable.')
        if id=='reproductive':limitations.append('Standalone published gonadotropin dynamics use prescribed ovarian inputs; no autonomous menstrual cycle, pregnancy dynamics or matching to male reference anatomy.')
        if id=='bioelectric':limitations.append('BETSE generic-tissue solver and human wound-field fit exist; human channel/pump kinetics are not identified.')
        if id=='csf':limitations.append('Published CSF model executes with optional native MAP one-way input; no matched ICP calibration, feedback to native physiology, glymphatic transport or posture-specific venous dynamics.')
        status='native_model_partial' if systems or comps else 'declared_not_executed'
        if id=='csf' and (out/'csf/index.json').is_file():status='published_csf_model_with_native_pressure_input'
        if id=='thermal' and (out/'thermal/index.json').is_file():
            status='native_energy_and_published_supine_thermoregulation'
            limitations.append('JOS-3 executes 85 thermal nodes and measured whole-body bedding boundaries separately from BioGears; no matched cross-engine heat/perfusion calibration.')
        if id=='bioelectric':status='tissue_solver_and_human_field_fit'
        if id=='reproductive' and (out/'reproductive/trajectory.json').is_file():status='published_standalone_reproductive_endocrine_model'
        if id=='musculoskeletal':status='native_tissue_and_opensim_mechanics'
        rows.append({'id':id,'name':name,'components':physical,'components_count':len(physical),'native_systems':systems,
                     'native_compartments':comps,'native_compartment_match':'literal name keyword; does not establish a coordinate transform',
                     'measured_components':observations,'geometry_count':len(meshes),'geometry_ids':meshes,'geometry_by_family':geometry_by_family,'geometry_families_are_independent':True,'status':status,'limitations':limitations})
    declared={c['id'] for row in rows for c in row['components']}
    missing=set(registry.components)-declared
    if missing:raise ValueError(f'Unclassified declared quantities: {sorted(missing)}')
    result={'schema_version':1,'systems':rows,'summary':{'system_domains':len(rows),'declared_components':len(components),
            'population_measured_components':len(measured),'native_circuit_nodes':g['summary']['unique_nodes'],
            'native_circuit_paths':g['summary']['unique_paths'],'native_compartments':len(g['compartments']),
            'native_properties':g['summary']['circuit_scalar_fields']+g['summary']['system_scalar_fields'],
            'all_interactions_calibrated':False},'native_source':g['source'],'native_graph_sha256':sha256(graph_path),
            'limitations':['Coverage means a mechanism or observation is represented, not complete biological validation.',
                          'Model families and reference anatomies do not describe one matched subject.',
                          'Measured cross-sectional population states do not identify temporal interaction coefficients.']}
    (out/'system-coverage.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['summary']));return result

if __name__=='__main__':build()
