"""unify searchable provenance without merging incompatible biological models."""
from collections import Counter
import json
from pathlib import Path
import sqlite3
from ihm.forge.acquisition import sha256

ROOT=Path('data/derived')


def read(path): return json.loads(Path(path).read_text())
def rows(path):
    with Path(path).open() as stream:
        for line in stream:
            if line.strip(): yield json.loads(line)


def main():
    path=ROOT/'evidence-catalog.sqlite'; temporary=path.with_suffix('.partial.sqlite')
    temporary.unlink(missing_ok=True);db=sqlite3.connect(temporary)
    db.executescript('CREATE TABLE parameters(id INTEGER PRIMARY KEY,source TEXT,name TEXT,value TEXT,units TEXT,basis TEXT,source_file TEXT,reference TEXT,context TEXT,raw_json TEXT);CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);')
    counts=Counter();bases=Counter();inputs=[]
    def add(source,name,value,units,basis,file,reference,raw):
        context=raw.get('equation_context') or raw.get('source_label') or raw.get('calibration_context') or ''
        db.execute('INSERT INTO parameters(source,name,value,units,basis,source_file,reference,context,raw_json) VALUES (?,?,?,?,?,?,?,?,?)',
           (source,name,json.dumps(value,allow_nan=False),units,basis,file,reference,context,json.dumps(raw,allow_nan=False)))
        counts[source]+=1;bases[basis]+=1
    physiology=ROOT/'physiology'
    for filename,basis in [('coefficients.jsonl','upstream_model_candidate'),('reported_fitted_coefficients.jsonl','upstream_reported_fit'),('validation_targets.jsonl','literature_validation_target')]:
        p=physiology/filename;inputs.append(p)
        for row in rows(p):
            add(row.get('source','biogears'),row.get('name') or row.get('property') or row.get('target') or row.get('property_name') or filename,
                row.get('reference_value',row.get('value') if row.get('value') is not None else row.get('expression')),row.get('units'),basis,row.get('source_file',''),
                json.dumps(row.get('citation_keys',row.get('reference',''))),row)
    for p in sorted((ROOT/'anatomy').glob('opensim__*.json')):
        inputs.append(p);model=read(p)
        for category in ('bodies','muscles'):
            for entity in model[category]:
                for name,parameter in entity.get('parameters',{}).items():
                    value=parameter.get('value')
                    if value in ('',None): continue
                    raw={'model':model['name'],'entity':entity['name'],'category':category,
                         'name':name,**parameter,'sha256':model['sha256'],'revision':model['source_revision'],
                         'subject_specific_calibration':False}
                    add('opensim',model['name']+'/'+entity['name']+'/'+name,value,parameter.get('unit'),
                        'published_musculoskeletal_model_parameter',model['source_path'],str(model['publications']),raw)
    p=ROOT/'population/nhanes-2017-2018/population-priors.jsonl';inputs.append(p)
    for row in rows(p):
        add('nhanes-2017-2018',row['population']+'/'+row['component'],row['estimate']['mean'],row['unit'],
            'empirical_population_state_prior',row['source_file'],row['reference'],row)
    p=ROOT/'semantics/vascular-geometry-measurements.jsonl';inputs.append(p)
    for row in rows(p):
        add('hra-vccf',row['entity']+'/'+row['property']+'/'+row['sex'],row['value'],row['units'],
            'literature_geometry_summary',row['source_file'],row['reference'],row)
    p=ROOT/'integumentary/validation-targets.jsonl'
    if p.exists():
        inputs.append(p)
        for row in rows(p):
            add('human-wound-study-2011',row['name'],row['value'],row['units'],row['basis'],row['source_file'],row['reference'],row)
    # Executed native scalar fields retain state/parameter status; not calibration.
    p=ROOT/'native-circuits/graph.json'
    if p.exists():
        inputs.append(p);graph=read(p)
        for category in ('nodes','paths','systems'):
            for owner in graph[category]:
                for name,value in owner.get('properties',{}).items():
                    if not isinstance(value,dict):continue
                    raw={'owner':owner.get('id',owner.get('type')),'quantity':name,**value,'equation_context':category,
                         'source_sha256':graph['source']['sha256'],'independently_calibrated':False}
                    add('biogears-native',str(raw['owner'])+'/'+name,value.get('value'),value.get('unit'),
                        'native_initialized_scalar_state_or_parameter',graph['source']['path'],json.dumps(graph['source']),raw)
    p=ROOT/'reproductive/trajectory.json'
    if p.exists():
        inputs.append(p);model=read(p)
        for parameter in model['parameters']:
            add('schlosser-selgrade-2000',parameter['component']+'/'+parameter['id'],parameter['value'],parameter['unit'],
                'published_reproductive_model_parameter',str(p),json.dumps(model['source']),parameter)
    p=ROOT/'csf/baseline.json'
    if p.exists():
        inputs.append(p);model=read(p)
        for parameter in model['parameters']:
            add('ursino-lodi-1997',parameter['id'],parameter['value'],parameter['unit'],
                'published_csf_model_parameter',str(p),json.dumps(model['source']),parameter)
    p=ROOT/'thermal/index.json'
    if p.exists():
        inputs.append(p);index=read(p)
        for run in index['runs']:
            cp=Path(run['coefficient_path']);inputs.append(cp);model=read(cp)
            source={'revision':index['source']['revision'],'coefficient_sha256':sha256(cp),'profile':run['id'],'independently_calibrated':False}
            for name,unit in [('capacity_J_K','J/K'),('body_surface_area_m2','m2')]:
                for i,value in enumerate(model[name]):
                    add('jos3',run['id']+'/'+name+'/'+str(i),value,unit,'published_thermal_model_parameter',str(cp),index['source']['url'],source)
            for i,row in enumerate(model['conductance_W_K']):
                for j,value in enumerate(row):
                    if value:
                        add('jos3',run['id']+f'/conductance/{i}/{j}',value,'W/K','published_thermal_model_parameter',str(cp),index['source']['url'],source)
            for name,unit in [('dry_resistance_m2K_W','m2 K/W'),('evaporative_resistance_m2kPa_W','m2 kPa/W')]:
                for i,value in enumerate(model['last_step'][name]):
                    add('jos3',run['id']+'/'+name+'/'+str(i),value,unit,'source_model_boundary_state',str(cp),index['source']['url'],source)
    p=ROOT/'calibration/skin-fit.json'
    if p.exists():
        inputs.append(p);model=read(p)
        for i,name in enumerate(model['feature_names']):
            raw={'calibration_context':model['model_equation'],'source_sha256':sha256(p),
                 'parameter_se':model['uncertainty']['parameter_se'][i],'metrics':model['metrics'],
                 'external_validation':False,'causal':False}
            add('human-wound-phenotype-fit',name,model['fit']['parameters'][i],'V/m',
                'human_observation_phenotypic_fit',str(p),json.dumps(model['provenance']['references']),raw)
    anatomy=read(ROOT/'anatomy/opensim_index.json');atlas=read(ROOT/'anatomy/bodyparts3d_index.json')
    vascular=read(ROOT/'vascular/vmr_index.json');population=read(ROOT/'population/nhanes-2017-2018/index.json')
    semantics=read(ROOT/'semantics/index.json');phys=read(ROOT/'physiology/collection.json')
    raw_files=[p for p in Path('data/raw').rglob('*') if p.is_file() and '.git' not in p.parts and '.venv' not in p.parts]
    summary={'schema_version':1,'raw_file_count_excluding_git':len(raw_files),'raw_bytes_excluding_git':sum(p.stat().st_size for p in raw_files),
      'searchable_parameter_constraint_records':sum(counts.values()),'record_counts_by_source':dict(counts),'record_counts_by_basis':dict(bases),
      'anatomy':{'opensim':anatomy['counts'],'bodyparts3d':atlas['counts']},
      'vascular':{'human_cases':len(vascular['cases']),'raw_bytes':vascular['local_raw_bytes'],
                  'case_ids':[c['case_id'] for c in vascular['cases']]},
      'population':{k:v for k,v in population.items() if k!='tables'},
      'semantics':{k:semantics[k] for k in ('asctb_rows','vessels','geometry_measurements','vessel_mesh_annotation_candidates','upstream_malformed_citation_rows')},
      'physiology':phys,'whole_body_calibrated':False,'whole_body_registered_3d':False,
      'status':'source-backed native, population, circuit, temporal and observation-level materializations; independent integrated calibration incomplete',
      'limits':['source models are not one subject or one jointly calibrated system',
                'upstream fitted parameters are author-reported, not independently refitted',
                'simulation fields are not raw in-vivo measurements',
                'anatomical annotations do not define a coordinate transform',
                'published model parameters may overlap across variants; counts are not unique human structures']}
    native_summary=ROOT/'physiology/biogears_native_run/summary.json'
    if native_summary.exists():
        inputs.append(native_summary)
        summary['native_physiology_execution']=read(native_summary)
    skin_index=ROOT/'integumentary/index.json'
    if skin_index.exists():
        inputs.append(skin_index)
        summary['human_skin_experiment']=read(skin_index)
    db.execute('INSERT INTO metadata VALUES (?,?)',('summary',json.dumps(summary)))
    db.execute('INSERT INTO metadata VALUES (?,?)',('inputs',json.dumps([{'path':str(p),'sha256':sha256(p)} for p in inputs])))
    db.executescript('CREATE INDEX parameter_source ON parameters(source); CREATE INDEX parameter_name ON parameters(name);')
    db.commit();db.close();temporary.replace(path)
    (ROOT/'collection-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ('raw_file_count_excluding_git','raw_bytes_excluding_git','searchable_parameter_constraint_records','record_counts_by_source','record_counts_by_basis')}))


if __name__=='__main__':main()
