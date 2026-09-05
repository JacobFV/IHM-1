"""preserve atlas cell/biomarker relations and vascular anatomical measurements.

Anatomical relation edges are not flow directions or electrical coupling weights.
Geometry.csv has malformed trailing citation fields in upstream: numeric prefixes
are validated and retained, while unresolved citation tails remain explicitly raw.
"""
import csv
import hashlib
import io
import json
from pathlib import Path
import sqlite3
import subprocess
import zipfile
from ihm.forge.acquisition import sha256

ROOT=Path('data/raw/semantics'); OUT=Path('data/derived/semantics')


def canon(value):
    if value.startswith('FMA') and value[3:].isdigit(): return 'FMA:'+value[3:]
    for prefix in ('http://purl.obolibrary.org/obo/','https://purl.obolibrary.org/obo/'):
        if value.startswith(prefix): return value[len(prefix):].replace('_',':',1)
    return value


def main():
    OUT.mkdir(parents=True,exist_ok=True); output=OUT/'human-reference-atlas.sqlite'
    db=sqlite3.connect(output)
    db.executescript('DROP TABLE IF EXISTS asctb; DROP TABLE IF EXISTS vessel; DROP TABLE IF EXISTS mesh_annotation; CREATE TABLE asctb(anatomy_label TEXT,cell_label TEXT,biomarker_label TEXT,anatomy_id TEXT,cell_id TEXT,biomarker_id TEXT,biomarker_type TEXT); CREATE TABLE vessel(name TEXT,parent TEXT,kind TEXT,anatomy_id TEXT,raw_json TEXT); CREATE TABLE mesh_annotation(mesh_id TEXT,anatomy_id TEXT,label TEXT,source_path TEXT);')
    archive=ROOT/'asctb-records-v0.11.1.csv.zip'; count=0
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None: raise ValueError('ASCT+B archive CRC mismatch')
        with z.open('asctb-records.csv') as stream:
            reader=csv.DictReader(io.TextIOWrapper(stream,encoding='utf-8-sig')); batch=[]
            for row in reader:
                batch.append((row['as_label'],row['ct_label'],row['bm_label'],canon(row['as']),canon(row['ct']),row['bm'],row['bmType']))
                if len(batch)==10000:
                    db.executemany('INSERT INTO asctb VALUES (?,?,?,?,?,?,?)',batch); count+=len(batch);batch=[]
            db.executemany('INSERT INTO asctb VALUES (?,?,?,?,?,?,?)',batch);count+=len(batch)
    with (ROOT/'hra-vccf/Vessel.csv').open(newline='',encoding='utf-8-sig') as stream:
        vessels=list(csv.DictReader(stream))
    db.executemany('INSERT INTO vessel VALUES (?,?,?,?,?)',[(r['Vessel'],r['BranchesFrom'],r['VesselType'],canon(r['ASID']),json.dumps(r)) for r in vessels])
    atlas=Path('data/derived/anatomy/bodyparts3d_index.json'); annotations=0
    if atlas.exists():
        for mesh in json.loads(atlas.read_text())['meshes']:
            rows=[(mesh['element_id'],canon(c['concept_id']),c['name'],mesh['source_path']) for c in mesh['concepts']]
            db.executemany('INSERT INTO mesh_annotation VALUES (?,?,?,?)',rows);annotations+=len(rows)
    db.executescript('CREATE INDEX asctb_anatomy ON asctb(anatomy_id); CREATE INDEX vessel_anatomy ON vessel(anatomy_id); CREATE INDEX mesh_anatomy ON mesh_annotation(anatomy_id);')
    candidate_count=db.execute('SELECT COUNT(*) FROM vessel v JOIN mesh_annotation m ON v.anatomy_id=m.anatomy_id WHERE v.anatomy_id!=""').fetchone()[0]
    candidates=[dict(zip(('vessel','anatomy_id','mesh','source_path'),row)) for row in db.execute('SELECT v.name,v.anatomy_id,m.mesh_id,m.source_path FROM vessel v JOIN mesh_annotation m ON v.anatomy_id=m.anatomy_id WHERE v.anatomy_id!=""')]
    (OUT/'vascular-mesh-annotation-candidates.json').write_text(json.dumps({'status':'semantic_association_only_not_coordinate_registration','includes_ancestor_annotations':True,'candidates':candidates},indent=2)+'\n')
    db.commit();db.close()
    geometry=ROOT/'hra-vccf/Geometry.csv'; measurements=[]; malformed=0
    with geometry.open(newline='',encoding='utf-8-sig') as stream:
        reader=csv.reader(stream);header=next(reader)
        for line,row in enumerate(reader,start=2):
            if len(row)<13: raise ValueError(f'geometry numeric/identity prefix truncated at {line}')
            prefix=dict(zip(header[:13],row[:13])); bad=len(row)!=len(header); malformed+=int(bad)
            def number(key):
                value=prefix[key].strip()
                if value in ('','-1'): return None
                return float(value)
            measurements.append({'source':'HRA-VCCF','source_file':str(geometry),'source_sha256':sha256(geometry),
                'line':line,'entity':prefix['Vessel'],'anatomy_id':canon(prefix['VesselID']),
                'property':prefix['Property'],'sex':prefix['Sex'],'value':number('Value'),
                'standard_deviation':number('StandardDeviation'),'range_low':number('RangeLow'),
                'range_high':number('RangeHigh'),'units':prefix['Units'],'population':prefix['Population'],
                'sample_size':number('SampleSize'),'method':prefix['Method'],'reference':prefix['ReferenceURL'],
                'basis':'atlas_curated_literature_geometry_summary','primary_study_independently_verified':False,
                'citation_columns_malformed':bad,'citation_tail_raw':row[13:]})
    with (OUT/'vascular-geometry-measurements.jsonl').open('w') as f:
        for row in measurements:f.write(json.dumps(row)+'\n')
    revision=subprocess.check_output(['git','-C',str(ROOT/'hra-vccf'),'rev-parse','HEAD'],text=True).strip()
    inventory=[]
    for path in sorted((ROOT/'hra-vccf').glob('*')):
        if path.is_file(): inventory.append({'path':str(path),'sha256':sha256(path),'bytes':path.stat().st_size})
    report={'asctb_rows':count,'asctb_archive_sha256':sha256(archive),'vessels':len(vessels),
        'mesh_annotations_including_ancestors':annotations,'vessel_mesh_annotation_candidates':candidate_count,
        'geometry_measurements':len(measurements),'upstream_malformed_citation_rows':malformed,
        'vccf_revision':revision,'files':inventory,'vccf_data_license':'CC-BY-4.0','vccf_code_license':'MIT',
        'source_references':['https://github.com/hubmapconsortium/hra-vccf','https://github.com/x-atlas-consortia/hra-pop/blob/main/output-data/v0.11.1/reports/hra/asctb-records.csv.zip'],
        'limitations':['semantic crosswalk is not geometric registration','branch relations are not validated blood flow directions',
          'upstream malformed citation suffixes preserved; validated numeric prefixes retained','atlas populations differ from NHANES and mesh donors']}
    (OUT/'index.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('files','limitations','source_references')}))


if __name__=='__main__': main()
