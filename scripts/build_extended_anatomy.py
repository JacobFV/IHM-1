#!/usr/bin/env python3
"""Export actual source anatomy in its own frame; optionally append app manifest.

Run with the project venv. Blender executes the --extract stage with autoexec off.
Default writes an isolated manifest fragment. --append-manifest is explicit integration.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
if '--extract' in sys.argv:
    sys.path.insert(0, str(ROOT / '.cache/blender-python'))
import numpy as np
from collect_extended_anatomy import REVISION as PINNED_REVISION, EXPECTED_SHA256 as PINNED_FILES
PINNED_BLEND_SHA256 = "9f08a17ea0115fed80b2a73ecdf0a1bc2ab2f6956f37c593ce23d513ea35afcd"
EXTRACTOR_SEMANTIC_VERSION = 1
EXTRACTION_DESCRIPTION = "Original datablocks and source-authored evaluated surfaces retained separately. Original viewport modifiers, levels and curve bevels evaluated unchanged. Triangulated without decimation. Object matrix_world applied; reflected face winding reversed only for negative determinant."
RAW = ROOT / 'data/raw/anatomy/extended'
OUT = ROOT / 'data/derived/anatomy/extended'
MODEL_ID = 'z-anatomy'
COLORS = {'skeletal':'#d9c6a6','muscular':'#b95864','arterial':'#ef725f','venous':'#598ed4','cardiac':'#cd6474',
    'lymphatic':'#78b991','respiratory':'#99bad1','digestive':'#c49a67','urinary':'#a47fad',
    'reproductive':'#d58cba','endocrine':'#dab25c','integumentary':'#d4ac98','connective':'#8dafa6','other':'#9badb4'}

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def require_digest(path, expected, label):
    path = Path(path)
    if not path.is_file() or sha256(path) != expected:
        raise ValueError(f'{label} checksum mismatch or missing file: {path}')


def validate_provenance():
    """Check actual pinned raw bytes before extraction, cached builds or integration."""
    provenance=json.loads((RAW/'provenance.json').read_text())
    if provenance.get('revision') != PINNED_REVISION:
        raise ValueError('Raw source revision differs from pinned revision')
    records=provenance.get('files',[])
    if len(records)!=len(PINNED_FILES) or {Path(r['path']).name for r in records}!=set(PINNED_FILES):
        raise ValueError('Incomplete or duplicate raw source provenance records')
    for record in records:
        path=ROOT/record['path'];expected=PINNED_FILES[path.name]
        if record.get('sha256')!=expected:
            raise ValueError(f'Pinned raw source digest mismatch: {path}')
        require_digest(path,expected,'raw source')
        if path.stat().st_size != record.get('bytes'):
            raise ValueError(f'Raw source size mismatch: {path}')
    if provenance.get('blend_sha256')!=PINNED_BLEND_SHA256:
        raise ValueError('Pinned blend digest mismatch')
    require_digest(ROOT/provenance['blend_path'],PINNED_BLEND_SHA256,'blend')
    return provenance


def validate_cached_sources(provenance):
    """Fail closed on stale, incompatible or altered intermediate geometry.

    Existing schema-2 exports predate the explicit semantic version. Their exact
    recorded extraction description identifies the same v1 algorithm; runtime
    validation/display-only edits do not invalidate these historical exports.
    """
    index=json.loads((OUT/'source_index.json').read_text())
    version=index.get('extractor_semantic_version')
    if index.get('schema_version')!=2 or (version is None and index.get('extraction')!=EXTRACTION_DESCRIPTION) or (version is not None and version!=EXTRACTOR_SEMANTIC_VERSION):
        raise ValueError('Unsupported cached extractor semantics; run --force-extract')
    if index.get('source_blend_sha256')!=provenance['blend_sha256']:
        raise ValueError('Source index differs from pinned blend')
    entries=index.get('meshes',[])
    if not entries or len({e['id'] for e in entries})!=len(entries):
        raise ValueError('Empty or duplicate cached source identities')
    for entry in entries:
        for stage in ('source','base'):
            require_digest(ROOT/entry[f'{stage}_geometry_path'],entry[f'{stage}_geometry_sha256'],f'{stage} geometry')
    return index


def validate_fragment(fragment, provenance, index):
    entries={e['id']:e for e in index['meshes']}
    structures=fragment.get('structures',[])
    if len(structures)!=len(entries) or {s['id'] for s in structures}!=set(entries):
        raise ValueError('Fragment identities differ from validated source index')
    for structure in structures:
        entry=entries[structure['id']];source=structure['source']
        if source.get('sha256')!=provenance['blend_sha256'] or source.get('geometry_sha256')!=entry['source_geometry_sha256'] or source.get('base_geometry_sha256')!=entry['base_geometry_sha256']:
            raise ValueError(f'Fragment source provenance differs: {structure["id"]}')
        require_digest(ROOT/structure['geometry_path'],structure['geometry_sha256'],'display geometry')


def write(path, data):
    path.parent.mkdir(parents=True,exist_ok=True)
    raw=json.dumps(data,separators=(',',':'),allow_nan=False).encode()
    if str(path).endswith('.gz'):
        with path.open('wb') as f:
            with gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as z:z.write(raw)
    else:path.write_bytes(raw)

def category(name, collections):
    c=set(collections);text=name.lower()
    if '6: Lymphoid organs' in c:return 'lymphatic'
    if 'parathyroid' in text:return 'endocrine'
    if any(w in text for w in ['pharynx','nasal cavity']):return 'respiratory'
    if any(w in text for w in ['omentum','meso-appendix','mesocolon','soft palate']):return 'digestive'
    if '9: Regions of human body' in c:return 'integumentary'
    if '1: Skeletal system' in c:return 'skeletal'
    if '3: Joints' in c:return 'connective'
    if '4: Muscular system' in c:
        return 'connective' if any(w in text for w in ['bursa','fascia','tendon','aponeurosis','sheath','retinacul']) else 'muscular'
    if '5: Cardiovascular system' in c:
        if any('vein' in v.lower() for v in c) or any(w in text for w in ['vein','venous','sinus']):return 'venous'
        if any('arter' in v.lower() or v=='Aorta' for v in c):return 'arterial'
        return 'cardiac'
    for label,system in [('Endocrine glands','endocrine'),('Respiratory system','respiratory'),('Urinary system','urinary'),("Genital systems'",'reproductive'),('Digestive system','digestive')]:
        if label in c:return system
    return 'other'

def rejection(obj):
    if obj.type not in {'MESH','CURVE'}:return 'non-geometric label/camera/light'
    if re.search(r'\.(j|g|t|i)(\.\d+)?$',obj.name):return 'source annotation/leader/group-label geometry'
    cs={c.name for c in obj.users_collection}
    if '7: Nervous system & Sense organs' in cs:return 'sense-organ/nervous collection excluded: third-party license scope'
    if 'kidney' in obj.name.lower():return 'kidney excluded: third-party NC license scope'
    if not cs.intersection({'1: Skeletal system','3: Joints','4: Muscular system','5: Cardiovascular system','6: Lymphoid organs','8: Visceral systems','9: Regions of human body'}):return 'outside selected anatomical collections (including insertion annotations and unlinked objects)'
    return None

def modifier_metadata(obj):
    result=[]
    for modifier in obj.modifiers:
        record={}
        for prop in modifier.bl_rna.properties:
            if prop.identifier=='rna_type':continue
            try:value=getattr(modifier,prop.identifier)
            except Exception:continue
            if isinstance(value,(str,bool,int,float)) or value is None:record[prop.identifier]=value
            elif hasattr(value,'name'):record[prop.identifier]={'datablock_name':value.name}
        result.append(record)
    return result

def extract():
    import bpy
    provenance=validate_provenance()
    if Path(bpy.data.filepath).resolve() != (ROOT/provenance['blend_path']).resolve():
        raise ValueError('Loaded Blender document is not the validated pinned source')
    destination=OUT/'source_geometry';destination.mkdir(parents=True,exist_ok=True)
    base_destination=OUT/'base_geometry';base_destination.mkdir(parents=True,exist_ok=True)
    entries=[];excluded=[]
    dg=bpy.context.evaluated_depsgraph_get()
    def arrays(mesh,matrix):
        mesh.calc_loop_triangles()
        local=np.empty((len(mesh.vertices),3),dtype=np.float64)
        mesh.vertices.foreach_get('co',local.ravel())
        world=local@np.array(matrix,dtype=np.float64)[:3,:3].T+np.array(matrix,dtype=np.float64)[:3,3]
        faces=np.empty((len(mesh.loop_triangles),3),dtype=np.int32)
        mesh.loop_triangles.foreach_get('vertices',faces.ravel())
        # Preserve original winding in local arrays; reflect world winding only for negative determinants.
        world_faces=faces[:,::-1] if matrix.to_3x3().determinant()<0 else faces
        return local,world,faces,world_faces
    for obj in sorted(bpy.data.objects,key=lambda o:o.name):
        reason=rejection(obj)
        if reason:
            excluded.append({'name':obj.name,'type':obj.type,'reason':reason});continue
        identity='za-'+hashlib.sha256(obj.name.encode()).hexdigest()[:16]
        base_mesh=obj.data if obj.type=='MESH' else obj.to_mesh()
        base_local,base_v,base_f,base_world_f=arrays(base_mesh,obj.matrix_world)
        base_path=base_destination/(identity+'.npz')
        np.savez_compressed(base_path,local_vertices=base_local,vertices=base_v,local_faces=base_f,faces=base_world_f)
        if obj.type=='CURVE':obj.to_mesh_clear()
        evaluated=obj.evaluated_get(dg)
        mesh=evaluated.to_mesh()
        local,vertices,local_faces,faces=arrays(mesh,evaluated.matrix_world)
        if not len(vertices) or not len(faces):
            excluded.append({'name':obj.name,'type':obj.type,'reason':'empty evaluated surface','base_geometry_path':str(base_path.relative_to(ROOT))})
            evaluated.to_mesh_clear();continue
        path=destination/(identity+'.npz')
        np.savez_compressed(path,local_vertices=local,vertices=vertices,local_faces=local_faces,faces=faces)
        collections=sorted(c.name for c in obj.users_collection)
        entries.append({'id':identity,'name':obj.name,'source_object':obj.name,'source_type':obj.type,
            'system':category(obj.name,collections),'source_collections':collections,
            'source_geometry_path':str(path.relative_to(ROOT)),'source_geometry_sha256':sha256(path),
            'base_geometry_path':str(base_path.relative_to(ROOT)),'base_geometry_sha256':sha256(base_path),
            'base_vertices':len(base_v),'base_triangles':len(base_f),
            'base_bounds':[base_v.min(0).tolist(),base_v.max(0).tolist()] if len(base_v) else None,
            'source_vertices':len(vertices),'source_triangles':len(faces),'bounds':[vertices.min(0).tolist(),vertices.max(0).tolist()],
            'matrix_world':[list(row) for row in evaluated.matrix_world],
            'original_matrix_world':[list(row) for row in obj.matrix_world],
            'evaluated_modifiers':modifier_metadata(obj)})
        evaluated.to_mesh_clear()
        if len(entries)%300==0:print('evaluated source geometry',len(entries),flush=True)
    inventory=[{'name':o.name,'type':o.type,'source_collections':sorted(c.name for c in o.users_collection),
        'matrix_world':[list(row) for row in o.matrix_world],'mesh_vertices':len(o.data.vertices) if o.type=='MESH' else None,
        'mesh_polygons':len(o.data.polygons) if o.type=='MESH' else None,
        'modifiers':modifier_metadata(o), 'export_exclusion':rejection(o)} for o in bpy.data.objects]
    write(OUT/'raw_object_inventory.json',{'objects':inventory})
    write(OUT/'source_index.json',{'schema_version':2,'extractor_semantic_version':EXTRACTOR_SEMANTIC_VERSION,'blender_version':bpy.app.version_string,
        'source_blend_version':list(bpy.data.version),'dependency_graph_mode':dg.mode,
        'evaluation_log':'data/derived/anatomy/extended/blender_evaluation.log',
        'source_blend_sha256':sha256(bpy.data.filepath),
        'unit_settings':{'system':bpy.context.scene.unit_settings.system,'scale_length':bpy.context.scene.unit_settings.scale_length},
        'extraction':EXTRACTION_DESCRIPTION,
        'meshes':entries,'excluded':excluded})
    print('source complete',len(entries),flush=True)

def build(force_extract=False):
    import trimesh
    OUT.mkdir(parents=True,exist_ok=True)
    provenance=validate_provenance()
    if force_extract or not (OUT/'source_index.json').exists():
        with (OUT/'blender_evaluation.log').open('w') as evaluation_log:
            subprocess.run(['blender','--background','--disable-autoexec',str(ROOT/provenance['blend_path']),
            '--python-exit-code','1','--python',str(Path(__file__).resolve()),'--','--extract'],check=True,
            env={**os.environ,'PYTHONHOME':os.environ.get('EXTENDED_BLENDER_PYTHONHOME','/usr')},
            stdout=evaluation_log,stderr=subprocess.STDOUT)
    index=validate_cached_sources(provenance)
    entries=index['meshes'];rotation=np.array([[1,0,0],[0,0,1],[0,-1,0.]])
    bounds=np.array([e['bounds'] for e in entries]).reshape(-1,3)@rotation.T
    center=(bounds.min(0)+bounds.max(0))/2
    # Do not assume Blender scene METRIC is independently measured SI. Normalize one atlas only.
    scale=2/np.ptp(bounds,axis=0).max()
    model={'id':MODEL_ID,'name':'Z-Anatomy · extended reference atlas',
        'description':'Authored reference anatomy derived partly from BodyParts3D; includes regional lymph nodes, muscle surfaces and body surface regions. Separate source frame; no subject calibration.',
        'frame':'z-anatomy-display-normalized','source_frame':'z-anatomy-blender-world','source_units':'Blender units (scene METRIC, scale_length 1; physical dimensions unverified)',
        'display_units':'normalized','bounds':{'min':((bounds.min(0)-center)*scale).tolist(),'max':((bounds.max(0)-center)*scale).tolist()},
        'display_transform':{'rotation':rotation.tolist(),'scale':float(scale),'translation':(-center*scale).tolist()},
        'calibration_status':'reference geometry; display normalized; no independent human calibration',
        'attribution':provenance['attribution'],'license':'CC-BY-SA-4.0','license_url':provenance['license_url'],
        'source_revision':provenance['revision'],'evaluation':{'blender_version':index['blender_version'],'source_blend_version':index['source_blend_version'],'mode':index['dependency_graph_mode'],'log':index['evaluation_log']},'independent_subject_count':0,'cross_family_registration':False}
    structures=[]
    for i,e in enumerate(entries):
        with np.load(ROOT/e['source_geometry_path'],allow_pickle=False) as data:
            v=data['vertices'];f=data['faces']
        original_faces=len(f);m=trimesh.Trimesh(vertices=(v@rotation.T-center)*scale,faces=f,process=False)
        path=OUT/'geometry'/(e['id']+'.json.gz')
        write(path,{'positions':np.round(m.vertices,8).ravel().tolist(),'indices':m.faces.ravel().tolist(),
            'normals':np.round(m.vertex_normals,6).ravel().tolist(),'display_decimation':False,'original_faces':original_faces,'source_triangles':original_faces,'base_triangles':e['base_triangles'],'output_triangles':len(m.faces),'full_resolution':True})
        is_region=e['system']=='integumentary'
        structures.append({'id':e['id'],'name':e['name'],'model_id':MODEL_ID,'system':e['system'],'kind':'mesh',
            'geometry_url':'/api/geometry/'+e['id'],'geometry_path':str(path.relative_to(ROOT)),'geometry_sha256':sha256(path),
            'display_geometry':{'resolution':'full-source-evaluated','source_faces':e['source_triangles'],'base_faces':e['base_triangles'],'target_faces':e['source_triangles']},
            'color':COLORS[e['system']],'default_visible':e['system']=='lymphatic',
            'calibration_status':'authored body surface region; no skin thickness or layers' if is_region else 'authored reference anatomical surface',
            'concepts':[],'source_collections':e['source_collections'],
            'evaluation_warnings':['Source oesophagus/profile dependency cycle; see blender_evaluation.log'] if e['name']=='Oesophagus' else [],
            'source':{'label':'Z-Anatomy','url':provenance['repository']+'/tree/'+provenance['revision'],
                'sha256':provenance['blend_sha256'],'revision':provenance['revision'],'object_name':e['source_object'],
                'geometry_sha256':e['source_geometry_sha256'],'base_geometry_sha256':e['base_geometry_sha256'],'geometry_stage':'source-authored evaluated viewport surface','frame':model['source_frame'],'units':model['source_units'],
                'specimen':'authored reference atlas; partly derived from BodyParts3D, not an additional measured human',
                'status':'acquired','source_category':provenance['source_category'],'license':'CC-BY-SA-4.0',
                'license_url':provenance['license_url'],'attribution':provenance['attribution']}})
        if i%500==0:print('display geometry',i,flush=True)
    counts=dict(sorted(Counter(s['system'] for s in structures).items()))
    limitations=['Z-Anatomy partly derives from BodyParts3D: overlapping anatomy is not an independent human or independent validation.',
        'No cross-family registration is assumed. Display normalization does not establish physical dimensions.',
        'Regional body surface meshes do not provide skin thickness, histological layers, hair follicles or sweat-gland geometry.',
        'This archive has lymph-node group surfaces, spleen, thymus and tonsils, but no named thoracic duct, cisterna chyli or lymphatic-vessel geometry.',
        'Node-group objects may contain several nodes; object counts are not counts of individual anatomical nodes.',
        'Blender reports a source oesophagus/profile dependency cycle; affected evaluation is retained with a warning and is not independently resolved.',
        'Source meshes are triangulated without decimation; original Blender archive and full source mesh exports are retained. Source-authored viewport modifiers are evaluated at their original settings; pre-modifier geometry is retained separately.',
        provenance['third_party_license_caveat']]
    write(OUT/'manifest_fragment.json',{'schema_version':1,'models':[model],'structures':structures,
        'systems':[{'id':k,'name':k.title(),'color':v} for k,v in COLORS.items()],
        'cross_family_registration':False,'limitations':limitations})
    write(OUT/'coverage.json',{'structure_count':len(structures),'counts_by_system':counts,
        'base_triangles':sum(e['base_triangles'] for e in entries),'source_triangles':sum(e['source_triangles'] for e in entries),'output_triangles':sum(e['source_triangles'] for e in entries),'display_decimation':False,'source_vertices':sum(e['source_vertices'] for e in entries),
        'independent_subject_count':0,'lymph_node_group_objects':sum(e['system']=='lymphatic' and 'node' in e['name'].lower() for e in entries),
        'lymphatic_vessel_objects':0,'excluded_objects':len(index['excluded']),'limitations':limitations})
    print(json.dumps({'structure_count':len(structures),'counts_by_system':counts},indent=2))

def append_manifest(path):
    """Explicit integration: copy geometry and idempotently replace only this family."""
    manifest=json.loads(path.read_text());fragment=json.loads((OUT/'manifest_fragment.json').read_text())
    provenance=validate_provenance()
    index=validate_cached_sources(provenance)
    validate_fragment(fragment,provenance,index)
    geometry=path.parent/'geometry';geometry.mkdir(parents=True,exist_ok=True)
    for s in fragment['structures']:shutil.copy2(ROOT/s['geometry_path'],geometry/(s['id']+'.json.gz'))
    manifest['models']=[m for m in manifest['models'] if m['id']!=MODEL_ID]+fragment['models']
    manifest['structures']=[s for s in manifest['structures'] if s['model_id']!=MODEL_ID]+fragment['structures']
    manifest['limitations']=list(dict.fromkeys(manifest.get('limitations',[])+fragment['limitations']))
    manifest['cross_family_registration']=False
    write(path,manifest)
    print('Integrated Z-Anatomy:',len(fragment['structures']))

if __name__=='__main__':
    if '--extract' in sys.argv:extract()
    else:
        parser=argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--force-extract',action='store_true');parser.add_argument('--append-manifest',type=Path)
        parser.add_argument('--integrate-only',action='store_true');args=parser.parse_args()
        if not args.integrate_only:build(args.force_extract)
        if args.append_manifest:append_manifest(args.append_manifest)
