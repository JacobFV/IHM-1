"""Build real, decimated 3D display assets; physical source identities remain separate."""
import csv
import gzip
import json
from pathlib import Path
from collections import Counter
import numpy as np
import trimesh
from ihm.forge.acquisition import sha256
from ihm.spatial.vtk import surface,read_arrays

OUT=Path('data/derived/app');GEOM=OUT/'geometry'
COLORS={'cardiac':'#cd6474','skeletal':'#d9c6a6','muscular':'#b95864','arterial':'#ef725f','venous':'#598ed4','nervous':'#e8c568',
'lymphatic':'#78b991','respiratory':'#99bad1','digestive':'#c49a67','urinary':'#a47fad','reproductive':'#d58cba',
'endocrine':'#dab25c','integumentary':'#d4ac98','connective':'#8dafa6','other':'#9badb4'}

def dump(path,data):
    data=json.dumps(data,separators=(',',':'),allow_nan=False).encode()
    if str(path).endswith('.gz'):
        with Path(path).open('wb') as f:
            with gzip.GzipFile(fileobj=f,mode='wb',mtime=0) as z:z.write(data)
    else:Path(path).write_bytes(data)

def classify(names):
    text=' | '.join(names).lower()
    for system,words in [('cardiac',['heart','myocardium','cardiac','atrioventricular','wall of atrium','wall of ventricle']),('skeletal',['bone organ','tooth','cartilage']),('muscular',['muscle organ','muscle tissue']),
      ('arterial',['artery','arterial']),('venous',['vein','venous']),('lymphatic',['lymph','spleen','thymus']),
      ('nervous',['nerve','brain','spinal cord','ganglion']),('respiratory',['lung','bronch','trachea','larynx','nasal']),
      ('urinary',['kidney','ureter','urinary','renal pelvis']),('reproductive',['testis','penis','prostate','seminal','epididym','ductus deferens']),
      ('endocrine',['thyroid','adrenal','pituitary','pineal']),('digestive',['liver','pancreas','stomach','intestin','colon','rectum','esophag','gallbladder','duoden','jejun','ileum']),
      ('integumentary',['skin','nail','hair','mammary']),('connective',['ligament','tendon','fascia','aponeurosis'])]:
        if any(word in text for word in words):return system
    return 'other'

def mesh_payload(vertices,faces,max_faces=700):
    m=trimesh.Trimesh(vertices=vertices,faces=faces,process=False)
    if len(faces)>max_faces:m=m.simplify_quadric_decimation(face_count=max_faces)
    return {'positions':np.round(m.vertices,7).ravel().tolist(),'indices':m.faces.ravel().tolist(),
            'normals':np.round(m.vertex_normals,6).ravel().tolist(),'display_decimation':True,'original_faces':len(faces)}

def build():
    GEOM.mkdir(parents=True,exist_ok=True);structures=[];models=[]
    old=json.loads((OUT/'manifest.json').read_text()) if (OUT/'manifest.json').exists() else {}
    old_structures={s['id']:s for s in old.get('structures',[])}
    old_models={m['id']:m for m in old.get('models',[])}
    atlas=json.loads(Path('data/derived/anatomy/bodyparts3d_index.json').read_text())
    counts=Counter(c['concept_id'] for m in atlas['meshes'] for c in m['concepts'])
    rotation=np.array([[1,0,0],[0,0,1],[0,-1,0.]])
    all_bounds=np.array([m['bounds_in_source_coordinates'] for m in atlas['meshes']]).reshape(-1,3)@rotation.T*.001
    center=(all_bounds.min(0)+all_bounds.max(0))/2
    model={'id':'bodyparts3d','name':'BodyParts3D · reference human','description':'Original coherent adult male atlas. Millimeter source coordinates verified from official diagram.',
           'frame':'bodyparts3d-display-m','source_frame':'bodyparts3d-mm-left-posterior-superior','source_units':'mm','display_units':'m',
           'calibration_status':'reference anatomy; not subject-specific','bounds':{'min':(all_bounds.min(0)-center).tolist(),'max':(all_bounds.max(0)-center).tolist()},
           'display_transform':{'rotation':rotation.tolist(),'scale':.001,'translation':(-center).tolist()},
           'coordinate_reference':'https://dbarchive.biosciencedbc.jp/data/bodyparts3d/20130619/coordinate_system.png',
           'attribution':'BodyParts3D, © The Database Center for Life Science, CC BY 4.0'}
    models.append(model)
    for i,entry in enumerate(atlas['meshes']):
        name=min(entry['concepts'],key=lambda c:(counts[c['concept_id']],-len(c['name'])))['name'] if entry['concepts'] else entry['element_id']
        system=classify([c['name'] for c in entry['concepts']]);id='bp3d-'+entry['element_id']
        p=GEOM/(id+'.json.gz')
        cached=old_structures.get(id,{})
        cache_valid=(p.exists() and cached.get('source',{}).get('sha256')==entry['sha256']
                     and cached.get('geometry_sha256')==sha256(p)
                     and old_models.get(model['id'],{}).get('display_transform')==model['display_transform'])
        if not cache_valid:
            m=trimesh.load(entry['source_path'],force='mesh',process=False)
            dump(p,mesh_payload(np.asarray(m.vertices)@rotation.T*.001-center,np.asarray(m.faces)))
        structures.append({'id':id,'name':name,'system':system,'model_id':model['id'],'kind':'mesh','geometry_url':'/api/geometry/'+id,
            'color':COLORS[system],'calibration_status':'reference anatomical surface','concepts':entry['concepts'],
            'source':{'label':'BodyParts3D 4.0','url':'https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html',
                'sha256':entry['sha256'],'frame':model['source_frame'],'specimen':'adult male reference atlas','units':'mm','status':'acquired'},
            'geometry_sha256':sha256(p),'default_visible':system in ('cardiac','skeletal','respiratory','digestive','urinary')})
        if i%400==0:print('atlas',i,flush=True)
    vascular=json.loads(Path('data/derived/vascular/vmr_index.json').read_text())
    for case,modelid,title in [('0001_H_AO_SVD','vascular-aorta','Aorta · pediatric Fontan case'),('0050_H_CERE_H','vascular-cerebral','Cerebral vessels · VMR'),('0077_H_PULM_H','vascular-pulmonary','Pulmonary vessels · VMR')]:
        raw=Path('data/raw/vascular/vmr/extracted')/case;file=next(raw.rglob('walls_combined.vtp'))
        v,f=surface(file);center=(v.max(0)+v.min(0))/2;scale=2/np.ptp(v,axis=0).max()
        vm={'id':modelid,'name':title,'description':'Separate case-specific vascular surface; original CFD where available. Display normalized; physical units require case confirmation.',
            'frame':case+'-display','source_frame':case,'source_units':'unconfirmed source units','display_units':'normalized',
            'calibration_status':'archived source case; no independent validation',
            'bounds':{'min':((v.min(0)-center)*scale).tolist(),'max':((v.max(0)-center)*scale).tolist()},
            'display_transform':{'scale':float(scale),'translation':(-center*scale).tolist()},'source_metadata':next(c['source_metadata'] for c in vascular['cases'] if c['case_id']==case)}
        models.append(vm);id=modelid+'-wall';p=GEOM/(id+'.json.gz');dump(p,mesh_payload((v-center)*scale,f,12000))
        source={'label':'Vascular Model Repository '+case,'url':'https://www.vascularmodel.com/','sha256':sha256(file),'frame':case,
                'specimen':case,'units':'unconfirmed source units','status':'archived; saved job status Simulation failed'}
        structures.append({'id':id,'name':title+' wall','system':'arterial','model_id':modelid,'kind':'mesh','geometry_url':'/api/geometry/'+id,
                           'color':'#be6a76','source':source,'calibration_status':'case anatomy; archived job failure retained','geometry_sha256':sha256(p),'default_visible':True})
        if case=='0050_H_CERE_H':
            with np.load('data/derived/vascular/0050_H_CERE_H_first_cfd_frame.npz',allow_pickle=False) as d:
                points=d['Points/Points'];velocity=d['PointData/velocity_01010'];pressure=d['PointData/pressure_01010']
            valid=np.flatnonzero(np.linalg.norm(velocity,axis=1)>np.percentile(np.linalg.norm(velocity,axis=1),35))
            sample=valid[np.linspace(0,len(valid)-1,700,dtype=int)]
            id=modelid+'-flow';flow={'positions':np.round((points[sample]-center)*scale,7).ravel().tolist(),
                'velocities':velocity[sample].ravel().tolist(),'pressure':pressure[sample].tolist(),'vector_scale':.04,
                'times':[1.01],'frames':[],'source_node_indices':sample.tolist(),'units':'source velocity/pressure; not confirmed SI','time_basis':'solver.inp dt=.001 seconds × archived step number'}
            # Decode original archived states rather than synthesize vector animation.
            framefile=next(p for p in Path('data/raw/vascular').rglob('*.vtu') if p.stat().st_size>1_000_000_000)
            for step in range(1010,2011,5):
                values=read_arrays(framefile,names={f'PointData/velocity_{step:05}',f'PointData/pressure_{step:05}'})
                flow['frames'].append({'velocities':values[f'PointData/velocity_{step:05}'][sample].ravel().tolist(),
                                       'pressure':values[f'PointData/pressure_{step:05}'][sample].tolist()})
                if step%200==10:print('CFD original step',step,flush=True)
            flow['times']=[step*.001 for step in range(1010,2011,5)]
            p=GEOM/(id+'.json.gz');dump(p,flow)
            structures.append({'id':id,'name':'Archived cerebral flow · 201 states','system':'arterial','model_id':modelid,'kind':'vectors',
                'geometry_url':'/api/geometry/'+id,'color':'#73dfd4','source':{**source,'sha256':sha256(framefile)},
                'calibration_status':'source CFD output; not in-vivo measurement','geometry_sha256':sha256(p),'default_visible':True})
    manifest={'schema_version':1,'models':models,'structures':structures,'systems':[{'id':k,'name':k.title(),'color':v} for k,v in COLORS.items()],
        'bounds':models[0]['bounds'],'whole_body_calibrated':False,'cross_family_registration':False,
        'limitations':['Different model families have different specimens; selecting a common display does not register them.',
                      'Display decimation never replaces original physical geometry.','Vascular archived jobs contain failure flags; flow is not certified converged.']}
    dump(OUT/'manifest.json',manifest);print('manifest',len(structures),'structures',flush=True)

if __name__=='__main__':build()
