"""Append the coherent source Rajagopal default pose to the local app manifest."""
import sys,json,gzip,hashlib,os
from pathlib import Path
import numpy as np
import trimesh
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.spatial.opensim import load_model
from ihm.spatial.vtk import surface,read_arrays
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/derived/app';MODEL='opensim-rajagopal'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(path,value):
    b=json.dumps(value,separators=(',',':'),allow_nan=False).encode();path=Path(path)
    if str(path).endswith('.gz'):
        with path.open('wb') as f:
            with gzip.GzipFile(fileobj=f,mode='wb',mtime=0) as z:z.write(b)
    else:path.write_bytes(b)
def source_surface(path):
    # VTK point arrays may carry arbitrary names; one source mesh does so.
    arrays=read_arrays(path);point_keys=[k for k in arrays if k.startswith('Points/')]
    if len(point_keys)!=1:raise ValueError('Expected exactly one point array')
    points=arrays[point_keys[0]];conn=arrays['Polys/connectivity'];offsets=arrays['Polys/offsets'];faces=[];start=0
    for stop in offsets:
        poly=conn[start:int(stop)];start=int(stop)
        if len(poly)<3:raise ValueError('Invalid source polygon')
        faces.extend((poly[0],poly[i],poly[i+1]) for i in range(1,len(poly)-1))
    faces=np.asarray(faces,dtype=int)
    if faces.min()<0 or faces.max()>=len(points):raise ValueError('Invalid mesh index')
    return points,faces

def build():
    source=ROOT/'data/raw/anatomy/opensim-models/source/Models/Rajagopal/Rajagopal2016.osim';m=load_model(source)
    source_meta=json.loads((ROOT/'data/derived/anatomy/opensim__Rajagopal__Rajagopal2016.json').read_text())
    if sha(source)!=source_meta['sha256']:raise ValueError('Source model identity changed')
    native_path=ROOT/'data/derived/opensim/native_corrected/baseline/summary.json'
    native=json.loads(native_path.read_text()) if native_path.exists() else None
    native_muscles={}
    if native:
        verification=json.loads((native_path.parent.parent/'verification.json').read_text())
        if not verification['moment_arm_consistency_passed'] or native['execution']['model_sha256']!=sha(source):raise ValueError('Native mechanics verification or identity mismatch')
        native_muscles={x['name']:x for x in native['muscles']}
        for name,t in native['body_transforms_ground'].items():np.testing.assert_allclose(t,m['frames'][name],atol=1e-9,rtol=0)
    rotation=np.array([[0,0,-1],[0,1,0],[1,0,0.]])
    meshes=[]
    for entry in m['meshes']:
        points,faces=source_surface(entry['path']);t=m['frames'][entry['body']]
        ground=(points*entry['scale'])@t[:3,:3].T+t[:3,3]
        meshes.append((entry,ground@rotation.T,faces))
    all_points=np.vstack([v for _,v,_ in meshes]);center=(all_points.max(0)+all_points.min(0))/2
    (OUT/'geometry').mkdir(parents=True,exist_ok=True);structures=[]
    source_record={'label':'Rajagopal2016 · published OpenSim model','url':'https://github.com/opensim-org/opensim-models/tree/master/Models/Rajagopal','sha256':sha(source),'source_revision':source_meta['source_revision'],'frame':'opensim-ground-m','specimen':'Rajagopal2016 source model; not the BodyParts3D subject','units':'m','status':'source default pose; no independent calibration'}
    def structure(id,name,system,kind,p,extra):
        structures.append({'id':id,'name':name,'model_id':MODEL,'system':system,'kind':kind,'geometry_url':'/api/geometry/'+id,'geometry_sha256':sha(p),'color':'#dcc9a7' if system=='skeletal' else '#c56272','default_visible':True,'source':source_record,'calibration_status':'published source geometry; not subject calibrated',**extra})
    for entry,vertices,faces in meshes:
        id=MODEL+'-bone-'+entry['name'];p=OUT/'geometry'/(id+'.json.gz')
        mesh=trimesh.Trimesh(vertices=vertices-center,faces=faces,process=False)
        original_faces=len(faces)
        if len(faces)>2500:mesh=mesh.simplify_quadric_decimation(face_count=2500)
        dump(p,{'positions':np.asarray(mesh.vertices).ravel().tolist(),'indices':np.asarray(mesh.faces).ravel().tolist(),'normals':np.asarray(mesh.vertex_normals).ravel().tolist(),'units':'m','original_faces':original_faces,'display_decimation':len(mesh.faces)<original_faces})
        structure(id,entry['path'].stem+' · '+entry['body'],'skeletal','mesh',p,{'body_frame':entry['body'],'geometry_source':str(entry['path'].relative_to(ROOT)),'geometry_source_sha256':sha(entry['path'])})
    for muscle in m['muscles']:
        points=np.asarray([p['ground_m'] for p in muscle['points']])@rotation.T-center
        pairs=np.asarray([[a,b] for a,b in zip(points,points[1:])]).reshape(-1,3)
        id=MODEL+'-muscle-'+muscle['name'];p=OUT/'geometry'/(id+'.json.gz')
        payload={'positions':pairs.ravel().tolist(),'units':'m','attachment_points':muscle['points'],'wrap_objects':muscle['wraps'],'wrap_solved':False,'path_interpretation':'straight segments between source attachment/via points; wrapping not solved','activation':None,'force_N':None,'moment_arm_m':None}
        extra={'path_interpretation':'attachment/via-point chords; no solved wrapping or force','wrap_solved':False,'wrap_count':len(muscle['wraps']),'source_parameters':{k:v for k,v in muscle.items() if k.endswith('_N') or k.endswith('_m')}}
        label='attachment path'
        if native:
            dump(p.with_name(id+'-attachment-fallback.json.gz'),payload)
            result=native_muscles[muscle['name']];points=np.asarray(result['path_ground_m'])@rotation.T-center
            pairs=np.asarray([[a,b] for a,b in zip(points,points[1:])]).reshape(-1,3)
            mechanics={k:result[k] for k in ['length_m','tendon_force_N','active_fiber_force_N','passive_fiber_force_N','fiber_length_m','moment_arms_m','moment_arm_check']}
            mechanics.update(activation=native['activation'],engine_variant=native['execution']['engine_variant'],evidence_kind='source_simulation',external_force_balance_solved=False)
            payload.update(positions=pairs.ravel().tolist(),wrap_solved=True,path_interpretation='native OpenSim wrapped path; cache invalidation correction applied',native_mechanics=mechanics,force_N=result['tendon_force_N'],activation=native['activation'])
            extra.update(wrap_solved=True,path_interpretation=payload['path_interpretation'],native_mechanics=mechanics,native_source={'engine_revision':native['execution']['engine_source_revision'],'simbody_revision':native['execution']['simbody_source_revision'],'library_sha256':native['execution']['simulation_library_sha256'],'export_sha256':sha(native_path)},calibration_status='source muscle equilibrium at explicit activation; not subject calibrated')
            label='native wrapped path'
        dump(p,payload)
        structure(id,muscle['name']+' · '+label,'muscular','lines',p,extra)
    model={'id':MODEL,'name':'Rajagopal2016 · coherent musculoskeletal model','description':'Source default-coordinate skeleton and muscle attachment paths, in meters. Wrap paths and dynamics are unsolved. Separate specimen from reference anatomy.',
        'frame':'opensim-rajagopal-display-m','source_frame':'opensim-ground-m','source_units':'m','display_units':'m','calibration_status':'published model; no subject-specific calibration',
        'bounds':{'min':(all_points.min(0)-center).tolist(),'max':(all_points.max(0)-center).tolist()},'display_transform':{'rotation':rotation.tolist(),'scale':1,'translation':(-center).tolist()},'attribution':source_meta['publications'],'source':source_record}
    if native:
        model['description']='Native OpenSim default pose with wrapped muscle paths, equilibrium forces and moment arms; explicit cache correction independently tested. Separate source specimen.'
        model['mechanics_source']=native['execution']
    manifest_path=OUT/'manifest.json';manifest=json.loads(manifest_path.read_text());manifest['models']=[x for x in manifest['models'] if x['id']!=MODEL]+[model]
    manifest['structures']=[x for x in manifest['structures'] if x['model_id']!=MODEL]+structures
    temporary=manifest_path.with_suffix('.opensim.tmp');dump(temporary,manifest);os.replace(temporary,manifest_path)
    report={'schema_version':1,'source':source_record,'default_coordinates':m['coordinates'],'body_transforms_ground':{k:v.tolist() for k,v in m['frames'].items()},'joint_count':m['joint_count'],'mesh_count':len(meshes),'muscle_count':len(m['muscles']),'attachment_point_count':sum(len(x['points']) for x in m['muscles']),'wrap_reference_count':sum(len(x['wraps']) for x in m['muscles']),'unresolved_default_frames':m['unresolved'],'wraps_solved':bool(native),'muscle_equilibrium_solved':bool(native),'dynamics_solved':False,'moment_arms_computed':bool(native),'native_export':str(native_path.relative_to(ROOT)) if native else None,'limitations':['SimmSpline evaluation restricted to exact source knots; non-knot defaults fail.','UniversalJoint supported only at zero default coordinates.','Source constraints checked for default consistency; no nonlinear assembly.',('Native wrapped paths and static muscle equilibrium; skeletal external force balance and time integration not solved.' if native else 'Muscle lines are attachment/via-point chords; wrapping and dynamics are not solved.'),'Display pose is source standing/default pose, not supine physiology.'],'reference_source_sha256':{p.name:sha(p) for p in (ROOT/'data/derived/opensim/reference').glob('*') if p.is_file()}}
    dump(ROOT/'data/derived/opensim/default_pose.json',report)
    print(json.dumps({k:report[k] for k in ['joint_count','mesh_count','muscle_count','attachment_point_count','wrap_reference_count','unresolved_default_frames']},indent=2))
if __name__=='__main__':build()
