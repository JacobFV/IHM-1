"""Register retained Arm26 paths and wraps onto existing native arm segments.

Source-only XML engineering: no engine initialization, optimization or new body
mass. Full native moment-arm and force verification is a separate acceptance.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from scipy.spatial.transform import Rotation
from ihm.assembly.sensorimotor_catalog import native_muscle_catalog
ARM_SOURCE='data/raw/anatomy/opensim-models/source/Models/Arm26/arm26.osim'
TARGET_SOURCE='data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled.osim'
SCHEMA_SOURCE='data/raw/mechanics/opensim-core/OpenSim/Simulation/Model/Force.cpp'


def vector(element,tag):
    result=np.fromstring(element.findtext(tag,''),sep=' ')
    if result.shape!=(3,) or not np.isfinite(result).all():raise ValueError('Invalid source vector '+tag)
    return result


def setvector(element,tag,values):
    child=element.find(tag)
    if child is None:child=ET.SubElement(element,tag)
    child.text=' '.join(format(float(x),'.17g') for x in values)


def frame(model,joint,body):
    j=model.find(f'.//JointSet/objects/*[@name="{joint}"]')
    if j is None:raise ValueError('Missing source joint '+joint)
    matches=[f for f in j.findall('./frames/PhysicalOffsetFrame') if f.findtext('socket_parent')=='/bodyset/'+body]
    if len(matches)!=1:raise ValueError('Ambiguous source joint frame')
    return matches[0]


def align(a,b):
    a=a/np.linalg.norm(a);b=b/np.linalg.norm(b)
    cross=np.cross(a,b);cos=float(a@b)
    if cos<-1+1e-10:raise ValueError('Antiparallel source axes need explicit registration')
    skew=np.array([[0,-cross[2],cross[1]],[cross[2],0,-cross[0]],[-cross[1],cross[0],0]])
    return np.eye(3)+skew+skew@skew/(1+cos)


def build_upperbody_model(root,output):
    root=Path(root).resolve();output=Path(output).resolve()
    if not output.is_relative_to(root):raise ValueError('Output must be retained within root')
    output.mkdir(parents=True,exist_ok=True)
    if (output/'registration.json').exists():raise ValueError('Preserve existing upperbody materialization')
    donor=ET.parse(root/ARM_SOURCE).getroot();target=ET.parse(root/TARGET_SOURCE).getroot()
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    sources={p:sha(root/p) for p in (ARM_SOURCE,TARGET_SOURCE,SCHEMA_SOURCE)}
    if donor.get('Version')!='40000' or target.get('Version')!='40500':
        raise ValueError('Unreviewed donor/target OpenSim XML migration')
    shoulder=vector(frame(donor,'r_shoulder','base'),'translation')
    elbow=vector(frame(donor,'r_elbow','r_humerus'),'translation')
    donor_muscles=donor.findall('.//ForceSet/objects/Thelen2003Muscle')
    if len(donor_muscles)!=6:raise ValueError('Unexpected donor muscle coverage')
    bodies={b.get('name'):b for b in target.findall('.//BodySet/objects/Body')}
    forces=target.find('.//ForceSet/objects');catalog=native_muscle_catalog(root);registrations=[]
    for side in ('r','l'):
        target_shoulder=vector(frame(target,'acromial_'+side,'torso'),'translation')
        target_elbow=vector(frame(target,'elbow_'+side,'humerus_'+side),'translation')
        mirror=np.diag([1.,1.,1. if side=='r' else -1.])
        scale=float(np.linalg.norm(target_elbow)/np.linalg.norm(elbow))
        rotation=align(elbow,mirror@target_elbow)
        linear=scale*mirror@rotation
        radius_offset=vector(frame(target,'radioulnar_'+side,'ulna_'+side),'translation')
        # Parent/child offset orientations cancel at q=0 for this target source.
        for joint,parent,child in [('elbow_'+side,'humerus_'+side,'ulna_'+side),('radioulnar_'+side,'ulna_'+side,'radius_'+side)]:
            pf=frame(target,joint,parent);cf=frame(target,joint,child)
            if not np.allclose(vector(pf,'orientation'),vector(cf,'orientation'),atol=1e-12) or not np.allclose(vector(cf,'translation'),0):
                raise ValueError('Unresolved neutral source joint transform')
        def map_point(body,point,muscle_name=None):
            if body=='base':return 'torso',target_shoulder+linear@(point-shoulder)
            if body=='r_humerus':return 'humerus_'+side,linear@point
            if body=='r_ulna_radius_hand':
                # The donor fuses forearm/hand. Biceps radial insertion is an
                # anatomical assignment, not recovered from that fused label.
                radial=muscle_name in ('BIClong','BICshort')
                return ('radius_' if radial else 'ulna_')+side,linear@point-(radius_offset if radial else 0.)
            raise ValueError('Unresolved donor attachment '+body)
        wraps=[]
        for body in donor.findall('.//BodySet/objects/Body'):
            for original in body.findall('./WrapObjectSet/objects/*'):
                w=deepcopy(original);name=original.get('name');w.set('name','arm26_'+name+'_'+side)
                target_body,position=map_point(body.get('name'),vector(original,'translation'))
                setvector(w,'translation',position)
                original_rotation=Rotation.from_euler('XYZ',vector(original,'xyz_body_rotation')).as_matrix()
                # Reflected local geometry restores determinant +1; dimensions
                # are invariant, while the local z quadrant reverses sign.
                wrap_rotation=mirror@rotation@original_rotation@mirror
                if not np.allclose(wrap_rotation.T@wrap_rotation,np.eye(3),atol=1e-12) or np.linalg.det(wrap_rotation)<0:
                    raise ValueError('Improper wrap orientation')
                setvector(w,'xyz_body_rotation',Rotation.from_matrix(wrap_rotation).as_euler('XYZ'))
                if side=='l' and w.findtext('quadrant') in ('z','-z'):
                    w.find('quadrant').text='-z' if w.findtext('quadrant')=='z' else 'z'
                for field in ('radius','length'):
                    if w.find(field) is not None:w.find(field).text=format(float(w.findtext(field))*scale,'.17g')
                if w.find('dimensions') is not None:setvector(w,'dimensions',vector(w,'dimensions')*scale)
                parent=bodies[target_body].find('WrapObjectSet')
                if parent is None:parent=ET.SubElement(bodies[target_body],'WrapObjectSet')
                objects=parent.find('objects')
                if objects is None:objects=ET.SubElement(parent,'objects')
                objects.append(w)
                wraps.append({'name':w.get('name'),'source_name':name,'body':target_body,'type':w.tag,
                              'rotation_matrix':wrap_rotation.tolist(),'source_quadrant':original.findtext('quadrant'),'quadrant':w.findtext('quadrant')})
        for original in donor_muscles:
            m=deepcopy(original);name=original.get('name');key='arm26_'+name+'_'+side;m.set('name',key)
            geometry_paths=m.findall('./GeometryPath')
            if len(geometry_paths)!=1 or geometry_paths[0].get('name')!='geometrypath':
                raise ValueError('Unexpected donor path property schema')
            # Force::updateFromXMLNode performs this exact migration for a
            # pre-40500 document. Embedding the donor into the 40500 target
            # bypasses that document-version migration, so apply it here.
            # It is a concrete one-object property; do not add a <path> node.
            geometry_paths[0].set('name','path')
            path=[]
            for point in m.findall('.//PathPointSet/objects/*'):
                if point.tag!='PathPoint':raise ValueError('Moving/conditional donor point needs explicit transform')
                source_body=point.findtext('socket_parent_frame','').rsplit('/',1)[-1]
                source_point=vector(point,'location');body,position=map_point(source_body,source_point,name)
                point.set('name',key+'_'+point.get('name','point'))
                point.find('socket_parent_frame').text='/bodyset/'+body;setvector(point,'location',position)
                path.append({'source_body':source_body,'source_point_m':source_point.tolist(),'body':body,'point_m':position.tolist()})
            for wrap in m.findall('.//PathWrapSet/objects/*'):
                wrap.find('wrap_object').text='arm26_'+wrap.findtext('wrap_object')+'_'+side
            for field in ('optimal_fiber_length','tendon_slack_length'):
                m.find(field).text=format(float(m.findtext(field))*scale,'.17g')
            forces.append(m)
            hemi='lh' if side=='r' else 'rh'
            catalog.append({'id':key,'side':side,'body_group':'shoulder_elbow','attachment_bodies':sorted({a['body'] for a in path}),
                'sensory_region':f'brain-{hemi}-postcentral','motor_region':f'brain-{hemi}-precentral',
                'max_isometric_force_n':float(m.findtext('max_isometric_force')),
                'optimal_fiber_length_m':float(m.findtext('optimal_fiber_length')),'tendon_slack_length_m':float(m.findtext('tendon_slack_length')),
                'source_path':ARM_SOURCE,'source_sha256':sources[ARM_SOURCE],'source_muscle_name':name,
                'source_muscle_law':m.tag,'path_points':path,
                'assignment_basis':'contralateral regional cortical engineering prior; no identified recruitment',
                'geometry_basis':'donor similarity transfer using shoulder/elbow landmarks; left mirrored prior; biceps radius insertion assigned from named anatomy',
                'parameter_basis':'source Thelen coefficients and Fmax preserved; fiber/slack lengths scaled by humeral landmark ratio',
                'native_force_validation':False})
        registrations.append({'side':side,'length_scale':scale,'linear_map':linear.tolist(),
            'source_shoulder_m':shoulder.tolist(),'target_shoulder_m':target_shoulder.tolist(),
            'source_elbow_m':elbow.tolist(),'target_elbow_m':target_elbow.tolist(),
            'elbow_landmark_residual_m':float(np.linalg.norm(linear@elbow-target_elbow)),
            'radius_offset_from_ulna_m':radius_offset.tolist(),'wraps':wraps,'reference_pose':'all target arm coordinates zero',
            'reflection_prior':side=='l','independent_registration_landmarks':2,
            'unresolved':'3D moment-arm validity, independent scapula motion, regional morphology and intersubject strength scaling'})
    model_path=output/'subject_with_arms.osim';ET.ElementTree(target).write(model_path,encoding='utf-8',xml_declaration=True)
    catpath=output/'catalog.json';catpath.write_text(json.dumps(catalog,indent=2)+'\n')
    manifest={'schema':'ihm.upperbody-registration.v1','model_path':str(model_path.relative_to(root)),'model_sha256':sha(model_path),
        'catalog_path':str(catpath.relative_to(root)),'catalog_sha256':sha(catpath),'sources':sources,
        'builder_sha256':sha(__file__),'registrations':registrations,'muscle_count':92,'added_upperbody_muscles':12,
        'xml_migration':{'source_document_version':40000,'target_document_version':40500,
            'law_source':SCHEMA_SOURCE,'law_source_sha256':sources[SCHEMA_SOURCE],
            'operation':'Force::updateFromXMLNode: GeometryPath name geometrypath -> path; concrete element retained without wrapper',
            'applied_to':12,'numerical_parameters_changed_by_migration':False},
        'body_count':len(bodies),'body_mass_modified':False,'joint_tree_modified':False,
        'source_urls':{'Arm26':'https://github.com/opensim-org/opensim-models/blob/master/Models/Arm26/arm26.osim',
                       'law':'https://github.com/opensim-org/opensim-core/blob/main/OpenSim/Actuators/Thelen2003Muscle.cpp'},
        'license':'Arm26: CC BY 3.0; Holzbaur, Murray, Delp, and OpenSim adaptation team credited in source',
        'scope':'source-derived bilateral shoulder/elbow actuators; XML validity only until native mechanical tests; no new mass, no torque surrogate',
        'native_verified':False}
    (output/'registration.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    out=root/'data/derived/mechanics/whole_body_arm26_v2'
    report=build_upperbody_model(root,out)
    print(json.dumps({'path':str(out),'muscles':report['muscle_count'],'native_verified':False}))
