"""Restricted exact XML default-pose evaluator; never approximates spline/wrap solutions.

Frame convention X_GB maps homogeneous coordinates in B to Ground. Supported
SimmSpline values must be exact source knots; this suffices for Rajagopal2016's
complete default pose. This is not an OpenSim dynamics or muscle wrapping engine.
"""
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from scipy.spatial.transform import Rotation

def vector(text,n=3):
    a=np.asarray([float(x) for x in (text or '').split()])
    if a.shape!=(n,) or not np.isfinite(a).all():raise ValueError('Invalid source vector')
    return a

def offset_transform(translation,orientation):
    t=np.asarray(translation,float);r=np.asarray(orientation,float)
    if t.shape!=(3,) or r.shape!=(3,) or not np.isfinite([t,r]).all():raise ValueError('Invalid offset transform')
    out=np.eye(4);out[:3,:3]=Rotation.from_euler('XYZ',r).as_matrix();out[:3,3]=t
    return out

def evaluate_function(f,args):
    args=np.asarray(args,float)
    if not np.isfinite(args).all():raise ValueError('Nonfinite coordinate')
    if f.tag=='Constant':return float(f.findtext('value'))
    if f.tag=='LinearFunction':
        c=vector(f.findtext('coefficients'),len(args)+1);return float(c[:-1]@args+c[-1])
    if f.tag=='SimmSpline':
        x=np.asarray([float(v) for v in f.findtext('x').split()]);y=np.asarray([float(v) for v in f.findtext('y').split()])
        if len(args)!=1 or x.shape!=y.shape or len(x)<2 or not np.isfinite([x,y]).all() or not np.all(np.diff(x)>0):raise ValueError('Invalid source spline')
        # Exact interpolation-node reproduction needs no assumption about endpoint derivatives.
        at=np.flatnonzero(x==args[0])
        if len(at)!=1:raise ValueError('SimmSpline off-knot evaluation requires upstream engine; unresolved')
        return float(y[at[0]])
    raise ValueError(f'Unsupported function {f.tag}; unresolved')

def _function(node):
    functions=[c for c in node if c.tag not in ('coordinates','axis')]
    if len(functions)!=1:raise ValueError('Expected one TransformAxis function')
    return functions[0]

def solve_frames(bodies,joints):
    frames={'ground':np.eye(4)};pending=list(joints)
    if len(set(bodies))!=len(bodies):raise ValueError('Duplicate body')
    while pending:
        progress=False
        for j in pending[:]:
            parent,child,xpf,xfm,xbm=j
            if parent not in frames:continue
            if child in frames or child not in bodies:raise ValueError('Duplicate/unknown child body')
            frames[child]=frames[parent]@xpf@xfm@np.linalg.inv(xbm)
            pending.remove(j);progress=True
        if not progress:raise ValueError('Unresolved cyclic or disconnected joint frames')
    if set(frames)!={'ground',*bodies}:raise ValueError('Unresolved body frames')
    return frames

def load_model(path):
    path=Path(path);model=ET.parse(path).getroot().find('Model')
    if model.findtext('length_units')!='meters':raise ValueError('Only explicit meter source units supported')
    bodies=list(model.find('BodySet/objects'));joints=list(model.find('JointSet/objects'))
    if any(b.tag!='Body' for b in bodies):raise ValueError('Unsupported BodySet component')
    coordinates={c.attrib['name']:float(c.findtext('default_value')) for j in joints for c in j.findall('coordinates/Coordinate')}
    if not np.isfinite(list(coordinates.values())).all():raise ValueError('Nonfinite coordinate default')
    for c in model.findall('ConstraintSet/objects/*'):
        if c.findtext('isEnforced','true')=='false':continue
        if c.tag!='CoordinateCouplerConstraint':raise ValueError('Unsupported constraint')
        names=c.findtext('independent_coordinate_names').split();dependent=c.findtext('dependent_coordinate_name')
        value=evaluate_function(list(c.find('coupled_coordinates_function'))[0],[coordinates[n] for n in names])*float(c.findtext('scale_factor','1'))
        if not np.isclose(coordinates[dependent],value,atol=1e-12,rtol=0):raise ValueError('Default coordinates violate source constraint; requires assembly')
    connections=[]
    def bodyname(socket):
        if socket=='/ground':return 'ground'
        if not socket.startswith('/bodyset/') or socket.count('/')!=2:raise ValueError('Unsupported frame socket '+socket)
        return socket.rsplit('/',1)[1]
    for joint in joints:
        offsets={f.attrib['name']:f for f in joint.findall('frames/PhysicalOffsetFrame')}
        pf=offsets[joint.findtext('socket_parent_frame')];cf=offsets[joint.findtext('socket_child_frame')]
        transform=lambda f:offset_transform(vector(f.findtext('translation')),vector(f.findtext('orientation')))
        motion=np.eye(4);q=[coordinates[c.attrib['name']] for c in joint.findall('coordinates/Coordinate')]
        if joint.tag=='CustomJoint':
            axes=joint.findall('SpatialTransform/TransformAxis')
            if [a.attrib['name'] for a in axes]!=['rotation1','rotation2','rotation3','translation1','translation2','translation3']:raise ValueError('Unsupported spatial axis order')
            for i,a in enumerate(axes):
                axis=vector(a.findtext('axis'));length=np.linalg.norm(axis)
                if length==0:raise ValueError('Zero transform axis')
                axis=axis/length;names=(a.findtext('coordinates') or '').split();v=evaluate_function(_function(a),[coordinates[n] for n in names])
                if i<3:motion[:3,:3]=motion[:3,:3]@Rotation.from_rotvec(axis*v).as_matrix()
                else:motion[:3,3]+=axis*v
        elif joint.tag=='PinJoint' and len(q)==1:
            motion[:3,:3]=Rotation.from_rotvec([0,0,q[0]]).as_matrix()
        elif joint.tag=='UniversalJoint' and len(q)==2 and np.all(np.asarray(q)==0):pass
        else:raise ValueError('Unsupported joint/default coordinate '+joint.tag)
        connections.append((bodyname(pf.findtext('socket_parent')),bodyname(cf.findtext('socket_parent')),transform(pf),motion,transform(cf)))
    frames=solve_frames([b.attrib['name'] for b in bodies],connections)
    meshes=[]
    for body in bodies:
        if body.findall('components/*'):raise ValueError('Body component frames unresolved')
        for mesh in body.findall('attached_geometry/*'):
            if mesh.tag!='Mesh' or mesh.findtext('socket_frame')!='..':raise ValueError('Unsupported attached geometry frame')
            file=path.parent/'Geometry'/mesh.findtext('mesh_file')
            if not file.is_file():file=path.parents[2]/'Geometry'/mesh.findtext('mesh_file')
            if not file.is_file():raise ValueError('Missing referenced geometry '+str(file))
            meshes.append({'name':mesh.attrib['name'],'body':body.attrib['name'],'path':file,'scale':vector(mesh.findtext('scale_factors'))})
    muscles=[]
    for force in model.findall('ForceSet/objects/*'):
        geometry=force.find('GeometryPath')
        if geometry is None:continue
        points=[]
        for point in geometry.findall('PathPointSet/objects/*'):
            if point.tag!='PathPoint':raise ValueError('Moving/conditional path point requires supported coordinate evaluator')
            body=bodyname(point.findtext('socket_parent_frame'));local=vector(point.findtext('location'))
            world=frames[body]@np.r_[local,1]
            points.append({'name':point.attrib['name'],'body':body,'local_m':local.tolist(),'ground_m':world[:3].tolist()})
        wraps=[{'name':w.attrib['name'],'wrap_object':w.findtext('wrap_object'),'method':w.findtext('method'),'range':w.findtext('range')} for w in geometry.findall('PathWrapSet/objects/*')]
        muscles.append({'name':force.attrib['name'],'source_type':force.tag,'points':points,'wraps':wraps,'wrap_solved':False,'max_isometric_force_N':float(force.findtext('max_isometric_force')),'optimal_fiber_length_m':float(force.findtext('optimal_fiber_length')),'tendon_slack_length_m':float(force.findtext('tendon_slack_length'))})
    return {'frames':frames,'coordinates':coordinates,'meshes':meshes,'muscles':muscles,'model_name':model.attrib['name'],'joint_count':len(joints),'unresolved':[]}
