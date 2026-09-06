"""Explicit source-only cervical reference registration; never activates a model.

Only affine scalar functions and fixed translations are evaluated. Unsupported
source joints/functions fail closed instead of approximating native kinematics.
"""
import numpy as np
from scipy.spatial.transform import Rotation


def numbers(text):
    values=np.asarray([float(x) for x in (text or '').split()])
    if not np.isfinite(values).all():raise ValueError('Nonfinite source values')
    return values


def proper_fit(source,target):
    x=np.asarray(source,dtype=float);y=np.asarray(target,dtype=float)
    if x.shape!=y.shape or x.ndim!=2 or x.shape[1]!=3 or len(x)<3 or not np.isfinite([x,y]).all():
        raise ValueError('Finite paired 3D correspondences required')
    a=x-x.mean(0);b=y-y.mean(0)
    singular=np.linalg.svd(a,compute_uv=False)
    if singular[1]<max(singular[0]*1e-6,1e-12) or np.linalg.svd(b,compute_uv=False)[1]<1e-12:
        raise ValueError('Collinear registration cannot identify a rotation')
    u,s,vh=np.linalg.svd(a.T@b)
    sign=np.linalg.det(vh.T@u.T)
    r=vh.T@np.diag([1,1,sign])@u.T
    t=np.eye(4);t[:3,:3]=r;t[:3,3]=y.mean(0)-r@x.mean(0)
    residual=np.linalg.norm(x@r.T+t[:3,3]-y,axis=1)
    return t,{'rms_m':float(np.sqrt(np.mean(residual**2))),'max_m':float(residual.max()),
              'residuals_m':residual.tolist(),'source_singular_values_m':singular.tolist(),
              'cross_covariance_singular_values_m2':s.tolist(),'determinant':float(np.linalg.det(r)),
              'reflection_correction_applied':bool(sign<0)}


def scalar_function(function,x):
    if function is None or not np.isfinite(x):raise ValueError('Missing or nonfinite function')
    if function.tag=='Constant':return float(numbers(function.findtext('value'))[0])
    if function.tag=='LinearFunction':
        c=numbers(function.findtext('coefficients'))
        if len(c)==1:return float(c[0])  # zero-argument linear function: constant term
        if len(c)==2:return float(c[0]*x+c[1])
    if function.tag=='SimmSpline':
        knots=numbers(function.findtext('x'));values=numbers(function.findtext('y'))
        if len(knots)==len(values)==2 and knots[1]>knots[0]:
            return float(values[0]+(x-knots[0])*(values[1]-values[0])/(knots[1]-knots[0]))
    raise ValueError('Unsupported function; only Constant, affine LinearFunction and two-knot SimmSpline: '+function.tag)


def reference_joint(joint,coordinates):
    if joint.tag not in ('CustomJoint','WeldJoint') or joint.findtext('reverse','false')!='false':
        raise ValueError('Unsupported joint/reverse: '+joint.tag)
    def frame(location,orientation):
        p=numbers(joint.findtext(location));a=numbers(joint.findtext(orientation))
        if p.shape!=(3,) or a.shape!=(3,):raise ValueError('Missing joint frame')
        t=np.eye(4);t[:3,:3]=Rotation.from_euler('XYZ',a).as_matrix();t[:3,3]=p;return t
    mobility=np.eye(4)
    axes=joint.findall('./SpatialTransform/TransformAxis')
    if joint.tag=='CustomJoint':
        if [a.get('name') for a in axes]!=['rotation1','rotation2','rotation3','translation1','translation2','translation3']:
            raise ValueError('Unsupported transform axis ordering')
        for index,axis in enumerate(axes):
            names=(axis.findtext('coordinates') or '').split()
            if len(names)>1 or any(n not in coordinates for n in names):raise ValueError('Unresolved axis coordinate')
            function=axis.find('function')
            if function is None or len(function)!=1:raise ValueError('Ambiguous axis function')
            value=scalar_function(function[0],coordinates[names[0]] if names else 0.)
            vector=numbers(axis.findtext('axis'))
            if vector.shape!=(3,) or not np.isclose(np.linalg.norm(vector),1):raise ValueError('Non-unit transform axis')
            if index<3:mobility[:3,:3]=mobility[:3,:3]@Rotation.from_rotvec(vector*value).as_matrix()
            elif value!=0:raise ValueError('Nonzero mobility translation unsupported')
    elif axes:raise ValueError('Weld with unexpected transform')
    return frame('location_in_parent','orientation_in_parent')@mobility@np.linalg.inv(frame('location','orientation'))


def reference_frames(model,names,coordinates):
    bodies={b.get('name'):b for b in model.findall('./BodySet/objects/Body')}
    frames={'spine':np.eye(4)};pending=list(names)
    while pending:
        progress=False
        for name in list(pending):
            body=bodies[name];joint=body.find('Joint')
            if joint is None or len(joint)!=1:raise ValueError('Expected one source joint')
            parent=joint[0].findtext('parent_body')
            if parent not in frames:continue
            frames[name]=frames[parent]@reference_joint(joint[0],coordinates)
            pending.remove(name);progress=True
        if not progress:raise ValueError('Unresolved source joint graph: '+str(pending))
    return frames


def validate_recipe(recipe,target_model_bytes):
    """Reject a recipe for any other native target, and independently check ledger."""
    import hashlib
    from ihm.assembly.cervical_inertia import transform_body,combine_bodies
    if hashlib.sha256(target_model_bytes).hexdigest()!=recipe['target_model_identity']['sha256']:
        raise ValueError('Target model identity mismatch; no silent mass normalization')
    if recipe['native_activation_allowed'] is not False:raise ValueError('Source recipe cannot authorize native activation')
    names={b['name'] for b in recipe['bodies'].values()}
    if len(names)!=9 or len(recipe['coordinates'])!=24 or len(recipe['constraints'])!=18 or len(recipe['muscles'])!=78:
        raise ValueError('Unexpected cervical component inventory')
    total=combine_bodies([recipe['replace_target_torso_inertia']]+[
        transform_body(b['mass_properties_in_new_body_frame'],b['body_to_target_torso_reference']) for b in recipe['bodies'].values()])
    import xml.etree.ElementTree as ET
    torso=ET.fromstring(target_model_bytes).find('./Model/BodySet/objects/Body[@name="torso"]')
    xx,yy,zz,xy,xz,yz=numbers(torso.findtext('inertia'));inertia=np.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])
    if not np.isclose(total['mass_kg'],float(torso.findtext('mass')),rtol=0,atol=1e-10) or not np.allclose(total['center_m'],numbers(torso.findtext('mass_center')),rtol=0,atol=1e-10) or not np.allclose(total['inertia_kg_m2'],inertia,rtol=0,atol=1e-10):
        raise ValueError('Recipe violates original torso mass, first moment or inertia')
    count=0
    for muscle in recipe['muscles']:
        if muscle['name'] in names:raise ValueError('Duplicate component identity')
        names.add(muscle['name'])
        for p in muscle['path_points']:
            if p['target_body'] not in {'torso'}|{b['name'] for b in recipe['bodies'].values()}:
                raise ValueError('Unresolved muscle body')
            if numbers(' '.join(map(str,p['target_location_m']))).shape!=(3,):raise ValueError('Invalid attachment')
            count+=1
    if count!=198:raise ValueError('Unexpected attachment inventory')
    return {'identity_bound':True,'conserved_torso_mass_kg':total['mass_kg'],'attachments':count}
