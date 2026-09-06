"""Bounded XML-only donor audit; does not import OpenSim or run a model."""
from collections import Counter
import hashlib,json,math,sys
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
NECK={f'cerv{i}' for i in range(1,8)}|{'skull','jaw'}


def identity(path):
    raw=path.read_bytes()
    return {'path':str(path.relative_to(ROOT)),'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}


def vector(text):return [float(x) for x in text.split()]


def inertia(body):
    entries=vector(body.findtext('inertia')) if body.find('inertia') is not None else [
        float(body.findtext('inertia_'+key)) for key in ('xx','yy','zz','xy','xz','yz')]
    xx,yy,zz,xy,xz,yz=entries
    matrix=np.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])
    eig=np.linalg.eigvalsh(matrix)
    # I=trace(S)1-S, where S is the nonnegative second mass moment.
    second=np.trace(matrix)/2*np.eye(3)-matrix
    valid=bool(np.linalg.eigvalsh(second).min()>=-1e-12)
    return {'tensor_kg_m2_xx_yy_zz_xy_xz_yz':entries,'principal_moments_kg_m2':eig.tolist(),
            'positive_definite':bool(eig.min()>0),'physical_mass_distribution_possible':valid,
            'second_mass_moment_eigenvalues_kg_m2':np.linalg.eigvalsh(second).tolist()}


def model(path):
    root=ET.parse(path).getroot();m=root.find('Model')
    bodies=m.findall('./BodySet/objects/Body');forces=m.findall('./ForceSet/objects/*')
    couplers=m.findall('./ConstraintSet/objects/CoordinateCouplerConstraint')
    cervical=[];coordinates=[]
    for body in bodies:
        if body.get('name') not in NECK:continue
        joint=list(body.find('Joint'))[0]
        coords=[{'name':c.get('name'),'default':float(c.findtext('default_value')),
                 'range_rad':vector(c.findtext('range')),'locked':c.findtext('locked')}
                for c in joint.findall('./CoordinateSet/objects/Coordinate')]
        coordinates.extend(c['name'] for c in coords)
        cervical.append({'name':body.get('name'),'mass_kg':float(body.findtext('mass')),
            'mass_center_body_m':vector(body.findtext('mass_center')),'inertia':inertia(body),
            'joint':{'name':joint.get('name'),'type':joint.tag,'parent':joint.findtext('parent_body'),
                     'frames':{k:vector(joint.findtext(k)) for k in ('location_in_parent','orientation_in_parent','location','orientation')},
                     'coordinates':coords,
                     'axes':[{'name':a.get('name'),'coordinate':a.findtext('coordinates'),
                              'axis':vector(a.findtext('axis')),
                              'function_xml':ET.tostring(a.find('function'),encoding='unicode').strip()}
                             for a in joint.findall('./SpatialTransform/TransformAxis')]}})
    neck_couplers=[{'name':c.get('name'),'independent':c.findtext('independent_coordinate_names'),
                    'dependent':c.findtext('dependent_coordinate_name'),
                    'function_xml':ET.tostring(c.find('coupled_coordinates_function'),encoding='unicode').strip()}
                   for c in couplers if c.findtext('dependent_coordinate_name') in coordinates]
    dependent={c['dependent'] for c in neck_couplers}
    muscles=[]
    for muscle in forces:
        if 'Muscle' not in muscle.tag:continue
        muscles.append({'name':muscle.get('name'),'type':muscle.tag,
            'max_isometric_force_n':float(muscle.findtext('max_isometric_force')),
            'optimal_fiber_length_m':float(muscle.findtext('optimal_fiber_length')),
            'tendon_slack_length_m':float(muscle.findtext('tendon_slack_length')),
            'path_bodies':[p.text.strip() for p in muscle.findall('./GeometryPath/PathPointSet/objects/*/body')]})
    return {'identity':identity(path),'xml_version':root.get('Version'),'name':m.get('name'),
            'credits':m.findtext('credits'),'publications':m.findtext('publications'),
            'units':{'length':m.findtext('length_units'),'force':m.findtext('force_units'),'gravity_m_s2':vector(m.findtext('gravity'))},
            'body_count_including_ground':len(bodies),'forces_by_type':dict(Counter(f.tag for f in forces)),
            'neck_head_bodies':cervical,'neck_head_mass_kg':sum(b['mass_kg'] for b in cervical),
            'neck_coordinate_count':len(coordinates),'neck_couplers':neck_couplers,
            'neck_independent_coordinates':[c for c in coordinates if c not in dependent],
            'muscles':muscles,'muscle_path_body_dependencies':sorted({b for x in muscles for b in x['path_bodies']}),
            'bushings':[{c.tag:c.text.strip() for c in b} | {'name':b.get('name')} for b in forces if b.tag=='BushingForce'],
            'visual_geometry_not_acquired':sorted({e.text for e in m.findall('.//geometry_file')}),
            'physical_inertia_failures':[b['name'] for b in cervical if not b['inertia']['physical_mass_distribution_possible']]}


def audit():
    manifest=json.loads((HERE/'sources.json').read_text())
    for entry in manifest['retained_files']:
        assert identity(ROOT/entry['path'])=={k:entry[k] for k in ('path','bytes','sha256')},entry['path']
    metadata=json.loads((HERE/'figshare_masi_4249808.json').read_text())
    payload=(HERE/'MASI_HMaleMuscle_HMaleMassDistr.osim').read_bytes()
    assert hashlib.md5(payload).hexdigest()==metadata['files'][0]['computed_md5']
    payload=(HERE/'Mortensen2018_reference.osim').read_bytes()
    assert hashlib.sha1(b'blob '+str(len(payload)).encode()+b'\0'+payload).hexdigest()=='41fb3623066aff0f7ce9fc4f5c5f13ae07589fe1'
    models=[model(HERE/name) for name in ('MASI_HMaleMuscle_HMaleMassDistr.osim','Mortensen2018_reference.osim')]
    for m in models:
        assert m['neck_coordinate_count']==24 and len(m['neck_couplers'])==18
        assert set(m['neck_independent_coordinates'])=={'pitch1','roll1','yaw1','pitch2','roll2','yaw2'}
        assert set(m['physical_inertia_failures'])=={f'cerv{i}' for i in range(1,8)}
        assert len(m['neck_head_bodies'])==9
    assert len(models[0]['muscles'])==78 and len(models[1]['muscles'])==72
    target=ROOT/'data/derived/mechanics/whole_body_arm26_v2/subject_with_arms.osim'
    m=ET.parse(target).getroot().find('Model');torso=m.find('./BodySet/objects/Body[@name="torso"]')
    held=ROOT/'data/raw/mechanics/opensim-core/OpenSim/Simulation/SimbodyEngine/Body.cpp'
    return {'schema':'ihm.cervical-source-audit.v1','scope':'XML and provenance audit only; no native execution or accepted donor inertia',
        'models':models,'target':{'identity':identity(target),'body_count':len(m.findall('./BodySet/objects/Body')),
        'torso_mass_kg':float(torso.findtext('mass')),'torso_mass_center_m':vector(torso.findtext('mass_center')),'torso_inertia':inertia(torso)},
        'opensim_invalid_inertia_handler':{'identity':identity(held),'lines':'146-163',
        'action':'Catches invalid inertia, logs warning/error, substitutes spherical inertia',
        'cervical_fallback_diagonal_kg_m2':math.sqrt(.02**2+.08**2+.02**2)/math.sqrt(3)}}


if __name__=='__main__':
    result=audit();destination=HERE/'audit.json'
    if '--check' in sys.argv:assert json.loads(destination.read_text())==result
    else:destination.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'verified':True,'models':len(result['models']),'native_runs':0,
                      'failed_inertia_bodies_per_model':[m['physical_inertia_failures'] for m in result['models']]}))
