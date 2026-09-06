"""Assemble an opt-in cervical JSON recipe, preserving source XML and identities.

No native initialization, default model mutation, coordinate atlas, or tissue
attachment accuracy claim. The output is not an OpenSim-loadable model.
"""
import argparse,gzip,hashlib,json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from ihm.assembly.cervical_registration import numbers,proper_fit,scalar_function,reference_frames
from ihm.assembly.cervical_inertia import convex_geometry_prior,transform_body,combine_bodies
from ihm.spatial.vtk import surface

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'data/research/cervical_inertia/v2/manifest.json'
DONOR=ROOT/'data/research/cervical/MASI_HMaleMuscle_HMaleMassDistr.osim'
GEOMETRY=ROOT/'data/raw/anatomy/opensim-models/source/Geometry'
NAMES=[f'cerv{i}' for i in range(7,0,-1)]+['skull','jaw']
ANCHORS=['torso','lclavicle','lscapula','rclavicle','rscapula']


def sha(raw):return hashlib.sha256(raw).hexdigest()
def encoded(element):return ET.tostring(element,encoding='unicode')
def point(t,p):return (np.asarray(t)@np.r_[p,1])[:3]


def build(output):
    output=Path(output).resolve()
    if output.exists() or not output.is_relative_to(ROOT):raise ValueError('Fresh workspace output required')
    captured={}
    def receipt(path):
        raw=path.read_bytes();digest=sha(raw);captured[digest]=raw
        return {'path':str(path.relative_to(ROOT)),'sha256':digest,'bytes':len(raw),'retained_copy':'inputs/'+digest+'.gz'}
    inputs=[receipt(p) for p in (PRIOR,DONOR,Path(__file__).resolve(),ROOT/'ihm/assembly/cervical_registration.py',ROOT/'ihm/assembly/cervical_inertia.py',ROOT/'ihm/spatial/vtk.py')]
    prior=json.loads(PRIOR.read_bytes())
    # Bind the entire prior dependency chain, including its original native model.
    for entry in prior['inputs']:
        raw=gzip.decompress((PRIOR.parent/entry['retained_copy']).read_bytes())
        if sha(raw)!=entry['sha256']:raise ValueError('Prior source archive mismatch')
    target_receipt=next(e for e in prior['inputs'] if e['path'].endswith('/native/assembled_model.osim'))
    model=ET.parse(DONOR).getroot().find('Model');bodies={b.get('name'):b for b in model.findall('./BodySet/objects/Body')}
    coordinates={c.get('name'):float(c.findtext('default_value')) for c in model.findall('.//CoordinateSet/objects/Coordinate')}
    neck_coordinates={c.get('name'):c for name in NAMES for c in bodies[name].findall('./Joint/*/CoordinateSet/objects/Coordinate')}
    coordinates.update({n:0. for n in neck_coordinates})
    constraints=[]
    for con in model.findall('./ConstraintSet/objects/*'):
        dependent=con.findtext('dependent_coordinate_name')
        if dependent not in neck_coordinates:continue
        independent=(con.findtext('independent_coordinate_names') or '').split()
        if con.tag!='CoordinateCouplerConstraint' or con.findtext('isDisabled')!='false' or len(independent)!=1 or independent[0] not in neck_coordinates:
            raise ValueError('Unsupported neck constraint')
        value=scalar_function(con.find('coupled_coordinates_function')[0],0.)
        if abs(value)>1e-14:raise ValueError('Zero neck pose violates donor coupling')
        constraints.append({'name':'cervical_'+con.get('name'),'source_xml':encoded(con),'dependent':'cervical_'+dependent,'independent':'cervical_'+independent[0]})
    if len(neck_coordinates)!=24 or len(constraints)!=18:raise ValueError('Expected 24 coordinates and 18 neck couplers')
    frames=reference_frames(model,NAMES+ANCHORS,coordinates)
    correspondences=[];geometry_records={};unsupported=[]
    for name in NAMES:
        visual=bodies[name].find('VisibleObject');geometries=visual.findall('./GeometrySet/objects/DisplayGeometry')
        for element in [visual]+geometries:
            if not np.array_equal(numbers(element.findtext('scale_factors')),np.ones(3)) or not np.array_equal(numbers(element.findtext('transform')),np.zeros(6)):
                raise ValueError('Unsupported donor geometry scaling/transform')
        if len(geometries)!=1:raise ValueError('Expected one source mesh per body')
        filename=geometries[0].findtext('geometry_file');path=GEOMETRY/filename
        if not path.exists():
            unsupported.append({'body':name,'missing_exact_geometry':filename,'decision':'No substitute or correspondence inferred'});continue
        identity=receipt(path);vertices,faces=surface(path)
        geom=convex_geometry_prior(vertices,faces)
        source=point(frames[name],geom['centroid_m'])
        target=point(prior['canonical_to_current_torso'],prior['bodies'][name]['geometry']['centroid_m'])
        correspondences.append({'body':name,'source_spine_m':source.tolist(),'target_torso_m':target.tolist()})
        geometry_records[name]={'source':identity,'convex_prior':geom,'provenance_limit':'Held OpenSim geometry library, exact filename match only; compatibility with MASI v2 mesh revision unverified; no inherited MIT license asserted'}
    transform,fit=proper_fit([e['source_spine_m'] for e in correspondences],[e['target_torso_m'] for e in correspondences])
    for record,error in zip(correspondences,fit['residuals_m']):record['residual_m']=error
    for name,record in geometry_records.items():
        record['body_origin_in_target_torso_m']=point(transform@frames[name],[0,0,0]).tolist()
    assembled={};roundtrip=[]
    for name in NAMES:
        t=transform@frames[name];local=transform_body(prior['bodies'][name]['current_torso_frame_prior'],np.linalg.inv(t))
        recovered=transform_body(local,t);roundtrip.append(recovered)
        joint=bodies[name].find('Joint')[0]
        assembled[name]={'name':'cervical_'+name,'source_body':name,'body_to_target_torso_reference':t.tolist(),
            'mass_properties_in_new_body_frame':local,'source_joint_xml':encoded(joint),
            'parent_body':'torso' if name=='cerv7' else 'cervical_'+joint.findtext('parent_body'),
            'root_parent_frame_to_target_torso':transform.tolist() if name=='cerv7' else None,
            'canonical_bone_ids':[r['id'] for r in prior['bodies'][name]['sources']],
            'canonical_geometry_to_new_body':(np.linalg.inv(t)@np.array(prior['canonical_to_current_torso'])).tolist()}
    combined=combine_bodies([prior['residual_torso']]+roundtrip)
    original=prior['current_torso']
    ledger={'mass_error_kg':combined['mass_kg']-original['mass_kg'],
        'com_error_m':float(np.linalg.norm(np.array(combined['center_m'])-original['center_m'])),
        'com_inertia_frobenius_error_kg_m2':float(np.linalg.norm(np.array(combined['inertia_kg_m2'])-original['inertia_kg_m2']))}
    if max(abs(v) for v in ledger.values())>1e-10:raise ValueError('Local-frame conversion failed mass/inertia conservation')
    muscles=[];external_counts={};point_count=0
    for muscle in model.findall('./ForceSet/objects/*'):
        if 'Muscle' not in muscle.tag:continue
        if muscle.tag!='Thelen2003Muscle' or muscle.findall('.//PathWrap') or muscle.findall('.//WrapObject'):raise ValueError('Unsupported donor muscle type/wrapping')
        points=[]
        for p in muscle.findall('./GeometryPath/PathPointSet/objects/*'):
            if p.tag!='PathPoint':raise ValueError('Unsupported moving/conditional attachment')
            owner=p.findtext('body');location=numbers(p.findtext('location'))
            if owner not in frames or location.shape!=(3,):raise ValueError('Unresolved muscle attachment')
            internal=owner in NAMES
            mapped=location if internal else point(transform@frames[owner],location)
            external_counts[owner]=external_counts.get(owner,0)+(not internal)
            points.append({'name':p.get('name'),'source_body':owner,'source_location_m':location.tolist(),
                'target_body':'cervical_'+owner if internal else 'torso','target_location_m':mapped.tolist(),
                'mapping':'Unchanged donor body attachment' if internal else 'Rigidly collapsed donor reference shoulder/spine frame; engineered fixed anchor'})
        if len(points)<2:raise ValueError('Missing muscle route')
        point_count+=len(points)
        muscles.append({'name':'cervical_'+muscle.get('name'),'source_xml':encoded(muscle),'path_points':points,
            'force_parameter_policy':'Retain donor parameters without strength or length scaling; calibration unresolved'})
    if len(muscles)!=78 or point_count!=198:raise ValueError('Unexpected source muscle inventory')
    result={'schema':'ihm.cervical-registration-recipe.v1','native_activation_allowed':False,
        'status':'Source-only assembled recipe; requires explicit native conversion and acceptance',
        'inputs':inputs,'target_model_identity':target_receipt,'inertial_prior_manifest_identity':inputs[0],
        'coordinate_conventions':{'donor':'meters, radians, Y superior, X anterior, Z right; body-fixed Z/X/Y cervical rotations',
            'target':'Frozen current torso body frame, common global registration retained',
            'fit':'One proper-rigid donor-spine to target-torso fit; no scale, no per-body joint repositioning',
            'reference':'All 24 neck coordinates explicitly zero, satisfying 18 cervical couplers. Shoulder coordinates frozen at source defaults; shoulder coupler assembly not claimed.'},
        'donor_spine_to_target_torso':transform.tolist(),'fit':fit,'correspondences':correspondences,
        'fit_evidence_scope':'Corresponding convex bone-envelope centroids are engineering proxies, not measured homologous landmarks or a validated anatomical atlas.',
        'unsupported_geometry':unsupported,'donor_geometry':geometry_records,
        'coordinates':{n:{'name':'cervical_'+n,'reference_value':0.,'source_xml':encoded(c)} for n,c in neck_coordinates.items()},
        'constraints':constraints,'bodies':assembled,'replace_target_torso_inertia':prior['residual_torso'],'conservation_check':ledger,
        'fixed_anchor_frames':{n:{'source_to_target_torso':(transform@frames[n]).tolist(),
            'source_default_coordinates':{c.get('name'):coordinates[c.get('name')] for c in bodies[n].findall('./Joint/*/CoordinateSet/objects/Coordinate')}} for n in ['spine']+ANCHORS},
        'muscles':muscles,'muscle_external_attachment_counts':{k:v for k,v in external_counts.items() if v},
        'recipe_inventory':{'new_bodies':9,'target_total_bodies':31,'new_muscles':78,'target_total_muscles':170,'new_coordinates':24,'new_constraints':18,'source_path_points':198},
        'required_before_native_acceptance':['Convert retained old XML components explicitly to target schema and unique names; resolve all frame sockets.',
            'Verify body-fixed transform convention with native coordinate perturbations, couplers and muscle moment arms.',
            'Resolve source C7 mesh absence and quantify surface alignment; centroid residual alone is insufficient.',
            'Calibrate fixed shoulder attachments or implement explicit shoulder mechanics; source defaults are not validated target anchors.',
            'Replace current torso-owned head/neck skin/contact and display ownership with these nine body mappings; avoid duplicate head geometry or contacts.',
            'Audit muscle overlap with existing 92 effectors, force-length operating points and strength; no donor whole-torso or shoulder mass is added.',
            'Recheck inertia physicality, total mass/COM/inertia, contact acceptance and supine equilibrium in the opt-in native model.'],
        'uncertainty':['Canonical global registration inherits 62.356mm RMS COM/envelope proxy residual.',
            'Held donor geometry filename matches do not authenticate MASI v2 geometry revision; mesh license/provenance needs resolution before redistribution or anatomical acceptance.',
            'Bone hull centroids and source neutral pose do not establish anatomical joint centers. Preserve residuals; no numerical correction of relative cervical joints.',
            'Fixed clavicle/scapula anchors omit shoulder motion and require experimental validation.',
            'Generic convex soft-tissue inertial priors retain all limitations of the bound prior manifest; no duplicate tissue or empirical density claim.']}
    output.mkdir(parents=True);(output/'inputs').mkdir()
    for digest,raw in captured.items():(output/'inputs'/f'{digest}.gz').write_bytes(gzip.compress(raw,mtime=0))
    (output/'recipe.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return {'path':str(output),'fit':fit,'conservation':ledger,'unsupported':unsupported,'inventory':result['recipe_inventory']}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True)
    print(json.dumps(build(parser.parse_args().output),indent=2))
