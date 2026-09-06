#!/usr/bin/env python3
"""Preserve a shoulder donor complement and its dependencies, without registration."""
from copy import deepcopy
from collections import Counter
import hashlib,json
from pathlib import Path
import xml.etree.ElementTree as E
ROOT=Path(__file__).resolve().parents[1]
BASE='data/research/shoulder_complement'
SELECTED=('DELT1','DELT2','DELT3','SUPSP','INFSP','SUBSC','TMIN','TMAJ','PECM1','PECM2','PECM3','LAT1','LAT2','LAT3','CORB')
UNITS={'max_isometric_force':'N','optimal_fiber_length':'m','tendon_slack_length':'m','pennation_angle_at_optimal':'rad','max_contraction_velocity':'optimal_fiber_lengths/s','default_activation':'dimensionless','default_fiber_length':'m'}

def xml(e):return E.tostring(e,encoding='unicode')

def generate(root=ROOT):
 root=Path(root);base=root/BASE;acquisition=json.loads((base/'acquisition.json').read_text())
 for row in acquisition['sources']:
  b=(root/row['path']).read_bytes()
  if len(b)!=row['bytes'] or hashlib.sha256(b).hexdigest()!=row['sha256']:raise ValueError('Changed donor source: '+row['path'])
 receipt=next(x for x in acquisition['sources'] if x['path'].endswith('/MOBL_ARMS_41.osim'))
 tree=E.fromstring((root/receipt['path']).read_bytes());model=tree.find('Model')
 if model.findtext('length_units')!='meters' or model.findtext('force_units')!='N':raise ValueError('Unexpected donor units')
 muscles={x.get('name'):x for x in model.findall('./ForceSet/objects/*') if 'Muscle' in x.tag}
 bodies=model.findall('./BodySet/objects/Body');coordinates={x.get('name') for x in model.findall('.//Coordinate')}
 wraps={w.get('name'):(b.get('name'),w) for b in bodies for w in b.findall('./WrapObjectSet/objects/*')}
 force_doc=E.Element('OpenSimDocument',Version=tree.get('Version'));force_set=E.SubElement(force_doc,'ForceSet',name='unregistered_source_shoulder_complement');objects=E.SubElement(force_set,'objects')
 rows=[];point_types=Counter();referenced_wraps=set();attachment_bodies=set()
 for name in SELECTED:
  muscle=muscles[name];objects.append(deepcopy(muscle));points=[]
  for point in muscle.findall('.//PathPointSet/objects/*'):
   frame=point.findtext('socket_parent_frame');body=frame.split('/')[-1]
   if body not in {b.get('name') for b in bodies}:raise ValueError('Unresolved source path body')
   attachment_bodies.add(body);point_types[point.tag]+=1
   deps={x.text.strip().split('/')[-1] for x in point.iter() if 'coordinate' in x.tag and x.text and x.text.strip()}
   if not deps<=coordinates:raise ValueError('Unresolved source path coordinate')
   points.append({'name':point.get('name'),'type':point.tag,'source_frame':frame,
    'fixed_location_m':[float(x) for x in point.findtext('location').split()] if point.find('location') is not None else None,
    'coordinate_dependencies':sorted(deps),'complete_source_xml':xml(point)})
  path_wraps=[]
  for wrap in muscle.findall('.//PathWrap'):
   ref=wrap.findtext('wrap_object')
   if ref not in wraps:raise ValueError('Unresolved source wrap')
   referenced_wraps.add(ref);path_wraps.append({'name':wrap.get('name'),'source_wrap_name':ref,'source_body':wraps[ref][0],'method':wrap.findtext('method'),'range':wrap.findtext('range'),'complete_source_xml':xml(wrap)})
  rows.append({'name':name,'law_type':muscle.tag,'parameter_evidence':'exact supplied donor model values; population/model prior, not this body measurement',
   'parameters':{k:{'value':float(muscle.findtext(k)),'unit':unit} for k,unit in UNITS.items() if muscle.find(k) is not None},
   'path_points':points,'path_wraps':path_wraps,'complete_muscle_subtree_sha256':hashlib.sha256(E.tostring(muscle)).hexdigest()})
 dependencies=E.Element('SourceDependencies',native_installation='false',registration='unbound')
 # Full source sets preserve function coefficients, inertia, transforms and sockets.
 # Their presence is a donor inventory, not a request to add their masses to IHM.
 for tag in ['Ground','BodySet','JointSet','ConstraintSet']:dependencies.append(deepcopy(model.find(tag)))
 coordrows=[{'name':c.get('name'),'default_value':float(c.findtext('default_value','0')),'locked':c.findtext('locked'),'range':c.findtext('range'),'complete_source_xml':xml(c)} for c in model.findall('.//Coordinate')]
 couplers=[{'name':c.get('name'),'dependent':c.findtext('dependent_coordinate_name'),'independent':c.findtext('independent_coordinate_names'),'complete_source_xml':xml(c)} for c in model.findall('./ConstraintSet/objects/*')]
 report={'schema':'ihm.shoulder-complement-source.v1','donor':dict(receipt,document_version=tree.get('Version'),model_name=model.get('name'),muscle_count=len(muscles),body_count=len(bodies),credits=model.findtext('credits')),
  'acquisition_bytes':acquisition['new_source_bytes'],'mirror_equals_official_release_verified':False,'official_model_zip_held':False,
  'license':{'mirror_repository':'Apache-2.0 file retained; does not establish upstream rights','current_primary_project':'noncommercial use conditions plus BSD-3-clause wording; citation required','older_catalog':'MIT listed in older OpenSim model catalog; not relied upon','scope':'local research candidate; no public publication or deployment authorization inferred'},
  'units':{'position':'m','force':'N','mass':'kg','inertia':'kg m2','angle':'rad'},'muscles':rows,
  'path_point_type_counts':dict(point_types),'attachment_body_names':sorted(attachment_bodies),'referenced_wrap_names':sorted(referenced_wraps),
  'referenced_wraps':[{'name':name,'body':wraps[name][0],'type':wraps[name][1].tag,'complete_source_xml':xml(wraps[name][1])} for name in sorted(referenced_wraps)],
  'coordinates':coordrows,'couplers':couplers,
  'bodies':[{'name':b.get('name'),'mass_kg':float(b.findtext('mass')),'mass_center_m':b.findtext('mass_center'),'inertia_kg_m2':b.findtext('inertia')} for b in bodies],
  'unrepresented_girdle_muscles':['trapezius','serratus anterior','rhomboids','levator scapulae','pectoralis minor'],
  'registration':{'body_transforms':None,'coordinate_mapping':None,'bilateral_reflection':None,'girdle_mass_partition':None,'added_native_mass_kg':0,
   'required_target_owners':['existing torso residual after cervical/thoracic partition','exclusive clavicle/scapula per side','existing humerus per side'],
   'gates':['resolve donor provenance and reuse scope','register anatomical girdle and humeral frames','coherent bilateral coordinate/coupler mapping','exclusive mass/COM/inertia partition with cervical/thoracic candidate','preserve moving/conditional path functions and source wrap ranges','reconcile existing Arm26 shoulder origins without duplicated muscle force','source-range geometry and tension-feasible moment-arm validation in native engine']},
  'native_installation':False,'native_load_test_performed':False,'donor_geometry_meshes_held':False,'parameter_tuning_performed':False,
  'source_primary_basis':{'species':'Homo sapiens','anthropometry':'50th-percentile male model','force_generation':'Saul2015 reports32muscle volumes and joint moments from5healthy young male subjects; fitted specific tension50.8N/cm2','uncertainty':'No per-parameter covariance or individual donor sample distribution supplied here; no target calibration inferred'}}
 return report,E.tostring(force_doc,encoding='unicode'),E.tostring(dependencies,encoding='unicode')

if __name__=='__main__':
 report,forces,dependencies=generate();base=ROOT/BASE
 (ROOT/'data/sources/shoulder_complement.json').write_text(json.dumps(report,indent=2)+'\n')
 (base/'shoulder_forces.xml').write_text(forces+'\n');(base/'source_dependencies.xml').write_text(dependencies+'\n')
 print(json.dumps({'selected_muscles':len(report['muscles']),'wraps':len(report['referenced_wrap_names']),'point_types':report['path_point_type_counts'],'source_bytes':report['acquisition_bytes']}))
