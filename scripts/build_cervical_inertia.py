"""Materialize opt-in inertial inputs from bounded retained geometry; no OpenSim."""
import argparse,gzip,hashlib,json,sys
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from ihm.assembly.cervical_inertia import convex_geometry_prior,segment_prior,homogeneous_prior,transform_body,partition_body

ROOT=Path(__file__).resolve().parents[1]
REFERENCE='data/derived/audits/cutaneous-factory-ylro2d66/body/mechanics'
CERVICAL={'cerv1':'FJ3176','cerv2':'FJ3177','cerv3':'FJ3161','cerv4':'FJ3164','cerv5':'FJ3167','cerv6':'FJ3170','cerv7':'FJ3172'}
CRANIUM=['FJ3199','FJ3200','FJ3269','FJ3273','FJ3274','FJ3281','FJ3287','FJ3309','FJ3375','FJ3379','FJ3380','FJ3386','FJ3392','FJ3394','FJ3395']


def identity(path):
    raw=path.read_bytes()
    return {'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}


def build(output):
    output=Path(output).resolve()
    if output.exists() or not output.is_relative_to(ROOT):raise ValueError('Fresh workspace-owned output required')
    reference=ROOT/REFERENCE
    paths=[reference/'canonical_mechanics.json',reference/'registration.json',reference/'identity.json',reference/'native/assembled_model.osim',
        reference/'native/inputs/subject_walk_scaled.osim',ROOT/'data/research/cervical/MASI_HMaleMuscle_HMaleMassDistr.osim',
        ROOT/'ihm/assembly/cervical_inertia.py',Path(__file__).resolve()]
    snapshots=[p.read_bytes() for p in paths]
    identities=[identity(p)|{'retained_copy':f'inputs/{index:02d}-{p.name}.gz','compression':'gzip'} for index,p in enumerate(paths)]
    if any(hashlib.sha256(raw).hexdigest()!=entry['sha256'] for raw,entry in zip(snapshots,identities)):raise ValueError('Inputs changed during capture')
    specs={e['id']:e for e in json.loads(paths[0].read_bytes())['entities']}
    registration=json.loads(paths[1].read_bytes())
    target=ET.parse(paths[3]).getroot().find('Model');torso=target.find('./BodySet/objects/Body[@name="torso"]')
    if len(target.findall('./BodySet/objects/Body'))!=22 or sum('Muscle' in f.tag for f in target.findall('./ForceSet/objects/*'))!=92:raise ValueError('Expected frozen 22-body/92-muscle reference')
    declared=json.loads(paths[2].read_bytes())
    if declared['canonical_mechanics_sha256']!=identities[0]['sha256']:raise ValueError('Reference mechanics identity mismatch')
    transform=np.array(registration['groups']['torso']['canonical_reference_embedding_in_source_body'])
    # Validate the embedding is from the frozen initial registration itself.
    recomputed=np.linalg.inv(np.array(registration['groups']['torso']['native_reference_transform']))@np.linalg.inv(np.array(registration['global_rigid_fit']['source_to_canonical_ground']))
    if not np.allclose(transform,recomputed,atol=1e-12):raise ValueError('Inconsistent canonical/native initial frame')
    xx,yy,zz,xy,xz,yz=map(float,torso.findtext('inertia').split())
    parent={'mass_kg':float(torso.findtext('mass')),'center_m':list(map(float,torso.findtext('mass_center').split())),
            'inertia_kg_m2':[[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]]}
    donor=ET.parse(paths[5]).getroot().find('Model')
    masses={b.get('name'):float(b.findtext('mass')) for b in donor.findall('./BodySet/objects/Body')}
    groups={name:[ident] for name,ident in CERVICAL.items()}|{'skull':CRANIUM,'jaw':['FJ3289']}
    geometries={};captured={};records={};bodies=[]
    for name,ids in groups.items():
        vertices=[];faces=[];offset=0;sources=[]
        for ident in ids:
            spec=specs['body-bp3d-'+ident];ref=spec['reference_geometry'];path=ROOT/ref['path'];raw=path.read_bytes()
            if len(raw)>2000000 or hashlib.sha256(raw).hexdigest()!=ref['sha256']:raise ValueError('Geometry size or identity mismatch: '+ident)
            payload=json.loads(gzip.decompress(raw))
            if payload['display_decimation']:raise ValueError('Full source geometry required')
            v=np.asarray(payload['positions'],dtype=float).reshape(-1,3);f=np.asarray(payload['indices'],dtype=int).reshape(-1,3)
            if len(v)>50000 or len(f)>100000:raise ValueError('Bounded geometry limit exceeded')
            vertices.append(v);faces.append(f+offset);offset+=len(v)
            sources.append({'id':spec['id'],'name':spec['name'],**identity(path)})
            captured[path.name]=raw
        geometry=convex_geometry_prior(np.vstack(vertices),np.vstack(faces));geometries[name]=geometry
        body=segment_prior(geometry,masses[name]) if name.startswith('cerv') else homogeneous_prior(geometry,masses[name])
        converted=transform_body(body,transform);converted['name']=name;bodies.append(converted)
        original=donor.find('./BodySet/objects/Body[@name="'+name+'"]')
        records[name]={'geometry':geometry,'sources':sources,'canonical_body_prior':body,'current_torso_frame_prior':converted,
                       'unchanged_donor_mass_kg':masses[name],
                       'original_donor_inertia_kg_m2':[float(original.findtext('inertia_'+k)) for k in ('xx','yy','zz','xy','xz','yz')],
                       'decision':'Replace donor tensor with geometry-conditioned prior; no eigenvalue clipping and no Body.cpp fallback'}
    try:residual=partition_body(parent,bodies);physical=True;partition_error=None
    except ValueError as error:residual=None;physical=False;partition_error=str(error)
    sensitivities=[]
    for bone_density in (1500.,2200.):
        for soft_density in (900.,1100.):
            scenario=[transform_body(segment_prior(geometries[name],masses[name],bone_density_kg_m3=bone_density,soft_density_kg_m3=soft_density),transform)
                      for name in CERVICAL]+[records[name]['current_torso_frame_prior'] for name in ('skull','jaw')]
            try:
                rem=partition_body(parent,scenario);summary={'physical':True,'torso_second_moment_min_kg_m2':min(rem['second_moment_eigenvalues_kg_m2'])}
            except ValueError as error:summary={'physical':False,'error':str(error)}
            sensitivities.append({'bone_density_kg_m3':bone_density,'soft_density_kg_m3':soft_density,**summary})
    overlap=[]
    for index,(a,_) in enumerate(CERVICAL.items()):
        for b in list(CERVICAL)[index+1:]:
            ca=records[a]['canonical_body_prior']['components'];cb=records[b]['canonical_body_prior']['components']
            lower=np.maximum(ca['envelope_bounds_min_m'],cb['envelope_bounds_min_m']);upper=np.minimum(ca['envelope_bounds_max_m'],cb['envelope_bounds_max_m'])
            volume=float(np.prod(np.maximum(upper-lower,0)))
            if volume:overlap.append({'a':a,'b':b,'aabb_overlap_upper_bound_m3':volume})
    result={'schema':'ihm.cervical-inertial-prior.v1','inputs':identities,'reference_kind':'Frozen initial 22-body/92-muscle factory registration, no moving or settled pose',
        'canonical_to_current_torso':transform.tolist(),'registration_uncertainty':registration['global_rigid_fit'],
        'current_torso':parent,'donor_mass_scaling':'Original MASI neck/head masses retained; debited from actual mass-scaled torso, not added to total mass',
        'bodies':records,'residual_torso':residual,'physical_partition':physical,'partition_error':partition_error,
        'density_sensitivity':sensitivities,'envelope_overlap_bounds':overlap,
        'limitations':['Convex hull fills source concavities and canals; source closed-volume validity is not asserted.',
            '1900/1000kg/m3 bone/soft densities are generic engineering priors; intervals1500–2200/900–1100 are sensitivity scenarios, not empirical confidence.',
            'Segment soft shells can overlap. They define positive lumped mass distributions, not a certified spatial segmentation or continuum density field.',
            'Head composite is homogeneous inside a craniofacial bone hull; skull/brain/scalp density heterogeneity and soft outer surface are unresolved.',
            'Jaw convex hull fills its open mouth concavity; source0.2kg mass retained as a separate proxy.',
            'Global reference registration has62.356mm RMS landmark residual; inertia feasibility does not validate cervical joint centers.',
            'No native model, inertia fallback, activation, equilibrium or collision run performed.'],
        'opt_in_augmentation_recipe':{'status':'inertial_inputs_only' if physical else 'blocked_by_torso_inertia',
            'native_activation_allowed':False,'proposed_body_count':31,'masi_proposed_muscle_count':170,
            'new_body_inertias_in_current_torso_frame':bodies,'replace_current_torso_with':residual,
            'required_before_native':['Register donor cervical joint root and each source body frame to these common-frame priors; express each tensor about its COM in that body frame.',
                'Retain source24coordinates/18couplers, all78MASI paths with unique identities; no duplicate whole-body/shoulder inertia.',
                'Bind source clavicle/scapula/thoracic attachment frames explicitly; do not substitute uncalibrated parent names.',
                'Reassign skull/cervical skin/contact ownership from torso, retain one world transform, and rebuild posed surface contact.',
                'Verify instantiated native tensors exactly match the proposed tensors and reject any automatic inertia repair.']}}
    # Compute and validate everything before publishing the fresh artifact.
    output.mkdir(parents=True);(output/'geometry').mkdir();(output/'inputs').mkdir()
    for entry,raw in zip(identities,snapshots):(output/entry['retained_copy']).write_bytes(gzip.compress(raw,mtime=0))
    for name,raw in captured.items():(output/'geometry'/name).write_bytes(raw)
    (output/'manifest.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


def verify(manifest_path):
    manifest_path=Path(manifest_path);output=manifest_path.parent;result=json.loads(manifest_path.read_bytes())
    for entry in result['inputs']:
        raw=gzip.decompress((output/entry['retained_copy']).read_bytes())
        if len(raw)!=entry['bytes'] or hashlib.sha256(raw).hexdigest()!=entry['sha256']:raise ValueError('Retained input identity changed')
    for name,record in result['bodies'].items():
        vertices=[];faces=[];offset=0
        for source in record['sources']:
            raw=(output/'geometry'/Path(source['path']).name).read_bytes()
            if len(raw)!=source['bytes'] or hashlib.sha256(raw).hexdigest()!=source['sha256']:raise ValueError('Retained geometry identity changed')
            payload=json.loads(gzip.decompress(raw));v=np.asarray(payload['positions']).reshape(-1,3);f=np.asarray(payload['indices']).reshape(-1,3)
            vertices.append(v);faces.append(f+offset);offset+=len(v)
        geometry=convex_geometry_prior(np.vstack(vertices),np.vstack(faces))
        if geometry!=record['geometry']:raise ValueError('Geometry moments or approximation audit changed')
        body=segment_prior(geometry,record['unchanged_donor_mass_kg']) if name.startswith('cerv') else homogeneous_prior(geometry,record['unchanged_donor_mass_kg'])
        reconstructed=transform_body(body,result['canonical_to_current_torso'])
        for key in ('mass_kg','center_m','inertia_kg_m2'):
            if not np.allclose(reconstructed[key],record['current_torso_frame_prior'][key],rtol=1e-12,atol=1e-14):raise ValueError('Materialized body prior changed')
    residual=partition_body(result['current_torso'],[b['current_torso_frame_prior'] for b in result['bodies'].values()])
    if residual!=result['residual_torso']:raise ValueError('Exclusive torso partition changed')
    if result['opt_in_augmentation_recipe']['native_activation_allowed']:raise ValueError('Inertial recipe cannot authorize native activation')
    return {'verified':True,'retained_input_count':len(result['inputs']),'body_count':len(result['bodies']),'mass_conservation_residual_kg':residual['reconstruction_residual']['mass_kg'],
            'inertia_conservation_residual_kg_m2':residual['reconstruction_residual']['inertia_frobenius_kg_m2'],'native_runs':0}


if __name__=='__main__':
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True);group.add_argument('--output');group.add_argument('--verify');args=parser.parse_args()
    if args.verify:print(json.dumps(verify(ROOT/args.verify)));sys.exit(0)
    result=build(ROOT/args.output)
    print(json.dumps({'output':args.output,'physical_partition':result['physical_partition'],'error':result['partition_error'],
        'new_segment_mass_kg':sum(b['unchanged_donor_mass_kg'] for b in result['bodies'].values()),
        'residual_torso':result['residual_torso'],'density_sensitivity':result['density_sensitivity'],'native_runs':0}))
