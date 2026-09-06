#!/usr/bin/env python3
"""Freeze a candidate-root manifest and registration, without native execution."""
from copy import deepcopy
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from prepare_skin_layer_migration import ANATOMY, MECHANICS, EVIDENCE, digest, encode
from ihm.assembly.body import build as build_body, CanonicalBody
from ihm.assembly.respiration import build as build_respiration
from ihm.assembly.articulated import CanonicalRegistration
from ihm.assembly.skin_layers import physical_skin_support

RUNTIME_NAMES=('body','body_protocol','body_runtime','body_states','contracts','interfaces','cosimulation','evidence','mechanics','brain','respiration','peripheral','certainty','temporal')
REFERENCE='data/derived/supine-support-5ma720yd/initial_native.json'


def support_equivalence(old, new, evidence):
    """Mask receipt stays historical; all geometry-bearing inputs must be equal."""
    oldids={e['id']:e for e in old['entities']};newids={e['id']:e for e in new['entities']}
    receipts=evidence['source_receipts']
    for receipt in receipts:
        ident=receipt['entity_id'];a=oldids[ident];b=newids[ident]
        for key in ('reference_geometry','centroid_m','bounds_m'):
            if a[key]!=b[key]:raise ValueError('Mask geometric dependency changed: '+ident)
        if b['reference_geometry']['sha256']!=receipt['sha256'] or b['reference_geometry']['path']!=receipt['path']:
            raise ValueError('Mask geometry receipt differs: '+ident)
    return {'accepted_for':'historical geometric component-mask reuse only',
        'source_entities':[r['entity_id'] for r in receipts],
        'geometry_metadata_equal':True,'anatomy_epoch_equal':False,
        'historical_evidence_rewritten':False,
        'scope':'Retained geometry references, bounds and centroids equal; exact skin bytes separately checked; other source bytes not rescanned'}


def entity_diff(old,new):
    a={e['id']:e for e in old['entities']};b={e['id']:e for e in new['entities']}
    if set(a)!=set(b):raise ValueError('Entity IDs changed')
    return [{'id':key,'fields':sorted(k for k in set(a[key])|set(b[key]) if a[key].get(k)!=b[key].get(k))}
            for key in a if a[key]!=b[key]]


def stage(root,candidate,output):
    root=Path(root).resolve();candidate=Path(candidate).resolve();output=Path(output).resolve()
    if output.exists() or not output.is_relative_to(root/'data/derived') or output.is_relative_to(root/'data/derived/canonical'):
        raise ValueError('Fresh derived epoch outside canonical required')
    manifest=json.loads((candidate/'manifest.json').read_bytes())
    inputs={p:(root/p).read_bytes() for p in manifest['inputs']}
    for p,raw in inputs.items():
        if digest(raw)!=manifest['inputs'][p]:raise ValueError('Migration input changed: '+p)
    migrated={name:(candidate/name).read_bytes() for name in manifest['outputs']}
    for name,raw in migrated.items():
        if digest(raw)!=manifest['outputs'][name]:raise ValueError('Migration candidate changed: '+name)
    olda=json.loads(inputs[ANATOMY]);olda_mech=json.loads(inputs[MECHANICS]);newa=json.loads(migrated['anatomy.json']);newm=json.loads(migrated['mechanics.json'])
    equivalent=support_equivalence(olda,newa,json.loads(inputs[EVIDENCE]))
    old_entities={e['id']:e for e in olda['entities']}
    for entity in newa['entities']:
        for key in ('reference_geometry','bounds_m','centroid_m','principal_axis'):
            if entity.get(key)!=old_entities[entity['id']].get(key):raise ValueError('Source geometry metadata changed: '+entity['id'])
    skin,=[e for e in newa['entities'] if e['role']=='skin']
    physical_skin_support(skin['reference_geometry'],inputs[skin['reference_geometry']['path']],inputs[EVIDENCE])
    captured={p:(root/p).read_bytes() for p in [REFERENCE,*('ihm/assembly/'+n+'.py' for n in RUNTIME_NAMES),
        *('data/derived/canonical/'+n+'.json' for n in ('profile','brain','peripheral','respiration','body'))]}
    profile=json.loads(captured['data/derived/canonical/profile.json']);patient=profile['native_patient_path'];captured[patient]=(root/patient).read_bytes()
    if digest(captured[patient])!=profile['native_patient_sha256']:raise ValueError('Patient receipt differs')
    # This staging deliberately does not instantiate ArticulatedBodyPlant/native.
    reference=json.loads(captured[REFERENCE])
    oldreg=CanonicalRegistration(olda_mech,reference).manifest();newreg=CanonicalRegistration(newm,reference).manifest()
    if oldreg!=newreg:raise ValueError('Skin metadata migration unexpectedly changes registration')
    output.mkdir(parents=True);staged=output/'root';canonical=staged/'data/derived/canonical'
    def write(path,raw):
        target=staged/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
    for path,raw in captured.items():
        if path not in ('data/derived/canonical/body.json','data/derived/canonical/respiration.json'):write(path,raw)
    for name,raw in migrated.items():write('data/derived/canonical/'+name,raw)
    write(EVIDENCE,inputs[EVIDENCE]);write(skin['reference_geometry']['path'],inputs[skin['reference_geometry']['path']])
    respiration=build_respiration(staged);write('data/derived/canonical/respiration.json',encode(respiration))
    oldresp=json.loads(captured['data/derived/canonical/respiration.json'])
    respfields=sorted(k for k in set(oldresp)|set(respiration) if oldresp.get(k)!=respiration.get(k))
    if respfields!=['anatomy_sha256']:raise ValueError('Respiration changed beyond anatomy receipt')
    body=build_body(staged);loaded=CanonicalBody.from_workspace(staged)
    (output/'registration.json').write_bytes(encode(newreg))
    (output/'migration_manifest.json').write_bytes((candidate/'manifest.json').read_bytes())
    oldbody=json.loads(captured['data/derived/canonical/body.json'])
    top=lambda a,b:sorted(k for k in set(a)|set(b) if a.get(k)!=b.get(k))
    diffs={'anatomy_entities':entity_diff(olda,newa),'anatomy_top_level':top(olda,newa),
        'mechanics_entities':entity_diff(olda_mech,newm),'mechanics_top_level':top(olda_mech,newm),
        'changed_link_damping_count':sum(a['damping_ns_m']!=b['damping_ns_m'] for a,b in zip(olda_mech['links'],newm['links'])),
        'respiration_fields':respfields,'body_manifest_fields':top(oldbody,body),'registration_equal':True}
    (output/'differences.json').write_bytes(encode(diffs))
    checks={
        'body_manifest':{'status':'passed','consumer':'CanonicalBody.from_workspace','entity_count':len(loaded.entities)},
        'articulated_registration':{'status':'passed','result':'Recomputed old/new against same retained native reference; exactly equal'},
        'supine_support':{'status':'not_rebuilt','reason':'Current native equilibrium keeps its own input epoch; shell thickness available for next owner'},
        'territory_mask':{'status':'geometry_equivalent','receipt':equivalent},
        'respiration_systemic_projection':{'status':'respiration_rebuilt','result':'Only anatomy hash changes; systemic/native trajectory not run'},
        'hair_manifest':{'status':'historical_receipt_not_rewritten','reason':'All entity geometry references preserved; future full hair artifact must receive its own equivalence receipt'},
        'material_domains':{'status':'not_rebuilt','reason':'Canonical proxy allocation changed; physical domain ownership/normalization needs dedicated builder acceptance'},
        'body_runtime_touch':{'status':'new_body_manifest_loads','reason':'No state registry/touch materialization or runtime stepping performed'},
        'viewer_geometry':{'status':'partial_staging','reason':'Only skin geometry copied; remaining geometry references preserved and served from original source only after explicit publication'},
    }
    hashes={str(p.relative_to(output)):digest(p.read_bytes()) for p in output.rglob('*') if p.is_file()}
    report={'schema':'ihm.skin-candidate-staging.v1','published':False,'candidate_root':str(staged.relative_to(root)),
        'inputs':{**manifest['inputs'],**{p:digest(raw) for p,raw in captured.items()}},'outputs':hashes,
        'implementation_sources':{str(Path(__file__).relative_to(root)):digest(Path(__file__).read_bytes()),'ihm/assembly/articulated.py':digest((root/'ihm/assembly/articulated.py').read_bytes())},
        'acceptance':checks,'native_executed':False,'native_mass_migrated':False,
        'scope':'Metadata-only candidate root and retained-reference registration; no native build, launch, full geometry scan, or viewer deployment'}
    # Detect source changes during this bounded build rather than accepting mixed epochs.
    for p,raw in {**inputs,**captured}.items():
        if (root/p).read_bytes()!=raw:raise ValueError('Source changed during staging: '+p)
    (output/'acceptance.json').write_bytes(encode(report))
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);parser.add_argument('--candidate',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    report=stage(args.root,args.candidate,args.output);print(json.dumps({'candidate_root':report['candidate_root'],'body_load':report['acceptance']['body_manifest'],'native_mass_migrated':False}))
