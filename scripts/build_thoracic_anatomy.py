"""Bounded source thorax recipe: retained geometry, candidate joints and mass debit."""
import argparse,gzip,hashlib,json,re
from pathlib import Path
import numpy as np
from ihm.assembly.thoracic_anatomy import surface_mass_prior,hinge_candidate,nearest_node_pair,descent_weights
from ihm.assembly.cervical_inertia import transform_body,partition_body,combine_bodies

ROOT=Path(__file__).resolve().parents[1]
PRIOR=ROOT/'data/research/cervical_inertia/v2/manifest.json'
ORDINALS=['first','second','third','fourth','fifth','sixth','seventh','eighth','ninth','tenth','eleventh','twelfth']


def build(output):
    output=Path(output).resolve()
    if output.exists() or not output.is_relative_to(ROOT):raise ValueError('Fresh workspace-owned output required')
    prior=json.loads(PRIOR.read_bytes());entry=next(x for x in prior['inputs'] if x['path'].endswith('/canonical_mechanics.json'))
    raw=gzip.decompress((PRIOR.parent/entry['retained_copy']).read_bytes())
    if hashlib.sha256(raw).hexdigest()!=entry['sha256']:raise ValueError('Frozen canonical identity mismatch')
    all_entities=json.loads(raw)['entities'];embedding=np.asarray(prior['canonical_to_current_torso'])
    entities={};captured={};geometry={};rows={};names={}
    def receipt(path,raw=None):
        raw=path.read_bytes() if raw is None else raw
        return {'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
    for e in all_entities:
        name=e['name'];kind=None
        if re.fullmatch(r'(left|right) \w+ rib',name):kind='rib'
        elif re.fullmatch(r'(left|right) \w+ costal cartilage',name):kind='cartilage'
        elif name in ['manubrium','body of sternum','xiphoid process']:kind='sternum'
        elif name=='diaphragm':kind='diaphragm'
        elif name in ['external intercostal muscle','internal intercostal muscle','innermost intercostal muscle']:kind='intercostal'
        elif re.fullmatch(r'\w+ thoracic vertebra',name) or name in [f'{s} lumbar vertebra' for s in ORDINALS[:3]]:kind='posterior_support'
        if kind is None:continue
        ref=e['reference_geometry'];path=ROOT/ref['path'];raw=path.read_bytes()
        if len(raw)>8000000 or hashlib.sha256(raw).hexdigest()!=ref['sha256'] or ref['units']!='m':raise ValueError('Geometry byte/unit mismatch')
        p=json.loads(gzip.decompress(raw))
        if p['display_decimation']:raise ValueError('Source geometry required')
        v=np.asarray(p['positions'],float).reshape(-1,3);f=np.asarray(p['indices'],int).reshape(-1,3)
        if len(v)>250000 or len(f)>500000:raise ValueError('Geometry exceeds bounded scope')
        local=v@embedding[:3,:3].T+embedding[:3,3]
        shell=surface_mass_prior(local,f,e['mass_kg'])
        ident=e['id'];entities[ident]=e;geometry[ident]=(v,f,local);names[name]=ident
        geometry_copy='geometry/'+path.name;captured[geometry_copy]=raw
        rows[ident]={'name':name,'kind':kind,'source_geometry':receipt(path,raw)|{'retained_copy':geometry_copy},
            'source_vertex_count':len(v),'source_face_count':len(f),'source_geometry_to_torso':embedding.tolist(),
            'source_proxy_mass_kg':e['mass_kg'],'source_mass_role':e.get('mass_role'),'source_volume_basis':e.get('volume_basis'),
            'torso_frame_mass_prior':shell,'mass_decision':'Retain inside parent torso' if kind=='posterior_support' else 'Debit exactly once from cervical-reduced torso; source proxy, not measured tissue mass',
            'source_material_parameters_used':False}
    counts={k:sum(r['kind']==k for r in rows.values()) for k in ['rib','cartilage','sternum','diaphragm','intercostal','posterior_support']}
    if counts!={'rib':24,'cartilage':14,'sternum':3,'diaphragm':1,'intercostal':6,'posterior_support':15}:raise ValueError('Unexpected anatomical inventory: '+str(counts))
    def binding(ident,vertex):
        f=geometry[ident][1];hits=np.argwhere(f==vertex)
        if len(hits)==0:raise ValueError('Attachment node has no source face')
        face,slot=map(int,hits[0]);bary=[0.,0.,0.];bary[slot]=1.
        return {'entity_id':ident,'vertex_index':int(vertex),'face_index':face,'barycentric':bary,'point_in_torso_m':geometry[ident][2][vertex].tolist()}
    def pair(a,b):
        p=nearest_node_pair(geometry[a][2],geometry[b][2])
        return {'a':binding(a,p['source_vertex']),'b':binding(b,p['target_vertex']),'gap_m':p['gap_m'],
            'basis':'Geometric nearest-node candidate; source contains no measured enthesis correspondence'}
    hinges={};attachments=[];joint_index={}
    for side in ['left','right']:
        for order in ORDINALS:
            rib=names[f'{side} {order} rib'];vertebra=names[f'{order} thoracic vertebra']
            h=hinge_candidate(geometry[rib][2],geometry[vertebra][2]);joint_index[rib]=len(hinges)
            body_frame=np.eye(4);body_frame[:3,3]=h['origin_m']
            rows[rib]['proposed_body_frame_to_torso']=body_frame.tolist()
            rows[rib]['proposed_body_local_mass_prior']=transform_body(rows[rib]['torso_frame_mass_prior'],np.linalg.inv(body_frame))
            hinges[rib]={'name':'thorax_'+rib,'coordinate_index':joint_index[rib],'joint_type':'PinJoint candidate','parent':'residual_torso',
                **h,'axis_conditioning_status':'ambiguous_geometric_axis' if h['first_to_second_axis_ratio']<2 else 'geometric_candidate_only','conditioning_ratio_threshold':2.,'posterior_support_entity':vertebra,'coordinate_units':'rad','reference_value':0.,
                'calibrated_range_rad':None,'positive_direction':'Deterministic PCA-axis sign; not classified as inspiratory',
                'attachment':pair(rib,vertebra)}
            cartilage=names.get(f'{side} {order} costal cartilage')
            if cartilage:attachments.append({'kind':'costochondral_candidate',**pair(rib,cartilage)})
    sternum_ids=[names[n] for n in ['manubrium','body of sternum','xiphoid process']]
    for a,b in zip(sternum_ids,sternum_ids[1:]):attachments.append({'kind':'sternum_piece_candidate',**pair(a,b)})
    for ident,row in rows.items():
        if row['kind']=='cartilage':
            candidates=[pair(ident,s) for s in sternum_ids]
            attachments.append({'kind':'cartilage_sternum_candidate',**min(candidates,key=lambda x:x['gap_m'])})
    # Intercostals are six aggregate surfaces, not individual numbered muscles.
    muscles=[]
    for ident,row in rows.items():
        if row['kind']!='intercostal':continue
        side='left' if geometry[ident][0][:,0].mean()>0 else 'right'
        for upper,lower in zip(ORDINALS,ORDINALS[1:]):
            a=names[f'{side} {upper} rib'];b=names[f'{side} {lower} rib']
            pa=pair(ident,a);pb=pair(ident,b)
            span=float(np.linalg.norm(np.array(pa['b']['point_in_torso_m'])-pb['b']['point_in_torso_m']))
            poor=max(pa['gap_m'],pb['gap_m'])>.5*span
            muscles.append({'source_entity':ident,'side_inferred_from_canonical_x':side,'rib_pair':[a,b],
                'upper_candidate':pa,'lower_candidate':pb,
                'candidate_span_m':span,'geometry_status':'unsupported_by_source_proximity' if poor else 'unvalidated_geometric_candidate',
                'proximity_screen':'Reject if either muscle-to-rib node gap exceeds half rib-to-rib candidate span; engineering screen, not measured enthesis tolerance',
                'fiber_path_status':'Two attachment candidates only; no fiber orientation, pennation, physiological cross section or optimal length identified',
                'active_force_parameters':None})
    diaphragm=names['diaphragm'];diaphragm_attachments=[]
    targets=[names[f'{side} {o} rib'] for side in ['left','right'] for o in ORDINALS[6:]]+[names['xiphoid process']]+[names[f'{o} lumbar vertebra'] for o in ORDINALS[:3]]
    for target in targets:diaphragm_attachments.append({'kind':'inferred_costal_sternal_crural_candidate',**pair(diaphragm,target)})
    anchors=sorted(set(x['a']['vertex_index'] for x in diaphragm_attachments))
    weights=descent_weights(geometry[diaphragm][2],anchors)
    weight_payload={'entity_id':diaphragm,'vertex_count':len(weights),'weights':weights.tolist(),
        'anchor_vertices':anchors,'axis_in_torso':(-embedding[:3,1]).tolist(),
        'formula':'Squared nearest-anchor-node distance normalized by maximum; zero at candidate anchors, max1. Geometry prior, not measured diaphragm fibers.',
        'motion':'x(q)=x_reference+weight*q_descent*axis; q_descent in meters. Moving rib-anchor coupling remains required.'}
    # Positive mass for the candidate diaphragm coordinate, from exact linear
    # nodal interpolation of w over each source triangle, not a generic1.5kg mode.
    v,f,_=geometry[diaphragm];tri=v[f];area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    w=weights[f];mean_w2=((w.sum(1))**2+(w*w).sum(1))/12
    effective=rows[diaphragm]['source_proxy_mass_kg']*float(area@mean_w2/area.sum())
    # All moved material shares enter this debit. Posterior supports stay owned
    # by the parent; never subtract them merely because their meshes were read.
    moving=[r['torso_frame_mass_prior'] for r in rows.values() if r['kind']!='posterior_support']
    reduced=prior['residual_torso'];remaining=partition_body(reduced,moving)
    total=combine_bodies([remaining]+moving)
    ledger={'debited_mass_kg':sum(b['mass_kg'] for b in moving),'parent_before_kg':reduced['mass_kg'],'parent_after_kg':remaining['mass_kg'],
        'mass_error_kg':total['mass_kg']-reduced['mass_kg'],
        'com_error_m':float(np.linalg.norm(np.asarray(total['center_m'])-reduced['center_m'])),
        'inertia_frobenius_error_kg_m2':float(np.linalg.norm(np.asarray(total['inertia_kg_m2'])-reduced['inertia_kg_m2']))}
    weight_bytes=gzip.compress(json.dumps(weight_payload,separators=(',',':')).encode(),mtime=0)
    inputs=[receipt(PRIOR),entry,receipt(Path(__file__).resolve()),receipt(ROOT/'ihm/assembly/thoracic_anatomy.py')]
    recipe={'schema':'ihm.thoracic-anatomy-recipe.v1','native_activation_allowed':False,'inputs':inputs,
        'target_model_identity':next(x for x in prior['inputs'] if x['path'].endswith('/native/assembled_model.osim')),
        'required_prior_partition':'Apply exact cervical_inertia/v2 residual ledger first; this is an alternative composed recipe, not a claim of current native ownership',
        'coordinate_frame':'Frozen initial actual native torso body frame; source meshes retained in canonical meters with one proper-rigid embedding',
        'source_license':'BodyParts3D ©DBCLS CC BY4.0; source/material priors remain differentiated',
        'counts':counts,'entities':rows,'rib_hinges':hinges,'attachment_candidates':attachments,
        'intercostal_candidates':muscles,'diaphragm_attachment_candidates':diaphragm_attachments,
        'diaphragm_mode':{'weights_file':'diaphragm_mode.json.gz','weights_sha256':hashlib.sha256(weight_bytes).hexdigest(),'coordinate_index':25,'coordinate_units':'m','effective_mass_kg':effective,
            'total_material_mass_kg':rows[diaphragm]['source_proxy_mass_kg'],'effective_mass_scope':'Derived kinetic coefficient for fixed-anchor descent alone, not extra mass; torso cross terms and moving-anchor modes not omitted by claiming closure'},
        'sternum_mode':{'coordinate_index':24,'coordinate_units':'m','axis_in_torso':embedding[:3,2].tolist(),'source_entities':sternum_ids,
            'kind':'Anterior translation candidate; cartilage deformation/closed-chain constraints unresolved'},
        'mass_partition':ledger,'replace_reduced_torso_with':remaining,
        'mechanical_recipe':{'coordinates':26,'rib_bodies':24,'sternum_compound_bodies':1,'diaphragm_kind':'Surface material mode requiring custom generalized mechanics; cannot substitute a free rigid mass',
            'cartilage_role':'Source-registered deformable links; no imported stiffness',
            'intercostal_role':'Source-registered material and candidate actuators; no imported active law',
            'loadable_next_step':'Create24PinJoint bodies and compound sternum material carrier from supplied local moments; convert diaphragm surface mode with consistent mass/Jacobian and close moving attachment constraints before coupling.',
            'no_extra_mass':True,'stiffness_or_damping':None,'muscle_activation_parameters':None},
        'limitations':['Source proxy masses inherit generic allocation and open-volume uncertainties; exact lamina moments do not make those masses measured.',
            'Nearest-node candidates are not verified anatomical insertions; report gaps and source indices, no hidden snapping.',
            'Rib PCA hinges are uncertain geometric hypotheses, including floating-rib joints; no physiological rotation sign/range inferred.',
            'Diaphragm fixed-anchor shape is a kinematic starting map; moving rib/abdominal support coupling and actual central tendon remain unresolved.',
            'Cavity V(q),J_V and native gas reference-volume calibration are not supplied by open tissue surfaces; no fake enclosed lung volume.',
            'No stiffness or recoil added; native bilateral chest-wall compliance and ideal pressure driver remain owners until explicitly transferred.',
            'No closed work ledger, native acceptance or patient-specific anatomical validation claimed.']}
    output.mkdir(parents=True);(output/'geometry').mkdir()
    for path,raw in captured.items():(output/path).write_bytes(raw)
    (output/'diaphragm_mode.json.gz').write_bytes(weight_bytes)
    (output/'manifest.json').write_text(json.dumps(recipe,indent=2,allow_nan=False)+'\n')
    return {'output':str(output),'counts':counts,'mass_partition':ledger,'diaphragm_effective_mass_kg':effective,
        'max_hinge_gap_m':max(h['attachment']['gap_m'] for h in hinges.values()),
        'max_diaphragm_candidate_gap_m':max(x['gap_m'] for x in diaphragm_attachments)}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True)
    print(json.dumps(build(p.parse_args().output),indent=2))
