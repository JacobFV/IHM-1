"""Prepare source-only moving material maps and an engineered lung envelope."""
import argparse,gzip,hashlib,io,json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree,ConvexHull
from ihm.assembly.thoracic_mechanism import closed_volume

ROOT=Path(__file__).resolve().parents[1]
ANATOMY=ROOT/'data/research/thoracic_anatomy/v3/manifest.json'


def packed(**arrays):
    b=io.BytesIO();np.savez_compressed(b,**arrays);return b.getvalue()
def receipt(path):
    raw=path.read_bytes();return {'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
def cardinal(points,anchors):
    k=min(3,len(anchors));distance,index=cKDTree(anchors).query(points,k=k)
    if k==1:distance=distance[:,None];index=index[:,None]
    weights=1/np.maximum(distance,1e-12)**2;weights/=weights.sum(1)[:,None]
    exact=distance[:,0]<1e-12;weights[exact]=0;weights[exact,0]=1
    return index,weights


def build(output):
    output=Path(output).resolve()
    if output.exists() or not output.is_relative_to(ROOT):raise ValueError('Fresh workspace artifact required')
    anatomy=json.loads(ANATOMY.read_bytes());entities=anatomy['entities'];geometry={};drivers={}
    for ident,h in anatomy['rib_hinges'].items():drivers[ident]=h['coordinate_index']
    for ident in anatomy['sternum_mode']['source_entities']:drivers[ident]=24
    for ident,row in entities.items():
        if row['kind']=='posterior_support':continue
        raw=(ANATOMY.parent/row['source_geometry']['retained_copy']).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=row['source_geometry']['sha256']:raise ValueError('Anatomy geometry identity mismatch')
        p=json.loads(gzip.decompress(raw));v=np.asarray(p['positions'],float).reshape(-1,3)
        t=np.asarray(row['source_geometry_to_torso']);geometry[ident]=v@t[:3,:3].T+t[:3,3]
    maps={};materials={};audit=[]
    raw=(ANATOMY.parent/anatomy['diaphragm_mode']['weights_file']).read_bytes();diaphragm=json.loads(gzip.decompress(raw))
    for ident,v in geometry.items():
        row=entities[ident];candidates=[]
        if row['kind'] in ['rib','sternum']:
            arrays={};entry={'map_kind':'rigid','driver':drivers[ident]}
        else:
            if row['kind']=='cartilage':
                for pair in anatomy['attachment_candidates']:
                    if pair['a']['entity_id']==ident and pair['b']['entity_id'] in drivers:candidates.append((pair['a'],pair['b'],pair['gap_m']))
                    elif pair['b']['entity_id']==ident and pair['a']['entity_id'] in drivers:candidates.append((pair['b'],pair['a'],pair['gap_m']))
            elif row['kind']=='intercostal':
                for route in anatomy['intercostal_candidates']:
                    if route['source_entity']!=ident:continue
                    if route['geometry_status']=='unsupported_by_source_proximity':continue
                    for key in ['upper_candidate','lower_candidate']:
                        pair=route[key];candidates.append((pair['a'],pair['b'],pair['gap_m']))
            elif row['kind']=='diaphragm':
                for pair in anatomy['diaphragm_attachment_candidates']:candidates.append((pair['a'],pair['b'],pair['gap_m']))
            else:raise ValueError('Unowned moving material')
            unique={};discard=[]
            for source,target,gap in sorted(candidates,key=lambda c:c[2]):
                vertex=source['vertex_index']
                if vertex in unique:discard.append({'source_vertex':vertex,'target_entity':target['entity_id'],'gap_m':gap});continue
                unique[vertex]=(source,target,gap)
            if len(unique)<2:raise ValueError('Insufficient distinct source anchors: '+ident)
            nodes=np.array(list(unique));reference=v[nodes]
            codes=np.array([drivers.get(target['entity_id'],-1) for source,target,gap in unique.values()])
            indices,weights=cardinal(v,reference)
            arrays={'anchor_reference':reference,'anchor_driver':codes,'anchor_source_vertices':nodes,'indices':indices,'weights':weights}
            if row['kind']=='diaphragm':arrays['descent_weights']=np.asarray(diaphragm['weights'])
            entry={'map_kind':'moving_anchors','anchor_count':len(unique),'source_anchor_nodes':nodes.tolist(),
                'target_owners':[target['entity_id'] for source,target,gap in unique.values()],
                'reference_gap_m':[gap for source,target,gap in unique.values()],
                'discarded_duplicate_node_candidates':discard,
                'interpolation':'Cardinal three-nearest-anchor inverse-square displacement interpolation; source-near attachment offsets rotate with target rigid owner',
                'scope':'Engineering material map; source topology retained, no measured fiber or strain calibration'}
            audit.append({'entity':ident,'anchors':len(unique),'maximum_node_to_anchor_distance_m':float(cKDTree(reference).query(v)[0].max())})
        blob=packed(**arrays);filename='maps/'+ident+'.npz';maps[filename]=blob
        materials[ident]=entry|{'map_file':filename,'map_sha256':hashlib.sha256(blob).hexdigest(),'mass_kg':row['source_proxy_mass_kg']}
    # The convex solid is only an engineered cavity geometry. It is NOT native
    # lung gas volume, not a tissue mass, and includes the interlobar/mediastinal gap.
    prior=json.loads((ROOT/'data/research/cervical_inertia/v2/manifest.json').read_bytes())
    original=next(e for e in prior['inputs'] if e['path'].endswith('/canonical_mechanics.json'))
    raw=gzip.decompress((ROOT/'data/research/cervical_inertia/v2'/original['retained_copy']).read_bytes())
    if hashlib.sha256(raw).hexdigest()!=original['sha256']:raise ValueError('Frozen anatomy source changed')
    canonical=json.loads(raw);lung=[];lung_receipts=[];lung_archives={};t=np.asarray(prior['canonical_to_current_torso'])
    for e in canonical['entities']:
        if not ('lobe of' in e['name'] and 'lung' in e['name']):continue
        path=ROOT/e['reference_geometry']['path'];raw=path.read_bytes()
        if hashlib.sha256(raw).hexdigest()!=e['reference_geometry']['sha256'] or len(raw)>8000000:raise ValueError('Lung source size/identity mismatch')
        p=json.loads(gzip.decompress(raw))
        if 'display_decimation' in p:
            if p['display_decimation']:raise ValueError('Decimated lung geometry unsupported')
        elif not (p.get('registration_id')=='z_anatomy' and p.get('units')=='m' and len(p.get('source_geometry_sha256',''))==64):
            raise ValueError('Unknown lung geometry receipt schema')
        v=np.asarray(p['positions'],float).reshape(-1,3);lung.append(v@t[:3,:3].T+t[:3,3])
        filename='lung_sources/'+path.name;lung_archives[filename]=raw
        lung_receipts.append(receipt(path)|{'retained_copy':filename,'entity_id':e['id'],'name':e['name'],'registered_source_geometry_sha256':p.get('source_geometry_sha256'),'source_schema':'Registered Z-Anatomy topology; this snapshot has no display-decimation field' if 'display_decimation' not in p else 'Canonical topology with explicit decimation flag'})
    if len(lung)!=5:raise ValueError('Expected five retained lobe geometries')
    cloud=np.vstack(lung);hull=ConvexHull(cloud-cloud.mean(0));faces=hull.simplices.copy()
    tri=cloud[faces];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);wrong=np.einsum('ij,ij->i',normal,hull.equations[:,:3])<0
    faces[wrong]=faces[wrong][:,[0,2,1]];keep=np.unique(faces);lookup=np.full(len(cloud),-1,int);lookup[keep]=np.arange(len(keep));cavity=cloud[keep];faces=lookup[faces]
    volume,_=closed_volume(cavity,faces)
    if not np.isclose(volume,hull.volume,rtol=1e-10):raise ValueError('Closed cavity integration mismatch')
    owner_ids=[ident for ident,row in entities.items() if row['kind'] in ['rib','sternum','diaphragm']]
    all_nodes=np.vstack([geometry[i] for i in owner_ids]);node_codes=np.concatenate([np.full(len(geometry[i]),code,int) for code,i in enumerate(owner_ids)])
    node_indices=np.concatenate([np.arange(len(geometry[i])) for i in owner_ids])
    indices,weights=cardinal(cavity,all_nodes)
    cavity_blob=packed(reference=cavity,faces=faces,entity_codes=node_codes[indices],vertex_indices=node_indices[indices],weights=weights)
    locked=sorted(h['coordinate_index'] for h in anatomy['rib_hinges'].values() if h['axis_conditioning_status']=='ambiguous_geometric_axis')
    manifest={'schema':'ihm.thoracic-mechanism.v1','native_activation_allowed':False,'anatomy_manifest':receipt(ANATOMY),
        'code_inputs':[receipt(Path(__file__).resolve()),receipt(ROOT/'ihm/assembly/thoracic_mechanism.py')],
        'materials':materials,'moving_material_count':len(materials),'locked_internal_coordinates':locked,
        'independent_generalized_coordinates':32-len(locked),'coordinate_convention':'v_torso(3),omega_torso(3),24ribangles(rad),sternum_anterior(m),diaphragm_inferior(m); all spatial values in frozen torso axes',
        'engineering_domain':{'rib_angles_abs_rad':.05,'sternum_abs_m':.005,'diaphragm_abs_m':.02,'scope':'Declared kinematic exploration bound, not physiological ROM or stiffness calibration; reject local cavity facet flips'},
        'diaphragm_axis_in_torso':diaphragm['axis_in_torso'],'material_map_audit':audit,
        'cavity':{'file':'cavity.npz','sha256':hashlib.sha256(cavity_blob).hexdigest(),'reference_geometric_volume_m3':volume,
            'vertices':len(cavity),'faces':len(faces),'source_lobes':lung_receipts,'material_entity_order':owner_ids,
            'maximum_nearest_material_node_gap_m':float(cKDTree(all_nodes).query(cavity)[0].max()),
            'derived_data_license':'CC BY-SA4.0; retained Z-Anatomy upstream attribution applies to lobe-derived cavity data',
            'upstream_license':receipt(ROOT/'data/raw/anatomy/extended/License.txt')|{'retained_copy':'Z-Anatomy-License.txt'},
            'basis':'Closed oriented convex envelope of five registered lobe surfaces, displacement-bound to nearest mechanical material nodes',
            'gas_mapping':'Optional explicit native reference supplies constant offset Vgas0−Vgeom0. No native gas store or gas-exchange area introduced.',
            'limitations':['Convex envelope fills interlobar gaps, mediastinum and concavities; it is not actual gas-exchange volume.',
                'Lobe-to-thoracic-node interpolation is engineered; nearest-node gaps retained, not calibrated pleural contact.',
                'Fixed topology/orientation checks do not prove global absence of self-intersection after deformation.']},
        'mass_ownership':'Every one of48debited material shares is represented with positive degree2 triangle quadrature; residual17.937704kg core supplies rigid6DOF inertia. No separate modal mass added.',
        'native_recoil_ownership':'BioGears left/right PleuralCavityToRespiratoryMuscle compliances remain native; mechanism has noK/damping. Its ideal pressure driver remains owner until an explicit port replaces it.',
        'pressure_work':'Qinternal=p*dV/dq. Closed-envelope uniform pressure resultant/moment should cancel. Pointloads expose6parent+26internal generalized force and fixed-parent reaction.',
        'not_implemented':['Passive constitutive laws, damping, chemical activation, native integration or coupled clock acceptance',
            'Subject-specific strain/attachment calibration and globally certified cavity/contact geometry',
            'Other whole-body segments and cervical mechanics: this is the explicitly composed reduced-torso subsystem only']}
    output.mkdir(parents=True);(output/'maps').mkdir();(output/'lung_sources').mkdir()
    for filename,blob in maps.items():(output/filename).write_bytes(blob)
    for filename,blob in lung_archives.items():(output/filename).write_bytes(blob)
    (output/'Z-Anatomy-License.txt').write_bytes((ROOT/'data/raw/anatomy/extended/License.txt').read_bytes())
    (output/'ATTRIBUTION.md').write_text('# Derived geometry attribution\n\nLobe sources and derived cavity geometry: Z-Anatomy — The libre3D atlas of anatomy, CC BY-SA4.0; BodyParts3D — Database Center for Life Science, upstream attribution retained in Z-Anatomy-License.txt. The generated lobe-derived cavity data are shared under CC BY-SA4.0.\n\nThoracic material maps derive from the identity-bound BodyParts3D source recipe (DBCLS, CC BY4.0). Source licenses concern retained/derived anatomical data; they do not assert a new license for unrelated repository code.\n')
    (output/'cavity.npz').write_bytes(cavity_blob);(output/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    return {'output':str(output),'materials':len(materials),'independent_dofs':32-len(locked),'locked':locked,'cavity_m3':volume,'cavity_vertices':len(cavity)}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);print(json.dumps(build(p.parse_args().output),indent=2))
