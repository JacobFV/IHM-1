#!/usr/bin/env python3
"""Pinned bone-frame quadrant masks: engineered material regions, not lymph watersheds."""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.lymph_registration import SkinTerritoryInventory

BONES = {'left': ('body-bp3d-FJ3282', 'body-bp3d-FJ3260'),
         'right': ('body-bp3d-FJ3387', 'body-bp3d-FJ3366')}
DEFAULTS = {'end_quantile': .05, 'boundary_band_m': .003, 'maximum_radius_m': .13}


def surface_diagnostic(xyz, triangles):
    """Detect duplicate facets and index-connected shells without deleting data."""
    directed_edges=np.concatenate([triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]]])
    edges=np.sort(directed_edges,axis=1)
    edge,inverse_edges,counts=np.unique(edges,axis=0,return_inverse=True,return_counts=True)
    orientation_sum=np.bincount(inverse_edges,weights=np.where(directed_edges[:,0]<directed_edges[:,1],1,-1))
    graph=coo_matrix((np.ones(2*len(edge)),
                     (np.concatenate([edge[:,0],edge[:,1]]),np.concatenate([edge[:,1],edge[:,0]]))),
                    shape=(len(xyz),len(xyz))).tocsr()
    count,components=connected_components(graph,directed=False)
    face_components=components[triangles[:,0]]
    area=np.linalg.norm(np.cross(xyz[triangles[:,1]]-xyz[triangles[:,0]],
                                 xyz[triangles[:,2]]-xyz[triangles[:,0]]),axis=1)/2
    duplicates={}
    for decimals in [8,6]:
        unique,inverse=np.unique(xyz.round(decimals),axis=0,return_inverse=True)
        facets=np.sort(inverse[triangles],axis=1)
        duplicates[str(decimals)]={'coincident_vertices':len(xyz)-len(unique),
                                  'coincident_facets_either_winding':len(facets)-len(np.unique(facets,axis=0))}
    rows=[]
    for component in range(count):
        f=triangles[face_components==component]; points=xyz[components==component]
        if not len(f):continue
        rows.append({'id':component,'triangle_count':len(f),'area_m2':float(area[face_components==component].sum()),
                     'boundary_edges':int(np.sum((counts==1)&(components[edge[:,0]]==component))),
                     'contact_eligibility':'candidate_until_largest_component_selection',
                     'bounds_m':[points.min(0).tolist(),points.max(0).tolist()],
                     'signed_volume_integral_m3':float(np.einsum('ij,ij->i',xyz[f[:,0]],np.cross(xyz[f[:,1]],xyz[f[:,2]])).sum()/6)})
    rows.sort(key=lambda r:(-r['area_m2'],r['id']))
    outer=rows[0]
    if outer['signed_volume_integral_m3']<=0 or any(r['coincident_facets_either_winding'] for r in duplicates.values()):
        raise ValueError('Explicit outer-surface review required')
    for row in rows:
        row['contact_eligibility']='engineered_exterior_proxy' if row['id']==outer['id'] else 'excluded_inner_or_seam_component'
    eligible=face_components==outer['id']
    return {'components':rows,'duplicate_coordinate_diagnostics':duplicates,
            'face_component_ids':face_components.tolist(),
            'inconsistent_interior_edge_winding':int(np.sum((counts==2)&(orientation_sum!=0))),
            'boundary_edges':int(np.sum(counts==1)),'nonmanifold_edges':int(np.sum(counts>2)),
            'selected_component_id':outer['id'],'selected_area_m2':outer['area_m2'],
            'selection_evidence':'inferred_largest_component_positive_signed_volume_integral',
            'signed_integral_is_closed_volume':False,'self_intersections_tested':False,
            'physical_surface_exclusivity_validated':False},eligible


def generate(root=ROOT, parameters=None):
    parameters = dict(DEFAULTS if parameters is None else parameters)
    if set(parameters) != set(DEFAULTS):
        raise ValueError('Unknown or missing engineering parameter')
    for value in parameters.values():
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not np.isfinite(value):
            raise ValueError('Invalid engineering parameter')
    if not 0 < parameters['end_quantile'] < .25 or not 0 <= parameters['boundary_band_m'] < .03 or not .03 < parameters['maximum_radius_m'] < .2:
        raise ValueError('Engineering parameter outside bounded recipe')
    anatomy_path = root / 'data/derived/canonical/anatomy.json'
    anatomy_raw = anatomy_path.read_bytes()
    anatomy = json.loads(anatomy_raw)
    if anatomy['frame']['id'] != 'bodyparts3d-display-m':
        raise ValueError('Unrecognized canonical coordinate frame')
    entities = {e['id']: e for e in anatomy['entities']}
    receipts = []

    def mesh(entity_id):
        entity = entities[entity_id]
        ref = entity['reference_geometry']; raw = (root / ref['path']).read_bytes()
        if hashlib.sha256(raw).hexdigest() != ref['sha256'] or ref['units'] != 'm' or ref['frame'] != anatomy['frame']['id']:
            raise ValueError('Changed or incompatible geometry')
        receipts.append({'entity_id': entity_id, 'name': entity['name'], 'path': ref['path'],
                         'sha256': ref['sha256'], 'bytes': len(raw), 'provenance': entity['provenance']})
        d = json.loads(gzip.decompress(raw) if ref['path'].endswith('.gz') else raw)
        return np.asarray(d['positions'], dtype=float).reshape(-1,3), np.asarray(d['indices'], dtype=int).reshape(-1,3)

    skin_id = 'body-bp3d-FJ2810'; xyz, triangles = mesh(skin_id)
    diagnostic, exterior = surface_diagnostic(xyz, triangles)
    vertices = xyz[triangles]; labels = np.zeros(len(triangles), dtype=np.uint8)
    regions, frames = [], {}
    for side, (tibia_id, fibula_id) in BONES.items():
        tibia, _ = mesh(tibia_id); fibula, _ = mesh(fibula_id)
        center = tibia.mean(axis=0)
        _, _, vh = np.linalg.svd(tibia-center, full_matrices=False)
        superior = vh[0] * (1 if vh[0,1] > 0 else -1)
        lateral = fibula.mean(axis=0) - center
        lateral -= np.dot(lateral, superior)*superior; lateral /= np.linalg.norm(lateral)
        anterior = np.cross(lateral, superior)
        if anterior[2] < 0: anterior *= -1
        anterior /= np.linalg.norm(anterior)
        projected = (tibia-center) @ superior
        lo, hi = np.quantile(projected, [parameters['end_quantile'],1-parameters['end_quantile']])
        rel = vertices-center; along = rel @ superior; lat = rel @ lateral; ant = rel @ anterior
        radius = np.sqrt(lat*lat+ant*ant)
        side_gate = np.all(vertices[:,:,0] > 0, axis=1) if side == 'left' else np.all(vertices[:,:,0] < 0, axis=1)
        eligible = (np.all((along>lo)&(along<hi)&(radius<parameters['maximum_radius_m']),axis=1) & side_gate & exterior)
        frames[side] = {'origin_m':center.tolist(), 'superior':superior.tolist(), 'lateral':lateral.tolist(),
                        'anterior':anterior.tolist(), 'axial_interval_m':[float(lo),float(hi)],
                        'definition':'Tibia vertex-PCA axis; fibula centroid defines lateral; +canonical z sets anterior sign'}
        band = parameters['boundary_band_m']
        for name, lat_sign, ant_sign in [('anteromedial',-1,1),('anterolateral',1,1),('posteromedial',-1,-1),('posterolateral',1,-1)]:
            mask = eligible & np.all(lat*lat_sign>band,axis=1) & np.all(ant*ant_sign>band,axis=1)
            if np.any(labels[mask]) or not np.any(mask):
                raise ValueError('Overlap or empty engineered quadrant')
            rid = side+'_lower_leg_'+name; label = len(regions)+1; labels[mask]=label
            regions.append({'id':rid,'label':label,'triangle_ids':np.flatnonzero(mask).tolist(),
                            'geometry_sha256':receipts[0]['sha256'],
                            'annotation_receipt':{'kind':'inferred_bone_frame_geometric_partition',
                                                  'recipe':'materialize_engineered_skin_territories_v1',
                                                  'anatomical_drainage_mask':False},'native_owner':None})
    inventory = SkinTerritoryInventory(root/receipts[0]['path'],receipts[0]['sha256'],regions).report()
    return {'schema':'engineered_skin_territories_v1','parameters':parameters,'frames':frames,
            'anatomy_sha256':hashlib.sha256(anatomy_raw).hexdigest(),'source_receipts':receipts,
            'regions':regions,'inventory':inventory,'unresolved_label':0,
            'surface_diagnostic':diagnostic,
            'contact_eligible_triangle_ids':np.flatnonzero(exterior).tolist(),
            'contact_excluded_triangle_count':int(np.sum(~exterior)),
            'contact_unresolved_area_m2':diagnostic['selected_area_m2']-sum(r['area_m2'] for r in inventory['regions'] if r['id']!='unresolved'),
            'evidence_kind':'inferred_geometry_partitions_not_measured_lymph_or_venous_watersheds',
            'anatomical_basis':'Lower-leg four-pathway terminology: PMID31746690; geometric boundaries are engineering choices.',
            'native_configuration':None,'native_commands':[]}


def three_region_proposal(materialization, first_ids, second_ids):
    """Explicit caller-selected unions and complement; never emit a/b commands."""
    if not first_ids or not second_ids or len(set(first_ids))!=len(first_ids) or len(set(second_ids))!=len(second_ids) or set(first_ids)&set(second_ids):
        raise ValueError('Two nonempty exclusive territory unions required')
    rows = {r['id']:r for r in materialization['inventory']['regions']}
    if any(i not in rows or i=='unresolved' for i in first_ids+second_ids):
        raise ValueError('Unknown or unresolved selected territory')
    total = materialization['surface_diagnostic']['selected_area_m2']
    a = sum(rows[i]['area_m2'] for i in first_ids)/total
    b = sum(rows[i]['area_m2'] for i in second_ids)/total
    if a+b>=1:raise ValueError('Positive complement required')
    return {'selected_territory_unions':[list(first_ids),list(second_ids)],
            'initial_area_weights':[a,b,1-a-b], 'area_denominator_m2':total,
            'surface_evidence':'inferred_outer_component_proxy_not_clinical_body_surface_area', 'weight_evidence':'uniform_thickness_and_extracellular_fraction_engineering_prior',
            'installed_native_binding':None,'requires_new_native_variant':True,'native_commands':[]}


def reduce_contacts(materialization, contacts, root=ROOT):
    """Reject inner/seam contacts, preserve full-source ownership and outer area."""
    eligible=set(materialization['contact_eligible_triangle_ids'])
    if any(c['triangle_id'] not in eligible for c in contacts):
        raise ValueError('Contact lies on excluded inner/seam source geometry')
    source=materialization['source_receipts'][0]
    model=SkinTerritoryInventory(root/source['path'],source['sha256'],materialization['regions'])
    result=model.reduce_contacts(contacts)
    unresolved=result['regions']['unresolved']
    unresolved['whole_territory_mean_pressure_pa']=unresolved['normal_force_n']/materialization['contact_unresolved_area_m2']
    result['contact_surface_evidence']='inferred_outer_component_proxy'
    return result


if __name__ == '__main__':
    output = ROOT/'data/research/engineered_skin_territories'
    output.mkdir(parents=True,exist_ok=True)
    result = generate()
    (output/'materialization.json').write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps({'source_triangles':result['inventory']['triangle_count'],
                      'contact_excluded_triangles':result['contact_excluded_triangle_count'],
                      'outer_proxy_area_m2':result['surface_diagnostic']['selected_area_m2'],
                      'assigned_triangles':sum(len(r['triangle_ids']) for r in result['regions'])}))
