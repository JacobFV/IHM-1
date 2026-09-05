#!/usr/bin/env python3
"""Build a source-constrained canonical body; integration is an explicit hook."""
import argparse
from collections import Counter
import copy
import json
from pathlib import Path
import shutil
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'scripts'))
from ihm.assembly.anatomy import (FRAME, MODEL_ID, ROTATION, LandmarkRegistration,
    mesh_properties, normalized_name, physical_role, read_geometry, sha256, verify_assembly, write_json)
OUT = ROOT/'data/derived/canonical'


def file_record(path):
    path = Path(path)
    return {'path': str(path.relative_to(ROOT)), 'sha256': sha256(path)}


def uncertainty(kind, registration=None):
    return {'biological': {'status': 'generic authored reference; population variability not estimated', 'confidence_percent': None,
                           'independent_subject_count': 0, 'subject_calibrated': False},
            'registration_or_synthesis': {'kind': kind, 'registration_id': registration, 'calibrated_probability': None},
            'display_numerics': {'surface_only': True, 'display_reduction_applied': False,
                                'position_units': 'm', 'rounding': 'JSON float serialization'}}


def entity_from_geometry(identity, source_id, name, system, path, vertices, faces, evidence_kind, provenance, assumptions):
    return {'id': identity, 'model_id': MODEL_ID, 'source_id': source_id, 'name': name,
            'system': system, 'role': physical_role(name, system), 'evidence_kind': evidence_kind,
            'reference_geometry': {**file_record(path), 'frame': FRAME, 'units': 'm', 'representation': 'triangular_surface'},
            **mesh_properties(vertices, faces), 'provenance': provenance, 'assumptions': assumptions,
            'uncertainty': uncertainty(evidence_kind, 'z_anatomy' if evidence_kind == 'registered_geometry' else None),
            'connections': []}


def register_lymphatic_graph(bp_by_name, entities):
    """Articulate publisher T-pose graph into the body as a structural network.

    Graph coordinates have no anatomical landmark labels. The source pose joint
    centers are engineering priors, explicitly separate from named atlas fits.
    """
    graph_path = ROOT/'data/derived/lymphatic/graph.json'
    graph = json.loads(graph_path.read_text())
    original = np.asarray([n['position_mm'] for n in graph['nodes']])*.001
    rotated = original @ ROTATION.T
    skin = bp_by_name['skin']
    low, high = np.array(skin['bounds_m']['min']), np.array(skin['bounds_m']['max'])
    height_scale = (high[1]-low[1])/(rotated[:, 1].max()-rotated[:, 1].min())
    transformed = rotated*np.array([.65, height_scale, .8]) + np.array([0., low[1], -.045])
    articulated = np.zeros(len(rotated), dtype=bool)
    control_records = []
    for side, sign in [('left', 1), ('right', -1)]:
        humerus, radius, hand = (bp_by_name[side+' '+name] for name in ('humerus', 'radius', 'third metacarpal bone'))
        hb, rb, mb = humerus['bounds_m'], radius['bounds_m'], hand['bounds_m']
        target = np.array([[sign*.165, hb['max'][1], humerus['centroid_m'][2]],
                           [humerus['centroid_m'][0]+sign*.025, hb['min'][1], humerus['centroid_m'][2]],
                           [radius['centroid_m'][0], rb['min'][1], radius['centroid_m'][2]],
                           [hand['centroid_m'][0]+sign*.025, mb['min'][1]-.065, hand['centroid_m'][2]]])
        source = np.array([[.20, 1.43, -.035], [.48, 1.39, -.045], [.72, 1.37, -.035], [.85, 1.36, -.015]])
        source[:, 0] *= sign
        control_records.append({'side': side, 'source_control_points_m': source.tolist(), 'target_control_points_m': target.tolist(),
                                'target_support_ids': [humerus['id'], radius['id'], hand['id']]})
        selected = np.where((rotated[:, 0]*sign > .20) & (rotated[:, 1] > 1.25))[0]
        for index in selected:
            point = rotated[index]
            segment = int(np.clip(np.searchsorted(source[:, 0]*sign, point[0]*sign)-1, 0, 2))
            t = np.clip((point[0]-source[segment, 0])/(source[segment+1, 0]-source[segment, 0]), 0., 1.)
            center = (1-t)*source[segment]+t*source[segment+1]
            destination = (1-t)*target[segment]+t*target[segment+1]
            tangent = target[segment+1]-target[segment]
            tangent /= np.linalg.norm(tangent)
            outward = np.array([-tangent[1]*sign, tangent[0]*sign, 0.])
            outward /= np.linalg.norm(outward)
            delta = point-center
            transformed[index] = destination + outward*delta[1]*.65 + np.array([0., 0., delta[2]*.8])
            articulated[index] = True
    groups = [e for e in entities if e['role'] == 'lymph_node_group']
    centers = np.array([e['centroid_m'] for e in groups])
    assignments = []
    for i, node in enumerate(graph['nodes']):
        if node['is_lymph_node']:
            distances = np.linalg.norm(centers-transformed[i], axis=1)
            nearest = int(np.argmin(distances))
            assignments.append({'source_node_id': node['id'], 'candidate_entity_id': groups[nearest]['id'],
                                'distance_m': float(distances[nearest]), 'accepted': bool(distances[nearest] <= .06),
                                'interpretation': 'regional nearest-group association; not anatomical node identity'})
    identity = 'body-published-lymphatic-network'
    edges = copy.deepcopy(graph['edges'])
    for edge in edges:
        a, b = edge['source_from']-1, edge['source_to']-1
        edge['canonical_chord_length_m'] = float(np.linalg.norm(transformed[a]-transformed[b]))
    canonical_graph = {'model_id': MODEL_ID, 'frame': FRAME, 'units': 'm', 'directed': False,
                       'nodes': [{**n, 'position_m': transformed[i].tolist(), 'arm_pose_transform': bool(articulated[i])} for i,n in enumerate(graph['nodes'])],
                       'edges': edges, 'node_group_associations': assignments, 'source_graph': file_record(graph_path),
                       'interpretation': 'Registered publisher topology. Chords do not resolve vessel lumens, flow, valves or proven anatomical connections.'}
    registered_path = OUT/'lymphatic_graph.json'
    write_json(registered_path, canonical_graph)
    geometry_path = OUT/'geometry'/f'{identity}.json.gz'
    positions = [v for edge in edges for i in (edge['source_from']-1, edge['source_to']-1) for v in transformed[i].tolist()]
    write_json(geometry_path, {'positions': positions, 'units': 'm', 'frame': FRAME, 'primitive': 'line_segments', 'source_edge_ids': [e['id'] for e in edges]})
    e = {'id': identity, 'model_id': MODEL_ID, 'source_id': graph['model_id'], 'name': 'Lymphatic structural network',
         'system': 'lymphatic', 'role': 'lymphatic_network', 'evidence_kind': 'registered_structural_graph',
         'bounds_m': {'min': transformed.min(0).tolist(), 'max': transformed.max(0).tolist()},
         'centroid_m': transformed.mean(0).tolist(), 'volume_m3': None, 'volume_method': 'unknown lumen radii',
         'reference_geometry': {**file_record(geometry_path), 'frame': FRAME, 'units': 'm', 'representation': 'structural_graph_chords'},
         'graph': file_record(registered_path), 'assumptions': ['LYMPHATIC-POSE-PRIOR'],
         'provenance': {'source_ids': [graph['model_id']], 'files': [file_record(graph_path)], 'source': graph['provenance'], 'dependency_group': 'PlasticBoy-derived-published-estimate'},
         'uncertainty': uncertainty('registered_structural_graph', 'lymphatic_pose'),
         'connections': [{'entity_id': identity, 'relation': 'regional_network_association', 'evidence': 'nearest registered source lymph-node group within 60 mm; no node identity asserted'} for identity in sorted({a['candidate_entity_id'] for a in assignments if a['accepted']})]}
    report = {'method': 'body extent alignment plus articulated arm centerline pose map', 'source_graph': file_record(graph_path),
              'source_rotation': ROTATION.tolist(), 'body_scale_xyz': [.65, float(height_scale), .8], 'body_translation_m': [0., float(low[1]), -.045],
              'arm_controls': control_records, 'arm_vertex_count': int(articulated.sum()), 'source_vertex_count': len(graph['nodes']), 'source_edge_count': len(edges),
              'node_group_association_threshold_m': .06, 'accepted_node_group_associations': sum(a['accepted'] for a in assignments),
              'node_group_association_rms_m': float(np.sqrt(np.mean([a['distance_m']**2 for a in assignments]))),
              'anatomically_labeled_source_landmarks': 0, 'held_out_rms_m': None,
              'interpretation': 'Pose prior, not landmark-validated registration. Regional group distances are diagnostics, not calibrated measurement uncertainty. Source graph topology retained exactly; arm source controls are assumed from its T-pose extent.'}
    model_source = json.loads((ROOT/'data/derived/lymphatic/manifest_fragment.json').read_text())['structures'][0]['source']
    s = {'id': identity, 'model_id': MODEL_ID, 'name': e['name'], 'system': 'lymphatic', 'kind': 'lines',
         'geometry_url': '/api/geometry/'+identity, 'geometry_path': str(geometry_path.relative_to(ROOT)),
         'geometry_sha256': e['reference_geometry']['sha256'], 'color': '#78b991', 'default_visible': True,
         'source': model_source, 'calibration_status': 'registered structural topology with explicit pose priors',
         'evidence_kind': e['evidence_kind'], 'canonical_entity_id': identity, 'assumptions': e['assumptions'], 'uncertainty': e['uncertainty']}
    return e, s, report


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((ROOT/'data/derived/app/manifest.json').read_text())
    bp_model = next(x for x in manifest['models'] if x['id'] == 'bodyparts3d')
    bp_structures = [x for x in manifest['structures'] if x['model_id'] == 'bodyparts3d']
    if len(bp_structures) != 2234:
        raise ValueError('Canonical build requires the complete 2,234-mesh BP scaffold')
    bp_index = {x['element_id']: x for x in json.loads((ROOT/'data/derived/anatomy/bodyparts3d_index.json').read_text())['meshes']}
    from build_extended_anatomy import validate_provenance, validate_cached_sources
    za_provenance = validate_provenance()
    za_index = validate_cached_sources(za_provenance)
    za_fragment = json.loads((ROOT/'data/derived/anatomy/extended/manifest_fragment.json').read_text())
    za_structures = {x['id']: x for x in za_fragment['structures']}
    entities, structures = [], []
    bp_by_name = {}
    for s in bp_structures:
        identity = 'body-' + s['id']
        source_path = ROOT/'data/derived/app/geometry'/f'{s["id"]}.json.gz'
        raw = bp_index[s['id'][5:]]
        if sha256(ROOT/raw['source_path']) != raw['sha256']:
            raise ValueError('BP raw checksum mismatch: ' + s['id'])
        if sha256(source_path) != s['geometry_sha256']:
            raise ValueError('BP display checksum mismatch: ' + s['id'])
        geometry = read_geometry(source_path)
        vertices = np.asarray(geometry['positions']).reshape(-1, 3)
        faces = np.asarray(geometry['indices']).reshape(-1, 3)
        transformed_bounds = np.asarray(raw['bounds_in_source_coordinates']) @ ROTATION.T * .001 + np.asarray(bp_model['display_transform']['translation'])
        if not np.allclose(vertices.min(0), transformed_bounds.min(0), atol=1e-8, rtol=0) or not np.allclose(vertices.max(0), transformed_bounds.max(0), atol=1e-8, rtol=0):
            raise ValueError('BP canonical coordinates disagree with acquired source bounds: '+s['id'])
        if len(vertices) != raw['vertices'] or len(faces) != raw['faces']:
            raise ValueError('Canonical scaffold requires full source geometry; run scripts/build_spatial_atlas.py to rebuild full-source surfaces')
        path = OUT/'geometry'/f'{identity}.json.gz'
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_path, path)
        e = entity_from_geometry(identity, s['id'], s['name'], s['system'], path, vertices, faces,
            'source_geometry', {'source_ids': [s['id']], 'files': [file_record(ROOT/raw['source_path'])],
                'source': copy.deepcopy(s['source']), 'source_to_canonical': bp_model['display_transform'],
                'dependency_group': 'bodyparts3d-derived-reference', 'license': s['source']['license']},
            ['CANONICAL-GENERIC-REFERENCE'])
        e['concepts'] = s.get('concepts', [])
        e['classification'] = s.get('classification', {})
        entities.append(e)
        bp_by_name.setdefault(normalized_name(s['name']), e)
        display = copy.deepcopy(s)
        display.update(id=identity, model_id=MODEL_ID, geometry_url='/api/geometry/'+identity,
                       geometry_path=str(path.relative_to(ROOT)), geometry_sha256=e['reference_geometry']['sha256'],
                       evidence_kind=e['evidence_kind'], canonical_entity_id=identity,
                       uncertainty=e['uncertainty'], assumptions=e['assumptions'])
        structures.append(display)
        if len(entities) % 500 == 0:
            print(f'Canonical scaffold: {len(entities)}/2234', flush=True)

    source, target, anchors = [], [], []
    for z in za_index['meshes']:
        name = normalized_name(z['name'])
        if z['system'] == 'skeletal' and 'cartilage' not in name and name in bp_by_name:
            a = np.mean(z['bounds'], axis=0) @ ROTATION.T
            b = np.asarray(bp_by_name[name]['centroid_m'])
            source.append(a); target.append(b)
            anchors.append({'name': name, 'source_id': z['id'], 'target_id': bp_by_name[name]['id'],
                            'source_rotated_m': a.tolist(), 'target_m': b.tolist(),
                            'landmark_definition': 'axis-aligned surface bounding-box center'})
    source, target = np.array(source), np.array(target)
    registration = LandmarkRegistration(source, target)
    reg_report = registration.report()
    held_out = []
    # Deterministic 5-fold holdout, retaining a large spatially distributed support set.
    for fold in range(5):
        test = np.arange(len(source)) % 5 == fold
        fit = LandmarkRegistration(source[~test], target[~test])
        held_out.extend(np.linalg.norm(fit.transform(source[test]) - target[test], axis=1).tolist())
    reg_report.update(landmarks=anchors, held_out_method='deterministic five-fold centroid holdout',
                      held_out_rms_m=float(np.sqrt(np.mean(np.square(held_out)))), held_out_max_m=float(max(held_out)),
                      held_out_residuals_m=held_out,
                      source_rotation=ROTATION.tolist(), source_units_assumption='Blender metric scene scale_length=1; atlas registration determines final physical scale',
                      independent_subject_count=0, dependent_source='Z-Anatomy derives in part from BodyParts3D',
                      source_index=file_record(ROOT/'data/derived/anatomy/extended/source_index.json'))
    exclusions, jacobians = [], []
    missing_essential_organs = {'za-b0e9cdd77cc7f9fe', 'za-2de3d6737090b6f5', 'za-984141fb35e90177', 'za-915a592641ffe9b0', 'za-0667ce7832fc9367', 'za-05ffc920dd6517bb'}
    for z in za_index['meshes']:
        if z['id'] not in missing_essential_organs and z['system'] != 'lymphatic' and 'parathyroid' not in z['name'].lower():
            continue
        name = normalized_name(z['name'])
        if name in bp_by_name:
            exclusions.append({'source_id': z['id'], 'canonical_id': bp_by_name[name]['id'], 'reason': 'existing canonical organ; do not duplicate'})
            continue
        path_source = ROOT/z['source_geometry_path']
        if sha256(path_source) != z['source_geometry_sha256']:
            raise ValueError('Z source geometry checksum mismatch: '+z['id'])
        with np.load(path_source) as mesh:
            vertices = mesh['vertices'] @ ROTATION.T
            faces = mesh['faces']
        transformed = registration.transform(vertices)
        jacobians.extend(registration.jacobian_determinants(vertices).tolist())
        identity = 'body-' + z['id']
        path = OUT/'geometry'/f'{identity}.json.gz'
        write_json(path, {'positions': transformed.ravel().tolist(), 'indices': faces.ravel().tolist(), 'units': 'm',
                         'frame': FRAME, 'registration_id': 'z_anatomy', 'source_geometry_sha256': z['source_geometry_sha256']})
        e = entity_from_geometry(identity, z['id'], name, z['system'], path, transformed, faces, 'registered_geometry',
            {'source_ids': [z['id']], 'files': [file_record(path_source)],
             'source': copy.deepcopy(za_structures[z['id']]['source']), 'registration_id': 'z_anatomy',
             'dependency_group': 'bodyparts3d-derived-reference', 'license': 'CC-BY-SA-4.0'},
            ['CANONICAL-GENERIC-REFERENCE', 'Z-LANDMARK-REGISTRATION'])
        entities.append(e)
        s = copy.deepcopy(za_structures[z['id']])
        s.update(id=identity, name=name, model_id=MODEL_ID, geometry_url='/api/geometry/'+identity,
                 geometry_path=str(path.relative_to(ROOT)), geometry_sha256=e['reference_geometry']['sha256'],
                 default_visible=True, evidence_kind=e['evidence_kind'], canonical_entity_id=identity,
                 calibration_status='registered generic anatomy; not independently measured',
                 uncertainty=e['uncertainty'], assumptions=e['assumptions'], display_geometry={'resolution': 'full-source registered', 'source_vertices': len(vertices), 'source_faces': len(faces), 'surface_only': True})
        structures.append(s)
    reg_report['jacobian_determinant_min_at_vertices'] = float(min(jacobians))
    reg_report['jacobian_determinant_max_at_vertices'] = float(max(jacobians))
    reg_report['jacobian_vertex_samples'] = len(jacobians)
    reg_report['duplicate_sources_excluded'] = exclusions

    # Layer quadrature uses the acquired skin surface and depth intervals, not
    # offset meshes that could falsely imply validated collision-free volumes.
    skin = bp_by_name['skin']
    layers = [('epidermis', .0001, [.00005, .0002]), ('dermis', .0015, [.0005, .003]), ('hypodermis', .005, [.001, .02])]
    depth = 0.
    for name, thickness, prior_range in layers:
        e = copy.deepcopy(skin)
        e.update(id='body-skin-'+name, source_id=None, name=name, role='skin_layer', evidence_kind='synthesized_layer',
                 concepts=[], classification={}, assumptions=['SKIN-LAYER-PRIOR'],
                 provenance={'source_ids': [skin['id']], 'files': [skin['reference_geometry']], 'dependency_group': 'bodyparts3d-derived-reference'},
                 uncertainty=uncertainty('synthesized_layer'),
                 connections=[{'entity_id': skin['id'], 'relation': 'layer_of', 'evidence': 'model partition; acquired skin surface supplies quadrature support'}])
        e['reference_geometry']['representation'] = 'surface_shell_quadrature_layer'
        e['shell'] = {'thickness_m': thickness, 'prior_range_m': prior_range, 'depth_interval_m': [depth, depth+thickness],
                      'depth_direction': 'inward along surface normal; integration coordinate only',
                      'geometry_materialized': False, 'prior_source': 'explicit engineering modeling assumption; not a measured or fitted thickness'}
        e['volume_m3'] = skin['surface_area_m2']*thickness
        e['volume_method'] = 'thin-shell area times assumed thickness; not measured volume'
        entities.append(e)
        skin['connections'].append({'entity_id': e['id'], 'relation': 'has_layer', 'evidence': 'explicit shell partition'})
        depth += thickness
    skin['assumptions'].append('SKIN-LAYER-PRIOR')
    skin['mechanical_representation'] = 'boundary support for three shell layers; avoid double-counting parent mass'
    # A support relation establishes a spatial substrate without claiming a
    # histological attachment or a physiological flow connection.
    bones = [e for e in entities if e['role'] == 'rigid_bone']
    centers = np.array([e['centroid_m'] for e in bones])
    for e in entities:
        if e['role'] == 'skin_layer':
            continue
        distance = np.linalg.norm(centers - np.array(e['centroid_m']), axis=1)
        if e['role'] == 'rigid_bone':
            distance[[i for i, b in enumerate(bones) if b['id'] == e['id']]] = np.inf
        index = int(np.argmin(distance))
        e['connections'].append({'entity_id': bones[index]['id'], 'relation': 'regional_spatial_support',
            'distance_m': float(distance[index]), 'evidence': 'inferred nearest bone bounding-box center; not an attachment or contact constraint'})
    graph_entity, graph_structure, graph_registration = register_lymphatic_graph(bp_by_name, entities)
    entities.append(graph_entity)
    structures.append(graph_structure)
    assumptions = [
        {'id': 'LYMPHATIC-POSE-PRIOR', 'kind': 'registered structural topology', 'statement': 'Publisher graph lacks anatomical labels; source T-pose joint centers, transverse body scale 0.65, anterior scale 0.8 and -45 mm anterior translation are explicit pose priors. Body height follows BP skin; arm controls follow BP humerus/radius/hand. Nearest node-group associations within 60 mm are regional candidates, not anatomical identity or proven drainage. No source vertices or edges are removed.'},
        {'id': 'CANONICAL-GENERIC-REFERENCE', 'kind': 'reference selection', 'statement': 'The adult male BP atlas defines one generic body and its size. Dependent derivative sources add missing structures without adding independent subjects.'},
        {'id': 'Z-LANDMARK-REGISTRATION', 'kind': 'inferred geometry', 'statement': 'Shared named bone bounding-box centers constrain affine plus smooth residual deformation. Between anchors, smoothness is a model assumption; registration does not establish measured tissue boundaries.', 'smoothing_prior': .002, 'held_out_rms_m': reg_report['held_out_rms_m']},
        {'id': 'SKIN-LAYER-PRIOR', 'kind': 'constitutive geometry prior', 'statement': 'Epidermis/dermis/hypodermis shell quadrature uses one acquired skin surface and nonoverlapping inward depth intervals. Uniform thickness and ranges are explicit engineering assumptions; no patient measurement or regional thickness field is asserted.', 'values': [{'layer': n, 'thickness_m': t, 'prior_range_m': r} for n,t,r in layers]},
        {'id': 'REGIONAL-SUPPORT', 'kind': 'spatial relation prior', 'statement': 'Nearest-bone relations provide a regional indexing substrate. They are not insertion sites, contacts, joints, or physiological exchange pathways.'}]
    assembly = {'schema_version': 1, 'model_id': MODEL_ID, 'name': 'IHM · canonical generic human',
        'frame': {'id': FRAME, 'units': 'm', 'axes': {'x': 'left', 'y': 'superior', 'z': 'anterior'},
                  'origin': 'center of acquired BP full-body axis-aligned bounds', 'source_transform': bp_model['display_transform']},
        'entities': entities, 'registrations': {'z_anatomy': reg_report, 'lymphatic_pose': graph_registration}, 'assumption_ledger': assumptions,
        'counts': {'entities': len(entities), 'systems': dict(Counter(e['system'] for e in entities)), 'roles': dict(Counter(e['role'] for e in entities)), 'evidence_kinds': dict(Counter(e['evidence_kind'] for e in entities))},
        'limits': ['Reference surface assembly, not a validated volumetric anatomical mesh.', 'Source meshes may contain overlaps, open surfaces and interfaces; support relations do not certify contact.', 'Bone landmarks and derivative atlases are not independent subject measurements.', 'Lymph node groups preserve authored group meshes; their number is not an individual lymph-node count.', 'Anatomical support connections alone do not establish physiological causality.']}
    write_json(OUT/'anatomy.json', assembly)
    model = copy.deepcopy(bp_model)
    model.update(id=MODEL_ID, name='IHM · one generic human', description='A single body with acquired skeletal, muscle, organ, nerve and vascular surfaces; registered lung parenchyma, regional lymphatic additions and explicit skin-layer priors.',
                 anatomy_path=str((OUT/'anatomy.json').relative_to(ROOT)), source_family_ids=['bodyparts3d', 'z-anatomy', 'published-lymphatic-network'],
                 cross_family_registration=True, calibration_status='generic anatomical assembly with recorded priors; not subject calibrated',
                 assumption_ledger=assumptions, registrations={'z_anatomy': {k:v for k,v in reg_report.items() if k not in ('landmarks', 'held_out_residuals_m')}},
                 attribution=[bp_model['attribution'], za_fragment['models'][0]['attribution']], license='Composite: BP CC BY 4.0; registered Z-Anatomy additions CC BY-SA 4.0')
    write_json(OUT/'manifest_fragment.json', {'models': [model], 'structures': structures})
    report = verify_assembly(assembly, ROOT)
    write_json(OUT/'verification.json', report)
    print(json.dumps(report, indent=2), flush=True)
    return assembly


def append_to_manifest(manifest_path=ROOT/'data/derived/app/manifest.json'):
    fragment = json.loads((OUT/'manifest_fragment.json').read_text())
    verify_assembly(json.loads((OUT/'anatomy.json').read_text()), ROOT)
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    for key in ('models', 'structures'):
        manifest[key] = [x for x in manifest[key] if x.get('model_id', x['id']) != MODEL_ID] + fragment[key]
    for s in fragment['structures']:
        source = ROOT/s['geometry_path']
        if sha256(source) != s['geometry_sha256']:
            raise ValueError('Canonical geometry changed: '+s['id'])
        destination = manifest_path.parent/'geometry'/f'{s["id"]}.json.gz'
        if source.resolve() != destination.resolve():
            shutil.copyfile(source, destination)
    temporary = manifest_path.with_suffix('.canonical.tmp')
    write_json(temporary, manifest)
    temporary.replace(manifest_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--append', action='store_true')
    args = parser.parse_args()
    build()
    if args.append:
        append_to_manifest()
