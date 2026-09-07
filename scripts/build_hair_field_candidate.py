#!/usr/bin/env python3
"""Promotion CANDIDATE: hair as named anatomical fields on the measured outer envelope.

Emits one canonical entity record per hair field, its literature record, a display
manifest fragment and an input-hashing manifest. Writes only under
data/derived/hair-field-candidate-v1/. No canonical asset is modified.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import platform
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.anatomy import FRAME, MODEL_ID, ROTATION, LandmarkRegistration  # noqa: E402
from ihm.assembly.hair_fields import (APPENDAGE_SHELLS, FIELDS, LEXICON, ROLE_VOCABULARY, SHELL_FIELDS,
                                      SHELL_STANDOFF_M, area_sample, assign_regions, face_areas, lexicon_index,
                                      region_base, shell_footprint, strand_mass_kg)  # noqa: E402
from hair_field_sources import ASSUMPTION_LEDGER, FIELD_EVIDENCE, MORPHOLOGY, NOTES, SOURCES  # noqa: E402

OUT = ROOT/'data/derived/hair-field-candidate-v1'
SCHEMA = 'ihm.hair-field-promotion-candidate.v1'
SAMPLE_SPACING_M = .003
SHELL_SPACING_M = .002
SEED = 20260907
INPUTS = ['data/derived/canonical/anatomy.json', 'data/derived/outer-envelope/outer-envelope.npz',
          'data/derived/outer-envelope/manifest.json', 'data/derived/anatomy/extended/manifest_fragment.json',
          'data/derived/anatomy/extended/source_index.json', 'data/derived/hair/elastic_v3/manifest_fragment.json',
          'data/derived/interstitial-composition-prior-v1/ledger.json', 'ihm/assembly/hair_fields.py',
          'ihm/assembly/anatomy.py', 'scripts/hair_field_sources.py', 'scripts/build_hair_field_candidate.py']


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json.dumps(value, separators=(',', ':'), allow_nan=False, sort_keys=False).encode()
    if path.suffix == '.gz':
        temporary = path.with_name(path.name + '.writing')
        with temporary.open('wb') as stream, gzip.GzipFile(filename='', fileobj=stream, mode='wb', mtime=0) as zipped:
            zipped.write(raw)
        temporary.replace(path)
    else:
        path.write_bytes(raw)
    return sha256(path)


def registration(anatomy):
    marks = anatomy['registrations']['z_anatomy']['landmarks']
    fit = LandmarkRegistration(np.array([m['source_rotated_m'] for m in marks]),
                               np.array([m['target_m'] for m in marks]))
    if not np.allclose(fit.report()['affine_4x3'], anatomy['registrations']['z_anatomy']['affine_4x3']):
        raise ValueError('z_anatomy registration did not reproduce from the stored landmarks')
    return fit


def extended_surfaces(fit):
    index = json.loads((ROOT/'data/derived/anatomy/extended/source_index.json').read_bytes())
    fragment = json.loads((ROOT/'data/derived/anatomy/extended/manifest_fragment.json').read_bytes())
    by_id = {m['id']: m for m in index['meshes']}
    surfaces, hashes = {}, {}
    for structure in fragment['structures']:
        if structure.get('system') != 'integumentary':
            continue
        entry = by_id[structure['id']]
        path = ROOT/entry['source_geometry_path']
        digest = sha256(path)
        if digest != entry['source_geometry_sha256']:
            raise ValueError('Extended source geometry changed: ' + structure['id'])
        with np.load(path) as mesh:
            vertices, faces = mesh['vertices'] @ ROTATION.T, mesh['faces']
        surfaces[structure['name']] = (structure['id'], fit.transform(vertices), faces)
        hashes[entry['source_geometry_path']] = digest
    return surfaces, hashes


def partition(surfaces, envelope_centroids, rng):
    labels, points, names = [], [], []
    for name, (_, vertices, faces) in sorted(surfaces.items()):
        if name in APPENDAGE_SHELLS:
            continue
        sample, _ = area_sample(vertices, faces, SAMPLE_SPACING_M, rng)
        names.append(name)
        points.append(sample)
        labels.append(np.full(len(sample), len(names) - 1))
    region, distance = assign_regions(envelope_centroids, np.concatenate(points), np.concatenate(labels))
    return np.array(names), region, distance


def build(strict=True):
    anatomy = json.loads((ROOT/'data/derived/canonical/anatomy.json').read_bytes())
    fit = registration(anatomy)
    rng = np.random.default_rng(SEED)
    surfaces, source_hashes = extended_surfaces(fit)
    with np.load(ROOT/'data/derived/outer-envelope/outer-envelope.npz') as mesh:
        vertices, faces = mesh['positions'], mesh['indices']
    areas = face_areas(vertices, faces)
    centroids = vertices[faces].mean(axis=1)
    names, region, region_distance = partition(surfaces, centroids, rng)
    index = lexicon_index()
    missing = sorted({region_base(n) for n in names} - set(index))
    if missing:
        raise ValueError('Surface regions absent from the hair-field lexicon: ' + ', '.join(missing))
    field = np.array([index[region_base(n)] for n in names])[region]

    shell_report = {}
    for identity in ('scalp', 'eyebrow', 'eyelash', 'pubic'):
        sample = np.concatenate([area_sample(surfaces[n][1], surfaces[n][2], SHELL_SPACING_M, rng)[0]
                                 for n in SHELL_FIELDS[identity]])
        inside, distance = shell_footprint(centroids, sample, SHELL_STANDOFF_M)
        sweep = {}
        for threshold in (.004, .006, .008, .012, .020):
            sweep[f'{threshold:.3f}'] = {'faces': int((distance < threshold).sum()),
                                         'area_m2': float(areas[distance < threshold].sum())}
        displaced = sorted({f: float(areas[inside & (field == f)].sum()) for f in set(field[inside])}.items(),
                           key=lambda kv: -kv[1])
        field[inside] = identity
        shell_report[identity] = {'shells': list(SHELL_FIELDS[identity]), 'standoff_m': SHELL_STANDOFF_M,
                                  'faces': int(inside.sum()), 'area_m2': float(areas[inside].sum()),
                                  'standoff_sweep_area_m2': sweep,
                                  'lexicon_area_absorbed_m2': [{'field': f, 'area_m2': a} for f, a in displaced]}

    records, entities, structures, geometry_hashes = [], [], [], {}
    total_area = float(areas.sum())
    for identity, (label, kind, laterality) in FIELDS.items():
        mask = field == identity
        evidence = FIELD_EVIDENCE[identity]
        area = float(areas[mask].sum())
        resolved = bool(mask.any())
        entity_id = 'body-hair-' + identity.replace('_', '-')
        geometry = None
        if resolved:
            patch_vertices, patch_faces = _submesh(vertices, faces, mask)
            path = OUT/'geometry'/(entity_id + '.json.gz')
            digest = write_json(path, {'positions': patch_vertices.ravel().tolist(),
                                       'indices': patch_faces.ravel().tolist(), 'units': 'm', 'frame': FRAME,
                                       'registration_id': 'z_anatomy',
                                       'source': 'data/derived/outer-envelope/outer-envelope.npz'})
            geometry_hashes[str(path.relative_to(ROOT))] = digest
            geometry = {'path': str(path.relative_to(ROOT)), 'sha256': digest, 'frame': FRAME, 'units': 'm',
                        'representation': 'triangular_surface',
                        'basis': 'outer body envelope faces carrying this field; a patch of the measured envelope, '
                                 'not an offset or authored shell'}
        record = _field_record(identity, label, kind, laterality, evidence, area, resolved,
                               names, region, field, mask, areas, total_area)
        records.append(record)
        entities.append(_entity(entity_id, label, identity, kind, laterality, record, geometry,
                                patch_vertices if resolved else None, patch_faces if resolved else None,
                                source_hashes))
        structures.append(_structure(entity_id, label, identity, record, geometry))

    summary = _summary(records, anatomy, total_area, shell_report, region_distance, names, field, areas)
    if strict:
        _check(records, summary, total_area)
    return {'anatomy': anatomy, 'records': records, 'entities': entities, 'structures': structures,
            'summary': summary, 'shell_report': shell_report, 'source_hashes': source_hashes,
            'geometry_hashes': geometry_hashes,
            'region_distance_m': {'median': float(np.median(region_distance)),
                                  'p95': float(np.percentile(region_distance, 95)),
                                  'max': float(region_distance.max())}}


def _submesh(vertices, faces, mask):
    from ihm.assembly.hair_fields import submesh
    return submesh(vertices, faces, mask)


def _field_record(identity, label, kind, laterality, evidence, area, resolved, names, region, field, mask,
                  areas, total_area):
    density = copy.deepcopy(evidence['follicle_density_cm2'])
    counted = evidence.get('follicle_count')
    composition = sorted({names[r]: float(areas[(region == r) & mask].sum())
                          for r in np.unique(region[mask])}.items(), key=lambda kv: -kv[1]) if resolved else []
    if density.get('subregion_cm2') and resolved:
        table = density['subregion_cm2']
        weighted = [(table[region_base(n)], a) for n, a in composition if region_base(n) in table]
        covered = sum(a for _, a in weighted)
        density['value'] = float(sum(v*a for v, a in weighted)/covered) if covered else None
        density['weighting'] = {'basis': 'footprint-area weighted over the surface regions the source measured',
                                'covered_area_fraction': covered/area if area else None,
                                'weights_cm2': [{'region': n, 'area_cm2': a*1e4, 'density_cm2': table[region_base(n)]}
                                                for n, a in composition if region_base(n) in table]}
    follicles = None
    if counted is not None:
        follicles = None if counted['value'] is None else float(counted['value'])
    elif density['value'] is not None and resolved:
        follicles = float(density['value'])*area*1e4
    morphology = MORPHOLOGY[identity]
    rho = SOURCES[MORPHOLOGY['_fibre_density_source']]['fibre_density_kg_m3']
    radius, length, fraction = morphology['shaft_radius_m'], morphology['shaft_length_m'], \
        morphology['shaft_bearing_fraction']
    mass = band = None
    if follicles is not None and None not in (radius['value'], length['value'], fraction['value']):
        mass = strand_mass_kg(follicles*fraction['value'], radius['value'], length['value'], rho)
        if radius.get('band') and length.get('band'):
            band = [strand_mass_kg(follicles*fraction['value'], radius['band'][i], length['band'][i], rho)
                    for i in (0, 1)]
    return {
        'id': identity, 'name': label, 'hair_kind': kind, 'laterality': laterality,
        'footprint': {'resolved': resolved,
                      'area_m2': area if resolved else None,
                      'area_cm2': area*1e4 if resolved else None,
                      'envelope_face_count': int(mask.sum()),
                      'envelope_area_fraction': area/total_area if resolved else None,
                      'basis': 'outer body envelope faces, area weighted',
                      'surface_regions_cm2': [{'region': n, 'area_cm2': a*1e4} for n, a in composition],
                      'unresolved_reason': None if resolved else evidence['unresolved_reason']},
        'follicle_density_cm2': density,
        'follicle_count': {'value': follicles,
                           'range': None if counted is None else counted.get('range'),
                           'tier': 'derived' if counted is None else counted['tier'],
                           'basis': (counted['basis'] if counted is not None else
                                     'follicle density times measured envelope footprint area'),
                           'source': counted['source'] if counted is not None else density['source']},
        'morphology': morphology,
        'hair_mass_kg': {'value': mass, 'band': band, 'tier': 'derived',
                         'fibre_density_kg_m3': rho, 'fibre_density_source': MORPHOLOGY['_fibre_density_source'],
                         'basis': 'follicle count times shaft-bearing fraction times pi r^2 L rho; the band is '
                                  'the same product at the low and high ends of the radius and length bands',
                         'cross_section_model': MORPHOLOGY['_cross_section_model'],
                         'absent_reason': None if mass is not None else
                         ('no follicle count' if follicles is None else
                          '; '.join(q['basis'] for q, k in ((radius, 'radius'), (length, 'length'),
                                                            (fraction, 'fraction')) if q['value'] is None))},
        'notes': evidence.get('notes', []),
    }


def _entity(entity_id, label, identity, kind, laterality, record, geometry, vertices, faces, source_hashes):
    bounds = centroid = principal = None
    if vertices is not None:
        low, high = vertices.min(axis=0), vertices.max(axis=0)
        bounds = {'min': low.tolist(), 'max': high.tolist()}
        centroid = ((low + high)/2).tolist()
        centred = vertices - vertices.mean(axis=0)
        principal = np.linalg.svd(centred, full_matrices=False)[2][0].tolist()
    return {
        'id': entity_id, 'model_id': MODEL_ID, 'source_id': None, 'name': label, 'system': 'integumentary',
        'role': ROLE_VOCABULARY['role'], 'role_vocabulary': ROLE_VOCABULARY,
        'evidence_kind': 'derived_surface_field',
        'display_group': 'hair',
        'hair_field': {'field_id': identity, 'hair_kind': kind, 'laterality': laterality,
                       'follicle_density_cm2': record['follicle_density_cm2']['value'],
                       'follicle_count': record['follicle_count']['value'],
                       'hair_mass_kg': record['hair_mass_kg']['value']},
        'reference_geometry': geometry,
        'bounds_m': bounds, 'centroid_m': centroid,
        'centroid_definition': 'axis-aligned surface bounding-box center; not center of mass',
        'surface_area_m2': record['footprint']['area_m2'],
        'volume_m3': None,
        'volume_method': 'not applicable: a hair field is a surface population, not a closed volume',
        'watertight_edge_incidence': False,
        'principal_axis': principal,
        'source_vertex_count': None if vertices is None else int(len(vertices)),
        'source_face_count': None if faces is None else int(len(faces)),
        'provenance': {
            'source_ids': [],
            'dataset': 'BodyParts3D 4.0 skin surface, Z-Anatomy Terminologia Anatomica surface regions, '
                       'and the hair literature in sources.json',
            'script': 'scripts/build_hair_field_candidate.py',
            'transform': 'Z-Anatomy source vertices rotated by the canonical source rotation, registered by the '
                         'stored z_anatomy affine plus thin-plate-spline landmark fit, area sampled at '
                         f'{SAMPLE_SPACING_M} m, then each outer-envelope face assigned to its nearest surface '
                         'region sample; authored hair shells override the region lexicon within '
                         f'{SHELL_STANDOFF_M} m',
            'registration_id': 'z_anatomy',
            'files': [{'path': p, 'sha256': sha256(ROOT/p)} for p in
                      ('data/derived/outer-envelope/outer-envelope.npz',
                       'data/derived/anatomy/extended/source_index.json',
                       'data/derived/canonical/anatomy.json')],
            'surface_region_source_geometry_count': len(source_hashes),
            'dependency_group': 'bodyparts3d-derived-reference',
            'license': 'CC BY 4.0 (BodyParts3D) and CC-BY-SA-4.0 (Z-Anatomy)',
        },
        'assumptions': ['CANONICAL-GENERIC-REFERENCE', 'Z-LANDMARK-REGISTRATION', 'HAIR-FIELD-LEXICON',
                        'HAIR-DENSITY-TRANSFER'],
        'uncertainty': {
            'biological': {'status': 'one generic authored body; follicle density and morphology are transferred '
                                     'from published cohorts that are not this specimen',
                           'confidence_percent': None, 'independent_subject_count': 0, 'subject_calibrated': False},
            'registration_or_synthesis': {'kind': 'derived_surface_field', 'registration_id': 'z_anatomy',
                                          'calibrated_probability': None},
            'display_numerics': {'surface_only': True, 'display_reduction_applied': False, 'position_units': 'm',
                                 'rounding': 'JSON float serialization'},
            'evidence_tiers': {'follicle_density': record['follicle_density_cm2']['tier'],
                               'shaft_radius': record['morphology']['shaft_radius_m']['tier'],
                               'shaft_length': record['morphology']['shaft_length_m']['tier'],
                               'footprint': 'measured' if record['footprint']['resolved'] else 'unresolved'},
        },
        'connections': [{'entity_id': 'body-bp3d-FJ2810', 'relation': 'field_on',
                         'evidence': 'the outer sheet of the acquired skin slab carries this field'}],
        'concepts': [], 'classification': {'terminologia_anatomica_regions':
                                           [r['region'] for r in record['footprint']['surface_regions_cm2']]},
    }


def _structure(entity_id, label, identity, record, geometry):
    return {'id': entity_id, 'name': label, 'model_id': MODEL_ID, 'system': 'hair', 'kind': 'mesh',
            'color': '#382b20', 'default_visible': False, 'canonical_entity_id': entity_id,
            'canonical_system': 'integumentary',
            'evidence_kind': 'derived_surface_field',
            'calibration_status': 'registered surface region footprint with transferred follicle density',
            'geometry_url': '/api/geometry/' + entity_id,
            'geometry_path': None if geometry is None else geometry['path'],
            'geometry_sha256': None if geometry is None else geometry['sha256'],
            'hair_field': {'field_id': identity, 'follicle_count': record['follicle_count']['value'],
                           'area_cm2': record['footprint']['area_cm2']}}


def _summary(records, anatomy, total_area, shell_report, region_distance, names, field, areas):
    resolved = [r for r in records if r['footprint']['resolved']]
    follicles = [r['follicle_count']['value'] for r in records if r['follicle_count']['value'] is not None]
    masses = [r['hair_mass_kg']['value'] for r in records if r['hair_mass_kg']['value'] is not None]
    ledger = json.loads((ROOT/'data/derived/interstitial-composition-prior-v1/ledger.json').read_bytes())
    body = ledger['cases']['voxel_8mm']['implied_total_body_mass_kg']
    verified = sum(1 for r in records if r['follicle_density_cm2']['tier'] == 'measured')
    return {
        'entities_before': {'total': anatomy['counts']['entities'],
                            'integumentary': anatomy['counts']['systems']['integumentary'],
                            'named_hair_entities': 3,
                            'named_hair_entity_ids': ['body-bp3d-FJ2813', 'body-bp3d-FJ2815', 'body-bp3d-FJ2812'],
                            'dynamic_display_hair_structures': 2},
        'entities_after': {'total': anatomy['counts']['entities'] + len(records),
                           'integumentary': anatomy['counts']['systems']['integumentary'] + len(records),
                           'hair_field_entities': len(records)},
        'footprint': {'envelope_area_m2': total_area,
                      'partitioned_area_m2': float(sum(r['footprint']['area_m2'] for r in resolved)),
                      'fields_resolved': len(resolved), 'fields_unresolved': len(records) - len(resolved),
                      'unresolved_fields': [r['id'] for r in records if not r['footprint']['resolved']],
                      'surface_regions_used': int(len(names)),
                      'region_assignment_distance_m': {'median': float(np.median(region_distance)),
                                                       'p95': float(np.percentile(region_distance, 95)),
                                                       'max': float(region_distance.max())}},
        'population': {'total_follicles': float(sum(follicles)),
                       'total_hair_mass_kg': float(sum(masses)),
                       'total_hair_mass_band_kg': [float(sum(r['hair_mass_kg']['band'][i] for r in records
                                                             if r['hair_mass_kg'].get('band'))) for i in (0, 1)],
                       'fields_with_follicle_count': len(follicles),
                       'fields_without_follicle_count': len(records) - len(follicles),
                       'fields_without_follicle_count_ids': [r['id'] for r in records
                                                             if r['follicle_count']['value'] is None],
                       'fields_with_mass': len(masses),
                       'fields_without_mass_ids': [r['id'] for r in records
                                                   if r['hair_mass_kg']['value'] is None],
                       'scalp_share_of_mass': (next(r['hair_mass_kg']['value'] for r in records
                                                    if r['id'] == 'scalp')/float(sum(masses))) if masses else None,
                       'area_weighted_mean_density_cm2': float(sum(follicles))/(total_area*1e4)},
        'cross_checks': {
            'pooled_density': {'model_area_weighted_cm2': float(sum(follicles))/(total_area*1e4),
                               'xu2017_pooled_ex_vivo_cm2': 140, 'xu2017_ci_cm2': [91, 189],
                               'reading': 'the model sits about four times below the pooled ex vivo '
                                          'meta-analysis. The gap is method: cyanoacrylate stripping counts '
                                          'orifices emerging at the surface, histology counts follicles in '
                                          'section. Nothing here reconciles them and the model does not adopt '
                                          'the higher number, because no cited cohort measured it site by site.'},
            'scalp_area': {'model_cm2': next(r['footprint']['area_cm2'] for r in records if r['id'] == 'scalp'),
                           'patel2026_scalp_surface_area_cm2': 492.36,
                           'patel2026_range_cm2': [447.29, 592.52],
                           'reading': 'their measured region is the head crown including frontal and mid-scalp; '
                                      'this footprint is the whole hair-bearing scalp including occipital and '
                                      'temporal, so it should be larger. It is, by 1.6 times.'},
            'scalp_count': {'model': next(r['follicle_count']['value'] for r in records if r['id'] == 'scalp'),
                            'patel2026_crown_hair_count': [63996, 66426],
                            'reading': 'same region caveat. No primary whole-scalp follicle count was retrieved '
                                       'in this sweep, so the familiar 100000 figure is not used as a target.'},
            'scalp_mass_independent_route': {
                'method': 'wikramanayake2012 hair mass index, mm2 of fibre cross-section per cm2 of scalp, times '
                          'footprint area times assumed length times fibre density. It shares no input with the '
                          'count route except the footprint and the fibre density.',
                'hmi_mm2_per_cm2': [.75, 1.00], 'hmi_average_mm2_per_cm2': .87,
                'kg': [float(h*1e-6*1e4*next(r['footprint']['area_m2'] for r in records if r['id'] == 'scalp')
                             * .03*1312.) for h in (.75, .87, 1.00)],
                'count_route_kg': next(r['hair_mass_kg']['value'] for r in records if r['id'] == 'scalp'),
                'reading': 'both routes use the 30 mm assumed length, so they agree on nothing about length; '
                           'what they check independently is fibre cross-section per unit scalp area.'}},
        'mass_ledger': {'body_mass_kg': body,
                        'source': 'data/derived/interstitial-composition-prior-v1/ledger.json voxel_8mm '
                                  'central_assumed_1030',
                        'hair_mass_kg': float(sum(masses)),
                        'hair_fraction_of_body_mass': float(sum(masses))/body,
                        'hair_mass_ppm_of_body_mass': float(sum(masses))/body*1e6,
                        'accounting': 'hair sits outside the envelope and outside the 8 mm interior occupancy that '
                                      'the ledger partitions, so it is an addition to the 70.7713 kg total, not a '
                                      'reallocation within it'},
        'evidence_split': {'fields': len(records), 'density_measured': verified,
                           'density_transferred': sum(1 for r in records
                                                      if r['follicle_density_cm2']['tier'] == 'transferred'),
                           'density_assumed': sum(1 for r in records
                                                  if r['follicle_density_cm2']['tier'] == 'assumed'),
                           'density_absent': sum(1 for r in records
                                                 if r['follicle_density_cm2']['value'] is None)},
        'shell_footprints': shell_report,
        'reconciliation': NOTES['reconciliation'],
        'limits': NOTES['limits'],
    }


def _check(records, summary, total_area):
    partitioned = summary['footprint']['partitioned_area_m2']
    if abs(partitioned - total_area) > 1e-9:
        raise ValueError(f'Field footprints do not partition the envelope: {partitioned} against {total_area}')
    for record in records:
        density = record['follicle_density_cm2']
        if density['value'] is not None and density['tier'] not in ('measured', 'transferred', 'assumed'):
            raise ValueError('Bad evidence tier on ' + record['id'])
        if density['value'] is not None and not density['source']:
            raise ValueError('Density without a source on ' + record['id'])
    if summary['entities_after']['total'] - summary['entities_before']['total'] != len(records):
        raise ValueError('Entity count arithmetic disagrees')


def emit(result):
    OUT.mkdir(parents=True, exist_ok=True)
    outputs = {}
    (OUT/'entities.jsonl').write_bytes(b''.join(json.dumps(e, separators=(',', ':'), allow_nan=False).encode() + b'\n'
                                                for e in result['entities']))
    outputs['entities.jsonl'] = sha256(OUT/'entities.jsonl')
    for name, value in (('fields.json', {'schema': SCHEMA, 'fields': result['records']}),
                        ('sources.json', {'schema': SCHEMA, 'sources': SOURCES, 'notes': NOTES}),
                        ('summary.json', {'schema': SCHEMA, **result['summary']}),
                        ('assumption_ledger.json', {'schema': SCHEMA, 'entries': NOTES['assumption_ledger'],
                                                    'target': 'data/derived/canonical/anatomy.json '
                                                              'assumption_ledger, on promotion'}),
                        ('manifest_fragment.json', {'schema': SCHEMA, 'structures': result['structures'],
                                                    'display_group': {'id': 'hair', 'name': 'Hair',
                                                                      'canonical_system': 'integumentary'},
                                                    'supersedes': NOTES['reconciliation']['supersedes'],
                                                    'retains': NOTES['reconciliation']['retains']})):
        outputs[name] = write_json(OUT/name, value)
    manifest = {'schema': SCHEMA, 'python': sys.version, 'platform': platform.platform(),
                'packages': {'numpy': np.__version__, 'scipy': __import__('scipy').__version__},
                'seed': SEED, 'sample_spacing_m': SAMPLE_SPACING_M, 'shell_sample_spacing_m': SHELL_SPACING_M,
                'shell_standoff_m': SHELL_STANDOFF_M,
                'inputs_sha256': {p: sha256(ROOT/p) for p in INPUTS},
                'extended_source_geometry_sha256': result['source_hashes'],
                'outputs_sha256': {**outputs, **{Path(p).name: h for p, h in result['geometry_hashes'].items()}},
                'geometry_sha256': result['geometry_hashes'],
                'canonical_assets_modified': False,
                'promotion': 'CANDIDATE. Nothing under data/derived/canonical or data/derived/app is written.'}
    write_json(OUT/'manifest.json', manifest)
    return manifest


def self_test():
    checks = []

    def check(name, ok, detail=''):
        checks.append({'check': name, 'pass': bool(ok), 'detail': detail})

    index = lexicon_index()
    check('lexicon claims each region once', len(index) == sum(len(v) for v in LEXICON.values()),
          f'{len(index)} regions')
    check('every field has evidence and morphology',
          set(FIELDS) == set(FIELD_EVIDENCE) == set(k for k in MORPHOLOGY if not k.startswith('_')),
          f'{len(FIELDS)} fields')
    for identity, evidence in FIELD_EVIDENCE.items():
        density = evidence['follicle_density_cm2']
        if density['value'] is not None:
            check(f'{identity} density cites a known source', density['source'] and
                  all(s in SOURCES for s in density['source']), str(density['source']))
            check(f'{identity} density tier is declared', density['tier'] in ('measured', 'transferred', 'assumed'),
                  density['tier'])
    for key, source in SOURCES.items():
        if source.get('pmid'):
            check(f'{key} PMID was verified against PubMed',
                  source.get('identifier_check', {}).get('status') == 'resolved',
                  str(source.get('identifier_check')))
    radius = np.array([1e-5, 3.5e-5])
    mass = strand_mass_kg(1e5, 3.5e-5, .3, 1300.)
    check('strand mass arithmetic', abs(mass - 1e5*np.pi*3.5e-5**2*.3*1300.) < 1e-15, f'{mass:.6f} kg')
    check('vellus is lighter than terminal per strand',
          strand_mass_kg(1, radius[0], .001, 1300.) < strand_mass_kg(1, radius[1], .3, 1300.))
    result = build(strict=False)
    total = sum(r['footprint']['area_m2'] for r in result['records'] if r['footprint']['resolved'])
    check('fields partition the envelope area', abs(total - 1.7812538727209595) < 1e-9, f'{total:.10f} m2')
    check('no field claims a face twice',
          sum(r['footprint']['envelope_face_count'] for r in result['records']) == 108716,
          str(sum(r['footprint']['envelope_face_count'] for r in result['records'])))
    check('scalp footprint is plausible for a scalp',
          .04 < next(r for r in result['records'] if r['id'] == 'scalp')['footprint']['area_m2'] < .12,
          f"{next(r for r in result['records'] if r['id'] == 'scalp')['footprint']['area_cm2']:.1f} cm2")
    check('glabrous field carries zero follicles',
          next(r for r in result['records'] if r['id'] == 'glabrous')['follicle_count']['value'] == 0)
    passed = sum(c['pass'] for c in checks)
    for entry in checks:
        print(('PASS ' if entry['pass'] else 'FAIL ') + entry['check'] + ('  ' + entry['detail'] if entry['detail'] else ''))
    print(f'self-test {passed}/{len(checks)}')
    return passed == len(checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        raise SystemExit(0 if self_test() else 1)
    result = build()
    manifest = emit(result)
    summary = result['summary']
    print(f"entities {summary['entities_before']['total']} -> {summary['entities_after']['total']}; "
          f"integumentary {summary['entities_before']['integumentary']} -> {summary['entities_after']['integumentary']}")
    print(f"envelope {summary['footprint']['envelope_area_m2']:.6f} m2 partitioned over "
          f"{summary['footprint']['fields_resolved']} fields, {summary['footprint']['fields_unresolved']} unresolved")
    for record in result['records']:
        area = record['footprint']['area_cm2']
        count = record['follicle_count']['value']
        mass = record['hair_mass_kg']['value']
        print(f"  {record['id']:14s} {('%9.1f cm2' % area) if area is not None else '  unresolved':>13s} "
              f"{('%11.0f' % count) if count is not None else '          -':>12s} follicles "
              f"{('%9.6f kg' % mass) if mass is not None else '          -':>13s} "
              f"[{record['follicle_density_cm2']['tier']}]")
    print(f"follicles {summary['population']['total_follicles']:.0f}; hair mass "
          f"{summary['population']['total_hair_mass_kg']:.6f} kg = "
          f"{summary['mass_ledger']['hair_mass_ppm_of_body_mass']:.1f} ppm of "
          f"{summary['mass_ledger']['body_mass_kg']:.4f} kg")
    print(f"density evidence measured {summary['evidence_split']['density_measured']} / transferred "
          f"{summary['evidence_split']['density_transferred']} / assumed {summary['evidence_split']['density_assumed']} "
          f"/ absent {summary['evidence_split']['density_absent']}")
    print('manifest ' + str((OUT/'manifest.json').relative_to(ROOT)) + ' inputs ' + str(len(manifest['inputs_sha256'])))


if __name__ == '__main__':
    main()
