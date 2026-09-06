#!/usr/bin/env python3
"""Prepare a fresh isolated skin-volume migration; never publish canonical data."""
from copy import deepcopy
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ihm.assembly.skin_layers import physical_skin_support

ANATOMY = 'data/derived/canonical/anatomy.json'
MECHANICS = 'data/derived/canonical/mechanics.json'
EVIDENCE = 'data/research/engineered_skin_territories/materialization.json'


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def encode(value):
    return (json.dumps(value, separators=(',', ':'), allow_nan=False)+'\n').encode()


def positive(value):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
        raise ValueError('Expected finite positive material quantity')
    return value


def close(a, b):
    if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12):
        raise ValueError('Retained mechanics differs from audited builder law')


def normalized(mechanics):
    rows = mechanics['entities']
    carriers = [e for e in rows if e['mass_role'] == 'numerical_boundary_carrier']
    raw = {e['id']: 0. if e in carriers else positive(e['volume_m3'])*positive(e['material']['density']['value']) for e in rows}
    total = positive(sum(raw.values()))
    target = positive(mechanics['mass_allocation']['target_mass_kg'])
    factor = positive((target-len(carriers)*1e-6)/total)
    for e in rows:
        if e['mass_role'] == 'numerical_boundary_carrier':
            e.update(mass_kg=1e-6, inertia_diagonal_kg_m2=[1e-8]*3, material_volume_m3=0.)
        else:
            mass = raw[e['id']]*factor
            ext = [max(b-a, .0001) for a, b in zip(e['bounds_m']['min'], e['bounds_m']['max'])]
            if len(ext) != 3 or not all(math.isfinite(x) for x in ext):
                raise ValueError('Invalid retained bounds')
            e.update(mass_kg=mass, inertia_diagonal_kg_m2=[mass*sum(x*x for x in ext)/18]*3, material_volume_m3=e['volume_m3'])
    byid = {e['id']: e for e in rows}
    for link in mechanics['links']:
        basis = link['parameter_basis']
        if basis['stiffness'] != 'EA/L' or basis['damping_ratio'] != .25:
            raise ValueError('Unsupported link law')
        close(link['stiffness_n_m'], positive(basis['young_modulus_pa'])*positive(basis['effective_area_m2'])/positive(basis['effective_length_m']))
        ma, mb = byid[link['a']]['mass_kg'], byid[link['b']]['mass_kg']
        link['damping_ns_m'] = .5*math.sqrt(link['stiffness_n_m']*ma*mb/(ma+mb))
    mechanics['mass_allocation'].update(unscaled_proxy_mass_kg=total, uniform_scale=factor,
        numerical_carrier_count=len(carriers), numerical_carrier_mass_kg=len(carriers)*1e-6)


def migrate_copy(anatomy_raw, mechanics_raw, geometry_raw, evidence_raw):
    anatomy, mechanics = json.loads(anatomy_raw), json.loads(mechanics_raw)
    if mechanics['source_files'].get(ANATOMY) != digest(anatomy_raw):
        raise ValueError('Retained mechanics anatomy receipt is stale')
    if json.loads(evidence_raw).get('anatomy_sha256') != digest(anatomy_raw):
        raise ValueError('Component evidence must describe the retained input anatomy')
    anatomy_ids = [e['id'] for e in anatomy['entities']]
    mechanics_ids = [e['id'] for e in mechanics['entities']]
    if len(set(anatomy_ids)) != len(anatomy_ids) or len(set(mechanics_ids)) != len(mechanics_ids) or set(anatomy_ids) != set(mechanics_ids):
        raise ValueError('Anatomy/mechanics entity identities differ')
    expected = deepcopy(mechanics)
    normalized(expected)
    for old, check in zip(mechanics['entities'], expected['entities']):
        close(old['mass_kg'], check['mass_kg'])
        if len(old['inertia_diagonal_kg_m2']) != 3:raise ValueError('Invalid inertia')
        for a, b in zip(old['inertia_diagonal_kg_m2'], check['inertia_diagonal_kg_m2']):close(a, b)
    for old, check in zip(mechanics['links'], expected['links']):close(old['damping_ns_m'], check['damping_ns_m'])
    for key in ('unscaled_proxy_mass_kg', 'uniform_scale', 'numerical_carrier_mass_kg', 'numerical_carrier_count'):
        close(mechanics['mass_allocation'][key], expected['mass_allocation'][key])
    skin, = [e for e in anatomy['entities'] if e['role'] == 'skin']
    support = physical_skin_support(skin['reference_geometry'], geometry_raw, evidence_raw)
    support['component_evidence_path'] = EVIDENCE
    skin['physical_surface_support'] = deepcopy(support)
    layers = sorted((e for e in anatomy['entities'] if e['role'] == 'skin_layer'), key=lambda e:e['shell']['depth_interval_m'][0])
    if {e['id'] for e in layers} != {'body-skin-epidermis', 'body-skin-dermis', 'body-skin-hypodermis'}:
        raise ValueError('Expected the three retained skin layers')
    byid = {e['id']: e for e in mechanics['entities']}
    depth = 0.
    changes = []
    for e in layers:
        shell = e['shell']; thickness = positive(shell['thickness_m'])
        close(shell['depth_interval_m'][0], depth)
        close(shell['depth_interval_m'][1], depth+thickness)
        if shell['geometry_materialized'] is not False or e['reference_geometry']['sha256'] != skin['reference_geometry']['sha256']:
            raise ValueError('Unsupported layer geometry')
        spec = byid[e['id']]
        close(spec['volume_m3'], e['volume_m3'])
        changes.append({'id':e['id'], 'old_volume_m3':e['volume_m3'], 'new_volume_m3':support['area_m2']*thickness, 'thickness_m':thickness})
        e.update(volume_m3=support['area_m2']*thickness,
                 volume_method='inferred exterior component area times assumed thickness; open-shell quadrature prior, not measured volume',
                 physical_surface_support=deepcopy(support))
        spec.update(volume_m3=e['volume_m3'], volume_basis=e['volume_method'], shell=deepcopy(shell), physical_surface_support=deepcopy(support))
        depth += thickness
    byid[skin['id']]['physical_surface_support'] = deepcopy(support)
    ledger, = [e for e in anatomy['assumption_ledger'] if e['id']=='SKIN-LAYER-PRIOR']
    ledger['statement'] = ('Layer quadrature uses an inferred exterior component of retained source skin and explicit nonoverlapping inward depth intervals. '
        'Uniform thickness is an engineering prior; exterior exclusivity and self-intersections remain unvalidated. Source geometry is preserved.')
    normalized(mechanics)
    candidate_anatomy = encode(anatomy)
    mechanics['source_files'][ANATOMY] = digest(candidate_anatomy)
    candidate_mechanics = encode(mechanics)
    manifest = dict(schema='ihm.skin-layer-migration-candidate.v1', published=False,
        inputs={ANATOMY:digest(anatomy_raw), MECHANICS:digest(mechanics_raw), EVIDENCE:digest(evidence_raw), skin['reference_geometry']['path']:digest(geometry_raw)},
        outputs={'anatomy.json':digest(candidate_anatomy), 'mechanics.json':digest(candidate_mechanics)},
        layers=changes, physical_surface_support=support, total_thickness_m=depth,
        old_mass_allocation=json.loads(mechanics_raw)['mass_allocation'], new_mass_allocation=mechanics['mass_allocation'],
        native_22body_mass_modified=False, geometry_modified=False,
        validation_scope='Retained metadata, exact skin geometry, builder normalization/inertia/link laws; no full geometry scan or native acceptance',
        logical_source_root='Candidate anatomy.json replaces the logical canonical anatomy path only within this candidate; it is not installed')
    return candidate_anatomy, candidate_mechanics, encode(manifest)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, required=True)
    args=parser.parse_args(); root=args.root.resolve(); output=args.output.resolve()
    if output == root/'data/derived/canonical' or (root/'data/derived/canonical') in output.parents:
        raise ValueError('Canonical output is forbidden')
    if output.exists():raise ValueError('Fresh candidate output directory required')
    a=(root/ANATOMY).read_bytes();m=(root/MECHANICS).read_bytes();e=(root/EVIDENCE).read_bytes()
    skin, = [r for r in json.loads(a)['entities'] if r['role']=='skin']
    for path, sha in json.loads(m)['source_files'].items():
        if digest((root/path).read_bytes()) != sha:raise ValueError('Stale retained mechanics source: '+path)
    g=(root/skin['reference_geometry']['path']).read_bytes()
    result=list(migrate_copy(a,m,g,e))
    manifest=json.loads(result[2])
    manifest['implementation_sources']={path:digest((root/path).read_bytes()) for path in ('scripts/prepare_skin_layer_migration.py','scripts/build_body_mechanics.py','ihm/assembly/skin_layers.py')}
    result[2]=encode(manifest)
    output.mkdir(parents=True)
    for name, raw in zip(('anatomy.json','mechanics.json','manifest.json'),result):(output/name).write_bytes(raw)
    print(json.dumps({'candidate':str(output),'published':False,'outputs':json.loads(result[2])['outputs']}))

if __name__=='__main__':main()
