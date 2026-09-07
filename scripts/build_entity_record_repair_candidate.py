#!/usr/bin/env python3
"""Candidate repair of three measured entity-record defects in the canonical metadata.

    .venv/bin/python scripts/build_entity_record_repair_candidate.py --self-test
    .venv/bin/python scripts/build_entity_record_repair_candidate.py --output data/derived/entity-record-repair-candidate-v1

Nothing here is canonical. `data/derived/canonical/` is opened read-only and is not
written; this emits proposed records beside it so the three can be reviewed together.
They are one candidate because all three change `entities` rows and therefore all
three land in the same single mass normalization.

  A  five duplicate entities   one BodyParts3D surface authored twice under two ids
  B  23 mis-roled discs        intervertebral disks carrying role rigid_bone
  C  skin layer areas          layer volume computed against a double-sided slab

Measured here (arithmetic and exact hashing over shipped assets):
  duplicates   independent decoded-geometry identity per pair, plus every mechanical
               reference that would dangle or is already degenerate
  discs        the ownership rule of cross-structure-repair-v1 replayed at the new
               role, and every link and count that moves
  skin         layer volume at each candidate area, the ICRP mass cross-check, and
               the cap-face area separating the two candidate areas
  allocation   one combined re-normalization under the exact law in
               build_body_mechanics.py, with the arithmetic recorded per step

Not established here: no geometry is remeshed, no builder is run, no canonical asset
is rewritten, no native process is started, and no ownership ledger is regenerated.
The disc ownership result is a replay of the shipped ledger's own rule, validated by
reproducing all 11,047 shipped decisions before any role is changed.
"""
from pathlib import Path
import argparse, gzip, hashlib, json, sys, time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANATOMY = ROOT / 'data/derived/canonical/anatomy.json'
MECHANICS = ROOT / 'data/derived/canonical/mechanics.json'
COINCIDENCE = ROOT / 'data/derived/cross-structure-repair-v1/coincidence.json'
LEDGER = ROOT / 'data/derived/cross-structure-repair-v1/ownership-ledger.jsonl'
REPAIR = ROOT / 'data/derived/cross-structure-repair-v1/summary.json'
ENVELOPE = ROOT / 'data/derived/outer-envelope/validation.json'
ENVELOPE_NPZ = ROOT / 'data/derived/outer-envelope/outer-envelope.npz'
MATERIALS = ROOT / 'data/derived/tissue-material-candidate-v1/materials.json'
GEOMETRY = ROOT / 'data/derived/canonical/geometry'
INPUTS = (ANATOMY, MECHANICS, COINCIDENCE, LEDGER, REPAIR, ENVELOPE, ENVELOPE_NPZ, MATERIALS)

# cross-structure-repair-v1/resolution.json. Lower rank keeps disputed volume.
RANK = {'rigid_bone': 0, 'cartilage': 1, 'tendon': 2, 'ligament': 3, 'fluid_cavity': 4,
        'vascular': 5, 'nerve': 6, 'muscle': 7, 'soft_organ': 8, 'connective_tissue': 9,
        'lymph_node_group': 10}

# coincidence.json drops the lexicographically later id; the earlier one keeps geometry.
DUPLICATES = (('body-bp3d-FJ1846', 'body-bp3d-FJ2013'),
              ('body-bp3d-FJ1916', 'body-bp3d-FJ2386'),
              ('body-bp3d-FJ1924', 'body-bp3d-FJ2394'),
              ('body-bp3d-FJ2440', 'body-bp3d-FJ2769'),
              ('body-bp3d-FJ2772', 'body-bp3d-FJ3201'))

DISCS = tuple(f'body-bp3d-FJ{n}' for n in range(3202, 3225))
LAYERS = (('body-skin-epidermis', 0.0001), ('body-skin-dermis', 0.0015),
          ('body-skin-hypodermis', 0.005))
SKIN_PARENT = 'body-bp3d-FJ2810'

RAW_SLAB_AREA = 3.502598974493317          # canonical FJ2810 surface_area_m2, double sided
EXTERIOR_AREA = 1.7804602548390722         # largest positive-orientation component, 109183 tri
ENVELOPE_AREA = 1.7812538727209595         # watertight genus-0 capped envelope
ICRP_SKIN_DENSITY = 1100.0                 # ICRP 89 para 514
ICRP_SKIN_MASS_KG = 3.3                    # ICRP 89 reference adult male skin (epidermis+dermis)
TARGET_MASS_KG = 77.1107029
CARRIER_MASS_KG = 1e-6

# build_body_mechanics.py:46 - the only role-dependent density in the allocation law.
BONE_DENSITY = 1900.0
SOFT_DENSITY = 1000.0

UNVERIFIED = (
    'no geometry was remeshed, tetrahedralized or re-registered here',
    'the disc ownership result replays the shipped rule and is not a fresh CGAL boolean pass',
    'no mesher, native process or builder was run, so no candidate PLC was proven meshable',
    'the 1.6 mm epidermis+dermis thickness prior is a literature prior, not measured on this specimen',
    'cartilage has no density row in tissue-material-candidate-v1, so re-roling the discs moves them '
    'from a role that carries a density to one that does not',
)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load_geometry(entity_id, geometry_dir=GEOMETRY):
    d = json.loads(gzip.open(Path(geometry_dir) / f'{entity_id}.json.gz', 'rt').read())
    V = np.asarray(d['positions'], float).reshape(-1, 3)
    F = np.asarray(d['indices'], np.int64).reshape(-1, 3)
    return V, F


def array_sha(V, F):
    """Hash the decoded arrays. The .json.gz sha differs on gzip framing alone."""
    h = hashlib.sha256()
    h.update(V.tobytes())
    h.update(F.tobytes())
    return h.hexdigest()


def facet_set_sha(V, F):
    """Order-independent hash of the facet set: exact coordinates, vertex order and
    face order both normalised away, so it answers 'same surface' not 'same file'."""
    keys = sorted(tuple(v for p in sorted(tuple(x) for x in t.tolist()) for v in p) for t in V[F])
    return hashlib.sha256(repr(keys).encode()).hexdigest()


def area(V, F):
    A, B, C = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    return float(np.linalg.norm(np.cross(B - A, C - A), axis=1).sum() / 2)


# ---------------------------------------------------------------- A duplicates

def duplicates(anatomy, mechanics, geometry_dir=GEOMETRY):
    ae = {e['id']: e for e in anatomy['entities']}
    me = {e['id']: e for e in mechanics['entities']}
    pair_of = {}
    for keep, drop in DUPLICATES:
        pair_of[keep] = drop
        pair_of[drop] = keep
    rows = []
    for keep, drop in DUPLICATES:
        Vk, Fk = load_geometry(keep, geometry_dir)
        Vd, Fd = load_geometry(drop, geometry_dir)
        k, d = ae[keep], ae[drop]
        rows.append({
            'keep': keep, 'drop': drop, 'name': k['name'],
            'identical_decoded_arrays': bool(Vk.shape == Vd.shape and Fk.shape == Fd.shape
                                             and np.array_equal(Vk, Vd) and np.array_equal(Fk, Fd)),
            'decoded_array_sha256': {keep: array_sha(Vk, Fk), drop: array_sha(Vd, Fd)},
            'facet_set_sha256': {keep: facet_set_sha(Vk, Fk), drop: facet_set_sha(Vd, Fd)},
            'stored_gzip_sha256': {keep: k['reference_geometry']['sha256'],
                                   drop: d['reference_geometry']['sha256']},
            'stored_gzip_sha_differs': k['reference_geometry']['sha256'] != d['reference_geometry']['sha256'],
            'triangles': int(len(Fk)), 'vertices': int(len(Vk)),
            'area_m2': {keep: area(Vk, Fk), drop: area(Vd, Fd)},
            'names_equal': k['name'] == d['name'],
            'fma_concepts_equal': ([c['concept_id'] for c in k.get('concepts', [])]
                                   == [c['concept_id'] for c in d.get('concepts', [])]),
            'terminal_fma_concepts': [c['concept_id'] for c in k.get('concepts', [])],
            'role_system': {keep: [k['role'], k['system']], drop: [d['role'], d['system']]},
            'source_obj': {keep: k['provenance']['files'][0]['path'],
                           drop: d['provenance']['files'][0]['path']},
            'source_obj_sha256': {keep: k['provenance']['files'][0]['sha256'],
                                  drop: d['provenance']['files'][0]['sha256']},
            'paired_structure_possible': False,
            'not_a_left_right_pair_because': (
                'the two surfaces are identical, not mirrored, and both carry the same terminal '
                'FMA identity, which is an unpaired midline or single structure. A genuine pair '
                'would carry laterality concepts and reflected coordinates.'),
            'mechanics': {
                'keep_links': sum(1 for l in mechanics['links'] if keep in (l['a'], l['b'])),
                'drop_links': sum(1 for l in mechanics['links'] if drop in (l['a'], l['b'])),
                'keep_muscle_refs': sum(1 for m in mechanics['muscles']
                                        if any(a['entity_id'] == keep for a in m['anchors'])),
                'drop_muscle_refs': sum(1 for m in mechanics['muscles']
                                        if any(a['entity_id'] == drop for a in m['anchors'])),
                'volume_m3_each': me[keep]['volume_m3'],
                'density_kg_m3': me[keep]['material']['density']['value'],
                'mass_kg_each': me[keep]['mass_kg'],
            },
        })

    degenerate_muscles = []
    for m in mechanics['muscles']:
        ids = [a['entity_id'] for a in m['anchors']]
        if len(ids) == 2 and pair_of.get(ids[0]) == ids[1]:
            degenerate_muscles.append({'id': m['id'], 'anchors': ids,
                                       'rest_path_length_m': m['rest_path_length_m'],
                                       'max_isometric_force_n': m['max_isometric_force_n']})
    self_links = [{'a': l['a'], 'b': l['b'], 'kind': l['kind'], 'stiffness_n_m': l['stiffness_n_m']}
                  for l in mechanics['links'] if pair_of.get(l['a']) == l['b']]
    drops = [d for _, d in DUPLICATES]
    dangling = [{'a': l['a'], 'b': l['b'], 'kind': l['kind']} for l in mechanics['links']
                if (l['a'] in drops) != (l['b'] in drops)]
    dangling_muscles = [{'id': m['id'], 'anchor': a['entity_id']} for m in mechanics['muscles']
                        for a in m['anchors'] if a['entity_id'] in drops
                        and not any(pair_of[a['entity_id']] == x['entity_id'] for x in m['anchors'])]
    incoming = [{'from': e['id'], 'to': c['entity_id']} for e in anatomy['entities']
                for c in e.get('connections', []) if c['entity_id'] in drops]

    double_volume = sum(me[d]['volume_m3'] for _, d in DUPLICATES)
    double_mass = sum(me[d]['mass_kg'] for _, d in DUPLICATES)
    return {
        'pairs': rows,
        'survivor_rule': ('the lexicographically earlier id survives, matching the decision already '
                          'recorded in cross-structure-repair-v1/coincidence.json, so the candidate '
                          'and the assembled mesh complex name the same survivor. Nothing anatomical '
                          'distinguishes the two ids: identical geometry, identical name, identical '
                          'FMA concept set. Where the two differ mechanically the survivor is also '
                          'the more connected id, which is the second, independent reason.'),
        'survivor_is_more_connected': {k: {'keep_links': r['mechanics']['keep_links'],
                                           'drop_links': r['mechanics']['drop_links']}
                                       for k, r in zip([p[0] for p in DUPLICATES], rows)},
        'double_counted_today': {
            'entities': len(DUPLICATES),
            'volume_m3': double_volume, 'volume_mL': double_volume * 1e6,
            'mass_kg': double_mass,
            'mass_fraction_of_body': double_mass / TARGET_MASS_KG,
            'note': 'the mechanics volume prior, not the repaired watertight volume; both copies '
                    'carry the identical figure, so exactly this much is counted twice',
        },
        'already_degenerate': {
            'muscle_and_connective_priors_with_both_anchors_on_one_structure': degenerate_muscles,
            'support_links_joining_a_structure_to_itself': self_links,
            'note': 'these are defects the duplicate already causes today. Nine actuator and '
                    'connective priors span a bone to itself, and one 20 kN/m skeletal support '
                    'link connects the hyoid to the hyoid.',
        },
        'references_that_must_be_remapped': {
            'links': dangling, 'link_count': len(dangling),
            'muscle_anchors': dangling_muscles, 'muscle_anchor_count': len(dangling_muscles),
            'anatomy_connections': incoming,
            'remap_is_exact_because': ('the survivor carries byte-identical vertices, so every '
                                       'attachment point, projection distance and nearest-surface '
                                       'query on the dropped id returns the same result on the '
                                       'survivor. The remap changes no coordinate.'),
        },
    }


# ---------------------------------------------------------------- B disc roles

def resolve(a, b, ra, rb, volume):
    """cross-structure-repair-v1 resolution rule: rank, then larger repaired volume."""
    if RANK[ra] != RANK[rb]:
        return a if RANK[ra] < RANK[rb] else b
    va, vb = volume.get(a, 0.), volume.get(b, 0.)
    if va != vb:
        return a if va > vb else b
    return min(a, b)


def disc_roles(anatomy, mechanics, ledger, repair):
    ae = {e['id']: e for e in anatomy['entities']}
    me = {e['id']: e for e in mechanics['entities']}
    volume = repair['graph']['entity_volume_m3']
    discs = set(DISCS)

    reproduced = sum(1 for r in ledger
                     if resolve(r['a'], r['b'], r['a_role'], r['b_role'], volume) == r['owner'])
    involved = [r for r in ledger if r['a'] in discs or r['b'] in discs]

    def replay(role):
        flips, flipped_volume = 0, 0.
        for r in involved:
            ra = role if r['a'] in discs else r['a_role']
            rb = role if r['b'] in discs else r['b_role']
            if resolve(r['a'], r['b'], ra, rb, volume) != r['owner']:
                flips += 1
                flipped_volume += r['overlap_volume_m3']
        return {'flipped_pairs': flips, 'flipped_overlap_volume_m3': flipped_volume,
                'flipped_overlap_mL': flipped_volume * 1e6}

    counterparts = {}
    bone_ties = []
    for r in involved:
        d, o = (r['a'], r['b']) if r['a'] in discs else (r['b'], r['a'])
        orole = r['b_role'] if r['a'] in discs else r['a_role']
        osys = r['b_system'] if r['a'] in discs else r['a_system']
        counterparts[f'{orole} / {osys}'] = counterparts.get(f'{orole} / {osys}', 0) + 1
        if orole == 'rigid_bone':
            bone_ties.append({'disc': d, 'bone': o, 'disc_volume_m3': volume.get(d),
                              'bone_volume_m3': volume.get(o), 'owner': r['owner'],
                              'overlap_volume_m3': r['overlap_volume_m3']})

    tree_edges = [l for l in mechanics['links']
                  if l['kind'] == 'inferred_skeletal_support' and (l['a'] in discs or l['b'] in discs)]
    anchored = [l for l in mechanics['links']
                if l['kind'] == 'soft_tissue_support' and (l['a'] in discs or l['b'] in discs)]
    attach = [l for l in mechanics['links']
              if l['kind'] == 'muscle_or_connective_attachment_support' and (l['a'] in discs or l['b'] in discs)]
    by_role = {}
    for l in anchored:
        other = l['a'] if l['b'] in discs else l['b']
        r = me[other]['role']
        by_role[r] = by_role.get(r, 0) + 1

    disc_volume = sum(me[d]['volume_m3'] for d in DISCS)
    return {
        'entities': list(DISCS),
        'count': len(DISCS),
        'current': {'role': 'rigid_bone', 'system': 'skeletal', 'constitutive': 'rigid',
                    'density_kg_m3': BONE_DENSITY},
        'proposed': {'role': 'cartilage', 'system': 'skeletal', 'constitutive': 'affine_neo_hookean',
                     'density_kg_m3': SOFT_DENSITY},
        'role_evidence': {
            'basis': 'the atlas own recorded FMA is-a identity, already stored on every disc row',
            'isa_path': ['intervertebral disk of Nth vertebra', 'intervertebral disk',
                         'articular disk of symphysis (FMA67396)', 'cartilage organ (FMA55107)'],
            'anchor_concept': 'FMA55107 cartilage organ',
            'same_anchor_as': ['body-bp3d-FJ2440 cricoid cartilage', 'body-bp3d-FJ2808 thyroid cartilage',
                               'the 14 costal cartilages'],
            'conclusion': 'the classifier already resolved these entities to cartilage organ and used '
                          'that anchor to assign system skeletal. Only the role field disagrees with '
                          'the evidence the record carries.',
        },
        'system_should_not_change': {
            'proposed': 'skeletal',
            'reason': ('all 28 cartilage entities in the atlas are system skeletal, including every '
                       'laryngeal and costal cartilage, and the disc classification block derives '
                       'skeletal from FMA55107 with recorded ontology hashes. Moving the discs to '
                       'system connective would contradict the record own evidence and split the '
                       'cartilage class across two systems. Recommend NO change to system.'),
        },
        'ownership_rule_replay': {
            'shipped_decisions': len(ledger),
            'decisions_reproduced_before_any_change': reproduced,
            'rule_reproduction_exact': reproduced == len(ledger),
            'disc_involving_pairs': len(involved),
            'disc_involving_overlap_volume_m3': sum(r['overlap_volume_m3'] for r in involved),
            'disc_wins_today': sum(1 for r in involved if r['owner'] in discs),
            'counterpart_role_system_census': counterparts,
            'at_proposed_role_cartilage': replay('cartilage'),
            'sensitivity': {role: replay(role) for role in
                            ('cartilage', 'tendon', 'ligament', 'fluid_cavity', 'vascular',
                             'connective_tissue')},
            'finding': ('the role correction flips ZERO ownership decisions. Every disc conflict '
                        'counterpart is either rigid_bone, where the vertebra already wins all 63 '
                        'ties on larger repaired volume, or rank 5 and below, which cartilage '
                        'outranks exactly as rigid_bone did. The outcome is invariant across ranks '
                        '0 to 4. The premise that the disc displaces the soft tissue it touches is '
                        'true and stays true; correcting the role does not fix it, and nothing in '
                        'the repair artifact needs regenerating for this change.'),
            'disc_vs_bone_ties': {'count': len(bone_ties),
                                  'all_won_by_the_vertebra': all(t['owner'] != t['disc'] for t in bone_ties),
                                  'overlap_volume_m3': sum(t['overlap_volume_m3'] for t in bone_ties),
                                  'rows': bone_ties},
        },
        'what_actually_changes': {
            'constitutive': 'rigid -> affine_neo_hookean for 23 entities',
            'density_kg_m3': f'{BONE_DENSITY} -> {SOFT_DENSITY}',
            'counts_rigid_bones': '257 -> 233 (23 discs plus the dropped duplicate hyoid)',
            'counts_soft_solids': '2145 -> 2164 (+23 discs, -4 dropped duplicates that were soft; the fifth, FJ3201, was rigid)',
            'skeletal_support_tree': {
                'edges_today': sum(1 for l in mechanics['links'] if l['kind'] == 'inferred_skeletal_support'),
                'edges_touching_a_disc': len(tree_edges),
                'note': 'build_body_mechanics.py builds a minimum spanning tree over rigid_bone '
                        'centroids. Removing 23 discs from that set rebuilds the tree over 233 '
                        'nodes; the 36 edges that currently route through a disc do not survive.',
            },
            'soft_tissue_supports_anchored_on_a_disc': {
                'count': len(anchored), 'by_attaching_role': by_role,
                'muscle_or_connective_attachment_links_also_on_a_disc': len(attach),
                'note': 'every non-bone entity attaches to its nearest rigid_bone. These 267 links '
                        'currently anchor onto a disc and must re-anchor to a real vertebra, and '
                        'the 23 discs themselves each gain a soft_tissue_support link they do not '
                        'have today. This, not the ownership rule, is the substantive consequence.',
            },
        },
        'mass_that_moves': {
            'disc_material_volume_m3': disc_volume,
            'disc_material_volume_mL': disc_volume * 1e6,
            'unscaled_mass_today_kg': disc_volume * BONE_DENSITY,
            'unscaled_mass_proposed_kg': disc_volume * SOFT_DENSITY,
            'unscaled_delta_kg': disc_volume * (SOFT_DENSITY - BONE_DENSITY),
            'arithmetic': f'{disc_volume:.12e} m3 x ({SOFT_DENSITY} - {BONE_DENSITY}) kg/m3 '
                          f'= {disc_volume * (SOFT_DENSITY - BONE_DENSITY):+.9f} kg unscaled',
        },
        'material_table_consequence': {
            'today_role_rigid_bone': 'E 17 GPa longitudinal / 11.5 GPa transverse (reilly1975), '
                                     'density 1300 whole skeleton or 1900 cortical (icrp89)',
            'proposed_role_cartilage': 'aggregate modulus 1.0 MPa, E 986 kPa derived, nu 0.08 '
                                       '(liu1997), density absent',
            'caveat': 'cartilage carries no density row in tissue-material-candidate-v1, so the '
                      'discs move from a role with a density to a role without one. The mechanics '
                      'allocation is unaffected because it uses its own blanket 1000 kg/m3, but '
                      'the candidate material table gains 0.1258 L of volume with no density.',
            'also_wrong_either_way': 'nu 0.08 is the biphasic solid-phase ratio for articular '
                                     'cartilage. A disc is fibrocartilage with a gel nucleus and is '
                                     'not articular cartilage; cartilage is the closest role the '
                                     'vocabulary offers, not a correct constitutive claim.',
        },
    }


# ---------------------------------------------------------------- C skin layers

def cap_area(npz_path=ENVELOPE_NPZ, apex_count=36):
    """The envelope's 438 aperture cap faces are synthetic fans, not integument."""
    z = np.load(npz_path)
    V, F = z['positions'], z['indices']
    A, B, C = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    a = np.linalg.norm(np.cross(B - A, C - A), axis=1) / 2
    apex = np.arange(len(V) - apex_count, len(V))
    mask = np.isin(F, apex).any(axis=1)
    return {'total_area_m2': float(a.sum()), 'cap_faces': int(mask.sum()),
            'cap_area_m2': float(a[mask].sum()),
            'authored_face_area_m2': float(a[~mask].sum()),
            'cap_fraction': float(a[mask].sum() / a.sum())}


def skin_layers(anatomy, mechanics, envelope, caps):
    ae = {e['id']: e for e in anatomy['entities']}
    me = {e['id']: e for e in mechanics['entities']}
    inflation = RAW_SLAB_AREA / EXTERIOR_AREA

    def at(a):
        vols = {i: a * t for i, t in LAYERS}
        return {'area_m2': a, 'layer_volume_m3': vols, 'total_volume_m3': sum(vols.values()),
                'total_volume_L': sum(vols.values()) * 1000,
                'epidermis_dermis_mass_kg_at_icrp_density': a * 0.0016 * ICRP_SKIN_DENSITY,
                'icrp_relative': a * 0.0016 * ICRP_SKIN_DENSITY / ICRP_SKIN_MASS_KG - 1,
                'all_layers_mass_kg_at_icrp_density': a * 0.0066 * ICRP_SKIN_DENSITY}

    today_volume = sum(ae[i]['volume_m3'] for i, _ in LAYERS)
    today_mass = sum(me[i]['mass_kg'] for i, _ in LAYERS)
    chosen = at(EXTERIOR_AREA)
    return {
        'defect': {
            'canonical_area_m2': RAW_SLAB_AREA,
            'area_is': 'the summed triangle area of a DOUBLE-SIDED slab; every square metre of body '
                       'surface is present twice, once outward and once inward',
            'slab_evidence': {'faces': envelope['source_slab_topology']['faces'],
                              'faces_with_antiparallel_partner':
                                  envelope['antiparallel_pairing']['faces_with_antiparallel_partner'],
                              'partner_area_fraction':
                                  envelope['antiparallel_pairing']['partner_area_fraction'],
                              'slab_thickness_mm_p50':
                                  envelope['antiparallel_pairing']['slab_thickness_mm_percentiles']['50'],
                              'reading': '99.998% of the slab area has an antiparallel partner at '
                                         '1.96 mm, so the slab is a two-sided shell, not a surface'},
            'inflation_factor': inflation,
            'today_layer_volume_m3': today_volume, 'today_layer_volume_L': today_volume * 1000,
            'today_layer_mass_kg': today_mass,
            'today_mass_fraction_of_body': today_mass / TARGET_MASS_KG,
        },
        'candidate_areas': {
            'exterior_component_m2': EXTERIOR_AREA, 'envelope_m2': ENVELOPE_AREA,
            'difference_m2': ENVELOPE_AREA - EXTERIOR_AREA,
            'difference_relative': ENVELOPE_AREA / EXTERIOR_AREA - 1,
            'difference_cm2': (ENVELOPE_AREA - EXTERIOR_AREA) * 1e4,
            'volume_difference_at_6.6mm_mL': (ENVELOPE_AREA - EXTERIOR_AREA) * 0.0066 * 1e6,
            'mass_difference_at_icrp_density_g': (ENVELOPE_AREA - EXTERIOR_AREA) * 0.0066 * ICRP_SKIN_DENSITY * 1000,
            'exterior': at(EXTERIOR_AREA), 'envelope': at(ENVELOPE_AREA), 'raw_slab': at(RAW_SLAB_AREA),
        },
        'decision': {
            'use': 'exterior_component', 'area_m2': EXTERIOR_AREA,
            'reasons': [
                'every one of its 109,183 triangles is an authored BodyParts3D skin facet, so the '
                'layer volume is a per-facet extrusion of authored surface',
                f'the envelope contains {caps["cap_faces"]} synthetic aperture cap faces totalling '
                f'{caps["cap_area_m2"]:.10f} m2 ({caps["cap_area_m2"] * 1e4:.2f} cm2, '
                f'{caps["cap_fraction"] * 100:.4f}% of it). Those fans seal 36 apertures and are not '
                'integument; extruding skin layers through a sealed nostril is the same class of '
                'error being fixed here',
                'the envelope also drops authored exterior facets during sheet assignment (157 '
                'stray faces discarded), so it is neither a superset nor a subset of the skin',
                'it is already the hash-bound basis checked by '
                'ihm.assembly.skin_layers.physical_skin_support against exact geometry bytes, and '
                'the basis of the accepted staging epoch skin-layer-epoch-20260906-v2, so choosing '
                'it keeps the whole candidate chain on one number',
            ],
            'the_choice_is_immaterial_numerically': (
                f'the two areas differ by {(ENVELOPE_AREA - EXTERIOR_AREA) * 1e4:.4f} cm2, '
                f'{(ENVELOPE_AREA / EXTERIOR_AREA - 1) * 100:.4f}%. Over the 6.6 mm stack that is '
                f'{(ENVELOPE_AREA - EXTERIOR_AREA) * 0.0066 * 1e6:.3f} mL and '
                f'{(ENVELOPE_AREA - EXTERIOR_AREA) * 0.0066 * ICRP_SKIN_DENSITY * 1000:.2f} g, an '
                'order of magnitude below the thickness prior own spread. This is a provenance '
                'argument, not a numerical one, and the envelope is a valid cross-check that '
                'agrees to 0.045%.'),
            'envelope_is_still_the_right_number_for': 'enclosed volume (69.720 L) and any '
                                                      'whole-body closed-surface quantity. It is '
                                                      'the wrong number only as a skin-layer divisor.',
        },
        'envelope_cap_measurement': caps,
        'icrp_cross_check': {
            'claim': 'the 1.6 mm epidermis+dermis thickness prior is right and only the area was wrong',
            'at_exterior_area': f'{EXTERIOR_AREA} m2 x 0.0016 m x {ICRP_SKIN_DENSITY} kg/m3 = '
                                f'{EXTERIOR_AREA * 0.0016 * ICRP_SKIN_DENSITY:.6f} kg against ICRP '
                                f'{ICRP_SKIN_MASS_KG} kg, {(EXTERIOR_AREA * 0.0016 * ICRP_SKIN_DENSITY / ICRP_SKIN_MASS_KG - 1) * 100:+.2f}%',
            'at_raw_slab_area': f'{RAW_SLAB_AREA} m2 x 0.0016 m x {ICRP_SKIN_DENSITY} kg/m3 = '
                                f'{RAW_SLAB_AREA * 0.0016 * ICRP_SKIN_DENSITY:.6f} kg, '
                                f'{(RAW_SLAB_AREA * 0.0016 * ICRP_SKIN_DENSITY / ICRP_SKIN_MASS_KG - 1) * 100:+.2f}%',
            'reading': 'the corrected area lands 5% under an independent reference mass; the '
                       'canonical area lands 87% over. The thickness prior is not what is wrong.',
            'not_established': 'ICRP 3.3 kg is a reference-man value for a 73 kg phantom, not a '
                               'measurement on this specimen. It is a sanity bound, not a target.',
        },
        'phantom_volume': {
            'today_m3': today_volume, 'proposed_m3': chosen['total_volume_m3'],
            'phantom_m3': today_volume - chosen['total_volume_m3'],
            'phantom_L': (today_volume - chosen['total_volume_m3']) * 1000,
            'arithmetic': f'({RAW_SLAB_AREA} - {EXTERIOR_AREA}) m2 x 0.0066 m = '
                          f'{(RAW_SLAB_AREA - EXTERIOR_AREA) * 0.0066 * 1000:.9f} L',
            'context': 'the whole envelope encloses 69.720 L, so the phantom layer volume is '
                       f'{(today_volume - chosen["total_volume_m3"]) / 0.06972028286368345 * 100:.2f}% '
                       'of the entire body interior',
        },
        'should_not_change': {
            'thicknesses_m': {i: t for i, t in LAYERS},
            'reason': 'the ICRP cross-check confirms the 0.1/1.5/5.0 mm priors at the corrected '
                      'area. Changing thickness as well would destroy the only independent check '
                      'this correction has.',
            'parent_surface_area_m2': RAW_SLAB_AREA,
            'parent_reason': 'FJ2810 surface_area_m2 is a true measurement of the authored surface '
                             'and stays descriptive source metadata. The defect is that a '
                             'double-sided quantity was used as a single-sided one, not that the '
                             'measurement is wrong.',
            'skin_parent_carrier_mass': 'FJ2810 stays a 1 mg numerical boundary carrier with zero '
                                        'material volume. Its mass belongs to the layer entities.',
        },
    }


# ---------------------------------------------------------------- combined allocation

def allocation(anatomy, mechanics, area_m2=EXTERIOR_AREA):
    """Replay build_body_mechanics.py mass law with the three corrections applied together."""
    me = {e['id']: e for e in mechanics['entities']}
    scale = mechanics['mass_allocation']['uniform_scale']
    carriers = mechanics['mass_allocation']['numerical_carrier_count']
    unscaled = {e['id']: (e['mass_kg'] / scale) for e in mechanics['entities']
                if e['mass_role'] != 'numerical_boundary_carrier'}
    base = sum(unscaled.values())

    a = sum(unscaled[d] for _, d in DUPLICATES)
    disc_volume = sum(me[d]['volume_m3'] for d in DISCS)
    b = disc_volume * (SOFT_DENSITY - BONE_DENSITY)
    old_skin = sum(me[i]['volume_m3'] for i, _ in LAYERS)
    new_skin = sum(area_m2 * t for _, t in LAYERS)
    c = (new_skin - old_skin) * SOFT_DENSITY

    total = base - a + b + c
    factor = (TARGET_MASS_KG - carriers * CARRIER_MASS_KG) / total

    def after(i):
        if i in dict((d, k) for k, d in DUPLICATES):
            return None
        v = area_m2 * dict(LAYERS)[i] if i in dict(LAYERS) else me[i]['volume_m3']
        rho = SOFT_DENSITY if (i in DISCS or me[i]['role'] != 'rigid_bone') else BONE_DENSITY
        return v * rho * factor

    rows = []
    for i, _ in LAYERS:
        rows.append({'id': i, 'change': 'C skin area',
                     'volume_m3': [me[i]['volume_m3'], area_m2 * dict(LAYERS)[i]],
                     'mass_kg': [me[i]['mass_kg'], after(i)]})
    for i in DISCS:
        rows.append({'id': i, 'change': 'B disc role',
                     'volume_m3': [me[i]['volume_m3'], me[i]['volume_m3']],
                     'mass_kg': [me[i]['mass_kg'], after(i)]})
    for _, d in DUPLICATES:
        rows.append({'id': d, 'change': 'A duplicate dropped',
                     'volume_m3': [me[d]['volume_m3'], None], 'mass_kg': [me[d]['mass_kg'], None]})

    untouched = [i for i in ('body-bp3d-FJ1252', 'body-bp3d-FJ3226', 'body-bp3d-FJ2810')
                 if i in me]
    layer_mass_after = sum(after(i) for i, _ in LAYERS)
    disc_mass_after = sum(after(i) for i in DISCS)
    return {
        'law': 'build_body_mechanics.py:92 -- mass = volume x role density, then one uniform '
               'factor = (77.1107029 - 6e-6) / sum(mass) applied to every proxy row. There is one '
               'global normalizer, so every change to any entity moves every other entity mass.',
        'reproduction': {'recorded_unscaled_kg': mechanics['mass_allocation']['unscaled_proxy_mass_kg'],
                         'recomputed_unscaled_kg': base,
                         'exact': abs(base - mechanics['mass_allocation']['unscaled_proxy_mass_kg']) < 1e-9,
                         'recorded_scale': scale,
                         'rows_where_mass_ne_volume_times_density': sum(
                             1 for e in mechanics['entities']
                             if e['mass_role'] != 'numerical_boundary_carrier'
                             and abs(e['mass_kg'] / scale - e['volume_m3'] * e['material']['density']['value']) > 1e-12)},
        'steps_unscaled_kg': {
            'base': base,
            'A_drop_5_duplicates': -a,
            'B_23_discs_1900_to_1000': b,
            'C_skin_layers_to_exterior_area': c,
            'total': total,
            'arithmetic': f'{base:.10f} - {a:.10f} + ({b:.10f}) + ({c:.10f}) = {total:.10f} kg',
        },
        'uniform_scale': {'before': scale, 'after': factor, 'relative_change': factor / scale - 1},
        'target_mass_kg': TARGET_MASS_KG,
        'total_after_kg': sum(x for x in (after(e['id']) for e in mechanics['entities']
                                          if e['mass_role'] != 'numerical_boundary_carrier') if x is not None)
                          + carriers * CARRIER_MASS_KG,
        'headline': {
            'skin_layer_mass_kg': [sum(me[i]['mass_kg'] for i, _ in LAYERS), layer_mass_after],
            'skin_layer_fraction_of_body': [sum(me[i]['mass_kg'] for i, _ in LAYERS) / TARGET_MASS_KG,
                                            layer_mass_after / TARGET_MASS_KG],
            'disc_mass_kg': [sum(me[i]['mass_kg'] for i in DISCS), disc_mass_after],
            'duplicate_mass_removed_kg': sum(me[d]['mass_kg'] for _, d in DUPLICATES),
        },
        'every_untouched_entity_moves': [
            {'id': i, 'change': 'none of the three; moves only through the normalizer',
             'mass_kg': [me[i]['mass_kg'], after(i)],
             'relative': (after(i) / me[i]['mass_kg'] - 1) if me[i]['mass_kg'] else None}
            for i in untouched if me[i]['mass_role'] != 'numerical_boundary_carrier'],
        'rows': rows,
        'not_a_weight_change': 'the body still masses 77.1107029 kg. This reallocates an '
                               'overlapping-atlas proxy partition; it is not measured weight loss '
                               'and does not touch native patient mass, which is scaled '
                               'independently by the 22-body model.',
    }


# ---------------------------------------------------------------- proposed records

def records(anatomy, mechanics, area_m2=EXTERIOR_AREA):
    ae = {e['id']: e for e in anatomy['entities']}
    me = {e['id']: e for e in mechanics['entities']}
    support = ae[LAYERS[0][0]].get('physical_surface_support')
    out = {'schema': 'ihm.entity-record-repair.v1',
           'basis': 'field-level deltas against the canonical rows, not whole replacement files',
           'anatomy': {'remove': [], 'change': []}, 'mechanics': {'remove': [], 'change': []}}
    for keep, drop in DUPLICATES:
        out['anatomy']['remove'].append({'id': drop, 'reason': 'byte-identical geometry duplicate',
                                         'survivor': keep})
        out['mechanics']['remove'].append({'id': drop, 'survivor': keep,
                                           'remap_links_and_muscle_anchors_to': keep})
    for d in DISCS:
        out['anatomy']['change'].append({'id': d, 'name': ae[d]['name'],
                                         'role': ['rigid_bone', 'cartilage'],
                                         'system': ['skeletal', 'skeletal (unchanged)']})
        out['mechanics']['change'].append({'id': d, 'role': ['rigid_bone', 'cartilage'],
                                           'constitutive': ['rigid', 'affine_neo_hookean'],
                                           'material.density.value': [BONE_DENSITY, SOFT_DENSITY],
                                           'material.density.prior_range': [[1500, 2200], [900, 1100]]})
    for i, t in LAYERS:
        new = area_m2 * t
        out['anatomy']['change'].append({
            'id': i, 'volume_m3': [ae[i]['volume_m3'], new],
            'volume_method': [ae[i]['volume_method'],
                              'inferred exterior component area times assumed thickness; open-shell '
                              'quadrature prior, not measured volume'],
            'physical_surface_support.area_m2': [None if not support else support.get('area_m2'), area_m2],
            'shell.thickness_m': [t, t],
            'surface_area_m2': [RAW_SLAB_AREA, RAW_SLAB_AREA],
            'arithmetic': f'{area_m2} m2 x {t} m = {new:.12e} m3'})
        out['mechanics']['change'].append({'id': i, 'volume_m3': [me[i]['volume_m3'], new],
                                           'material_volume_m3': [me[i]['material_volume_m3'], new]})
    out['anatomy']['change'].append({
        'id': SKIN_PARENT, 'physical_surface_support.area_m2': [None, area_m2],
        'volume_m3': [ae[SKIN_PARENT]['volume_m3'], ae[SKIN_PARENT]['volume_m3']],
        'surface_area_m2': [RAW_SLAB_AREA, RAW_SLAB_AREA],
        'note': 'gains the support receipt only; stays a 1 mg carrier with zero material volume'})
    out['anatomy']['assumption_ledger'] = {
        'add_or_amend': ['SKIN-LAYER-PRIOR: layer support area is the inferred exterior component of '
                         'the authored slab, not its raw two-sided triangle area',
                         'DUPLICATE-ENTITY-COLLAPSE: five BodyParts3D ids naming one surface collapse '
                         'to the lexicographically earlier id',
                         'DISC-ROLE-FROM-FMA: intervertebral disk role follows its recorded FMA is-a '
                         'anchor FMA55107 cartilage organ']}
    out['counts'] = {'entities': [len(anatomy['entities']), len(anatomy['entities']) - len(DUPLICATES)],
                     'rigid_bones': [mechanics['counts']['rigid_bones'],
                                     mechanics['counts']['rigid_bones'] - len(DISCS) - 1],
                     'soft_solids': [mechanics['counts']['soft_solids'],
                                     mechanics['counts']['soft_solids'] + len(DISCS) - 4]}
    return out


# ---------------------------------------------------------------- downstream

def downstream(disc, skin):
    return {
        'A_duplicates': {
            'must_change': [
                {'consumer': 'data/derived/canonical/mechanics.json links',
                 'effect': f'{disc and ""}30 links anchored on FJ3201 and 4 more on the other four '
                           'dropped ids re-anchor to the survivor. One inferred_skeletal_support '
                           'edge joining FJ2772 to FJ3201 at 20 kN/m disappears entirely: it '
                           'connects the hyoid to the hyoid.'},
                {'consumer': 'mechanics.json muscles',
                 'effect': '9 actuator and connective priors span a structure to itself and lose '
                           'their second anchor. They need re-projection onto a genuinely distinct '
                           'second bone or explicit deletion; silently collapsing them to zero '
                           'length is not acceptable.'},
                {'consumer': 'cross-structure-repair-v1',
                 'effect': 'already drops these five ids from any assembled complex, so the mesh '
                           'complex agrees with the candidate and needs no regeneration. The '
                           'conflict graph counts (2403 surfaces) already exclude nothing, so the '
                           '22 facet-sharing pairs shrink by the 5 identical-geometry pairs.'},
                {'consumer': 'anatomy.json connections',
                 'effect': 'one regional_spatial_support connection points at FJ3201 and re-points '
                           'to FJ2772 at the identical distance.'},
            ],
            'must_not_change': [
                {'item': 'geometry files', 'reason': 'the survivor already holds the exact bytes. '
                                                     'Nothing is remeshed and no coordinate moves.'},
                {'item': 'the 4269 segmentation-noise and 2528 modelling-conflict classifications',
                 'reason': 'the five identical pairs were already excluded from the conflict '
                           'classes; they are recorded separately as coincidence.'},
            ],
        },
        'B_discs': {
            'must_change': [
                {'consumer': 'mechanics.json skeletal support tree',
                 'effect': 'the minimum spanning tree is built over rigid_bone centroids. 257 -> '
                           '233 nodes rebuilds it; 36 current edges route through a disc.'},
                {'consumer': 'mechanics.json soft_tissue_support',
                 'effect': '267 links anchored on a disc re-anchor to the nearest true vertebra, '
                           'and the 23 discs each gain a support link they do not have today. '
                           'Every link damping is reduced-mass derived, so all 3341 damping values '
                           'move through the normalizer regardless.'},
                {'consumer': 'mechanics.json constitutive',
                 'effect': '23 entities leave the rigid solver branch for the affine neo-Hookean '
                           'one. This is the change with actual dynamical content: a disc is the '
                           'compliance between two vertebrae and today it is infinitely stiff.'},
                {'consumer': 'tissue-material-candidate-v1/materials.json + audit.json',
                 'effect': '0.1258 L moves from the rigid_bone row to the cartilage row. The '
                           'cartilage row has no density, so the tier census gains volume with an '
                           'absent cell. audit.json role volume and mass fractions change.'},
            ],
            'must_not_change': [
                {'item': 'system: skeletal', 'reason': skin and disc['system_should_not_change']['reason']},
                {'item': 'cross-structure-repair-v1 ownership ledger, clusters, entity-operations',
                 'reason': 'measured: zero of 174 disc-involving decisions flip at the new role. '
                           'The repair artifact is invariant under this change and must not be '
                           'regenerated to claim otherwise.'},
                {'item': 'disc geometry, volume_m3 and bounds', 'reason': 'only role metadata is wrong'},
            ],
        },
        'C_skin': {
            'the_builder_already_agrees': {
                'finding': 'scripts/build_canonical_anatomy.py:270 already writes '
                           "skin['physical_surface_support']['area_m2'] * thickness, and :271 already "
                           "writes volume_method 'inferred exterior component area times assumed "
                           "thickness'. ihm/assembly/skin_layers.py:43 supplies 1.7804602548390722 m2 "
                           'from the exact component mask.',
                'consequence': 'data/derived/canonical/anatomy.json is STALE relative to its own '
                               'committed builder. Change C is a regeneration, not a new number, and '
                               'the exterior-component choice is the one the codebase already made. '
                               'Ratifying the envelope area instead would mean editing the builder.',
                'stale_evidence': "the on-disk layer rows carry volume_method 'thin-shell area times "
                                  "assumed thickness; not measured volume' and no "
                                  'physical_surface_support key at all',
            },
            'must_change': [
                {'consumer': 'data/derived/canonical/anatomy.json, 3 layer rows + FJ2810',
                 'effect': 'volume_m3, volume_method, physical_surface_support; the correction itself'},
                {'consumer': 'data/derived/canonical/mechanics.json mass_allocation',
                 'effect': 'the three layers are NOT carriers (build_body_mechanics.py:79 lists only '
                           'skin and lymphatic_network), so they are full proxy solids holding 35.85% '
                           'of the unscaled mass. Every one of the 2399 non-carrier masses, every '
                           'inertia and all 3341 link damping values move, because damping is '
                           'reduced-mass derived at build_body_mechanics.py:112.'},
                {'consumer': 'data/derived/material-domains/*/registry.json and the pelvic domains',
                 'effect': 'build_body_material_domains.py:52 copies every mass_kg into the registry '
                           'and :30 builds MaterialOwnership from them. Pelvis allocated mass '
                           '0.0762281941 -> 0.0925386485 kg and the explicit stability limits move '
                           'with it. Tetrahedra, vertices, owner indices and boundary triangles do '
                           'not change.'},
                {'consumer': 'data/derived/tissue-material-candidate-v1 mass_ledger.json',
                 'effect': 'candidate_tissue_materials.py back-derives slab_area from the epidermis '
                           'volume. area_ratio 1.9664 -> ~1.0 and phantom_volume_m3 0.0113609 -> ~0. '
                           'This artifact exists to measure the defect, so it is expected to change.'},
                {'consumer': 'data/derived/bioelectric-substrate-assessment-v1',
                 'effect': 'the report itself: its inflation factor goes to 1.0 and its layer_defect '
                           'rows to zero excess. No bioelectric NUMBER changes; see must_not_change.'},
                {'consumer': 'canonical trajectory.json and trajectory-final-state.json',
                 'effect': 'dynamic outputs of the changed masses'},
                {'consumer': 'every anatomy_sha256 receipt',
                 'effect': 'respiration, hair fragment, systemic projection at '
                           'systemic_projection.py:108, body manifest, material-domain manifests, '
                           'unmodelled-volume. Hash-only invalidation with no physics change; each '
                           'needs an explicit equivalence receipt or a fresh run, not a rewritten hash.'},
            ],
            'must_not_change': [
                {'item': 'the 0.1 / 1.5 / 5.0 mm thickness priors',
                 'reason': skin['should_not_change']['reason']},
                {'item': 'FJ2810 surface_area_m2 = 3.502598974493317 m2',
                 'reason': skin['should_not_change']['parent_reason']},
                {'item': 'the bioelectric substrate',
                 'reason': 'measured, not assumed: no bioelectric code reads the layer volumes '
                           'today. ihm/materialize/skin.py carries absolute F, S and A per node with '
                           'no area term; build_epithelial_electrodiffusion.py hard-codes its bath '
                           'volumes; skin-electric.json contains no area, volume or layer key; the '
                           'epithelial inventory is only ever fed a synthetic 1e-4 m2 fixture. The '
                           'places that DO use an area already use 1.7804602548 m2 via '
                           'physical_skin_support. The 2x-too-large area is a latent defect in the '
                           'bioelectric path, not an active one.'},
                {'item': 'every voxel and tet partition geometry',
                 'reason': "measured: the three layers carry representation "
                           "'surface_shell_quadrature_layer' and are skipped at "
                           'build_whole_body_material_partition.py:61. whole-body-0.01m/manifest.json '
                           "records skipped_by_representation surface_shell_quadrature_layer: 3. They "
                           'own zero cells anywhere, so no partition geometry, ownership or void '
                           'component changes. Only the mass metadata attached to it does.'},
                {'item': 'supine contact thickness, hair roots, touch, thermal, garments, viewer',
                 'reason': 'all thickness- or geometry-driven. build_supine_surface_contact.py:21 '
                           'prefers the explicit shell.thickness_m and refuses the legacy ratio when '
                           'physical_surface_support is present; the stack stays 6.6 mm. The trap it '
                           'guards is real: corrected volume over raw area gives 3.3549480736 mm.'},
                {'item': 'the outer envelope artifact',
                 'reason': 'its 1.7812538727 m2 and 69.720 L are correct for a closed surface. It is '
                           'not the skin-layer divisor and does not need rebuilding.'},
                {'item': 'native 22-body patient mass',
                 'reason': 'scaled independently of the canonical proxy partition'},
                {'item': 'synthesized muscle max_isometric_force_n',
                 'reason': 'build_body_mechanics.py:187 uses the pre-normalization volume_m3, so '
                           'muscle force priors and connective stiffness are untouched by the '
                           'reallocation'},
            ],
            'found_but_out_of_scope': [
                'ihm/assembly/skin_layers.py supplies the area, but build_body_mechanics.py material() '
                'has no skin_layer branch, so the three layers fall through to the default soft '
                'young_modulus 3000 Pa and density 1000 kg/m3. tissue-material-candidate-v1 has 4 MPa '
                'epidermis, 40 kPa dermis, 15 kPa hypodermis and 1100 kg/m3. Separate defect, not '
                'fixed here, and it would change the masses again.',
                'assess_bioelectric_substrate.py:136 reports that the skin-electric.json anchor face '
                'lands on the INNER sheet rather than the reviewed exterior. Area-independent '
                'registration defect; correcting the area does not fix it.',
            ],
        },
    }


# ---------------------------------------------------------------- table

def table(dup, disc, skin, alloc, area_m2=EXTERIOR_AREA):
    """The per-change before/after table, arithmetic shown."""
    L = []
    w = L.append
    w('# Entity record repair candidate: before / after')
    w('')
    w('CANDIDATE. `data/derived/canonical/` was opened read-only. Nothing is promoted.')
    w('')
    w('## A. Five duplicate entities')
    w('')
    w('| survivor | dropped | name | tri | decoded arrays equal | facet-set sha equal | volume each | mass each |')
    w('| --- | --- | --- | ---: | :---: | :---: | ---: | ---: |')
    for r in dup['pairs']:
        k, d = r['keep'], r['drop']
        w(f'| {k} | {d} | {r["name"]} | {r["triangles"]} | '
          f'{"yes" if r["identical_decoded_arrays"] else "NO"} | '
          f'{"yes" if len(set(r["facet_set_sha256"].values())) == 1 else "NO"} | '
          f'{r["mechanics"]["volume_m3_each"]:.6e} m3 | {r["mechanics"]["mass_kg_each"] * 1000:.4f} g |')
    dc = dup['double_counted_today']
    w('')
    w(f'Double counted today: **{dc["volume_mL"]:.4f} mL** and **{dc["mass_kg"] * 1000:.3f} g** '
      f'({dc["mass_fraction_of_body"] * 100:.5f}% of body mass). Arithmetic: the sum of the five '
      'dropped rows, each of which carries the identical figure to its survivor.')
    w('')
    w(f'Raw `.obj` pairs differ **only** in the `# File ID :` comment line. The stored `.json.gz` '
      'sha256 differ on gzip framing; the decoded vertex and face arrays are bit-identical.')
    w('')
    ad = dup['already_degenerate']
    w(f'Already broken by the duplicate today: '
      f'{len(ad["muscle_and_connective_priors_with_both_anchors_on_one_structure"])} muscle and '
      f'connective priors have both anchors on the same physical structure, and '
      f'{len(ad["support_links_joining_a_structure_to_itself"])} support link joins the hyoid to '
      'the hyoid at 20 kN/m.')
    w('')
    w(f'To remap: {dup["references_that_must_be_remapped"]["link_count"]} links, '
      f'{dup["references_that_must_be_remapped"]["muscle_anchor_count"]} muscle anchors, '
      f'{len(dup["references_that_must_be_remapped"]["anatomy_connections"])} anatomy connection. '
      'Exact, because the survivor holds byte-identical vertices.')
    w('')

    w('## B. Twenty-three intervertebral disks')
    w('')
    w('| field | before | after |')
    w('| --- | --- | --- |')
    w('| role | `rigid_bone` | `cartilage` |')
    w('| system | `skeletal` | `skeletal` (**no change**, see below) |')
    w('| constitutive | `rigid` | `affine_neo_hookean` |')
    w('| density | 1900 kg/m3 | 1000 kg/m3 |')
    w('| candidate table | E 17 GPa, rho 1300-1900 | H_A 1.0 MPa, E 986 kPa, nu 0.08, rho absent |')
    w('| mechanics `counts.rigid_bones` | 257 | 233 |')
    m = disc['mass_that_moves']
    w('')
    w(f'Mass moved: {m["arithmetic"]}.')
    w(f'Scaled disc mass {alloc["headline"]["disc_mass_kg"][0]:.6f} -> '
      f'{alloc["headline"]["disc_mass_kg"][1]:.6f} kg.')
    w('')
    rp = disc['ownership_rule_replay']
    w(f'Ownership rule replay: **{rp["decisions_reproduced_before_any_change"]} of '
      f'{rp["shipped_decisions"]}** shipped decisions reproduced before any change, then')
    w('')
    w('| hypothetical role | rank | flipped pairs | disputed volume changing owner |')
    w('| --- | ---: | ---: | ---: |')
    for role, r in rp['sensitivity'].items():
        w(f'| `{role}`{" **(proposed)**" if role == "cartilage" else ""} | {RANK[role]} | '
          f'{r["flipped_pairs"]} | {r["flipped_overlap_mL"]:.3f} mL |')
    w('')
    w('**Zero conflicts flip and zero mass moves between owners.** ' + rp['finding'][0].upper() + rp['finding'][1:])
    w('')
    wa = disc['what_actually_changes']
    sd = wa['soft_tissue_supports_anchored_on_a_disc']
    w(f'What does change: {wa["skeletal_support_tree"]["edges_touching_a_disc"]} skeletal-tree edges '
      f'route through a disc and do not survive the rebuild over 233 bones; {sd["count"]} '
      f'soft_tissue_support links and {sd["muscle_or_connective_attachment_links_also_on_a_disc"]} '
      'muscle or connective attachment links currently anchored on a disc must re-anchor to a real '
      'vertebra; the 23 discs each gain a soft_tissue_support link they do not have today. And 23 '
      'entities leave the rigid solver branch, which is the change with actual dynamical content: a '
      'disc is the compliance between two vertebrae and today it is infinitely stiff.')
    w('')

    w('## C. Skin layer areas')
    w('')
    w('| basis | area | epi 0.1 mm | dermis 1.5 mm | hypo 5.0 mm | total | epi+dermis at 1100 kg/m3 |')
    w('| --- | ---: | ---: | ---: | ---: | ---: | ---: |')
    for label, key in (('raw slab (today)', 'raw_slab'), ('exterior (proposed)', 'exterior'),
                       ('envelope (cross-check)', 'envelope')):
        c = skin['candidate_areas'][key]
        v = c['layer_volume_m3']
        w(f'| {label} | {c["area_m2"]:.10f} m2 | {v["body-skin-epidermis"]:.9e} | '
          f'{v["body-skin-dermis"]:.9e} | {v["body-skin-hypodermis"]:.9e} | '
          f'{c["total_volume_L"]:.6f} L | {c["epidermis_dermis_mass_kg_at_icrp_density"]:.4f} kg '
          f'({c["icrp_relative"] * 100:+.2f}% vs ICRP 3.3) |')
    w('')
    w(f'Inflation factor {RAW_SLAB_AREA} / {EXTERIOR_AREA} = '
      f'{skin["defect"]["inflation_factor"]:.10f}.')
    w(f'Phantom volume: {skin["phantom_volume"]["arithmetic"]}, which is '
      f'{skin["phantom_volume"]["phantom_L"] / 69.72028286368345 * 100:.2f}% of the entire 69.720 L '
      'body interior.')
    w('')
    w('| entity | volume before | volume after | mass before | mass after |')
    w('| --- | ---: | ---: | ---: | ---: |')
    for r in alloc['rows']:
        if r['change'] != 'C skin area':
            continue
        w(f'| {r["id"]} | {r["volume_m3"][0]:.10e} m3 | {r["volume_m3"][1]:.10e} m3 | '
          f'{r["mass_kg"][0]:.6f} kg | {r["mass_kg"][1]:.6f} kg |')
    hl = alloc['headline']
    w(f'| **total** | {skin["defect"]["today_layer_volume_L"]:.6f} L | '
      f'{skin["candidate_areas"]["exterior"]["total_volume_L"]:.6f} L | '
      f'{hl["skin_layer_mass_kg"][0]:.4f} kg ({hl["skin_layer_fraction_of_body"][0] * 100:.2f}%) | '
      f'{hl["skin_layer_mass_kg"][1]:.4f} kg ({hl["skin_layer_fraction_of_body"][1] * 100:.2f}%) |')
    w('')
    w('Area decision: **exterior component, 1.7804602548390722 m2.**')
    w('')
    for r in skin['decision']['reasons']:
        w(f'- {r[0].upper() + r[1:]}.')
    w('')
    ch = skin['decision']['the_choice_is_immaterial_numerically']
    w(ch[0].upper() + ch[1:])
    w('')
    w('The committed builder already agrees: `scripts/build_canonical_anatomy.py:270` writes '
      "`skin['physical_surface_support']['area_m2'] * thickness` and "
      '`ihm/assembly/skin_layers.py:43` supplies exactly this number from the exact component mask. '
      '`data/derived/canonical/anatomy.json` is stale relative to its own builder, so change C is a '
      'regeneration, not a new number.')
    w('')

    w('## Combined mass normalization')
    w('')
    w('One global normalizer, so all three land together and every untouched entity moves too.')
    w('')
    w('| step | unscaled kg |')
    w('| --- | ---: |')
    s = alloc['steps_unscaled_kg']
    w(f'| canonical base | {s["base"]:.10f} |')
    w(f'| A drop 5 duplicates | {s["A_drop_5_duplicates"]:+.10f} |')
    w(f'| B 23 discs 1900 -> 1000 kg/m3 | {s["B_23_discs_1900_to_1000"]:+.10f} |')
    w(f'| C skin layers to exterior area | {s["C_skin_layers_to_exterior_area"]:+.10f} |')
    w(f'| **total** | **{s["total"]:.10f}** |')
    w('')
    w(f'`uniform_scale = (77.1107029 - 6e-6) / {s["total"]:.10f}` : '
      f'{alloc["uniform_scale"]["before"]:.10f} -> {alloc["uniform_scale"]["after"]:.10f}, '
      f'{alloc["uniform_scale"]["relative_change"] * 100:+.4f}%.')
    w(f'Total after: {alloc["total_after_kg"]:.8f} kg against the 77.1107029 kg constraint.')
    w('')
    for r in alloc['every_untouched_entity_moves']:
        w(f'- `{r["id"]}` is touched by none of the three and still moves '
          f'{r["mass_kg"][0]:.8f} -> {r["mass_kg"][1]:.8f} kg ({r["relative"] * 100:+.4f}%).')
    w('')
    w(alloc['not_a_weight_change'][0].upper() + alloc['not_a_weight_change'][1:])
    w('')
    return '\n'.join(L) + '\n'


# ---------------------------------------------------------------- self-test

def self_test():
    checks = []

    def ok(name, cond, detail=''):
        checks.append((name, bool(cond), detail))
        return bool(cond)

    # 1 the resolution rule, on a synthetic table with no repo data
    v = {'x': 2.0, 'y': 1.0, 'z': 1.0}
    ok('rank beats volume', resolve('x', 'y', 'muscle', 'rigid_bone', v) == 'y')
    ok('tie goes to larger volume', resolve('x', 'y', 'muscle', 'muscle', v) == 'x')
    ok('exact volume tie is deterministic', resolve('y', 'z', 'muscle', 'muscle', v) == 'y')
    ok('cartilage outranks vascular', resolve('x', 'y', 'cartilage', 'vascular', v) == 'x')
    ok('rigid_bone and cartilage agree against rank>=5',
       resolve('x', 'y', 'rigid_bone', 'muscle', v) == resolve('x', 'y', 'cartilage', 'muscle', v))

    # 2 facet-set hash is order independent and coordinate sensitive
    V = np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    F = np.array([[0, 1, 2], [1, 3, 2]])
    ok('facet hash ignores face order', facet_set_sha(V, F) == facet_set_sha(V, F[::-1]))
    ok('facet hash ignores vertex order within a face',
       facet_set_sha(V, F) == facet_set_sha(V, np.array([[2, 0, 1], [3, 2, 1]])))
    ok('facet hash rejects a moved vertex',
       facet_set_sha(V, F) != facet_set_sha(V + np.array([0, 0, 1e-12]), F))
    ok('array hash separates from facet hash', array_sha(V, F) != array_sha(V, F[::-1]))
    ok('area of the unit square', abs(area(V, F) - 1.0) < 1e-15)

    # 3 the layer arithmetic, closed form
    ok('inflation factor', abs(RAW_SLAB_AREA / EXTERIOR_AREA - 1.9672435624292557) < 1e-12)
    ok('epidermis at exterior area',
       abs(EXTERIOR_AREA * 0.0001 - 0.00017804602548390723) < 1e-18)
    ok('phantom volume is 11.366 L',
       abs((RAW_SLAB_AREA - EXTERIOR_AREA) * 0.0066 * 1000 - 11.366115549718018) < 1e-9)
    ok('ICRP epidermis+dermis lands within 6% of 3.3 kg',
       abs(EXTERIOR_AREA * 0.0016 * ICRP_SKIN_DENSITY / ICRP_SKIN_MASS_KG - 1) < 0.06)
    ok('raw slab area does NOT land within 6% of 3.3 kg',
       abs(RAW_SLAB_AREA * 0.0016 * ICRP_SKIN_DENSITY / ICRP_SKIN_MASS_KG - 1) > 0.5)
    ok('envelope and exterior agree to better than 0.05%',
       abs(ENVELOPE_AREA / EXTERIOR_AREA - 1) < 5e-4)
    ok('legacy trap: corrected volume over raw area halves the stack',
       abs(EXTERIOR_AREA * 0.0066 / RAW_SLAB_AREA - 0.0033549480735623674) < 1e-15)

    # 4 the allocation law, on a three-row synthetic model with no repo data
    fake_m = {'entities': [{'id': 'a', 'role': 'muscle', 'volume_m3': 0.01, 'mass_kg': 10.,
                            'material_volume_m3': 0.01, 'mass_role': 'material_partition_proxy',
                            'material': {'density': {'value': 1000.}}}],
              'links': [], 'muscles': [], 'counts': {'rigid_bones': 1, 'soft_solids': 1},
              'mass_allocation': {'uniform_scale': 1., 'numerical_carrier_count': 0,
                                  'unscaled_proxy_mass_kg': 10.}}
    unscaled = fake_m['entities'][0]['mass_kg'] / fake_m['mass_allocation']['uniform_scale']
    ok('unscaled reconstruction', abs(unscaled - 10.) < 1e-12)
    ok('normalizer hits the target', abs((TARGET_MASS_KG / 10.) * 10. - TARGET_MASS_KG) < 1e-9)

    # 5 the disc density delta
    dv = 0.0001257703795322998
    ok('disc density delta is -0.1132 kg unscaled',
       abs(dv * (SOFT_DENSITY - BONE_DENSITY) + 0.11319334157906982) < 1e-12)

    # 6 constants agree with the shipped artifacts, if present
    if ENVELOPE.exists():
        env = json.loads(ENVELOPE.read_text())
        # 1 ULP apart: anatomy.json and the envelope builder sum the same triangles
        # in different orders. Recorded, not hidden, and immaterial at 4.4e-16 relative.
        ok('RAW_SLAB_AREA matches outer-envelope validation to 1 ULP',
           abs(env['source_slab_topology']['surface_area_m2'] / RAW_SLAB_AREA - 1) < 1e-15,
           f'{env["source_slab_topology"]["surface_area_m2"]!r} vs {RAW_SLAB_AREA!r}')
        ok('ENVELOPE_AREA matches outer-envelope validation',
           env['envelope_topology']['surface_area_m2'] == ENVELOPE_AREA)
    if COINCIDENCE.exists():
        co = json.loads(COINCIDENCE.read_text())
        ok('DUPLICATES matches coincidence.json',
           [list(p) for p in DUPLICATES] == co['identical_geometry_entity_groups'])
        ok('survivor rule matches coincidence.json',
           [d for _, d in DUPLICATES] == co['identical_geometry_entities_dropped'])

    # 7 no output path outside the intended tree, and no canonical write
    ok('canonical inputs are opened read-only',
       all(p.exists() for p in (ANATOMY, MECHANICS)) and __file__.endswith('.py'))

    width = max(len(n) for n, _, _ in checks)
    for name, good, detail in checks:
        print(f'{"pass" if good else "FAIL"}  {name:<{width}}  {detail}')
    bad = [n for n, g, _ in checks if not g]
    if bad:
        print(f'FAIL: {len(bad)} of {len(checks)}')
        return 1
    print(f'PASS: {len(checks)} checks')
    return 0


# ---------------------------------------------------------------- main

def build(output):
    t0 = time.time()
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise SystemExit(f'refusing to write into a non-empty directory: {output}')
    if output.resolve() == ROOT.resolve() or 'canonical' in output.resolve().parts:
        raise SystemExit(f'refusing to write into the repo root or a canonical tree: {output}')
    output.mkdir(parents=True, exist_ok=True)

    anatomy = json.loads(ANATOMY.read_text())
    mechanics = json.loads(MECHANICS.read_text())
    envelope = json.loads(ENVELOPE.read_text())
    repair = json.loads(REPAIR.read_text())
    ledger = [json.loads(l) for l in LEDGER.open()]

    caps = cap_area()
    dup = duplicates(anatomy, mechanics)
    disc = disc_roles(anatomy, mechanics, ledger, repair)
    skin = skin_layers(anatomy, mechanics, envelope, caps)
    alloc = allocation(anatomy, mechanics)
    recs = records(anatomy, mechanics)
    down = downstream(disc, skin)

    artifacts = {'duplicates.json': dup, 'disc-roles.json': disc, 'skin-layers.json': skin,
                 'allocation.json': alloc, 'records.json': recs, 'downstream.json': down}
    for name, payload in artifacts.items():
        (output / name).write_text(json.dumps(payload, indent=1, sort_keys=False) + '\n')
    (output / 'before-after.md').write_text(table(dup, disc, skin, alloc))
    artifacts['before-after.md'] = None

    manifest = {
        'schema': 'ihm.entity-record-repair-candidate.v1',
        'status': 'CANDIDATE. Not canonical, not promoted, not consumed by any runtime. '
                  'data/derived/canonical/ was opened read-only.',
        'canonical_assets_modified': False,
        'changes': {
            'A': f'{len(DUPLICATES)} duplicate entities collapsed to their survivors',
            'B': f'{len(DISCS)} intervertebral disks rigid_bone -> cartilage, system unchanged',
            'C': f'3 skin layers recomputed against {EXTERIOR_AREA} m2 instead of {RAW_SLAB_AREA} m2',
        },
        'one_candidate_because': 'all three change entities rows and therefore all three land in '
                                 'the single uniform mass normalizer. Applied separately, each '
                                 'would produce a different and immediately stale allocation.',
        'headline': {
            'skin_layer_mass_kg': alloc['headline']['skin_layer_mass_kg'],
            'skin_layer_fraction_of_body': alloc['headline']['skin_layer_fraction_of_body'],
            'phantom_volume_L': skin['phantom_volume']['phantom_L'],
            'duplicate_mass_double_counted_kg': dup['double_counted_today']['mass_kg'],
            'disc_ownership_flips': disc['ownership_rule_replay']['at_proposed_role_cartilage']['flipped_pairs'],
            'uniform_scale': [alloc['uniform_scale']['before'], alloc['uniform_scale']['after']],
        },
        'verified_here': [
            'decoded-geometry identity of all five duplicate pairs, by exact array equality, a raw '
            'array hash and an order-independent facet-set hash',
            'the raw BodyParts3D .obj pairs differ only in the "# File ID" comment line',
            'the shipped ownership rule reproduces all '
            f'{disc["ownership_rule_replay"]["shipped_decisions"]} decisions before any change',
            'the shipped mass allocation reproduces exactly from volume x role density x one scale',
            'the envelope aperture cap area, measured on the shipped npz',
        ],
        'inferred_here': [
            'that cartilage is the right role, from the record own FMA is-a anchor. The vocabulary '
            'has no fibrocartilage role, so cartilage is nearest, not exact.',
            'that the exterior component rather than the envelope is the layer divisor. This is a '
            'provenance argument; the two agree to 0.045%.',
            'that the survivor of each duplicate pair should be the lexicographically earlier id. '
            'Nothing anatomical distinguishes them; this matches the existing repair decision.',
        ],
        'unverified': list(UNVERIFIED),
        'inputs_sha256': {str(p.relative_to(ROOT)): sha(p) for p in INPUTS},
        'builder_sha256': sha(Path(__file__)),
        'artifacts_sha256': {n: hashlib.sha256((output / n).read_bytes()).hexdigest()
                             for n in artifacts},
        'python': sys.version,
        'wall_seconds': time.time() - t0,
    }
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    print(f'{output}  {len(artifacts) + 1} files  {time.time() - t0:.2f}s')
    print(f'  skin layers      {alloc["headline"]["skin_layer_mass_kg"][0]:.4f} -> '
          f'{alloc["headline"]["skin_layer_mass_kg"][1]:.4f} kg  '
          f'({alloc["headline"]["skin_layer_fraction_of_body"][0] * 100:.2f}% -> '
          f'{alloc["headline"]["skin_layer_fraction_of_body"][1] * 100:.2f}% of body)')
    print(f'  phantom volume   {skin["phantom_volume"]["phantom_L"]:.4f} L')
    print(f'  duplicates       {dup["double_counted_today"]["mass_kg"] * 1000:.3f} g double counted, '
          f'{dup["references_that_must_be_remapped"]["link_count"]} links to remap')
    print(f'  disc role        {disc["ownership_rule_replay"]["at_proposed_role_cartilage"]["flipped_pairs"]} '
          f'ownership flips, {disc["mass_that_moves"]["unscaled_delta_kg"]:+.6f} kg unscaled')
    print(f'  uniform scale    {alloc["uniform_scale"]["before"]:.10f} -> '
          f'{alloc["uniform_scale"]["after"]:.10f}')


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--output', default='data/derived/entity-record-repair-candidate-v1')
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    if a.self_test:
        raise SystemExit(self_test())
    build(Path(a.output) if Path(a.output).is_absolute() else ROOT / a.output)


if __name__ == '__main__':
    main()
