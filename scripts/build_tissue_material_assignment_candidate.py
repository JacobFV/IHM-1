#!/usr/bin/env python3
"""Measure what the 1,597 promoted display structures did to the mass ledger, and
what a sourced material assignment does to it.

    .venv/bin/python scripts/build_tissue_material_assignment_candidate.py --self-test
    .venv/bin/python scripts/build_tissue_material_assignment_candidate.py \
        --output data/derived/tissue-material-assignment-candidate-v1

Nothing here is canonical. Every canonical asset is opened read-only. The
assignment itself lives in `ihm/assembly/tissue_materials.py`, which
`scripts/build_body_mechanics.py` imports, so a rebuild reaches the same answer;
this script measures that rule against the shipped one and writes the receipts.

Measured here (arithmetic over shipped assets, all of it re-derived by importing
the shipped rule rather than restating it):

  decomposition   the unscaled proxy excess over the declared body mass, split
                  into the part that is wrong DENSITY and the part that is
                  wrong VOLUME, with the volume part broken down by which
                  representation defect produced it
  assignment      one row per entity: tissue class, density, density source and
                  tier, volume rule, thickness and thickness tier
  sensitivity     the assignment re-run over the full sweep bounds of every
                  assumed thickness, so no closure rests on an untested constant
  colocation      a name-key pass over the 291 promoted structures that a
                  geometry-only key would have called aliases, saying which
                  are resolvable from naming alone and which are not
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.assembly import tissue_materials as tm

ANATOMY = ROOT / 'data/derived/canonical/anatomy.json'
MECHANICS = ROOT / 'data/derived/canonical/mechanics.json'
PROFILE = ROOT / 'data/derived/canonical/profile.json'
COLOCATION = ROOT / 'data/derived/display-structure-promotion-v1/surface-colocation.json'
LUMEN = ROOT / 'data/derived/interstitial-matrix-v1/lumen.json'
COMPOSITION = ROOT / 'data/derived/interstitial-matrix-v1/composition.json'
INPUTS = (ANATOMY, MECHANICS, PROFILE, COLOCATION, LUMEN, COMPOSITION)

# The densities the shipped builder used before this assignment. Not a tissue
# density in either case: 1000 is water and 1900 is marrow-free cortical bone.
SUPERSEDED_DENSITY = {'rigid_bone': 1900., 'other': 1000.}
# The allocation those densities produced, frozen here rather than read back from
# mechanics.json. A receipt that recomputes its own "before" from an artifact the
# promotion overwrites stops being a receipt the moment the promotion runs; this
# script therefore states the superseded numbers and reports whether the
# canonical asset still carries them.
SUPERSEDED_ALLOCATION = {'unscaled_proxy_mass_kg': 93.47748928171038,
                         'uniform_scale': 0.7570923710490256,
                         'entity_volume_m3': 0.07473798928171045,
                         'measured_on': 'the 4000-entity canonical model at anatomy.json sha256 '
                                        '874479e434e2637c08a1f5153de88fdaf138d33f55ae91a3bb7b1578ed891d5f, '
                                        'immediately after the 1,597-structure display promotion'}
PROMOTION_ASSUMPTION = 'DISPLAY-STRUCTURE-PROMOTION'


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def carrier(entity):
    """The shipped builder's carrier rule, restated so this script sees it."""
    return (entity['role'] in ('skin', 'lymphatic_network', 'surface_region')
            or 'cavity of' in entity['name'].lower())


# ---------------------------------------------------------------- assignment

def assign(entities, prior_volume, thickness=None):
    """One row per entity under the shipped rule, at the given thickness priors."""
    rows = []
    for e in entities:
        name, role = e['name'], e['role']
        volume, representation = tm.material_volume(name, role, prior_volume[e['id']],
                                                    e.get('surface_area_m2'),
                                                    bool(e.get('watertight_edge_incidence')),
                                                    thickness)
        rho, record = tm.density(name, role)
        is_carrier = carrier(e)
        rows.append({
            'entity_id': e['id'], 'name': name, 'role': role, 'system': e.get('system'),
            'evidence_kind': e.get('evidence_kind'),
            'promoted_display_structure': PROMOTION_ASSUMPTION in (e.get('assumptions') or []),
            'watertight_edge_incidence': bool(e.get('watertight_edge_incidence')),
            'surface_area_m2': e.get('surface_area_m2'),
            'superseded_volume_m3': prior_volume[e['id']],
            'volume_m3': 0. if is_carrier else volume,
            'volume_rule': 'numerical_boundary_carrier' if is_carrier else representation['rule'],
            'volume_thickness_m': representation.get('thickness_m'),
            'volume_thickness_tier': representation.get('thickness_tier'),
            'volume_basis': representation['basis'],
            'density_kg_m3': record['value_kg_m3'], 'tissue_class': record['tissue_class'],
            'density_source': record['source'], 'density_tier': record['tier'],
            'density_note': record['note'],
            'superseded_density_kg_m3': (SUPERSEDED_DENSITY['rigid_bone'] if role == 'rigid_bone'
                                         else SUPERSEDED_DENSITY['other']),
            'mass_kg': 0. if is_carrier else volume * rho,
        })
    return rows


# ---------------------------------------------------------------- decomposition

def decompose(entities, prior_volume, rows, target_mass_kg):
    """Split the unscaled proxy excess into density, volume and residual.

    The split is done by evaluating the SAME proxy sum four times and taking
    differences, so no term is estimated: superseded density on superseded
    volume, sourced density on superseded volume, sourced density on assigned
    volume, and the target.
    """
    by_id = {e['id']: e for e in entities}
    old_d_old_v = new_d_old_v = new_d_new_v = 0.
    for r in rows:
        if r['volume_rule'] == 'numerical_boundary_carrier':
            continue
        old_d_old_v += r['superseded_volume_m3'] * r['superseded_density_kg_m3']
        new_d_old_v += r['superseded_volume_m3'] * r['density_kg_m3']
        new_d_new_v += r['volume_m3'] * r['density_kg_m3']
    excess = old_d_old_v - target_mass_kg
    density_term = new_d_old_v - old_d_old_v
    volume_term = new_d_new_v - new_d_old_v

    # Where the volume term came from, by representation rule, priced at the
    # sourced density so the parts sum to the whole.
    by_rule = {}
    for r in rows:
        if r['volume_rule'] == 'numerical_boundary_carrier':
            continue
        released = (r['superseded_volume_m3'] - r['volume_m3'])
        if abs(released) < 1e-15:
            continue
        slot = by_rule.setdefault(r['volume_rule'], {'entities': 0, 'volume_released_m3': 0.,
                                                     'mass_released_kg': 0., 'examples': []})
        slot['entities'] += 1
        slot['volume_released_m3'] += released
        slot['mass_released_kg'] += released * r['density_kg_m3']
        if len(slot['examples']) < 8:
            slot['examples'].append({'name': r['name'], 'role': r['role'],
                                     'superseded_volume_m3': r['superseded_volume_m3'],
                                     'volume_m3': r['volume_m3'],
                                     'surface_area_m2': r['surface_area_m2']})
    for slot in by_rule.values():
        slot['examples'].sort(key=lambda x: -x['superseded_volume_m3'])

    # The colocation term is measured, not modelled: the promoted structures a
    # geometry-only key would have called aliases, priced under the assignment.
    coloc = json.loads(COLOCATION.read_text())
    coloc_ids = {r['promoted_id'] for r in coloc['best_match_per_structure']}
    row_by_id = {r['entity_id']: r for r in rows}
    coloc_mass_superseded = sum(row_by_id[i]['superseded_volume_m3']
                                * row_by_id[i]['superseded_density_kg_m3']
                                for i in coloc_ids if i in row_by_id)
    coloc_mass_assigned = sum(row_by_id[i]['mass_kg'] for i in coloc_ids if i in row_by_id)

    promoted = [r for r in rows if r['promoted_display_structure']]
    return {
        'question': 'the unscaled proxy mass is %.4f kg against a declared %.4f kg body. How much '
                    'of the %.4f kg excess is wrong density and how much is wrong volume?'
                    % (old_d_old_v, target_mass_kg, excess),
        'superseded_density_on_superseded_volume_kg': old_d_old_v,
        'recomputed_superseded_matches_the_frozen_receipt': abs(
            old_d_old_v - SUPERSEDED_ALLOCATION['unscaled_proxy_mass_kg']) < 1e-6,
        'declared_body_mass_kg': target_mass_kg,
        'unscaled_proxy_excess_kg': excess,
        'attributable_to_density_kg': density_term,
        'attributable_to_volume_kg': volume_term,
        'remaining_after_both_kg': new_d_new_v - target_mass_kg,
        'reading': 'moving every entity from the superseded 1000/1900 kg/m3 pair to its sourced '
                   'tissue density changes the proxy by %.4f kg, which is %.1f%% of the excess. '
                   'Correcting the volume each surface actually holds as tissue changes it by '
                   '%.4f kg, which is %.1f%%. The excess is a VOLUME defect, not a density one, '
                   'and a density fudge large enough to absorb it would have to take every soft '
                   'tissue below adipose.'
                   % (density_term, 100. * abs(density_term) / abs(excess),
                      volume_term, 100. * abs(volume_term) / abs(excess)),
        'volume_term_by_representation_rule': by_rule,
        'volume_double_count_by_colocation': {
            'promoted_structures_clearing_loose_alias_bounds': len(coloc_ids),
            'source': 'data/derived/display-structure-promotion-v1/surface-colocation.json',
            'unscaled_proxy_mass_at_superseded_density_kg': coloc_mass_superseded,
            'mass_under_this_assignment_kg': coloc_mass_assigned,
            'statement': 'this mass is NOT part of the volume term above and is not removed by '
                         'this assignment. It is the mass of promoted structures that a geometry-'
                         'only key would have called aliases of a structure the model already '
                         'carried. Where such a pair is one structure under two names the model '
                         'carries it twice, and no representation rule can see that: it needs a '
                         'concept key. See colocation.json for what naming alone can settle.'},
        'promoted_structures': {
            'count': len(promoted),
            'superseded_proxy_mass_kg': sum(r['superseded_volume_m3'] * r['superseded_density_kg_m3']
                                            for r in promoted),
            'assigned_mass_kg': sum(r['mass_kg'] for r in promoted),
        },
        'not_addressed_here': [
            'part-of nesting: `mucosa of stomach` (0.425 L) is a layer of `stomach` (0.567 L) and '
            'both are entities. The BodyParts3D part-of graph would settle these; this assignment '
            'does not read it.',
            '69.4 L of the superseded 74.7 L proxy volume rested on an OPEN-surface signed '
            'integral. This assignment only overrides that prior where the surface topology says '
            'the enclosed volume is categorically not the entity (a sheet, a lumen); for a viscus '
            'with a small hole the open integral is still the estimate.',
        ],
    }


# ---------------------------------------------------------------- sensitivity

def sensitivity(entities, prior_volume, target_mass_kg):
    """Re-run the assignment at every corner of every assumed thickness prior."""
    out = {}
    base = {k: v[0] for k, v in tm.THICKNESS.items()}
    for key, (value, tier, bounds, note) in tm.THICKNESS.items():
        if tier != 'assumed':
            continue
        row = {'central_m': value, 'tier': tier, 'sweep_m': list(bounds), 'note': note, 'cases': {}}
        for label, t in (('lower', bounds[0]), ('central', value), ('upper', bounds[1])):
            override = {k: (t if k == key else base[k], tm.THICKNESS[k][1], tm.THICKNESS[k][2],
                            tm.THICKNESS[k][3]) for k in tm.THICKNESS}
            rows = assign(entities, prior_volume, override)
            proxy = sum(r['mass_kg'] for r in rows)
            row['cases'][label] = {'thickness_m': t, 'unscaled_proxy_mass_kg': proxy,
                                   'uniform_scale': target_mass_kg / proxy,
                                   'entity_volume_m3': sum(r['volume_m3'] for r in rows)}
        span = (row['cases']['upper']['unscaled_proxy_mass_kg']
                - row['cases']['lower']['unscaled_proxy_mass_kg'])
        row['proxy_mass_span_over_sweep_kg'] = span
        out[key] = row
    return out


# ---------------------------------------------------------------- colocation

ORDINAL = {'first': 1, 'second': 2, 'third': 3, 'fourth': 4, 'fifth': 5, 'sixth': 6, 'seventh': 7,
           'eighth': 8, 'ninth': 9, 'tenth': 10, 'eleventh': 11, 'twelfth': 12}
SPINE = {'cervical': 'c', 'thoracic': 't', 'lumbar': 'l', 'sacral': 's'}
# Token substitutions that standard anatomical nomenclature settles without a
# concept key. Each is an accepted alternative NAME for one structure in
# Terminologia Anatomica or in standard dental and vascular usage; none is a
# judgement about two neighbouring structures. Listed one per line so the
# reasoning is auditable rather than buried in a regex.
SYNONYM_TOKENS = [
    ('flexor accessorius', 'quadratus plantae'),   # TA: m. quadratus plantae, syn. flexor accessorius
    ('colli', 'cervicis'),                         # Latin genitive of "neck"; same muscle group
    ('little finger', 'fifth finger'),             # digitus minimus manus
    ('ring finger', 'fourth finger'),
    ('middle finger', 'third finger'),
    ('index finger', 'second finger'),
    ('little toe', 'fifth toe'),
    ('central incisor', 'medial incisor'),         # dens incisivus medialis
    ('hemi azygos', 'hemiazygos'),
]
# Words that are number, laterality or grouping conventions rather than identity.
DROP_WORDS = {'secondary', 'muscle', 'bone', 'tooth', 'set', 'group', 'of', 'the'}
LATERALITY = {'left', 'right'}
# The only extra words a subset match may differ by. Laterality, because an atlas
# may or may not lateralise a name; hand and foot, because a name may or may not
# say which limb a digit belongs to and two structures that differ only by hand
# versus foot cannot be geometrically colocated. Every other extra word is an
# anatomical qualifier -- `perforating`, `mucosa`, `branch`, `segment` -- and
# distinguishes a structure from its parent or its neighbour, so a subset match
# on one of those would collapse a part into its whole.
CONVENTION_EXTRAS = {'left', 'right', 'hand', 'foot'}


def _singular(word):
    if len(word) > 3 and word.endswith('ies'):
        return word[:-3] + 'y'
    if len(word) > 3 and word.endswith(('ches', 'shes', 'sses', 'xes')):
        return word[:-2]
    if len(word) > 3 and word.endswith('s') and not word.endswith('ss'):
        return word[:-1]
    return word


def normal_name(name):
    """A comparison key that absorbs the naming conventions the two atlases differ on."""
    n = name.lower().strip().replace('-', '')
    for canonical, alternative in SYNONYM_TOKENS:
        n = n.replace(alternative, canonical)
    n = re.sub(r'\bdisk\b', 'disc', n)
    for word, number in ORDINAL.items():         # 'eleventh thoracic vertebra' -> 't11'
        for latin, letter in SPINE.items():
            n = n.replace('%s %s vertebra' % (word, latin), 'vertebra %s%d' % (letter, number))
            n = n.replace('intervertebral disc of %s %s vertebra' % (word, latin),
                          'intervertebral disc %s%d' % (letter, number))
    # the vertebra substitution above rewrites the disc's own qualifier too
    n = re.sub(r'\bdisc of vertebra\b', 'disc', n)
    # An intervertebral disc is named either for the vertebra above it or for the
    # pair it lies between; `t4t5` and `t4` are the same disc.
    n = re.sub(r'\b([ctls])(\d+)(\1)(\d+)\b',
               lambda m: '%s%s' % (m.group(1), m.group(2))
               if int(m.group(4)) == int(m.group(2)) + 1 else m.group(0), n.replace('-', ''))
    # 'navicular bone of left foot' and 'left navicular bone' are one convention
    # apart, not two structures.
    n = re.sub(r'\bof (left|right) (foot|hand|leg|arm|thigh|forearm|shoulder)\b', r'\1', n)
    n = re.sub(r'[^a-z0-9]+', ' ', n).strip()
    words = [_singular(w) for w in n.split()]
    return ' '.join(sorted(w for w in words if w not in DROP_WORDS))


def _word_set(key):
    return frozenset(key.split())


def resolve_colocations():
    """What a name key alone can settle among the 291, and what it cannot."""
    report = json.loads(COLOCATION.read_text())
    resolved, unresolved = [], []
    for r in report['best_match_per_structure']:
        a, b = normal_name(r['promoted_name']), normal_name(r['canonical_name'])
        sa, sb = _word_set(a), _word_set(b)
        why = None
        if a == b:
            why = 'names normalise to the same key under the synonym, ordinal/numeral, singular, '\
                  '"set of"/"head of"/"part of" and secondary-dentition conventions'
        elif (sa and sb and (sa ^ sb) and (sa ^ sb) <= CONVENTION_EXTRAS
              and r['agreement']['bbox_iou'] >= 0.5):
            why = 'the two keys differ only by %s, which is a naming convention and not an '\
                  'anatomical qualifier: one name is the unlateralised or unregionalised form of '\
                  'the other, and the bounding boxes agree at IoU %.3f' \
                  % (' and '.join(sorted(sa ^ sb)), r['agreement']['bbox_iou'])
        row = {'promoted_id': r['promoted_id'], 'promoted_name': r['promoted_name'],
               'promoted_role': r['promoted_role'], 'canonical_id': r['canonical_id'],
               'canonical_name': r['canonical_name'], 'canonical_role': r['canonical_role'],
               'bbox_iou': r['agreement']['bbox_iou'],
               'clears_strict_duplicate_bounds': r['clears_duplicate_bounds'],
               'normalised_promoted': a, 'normalised_canonical': b}
        if why:
            row['verdict'] = 'same_structure_under_two_names'
            row['grounds'] = why
            resolved.append(row)
        else:
            row['verdict'] = 'unresolved'
            row['grounds'] = 'naming does not settle it and the repository holds no concept key: '\
                             'the Z-Anatomy display records carry no FMA or UBERON identifier, and '\
                             'no external synonym source was obtainable in the session that wrote '\
                             'this script. Geometry cannot settle it either -- the same bounds '\
                             'contain true cross-naming duplicates and true distinct neighbours.'
            unresolved.append(row)
    return {
        'schema': 'ihm.promoted-colocation-name-resolution.v1',
        'input': 'data/derived/display-structure-promotion-v1/surface-colocation.json',
        'question': 'of the promoted structures a geometry-only key would have called aliases, '
                    'which does NAMING settle?',
        'considered': len(resolved) + len(unresolved),
        'resolved_as_duplicate': len(resolved),
        'unresolved': len(unresolved),
        'action_taken': 'NONE. This is a report. Collapsing a promoted structure changes the '
                        'entity set and the reference geometry, which is a different lane from a '
                        'material assignment; the pairs are recorded so the double count is '
                        'auditable rather than invisible, and so the lane that does collapse them '
                        'has a defended starting list.',
        'resolved': sorted(resolved, key=lambda r: -r['bbox_iou']),
        'unresolved_pairs': sorted(unresolved, key=lambda r: -r['bbox_iou']),
    }


# ---------------------------------------------------------------- fat budget

def fat_budget(rows, composition):
    """Where body fat lives once the hypodermis shell carries an adipose density.

    The interstitial composition build can only account for fat in the FILL -- the
    void outside every entity. The declared fraction is whole-body, so the
    remainder has to be inside entity geometry, and the build's own ICRP 89
    cross-check says how much of that is defensible. Until this assignment there
    was no adipose-density material anywhere in the entity set, so none of it
    could be located.
    """
    closure = composition['body_fat']['closure']
    cross = closure['icrp89_cross_check']
    fat_of_adipose = 0.8
    adipose_rows = [r for r in rows if r['tissue_class'] == 'adipose']
    entity_adipose = sum(r['mass_kg'] for r in adipose_rows)
    entity_fat = entity_adipose * fat_of_adipose
    required = closure['required_total_fat_kg']
    fill = closure['fat_carried_by_the_fill_kg']
    non_separable = cross['scaled_non_separable_fat_kg']
    return {
        'declared_fraction': composition['body_fat']['declared_fraction'],
        'required_total_fat_kg': required,
        'fat_carried_by_the_interstitial_fill_kg': fill,
        'before': {
            'adipose_density_entities': 0,
            'fat_locatable_in_entity_geometry_kg': 0.,
            'unlocated_fat_kg': required - fill,
            'icrp89_scaled_non_separable_fat_kg': non_separable,
            'excess_over_icrp_non_separable_kg': cross['residual_minus_scaled_non_separable_fat_kg'],
        },
        'after': {
            'adipose_density_entities': len(adipose_rows),
            'adipose_entities': [{'name': r['name'], 'volume_m3': r['volume_m3'],
                                  'density_kg_m3': r['density_kg_m3'], 'mass_kg': r['mass_kg']}
                                 for r in adipose_rows],
            'entity_adipose_tissue_kg': entity_adipose,
            'fat_locatable_in_entity_geometry_kg': entity_fat,
            'unlocated_fat_kg': required - fill - entity_fat,
            'icrp89_scaled_non_separable_fat_kg': non_separable,
            'excess_over_icrp_non_separable_kg': required - fill - entity_fat - non_separable,
        },
        'fat_fraction_of_adipose_tissue': fat_of_adipose,
        'double_count_disclosure': {
            'statement': 'the hypodermis shell has no volumetric geometry, is skipped by the '
                         'occupancy pass, and describes space INSIDE the void the interstitial '
                         'matrix fills. Crediting its adipose here and keeping the fill\'s '
                         'subcutaneous adipose counts part of the same space twice.',
            'shell_volume_m3': sum(r['volume_m3'] for r in adipose_rows),
            'fill_subcutaneous_compartment_m3': sum(
                c['volume_m3'] for c in composition['constituents']
                if c['compartment'] == 'subcutaneous'),
            'decision': 'the assignment gives the shell the adipose DENSITY, which is what it is '
                        'made of, and does not change its volume or add it to the fill. The '
                        'overlap is a property of the interstitial lane\'s ledger, not of this '
                        'material table, and is reported there and here rather than netted out '
                        'silently in either.',
        },
    }


# ---------------------------------------------------------------- cross-check

# ICRP 89 reference adult male, Tables 2.8 and 2.9, scaled to this body by mass.
# Not a validation: the reference male is 73 kg and 176 cm and is not this
# specimen. It is the only independent statement of what a body of this mass is
# made of that the repository holds, so it is the check the assignment has to
# survive.
ICRP89_REFERENCE = {
    'total_body_kg': 73.0,
    'skeletal_muscle_kg': 29.0,
    'total_skeleton_kg': 10.5,
    'skin_kg': 3.3,
    'lung_with_blood_kg': 1.2,
    'adipose_tissue_kg': 18.2,
}


def icrp_cross_check(rows, target_mass_kg, uniform_scale):
    """Weigh the assignment against the only independent composition the repo holds."""
    factor = target_mass_kg / ICRP89_REFERENCE['total_body_kg']

    def mass(predicate):
        return sum(r['mass_kg'] for r in rows if predicate(r)) * uniform_scale

    checks = {
        'skeletal_muscle': (mass(lambda r: r['tissue_class'] == 'skeletal_muscle'),
                            ICRP89_REFERENCE['skeletal_muscle_kg'] * factor,
                            'ICRP 89 skeletal muscle'),
        'skeleton': (mass(lambda r: r['tissue_class'] == 'skeleton_with_marrow'),
                     ICRP89_REFERENCE['total_skeleton_kg'] * factor,
                     'ICRP 89 total skeleton; the atlas geometry is whole bones including the '
                     'marrow cavity, which is what the 1300 kg/m3 row describes'),
        'skin_epidermis_and_dermis': (mass(lambda r: r['tissue_class'] == 'skin'),
                                      ICRP89_REFERENCE['skin_kg'] * factor,
                                      'ICRP 89 skin. The hypodermis is deliberately excluded: ICRP '
                                      'counts it as adipose tissue, not as skin, and so does this '
                                      'assignment'),
        'lung': (mass(lambda r: r['tissue_class'] == 'lung'),
                 ICRP89_REFERENCE['lung_with_blood_kg'] * factor,
                 'ICRP 89 lung with blood'),
        'adipose_in_entity_geometry': (mass(lambda r: r['tissue_class'] == 'adipose'),
                                       ICRP89_REFERENCE['adipose_tissue_kg'] * factor,
                                       'ICRP 89 TOTAL adipose tissue. The entity set holds only the '
                                       'hypodermis shell, so this row is expected to fall short by '
                                       'whatever the interstitial fill carries; it is here to say '
                                       'by how much'),
    }
    return {
        'reference': dict(ICRP89_REFERENCE, source='icrp89', tier='transferred',
                          note='reference adult male, 73 kg, 176 cm. Not this specimen.'),
        'scale_to_this_body': factor,
        'rows': {k: {'assigned_kg': a, 'icrp89_scaled_kg': b, 'ratio': a / b if b else None,
                     'difference_kg': a - b, 'basis': note}
                 for k, (a, b, note) in checks.items()},
        'reading': 'the assigned masses are AFTER the uniform normalizer, which is the number the '
                   'model actually simulates. Nothing here was tuned to these targets: the '
                   'densities come from ICRU-44 and ICRP 89 density tables and the volumes come '
                   'from this specimen\'s own surfaces, so the agreement is a check, not a fit.',
    }


# ---------------------------------------------------------------- provenance

# How a density tier and a volume rule map onto the four tiers of
# ihm.structure-provenance.v1. `measured` is reserved for a value measured on
# THIS specimen: the surface geometry is, so a closed-surface volume reaches it;
# no density does.
VOLUME_TIER = {
    'solid': ('measured', 'the signed divergence integral of this specimen\'s own closed surface'),
    'closed_membrane_measured': ('measured', 'the signed divergence integral of this specimen\'s '
                                             'own closed surface'),
    'hollow_viscus_thinner_than_its_wall': ('measured', 'the enclosed volume of this specimen\'s '
                                                        'own surface, already below one wall '
                                                        'thickness so no lumen was subtracted'),
    'open_membrane_sheet': ('synthesized', 'this specimen\'s measured surface AREA times an assumed '
                                           'membrane thickness that no retrieved source constrains'),
    'hollow_viscus_wall': ('synthesized', 'this specimen\'s measured surface AREA times an assumed '
                                          'wall thickness that no retrieved source constrains'),
    'vessel_wall': ('synthesized', 'this specimen\'s measured surface AREA times the assumed 0.3 mm '
                                   'wall the builder has always used'),
    'skin_boundary_carrier': ('synthesized', 'a 1 micrometre numerical carrier, not a tissue volume'),
    'numerical_boundary_carrier': ('synthesized', 'zero material volume by the carrier rule; this '
                                                  'structure is a boundary, not a tissue'),
}
DENSITY_TIER = {'transferred': 'transferred', 'assumed': 'synthesized'}
WEAKEST = ('measured', 'transferred', 'derived', 'synthesized')


def provenance_records(rows, module_sha, builder_sha, inputs):
    """One ihm.structure-provenance.v1 record per entity, for its MATERIAL.

    The canonical structure-provenance index answers where a structure's geometry
    came from. This answers the same question for the two numbers that turn that
    geometry into mass, and it is emitted for every entity, not only the ones
    whose value changed, so no assignment is traceable only by absence.
    """
    for r in rows:
        volume_tier, volume_basis = VOLUME_TIER[r['volume_rule']]
        density_tier = DENSITY_TIER[r['density_tier']]
        tier = max(volume_tier, density_tier, key=WEAKEST.index)
        source = tm.SOURCES[r['density_source']]
        yield {
            'schema': 'ihm.structure-provenance.v1',
            'record_kind': 'material_assignment',
            'structure_id': r['entity_id'],
            'model_id': 'ihm-body',
            'canonical_entity_id': r['entity_id'],
            'name': r['name'], 'system': r['system'], 'role': r['role'],
            'evidence_kind': r['evidence_kind'],
            'dataset': {
                'id': r['density_source'],
                'label': source['finding'].split(':')[0],
                'version': None, 'revision': None,
                'url': source['url'],
                'specimen': 'reference-man tabulation. Not this specimen, and not any single '
                            'measured subject.',
                'units': 'kg/m3', 'frame': None,
                'attribution': source['finding'],
                'acquisition_status': 'published reference table read into this repository',
                'license': None,
                'license_absent_reason': 'a numeric density from a published reference table; the '
                                         'table itself is cited, no file of it is retained here',
            },
            'source_file': {
                'path': 'ihm/assembly/tissue_materials.py', 'sha256': module_sha,
                'sha256_verified': True, 'bytes': None,
                'note': 'the shipped module that owns the density table and the representation '
                        'rules. The value below is read from it, not restated.',
            },
            'build': {
                'script': 'scripts/build_body_mechanics.py',
                'script_sha256': sha(ROOT / 'scripts/build_body_mechanics.py'),
                'candidate_script': 'scripts/build_tissue_material_assignment_candidate.py',
                'candidate_script_sha256': builder_sha,
                'commit': None, 'commit_covers_working_tree': False, 'uncommitted_changes': True,
                'inputs_sha256': inputs,
            },
            'assignment': {
                'tissue_class': r['tissue_class'],
                'density_kg_m3': r['density_kg_m3'],
                'density_source': r['density_source'],
                'density_tier': r['density_tier'],
                'density_note': r['density_note'],
                'superseded_density_kg_m3': r['superseded_density_kg_m3'],
                'volume_m3': r['volume_m3'],
                'volume_rule': r['volume_rule'],
                'volume_basis': r['volume_basis'],
                'volume_thickness_m': r['volume_thickness_m'],
                'volume_thickness_tier': r['volume_thickness_tier'],
                'superseded_volume_m3': r['superseded_volume_m3'],
                'surface_area_m2': r['surface_area_m2'],
                'watertight_edge_incidence': r['watertight_edge_incidence'],
                'unscaled_mass_kg': r['mass_kg'],
            },
            'transforms': [
                {'step': 'tissue_class_from_name_and_role', 'residual': None,
                 'residual_reason': 'a classification, not a fit'},
                {'step': 'density_lookup', 'residual': None,
                 'residual_reason': 'a table read; no fit to this specimen was performed, and none '
                                    'is possible without weighing it'},
                {'step': 'volume_rule:' + r['volume_rule'], 'residual': None,
                 'residual_reason': volume_basis},
            ],
            'frame_relation': 'canonical',
            'tier': tier,
            'tier_basis': 'the weaker of the volume tier (%s) and the density tier (%s). A density '
                          'from a reference table is transferred, never measured; a value no '
                          'retrieved source constrains is synthesized.' % (volume_tier, density_tier),
            'tier_evidence': {'volume': volume_basis, 'density': r['density_note']},
            'assumptions': ['TISSUE-MATERIAL-ASSIGNMENT'],
            'completeness': {
                'answerable': True,
                'measured_on_this_specimen': {'geometry': True, 'density': False,
                                              'thickness': False},
            },
        }


# ---------------------------------------------------------------- closure

def mass_closure(rows, composition, target_mass_kg):
    """Does the assigned tissue plus the space it leaves weigh what the body weighs?

    The entity set does not tile the body: the interstitial matrix build measured
    a void that belongs to no entity. So a correct entity proxy must come out
    BELOW the declared mass, and the interesting question is whether the density
    the remainder then needs is a density any tissue has. Under the superseded
    numbers the remainder was NEGATIVE 22.7 kg, which is not a body.
    """
    shells = [r for r in rows if r['role'] == 'skin_layer']
    shell_volume = sum(r['volume_m3'] for r in shells)
    proxy = sum(r['mass_kg'] for r in rows)
    volume = sum(r['volume_m3'] for r in rows)
    envelope = composition['entity_side']['volume_m3'] + composition['fill']['volume_m3']
    # The three skin shells have no volumetric geometry, so they are not part of
    # the voxel occupancy the envelope is partitioned by; the space they describe
    # is inside the void. They are therefore excluded from the volume that is
    # subtracted from the envelope, and named as the overlap they are.
    solid_volume = volume - shell_volume
    remainder_mass = target_mass_kg - proxy
    # Two readings, because the shells are ambiguous by construction. The void
    # is what the occupancy pass leaves, and the shells are skipped by it, so
    # the shells sit INSIDE that void. Either the shells are the model's
    # description of the outer part of that void, in which case the void free of
    # them is what the remainder mass has to fill; or they are treated as
    # separate, in which case they overlap it. Both are stated.
    void = envelope - solid_volume
    void_free_of_shells = void - shell_volume
    return {
        'envelope_interior_m3': envelope,
        'envelope_basis': 'data/derived/interstitial-matrix-v1/composition.json entity_side plus '
                          'fill, the 1 mm partition of the acquired envelope',
        'assigned_entity_volume_m3': volume,
        'assigned_entity_volume_excluding_skin_shells_m3': solid_volume,
        'skin_shell_volume_m3': shell_volume,
        'assigned_entity_mass_kg': proxy,
        'declared_body_mass_kg': target_mass_kg,
        'remainder_mass_kg': remainder_mass,
        'void_outside_every_occupying_entity_m3': void,
        'void_free_of_the_skin_shells_m3': void_free_of_shells,
        'remainder_density_if_shells_are_part_of_the_void_kg_m3': (
            remainder_mass / void_free_of_shells if void_free_of_shells > 0 else None),
        'remainder_density_if_shells_are_counted_separately_kg_m3': (
            remainder_mass / void if void > 0 else None),
        'reading': 'the assigned entities weigh %.4f kg. The declared body needs a further '
                   '%.4f kg from space that belongs to no occupying entity. The three skin '
                   'quadrature shells have no volumetric geometry and are skipped by the occupancy '
                   'pass, so their %.6f m3 lies inside that space and the answer depends on how '
                   'they are read: taking the shells AS the outer part of the void leaves '
                   '%.6f m3 needing %.1f kg/m3, which is within a percent of ICRU-44 adipose '
                   'tissue at 950 and below loose connective tissue and interstitial fluid at '
                   '1030; counting them separately leaves %.6f m3 '
                   'needing %.1f kg/m3, which no tissue reaches and which measures the overlap '
                   'rather than a material. Under the superseded densities and volumes the same '
                   'arithmetic asked the remainder to weigh %.4f kg, which is negative mass in a '
                   'positive volume and is not a body.'
                   % (proxy, remainder_mass, shell_volume, void_free_of_shells,
                      remainder_mass / void_free_of_shells if void_free_of_shells > 0 else float('nan'),
                      void, remainder_mass / void if void > 0 else float('nan'),
                      target_mass_kg - SUPERSEDED_ALLOCATION['unscaled_proxy_mass_kg']),
        'unresolved': 'this assignment does not decide which reading is right. Deciding it is the '
                      'act of materialising the interstitial matrix as a first-class structure, '
                      'which retires the shells; that is a separate lane and it is the lane the '
                      'double_count_hazard block of the interstitial composition build is written '
                      'for.',
    }


# ---------------------------------------------------------------- build

def build(output):
    anatomy = json.loads(ANATOMY.read_text())
    mech = json.loads(MECHANICS.read_text())
    profile = json.loads(PROFILE.read_text())
    composition = json.loads(COMPOSITION.read_text())
    entities = anatomy['entities']
    # The superseded volume is read from the geometric fallback the builder
    # produced BEFORE any representation rule: mechanics.json carries the volume
    # after the rule, so the pre-rule number is recovered from the record the
    # rule writes, and is the mechanics volume where no rule fired.
    prior_volume = {}
    for s in mech['entities']:
        rep = s.get('volume_representation') or {}
        prior_volume[s['id']] = float(rep.get('open_surface_signed_integral_m3')
                                      or rep.get('enclosed_volume_including_lumen_m3')
                                      or s['volume_m3'])
    for e in entities:
        prior_volume.setdefault(e['id'], e.get('volume_m3') or 0.)
    target = float(profile['mass_kg'])

    rows = assign(entities, prior_volume)
    proxy = sum(r['mass_kg'] for r in rows)
    doc = {
        'schema': 'ihm.tissue-material-assignment-candidate.v1',
        'status': 'CANDIDATE. Every canonical asset is opened read-only. The assignment itself is '
                  'shipped in ihm/assembly/tissue_materials.py and is applied to the canonical '
                  'model by scripts/promote_tissue_material_assignment.py.',
        'rule_owner': 'ihm/assembly/tissue_materials.py',
        'rule_owner_sha256': sha(ROOT / 'ihm/assembly/tissue_materials.py'),
        'builder': 'scripts/build_tissue_material_assignment_candidate.py',
        'builder_sha256': sha(Path(__file__)),
        'python': sys.version,
        'inputs_sha256': {str(p.relative_to(ROOT)): sha(p) for p in INPUTS},
        'tier_vocabulary': {
            'measured': 'measured on THIS specimen. No density or thickness in this table reaches '
                        'it; the geometry does, and is marked so.',
            'transferred': 'a published value or reference-table row carried in from a donor '
                           'cohort. A density taken from ICRU-44 or ICRP 89 is transferred.',
            'derived': 'computed from other rows in this table by a recorded step.',
            'assumed': 'no retrieved source constrains it. Every assumed value is swept.',
        },
        'sources': tm.SOURCES,
        'density_table': {k: {'value_kg_m3': v[0], 'source': v[1], 'tier': v[2], 'note': v[3]}
                          for k, v in tm.DENSITY.items()},
        'thickness_table': {k: {'value_m': v[0], 'tier': v[1], 'sweep_m': list(v[2]), 'note': v[3]}
                            for k, v in tm.THICKNESS.items()},
        'allocation': {
            'declared_body_mass_kg': target,
            'superseded': dict(SUPERSEDED_ALLOCATION,
                               canonical_mechanics_still_carries_it=abs(
                                   mech['mass_allocation']['unscaled_proxy_mass_kg']
                                   - SUPERSEDED_ALLOCATION['unscaled_proxy_mass_kg']) < 1e-6,
                               canonical_mechanics_live_value_kg=mech['mass_allocation']['unscaled_proxy_mass_kg']),
            'assigned': {'unscaled_proxy_mass_kg': proxy, 'uniform_scale': target / proxy,
                         'entity_volume_m3': sum(r['volume_m3'] for r in rows)},
        },
        'icrp89_cross_check': icrp_cross_check(rows, target, target / proxy),
        'mass_closure': mass_closure(rows, composition, target),
        'decomposition': decompose(entities, prior_volume, rows, target),
        'fat_budget': fat_budget(rows, composition),
        'sensitivity': sensitivity(entities, prior_volume, target),
        'counts_by_volume_rule': _census(rows, 'volume_rule'),
        'counts_by_tissue_class': _census(rows, 'tissue_class'),
        'counts_by_density_tier': _census(rows, 'density_tier'),
        'registered_entities': _registered(rows),
        'not_established': [
            'no density or thickness here was measured on this specimen. Every density is a '
            'reference-table row and every thickness is an assumed prior.',
            'the membrane thickness and the hollow-viscus wall thickness are `assumed`: no source '
            'in this repository constrains either, and no external measurement was obtainable in '
            'the session that wrote this table. Both are swept and the span is reported.',
            'the assignment reduces volume only where surface topology says the enclosed volume is '
            'categorically not the entity. It does not resolve part-of nesting, and it does not '
            'resolve the 291 colocated promoted structures, which need a concept key.',
        ],
    }
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / 'assignment.json').write_text(json.dumps(doc, indent=1) + '\n')
    (output / 'entities.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in rows))
    (output / 'colocation.json').write_text(json.dumps(resolve_colocations(), indent=1) + '\n')
    records = list(provenance_records(rows, doc['rule_owner_sha256'], doc['builder_sha256'],
                                      doc['inputs_sha256']))
    (output / 'material-provenance.jsonl').write_text(
        ''.join(json.dumps(r) + '\n' for r in records))
    tiers = {}
    for r in records:
        tiers[r['tier']] = tiers.get(r['tier'], 0) + 1
    (output / 'provenance-index.json').write_text(json.dumps({
        'schema': 'ihm.tissue-material-provenance-index.v1',
        'record_schema': 'ihm.structure-provenance.v1',
        'record_kind': 'material_assignment',
        'records': len(records),
        'records_path': 'material-provenance.jsonl',
        'records_by_tier': dict(sorted(tiers.items(), key=lambda kv: -kv[1])),
        'statement': 'one record per canonical entity, emitted whether or not the entity\'s value '
                     'changed, so no assignment is traceable only by absence. `measured` is never '
                     'reached by a density: it is reserved for a value measured on this specimen, '
                     'and only the surface geometry qualifies.',
        'ids': [r['structure_id'] for r in records],
    }, indent=1) + '\n')
    doc['provenance'] = {'records': len(records), 'records_by_tier': tiers,
                         'path': 'material-provenance.jsonl'}
    (output / 'assignment.json').write_text(json.dumps(doc, indent=1) + '\n')
    manifest = {
        'schema': 'ihm.tissue-material-assignment-candidate-manifest.v1',
        'builder': 'scripts/build_tissue_material_assignment_candidate.py',
        'builder_sha256': sha(Path(__file__)),
        'rule_owner_sha256': sha(ROOT / 'ihm/assembly/tissue_materials.py'),
        'inputs_sha256': doc['inputs_sha256'],
        'canonical_assets_modified': False,
        'artifacts_sha256': {p.name: sha(p) for p in sorted(output.iterdir())
                             if p.is_file() and p.name != 'manifest.json'},
    }
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    print(json.dumps({'allocation': doc['allocation'],
                      'decomposition': {k: v for k, v in doc['decomposition'].items()
                                        if not isinstance(v, dict)},
                      'fat_budget_after': doc['fat_budget']['after']['excess_over_icrp_non_separable_kg'],
                      'counts_by_volume_rule': doc['counts_by_volume_rule']}, indent=1))
    return doc


def _census(rows, key):
    out = {}
    for r in rows:
        slot = out.setdefault(r[key], {'entities': 0, 'volume_m3': 0., 'mass_kg': 0.})
        slot['entities'] += 1
        slot['volume_m3'] += r['volume_m3']
        slot['mass_kg'] += r['mass_kg']
    return dict(sorted(out.items(), key=lambda kv: -kv[1]['mass_kg']))


def _registered(rows):
    """The 1,767 entities whose geometry was fitted in from another dataset."""
    reg = [r for r in rows if r['evidence_kind'] == 'registered_geometry']
    return {'entities': len(reg),
            'by_role': _census(reg, 'role'),
            'by_tissue_class': _census(reg, 'tissue_class'),
            'by_volume_rule': _census(reg, 'volume_rule'),
            'by_density_tier': _census(reg, 'density_tier')}


# ---------------------------------------------------------------- self-test

class Tests(unittest.TestCase):
    def test_open_sheet_loses_its_enclosure_and_a_closed_one_does_not(self):
        # An open fascia sheet of 0.17 m2 wrapping 4 L keeps 0.17 L, not 4 L.
        v, rep = tm.material_volume('left fascia lata', 'connective_tissue', 0.004086, 0.16932, False)
        self.assertEqual(rep['rule'], 'open_membrane_sheet')
        self.assertAlmostEqual(v, 0.00016932, places=9)
        # The same sheet authored closed measures itself and is left alone.
        v, rep = tm.material_volume('left piriformis fascia', 'connective_tissue', 0.000034, 0.0092, True)
        self.assertEqual(rep['rule'], 'closed_membrane_measured')
        self.assertAlmostEqual(v, 0.000034, places=9)

    def test_a_muscle_named_after_a_fascia_is_not_a_membrane(self):
        v, rep = tm.material_volume('left tensor fasciae latae', 'muscle', 0.0000493, 0.0116, False)
        self.assertEqual(rep['rule'], 'solid')
        self.assertAlmostEqual(v, 0.0000493, places=9)

    def test_hollow_viscus_drops_its_lumen_but_never_below_its_wall(self):
        v, rep = tm.material_volume('stomach', 'soft_organ', 0.0005665, 0.0383, False)
        self.assertEqual(rep['rule'], 'hollow_viscus_wall')
        self.assertAlmostEqual(v, 0.0383 * 0.003, places=9)
        self.assertGreater(rep['lumen_volume_m3'], 0.)
        thin, rep = tm.material_volume('ureter', 'soft_organ', 1e-9, 0.0005, False)
        self.assertEqual(rep['rule'], 'hollow_viscus_thinner_than_its_wall')
        self.assertEqual(thin, 1e-9)

    def test_hypodermis_is_the_only_adipose_material_and_it_is_950(self):
        rho, record = tm.density('hypodermis', 'skin_layer')
        self.assertEqual(rho, 950.)
        self.assertEqual(record['tissue_class'], 'adipose')
        self.assertEqual(tm.density('dermis', 'skin_layer')[0], 1100.)

    def test_bone_leaves_the_marrow_free_cortical_density(self):
        rho, record = tm.density('left femur', 'rigid_bone')
        self.assertEqual(rho, 1300.)
        self.assertEqual(record['tier'], 'transferred')
        self.assertNotEqual(rho, SUPERSEDED_DENSITY['rigid_bone'])

    def test_every_assumed_row_is_swept_and_every_transferred_row_names_a_source(self):
        for key, (value, source, tier, note) in tm.DENSITY.items():
            self.assertIn(tier, ('transferred', 'assumed'), key)
            self.assertIn(source, tm.SOURCES, key)
            self.assertTrue(note.strip(), key)
        for key, (value, tier, bounds, note) in tm.THICKNESS.items():
            if tier == 'assumed':
                self.assertLess(bounds[0], bounds[1], key)
                self.assertTrue(bounds[0] <= value <= bounds[1], key)

    def test_decomposition_terms_sum_to_the_measured_excess(self):
        anatomy = json.loads(ANATOMY.read_text())
        mech = json.loads(MECHANICS.read_text())
        prior = {}
        for s in mech['entities']:
            rep = s.get('volume_representation') or {}
            prior[s['id']] = float(rep.get('open_surface_signed_integral_m3')
                                   or rep.get('enclosed_volume_including_lumen_m3')
                                   or s['volume_m3'])
        target = float(json.loads(PROFILE.read_text())['mass_kg'])
        rows = assign(anatomy['entities'], prior)
        d = decompose(anatomy['entities'], prior, rows, target)
        self.assertAlmostEqual(d['unscaled_proxy_excess_kg'] + d['attributable_to_density_kg']
                               + d['attributable_to_volume_kg'], d['remaining_after_both_kg'],
                               places=6)

    def test_name_key_separates_a_cross_naming_duplicate_from_a_true_neighbour(self):
        self.assertEqual(normal_name('vertebra t11'), normal_name('eleventh thoracic vertebra'))
        self.assertEqual(normal_name('left navicular bone'), normal_name('navicular bone of left foot'))
        self.assertEqual(normal_name('right upper first molar tooth'),
                         normal_name('right upper first secondary molar tooth'))
        self.assertNotEqual(normal_name('left internal abdominal oblique muscle'),
                            normal_name('left external oblique'))
        self.assertNotEqual(normal_name('posterior lateral segment of liver vii'),
                            normal_name('hepatovenous segment ix'))
        self.assertEqual(normal_name('hemi-azygos vein'), normal_name('hemiazygos vein'))
        self.assertEqual(normal_name('intervertebral disc t4-t5'),
                         normal_name('intervertebral disk of fourth thoracic vertebra'))

    def test_a_part_is_never_collapsed_into_its_whole(self):
        # A subset match may only differ by a naming convention. `perforating
        # branches of X` and `X` are a branch and its parent, and `mucosa of
        # stomach` is a layer of the stomach; neither may resolve as a duplicate.
        for promoted, canonical in (
                ('right perforating branches of plantar metatarsal artery', 'plantar metatarsal artery'),
                ('mucosa of stomach', 'stomach'),
                ('left atrium', 'wall of left atrium'),
                ('left transversus abdominis muscle', 'left external oblique')):
            a, b = _word_set(normal_name(promoted)), _word_set(normal_name(canonical))
            self.assertTrue(a != b and not (a ^ b) <= CONVENTION_EXTRAS, (promoted, canonical))

    def test_a_naming_convention_is_collapsed(self):
        for promoted, canonical in (
                ('left external intercostal muscles', 'external intercostal muscle'),
                ('right dorsal metacarpal arteries', 'set of dorsal metacarpal arteries'),
                ('right middle phalanx of fifth finger of hand', 'middle phalanx of right little finger')):
            a, b = _word_set(normal_name(promoted)), _word_set(normal_name(canonical))
            self.assertTrue(a == b or (a ^ b) <= CONVENTION_EXTRAS, (promoted, canonical))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', default='data/derived/tissue-material-assignment-candidate-v1')
    p.add_argument('--self-test', action='store_true')
    a = p.parse_args()
    if a.self_test:
        result = unittest.TextTestRunner(verbosity=2).run(
            unittest.TestLoader().loadTestsFromTestCase(Tests))
        raise SystemExit(0 if result.wasSuccessful() else 1)
    build(ROOT / a.output if not Path(a.output).is_absolute() else Path(a.output))


if __name__ == '__main__':
    main()
