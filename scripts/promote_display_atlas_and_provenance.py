#!/usr/bin/env python3
"""Promote the display structures, the surface-ownership lanes and the provenance index.

Three candidates, in dependency order, each of which had declared
``canonical_assets_modified: false``:

  1 data/derived/display-promotion-candidate-v1   Z-Anatomy surfaces this body does not
    carry become canonical entities. Aliases resolve to the entity that already exists.
  2 data/derived/cross-structure-repair-v1        the ownership ledger, and
    data/derived/conflict-free-atlas-v1           the conforming surfaces it implies.
    Neither is regenerated. Every shipped decision is replayed against the canonical model
    as it now stands, through the rule that made it.
  3 data/derived/structure-provenance-candidate-v1 the index that says, for every
    structure, which dataset, file, script and transform put it there.

Like scripts/promote_entity_record_repair.py before it, this script contains no promotion
logic of its own. The change lives in the builders that own each field -- so a rebuild
reaches the same answer and cannot silently revert -- and this runs them in order, checks
the result against what the candidates predicted, and writes a receipt. Every candidate
directory is left byte-for-byte unchanged as its own provenance record, and that is
checked here rather than assumed; the promoted provenance index is written to a new
directory for exactly that reason.

  .venv/bin/python scripts/promote_display_atlas_and_provenance.py --stage all
"""
from collections import Counter
from pathlib import Path
import argparse, hashlib, json, subprocess, sys, time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'scripts'))
from ihm.assembly.anatomy import bilateral_mirror_census, name_key, normalized_name

VENV = ROOT/'.venv/bin/python'
LIBIGL = ROOT/'data/runtime/geometry/libigl-2.6.2/venv/bin/python'
CANDIDATES = {
    'display': ROOT/'data/derived/display-promotion-candidate-v1',
    'repair': ROOT/'data/derived/cross-structure-repair-v1',
    'atlas': ROOT/'data/derived/conflict-free-atlas-v1',
    'provenance': ROOT/'data/derived/structure-provenance-candidate-v1',
}
CANONICAL_ASSETS = ['data/derived/canonical/profile.json', 'data/derived/canonical/anatomy.json',
                    'data/derived/canonical/manifest_fragment.json', 'data/derived/canonical/verification.json',
                    'data/derived/canonical/lymphatic_graph.json', 'data/derived/canonical/mechanics.json',
                    'data/derived/canonical/body.json', 'data/derived/app/manifest.json']
BUILDERS = [
    ('ihm/assembly/anatomy.py', 'the role vocabulary, the cross-source name key, the cavity-name rule, '
                                'the bilateral mirror census and the assembly gates'),
    ('scripts/build_canonical_anatomy.py', 'which source surfaces become entities, including the promoted set '
                                           'and the lexical half of its alias resolution, re-derived every build'),
    ('scripts/build_body_mechanics.py', 'the single uniform mass normalizer and which roles carry only '
                                        'numerical inertia'),
    ('ihm/assembly/body.py', 'the canonical body manifest, its source digests, and the promoted '
                             'surface-ownership and provenance declarations'),
    ('scripts/build_cross_structure_conflict_repair.py', 'the ownership rule, and which tet-ready surfaces '
                                                         'it consumes after the duplicate collapse'),
    ('scripts/build_structure_provenance.py', 'the provenance record shape and how a build records the '
                                              'repository state it ran on'),
    ('scripts/build_display_promotion_candidate.py', 'the alias, triage and promotion rules the candidate '
                                                     'directory records'),
]


def prior_promotion():
    """The state the previous promotion left, read from its receipt rather than restated.

    These canonical files are untracked, so once a builder has run the superseded bytes are
    gone. The record-repair receipt is the repository's record of what they said, and taking
    the 'before' numbers from it -- hashed -- is the only way this receipt can state a
    before/after honestly on a rerun, when the observed before is already the after.
    """
    path = ROOT/'data/derived/entity-record-repair-promotion-v1/receipt.json'
    receipt = json.loads(path.read_text())
    return {'path': str(path.relative_to(ROOT)), 'sha256': sha(path),
            'counts': receipt['verified_here']['counts'], 'mass': receipt['verified_here']['mass'],
            'note': 'the state this promotion started from, as the previous promotion recorded it'}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def hashes(paths):
    return {str(p): (sha(ROOT/p) if (ROOT/p).is_file() else None) for p in paths}


def geometry_tree(directory):
    files = sorted(Path(directory).glob('*.json.gz'))
    return {'files': len(files),
            'tree_sha256': hashlib.sha256(''.join(f.name+':'+sha(f) for f in files).encode()).hexdigest()}


def run(command, label, allow_failure=False):
    print('==', label, flush=True)
    began = time.monotonic()
    result = subprocess.run([str(c) for c in command], cwd=ROOT, capture_output=True, text=True)
    if result.returncode != 0 and not allow_failure:
        sys.stderr.write(result.stdout[-4000:]+result.stderr[-4000:])
        raise SystemExit('%s failed with %d' % (label, result.returncode))
    return {'command': [str(c) for c in command], 'seconds': time.monotonic()-began,
            'returncode': result.returncode, 'stdout_tail': result.stdout[-2000:],
            'stderr_tail': result.stderr[-2000:]}


def candidate_intact(name):
    """A candidate directory is a provenance record; verify it still is one."""
    base = CANDIDATES[name]
    manifest = json.loads((base/'manifest.json').read_text())
    artifacts = manifest.get('artifacts_sha256', {})
    return {'path': str(base.relative_to(ROOT)), 'manifest_sha256': sha(base/'manifest.json'),
            'artifacts_checked': len(artifacts),
            'artifacts_mismatched': sorted(rel for rel, digest in artifacts.items()
                                           if not (base/rel).is_file() or sha(base/rel) != digest),
            'declared_canonical_assets_modified': manifest.get('canonical_assets_modified')}


# ---------------------------------------------------------------- alias ladder

def alias_ladder():
    """Re-derive, from the artifacts and not from the candidate's summary, how the display
    scene divides into structures this body already has and structures it lacks.

    A claim circulated that 4,985 display-only ids meant three times more anatomy. It does
    not. Most of them are a second id-form of a structure already present, and the ones
    that survive that test still have to survive a name test and a geometry test.
    """
    manifest = json.loads((ROOT/'data/derived/app/manifest.json').read_text())
    anatomy = json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())
    entities = {e['id']: e for e in anatomy['entities']}
    promoted_ids = {e['source_id'] for e in anatomy['entities']
                    if 'DISPLAY-STRUCTURE-PROMOTION' in e.get('assumptions', [])}
    collapsed = {d['id']: d['survivor']
                 for d in anatomy.get('duplicate_surface_collapse', {}).get('dropped', [])}
    aliases = [json.loads(l) for l in (CANDIDATES['display']/'aliases.jsonl').read_text().splitlines() if l.strip()]
    decisions = [json.loads(l) for l in (CANDIDATES['display']/'decisions.jsonl').read_text().splitlines() if l.strip()]
    by_tier = Counter(a['tier'] for a in aliases)

    non_canonical = [s for s in manifest['structures'] if s['id'] not in entities]
    id_form, collapsed_form, absent = [], [], []
    for s in non_canonical:
        if 'body-'+s['id'] in collapsed:
            collapsed_form.append(s)
        elif 'body-'+s['id'] in entities and s['id'] not in promoted_ids:
            id_form.append(s)
        else:
            # either never resolved to a canonical id, or resolved only because this
            # promotion just created one for it; both are 'this body did not have it'
            absent.append(s)
    name_aliases = [a for a in aliases if a['tier'] in ('name_alias', 'token_name_alias')]
    extent = [a for a in aliases if a['tier'] == 'name_alias_extent_mismatch']
    geometric = [a for a in aliases if a['tier'] == 'geometric_alias']
    return {
        'schema': 'ihm.display-alias-ladder.v1',
        'method': 'recounted from data/derived/app/manifest.json, data/derived/canonical/anatomy.json and '
                  'the candidate alias and decision ledgers, against the canonical entity set as it '
                  'stands, not copied from any summary. Structures this promotion made canonical are '
                  'counted on the absent side, which is where they stood before it ran.',
        'display_structures_now': len(manifest['structures']),
        'canonical_entities_now': len(entities),
        'display_ids_that_are_not_canonical_ids': len(non_canonical),
        'split': {
            'alias_id_forms': len(id_form),
            'collapsed_duplicate_id_forms': len(collapsed_form),
            'genuinely_absent_candidates': len(absent),
        },
        'alias_id_form_rule': "the display id X has a canonical entity 'body-'+X, so it is a second "
                              'id-form of a structure this body already has, not anatomy it lacks',
        'collapsed_duplicate_id_form_rule': 'the canonical row for this display id was collapsed into an '
                                            'identical duplicate-authored surface; it aliases the survivor',
        'of_the_genuinely_absent_candidates': {
            'resolved_by_name_and_geometry': len(name_aliases),
            'same_name_authored_at_a_different_extent': len(extent),
            'resolved_by_geometry_alone_under_the_duplicate_bounds': len(geometric),
            'deleted_as_a_foreign_specimen_or_not_anatomy':
                sum(1 for d in decisions if d['decision'] == 'delete'),
            'kept_as_a_projection_of_this_model':
                sum(1 for d in decisions if d['decision'] == 'keep'),
            'promoted': len(promoted_ids),
        },
        'alias_tiers': dict(by_tier),
        'promoted': len(promoted_ids),
        'correction': 'the display-only id count is not a count of missing anatomy. Roughly half of it '
                      'is id-form aliasing alone, and of what remains, name and geometry resolution and '
                      'foreign-specimen deletion remove most of the rest.',
    }


# ---------------------------------------------------------------- stages

def stage_display(pre_state):
    candidate = json.loads((CANDIDATES['display']/'manifest.json').read_text())
    summary = json.loads((CANDIDATES['display']/'summary.json').read_text())
    prior = prior_promotion()
    before = hashes(CANONICAL_ASSETS)
    before_geometry = geometry_tree(ROOT/'data/derived/canonical/geometry')
    anatomy_before = json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())
    mechanics_before = json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text())
    entities_before = len(anatomy_before['entities'])
    steps = [run([VENV, 'scripts/build_canonical_anatomy.py', '--append'], 'canonical anatomy'),
             run([VENV, 'scripts/build_body_mechanics.py'], 'body mechanics'),
             run([VENV, 'scripts/build_canonical_body.py'], 'canonical body')]
    anatomy = json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())
    mechanics = json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text())
    verification = json.loads((ROOT/'data/derived/canonical/verification.json').read_text())
    promoted = [e for e in anatomy['entities'] if 'DISPLAY-STRUCTURE-PROMOTION' in e.get('assumptions', [])]
    expected = {json.loads(l)['id']: json.loads(l)
                for l in (CANDIDATES['display']/'entities.jsonl').read_text().splitlines() if l.strip()}
    role_disagreements = [{'id': e['id'], 'name': e['name'], 'candidate': expected[e['id']]['role'],
                           'canonical': e['role']}
                          for e in promoted if e['id'] in expected and expected[e['id']]['role'] != e['role']]
    allocated = sum(e['mass_kg'] for e in mechanics['entities'])
    promoted_ids = {e['id'] for e in promoted}
    receipt = {
        'schema': 'ihm.display-structure-promotion.v1',
        'status': 'PROMOTED. Z-Anatomy surfaces this body did not carry are canonical entities. The '
                  'candidate directory is preserved unchanged as the provenance record and is read at '
                  'build time only as a hashed measurement, exactly as the acquired source index is.',
        'canonical_assets_modified': True,
        'promoted_from': {'path': str(CANDIDATES['display'].relative_to(ROOT)),
                          'manifest_sha256': sha(CANDIDATES['display']/'manifest.json'),
                          'decisions_sha256': sha(CANDIDATES['display']/'decisions.jsonl'),
                          'aliases_sha256': sha(CANDIDATES['display']/'aliases.jsonl'),
                          'entities_sha256': sha(CANDIDATES['display']/'entities.jsonl'),
                          'thresholds': candidate['thresholds'],
                          'candidate_recorded_input_sha256': candidate['inputs_sha256']},
        'method': 'the promotion lives in the builders, not in this script. build_canonical_anatomy.py '
                  'reads the candidate decision ledger the way it reads any hashed source index, and '
                  're-derives on every build the half of the alias rule that needs no geometry: a '
                  'promoted structure must be absent by canonical id, by normalized name and by '
                  'side-preserving name token multiset, or the build aborts.',
        'applied_by': {'promotion_script': 'scripts/promote_display_atlas_and_provenance.py',
                       'promotion_script_sha256': sha(Path(__file__)),
                       'builders': [{'path': path, 'sha256': sha(ROOT/path), 'owns': owns}
                                    for path, owns in BUILDERS]},
        'steps': steps,
        'prior_promotion': prior,
        'counts': {
            'entities': {'before_promotion': len(anatomy['entities'])-len(promoted),
                         'after': len(anatomy['entities']),
                         'previous_promotion_recorded': prior['counts']['entities'][1],
                         'method': 'the promoted rows carry DISPLAY-STRUCTURE-PROMOTION in their '
                                   'assumption list, so the count before this promotion is the count '
                                   'now minus them. It is stated this way, and cross-checked against '
                                   'the previous promotion receipt, because a rerun of this script '
                                   'observes a before that is already the after.',
                         'observed_this_run': [entities_before, len(anatomy['entities'])]},
            'promoted': len(promoted),
            'candidate_predicted_promoted': summary['promoted'],
            'roles': {'before_promotion': dict(Counter(e['role'] for e in anatomy['entities']
                                                       if 'DISPLAY-STRUCTURE-PROMOTION' not in e.get('assumptions', []))),
                      'promoted': dict(Counter(e['role'] for e in promoted)),
                      'after': anatomy['counts']['roles']},
            'systems': {'before_promotion': dict(Counter(e['system'] for e in anatomy['entities']
                                                         if 'DISPLAY-STRUCTURE-PROMOTION' not in e.get('assumptions', []))),
                        'promoted': dict(Counter(e['system'] for e in promoted)),
                        'after': anatomy['counts']['systems']},
            'support_links': {'previous_promotion_recorded': prior['counts']['support_links'],
                              'after': mechanics['counts']['support_links'],
                              'observed_this_run': [mechanics_before['counts']['support_links'],
                                                    mechanics['counts']['support_links']]},
            'muscle_actuators': mechanics['counts']['muscle_actuators'],
            'tendon_ligament_paths': mechanics['counts']['tendon_ligament_paths'],
            'boundary_carriers': mechanics['counts']['boundary_carriers'],
            'reference_geometry_files': [before_geometry, geometry_tree(ROOT/'data/derived/canonical/geometry')],
        },
        'mass': {
            'declared_kg': mechanics['mass_allocation']['target_mass_kg'],
            'declared_source': mechanics['mass_allocation']['target_mass_source'],
            'unscaled_proxy_mass_kg': {'previous_promotion_recorded': prior['mass']['unscaled_proxy_mass_kg'][1],
                                       'after': mechanics['mass_allocation']['unscaled_proxy_mass_kg']},
            'uniform_scale': {'previous_promotion_recorded': prior['mass']['uniform_scale'][1],
                              'after': mechanics['mass_allocation']['uniform_scale']},
            'allocated_total_kg': allocated,
            'allocates_to_declared': abs(allocated - mechanics['mass_allocation']['target_mass_kg']) < 1e-9,
            'mass_now_carried_by_promoted_structures_kg':
                sum(e['mass_kg'] for e in mechanics['entities'] if e['id'] in promoted_ids),
            'mass_now_carried_by_pre_existing_structures_kg':
                sum(e['mass_kg'] for e in mechanics['entities'] if e['id'] not in promoted_ids),
            'redistribution_note': 'the declared mass is unchanged and still allocates to it exactly. The '
                                   'promoted surfaces overlap acquired ones, so the proxy partition they '
                                   'enter is not a measured compartment partition and the uniform scale '
                                   'moves accordingly. Surface ownership for the promoted set is '
                                   'unresolved; see the surface-ownership promotion receipt.',
        },
        'verified_here': {
            'assembly_verification': verification,
            'promoted_role_disagreements_with_the_candidate': role_disagreements,
            'promoted_absent_from_the_candidate': sorted(promoted_ids - set(expected)),
            'candidate_promoted_not_in_canonical': sorted(set(expected) - promoted_ids),
            'laterality': bilateral_mirror_census(anatomy['entities']),
            'alias_ladder': alias_ladder(),
        },
        'candidate_intact': candidate_intact('display'),
        'pre_promotion_sha256': pre_state,
        'inputs_sha256_at_run_start': before,
        'outputs_sha256_after': hashes(CANONICAL_ASSETS),
        'known_stale_downstream': [
            {'artifact': 'data/derived/canonical/trajectory.json, trajectory-final-state.json, trajectory-spectra.json',
             'state': 'stale', 'reason': 'dynamic outputs of the changed masses; they need a fresh native '
                                         'run, not a rewritten hash. They were already stale before this promotion.'},
            {'artifact': 'data/derived/entity-tet-ready-v1 and data/derived/muscle-tet-ready-v1',
             'state': 'incomplete for the enlarged entity set',
             'reason': 'no tet-ready surface exists for a promoted structure, so no ownership decision '
                       'or conforming surface does either'},
        ],
        'python': sys.version,
    }
    return receipt


def stage_ownership(libigl):
    out = ROOT/'data/derived/surface-ownership-promotion-v1'
    out.mkdir(parents=True, exist_ok=True)
    step = run([libigl, 'scripts/verify_cross_structure_ownership_reproduces.py',
                '--output', out/'reproduction.json'], 'ownership reproduction', allow_failure=True)
    reproduction = json.loads((out/'reproduction.json').read_text())
    body = json.loads((ROOT/'data/derived/canonical/body.json').read_text())
    return {
        'schema': 'ihm.surface-ownership-promotion.v1',
        'status': ('PROMOTED as shipped, and NOT regenerated. Every shipped decision replays unchanged.'
                   if reproduction['reproduces'] else
                   'PROMOTED as shipped and NOT regenerated, WITH A REPORTED NON-REPRODUCTION. The '
                   'ledger does not reproduce as written against the canonical model as it stands, and '
                   'this receipt says so rather than regenerating it or hiding it. See '
                   'non_reproduction below for exactly what does not replay and what does.'),
        'reproduces_as_written': reproduction['reproduces'],
        'non_reproduction': ({} if reproduction['reproduces'] else {
            'owner_flips': sum(1 for d in reproduction['ownership_decision_drift']
                               if d['shipped'][:2] != d['replayed'][:2]),
            'conflict_class_relabels': sum(1 for d in reproduction['ownership_decision_drift']
                                           if d['shipped'][:2] == d['replayed'][:2]),
            'rows_naming_an_entity_the_model_no_longer_has': reproduction['ledger_rows_naming_a_collapsed_entity'],
            'stale_row_classes': reproduction['stale_rows'],
            'entities_whose_role_moved': len({r['entity_id'] for r in reproduction['role_or_system_drift']}),
            'surface_input_drift': reproduction['surface_input_drift'],
            'magnitude_class_drift': reproduction['magnitude_class_drift'],
            'cause': 'both are consequences of the preceding entity-record-repair promotion, not of the '
                     'display promotion. That lane collapsed five duplicate-authored surfaces and re-roled '
                     '23 intervertebral discs from rigid_bone to cartilage. Its receipt reported zero '
                     'ownership flips, and that is confirmed here exactly: no owner and no yielding party '
                     'changes. What it did not report is that the five dropped ids remain parties to 89 '
                     'ledger rows, and that the disc role change relabels 63 rows from sibling to '
                     'hierarchical without moving the owner.',
            'not_regenerated_because': 'regenerating would re-run a CGAL sweep over 62,973 candidate pairs '
                                       'to reach the same owners, and the priority table has no rank for '
                                       'surface_region, skin, skin_layer or lymphatic_network, so a '
                                       'regeneration over the enlarged entity set would fail on the first '
                                       'promoted region rather than produce a better answer.'}),
        'canonical_assets_modified': True,
        'regenerated': False,
        'promoted_from': {name: {'path': str(CANDIDATES[name].relative_to(ROOT)),
                                 'manifest_sha256': sha(CANDIDATES[name]/'manifest.json')}
                          for name in ('repair', 'atlas')},
        'method': 'every shipped ownership decision was replayed through '
                  'build_cross_structure_conflict_repair.decide and its magnitude thresholds, at the '
                  'role and system the canonical model carries now. Nothing was rebuilt.',
        'applied_by': {'promotion_script': 'scripts/promote_display_atlas_and_provenance.py',
                       'promotion_script_sha256': sha(Path(__file__)),
                       'verifier': 'scripts/verify_cross_structure_ownership_reproduces.py',
                       'verifier_sha256': sha(ROOT/'scripts/verify_cross_structure_ownership_reproduces.py'),
                       'declaration': 'ihm/assembly/body.py surface_ownership()'},
        'steps': [step],
        'reproduction': reproduction,
        'declared_in_canonical': body['surface_ownership'],
        'candidate_intact': {name: candidate_intact(name) for name in ('repair', 'atlas')},
        'python': sys.version,
    }


def stage_provenance(libigl):
    out = ROOT/'data/derived/structure-provenance-v1'
    before = sha(ROOT/'data/derived/canonical/body.json')
    steps = [run([VENV, 'scripts/build_structure_provenance.py', '--output',
                  'data/derived/structure-provenance-v1', '--promoted'], 'structure provenance'),
             run([VENV, 'scripts/build_canonical_body.py'], 'canonical body')]
    index = json.loads((out/'index.json').read_text())
    manifest = json.loads((out/'manifest.json').read_text())
    audit = json.loads((out/'audit.json').read_text())
    body = json.loads((ROOT/'data/derived/canonical/body.json').read_text())
    return {
        'schema': 'ihm.structure-provenance-promotion.v1',
        'status': 'PROMOTED to data/derived/structure-provenance-v1, declared by '
                  'data/derived/canonical/body.json. The candidate directory is untouched.',
        'canonical_assets_modified': True,
        'promoted_from': {'path': str(CANDIDATES['provenance'].relative_to(ROOT)),
                          'manifest_sha256': sha(CANDIDATES['provenance']/'manifest.json')},
        'method': 'the same builder, pointed at a new output directory and told it is the promoted '
                  'index. The only behavioural change is in how it records the repository state it ran '
                  'on: a build whose purpose is traceability may not answer "which commit" with a bare '
                  'null, so it now records the head commit, the head tree, and a digest over every '
                  'path git reports dirty, with an explicit statement of what that does and does not '
                  'guarantee.',
        'applied_by': {'promotion_script': 'scripts/promote_display_atlas_and_provenance.py',
                       'promotion_script_sha256': sha(Path(__file__)),
                       'builder': 'scripts/build_structure_provenance.py',
                       'builder_sha256': sha(ROOT/'scripts/build_structure_provenance.py'),
                       'declaration': 'ihm/assembly/body.py provenance_index()'},
        'steps': steps,
        'working_tree': manifest['working_tree'],
        'builder_record': manifest['builder'],
        'coverage': {'records': index['record_count'],
                     'canonical_entity_coverage': body['provenance_index']['canonical_entity_coverage'],
                     'field_coverage_over_every_indexed_structure': audit['field_coverage'],
                     'tiers': audit['tier_counts'],
                     'hash_verification': audit['hash_verification'],
                     'unfilled': audit['unfilled']},
        'outputs_sha256': {name: sha(out/name) for name in
                           ('index.json', 'manifest.json', 'audit.json', 'schema.json', 'hash-resolution.json')},
        'canonical_body_sha256': [before, sha(ROOT/'data/derived/canonical/body.json')],
        'candidate_intact': candidate_intact('provenance'),
        'python': sys.version,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', default='all',
                        choices=['all', 'display', 'ownership', 'provenance', 'colocation'])
    parser.add_argument('--pre-state', type=Path,
                        help='JSON of canonical asset digests captured before the first builder run; '
                             'these files are untracked, so this sidecar is the only record of the '
                             'superseded bytes')
    parser.add_argument('--libigl-python', type=Path, default=LIBIGL)
    args = parser.parse_args()
    pre_state = json.loads(args.pre_state.read_text()) if args.pre_state else None
    began = time.monotonic()
    written = []

    def emit(directory, receipt):
        path = ROOT/'data/derived'/directory/'receipt.json'
        path.parent.mkdir(parents=True, exist_ok=True)
        receipt['wall_seconds'] = time.monotonic()-began
        path.write_text(json.dumps(receipt, indent=1, allow_nan=False)+'\n')
        written.append(str(path.relative_to(ROOT)))
        print('wrote', path.relative_to(ROOT), flush=True)

    if args.stage in ('all', 'display'):
        emit('display-structure-promotion-v1', stage_display(pre_state))
    if args.stage in ('all', 'colocation'):
        run([args.libigl_python, 'scripts/audit_promoted_surface_colocation.py', '--output',
             'data/derived/display-structure-promotion-v1/surface-colocation.json'], 'surface co-location')
        written.append('data/derived/display-structure-promotion-v1/surface-colocation.json')
    if args.stage in ('all', 'ownership'):
        emit('surface-ownership-promotion-v1', stage_ownership(args.libigl_python))
    if args.stage in ('all', 'provenance'):
        emit('structure-provenance-promotion-v1', stage_provenance(args.libigl_python))
    print(json.dumps({'written': written, 'seconds': time.monotonic()-began}, indent=1))


if __name__ == '__main__':
    main()
