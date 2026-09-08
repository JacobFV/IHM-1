#!/usr/bin/env python3
"""Replay every shipped ownership decision against the canonical model as it now stands.

data/derived/cross-structure-repair-v1 is not regenerated. Regeneration would re-run a
CGAL sweep over 62,973 bounding-box-overlapping pairs to reach an answer that is only
allowed to change if one of its actual inputs changed, and it would be a fresh chance to
introduce drift. So this asserts the cheaper and stronger thing instead: that the shipped
ledger is exactly what the shipped rule produces from the current canonical model.

Every decision is replayed through build_cross_structure_conflict_repair.decide and the
magnitude thresholds in that module -- not through a restatement of them here -- using
the role and system each entity carries in data/derived/canonical/anatomy.json today and
the overlap geometry the sweep measured. Any change to a role, a system, a volume or an
identity moves a decision and fails this. Entities the ledger does not cover are counted
and reported rather than assumed to be conflict-free.

Run with the isolated libigl environment, which is what the repair lane imports:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/verify_cross_structure_ownership_reproduces.py
"""
from pathlib import Path
import argparse, hashlib, json, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'scripts'))
from build_cross_structure_conflict_repair import (CONFLICT_FRACTION, CONFLICT_VOLUME_M3, MUSCLE,
                                                   NOISE_FRACTION, NOISE_THICKNESS_M, PRIORITY,
                                                   decide, sha, sources)

REPAIR = ROOT/'data/derived/cross-structure-repair-v1'
ATLAS = ROOT/'data/derived/conflict-free-atlas-v1'
ANATOMY = ROOT/'data/derived/canonical/anatomy.json'


def magnitude(volume, fraction, thickness):
    """The shipped magnitude rule, at the shipped thresholds imported from the builder."""
    if volume is None:
        return 'unmeasured'
    if fraction is not None and (fraction > CONFLICT_FRACTION or volume > CONFLICT_VOLUME_M3):
        return 'modelling_conflict'
    if fraction is not None and fraction <= NOISE_FRACTION and (thickness is None or thickness <= NOISE_THICKNESS_M):
        return 'segmentation_noise'
    return 'substantive'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    anatomy = json.loads(ANATOMY.read_text())
    entities = {e['id']: e for e in anatomy['entities']}
    skipped = []
    rows = sources(skipped)
    by = {r['entity_id']: r for r in rows}
    collapsed = {d['id']: d['survivor']
                 for d in anatomy.get('duplicate_surface_collapse', {}).get('dropped', [])}
    graph = json.loads((REPAIR/'summary.json').read_text())['graph']
    volumes = graph['entity_volume_m3']
    manifest = json.loads((REPAIR/'manifest.json').read_text())

    # 1. the surfaces the lane consumed are byte-identical to what it recorded
    surface_drift = [r['entity_id'] for r in rows
                     if manifest['per_entity_input_sha256'].get(r['entity_id']) != sha(r['path'])]
    surfaces_missing = sorted(set(manifest['per_entity_input_sha256']) - set(by) - set(collapsed))
    surfaces_collapsed_away = sorted(set(manifest['per_entity_input_sha256']) & set(collapsed))

    # 2. every entity the ledger names still exists and still carries the role and system
    #    the decision was taken at
    ledger = [json.loads(line) for line in (REPAIR/'ownership-ledger.jsonl').read_text().splitlines() if line.strip()]
    role_drift, decision_drift, magnitude_drift, missing = [], [], [], []
    # Rows naming a collapsed id are not replayable: the entity is gone. Each is either an
    # exact restatement of a row the survivor already carries -- the dropped surface is
    # byte-identical, so it conflicts with exactly the same partners -- or the
    # survivor-versus-dropped pair itself, which after the collapse is an entity against
    # itself. Both are classified, not skipped.
    pair_key = lambda x, y: tuple(sorted((x, y)))
    present = {pair_key(r['a'], r['b']) for r in ledger}
    stale_rows = {'restated_by_survivor': [], 'self_pair_after_collapse': [], 'unresolved': []}
    for record in ledger:
        a, b = record['a'], record['b']
        if a in collapsed or b in collapsed:
            sa, sb = collapsed.get(a, a), collapsed.get(b, b)
            row = {'pair': [a, b], 'as_survivors': [sa, sb], 'overlap_volume_m3': record['overlap_volume_m3'],
                   'owner': record['owner'], 'conflict_class': record['conflict_class']}
            if sa == sb:
                stale_rows['self_pair_after_collapse'].append(row)
            elif pair_key(sa, sb) in present:
                stale_rows['restated_by_survivor'].append(row)
            else:
                stale_rows['unresolved'].append(row)
            continue
        if a not in entities or b not in entities:
            missing.append((a, b)); continue
        for ident, side in ((a, 'a'), (b, 'b')):
            e = entities[ident]
            if e['role'] != record[side+'_role'] or e['system'] != record[side+'_system']:
                role_drift.append({'entity_id': ident, 'ledger_role': record[side+'_role'],
                                   'ledger_system': record[side+'_system'],
                                   'canonical_role': e['role'], 'canonical_system': e['system']})
        owner, loser, kind = decide({'entity_id': a, 'role': entities[a]['role'], 'system': entities[a]['system']},
                                    {'entity_id': b, 'role': entities[b]['role'], 'system': entities[b]['system']},
                                    volumes)
        if (owner, loser, kind) != (record['owner'], record['yields'], record['conflict_class']):
            decision_drift.append({'pair': [a, b], 'shipped': [record['owner'], record['yields'], record['conflict_class']],
                                   'replayed': [owner, loser, kind]})
        replayed = magnitude(record['overlap_volume_m3'], record['overlap_fraction_of_smaller'],
                             record['overlap_thickness_proxy_m'])
        if replayed != record['magnitude_class']:
            magnitude_drift.append({'pair': [a, b], 'shipped': record['magnitude_class'], 'replayed': replayed})

    # 3. what the ledger does not cover. The lane resolves the surfaces its two tet-ready
    #    parents produced; any canonical entity outside those parents has no ownership
    #    decision at all, and saying so is the point of this number.
    covered = set(by)
    uncovered = sorted(set(entities) - covered)
    uncovered_by_role = {}
    for ident in uncovered:
        uncovered_by_role[entities[ident]['role']] = uncovered_by_role.get(entities[ident]['role'], 0) + 1

    # 4. the atlas built on top of this ledger still matches its own recorded artifacts
    atlas_manifest = json.loads((ATLAS/'manifest.json').read_text())
    atlas_drift = [rel for rel, digest in atlas_manifest['artifacts_sha256'].items()
                   if not (ATLAS/rel).exists() or sha(ATLAS/rel) != digest]
    atlas_input_drift = {rel: {'recorded': digest, 'now': sha(ROOT/rel) if (ROOT/rel).exists() else None}
                         for rel, digest in atlas_manifest['inputs_sha256'].items()
                         if not (ROOT/rel).exists() or sha(ROOT/rel) != digest}

    report = {
        'schema': 'ihm.cross-structure-ownership-reproduction.v1',
        'method': 'every shipped decision replayed through build_cross_structure_conflict_repair.decide '
                  'and its magnitude thresholds, at the role and system canonical carries now',
        'regenerated': False,
        'ledger_decisions': len(ledger),
        'ledger_rows_naming_a_collapsed_entity': sum(len(v) for v in stale_rows.values()),
        'stale_rows': {k: len(v) for k, v in stale_rows.items()},
        'stale_rows_detail': stale_rows,
        'tet_ready_surfaces_skipped_as_collapsed': skipped,
        'surfaces_collapsed_away': surfaces_collapsed_away,
        'entities_covered': len(covered),
        'canonical_entities': len(entities),
        'priority_roles': sorted(PRIORITY),
        'surface_input_drift': surface_drift,
        'surface_inputs_missing': surfaces_missing,
        'entities_missing_from_canonical': missing,
        'role_or_system_drift': role_drift,
        'ownership_decision_drift': decision_drift,
        'magnitude_class_drift': magnitude_drift,
        'uncovered_canonical_entities': len(uncovered),
        'uncovered_by_role': uncovered_by_role,
        'uncovered_note': 'canonical entities with no tet-ready surface in either parent lane, so no '
                          'ownership decision exists for them. Promoted display structures land here '
                          'until the tet-ready lanes are re-run over the enlarged entity set.',
        'roles_absent_from_priority_table': sorted({entities[i]['role'] for i in uncovered} - set(PRIORITY)),
        'atlas_artifact_drift': atlas_drift,
        'atlas_recorded_input_drift': atlas_input_drift,
        'atlas_input_drift_note': 'the atlas records the digest of the whole anatomy.json it read. That '
                                  'file changes whenever any entity anywhere changes, so a digest '
                                  'mismatch is not by itself a change to anything the lane consumed; '
                                  'role_or_system_drift and surface_input_drift are the fields that say '
                                  'whether a consumed value moved.',
    }
    # Reproduction is claimed only if every shipped row replays. Rows naming an entity the
    # canonical model no longer has do not replay, whatever their content, so their presence
    # is a failure to reproduce and is reported as one.
    report['reproduces'] = not (surface_drift or surfaces_missing or missing or role_drift
                                or decision_drift or magnitude_drift or atlas_drift
                                or report['ledger_rows_naming_a_collapsed_entity'])
    report['surviving_decisions_reproduce'] = not (surface_drift or surfaces_missing or missing
                                                   or role_drift or decision_drift or magnitude_drift)
    report['verdict'] = ('every shipped decision replays unchanged'
                         if report['reproduces'] else
                         'the decisions over entities the canonical model still has all replay '
                         'unchanged, but the ledger file itself names entities that no longer '
                         'exist and therefore does not reproduce as written'
                         if report['surviving_decisions_reproduce'] else
                         'a shipped decision does not replay; see the drift fields')
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=1)+'\n')
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ('priority_roles', 'stale_rows_detail')}, indent=1))
    return 0 if report['reproduces'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
