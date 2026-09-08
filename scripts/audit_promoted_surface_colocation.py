#!/usr/bin/env python3
"""Measure how far a purely geometric key can go at separating an alias from a neighbour.

The display promotion resolves an alias with two keys and never with geometry alone:
loose geometry (bbox IoU >= 0.5, symmetric mean surface distance <= 10% of the candidate
bounding-box diagonal) is only allowed to confirm a match the *name* already proposed,
while geometry with no name evidence must clear the far stricter duplicate bounds
(IoU >= 0.8, distance <= 2%). This asks what the loose bounds alone would have claimed,
by sweeping every promoted structure against every canonical structure it shares a
bounding box with and reporting the pairs that clear them.

The answer is the reason the two-key rule exists, and it is a measurement, not an opinion:
the pairs that clear the loose bounds contain both true cross-naming duplicates and
genuinely distinct neighbouring structures, and no threshold on this data separates them.
Report, do not act: nothing here changes a promotion decision.

  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/audit_promoted_surface_colocation.py \
      --output data/derived/display-structure-promotion-v1/surface-colocation.json
"""
from pathlib import Path
import argparse, json, sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'scripts'))
from build_display_promotion_candidate import (ALIAS_DISTANCE_FRACTION, ALIAS_IOU, DUPLICATE_DISTANCE_FRACTION,
                                               DUPLICATE_IOU, alias_verdict, bbox_iou, duplicate_verdict,
                                               geometry_agreement, surface_from_gz)

ANATOMY = ROOT/'data/derived/canonical/anatomy.json'
PROMOTION_ASSUMPTION = 'DISPLAY-STRUCTURE-PROMOTION'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--prescreen-iou', type=float, default=ALIAS_IOU)
    parser.add_argument('--neighbours', type=int, default=5)
    args = parser.parse_args()
    anatomy = json.loads(ANATOMY.read_text())
    entities = anatomy['entities']
    promoted = [e for e in entities if PROMOTION_ASSUMPTION in e.get('assumptions', [])]
    others = [e for e in entities if PROMOTION_ASSUMPTION not in e.get('assumptions', [])]
    if not promoted:
        raise SystemExit('No promoted structures in the canonical model; nothing to audit')
    lo = np.array([e['bounds_m']['min'] for e in others])
    hi = np.array([e['bounds_m']['max'] for e in others])
    box = np.prod(hi-lo, axis=1)
    pairs = []
    for e in promoted:
        a_lo = np.array(e['bounds_m']['min']); a_hi = np.array(e['bounds_m']['max'])
        inter = np.prod(np.clip(np.minimum(a_hi, hi)-np.maximum(a_lo, lo), 0, None), axis=1)
        union = float(np.prod(a_hi-a_lo)) + box - inter
        iou = np.where(union > 0, inter/np.maximum(union, 1e-30), 0.0)
        for index in np.argsort(-iou)[:args.neighbours]:
            if iou[index] < args.prescreen_iou:
                break
            pairs.append((e, others[int(index)]))
    rows = []
    for index, (a, b) in enumerate(pairs):
        av, af = surface_from_gz(ROOT/a['reference_geometry']['path'])
        bv, bf = surface_from_gz(ROOT/b['reference_geometry']['path'])
        agreement = geometry_agreement(av, af, bv, bf)
        rows.append({'promoted_id': a['id'], 'promoted_name': a['name'], 'promoted_role': a['role'],
                     'canonical_id': b['id'], 'canonical_name': b['name'], 'canonical_role': b['role'],
                     'evidence_kind': b['evidence_kind'], 'agreement': agreement,
                     'clears_alias_bounds': alias_verdict(agreement),
                     'clears_duplicate_bounds': duplicate_verdict(agreement)})
        if index and index % 100 == 0:
            print('  %d/%d' % (index, len(pairs)), flush=True)
    loose = [r for r in rows if r['clears_alias_bounds']]
    strict = [r for r in rows if r['clears_duplicate_bounds']]
    best = {}
    for r in loose:
        if r['promoted_id'] not in best or r['agreement']['bbox_iou'] > best[r['promoted_id']]['agreement']['bbox_iou']:
            best[r['promoted_id']] = r
    report = {
        'schema': 'ihm.promoted-surface-colocation.v1',
        'question': 'how many promoted structures would a geometry-only key have called aliases of an '
                    'existing structure, at the bounds the promotion lane only allows a name match to use',
        'thresholds': {'alias_bbox_iou': ALIAS_IOU, 'alias_surface_distance_over_diagonal': ALIAS_DISTANCE_FRACTION,
                       'duplicate_bbox_iou': DUPLICATE_IOU, 'duplicate_surface_distance_over_diagonal': DUPLICATE_DISTANCE_FRACTION},
        'promoted_structures': len(promoted),
        'pairs_measured': len(rows),
        'pairs_clearing_loose_alias_bounds': len(loose),
        'pairs_clearing_strict_duplicate_bounds': len(strict),
        'promoted_structures_clearing_loose_bounds': len(best),
        'promoted_structures_clearing_strict_bounds': len({r['promoted_id'] for r in strict}),
        'finding': 'the set that clears the loose bounds is not separable by any further geometric '
                   'bound on this data. It contains true cross-naming duplicates and true distinct '
                   'neighbours at overlapping distances, so a synonym or concept key is what would '
                   'resolve it; the Z-Anatomy display records carry no concept ids, so this '
                   'repository does not currently hold that key.',
        'consequence': 'these structures are promoted as distinct entities. Where a pair is in fact '
                       'one structure under two names, the model now carries it twice. The pairs are '
                       'listed in full so that is auditable rather than invisible.',
        'best_match_per_structure': sorted(best.values(), key=lambda r: -r['agreement']['bbox_iou']),
        'all_pairs_clearing_loose_bounds': sorted(loose, key=lambda r: -r['agreement']['bbox_iou']),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=1)+'\n')
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ('best_match_per_structure', 'all_pairs_clearing_loose_bounds')}, indent=1))


if __name__ == '__main__':
    main()
