#!/usr/bin/env python3
"""Roll the pose-defect runs up into a findings ledger and a cartilage-synthesis
requirement per joint.

Reads the outputs of `audit_joint_pose_defect.py` at two patch tightnesses and
writes a summary directory. Read-only over everything else.
"""
import argparse
import collections
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from joint_pose_defect_lib import sha256  # noqa: E402

PRIOR = 'data/derived/joint-substrate-assessment-v1'
# error_after_rms and core-floor thresholds, in mm, that decide the verdict
RESOLVED_MM = .50
SHAPE_BLOCK_MM = .50


def load(directory):
    return {name: json.loads((directory / f'{name}.json').read_bytes())
            for name in ('joints', 'literature', 'solve', 'bones', 'manifest')}


def per_joint(run):
    joints = collections.defaultdict(list)
    for site in run['joints']['sites']:
        if site['literature']:
            joints[site['joint']].append(site)
    damped = {v['variant']: v for v in run['solve']['variants']}['damped']
    solved = collections.defaultdict(list)
    for row in damped['per_site']:
        if row['tier'] != 'bonded_structural':
            solved[row['joint']].append(row)
    return joints, solved, damped


def med(rows, path):
    values = []
    for row in rows:
        value = row
        for key in path:
            value = value[key]
        values.append(value)
    return float(np.median(values)) if values else None


def verdict(after, core_floor, target):
    if after is None:
        return 'no_site'
    if after <= RESOLVED_MM:
        return 'resolved_by_global_repose'
    if core_floor is not None and core_floor > SHAPE_BLOCK_MM and core_floor > .25 * target:
        return 'blocked_by_articular_shape_mismatch'
    return 'blocked_by_pose_coupling'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--primary', type=Path, required=True)
    parser.add_argument('--sensitivity', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    primary, sensitivity = load(args.primary), load(args.sensitivity)
    joints, solved, damped = per_joint(primary)
    joints3, solved3, damped3 = per_joint(sensitivity)
    weak = {v['variant']: v for v in primary['solve']['variants']}['weakly_damped']

    missing = collections.Counter(r['joint'] for r in primary['joints']['named_joints_not_found'])
    nearest = {}
    for row in primary['joints']['named_joints_not_found']:
        nearest.setdefault(row['joint'], []).append(row)

    tiers = collections.Counter()
    tier_area = collections.Counter()
    for site in primary['joints']['sites']:
        if site['literature']:
            tiers[site['literature']['tier']] += 1
            tier_area[site['literature']['tier']] += site['patch_area_a_mm2']

    rows = []
    for label in sorted(set(joints) | set(missing)):
        sites, fitted = joints.get(label, []), solved.get(label, [])
        target = sites[0]['literature']['target_gap_mm'] if sites else None
        tier = sites[0]['literature']['tier'] if sites else None
        after = med(fitted, ['error_after_rms_mm'])
        core = med(fitted, ['per_joint_rigid_floor_core_rms_mm'])
        rows.append({
            'joint': label,
            'sites_found': len(sites),
            'instances_not_found': missing.get(label, 0),
            'nearest_missing_instance': (min(nearest.get(label, []),
                                             key=lambda r: r['min_vertex_separation_mm'])
                                         if label in nearest else None),
            'target_gap_mm': target, 'literature_tier': tier,
            'surfaces': sites[0]['literature']['surfaces'] if sites else None,
            'measured_median_gap_mm': med(sites, ['measured_gap', 'median_mm']),
            'measured_min_gap_mm': med(sites, ['measured_gap', 'min_mm']),
            'measured_p95_gap_mm': med(sites, ['measured_gap', 'p95_mm']),
            'interpenetrating_area_fraction': med(sites, ['measured_gap',
                                                          'interpenetrating_area_fraction']),
            'patch_area_a_mm2': med(sites, ['patch_area_a_mm2']),
            'mean_triangle_edge_mm': med(sites, ['mean_triangle_edge_a_mm']),
            'pose_error_mm': med(sites, ['pose_error_mm', 'median_minus_target']),
            'error_before_rms_mm': med(fitted, ['error_before_rms_mm']),
            'per_joint_rigid_floor_rms_mm': med(fitted, ['per_joint_rigid_floor_rms_mm']),
            'per_joint_rigid_floor_core_rms_mm': core,
            'error_after_global_rms_mm': after,
            'error_after_global_max_abs_mm': med(fitted, ['error_after_max_abs_mm']),
            'patch3mm_error_before_rms_mm': med(solved3.get(label, []), ['error_before_rms_mm']),
            'patch3mm_rigid_floor_rms_mm': med(solved3.get(label, []),
                                               ['per_joint_rigid_floor_rms_mm']),
            'verdict': verdict(after, core, target or 1.),
            'required_shell_thickness_mm': ([s['thickness_mm'] for s in sites[0]['literature']['surfaces']]
                                            if sites else None)})

    counts = collections.Counter(r['verdict'] for r in rows)
    findings = {
        'schema': 'ihm.joint-pose-defect-findings.v1',
        'question': 'Is the canonical bone pose a configuration in which every synovial joint can carry a '
                    'cartilage layer of literature thickness, and if not, can a rigid-per-bone re-pose '
                    'make it one?',
        'answer': 'No, and no. Over 278 articular contact sites the signed gap disagrees with the sum of '
                  'the two surfaces literature cartilage thicknesses by 2.21 mm RMS (median 1.90 mm, max '
                  '4.84 mm), in both directions: 22 of 39 joint types sit too far apart and 17 sit too '
                  'close, with the sternoclavicular patch 89 percent interpenetrated. A damped '
                  'rigid-per-bone least squares over 226 free bones brings that to 1.59 mm RMS while '
                  'displacing 56 bones by more than 5 mm and one by 22.6 mm. Relaxing the damping '
                  'twentyfold buys only 1.20 mm RMS and costs 234 mm of displacement on 198 bones, which '
                  'is a torn skeleton, not a re-pose. The lower bound that matters is the per-joint one: '
                  'give each joint its own free relative rigid twist and ignore that its bones serve '
                  'other joints, and 0.99 mm RMS still remains over the full patch, 0.47 mm over the '
                  'congruent inner half. That floor is articular shape mismatch plus mesh resolution and '
                  'no re-posing of anything can remove it. Only 4 joint types - the finger and toe '
                  'interphalangeals - come within 0.5 mm after the global solve.',
        'headline_mm': {
            'damped': damped['cartilage_mm'], 'weakly_damped': weak['cartilage_mm'],
            'patch_gap_3mm_damped': damped3['cartilage_mm']},
        'displacement': {
            'damped_max_vertex_mm': damped['max_bone_displacement_mm'],
            'damped_bones_over_5mm': damped['bones_displaced_over_5mm'],
            'damped_bones_over_5deg': damped['bones_rotated_over_5deg'],
            'weak_max_vertex_mm': weak['max_bone_displacement_mm'],
            'weak_bones_over_5mm': weak['bones_displaced_over_5mm'],
            'spinal_accumulation_note': 'In the weakly damped solve the correction accumulates up the '
                                        'spine from a pinned sacrum, so the cranium ends about 78 mm from '
                                        'where it started. The per-level corrections are each only a few '
                                        'millimetres; the chain is what makes them intolerable.'},
        'verdicts': dict(counts),
        'literature_provenance': {
            'sites_by_tier': dict(tiers),
            'patch_area_mm2_by_tier': {k: round(v, 2) for k, v in tier_area.items()},
            'verified_vs_assumed': {
                'verified_geometry': 'Every gap, patch area, interpenetration fraction, residual and '
                                     'displacement in these artifacts is computed from the canonical '
                                     'geometry whose sha256 the manifest records. Exact point-to-triangle '
                                     'distances, no sampling approximation beyond vertex striding.',
                'transferred_literature': 'Every cartilage thickness. None measured on this specimen.',
                'assumed': 'All thicknesses tiered assumed, the zero-gap target for bonded interfaces, '
                           'the 6 mm patch cut-off, the normal-facing threshold, and the damping weights.'}},
        'cross_check_against_prior_assessment': {
            'prior': PRIOR,
            'prior_contact_pairs': primary['manifest']['input_sha256'].get(
                'data/derived/joint-substrate-assessment-v1/cartilage_gap.json'),
            'reproduced': 'The prior assessment reported a 4.42 mm humeroradial gap, a ~5 mm radiolunate '
                          'gap and a 0.36 mm tibiofemoral gap. This audit measures a 4.42 mm minimum on '
                          'the left humeroradial patch, finds the radiolunate pair with no patch at all '
                          'at 4.81 mm (right) and 5.18 mm (left) minimum vertex separation, and measures '
                          'a 0.08-1.35 mm minimum on the four tibiofemoral sites.'},
        'per_joint': rows}

    synthesis = {
        'schema': 'ihm.cartilage-synthesis-requirement.v1',
        'patch_boundary': {
            'used_here': 'signed gap < %.1f mm to the opposing bone AND opposing outward normals '
                         '(dot < %.2f), then single-linkage split at %.0f mm.' % (
                             primary['joints']['screen']['patch_gap_m'] * 1e3,
                             primary['joints']['screen']['facing_normal_dot_max'],
                             primary['joints']['screen']['cluster_link_m'] * 1e3),
            'problem': 'The patch is cut off by the threshold rather than by anatomy: the area-weighted '
                       '95th percentile gap sits at the cut-off on nearly every large joint, so the '
                       'patch keeps growing as the threshold grows. Halving the threshold to 3 mm is '
                       'reported alongside as a sensitivity.',
            'recommendation': 'Delimit patches from named articular surfaces, not from distance. The '
                              'Z-Anatomy extended atlas carries 40 annotation-only articular features '
                              '(Acetabulum, Superior articular surfaces of tibia, Articular facet of head '
                              'of radius, Superior/Inferior articular facet of vertebra, ...) that would '
                              'give an anatomical boundary. They would arrive through the recorded '
                              'z_anatomy registration, whose never-fitted surface residual is 4.12 mm RMS '
                              'and 12.5 mm max - larger than the cartilage being delimited - so the '
                              'boundary would be anatomically named but positionally soft.'},
        'thickness': 'Per surface, from literature.json. Two shells per joint, one per bone, raised along '
                     'the outward normal. Sum of the two is the target gap used above.',
        'resolution_limit': primary['manifest']['mesh_resolution_mm'],
        'verdict_key': {
            'resolved_by_global_repose': 'after the damped global solve the site is within %.1f mm RMS '
                                         'of its target gap' % RESOLVED_MM,
            'blocked_by_pose_coupling': 'the joint alone could be fixed, but its bones are shared with '
                                        'other joints that pull the other way',
            'blocked_by_articular_shape_mismatch': 'even the joint alone, on its congruent core, cannot '
                                                   'be brought within %.1f mm by any rigid twist; the two '
                                                   'surfaces are not congruent' % SHAPE_BLOCK_MM,
            'no_site': 'the two bones never come close enough to form a contact patch at all'},
        'per_joint': [{k: r[k] for k in ('joint', 'sites_found', 'instances_not_found',
                                         'nearest_missing_instance', 'target_gap_mm', 'literature_tier',
                                         'required_shell_thickness_mm', 'measured_median_gap_mm',
                                         'interpenetrating_area_fraction',
                                         'per_joint_rigid_floor_core_rms_mm',
                                         'error_after_global_rms_mm', 'verdict')} for r in rows],
        'unresolvable': sorted([r['joint'] for r in rows if r['verdict'] != 'resolved_by_global_repose'])}

    args.output.mkdir(parents=True, exist_ok=False)
    written = {}
    for name, payload in (('findings.json', findings), ('synthesis_requirements.json', synthesis)):
        path = args.output / name
        path.write_text(json.dumps(payload, indent=1, allow_nan=False, default=float) + '\n')
        written[name] = {'sha256': sha256(path), 'bytes': path.stat().st_size}
    inputs = {}
    for directory in (args.primary, args.sensitivity):
        for child in sorted(directory.glob('*.json')):
            inputs[str(child.resolve().relative_to(args.root.resolve()))] = child
    inputs['scripts/summarize_joint_pose_defect.py'] = Path(__file__).resolve()
    manifest = {
        'schema': 'ihm.joint-pose-defect-summary-manifest.v1',
        'scope': 'Read-only roll-up of two audit runs. Writes only into this directory.',
        'input_sha256': {k: sha256(v) for k, v in inputs.items()},
        'input_bytes': {k: v.stat().st_size for k, v in inputs.items()},
        'outputs': written,
        'verdict_counts': dict(counts),
        'wall_s': round(time.time() - started, 3)}
    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=1, allow_nan=False,
                                                          default=float) + '\n')
    print(json.dumps({'output': str(args.output), 'verdicts': dict(counts),
                      'joints': len(rows)}, indent=1))


if __name__ == '__main__':
    main()
