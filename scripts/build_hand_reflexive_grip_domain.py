"""One hand as a single conforming soft-body domain, and what a reflexive grip would still need.

The OpenSim plant welds each hand to the radius as one rigid body: no finger joints, no thumb, no
intrinsic hand muscle. The implicit model has the anatomy anyway, as geometry. This module asks
whether that geometry can be materialised as ONE conforming tetrahedral domain in which flexion
could emerge from tissue and tendon contraction rather than from a declared degree of freedom.

  select    every canonical and joint-substrate structure belonging to one hand, chosen by a
            measured side criterion, clipped exactly to the hand box with a CGAL boolean, with an
            exhaustive conflict census over the clipped members.
  mesh      the members are put in one soup, meshed as a single domain with the interstitial
            complement, tets labelled by winding number, disputed tets awarded to the declared
            role-priority owner, and validated exactly as the prior lanes validated theirs:
            positive volume, shared nodes at interfaces, per-structure volume against each
            surface's divergence integral, residual pairwise overlap, element quality.
  audit     the ownership and undisplaced-fidelity receipts, recomputed from what the mesh stage
            already wrote. Exact, free, and it never re-runs the mesher.
  actuate   material assignment from the sourced candidate table, per-tet fibre transfer from
            data/derived/muscle-fibre-field-v1, the explicit stability limit of
            ihm/assembly/contact_dynamics.py, and, when it fits, a passive settle and a single
            active contraction with conservation receipts.

Requires the isolated libigl environment:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_hand_reflexive_grip_domain.py --self-test
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_hand_reflexive_grip_domain.py \
      --out data/derived/hand-reflexive-grip-domain-v1 --stage all

Nothing outside the --out directory is written. scripts/build_cross_structure_conflict_repair.py,
scripts/build_meshing_route_study.py and every data/derived input are imported and read, never
modified.
"""
from pathlib import Path
import argparse, gzip, hashlib, json, multiprocessing as mp, os, shutil, sys, tempfile, time
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_cross_structure_conflict_repair as R   # noqa: E402  imported, never modified
import build_meshing_route_study as S               # noqa: E402  imported, never modified
import igl                                          # noqa: E402
import igl.copyleft.cgal as cgal                    # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REPAIR = ROOT / 'data/derived/cross-structure-repair-v1'
JOINTS = ROOT / 'data/derived/joint-substrate-candidate-v1'
FIBRE = ROOT / 'data/derived/muscle-fibre-field-v1'
MATERIALS = ROOT / 'data/derived/tissue-material-candidate-v1'
ENVELOPE = ROOT / 'data/derived/outer-envelope/outer-envelope.json.gz'

# ------------------------------------------------------------------ selection policy
#
# The hand box is the axis-aligned bounding box of the 27 bones of one hand, padded. 10 mm is one
# metacarpal diameter: enough to admit the intrinsic muscles and the retinacula that lie over the
# skeleton, short enough that the cut through the extrinsic tendons lands in the distal forearm
# rather than in a muscle belly.
HAND_PAD_M = 0.010
BONE_TOKENS = ('metacarpal bone', 'capitate', 'hamate', 'lunate', 'pisiform', 'scaphoid',
               'trapezium', 'trapezoid', 'triquetral')
# Roles admitted to the domain. Vascular and nerve are excluded by declaration, not by accident:
# they carry no active stress, they are the two roles the priority table itself flags as
# unjustifiable to displace, and admitting them would add ~43 entities of sub-millimetre tubing to
# a domain whose question is whether tendon traction curls a finger. The count that WOULD qualify
# is measured and reported.
ADMITTED_ROLES = ('rigid_bone', 'cartilage', 'tendon', 'ligament', 'muscle', 'connective_tissue')
EXCLUDED_ROLES = ('vascular', 'nerve', 'soft_organ', 'fluid_cavity', 'lymph_node_group')
# A structure enters only if its exact clipped volume clears this. 1 mm3 removes grazing corner
# contact with the box without removing any real hand structure; the smallest retained ligament is
# 1.9 mm3.
MIN_CLIPPED_VOLUME_M3 = 1e-9
# A candidate on the far side of the midline is rejected unless it straddles it.
MIDLINE_TOLERANCE_M = 0.05

# joint_class of the joint-substrate lane mapped onto the priority table of the repair module.
ZA_ROLE = {'ligament': 'ligament', 'retinaculum': 'ligament', 'articular_disc': 'cartilage',
           'tendon_sheath': 'connective_tissue', 'interosseous_membrane': 'connective_tissue',
           'fascia_aponeurosis': 'connective_tissue', 'synovial_bursa': 'fluid_cavity',
           'articular_capsule': 'connective_tissue', 'tendon': 'tendon', 'meniscus': 'cartilage',
           'articular_cartilage': 'cartilage', 'labrum': 'cartilage', 'disc': 'cartilage'}
ENVELOPE_ID = 'hand-outer-envelope'
ENVELOPE_ROLE = 'body_envelope'
# Rank below every tissue: the envelope owns what no structure claims, which is the interstitium.
ENVELOPE_RANK = 11

# ------------------------------------------------------------------ mesher settings
#
# route-decision.json: the exact CGAL arrangement + TetGen -pY route is reproducible only to 48
# entities and its elements are unusable (min dihedral 0.000 deg, one negative-volume tet at 48).
# This domain has 109 members, so route B is the only route that can carry it. epsr is fTetWild's
# envelope, relative to the bounding box diagonal; lr is its target edge length, same normalisation.
# -l is an interior target only: 82 288 input facets whose median edge is far below it force the
# mesh finer than -l wherever a surface runs, so the element count is set by the anatomy, not by -l.
#
# The route study's tightest measured setting, epsr 1e-4 with lr 0.01, was tried first on this
# domain and did not return: FloatTetwild_bin ran 3600 s on the 82 288-facet hand soup and was
# killed by the subprocess timeout with no mesh written. That is a measured ceiling for this host
# at this facet count, not a guess, and it is recorded in domain.json.settings_history.
FINE = {'epsr': 3e-4, 'lr': 0.02}
COARSE = {'epsr': 1e-3, 'lr': 0.05}
SETTINGS_HISTORY = [
    {'epsilon_relative': 1e-3, 'edge_length_relative': 0.05, 'outcome': 'meshed',
     'ftetwild_seconds': 512,
     'note': 'the route study\'s validated 100-entity setting. This is the domain that was built '
             'and validated, on an otherwise quiet host.'},
    {'epsilon_relative': 1e-4, 'edge_length_relative': 0.01, 'outcome': 'timed out',
     'seconds_before_kill': 3600,
     'note': 'the route study\'s tightest thigh setting. On this 109-member, 82 288-facet hand '
             'soup fTetWild had not returned after one hour and was killed; no mesh was written.'},
    {'epsilon_relative': 3e-4, 'edge_length_relative': 0.02, 'outcome': 'timed out',
     'seconds_before_kill': 5403,
     'note': 'the relaxed retry, run while another lane held this 20-core host at load 100-200 so '
             'FloatTetwild_bin averaged 20-27 percent of one core. 5403 s of wall at that share is '
             'roughly 22 min of exclusive CPU, only about 2.5x the coarse run, so this is a '
             'measured wall-clock failure under contention and NOT evidence that the setting is '
             'unreachable on an idle host.'}]
# The route study measured the exact route's interface at 8.96e-16 m and fTetWild's at 2.08e-04 m
# for epsr 1e-4 / lr 0.01 on a thigh box. The hand box diagonal is ~0.24 m, so the same relative
# envelope is ~2.4e-05 m here, and that is the number this lane must beat, not the thigh figure.
SURFACE_DEVIATION_BUDGET_M = S.SURFACE_DEVIATION_BUDGET_M
VOLUME_ERROR_BUDGET = S.VOLUME_ERROR_BUDGET

# ------------------------------------------------------------------ actuation policy
#
# Bone has no modulus anywhere in the sourced table. The soft-body direction document says bones
# are "rigid inclusions or very stiff soft material"; a literature cortical modulus makes the
# explicit stability limit unusable, so both are computed and both are reported.
BONE_E_LITERATURE_PA = 17.0e9      # cortical bone, longitudinal; NOT in the candidate table
BONE_E_SURROGATE_PA = 1.0e7        # declared stiff-soft surrogate, 850x the muscle modulus
POISSON_ASSUMED = 0.45             # blanket assumption recorded in the candidate table
DENSITY_FALLBACK_KG_M3 = 1000.0    # blanket assumption recorded in the candidate table
# Peak active stress of skeletal muscle. Transferred, not measured on this specimen.
ACTIVE_STRESS_PA = 3.0e5


def sha(path):
    return R.sha(path)


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def box_surface(lo, hi):
    """Outward-oriented axis-aligned box. R.BOX is wound inward, which a boolean reads as empty."""
    V = np.ascontiguousarray(np.asarray(lo, float) + R.BOX_CORNERS * (np.asarray(hi, float) - np.asarray(lo, float)))
    F = np.ascontiguousarray(R.BOX[:, [0, 2, 1]].copy())
    return V, F


def load_za(structure_id):
    d = json.loads(gzip.decompress((JOINTS / 'geometry' / (structure_id + '.json.gz')).read_bytes()))
    return (np.ascontiguousarray(np.asarray(d['positions'], float).reshape(-1, 3)),
            np.ascontiguousarray(np.asarray(d['indices'], np.int64).reshape(-1, 3)))


def load_envelope():
    d = json.loads(gzip.decompress(ENVELOPE.read_bytes()))
    return (np.ascontiguousarray(np.asarray(d['positions'], float).reshape(-1, 3)),
            np.ascontiguousarray(np.asarray(d['indices'], np.int64).reshape(-1, 3)))


# --------------------------------------------------------------------- phase: select

_CTX = {}


def _init(rows_by, lo, hi):
    _CTX['by'] = rows_by
    _CTX['lo'] = np.asarray(lo, float)
    _CTX['hi'] = np.asarray(hi, float)


def _source_mesh(kind, key):
    if kind == 'c':
        return R.load_mesh(_CTX['by'][key]['path'])
    if kind == 'z':
        return load_za(key)
    return load_envelope()


def _clip_worker(item):
    kind, key = item
    lo, hi = _CTX['lo'], _CTX['hi']
    V, F = _source_mesh(kind, key)
    full = abs(R.divergence(V, F))
    if (V.min(0) > hi).any() or (V.max(0) < lo).any():
        return {'kind': kind, 'key': key, 'clipped_volume_m3': 0.0, 'full_volume_m3': full,
                'clipped_faces': 0, 'full_faces': int(len(F)), 'outcome': 'outside_box'}
    BV, BF = box_surface(lo, hi)
    try:
        b = cgal.mesh_boolean(V, F, BV, BF, type_str='intersect')
    except BaseException as error:
        return {'kind': kind, 'key': key, 'clipped_volume_m3': None, 'full_volume_m3': full,
                'clipped_faces': 0, 'full_faces': int(len(F)),
                'outcome': 'boolean_failed:' + type(error).__name__}
    CV = np.ascontiguousarray(np.asarray(b[0], float))
    CF = np.ascontiguousarray(np.asarray(b[1], np.int64))
    return {'kind': kind, 'key': key, 'clipped_volume_m3': abs(R.divergence(CV, CF)) if len(CF) else 0.0,
            'full_volume_m3': full, 'clipped_faces': int(len(CF)), 'full_faces': int(len(F)),
            'clipped_vertices': int(len(CV)), 'outcome': 'clipped'}


def hand_bones(anatomy, side):
    sign = 1.0 if side == 'left' else -1.0
    out = []
    for e in anatomy['entities']:
        if e['role'] != 'rigid_bone':
            continue
        n = e['name'].lower()
        if 'phalanx' in n and ('finger' in n or 'thumb' in n):
            pass
        elif any(t in n for t in BONE_TOKENS) and 'toe' not in n and 'metatarsal' not in n:
            pass
        else:
            continue
        if e['centroid_m'][0] * sign > 0:
            out.append(e['id'])
    return sorted(out)


def side_census(rows, rows_by, meta, za_rows, side, workers):
    bones = hand_bones({'entities': list(meta.values())}, side)
    corners = [R.load_mesh(rows_by[b]['path'])[0] for b in bones if b in rows_by]
    lo = np.min([v.min(0) for v in corners], 0) - HAND_PAD_M
    hi = np.max([v.max(0) for v in corners], 0) + HAND_PAD_M
    sign = 1.0 if side == 'left' else -1.0
    items, excluded_roles = [], {}
    for r in rows:
        e = meta[r['entity_id']]
        if e['role'] not in ADMITTED_ROLES:
            excluded_roles.setdefault(e['role'], 0)
            continue
        if e['centroid_m'][0] * sign < 0 and abs(e['centroid_m'][0]) > MIDLINE_TOLERANCE_M:
            continue
        items.append(('c', r['entity_id']))
    for z in za_rows:
        if not (z['promotable'] and z['tet_ready']):
            continue
        c = z['placement']['centroid_m']
        if c[0] * sign < 0 and abs(c[0]) > MIDLINE_TOLERANCE_M:
            continue
        items.append(('z', z['structure_id']))
    items.append(('e', ENVELOPE_ID))
    with mp.get_context('fork').Pool(workers, initializer=_init, initargs=(rows_by, lo, hi)) as pool:
        res = pool.map(_clip_worker, items, chunksize=2)
    kept = [r for r in res if (r['clipped_volume_m3'] or 0.0) > MIN_CLIPPED_VOLUME_M3]
    # the vascular and nerve entities that WOULD have qualified, measured rather than assumed away
    vitems = [('c', r['entity_id']) for r in rows
              if meta[r['entity_id']]['role'] in EXCLUDED_ROLES
              and not (meta[r['entity_id']]['centroid_m'][0] * sign < 0
                       and abs(meta[r['entity_id']]['centroid_m'][0]) > MIDLINE_TOLERANCE_M)]
    with mp.get_context('fork').Pool(workers, initializer=_init, initargs=(rows_by, lo, hi)) as pool:
        vres = pool.map(_clip_worker, vitems, chunksize=4)
    vkept = [r for r in vres if (r['clipped_volume_m3'] or 0.0) > MIN_CLIPPED_VOLUME_M3]
    bone_faces = {b: int(len(R.load_mesh(rows_by[b]['path'])[1])) for b in bones}
    distal = [b for b in bones if 'distal phalanx' in meta[b]['name']]
    return {'side': side, 'box_lo_m': lo.tolist(), 'box_hi_m': hi.tolist(),
            'box_volume_m3': float(np.prod(hi - lo)),
            'box_diagonal_m': float(np.linalg.norm(hi - lo)),
            'hand_bones': bones, 'hand_bone_count': len(bones),
            'hand_bone_faces_total': int(sum(bone_faces.values())),
            'hand_bone_faces_min': int(min(bone_faces.values())),
            'distal_phalanx_faces_total': int(sum(bone_faces[b] for b in distal)),
            'candidates_tested': len(items), 'members': kept,
            'member_count': len(kept),
            'member_clipped_faces': int(sum(r['clipped_faces'] for r in kept)),
            'excluded_by_role_present_in_atlas': sorted(excluded_roles),
            'excluded_role_entities_that_would_qualify': len(vkept),
            'excluded_role_clipped_volume_m3': float(sum(r['clipped_volume_m3'] for r in vkept)),
            'excluded_role_members': [{'entity_id': r['key'], 'name': meta[r['key']]['name'],
                                       'role': meta[r['key']]['role'],
                                       'clipped_volume_m3': r['clipped_volume_m3']} for r in vkept],
            'all_results': res}


def choose_side(left, right):
    """Declared, measured criterion. Both hands carry the same 27 bones and the same 109 candidate
    structures, so the discriminator has to be resolution and recorded conflict, not membership."""
    score = {}
    for c in (left, right):
        score[c['side']] = {'hand_bone_faces_total': c['hand_bone_faces_total'],
                            'hand_bone_faces_min': c['hand_bone_faces_min'],
                            'distal_phalanx_faces_total': c['distal_phalanx_faces_total'],
                            'member_count': c['member_count'],
                            'recorded_conflicts': c.get('recorded_conflicts')}
    a, b = left, right
    better = a if (a['distal_phalanx_faces_total'], a['hand_bone_faces_min']) >= \
                  (b['distal_phalanx_faces_total'], b['hand_bone_faces_min']) else b
    return better['side'], {
        'criterion': 'the more finely sampled distal skeleton wins, tie-broken by the coarsest '
                     'single bone in the hand. Distal phalanges are where a curl has to be '
                     'resolved and they are the coarsest bones in the atlas, so their face count '
                     'sets the resolution ceiling of the whole domain.',
        'scores': score, 'chosen': better['side']}


# --------------------------------------------------------------- phase: conflict census

_PARTS = None


def _conflict_worker(task):
    out = []
    for i, j in task:
        VA, FA = _PARTS[i]
        VB, FB = _PARTS[j]
        rec = {'a': int(i), 'b': int(j)}
        try:
            IF = cgal.intersect_other(VA, FA, VB, FB, True, False, False, False, 2_000_000)[0]
            rec['intersecting_face_pairs'] = int(len(np.asarray(IF).reshape(-1, 2)))
        except BaseException as error:
            rec['error'] = 'intersect_other:' + type(error).__name__
            out.append(rec)
            continue
        if not rec['intersecting_face_pairs']:
            continue
        try:
            b = cgal.mesh_boolean(VA, FA, VB, FB, type_str='intersect')
            IV = np.asarray(b[0], float)
            IFc = np.asarray(b[1], np.int64)
            v = abs(R.divergence(IV, IFc)) if len(IFc) else 0.0
            a = R.area(IV, IFc) if len(IFc) else 0.0
            rec.update(overlap_volume_m3=v, overlap_area_m2=a,
                       overlap_thickness_proxy_m=(2.0 * v / a) if a > 0 else 0.0)
        except BaseException as error:
            rec['overlap_error'] = type(error).__name__
        out.append(rec)
    return out


def conflict_census(parts, priority_of, workers, chunk=48):
    global _PARTS
    _PARTS = [(p['V'], p['F']) for p in parts]
    began = time.monotonic()
    lo = np.stack([V.min(0) for V, _ in _PARTS])
    hi = np.stack([V.max(0) for V, _ in _PARTS])
    ov = (lo[:, None, :] <= hi[None, :, :]).all(2) & (hi[:, None, :] >= lo[None, :, :]).all(2)
    np.fill_diagonal(ov, False)
    pairs = np.argwhere(np.triu(ov))
    tasks = [[(int(a), int(b)) for a, b in pairs[k:k + chunk]] for k in range(0, len(pairs), chunk)]
    found = []
    with mp.get_context('fork').Pool(workers) as pool:
        for batch in pool.imap_unordered(_conflict_worker, tasks):
            found.extend(batch)
    ledger = []
    for rec in found:
        a, b = parts[rec['a']], parts[rec['b']]
        owner, loser, kind = R.decide({'entity_id': a['entity_id'], 'role': a['role'], 'system': a['system']},
                                      {'entity_id': b['entity_id'], 'role': b['role'], 'system': b['system']},
                                      {a['entity_id']: a['surface_volume_m3'],
                                       b['entity_id']: b['surface_volume_m3']}) \
            if a['role'] in R.PRIORITY and b['role'] in R.PRIORITY else (
                (a['entity_id'], b['entity_id'], 'envelope') if priority_of[a['entity_id']] <= priority_of[b['entity_id']]
                else (b['entity_id'], a['entity_id'], 'envelope'))
        ledger.append({'a': a['entity_id'], 'b': b['entity_id'], 'a_name': a['name'], 'b_name': b['name'],
                       'a_role': a['role'], 'b_role': b['role'],
                       'intersecting_face_pairs': rec.get('intersecting_face_pairs', 0),
                       'overlap_volume_m3': rec.get('overlap_volume_m3'),
                       'overlap_thickness_proxy_m': rec.get('overlap_thickness_proxy_m'),
                       'owner': owner, 'yields': loser, 'conflict_class': kind,
                       'error': rec.get('error') or rec.get('overlap_error')})
    ledger.sort(key=lambda r: -(r['overlap_volume_m3'] or 0.0))
    # The envelope encloses every other member by construction, so an envelope pair is a
    # containment relation, not a conflict. It is counted separately and never mixed in.
    def stats(rs):
        v = np.array([r['overlap_volume_m3'] for r in rs if r['overlap_volume_m3'] is not None])
        t = np.array([r['overlap_thickness_proxy_m'] for r in rs
                      if r['overlap_thickness_proxy_m'] is not None])
        seen = {e for r in rs for e in (r['a'], r['b'])}
        return {'pairs': len(rs), 'entities_involved': len(seen),
                'intersecting_face_pairs': int(sum(r['intersecting_face_pairs'] for r in rs)),
                'total_overlap_volume_m3': float(v.sum()) if len(v) else 0.0,
                'overlap_volume_median_m3': float(np.median(v)) if len(v) else None,
                'overlap_volume_max_m3': float(v.max()) if len(v) else None,
                'overlap_thickness_proxy_median_m': float(np.median(t)) if len(t) else None}
    tissue = [r for r in ledger if ENVELOPE_ROLE not in (r['a_role'], r['b_role'])]
    envelope = [r for r in ledger if ENVELOPE_ROLE in (r['a_role'], r['b_role'])]
    inconflict = {e for r in tissue for e in (r['a'], r['b'])}
    return ledger, {
        'members': len(parts), 'all_pairs': len(parts) * (len(parts) - 1) // 2,
        'bounding_box_overlapping_pairs': int(len(pairs)),
        'exactly_tested_pairs': int(len(pairs)),
        'coverage': 'exhaustive over the clipped members: every bounding-box-overlapping pair was '
                    'tested with CGAL intersect_other and, where it intersects, an exact CGAL '
                    'boolean gave the overlap volume. No face cap.',
        'intersecting_pairs': len(ledger),
        'tissue_tissue': stats(tissue),
        'envelope_containment': dict(stats(envelope),
                                     meaning='the outer body envelope encloses every hand '
                                             'structure, so an envelope pair is a containment, not '
                                             'a conflict, and the priority table gives the envelope '
                                             'the lowest rank so it keeps only what no structure '
                                             'claims, which is the interstitium. Only the '
                                             'structures whose SURFACE crosses the envelope surface '
                                             'appear here at all, most of them because they share '
                                             'the flat cut face of the hand box with it; a '
                                             'structure strictly inside the body registers no '
                                             'crossing and is absent from this count.'),
        'entities_in_tissue_conflict': len(inconflict),
        'entities_free_of_tissue_conflict': len(parts) - 1 - len(inconflict),
        'errored_pairs': int(sum(1 for r in ledger if r['error'])),
        'wall_seconds': time.monotonic() - began}


# ------------------------------------------------------------------------ phase: mesh

def priority_table(parts):
    out = {}
    for p in parts:
        out[p['entity_id']] = ENVELOPE_RANK if p['role'] == ENVELOPE_ROLE else R.rank(p['role'])
    return out


def _overlap_worker(task):
    out = []
    for i, j in task:
        VA, FA = _PARTS[i]
        VB, FB = _PARTS[j]
        if not len(FA) or not len(FB):
            continue
        try:
            n = int(len(np.asarray(cgal.intersect_other(VA, FA, VB, FB, True, False, False, False,
                                                        2_000_000)[0]).reshape(-1, 2)))
            b = cgal.mesh_boolean(VA, FA, VB, FB, type_str='intersect')
            v = abs(R.divergence(np.asarray(b[0], float), np.asarray(b[1], np.int64))) if len(b[1]) else 0.0
        except BaseException as error:
            out.append({'a': int(i), 'b': int(j), 'error': type(error).__name__})
            continue
        if n or v:
            out.append({'a': int(i), 'b': int(j), 'shared_interface_face_contacts': n,
                        'residual_overlap_volume_m3': v})
    return out


def residual_overlap_parallel(shells, parts, workers, chunk=32, cap=40):
    global _PARTS
    _PARTS = [(np.ascontiguousarray(V), np.ascontiguousarray(F)) for V, F in shells]
    live = [i for i, (V, F) in enumerate(shells) if len(F)]
    lo = {i: shells[i][0].min(0) for i in live}
    hi = {i: shells[i][0].max(0) for i in live}
    pairs = [(i, j) for a, i in enumerate(live) for j in live[a + 1:]
             if not ((lo[i] > hi[j]).any() or (lo[j] > hi[i]).any())]
    tasks = [pairs[k:k + chunk] for k in range(0, len(pairs), chunk)]
    found = []
    with mp.get_context('fork').Pool(workers) as pool:
        for batch in pool.imap_unordered(_overlap_worker, tasks):
            found.extend(batch)
    total = sum(c.get('residual_overlap_volume_m3', 0.0) for c in found)
    for c in found:
        c['a'] = parts[c['a']]['entity_id']
        c['b'] = parts[c['b']]['entity_id']
    return {'pairs_tested': len(pairs), 'pairs_with_empty_shell': len(parts) - len(live),
            'pairs_in_contact': len(found),
            'total_residual_overlap_volume_m3': float(total),
            'errored_pairs': int(sum(1 for c in found if 'error' in c)),
            'contacts': sorted(found, key=lambda c: -c.get('residual_overlap_volume_m3', 0.0))[:cap],
            'basis': 'exact CGAL boolean intersection of every bounding-box-overlapping pair of '
                     'extracted owned surfaces. A nonzero face-contact count with zero volume is '
                     'the intended outcome: the regions meet on one shared conforming interface.'}


def mesh_domain(parts, priority_of, workdir, epsr, lr, workers, label, full_audit=True,
                timeout=7200):
    workdir = Path(workdir)
    if workdir.exists():
        shutil.rmtree(workdir)
    began = time.monotonic()
    report = {'label': label, 'mesher': 'fTetWild', 'epsilon_relative': epsr,
              'edge_length_relative': lr, 'member_count': len(parts)}
    try:
        rec, TV, TT = S.run_ftetwild(parts, workdir, epsr, lr, timeout=timeout)
    except BaseException as error:
        report['ftetwild'] = {'succeeded': False, 'failed': type(error).__name__,
                              'timeout_s': timeout, 'seconds': time.monotonic() - began,
                              'message': str(error)[:400]}
        report['wall_seconds'] = time.monotonic() - began
        return report, None
    report['ftetwild'] = rec
    if not rec['succeeded']:
        report['wall_seconds'] = time.monotonic() - began
        return report, None
    V, _, _ = S.soup(parts)
    diag = float(np.linalg.norm(V.max(0) - V.min(0)))
    report['envelope'] = {'input_bbox_diagonal_m': diag,
                          'epsilon_absolute_m': epsr * diag,
                          'target_edge_length_m': lr * diag,
                          'note': 'fTetWild normalises both by the bounding box diagonal'}
    mesh, TT, owner, det = S.label_and_measure(TV, TT, parts, priority_of)
    lo = TV.min(0)
    hi = TV.max(0)
    mesh['mesher_box_volume_m3'] = float(np.prod(hi - lo))
    mesh['mesher_box_closure_relative_error'] = float(np.abs(det).sum() / np.prod(hi - lo) - 1)
    report['conforming_mesh'] = mesh
    report['shared_node_audit'] = S.shared_node_audit(TV, TT, owner, parts)
    per = {r['entity_id']: r for r in mesh['per_structure']}
    env = [p['entity_id'] for p in parts if p['role'] == ENVELOPE_ROLE]
    tissue = [r for r in mesh['per_structure'] if r['entity_id'] not in env]
    for r in mesh['per_structure']:
        r['volume_conceded_m3'] = r['claimed_volume_m3'] - r['owned_volume_m3']
        r['fraction_conceded'] = (r['volume_conceded_m3'] / r['claimed_volume_m3']
                                  if r['claimed_volume_m3'] > 0 else None)
    report['ownership'] = {
        'note': 'the outer body envelope claims every structure tet by construction, so '
                'conforming_mesh.disputed_tets counts containment as well as conflict. The '
                'tissue-level numbers below exclude it.',
        'tissue_volume_conceded_m3': float(sum(r['volume_conceded_m3'] for r in tissue)),
        'envelope_volume_conceded_m3': float(sum(per[e]['volume_conceded_m3'] for e in env)),
        'envelope_owned_volume_m3': float(sum(per[e]['owned_volume_m3'] for e in env)),
        'entities_fully_displaced': [r['entity_id'] for r in tissue if r['owned_tets'] == 0],
        'worst_conceders': sorted([{k: r[k] for k in ('entity_id', 'name', 'role',
                                                      'surface_volume_m3', 'claimed_volume_m3',
                                                      'owned_volume_m3', 'fraction_conceded')}
                                   for r in tissue if r['fraction_conceded']],
                                  key=lambda r: -(r['fraction_conceded'] or 0))[:15]}
    if full_audit:
        shells, fid = S.interface_fidelity(TV, TT, owner, parts, det)
        # An entity that conceded nothing has an owned shell that IS its input surface, so its
        # deviation is the mesher's envelope error and nothing else. An entity that conceded volume
        # has an owned shell that runs along the interface with whatever displaced it, which lies
        # inside its input surface by the depth of the overlap, so its deviation measures the
        # conflict, not the mesher. Only the first class is a fidelity measurement.
        undisplaced = {r['entity_id'] for r in tissue
                       if r['owned_tets'] == r['claimed_tets'] and r['claimed_tets'] > 0}
        clean = [r for r in fid['per_structure'] if r.get('faces') and r['entity_id'] in undisplaced]
        fid['undisplaced_reference'] = {
            'entities': len(clean),
            'basis': 'entities whose owned tet set equals their claimed tet set, so no other '
                     'structure displaced them and the deviation is the mesher envelope alone',
            'max_deviation_m': max((r['deviation_max_m'] for r in clean), default=None),
            'median_of_per_entity_max_deviation_m':
                float(np.median([r['deviation_max_m'] for r in clean])) if clean else None,
            'max_abs_volume_relative_error':
                max((abs(r['volume_relative_error'] or 0.0) for r in clean), default=None),
            'median_abs_volume_relative_error':
                float(np.median([abs(r['volume_relative_error'] or 0.0) for r in clean])) if clean else None,
            'worst': sorted([{k: r[k] for k in ('entity_id', 'name', 'role', 'faces',
                                                'deviation_max_m', 'volume_relative_error')}
                             for r in clean], key=lambda r: -r['deviation_max_m'])[:10]}
        report['interface_fidelity'] = fid
        report['resolved_disjointness'] = residual_overlap_parallel(shells, parts, workers)
    report['acceptance'] = {
        'all_tets_positive_volume': bool(mesh['all_positive_volume']),
        'shared_nodes_at_interfaces': bool(report['shared_node_audit']['shared_nodes_at_interfaces']),
        'max_abs_per_structure_volume_error': mesh['max_abs_claimed_volume_error'],
        'volume_error_budget': VOLUME_ERROR_BUDGET,
        'volume_error_within_budget': bool(mesh['max_abs_claimed_volume_error'] < VOLUME_ERROR_BUDGET),
        'max_interface_deviation_m': (report.get('interface_fidelity') or {}).get('max_deviation_m'),
        'max_undisplaced_deviation_m': ((report.get('interface_fidelity') or {})
                                        .get('undisplaced_reference') or {}).get('max_deviation_m'),
        'max_undisplaced_volume_error': ((report.get('interface_fidelity') or {})
                                         .get('undisplaced_reference') or {}).get('max_abs_volume_relative_error'),
        'surface_deviation_budget_m': SURFACE_DEVIATION_BUDGET_M,
        'residual_pairwise_overlap_m3': (report.get('resolved_disjointness') or {}).get('total_residual_overlap_volume_m3'),
        'min_dihedral_deg': mesh['quality']['min_dihedral_deg'],
        'tets_with_min_dihedral_under_5deg': mesh['quality']['tets_with_min_dihedral_under_5deg']}
    report['wall_seconds'] = time.monotonic() - began
    return report, {'TV': TV, 'TT': TT, 'owner': owner, 'det': det}


def exact_route_control(parts, priority_of, pad, flags, timeout_note):
    """The exact CGAL arrangement + TetGen -pY route, run once so its failure is a receipt.

    TetGen writes tetgen-tmpfile_skipped.* into the working directory when it skips a facet, so the
    call runs from a scratch directory and leaves nothing in the repository."""
    began = time.monotonic()
    cwd = os.getcwd()
    scratch = tempfile.mkdtemp(prefix='hand-exact-')
    try:
        os.chdir(scratch)
        report, mesh = R.mesh_cluster(parts, pad, flags, priority_of, emit=None,
                                      label='hand-exact-control', snap_below=0.0)
    except BaseException as error:
        return {'route': '0-exact-arrangement-control', 'member_count': len(parts),
                'failed': type(error).__name__, 'message': str(error),
                'wall_seconds': time.monotonic() - began, 'note': timeout_note}
    finally:
        os.chdir(cwd)
        shutil.rmtree(scratch, ignore_errors=True)
    report['route'] = '0-exact-arrangement-control'
    report['wall_seconds'] = time.monotonic() - began
    report.pop('resolved_surfaces', None)
    return report


# --------------------------------------------------------------------- phase: actuate

def material_table():
    m = json.loads((MATERIALS / 'materials.json').read_text())
    return m['materials'], m.get('entity_overrides', {})


def lame(young_pa, nu):
    mu = young_pa / (2.0 * (1.0 + nu))
    lam = young_pa * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return mu, lam


def assign_materials(parts, owner, materials, bone_e_pa):
    """Per-tet mu, lambda, rho with an explicit sourced/assumed tier for every number used."""
    rows = []
    mu = np.zeros(len(owner))
    lam = np.zeros(len(owner))
    rho = np.zeros(len(owner))
    for i, p in enumerate(parts):
        role = p['role']
        entry = materials.get(role if role != ENVELOPE_ROLE else 'adipose', {})
        lin = entry.get('linear', {})
        e = (lin.get('young_modulus') or {}).get('value')
        tier = (lin.get('young_modulus') or {}).get('tier')
        source = (lin.get('young_modulus') or {}).get('source')
        if role == 'rigid_bone':
            e, tier, source = bone_e_pa, 'assumed', 'declared in this module; the candidate table has no bone modulus'
        if e is None:
            g = (lin.get('shear_modulus') or {}).get('value')
            if g is not None:
                e = 3.0 * g
                tier, source = 'derived', (lin.get('shear_modulus') or {}).get('source')
        if e is None:
            e = materials['muscle']['linear']['young_modulus']['value']
            tier, source = 'assumed', 'no value in the candidate table for this role; muscle modulus substituted'
        d = (lin.get('density') or {}).get('value') or DENSITY_FALLBACK_KG_M3
        dtier = (lin.get('density') or {}).get('tier') or 'assumed'
        m_, l_ = lame(float(e), POISSON_ASSUMED)
        sel = owner == i
        mu[sel] = m_
        lam[sel] = l_
        rho[sel] = float(d)
        rows.append({'entity_id': p['entity_id'], 'name': p['name'], 'role': role,
                     'tets': int(sel.sum()), 'young_modulus_pa': float(e),
                     'young_modulus_tier': tier, 'young_modulus_source': source,
                     'density_kg_m3': float(d), 'density_tier': dtier,
                     'poisson_ratio': POISSON_ASSUMED, 'poisson_tier': 'assumed'})
    return mu, lam, rho, rows


def transfer_fibres(parts, TV, TT, owner):
    """Nearest-source-tet transfer of the per-tet fibre direction field onto this domain."""
    from scipy.spatial import cKDTree
    index = {json.loads(l)['entity_id']: json.loads(l)
             for l in (FIBRE / 'entities.jsonl').read_text().splitlines()}
    centre = TV[TT].mean(1)
    direction = np.zeros((len(TT), 3))
    rows = []
    for i, p in enumerate(parts):
        if p['role'] != 'muscle':
            continue
        eid = p['entity_id']
        sel = np.flatnonzero(owner == i)
        rec = index.get(eid)
        if rec is None or not len(sel):
            rows.append({'entity_id': eid, 'name': p['name'], 'tets': int(len(sel)),
                         'provenance': None, 'transferred': False,
                         'reason': 'no fibre field for this entity' if rec is None else 'owns no tet'})
            continue
        d = np.load(FIBRE / rec['fibre_path'])
        src = d['tet_vertices_m'][d['tets']].mean(1)
        tree = cKDTree(src)
        dist, idx = tree.query(centre[sel])
        f = np.asarray(d['fibre'], float)[idx]
        f /= np.maximum(np.linalg.norm(f, axis=1, keepdims=True), 1e-300)
        direction[sel] = f
        h = np.cbrt(np.abs(np.linalg.det(np.swapaxes(TV[TT[sel, 1:]] - TV[TT[sel, 0, None]], 1, 2)) / 6.0))
        rows.append({'entity_id': eid, 'name': p['name'], 'tets': int(len(sel)),
                     'provenance': rec['provenance'], 'transferred': True,
                     'source_tets': int(len(src)),
                     'transfer_distance_median_m': float(np.median(dist)),
                     'transfer_distance_max_m': float(dist.max()),
                     'transfer_distance_over_element_size_median': float(np.median(dist / np.maximum(h, 1e-300))),
                     'fibre_coherence_source': rec.get('fibre_coherence')})
    return direction, rows


def active_force(region, direction, sigma_pa, positions):
    """Nodal forces of an active fibre stress P_act = sigma * (F a0) (x) a0, assembled exactly as
    DeformableRegion assembles its passive first Piola term."""
    f = region.deformation(positions)
    fa = np.einsum('nij,nj->ni', f, direction)
    piola = sigma_pa[:, None, None] * fa[:, :, None] * direction[:, None, :]
    h = region.volumes[:, None, None] * piola @ np.swapaxes(region.inverse, 1, 2)
    local = np.concatenate((-h.sum(axis=2)[:, None, :], np.swapaxes(h, 1, 2)), axis=1)
    out = np.zeros_like(positions)
    for k in range(4):
        np.add.at(out, region.tets[:, k], local[:, k])
    return -out


def build_body(TV, TT, owner, mu, lam, rho, direction, fixed_mask):
    """Only the tets inside the body envelope enter the solid; the exterior box tets are air."""
    sys.path.insert(0, str(ROOT))
    from ihm.assembly.contact_dynamics import DynamicTetrahedra
    keep = owner >= 0
    tets = TT[keep]
    used = np.unique(tets)
    remap = np.full(len(TV), -1, np.int64)
    remap[used] = np.arange(len(used))
    V = np.ascontiguousarray(TV[used])
    T = np.ascontiguousarray(remap[tets.ravel()].reshape(-1, 4))
    det = np.linalg.det(np.swapaxes(V[T[:, 1:]] - V[T[:, 0, None]], 1, 2))
    flip = det < 0
    if flip.any():
        T[flip] = T[flip][:, [0, 2, 1, 3]]
    body = DynamicTetrahedra(V, T, mu_pa=mu[keep], lambda_pa=lam[keep], density_kg_m3=rho[keep],
                             fixed_nodes=np.flatnonzero(fixed_mask[used]).astype(int))
    return body, keep, used, T, direction[keep]


def run_dynamics(body, direction, active_sigma, steps, dt, gravity):
    sys.path.insert(0, str(ROOT))
    from ihm.assembly.contact_dynamics import step_coupled
    tip = np.argmax(-body.reference_position_m[:, 1])
    began = time.monotonic()
    trace = []
    x0 = body.position_m.copy()
    for k in range(steps):
        ext = None
        if active_sigma is not None:
            ext = {'hand': active_force(body.region, direction, active_sigma, body.position_m)}
        rec = step_coupled({'hand': body}, dt, external_forces_n=ext, gravity_m_s2=gravity)
        if k % max(1, steps // 10) == 0 or k == steps - 1:
            trace.append({'step': k, 'time_s': rec['time_s'],
                          'elastic_energy_j': rec['elastic_energy_j'],
                          'kinetic_energy_j': rec['kinetic_energy_j'],
                          'total_energy_j': rec['total_energy_j'],
                          'external_work_j': rec['external_work_j'],
                          'numerical_energy_defect_j': rec['numerical_energy_defect_j'],
                          'momentum_residual_ns': [float(v) for v in rec['momentum_residual_ns']],
                          'max_displacement_m': float(np.abs(body.position_m - x0).max()),
                          'min_jacobian': body.minimum_jacobian()})
    d = body.position_m - x0
    return {'steps': steps, 'dt_s': dt, 'simulated_time_s': steps * dt,
            'wall_seconds': time.monotonic() - began,
            'wall_seconds_per_step': (time.monotonic() - began) / steps,
            'max_node_displacement_m': float(np.linalg.norm(d, axis=1).max()),
            'mean_node_displacement_m': float(np.linalg.norm(d, axis=1).mean()),
            'most_distal_node_displacement_m': float(np.linalg.norm(d[tip])),
            'min_jacobian': body.minimum_jacobian(),
            'trace': trace}


# ------------------------------------------------------------------------ self-test

def self_test():
    checks = []

    def check(name, ok, detail):
        checks.append({'check': name, 'passed': bool(ok), 'detail': detail})

    lo = np.zeros(3)
    hi = np.ones(3)
    BV, BF = box_surface(lo, hi)
    check('box surface is outward oriented', abs(R.divergence(BV, BF) - 1.0) < 1e-12,
          R.divergence(BV, BF))

    # a unit cube clipped by a half-box keeps exactly half its volume
    CV, CF = box_surface([0, 0, 0], [1, 1, 1])
    HV, HF = box_surface([0, 0, 0], [0.5, 1, 1])
    b = cgal.mesh_boolean(CV, CF, HV, HF, type_str='intersect')
    v = abs(R.divergence(np.asarray(b[0], float), np.asarray(b[1], np.int64)))
    check('exact clip keeps exactly half the cube', abs(v - 0.5) < 1e-12, v)

    # the envelope rank is below every tissue rank in the imported priority table
    check('envelope ranks below every tissue role',
          ENVELOPE_RANK > max(v[0] for v in R.PRIORITY.values()),
          (ENVELOPE_RANK, max(v[0] for v in R.PRIORITY.values())))
    check('every joint_class maps to a declared priority role',
          set(ZA_ROLE.values()) <= set(R.PRIORITY), sorted(set(ZA_ROLE.values()) - set(R.PRIORITY)))

    # Lame conversion round trip
    mu, lam = lame(11730.0, 0.45)
    e = mu * (3 * lam + 2 * mu) / (lam + mu)
    check('Lame conversion round trips Young', abs(e - 11730.0) < 1e-6, (mu, lam, e))

    # ownership, labelling and validation on two overlapping cubes meshed by the real mesher
    parts = []
    for eid, role, box in (('A', 'rigid_bone', ([0, 0, 0], [1, 1, 1])),
                           ('B', 'muscle', ([0.75, 0, 0], [1.75, 1, 1]))):
        V, F = box_surface(*box)
        parts.append({'entity_id': eid, 'name': eid, 'role': role, 'system': 'test',
                      'frame': 'test', 'source_path': 'test', 'source_sha256': '0' * 64,
                      'surface_volume_m3': abs(R.divergence(V, F)), 'V': V, 'F': F})
    priority = priority_table(parts)
    check('bone outranks muscle in the domain priority', priority['A'] < priority['B'], priority)
    ledger, census = conflict_census(parts, priority, 2)
    check('census finds the one overlapping pair', census['tissue_tissue']['pairs'] == 1, census)
    check('census measures the 0.25 m3 slab',
          abs((ledger[0]['overlap_volume_m3'] or 0) - 0.25) < 1e-9, ledger[0])
    check('bone is awarded the disputed slab', ledger[0]['owner'] == 'A', ledger[0])

    with tempfile.TemporaryDirectory() as scratch:
        report, mesh = mesh_domain(parts, priority, Path(scratch) / 'ftw', 1e-3, 0.05, 2,
                                   'self-test')
        check('fTetWild meshed the self-test domain', report['ftetwild']['succeeded'],
              report['ftetwild'].get('stderr_tail'))
        if mesh is not None:
            m = report['conforming_mesh']
            check('all tets positive', m['all_positive_volume'], m['min_tet_volume_m3'])
            check('shared nodes at interfaces',
                  report['shared_node_audit']['shared_nodes_at_interfaces'],
                  report['shared_node_audit'])
            a = [p for p in m['per_structure'] if p['entity_id'] == 'A'][0]
            b2 = [p for p in m['per_structure'] if p['entity_id'] == 'B'][0]
            check('bone keeps its whole 1 m3', abs(a['owned_volume_m3'] - 1.0) < 5e-2, a)
            check('muscle concedes the slab and keeps 0.75 m3',
                  abs(b2['owned_volume_m3'] - 0.75) < 5e-2, b2)
            check('residual pairwise overlap is zero',
                  report['resolved_disjointness']['total_residual_overlap_volume_m3'] < 1e-12,
                  report['resolved_disjointness']['total_residual_overlap_volume_m3'])

            # a stiff cube under gravity, one explicit step, must conserve what step_coupled claims
            owner = mesh['owner']
            mu_a, lam_a, rho_a, _ = assign_materials(parts, owner, material_table()[0], 1e6)
            direction = np.tile(np.array([1.0, 0.0, 0.0]), (len(mesh['TT']), 1))
            inside = np.unique(mesh['TT'][owner >= 0])
            fixed = mesh['TV'][:, 1] < mesh['TV'][inside, 1].min() + 1e-6
            body, keep, used, T, dirk = build_body(mesh['TV'], mesh['TT'], owner, mu_a, lam_a,
                                                   rho_a, direction, fixed)
            dt = 0.5 * body.max_explicit_dt_s
            passive = run_dynamics(body, dirk, None, 3, dt, (0.0, -9.81, 0.0))
            last = passive['trace'][-1]
            check('passive step conserves energy to the reported defect',
                  abs(last['numerical_energy_defect_j']) < 1e-6 * max(abs(last['total_energy_j']), 1.0),
                  last)
            check('passive step conserves momentum',
                  float(np.abs(last['momentum_residual_ns']).max()) < 1e-9,
                  last['momentum_residual_ns'])
            # An active fibre stress is not derived from a potential, so step_coupled's energy
            # defect absorbs the work the scheme mis-integrates rather than going to zero. The
            # substantive test is that it converges: halving dt must shrink the per-step defect.
            sigma = np.full(len(dirk), 1e4)
            saved = body.checkpoint()
            before = body.position_m.copy()
            active = run_dynamics(body, dirk, sigma, 4, dt, (0.0, 0.0, 0.0))
            check('active fibre stress moves the solid',
                  float(np.abs(body.position_m - before).max()) > 0,
                  active['max_node_displacement_m'])
            coarse_defect = abs(active['trace'][-1]['numerical_energy_defect_j'])
            body.restore(saved)
            refined = run_dynamics(body, dirk, sigma, 8, dt / 2, (0.0, 0.0, 0.0))
            fine_defect = abs(refined['trace'][-1]['numerical_energy_defect_j'])
            ratio = fine_defect / coarse_defect if coarse_defect else 0.0
            check('active energy defect converges when dt is halved', ratio < 0.8,
                  {'defect_at_dt_j': coarse_defect, 'defect_at_half_dt_j': fine_defect,
                   'ratio': ratio})
            check('active step conserves momentum',
                  float(np.abs(refined['trace'][-1]['momentum_residual_ns']).max()) < 1e-9,
                  refined['trace'][-1]['momentum_residual_ns'])

    passed = all(c['passed'] for c in checks)
    print(json.dumps({'self_test': 'ihm.hand-reflexive-grip-domain', 'passed': passed,
                      'checks': checks}, indent=2))
    return passed


# ------------------------------------------------------------------------------ run

def build_parts(side_report, rows_by, meta, za_by, keep_ids=None):
    parts = []
    lo = np.asarray(side_report['box_lo_m'])
    hi = np.asarray(side_report['box_hi_m'])
    BV, BF = box_surface(lo, hi)
    for m in side_report['members']:
        key = m['key']
        if keep_ids is not None and key not in keep_ids:
            continue
        if m['kind'] == 'c':
            V, F = R.load_mesh(rows_by[key]['path'])
            name, role, system = meta[key]['name'], meta[key]['role'], meta[key]['system']
            src, digest = str(rows_by[key]['path'].relative_to(ROOT)), rows_by[key]['sha256']
        elif m['kind'] == 'z':
            V, F = load_za(key)
            z = za_by[key]
            name, role, system = z['name'], ZA_ROLE[z['joint_class']], z['system']
            src = 'data/derived/joint-substrate-candidate-v1/' + z['output_path']
            digest = z['output_sha256']
        else:
            V, F = load_envelope()
            name, role, system = 'hand outer body envelope', ENVELOPE_ROLE, 'integumentary'
            src, digest = str(ENVELOPE.relative_to(ROOT)), sha(ENVELOPE)
        b = cgal.mesh_boolean(V, F, BV, BF, type_str='intersect')
        CV = np.ascontiguousarray(np.asarray(b[0], float))
        CF = np.ascontiguousarray(np.asarray(b[1], np.int64))
        parts.append({'entity_id': key, 'name': name, 'role': role, 'system': system,
                      'lane': {'c': 'canonical', 'z': 'joint-substrate-candidate-v1',
                               'e': 'outer-envelope'}[m['kind']],
                      'frame': 'bodyparts3d-display-m', 'source_path': src, 'source_sha256': digest,
                      'full_volume_m3': m['full_volume_m3'], 'full_faces': m['full_faces'],
                      'clipped_faces': int(len(CF)), 'clipped_vertices': int(len(CV)),
                      'surface_volume_m3': abs(R.divergence(CV, CF)),
                      'clipped': m['full_faces'] != int(len(CF)) or m['clipped_volume_m3'] < m['full_volume_m3'] * (1 - 1e-9),
                      'V': CV, 'F': CF})
    parts.sort(key=lambda p: p['entity_id'])
    return parts


def stage_select(out, args):
    began = time.monotonic()
    rows = R.sources()
    rows_by = {r['entity_id']: r for r in rows}
    anatomy = json.loads((ROOT / 'data/derived/canonical/anatomy.json').read_text())
    meta = {e['id']: e for e in anatomy['entities']}
    za_rows = [json.loads(l) for l in (JOINTS / 'entities.jsonl').read_text().splitlines()]
    za_by = {z['structure_id']: z for z in za_rows}
    ledger = [json.loads(l) for l in (REPAIR / 'ownership-ledger.jsonl').read_text().splitlines()]
    tetready = {}
    for lane in ('muscle-tet-ready-v1', 'entity-tet-ready-v1'):
        for l in (ROOT / 'data/derived' / lane / 'entities.jsonl').read_text().splitlines():
            r = json.loads(l)
            tetready[r['entity_id']] = r
    sides = {}
    for side in ('left', 'right'):
        c = side_census(rows, rows_by, meta, za_rows, side, args.workers)
        ids = {m['key'] for m in c['members'] if m['kind'] == 'c'}
        edges = [e for e in ledger if e['a'] in ids and e['b'] in ids]
        c['recorded_conflicts'] = len(edges)
        c['recorded_conflict_overlap_volume_m3'] = float(sum(e['overlap_volume_m3'] or 0.0 for e in edges))
        c['recorded_conflicts_flagged'] = int(sum(e['flagged_for_review'] for e in edges))
        c['recorded_conflict_magnitude'] = {k: sum(1 for e in edges if e['magnitude_class'] == k)
                                            for k in ('segmentation_noise', 'substantive',
                                                      'modelling_conflict', 'unmeasured')}
        c['members_not_tet_ready'] = [{'entity_id': m['key'], 'name': meta[m['key']]['name'],
                                       'closed_after_repair': tetready[m['key']]['after']['closed']}
                                      for m in c['members']
                                      if m['kind'] == 'c' and not tetready[m['key']]['tet_ready']]
        sides[side] = c
    side, why = choose_side(sides['left'], sides['right'])
    chosen = sides[side]
    # drop members that no volumetric mesher can take, and say what replaced them
    drop = {}
    for m in list(chosen['members']):
        if m['kind'] == 'c' and not tetready[m['key']]['tet_ready']:
            drop[m['key']] = {'entity_id': m['key'], 'name': meta[m['key']]['name'],
                              'reason': 'not tet-ready in its lane: the repaired surface is still '
                                        'open, so it bounds no volume',
                              'closed_after_repair': tetready[m['key']]['after']['closed']}
    replacement = {}
    for k in drop:
        n = meta[k]['name'].lower()
        for z in chosen['members']:
            if z['kind'] != 'z':
                continue
            zn = za_by[z['key']]['name'].lower()
            if 'retinaculum' in n and 'retinaculum' in zn and 'flexor' in n and 'flexor' in zn:
                replacement[k] = {'replaced_by': z['key'], 'replacement_name': za_by[z['key']]['name'],
                                  'basis': 'the joint-substrate lane carries the same structure as '
                                           'a closed, tet-ready surface'}
    keep_ids = {m['key'] for m in chosen['members'] if m['key'] not in drop}
    parts = build_parts(chosen, rows_by, meta, za_by, keep_ids)
    priority = priority_table(parts)
    conflicts, census = conflict_census(parts, priority, args.workers)
    # capsule and aponeurosis candidates the joint-substrate lane has but cannot promote
    unpromotable = []
    lo = np.asarray(chosen['box_lo_m'])
    hi = np.asarray(chosen['box_hi_m'])
    sign = 1.0 if side == 'left' else -1.0
    for z in za_rows:
        if z['promotable'] and z['tet_ready']:
            continue
        c = z['placement']['centroid_m']
        if c[0] * sign < 0 or not ((np.asarray(c) >= lo).all() and (np.asarray(c) <= hi).all()):
            continue
        unpromotable.append({'structure_id': z['structure_id'], 'name': z['name'],
                             'joint_class': z['joint_class'], 'promotable': z['promotable'],
                             'tet_ready': z['tet_ready'],
                             'closed_after_repair': z['after']['closed'],
                             'boundary_edges_after_repair': z['after']['boundary_edges']})
    by_role = {}
    for p in parts:
        by_role.setdefault(p['role'], []).append(p['entity_id'])
    by_lane = {}
    for p in parts:
        by_lane[p['lane']] = by_lane.get(p['lane'], 0) + 1
    selection = {
        'schema': 'ihm.hand-reflexive-grip-selection.v1',
        'side_chosen': side, 'side_criterion': why,
        'side_census': {s: {k: v for k, v in c.items() if k not in ('all_results', 'members',
                                                                    'excluded_role_members')}
                        for s, c in sides.items()},
        'box_lo_m': chosen['box_lo_m'], 'box_hi_m': chosen['box_hi_m'],
        'box_volume_m3': chosen['box_volume_m3'], 'box_diagonal_m': chosen['box_diagonal_m'],
        'box_padding_m': HAND_PAD_M,
        'entities': len(parts),
        'entities_by_role': {k: len(v) for k, v in sorted(by_role.items())},
        'entities_by_lane': by_lane,
        'total_clipped_faces': int(sum(p['clipped_faces'] for p in parts)),
        'total_unclipped_faces': int(sum(p['full_faces'] for p in parts)),
        'entities_clipped_by_the_box': int(sum(1 for p in parts if p['clipped'])),
        'total_clipped_volume_m3': float(sum(p['surface_volume_m3'] for p in parts)),
        'dropped_members': list(drop.values()),
        'dropped_member_replacements': replacement,
        'excluded_roles': list(EXCLUDED_ROLES),
        'excluded_role_entities_that_would_qualify': chosen['excluded_role_entities_that_would_qualify'],
        'excluded_role_clipped_volume_m3': chosen['excluded_role_clipped_volume_m3'],
        'joint_substrate_candidates_in_the_box_that_cannot_be_promoted': unpromotable,
        'recorded_conflicts_among_canonical_members': chosen['recorded_conflicts'],
        'recorded_conflict_source': 'data/derived/cross-structure-repair-v1/ownership-ledger.jsonl',
        'recorded_conflict_magnitude': chosen['recorded_conflict_magnitude'],
        'recorded_conflicts_flagged': chosen['recorded_conflicts_flagged'],
        'recorded_conflict_overlap_volume_m3': chosen['recorded_conflict_overlap_volume_m3'],
        'measured_conflicts_over_all_clipped_members': census,
        'selection_policy': {
            'hand_pad_m': HAND_PAD_M, 'min_clipped_volume_m3': MIN_CLIPPED_VOLUME_M3,
            'admitted_roles': list(ADMITTED_ROLES),
            'clip': 'every member is intersected with the hand box by an exact CGAL boolean, so a '
                    'structure crossing the wrist enters as exactly the part of itself inside the '
                    'domain and the cut face is a real boundary where proximal traction enters',
            'za_role_map': ZA_ROLE, 'envelope_rank': ENVELOPE_RANK},
        'wall_seconds': time.monotonic() - began}
    write_json(out / 'selection.json', selection)
    with (out / 'members.jsonl').open('w') as f:
        for p in parts:
            f.write(json.dumps({k: v for k, v in p.items() if k not in ('V', 'F')},
                               allow_nan=False) + '\n')
    with (out / 'conflicts.jsonl').open('w') as f:
        for c in conflicts:
            f.write(json.dumps(c, allow_nan=False) + '\n')
    np.savez_compressed(out / 'clipped-members.npz',
                        **{'V%d' % i: p['V'] for i, p in enumerate(parts)},
                        **{'F%d' % i: p['F'] for i, p in enumerate(parts)},
                        ids=np.array([p['entity_id'] for p in parts]))
    print('select: side=%s entities=%d faces=%d conflicts=%d (%.0f s)'
          % (side, len(parts), selection['total_clipped_faces'], census['tissue_tissue']['pairs'],
             selection['wall_seconds']), flush=True)
    return selection, parts


def load_parts(out):
    rows = [json.loads(l) for l in (out / 'members.jsonl').read_text().splitlines()]
    d = np.load(out / 'clipped-members.npz', allow_pickle=False)
    parts = []
    for i, r in enumerate(rows):
        p = dict(r)
        p['V'] = np.ascontiguousarray(d['V%d' % i])
        p['F'] = np.ascontiguousarray(d['F%d' % i])
        parts.append(p)
    return parts


def stage_mesh(out, args, parts):
    priority = priority_table(parts)
    prior = out / 'domain.json'
    results = json.loads(prior.read_text())['resolutions'] if prior.exists() else {}
    for label, setting in (('fine', FINE), ('coarse', COARSE)):
        if args.resolution not in ('both', label):
            continue
        print('mesh[%s]: epsr=%g lr=%g' % (label, setting['epsr'], setting['lr']), flush=True)
        report, mesh = mesh_domain(parts, priority, out / 'ftetwild' / label,
                                   setting['epsr'], setting['lr'], args.workers, label,
                                   full_audit=not args.fast, timeout=args.ftw_timeout)
        results[label] = report
        if mesh is not None:
            np.savez_compressed(out / ('tetmesh-%s.npz' % label), TV=mesh['TV'], TT=mesh['TT'],
                                owner=mesh['owner'])
            print('  tets=%d verts=%d min_dihedral=%.2f deg'
                  % (report['conforming_mesh']['tets'], report['conforming_mesh']['tet_vertices'],
                     report['conforming_mesh']['quality']['min_dihedral_deg']), flush=True)
        write_json(out / 'domain.json', {'schema': 'ihm.hand-reflexive-grip-domain.v1',
                                         'route_decision': 'data/derived/meshing-route-study-v1/route-decision.json',
                                         'mesher_choice': MESHER_CHOICE,
                                         'settings_history': SETTINGS_HISTORY,
                                         'resolutions': results})
    if args.exact_control:
        print('mesh: exact-arrangement control', flush=True)
        results['exact_control'] = exact_route_control(
            parts, priority, 0.005, args.flags,
            'the exact route is recorded to be reproducible only to 48 entities')
        write_json(out / 'domain.json', {'schema': 'ihm.hand-reflexive-grip-domain.v1',
                                         'route_decision': 'data/derived/meshing-route-study-v1/route-decision.json',
                                         'mesher_choice': MESHER_CHOICE,
                                         'settings_history': SETTINGS_HISTORY,
                                         'resolutions': results})
    return results


MESHER_CHOICE = {
    'chosen': 'fTetWild (route B)',
    'why': 'the domain has 109 members. data/derived/meshing-route-study-v1/route-decision.json '
           'records that the exact CGAL arrangement + TetGen -pY route is reproducible only to 48 '
           'entities, produced one negative-volume tet and a 0.000 deg minimum dihedral at that '
           'rung, and failed 6 of 6 controlled repeats at 60 on a byte-identical PLC. fTetWild '
           'reached 100 entities with a worst element near 10 deg and preserves shared nodes.',
    'price': 'fTetWild does not preserve the position of the interface. The route study measured '
             '2.08e-04 m displacement at epsr 1e-4 / lr 0.01 on a thigh box; the number this lane '
             'actually paid on the hand box is in interface_fidelity.max_deviation_m.'}


def stage_actuate(out, args, parts):
    began = time.monotonic()
    materials, _ = material_table()
    result = {'schema': 'ihm.hand-reflexive-grip-actuation.v1'}
    fibre_index = {json.loads(l)['entity_id']: json.loads(l)
                   for l in (FIBRE / 'entities.jsonl').read_text().splitlines()}
    muscles = [p for p in parts if p['role'] == 'muscle']
    tendon_sheaths = [p for p in parts if p['lane'] == 'joint-substrate-candidate-v1'
                      and p['role'] == 'connective_tissue']
    result['actuator_inventory'] = {
        'muscle_entities_in_domain': len(muscles),
        'with_fibre_field': sum(1 for p in muscles if p['entity_id'] in fibre_index),
        'derived_fibre_entities': sum(1 for p in muscles if fibre_index.get(p['entity_id'], {}).get('provenance') == 'derived'),
        'inferred_fibre_entities': sum(1 for p in muscles if fibre_index.get(p['entity_id'], {}).get('provenance') == 'inferred'),
        'without_fibre_field': [p['entity_id'] for p in muscles if p['entity_id'] not in fibre_index],
        'grouped_entities': [{'entity_id': p['entity_id'], 'name': p['name']}
                             for p in muscles if p['name'].lower().startswith('set of')],
        'extrinsic_clipped_by_the_wrist': [
            {'entity_id': p['entity_id'], 'name': p['name'],
             'fraction_of_the_muscle_inside_the_domain': p['surface_volume_m3'] / p['full_volume_m3']
             if p['full_volume_m3'] else None}
            for p in sorted(muscles, key=lambda q: q['entity_id'])
            if p['full_volume_m3'] and p['surface_volume_m3'] < 0.95 * p['full_volume_m3']],
        'intrinsic_whole_in_the_domain': [
            {'entity_id': p['entity_id'], 'name': p['name']}
            for p in sorted(muscles, key=lambda q: q['entity_id'])
            if p['full_volume_m3'] and p['surface_volume_m3'] >= 0.95 * p['full_volume_m3']],
        'tendon_sheath_entities': len(tendon_sheaths)}
    materials_used = {}
    for role in sorted({p['role'] for p in parts}):
        entry = materials.get(role if role != ENVELOPE_ROLE else 'adipose', {})
        lin = entry.get('linear', {})
        materials_used[role] = {
            'young_modulus': lin.get('young_modulus'), 'shear_modulus': lin.get('shear_modulus'),
            'density': lin.get('density'), 'poisson_ratio': lin.get('poisson_ratio'),
            'entities': sum(1 for p in parts if p['role'] == role)}
    result['material_table'] = {
        'source': 'data/derived/tissue-material-candidate-v1/materials.json',
        'per_role': materials_used,
        'roles_with_no_modulus_anywhere': sorted(
            r for r, v in materials_used.items()
            if not (v['young_modulus'] or {}).get('value') and not (v['shear_modulus'] or {}).get('value')),
        'active_stress_pa': ACTIVE_STRESS_PA,
        'active_stress_tier': 'transferred; peak isometric skeletal muscle stress is not in the '
                              'candidate table and is not measured on this specimen',
        'poisson_ratio_assumed': POISSON_ASSUMED,
        'bone_modulus_literature_pa': BONE_E_LITERATURE_PA,
        'bone_modulus_surrogate_pa': BONE_E_SURROGATE_PA,
        'bone_note': 'the candidate table has no bone modulus at any tier. Both a literature '
                     'cortical value and a declared stiff-soft surrogate are carried through the '
                     'stability estimate so the cost of the choice is visible.'}
    for name in (args.dynamics_on,) if args.dynamics_on != 'both' else ('coarse', 'fine'):
        path = out / ('tetmesh-%s.npz' % name)
        if not path.exists():
            result.setdefault('resolutions', {})[name] = {'available': False}
            continue
        d = np.load(path)
        TV, TT, owner = d['TV'], d['TT'], d['owner']
        entry = {'available': True, 'tets': int(len(TT)), 'tet_vertices': int(len(TV)),
                 'structure_tets': int((owner >= 0).sum())}
        direction, fibre_rows = transfer_fibres(parts, TV, TT, owner)
        entry['fibre_transfer'] = {
            'entities': len(fibre_rows),
            'transferred': sum(1 for r in fibre_rows if r['transferred']),
            'derived': sum(1 for r in fibre_rows if r.get('provenance') == 'derived'),
            'inferred': sum(1 for r in fibre_rows if r.get('provenance') == 'inferred'),
            'muscle_tets_with_a_direction': int((np.linalg.norm(direction, axis=1) > 0.5).sum()),
            'max_transfer_distance_m': max((r.get('transfer_distance_max_m', 0.0) for r in fibre_rows), default=0.0),
            'per_entity': fibre_rows,
            'basis': 'nearest source tet centroid in data/derived/muscle-fibre-field-v1; the '
                     'source field lives on each muscle\'s own per-entity TetGen mesh, not on this '
                     'domain, so every direction is a transfer and the distance is reported'}
        for bone_label, bone_e in (('literature_cortical', BONE_E_LITERATURE_PA),
                                   ('stiff_soft_surrogate', BONE_E_SURROGATE_PA)):
            mu, lam, rho, rows = assign_materials(parts, owner, materials, bone_e)
            keep = owner >= 0
            det = np.abs(np.linalg.det(np.swapaxes(TV[TT[keep][:, 1:]] - TV[TT[keep][:, 0, None]], 1, 2)) / 6.0)
            x = TV[TT[keep]]
            areas = np.stack([np.linalg.norm(np.cross(x[:, b] - x[:, a], x[:, c] - x[:, a]), axis=1) / 2
                              for a, b, c in ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))], axis=1)
            height = 3 * det / areas.max(axis=1)
            wave = np.sqrt((lam[keep] + 2 * mu[keep]) / rho[keep])
            dt = float(0.1 * np.min(height / wave))
            entry.setdefault('stability', {})[bone_label] = {
                'bone_young_modulus_pa': bone_e,
                'max_explicit_dt_s': dt,
                'limiting_wave_speed_m_s': float(wave.max()),
                'min_element_height_m': float(height.min()),
                'steps_for_100_ms_contraction': int(np.ceil(0.1 / dt)),
                'materials': rows if bone_label == 'stiff_soft_surrogate' else None}
        if args.run_dynamics and name == args.dynamics_on.replace('both', 'coarse'):
            mu, lam, rho, _ = assign_materials(parts, owner, materials, BONE_E_SURROGATE_PA)
            structure_nodes = np.unique(TT[owner >= 0])
            ymax = float(TV[structure_nodes, 1].max())
            proximal = TV[:, 1] >= ymax - args.grip_clamp_m
            body, keep, used, T, dirk = build_body(TV, TT, owner, mu, lam, rho, direction, proximal)
            entry['solid'] = {'nodes': int(len(body.position_m)), 'tets': int(len(T)),
                              'fixed_nodes': int(len(body.fixed_nodes)),
                              'fixed_node_basis': 'the proximal cut face of the domain: every node '
                                                  'within %g m of the most proximal structure node '
                                                  '(y = %.6f m), which is where the box cut the '
                                                  'forearm' % (args.grip_clamp_m, ymax),
                              'total_mass_kg': float(body.mass_kg.sum()),
                              'max_explicit_dt_s': body.max_explicit_dt_s}
            names = {p['entity_id']: p['name'].lower() for p in parts}
            group = {'flexor': [i for i, p in enumerate(parts) if p['role'] == 'muscle'
                                and ('flexor' in names[p['entity_id']]
                                     or 'lumbrical' in names[p['entity_id']]
                                     or 'interosse' in names[p['entity_id']])],
                     'extensor': [i for i, p in enumerate(parts) if p['role'] == 'muscle'
                                  and 'extensor' in names[p['entity_id']]]}
            owner_kept = owner[keep]
            # A frame-free actuation probe, one force assembly and no time stepping: with the
            # active stress switched on at the reference configuration, what resultant does each
            # distal phalanx feel, and do the flexor and extensor groups oppose each other on it?
            distal = [i for i, p in enumerate(parts)
                      if p['role'] == 'rigid_bone' and 'distal phalanx' in p['name'].lower()]
            resultant = {}
            for gname, idx in group.items():
                sigma = np.zeros(len(dirk))
                sigma[np.isin(owner_kept, idx)] = ACTIVE_STRESS_PA
                f = active_force(body.region, dirk, sigma, body.reference_position_m)
                resultant[gname] = {'actuated_entities': [parts[i]['name'] for i in idx],
                                    'actuated_tets': int(np.isin(owner_kept, idx).sum()),
                                    'total_force_norm_n': float(np.linalg.norm(f, axis=1).sum()),
                                    'net_force_n': [float(v) for v in f.sum(0)],
                                    '_f': f}
            probe = []
            for i in distal:
                nodes = np.unique(T[owner_kept == i])
                if not len(nodes):
                    continue
                row = {'entity_id': parts[i]['entity_id'], 'name': parts[i]['name'],
                       'nodes': int(len(nodes))}
                vec = {}
                for gname in group:
                    r = resultant[gname]['_f'][nodes].sum(0)
                    vec[gname] = r
                    row[gname + '_resultant_n'] = [float(v) for v in r]
                    row[gname + '_resultant_magnitude_n'] = float(np.linalg.norm(r))
                a, b = vec['flexor'], vec['extensor']
                na, nb = np.linalg.norm(a), np.linalg.norm(b)
                row['flexor_extensor_angle_deg'] = float(np.degrees(np.arccos(
                    np.clip(a @ b / (na * nb), -1, 1)))) if na > 0 and nb > 0 else None
                probe.append(row)
            for g in resultant.values():
                g.pop('_f')
            entry['actuation_probe'] = {
                'basis': 'active fibre stress switched on at the reference configuration; the '
                         'nodal force resultant it puts on each distal phalanx, with no time '
                         'integration. An antagonist pair that pulls a fingertip in opposite '
                         'directions is the emergent-articulation signal: nothing declares a '
                         'joint, the direction comes from the tissue and the tendon geometry.',
                'active_stress_pa': ACTIVE_STRESS_PA,
                'groups': resultant, 'distal_phalanges': probe,
                'flexor_extensor_angle_median_deg':
                    float(np.median([r['flexor_extensor_angle_deg'] for r in probe
                                     if r['flexor_extensor_angle_deg'] is not None])) if probe else None,
                'caveat': 'BodyParts3D authors each long flexor and extensor as ONE surface from '
                          'belly to insertion, with no separate tendon entity, so the part of it '
                          'inside this domain is anatomically tendon and is being given a '
                          'contractile stress it does not have. The physically correct load is a '
                          'traction on the proximal cut face; this probe over-reads distal force '
                          'and the fraction of each muscle inside the domain says by how much.'}
            dt = args.dt_fraction * body.max_explicit_dt_s
            entry['passive_settle'] = run_dynamics(body, dirk, None, args.steps, dt,
                                                   (0.0, -9.81, 0.0))
            flex = np.zeros(len(dirk))
            sel = np.isin(owner_kept, group['flexor'])
            flex[sel] = ACTIVE_STRESS_PA
            entry['contraction'] = run_dynamics(body, dirk, flex, args.steps, dt,
                                                (0.0, -9.81, 0.0))
            entry['contraction']['actuated_entities'] = [parts[i]['name'] for i in group['flexor']]
            entry['contraction']['actuated_tets'] = int(sel.sum())
            entry['contraction']['active_stress_pa'] = ACTIVE_STRESS_PA
            per_step = entry['contraction']['wall_seconds_per_step']
            entry['contraction_cost_extrapolation'] = {
                'wall_seconds_per_step': per_step,
                'dt_s': dt,
                'target_contraction_s': 0.1,
                'steps_required': int(np.ceil(0.1 / dt)),
                'wall_hours_required': float(np.ceil(0.1 / dt) * per_step / 3600.0),
                'cpu_share_of_one_core_percent': None,
                'basis': 'explicit symplectic Euler at half the stability limit of '
                         'ihm/assembly/contact_dynamics.py with the stiff-soft bone surrogate, '
                         'measured on this host at this element count. wall_seconds_per_step is '
                         'whatever CPU share this process actually got; on a contended host it '
                         'over-reads. The step COUNT is a property of the discretisation and the '
                         'material and does not move.'}
        result.setdefault('resolutions', {})[name] = entry
    result['wall_seconds'] = time.monotonic() - began
    write_json(out / 'actuation.json', result)
    return result


def stage_audit(out, parts):
    """Recompute the ownership and undisplaced-fidelity blocks from an existing domain.json.

    Both are exact functions of numbers the mesh stage already wrote, so a receipt added later
    costs nothing and never re-runs the mesher."""
    path = out / 'domain.json'
    doc = json.loads(path.read_text())
    env = {p['entity_id'] for p in parts if p['role'] == ENVELOPE_ROLE}
    for label, report in doc['resolutions'].items():
        mesh = report.get('conforming_mesh')
        if not mesh:
            continue
        for r in mesh['per_structure']:
            r['volume_conceded_m3'] = r['claimed_volume_m3'] - r['owned_volume_m3']
            r['fraction_conceded'] = (r['volume_conceded_m3'] / r['claimed_volume_m3']
                                      if r['claimed_volume_m3'] > 0 else None)
        per = {r['entity_id']: r for r in mesh['per_structure']}
        tissue = [r for r in mesh['per_structure'] if r['entity_id'] not in env]
        report['ownership'] = {
            'note': 'the outer body envelope claims every structure tet by construction, so '
                    'conforming_mesh.disputed_tets counts containment as well as conflict. The '
                    'tissue-level numbers below exclude it.',
            'tissue_volume_conceded_m3': float(sum(r['volume_conceded_m3'] for r in tissue)),
            'envelope_volume_conceded_m3': float(sum(per[e]['volume_conceded_m3'] for e in env if e in per)),
            'envelope_owned_volume_m3': float(sum(per[e]['owned_volume_m3'] for e in env if e in per)),
            'entities_fully_displaced': [r['entity_id'] for r in tissue if r['owned_tets'] == 0],
            'worst_conceders': sorted([{k: r[k] for k in ('entity_id', 'name', 'role',
                                                          'surface_volume_m3', 'claimed_volume_m3',
                                                          'owned_volume_m3', 'fraction_conceded')}
                                       for r in tissue if r['fraction_conceded']],
                                      key=lambda r: -(r['fraction_conceded'] or 0))[:15]}
        fid = report.get('interface_fidelity')
        if fid:
            undisplaced = {r['entity_id'] for r in tissue
                           if r['owned_tets'] == r['claimed_tets'] and r['claimed_tets'] > 0}
            clean = [r for r in fid['per_structure']
                     if r.get('faces') and r['entity_id'] in undisplaced]
            fid['undisplaced_reference'] = {
                'entities': len(clean),
                'basis': 'entities whose owned tet set equals their claimed tet set, so no other '
                         'structure displaced them and the deviation is the mesher envelope alone',
                'max_deviation_m': max((r['deviation_max_m'] for r in clean), default=None),
                'median_of_per_entity_max_deviation_m':
                    float(np.median([r['deviation_max_m'] for r in clean])) if clean else None,
                'max_abs_volume_relative_error':
                    max((abs(r['volume_relative_error'] or 0.0) for r in clean), default=None),
                'median_abs_volume_relative_error':
                    float(np.median([abs(r['volume_relative_error'] or 0.0) for r in clean])) if clean else None,
                'worst': sorted([{k: r[k] for k in ('entity_id', 'name', 'role', 'faces',
                                                    'deviation_max_m', 'volume_relative_error')}
                                 for r in clean], key=lambda r: -r['deviation_max_m'])[:10]}
            a = report.setdefault('acceptance', {})
            a['max_undisplaced_deviation_m'] = fid['undisplaced_reference']['max_deviation_m']
            a['max_undisplaced_volume_error'] = fid['undisplaced_reference']['max_abs_volume_relative_error']
        e = np.array([abs(r['claimed_relative_volume_error'] or 0.0) for r in mesh['per_structure']])
        report.setdefault('acceptance', {})['per_structure_volume_error_distribution'] = {
            'entities': int(len(e)), 'median': float(np.median(e)),
            'p90': float(np.percentile(e, 90)), 'max': float(e.max()),
            'entities_over_1e-3': int((e > 1e-3).sum()), 'entities_over_1e-2': int((e > 1e-2).sum()),
            'entities_over_1e-1': int((e > 1e-1).sum())}
    doc['settings_history'] = SETTINGS_HISTORY
    write_json(path, doc)
    return doc


def manifest(out, extra_inputs=()):
    inputs = {'data/derived/canonical/anatomy.json': sha(R.ANATOMY),
              'data/derived/muscle-tet-ready-v1/manifest.json': sha(R.MUSCLE / 'manifest.json'),
              'data/derived/entity-tet-ready-v1/manifest.json': sha(R.ENTITY / 'manifest.json'),
              'data/derived/cross-structure-repair-v1/summary.json': sha(REPAIR / 'summary.json'),
              'data/derived/joint-substrate-candidate-v1/manifest.json': sha(JOINTS / 'manifest.json'),
              'data/derived/muscle-fibre-field-v1/manifest.json': sha(FIBRE / 'manifest.json'),
              'data/derived/tissue-material-candidate-v1/manifest.json': sha(MATERIALS / 'manifest.json'),
              'data/derived/outer-envelope/outer-envelope.json.gz': sha(ENVELOPE),
              'data/derived/meshing-route-study-v1/route-decision.json':
                  sha(ROOT / 'data/derived/meshing-route-study-v1/route-decision.json'),
              'scripts/build_cross_structure_conflict_repair.py': sha(R.__file__),
              'scripts/build_meshing_route_study.py': sha(S.__file__),
              'scripts/build_hand_reflexive_grip_domain.py': sha(__file__),
              'ihm/assembly/contact_dynamics.py': sha(ROOT / 'ihm/assembly/contact_dynamics.py')}
    if S.FTETWILD.exists():
        inputs['data/runtime/tolerant-mesher/build/FloatTetwild_bin'] = sha(S.FTETWILD)
    for k in extra_inputs:
        inputs[k] = sha(ROOT / k)
    artifacts = {p.relative_to(out).as_posix(): sha(p) for p in sorted(out.rglob('*'))
                 if p.is_file() and p.name != 'manifest.json'
                 and p.suffix not in ('.obj', '.msh', '.csv')}
    write_json(out / 'manifest.json',
               {'schema': 'ihm.hand-reflexive-grip-domain-manifest.v1',
                'inputs_sha256': inputs, 'artifacts_sha256': artifacts,
                'canonical_assets_modified': False, 'existing_repo_files_modified': False,
                'python': sys.version})


def run(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    began = time.monotonic()
    if args.stage in ('select', 'all'):
        stage_select(out, args)
    parts = load_parts(out)
    if args.stage in ('mesh', 'all'):
        stage_mesh(out, args, parts)
    if args.stage in ('audit', 'mesh', 'all') and (out / 'domain.json').exists():
        stage_audit(out, parts)
    if args.stage in ('actuate', 'all'):
        stage_actuate(out, args, parts)
    manifest(out)
    print('done in %.0f s' % (time.monotonic() - began), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--out', default=str(ROOT / 'data/derived/hand-reflexive-grip-domain-v1'))
    p.add_argument('--stage', default='all', choices=['select', 'mesh', 'audit', 'actuate', 'all'])
    p.add_argument('--resolution', default='both', choices=['fine', 'coarse', 'both', 'none'])
    p.add_argument('--workers', type=int, default=12)
    p.add_argument('--fast', action='store_true', help='skip the exact pairwise overlap audit')
    p.add_argument('--exact-control', action='store_true',
                   help='also run the exact CGAL arrangement + TetGen -pY route as a control')
    p.add_argument('--flags', default='pY')
    p.add_argument('--ftw-timeout', type=int, default=7200)
    p.add_argument('--run-dynamics', action='store_true')
    p.add_argument('--dynamics-on', default='coarse', choices=['coarse', 'fine', 'both'])
    p.add_argument('--steps', type=int, default=200)
    p.add_argument('--dt-fraction', type=float, default=0.5)
    p.add_argument('--grip-clamp-m', type=float, default=0.002)
    p.add_argument('--self-test', action='store_true')
    a = p.parse_args()
    if a.self_test:
        sys.exit(0 if self_test() else 1)
    run(a)
