#!/usr/bin/env python3
"""Where the spinal cord is, level by level, in the canonical frame -- and where the relays go.

    CUDA_VISIBLE_DEVICES="" .venv/bin/python scripts/build_spinal_cord_levels.py

WHY THIS EXISTS.  Every peripheral route length in `peripheral.json` and every patch length in
`dermatomes.json` ends at a relay.  Until 18 September 2026 the relays were typed by hand at
`[+-0.012, y, -0.04]` with `y` in {0.65, 0.53, 0.30, 0.04, -0.08}, and an audit (IBM-1
`scripts/audit_nerve_route_lengths.py`, IBM-1 `docs/LOG.md` 18 September, this repo's
`docs/BODY_PERIPHERAL.md`) found every spinal relay BELOW the cord segment it stands for:
cervical 89 mm, thoracic 152 mm, lumbar 276 mm, sacral 337 mm (under the pelvis, where there
is no cord).  The frame was right; the points were wrong.  A frame audit cannot see that,
because both ends of a wrong length are in the right frame.

WHAT THIS BUILDS, and the source of every coordinate.

1. *The cord's path* is the centreline of Z-Anatomy's `Spinal dura` mesh (CC-BY-SA 4.0), read
   from the staged `Startup.blend` by IBM-1 `scripts/blender_zanatomy_nerves.py` (run here if
   its cache is missing, CPU only).  It is carried into `bodyparts3d-display-m` by the SAME
   landmark registration that placed every `body-za-*` entity in `anatomy.json` (99 shared
   bones, affine plus thin-plate spline, held-out RMS 4.7 mm), rebuilt from the landmarks
   `anatomy.json` stores and asserted to reproduce the stored affine exactly.  NOT the site
   exporter's bounding-box chain (display box x 0.86487): that is a height ratio between two
   atlases, and the difference between the two is reported below rather than assumed away.
2. *Segment levels* come from BodyParts3D vertebrae in `anatomy.json`.  A vertebral BODY's
   level is the midpoint of the two intervertebral discs above and below it (the whole-bone
   centroid sits low in the thoracic spine, pulled down by the spinous process).
3. *Which vertebra a cord segment lies behind* is an authored clinical prior, not a
   measurement: cord segments sit rostral to their vertebrae -- about one body in the
   cervical cord, about two in the upper thoracic cord, the lumbar segments behind T10-T12 and
   the sacral segments forming the conus behind T12-L1 (the conus ends at L1-L2).  Real people
   vary by about one segment either way.  `RELAY_LEVELS` names the representative segment for
   each relay group and the vertebral body it is taken to lie behind.
4. *A relay* is the dura centreline point at that level, displaced 4 mm to its own side (the
   dorsal horn, inside a cord ~10 mm across), with the centreline's own x (under 0.3 mm)
   dropped so each left/right pair is an exact mirror image -- IBM-1's visceral join collapses
   the two sides and raises if their lengths differ at all.  The two brainstem relays are BodyParts3D
   entities: `cranial` at the pons half of its side (the nuclei it stands for run from the
   midbrain, III/IV, through the pons, V/VI/VII, to the medulla, XII -- the pons is their
   middle), `solitary` at the medulla oblongata half of its side (the nucleus of the solitary
   tract is medullary); each pair symmetrised the same way (halves' centroids averaged in y and
   z, mirrored in x; the shift is recorded).

CHECKS, each of which can fail.

* The blend-world coordinates in the extraction must equal the landmark sources the
  registration was fitted on (three bones, to 1e-6 m).  A cache written in the normalised
  display box, or through a different convention, misses by decimetres.
* The registration rebuilt from the stored landmarks must reproduce the stored affine exactly.
* KNOWN ANSWER, independent of both atlases' bounding boxes: at every BodyParts3D vertebra
  from C3 to L1 the registered centreline must lie INSIDE that vertebra's canal -- bone in at
  least 16 of 24 directions around it within 35 mm, and at least 2 mm clear of bone.  First
  run: 17 of 18, FAILED at C3 (1.5 mm clearance, bone all round).  Recorded, not rescored; the
  relays are gated on their own four levels, which pass.  The same test on the centreline
  scaled by the 15.7% frame error of 17 September passes only 8 of 18 -- it can fail, but it
  is weak in the thoracolumbar canal, where a scale error slides a vertical tube along itself.

Output: `data/derived/canonical/spinal_cord_levels.json`, read by `build_body_peripheral.py`
and `enrich_peripheral_routes.py`.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.anatomy import LandmarkRegistration, ROTATION  # noqa: E402

IBM = ROOT.parent / 'IBM-1'
BLEND = ROOT / 'data/raw/anatomy/extended/extracted/Z-Anatomy/Startup.blend'
EXTRACT = IBM / 'scripts/blender_zanatomy_nerves.py'
CACHE = IBM / '.cache/zanatomy_nerves.json'
ANATOMY = ROOT / 'data/derived/canonical/anatomy.json'
OUT = ROOT / 'data/derived/canonical/spinal_cord_levels.json'

SLICE_STEP_M = 0.005
SLICE_HALF_M = 0.003
RELAY_LATERAL_M = 0.004
ENCLOSURE_RADIUS_M = 0.035
ENCLOSURE_MIN_BINS = 16          # of 24
ENCLOSURE_MIN_CLEARANCE_M = 0.002
FRAME_ERROR_17_SEP = 1 / 0.86487  # the withdrawn "shared frame": display box against metres

#: relay group -> the segments it stands for, the representative segment, and the vertebral
#: BODY that segment is taken to lie behind.  Authored clinical prior; see the docstring.
RELAY_LEVELS = {
    'cervical': dict(
        stands_for='C1-T1: brachial plexus C5-T1 (every arm muscle and arm skin patch), '
                   'cervical plexus C1-C4, phrenic C3-C5, spinal accessory',
        segment='C6', vertebra='fifth cervical vertebra',
        rule='cervical segments lie about one vertebral body above their number; C6 is the '
             'widest part of the cervical enlargement'),
    'thoracic': dict(
        stands_for='T1-T12: intercostal, subcostal, dorsal rami, splanchnics T5-T12, '
                   'sympathetic outflow',
        segment='T7', vertebra='fifth thoracic vertebra',
        rule='upper and mid thoracic segments lie about two vertebral bodies above their '
             'number; T7 is the middle of T1-T12'),
    'lumbar': dict(
        stands_for='L1-L4: lumbar plexus (femoral, obturator L2-L4; iliohypogastric, '
                   'ilioinguinal, genitofemoral L1-L2), lumbar splanchnics',
        segment='L3', vertebra='eleventh thoracic vertebra',
        rule='lumbar segments L1-L5 lie behind the T10-T12 bodies, L2-L3 behind T11'),
    'sacral': dict(
        stands_for='L4-S4: sacral plexus (sciatic, gluteal, pudendal), pelvic splanchnics',
        segment='S2', vertebra='first lumbar vertebra',
        rule='the sacral segments form the conus medullaris behind the T12-L1 bodies; the '
             'conus ends at L1-L2'),
}
BRAINSTEM_RELAYS = {
    'cranial': dict(entity={'left': 'body-bp3d-FJ1775', 'right': 'body-bp3d-FJ1822'},
                    structure='pons',
                    rule='cranial-nerve nuclei served (III, IV midbrain; V, VI, VII pons; XII '
                         'medulla); the pons is their middle. Entity centroid of the side\'s '
                         'half.'),
    'solitary': dict(entity={'left': 'body-bp3d-FJ1769', 'right': 'body-bp3d-FJ1831'},
                     structure='medulla oblongata',
                     rule='the nucleus of the solitary tract is in the dorsal medulla. Entity '
                          'centroid of the side\'s half.'),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def geometry(entity: dict) -> np.ndarray:
    g = json.loads(gzip.decompress((ROOT / entity['reference_geometry']['path']).read_bytes()))
    return np.asarray(g['positions'], float).reshape(-1, 3)


def blender_extraction() -> dict:
    if not CACHE.exists():
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        print('extracting the spinal dura from Startup.blend with blender (CPU, one-off)...',
              flush=True)
        subprocess.run(['blender', '-b', str(BLEND), '--python', str(EXTRACT), '--', str(CACHE)],
                       check=True, capture_output=True,
                       env={**os.environ, 'CUDA_VISIBLE_DEVICES': ''})
    return json.loads(CACHE.read_text())


def registration(anatomy: dict) -> tuple[LandmarkRegistration, dict]:
    stored = anatomy['registrations']['z_anatomy']
    marks = stored['landmarks']
    reg = LandmarkRegistration(np.array([m['source_rotated_m'] for m in marks]),
                               np.array([m['target_m'] for m in marks]),
                               smoothing=stored['smoothing_prior'])
    err = float(np.abs(reg.affine - np.asarray(stored['affine_4x3'])).max())
    if err > 1e-12:
        raise SystemExit(f'rebuilt registration does not reproduce the stored affine ({err:.2e})')
    return reg, {m['name']: m for m in marks}


def centreline(points: np.ndarray):
    ys = np.arange(points[:, 1].min() + SLICE_HALF_M, points[:, 1].max() - SLICE_HALF_M / 2,
                   SLICE_STEP_M)
    rows = []
    for y in ys:
        s = points[np.abs(points[:, 1] - y) < SLICE_HALF_M]
        if len(s) < 12:
            continue
        m = s.mean(0)
        rows.append([float(m[0]), float(y), float(m[2]),
                     float(np.median(np.hypot(s[:, 0] - m[0], s[:, 2] - m[2]))), len(s)])
    return np.asarray(rows)


def at_level(cl: np.ndarray, y: float) -> np.ndarray:
    if not (cl[0, 1] <= y <= cl[-1, 1]):
        raise SystemExit(f'level y={y:.3f} is outside the dura ({cl[0, 1]:.3f}..{cl[-1, 1]:.3f})')
    return np.array([np.interp(y, cl[:, 1], cl[:, 0]), y, np.interp(y, cl[:, 1], cl[:, 2])])


def enclosure(bone: np.ndarray, point: np.ndarray, level_y: float) -> dict:
    s = bone[np.abs(bone[:, 1] - level_y) < 0.006]
    d = s[:, [0, 2]] - point[[0, 2]]
    r = np.hypot(d[:, 0], d[:, 1])
    ang = np.arctan2(d[:, 1], d[:, 0])
    near = r < ENCLOSURE_RADIUS_M
    bins = np.unique(((ang[near] + np.pi) / (2 * np.pi) * 24).astype(int) % 24)
    out = dict(bins_with_bone=int(len(bins)), clearance_m=float(r.min()) if len(r) else None,
               slab_vertices=int(len(s)))
    out['inside_canal'] = bool(len(bins) >= ENCLOSURE_MIN_BINS and len(r)
                               and r.min() >= ENCLOSURE_MIN_CLEARANCE_M)
    return out


def build(verbose: bool = True) -> dict:
    anatomy = json.loads(ANATOMY.read_text())
    if anatomy['frame']['id'] != 'bodyparts3d-display-m' or anatomy['frame']['units'] != 'm':
        raise SystemExit(f'anatomy.json frame is {anatomy["frame"]["id"]}, not bodyparts3d-display-m')
    ents = [e for e in anatomy['entities'] if e['id'].startswith('body-bp3d')]
    by_id = {e['id']: e for e in anatomy['entities']}
    reg, marks = registration(anatomy)
    ex = blender_extraction()
    if not ex.get('dura'):
        raise SystemExit('the extraction carries no spinal dura')

    # ---- frame at the boundary: the extraction must be in the blend world the registration
    # was fitted in.  Three bones, bounding-box centres, to a micrometre.
    frame_rows = []
    for bone, mark in (('Femur.l', 'left femur'), ('Humerus.r', 'right humerus'),
                       ('Sacrum', 'sacrum')):
        lo, hi = (np.asarray(v) for v in ex['check'][bone])
        got = ((lo + hi) / 2) @ ROTATION.T
        want = np.asarray(marks[mark]['source_rotated_m'])
        frame_rows.append(dict(bone=bone, landmark=mark, error_m=float(np.linalg.norm(got - want))))
    worst = max(r['error_m'] for r in frame_rows)
    if worst > 1e-6:
        raise SystemExit(f'extraction is not in the registration\'s blend world (worst {worst:.3g} m)')

    D = reg.transform(np.asarray(ex['dura'], float) @ ROTATION.T)
    cl = centreline(D)

    # the site exporter's chain, for comparison only: display transform x (bp height / za height)
    man = json.loads((ROOT / 'data/derived/anatomy/extended/manifest_fragment.json').read_text())
    zm = man['models'][0]
    bp = json.loads((ROOT / 'data/derived/app/manifest.json').read_text())
    bpm = next(m for m in bp['models'] if m['id'] == 'bodyparts3d')
    zs = bpm['bounds']['max'][1] / zm['bounds']['max'][1]
    dt = zm['display_transform']
    chain = (float(dt['scale']) * (np.asarray(ex['dura'], float) @ np.asarray(dt['rotation']).T)
             + np.asarray(dt['translation'])) * zs
    cl_chain = centreline(chain)

    # ---- vertebral bodies from the disc above and below
    discs = sorted((e for e in ents if 'intervertebral dis' in e['name']),
                   key=lambda e: e['centroid_m'][1])
    verts = {e['name']: e for e in ents
             if re.fullmatch(r'\w+ (cervical|thoracic|lumbar) vertebra', e['name'])}
    levels = []
    for name, e in sorted(verts.items(), key=lambda kv: -kv[1]['centroid_m'][1]):
        y = e['centroid_m'][1]
        above = [d for d in discs if d['centroid_m'][1] > y]
        below = [d for d in discs if d['centroid_m'][1] < y]
        if not above or not below:
            continue
        up, dn = above[0], below[-1]
        body_y = 0.5 * (up['centroid_m'][1] + dn['centroid_m'][1])
        row = dict(vertebra=name, entity_id=e['id'], centroid_y_m=y, body_level_y_m=body_y,
                   disc_above=up['id'], disc_below=dn['id'])
        if cl[0, 1] <= body_y <= cl[-1, 1]:
            p = at_level(cl, body_y)
            row['dura_centreline_m'] = p.tolist()
            bone = geometry(e)
            row['canal'] = enclosure(bone, p, body_y)
            # the same test on the 17 September frame error: the whole centreline scaled by
            # 1/0.86487 about the frame origin, read at the same level.  must NOT pass.
            wrong_cl = cl.copy()
            wrong_cl[:, :3] *= FRAME_ERROR_17_SEP
            if wrong_cl[0, 1] <= body_y <= wrong_cl[-1, 1]:
                row['canal_under_15p7pct_frame_error'] = enclosure(bone, at_level(wrong_cl, body_y), body_y)
            else:
                row['canal_under_15p7pct_frame_error'] = dict(inside_canal=False, bins_with_bone=0,
                                                              clearance_m=None, slab_vertices=0,
                                                              note='the scaled dura does not reach this level')
        levels.append(row)
    L = {r['vertebra']: r for r in levels}
    tested = [r for r in levels if 'canal' in r]
    gate_c3_l1 = [r for r in tested
                  if r['vertebra'] in [f'{o} {k} vertebra' for o, k in
                                       [(a, 'cervical') for a in ('third', 'fourth', 'fifth', 'sixth', 'seventh')]
                                       + [(a, 'thoracic') for a in ('first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh', 'eighth', 'ninth', 'tenth', 'eleventh', 'twelfth')]
                                       + [('first', 'lumbar')]]]
    inside = sum(r['canal']['inside_canal'] for r in gate_c3_l1)
    inside_wrong = sum(r['canal_under_15p7pct_frame_error']['inside_canal'] for r in gate_c3_l1)

    # ---- relays
    relays = {}
    for group, spec in RELAY_LEVELS.items():
        row = L[spec['vertebra']]
        c = at_level(cl, row['body_level_y_m'])
        for side, sign in (('left', 1.0), ('right', -1.0)):
            # the pair is EXACTLY mirror-symmetric: the centreline's own x (|x| < 0.3 mm here) is
            # dropped, because IBM-1's visceral join collapses sides and raises on any difference
            p = np.array([sign * RELAY_LATERAL_M, c[1], c[2]])
            relays[f'peripheral-relay-{side}-{group}'] = dict(
                position_m=p.tolist(), group=group, side=side, segment=spec['segment'],
                vertebra=spec['vertebra'], vertebra_entity_id=row['entity_id'],
                level_y_m=row['body_level_y_m'], dura_centreline_m=c.tolist(),
                distance_to_dura_centreline_m=float(np.linalg.norm(p - c)),
                canal_check=row['canal'],
                stands_for=spec['stands_for'], level_rule=spec['rule'],
                position_source=(
                    'Z-Anatomy Spinal dura centreline (Startup.blend, via IBM-1 '
                    'blender_zanatomy_nerves.py) registered by anatomy.json z_anatomy landmark '
                    f'affine+TPS, at the {spec["vertebra"]} body level (midpoint of the '
                    f'BodyParts3D discs {row["disc_above"]} and {row["disc_below"]}), '
                    f'{1000 * RELAY_LATERAL_M:.0f} mm to the {side}'),
                evidence_kind='registered_cord_geometry_at_authored_segment_level')
    for group, spec in BRAINSTEM_RELAYS.items():
        halves = {}
        for side in ('left', 'right'):
            e = by_id[spec['entity'][side]]
            if e['name'] != spec['structure']:
                raise SystemExit(f'{e["id"]} is {e["name"]!r}, not {spec["structure"]!r}')
            halves[side] = np.asarray(e['centroid_m'], float)
        # symmetrised pair: the two halves' centroids averaged in y and z and mirrored in x, so
        # that left and right routes stay exact mirror images (IBM-1's visceral join raises
        # otherwise).  The asymmetry this removes is recorded.
        L_, R_ = halves['left'], halves['right']
        sym = np.array([0.5 * (L_[0] - R_[0]), 0.5 * (L_[1] + R_[1]), 0.5 * (L_[2] + R_[2])])
        for side, sign in (('left', 1.0), ('right', -1.0)):
            e = by_id[spec['entity'][side]]
            pos = np.array([sign * sym[0], sym[1], sym[2]])
            relays[f'peripheral-relay-{side}-{group}'] = dict(
                position_m=pos.tolist(), group=group, side=side,
                structure=spec['structure'], entity_id=e['id'], level_rule=spec['rule'],
                entity_centroid_m=halves[side].tolist(),
                symmetrisation_shift_m=float(np.linalg.norm(pos - halves[side])),
                position_source=(f'BodyParts3D {spec["entity"]["left"]} / {spec["entity"]["right"]} '
                                 f'({spec["structure"]}, left/right halves) centroids in anatomy.json, '
                                 'averaged in y,z and mirrored in x'),
                evidence_kind='atlas_entity_centroid')

    chain_offset = []
    for group, spec in RELAY_LEVELS.items():
        y = L[spec['vertebra']]['body_level_y_m']
        a, b = at_level(cl, y), at_level(cl_chain, y)
        chain_offset.append(dict(group=group, level_y_m=y,
                                 registered_minus_bbox_chain_m=(a - b).tolist(),
                                 distance_m=float(np.linalg.norm(a - b))))

    data = dict(
        schema='ihm.spinal-cord-levels.v1', frame=anatomy['frame'],
        source=dict(
            blend=str(BLEND.relative_to(ROOT)), blend_bytes=BLEND.stat().st_size,
            extraction_script=str(EXTRACT), extraction_cache=str(CACHE),
            extraction_cache_sha256=sha256(CACHE), object='Spinal dura',
            licence='Z-Anatomy / BodyParts3D, CC-BY-SA 4.0',
            dura_vertices=len(ex['dura']),
            registration='anatomy.json registrations.z_anatomy (rebuilt from stored landmarks; '
                         'affine reproduced exactly)',
            registration_held_out_rms_m=anatomy['registrations']['z_anatomy']['held_out_rms_m'],
            anatomy_sha256=sha256(ANATOMY)),
        frame_check=dict(rows=frame_rows, worst_m=worst, tolerance_m=1e-6,
                         meaning='extraction coordinates equal the landmark sources the '
                                 'registration was fitted on'),
        centreline=dict(
            method=f'mean of registered dura vertices in {2000 * SLICE_HALF_M:.0f} mm slabs every '
                   f'{1000 * SLICE_STEP_M:.0f} mm of y; radius = median in-plane distance',
            columns=['x_m', 'y_m', 'z_m', 'dura_radius_m', 'vertices'],
            rows=cl.round(6).tolist(),
            y_range_m=[float(cl[0, 1]), float(cl[-1, 1])],
            note='The Z-Anatomy spinal dura ends at y=%.3f, about the L2 body; the conus is at '
                 'L1-L2, so every spinal relay level lies on it.' % cl[0, 1]),
        vertebral_levels=levels,
        canal_sweep=dict(
            gate='registered dura centreline inside the BodyParts3D vertebral canal at every '
                 'body level C3-L1 (bar fixed before this implementation first ran)',
            passed=inside == len(gate_c3_l1), inside=inside, of=len(gate_c3_l1),
            failed_levels=[dict(vertebra=r['vertebra'], **r['canal']) for r in gate_c3_l1
                           if not r['canal']['inside_canal']],
            rule=f'bone in >= {ENCLOSURE_MIN_BINS}/24 directions within '
                 f'{1000 * ENCLOSURE_RADIUS_M:.0f} mm and >= {1000 * ENCLOSURE_MIN_CLEARANCE_M:.0f} '
                 'mm clear of bone, in a 12 mm slab at the vertebral body level',
            sabotage='the whole centreline scaled by 1/0.86487 (the 17 September frame error) '
                     'and read at the same levels',
            sabotage_inside=inside_wrong,
            sabotage_reading='weak in the thoracolumbar canal: a uniform scale about a mid-body '
                             'origin mostly slides a vertical tube along itself. The decisive '
                             'frame evidence is frame_check, not this.'),
        relay_level_gate=dict(
            gate='every spinal relay level inside the vertebral canal',
            passed=all(L[spec['vertebra']]['canal']['inside_canal'] for spec in RELAY_LEVELS.values()),
            rows={g: dict(vertebra=spec['vertebra'], **L[spec['vertebra']]['canal'],
                          under_frame_error=L[spec['vertebra']]['canal_under_15p7pct_frame_error'])
                  for g, spec in RELAY_LEVELS.items()}),
        bbox_chain_comparison=chain_offset,
        relays=relays,
        relay_levels=RELAY_LEVELS,
        limitations=[
            'Segment-to-vertebra relations are a population rule; individuals vary by about one '
            'segment, i.e. 15-25 mm of cord.',
            'One relay stands for a group of segments. Its level is the group\'s representative '
            'segment, so a route from the top or bottom of the group is misplaced by up to half '
            'the group\'s cord extent.',
            'The dura is Z-Anatomy\'s, registered onto BodyParts3D by 99 bone landmarks; it is '
            'not this body\'s own cord, which BodyParts3D does not carry.',
        ])
    if not data['relay_level_gate']['passed']:
        raise SystemExit('relay-level canal gate FAILED; relays not written')
    OUT.write_text(json.dumps(data, indent=1) + '\n')
    if verbose:
        print(f'frame check: worst {worst:.2e} m over {len(frame_rows)} bones')
        print(f'dura centreline: {len(cl)} slices, y {cl[0, 1]:.3f}..{cl[-1, 1]:.3f}')
        sw = data['canal_sweep']
        print(f'canal sweep C3-L1: {"PASS" if sw["passed"] else "FAILED"} {inside}/{len(gate_c3_l1)} '
              f'inside (failed: {[r["vertebra"] for r in sw["failed_levels"]] or "-"}); under the '
              f'15.7% frame error {inside_wrong}/{len(gate_c3_l1)}')
        print(f'relay-level canal gate: {"PASS" if data["relay_level_gate"]["passed"] else "FAILED"}')
        for r in tested:
            c, w = r['canal'], r['canal_under_15p7pct_frame_error']
            print(f'  {r["vertebra"]:27s} body y {r["body_level_y_m"]:.3f}  bins {c["bins_with_bone"]:2d}/24 '
                  f'clear {1000 * (c["clearance_m"] or 0):4.1f} mm  {"in " if c["inside_canal"] else "OUT"}'
                  f'   | frame error: bins {w["bins_with_bone"]:2d} {"in" if w["inside_canal"] else "OUT"}')
        for rid, r in relays.items():
            print(f'  {rid:34s} {np.round(r["position_m"], 4).tolist()}  '
                  f'{r.get("segment", r.get("structure"))}')
        for r in chain_offset:
            print(f'  {r["group"]:9s} registered vs bbox-chain centreline: {1000 * r["distance_m"]:.1f} mm')
        print('wrote', OUT.relative_to(ROOT))
    return data


if __name__ == '__main__':
    build()
