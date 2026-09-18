#!/usr/bin/env python3
"""Known answers for the runtime anatomy pose (ihm/assembly/anatomy_pose.py).

    CUDA_VISIBLE_DEVICES="" .venv/bin/python scripts/verify_anatomy_pose.py

Writes `data/derived/anatomy-segment-binding/runtime_pose_report.json` and exits non-zero if a
gate fails.  Every gate here can fail; the ones that test a guard are run on input built to
trip it, because a guard nobody has watched fire is not a guard.

1. FK against Simbody: the pure-python forward kinematics must reproduce `transform_ground`
   in stored native frames (every body, every frame found).
2. Rest pose: at the binding's registered reference pose every entity is where the atlas put it
   (max displacement < 1e-9 m).  The model's DEFAULT pose is also reported, and it is NOT the
   atlas pose -- the atlas was registered at a fitted pose, not at the .osim defaults.
3. Distal only: each of the 31 independent coordinates is driven through its declared range
   from the reference pose, alone.  Exactly the entities on segments kinematically downstream
   of that coordinate may move, and all of them must.  Both pivot modes.
4. Idempotence: the same state twice, and through a second poser, gives bit-identical output.
5. Frame guards: millimetres, the atlas frame, the 15.7% rescale, a missing body, and vertices
   declared in the normalised display frame must each be REFUSED.
6. Coverage: every anatomy.json entity is posed rigidly, blended with the skin, or listed with
   its reason.
7. Joint opening, measured rather than assumed: the 64 closest cross-joint bone point pairs at
   rest, and how far their separation changes through each joint's range, for OpenSim pivots
   against anatomical (closest-bone-surface) pivots.  Reported; not a gate.
8. A real motion: every frame of `data/derived/gait-best/trajectory.json` and the stored native
   frames, posed; finite and rigid; time per call.
"""
from __future__ import annotations

import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.anatomy_pose import (AnatomyPoser, FrameError, FRAME, SKIN_ID)  # noqa: E402

OUT = ROOT / 'data/derived/anatomy-segment-binding/runtime_pose_report.json'
NATIVE = ['data/derived/supine-equilibrium-initial-xfhqpo4e/advanced_1us.json',
          'data/derived/native-stream-smoke-n6e0pvqi/supine/smoke.json',
          'data/derived/opensim-instance-mass-7lmt4max/baseline_frames.json',
          'data/derived/opensim-instance-mass-7lmt4max/variant_frames.json',
          'data/derived/native-surface-foundation-cats536k/surface_initial.json',
          'data/derived/supine-supported-tones-db6xdj5m/initialized.json']
GAIT = 'data/derived/gait-best/trajectory.json'
MOVE_TOL = 1e-9


def native_frames(x):
    if isinstance(x, dict):
        if 'bodies' in x and 'coordinates' in x:
            yield x
        for v in x.values():
            yield from native_frames(v)
    elif isinstance(x, list):
        for v in x:
            yield from native_frames(v)


def main() -> int:
    t0 = time.time()
    P = AnatomyPoser.from_workspace(ROOT)
    Pa = AnatomyPoser.from_workspace(ROOT, pivot='anatomical')
    rec: dict = {'frame': FRAME, 'load_s': time.time() - t0, 'gates': []}
    kin = P.kin

    def gate(name, passed, **extra):
        rec['gates'].append(dict(gate=name, passed=bool(passed), **extra))
        print(f'{"PASS" if passed else "FAIL"}  {name}  ' +
              '  '.join(f'{k}={v}' for k, v in extra.items() if not isinstance(v, (list, dict))))

    # ---- 1. FK against Simbody
    worst, per_body, nframes, files = 0.0, {}, 0, []
    for f in NATIVE:
        p = ROOT / f
        if not p.exists():
            continue
        files.append(f)
        for fr in native_frames(json.loads(p.read_text())):
            q = {k: (v['value'] if isinstance(v, dict) else v) for k, v in fr['coordinates'].items()}
            T = kin.forward(kin.complete(q))
            for b, r in fr['bodies'].items():
                e = float(np.abs(np.asarray(r['transform_ground']) - T[b]).max())
                per_body[b] = max(per_body.get(b, 0.0), e)
                worst = max(worst, e)
            nframes += 1
    gate('FK reproduces Simbody transform_ground', worst < 1e-6, frames=nframes, files=len(files),
         worst_abs=f'{worst:.2e}', worst_body=max(per_body, key=per_body.get))
    rec['fk_per_body_worst'] = per_body

    # ---- 2. rest pose
    ref = P.pose_from_coordinates({})
    d = np.linalg.norm(ref.centroid_m - P.rest_centroid, axis=1)
    rot = np.abs(ref.rotation - np.eye(3)).max()
    gate('reference pose reproduces the atlas rest pose', d.max() < 1e-9 and rot < 1e-9,
         max_m=f'{d.max():.2e}', median_m=f'{np.median(d):.2e}', max_rotation_dev=f'{rot:.1e}')
    # also every VERTEX of a spread of entities, not only centroids
    ents = {e['id']: e for e in json.loads((ROOT / 'data/derived/canonical/anatomy.json').read_text())['entities']}
    sample = P.entity_ids[::200]
    vmax = 0.0
    for eid in sample:
        with gzip.open(ROOT / ents[eid]['reference_geometry']['path']) as f:
            V = np.asarray(json.load(f)['positions'], float).reshape(-1, 3)
        vmax = max(vmax, float(np.abs(ref.transform_vertices(eid, V) - V).max()))
    gate('reference pose: vertices of 20 entities unmoved', vmax < 1e-9, entities=len(sample), max_m=f'{vmax:.2e}')
    dflt = P.pose_from_coordinates({}, fill='default')
    dd = np.linalg.norm(dflt.centroid_m - P.rest_centroid, axis=1)
    rec['model_default_pose'] = dict(
        max_m=float(dd.max()), median_m=float(np.median(dd)),
        reading='the .osim default pose is not the pose the atlas was registered at (binding '
                'reference_pose_rad: pelvis_ty 1.0185 vs 0.93 default, lumbar_rotation 0.236 rad, '
                'hip_rotation -0.14..-0.18 rad, ...); this is the anatomy moved into the model '
                'default, not an error')
    print(f'      model DEFAULT pose (not the atlas pose): max {1000 * dd.max():.1f} mm, median '
          f'{1000 * np.median(dd):.1f} mm from rest')

    # ---- 3. distal only, both pivot modes
    coupled = {c['independent'][0]: c['dependent'] for c in kin.couplers}
    dependents = set(coupled.values())
    rows = []
    for mode, poser in (('opensim', P), ('anatomical', Pa)):
        base = poser.pose_from_coordinates({})
        bad = 0
        for cname, spec in kin.coordinates.items():
            if cname in dependents:
                continue
            joints = [j for j in kin.order if cname in j['coords']
                      or (cname in coupled and coupled[cname] in j['coords'])]
            affected = set().union(*(kin.descendants(j['child']) for j in joints))
            expect = np.isin(np.array(poser.segments)[poser.segment_of], sorted(affected))
            lo, hi = spec['range']
            if not np.isfinite(lo) or hi - lo > 2 * np.pi:
                lo, hi = -np.pi / 2, np.pi / 2           # the shoulders declare [-10, 10] rad
            ref_v = poser.reference_pose[cname]
            moved = np.zeros(len(poser.entity_ids), bool)
            for v in np.linspace(lo, hi, 7):
                if abs(v - ref_v) < 1e-6:
                    continue
                p = poser.pose_from_coordinates({cname: float(v)})
                moved |= (np.linalg.norm(p.centroid_m - base.centroid_m, axis=1) > MOVE_TOL) | \
                         (np.abs(p.rotation - base.rotation).max(axis=(1, 2)) > MOVE_TOL)
            wrong = int((moved & ~expect).sum())
            still = int((~moved & expect).sum())
            bad += wrong + still
            rows.append(dict(mode=mode, coordinate=cname, segments=sorted(affected),
                             expected=int(expect.sum()), moved=int(moved.sum()),
                             moved_but_not_distal=wrong, distal_but_unmoved=still))
        gate(f'single coordinate moves exactly its distal entities ({mode} pivots)', bad == 0,
             coordinates=sum(r['mode'] == mode for r in rows), violations=bad)
    rec['distal_only'] = rows
    for r in rows:
        if r['mode'] == 'opensim' and r['coordinate'] in ('knee_angle_r', 'elbow_flex_l', 'ankle_angle_r',
                                                          'hip_flexion_l', 'lumbar_extension', 'pelvis_tilt'):
            print(f'      {r["coordinate"]:18s} moved {r["moved"]:4d}, expected {r["expected"]:4d} '
                  f'({", ".join(r["segments"]) if len(r["segments"]) < 6 else str(len(r["segments"])) + " segments"})')

    # ---- 4. idempotence
    q = {'knee_angle_r': 0.9, 'hip_flexion_l': 0.4, 'elbow_flex_r': 1.1, 'lumbar_bending': 0.1}
    a, b = P.pose_from_coordinates(q), P.pose_from_coordinates(q)
    P2 = AnatomyPoser.from_workspace(ROOT)
    c = P2.pose_from_coordinates(q)
    same = all(np.array_equal(getattr(a, k), getattr(x, k)) for x in (b, c)
               for k in ('rotation', 'translation', 'centroid_m', 'segment_motion'))
    sa, sb = P.skin_vertices(a), P.skin_vertices(b)
    aa, ab = Pa.pose_from_coordinates(q), Pa.pose_from_coordinates(q)
    same_a = all(np.array_equal(getattr(aa, k), getattr(ab, k)) for k in ('rotation', 'translation'))
    gate('idempotent: same state twice, and via a fresh poser, bit-identical',
         same and np.array_equal(sa, sb) and same_a and a.entity_ids == c.entity_ids)

    # ---- 5. frame guards, each driven to fire
    T = kin.forward(kin.complete({}, P.reference_pose))
    trials = {
        'millimetres': {k: np.diag([1, 1, 1, 1]) @ np.block([[v[:3, :3], 1000 * v[:3, 3:]], [v[3:]]]) for k, v in T.items()},
        'atlas frame (A applied)': {k: P.A @ v for k, v in T.items()},
        'atlas frame, rotation part only rigid': {k: np.block([[v[:3, :3], P.scale * v[:3, 3:]], [v[3:]]]) for k, v in T.items()},
        '15.7% rescale (display box)': {k: np.block([[v[:3, :3], v[:3, 3:] / 0.86487], [v[3:]]]) for k, v in T.items()},
        'missing body': {k: v for k, v in T.items() if k != 'hand_l'},
    }
    fired = {}
    for label, bodies in trials.items():
        try:
            P.pose(bodies)
            fired[label] = 'NOT REFUSED'
        except (FrameError, ValueError) as e:
            fired[label] = f'refused: {str(e)[:90]}'
    try:
        ref.transform_vertices(P.entity_ids[0], np.zeros((1, 3)), frame='z-anatomy-display-normalized')
        fired['vertices in the normalised display frame'] = 'NOT REFUSED'
    except FrameError as e:
        fired['vertices in the normalised display frame'] = f'refused: {str(e)[:90]}'
    try:
        ref.transform_vertices(SKIN_ID, np.zeros((1, 3)))
        fired['rigid pose asked of the skin'] = 'NOT REFUSED'
    except KeyError as e:
        fired['rigid pose asked of the skin'] = f'refused: {str(e)[:90]}'
    ok_native = True
    try:
        P.pose(T)
    except FrameError:
        ok_native = False
    gate('frame guards refuse every wrong input and accept the right one',
         ok_native and all(v.startswith('refused') for v in fired.values()), trials=len(fired))
    for k, v in fired.items():
        print(f'      {k:42s} {v}')
    rec['frame_guards'] = fired

    # ---- 6. coverage
    cov = dict(P.coverage)
    cov['unbound_reasons'] = {k: P.unbound[k] for k in sorted(P.unbound)}
    cov['symmetry_overrides'] = P.symmetry_overrides
    rec['coverage'] = cov
    accounted = cov['rigid'] + cov['skin_blend'] + cov['skin_layers_follow'] + len(cov['not_posed'])
    gate('every anatomy.json entity posed or listed with a reason', accounted == cov['anatomy_entities'],
         entities=cov['anatomy_entities'], rigid=cov['rigid'], skin=cov['skin_blend'],
         layers=cov['skin_layers_follow'], not_posed=len(cov['not_posed']),
         symmetry_overrides=len(P.symmetry_overrides))
    skin = P.skin_vertices(ref)
    rest = np.asarray(json.loads(gzip.decompress((ROOT / ents[SKIN_ID]['reference_geometry']['path']).read_bytes()))['positions']).reshape(-1, 3)
    gate('skin blend reproduces the rest skin at the reference pose', np.abs(skin - rest).max() < 1e-9,
         vertices=len(rest), max_m=f'{np.abs(skin - rest).max():.1e}')

    # ---- 7. joint opening, OpenSim vs anatomical pivots
    from scipy.spatial import cKDTree
    piv = Pa.pivots()
    bones = {}
    for s in P.segments:
        pts = []
        for bid in P._binding['segment_named_bones'][s]:
            with gzip.open(ROOT / ents[bid]['reference_geometry']['path']) as f:
                pts.append(np.asarray(json.load(f)['positions'], float).reshape(-1, 3))
        Pp = np.concatenate(pts)
        bones[s] = Pp[:: max(1, len(Pp) // 40000)]
    opening = []
    for jname, cname in (('walker_knee_r', 'knee_angle_r'), ('hip_l', 'hip_flexion_l'),
                         ('elbow_r', 'elbow_flex_r'), ('acromial_l', 'arm_flex_l'),
                         ('ankle_r', 'ankle_angle_r'), ('back', 'lumbar_extension'),
                         ('radioulnar_r', 'pro_sup_r')):
        pa, ch = piv[jname]['parent'], piv[jname]['child']
        a, b = bones[pa], bones[ch]
        dist, idx = cKDTree(b).query(a)
        sel = np.argsort(dist)[:64]
        xa, xb, d0 = a[sel], b[idx[sel]], dist[sel]
        lo, hi = kin.coordinates[cname]['range']
        if hi - lo > 2 * np.pi:
            lo, hi = -np.pi / 2, np.pi / 2
        row = dict(joint=jname, coordinate=cname, pivot_offset_m=float(np.linalg.norm(
            piv[jname]['point_m'] - (P.A @ (kin.forward(kin.complete({}, P.reference_pose))[ch]
                                             @ kin.joint_of[ch]['Xcf']))[:3, 3])))
        for mode, poser in (('opensim', P), ('anatomical', Pa)):
            worst = 0.0
            for v in np.linspace(lo, hi, 9):
                M = poser.pose_from_coordinates({cname: float(v)}).segment_motion
                Ma, Mb = M[poser.seg_index[pa]], M[poser.seg_index[ch]]
                ya = xa @ Ma[:3, :3].T + Ma[:3, 3]
                yb = xb @ Mb[:3, :3].T + Mb[:3, 3]
                worst = max(worst, float(np.abs(np.linalg.norm(ya - yb, axis=1) - d0).max()))
            row[f'max_opening_{mode}_m'] = worst
        opening.append(row)
        print(f'      {jname:14s} pivot offset {1000 * row["pivot_offset_m"]:5.1f} mm   worst opening: '
              f'OpenSim pivot {1000 * row["max_opening_opensim_m"]:6.1f} mm, anatomical '
              f'{1000 * row["max_opening_anatomical_m"]:6.1f} mm')
    rec['joint_opening'] = opening

    # ---- 8. real motion + timing
    g = json.loads((ROOT / GAIT).read_text())
    frames = [{k: float(v['value']) for k, v in f['joints'].items()} for f in g['frames']]
    t1 = time.time()
    ok = True
    for q in frames:
        p = P.pose_from_coordinates(q)
        R = p.segment_motion[:, :3, :3]
        ok &= bool(np.isfinite(p.centroid_m).all()) and bool(
            np.allclose(np.einsum('sij,sik->sjk', R, R), np.eye(3), atol=1e-9))
    per = (time.time() - t1) / len(frames)
    # what anatomical pivots cost: how far they carry the anatomy off the simulated segments
    dev = [np.linalg.norm(Pa.pose_from_coordinates(q).centroid_m - P.pose_from_coordinates(q).centroid_m,
                          axis=1) for q in frames[::5]]
    dev = np.concatenate(dev)
    rec['anatomical_pivot_departure_from_simulated_segments_m'] = dict(
        max=float(dev.max()), median=float(np.median(dev)), p95=float(np.percentile(dev, 95)),
        frames=len(frames[::5]), note='entity centroids, anatomical minus OpenSim pivots, gait-best')
    print(f'      anatomical pivots carry entities off the simulated segments by median '
          f'{1000 * np.median(dev):.1f} mm, 95th pct {1000 * np.percentile(dev, 95):.1f}, max {1000 * dev.max():.1f} mm (gait-best)')
    nat = 0
    t2 = time.time()
    for f in NATIVE:
        if (ROOT / f).exists():
            for fr in native_frames(json.loads((ROOT / f).read_text())):
                P.pose_from_native(fr)
                nat += 1
    t3 = time.time()
    P.skin_vertices(p)
    skin_s = time.time() - t3
    gate('gait-best trajectory and native frames posed: finite, rigid', ok and nat > 0,
         gait_frames=len(frames), native_frames=nat, ms_per_pose=f'{1000 * per:.2f}',
         skin_ms=f'{1000 * skin_s:.0f}')
    rec['timing'] = dict(pose_from_coordinates_s=per, skin_vertices_s=skin_s, gait_frames=len(frames),
                         native_frames=nat)
    OUT.write_text(json.dumps(rec, indent=1, default=lambda o: getattr(o, 'tolist', str)()) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    return 0 if all(g['passed'] for g in rec['gates']) else 1


if __name__ == '__main__':
    sys.exit(main())
