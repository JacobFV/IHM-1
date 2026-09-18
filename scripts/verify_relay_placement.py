#!/usr/bin/env python3
"""Known answers for the relays on the cord: route lengths against Z-Anatomy centrelines.

    CUDA_VISIBLE_DEVICES="" .venv/bin/python scripts/verify_relay_placement.py

Run after `build_spinal_cord_levels.py`, `build_body_peripheral.py` and
`build_dermatome_patches.py`.  Writes `data/derived/canonical/relay_placement_report.json`.

WHAT IT CHECKS, in the order a reader should believe it.

1. INSTRUMENT.  IHM's length method (`1.15 x` straight line to the relay for bindings and
   patches, upper median of those for a binding-backed nerve, sum of the polyline otherwise)
   is re-implemented here and must reproduce every declared nerve length in the NEW
   `peripheral.json` exactly, and -- with the pre-18-Sep relays and authored endpoints put
   back -- every length in the OLD one (146 of each; pass the old file as the first
   argument to check that direction).  Only then is the before/after table
   a statement about the relays rather than about this script.
2. RELAYS.  Old and new position; each new spinal relay's distance to the registered dura
   centreline at its level; the vertebral-canal enclosure at that level (from
   `spinal_cord_levels.json`).
3. KNOWN ANSWERS.  Each route with a Z-Anatomy counterpart against the length along that
   atlas nerve, from the point nearest the route's own endpoint to the cord at the relay's
   level (target: within ~10%).  The Z-Anatomy centrelines are the ones IBM-1's site exporter
   routed (`IBM-1/site/data/body.js`, read only), carried back to blender world through the
   exporter's own chain and then into this frame by the SAME landmark registration the cord
   used; the vagus and plantar curves come from IBM-1's audit extraction.  The cord chain of
   the site payload, carried the same way, must land on this repo's dura centreline.
4. The sciatic-tibial path, conus to sole, must exceed 1 m ("can exceed one metre").
"""
from __future__ import annotations

import gzip
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
from ihm.assembly.anatomy import LandmarkRegistration, ROTATION  # noqa: E402
from peripheral_route_catalog import SOMATIC, VISCERAL, ANCHORED   # noqa: E402

IBM = ROOT.parent / 'IBM-1'
CAN = ROOT / 'data/derived/canonical'
SITE = IBM / 'site/data/body.js'
AUDIT_CACHE = IBM / '.cache/zanatomy_nerve_audit.json'
OUT = CAN / 'relay_placement_report.json'

#: the relays as typed in build_body_peripheral.py until 18 Sep 2026, and the solitary default
OLD_RELAY_Y = {'cranial': .65, 'cervical': .53, 'thoracic': .30, 'lumbar': .04, 'sacral': -.08}
OLD_SOLITARY = [.012, .655, -.035]
TOLERANCE = 0.10


def arclen(p):
    return float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum())


def chains(v, edges):
    adj: dict = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    used, out = set(), []
    for st in [k for k, n in adj.items() if len(n) != 2] + list(adj):
        for nb in adj[st]:
            k = (min(st, nb), max(st, nb))
            if k in used:
                continue
            used.add(k)
            ch, prev, cur = [st, nb], st, nb
            while len(adj[cur]) == 2:
                nx = [n for n in adj[cur] if n != prev][0]
                k = (min(cur, nx), max(cur, nx))
                if k in used:
                    break
                used.add(k)
                ch.append(nx)
                prev, cur = cur, nx
            out.append(v[ch])
    return out


def old_relays(per):
    out = {}
    for r in per['relays']:
        side, level = r['side'], r['id'].split(f"-{r['side']}-")[1]
        sign = 1.0 if side == 'left' else -1.0
        if level in OLD_RELAY_Y:
            out[r['id']] = np.array([sign * .012, OLD_RELAY_Y[level], -.04])
        elif level == 'solitary':
            out[r['id']] = np.array([sign * OLD_SOLITARY[0], *OLD_SOLITARY[1:]])
        else:
            out[r['id']] = np.asarray(r['position_m'], float)
    return out


def ihm_lengths(per, E, relays, first_points):
    """IHM's own method (build_body_peripheral / enrich_peripheral_routes), relays passed in."""
    per_nerve: dict = {}
    binding_len = {}
    for b in per['muscle_bindings']:
        d = 1.15 * float(np.linalg.norm(np.asarray(E[b['canonical_entity_id']]['centroid_m'])
                                        - relays[b['relay_id']]))
        binding_len[b['muscle_id']] = max(.03, d)
        per_nerve.setdefault(b['nerve_id'], []).append((b['path_length_m'], binding_len[b['muscle_id']],
                                                        b['muscle_id']))
    for p in per['receptor_patches']:
        d = 1.15 * float(np.linalg.norm(np.asarray(p['position_m']) - relays[p['relay_id']]))
        per_nerve.setdefault(p['nerve_id'], []).append((p['path_length_m'], max(.03, d), p['id']))
    out = {}
    for n in per['nerves']:
        if n['id'] in per_nerve:
            # the builder sorts by its OWN (new) lengths; re-sort by the lengths computed here
            c = sorted(per_nerve[n['id']], key=lambda t: (t[1], t[2]))
            out[n['id']] = c[len(c) // 2][1]
        else:
            r = relays[n['relay_id']]
            p0 = np.asarray(first_points[n['id']], float)
            sgn = 1.0 if n['side'] == 'left' else -1.0
            q = [p0, r] if n.get('route_kind') == 'special_sense' else [p0, np.array([sgn * .04, r[1], r[2]]), r]
            out[n['id']] = float(sum(np.linalg.norm(a - b) for a, b in zip(q, q[1:])))
    return out


def main() -> int:
    per = json.loads((CAN / 'peripheral.json').read_text())
    der = json.loads((CAN / 'dermatomes.json').read_text())
    lev = json.loads((CAN / 'spinal_cord_levels.json').read_text())
    anat = json.loads((CAN / 'anatomy.json').read_text())
    for label, d in (('peripheral.json', per), ('dermatomes.json', der), ('spinal_cord_levels.json', lev)):
        if d['frame']['id'] != anat['frame']['id']:
            raise SystemExit(f'{label} is in {d["frame"]["id"]}, anatomy.json in {anat["frame"]["id"]}')
    E = {e['id']: e for e in anat['entities']}
    N = {n['id']: n for n in per['nerves']}
    R_new = {r['id']: np.asarray(r['position_m'], float) for r in per['relays']}
    R_old = old_relays(per)
    rec: dict = {'frame': anat['frame']['id']}

    # ---- 1. instrument, both directions
    first_new = {n['id']: n['points_m'][0] for n in per['nerves'] if n.get('points_m')}
    first_old = dict(first_new)
    catalog = {name: pos for name, _, _, pos in SOMATIC + VISCERAL}
    for n in per['nerves']:
        nm = n['id'].split(f"-{n['side']}-", 1)[1]
        if nm in ANCHORED:
            sgn = 1.0 if n['side'] == 'left' else -1.0
            first_old[n['id']] = [sgn * catalog[nm][0], *catalog[nm][1:]]
    new = ihm_lengths(per, E, R_new, first_new)
    old = ihm_lengths(per, E, R_old, first_old)
    worst_new = max(abs(new[k] - N[k]['path_length_m']) for k in new)
    rec['instrument'] = dict(nerves=len(new), worst_abs_m_vs_declared=worst_new,
                             reproduces_declared=worst_new < 1e-9)
    print(f'instrument: IHM method reproduces {sum(abs(new[k] - N[k]["path_length_m"]) < 1e-9 for k in new)}'
          f'/{len(new)} declared lengths (worst {worst_new:.1e} m)')
    if worst_new >= 1e-9:
        raise SystemExit('the re-implementation does not reproduce peripheral.json; nothing below is trustworthy')
    if len(sys.argv) > 1:            # the pre-change peripheral.json, to check the OLD direction
        prev = {n['id']: n['path_length_m'] for n in json.loads(Path(sys.argv[1]).read_text())['nerves']}
        worst_old = max(abs(old[k] - prev[k]) for k in old)
        rec['instrument']['worst_abs_m_vs_previous_file'] = worst_old
        rec['instrument']['previous_file'] = sys.argv[1]
        print(f'instrument: old relays + authored endpoints reproduce {sum(abs(old[k] - prev[k]) < 1e-9 for k in old)}'
              f'/{len(old)} lengths of {sys.argv[1]} (worst {worst_old:.1e} m)')
        if worst_old >= 1e-9:
            raise SystemExit('the old-relay reconstruction does not reproduce the previous file')
    rec['old_lengths_m'] = old
    rec['new_lengths_m'] = new

    # ---- 2. relays
    print('\nrelays                               old y     new position (m)            off-cord  canal')
    rows = []
    for rid in sorted(R_new):
        lvl = rid.split('-', 3)[3]
        if lvl not in OLD_RELAY_Y and lvl != 'solitary':
            continue
        placed = lev['relays'][rid]
        c = placed.get('dura_centreline_m')
        row = dict(relay=rid, old_m=R_old[rid].tolist(), new_m=R_new[rid].tolist(),
                   moved_m=float(np.linalg.norm(R_new[rid] - R_old[rid])),
                   level=placed.get('segment') or placed.get('structure'),
                   vertebra=placed.get('vertebra'),
                   distance_to_dura_centreline_m=(float(np.linalg.norm(R_new[rid] - np.asarray(c)))
                                                  if c else None),
                   canal=placed.get('canal_check'), source=placed['position_source'])
        rows.append(row)
        if rid.startswith('peripheral-relay-left'):
            cc = row['canal']
            off = '' if c is None else f"{1000 * row['distance_to_dura_centreline_m']:.1f} mm"
            canal = '' if not cc else (f"inside, {cc['bins_with_bone']}/24, "
                                       f"{1000 * cc['clearance_m']:.1f} mm clear")
            where = row['level'] + (', ' + row['vertebra'] if row['vertebra'] else '')
            print(f'  {rid:34s} {R_old[rid][1]:+.3f}   {str(np.round(R_new[rid], 4).tolist()):28s} '
                  f'{off:9s} {canal}  [{where}]')
    rec['relays'] = rows

    # ---- 3. Z-Anatomy centrelines in THIS frame
    stored = anat['registrations']['z_anatomy']
    reg = LandmarkRegistration(np.array([m['source_rotated_m'] for m in stored['landmarks']]),
                               np.array([m['target_m'] for m in stored['landmarks']]),
                               smoothing=stored['smoothing_prior'])
    zm = json.loads((ROOT / 'data/derived/anatomy/extended/manifest_fragment.json').read_text())['models'][0]
    bpm = json.loads((CAN / 'manifest_fragment.json').read_text())['models'][0]
    ZS = bpm['bounds']['max'][1] / zm['bounds']['max'][1]
    dt = zm['display_transform']
    Rz, sz, tz = np.asarray(dt['rotation']), float(dt['scale']), np.asarray(dt['translation'])
    if not np.allclose(Rz @ Rz.T, np.eye(3), atol=1e-9):
        raise SystemExit('display rotation is not orthonormal; the chain cannot be inverted')

    def site_to_here(p):                  # site metres -> blend world -> registered canonical
        blend = ((np.asarray(p, float) / ZS - tz) / sz) @ Rz
        return reg.transform(blend @ ROTATION.T)

    def blend_to_here(p):
        return reg.transform(np.asarray(p, float) @ ROTATION.T)

    s = SITE.read_text()
    site = json.loads(s[s.index('{'): s.rindex('}') + 1])['nerves']
    P = np.asarray(site['p'], float).reshape(-1, 3)
    owner = np.zeros(len(P), int)
    for i, (o, n, *_) in enumerate(site['chains']):
        owner[o:o + n] = i
    Ph = site_to_here(P)
    cord_rows = np.asarray(lev['centreline']['rows'])
    c0 = site['chains'][0]
    assert site['names'][0] == 'Spinal cord'
    sc = Ph[c0[0]: c0[0] + c0[1]]
    inside = sc[(sc[:, 1] > cord_rows[0, 1]) & (sc[:, 1] < cord_rows[-1, 1])]
    cx = np.interp(inside[:, 1], cord_rows[:, 1], cord_rows[:, 0])
    cz = np.interp(inside[:, 1], cord_rows[:, 1], cord_rows[:, 2])
    cord_gap = np.hypot(inside[:, 0] - cx, inside[:, 2] - cz)
    rec['site_cord_vs_levels_centreline_m'] = dict(median=float(np.median(cord_gap)), max=float(cord_gap.max()))
    print(f'\nsite cord chain carried into this frame vs the dura centreline: median '
          f'{1000 * np.median(cord_gap):.1f} mm, max {1000 * cord_gap.max():.1f} mm')

    paths = {p['name']: np.asarray(p['idx']) for p in site['paths']}

    def z_to_level(name, endpoint, level_y):
        """along the site path from the point nearest `endpoint` to the cord at `level_y`."""
        idx = paths[name]
        pts = Ph[idx]
        seg = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
        kc = int(np.argmax([site['names'][owner[i]] == 'Spinal cord' for i in idx]))
        k = int(np.argmin(np.linalg.norm(pts[:kc] - endpoint, axis=1)))
        kl = kc + int(np.argmin(np.abs(pts[kc:, 1] - level_y)))
        return dict(length_m=float(seg[kl] - seg[k]), gap_m=float(np.linalg.norm(pts[k] - endpoint)),
                    cord_entry_y_m=float(pts[kc, 1]), level_y_m=float(level_y), pts=pts, seg=seg, kc=kc,
                    stops_short=bool(k == 0 and np.linalg.norm(pts[k] - endpoint) > .01))

    ka = []

    def compare(route, zname, ihm_len, endpoint, level_y, note=''):
        z = z_to_level(zname, endpoint, level_y)
        row = dict(route=route, z_anatomy=zname, ihm_old_m=old[route], ihm_new_m=ihm_len,
                   z_m=z['length_m'], endpoint_gap_m=z['gap_m'], cord_entry_y_m=z['cord_entry_y_m'],
                   relay_level_y_m=z['level_y_m'], note=note, z_stops_short=z['stops_short'])
        row['ratio_old'] = row['ihm_old_m'] / row['z_m']
        row['ratio_new'] = row['ihm_new_m'] / row['z_m']
        row['within_tolerance'] = abs(row['ratio_new'] - 1) <= TOLERANCE
        ka.append(row)
        return z

    for side, sfx in (('left', 'l'), ('right', 'r')):
        for ihm, zn in (('tibial', 'Tibial nerve'), ('deep_fibular', 'Deep fibular nerve'),
                        ('superficial_fibular', 'Superficial fibular nerve'),
                        ('femoral', 'Femoral nerve'), ('median', 'Median nerve'),
                        ('radial', 'Radial nerve')):
            n = N[f'peripheral-nerve-{side}-{ihm}']
            ep = np.asarray(E[n['endpoint_id']]['centroid_m'])
            compare(n['id'], f'{zn}.{sfx}', n['path_length_m'], ep, R_new[n['relay_id']][1],
                    note=f'endpoint {E[n["endpoint_id"]]["name"]}')

    # medial plantar: the site's tibial path (distal end -> cord at the sacral level) plus the
    # medial plantar curve from the audit extraction, from its proximal end to the point
    # nearest the IHM endpoint.  This is the sciatic-tibial path, conus to sole.
    aud = json.loads(AUDIT_CACHE.read_text())
    curves = {c['name']: c for c in aud['curves']}

    def longest(nm):
        c = curves[nm]
        cs = [blend_to_here(np.asarray(c['v'])[ch_idx]) for ch_idx in
              [np.asarray(x) for x in chains(np.arange(len(c['v'])), c['e'])]]
        m = max(cs, key=arclen)
        return m if m[0, 1] >= m[-1, 1] else m[::-1]     # top first

    for side, sfx in (('left', 'l'), ('right', 'r')):
        n = N[f'peripheral-nerve-{side}-medial_plantar']
        ep = np.asarray(E[n['endpoint_id']]['centroid_m'])
        tib = z_to_level(f'Tibial nerve.{sfx}', np.asarray(Ph[paths[f'Tibial nerve.{sfx}']][0]),
                         R_new[n['relay_id']][1])
        mp = longest(f'Medial plantar nerve.{sfx}')
        ms = np.r_[0, np.cumsum(np.linalg.norm(np.diff(mp, axis=0), axis=1))]
        k = int(np.argmin(np.linalg.norm(mp - ep, axis=1)))
        full = tib['length_m'] + float(ms[-1])
        to_ep = tib['length_m'] + float(ms[k])
        row = dict(route=n['id'], z_anatomy=f'Tibial nerve.{sfx} (site path) + Medial plantar nerve.{sfx}',
                   ihm_old_m=old[n['id']], ihm_new_m=n['path_length_m'], z_m=to_ep,
                   z_conus_to_sole_m=full, endpoint_gap_m=float(np.linalg.norm(mp[k] - ep)),
                   relay_level_y_m=float(R_new[n['relay_id']][1]),
                   note=f'endpoint {E[n["endpoint_id"]]["name"]}',
                   z_stops_short=bool(k == len(mp) - 1 and np.linalg.norm(mp[k] - ep) > .01))
        row['ratio_old'] = row['ihm_old_m'] / row['z_m']
        row['ratio_new'] = row['ihm_new_m'] / row['z_m']
        row['within_tolerance'] = abs(row['ratio_new'] - 1) <= TOLERANCE
        ka.append(row)

    # vagus: Z curve from its top (brainstem) to the point nearest the gastric endpoint.  The
    # declared right route is the MIRROR of the left (IBM-1's bilateral visceral contract), so its
    # endpoint is not on the stomach; the right vagus is also computed with its own per-side
    # anchor (posterior gastric wall) by the same polyline method, which is what the Z right
    # vagus -- the only one that reaches the abdomen -- can be compared with.
    from enrich_peripheral_routes import anchor_point
    for side, sfx in (('left', 'l'), ('right', 'r')):
        n = N[f'peripheral-nerve-{side}-vagus']
        if side == 'left':
            ep, length, label = np.asarray(n['points_m'][0]), n['path_length_m'], n['id']
            how = f'declared; endpoint on {n["endpoint_source"]["entity_name"]} by {n["endpoint_source"]["rule"]}'
        else:
            ep = np.asarray(anchor_point(ROOT, E, 'vagus', 'right', catalog['vagus'])[0])
            r = R_new[n['relay_id']]
            q = [ep, np.array([-.04, r[1], r[2]]), r]
            length = float(sum(np.linalg.norm(a - b) for a, b in zip(q, q[1:])))
            label = n['id'] + ' (per-side anchor, posterior wall; NOT the declared route)'
            how = (f'declared right route is the mirror of the left, {1000 * n["path_length_m"]:.0f} mm, '
                   'endpoint off the stomach; this row re-anchors it on the posterior wall')
        vg = longest(f'Vagus nerve (X).{sfx}')
        vs = np.r_[0, np.cumsum(np.linalg.norm(np.diff(vg, axis=0), axis=1))]
        k = int(np.argmin(np.linalg.norm(vg - ep, axis=1)))
        row = dict(route=label, z_anatomy=f'Vagus nerve (X).{sfx}', ihm_old_m=old[n['id']],
                   ihm_new_m=length, z_m=float(vs[k]),
                   endpoint_gap_m=float(np.linalg.norm(vg[k] - ep)), z_top_y_m=float(vg[0, 1]),
                   z_stops_short=bool((k == len(vg) - 1 and np.linalg.norm(vg[k] - ep) > .01)
                                      or vg[:, 1].min() > ep[1] + .01),
                   note='Z top of curve to the point nearest the gastric endpoint; ' + how)
        row['ratio_old'] = row['ihm_old_m'] / row['z_m']
        row['ratio_new'] = row['ihm_new_m'] / row['z_m']
        row['within_tolerance'] = abs(row['ratio_new'] - 1) <= TOLERANCE
        ka.append(row)
    for r in ka:
        # a Z curve whose nearest point to the IHM endpoint is its own terminal end STOPS SHORT of
        # that endpoint: its length is a lower bound, not a same-scope known answer.  (Rule written
        # after the first table, where it marked radial and left vagus; both rows stay printed.)
        r['comparable'] = not r['z_stops_short']
    rec['known_answers'] = ka
    comp = [r for r in ka if r['comparable']]
    rec['known_answer_summary'] = dict(
        comparable=len(comp), within_tolerance=sum(r['within_tolerance'] for r in comp),
        outside=[r['route'] for r in comp if not r['within_tolerance']],
        not_comparable=[dict(route=r['route'], reason='Z curve ends %.0f mm short of the IHM endpoint'
                             % (1000 * r['endpoint_gap_m'])) for r in ka if not r['comparable']],
        old_within_tolerance=sum(abs(r['ratio_old'] - 1) <= TOLERANCE for r in comp))

    print(f'\nknown answers (target within {100 * TOLERANCE:.0f}%)')
    print(f'  {"route":40s} {"old":>6s} {"new":>6s} {"Z":>6s} {"old/Z":>6s} {"new/Z":>6s} {"gap":>6s}')
    for r in ka:
        print(f'  {r["route"][17:57]:40s} {1000 * r["ihm_old_m"]:6.0f} {1000 * r["ihm_new_m"]:6.0f} '
              f'{1000 * r["z_m"]:6.0f} {r["ratio_old"]:6.2f} {r["ratio_new"]:6.2f} '
              f'{1000 * r["endpoint_gap_m"]:4.0f}mm  '
              f'{"Z STOPS SHORT (lower bound)" if not r["comparable"] else "ok" if r["within_tolerance"] else "OUTSIDE"}')
    sm = rec['known_answer_summary']
    print(f'  comparable {sm["comparable"]}: {sm["within_tolerance"]} within {100 * TOLERANCE:.0f}% now, '
          f'{sm["old_within_tolerance"]} before; outside now: {sm["outside"] or "-"}')

    # ---- 4. over a metre, conus to sole
    sole = [r for r in ka if r['route'].endswith('medial_plantar')]
    der_foot = max((p for p in der['patches'] if p['region'] == 'foot'), key=lambda p: p['path_length_m'])
    foot_patch = next(p for p in per['receptor_patches'] if p['id'] == 'peripheral-skin-left-foot')
    longest_route = max(N.values(), key=lambda n: n['path_length_m'])
    rec['over_a_metre'] = dict(
        z_conus_to_sole_m=[r['z_conus_to_sole_m'] for r in sole],
        ihm_medial_plantar_m=[r['ihm_new_m'] for r in sole],
        ihm_foot_receptor_patch_m=foot_patch['path_length_m'],
        longest_dermatome_patch=dict(id=der_foot['id'], m=der_foot['path_length_m']),
        longest_nerve=dict(id=longest_route['id'], m=longest_route['path_length_m']))
    rec['over_a_metre']['passed'] = bool(all(r['ihm_new_m'] > 1.0 for r in sole)
                                         and foot_patch['path_length_m'] > 1.0)
    print(f'\nconus -> sole: Z {", ".join(f"{1000 * v:.0f}" for v in rec["over_a_metre"]["z_conus_to_sole_m"])} mm; '
          f'IHM medial plantar {", ".join(f"{1000 * r["ihm_new_m"]:.0f}" for r in sole)} mm '
          f'(was {", ".join(f"{1000 * r["ihm_old_m"]:.0f}" for r in sole)}); '
          f'sole receptor patch {1000 * foot_patch["path_length_m"]:.0f} mm; longest dermatome patch '
          f'{1000 * der_foot["path_length_m"]:.0f} mm -> {"PASS" if rec["over_a_metre"]["passed"] else "FAIL"} (> 1 m)')

    # ---- every route, old -> new
    ratios = np.array([new[k] / old[k] for k in new])
    rec['all_routes'] = dict(ratio_min=float(ratios.min()), ratio_median=float(np.median(ratios)),
                             ratio_max=float(ratios.max()))
    print(f'\nall 146 routes: new/old x{ratios.min():.2f} .. x{ratios.max():.2f}, median x{np.median(ratios):.2f}')
    OUT.write_text(json.dumps(rec, indent=1, default=lambda o: getattr(o, 'tolist', str)()) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
