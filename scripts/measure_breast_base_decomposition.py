#!/usr/bin/env python3
"""The shape-mismatch measurement RESTRICTED TO BASE NODES -- checking my own finding.

Last hour I concluded "no placement seats these breasts" from the spread of signed distance over
ALL breast vertices. But only 31% of a prepared breast's nodes are base nodes (4,097 of 13,216 for
s1159-left); the other 69% are the front surface, which sits far from the chest wall and never has
to seat. Including them can only inflate the spread, and could have manufactured the conclusion.

`prepared.npz` carries the repo's own `base` index and `gap0`, the signed gap at those nodes. This
redoes the measurement on the base alone, with the same objective (minimise the VARIANCE of signed
distance under a bounded translation) and the same known answer.
"""
import numpy as np, json, gzip, sys
from pathlib import Path
from scipy.spatial import cKDTree
ROOT = Path('/home/brandonin/Documents/IHM-1')
LIMIT_MM = 25.0

ents = json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())['entities']
WV = []
for e in ents:
    n = e.get('name','').lower()
    if e.get('reference_geometry') and ('pectoralis major' in n or 'rib' in n or 'sternum' in n):
        g = json.loads(gzip.decompress((ROOT/e['reference_geometry']['path']).read_bytes()))
        WV.append(np.asarray(g['positions'],float).reshape(-1,3))
WV = np.vstack(WV); tree = cKDTree(WV); WVz = WV[:,2]
print(f"chest wall: {len(WV):,} vertices")

def depth(V, t=np.zeros(3)):
    P = V + t
    _, i = tree.query(P, k=1, workers=-1)
    return P[:,2] - WVz[i]

def fit(V, limit=LIMIT_MM):
    t = np.zeros(3); cost = lambda tt: float(depth(V,tt).var())
    c, step = cost(t), 0.008
    while step > 0.0002:
        imp = False
        for k in range(3):
            for s in (+1.,-1.):
                tt = t.copy(); tt[k] += s*step
                if np.linalg.norm(tt)*1e3 > limit: continue
                cc = cost(tt)
                if cc < c - 1e-16: t, c, imp = tt, cc, True
        if not imp: step *= 0.5
    return t, float(np.sqrt(c))*1e3

# KNOWN ANSWER on base nodes
P0 = np.load(ROOT/'data/derived/female-breast-seated-v1/s1159/left/prepared.npz')
B0 = P0['X'][P0['base']]
t0, sd0 = fit(B0); t1, sd1 = fit(B0 + np.array([0,0,-0.010]))
miss = float(np.linalg.norm((t1 - np.array([0,0,0.010])) - t0)*1e3)
ok = miss < 2.0 and abs(sd1-sd0) < 0.3
print(f"\nKNOWN ANSWER (base nodes): 10 mm posterior displacement must be compensated exactly.")
print(f"  miss {miss:.2f} mm, residual sd moves {sd1-sd0:+.3f} mm   {'PASS' if ok else 'FAIL'}")
if not ok: sys.exit("known answer FAILED")

print(f"\n{'subject':8s} {'side':6s} {'n base':>7s} {'sd ALL':>7s} {'sd BASE':>8s} {'|t| mm':>7s} {'residual':>9s} {'removed':>8s}")
ALLV = {('s0790','left'):4.86,('s0790','right'):4.67,('s1067','left'):5.18,('s1067','right'):6.40,
        ('s1159','left'):4.10,('s1159','right'):4.44,('s0970','left'):6.49,('s0970','right'):6.14}
rows=[]
for sid in ('s0790','s1067','s1159','s0970'):
    for side in ('left','right'):
        f = ROOT/f'data/derived/female-breast-seated-v1/{sid}/{side}/prepared.npz'
        if not f.exists(): continue
        P = np.load(f); B = P['X'][P['base']]
        sdb = float(depth(B).std())*1e3
        t, sda = fit(B); tm = float(np.linalg.norm(t)*1e3)
        rem = 1.0 - sda/max(sdb,1e-12)
        print(f"{sid:8s} {side:6s} {len(B):7d} {ALLV[(sid,side)]:7.2f} {sdb:8.2f} {tm:7.2f} {sda:9.2f} {rem:7.1%}")
        rows.append(dict(subject=sid,side=side,n_base=int(len(B)),sd_all_mm=ALLV[(sid,side)],
                         sd_base_mm=sdb,translation_mm=tm,residual_sd_mm=sda,fraction_removed=float(rem)))
rem = np.array([r['fraction_removed'] for r in rows]); res = np.array([r['residual_sd_mm'] for r in rows])
sdb = np.array([r['sd_base_mm'] for r in rows]); sda_all = np.array([r['sd_all_mm'] for r in rows])
print(f"\n  base-node spread is {np.median(sdb/sda_all):.2f}x the all-vertex spread (median)")
print(f"  best translation removes a median of {np.median(rem):.1%} of the BASE spread")
print(f"  residual after translation: {res.min():.1f}-{res.max():.1f} mm sd")
print("\n  VERDICT: " + ("mostly RIGID on the base -- my all-vertex conclusion was wrong."
      if np.median(rem) > 0.5 else
      "still mostly a SHAPE MISMATCH on the base alone. The all-vertex conclusion holds."))
Path(ROOT/'out/breast_base_only_decomposition.json').write_text(json.dumps(rows, indent=2)+"\n")
print(f"\nwrote out/breast_base_only_decomposition.json")
