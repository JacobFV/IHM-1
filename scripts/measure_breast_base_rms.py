#!/usr/bin/env python3
"""Closing a gap in my own measure: variance ignores the MEAN.

The base-node decomposition minimised the VARIANCE of signed distance and reported the residual
spread, 2.3-3.7 mm sd. But variance is invariant to a uniform offset: a base sitting uniformly
40 mm inside the chest wall has near-zero spread and would score as perfectly seatable while
still being 40 mm inside. The measure answered "can a translation make the base PARALLEL to the
wall", not "can a translation get the base ONTO it".

That matters because it is the gap between the two numbers that has never added up: the residual
spread is ~3 mm, and the drive that reached fraction 1.0000 deformed tissue by 67.96 mm. Three
millimetres of shape disagreement does not cost sixty-eight millimetres of strain.

The objective that asks the real question is the mean SQUARED SIGNED DISTANCE, since

    RMS^2  =  mean^2  +  variance

so it penalises being buried and being mis-shaped together, and is not degenerate -- moving the
breast away from the chest increases it, unlike minimising penetration. Reporting the residual
RMS with its mean and sd parts separated says which of the two is actually fatal.

KNOWN ANSWER, same as before: a 10 mm posterior displacement must be compensated exactly, and the
residual RMS it achieves must be unchanged.
"""
import numpy as np, json, gzip, sys
from pathlib import Path
from scipy.spatial import cKDTree
ROOT = Path('/home/brandonin/Documents/IHM-1'); LIMIT_MM = 25.0

ents = json.loads((ROOT/'data/derived/canonical/anatomy.json').read_text())['entities']
WV=[]
for e in ents:
    n=e.get('name','').lower()
    if e.get('reference_geometry') and ('pectoralis major' in n or 'rib' in n or 'sternum' in n):
        g=json.loads(gzip.decompress((ROOT/e['reference_geometry']['path']).read_bytes()))
        WV.append(np.asarray(g['positions'],float).reshape(-1,3))
WV=np.vstack(WV); tree=cKDTree(WV); WVz=WV[:,2]

def depth(V,t=np.zeros(3)):
    P=V+t; _,i=tree.query(P,k=1,workers=-1); return P[:,2]-WVz[i]

def fit_rms(V, limit=LIMIT_MM):
    t=np.zeros(3); cost=lambda tt: float((depth(V,tt)**2).mean())
    c,step=cost(t),0.008
    while step>0.0002:
        imp=False
        for k in range(3):
            for s in (+1.,-1.):
                tt=t.copy(); tt[k]+=s*step
                if np.linalg.norm(tt)*1e3>limit: continue
                cc=cost(tt)
                if cc<c-1e-16: t,c,imp=tt,cc,True
        if not imp: step*=0.5
    d=depth(V,t)*1e3
    return t, float(np.sqrt(c))*1e3, float(d.mean()), float(d.std())

P0=np.load(ROOT/'data/derived/female-breast-seated-v1/s1159/left/prepared.npz')
B0=P0['X'][P0['base']]
t0,r0,_,_=fit_rms(B0); t1,r1,_,_=fit_rms(B0+np.array([0,0,-0.010]))
miss=float(np.linalg.norm((t1-np.array([0,0,0.010]))-t0)*1e3)
ok = miss<2.0 and abs(r1-r0)<0.3
print(f"KNOWN ANSWER: 10 mm posterior displacement compensated exactly, residual RMS unchanged.")
print(f"  miss {miss:.2f} mm, residual RMS moves {r1-r0:+.3f} mm   {'PASS' if ok else 'FAIL'}")
if not ok: sys.exit("known answer FAILED")

print(f"\n{'subject':8s} {'side':6s} {'RMS before':>10s} {'|t|':>6s} {'RMS after':>10s} "
      f"{'  = mean':>9s} {'+ sd':>7s} {'deepest in':>11s}")
rows=[]
for sid in ('s0790','s1067','s1159','s0970'):
    for side in ('left','right'):
        f=ROOT/f'data/derived/female-breast-seated-v1/{sid}/{side}/prepared.npz'
        if not f.exists(): continue
        P=np.load(f); B=P['X'][P['base']]
        d0=depth(B)*1e3; rms0=float(np.sqrt((d0**2).mean()))
        t,rms,mean,sd=fit_rms(B); tm=float(np.linalg.norm(t)*1e3)
        deepest=float(-(depth(B,t)*1e3).min())
        bounded = tm > LIMIT_MM-0.5
        print(f"{sid:8s} {side:6s} {rms0:10.2f} {tm:6.2f}{'*' if bounded else ' '} {rms:10.2f} "
              f"{mean:9.2f} {sd:7.2f} {deepest:10.2f}")
        rows.append(dict(subject=sid,side=side,rms_before_mm=rms0,translation_mm=tm,
                         rms_after_mm=rms,residual_mean_mm=mean,residual_sd_mm=sd,
                         deepest_behind_after_mm=deepest,at_bound=bool(bounded)))
R=np.array([r['rms_after_mm'] for r in rows]); M=np.array([r['residual_mean_mm'] for r in rows])
S=np.array([r['residual_sd_mm'] for r in rows]); D=np.array([r['deepest_behind_after_mm'] for r in rows])
print(f"\n  residual RMS after the best translation: {R.min():.2f}-{R.max():.2f} mm")
print(f"    of which mean (still buried): {np.abs(M).min():.2f}-{np.abs(M).max():.2f} mm")
print(f"    of which sd  (mis-shaped)   : {S.min():.2f}-{S.max():.2f} mm")
print(f"  deepest node still behind the wall after the best translation: {D.min():.1f}-{D.max():.1f} mm")
dominant = 'MEAN (still buried)' if np.median(np.abs(M)) > np.median(S) else 'SPREAD (mis-shaped)'
print(f"\n  the residual is dominated by the {dominant}")
Path(ROOT/'out/breast_base_rms_decomposition.json').write_text(json.dumps(rows,indent=2)+"\n")
print("wrote out/breast_base_rms_decomposition.json")
