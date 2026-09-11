"""Symmetric normal shooting: a surface correspondence whose direction comes from the surfaces.

WHY THIS EXISTS. Matching a source point to the NEAREST point on a target surface is biased on
curvature by order d^2/R for a displacement d and a radius R. On this body's ribs, a few millimetres
thick, that bias was measured against a known preimage at mean 0.671 mm, p90 2.281 mm and max
7.152 mm for displacements averaging 1.53 mm -- large enough that a thin-plate spline fitted to
those targets reproduced them to 0.01 mm and still missed the true surface by 1.27 mm RMS
(docs/BODY_PARAMETERS.md, "The deformable chest-wall fit"). Two warp lines in this repo are fitted
on nearest-point targets, so this is written as a reusable instrument rather than inline.

THE RULE. For each source sample, shoot along the SOURCE surface's own outward normal (both ways)
and take the first intersection with the target within a cap. Shoot back from that hit along the
TARGET's normal and keep the pair only where the return lands within a tolerance of the source
point. Direction comes from the surfaces; proximity never enters.

Normals are interpolated from area-weighted vertex normals, not taken from the facet: a
marching-cubes or atlas surface has facet normals that jump between neighbours, and a direction that
jumps gives a correspondence that jumps with it.
"""
import numpy as np
from scipy.spatial import cKDTree

__all__ = ["vertex_normals", "surface_samples", "build_index", "first_hit", "shoot_pairs"]


def vertex_normals(V, F):
    n = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    vn = np.zeros_like(V)
    for k in range(3): np.add.at(vn, F[:, k], n)
    return vn / np.maximum(np.linalg.norm(vn, axis=1, keepdims=True), 1e-30)


def surface_samples(V, F, n, seed=0):
    """Area-weighted samples with their interpolated normals, and the (face, barycentric) they came
    from -- so the same material point can be located on a second mesh with the same connectivity,
    which is what a known-answer test needs."""
    rng = np.random.default_rng(seed)
    tri = V[F]
    a = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    f = rng.choice(len(F), n, p=a / a.sum())
    u = rng.random((n, 1)); v = rng.random((n, 1))
    over = (u + v > 1); u[over] = 1 - u[over]; v[over] = 1 - v[over]
    w = np.hstack([1 - u - v, u, v])
    P = np.einsum('nk,nkj->nj', w, V[F[f]])
    N = np.einsum('nk,nkj->nj', w, vertex_normals(V, F)[F[f]])
    return P, N / np.maximum(np.linalg.norm(N, axis=1, keepdims=True), 1e-30), f, w


def build_index(V, F, tier=99):
    """Two tiers, because one oversized triangle would set the query radius for every ray: face radii
    on these surfaces run ~1.6 mm median against a 28.9 mm max."""
    cent = V[F].mean(1)
    rad = np.linalg.norm(V[F] - cent[:, None, :], axis=2).max(1)
    cut = float(np.percentile(rad, tier))
    small = np.flatnonzero(rad <= cut); large = np.flatnonzero(rad > cut)
    return dict(tree=cKDTree(cent[small]), small=small, large=large, radius=cut)


def _moller_trumbore(P, D, v0, e1, e2):
    h = np.cross(D, e2); a = np.einsum('fj,fj->f', e1, h)
    usable = np.abs(a) > 1e-16
    inv = 1.0 / np.where(usable, a, 1.0)
    s = P - v0
    u = inv * np.einsum('fj,fj->f', s, h)
    q = np.cross(s, e1)
    v = inv * (q @ D)
    t = inv * np.einsum('fj,fj->f', e2, q)
    return t, usable & (u >= -1e-9) & (v >= -1e-9) & (u + v <= 1 + 1e-9)


def first_hit(P, D, V, F, index, cap):
    """Nearest intersection of each ray P + tD with |t| <= cap, shooting BOTH ways along D."""
    v0 = V[F[:, 0]]; e1 = V[F[:, 1]] - v0; e2 = V[F[:, 2]] - v0
    dist = np.full(len(P), np.inf); face = np.full(len(P), -1, np.int64); sign = np.zeros(len(P))
    small, large, tree = index["small"], index["large"], index["tree"]
    for i, cand in enumerate(tree.query_ball_point(P, r=cap + index["radius"], workers=-1)):
        c = small[np.asarray(cand)] if len(cand) else np.empty(0, np.int64)
        if len(large): c = np.concatenate([c, large])
        if not len(c): continue
        t, ok = _moller_trumbore(P[i], D[i], v0[c], e1[c], e2[c])
        ok &= (np.abs(t) <= cap) & (np.abs(t) > 1e-12)
        if not ok.any(): continue
        j = int(np.argmin(np.where(ok, np.abs(t), np.inf)))
        dist[i] = abs(t[j]); face[i] = c[j]; sign[i] = np.sign(t[j])
    return dist, face, sign


def shoot_pairs(src_V, src_F, tgt_V, tgt_F, *, n=2000, cap_m=0.02, return_tol_m=1e-3, seed=0,
                samples=None, src_index=None, tgt_index=None, min_normal_agreement=None):
    """Symmetric normal shooting from source to target.

    Returns the kept source and target points, the mask over the samples, and why the rest were
    dropped: no hit along the source normal within cap_m, a return that landed further than
    return_tol_m from where it started, or -- when min_normal_agreement is set -- a target whose
    normal does not agree with the source's.

    On min_normal_agreement: a normal shot at a THIN sheet can cross it and hit the far side, where
    the surface faces the other way and the return test still passes because the two walls are
    roughly parallel. Measured on this body's right rib 2, those hits are 99 of 362 kept pairs, with
    a normal agreement of -0.850 and a median error of 4.52 mm against the known preimage; requiring
    agreement > 0 takes that rib from 1.215 mm mean error to 0.125 mm. It defaults to OFF because
    the gate it is judged by (docs/BODY_PARAMETERS.md) specifies shooting and the return test only.
    """
    P, N, face, bary = surface_samples(src_V, src_F, n, seed) if samples is None else samples
    tgt_index = tgt_index or build_index(tgt_V, tgt_F)
    src_index = src_index or build_index(src_V, src_F)
    d, f, s = first_hit(P, N, tgt_V, tgt_F, tgt_index, cap_m)
    hit = f >= 0
    Q = np.full_like(P, np.nan)
    Q[hit] = P[hit] + (s[hit] * d[hit])[:, None] * N[hit]
    back = np.zeros(len(P), bool); ret = np.full(len(P), np.inf)
    if hit.any():
        tgt_vn = vertex_normals(tgt_V, tgt_F)
        tri = tgt_V[tgt_F[f[hit]]]
        w = _barycentric(Q[hit], tri)
        M = np.einsum('nk,nkj->nj', w, tgt_vn[tgt_F[f[hit]]])
        M /= np.maximum(np.linalg.norm(M, axis=1, keepdims=True), 1e-30)
        d2, f2, s2 = first_hit(Q[hit], M, src_V, src_F, src_index, cap_m)
        R = Q[hit] + (s2 * d2)[:, None] * M
        good = f2 >= 0
        r = np.full(len(R), np.inf); r[good] = np.linalg.norm(R[good] - P[hit][good], axis=1)
        ret[hit] = r; back[hit] = good
    keep = hit & back & (ret <= return_tol_m)
    disagree = np.zeros(len(P), bool)
    if min_normal_agreement is not None and keep.any():
        tgt_vn = vertex_normals(tgt_V, tgt_F)
        idx = np.flatnonzero(keep)
        tri = tgt_V[tgt_F[f[idx]]]
        M = np.einsum('nk,nkj->nj', _barycentric(Q[idx], tri), tgt_vn[tgt_F[f[idx]]])
        M /= np.maximum(np.linalg.norm(M, axis=1, keepdims=True), 1e-30)
        bad = np.einsum('ij,ij->i', N[idx], M) <= min_normal_agreement
        disagree[idx[bad]] = True; keep[idx[bad]] = False
    return dict(source=P[keep], target=Q[keep], normal=N[keep], keep=keep, sampled=len(P),
                no_hit=int((~hit).sum()), no_return=int((hit & ~back).sum()),
                return_too_far=int((hit & back & (ret > return_tol_m)).sum()),
                normal_disagreed=int(disagree.sum()),
                return_distance_m=ret, all_source=P, all_target=Q, face=face, bary=bary)


def _barycentric(P, tri):
    v0 = tri[:, 1] - tri[:, 0]; v1 = tri[:, 2] - tri[:, 0]; v2 = P - tri[:, 0]
    d00 = np.einsum('ij,ij->i', v0, v0); d01 = np.einsum('ij,ij->i', v0, v1); d11 = np.einsum('ij,ij->i', v1, v1)
    d20 = np.einsum('ij,ij->i', v2, v0); d21 = np.einsum('ij,ij->i', v2, v1)
    den = np.maximum(d00 * d11 - d01 * d01, 1e-30)
    v = (d11 * d20 - d01 * d21) / den; w = (d00 * d21 - d01 * d20) / den
    return np.stack([1 - v - w, v, w], axis=1)
