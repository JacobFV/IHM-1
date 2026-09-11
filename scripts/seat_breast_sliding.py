"""Seat each registered female breast on this body's chest wall with a SLIDING base.

WHY THIS BOUNDARY CONDITION. Prescribing each penetrating base node onto its closest bed point
(scripts/seat_breast_on_chest_wall.py) is ill-posed: it collapsed 69.4% of s1159-left's base
triangles below a fifth of their area, so no solver could keep J > 0.2. The anatomical condition is
different -- the breast is attached to the pectoral fascia along the NORMAL and slides on it.

THE BOUNDARY CONDITION (docs/BODY_PARAMETERS.md, "The sliding base", fixed before this was built):
every base node that penetrates the bed in the registered position is held ON the bed along the
bed's normal, bilateral and frictionless tangentially; every other base node may not enter the bed
(unilateral); the rest of the surface is free. FEBio's counterpart is sliding contact against a
rigid bed with tension allowed on the held nodes.

THE BED is the anterior (forward-facing) surface of this body's pectoralis major -- clavicular,
sternocostal and abdominal parts, BodyParts3D only, the breast's anatomical bed, not the ribs.

WHICH NODES ARE BASE, AND WHICH SHEET IS THEIR BED (implementation, not fixed by the judge). Each
posterior-facing boundary node (outward normal z < 0) CASTS A RAY along its own outward normal, held
fixed from the registered position. A hit ANTERIOR to the node means the muscle's surface is in
front of it -- the node is inside the chest, so it is HELD; a hit POSTERIOR means the bed is behind
it, so it is UNILATERAL; no hit either way means there is no muscle along that line -- the lateral
breast over serratus anterior -- and the node is FREE, which is what "the rest of the surface is
free" means for it.

Not "the nearest bed point": the bed is three overlapping pectoralis parts, and for 34% of the
candidates on s1159-left the NEAREST sheet was on the wrong side of the node (a median 24.6 mm
away), so their measured gap flipped sign as they slid and the re-linearisation never converged --
tens of millimetres of apparent penetration through a constraint that forbids any. A ray cannot
pick the wrong sheet, and its hit point travels continuously as the node slides.

SEAT IT IN TWO STEPS (docs/BODY_PARAMETERS.md, fixed 2026-09-10 before it ran). 45 mm of overlap
is placement, not tissue deformation: the breast is another woman's tissue where a similarity
registration put it, and neither solver will push that out (FEBio 582 s without finishing a step;
the in-repo solve 0.1% of the seating per step). So step 1 PLACES the breast rigidly -- translation
and rotation, no scale, so volume cannot change -- minimising the sum of squared penetration depths
of the base nodes along their own outward normals, by L-BFGS over the six degrees of freedom. Step 2
CONFORMS it with the sliding boundary condition above, from the placed state. A breast needing more
than PLACEMENT_LIMIT_M of translation is a REGISTRATION failure for that subject, reported unseated
rather than moved until it fits.

THE SOLVERS. In-repo: SlidingRegion (ihm/assembly/sliding_contact.py), projected Newton on
DeformableRegion's energy with each base node's frame rotated to the bed normal, so the normal is
a bound; frames re-linearised as nodes slide. FEBio 4.13: ihm/assembly/febio_sliding.py, rigid
shell bed, sliding-elastic contact, augmented Lagrangian, tension on the held faces.

GATES per breast: (a) volume within 1%; (b) every tet J > 0.2; (c) the two solvers within 5% of the
maximum displacement; (d) no base triangle flipped in the solved state. Reported, not judged: the
held nodes' normal gap, tangential slide, this body's rib points inside the breast, and whether
Young's modulus cancels. A breast that fails is recorded; the boundary condition is not changed and
re-run on the same subjects.

CAVEATS: one clinical subject per breast, as a segmentation model drew it; 'breast' is one
soft-tissue label with no gland, ducts or nipple; nu = 0.49 is assumed; gravity and the unloaded
supine reference shape are out of scope; the mapped female TRUNK is unchanged, so this body's
sternum still sits outside it (the separate skin-envelope problem).
"""
import argparse, gzip, importlib.util, json, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.sliding_contact import SlidingRegion, seat_on_bed                       # noqa: E402
from ihm.assembly.prescribed_deformation import lame, tet_volumes                          # noqa: E402
from ihm.assembly.febio_crosscheck import read_msh_tets, run_febio, rcm_order              # noqa: E402
from ihm.assembly.febio_sliding import write_sliding_feb, split_base_faces                 # noqa: E402

OUT = ROOT / "data/derived/female-breast-sliding-v1"
REG = {"s0790": "female-torso-registered-v1", "s1067": "female-torso-registered-s1067-v1",
       "s1159": "female-torso-registered-s1159-v1", "s0970": "female-torso-registered-s0970-v1"}
# THE MUSCULAR CHEST WALL (docs/BODY_PARAMETERS.md, corrected 2026-09-10): pectoralis major's three
# parts, pectoralis minor, serratus anterior, the external oblique and rectus abdominis -- and still
# NOT the ribs. Pectoralis major alone is not the bed under 13-29% of a breast; those nodes were
# measured against the muscle's EDGE. Rectus abdominis has no BodyParts3D entity, so its Z-Anatomy
# mesh is the only source for it (the duplicate-copy caution elsewhere does not apply: there is
# nothing to duplicate).
BED_MUSCLES = {"right": ["body-bp3d-FJ1447", "body-bp3d-FJ1464", "body-bp3d-FJ1446", "body-bp3d-FJ1456",
                         "body-bp3d-FJ1459", "body-bp3d-FJ1452", "body-za-f8e80e0f249b354b"],
               "left": ["body-bp3d-FJ1447M", "body-bp3d-FJ1464M", "body-bp3d-FJ1446M", "body-bp3d-FJ1456M",
                        "body-bp3d-FJ1459M", "body-bp3d-FJ1452M", "body-za-a45e813cd3661ae2"]}
FTETWILD = ROOT / "data/runtime/tolerant-mesher/build/FloatTetwild_bin"
NU, E_PA, E_CHECK_PA = 0.49, 1000.0, 10000.0
DECIMATE_FACES, DECIMATE_VOLUME_TOL, MESH_LR = 8000, 0.005, 0.08
# A ray finds the bed only NEARBY: uncapped, a ray grazes off to a distant fold of the muscle and
# calls it the bed (175 mm of 'penetration' on a breast 13 cm deep). 60 mm is well beyond the
# deepest real penetration measured on these four subjects (43 mm).
# A ray that meets the bed at a shallow angle has no well-defined bed direction: a micrometre of
# motion slides its hit point millimetres along the surface, which reads as the constraint moving.
# Those nodes are free, like the ones with no muscle on their line.
# Both convergence tolerances are set to the BED's facet scale, ~1 mm, and declared together: the
# muscular wall is a faceted anatomical mesh whose normals turn several degrees over 2-4 mm, so a
# sliding node's hit point moves ~0.5 mm per pass and its measured gap wanders by ~0.2 mm. Asking
# either to settle to 0.05 mm asks the solver to resolve the bed finer than the bed is defined.
# These are convergence criteria, not gates: (a)-(d) are unchanged and the contact gap is reported.
RAY_REACH_M, GRAZING_COS, LOAD_STEPS, GAP_TOL_M = 0.060, 0.3, 8, 1e-3
# The re-linearisation is called self-consistent when the association stops moving. 0.1 mm was
# below the BED's own facet scale -- the muscular wall's normals turn several degrees over
# 2-4 mm, so a sliding node's hit point genuinely moves ~0.5 mm per pass and the criterion
# could never be met. It is set to 1 mm, the scale of the surface the tissue slides on, and
# declared here rather than tuned quietly.
ASSOC_TOL_M = 1e-3
CONTROL_JUMP_LIMIT_M = 5e-4   # raising this to 5e-3 made R-prime worse, not better: see the doc
# The objective is penetration only, so it has a degenerate minimum: fly the breast away and every
# ray misses the muscle. Unbounded L-BFGS took it in one step (999.7 mm). The search is therefore
# bounded to the region where 'placement' means anything -- the 25 mm honesty limit itself, and
# 15 deg -- so a breast that wants more sits at the bound and is reported a registration failure.
# The finite-difference step is 0.5 mm, not 0.1: a ray-cast objective changes in facet-sized jumps
# and a finer step reads as noise (the search stalled after 2 iterations). Several starts are tried,
# including pure anterior offsets -- lifting the breast off the chest is the obvious direction -- and
# the best is kept.
# THE TRIM (docs/BODY_PARAMETERS.md, fourth attempt, fixed before it ran): breast tissue does not
# lie behind pectoralis major, so a label that does is a segmentation error -- these labels come
# from a model run on a clinical scan. Every boundary vertex more than TRIM_DEPTH_M behind the
# muscular wall along its own outward ray is removed with the tets it belongs to, and the rest is
# re-meshed. 20 mm sits below the measured p99 of 24.7 mm, so the trim is a real test. If it takes
# more than TRIM_VOLUME_TOL of the breast, the label is not locally wrong -- the breast is in the
# wrong place -- and the subject is a registration failure with no seating attempted.
TRIM_DEPTH_M, TRIM_VOLUME_TOL = 0.020, 0.01
TRIM_ENABLED = False
PLACEMENT_LIMIT_M, PLACEMENT_SAMPLE = 0.025, 1200
PLACEMENT_CAP_M, PLACEMENT_TURN_CAP_RAD = 0.100, 0.524
PLACEMENT_EPS_M, PLACEMENT_STARTS_MM = 5e-4, (0.0, 8.0, 16.0, 22.0)
CAVEATS = ["one clinical subject per breast, as a segmentation model drew it",
           "'breast' is a single soft-tissue label: no gland, ducts or nipple",
           "nu = 0.49 is assumed (adipose is nearly incompressible)",
           "gravity and the unloaded supine reference shape are out of scope",
           "the mapped female TRUNK is not changed; this body's sternum still sits outside it (skin-envelope problem)"]


def say(*a): print(*a, flush=True)


def read_obj(p):
    V, F = [], []
    for l in open(p):
        if l.startswith("v "): V.append([float(x) for x in l.split()[1:4]])
        elif l.startswith("f "): F.append([int(x.split("/")[0]) - 1 for x in l.split()[1:4]])
    return np.asarray(V), np.asarray(F, np.int64)


def write_obj(p, V, F, header):
    with open(p, "w") as h:
        for line in header: h.write(f"# {line}\n")
        for x in V: h.write("v %.9g %.9g %.9g\n" % tuple(x))
        for t in F: h.write("f %d %d %d\n" % (t[0] + 1, t[1] + 1, t[2] + 1))


def boundary_faces(T):
    f = np.vstack([T[:, [1, 2, 3]], T[:, [0, 3, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 1]]])
    _, idx, cnt = np.unique(np.sort(f, 1), axis=0, return_index=True, return_counts=True)
    return f[idx[cnt == 1]]


def face_normals(V, F):
    n = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    return n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-30)


def bed(side, near=None, margin_m=0.08):
    """The muscular chest wall on this side: whole muscle surfaces, NOT filtered to anterior-facing
    faces. The filter existed because a closest-point association needs a single-valued patch; a ray
    hits the superficial surface first by construction, and whole surfaces have no rim holes to
    catch nodes. It also matters for the added muscles: serratus and the obliques face laterally, so
    an anterior-only filter would discard most of the coverage this correction adds. `near` crops
    the bed to that point cloud's bounding box plus margin_m, which is what makes ray casting
    affordable."""
    ents = {e["id"]: e for e in json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]}
    Vs, Fs, off = [], [], 0
    for i in BED_MUSCLES[side]:
        g = json.loads(gzip.decompress((ROOT / ents[i]["reference_geometry"]["path"]).read_bytes()))
        V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
        Vs.append(V); Fs.append(F + off); off += len(V)
    V = np.vstack(Vs); F = np.vstack(Fs)
    if near is not None:
        lo = np.asarray(near).min(0) - margin_m; hi = np.asarray(near).max(0) + margin_m
        keep = ((V[F] >= lo).all(2) & (V[F] <= hi).all(2)).any(1)
        F = F[keep]
    used = np.unique(F); remap = -np.ones(len(V), np.int64); remap[used] = np.arange(len(used))
    V, F = V[used], remap[F]
    e = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), 1)
    uniq, cnt = np.unique(e, axis=0, return_counts=True)
    return V, F, np.unique(uniq[cnt == 1])                      # rim vertices of the open patch


def ray_hits(P, D, V, F, chunk=192):
    """Nearest hit of the rays P + t D (t > 0) on the triangles (V, F): distances and face indices."""
    v0 = V[F[:, 0]]; e1 = V[F[:, 1]] - v0; e2 = V[F[:, 2]] - v0
    dist = np.full(len(P), np.inf); face = np.full(len(P), -1, np.int64)
    for s in range(0, len(P), chunk):
        p, d = P[s:s + chunk], D[s:s + chunk]
        h = np.cross(d[:, None, :], e2[None, :, :])
        a = np.einsum('fj,mfj->mf', e1, h); usable = np.abs(a) > 1e-16
        inv = 1.0 / np.where(usable, a, 1.0)
        sv = p[:, None, :] - v0[None, :, :]
        u = inv * np.einsum('mfj,mfj->mf', sv, h)
        q = np.cross(sv, e1[None, :, :])
        v = inv * np.einsum('mj,mfj->mf', d, q)
        t = inv * np.einsum('fj,mfj->mf', e2, q)
        ok = usable & (u >= -1e-9) & (v >= -1e-9) & (u + v <= 1 + 1e-9) & (t > 1e-9)
        t = np.where(ok, t, np.inf)
        j = np.argmin(t, axis=1); tm = t[np.arange(len(p)), j]
        dist[s:s + chunk] = tm; face[s:s + chunk] = np.where(np.isfinite(tm), j, -1)
    return dist, face


def build_bed_index(V, F):
    """Two tiers, because one oversized triangle would set the query radius for every ray: the bed's
    face radii run 1.6 mm median but 28.9 mm max, and querying with the max pulls 3,079 candidates
    per ray instead of 374. Small faces go in a KD-tree; the 1% oversized ones are tested in bulk."""
    cent = V[F].mean(1); rad = np.linalg.norm(V[F] - cent[:, None, :], axis=2).max(1)
    cut = float(np.percentile(rad, 99))
    small = np.flatnonzero(rad <= cut); large = np.flatnonzero(rad > cut)
    return dict(tree=cKDTree(cent[small]), small=small, large=large, radius=cut)


def ray_hits_indexed(P, D, V, F, index, reach_m):
    """Same answer as ray_hits, but each ray is tested only against faces near its own segment: the
    bed of the whole muscular chest wall is 125k faces, and testing all of them for every ray costs
    2.45 s per 200 nodes. Any face the segment can hit has its centroid within reach/2 + the face
    radius of the segment's midpoint."""
    small, large, tree = index["small"], index["large"], index["tree"]
    v0 = V[F[:, 0]]; e1 = V[F[:, 1]] - v0; e2 = V[F[:, 2]] - v0
    dist = np.full(len(P), np.inf); face = np.full(len(P), -1, np.int64)
    if len(large):                                   # the oversized few, against every ray at once
        tl, fl = ray_hits(P, D, V, F[large])
        take = (tl <= reach_m) & (tl < dist)
        dist[take] = tl[take]; face[take] = large[fl[take]]
    mid = P + D * (reach_m / 2)
    for i, cand in enumerate(tree.query_ball_point(mid, r=reach_m / 2 + index["radius"], workers=-1)):
        if not cand: continue
        c = small[np.asarray(cand)]
        h = np.cross(D[i], e2[c]); a = np.einsum('fj,fj->f', e1[c], h)
        usable = np.abs(a) > 1e-16
        inv = 1.0 / np.where(usable, a, 1.0)
        sv = P[i] - v0[c]
        u = inv * np.einsum('fj,fj->f', sv, h)
        q = np.cross(sv, e1[c])
        v = inv * (q @ D[i])
        t = inv * np.einsum('fj,fj->f', e2[c], q)
        ok = usable & (u >= -1e-9) & (v >= -1e-9) & (u + v <= 1 + 1e-9) & (t > 1e-9) & (t <= reach_m)
        if not ok.any(): continue
        t = np.where(ok, t, np.inf); j = int(np.argmin(t))
        if t[j] < dist[i]: dist[i] = t[j]; face[i] = c[j]
    return dist, face


def bed_rays(V, F, directions, reach_m=0.060):
    """association(points) -> (c, n): where each node's fixed ray meets the bed, and the bed's
    anterior-facing normal there. The ray is cast both ways and the nearer hit wins, so a node
    inside the chest associates with the sheet in FRONT of it and one outside with the sheet
    behind it -- a classification the nearest-point rule gets wrong wherever the muscle's three
    parts overlap."""
    fn = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    # SMOOTH vertex normals, interpolated at the hit point. The bed is a marching-cubes surface, so
    # its facet normals are noisy: neighbouring base nodes given their own raw facet normal are
    # pushed in measurably different directions, which distorts the elements between them and
    # inverts them even at micrometre steps.
    vn = np.zeros_like(V)
    for k in range(3): np.add.at(vn, F[:, k], fn)
    vn /= np.maximum(np.linalg.norm(vn, axis=1, keepdims=True), 1e-30)
    fn /= np.maximum(np.linalg.norm(fn, axis=1, keepdims=True), 1e-30)
    D = np.asarray(directions, float)
    D = D / np.maximum(np.linalg.norm(D, axis=1, keepdims=True), 1e-30)
    index = build_bed_index(V, F)
    def association(P):
        P = np.ascontiguousarray(np.asarray(P, float))
        t_back, f_back = ray_hits_indexed(P, D, V, F, index, reach_m)     # bed behind
        t_front, f_front = ray_hits_indexed(P, -D, V, F, index, reach_m)  # bed in front
        front = t_front < t_back
        t = np.where(front, t_front, t_back); f = np.where(front, f_front, f_back)
        hit = np.isfinite(t)
        c = P + np.where(front, -1.0, 1.0)[:, None] * np.where(hit, t, 0.0)[:, None] * D
        tri = V[F[np.where(f >= 0, f, 0)]]
        w = np.stack([np.linalg.norm(np.cross(tri[:, 1] - c, tri[:, 2] - c), axis=1),
                      np.linalg.norm(np.cross(tri[:, 2] - c, tri[:, 0] - c), axis=1),
                      np.linalg.norm(np.cross(tri[:, 0] - c, tri[:, 1] - c), axis=1)], 1)
        w /= np.maximum(w.sum(1, keepdims=True), 1e-30)
        n = np.einsum('ij,ijk->ik', w, vn[F[np.where(f >= 0, f, 0)]])
        n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-30)
        n[np.einsum('ij,ij->i', n, D) > 0] *= -1        # orient anteriorly, opposite the outward normal
        # No hit within reach: INVALID, not a plane through the node. Returning the node's own
        # position reads as a zero gap, so a held node 14 mm deep reports its error as its own depth
        # and the re-linearisation never converges. The caller keeps the node's previous association.
        c[~hit] = np.nan; n[~hit] = np.nan
        return c, n
    def has_bed(P):
        P = np.ascontiguousarray(np.asarray(P, float))
        return (ray_hits_indexed(P, D, V, F, index, reach_m)[0] <= reach_m) | \
               (ray_hits_indexed(P, -D, V, F, index, reach_m)[0] <= reach_m)
    return association, has_bed


def largest_component(T):
    """Keep only the largest tet component: trimming can leave islands."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    f = np.vstack([T[:, [1, 2, 3]], T[:, [0, 3, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 1]]])
    _, inv = np.unique(np.sort(f, 1), axis=0, return_inverse=True)
    tet_of = np.tile(np.arange(len(T)), 4)
    order = np.argsort(inv, kind="stable"); inv_s, tet_s = inv[order], tet_of[order]
    shared = np.flatnonzero(np.bincount(inv) == 2)
    if not len(shared): return np.ones(len(T), bool)
    start = np.searchsorted(inv_s, shared)
    A = coo_matrix((np.ones(len(shared)), (tet_s[start], tet_s[start + 1])), shape=(len(T), len(T)))
    ncomp, label = connected_components(A, directed=False)
    return np.ones(len(T), bool) if ncomp == 1 else label == np.argmax(np.bincount(label))


def penetration_of(X, B, bV, bF, nodes=None):
    """Depth behind the muscular wall along each boundary vertex's own outward ray."""
    vn = np.zeros_like(X)
    for k in range(3): np.add.at(vn, B[:, k], np.cross(X[B[:, 1]] - X[B[:, 0]], X[B[:, 2]] - X[B[:, 0]]))
    idx = np.unique(B) if nodes is None else np.asarray(nodes)
    dirs = vn[idx] / np.maximum(np.linalg.norm(vn[idx], axis=1, keepdims=True), 1e-30)
    association, has = bed_rays(bV, bF, dirs, RAY_REACH_M)
    c, n = association(X[idx])
    gap = np.einsum('ij,ij->i', n, X[idx] - c)
    return idx, np.where(np.isfinite(gap), np.maximum(0.0, -gap), 0.0), dirs, has


def stage_prepare(sid, side, d):
    import igl
    src = ROOT / "data/derived" / REG[sid] / f"breast_{side}.obj"
    md = d / f"mesh_dec{DECIMATE_FACES}"; md.mkdir(parents=True, exist_ok=True)
    Vo, Fo = read_obj(src); to = Vo[Fo]
    vol_o = abs(np.einsum("ij,ij->i", to[:, 0], np.cross(to[:, 1], to[:, 2])).sum() / 6)
    dec, msh = md / "decimated.obj", md / "out.msh"
    if not dec.exists():
        out = [np.asarray(o) for o in igl.qslim(np.ascontiguousarray(Vo), np.ascontiguousarray(Fo), DECIMATE_FACES) if hasattr(o, "shape")]
        U = next(o for o in out if o.dtype.kind == "f" and o.ndim == 2 and o.shape[1] == 3)
        G = next(o for o in out if o.dtype.kind in "iu" and o.ndim == 2 and o.shape[1] == 3)
        tu = U[G]; vol_d = abs(np.einsum("ij,ij->i", tu[:, 0], np.cross(tu[:, 1], tu[:, 2])).sum() / 6)
        change = vol_d / vol_o - 1
        if abs(change) > DECIMATE_VOLUME_TOL: raise SystemExit(f"decimation changed the volume by {100*change:+.3f}%")
        write_obj(dec, U, G, [f"{sid} {side} breast, qslim to {len(G)} faces ({100*change:+.3f}% volume)"])
        print(f"  decimated {len(Fo)} -> {len(G)} faces ({100*change:+.3f}% volume)")
    if not msh.exists():
        cmd = [str(FTETWILD), "-i", str(dec), "-o", str(msh), "--no-binary", "-e", "1e-3", "-l", str(MESH_LR), "--max-threads", "12"]
        t0 = time.time(); p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        (md / "ftetwild.log").write_text(p.stdout[-4000:] + p.stderr[-4000:])
        if p.returncode != 0 or not msh.exists(): raise SystemExit(f"fTetWild failed on {sid} {side}")
        print(f"  fTetWild {time.time()-t0:.0f} s")
    X, T = read_msh_tets(msh)
    used = np.unique(T); remap = -np.ones(len(X), np.int64); remap[used] = np.arange(len(used)); X, T = X[used], remap[T]
    v = tet_volumes(X, T); T[v < 0] = T[v < 0][:, [0, 2, 1, 3]]
    B = boundary_faces(T)
    bV, bF, rim = bed(side, near=X)
    # --- the trim, gated before any seating (the fourth attempt; superseded) ---
    # The label trim answered "is the tissue behind the wall mislabelled or misplaced". It answered
    # MISPLACED on all eight breasts (2.0-7.7% removed against a 1% gate), and the registration line
    # that followed has now ended: a rib cage is struts with air between them, and pairing, closing
    # and field-reading all failed on that one cause. The registration error is therefore ABSORBED by
    # the seating rather than corrected before it, so the trim is off unless asked for.
    if not TRIM_ENABLED:
        trim = dict(enabled=False, note="registration error absorbed by the seating, not trimmed")
        (d / "trim.json").write_text(json.dumps(trim, indent=2) + "\n")
        nb = face_normals(X, B)
    else:
        idx, pen, _, _ = penetration_of(X, B, bV, bF)
        deep = idx[pen > TRIM_DEPTH_M]
        volume_all = tet_volumes(X, T).sum()
        trim = dict(depth_mm=TRIM_DEPTH_M * 1e3, deep_vertices=int(len(deep)),
                    penetration_before_mm=dict(max=float(pen.max() * 1e3), p99=float(np.percentile(pen, 99) * 1e3)))
        if len(deep):
            isdeep = np.zeros(len(X), bool); isdeep[deep] = True
            keep = ~isdeep[T].any(1)
            keep &= largest_component(T[keep])[np.cumsum(keep) - 1] if keep.any() else keep
            removed = 1.0 - tet_volumes(X, T[keep]).sum() / volume_all
            trim.update(tets_removed=int((~keep).sum()), removed_volume_fraction=float(removed))
            print(f"  trim: {len(deep)} vertices deeper than {TRIM_DEPTH_M*1e3:.0f} mm (max {pen.max()*1e3:.1f} mm), "
                  f"{int((~keep).sum())} tets, {100*removed:.3f}% of the breast's volume")
            if removed > TRIM_VOLUME_TOL:
                (d / "trim.json").write_text(json.dumps({**trim, "registration_failure": True}, indent=2) + "\n")
                raise SystemExit(f"REGISTRATION FAILURE for {sid} {side}: the trim takes {100*removed:.2f}% > "
                                 f"{100*TRIM_VOLUME_TOL:.0f}%; the breast is misplaced, not mislabelled; no seating attempted")
            # re-mesh the trimmed body
            Bk = boundary_faces(T[keep]); kv = np.unique(Bk); kmap = -np.ones(len(X), np.int64); kmap[kv] = np.arange(len(kv))
            write_obj(md / "trimmed.obj", X[kv], kmap[Bk], [f"{sid} {side} breast, trimmed {100*removed:.3f}% behind the muscular wall"])
            tmsh = md / "trimmed.msh"
            if not tmsh.exists():
                cmd = [str(FTETWILD), "-i", str(md / "trimmed.obj"), "-o", str(tmsh), "--no-binary", "-e", "1e-3",
                       "-l", str(MESH_LR), "--max-threads", "12"]
                t0 = time.time(); pr = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
                (md / "ftetwild_trimmed.log").write_text(pr.stdout[-4000:] + pr.stderr[-4000:])
                if pr.returncode != 0 or not tmsh.exists(): raise SystemExit(f"fTetWild failed on the trimmed {sid} {side}")
                print(f"  re-meshed the trimmed breast in {time.time()-t0:.0f} s")
            X, T = read_msh_tets(tmsh)
            used = np.unique(T); remap = -np.ones(len(X), np.int64); remap[used] = np.arange(len(used)); X, T = X[used], remap[T]
            v = tet_volumes(X, T); T[v < 0] = T[v < 0][:, [0, 2, 1, 3]]
            B = boundary_faces(T)
            trim["volume_after_remesh_ml"] = float(tet_volumes(X, T).sum() * 1e6)
            trim["removed_volume_fraction_after_remesh"] = float(1.0 - tet_volumes(X, T).sum() / volume_all)
            _, pen2, _, _ = penetration_of(X, B, bV, bF)
            trim["penetration_after_mm"] = dict(max=float(pen2.max() * 1e3), p99=float(np.percentile(pen2, 99) * 1e3))
            print(f"  after the trim and re-mesh: {len(X)} nodes, {len(T)} tets, penetration max {pen2.max()*1e3:.1f} mm")
    (d / "trim.json").write_text(json.dumps(trim, indent=2) + "\n")
    nb = face_normals(X, B)
    posterior = np.unique(B[nb[:, 2] < 0]); anterior = np.unique(B[nb[:, 2] > 0])
    vn = np.zeros_like(X)
    for k in range(3): np.add.at(vn, B[:, k], np.cross(X[B[:, 1]] - X[B[:, 0]], X[B[:, 2]] - X[B[:, 0]]))
    dirs = vn[posterior] / np.maximum(np.linalg.norm(vn[posterior], axis=1, keepdims=True), 1e-30)
    association, has_bed = bed_rays(bV, bF, dirs, RAY_REACH_M)
    c, n = association(X[posterior])
    facing = np.abs(np.einsum('ij,ij->i', np.nan_to_num(n), dirs)) >= GRAZING_COS
    on_line = has_bed(X[posterior]) & facing
    base = posterior[on_line]; dirs = dirs[on_line]
    gap0 = np.einsum('ij,ij->i', n[on_line], X[base] - c[on_line])
    np.savez(d / "prepared.npz", X=X, T=T, base=base, posterior=posterior, anterior=anterior, boundary=B,
             gap0=gap0, bedV=bV, bedF=bF, rim=rim, ray_directions=dirs)
    info = dict(nodes=len(X), tets=len(T), mesh_ml=tet_volumes(X, T).sum() * 1e6, surface_ml=vol_o * 1e6,
                posterior_nodes=len(posterior), base_nodes=len(base), held=int((gap0 < 0).sum()),
                unilateral=int((gap0 >= 0).sum()), free_no_bed_on_the_line=int((~on_line).sum()), grazing_excluded=int((~facing).sum()),
                deepest_penetration_mm=float(-gap0.min() * 1e3), trim=trim)
    (d / "prepared.json").write_text(json.dumps(info, indent=2) + "\n")
    print(f"  {len(X)} nodes, {len(T)} tets; base {len(base)} of {len(posterior)} posterior nodes "
          f"({int((~on_line).sum())} free: no muscle along their ray): {info['held']} held, {info['unilateral']} unilateral; "
          f"deepest penetration {info['deepest_penetration_mm']:.1f} mm")


def stage_place(sid, side, d):
    """Step 1: the rigid motion (no scale) minimising the summed squared penetration of the base
    nodes along their own outward normals. Rewrites the prepared state in the PLACED pose, keeping
    the registered one beside it."""
    from scipy.optimize import minimize
    from scipy.spatial.transform import Rotation
    P = np.load(d / "prepared.npz")
    X0 = P["X_registered"] if "X_registered" in P.files else P["X"]
    T, B, posterior, bV, bF = P["T"], P["boundary"], P["posterior"], P["bedV"], P["bedF"]
    vn = np.zeros_like(X0)
    for k in range(3): np.add.at(vn, B[:, k], np.cross(X0[B[:, 1]] - X0[B[:, 0]], X0[B[:, 2]] - X0[B[:, 0]]))
    dirs0 = vn[posterior] / np.maximum(np.linalg.norm(vn[posterior], axis=1, keepdims=True), 1e-30)
    centroid = X0.mean(0)

    # Translation and rotation MAGNITUDES are what the rules are about, so the six free parameters
    # are squashed onto magnitude caps rather than boxed per component: |t| < PLACEMENT_CAP_M and
    # |rotvec| < PLACEMENT_TURN_CAP_RAD, smoothly, so L-BFGS stays unconstrained. Boxing components
    # let |t| reach 43 mm under a 25 mm rule; the 25 mm rule then judges the result.
    squash = lambda v, cap: cap * v / np.sqrt(1.0 + float(v @ v))

    def place(params):
        t = squash(params[:3], PLACEMENT_CAP_M); rv = squash(params[3:], PLACEMENT_TURN_CAP_RAD)
        R = Rotation.from_rotvec(rv).as_matrix()
        return (X0 - centroid) @ R.T + centroid + t, dirs0 @ R.T, t, rv

    def raw_penetration(params, idx):
        Y, D, _, _ = place(params)
        association, _ = bed_rays(bV, bF, D[idx], RAY_REACH_M)
        c, n = association(Y[posterior[idx]])
        gap = np.einsum('ij,ij->i', n, Y[posterior[idx]] - c)
        return np.where(np.isfinite(gap), np.maximum(0.0, -gap), np.nan)

    pen0 = np.nan_to_num(raw_penetration(np.zeros(6), np.arange(len(posterior))))

    def penetration(params, idx):
        """A node that LOSES its bed is counted at its PRE-placement penetration, so carrying tissue
        off the muscle cannot pay -- the degenerate minimum of a penetration-only objective, which
        flew the breast 999.7 mm unbounded and to the corner of the box bounded."""
        pen = raw_penetration(params, idx)
        return np.where(np.isfinite(pen), pen, pen0[idx])

    sample = np.random.default_rng(0).choice(len(posterior), min(PLACEMENT_SAMPLE, len(posterior)), replace=False)
    before = pen0.copy()
    t0 = time.time()
    objective = lambda q: float((penetration(q, sample) ** 2).sum())
    anterior = -dirs0.mean(0); anterior /= max(np.linalg.norm(anterior), 1e-30)
    tries = []
    for mm in PLACEMENT_STARTS_MM:
        q0 = np.zeros(6); q0[:3] = anterior * (mm * 1e-3) / PLACEMENT_CAP_M   # in squashed coordinates
        ri = minimize(objective, q0, method="L-BFGS-B",
                      options={"eps": PLACEMENT_EPS_M, "maxiter": 60, "ftol": 1e-14, "gtol": 1e-14})
        tries.append(ri)
        print(f"    start {mm:.0f} mm anterior: objective {objective(q0):.3e} -> {ri.fun:.3e}, "
              f"|t| {np.linalg.norm(place(ri.x)[2])*1e3:.2f} mm, {ri.nit} iterations", flush=True)
    r = min(tries, key=lambda ri: ri.fun)
    after = penetration(r.x, np.arange(len(posterior)))

    def with_bed(params):
        Y, D, _, _ = place(params); _, has = bed_rays(bV, bF, D, RAY_REACH_M); return int(has(Y[posterior]).sum())
    bed_before, bed_after = with_bed(np.zeros(6)), with_bed(r.x)
    _, _, t_opt, rv_opt = place(r.x)
    shift = float(np.linalg.norm(t_opt)); turn = float(np.degrees(np.linalg.norm(rv_opt)))
    rec = dict(translation_mm=(t_opt * 1e3).tolist(), translation_magnitude_mm=shift * 1e3,
               rotation_deg=turn, rotation_axis=(rv_opt / max(np.linalg.norm(rv_opt), 1e-30)).tolist(),
               penetration_before_mm=dict(max=float(before.max() * 1e3), median_over_penetrating=float(np.median(before[before > 0]) * 1e3) if (before > 0).any() else 0.0,
                                          nodes=int((before > 0).sum())),
               penetration_after_mm=dict(max=float(after.max() * 1e3), median_over_penetrating=float(np.median(after[after > 0]) * 1e3) if (after > 0).any() else 0.0,
                                         nodes=int((after > 0).sum())),
               objective_before=float((before ** 2).sum()), objective_after=float((after ** 2).sum()),
               iterations=int(r.nit), starts_mm=list(PLACEMENT_STARTS_MM),
               start_objectives=[float(ri.fun) for ri in tries], wall_seconds=time.time() - t0,
               nodes_with_bed_on_their_ray=dict(before=bed_before, after=bed_after),
               registration_failure=bool(shift > PLACEMENT_LIMIT_M),
               limit_mm=PLACEMENT_LIMIT_M * 1e3)
    (d / "place.json").write_text(json.dumps(rec, indent=2) + "\n")
    print(f"  rigid placement: translation {shift*1e3:.2f} mm {np.round(t_opt*1e3,2).tolist()}, rotation {turn:.2f} deg, "
          f"{r.nit} iterations, {time.time()-t0:.0f} s")
    print(f"  penetration: max {before.max()*1e3:.1f} -> {after.max()*1e3:.1f} mm; nodes penetrating {int((before>0).sum())} -> {int((after>0).sum())}; "
          f"nodes with muscle on their ray {bed_before} -> {bed_after}")
    if rec["registration_failure"]:
        raise SystemExit(f"REGISTRATION FAILURE for {sid} {side}: placement wants {shift*1e3:.1f} mm "
                         f"(bound {PLACEMENT_LIMIT_M*1e3:.0f} mm); reported unseated")
    # re-classify in the placed pose: step 2 starts here
    X, D, t_best, rv_best = place(r.x)
    association, has_bed = bed_rays(bV, bF, D, RAY_REACH_M)
    c, n = association(X[posterior])
    facing = np.abs(np.einsum('ij,ij->i', np.nan_to_num(n), D)) >= GRAZING_COS
    on_line = has_bed(X[posterior]) & facing
    base = posterior[on_line]; dirs = D[on_line]
    gap0 = np.einsum('ij,ij->i', n[on_line], X[base] - c[on_line])
    np.savez(d / "prepared.npz", X=X, X_registered=X0, T=T, base=base, posterior=posterior, anterior=P["anterior"],
             boundary=B, gap0=gap0, bedV=bV, bedF=bF, rim=P["rim"], ray_directions=dirs,
             rigid_translation=r.x[:3], rigid_rotvec=r.x[3:], rigid_centroid=centroid)
    print(f"  placed: base {len(base)} of {len(posterior)} posterior nodes, {int((gap0 < 0).sum())} held, "
          f"{int((gap0 >= 0).sum())} unilateral; deepest penetration {-gap0.min()*1e3:.1f} mm")


def stage_dr(sid, side, d, young, tag=""):
    P = np.load(d / "prepared.npz")
    mu, lam = lame(young, NU)
    region = SlidingRegion(P["X"], P["T"], mu_pa=mu, lambda_pa=lam, density_kg_m3=950.0)
    closest, _ = bed_rays(P["bedV"], P["bedF"], P["ray_directions"], RAY_REACH_M)
    t0 = time.time()
    smoothed = d / "smoothed_move.npy"
    move = np.load(smoothed) if smoothed.exists() else None
    if move is not None:
        print(f"  seating with the SMOOTHED depth field (declared modelling choice): "
              f"median travel {1000*np.median(move[move>0]):.2f} mm", flush=True)
    r = seat_on_bed(region, P["base"], closest, load_steps=LOAD_STEPS, gap_tol_m=GAP_TOL_M, jump_limit_m=5e-4, move=move, assoc_tol_m=ASSOC_TOL_M,
                    log=lambda m, flush=True: print(m, flush=True))
    np.save(d / f"dr_displacement{tag}.npy", r["displacement"])
    (d / f"dr{tag}.json").write_text(json.dumps(dict(
        steps=r["steps"], cutbacks=r["cutbacks"], minimum_jacobian=r["minimum_jacobian"], converged=bool(r["converged"]),
        held=int(r["held"].sum()), young_pa=young, wall_seconds=time.time() - t0), indent=2) + "\n")
    print(f"  in-repo E={young:g}: min J {r['minimum_jacobian']:.4f}, held gap max "
          f"{np.abs(r['gap_m'][r['held']]).max()*1e3:.4f} mm, {time.time()-t0:.0f} s")


def stage_febio(sid, side, d):
    P = np.load(d / "prepared.npz"); X, T, B = P["X"], P["T"], P["boundary"]
    base, gap0 = P["base"], P["gap0"]
    hf, ff, oneway = split_base_faces(B, base, base[gap0 < 0])
    fd = d / "febio"; fd.mkdir(exist_ok=True)
    # the settings the cylinder gate was re-passed with: a soft start eases tens of millimetres of
    # initial penetration out, a stiff finish leaves a small gap
    write_sliding_feb(fd / "seat.feb", X, T, hf, ff, P["bedV"], P["bedF"], young_pa=E_PA, nu=NU,
                      steps=20, search_radius=float(-gap0.min() * 2 + 0.01), reorder=rcm_order(len(X), T),
                      penalty=50.0, ramp_from=0.001)
    print(f"  FEBio: {len(hf)} held faces, {len(ff)} unilateral faces, {len(oneway)} held nodes one-way only")
    t0 = time.time(); r = run_febio(fd / "seat.feb", threads=2, timeout=86400)
    np.save(d / "febio_displacement.npy", r["displacement"][:len(X)])
    (d / "febio.json").write_text(json.dumps({**{k: v for k, v in r.items() if k != "displacement"},
                                              "held_faces": len(hf), "unilateral_faces": len(ff),
                                              "held_nodes_one_way_only": len(oneway),
                                              "wall_seconds": time.time() - t0}, indent=2) + "\n")
    print(f"  FEBio: final time {r['final_time']:.3f}, normal termination {r['normal_termination']}, {time.time()-t0:.0f} s")


def ribs_inside(side, V, F):
    spec = importlib.util.spec_from_file_location("cwo", ROOT / "scripts/measure_female_chest_wall_offset.py")
    M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
    ents = {e["name"]: e for e in json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
            if e["role"] == "rigid_bone" and e["id"].startswith("body-bp3d-")}
    pts = []
    for o in ("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth"):
        g = json.loads(gzip.decompress((ROOT / ents[f"{side} {o} rib"]["reference_geometry"]["path"]).read_bytes()))
        p = np.asarray(g["positions"], float).reshape(-1, 3)
        pts.append(p[np.random.default_rng(0).choice(len(p), min(len(p), 1500), replace=False)])
    return float(M.inside(V, F, np.vstack(pts)).mean())


def stage_judge(sid, side, d):
    P = np.load(d / "prepared.npz"); X, T, B, base, gap0 = P["X"], P["T"], P["boundary"], P["base"], P["gap0"]
    u = np.load(d / "dr_displacement.npy"); Y = X + u
    fe = np.load(d / "febio_displacement.npy") if (d / "febio_displacement.npy").exists() else None
    J = np.linalg.det(np.swapaxes(Y[T[:, 1:]] - Y[T[:, 0, None]], 1, 2) @ np.linalg.inv(np.swapaxes(X[T[:, 1:]] - X[T[:, 0, None]], 1, 2)))
    vr = float(tet_volumes(Y, T).sum() / tet_volumes(X, T).sum())
    isb = np.zeros(len(X), bool); isb[base] = True
    tri = B[isb[B].all(1)]
    n0 = np.cross(X[tri[:, 1]] - X[tri[:, 0]], X[tri[:, 2]] - X[tri[:, 0]])
    n1 = np.cross(Y[tri[:, 1]] - Y[tri[:, 0]], Y[tri[:, 2]] - Y[tri[:, 0]])
    flipped = int((np.einsum('ij,ij->i', n0, n1) < 0).sum())
    umax = float(np.linalg.norm(u, axis=1).max())
    rms = float(np.sqrt(np.mean(np.sum((u - fe) ** 2, 1)))) if fe is not None else float("nan")
    gates = {"a_volume_within_1pct": abs(vr - 1) <= 0.01, "b_every_J_above_0_2": bool(J.min() > 0.2),
             "c_solvers_within_5pct": bool(fe is not None and rms / umax <= 0.05), "d_no_flipped_base_triangle": flipped == 0}
    closest, _ = bed_rays(P["bedV"], P["bedF"], P["ray_directions"], RAY_REACH_M)
    c, n = closest(Y[base]); gap = np.einsum('ij,ij->i', n, Y[base] - c)
    held = gap0 < 0
    ub = u[base]; slide = np.linalg.norm(ub - np.einsum('ij,ij->i', ub, n)[:, None] * n, axis=1)
    ant = P["anterior"]; ua = np.linalg.norm(u[ant], axis=1)
    rec = dict(subject=sid, side=side, nodes=int(len(X)), tets=int(len(T)), volume_ratio=vr, min_J=float(J.min()),
               max_displacement_mm=umax * 1e3, rms_dr_vs_febio_mm=rms * 1e3, rms_over_max=rms / umax if fe is not None else None,
               flipped_base_triangles=flipped, base_triangles=int(len(tri)), gates=gates, passes=all(gates.values()),
               reported=dict(held_nodes=int(held.sum()), unilateral_nodes=int((~held).sum()),
                             held_gap_median_mm=float(np.median(np.abs(gap[held])) * 1e3) if held.any() else None,
                             held_gap_max_mm=float(np.abs(gap[held]).max() * 1e3) if held.any() else None,
                             unilateral_penetration_max_mm=float(max(0.0, -gap[~held].min()) * 1e3) if (~held).any() else None,
                             tangential_slide_median_mm=float(np.median(slide) * 1e3), tangential_slide_max_mm=float(slide.max() * 1e3),
                             base_displacement_median_mm=float(np.median(np.linalg.norm(ub, axis=1)) * 1e3),
                             base_displacement_max_mm=float(np.linalg.norm(ub, axis=1).max() * 1e3),
                             anterior_displacement_median_mm=float(np.median(ua) * 1e3), anterior_displacement_max_mm=float(ua.max() * 1e3),
                             # The registration error is absorbed here rather than corrected, so how far
                             # the tissue is moved has to be visible: a seating that moves tissue further
                             # than the breast's own dimension is not that breast any more.
                             deformation=dict(
                                 median_mm=float(np.median(np.linalg.norm(u, axis=1)) * 1e3),
                                 p90_mm=float(np.percentile(np.linalg.norm(u, axis=1), 90) * 1e3),
                                 max_mm=float(np.linalg.norm(u, axis=1).max() * 1e3),
                                 breast_cube_root_volume_mm=float((tet_volumes(X, T).sum() ** (1 / 3)) * 1e3),
                                 breast_bbox_diagonal_mm=float(np.linalg.norm(np.ptp(X, axis=0)) * 1e3),
                                 max_over_dimension=float(np.linalg.norm(u, axis=1).max() / np.linalg.norm(np.ptp(X, axis=0)))),
                             rib_points_inside_before=ribs_inside(side, X, B), rib_points_inside_after=ribs_inside(side, Y, B)),
               caveats=CAVEATS)
    e2 = d / "dr_displacement_E10000.npy"
    if e2.exists():
        u2 = np.load(e2)
        rec["reported"]["young_modulus_cancels"] = dict(E_pa=[E_PA, E_CHECK_PA], max_abs_difference_mm=float(np.abs(u2 - u).max() * 1e3),
                                                        relative_to_max_displacement=float(np.abs(u2 - u).max() / umax))
    write_obj(d / f"breast_{side}_seated.obj", Y, B, [f"{sid} {side} breast seated on this body's pectoralis major, sliding base "
                                                      f"(canonical frame, metres)"] + CAVEATS)
    (d / "judge.json").write_text(json.dumps(rec, indent=2) + "\n")
    g = rec["gates"]
    print(f"  JUDGE {sid} {side}: volume {100*(vr-1):+.3f}% {'PASS' if g['a_volume_within_1pct'] else 'FAIL'} | min J {J.min():.3f} "
          f"{'PASS' if g['b_every_J_above_0_2'] else 'FAIL'} | solvers {100*rms/umax:.2f}% of {umax*1e3:.1f} mm "
          f"{'PASS' if g['c_solvers_within_5pct'] else 'FAIL'} | flipped base triangles {flipped} "
          f"{'PASS' if g['d_no_flipped_base_triangle'] else 'FAIL'} -> {'PASS' if rec['passes'] else 'FAIL'}")
    dfm = rec['reported']['deformation']
    print(f"  deformation: median {dfm['median_mm']:.2f} mm, p90 {dfm['p90_mm']:.2f}, max {dfm['max_mm']:.2f} "
          f"against a breast {dfm['breast_bbox_diagonal_mm']:.0f} mm across ({100*dfm['max_over_dimension']:.1f}% of it)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subject", required=True, choices=sorted(REG)); ap.add_argument("--side", required=True, choices=("left", "right"))
    ap.add_argument("--no-project", action="store_true")
    # DECLARED, not chosen silently: what happens to a node whose ray leaves the bed. "hold" keeps its
    # last valid association and keeps driving it; "release" drops its constraint for that step. The
    # control is run with "hold" because releasing would change WHICH nodes are driven midway, so the
    # two arms would differ in more than staleness -- and the count is reported per step either way.
    ap.add_argument("--lost-bed", choices=("hold", "release"), default="hold")
    ap.add_argument("--stage", required=True, choices=("prepare", "place", "smooth", "control-r", "control-r-prime", "control-r2", "gate-s", "gate-t-none", "gate-t-all", "gate-u", "dr", "dr-check-E", "febio", "judge"))
    a = ap.parse_args(); d = OUT / a.subject / a.side; d.mkdir(parents=True, exist_ok=True)
    print(f"{a.subject} {a.side}: {a.stage}", flush=True)
    {"prepare": lambda: stage_prepare(a.subject, a.side, d), "place": lambda: stage_place(a.subject, a.side, d),
     "smooth": lambda: stage_smooth(a.subject, a.side, d),
     "control-r": lambda: stage_control_r(a.subject, a.side, d),
     "control-r-prime": lambda: stage_control_r(a.subject, a.side, d, bed_constraint=False),
     "control-r2": lambda: stage_control_r(a.subject, a.side, d, bed_constraint=False, freeze_frames=True),
     "gate-u": lambda: stage_control_r(a.subject, a.side, d, bed_constraint=False, freeze_frames=True,
                                       association="allornothing", lost_bed=a.lost_bed),
     "gate-t-none": lambda: stage_control_r(a.subject, a.side, d, bed_constraint=False, freeze_frames=True,
                                            association="none", lost_bed=a.lost_bed),
     "gate-t-all": lambda: stage_control_r(a.subject, a.side, d, bed_constraint=False, freeze_frames=True,
                                           association="all", lost_bed=a.lost_bed),
     "gate-s": lambda: stage_control_r(a.subject, a.side, d, bed_constraint=False, freeze_frames=True,
                                       project=not a.no_project, stop_fraction=0.25,
                                       solver_log=(lambda m: print(m, flush=True))),
     "dr": lambda: stage_dr(a.subject, a.side, d, E_PA),
     "dr-check-E": lambda: stage_dr(a.subject, a.side, d, E_CHECK_PA, "_E10000"),
     "febio": lambda: stage_febio(a.subject, a.side, d), "judge": lambda: stage_judge(a.subject, a.side, d)}[a.stage]()


# ---- the declared modelling choice: smooth the DEPTH FIELD, bounded by invertibility ------------
# (docs/BODY_PARAMETERS.md, 2a2ac9c.) 11% of base edges require their ends to slide past each other
# by more than the edge's own length, so the field the seating must apply is not a deformation. The
# field is therefore smoothed ON the base surface -- never the geometry -- with the bandwidth raised
# from zero and stopped at the FIRST value where the worst neighbour gradient reaches GRADIENT_BAR.
# Nothing here is chosen to make a number come out: the stopping point is what the solver needs to
# keep elements from inverting, and if reaching it costs more than half the field's magnitude the
# breast is being reshaped by the chest rather than seated on it, and the run stops.
GRADIENT_BAR = 0.5
SMOOTH_BANDWIDTHS_MM = (0., 1., 2., 3., 4., 5., 6., 8., 10., 12., 15., 20., 25., 30., 40., 50.)
MAGNITUDE_STOP = 0.5
ADAPTED_LABEL = ("a female breast adapted to this body's chest wall, derived from subject {sid}, "
                 "not a model of her anatomy")


def smooth_on_surface(X, nodes, values, sigma_m):
    """Gaussian averaging over the base surface: one number per node, geometry untouched."""
    if sigma_m <= 0: return values.copy()
    P = X[nodes]; tree = cKDTree(P)
    out = np.empty_like(values)
    for a, neigh in enumerate(tree.query_ball_point(P, r=3 * sigma_m, workers=-1)):
        j = np.asarray(neigh)
        w = np.exp(-((np.linalg.norm(P[j] - P[a], axis=1) / sigma_m) ** 2) / 2)
        out[a] = float((w * values[j]).sum() / w.sum())
    return out


def base_edges(X, T, base):
    isb = np.zeros(len(X), bool); isb[base] = True
    E = np.vstack([T[:, [a, b]] for a in range(4) for b in range(a + 1, 4)])
    E = np.unique(np.sort(E, 1), axis=0)
    E = E[isb[E[:, 0]] & isb[E[:, 1]]]
    order = -np.ones(len(X), np.int64); order[base] = np.arange(len(base))
    return order[E], np.linalg.norm(X[E[:, 0]] - X[E[:, 1]], axis=1)


def bandwidth_sweep(X, T, base, depth):
    """Raise the bandwidth until the worst neighbour gradient reaches the bar; report the cost."""
    E, L = base_edges(X, T, base)
    raw_mag = float(np.abs(depth).mean())
    rows = []
    chosen = None
    for sig in SMOOTH_BANDWIDTHS_MM:
        d = smooth_on_surface(X, base, depth, sig * 1e-3)
        g = np.abs(d[E[:, 0]] - d[E[:, 1]]) / np.maximum(L, 1e-9)
        removed = 1.0 - float(np.abs(d).mean()) / raw_mag
        rows.append(dict(bandwidth_mm=sig, max_gradient=float(g.max()), p90_gradient=float(np.percentile(g, 90)),
                         over_one=int((g > 1).sum()), magnitude_removed=removed,
                         median_mm=float(1000 * np.median(d[depth > 0])), p90_mm=float(1000 * np.percentile(d, 90))))
        say(f"  bandwidth {sig:5.1f} mm: worst gradient {g.max():8.2f}, p90 {np.percentile(g,90):5.2f}, "
            f"edges over 1.0: {int((g>1).sum()):5d}, magnitude removed {100*removed:5.1f}%, "
            f"depth median {1000*np.median(d[depth>0]):5.2f} mm")
        if chosen is None and g.max() <= GRADIENT_BAR:
            chosen = dict(rows[-1]); chosen["field"] = d
            break
    return chosen, rows, raw_mag


def stage_smooth(sid, side, d_dir):
    P = np.load(d_dir / "prepared.npz")
    X, T, base, gap0 = P["X"], P["T"], P["base"], P["gap0"]
    depth = np.where(gap0 < 0, -gap0, 0.0)
    say(f"{sid} {side}: raw field over {len(base)} base nodes, median {1000*np.median(depth[depth>0]):.2f} mm, "
        f"p90 {1000*np.percentile(depth,90):.2f}, max {1000*depth.max():.2f}")
    chosen, rows, raw_mag = bandwidth_sweep(X, T, base, depth)
    rec = dict(subject=sid, side=side, sweep=rows, gradient_bar=GRADIENT_BAR,
               raw=dict(median_mm=float(1000 * np.median(depth[depth > 0])),
                        p90_mm=float(1000 * np.percentile(depth, 90)), max_mm=float(1000 * depth.max()),
                        mean_magnitude_mm=1000 * raw_mag))
    if chosen is None:
        say(f"NO BANDWIDTH within {SMOOTH_BANDWIDTHS_MM[-1]:.0f} mm brings the worst gradient to {GRADIENT_BAR}")
        rec["stopped"] = "no bandwidth reached the bar"
        (d_dir / "smoothing.json").write_text(json.dumps(rec, indent=2) + "\n")
        return None
    removed = chosen["magnitude_removed"]
    say(f"BANDWIDTH REQUIRED: {chosen['bandwidth_mm']:.0f} mm; it removes {100*removed:.1f}% of the field's "
        f"magnitude; depth median {rec['raw']['median_mm']:.2f} -> {chosen['median_mm']:.2f} mm, "
        f"p90 {rec['raw']['p90_mm']:.2f} -> {chosen['p90_mm']:.2f} mm")
    rec["chosen"] = {k: v for k, v in chosen.items() if k != "field"}
    if removed > MAGNITUDE_STOP:
        say(f"STOP: reaching the bar removes {100*removed:.1f}% of the field, more than half. The breast is being "
            f"reshaped by the chest rather than seated on it; no seating is run.")
        rec["stopped"] = "more than half the field's magnitude removed"
        (d_dir / "smoothing.json").write_text(json.dumps(rec, indent=2) + "\n")
        return None
    np.save(d_dir / "smoothed_move.npy", chosen["field"])
    (d_dir / "smoothing.json").write_text(json.dumps(rec, indent=2) + "\n")
    return chosen


def stage_control_r(sid, side, d_dir, bed_constraint=True, freeze_frames=False, project=True,
                    stop_fraction=1.0, solver_log=None, association='persistent', lost_bed='hold'):
    """CONTROL R: the same 3,123 held nodes, the same bed, the same solver and the same J > 0.2
    floor, driven by a RIGID TRANSLATION of the whole base equal to the smoothed field's median
    displacement. Gate R: completes to fraction 1.0 with zero inversions."""
    P = np.load(d_dir / "prepared.npz")
    X, T, base, gap0 = P["X"], P["T"], P["base"], P["gap0"]
    held = gap0 < 0
    # R' (docs/BODY_PARAMETERS.md, 73baa16): R's premise was wrong. Only the held nodes are driven,
    # while the interior is free and the BED holds every other base node out of itself, so a body
    # whose held surface translates into a fixed bed cannot translate rigidly as a whole -- the
    # rigid motion is not admissible there and elements may legitimately deform. With the bed
    # constraint removed, the whole-body translation IS admissible and zero-energy, so nothing may
    # invert for any reason of physics or mesh quality. That is the control R should have been.
    if not bed_constraint:
        base = base[held]; gap0 = gap0[held]; held = gap0 < 0
    mu, lam = lame(E_PA, NU)
    region = SlidingRegion(X, T, mu_pa=mu, lambda_pa=lam, density_kg_m3=950.0)
    dirs = P["ray_directions"] if bed_constraint else P["ray_directions"][np.load(d_dir / "prepared.npz")["gap0"] < 0]
    closest, _ = bed_rays(P["bedV"], P["bedF"], dirs, RAY_REACH_M)
    c0, n0 = closest(X[base])
    smoothed = d_dir / "smoothed_move.npy"
    travel = np.load(smoothed) if smoothed.exists() else np.where(gap0 < 0, -gap0, 0.0)
    if not bed_constraint and len(travel) != len(base):
        travel = travel[np.load(d_dir / "prepared.npz")["gap0"] < 0]
    magnitude = float(np.median(travel[held]))
    direction = n0[held].mean(0); direction /= np.linalg.norm(direction)
    rigid = magnitude * direction
    say(f"{sid} {side} CONTROL {'R' if bed_constraint else 'R-prime (bed constraint removed)'}: rigid translation of {1000*magnitude:.2f} mm along "
        f"[{direction[0]:+.3f}, {direction[1]:+.3f}, {direction[2]:+.3f}] "
        f"({int(held.sum())} held nodes, same bed, same solver, same J > 0.2 floor)")
    t0 = time.time()
    try:
        # The association's jump limit must exceed how far a node moves in one increment, or the
        # update is rejected, the constraint normal goes stale, and the target stops being consistent
        # with the motion being driven. At 5e-4 m the log pinned at 0.4994 mm every pass while nodes
        # moved 0.94 mm.
        r = seat_on_bed(region, base, closest, load_steps=LOAD_STEPS, gap_tol_m=GAP_TOL_M,
                        jump_limit_m=CONTROL_JUMP_LIMIT_M, assoc_tol_m=ASSOC_TOL_M, rigid_m=rigid,
                        freeze_frames=freeze_frames, project=project, stop_fraction=stop_fraction,
                        solver_log=solver_log, association=association, lost_bed=lost_bed,
                        facet_m=GAP_TOL_M,
                        log=lambda m, flush=True: print(m, flush=True))
    except Exception as failure:
        say(f"GATE {'R' if bed_constraint else 'R-prime'}: FAILED -- {type(failure).__name__}: {failure}")
        say("  Reported as observed. R's premise was withdrawn (73baa16): only the held nodes are driven "
            "while the bed holds the others out of itself, so a whole-body rigid translation is NOT "
            "admissible here and elements may legitimately deform. R-prime, with the bed removed, is the "
            "control that can carry that claim.")
        (d_dir / ("control_r.json" if bed_constraint else "control_r_prime.json")).write_text(json.dumps(dict(passes=False, error=str(failure),
            rigid_mm=1000 * magnitude, direction=direction.tolist(), seconds=time.time() - t0), indent=2) + "\n")
        return False
    steps = r["steps"]
    if r.get("refusal_rate") is not None:
        say(f"  REFUSAL RATE: {100*r['refusal_rate']:.1f}% of association updates refused "
            f"({sum(r['refusals'])} of {len(r['refusals'])}); threshold = the step's node motion + the "
            f"{1000*GAP_TOL_M:.0f} mm facet scale")
    if r.get("lost_bed_per_association"):
        say(f"  nodes whose ray left the bed, per association: {r['lost_bed_per_association']} "
            f"(rule: {r['lost_bed_rule']}, association: {r['association']})")
    if freeze_frames:
        say("  min J across the run, with the association frozen inside each step:")
        say(f"    {'step':>5} {'fraction':>9} {'min J at entry':>15} {'min J at exit':>14} {'drop in step':>13} "
            f"{'association moved before it':>28}")
        prev_exit = None
        for k, st in enumerate(steps, 1):
            drop_in = st['min_J_entry'] - st['min_J']
            across = "" if prev_exit is None else f"{prev_exit - st['min_J_entry']:+.4f}"
            say(f"    {k:5d} {st['fraction']:9.4f} {st['min_J_entry']:15.4f} {st['min_J']:14.4f} {drop_in:+13.4f} "
                f"{1000*st['association_moved_m']:9.3f} mm, J across it {across:>9}")
            prev_exit = st['min_J']
        within = sum(st['min_J_entry'] - st['min_J'] for st in steps)
        across_total = sum((steps[i - 1]['min_J'] - steps[i]['min_J_entry']) for i in range(1, len(steps)))
        say(f"  total fall WITHIN steps (constraints fixed): {within:+.4f}; ACROSS re-associations: {across_total:+.4f}")
        if abs(within) < 1e-6 and abs(across_total) < 1e-6:
            say("  -> no min J is lost anywhere: every step holds the exact solution")
        elif across_total > within:
            say("  -> min J is lost across the re-associations")
        else:
            say("  -> min J is lost inside steps; note that a re-association moves no nodes, so its effect "
                "appears in the NEXT step's start, which is inside a step (gate S)")
    u = r["displacement"]
    Y = X + u
    J = np.linalg.det(np.swapaxes(Y[T[:, 1:]] - Y[T[:, 0, None]], 1, 2)
                      @ np.linalg.inv(np.swapaxes(X[T[:, 1:]] - X[T[:, 0, None]], 1, 2)))
    drift = np.linalg.norm(u - rigid, axis=1)
    ok = bool(J.min() > 0.2)
    say(f"GATE {'R' if bed_constraint else 'R-prime'}: completed to fraction 1.0 in {time.time()-t0:.0f} s, min J {J.min():.4f}, "
        f"inversions {int((J <= 0).sum())} -> {'PASS' if ok else 'FAIL'}")
    say(f"  the solution against the ideal rigid translation: median |u - d| {1000*np.median(drift):.3f} mm, "
        f"max {1000*drift.max():.3f} mm; volume ratio {tet_volumes(Y, T).sum()/tet_volumes(X, T).sum():.6f}")
    (d_dir / ("control_r.json" if bed_constraint else ("control_r_double_prime.json" if freeze_frames else "control_r_prime.json"))).write_text(json.dumps(dict(passes=ok, steps=steps, rigid_mm=1000 * magnitude,
        direction=direction.tolist(), min_J=float(J.min()), inversions=int((J <= 0).sum()),
        drift_median_mm=float(1000 * np.median(drift)), drift_max_mm=float(1000 * drift.max()),
        volume_ratio=float(tet_volumes(Y, T).sum() / tet_volumes(X, T).sum()),
        cutbacks=r["cutbacks"], seconds=time.time() - t0), indent=2) + "\n")
    return ok


if __name__ == "__main__": main()
