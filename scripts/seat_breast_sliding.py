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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.sliding_contact import SlidingRegion, seat_on_bed                       # noqa: E402
from ihm.assembly.prescribed_deformation import lame, tet_volumes                          # noqa: E402
from ihm.assembly.febio_crosscheck import read_msh_tets, run_febio, rcm_order              # noqa: E402
from ihm.assembly.febio_sliding import write_sliding_feb, split_base_faces                 # noqa: E402

OUT = ROOT / "data/derived/female-breast-sliding-v1"
REG = {"s0790": "female-torso-registered-v1", "s1067": "female-torso-registered-s1067-v1",
       "s1159": "female-torso-registered-s1159-v1", "s0970": "female-torso-registered-s0970-v1"}
PECTORALIS = {"right": ["body-bp3d-FJ1447", "body-bp3d-FJ1464", "body-bp3d-FJ1446"],
              "left": ["body-bp3d-FJ1447M", "body-bp3d-FJ1464M", "body-bp3d-FJ1446M"]}
FTETWILD = ROOT / "data/runtime/tolerant-mesher/build/FloatTetwild_bin"
NU, E_PA, E_CHECK_PA = 0.49, 1000.0, 10000.0
DECIMATE_FACES, DECIMATE_VOLUME_TOL, MESH_LR = 8000, 0.005, 0.08
# A ray finds the bed only NEARBY: uncapped, a ray grazes off to a distant fold of the muscle and
# calls it the bed (175 mm of 'penetration' on a breast 13 cm deep). 60 mm is well beyond the
# deepest real penetration measured on these four subjects (43 mm).
# A ray that meets the bed at a shallow angle has no well-defined bed direction: a micrometre of
# motion slides its hit point millimetres along the surface, which reads as the constraint moving.
# Those nodes are free, like the ones with no muscle on their line.
RAY_REACH_M, GRAZING_COS, LOAD_STEPS, GAP_TOL_M = 0.060, 0.3, 8, 5e-5
CAVEATS = ["one clinical subject per breast, as a segmentation model drew it",
           "'breast' is a single soft-tissue label: no gland, ducts or nipple",
           "nu = 0.49 is assumed (adipose is nearly incompressible)",
           "gravity and the unloaded supine reference shape are out of scope",
           "the mapped female TRUNK is not changed; this body's sternum still sits outside it (skin-envelope problem)"]


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


def bed(side):
    """The anterior-facing surface of this body's pectoralis major, and its rim vertices."""
    ents = {e["id"]: e for e in json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]}
    Vs, Fs, off = [], [], 0
    for i in PECTORALIS[side]:
        g = json.loads(gzip.decompress((ROOT / ents[i]["reference_geometry"]["path"]).read_bytes()))
        V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
        t = V[F]
        if np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum() < 0: F = F[:, [0, 2, 1]]
        Vs.append(V); Fs.append(F + off); off += len(V)
    V = np.vstack(Vs); F = np.vstack(Fs)[face_normals(np.vstack(Vs), np.vstack(Fs))[:, 2] > 0]
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
    def association(P):
        P = np.ascontiguousarray(np.asarray(P, float))
        t_back, f_back = ray_hits(P, D, V, F)          # along the outward normal: bed behind
        t_front, f_front = ray_hits(P, -D, V, F)       # against it: bed in front
        t_back[t_back > reach_m] = np.inf; t_front[t_front > reach_m] = np.inf
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
        return (ray_hits(P, D, V, F)[0] <= reach_m) | (ray_hits(P, -D, V, F)[0] <= reach_m)
    return association, has_bed


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
    B = boundary_faces(T); nb = face_normals(X, B)
    posterior = np.unique(B[nb[:, 2] < 0]); anterior = np.unique(B[nb[:, 2] > 0])
    bV, bF, rim = bed(side)
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
                deepest_penetration_mm=float(-gap0.min() * 1e3))
    (d / "prepared.json").write_text(json.dumps(info, indent=2) + "\n")
    print(f"  {len(X)} nodes, {len(T)} tets; base {len(base)} of {len(posterior)} posterior nodes "
          f"({int((~on_line).sum())} free: no muscle along their ray): {info['held']} held, {info['unilateral']} unilateral; "
          f"deepest penetration {info['deepest_penetration_mm']:.1f} mm")


def stage_dr(sid, side, d, young, tag=""):
    P = np.load(d / "prepared.npz")
    mu, lam = lame(young, NU)
    region = SlidingRegion(P["X"], P["T"], mu_pa=mu, lambda_pa=lam, density_kg_m3=950.0)
    closest, _ = bed_rays(P["bedV"], P["bedF"], P["ray_directions"], RAY_REACH_M)
    t0 = time.time()
    r = seat_on_bed(region, P["base"], closest, load_steps=LOAD_STEPS, gap_tol_m=GAP_TOL_M, jump_limit_m=5e-4,
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subject", required=True, choices=sorted(REG)); ap.add_argument("--side", required=True, choices=("left", "right"))
    ap.add_argument("--stage", required=True, choices=("prepare", "dr", "dr-check-E", "febio", "judge"))
    a = ap.parse_args(); d = OUT / a.subject / a.side; d.mkdir(parents=True, exist_ok=True)
    print(f"{a.subject} {a.side}: {a.stage}", flush=True)
    {"prepare": lambda: stage_prepare(a.subject, a.side, d), "dr": lambda: stage_dr(a.subject, a.side, d, E_PA),
     "dr-check-E": lambda: stage_dr(a.subject, a.side, d, E_CHECK_PA, "_E10000"),
     "febio": lambda: stage_febio(a.subject, a.side, d), "judge": lambda: stage_judge(a.subject, a.side, d)}[a.stage]()


if __name__ == "__main__": main()
