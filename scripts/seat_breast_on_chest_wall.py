"""seat each registered female breast on this body's chest wall by a soft-tissue deformation.

WHY. Four women registered into this body (s0790, s1067, s1159, s0970) all show this body's
anterior chest wall in front of theirs: its ribs 2-7 sit inside their mapped breasts by up to
~2 cm. No correction that moves the WHOLE breast can seat it: clearing the rib-5 overlap lifts the
rest of the base off the chest wall (scripts/fit_breast_seating_correction.py). This deforms the
breast instead: its BASE is placed on this body's bed and the tissue is solved as a
near-incompressible hyperelastic solid.

THE BED. The anterior (forward-facing) surface of this body's pectoralis major -- clavicular,
sternocostal and abdominal parts, BodyParts3D only -- the breast's anatomical bed, not the ribs.

THE BASE RULE, fixed before any breast was judged. Candidates are the tet mesh's
posterior-facing boundary nodes (outward normal z < 0; +z is anterior in this frame) within
BASE_REACH_M of the bed. A candidate that PENETRATES the bed -- lies behind the outward normal of
its closest bed face -- is placed on its closest bed point, at ANY depth. A candidate with a GAP (in
front of the bed) is left FREE, and only those are subject to BASE_REACH_M, which no longer affects
the result since gapped nodes are free either way.
(Why no depth cap. A 30 mm cap was first justified from the rib penetrations, ~21 mm at worst. But
the bed is the pectoralis, which sits a measured median 20.2 mm in front of this body's ribs 3-6, so
the breast penetrates the BED by up to ~42.6 mm on s1159-left -- 281 penetrating nodes lay deeper
than 30 mm and would have been left free INSIDE the muscle. Their closest bed points sit well within
the muscle's height, so they are not edge artefacts. A node behind the bed can only be seated by
being pushed out; capping that depth was the error. Changed before any solve finished.) (A first version placed EVERY candidate on the bed, pulling gapped nodes in: on s1159-left
that moved the base a median 12.8 mm and up to the full 30 mm reach, far beyond the 5-8 mm the rib
penetrations call for. It was stopped before any solve finished and replaced by this rule, which
is the one the task specified.)

THE SOLVE. PrescribedRegion (ihm/assembly/prescribed_deformation.py, a subclass of
DeformableRegion) at nu = 0.49, cross-checked by FEBio 4.13 on the same mesh and boundary
conditions (ihm/assembly/febio_crosscheck.py). With displacement-only boundary conditions and no
body force, Young's modulus cancels from the resting shape; that is shown, not assumed, by solving
one breast at two values of E.

THE JUDGE (docs/BODY_PARAMETERS.md, committed before any code):
  solver first  nu=0.49 block volume within 1%; DeformableRegion vs FEBio within 5% of max
                displacement, RMS over nodes, on the same mesh and BCs
  per breast    (a) volume within 1%; (b) every tet J > 0.2; (c) the two solvers within 5% of
                the deformed breast's maximum displacement
  reported, NOT judged (the BCs enforce them): base contact median <= 3 mm; this body's rib points
                inside the breast <= 1%

CAVEATS: one clinical subject per breast, as a segmentation model drew it; 'breast' is one
soft-tissue label with no gland, ducts or nipple; nu = 0.49 is assumed; gravity and the unloaded
supine reference shape are out of scope; the mapped female TRUNK is not changed, so this body's
sternum still sits outside it (the separate skin-envelope problem).

STATUS (2026-09-10): the solver-first gate PASSES and the per-breast gate STOPS on s1159-left.
The block known-answer case passes -- nu=0.49 compression -0.221% volume; DeformableRegion and FEBio
agree to 1e-5 / 9e-5 % of max displacement, and a mismatched-case control reads 124.8%. But the base
rule below, closest-point projection onto the bed, is ILL-POSED: on s1159-left it collapses 69.4% of
the base triangles below 0.2 of their area (20.4% flip), so a tet on such a face cannot reach J > 0.2
whatever its fourth node does. FEBio fails at t = 0.227 with 53 negative-Jacobian reports. No seated
breast is written. The fix is a different boundary condition -- constrain the base only along the
bed's normal and let it slide tangentially, cross-checked by FEBio sliding-elastic contact with
tension -- which the committed judge does not yet cover.

Stages save their arrays under OUT/<subject>/<side>/ so each can run, be killed and resume:
  prepare  tet-mesh the breast, build the bed, select and project the base
  dr       solve with PrescribedRegion (optionally at a second E)
  febio    solve the same problem in FEBio
  judge    apply the judge and write the deformed surface and the manifest entry
"""
import argparse, gzip, importlib.util, json, subprocess, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.prescribed_deformation import PrescribedRegion, lame, tet_volumes   # noqa: E402
from ihm.assembly.febio_crosscheck import write_feb, run_febio, read_msh_tets        # noqa: E402

OUT = ROOT / "data/derived/female-breast-seated-v1"
REG = {"s0790": "female-torso-registered-v1", "s1067": "female-torso-registered-s1067-v1",
       "s1159": "female-torso-registered-s1159-v1", "s0970": "female-torso-registered-s0970-v1"}
PECTORALIS = {"right": ["body-bp3d-FJ1447", "body-bp3d-FJ1464", "body-bp3d-FJ1446"],
              "left": ["body-bp3d-FJ1447M", "body-bp3d-FJ1464M", "body-bp3d-FJ1446M"]}
FTETWILD = ROOT / "data/runtime/tolerant-mesher/build/FloatTetwild_bin"
BASE_REACH_M, NU, E_PA, E_CHECK_PA = 0.030, 0.49, 1000.0, 10000.0
# fTetWild's relative target edge length. 0.05 gave ~19k nodes per breast, for which FEBio's skyline
# solver stored 870 million entries and had not finished one factorisation in 12 minutes. 0.08 is
# ~4x fewer nodes; mesh resolution is not a judged gate, and both solvers share the discretisation.
MESH_LR = 0.08
# ...and -l alone does not coarsen: with a tight envelope (-e 1e-3) fTetWild must resolve every input
# triangle, and the breast surfaces are fine marching-cubes meshes (~67k faces). -l 0.08 gave 19,009
# nodes, the same as -l 0.05's 19,151. So the surface is DECIMATED first, to DECIMATE_FACES, gated:
# the decimated surface must keep the original breast volume within DECIMATE_VOLUME_TOL.
# Chosen ONCE for all eight breasts, before any solve, from a survey of every breast at 8k/12k/16k/24k
# faces: igl.decimate (shortest-edge collapse) fails the gate at 8k (worst -0.623%, s1159 right);
# igl.qslim (quadric error) passes at 8k on every breast, worst -0.105% -- five times inside the gate.
DECIMATE_FACES, DECIMATE_VOLUME_TOL = 8000, 0.005
ANTERIOR_NZ = 0.2          # a bed face is anterior-facing if its outward normal has z > this
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
    """outward-oriented boundary triangles of positively oriented tets."""
    f = np.vstack([T[:, [1, 2, 3]], T[:, [0, 3, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 1]]])
    key = np.sort(f, 1); _, idx, cnt = np.unique(key, axis=0, return_index=True, return_counts=True)
    return f[idx[cnt == 1]]

def normals(V, F):
    n = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    return n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-30)

def bed(side):
    ents = {e["id"]: e for e in json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]}
    Vs, Fs, off = [], [], 0
    for i in PECTORALIS[side]:
        g = json.loads(gzip.decompress((ROOT / ents[i]["reference_geometry"]["path"]).read_bytes()))
        V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
        t = V[F]
        if np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum() < 0: F = F[:, [0, 2, 1]]   # make outward
        Vs.append(V); Fs.append(F + off); off += len(V)
    V = np.vstack(Vs); F = np.vstack(Fs)
    return V, F[normals(V, F)[:, 2] > ANTERIOR_NZ]

def stage_prepare(sid, side, d):
    import igl
    src = ROOT / "data/derived" / REG[sid] / f"breast_{side}.obj"
    md = d / f"mesh_dec{DECIMATE_FACES}"; md.mkdir(parents=True, exist_ok=True)
    msh = md / "out.msh"
    Vo, Fo = read_obj(src); to = Vo[Fo]; vol_o = abs(np.einsum("ij,ij->i", to[:, 0], np.cross(to[:, 1], to[:, 2])).sum() / 6)
    dec = md / "decimated.obj"
    if not dec.exists():
        # this igl build returns (V, F, J, I) with no success flag; pick the arrays by type rather than position
        out = [np.asarray(o) for o in igl.qslim(np.ascontiguousarray(Vo), np.ascontiguousarray(Fo), DECIMATE_FACES)
               if hasattr(o, "shape")]
        U = next(o for o in out if o.dtype.kind == "f" and o.ndim == 2 and o.shape[1] == 3)
        G = next(o for o in out if o.dtype.kind in "iu" and o.ndim == 2 and o.shape[1] == 3)
        print(f"  qslim returned {[(o.dtype.kind, o.shape) for o in out]}; vertices {U.shape}, faces {G.shape}")
        if G.max() >= len(U): raise SystemExit("decimated faces index past the decimated vertices")
        tu = U[G]; vol_d = abs(np.einsum("ij,ij->i", tu[:, 0], np.cross(tu[:, 1], tu[:, 2])).sum() / 6)
        change = vol_d / vol_o - 1
        print(f"  decimated {len(Fo)} -> {len(G)} faces; volume {vol_o*1e6:.1f} -> {vol_d*1e6:.1f} mL ({100*change:+.3f}%)")
        if abs(change) > DECIMATE_VOLUME_TOL: raise SystemExit(f"decimation changed the breast volume by {100*change:+.3f}% (> {100*DECIMATE_VOLUME_TOL}%)")
        write_obj(dec, U, G, [f"{sid} {side} breast, qslim-decimated to {len(G)} faces for meshing ({100*change:+.3f}% volume)"])
    if not msh.exists():
        cmd = [str(FTETWILD), "-i", str(dec), "-o", str(msh), "--no-binary", "-e", "1e-3", "-l", str(MESH_LR), "--max-threads", "12"]
        t0 = time.time(); p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        (md / "ftetwild.log").write_text(p.stdout[-4000:] + p.stderr[-4000:])
        if p.returncode != 0 or not msh.exists(): raise SystemExit(f"fTetWild failed on {sid} {side}")
        print(f"  fTetWild {time.time()-t0:.0f} s")
    X, T = read_msh_tets(msh)
    used = np.unique(T); remap = -np.ones(len(X), int); remap[used] = np.arange(len(used)); X = X[used]; T = remap[T]
    vol = tet_volumes(X, T); T[vol < 0] = T[vol < 0][:, [0, 2, 1, 3]]
    Vs, Fs = read_obj(src); t = Vs[Fs]; surf_ml = abs(np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum() / 6) * 1e6
    mesh_ml = tet_volumes(X, T).sum() * 1e6
    B = boundary_faces(T); nb = normals(X, B)
    post = np.unique(B[nb[:, 2] < 0]); ant = np.unique(B[nb[:, 2] > 0])
    bV, bF = bed(side)
    d2, I, C = igl.point_mesh_squared_distance(X[post], bV, bF)
    behind = np.einsum("ij,ij->i", X[post] - C, normals(bV, bF)[I]) < 0          # penetrating the bed
    take = behind                                   # every penetrating node, at any depth
    base = post[take]; target = C[take]
    gapped = int(((np.sqrt(d2) <= BASE_REACH_M) & ~behind).sum())
    mask = np.zeros_like(X, bool); disp = np.zeros_like(X)
    mask[base] = True; disp[base] = target - X[base]
    np.savez(d / "prepared.npz", X=X, T=T, mask=mask, disp=disp, base=base, post=post, ant=ant, boundary=B)
    m = np.linalg.norm(disp[base], axis=1)
    info = dict(nodes=len(X), tets=len(T), mesh_ml=mesh_ml, surface_ml=surf_ml, posterior_nodes=len(post),
                base_nodes=len(base), gapped_candidates_left_free=gapped, base_move_median_mm=float(np.median(m) * 1e3), base_move_max_mm=float(m.max() * 1e3))
    (d / "prepared.json").write_text(json.dumps(info, indent=2) + "\n")
    print(f"  {len(X)} nodes, {len(T)} tets; tet volume {mesh_ml:.1f} mL vs surface {surf_ml:.1f} mL ({100*(mesh_ml/surf_ml-1):+.2f}%)")
    print(f"  base: {len(base)} penetrating nodes projected, {gapped} gapped candidates left free (of {len(post)} posterior boundary nodes); "
          f"moved median {info['base_move_median_mm']:.1f} mm, max {info['base_move_max_mm']:.1f} mm")

def stage_dr(sid, side, d, young):
    P = np.load(d / "prepared.npz")
    mu, lam = lame(young, NU)
    body = PrescribedRegion(P["X"], P["T"], mu_pa=mu, lambda_pa=lam, density_kg_m3=950.0)
    r = body.solve_prescribed(P["mask"], P["disp"])
    tag = "" if young == E_PA else f"_E{int(young)}"
    np.save(d / f"dr_displacement{tag}.npy", r["displacement"])
    keep = {k: v for k, v in r.items() if k not in ("positions", "displacement")}
    (d / f"dr{tag}.json").write_text(json.dumps(keep, indent=2) + "\n")
    print(f"  DR E={young:g}: volume ratio {r['volume_ratio']:.5f}, min J {r['minimum_jacobian']:.4f}, residual/reaction "
          f"{r['residual_relative']:.2e}, {r['iterations']} iterations, {r['wall_seconds']:.0f} s")

def stage_febio(sid, side, d):
    P = np.load(d / "prepared.npz"); fd = d / "febio"; fd.mkdir(exist_ok=True)
    write_feb(fd / "seat.feb", P["X"], P["T"], P["mask"], P["disp"], young_pa=E_PA, nu=NU, steps=20, auto_stepper=True, reorder=True)
    t0 = time.time(); r = run_febio(fd / "seat.feb", threads=12)
    np.save(d / "febio_displacement.npy", r["displacement"])
    (d / "febio.json").write_text(json.dumps({k: v for k, v in r.items() if k != "displacement"}, indent=2) + "\n")
    print(f"  FEBio: final time {r['final_time']:.3f}, normal termination {r['normal_termination']}, {time.time()-t0:.0f} s")

def ribs_inside(side, V, F):
    spec = importlib.util.spec_from_file_location("cwo", ROOT / "scripts/measure_female_chest_wall_offset.py")
    M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
    ents = {e["name"]: e for e in json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
            if e["role"] == "rigid_bone" and e["id"].startswith("body-bp3d-")}
    pts = []
    for o in ("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth"):
        g = json.loads(gzip.decompress((ROOT / ents[f"{side} {o} rib"]["reference_geometry"]["path"]).read_bytes()))
        P = np.asarray(g["positions"], float).reshape(-1, 3)
        pts.append(P[np.random.default_rng(0).choice(len(P), min(len(P), 1500), replace=False)])
    return float(M.inside(V, F, np.vstack(pts)).mean())

def stage_judge(sid, side, d):
    import igl
    P = np.load(d / "prepared.npz"); X, T, mask, B = P["X"], P["T"], P["mask"], P["boundary"]
    u = np.load(d / "dr_displacement.npy"); f = np.load(d / "febio_displacement.npy")
    Y = X + u; J = np.linalg.det(np.swapaxes(Y[T[:, 1:]] - Y[T[:, 0, None]], 1, 2) @
                                 np.linalg.inv(np.swapaxes(X[T[:, 1:]] - X[T[:, 0, None]], 1, 2)))
    vr = tet_volumes(Y, T).sum() / tet_volumes(X, T).sum()
    umax = float(np.linalg.norm(u, axis=1).max()); rms = float(np.sqrt(np.mean(np.sum((u - f) ** 2, 1))))
    g = dict(volume_within_1pct=abs(vr - 1) <= 0.01, every_J_above_0_2=bool(J.min() > 0.2), solvers_within_5pct=rms / umax <= 0.05)
    bV, bF = bed(side)
    post = P["post"]; dcont = np.sqrt(igl.point_mesh_squared_distance(Y[post], bV, bF)[0])
    ant = P["ant"]; ua = np.linalg.norm(u[ant], axis=1); ub = np.linalg.norm(u[P["base"]], axis=1)
    rib_before, rib_after = ribs_inside(side, X, B), ribs_inside(side, Y, B)
    rec = dict(volume_ratio=float(vr), min_J=float(J.min()), rms_dr_vs_febio_mm=rms * 1e3, max_displacement_mm=umax * 1e3,
               rms_over_max=rms / umax, gates=g, passes=all(g.values()),
               reported=dict(base_contact_median_mm=float(np.median(dcont) * 1e3), posterior_nodes_to_bed_p95_mm=float(np.percentile(dcont, 95) * 1e3),
                             rib_points_inside_before=rib_before, rib_points_inside_after=rib_after,
                             base_displacement_median_mm=float(np.median(ub) * 1e3), base_displacement_max_mm=float(ub.max() * 1e3),
                             anterior_displacement_median_mm=float(np.median(ua) * 1e3), anterior_displacement_max_mm=float(ua.max() * 1e3)))
    E2 = d / "dr_displacement_E10000.npy"
    if E2.exists():
        u2 = np.load(E2); rec["young_modulus_cancels"] = dict(E_pa=[E_PA, E_CHECK_PA], max_abs_difference_mm=float(np.abs(u2 - u).max() * 1e3),
                                                               relative=float(np.abs(u2 - u).max() / umax))
    write_obj(d / f"breast_{side}_seated.obj", Y, B[:, [0, 2, 1]] if False else B,
              [f"{sid} {side} breast seated on this body's pectoralis major (canonical frame, metres)"] + CAVEATS)
    (d / "judge.json").write_text(json.dumps(rec, indent=2) + "\n")
    print(f"  JUDGE {sid} {side}: volume {100*(vr-1):+.3f}% | min J {J.min():.3f} | DR vs FEBio RMS {100*rms/umax:.3f}% of {umax*1e3:.1f} mm "
          f"-> {'PASS' if rec['passes'] else 'FAIL'} | contact median {rec['reported']['base_contact_median_mm']:.2f} mm | "
          f"ribs inside {100*rib_before:.2f}% -> {100*rib_after:.2f}%")

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subject", required=True, choices=sorted(REG)); ap.add_argument("--side", required=True, choices=("left", "right"))
    ap.add_argument("--stage", required=True, choices=("prepare", "dr", "dr-check-E", "febio", "judge"))
    a = ap.parse_args(); d = OUT / a.subject / a.side; d.mkdir(parents=True, exist_ok=True)
    print(f"{a.subject} {a.side}: {a.stage}")
    {"prepare": lambda: stage_prepare(a.subject, a.side, d), "dr": lambda: stage_dr(a.subject, a.side, d, E_PA),
     "dr-check-E": lambda: stage_dr(a.subject, a.side, d, E_CHECK_PA), "febio": lambda: stage_febio(a.subject, a.side, d),
     "judge": lambda: stage_judge(a.subject, a.side, d)}[a.stage]()

if __name__ == "__main__": main()
