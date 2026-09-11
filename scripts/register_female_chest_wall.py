"""Register each subject's CHEST WALL onto this body, and judge the registration.

WHY. Four seating instruments were built and judged against the whole-torso similarity's placement
of these breasts, and all four failed on the same thing: 2.0-7.7% of every breast sits more than
20 mm behind this body's muscular chest wall, by an amount that scales with each subject's own
chest-wall offset (docs/BODY_PARAMETERS.md). That is not a label error and not a solver error; one
similarity fitted to a whole torso cannot place a breast against this body's chest wall. Her CT
carries ribs 1-12 per side, sternum, clavicles and T1-T12 as separate labels, so a chest-local fit
needs none of her soft tissue.

THE INSTRUMENT. One similarity per subject fitted to her ribs 2-7, both sides, onto this body's, by
the one-way trimmed ICP the pelvic registration uses (scripts/register_pelvic_organs.py): her
PARTIAL scan surfaces are sampled and matched onto this body's COMPLETE bones, never the reverse,
so the missing parts of a cut rib pull nothing. Ribs 2-7 because that is the breast's measured base.

GATES (docs/BODY_PARAMETERS.md, "The chest-wall registration", fixed before this was built):
  1 known answer  this body's own chest bones, truncated the way her scan cuts them -- to the part
                  her OWN labels cover, per subject, which is a median 0.80-0.89 of each rib's
                  surface -- moved by a known similarity, recovered within 2 mm, 2 deg, 1%
  2 laterality    her left ribs map nearer this body's left ribs than its right
  3 held out      ribs 2-7 are FITTED; sternum, clavicles, rib 1, ribs 8-12 and T1-T12 are held
                  out, and the median nearest-surface distance there is no worse than the
                  whole-torso similarity achieves on the same structures, subject by subject
  4 consequence   breast volume more than 20 mm behind the muscular chest wall <= 1% per subject
Reported, not judged: how far each chest-local transform differs from its whole-torso one, and the
chest-wall offset recomputed under it.
"""
import argparse, gzip, importlib.util, json, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location("seat", ROOT / "scripts/seat_breast_sliding.py")
SEAT = importlib.util.module_from_spec(spec); spec.loader.exec_module(SEAT)   # bed, rays, trim machinery

SUBJECTS = ("s0790", "s1067", "s1159", "s0970")
SRC = ROOT / "data/derived/female-torso-totalsegmentator-v1"
REGISTERED = {"s0790": "female-torso-registered-v1", "s1067": "female-torso-registered-s1067-v1",
              "s1159": "female-torso-registered-s1159-v1", "s0970": "female-torso-registered-s0970-v1"}
OUT = ROOT / "data/derived/female-chest-wall-registered-v1"
ORD = ("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth")
STERNUM_PARTS = ["manubrium", "body of sternum", "xiphoid process"]
FIT_LABELS = [f"rib_{s}_{i}" for s in ("left", "right") for i in range(2, 8)]
HELD_OUT = (["sternum", "clavicula_left", "clavicula_right"] +
            [f"rib_{s}_{i}" for s in ("left", "right") for i in [1] + list(range(8, 13))] +
            [f"vertebrae_T{i}" for i in range(1, 13)])
KA_T_MM, KA_R_DEG, KA_S = 2.0, 2.0, 0.01
TRIM_DEPTH_M, TRIM_VOLUME_TOL = 0.020, 0.01
N_SRC, N_TGT, TRIM, ITERS, TOL = 4000, 40000, 0.2, 60, 1e-9


def say(*a): print(*a, flush=True)
def sim(s, R, t): M = np.eye(4); M[:3, :3] = s * R; M[:3, 3] = t; return M
def apply(M, p): return np.asarray(p) @ M[:3, :3].T + M[:3, 3]
def scale_of(M): return float(np.cbrt(np.linalg.det(M[:3, :3])))


def umeyama(a, b):
    ca, cb = a.mean(0), b.mean(0); x, y = a - ca, b - cb
    u, s, vt = np.linalg.svd(x.T @ y / len(a))
    d = np.diag([1, 1, np.sign(np.linalg.det(vt.T @ u.T))])
    r = vt.T @ d @ u.T; scale = float((s @ np.diag(d)).sum() / (x ** 2).sum() * len(a))
    return scale, r, cb - scale * r @ ca


def areas(V, F): t = V[F]; return np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2


def area_samples(V, F, n, seed=0):
    rng = np.random.default_rng(seed); a = areas(V, F)
    f = rng.choice(len(F), n, p=a / a.sum()); u = rng.random((n, 1)); v = rng.random((n, 1))
    over = (u + v > 1); u[over] = 1 - u[over]; v[over] = 1 - v[over]
    t = V[F[f]]
    return t[:, 0] + u * (t[:, 1] - t[:, 0]) + v * (t[:, 2] - t[:, 0])


def body_meshes():
    """this body's chest bones, keyed by her label names."""
    ents = {e["name"]: e for e in json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
            if e["role"] == "rigid_bone" and e["id"].startswith("body-bp3d-")}
    def mesh(name):
        g = json.loads(gzip.decompress((ROOT / ents[name]["reference_geometry"]["path"]).read_bytes()))
        return np.asarray(g["positions"], float).reshape(-1, 3), np.asarray(g["indices"], np.int64).reshape(-1, 3)
    out = {}
    for side in ("left", "right"):
        for i, o in enumerate(ORD, start=1): out[f"rib_{side}_{i}"] = mesh(f"{side} {o} rib")
        out[f"clavicula_{side}"] = mesh(f"{side} clavicle")
    for i, o in enumerate(ORD, start=1): out[f"vertebrae_T{i}"] = mesh(f"{o} thoracic vertebra")
    Vs, Fs, off = [], [], 0
    for part in STERNUM_PARTS:
        V, F = mesh(part); Vs.append(V); Fs.append(F + off); off += len(V)
    out["sternum"] = (np.vstack(Vs), np.vstack(Fs))
    return out


def subject_meshes(subject):
    out = {}
    for p in sorted((SRC / subject / "meshes").glob("*.obj")):
        if p.stem in ("skin", "breast_left", "breast_right"): continue
        V, F = SEAT.read_obj(p)
        if len(F): out[p.stem] = (V, F)
    return out


def fit(src, body, labels, seed=0):
    """one-way trimmed ICP: her PARTIAL surfaces sampled onto this body's COMPLETE bones."""
    labs = [l for l in labels if l in src and l in body]
    ca = np.array([area_samples(*src[l], 3000).mean(0) for l in labs])
    cb = np.array([area_samples(*body[l], 3000).mean(0) for l in labs])
    M = sim(*umeyama(ca, cb))
    S = {l: area_samples(*src[l], N_SRC, seed + i) for i, l in enumerate(labs)}
    T = {l: cKDTree(area_samples(*body[l], N_TGT, 100 + i)) for i, l in enumerate(labs)}
    it = 0
    for it in range(1, ITERS + 1):
        A_, B_, D_ = [], [], []
        for l in labs:
            P = apply(M, S[l]); d, j = T[l].query(P)
            A_.append(S[l]); B_.append(T[l].data[j]); D_.append(d)
        A_, B_, D_ = np.vstack(A_), np.vstack(B_), np.concatenate(D_)
        keep = D_ <= np.quantile(D_, 1 - TRIM)
        Mn = sim(*umeyama(A_[keep], B_[keep]))
        if np.abs(Mn - M).max() < TOL: M = Mn; break
        M = Mn
    d = np.concatenate([T[l].query(apply(M, S[l]))[0] for l in labs])
    rms = float(np.sqrt((np.sort(d)[: int(len(d) * (1 - TRIM))] ** 2).mean()))
    return M, it, rms, labs


def known_answer(body, src, M_whole, tag, cover_m=0.005):
    """This body's own chest bones, truncated to the part HER labels actually cover, moved by a known
    similarity: the one-way fit must recover it. A superior-inferior band cut nothing -- her 369 mm
    field of view spans ribs 2-7 whole -- but her rib labels carry only 0.80-0.89 of each rib's
    surface, so the cut that matters is the one her segmentation makes, measured here per rib."""
    chest = {l: body[l] for l in FIT_LABELS}
    allv = np.vstack([v for v, _ in chest.values()]); c = allv.mean(0)
    # How much of each rib her label carries is measured by AREA RATIO, which depends only on the
    # scale and not on the pose. Testing proximity to her mapped label instead would cut away the
    # whole-torso misalignment itself (it kept a median 0.27-0.46 of each rib, against a measured
    # 0.80-0.89) and would ask the instrument to recover a similarity from a quarter of a rib.
    # A CT rib label stops short at the anterior costal-cartilage end, so the cut is taken there.
    scale = scale_of(M_whole)
    cut = {}
    for lab, (V, F) in chest.items():
        if lab not in src: continue
        share = float(np.clip(areas(*src[lab]).sum() * scale ** 2 / areas(V, F).sum(), 0.5, 1.0))
        a = areas(V, F); order = np.argsort(-V[F].mean(1)[:, 2])        # anterior (+z) end first
        drop = np.cumsum(a[order]) <= (1 - share) * a.sum()
        keepf = np.ones(len(F), bool); keepf[order[drop]] = False
        if keepf.sum() < 50: continue
        cut[lab] = (V, F[keepf])
    frac = {lab: round(float(areas(V, F).sum() / areas(*chest[lab]).sum()), 3) for lab, (V, F) in cut.items()}
    axis = np.array([1, -2, 1.5]); R = Rotation.from_rotvec(np.deg2rad(9) * axis / np.linalg.norm(axis)).as_matrix()
    Mt = sim(1.06, R, np.array([0.04, -0.03, 0.05]))                     # body -> "scan"
    scan = {lab: (apply(Mt, V), F) for lab, (V, F) in cut.items()}
    M, it, rms, _ = fit(scan, chest, list(cut))
    Minv = np.linalg.inv(Mt)
    terr = 1000 * float(np.linalg.norm(apply(M, apply(Mt, c[None]))[0] - c))
    Rf = M[:3, :3] / scale_of(M); Ri = Minv[:3, :3] / scale_of(Minv)
    rerr = float(np.degrees(np.arccos(np.clip((np.trace(Rf.T @ Ri) - 1) / 2, -1, 1))))
    serr = abs(scale_of(M) / scale_of(Minv) - 1)
    ok = terr <= KA_T_MM and rerr <= KA_R_DEG and serr <= KA_S
    say(f"  kept surface share after the cut: min {min(frac.values()):.2f}, median {np.median(list(frac.values())):.2f}")
    say(f"GATE 1 known answer: translation {terr:.3f} mm (<= {KA_T_MM}), rotation {rerr:.3f} deg (<= {KA_R_DEG}), "
        f"scale {100*serr:.3f}% (<= {100*KA_S:.0f}%), {it} iterations, trimmed residual {1000*rms:.2f} mm -> {'PASS' if ok else 'FAIL'}")
    return ok, dict(subject_defining_the_cut=tag, translation_mm=terr, rotation_deg=rerr, scale_rel=serr,
                    iterations=it, trimmed_residual_mm=1000 * rms, kept_share=frac)


def laterality(M, src, body):
    out = {}
    for side, other in (("left", "right"), ("right", "left")):
        P = np.vstack([apply(M, area_samples(*src[f"rib_{side}_{i}"], 2000)) for i in range(2, 8) if f"rib_{side}_{i}" in src])
        own = cKDTree(np.vstack([area_samples(*body[f"rib_{side}_{i}"], 20000) for i in range(2, 8)]))
        opp = cKDTree(np.vstack([area_samples(*body[f"rib_{other}_{i}"], 20000) for i in range(2, 8)]))
        out[side] = (float(np.median(own.query(P)[0])), float(np.median(opp.query(P)[0])))
    ok = all(a < b for a, b in out.values())
    say(f"GATE 2 laterality: left {1000*out['left'][0]:.1f} vs {1000*out['left'][1]:.1f} mm, "
        f"right {1000*out['right'][0]:.1f} vs {1000*out['right'][1]:.1f} mm -> {'PASS' if ok else 'FAIL'}")
    return ok, {k: dict(own_mm=1000 * v[0], opposite_mm=1000 * v[1]) for k, v in out.items()}


def held_out(M, M_whole, src, body):
    """median nearest-surface distance on the structures NOT fitted, chest-local vs whole-torso."""
    rows = {}
    for lab in HELD_OUT:
        if lab not in src or lab not in body: continue
        P = area_samples(*src[lab], 3000, 7)
        tree = cKDTree(area_samples(*body[lab], 30000, 11))
        rows[lab] = (float(np.median(tree.query(apply(M, P))[0])), float(np.median(tree.query(apply(M_whole, P))[0])))
    chest = float(np.median([v[0] for v in rows.values()])); whole = float(np.median([v[1] for v in rows.values()]))
    ok = chest <= whole
    say(f"GATE 3 held out ({len(rows)} structures): median nearest-surface distance chest-local {1000*chest:.2f} mm "
        f"vs whole-torso {1000*whole:.2f} mm -> {'PASS' if ok else 'FAIL'}")
    worse = sorted(((v[0] - v[1]) * 1000, l) for l, v in rows.items())[-3:]
    say(f"  worst three relative to whole-torso: " + ", ".join(f"{l} {d:+.1f} mm" for d, l in reversed(worse)))
    return ok, dict(chest_local_median_mm=1000 * chest, whole_torso_median_mm=1000 * whole,
                    per_structure_mm={l: dict(chest_local=1000 * v[0], whole_torso=1000 * v[1]) for l, v in rows.items()})


def breast_behind_wall(subject, side, M, workdir):
    """GATE 4: map her breast label with the new transform, mesh it, and measure the volume more
    than TRIM_DEPTH_M behind the muscular chest wall -- the trim gate, now testing the transform."""
    V, F = SEAT.read_obj(SRC / subject / "meshes" / f"breast_{side}.obj")
    Y = apply(M, V)
    workdir.mkdir(parents=True, exist_ok=True)
    import igl
    out = [np.asarray(o) for o in igl.qslim(np.ascontiguousarray(Y), np.ascontiguousarray(F), SEAT.DECIMATE_FACES) if hasattr(o, "shape")]
    U = next(o for o in out if o.dtype.kind == "f" and o.ndim == 2 and o.shape[1] == 3)
    G = next(o for o in out if o.dtype.kind in "iu" and o.ndim == 2 and o.shape[1] == 3)
    SEAT.write_obj(workdir / "breast.obj", U, G, [f"{subject} {side} breast under the chest-local transform"])
    msh = workdir / "breast.msh"
    if not msh.exists():
        cmd = [str(SEAT.FTETWILD), "-i", str(workdir / "breast.obj"), "-o", str(msh), "--no-binary",
               "-e", "1e-3", "-l", str(SEAT.MESH_LR), "--max-threads", "12"]
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        (workdir / "ftetwild.log").write_text(p.stdout[-3000:] + p.stderr[-3000:])
        if p.returncode != 0 or not msh.exists(): raise SystemExit(f"fTetWild failed on {subject} {side}")
    X, T = SEAT.read_msh_tets(msh)
    used = np.unique(T); remap = -np.ones(len(X), np.int64); remap[used] = np.arange(len(used)); X, T = X[used], remap[T]
    v = SEAT.tet_volumes(X, T); T[v < 0] = T[v < 0][:, [0, 2, 1, 3]]
    B = SEAT.boundary_faces(T)
    bV, bF, _ = SEAT.bed(side, near=X)
    idx, pen, _, _ = SEAT.penetration_of(X, B, bV, bF)
    deep = idx[pen > TRIM_DEPTH_M]
    total = SEAT.tet_volumes(X, T).sum()
    if not len(deep): return 0.0, float(pen.max() * 1e3), 0
    isdeep = np.zeros(len(X), bool); isdeep[deep] = True
    keep = ~isdeep[T].any(1)
    return float(1.0 - SEAT.tet_volumes(X, T[keep]).sum() / total), float(pen.max() * 1e3), int(len(deep))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subject", choices=SUBJECTS); ap.add_argument("--stage", default="all", choices=("all", "known-answer"))
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    body = body_meshes()
    subjects = [a.subject] if a.subject else list(SUBJECTS)
    # the known answer's cut band is her field of view, measured from her own ribs under the existing map
    src0 = subject_meshes(subjects[0])
    M_whole0 = np.array(json.loads((ROOT / "data/derived" / REGISTERED[subjects[0]] / "manifest.json").read_text())["transform"])
    ribs = np.vstack([apply(M_whole0, v) for l, (v, f) in src0.items() if l.startswith("rib_")])
    band = (float(ribs[:, 1].min()), float(ribs[:, 1].max()))
    say(f"her field of view, superior-inferior: {1000*(band[1]-band[0]):.0f} mm band "
        f"(it cuts none of ribs 2-7; her SEGMENTATION is what leaves them partial)")
    ka_ok, ka = known_answer(body, src0, M_whole0, subjects[0])
    (OUT / "known_answer.json").write_text(json.dumps(ka, indent=2) + "\n")
    if a.stage == "known-answer" or not ka_ok:
        if not ka_ok: say("GATE 1 FAILED: the instrument does not recover a known similarity; no subject is fitted.")
        return
    results = {}
    for subject in subjects:
        say(f"\n=== {subject} ===")
        src = subject_meshes(subject)
        M_whole = np.array(json.loads((ROOT / "data/derived" / REGISTERED[subject] / "manifest.json").read_text())["transform"])
        t0 = time.time(); M, it, rms, labs = fit(src, body, FIT_LABELS); secs = time.time() - t0
        say(f"  chest-local fit on {len(labs)} rib labels: {it} iterations, trimmed residual {1000*rms:.2f} mm, {secs:.0f} s")
        g2, lat = laterality(M, src, body)
        g3, ho = held_out(M, M_whole, src, body)
        dM = np.linalg.inv(M_whole) @ M
        dt = 1000 * float(np.linalg.norm(dM[:3, 3])); dr = float(np.degrees(np.arccos(np.clip((np.trace(dM[:3, :3] / scale_of(dM)) - 1) / 2, -1, 1))))
        say(f"  reported: chest-local differs from whole-torso by {dt:.1f} mm, {dr:.2f} deg, scale {scale_of(M):.4f} vs {scale_of(M_whole):.4f}")
        g4 = {}
        for side in ("left", "right"):
            frac, pmax, ndeep = breast_behind_wall(subject, side, M, OUT / subject / side)
            g4[side] = dict(volume_behind_wall=frac, deepest_mm=pmax, deep_vertices=ndeep)
            say(f"GATE 4 {side}: {100*frac:.3f}% of the breast more than {TRIM_DEPTH_M*1e3:.0f} mm behind the wall "
                f"(deepest {pmax:.1f} mm, {ndeep} vertices) -> {'PASS' if frac <= TRIM_VOLUME_TOL else 'FAIL'}")
        g4_ok = all(v["volume_behind_wall"] <= TRIM_VOLUME_TOL for v in g4.values())
        results[subject] = dict(transform=M.tolist(), iterations=it, trimmed_residual_mm=1000 * rms,
                                gate2_laterality=dict(passes=g2, **lat), gate3_held_out=dict(passes=g3, **ho),
                                gate4_behind_wall=dict(passes=g4_ok, **g4),
                                differs_from_whole_torso=dict(translation_mm=dt, rotation_deg=dr,
                                                              scale=scale_of(M), whole_torso_scale=scale_of(M_whole)),
                                passes=bool(g2 and g3 and g4_ok))
        (OUT / f"{subject}.json").write_text(json.dumps(results[subject], indent=2) + "\n")
    (OUT / "summary.json").write_text(json.dumps({"known_answer": ka, "subjects": results}, indent=2) + "\n")
    say("\n=== gate table ===")
    for s, r in results.items():
        say(f"  {s}: laterality {'pass' if r['gate2_laterality']['passes'] else 'FAIL'} | held out "
            f"{'pass' if r['gate3_held_out']['passes'] else 'FAIL'} | behind the wall "
            f"{100*r['gate4_behind_wall']['left']['volume_behind_wall']:.2f}% / "
            f"{100*r['gate4_behind_wall']['right']['volume_behind_wall']:.2f}% "
            f"{'pass' if r['gate4_behind_wall']['passes'] else 'FAIL'} -> {'PASS' if r['passes'] else 'FAIL'}")


if __name__ == "__main__": main()
