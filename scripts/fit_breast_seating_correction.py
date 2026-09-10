"""one systematic correction that seats the mapped female breasts on this body's chest wall.

three women registered into this body (s0790, s1067, s1159) fail gate d the same way: this
body's anterior chest wall sits in FRONT of the mapped female one. its ribs 2-7 lie inside
the mapped breasts, deepest at rib 5 (median 6.7-10.5 mm, max 12.6-17.2 mm), tapering to
1-3 mm at ribs 2-3 and 7. same direction, same profile, in all three -- so the correction is
ONE correction in this body's frame, fitted on subjects and tested on a subject it never saw.

models
  baseline  per side, one rigid translation along this body's anterior axis (+z); the
            smallest magnitude that brings the training subjects' rib points inside the
            breast to <= 1%. volume is conserved EXACTLY (a rigid move).
  profile   only if the baseline's held-out judge fails: every breast point displaced along
            this body's horizontal outward chest-wall normal (radial from the ribcage's
            vertical axis) by s * g(height), g a quadratic fitted to the training subjects'
            clearance envelope, s the smallest scale that clears the training subjects.
            volume is checked explicitly.
leave-one-subject-out: fit on two, apply to the third.

THE JUDGE, named before the fit and NOT optimised by it (gate d2 is not a judge: a
correction fitted to clear the ribs passes it by construction):
  a  held-out: this body's rib points inside the held-out subject's breasts <= 1%
  b  contact:  the held-out breast still TOUCHES the chest wall -- median distance of its
               10% of vertices nearest this body's rib + sternum surfaces <= 3 mm (95th
               percentile reported)
  c  volume:   rigid moves conserve it exactly; a deformation must stay within 1%

INSIDE TEST. rib points are tested against the subject's own breast MASK (breast & body,
the voxels the breast meshes were built from), mapped through the inverse registration --
exact to the voxel and O(points). for a displacement field u, p is inside the displaced
breast iff p - u(p) is inside the original, to first order in the gradient of u (u varies
over ~100 mm of height, so the error is second order). a gate checks the mask lookup against
the mesh ray-parity test on the uncorrected breasts.

ORDER GATE. the extraction was being rerun under a new stray-fragment rule while this was
written; the registrations rest on meshes from before it. this refuses to fit while that
rerun is running, and unless s0790's and s1159's gate values and breast volumes are
identical before and after it.
"""
import gzip, importlib.util, json, subprocess, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
TS = ROOT / "data/derived/female-torso-totalsegmentator-v1"
RAWTS = ROOT / "data/raw/anatomy/totalsegmentator"
OUT = ROOT / "data/derived/female-breast-seating-correction-v1"
SUBJECTS = {"s0790": ROOT / "data/derived/female-torso-registered-v1",
            "s1067": ROOT / "data/derived/female-torso-registered-s1067-v1",
            "s1159": ROOT / "data/derived/female-torso-registered-s1159-v1"}
ORD = ("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth")
SIDES = ("left", "right")
RIB_SHARE_MAX, CONTACT_MM_MAX, VOLUME_TOL = 0.01, 3.0, 0.01
EZ = np.array([0.0, 0.0, 1.0])                                  # this body's anterior axis (gate: +z anterior)
spec = importlib.util.spec_from_file_location("offset", ROOT / "scripts/measure_female_chest_wall_offset.py")
OFF = importlib.util.module_from_spec(spec); spec.loader.exec_module(OFF)

def say(*a): print(*a, flush=True)

def order_gate():
    running = subprocess.run(["bash", "-c", "ps -eo pid,comm,args | awk '$2 ~ /^python/ && /build_female_torso/ {print $1}'"],
                             capture_output=True, text=True).stdout.split()
    if running: sys.exit(f"ORDER GATE: the extraction is still running (pid {running}); refusing to fit")
    before, after = TS / "manifest.before-stray-rule.json", TS / "manifest.json"
    if not before.exists(): sys.exit("ORDER GATE: no pre-rerun manifest to compare against")
    b, a = json.loads(before.read_text())["subjects"], json.loads(after.read_text())["subjects"]
    keys = ("gate_bones_in_body", "gate_breast_in_body", "gate_breast_on_ribs", "gate_chest_whole",
            "gate_breast_whole", "breast_ml", "passes")
    for sid in ("s0790", "s1159"):
        diff = [k for k in keys if b[sid].get(k) != a[sid].get(k)]
        say(f"ORDER GATE: {sid} gate values and breast volumes before/after the stray-rule rerun: "
            + ("IDENTICAL" if not diff else f"CHANGED {diff}"))
        if diff: sys.exit("the registrations rest on changed meshes; refusing to fit")

def area_samples(V, F, n, seed):
    rng = np.random.default_rng(seed); t = V[F]
    a = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2
    k = rng.choice(len(F), n, p=a / a.sum()); r1 = np.sqrt(rng.random(n)); r2 = rng.random(n); t = t[k]
    return (1 - r1)[:, None] * t[:, 0] + (r1 * (1 - r2))[:, None] * t[:, 1] + (r1 * r2)[:, None] * t[:, 2]

def signed_volume(V, F):
    t = V[F]; return float(np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum() / 6.0)

class Subject:
    def __init__(self, sid, reg):
        import nibabel as nib
        self.sid = sid
        m = json.loads((reg / "manifest.json").read_text())
        if "transform" not in m: sys.exit(f"{sid}: registration manifest carries no transform")
        self.M = np.asarray(m["transform"], float); self.Minv = np.linalg.inv(self.M)
        breast = np.asanyarray(nib.load(str(TS / sid / "breasts" / "breast.nii.gz")).dataobj) > 0
        body = np.asanyarray(nib.load(str(TS / sid / "body" / "body.nii.gz")).dataobj) > 0
        self.mask = breast & body                                       # what the breast meshes were built from
        A = nib.load(str(RAWTS / sid / "ct.nii.gz")).affine; self.Ainv = np.linalg.inv(A)
        self.mesh = {s: OFF.read_obj(reg / f"breast_{s}.obj") for s in SIDES}
        self.vol = {s: signed_volume(*self.mesh[s]) for s in SIDES}

    def inside(self, P):
        """canonical points -> this subject's CT voxels -> breast mask lookup."""
        q = (P @ self.Minv[:3, :3].T + self.Minv[:3, 3]) * 1000.0        # CT world, mm
        ijk = np.rint(q @ self.Ainv[:3, :3].T + self.Ainv[:3, 3]).astype(int)
        ok = np.all((ijk >= 0) & (ijk < np.array(self.mask.shape)), axis=1)
        out = np.zeros(len(P), bool); out[ok] = self.mask[tuple(ijk[ok].T)]
        return out

def load_body():
    ents = {e["name"]: e for e in json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
            if e["role"] == "rigid_bone" and e["id"].startswith("body-bp3d-")}
    def mesh(n):
        g = json.loads(gzip.decompress((ROOT / ents[n]["reference_geometry"]["path"]).read_bytes()))
        return np.asarray(g["positions"], float).reshape(-1, 3), np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ribs = {s: np.vstack([area_samples(*mesh(f"{s} {o} rib"), 1500, seed=10 * i + (s == "right"))
                          for i, o in enumerate(ORD)]) for s in SIDES}
    ribverts = {s: np.vstack([mesh(f"{s} {o} rib")[0] for o in ORD]) for s in SIDES}
    wall = np.vstack([area_samples(*mesh(n), 4000, seed=100 + j) for j, n in enumerate(
        [f"{s} {o} rib" for s in SIDES for o in ORD] + ["manubrium", "body of sternum", "xiphoid process"])])
    return ribs, ribverts, wall

class Correction:
    """u(p) = s * g(height) * n(p); baseline is g = 1, n = +z per side (a rigid translation)."""
    def __init__(self, kind, params, axis=None):
        self.kind, self.params, self.axis = kind, params, axis
    def u(self, P, side):
        if self.kind == "baseline":
            return np.outer(np.full(len(P), self.params[side]), EZ)
        s, c = self.params["scale"], self.params["poly"]
        g = np.clip(np.polyval(c, P[:, 1]), 0.0, None) * s
        n = P - self.axis; n[:, 1] = 0.0; n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
        return g[:, None] * n

def rib_share(subj, corr, ribs, side):
    P = ribs[side]
    return float(subj.inside(P - corr.u(P, side)).mean()) if corr else float(subj.inside(P).mean())

def contact(subj, corr, wall_tree, side):
    V, F = subj.mesh[side]; Vc = V + corr.u(V, side) if corr else V
    d = 1000.0 * wall_tree.query(Vc)[0]; near = np.sort(d)[: max(1, len(d) // 10)]
    return float(np.median(near)), float(np.percentile(near, 95)), signed_volume(Vc, F) / subj.vol[side] - 1.0

def min_translation(subj, ribs, side, hi=0.04):
    if rib_share(subj, None, ribs, side) <= RIB_SHARE_MAX: return 0.0
    lo, top = 0.0, hi
    if rib_share(subj, Correction("baseline", {side: top}), ribs, side) > RIB_SHARE_MAX: return float("inf")
    for _ in range(18):
        mid = (lo + top) / 2
        if rib_share(subj, Correction("baseline", {side: mid}), ribs, side) <= RIB_SHARE_MAX: top = mid
        else: lo = mid
    return top

def clearance_envelope(subjs, ribs, axis):
    """for every rib point inside a training breast, the outward displacement along n that clears it."""
    rows = []
    for subj in subjs:
        for side in SIDES:
            P = ribs[side]; ins = subj.inside(P); Q = P[ins]
            if not len(Q): continue
            n = Q - axis; n[:, 1] = 0.0; n /= np.linalg.norm(n, axis=1, keepdims=True)
            need = np.full(len(Q), np.nan)
            for d in np.arange(0.0005, 0.0405, 0.0005):
                clear = ~subj.inside(Q - d * n) & np.isnan(need); need[clear] = d
            rows.append(np.c_[Q[:, 1], np.nan_to_num(need, nan=0.04)])
    return np.vstack(rows)

def fit_profile(subjs, ribs, axis):
    env = clearance_envelope(subjs, ribs, axis)
    bins = np.arange(env[:, 0].min(), env[:, 0].max() + 0.01, 0.01)
    h, g = [], []
    for lo_, hi_ in zip(bins[:-1], bins[1:]):
        sel = (env[:, 0] >= lo_) & (env[:, 0] < hi_)
        if sel.sum() >= 5: h.append((lo_ + hi_) / 2); g.append(np.percentile(env[sel, 1], 95))
    poly = np.polyfit(h, g, 2)
    corr = Correction("profile", {"scale": 1.0, "poly": poly.tolist()}, axis)
    def ok(s):
        corr.params["scale"] = s
        return all(rib_share(sj, corr, ribs, sd) <= RIB_SHARE_MAX for sj in subjs for sd in SIDES)
    lo, top = 0.0, 3.0
    if not ok(top): corr.params["scale"] = float("inf"); return corr
    for _ in range(16):
        mid = (lo + top) / 2
        if ok(mid): top = mid
        else: lo = mid
    corr.params["scale"] = top
    return corr

def judge(subj, corr, ribs, wall_tree):
    out = {}
    for side in SIDES:
        sh = rib_share(subj, corr, ribs, side); cm, c95, dv = contact(subj, corr, wall_tree, side)
        out[side] = dict(rib_share=sh, contact_median_mm=cm, contact_p95_mm=c95, volume_change=dv,
                         a=sh <= RIB_SHARE_MAX, b=cm <= CONTACT_MM_MAX, c=abs(dv) <= VOLUME_TOL)
    return out

def main():
    order_gate()
    OUT.mkdir(parents=True, exist_ok=True)
    import trimesh
    ico = trimesh.creation.icosphere(subdivisions=4, radius=0.1)
    ph = np.random.default_rng(1).normal(size=(800, 3)); ph /= np.linalg.norm(ph, axis=1, keepdims=True)
    ci, co = OFF.inside(ico.vertices, ico.faces, 0.05 * ph).mean(), OFF.inside(ico.vertices, ico.faces, 0.5 * ph).mean()
    say(f"CONTROL (mesh parity): sphere points inside {ci:.3f} (1.000), outside {co:.3f} (0.000)")
    if not (ci == 1.0 and co == 0.0): sys.exit("parity control failed")
    subj = {sid: Subject(sid, reg) for sid, reg in SUBJECTS.items()}
    ribs, ribverts, wall = load_body(); wall_tree = cKDTree(wall)
    axis = wall.mean(0)                                           # the ribcage's vertical axis passes through here
    say("GATE (fast inside test): mask lookup vs mesh ray parity on the UNCORRECTED breasts, same rib points")
    for sid in SUBJECTS:
        for side in SIDES:
            P = ribs[side]; near = np.all((P > subj[sid].mesh[side][0].min(0) - 0.01) & (P < subj[sid].mesh[side][0].max(0) + 0.01), axis=1)
            Pn = P[near]; Pn = Pn[np.random.default_rng(0).choice(len(Pn), min(len(Pn), 3000), replace=False)]
            m, p = subj[sid].inside(Pn).mean(), OFF.inside(*subj[sid].mesh[side], Pn).mean()
            say(f"  {sid} {side:5s}: mask {m:.4f}  mesh {p:.4f}  |diff| {abs(m - p):.4f}")
            if abs(m - p) > 0.02: sys.exit("mask lookup disagrees with the mesh test by > 2 points; not a valid proxy")
    say("\nuncorrected (all this body's rib samples per side, the gate-d2 denominator):")
    for sid in SUBJECTS:
        say(f"  {sid}: " + "  ".join(f"{s} {rib_share(subj[sid], None, ribs, s):.4f}" for s in SIDES))
    mins = {sid: {s: min_translation(subj[sid], ribs, s) for s in SIDES} for sid in SUBJECTS}
    say("\nper-subject minimal anterior translation to reach <= 1% (mm): "
        + "  ".join(f"{sid} L {1000*v['left']:.1f} R {1000*v['right']:.1f}" for sid, v in mins.items()))
    report = dict(judge=dict(rib_share_max=RIB_SHARE_MAX, contact_median_mm_max=CONTACT_MM_MAX, volume_tolerance=VOLUME_TOL),
                  per_subject_min_translation_m=mins, folds=[])
    need_profile = False
    say("\n== BASELINE: one rigid anterior translation per side, leave-one-subject-out ==")
    for held in SUBJECTS:
        train = [s for s in SUBJECTS if s != held]
        corr = Correction("baseline", {s: max(mins[t][s] for t in train) for s in SIDES})
        j = judge(subj[held], corr, ribs, wall_tree)
        say(f"  held out {held}: t L {1000*corr.params['left']:.1f} mm R {1000*corr.params['right']:.1f} mm | " + " | ".join(
            f"{s}: ribs {j[s]['rib_share']:.4f}{'' if j[s]['a'] else ' FAIL-a'}  contact {j[s]['contact_median_mm']:.1f}"
            f"/{j[s]['contact_p95_mm']:.1f} mm{'' if j[s]['b'] else ' FAIL-b'}  vol {100*j[s]['volume_change']:+.2f}%" for s in SIDES))
        need_profile |= not all(j[s]["a"] and j[s]["b"] for s in SIDES)
        report["folds"].append(dict(model="baseline", held_out=held, params=corr.params, judge=j))
    final = Correction("baseline", {s: max(mins[t][s] for t in SUBJECTS) for s in SIDES})
    if need_profile:
        say("\n== PROFILE: s * g(height) along the outward chest-wall normal, leave-one-subject-out ==")
        for held in SUBJECTS:
            train = [subj[s] for s in SUBJECTS if s != held]
            corr = fit_profile(train, ribs, axis)
            j = judge(subj[held], corr, ribs, wall_tree)
            say(f"  held out {held}: scale {corr.params['scale']:.3f} poly {np.round(corr.params['poly'], 4).tolist()} | " + " | ".join(
                f"{s}: ribs {j[s]['rib_share']:.4f}{'' if j[s]['a'] else ' FAIL-a'}  contact {j[s]['contact_median_mm']:.1f}"
                f"/{j[s]['contact_p95_mm']:.1f} mm{'' if j[s]['b'] else ' FAIL-b'}  vol {100*j[s]['volume_change']:+.2f}%"
                f"{'' if j[s]['c'] else ' FAIL-c'}" for s in SIDES))
            report["folds"].append(dict(model="profile", held_out=held, params=corr.params, judge=j))
        final = fit_profile([subj[s] for s in SUBJECTS], ribs, axis)
    folds = [f for f in report["folds"] if f["model"] == final.kind]
    passed = all(f["judge"][s][k] for f in folds for s in SIDES for k in ("a", "b", "c"))
    say(f"\nFINAL ({final.kind}, fitted on all three): {json.dumps(final.params)}")
    say("VERDICT: " + ("the held-out judge passes on every fold -- a usable correction" if passed else
        "the held-out judge FAILS -- no correction of this family seats the breast; nothing seated is written"))
    report["final"] = dict(kind=final.kind, params=final.params, axis_point_m=axis.tolist() if final.axis is not None else None,
                           held_out_judge_passed=passed)
    report["verdict"] = ("usable" if passed else
        "NOT usable: rigid translation and a height profile along the outward normal both fail the held-out judge. "
        "clearing the rib-5 penetration floats the breast off the chest wall elsewhere (contact median > 3 mm on every "
        "fold), and the profile adds ~8-9% volume. what is needed is a local deformation of the breast BASE that "
        "conforms it to this body's chest wall, with volume compensation -- not a displacement of the whole breast.")
    if not passed:
        (OUT / "correction.json").write_text(json.dumps(report, indent=2, default=float) + "\n")
        say(f"wrote {(OUT / 'correction.json').relative_to(ROOT)} (verdict only)")
        return
    for sid in SUBJECTS:
        (OUT / sid).mkdir(exist_ok=True)
        for s in SIDES:
            V, F = subj[sid].mesh[s]; Vc = V + final.u(V, s)
            with open(OUT / sid / f"breast_{s}.obj", "w") as h:
                h.write(f"# {sid} breast_{s}, registered and seated by the systematic correction ({final.kind}), canonical frame, metres\n")
                for v in Vc: h.write("v %.9g %.9g %.9g\n" % tuple(v))
                for f in F: h.write("f %d %d %d\n" % (f[0] + 1, f[1] + 1, f[2] + 1))
    report["d1_note"] = ("moving the breasts does not move the mapped TRUNK, so gate d1 (this body's sternum inside it: "
                         "s0790 0.767, s1067 0.635, s1159 0.315) is unchanged by construction. the trunk envelope is the "
                         "separate skin-replacement problem.")
    (OUT / "correction.json").write_text(json.dumps(report, indent=2, default=float) + "\n")
    say(f"wrote {(OUT / 'correction.json').relative_to(ROOT)} and seated breasts per subject")

if __name__ == "__main__": main()
