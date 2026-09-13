#!/usr/bin/env python3
"""Is the ~19 mm anterior deficit the similarity fit's SINGLE SCALE?

Pre-registered in `docs/BODY_PARAMETERS.md`, 2026-09-13, and committed before this script ran.
The gate is quoted below verbatim and is not moved by anything printed here.

`register_female_torso.py` fits a SIMILARITY -- rotation + one uniform scale + translation
(`umeyama`). One scale cannot represent a thorax that differs from this body's by a different
ratio in AP depth than in ML width or SI height. The 34 correspondences are 24 ribs, 2
clavicles, 12 thoracic vertebrae and the sternum; that population's spread is dominated by ML
and SI, so the fitted scale is set almost entirely by those two axes. If the CT subjects'
thoraces are proportionally shallower in AP, the scale overshoots AP, the mapped chest wall
lands too far posterior, and the breasts land with it.

THE GATE, fixed in advance:
  CONFIRMED if s_AP is the SMALLEST of the three scales in >= 3 of 4 subjects AND the predicted
    anterior shortfall at the sternum is >= 50% of that subject's measured correction in >= 3
    of 4 subjects.
  REFUTED otherwise. A partial result is recorded as REFUTED with the fraction stated, never as
    "partly confirmed".

THREE KNOWN ANSWERS, and the third exists because the first two can both pass for a broken
instrument.

  0. THE REBUILD. The uniform scale computed here must reproduce the scale `register.log` wrote
     for the same subject's centroid fit, to < 0.5%. Without this the fit is on some other set
     of points and every number below is about a population the registration never saw. The
     dropped-label list is read from that same log rather than recomputed, so the population is
     the registration's by construction.

  1. AP-ONLY SQUASH. Apply a known anisotropic map -- AP 0.900, ML and SI 1.000, composed with a
     random rotation and translation -- to a subject's own centroids. The anisotropic fit must
     recover 0.900 to < 0.5%; the uniform fit must NOT, returning one scale near the geometric
     mean (~0.965).

  2. ISOTROPIC INPUT, and this is the one that matters. An anisotropic fit has three parameters
     where the similarity has one, so it will always report SOME anisotropy. Under an isotropic
     synthetic map (all three scales 1.070) it must return three scales equal to within 0.5%.
     If it manufactures anisotropy from an isotropic input then every anisotropy it reports on
     real data is its own, and the measurement is void. Known answer 1 alone would pass for a
     fit that merely has more freedom.

  3. THE AXIS CONVENTION. Every number here is "the AP scale", so the identity of the AP axis is
     load-bearing and is asserted, not assumed: T1 -> T12 must displace predominantly along the
     axis called superoinferior, and left clavicle -> right clavicle predominantly along the one
     called lateral. Whichever remains is anterior. This repository has burned days on AP/ML
     selection failures where the instrument passed and the selection failed.

THE MODEL. The anisotropic map is  b ~= S R a + t  with S = diag(s_x, s_y, s_z) in THIS BODY's
frame -- scales applied after the rotation, so they are scales along anatomical axes rather than
along the subject's scanner axes. Fitted by Levenberg-Marquardt over 9 parameters (rotation
vector, three log-scales, translation), initialised at the umeyama solution, so it can never do
worse than the similarity.

CAVEAT that travels with any use of this: four CT subjects registered onto a male-derived body.
This measures the disagreement between that body's thorax and a registered subject's, which is a
statement about the registration, not about any subject's anatomy.
"""
import argparse, ast, gzip, json, re, sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = ROOT / "data/derived/female-torso-totalsegmentator-v1"
REG = {"s0790": "female-torso-registered-v1", "s1067": "female-torso-registered-s1067-v1",
       "s1159": "female-torso-registered-s1159-v1", "s0970": "female-torso-registered-s0970-v1"}
# the measured anterior correction per subject, from the eight-breast table in
# docs/BODY_PARAMETERS.md (mean of the two sides). Quoted, not recomputed here.
MEASURED_ANTERIOR_MM = {"s0790": 19.96, "s0970": 24.73, "s1067": 13.70, "s1159": 16.11}
ORD = ["first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth",
       "tenth", "eleventh", "twelfth"]
SHORTFALL_FRACTION_MIN = 0.50          # the gate, fixed in advance
N_SUBJECTS_MIN = 3                     # of 4, both halves of the gate


def read_obj(p):
    V, F = [], []
    for line in open(p):
        if line.startswith("v "):
            V.append([float(x) for x in line.split()[1:4]])
        elif line.startswith("f "):
            F.append([int(x.split("/")[0]) - 1 for x in line.split()[1:4]])
    return np.asarray(V, float), np.asarray(F, np.int64)


def areas(V, F):
    t = V[F]
    return np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2


def area_centroid(V, F):
    a = areas(V, F)
    return (a[:, None] * V[F].mean(1)).sum(0) / a.sum()


def entity_mesh(e):
    g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
    return np.asarray(g["positions"], float).reshape(-1, 3), np.asarray(g["indices"], np.int64).reshape(-1, 3)


def umeyama(a, b):
    """Exactly `register_female_torso.py`'s, copied so the comparison is against the same code."""
    ca, cb = a.mean(0), b.mean(0)
    x, y = a - ca, b - cb
    u, s, vt = np.linalg.svd(x.T @ y)
    d = np.diag([1.0, 1.0, np.sign(np.linalg.det(vt.T @ u.T))])
    r = vt.T @ d @ u.T
    scale = float((s @ np.diag(d)).sum() / (x ** 2).sum())
    return scale, r, cb - scale * r @ ca


def aniso_fit(a, b, s0, R0, t0):
    """b ~= diag(s) R a + t, nine parameters, initialised at the similarity solution."""
    r0 = Rotation.from_matrix(R0).as_rotvec()
    p0 = np.concatenate([r0, np.log([s0, s0, s0]), t0])

    def resid(p):
        R = Rotation.from_rotvec(p[:3]).as_matrix()
        S = np.exp(p[3:6])
        return ((a @ R.T) * S + p[6:9] - b).ravel()

    sol = least_squares(resid, p0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    R = Rotation.from_rotvec(sol.x[:3]).as_matrix()
    return np.exp(sol.x[3:6]), R, sol.x[6:9], float(np.sqrt((sol.fun ** 2).reshape(-1, 3).sum(1).mean()))


def sim_rms(a, b, s, R, t):
    return float(np.sqrt((((a @ R.T) * s + t - b) ** 2).sum(1).mean()))


def dropped_from_log(sid):
    """The population the registration actually fitted, taken from the log it wrote.

    Rebuilding the drop rule here would be a reconstruction that happens to be the right
    length -- the failure mode CLAUDE.md records for the THINGS-EEG2 image order. The log is
    the dataset that defines it.
    """
    log = (ROOT / "data/derived" / REG[sid] / "register.log").read_text()
    n_kept = int(re.search(r"(\d+) bone correspondences", log).group(1))
    m = re.search(r"dropped (\[[^\]]*\])", log)
    # s0790 was registered WITHOUT --whole-labels-only and has no dropped line: it kept all 39.
    # The count check below is what catches a mis-read either way, so this is not an assumption.
    return (set() if m is None else set(ast.literal_eval(m.group(1)))), n_kept


def logged_centroid_scale(sid):
    log = (ROOT / "data/derived" / REG[sid] / "register.log").read_text()
    return float(re.search(r"centroid fit\s*:\s*scale\s*([0-9.]+)", log).group(1))


def build_correspondences(sid):
    pairs = {}
    for side in ("left", "right"):
        for i in range(12):
            pairs[f"rib_{side}_{i + 1}"] = [f"{side} {ORD[i]} rib"]
        pairs[f"clavicula_{side}"] = [f"{side} clavicle"]
    for i in range(12):
        pairs[f"vertebrae_T{i + 1}"] = [f"{ORD[i]} thoracic vertebra"]
    pairs["sternum"] = ["manubrium", "body of sternum", "xiphoid process"]

    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    by_name = {}
    for e in ents:
        if e["role"] == "rigid_bone" and e["id"].startswith("body-bp3d-"):
            by_name.setdefault(e["name"], []).append(e)

    dropped, n_kept = dropped_from_log(sid)
    labels = [l for l in sorted(pairs) if l not in dropped]
    if len(labels) != n_kept:
        raise SystemExit(f"{sid}: rebuilt {len(labels)} labels, register.log fitted {n_kept}")

    src, dst, named = [], [], {}
    for lab in labels:
        Vs, Fs = read_obj(SRC_ROOT / sid / "meshes" / f"{lab}.obj")
        parts = []
        for nm in pairs[lab]:
            cand = by_name.get(nm, [])
            if len(cand) != 1:
                raise SystemExit(f"{lab}: expected one entity named {nm!r}, found {len(cand)}")
            parts.append(entity_mesh(cand[0]))
        off = np.cumsum([0] + [len(v) for v, _ in parts[:-1]])
        Vb = np.vstack([v for v, _ in parts])
        Fb = np.vstack([f + o for (_, f), o in zip(parts, off)])
        src.append(area_centroid(Vs, Fs))
        dst.append(area_centroid(Vb, Fb))
        named[lab] = (src[-1], dst[-1])
    return np.array(src), np.array(dst), labels, named


def identify_axes(named):
    """Assert which column is lateral, superoinferior, anterior -- in THIS BODY's frame.

    Not assumed. T1 -> T12 is a superoinferior displacement; left clavicle -> right clavicle is
    a lateral one. Whichever axis neither claims is anterior. Every 'AP scale' below depends on
    this being right, and AP/ML selection has failed five times on this line.
    """
    def body(lab):
        return named[lab][1]
    si = np.abs(body("vertebrae_T12") - body("vertebrae_T1"))
    ml = np.abs(body("clavicula_right") - body("clavicula_left"))
    ax_si, ax_ml = int(np.argmax(si)), int(np.argmax(ml))
    ok_si = si[ax_si] > 2.0 * np.sort(si)[-2]
    ok_ml = ml[ax_ml] > 2.0 * np.sort(ml)[-2]
    if ax_si == ax_ml or not (ok_si and ok_ml):
        raise SystemExit(f"axis identification FAILED: SI {si} -> {ax_si}, ML {ml} -> {ax_ml}")
    ax_ap = ({0, 1, 2} - {ax_si, ax_ml}).pop()
    return ax_ml, ax_si, ax_ap, si, ml


def known_answers(src, dst):
    rng = np.random.default_rng(0)
    print("KNOWN ANSWERS")
    ok = True

    # 0 is checked per subject in main(); 3 is checked in identify_axes().
    for name, S_true in (("1  AP-ONLY SQUASH   ", np.array([1.0, 1.0, 0.900])),
                         ("2  ISOTROPIC INPUT  ", np.array([1.070, 1.070, 1.070]))):
        R_t = Rotation.from_rotvec(np.deg2rad(11) * np.array([1, -2, 3.0]) / np.sqrt(14)).as_matrix()
        t_t = np.array([0.12, -0.31, 0.05])
        b = (src @ R_t.T) * S_true + t_t
        s_u, R_u, t_u = umeyama(src, b)
        S_a, R_a, t_a, rms_a = aniso_fit(src, b, s_u, R_u, t_u)
        rms_u = sim_rms(src, b, s_u, R_u, t_u)
        ratio = S_a / S_a.max()
        print(f"  {name}true scales {S_true}  ->  anisotropic fit {np.round(S_a, 5)}")
        print(f"     uniform fit one scale {s_u:.5f} (geometric mean of the truth "
              f"{float(np.prod(S_true) ** (1/3)):.5f}); "
              f"rms uniform {rms_u*1e3:9.4f} mm, anisotropic {rms_a*1e3:9.4f} mm")
        if S_true[0] != S_true[2]:                       # known answer 1
            got = float(S_a[2] / S_a[0])
            hit = abs(got - 0.900) / 0.900 < 0.005
            uni_blind = abs(s_u - 0.900) / 0.900 > 0.005
            print(f"     anisotropic recovers AP/ML {got:.5f} vs 0.900 -> {'PASS' if hit else 'FAIL'}; "
                  f"uniform does NOT recover it -> {'PASS' if uni_blind else 'FAIL'}")
            ok &= hit and uni_blind
        else:                                            # known answer 2 -- the one that matters
            spread = float(S_a.max() / S_a.min() - 1.0)
            hit = spread < 0.005
            print(f"     spread across the three fitted scales {spread:.2e} (< 5e-3) -> "
                  f"{'PASS' if hit else 'FAIL'}   [an isotropic input must NOT yield anisotropy]")
            ok &= hit
    print()
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out/registration_anisotropy.json")
    a = ap.parse_args()

    print(__doc__.split("CAVEAT")[0].strip().splitlines()[0])
    print(f"pre-registered gate: s_AP smallest in >= {N_SUBJECTS_MIN} of 4 AND predicted shortfall "
          f">= {SHORTFALL_FRACTION_MIN:.0%} of the measured correction in >= {N_SUBJECTS_MIN} of 4\n")

    rows, ka_done = [], False
    for sid in sorted(REG):
        src, dst, labels, named = build_correspondences(sid)
        ax_ml, ax_si, ax_ap, si, ml = identify_axes(named)
        if not ka_done:
            print(f"KNOWN ANSWER 3  AXIS CONVENTION (from {sid}'s body-side centroids, identical "
                  f"for all four)")
            print(f"  T1 -> T12 displacement |{np.round(si*1e3, 1)}| mm -> axis {ax_si} is SUPEROINFERIOR")
            print(f"  L -> R clavicle        |{np.round(ml*1e3, 1)}| mm -> axis {ax_ml} is LATERAL")
            print(f"  by elimination axis {ax_ap} is ANTERIOR                                  PASS\n")
            if not known_answers(src, dst):
                sys.exit("a known answer FAILED; the anisotropy measurement is void")
            ka_done = True

        s_u, R_u, t_u = umeyama(src, dst)
        logged = logged_centroid_scale(sid)
        rel = abs(s_u - logged) / logged
        print(f"KNOWN ANSWER 0  REBUILD {sid}: {len(labels)} correspondences, uniform scale "
              f"{s_u:.4f} vs register.log {logged:.4f}, {rel:.2%} "
              f"-> {'PASS' if rel < 0.005 else 'FAIL'}")
        if rel >= 0.005:
            sys.exit(f"{sid}: the rebuilt correspondences are not the registration's")

        S_a, R_a, t_a, rms_a = aniso_fit(src, dst, s_u, R_u, t_u)
        rms_u = sim_rms(src, dst, s_u, R_u, t_u)

        # the pre-registered predictor: (s_AP - s_uniform) x the AP distance from the fit's
        # centroid to the sternum, in the body frame.
        src_stern, dst_stern = named["sternum"]
        d_ap = float(abs((dst_stern - dst.mean(0))[ax_ap]))
        pred = float((S_a[ax_ap] - s_u) * d_ap * 1e3)
        # and the direct version, reported beside it as a check on the formula: where each map
        # actually puts the source sternum centroid along the anterior axis.
        z_u = float(((src_stern @ R_u.T) * s_u + t_u)[ax_ap])
        z_a = float(((src_stern @ R_a.T) * S_a + t_a)[ax_ap])
        direct = (z_a - z_u) * 1e3
        rows.append(dict(subject=sid, n=len(labels), uniform_scale=s_u, logged_scale=logged,
                         s_ml=float(S_a[ax_ml]), s_si=float(S_a[ax_si]), s_ap=float(S_a[ax_ap]),
                         rms_uniform_mm=rms_u * 1e3, rms_aniso_mm=rms_a * 1e3,
                         ap_is_smallest=bool(S_a[ax_ap] == S_a.min()),
                         sternum_ap_lever_mm=d_ap * 1e3,
                         predicted_shortfall_mm=pred, direct_shortfall_mm=direct,
                         measured_anterior_mm=MEASURED_ANTERIOR_MM[sid],
                         fraction=pred / MEASURED_ANTERIOR_MM[sid]))

    print(f"\n{'subject':8s} {'n':>3s} {'s_uniform':>10s} {'s_ML':>8s} {'s_SI':>8s} {'s_AP':>8s} "
          f"{'AP min?':>8s} {'rms uni':>8s} {'rms ani':>8s}")
    for r in rows:
        print(f"{r['subject']:8s} {r['n']:3d} {r['uniform_scale']:10.4f} {r['s_ml']:8.4f} "
              f"{r['s_si']:8.4f} {r['s_ap']:8.4f} {'YES' if r['ap_is_smallest'] else 'no':>8s} "
              f"{r['rms_uniform_mm']:7.2f}m {r['rms_aniso_mm']:7.2f}m")

    print(f"\n{'subject':8s} {'lever mm':>9s} {'predicted':>10s} {'direct':>8s} {'measured':>9s} "
          f"{'fraction':>9s}")
    for r in rows:
        print(f"{r['subject']:8s} {r['sternum_ap_lever_mm']:9.1f} {r['predicted_shortfall_mm']:+10.2f} "
              f"{r['direct_shortfall_mm']:+8.2f} {r['measured_anterior_mm']:+9.2f} "
              f"{r['fraction']:9.1%}")

    n_ap_min = sum(r["ap_is_smallest"] for r in rows)
    n_frac = sum(r["fraction"] >= SHORTFALL_FRACTION_MIN for r in rows)
    confirmed = n_ap_min >= N_SUBJECTS_MIN and n_frac >= N_SUBJECTS_MIN
    print(f"\n  s_AP is the smallest of the three in {n_ap_min} of 4 (need >= {N_SUBJECTS_MIN})")
    print(f"  the predicted shortfall reaches {SHORTFALL_FRACTION_MIN:.0%} of the measured "
          f"correction in {n_frac} of 4 (need >= {N_SUBJECTS_MIN})")
    print(f"\n  VERDICT: " + ("CONFIRMED -- the single scale is what puts the breasts too deep, and "
                              "an\n  anisotropic registration is the remedy. The 2.61-4.03 mm residual "
                              "shape mismatch\n  is untouched by this and stays exactly where it is."
                              if confirmed else
                              "REFUTED on the gate as written. The uniform scale is NOT what puts the\n"
                              "  breasts ~19 mm too deep; the next candidate must be named rather than "
                              "this one\n  rescued, and the gate is not moved to fit what came back."))

    p = ROOT / a.out
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(dict(
        gate=dict(n_subjects_min=N_SUBJECTS_MIN, shortfall_fraction_min=SHORTFALL_FRACTION_MIN),
        confirmed=bool(confirmed), n_ap_smallest=int(n_ap_min), n_fraction_met=int(n_frac),
        rows=rows,
        caveat="four CT subjects registered onto a male-derived body; a statement about the "
               "registration, not about any subject's anatomy"), indent=2) + "\n")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
