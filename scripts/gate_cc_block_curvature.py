"""GATE CC: what does persistent association COST, and is the cost a bug or the price of curvature?

R2's rigid-translation control degrades under persistent association -- min J 1.000 -> 0.636 by
fraction 0.75 -- and the pre-repair code degrades too (min J 0.710 at fraction 0.25 against the
repaired run's 0.920), so the repair is not the cause. R2's earlier pass, platens 0.00000% and
cylinder 6.382 um, was measured WITHOUT persistent association and stands only for that
configuration; it is re-scoped, not retired, and is NOT quoted here.

THE DECOMPOSITION, which is what makes this a known answer rather than a tolerance:

  CC-flat    a PLANE bed. A plane's normal never rotates, so a stale association is not stale in any
             way that matters: the constraint direction at the new position is the constraint
             direction at the old one. The cost must be EXACTLY ZERO. No tolerance attaches. Any
             deviation at all is a defect, not a limitation.

  CC-curved  a CYLINDER of radius R, the same block, the same drive. Moving d along the surface
             rotates the closest point's normal by d/R, so a step of d_step leaves the node off the
             true surface by the sagitta R(1 - cos(d_step/R)) ~ d_step^2 / 2R. Over N steps that
             accumulates to about N * d_step^2 / 2R = d^2 / 2RN. Measured against predicted:
             matching means persistence has a known, quantified price on curved geometry;
             much larger means something beyond staleness is happening.

DEVIATION FROM THE GATE AS WORDED, declared rather than quietly taken: the coordinator's CC-flat is
"the platens", whose 0.00000% was a COMPRESSION compared against the exact homogeneous field. I have
built the flat arm as the same rigid tangential drive as the curved arm, on a plane, so that the two
arms differ in CURVATURE AND NOTHING ELSE. A pair that differs in two things cannot decompose
anything. The old platens number is therefore not comparable to this one and is not quoted.

WHY THE DRIVE IS TANGENT, AND WHY THAT TOOK FOUR TRIES. With travel = 0 the contact face is pinned
normally, n.u = 0, so the rigid translation is admissible ONLY if d is tangent. Tilting d out of the
tangent -- which was the fix for attempt 2, where a tangential drive applied to the contact face
itself moved nothing -- stretches the block instead, and the plane read 2.602647 mm. The tangential
drive now comes from the PRESCRIBED FAR FACE rather than from the normal constraint, so the reason
for the tilt is gone with it.

On a plane, u = d is then the unique zero-energy admissible solution and the error must be exactly
zero. On a cylinder the contact face must hold its gap while following a surface that curves away,
so the same translation forces a real deformation of order the sagitta.

THE SWEEP is the real test, not any single comparison: R is varied over a decade, and the error must
fall as 1/R. One number matching one prediction can be a coincidence; a slope cannot.
"""
import sys, time
sys.path.insert(0, "/home/brandonin/Documents/IHM-1")
import numpy as np
from ihm.assembly.mechanics_backend import tetra_box
from ihm.assembly.sliding_contact import SlidingRegion, seat_on_bed
from ihm.assembly.prescribed_deformation import lame, tet_volumes

DELTA_M = 1e-3            # how far the bed sits inside the block's face, so those nodes are held
DRIVE_M = 0.00748         # the same 7.48 mm every control on this line uses
DIRECTION = np.array([1.0, 0.0, 0.0])   # TANGENT to the bed at the contact: see the note below
STEPS = 8
RADII = (0.40, 0.20, 0.10, 0.05)


def plane_bed(z0):
    n = np.array([0.0, 0.0, 1.0])
    def f(Q):
        Q = np.asarray(Q, float)
        c = Q.copy(); c[:, 2] = z0
        return c, np.broadcast_to(n, Q.shape).copy()
    return f


def cylinder_bed(z0, R):
    """Axis along y at z = z0 - R, so the surface's top is at z = z0 with outward normal +z there --
    the R -> infinity limit of this is exactly plane_bed(z0)."""
    axis_z = z0 - R
    def f(Q):
        Q = np.asarray(Q, float)
        w = np.stack([Q[:, 0], Q[:, 2] - axis_z], 1)          # radial part, y is along the axis
        r = np.linalg.norm(w, axis=1, keepdims=True)
        u = w / np.maximum(r, 1e-30)
        c = np.stack([R * u[:, 0], Q[:, 1], axis_z + R * u[:, 1]], 1)
        n = np.stack([u[:, 0], np.zeros(len(Q)), u[:, 1]], 1)
        return c, n
    return f


def run(bed, label, R, frozen=True):
    """The contact face slides along the bed, dragged TANGENTIALLY THROUGH THE ELASTIC BLOCK by the
    far face, which is prescribed the translation d. That detail is the whole construction: under a
    normal-only constraint the tangential direction is undetermined, so a tangential drive applied to
    the contact face itself moves nothing (attempt 1) and a tilted one moves only along the normal
    (attempt 2, |u - d| came back as exactly the tangential component of d). The far face gives the
    contact nodes something to be pulled by, and the bed then has to be followed."""
    X, T = tetra_box((0.060, 0.040, 0.030), (12, 8, 6))
    X = X - X.mean(0)
    mu, lam = lame(1000.0, 0.49)
    region = SlidingRegion(X, T, mu_pa=mu, lambda_pa=lam, density_kg_m3=950.0)
    bottom = np.flatnonzero(X[:, 2] < X[:, 2].min() + 1e-9)     # the contact face, on the bed
    top = np.flatnonzero(X[:, 2] > X[:, 2].max() - 1e-9)        # the driven face
    z0 = X[:, 2].min() + DELTA_M                                # the bed sits INSIDE the contact face
    d = DRIVE_M * DIRECTION / np.linalg.norm(DIRECTION)
    t0 = time.time()
    try:
        r = seat_on_bed(region, bottom, bed(z0) if R is None else bed(z0, R), load_steps=STEPS,
                        # THE JUMP LIMIT MUST NOT BIND, or this gate measures nothing. At 0.5 mm
                        # against a per-step slide of 0.935 mm the association was REFUSED every
                        # step and reported 'moved 0.0000 mm' -- which reads as settled, not as
                        # refused. Persistent association was not being exercised at all, and
                        # frozen and re-linearised were bit-identical for that reason.
                        gap_tol_m=1e-3, jump_limit_m=0.01,
                        prescribe=(top, np.broadcast_to(d, (len(top), 3)).copy()),
                        # travel ZERO. Left to its default it is -gap0, so the contact face would
                        # also be seated 1 mm while the far face translates -- a real deformation
                        # that the rigid translation cannot match, which is what made the plane
                        # read 1.402927 mm on attempt 3 in BOTH the repaired and the pre-repair
                        # code. The reference was wrong, not the solver.
                        move=np.zeros(len(bottom)),
                        freeze_frames=frozen, association="persistent", lost_bed="hold",
                        # a re-linearised pass only happens if the step is NOT declared settled,
                        # and at assoc_tol_m = 1 mm it always is: frozen and re-linearised came
                        # back bit-identical, which is how the void A/B announced itself.
                        assoc_tol_m=1e-3 if frozen else 1e-9,
                        facet_m=1e-3, log=None)
    except Exception as e:
        print(f"  {label}: DID NOT COMPLETE -- {type(e).__name__}: {e}", flush=True)
        return None
    u = r["displacement"]; Y = X + u
    if np.abs(u).max() < 1e-9:
        print(f"  {label}: NOTHING MOVED -- the drive is not driving; this case is void", flush=True)
        return None
    J = np.linalg.det(np.swapaxes(Y[T[:, 1:]] - Y[T[:, 0, None]], 1, 2)
                      @ np.linalg.inv(np.swapaxes(X[T[:, 1:]] - X[T[:, 0, None]], 1, 2)))
    err = np.linalg.norm(u - d, axis=1)
    print(f"  {label}: min J {J.min():.6f}, inversions {int((J <= 0).sum())}, "
          f"max |u - d| {1000*err.max():.6f} mm, median {1000*np.median(err):.6f} mm, "
          f"passes/step {min(st['passes'] for st in r['steps'])}-{max(st['passes'] for st in r['steps'])}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    passes = [st["passes"] for st in r["steps"]]
    return dict(label=label, R=R, minJ=float(J.min()), inv=int((J <= 0).sum()),
                err=float(err.max()), u=u, passes=passes)


print("GATE CC: the cost of persistent association, decomposed into flat and curved")
print(f"  block 60 x 40 x 30 mm, {STEPS} steps, far face translated {1000*DRIVE_M:.2f} mm, "
      f"contact face on a bed {1000*DELTA_M:.1f} mm inside it")
print("\nCC-flat -- a plane. The rigid translation satisfies every constraint exactly (n.u = n.d for")
print("any n, and a plane's normal never rotates), so the cost must be EXACTLY zero:")
flat = run(plane_bed, "plane (R = infinity)", None)

print("\nCC-curved -- a cylinder, FROZEN frames against RE-LINEARISED ones. Their difference is the")
print("staleness and nothing else: same geometry, same drive, same association mode.")
rows = []
for R in RADII:
    a = run(cylinder_bed, f"cylinder R = {1000*R:.0f} mm, frozen        ", R, frozen=True)
    b = run(cylinder_bed, f"cylinder R = {1000*R:.0f} mm, re-linearised ", R, frozen=False)
    if a and b:
        diff = float(np.abs(a["u"] - b["u"]).max())
        pred = DRIVE_M ** 2 / (2 * R * STEPS)
        rows.append(dict(R=R, diff=diff, pred=pred))
        print(f"    -> staleness = |frozen - re-linearised| {1000*diff:.6f} mm against the predicted "
              f"d^2/2RN {1000*pred:.6f} mm, ratio {diff/pred:.3f}", flush=True)

print("\n" + "=" * 100)
if flat is None or len(rows) < 2:
    print("GATE CC: INCOMPLETE -- an arm did not finish, so neither half is readable.")
    raise SystemExit(1)
flat_ok = flat["err"] == 0.0 and flat["inv"] == 0
print(f"CC-flat: max |u - d| {1000*flat['err']:.6f} mm, min J {flat['minJ']:.6f} -> "
      + ("PASS, exactly zero as a plane requires" if flat_ok else
         "FAIL -- a plane's normal never rotates, so any deviation here is a defect and not a cost"))
ratios = [r["diff"] / r["pred"] for r in rows]
sl = np.polyfit(np.log([1.0 / r["R"] for r in rows]), np.log([r["diff"] for r in rows]), 1)[0]
print(f"CC-curved: measured/predicted " + ", ".join(f"{x:.3f}" for x in ratios)
      + f"; staleness scales as (1/R)^{sl:.3f}, and the sagitta argument says 1.000")
print("  -> " + ("the price of persistence on curved geometry, quantified: a limitation to carry"
                 if 0.5 <= min(ratios) and max(ratios) <= 2.0 and abs(sl - 1.0) <= 0.25 else
                 "NOT explained by staleness alone -- the measured cost does not track d^2/2RN"))
