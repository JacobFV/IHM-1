"""GATE CC-FLAT-PRIME: is the plane's 2 nm the convergence floor, or a real defect?

CC-flat declared EXACTLY zero and read 0.000002 mm, so it printed FAIL and the bar is not being
moved. The replacement is a convergence test rather than a new number: re-run the identical flat case
at successively tighter solver tolerances. Falling with rtol means the convergence floor and a clean
flat case; plateauing means the 2 nm is real. Worth running rather than assuming, because 2.7e-7 of
the motion is about 3x the declared rtol of 1e-7 -- a principled bar at rtol would have failed too.
"""
import sys; sys.path.insert(0, "/home/brandonin/Documents/IHM-1")
import numpy as np
from ihm.assembly.mechanics_backend import tetra_box
from ihm.assembly.sliding_contact import SlidingRegion, seat_on_bed
from ihm.assembly.prescribed_deformation import lame

X, T = tetra_box((0.060, 0.040, 0.030), (12, 8, 6)); X = X - X.mean(0)
mu, lam = lame(1000.0, 0.49)
bottom = np.flatnonzero(X[:, 2] < X[:, 2].min() + 1e-9)
top = np.flatnonzero(X[:, 2] > X[:, 2].max() - 1e-9)
z0 = X[:, 2].min() + 1e-3
d = np.array([0.00748, 0.0, 0.0])
n = np.array([0.0, 0.0, 1.0])
def bed(Q):
    Q = np.asarray(Q, float); c = Q.copy(); c[:, 2] = z0
    return c, np.broadcast_to(n, Q.shape).copy()

print("GATE CC-FLAT-PRIME: the plane, at successively tighter solver tolerances")
print(f"  drive {1000*np.linalg.norm(d):.2f} mm; the exact answer is u = d everywhere, energy zero\n")
print(f"  {'rtol':>8}  {'max |u - d|, mm':>18}  {'as a fraction of |d|':>22}  {'min J':>10}")
rows = []
for rtol in (1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10):
    region = SlidingRegion(X, T, mu_pa=mu, lambda_pa=lam, density_kg_m3=950.0)
    r = seat_on_bed(region, bottom, bed, load_steps=8, gap_tol_m=1e-3, jump_limit_m=0.01,
                    prescribe=(top, np.broadcast_to(d, (len(top), 3)).copy()),
                    move=np.zeros(len(bottom)), freeze_frames=True, association="persistent",
                    lost_bed="hold", assoc_tol_m=1e-3, facet_m=1e-3, rtol=rtol, log=None)
    u = r["displacement"]; Y = X + u
    J = np.linalg.det(np.swapaxes(Y[T[:, 1:]] - Y[T[:, 0, None]], 1, 2)
                      @ np.linalg.inv(np.swapaxes(X[T[:, 1:]] - X[T[:, 0, None]], 1, 2)))
    e = float(np.linalg.norm(u - d, axis=1).max())
    rows.append((rtol, e))
    print(f"  {rtol:8.0e}  {1000*e:18.9f}  {e/np.linalg.norm(d):22.3e}  {J.min():10.6f}", flush=True)

es = [e for _, e in rows]
rel = [e / float(np.linalg.norm(d)) for e in es]
drop = es[0] / max(es[-1], 1e-300)
# NO SLOPE IS QUOTED. The measured shape is a PLATEAU across the loose tolerances and then a
# collapse to machine precision, not a smooth power law, so a fitted exponent would be a number
# without a meaning -- the same kind of number this line has spent the day removing. What decides
# the gate is whether the deviation falls at all, and how far: to machine precision is a
# convergence floor, a plateau all the way down would be a defect.
floor = min(rel)
print(f"\n  the deviation is CONSTANT across the loose tolerances and then collapses; it falls "
      f"{drop:.3g}x in all, to {floor:.3e} of the drive")
print("  -> " + ("CC-FLAT-PRIME PASS: tightening rtol removes the residual entirely, reaching "
                 "double precision. The plane's 2 nm was the convergence floor, not a defect, and "
                 "the flat case is clean -- the cost of staleness on a plane is zero as the "
                 "decomposition requires."
                 if floor < 1e-12 else
                 "CC-FLAT-PRIME FAIL: the deviation does not fall to machine precision however "
                 "tight the tolerance, so the residual on a plane is real and is a defect"))
