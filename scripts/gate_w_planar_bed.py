"""GATE W: an analytic PLANAR bed, driven by a rigid translation TANGENT to it.

The plane makes the answer closed-form: a tangential rigid translation is zero-energy, the body does
not deform, every association slides by exactly the translation distance, and no association can
teleport because the bed has no features to jump between. It is the minimal control that exercises
BED-FOLLOWING -- the regime every rigid-drive control on this line has missed.
"""
import sys, importlib.util; sys.path.insert(0, '/home/brandonin/Documents/IHM-1')
import numpy as np
from ihm.assembly.sliding_contact import SlidingRegion, seat_on_bed
from ihm.assembly.prescribed_deformation import lame, tet_volumes

P = np.load("/home/brandonin/Documents/IHM-1/data/derived/female-breast-sliding-v1/s1159/left/prepared.npz")
X, T, base_all, gap0_all = P["X"], P["T"], P["base"], P["gap0"]
base = base_all[gap0_all < 0]                      # the same held set the other controls drive
mu, lam = lame(1000.0, 0.49)
region = SlidingRegion(X, T, mu_pa=mu, lambda_pa=lam, density_kg_m3=950.0)

# the plane: through 5 mm in FRONT of the base's centroid, so every base node sits behind it (held)
spec = importlib.util.spec_from_file_location("seat", "/home/brandonin/Documents/IHM-1/scripts/seat_breast_sliding.py")
sys.argv = ["x"]; S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
closest0, _ = S.bed_rays(P["bedV"], P["bedF"], P["ray_directions"][gap0_all < 0], S.RAY_REACH_M)
n_plane = closest0(X[base])[1].mean(0); n_plane /= np.linalg.norm(n_plane)
p0 = X[base].mean(0) + 0.005 * n_plane

def plane_closest(Q):
    Q = np.asarray(Q, float)
    d = (Q - p0) @ n_plane
    return Q - d[:, None] * n_plane[None, :], np.broadcast_to(n_plane, Q.shape).copy()

# a translation TANGENT to the plane, of the same 7.48 mm the other controls use
t = np.cross(n_plane, [0.0, 0.0, 1.0]); t /= np.linalg.norm(t)
d = 0.00748 * t
print(f"plane normal [{n_plane[0]:+.3f},{n_plane[1]:+.3f},{n_plane[2]:+.3f}], tangential drive "
      f"{1000*np.linalg.norm(d):.2f} mm along [{t[0]:+.3f},{t[1]:+.3f},{t[2]:+.3f}]; "
      f"n.d = {1000*abs(n_plane @ d):.2e} mm (tangent to machine precision)")
r = seat_on_bed(region, base, plane_closest, load_steps=8, gap_tol_m=1e-3, jump_limit_m=5e-4,
                assoc_tol_m=1e-3, rigid_m=d, drive_full=True, freeze_frames=True,
                association="allornothing", lost_bed="hold", facet_m=1e-3,
                log=lambda m, flush=True: print(m, flush=True))
u = r["displacement"]; Y = X + u
J = np.linalg.det(np.swapaxes(Y[T[:, 1:]] - Y[T[:, 0, None]], 1, 2)
                  @ np.linalg.inv(np.swapaxes(X[T[:, 1:]] - X[T[:, 0, None]], 1, 2)))
err = np.linalg.norm(u - d, axis=1)
per_step = np.linalg.norm(d) / 8
moves = [st.get("association_moved_m", float("nan")) for st in r["steps"]]
print(f"\nGATE W: completed to 1.0, min J {J.min():.4f}, inversions {int((J<=0).sum())}")
print(f"  refusal rate {100*r['refusal_rate']:.1f}% ({sum(r['refusals'])} of {len(r['refusals'])})")
print(f"  association motion per step: {[f'{1000*m:.4f}' for m in moves]} mm against the translation "
      f"distance {1000*per_step:.4f} mm")
print(f"  solution vs the exact tangential translation: median {1000*np.median(err):.6f} mm, max {1000*err.max():.6f} mm")
print(f"  volume ratio {tet_volumes(Y, T).sum()/tet_volumes(X, T).sum():.8f}")
ok = (J.min() >= 0.99 and r['refusal_rate'] == 0 and
      all(abs(m - per_step) <= 1e-5 for m in moves if np.isfinite(m)))
print(f"  -> {'PASS' if ok else 'FAIL'} (min J >= 0.99, refusal rate 0, association motion = translation distance)")
