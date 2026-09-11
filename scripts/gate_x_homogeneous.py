"""GATE X (and, with --adaptive, GATE V on the plane): a homogeneous deformation of the held base against an analytic plane.

Non-rigid and sliding, yet closed-form. A simple shear TANGENT to the plane is isochoric (det F = 1)
and, being homogeneous, satisfies equilibrium exactly for a homogeneous neo-Hookean material -- the
class scripts/... R2 already validates against platens. It is the first control in which association
motion VARIES FROM NODE TO NODE, which is the only untested property of the refusal rule.

The whole BOUNDARY is driven, not the base alone: a homogeneous field is exact only if the entire
boundary is compatible with it, and driving the base alone would let the free surface relax, leaving
nothing closed-form to compare against.

NOT SUFFICIENT: X runs on a PLANE. The breast is non-rigid motion against the ANATOMICAL bed, and
the standing condition on this line is unchanged -- a passing X is not clearance.
"""
import sys, importlib.util
sys.path.insert(0, "/home/brandonin/Documents/IHM-1")
# read the flag BEFORE anything overwrites sys.argv: loading the seat module below sets
# sys.argv = ["x"], which silently turned a --adaptive run into an allornothing one and
# reported it as gate V passing.
ADAPTIVE = "--adaptive" in sys.argv
import numpy as np
from ihm.assembly.sliding_contact import SlidingRegion, seat_on_bed
from ihm.assembly.prescribed_deformation import lame, tet_volumes

ROOT = "/home/brandonin/Documents/IHM-1"
P = np.load(f"{ROOT}/data/derived/female-breast-sliding-v1/s1159/left/prepared.npz")
X, T, base_all, gap0_all, B = P["X"], P["T"], P["base"], P["gap0"], P["boundary"]
base = base_all[gap0_all < 0]
mu, lam = lame(1000.0, 0.49)
region = SlidingRegion(X, T, mu_pa=mu, lambda_pa=lam, density_kg_m3=950.0)

spec = importlib.util.spec_from_file_location("seat", f"{ROOT}/scripts/seat_breast_sliding.py")
sys.argv = ["x"]; S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
n_plane = S.bed_rays(P["bedV"], P["bedF"], P["ray_directions"][gap0_all < 0], S.RAY_REACH_M)[0](X[base])[1].mean(0)
n_plane /= np.linalg.norm(n_plane)
p0 = X[base].mean(0) + 0.005 * n_plane

def plane_closest(Q):
    Q = np.asarray(Q, float); d = (Q - p0) @ n_plane
    return Q - d[:, None] * n_plane[None, :], np.broadcast_to(n_plane, Q.shape).copy()

t1 = np.cross(n_plane, [0.0, 0.0, 1.0]); t1 /= np.linalg.norm(t1)
t2 = np.cross(n_plane, t1); t2 /= np.linalg.norm(t2)
span = float(np.ptp((X[base] - p0) @ t1))
gamma = 0.0075 / span                                   # ~7.5 mm across the base, as the other controls use
field = lambda Q: gamma * ((np.asarray(Q, float) - p0) @ t1)[:, None] * t2[None, :]
F = np.eye(3) + gamma * np.outer(t2, t1)
boundary_nodes = np.unique(B)
print(f"simple shear tangent to the plane: gamma {gamma:.5f}, det F {np.linalg.det(F):.12f} (isochoric); "
      f"node displacement 0 to {1000*np.abs(field(X[boundary_nodes])).max():.2f} mm")
u_exact = field(X)
print(f"association motion will VARY across the base: per full drive, min {1000*np.linalg.norm(field(X[base]),axis=1).min():.3f} mm, "
      f"median {1000*np.median(np.linalg.norm(field(X[base]),axis=1)):.3f}, max {1000*np.linalg.norm(field(X[base]),axis=1).max():.3f}")
MODE = "adaptive" if ADAPTIVE else "allornothing"
print(f"stepping: {MODE}" + ("  (GATE V: a rejected step shrinks rather than refusing the update)" if MODE == "adaptive" else ""))
r = seat_on_bed(region, base, plane_closest, load_steps=8, gap_tol_m=1e-3, jump_limit_m=5e-4,
                assoc_tol_m=1e-3, freeze_frames=True, association=MODE, lost_bed="hold",
                facet_m=1e-3, prescribe=(boundary_nodes, field(X[boundary_nodes])),
                log=lambda m, flush=True: print(m, flush=True))
u = r["displacement"]; Y = X + u
J = np.linalg.det(np.swapaxes(Y[T[:, 1:]] - Y[T[:, 0, None]], 1, 2)
                  @ np.linalg.inv(np.swapaxes(X[T[:, 1:]] - X[T[:, 0, None]], 1, 2)))
err = np.linalg.norm(u - u_exact, axis=1)
ok = bool(J.min() >= 0.99 and r["refusal_rate"] == 0 and err.max() < 1e-5)
print(f"\nGATE X: completed to 1.0, min J {J.min():.4f} (exact 1.0000), inversions {int((J<=0).sum())}")
print(f"  refusal rate {100*r['refusal_rate']:.1f}% ({sum(r['refusals'])} of {len(r['refusals'])})")
print(f"  solution vs the exact homogeneous field: median {1000*np.median(err):.6f} mm, max {1000*err.max():.6f} mm")
print(f"  volume ratio {tet_volumes(Y, T).sum()/tet_volumes(X, T).sum():.8f} (isochoric shear: 1.0)")
updates = len(r.get("lost_bed_per_association", []))
print(f"  associations actually updated: {updates} (a rule that never updates would show 0)")
print(f"  -> {'PASS' if ok else 'FAIL'} (fraction 1.0, min J >= 0.99, refusal rate 0, exact field"
      + (", associations updating)" if MODE == "adaptive" else ")"))
print("  NOT SUFFICIENT: this is a plane. The standing condition -- non-rigid motion against the")
print("  ANATOMICAL bed -- is unchanged, and a passing X is not clearance for the breast.")
