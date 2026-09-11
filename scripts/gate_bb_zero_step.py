"""GATE BB-PRIME: A ZERO-SIZED STEP, REPEATED AT ONE FIXED FRACTION, MUST REACH A FIXED POINT.

"A step of size zero must move nothing" was too strong, and the coordinator replaced it (f3b7c08)
after this script's first run: where re-association has moved the plane AWAY from a node, the node's
unilateral lower bound genuinely LOOSENS, and the body ought to relax into it. Forbidding that would
forbid contact from ever releasing. What must be true is not that the motion is zero but that it
CONVERGES -- repeated zero-sized steps at ONE fraction must reach a fixed point.

  decays  -> a one-time relaxation into a legitimately loosened feasible set
  persists -> the association and the solve are chasing each other at a FIXED load, which is the same
              non-convergence in a new place and is NOT fixed by the repair

THE BOUNDS half is unambiguous and is unchanged by the rewording: 0.000e+00 mm over 0 nodes with no
inversions, against 7.19 mm over 26 nodes and four inverted elements before the repair, repeatable
to the last digit. That is the defect closed as a test.

DECLARED BEFORE THE RUN, so it cannot be fitted to the numbers: the fixed point counts as reached
when a repeat moves less than 1 micrometre. That is a thousandth of the bed's 1 mm facet scale --
the smallest length at which this geometry means anything -- and it is not a tolerance chosen to be
passed. The whole sequence is printed either way.

NOTE ON READING THE FIRST RUN'S NUMBERS: 0.335 mm at fraction 0.0625 and 0.477 mm at 0.1094 are at
DIFFERENT fractions. They are not a convergence sequence and the rise between them says nothing
about whether either converges.

THE DEFECT. A held node's bound was set absolutely to `on_plane + gap0 + fraction*travel`, with
`on_plane = n.(c - X[base])` the displacement that puts the node on its association plane. The load
fraction scales `travel` and does not scale `on_plane`, so once the association had drifted, the
head of the next increment demanded that drift instantly however small the step. That single fact
accounts for gates R, R', R'', S, T, U, V and AA: every one of them was adjusting the term that was
already under control.

THE REPAIR. The bound is interpolated from where the node IS to where the constraint wants it, by
theta, the share of the REMAINING drive the increment consumes. At theta = 0 the bound is the node's
own coordinate, computed with the same einsum over the same frames the solver uses to form v_prev,
so the two agree bit for bit and the displacement is exactly zero rather than float dust. At
theta = 1 the bound is the constraint itself, so the endpoint of the drive is unchanged.

WHAT THIS GATE CANNOT TELL ME, said plainly: after the repair BB passes BY CONSTRUCTION -- theta = 0
multiplies the whole bracket. So BB on its own is no evidence that the repair is any good. It fixes
the statement of the defect so it cannot come back silently. The evidence is elsewhere: the
regressions must still pass, and the drive must get further than it did.
"""
import sys, importlib.util, time
sys.path.insert(0, "/home/brandonin/Documents/IHM-1")
CAP_S = float(sys.argv[sys.argv.index("--cap") + 1]) if "--cap" in sys.argv else 3000.0
import numpy as np
from ihm.assembly.sliding_contact import SlidingRegion, seat_on_bed
from ihm.assembly.prescribed_deformation import lame

ROOT = "/home/brandonin/Documents/IHM-1"
D = f"{ROOT}/data/derived/female-breast-sliding-v1/s1159/left"
spec = importlib.util.spec_from_file_location("seat", f"{ROOT}/scripts/seat_breast_sliding.py")
sys.argv = ["x"]; S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)

P = np.load(f"{D}/prepared.npz")
X, T, base = P["X"], P["T"], P["base"]
move = np.load(f"{D}/smoothed_move.npy")
mu, lam = lame(1000.0, 0.49)
region = SlidingRegion(X, T, mu_pa=mu, lambda_pa=lam, density_kg_m3=950.0)
closest, _ = S.bed_rays(P["bedV"], P["bedF"], P["ray_directions"], S.RAY_REACH_M)

start = time.time()
class Cap(RuntimeError): pass
def logger(m, flush=True):
    print(f"    [{time.time()-start:6.0f}s] {m}", flush=True)
    if time.time() - start > CAP_S: raise Cap(f"wall-clock cap {CAP_S:.0f}s")

rows = []
REPEATS = int(sys.argv[sys.argv.index("--repeats") + 1]) if "--repeats" in sys.argv else 6
FIXED_POINT_M = 1e-6                      # declared above, before the run: a thousandth of a facet
class Done(RuntimeError): pass


def zero_step(fraction, advance, u, aim):
    """REPEATS zero-sized steps at ONE fraction, each starting from the last one's answer. Declared:
    each call re-associates, which is the point -- the test is whether the association and the solve
    settle against each other at a fixed load, not whether one call is quiet."""
    print(f"  GATE BB-PRIME at fraction {fraction:.4f}: {REPEATS} zero-sized steps, each from the "
          "last one's answer", flush=True)
    seq, u_i = [], u
    for k in range(REPEATS):
        u_next, info, _, _, _ = advance(fraction, u_i, 0.0, aim)
        d = np.abs(u_next - u_i)
        interior = np.ones(len(d), bool); interior[base] = False
        seq.append(dict(k=k + 1, bound_m=info["bound_displacement_m"],
                        nodes=info["bound_nodes_over_1mm"], moved_m=float(d.max()),
                        base_m=float(d[base].max()), interior_m=float(d[interior].max()),
                        min_J=info["min_J"]))
        r = seq[-1]
        print(f"    repeat {r['k']}: bounds displaced {1000*r['bound_m']:.3e} mm over {r['nodes']} "
              f"nodes; the solve moved {1000*r['moved_m']:.6f} mm (base {1000*r['base_m']:.6f}, "
              f"interior {1000*r['interior_m']:.6f}); min J {r['min_J']:.4f}", flush=True)
        u_i = u_next
        if r["moved_m"] < FIXED_POINT_M: break
    rows.append(dict(fraction=fraction, seq=seq))
    raise Done("the gate has its sequence; the drive is not needed past this point")


print(f"GATE BB-PRIME: a zero-sized step repeated at one fixed fraction must reach a fixed point "
      f"(declared: under {1000*FIXED_POINT_M:g} mm)")
print("GATE BB: a step of size zero must move nothing, asked after every accepted step")
err = None
try:
    seat_on_bed(region, base, closest, load_steps=S.LOAD_STEPS, gap_tol_m=S.GAP_TOL_M,
                jump_limit_m=5e-4, move=move, assoc_tol_m=S.ASSOC_TOL_M,
                on_step=zero_step, log=logger)
except Done: pass
except Cap as e: err = str(e)
except Exception as e: err = f"{type(e).__name__}: {e}"

if not rows:
    print(f"\nGATE BB-PRIME: NOT TESTED -- the drive accepted no step" + (f": {err}" if err else ""))
    raise SystemExit(1)
seq = rows[0]["seq"]
worst_bound = max(r["bound_m"] for r in seq)
last = seq[-1]["moved_m"]
print(f"\nthe sequence at fraction {rows[0]['fraction']:.4f}, in mm: "
      + " -> ".join(f"{1000*r['moved_m']:.6f}" for r in seq))
if worst_bound > 0.0:
    print(f"GATE BB-PRIME: FAIL on the bounds -- a zero-sized step displaced the configuration by "
          f"up to {1000*worst_bound:.3e} mm; the requirement is exactly zero")
elif last < FIXED_POINT_M:
    print(f"GATE BB-PRIME: PASS -- the bounds displaced exactly 0.0 mm at every repeat, and the "
          f"motion fell to {1000*last:.6f} mm, under the declared {1000*FIXED_POINT_M:g} mm, in "
          f"{len(seq)} repeats. The relaxation is one-time: the association and the solve settle "
          "against each other at a fixed load.")
else:
    print(f"GATE BB-PRIME: FAIL -- the bounds displaced exactly 0.0 mm, but after {len(seq)} repeats "
          f"at ONE fraction the motion is still {1000*last:.6f} mm and has not reached a fixed "
          "point. The association and the solve are chasing each other at a fixed load, which the "
          "repair does not address.")
