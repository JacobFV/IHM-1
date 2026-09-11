"""GATES Z and Y: known answers for the ANATOMICAL bed, which has no closed form.

On the anatomical bed with a non-rigid motion there is no analytic solution, and prescribing the whole
boundary as gate X does would remove the bed's role entirely. Two known answers need neither:

  GATE Z, PATH INDEPENDENCE: the same final configuration reached in 16 increments and in 64 must
  agree. Step count is exactly what varies association STALENESS while holding everything else fixed,
  and a stale-association error is path-dependent by construction, so Z points straight at it.

  GATE Y, REVERSIBILITY: driven forward by a non-rigid motion and then back, the return must reproduce
  the start. For frictionless unilateral contact on a hyperelastic body the quasi-static solution is a
  property of the configuration, not of the path to it.

TOLERANCE: the bed's facet scale, 1 mm -- the same measured quantity the threshold uses. The
association is discrete, so two step counts may legitimately select different facets and differ at
that scale; anything larger is path dependence the discretisation cannot excuse.

THE LIMIT, stated so a pass is not over-read: reversibility and path independence are NECESSARY, NOT
SUFFICIENT. A consistently wrong answer is reversible and path-independent too. They bound the
SOLVER's contribution to the error and say nothing about whether the seated breast is anatomically
right, which is judged by gates (a)-(d) and by nothing else.
"""
import sys, importlib.util, time
sys.path.insert(0, "/home/brandonin/Documents/IHM-1")
import numpy as np
from ihm.assembly.sliding_contact import SlidingRegion, seat_on_bed
from ihm.assembly.prescribed_deformation import lame

ROOT = "/home/brandonin/Documents/IHM-1"
D = f"{ROOT}/data/derived/female-breast-sliding-v1/s1159/left"
TOL_M = 1e-3

spec = importlib.util.spec_from_file_location("seat", f"{ROOT}/scripts/seat_breast_sliding.py")
sys.argv = ["x"]; S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
P = np.load(f"{D}/prepared.npz")
X, T, base = P["X"], P["T"], P["base"]
move = np.load(f"{D}/smoothed_move.npy")
mu, lam = lame(1000.0, 0.49)

def run(steps, phases, label):
    region = SlidingRegion(X, T, mu_pa=mu, lambda_pa=lam, density_kg_m3=950.0)
    closest, _ = S.bed_rays(P["bedV"], P["bedF"], P["ray_directions"], S.RAY_REACH_M)
    t0 = time.time()
    try:
        r = seat_on_bed(region, base, closest, load_steps=steps, gap_tol_m=S.GAP_TOL_M,
                        jump_limit_m=5e-4, assoc_tol_m=S.ASSOC_TOL_M, move=move, freeze_frames=True,
                        association="allornothing", lost_bed="hold", facet_m=S.GAP_TOL_M, phases=phases,
                        log=lambda m, flush=True: print(f"    [{label}] {m}", flush=True))
        u = r["displacement"]
        print(f"  {label}: completed, {len(r['steps'])} increments, refusals "
              f"{100*r['refusal_rate']:.1f}%, min J {r['steps'][-1]['min_J']:.4f}, {time.time()-t0:.0f} s")
        return u, r
    except Exception as e:
        print(f"  {label}: STALLED -- {type(e).__name__}: {str(e)[:90]} ({time.time()-t0:.0f} s)")
        return None, None

print("GATE Z: the same drive in 16 increments and in 64, on the anatomical bed")
u16, r16 = run(16, (1.0,), "16 increments")
u64, r64 = run(64, (1.0,), "64 increments")
if u16 is not None and u64 is not None:
    d = np.linalg.norm(u16 - u64, axis=1)
    print(f"  16 against 64: median {1000*np.median(d):.4f} mm, max {1000*d.max():.4f} mm "
          f"-> {'PASS' if d.max() <= TOL_M else 'FAIL'} (tolerance {1000*TOL_M:.0f} mm, the facet scale)")
else:
    print("  -> Z INCONCLUSIVE: a run did not complete, so there is no common final configuration.")

print("\nGATE Y: driven forward and back in one schedule")
uy, ry = run(16, (1.0, 0.0), "forward then back")
if uy is not None:
    d = np.linalg.norm(uy, axis=1)
    print(f"  return against the start: median {1000*np.median(d):.4f} mm, max {1000*d.max():.4f} mm "
          f"-> {'PASS' if d.max() <= TOL_M else 'FAIL'} (tolerance {1000*TOL_M:.0f} mm)")
else:
    print("  -> Y INCONCLUSIVE: the schedule did not complete.")

print("\nNECESSARY, NOT SUFFICIENT: a consistently wrong answer is reversible and path-independent too.")
print("These bound the solver's contribution to the error; whether the seated breast is anatomically")
print("right is judged by gates (a)-(d) and by nothing else.")
