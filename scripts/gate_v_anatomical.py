"""GATE V on the ANATOMICAL bed: is association-driven adaptive stepping what makes this regime
affordable?

On the plane V passes, but refusals were already zero there under the old rule, so that pass shows
the machinery is correct and discriminating -- not that it helps. This is the affordability question:
the same non-rigid seating drive, on the bed that refused 93.8% of updates under a RIGID drive and on
which gate Z's first leg has run past 42 minutes unfinished.

Reported beyond the gate: the rejection rate and how far the step must shrink to get a consistent
update; the wall-clock cost against Z's leg; and whether associations ACTUALLY UPDATE -- a run that
completes by shrinking to nothing and never updating is T-none wearing a different name, which is
exactly what the conjunction exists to catch.
"""
import sys, importlib.util, time, json
sys.path.insert(0, "/home/brandonin/Documents/IHM-1")
CAP_S = float(sys.argv[sys.argv.index("--cap") + 1]) if "--cap" in sys.argv else 1800.0
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
    if time.time() - start > CAP_S: raise Cap(f"wall-clock cap {CAP_S:.0f}s reached")

print(f"GATE V on the anatomical bed: {int((P['gap0'] < 0).sum())} held nodes, "
      f"the smoothed seating field (median travel {1000*np.median(move[P['gap0'] < 0]):.2f} mm), cap {CAP_S:.0f}s")
r, err = None, None
try:
    r = seat_on_bed(region, base, closest, load_steps=8, gap_tol_m=S.GAP_TOL_M, jump_limit_m=5e-4,
                    assoc_tol_m=S.ASSOC_TOL_M, move=move, freeze_frames=True, association="adaptive",
                    lost_bed="hold", facet_m=S.GAP_TOL_M, log=logger)
except Cap as e:
    err = str(e)
except Exception as e:
    err = f"{type(e).__name__}: {e}"
elapsed = time.time() - start
if r is not None:
    steps = r["steps"]; rej = sum(r["refusals"]); total = len(r["refusals"])
    print(f"\nGATE V (anatomical): COMPLETED to 1.0 in {elapsed:.0f}s, {len(steps)} accepted increments")
    print(f"  steps rejected {rej} of {total} ({100*rej/max(total,1):.1f}%), cutbacks {r['cutbacks']}")
    print(f"  associations actually updated: {len(r['lost_bed_per_association'])}")
    print(f"  smallest accepted increment: {min(s['fraction'] for s in steps[1:]) if len(steps)>1 else 1.0:.6f} of the drive")
    print(f"  min J at the end {steps[-1]['min_J']:.4f}")
else:
    print(f"\nGATE V (anatomical): DID NOT COMPLETE in {elapsed:.0f}s -- {err}")
    print("  the fraction reached and the rejection pattern are in the log above; a run that shrinks")
    print("  without updating its associations is T-none under another name, which is what the")
    print("  conjunction exists to catch.")
print(f"  against gate Z's current-stepping leg, which passed 42 minutes unfinished on the same drive")
