"""THE ZERO-INCREMENT PROBE on the GLOBAL association -- the one thing today's work did not touch.

The seating drive stalls with 13 and 2 inverted elements at the START of an increment, reached from
states that had already been ACCEPTED. Today established only where those inversions do not come
from: not the mesh, not the initial bound configuration, not the first increment, not the clip, not
the linear response -- all of which read zero from the pristine state.

What was never separated is the increment from the RE-ASSOCIATION. Every cut-back re-associates
before it re-solves, so shrinking the step never tested the step alone. This asks the drive for a
step of size ZERO at the fraction it has already accepted. A zero-sized step asks for no new travel,
so if it still inverts, the increment is not involved at all and the re-association is; if it
succeeds, the increment is involved after all and the picture I am carrying is wrong.

PREDICTION, on the record before the run: the zero-sized step fails, with a finite bound
displacement, because re-association moves the bounds under an already-deformed configuration.
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

out = {}


def on_stall(fraction, advance, u):
    """The whole point of the script. advance(fraction, u) at the fraction ALREADY ACCEPTED is a
    step of size zero: the held nodes' targets are unchanged, so nothing is asked to travel."""
    out["fraction"] = fraction
    print(f"\n  STALLED at fraction {fraction:.4f}; the body has moved at most "
          f"{1000 * np.abs(u).max():.2f} mm", flush=True)
    try:
        advance(fraction, u)
        out["zero"] = "SUCCEEDS -- a zero-sized step is fine, so the failure DOES need a finite increment"
    except Exception as e:
        out["zero"] = f"FAILS TOO -- {type(e).__name__}: {e}"
    print(f"  ZERO-INCREMENT PROBE: {out['zero']}", flush=True)
    # And again, to see whether the zero-sized step is repeatable or drifts: two identical requests
    # from the same accepted state must give the same answer, for the same reason a query must.
    try:
        advance(fraction, u)
        out["zero2"] = "SUCCEEDS"
    except Exception as e:
        out["zero2"] = f"FAILS -- {e}"
    print(f"  ZERO-INCREMENT PROBE, repeated from the same accepted state: {out['zero2']}", flush=True)


print("THE ZERO-INCREMENT PROBE on the global association, production settings "
      f"(load_steps={S.LOAD_STEPS}, jump_limit 0.5 mm)")
try:
    r = seat_on_bed(region, base, closest, load_steps=S.LOAD_STEPS, gap_tol_m=S.GAP_TOL_M,
                    jump_limit_m=5e-4, move=move, assoc_tol_m=S.ASSOC_TOL_M,
                    on_stall=on_stall, log=logger)
    print(f"\nthe drive COMPLETED in {time.time()-start:.0f}s, min J {r['minimum_jacobian']:.4f} "
          "-- there was no stall to probe")
except Cap as e:
    print(f"\nno stall reached before the cap: {e}")
except Exception as e:
    print(f"\nthe drive stopped: {type(e).__name__}: {e}")
if "zero" in out:
    print(f"\nANSWER at fraction {out['fraction']:.4f}: a step of size zero {out['zero']}")
    print(f"  repeated: {out['zero2']}")
else:
    print("\nANSWER: none -- the drive never reached a stall, so there was nothing to probe.")
