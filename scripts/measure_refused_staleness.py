"""THE STALENESS, IN MILLIMETRES, WITH NOTHING ASSUMED.

The per-node jump limit freezes 31-37% of held nodes at every step, and the log never showed it:
"association moved 0.4999 mm" is the maximum over ACCEPTED updates, so a refused node contributes
nothing to it. This accumulates, per node, exactly how far each association WOULD have moved and did
not -- summed over the drive, in millimetres. No curvature model, no d^2/2R, nothing fitted.

  a few tenths over the whole drive -> the 31-37% is a bookkeeping fact and the anatomical seating
                                      results survive with a caveat
  millimetres                      -> those results were computed against constraints pointing where
                                      the bed no longer is, and are re-run rather than re-read

The refused nodes are refused BECAUSE their associations want to move furthest, so they are the
worst-case population by construction. That is a reason the number COULD be large, not a prediction
that it is, and it is why the distribution is reported rather than a single figure: the median over
held nodes says what the typical constraint suffers, the median over REFUSED nodes says what the
affected population suffers, and the maximum says how bad one constraint can get. Quoting only one
of the three would misrepresent the other two.

Accounting, stated so it can be checked: a refusal inside a step that was later CUT BACK is not
counted -- the accumulator is saved and restored with the association, so only refusals on accepted
steps contribute. Nodes whose ray finds no bed at all are counted separately as a count, never as a
distance, because "no bed to move to" is not a displacement.
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
c0, n0 = closest(X[base]); gap0 = np.einsum('ij,ij->i', n0, X[base] - c0)
held = gap0 < 0
NH = int(held.sum())

start_t = time.time()
class Cap(RuntimeError): pass
def logger(m, flush=True):
    if time.time() - start_t > CAP_S: raise Cap(f"wall-clock cap {CAP_S:.0f}s")

seen = []


def report(fraction, advance, u, aim, info):
    t = info["refused_travel_m"][held] * 1000.0          # mm, accumulated, held nodes only
    k = info["refused_steps"][held]
    ever = k > 0
    seen.append(dict(fraction=fraction, t=t.copy(), k=k.copy()))
    print(f"  fraction {fraction:.4f}: accumulated refused association travel over held nodes -- "
          f"median {np.median(t):.4f} mm, p90 {np.percentile(t, 90):.4f} mm, max {t.max():.4f} mm",
          flush=True)
    print(f"      {int(ever.sum())} of {NH} held nodes ({100*ever.mean():.1f}%) have been refused at "
          f"least once, up to {int(k.max())} times; over those nodes alone the median accumulated "
          f"travel is {np.median(t[ever]) if ever.any() else 0.0:.4f} mm", flush=True)
    # THE NUMBER THAT ACTUALLY MEASURES STALENESS: how far each constraint is from the bed NOW.
    w = info["wanted_now_m"][held] * 1000.0
    live = w > 0
    seen[-1]["w"] = w.copy()
    print(f"      INSTANTANEOUS, the distance each refused constraint is from the bed right now: "
          f"{int(live.sum())} nodes refused this step, median {np.median(w[live]) if live.any() else 0.0:.4f} mm, "
          f"p90 {np.percentile(w[live], 90) if live.any() else 0.0:.4f} mm, "
          f"max {w.max():.4f} mm", flush=True)


print("THE STALENESS, IN MILLIMETRES: accumulated association displacement that was refused")
print(f"  {NH} held nodes; the jump limit is {1000*5e-4:.1f} mm per node per step\n")
err = None
try:
    seat_on_bed(region, base, closest, load_steps=S.LOAD_STEPS, gap_tol_m=S.GAP_TOL_M,
                jump_limit_m=5e-4, move=move, assoc_tol_m=S.ASSOC_TOL_M,
                on_step=report, log=logger)
except Cap as e: err = str(e)
except Exception as e: err = f"{type(e).__name__}: {e}"

print("\n" + "=" * 100)
if not seen:
    print("NOT MEASURED -- the drive accepted no step" + (f": {err}" if err else ""))
    raise SystemExit(1)
last = seen[-1]; t, k = last["t"], last["k"]
ever = k > 0
print(f"at fraction {last['fraction']:.4f}, after {len(seen)} accepted steps, the accumulated "
      "refused association travel:")
print(f"  over ALL {NH} held nodes:      median {np.median(t):.4f} mm, p90 {np.percentile(t,90):.4f} mm, max {t.max():.4f} mm")
if ever.any():
    print(f"  over the {int(ever.sum())} ever refused:  median {np.median(t[ever]):.4f} mm, "
          f"p90 {np.percentile(t[ever],90):.4f} mm, max {t[ever].max():.4f} mm")
print(f"  nodes over 1 mm: {int((t > 1.0).sum())} of {NH} ({100*(t > 1.0).mean():.1f}%); "
      f"over 0.5 mm: {int((t > 0.5).sum())} ({100*(t > 0.5).mean():.1f}%)")
# THE BOUNDARY BETWEEN THE TWO BRANCHES, declared here before the run: 1 mm. The coordinator's
# branches are "a few tenths" and "millimetres", and 1 mm is where those meet. It is also the bed's
# own facet scale, the length below which this surface cannot distinguish anything -- so a staleness
# under it cannot point a constraint at a different feature, and one over it can.
BOUNDARY_MM = 1.0
w = last.get("w"); live = w > 0 if w is not None else None
if w is not None and live.any():
    print(f"  INSTANTANEOUS staleness, the distance a refused constraint is from the bed NOW, over "
          f"the {int(live.sum())} refused this step:")
    print(f"    median {np.median(w[live]):.4f} mm, p90 {np.percentile(w[live],90):.4f} mm, "
          f"max {w.max():.4f} mm; over 1 mm: {int((w > 1.0).sum())} of {NH} "
          f"({100*(w > 1.0).mean():.1f}% of all held nodes)")
# THE VERDICT IS TAKEN ON THE INSTANTANEOUS NUMBER, not the accumulated one, because the
# accumulated one counts a repeatedly-refused teleport once per refusal.
affected = float(np.median(w[live])) if (w is not None and live.any()) else 0.0
if affected >= BOUNDARY_MM:
    print(f"  -> MILLIMETRES ({affected:.4f} mm median over the affected nodes, against the declared "
          f"{BOUNDARY_MM:g} mm): the anatomical results on this line were computed against "
          "constraints pointing where the bed no longer is, and are re-run rather than re-read")
else:
    print(f"  -> A FEW TENTHS ({affected:.4f} mm median over the affected nodes, under the declared "
          f"{BOUNDARY_MM:g} mm): the 31-37% refusal is a bookkeeping fact, and the anatomical "
          "seating results survive with a caveat naming it")
print("  NOTE: this is the state at the LAST accepted step, not the whole drive -- the drive was "
      "still advancing" + (f" when it stopped: {err}" if err else ""))
