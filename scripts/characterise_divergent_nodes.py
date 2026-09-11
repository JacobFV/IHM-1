"""SIZING THE 39: which held nodes move AWAY from their target, and why.

At fraction 0.0625 the drive closes +0.4540 mm at the median -- 97% of the scheduled 0.468 mm -- but
39 of 3,123 held nodes are FURTHER from full seating than they started. 1.25% sits close enough to
this line's own 1% gates to be a number rather than an impression, and the maximum was one of these.

FOUR QUESTIONS, each answerable and none of them a threshold:
  WHERE      are they on the base -- clustered on one patch, or scattered?
  THE SAME?  are they the same nodes at every fraction, or a rotating cast?
  REFUSED?   do they coincide with the nodes whose association update was refused at the 0.5 mm jump
             limit? The refusal is per-node and was invisible in the log, because 'association moved
             0.4999 mm' is the maximum over ACCEPTED updates and a refused node contributes nothing
             to it. A coincidence here connects them to a cause already on the record.
  HOW FAR    have they gone backwards, in millimetres, against the travel they were asked for?

Chance agreement is reported alongside every overlap: with 39 divergent nodes and R refused ones out
of 3,123, the overlap expected at random is 39*R/3123. An overlap is only evidence if it beats that,
and quoting an overlap without it would be the same error as an inversion count without a
denominator.
"""
import sys, importlib.util, time
sys.path.insert(0, "/home/brandonin/Documents/IHM-1")
CAP_S = float(sys.argv[sys.argv.index("--cap") + 1]) if "--cap" in sys.argv else 2400.0
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
held_idx = np.flatnonzero(held)                      # indices INTO base
start = np.abs(move[held])
NH = int(held.sum())

start_t = time.time()
class Cap(RuntimeError): pass
def logger(m, flush=True):
    if time.time() - start_t > CAP_S: raise Cap(f"wall-clock cap {CAP_S:.0f}s")

seen = []


def record(fraction, advance, u, aim, info):
    to_aim = info["to_aim_m"]
    bad = to_aim > start + 1e-12
    refused = info["refused_mask"][held_idx]
    both = int((bad & refused).sum())
    chance = float(bad.sum()) * float(refused.sum()) / NH
    pts = (X + u)[base][held_idx][bad]
    spread = np.linalg.norm(pts - pts.mean(0), axis=1) if len(pts) else np.zeros(1)
    seen.append(dict(fraction=fraction, bad=np.flatnonzero(bad), refused_n=int(refused.sum()),
                     both=both, chance=chance,
                     backwards_mm=1000 * float((to_aim - start)[bad].max()) if bad.any() else 0.0,
                     travel_mm=1000 * float(start[bad].mean()) if bad.any() else 0.0,
                     spread_mm=1000 * float(spread.max())))
    r = seen[-1]
    print(f"  fraction {fraction:.4f}: {len(r['bad'])} of {NH} held nodes moving away, worst by "
          f"{r['backwards_mm']:.3f} mm; their own asked travel averages {r['travel_mm']:.3f} mm; "
          f"they span {r['spread_mm']:.1f} mm on the base", flush=True)
    print(f"      association refused for {r['refused_n']} of {NH} nodes this step; "
          f"{r['both']} of the {len(r['bad'])} diverging are among them, against {r['chance']:.1f} "
          f"expected by chance", flush=True)
    if len(seen) > 1:
        prev = set(seen[-2]["bad"].tolist()); now = set(r["bad"].tolist())
        keep = len(prev & now)
        print(f"      {keep} of {len(prev)} are the same nodes as at the previous fraction "
              f"({100*keep/max(len(prev),1):.0f}% persistent)", flush=True)


print("SIZING THE 39: which held nodes move away from their target, and whether they are the "
      "nodes whose association was refused")
print(f"  {NH} held nodes of {len(base)} base nodes; median asked travel "
      f"{1000*np.median(start):.3f} mm\n")
err = None
try:
    seat_on_bed(region, base, closest, load_steps=S.LOAD_STEPS, gap_tol_m=S.GAP_TOL_M,
                jump_limit_m=5e-4, move=move, assoc_tol_m=S.ASSOC_TOL_M,
                on_step=record, log=logger)
except Cap as e: err = str(e)
except Exception as e: err = f"{type(e).__name__}: {e}"

print("\n" + "=" * 100)
if not seen:
    print(f"NOT CHARACTERISED -- the drive accepted no step" + (f": {err}" if err else ""))
    raise SystemExit(1)
core = set(seen[0]["bad"].tolist())
for r in seen[1:]: core &= set(r["bad"].tolist())
print(f"{len(core)} nodes diverge at EVERY accepted fraction -- the persistent core of a set that "
      f"runs {min(len(r['bad']) for r in seen)} to {max(len(r['bad']) for r in seen)} nodes "
      f"over {len(seen)} accepted steps")
tot_both = sum(r["both"] for r in seen); tot_chance = sum(r["chance"] for r in seen)
print(f"overlap with the refused-association nodes, summed over steps: {tot_both} against "
      f"{tot_chance:.1f} expected by chance"
      + ("  -> the refusal is implicated" if tot_both > 2 * tot_chance else
         "  -> NOT distinguishable from chance; the refusal is not their cause"))
if err: print(f"\n(the drive stopped: {err})")
