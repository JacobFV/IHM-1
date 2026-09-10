"""the registration gate WITHOUT the partition confound.

per-segment enclosure cannot pass for radius/ulna/patella under ANY registration,
because those segments own a strip or a patch of skin, not a closed region.  so
test each segment's registered OpenSim bones against the WHOLE capped skin.
ceiling: the body's own anatomical bones against the same whole skin, same frame.
"""
import argparse, importlib.util, json, gzip, sys
import numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
_ap = argparse.ArgumentParser()
_ap.add_argument('--warp', type=Path, default=None,
                 help='a scripts/skin_warp.py warp (.npz) whose base must be the binding map; the warped skin is '
                      'measured beside the own-bones ceiling and the binding map. Without it the output is unchanged.')
_ap.add_argument('--json', type=Path, default=None, help='with --warp: write the per-segment columns here')
ARGS = _ap.parse_args()
spec = importlib.util.spec_from_file_location('bscm', ROOT/'scripts/build_skin_contact_meshes.py')
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
mech = json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text())
ents = {e['id']: e for e in mech['entities']}
skin = next(e for e in mech['entities'] if e['role'] == 'skin')
g = json.loads(gzip.decompress((ROOT/skin['reference_geometry']['path']).read_bytes()))
V = np.asarray(g['positions'], float).reshape(-1, 3); F = np.asarray(g['indices'], np.int64).reshape(-1, 3)
ext = np.asarray(json.loads((ROOT/B.EVIDENCE).read_text())['contact_eligible_triangle_ids'], np.int64)
used, inv = np.unique(F[ext], return_inverse=True); sv, sf = V[used], inv.reshape(-1, 3)
if B.simtk_precondition(sv, sf) is not None:
    cv, cf = B.repair(sv, sf); sv, sf, cap = B.cap_boundaries(cv, cf)
    print(f"whole exterior skin capped: {len(sf):,} triangles, {cap:.4f} m2 of cap added")
bind = json.loads((ROOT/'data/derived/anatomy-segment-binding/binding.json').read_text())
rng = np.random.default_rng(0); own = {}
for eid, rec in bind['entities'].items():
    if ents.get(eid, {}).get('role') != 'rigid_bone': continue
    p = np.asarray(json.loads(gzip.decompress((ROOT/ents[eid]['reference_geometry']['path']).read_bytes()))['positions'], float).reshape(-1, 3)
    own.setdefault(rec['segment'], []).append(p[rng.choice(len(p), min(len(p), 200), replace=False)])
own = {k: np.concatenate(v) for k, v in own.items()}
osim = B.bone_clouds()

# map 1: binding similarity, at the pose it was fitted in
Tb, frames_b, _ = B.binding_registration()
# map 2: CanonicalRegistration.global_map, at the reference run's t=0 pose
ref = json.loads((ROOT/B.DEFAULT_REFERENCE).read_text())
Tc = np.linalg.inv(B.CanonicalRegistration(mech, ref).global_map)
frames_c = {n: np.asarray(v['transform_ground'], float) for n, v in ref['bodies'].items()}

def skin_in(T): return sv @ T[:3, :3].T + T[:3, 3]
def placed(seg, frames): M = frames[seg]; return osim[seg] @ M[:3, :3].T + M[:3, 3]
skin_b, skin_c = skin_in(Tb), skin_in(Tc)

segs = sorted(s for s in own if s in osim and s in frames_b and s in frames_c)
if ARGS.warp is not None:
    # the warped skin, beside the ceiling and the map it is built on.  Same capped surface, same
    # frames, same enclosure() and samples: only the skin's vertex positions differ.
    sys.path.insert(0, str(ROOT / 'scripts'))
    from skin_warp import Warp
    W = Warp.load(ARGS.warp)
    if not np.array_equal(W.base, Tb): sys.exit('the warp\'s base is not the binding map; nothing measured')
    skin_w = W.apply(sv)
    print(f"warp {ARGS.warp}: {len(W.centres)} centres, max |displacement| at skin {1000*np.abs(skin_w-skin_b).max():.2f} mm")
    print(f"\n{'segment':12s} {'own bones (ceiling)':>20s} {'binding map':>12s} {'warped':>8s}")
    rows = []
    for s in segs:
        a = B.enclosure(sv, sf, own[s], samples=500)
        b = B.enclosure(skin_b, sf, placed(s, frames_b), samples=500)
        w = B.enclosure(skin_w, sf, placed(s, frames_b), samples=500)
        rows.append((a, b, w)); print(f"{s:12s} {a:20.3f} {b:12.3f} {w:8.3f}", flush=True)
    r = np.array(rows)
    print(f"\n{'mean':12s} {r[:,0].mean():20.3f} {r[:,1].mean():12.3f} {r[:,2].mean():8.3f}")
    print(f"{'>= 0.99':12s} {int((r[:,0]>=.99).sum()):>17d}/{len(r)} {int((r[:,1]>=.99).sum()):>9d}/{len(r)} {int((r[:,2]>=.99).sum()):>5d}/{len(r)}")
    if ARGS.json is not None:
        ARGS.json.write_text(json.dumps(dict(warp=str(ARGS.warp), segments=segs,
            own_bones_ceiling={s: x[0] for s, x in zip(segs, rows)}, binding_map={s: x[1] for s, x in zip(segs, rows)},
            warped={s: x[2] for s, x in zip(segs, rows)},
            mean=dict(own_bones_ceiling=float(r[:,0].mean()), binding_map=float(r[:,1].mean()), warped=float(r[:,2].mean()))), indent=2) + "\n")
    sys.exit(0)
print(f"\n{'segment':12s} {'own bones (ceiling)':>20s} {'canonical map':>14s} {'binding map':>12s}")
rows = []
for s in segs:
    a = B.enclosure(sv, sf, own[s], samples=500)
    c = B.enclosure(skin_c, sf, placed(s, frames_c), samples=500)
    b = B.enclosure(skin_b, sf, placed(s, frames_b), samples=500)
    rows.append((a, c, b)); print(f"{s:12s} {a:20.3f} {c:14.3f} {b:12.3f}", flush=True)
r = np.array(rows)
print(f"\n{'mean':12s} {r[:,0].mean():20.3f} {r[:,1].mean():14.3f} {r[:,2].mean():12.3f}")
print(f"{'>= 0.99':12s} {int((r[:,0]>=.99).sum()):>17d}/{len(r)} {int((r[:,1]>=.99).sum()):>11d}/{len(r)} {int((r[:,2]>=.99).sum()):>9d}/{len(r)}")
print("per-segment (partition-confounded) gate, for comparison: canonical 0.273, binding 0.445")
