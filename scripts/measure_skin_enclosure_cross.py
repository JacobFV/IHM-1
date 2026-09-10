"""which skin piece DOES contain each segment's own bones?

prediction if the failure is the PARTITION: a bone that fails its own piece sits
inside a NEIGHBOUR's piece (radius -> hand or humerus), and reads ~0 in pieces far
away.  the far pieces are the negative control; hand_l is the positive one.
"""
import importlib.util, json, gzip, sys
import numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location('bscm', ROOT/'scripts/build_skin_contact_meshes.py')
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
mech = json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text())
ents = {e['id']: e for e in mech['entities']}
skin = next(e for e in mech['entities'] if e['role'] == 'skin')
g = json.loads(gzip.decompress((ROOT/skin['reference_geometry']['path']).read_bytes()))
V = np.asarray(g['positions'], float).reshape(-1, 3); F = np.asarray(g['indices'], np.int64).reshape(-1, 3)
ext = np.asarray(json.loads((ROOT/B.EVIDENCE).read_text())['contact_eligible_triangle_ids'], np.int64)
sb = json.loads(gzip.decompress((ROOT/B.BINDING).read_bytes()))
segs = [s['id'] for s in sb['segments']]; W = np.asarray(sb['weights'], np.float32)
owner = ((W[F[:, 0]] + W[F[:, 1]] + W[F[:, 2]]) / 3).argmax(1)
bind = json.loads((ROOT/'data/derived/anatomy-segment-binding/binding.json').read_text())
rng = np.random.default_rng(0); bones = {}
for eid, rec in bind['entities'].items():
    if ents.get(eid, {}).get('role') != 'rigid_bone': continue
    p = np.asarray(json.loads(gzip.decompress((ROOT/ents[eid]['reference_geometry']['path']).read_bytes()))['positions'], float).reshape(-1, 3)
    bones.setdefault(rec['segment'], []).append(p[rng.choice(len(p), min(len(p), 400), replace=False)])
bones = {k: np.concatenate(v) for k, v in bones.items()}

pieces = {}
print(f"{'piece':10s} {'tris':>7s} {'area m2':>8s}  y-extent of skin    y-extent of own bones")
for i, s in enumerate(segs):
    sel = ext[owner[ext] == i]
    if len(sel) < 64: continue
    used, inv = np.unique(F[sel], return_inverse=True); lv, lf = V[used], inv.reshape(-1, 3)
    t = lv[lf]; area = float(np.linalg.norm(np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]), axis=1).sum()/2)
    if B.simtk_precondition(lv, lf) is not None:
        cv, cf = B.repair(lv, lf); lv, lf, _ = B.cap_boundaries(cv, cf)
    pieces[s] = (lv, lf)
    by = bones.get(s); bs = f"{by[:,1].min():+.3f}..{by[:,1].max():+.3f}" if by is not None else "-"
    print(f"{s:10s} {len(sel):7d} {area:8.4f}  {lv[:,1].min():+.3f}..{lv[:,1].max():+.3f}   {bs}", flush=True)

probe = ['hand_l', 'radius_l', 'ulna_l', 'patella_l', 'calcn_l', 'femur_l']
cols = [c for c in segs if c in pieces]
print("\nshare of ROW bones inside COLUMN piece (only cells > 0.02 shown; '.' = ~0)")
print(f"{'bones':10s} " + " ".join(f"{c[:7]:>7s}" for c in cols))
for r in probe:
    vals = [B.enclosure(*pieces[c], bones[r], samples=600) for c in cols]
    print(f"{r:10s} " + " ".join(f"{v:7.2f}" if v > .02 else f"{'.':>7s}" for v in vals), flush=True)
