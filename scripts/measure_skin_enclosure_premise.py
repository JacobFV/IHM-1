"""premise gate: are this body's OWN anatomical bones inside its OWN skin?

both come from one acquired body in one canonical frame, so no registration is
involved and the answer should be ~1.0.  if it is, fitting the OpenSim bones onto
these per segment carries enclosure across; if it is not, that plan cannot work.
same cut, same caps, same enclosure() as build_skin_contact_meshes.py, so the
number is directly comparable to the 0.273 / 0.445 the registered gate prints.
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
V = np.asarray(g['positions'], float).reshape(-1, 3)
F = np.asarray(g['indices'], np.int64).reshape(-1, 3)
exterior = np.asarray(json.loads((ROOT/B.EVIDENCE).read_text())['contact_eligible_triangle_ids'], np.int64)
sb = json.loads(gzip.decompress((ROOT/B.BINDING).read_bytes()))
segments = [s['id'] for s in sb['segments']]
W = np.asarray(sb['weights'], np.float32)
owner = ((W[F[:, 0]] + W[F[:, 1]] + W[F[:, 2]]) / 3).argmax(1)

bind = json.loads((ROOT/'data/derived/anatomy-segment-binding/binding.json').read_text())
rng = np.random.default_rng(0)
bones = {}
for eid, rec in bind['entities'].items():
    if ents.get(eid, {}).get('role') != 'rigid_bone': continue
    p = np.asarray(json.loads(gzip.decompress((ROOT/ents[eid]['reference_geometry']['path']).read_bytes()))['positions'], float).reshape(-1, 3)
    if len(p) > 400: p = p[rng.choice(len(p), 400, replace=False)]
    bones.setdefault(rec['segment'], []).append(p)
bones = {k: np.concatenate(v) for k, v in bones.items()}

# positive control in the output itself: hand, pelvis and torso must read ~1.0
print(f"{'segment':12s} {'bone pts':>9s} {'own bones inside own skin':>26s}")
rows = []
for i, seg in enumerate(segments):
    sel = exterior[owner[exterior] == i]
    if len(sel) < 64 or seg not in bones: continue
    used, inv = np.unique(F[sel], return_inverse=True)
    lv, lf = V[used], inv.reshape(-1, 3)
    if B.simtk_precondition(lv, lf) is not None:
        cv, cf = B.repair(lv, lf)
        try: lv, lf, _ = B.cap_boundaries(cv, cf)
        except Exception as e: print(f"{seg:12s} capping failed: {e}"); continue
    e = B.enclosure(lv, lf, bones[seg])
    rows.append(e); print(f"{seg:12s} {len(bones[seg]):9d} {e:26.3f}", flush=True)
print(f"\nmean {np.mean(rows):.3f} over {len(rows)} segments; {sum(r >= .99 for r in rows)} enclose >= 0.99")
print("registered OpenSim bones for comparison: canonical map 0.273, binding map 0.445, 0/20 each")
