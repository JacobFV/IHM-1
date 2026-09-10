"""decisive test of the partition explanation.

if radius/ulna/patella fail because each owns only a STRIP or PATCH of a shared
skin region, then giving the shared region to both -- one forearm, one foot, one
knee -- must bring their bones inside.  if merged pieces still fail, the partition
is not the explanation and something about the skin itself is.
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

def piece(members):
    idx = [segs.index(m) for m in members]
    sel = ext[np.isin(owner[ext], idx)]
    used, inv = np.unique(F[sel], return_inverse=True); lv, lf = V[used], inv.reshape(-1, 3)
    if B.simtk_precondition(lv, lf) is not None:
        cv, cf = B.repair(lv, lf); lv, lf, _ = B.cap_boundaries(cv, cf)
    return lv, lf

tests = [  # (skin pieces merged, bones tested)
    (['radius_l'], 'radius_l'), (['ulna_l'], 'ulna_l'),
    (['radius_l', 'ulna_l'], 'radius_l'), (['radius_l', 'ulna_l'], 'ulna_l'),
    (['radius_l', 'ulna_l', 'hand_l'], 'radius_l'),
    (['patella_l'], 'patella_l'), (['patella_l', 'tibia_l'], 'patella_l'),
    (['patella_l', 'tibia_l', 'femur_l'], 'patella_l'),
    (['calcn_l'], 'calcn_l'), (['calcn_l', 'toes_l'], 'calcn_l'),
    (['femur_l'], 'femur_l'), (['femur_l', 'pelvis'], 'femur_l'),
]
cache = {}
print(f"{'skin pieces merged':34s} {'bones':10s} {'inside':>7s}")
for members, b in tests:
    k = tuple(members)
    if k not in cache: cache[k] = piece(members)
    print(f"{'+'.join(members):34s} {b:10s} {B.enclosure(*cache[k], bones[b], samples=800):7.3f}", flush=True)
