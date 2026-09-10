"""causal probe: does correcting the JOINT registration bring failing ligaments inside
their own ultimate strain?

measure_ligament_attachment_vs_path.py showed that moving attachments within a
ligament's own footprint rescues at most 13 of 51.  measure_joint_centre_registration.py
showed the global atlas->scaffold map displaces the atlas's joints from the
scaffold's -- by comparable amounts to a cruciate's own length.  correlation is not
cause, so intervene: translate each failing element, rigidly, by the displacement
between the atlas's articular centre and the scaffold's articular centre at the
joint it crosses (the SAME closest-approach rule on both skeletons), re-derive its
tip attachments exactly as the builder does, and re-sweep.

a pure translation per joint is the crudest possible registration correction.  if
it rescues the element, registration is the cause and a proper per-segment
registration is the fix.  if it does not, the straight path is, and the element
needs a wrap surface.  gate: a ZERO translation must reproduce the stored peak
strain to machine precision.
"""
import gzip, importlib.util, json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
def load(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
B = load("btfe", "scripts/build_tissue_force_elements.py")
bscm = load("bscm", "scripts/build_skin_contact_meshes.py")
bind = B.load_module("bind_anatomy", ROOT / "scripts/bind_anatomy_to_segments.py")
render = B.load_module("render_body_3d", ROOT / "scripts/render_body_3d.py")
crawl = B.load_module("crawl", ROOT / "scripts/crawl.py")
entities = json.loads(B.ANATOMY.read_text())["entities"]
gpath = {e["id"]: ROOT / e["reference_geometry"]["path"] for e in entities if e.get("reference_geometry")}
bind.GEOMETRY_PATH.update(gpath)
groups = {}
for e in entities:
    if e["role"] == "rigid_bone": groups.setdefault(bind.named_segment(e["name"]), []).append(e["id"])
segments = sorted(groups)
atlas_bone = {s: np.concatenate([bind.read_geometry(i) for i in groups[s]]) for s in segments}
trees = [cKDTree(atlas_bone[s]) for s in segments]
binding = json.loads(B.BINDING.read_text())
inv = np.linalg.inv(np.asarray(binding["similarity_atlas_from_opensim_ground"], float))
ref = binding["reference_pose_rad"]; model = render.OsimModel(B.MODEL); rest = model.forward(ref)
parent = {j["child"]: j["parent"] for j in model.joints}
ranges = crawl.declared_ranges(B.MODEL)
jc = {}
for j in model.joints:
    for c in j["coords"]:
        if c in ranges: jc.setdefault(j["child"], []).append(c)

# ---- per-joint displacement, same closest-approach rule on both skeletons
def ground(p): return p @ inv[:3, :3].T + inv[:3, 3]
atlas_g = {s: ground(p) for s, p in atlas_bone.items()}
osim_g = {s: p @ rest[s][:3, :3].T + rest[s][:3, 3] for s, p in bscm.bone_clouds().items() if s in rest}
def centre(pa, pc, margin=0.008):
    d, j = cKDTree(pc).query(pa); ok = d <= d.min() + margin
    return ((pa[ok] + pc[j[ok]]) / 2).mean(0)
JOINTS = {}  # keyed by child body of the joint
for par, ch in [("femur_l", "tibia_l"), ("femur_r", "tibia_r"), ("tibia_l", "talus_l"), ("tibia_r", "talus_r"),
                ("humerus_l", "ulna_l"), ("humerus_r", "ulna_r"), ("pelvis", "femur_l"), ("pelvis", "femur_r")]:
    JOINTS[ch] = centre(osim_g[par], osim_g[ch]) - centre(atlas_g[par], atlas_g[ch])
print("atlas articular centre -> scaffold articular centre (same rule both sides):")
for ch, t in JOINTS.items(): print(f"  joint above {ch:9s} {1000*np.linalg.norm(t):6.1f} mm   ({', '.join(f'{1000*x:+.1f}' for x in t)})")

def chain(x):
    out = [x]
    while x in parent: x = parent[x]; out.append(x)
    return out
def crossed(a, b):
    ca, cb = chain(a), chain(b); common = next(x for x in ca if x in cb)
    return ca[:ca.index(common)] + cb[:cb.index(common)]
def spanning(a, b): return [c for body in crossed(a, b) for c in jc.get(body, [])]
cache = {}
def T(c, v):
    k = (c, round(float(v), 12))
    if k not in cache: p = dict(ref); p[c] = float(v); cache[k] = model.forward(p)
    return cache[k]
def peak(pa_atlas, pb_atlas, b1, b2, t):
    la = np.linalg.solve(rest[b1][:3, :3], ground(pa_atlas) + t - rest[b1][:3, 3])
    lb = np.linalg.solve(rest[b2][:3, :3], ground(pb_atlas) + t - rest[b2][:3, 3])
    def L(tr): return np.linalg.norm((tr[b1][:3, :3] @ la + tr[b1][:3, 3]) - (tr[b2][:3, :3] @ lb + tr[b2][:3, 3]))
    s = L(rest)
    return max((L(T(c, v)) - s) / s for c in spanning(b1, b2) for v in np.linspace(*ranges[c], B.ADMISSIBILITY_SAMPLES))

rows = json.loads((B.OUT / "ligaments.json").read_text())["elements"]
bad = [r for r in rows if r.get("status") == "two_segment" and not r["kinematically_admissible"]]
u = B.ULTIMATE_STRAIN; out, gate = [], []
print(f"\n{'element':44s} {'joint':9s} {'stored':>7s} {'t=0':>7s} {'corrected':>9s}")
for r in bad:
    g = json.loads(gzip.decompress(gpath[r["id"]].read_bytes())); V = np.asarray(g["positions"], float).reshape(-1, 3)
    nearest = np.argmin(np.stack([tr.query(V)[0] for tr in trees], 1), 1)
    va, vb = V[nearest == segments.index(r["body1"])], V[nearest == segments.index(r["body2"])]
    def tip(side, other):
        rad = np.linalg.norm(side - other, axis=1); ch = side[rad >= np.percentile(rad, 100 - B.TIP_QUANTILE)]
        return (ch if len(ch) else side).mean(0)
    ta, tb = tip(va, vb.mean(0)), tip(vb, va.mean(0))
    joint = next((c for c in crossed(r["body1"], r["body2"]) if c in JOINTS), None)
    zero = peak(ta, tb, r["body1"], r["body2"], np.zeros(3)); gate.append(abs(zero - r["peak_strain_over_declared_range"]))
    corr = peak(ta, tb, r["body1"], r["body2"], JOINTS[joint]) if joint else None
    out.append(dict(name=r["name"], joint=joint, stored=r["peak_strain_over_declared_range"], corrected=corr))
    print(f"{r['name'][:44]:44s} {joint or '-':9s} {r['peak_strain_over_declared_range']:7.3f} {zero:7.3f} "
          + (f"{corr:9.3f}" if corr is not None else f"{'no joint':>9s}"), flush=True)
print(f"\nGATE: max |t=0 - stored| = {max(gate):.2e} " + ("PASS" if max(gate) < 1e-6 else "FAIL -- nothing below is interpretable"))
covered = [o for o in out if o["corrected"] is not None]
print(f"covered by a measured joint: {len(covered)}/{len(out)}")
print(f"admissible after the per-joint translation: {sum(o['corrected'] <= u for o in covered)}/{len(covered)}"
      f"   (before: 0/{len(covered)})")
print(f"peak strain fell for {sum(o['corrected'] < o['stored'] for o in covered)}/{len(covered)};"
      f" median {np.median([o['stored'] for o in covered]):.3f} -> {np.median([o['corrected'] for o in covered]):.3f}")
(B.OUT / "registration_probe.json").write_text(json.dumps(dict(ultimate_strain=u, gate_max_abs=max(gate),
    joint_translation_m={k: v.tolist() for k, v in JOINTS.items()}, rows=out), indent=2) + "\n")
