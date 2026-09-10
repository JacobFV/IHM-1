"""is a failing ligament's ATTACHMENT wrong, or its PATH?

51 of 117 tensile elements pass ligament ultimate strain (17.1%) inside the
declared range of a joint they span -- the cruciates at 118-123%.  Two causes are
named in docs/TISSUE_MECHANICS.md and they call for different fixes:

  attachment  the builder puts each end at a TIP centroid (the farthest 25% of
              the vertices voting for that bone), which for a cruciate sits
              ~22 mm off the flexion axis.  fix: a different attachment rule.
  path        a straight line is not a ligament's path; it wraps.  fix: a wrap
              surface per joint.

This separates them.  For every inadmissible element, peak strain over the same
sweep the builder runs, under three attachment choices:

  (i)   the builder's tip centroids.  MUST reproduce the stored peak strain, or
        this replication differs from the builder and nothing below is
        interpretable.  That is the gate.
  (ii)  the whole-footprint centroid of each end.
  (iii) the BEST pair of vertices inside the ligament's own footprint -- the most
        isometric fibre the ligament's own geometry contains.  An optimistic
        bound on what moving attachments alone can do, and a selected one, so it
        is reported with the share of the range over which that fibre is taut:
        a fibre that is admissible because it is slack everywhere restrains
        nothing, and the doc already says the 66 survivors largely are that.

If (iii) clears 17.1% for an element, its attachment rule is the defect.  If even
the best fibre in its own footprint cannot, the straight path is, and it needs a
wrap surface.
"""
import gzip, importlib.util, json, sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("btfe", ROOT / "scripts/build_tissue_force_elements.py")
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)
CANDIDATES = 60

bind = B.load_module("bind_anatomy", ROOT / "scripts/bind_anatomy_to_segments.py")
render = B.load_module("render_body_3d", ROOT / "scripts/render_body_3d.py")
crawl = B.load_module("crawl", ROOT / "scripts/crawl.py")
entities = json.loads(B.ANATOMY.read_text())["entities"]
by_id = {e["id"]: e for e in entities}
geometry_path = {e["id"]: ROOT / e["reference_geometry"]["path"] for e in entities if e.get("reference_geometry")}
bind.GEOMETRY_PATH.update(geometry_path)
groups = {}
for e in entities:
    if e["role"] == "rigid_bone":
        groups.setdefault(bind.named_segment(e["name"]), []).append(e["id"])
segments = sorted(groups)
trees = [cKDTree(np.concatenate([bind.read_geometry(i) for i in groups[s]])) for s in segments]

binding = json.loads(B.BINDING.read_text())
inverse = np.linalg.inv(np.asarray(binding["similarity_atlas_from_opensim_ground"], float))
reference_pose = binding["reference_pose_rad"]
model = render.OsimModel(B.MODEL)
rest = model.forward(reference_pose)
parent = {j["child"]: j["parent"] for j in model.joints}
ranges = crawl.declared_ranges(B.MODEL)
joint_coordinates = {}
for j in model.joints:
    for c in j["coords"]:
        if c in ranges: joint_coordinates.setdefault(j["child"], []).append(c)

def to_local(p, seg):
    w = p @ inverse[:3, :3].T + inverse[:3, 3]; T = rest[seg]
    return (w - T[:3, 3]) @ np.linalg.inv(T[:3, :3]).T

def spanning(a, b):
    def chain(x):
        out = [x]
        while x in parent: x = parent[x]; out.append(x)
        return out
    ca, cb = chain(a), chain(b); common = next(x for x in ca if x in cb)
    out = []
    for body in ca[:ca.index(common)] + cb[:cb.index(common)]: out += joint_coordinates.get(body, [])
    return out

pose_cache = {}
def transforms(coord, value):
    k = (coord, round(float(value), 12))
    if k not in pose_cache:
        pose = dict(reference_pose); pose[coord] = float(value); pose_cache[k] = model.forward(pose)
    return pose_cache[k]

def sweep(pa, pb, b1, b2):
    """strain of every fibre pa[i] -> pb[j] at every sweep pose: (poses, na, nb)."""
    def dist(T):
        A = pa @ T[b1][:3, :3].T + T[b1][:3, 3]; Bp = pb @ T[b2][:3, :3].T + T[b2][:3, 3]
        return np.linalg.norm(A[:, None] - Bp[None], axis=2)
    slack = dist(rest)
    out = [(dist(transforms(c, v)) - slack) / slack
           for c in spanning(b1, b2) for v in np.linspace(*ranges[c], B.ADMISSIBILITY_SAMPLES)]
    return np.stack(out)

rows = json.loads((B.OUT / "ligaments.json").read_text())["elements"]
bad = [r for r in rows if r.get("status") == "two_segment" and not r["kinematically_admissible"]]
rng = np.random.default_rng(0)
print(f"{len(bad)} inadmissible elements; ultimate strain {B.ULTIMATE_STRAIN}\n")
print(f"{'element':44s} {'stored':>7s} {'(i) tip':>8s} {'(ii) cen':>8s} {'(iii) best':>10s} {'taut':>5s}")
res, gate = [], []
for r in bad:
    g = json.loads(gzip.decompress(geometry_path[r["id"]].read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3)
    nearest = np.argmin(np.stack([t.query(V)[0] for t in trees], 1), 1)
    a, b = segments.index(r["body1"]), segments.index(r["body2"])
    va, vb = V[nearest == a], V[nearest == b]
    def tip(side, other):
        rad = np.linalg.norm(side - other, axis=1); ch = side[rad >= np.percentile(rad, 100 - B.TIP_QUANTILE)]
        return (ch if len(ch) else side).mean(0)
    ta, tb = tip(va, vb.mean(0)), tip(vb, va.mean(0))
    def cands(side, extra):
        pick = side[rng.choice(len(side), min(len(side), CANDIDATES), replace=False)]
        return np.vstack([extra, pick])
    ca, cb = cands(va, [ta, va.mean(0)]), cands(vb, [tb, vb.mean(0)])
    S = sweep(to_local(ca, r["body1"]), to_local(cb, r["body2"]), r["body1"], r["body2"])
    peak = S.max(0)
    i_tip, ii_cen = float(peak[0, 0]), float(peak[1, 1])
    k = np.unravel_index(np.argmin(peak), peak.shape); iii = float(peak[k])
    taut = float((S[:, k[0], k[1]] > 1e-9).mean())
    gate.append(abs(i_tip - r["peak_strain_over_declared_range"]))
    res.append((r["name"], r["peak_strain_over_declared_range"], i_tip, ii_cen, iii, taut))
    print(f"{r['name'][:44]:44s} {res[-1][1]:7.3f} {i_tip:8.3f} {ii_cen:8.3f} {iii:10.3f} {taut:5.2f}", flush=True)

u = B.ULTIMATE_STRAIN; R = np.array([x[1:] for x in res])
print(f"\nGATE: max |recomputed tip peak - stored| = {max(gate):.2e}  "
      + ("PASS" if max(gate) < 1e-6 else "FAIL -- this replication is not the builder's; nothing below is interpretable"))
print(f"admissible at {u}:  (i) tip {int((R[:,1]<=u).sum())}/{len(R)}   (ii) centroid {int((R[:,2]<=u).sum())}/{len(R)}   "
      f"(iii) best fibre {int((R[:,3]<=u).sum())}/{len(R)}")
print(f"  of those rescued by (iii), taut over >25% of the sweep: {int(((R[:,3]<=u)&(R[:,4]>.25)).sum())}")
out = ROOT / "data/derived/tissue-force-elements-v1/attachment_vs_path.json"
out.write_text(json.dumps(dict(ultimate_strain=u, candidates_per_end=CANDIDATES, gate_max_abs=max(gate),
    rows=[dict(name=n, stored=s, tip=t, centroid=c, best_fibre=bf, best_fibre_taut_share=tt)
          for n, s, t, c, bf, tt in res]), indent=2) + "\n")
print("wrote", out.relative_to(ROOT))
